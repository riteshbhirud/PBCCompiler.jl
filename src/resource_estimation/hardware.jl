"""
    ResourceEstimation.grid(rows, cols) -> Graphs.SimpleGraph

Square-lattice connectivity with node `(r, c)` at index `(r - 1) * cols + c`, coupled
to its horizontal and vertical neighbours. Not exported, to avoid clashing with `Graphs.grid`.
"""
function grid(rows::Int, cols::Int)
    g = Graphs.SimpleGraph(rows * cols)
    id(r, c) = (r - 1) * cols + c
    for r in 1:rows, c in 1:cols
        c < cols && Graphs.add_edge!(g, id(r, c), id(r, c + 1))
        r < rows && Graphs.add_edge!(g, id(r, c), id(r + 1, c))
    end
    g
end

"""
    heavy_hex(n) -> Graphs.SimpleGraph

Heavy-hexagon connectivity, as on IBM devices: an `n`-by-`n` patch of the hexagonal
lattice with an extra qubit on every edge. Every node has degree at most 3.
"""
function heavy_hex(n::Int)
    n >= 1 || throw(ArgumentError("n must be ≥ 1"))
    rows = n + 1; cols = 2n + 2
    id(r, c) = (r - 1) * cols + c
    edges = Tuple{Int,Int}[]
    for r in 1:rows, c in 1:cols
        c < cols && push!(edges, (id(r, c), id(r, c + 1)))
        r < rows && iseven(c + r) && push!(edges, (id(r, c), id(r + 1, c)))
    end
    nv = rows * cols
    g = Graphs.SimpleGraph(nv + length(edges))
    for (k, (a, b)) in enumerate(edges)
        m = nv + k
        Graphs.add_edge!(g, a, m); Graphs.add_edge!(g, m, b)
    end
    g
end

function hardware_graph(conn, nodes::Int)
    conn isa Graphs.SimpleGraph && return conn
    if conn == :grid
        s = ceil(Int, sqrt(nodes))
        return grid(s, s)
    elseif conn == :heavy_hex
        n = 1
        while Graphs.nv(heavy_hex(n)) < nodes
            n += 1
        end
        return heavy_hex(n)
    end
    throw(ArgumentError("no hardware graph for connectivity $conn"))
end

"""
    tanner_graph(m::CodeModel) -> Graphs.SimpleGraph

The code's Tanner graph: data qubits are nodes `1:n`, checks are nodes `n+1:n+s`,
and each check is joined to the data qubits it acts on.
"""
function tanner_graph(m::CodeModel)
    Hx, Hz = parity_check_matrices(m)
    tanner_graph(vcat(Hx, Hz))
end

function tanner_graph(H::AbstractMatrix{Bool})
    s, n = size(H)
    g = Graphs.SimpleGraph(n + s)
    for r in 1:s, c in 1:n
        H[r, c] && Graphs.add_edge!(g, c, n + r)
    end
    g
end

"""
    EmbeddingResult

Placement of a code's Tanner graph onto hardware: `placement[v]` is the hardware
node of Tanner node `v`, `cost` is the total shortest-path length of the Tanner
edges, and `edge_lengths` lists the hardware distance of every Tanner edge.
"""
struct EmbeddingResult
    placement::Vector{Int}
    cost::Int
    edge_lengths::Vector{Int}
end

function _all_pairs_distances(hw::Graphs.SimpleGraph)
    H = Graphs.nv(hw)
    D = Matrix{Int}(undef, H, H)
    for s in 1:H
        D[s, :] .= Graphs.gdistances(hw, s)
    end
    D
end

"""
    embed_tanner_graph(code_graph, hw_graph; iterations = 20_000, rng, temperature = 2.0) -> EmbeddingResult

Place the Tanner graph injectively on hardware nodes, minimizing the total
shortest-path length of its edges by simulated annealing. The result is
deterministic for a given `rng` seed and iteration count.
"""
function embed_tanner_graph(G::Graphs.SimpleGraph, hw::Graphs.SimpleGraph; iterations::Int = 20_000,
                            rng::AbstractRNG = Random.MersenneTwister(0), temperature::Float64 = 2.0)
    V = Graphs.nv(G); H = Graphs.nv(hw)
    V <= H || throw(ArgumentError("hardware has $H qubits, the code needs $V"))
    D = _all_pairs_distances(hw)
    any(==(typemax(Int)), D) && throw(ArgumentError("hardware graph is disconnected"))
    seen = falses(V); oc = Int[]
    for s in 1:V
        seen[s] && continue
        seen[s] = true; queue = [s]
        while !isempty(queue)
            v = popfirst!(queue); push!(oc, v)
            for u in Graphs.neighbors(G, v)
                seen[u] || (seen[u] = true; push!(queue, u))
            end
        end
    end
    center = argmin([sum(D[h, :]) for h in 1:H])
    oh = sortperm(D[center, :])
    place = zeros(Int, V); occupant = zeros(Int, H)
    for (k, v) in enumerate(oc)
        place[v] = oh[k]; occupant[oh[k]] = v
    end
    edges = [(Graphs.src(e), Graphs.dst(e)) for e in Graphs.edges(G)]
    nbrs = [Graphs.neighbors(G, v) for v in 1:V]
    cost = sum(D[place[a], place[b]] for (a, b) in edges; init = 0)
    local_cost(v, at) = sum(D[at, place[u]] for u in nbrs[v]; init = 0)
    for it in 1:iterations
        T = temperature * 0.01^(it / iterations)
        v = rand(rng, 1:V)
        h = rand(rng, 1:H)
        h == place[v] && continue
        w = occupant[h]
        old = place[v]
        if w == 0
            Δ = local_cost(v, h) - local_cost(v, old)
        else
            # an edge between v and w keeps its length under a swap; correct the local costs for it
            Δ = local_cost(v, h) - local_cost(v, old) + local_cost(w, old) - local_cost(w, h)
            w in nbrs[v] && (Δ += 2 * D[h, old])
        end
        if Δ <= 0 || rand(rng) < exp(-Δ / T)
            if w == 0
                occupant[old] = 0
            else
                place[w] = old; occupant[old] = w
            end
            place[v] = h; occupant[h] = v
            cost += Δ
        end
    end
    lens = [D[place[a], place[b]] for (a, b) in edges]
    sum(lens; init = 0) == cost || error("internal: annealing cost drifted")
    EmbeddingResult(place, cost, lens)
end

"""
    EmbeddingOverhead

Routing overhead of running a code on non-all-to-all hardware: extra routing qubits,
extra two-qubit-gate layers per syndrome round, the number of non-native Tanner
edges, the routing `model` (`:bridge` or `:swap`), and the embedding itself.
"""
struct EmbeddingOverhead
    extra_qubits::Int
    extra_2q_layers::Int
    non_native_edges::Int
    model::Symbol
    embedding::Union{Nothing,EmbeddingResult}
end

"""
    embedding_overhead(m::CodeModel, hw::HardwareModel; routing = :bridge, iterations, rng) -> EmbeddingOverhead

Overhead of embedding the code's Tanner graph into the hardware connectivity. A
Tanner edge of hardware length `L > 1` is non-native.

- `:bridge`: `L - 1` bridge qubits per non-native edge and `2 * (max L - 1)` extra
  CNOT layers per round.
- `:swap`: no extra qubits and `6 * (max L - 1)` extra CNOT layers per round.

`:all_to_all` hardware has no overhead.
"""
function embedding_overhead(m::CodeModel, hw::HardwareModel; routing::Symbol = :bridge, iterations::Int = 20_000,
                            rng::AbstractRNG = Random.MersenneTwister(0))
    routing in (:bridge, :swap) || throw(ArgumentError("routing must be :bridge or :swap"))
    hw.connectivity === :all_to_all && return EmbeddingOverhead(0, 0, 0, routing, nothing)
    G = tanner_graph(m)
    hwg = hardware_graph(hw.connectivity, Graphs.nv(G))
    emb = embed_tanner_graph(G, hwg; iterations, rng)
    long = filter(>(1), emb.edge_lengths)
    maxL = isempty(long) ? 1 : maximum(long)
    if routing == :bridge
        EmbeddingOverhead(sum(L - 1 for L in long; init = 0), 2 * (maxL - 1), length(long), routing, emb)
    else
        EmbeddingOverhead(0, 6 * (maxL - 1), length(long), routing, emb)
    end
end
