# E1: logical layout of data qubits onto code blocks.

"""
    LogicalLayout(code_blocks, qubit_to_block, qubit_within_block)

Which code block each logical (data) qubit lives in, and its logical index inside
that block. `code_blocks[b]` is the [`CodeModel`](@ref) of block `b`.
"""
struct LogicalLayout
    code_blocks::Vector{CodeModel}
    qubit_to_block::Dict{Int,Int}
    qubit_within_block::Dict{Int,Int}
end

LogicalLayout(blocks::AbstractVector{<:CodeModel}, q2b::AbstractDict, q2i::AbstractDict) =
    LogicalLayout(collect(CodeModel, blocks), Dict{Int,Int}(q2b), Dict{Int,Int}(q2i))

"""
    validate_layout(layout, qubits) -> Vector{String}

Problems with a layout for the given logical qubits (empty if it is valid):
- every qubit is mapped to an existing block;
- no two qubits share a `(block, index)` slot;
- each index lies in `1:k` of its block.
"""
function validate_layout(layout::LogicalLayout, qubits)
    errs = String[]
    nb = length(layout.code_blocks)
    used = Dict{Tuple{Int,Int},Int}()
    for q in sort!(collect(qubits))
        if !haskey(layout.qubit_to_block, q) || !haskey(layout.qubit_within_block, q)
            push!(errs, "logical qubit $q is not mapped to a block")
            continue
        end
        b = layout.qubit_to_block[q]; i = layout.qubit_within_block[q]
        if !(1 <= b <= nb)
            push!(errs, "logical qubit $q is mapped to block $b, but there are $nb blocks")
            continue
        end
        kb = code_k(layout.code_blocks[b])
        (1 <= i <= kb) || push!(errs, "logical qubit $q has index $i in block $b, which encodes only $kb logical qubits")
        if haskey(used, (b, i))
            push!(errs, "logical qubits $(used[(b, i)]) and $q share slot $i of block $b")
        else
            used[(b, i)] = q
        end
    end
    errs
end

"""
    pack_sequential(nqubits, models) -> LogicalLayout

Fill the blocks in order: qubit 1 goes to slot 1 of block 1, and so on, moving to the
next block when a block's `k` slots are full. Throws if the blocks cannot hold
`nqubits`.
"""
function pack_sequential(nqubits::Int, models::AbstractVector{<:CodeModel})
    cap = sum(code_k, models; init = 0)
    cap >= nqubits || throw(ArgumentError("blocks hold $cap logical qubits, need $nqubits"))
    q2b = Dict{Int,Int}(); q2i = Dict{Int,Int}()
    b = 1; i = 0
    for q in 1:nqubits
        if i == code_k(models[b])
            b += 1; i = 0
        end
        i += 1
        q2b[q] = b; q2i[q] = i
    end
    LogicalLayout(models, q2b, q2i)
end

"""
    interaction_graph(circuit, n_data) -> Dict{Tuple{Int,Int},Int}

Logical interaction weights: for every measurement, every pair of data qubits
(`≤ n_data`) in its Pauli support gets +1. Keys are `(a, b)` with `a < b`.
"""
function interaction_graph(c::AbstractVector{CircuitOp.Type}, n_data::Int)
    w = Dict{Tuple{Int,Int},Int}()
    for op in c
        m = _measurement_parts(op)
        m === nothing && continue
        s = [q for q in _pauli_support(m[1], m[3]) if q <= n_data]
        for x in eachindex(s), y in x+1:length(s)
            key = minmax(s[x], s[y])
            w[key] = get(w, key, 0) + 1
        end
    end
    w
end

"""Total interaction weight between qubits in different blocks."""
function cross_block_weight(w::Dict{Tuple{Int,Int},Int}, q2b::Dict{Int,Int})
    sum((v for ((a, b), v) in w if q2b[a] != q2b[b]); init = 0)
end

"""
    pack_by_interaction(circuit, models; n_data = width, refine = true, max_passes = 50) -> (layout, cross_weight)

Assign the circuit's data qubits to blocks so as to keep strongly interacting
qubits together. The interaction graph comes from the compiled circuit's
measurements (edge weight = number of measurements touching both qubits).

The assignment is greedy: qubits are placed in order of decreasing total weight,
each into the block (with a free slot) where it has the most weight already. A
Kernighan–Lin-style refinement then performs the best improving swap or move until
none remains. Returns the layout and the remaining cross-block weight.
"""
function pack_by_interaction(c::AbstractVector{CircuitOp.Type}, models::AbstractVector{<:CodeModel};
                             n_data::Int = _est_width(c), refine::Bool = true, max_passes::Int = 50)
    cap = [code_k(m) for m in models]
    sum(cap) >= n_data || throw(ArgumentError("blocks hold $(sum(cap)) logical qubits, need $n_data"))
    w = interaction_graph(c, n_data)
    nbr = [Dict{Int,Int}() for _ in 1:n_data]
    for ((a, b), v) in w
        nbr[a][b] = v; nbr[b][a] = v
    end
    deg = [sum(values(nbr[q]); init = 0) for q in 1:n_data]
    order = sort(1:n_data; by = q -> (-deg[q], q))
    q2b = Dict{Int,Int}()
    load = zeros(Int, length(models))
    for q in order
        best = 0; bestaff = -1
        for b in eachindex(models)
            load[b] < cap[b] || continue
            aff = sum((v for (r, v) in nbr[q] if get(q2b, r, 0) == b); init = 0)
            if aff > bestaff
                best = b; bestaff = aff
            end
        end
        q2b[q] = best; load[best] += 1
    end
    if refine
        aff(q, b) = sum((v for (r, v) in nbr[q] if q2b[r] == b && r != q); init = 0)
        for _ in 1:max_passes
            bestgain = 0; bestmove = nothing
            for q in 1:n_data
                bq = q2b[q]
                for b in eachindex(models)            # move into a free slot
                    (b == bq || load[b] >= cap[b]) && continue
                    g = aff(q, b) - aff(q, bq)
                    g > bestgain && (bestgain = g; bestmove = (:move, q, b))
                end
                for r in q+1:n_data                   # swap with a qubit in another block
                    br = q2b[r]
                    br == bq && continue
                    wqr = get(nbr[q], r, 0)
                    g = aff(q, br) - aff(q, bq) + aff(r, bq) - aff(r, br) - 2wqr
                    g > bestgain && (bestgain = g; bestmove = (:swap, q, r))
                end
            end
            bestmove === nothing && break
            if bestmove[1] == :move
                _, q, b = bestmove
                load[q2b[q]] -= 1; load[b] += 1; q2b[q] = b
            else
                _, q, r = bestmove
                q2b[q], q2b[r] = q2b[r], q2b[q]
            end
        end
    end
    q2i = Dict{Int,Int}(); next = ones(Int, length(models))
    for q in 1:n_data
        b = q2b[q]; q2i[q] = next[b]; next[b] += 1
    end
    (LogicalLayout(models, q2b, q2i), cross_block_weight(w, q2b))
end
