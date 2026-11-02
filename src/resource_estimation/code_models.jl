# E2: QEC code models.
#
# Sources for every default below:
#  [GoSC]   D. Litinski, "A Game of Surface Codes", Quantum 3, 128 (2019), arXiv:1808.02892.
#  [Fowler] A. G. Fowler et al., "Surface codes: Towards practical large-scale quantum
#           computation", PRA 86, 032324 (2012), arXiv:1208.0928.
#  [BB]     S. Bravyi et al., "High-threshold and low-overhead fault-tolerant quantum
#           memory", Nature 627, 778 (2024), arXiv:2308.07915.
#  [Adapt]  E. Swaroop, T. Jochym-O'Connor, T. J. Yoder, "Universal adapters between
#           quantum LDPC codes", arXiv:2410.03628.
# Anything without a source is a required argument.

"""
    HardwareModel(connectivity; t_1q, t_2q, t_meas, t_reaction, p_phys)

Physical hardware: a connectivity (`:all_to_all`, `:grid`, `:heavy_hex`, or a
`Graphs.SimpleGraph`), gate/measurement times in seconds, the classical
decode-and-feed-forward latency `t_reaction` (seconds), and the physical error
rate `p_phys` used by the logical-error models. All times and rates are required:
they are hardware properties and have no source-backed defaults.
"""
struct HardwareModel
    connectivity::Union{Symbol,Graphs.SimpleGraph{Int}}
    t_1q::Float64
    t_2q::Float64
    t_meas::Float64
    t_reaction::Float64
    p_phys::Float64
end

function HardwareModel(connectivity = :all_to_all; t_1q::Real, t_2q::Real, t_meas::Real, t_reaction::Real, p_phys::Real)
    connectivity isa Symbol && !(connectivity in (:all_to_all, :grid, :heavy_hex)) &&
        throw(ArgumentError("connectivity must be :all_to_all, :grid, :heavy_hex or a SimpleGraph"))
    all(>=(0), (t_1q, t_2q, t_meas, t_reaction)) || throw(ArgumentError("times must be non-negative"))
    0 < p_phys < 1 || throw(ArgumentError("p_phys must be in (0, 1)"))
    HardwareModel(connectivity, t_1q, t_2q, t_meas, t_reaction, p_phys)
end

# ---------------------------------------------------------------- logical error models

"""Logical error rate per code cycle as a function of physical error rate and distance."""
abstract type PLModel end

"""
    AnalyticPL(; A = 0.1, p_th = 0.01)

`p_L(p, d) = A·(p/p_th)^((d+1)/2)` per code cycle. The defaults are [GoSC] Eq. (10),
`p_L = 0.1(100p)^((d+1)/2)` for circuit-level noise in the surface code.
"""
Base.@kwdef struct AnalyticPL <: PLModel
    A::Float64 = 0.1
    p_th::Float64 = 0.01
end
pl(m::AnalyticPL, p::Real, d::Int) = m.A * (p / m.p_th)^((d + 1) / 2)

"""
    FitPL(f)

A user-supplied fit: `p_L = f(p, d)` per code cycle.
"""
struct FitPL{F} <: PLModel
    f::F
end
pl(m::FitPL, p::Real, d::Int) = Float64(m.f(p, d))

"""
    BBFitPL(d_circ, c0, c1, c2)

The bivariate-bicycle fitting formula of [BB] (Sec. 8, Table 6),
`p_L(p) = p^(d_circ/2)·exp(c0 + c1·p + c2·p²)` per syndrome cycle. The distance
argument is ignored (the fit is for one specific code).
"""
struct BBFitPL <: PLModel
    d_circ::Int
    c0::Float64
    c1::Float64
    c2::Float64
end
pl(m::BBFitPL, p::Real, d::Int) = p^(m.d_circ / 2) * exp(m.c0 + m.c1 * p + m.c2 * p^2)

"""
    TablePL(rows; code)
    TablePL(csvpath; code)

Logical error rates from the user's own threshold simulations: rows of
`(code, d, p, p_L, ci_low, ci_high)` (a CSV with exactly these columns and a header),
interpolated linearly in log(p)–log(p_L) at fixed `d`. Querying a `d` that is not
in the table, or a `p` outside the tabulated range, throws.
"""
struct TablePL <: PLModel
    code::String
    points::Dict{Int,Vector{Tuple{Float64,Float64}}}   # d => sorted (p, p_L)
end

function TablePL(rows::AbstractVector; code::AbstractString)
    pts = Dict{Int,Vector{Tuple{Float64,Float64}}}()
    for r in rows
        string(r[1]) == code || continue
        push!(get!(pts, Int(r[2]), Tuple{Float64,Float64}[]), (Float64(r[3]), Float64(r[4])))
    end
    isempty(pts) && throw(ArgumentError("no rows for code $code"))
    foreach(sort!, values(pts))
    TablePL(String(code), pts)
end

function TablePL(path::AbstractString; code::AbstractString)
    lines = filter(!isempty, strip.(readlines(path)))
    hdr = strip.(split(lines[1], ","))
    hdr == ["code", "d", "p", "p_L", "ci_low", "ci_high"] ||
        throw(ArgumentError("CSV header must be code,d,p,p_L,ci_low,ci_high (got $(join(hdr, ",")))"))
    rows = [(f = strip.(split(l, ",")); (f[1], parse(Int, f[2]), parse(Float64, f[3]), parse(Float64, f[4]),
             parse(Float64, f[5]), parse(Float64, f[6]))) for l in lines[2:end]]
    TablePL(rows; code)
end

function pl(m::TablePL, p::Real, d::Int)
    haskey(m.points, d) || throw(ArgumentError("TablePL($(m.code)) has no data for d = $d"))
    pts = m.points[d]
    length(pts) == 1 && pts[1][1] == p && return pts[1][2]
    (pts[1][1] <= p <= pts[end][1]) ||
        throw(ArgumentError("p = $p outside the tabulated range [$(pts[1][1]), $(pts[end][1])] for d = $d"))
    i = findlast(x -> x[1] <= p, pts)
    i === nothing && throw(ArgumentError("p = $p below the tabulated range"))
    i == length(pts) && return pts[end][2]
    (p1, y1), (p2, y2) = pts[i], pts[i + 1]
    t = (log(p) - log(p1)) / (log(p2) - log(p1))
    exp(log(y1) + t * (log(y2) - log(y1)))
end

# ---------------------------------------------------------------- code models

"""
    CodeModel

A QEC code block. Interface (implemented by every concrete model):
`code_n`, `code_k`, `code_d`, `syndrome_ancillas`, `cycle_time(m, hw)`,
`logical_error_per_cycle(m, p)`, `intra_measurement_cost(m, support)`,
`max_parallel_surgeries(m)`, `aux_graph_qubits(m, weight)`, `tanner_graph(m)`,
plus `distance_family(m)` and `with_distance(m, d)` for distance optimization.

The spec's `n`, `k`, `d` are named `code_n`, `code_k`, `code_d` here (the QECCore
convention), to avoid exporting one-letter names.
"""
abstract type CodeModel end

max_parallel_surgeries(::CodeModel) = 1
logical_error_per_cycle(m::CodeModel, p::Real) = pl(m.pl, p, code_d(m))
distance_family(::CodeModel) = nothing        # not optimizable unless a model says so

"""
    SurfaceCodeModel(d; pl = AnalyticPL(), ancilla_tiles_per_qubit = 1.0, max_parallel = 1)

Rotated surface code patch, `[[d², 1, d]]` with `d² − 1` syndrome qubits (one per
stabilizer), so a tile is `2d² − 1` physical qubits.

- Cycle: [Fowler] Fig. 1, eight steps: initialization, Hadamard, four CNOTs,
  Hadamard, measurement, so `cycle_time = 3·t_1q + 4·t_2q + t_meas`.
- `pl`: [GoSC] Eq. (10).
- Pauli product measurements take `d` code cycles of lattice surgery ([GoSC] Sec. 2:
  "d code cycles to account for measurement errors").
- `ancilla_tiles_per_qubit`: lattice-surgery routing tiles per logical qubit in the
  measured product. Default 1, a **modelling assumption**: one adjacent routing
  tile per patch, as in [GoSC]'s data-block layouts (the compact block uses 1.53
  tiles per stored qubit including its routing region). Set it for your layout.
"""
struct SurfaceCodeModel <: CodeModel
    d::Int
    pl::PLModel
    ancilla_tiles_per_qubit::Float64
    max_parallel::Int
    function SurfaceCodeModel(d::Int; pl::PLModel = AnalyticPL(), ancilla_tiles_per_qubit::Real = 1.0, max_parallel::Int = 1)
        (d >= 3 && isodd(d)) || throw(ArgumentError("surface code distance must be odd and ≥ 3 (got $d)"))
        new(d, pl, ancilla_tiles_per_qubit, max_parallel)
    end
end

code_n(m::SurfaceCodeModel) = m.d^2
code_k(::SurfaceCodeModel) = 1
code_d(m::SurfaceCodeModel) = m.d
syndrome_ancillas(m::SurfaceCodeModel) = m.d^2 - 1
cycle_time(::SurfaceCodeModel, hw::HardwareModel) = 3hw.t_1q + 4hw.t_2q + hw.t_meas
max_parallel_surgeries(m::SurfaceCodeModel) = m.max_parallel
_tile(m::SurfaceCodeModel) = 2m.d^2 - 1
function intra_measurement_cost(m::SurfaceCodeModel, support::Int)
    (ceil(Int, m.ancilla_tiles_per_qubit * support) * _tile(m), m.d)
end
aux_graph_qubits(m::SurfaceCodeModel, weight::Int) = ceil(Int, m.ancilla_tiles_per_qubit) * _tile(m)
distance_family(::SurfaceCodeModel) = :surface
with_distance(m::SurfaceCodeModel, d::Int) = SurfaceCodeModel(d; pl = m.pl, ancilla_tiles_per_qubit = m.ancilla_tiles_per_qubit,
                                                               max_parallel = m.max_parallel)
next_distance(::SurfaceCodeModel, d::Int) = d + 2
min_distance(::SurfaceCodeModel) = 3

"""
    BBCodeModel(name; intra_ancilla, intra_rounds, aux_graph_qubits, pl = preset fit, max_parallel = 1)

Bivariate bicycle code presets from [BB] Table 3:
- `:gross`: `[[144, 12, 12]]`, l = 12, m = 6, A = x³ + y + y², B = y³ + x + x²
- `:two_gross`: `[[288, 12, 18]]`, l = 12, m = 12, A = x³ + y² + y⁷, B = y³ + x + x²

From [BB]: `n` check qubits (2n physical qubits in total); a syndrome cycle with 7
CNOT layers plus initialization and measurement (`cycle_time = 7·t_2q + t_1q +
t_meas`, Sec. 5); and the `p_L` fit of Sec. 8 / Table 6 with `d_circ` the Table 1
upper bound (10 and 18).

Required (no published source for a generic value): `intra_ancilla` and
`intra_rounds`, the ancilla qubits and rounds of an in-block logical measurement
(for example from an auxiliary-graph construction), and `aux_graph_qubits`, a
function `weight -> qubits` for the auxiliary graph used when this block takes part
in an adapter join ([Adapt] Table II gives 35–47 qubits for weight 12–14 logicals
of a [[98,6,12]] code).
"""
struct BBCodeModel <: CodeModel
    name::Symbol
    n::Int
    k::Int
    d::Int
    l::Int
    m::Int
    A::NTuple{3,Tuple{Int,Int}}       # monomials (x-power, y-power)
    B::NTuple{3,Tuple{Int,Int}}
    pl::PLModel
    intra_ancilla::Int
    intra_rounds::Int
    aux_graph::Function
    max_parallel::Int
end

const BB_PRESETS = Dict(
    :gross => (n = 144, k = 12, d = 12, l = 12, m = 6, A = ((3, 0), (0, 1), (0, 2)), B = ((0, 3), (1, 0), (2, 0)),
               pl = BBFitPL(10, 18.04, 1337.0, -96007.0)),
    :two_gross => (n = 288, k = 12, d = 18, l = 12, m = 12, A = ((3, 0), (0, 2), (0, 7)), B = ((0, 3), (1, 0), (2, 0)),
                   pl = BBFitPL(18, 32.04, 3522.0, -294482.0)),
)

function BBCodeModel(name::Symbol; intra_ancilla::Int, intra_rounds::Int, aux_graph_qubits::Function,
                     pl::Union{Nothing,PLModel} = nothing, max_parallel::Int = 1)
    haskey(BB_PRESETS, name) || throw(ArgumentError("unknown BB preset $name (have $(collect(keys(BB_PRESETS))))"))
    p = BB_PRESETS[name]
    BBCodeModel(name, p.n, p.k, p.d, p.l, p.m, p.A, p.B, pl === nothing ? p.pl : pl, intra_ancilla, intra_rounds,
                aux_graph_qubits, max_parallel)
end

code_n(m::BBCodeModel) = m.n
code_k(m::BBCodeModel) = m.k
code_d(m::BBCodeModel) = m.d
syndrome_ancillas(m::BBCodeModel) = m.n
cycle_time(::BBCodeModel, hw::HardwareModel) = 7hw.t_2q + hw.t_1q + hw.t_meas
max_parallel_surgeries(m::BBCodeModel) = m.max_parallel
intra_measurement_cost(m::BBCodeModel, support::Int) = (m.intra_ancilla, m.intra_rounds)
aux_graph_qubits(m::BBCodeModel, weight::Int) = Int(m.aux_graph(weight))
distance_family(::BBCodeModel) = :bivariate_bicycle
const _BB_ORDER = (:gross, :two_gross)
function with_distance(m::BBCodeModel, d::Int)
    name = d == 12 ? :gross : d == 18 ? :two_gross : throw(ArgumentError("no BB preset with d = $d"))
    p = BB_PRESETS[name]
    pl = m.pl isa BBFitPL ? p.pl : m.pl
    BBCodeModel(name, p.n, p.k, p.d, p.l, p.m, p.A, p.B, pl, m.intra_ancilla, m.intra_rounds, m.aux_graph, m.max_parallel)
end
next_distance(::BBCodeModel, d::Int) = d == 12 ? 18 : throw(ArgumentError("no larger BB preset than d = $d"))
min_distance(::BBCodeModel) = 12

"""
    GenericCSSModel(code; d = nothing, compute_distance = false, pl, cycle_layers,
                    intra_ancilla, intra_rounds, aux_graph_qubits, syndrome_ancillas = nothing)

Wraps any QuantumClifford stabilizer code. `n` and `k` are read from the code, and
`syndrome_ancillas` defaults to the number of checks (one ancilla per check). `d`
must be given, unless `compute_distance = true`, which calls QuantumClifford's
`distance` (this may need its optional solver extensions). Required:
- the error model `pl`;
- `cycle_layers = (n_1q, n_2q, n_meas)`, the depth of one syndrome cycle in each
  kind of layer;
- the in-block measurement cost (`intra_ancilla`, `intra_rounds`);
- `aux_graph_qubits(weight)`.
"""
struct GenericCSSModel{C} <: CodeModel
    code::C
    n::Int
    k::Int
    d::Int
    nsyn::Int
    pl::PLModel
    cycle_layers::NTuple{3,Int}
    intra_ancilla::Int
    intra_rounds::Int
    aux_graph::Function
end

function GenericCSSModel(code; d::Union{Nothing,Int} = nothing, compute_distance::Bool = false, pl::PLModel,
                         cycle_layers::NTuple{3,Int}, intra_ancilla::Int, intra_rounds::Int, aux_graph_qubits::Function,
                         syndrome_ancillas::Union{Nothing,Int} = nothing)
    n = QuantumClifford.ECC.code_n(code)
    k = QuantumClifford.ECC.code_k(code)
    if d === nothing
        compute_distance || throw(ArgumentError("give d, or set compute_distance = true"))
        d = QuantumClifford.ECC.distance(code)
    end
    ns = syndrome_ancillas === nothing ? size(QuantumClifford.ECC.parity_checks(code), 1) : syndrome_ancillas
    GenericCSSModel(code, n, k, d, ns, pl, cycle_layers, intra_ancilla, intra_rounds, aux_graph_qubits)
end

code_n(m::GenericCSSModel) = m.n
code_k(m::GenericCSSModel) = m.k
code_d(m::GenericCSSModel) = m.d
syndrome_ancillas(m::GenericCSSModel) = m.nsyn
cycle_time(m::GenericCSSModel, hw::HardwareModel) = m.cycle_layers[1] * hw.t_1q + m.cycle_layers[2] * hw.t_2q + m.cycle_layers[3] * hw.t_meas
intra_measurement_cost(m::GenericCSSModel, support::Int) = (m.intra_ancilla, m.intra_rounds)
aux_graph_qubits(m::GenericCSSModel, weight::Int) = Int(m.aux_graph(weight))

Base.show(io::IO, m::SurfaceCodeModel) = print(io, "SurfaceCodeModel(d=", m.d, ")")
Base.show(io::IO, m::BBCodeModel) = print(io, "BBCodeModel(:", m.name, ", [[", m.n, ",", m.k, ",", m.d, "]])")
Base.show(io::IO, m::GenericCSSModel) = print(io, "GenericCSSModel([[", m.n, ",", m.k, ",", m.d, "]])")
