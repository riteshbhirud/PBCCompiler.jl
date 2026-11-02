# E4: adapter (cross-block joint measurement) cost models.

"""Cost model for a Pauli measurement spanning several code blocks."""
abstract type AdapterCostModel end

"""
    AnalyticAdapterModel(; cycle_time_factor, aux_qubits_per_block = aux_graph_qubits,
                           adapter_qubits_per_join = d -> 2d - 1, rounds = d -> d)

Closed-form cost of a joint measurement touching blocks `B₁ … B_a` with per-block
logical weights `wᵢ`:

- `d_joint = minᵢ d(Bᵢ)`;
- ancilla qubits `= Σᵢ aux_qubits_per_block(Bᵢ, wᵢ) + (a − 1)·adapter_qubits_per_join(d_joint)`;
- rounds `= rounds(d_joint)`;
- time per round `= cycle_time_factor · maxᵢ cycle_time(Bᵢ)`.

Fields and the source of each default:

- `aux_qubits_per_block`: the auxiliary-graph (or lattice-surgery routing) qubits each
  block needs. The default delegates to the block's own `aux_graph_qubits`, which
  is a required argument for qLDPC models and one routing tile for surface codes.
- `adapter_qubits_per_join`: default `2d − 1`. Swaroop, Jochym-O'Connor, Yoder,
  arXiv:2410.03628: an adapter of width `w = d` uses `w` qubits and `w − 1` checks
  (Table III: "27 (for adapter)" at `w = 14`). The spec suggested `3d`; the sourced
  value is used instead.
- `rounds`: default `d`. Lattice surgery and auxiliary-graph surgery need `O(d)`
  rounds to handle measurement errors (Litinski, GoSC, Sec. 2; arXiv:2410.03628
  Sec. I).
- `cycle_time_factor`: **required**. Merged codes have higher-weight checks (up to
  8 vs 6 in arXiv:2410.03628 Table III), so their syndrome cycle is longer, but no
  published value converts that into time.
"""
struct AnalyticAdapterModel{F1,F2,F3} <: AdapterCostModel
    aux_qubits_per_block::F1
    adapter_qubits_per_join::F2
    rounds::F3
    cycle_time_factor::Float64
end

function AnalyticAdapterModel(; cycle_time_factor::Real, aux_qubits_per_block = aux_graph_qubits,
                              adapter_qubits_per_join = d -> 2d - 1, rounds = d -> d)
    cycle_time_factor > 0 || throw(ArgumentError("cycle_time_factor must be positive"))
    AnalyticAdapterModel(aux_qubits_per_block, adapter_qubits_per_join, rounds, Float64(cycle_time_factor))
end

"""
    AdapterCost

Cost of one cross-block joint measurement: ancilla qubits, rounds, seconds per
round, the joint distance, a `source` (`:analytic` or `:constructed`), an optional
note (for example a recorded fallback), and `alt_ancilla`, the other model's
ancilla count for cross-checking (−1 if not computed).
"""
struct AdapterCost
    ancilla::Int
    rounds::Int
    round_time::Float64
    d_joint::Int
    source::Symbol
    note::String
    alt_ancilla::Int
end

"""
    adapter_cost(model, blocks::Vector{CodeModel}, weights::Vector{Int}, hw; kwargs...) -> AdapterCost

Cost of a joint measurement of logical weight `weights[i]` on block `blocks[i]`.
"""
function adapter_cost(am::AnalyticAdapterModel, blocks::AbstractVector{<:CodeModel}, weights::AbstractVector{Int},
                      hw::HardwareModel; kwargs...)
    a = length(blocks)
    a >= 2 || throw(ArgumentError("an adapter joins at least two blocks"))
    dj = minimum(code_d, blocks)
    anc = sum(Int(am.aux_qubits_per_block(b, w)) for (b, w) in zip(blocks, weights)) +
          (a - 1) * Int(am.adapter_qubits_per_join(dj))
    AdapterCost(anc, Int(am.rounds(dj)), am.cycle_time_factor * maximum(b -> cycle_time(b, hw), blocks), dj, :analytic, "", -1)
end

# ---------------------------------------------------------------- E5: constructed adapters

"""
    adapters_available() -> Bool

Whether the loaded QuantumClifford provides the universal-adapter construction
(`QuantumClifford.ECC.build_adapter`, from the `riteshbhirud/QuantumClifford.jl`
fork, PR #742).
"""
adapters_available() = isdefined(QuantumClifford.ECC, :build_adapter)

"""
    ConstructedAdapterModel(; cycle_time_factor, fallback = AnalyticAdapterModel(; cycle_time_factor))

Exact adapter costs from the universal-adapter construction of Swaroop,
Jochym-O'Connor and Yoder (arXiv:2410.03628), as implemented by
`QuantumClifford.ECC.build_adapter` (fork `riteshbhirud/QuantumClifford.jl`,
`master`, PR #742).

For a measurement on exactly two blocks, whose part on each block is a single X- or
Z-type logical, the merged code is built from the blocks' codes and the logical
representatives `logz_ops(code)[i]` (or `logx_ops(code)[i]` for X, used as a Z-type
logical of the code with its X and Z checks swapped), where `i` is the logical index from the layout. The costs are read
off the construction:
- ancilla = new qubits + new checks (one syndrome ancilla per check);
- rounds = `d_joint`;
- round time = `cycle_time_factor` × the slowest block's cycle.
The merged code's maximum check weights are recorded in the note, and the analytic
model's ancilla count is kept alongside for cross-checking. Constructions are
cached per (block code, logical index, letter) pair.

Anything else (more than two blocks, Y-type parts, several logicals of one block,
the fork not loaded, or a construction that throws, for example because the fork's
inter-code construction needs logical supports of equal length) uses `fallback`, and the fallback is **recorded** in the
estimate's notes. The merged-code distance is not computed (an ILP); `d_joint` is
reported as a bound.
"""
struct ConstructedAdapterModel <: AdapterCostModel
    cycle_time_factor::Float64
    fallback::AnalyticAdapterModel
    cache::Dict{Any,Tuple{Int,Int,String}}
end

function ConstructedAdapterModel(; cycle_time_factor::Real, fallback::Union{Nothing,AnalyticAdapterModel} = nothing)
    fb = fallback === nothing ? AnalyticAdapterModel(; cycle_time_factor) : fallback
    ConstructedAdapterModel(Float64(cycle_time_factor), fb, Dict{Any,Tuple{Int,Int,String}}())
end

"""The block's code as a QuantumClifford CSS code."""
block_code(m::GenericCSSModel) = m.code
function block_code(m::CodeModel)
    Hx, Hz = parity_check_matrices(m)
    QuantumClifford.ECC.CSS(Matrix{Bool}(Hx), Matrix{Bool}(Hz))
end

_model_key(m::CodeModel) = (typeof(m), code_n(m), code_k(m), code_d(m), m isa BBCodeModel ? m.name : :none)

function adapter_cost(cm::ConstructedAdapterModel, blocks::AbstractVector{<:CodeModel}, weights::AbstractVector{Int},
                      hw::HardwareModel; layout = nothing, qubits = Int[], pauli = nothing, n_data = 0, kwargs...)
    fb = adapter_cost(cm.fallback, blocks, weights, hw)
    reason = if !adapters_available()
        "the universal-adapter construction is not loaded (needs the riteshbhirud/QuantumClifford.jl fork)"
    elseif length(blocks) != 2
        "constructed adapters join exactly two blocks ($(length(blocks)) here)"
    elseif any(!=(1), weights)
        "a block contributes a product of several logical qubits"
    elseif layout === nothing || pauli === nothing
        "missing layout or Pauli"
    else
        ""
    end
    if isempty(reason)
        parts = Tuple{Int,Char}[]          # (logical index, letter) per block
        for b in sort!(unique([layout.qubit_to_block[q] for q in qubits if q <= n_data]))
            q = only(qq for qq in qubits if qq <= n_data && layout.qubit_to_block[qq] == b)
            x, z = pauli[q]
            push!(parts, (layout.qubit_within_block[q], x ? (z ? 'Y' : 'X') : 'Z'))
        end
        any(p -> p[2] == 'Y', parts) && (reason = "Y-type logical parts are not supported by the construction")
    end
    if !isempty(reason)
        return AdapterCost(fb.ancilla, fb.rounds, fb.round_time, fb.d_joint, :analytic, "constructed adapter fallback: " * reason, -1)
    end
    key = (_model_key(blocks[1]), parts[1], _model_key(blocks[2]), parts[2])
    anc, nq, note = get!(cm.cache, key) do
        try
            _construct_adapter(blocks, parts)
        catch err
            err isa InterruptException && rethrow()
            (-1, -1, sprint(showerror, err))          # cached failure: never retried
        end
    end
    if anc < 0
        return AdapterCost(fb.ancilla, fb.rounds, fb.round_time, fb.d_joint, :analytic,
                           "constructed adapter fallback: the construction failed: " * note, -1)
    end
    dj = minimum(code_d, blocks)
    AdapterCost(anc, dj, cm.cycle_time_factor * maximum(b -> cycle_time(b, hw), blocks), dj, :constructed,
                note * "; analytic model: $(fb.ancilla) ancilla", fb.ancilla)
end

function _construct_adapter(blocks, parts)
    ECC = QuantumClifford.ECC
    codes = Any[]; supports = Vector{Int}[]
    for ((idx, letter), m) in zip(parts, blocks)
        code = block_code(m)
        if letter == 'X'
            # logical i's X representative (paired with logz_ops[i]), as a Z-type operator
            # of the Hadamard-conjugated code (X and Z checks swapped)
            p = getfield(ECC, :logx_ops)(code)[idx]
            push!(supports, [j for j in 1:nqubits(p) if p[j][1]])
            code = ECC.CSS(Matrix{Bool}(ECC.parity_matrix_z(code)), Matrix{Bool}(ECC.parity_matrix_x(code)))
        else
            p = getfield(ECC, :logz_ops)(code)[idx]
            push!(supports, [j for j in 1:nqubits(p) if p[j][2]])
        end
        push!(codes, code)
    end
    a = getfield(ECC, :build_adapter)(codes[1], codes[2], supports[1], supports[2])
    rows(c) = size(ECC.parity_matrix_x(c), 1) + size(ECC.parity_matrix_z(c), 1)
    new_qubits = code_n(a) - code_n(codes[1]) - code_n(codes[2])
    new_checks = rows(a) - rows(codes[1]) - rows(codes[2])
    kexp = code_k(codes[1]) + code_k(codes[2]) - 1
    code_k(a) == kexp || error("internal: merged code has k = $(code_k(a)), expected $kexp")
    wx = maximum(sum(ECC.parity_matrix_x(a); dims = 2)); wz = maximum(sum(ECC.parity_matrix_z(a); dims = 2))
    (new_qubits + new_checks, new_qubits, "merged code: +$new_qubits qubits, +$new_checks checks, max check weights X=$wx Z=$wz")
end
