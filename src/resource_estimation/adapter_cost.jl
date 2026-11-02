"""Cost model for a Pauli measurement spanning several code blocks."""
abstract type AdapterCostModel end

"""
    AnalyticAdapterModel(; cycle_time_factor, aux_qubits_per_block = aux_graph_qubits,
                           adapter_qubits_per_join = d -> 2d - 1, rounds = d -> d)

Closed-form cost of a joint measurement on blocks `B_1, ..., B_a` with per-block
logical weights `w_i`:

- `d_joint = min_i d(B_i)`
- ancilla `= sum_i aux_qubits_per_block(B_i, w_i) + (a - 1) * adapter_qubits_per_join(d_joint)`
- rounds `= rounds(d_joint)`
- seconds per round `= cycle_time_factor * max_i cycle_time(B_i)`

Keyword arguments:

- `cycle_time_factor`: required, since no published value gives the slowdown of the merged code's syndrome cycle.
- `aux_qubits_per_block`: per-block auxiliary qubits; defaults to the block's `aux_graph_qubits`.
- `adapter_qubits_per_join`: default `2d - 1` (arXiv:2410.03628, Table III).
- `rounds`: default `d` (Litinski 2019, Sec. 2; arXiv:2410.03628, Sec. I).
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
round, the joint distance, a `source` (`:analytic` or `:constructed`), a note
(for example a recorded fallback), and `alt_ancilla`, the analytic model's
ancilla count for comparison (-1 if not computed).
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

"""
    adapters_available() -> Bool

Whether the loaded QuantumClifford provides `QuantumClifford.ECC.build_adapter`
(currently only in the `riteshbhirud/QuantumClifford.jl` fork, PR #742).
"""
adapters_available() = isdefined(QuantumClifford.ECC, :build_adapter)

"""
    ConstructedAdapterModel(; cycle_time_factor, fallback = AnalyticAdapterModel(; cycle_time_factor))

Adapter costs from an explicit construction of the merged code with the universal
adapters of Swaroop, Jochym-O'Connor and Yoder (arXiv:2410.03628), via
`QuantumClifford.ECC.build_adapter`.

Applies to measurements on exactly two blocks whose part on each block is a single
X- or Z-type logical. The ancilla count is the number of new qubits plus new checks,
rounds are `d_joint`, and seconds per round are `cycle_time_factor` times the slowest
block's cycle. The note records the merged code's maximum check weights. All other
measurements, or a failed construction, use `fallback`, and the fallback is recorded
in the estimate's notes. `d_joint` is a bound; the merged-code distance is not computed.
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

"""
    block_code(m::CodeModel)

The block's code as a QuantumClifford CSS code.
"""
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
        parts = Tuple{Int,Char}[]
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
            (-1, -1, sprint(showerror, err))
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
            # treat the X logical as a Z logical of the code with X and Z checks swapped
            p = getfield(ECC, :logx_ops)(code)[idx]
            push!(supports, [j for j in 1:(nqubits(p)::Int) if p[j][1]])
            code = ECC.CSS(Matrix{Bool}(ECC.parity_matrix_z(code)), Matrix{Bool}(ECC.parity_matrix_x(code)))
        else
            p = getfield(ECC, :logz_ops)(code)[idx]
            push!(supports, [j for j in 1:(nqubits(p)::Int) if p[j][2]])
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
