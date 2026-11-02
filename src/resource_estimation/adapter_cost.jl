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
round, the joint distance, a `source` (`:analytic` or `:constructed`), and an
optional note (for example a recorded fallback).
"""
struct AdapterCost
    ancilla::Int
    rounds::Int
    round_time::Float64
    d_joint::Int
    source::Symbol
    note::String
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
    AdapterCost(anc, Int(am.rounds(dj)), am.cycle_time_factor * maximum(b -> cycle_time(b, hw), blocks), dj, :analytic, "")
end
