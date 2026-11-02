"""
    PBCCompiler.ResourceEstimation

Physical resource estimates for a compiled PBC program: physical qubits (data,
syndrome, surgery/adapter ancilla, factories, routing), time, and logical failure
probability, given a [`LogicalLayout`](@ref) of logical qubits onto QEC code blocks,
a [`HardwareModel`](@ref), an adapter cost model and a magic-state factory model.
The entry point is [`estimate_resources`](@ref).
"""
module ResourceEstimation

import ..PBCCompiler
using ..PBCCompiler: Circuit, CircuitOp, affectedqubits
using ..GadgetShapes: count_gadgets, find_gadgets
using Moshi.Match: @match
import QuantumClifford
using QuantumClifford: PauliOperator, nqubits, embed, comm
import Graphs
import LinearAlgebra
import Random
using Random: AbstractRNG

include("code_models.jl")
include("checks.jl")
include("layout.jl")
include("adapter_cost.jl")
include("classify.jl")
include("hardware.jl")
include("scheduler.jl")
include("magic.jl")
include("estimate.jl")
include("report.jl")

export HardwareModel, PLModel, AnalyticPL, FitPL, BBFitPL, TablePL,
       CodeModel, SurfaceCodeModel, BBCodeModel, GenericCSSModel,
       code_n, code_k, code_d, syndrome_ancillas, cycle_time, logical_error_per_cycle,
       intra_measurement_cost, max_parallel_surgeries, aux_graph_qubits,
       parity_check_matrices, gf2_rank,
       LogicalLayout, validate_layout, pack_sequential, pack_by_interaction, interaction_graph,
       AdapterCostModel, AnalyticAdapterModel, AdapterCost, adapter_cost,
       EstOp, classify_operations, weight_histogram,
       heavy_hex, tanner_graph, EmbeddingResult, embed_tanner_graph, EmbeddingOverhead, embedding_overhead,
       ScheduleResult, dependency_dag, schedule_operations,
       FactoryModel, FifteenToOneFactory, TableFactory, litinski_factory, FactorySpec, factory_spec,
       magic_state_counts, choose_factory_count,
       ResourceEstimate, estimate_resources, to_dict

end # module ResourceEstimation
