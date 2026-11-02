"""
    PBCCompiler.Verification

Checks that a compiled PBC program implements the same logical operation as its
input circuit, and finds counterexamples when it does not.

The entry point is [`verify_compilation`](@ref). The verifier deliberately
reuses none of the compiler's transformation code (only the `CircuitOp` data
types and QuantumClifford.jl primitives), so a bug in a compiler pass cannot hide
itself from the check.
"""
module Verification

using ..PBCCompiler: Circuit, CircuitOp
using Moshi.Match: @match
import QuantumClifford
import LinearAlgebra
using Random: AbstractRNG, randperm
using QuantumClifford: PauliOperator, CliffordOperator, nqubits, embed, comm

include("affine.jl")
include("ops.jl")
include("wellformed.jl")
include("frame.jl")
include("paulimap.jl")
include("lowering.jl")
include("normal_form.jl")
include("dense.jl")
include("random_circuits.jl")
include("cyclotomic.jl")
include("gadgets.jl")

export AffineBit, affine_var, isconstant, evaluate, substitute,
       check_wellformed, wellformedness_issues, WellformednessIssue,
       SymbolicPauli, mul_conditional!, mul_frame!, conjugate!, commutation_parity, frame_variables,
       PauliMap, identity_map, apply_map, rotate!, compose, to_clifford,
       VOp, LoweredProgram, lower, NFItem, NormalForm, normal_form, frame_clifford,
       random_clifford_t_circuit,
       GadgetProof, prove_reduced_gadget, synthesize_z1_clifford, LiftedGadget, GadgetFailure,
       LiftingResult, lift_gadgets, lower_lifted

end # module Verification
