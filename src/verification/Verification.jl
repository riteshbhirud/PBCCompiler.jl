"""
    PBCCompiler.Verification

Checks that a compiled PBC program implements the same logical operation as its
input circuit, and finds counterexamples when it does not.

The entry point is [`verify_compilation`](@ref).
"""
module Verification

import ..PBCCompiler
using ..PBCCompiler: Circuit, CircuitOp
using Moshi.Match: @match
import QuantumClifford
import LinearAlgebra
import Random
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
include("canonical.jl")
include("compare.jl")
include("types.jl")
include("oracle_dense.jl")
include("oracle_stabilizer.jl")
include("oracle_dispatch.jl")
include("counterexample.jl")
include("passwise.jl")
include("verify.jl")

export AffineBit, affine_var, isconstant, evaluate, substitute,
       check_wellformed, wellformedness_issues, WellformednessIssue,
       SymbolicPauli, mul_conditional!, mul_frame!, conjugate!, commutation_parity, frame_variables,
       PauliMap, identity_map, apply_map, rotate!, compose, to_clifford,
       VOp, LoweredProgram, lower, NFItem, NormalForm, normal_form, frame_clifford,
       random_clifford_t_circuit,
       GadgetProof, prove_reduced_gadget, synthesize_z1_clifford, LiftedGadget, GadgetFailure,
       LiftingResult, lift_gadgets, lower_lifted,
       CanonicalForm, canonicalize, ComparisonResult, compare_forms,
       VerificationStatus, EQUIVALENT, NOT_EQUIVALENT, UNKNOWN,
       Certificate, Counterexample, VerificationResult, verify_compilation,
       DenseVerdict, dense_oracle, StabilizerVerdict, stabilizer_oracle,
       minimize_counterexample, minimize_pair,
       PREPROCESS_PASSES, PassStep, run_passes, verify_passes, verify_passes_detailed

end # module Verification
