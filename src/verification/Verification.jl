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
using QuantumClifford: PauliOperator, CliffordOperator, nqubits, embed, comm

include("affine.jl")
include("ops.jl")
include("wellformed.jl")
include("frame.jl")

export AffineBit, affine_var, isconstant, evaluate, substitute,
       check_wellformed, wellformedness_issues, WellformednessIssue,
       SymbolicPauli, mul_conditional!, mul_frame!, conjugate!, commutation_parity, frame_variables

end # module Verification
