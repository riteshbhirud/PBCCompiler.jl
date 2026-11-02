module PBCCompiler

include("type.jl")
include("io.jl")
include("traversal.jl")
include("affectedqubits.jl")
include("plotting.jl")
include("pair_transformation.jl")
include("preprocess.jl")
include("random_circuit.jl")
include("joint_measurement_check.jl")
include("logic.jl")

include("gadget_shapes.jl")
include("verification/Verification.jl")
using .Verification: check_wellformed, verify_compilation, verify_passes, VerificationResult, VerificationStatus,
    EQUIVALENT, NOT_EQUIVALENT, UNKNOWN
export Verification, check_wellformed, verify_compilation, verify_passes, VerificationResult, VerificationStatus,
    EQUIVALENT, NOT_EQUIVALENT, UNKNOWN

end # module PBCCompiler
