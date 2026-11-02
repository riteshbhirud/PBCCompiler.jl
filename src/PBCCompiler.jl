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
include("resource_estimation/ResourceEstimation.jl")
using .ResourceEstimation: estimate_resources, LogicalLayout, HardwareModel
export Verification, ResourceEstimation, check_wellformed, estimate_resources, LogicalLayout, HardwareModel, verify_compilation, verify_passes, VerificationResult, VerificationStatus,
    EQUIVALENT, NOT_EQUIVALENT, UNKNOWN

end # module PBCCompiler
