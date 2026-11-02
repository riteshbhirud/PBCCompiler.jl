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

include("verification/Verification.jl")
using .Verification: check_wellformed
export Verification, check_wellformed

end # module PBCCompiler
