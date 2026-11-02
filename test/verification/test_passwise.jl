@testitem "Verification V8: per-pass translation validation" tags=[:verification, :v8] begin
using PBCCompiler
using PBCCompiler: CircuitOp, parse_input, preprocess_circuit, random_test_circuit
using PBCCompiler.Verification: verify_passes_detailed, run_passes, PREPROCESS_PASSES, random_clifford_t_circuit
using QuantumClifford: @P_str
using Random
const C = CircuitOp

first_bad(steps) = (i = findfirst(s -> s.result.status != EQUIVALENT, steps); i === nothing ? nothing : steps[i].pass)

# drift guard: the mirrored pass list reproduces preprocess_circuit exactly
rng = MersenneTwister(88)
for _ in 1:30
    c = random_clifford_t_circuit(rng, rand(rng, 1:5); depth=rand(rng, 2:12), meas_density=0.1, cc_density=0.1)
    trace = run_passes(c)
    @test [t[1] for t in trace] == collect(PREPROCESS_PASSES)
    @test trace[end][2] == preprocess_circuit(copy(c))
end
@test isempty(run_passes([C.Measurement(P"Z", 1, [1])]))      # preprocess_circuit's early return

# every pass is validated on the compiler's own test circuits and the fixture
for c in ([C.PauliConditional(P"X", [1], P"Z", [2]), C.ExpHalfPiPauli(P"YZ", [1, 2]), C.ExpEighPiPauli(P"Z", [2]),
           C.Measurement(P"Z", 1, [1]), C.Measurement(P"Z", 2, [2])],
          parse_input(joinpath(pkgdir(PBCCompiler), "test", "fixtures", "sample_input.qasm")))
    steps = verify_passes_detailed(c)
    @test length(steps) == 6 && first_bad(steps) === nothing
end
Random.seed!(3)
for _ in 1:10
    @test all(r -> r.status == EQUIVALENT, verify_passes(random_test_circuit(rand(4:10), rand(2:4))))
end

# bugs are localized to the pass that introduces them
b1 = [C.ExpQuatPiPauli(P"X", [1]), C.PauliConditional(P"Z", [2], P"X", [1]), C.Measurement(P"Z", 1, [1]), C.Measurement(P"Z", 2, [2])]
steps = verify_passes_detailed(b1)
@test first_bad(steps) == :remove_pauliconditional
@test steps[1].result.status == NOT_EQUIVALENT
b2 = [C.ExpEighPiPauli(P"XZ", [2, 1]), C.Measurement(P"Z", 1, [1]), C.Measurement(P"Z", 2, [2])]
steps = verify_passes_detailed(b2)
@test first_bad(steps) == :remove_nonclifford
@test steps[5].result.status == NOT_EQUIVALENT
# a compiler exception is reported, not thrown
steps = verify_passes_detailed([C.ExpQuatPiPauli(P"Z", [1]), C.ExpQuatPiPauli(P"X", [1])])
@test steps[end].pass == :error && steps[end].result.status == UNKNOWN
end
