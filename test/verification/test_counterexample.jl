@testitem "Verification V7: minimized, confirmed counterexamples for compiler bugs" tags=[:verification, :v7] begin
using PBCCompiler
using PBCCompiler: CircuitOp
using PBCCompiler.Verification: random_clifford_t_circuit, dense_oracle, _width
using QuantumClifford: @P_str
using Random
const C = CircuitOp

# Bug B1 buried in a larger random circuit: the minimized input isolates the reversed CNOT.
rng = MersenneTwister(808)
for trial in 1:5
    n = 4
    pre = random_clifford_t_circuit(rng, n; depth=10, t_density=0.2, sorted_qubits=true, final_measure=false)
    a = vcat(pre, [C.PauliConditional(P"Z", [4], P"X", [1])], [C.Measurement(P"Z", q, [q]) for q in 1:n])
    r = verify_compilation(a)
    @test r.status == NOT_EQUIVALENT
    cx = r.counterexample
    @test cx.confirmed
    @test cx.minimized_input !== nothing && length(cx.minimized_input) < length(a)
    @test any(op -> PBCCompiler.Verification._opkind(op) == PBCCompiler.Verification.OP_PAULICOND, cx.minimized_input)
    # the minimized pair itself is still inequivalent (independent dense re-check)
    mi, mc = cx.minimized_input, cx.minimized_compiled
    @test !dense_oracle(mi, mc; n_data=_width(mi)).agree
    # the counterexample carries a concrete distinguishing input state and record
    @test length(cx.input_state) == 1 << n
    @test abs(cx.prob_input - cx.prob_compiled) > 1e-6
end

# Bug B2 (unsorted T) only bites when the T is not moved past an earlier Clifford
# (moving it re-embeds the Pauli correctly), so a preceding Clifford hides it:
hidden = [C.ExpQuatPiPauli(P"X", [3]), C.ExpEighPiPauli(P"XZ", [2, 1]), C.ExpQuatPiPauli(P"Z", [3]),
          C.Measurement(P"Z", 1, [1]), C.Measurement(P"Z", 2, [2]), C.Measurement(P"Z", 3, [3])]
@test verify_compilation(hidden).status == EQUIVALENT
# with the T first it is exposed, and the counterexample is minimized to the T gate plus measurements
a = [C.ExpEighPiPauli(P"XZ", [2, 1]), C.ExpQuatPiPauli(P"Z", [3]),
     C.Measurement(P"Z", 1, [1]), C.Measurement(P"Z", 2, [2]), C.Measurement(P"Z", 3, [3])]
r = verify_compilation(a)
@test r.status == NOT_EQUIVALENT && r.counterexample.confirmed
@test count(op -> PBCCompiler.Verification._opkind(op) == PBCCompiler.Verification.OP_EIGHPI, r.counterexample.minimized_input) == 1
@test length(r.counterexample.minimized_input) <= 4
end

@testitem "Verification V7: pair minimization for mutated compiled outputs" tags=[:verification, :v7] begin
using PBCCompiler
using PBCCompiler: CircuitOp
using PBCCompiler.Verification: dense_oracle, _width, random_clifford_t_circuit
using Random
include(joinpath(pkgdir(PBCCompiler), "test", "verification", "mutations.jl"))

rng = MersenneTwister(4242)
nfound = Ref(0)
for trial in 1:30
    n = rand(rng, 2:4)
    a = random_clifford_t_circuit(rng, n; depth=3n, t_density=0.25, sorted_qubits=true)
    b = PBCCompiler.preprocess_circuit(copy(a))
    bm, ok = mutate(rng, b, rand(rng, (:flip_sign, :change_letter, :drop_correction)); n_data=n)
    ok || continue
    r = verify_compilation(a, bm; outcome_map=:identity)
    r.status == NOT_EQUIVALENT || continue
    nfound[] += 1
    cx = r.counterexample
    @test cx.confirmed
    mi, mc = cx.minimized_input, cx.minimized_compiled
    @test mi !== nothing
    @test count(op -> PBCCompiler.Verification._opkind(op) == PBCCompiler.Verification.OP_MEASUREMENT, mi) <=
          count(op -> PBCCompiler.Verification._opkind(op) == PBCCompiler.Verification.OP_MEASUREMENT, a)
    @test !dense_oracle(mi, mc; n_data=n, outcome_map=:identity).agree
end
@test nfound[] >= 10
end

@testitem "Verification V7: divergence and stabilizer witness" tags=[:verification, :v7] begin
using PBCCompiler
using PBCCompiler: CircuitOp
using QuantumClifford: @P_str
using Random
const C = CircuitOp
# the first diverging canonical item is reported
r = verify_compilation([C.ExpEighPiPauli(P"Z", [1]), C.Measurement(P"X", 1, [1])],
                       [C.ExpEighPiPauli(P"-Z", [1]), C.Measurement(P"X", 1, [1])]; outcome_map=:identity)
@test r.status == NOT_EQUIVALENT
# T(−Z) normalizes to T(Z)·S†, so the rotations agree and the measurements (X vs ±Y) diverge
k, xi, yi = r.counterexample.divergence
@test k >= 1 && xi !== nothing && yi !== nothing && xi != yi
# beyond the dense limit the witness is a stabilizer state; marked unconfirmed unless minimization
# makes it small enough for the dense oracle
n = 12
a = vcat([C.ExpQuatPiPauli(P"X", [q]) for q in 1:n], [C.ExpEighPiPauli(P"Z", [5]), C.PauliConditional(P"Z", [9], P"X", [3])],
         [C.Measurement(P"Z", q, [q]) for q in 1:n])
r = verify_compilation(a; max_dense_qubits=10, rng=MersenneTwister(2))
@test r.status == NOT_EQUIVALENT && r.counterexample.method == :stabilizer_oracle
@test r.counterexample.input_stabilizer !== nothing
@test r.counterexample.confirmed == (PBCCompiler.Verification._width(r.counterexample.minimized_compiled) <= 10)
end
