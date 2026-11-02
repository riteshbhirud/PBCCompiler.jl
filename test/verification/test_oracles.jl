@testmodule OracleTestUtils begin
    using PBCCompiler
    using PBCCompiler: CircuitOp
    using QuantumClifford: PauliOperator, nqubits
    using Random
    const V = PBCCompiler.Verification

    # random mutation of a compiled program, returns (mutated, kind)
    function mutate_compiled(rng, c)
        c = copy(c)
        idx = [i for (i, op) in enumerate(c)]
        i = rand(rng, idx)
        op = c[i]
        k = V._opkind(op)
        kind = rand(rng, (:sign, :letter, :delete, :swap))
        if kind == :delete || !(k == V.OP_MEASUREMENT || V._isrotation(k) || k == V.OP_BITCOND)
            deleteat!(c, i); return (c, :delete)
        elseif kind == :swap && i < length(c)
            c[i], c[i + 1] = c[i + 1], c[i]; return (c, :swap)
        end
        inner, b = k == V.OP_BITCOND ? V._bitcond_parts(op) : (op, 0)
        ik = V._opkind(inner)
        (ik == V.OP_MEASUREMENT || V._isrotation(ik)) || (deleteat!(c, i); return (c, :delete))
        p, q = V._pauli_and_qubits(inner)
        p = copy(p)
        if kind == :sign
            p.phase[] = (p.phase[] + 0x02) & 0x03
        else
            j = rand(rng, 1:nqubits(p)); x, z = p[j]
            p[j] = (x ⊻ rand(rng, Bool), z ⊻ true)
        end
        newinner = ik == V.OP_MEASUREMENT ? CircuitOp.Measurement(p, V._measurement_bit(inner), q) :
                   ik == V.OP_EIGHPI ? CircuitOp.ExpEighPiPauli(p, q) :
                   ik == V.OP_QUATPI ? CircuitOp.ExpQuatPiPauli(p, q) :
                   ik == V.OP_HALFPI ? CircuitOp.ExpHalfPiPauli(p, q) : CircuitOp.Pauli(p, q)
        c[i] = k == V.OP_BITCOND ? CircuitOp.BitConditional(newinner, b) : newinner
        (c, kind)
    end
end

@testitem "Verification V6: normal form never contradicts the dense oracle" tags=[:verification, :v6] setup=[OracleTestUtils] begin
using PBCCompiler
using PBCCompiler.Verification: random_clifford_t_circuit, dense_oracle, _width
using Random
using .OracleTestUtils: mutate_compiled

rng = MersenneTwister(6060)
ncases = Ref(0); nnf_eq = Ref(0); ndense_neq = Ref(0)
for trial in 1:300
    n = rand(rng, 1:4)
    a = random_clifford_t_circuit(rng, n; depth=rand(rng, 1:3n+2), t_density=0.2, multiqubit_density=0.2,
                                  meas_density=0.1, cc_density=0.1, sorted_qubits=rand(rng, Bool))
    b = try PBCCompiler.preprocess_circuit(copy(a)) catch; continue end
    _width(b) <= 8 || continue
    rand(rng) < 0.5 && (b = first(mutate_compiled(rng, b)))
    nf = verify_compilation(a, b; oracle=:none)
    dv = dense_oracle(a, b; n_data=n, outcome_map=:infer, max_qubits=8)
    dv.decided || continue
    ncases[] += 1
    nf.status == EQUIVALENT && (nnf_eq[] += 1)
    dv.agree || (ndense_neq[] += 1)
    # the normal form never claims equivalence that the dense ground truth refutes
    @test !(nf.status == EQUIVALENT && !dv.agree)
    if nf.status == EQUIVALENT && nf.certificate.identity_map
        @test dense_oracle(a, b; n_data=n, outcome_map=:identity, max_qubits=8).agree
    end
end
@test ncases[] >= 200 && nnf_eq[] >= 50 && ndense_neq[] >= 30
end

@testitem "Verification V6: dense oracle refutes and certifies nothing by itself" tags=[:verification, :v6] begin
using PBCCompiler
using PBCCompiler: CircuitOp
using PBCCompiler.Verification: dense_oracle
using QuantumClifford: @P_str
const C = CircuitOp
# equivalent programs: the oracle agrees, but verify_compilation never turns that into EQUIVALENT
a = [C.Measurement(P"Z", 1, [1]), C.Measurement(P"Z", 2, [1])]
b = [C.Measurement(P"Z", 1, [1]), C.Measurement(P"Z", 2, [1])]
@test dense_oracle(a, b).agree
# equivalent only under a parity map: the normal form cannot find it, the dense oracle can
a2 = [C.Measurement(P"Z", 1, [1]), C.Measurement(P"Z", 2, [2])]
b2 = [C.Measurement(P"Z", 1, [1]), C.Measurement(P"ZZ", 2, [1, 2])]
r = verify_compilation(a2, b2)
@test r.status == UNKNOWN && r.method == :dense_oracle
@test dense_oracle(a2, b2; outcome_map=:infer).agree
@test !dense_oracle(a2, b2; outcome_map=:identity).agree
# with the identity map the same pair is refuted
r = verify_compilation(a2, b2; outcome_map=:identity)
@test r.status == NOT_EQUIVALENT && r.counterexample.confirmed
@test abs(r.counterexample.prob_input - r.counterexample.prob_compiled) > 0.1
# channel semantics: a dropped trailing Clifford changes the post-measurement state
c1 = [C.Measurement(P"Z", 1, [1]), C.ExpQuatPiPauli(P"X", [1])]
c2 = [C.Measurement(P"Z", 1, [1])]
@test dense_oracle(c1, c2).agree
@test !dense_oracle(c1, c2; semantics=:channel).agree
@test verify_compilation(c1, c2; semantics=:channel).status == NOT_EQUIVALENT
end

@testitem "Verification V6: stabilizer oracle finds injected bugs, never refutes correct output" tags=[:verification, :v6] setup=[OracleTestUtils] begin
using PBCCompiler
using PBCCompiler: CircuitOp
using PBCCompiler.Verification: random_clifford_t_circuit, stabilizer_oracle, dense_oracle, _width
using QuantumClifford: @P_str
using Random
using .OracleTestUtils: mutate_compiled

rng = MersenneTwister(7070)
false_pos = Ref(0); found = Ref(0); changing = Ref(0)
for trial in 1:40
    n = rand(rng, 3:5)
    a = random_clifford_t_circuit(rng, n; depth=3n, t_density=0.1, multiqubit_density=0.15, sorted_qubits=true)
    b = PBCCompiler.preprocess_circuit(copy(a))
    # correct compilation: must never be refuted
    sv = stabilizer_oracle(a, b; n_data=n, trials=5, rng)
    sv.applicable && sv.refuted && (false_pos[] += 1)
    # injected bug, with dense ground truth when small enough
    bm, kind = mutate_compiled(rng, b)
    _width(bm) <= 9 || continue
    dv = dense_oracle(a, bm; n_data=n, outcome_map=:identity, max_qubits=9)
    (dv.decided && !dv.agree) || continue          # only channel-changing mutations count
    changing[] += 1
    svm = stabilizer_oracle(a, bm; n_data=n, trials=20, rng)
    svm.applicable && svm.refuted && (found[] += 1)
end
@test false_pos[] == 0
@test changing[] >= 10
@test found[] >= 0.8 * changing[]
# beyond the dense limit: bug B1 (reversed CNOT) on 12 qubits
n = 12
a = vcat([CircuitOp.ExpQuatPiPauli(P"Z", [q]) for q in 1:n], [CircuitOp.ExpQuatPiPauli(P"X", [q]) for q in 1:n],
         [CircuitOp.ExpEighPiPauli(P"Z", [3]), CircuitOp.PauliConditional(P"Z", [7], P"X", [2])],
         [CircuitOp.Measurement(P"Z", q, [q]) for q in 1:n])
r = verify_compilation(a; max_dense_qubits=10, stabilizer_trials=30, rng=MersenneTwister(1))
@test r.status == NOT_EQUIVALENT && r.method == :stabilizer_oracle
@test r.counterexample.input_stabilizer !== nothing
end
