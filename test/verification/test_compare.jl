@testitem "Verification V5: compiler test circuits verify EQUIVALENT via the normal form" tags=[:verification, :v5] begin
using PBCCompiler
using PBCCompiler: CircuitOp, parse_input, random_test_circuit
using QuantumClifford: @P_str
using Random
const C = CircuitOp

isnf(r) = r.status == EQUIVALENT && r.method == :normal_form && r.certificate.identity_map

# circuits from the existing compiler tests (test_compilation_pass.jl)
@test isnf(verify_compilation([C.ExpQuatPiPauli(P"Z", [1]), C.ExpQuatPiPauli(P"X", [1]), C.ExpQuatPiPauli(P"Z", [1]),
    C.PauliConditional(P"Z", [1], P"X", [2]), C.ExpEighPiPauli(P"Z", [2]), C.Measurement(P"Z", 1, [1]), C.Measurement(P"Z", 2, [2])]))
@test isnf(verify_compilation([C.PauliConditional(P"X", [1], P"Z", [2]), C.ExpHalfPiPauli(P"YZ", [1, 2]),
    C.ExpEighPiPauli(P"Z", [2]), C.Measurement(P"Z", 1, [1]), C.Measurement(P"Z", 2, [2])]))
# repo fixture (121 T gates) and QFT-2 through the QASM parser
@test isnf(verify_compilation(parse_input(joinpath(pkgdir(PBCCompiler), "test", "fixtures", "sample_input.qasm"))))
qasm = "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[2];\ncreg c[2];\nh q[0];\nt q[0];\nt q[1];\ncx q[0],q[1];\ntdg q[1];\ncx q[0],q[1];\nh q[1];\nmeasure q[0] -> c[1];\nmeasure q[1] -> c[0];\n"
path = tempname() * ".qasm"; write(path, qasm)
@test isnf(verify_compilation(parse_input(path)))
rm(path)
# the compiler's own random generator (test_random_circuit.jl)
Random.seed!(11)
for _ in 1:20
    @test isnf(verify_compilation(random_test_circuit(rand(4:12), rand(1:4))))
end
# random Clifford+T circuits that avoid the known compiler bugs B1/B2
rng = MersenneTwister(2027)
for _ in 1:60
    n = rand(rng, 1:6)
    c = PBCCompiler.Verification.random_clifford_t_circuit(rng, n; depth=rand(rng, 1:6n), sorted_qubits=true,
                                                            t_density=0.35, multiqubit_density=0.25)
    @test isnf(verify_compilation(c))
end
end

@testitem "Verification V5: hand-built equivalent pairs" tags=[:verification, :v5] begin
using PBCCompiler
using PBCCompiler: CircuitOp
using QuantumClifford: @P_str
const C = CircuitOp
v(a, b; kw...) = verify_compilation(a, b; oracle=:none, kw...)
M1 = C.Measurement(P"Z", 1, [1]); M2 = C.Measurement(P"Z", 2, [2])
T(p, q) = C.ExpEighPiPauli(p, q); S(p, q) = C.ExpQuatPiPauli(p, q); X(p, q) = C.ExpHalfPiPauli(p, q)

# reordering commuting non-Clifford rotations
@test v([T(P"ZZ", [1, 2]), T(P"Z", [1]), C.Measurement(P"X", 1, [1]), C.Measurement(P"X", 2, [2])],
        [T(P"Z", [1]), T(P"ZZ", [1, 2]), C.Measurement(P"X", 1, [1]), C.Measurement(P"X", 2, [2])]).status == EQUIVALENT
# merging T·T = S, and splitting T† = T·S†
@test v([T(P"X", [1]), T(P"X", [1]), M1], [S(P"X", [1]), M1]).status == EQUIVALENT
@test v([T(P"-X", [1]), M1], [T(P"X", [1]), S(P"-X", [1]), M1]).status == EQUIVALENT
# rotations about the same Pauli separated by commuting operations merge
@test v([T(P"X", [1]), S(P"Z", [2]), T(P"X", [1]), M1, M2], [S(P"X", [1]), S(P"Z", [2]), M1, M2]).status == EQUIVALENT
# pushing a Clifford through a measurement (what remove_clifford does)
@test v([S(P"X", [1]), M1], [C.Measurement(P"Y", 1, [1])]).status == EQUIVALENT
# outcome-flipping Pauli: equivalent only under the parity map bit1 = 1 ⊕ m
r = v([M1], [X(P"X", [1]), M1])
@test r.status == EQUIVALENT && !r.certificate.identity_map
@test v([M1], [X(P"X", [1]), M1]; outcome_map=:identity).status != EQUIVALENT
# explicit outcome map: input bit 1 = 1 ⊕ compiled bit 1
using PBCCompiler.Verification: AffineBit
@test v([M1], [X(P"X", [1]), M1]; outcome_map=Dict(1 => AffineBit(true, [1]))).status == EQUIVALENT
@test v([M1], [X(P"X", [1]), M1]; outcome_map=Dict(1 => AffineBit(false, [1]))).status != EQUIVALENT
# classical control: a conditional correction in both programs
cc_in = [C.Measurement(P"X", 1, [1]), C.BitConditional(C.ExpHalfPiPauli(P"Z", [1]), 1), T(P"Z", [2]), C.Measurement(P"X", 3, [1]), M2]
cc_out = [C.Measurement(P"X", 1, [1]), T(P"Z", [2]), C.BitConditional(C.ExpHalfPiPauli(P"Z", [1]), 1), C.Measurement(P"X", 3, [1]), M2]
@test v(cc_in, cc_out).status == EQUIVALENT
# channel semantics: identical programs are equivalent, a dropped trailing Clifford is not
@test v([S(P"X", [1]), M1, S(P"Z", [1])], [S(P"X", [1]), M1, S(P"Z", [1])]; semantics=:channel).status == EQUIVALENT
@test v([S(P"X", [1]), M1, S(P"Z", [1])], [S(P"X", [1]), M1]; semantics=:channel).status != EQUIVALENT
@test v([S(P"X", [1]), M1, S(P"Z", [1])], [S(P"X", [1]), M1]).status == EQUIVALENT     # outcomes agree
end

@testitem "Verification V5: never EQUIVALENT on inequivalent pairs (dense ground truth)" tags=[:verification, :v5] setup=[NFTestUtils] begin
using PBCCompiler
using PBCCompiler: CircuitOp
using PBCCompiler.Verification: random_clifford_t_circuit
using QuantumClifford: PauliOperator
using Random
using .NFTestUtils: circuit_povm, povm_equal
const C = CircuitOp

# a random single-op mutation of a circuit
function mutate(rng, c)
    c = copy(c)
    k = rand(rng, eachindex(c))
    op = c[k]
    kind = rand(rng, 1:4)
    if kind == 1
        deleteat!(c, k)
    elseif kind == 2 && k < length(c)
        c[k], c[k + 1] = c[k + 1], c[k]
    elseif kind == 3
        q = rand(rng, 1:PBCCompiler.Verification._width(c))
        insert!(c, k, C.ExpEighPiPauli(PauliOperator(0x00, [rand(rng, Bool)], [true]), [q]))
    else
        q = rand(rng, 1:PBCCompiler.Verification._width(c))
        insert!(c, k, C.ExpQuatPiPauli(PauliOperator(rand(rng, Bool) ? 0x02 : 0x00, [true], [rand(rng, Bool)]), [q]))
    end
    c
end

rng = MersenneTwister(99)
n_ineq = Ref(0); n_eq_proved = Ref(0)
for _ in 1:200
    n = rand(rng, 1:3)
    a = random_clifford_t_circuit(rng, n; depth=rand(rng, 1:3n+2), meas_density=0.15, cc_density=0.15, t_density=0.35)
    b = mutate(rng, a)
    PBCCompiler.Verification._width(b) == n || continue
    any(op -> PBCCompiler.Verification._opkind(op) == PBCCompiler.Verification.OP_MEASUREMENT, b) || continue
    isempty(check_wellformed(b)) || continue
    same = povm_equal(circuit_povm(a, n), circuit_povm(b, n))
    r = verify_compilation(a, b; oracle=:none, outcome_map=:identity)
    same || (n_ineq[] += 1)
    r.status == EQUIVALENT && (n_eq_proved[] += 1)
    @test !(r.status == EQUIVALENT && !same)          # soundness
end
@test n_ineq[] > 50 && n_eq_proved[] > 10              # both cases really occur
end

@testitem "Verification V5: known compiler bugs are never certified" tags=[:verification, :v5] begin
using PBCCompiler
using PBCCompiler: CircuitOp
using QuantumClifford: @P_str
const C = CircuitOp
# B1: CNOT with control index > target index
@test verify_compilation([C.PauliConditional(P"Z", [2], P"X", [1]), C.Measurement(P"Z", 1, [1]), C.Measurement(P"Z", 2, [2])]).status != EQUIVALENT
# B2: T gate on an unsorted qubit list
@test verify_compilation([C.ExpEighPiPauli(P"XZ", [2, 1]), C.Measurement(P"Z", 1, [1]), C.Measurement(P"Z", 2, [2])]).status != EQUIVALENT
# B3: gadget bit collides with a program bit (contract violation, never certified)
r = verify_compilation([C.Measurement(P"Z", 1, [1]), C.ExpEighPiPauli(P"X", [1]), C.Measurement(P"Z", 2, [1])])
@test r.status != EQUIVALENT && !isempty(r.wellformedness_errors)
# B4: no measurement, the compiler throws and the error is reported
r = verify_compilation([C.ExpQuatPiPauli(P"Z", [1]), C.ExpQuatPiPauli(P"X", [1])])
@test r.status == UNKNOWN && occursin("compiler threw", only(r.stats[:reasons]))
end
