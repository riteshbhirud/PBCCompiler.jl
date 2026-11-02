@testitem "Verification V4: exact arithmetic in Q(ζ8)" tags=[:verification, :v4] begin
using PBCCompiler.Verification: Cyc8, zeta, INV_SQRT2, fieldnorm, abs2c, zeta_exponent, galois
using Random
@test zeta(1) * zeta(1) == zeta(2)
@test zeta(4) == Cyc8(-1) && zeta(8) == Cyc8(1) && zeta(-1) == zeta(7)
@test INV_SQRT2 * INV_SQRT2 == Cyc8(1 // 2)
@test conj(zeta(3)) == zeta(5)
@test zeta_exponent(zeta(6)) == 6 && zeta_exponent(Cyc8(2)) === nothing
rng = MersenneTwister(1)
for _ in 1:30
    a = Cyc8(rand(rng, -5:5) // rand(rng, 1:4), rand(rng, -5:5), rand(rng, -5:5) // 3, rand(rng, -5:5))
    iszero(a) && continue
    @test a * inv(a) == Cyc8(1)
    @test conj(abs2c(a)) == abs2c(a)                          # |a|² is real
    @test fieldnorm(a) == fieldnorm(galois(a, 3))             # the norm is Galois invariant
end
end

@testitem "Verification V4: Clifford synthesis V·P·V† = Z₁ (dense)" tags=[:verification, :v4] begin
using PBCCompiler.Verification: synthesize_z1_clifford, DensePauli, apply_rotation!, dense_matrix
using QuantumClifford: random_pauli, PauliOperator
using LinearAlgebra, Random
rng = MersenneTwister(9)
for _ in 1:40
    n = rand(rng, 1:5)
    p = random_pauli(rng, n; nophase=false, realphase=true)
    all(iszero, p.xz) && continue
    rots, _ = synthesize_z1_clifford(p)
    U = Matrix{ComplexF64}(I, 1 << n, 1 << n)
    for (q, k) in rots
        apply_rotation!(U, DensePauli(q), k // 4)
    end
    z1 = zero(PauliOperator, n); z1[1] = (false, true)
    @test U * dense_matrix(p) * U' ≈ dense_matrix(z1)
end
end

@testitem "Verification V4: every compiler-emitted gadget lifts" tags=[:verification, :v4] begin
using PBCCompiler
using PBCCompiler: CircuitOp, preprocess_circuit, parse_input
using PBCCompiler.Verification: lift_gadgets, lower_lifted, random_clifford_t_circuit
using PBCCompiler.GadgetShapes: find_gadgets, count_gadgets
using QuantumClifford: @P_str
using Random

function check_all_lift(input)
    n = PBCCompiler.Verification._width(input)
    out = preprocess_circuit(copy(input))
    lr = lift_gadgets(out, n)
    nmagic = count(op -> PBCCompiler.Verification._opkind(op) == PBCCompiler.Verification.OP_MEASUREMENT &&
                         all(>(n), PBCCompiler.Verification._support(op)), out)
    ok = isempty(lr.failed) && length(lr.lifted) == nmagic &&
         all(g -> g.angle == 1 // 8 && iszero(g.residual), lr.lifted) &&
         isempty(lower_lifted(out, lr).unsupported)
    ok, length(lr.lifted)
end

ok, k = check_all_lift([CircuitOp.ExpEighPiPauli(P"Z", [1]), CircuitOp.Measurement(P"Z", 1, [1])])
@test ok && k == 1
ok, k = check_all_lift([CircuitOp.ExpEighPiPauli(P"-Z", [1]), CircuitOp.ExpQuatPiPauli(P"X", [1]),
                        CircuitOp.ExpEighPiPauli(P"YX", [1, 2]), CircuitOp.Measurement(P"Z", 1, [1]),
                        CircuitOp.Measurement(P"Z", 2, [2])])
@test ok && k == 2
fx = parse_input(joinpath(pkgdir(PBCCompiler), "test", "fixtures", "sample_input.qasm"))
ok, k = check_all_lift(fx)
@test ok && k == 111
rng = MersenneTwister(404)
total = Ref(0)
for _ in 1:40
    n = rand(rng, 1:5)
    c = random_clifford_t_circuit(rng, n; depth=rand(rng, 1:5n), t_density=0.4, multiqubit_density=0.25,
                                  sorted_qubits=true)
    ok, k = check_all_lift(c)
    @test ok
    total[] += k
end
@test total[] > 40
# the shared recognizer counts magic states by template
out = preprocess_circuit(copy(fx))
@test count_gadgets(out, 2) == Dict(:pi8_autocorrected => 111)
end

@testitem "Verification V4: corrupted gadgets are rejected with a reason" tags=[:verification, :v4] begin
using PBCCompiler: CircuitOp, preprocess_circuit
using PBCCompiler.Verification: lift_gadgets
using QuantumClifford: @P_str
const C = CircuitOp

good = preprocess_circuit([C.ExpEighPiPauli(P"Z", [1]), C.Measurement(P"Z", 1, [1])])
# good = [M(ZZ,2,[1,2]), M(X,3,[2]), BC(Quat Z,2), BC(Half Z,3), M(Z,1,[1])]
@test length(lift_gadgets(good, 1).lifted) == 1

function lift_one(c; magic_k=1)
    lr = lift_gadgets(c, 1; magic_k)
    @test length(lr.lifted) + length(lr.failed) <= 1
    lr
end
reason(lr) = isempty(lr.failed) ? "" : lr.failed[1].reason

# wrong correction: π/4 replaced by π/2
c = copy(good); c[3] = C.BitConditional(C.ExpHalfPiPauli(P"Z", [1]), 2)
lr = lift_one(c); @test isempty(lr.lifted) && occursin("different rotations", reason(lr))
# wrong conditioning bits (swapped)
c = copy(good); c[3] = C.BitConditional(C.ExpQuatPiPauli(P"Z", [1]), 3); c[4] = C.BitConditional(C.ExpHalfPiPauli(P"Z", [1]), 2)
lr = lift_one(c); @test isempty(lr.lifted) && occursin("different rotations", reason(lr))
# missing π/4 correction
c = deleteat!(copy(good), 3)
lr = lift_one(c); @test isempty(lr.lifted) && occursin("different rotations", reason(lr))
# correction about a different Pauli
c = copy(good); c[3] = C.BitConditional(C.ExpQuatPiPauli(P"X", [1]), 2)
lr = lift_one(c); @test isempty(lr.lifted) && occursin("not ±", reason(lr))
# wrong magic state |+⟩ (k = 0): the bare gadget only applies Z^{b2}, so after the
# π/4 correction the branches implement different rotations
lr = lift_one(good; magic_k=0); @test isempty(lr.lifted) && occursin("different rotations", reason(lr))
# wrong magic state S|+⟩ (k = 2): a valid gadget, but for a different (Clifford) angle
lr = lift_one(good; magic_k=2); @test isempty(lr.lifted) || lr.lifted[1].angle != 1 // 8
# magic qubit measured in Z instead of X: the shape is not recognized at all
c = copy(good); c[2] = C.Measurement(P"Z", 3, [2])
lr = lift_one(c); @test isempty(lr.lifted) && isempty(lr.failed)
# missing π/2 correction: lifts, but leaves the residual Pauli frame P^{b2}
c = deleteat!(copy(good), 4)
lr = lift_one(c); @test length(lr.lifted) == 1 && lr.lifted[1].residual.vars == [2]
# π/4 correction with the opposite sign: equivalent up to the Pauli frame P^{b1}
c = copy(good); c[3] = C.BitConditional(C.ExpQuatPiPauli(P"-Z", [1]), 2)
lr = lift_one(c); @test length(lr.lifted) == 1 && lr.lifted[1].angle == 1 // 8 && lr.lifted[1].residual.vars == [1]
end

@testitem "Verification V4: exact gadget prover vs independent floating-point check" tags=[:verification, :v4] begin
using PBCCompiler.Verification: prove_reduced_gadget
using LinearAlgebra, Random

L = Dict('X' => ComplexF64[0 1; 1 0], 'Y' => ComplexF64[0 -im; im 0], 'Z' => ComplexF64[1 0; 0 -1])
I2 = Matrix{ComplexF64}(I, 2, 2)
eig(l, b) = (v = eigen(Hermitian(L[l])); v.vectors[:, findfirst(x -> isapprox(x, b ? -1 : 1), v.values)])
rot(θ) = cispi(-θ) * [1 0; 0 0] + cispi(θ) * [0 0; 0 1]        # exp(-iπθZ)
# float data Kraus operator of branch (b1, b2); index = d + 2m, i.e. kron(magic, data)
function kraus(jl, ml, ms, corr, k, b1, b2)
    J = kron(L[jl], L['Z']); T = kron(ms ? -L[ml] : L[ml], I2)
    m = ComplexF64[1, cispi(k / 4)] / sqrt(2)
    M = (I + (b2 ? -1 : 1) * T) / 2 * (I + (b1 ? -1 : 1) * J) / 2 * kron(m, I2)
    for (θ, role) in corr
        (role == 1 ? b1 : b2) && (M = kron(I2, rot(θ)) * M)
    end
    kron(eig(ml, b2 ⊻ ms)', I2) * M
end
proportional(A, B) = (c = tr(B' * A) / tr(B' * B); norm(A - c * B) < 1e-10)

rng = MersenneTwister(12)
nok = Ref(0); nfail = Ref(0)
for jl in "XYZ", ml in "XYZ", ms in (false, true), k in 0:7
    corr = Tuple{Rational{Int},Int}[(rand(rng, -3:4) // 8, rand(rng, 1:2)) for _ in 1:rand(rng, 0:2)]
    pr = prove_reduced_gadget(jl, ml, ms, corr, k)
    Ks = [kraus(jl, ml, ms, corr, k, b1, b2) for b1 in (false, true), b2 in (false, true)]
    if pr.ok
        nok[] += 1
        c0, c1, c2 = pr.residual
        good = true
        for b1 in (false, true), b2 in (false, true)
            K = Ks[b1 + 1, b2 + 1]
            norm(K) < 1e-12 && continue
            r = c0 ⊻ (c1 & b1) ⊻ (c2 & b2)
            good &= proportional(K, (r ? L['Z'] : I2) * rot(pr.angle))
        end
        @test good
        @test sum(K' * K for K in Ks) ≈ I                   # trace preserving
    else
        nfail[] += 1
        # some branch is not a common rotation up to Z: there is no θ ∈ π/8·Z that works
        works = any(0:15) do j
            θ = j // 16
            all(norm(K) < 1e-12 || proportional(K, rot(θ)) || proportional(K, L['Z'] * rot(θ)) for K in Ks)
        end
        @test !works || !(sum(K' * K for K in Ks) ≈ I)
    end
end
@test nok[] > 10 && nfail[] > 10
end
