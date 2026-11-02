@testitem "Verification V2: AffineBit" tags=[:verification, :v2] begin

using PBCCompiler.Verification: AffineBit, affine_var, isconstant, evaluate, substitute

@test AffineBit(true, [3, 1, 3, 2]) == AffineBit(true, [1, 2])          # repeated var cancels
@test affine_var(1) ⊻ affine_var(1) == AffineBit(false)
@test (affine_var(2) ⊻ true) ⊻ affine_var(5) == AffineBit(true, [2, 5])
@test isconstant(AffineBit(true)) && !isconstant(affine_var(1))
@test evaluate(AffineBit(true, [1, 3]), [true, false, true]) == true
@test evaluate(AffineBit(false, [1, 3]), Dict(1 => true, 3 => false)) == true
# substitution: m1 -> 1 ⊕ m4, m3 -> m4
s = substitute(AffineBit(false, [1, 3]), v -> v == 1 ? AffineBit(true, [4]) : v == 3 ? affine_var(4) : affine_var(v))
@test s == AffineBit(true)
@test hash(AffineBit(true, [2, 1])) == hash(AffineBit(true, [1, 2]))
@test sprint(show, AffineBit(true, [1, 3])) == "1 ⊕ m1 ⊕ m3"

end

@testitem "Verification V2: symbolic Pauli frame vs explicit branches" tags=[:verification, :v2] begin

using PBCCompiler.Verification: AffineBit, SymbolicPauli, mul_conditional!, mul_frame!, conjugate!,
    commutation_parity, evaluate, frame_variables
using QuantumClifford: PauliOperator, Stabilizer, apply!, comm, random_clifford, random_pauli, nqubits
using Random

const M = 12              # outcome variables
const NB = 1 << M         # all assignments

assignment(k) = BitVector(((k >> (i - 1)) & 1) == 1 for i in 1:M)
rand_affine(rng; p=0.3) = AffineBit(rand(rng, Bool), [v for v in 1:M if rand(rng) < p])
samebits(a::PauliOperator, b::PauliOperator) = a.xz == b.xz          # equal up to phase

rng = MersenneTwister(20261102)
for trial in 1:8
    n = rand(rng, 1:5)
    f = SymbolicPauli(n)
    # random symbolic start
    for j in 1:n
        f.x[j] = rand_affine(rng); f.z[j] = rand_affine(rng)
    end
    start = copy(f)
    # random operation sequence, recorded for explicit replay
    ops = Any[]
    for _ in 1:6
        if rand(rng, Bool)
            q = random_pauli(rng, n; nophase=true)
            c = rand_affine(rng)
            mul_conditional!(f, q, c)
            push!(ops, (:mul, q, c))
        else
            C = random_clifford(rng, n)
            conjugate!(f, C)
            push!(ops, (:conj, C))
        end
    end
    probe = random_pauli(rng, n; nophase=true)
    par = commutation_parity(f, probe)
    ok_frame = true; ok_parity = true
    for k in 0:NB-1
        a = assignment(k)
        p = evaluate(start, a)
        for op in ops
            if op[1] == :mul
                evaluate(op[3], a) && (p = op[2] * p)
            else
                p = apply!(Stabilizer([p]), op[2])[1]
            end
        end
        ok_frame &= samebits(p, evaluate(f, a))
        ok_parity &= (comm(p, probe) == 0x01) == evaluate(par, a)
    end
    @test ok_frame
    @test ok_parity
end

@testset "equality as functions" begin
    rng = MersenneTwister(7)
    for _ in 1:20
        n = rand(rng, 1:3)
        f = SymbolicPauli(n); g = SymbolicPauli(n)
        for j in 1:n
            f.x[j] = rand_affine(rng; p=0.15); f.z[j] = rand_affine(rng; p=0.15)
        end
        copyto!(g.x, f.x); copyto!(g.z, f.z)
        rand(rng, Bool) && (g.z[rand(rng, 1:n)] ⊻= rand_affine(rng; p=0.15))
        exhaustive = all(evaluate(f, assignment(k)).xz == evaluate(g, assignment(k)).xz for k in 0:NB-1)
        @test (f == g) == exhaustive
    end
    # building the same function two ways
    f = SymbolicPauli(2); g = SymbolicPauli(2)
    c1 = AffineBit(false, [1]); c2 = AffineBit(true, [2])
    mul_conditional!(f, PauliOperator(0x00, [true, false], [false, true]), c1)
    mul_conditional!(f, PauliOperator(0x00, [true, false], [false, true]), c2)
    mul_conditional!(g, PauliOperator(0x00, [true, false], [false, true]), c1 ⊻ c2)
    @test f == g
    h = SymbolicPauli(2); mul_frame!(h, f)
    @test h == f
    @test frame_variables(f) == [1, 2]
end

end
