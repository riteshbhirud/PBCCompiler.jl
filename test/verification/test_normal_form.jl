@testmodule NFTestUtils begin
    using PBCCompiler.Verification
    using PBCCompiler.Verification: DensePauli, apply_rotation!, apply_projector!, dense_matrix, dense_unitary,
        _dense_dfs, ITEM_MEASURE, ITEM_ROTATION
    using LinearAlgebra

    # Global-phase-insensitive equality of two unitaries.
    phase_equal(A, B; atol=1e-9) = abs(abs(tr(A' * B)) - size(A, 1)) < atol * size(A, 1)

    # POVM of a circuit on n data qubits, keyed by the final values of all written bits.
    function circuit_povm(c, n)
        E = Dict{Vector{Bool},Matrix{ComplexF64}}()
        K0 = Matrix{ComplexF64}(I, 1 << n, 1 << n)
        _dense_dfs(c, n, K0) do K, reg, outs
            key = [reg[b] for b in sort!(collect(keys(reg)))]
            E[key] = get(E, key, zeros(ComplexF64, 1 << n, 1 << n)) + K' * K
        end
        E
    end

    # POVM of the items of a normal form (frames are irrelevant for outcomes).
    function items_povm(nf::NormalForm)
        n = nf.n
        E = Dict{Vector{Bool},Matrix{ComplexF64}}()
        bits = sort!(collect(keys(nf.final_var)))
        function rec(k, K, assign)
            if k > length(nf.items)
                key = [assign[nf.final_var[b]] for b in bits]
                E[key] = get(E, key, zeros(ComplexF64, 1 << n, 1 << n)) + K' * K
                return
            end
            it = nf.items[k]
            if it.kind == ITEM_MEASURE
                s = evaluate(it.sign, assign)
                for v in (false, true)
                    K2 = copy(K)
                    apply_projector!(K2, DensePauli(it.pauli), v ⊻ s)
                    sum(abs2, K2) < 1e-28 && continue
                    a2 = copy(assign); a2[it.var] = v
                    rec(k + 1, K2, a2)
                end
            else
                K2 = copy(K)
                if evaluate(it.cond, assign)
                    θ = evaluate(it.sign, assign) ? -it.angle : it.angle
                    apply_rotation!(K2, DensePauli(it.pauli), θ)
                end
                rec(k + 1, K2, assign)
            end
        end
        rec(1, Matrix{ComplexF64}(I, 1 << n, 1 << n), Dict{Int,Bool}())
        E
    end

    function povm_equal(A, B; atol=1e-9)
        for k in union(keys(A), keys(B))
            a = get(A, k, nothing); b = get(B, k, nothing)
            a === nothing && (a = zero(b)); b === nothing && (b = zero(a))
            isapprox(a, b; atol=atol) || return false
        end
        true
    end
end

@testitem "Verification V3: PauliMap vs dense conjugation" tags=[:verification, :v3] begin
using PBCCompiler.Verification: identity_map, rotate!, apply_map, compose, to_clifford, DensePauli,
    apply_rotation!, dense_matrix
using QuantumClifford: random_pauli, Stabilizer, apply!, PauliOperator, nqubits
using LinearAlgebra, Random

rng = MersenneTwister(31)
gens(n) = vcat([(p = zero(PauliOperator, n); p[j] = (true, false); p) for j in 1:n],
               [(p = zero(PauliOperator, n); p[j] = (false, true); p) for j in 1:n])
for trial in 1:40
    n = rand(rng, 1:4)
    fwd = identity_map(n); inv = identity_map(n)
    U = Matrix{ComplexF64}(I, 1 << n, 1 << n)
    for _ in 1:8
        q = random_pauli(rng, n; nophase=false, realphase=true)
        k = rand(rng, 1:3)
        rotate!(fwd, inv, q, k)
        apply_rotation!(U, DensePauli(q), k // 4)
    end
    probes = vcat(gens(n), [random_pauli(rng, n; nophase=false) for _ in 1:4])
    @test all(dense_matrix(apply_map(fwd, P)) ≈ U * dense_matrix(P) * U' for P in probes)
    @test all(dense_matrix(apply_map(inv, P)) ≈ U' * dense_matrix(P) * U for P in probes)
    # cross-check against QuantumClifford's own Clifford action (phases included)
    C = to_clifford(fwd)
    @test all(apply!(Stabilizer([P]), C)[1] == apply_map(fwd, P) for P in probes if iseven(P.phase[]))
    # composition
    fwd2 = identity_map(n); inv2 = identity_map(n); U2 = Matrix{ComplexF64}(I, 1 << n, 1 << n)
    for _ in 1:3
        q = random_pauli(rng, n; nophase=true); k = rand(rng, 1:3)
        rotate!(fwd2, inv2, q, k); apply_rotation!(U2, DensePauli(q), k // 4)
    end
    W = U * U2
    @test all(dense_matrix(apply_map(compose(fwd, fwd2), P)) ≈ W * dense_matrix(P) * W' for P in probes)
end
end

@testitem "Verification V3: PauliConditional lowering vs projector definition" tags=[:verification, :v3] setup=[NFTestUtils] begin
using PBCCompiler: CircuitOp
using PBCCompiler.Verification: lower, DensePauli, apply_rotation!, dense_unitary
using QuantumClifford: random_pauli
using LinearAlgebra, Random
using .NFTestUtils: phase_equal

rng = MersenneTwister(5)
for trial in 1:25
    n = rand(rng, 2:4)
    perm = randperm(rng, n)
    k = rand(rng, 1:n-1)
    cq = perm[1:k]; tq = perm[k+1:end]
    cp = random_pauli(rng, length(cq); nophase=false, realphase=true)
    tp = random_pauli(rng, length(tq); nophase=false, realphase=true)
    c = [CircuitOp.PauliConditional(cp, cq, tp, tq)]
    Ud = dense_unitary(c, n)
    prog = lower(c, n)
    @test isempty(prog.unsupported)
    V = Matrix{ComplexF64}(I, 1 << n, 1 << n)
    for v in prog.ops
        apply_rotation!(V, DensePauli(v.pauli), v.angle)
    end
    @test phase_equal(Ud, V)
end
end

@testitem "Verification V3: normal form of unitary Clifford+T circuits (dense)" tags=[:verification, :v3] setup=[NFTestUtils] begin
using PBCCompiler.Verification: normal_form, random_clifford_t_circuit, apply_map, dense_unitary, dense_matrix,
    DensePauli, apply_rotation!, evaluate, isconstant
using QuantumClifford: PauliOperator, random_pauli
using LinearAlgebra, Random

rng = MersenneTwister(2026)
gens(n) = vcat([(p = zero(PauliOperator, n); p[j] = (true, false); p) for j in 1:n],
               [(p = zero(PauliOperator, n); p[j] = (false, true); p) for j in 1:n])
for trial in 1:60
    n = rand(rng, 1:6)
    c = random_clifford_t_circuit(rng, n; depth=rand(rng, 1:4n+4), final_measure=false,
                                  t_density=0.3, multiqubit_density=0.3)
    nf = normal_form(c; n=n)
    U = dense_unitary(c, n)
    V = Matrix{ComplexF64}(I, 1 << n, 1 << n)
    for it in nf.items
        @test isconstant(it.sign) && isone(it.cond)
        θ = evaluate(it.sign, Bool[]) ? -it.angle : it.angle
        apply_rotation!(V, DensePauli(it.pauli), θ)
    end
    @test isone(nf.pauli_frame)
    W = U * V'                         # must be the Clifford frame F (up to phase)
    @test all(W * dense_matrix(e) * W' ≈ dense_matrix(apply_map(nf.fwd, e)) for e in gens(n))
    @test all(W' * dense_matrix(e) * W ≈ dense_matrix(apply_map(nf.inv, e)) for e in gens(n))
    # every item is non-Clifford (Clifford rotations are always absorbed)
    @test all(denominator(it.angle) == 8 for it in nf.items)
end
end

@testitem "Verification V3: normal-form items reproduce the outcome POVM" tags=[:verification, :v3] setup=[NFTestUtils] begin
using PBCCompiler.Verification: normal_form, random_clifford_t_circuit
using Random
using .NFTestUtils: circuit_povm, items_povm, povm_equal

rng = MersenneTwister(77)
nconditional = Ref(0)
for trial in 1:60
    n = rand(rng, 1:4)
    c = random_clifford_t_circuit(rng, n; depth=rand(rng, 2:3n+3), meas_density=0.2, cc_density=0.25,
                                  t_density=0.3, multiqubit_density=0.25)
    nf = normal_form(c; n=n)
    nconditional[] += count(it -> !isone(it.cond) || !isempty(it.sign.vars), nf.items)
    @test povm_equal(circuit_povm(c, n), items_povm(nf))
end
@test nconditional[] > 0      # the suite really exercises outcome-dependent items
end

@testitem "Verification V3: unsupported constructs" tags=[:verification, :v3] begin
using PBCCompiler: CircuitOp
using PBCCompiler.Verification: normal_form, lower
using QuantumClifford: @P_str
@test_throws ArgumentError normal_form([CircuitOp.PrepMagic(1, [2]), CircuitOp.Measurement(P"Z", 1, [1])]; n=2)
prog = lower([CircuitOp.BitConditional(CircuitOp.BitConditional(CircuitOp.ExpHalfPiPauli(P"X", [1]), 1), 1)], 1)
@test length(prog.unsupported) == 1
end
