@testitem "V9: benchmark circuit generators are correct (dense)" tags=[:verification, :v9] begin
using PBCCompiler
using PBCCompiler.Verification: dense_unitary
using LinearAlgebra
include(joinpath(pkgdir(PBCCompiler), "benchmark", "verification", "circuits.jl"))

phase_equal(A, B) = abs(abs(tr(A' * B)) - size(A, 1)) < 1e-9 * size(A, 1)
# little-endian basis: qubit j is bit j-1
function perm_unitary(n, f)
    U = zeros(ComplexF64, 1 << n, 1 << n)
    for x in 0:(1 << n)-1
        U[f(x) + 1, x + 1] = 1
    end
    U
end
bit(x, j) = (x >> (j - 1)) & 1

# CNOT in both orientations
for (c, t) in ((1, 2), (2, 1), (3, 1))
    U = perm_unitary(3, x -> bit(x, c) == 1 ? x ⊻ (1 << (t - 1)) : x)
    @test phase_equal(dense_unitary(bm_cnot(c, t), 3), U)
end
# controlled-S: phase i on |11⟩
CS = Diagonal(ComplexF64[1, 1, 1, im])
@test phase_equal(dense_unitary(bm_cs(1, 2), 2), Matrix(CS))
@test phase_equal(dense_unitary(bm_cs(2, 1), 2), Matrix(CS))
# Toffoli
for (a, b, c) in ((1, 2, 3), (3, 1, 2))
    U = perm_unitary(3, x -> (bit(x, a) & bit(x, b)) == 1 ? x ⊻ (1 << (c - 1)) : x)
    @test phase_equal(dense_unitary(bm_toffoli(a, b, c), 3), U)
end
# Cuccaro adder: b = a + b + cin, carry in cout
for k in 1:2
    c = filter(op -> PBCCompiler.Verification._opkind(op) != PBCCompiler.Verification.OP_MEASUREMENT, cuccaro_adder(k))
    n = 2k + 2
    U = dense_unitary(c, n)
    ok = true
    for A in 0:(1 << k)-1, B in 0:(1 << k)-1, cin in 0:1
        x = cin
        for i in 1:k
            x |= ((A >> (i - 1)) & 1) << (2i - 1)       # a_i on qubit 2i
            x |= ((B >> (i - 1)) & 1) << (2i)           # b_i on qubit 2i+1
        end
        y = argmax(abs.(U[:, x + 1])) - 1
        S = A + B + cin
        got = sum(((y >> (2i)) & 1) << (i - 1) for i in 1:k) + (((y >> (2k + 1)) & 1) << k)
        gotA = sum(((y >> (2i - 1)) & 1) << (i - 1) for i in 1:k)
        ok &= abs(abs(U[y + 1, x + 1]) - 1) < 1e-9 && got == S && gotA == A && (y & 1) == cin
    end
    @test ok
end
# approximate QFT-2 is H, CS, H (exact QFT-2 without the final swap)
H = [1 1; 1 -1] / sqrt(2)
qft2 = kron(H, I(2)) * Matrix(CS) * kron(I(2), H)       # qubit 2 is the most significant factor
c = filter(op -> PBCCompiler.Verification._opkind(op) != PBCCompiler.Verification.OP_MEASUREMENT, approx_qft(2))
@test phase_equal(dense_unitary(c, 2), qft2)
end
