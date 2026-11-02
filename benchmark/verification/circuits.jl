# Benchmark circuit families in PBCCompiler's IR (shared by the verification and
# estimation benchmark scripts).
#
# CNOTs whose control index exceeds the target index are emitted as
# (H⊗H)·CNOT(t→c)·(H⊗H), an exact identity, so that compiler bug B1 (see
# docs/NOTES_FOR_COMPILER.md) does not turn every benchmark into a refutation.

using PBCCompiler: CircuitOp
using QuantumClifford: PauliOperator

const _Zp = PauliOperator(0x00, [false], [true])
const _mZp = PauliOperator(0x02, [false], [true])
const _Xp = PauliOperator(0x00, [true], [false])

bm_h(q) = [CircuitOp.ExpQuatPiPauli(_Zp, [q]), CircuitOp.ExpQuatPiPauli(_Xp, [q]), CircuitOp.ExpQuatPiPauli(_Zp, [q])]
bm_t(q) = [CircuitOp.ExpEighPiPauli(_Zp, [q])]
bm_tdg(q) = [CircuitOp.ExpEighPiPauli(_mZp, [q])]
function bm_cnot(c, t)
    c < t && return [CircuitOp.PauliConditional(_Zp, [c], _Xp, [t])]
    vcat(bm_h(c), bm_h(t), [CircuitOp.PauliConditional(_Zp, [t], _Xp, [c])], bm_h(c), bm_h(t))
end
bm_measure_all(n; bit0 = 0) = [CircuitOp.Measurement(_Zp, bit0 + q, [q]) for q in 1:n]

"""Controlled-S between qubits a and b: T_a T_b CNOT(a,b) T†_b CNOT(a,b)."""
bm_cs(a, b) = vcat(bm_t(a), bm_t(b), bm_cnot(min(a, b), max(a, b)), bm_tdg(max(a, b)), bm_cnot(min(a, b), max(a, b)))

"""Standard 7-T Toffoli decomposition (controls a, b; target c)."""
function bm_toffoli(a, b, c)
    vcat(bm_h(c), bm_cnot(b, c), bm_tdg(c), bm_cnot(a, c), bm_t(c), bm_cnot(b, c), bm_tdg(c), bm_cnot(a, c),
         bm_t(b), bm_t(c), bm_h(c), bm_cnot(a, b), bm_t(a), bm_tdg(b), bm_cnot(a, b))
end

"""Approximate QFT on n qubits keeping only the controlled-S (R₂) rotations (exact in Clifford+T)."""
function approx_qft(n)
    c = CircuitOp.Type[]
    for j in 1:n
        append!(c, bm_h(j))
        j < n && append!(c, bm_cs(j + 1, j))
    end
    append!(c, bm_measure_all(n))
    c
end

"""
Cuccaro ripple-carry adder on two k-bit registers plus a carry-in and a carry-out
(n = 2k + 2 qubits), with Toffolis decomposed into Clifford+T. All qubits are
measured at the end.
"""
function cuccaro_adder(k)
    # layout: cin = 1, a_i = 2i, b_i = 2i + 1 (i = 1..k), cout = 2k + 2
    a(i) = 2i; b(i) = 2i + 1; cin = 1; cout = 2k + 2
    c = CircuitOp.Type[]
    maj(x, y, z) = vcat(bm_cnot(z, y), bm_cnot(z, x), bm_toffoli(x, y, z))
    uma(x, y, z) = vcat(bm_toffoli(x, y, z), bm_cnot(z, x), bm_cnot(x, y))
    prev = cin
    for i in 1:k
        append!(c, maj(prev, b(i), a(i))); prev = a(i)
    end
    append!(c, bm_cnot(a(k), cout))
    for i in k:-1:1
        append!(c, uma(i == 1 ? cin : a(i - 1), b(i), a(i)))
    end
    append!(c, bm_measure_all(2k + 2))
    c
end
