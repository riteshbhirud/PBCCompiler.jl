# Dense simulation primitives (ComplexF64), shared by the dense oracle and the tests.
#
# Qubit j corresponds to bit (j-1) of the row index (little-endian). Operators act
# on the rows of a matrix K whose columns are states, so simulating the operator
# itself is K = I, and simulating a set of input states is K = [ψ₁ ψ₂ …].
#
# This simulator gives every CircuitOp its semantics *directly* (for example
# PauliConditional via its projector definition), independently of the lowering
# used by the normal form. It is therefore a genuine cross-check.

const DENSE_MAX_QUBITS = 30

"""A Pauli operator as bit masks: P|r⟩ = coeff·(−1)^{popcount(r & zmask)}·|r ⊻ xmask⟩."""
struct DensePauli
    xmask::Int
    zmask::Int
    coeff::ComplexF64
end

function DensePauli(p::PauliOperator)
    n = nqubits(p)
    n <= DENSE_MAX_QUBITS || throw(ArgumentError("dense simulation supports at most $DENSE_MAX_QUBITS qubits"))
    xm = 0; zm = 0; ny = 0
    for j in 1:n
        x, z = p[j]
        x && (xm |= 1 << (j - 1))
        z && (zm |= 1 << (j - 1))
        (x && z) && (ny += 1)
    end
    DensePauli(xm, zm, ComplexF64(im)^mod(Int(p.phase[]) + ny, 4))
end

"""out ← P·K (out and K must not alias)."""
function pauli_mul!(out::AbstractMatrix{ComplexF64}, P::DensePauli, K::AbstractMatrix{ComplexF64})
    nr = size(K, 1)
    @inbounds for col in axes(K, 2), r in 0:nr-1
        c = isodd(count_ones(r & P.zmask)) ? -P.coeff : P.coeff
        out[(r ⊻ P.xmask) + 1, col] = c * K[r + 1, col]
    end
    out
end

"""K ← exp(-iπθ·P)·K."""
function apply_rotation!(K::AbstractMatrix{ComplexF64}, P::DensePauli, θ::Rational{Int})
    tmp = similar(K)
    pauli_mul!(tmp, P, K)
    c = cospi(Float64(θ)); s = sinpi(Float64(θ))
    @. K = c * K - im * s * tmp
    K
end

"""K ← ½(I + (−1)^outcome P)·K, the projector onto the outcome-`outcome` eigenspace."""
function apply_projector!(K::AbstractMatrix{ComplexF64}, P::DensePauli, outcome::Bool)
    tmp = similar(K)
    pauli_mul!(tmp, P, K)
    s = outcome ? -0.5 : 0.5
    @. K = 0.5 * K + s * tmp
    K
end

"""K ← C·K with C = (I + P₁ + P₂ − P₁P₂)/2 (apply P₂ when P₁ = −1)."""
function apply_pauliconditional!(K::AbstractMatrix{ComplexF64}, P1::DensePauli, P2::DensePauli)
    t1 = similar(K); t2 = similar(K); t12 = similar(K)
    pauli_mul!(t1, P1, K)
    pauli_mul!(t2, P2, K)
    pauli_mul!(t12, P1, t2)
    @. K = 0.5 * (K + t1 + t2 - t12)
    K
end

"""Dense matrix of a Pauli operator (for tests)."""
function dense_matrix(p::PauliOperator)
    d = 1 << nqubits(p)
    K = Matrix{ComplexF64}(LinearAlgebra.I, d, d)
    out = similar(K)
    pauli_mul!(out, DensePauli(p), K)
end

"""The magic state (|0⟩ + ζ₈^k |1⟩)/√2 as a vector (k = 1 is T|+⟩)."""
_magic_vector(k::Int) = ComplexF64[1, cispi(k / 4)] / sqrt(2)

"""
Initial operator: columns are the 2^n_data computational basis inputs, tensored with
the magic state on each qubit n_data+1 … N.
"""
function _initial_operator(n_data::Int, N::Int, magic_k::Int)
    N <= DENSE_MAX_QUBITS || throw(ArgumentError("too many qubits for dense simulation"))
    D = 1 << n_data
    K = zeros(ComplexF64, 1 << N, D)
    mv = _magic_vector(magic_k)
    nm = N - n_data
    for m in 0:(1 << nm)-1
        amp = one(ComplexF64)
        for t in 0:nm-1
            amp *= mv[((m >> t) & 1) + 1]
        end
        for j in 0:D-1
            K[(j | (m << n_data)) + 1, j + 1] = amp
        end
    end
    K
end

_dense_embed(p::PauliT, q::Vector{Int}, N::Int) = DensePauli(embed(N, q, p))

# Apply a non-measurement op given the current register; errors on unsupported ops.
function _dense_apply_op!(K, op, N::Int, reg::Dict{Int,Bool})
    k = _opkind(op)
    if _isrotation(k)
        p, q = _pauli_and_qubits(op)
        apply_rotation!(K, _dense_embed(p, q, N), _rotation_angle(k))
    elseif k == OP_PAULICOND
        cp, cq, tp, tq = _paulicond_parts(op)
        apply_pauliconditional!(K, _dense_embed(cp, cq, N), _dense_embed(tp, tq, N))
    elseif k == OP_BITCOND
        inner, b = _bitcond_parts(op)
        haskey(reg, b) || throw(ArgumentError("dense simulation: bit $b read before write"))
        _opkind(inner) == OP_MEASUREMENT && throw(ArgumentError("dense simulation: conditional measurement unsupported"))
        reg[b] && _dense_apply_op!(K, inner, N, reg)
    elseif k == OP_MEASUREMENT
        throw(ArgumentError("internal: measurement passed to _dense_apply_op!"))
    else
        throw(ArgumentError("dense simulation: PrepMagic has no defined semantics"))
    end
    K
end

"""
    dense_unitary(c, N::Int) -> Matrix{ComplexF64}

The 2^N × 2^N unitary of a measurement-free circuit (for tests).
"""
function dense_unitary(c::AbstractVector{CircuitOp.Type}, N::Int)
    K = Matrix{ComplexF64}(LinearAlgebra.I, 1 << N, 1 << N)
    reg = Dict{Int,Bool}()
    for op in c
        _opkind(op) == OP_MEASUREMENT && throw(ArgumentError("dense_unitary: circuit has a measurement"))
        _dense_apply_op!(K, op, N, reg)
    end
    K
end

"""
    _dense_dfs(visit, c, N, K0; allow, max_leaves, prune_tol)

Depth-first enumeration of all measurement-outcome branches of `c`, starting
from the operator `K0` (rows: 2^N). At each measurement of bit `b` with candidate
outcome `v`, the branch is explored only if `allow(op_index, b, v, reg)` is
true. Branches whose operator has squared Frobenius norm `≤ prune_tol` are
dropped (they contribute nothing). At each leaf, `visit(K, reg, outcomes)` is
called, where `outcomes` lists the `(op index, bit, value)` triples of that
branch. Returns `(leaves_visited, complete)`; `complete` is false if
`max_leaves` was hit.
"""
function _dense_dfs(visit, c::AbstractVector{CircuitOp.Type}, N::Int, K0::Matrix{ComplexF64};
                    allow = (i, b, v, reg) -> true, max_leaves::Int = typemax(Int), prune_tol::Float64 = 1e-28)
    counter = Ref(0)
    complete = Ref(true)
    _dense_dfs_rec(visit, c, N, 1, copy(K0), Dict{Int,Bool}(), Tuple{Int,Int,Bool}[], allow, max_leaves, prune_tol, counter, complete)
    (counter[], complete[])
end

function _dense_dfs_rec(visit, c, N, i, K, reg, outs, allow, max_leaves, prune_tol, counter, complete)
    complete[] || return
    while i <= length(c)
        op = c[i]
        if _opkind(op) == OP_MEASUREMENT
            p, q = _pauli_and_qubits(op)
            b = _measurement_bit(op)
            P = _dense_embed(p, q, N)
            for v in (false, true)
                allow(i, b, v, reg) || continue
                K2 = copy(K)
                apply_projector!(K2, P, v)
                sum(abs2, K2) <= prune_tol && continue
                reg2 = copy(reg); reg2[b] = v
                outs2 = copy(outs); push!(outs2, (i, b, v))
                _dense_dfs_rec(visit, c, N, i + 1, K2, reg2, outs2, allow, max_leaves, prune_tol, counter, complete)
                complete[] || return
            end
            return
        end
        _dense_apply_op!(K, op, N, reg)
        i += 1
    end
    if counter[] >= max_leaves
        complete[] = false
        return
    end
    counter[] += 1
    visit(K, reg, outs)
    nothing
end
