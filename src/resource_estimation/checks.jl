const EPauli = PauliOperator{Array{UInt8,0},Vector{UInt64}}

_measurement_parts(op) = @match op begin
    CircuitOp.Measurement(p, b, q) => (p, b, q)
    _ => nothing
end

_bitcond_parts_est(op) = @match op begin
    CircuitOp.BitConditional(inner, b) => (inner, b)
    _ => nothing
end

_rotation_parts(op) = @match op begin
    CircuitOp.Pauli(p, q) => (p, q, 1 // 2)
    CircuitOp.ExpHalfPiPauli(p, q) => (p, q, 1 // 2)
    CircuitOp.ExpQuatPiPauli(p, q) => (p, q, 1 // 4)
    CircuitOp.ExpEighPiPauli(p, q) => (p, q, 1 // 8)
    _ => nothing
end

_paulicond_parts_est(op) = @match op begin
    CircuitOp.PauliConditional(cp, cq, tp, tq) => (cp, cq, tp, tq)
    _ => nothing
end

_pauli_support(p::PauliOperator, qubits::Vector{Int}) = sort!([qubits[j] for j in eachindex(qubits) if p[j][1] || p[j][2]])

_est_width(c) = maximum((isempty(affectedqubits(op)) ? 0 : maximum(affectedqubits(op)) for op in c); init = 0)

"""
    parity_check_matrices(m::CodeModel) -> (Hx::BitMatrix, Hz::BitMatrix)

Explicit X- and Z-check matrices of the model's code.

- Surface: the rotated surface code of distance `d`, with qubit `(i, j)` at index
  `i*d + j + 1` (row-major); X boundaries on the top and bottom, Z on the left and right.
- BB: `Hx = [A | B]`, `Hz = [B' | A']` (Bravyi et al. 2024, Eq. 1).
- Generic: taken from the wrapped QuantumClifford code.
"""
function parity_check_matrices(m::SurfaceCodeModel)
    d = m.d
    q(i, j) = i * d + j + 1
    X = Vector{Vector{Int}}(); Z = Vector{Vector{Int}}()
    for i in 0:d-2, j in 0:d-2
        s = [q(i, j), q(i + 1, j), q(i, j + 1), q(i + 1, j + 1)]
        iseven(i + j) ? push!(X, s) : push!(Z, s)
    end
    for j in 0:d-2
        isodd(j) && push!(X, [q(0, j), q(0, j + 1)])
        iseven(j) && push!(X, [q(d - 1, j), q(d - 1, j + 1)])
    end
    for i in 0:d-2
        iseven(i) && push!(Z, [q(i, 0), q(i + 1, 0)])
        isodd(i) && push!(Z, [q(i, d - 1), q(i + 1, d - 1)])
    end
    (_rows_to_matrix(X, d^2), _rows_to_matrix(Z, d^2))
end

function _rows_to_matrix(rows, n::Int)
    H = falses(length(rows), n)
    for (r, s) in enumerate(rows), c in s
        H[r, c] = true
    end
    H
end

function _shift(n::Int)
    S = falses(n, n)
    for i in 1:n
        S[i, mod1(i + 1, n)] = true
    end
    S
end

_gf2mul(A::AbstractMatrix{Bool}, B::AbstractMatrix{Bool}) = BitMatrix(isodd.(Int.(A) * Int.(B)))

function parity_check_matrices(m::BBCodeModel)
    x = BitMatrix(kron(_shift(m.l), Matrix{Bool}(LinearAlgebra.I, m.m, m.m)))
    y = BitMatrix(kron(Matrix{Bool}(LinearAlgebra.I, m.l, m.l), _shift(m.m)))
    mono(a, b) = _gf2mul(_matpow(x, a), _matpow(y, b))
    poly(t) = reduce((P, Q) -> P .⊻ Q, (mono(a, b) for (a, b) in t))
    A = poly(m.A); B = poly(m.B)
    (hcat(A, B), hcat(permutedims(B), permutedims(A)))
end

function _matpow(M::BitMatrix, k::Int)
    R = BitMatrix(Matrix{Bool}(LinearAlgebra.I, size(M)...))
    for _ in 1:k
        R = _gf2mul(R, M)
    end
    R
end

function parity_check_matrices(m::GenericCSSModel)
    S = QuantumClifford.ECC.parity_checks(m.code)
    n::Int = QuantumClifford.nqubits(S)
    rowsx = Vector{Vector{Int}}(); rowsz = Vector{Vector{Int}}()
    for r in 1:length(S)
        p = S[r]
        xs = [j for j in 1:n if p[j][1]]; zs = [j for j in 1:n if p[j][2]]
        isempty(xs) || push!(rowsx, xs)
        isempty(zs) || push!(rowsz, zs)
    end
    (_rows_to_matrix(rowsx, n), _rows_to_matrix(rowsz, n))
end

"""
    gf2_rank(H::AbstractMatrix{Bool}) -> Int

Rank of a matrix over GF(2).
"""
function gf2_rank(H::AbstractMatrix{Bool})
    M = BitMatrix(H); r = 0
    nr, nc = size(M)
    for c in 1:nc
        piv = findfirst(i -> M[i, c], r+1:nr)
        piv === nothing && continue
        piv += r
        r += 1
        M[r, :], M[piv, :] = M[piv, :], M[r, :]
        for i in 1:nr
            (i != r && M[i, c]) && (M[i, :] .⊻= M[r, :])
        end
        r == nr && break
    end
    r
end
