# Exact arithmetic in the cyclotomic field Q(ζ₈), ζ₈ = e^{iπ/4}.
#
# Every amplitude in a magic-state gadget proof lies in Q(ζ₈): projectors (factors
# of 1/2), Paulis, phases diag(1, ζ^k), and the magic state (|0⟩ + ζ^k|1⟩)/√2 with
# 1/√2 = (ζ − ζ³)/2. Doing the proof in this field makes it exact rather than a
# floating-point check.
#
# An element is a₀ + a₁ζ + a₂ζ² + a₃ζ³ with rational aᵢ, using ζ⁴ = −1.

"""
    Cyc8(a0, a1, a2, a3)

The element `a0 + a1·ζ + a2·ζ² + a3·ζ³` of Q(ζ₈), with `ζ = e^{iπ/4}` and exact
rational coefficients.
"""
struct Cyc8
    c::NTuple{4,Rational{BigInt}}
    Cyc8(c::NTuple{4,Rational{BigInt}}) = new(c)
end

Cyc8(a0, a1, a2, a3) = Cyc8((Rational{BigInt}(a0), Rational{BigInt}(a1), Rational{BigInt}(a2), Rational{BigInt}(a3)))
Cyc8(a::Integer) = Cyc8(a, 0, 0, 0)
Cyc8(a::Rational) = Cyc8(a, 0, 0, 0)

Base.zero(::Type{Cyc8}) = Cyc8(0)
Base.one(::Type{Cyc8}) = Cyc8(1)
Base.zero(::Cyc8) = Cyc8(0)
Base.one(::Cyc8) = Cyc8(1)
Base.:(==)(a::Cyc8, b::Cyc8) = a.c == b.c
Base.hash(a::Cyc8, h::UInt) = hash(a.c, hash(:Cyc8, h))
Base.iszero(a::Cyc8) = all(iszero, a.c)
Base.:+(a::Cyc8, b::Cyc8) = Cyc8(a.c .+ b.c)
Base.:-(a::Cyc8, b::Cyc8) = Cyc8(a.c .- b.c)
Base.:-(a::Cyc8) = Cyc8(.-(a.c))

function Base.:*(a::Cyc8, b::Cyc8)
    r = zeros(Rational{BigInt}, 4)
    for i in 0:3, j in 0:3
        t = a.c[i + 1] * b.c[j + 1]
        iszero(t) && continue
        k = i + j
        if k < 4
            r[k + 1] += t
        else
            r[k - 3] -= t          # ζ^k = −ζ^{k−4}
        end
    end
    Cyc8(r[1], r[2], r[3], r[4])
end

Base.:*(a::Cyc8, b::Union{Integer,Rational}) = Cyc8(a.c .* Rational{BigInt}(b))
Base.:*(b::Union{Integer,Rational}, a::Cyc8) = a * b

"""`ζ₈^k` for any integer `k`."""
function zeta(k::Integer)
    k = mod(k, 8)
    s = k >= 4 ? -1 : 1
    k = mod(k, 4)
    Cyc8(ntuple(i -> Rational{BigInt}(i == k + 1 ? s : 0), 4))
end

"""The Galois automorphism ζ ↦ ζ^m (m odd)."""
function galois(a::Cyc8, m::Int)
    r = Cyc8(0)
    for i in 0:3
        iszero(a.c[i + 1]) && continue
        r = r + a.c[i + 1] * zeta(m * i)
    end
    r
end

"""Complex conjugate (the automorphism ζ ↦ ζ⁻¹ = ζ⁷)."""
Base.conj(a::Cyc8) = galois(a, 7)

"""Field norm, a rational number: the product of the four Galois conjugates."""
function fieldnorm(a::Cyc8)
    p = a * galois(a, 3) * galois(a, 5) * galois(a, 7)
    all(iszero, p.c[2:4]) || error("internal: norm is not rational")
    p.c[1]
end

function Base.inv(a::Cyc8)
    iszero(a) && throw(DivideError())
    num = galois(a, 3) * galois(a, 5) * galois(a, 7)
    num * inv(fieldnorm(a))
end

Base.:/(a::Cyc8, b::Cyc8) = a * inv(b)

"""`|a|²` as an element of the field (it is real)."""
abs2c(a::Cyc8) = a * conj(a)

"""`1/√2 = (ζ − ζ³)/2`."""
const INV_SQRT2 = Cyc8(0, 1 // 2, 0, -1 // 2)

"""
    zeta_exponent(a::Cyc8) -> Union{Int,Nothing}

The `k ∈ 0:7` with `a == ζ^k`, or `nothing` if `a` is not an 8th root of unity.
"""
function zeta_exponent(a::Cyc8)
    for k in 0:7
        a == zeta(k) && return k
    end
    nothing
end

Base.show(io::IO, a::Cyc8) = print(io, "Cyc8", string.(a.c))

# Small exact matrix helpers (dense, Cyc8 entries).
function _cmatmul(A::Matrix{Cyc8}, B::Matrix{Cyc8})
    size(A, 2) == size(B, 1) || throw(DimensionMismatch())
    C = fill(Cyc8(0), size(A, 1), size(B, 2))
    for i in axes(A, 1), j in axes(B, 2)
        s = Cyc8(0)
        for k in axes(A, 2)
            (iszero(A[i, k]) || iszero(B[k, j])) && continue
            s = s + A[i, k] * B[k, j]
        end
        C[i, j] = s
    end
    C
end

_cidentity(d::Int) = [i == j ? Cyc8(1) : Cyc8(0) for i in 1:d, j in 1:d]
_cadjoint(A::Matrix{Cyc8}) = [conj(A[j, i]) for i in axes(A, 2), j in axes(A, 1)]
