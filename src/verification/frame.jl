# V2: symbolic Pauli frame.
#
# A Pauli operator whose X and Z bits on each qubit are affine functions (over F2)
# of the measurement-outcome variables. It represents all 2^m outcome branches at
# once.
#
# The frame is tracked up to an outcome-dependent global phase. That phase is
# unobservable within a branch, and is in general not affine in the outcomes: the
# sign of a conjugated frame can be quadratic, e.g. S maps X→Y and Y→−X. The
# observable consequences of the frame (sign flips of later measurements and
# rotations) are commutation parities, which are affine; see
# `commutation_parity`.

"""
    SymbolicPauli(n::Int)

An `n`-qubit Pauli operator whose X and Z bits are [`AffineBit`](@ref)s, i.e.
affine functions over F2 of measurement-outcome variables. Constructed as the
identity. It is tracked up to a global phase per outcome branch.

Operations: [`mul_conditional!`](@ref) (multiply by `Q^c` for an affine
condition `c`), [`conjugate!`](@ref) (conjugate by a Clifford),
[`commutation_parity`](@ref), [`evaluate`](@ref), and `==` (equality as
functions of the outcomes).
"""
struct SymbolicPauli
    x::Vector{AffineBit}
    z::Vector{AffineBit}
end

SymbolicPauli(n::Int) = SymbolicPauli([AffineBit(false) for _ in 1:n], [AffineBit(false) for _ in 1:n])

QuantumClifford.nqubits(f::SymbolicPauli) = length(f.x)
Base.copy(f::SymbolicPauli) = SymbolicPauli(copy(f.x), copy(f.z))
Base.:(==)(a::SymbolicPauli, b::SymbolicPauli) = a.x == b.x && a.z == b.z
Base.hash(a::SymbolicPauli, h::UInt) = hash(a.z, hash(a.x, hash(:SymbolicPauli, h)))

"""Whether the frame is the identity for every outcome assignment."""
Base.isone(f::SymbolicPauli) = all(iszero, f.x) && all(iszero, f.z)

"""
    mul_conditional!(f::SymbolicPauli, q::PauliOperator, c::AffineBit) -> f

Multiply the frame by `q^c`: apply `q` in exactly the branches where `c` evaluates
to 1. The phase of `q` is ignored, since the frame is tracked up to a global
phase.
"""
function mul_conditional!(f::SymbolicPauli, q::PauliOperator, c::AffineBit)
    nqubits(q) == length(f.x) || throw(DimensionMismatch("frame has $(length(f.x)) qubits, Pauli has $(nqubits(q))"))
    iszero(c) && return f
    @inbounds for j in 1:length(f.x)
        xj, zj = q[j]
        xj && (f.x[j] = f.x[j] ⊻ c)
        zj && (f.z[j] = f.z[j] ⊻ c)
    end
    f
end

mul_conditional!(f::SymbolicPauli, q::PauliOperator, c::Bool) = mul_conditional!(f, q, AffineBit(c))

"""
    mul_frame!(f::SymbolicPauli, g::SymbolicPauli) -> f

Multiply the frame `f` by the frame `g` in place (bitwise XOR, up to phase).
"""
function mul_frame!(f::SymbolicPauli, g::SymbolicPauli)
    length(f.x) == length(g.x) || throw(DimensionMismatch("frames have different sizes"))
    @inbounds for j in 1:length(f.x)
        f.x[j] = f.x[j] ⊻ g.x[j]
        f.z[j] = f.z[j] ⊻ g.z[j]
    end
    f
end

"""
    commutation_parity(f::SymbolicPauli, p::PauliOperator) -> AffineBit

The affine function that is 1 exactly in the branches where the frame
anticommutes with `p`, i.e. the symplectic product of the frame's bits with `p`.
"""
function commutation_parity(f::SymbolicPauli, p::PauliOperator)
    nqubits(p) == length(f.x) || throw(DimensionMismatch("frame has $(length(f.x)) qubits, Pauli has $(nqubits(p))"))
    r = AffineBit(false)
    @inbounds for j in 1:length(f.x)
        xj, zj = p[j]
        zj && (r = r ⊻ f.x[j])
        xj && (r = r ⊻ f.z[j])
    end
    r
end

"""
    conjugate!(f::SymbolicPauli, images::AbstractVector{<:PauliOperator}) -> f
    conjugate!(f::SymbolicPauli, c::CliffordOperator) -> f

Replace the frame `F` by `C F C†` for the (unconditional) Clifford `C`, given
either as a `CliffordOperator` or as the vector of images `[C X₁ C†, …, C Xₙ C†,
C Z₁ C†, …, C Zₙ C†]`. The images' bits enter linearly, so the result stays affine
in the outcomes.
"""
function conjugate!(f::SymbolicPauli, images::AbstractVector{<:PauliOperator})
    n = length(f.x)
    length(images) == 2n || throw(DimensionMismatch("need 2n = $(2n) images, got $(length(images))"))
    nx = [AffineBit(false) for _ in 1:n]
    nz = [AffineBit(false) for _ in 1:n]
    @inbounds for j in 1:n
        for (src, img) in ((f.x[j], images[j]), (f.z[j], images[n + j]))
            iszero(src) && continue
            for k in 1:n
                xk, zk = img[k]
                xk && (nx[k] = nx[k] ⊻ src)
                zk && (nz[k] = nz[k] ⊻ src)
            end
        end
    end
    copyto!(f.x, nx)
    copyto!(f.z, nz)
    f
end

conjugate!(f::SymbolicPauli, c::CliffordOperator) = conjugate!(f, [c[i] for i in 1:2 * length(f.x)])

"""
    evaluate(f::SymbolicPauli, assignment) -> PauliOperator

The concrete Pauli (with phase `+1`) that the frame equals when the outcome
variables take the values in `assignment`.
"""
function evaluate(f::SymbolicPauli, assignment)
    n = length(f.x)
    PauliOperator(0x00, [evaluate(f.x[j], assignment) for j in 1:n], [evaluate(f.z[j], assignment) for j in 1:n])
end

"""
    substitute(f::SymbolicPauli, g) -> SymbolicPauli

Replace every outcome variable `v` by the affine function `g(v)` in all bits.
"""
substitute(f::SymbolicPauli, g) = SymbolicPauli([substitute(a, g) for a in f.x], [substitute(a, g) for a in f.z])

"""All outcome variables the frame depends on, sorted."""
function frame_variables(f::SymbolicPauli)
    s = Set{Int}()
    for a in f.x, v in a.vars; push!(s, v); end
    for a in f.z, v in a.vars; push!(s, v); end
    sort!(collect(s))
end

function Base.show(io::IO, f::SymbolicPauli)
    print(io, "SymbolicPauli(")
    first = true
    for j in 1:length(f.x)
        (iszero(f.x[j]) && iszero(f.z[j])) && continue
        first || print(io, ", ")
        first = false
        print(io, "q$j: X^(", f.x[j], ") Z^(", f.z[j], ")")
    end
    first && print(io, "I")
    print(io, ")")
end
