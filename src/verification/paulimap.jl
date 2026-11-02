"""
    PauliMap

The conjugation action of a Clifford unitary, stored as the images of the 2n
generators: `img[j]` is the image of `Xⱼ` and `img[n+j]` the image of `Zⱼ`, signs
included. It determines the Clifford up to a global phase.
"""
struct PauliMap
    n::Int
    img::Vector{PauliT}
end

"""The identity map on `n` qubits."""
function identity_map(n::Int)
    img = Vector{PauliT}(undef, 2n)
    for j in 1:n
        x = zero(PauliOperator, n); x[j] = (true, false); img[j] = x
        z = zero(PauliOperator, n); z[j] = (false, true); img[n + j] = z
    end
    PauliMap(n, img)
end

Base.copy(m::PauliMap) = PauliMap(m.n, [copy(p) for p in m.img])
Base.:(==)(a::PauliMap, b::PauliMap) = a.n == b.n && a.img == b.img

_set_phase(p::PauliT, ph::Integer) = (q = copy(p); q.phase[] = UInt8(mod(ph, 4)); q)

"""
    apply_map(m::PauliMap, p::PauliOperator) -> PauliOperator

The image of `p` under the Clifford map, including its phase.
"""
function apply_map(m::PauliMap, p::PauliOperator)
    n = m.n
    nqubits(p) == n || throw(DimensionMismatch("map acts on $n qubits, Pauli has $(nqubits(p))"))
    r = zero(PauliOperator, n)
    extra = Int(p.phase[])
    @inbounds for j in 1:n
        xj, zj = p[j]
        if xj
            r = r * m.img[j]
        end
        if zj
            r = r * m.img[n + j]
        end
        (xj && zj) && (extra += 1)
    end
    r.phase[] = UInt8(mod(Int(r.phase[]) + extra, 4))
    r
end

# generator idx is Xⱼ for idx ≤ n, Zⱼ otherwise
function _generator_anticommutes(n::Int, idx::Int, q::PauliOperator)
    if idx <= n
        return q[idx][2]          # Xⱼ anticommutes with q iff q has Z on j
    else
        return q[idx - n][1]      # Zⱼ anticommutes with q iff q has X on j
    end
end

# if {r, q} = 0 then G r G† = r exp(2iθq) for G = exp(-iθq)
"""
    rotate!(fwd::PauliMap, inv::PauliMap, q::PauliOperator, k::Integer)

Update the pair of maps for `F ← exp(-i k π/4 · q) F`, with `k ∈ {1, 2, 3}` and
`q` Hermitian (its sign is respected).
"""
function rotate!(fwd::PauliMap, inv::PauliMap, q::PauliOperator, k::Integer)
    (k == 1 || k == 2 || k == 3) || throw(ArgumentError("k must be 1, 2 or 3 (got $k)"))
    n = fwd.n
    invq = apply_map(inv, q)
    @inbounds for idx in 1:2n
        r = fwd.img[idx]
        if comm(r, q) == 0x01
            fwd.img[idx] = k == 1 ? 1im * (r * q) : k == 3 ? -1im * (r * q) : -r
        end
        if _generator_anticommutes(n, idx, q)
            e = inv.img[idx]
            inv.img[idx] = k == 1 ? -1im * (e * invq) : k == 3 ? 1im * (e * invq) : -e
        end
    end
    nothing
end

"""
    compose(a::PauliMap, b::PauliMap) -> PauliMap

The map `P ↦ a(b(P))`.
"""
compose(a::PauliMap, b::PauliMap) = PauliMap(a.n, [apply_map(a, p) for p in b.img])

"""
    to_clifford(m::PauliMap) -> CliffordOperator

The QuantumClifford `CliffordOperator` with the same generator images.
"""
to_clifford(m::PauliMap) = CliffordOperator(m.img)

"""
    pauli_sign_conjugate_fwd!(fwd::PauliMap, f::PauliOperator)
    pauli_sign_conjugate_inv!(inv::PauliMap, f::PauliOperator)

Update the forward map `P ↦ F P F†` (resp. the inverse map `P ↦ F† P F`) for
`F ← f F` with `f` a Pauli.
"""
function pauli_sign_conjugate_fwd!(m::PauliMap, f::PauliOperator)
    for idx in 1:2m.n
        comm(m.img[idx], f) == 0x01 && (m.img[idx] = -m.img[idx])
    end
    m
end

function pauli_sign_conjugate_inv!(m::PauliMap, f::PauliOperator)
    for idx in 1:2m.n
        _generator_anticommutes(m.n, idx, f) && (m.img[idx] = -m.img[idx])
    end
    m
end
