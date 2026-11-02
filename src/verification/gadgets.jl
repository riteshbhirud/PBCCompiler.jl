# V4: gadget recognition (via GadgetShapes) and *semantic* proof ("lifting").
#
# For each recognized gadget with joint measurement (−1)^s·P_d⊗σ_a, magic
# measurement (−1)^{s'}τ_a and corrections exp(-iπθᵢ Qᵢ) conditioned on the two
# outcomes:
#
# 1. Let P = (−1)^s P_d (the signed data part). Synthesize a Clifford V with
#    V P V† = +Z₁, as an explicit rotation sequence, and check it exactly.
# 2. Conjugate the gadget by V: the joint measurement becomes Z₁⊗σ, and each
#    correction Pauli Qᵢ must map to εᵢZ₁ (otherwise the gadget does not rotate
#    about P and is rejected).
# 3. The reduced gadget acts on data qubit 1 and the magic qubit only. Prove it in
#    exact Q(ζ₈) arithmetic over all four outcome branches: every branch's data Kraus
#    operator must be (scalar)·Z^{r_b}·exp(-iπθZ) with one θ for all branches, the
#    branch weights must sum to the identity, and r_b must be affine in (b1, b2).
# 4. Emit exp(-iπθ P) and the residual Pauli frame P^{r(b1,b2)}. The gadget's own
#    outcome variables stay as free variables of that frame, and are otherwise
#    hidden.
#
# Proofs are cached by (template, reduced shape, magic state). Only step 1 is per
# instance.

using ..GadgetShapes: GadgetTemplate, GADGET_TEMPLATES, GadgetMatch, find_gadgets

"""
    GadgetProof

Result of the exact proof of a reduced gadget: whether it succeeded (`ok`), the
reason when it did not, the rotation angle θ (units of π, in [0, 1/2)), the
residual frame `r = c₀ ⊕ c₁·b1 ⊕ c₂·b2` as `(c₀, c₁, c₂)`, and the
state-independent probability of each outcome branch `(b1, b2)` (in the order
00, 01, 10, 11), as exact field elements.
"""
struct GadgetProof
    ok::Bool
    reason::String
    angle::Rational{Int}
    residual::NTuple{3,Bool}
    branch_probs::NTuple{4,Cyc8}
end

_fail(reason) = GadgetProof(false, reason, 0 // 1, (false, false, false), ntuple(_ -> Cyc8(0), 4))

const _I2 = [Cyc8(1) Cyc8(0); Cyc8(0) Cyc8(1)]
const _X2 = [Cyc8(0) Cyc8(1); Cyc8(1) Cyc8(0)]
const _Y2 = [Cyc8(0) -zeta(2); zeta(2) Cyc8(0)]
const _Z2 = [Cyc8(1) Cyc8(0); Cyc8(0) Cyc8(-1)]
_letter_matrix(l::Char) = l == 'X' ? _X2 : l == 'Y' ? _Y2 : l == 'Z' ? _Z2 : _I2

# Kronecker product with `A` as the outer (most significant) factor.
function _ckron(A::Matrix{Cyc8}, B::Matrix{Cyc8})
    ra, ca = size(A); rb, cb = size(B)
    C = fill(Cyc8(0), ra * rb, ca * cb)
    for i1 in 1:ra, j1 in 1:ca, i2 in 1:rb, j2 in 1:cb
        C[(i1 - 1) * rb + i2, (j1 - 1) * cb + j2] = A[i1, j1] * B[i2, j2]
    end
    C
end

# Normalized eigenvector of the single-qubit Pauli `l` with eigenvalue (−1)^b, exactly.
function _eigvec(l::Char, b::Bool)
    if l == 'Z'
        return b ? reshape([Cyc8(0), Cyc8(1)], 2, 1) : reshape([Cyc8(1), Cyc8(0)], 2, 1)
    elseif l == 'X'
        return reshape([INV_SQRT2, b ? -INV_SQRT2 : INV_SQRT2], 2, 1)
    else # 'Y': eigenvalue +1 ↔ (1, i)/√2
        return reshape([INV_SQRT2, (b ? -zeta(2) : zeta(2)) * INV_SQRT2], 2, 1)
    end
end

_projector(O::Matrix{Cyc8}, b::Bool) = [(i == j ? Cyc8(1 // 2) : Cyc8(0)) + (b ? -1 // 2 : 1 // 2) * O[i, j] for i in axes(O, 1), j in axes(O, 2)]

"""
    prove_reduced_gadget(joint_letter, magic_letter, magic_sign, corrections, magic_k) -> GadgetProof

Exact proof of a gadget reduced to one data qubit (`Z ≡ P`) and one magic qubit,
prepared in `(|0⟩ + ζ₈^magic_k |1⟩)/√2`. The joint measurement is `σ⊗Z`, with
`σ = joint_letter` on the magic qubit. The magic measurement is
`(−1)^magic_sign · magic_letter`. `corrections` lists `(θ, role)`: apply
`exp(-iπθ Z)` to the data qubit if outcome `role` (1 = joint, 2 = magic) is 1.
"""
function prove_reduced_gadget(joint_letter::Char, magic_letter::Char, magic_sign::Bool,
                              corrections::Vector{Tuple{Rational{Int},Int}}, magic_k::Int)
    for (θ, _) in corrections
        isinteger(8θ) || return _fail("correction angle $(θ)π is not a multiple of π/8")
    end
    J = _ckron(_letter_matrix(joint_letter), _Z2)                       # index = d + 2m
    Tm = _ckron(magic_sign ? -1 * _letter_matrix(magic_letter) : _letter_matrix(magic_letter), _I2)
    ketm = reshape([INV_SQRT2, zeta(magic_k) * INV_SQRT2], 2, 1)
    M0 = _ckron(ketm, _I2)                                             # |d⟩ ↦ |m⟩⊗|d⟩
    θs = Dict{Tuple{Bool,Bool},Rational{Int}}()
    probs = Cyc8[]
    sumα = Cyc8(0); sumβ = Cyc8(0)
    for b1 in (false, true), b2 in (false, true)
        M = _cmatmul(_projector(J, b1), M0)
        M = _cmatmul(_projector(Tm, b2), M)
        for (θ, role) in corrections
            (role == 1 ? b1 : b2) || continue
            D = [Cyc8(1) Cyc8(0); Cyc8(0) zeta(Int(8θ))]             # ∝ exp(-iπθZ)
            M = _cmatmul(_ckron(_I2, D), M)
        end
        e = _eigvec(magic_letter, b2 ⊻ magic_sign)
        Kd = _cmatmul(_ckron(_cadjoint(e), _I2), M)
        _cmatmul(_ckron(e, _I2), Kd) == M || return _fail("internal: magic qubit not left in the measured eigenstate")
        (iszero(Kd[1, 2]) && iszero(Kd[2, 1])) || return _fail("branch ($(Int(b1)),$(Int(b2))) is not diagonal in the eigenbasis of P")
        α, β = Kd[1, 1], Kd[2, 2]
        pα, pβ = abs2c(α), abs2c(β)
        sumα = sumα + pα; sumβ = sumβ + pβ
        push!(probs, pα)
        (iszero(α) && iszero(β)) && continue
        pα == pβ || return _fail("branch ($(Int(b1)),$(Int(b2))) is not unitary on the data (|⟨0|K|0⟩|² ≠ |⟨1|K|1⟩|²), e.g. a wrong magic state or measurement basis")
        j = zeta_exponent(β / α)
        j === nothing && return _fail("branch ($(Int(b1)),$(Int(b2))) applies a relative phase that is not a multiple of π/4")
        θs[(b1, b2)] = _normalize_angle(j // 8)                            # diag(1, ζʲ) ∝ exp(-iπ(j/8)Z)
    end
    (sumα == Cyc8(1) && sumβ == Cyc8(1)) || return _fail("the branch Kraus operators do not sum to a trace-preserving map")
    isempty(θs) && return _fail("every branch has probability zero")
    ks = sort!(collect(keys(θs)))
    θ0 = mod(θs[ks[1]], 1 // 2)
    for k in ks
        iszero(mod(θs[k] - θ0, 1 // 2)) || return _fail("branches implement different rotations ($(θs[ks[1]])π vs $(θs[k])π): an outcome-dependent non-Pauli correction is missing or wrong")
    end
    r = Dict(k => isodd(Int(mod(θs[k] - θ0, 1) * 2)) for k in ks)
    fit = nothing
    for c0 in (false, true), c1 in (false, true), c2 in (false, true)
        if all(r[k] == (c0 ⊻ (c1 & k[1]) ⊻ (c2 & k[2])) for k in ks)
            fit = (c0, c1, c2)
            break
        end
    end
    fit === nothing && return _fail("the residual Pauli frame is not an affine function of the gadget outcomes")
    GadgetProof(true, "", θ0, fit, Tuple(probs))
end

const _GADGET_PROOF_CACHE = Dict{Tuple{Symbol,Char,Char,Bool,Vector{Tuple{Rational{Int},Int}},Int},GadgetProof}()
const _GADGET_PROOF_LOCK = ReentrantLock()

function _cached_proof(key)
    lock(_GADGET_PROOF_LOCK) do
        hit = haskey(_GADGET_PROOF_CACHE, key)
        proof = hit ? _GADGET_PROOF_CACHE[key] : (_GADGET_PROOF_CACHE[key] = prove_reduced_gadget(key[2], key[3], key[4], key[5], key[6]))
        (proof, hit)
    end
end

"""
    synthesize_z1_clifford(p::PauliOperator) -> (rotations, map)

A Clifford `V` with `V p V† = +Z₁`, for a Hermitian non-identity Pauli `p`. It is
returned as a list of rotations `(Q, k)` meaning `exp(-ikπ/4·Q)` (applied in
order), together with its forward `PauliMap`. The construction maps every X or Y
letter to Z, folds the Z support onto one qubit with CNOTs, moves that qubit to
position 1 with two CNOTs, and fixes the sign with an X.
"""
function synthesize_z1_clifford(p::PauliOperator)
    n = nqubits(p)
    (_ishermitian(p) && !_isidentity(p)) || throw(ArgumentError("need a Hermitian non-identity Pauli"))
    rots = Tuple{PauliT,Int}[]
    fwd = identity_map(n); inv = identity_map(n)
    single(l, j) = (q = zero(PauliOperator, n); q[j] = (l == 'X' || l == 'Y', l == 'Z' || l == 'Y'); q)
    function rot!(q, k)
        push!(rots, (q, k)); rotate!(fwd, inv, q, k)
    end
    # CNOT(c, t) ∝ (Z_c)_{π/4} (X_t)_{π/4} (Z_c X_t)_{−π/4}
    cnot!(c, t) = (rot!(single('Z', c), 1); rot!(single('X', t), 1); rot!(single('Z', c) * single('X', t), 3))
    for j in 1:n
        x, z = p[j]
        x || continue
        rot!(z ? single('X', j) : single('Y', j), 1)       # Y ↦ ±Z via X-rotation, X ↦ ±Z via Y-rotation
    end
    cur = apply_map(fwd, p)
    S = [j for j in 1:n if cur[j][2]]
    t = S[1]
    for a in S[2:end]
        cnot!(a, t)                                         # Z_a Z_t ↦ Z_t
    end
    if t != 1
        cnot!(1, t)                                         # Z_t ↦ Z₁Z_t
        cnot!(t, 1)                                         # Z₁Z_t ↦ Z₁
    end
    cur = apply_map(fwd, p)
    cur.phase[] == 0x02 && rot!(single('X', 1), 2)           # −Z₁ ↦ +Z₁
    apply_map(fwd, p) == single('Z', 1) || error("internal: Clifford synthesis failed for $p")
    (rots, fwd)
end

"""
    LiftedGadget

A gadget whose action was proved: it implements `exp(-iπ·angle·data_pauli)`
followed by the residual frame `data_pauli^residual` (`residual` is affine in the
gadget's outcome variables). `clifford_V` is the rotation sequence of the
synthesized Clifford `V` with `V·data_pauli·V† = Z₁`, and `proof_key` is the
key of the cached exact proof.
"""
struct LiftedGadget
    match::GadgetMatch
    data_pauli::PauliT
    angle::Rational{Int}
    residual::AffineBit
    clifford_V::Vector{Tuple{PauliT,Int}}
    proof_key::Tuple{Symbol,Char,Char,Bool,Vector{Tuple{Rational{Int},Int}},Int}
    proof::GadgetProof
end

"""A recognized gadget that could not be proved, with the reason."""
struct GadgetFailure
    match::GadgetMatch
    reason::String
end

"""
    LiftingResult

Output of [`lift_gadgets`](@ref): the proved gadgets, the gadgets that were
recognized but failed their proof, and cache statistics.
"""
struct LiftingResult
    n_data::Int
    lifted::Vector{LiftedGadget}
    failed::Vector{GadgetFailure}
    cache_hits::Int
    cache_misses::Int
end

_letter_of(p::PauliOperator, j::Int) = (x = p[j][1]; z = p[j][2]; x ? (z ? 'Y' : 'X') : (z ? 'Z' : 'I'))

"""
    lift_gadgets(c, n_data; magic_k = 1, templates = GADGET_TEMPLATES) -> LiftingResult

Recognize the magic-state gadgets of a compiled program (qubits `> n_data` are
magic qubits in `(|0⟩ + ζ₈^magic_k|1⟩)/√2`; `magic_k = 1` is PBCCompiler's `T|+⟩`)
and prove each one exactly. See the file header for the method.
"""
function lift_gadgets(c::AbstractVector{CircuitOp.Type}, n_data::Int; magic_k::Int = 1, templates = GADGET_TEMPLATES)
    lifted = LiftedGadget[]; failed = GadgetFailure[]
    hits = 0; misses = 0
    for g in find_gadgets(c, n_data; templates)
        r = _lift_one(c, g, n_data, magic_k)
        if r isa GadgetFailure
            push!(failed, r)
        else
            lg, hit = r
            push!(lifted, lg)
            hit ? (hits += 1) : (misses += 1)
        end
    end
    LiftingResult(n_data, lifted, failed, hits, misses)
end

function _lift_one(c, g::GadgetMatch, n_data::Int, magic_k::Int)
    g.joint_bit == g.magic_bit && return GadgetFailure(g, "both gadget measurements write bit $(g.joint_bit)")
    jp, jq = _pauli_and_qubits(c[g.joint_index])
    _ishermitian(jp) || return GadgetFailure(g, "joint measurement Pauli is not Hermitian")
    apos = findfirst(==(g.magic_qubit), jq)
    apos === nothing && return GadgetFailure(g, "internal: magic qubit not in the joint measurement")
    joint_letter = _letter_of(jp, apos)
    # signed data part: everything except the magic position, carrying the overall sign
    dpos = [j for j in eachindex(jq) if j != apos]
    dp = PauliOperator(jp.phase[], [jp[j][1] for j in dpos], [jp[j][2] for j in dpos])
    P = embed(n_data, jq[dpos], dp)
    mp, _ = _pauli_and_qubits(c[g.magic_index])
    _ishermitian(mp) || return GadgetFailure(g, "magic measurement Pauli is not Hermitian")
    magic_letter = _letter_of(mp, 1)
    magic_sign = mp.phase[] == 0x02
    rots, Vmap = synthesize_z1_clifford(P)
    Z1 = (z = zero(PauliOperator, n_data); z[1] = (false, true); z)
    corrections = Tuple{Rational{Int},Int}[]
    for ci in g.correction_indices
        inner, b = _bitcond_parts(c[ci])
        k = _opkind(inner)
        q, qs = _pauli_and_qubits(inner)
        Q = _embed_checked(q, qs, n_data, "gadget correction")
        Q isa String && return GadgetFailure(g, Q)
        VQ = apply_map(Vmap, Q)
        ε = if VQ == Z1
            1
        elseif VQ == -Z1
            -1
        else
            return GadgetFailure(g, "correction at op $ci rotates about $q on $qs, which is not ±(the measured data Pauli $P)")
        end
        role = b == g.joint_bit ? 1 : 2
        push!(corrections, (ε * _rotation_angle(k), role))
    end
    key = (g.template, joint_letter, magic_letter, magic_sign, corrections, magic_k)
    proof, hit = _cached_proof(key)
    proof.ok || return GadgetFailure(g, proof.reason)
    c0, c1, c2 = proof.residual
    vars = Int[]
    c1 && push!(vars, g.joint_index)
    c2 && push!(vars, g.magic_index)
    (LiftedGadget(g, P, proof.angle, AffineBit(c0, vars), rots, key, proof), hit)
end

"""
    lower_lifted(c, lr::LiftingResult) -> LoweredProgram

Lower a compiled program on its data qubits, replacing each proved gadget by the
rotation and residual frame it implements. The gadget's two measurements still
count as writes of their bits, with their own outcome variables. Any remaining op
on a magic qubit (for example inside a gadget that failed its proof) is reported
as unsupported.
"""
function lower_lifted(c::AbstractVector{CircuitOp.Type}, lr::LiftingResult)
    repl = Dict{Int,Tuple{Int,Vector{VOp},Vector{Tuple{Int,Int}}}}()
    for lg in lr.lifted
        g = lg.match
        vops = VOp[_vrot(lg.data_pauli, lg.angle, AffineBit(true), g.joint_index)]
        iszero(lg.residual) || push!(vops, _vrot(lg.data_pauli, 1 // 2, lg.residual, g.joint_index))
        repl[g.joint_index] = (last(g.ops), vops, [(g.joint_bit, g.joint_index), (g.magic_bit, g.magic_index)])
    end
    lower(c, lr.n_data; replacements = repl)
end
