"""
    StabilizerVerdict

Result of [`stabilizer_oracle`](@ref). `applicable` is false when the oracle
could not run on these programs, `refuted` is true when it found a difference,
`trials` is the number of trials run and `max_deviation` the largest probability
difference seen. When refuted, `witness_state` (a stabilizer state on the data
qubits), `witness_record`, `p_input` and `p_compiled` describe the difference.
"""
struct StabilizerVerdict
    applicable::Bool
    refuted::Bool
    trials::Int
    max_deviation::Float64
    witness_state::Union{Nothing,QuantumClifford.Stabilizer}
    witness_record::Dict{Int,Bool}
    p_input::Float64
    p_compiled::Float64
    note::String
end

_gs_notapplicable(note) = StabilizerVerdict(false, false, 0, NaN, nothing, Dict{Int,Bool}(), NaN, NaN, note)

function _gs_clifford_rotation(P::PauliT, θ::Rational{Int}, N::Int)
    fwd = identity_map(N); inv = identity_map(N)
    rotate!(fwd, inv, P, θ == 1 // 4 ? 1 : θ == -1 // 4 ? 3 : 2)
    to_clifford(fwd)
end

# QuantumClifford 0.11.8 throws a DomainError in `apply!(::GeneralizedStabilizer, ::PauliChannel)`
# when a channel Pauli has a -1 sign, so channel Paulis are kept unsigned and the sign goes into the weights.
_unsigned(P::PauliT) = (P.phase[] == 0x00 ? (P, 1.0) : (_set_phase(P, 0), P.phase[] == 0x02 ? -1.0 : error("non-Hermitian Pauli")))

function _gs_apply_rotation!(st, P::PauliT, θ::Rational{Int}, N::Int)
    t = _normalize_angle(θ)
    (iszero(t) || _isidentity(P)) && return st
    if denominator(t) <= 4
        QuantumClifford.apply!(st, _gs_clifford_rotation(P, t, N))
    else
        P0, sg = _unsigned(P)
        ch = QuantumClifford.UnitaryPauliChannel([zero(PauliOperator, N), P0], ComplexF64[cospi(Float64(t)), -im * sg * sinpi(Float64(t))])
        QuantumClifford.apply!(st, ch; prune_threshold = 0.0)
    end
    st
end

function _gs_apply_op!(st, op, N::Int, reg::Dict{Int,Bool})
    k = _opkind(op)
    if _isrotation(k)
        p, q = _pauli_and_qubits(op)
        _gs_apply_rotation!(st, embed(N, q, p), _rotation_angle(k), N)
    elseif k == OP_PAULICOND
        cp, cq, tp, tq = _paulicond_parts(op)
        P1 = embed(N, cq, cp); P2 = embed(N, tq, tp)
        _gs_apply_rotation!(st, P1, 1 // 4, N)
        _gs_apply_rotation!(st, P2, 1 // 4, N)
        _gs_apply_rotation!(st, P1 * P2, -1 // 4, N)
    elseif k == OP_BITCOND
        inner, b = _bitcond_parts(op)
        haskey(reg, b) || throw(ArgumentError("bit $b read before write"))
        _opkind(inner) == OP_MEASUREMENT && throw(ArgumentError("conditional measurement"))
        reg[b] && _gs_apply_op!(st, inner, N, reg)
    else
        throw(ArgumentError("unsupported op for the stabilizer oracle: $(_kindname(k))"))
    end
    st
end

# Unnormalized projection; its trace is the outcome probability. The internal QuantumClifford
# routine updates the tableau, whereas applying (I ± P)/2 as a channel would grow χ up to 4x per measurement.
function _gs_project(st, P::PauliT, outcome::Bool, N::Int)
    P0, sg = _unsigned(P)
    res = (outcome ⊻ (sg < 0)) ? 1 : 0
    first(QuantumClifford.PauliChannelNonClifford._projectrand_notnorm(copy(st), P0, res))
end

_gs_trace(st) = real(LinearAlgebra.tr(st))

function _gs_initial(s_data, n::Int, N::Int, magic_k::Int)
    N == n && return QuantumClifford.GeneralizedStabilizer(s_data)
    magic = one(QuantumClifford.Stabilizer, N - n; basis = :X)
    g = QuantumClifford.GeneralizedStabilizer(QuantumClifford.tensor(s_data, magic))
    for a in n+1:N
        QuantumClifford.apply!(g, embed(N, [a], QuantumClifford.pcPhase(magic_k * π / 4)); prune_threshold = 0.0)
    end
    g
end

# Returns the probability of `record`. Measurements not in `visible` are summed over;
# with `sample`, missing record bits are sampled first.
function _gs_run(c, N::Int, st0, visible::Dict{Int,Int}, flip::Dict{Int,Bool}, record::Dict{Int,Bool},
                 sample::Bool, rng, max_branches::Int)
    branches = [(st0, Dict{Int,Bool}())]
    for (i, op) in enumerate(c)
        if _opkind(op) == OP_MEASUREMENT
            p, q = _pauli_and_qubits(op)
            P = embed(N, q, p); b = _measurement_bit(op)
            if haskey(visible, i)
                vb = visible[i]
                if !haskey(record, vb)
                    sample || throw(ArgumentError("record lacks bit $vb"))
                    p1 = sum(_gs_trace(_gs_project(st, P, true, N)) for (st, _) in branches)
                    tot = sum(_gs_trace(st) for (st, _) in branches)
                    record[vb] = rand(rng) * tot < p1
                end
                val = record[vb] ⊻ get(flip, i, false)
                newb = similar(branches, 0)
                for (st, reg) in branches
                    s2 = _gs_project(st, P, val, N)
                    _gs_trace(s2) > 1e-300 || continue
                    reg2 = copy(reg); reg2[b] = val
                    push!(newb, (s2, reg2))
                end
                branches = newb
            else
                newb = similar(branches, 0)
                for (st, reg) in branches, val in (false, true)
                    s2 = _gs_project(st, P, val, N)
                    _gs_trace(s2) > 1e-300 || continue
                    reg2 = copy(reg); reg2[b] = val
                    push!(newb, (s2, reg2))
                end
                branches = newb
                length(branches) <= max_branches || throw(ArgumentError("more than $max_branches hidden branches"))
            end
            isempty(branches) && return 0.0
        else
            for (st, reg) in branches
                _gs_apply_op!(st, op, N, reg)
            end
        end
    end
    sum((_gs_trace(st) for (st, _) in branches); init = 0.0)
end

"""
    stabilizer_oracle(input, compiled; n_data, outcome_map = :identity, magic_k = 1, trials = 20,
                      rng = Random.default_rng(), max_branches = 256, tol = 1e-6) -> StabilizerVerdict

Look for a difference between two programs. For random stabilizer input states,
an outcome record is sampled from the input program and its probability is
compared under both programs. Returns a [`StabilizerVerdict`](@ref). This oracle
can only refute equivalence, never certify it.

`outcome_map` must be `:identity` or a `Dict` whose values are single compiled
bits (an `Int`, or an `AffineBit` with one variable and an optional constant).
Parity maps are not supported.
"""
function stabilizer_oracle(input::AbstractVector{CircuitOp.Type}, compiled::AbstractVector{CircuitOp.Type};
                           n_data::Int = _width(input), outcome_map = :identity, magic_k::Int = 1, trials::Int = 20,
                           rng::AbstractRNG = Random.default_rng(), max_branches::Int = 256, tol::Float64 = 1e-6,
                           visible_bits::Union{Nothing,AbstractVector{Int}} = nothing)
    n = n_data
    N = max(n, _width(compiled))
    Nin = max(n, _width(input))
    fw_in = _final_writes(input); fw_out = _final_writes(compiled)
    vbits = visible_bits === nothing ? collect(keys(fw_in)) : collect(Int, visible_bits)
    all(b -> haskey(fw_in, b), vbits) || return _gs_notapplicable("a visible bit is never written by the input")
    vis_in = Dict(fw_in[b] => b for b in vbits)
    vis_out = Dict{Int,Int}(); flip = Dict{Int,Bool}()
    for b in vbits
        e = outcome_map isa AbstractDict ? get(outcome_map, b, nothing) : affine_var(b)
        e === nothing && return _gs_notapplicable("outcome map has no entry for input bit $b")
        e isa Integer && (e = affine_var(e))
        length(e.vars) == 1 || return _gs_notapplicable("parity outcome maps are not supported by the stabilizer oracle")
        cb = e.vars[1]
        haskey(fw_out, cb) || return StabilizerVerdict(true, true, 0, Inf, nothing, Dict{Int,Bool}(), NaN, NaN,
                                                       "the compiled program never writes bit $cb, which input bit $b maps to")
        vis_out[fw_out[cb]] = b
        flip[fw_out[cb]] = e.constant
    end
    maxdev = 0.0
    for t in 1:trials
        s_data = QuantumClifford.random_stabilizer(rng, n)
        st_in = _gs_initial(s_data, n, Nin, magic_k)
        st_out = _gs_initial(s_data, n, N, magic_k)
        record = Dict{Int,Bool}()
        pin, pout = try
            (_gs_run(input, Nin, st_in, vis_in, Dict{Int,Bool}(), record, true, rng, max_branches),
             _gs_run(compiled, N, st_out, vis_out, flip, record, false, rng, max_branches))
        catch err
            err isa ArgumentError || rethrow()
            return _gs_notapplicable("not simulable: " * err.msg)
        end
        d = abs(pin - pout)
        maxdev = max(maxdev, d)
        if d > tol * max(1.0, pin, pout)
            return StabilizerVerdict(true, true, t, maxdev, s_data, record, pin, pout,
                                     "the probability of an outcome record differs for a stabilizer input state")
        end
    end
    StabilizerVerdict(true, false, trials, maxdev, nothing, Dict{Int,Bool}(), NaN, NaN,
                      "no difference found in $trials trials (not a proof)")
end
