"""
    ComparisonResult

Result of [`compare_forms`](@ref). `equal` says whether the forms match. When
they do not, `reason` explains why and `position` is the first differing item
(0 if the difference is not tied to an item). When they do, `outcome_map` maps
each input bit to an affine function of the compiled outcome variables.
"""
struct ComparisonResult
    equal::Bool
    reason::String
    position::Int
    outcome_map::Dict{Int,AffineBit}
end

_mismatch(reason, k) = ComparisonResult(false, reason, k, Dict{Int,AffineBit}())

"""
    compare_forms(a::CanonicalForm, b::CanonicalForm; outcome_map = :infer) -> ComparisonResult

Check whether the canonical forms `a` (input) and `b` (compiled) describe
equivalent programs.

- `outcome_map = :infer`: accept any affine (parity) map between outcomes and return it.
- `:identity`: each input bit must equal the compiled program's final value of the same bit.
- a `Dict{Int,AffineBit}`: maps each input bit to a parity of compiled bits (read at their final write); this map is checked.

`visible_bits` restricts the comparison to those input bits (default: every bit
the input writes).

A mismatch means the forms could not be proved equal, not that the programs are
different.
"""
function compare_forms(a::CanonicalForm, b::CanonicalForm; outcome_map = :infer,
                       visible_bits::Union{Nothing,AbstractVector{Int}} = nothing)
    a.semantics == b.semantics || throw(ArgumentError("forms use different semantics"))
    A = a.nf.items; B = b.nf.items
    vmap = Dict{Int,AffineBit}()
    translate(x::AffineBit) = all(v -> haskey(vmap, v), x.vars) ? substitute(x, v -> vmap[v]) : nothing
    for k in 1:min(length(A), length(B))
        x = A[k]; y = B[k]
        x.kind == y.kind || return _mismatch("item kinds differ ($(x.kind) vs $(y.kind))", k)
        x.pauli == y.pauli || return _mismatch("Paulis differ ($(x.pauli) vs $(y.pauli))", k)
        tc = translate(x.cond)
        tc === nothing && return _mismatch("input condition uses an unmapped outcome", k)
        tc == y.cond || return _mismatch("conditions differ ($(tc) vs $(y.cond))", k)
        ts = translate(x.sign)
        ts === nothing && return _mismatch("input sign uses an unmapped outcome", k)
        if x.kind == ITEM_ROTATION
            x.angle == y.angle || return _mismatch("rotation angles differ ($(x.angle)π vs $(y.angle)π)", k)
            ts == y.sign || return _mismatch("rotation signs differ ($(ts) vs $(y.sign))", k)
        else
            vmap[x.var] = affine_var(y.var) ⊻ ts ⊻ y.sign
        end
    end
    length(A) == length(B) || return _mismatch("different numbers of items ($(length(A)) vs $(length(B)))",
                                                 min(length(A), length(B)) + 1)
    omap = Dict{Int,AffineBit}()
    for bit in (visible_bits === nothing ? keys(a.nf.final_var) : visible_bits)
        haskey(a.nf.final_var, bit) || return _mismatch("visible bit $bit is never written by the input", 0)
        v = a.nf.final_var[bit]
        haskey(vmap, v) || return _mismatch("input bit $bit is not produced by a matched measurement", 0)
        omap[bit] = vmap[v]
    end
    if outcome_map === :identity
        for (bit, e) in omap
            haskey(b.nf.final_var, bit) || return _mismatch("compiled program never writes bit $bit", 0)
            e == affine_var(b.nf.final_var[bit]) ||
                return _mismatch("bit $bit equals $(e) in the compiled program, not its own final value (identity outcome map required)", 0)
        end
    elseif outcome_map isa AbstractDict
        for (bit, want) in outcome_map
            haskey(omap, bit) || return _mismatch("explicit outcome map mentions unknown input bit $bit", 0)
            w = _bits_to_vars(want, b.nf)
            w === nothing && return _mismatch("explicit outcome map for bit $bit uses a compiled bit that is never written", 0)
            omap[bit] == w || return _mismatch("bit $bit equals $(omap[bit]), not the given $(want)", 0)
        end
    elseif outcome_map !== :infer
        throw(ArgumentError("outcome_map must be :infer, :identity, or a Dict"))
    end
    if a.semantics == :channel
        a.nf.fwd == b.nf.fwd || return _mismatch("final Clifford frames differ", 0)
        tf = substitute(a.nf.pauli_frame, v -> haskey(vmap, v) ? vmap[v] : affine_var(-v))
        tf == b.nf.pauli_frame || return _mismatch("final Pauli frames differ", 0)
    end
    ComparisonResult(true, "", 0, omap)
end

function _bits_to_vars(e, nf::NormalForm)
    e isa Integer && (e = affine_var(e))
    e isa AffineBit || throw(ArgumentError("outcome map values must be Int or AffineBit"))
    all(bit -> haskey(nf.final_var, bit), e.vars) || return nothing
    substitute(e, bit -> affine_var(nf.final_var[bit]))
end
