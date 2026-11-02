# Affine functions over F2 of measurement-outcome variables.
#
# A value `c ⊕ x_{v1} ⊕ x_{v2} ⊕ …` is stored as the constant `c` plus the sorted,
# duplicate-free list of variable ids. Variable ids are opaque integers; the normal
# form uses one id per measurement *write* (so a classical bit that is overwritten
# gets a fresh id for each write).

"""
    AffineBit(constant::Bool, vars)

An affine function over F2 of outcome variables: `constant ⊕ (⊕ᵢ x_{vars[i]})`.

`vars` is normalized to a sorted list with every variable appearing at most once
(a variable listed twice cancels). Supports `⊻` (addition over F2), evaluation on
an assignment, and structural equality, which coincides with equality as functions.

```jldoctest
julia> using PBCCompiler.Verification: AffineBit, affine_var

julia> a = affine_var(3) ⊻ affine_var(1) ⊻ AffineBit(true)
1 ⊕ m1 ⊕ m3

julia> a ⊻ affine_var(3)
1 ⊕ m1
```
"""
struct AffineBit
    constant::Bool
    vars::Vector{Int}
    function AffineBit(constant::Bool, vars::AbstractVector{<:Integer})
        v = sort!(collect(Int, vars))
        # drop pairs (x ⊕ x = 0)
        out = Int[]
        i = 1
        while i <= length(v)
            j = i
            while j < length(v) && v[j + 1] == v[i]
                j += 1
            end
            isodd(j - i + 1) && push!(out, v[i])
            i = j + 1
        end
        new(constant, out)
    end
end

AffineBit(constant::Bool) = AffineBit(constant, Int[])

"""
    affine_var(i::Integer) -> AffineBit

The affine function consisting of the single outcome variable `i`.
"""
affine_var(i::Integer) = AffineBit(false, [Int(i)])

Base.:(==)(a::AffineBit, b::AffineBit) = a.constant == b.constant && a.vars == b.vars
Base.hash(a::AffineBit, h::UInt) = hash(a.vars, hash(a.constant, hash(:AffineBit, h)))
Base.iszero(a::AffineBit) = !a.constant && isempty(a.vars)
Base.isone(a::AffineBit) = a.constant && isempty(a.vars)

"""
    isconstant(a::AffineBit) -> Bool

Whether `a` does not depend on any outcome variable.
"""
isconstant(a::AffineBit) = isempty(a.vars)

# Symmetric difference of two sorted, duplicate-free vectors.
function _symdiff_sorted(a::Vector{Int}, b::Vector{Int})
    out = Vector{Int}(undef, 0)
    sizehint!(out, length(a) + length(b))
    i = 1; j = 1
    @inbounds while i <= length(a) && j <= length(b)
        if a[i] < b[j]
            push!(out, a[i]); i += 1
        elseif b[j] < a[i]
            push!(out, b[j]); j += 1
        else
            i += 1; j += 1
        end
    end
    @inbounds while i <= length(a)
        push!(out, a[i]); i += 1
    end
    @inbounds while j <= length(b)
        push!(out, b[j]); j += 1
    end
    out
end

# Internal constructor that trusts its input to be sorted and duplicate-free.
_affine_unchecked(c::Bool, v::Vector{Int}) = (r = AffineBit(c); append!(r.vars, v); r)

Base.xor(a::AffineBit, b::AffineBit) = _affine_unchecked(a.constant ⊻ b.constant, _symdiff_sorted(a.vars, b.vars))
Base.xor(a::AffineBit, b::Bool) = _affine_unchecked(a.constant ⊻ b, copy(a.vars))
Base.xor(b::Bool, a::AffineBit) = xor(a, b)

"""
    affine_and_const(a::AffineBit, b::Bool) -> AffineBit

`a` if `b` is true, the zero function otherwise (multiplication by a constant bit).
"""
affine_and_const(a::AffineBit, b::Bool) = b ? a : AffineBit(false)

"""
    evaluate(a::AffineBit, assignment) -> Bool

Value of `a` when variable `v` takes the value `assignment[v]` (any indexable
collection or `AbstractDict` of `Bool`s).
"""
function evaluate(a::AffineBit, assignment)
    r = a.constant
    for v in a.vars
        r ⊻= Bool(assignment[v])
    end
    r
end

"""
    substitute(a::AffineBit, f) -> AffineBit

Replace every variable `v` of `a` by the affine function `f(v)::AffineBit`
(variables for which `f` returns `affine_var(v)` are kept).
"""
function substitute(a::AffineBit, f)
    r = AffineBit(a.constant)
    for v in a.vars
        r = r ⊻ (f(v)::AffineBit)
    end
    r
end

function Base.show(io::IO, a::AffineBit)
    if isempty(a.vars)
        print(io, a.constant ? "1" : "0")
        return
    end
    parts = String[]
    a.constant && push!(parts, "1")
    for v in a.vars
        push!(parts, "m$v")
    end
    print(io, join(parts, " ⊕ "))
end
