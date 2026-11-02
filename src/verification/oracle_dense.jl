# V6: exact dense oracle for small instances (floating point, tolerance-based).
#
# The input program's visible POVM {E_in[b]} (b = final values of the bits it
# writes) is compared with the compiled program's raw POVM {E_out[r]} (r = the
# outcomes of *all* its measurements, including hidden gadget outcomes), through an
# outcome map M: input bit j = cⱼ ⊕ ⟨Sⱼ, r⟩. Equivalence under M means
#     E_in[b] = Σ_{r : M(r) = b} E_out[r]   for every b,
# or equivalently, in the Walsh basis (V_T = Σ_b (−1)^{⟨T,b⟩} E_in[b],
# W_S = Σ_r (−1)^{⟨S,r⟩} E_out[r]):
#     V_T = (−1)^{⟨c,T⟩} W_{⊕_{j∈T} Sⱼ}   for every subset T of input bits.
# With `outcome_map = :infer`, candidate (Sⱼ, cⱼ) are found from the single-bit
# conditions, then checked on every T. A refutation therefore means that *no*
# affine outcome map works. Being floating point, this oracle can refute but never
# certify.

"""
    DenseVerdict

Result of [`dense_oracle`](@ref):

- `decided`: every branch was enumerated within the limits;
- `agree`: some admissible outcome map reproduces the input's outcome
  distribution, and the post-measurement channel too for `:channel`;
- `outcome_map`: that map (input bit ⇒ affine function of compiled measurement-op
  variables);
- `max_deviation`: largest entrywise deviation under the best map tried;
- branch counts;
- a witness when they disagree: the input record, a dense input state, and the
  record's probability under each program;
- `note`.
"""
struct DenseVerdict
    decided::Bool
    agree::Bool
    outcome_map::Union{Nothing,Dict{Int,AffineBit}}
    max_deviation::Float64
    leaves_input::Int
    leaves_compiled::Int
    witness_record::Dict{Int,Bool}
    witness_state::Vector{ComplexF64}
    p_input::Float64
    p_compiled::Float64
    note::String
end

_undecided(note) = DenseVerdict(false, false, nothing, NaN, 0, 0, Dict{Int,Bool}(), ComplexF64[], NaN, NaN, note)

# op index of the last top-level measurement writing each bit
function _final_writes(c)
    fw = Dict{Int,Int}()
    for (i, op) in enumerate(c)
        _opkind(op) == OP_MEASUREMENT && (fw[_measurement_bit(op)] = i)
    end
    fw
end

_measurement_ops(c) = [i for (i, op) in enumerate(c) if _opkind(op) == OP_MEASUREMENT]

"""
    dense_oracle(input, compiled; n_data, semantics = :outcomes, outcome_map = :identity,
                 magic_k = 1, max_qubits = 10, max_branches = 4096, tol = 1e-9) -> DenseVerdict

Compare two programs exactly (up to floating point) by dense simulation of every
outcome branch. Compiled qubits above `n_data` are magic qubits, prepared in
`(|0⟩ + ζ₈^magic_k|1⟩)/√2` (in either program). `outcome_map` is `:identity`,
`:infer` (search all affine maps), or a `Dict` input bit ⇒ `AffineBit` over compiled
bits. `visible_bits` selects the input bits forming the record (default: all).
"""
function dense_oracle(input::AbstractVector{CircuitOp.Type}, compiled::AbstractVector{CircuitOp.Type};
                      n_data::Int = _width(input), semantics::Symbol = :outcomes, outcome_map = :identity,
                      magic_k::Int = 1, max_qubits::Int = 10, max_branches::Int = 4096, tol::Float64 = 1e-9,
                      visible_bits::Union{Nothing,AbstractVector{Int}} = nothing)
    n = n_data
    N = max(n, _width(compiled))
    Nin = max(n, _width(input))
    max(N, Nin) <= max_qubits || return _undecided("$(max(N, Nin)) qubits exceed the dense limit $max_qubits")
    D = 1 << n
    # ---- input visible POVM (and Kraus lists for the channel semantics)
    fw_in = _final_writes(input)
    bits = visible_bits === nothing ? sort!(collect(keys(fw_in))) : sort!(collect(Int, visible_bits))
    all(b -> haskey(fw_in, b), bits) || return _undecided("a visible bit is never written by the input")
    v = length(bits)
    v <= 16 || return _undecided("too many visible bits ($v)")
    Ein = [zeros(ComplexF64, D, D) for _ in 1:(1 << v)]
    Kin = [Matrix{ComplexF64}[] for _ in 1:(1 << v)]
    nin, done_in = try
        _dense_dfs(input, Nin, _initial_operator(n, Nin, magic_k); max_leaves = max_branches) do K, reg, outs
            b = sum((reg[bits[j]] ? 1 : 0) << (j - 1) for j in 1:v; init = 0)
            Ein[b + 1] .+= K' * K
            if semantics == :channel
                for m in 0:(1 << (Nin - n))-1           # trace out magic qubits
                    push!(Kin[b + 1], K[[(m << n) + j + 1 for j in 0:D-1], :])
                end
            end
        end
    catch err
        err isa ArgumentError || rethrow()
        return _undecided("input not simulable: " * err.msg)
    end
    done_in || return _undecided("input has more than $max_branches outcome branches")
    # ---- compiled raw POVM
    mops = _measurement_ops(compiled)
    R = length(mops)
    R <= 20 || return _undecided("compiled program has too many measurements ($R) for the dense oracle")
    pos = Dict(op => k for (k, op) in enumerate(mops))
    Eout = Dict{Int,Matrix{ComplexF64}}()
    Kout = Dict{Int,Vector{Matrix{ComplexF64}}}()
    K0 = _initial_operator(n, N, magic_k)
    nout, done_out = try
        _dense_dfs(compiled, N, K0; max_leaves = max_branches) do K, reg, outs
            r = 0
            for (i, _, val) in outs
                val && (r |= 1 << (pos[i] - 1))
            end
            Eout[r] = get(Eout, r, zeros(ComplexF64, D, D)) + K' * K
            if semantics == :channel
                lst = get!(Kout, r, Matrix{ComplexF64}[])
                for m in 0:(1 << (N - n))-1         # trace out the magic qubits
                    rows = [(m << n) + j + 1 for j in 0:D-1]
                    push!(lst, K[rows, :])
                end
            end
        end
    catch err
        err isa ArgumentError || rethrow()
        return _undecided("compiled program not simulable: " * err.msg)
    end
    done_out || return _undecided("compiled program has more than $max_branches outcome branches")

    fw_out = _final_writes(compiled)
    # candidate maps: (S, c) per input bit, with S a bitmask over raw positions
    function fixed_map(spec)
        Ss = zeros(Int, v); cs = falses(v)
        for (j, bit) in enumerate(bits)
            e = if spec === :identity
                haskey(fw_out, bit) || return nothing
                affine_var(bit)
            elseif !(spec isa AbstractDict)
                throw(ArgumentError("outcome_map must be :identity, :infer or a Dict (got $spec)"))
            else
                haskey(spec, bit) || return nothing
                want = spec[bit]
                want isa Integer ? affine_var(want) : want
            end
            cs[j] = e.constant
            for b in e.vars
                haskey(fw_out, b) || return nothing
                Ss[j] |= 1 << (pos[fw_out[b]] - 1)
            end
        end
        (Ss, cs)
    end
    mapped_povm(Ss, cs) = begin
        G = [zeros(ComplexF64, D, D) for _ in 1:(1 << v)]
        for (r, E) in Eout
            b = 0
            for j in 1:v
                (cs[j] ⊻ isodd(count_ones(r & Ss[j]))) && (b |= 1 << (j - 1))
            end
            G[b + 1] .+= E
        end
        G
    end
    deviation(G) = maximum(maximum(abs, Ein[b] - G[b]) for b in eachindex(Ein))

    candidates = Tuple{Vector{Int},BitVector}[]
    complete = true
    if outcome_map === :infer
        # try the identity map first; search the affine maps only if it fails
        fm = fixed_map(:identity)
        if fm !== nothing
            G = mapped_povm(fm...)
            if deviation(G) <= tol && (semantics == :outcomes || _channels_match(Kin, Kout, fm[1], fm[2], v, tol))
                push!(candidates, fm)
            end
        end
        if isempty(candidates)
            complete = _infer_maps!(candidates, Ein, Eout, v, R, D, tol)
            fm === nothing || push!(candidates, fm)        # used for the witness if nothing matches
        end
    else
        fm = fixed_map(outcome_map)
        fm === nothing && return _dense_verdict_nomap(nin, nout, bits, outcome_map)
        push!(candidates, fm)
    end
    best = nothing; bestdev = Inf
    for (Ss, cs) in candidates
        G = mapped_povm(Ss, cs)
        d = deviation(G)
        if d <= tol && (semantics == :outcomes || _channels_match(Kin, Kout, Ss, cs, v, tol))
            omap = Dict(bits[j] => AffineBit(cs[j], [mops[k] for k in 1:R if (Ss[j] >> (k - 1)) & 1 == 1]) for j in 1:v)
            return DenseVerdict(true, true, omap, d, nin, nout, Dict{Int,Bool}(), ComplexF64[], NaN, NaN,
                                outcome_map === :infer ? "an affine outcome map reproduces the input" : "")
        end
        d < bestdev && (bestdev = d; best = (Ss, cs, G))
    end
    complete || return _undecided("the identity map does not match and the affine outcome-map search exceeded its limits")
    # disagreement: witness under the best (or identity) map
    if best === nothing
        fm = fixed_map(:identity)
        fm === nothing && return _dense_verdict_nomap(nin, nout, bits, outcome_map)
        best = (fm[1], fm[2], mapped_povm(fm...))
        bestdev = deviation(best[3])
    end
    G = best[3]
    bw = argmax(b -> opnorm(Ein[b] - G[b]), eachindex(Ein))
    ev = LinearAlgebra.eigen(LinearAlgebra.Hermitian((Ein[bw] - G[bw] + (Ein[bw] - G[bw])') / 2))
    k = argmax(abs.(ev.values))
    ψ = ev.vectors[:, k]
    pin = real(ψ' * Ein[bw] * ψ); pout = real(ψ' * G[bw] * ψ)
    rec = Dict(bits[j] => isodd((bw - 1) >> (j - 1)) for j in 1:v)
    note = if bestdev <= tol
        "outcome distributions agree but the post-measurement states differ"
    elseif outcome_map === :infer
        "no affine outcome map makes the outcome distributions agree"
    else
        "outcome distributions differ under the given outcome map"
    end
    DenseVerdict(true, false, nothing, bestdev, nin, nout, rec, ψ, pin, pout, note)
end

opnorm(A) = LinearAlgebra.opnorm(A)

function _dense_verdict_nomap(nin, nout, bits, outcome_map)
    DenseVerdict(true, false, nothing, Inf, nin, nout, Dict{Int,Bool}(), ComplexF64[], NaN, NaN,
                 "the compiled program does not write every bit the outcome map ($(outcome_map)) refers to")
end

# Search all affine outcome maps through the Walsh transforms (see the file header).
# Appends every map that satisfies all Walsh conditions to `candidates`. Returns
# whether the search was exhaustive (false if a size limit cut it short).
function _infer_maps!(candidates, Ein, Eout, v, R, D, tol; max_combos::Int = 4096)
    (R <= 14 && (1 << R) * D * D <= (1 << 23)) || return false
    W = zeros(ComplexF64, D, D, 1 << R)
    for (r, E) in Eout
        W[:, :, r + 1] .= E
    end
    # fast Walsh–Hadamard transform along the third dimension
    h = 1
    while h < (1 << R)
        for i in 0:2h:(1 << R)-1, j in i:i+h-1
            a = W[:, :, j + 1]; b = W[:, :, j + h + 1]
            W[:, :, j + 1] .= a .+ b
            W[:, :, j + h + 1] .= a .- b
        end
        h *= 2
    end
    V(T) = sum(((-1)^count_ones(b & T)) * Ein[b + 1] for b in 0:(1 << v)-1)
    perbit = Vector{Vector{Tuple{Int,Bool}}}(undef, v)
    for j in 1:v
        Vj = V(1 << (j - 1))
        perbit[j] = Tuple{Int,Bool}[]
        for S in 0:(1 << R)-1
            WS = view(W, :, :, S + 1)
            maximum(abs, Vj - WS) <= tol && push!(perbit[j], (S, false))
            maximum(abs, Vj + WS) <= tol && push!(perbit[j], (S, true))
        end
        isempty(perbit[j]) && return true          # no map can exist
    end
    VT = [V(T) for T in 0:(1 << v)-1]
    # backtracking over bits: after fixing bit j, check every T ⊆ {1..j} that contains j
    chosen = Vector{Tuple{Int,Bool}}(undef, v)
    nodes = Ref(0)
    exhausted = Ref(false)
    function consistent(j)
        for T0 in 0:(1 << (j - 1))-1
            T = T0 | (1 << (j - 1))
            S = 0; c = false
            for k in 1:j
                if (T >> (k - 1)) & 1 == 1
                    S ⊻= chosen[k][1]; c ⊻= chosen[k][2]
                end
            end
            maximum(abs, VT[T + 1] - (c ? -1 : 1) * view(W, :, :, S + 1)) > tol && return false
        end
        true
    end
    function search(j)
        exhausted[] && return
        if j > v
            push!(candidates, ([x[1] for x in chosen], BitVector([x[2] for x in chosen])))
            length(candidates) >= max_combos && (exhausted[] = true)
            return
        end
        for cand in perbit[j]
            nodes[] += 1
            nodes[] > 200_000 && (exhausted[] = true; return)
            chosen[j] = cand
            consistent(j) && search(j + 1)
            exhausted[] && return
        end
    end
    search(1)
    !exhausted[] || !isempty(candidates)
end

# Channel semantics: per visible record, compare the Choi matrices of the two Kraus
# sets via ‖J_in − J_out‖²_F = Σ|⟨aᵢ,aⱼ⟩|² + Σ|⟨bᵢ,bⱼ⟩|² − 2Σ|⟨aᵢ,bⱼ⟩|².
function _channels_match(Kin, Kout, Ss, cs, v, tol)
    groups = [Matrix{ComplexF64}[] for _ in 1:(1 << v)]
    for (r, lst) in Kout
        b = 0
        for j in 1:v
            (cs[j] ⊻ isodd(count_ones(r & Ss[j]))) && (b |= 1 << (j - 1))
        end
        append!(groups[b + 1], lst)
    end
    for b in eachindex(groups)
        A = [vec(K) for K in Kin[b]]; B = [vec(K) for K in groups[b]]
        s = sum((abs2(LinearAlgebra.dot(x, y)) for x in A for y in A); init = 0.0) +
            sum((abs2(LinearAlgebra.dot(x, y)) for x in B for y in B); init = 0.0) -
            2 * sum((abs2(LinearAlgebra.dot(x, y)) for x in A for y in B); init = 0.0)
        s > tol * max(1.0, length(A) + length(B)) && return false
    end
    true
end
