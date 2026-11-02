"""
    random_clifford_t_circuit(rng, n; depth = 10n, t_density = 0.25, meas_density = 0.0,
                              cc_density = 0.0, multiqubit_density = 0.2,
                              final_measure = true, sorted_qubits = false) -> Circuit

A random circuit on qubits `1:n` built from the IR's own op types:

- single-qubit Clifford rotations (`ExpQuatPiPauli` / `ExpHalfPiPauli` on X, Y or Z,
  either sign), Hadamards (as `Z,X,Z` quarter turns, like the parser), and CNOT-like
  `PauliConditional(Z, [c], X, [t])`;
- T-type gates (`ExpEighPiPauli` on ±Z), with probability `t_density` per gate;
- with probability `multiqubit_density`, a rotation (π/4 or π/8) about a random
  multi-qubit Pauli;
- with probability `meas_density`, a mid-circuit Pauli measurement into a fresh bit;
- with probability `cc_density` (once some bit exists), a `BitConditional`
  single-qubit Clifford conditioned on a random earlier bit;
- if `final_measure`, a Z measurement of every qubit at the end, into fresh bits.

`sorted_qubits = true` keeps every qubit list sorted and every CNOT control below
its target, which avoids the known compiler bugs B1/B2 (see `docs/NOTES_FOR_COMPILER.md`).
"""
function random_clifford_t_circuit(rng::AbstractRNG, n::Int; depth::Int = 10n, t_density::Real = 0.25,
                                   meas_density::Real = 0.0, cc_density::Real = 0.0,
                                   multiqubit_density::Real = 0.2, final_measure::Bool = true,
                                   sorted_qubits::Bool = false)
    n >= 1 || throw(ArgumentError("n must be ≥ 1"))
    c = CircuitOp.Type[]
    nbits = 0
    letter(rng) = rand(rng, ('X', 'Y', 'Z'))
    single(l, s) = PauliOperator(s ? 0x02 : 0x00, [l == 'X' || l == 'Y'], [l == 'Z' || l == 'Y'])
    for _ in 1:depth
        u = rand(rng)
        q = rand(rng, 1:n)
        if nbits > 0 && rand(rng) < cc_density
            inner = rand(rng, Bool) ? CircuitOp.ExpHalfPiPauli(single(letter(rng), false), [q]) :
                                      CircuitOp.ExpQuatPiPauli(single(letter(rng), rand(rng, Bool)), [q])
            push!(c, CircuitOp.BitConditional(inner, rand(rng, 1:nbits)))
        elseif rand(rng) < meas_density
            nbits += 1
            qs = _random_qubits(rng, n, rand(rng, 1:min(n, 3)), sorted_qubits)
            push!(c, CircuitOp.Measurement(_random_hermitian(rng, length(qs)), nbits, qs))
        elseif n >= 2 && rand(rng) < multiqubit_density
            qs = _random_qubits(rng, n, rand(rng, 2:min(n, 4)), sorted_qubits)
            p = _random_hermitian(rng, length(qs))
            push!(c, rand(rng, Bool) ? CircuitOp.ExpEighPiPauli(p, qs) : CircuitOp.ExpQuatPiPauli(p, qs))
        elseif rand(rng) < t_density
            push!(c, CircuitOp.ExpEighPiPauli(single('Z', rand(rng, Bool)), [q]))
        elseif n >= 2 && u < 0.3
            a, b = _random_qubits(rng, n, 2, sorted_qubits)
            push!(c, CircuitOp.PauliConditional(single('Z', false), [a], single('X', false), [b]))
        elseif u < 0.5
            push!(c, CircuitOp.ExpQuatPiPauli(single('Z', false), [q]))
            push!(c, CircuitOp.ExpQuatPiPauli(single('X', false), [q]))
            push!(c, CircuitOp.ExpQuatPiPauli(single('Z', false), [q]))
        elseif u < 0.85
            push!(c, CircuitOp.ExpQuatPiPauli(single(letter(rng), rand(rng, Bool)), [q]))
        else
            push!(c, CircuitOp.ExpHalfPiPauli(single(letter(rng), false), [q]))
        end
    end
    if final_measure
        for q in 1:n
            nbits += 1
            push!(c, CircuitOp.Measurement(single('Z', false), nbits, [q]))
        end
    end
    c
end

function _random_qubits(rng::AbstractRNG, n::Int, k::Int, sorted::Bool)
    qs = randperm(rng, n)[1:k]
    sorted ? sort!(qs) : qs
end

function _random_hermitian(rng::AbstractRNG, k::Int)
    while true
        x = rand(rng, Bool, k); z = rand(rng, Bool, k)
        (any(x) || any(z)) && return PauliOperator(rand(rng, Bool) ? 0x02 : 0x00, x, z)
    end
end
