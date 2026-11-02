# [Verifying compilation](@id verification)

```@meta
CurrentModule = PBCCompiler.Verification
```

`PBCCompiler.Verification` checks that a compiled Pauli-based-computation program
implements the same logical operation as the input circuit, and produces a concrete
counterexample when it does not.

```julia
using PBCCompiler, QuantumClifford
using PBCCompiler: CircuitOp

input = [CircuitOp.PauliConditional(P"Z", [1], P"X", [2]), CircuitOp.ExpEighPiPauli(P"Z", [2]),
         CircuitOp.Measurement(P"Z", 1, [1]), CircuitOp.Measurement(P"Z", 2, [2])]

r = verify_compilation(input)           # compiles with preprocess_circuit, then verifies
r.status                                # EQUIVALENT
r.method                                # :normal_form

compiled = PBCCompiler.preprocess_circuit(copy(input))
verify_compilation(input, compiled)     # verify a given pair
verify_passes(input)                    # one result per compiler pass
```

The verifier does not reuse any of the compiler's transformation code. It relies
only on the `CircuitOp` data types and QuantumClifford.jl primitives, so a bug in a
compiler pass cannot hide itself.

## What is checked

**Programs.** The input is a Clifford+T circuit with measurements and possibly
classical control. The compiled program consists of Pauli product measurements,
magic-state gadgets, Clifford rotations and outcome-conditioned corrections. The
conventions (verified numerically, see `docs/dev/codebase_map.md`) are:

- rotations are `exp(-iθP)`: `ExpQuatPiPauli` is θ = π/4, `ExpEighPiPauli` is θ = π/8;
- a measurement writes 1 to its bit when it observes the eigenvalue −1;
- `BitConditional` applies its op when its bit is 1;
- a `PauliConditional` applies the target Pauli when the control Pauli is −1;
- qubits beyond the input's width are magic qubits, implicitly prepared in
  `T|+⟩ = (|0⟩ + e^{iπ/4}|1⟩)/√2`.

**Semantics.** Two semantics are available (`semantics` keyword):

- `:outcomes` (default): for every input state, the joint distribution of the
  input's classical outcomes (the final value of every bit it writes) is the same.
  This is what PBCCompiler preserves: it drops everything after the last
  measurement and conjugates measurements, so the post-measurement quantum state is
  *not* preserved.
- `:channel`: additionally, for every outcome record, the post-measurement quantum
  operation on the data qubits is the same (magic qubits are traced out).

Both are taken up to a global phase per outcome branch.

**Outcome map.** The compiled program may report an input outcome as a parity of
its own raw outcomes, for example when it tracks a Pauli frame. `outcome_map`
controls this:

- `:identity`: input bit `b` is the compiled program's final value of bit `b`. This
  is PBCCompiler's runtime contract, and the default of `verify_compilation(input)`.
- `:infer`: any affine map `bᵢ = cᵢ ⊕ ⊕ⱼ rⱼ` over compiled outcomes is allowed. The
  map found is reported in the certificate. This is the default of the two-argument
  form.
- an explicit `Dict` from input bit to [`AffineBit`](@ref) over compiled bits.

Outcomes that the compiled program uses only internally (gadget outcomes) are
hidden. They must not change the result.

## Method layers

| layer | can return | sound? | complete? |
|---|---|---|---|
| well-formedness ([`check_wellformed`](@ref)) | blocks `EQUIVALENT` | — | — |
| gadget lifting ([`lift_gadgets`](@ref)) | (input to the normal form) | exact proof | for the recognized templates |
| normal form ([`normal_form`](@ref), [`canonicalize`](@ref), [`compare_forms`](@ref)) | `EQUIVALENT` | **yes** | for the rewrites listed below |
| dense oracle ([`dense_oracle`](@ref)) | `NOT_EQUIVALENT` | yes (tolerance 1e-9) | yes, up to its size limit |
| stabilizer oracle ([`stabilizer_oracle`](@ref)) | `NOT_EQUIVALENT` | yes (tolerance 1e-6) | no (randomized) |

1. **Well-formedness.** Static contract checks: qubit ranges, Hermitian Paulis,
   read-before-write, bit overwrites, and magic-qubit discipline. A program that
   violates the contract is never certified, even if it happens to be equivalent.

2. **Gadget lifting.** Each magic-state gadget is recognized structurally from a
   template table (`GadgetShapes.GADGET_TEMPLATES`) and then *proved*:
   - a Clifford `V` with `V P V† = Z₁` is synthesized and checked exactly;
   - the gadget is reduced to one data qubit plus the magic qubit;
   - the reduced gadget is checked over all four outcome branches in exact
     arithmetic in the cyclotomic field `Q(ζ₈)`. Each branch must be a common
     rotation `exp(-iθP)` up to a Pauli-frame residual that is affine in the
     gadget's outcomes, and the branches must form a trace-preserving map.
   Proofs are cached by shape.

3. **Normal form.** Both programs are rewritten as
   `PF · F · Iₖ ⋯ I₁`: non-Clifford and conditional rotations and measurements
   (items) in input coordinates, a final Clifford frame `F`, and an outcome-dependent
   Pauli frame `PF` whose bits are affine functions of the outcomes. Angles are exact
   rationals. The canonical form then, until nothing changes:
   - drops items that cannot affect any outcome (outcome semantics only);
   - merges same-Pauli rotations across commuting items;
   - normalizes each non-Clifford angle into (0, π/4), absorbing the Clifford
     remainder;
   - orders items by Foata layers of their dependency DAG.
   Equal canonical forms (and equal frames, for `:channel`) prove equivalence.
   - *Sound:* every rewrite is an exact identity, so equality implies equivalence
     under the reported outcome map.
   - *Complete* for compilers that only push Cliffords through, reorder commuting
     operations, merge or split rotations about the same Pauli, and use recognized
     gadgets. This is exactly what PBCCompiler's passes do.
   - *Cannot decide* equivalences that need other identities: rewriting a
     measurement by stabilizers of the current state (for example measuring `Z₁Z₂`
     after `Z₁` in place of `Z₂`), outcome-dependent rotations that only cancel in
     combination, or programs whose canonical orders differ because items with equal
     keys appear in different orders. These come out `UNKNOWN` (or go to the
     oracles), never `EQUIVALENT`.

4. **Oracles** (`oracle = :auto | :dense | :stabilizer | :none`). With `:auto`, the
   dense oracle runs when all qubits fit (`max_dense_qubits`, default 10), otherwise
   the stabilizer oracle. Oracles only *refute*. A numerical agreement is reported as
   `UNKNOWN` with the deviation in `stats`.
   - The dense oracle enumerates every outcome branch, builds both POVMs, and
     compares them under the outcome map. With `:infer` it searches *all* affine maps
     exactly, through Walsh–Hadamard transforms of the POVMs. A refutation then means
     that no affine map works. If the search would exceed its limits, the result is
     undecided rather than a refutation.
   - The stabilizer oracle uses QuantumClifford's `GeneralizedStabilizer`. It draws
     random stabilizer inputs, samples an outcome record from the input program, and
     computes that record's probability exactly under both programs (hidden branches
     enumerated). Under `:infer`, a difference it finds with the identity map is only
     reported as evidence (`UNKNOWN`), because it does not exclude a parity map.

## Reading a result

```julia
r = verify_compilation(input)
r.status                    # EQUIVALENT | NOT_EQUIVALENT | UNKNOWN
r.method                    # what decided it
r.stats[:reasons]           # why it is not EQUIVALENT, when it is not
r.unverified_components     # compiled op indices that could not be proved
r.wellformedness_errors
r.certificate               # EQUIVALENT only
r.counterexample            # NOT_EQUIVALENT only
```

The certificate holds both canonical forms, the outcome map (input bit ⇒
`AffineBit` over compiled outcome variables, where a variable is the index of a
compiled measurement op and `variable_bits` gives its bit), whether that map is the
identity, and the lifted gadgets with their proofs.

### Counterexamples

A [`Counterexample`](@ref) contains:

- `divergence`: the first canonical position where the normal forms differ, with the
  input and compiled items there. This is where to start reading.
- `input_state` (a dense vector on the data qubits, from the dense oracle) or
  `input_stabilizer` (from the stabilizer oracle): an input on which the programs
  behave differently.
- `branch`: an outcome record (input bit ⇒ value), with its probability under the
  input (`prob_input`) and the compiled program (`prob_compiled`). These differ.
- `minimized_input`, `minimized_compiled`: a small failing pair. For
  `verify_compilation(input)` this is a 1-minimal input found by delta debugging,
  recompiling every candidate through the public compiler API
  ([`minimize_counterexample`](@ref)). For a given pair, the visible bits are
  minimized ([`minimize_pair`](@ref)).
- `confirmed`: whether the dense oracle re-confirmed the difference. Large
  stabilizer-oracle witnesses that cannot be minimized below the dense limit stay
  unconfirmed.

For example, a CNOT whose control index is larger than its target index exposes a
compiler bug (B1 in `docs/NOTES_FOR_COMPILER.md`):

```julia
r = verify_compilation([CircuitOp.PauliConditional(P"Z", [2], P"X", [1]),
                        CircuitOp.Measurement(P"Z", 1, [1]), CircuitOp.Measurement(P"Z", 2, [2])])
r.status                           # NOT_EQUIVALENT
r.counterexample.branch            # Dict(1 => true, 2 => false)
r.counterexample.prob_input        # ≈ 0.93
r.counterexample.prob_compiled     # ≈ 0.07
```

and [`verify_passes`](@ref) localizes it to `remove_pauliconditional`.

## Limitations

- `EQUIVALENT` is only ever produced by the normal form (exact). The dense oracle
  uses floating point, so its agreement is never upgraded to a proof.
- The normal form is incomplete outside the rewrite system described above.
- Only gadget shapes in the template table are lifted. A gadget that does not match
  the template, or fails its proof, leaves its ops unverified (they are reported in
  `unverified_components`), and the result then depends on the oracles.
- The magic state is an assumption (`magic_k`, default `T|+⟩`), since the IR has no
  preparation op. `PrepMagic` has no defined semantics and is unsupported, as are
  nested `BitConditional`s and conditional measurements.
- The dense oracle is limited to `max_dense_qubits` total qubits (data plus magic)
  and `max_branches` outcome branches. The `:infer` map search has its own size
  limits.
- The stabilizer oracle cost grows exponentially with the number of non-Clifford
  gates. It needs single-bit outcome maps (no parity maps), and it can miss
  differences.
- Per-pass validation mirrors `preprocess_circuit`'s pass order in
  [`PREPROCESS_PASSES`](@ref). A test detects drift, but the compiler does not yet
  expose the order as data.

## API

```@autodocs
Modules = [PBCCompiler.Verification, PBCCompiler.GadgetShapes]
Private = false
```
