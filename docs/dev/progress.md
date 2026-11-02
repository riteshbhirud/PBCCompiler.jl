# Progress log: verification and resource estimation

Work is on the local `dev` branch, created from `riteshbhirud/PBCCompiler.jl`
`master` at `ed0efec`. Nothing has been pushed.

## M0: Codebase map

- Wrote `docs/dev/codebase_map.md`. It covers every `CircuitOp` variant with exact
  semantics, the entry points and passes, printed compiled outputs (single T,
  CNOT+T, QFT-2, and stats for the repo fixture), the utilities, and the exact
  gadget shape.
- Every convention was checked numerically with a dense-matrix script, not read off
  docstrings:
  - rotations are `exp(-iθP)`;
  - bit = 1 means eigenvalue −1;
  - the gadget implements `exp(-iπ/8 P)` with `|m⟩ = T|+⟩` in all 4 branches;
  - `(A*B)(p) = A(B(p))` for `CliffordOperator`.
- Found 6 compiler issues (B1–B6). Each was reproduced before being written down,
  and they are recorded in `docs/NOTES_FOR_COMPILER.md`.
- Deviations recorded:
  - the default semantics is `:outcomes`, not the full channel;
  - there is no QFT-4 in the repo and it is not exact in Clifford+T;
  - the Pauli frame is tracked without a sign;
  - `PrepMagic` is unsupported (no defined semantics).
- Tests: the existing suite was run before starting and passed 178/178. M0 adds no
  code.

## V1: Well-formedness checks

- Added the `PBCCompiler.Verification` submodule (`src/verification/Verification.jl`),
  `affine.jl` (`AffineBit`, the F2-affine outcome expression used everywhere), `ops.jl`
  (op accessors written independently of the compiler's helpers), and `wellformed.jl`
  (`check_wellformed`, `wellformedness_issues`, `WellformednessIssue`).
- Rules: qubit range; duplicate qubits; Pauli length matching the qubit list;
  Hermitian measured and rotated Paulis; disjoint control and target qubits;
  read-before-write and dangling bits; bit rewrites; conditional measurements; and
  magic-qubit discipline (prepared by exactly one joint measurement, consumed by one
  single-qubit measurement, no other use).
- Tests: `test/verification/test_wellformed.jl` violates each rule once and checks the
  reported index and rule (37 assertions). Real compiler output is clean, and compiler
  bug B3 (gadget bit collision) is flagged as `:bit_rewrite`.
- The full suite passes, 215/215. The first full run caught an `affine_var` bug
  (Tuple vs Vector) through its doctest; it was fixed before this commit.

## V2: Symbolic Pauli frame

- `src/verification/frame.jl`: `SymbolicPauli`, with X/Z bits that are `AffineBit`s of
  the outcome variables. Operations: `mul_conditional!` (multiply by `Q^c`),
  `mul_frame!`, `conjugate!` (by a `CliffordOperator` or a vector of generator
  images), `commutation_parity`, `evaluate`, `substitute`, and `==` (equality as
  functions, which is exact because the affine representation is canonical).
- The frame is tracked without a sign (deviation 3 in the codebase map: the sign is
  not affine in general, and only commutation parities are observable).
- Tests (`test/verification/test_frame.jl`, 48 assertions):
  - `AffineBit` algebra;
  - 8 random frames over **12 outcome variables**, each put through 6 random
    conditional multiplications and Clifford conjugations, compared bit for bit
    against explicit per-branch replay for **all 2^12 assignments**, including
    commutation parities;
  - equality-as-functions agrees with exhaustive comparison on 20 random pairs.
- Test power check: a planted bug (X/Z swapped in the parity) is caught in 17 of 20
  random trials by the same exhaustive comparison.
- Added `Random` and `LinearAlgebra` to `test/Project.toml`. The `Pkg.test` sandbox
  cannot load undeclared stdlibs.
- Full suite: 263/263 pass.

## V3: Rotation normal form

- `paulimap.jl`: `PauliMap`, a Clifford frame stored as the 2n generator images
  (signs included), with `rotate!` (O(n) Pauli products per Clifford rotation),
  `apply_map`, `compose`, and `to_clifford`. The update rules for `exp(-ikπ/4·Q)`
  are derived in the file header.
- `lowering.jl`: `lower` turns a `Circuit` into `VOp`s. Each measurement write gets
  its own SSA variable, so bit overwrites are sound. `PauliConditional` becomes
  `(P_c)_{π/4}(P_t)_{π/4}(P_cP_t)_{−π/4}` (derived from `C = (I+P_c+P_t−P_cP_t)/2`).
  Unsupported constructs (`PrepMagic`, nested conditionals, conditional
  measurements, reads before writes) are collected rather than guessed at.
- `normal_form.jl`: `NormalForm` / `NFItem`, the "push Cliffords to the end"
  builder. Angles are exact `Rational`s in units of π. Unconditional Clifford
  rotations are absorbed into the frame, with a sign correction folded into the
  Pauli frame. Paulis (including conditional ones) go into the symbolic Pauli frame.
  Everything else becomes an item with an affine sign or condition. The constant part
  of the Pauli frame is folded into the Clifford frame at the end.
- `dense.jl`: a dense simulator that gives each `CircuitOp` its semantics directly
  (`PauliConditional` via the projector definition, not via the lowering), so it is
  an independent cross-check.
- `random_circuits.jl`: the seeded `random_clifford_t_circuit` generator (qubits,
  depth, T density, measurement density, classical-control density, multi-qubit
  rotations, and an option to avoid compiler bugs B1/B2).
- Tests (`test_normal_form.jl`, 724 assertions):
  - `PauliMap` against dense `U P U†` and `U† P U` (40 random sequences, n ≤ 4),
    against QuantumClifford's own `apply!` (exact phases), and `compose`;
  - the `PauliConditional` decomposition against the projector definition (25 cases);
  - **60 random unitary Clifford+T circuits (n ≤ 6)**: the input's dense unitary
    equals `F · items` up to phase, checked generator by generator against both
    `fwd` and `inv`;
  - **60 random circuits with mid-circuit measurements and classical control**: the
    per-outcome POVM of the normal-form items equals the circuit's POVM.
- Test power check: flipping the sign of one measurement item is detected in 50 of
  60 random circuits. The other 10 flip outcomes that are uniformly random for every
  input, where the flip is semantically a no-op.
- `LinearAlgebra` was added to `[deps]`, with compat `1.11.0` mirroring `Random`.
- Full suite: 987/987 pass.

## V4: Gadget recognition and lifting

- `src/gadget_shapes.jl`: the `PBCCompiler.GadgetShapes` module. It does template-driven
  structural recognition (`GADGET_TEMPLATES`, `find_gadgets`, `count_gadgets`) and
  depends only on the data types and QuantumClifford, so the resource estimator (E8)
  can share it without depending on the verifier.
- `cyclotomic.jl`: exact arithmetic in Q(ζ₈) (`Cyc8`), with field inverse via the
  Galois norm.
- `gadgets.jl`:
  - `synthesize_z1_clifford` builds V with V·P·V† = +Z₁ as an explicit rotation
    sequence and checks it exactly;
  - each correction must map to ±Z₁ under V;
  - `prove_reduced_gadget` proves the 2-qubit reduced gadget **exactly in Q(ζ₈)**
    over all four branches: diagonal and unitary per branch, trace preserving, angle
    a multiple of π/8, residual Pauli frame affine in the gadget outcomes;
  - proofs are cached by (template, reduced shape, magic state);
  - `lift_gadgets` and `lower_lifted` splice the proved rotation, plus any residual
    frame, into the lowered program.
- `lower` gained a `replacements` hook, so gadgets are spliced in during lowering.
- Results:
  - **all 111 gadgets of the repo fixture lift**, with angle 1/8, no residual, and
    every branch probability exactly 1/4;
  - all gadgets of 40 random compiled circuits lift;
  - the second time the same shape is seen, the proof is a cache hit.
- Corrupted gadgets each fail with an explicit reason:
  - π/4 correction replaced by π/2: rejected;
  - swapped conditioning bits: rejected;
  - missing π/4 correction: rejected;
  - correction about the wrong Pauli: rejected;
  - magic state |+⟩: rejected;
  - magic state S|+⟩: not a T;
  - magic qubit measured in Z: not recognized.
  - The missing-π/2 and negated-π/4 cases lift with an explicit residual Pauli frame,
    which is mathematically correct. They are caught at comparison under the identity
    outcome map.
- An independent floating-point cross-check of the exact prover over all 144
  (joint letter, magic letter, sign, magic state) combinations with random
  corrections agrees on every success and every failure.
- Correction during the work: I first predicted that the |+⟩ magic state gives
  non-unitary branches. The prover showed (and a hand recomputation confirmed) that
  the X measurement makes them unitary, and the rejection reason is "different
  rotations". The test was fixed to match the mathematics.
- Aqua flagged constructor ambiguities from `Cyc8 <: Number`. Fixed by dropping the
  supertype and adding a typed inner constructor. Aqua is clean again.
- Full suite: 1344/1344 pass.

## V5: Canonicalization and comparison

- `canonical.jl`: `canonicalize` does four things, repeated to a fixed point:
  1. prunes items that cannot affect any outcome (`:outcomes` semantics only);
  2. merges same-Pauli rotations across commuting items, re-absorbing Clifford results
     into the frame through an exact rebuild;
  3. **normalizes angles modulo Cliffords** (each unconditional constant-sign
     non-Clifford rotation becomes angle r ∈ (0, π/4), with the Clifford remainder
     absorbed into the frame);
  4. orders the items by Foata layers with a key that does not depend on variable names.
- `compare.jl`: `compare_forms` walks both canonical lists and builds the outcome map
  (`r_in = r_out ⊕ s_in ⊕ s_out` for matched measurements). Supports `:infer`,
  `:identity`, and explicit maps. For `:channel`, the frames must also match.
- `types.jl`: `VerificationStatus`, `Certificate`, `Counterexample` (fully defined
  now, filled in V7), and `VerificationResult`. `verify.jl`: `verify_compilation`
  (two-argument and compile-it-yourself forms), with a readable `show`. A program
  that violates the IR contract is never certified (status `UNKNOWN`, method
  `:wellformedness`).
- **Found and fixed an incompleteness:**
  - The repo fixture first came out `UNKNOWN`: its normal forms differed at
    canonical position 94.
  - Dense POVMs (2 data qubits) showed the compiled program *is* correct (max diff
    1.3e-14), so the gap was in the verifier.
  - Delta debugging reduced it to `T†·T·T·S_X; M_Z`. The verifier merged T†·T and
    kept a +T; the compiler merged T·T into S and kept the T†. They are equal because
    T = T†·S, but the two normal forms differed.
  - Fix: step 3 above, angle normalization modulo Cliffords. The fixture now verifies
    `EQUIVALENT` with the identity outcome map in about 0.03 s.
- Tests (`test_compare.jl`, 293 assertions):
  - every compiler test circuit verifies `EQUIVALENT` via `:normal_form` with the
    identity map: the existing compiler tests, the fixture, QFT-2 through the parser,
    20 circuits from the compiler's own `random_test_circuit`, and 60 random
    Clifford+T circuits;
  - hand-built equivalent pairs: reordering, merging, splitting, merging across
    commuting ops, Clifford through measurement, parity outcome maps
    (`:infer` / `:identity` / explicit), classical control, `:channel` semantics;
  - **soundness:** 200 random (circuit, mutated circuit) pairs with dense POVM ground
    truth; `EQUIVALENT` is never returned when the POVMs differ (more than 50 pairs
    were inequivalent, more than 10 proved equivalent);
  - B1, B2, B3 and B4 are never certified.
- Full suite: 1637/1637 pass.

## V6: Dense and stabilizer oracles

- `oracle_dense.jl`: `dense_oracle` builds the input's visible POVM and the compiled
  program's raw POVM (all outcomes, gadget ones included) by exact enumeration of
  every branch, with magic qubits prepared in `T|+⟩`.
  - `:identity` and explicit maps are checked directly.
  - **`:infer` searches every affine outcome map exactly** via Walsh–Hadamard
    transforms (`V_T = ±W_{S(T)}` for all T), trying the identity map first and then
    backtracking with a node budget. An exhausted budget is reported as undecided,
    never as a refutation.
  - `:channel` additionally compares per-record Choi matrices via Kraus Gram sums.
  - It uses a tolerance of 1e-9, so it refutes and never certifies.
- `oracle_stabilizer.jl`: `stabilizer_oracle` (refutation only) uses QuantumClifford's
  `GeneralizedStabilizer`. It draws random stabilizer inputs, samples a record from
  the input, and computes the record's exact probability under both programs, with
  hidden branches enumerated.
- `oracle_dispatch.jl`: `:auto` tries the normal form first, then the dense oracle if
  the instance is small enough, otherwise the stabilizer oracle. A numerical
  agreement is reported as `UNKNOWN`. Under `:infer`, a stabilizer difference found
  with the identity map is reported as evidence (`UNKNOWN`), because it does not rule
  out a parity map.
- The one-argument `verify_compilation(input)` now defaults to
  `outcome_map = :identity`, which is PBCCompiler's runtime contract.
- **Upstream QuantumClifford 0.11.8 bug found:** `apply!(::GeneralizedStabilizer,
  ::PauliChannel)` evaluates `im^(phaseₗ − phaseᵣ + 4)`, which throws `DomainError` when
  a channel Pauli has a −1 sign. Worked around by passing unsigned Paulis. Worth
  reporting upstream.
- Performance fixes made during this milestone:
  - Forced projection as a non-unitary Pauli channel grew χ 4× per measurement (a
    test took more than 10 minutes). It now uses QuantumClifford's own tableau-updating
    measurement routine (internal `_projectrand_notnorm`, wrapped in one place).
  - The `:infer` map search hung on one case (8 visible bits, 10 raw bits) because it
    enumerated candidate combinations before trying the identity map. It now tries
    the identity first and backtracks with pruning.
- Tests (`test_oracles.jl`):
  - **300 seeded random cases: the normal form never contradicts the dense oracle**
    (≥ 200 decided, ≥ 50 proved equivalent, ≥ 30 refuted);
  - the dense oracle finds parity maps that the normal form cannot;
  - the channel-semantics refutation works;
  - the stabilizer oracle never refutes a correct compilation (40 programs), and
    catches ≥ 80% of dense-confirmed channel-changing mutants;
  - B1 on 12 qubits (beyond the dense limit) is refuted by the stabilizer oracle.
- **Long run** (`benchmark/verification/long_crosscheck.jl 3000 1`):
  3000 cases, 1434 proved `EQUIVALENT` by the normal form, 853 refuted by the dense
  oracle, **0 contradictions**, 188 s.

## V7: Counterexamples and minimization

- `counterexample.jl`:
  - `minimize_counterexample` delta-debugs the *input* and recompiles each candidate
    through `preprocess_circuit` (the compiler is deterministic; checked in M0);
  - `minimize_pair` delta-debugs the visible bits of a given (input, compiled) pair;
  - every counterexample is re-confirmed by the dense oracle when small enough and
    marked confirmed or unconfirmed.
- A `Counterexample` reports:
  - the first diverging canonical items;
  - a dense input state, or a stabilizer state from the stabilizer oracle;
  - the visible record and its probability under each program;
  - the minimized programs.
- Tests (`test_counterexample.jl`, 121 assertions):
  - B1 buried in random 4-qubit circuits (5 cases) is minimized to a few ops
    containing the reversed CNOT, confirmed, and independently re-checked;
  - B2 is minimized to the T gate plus measurements;
  - 30 mutated compiled outputs give confirmed, pair-minimized counterexamples;
  - beyond the dense limit the witness is a stabilizer state.
- Finding: **B2 only bites when the unsorted T is not moved past an earlier
  Clifford**, because moving it re-embeds the Pauli correctly. The verifier correctly
  certified such a circuit, so my test expectation was wrong. The notes now say
  exactly when B2 occurs.

## V8: Per-pass translation validation

- `passwise.jl`: `PREPROCESS_PASSES` mirrors `preprocess_circuit`'s pass order (hook
  request H1 in `docs/NOTES_FOR_COMPILER.md`), plus `run_passes`, `verify_passes` and
  `verify_passes_detailed`. Each pass's output is verified against its input, with the
  original input's bits as the visible record.
- To support intermediate programs, `verify_compilation` gained `visible_bits`. It also
  lifts gadgets on the input side, and both oracles now prepare magic qubits in either
  program.
- Results:
  - **B1 is localized to `remove_pauliconditional`, B2 to `remove_nonclifford`**;
  - the fixture and the compiler's own test circuits are `EQUIVALENT` at every pass.
- Drift guard: running `PREPROCESS_PASSES` in sequence equals `preprocess_circuit` on
  30 random programs, and the early return for programs with fewer than two ops is
  mirrored.
- Full suite at this point: 2273/2273 pass. That includes the V9 test files written
  so far; those are committed with V9.
