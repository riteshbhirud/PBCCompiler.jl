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
