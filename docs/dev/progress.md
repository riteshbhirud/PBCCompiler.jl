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
