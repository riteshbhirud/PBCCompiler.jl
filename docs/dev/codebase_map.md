# PBCCompiler.jl codebase map (M0)

This is the reference that the verifier (`src/verification/`) and the resource
estimator (`src/resource_estimation/`) are built against. Every convention below
was checked numerically against dense matrices (script output reproduced where
relevant), not just read off docstrings. When this file and a docstring in
`src/` disagree, this file describes what the code actually does.

Snapshot: upstream `master` at `ed0efec` (identical to `riteshbhirud/PBCCompiler.jl`
`master` when the `dev` branch was created).

## 1. `CircuitOp` variants (`src/type.jl`)

`CircuitOp` is a Moshi `@data` type. Pauli fields have type
`P = PauliOperator{Array{UInt8,0},Vector{UInt64}}`, written in QuantumClifford's
convention: Hermitian letters `I,X,Y,Z` times a phase `i^p` (`P"Y"` has phase 0).
`pauli[j]` acts on qubit `qubits[j]`, so the qubit order in `qubits` matters and
need not be sorted (`embed(n, qubits, pauli)` places letter `j` on `qubits[j]`).
The **rotation convention is `exp(-iθP)`**, which is Litinski's `P_θ`.

| Variant | Fields | Semantics (exact) |
|---|---|---|
| `Measurement` | `pauli, bit, qubits` | Projective measurement of the Hermitian Pauli `pauli` (sign included) on `qubits`. The outcome is written to classical bit `bit`: **bit = 1 iff eigenvalue −1** (runtime: `Bool(phase >> 1)` on the `project!` phase). The state collapses onto the eigenspace. |
| `Pauli` | `pauli, qubits` | Apply the Pauli operator `pauli`. Equal to `ExpHalfPiPauli` up to a global phase. Never produced by the compiler; the compiler passes cannot commute anything past it. |
| `ExpHalfPiPauli` | `pauli, qubits` | `exp(-iπ/2·P) = -iP`, i.e. a Pauli gate up to global phase. Parser: `x,y,z`. |
| `ExpQuatPiPauli` | `pauli, qubits` | `exp(-iπ/4·P)`. Parser: `s → +Z` (S = e^{iπ/4}·exp(-iπ/4 Z)), `sdg → −Z`, `h → Z,X,Z` (H ∝ Z_{π/4}X_{π/4}Z_{π/4}). |
| `ExpEighPiPauli` | `pauli, qubits` | `exp(-iπ/8·P)`. Parser: `t → +Z` (T = e^{iπ/8}·exp(-iπ/8 Z)), `tdg → −Z`. A negative sign on `P` is the inverse rotation. |
| `PrepMagic` | `qubit, qubits` | **Undefined.** Docstring is `TODO`, no pass, runtime or test gives it semantics (it appears only in `affectedqubits` and plotting). The verifier treats it as unsupported (result `UNKNOWN`, with a reason). The compiler never emits it. |
| `PauliConditional` | `control_pauli, control_qubits, target_pauli, target_qubits` | Unitary `C = Π₊ ⊗ I + Π₋ ⊗ P_t` with `Π± = (I ± P_c)/2`: apply the Pauli `P_t` when the control is in the −1 eigenspace of `P_c`. Equivalently `C = (I + P_c + P_t − P_cP_t)/2`. Parser: `cx q[c],q[t] → PauliConditional(Z,[c],X,[t])` (CNOT). See note (a). |
| `BitConditional` | `op, bit` | Apply `op` iff classical bit `bit` is 1. **Single bit, not a parity.** The runtime splices `op` in when the bit is `true` and deletes it when `false`. |

Notes:

- (a) The `PauliConditional` docstring says "(pi/2) Pauli rotation", which taken
  literally (`exp(-iπ/2 P_t) = -iP_t` under control) would add a controlled phase.
  The code does not do that. `remove_pauliconditional` uses
  `(P_c⊗P_t)_{π/4}·(P_t)_{−π/4}·(P_c)_{−π/4}`, which is `∝ C` above (C is Hermitian
  and unitary, so `C† = C`), and the parser uses it for CNOT. The verifier uses
  the projector definition above, cross-checked against dense matrices.
- The control and target qubits must be disjoint for `C` to be unitary. The verifier
  reports overlap as a well-formedness error.
- Hermiticity: `Measurement` and rotation Paulis must have phase ±1. A phase of ±i
  makes `exp(-iθP)` non-unitary, so the verifier reports it.

## 2. Public entry points

None of these are exported (the module had no `export`s before this work), so they
are reached as `PBCCompiler.<name>`.

| Function | File | Behaviour |
|---|---|---|
| `preprocess_circuit(c::Circuit)` | `preprocess.jl` | **Mutates `c` in place** and returns it. Returns early, unchanged, if `length(c) < 2`. Throws `ArgumentError` (reduction over an empty collection) if `c` has no `Measurement`. This is the static compiler. |
| `remove_pauliconditional(c)` | `preprocess.jl` | Pass 1: each `PauliConditional` becomes three `ExpQuatPiPauli`. |
| `group_nonclifford(c)` | `preprocess.jl` | Pass 2: bubbles each `ExpEighPiPauli` to the left past Clifford rotations (conjugating them). Stops at measurements, `Pauli`, and conditionals. |
| `merge_ops(c)` | `preprocess.jl` | Pass 3: merges adjacent same-Pauli rotations. |
| `remove_clifford(c)` | `preprocess.jl` | Pass 4: bubbles each `Measurement` to the left past Clifford rotations (conjugating them). The Cliffords are **not deleted here**; they end up after the measurements. |
| `remove_nonclifford(c)` | `preprocess.jl` | Pass 5: replaces each `ExpEighPiPauli` with a 4-op magic-state gadget (Section 5). |
| `remove_post_measurement(c)` | `preprocess.jl` | Pass 6: truncates everything after the last `Measurement`. This is where the Cliffords pushed past the final measurements are dropped. |
| `run(c, rt, input_state=nothing)` | `logic.jl` | Builds a `CompilerState` (preprocessing a copy) and executes it measurement by measurement, re-running `preprocess_circuit` on the *live* circuit whenever a conditional resolves or a random outcome injects a Clifford fix-up. `state.circuit` is therefore the dynamically rewritten circuit, **not** a static compiled program. |
| `parse_input(path)` | `io.jl` | OpenQASM 2 subset: `h s sdg t tdg x y z cx ccx measure` (`ccx` maps to a parity-controlled X, see Section 7). Qubit and bit indices are shifted to 1-based. |

**Where the compiled program lives:** `out = preprocess_circuit(copy(input))`. The
passes are callable one at a time, in the order above (this is exactly what
`preprocess_circuit` does after its length check), so per-pass validation (V8) is
possible without a hook. The order is duplicated rather than exposed as data,
though; see `docs/NOTES_FOR_COMPILER.md`.

**Semantic contract of the compiled program.** The compiled program is a
straight-line program of measurements, gadget corrections (`BitConditional`) and
possibly some unconditional Clifford rotations (merged T pairs) between gadgets.
Because pass 6 drops every operation after the last measurement, and because
measurements are conjugated (the measured *observable* changes even though the
outcome does not), the compiled program reproduces the **joint distribution of the
classical outcomes** of the input for every input state. It does **not** reproduce
the post-measurement quantum state. The verifier therefore defaults to
`semantics = :outcomes` (POVM equality per outcome record). It also offers the full
channel check `semantics = :channel`, which the compiler is not expected to pass.
This is a deliberate, documented deviation from spec A.1 ("same quantum channel").

**Classical register semantics.** A bit holds its most recent write
(`classical_register[bit] = result` in `do_quantum_step`). The visible outcome
record is the final value of each bit written by the input circuit.

## 3. Compiled output, op by op

(Output of `preprocess_circuit(copy(c))`, printed verbatim.)

### 3.1 Single T: `[ExpEighPiPauli(Z,[1]), Measurement(Z,1,[1])]`
```
[1] Measurement(pauli=+ ZZ, bit=2, qubits=[1, 2])
[2] Measurement(pauli=+ X, bit=3, qubits=[2])
[3] BitConditional(op=ExpQuatPiPauli(pauli=+ Z, qubits=[1]), bit=2)
[4] BitConditional(op=ExpHalfPiPauli(pauli=+ Z, qubits=[1]), bit=3)
[5] Measurement(pauli=+ Z, bit=1, qubits=[1])
```

### 3.2 CNOT then T: `[PauliConditional(Z,[1],X,[2]), ExpEighPiPauli(Z,[2]), Measurement(Z,1,[1]), Measurement(Z,2,[2])]`
```
[1] Measurement(pauli=+ ZZZ, bit=3, qubits=[1, 2, 3])
[2] Measurement(pauli=+ X, bit=4, qubits=[3])
[3] BitConditional(op=ExpQuatPiPauli(pauli=+ ZZ, qubits=[1, 2]), bit=3)
[4] BitConditional(op=ExpHalfPiPauli(pauli=+ ZZ, qubits=[1, 2]), bit=4)
[5] Measurement(pauli=+ Z_, bit=1, qubits=[1, 2])
[6] Measurement(pauli=+ ZZ, bit=2, qubits=[1, 2])
```

### 3.3 QFT example

**Deviation:** the repo has no QFT-4 example (no file or test mentions QFT). A
QFT-4 is also not exactly expressible in Clifford+T: it needs controlled-R₃ and R₄
phases, i.e. π/16 and π/32 rotations. In its place, here is the exact QFT-2 in the
parser's gate set (H, controlled-S, H, with the final swap done by relabelling the
measurements):

```
h q[0]; t q[0]; t q[1]; cx q[0],q[1]; tdg q[1]; cx q[0],q[1]; h q[1];
measure q[0] -> c[1]; measure q[1] -> c[0];
```
Compiled:
```
[1]  Measurement(pauli=+ XZ, bit=7, qubits=[1, 5])
[2]  Measurement(pauli=+ X, bit=8, qubits=[5])
[3]  BitConditional(op=ExpQuatPiPauli(pauli=+ X, qubits=[1]), bit=7)
[4]  BitConditional(op=ExpHalfPiPauli(pauli=+ X, qubits=[1]), bit=8)
[5]  Measurement(pauli=+ _ZZ, bit=5, qubits=[1, 2, 4])
[6]  Measurement(pauli=+ X, bit=6, qubits=[4])
[7]  BitConditional(op=ExpQuatPiPauli(pauli=+ _Z, qubits=[1, 2]), bit=5)
[8]  BitConditional(op=ExpHalfPiPauli(pauli=+ _Z, qubits=[1, 2]), bit=6)
[9]  Measurement(pauli=- XZZ, bit=3, qubits=[1, 2, 3])
[10] Measurement(pauli=+ X, bit=4, qubits=[3])
[11] BitConditional(op=ExpQuatPiPauli(pauli=- XZ, qubits=[1, 2]), bit=3)
[12] BitConditional(op=ExpHalfPiPauli(pauli=- XZ, qubits=[1, 2]), bit=4)
[13] Measurement(pauli=+ X_, bit=2, qubits=[1, 2])
[14] Measurement(pauli=+ _X, bit=1, qubits=[1, 2])
```
The benchmarks (V9) use an *approximate* QFT-n that keeps only the controlled-S
(R₂) rotations, which is exact in Clifford+T.

### 3.4 Repo fixture `test/fixtures/sample_input.qasm`
562 parsed ops (121 T/T†, 3 CNOT, 2 measurements) compile to 451 ops: 224
measurements (111 gadgets × 2 + 2 data), 222 `BitConditional`, and **5 unconditional
`ExpQuatPiPauli` left between gadgets** (merged T·T pairs that pass 4 cannot move
past the neighbouring T gates). The largest qubit index is 113 (2 data + 111 magic).

## 4. Utilities already available

- `affectedqubits(op)` / `affectedqubits(circuit)` (`affectedqubits.jl`): sorted unique
  qubit indices.
- `find_variant_indices(vec, T)`, `get_circuit_width(c)` (= max qubit index),
  `get_bit_number(c)` (= max measured bit), `paulis(op)`, `complete_paulis`,
  `check_commutation` (`preprocess.jl`, `pair_transformation.jl`).
- `traversal` (`traversal.jl`).
- `random_test_circuit(num_ops, num_qubits)` (`random_circuit.jl`): random ops on
  *all* qubits with control qubits always below target qubits.
- Existing tests: `test_affectedqubits`, `test_traversal`, `test_commutation`,
  `test_conjugation`, `test_compilation_pass`, `test_io`, `test_random_circuit`,
  `test_plotting`, Aqua, doctests, and JET (`@test_broken`).

Per spec 0.2 the verifier uses **none** of the transformation or commutation helpers.
It reimplements width, support, and Pauli conjugation from QuantumClifford
primitives. The resource estimator may use `affectedqubits`.

## 5. How T gates appear in the output: the gadget

Each `ExpEighPiPauli(P, Q)` becomes, contiguously (`gadgetize`):
```
Measurement(P ⊗ Z, b1, [Q; a])                     # joint measurement with magic qubit a
Measurement(X, b2, [a])                             # magic qubit measured in X
BitConditional(ExpQuatPiPauli(P, Q), b1)            # P_{π/4} correction if b1 = 1
BitConditional(ExpHalfPiPauli(P, Q), b2)            # P_{π/2} (Pauli) correction if b2 = 1
```
- Magic qubits are numbered `w+1, w+2, …` where `w = get_circuit_width` at pass 5.
  They are assigned in **reverse** program order (the last T gets `w+1`).
- Bits are `w + 2j − 1`, `w + 2j` for the j-th gadget counted from the end.
- **Magic state:** there is no preparation op in the IR. Magic qubits are implicitly
  prepared in `|m⟩ = T|+⟩ = (|0⟩ + e^{iπ/4}|1⟩)/√2` (`create_magic_state`: `|+⟩` then
  `pcT`), and are never reused.
- Verified exactly (dense, every branch): with `|m⟩` as above, the gadget implements
  `exp(-iπ/8·P)` on the data in every branch `(b1,b2)`, each with probability 1/4,
  for an arbitrary data state. Script output:
  ```
  b1=false b2=false  |c|^2=0.25  residual=1.0e-16
  b1=false b2=true   |c|^2=0.25  residual=1.15e-16
  b1=true  b2=false  |c|^2=0.25  residual=1.0e-16
  b1=true  b2=true   |c|^2=0.25  residual=1.04e-16
  ```
  The opposite convention `exp(+iπ/8 Z)` gives residual 0.5, so the convention is
  unambiguous.

## 6. Deviations from the spec, and design decisions forced by the code

1. **Semantics:** default `:outcomes` (Section 2), with `:channel` available.
2. **Magic qubits:** they are identified as the qubits with index greater than the
   input circuit's width, implicitly prepared in `|m⟩`. "Prepared before use" (V1)
   becomes "each magic qubit is first touched by exactly one joint measurement and
   then consumed by exactly one single-qubit measurement, with no other use".
3. **Symbolic Pauli frame sign:** it is tracked as X/Z bit-vectors that are affine
   over the outcome bits, **without a sign**. A per-branch global phase is
   unobservable (spec A.1), and the sign of a conjugated frame is in general
   *quadratic* in the outcome bits (for example S maps X→Y and Y→−X, so the sign
   depends on x·z). It is therefore not an affine function. The observable effect of
   the frame, sign flips on later measurements and rotations, is the commutation
   parity, which *is* affine, and that is what the normal form records.
4. **Conditional Paulis** are always folded into the symbolic Pauli frame (exact), so
   no `CondPauliItem` is ever emitted.
5. **Gadget proofs** use exact arithmetic in `Q(ζ₈)` (cyclotomic field), not
   floating point, so a lifted gadget is a proof, not a numerical check. The
   magic state is parameterised as `(|0⟩ + ζ₈^k|1⟩)/√2` with `k = 1` by default.
6. **Dense oracle:** it uses floating point with tolerance 1e-9, so by spec Section 3
   it can only refute (`NOT_EQUIVALENT`). Agreement is reported as `UNKNOWN` with the
   deviation in `stats`. `EQUIVALENT` comes only from the normal form.
7. **QFT-4** replaced by QFT-2 and approximate QFT (Section 3.3).

## 7. Compiler bugs found while mapping (owned by Yuhan; not fixed here)

All of these reproduce on `ed0efec`. They are also listed in
`docs/NOTES_FOR_COMPILER.md`, and the verifier test suite has a test for each
showing that it is detected.

1. **`PauliConditional` with control qubit index > target index** (`remove_pauliconditional`):
   the joint rotation is built as `ExpQuatPiPauli(cp⊗tp, sort(union(cq,tq)))`, which
   pairs the letters with the wrong qubits. `cx q[1],q[0]` compiles to measurements
   of `−XX` and `ZY` instead of `ZZ` and `Z`.
2. **`ExpEighPiPauli` with unsorted `qubits`** (`gadgetize`): it uses `paulis(op)` (in
   `op.qubits` order) together with `affectedqubits(op)` (sorted). A T on `XZ` over
   `[2,1]` compiles to a gadget measuring `X` on qubit 1 and `Z` on qubit 2.
3. **Gadget bit numbering** starts at the qubit width, not the bit count. An input
   that writes more bits than it has qubits (for example measuring one qubit twice)
   gets gadget bits that collide with its own bits.
4. **No measurements:** `preprocess_circuit` throws (`maximum` of an empty collection).
5. **`validate_CircuitOp`** throws `UndefVarError` (it interpolates undefined `$name`,
   `$p`, `$q`) instead of `PauliQubitMismatchError` for an empty `PauliConditional`.
6. **`ccx`** is parsed as `PauliConditional(ZZ, [c1,c2], X, [t])`, which is a
   parity-controlled X (it applies X when Z⊗Z = −1, i.e. c1 ⊕ c2 = 1), not a Toffoli.
