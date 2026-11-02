# Notes for the compiler owner

The verification and resource-estimation work does not modify anything in the
compilation pipeline. This file collects (1) hooks that would make verification
easier, and (2) compiler bugs the verifier found, each with a minimal reproduction.
Every bug below reproduces on upstream `ed0efec`. The verifier's test suite has a
test for each one asserting that it is detected, so when a bug is fixed that test
will start failing and should be flipped.

## Requested hooks

### H1. Expose the pass pipeline as data
`preprocess_circuit` hard-codes its six passes. `verify_passes` (per-pass
translation validation) currently duplicates that list in
`src/verification/passwise.jl`. If the pipeline changes, that copy silently
drifts. A guard test (`verify_passes` reproduces `preprocess_circuit` exactly)
fails if that happens.

Requested: something like
```julia
const PREPROCESS_PASSES = (remove_pauliconditional, group_nonclifford, merge_ops,
                           remove_clifford, remove_nonclifford, remove_post_measurement)
preprocess_circuit(c) = (length(c) < 2 ? c : (foreach(p -> p(c), PREPROCESS_PASSES); c))
```
so the verifier can iterate over `PBCCompiler.PREPROCESS_PASSES`.

### H2. Make magic-state preparation explicit (optional)
Magic qubits are implicit: any qubit beyond the input width is assumed to start in
`T|+⟩`. An explicit preparation op (the unused `PrepMagic`, once it has defined
semantics) would let the verifier check the magic state instead of assuming it.

## Bugs found

### B1. `PauliConditional` with control index > target index compiles wrong
`remove_pauliconditional` builds the joint rotation as
`ExpQuatPiPauli(cp ⊗ tp, sort(union(cq, tq)))`. `cp ⊗ tp` is ordered (control,
target), but the qubit list is sorted, so the letters land on the wrong qubits
whenever a control qubit has a larger index than a target qubit.
```julia
c = Circuit([PauliConditional(P"Z",[2],P"X",[1]), Measurement(P"Z",1,[1]), Measurement(P"Z",2,[2])])
preprocess_circuit(c)   # measures -XX (bit 1) and +ZY (bit 2); correct is Z₁Z₂ (bit 1) and Z₂ (bit 2)
```
Any QASM file with `cx q[j],q[i]` for `j > i` hits this.
Fix: build the joint Pauli over `vcat(cq, tq)` (unsorted) and let `embed` order it.

### B2. `ExpEighPiPauli` with unsorted `qubits` produces a wrong gadget
`gadgetize` takes `P = paulis(op)` (in `op.qubits` order) and
`Q = affectedqubits(op)` (sorted), then pairs them.
```julia
c = Circuit([ExpEighPiPauli(P"XZ",[2,1]), Measurement(P"Z",1,[1]), Measurement(P"Z",2,[2])])
preprocess_circuit(c)   # gadget measures X on q1 and Z on q2 (should be Z on q1, X on q2)
```
Fix: use `op.qubits` instead of `affectedqubits(op)`.

### B3. Gadget bit indices can collide with the program's own bits
`remove_nonclifford` numbers gadget bits from the *qubit* width, not from the
largest bit index already in use.
```julia
c = Circuit([Measurement(P"Z",1,[1]), ExpEighPiPauli(P"X",[1]), Measurement(P"Z",2,[1])])
preprocess_circuit(c)   # the gadget's joint measurement writes bit 2, the same bit as the final data measurement
```
Fix: start gadget bits at `max(get_bit_number(c), width) + 1`.

### B4. Circuits with no measurement crash
`remove_post_measurement` calls `maximum` on an empty index list, so
`preprocess_circuit([ExpQuatPiPauli(P"Z",[1]), ExpQuatPiPauli(P"X",[1])])` throws
`ArgumentError`.

### B5. `validate_CircuitOp` throws `UndefVarError`
The empty-`PauliConditional` branch interpolates `$name`, `$p`, `$q`, which are not
defined there, so the intended `PauliQubitMismatchError` is never thrown.

### B6. `ccx` is not a Toffoli
`parse_input` maps `ccx a,b,t` to `PauliConditional(ZZ, [a,b], X, [t])`. That
applies X when `Z_a Z_b = −1`, i.e. when `a ⊕ b = 1` (a parity-controlled X), not
when `a ∧ b = 1`. A Toffoli is not Clifford, so it needs a T-gate decomposition.
