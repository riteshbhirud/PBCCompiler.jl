# PBCCompiler.jl

Tools for Pauli-based computation (PBC). The package compiles Clifford+T circuits
into sequences of Pauli product measurements with magic-state injection and
classical feed-forward, and simulates running them.

```julia
using PBCCompiler, QuantumClifford
using PBCCompiler: CircuitOp, preprocess_circuit

input = [CircuitOp.PauliConditional(P"Z", [1], P"X", [2]), CircuitOp.ExpEighPiPauli(P"Z", [2]),
         CircuitOp.Measurement(P"Z", 1, [1]), CircuitOp.Measurement(P"Z", 2, [2])]
compiled = preprocess_circuit(copy(input))
```

## Verifying compilation

`PBCCompiler.Verification` checks that a compiled program implements the same
logical operation as its input, without reusing any of the compiler's own
transformation code.

```julia
r = verify_compilation(input)          # compile, then verify
r.status                               # EQUIVALENT / NOT_EQUIVALENT / UNKNOWN
verify_passes(input)                   # translation validation of each compiler pass
```

The method layers are:
- static well-formedness checks;
- exact proofs of every magic-state gadget (arithmetic in Q(ζ₈));
- a canonical rotation normal form that *proves* equivalence;
- dense and stabilizer-state oracles that *refute* it.

A refutation comes with a minimized, dense-confirmed counterexample. See the
[verification docs page](docs/src/verification.md).

## Resource estimation

`PBCCompiler.ResourceEstimation` turns a compiled program into:
- physical qubits (data, syndrome, surgery/adapter ancilla, factories, routing);
- time;
- logical failure probability, with per-block and per-measurement breakdowns.

```julia
using PBCCompiler.ResourceEstimation
hw = HardwareModel(:all_to_all; t_1q = 25e-9, t_2q = 50e-9, t_meas = 500e-9, t_reaction = 1e-6, p_phys = 1e-3)
layout = pack_sequential(2, [SurfaceCodeModel(11), SurfaceCodeModel(11)])
est = estimate_resources(compiled, layout; hardware = hw,
                         adapter_model = AnalyticAdapterModel(cycle_time_factor = 1.0), optimize_distances = true)
```

It includes:
- **code models:** surface codes, bivariate-bicycle codes (gross and two-gross), and
  any QuantumClifford code;
- **adapters:** analytic cross-block costs, or exact ones from the universal-adapter
  construction (`ConstructedAdapterModel`, which needs a QuantumClifford with
  `build_adapter`);
- **hardware:** connectivity embedding;
- **scheduling:** a dependency-aware scheduler, magic-state factories, and an
  error-budget-driven distance search.

Every default constant has a published source; the rest are required arguments.
See the [resource estimation docs page](docs/src/resource_estimation.md).

## Development notes

- `docs/dev/codebase_map.md`: exact semantics of the IR.
- `docs/NOTES_FOR_COMPILER.md`: compiler bugs found by the verifier, and requested
  hooks.
- `docs/dev/progress.md`: the milestone log.
- Run the tests with `julia -tauto --project -e 'using Pkg; Pkg.test("PBCCompiler")'`.
