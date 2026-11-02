# [Resource estimation](@id resource_estimation)

```@meta
CurrentModule = PBCCompiler.ResourceEstimation
```

`PBCCompiler.ResourceEstimation` turns a compiled PBC program into physical
resources:

- physical qubits (data, syndrome, surgery/adapter ancilla, factories, routing);
- time (in code cycles and seconds);
- logical failure probability, with per-block and per-measurement breakdowns;
- the chosen code distances.

```julia
using PBCCompiler, QuantumClifford
using PBCCompiler: CircuitOp, preprocess_circuit
using PBCCompiler.ResourceEstimation

compiled = preprocess_circuit([CircuitOp.ExpEighPiPauli(P"ZZ", [1, 2]),
                               CircuitOp.Measurement(P"Z", 1, [1]), CircuitOp.Measurement(P"Z", 2, [2])])
layout = pack_sequential(2, [SurfaceCodeModel(11), SurfaceCodeModel(11)])
hw = HardwareModel(:all_to_all; t_1q = 25e-9, t_2q = 50e-9, t_meas = 500e-9, t_reaction = 1e-6, p_phys = 1e-3)
est = estimate_resources(compiled, layout; hardware = hw,
                         adapter_model = AnalyticAdapterModel(cycle_time_factor = 1.0),
                         optimize_distances = true, error_budget = 1e-2)
```

Every default constant comes with a published source. Anything that has no source
(all hardware times and rates, the qLDPC in-block measurement costs, the
auxiliary-graph sizes, the adapter cycle-time factor) is a **required** argument.

## Model of computation

The compiled program is a sequence of Pauli product measurements, magic-state
gadgets, Clifford corrections and Pauli-frame updates. Each operation is costed as
follows:

| operation | treatment |
|---|---|
| measurement on one block | `intra_measurement_cost(block, weight)`: ancilla and rounds |
| measurement on ≥ 2 blocks | the adapter model (analytic or constructed) |
| measurement on magic qubits only | one code cycle, no extra qubits (the factory model covers the rest) |
| π/4 rotation (conditional or not), `PauliConditional` | one product measurement on the same support (via a \|Y⟩ resource, GoSC Sec. 1) |
| π/2 rotation, Pauli (conditional or not) | Pauli frame: tracked classically, free |
| joint measurement touching a magic qubit, or an unlifted π/8 | consumes one magic state |

Operations are scheduled on a dependency DAG. Two operations depend on each other
when their Paulis anticommute, or classically when one reads a bit the other wrote,
which adds `t_reaction`. Scheduling is list scheduling:
- each block runs at most `max_parallel_surgeries` operations at once;
- a cross-block operation holds every block it touches;
- a magic-consuming operation waits for the next factory output (storage is
  unlimited).

## Cost formulas and parameters

**Surface code** ([`SurfaceCodeModel`](@ref)): `[[d², 1, d]]` with `d² − 1` syndrome
qubits.
- Cycle time `3·t_1q + 4·t_2q + t_meas`: the 8-step cycle of Fowler et al., PRA 86,
  032324 (2012).
- `p_L(p, d) = 0.1·(100p)^((d+1)/2)` per cycle: Litinski, *A Game of Surface Codes*
  (GoSC), Eq. (10).
- In-block (lattice surgery) measurement: `d` rounds. Ancilla
  `ancilla_tiles_per_qubit · weight · (2d² − 1)`; the default of one routing tile per
  patch is a documented modelling assumption.

**Bivariate bicycle codes** ([`BBCodeModel`](@ref)), from Bravyi et al., Nature 627,
778 (2024):
- presets `:gross` `[[144,12,12]]` and `:two_gross` `[[288,12,18]]` (polynomials from
  their Table 3);
- `n` check qubits;
- cycle time `7·t_2q + t_1q + t_meas` (a depth-7 CNOT circuit);
- `p_L = p^(d_circ/2)·exp(c₀ + c₁p + c₂p²)` with their Table 6 coefficients. It
  reproduces their Table 1.
The in-block measurement cost and `aux_graph_qubits(w)` are required.

**Any QuantumClifford code** ([`GenericCSSModel`](@ref)): `n` and `k` are read from
the code, and syndrome ancillas default to one per check. `d` (or
`compute_distance = true`), `pl`, `cycle_layers`, and the measurement costs are
required.

**Logical error models** ([`PLModel`](@ref)): [`AnalyticPL`](@ref) (`A·(p/p_th)^((d+1)/2)`),
[`FitPL`](@ref) (any closure `(p, d) -> p_L`), [`BBFitPL`](@ref), and [`TablePL`](@ref).

**Analytic adapters** ([`AnalyticAdapterModel`](@ref)) for a measurement on blocks
`B₁…B_a` with weights `wᵢ`:

    d_joint  = min d(Bᵢ)
    ancilla  = Σ aux_qubits_per_block(Bᵢ, wᵢ) + (a − 1)·adapter_qubits_per_join(d_joint)
    rounds   = rounds(d_joint)                                 # default d_joint
    round    = cycle_time_factor · max cycle_time(Bᵢ)          # factor required

The default `adapter_qubits_per_join(d) = 2d − 1` is the `w` qubits plus `w − 1`
checks of a width-`d` adapter (Swaroop, Jochym-O'Connor, Yoder, arXiv:2410.03628,
Table III).

**Constructed adapters** ([`ConstructedAdapterModel`](@ref)) build the actual merged
code with the universal-adapter construction (`QuantumClifford.ECC.build_adapter`,
from the `riteshbhirud/QuantumClifford.jl` fork, PR #742), when it is loaded. The
result is: ancilla = new qubits + new checks, rounds `d_joint`, cached per logical
pair, and reported next to the analytic number. Measurements it cannot build fall
back to the analytic model, and **the fallback and its reason are recorded** in
`est.notes`.

**Factories** ([`FifteenToOneFactory`](@ref), [`litinski_factory`](@ref)):
- GoSC's 15-to-1 block uses `11 · 2d²` qubits, takes `11d / (1 − p)¹⁵` code cycles
  per state, and has output error `35p³`;
- or fixed numbers from Litinski, *Magic State Distillation: Not as Costly as You
  Think*, Table 1. A note is recorded if the hardware `p` differs from the table's.
The number of factories is chosen by scheduling each candidate and minimizing
makespan (or space-time volume, `factory_objective = :spacetime`).

**Hardware** ([`HardwareModel`](@ref)): connectivity `:all_to_all`, `:grid`,
`:heavy_hex` or a custom `Graphs.SimpleGraph`. For anything other than all-to-all,
each block's Tanner graph is embedded by seeded simulated annealing
([`embed_tanner_graph`](@ref)). The non-native edges then cost either
- `:bridge`: `L − 1` routing qubits each, plus `2(max L − 1)` CNOT layers per round;
  or
- `:swap`: `6(max L − 1)` layers per round.
The overhead is reported separately (`qubit_breakdown[:routing]`).

**Error budget.** The failure probability is

    Σ_blocks p_L(block) · cycles alive  +  Σ_operations Σ_their blocks p_L · rounds  +  #magic · factory error

With `optimize_distances = true`, every code family starts at its smallest distance
(surface codes share one distance; BB moves gross → two-gross). The family
contributing the most failure is then increased, and everything is re-scheduled,
until the budget is met. If it cannot be met (for example the magic-state error alone
exceeds the budget), `meets_budget` is false and a note explains why.

## Building a layout

- [`pack_sequential`](@ref)`(nqubits, models)` fills blocks in qubit order.
- [`pack_by_interaction`](@ref)`(compiled, models)` groups qubits that are measured
  together. It places qubits greedily by affinity, then refines with Kernighan–Lin
  style moves and swaps, and returns the layout and the remaining cross-block
  measurement weight.
- Or build `LogicalLayout(blocks, qubit_to_block, qubit_within_block)` yourself.
  [`validate_layout`](@ref) lists every problem.

## Plugging in threshold data

Write your simulation results as CSV with the header
`code,d,p,p_L,ci_low,ci_high`, and pass `pl = TablePL("results.csv"; code = "mycode")`
to a code model. Values are interpolated in log(p)–log(p_L) at each tabulated
distance. Querying an untabulated distance, or a `p` outside the tabulated range,
throws instead of extrapolating. For a fitted formula use `FitPL((p, d) -> ...)`.

## Reading the report

`show(est)` prints:
- the qubit breakdown;
- time in seconds and in cycles of the slowest block's syndrome cycle;
- the logical depth;
- the failure probability split into blocks, measurements and magic;
- the chosen distances, the magic-state counts by gadget type, the number of
  factories, and the measurement-weight histogram;
- the critical path (if the makespan exceeds the critical path's time, resources
  were the bottleneck);
- per-block p_L, failure and utilization;
- every note.

`est.per_measurement` lists, per operation: class, blocks, weights, ancilla, rounds,
start time, failure contribution, cost source, and (with constructed adapters) the
analytic ancilla for comparison. [`to_dict`](@ref) gives the whole estimate as
JSON-friendly `Dict`s.

## Limitations

- These are first-order models. The surface-code routing (one tile per patch per
  product) is pessimistic for compact layouts. The GoSC comparison in
  `docs/dev/estimation_validation.md` gives 25% more qubits than GoSC's compact
  block, with the same time.
- Measurement failure is `p_L × rounds` per involved block. Merged-code distances are
  not computed (`d_joint` is used as a bound).
- The constructed adapter model handles two blocks with one X or Z logical each, and
  needs the fork. Other cases fall back (recorded).
- The magic-state storage is unlimited and free, and magic-only measurements take
  one cycle.
- Routing overhead is computed for code blocks only, not for merged (adapter) codes.
- The scheduler is greedy list scheduling (no backfilling), and its dependency DAG is
  quadratic in the number of operations.

## API

```@autodocs
Modules = [PBCCompiler.ResourceEstimation]
Private = false
```
