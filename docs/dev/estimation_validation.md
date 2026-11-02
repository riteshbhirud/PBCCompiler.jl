# Resource estimator validation (E10)

## 1. Hand-computable unit tests

These are in `test/resource_estimation/test_estimate.jl` (tags `:e9`, `:e10`), with
the expected numbers written out in the test:

- **two blocks, one cross-block and one intra-block measurement:** 73 qubits
  (18 data + 16 syndrome + 39 peak ancilla), 5.49 µs, failure 0.027;
- **a single block:** a single d = 5 measurement;
- **a classical dependency:** `t_reaction` lands on the critical path.

## 2. Property tests

These are in the same file:
- the total qubit count equals the sum of the breakdown;
- increasing d never decreases the qubit count, and memory plus surgery failure
  decreases;
- adding a measurement never decreases the makespan;
- routing is zero on all-to-all hardware;
- distance optimization returns the smallest distance that meets the budget.

## 3. Reference comparison: Litinski, "A Game of Surface Codes", Sec. 4

Script: `benchmark/estimation/gosc_comparison.jl`. Raw output:
`benchmark/estimation/results/gosc_comparison.txt`.

The reference is GoSC's minimal setup at p = 10⁻⁴: 100 logical qubits, 10⁸ T gates,
a compact data block of 153 tiles plus one 15-to-1 block of 11 tiles. With d = 13
this is 164 · 2 · 13² ≈ **55,400 qubits**, 13·11·10⁸ code cycles = **3.97 h** at 1 µs,
and a final (memory) error of **0.2%** (19.8% at d = 11).

**Assumptions made to match GoSC.** None of these is tuned to the result.
- One surface-code patch per logical qubit, and a 1 µs code cycle (each of the
  8 cycle steps takes 1/8 µs).
- No reaction time.
- Lattice surgery with one routing tile per involved patch, and no adapter qubits.
- One GoSC 15-to-1 factory.
- Random π/8 Pauli-product rotations on all 100 qubits, compiled by PBCCompiler.
  Results are extrapolated linearly from N = 400 T gates. The per-T time and failure
  are flat for N = 100, 200 and 400, which the script checks.

    Scaling check (per-T quantities should not depend on N):
      N_T= 100 d=11  qubits= 47247  cycles/T=  121.5  failure/T=1.416e-09
      N_T= 100 d=13  qubits= 66063  cycles/T=  143.6  failure/T=5.133e-11
      N_T= 200 d=11  qubits= 47729  cycles/T=  121.3  failure/T=1.415e-09
      N_T= 200 d=13  qubits= 66737  cycles/T=  143.4  failure/T=5.130e-11
      N_T= 400 d=11  qubits= 49416  cycles/T=  121.3  failure/T=1.414e-09
      N_T= 400 d=13  qubits= 69096  cycles/T=  143.3  failure/T=5.130e-11
    
    Extrapolation to 10^8 T gates (from the largest N):
      d=11: 49416 physical qubits (data 12100, syndrome 12000, surgery ancilla peak 22654, factory 2662), 1.21e+10 cycles = 3.37 h, failure 0.141
      d=13: 69096 physical qubits (data 16900, syndrome 16800, surgery ancilla peak 31678, factory 3718), 1.43e+10 cycles = 3.98 h, failure 0.00513
    
    Reference (GoSC Sec. 4.3, d = 13): 55,400 qubits, 13·11·10^8 = 1.43e10 cycles = 3.97 h, failure 0.2% (d = 11: 19.8%; their failure counts memory errors on all 164 tiles; the T-gate error budget is separate)

| quantity (d = 13) | this estimator | GoSC | gap |
|---|---|---|---|
| time | 1.43·10¹⁰ cycles = 3.98 h | 1.43·10¹⁰ cycles = 3.97 h | +0.2%: our factory time includes the (1−p)¹⁵ success probability |
| physical qubits | 69,096 | 55,400 | +25%, see below |
| failure probability | 0.51% (of which magic 0.35%) | 0.2% memory, plus a separate ≤1% T-gate budget | memory + surgery ≈ 0.16% |
| failure at d = 11 | 14.1% | 19.8% | both choose d = 13 for a 1% memory budget |

**Gaps, explained (not tuned away).**

- *Qubits.* Data and syndrome (2·100·d²) plus the factory (11·2d²) give 37,418. The rest
  is the peak lattice-surgery routing (31,678): our model charges one tile per patch
  in each product. A random product touches about 94 of the 100 patches, and the peak
  is a maximum over all products, so it creeps up slowly with N (66k at N = 100, 69k
  at N = 400). GoSC's compact block uses a fixed 53-tile routing region for 100
  qubits. Our per-product model is pessimistic for this layout. A
  `ancilla_tiles_per_qubit` of about 0.56 would reproduce the compact block, but it is
  left at its documented default rather than fitted.
- *Failure.* Our total includes the magic-state error, 35p³ · 10⁸ = 0.35%, which GoSC
  budgets separately. Without it, memory plus surgery is about 0.16%, against GoSC's
  0.2%. GoSC charges all 164 tiles as idle memory for the whole run; we charge 100
  patches all the time plus the routing tiles only while they are in use.

## 4. Analytic (E4) vs constructed (E5) adapter costs

Script: `benchmark/estimation/adapter_agreement.jl`, run in an environment with the
`riteshbhirud/QuantumClifford.jl` fork (master, PR #742).

"Analytic join" is the `2d − 1` adapter term of arXiv:2410.03628 alone. "Analytic
total" adds the block's auxiliary-graph assumption: one routing tile per patch for
surface codes, and 0 for BB here, so that only the join term is compared. A value of
−1 means "not computed (fallback)".

| pair | d_joint | constructed ancilla | analytic join (2d−1) | analytic total | construction note |
|---|---|---|---|---|---|
| surface d=3 × d=3, ZZ | 3 | 15 | 5 | 39 | merged code: +7 qubits, +8 checks, max check weights X=5 Z=4 (1.70 s) |
| surface d=3 × d=3, XZ | 3 | 15 | 5 | 39 | merged code: +7 qubits, +8 checks, max check weights X=5 Z=4 (0.02 s) |
| surface d=5 × d=5, ZZ | 5 | 27 | 9 | 107 | merged code: +13 qubits, +14 checks, max check weights X=5 Z=4 (0.02 s) |
| surface d=5 × d=5, XZ | 5 | 27 | 9 | 107 | merged code: +13 qubits, +14 checks, max check weights X=5 Z=4 (0.00 s) |
| surface d=7 × d=7, ZZ | 7 | 39 | 13 | 207 | merged code: +19 qubits, +20 checks, max check weights X=5 Z=4 (0.04 s) |
| surface d=7 × d=7, XZ | 7 | 39 | 13 | 207 | merged code: +19 qubits, +20 checks, max check weights X=5 Z=4 (0.00 s) |
| surface d=9 × d=9, ZZ | 9 | 51 | 17 | 339 | merged code: +25 qubits, +26 checks, max check weights X=5 Z=4 (0.08 s) |
| surface d=9 × d=9, XZ | 9 | 51 | 17 | 339 | merged code: +25 qubits, +26 checks, max check weights X=5 Z=4 (0.01 s) |
| surface d=11 × d=11, ZZ | 11 | 63 | 21 | 503 | merged code: +31 qubits, +32 checks, max check weights X=5 Z=4 (0.10 s) |
| surface d=11 × d=11, XZ | 11 | 63 | 21 | 503 | merged code: +31 qubits, +32 checks, max check weights X=5 Z=4 (0.01 s) |
| gross × gross, Z̄₁Z̄₁ | 12 | 241 | 23 | 23 | merged code: +120 qubits, +121 checks, max check weights X=8 Z=8 (0.62 s) |
| gross × gross, Z̄₁Z̄₂ | 12 | 23 | 23 | -1 | constructed adapter fallback: the construction failed: DimensionMismatch: build_adapter: logical supports have unequal lengths (26 vs 38) (0.01 s) |
| gross × two-gross, Z̄₁Z̄₁ | 12 | 23 | 23 | -1 | constructed adapter fallback: the construction failed: DimensionMismatch: build_adapter: logical supports have unequal lengths (26 vs 50) (0.23 s) |
| two-gross × two-gross, Z̄₁Z̄₁ | 18 | 557 | 35 | 35 | merged code: +278 qubits, +279 checks, max check weights X=8 Z=12 (0.30 s) |

Observations:

- **Surface code, constructed.** The constructed adapter costs exactly 6d − 3 ancilla.
  Of those, 3d − 2 are new qubits. That is consistent with (d − 1) + (d − 1)
  auxiliary-graph path edges plus a width-d adapter, an inferred split that was not
  inspected edge by edge. The other 3d − 1 are new checks (one ancilla each). The analytic
  join term (2d − 1) covers only the adapter itself.
- **Surface code, analytic total.** The analytic total (one routing tile per patch,
  2d² − 1 each) is far more conservative. That tile is a lattice-surgery routing
  assumption, not an adapter construction.
- **BB codes.** Using `logz_ops`' representatives (weights 26 and 50, not minimal),
  the constructed merged codes add 241 ancilla for gross × gross and 557 for
  two-gross × two-gross, an order of magnitude above `2d − 1`. The auxiliary graph
  dominates. The analytic BB costs need a realistic `aux_graph_qubits(w)` from the
  user.
- **Findings in the adapter code (PR #742):**
  - its *inter-code* `build_adapter` rejects logical supports of unequal length
    (26 vs 38 for gross Z̄₁ × Z̄₂). The paper's Sec. VII.C uses the minimum width for
    unequal logicals, and the PR does this on its intra-code path only. The estimator
    records the failure and falls back to the analytic model.
  - two-gross × two-gross reaches Z-check weight 12, above the ≤ 8 of the paper's
    examples. That may be due to the non-minimal representative and was not
    investigated further.
