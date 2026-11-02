# Verifier evaluation (V9)

All numbers below come from scripts in `benchmark/verification/`, run on this
machine (Apple Silicon, Julia 1.12.3, `-tauto`). Re-run them to regenerate.

## Mutation testing

`julia -tauto --project=test benchmark/verification/mutation_report.jl 300 2026`

- 300 random programs (1–4 data qubits), compiled by PBCCompiler, then each mutation
  kind applied to a copy of the compiled output.
- Ground truth for "channel-changing" is the dense oracle under the identity outcome
  map, which is PBCCompiler's contract.
- "detected" means `verify_compilation` returned `NOT_EQUIVALENT` with the default
  oracle stack.
- "stabilizer-only rate" is the same with the dense oracle disabled, so the refutation
  comes from the independent `GeneralizedStabilizer` oracle.
- "false refutations" counts `NOT_EQUIVALENT` on a mutation that does not change the
  behaviour, from either path.

| mutation | applied | channel-changing | detected | rate | by dense | by stabilizer | stabilizer-only rate | false refutations |
|---|---|---|---|---|---|---|---|---|
| flip_condition_bit | 192 | 98 | 98 | 100.0% | 98 | 0 | 100.0% | 0 |
| drop_correction | 198 | 113 | 113 | 100.0% | 113 | 0 | 100.0% | 0 |
| flip_sign | 283 | 188 | 188 | 100.0% | 188 | 0 | 100.0% | 0 |
| swap_anticommuting | 153 | 92 | 92 | 100.0% | 92 | 0 | 100.0% | 0 |
| change_letter | 283 | 246 | 246 | 100.0% | 246 | 0 | 100.0% | 0 |
| change_angle | 198 | 106 | 106 | 100.0% | 106 | 0 | 100.0% | 0 |
| drop_gadget | 199 | 115 | 115 | 100.0% | 115 | 0 | 100.0% | 0 |

**Misses: none.** All channel-changing mutants are detected (100% for every mutation
kind) and there are no false refutations. In these small instances the dense oracle
is reached first, so "by dense" covers everything. The stabilizer-only column shows
that the independent oracle alone also catches every one. The mutations that do not
change behaviour (applied minus channel-changing) are excluded from the rates. One
verified example is a sign flip on a measurement whose outcome is uniformly random
for every input (see the V3 notes in `progress.md`).

The test suite runs a smaller version of this (`test_mutation.jl`, 60 programs). It
asserts 100% detection, zero false refutations, at least 90% stabilizer-only
detection, and at least one mutant of every kind.

## Timing

`julia -tauto --project=test benchmark/verification/timing.jl`, then `plot_timing.jl`
for the PNGs (`timing_vs_{qubits,t_count,outcomes}.png`, git-ignored). The CSV is
`benchmark/verification/results/timing.csv`. These are normal-form verification times
(oracles off), best of 3. Every program below verified `EQUIVALENT`.

| family | param | qubits | T count | compiled outcomes | time |
|---|---|---|---|---|---|
| approx_qft | 2 | 2 | 3 | 13 | 0.05 ms |
| approx_qft | 4 | 4 | 9 | 31 | 0.18 ms |
| approx_qft | 6 | 6 | 15 | 49 | 0.35 ms |
| approx_qft | 8 | 8 | 21 | 67 | 0.56 ms |
| approx_qft | 10 | 10 | 27 | 85 | 0.99 ms |
| approx_qft | 12 | 12 | 33 | 103 | 1.37 ms |
| approx_qft | 14 | 14 | 39 | 121 | 2.07 ms |
| approx_qft | 16 | 16 | 45 | 139 | 2.71 ms |
| cuccaro_adder | 1 | 4 | 14 | 119 | 0.25 ms |
| cuccaro_adder | 2 | 6 | 28 | 235 | 1.16 ms |
| cuccaro_adder | 3 | 8 | 42 | 351 | 3.21 ms |
| cuccaro_adder | 4 | 10 | 56 | 467 | 5.47 ms |
| cuccaro_adder | 5 | 12 | 70 | 583 | 10.61 ms |
| random_t0.1 | 2 | 2 | 2 | 24 | 0.05 ms |
| random_t0.3 | 2 | 2 | 8 | 20 | 0.10 ms |
| random_t0.1 | 4 | 4 | 5 | 46 | 0.15 ms |
| random_t0.3 | 4 | 4 | 6 | 46 | 0.20 ms |
| random_t0.1 | 8 | 8 | 9 | 84 | 0.44 ms |
| random_t0.3 | 8 | 8 | 18 | 80 | 0.66 ms |
| random_t0.1 | 12 | 12 | 20 | 138 | 0.99 ms |
| random_t0.3 | 12 | 12 | 33 | 124 | 2.31 ms |
| random_t0.1 | 16 | 16 | 27 | 172 | 2.28 ms |
| random_t0.3 | 16 | 16 | 35 | 178 | 3.02 ms |
| fixture | 0 | 2 | 121 | 562 | 32.58 ms |

Empirically, time grows like a power of the T count: about T^1.5 for the QFT family
(T 3→45, 0.05→2.71 ms) and about T^2.3 for the adders (T 14→70, 0.25→10.6 ms).
The canonicalization's pairwise DAG and merge scans are quadratic in the number of
items. The 121-T
fixture (224 compiled measurements on 2 data qubits) is the slowest, at about 30 ms,
because canonicalization re-absorbs the merged Clifford pairs. BenchmarkTools
numbers for three of these are in `benchmark/benchmarks.jl` (group
`"verification"`): fixture 22.2 ms, approx QFT-12 1.20 ms, Cuccaro adder k=3 2.91 ms.

## Long cross-check

`julia -tauto --project=test benchmark/verification/long_crosscheck.jl 3000 1`:
3000 random (input, compiled-or-mutated) pairs with up to 8 total qubits. 1434 were
proved `EQUIVALENT` by the normal form and 853 refuted by the dense oracle, with
**0 contradictions**, in 188 s.
