# Simulation and numerical verification report

Generated 2026-09-17T03:36:59.469280+00:00 using `python -m plasma_surface.validation`.

**Numerical checks: PASS.** This verifies implementation and internal consistency, not agreement with measured plasma etching.

- Regression test output: [tests.txt](tests.txt).
- 4,096 recipes evaluated by both backends; maximum EPC difference: 2.22e-16 nm/cycle.
- Six stochastic ensembles: all means within five standard errors of the exact mean-field expectation: True.
- Subthreshold/chemical-only controls and dose saturation: True.
- Synthetic surrogate held-out MAE: 0.00284 nm/cycle; constant baseline: 0.05223; high-energy holdout: 0.02490.
- Median kMC times: Python 2.640481 s, C++ 0.002218 s; workload-specific speedup 1190.5×, including the improved site-selection algorithm.
- Public DFT archive: checksum_verified. Structures are not experimental etch-rate labels.

![Convergence and saturation](verification.png)

## Results to inspect

- [Python ALE figure](python/overview.png) and [C++ ALE figure](cpp/overview.png).
- [Si/SiNx hypothetical scenarios](materials/comparison.png).
- [Stochastic convergence](convergence.csv), [raw ensemble samples](ensemble_samples.csv), and [half-cycle controls](half_cycle_controls.csv).
- [Surrogate metrics](ml/surrogate_metrics.json) and [predictions](ml/surrogate_predictions.csv).
- [Benchmark](benchmark.json), [summary](summary.json), and [provenance](manifest.json).

SiNx composition cards remain hypothetical. Default rates are uncalibrated. Statistical tolerances were fixed at five standard errors before evaluating ensembles; no seeds were discarded. Finite-site scatter should scale approximately as N^(-1/2), but this is not proof of experimental accuracy. Public geometry checks do not validate reaction barriers. No experimental validation or HPC scaling was performed.
