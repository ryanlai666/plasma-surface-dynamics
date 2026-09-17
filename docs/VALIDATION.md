# Local validation record

Validated on Windows 11 using Python 3.13.5, NumPy 2.1.3, SciPy 1.15.3, scikit-learn 1.6.1, Matplotlib 3.10.0, and MSYS2 GCC 15.2.0.

## Numerical and software checks

19 tests passed with the native library built. Checks include adsorption's analytical limit, no-flux/subthreshold limits, deposition sign and material accounting, seeded repeatability, event-budget enforcement, deterministic campaign sharding, nonoverlapping surrogate recipe splits, synthetic parameter recovery, Python/C++ mean-field agreement, and stochastic ensembles against exact expectations. Ensemble comparisons use a five-standard-error tolerance; they do not establish physical correctness against experiments.

The C++ tests exercise ALE, mixed etching/growth, pure deposition, and idle phases. Five parameterized deterministic cases agree within 1e-12 absolute tolerance. Native tests are skipped if the library has not been built; a passing Python-only run is not evidence that the C++ backend was tested.

## CPU benchmark

Workload: 2,048 sites, 10 four-phase ALE cycles, three seeds per backend, approximately 26,700 events per run. Both paths were warmed up before timing. Median complete-call wall times:

| Backend | Median time |
|---|---:|
| Python reference | 0.54815 s |
| C++17 | 0.001751 s |

Observed speedup: approximately **313×** for this workload on this machine. The C++ implementation also replaces array scans with constant-time eligible-site selection; this is not a compiler-only comparison. Random streams differ. The shortest measurements are sensitive to timing noise, and performance should be remeasured for the intended site count and event density. No universal speedup, GPU measurement, or cluster scaling claim is made.

Raw repeat timings, event counts, outputs, compiler flags, package versions, and hashes are in `outputs/benchmark.json`. Reproduce with `python -m plasma_surface.cli benchmark`. Rebuilding requires `python scripts/build_cpp.py`.

## ML and data checks

The 128-recipe C++ demo gave held-out synthetic EPC MAE approximately 0.0131 nm/cycle, versus 0.0548 for the training-mean baseline. High-energy holdout MAE was approximately 0.0355 nm/cycle. These numbers measure approximation of an uncalibrated model, not experimental predictive accuracy. The high-energy audit uses a separate training fit restricted below 65 eV.

The Si–HCl public ZIP was verified against its provider MD5; 22 structures were inventoried, with member SHA-256 hashes and original attribution. Downloaded structures do not calibrate the ALE model. SiNₓ curves remain hypothetical parameter studies.

Both generated overview and composition-comparison figures were visually inspected. The SLURM template and hosted CI configuration have not been executed on a cluster or remote CI service.
