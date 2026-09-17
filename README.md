# Plasma Surface Lab

A local-computer research and portfolio project for **silicon atomic layer etching (ALE)**, with exploratory **SiNₓ composition scenarios**, Python analysis, and an optional C++17 solver.

Research question: **How do surface modification, competing removal pathways, and uncertain reaction rates determine etch-per-cycle saturation and the usable ALE energy window?**

This is a working scientific software foundation, not an experimentally validated process simulator. Default rates are illustrative. SiNₓ cards contain explicitly invented sensitivity hypotheses, not fitted composition-dependent material properties. Public DFT structures are kept separate from synthetic model outputs.

## Run locally

Python 3.11+; no GPU, paid data service, or cluster required. From this directory:

```powershell
python -m venv .venv
.venv\Scripts\python -m pip install -e ".[dev]"
.venv\Scripts\python -m plasma_surface.cli demo
```

If your active Python already has NumPy, SciPy, scikit-learn, Matplotlib, and pytest:

```powershell
python -m plasma_surface.cli demo
python -m pytest -q --basetemp=outputs/pytest-run
```

`outputs/demo/` contains phase-resolved mean-field and kMC results, a 128-recipe synthetic campaign, held-out ML predictions, metrics, a figure, and a provenance manifest.

## C++ acceleration

A C++17 GCC/Clang compiler is optional. The current Windows machine has MSYS2 GCC. Build locally; no Python binding package is needed:

```powershell
python scripts/build_cpp.py
python -m plasma_surface.cli demo --backend cpp --output outputs/demo_cpp
python -m plasma_surface.cli benchmark
python -m pytest -q --basetemp=outputs/pytest-cpp
```

The C++ library implements both exact mean-field integration and continuous-time Gillespie kMC. A permutation partition enables constant-time selection of bare or modified sites. Python's reference kMC scans site arrays, so the measured improvement includes both compilation and an improved selection algorithm. Python still handles rate construction, analysis, plotting, and ML.

Benchmark results are saved to `outputs/benchmark.json`, including workload, repeat timings, event counts, compiler options, and source/binary hashes. Stochastic trajectories differ between backends; agreement is statistical. C++ mean-field outputs are checked against Python to numerical precision. Small mean-field problems may not speed up because binding overhead dominates. Linux/macOS builds are supported by the script but have not been run on this machine. The backend is designed for a source checkout/editable install; a portable binary wheel is not provided.

## Silicon and silicon nitride

```powershell
python -m plasma_surface.cli materials --backend cpp
```

This compares Si, SiN₀.₈, SiN₁.₀, and Si₃N₄ tags using explicit parameter cards in [configs/materials.json](configs/materials.json). Here x=N/Si; stoichiometric Si₃N₄ has x=4/3. Outputs include a comparison figure, CSV, and all effective parameters.

**Composition is metadata in this first version, not an evolving surface state.** Differences between curves come from the chosen hypothetical rates. Hydrogen fraction is recorded but does not drive kinetics. The generic phase sequence is not a validated nitride etch chemistry. The Si thickness conversion is also inherited in these exploratory cards; real SiNₓ calibration must replace it. Preferential N removal, composition enrichment, hydrogen chemistry, film density, and product speciation require a richer reaction network. See the [research plan](docs/RESEARCH_PLAN.md).

## Public data

```powershell
python -m plasma_surface.cli fetch-data
```

The verified public [Si–HCl DFT archive](https://zenodo.org/records/10211009) is ~75 kB and licensed CC BY 4.0. The downloader preserves attribution, source metadata, and checksums; inventory reads the ZIP without extraction. This workspace contains 22 inventoried structures in `data/raw/si_hcl/`. These geometries do not provide a calibrated chlorine/argon ALE model.

For SiNₓ, the [SAIT MLFF benchmark](https://github.com/SAITPublic/MLFF-Framework) is a strong public atomistic starting point. Its Si/N data must not be treated as a reactive halogen force field. Download links, suitability, limitations, and other sources are in [DATASETS.md](docs/DATASETS.md). Large archives are not downloaded automatically.

## Experimental calibration interface

```powershell
python -m plasma_surface.cli calibrate path/to/measurements.csv
```

CSV columns: `energy_ev,dose_s,ion_s,temperature_k,cycle,epc_nm,sigma_nm,source_id`. Supply at least four measurements, positive standard uncertainties, and enough distinct doses/energies to constrain both fitted parameters. The current fitter assumes the exact fluxes, purge durations, initial state, and base silicon parameters in `ale_recipe`; do not feed unrelated chemistries or mixed-material data into it. It fits sticking and chemical yield, reports Jacobian rank/singular values, and records the input checksum. No experimental measurements are bundled or fabricated; the parameter-recovery test uses explicitly synthetic data.

## Etching and deposition extensions

The same kernel supports simultaneous neutral and ion exposure for continuous etching and a generic precursor-arrival channel for net growth:

```python
from plasma_surface.model import Phase, mean_field

continuous_etch = Phase("etch", 10, radical_flux_m2_s=2e19,
                       ion_flux_m2_s=2e19, ion_energy_ev=60)
generic_growth = Phase("growth", 10, precursor_flux_m2_s=1e19)
print(mean_field([continuous_etch], cycles=1)[-1])
print(mean_field([generic_growth], cycles=1)[-1])
```

Growth is a minimal deposition channel, not a PECVD or PEALD mechanism. Positive net removal means etching; negative means growth.

## Simulation and numerical verification campaign

Run `python -m plasma_surface.validation` after building the C++ backend. This executes regression tests, both demos, Si/SiNx scenarios, 4,096 recipes with both backends, stochastic convergence ensembles, half-cycle controls, dose saturation, ML evaluation, and a five-repeat benchmark.

The generated [simulation report](docs/results/REPORT.md) and supporting results are tracked under `docs/results/`. Raw archives and compiled binaries remain local artifacts. The campaign exits with an error if its predefined numerical checks fail. Experimental validation still requires compatible measured data.

For an optional reaction-discovery extension, see the [HiPRGen assessment and proposed connection](docs/HIPRGEN.md). It is an upstream candidate generator, not a calibrated plasma rate source.

## Project map

| Location | Purpose |
|---|---|
| `plasma_surface/model.py` | Reference kinetics, exact solver, Gillespie kMC |
| `cpp/surface.cpp`, `plasma_surface/native.py` | C++ kernel and Python interface |
| `plasma_surface/workflows.py` | Campaign, ML evaluation, two-parameter calibration |
| `plasma_surface/datasets.py` | Public-data download, checksum, geometry inventory |
| `configs/materials.json` | Explicit Si/SiNₓ hypothetical parameter cards |
| `tests/` | Analytical limits, stochastic agreement, calibration, backend checks |
| `hpc/sweep.slurm` | Optional future cluster array template; not required locally |
| `docs/` | Equations, research milestones, datasets, validation guidance |

The SLURM array shards a deterministic design by recipe ID; merging and sorting the shards reproduces serial output. It is a portability template, not evidence of a completed HPC campaign. Source code is currently local and has not been published.
