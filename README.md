# Plasma Surface Dynamics

**Physics-based and machine-learning workflows for silicon and silicon-nitride atomic layer etching.**

Combining surface reaction kinetics, Python/C++ simulation, public atomistic data, and pretrained MACE force-field comparisons. 

**Research question:** How do surface modification, competing removal pathways, and uncertain reaction rates determine etch-per-cycle saturation and the usable ALE energy window?

> **Scientific scope:** This project verifies numerical implementations and compares atomistic models. Default ALE rates are uncalibrated, and SiNx composition cards are hypothetical sensitivity scenarios. The results do not establish experimental etching accuracy.

[Animations](#animation-gallery) | [Quick start](#quick-start) | [Results](#recorded-results) | [DFT/TS](#dft-data-and-transition-state-results) | [Equations and literature](#physical-hypotheses-equations-and-literature-support) | [Documentation](#documentation-and-project-map)

## Animation gallery

Each movie follows three ALE cycles on 256 independent surface columns at 35 eV. **Teal** marks bare sites, **gold** marks modified sites, and **pink** marks removal since the preceding frame. The 3D perspective and the top-view/cross-section movies use the same seeded trajectory. Each 21 s simulation plays in 12.1 s.

### Silicon

![Silicon kMC in 3D perspective](docs/animations/si_3d.gif)

<details>
<summary>Si: top view and linked height cross-section</summary>

![Silicon site map and height cross-section](docs/animations/si.gif)

</details>

### Silicon nitride composition scenarios

<details>
<summary><strong>SiN0.8: 3D perspective and top view</strong></summary>

![SiN0.8 kMC in 3D perspective](docs/animations/sin0p8_3d.gif)

![SiN0.8 site map and height cross-section](docs/animations/sin0p8.gif)

</details>

<details>
<summary><strong>SiN1.0: 3D perspective and top view</strong></summary>

![SiN1.0 kMC in 3D perspective](docs/animations/sin1p0_3d.gif)

![SiN1.0 site map and height cross-section](docs/animations/sin1p0.gif)

</details>

<details>
<summary><strong>Si3N4: 3D perspective and top view</strong></summary>

![Si3N4 kMC in 3D perspective](docs/animations/si3n4_3d.gif)

![Si3N4 site map and height cross-section](docs/animations/si3n4.gif)

</details>

The columns represent model states and removal increments. Their spacing is schematic, and the 3D vertical scale is exaggerated 1.5x; neither view resolves atom identities, chemical bonds, or a crystal lattice. Nitride compositions label hypothetical rate cards.

Reproduce all eight GIFs with `python -m plasma_surface.animate`. [Visualization settings and fidelity checks](docs/ANIMATIONS.md) | [Trajectory and rendering provenance](docs/animations/manifest.json).

## Quick start

Python 3.11+ and a local CPU are sufficient. From a source checkout on Windows:

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install -e ".[dev]"
python -m plasma_surface.cli demo
```

On Linux/macOS, activate with `source .venv/bin/activate`. The demo writes phase histories, a synthetic ML campaign, metrics, plots, and provenance to `outputs/demo/`.

| Task | Command |
|---|---|
| Generate the 2D and 3D animations | `python -m plasma_surface.animate` |
| Compare the four material cards | `python -m plasma_surface.cli materials` |
| Fetch the public DFT geometries | `python -m plasma_surface.cli fetch-data` |
| Run regression tests | `python -m pytest -q --basetemp=outputs/pytest-run` |
| Build the optional C++ backend | `python scripts/build_cpp.py` |
| Run the C++ demo | `python -m plasma_surface.cli demo --backend cpp` |
| Run the full numerical campaign after building C++ | `python -m plasma_surface.validation` |

The C++ build requires a C++17 GCC/Clang compiler. [Detailed setup, benchmarking, calibration, and deposition examples](docs/USAGE.md) | [Optional MACE environment](docs/MACE.md).

## Recorded results

| Check or comparison | Recorded result | Interpretation |
|---|---:|---|
| Regression suite | 25 tests passed | Analytical limits, kinetics, calibration recovery, backend checks, and animation-state fidelity |
| Python/C++ recipe comparison | 4,096 recipes | Same deterministic model evaluated by both implementations |
| Maximum backend EPC difference | 2.22e-16 nm/cycle | Numerical agreement |
| Stochastic ensembles | 6 of 6 within five standard errors | Agreement with the exact mean-field expectation |
| Synthetic surrogate held-out MAE | 0.00284 nm/cycle | Approximation of the uncalibrated model |
| Constant-baseline MAE | 0.05223 nm/cycle | Reference for the synthetic surrogate |
| MACE structural evaluations | 24 across two checkpoints | 12 structures evaluated by each model |
| MACE energy/force consistency | Both checks passed | Finite-difference tolerance: 1e-4 eV/angstrom |
| Public Si-HCl DFT geometry inventory | 22 structures | Provenance and checksums retained |

[Full simulation report](docs/results/REPORT.md) | [Raw benchmark timings](docs/results/benchmark.json) | [Research roadmap](docs/RESEARCH_PLAN.md)

The C++ kMC backend uses constant-time eligible-site selection. The recorded 2,048-site, 10-cycle benchmark had median times of 2.64048 s in Python and 0.002218 s in C++. This workload-specific comparison includes the algorithm change as well as compilation; timings vary with machine load and do not imply a universal speedup.

## Simulation and atomistic comparisons

<details>
<summary><strong>Si ALE: phase histories and ion-energy dependence</strong></summary>

![Silicon ALE phase histories and ion-energy sweep](docs/results/cpp/overview.png)

Exact mean-field integration and stochastic Gillespie simulation describe the same illustrative kinetics. The energy sweep explores modification-dose dependence.

</details>

<details>
<summary><strong>Si/SiNx: comparison of two pretrained MACE models</strong></summary>

![MACE energy-strain comparison](docs/mace_results/comparison.png)

MACE-MP-0b2 small and MACE-MPA-0 medium are evaluated on diamond Si, beta-Si3N4, and two unrelaxed nitrogen-vacancy probes. Energies are referenced separately within each composition. These are structural comparisons, not an etch-rate validation or accuracy ranking.

[Atomistic report](docs/mace_results/REPORT.md) | [Hypothetical composition-rate comparison](docs/results/materials/comparison.png).

</details>

<details>
<summary><strong>Numerical verification: finite-site convergence and dose saturation</strong></summary>

![Finite-site convergence and dose saturation](docs/results/verification.png)

Ensemble scatter decreases as the number of sites increases. The dose sweep checks saturation within the assumed mechanism.

</details>

<!-- BEGIN ATOMISTIC RESULTS -->
## DFT data and transition-state results

A small NH3 inversion calculation tests the molecular saddle workflow using **direct DFT** and the existing **OMol25-trained DPA model**:

| Method | Electronic barrier (eV) | TS max force (eV/A) | Saddle checks |
|---|---:|---:|---|
| DPA-3.3 / OMol25 | 0.196504 | 7.86e-07 | Passed |
| DFT PBE / def2-SVP | 0.277387 | 1.88e-05 | Passed |
| DFT PBE / def2-TZVP | 0.215378 | 3.15e-05 | Passed |

Changing PBE from def2-SVP to def2-TZVP changes this electronic barrier by -0.062009 eV. This two-basis comparison does not establish the complete-basis limit.


![Direct DFT and OMol25 molecular inversion paths](docs/barrier_results/comparison.png)

All reported saddles pass full-force, negative-curvature, and two-sided relaxation checks. NH3 inversion is a molecular workflow diagnostic, **not a SiN surface etching barrier**. The comparison uses different reference methods and excludes zero-point/free-energy corrections.

| Collected reference data | Contents |
|---|---|
| Published Si/O/C/F DFT | 70 energy/force frames, including 25 quasi-static drag configurations; drag peaks are not validated TS barriers |
| Si-H-Cl source geometries | 22 structures, including 15 source-named TS; the archive lacks energy/force/Hessian labels |
| Published Cl diffusion on Si(111)-(5x5) | 1.73 eV hopping and 1.34 eV SiCl-complex diffusion, with sources and applicability limits |
| Fragment/side-reaction catalog | 35 species/bookkeeping units, 21 balanced candidate reactions, and eight proposed experiment sets |

[Computed results and reproduction](docs/barrier_results/REPORT.md) | [Data, diffusion equations, and sources](docs/ATOMISTIC_DATA.md) | [Fragments and experiment sets](docs/REACTION_CANDIDATES.md).

<!-- END ATOMISTIC RESULTS -->

## Physical hypotheses, equations, and literature support

The implemented model resolves bare/modified surface states. With modified fraction $\theta$, modification hazard $a$, desorption $d$, chemical removal $c$, physical removal $s$, and growth $g$ (all in s$^{-1}$):

$$
\frac{d\theta}{dt}=a(1-\theta)-(d+c+s+g)\theta,
\qquad \frac{dD}{dt}=\ell(c\theta+s-g).
$$

Here $D$ is net removed thickness and $\ell$ is the thickness increment. The rates use arrival fluxes, an Arrhenius desorption term, and an illustrative threshold yield:

$$
a=\frac{S\Gamma_r}{n_s},\quad d=\nu\exp[-E_d/(k_BT)],\quad
c=\frac{\Gamma_i}{n_s}Y_c(E),\quad
Y_j(E)=A_j\max\!\left(\sqrt{E/E_j}-1,0\right).
$$

**Hypotheses:** modification saturates available sites; ion-assisted removal depends on modification; newly exposed columns are bare; columns are independent. Increasing nitrogen content is assigned lower sticking and higher thresholds only as a sensitivity hypothesis, not an established composition law.

| Atomistic evidence | Implication for this project |
|---|---|
| [Si chlorination DFT/TDDFT](https://doi.org/10.1016/j.mssp.2022.107169) reports lower removal energy after modification | Supports separate bare/modified states; does not calibrate incident-ion thresholds |
| [SiN + CH3F first-principles study](https://doi.org/10.1016/j.apsusc.2020.148557) finds very low molecular sticking on its H-terminated surface | Current sticking values describe a generic effective radical channel, not measured CH3F sticking |
| [SiN ion-impact AIMD/tight-binding MD](https://doi.org/10.1063/5.0155929) examines adsorbate decomposition | Motivates explicit precursor fragments and collision-assisted reactions in a future mechanism |
| [Amorphous hydrogenated SiN + HF DFT](https://doi.org/10.1016/j.apsusc.2024.159414) resolves bond-cleavage pathways | Thermal barriers need specific states, hydrogen content, and chemistry before transfer to kMC |

**Parameter status:** the current values remain uncalibrated. A DFT reaction barrier in eV is not an ion-impact threshold in eV. The 35 eV animations illustrate the chosen cards; they are not a literature-validated ALE window.

[Detailed evidence, numerical parameters, hypotheses, and DFT/AIMD-to-kMC equations](docs/LITERATURE.md) | [Exact solver equations](docs/MODEL.md) | [Small public DFT barrier extract](data/literature/sin_hf_barriers.csv).

## Modeling workflow

```mermaid
flowchart LR
    A[Public structures and literature] --> B[Atomistic model comparison]
    B -. Future validated rates .-> C[Surface reaction kinetics]
    C --> D[Python and C++ solvers]
    D --> E[Verification and parameter sweeps]
    E --> F[ML surrogate evaluation]
    G[HiPRGen candidate reactions] -. Proposed integration .-> C
```

Solid arrows describe implemented workflows; dashed arrows require further reaction data and model development. [HiPRGen assessment](docs/HIPRGEN.md).

## Public data and scientific scope

The project combines public atomistic references, local molecular DFT/ML calculations, and explicitly synthetic kinetic campaigns. The [Si-HCl DFT archive](https://zenodo.org/records/10211009) contains 22 inventoried geometries with CC BY 4.0 attribution and checksums. Two [public HF/SiN DFT barrier entries](data/literature/sin_hf_barriers.csv) are retained as reference data and are not loaded into the solver. See [dataset sources and suitability](docs/DATASETS.md).

Composition is metadata in the current [Si/SiNx cards](configs/materials.json). Hydrogen content, preferential Si/N removal, film density, product speciation, and lateral interactions do not drive the kinetics. The inherited thickness conversion is illustrative for nitrides. The generic growth channel is not a validated PECVD/PEALD mechanism. Experimental calibration requires measurements matched to the implemented recipe and chemistry; none are fabricated or bundled.

## Documentation and project map

| Topic | Location |
|---|---|
| Equations and exact integration | [Model](docs/MODEL.md) |
| DFT/TS data and computed molecular barriers | [Acquisition](docs/ATOMISTIC_DATA.md), [results](docs/barrier_results/REPORT.md) |
| Fragments, side reactions, and experiment sets | [Candidate mechanisms](docs/REACTION_CANDIDATES.md) |
| DFT/AIMD evidence, parameter audit, and hypotheses | [Literature](docs/LITERATURE.md) |
| Reproduction and extended command examples | [Usage guide](docs/USAGE.md) |
| Animation interpretation and settings | [Animations](docs/ANIMATIONS.md) |
| Validation method and recorded campaign | [Validation](docs/VALIDATION.md), [results](docs/results/REPORT.md) |
| Public atomistic data and MACE comparisons | [Datasets](docs/DATASETS.md), [MACE](docs/MACE.md) |
| Future reaction discovery and research milestones | [HiPRGen](docs/HIPRGEN.md), [research plan](docs/RESEARCH_PLAN.md) |
| Python kinetics and animation rendering | `plasma_surface/model.py`, `plasma_surface/animate.py` |
| C++ kernel and Python interface | `cpp/surface.cpp`, `plasma_surface/native.py` |
| Campaigns, ML, and calibration | `plasma_surface/workflows.py` |
| Material assumptions and regression checks | `configs/materials.json`, `tests/` |
| Optional cluster array template | `hpc/sweep.slurm` |

The SLURM array shards a deterministic design by recipe ID; merging and sorting the shards reproduces serial output.
