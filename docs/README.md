# Research and repository guide

Start with the [project overview](../README.md), then use the relevant model and evidence track below. These models have different scopes; their rates and validation claims are not interchangeable.

## Models and results

| Model | Implementation | Results | What it represents |
|---|---|---|---|
| Species-resolved HF/SiN:H | [Python](../plasma_surface/species_kmc.py), [C++](../cpp/species.cpp) | [45-state results](species_kmc_results/REPORT.md) | Conditional motif kinetics; one initial inventory |
| Explicit multilayer graph | [Python graph engine](../plasma_surface/multilayer.py) | [Multilayer report](multilayer_results/REPORT.md) | Bonds, terminations and accessibility across six depth bands; unvalidated rates |
| F2/Si first event | [Source rate law](../plasma_surface/dry_etch.py) | [Dry-etch results](dry_etch_results/REPORT.md) | Published facet-specific first-reaction kinetics |
| Generic two-state baseline | [Model](../plasma_surface/model.py) | [Model equations](MODEL.md) | Software/sensitivity example with illustrative parameters |

## Evidence to kinetics

1. [Public structures and datasets](DATASETS.md), [atomistic inventory](ATOMISTIC_DATA.md), [literature](LITERATURE.md).
2. [Molecular/surface calculations](dry_etch_results/MOLECULAR_SURFACES.md), [relaxed adsorption](dry_etch_results/RELAXED_ADSORPTION.md), [transition states](dry_etch_results/TRANSITION_STATES.md).
3. [Parameter inventory and thermal-rate equations](KINETIC_PARAMETERS.md), [source energetic mapping](DRY_ETCH_PARAMETERS.md), [rate qualification gate](../plasma_surface/rate_evidence.py).
4. [Vertical kMC algorithm](../README.md#kmc-sampling-flowchart), [sampling details](KMC_ALGORITHM.md), [validation](VALIDATION.md).

## Unresolved chemistry and next calculations

- [New coadsorbate structures and calculations](intermediate_results/REPORT.md).
- [Intermediate gap register](../data/reaction_network/intermediate_gaps.csv) and [HiPRGen audit](HIPRGEN.md).
- [Cycle-depth stall diagnosis](multilayer_results/cycle_diagnosis.json): 13 bare-N and eight NH2 final-cleavage sites block continued recession.
- [Environment-specific calculation queue](../data/multilayer/rate_requests/index.json).
- [Scientific critique](SCIENTIFIC_REVIEW.md), [fragment and side-reaction hypotheses](REACTION_CANDIDATES.md), [research plan](RESEARCH_PLAN.md).

## New assumption checks

- [Full coadsorbate stability](intermediate_results/FULL_STABILITY.md): mobile-host Hessians and unstable-reference exclusions.
- [Final-cleavage DFT comparison](final_cleavage_results/REPORT.md): two basis sets, failed NEBs and a checked constrained scan.
- [Gas thermochemistry](GAS_THERMOCHEMISTRY.md): NIST coefficients, pressure corrections and reference limits.
- [Material/model decisions](MATERIAL_MODEL_DECISIONS.md): independent literature and discriminating tests.

## Folder map

| Folder | Purpose |
|---|---|
| `plasma_surface/` | Importable simulation, rate and analysis code |
| `cpp/` | Compiled backends for the original models; not the new graph engine |
| `configs/` | Model topology, assumptions and qualified-rate library |
| `data/literature/`, `data/reference/` | Source extracts and public structures |
| `data/surface_paths/` | Per-system evaluated atomistic structures and energies |
| `data/reaction_network/` | Candidate-network evidence and intermediate gaps |
| `data/multilayer/` | Explicit graph, trajectories and rate requests |
| `data/intermediate_campaign/` | Per-surface coadsorbate calculations, references and candidate connectivity |
| `data/final_cleavage/` | Molecular-proxy paths, DFT diagnostics and explicit invalidation records |
| `data/kinetic_audit/` | Units, provenance, parameters and evaluated conditional rates |
| `docs/*_results/`, `docs/reaction_network/` | Reports, plots, animations and rendering provenance |
| `scripts/` | [Reproducible calculation/report workflows](../scripts/README.md) |
| `tests/` | Conservation, numerical, backend and provenance checks |
| `third_party/` | Licensed pinned source snapshots |
| `archive/` | Historical/unrelated diagnostics, excluded from active evidence |

Existing source paths are retained because calculation records and hashes reference them. This index organizes their purpose without breaking provenance. [Earlier extended project notes](PROJECT_DETAILS.md) remain available as historical context.
