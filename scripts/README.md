# Script catalog

Run every script from the repository root. The **Env** column says which environment a script needs:

| Env | Meaning | Setup |
|---|---|---|
| core | `pip install -e ".[dev]"` (NumPy, SciPy, Matplotlib) | [usage](../docs/USAGE.md) |
| C++ | core plus a C++17 compiler; run the matching `build_*` script first | [usage](../docs/USAGE.md#c-acceleration) |
| ASE | core plus `pip install -e ".[atomistic]"` | |
| MACE | ASE plus PyTorch and a MACE model checkpoint (not bundled) | [MACE](../docs/MACE.md) |
| OMol25 | ASE plus DeepMD/PyTorch with the OMol25 model | [transition states](../docs/dry_etch_results/TRANSITION_STATES.md) |
| DFT | ASE plus PySCF (the author used WSL) | [transition states](../docs/dry_etch_results/TRANSITION_STATES.md) |

**Ground rules.** *Calculation* scripts produce evidence, and most refuse to overwrite a completed calculation directory. *Report/render* scripts only present saved evidence and never create barriers. Saved results record the path and SHA-256 of the script that produced them, so **do not rename or move these files**. Formatting-only edits are fine if recorded with `record_format_equivalence.py`. Any logic change to a calculation script needs fresh artifacts. Editing a report must never relabel an unconverged path as a validated transition state.

## 1. Build and maintenance

| Script | Env | Purpose |
|---|---|---|
| `build_cpp.py` | C++ | Build the two-state model backend (`build/surface.*`) |
| `build_species_cpp.py` | C++ | Build the species-kMC backend (`build/species.*`) |
| `make_figures.py` | core | Draw the presentation figures in `docs/figures/` from saved results ([guide](../docs/figures/README.md)) |
| `update_research_readme.py` | core | Refresh the generated status blocks in the top-level README |
| `record_format_equivalence.py` | core | Record formatting-only code changes in `provenance/code_format_ledger.json` |

## 2. Species-resolved kMC (45-state HF/SiN:H model)

| Script | Env | Purpose | Output |
|---|---|---|---|
| `run_species_kmc.py` | C++ | Python + C++ ensembles over temperature, verification and sensitivity runs | `docs/species_kmc_results/` |
| `report_species_kmc.py` | core | Report and plots from the saved campaign | `docs/species_kmc_results/REPORT.md` |
| `render_species_lattice.py` | core | Animation built only from recorded snapshots | `species_kmc.gif` |
| `analyze_kinetic_priorities.py` | core | One-at-a-time sensitivities and motif-mixture identifiability | `kinetic_priorities.*` |
| `compare_literature_channels.py` | core | Evaluate published coadsorbate rate laws (no fitting) | `docs/species_kmc_results/` |
| `audit_kinetic_parameters.py` | core | Inventory every kinetic input and reproduce conditional rates | `data/kinetic_audit/`, `docs/KINETIC_PARAMETERS.md` |
| `report_gas_thermochemistry.py` | core | NIST gas chemical potentials | `docs/GAS_THERMOCHEMISTRY.md` |
| `plot_reaction_networks.py` | core | Full enabled/disabled kMC network and HiPRGen ladders | `docs/reaction_network/` |
| `draw_species_network.py` | core | Per-family network panels | `docs/reaction_network/` |
| `draw_surface_states.py` | core | Bookkeeping sketch of every state with exact atom inventory | `docs/reaction_network/surface_states/` |

## 3. Multilayer bond-graph kMC

| Script | Env | Purpose | Output |
|---|---|---|---|
| `build_multilayer_graph.py` | ASE | Build a terminated Si3N4 bond graph from the public CIF | `data/multilayer/sin_graph.json` |
| `run_multilayer.py` | core | Demonstration run, accessibility controls and strict evidence-gated control | `data/multilayer/`, `docs/multilayer_results/summary.json` |
| `analyze_multilayer_stall.py` | core | 12-cycle run explaining why depth progression stops | `docs/multilayer_results/cycle_diagnosis.json` |
| `queue_multilayer_rates.py` | core | List every missing-rate environment as a calculation request | `data/multilayer/rate_requests/` |
| `render_multilayer.py` | core | GIF and storyboard from saved graph states | `docs/multilayer_results/` |
| `report_multilayer.py` | core | Multilayer report and rate-evidence work queue | `docs/multilayer_results/REPORT.md` |

## 4. Reaction-network discovery

| Script | Env | Purpose |
|---|---|---|
| `hiprgen_intermediate_audit.py` | core | Limited HiPRGen bucketing/decision-tree pilot on capped motifs; no rates assigned |
| `check_reaction_catalog.py` | core | Element and charge balance of proposed reactions |

## 5. Public data collection

| Script | Env | Purpose |
|---|---|---|
| `collect_atomistic.py` | ASE | Curate small attributed public DFT subsets |
| `collect_surface_paths.py` | ASE | Fetch published SiCl4/Si(100) stationary points |
| `animate_surface_paths.py` | ASE | Render published IS → TS → FS points (no interpolation) |

## 6. Surface screening and adsorption (MACE)

| Script | Env | Purpose |
|---|---|---|
| `molecular_surface_campaign.py` | MACE | Molecule × surface matrix, keeping every geometry and energy |
| `screen_molecular_sites.py` | MACE | Three sites × distinct orientations per molecule/surface |
| `screen_molecular_sites_energy.py` | MACE | Energy-only variant of the site screen |
| `check_molecular_energies.py` | MACE | Check energy-only evaluation against ASE energy/force evaluation |
| `relax_adsorption_multistart.py` | MACE | Relax HF plus the upper substrate from two screened starts per surface |
| `audit_adsorption_references.py` | MACE | Detect metastable clean-slab references in adsorption energies |
| `report_relaxed_adsorption.py` | ASE | Relaxed-adsorption report (`docs/dry_etch_results/RELAXED_ADSORPTION.md`) |
| `render_molecular_campaign.py` | ASE | Render evaluated geometries only |
| `report_molecular_campaign.py` | core | Summarize saved calculations without reclassifying failed paths |

## 7. Surface paths and saddles (MACE)

| Script | Env | Purpose |
|---|---|---|
| `run_surface_neb.py` | MACE | Rigid-surface halogen diffusion, climbing-image NEB |
| `run_surface_neb_bfgs.py` | MACE | Compatibility entry point for the BFGS runner |
| `refine_surface_neb.py` | MACE | Densify a saved NEB and re-optimize inserted images |
| `validate_surface_neb.py` | MACE | Path energies, endpoint/peak curvature and second-model check |
| `run_molecular_reaction_neb.py` | MACE | Molecular dissociation with a locally mobile surface pair |
| `refine_molecular_surface_paths.py` | MACE | Resume saved molecular NEBs with FIRE |
| `trace_surface_saddle.py` | MACE | Trace both sides of a saddle to the minima actually reached |
| `audit_surface_saddle_identity.py` | ASE | Identify the connected saddle by bond distances |
| `audit_molecular_paths.py` | ASE | Reject geometrically pathological paths before reading peaks |
| `render_calculated_paths.py` | ASE | Render only evaluated NEB images |
| `report_surface_paths.py` | core | Path-report tables from checked outputs |

## 8. Local cleavage proxies (OMol25 / DFT)

| Script | Env | Purpose |
|---|---|---|
| `cluster_reaction_paths.py` | OMol25, DFT | Si-N/Si-O cleavage proxies: NEB plus optional DFT |
| `cluster_reaction_paths_refined.py` | OMol25, DFT | Refined version of the above |
| `dft_molecular_singlepoints.py` | DFT | PBE/def2-SVP energies and forces on saved geometries |
| `resume_dft_singlepoints.py` | DFT | Resume failed single points with a Newton SCF fallback |
| `finish_molecular_dft.py` | DFT | Evaluate completed path files as they appear |
| `check_dft_scf_root.py` | DFT | Alternative SCF starts and stability for a problematic image |

## 9. Follow-up validation campaigns

| Workflow | Calculation scripts | Report |
|---|---|---|
| HF/HF and HF/H2O coadsorbates (MACE) | `relax_coadsorbate_intermediates.py`, `refine_coadsorbate_intermediates.py`, `validate_coadsorbate_minima.py` | `report_coadsorbate_intermediates.py`, `report_followup_validation.py` |
| Final-cleavage molecular proxy (OMol25/DFT) | `final_cleavage_proxy.py`, `refine_final_cleavage_neb.py`, `check_proxy_dft.py`, `check_proxy_force_consistency.py` | `report_followup_validation.py` |
| Fixed-frame local coordinate scan (OMol25) | `scan_final_cleavage_fixed_frame.py`, with helper `fixed_frame_constraints.py` | same report; not a TS calculation |

`scan_final_cleavage_coordinates.py` is kept only to reproduce an invalidated constraint-ordering diagnostic (see its `INVALIDATION.json`). Use the fixed-frame replacement for new work.
