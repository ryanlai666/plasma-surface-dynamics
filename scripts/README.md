# Calculation and reporting workflows

Run scripts from the repository root, using the appropriate installed environment. The [usage guide](../docs/USAGE.md) covers setup. Calculation scripts produce evidence; report/render scripts present saved evidence and do not create new barriers.

| Task | Run or inspect | Principal output |
|---|---|---|
| Audit every active kinetic input | `audit_kinetic_parameters.py` | `data/kinetic_audit/`, `docs/KINETIC_PARAMETERS.md` |
| Reproduce species kMC | `run_species_kmc.py`, `report_species_kmc.py`, `render_species_lattice.py` | `docs/species_kmc_results/` |
| Build/run multilayer graph | `build_multilayer_graph.py`, `run_multilayer.py` | `data/multilayer/`, `docs/multilayer_results/summary.json` |
| Diagnose stalled depth progression | `analyze_multilayer_stall.py` | 12-cycle trajectory and cycle-depth diagnosis |
| Generate missing-rate requests | `queue_multilayer_rates.py` | `data/multilayer/rate_requests/` |
| Present multilayer results | `render_multilayer.py`, `report_multilayer.py` | GIF, depth controls, report |
| Enumerate reaction candidates | `hiprgen_intermediate_audit.py` | Limited HiPRGen pilot; no automatic rates |
| Draw species chemistry | `draw_surface_states.py`, `draw_species_network.py` | `docs/reaction_network/` |
| Screen and relax surfaces | `molecular_surface_campaign.py`, `relax_adsorption_multistart.py` | Per-surface/per-species MACE calculations |
| Search/check local paths | `run_molecular_reaction_neb.py`, `cluster_reaction_paths_refined.py`, `validate_surface_neb.py` | Evaluated paths with explicit convergence status |
| Direct molecular DFT diagnostics | `dft_molecular_singlepoints.py`, `check_dft_scf_root.py` | Energies/forces and SCF checks; not automatically DFT-optimized paths |
| Refresh generated overview blocks | `update_research_readme.py` | Status sections only |

Review each script's arguments and output policy before rerunning expensive work. Many calculations intentionally preserve existing outputs. Changes to calculation code require fresh artifacts or explicit provenance; editing a report must not relabel an unconverged path as a validated transition state.
