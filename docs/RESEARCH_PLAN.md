# Research portfolio plan: Si and SiNₓ plasma-surface modeling

## Objective and scope

Build a defensible multiscale study of how surface modification and ion-driven removal compete in silicon ALE, then investigate which additional states are required for non-stoichiometric silicon nitride. Work on a local CPU and reuse public data. A strong outcome is a reproducible analysis that explains where a model succeeds, where it fails, and which new measurements or calculations would resolve the failure.

The working software already supplies deterministic and stochastic solvers, an accelerated C++ implementation, surrogate baselines, data provenance, and tests. New DFT/MD calculations, reactive force-field training, and experimental validation remain research work; they have not been performed by creating this repository.

## Suggested milestones

| Stage | Work | Reviewable deliverable / acceptance criterion |
|---|---|---|
| 1: weeks 1–2 | Reproduce local examples, map phase histories, vary seeds/site count, inspect public structures | Equations, analytical limits, stochastic convergence plots, provenance record |
| 2: weeks 3–4 | Curate a small compatible published Si dataset; record digitization uncertainty if only figures exist | Source-linked table, units, conditions, independent validation split; no invented labels |
| 3: weeks 5–6 | Fit low-dimensional reaction model; run chemical-only/ion-only controls and dose sweeps | Identifiability study, residual plots, credible domain limits, parameter intervals |
| 4: weeks 7–8 | Introduce SiNₓ structure data and composition descriptors; hold out a composition | Composition-resolved atomistic benchmark; distinguish structural from reaction accuracy |
| 5: weeks 9–10 | Replace hypothetical nitride cards with a species-resolved reaction network where data permit | Surface N/Si evolution and species balance; calibrated thickness conversion |
| 6: weeks 11–12 | Compare direct simulation and surrogates; publish a reproducible research report locally | Accuracy/cost comparison, ablations, failure cases, documented CPU runtime |

Timing is an organizational estimate, not a guarantee of publishable new physics. If compatible reaction data cannot be obtained, finish an explicit model-identifiability and public-data-gap study rather than claim calibrated process predictions.

## Atomistic-to-mesoscale bridge

1. Inspect public relaxed and transition-state geometries. Preserve cell, orientation, termination, composition, charge, and calculation method. Geometry alone does not determine a reaction barrier.
2. For a future DFT study, reproduce reference energies with converged slab thickness, vacuum, k-points, cutoff, and electronic treatment. Establish adsorption energies and transition states, checking the transition-state mode and reference chemical potentials. Use an accessible code only after confirming resources and pseudopotential compatibility.
3. For MD, validate a potential against the actual Si/N/H/F/Cl/Ar species and collision energies needed. A Si–N model that reproduces elastic properties is not automatically suitable for halogen etching. Report energy/force errors by environment, trajectory stability, and collision outcomes.
4. Estimate removal probabilities and products across ion energy, angle, dose, coverage, temperature, and composition. Use repeated independent impacts; report binomial/count uncertainties, not only means. Keep trajectories from the same slab/impact family in the same split.
5. Translate those statistics into reaction-specific hazards, then test the mesoscale model against held-out MD or published measurements. Preserve the link from each rate to its source, fitting range, and uncertainty.

These steps are a research protocol, not an implemented DFT/MD driver. The current public archive supports step 1. The SAIT SiN dataset supports structure/energy/force benchmarking but requires additional reactive configurations for etching.

## SiNₓ questions

Study x=N/Si at 0.8, 1.0, and 4/3, but treat hydrogen content, density, deposition history, and surface termination as separate covariates. Do not assume composition alone determines etch response. Compare at fixed incident conditions and separately at fixed target EPC.

A useful next model would track exposed Si and N sites, their chemical modification/passivation states, and a finite subsurface reservoir. Preferential removal changes the exposed composition; replenishment draws from the bulk composition. Every event must conserve Si/N counts, account for outgoing products, and use material-specific volume conversion. Test zero-yield limits and species balances before interpreting selectivity. The current composition cards deliberately do not claim these capabilities.

## Physics + ML study design

The implemented ExtraTrees model predicts synthetic ROM EPC. It uses recipe-group holdouts, a constant baseline, and a separate high-energy holdout audit. This is a workflow benchmark, not an ML interatomic potential or experimental virtual-metrology model. Tree spread is only a heuristic diagnostic.

For a research result, compare the physical model, a data-only surrogate, and a residual correction to the physical model. Reserve complete process runs/material compositions for final testing. Fit scaling and hyperparameters only within training folds. Audit uncertainty coverage on held-out data and separate stochastic noise from uncertain physics. Active sampling should seek both model disagreement and meaningful physical regimes; avoid choosing points exclusively where an uncalibrated surrogate is confident.

For public optical-emission data, start with run-level grouping and simple spectral baselines. Do not infer etch rates without matching thickness measurements, or pool optical spectra and DFT energies as if they were the same training target.

## Experimental collaboration artifact

Even without laboratory access, prepare a measurement request table: material N/Si/H and density, chamber/wafer ID, gas and wall-conditioning history, delivered ion energy distribution, flux, pulse timing, temperature, thickness versus cycle, and uncertainty. Request single-half-cycle controls, saturation sweeps, replicates, and an untouched validation campaign. Interpret OES with potential wall contributions and instrument response in mind.

## Portfolio claims supported today

- Implemented and tested mean-field and stochastic surface kinetics in Python and C++.
- Benchmarked local CPU performance with reproducible workload and build metadata.
- Curated and inventoried a real public DFT geometry dataset with attribution and checksums.
- Built synthetic surrogate evaluation and a tested calibration interface.
- Defined SiNₓ sensitivity scenarios and a public-data path to an eventual composition-resolved model.

Do not claim completed ab initio reaction calculations, a validated reactive ML force field, measured SiN selectivity, experimental agreement, or production HPC deployment until that work exists.
