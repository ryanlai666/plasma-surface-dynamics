# Dry-etch reaction parameter results

This run reproduces published silicon/F2 first-reaction kinetics and explores conditional sensitivity to SiN:H/HF barriers. It performs **no new DFT** and does not establish a calibrated plasma ALE etch-per-cycle prediction.

![Published F2 surface rates and first-event kinetics](f2_surface_kinetics.png)

At 500 K and 1 Torr F2, with equal gas and surface temperatures:

| Si facet | Source barrier (eV) | First-reaction hazard (s^-1) |
|---|---:|---:|
| r100 | 0.13 | 25.9988 |
| u100 | 0.31 | 0.555994 |
| 110 | 0.35 | 1.56878 |
| 111 | 0.57 | 0.00493066 |

Source: [Dwivedi et al., Table 1 and equations 5a/5b](https://arxiv.org/abs/2305.09037v2). The first-event Gillespie calculation uses 512 sites and 128 replicates per facet. All ensembles remain within 2.01 analytic standard errors of `1-exp(-r*t)` at the sampled times. Common random numbers are used across facets, so those four numerical checks are correlated, not four independent experiments.

The temperature/facet sweep contains 484 rows in [f2_rates.csv](f2_rates.csv). A separate [SiN:H/HF sensitivity table](sin_hf_sensitivity.csv) contains 192 conditional rates for 16 pathways, four temperatures and three **assumed** prefactors. It does not predict relative product yields across different initial motifs.

[Verification and code/source hashes](verification.json) | [Equations, parameter mapping and applicability](../DRY_ETCH_PARAMETERS.md) | [IS/TS/FS animation gallery](TRANSITION_STATES.md).

Run `python -m plasma_surface.dry_etch` from the repository root. The main regression suite passed 29 tests after this implementation. Published HF-impact yields and salt energetics are reference data, with their distinct units retained; they are not fitted to the generic Si/SiNx cards.
