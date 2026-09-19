# Relaxed HF adsorption: four surfaces, two distinct starts each

This is a force-driven follow-up to rigid approach screening. HF and the **upper half of the substrate atoms by height** can move; atoms at or below the median substrate height remain fixed. The bare slab uses exactly the same fixed indices. Cell vectors stay fixed. MACE-MP-0b2 small and FIRE are used with a 0.04 eV/angstrom force criterion and 300-step budget.

Starts are selected from evaluated scan heights >=1.9 angstrom: first the lowest interaction-energy point, then the lowest point with both a different site and a different orientation. Thus two different starts are actually tested, but this does not exhaust orientations, reconstructions or amorphous environments.

![Energy during relaxation](relaxed_adsorption.png)

| Surface / start | Initial site / orientation | Combined / bare slab converged | Final mobile fmax (eV/A) | Relaxation energy (eV) | Apparent adsorption difference (eV) | HF distance (A) | Artifacts |
|---|---|---|---:|---:|---:|---:|---|
| Si100 / 1 | hollow / upright | yes / yes | 0.0390 | -3.041 | -3.678 | 0.952 | [data + GIF](../../data/surface_paths/Si100/HF/relaxed_adsorption/start_1/README.md) |
| Si100 / 2 | bridge_Si_Si / parallel | yes / yes | 0.0384 | -0.075 | -0.440 | 0.954 | [data + GIF](../../data/surface_paths/Si100/HF/relaxed_adsorption/start_2/README.md) |
| Si111 / 1 | hollow / parallel | yes / yes | 0.0361 | -9.108 | -10.193 | 0.951 | [data + GIF](../../data/surface_paths/Si111/HF/relaxed_adsorption/start_1/README.md) |
| Si111 / 2 | bridge_Si_Si / upright | no / yes | 0.1398 | -8.025 | not reported | 0.947 | [data + GIF](../../data/surface_paths/Si111/HF/relaxed_adsorption/start_2/README.md) |
| beta_Si3N4_001 / 1 | atop_N / upright | yes / yes | 0.0399 | -1.415 | -0.246 | 0.949 | [data + GIF](../../data/surface_paths/beta_Si3N4_001/HF/relaxed_adsorption/start_1/README.md) |
| beta_Si3N4_001 / 2 | bridge_Si_N / parallel | yes / yes | 0.0380 | -1.545 | -0.303 | 0.946 | [data + GIF](../../data/surface_paths/beta_Si3N4_001/HF/relaxed_adsorption/start_2/README.md) |
| alpha_quartz_001 / 1 | atop_O / flipped | yes / yes | 0.0390 | -1.123 | -0.855 | 1.007 | [data + GIF](../../data/surface_paths/alpha_quartz_001/HF/relaxed_adsorption/start_1/README.md) |
| alpha_quartz_001 / 2 | bridge_Si_O / parallel | yes / yes | 0.0386 | -1.487 | -0.897 | 1.026 | [data + GIF](../../data/surface_paths/alpha_quartz_001/HF/relaxed_adsorption/start_2/README.md) |

`Relaxation energy = E_final - E_start` follows one composition and one starting geometry. It can include substrate reconstruction and does not measure an activation barrier. FIRE need not decrease the energy monotonically at every step.

`Apparent adsorption difference = E_relaxed(slab+HF) - E_relaxed(slab) - E_isolated(HF)` is tabulated only if the combined and bare-slab optimizations both meet the force criterion. The gas geometry is the previously optimized isolated HF using the same model. The combined system has finite periodic coverage. A dissociated HF geometry still uses intact HF as the reference. Reconstruction into different local basins can influence this difference; it is not a converged experimental adsorption enthalpy.

HF distance >1.4 angstrom is only a screening label for bond separation. See the reference decomposition below before interpreting this number as HF binding. Stable minima require Hessian checks; elementary reactions need connected saddle searches. Force convergence alone supplies neither. No rates were changed using these results.

Remaining requirements: realistic termination and amorphous SiNx:H structures, coverage and cell/thickness convergence, more starts, curvature and connectivity checks, then matched DFT and experimental validation. The calculations retain ideal unpassivated surfaces and the MACE domain limitations described in the [molecular campaign](MOLECULAR_SURFACES.md).

Reproduce in the optional MACE environment:

```sh
python scripts/relax_adsorption_multistart.py
python scripts/report_relaxed_adsorption.py
python scripts/update_research_readme.py
```

Outputs are per surface/species under `data/surface_paths/<surface>/HF/relaxed_adsorption/`. Every GIF uses sampled, evaluated optimization steps; none uses interpolated geometries.

## Reference-state audit: separate binding from reconstruction

The apparent adsorption difference decomposes exactly as `E_interaction(frozen fragments) + [E_slab(final geometry) - E_slab(reference)] + [E_HF(periodic final geometry) - E_HF(isolated reference)]`. The molecular term includes both distortion and periodic-layer reference effects. A negative slab term demonstrates that the original clean reference was not the lowest accessible geometry, even if its force criterion passed.

| Surface / start | Frozen-fragment interaction (eV) | Slab-reference shift (eV) | Molecular reference shift (eV) | Lower clean-slab geometry found? |
|---|---:|---:|---:|---|
| alpha_quartz_001 / 1 | -1.000 | 0.098 | 0.046 | not flagged |
| alpha_quartz_001 / 2 | -0.973 | -0.012 | 0.088 | not flagged |
| beta_Si3N4_001 / 1 | -0.261 | 0.011 | 0.004 | not flagged |
| beta_Si3N4_001 / 2 | -0.316 | 0.007 | 0.006 | not flagged |
| Si100 / 1 | -0.627 | -3.054 | 0.002 | yes |
| Si100 / 2 | -0.440 | -0.001 | 0.001 | not flagged |
| Si111 / 1 | -0.456 | -9.740 | 0.003 | yes |
| Si111 / 2 | -0.630 | -8.386 | 0.005 | yes |

The -0.1 eV flag is a diagnostic threshold, not a validated accuracy tolerance. Flagged apparent differences are excluded from quantitative adsorption claims and kinetic parameterization. Remedy: independently reconstruct/passivate bare slabs, perform multistart reference searches, then compare identical structural basins and DFT forces. Merely allowing more atoms to move does not solve reference-state bias.

[All numerical decompositions](adsorption_reference_audit.json). Run `python scripts/audit_adsorption_references.py` before regenerating this report.

Calculation metadata retain the original Windows relative-path strings for provenance; reporting and audit readers normalize their separators.
