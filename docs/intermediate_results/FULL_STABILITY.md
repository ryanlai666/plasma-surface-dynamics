# Full mobile-coordinate stability checks

All four single-HF references and eight coadsorbates were checked using every mobile host and adsorbate coordinate. Lower-host atoms remain fixed. Central differences use 0.01 angstrom. Negative-mode displacements and independent BFGS escapes are retained, including failures.

| Surface / added molecule / start | Parent minimum curvature (eV/A2) | Candidate minimum curvature (eV/A2) | Candidate force (eV/A) | Raw energy difference (eV) | Screen |
|---|---:|---:|---:|---:|---|
| beta_Si3N4_001 / HF / 1 | -0.3540 | +0.0031 | 0.0164 | -1.556 | exclude |
| beta_Si3N4_001 / H2O / 1 | -0.3540 | +0.0445 | 0.0162 | -2.444 | exclude |
| beta_Si3N4_001 / HF / 2 | +0.3693 | +0.3497 | 0.0189 | -3.742 | pass |
| beta_Si3N4_001 / H2O / 2 | +0.3693 | +0.3646 | 0.0189 | -1.025 | pass |
| alpha_quartz_001 / HF / 1 | +0.0449 | +0.0187 | 0.0099 | -0.826 | pass |
| alpha_quartz_001 / H2O / 1 | +0.0449 | +0.0065 | 0.0083 | -0.754 | pass |
| alpha_quartz_001 / HF / 2 | -0.0157 | +0.0386 | 0.0184 | -0.678 | pass |
| alpha_quartz_001 / H2O / 2 | -0.0157 | -0.0006 | 0.0198 | -0.773 | pass |

A pass requires both parent and candidate forces <=0.02 eV/A and no eigenvalue below -0.02 eV/A2. Small negative modes are not automatically physical zero modes. This is provisional MACE screening, not proof of a DFT minimum. Unweighted Hessian curvatures are not vibrational frequencies.

The first nitride parent has appreciable negative curvature despite a small force. Its two coadsorbate energy differences retain an unstable reference and are **excluded from screened association energies**. The JSON stores null for excluded values while preserving raw diagnostics. Force convergence alone is insufficient.

The oxide/water start-1 negative mode was explicitly followed. Compare the original and follow-up structures to inspect changes under identical constraints. No activation barrier or kMC rate is inferred.

[Full Hessians and escape metadata](../../data/intermediate_campaign/full_stability.json) | [Eligibility records](../../data/intermediate_campaign/stability_qualification.json). Each system has a separate `refinement/full_stability` folder. Reproduce the calculator using `scripts/validate_coadsorbate_minima.py` in a fresh output location, then this reporter.
