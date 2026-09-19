# MACE comparison plan and execution

## Selected checkpoints

The [official MACE model registry](https://github.com/ACEsuit/mace-foundations) lists MACE-MP-0b2 and MACE-MPA-0 as materials models covering Si and N, with MIT-licensed checkpoints. This project compares the small MP-0b2 model with medium MPA-0. The pair offers a practical CPU baseline and a larger alternative with different training coverage. The size and training-set differences are confounded; this is not a controlled architecture ablation.

Other candidates include MACE-MATPES-PBE/r2SCAN and MACE-MH-1. Their model cards and licenses should be reviewed for the intended use before adding them. MACE-OFF23 is unsuitable for this Si-based test because its published element set excludes Si. The [foundation-model documentation](https://mace-docs.readthedocs.io/en/latest/guide/foundation_models.html) describes these model families and ASE usage.

## Reproduce locally

Use a separate environment so MACE dependencies do not change the original simulator environment:

```powershell
python -m venv .venv-mace
.\.venv-mace\Scripts\python.exe -m pip install -r requirements-mace.txt
.\.venv-mace\Scripts\python.exe -m plasma_surface.mace_compare
```

Use an environment with its own numerical libraries. Sharing Anaconda MKL packages with the PyTorch wheel caused an OpenMP runtime collision during the first attempt; a fully isolated environment avoids that combination. The result manifest records actual versions. The script downloads official checkpoint files into `.cache/mace/`, records hashes, and evaluates on CPU in float64 with two Torch threads. The models and environment are excluded from Git.

Results are written to [mace_results/REPORT.md](mace_results/REPORT.md). Each run compares 12 structures across two checkpoints, followed by a central finite-difference force check for each model. Outputs include model disagreement, shifted energy-strain curves, all input structures, model hashes, and environment metadata.

## Interpretation

Public beta-Si3N4 coordinates come from [COD 2102550](https://www.crystallography.net/cod/2102550.html), CC0, attributed to du Boulay et al. (2004). Diamond Si is generated at a chosen 5.43 angstrom lattice parameter. Si6N6 and Si6N5 are derived by removing selected nitrogen atoms from the nitride cell. They are unrelaxed defect probes, not realistic equilibrated amorphous films at those compositions.

Compare forces and energy changes within fixed composition. Raw total energies across different atom counts do not establish relative stability; proper formation energies need consistent chemical potentials. Neither checkpoint supplies an independent reference for the other. Gradient checks establish numerical consistency only, and no energy/force accuracy MAE is reported without reference labels.

A follow-up accuracy benchmark should use the public SAIT SiN dataset with its held-out trajectories/compositions and verify DFT energy conventions. Plasma etching needs additional surface, halogen, collision, and reaction-path data. Do not replace the ALE model rates with bulk MACE energies or turn reaction free energies directly into activation barriers.

## Completed local comparison

Both checkpoints ran successfully on CPU: 12 structures each, 24 energy/force evaluations total. Both gradient checks passed, with absolute differences of 3.42e-8 and 6.83e-8 eV/angstrom at the tested displaced coordinate.

At zero strain, the RMS force-component disagreement is 0.0474 eV/angstrom for beta-Si3N4, 0.2640 for the Si6N6 vacancy probe, and 0.3101 for Si6N5. The largest disagreement in these probes is 0.4910 eV/angstrom for compressed Si6N5. This suggests prioritizing defect environments for reference calculations; it does not establish which model is correct. Nearly zero forces in uniformly strained perfect diamond Si follow from symmetry and do not demonstrate force accuracy for distorted silicon environments.

The ASE CIF reader emitted a generic hexagonal-setting warning. The loaded cell contains 6 Si and 8 N atoms and preserves the source cell dimensions; the exact expanded structures are saved for inspection. No relaxation or atomistic etching trajectory was performed.

## Surface reaction paths

A separate [surface-path workflow](dry_etch_results/SURFACE_PATHS.md) now computes rigid-surface F/Cl migration with MACE CI-NEB, checks adsorbate-only saddle curvature, and evaluates a second model on the same images. These calculations include relaxation and reaction-path energies; they do not change the interpretation of the earlier bulk-strain comparison.
