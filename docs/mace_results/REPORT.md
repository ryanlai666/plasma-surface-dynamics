# MACE force-field comparison

Evaluated 12 structural probes using two official checkpoints on CPU (24 energy/force evaluations).

- MACE-MP-0b2 small and MACE-MPA-0 medium.
- Diamond Si and public beta-Si3N4 (COD 2102550).
- Derived Si6N6 and Si6N5 cells with unrelaxed nitrogen vacancies.
- Each structure at -3%, 0%, and +3% isotropic linear strain.
- Energy-gradient finite-difference checks passed: True (absolute tolerance 1e-4 eV/angstrom).

![Comparison](comparison.png)

[Predictions](predictions.csv), [force disagreement](model_disagreement.csv), [gradient checks](force_gradient_checks.json), and [checkpoint/environment provenance](manifest.json).

Energy offsets are removed separately for each structure and model. Do not use these plots to rank compositions thermodynamically. No reference DFT forces or energies are available for these probes, so model disagreement is not an accuracy score. The finite-difference check verifies energy/force consistency, not agreement with experiments. Vacancy cells are deliberately unrelaxed structural tests, not synthesized or equilibrated SiNx films. This comparison does not validate plasma-impact chemistry or ALE rates.

The public crystal data are CC0: du Boulay et al., Acta Crystallographica B 60 (2004), 388-405, [COD 2102550](https://www.crystallography.net/cod/2102550.html). Models and citations: [MACE foundation-model repository](https://github.com/ACEsuit/mace-foundations). Model checkpoints remain cached locally and are not committed.
