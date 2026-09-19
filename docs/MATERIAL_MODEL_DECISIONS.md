# Material and mechanism decisions

The current models answer different questions. A finite SiF3NH2/HF proxy isolates terminal-NH2 chemistry, an ideal periodic slab screens adsorption configurations, and the multilayer bond graph tests stochastic bookkeeping and exposure. None is an interchangeable representation of an experimental amorphous film.

| Independent evidence | Implication for this repository | Concrete next calculation or measurement |
|---|---|---|
| [Khumaini et al., 2024](https://doi.org/10.1016/j.apsusc.2024.159414): hydrogenated amorphous nitride, competing hydrogenated volatile products and retained salt | SiF4-only removal and crystalline, hydrogen-poor environments cannot cover the reported chemistry | Match hydrogen inventory and remaining backbonds; compare SiF4, SiHF3 and SiH2F2 pathways on multiple amorphous local environments |
| [Jung et al., 2020](https://doi.org/10.1116/1.5125569): coadsorbates on fluorinated clusters and distinct treatment of vibrationally excited HF | A single HF occupancy flag omits assisted proton transfer; an excited-reactant rate cannot be treated as a ground-state barrier | Compare single HF, HF/HF and HF/H2O paths at matched fluorination, retaining vibrational-state labels |
| [Lill et al., 2024](https://doi.org/10.1116/6.0004019): water effects on low-temperature etching and retained salt | Faster etching need not mean a lower local Si-N cleavage barrier | Separate cleavage, salt formation and salt removal; correlate water dose with retained N/F and volatile products |
| [First-principles OH/O-terminated nitride study, 2012](https://doi.org/10.1021/jp3041605): termination and ring environment affect HF reactions | A bare slab is not a general reference for oxidized nitride; the aqueous context is not a dry-plasma rate source | Build matched O, OH and NHx terminations; use the paper as a mechanistic contrast, not a transferable rate table |
| [Zhang, Zhao and Teo, 2004](https://doi.org/10.1103/PhysRevB.69.125319): fluorination changes silicon backbond strength | A bond-energy trend does not specify an activation barrier | Compute matched fluorination-stage IS/TS/FS; avoid converting bond energies directly into Arrhenius barriers |

These assessments use the accessible primary article abstracts or author manuscript; they do not imply access to unavailable supplementary structures. No new rates are assigned from this table.

## What the present calculations can decide

1. Test whether a candidate is stationary under its stated constraints, then test the complete mobile-coordinate Hessian. The new coadsorbate campaign found a host-coupled instability missed by an adsorbate-only Hessian.
2. Compare energies and forces against direct DFT at identical geometries before using a pretrained model for TS searches. Agreement in one reaction energy is insufficient if endpoint forces remain large.
3. Separate local cleavage from product diffusion/desorption. A band ending in a distant gas fragment can mix processes and produce poorly distributed images; a constrained coordinate scan diagnoses this but is not a TS proof.
4. Keep failed searches as evidence about method reliability, without enabling their peaks in kMC. Then prioritize embedded final-cleavage environments responsible for the observed graph stall.

The next experimental discriminants are hydrogen content/composition, termination-resolved surface spectroscopy, volatile-product ratios, retained salt versus temperature, and EPC versus dose and purge duration. Those observables can distinguish missing chemistry from fitted-rate compensation.
