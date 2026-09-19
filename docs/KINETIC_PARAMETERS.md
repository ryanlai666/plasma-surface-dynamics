# Kinetic parameters, free energies and evidence

**Current implementation:** published DFT barriers with assumed Arrhenius prefactors drive the conditional HF models. A separate F2 first-event model uses a published cluster-DFT/TST fit. No AIMD-derived rates, full activation free energies or diffusion hops are enabled in the active kMC engines. The strict multilayer qualified-rate library is empty.

The reproducible audit contains **183 parameter records**, **220 event/temperature rates** (55 channels at four temperatures), **18 multilayer thermal rules**, and **16 source-specific coadsorbate rate evaluations**. Numerical values and missing entries are distinguished explicitly; a missing free energy is not zero.

[All parameter values, units, status and source](../data/kinetic_audit/parameters.csv) | [Every motif event rate](../data/kinetic_audit/species_event_rates.csv) | [Multilayer rules](../data/kinetic_audit/multilayer_rate_rules.csv) | [Input hashes](../data/kinetic_audit/manifest.json).

## Active controls and assumptions

| Quantity | Motif model | Multilayer demonstration | Evidence |
|---|---|---|---|
| Temperature | 325, 350, 400, 450 K | 450 K | Chosen simulation controls |
| HF arrival | 5 s^-1 per available motif | 5 exp(-depth / 1 A) s^-1 in dose | Assumed; not derived from measured flux |
| HF complex desorption | 100 s^-1 | 100 s^-1 | Assumed |
| Thermal prefactor | 1e12 s^-1 | 1e12 s^-1 | Assumed; no event-specific vibrational entropy |
| Dose / purge | 2 s / 1 s | 2 s / 1 s repeated | Chosen protocol |
| Access depth / attenuation | Not spatial | 2 A / 1 A | Assumed geometric rule |
| Detached Si-product release Ea | Literature motif event where specified | 0.53 eV | Multilayer desorption value is assumed, not source P4 |
| Detached NH3 release Ea | Included in source motif channel | 0.45 eV | Multilayer desorption assumed |
| Si-Cl exchange Ea | Not included | 0.75 eV | Assumed |
| Final SiF3-NH2 / SiF3-N cleavage | Disabled | Disabled | Missing matched final-step barriers |
| Lateral diffusion | Absent | Absent | No calibrated hop rates |
| Activation free energy | Not supplied | Not supplied | No matched full IS/TS thermochemistry |

The audit also lists all ten generic two-state model defaults, recipe fluxes/energies, composition-card overrides, initial motif fractions, graph cap probabilities, F2 law parameters and salt-reference energies. The generic model's 15/40 eV ion thresholds are not thermal activation barriers; its 0.65 eV modifier loss barrier is illustrative.

## Source SiN:H/HF energies

| Path | Ea (eV) | Ephy (eV) | Reaction DeltaE (eV) |
|---|---:|---:|---:|
| P1a | 0.45 | -1.09 | -0.24 |
| P1b | 2.10 | -0.83 | -0.63 |
| P1c | 0.02 | -0.91 | -1.40 |
| P1d | 0.72 | -0.90 | -0.97 |
| P2a | 0.68 | -0.54 | -0.76 |
| P2b | 1.80 | -0.42 | -0.65 |
| P2c | 0.79 | -0.35 | -0.24 |
| P2d | 0.75 | -0.51 | -0.96 |
| P2e | 0.81 | -0.51 | -1.29 |
| P2f | 0.79 | -0.25 | -1.31 |
| P3a | 0.88 | -0.64 | -0.42 |
| P3b | 1.54 | -0.29 | -1.08 |
| P3c | 0.74 | -0.37 | -0.67 |
| P3d | 0.90 | -0.30 | -0.56 |
| P3e | 0.96 | -0.39 | -1.13 |
| P4 | 0.53 | -0.32 | -0.77 |

[Source extraction and reference caveats](DRY_ETCH_PARAMETERS.md#sinhhf-all-fluorination-paths-and-competing-fragments). These are different source environments, not a common branching state. `Ea`, physisorption energy and reaction energy have separate meanings. In the current conditional kinetics, only `Ea` enters `nu exp(-Ea/kBT)`; Ephy and DeltaE do not silently become rate barriers.

## Thermal-rate workflow

The new [thermal kinetics module](../plasma_surface/thermal_kinetics.py) implements the following calculations with explicit units. It produces numerical candidates and does not automatically enable a kMC event.

1. **Arrhenius:** `k = nu exp(-Ea / kBT)`, with `Ea` in eV, `nu` in s^-1 and T in K. A 0.1 eV barrier change is exponentially significant; uncertainty should propagate through the rate rather than disappear into a fitted prefactor.
2. **Harmonic adsorbate thermochemistry:** `F_vib = sum[epsilon/2 + kBT ln(1-exp(-epsilon/kBT))]`, using positive mode energies `epsilon = h nu`. The first term is ZPE. The input IS and TS must represent the same environment and constrained degrees of freedom. Exactly one TS unstable mode is excluded; zero/imaginary stable-mode inputs are rejected. Soft translations or hindered rotations require an explicit physical treatment, not arbitrary deletion. [ASE thermochemistry implementation](https://docs.ase-lib.org/_modules/ase/thermochemistry.html).
3. **Surface HTST:** `DeltaF_dagger = (E_TS-E_IS) + F_vib_TS - F_vib_IS`, then `k = kappa (kBT/h) exp(-DeltaF_dagger/kBT)`. This is a fixed-volume surface Helmholtz approximation. It is not automatically a Gibbs activation free energy. A Gibbs treatment also needs consistent standard states, pressure/chemical potentials, and applicable gas translations/rotations. Do not multiply an independent Arrhenius prefactor onto the Eyring result.
4. **Gas arrival:** ideal-gas flux `J = p/sqrt(2 pi m kBT_g)` and per-site hazard `a_ads = J S A_site`. Inputs are pressure (Pa), gas temperature, molecular mass, site area and sticking. Surface temperature belongs in thermal surface hazards; it need not equal gas temperature.
5. **Surface diffusion:** `k_hop = nu_hop exp(-E_m/kBT)` and `D = f z k_hop l^2/(2d)` for equivalent isotropic neighbor hops, with per-neighbor rate, hop length l, coordination z, dimension d and correlation factor f. Current Cl/Si(111) reference barriers are 1.73 eV for Cl and 1.34 eV for SiCl; missing prefactors and differing environments prevent a numerical current-model diffusion rate. [Recorded sources](../data/literature/diffusion_barriers.csv).

The module requires actual matched mode data for HTST. The coadsorbate campaign's partial, unweighted adsorbate Hessian is not a complete mass-weighted IS/TS spectrum and is therefore not used as thermochemical input.

## AIMD and thermal kinetics: what is still needed

AIMD can test precursor stability, proton-transfer events and temperature-dependent motion. A short trajectory with no reaction does not measure a zero rate or an activation barrier. Define reactant/product basins, exclude recrossings and accumulate eligible-site residence time. Under a stationary Poisson model, `k_hat = N_committed_events / total_reactant_site_time`. The helper `md_count_rate` supplies the estimate and a one-sided Poisson upper bound; these assumptions require checking with replicas and waiting-time statistics.

For illustration only, zero events during 100 ps at one continuously eligible site gives a 95% upper bound of approximately **3.0e10 s^-1**, not useful evidence for slow etch kinetics. This is an analytical observation-window example, not a new AIMD result. Use longer/replicated sampling, validated rare-event methods or matched DFT/NEB+thermochemistry for activated events. Do not infer a barrier from one temperature without a justified prefactor.

Diffusion from AIMD/MD requires unwrapped coordinates and a diffusive fitting interval: `D = slope(MSD)/(2d)`, or the corresponding velocity-autocorrelation integral. Reject ballistic, caged or drifting regimes and estimate uncertainty across independent samples. [LAMMPS diffusion guidance](https://docs.lammps.org/Howto_diffusion.html).

## Calculated rates and next gate

[Source coadsorbate thermal rates](../data/kinetic_audit/published_coadsorbate_rates.csv) evaluate the existing Jung 2020 R6/R7 parameterization at 300/400/450/500 K. They describe its fluorinated clusters and vibrational-assistance assumptions, not our new ideal-slab candidates. [Primary paper](https://doi.org/10.1116/1.5125569).

The next critical graph rates are the **13 bare-N and eight NH2 final-cleavage environments** responsible for the [12-cycle stall](multilayer_results/cycle_diagnosis.json). For each, obtain matched relaxed endpoints, an index-one saddle connected to both basins, consistent energy/force accuracy and full vibrational treatment; then record geometry applicability, uncertainty and artifact hashes in the [evidence gate](../plasma_surface/rate_evidence.py). No energy from a rigid scan or an unconverged optimizer is promoted to a rate.

```sh
python scripts/audit_kinetic_parameters.py
python -m pytest tests/test_thermal_kinetics.py -q
```
