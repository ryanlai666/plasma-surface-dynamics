# Scientific review and a defensible route to publication

## Verdict and narrower research question

This repository is a reproducible research prototype, **not a publication-ready predictive ALE model**. Numerical verification, a large reaction graph and animated trajectories do not establish chemical accuracy. The present scientifically useful result is an audit of which assumptions dominate a bounded HF/SiN:H model, together with explicit failures of atomistic reference choices.

A defensible first paper would ask: **Which surface motifs and competing precursor/retention processes determine fluorination versus volatile Si release in a specified hydrogenated SiN film?** This is narrower than predicting all Si/SiNx ALE, deposition and plasma recipes. Si/Cl2/Ar+ energy-window physics should be a separate validated branch.

## Independent literature, different regimes

| Evidence | Insight | Consequence for this project |
|---|---|---|
| [Khumaini et al., 2024](https://doi.org/10.1016/j.apsusc.2024.159414): DFT on hydrogenated amorphous SiN | Film structure and H content matter; cleavage competes with salt formation/desorption. | Our crystalline bare nitride is not a substitute for that film. The enabled source pathways omit retention, so inventory exhaustion cannot establish physical self-limitation. |
| [Jung et al., 2020](https://doi.org/10.1116/1.5125569): strongly fluorinated nitride/oxide clusters in a remote-plasma model | HF/H2O and HF/HF(v=1) coadsorbates introduce distinct removal channels. | Include precursor populations and excitation/quenching in a separate matching regime; do not splice its barriers into the amorphous-film chain. |
| [Cheng and Hwang, 2021](https://doi.org/10.1016/j.apsusc.2020.148557): CH3F on terminated nitride | Their first-principles study questions intact-molecule chemisorption as the explanation for full-layer removal. | Neutral CH3F scans need radical-fragment and ion-created-site controls. A low scan minimum is not a sticking probability. |
| [Dwivedi et al., 2023](https://arxiv.org/abs/2305.09037): F2 on fluorine-terminated Si | Orientation trends involve termination and Si-Si cleavage; the tested reactive potential misses important barriers. | Compare matched coverage and bond-breaking paths before claiming Si(100)/Si(111) selectivity. |
| [Vella et al., 2022](https://graves.princeton.edu/publications/molecular-dynamics-study-silicon-atomic-layer-etching-chlorine-gas-and-argon-ions), [2025](https://doi.org/10.1021/acs.jpcb.5c01378): Si/Cl2/Ar+ | Ion fluence and competing chemical/physical sputtering govern ALE behavior. | Thermal motif kinetics cannot establish an ion-energy window; validate incident-energy, angle and fluence distributions separately. |
| [Gil et al., 2023](https://www.nature.com/articles/s41598-023-38359-4): remote plasma, methanol and heating | An experimental recipe produces oxide-selective removal through retained surface products. | Nitride/oxide selectivity is not a fixed material property; recipe and product retention can change its direction. |
| [Reuter and Scheffler, 2006](https://doi.org/10.1103/PhysRevB.73.045433): first-principles kMC | Spatial correlations are part of the modeled kinetics, beyond a lattice picture. | Our independent motif grid has no neighbor interactions. It must not be presented as spatially predictive. |
| [Meskine et al., 2009](https://arxiv.org/abs/0805.4356): rate-control analysis | Sensitivity can identify kinetically influential processes. | We now compute pulse-endpoint sensitivities; these are not steady-state degree-of-rate-control coefficients. |
| [Official MACE fine-tuning guidance](https://mace-docs.readthedocs.io/en/latest/guide/finetuning_guidance.html) | Validate the target observable, including barrier profiles and dynamics. | Hold out whole reactions/surface environments; near-duplicate path-frame splits would overstate generalization. |

These sources provide independent mechanistic and methodological evidence, not interchangeable numbers. [Machine-readable evidence map and access scope](../data/literature/research_evidence_map.json). Some publisher results were accessible as abstracts only; no inaccessible full text is claimed to have been reviewed.

## Improvements actually computed in this revision

### 1. Test the rigid-slab assumption

Eight HF/surface starts were relaxed with the upper substrate mobile. Seven meet the force threshold; one Si(111) start reaches the 300-step limit without convergence. **All retain a short HF bond** under the stated distance criterion. None establishes a dry-etch transition state.

A stronger diagnostic removes HF from each final geometry and evaluates the resulting frozen clean slab and molecular layer. The exact energy decomposition separates interaction from substrate reconstruction and reference effects. On Si(111), the very negative apparent adsorption difference is dominated by a lower-energy substrate configuration relative to the original clean-slab local minimum. It is not evidence for exceptionally strong HF chemisorption. Force convergence of the original reference was insufficient.

This changes how the results are used: flagged reference differences are excluded from quantitative adsorption claims and all kinetic parameterization. The next remedy is multistart reconstructed/passivated reference slabs and matching DFT checks, not treating these large negative values as a discovery. [Numerical audit, trajectories and decomposition](dry_etch_results/RELAXED_ADSORPTION.md).

### 2. Determine what actually controls the current kinetics

![Transient rate sensitivity and initial-motif dependence](species_kmc_results/kinetic_priorities.png)

At 400 K after a 2 s dose plus 1 s purge, the deterministic expected Si-release fraction is **0.3884**; the finite C++ ensemble gives **0.3871**. Local logarithmic sensitivity is **+0.242 for HF arrival**, **-0.081 for HF desorption**, and **+0.044 for P2c**, the largest individual chemical-pathway group. These values are calculated here, not measured or taken from a paper.

For each rate group, all matching events are multiplied by 1.1 or its reciprocal:

$$S_j(t_f)=\frac{\ln Y(k_j\times1.1,t_f)-\ln Y(k_j/1.1,t_f)}{2\ln1.1}.$$

P4 is perturbed jointly wherever the same source parameter is reused. Other rates and initial populations are held fixed. This is a local transient diagnostic, not a global uncertainty interval, a universal ranking or a thermodynamically constrained steady-state rate-control calculation.

**Actionable inference:** quantify HF delivery/sticking and precursor residence first; among the current enabled chemical groups, refine P2c before spending all the DFT budget on a numerically insensitive path. A missing pathway cannot appear in this sensitivity ranking, so the structural gap audit still matters.

### 3. Expose mixture confounding and false saturation

With the same rate constants, a pure NH-bridge inventory releases about **94.5%** of its Si centers, a pure pre-SiH2F inventory **97.7%**, and pre-SiHF2 **7.3%** by 3 s. Terminal-NH2, bare-N and Si-Si inventories release zero because their final pathways are disabled, not because selectivity has been experimentally established.

The weighted pure-motif predictions reproduce the mixture calculation numerically. This follows from the linear generator:

$$\dot{\mathbf n}=Q(t)\mathbf n,\qquad Y(t)=\sum_m w_mY_m(t).$$

There are no lateral interactions, no motif interconversion and no subsurface refill. The currently assumed initial fractions therefore set a **60% eventual topological removal ceiling**, even with unlimited exposure and all enabled barriers crossed. A plateau can be imposed by missing graph edges or finite inventory; it is not proof of ALE self-limitation. Spatial kMC adds trajectory fluctuations here, but its mean is already represented by the linear master equation.

**Actionable inference:** constrain motif weights independently using film composition/termination evidence and product or coverage measurements. Fitting a single removal curve cannot identify all motif weights and rates. Model new-site exposure and retention before interpreting cycle-to-cycle saturation.

[All finite differences, pure-motif results and provenance](species_kmc_results/kinetic_priorities.json) | [Runner](../scripts/analyze_kinetic_priorities.py).

### 4. Quantify why coadsorbate identity matters

The following rates are recomputed at 400 K from **Jung Table I, R6/R7**, using $k=A\exp[-(E_a/R)/T]$. They apply per occupied coadsorbate state, not per arbitrary bare site. [Author-hosted source](https://cpseg.eecs.umich.edu/pub/articles/JVSTA_38_023008_2020.pdf).

| Occupied state / channel | Si3N4 rate (s^-1) | SiO2 rate (s^-1) |
|---|---:|---:|
| HF + H2O / R6 | 34.34 | 135204.35 |
| HF + HF(v=1) / R7 | 4.985e9 | 1.159 |

The contrasting intrinsic channel preference shows why a single neutral HF pathway cannot determine selectivity. Occupancy, incident flux, desorption and retained products still determine the net process. These values are **not inserted into the 55-event model**. [Transcribed parameters](../data/literature/jung2020_coadsorption_rates.json), [calculated rates](species_kmc_results/literature_channel_rates.csv), [reproduction](../scripts/compare_literature_channels.py).

## Concrete next calculations and acceptance criteria

These are proposed requirements, not completed results or universal standards.

| Priority | Calculation or data needed | Decision it resolves | Acceptance requirement |
|---|---|---|---|
| 1 | Independently relaxed reconstructed Si; passivated and amorphous SiNx:H reference ensembles | Are apparent adsorption energies chemical binding or reference artifacts? | Multiple starts, cell/thickness convergence, consistent reference basins; DFT forces on chosen structures |
| 2 | Matched HF adsorption/desorption and P2c IS/TS/FS on a specified a-SiN:H motif | Which measured flux and precursor occupancy produce reaction? | Explicit energy zero, curvature/endpoint checks, thermal corrections and site-density/sticking model |
| 3 | HF-HF and HF-H2O complexes on matching fluorinated N/O sites; retained fluoride/salt states | Can omitted channels change branching or self-limitation? | Stable intermediates, balanced site/atom inventory, connected paths; keep excitation and salt-removal regimes distinct |
| 4 | Paired trajectories with intact CH3F, relevant fragments and ion-created defects | Is molecular dosing sufficient for modification? | Matched initial surfaces and controlled incident distributions; charge/spin treatment for radicals |
| 5 | Multilayer exposure and neighboring-site events | Is a plateau finite-inventory exhaustion or physical passivation? | Explicit connectivity, surface-only accessibility, ledgers including newly exposed material; no automatic unvalidated site refill |
| 6 | Matched public etch/coverage/product data across dose, temperature and purge | Does the model predict held-out observations? | Reserve whole conditions for validation; compare against the simpler independent-site model and quantify uncertainty |

An adsorption hazard should ultimately be tied to gas delivery, for example $k_{ads}=S A_{site}p/\sqrt{2\pi m k_BT_g}$ for a thermal gas with appropriate assumptions. A thermal reaction rate needs a consistent transition-state treatment, e.g. $k=(k_BT/h)\exp[-\Delta G^\ddagger/(k_BT)]$. Reverse rates require the same thermochemical references; an irreversible product sink is justified only for the specified dilute-flow regime. These are modeling requirements, not fitted replacements for the present assumptions.

At 400 K, a **0.1 eV barrier error changes an Arrhenius rate by about 18.2-fold**. A factor-of-three rate target corresponds to roughly 0.038 eV if barrier error were the only uncertainty. Therefore generic energy/force RMSE or a visually smooth NEB cannot justify precise kinetic predictions. Use reaction-specific DFT energy/force tests, source-model disagreement, correlated uncertainty propagation and validation against the intended observable.

The proposed publication contribution is a tested mechanistic discrimination with uncertainty and observable consequences. It is not the number of molecules screened, number of GIFs, network size or local C++ speedup.
