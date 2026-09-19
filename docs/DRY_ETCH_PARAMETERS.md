# Dry-etch transition states and parameter mapping

Current assessment: [multi-source scientific critique and calculated priorities](SCIENTIFIC_REVIEW.md). The current species model has no spatial interactions or subsurface refill; candidate enumeration does not establish completeness.

The active atomistic evidence concerns surface reactions. This update **reuses published DFT transition states and rate laws**; it does not claim new periodic DFT or a calibrated Cl2/Ar ALE window. The unrelated NH3 inversion example has been moved to `archive/molecular_diagnostics/` and is excluded from dry-etch evidence.

## What can constrain which model

| Chemistry | Acquired evidence | Appropriate parameter use | Missing |
|---|---|---|---|
| F2 + F-terminated Si, four facets | Four DFT/TST barriers and fitted prefactors/probabilities | Implemented first-event gas-surface hazard, within 298.15-1000 K | Complete removal sequence, F-escape branching and new-layer exposure |
| HF + amorphous SiN:H | 16 fluorination Ea/Ephy/DeltaE entries; two salt-formation barriers | State-specific thermal reaction candidates; conditional prefactor sensitivity | Exact energy-reference audit, vibrational prefactors and coordinates |
| SiCl4 on reconstructed Si(100) | 17 stationary structures, six matched IS/TS/FS paths and energies | Readsorption, recombination and dimer-flip channels | Full NEB images, Hessians, raw forces and cell vectors |
| HF impacts on amorphous Si3N4 | Six published ML-MD yields and seven state-resolved area coefficients | Chemistry-specific impact-model reference | Complete Table S4 operating conditions; transferable uncertainty; raw model/training files |
| Cl2/Ar ALE and Cl diffusion | Previously collected MD context and facet-specific diffusion barriers | Mechanism constraints and future explicit hops | Matched ion-yield ensembles and modifier-specific TS/prefactors |

The [machine-readable audit](../data/literature/ale_parameter_audit.json) checks **all ten current baseline parameters**. None has been silently replaced. In particular, the illustrative 0.65 eV modifier-desorption barrier has no identified elementary event; it is not validated by an HF cleavage barrier. Nitrogen composition alone does not determine any of these rates.

## Silicon/F2: implemented reaction-specific rates

[Dwivedi et al., author manuscript, Table 1 and equations 5a/5b](https://arxiv.org/abs/2305.09037v2) report first-step dissociative chemisorption on F-terminated Si cluster models. The parameterization is broken-symmetry unrestricted B3LYP, with a Si32F32 cluster. This is chemical dry etching rather than a Cl2/Ar plasma ALE cycle.

| Facet | Ea (eV) | A0 (m3/s) | n1 | gamma0 | n2 |
|---|---:|---:|---:|---:|---:|
| reconstructed (100) | 0.13 | 8.82e-21 | 2.20 | 0.0011 | 1.60 |
| unreconstructed (100) | 0.31 | 1.23e-20 | 2.20 | 0.0016 | 1.70 |
| (110) | 0.35 | 7.52e-20 | 2.50 | 0.0096 | 1.96 |
| (111) | 0.57 | 3.90e-20 | 2.50 | 0.0050 | 1.96 |

$$k(T)=A_0(T/298.15)^{n_1}\exp[-E_a/(k_BT)],\qquad r=k(T)n_{F_2},\qquad n_{F_2}=p_{F_2}/(k_B^{SI}T_g).$$

The reproduction assumes gas and surface temperatures are equal. `k` is bimolecular (m3/s), while `r` is the per-site hazard (s^-1). A separate published fit gives the per-impact probability `gamma=gamma0*(T/298.15)^n2*exp(-Ea/kBT)`; it is not multiplied into `k` again. The implementation rejects extrapolation outside the source fit range. At fixed conditions, each initially susceptible site reacts at most once, so `f(t)=1-exp(-r*t)`. It does not remove a column, reset a fresh site or assign nm/cycle.

Reproduce with `python -m plasma_surface.dry_etch`. The [results](dry_etch_results/REPORT.md) contain 484 temperature/facet points and Gillespie verification. The paper discusses F escape versus two-F retention after the first TS, which can change isotropy; its branching is not determined by the four barriers.

## SiN:H/HF: all fluorination paths and competing fragments

The [Khumaini et al. source poster](https://avssymposium.org/ALD2024/Sessions/SupplementalDocumentDownload/80689?sessionId=79078), associated with [the DFT paper](https://doi.org/10.1016/j.apsusc.2024.159414), contains 16 pathways. The [full numerical extract](../data/literature/sin_hf_pathways.csv) retains the source motif notation, stage, Ephy, DeltaE, Ea and product separately; [provenance](../data/literature/sin_hf_pathways_provenance.json) records the PDF hash.

| Event family | Path IDs | Ea (eV) |
|---|---|---|
| NH3 release from Si-N cleavage | P1a, P2a, P3a | 0.45, 0.68, 0.88 |
| Competing Si-H cleavage / H2 release | P1b, P2b, P3b | 2.10, 1.80, 1.54 |
| Bridge cleavage, first fluorination | P1c, P1d | 0.02, 0.72 |
| Bridge/Si-Si cleavage, second stage | P2c, P2d, P2f | 0.79, 0.75, 0.79 |
| Bridge cleavage, third stage | P3c, P3d | 0.74, 0.90 |
| SiH2F2, SiHF3, SiF4 release | P2e, P3e, P4 | 0.81, 0.96, 0.53 |

These local motifs are not mutually exclusive exits from one identical initial state. They must not be normalized into product fractions without motif populations, adsorbed-HF availability and connectivity. NH3 remains a relevant **product of Si-N cleavage**; its gas-phase umbrella motion supplies no etch parameter.

The run also outputs 192 conditional sensitivity values using `nu=1e11, 1e12, 1e13 s^-1` and 300/400/500/600 K. These prefactors are assumptions, not measurements. `nu exp(-Ea/kBT)` is a conditional hazard for the occupied source reaction complex. Full energy-reference conventions, entropy and HF adsorption equilibrium remain unresolved, so this table is not a forecast of film etching.

The [salt extract](../data/literature/sin_salt_energetics.json) adds formation barriers 0.40 and 0.07 eV and separate desorption energetics. AFS can release directly, decompose to NH3/HF/SiF4, or pass through retained NH5F2. The reported 1.97-4.79 eV desorption energies are **not automatically TS barriers**. Source free-energy crossover temperatures are condition dependent, not recommended process setpoints.

## Chlorinated silicon: source-based stationary-point animations

[Zhang, Zhu and Li](https://doi.org/10.3390/sym15010213) supply 117-atom stationary structures for SiCl4 dissociation on Si(100)-c(4x2). We parsed their coordinates, preserved atom order and reconstructed the published pathway assignments. The energy method is periodic CP2K PBE-D3(BJ)/TZVP-MOLOPT-GTH (550 Ry); the separately described B3LYP cluster analysis is not the source of the plotted reaction energies.

The [six GIFs and barrier table](dry_etch_results/TRANSITION_STATES.md) run in the reverse direction: adsorbed SiCl3 + Cl recombine into a physisorbed SiCl4 complex (plus one dimer-flip path). This is a relevant competing surface channel, **not a demonstrated net removal of substrate Si**. The initial, transition and final structures are source data; the animations contain only the published IS, TS and FS, without invented intermediate frames. No unreported saddle Hessian check is claimed.

For each source triplet, `Ea_forward=E_TS-E_IS`, `Ea_reverse=E_TS-E_FS` and `DeltaE=E_FS-E_IS`. We check `Ea_forward-Ea_reverse=DeltaE`. Source energies are referenced to separated SiCl4 plus clean surface, in kcal/mol; subtract first, then divide by 23.060547830619 to obtain eV. A submerged TS relative to gas is not a negative barrier from the adsorbed initial state.

Actual substrate-removal evidence is kept distinct: [the Si(100) chlorine-etch study](https://doi.org/10.1016/S0039-6028(99)00610-X) discusses Cl transfer, SiCl2 back-bond breaking and higher chlorination, while [the Si(111) photo-ALE study](https://doi.org/10.1016/j.mssp.2022.107169) reports 2.5159 eV chlorinated-Si removal versus 6.1832 eV bare-Si removal. Facet, coverage and activation method differ; these cannot calibrate a universal ion threshold.

## Impact yields and the missing connection to plasma ALE

[Hong et al., Supporting Information Tables S3-S4](https://doi.org/10.1021/acsami.4c07949.s001) provide six HF-impact yields and seven species-resolved reaction/desorption coefficients. We inspected the public 11-page supplement. It does **not** provide a model checkpoint, raw training set, or matched TS coordinates. [Extract and availability audit](../data/literature/sin_hf_mlip_reference.json), [yield table](../data/literature/sin_hf_impact_yields.csv).

The yields are in **Si3N4/HF**, not Si atoms/Ar ion. The S4 coefficients are **area per incident HF**, despite being called probabilities; multiplication by a matched flux gives s^-1. They cannot be copied into dimensionless `chemical_yield_scale`. The angular and termination dependence makes them useful test cases for a dedicated HF-impact extension once all source conditions are recovered.

For plasma-driven removal, the required closure is channel-specific:

$$r_j=\frac{\Gamma_i}{n_s}\int Y_j(E,\alpha,\mathbf{s})\,f_i(E,\alpha)\,dE\,d\alpha,$$

where `f_i` is a normalized incident energy/angle distribution and `s` contains coverage, composition, H content and defects. MD/beam counts must distinguish modifier stripping, substrate removal, implantation and retained fluorocarbon. Thermal rates require matched minima, saddle and vibrational free energies. These are different measurements/calculations, even when both energies are expressed in eV.

## Next missing atomistic calculations, in order

1. On chlorinated Si, distinguish Cl loss from SiCl2 detachment at specified coverage; include the neighboring Cl-transfer step and defect sites. Obtain matched slabs and cell/constraints before CI-NEB.
2. On SiNx:H, refine HF cleavage at terminal NH2 and bridging NH/N sites, then SiHF3/SiH2F2/SiF4 release for several amorphous environments. Reuse the source motif definitions rather than extrapolating a single zero-H composition card.
3. Add CFx retention and N removal through HCN/FCN/CNx, with O-assisted carbon removal as a separate chemistry. Sources and proposed discriminating experiments are in the [fragment guide](REACTION_CANDIDATES.md).
4. Only after checking endpoint connectivity, force convergence and the unstable mode, compute prefactors/free-energy corrections and fit a state-resolved kinetic network. The acquired SiCl4 XYZ files omit cell vectors; they are not ready-to-run periodic DFT inputs.
