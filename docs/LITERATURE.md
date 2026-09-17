# Literature, hypotheses, and parameter provenance

This document separates published atomistic evidence from this repository's coarse-grained assumptions. Literature checked 2026-09-16. No new DFT, AIMD, or barrier calculations were performed for this update. The existing MACE calculations are structural probes, not DFT reaction validation. Publisher abstracts/previews and the public author poster below support the extracted claims; inaccessible full texts were not treated as reviewed.

## Evidence and applicability

| ID and primary source | Method and system | Reported evidence | Appropriate use and limitation |
|---|---|---|---|
| L1: [Mechanism of photo-assisted atomic layer etching of chlorinated Si(111) surfaces: Insights from DFT/TDDFT calculations (2023)](https://doi.org/10.1016/j.mssp.2022.107169) | DFT/TDDFT; chlorinated Si(111), photo-assisted removal | Reported removal/desorption energy changes from 6.1832 to 2.5159 eV after chlorination | Supports modification-dependent removal energetics. These values concern the paper's reaction coordinate; neither is the incident Ar-ion threshold or the kernel's modifier-desorption barrier. |
| L2: Wang and Fang, [Ab initio simulations of ultrashort laser pulse interaction with Cl-Si(100) (2023)](https://doi.org/10.1039/D3CP02388E) | Real-time TDDFT coupled to nuclear dynamics under laser excitation | Electronic excitation changes SiCl desorption dynamics | An example of excited-state first-principles dynamics. Laser-driven Ehrenfest dynamics should not be presented as thermal ground-state AIMD or as an Ar sputtering calibration. |
| L3: Cheng and Hwang, [Dissociative chemisorption of methyl fluoride and its implications for atomic layer etching of silicon nitride (2021)](https://doi.org/10.1016/j.apsusc.2020.148557) | PBE/PAW DFT; H-terminated N-rich beta-Si3N4 surface | Initial CH3F sticking below 1e-7 in the absence of dangling bonds; direct chemisorption alone does not explain the proposed complete-layer removal | Challenges interpreting the present S=0.18-0.25 as molecular CH3F sticking. The current radical species is unspecified, so these are different channels. Surface activation and precursor fragmentation need explicit states. |
| L4: Cheng and Hwang, [Low-energy argon ion bombardment-induced decomposition of physisorbed hydrofluorocarbons on silicon nitride surfaces (2023)](https://doi.org/10.1063/5.0155929), [author abstract in PubMed](https://pubmed.ncbi.nlm.nih.gov/37417756/) | Combined tight-binding MD and AIMD; adsorbed HFCs, ion energies up to 35 eV | Collision-assisted surface reactions and direct decomposition depend on impact conditions; collision-assisted routes dominate near 11 eV in the studied cases | Supports adding adsorbate and impact-conditioned pathways. Adsorbate decomposition onset is not a substrate etch threshold; the 35 eV study range does not validate our 35 eV animation. |
| L5: Khumaini et al., [Etching mechanism of amorphous hydrogenated silicon nitride by hydrogen fluoride (2024)](https://doi.org/10.1016/j.apsusc.2024.159414), [author institutional abstract](https://sejong.elsevierpure.com/en/publications/etching-mechanism-of-amorphous-hydrogenated-silicon-nitride-by-hy/) | DFT reaction pathways on an amorphous Si-rich nitride with about 25 at.% H, prepared with MD and DFT | Si-N/Si-Si cleavage pathways have barriers at or below 0.90 eV; products and salt formation matter | Supports an HF-specific, hydrogen- and fluorination-resolved network. Current cards have zero H and no salt/product states; transferring these barriers into the generic kernel is unjustified. The source's use of MD plus DFT is not by itself evidence that every trajectory was AIMD. |
| L6: [Si-Cl2-Ar+ Atomic Layer Etching Window: A Fundamental Study Using Molecular Dynamics Simulations and a Reduced Order Model (2025)](https://doi.org/10.1021/acs.jpcb.5c01378) | Classical MD and reduced-order modeling; normal-incidence Ar ions | A narrow window around 15-20 eV in that study; selective Cl loss limits EPC | A useful challenge for this model's broad 15/40 eV threshold separation. This is MD, not AIMD. It does not establish universal thresholds across ion angles, force fields, and surface conditions. |
| L7: Kroll, [Structure and reactivity of amorphous silicon nitride investigated with density-functional methods (2001)](https://doi.org/10.1016/S0022-3093(01)00676-7) | Constructed amorphous networks followed by density-functional relaxation and ab initio MD | Chemically ordered networks and hydrogenated structures provide atomistic structural models | Supports using representative amorphous structures when extending the project. Does not supply a composition-dependent etch-rate law or validate crystalline vacancy probes as amorphous films. |

The two-state mechanism is a useful research baseline, but L3-L5 identify chemistry that it cannot resolve. This evidence motivates model extensions rather than retrospective claims that the current constants were derived from first principles.

## Small public numerical extract

The authors' [ALD 2024 poster, page 1](https://avssymposium.org/ALD2024/Sessions/SupplementalDocumentDownload/80689?sessionId=79078) associated with L5 publicly tabulates pathway energies. Two factual entries are transcribed in [sin_hf_barriers.csv](../data/literature/sin_hf_barriers.csv), with [provenance and interpretation](../data/literature/provenance.json):

| Source pathway | Context | Activation energy Ea (eV) | Reaction energy Delta E (eV) |
|---|---|---:|---:|
| P1d | First fluorination; bridging Si-N-Si + HF reaction | 0.72 | -0.97 |
| P4 | Fourth fluorination; SiF4 release from an NH-bridged fluorinated site | 0.53 | -0.77 |

These entries are examples for a future reaction catalog. They are not statistical uncertainty bounds, a complete dataset, or calibrated inputs. Preserve the poster's reference convention when reusing Ea; do not reinterpret energies relative to physisorbed reactants as gas-to-surface sticking barriers without reconstructing the full pathway. Prefactors and free-energy corrections are not supplied in this extract. The CSV is not loaded by the solver. Attribution is retained; the paper and poster are not redistributed or assigned a new license.

## Current parameters and what would constrain them

All numbers in this table are implementation values, not measurements or literature fits. SiNx cards inherit omitted entries from `Parameters`.

| Parameter | Si | SiN0.8 | SiN1.0 | Si3N4 | Evidence required to replace the assumption |
|---|---:|---:|---:|---:|---|
| Effective sticking S | 0.25 | 0.22 | 0.20 | 0.18 | Species- and surface-specific uptake versus fluence; impact or adsorption calculations. L3 rules out a naive CH3F interpretation. |
| Chemical threshold Ec (eV) | 15 | 22 | 25 | 30 | Energy-resolved ion-impact yields on modified surfaces, matched to beam and film conditions |
| Physical threshold Es (eV) | 40 | 45 | 48 | 50 | Bare-surface sputter yields; L6 cautions against treating 15/40 eV as a validated window |
| Chemical yield scale Ac | 0.8 | 0.65 | 0.60 | 0.55 | Removal counts per incident ion at controlled modification coverage |
| Physical yield scale As | 0.15 | 0.15 | 0.15 | 0.15 | Bare-film ion-only controls and matched MD ensembles |
| Modifier desorption Ed (eV) | 0.65 | 0.65 | 0.65 | 0.65 | Identify the departing species and initial/final states, then calculate or measure its barrier |
| Desorption prefactor nu (s^-1) | 1e6 | 1e6 | 1e6 | 1e6 | Vibrational/free-energy analysis or temperature-dependent kinetics; do not fit Ed independently of nu without identifiability checks |
| Site density ns (m^-2) | 7e18 | 7e18 | 7e18 | 7e18 | Surface-site definition, exposed area, film structure and density |
| Thickness increment ell (nm) | 0.136 | 0.136 | 0.136 | 0.136 | Event stoichiometry and removed volume per area; inherited Si-like scale is illustrative for nitride |
| Deposition sticking | 0.2 | 0.2 | 0.2 | 0.2 | A separate deposition chemistry; inactive when precursor flux is zero |

The assumed ordering with N/Si is not fitted, interpolated from literature, or inferred from L1-L7. There is no established monotonic relationship here. In particular, the 0.65 eV desorption barrier must not be described as validated simply because it is numerically near an unrelated HF bond-cleavage barrier.

## Implemented equations

Let Gamma_r, Gamma_i and Gamma_p be incident radical, ion and precursor fluxes in m^-2 s^-1, ns the column density, T the surface temperature, and E the monoenergetic ion energy. The five per-eligible-column hazards in s^-1 are

$$
a=\frac{S\Gamma_r}{n_s},\quad
 d=\nu e^{-E_d/(k_BT)},\quad
 c=\frac{\Gamma_iY_c(E)}{n_s},\quad
 s=\frac{\Gamma_iY_s(E)}{n_s},\quad
 g=\frac{S_p\Gamma_p}{n_s}.
$$

$$
Y_j(E)=A_j\left[\sqrt{E/E_j}-1\right]_+,\qquad j\in\{c,s\},\qquad [x]_+=\max(x,0).
$$

The yield law is an illustrative closure, not a DFT-derived law. The model does not simulate individual ion trajectories. Y is an effective removal-count factor; it is not constrained to be a Bernoulli probability. The conversion from atom-resolved yields to column-removal units must be defined before calibration.

For N columns and n modified columns, total event propensities are

$$
(\alpha_a,\alpha_d,\alpha_c,\alpha_s,\alpha_g)
=(a(N-n),dn,cn,sN,gN).
$$

With independent uniform draws u1,u2, Gillespie sampling uses

$$
\Delta t=-\ln(u_1)/\alpha_0,\qquad
P(j)=\alpha_j/\alpha_0,\qquad \alpha_0=\sum_j\alpha_j.
$$

A cumulative propensity comparison with u2 selects j. Events beyond a phase boundary are discarded and rates are updated. For theta=n/N and positive net removed thickness D, the exact expectations satisfy

$$
\dot\theta=a-K\theta,\quad K=a+d+c+s+g,\qquad
\dot D=\ell(c\theta+s-g).
$$

During a constant-rate phase of length t with K>0,

$$
\theta(t)=\frac aK+\left(\theta_0-\frac aK\right)e^{-Kt},\qquad
I_\theta(t)=\frac{at}{K}+\left(\theta_0-\frac aK\right)\frac{1-e^{-Kt}}{K},
$$

$$
\Delta D=\ell\left[cI_\theta(t)+(s-g)t\right].
$$

For K=0, theta is constant and no events occur. EPC is the difference in cumulative D between consecutive cycle ends, not necessarily the steady-cycle limit. These equations match [model.py](../plasma_surface/model.py) and the C++ solver.

## Falsifiable hypotheses

| Hypothesis | Consequence within the model | Test or observation that could require revision |
|---|---|---|
| H1: One exposed site can be modified once before removal | With adsorption alone, theta=1-(1-theta0)exp(-at); uptake saturates | Multilayer penetration, subsurface modification, or several distinct uptake timescales |
| H2: Modification gates chemical removal | With radicals off and d=s=g=0, removal is ell*theta0*(1-exp(-ct)); it saturates at ell*theta0 | Continued removal after exhausting initial coverage, damage-assisted reactivation, or substantial ion-only etch |
| H3: Columns have identical, independent hazards | Spatial arrangement and height do not influence event rates | Correlated roughening, facet dependence, neighbor effects, transport limitations |
| H4: Nitrogen content decreases S and increases thresholds | The selected cards tend to remove less at the shared 35 eV recipe | Matched composition experiments or reactive trajectories give opposite/nonmonotonic trends |
| H5: A single Arrhenius modifier-loss channel suffices | Log(d) is linear in inverse T for fixed nu and Ed | Multiple products, coverage-dependent barriers, or nonthermal loss |

H1-H3 are mathematical/model assumptions; H4 is an invented material sensitivity scenario; H5 is an uncalibrated kinetic closure. Current tests verify consequences of these assumptions, not their experimental truth.

## How DFT and AIMD would supply rates

The following workflow is proposed, not implemented by the present animation or MACE demo.

1. **Define chemistry and structures.** Separate Si/Cl2/Ar from SiNx/HFC/Ar and SiNx:H/HF. For nitrides, use several amorphous slabs for each composition, hydrogen fraction, density, and termination. Retain slab area, active-site counts, preparation protocol and uncertainties.
2. **DFT reaction paths.** Relax reactants/products and transition states (for example with NEB and a saddle-point check). Record reaction energy Delta E=E_products-E_reactants and barrier E_dagger=E_TS-E_reactants with a common energy reference. Adsorption energy alone is not an activation barrier. Check convergence, slab thickness, vacuum, functional and spin/charge assumptions.
3. **Thermal kinetics.** For an equilibrated activated step, use transition-state theory with transmission coefficient kappa and activation free energy:

$$
k_r(T)=\kappa_r(T)\frac{k_BT}{h}\exp\!\left[-\frac{\Delta G_r^\ddagger(T)}{k_BT}\right].
$$

Use consistent energy units (eV or J). Alternatively use harmonic transition-state prefactors and potential-energy barriers with the matching vibrational corrections. Do not mix a free-energy barrier with an independently entropy-corrected prefactor. For two reversible equilibrated states at the same reservoir conditions, forward/reverse rates should reproduce their free-energy difference. Plasma-driven channels need not obey equilibrium detailed balance.

4. **AIMD/validated reactive MD impacts.** Sample incident energies, angles, surface configurations and independent impact locations. Estimate adsorption probability S=N_ads/N_inc and etch yield Y=N_removed/N_inc with sampling uncertainty. Short trajectories alone do not measure slow thermal rates; missing rare events are not evidence of a zero barrier or zero long-time rate. Check whether energy loss, electronic excitation and charged impacts exceed the force-field approximation.
5. **Fold in incident distributions.** For a normalized joint energy/solid-angle distribution f(E,Omega), a future impact rate is

$$
k_{r,\mathrm{ion}}(q)=\frac{\Gamma_i}{n_s}\int Y_r(E,\Omega,q)f(E,\Omega)\,dE\,d\Omega.
$$

Here q is the local surface state. The current kernel replaces the integral by one energy and omits q beyond bare/modified. Impact energies are not thermal activation energies, even though both may be expressed in eV.

6. **Build and validate the catalog.** Add explicit Si/N/H/F/Cl states and product/element balances before importing the public HF barriers. Include modification-depth or damage states if the selected mechanism needs them. Fit only identifiable remaining parameters against matching public measurements, then validate held-out energies, fluences, temperatures and compositions. Propagate structural, DFT, impact-sampling, fitted-parameter and model-form uncertainties separately.

The public geometry archive and current MACE checkpoints provide starting structures and comparative predictions. They do not establish reactive halogen/ion accuracy. New DFT reference energies, forces, and pathways would be required to train or validate a reactive ML potential. HiPRGen may propose candidate reactions but cannot assign surface-specific barriers by itself.
