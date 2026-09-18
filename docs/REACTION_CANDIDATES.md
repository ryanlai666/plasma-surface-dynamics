# Fragments, side reactions, and discriminating experiment sets

The [machine-readable catalog](../configs/reaction_candidates.json) contains 35 species/bookkeeping units and 21 atom- and charge-balanced candidate reactions. These are a research shortlist, not a validated kinetic mechanism. No channel probabilities or rate constants are assigned. Literature-supported products and hypothetical side reactions are explicitly labeled per entry. The [checker](../scripts/check_reaction_catalog.py) verifies conservation; it does not prove a reaction occurs.

## Candidate fragments by chemistry

| Process family | Candidate volatile products or intermediates | Competing surface/wall behavior | Evidence |
|---|---|---|---|
| Si + HCl or chlorine chemistry | SiCl2, SiCl4, chlorohydrosilanes such as SiHCl3/SiH2Cl2; HCl, H2 and Cl-containing radicals | Hydrogen/chlorine passivation; loss of Cl before Si removal; radical recombination and wall memory | [Si-HCl DFT/kMC study](https://doi.org/10.1016/j.apsusc.2024.159836); [Cl2/Ar ALE MD study](https://doi.org/10.1021/acs.jpcb.5c01378). Product identities and relative yields depend on the selected chemistry |
| Hydrogenated SiN + HF | NH3, SiF4, SiHF3, SiH2F2; stepwise fluorinated surface sites | NH4F and ammonium fluorosilicate (AFS) retention; delayed decomposition/desorption | [Amorphous SiN:H/HF DFT study](https://doi.org/10.1016/j.apsusc.2024.159414). Existing zero-H parameter cards do not represent this mechanism |
| SiN + CH3F/CF4 with ion activation | CH3/F from CH3F; CF2/F from CF4; candidate secondary CF3 and NHx fragments | Fluorocarbon retention/passivation versus removal; precursor activation and ion-induced damage | [AIMD/tight-binding impact study](https://doi.org/10.1063/5.0155929). Secondary channels remain hypotheses; [DFT CH3F study](https://doi.org/10.1016/j.apsusc.2020.148557) challenges simple direct molecular adsorption |
| Background/measurement controls | H2O, HF/HCl, NHx and SiFx fragments in the measurement chain | Wall recombination, moisture-sensitive residue, ionizer fragmentation, delayed chamber release | Candidate confounders to test; no reactor-specific rates established |

An illustrative balanced nitride removal equation is `Si3N4 + 12 HF -> 3 SiF4 + 4 NH3`. It is a **net balance**, not one elementary reaction, and cannot predict preferential removal or hydrogen dependence. Salt retention can be represented by `SiF4 + 2 NH4F -> (NH4)2SiF6`. Surface atoms and salts are bookkeeping units in the catalog, not isolated molecules to pass blindly to an ML calculator.

For initial molecular calculations, stable closed-shell fragments use neutral singlet inputs; H, F, Cl, CH3, CF3, NH2 and SiF3 use neutral doublet candidates. These reference-state choices do not exclude excited or charged plasma species. Ar+ is a separate charged doublet entry. Surface/slab charge and spin require their own definition.

## Proposed experimental sets

Choose actual energy, temperature, dose and pressure levels from the apparatus's characterized operating range. The project's illustrative 35 eV recipe and hypothetical nitride thresholds are not experimental setpoints. Record measured ion energy/angle distribution and fluence, neutral dose, surface temperature, composition, hydrogen content, film thickness and chamber history.

| Set | Comparison | Measurements to correlate | Mechanism it can distinguish |
|---|---|---|---|
| E0: half-cycle controls | Untreated, neutral-only, ion-only, full sequential cycle; include purge-only and a matched blank substrate | Thickness/mass change, surface coverage, phase-resolved products | Synergistic removal versus spontaneous chemistry, sputtering or instrument background |
| E1: modification saturation | Short/intermediate/long neutral dose at fixed removal fluence; then vary removal fluence at saturated dose | Uptake plateau, EPC, residual Cl/F, Si-containing product evolution | Finite modified-layer removal versus continued penetration, damage or precursor carryover |
| E2: Si chlorination selectivity | Compare chlorine uptake and Si loss through the removal half-cycle; repeat with consistent H termination | Calibrated SiClx product patterns; surface Cl; Si removal | Preferential Cl stripping versus useful substrate removal; missing channels in the two-state model |
| E3: nitride composition and H | SiN0.8, SiN1.0 and Si3N4 where available; independently characterize and vary H content rather than attributing it to N/Si | Si/N/F ratios, NH3/SiF4-family products, film density and thickness | Composition effects versus hydrogen-assisted chemistry and conversion-factor artifacts |
| E4: salt retention | Compare post-modification, post-purge and post-removal surfaces; inspect controlled thermal-release behavior | Mass gain/loss, N-H/F and AFS-compatible signatures, correlated delayed NH3/HF/SiF4 | A nonvolatile intermediate or residue versus immediate volatilization |
| E5: fluorocarbon competition | Matched-dose precursor-only versus ion-assisted exposure; compare CH3F/CF4 with composition tracked | Surface C-F, C/Si, F/N, net mass and removal; precursor fragments | Passivation/polymer buildup versus chemical activation and etching |
| E6: diffusion/redistribution | Same uptake followed by different dwell times and temperatures with no new incident flux | Spatial or chemical redistribution and unchanged total uptake where applicable | Diffusion versus desorption; distinguish a redistribution timescale from removal |
| E7: background and memory | Blank substrate/chamber exposures; repeated purge durations and run order | Product tails, H2O/O-containing signals, wall-history dependence | Chamber residence time, contamination, and wall reactions masquerading as surface kinetics |

For E6, an Arrhenius slope alone is not proof of a specific diffusion pathway: a transport or desorption process can produce an apparent activation energy. Use the specific reconstructed surface before comparing to the [Cl diffusion values](../data/literature/diffusion_barriers.csv).

Use repeated reference runs to assess drift, randomize condition order where feasible, and retain whole recipes as held-out validation sets. Begin with a small discriminating design rather than a full factorial across every variable. Matched half-cycle controls are more informative for mechanism selection than merely collecting more EPC points.

## Product detection and atom balances

Gas analysis should use calibrated parent/fragment patterns and phase timing. A SiF3 signal can come from ionization of SiF4; it is not sufficient evidence that SiF3 desorbed from the film. Likewise, NHx channels can overlap other backgrounds. Chlorine isotope patterns help with assignments but do not remove the need for calibration. Couple gas signals with surface spectroscopy and independent thickness or mass changes; do not infer stoichiometry from one peak.

Track Si, N, H, F/Cl and C balances separately, including retained film, volatile products and unmeasured channels. A persistent C-F layer with mass gain suggests that a deposition/passivation channel is needed; delayed N/H/F release suggests an additional retained-state pool. Neither can be represented just by lowering a single sticking number.

## Computational sequence for these sets

1. Use the NH3 inversion result to verify the saddle/force/Hessian workflow on a small nitride-related gas product. It does not estimate NH3 release from a nitride surface.
2. Rank **candidate calculations**, not predicted reactor yield: adsorbed HF attack on a hydrogenated Si-N bridge; fluorination of Si-H groups; Cl migration on a defined reconstructed Si surface; fragment-assisted Si-N cleavage; salt formation/decomposition where relevant.
3. For cheap initial screening, the catalog includes successive `SiH4 + HF -> SiH3F + H2`-type molecular exchanges. These are atom-balanced surrogate hypotheses, not established surface pathways. Their TS and reaction energies are not yet calculated.
4. For each selected real surface event, compute matched reactant, saddle and product energies, verify the saddle and connectivity, then evaluate local-environment and charge/spin dependence. Check ML candidates with DFT before constructing rates.
5. Fit a richer state model only after observations discriminate bare, modified, damaged, fluorocarbon-covered and salt-retaining states. Keep diffusion as explicit hops rather than folding it into etch yield.

Reproduce the catalog check with `python scripts/check_reaction_catalog.py`. Links: [data acquisition and diffusion theory](ATOMISTIC_DATA.md), [computed molecular results](barrier_results/REPORT.md), [earlier mechanism evidence](LITERATURE.md).
