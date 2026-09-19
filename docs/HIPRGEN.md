# HiPRGen pilot and intermediate coverage

Update: a separate [multilayer bond-graph prototype and evidence-gated rate queue](multilayer_results/REPORT.md) has now been implemented. Statements below about independent motifs and no subsurface refill refer to the original 45-state model. The new prototype does not resolve the energetics or experimental-validation gaps.

Current assessment: [multi-source scientific critique and calculated priorities](SCIENTIFIC_REVIEW.md). The current species model has no spatial interactions or subsurface refill; candidate enumeration does not establish completeness.

A **limited HiPRGen function-level pilot has now been run**. It executes the upstream composition-bucketing and decision-tree functions on a supplied library of capped Si-N/Si-O molecular motifs. It is not the full MPI/species-filtering/thermochemical pipeline and does not verify barriers, stable surface structures or mechanism completeness.

![HiPRGen candidate network](reaction_network/hiprgen_candidates.png)

## What was actually executed

The pinned upstream revision is [`a0dddfedc21be0121745e5f33f27ad8aafe796ea`](https://github.com/BlauGroup/HiPRGen/tree/a0dddfedc21be0121745e5f33f27ad8aafe796ea). Unmodified `bucket()` and `run_decision_tree()` definitions are loaded from the [licensed source snapshot](../third_party/hiprgen_snapshot/README.md). This avoids installing unrelated MPI/OpenBabel/reporting dependencies for a small local pilot. A compatible terminal enum and project-specific filters are supplied; this execution boundary is recorded in the output.

Each of four families contains **15 Si-centered motif identities plus HX and its coproduct**: NH2/F, NH2/Cl, OH/F, and OH/Cl. These are 68 family-specific records, with duplicates across families; not 68 unique chemical species. The supplied molecules have form `SiH_a X_b L_c`, with `a+b+c=4`, `L=NH2 or OH`, and `X=F or Cl`. H caps remain fixed during a reaction. The prospective step replaces one Si-L bond with Si-X and transfers the HX proton to the departing ligand:

`SiH_a X_b L_c + HX -> SiH_a X_(b+1) L_(c-1) + HL`.

Here `HL` is NH3 or H2O. This is a capped molecular proxy, not a periodic surface formula. In particular, NH2/OH ligands do not represent every bridging N/O connectivity in a real film.

| Family | Supplied records | Composition buckets | Directed pairs examined | Retained directions |
|---|---:|---:|---:|---:|
| NH2 / F | 17 | 85 | 296 | 20 |
| NH2 / Cl | 17 | 85 | 296 | 20 |
| OH / F | 17 | 85 | 296 | 20 |
| OH / Cl | 17 | 85 | 296 | 20 |

The result is **40 forward substitutions and 40 reverse candidates**. Filters retain one Si center, unchanged caps and one ligand substitution, and reject trivial/shared-species reactions. Atomic composition is checked by native HiPRGen bucketing. No free energies, synthetic barriers or rates are assigned. Molecular geometry generation, graph-isomorphism filtering, spin-state validation and TS searches are not performed by this pilot.

[Full species, reactions, filter counts and source hashes](../data/reaction_network/hiprgen/network.json) | [SQLite buckets](../data/reaction_network/hiprgen/) | [Runner](../scripts/hiprgen_intermediate_audit.py).

## Do we have enough intermediates?

**Not for a complete predictive dry-etch/ALE mechanism.** We have a bounded candidate library and a focused conditional kinetics model. The [intermediate gap table](../data/reaction_network/intermediate_gaps.csv) names the missing structural and kinetic evidence. A large node count is not evidence of completeness.

The highest-priority additions are explicit HF precursor complexes and protonation states; successive fluorination with remaining backbonds retained; separate final-release states for bridging NH, bare N and terminal NH2; hydrogenated Si product release; coadsorbed HF/H2O; and retained NH3/fluoride salts. Chlorine, fluorocarbon, defective/amorphous surfaces and ion-assisted events require separate branches and matching data.

An intermediate needs a chemically specified geometry/charge/spin and a stable minimum. A connected reaction needs a converged path, saddle curvature and endpoint connectivity. A usable thermal rate additionally needs a consistent energy reference and a justified prefactor/free-energy treatment. Surface site/anchor balance and gas reservoirs must remain explicit. Species with the same formula but different site, orientation, termination or bonding cannot be merged merely because their compositions match.

A network should be expanded until product branching and observable predictions are insensitive to plausible omitted pathways within the intended experimental domain, and then tested against matched measurements. Unknown barriers must not be replaced with arbitrary constants and labeled as verified kinetics. HiPRGen begins from the species library supplied to it; it cannot discover missing species outside that library or prove that all important intermediates were included. See the [HiPRGen method](https://doi.org/10.1039/D2DD00117A).

## Connection to the new kMC model

A separate [species-resolved HF/SiN:H model](species_kmc_results/REPORT.md) now has **45 named states and 55 enabled events**, representing all 16 published source pathways. It includes explicit HF complexes, successive fluorination states, surface nitrogen/hydrogen bookkeeping, and named gas products. Two future removal states remain unreachable because matching barriers are missing. Connecting source motifs into chains, arrival/desorption rates, prefactors and initial motif populations are assumptions; this is a conditional sensitivity calculation, not calibrated ALE.

The HiPRGen molecular candidate library is **not automatically imported as rates**. The new Python/C++ engine supports count-based single-motif events with external gas reservoirs. General bimolecular surface reactions, lateral diffusion and ion-impact propensities would require additional event handling and data. Candidate networks and enabled kinetics are plotted separately to make this boundary reviewable.

Reproduce:

```sh
python scripts/hiprgen_intermediate_audit.py --output outputs/hiprgen_reproduction
python scripts/build_species_cpp.py
python scripts/run_species_kmc.py
python scripts/plot_reaction_networks.py
```

The HiPRGen audit writes fresh SQLite buckets; use a new output directory for regeneration. Upstream third-party copyrights are retained; project commits remain authored by Ryan.
