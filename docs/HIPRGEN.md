# HiPRGen assessment for Si / SiNx plasma-surface reactions

## Recommendation

Use HiPRGen as an optional upstream **candidate-reaction generator**, followed by chemistry-specific filtering and rate validation. It is not a drop-in plasma-surface kinetics backend. The project has not installed or run HiPRGen, and no HiPRGen-generated rates have been used in the reported simulations.

The current upstream repository is [BlauGroup/HiPRGen](https://github.com/BlauGroup/HiPRGen). Its workflow enumerates composition-balanced molecular reactions, applies customizable filters, and uses MPI for network generation. RNMC is a separate network-simulation component. These capabilities make it a plausible research extension for a curated, small species set, rather than a source of ready-made Si/SiNx reaction data.

## What needs adaptation

The [method paper](https://pubs.rsc.org/en/content/articlepdf/2023/dd/d2dd00117a) starts from species with precomputed properties; it does not generate the required species dataset from scratch. Its demonstrated chemistry concerns electrochemical networks. For our application, prepare consistent Si/N/H/F/Cl molecular and surface-fragment calculations first. The downloaded HCl geometries alone do not provide those energies or a complete species library.

Inspection of [reaction_questions.py](https://github.com/BlauGroup/HiPRGen/blob/main/HiPRGen/reaction_questions.py) shows a default transition-state-like rate expression and a filter that can construct a barrier from reaction free energy plus a constant. Such assigned barriers are not measured or DFT transition-state barriers. The code also contains electrochemical electron-free-energy and solvent-related machinery. These are not substitutes for a nonthermal electron-energy distribution, ion impact, or surface charging model.

Our proposed adaptations are:

1. Represent gas molecules, radicals, surface-bound fragments, vacancies, and site types explicitly. Preserve substrate anchors and distinguish Si sites from N sites. A periodic slab cannot simply be treated as a free gas molecule.
2. Enforce elemental and site balance. Track electrons/ions or reservoirs explicitly whenever charge changes; document external energy delivery.
3. Replace battery-specific species/reaction filters. In particular, do not reject all uphill reactions merely because thermal free energy is positive: plasma excitation or impact may drive them. Conversely, a favorable free energy does not prove kinetic accessibility.
4. Obtain thermal barriers from compatible calculations or measurements. Parameterize ion-driven probabilities versus energy, angle, and local coverage from suitable MD/beam data. Derive electron-impact rates from collision data and the electron distribution.
5. Retain candidate provenance, unknown barriers, uncertainty, and applicability ranges. Reject conversion into executable kinetics when rate units, surface context, or required barriers are missing.

These are proposed design requirements, not functionality currently implemented in HiPRGen or this project.

## Proposed connection boundary

The inspected [reaction database schema](https://github.com/BlauGroup/HiPRGen/blob/main/HiPRGen/reaction_filter.py) exposes reactant/product IDs, rate, reaction free energy, barrier, and redox status. An adapter can read that SQLite database without loading molecule pickle files, then join an independently exported species table with composition, charge, phase, and surface-site metadata. Preserve upstream rates as unvalidated annotations until their physical meaning and units are established.

Suggested intermediate record:

```json
{
  "source_reaction_id": 0,
  "upstream_commit": "record-the-reviewed-commit",
  "reactants": [{"species_id": "surface_SiCl", "count": 1}],
  "products": [],
  "surface_context": "must-be-specified",
  "rate_law": null,
  "rate_units": null,
  "barrier_ev": null,
  "status": "candidate_not_simulatable"
}
```

This incomplete record illustrates the metadata contract only; it is not a balanced reaction or input accepted by the current solver. A future adapter must require balanced complete reactions before rate construction.

The current C++ solver has two site states and five fixed event classes. A general HiPRGen network needs a species-resolved propensity engine and Si/N inventory bookkeeping, or a documented reduction into those five event classes. Simply importing a reaction list would not create a correct coupled simulation.

## Local-computer pilot

Start with roughly 20-50 curated species/fragments for **one** chosen chemistry, not all Si/N/H/F/Cl combinations at once. Pin the upstream commit, review its license/dependencies, and use an isolated environment; test the upstream examples before changing filters. Prefer a small MPI job and inspect the resulting candidate count before expanding. Validate several reaction paths against independent barriers or published mechanisms before enabling any rates in the ALE solver. This pilot can establish whether the added network complexity improves the scientific question enough to justify further work.

Assessment based on the linked upstream README, database schema, reaction filters, and method paper inspected during this session. No third-party source code was copied into this repository.
