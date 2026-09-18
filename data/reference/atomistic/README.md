# Attributed public atomistic subset

Sources: [Oh et al., Zenodo 19491140](https://zenodo.org/records/19491140) and [Martinez i Diaz et al., Zenodo 10211009](https://zenodo.org/records/10211009), both CC BY 4.0. Full creator lists, source checksums, selection rules and output hashes are in `provenance.json`.

`dft_subset.extxyz` preserves 70 original frame blocks containing DFT energy/force labels. Twenty-five are quasi-static drag configurations, not verified transition states. The source etch configurations were ML-generated and subsequently DFT-labeled, not AIMD trajectories. `hcl_geometries.extxyz` is an ASE conversion of 22 VASP geometries, including 15 source-named transition states; source energies, forces and Hessians are unavailable in the CONTCAR files. Index JSON files retain source identities.

These datasets retain upstream attribution independently of the project's sole-contributor authorship. The source material is licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). The subset selection and format conversion are the changes made here; no labels were calculated or fitted during curation.

See [acquisition and interpretation](../../../docs/ATOMISTIC_DATA.md). Reproduce with `scripts/collect_atomistic.py` from the repository root.
