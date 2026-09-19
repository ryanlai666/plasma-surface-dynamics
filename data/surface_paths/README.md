# Surface / reactant path collection

| Surface | Species | Contents |
|---|---|---|
| [Si(100)](Si100/F/mace_neb/README.md) | F | Calculated rigid-surface migration path, energies, GIF and checks |
| [Si(100)](Si100/Cl/mace_neb/README.md) | Cl | Same protocol, different halogen |
| [Si(111)](Si111/F/mace_neb/README.md) | F | Different facet, same halogen |
| [Si(100)-c(4x2)](Si100_c4x2/SiCl4/published/README.md) | SiCl4 | 17 published geometries; six three-point slideshows |
| [Fluorinated Si3N4](Si3N4/HF_H2O/literature/README.md) | HF + H2O | Published R6 activation and prefactor |
| [Fluorinated SiO2](SiO2/HF_H2O/literature/README.md) | HF + H2O | Published R6 comparison |
| [Fluorinated Si3N4](Si3N4/HF_HFv1/literature/README.md) | HF + HF(v=1) | Source excitation-assisted R7 model |
| [Fluorinated SiO2](SiO2/HF_HFv1/literature/README.md) | HF + HF(v=1) | Source quenching/R7 comparison |
| [Amorphous SiNx:H](SiNxH_amorphous/HF/literature/README.md) | HF | 16 surface-motif pathways; unresolved energy zero |
| [SiO2 cluster/slab](SiO2/HF/literature/README.md) | HF | Seven DFT/ReaxFF reaction records; unit confirmation pending |
| [Reconstructed Si(100)](Si_r100/F2/literature/reference.json) | F2 | Published first-event rate |
| [Unreconstructed Si(100)](Si_u100/F2/literature/reference.json) | F2 | Published first-event rate |
| [Si(110)](Si_110/F2/literature/reference.json) | F2 | Published first-event rate |
| [Si(111)](Si_111/F2/literature/reference.json) | F2 | Published first-event rate |

`mace_neb` means a locally calculated ML path; `published` means source coordinates; `literature` contains cited parameters without a recovered Cartesian path. New methods belong in separate subdirectories so their energy zeros and provenance cannot be mixed.

[Comparison, equations and limitations](../../docs/dry_etch_results/SURFACE_PATHS.md).
