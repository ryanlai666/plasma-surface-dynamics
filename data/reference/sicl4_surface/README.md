# Published SiCl4 / Si(100) stationary structures

Source: Jianxun Zhang, Quan Zhu and Jun Li, [Symmetry 15, 213 (2023)](https://doi.org/10.3390/sym15010213), supporting information. **CC BY 4.0**; credit remains with the source authors.

Seventeen coordinate tables were converted to ASE extended XYZ, preserving their coordinates and atom order: clean surface, isolated SiCl4, four physisorbed complexes, four dissociated products, one intermediate and six source transition states. Each adsorbate/slab structure has 117 atoms (81 Si, 32 H, 4 Cl). The source DFT energies are in [the pathway table](../../literature/sicl4_surface_paths.json), not invented calculator labels in these XYZ files.

The coordinate tables do not provide cell vectors, raw forces, Hessians or all NEB images. No periodic cell was guessed. `pbc=False` in the conversion means absent cell metadata, not that the original calculation was nonperiodic. Recover the original cell and constraints before a periodic recalculation.

[Provenance](provenance.json) includes source SHA256 hashes; [inventory](inventory.json) has hashes of converted files. Regenerate with `python scripts/collect_surface_paths.py --archive PATH_TO_PUBLISHER_ZIP` from the repository root, using an environment with ASE and PyMuPDF. Omitting the argument downloads the public 27.8 MB supplement into ignored `data/raw/`.
