# SiO2 / HF: four cluster steps and three slab pathways

Sources: [Kim et al., ACS Omega (2021)](https://doi.org/10.1021/acsomega.1c01824), Figure 2; [supplement, Figure S1](https://doi.org/10.1021/acsomega.1c01824.s001).

The CSV transcribes printed numbers, not a digitized continuous curve. The retrieved Figure 2 and Figure S1 tables do not explicitly label energy units; the surrounding ReaxFF energy convention suggests kcal/mol, but this is **unconfirmed**. Values are therefore retained as reported, excluded from the eV barrier comparison and not used for kinetics. Obtain explicit author confirmation before conversion. Slab formulas are surface-motif notation with implicit substrate atoms, not free-molecule equations.

DFT path method: VASP PAW/PBE, 400 eV, NEB followed by climbing-image NEB and dimer refinement when needed. This differs from the ADF/PBE-D3/TZ2P bond scans discussed elsewhere in the same paper. The paper includes stationary-point pictures; the retrieved SI supplies ReaxFF parameters but no Cartesian path images or raw NEB energies. No atomic trajectory is reconstructed from pictures.

The four cluster steps successively fluorinate a surface proxy and finally release SiF4. The three slab paths distinguish first Si-F/O-H formation, second fluorination on the same Si, and fluorination on a neighboring Si with water release. The mismatch between DFT and ReaxFF, especially cluster step c, motivates independent path validation before rate fitting.

Source SI is CC BY-NC 4.0. Raw source figures/force-field tables are not bundled into the code license. Source PDFs and figures can be reacquired using their publisher/Europe PMC links.
