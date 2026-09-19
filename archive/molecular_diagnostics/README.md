# Archived molecular workflow diagnostic

The NH3 umbrella-inversion calculations were moved here because they do not validate dry etching, Si-N bond cleavage or ALE parameters. They are excluded from the current project results. Their raw outputs and source hashes are preserved for reproducibility.

- [Historical calculation report](barrier_results/REPORT.md)
- [Original calculation script](molecular_barrier.py)
- [Original report generator, preserved as text](report_barriers_legacy.py.txt)
- [Isolated PySCF requirements](requirements-dft.txt)

The original report generator is archived as text and must not be used to update the main README. Its source layout references describe the previous repository layout. Run the molecular calculation script only with a new explicit output directory if this diagnostic is needed separately. The main research workflow is documented in [dry-etch parameters](../../docs/DRY_ETCH_PARAMETERS.md).
