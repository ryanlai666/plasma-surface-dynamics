"""Evaluate completed local-motif path files with DFT as they become available."""

from pathlib import Path
import json, subprocess, sys, time

ROOT = Path(__file__).resolve().parents[1]
pending = ['SiN_capped_motif', 'SiO_capped_motif']
deadline = time.monotonic() + 2400
while pending:
    for motif in list(pending):
        parent = ROOT / 'data/surface_paths' / motif / 'HF'
        source = parent / 'omol25_refined'
        if not (source / 'summary.json').exists():
            continue
        d = json.loads((source / 'summary.json').read_text())
        if d['status'] in ['neb_converged', 'neb_unconverged']:
            subprocess.run(
                [
                    sys.executable,
                    str(ROOT / 'scripts/dft_molecular_singlepoints.py'),
                    str(source / 'images.extxyz'),
                    str(parent / 'dft_path'),
                ],
                check=True,
            )
            print('DFT complete', motif, 'ML geometry status', d['status'], flush=True)
            pending.remove(motif)
        elif d['status'] in ['endpoint_validation_failed', 'endpoints_collapsed']:
            raise RuntimeError('No reaction path available for ' + motif)
    if time.monotonic() > deadline:
        raise TimeoutError('ML path jobs did not finish within 40 minutes')
    if pending:
        time.sleep(10)
