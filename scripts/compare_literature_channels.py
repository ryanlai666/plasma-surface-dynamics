"""Evaluate published per-occupied-state rate laws; no etch/selectivity fitting."""

from pathlib import Path
import csv, json, math, hashlib

ROOT = Path(__file__).resolve().parents[1]


def main():
    p = ROOT / 'data/literature/jung2020_coadsorption_rates.json'
    d = json.loads(p.read_text())
    rows = []
    for t in [300, 350, 400, 450]:
        for r in d['rows']:
            rows.append(
                dict(
                    **r,
                    temperature_K=t,
                    barrier_eV=r['Ea_over_R_K'] * 8.617333262145e-5,
                    rate_s=r['A_s'] * math.exp(-r['Ea_over_R_K'] / t),
                )
            )
    out = ROOT / 'docs/species_kmc_results'
    with (out / 'literature_channel_rates.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=rows[0])
        w.writeheader()
        w.writerows(rows)
    (out / 'literature_channel_rates.json').write_text(
        json.dumps(
            dict(
                source_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
                runner='scripts/compare_literature_channels.py',
                runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                rows=rows,
                interpretation='Rates per occupied coadsorbate state; not net etch rate or material selectivity. Flux, coverage, desorption and retained-product kinetics are required.',
            ),
            indent=2,
        )
        + '\n'
    )
    print([r for r in rows if r['temperature_K'] == 400])


if __name__ == '__main__':
    main()
