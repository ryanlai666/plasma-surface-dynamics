"""Transient endpoint sensitivities and mixture identifiability for the conditional model."""

from pathlib import Path
from copy import deepcopy
import sys, csv, json, hashlib
import numpy as np
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from plasma_surface.species_kmc import build_network, rate_constants, expectation, KB


def main():
    n = build_network()
    rates = rate_constants(n)
    removed = np.array([s['Si_removed'] for s in n['states']])
    times = np.array([0.0, 2.0, 3.0])
    sites = 1000

    def outcome(net, r):
        c, g = expectation(net, r, sites, times)
        return float(c[-1, removed].sum() / sites)

    base = outcome(n, rates)
    rows = []
    groups = sorted({e['pathway'] for e in n['events'] if e['pathway']}) + [
        'HF_arrival',
        'HF_desorption',
    ]
    for name in groups:
        mask = np.array(
            [
                e['pathway'] == name
                if name not in ['HF_arrival', 'HF_desorption']
                else e['kind'] == ('adsorption' if name == 'HF_arrival' else 'desorption')
                for e in n['events']
            ]
        )
        values = []
        for factor in [1 / 1.1, 1.1]:
            r = rates.copy()
            r[mask] *= factor
            values.append(outcome(n, r))
        sensitivity = float((np.log(values[1]) - np.log(values[0])) / (2 * np.log(1.1)))
        rows.append(
            dict(
                group=name,
                relative_step=1.1,
                baseline_yield=base,
                lower_rate_yield=values[0],
                higher_rate_yield=values[1],
                log_yield_sensitivity=sensitivity,
                event_ids=[e['id'] for e, m in zip(n['events'], mask) if m],
            )
        )
    pure = []
    for initial, weight in n['initial_motif_fractions'].items():
        nn = deepcopy(n)
        nn['initial_motif_fractions'] = {initial: 1.0}
        pure.append(
            dict(
                initial_state=initial, assumed_fraction=weight, pure_motif_yield=outcome(nn, rates)
            )
        )
    reconstructed = sum(r['assumed_fraction'] * r['pure_motif_yield'] for r in pure)
    assert abs(reconstructed - base) < 1e-7
    out = ROOT / 'docs/species_kmc_results'
    d = dict(
        temperature_K=400,
        exposure_s=2,
        purge_s=1,
        observable='Si released per initial motif at 3 s; not steady-state rate',
        baseline_yield=base,
        sensitivities=rows,
        pure_motif_results=pure,
        mixture_reconstruction=reconstructed,
        mixture_linearity_error=abs(reconstructed - base),
        barrier_rate_factor_for_0p1_eV=float(np.exp(0.1 / (KB * 400))),
        network_sha256=hashlib.sha256(
            (ROOT / 'configs/species_kmc_network.json').read_bytes()
        ).hexdigest(),
        runner='scripts/analyze_kinetic_priorities.py',
        runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        interpretation='One-at-a-time transient sensitivities, not global uncertainty intervals or equilibrium degree of rate control. Forward rates changed without reverse reactions in this irreversible model.',
    )
    (out / 'kinetic_priorities.json').write_text(json.dumps(d, indent=2) + '\n')
    with (out / 'kinetic_priorities.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=rows[0])
        w.writeheader()
        w.writerows(rows)
    fig, axes = plt.subplots(1, 2, figsize=(12, 6), layout='constrained')
    order = sorted(rows, key=lambda r: abs(r['log_yield_sensitivity']), reverse=True)[:10]
    axes[0].barh(
        [r['group'] for r in order][::-1],
        [r['log_yield_sensitivity'] for r in order][::-1],
        color=['#C9754B' if r['log_yield_sensitivity'] < 0 else '#2B9B9E' for r in order][::-1],
    )
    axes[0].axvline(0, color='gray', lw=0.6)
    axes[0].set(
        title='What controls the 3 s endpoint?', xlabel='d ln(Si-release fraction) / d ln(rate)'
    )
    axes[1].barh(
        [r['initial_state'].replace('_F0', '').replace('_', ' ') for r in pure][::-1],
        [r['pure_motif_yield'] for r in pure][::-1],
        color='#A489CF',
    )
    axes[1].set(
        title='Same rates; different initial motifs',
        xlabel='Si release fraction for a pure motif inventory',
        xlim=(0, 1),
    )
    fig.suptitle('400 K conditional diagnostics: sensitivity and mixture confounding', fontsize=14)
    fig.savefig(out / 'kinetic_priorities.png', dpi=170)
    plt.close(fig)
    print(json.dumps(dict(baseline=base, top=order[:5], pure=pure)), flush=True)


if __name__ == '__main__':
    main()
