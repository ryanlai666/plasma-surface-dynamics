"""Plot candidate HiPRGen ladders and the actual enabled/disabled kMC network."""

from pathlib import Path
import json
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

ROOT = Path(__file__).resolve().parents[1]


def arrow(ax, a, b, color='#52687b', style='-', rad=0):
    ax.add_patch(
        FancyArrowPatch(
            a,
            b,
            arrowstyle='-|>',
            mutation_scale=9,
            linewidth=1.1,
            color=color,
            linestyle=style,
            connectionstyle=f'arc3,rad={rad}',
            shrinkA=16,
            shrinkB=16,
        )
    )


def main():
    d = json.loads((ROOT / 'data/reaction_network/hiprgen/network.json').read_text())
    out = ROOT / 'docs/reaction_network'
    out.mkdir(exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(15, 10), layout='constrained')
    for ax, family in zip(axes.flat, [f['family'] for f in d['families']]):
        states = [
            s
            for s in d['species']
            if s['family'] == family and s['phase'] == 'capped_molecular_proxy'
        ]
        positions = {s['id']: (s['halogens'], -s['cap_H']) for s in states}
        for r in d['reactions']:
            if r['family'] != family or r['direction'] != 'forward_substitution':
                continue
            a = next(x for x in r['reactants'] if x in positions)
            b = next(x for x in r['products'] if x in positions)
            arrow(ax, positions[a], positions[b])
        for s in states:
            x, y = positions[s['id']]
            ax.scatter(
                x,
                y,
                s=1100,
                marker='s',
                c='#d5e9ef' if s['bound_ligands'] else '#dfe7c5',
                edgecolor='#455e70',
                zorder=3,
            )
            ax.text(
                x, y, s['id'].replace('(', '\n('), ha='center', va='center', fontsize=7, zorder=4
            )
        ax.set(
            xlim=(-0.6, 4.6),
            ylim=(-4.6, 0.7),
            xticks=range(5),
            yticks=[0, -1, -2, -3, -4],
            yticklabels=[0, 1, 2, 3, 4],
            xlabel='Halogen ligands on central Si',
            ylabel='Fixed H caps',
            title=family
            + ': + H'
            + family.split('_')[1]
            + ' / coproduct '
            + ('NH3' if family.startswith('NH') else 'H2O'),
        )
        ax.grid(alpha=0.1)
    fig.suptitle(
        'HiPRGen function-level pilot: 40 forward candidate substitutions + 40 reverse\nSupplied capped-motif library; no energy screening, surface stability or TS verification',
        fontsize=13,
    )
    fig.savefig(out / 'hiprgen_candidates.png', dpi=170)
    fig.savefig(out / 'hiprgen_candidates.svg')
    plt.close(fig)
    from draw_species_network import main as draw_species

    draw_species()
    print('Saved reaction network PNG/SVG plots')


if __name__ == '__main__':
    main()
