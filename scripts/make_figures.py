"""Draw the presentation figures in docs/figures/ from saved results only.

This script runs no simulation and changes no result.  It reads saved
trajectories, summaries and energy tables, and writes redesigned figures with a
manifest recording the SHA-256 of every input.  The original diagnostic plots
written by the calculation scripts are left untouched.

    python scripts/make_figures.py            # all figures
    python scripts/make_figures.py kinetics   # one figure by name
"""

import csv
import hashlib
import io
import json
import sys
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from plasma_surface import figstyle as fs  # noqa: E402

OUT = ROOT / 'docs/figures'
INPUTS: set[str] = set()

FAMILY = {
    'bridge_NH': 'NH bridge',
    'terminal_NH2': 'terminal NH₂',
    'Si_H_rich': 'Si–H-rich',
    'bridge_N': 'bare N bridge',
    'pre_SiH2F': 'pre-existing SiH₂F',
    'pre_SiHF2': 'pre-existing SiHF₂',
    'Si_Si_backbond': 'Si–Si backbond',
}
SURFACE = {
    'Si100': 'Si(100)',
    'Si111': 'Si(111)',
    'beta_Si3N4_001': 'β-Si₃N₄(001)',
    'alpha_quartz_001': 'α-quartz(001)',
}
GAS = {'NH3': 'NH₃', 'SiF4': 'SiF₄', 'SiH2F2': 'SiH₂F₂', 'SiHF3': 'SiHF₃', 'H2': 'H₂'}


def use(path) -> Path:
    """Register a saved input so its hash lands in the manifest."""
    path = ROOT / path
    INPUTS.add(path.relative_to(ROOT).as_posix())
    return path


def load_json(path):
    return json.loads(use(path).read_text(encoding='utf-8'))


def load_csv(path):
    with use(path).open(encoding='utf-8') as stream:
        return list(csv.DictReader(stream))


def save(fig, name):
    fig.savefig(OUT / name)
    plt.close(fig)


def network_states():
    from plasma_surface.species_kmc import build_network

    use('data/literature/sin_hf_pathways.csv')
    return build_network()


# --------------------------------------------------------------------------- 1
def kinetics():
    """Si release vs temperature, and gas products at 400 K."""
    summary = load_json('docs/species_kmc_results/summary.json')
    rows = load_csv('docs/species_kmc_results/species_populations.csv')
    dose, purge = summary['exposure_s'], summary['purge_s']
    si_gases = ('SiF4', 'SiH2F2', 'SiHF3')

    def series(temperature, species):
        pts = sorted(
            (float(r['time_s']), float(r['fraction_per_initial_motif']))
            for r in rows
            if r['phase'] == 'gas_product'
            and float(r['temperature_K']) == temperature
            and r['species'] in species
        )
        times = sorted({t for t, _ in pts})
        return np.array(times), np.array([sum(v for t2, v in pts if t2 == t) for t in times])

    fig, (a, b) = plt.subplots(1, 2, figsize=(12, 5.2), gridspec_kw=dict(wspace=0.42))
    fig.subplots_adjust(top=0.78, bottom=0.14, left=0.07, right=0.9)
    temps = sorted(t['temperature_K'] for t in summary['temperatures'])
    ends, labels = [], []
    for temp, color in zip(temps, fs.RAMP4):
        t, y = series(temp, si_gases)
        a.plot(t, 100 * y, color=color)
        ends.append(100 * y[-1])
        labels.append(
            f'{temp:.0f} K  {100 * y[-1]:.0f}%' if y[-1] >= 0.01 else f'{temp:.0f} K  <1%'
        )
    fs.shade_phases(a, [(0, dose, 'HF dose'), (dose, dose + purge, 'purge')])
    fs.end_labels(a, dose + purge, ends, labels, fs.RAMP4, min_gap=4.5, dx=0.12)
    a.set(
        xlim=(0, dose + purge),
        ylim=(0, 55),
        xlabel='Time (s)',
        ylabel='Motifs that released Si (%)',
    )
    a.set_title('Si released, by temperature', pad=22)
    a.legend(
        handles=[Line2D([], [], color=c, label=f'{t:.0f} K') for t, c in zip(temps, fs.RAMP4)],
        loc='upper left',
        ncol=2,
        fontsize=9,
    )

    gases = [('NH3', fs.BLUE), ('SiF4', fs.ORANGE), ('SiH2F2', fs.AQUA), ('SiHF3', fs.YELLOW)]
    ends, labels = [], []
    for gas, color in gases:
        t, y = series(400.0, (gas,))
        b.plot(t, 100 * y, color=color)
        ends.append(100 * y[-1])
        labels.append(f'{GAS[gas]}  {100 * y[-1]:.1f}')
    fs.shade_phases(b, [(0, dose, 'HF dose'), (dose, dose + purge, 'purge')])
    fs.end_labels(b, dose + purge, ends, labels, [c for _, c in gases], min_gap=2.6, dx=0.12)
    b.set(
        xlim=(0, dose + purge),
        ylim=(0, 40),
        xlabel='Time (s)',
        ylabel='Molecules per 100 starting motifs',
    )
    b.set_title('Gas products at 400 K', pad=22)
    b.legend(
        handles=[Line2D([], [], color=c, label=GAS[g]) for g, c in gases],
        loc='upper left',
        ncol=2,
        fontsize=9,
    )
    fs.header(
        fig,
        'Si release switches on between 350 K and 400 K',
        'Species-resolved kMC: 2 s HF dose, then 1 s purge. Published barriers with assumed prefactors, '
        'HF flux and starting surface; conditional sensitivity, not a calibrated etch rate.',
    )
    fs.footnote(
        fig,
        'H₂ output is zero at every temperature. Source: docs/species_kmc_results/species_populations.csv',
    )
    save(fig, 'si_release_by_temperature.png')


# --------------------------------------------------------------------------- 2
def priorities():
    """What controls Si release at 400 K, and which motifs can release it."""
    d = load_json('docs/species_kmc_results/kinetic_priorities.json')
    network = network_states()
    states = network['states']
    events = {e['id']: e for e in network['events']}

    def describe(group, ids):
        if group == 'HF_arrival':
            return 'HF arrival rate (assumed)'
        if group == 'HF_desorption':
            return 'HF desorption rate (assumed)'
        steps = []
        for i in ids:
            e = events[i]
            src, dst = states[e['source']], states[e['target']]
            step = (
                'Si release'
                if dst['Si_removed']
                else f"F{src['fluorination_stage']}→F{dst['fluorination_stage']}"
            )
            steps.append(f"{FAMILY[src['family']]}, {step}")
        return f'{group}: ' + '; '.join(dict.fromkeys(steps))

    rows = [
        (describe(s['group'], s['event_ids']), s['log_yield_sensitivity'])
        for s in d['sensitivities']
    ]
    shown = sorted([r for r in rows if abs(r[1]) >= 1e-3], key=lambda r: abs(r[1]))
    hidden = len(rows) - len(shown)

    fig, (a, b) = plt.subplots(
        1, 2, figsize=(13, 5.6), gridspec_kw=dict(width_ratios=[1.15, 1], wspace=0.75)
    )
    fig.subplots_adjust(top=0.78, bottom=0.24, left=0.24, right=0.95)
    y = np.arange(len(shown))
    vals = [v for _, v in shown]
    bars = a.barh(y, vals, height=0.62, color=[fs.BLUE if v > 0 else fs.RED for v in vals])
    a.set_yticks(y, [n for n, _ in shown])
    a.axvline(0, color=fs.AXIS, lw=0.8)
    a.grid(axis='x')
    a.grid(axis='y', visible=False)
    for bar, v in zip(bars, vals):
        a.annotate(
            f'{v:+.3f}',
            xy=(v, bar.get_y() + bar.get_height() / 2),
            xytext=(4 if v > 0 else -4, 0),
            textcoords='offset points',
            ha='left' if v > 0 else 'right',
            va='center',
            fontsize=9,
            color=fs.INK2,
        )
    a.set_xlim(-0.14, 0.3)
    a.set_xlabel('Change in Si release per change in rate (log–log slope)')
    a.set_title('Which rate matters most?')
    a.text(
        0,
        -0.2,
        f'Blue: faster rate → more Si released. Red: faster rate → less.\n'
        f'The other {hidden} reaction steps have |slope| < 0.001 and are not shown.',
        transform=a.transAxes,
        fontsize=8.5,
        color=fs.MUTED,
        va='top',
    )

    pure = d['pure_motif_results']
    names = [
        FAMILY[states[[s['id'] for s in states].index(p['initial_state'])]['family']] for p in pure
    ]
    yb = np.arange(len(pure))[::-1]
    bars = b.barh(yb, [100 * p['pure_motif_yield'] for p in pure], height=0.62, color=fs.BLUE)
    b.set_yticks(
        yb, [f"{n}\n{100 * p['assumed_fraction']:.0f}% of start" for n, p in zip(names, pure)]
    )
    b.tick_params(axis='y', labelsize=9)
    for bar, p in zip(bars, pure):
        v = 100 * p['pure_motif_yield']
        b.annotate(
            f'{v:.0f}%' if v >= 0.5 else '0% (no release path)',
            xy=(v, bar.get_y() + bar.get_height() / 2),
            xytext=(4, 0),
            textcoords='offset points',
            va='center',
            fontsize=9,
            color=fs.INK2,
        )
    b.set_xlim(0, 125)
    b.set_xticks([0, 25, 50, 75, 100])
    b.grid(axis='x')
    b.grid(axis='y', visible=False)
    b.set_xlabel('Si released if the surface held only this motif (%)')
    b.set_title('Which starting motifs can release Si?')
    fs.header(
        fig,
        'HF supply, not any single reaction barrier, limits Si release at 400 K',
        f"One-at-a-time rate changes (±{100 * (d['sensitivities'][0]['relative_step'] - 1):.0f}%) at "
        f"{d['temperature_K']:.0f} K after a {d['exposure_s']:.0f} s dose and {d['purge_s']:.0f} s purge. "
        'Only 3 of the 7 motif types can release Si at all, so the assumed starting mixture sets the ceiling.',
    )
    save(fig, 'what_controls_si_release.png')


# --------------------------------------------------------------------------- 3
def multilayer_stall():
    d = load_json('docs/multilayer_results/cycle_diagnosis.json')
    cycles = d['cycles']
    fig, (a, b) = plt.subplots(
        1, 2, figsize=(12.5, 5), gridspec_kw=dict(width_ratios=[1.7, 1], wspace=0.35)
    )
    fig.subplots_adjust(top=0.76, bottom=0.14, left=0.06, right=0.97)
    x = np.array([c['cycle'] for c in cycles])
    b5 = np.array([c['removed_by_band'].get('5', 0) for c in cycles])
    b4 = np.array([c['removed_by_band'].get('4', 0) for c in cycles])
    other = np.array(
        [sum(v for k, v in c['removed_by_band'].items() if k not in ('4', '5')) for c in cycles]
    )
    assert not other.any(), 'only the top two bands are expected to lose atoms'
    a.bar(x, b5, width=0.62, color=fs.BLUE, label='Top band (B5)')
    a.bar(x, b4, width=0.62, bottom=b5, color=fs.ORANGE, label='Next band down (B4)')
    for xi, total in zip(x, b5 + b4):
        a.annotate(
            str(total),
            xy=(xi, total),
            xytext=(0, 3),
            textcoords='offset points',
            ha='center',
            fontsize=9,
            color=fs.INK2,
        )
    a.annotate(
        'No atoms removed\nin cycles 3–12',
        xy=(7.5, 1.2),
        ha='center',
        va='bottom',
        fontsize=10,
        color=fs.INK2,
    )
    a.set_xticks(x)
    a.set(
        xlabel='Dose/purge cycle (identical rates every cycle)',
        ylabel='Substrate atoms removed',
        ylim=(0, 25),
    )
    a.set_title('Atoms removed per cycle')
    a.legend(loc='upper right')

    blocked = Counter(
        'bare N' if r['partner_H'] == 0 else f"NH{r['partner_H']}"
        for r in d['blocked_final_cleavages']
    )
    labels = {'bare N': 'Bonded to bare N', 'NH2': 'Bonded to NH₂'}
    keys = sorted(blocked, key=blocked.get)
    bars = b.barh(range(len(keys)), [blocked[k] for k in keys], height=0.55, color=fs.BLUE)
    b.set_yticks(range(len(keys)), [labels.get(k, k) for k in keys])
    fs.bar_value_labels(b, bars)
    b.grid(axis='x')
    b.grid(axis='y', visible=False)
    b.set(xlabel='Stuck SiF₃ sites', xlim=(0, max(blocked.values()) * 1.25))
    b.set_title('Why removal stops')
    b.text(
        0,
        -0.3,
        'Each stuck Si has three F caps and one last Si–N bond.\n'
        'No barrier is available for breaking that bond,\nso the step is disabled, not slow.',
        transform=b.transAxes,
        fontsize=9,
        color=fs.INK2,
        va='top',
    )
    fs.header(
        fig,
        'Multilayer etching stalls after cycle 2 because a rate is missing',
        f"{sum(blocked.values())} reachable Si sites wait on a final Si–N cleavage step with no rate. "
        'This is a gap in the mechanism, not evidence of self-limiting ALE.',
    )
    fig.subplots_adjust(bottom=0.25)
    save(fig, 'multilayer_stall.png')


# --------------------------------------------------------------------------- 4
def relaxation():
    rows = load_json('docs/dry_etch_results/relaxed_adsorption.json')
    surfaces = list(dict.fromkeys(r['surface'] for r in rows))
    fig, axes = plt.subplots(1, 4, figsize=(14, 4.6), gridspec_kw=dict(wspace=0.32))
    fig.subplots_adjust(top=0.74, bottom=0.2, left=0.05, right=0.97)
    drops = []
    for ax, surface in zip(axes, surfaces):
        ends, labels, colors = [], [], []
        for r, color in zip([r for r in rows if r['surface'] == surface], (fs.BLUE, fs.ORANGE)):
            data = load_csv(r['folder'].replace('\\', '/') + '/energies.csv')
            steps = np.array([int(x['step']) for x in data])
            e = np.array([float(x['energy_eV']) for x in data])
            e -= e[0]
            ax.plot(steps, e, color=color)
            ax.plot(
                steps[-1],
                e[-1],
                'o',
                ms=7,
                color=color if r['converged'] else fs.SURFACE,
                mec=color,
                mew=2,
                zorder=4,
            )
            ends.append(e[-1])
            colors.append(color)
            labels.append(
                f"Start {r['start']}: {e[-1]:+.1f} eV"
                + ('' if r['converged'] else ' (not converged)')
            )
            drops.append(-e[-1])
        span = max(0.4, -min(ends))
        ax.set_ylim(-span * 1.12, span * 0.08)
        # Final values in the empty upper-right corner: coloured key, ink text.
        for row, (label, color) in enumerate(zip(labels, colors)):
            y = 0.9 - 0.09 * row
            ax.text(
                0.97,
                y,
                label,
                transform=ax.transAxes,
                ha='right',
                va='center',
                fontsize=9,
                color=fs.INK,
            )
            ax.plot([0.985, 1.02], [y, y], transform=ax.transAxes, color=color, lw=3, clip_on=False)
        ax.set_title(SURFACE.get(surface, surface))
        ax.set_xlabel('Optimizer step')
    axes[0].set_ylabel('Energy change from start (eV)')
    fig.legend(
        handles=[
            Line2D([], [], color=fs.BLUE, label='Start 1'),
            Line2D([], [], color=fs.ORANGE, label='Start 2'),
            Line2D(
                [],
                [],
                ls='',
                marker='o',
                color=fs.INK2,
                ms=7,
                label='converged (force < 0.04 eV/Å)',
            ),
            Line2D(
                [],
                [],
                ls='',
                marker='o',
                mfc=fs.SURFACE,
                mec=fs.INK2,
                mew=2,
                ms=7,
                label='not converged',
            ),
        ],
        loc='lower left',
        ncol=4,
        bbox_to_anchor=(0.04, 0.0),
    )
    converged = sum(r['converged'] for r in rows)
    fs.header(
        fig,
        f'HF plus the top of the slab relax downhill by {min(drops):.1f}–{max(drops):.1f} eV',
        f'MACE relaxation from two screened starting sites per surface; {converged} of {len(rows)} runs converge. '
        'The drop includes substrate reconstruction, so it is not an adsorption energy or a barrier.',
    )
    save(fig, 'hf_surface_relaxation.png')


# --------------------------------------------------------------------------- 5
def final_cleavage():
    base = 'data/final_cleavage/SiF3_NH2_HF'
    d = load_json(f'{base}/dft_initial_path/summary.json')
    scan = load_json(f'{base}/fixed_frame_three_coordinate_scan/summary.json')
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(15, 5), gridspec_kw=dict(wspace=0.3))
    fig.subplots_adjust(top=0.74, bottom=0.2, left=0.05, right=0.98)
    geoms = ['Reactant', 'Highest point\nof failed path', 'Product']
    methods = [('OMol25 (ML)', d['results'][0]['relative_ML_energies_eV'], fs.BLUE)]
    colors = {'def2-svp': fs.ORANGE, 'def2-tzvp': fs.AQUA}
    for r in d['results']:
        methods.append((f"DFT PBE/{r['basis']}", r['relative_DFT_energies_eV'], colors[r['basis']]))
    offsets = np.linspace(-0.18, 0.18, len(methods))
    for (label, energies, color), off in zip(methods, offsets):
        for i, e in enumerate(energies):
            a.hlines(
                e,
                i + off - 0.08,
                i + off + 0.08,
                color=color,
                lw=3,
                label=label if i == 0 else None,
            )
        a.plot(np.arange(3) + off, energies, color=color, lw=0.8, ls=':', alpha=0.8)
    peak = [m[1][1] for m in methods]
    a.annotate(
        f'{min(peak):.2f}–{max(peak):.2f} eV',
        xy=(1, max(peak)),
        xytext=(0, 8),
        textcoords='offset points',
        ha='center',
        fontsize=9.5,
        color=fs.INK2,
    )
    a.set_xticks(range(3), geoms)
    a.set(ylabel='Energy relative to reactant (eV)', xlim=(-0.5, 2.5), ylim=(-1, 3))
    a.set_title('Methods agree on the energies…')
    a.legend(loc='upper right', fontsize=9)

    x = np.arange(3)
    width = 0.32
    for k, r in enumerate(d['results']):
        forces = [f['max_mobile_DFT_force_eV_A'] for f in r['frames']]
        bars = b.bar(
            x + (k - 0.5) * width,
            forces,
            width * 0.9,
            color=colors[r['basis']],
            label=f"PBE/{r['basis']}",
        )
        for bar, f in zip(bars, forces):
            b.annotate(
                f'{f:.2f}',
                xy=(bar.get_x() + bar.get_width() / 2, f),
                xytext=(0, 3),
                textcoords='offset points',
                ha='center',
                fontsize=8.5,
                color=fs.INK2,
            )
    b.axhline(
        0.04, color=fs.INK, lw=1.2, zorder=5, label='< 0.04 eV/Å needed for a stationary point'
    )
    b.set_xticks(x, geoms)
    b.set(ylabel='Largest force on a mobile atom (eV/Å)', ylim=(0, 2.3))
    b.set_title('…but none of the points is stationary')
    b.legend(loc='upper left', fontsize=9)

    e0 = scan['results'][0]['energy_eV']
    for direction, color, label in (
        ('forward', fs.BLUE, 'Stretching the bond'),
        ('reverse', fs.ORANGE, 'Shortening it again'),
    ):
        rr = sorted(
            [r for r in scan['results'] if r['direction'] == direction], key=lambda r: r['SiN_A']
        )
        xs = [r['SiN_A'] for r in rr]
        ys = [r['energy_eV'] - e0 for r in rr]
        c.plot(xs, ys, color=color, label=label)
        ok = [r['converged'] for r in rr]
        c.plot(
            [v for v, o in zip(xs, ok) if o],
            [v for v, o in zip(ys, ok) if o],
            'o',
            color=color,
            ms=6,
        )
        c.plot(
            [v for v, o in zip(xs, ok) if not o],
            [v for v, o in zip(ys, ok) if not o],
            'o',
            mfc=fs.SURFACE,
            mec=color,
            mew=2,
            ms=7,
        )
    c.plot(
        [], [], 'o', mfc=fs.SURFACE, mec=fs.INK2, mew=2, ms=7, ls='', label='point not converged'
    )
    c.set(xlabel='Si–N distance held fixed (Å)', ylabel='Energy relative to start (eV)')
    c.set_title('Forward and reverse scans disagree')
    c.legend(loc='upper left', fontsize=9)
    fs.header(
        fig,
        'The final Si–N cleavage barrier is still unresolved',
        'Local SiF₃–NH₂ + HF molecular model. The DFT energies are real, but no point is a transition state, '
        'and the constrained scan shows hysteresis, so no rate enters the kMC.',
    )
    save(fig, 'final_cleavage_checks.png')


# ------------------------------------------------------------------ animations
def _gif(frames, name, ms=160, hold_ms=2200):
    durations = [ms] * len(frames)
    durations[-1] = hold_ms
    frames[0].save(
        OUT / name,
        save_all=True,
        append_images=frames[1:],
        duration=durations,
        loop=0,
        optimize=True,
    )
    return len(frames)


def _grab(fig):
    buffer = io.BytesIO()
    fig.savefig(buffer, format='png', dpi=100)
    plt.close(fig)
    buffer.seek(0)
    return Image.open(buffer).convert('P', palette=Image.ADAPTIVE, colors=128)


def species_animation():
    traj = np.load(use('docs/species_kmc_results/lattice_trajectory.npz'))
    summary = load_json('docs/species_kmc_results/summary.json')
    states = network_states()['states']
    dose, purge = summary['exposure_s'], summary['purge_s']
    removed = np.array([s['Si_removed'] for s in states])
    hf = np.array([s['HF_complex'] for s in states])
    stage = np.array([min(s['fluorination_stage'], 3) for s in states])
    times, grids, counts = traj['times'], traj['states'], traj['counts']
    n = grids.shape[1]
    side = int(round(np.sqrt(n)))
    fractions = np.stack(
        [counts[:, (~removed) & (stage == k)].sum(1) / n for k in range(4)]
        + [counts[:, removed].sum(1) / n],
        axis=1,
    )
    colors = fs.RAMP4 + [fs.GRID]
    names = ['F0: no F yet', 'F1', 'F2', 'F3: fully fluorinated', 'Si released']
    frames = []
    for k, t in enumerate(times):
        fig = plt.figure(figsize=(11, 5.4))
        fig.subplots_adjust(top=0.8, bottom=0.12, left=0.03, right=0.84, wspace=0.18)
        grid_ax = fig.add_subplot(1, 2, 1)
        s = grids[k].reshape(side, side)
        rgb = np.array([matplotlib.colors.to_rgb(c) for c in fs.RAMP4 + [fs.GRID]])
        img = rgb[np.where(removed[s], 4, stage[s])]
        grid_ax.imshow(img, interpolation='nearest')
        grid_ax.set_xticks(np.arange(-0.5, side, 1), minor=True)
        grid_ax.set_yticks(np.arange(-0.5, side, 1), minor=True)
        grid_ax.grid(which='minor', color=fs.SURFACE, lw=1.4)
        grid_ax.grid(which='major', visible=False)
        grid_ax.tick_params(which='both', length=0, labelbottom=False, labelleft=False)
        for sp in grid_ax.spines.values():
            sp.set_visible(False)
        rr, cc = np.nonzero(hf[s])
        grid_ax.scatter(cc, rr, s=28, color=fs.ORANGE, edgecolor=fs.SURFACE, lw=0.8, zorder=3)
        grid_ax.set_title(f'{n} surface sites (each square = one Si motif)', fontsize=10.5)
        area = fig.add_subplot(1, 2, 2)
        area.stackplot(times, 100 * fractions.T, colors=colors, lw=0)
        fs.shade_phases(area, [(0, dose, 'HF dose'), (dose, dose + purge, 'purge')])
        area.axvline(t, color=fs.INK, lw=1.2)
        area.set(xlim=(0, times[-1]), ylim=(0, 100), xlabel='Time (s)', ylabel='Share of sites (%)')
        area.set_title('How the surface changes over time', fontsize=10.5, pad=20)
        mids = np.cumsum(100 * fractions[-1]) - 50 * fractions[-1]
        for y, name, frac in zip(mids, names, fractions[-1]):
            if frac > 0.02:
                area.text(
                    times[-1] * 1.02, y, name, va='center', fontsize=9, color=fs.INK, clip_on=False
                )
        fig.legend(
            handles=[Patch(color=c, label=l) for c, l in zip(fs.RAMP4, names[:4])]
            + [Patch(color=fs.GRID, label='Si released (gone)')]
            + [Line2D([], [], ls='', marker='o', color=fs.ORANGE, label='HF adsorbed on site')],
            loc='lower center',
            ncol=6,
            fontsize=9,
            bbox_to_anchor=(0.45, 0.0),
        )
        phase = 'HF dose' if t < dose else 'purge (HF off)'
        released = 100 * fractions[k, 4]
        fs.header(
            fig,
            f'HF fluorinates surface sites until Si leaves as gas    t = {t:.2f} s, {phase}',
            '400 K. Blue'
            f" gets darker as F accumulates; a gray square means its Si "
            f'has left as gas ({released:.0f}% of sites so far). Recorded kMC states, no interpolation.',
        )
        frames.append(_grab(fig))
    return _gif(frames, 'species_kmc.gif')


def multilayer_animation():
    d = load_json('data/multilayer/demo_trajectory.json')
    nodes = d['initial']['nodes']
    params = d['parameters']
    pos = np.array([n['position_A'] for n in nodes])
    element = np.array([n['element'] for n in nodes])
    layer = np.array([n['layer'] for n in nodes])
    cut = (
        np.array([n['column'] for n in nodes]) // 6 == 2
    )  # same fixed row as the original renderer
    snaps = d['snapshots']
    cycle = params['dose_s'] + params['purge_s']
    horizon = snaps[-1]['time_s']
    spans = []
    t0 = 0.0
    while t0 < horizon - 1e-9:
        spans += [(t0, t0 + params['dose_s'], 'dose'), (t0 + params['dose_s'], t0 + cycle, 'purge')]
        t0 += cycle
    times = np.array([s['time_s'] for s in snaps])
    active = np.array([s['active'] for s in snaps])
    removed_si = [((~a) & (element == 'Si')).sum() for a in active]
    removed_n = [((~a) & (element == 'N')).sum() for a in active]
    per_band = Counter(layer)
    color = np.where(element == 'Si', fs.BLUE, fs.AQUA)
    zlo, zhi = pos[:, 2].min() - 1, pos[:, 2].max() + 1.5
    centers = [pos[layer == k, 2].mean() for k in range(6)]
    frames = []
    for k, snap in enumerate(snaps):
        act = active[k]
        exposed = np.array(snap['exposed'])
        hf = np.array(snap['precursors'])
        fig = plt.figure(figsize=(12, 6))
        gsp = fig.add_gridspec(
            2,
            2,
            width_ratios=[1.25, 1],
            hspace=0.65,
            wspace=0.28,
            left=0.07,
            right=0.97,
            top=0.8,
            bottom=0.16,
        )
        side = fig.add_subplot(gsp[:, 0])
        for i, j in snap['bonds']:
            if cut[i] and cut[j] and act[i] and act[j] and abs(pos[i, 0] - pos[j, 0]) < 4:
                side.plot(pos[[i, j], 0], pos[[i, j], 2], color=fs.AXIS, lw=1, zorder=1)
        m = cut & act
        side.scatter(pos[m, 0], pos[m, 2], s=60, c=color[m], edgecolor=fs.SURFACE, lw=0.8, zorder=3)
        m = cut & ~act
        side.scatter(
            pos[m, 0],
            pos[m, 2],
            s=60,
            facecolor='none',
            edgecolor=fs.AXIS,
            lw=1.2,
            ls=(0, (2, 2)),
            zorder=2,
        )
        m = cut & act & exposed
        side.scatter(
            pos[m, 0], pos[m, 2], s=150, facecolor='none', edgecolor=fs.ORANGE, lw=2, zorder=4
        )
        m = cut & act & hf
        side.scatter(pos[m, 0], pos[m, 2] + 0.9, s=36, marker='v', color=fs.INK, zorder=5)
        for band, z in enumerate(centers):
            side.text(
                pos[cut, 0].max() + 1.6,
                z,
                f'B{band}' + (' (fixed)' if band == 0 else ''),
                va='center',
                fontsize=8.5,
                color=fs.MUTED,
            )
        side.set(
            ylim=(zlo, zhi),
            xlim=(pos[cut, 0].min() - 1, pos[cut, 0].max() + 1),
            xlabel='x (Å)',
            ylabel='Height (Å)',
        )
        side.grid(False)
        side.set_title('Side view of one slice through the slab')

        line = fig.add_subplot(gsp[0, 1])
        line.plot(times, removed_n, color=fs.AQUA, label='N removed')
        line.plot(times, removed_si, color=fs.BLUE, label='Si removed')
        fs.shade_phases(line, spans)
        line.axvline(times[k], color=fs.INK, lw=1.2)
        line.set(
            xlim=(0, horizon),
            ylim=(0, max(removed_n) * 1.25 + 1),
            xlabel='Time (s)',
            ylabel='Atoms removed',
        )
        line.set_title('Removal over time (whole slab)', pad=18)
        line.legend(loc='lower right', fontsize=9)

        bars = fig.add_subplot(gsp[1, 1])
        gone = Counter(layer[~act])
        bands = list(range(5, -1, -1))
        y = np.arange(len(bands))[::-1]
        pct = [100 * gone.get(b, 0) / per_band[b] for b in bands]
        hb = bars.barh(y, pct, height=0.6, color=fs.BLUE)
        fs.bar_value_labels(bars, hb, fmt='{:.0f}%')
        bars.set_yticks(y, [f'B{b}' + (' top' if b == 5 else '') for b in bands])
        bars.set(xlim=(0, 60), xlabel='Share of band removed (%)')
        bars.grid(axis='x')
        bars.grid(axis='y', visible=False)
        bars.set_title('Where removal happened (depth band)')
        fig.legend(
            handles=[
                Line2D([], [], ls='', marker='o', color=fs.BLUE, ms=8, label='Si'),
                Line2D([], [], ls='', marker='o', color=fs.AQUA, ms=8, label='N'),
                Line2D(
                    [],
                    [],
                    ls='',
                    marker='o',
                    mfc='none',
                    mec=fs.ORANGE,
                    mew=2,
                    ms=11,
                    label='exposed to gas',
                ),
                Line2D([], [], ls='', marker='v', color=fs.INK, ms=6, label='HF adsorbed'),
                Line2D([], [], ls='', marker='o', mfc='none', mec=fs.AXIS, ms=8, label='removed'),
            ],
            loc='lower left',
            ncol=5,
            fontsize=9,
            bbox_to_anchor=(0.05, 0.0),
        )
        if k == len(snaps) - 1:
            phase = 'end of run'
        else:
            phase = 'HF dose' if times[k] % cycle < params['dose_s'] else 'purge'
        fs.header(
            fig,
            f'Etching removes the top layer, exposing the atoms below    t = {times[k]:.2f} s ({phase})',
            f"Multilayer bond-graph kMC, {len(nodes)} Si/N atoms in six depth bands, {params['temperature_K']:.0f} K. "
            'Demonstration rates, not validated; atoms stay at fixed crystal positions.',
        )
        frames.append(_grab(fig))
    return _gif(frames, 'multilayer_kmc.gif')


FIGURES = {
    'kinetics': kinetics,
    'priorities': priorities,
    'stall': multilayer_stall,
    'relaxation': relaxation,
    'cleavage': final_cleavage,
    'species_animation': species_animation,
    'multilayer_animation': multilayer_animation,
}


def main(argv):
    fs.apply()
    OUT.mkdir(parents=True, exist_ok=True)
    chosen = argv or list(FIGURES)
    for name in chosen:
        FIGURES[name]()
        print('drew', name)
    manifest_path = OUT / 'manifest.json'
    manifest = (
        json.loads(manifest_path.read_text(encoding='utf-8')) if manifest_path.exists() else {}
    )
    sha = lambda p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest()
    manifest.update(
        description='Presentation figures drawn from saved results only; no simulation is run.',
        runner='scripts/make_figures.py',
        runner_sha256=sha('scripts/make_figures.py'),
        style_sha256=sha('plasma_surface/figstyle.py'),
    )
    manifest.setdefault('inputs_sha256', {}).update({p: sha(p) for p in sorted(INPUTS)})
    manifest['figures_sha256'] = {
        f.name: hashlib.sha256(f.read_bytes()).hexdigest()
        for f in sorted(OUT.glob('*'))
        if f.suffix in ('.png', '.gif')
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')


if __name__ == '__main__':
    main(sys.argv[1:])
