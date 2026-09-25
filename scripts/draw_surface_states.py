"""Chemical bookkeeping sketches for every kMC state, with exact atom inventories.

Coordinates are diagram layouts, not calculated adsorption geometries.
"""

from pathlib import Path
from collections import Counter
import json, hashlib, math
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, FancyArrowPatch, Patch
from draw_species_network import COLORS, HF, REMOVED, MISSING, NAMES, save, ROOT

ELEMENT = {'Si': '#90A4BD', 'N': '#678BDD', 'F': '#71C98D', 'H': '#FFFFFF'}


def depiction(s):
    nodes = []
    bonds = []
    anchors = []

    def atom(el, x, y):
        nodes.append(dict(element=el, x=x, y=y))
        return len(nodes) - 1

    def bond(a, b, kind='bond'):
        bonds.append(dict(a=a, b=b, kind=kind))

    def group(h, x, y, central=None, anchor=False):
        i = atom('N', x, y)
        if central is not None:
            bond(central, i)
        if anchor:
            anchors.append(i)
        for k in range(h):
            angle = math.radians(60 + 120 * k) if h == 2 else math.radians(90)
            j = atom('H', x + 0.11 * math.cos(angle), y + 0.11 * math.sin(angle))
            bond(i, j)
        return i

    j = s['fluorination_stage']
    family = s['family']
    removed = s['Si_removed']
    unknown = removed and family in ['terminal_NH2', 'bridge_N']
    if unknown:
        # No validated residual connectivity: an inventory is not a chemical structure.
        return dict(
            state_id=s['id'],
            nodes=[],
            bonds=[],
            anchors=[],
            unknown=True,
            composition=s['composition'],
            scope='Unparameterized placeholder; structure deliberately not drawn',
        )
    slots = [(0.25, 0.57), (0.75, 0.57), (0.25, 0.31), (0.75, 0.31)]
    if family in ['bridge_NH', 'bridge_N', 'terminal_NH2']:
        si = None if removed else atom('Si', 0.50, 0.43)
        for k, (x, y) in enumerate(slots):
            if k < j:
                if si is not None:
                    bond(si, atom('F', x, y))
                if family in ['bridge_NH', 'bridge_N']:
                    group(2 if family == 'bridge_NH' else 1, 0.19 + 0.205 * k, 0.13, anchor=True)
            else:
                group(
                    1 if family == 'bridge_NH' else 0 if family == 'bridge_N' else 2,
                    x,
                    y,
                    si,
                    anchor=family != 'terminal_NH2',
                )
    elif family == 'Si_H_rich':
        if removed:
            group(2, 0.5, 0.23, anchor=True)
        else:
            si = atom('Si', 0.5, 0.43)
            group(1, 0.5, 0.20, si, True)
            for k, (x, y) in enumerate([(0.24, 0.43), (0.5, 0.66), (0.76, 0.43)]):
                bond(si, atom('F' if k < j else 'H', x, y))
    elif family.startswith('pre_'):
        if removed:
            group(1, 0.5, 0.23, anchor=True)
        else:
            si = atom('Si', 0.5, 0.43)
            group(0, 0.5, 0.20, si, True)
            for k, (x, y) in enumerate([(0.24, 0.43), (0.5, 0.66), (0.76, 0.43)]):
                bond(si, atom('F' if k < j else 'H', x, y))
    elif family == 'Si_Si_backbond':
        si = atom('Si', 0.44, 0.47)
        other = atom('Si', 0.70, 0.23)
        anchors.extend([si, other])
        bond(si, atom('F', 0.22, 0.56))
        if j == 1:
            bond(si, other)
        else:
            bond(si, atom('F', 0.63, 0.61))
            bond(other, atom('H', 0.83, 0.36))
    if s['HF_complex']:
        f = atom('F', 0.40, 0.84)
        h = atom('H', 0.60, 0.84)
        bond(f, h)
        bond(h, 0, 'encounter')
    inventory = dict(Counter(a['element'] for a in nodes))
    assert inventory == s['composition'], (s['id'], inventory, s['composition'])
    return dict(
        state_id=s['id'],
        nodes=nodes,
        bonds=bonds,
        anchors=anchors,
        unknown=False,
        composition=inventory,
        scope='Representative reduced-motif connectivity; not unique relaxed geometry; untracked substrate anchors implicit',
    )


def draw(ax, d, s, x=0, y=0, w=1, h=1, labels=True):
    def point(a):
        return x + w * a['x'], y + h * a['y']

    unknown = d['unknown']
    color = (
        'white'
        if unknown
        else HF
        if s['HF_complex']
        else REMOVED
        if s['Si_removed']
        else COLORS[min(s['fluorination_stage'], 3)]
    )
    ax.add_patch(
        Rectangle(
            (x, y - 0.10 * h),
            w,
            h * 1.14,
            facecolor='#FCFDFE',
            edgecolor=MISSING if unknown else '#D2DCE6',
            lw=0.7,
            linestyle='--' if unknown else '-',
        )
    )
    ax.add_patch(
        Rectangle(
            (x, y - 0.10 * h),
            w,
            0.12 * h,
            facecolor='#E4E9EE',
            hatch='///',
            edgecolor='#ABB7C3',
            lw=0.4,
        )
    )
    ax.add_patch(Rectangle((x, y + h * 0.95), w, 0.10 * h, facecolor=color, edgecolor='none'))
    title = (
        'unknown residual'
        if unknown
        else 'HF adsorbed'
        if s['HF_complex']
        else 'Si released'
        if s['Si_removed']
        else 'F' + str(s['fluorination_stage']) + ' motif'
    )
    ax.text(
        x + w / 2,
        y + h,
        title,
        ha='center',
        va='center',
        fontsize=7.3,
        color='white' if s['Si_removed'] and not unknown else '#193148',
    )
    if unknown:
        ax.text(x + w / 2, y + h * 0.50, '?', fontsize=27, ha='center', color=MISSING)
        ax.text(
            x + w / 2,
            y + h * 0.24,
            'No verified structure',
            fontsize=6.5,
            ha='center',
            color=MISSING,
        )
        return
    for i in d['anchors']:
        px, py = point(d['nodes'][i])
        ax.plot([px, px], [y, py], c='#A0ACB8', lw=0.8, ls=':', zorder=1)
    for b in d['bonds']:
        a, c = point(d['nodes'][b['a']]), point(d['nodes'][b['b']])
        ax.plot(
            [a[0], c[0]],
            [a[1], c[1]],
            color='#C68C23' if b['kind'] == 'encounter' else '#485C6E',
            lw=1.1,
            ls=':' if b['kind'] == 'encounter' else '-',
            zorder=2,
        )
    for a in d['nodes']:
        px, py = point(a)
        radius = w * (0.037 if a['element'] == 'H' else 0.050)
        ax.add_patch(
            Circle(
                (px, py),
                radius,
                facecolor=ELEMENT[a['element']],
                edgecolor='#42576B',
                lw=0.6,
                zorder=3,
            )
        )
        if labels:
            ax.text(
                px,
                py,
                a['element'],
                fontsize=5.1 if a['element'] == 'H' else 5.8,
                ha='center',
                va='center',
                zorder=4,
                color='#193148',
            )


def molecule(ax, name, x, y, scale=0.17):
    # Flat molecular connectivity icon; never a claimed optimized conformation.
    center = 'Si' if name.startswith('Si') else 'N' if name == 'NH3' else 'H'
    others = {
        'SiF4': ['F'] * 4,
        'SiH2F2': ['H', 'H', 'F', 'F'],
        'SiHF3': ['H', 'F', 'F', 'F'],
        'NH3': ['H'] * 3,
        'H2': ['H'],
        'HF': ['F'],
    }[name]
    ax.scatter([x], [y], s=48, c=ELEMENT[center], edgecolors='#43576A', zorder=4)
    ax.text(x, y, center, ha='center', va='center', fontsize=4.5, zorder=5)
    for k, el in enumerate(others):
        a = 2 * math.pi * k / len(others)
        xx = x + scale * math.cos(a)
        yy = y + scale * math.sin(a)
        ax.plot([x, xx], [y, yy], c='#485C6E', lw=0.7)
        ax.scatter([xx], [yy], s=30, c=ELEMENT[el], edgecolors='#43576A', linewidths=0.5, zorder=4)
        ax.text(xx, yy, el, ha='center', va='center', fontsize=4, zorder=5)


def family_panel(ax, n, family, dep):
    states = n['states']
    ids = [i for i, s in enumerate(states) if s['family'] == family]
    events = [e for e in n['events'] if e['source'] in ids]
    base = sorted(
        [i for i in ids if not states[i]['HF_complex']],
        key=lambda i: states[i]['fluorination_stage'],
    )
    order = []
    for i in base:
        order.append(i)
        order.extend(
            j
            for j in ids
            if states[j]['HF_complex']
            and states[j]['fluorination_stage'] == states[i]['fluorination_stage']
        )
    positions = {i: k * 1.60 for k, i in enumerate(order)}
    width = 1.12
    for i in order:
        draw(ax, dep[i], states[i], positions[i], 0.20, width, 1.15)
    for e in events:
        a = positions[e['source']]
        b = positions[e['target']]
        if e['kind'] == 'desorption':
            ax.add_patch(
                FancyArrowPatch(
                    (a + 0.3, 1.49),
                    (b + 0.8, 1.49),
                    arrowstyle='-|>',
                    connectionstyle='arc3,rad=.28',
                    mutation_scale=8,
                    color='#7A8A98',
                    lw=0.8,
                )
            )
        else:
            ax.add_patch(
                FancyArrowPatch(
                    (a + width, 0.80),
                    (b, 0.80),
                    arrowstyle='-|>',
                    mutation_scale=9,
                    color='#B45928' if e['kind'] == 'reaction' else '#7A8A98',
                    lw=1.1,
                )
            )
            x = (a + width + b) / 2
            if e['kind'] == 'reaction':
                ax.text(
                    x,
                    0.48,
                    e['pathway'] + '\n' + f"{e['barrier_eV']:.2f} eV",
                    ha='center',
                    va='top',
                    fontsize=6,
                    color='#843F1E',
                )
                products = [g for g, v in e['gas_delta'].items() if v > 0]
                for gas in products:
                    molecule(ax, gas, x, -0.12, scale=0.13)
                    ax.text(x, -0.43, gas + ' (g)', ha='center', fontsize=6, color='#843F1E')
            else:
                ax.text(x, 0.50, '+ HF', ha='center', fontsize=6, color='#637889')
    if family in ['terminal_NH2', 'bridge_N']:
        a = positions[order[-2]] + width
        b = positions[order[-1]]
        ax.add_patch(
            FancyArrowPatch(
                (a, 0.8), (b, 0.8), arrowstyle='-|>', mutation_scale=9, color=MISSING, ls='--', lw=1
            )
        )
        ax.text(
            (a + b) / 2, 0.46, 'missing\nbarrier', ha='center', va='top', fontsize=6, color=MISSING
        )
    if family == 'Si_Si_backbond':
        ax.text(5.2, 0.75, 'Further states and rates unresolved', fontsize=10, color=MISSING)
    ax.text(0, 1.94, NAMES[family], fontsize=12, weight='bold', color='#193148')
    ax.text(
        14,
        1.94,
        f'{len(ids)} states / {len(events)} events',
        ha='right',
        fontsize=9,
        color='#637889',
    )
    ax.set(xlim=(-0.05, 14.2), ylim=(-0.60, 2.1))
    ax.set_aspect('equal')
    ax.axis('off')


def main():
    source = ROOT / 'configs/species_kmc_network.json'
    n = json.loads(source.read_text())
    dep = [depiction(s) for s in n['states']]
    out = ROOT / 'docs/reaction_network'
    folder = out / 'surface_states'
    folder.mkdir(exist_ok=True)
    families = list(dict.fromkeys(s['family'] for s in n['states']))
    for s, d in zip(n['states'], dep):
        fig, ax = plt.subplots(figsize=(3.6, 3.8))
        draw(ax, d, s)
        ax.set(xlim=(-0.08, 1.08), ylim=(-0.16, 1.12))
        ax.set_aspect('equal')
        ax.axis('off')
        fig.suptitle(s['id'], fontsize=10)
        fig.text(0.5, 0.02, 'Schematic motif, not calculated geometry', ha='center', fontsize=8)
        save(fig, folder / s['id'])
    for family in families:
        fig, ax = plt.subplots(figsize=(20, 4.1))
        family_panel(ax, n, family, dep)
        fig.subplots_adjust(left=0.025, right=0.99, bottom=0.14, top=0.98)
        fig.text(
            0.5,
            0.035,
            'Hatched band: implicit substrate | dotted gold: HF encounter, not a covalent bond | gray return arrows: HF desorption',
            ha='center',
            fontsize=10,
        )
        save(fig, out / ('surface_network_' + family))
    fig, axes = plt.subplots(7, 1, figsize=(21, 25))
    fig.subplots_adjust(left=0.035, right=0.985, top=0.925, bottom=0.035, hspace=0.08)
    for ax, family in zip(axes, families):
        family_panel(ax, n, family, dep)
    fig.suptitle(
        'Surface motifs and adsorbed molecules in the kMC network',
        fontsize=22,
        weight='bold',
        y=0.989,
        color='#193148',
    )
    fig.text(
        0.5,
        0.970,
        '45 states / 55 enabled events | Representative bookkeeping connectivity, not relaxed atomic structures',
        ha='center',
        fontsize=11,
    )
    handles = [Patch(fc=ELEMENT[e], ec='#43576A', label=e) for e in ELEMENT] + [
        Patch(fc=HF, label='HF complex'),
        Patch(fc=REMOVED, label='Si released'),
        Patch(fc='white', ec=MISSING, ls='--', label='Structure / rate unresolved'),
    ]
    fig.legend(
        handles=handles,
        ncol=7,
        loc='upper center',
        bbox_to_anchor=(0.5, 0.962),
        frameon=False,
        fontsize=10,
    )
    fig.text(
        0.5,
        0.942,
        'State strips retain the kMC F0-F3 colors. Gray return arrows: HF desorption. Brown arrows: source Ea and gas products.',
        ha='center',
        fontsize=10,
    )
    fig.text(
        0.5,
        0.016,
        'Hatched substrate and dotted anchor lines are schematic context, not simulated lower layers. Dotted gold lines mark HF encounter complexes.\nOnly tracked atoms are drawn; untracked bonds/anchors are implicit. Missing residual structures are not invented.',
        ha='center',
        fontsize=10,
    )
    save(fig, out / 'species_kmc_network')
    manifest = json.loads((out / 'render_manifest.json').read_text())
    manifest.update(
        depictions=dep,
        depiction_renderer='scripts/draw_surface_states.py',
        depiction_renderer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        interpretation='Complete configured kMC graph with representative bookkeeping sketches; no calculated adsorption geometry or extra subsurface kinetics implied',
        network_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
    )
    (out / 'render_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    lines = [
        '# Surface and adsorbate reaction-network atlas',
        '',
        'These diagrams show representative reduced-motif connectivity, not atomistic structures. All explicitly depicted atoms match the state composition. Unparameterized residuals have no invented structure; their inventory remains in the configuration. Implicit anchors outside that inventory are drawn only as dotted connections to a hatched support. The support does not add a simulated subsurface layer.',
        '',
        'HF is drawn as H-F above the motif with a dotted encounter line. This line does not specify a binding atom, bond order, orientation or dissociation. Products beside chemical arrows are separate gas molecules.',
        '',
        '| Family | Full diagram |',
        '|---|---|',
    ]
    for family in families:
        lines.append(
            f'| {NAMES[family]} | [PNG](surface_network_{family}.png), [SVG](surface_network_{family}.svg) |'
        )
    lines += [
        '',
        '## Every state',
        '',
        '| State ID | Sketch | Explicit inventory |',
        '|---|---|---|',
    ]
    for s, d in zip(n['states'], dep):
        lines.append(
            f"| {s['id']} | [PNG](surface_states/{s['id']}.png), [SVG](surface_states/{s['id']}.svg) | {s['composition']}"
            + (' (unknown connectivity)' if d['unknown'] else '')
            + ' |'
        )
    lines += [
        '',
        'The same formula does not determine a unique surface structure. These sketches visualize the current model assumptions; obtaining validated connectivity requires relaxed structures, charge/spin information and reaction-path checks.',
        '',
        'Reproduce with `python scripts/plot_reaction_networks.py`. [Configuration](../../configs/species_kmc_network.json), [rendered atom inventories](render_manifest.json).',
    ]
    (out / 'SURFACE_ATLAS.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    print('Drew 45 state cards, 7 family panels and complete surface network')


if __name__ == '__main__':
    main()
