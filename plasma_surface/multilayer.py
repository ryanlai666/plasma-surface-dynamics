"""Bond-resolved multilayer HF etching prototype with explicit atom ledgers.

All transferred/local-environment rates are unvalidated demonstration parameters.
The strict policy blocks events without graph-environment-specific validated rates.
"""

from copy import deepcopy
from collections import Counter
import math
import numpy as np

KB = 8.617333262145e-5
GAS = {'HF': {'H': 1, 'F': 1}, 'H2': {'H': 2}, 'HCl': {'H': 1, 'Cl': 1}, 'NH3': {'N': 1, 'H': 3}}


def formula(counts):
    return ''.join(
        e + (str(counts[e]) if counts[e] > 1 else '')
        for e in ['Si', 'N', 'H', 'F', 'Cl']
        if counts.get(e, 0)
    )


def inventory(graph, precursors):
    c = Counter()
    for i, node in enumerate(graph['nodes']):
        if not node['active']:
            continue
        c[node['element']] += 1
        c.update(node['termination'])
        if precursors[i]:
            c.update({'H': 1, 'F': 1})
    return c


def neighbors(graph):
    result = [set() for _ in graph['nodes']]
    for a, b in graph['bonds']:
        if not graph['nodes'][a]['active'] or not graph['nodes'][b]['active']:
            raise ValueError('Bond to removed atom')
        result[a].add(b)
        result[b].add(a)
    return result


def validate(graph):
    nb = neighbors(graph)
    for i, node in enumerate(graph['nodes']):
        if not node['active']:
            continue
        target = 4 if node['element'] == 'Si' else 3
        if len(nb[i]) + len(node['termination']) != target:
            raise ValueError('Valence inventory violation at ' + str(i))
        if any(x not in ['H', 'F', 'Cl'] for x in node['termination']):
            raise ValueError('Unsupported termination')
        if node['element'] == 'N' and any(x != 'H' for x in node['termination']):
            raise ValueError('Only NHx anchors implemented')
    return nb


def accessibility(graph, penetration_A=2.0):
    tops = {}
    for node in graph['nodes']:
        if node['active']:
            tops[node['column']] = max(tops.get(node['column'], -math.inf), node['position_A'][2])
    depths = np.array(
        [
            max(0.0, tops.get(n['column'], n['position_A'][2]) - n['position_A'][2])
            if n['active']
            else math.inf
            for n in graph['nodes']
        ]
    )
    exposed = np.array(
        [n['active'] and not n['fixed'] and d < 1e-7 for n, d in zip(graph['nodes'], depths)]
    )
    reachable = np.array(
        [
            n['active'] and not n['fixed'] and d <= penetration_A + 1e-9
            for n, d in zip(graph['nodes'], depths)
        ]
    )
    return depths, exposed, reachable


def candidates(
    graph,
    precursors,
    temperature=450.0,
    arrival=5.0,
    desorption=100.0,
    penetration_A=2.0,
    attenuation_A=1.0,
    policy='validated_only',
    dose_on=True,
    rate_library=None,
):
    if policy not in ['demonstration', 'validated_only']:
        raise ValueError('Unknown rate policy')
    if min(temperature, attenuation_A) <= 0 or min(arrival, desorption, penetration_A) < 0:
        raise ValueError('Invalid physical parameter')
    nb = neighbors(graph)
    depth, exposed, reachable = accessibility(graph, penetration_A)
    events = []

    def add(kind, i, j, rate, barrier=None, reference=None):
        if rate >= 0:
            events.append(
                dict(
                    kind=kind,
                    site=i,
                    partner=j,
                    rate_s=float(rate),
                    barrier_eV=barrier,
                    reference_pathway=reference,
                    rate_status='unvalidated_local_environment_transfer'
                    if reference
                    else 'assumed_demonstration_rate',
                    depth_A=float(depth[i]),
                    exposed=bool(exposed[i]),
                )
            )

    def chem(kind, i, j, ea, ref=None):
        add(kind, i, j, 1e12 * math.exp(-ea / (KB * temperature)), ea, ref)

    for i, node in enumerate(graph['nodes']):
        if not node['active'] or node['fixed']:
            continue
        if precursors[i]:
            add('HF_desorption', i, None, desorption)
        if not reachable[i]:
            continue
        lig = Counter(node['termination'])
        if node['element'] == 'N':
            if not nb[i] and lig == {'H': 3}:
                chem('NH3_release', i, None, 0.45)
            continue
        if not nb[i] and len(node['termination']) == 4 and not precursors[i]:
            chem('Si_molecule_release', i, None, 0.53)
        if not precursors[i]:
            if dose_on:
                add('HF_adsorption', i, None, arrival * math.exp(-depth[i] / attenuation_A))
            continue
        fluor = lig['F']
        stage = min(fluor + 1, 3)
        for j in nb[i]:
            other = graph['nodes'][j]
            if other['fixed']:
                continue
            if other['element'] == 'N':
                h = other['termination'].count('H')
                if fluor >= 3:
                    if h == 1:
                        chem('SiN_cleavage', i, j, 0.53, 'P4')
                    else:
                        add('SiN_cleavage', i, j, 0.0)
                        events[-1]['rate_status'] = 'disabled_missing_final_cleavage_barrier'
                    continue
                suffix = 'a' if h >= 2 else 'c' if h == 1 else 'd'
                ref = f'P{stage}{suffix}'
                # Literature numbers only seed a hypothesis; graph environments differ.
                ea = {'a': [0.45, 0.68, 0.88], 'c': [0.02, 0.79, 0.74], 'd': [0.72, 0.75, 0.90]}[
                    suffix
                ][stage - 1]
                chem('SiN_cleavage', i, j, ea, ref)
            elif other['element'] == 'Si':
                chem('SiSi_cleavage', i, j, 0.79, 'P2f')
        if lig['H']:
            chem('SiH_exchange', i, None, [2.10, 1.80, 1.54][stage - 1], 'P' + str(stage) + 'b')
        if lig['Cl']:
            chem('SiCl_exchange', i, None, 0.75)
    if policy == 'validated_only':
        from .rate_evidence import environment, qualified_rate

        context = dict(
            temperature_K=float(temperature),
            arrival_s=float(arrival),
            desorption_s=float(desorption),
            penetration_A=float(penetration_A),
            attenuation_A=float(attenuation_A),
        )
        accepted = []
        for event in events:
            env = environment(graph, precursors, event, temperature)
            record = (rate_library or {}).get(env['key'])
            rate, reason = qualified_rate(record, env, event, context)
            if rate is not None:
                event.update(
                    rate_s=rate,
                    rate_status='qualified_record',
                    rate_record=record['id'],
                    environment_key=env['key'],
                )
                accepted.append(event)
        return accepted
    return events


def apply(graph, precursors, event, gas, gas_formulas):
    i = event['site']
    node = graph['nodes'][i]
    kind = event['kind']
    j = event['partner']

    def emit(name, counts):
        gas[name] = gas.get(name, 0) + 1
        gas_formulas[name] = counts

    if kind == 'HF_adsorption':
        if precursors[i] or not node['active']:
            raise ValueError('Occupied/inactive site')
        precursors[i] = True
        gas['HF'] = gas.get('HF', 0) - 1
    elif kind == 'HF_desorption':
        if not precursors[i]:
            raise ValueError('Missing precursor')
        precursors[i] = False
        emit('HF', GAS['HF'])
    elif kind in ['SiN_cleavage', 'SiSi_cleavage']:
        if not precursors[i]:
            raise ValueError('Missing precursor')
        pair = sorted([i, j])
        if pair not in graph['bonds']:
            raise ValueError('Missing bond')
        graph['bonds'].remove(pair)
        node['termination'].append('F')
        graph['nodes'][j]['termination'].append('H')
        precursors[i] = False
    elif kind in ['SiH_exchange', 'SiCl_exchange']:
        if not precursors[i]:
            raise ValueError('Missing precursor')
        old = 'H' if kind == 'SiH_exchange' else 'Cl'
        node['termination'].remove(old)
        node['termination'].append('F')
        precursors[i] = False
        name = 'H2' if old == 'H' else 'HCl'
        emit(name, GAS[name])
    elif kind in ['Si_molecule_release', 'NH3_release']:
        if neighbors(graph)[i] or precursors[i]:
            raise ValueError('Cannot release bonded/occupied atom')
        counts = dict(Counter([node['element']] + node['termination']))
        name = formula(counts)
        if kind == 'NH3_release' and counts != {'N': 1, 'H': 3}:
            raise ValueError('Invalid NH3 product')
        if kind == 'Si_molecule_release' and (counts.get('Si') != 1 or sum(counts.values()) != 5):
            raise ValueError('Incomplete Si molecule')
        emit(name, counts)
        node['active'] = False
        node['termination'] = []
    else:
        raise ValueError('Unknown event')


def simulate(
    initial,
    times,
    seed=71,
    temperature=450.0,
    arrival=5.0,
    desorption=100.0,
    penetration_A=2.0,
    attenuation_A=1.0,
    policy='validated_only',
    dose_s=2.0,
    purge_s=1.0,
    max_events=200000,
    rate_library=None,
):
    graph = deepcopy(initial)
    validate(graph)
    prec = np.zeros(len(graph['nodes']), bool)
    initial_ledger = inventory(graph, prec)
    gas = {}
    gas_formulas = deepcopy(GAS)
    rng = np.random.default_rng(seed)
    times = np.asarray(times, float)
    if not len(times) or times[0] < 0 or np.any(np.diff(times) < 0):
        raise ValueError('Invalid sample times')
    if dose_s < 0 or purge_s < 0 or dose_s + purge_s <= 0:
        raise ValueError('Invalid phases')
    history = []
    events = []
    t = 0.0
    sample = 0
    period = dose_s + purge_s

    def record(time):
        d, ex, re = accessibility(graph, penetration_A)
        history.append(
            dict(
                time_s=float(time),
                active=[n['active'] for n in graph['nodes']],
                terminations=[dict(Counter(n['termination'])) for n in graph['nodes']],
                precursors=prec.tolist(),
                gas=deepcopy(gas),
                bonds=deepcopy(graph['bonds']),
                depth_A=[None if not math.isfinite(x) else float(x) for x in d],
                exposed=ex.tolist(),
                reachable=re.tolist(),
            )
        )

    while sample < len(times):
        cycle = math.floor((t + 1e-12) / period)
        phase = t - cycle * period
        dose = phase < dose_s - 1e-10
        boundary = cycle * period + dose_s if dose else (cycle + 1) * period
        if boundary <= t + 1e-12:
            boundary = t + period
        options = candidates(
            graph,
            prec,
            temperature,
            arrival,
            desorption,
            penetration_A,
            attenuation_A,
            policy,
            dose,
            rate_library,
        )
        rates = np.array([e['rate_s'] for e in options])
        total = rates.sum()
        proposed = t + rng.exponential(1 / total) if total else math.inf
        next_t = min(proposed, boundary)
        while sample < len(times) and times[sample] <= next_t:
            record(times[sample])
            sample += 1
        if sample == len(times):
            break
        if boundary <= proposed:
            t = boundary
            continue
        idx = min(
            int(np.searchsorted(np.cumsum(rates), rng.random() * total, side='right')),
            len(options) - 1,
        )
        event = options[idx]
        before_ex = accessibility(graph, penetration_A)[1]
        apply(graph, prec, event, gas, gas_formulas)
        validate(graph)
        ledger = inventory(graph, prec)
        for name, amount in gas.items():
            for element, count in gas_formulas[name].items():
                ledger[element] += amount * count
        assert all(ledger[e] == initial_ledger[e] for e in set(ledger) | set(initial_ledger))
        after_ex = accessibility(graph, penetration_A)[1]
        event.update(
            time_s=float(proposed), newly_exposed=np.flatnonzero(after_ex & ~before_ex).tolist()
        )
        events.append(event)
        t = proposed
        if len(events) > max_events:
            raise RuntimeError('Event budget exceeded')
    return dict(
        policy=policy,
        parameters=dict(
            temperature_K=temperature,
            arrival_s=arrival,
            desorption_s=desorption,
            penetration_A=penetration_A,
            attenuation_A=attenuation_A,
            dose_s=dose_s,
            purge_s=purge_s,
        ),
        initial=initial,
        final=graph,
        snapshots=history,
        events=events,
        gas_formulas=gas_formulas,
        initial_inventory=dict(initial_ledger),
        conservation_exact=True,
        seed=seed,
        scope='Pulsed graph-reactive etching demonstration; rates unvalidated; no ALE self-limitation or EPC claim',
    )
