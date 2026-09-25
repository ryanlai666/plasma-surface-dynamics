"""Conservative rate-evidence matching; a descriptor match alone is insufficient."""

import hashlib, json, math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
from collections import Counter
import numpy as np


def environment(graph, precursors, event, temperature):
    from .multilayer import neighbors

    nb = neighbors(graph)
    i = event['site']
    j = event['partner']
    roots = {i} | ({j} if j is not None else set())
    selected = set(roots)
    for _ in range(2):
        selected |= {k for q in list(selected) for k in nb[q]}
    tags = {
        q: json.dumps(
            dict(
                element=graph['nodes'][q]['element'],
                termination=dict(sorted(Counter(graph['nodes'][q]['termination']).items())),
                precursor_HF=bool(precursors[q]),
                role='reactive_Si' if q == i else 'bond_partner' if q == j else 'neighbor',
                coordination=len(nb[q]),
            ),
            sort_keys=True,
        )
        for q in selected
    }
    labels = tags.copy()
    for _ in range(3):
        labels = {
            q: hashlib.sha256(
                (
                    tags[q] + '|' + '|'.join(sorted(labels[k] for k in nb[q] if k in selected))
                ).encode()
            ).hexdigest()
            for q in selected
        }
    # Geometry context consists of substrate pair lengths, with periodic MIC in x/y.
    cell = np.asarray(graph['cell_A'])
    inv = np.linalg.inv(cell)
    dist = []
    for a, b in graph['bonds']:
        if a in selected and b in selected:
            delta = np.array(graph['nodes'][a]['position_A']) - graph['nodes'][b]['position_A']
            frac = delta @ inv
            frac[:2] -= np.round(frac[:2])
            length = np.linalg.norm(frac @ cell)
            dist.append([*sorted([labels[a], labels[b]]), round(float(length), 3)])
    detail = dict(
        reaction=event['kind'],
        temperature_K=float(temperature),
        root_labels=sorted(labels[q] for q in roots),
        shell_labels=sorted(labels.values()),
        bond_geometry=sorted(dist),
        description='Two-shell labelled graph plus substrate bond lengths; not a unique relaxed adsorbate geometry or electronic state',
    )
    key = hashlib.sha256(json.dumps(detail, sort_keys=True).encode()).hexdigest()
    return dict(
        key=key,
        descriptor=detail,
        atom_ids=sorted(selected),
        central_site=i,
        partner=j,
        central_termination=dict(Counter(graph['nodes'][i]['termination'])),
        central_coordination=len(nb[i]),
        partner_termination=dict(Counter(graph['nodes'][j]['termination']))
        if j is not None
        else None,
    )


def qualified_rate(record, env, event, context):
    if record is None:
        return None, 'no_matched_record'
    required = [
        'id',
        'environment_key',
        'reaction',
        'rate_s',
        'source',
        'method',
        'context',
        'geometry_match_verified',
        'applicability_reviewed',
        'reference_state',
        'uncertainty',
        'artifact_sha256',
    ]
    if any(k not in record for k in required):
        return None, 'incomplete_evidence'
    if record['environment_key'] != env['key'] or record['reaction'] != event['kind']:
        return None, 'environment_mismatch'
    if record['context'] != context:
        return None, 'condition_mismatch'
    if not record['geometry_match_verified'] or not record['applicability_reviewed']:
        return None, 'unreviewed_transfer'
    if not record['artifact_sha256']:
        return None, 'missing_artifacts'
    for name, sha in record['artifact_sha256'].items():
        path = (ROOT / name).resolve()
        if (
            not path.is_relative_to(ROOT)
            or not path.is_file()
            or hashlib.sha256(path.read_bytes()).hexdigest() != sha
        ):
            return None, 'artifact_mismatch'
    rate = record['rate_s']
    if not isinstance(rate, (int, float)) or not math.isfinite(rate) or rate < 0:
        return None, 'invalid_rate'
    uncertainty = record['uncertainty']
    if not isinstance(uncertainty, dict) or not uncertainty.get('basis'):
        return None, 'missing_uncertainty_basis'
    bounds = uncertainty.get('rate_s_interval', [])
    if (
        len(bounds) != 2
        or not all(isinstance(x, (int, float)) and math.isfinite(x) for x in bounds)
        or not 0 <= bounds[0] <= rate <= bounds[1]
    ):
        return None, 'invalid_uncertainty_interval'
    if record['method'] == 'DFT_TST':
        checks = record.get('checks', {})
        if event['kind'] == 'HF_adsorption':
            return None, 'adsorption_requires_flux_and_sticking'
        if record['reference_state'] != 'same_adsorbed_or_surface_state':
            return None, 'inconsistent_reference'
        if checks.get('imaginary_modes') != 1 or not all(
            checks.get(k, False)
            for k in [
                'IS_minimum',
                'FS_minimum',
                'TS_force_converged',
                'connectivity',
                'prefactor_or_free_energy',
                'basis_cell_convergence',
                'charge_spin_reviewed',
            ]
        ):
            return None, 'transition_state_evidence_incomplete'
    elif record['method'] == 'matched_experiment':
        if not record.get('held_out_validation') or not record.get('matched_film_and_recipe'):
            return None, 'experiment_mismatch'
    else:
        return None, 'unsupported_evidence_method'
    return float(rate), 'qualified_record'
