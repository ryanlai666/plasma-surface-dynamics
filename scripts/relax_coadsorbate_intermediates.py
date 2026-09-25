"""Local coadsorbate screening to address A1/A5, not a barrier or rate fit."""

from pathlib import Path
import csv, hashlib, json, time
import numpy as np
import torch
from threadpoolctl import threadpool_limits
from ase import Atoms
from ase.io import read, write
from ase.optimize import FIRE
from ase.constraints import FixAtoms
from mace.calculators import MACECalculator
from molecular_surface_campaign import ROOT, independent, snapshot, molecule

OUT = ROOT / 'data/intermediate_campaign'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def optimize(a, calc, folder, steps=400, tolerance=0.025):
    folder.mkdir(parents=True, exist_ok=True)
    a.calc = independent(calc)
    frames = []
    opt = FIRE(a, logfile=str(folder / 'optimization.log'), dt=0.03, maxstep=0.07)

    def save():
        frames.append(snapshot(a))

    opt.attach(save, interval=1)
    start = time.perf_counter()
    ok = bool(opt.run(fmax=tolerance, steps=steps))
    write(folder / 'trajectory.extxyz', frames)
    write(folder / 'final.extxyz', frames[-1])
    result = dict(
        converged=ok,
        steps=opt.nsteps,
        energy_eV=float(a.get_potential_energy()),
        max_mobile_force_eV_A=float(np.linalg.norm(a.get_forces(), axis=1).max()),
        elapsed_s=time.perf_counter() - start,
        frames=len(frames),
    )
    with (folder / 'energies.csv').open('w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['optimization_step', 'energy_eV', 'energy_relative_to_start_eV'])
        energies = [float(x.get_potential_energy()) for x in frames]
        w.writerows((i, e, e - energies[0]) for i, e in enumerate(energies))
    return a, result


def adsorbate_curvature(a, first, delta=0.01):
    """Unweighted Hessian block; host frozen during finite differences, not full stability."""
    initial = a.positions.copy()
    dofs = [(i, j) for i in range(first, len(a)) for j in range(3)]
    h = np.zeros((len(dofs), len(dofs)))
    for col, (i, j) in enumerate(dofs):
        p = initial.copy()
        p[i, j] += delta
        a.set_positions(p)
        fp = a.get_forces(apply_constraint=False)[first:].ravel()
        p = initial.copy()
        p[i, j] -= delta
        a.set_positions(p)
        fm = a.get_forces(apply_constraint=False)[first:].ravel()
        h[:, col] = -(fp - fm) / (2 * delta)
    a.set_positions(initial)
    values = np.linalg.eigvalsh((h + h.T) / 2)
    return dict(
        displacement_A=delta,
        block='adsorbate Cartesian coordinates; every host coordinate held fixed',
        eigenvalues_eV_A2=values.tolist(),
        negative_eigenvalues_below_minus_0p02=int(sum(values < -0.02)),
        minimum_eigenvalue_eV_A2=float(values.min()),
        scope='Partial curvature screen only; no full mobile-host Hessian or vibrational thermochemistry',
    )


def placed_addition(base, species, start):
    # HF is last F,H in the parent. Approach along an upward/lateral donor axis,
    # with the second molecule accepting the parent's H via F or O.
    direction = np.array([0.65 if start == 1 else -0.65, 0.35, 0.68])
    direction /= np.linalg.norm(direction)
    center = base.positions[-1] + 1.8 * direction
    if species == 'HF':
        added = Atoms('FH', positions=[center, center + 0.94 * direction])
    else:
        tangent = np.cross(direction, [0.0, 1.0, 0.0])
        tangent /= np.linalg.norm(tangent)
        added = Atoms(
            'OH2',
            positions=[
                center,
                center + 0.586 * direction + 0.757 * tangent,
                center + 0.586 * direction - 0.757 * tangent,
            ],
        )
    a = base.copy()
    a += added
    # Reject overlapping starts under the actual periodic geometry.
    dist = a.get_all_distances(mic=True)
    cross = dist[len(base) :, : len(base)]
    if cross.min() < 1.0:
        raise ValueError('Overlapping coadsorbate start')
    return a, float(cross.min())


def main():
    OUT.mkdir(exist_ok=True)
    torch.set_num_threads(2)
    threadpool_limits(limits=1)
    model = ROOT / '.cache/mace/mp_0b2_small.model'
    calc = MACECalculator(model_paths=str(model), device='cpu', default_dtype='float64')
    refs = {}
    results = []
    for species in ('HF', 'H2O'):
        gas = molecule(species)
        gas, info = optimize(gas, calc, OUT / 'gas' / species, steps=200, tolerance=0.01)
        refs[species] = info
        (OUT / 'gas' / species / 'summary.json').write_text(json.dumps(info, indent=2) + '\n')
    for surface in ('beta_Si3N4_001', 'alpha_quartz_001'):
        for start in (1, 2):
            source = (
                ROOT
                / 'data/surface_paths'
                / surface
                / 'HF/relaxed_adsorption'
                / f'start_{start}/trajectory.extxyz'
            )
            parent = read(source, -1)
            nhost = len(parent) - 2
            base, mono = optimize(parent, calc, OUT / surface / f'parent_{start}')
            for species in ('HF', 'H2O'):
                folder = OUT / surface / ('HF_' + species) / f'start_{start}'
                if (folder / 'summary.json').exists():
                    results.append(json.loads((folder / 'summary.json').read_text()))
                    continue
                a, separation = placed_addition(base, species, start)
                fixed = base.constraints[0].get_indices().tolist()
                a.set_constraint(FixAtoms(indices=fixed))
                final, r = optimize(a, calc, folder)
                distances = final.get_all_distances(mic=True)
                structure = dict(
                    parent_HF_distance_A=float(distances[nhost, nhost + 1]),
                    added_bond_distances_A=[
                        float(distances[nhost + 2, i]) for i in range(nhost + 3, len(final))
                    ],
                    closest_added_heavy_to_parent_H_A=float(distances[nhost + 2, nhost + 1]),
                    minimum_added_to_host_A=float(distances[nhost + 2 :, :nhost].min()),
                )
                eadd = (
                    float(r['energy_eV'] - mono['energy_eV'] - refs[species]['energy_eV'])
                    if r['converged'] and mono['converged'] and refs[species]['converged']
                    else None
                )
                curvature = adsorbate_curvature(final, nhost) if r['converged'] else None
                r.update(
                    surface=surface,
                    coadsorbate=species,
                    start=start,
                    gap_ids=['A1', 'A5'],
                    source=str(source.relative_to(ROOT)).replace('\\', '/'),
                    source_sha256=sha(source),
                    fixed_indices=fixed,
                    host_atoms=nhost,
                    adsorbate_atoms=len(final) - nhost,
                    formula=final.get_chemical_formula(),
                    initial_min_cross_distance_A=separation,
                    incremental_association_energy_eV=eadd,
                    reference=mono,
                    gas_reference=refs[species],
                    geometry=structure,
                    adsorbate_curvature=curvature,
                    folder=str(folder.relative_to(ROOT)).replace('\\', '/'),
                    trajectory_sha256=sha(folder / 'trajectory.extxyz'),
                    final_sha256=sha(folder / 'final.extxyz'),
                    status='relaxed_candidate_not_TS_or_rate_evidence'
                    if r['converged']
                    else 'unconverged_candidate',
                    rate_enabled=False,
                )
                (folder / 'summary.json').write_text(json.dumps(r, indent=2) + '\n')
                results.append(r)
                print(
                    surface,
                    species,
                    start,
                    r['converged'],
                    r['max_mobile_force_eV_A'],
                    eadd,
                    flush=True,
                )
    report = dict(
        method='MACE-MP-0b2 small, float64 CPU; FIRE upper host and all adsorbates relaxed; lower host fixed',
        runner='scripts/relax_coadsorbate_intermediates.py',
        runner_sha256=sha(Path(__file__)),
        helper='scripts/molecular_surface_campaign.py',
        helper_sha256=sha(ROOT / 'scripts/molecular_surface_campaign.py'),
        model_sha256=sha(model),
        force_tolerance_eV_A=0.025,
        maximum_steps=400,
        results=results,
        limitations=[
            'Ideal crystal cuts, not amorphous SiNx:H or fluorinated clusters from source papers',
            'Ground-state potential: no vibrationally excited HF',
            'Only two local starts per surface and coadsorbate; periodic coverage fixed by cell',
            'Partial adsorbate curvature is not full minimum validation',
            'No TS, free energy, qualified rate, or experimental validation',
        ],
    )
    (OUT / 'summary.json').write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    main()
