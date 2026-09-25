"""Curate small attributed public DFT subsets without extracting archive paths."""

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import tarfile
import zipfile
import numpy as np
from ase.io import read, write


def xyz_blocks(text):
    lines = text.splitlines(keepends=True)
    i = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        n = int(lines[i])
        end = i + n + 2
        if n < 1 or end > len(lines):
            raise ValueError('Truncated or invalid XYZ frame')
        yield ''.join(lines[i:end])
        i = end


def digest(payload):
    return hashlib.sha256(payload).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--etch-archive', type=Path, required=True)
    p.add_argument('--hcl-archive', type=Path, default=Path('data/raw/si_hcl/data.zip'))
    p.add_argument('--output', type=Path, default=Path('data/reference/atomistic'))
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    payload = args.etch_archive.read_bytes()
    if hashlib.md5(payload).hexdigest() != 'c0b18d62a9ba8905b29fe3e90689c325':
        raise ValueError('Unexpected etching archive MD5')
    selected = []
    index = []
    inventory = []
    qsd_profiles = []
    with tarfile.open(fileobj=io.BytesIO(payload)) as tar:
        for member in sorted(tar.getmembers(), key=lambda m: m.name):
            if not member.isfile() or not member.name.endswith('.extxyz'):
                continue
            if member.size > 30_000_000:
                raise ValueError('Unexpected member size')
            raw = tar.extractfile(member).read()
            blocks = list(xyz_blocks(raw.decode()))
            is_qsd = '/qsd/' in member.name
            if is_qsd:
                # Keep all geometries from the first source configuration's six pair scans.
                import re

                groups = [re.search(r'group=(\S+)', b.splitlines()[1]).group(1) for b in blocks]
                chosen = sorted(set(groups))[:6]
                indices = [i for i, g in enumerate(groups) if g in chosen]
            else:
                indices = sorted(set(np.linspace(0, len(blocks) - 1, 5, dtype=int).tolist()))
            inventory.append(
                dict(
                    member=member.name,
                    frames=len(blocks),
                    selected=len(indices),
                    sha256=digest(raw),
                    interpretation='DFT quasi-static drag, not validated TS'
                    if is_qsd
                    else 'DFT single-point labels',
                )
            )
            for i in indices:
                atoms = read(io.StringIO(blocks[i]), format='extxyz')
                energy = float(atoms.get_potential_energy())
                forces = atoms.get_forces()
                if not np.isfinite(energy) or not np.isfinite(forces).all():
                    raise ValueError('Nonfinite labels')
                if forces.shape != (len(atoms), 3):
                    raise ValueError('Force shape')
                row = dict(
                    subset_frame=len(selected),
                    source_member=member.name,
                    source_frame=i,
                    source_frame_sha256=digest(blocks[i].encode()),
                    natoms=len(atoms),
                    formula=atoms.get_chemical_formula(),
                    energy_eV=energy,
                    has_forces=True,
                    label_kind='published_DFT',
                    source_doi='10.5281/zenodo.19491140',
                    split_group=atoms.info.get('group', member.name),
                    role='quasi_static_drag_not_TS' if is_qsd else 'snapshot',
                )
                if is_qsd:
                    row.update(
                        distance_A=float(atoms.info['dist']),
                        source_e_diff=float(atoms.info['e_diff']),
                    )
                    qsd_profiles.append(row.copy())
                index.append(row)
                selected.append(blocks[i])
    (args.output / 'dft_subset.extxyz').write_text(
        ''.join(selected), encoding='utf-8', newline='\n'
    )
    # The source frame blocks are unmodified. Do not invent missing energies for HCl geometries.
    hcl = args.hcl_archive.read_bytes()
    if hashlib.md5(hcl).hexdigest() != '34644b3c55b44441dc64492759cb6d31':
        raise ValueError('Unexpected HCl archive')
    geometry_rows = []
    geometries = []
    with zipfile.ZipFile(io.BytesIO(hcl)) as z:
        for name in sorted(z.namelist()):
            if not name.endswith('/CONTCAR') or '__MACOSX' in name:
                continue
            raw = z.read(name)
            atoms = read(io.StringIO(raw.decode()), format='vasp')
            role = 'source_named_TS_unverified_here' if '/TS/' in name else 'adsorbate_geometry'
            atoms.info.update(source_member=name, source_doi='10.5281/zenodo.10211009', role=role)
            geometries.append(atoms)
            geometry_rows.append(
                dict(
                    frame=len(geometries) - 1,
                    source_member=name,
                    sha256=digest(raw),
                    natoms=len(atoms),
                    formula=atoms.get_chemical_formula(),
                    role=role,
                    has_energy=False,
                    has_forces=False,
                )
            )
    write(args.output / 'hcl_geometries.extxyz', geometries, format='extxyz')

    def save(name, obj):
        (args.output / name).write_text(
            json.dumps(obj, indent=2, allow_nan=False) + '\n', encoding='utf-8'
        )

    save('frame_index.json', index)
    save('archive_inventory.json', inventory)
    save('hcl_geometry_index.json', geometry_rows)
    save('qsd_profiles.json', qsd_profiles)
    provenance = dict(
        created_utc=datetime.now(timezone.utc).isoformat(),
        script_sha256=digest(Path(__file__).read_bytes()),
        sources=[
            dict(
                doi='10.5281/zenodo.19491140',
                url='https://zenodo.org/records/19491140',
                license='CC BY 4.0',
                authors=[
                    'Sangmin Oh',
                    'Jinmu You',
                    'Jaesun Kim',
                    'Jiho Lee',
                    'Seungwu Han',
                    'Youngho Kang',
                ],
                archive_sha256=digest(payload),
            ),
            dict(
                doi='10.5281/zenodo.10211009',
                url='https://zenodo.org/records/10211009',
                license='CC BY 4.0',
                authors=['Biel Martinez i Diaz', 'Jing Li', 'Hector Prats', 'Benoit Sklenard'],
                archive_sha256=digest(hcl),
            ),
        ],
        selection='Five evenly spaced source frames per non-QSD member; all frames from six lexicographically first QSD groups; all HCl CONTCARs',
        published_DFT_frames=len(index),
        geometry_only_frames=len(geometries),
        source_named_TS=sum(r['role'] == 'source_named_TS_unverified_here' for r in geometry_rows),
        limitations=[
            'No new physical validation from subsampling',
            'QSD peaks are not transition-state barriers',
            'HCl archive lacks energies/forces in CONTCAR files',
            'No nitrogen in the Si/O/C/F DFT subset',
            'XC and pseudopotential details require paper audit before combining with other DFT sources',
            'Etch snapshots were generated with an ML potential then labeled with DFT; not AIMD',
            'Split by complete trajectory/configuration family to avoid train-test leakage',
        ],
        output_sha256={
            name: digest((args.output / name).read_bytes())
            for name in ['dft_subset.extxyz', 'hcl_geometries.extxyz']
        },
    )
    save('provenance.json', provenance)
    print(
        json.dumps(
            {
                k: provenance[k]
                for k in ['published_DFT_frames', 'geometry_only_frames', 'source_named_TS']
            },
            indent=2,
        )
    )


if __name__ == '__main__':
    main()
