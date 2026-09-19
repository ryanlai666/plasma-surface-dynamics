"""Summarize saved calculations without reclassifying failed paths as transition states."""
from pathlib import Path
import csv,json
ROOT=Path(__file__).resolve().parents[1]
def main():
    results=json.loads((ROOT/'docs/dry_etch_results/molecular_screening.json').read_text())
    assert len(results)==36
    lines=['## Calculated screening results','','![Molecular screening](molecular_screening.png)','','| Surface | Molecule | Site / orientation at lowest sampled energy | E_int (eV) | Height (angstrom) | Data |','|---|---|---|---:|---:|---|']
    for r in results:
        parent='data/surface_paths/'+r['surface']+'/'+r['species']
        lines.append(f"| {r['surface']} | {r['species']} | {r['site']} / {r['orientation']} | {r['minimum_interaction_energy_eV']:.4f} | {r['height_A']:.2f} | [curves + GIF](../../{parent}/README.md) |")
    reaction=[];rows=[]
    for method in ['molecular_neb','molecular_refined','mace_local_saddle','omol25_neb','omol25_refined']:
        for p in sorted((ROOT/'data/surface_paths').glob('*/*/'+method+'/summary.json')):
            d=json.loads(p.read_text());folder=p.parent;rel=folder.relative_to(ROOT).as_posix()
            peak=d.get('peak_above_IS_eV');valid=d.get('endpoint_connectivity_confirmed',False) and d.get('index_one_in_mobile_subspace',False) and d['status'] in ['neb_converged','constrained_saddle_connected']
            label='constrained saddle confirmed' if valid else d['status'].replace('_',' ')
            if method=='mace_local_saddle':label+='; water-derived OH + H; N-to-N H transfer; minimization branches, not IRC'
            elif d.get('index_one_in_mobile_subspace') and not valid:label+='; index-one candidate, connectivity unconfirmed'
            barrier=f'{peak:.4f}' if peak is not None else '--'
            if peak is not None and not valid:barrier+=' (unverified peak)'
            rows.append(f"| {folder.parent.parent.name} / {folder.parent.name} | {method} | {barrier} | {label} | [data](../../{rel}/README.md) |")
            text=f"# {folder.parent.parent.name} / {folder.parent.name}: {method}\n\nStatus: **{label}**.\n\n"
            if (folder/'path.gif').exists():text+='![Evaluated reaction images](path.gif)\n\n[Every image energy](energies.csv) | [Actual coordinates](images.extxyz).\n\n'
            text+=f"Peak relative to IS: **{barrier} eV**. A peak without convergence, curvature and connectivity is not a verified transition state.\n\n"
            text+='[Calculation provenance and full checks](summary.json) | [Scope, equations and sources](../../../../../docs/dry_etch_results/MOLECULAR_SURFACES.md).\n'
            if 'capped' in folder.parent.parent.name:text+='\nThis is a capped molecular Si-N/Si-O bonding proxy, not a periodic slab or an etch-rate calculation. SiH3 is fixed.\n'
            else:text+='\nTwo selected substrate atoms and all molecular atoms are mobile; remaining slab atoms are fixed.\n'
            (folder/'README.md').write_text(text,encoding='utf8')
            reaction.append(dict(folder=rel,**d))
    lines+=['','## Molecular reaction searches','','Failed searches are part of the evidence. `no distinct dissociation endpoints` means the initial and final breaking-bond distances did not identify intact versus broken states; it does not establish a barrierless reaction or identical endpoint structures.','', '| System | Method | Peak above IS (eV) | Numerical conclusion | Evidence |','|---|---|---:|---|---|']+rows
    lines+=['','## Direct DFT checks on saved ML geometries','','These compare electronic energies at the same saved geometries. Endpoint-only differences are reaction energies, not barriers. Full-path maxima remain single-point peaks, not DFT-optimized saddles.','','| Local motif | Geometry set | Images | DFT FS - IS (eV) | DFT max - IS (eV) | Evidence |','|---|---|---:|---:|---:|---|']
    dft=[]
    for p in sorted((ROOT/'data/surface_paths').glob('*_capped_motif/HF/dft_*/summary.json')):
        d=json.loads(p.read_text());assert d['status']=='complete';energies=[r['energy_eV'] for r in d['frames']];rel=p.parent.relative_to(ROOT).as_posix();n=len(energies)
        accepted=sum(r.get('accepted_for_comparison',True) for r in d['frames']);peak_text=f'{max(energies)-energies[0]:.4f}' if accepted==n else 'not reported: SCF unresolved'
        lines.append(f"| {p.parents[2].name} | {p.parent.name} | {accepted}/{n} accepted | {energies[-1]-energies[0]:.4f} | {peak_text} | [energies](../../{rel}/energies.csv), [SCF metadata](../../{rel}/summary.json) |")
        dft.append(dict(folder=rel,**d))
    doc=ROOT/'docs/dry_etch_results/MOLECULAR_SURFACES.md';s=doc.read_text().split('<!-- RESULTS -->')[0];doc.write_text(s+'<!-- RESULTS -->\n\n'+'\n'.join(lines)+'\n',encoding='utf8')
    (ROOT/'docs/dry_etch_results/molecular_reactions.json').write_text(json.dumps(reaction,indent=2)+'\n')
    print('Reported',len(results),'surface/species systems;',len(reaction),'reaction attempts;',len(dft),'DFT sets')
if __name__=='__main__':main()
