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
    section=['<!-- BEGIN MOLECULAR CAMPAIGN -->','### Molecules, surface orientations and adsorption sites','',f"**Nine molecules x four surfaces; {sum(r['combinations'] for r in results)} site/orientation curves and {sum(r['evaluations'] for r in results):,} evaluated geometries**, plus 324 baseline approach geometries. HF, HCl, F2, Cl2, H2, H2O, CH3F, SiF4 and SiCl4 each have separate surface/species folders.",'','| Surface orientation | Sites sampled | Molecular orientations |','|---|---|---|','| Si(100), Si(111) | atop Si, bridge, hollow | upright, parallel, flipped when distinct |','| beta-Si3N4(001) | Si, N, Si-N bridge projections | upright, parallel, flipped when distinct |','| alpha-quartz(001) | Si, O, Si-O bridge projections | upright, parallel, flipped when distinct |','','![Site and orientation energy comparison](docs/dry_etch_results/molecular_screening.png)','','These **rigid approach curves are not transition-state paths**. Each 3D GIF shows only calculated geometries with their energies. The table compares sampled interaction minima, not barriers or etch selectivity. Surfaces are ideal unpassivated cuts; coverage, reconstruction and amorphous composition effects remain unresolved.','','<details>','<summary>3D molecular approach examples and all 36 surface/species folders</summary>','']
    for surface,species in [('Si100','HCl'),('Si111','Cl2'),('beta_Si3N4_001','CH3F'),('alpha_quartz_001','HF')]:
        r=next(r for r in results if r['surface']==surface and r['species']==species);section.append(f"![{species} / {surface}: evaluated approach]({r['folder']}/path.gif)\n")
    section+=['[All curves, folders, energies and sources](docs/dry_etch_results/MOLECULAR_SURFACES.md#calculated-screening-results).','','</details>','','Molecular dissociation and local Si-N/Si-O cleavage are evaluated separately. Failed endpoint/NEB searches are retained with their convergence and Hessian checks. The OMol25 capped motifs and direct PBE/def2-SVP checks are **local molecular models**, not periodic-surface DFT validation.','','[Reaction search table, energy peaks, equations and hypotheses](docs/dry_etch_results/MOLECULAR_SURFACES.md#molecular-reaction-searches) | [Direct DFT results](docs/dry_etch_results/MOLECULAR_SURFACES.md#direct-dft-checks-on-saved-ml-geometries).','<!-- END MOLECULAR CAMPAIGN -->','']
    evidence=['','#### Molecular reaction-path results','','| System / event | ML peak above local IS (eV) | Interpretation |','|---|---:|---|']
    for d in reaction:
        method=d['folder'].split('/')[-1]
        if method not in ['mace_local_saddle','molecular_refined','omol25_refined']:continue
        peak=d.get('peak_above_IS_eV');confirmed=d.get('endpoint_connectivity_confirmed',False)
        label='Constrained saddle connected to checked local minima' if confirmed else d['status'].replace('_',' ')+'; TS not verified'
        if method=='mace_local_saddle':label+='; water-derived OH + H; N-to-N H transfer'
        system=d['folder'].split('/')[2]+' / '+d['folder'].split('/')[3]
        evidence.append(f"| [{system}]({d['folder']}/README.md) | {peak:.4f} | {label} |" if peak is not None else f"| {system} | -- | {label} |")
    evidence+=['','These values are model potential energies with fixed substrate/cap atoms. Capped Si-N/Si-O motifs are not full surfaces. A converged highest image alone does not establish a transition state.','','| Direct PBE/def2-SVP check | DFT FS - IS (eV) | DFT highest saved image - IS (eV) |','|---|---:|---:|']
    for d in dft:
        if not d['folder'].endswith('/dft_path'):continue
        es=[r['energy_eV'] for r in d['frames']];name=d['folder'].split('/')[2];accepted=sum(r.get('accepted_for_comparison',True) for r in d['frames']);peak_text=f'{max(es)-es[0]:.4f}' if accepted==len(es) else 'SCF unresolved; peak rejected'
        evidence.append(f"| [{name}: {len(es)} saved ML path geometries]({d['folder']}/energies.csv) | {es[-1]-es[0]:.4f} | {peak_text} |")
    evidence+=['','DFT values are single-point energies along the saved ML paths, not DFT-optimized transition states. See the full report for path convergence and endpoint checks.','']
    for d in reaction:
        if d['folder'].endswith('/mace_local_saddle') and (ROOT/d['folder']/'path.gif').exists():evidence.append(f"![Evaluated surface H-transfer saddle beside Si-OH]({d['folder']}/path.gif)\n")
    section[-2:-2]=evidence
    readme=ROOT/'README.md';raw=readme.read_bytes();newline='\r\n' if b'\r\n' in raw else '\n';s=raw.decode('utf8').replace('\r\n','\n')
    if '<!-- BEGIN MOLECULAR CAMPAIGN -->' in s:
        before,rest=s.split('<!-- BEGIN MOLECULAR CAMPAIGN -->',1);_,after=rest.split('<!-- END MOLECULAR CAMPAIGN -->',1);s=before+'\n'.join(section)+after
    else:s=s.replace('### Literature-parameterized dry-etch kinetics','\n'.join(section)+'\n### Literature-parameterized dry-etch kinetics')
    s=s.replace('Each new GIF shows **force-optimized, energy-evaluated NEB images**','Each migration GIF below shows **force-optimized, energy-evaluated NEB images**')
    readme.write_bytes(s.replace('\n',newline).encode('utf8'))
    (ROOT/'docs/dry_etch_results/molecular_reactions.json').write_text(json.dumps(reaction,indent=2)+'\n')
    print('Reported',len(results),'surface/species systems;',len(reaction),'reaction attempts;',len(dft),'DFT sets')
if __name__=='__main__':main()
