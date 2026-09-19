"""Refresh compact README status blocks from saved calculation evidence."""
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[1]
def block(text,name,body):
    start='<!-- BEGIN '+name+' -->';end='<!-- END '+name+' -->'
    before,rest=text.split(start,1);_,after=rest.split(end,1)
    return before+start+'\n'+body+'\n'+end+after

def main():
    d=json.loads((ROOT/'docs/species_kmc_results/summary.json').read_text());validation=ROOT/'docs/current_validation.json'
    tests=str(json.loads(validation.read_text())['passed'])+' tests passed' if validation.exists() else 'See recorded test results'
    lines=['| Component | Current status |','|---|---|',f"| Numerical checks | **{tests}**; atom conservation, independent master equation, Python/C++ statistical agreement and artifact provenance. |",f"| Species-resolved kMC | **{d['states']} states / {d['enabled_events']} enabled events**, 16 source pathways; **{len(d['temperatures'])*d['replicates']} C++ trajectories** across 325-450 K. Conditional kinetics, not calibrated ALE. |",'| Network discovery | HiPRGen pilot: **40 forward + 40 reverse candidates**, with missing-intermediate audit; no automatic rate assignment. |','| Atomistic evidence | MACE surface paths and molecular screening; OMol25 local models and direct DFT diagnostics. Convergence limits retained. |']
    p=ROOT/'README.md';s=p.read_text();r=ROOT/'docs/dry_etch_results/relaxed_adsorption.json'
    if r.exists():
        rows=json.loads(r.read_text());n=sum(x['converged'] for x in rows)
        lines.append(f'| Beyond rigid scans | **{n}/{len(rows)} HF/substrate relaxations force-converged**, two starts on each of four surfaces; see reference convergence below. |')
        body=f"**{n}/{len(rows)} adsorbate/substrate relaxations meet 0.04 eV/angstrom.** Lower-half atoms remain fixed. Force convergence does not establish a stable minimum or transition state. [Full table, energy traces and actual relaxation GIFs](docs/dry_etch_results/RELAXED_ADSORPTION.md).\n\n![Force-driven adsorption relaxation](docs/dry_etch_results/relaxed_adsorption.png)"
        s=block(s,'RELAXATION STATUS',body)
    multi=ROOT/'docs/multilayer_results/summary.json'
    if multi.exists():
        m=json.loads(multi.read_text());q=json.loads((ROOT/'data/multilayer/rate_requests/index.json').read_text())['requests']
        lines.append(f"| Multilayer prototype | {m['substrate_atoms']} Si/N atoms, six depth bands; explicit bonds and H/F/Cl termination. Unvalidated demonstration rates; strict mode blocks missing data. |")
        if '<!-- BEGIN MULTILAYER STATUS -->' in s:
            s=block(s,'MULTILAYER STATUS',f"The 9 s demonstration recorded **{m['baseline_events']} events**, **{m['removed'].get('Si',0)} Si + {m['removed'].get('N',0)} N removals**, and **{m['newly_exposed_atoms']} newly exposed atoms**. Twelve access-depth/seed controls accompany it. The rate audit identifies **{len(q)} distinct missing-rate environments**, with separate IS/FS connectivity requests for priority cases.")
    s=block(s,'CURRENT STATUS','\n'.join(lines));p.write_text(s,encoding='utf-8')
if __name__=='__main__':main()
