"""Report completed stability and DFT checks without assigning kinetic rates."""
from pathlib import Path
import json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]

def stability():
    source=ROOT/'data/intermediate_campaign/full_stability.json';d=json.loads(source.read_text())
    if d['status']!='complete':raise RuntimeError('Stability checks still running')
    rows=[];lines=[]
    for r in d['results']:
        p=d['parents'][r['parent_key']];ok=all(x['minimum_eigenvalue_eV_A2']>=-.02 and x['max_mobile_force_eV_A']<=.02 for x in (p,r))
        rows.append(dict(surface=r['surface'],coadsorbate=r['coadsorbate'],start=r['start'],parent_minimum_curvature_eV_A2=p['minimum_eigenvalue_eV_A2'],candidate_minimum_curvature_eV_A2=r['minimum_eigenvalue_eV_A2'],candidate_force_eV_A=r['max_mobile_force_eV_A'],raw_energy_difference_eV=r['incremental_association_energy_eV'],passes_force_and_appreciable_curvature_screen=ok,screened_association_energy_eV=r['incremental_association_energy_eV'] if ok else None,rate_enabled=False,structure=r['final_structure']))
        lines.append(f"| {r['surface']} / {r['coadsorbate']} / {r['start']} | {p['minimum_eigenvalue_eV_A2']:+.4f} | {r['minimum_eigenvalue_eV_A2']:+.4f} | {r['max_mobile_force_eV_A']:.4f} | {r['incremental_association_energy_eV']:+.3f} | {'pass' if ok else 'exclude'} |")
    (source.parent/'stability_qualification.json').write_text(json.dumps(dict(source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),threshold_eV_A2=-.02,force_tolerance_eV_A=.02,results=rows),indent=2)+'\n')
    text='''# Full mobile-coordinate stability checks

All four single-HF references and eight coadsorbates were checked using every mobile host and adsorbate coordinate. Lower-host atoms remain fixed. Central differences use 0.01 angstrom. Negative-mode displacements and independent BFGS escapes are retained, including failures.

| Surface / added molecule / start | Parent minimum curvature (eV/A2) | Candidate minimum curvature (eV/A2) | Candidate force (eV/A) | Raw energy difference (eV) | Screen |
|---|---:|---:|---:|---:|---|
'''+ '\n'.join(lines)+'''

A pass requires both parent and candidate forces <=0.02 eV/A and no eigenvalue below -0.02 eV/A2. Small negative modes are not automatically physical zero modes. This is provisional MACE screening, not proof of a DFT minimum. Unweighted Hessian curvatures are not vibrational frequencies.

The first nitride parent has appreciable negative curvature despite a small force. Its two coadsorbate energy differences retain an unstable reference and are **excluded from screened association energies**. The JSON stores null for excluded values while preserving raw diagnostics. Force convergence alone is insufficient.

The oxide/water start-1 negative mode was explicitly followed. Compare the original and follow-up structures to inspect changes under identical constraints. No activation barrier or kMC rate is inferred.

[Full Hessians and escape metadata](../../data/intermediate_campaign/full_stability.json) | [Eligibility records](../../data/intermediate_campaign/stability_qualification.json). Each system has a separate `refinement/full_stability` folder. Reproduce the calculator using `scripts/validate_coadsorbate_minima.py` in a fresh output location, then this reporter.
'''
    (ROOT/'docs/intermediate_results/FULL_STABILITY.md').write_text(text,encoding='utf-8')

def proxy():
    base=ROOT/'data/final_cleavage/SiF3_NH2_HF';d=json.loads((base/'dft_initial_path/summary.json').read_text());scan=json.loads((base/'fixed_frame_three_coordinate_scan/summary.json').read_text())
    if d.get('status')!='complete' or scan['status']!='complete':raise RuntimeError('Proxy checks still running')
    initial=json.loads((base/'omol25/summary.json').read_text());continued=json.loads((base/'omol25_bfgs/summary.json').read_text());consistency=json.loads((base/'force_consistency.json').read_text());out=ROOT/'docs/final_cleavage_results';out.mkdir(exist_ok=True)
    fig,ax=plt.subplots(1,3,figsize=(15,4.5),layout='constrained');x=np.arange(3);ax[0].plot(x,d['results'][0]['relative_ML_energies_eV'],'o--',label='OMol25');table=[]
    for r in d['results']:
        ax[0].plot(x,r['relative_DFT_energies_eV'],'o-',label='PBE/'+r['basis']);ax[1].plot(x,[f['max_mobile_DFT_force_eV_A'] for f in r['frames']],'o-',label=r['basis'])
        table.append(f"| PBE/{r['basis']} | {r['relative_DFT_energies_eV'][1]:.3f} | {r['relative_DFT_energies_eV'][2]:+.3f} | "+' / '.join(f"{f['max_mobile_DFT_force_eV_A']:.3f}" for f in r['frames'])+' |')
    for a in ax[:2]:a.set_xticks(x,['ML reactant','Failed-band peak','ML product']);a.tick_params(axis='x',labelsize=8);a.legend(fontsize=8);a.grid(alpha=.2)
    ax[0].set(ylabel='Energy relative to reactant (eV)',title='DFT checks at identical geometries');ax[1].set(ylabel='Maximum mobile force (eV/A)',title='DFT stationarity is not satisfied')
    for direction in ('forward','reverse'):
        rr=[r for r in scan['results'] if r['direction']==direction];ax[2].plot([r['SiN_A'] for r in rr],[r['energy_eV']-scan['results'][0]['energy_eV'] for r in rr],'o-',label=direction)
        bad=[r for r in rr if not r['converged']];ax[2].scatter([r['SiN_A'] for r in bad],[r['energy_eV']-scan['results'][0]['energy_eV'] for r in bad],marker='x',s=100,c='red')
    ax[2].set(xlabel='Constrained Si-N distance (A)',ylabel='Energy relative to initial geometry (eV)',title='Constrained scan; not a TS path');ax[2].legend();ax[2].grid(alpha=.2);fig.savefig(out/'dft_and_path_checks.png',dpi=150);plt.close(fig)
    maxerr=max(r['max_absolute_error_eV_A'] for r in consistency['results']);converged=sum(r['converged'] for r in scan['results']);end=next(r for r in scan['results'] if r['direction']=='reverse' and r['index']==10)
    text='''# Final Si-N cleavage: path searches and DFT diagnostics

**Target:** SiF3NH2 + HF -> SiF4 + NH3. NH3 is the nitrogen-containing removal product. This neutral-singlet molecular proxy has a fixed SiF3 frame; it is not a periodic or embedded surface and does not represent bare-N backbonds.

![Evaluated energies, DFT forces and scan hysteresis](dft_and_path_checks.png)

| Evaluation at original OMol25 geometries | Failed-band peak minus reactant (eV) | Product minus reactant (eV) | DFT maximum mobile forces: reactant / peak / product (eV/A) |
|---|---:|---:|---|
'''+ '\n'.join(table)+f'''

These fixed-geometry differences are **not activation barriers**. Density-fitted PBE uses grid level 3 and SCF tolerance 1e-9 hartree. All retained SCFs converged and passed internal orbital-stability checks; external spin stability was not tested. Nonzero DFT endpoint forces prevent a stationary-point or activation-free-energy claim.

| Path search | Peak above reactant (eV) | Maximum NEB force (eV/A) | Outcome |
|---|---:|---:|---|
| Initial FIRE band | {initial['peak_above_IS_eV']:.3f} | {initial['max_neb_force_eV_A']:.3f} | Unconverged; not a barrier |
| BFGS continuation | {continued['peak_above_IS_eV']:.3f} | {continued['max_neb_force_eV_A']:.3f} | Unconverged; not a barrier |

Neither peak has exactly one appreciable negative mobile-coordinate curvature. More iterations did not resolve the failure. Finite-difference energy derivatives agree with calculator forces to maximum discrepancy {maxerr:.4f} eV/A across the three original geometries and two displacements; consistency is not DFT accuracy.

## Changed search strategy and fixed-frame correction

The distant original NH3 endpoint mixes local chemistry with product separation. An initial two-coordinate scan also exposed a constraint-ordering defect: the bond projection moved nominally fixed Si by as much as 1.748 A. That run is explicitly [invalidated](../../data/final_cleavage/SiF3_NH2_HF/local_coordinate_scan/INVALIDATION.json) and excluded from this plot and all fixed-frame energy comparisons.

The corrected calculation enforces the fixed frame and three distances jointly: Si-N elongation, incoming H-to-N approach and incoming Si-F approach. Every step checks the original fixed coordinates. Analytic force projection is tested for tangency and virtual-work preservation. No atom masses are altered.

Every saved energy is evaluated with OMol25; none is interpolated. {converged}/22 points pass the projected-force tolerance of 0.04 eV/A. Red crosses mark failures. Physical forces before projection are retained separately. The local product-side Si-F target is {end['SiF_A']:.3f} A. This scan supplies candidate geometries and a hysteresis diagnostic, not a minimum-energy path or TS. Its imposed coordinates can conceal unstable directions and force a chosen mechanism. Unconstrained relaxation, saddle curvature and downhill connectivity remain required.

## Consequence for kMC

A3 remains open. These calculations do not repair the 13 bare-N and 8 NH2 final-cleavage blockers in the 12-cycle graph run. All new rates remain disabled. Even a qualified molecular-proxy saddle would require matched embedded-environment validation before graph transfer.

[Original atom map/path](../../data/final_cleavage/SiF3_NH2_HF/omol25/summary.json) | [Continuation](../../data/final_cleavage/SiF3_NH2_HF/omol25_bfgs/summary.json) | [DFT comparison](../../data/final_cleavage/SiF3_NH2_HF/dft_initial_path/summary.json) | [Coordinate scan](../../data/final_cleavage/SiF3_NH2_HF/fixed_frame_three_coordinate_scan/summary.json) | [Force consistency](../../data/final_cleavage/SiF3_NH2_HF/force_consistency.json).

Each method folder retains geometries, forces, optimizer logs, energies and provenance hashes. No trajectory here is labeled an IS-TS-FS animation because no TS was validated.
'''
    (out/'REPORT.md').write_text(text,encoding='utf-8')

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--section',choices=['all','coadsorbates','proxy'],default='all');args=parser.parse_args()
    if args.section in ('all','coadsorbates'):stability()
    if args.section in ('all','proxy'):proxy()
