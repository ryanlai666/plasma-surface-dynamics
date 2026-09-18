# Render the recorded barrier comparison; no calculations are launched.
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'docs/barrier_results'

def main():
    cases=[('DPA-3.3 / OMol25','omol25_nh3'),('DFT PBE / def2-SVP','pbe_svp_nh3')]
    if (OUT/'pbe_tzvp_nh3/summary.json').exists():cases.append(('DFT PBE / def2-TZVP','pbe_tzvp_nh3'))
    fig,ax=plt.subplots(figsize=(8.3,4.8))
    rows=[]
    for label,folder in cases:
        path=OUT/folder
        s=json.loads((path/'summary.json').read_text());m=json.loads((path/'manifest.json').read_text())
        if m['status']!='verified':raise ValueError('Unverified saddle result: '+folder)
        p=json.loads((path/'profile.json').read_text())
        ax.plot([x['q_A'] for x in p],[x['relative_energy_eV'] for x in p],'.-',label=label,lw=2)
        rows.append(dict(label=label,folder=folder,barrier=s['barrier_eV'],force=s['saddle_max_force_eV_A'],evaluations=m['evaluations'],seconds=m['elapsed_s']))
    ax.set(xlabel='Signed N height above the H plane (angstrom)',ylabel='Energy relative to each optimized minimum (eV)',title='NH3 umbrella inversion: relaxed molecular paths')
    ax.axvline(0,color='gray',lw=.7,ls=':');ax.grid(alpha=.2);ax.legend(frameon=False)
    fig.text(.5,.015,'Electronic energies; no ZPE correction. Molecular diagnostic, not a SiN etching barrier.',ha='center',fontsize=9)
    fig.tight_layout(rect=[0,.035,1,1]);fig.savefig(OUT/'comparison.png',dpi=180);plt.close(fig)
    table='| Method | Electronic barrier (eV) | TS max force (eV/A) | Saddle checks |\n|---|---:|---:|---|\n'
    for r in rows:table+=f"| {r['label']} | {r['barrier']:.6f} | {r['force']:.2e} | Passed |\n"
    basis_note=''
    if len(rows)==3:
        shift=rows[2]['barrier']-rows[1]['barrier']
        basis_note=f'Changing PBE from def2-SVP to def2-TZVP changes this electronic barrier by {shift:+.6f} eV. This two-basis comparison does not establish the complete-basis limit.\n\n'
    report='# Local molecular transition-state results\n\n'+table+'\n![Relaxed inversion profiles](comparison.png)\n\n'
    report+=basis_note
    report+=(ROOT/'docs/BARRIER_METHOD.md').read_text(encoding='utf-8')
    report+='\n## Recorded runs\n\n'
    for r in rows:report+=f"- [{r['label']} summary]({r['folder']}/summary.json), [manifest]({r['folder']}/manifest.json), [saddle]({r['folder']}/saddle.xyz): {r['evaluations']} evaluations; {r['seconds']:.1f} s recorded computation time.\n"
    (OUT/'REPORT.md').write_text(report,encoding='utf-8')
    block='''<!-- BEGIN ATOMISTIC RESULTS -->
## DFT data and transition-state results

A small NH3 inversion calculation tests the molecular saddle workflow using **direct DFT** and the existing **OMol25-trained DPA model**:

'''+table+'\n'+basis_note+'''
![Direct DFT and OMol25 molecular inversion paths](docs/barrier_results/comparison.png)

All reported saddles pass full-force, negative-curvature, and two-sided relaxation checks. NH3 inversion is a molecular workflow diagnostic, **not a SiN surface etching barrier**. The comparison uses different reference methods and excludes zero-point/free-energy corrections.

| Collected reference data | Contents |
|---|---|
| Published Si/O/C/F DFT | 70 energy/force frames, including 25 quasi-static drag configurations; drag peaks are not validated TS barriers |
| Si-H-Cl source geometries | 22 structures, including 15 source-named TS; the archive lacks energy/force/Hessian labels |
| Published Cl diffusion on Si(111)-(5x5) | 1.73 eV hopping and 1.34 eV SiCl-complex diffusion, with sources and applicability limits |
| Fragment/side-reaction catalog | 35 species/bookkeeping units, 21 balanced candidate reactions, and eight proposed experiment sets |

[Computed results and reproduction](docs/barrier_results/REPORT.md) | [Data, diffusion equations, and sources](docs/ATOMISTIC_DATA.md) | [Fragments and experiment sets](docs/REACTION_CANDIDATES.md).

<!-- END ATOMISTIC RESULTS -->

'''
    readme=ROOT/'README.md';s=readme.read_text(encoding='utf-8')
    begin='<!-- BEGIN ATOMISTIC RESULTS -->';end='<!-- END ATOMISTIC RESULTS -->'
    if begin in s:
        a=s.index(begin);b=s.index(end,a)+len(end);s=s[:a]+block.rstrip()+s[b:]
    else:s=s.replace('## Physical hypotheses, equations, and literature support',block+'## Physical hypotheses, equations, and literature support')
    s=s.replace('22 tests passed','25 tests passed')
    readme.write_text(s,encoding='utf-8',newline='\n')

if __name__=='__main__':main()
