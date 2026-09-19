"""Report calculated coadsorbate candidates without promoting them to rate evidence."""
from pathlib import Path
import csv, hashlib, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from ase.io import read
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/intermediate_campaign'
OUT=ROOT/'docs/intermediate_results'
COLORS={'Si':'#7C91AD','N':'#3478B9','O':'#D65B4E','F':'#269D78','H':'#E6CD83'}
CUTOFF={('F','H'):1.4,('H','N'):1.3,('H','O'):1.3,('H','Si'):1.8,('F','Si'):2.0,('N','Si'):2.15,('O','Si'):2.0}

def bonds(a):
    distances=a.get_all_distances(mic=True);sym=a.get_chemical_symbols()
    return [[i,j] for i in range(len(a)) for j in range(i) if distances[i,j]<CUTOFF.get(tuple(sorted((sym[i],sym[j]))),0)]

def main():
    OUT.mkdir(exist_ok=True);d=json.loads((DATA/'refined_summary.json').read_text());rr=d['results'];nodes=[];edges=[];table=[]
    fig=plt.figure(figsize=(16,9));energy,axes=plt.subplots(2,4,figsize=(15,7),layout='constrained')
    for k,r in enumerate(rr):
        folder=ROOT/r['folder'];a=read(folder/'final.extxyz');initial=read(ROOT/r['initial_attempt']/'trajectory.extxyz',0);pos=a.positions;sym=a.get_chemical_symbols();n=r['host_atoms'];bb=bonds(a)
        # Distances use MIC for classification. The drawing omits wraparound edges.
        ax=fig.add_subplot(2,4,k+1,projection='3d')
        seg=[pos[[i,j]] for i,j in bb if np.linalg.norm(pos[i]-pos[j])<2.2]
        ax.add_collection3d(Line3DCollection(seg,colors='#8895A3',linewidths=1.2))
        ax.scatter(*pos[:n].T,c=[COLORS[s] for s in sym[:n]],s=40,alpha=.55,depthshade=False)
        ax.scatter(*pos[n:].T,c=[COLORS[s] for s in sym[n:]],s=80,edgecolors='#324154',linewidths=.6,depthshade=False)
        for i in range(n,len(a)):ax.text(*pos[i],str(i),fontsize=7)
        ax.set_box_aspect(np.maximum(np.ptp(pos,axis=0),1));ax.view_init(18,-65);ax.axis('off')
        name=('SiN' if r['surface'].startswith('beta') else 'SiO2')+' / HF + '+r['coadsorbate']+f" / start {r['start']}"
        assoc=r['incremental_association_energy_eV'];value=f'{assoc:+.3f} eV' if assoc is not None else 'not assigned'
        ax.set_title(name+'\nRaw refinement energy: '+value,fontsize=10)
        rows=list(csv.DictReader((folder/'energies.csv').open()));steps=[int(x['optimization_step']) for x in rows];energies=[float(x['energy_relative_to_start_eV']) for x in rows]
        axes.flat[k].plot(steps,energies,c='#269D78');axes.flat[k].set(title=name,xlabel='BFGS step (not reaction coordinate)',ylabel='E - initial E (eV)');axes.flat[k].grid(alpha=.15)
        geometric='parent HF separated' if r['geometry']['parent_HF_distance_A']>1.4 else 'parent HF retained'
        threshold=1.4 if r['coadsorbate']=='HF' else 1.3
        geometric+='; added '+r['coadsorbate']+(' separated' if max(r['geometry']['added_bond_distances_A'])>threshold else ' retained')
        key=r['surface']+'_HF_'+r['coadsorbate']+'_start_'+str(r['start']);before={tuple(x) for x in bonds(initial)};after={tuple(x) for x in bb}
        node=dict(id=key,formula=a.get_chemical_formula(),host_atoms=n,adsorbate_atom_indices=list(range(n,len(a))),symbols=sym,positions_A=pos.tolist(),cell_A=a.cell.array.tolist(),pbc=a.pbc.tolist(),distance_graph=bb,formed_distance_edges=[list(x) for x in sorted(after-before)],lost_distance_edges=[list(x) for x in sorted(before-after)],classification=geometric,converged=r['converged'],curvature_screen=r['adsorbate_curvature'],structure=r['folder']+'/final.extxyz',structure_sha256=r['final_sha256'],rate_enabled=False)
        nodes.append(node);edges.append(dict(source=r['surface']+f"_parent_{r['start']}",gas_reactant=r['coadsorbate'],target=key,kind='coadsorption_candidate',incremental_association_energy_eV=assoc,barrier_eV=None,rate_s=None,rate_enabled=False))
        curve=r['adsorbate_curvature'];curv=f"{curve['minimum_eigenvalue_eV_A2']:+.3f}" if curve else 'not evaluated'
        table.append(f"| {name} | {'yes' if r['converged'] else 'no'} | {r['max_mobile_force_eV_A']:.3f} | {value} | {curv} | [structure](../../{r['folder']}/final.extxyz), [energies](../../{r['folder']}/energies.csv) |")
    fig.suptitle('Relaxed coadsorbate candidates: actual final MACE geometries',fontsize=17)
    fig.text(.5,.035,'Si gray-blue | N blue | O red | F green | H gold. Adsorbate atom indices shown. Faint atoms: substrate.\nDistance-cutoff bonds are a geometric diagnostic; periodic wrap bonds omitted visually. No transition states are claimed.',ha='center',fontsize=10)
    fig.subplots_adjust(top=.88,bottom=.11,wspace=.05,hspace=.16);fig.savefig(OUT/'coadsorbate_structures.png',dpi=150);plt.close(fig)
    energy.suptitle('Evaluated optimization energies; these are not minimum-energy reaction paths',fontsize=14);energy.savefig(OUT/'relaxation_energies.png',dpi=140);plt.close(energy)
    qualified=json.loads((DATA/'stability_qualification.json').read_text())
    lookup={(r['surface'],r['coadsorbate'],r['start']):r for r in qualified['results']}
    for node,edge,r in zip(nodes,edges,rr):
        q=lookup[(r['surface'],r['coadsorbate'],r['start'])]
        node['full_stability_screen']=q
        edge['incremental_energy_scope']='Raw original refinement; consult full-stability screened field'
        edge['screened_association_energy_eV']=q['screened_association_energy_eV']
        edge['screened_energy_structure']=q['structure']
        edge['DFT_validated']=False
    network=dict(scope='Eight calculated candidate endpoints and coadsorption bookkeeping edges; no kinetic network extension enabled',nodes=nodes,edges=edges,distance_cutoffs_A={','.join(k):v for k,v in CUTOFF.items()},runner='scripts/report_coadsorbate_intermediates.py',runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),summary_sha256=hashlib.sha256((DATA/'refined_summary.json').read_bytes()).hexdigest())
    (DATA/'candidate_network.json').write_text(json.dumps(network,indent=2)+'\n')
    count=sum(r['converged'] for r in rr)
    text=f'''# Intermediate-gap campaign: HF coadsorbates

**{count}/{len(rr)} local relaxations force-converged** at 0.02 eV/angstrom. These are MACE candidates addressing A1/A5, not DFT-validated intermediates, transition states or new kMC rates. Two inherited single-HF adsorption starts on each ideal nitride/oxide slab were extended with HF or H2O. All eight starts are retained, including failures.

![Calculated candidate structures](coadsorbate_structures.png)

## Original refinement results

**Use the [full mobile-coordinate follow-up](FULL_STABILITY.md) for current reference eligibility.** It includes substrate-coupled modes, additional escapes and explicitly excludes unstable parent references. The figure and table below preserve the original refinement stage, not the final stability-screened geometries.

| Surface / coadsorbate / start | Force converged | Max mobile force (eV/A) | Incremental association | Lowest adsorbate-block curvature (eV/A2) | Files |
|---|---|---:|---:|---:|---|
'''+ '\n'.join(table)+'''

The energy reference is

`Delta E_add = E_relaxed(slab + HF + X) - E_relaxed(slab + HF) - E_relaxed(X_gas)`, with `X = HF or H2O`.

Each parent was re-relaxed using the same fixed lower substrate and model, then BFGS-refined to the same 0.02 eV/angstrom tolerance as the coadsorbates. Gas references use the same model and a stricter 0.01 eV/angstrom tolerance. An association energy is assigned only if all three optimizations converged. It includes any rearrangement of the original HF and mobile substrate. It is neither a cleavage barrier nor a finite-temperature adsorption free energy. A negative value alone does not establish a stable adsorbed intermediate.

Finite differences of all adsorbate forces screen the adsorbate-only Cartesian Hessian block at 0.01 angstrom displacement. Values below -0.02 eV/A2 flag appreciable negative curvature; smaller values can still represent soft instabilities or numerical error. A nonnegative block does not exclude an instability involving mobile substrate atoms. These unweighted curvatures are not vibrational frequencies. Full host-coupled Hessians are now available in the linked follow-up; cell-size checks and DFT comparison remain required.

The [candidate network](../../data/intermediate_campaign/candidate_network.json) records atom indices, actual coordinates, distance-based connectivity and changes relative to each start. Distance cutoffs are diagnostics, not electronic bond orders. No bond-cutoff graph establishes a chemical TS or an elementary mechanism.

![Actual relaxation energies](relaxation_energies.png)

These curves use evaluated geometries at optimization steps. Their maxima are not TS peaks; optimization trajectories are not physical time or minimum-energy reaction paths.

## Why these gaps matter: independent evidence

- [Jung et al. (2020)](https://doi.org/10.1116/1.5125569) resolve HF/HF and HF/H2O coadsorbed states on fluorinated clusters. Their ground-state barriers must be distinguished from a model adjustment for vibrationally excited HF. Sticking and desorption prefactors include fitted quantities. The [six-row extraction](../../data/intermediate_campaign/literature_coadsorbate_parameters.csv) retains source location, reference environment and transfer restrictions. Kelvin barrier parameters are converted using `E_eV = (Ea/R)_K * kB_eV/K`; no values enter our kMC library.
- [Lill et al. (2024)](https://doi.org/10.1116/6.0004019) report water-enhanced thermal HF etching of SiN and weaker AFS binding with water substitution. This motivates two competing explanations: water-assisted local reaction and water-assisted removal of retained salt. Our clean-slab coadsorbates test neither explanation completely.
- [Khumaini et al. (2024)](https://doi.org/10.1016/j.apsusc.2024.159414) use hydrogenated amorphous nitride and consider salt formation. Their material state makes transfer to an ideal crystalline slab uncertain; hydrogen content and retained products must be varied explicitly.

[Source-access records](../../data/intermediate_campaign/sources.json). The latter two source assessments use publicly accessible abstracts; no unavailable supplementary geometries were reconstructed or assumed verified.

## What is closed, and what remains open

| Gap | Progress in this campaign | Required next evidence |
|---|---|---|
| A1: precursor library | Explicit two-molecule coordinates and relaxation records on two substrates | Full minimum tests, termination and coverage ensembles, adsorption/desorption kinetics |
| A5: coadsorbates | HF/HF and HF/H2O candidates with reference-consistent incremental energies | Matched stepwise/concerted proton-transfer paths and connected TS searches |
| A2/A3: sequential/final cleavage | Existing environment-specific requests remain available | Relaxed embedded IS/FS for NH, NH2 and bare-N backbonds; qualified barriers; do not borrow one stage's barrier |
| A6: retention | Salt-removal hypothesis now distinguished from coadsorption hypothesis | Surface AFS/NH4F structures, water substitution, decomposition/desorption free energies |
| B3: material realism | Crystalline nitride/oxide comparison only | Amorphous SiNx:H ensembles, defects, multiple terminations and density/composition controls |

No gap is declared fully closed. The current kMC uses a single HF occupancy flag and cannot represent two coadsorbates faithfully. Additional occupancy states, atom accounting and qualified rates are required before these nodes can become enabled events. The [246-environment calculation queue](../../data/multilayer/rate_requests/index.json) remains separate from this screening library; these new geometries are not matches to those embedded graph environments.

## Reproduce

Use the repository MACE/ASE environment with the pinned local model:

```sh
python scripts/relax_coadsorbate_intermediates.py
python scripts/refine_coadsorbate_intermediates.py
python scripts/validate_coadsorbate_minima.py
python scripts/report_followup_validation.py --section coadsorbates
python scripts/report_coadsorbate_intermediates.py
```

[Run metadata and hashes](../../data/intermediate_campaign/refined_summary.json). Each system has its own folder with the original FIRE attempt and a separate BFGS refinement, evaluated trajectories, energy tables, final coordinates, optimizer logs and summaries. The plotted energies are the BFGS continuation, not the original starting geometry. Initial FIRE files inherit old single-HF scan metadata; use the correctly labeled refinement exports and JSON metadata for this campaign. Gas and single-HF references are retained separately. This is local CPU screening on a general materials potential, with no vibrational excitation, charge-state control, entropy or electronic-structure validation.
'''
    (OUT/'REPORT.md').write_text(text,encoding='utf-8');print('Reported',len(rr),'coadsorbate candidates;',count,'force-converged')

if __name__=='__main__':main()
