"""Report relaxed geometry evidence without interpreting minimization as a TS path."""
from pathlib import Path
import json,csv,hashlib,io
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
from ase.io import read
from render_molecular_campaign import COLORS
ROOT=Path(__file__).resolve().parents[1]
def main():
    out=ROOT/'docs/dry_etch_results';rows=json.loads((out/'relaxed_adsorption.json').read_text())
    for r in rows:
        for key in ['source','folder']:r[key]=r[key].replace(chr(92),'/')
    surfaces=list(dict.fromkeys(r['surface'] for r in rows))
    fig,axes=plt.subplots(2,2,figsize=(10,7),layout='constrained');table=[]
    for ax,surface in zip(axes.flat,surfaces):
        for r in [r for r in rows if r['surface']==surface]:
            folder=ROOT/r['folder'];data=list(csv.DictReader((folder/'energies.csv').open()));steps=[int(x['step']) for x in data];es=np.array([float(x['energy_eV']) for x in data]);ax.plot(steps,es-es[0],label=f"Start {r['start']}: "+('converged' if r['converged'] else 'not converged'))
            reference=json.loads((folder.parent/'summary.json').read_text())['bare_slab'];ads='not reported' if r['adsorption_energy_eV'] is None else f"{r['adsorption_energy_eV']:.3f}"
            table.append(f"| {surface} / {r['start']} | {r['site']} / {r['orientation']} | {'yes' if r['converged'] else 'no'} / {'yes' if reference['converged'] else 'no'} | {r['max_mobile_force_eV_A']:.4f} | {r['relaxation_energy_eV']:.3f} | {ads} | {r['HF_distance_A']:.3f} | [data + GIF](../../{r['folder']}/README.md) |")
            if not (folder/'relaxation.gif').exists():
                images=read(folder/'trajectory.extxyz',':');indices=np.unique(np.linspace(0,len(images)-1,min(12,len(images)),dtype=int));frames=[]
                xyz=np.array([a.positions for a in images]);lo=xyz.min(axis=(0,1));hi=xyz.max(axis=(0,1));mid=(hi+lo)/2;span=max(hi-lo)
                for i in indices:
                    a=images[i];f=plt.figure(figsize=(9,4.5),dpi=100);view=f.add_subplot(121,projection='3d');ep=f.add_subplot(122)
                    view.scatter(*a.positions.T,s=[35 if x=='H' else 75 for x in a.get_chemical_symbols()],c=[COLORS[x] for x in a.get_chemical_symbols()],edgecolors='#394955',linewidths=.4)
                    for setter,c in zip([view.set_xlim,view.set_ylim,view.set_zlim],mid):setter(c-span*.55,c+span*.55)
                    view.set_box_aspect([1,1,1]);view.view_init(25,-55);view.axis('off');view.set_title(surface+' / HF')
                    ep.plot(steps,es-es[0],color='#288A91');ep.scatter(i,es[i]-es[0],color='#D15A30',s=60);ep.set(xlabel='FIRE optimization step',ylabel='Energy relative to start (eV)',title=f'Actual evaluated step {i}');ep.grid(alpha=.2)
                    f.suptitle('Upper substrate + HF relaxation; lower substrate fixed',fontsize=13);f.text(.5,.03,'Optimization trajectory, not a reaction coordinate or transition-state path. No interpolated frames.',ha='center',fontsize=9);f.subplots_adjust(bottom=.19,top=.81,wspace=.05)
                    b=io.BytesIO();f.savefig(b,format='png');b.seek(0);frames.append(Image.open(b).convert('RGB'));plt.close(f)
                frames[0].save(folder/'relaxation.gif',save_all=True,append_images=frames[1:],duration=400,loop=0)
                (folder/'render_manifest.json').write_text(json.dumps(dict(frame_indices=indices.tolist(),interpolated_frames=0,trajectory_sha256=hashlib.sha256((folder/'trajectory.extxyz').read_bytes()).hexdigest(),gif_sha256=hashlib.sha256((folder/'relaxation.gif').read_bytes()).hexdigest()),indent=2)+'\n')
            (folder/'README.md').write_text(f"# {surface} / HF - relaxed start {r['start']}\n\n![Actual relaxation frames](relaxation.gif)\n\nForce-converged: **{r['converged']}**. {r['geometry_classification']}; final HF separation {r['HF_distance_A']:.3f} angstrom. This distance classification is not a full bonding or stability analysis.\n\nUpper-substrate maximum displacement: {r['substrate_max_displacement_A']:.3f} angstrom. Mobile substrate atoms: {r['mobile_substrate_atoms']}; all HF atoms mobile.\n\n[Evaluated trajectory](trajectory.extxyz), [energies](energies.csv), [metadata](summary.json), [frame mapping](render_manifest.json).\n\nNo Hessian, TS, DFT or global-minimum validation is implied. [Campaign assumptions](../../../../../../docs/dry_etch_results/RELAXED_ADSORPTION.md).\n")
        ax.set(title=surface,xlabel='FIRE optimization step',ylabel='Energy relative to own start (eV)');ax.legend(fontsize=8);ax.grid(alpha=.2)
    fig.suptitle('Beyond rigid scans: HF and the upper substrate move',fontsize=15);fig.savefig(out/'relaxed_adsorption.png',dpi=160);plt.close(fig)
    text="""# Relaxed HF adsorption: four surfaces, two distinct starts each

This is a force-driven follow-up to rigid approach screening. HF and the **upper half of the substrate atoms by height** can move; atoms at or below the median substrate height remain fixed. The bare slab uses exactly the same fixed indices. Cell vectors stay fixed. MACE-MP-0b2 small and FIRE are used with a 0.04 eV/angstrom force criterion and 300-step budget.

Starts are selected from evaluated scan heights >=1.9 angstrom: first the lowest interaction-energy point, then the lowest point with both a different site and a different orientation. Thus two different starts are actually tested, but this does not exhaust orientations, reconstructions or amorphous environments.

![Energy during relaxation](relaxed_adsorption.png)

| Surface / start | Initial site / orientation | Combined / bare slab converged | Final mobile fmax (eV/A) | Relaxation energy (eV) | Apparent adsorption difference (eV) | HF distance (A) | Artifacts |
|---|---|---|---:|---:|---:|---:|---|
"""+'\n'.join(table)+"""

`Relaxation energy = E_final - E_start` follows one composition and one starting geometry. It can include substrate reconstruction and does not measure an activation barrier. FIRE need not decrease the energy monotonically at every step.

`Apparent adsorption difference = E_relaxed(slab+HF) - E_relaxed(slab) - E_isolated(HF)` is tabulated only if the combined and bare-slab optimizations both meet the force criterion. The gas geometry is the previously optimized isolated HF using the same model. The combined system has finite periodic coverage. A dissociated HF geometry still uses intact HF as the reference. Reconstruction into different local basins can influence this difference; it is not a converged experimental adsorption enthalpy.

HF distance >1.4 angstrom is only a screening label for bond separation. See the reference decomposition below before interpreting this number as HF binding. Stable minima require Hessian checks; elementary reactions need connected saddle searches. Force convergence alone supplies neither. No rates were changed using these results.

Remaining requirements: realistic termination and amorphous SiNx:H structures, coverage and cell/thickness convergence, more starts, curvature and connectivity checks, then matched DFT and experimental validation. The calculations retain ideal unpassivated surfaces and the MACE domain limitations described in the [molecular campaign](MOLECULAR_SURFACES.md).

Reproduce in the optional MACE environment:

```sh
python scripts/relax_adsorption_multistart.py
python scripts/report_relaxed_adsorption.py
python scripts/update_research_readme.py
```

Outputs are per surface/species under `data/surface_paths/<surface>/HF/relaxed_adsorption/`. Every GIF uses sampled, evaluated optimization steps; none uses interpolated geometries.
"""
    audit=json.loads((out/'adsorption_reference_audit.json').read_text())
    text+='\n## Reference-state audit: separate binding from reconstruction\n\nThe apparent adsorption difference decomposes exactly as `E_interaction(frozen fragments) + [E_slab(final geometry) - E_slab(reference)] + [E_HF(periodic final geometry) - E_HF(isolated reference)]`. The molecular term includes both distortion and periodic-layer reference effects. A negative slab term demonstrates that the original clean reference was not the lowest accessible geometry, even if its force criterion passed.\n\n| Surface / start | Frozen-fragment interaction (eV) | Slab-reference shift (eV) | Molecular reference shift (eV) | Lower clean-slab geometry found? |\n|---|---:|---:|---:|---|\n'
    for r in audit:text+=f"| {r['surface']} / {r['start']} | {r['frozen_fragment_interaction_eV']:.3f} | {r['substrate_reference_shift_eV']:.3f} | {r['molecular_layer_reference_shift_eV']:.3f} | {'yes' if r['lower_clean_slab_geometry_found'] else 'not flagged'} |\n"
    text+='\nThe -0.1 eV flag is a diagnostic threshold, not a validated accuracy tolerance. Flagged apparent differences are excluded from quantitative adsorption claims and kinetic parameterization. Remedy: independently reconstruct/passivate bare slabs, perform multistart reference searches, then compare identical structural basins and DFT forces. Merely allowing more atoms to move does not solve reference-state bias.\n\n[All numerical decompositions](adsorption_reference_audit.json). Run `python scripts/audit_adsorption_references.py` before regenerating this report.\n'
    (out/'RELAXED_ADSORPTION.md').write_text(text,encoding='utf-8');print('Reported and rendered',len(rows),'relaxed starts')
if __name__=='__main__':main()
