"""Inventory current kinetic inputs and reproduce conditional source rates."""

from pathlib import Path
from dataclasses import asdict
import csv, hashlib, inspect, json, sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from plasma_surface.model import Parameters, ale_recipe
from plasma_surface.species_kmc import rate_constants
from plasma_surface.thermal_kinetics import arrhenius, KB_EV, md_count_rate


def write_csv(path, rows):
    with path.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main():
    out = ROOT / 'data/kinetic_audit'
    out.mkdir(exist_ok=True)
    records = []

    def add(model, name, value, unit, status, source, meaning):
        records.append(
            dict(
                model=model,
                parameter=name,
                value=value,
                unit=unit,
                status=status,
                source=source,
                meaning=meaning,
            )
        )

    generic = asdict(Parameters())
    units = {
        'site_density_m2': 'm^-2',
        'layer_nm': 'nm',
        'sticking': '1',
        'desorption_prefactor_s': 's^-1',
        'desorption_barrier_ev': 'eV',
        'chemical_threshold_ev': 'eV',
        'physical_threshold_ev': 'eV',
        'chemical_yield_scale': '1',
        'physical_yield_scale': '1',
        'deposition_sticking': '1',
    }
    for k, v in generic.items():
        add(
            'generic_two_state',
            k,
            v,
            units[k],
            'illustrative_unfitted',
            'plasma_surface/model.py',
            'Generic baseline only; not graph chemistry',
        )
    for phase in ale_recipe():
        for k, v in asdict(phase).items():
            if k != 'name':
                add(
                    'generic_recipe_' + phase.name,
                    k,
                    v,
                    'K'
                    if k == 'temperature_k'
                    else 'eV'
                    if k == 'ion_energy_ev'
                    else 's'
                    if k == 'duration_s'
                    else 'm^-2 s^-1',
                    'assumed_recipe',
                    'plasma_surface/model.py',
                    'Default ale_recipe; purge appears twice',
                )
    cards = json.loads((ROOT / 'configs/materials.json').read_text())
    for card in cards['cards']:
        for k, v in {
            **{k: card[k] for k in ('n_si_ratio', 'hydrogen_atomic_fraction')},
            **card['parameters'],
        }.items():
            add(
                card['name'],
                k,
                v,
                units.get(k, '1'),
                'invented_composition_scenario',
                'configs/materials.json',
                'Not a fitted composition law',
            )
    network = json.loads((ROOT / 'configs/species_kmc_network.json').read_text())
    summary = json.loads((ROOT / 'docs/species_kmc_results/summary.json').read_text())
    temps = [r['temperature_K'] for r in summary['temperatures']]
    nu = summary['assumed_prefactor_s']
    arrival = summary['assumed_HF_arrival_s']
    des = summary['assumed_complex_desorption_s']
    for k, v, u in [
        ('prefactor', nu, 's^-1'),
        ('HF_arrival', arrival, 's^-1'),
        ('HF_desorption', des, 's^-1'),
        ('exposure', summary['exposure_s'], 's'),
        ('purge', summary['purge_s'], 's'),
    ]:
        add(
            '45_state',
            k,
            v,
            u,
            'assumed',
            'docs/species_kmc_results/summary.json',
            'No entropy or pressure-derived correction applied',
        )
    for k, v in network['initial_motif_fractions'].items():
        add(
            '45_state',
            'initial_fraction_' + k,
            v,
            '1',
            'assumed',
            'configs/species_kmc_network.json',
            'Fixed initial motif population',
        )
    for t in temps:
        add(
            '45_state',
            'temperature',
            t,
            'K',
            'chosen_control',
            'docs/species_kmc_results/summary.json',
            'Temperature ensemble',
        )
    paths = list(csv.DictReader((ROOT / 'data/literature/sin_hf_pathways.csv').open()))
    for r in paths:
        for k in ('activation_energy_eV', 'physisorption_energy_eV', 'reaction_energy_eV'):
            add(
                'SiN_HF_source',
                r['pathway'] + '_' + k,
                float(r[k]),
                'eV',
                'published_DFT; Ea conditionally used, other energies reference-only',
                'data/literature/sin_hf_pathways.csv',
                r['surface_motif_reaction'],
            )
    rates = []
    for t in temps:
        for e, k in zip(network['events'], rate_constants(network, t, nu, arrival, des)):
            rates.append(
                dict(
                    event_id=e['id'],
                    source_state=network['states'][e['source']]['id'],
                    kind=e['kind'],
                    pathway=e['pathway'],
                    temperature_K=t,
                    barrier_eV=e['barrier_eV'],
                    prefactor_s=nu if e['kind'] == 'reaction' else None,
                    rate_s=float(k),
                    activation_free_energy_eV=None,
                    status=e['rate_status'],
                )
            )
    write_csv(out / 'species_event_rates.csv', rates)
    multi = json.loads((ROOT / 'docs/multilayer_results/summary.json').read_text())
    for k, v in multi['parameters'].items():
        add(
            'multilayer',
            k,
            v,
            'K'
            if k.endswith('_K')
            else 'angstrom'
            if k.endswith('_A')
            else 's^-1'
            if k in ('arrival_s', 'desorption_s')
            else 's',
            'assumed_control',
            'docs/multilayer_results/summary.json',
            'Arrival attenuates with local depth; zero during purge',
        )
    for k, v, u, meaning in [
        ('thermal_prefactor', 1e12, 's^-1', 'All demonstration thermal events'),
        ('Si_molecule_release_barrier', 0.53, 'eV', 'Assumed desorption, not literature P4'),
        ('NH3_release_barrier', 0.45, 'eV', 'Assumed detached NH3 desorption'),
        ('SiCl_exchange_barrier', 0.75, 'eV', 'Assumed exchange'),
        ('Si_N_bond_cutoff', 2.05, 'angstrom', 'Initial graph construction'),
        ('lateral_columns', 36, '1', '6 by 6 moving accessibility columns'),
        ('depth_bands', 6, '1', 'Unit-cell bands, not atomic monolayers'),
        ('fixed_bottom_bands', 1, '1', '56 fixed host atoms'),
    ]:
        add(
            'multilayer',
            k,
            v,
            u,
            'assumed_model',
            'plasma_surface/multilayer.py; scripts/build_multilayer_graph.py',
            meaning,
        )
    for k, v in [('H', 0.2), ('F', 0.5), ('Cl', 0.3)]:
        add(
            'multilayer',
            'initial_Si_cap_probability_' + k,
            v,
            '1',
            'assumed',
            'scripts/build_multilayer_graph.py',
            'N dangling valences initially capped with H',
        )
    rules = []
    pathmap = {r['pathway']: float(r['activation_energy_eV']) for r in paths}
    for stage in (1, 2, 3):
        for suffix, condition in [
            ('a', 'N with at least 2 H'),
            ('c', 'N with 1 H'),
            ('d', 'N with 0 H'),
        ]:
            key = f'P{stage}{suffix}'
            rules.append(
                dict(
                    event='SiN_cleavage',
                    condition=f'Si F count {stage - 1}; ' + condition,
                    pathway=key,
                    barrier_eV=pathmap[key],
                    prefactor_s=1e12,
                    status='unvalidated environment transfer',
                )
            )
        key = f'P{stage}b'
        rules.append(
            dict(
                event='SiH_exchange',
                condition=f'Stage {stage}; stage capped at 3',
                pathway=key,
                barrier_eV=pathmap[key],
                prefactor_s=1e12,
                status='unvalidated environment transfer',
            )
        )
    for ev, condition, key, barrier, status in [
        ('SiSi_cleavage', 'eligible Si-Si bond', 'P2f', 0.79, 'unvalidated environment transfer'),
        (
            'SiN_cleavage',
            'Si F >=3 and partner N H=1',
            'P4',
            0.53,
            'unvalidated environment transfer',
        ),
        (
            'SiN_cleavage',
            'Si F >=3 and partner N H=0 or >=2',
            None,
            None,
            'disabled: missing matched final barrier',
        ),
        ('SiCl_exchange', 'adsorbed HF and Si-Cl cap', None, 0.75, 'assumed'),
        (
            'Si_molecule_release',
            'zero backbonds; four ligands; no HF precursor',
            None,
            0.53,
            'assumed',
        ),
        ('NH3_release', 'zero backbonds; N with 3 H', None, 0.45, 'assumed'),
    ]:
        rules.append(
            dict(
                event=ev,
                condition=condition,
                pathway=key,
                barrier_eV=barrier,
                prefactor_s=1e12 if barrier is not None else None,
                status=status,
            )
        )
    write_csv(out / 'multilayer_rate_rules.csv', rules)
    for r in csv.DictReader((ROOT / 'data/literature/diffusion_barriers.csv').open()):
        add(
            'diffusion_reference',
            r['id'],
            float(r['barrier_ev']),
            'eV',
            'not used in active kMC; prefactor missing',
            r['doi'],
            r['system'] + '; ' + r['process'],
        )
    f2 = json.loads((ROOT / 'data/literature/si_f2_kinetics.json').read_text())
    for law in f2['laws']:
        for k in ('barrier_eV', 'A0_m3_s', 'n1', 'gamma0', 'n2', 'reference_temperature_K'):
            add(
                'F2_first_event_' + law['facet'],
                k,
                law[k],
                'eV'
                if k == 'barrier_eV'
                else 'm^3 s^-1'
                if k == 'A0_m3_s'
                else 'K'
                if k.endswith('_K')
                else '1',
                'published cluster DFT/TST fit; separate first-event model',
                'data/literature/si_f2_kinetics.json',
                'Fit range 298.15-1000 K; not multilayer ALE',
            )
    salts = json.loads((ROOT / 'data/literature/sin_salt_energetics.json').read_text())
    for r in salts['formation'] + salts['desorption']:
        for k, v in r.items():
            if k.endswith('_eV') or k.endswith('_C'):
                add(
                    'salt_reference',
                    r['id'] + '_' + k,
                    v,
                    'eV' if k.endswith('_eV') else 'C',
                    'reference-only; desorption energies are not TS barriers',
                    salts['source'],
                    r['reaction'],
                )
    for name, unit, meaning in [
        ('activation_DeltaG', 'eV', 'No matched full IS/TS Gibbs free energies supplied'),
        ('activation_DeltaF', 'eV', 'No full matched IS/TS harmonic mode sets supplied'),
        ('ZPE_correction', 'eV', 'Not evaluated for active kMC barriers'),
        ('activation_entropy', 'eV/K', 'Not available; prefactor 1e12 is an assumption'),
        (
            'transmission_coefficient',
            '1',
            'Unknown; helper default 1 is an explicit TST assumption',
        ),
        ('diffusion_hop_rate', 's^-1', 'No enabled hop events in active engines'),
        ('diffusion_coefficient', 'm^2/s', 'No matched hop network or diffusive trajectory fit'),
        ('AIMD_rate', 's^-1', 'No AIMD event-time ensemble generated'),
        (
            'gas_chemical_potentials',
            'eV',
            'No pressure/standard-state thermochemistry in graph model',
        ),
    ]:
        add(
            'missing',
            name,
            None,
            unit,
            'not computed or not supplied',
            'plasma_surface/thermal_kinetics.py',
            meaning,
        )
    write_csv(out / 'parameters.csv', records)
    source_rates = []
    for row in json.loads((ROOT / 'data/literature/jung2020_coadsorption_rates.json').read_text())[
        'rows'
    ]:
        for t in (300.0, 400.0, 450.0, 500.0):
            source_rates.append(
                dict(
                    **row,
                    temperature_K=t,
                    barrier_eV=row['Ea_over_R_K'] * KB_EV,
                    conditional_rate_s=arrhenius(row['Ea_over_R_K'] * KB_EV, t, row['A_s']),
                    enabled_in_current_kmc=False,
                )
            )
    write_csv(out / 'published_coadsorbate_rates.csv', source_rates)
    files = [
        'plasma_surface/model.py',
        'plasma_surface/multilayer.py',
        'plasma_surface/species_kmc.py',
        'plasma_surface/thermal_kinetics.py',
        'configs/materials.json',
        'configs/species_kmc_network.json',
        'data/literature/sin_hf_pathways.csv',
        'data/literature/si_f2_kinetics.json',
        'data/literature/diffusion_barriers.csv',
        'data/literature/jung2020_coadsorption_rates.json',
        'data/literature/sin_salt_energetics.json',
        'docs/multilayer_results/summary.json',
        'docs/species_kmc_results/summary.json',
    ]
    manifest = dict(
        parameter_rows=len(records),
        species_event_rate_rows=len(rates),
        multilayer_rule_rows=len(rules),
        published_rate_rows=len(source_rates),
        temperatures_K=temps,
        runner='scripts/audit_kinetic_parameters.py',
        runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        input_sha256={p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in files},
        scope='Current physical kinetics parameters and retained literature energetic inputs; optimizer settings and full model weights remain in per-calculation manifests',
        zero_event_100ps_one_site_example=md_count_rate(0, 1e-10),
    )
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    doc = ROOT / 'docs/KINETIC_PARAMETERS.md'
    table = '\n'.join(
        f"| {r['pathway']} | {float(r['activation_energy_eV']):.2f} | {float(r['physisorption_energy_eV']):.2f} | {float(r['reaction_energy_eV']):.2f} |"
        for r in paths
    )
    text = f'''# Kinetic parameters, free energies and evidence

**Current implementation:** published DFT barriers with assumed Arrhenius prefactors drive the conditional HF models. A separate F2 first-event model uses a published cluster-DFT/TST fit. No AIMD-derived rates, full activation free energies or diffusion hops are enabled in the active kMC engines. The strict multilayer qualified-rate library is empty.

The reproducible audit contains **{len(records)} parameter records**, **{len(rates)} event/temperature rates** (55 channels at four temperatures), **{len(rules)} multilayer thermal rules**, and **{len(source_rates)} source-specific coadsorbate rate evaluations**. Numerical values and missing entries are distinguished explicitly; a missing free energy is not zero.

[All parameter values, units, status and source](../data/kinetic_audit/parameters.csv) | [Every motif event rate](../data/kinetic_audit/species_event_rates.csv) | [Multilayer rules](../data/kinetic_audit/multilayer_rate_rules.csv) | [Input hashes](../data/kinetic_audit/manifest.json).

## Active controls and assumptions

| Quantity | Motif model | Multilayer demonstration | Evidence |
|---|---|---|---|
| Temperature | 325, 350, 400, 450 K | 450 K | Chosen simulation controls |
| HF arrival | 5 s^-1 per available motif | 5 exp(-depth / 1 A) s^-1 in dose | Assumed; not derived from measured flux |
| HF complex desorption | 100 s^-1 | 100 s^-1 | Assumed |
| Thermal prefactor | 1e12 s^-1 | 1e12 s^-1 | Assumed; no event-specific vibrational entropy |
| Dose / purge | 2 s / 1 s | 2 s / 1 s repeated | Chosen protocol |
| Access depth / attenuation | Not spatial | 2 A / 1 A | Assumed geometric rule |
| Detached Si-product release Ea | Literature motif event where specified | 0.53 eV | Multilayer desorption value is assumed, not source P4 |
| Detached NH3 release Ea | Included in source motif channel | 0.45 eV | Multilayer desorption assumed |
| Si-Cl exchange Ea | Not included | 0.75 eV | Assumed |
| Final SiF3-NH2 / SiF3-N cleavage | Disabled | Disabled | Missing matched final-step barriers |
| Lateral diffusion | Absent | Absent | No calibrated hop rates |
| Activation free energy | Not supplied | Not supplied | No matched full IS/TS thermochemistry |

The audit also lists all ten generic two-state model defaults, recipe fluxes/energies, composition-card overrides, initial motif fractions, graph cap probabilities, F2 law parameters and salt-reference energies. The generic model's 15/40 eV ion thresholds are not thermal activation barriers; its 0.65 eV modifier loss barrier is illustrative.

## Source SiN:H/HF energies

| Path | Ea (eV) | Ephy (eV) | Reaction DeltaE (eV) |
|---|---:|---:|---:|
{table}

[Source extraction and reference caveats](DRY_ETCH_PARAMETERS.md#sinhhf-all-fluorination-paths-and-competing-fragments). These are different source environments, not a common branching state. `Ea`, physisorption energy and reaction energy have separate meanings. In the current conditional kinetics, only `Ea` enters `nu exp(-Ea/kBT)`; Ephy and DeltaE do not silently become rate barriers.

## Thermal-rate workflow

The new [thermal kinetics module](../plasma_surface/thermal_kinetics.py) implements the following calculations with explicit units. It produces numerical candidates and does not automatically enable a kMC event.

1. **Arrhenius:** `k = nu exp(-Ea / kBT)`, with `Ea` in eV, `nu` in s^-1 and T in K. A 0.1 eV barrier change is exponentially significant; uncertainty should propagate through the rate rather than disappear into a fitted prefactor.
2. **Harmonic adsorbate thermochemistry:** `F_vib = sum[epsilon/2 + kBT ln(1-exp(-epsilon/kBT))]`, using positive mode energies `epsilon = h nu`. The first term is ZPE. The input IS and TS must represent the same environment and constrained degrees of freedom. Exactly one TS unstable mode is excluded; zero/imaginary stable-mode inputs are rejected. Soft translations or hindered rotations require an explicit physical treatment, not arbitrary deletion. [ASE thermochemistry implementation](https://docs.ase-lib.org/_modules/ase/thermochemistry.html).
3. **Surface HTST:** `DeltaF_dagger = (E_TS-E_IS) + F_vib_TS - F_vib_IS`, then `k = kappa (kBT/h) exp(-DeltaF_dagger/kBT)`. This is a fixed-volume surface Helmholtz approximation. It is not automatically a Gibbs activation free energy. A Gibbs treatment also needs consistent standard states, pressure/chemical potentials, and applicable gas translations/rotations. Do not multiply an independent Arrhenius prefactor onto the Eyring result.
4. **Gas arrival:** ideal-gas flux `J = p/sqrt(2 pi m kBT_g)` and per-site hazard `a_ads = J S A_site`. Inputs are pressure (Pa), gas temperature, molecular mass, site area and sticking. Surface temperature belongs in thermal surface hazards; it need not equal gas temperature.
5. **Surface diffusion:** `k_hop = nu_hop exp(-E_m/kBT)` and `D = f z k_hop l^2/(2d)` for equivalent isotropic neighbor hops, with per-neighbor rate, hop length l, coordination z, dimension d and correlation factor f. Current Cl/Si(111) reference barriers are 1.73 eV for Cl and 1.34 eV for SiCl; missing prefactors and differing environments prevent a numerical current-model diffusion rate. [Recorded sources](../data/literature/diffusion_barriers.csv).

The module requires actual matched mode data for HTST. The coadsorbate campaign's partial, unweighted adsorbate Hessian is not a complete mass-weighted IS/TS spectrum and is therefore not used as thermochemical input.

## AIMD and thermal kinetics: what is still needed

AIMD can test precursor stability, proton-transfer events and temperature-dependent motion. A short trajectory with no reaction does not measure a zero rate or an activation barrier. Define reactant/product basins, exclude recrossings and accumulate eligible-site residence time. Under a stationary Poisson model, `k_hat = N_committed_events / total_reactant_site_time`. The helper `md_count_rate` supplies the estimate and a one-sided Poisson upper bound; these assumptions require checking with replicas and waiting-time statistics.

For illustration only, zero events during 100 ps at one continuously eligible site gives a 95% upper bound of approximately **3.0e10 s^-1**, not useful evidence for slow etch kinetics. This is an analytical observation-window example, not a new AIMD result. Use longer/replicated sampling, validated rare-event methods or matched DFT/NEB+thermochemistry for activated events. Do not infer a barrier from one temperature without a justified prefactor.

Diffusion from AIMD/MD requires unwrapped coordinates and a diffusive fitting interval: `D = slope(MSD)/(2d)`, or the corresponding velocity-autocorrelation integral. Reject ballistic, caged or drifting regimes and estimate uncertainty across independent samples. [LAMMPS diffusion guidance](https://docs.lammps.org/Howto_diffusion.html).

## Calculated rates and next gate

[Source coadsorbate thermal rates](../data/kinetic_audit/published_coadsorbate_rates.csv) evaluate the existing Jung 2020 R6/R7 parameterization at 300/400/450/500 K. They describe its fluorinated clusters and vibrational-assistance assumptions, not our new ideal-slab candidates. [Primary paper](https://doi.org/10.1116/1.5125569).

The next critical graph rates are the **13 bare-N and eight NH2 final-cleavage environments** responsible for the [12-cycle stall](multilayer_results/cycle_diagnosis.json). For each, obtain matched relaxed endpoints, an index-one saddle connected to both basins, consistent energy/force accuracy and full vibrational treatment; then record geometry applicability, uncertainty and artifact hashes in the [evidence gate](../plasma_surface/rate_evidence.py). No energy from a rigid scan or an unconverged optimizer is promoted to a rate.

```sh
python scripts/audit_kinetic_parameters.py
python -m pytest tests/test_thermal_kinetics.py -q
```
'''
    doc.write_text(text, encoding='utf-8')
    print(
        json.dumps(
            {
                k: manifest[k]
                for k in (
                    'parameter_rows',
                    'species_event_rate_rows',
                    'multilayer_rule_rows',
                    'published_rate_rows',
                )
            }
        )
    )


if __name__ == '__main__':
    main()
