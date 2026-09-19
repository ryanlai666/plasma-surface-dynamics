"""Update only the owned reaction-network block in the atomistic README section."""
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[1]
def main():
    d=json.loads((ROOT/'docs/species_kmc_results/summary.json').read_text())
    block="""<!-- BEGIN SPECIES NETWORK -->
### Reaction networks and species-resolved kMC

The [HiPRGen pilot](docs/HIPRGEN.md) retained **40 forward and 40 reverse molecular substitution candidates** across Si-N/Si-O motifs with F/Cl. These pass composition and fixed-cap filters; they are not verified surface reactions. The [intermediate gap inventory](data/reaction_network/intermediate_gaps.csv) identifies missing precursor, proton-transfer, backbond, product-retention and ion-driven states.

The HF/SiN:H kMC rerun now tracks **45 named states and 55 enabled events**, representing all 16 published pathway entries. It includes HF complexes, successive fluorination, NH/NH2 and Si-H environments, Si-Si cleavage, and named NH3/H2/SiF4/SiH2F2/SiHF3 products. Missing release barriers leave branches disabled.

![Enabled species-resolved reaction network](docs/reaction_network/species_kmc_network.png)

![Rerun species kinetics and products](docs/species_kmc_results/species_kinetics.png)

**Scientific status:** this is conditional sensitivity, not calibrated ALE. The source Ea values, assumed prefactors, arrival/desorption hazards and proposed connections between source motifs are recorded per event. Unconverged atomistic peaks are excluded. No EPC or film-composition prediction is claimed.

Python and C++ implementations were checked against the independent master equation, with exact Si/N/H/F accounting. The temperature campaign contains **512 C++ trajectories with 1,000 initial motifs each**; independent Python runs check backend agreement.

<details>
<summary>HiPRGen candidate network, intermediate populations, and the new 2D/3D kMC animation</summary>

![HiPRGen bounded candidate networks](docs/reaction_network/hiprgen_candidates.png)

![Named intermediate populations](docs/species_kmc_results/intermediate_populations.png)

![Actual species kMC snapshots, 2D and perspective](docs/species_kmc_results/species_kmc.gif)

The animation is an actual 400-motif kMC trajectory on a schematic site grid; it is not an atomistic crystal or physical height measurement. No frames are interpolated.

</details>

[Species kMC report, equations, assumptions and verification](docs/species_kmc_results/REPORT.md) | [Every state/event](configs/species_kmc_network.json) | [HiPRGen scope and completeness audit](docs/HIPRGEN.md) | [Event rates and sources](docs/species_kmc_results/event_rates.csv).
<!-- END SPECIES NETWORK -->
"""
    p=ROOT/'README.md';raw=p.read_bytes();nl='\r\n' if b'\r\n' in raw else '\n';s=raw.decode('utf8').replace('\r\n','\n');before=s.split('<!-- BEGIN ATOMISTIC RESULTS -->')[0];after=s.split('<!-- END ATOMISTIC RESULTS -->')[1]
    if '<!-- BEGIN SPECIES NETWORK -->' in s:
        a,b=s.split('<!-- BEGIN SPECIES NETWORK -->',1);_,c=b.split('<!-- END SPECIES NETWORK -->',1);s=a+block+c
    else:s=s.replace('### Literature-parameterized dry-etch kinetics',block+'\n### Literature-parameterized dry-etch kinetics')
    assert s.split('<!-- BEGIN ATOMISTIC RESULTS -->')[0]==before and s.split('<!-- END ATOMISTIC RESULTS -->')[1]==after
    p.write_bytes(s.replace('\n',nl).encode('utf8'))
if __name__=='__main__':main()
