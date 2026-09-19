# Gas chemical potentials: source-bounded thermal terms

The [NIST coefficient library](../data/literature/nist_gas_shomate.json) supplies Shomate thermochemistry for HF, SiF4, NH3, H2, HCl and water vapor. These gases cover reactants, volatile products and possible competing channels; their presence in this table does not establish an elementary surface reaction. Every entry records the source URL and temperature range. All coefficients use the source SI tables.

[81 evaluated gas states](../data/kinetic_audit/gas_chemical_potentials.csv) cover selected temperatures and partial pressures of 1 Pa, 1 Torr and 1 bar. **Water below 500 K is excluded**, because it falls outside the retrieved fit. The existing 450 K kMC therefore does not acquire a water free-energy correction from this dataset.

With t=T/1000, the source fits give:

```text
Cp = A + Bt + Ct^2 + Dt^3 + E/t^2                 [J mol^-1 K^-1]
H(T)-H(298.15) = At+Bt^2/2+Ct^3/3+Dt^4/4-E/t+F-H [kJ mol^-1]
S? = A ln(t)+Bt+Ct^2/2+Dt^3/3-E/(2t^2)+G         [J mol^-1 K^-1]
mu(T,p)-H(298.15) = DeltaH - T S? + RT ln(p/p?)
p? = 1 bar
```

The code converts all terms to consistent units before summation. The exported `mu_element_enthalpy_reference_eV` adds the source formation enthalpy at 298.15 K. This is an elemental **enthalpy** reference; it is neither an absolute DFT energy nor the standard Gibbs energy of formation. Element reference offsets cancel only in complete, balanced reactions with compatible references for every participant. Do not add these numbers directly to an unaligned DFT total energy.

Surface activation free energies still require matched IS/TS electronic energies, vibrational modes, constraints and standard states. Adsorption also requires a consistent gas-to-surface standard-state conversion and site convention. NIST gas data reduce a thermochemistry input gap; they do not validate the multilayer cleavage barriers, sticking probabilities or diffusion rates.

Reproduce with `python scripts/report_gas_thermochemistry.py`. Tests verify source table values, dH/dT=Cp, the ideal-gas pressure dependence and fit-range enforcement. [Provenance manifest](../data/kinetic_audit/gas_thermochemistry_manifest.json).
