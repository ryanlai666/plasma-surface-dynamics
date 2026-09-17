# Equations and interpretation

Each independent column has one exposed site with state bare B or modified M. Coverage θ is the fraction of modified columns. A removal or deposition event exposes a new bare site. Columns have integer height offsets but no lateral interactions.

| Event | Hazard per eligible site (s⁻¹) | State change | Height change |
|---|---|---|---|
| Adsorption/modification | a = sticking × radical flux / site density | B → M | 0 |
| Desorption | d = prefactor × exp(−barrier / kBT) | M → B | 0 |
| Ion-assisted chemical removal | c = ion flux × Ychem / site density | M → B | −1 |
| Physical removal | s = ion flux × Yphys / site density | B or M → B | −1 |
| Generic deposition | g = deposition sticking × precursor flux / site density | B or M → B | +1 |

Y(E) = scale × max(sqrt(E/Ethreshold) − 1, 0). This is an illustrative threshold law, not a fit to a cited sputter-yield database. Monoenergetic normal-incidence exposure is assumed. Rate hazards represent expected removal counts; the event model does not simulate individual incident ion trajectories or multi-atom collision cascades.

The exact mean-field equations during each constant-rate phase are:

```text
dθ/dt = a − (a+d+c+s+g) θ
dR/dt = layer_nm × (c θ + s)
dG/dt = layer_nm × g
net removal = R − G
```

For k=a+d+c+s+g > 0, θ(t)=a/k+(θ₀−a/k)exp(−kt). The integrated coverage is a t/k+(θ₀−a/k)(1−exp(−kt))/k. Both backends evaluate the relaxation term with `expm1` for small-step stability. Coverage persists across phases and cycles. Cycle-5 EPC is the difference between cumulative removal after cycles 5 and 4; it is not automatically a steady-state value.

The stochastic solver draws an exponential waiting time from total hazard, chooses a reaction proportional to its propensity, then uniformly chooses an eligible column. At a phase boundary, the previous waiting time is discarded and new hazards are used, consistent with the memoryless exponential process. Rates do not depend on column height. The reported roughness is the population standard deviation of independent column heights, not a prediction of AFM roughness or trench morphology.

Illustrative defaults use 7×10¹⁸ sites/m², 0.136 nm/event, chemical threshold 15 eV, and physical threshold 40 eV. These are software demonstration choices. In particular, the default threshold separation must not be presented as reproducing the much narrower window discussed in the public Si–Cl₂–Ar study listed in DATASETS.md.

## Required checks before interpreting material trends

- Calibrate species-specific thresholds, rates, densities, and thickness conversion for the selected chemistry and film.
- Measure or model the ion energy/angular distributions, neutral flux, temperature, and pulse transients. Reactor power is not ion energy.
- Compare chemical-only, ion-only, and combined cycles; inspect dose and removal-time saturation, residual coverage, and over-etch.
- Propagate both fitted-rate uncertainty and model discrepancy. ML tree disagreement currently measures neither reliably.
- Validate on independent runs and compositions. An accurate surrogate of this model does not establish physical accuracy.

No self-consistent plasma discharge, electron kinetics, reactor transport, evolving Si/N ratio, charging, redeposition, lattice orientation, or chemical product accounting is implemented. These exclusions define the scientific scope of the current kernel.
