# Paleozoic Seawater Lithium Models — BLAG F_sil Perturbation

Two-phase BLAG model (Berner, Lasaga & Garrels, 1983) with prescribed silicate weathering from lithium isotope data.

## Overview

This code implements a two-phase approach:

1. **Phase 1 (570–424 Ma):** Full BLAG feedback spinup — the system finds its own self-consistent equilibrium with CO₂-weathering feedback.
2. **Phase 2 (424–300 Ma):** Silicate weathering is prescribed from Li isotope data. The weathering intensity (WI) multiplier perturbs F_sil relative to the spinup end-state, while CO₂ and other variables respond via the full BLAG equations (Eqs. 59A–59H).

The silicate weathering forcing term is derived from δ⁷Li data converted to weathering/denudation (W/D). Since D (denudation) is held constant, changes in W can be directly related to changes in total silicate weathering.

## Requirements

- **MATLAB** (R2019b or later recommended)
  - No additional toolboxes required

## Files

| File | Description |
|------|-------------|
| `BLAG_model_Fsil_perturbation.m` | Main driver script — run this file |
| `BLAG_odes_Fsil.m` | ODE system function (Eqs. 59A–59H) |
| `data/GEOCARB_input_arrays_tMod.csv` | GEOCARB forcing data (land area fA, spreading rate fSR) |
| `data/Li_model_output_Ghosh2026.csv` | Li box model output (weathering intensity and pCO₂ estimates) |
| `BLAG_Fsil_perturbation_results.csv` | Example output results (424–300 Ma) |

## How to Run

1. Clone this repository:
   ```
   git clone https://github.com/ghosh17/PaleozoicSeawaterLithiumModels.git
   ```
2. Open MATLAB and navigate to the repository directory.
3. Run the main script:
   ```matlab
   BLAG_model_Fsil_perturbation
   ```

The script will:
- Load GEOCARB forcing and Li isotope data from the `data/` directory
- Run the Phase 1 spinup (570–424 Ma)
- Run the Phase 2 Li-forced simulation (424–300 Ma)
- Generate diagnostic plots
- Save results to `BLAG_Fsil_perturbation_results.csv`

## Input Data

### GEOCARB Forcing (`data/GEOCARB_input_arrays_tMod.csv`)
Contains land area factor (fA) and seafloor spreading rate factor (fSR) as functions of geological age, based on GEOCARB III (Berner & Kothavala, 2001).

### Li Model Output (`data/Li_model_output_Ghosh2026.csv`)
Contains weathering intensity (WI_high) and pCO₂ estimates derived from the Li box model (Ghosh, 2026). The WI values are normalized to the 424 Ma baseline to produce the cumulative weathering multiplier used in Phase 2.

## References

- Berner, R.A., Lasaga, A.C. & Garrels, R.M. (1983). The carbonate-silicate geochemical cycle and its effect on atmospheric carbon dioxide over the past 100 million years. *Am. J. Sci.* 283, 641–683.
- Berner, R.A. & Kothavala, Z. (2001). GEOCARB III: A revised model of atmospheric CO₂ over Phanerozoic time. *Am. J. Sci.* 301, 182–204.
- Ghosh (2026). Li box model output for Paleozoic seawater lithium isotope modeling.
