# Data Directory

This directory should contain the following data files required by the BLAG model:

1. **GEOCARB_input_arrays_tMod.csv** - GEOCARB forcing data (fA, fSR) used in the BLAG model spinup and forcing phases
2. **Li_model_output_Ghosh2026.xlsx** - Li box model output containing WI (weathering intensity) data from Ghosh (2026)

These files are referenced by `BLAG_model_Fsil_perturbation.m` and must be placed in this directory before running the model.
