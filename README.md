# Longitudinal-COMBAT-GAM-and-analysis-of-longitudinal-data
In this project, data is harmonized from participants who have been in ADNI for at least 10 years.  Following harmonization, we looked into the progression of cognitive decline by APOE groups by looking at hippocampal volume and meta-ROI


# Parallel Developtment for Speed Up

Professfor Tosun has indicated that 1000 iterations is prefered. The original code running in serial, requires 2 hours per bio marker
to compute 100 iteration boot strap.

## Parallel Library

First attempt to speed up will be the native R package `parallel`. A data set has been included in this repo for testing,
`lcmm_modelling_data`.

There is also a testing entry point script, `running_functions.R`.

