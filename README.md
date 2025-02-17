# Eros TPM
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[developer mail](mailto:jbeniyama@oca.eu)

## Overview
This is a repository for Eros TPM in prep.
Figures are made in /plot.

## data (in /data)
* 433.obj (shape model of Eros, downloaded from DAMIT, https://astro.troja.mff.cuni.cz/projects/damit/asteroid_models/view/3083, shape.obj)
* 433_spin.txt (spin file of Eros, downloaded from DAMIT, https://astro.troja.mff.cuni.cz/projects/damit/asteroid_models/view/3083, spin.txt)
* 433_obs_20250110.txt (thermal observations formatted for TPM, read the paper for the details)
* 433_eph_20250110.txt (ephemeris file paired with 433_obs.txt)

## Procedure
Skip and go to "Plot figures" if you just make figures in the paper.

1. Do TPM with brute-force method
```
bash 1_ErosTPM.sh 433.obj 433_spin.txt 433_obs_20250110.txt 433_eph_20250110.txt .
bash 1_ErosTPM.sh 433.obj 433_spin.txt 433_obs_20250110.txt 433_eph_20250110.txt .
```

2. Make neural-network model to predict thermal flux

3. Predict thermal fluxes with single-component TPM

4. Predict thermal fluxes with single-component TPM

## Plot figures
Do all commands in /plot.
```
```
