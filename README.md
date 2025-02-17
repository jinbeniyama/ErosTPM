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
* 433_akari_ukirt1998_ukirt2002_lim2005_3_SST_six_20250110.dat (observation files after preprocessing)

## Procedure
Skip and go to "Plot figures" if you just make figures in the paper.

1. Do TPM with brute-force method
```
bash 1_ErosTPM.sh (obj file) (spin file) (obs file) (ephemeris file) (output directory)
bash 1_ErosTPM.sh 433.obj 433_spin.txt 433_obs_20250110.txt 433_eph_20250110.txt TPMresult
```

2. Make neural-network (NN) model to predict thermal flux
```
bash 2_makeNNmodel.sh (directory with output files of 1.)
bash 2_makeNNmodel.sh TPMresult
```

3. Predict thermal fluxes with the NN model
```
bash 3_predictflux.sh (directory with output files of 2.)
bash 3_predictflux.sh 
```

## Plot figures
Do all commands in /plot.

``` 
# Location of Eros
python ../script/Eros_fig_loc.py
```


``` 
# Thermal flux of Eros
python ../script/Eros_fig_flux.py ../data/433_akari_ukirt1998_ukirt2002_lim2005_3_SST_six_20250110.dat
```
