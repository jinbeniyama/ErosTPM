# Eros TPM
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[developer mail](mailto:jbeniyama@oca.eu)

## Overview
This is a repository for Eros TPM in prep.
Figures are made in `./plot`.

## Data (in /data)
* `Eros_UKIRT_June_1998_modified_by_JB.txt` (provided by A. Harris and modified by J. Beniyama)
* `erosLim-2002sept22_Jy.flx` (thermal flux from Lim+2005, scanned from their Fig.7 on p. 398)
* `SPITZER_S0_4872960_0001_9_E7275827_tune.tbl`, `SPITZER_S0_4872960_0004_9_E7275831_tune.tbl`,
  `SPITZER_S0_4872960_0007_9_E7275841_tune.tbl`, `SPITZER_S0_4872960_0010_9_E7275844_tune.tbl`,
  `SPITZER_S2_4872960_0012_9_E7275703_tune.tbl`, `SPITZER_S2_4872960_0013_9_E7275707_tune.tbl`,
  `SPITZER_S2_4872960_0014_9_E7275704_tune.tbl`, `SPITZER_S2_4872960_0015_9_E7275696_tune.tbl`
  (SST/IRS spectra, downloaded from `https://pds-smallbodies.astro.umd.edu/data_other/Spitzer.shtm`)
* `pre_akari.dat` (AKARI observations downloaded from `https://darts.isas.jaxa.jp/astro/akari/data/AKARI-IRC_Catalogue_AllSky_ASTFLUX_1.0.html`)
* `433.obj` (shape model of Eros, downloaded from DAMIT, `https://astro.troja.mff.cuni.cz/projects/damit/asteroid_models/view/3083`, `shape.obj`)
* `433_spin.txt` (spin file of Eros, downloaded from DAMIT, `https://astro.troja.mff.cuni.cz/projects/damit/asteroid_models/view/3083`, `spin.txt`)
* `433_obs_N998.txt` (thermal observations formatted for TPM, read the paper for the details, made by J.B.)
* `433_eph_N998.txt` (ephemeris file paired with `433_obs_N998.txt`, made by J.B.)
* `Eros_flux_N998.txt` (observation files after preprocessing, which is made with `Eros_prepro_flux.py`)

## Preprocesses
Execute folloing commands in `./data` to make obs file and ephemeris file for the TPM.
``` 
# Preprocess of thermal infrared fluxes.
python ../script/Eros_prepro_flux.py ../data/Eros_UKIRT_June_1998_modified_by_JB.txt ../data/erosLim-2002sept22_Jy.flx ../data/pre_akari.dat
``` 

```
# Make obs and ephemeris files
make_obseph.py 433 Eros_flux_N998.txt --out_obs 433_obs_N998.txt --out_eph 433_eph_N998.txt
```

## Procedure
Skip and go to "Plot figures" if you just make figures in the paper.

1. Do TPM with brute-force method
```
bash 1_ErosTPM.sh (obj file) (spin file) (obs file) (ephemeris file) (output directory)
bash 1_ErosTPM.sh 433.obj 433_spin.txt 433_obs_N811.txt 433_eph_N811.txt TPMresult
```

2. Make neural-network (NN) model to predict thermal flux
```
bash 2_makeNNmodel.sh (directory with output files of 1.)
bash 2_makeNNmodel.sh TPMresult
```

3. Predict thermal fluxes with the NN model
```
bash 3_predictflux.sh (NN model made in process 2.)
bash 3_predictflux.sh NNmodel
```
The output files are in the format `.npy` such as `LUT_2450991.767627034.npy`.


## Plot figures
Do all commands in `/plot`.

``` 
# Location of Eros
python ../script/Eros_fig_loc.py
```

``` 
# Thermal flux of Eros
python ../script/Eros_fig_flux.py ../data/Eros_flux_N811.txt
```
