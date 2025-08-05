# Eros TPM
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[developer mail](mailto:jbeniyama@oca.eu)

## Overview
This is a repository for Eros TPM in prep.
Figures are made in `./plot`.

## Structure 
```
./
  data/
  script/
  .gitignored
  README.md
```

## Data (in /data)
* `433.obj` (shape model of Eros, downloaded from DAMIT, `https://astro.troja.mff.cuni.cz/projects/damit/asteroid_models/view/3083`, `shape.obj`)
* `433_eph_N448.txt` (ephemeris file paired with `433_obs_N448.txt`, made by J.B.)
* `433_obs_N448.txt` (thermal observations formatted for TPM, read the paper for the details, made by J.B.)
* `433_spin.txt` (spin file of Eros, downloaded from DAMIT, `https://astro.troja.mff.cuni.cz/projects/damit/asteroid_models/view/3083`, `spin.txt`)
* `Eros_UKIRT_June_1998_modified_by_JB.txt` (provided by A. Harris and modified by J. Beniyama)
* `Eros_flux_N448.txt` (observation files after preprocessing, which is made with `Eros_prepro_flux.py`)
* `SPITZER_S0_4872960_0001_9_E7275827_tune.tbl (not used)`, `SPITZER_S0_4872960_0004_9_E7275831_tune.tbl (not used)`,
  `SPITZER_S0_4872960_0007_9_E7275841_tune.tbl (not used)`, `SPITZER_S0_4872960_0010_9_E7275844_tune.tbl (not used)`,
  `SPITZER_S2_4872960_0012_9_E7275703_tune.tbl`, `SPITZER_S2_4872960_0013_9_E7275707_tune.tbl`,
  `SPITZER_S2_4872960_0014_9_E7275704_tune.tbl`, `SPITZER_S2_4872960_0015_9_E7275696_tune.tbl`
  (SST/IRS spectra, downloaded from `https://pds-smallbodies.astro.umd.edu/data_other/Spitzer.shtm`)
* `SSTinfo.txt` (SST pointing information)
* `erosLim-2002sept22_Jy.flx` (thermal flux from Lim+2005, scanned from their Fig.7 on p. 398)
* `lut_20250508.txt` (look-up table used to make NN model)
* `pre_akari.dat` (AKARI observations downloaded from `https://darts.isas.jaxa.jp/astro/akari/data/AKARI-IRC_Catalogue_AllSky_ASTFLUX_1.0.html`)

## Preprocesses
Execute folloing commands in `./data` to make obs file and ephemeris file for the TPM.

- Preprocess of thermal infrared fluxes.
``` 
# With only 2 averaged SST spectra (N=448)
python script/Eros_prepro_flux.py data/Eros_UKIRT_June_1998_modified_by_JB.txt data/erosLim-2002sept22_Jy.flx data/pre_akari.dat
``` 

- Make obs and ephemeris files (in `/data`)
```
# With only 2 averaged SST spectra (N=448)
make_obseph.py 433 Eros_flux_N448.txt --out_obs 433_obs_N448.txt --out_eph 433_eph_N448.txt
```

## Procedure
Skip and go to "Plot figures" if you just make figures in the paper.

1. **Do TPM with brute-force method**
```
# TPM with new grids on 2025-08-03
bash script/ErosTPM_newgrid.sh data/433.obj data/433_spin.txt data/433_obs_N448.txt data/433_eph_N448.txt TPMres_20250803
```

2. **Make look up table**
```
make_lut.py TPMresult_20250508_LagerrosApp --out data/lut_20250508.txt
```

3. **Make neural-network (NN) model to predict thermal flux**
```
DemoNN.ipynb
```
The results are saved in `saved_model`.

4. **Predict thermal fluxes with the NN model**
```
bash 3_predictflux.sh (NN model made in process 2.)
bash 3_predictflux.sh NNmodel
```
The output files are in the format `.npy` such as `LUT_2450991.767627034.npy`.


## Plot figures
Do all commands in `/plot`.

- Location of Eros (Figure X.)
``` 
python ../script/Eros_fig_loc.py
```

- Thermal flux of Eros (Figure X.)
``` 
python ../script/Eros_fig_flux.py ../data/Eros_flux_N448.txt
```

- Single-component fit wo/NN
```
# TI vs. chi-squared 
# SST 5%, Lim+ 10%, Wolters+2008 7%
plot_tpm_brute.py ../TPMres_20250803/tpmout_433_brute_ti* --fixscale -x TI --reduce --N_param 2 --logx --logy --ylim 1 100 --paper P14 --out chi2_single_20250803.pdf
#plot_tpm_brute.py ../TPMres_20250803/tpmout_433_brute_ti* --scale_per_obs -x TI --reduce --N_param 2 --logx --logy --ylim 1 4 --out chi2_single_scale_per_obs_20250803.pdf
#plot_tpm_brute.py ../TPMres_20250803/tpmout_433_brute_ti* --scale_per_obs -x TI --reduce --N_param 2 --logx --logy --ylim 1 4 --out chi2_single_scale_per_obs_20250803.pdf

# obs. vs model 
python ../script/Eros_fig_obsvsmodel.py /Users/beniyama/research/ErosTPM/paper/TPMres_20250803/tpmout_433_brute_ti60_ca90_cr0.4.dat --out single_obs_vs_model_20250801.pdf
```

- Single-component fit w/NN, TI vs. chi-squared 
```
(In prep.)
plot_tpm_brute_NN.py NN_TPMres.txt --fixscale --reduce --logx --logy --N_param 2 --out Eros_fig_chi2single_NN.pdf
```

- Dual-component fit wo/NN, TI vs. chi-squared (Figure X., in prep.)
```
blend_tpm_result_iter.py ../TPMres_20250803/tpmout* --resdir "../TPMres_20250803/" --fixscale --TI0 60 --out blended_flux_fixscale_iter.txt --T_typical 250
# -> TI_thresh = 34.80
plot_blended_result_RMS.py blended_flux_fixscale_iter.txt --reduce --dof 446 --TI_thresh 34.80
```
- Dual-component fit w/NN, TI vs. chi-squared (Figure X., in prep.)
```
(In prep.)
```

- Field of view at the time of SST observations (Figure X., in prep.)
```
# python ../script/Eros_fig_SSTFoV.py ../data/433_akari_ukirt1998_ukirt2002_lim2005_3_SST_six_20250110.dat
``` 

## Dependencies
This repository is depending on `Python`, `NumPy`, `pandas`, `SciPy`, `Astropy`, `Astroquery`, `TREM`.
The TPM code is downloaded from `https://www.oca.eu/images/LAGRANGE/pages_perso/delbo/thermops.tar.gz`.
