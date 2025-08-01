# Eros TPM
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[developer mail](mailto:jbeniyama@oca.eu)

## Overview
This is a repository for Eros TPM in prep.
Figures are made in `./plot`.

## Data (in /data)
* `433.obj` (shape model of Eros, downloaded from DAMIT, `https://astro.troja.mff.cuni.cz/projects/damit/asteroid_models/view/3083`, `shape.obj`)
* `433_eph_N998.txt` (ephemeris file paired with `433_obs_N998.txt`, made by J.B.)
* `433_obs_N998.txt` (thermal observations formatted for TPM, read the paper for the details, made by J.B.)
* `433_spin.txt` (spin file of Eros, downloaded from DAMIT, `https://astro.troja.mff.cuni.cz/projects/damit/asteroid_models/view/3083`, `spin.txt`)
* `Eros_UKIRT_June_1998_modified_by_JB.txt` (provided by A. Harris and modified by J. Beniyama)
* `Eros_flux_N998.txt` (observation files after preprocessing, which is made with `Eros_prepro_flux.py`)
* `SPITZER_S0_4872960_0001_9_E7275827_tune.tbl`, `SPITZER_S0_4872960_0004_9_E7275831_tune.tbl`,
  `SPITZER_S0_4872960_0007_9_E7275841_tune.tbl`, `SPITZER_S0_4872960_0010_9_E7275844_tune.tbl`,
  `SPITZER_S2_4872960_0012_9_E7275703_tune.tbl`, `SPITZER_S2_4872960_0013_9_E7275707_tune.tbl`,
  `SPITZER_S2_4872960_0014_9_E7275704_tune.tbl`, `SPITZER_S2_4872960_0015_9_E7275696_tune.tbl`
  (SST/IRS spectra, downloaded from `https://pds-smallbodies.astro.umd.edu/data_other/Spitzer.shtm`)
* SSTinfo.txt (SST pointing information)
* `erosLim-2002sept22_Jy.flx` (thermal flux from Lim+2005, scanned from their Fig.7 on p. 398)
* `pre_akari.dat` (AKARI observations downloaded from `https://darts.isas.jaxa.jp/astro/akari/data/AKARI-IRC_Catalogue_AllSky_ASTFLUX_1.0.html`)

## Preprocesses
Execute folloing commands in `./data` to make obs file and ephemeris file for the TPM.

- Preprocess of thermal infrared fluxes.
``` 
# With 8 SST spectra (N=998)
#python ../script/Eros_prepro_flux.py ../data/Eros_UKIRT_June_1998_modified_by_JB.txt ../data/erosLim-2002sept22_Jy.flx ../data/pre_akari.dat
# With only 4 SST spectra (N=602)
#python script/Eros_prepro_flux.py data/Eros_UKIRT_June_1998_modified_by_JB.txt data/erosLim-2002sept22_Jy.flx data/pre_akari.dat
# With only 2 averaged SST spectra (N=426)
python script/Eros_prepro_flux.py data/Eros_UKIRT_June_1998_modified_by_JB.txt data/erosLim-2002sept22_Jy.flx data/pre_akari.dat
``` 

- Make obs and ephemeris files (in `/data`)
```
# With only 8 SST spectra (N=998)
#make_obseph.py 433 Eros_flux_N998.txt --out_obs 433_obs_N998.txt --out_eph 433_eph_N998.txt
# With only 4 SST spectra (N=602)
#make_obseph.py 433 Eros_flux_N602.txt --out_obs 433_obs_N602.txt --out_eph 433_eph_N602.txt
# With only 2 averaged SST spectra (N=426)
## Nominal error
make_obseph.py 433 Eros_flux_N426.txt --out_obs 433_obs_N426_nominal.txt --out_eph 433_eph_N426.txt
## 5% error for Spitzer
make_obseph.py 433 Eros_flux_N426.txt --out_obs 433_obs_N426.txt --out_eph 433_eph_N426.txt
```

## Procedure
Skip and go to "Plot figures" if you just make figures in the paper.

1. Do TPM with brute-force method
```
# TPM N=998
#bash 1_ErosTPM.sh ../433.obj ../433_spin.txt 433_obs_N998.txt 433_eph_N998.txt TPMresult_20250508_LagerrosApp

# TPM with new grids on 2025-07-24
bash script/ErosTPM_newgrid.sh data/433.obj data/433_spin.txt data/433_obs_N602.txt data/433_eph_N602.txt TPMres_normal_20250724_re
# Test with new grids with landscape model on 2025-07-24
bash script/ErosTPM_landscape_newgrid.sh data/433.obj data/433_spin.txt data/433_obs_N602.txt data/433_eph_N602.txt TPMres_landscape_20250724

# TPM with new grids on 2025-08-01
bash script/ErosTPM_newgrid.sh data/433.obj data/433_spin.txt data/433_obs_N426.txt data/433_eph_N426.txt TPMres_normal_20250801
bash script/ErosTPM_newgrid.sh data/433.obj data/433_spin.txt data/433_obs_N426.txt data/433_eph_N426.txt TPMres_normal_20250801_SST_Lim_error_update

```

```
# Make look up table.
make_lut.py TPMresult_20250508_LagerrosApp --out data/lut_20250508.txt
# Test
make_lut.py TPMres_normal_20250724 --out data/lut_normal_20250725.txt
```

2. **Make neural-network (NN) model to predict thermal flux**
```
DemoNN.ipynb
```
The results are saved in `saved_model`.

3. **Predict thermal fluxes with the NN model**
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
#python ../script/Eros_fig_flux.py ../data/Eros_flux_N998.txt
python ../script/Eros_fig_flux.py ../data/Eros_flux_N426.txt
```

- Single-component fit, TI vs. chi-squared (Figure X., in prep.)
```
plot_tpm_brute.py ../TPMresult_2025Jan/tpmout* -x TI --reduce --N_param 3 --logx --logy --ylim 20 200 --out chi2_single_2025Jan.pdf
plot_tpm_brute.py ../TPMresult/tpmout* -x TI --reduce --N_param 3 --logx --logy --ylim 20 200 --out chi2_single_wo_Laggerros_app.pdf
plot_tpm_brute.py ../TPMresult_20250508_LagerrosApp/tpmout_433_brute_ti*  -x TI --reduce --N_param 2 --fixscale --logx --logy --ylim 20 200 --out chi2_single_w_Laggerros_app.pdf
# New grid 
plot_tpm_brute.py ../TPMres_normal_20250724/tpmout_433_brute_ti* --fixscale -x TI --reduce --N_param 2 --logx --logy --ylim 10 1500 --out chi2_normal_20250725.pdf
plot_tpm_brute.py ../TPMres_landscape_20250724/lstpmout_433_brute_ti* --fixscale  -x TI --reduce --N_param 2 --logx --logy --ylim 10 1500 --out chi2_landscape_20250725.pdf

plot_tpm_brute.py ../TPMres_normal_20250730/tpmout_433_brute_ti* --fixscale -x TI --reduce --N_param 2 --logx --logy --ylim 10 1500 --out chi2_normal_20250730.pdf
plot_tpm_brute.py ../TPMres_normal_20250730/tpmout_433_brute_ti* --scale_per_obs -x TI --reduce --N_param 2 --logx --logy --ylim 10 1500 --out chi2_normal_scale_per_obs_20250730.pdf

# SST 5% error
plot_tpm_brute.py ../TPMres_normal_20250801/tpmout_433_brute_ti* --fixscale -x TI --reduce --N_param 2 --logx --logy --ylim 10 1500 --out chi2_normal_20250801.pdf
plot_tpm_brute.py ../TPMres_normal_20250801/tpmout_433_brute_ti* --scale_per_obs -x TI --reduce --N_param 2 --logx --logy --ylim 10 1500 --out chi2_normal_scale_per_obs_20250801.pdf

# SST 5%, Lim 10%
plot_tpm_brute.py ../TPMres_normal_20250801_SST_Lim_error_update/tpmout_433_brute_ti* --fixscale -x TI --reduce --N_param 2 --logx --logy --ylim 4.5 400 --out chi2_SST_Lim_20250801.pdf
plot_tpm_brute.py ../TPMres_normal_20250801_SST_Lim_error_update/tpmout_433_brute_ti* --scale_per_obs -x TI --reduce --N_param 2 --logx --logy --ylim 1 10 --out chi2_SST_Lim_scale_per_obs_20250801.pdf
```

- Dual-component fit, TI vs. chi-squared (Figure X., in prep.)
```
plot_tpm_brute_NN.py NN_TPMres.txt --fixscale --reduce --logx --logy --N_param 2 --out Eros_fig_chi2single_NN.pdf
```

- Field of view at the time of SST observations (Figure X., in prep.)
```
python ../script/Eros_fig_SSTFoV.py ../data/433_akari_ukirt1998_ukirt2002_lim2005_3_SST_six_20250110.dat
``` 

## Dependencies
This repository is depending on `Python`, `NumPy`, `pandas`, `SciPy`, `Astropy`, `Astroquery`, `TREM`.
The tpm code is downloaded from `https://www.oca.eu/images/LAGRANGE/pages_perso/delbo/thermops.tar.gz`.
