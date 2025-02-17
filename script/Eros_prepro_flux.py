#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Preprocesses of thermal fluxes of Eros.

1. Eros in 1998 published in Harris+1999
2. Eros in 2002 sep. published in Wolters+2008
3. Eros in 2002 sep. published in Lim+2005
4. Eros in 2004 by SST
5. Eros in 2007 by AKARI
"""
import os 
from argparse import ArgumentParser as ap
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astroquery.jplhorizons import Horizons
from astropy.constants import c, au
from astropy.time import Time

from myplot import mycolor


def make_df_W08():
    """
    """

def remove_largevar(df, key, var_th):
    # Remove large variation
    df = df.reset_index(drop=True)
    idx_list = []

    for idx, row in df.iterrows():

        if idx == 0:
            val0 = row[key]
            valmin = val0*(1 - var_th)
            valmax = val0*(1 + var_th)
            idx_list.append(idx)
            print(valmin)
            print(valmax)
        else:
            val = row[key]

            if val > valmax:
                pass
            elif val < valmin:
                pass
            else:
                # Save index
                idx_list.append(idx)
                # Update
                val0 = val
                valmin = val0*(1 - var_th)
                valmax = val0*(1 + var_th)

    print(idx_list)
    df = df.loc[idx_list]
    return df


def remove_edge(df, n):
    # Remove large variation
    df = df.reset_index(drop=True)
    N = len(df)
    idx_list = df.index.values.tolist()
    idx_use = idx_list[n:N-n]

    df = df.loc[idx_use]
    return df


def SST_ltcor(target, utc):
    epoch = Time(str(utc), format='isot', scale='utc')
    epoch = epoch.jd

    ## Lighttime correction
    ast = Horizons(location='@sst',id=target, epochs=epoch)
    ast = ast.ephemerides()
    au_m = au.to("m").value
    lt_s = ast["delta"]*au_m/c
    lt_day = lt_s / 3600./24.
    epoch_ltcor = epoch - lt_day
    print(f"  epoch_jd       : {epoch}")
    print(f"  epoch_jd_ltcor : {epoch_ltcor.value[0]}")
    return epoch_ltcor.value[0]


if __name__ == "__main__":
    parser = ap(
        description="Preprocess of thermal fluxes of Eros.")
    parser.add_argument(
        "f_H99", type=str,
        help="Flux tables from Harris+1999.")
    parser.add_argument(
        "f_L05", type=str,
        help="Flux tables from Lim+2005.")
    parser.add_argument(
        "f_AKARI", type=str,
        help="Flux tables from JAXA website.")
    parser.add_argument(
        "--fdir", type=str, default="../data",
        help="input file directory")
    parser.add_argument(
        "--outdir", type=str, default="fig",
        help="output directory")
    parser.add_argument(
        "--outtype", default="pdf",
        help="format of output figure")
    args = parser.parse_args()


    outdir = args.outdir
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    # Read input files

    # 1. Harris+1999 ==========================================================
    # N = 
    f_H99 = args.f_H99
    #print(f" Original N={len(df_H99)} (H99)")
    # 1. Harris+1999 ==========================================================

    # 2. Lim+2005 =============================================================
    # N = 
    f_L05 = args.f_L05
    # READ BY EYE
    #print(f" Original N={len(df_L05)} (L05)")
    # 2. Lim+2005 =============================================================

    # 3. Wolters+2008 =========================================================
    # N = 
    #df_W08 = make_df_W08()
    #print(f" Original N={len(df_W08)} (W08)")
    # 3. Wolters+2008 =========================================================

    # 4. SST/IRS ==============================================================
    # TODO: Check
    # N = 565
    # # Note: The spectra are used in Vernazza+2010. The details are not written in the paper.
    # In Vernazza+2010: They used spectra obtained from 2004-09-30 00:52 to 01:01.
    # downloaded from https://pds-smallbodies.astro.umd.edu/data_other/sptz_02_INNER/a433.shtml#top
    # All spectra were taken on 2004-09-30

    # ch0: 2 out of 12 spectra
    #   SPITZER_S0_4872960_0007_9_E7275841_tune.tbl
    #   SPITZER_S0_4872960_0010_9_E7275844_tune.tbl
    # ch2: all 4 spectra
    #   SPITZER_S2_4872960_0012_9_E7275703_tune.tbl 
    #   SPITZER_S2_4872960_0013_9_E7275707_tune.tbl
    #   SPITZER_S2_4872960_0014_9_E7275704_tune.tbl
    #   SPITZER_S2_4872960_0015_9_E7275696_tune.tbl
    f_dir = args.fdir
    f_list = [
        "SPITZER_S0_4872960_0007_9_E7275841_tune.tbl",
        "SPITZER_S0_4872960_0010_9_E7275844_tune.tbl",
        "SPITZER_S2_4872960_0012_9_E7275703_tune.tbl", 
        "SPITZER_S2_4872960_0013_9_E7275707_tune.tbl",
        "SPITZER_S2_4872960_0014_9_E7275704_tune.tbl",
        "SPITZER_S2_4872960_0015_9_E7275696_tune.tbl",
    ]
    utc_list = [
        "2004-09-30T00:56:14.834",
        "2004-09-30T00:57:36.623",
        "2004-09-30T00:58:58.419", 
        "2004-09-30T00:59:27.618",
        "2004-09-30T01:00:00.021",
        "2004-09-30T01:00:29.216",
    ]
    specidx_list = [8, 11, 1, 2, 3, 4]

    key_skip = ["\\processing", "\\wavsamp", "\\char", "\\character", "\\COMMENT", "\\HISTORY", "\\int", "\\float"]
    df_list = []
    for f in f_list:
        f0 = os.path.join(f_dir, f)
        # To get header once
        idx_h = 0
        dftemp_list = []
        with open(f0, "r") as f:
            lines = f.readlines()
            for l in lines:
                sl = l.split(" ")

                if sl[0] in key_skip:
                    pass
                else:
                    # Header
                    if (idx_h == 0) & (l[0] == "|"):
                        l_use = l.split("|")
                        # Remove "\n", ""
                        l_use = [x for x in l_use if x != "\n"]
                        keys = [x.strip() for x in l_use if x != ""]
                        idx_h = 1
                    elif l[0] == "|":
                        pass
                    # Data
                    else:
                        l_data = l.split()
                        # Check dimension
                        assert len(keys) == len(l_data)
                        # Make df
                        d = dict(zip(keys, l_data))
                        
                        df = pd.DataFrame(d.values(), index=d.keys()).T

                        # Save data
                        dftemp_list.append(df)

                      
                      
        # Make a dataframe
        df_1spec = pd.concat(dftemp_list)
        df_list.append(df_1spec)

    # # Do light-time correction, save results, and plot
    # Obs time
    # DATE_OBS: starting time
    #   2004-09-30T00:52:07.827
    # SAMPTIME: integtime?
    #   1.0486 s
    # EXPTOT_T: w/deadtime?
    #   14.68 s
    
    # TODO: hyper parameters
    # Relative error
    err_th = 0.05
    # Variation (25%)
    var_th = 0.25
    
    key_w, key_flux, key_fluxerr = "wavelength", "flux_density", "error"
    df4out_list = []
    for idx, df in enumerate(df_list):
        
        # Remove data with flag
        df = df[df["bit-flag"] == "0"]
        df = df.astype(float)
        # Remove large error
        df = df[df[key_fluxerr]/df[key_flux] < err_th]
        # Remove large variation
        df = remove_largevar(df, key_flux, var_th)
        # Remove edge
        df = remove_edge(df, 3)
        
        # Save
        utc_obs = utc_list[idx]
        # ltcor
        jd_ltcor = SST_ltcor(433, utc_obs)
        df = df.rename(columns={"flux_density":"flux", "error":"fluxerr"})
    
        df["jd"] = jd_ltcor
        df["cflag"] = 999
        if idx < 2:
            ch = 0
        else:
            ch = 2
        df["memo"] = f"SSTch{ch}_{specidx_list[idx]}"
        df["code"] = "@sst"
    
        # Columns:
        # jd wavelength flux fluxerr code cflag memo
        col4out = ["jd", "wavelength", "flux", "fluxerr", "code", "cflag", "memo"]
        df4out = df[col4out]
        df4out_list.append(df4out)

    # Merge 6 spectra
    df_S = pd.concat(df4out_list)
    print(f" Original N={len(df_S)} (SST)")
    # 4. SST/IRS ==============================================================
 

    # 5. AKARI/IRC ============================================================
    # N = 5 
    # Light-time not corrected
    # Color-term corrected (thus we put -9 and -8 as flags)
    f_AKARI = args.f_AKARI
    df_A = pd.read_csv(f_AKARI, sep=" ")
    print(f" Original N={len(df_A)} (AKARI)")
    # Do light-time correction
    jd_ltcor_list = []
    for jd in [2454199.2663881136, 2454199.404411632, 2454199.818465822, 2454199.8874780326, 2454199.9564906135]:
        ast = Horizons(location='500',id=433, epochs=jd)
        ast = ast.ephemerides()
        au_m = au.to("m").value
        lt_s = ast["delta"]*au_m/c
        lt_day = lt_s / 3600./24.
        jd_ltcor = jd - lt_day
        jd_ltcor_list.append(jd_ltcor)
    df_A["jd"] = jd_ltcor_list
    # 5. AKARI/IRC ============================================================


    # TODO: SNR cut here?


    # Make a single merged file
    df = pd.concat([df_H99, df_L05, df_W08, df_S, df_A])
    print(df)

    # Save 
