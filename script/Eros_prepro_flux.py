#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Preprocesses of thermal fluxes of Eros.

N = 426 (= 175+13+53+180+5)

1. Eros in 1998 published in Harris+1999
   Q-band spectra are not used 
   due to the relatively large uncertainty in the absolute calibration.
2. Eros in 2002 Sep. published in Wolters+2008
3. Eros in 2002 Sep. published in Lim+2005
4. Eros in 2004 Sep. by SST
5. Eros in 2007 Apr. by AKARI (Usui+2011)

Note
----
- The time in the output files is light-time corrected.
- Q-band spectroscopy in Harris+1999 is not used in the project.
- Spitzer/IRAC photometry (at 3.6 & 4.5 micron, Trilling+2010, Mueller+2011) 
  was not used in the project. 
"""
import os 
from argparse import ArgumentParser as ap
import pandas as pd
import numpy as np
from astroquery.jplhorizons import Horizons
from astropy.constants import c, au
from astropy.time import Time

from Eros_common import (
    Eros_Harris1999, Eros_Lim2005_3, Eros_Wolters2008, remove_largevar, 
    remove_edge, SST_ltcor)


def make_ave_SST(df1, df2):
    """Make avareged spectrum for SST (of Eros).

    Parameters
    ----------
    df1 : pandas.DataFrame
        input data1
    df2 : pandas.DataFrame
        input data2

    Return
    ------
    df12 : pandas.DataFrame
        output merged data
    """

    # Most of the wavelengths are the same, but not always all of them.
    w1 = sorted(set(df1.wavelength))
    w2 = sorted(set(df2.wavelength))

    # Final wavelength grids
    w12 = sorted(set(w1 + w2))

    print(f"    Wavelength1  N = {len(w1)}")
    print(f"    Wavelength2  N = {len(w2)}")
    print(f"    -> combined  N = {len(w12)}\n")

    # Set wavelength as index
    df1_ = df1.set_index("wavelength")
    df2_ = df2.set_index("wavelength")

    # Make a new dataframe
    f_mean = []
    ferr_mean = []

    for w in w12:
        f1 = df1_["flux"][w] if w in df1_.index else np.nan
        f2 = df2_["flux"][w] if w in df2_.index else np.nan
        e1 = df1_["fluxerr"][w] if w in df1_.index else np.nan
        e2 = df2_["fluxerr"][w] if w in df2_.index else np.nan

        # Mean
        if not np.isnan(f1) and not np.isnan(f2):
            f = 0.5 * (f1 + f2)
            # Error
            ferr_ave = 0.5 * np.abs(f1 - f2) 
            # Original error
            ferr_ori = (0.5*(e1**2 + e2**2))**0.5
            ferr = (ferr_ave**2 + ferr_ori**2)**0.5
        elif not np.isnan(f1):
            f = f1
            ferr = e1
        elif not np.isnan(f2):
            f = f2
            ferr = e2
        else:
            f = np.nan
            ferr = np.nan
            assert False, "Check the code."

        f_mean.append(f)
        ferr_mean.append(ferr)

    df_12 = pd.DataFrame({
        "wavelength": w12,
        "flux": f_mean,
        "fluxerr": ferr_mean
    })
    return df_12



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
        "--fdir", type=str, default="data",
        help="input file directory")
    parser.add_argument(
        "--outdir", type=str, default="data",
        help="output directory")
    args = parser.parse_args()


    outdir = args.outdir
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    # Read input files
    # Columns used in the final output:
    # jd wavelength flux fluxerr code cflag memo
    col4out = ["jd", "wavelength", "flux", "fluxerr", "code", "cflag", "memo"]

    # 1. Harris+1999 ==========================================================
    # N = 175 (FIY, N = 303 incl. Q-band)
    f_H99 = args.f_H99
    df_H99 = Eros_Harris1999(f_H99)
    # Use only N-band
    df_H99 = df_H99[df_H99["wavelength"] < 15]
    df_H99 = df_H99[col4out]
    print(f"Original N={len(df_H99)} (H99)")
    print(f"    Columns: {df_H99.columns.tolist()}\n")
    # 1. Harris+1999 ==========================================================

    # 2. Lim+2005 =============================================================
    # N = 53
    # READ BY EYE
    f_L05 = args.f_L05
    df_L05 = Eros_Lim2005_3(f_L05)
    df_L05 = df_L05[col4out]
    # Set error as 10%
    df_L05["fluxerr"] = df_L05["flux"]*0.1

    print(f"Original N={len(df_L05)} (L05)")
    print(f"    Columns: {df_L05.columns.tolist()}\n")
    # 2. Lim+2005 =============================================================

    # 3. Wolters+2008 =========================================================
    # N = 13
    df_W08 = Eros_Wolters2008()
    df_W08 = df_W08[col4out]
    print(f"Original N={len(df_W08)} (W08)")
    print(f"    Columns: {df_W08.columns.tolist()}\n")
    # 3. Wolters+2008 =========================================================

    # 4. SST/IRS ==============================================================
    # N = 426
    # # Note: The spectra are used in Vernazza+2010. The details are not written in the paper.
    # In Vernazza+2010: They used spectra obtained from 2004-09-30 00:52 to 01:01.
    # downloaded from https://pds-smallbodies.astro.umd.edu/data_other/sptz_02_INNER/a433.shtml#top

    # ch0: -> Not use!
    # ch2: all 4 spectra
    f_dir = args.fdir
    f_list = [
        "SPITZER_S2_4872960_0012_9_E7275703_tune.tbl", 
        "SPITZER_S2_4872960_0013_9_E7275707_tune.tbl",
        "SPITZER_S2_4872960_0014_9_E7275704_tune.tbl",
        "SPITZER_S2_4872960_0015_9_E7275696_tune.tbl",
    ]
    utc_list = [
        "2004-09-30T00:58:58.419", 
        "2004-09-30T00:59:27.618",
        "2004-09-30T01:00:00.021",
        "2004-09-30T01:00:29.216",
    ]
    # Not used
    key_skip = [
        "\\processing", "\\wavsamp", "\\char", "\\character", 
        "\\COMMENT",  "\\HISTORY", "\\int", "\\float"]
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
    
    # Merge two spectra in each wavelength region
    assert len(df_list) == 4, "Check the code"

    df_SST1 = df_list[0]
    df_SST2 = df_list[1]
    df_SST3 = df_list[2]
    df_SST4 = df_list[3]

    print("Original (SST, four spectra)")
    print(f"    N1 = {len(df_SST1)}")
    print(f"    N2 = {len(df_SST2)}")
    print(f"    N3 = {len(df_SST3)}")
    print(f"    N4 = {len(df_SST4)}\n")

    # Remove data with flag
    df_SST1 = df_SST1[df_SST1["bit-flag"] == "0"]
    df_SST1 = df_SST1.astype(float)
    df_SST2 = df_SST2[df_SST2["bit-flag"] == "0"]
    df_SST2 = df_SST2.astype(float)
    df_SST3 = df_SST3[df_SST3["bit-flag"] == "0"]
    df_SST3 = df_SST3.astype(float)
    df_SST4 = df_SST4[df_SST4["bit-flag"] == "0"]
    df_SST4 = df_SST4.astype(float)

    print("After removal with flag ")
    print(f"    N1 = {len(df_SST1)}")
    print(f"    N2 = {len(df_SST2)}")
    print(f"    N3 = {len(df_SST3)}")
    print(f"    N4 = {len(df_SST4)}\n")


    # # Do light-time correction, save results, and plot
    # Obs time
    # DATE_OBS: starting time
    #   2004-09-30T00:52:07.827
    # SAMPTIME: integtime?
    #   1.0486 s
    # EXPTOT_T: w/deadtime?
    #   14.68 s
    
    # Relative error
    err_th = 0.05
    # Variation (25%)
    var_th = 0.25
    
    key_w, key_flux, key_fluxerr = "wavelength", "flux_density", "error"
    df_list_cor = []
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
        df = df.rename(columns={"flux_density":"flux", "error":"fluxerr"})
        df_list_cor.append(df)

    # Merge 2 spectra in each wavelength region
    df_SST1 = df_list_cor[0]
    df_SST2 = df_list_cor[1]
    df_SST3 = df_list_cor[2]
    df_SST4 = df_list_cor[3]
    
    # Get mean spectrum
    df_S_12 = make_ave_SST(df_SST1, df_SST2)
    df_S_34 = make_ave_SST(df_SST3, df_SST4)
    print(f" Merged N={len(df_S_12)} (SST 14--20 micron)")
    print(f" Merged N={len(df_S_34)} (SST 20--39 micron)\n")

    # Mean utc1 and utc2
    utc1 = utc_list[0]
    jd1 = Time(utc1, format="isot", scale="utc").jd
    utc2 = utc_list[1]
    jd2 = Time(utc2, format="isot", scale="utc").jd
    jd12 = np.mean([jd1, jd2])
    jd12_ltcor = SST_ltcor(433, jd12)
    df_S_12["jd"] = jd12_ltcor
    df_S_12["cflag"] = 999
    df_S_12["memo"] = f"SSTch2_1_and_2"
    df_S_12["code"] = "@sst"

    # Mean utc3 and utc4
    utc3 = utc_list[2]
    jd3 = Time(utc3, format="isot", scale="utc").jd
    utc4 = utc_list[3]
    jd4 = Time(utc4, format="isot", scale="utc").jd
    jd34 = np.mean([jd3, jd4])
    jd34_ltcor = SST_ltcor(433, jd34)
    df_S_34["jd"] = jd34_ltcor
    df_S_34["cflag"] = 999
    df_S_34["memo"] = f"SSTch2_3_and_4"
    df_S_34["code"] = "@sst"
    
    # Set relative error as 5% (see Hanus+2016)
    df_S_12["fluxerr"] = df_S_12["flux"]*0.05
    df_S_34["fluxerr"] = df_S_34["flux"]*0.05
    # 4. SST/IRS ==============================================================
 

    # 5. AKARI/IRC ============================================================
    # N = 5 
    # Light-time not corrected
    # Color-term corrected (thus we put -9 and -8 as flags)
    f_AKARI = args.f_AKARI
    df_A = pd.read_csv(f_AKARI, sep=" ")
    print(f"Original N={len(df_A)} (AKARI)")
    # Do light-time correction
    jd_ltcor_list = []
    for jd in [2454199.2663881136, 2454199.404411632, 2454199.818465822, 
               2454199.8874780326, 2454199.9564906135]:
        # Use code 500
        ast = Horizons(location='500',id=433, epochs=jd)
        ast = ast.ephemerides()
        au_m = au.to("m").value
        lt_s = ast["delta"]*au_m/c
        lt_day = lt_s / 3600./24.
        jd_ltcor = jd - lt_day
        jd_ltcor_list.append(jd_ltcor)
    df_A.loc[:, "jd"] = jd_ltcor_list
    print(f"    Columns: {df_A.columns.tolist()}\n")
    # 5. AKARI/IRC ============================================================

    # Make a single merged file
    df = pd.concat([df_H99, df_L05, df_W08, df_S_12, df_S_34, df_A])
    # Save 
    out = f"Eros_flux_N{len(df)}.txt"
    out = os.path.join(args.outdir, out)
    df.to_csv(out, sep=" ", index=False)
    
    N_total = len(df)
    print(f"Total number of obs. N={N_total}")
