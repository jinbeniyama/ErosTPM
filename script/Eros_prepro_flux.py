#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Preprocesses of thermal fluxes of Eros.
N = 811 (= 175+53+13+565+5)

1. Eros in 1998 published in Harris+1999
   Q-band spectra are not used 
   due to the relatively large uncertainty in the absolute calibration.
2. Eros in 2002 sep. published in Wolters+2008
3. Eros in 2002 sep. published in Lim+2005
4. Eros in 2004 by SST
5. Eros in 2007 by AKARI
"""
import os 
from argparse import ArgumentParser as ap
import pandas as pd
from astroquery.jplhorizons import Horizons
from astropy.constants import c, au
from astropy.time import Time


def Eros_Harris1999(f):
    """

    Parameter
    ---------
    f : str
        Spectra of Eros in Harris+1999

    Return
    ------
    df : pandas.DataFrame
        dataframe
    """
    ## UKIRT/Michelle in 1998 June
    # Originally provided "Eros_UKIRT_June_1998.txt" by Alan Harris. Thanks!
    # The given times refer (roughly!) to the mid-point of the set of exposures taken on each date.

    # Wavelength [micron], Flux [mJy], Flux error [mJy]
    # June 27 at 07:00: 2 spectra (8--13 micron)
    # June 29 at 06:00: 2 spectra (8--13 micron)
    # June 30 at 07:00: 3 spectra (8--13 micron)
    # June 28 at 06:00: 2 spectra (16--24 micron) -> these are used as well
    
    # A bit more searching has led me to the following timing information:
    # UT Date/time    Obs. no.
    # 1998 June                    
    # 27/0632      22  
    #    0733      34
    # 29/0549      27
    #    0557      28
    # 30/0553     25
    #    0708     38
    #    0808     45
    # 28/0553      17       
    #    0609      19
    idx_spec = 0
    df_list = []
    with open(f, "r") as f:
        lines = f.readlines()
        for idx_l, l in enumerate(lines):
            if l[0:4] == "1998":
                # make list
                jd_list = []
                w_list, f_list, ferr_list, idx_list = [], [], [], []
                idx_spec += 1
                t = Time(str(l), format='isot', scale='utc')
                epoch = t.jd
                
                # Light time correction
                ast = Horizons(location='568',id=433, epochs=epoch)
                ast = ast.ephemerides()
                au_m = au.to("m").value
                lt_s = ast["delta"]*au_m/c
                lt_day = lt_s / 3600./24.
                epoch_ltcor = epoch - lt_day
            else:
                w, f, ferr, _ = l.split()
                jd_list.append(float(epoch_ltcor))
                w_list.append(float(w))
                # mJy to Jy
                f_list.append(float(f)*1e-3)
                ferr_list.append(float(ferr)*1e-3)
                idx_list.append(int(idx_spec))
    
                # Check whether the last line or not
                # Last or     
                if (idx_l == len(lines)-1) or (lines[idx_l + 1][0:4] == "1998"):
                    df = pd.DataFrame(dict(jd=jd_list, wavelength=w_list, flux=f_list, fluxerr=ferr_list, n=idx_list))
                    df_list.append(df)
    df = pd.concat(df_list)
    df["code"] = 568
    df["cflag"] = 999
    df["memo"] = "UKIRT1998"
    return df


def Eros_Lim2005_3(f):
    """

    Parameter
    ---------
    f : str
        Spectra of Eros in Lim+2005 read by eye

    Return
    ------
    df : pandas.DataFrame
        dataframe
    """
    # (2002-09-21)
    # (2002-09-22 A 1st)
    # 2002-09-22 B 2nd
    # TODO: There are two other spectra.

    ## 2002-09-22 2nd spectrum
    ## 6:16:23--6:46:28
    ## 2452539.761377315--2452539.7822685186
    epoch = (2452539.761377315 + 2452539.7822685186)/2.
    
    ## Lighttime correction
    ast = Horizons(location='675',id=433, epochs=epoch)
    ast = ast.ephemerides()
    au_m = au.to("m").value
    lt_s = ast["delta"]*au_m/c
    lt_day = lt_s / 3600./24.
    epoch_ltcor = epoch - lt_day
    
    w_list, f_list, ferr_list = [], [], []
    with open(f, "r") as f:
        lines = f.readlines()[2:]
        for idx_l, l in enumerate(lines):
            w, f, ferr = l.split()
            w_list.append(float(w))
            f_list.append(float(f))
            ferr_list.append(float(ferr))
    
    df = pd.DataFrame(dict(wavelength=w_list, flux=f_list, fluxerr=ferr_list))
    
    df["jd"] = float(epoch_ltcor)
    df["code"] = 675
    df["cflag"] = 999
    df["memo"] = "Lim2005_3"
    return df


def Eros_Wolters2008():
    """
    """
    ## UKIRT/Michelle
    ## Obs time = 2002-09-28 09:39–09:58
    ## Central time 2002-09-28 09:50 = 2452545.909722222
    ## Wavelength [micron], Flux density [10^-13 W/m^2/micron], its uncertainty [10^-13 W/m^2/micron]
    w_list = [
        8.120, 8.369, 8.625, 8.884, 9.159, 10.166, 10.476, 10.783, 11.088, 11.393, 11.700, 12.011, 12.329]
    f_list = [
        2.86,  2.99,  3.07, 3.17, 3.15, 3.22, 3.19, 3.12, 3.08, 3.06, 3.05, 2.99, 2.84]
    ferr_list = [
        0.032, 0.020, 0.022, 0.012, 0.013, 0.0097, 0.0089, 0.014, 0.012, 0.0066, 0.010, 0.011, 0.021]
    df = pd.DataFrame(dict(wavelength=w_list, f=f_list, ferr=ferr_list))
    
    
    # Constant to convert W/m^2/m to Jy (see note on iPad)
    df["f"] *= 1e-13
    df["ferr"] *= 1e-13
    # x [W/m^2/m] = x*const*w**2 [Jy]
    const = 1e20/c.value
    # w, wavelength: micron
    # c: m/s
    df["flux"] = df["f"]*df["wavelength"]**2*const
    df["fluxerr"] = df["ferr"]*df["wavelength"]**2*const
    
    ## Lighttime correction
    epoch = 2452545.909722222
    ast = Horizons(location='568',id=433, epochs=epoch)
    
    ast = ast.ephemerides()
    au_m = au.to("m").value
    lt_s = ast["delta"]*au_m/c
    lt_day = lt_s / 3600./24.
    epoch_ltcor = epoch - lt_day
    
    df["jd"] = float(epoch_ltcor)
    df["code"] = 568
    df["cflag"] = 999
    df["memo"] = "UKIRT2002"
    return df

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
        "--outdir", type=str, default="../data",
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
    print(f" Original N={len(df_H99)} (H99)")
    print(f"  Columns: {df_H99.columns.tolist()}")
    # 1. Harris+1999 ==========================================================

    # 2. Lim+2005 =============================================================
    # N = 53
    # READ BY EYE
    f_L05 = args.f_L05
    df_L05 = Eros_Lim2005_3(f_L05)
    df_L05 = df_L05[col4out]
    print(f" Original N={len(df_L05)} (L05)")
    print(f"  Columns: {df_L05.columns.tolist()}")
    # 2. Lim+2005 =============================================================

    # 3. Wolters+2008 =========================================================
    # N = 13
    df_W08 = Eros_Wolters2008()
    df_W08 = df_W08[col4out]
    print(f" Original N={len(df_W08)} (W08)")
    print(f"  Columns: {df_W08.columns.tolist()}")
    # 3. Wolters+2008 =========================================================

    # 4. SST/IRS ==============================================================
    # TODO: Check
    # N = 565
    # # Note: The spectra are used in Vernazza+2010. The details are not written in the paper.
    # In Vernazza+2010: They used spectra obtained from 2004-09-30 00:52 to 01:01.
    # downloaded from https://pds-smallbodies.astro.umd.edu/data_other/sptz_02_INNER/a433.shtml#top

    # ch0: 4 out of 12 spectra
    # ch2: all 4 spectra
    f_dir = args.fdir
    f_list = [
        "SPITZER_S0_4872960_0001_9_E7275827_tune.tbl",
        "SPITZER_S0_4872960_0004_9_E7275831_tune.tbl",
        "SPITZER_S0_4872960_0007_9_E7275841_tune.tbl",
        "SPITZER_S0_4872960_0010_9_E7275844_tune.tbl",
        "SPITZER_S2_4872960_0012_9_E7275703_tune.tbl", 
        "SPITZER_S2_4872960_0013_9_E7275707_tune.tbl",
        "SPITZER_S2_4872960_0014_9_E7275704_tune.tbl",
        "SPITZER_S2_4872960_0015_9_E7275696_tune.tbl",
    ]
    utc_list = [
        "2004-09-30T00:52:43.624",
        "2004-09-30T00:54:32.420",
        "2004-09-30T00:56:14.834",
        "2004-09-30T00:57:36.623",
        "2004-09-30T00:58:58.419", 
        "2004-09-30T00:59:27.618",
        "2004-09-30T01:00:00.021",
        "2004-09-30T01:00:29.216",
    ]
    specidx_list = [2, 5, 8, 11, 1, 2, 3, 4]

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
        if idx < 4:
            ch = 0
        else:
            ch = 2
        df["memo"] = f"SSTch{ch}_{specidx_list[idx]}"
        df["code"] = "@sst"
    
        df4out = df[col4out]
        df4out_list.append(df4out)

    # Merge 6 spectra
    df_S = pd.concat(df4out_list)
    print(f" Original N={len(df_S)} (SST)")
    print(f"  Columns: {df_S.columns.tolist()}")
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
    df_A.loc[:, "jd"] = jd_ltcor_list
    print(f"  Columns: {df_A.columns.tolist()}")
    # 5. AKARI/IRC ============================================================


    # Make a single merged file
    df = pd.concat([df_H99, df_L05, df_W08, df_S, df_A])
    # Save 
    out = f"Eros_flux_N{len(df)}.txt"
    out = os.path.join(args.outdir, out)
    df.to_csv(out, sep=" ", index=False)
