#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Useful functions for Eros paper.
"""
import os
import pandas as pd
import numpy as np
import datetime
from astroquery.jplhorizons import Horizons
import matplotlib.pyplot as plt  
from astropy import units as u
from astropy.constants import c, au
from astropy.time import Time


mycolor = [
    "#AD002D", "#1e50a2", "#69821b", "#f055f0", "#afafb0", 
    "#0095b9", "#89c3eb", "#ec6800", "cyan", "gold", "magenta"
    ] 
mycolor = mycolor*10

mymark = ["o", "^", "s", "D", "*", "v", "<", ">", "h", "H"]
mymark = mymark*100


def tel2col(tel):
    """Define color for each telescope.

    Parameter
    ---------
    tel : str
        telescope name

    Return
    ------
    col : str
        color
    """
    if tel == "UKIRT1998":
        col = mycolor[0]
    elif tel == "UKIRT2002":
        col = mycolor[1]
    elif tel == "Lim2005_3":
        col = mycolor[2]
    elif tel in [
        "SSTch2_1_and_2", 
        "SSTch2_3_and_4"]:
        col = mycolor[3]
    elif tel == "akari":
        col = mycolor[5]
    else:
        col = "black"
    return col


def tel2lab(tel):
    """Define label for each telescope.

    Parameter
    ---------
    tel : str
        telescope name

    Return
    ------
    lab : str
        label
    """
    if tel == "UKIRT1998":
        lab = "1998-06-27,28,29,30 (Harris+1999)\nN=175"
    elif tel == "UKIRT2002":
        lab = "2002-09-28 (Wolters+2008)\nN=13"
    elif tel == "Lim2005_3":
        lab = "2002-09-21,22 (Lim+2005)\nN=53"
    elif tel in [
        "SSTch2_1_and_2", "SSTch2_3_and_4"]:
        # Only four, averaged
        lab = "2004-09-30 (SST)\nN=202"
    elif tel == "akari":
        lab = "2007-04-08(AKARI)\nN=5"
    else:
        lab = "label"
    return lab


def plot_spec_Eros(df, tel_list, out=None):
    """Plot spectra of Eros.

    Parameters
    ----------
    df : pandas.DataFrame
        data frame with spectra
    tel_list : array-like
        list of telescopes
    out : str, optional 
        output filename
    """
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_axes([0.10, 0.12, 0.85, 0.85])
    ax.set_xlabel("Wavelength [micron]")
    ax.set_ylabel("Flux density [Jy]")
    ax.set_xlim([5, 40])
    ax.set_ylim([0, 15])

    # Labels already used
    lab_used = []
    for tel in tel_list:
        df_tel = df[df["memo"] == tel]
        N_tel = len(df_tel)
        col = tel2col(tel)
        lab = tel2lab(tel)
        if lab in lab_used:
            lab = None
        else:
            lab_used.append(lab)
        ax.errorbar(
            df_tel["wavelength"], df_tel["flux"],  df_tel["fluxerr"], 
            marker="o", ls="None", ms=5, color=col, label=lab)
    ax.legend()
    if out:
        plt.savefig(out)
    plt.close()


def generate_latex_table(df):
    """Generate latex table of observations.

    Parameter
    ---------
    df : pandas.DataFrame
        dataframe with observations

    Return
    ------
    latex_table
        latex table
    """
    latex_table = r"""
\begin{longtable}{ccccc}
\centering
\caption{Observed fluxes at different epochs and wavelengths.} \label{tab:flux} \\
\hline\hline
Epoch & Wavelength & Flux [Jy] & Flux error [Jy] & Reference\\
(JD) & ($\mu$m) &  &  &\\
\hline
\endfirsthead

\hline\hline
Epoch & Wavelength & Flux [Jy] & Flux error [Jy] & Reference\\
(JD) & ($\mu$m) &  &  & \\
\hline
\endhead

\hline
\endfoot

\hline
"""
    for _, row in df.iterrows():
        latex_table += f"{row['jd']:.5f} & {row['wavelength']:.1f} & {row['flux']:.2f} & {row['fluxerr']:.2f} & {row['ref']}\\\\\n"
    latex_table += r"""\hline
\end{longtable}
\tablefoot{The errors represent 1$\sigma$ uncertainties.}
"""
    return latex_table


def Eros_Harris1999(f):
    """Handle Eros observations in Harris+1999.

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
                    df = pd.DataFrame(
                        dict(jd=jd_list, wavelength=w_list, flux=f_list, 
                        fluxerr=ferr_list, n=idx_list))
                    df_list.append(df)
    df = pd.concat(df_list)
    df["code"] = 568
    df["cflag"] = 999
    df["memo"] = "UKIRT1998"
    return df


def Eros_Lim2005_3(f):
    """Handle Eros observations in Lim+2005.

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

    # Just to check the consistency with Figure in the paper
    #df["flux_W"] = df["flux"]/1e19/(df["wavelength"]*1e-6)**2*c.value
    #df["fluxerr_W"] = df["fluxerr"]/1e19/(df["wavelength"]*1e-6)**2*c.value
    
    df["jd"] = float(epoch_ltcor)
    df["code"] = 675
    df["cflag"] = 999
    df["memo"] = "Lim2005_3"
    return df


def Eros_Wolters2008():
    """Handle Eros observations in Wolters+2008.
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

    # Add error
    
    
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
    """Remove ourliers.

    Parameters
    ----------
    df : pandas.DataFrame
        input dataframe with thermal fluxes
    key : str
        keyword to be checked
    var_th : float
        variation threshold.
        < 1-var_th and > 1-var_th are removed.

    Return
    ------
    df : pandas.DataFrame
        output dataframe without outliers
    """
    # Remove large variation
    df = df.reset_index(drop=True)
    idx_list = []

    for idx, row in df.iterrows():

        if idx == 0:
            val0 = row[key]
            valmin = val0*(1 - var_th)
            valmax = val0*(1 + var_th)
            idx_list.append(idx)
            #print(valmin)
            #print(valmax)
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

    #print(idx_list)
    df = df.loc[idx_list]
    return df


def remove_edge(df, n):
    """Remove data close to edge.

    Parameters
    ----------
    df : pandas.DataFrame
        input dataframe with thermal fluxes
    n : int
        number of data ignored on each side

    Return
    ------
    df : pandas.DataFrame
        output dataframe without outliers
    """
    df = df.reset_index(drop=True)
    N = len(df)
    idx_list = df.index.values.tolist()
    idx_use = idx_list[n:N-n]

    df = df.loc[idx_use]
    return df


def SST_ltcor(target, jd):
    """Light-time correction for SST data.

    Parameters
    ----------
    target : str
        target name
    jd : str
        julian day

    Return
    ------
    jd_ltcor : float 
        light-time corrected time in jd
    """
    ## Lighttime correction
    ast = Horizons(location='@sst',id=target, epochs=jd)
    ast = ast.ephemerides()
    au_m = au.to("m").value
    lt_s = ast["delta"]*au_m/c
    lt_day = lt_s / 3600./24.
    jd_ltcor = jd - lt_day
    #print(f"  epoch_jd       : {jd}")
    #print(f"  epoch_jd_ltcor : {jd_ltcor.value[0]}")
    jd_ltcor = jd_ltcor.value[0]
    return jd_ltcor
