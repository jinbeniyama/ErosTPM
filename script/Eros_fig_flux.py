#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plot thermal fluxes of Eros, and also print them.

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

from myplot import mycolor

def tel2col(tel):
    """
    Define color for a telescope.
    """
    if tel == "UKIRT1998":
        col = mycolor[0]
    elif tel == "UKIRT2002":
        col = mycolor[1]
    elif tel == "Lim2005_3":
        col = mycolor[2]
    elif tel in ["SSTch0_8", "SSTch0_11", "SSTch2_1",  "SSTch2_2", "SSTch2_3", "SSTch2_4"]:
        col = mycolor[3]
    elif tel == "akari":
        col = mycolor[5]
    else:
        col = "black"
    return col


def tel2lab(tel):
    """
    Define label for a telescope.
    """

    if tel == "UKIRT1998":
        lab = "1998 Jun. (Harris+1999)"
    elif tel == "UKIRT2002":
        lab = "2002 Sep. (Wolters+2009)"
    elif tel == "Lim2005_3":
        lab = "2002 Sep. (Lim+2005)"
    elif tel in ["SSTch0_8", "SSTch0_11", "SSTch2_1",  "SSTch2_2", "SSTch2_3", "SSTch2_4"]:
        lab = "2004 Sep. (SST)"
    elif tel == "akari":
        lab = "2007 Apr. (AKARI)"
    else:
        lab = "label"
        
    return lab



def plot_spec_Eros(df, tel_list, out=None):


    #set_gnuplot_style()

    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_axes([0.10, 0.12, 0.85, 0.85])
    ax.set_xlabel("Wavelength [micron]")
    ax.set_ylabel("Flux density [Jy]")

    # Labels already used
    lab_used = []
    for tel in tel_list:
        df_tel = df[df["memo"] == tel]
        col = tel2col(tel)
        lab = tel2lab(tel)
        if lab in lab_used:
            lab = None
        else:
            lab_used.append(lab)
            
        ax.errorbar(df_tel["wavelength"], df_tel["flux"],  df_tel["fluxerr"], marker="o", ls="None", ms=5, color=col, label=lab)

    ax.legend()
   
    if out:
        plt.savefig(out)


def generate_latex_table(df):
    latex_table = r"""
\begin{longtable}{cccc}
\centering
\caption{Observed fluxes at different epochs and wavelengths.} \label{tab:flux} \\
\hline\hline
Epoch & Wavelength & Flux [Jy] & Flux error [Jy] \\
(JD) & ($\mu$m) &  &  \\
\hline
\endfirsthead

\hline\hline
Epoch & Wavelength & Flux [Jy] & Flux error [Jy] \\
(JD) & ($\mu$m) &  &  \\
\hline
\endhead

\hline
\endfoot

\hline
"""

    for _, row in df.iterrows():
        latex_table += f"{row['jd']:.1f} & {row['wavelength']:.1f} & {row['flux']:.2f} & {row['fluxerr']:.2f} \\\\\n"

    latex_table += r"""\hline
\end{longtable}
\tablefoot{The errors represent 1$\sigma$ uncertainties.}
"""

    return latex_table

if __name__ == "__main__":
    parser = ap(
        description="Plot location of Eros.")
    parser.add_argument(
        "flux", type=str,
        help="Preprocessed thermal fluxes")
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

    # Read all data
    f = args.flux
    df = pd.read_csv(f, sep=" ")
    print(f" Original N={len(df)}")

    # TODO: Check data selection
    # Remove outliers
    df = df[df["flux"] > 0]
    print(f" Remove negative flux -> N={len(df)}")
    # Remove low SNR
    SNR_min = 5
    df = df[df["flux"]/df["fluxerr"] > SNR_min]
    print(f" Remove SNR < {SNR_min} -> N={len(df)}")
    
    
    # Plot fluxes

    tel_list = list(set(df.memo.values.tolist()))
    # Sort by year in ascending order
    tel_list = ["UKIRT1998", "UKIRT2002", "Lim2005_3", "SSTch2_3", "SSTch0_8", "SSTch0_11", "SSTch2_4", "SSTch2_1", "SSTch2_2", "akari"]

    out = "Eros_fig_flux.pdf"
    out = os.path.join(args.outdir, out)
    plot_spec_Eros(df, tel_list, out)





    latex_tab = generate_latex_table(df)
    print(latex_tab)


