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

from Eros_common import mycolor, tel2col, tel2lab, plot_spec_Eros, generate_latex_table

if __name__ == "__main__":
    parser = ap(
        description="Plot location of Eros.")
    parser.add_argument(
        "flux", type=str,
        help="Preprocessed thermal fluxes")
    parser.add_argument(
        "--outdir", type=str, default=".",
        help="output directory")
    parser.add_argument(
        "--outtype", default="pdf",
        help="format of output figure")
    args = parser.parse_args()


    outdir = args.outdir
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    # Read all data
    # Outliers (negative fluxes) rejections were already done in the preprocesses
    f = args.flux
    df = pd.read_csv(f, sep=" ")
    print(f" Original N={len(df)}")

    # Plot fluxes
    tel_list = list(set(df.memo.values.tolist()))
    # Sort by year in ascending order
    tel_list = ["UKIRT1998", "UKIRT2002", "Lim2005_3", "SSTch0_2", "SSTch0_5", "SSTch0_8", 
                "SSTch0_11", "SSTch2_3", "SSTch2_4", "SSTch2_1", "SSTch2_2", "akari"]


    # Make a reference column
    df["ref"] = 0
    df.loc[(df['memo'] == 'UKIRT1998'), 'ref'] = "Harris et al. (1999)"
    df.loc[(df['memo'] == 'UKIRT2002'), 'ref'] = "Wolters et al. (2008)"
    df.loc[(df['memo'] == 'Lim2005_3'), 'ref'] = "Lim et al. (2005)"
    df.loc[(df['memo'] == 'SSTch0_2'), 'ref'] = "SST/IRS"
    df.loc[(df['memo'] == 'SSTch0_5'), 'ref'] = "SST/IRS"
    df.loc[(df['memo'] == 'SSTch0_8'), 'ref'] = "SST/IRS"
    df.loc[(df['memo'] == 'SSTch0_11'), 'ref'] = "SST/IRS"
    df.loc[(df['memo'] == 'SSTch2_1'), 'ref'] = "SST/IRS"
    df.loc[(df['memo'] == 'SSTch2_2'), 'ref'] = "SST/IRS"
    df.loc[(df['memo'] == 'SSTch2_3'), 'ref'] = "SST/IRS"
    df.loc[(df['memo'] == 'SSTch2_4'), 'ref'] = "SST/IRS"
    df.loc[(df['memo'] == 'akari'), 'ref'] = "AKARI/IRC"


    out = "Eros_fig_flux.pdf"
    out = os.path.join(args.outdir, out)
    plot_spec_Eros(df, tel_list, out)

    latex_tab = generate_latex_table(df)
    print(latex_tab)
