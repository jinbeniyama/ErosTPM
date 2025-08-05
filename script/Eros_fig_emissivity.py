#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plot emissivity of Eros.
"""
import os 
from argparse import ArgumentParser as ap
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from Eros_common import mycolor, mymark



if __name__ == "__main__":
    parser = ap(
        description="Plot emissivity of Eros.")
    parser.add_argument(
        "eps", type=str,
        help="Preprocessed emissivity")
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
    df = pd.read_csv(args.eps, sep=" ")
    print(f" Original N={len(df)}")

    out = "Eros_fig_emissivity.pdf"
    out = os.path.join(args.outdir, out)

    fig = plt.figure(figsize=(8, 12))
    ax = fig.add_axes([0.15, 0.1, 0.8, 0.85])
    ax.set_xlabel("Wavelength [micron]")
    ax.set_ylabel("Emissivity (observation/model)")

    # Emissivity
    label = "Eros (Spitzer)"
    ax.errorbar(
        df["w"], df["emissivity"], df["emissivityerr"], fmt="o", color="black", 
        label=label, ms=10)
    # Running mean
    label = "Eros (Spitzer) running mean"
    ax.errorbar(
        df["w"], df["emissivity_smo"], df["emissivityerr_smo"], fmt="x", color="red", 
        label=label, ms=5)

    ax.grid()
    ax.legend()
    ax.set_ylim([0.8, 1.2])

    plt.savefig(out)
