#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plot thermal fluxes by SST.
"""
import os 
from argparse import ArgumentParser as ap
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time

from Eros_common import mycolor, mymark


def plot_SST_flux(df, label, color="black", out=None):
    fig = plt.figure(figsize=(16, 12))
    ax1 = fig.add_axes([0.12, 0.15, 0.8, 0.8])
    ax1.set_xlabel("Wavelength [micron]")
    ax1.set_ylabel("Flux density [Jy]")
    ax1.set_title(label)
    
    key_w, key_flux, key_fluxerr = "wavelength", "flux", "fluxerr"

    epoch_list = sorted(set(df.jd))
     
    for i, epoch in enumerate(epoch_list):
        df_epo = df[df["jd"] == epoch]
        utc = Time(str(epoch), format='jd', scale='utc').datetime.replace(microsecond=0)

        label = f"Epoch {epoch:.3f} {utc}"
        ax1.errorbar(
            df_epo[key_w], df_epo[key_flux], yerr=df_epo[key_fluxerr],
            fmt=mymark[i], ls="None", color=mycolor[i], ms=15, markerfacecolor='None', label=label)
    
    ax1.grid(which="both")
    ax1.set_yscale("log")
    
    ax1.legend(fontsize=10, ncol=2)
    if out:
        plt.savefig(out)


def plot_SST_emissivity(df, label, color="black", out=None):
    fig = plt.figure(figsize=(16, 12))
    ax1 = fig.add_axes([0.12, 0.15, 0.8, 0.8])
    ax1.set_xlabel("Wavelength [micron]")
    ax1.set_ylabel("Emissivity (obs/model)")
    ax1.set_title(label)
    
    key_w, key_flux, key_fluxerr = "wavelength", "flux", "fluxerr"

    epoch_list = sorted(set(df.jd))
     
    for i, epoch in enumerate(epoch_list):
        df_epo = df[df["jd"] == epoch]
        utc = Time(str(epoch), format='jd', scale='utc').datetime.replace(microsecond=0)
        label = f"Epoch {epoch:.3f} {utc}"
        
        # Fit by polynomial
        degree = 3
        wavelengths = df_epo[key_w].values
        fluxes = df_epo[key_flux].values
        coeffs = np.polyfit(wavelengths, fluxes, deg=degree)
        poly_model = np.poly1d(coeffs)
        df_epo["flux_model"] = poly_model(wavelengths)
        df_epo["emissivity"] = df_epo[key_flux]/df_epo["flux_model"]

        ax1.errorbar(
            df_epo[key_w], df_epo["emissivity"], yerr=df_epo[key_fluxerr]/df_epo[key_flux],
            fmt=mymark[i], ls="None", color=mycolor[i], ms=15, markerfacecolor='None', label=label)

    #ymin, ymax = ax1.get_ylim()
    #wid = np.max([abs(ymin), abs(ymax)])
    ax1.set_ylim([0.8, 1.2])
    ax1.grid(which="both")
    
    ax1.legend(fontsize=10, ncol=2)
    if out:
        plt.savefig(out)
    

if __name__ == "__main__":
    parser = ap(
        description="Plot SST spectra.")
    parser.add_argument(
        "--flux", type=str, nargs="*",
        help="SST fluxes")
    parser.add_argument(
        "--labels", type=str, nargs="*",
        help="Labels")
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

    assert len(args.flux) == len(args.labels), "Check input arguments."

    for (flux, label) in zip(args.flux, args.labels):
        df = pd.read_csv(flux, sep=" ")
        epoch_list = sorted(set(df.jd))
        N_epoch = len(epoch_list)
        print(f"{label}")
        print(f"  N_epoch = {N_epoch}")
        # Number of epochs
        
        # Plot fluxes
        out = f"SST_flux_{label}.pdf"
        out = os.path.join(outdir, out)
        plot_SST_flux(df, label, out=out)
        print(f"  Plot SST fluxes and save as {out}.")

        # Make emissivity
        out = f"SST_emissivity_{label}.pdf"
        out = os.path.join(outdir, out)
        plot_SST_emissivity(df, label, out=out)
        print(f"  Plot SST emissivities and save as {out}.")
