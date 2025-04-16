#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Plot FoV of SST at the time of observations of Eros.

The observations of Eros were performed with IRS Map mode.
"""
import os
from argparse import ArgumentParser as ap
from astroquery.jplhorizons import Horizons
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
from matplotlib.ticker import FormatStrFormatter
import matplotlib.patches as patches
from astropy.coordinates import Angle
from matplotlib.transforms import Affine2D

from myplot import mymark, mycolor, mygrid


def add_rectangle(ax, ra, dec, fltr, pa=0):
    if fltr == "SL":
        width_arcsec = 3.6
        height_arcsec = 57.0
        color = "blue"
    elif fltr == "LL":
        width_arcsec = 10.5
        height_arcsec = 57.0
        color = "green"

    width_deg = width_arcsec / 3600.0 / np.cos(np.radians(dec))
    height_deg = height_arcsec / 3600.0

    ra0 = ra - width_deg / 2
    dec0 = dec - height_deg / 2

    rect = patches.Rectangle(
        (ra0, dec0),
        width_deg,
        height_deg,
        linewidth=0.8,
        edgecolor=color,
        facecolor="none"
    )
    # Rotate
    transform = (
        Affine2D().rotate_deg_around(ra, dec, pa) + ax.transData)
    rect.set_transform(transform)
    ax.add_patch(rect)


def plot_SSTFoV(ast, df, t0, t1, tstep="1m", out="out.png"):
    """Plot location of asteroid with FoV of SST.

    Parameters
    ----------
    ast : str 
        asteroid name to query to JPL Horizons
    df : pandas.DataFrame
        dataframe with jd (julian day) etc.
    t0 : str
        starting time
    t1 : str
        ending time
    tstep : str
        time step
    out : str, optional
        file name if you make a plot
    """
    # Plot orbit
    res = Horizons(
        location='@sst', id=ast, epochs={'start':t0, 'stop':t1, 'step':tstep})
    eph = res.ephemerides()

    fig = plt.figure(figsize=(16, 6))
    ax = fig.add_axes([0.10, 0.05, 0.85, 0.9])
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))
    # Not equal
    cos_dec = np.cos(np.radians(eph["DEC"][0]))
    ax.set_aspect(1 / cos_dec)
    ax.set_xlabel("RA [deg]")
    ax.set_ylabel("Dec [deg]")
    
    ax.plot(eph["RA"], eph["DEC"], color="black", lw=1.5, ls="solid")
   
    # Plot Location at the time of observations
    for i, row in df.iterrows():
        date = row["date"]
        onslit = row["onslit"]
        fltr = row["fltr"]
        # Convert to jd
        jd = Time(date, format="fits").jd
        res = Horizons(
            location='@sst', id=ast, epochs=jd)
        eph = res.ephemerides()
        # Real position
        ra_obj, dec_obj = eph["RA"][0], eph["DEC"][0]
        ax.scatter(ra_obj, dec_obj, color="black", marker="x")
        
        # Pointing 
        ra, dec, pa = row["RA_RQST"], row["DEC_RQST"], row["PA_RQST"]
        # Important
        pa = -pa
        if onslit == 0:
            color = "gray"
        elif onslit == 1:
            color = "red"
        ax.scatter(ra, dec, color=color)

        add_rectangle(ax, ra, dec, fltr, pa=pa)
    ax.legend()
    plt.savefig(out)
    plt.close()


if __name__ == "__main__":
    parser = ap(
        description="Plot FoV of SST.")
    parser.add_argument(
        "sst", type=str,
        help="SST obs info.")
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
    df = pd.read_csv(args.sst, sep=" ")

    out = f"Eros_fig_SSTFoV.{args.outtype}"
    out = os.path.join(args.outdir, out)
    # 2004-09-30 00:48:16 to 2004-09-30 00:52:30.732617
    t0 = "2004-09-30 00:51"
    t1 = "2004-09-30 01:01"
    plot_SSTFoV("433", df, t0, t1, out=out)
