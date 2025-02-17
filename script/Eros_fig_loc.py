#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plot location of Eros and Earth.
"""
import os
from argparse import ArgumentParser as ap
from astroquery.jplhorizons import Horizons
import pandas as pd
import matplotlib.pyplot as plt
from astropy.time import Time
from scipy import interpolate
from astropy.constants import c, au

from myplot import mymark, mycolor, mygrid, mymathtext


def plot_ast(asteroid, df, date0="1974-01-01", date1="2024-12-31", wid=1.8, out=None):
    """
    Plot location of asteroid and the Earth.

    Parameters
    ----------
    asteroid : str 
        asteroid name to query to JPL Horizons
    df : pandas.DataFrame
        dataframe with jd (julian day), col (color), and lab (label)
    date0, date1 : str
        starting and ending date
    wid : float
        width of the plot in au
    out : str, optional
        file name if you make a plot
    """
    # 500 represents the center
    # @0 represents that coordinate origin is the solar barycenter, 
    # @10 represents that coordinate origin is the Sun-body center
    # location='500@0' seems bad for this purpose
    # The location is not overlapped.
    # Default is 500@10 (Sun-body center), not 500@0 (barucenter)

    
    # Just to plot circle
    ast1 = Horizons(
        location='500@10',
        id=asteroid, epochs={'start':"2024-01-01", 'stop':"2025-02-10", 'step':"1d"})
    ast1 = ast1.vectors()
    # Earth
    Earth1 = Horizons(
        location='500@10',
        id="399", epochs={'start':"2024-01-01", 'stop':"2025-01-01", 'step':"1d"})
    Earth1 = Earth1.vectors()
    
    
    # To plot the exact location at time time of observations
    ast = Horizons(
        location='500@10',
        id=asteroid, epochs={'start':date0, 'stop':date1, 'step':"1d"})

    ast = ast.vectors()
    
    # Earth
    Earth = Horizons(
        location='500@10',
        id="399", epochs={'start':date0, 'stop':date1, 'step':"1d"})
    Earth = Earth.vectors()

    # Prediction
    f_ast_x = interpolate.interp1d(ast["datetime_jd"], ast["x"], kind='linear', fill_value='extrapolate')
    f_ast_y = interpolate.interp1d(ast["datetime_jd"], ast["y"], kind='linear', fill_value='extrapolate')
  
    
    f_E_x = interpolate.interp1d(Earth["datetime_jd"], Earth["x"], kind='linear', fill_value='extrapolate')
    f_E_y = interpolate.interp1d(Earth["datetime_jd"], Earth["y"], kind='linear', fill_value='extrapolate')


    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_axes([0.15, 0.15, 0.8, 0.8])
    mymathtext()
    mygrid(ax)
    ax.set_xlim([-wid, wid])
    ax.set_ylim([-wid, wid])
    ax.set_xlabel("x [au]")
    ax.set_ylabel("y [au]")
    
    # asteroid
    lab = asteroid
    lab = None
    ax.plot(ast["x"], ast["y"], color="black", lw=1.5, ls="solid", label=lab)
    # Earth
    lab = "Earth"
    lab = None
    ax.plot(Earth1["x"], Earth1["y"], color="grey", lw=1.5, ls="dashed", label=lab)
    # Sun
    ax.scatter(0, 0, marker="x", color="black", lw=1.5, s=200)

    for i, row in df.iterrows():
        jd = row["jd"]
        col = row["col"]
        lab = row["lab"]
        # Predict location
        x_ast, y_ast = f_ast_x(jd), f_ast_y(jd)
        x_E, y_E = f_E_x(jd), f_E_y(jd)

        t = Time(str(jd), format='jd', scale='utc')
        year = t.datetime.year
        mark = mymark[i]
        ax.scatter(
            x_ast, y_ast,
            color=col, marker=mark, s=250, lw=1,
            label=lab, zorder=10, ec="black",)
        ax.scatter(
            x_E, y_E,
            color=col, marker=mark, s=250, lw=1,
            label=None, zorder=10, ec="black")

    ax.legend(loc="upper left")
    if out:
        plt.savefig(out)
    plt.close()


if __name__ == "__main__":
    parser = ap(
        description="Plot location of Eros.")
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

    # 2025-02-11 ver. 
    # Harris1999(UKIRT98) : 1998-06-27, 1998-06-28, 1998-06-29, 1998-06-30
    # Lim+2005            : (2002-09-21, 2002-09-22, )2002-09-23
    # Wolters2008         : 2002-09-28
    # SST                 : 2004-09-30
    # AKARI               : 2007-04-08
    utc_list = [
        "1998-06-27", 
        "2002-09-23", 
        "2002-09-28", 
        "2004-09-30",
        "2007-04-08",
        ]
    col_list = [
        mycolor[0],
        mycolor[1],
        mycolor[2],
        mycolor[3],
        mycolor[5]
        ]
    lab_list = [ 
        "1998 Jun. (Harris+1999)",
        "2002 Sep. (Wolters+2009)",        
        "2002 Sep. (Lim+2005)",
        "2004 Sep. (SST)",
        "2007 Apr. (AKARI)",
        ]
    jd_list = [Time(x, format='isot', scale='utc').jd for x in utc_list]
    df_E = pd.DataFrame(dict(jd=jd_list, col=col_list, lab=lab_list))

    out = f"Eros_fig_loc.{args.outtype}"
    out = os.path.join(args.outdir, out)
    plot_ast("433", df_E, out=out)
