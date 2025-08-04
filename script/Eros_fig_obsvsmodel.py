#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plot observations vs. model fluxes.

 Harris: 7 (spec)
   2450991.767627034, 2450991.809986014, 2450993.7376662632, 2450993.7432215353, 2450994.7403925494, 2450994.7924731895, 2450994.8341376856,
 Lim   : 1 (spec)
   2452539.768144042,
 Wolter: 1 (spec)
   2452545.9060276826
 SST   : 4 (spec)
   2453278.5354158804, 2453278.535753845, 2453278.5361288944, 2453278.5364668127
 AKARI : 5 (phot)
   2454199.2591583026, 2454199.397187427, 2454199.8112584595, 2454199.880273481, 2454199.949288874
"""
import os 
from argparse import ArgumentParser as ap
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from Eros_common import mycolor, tel2col, tel2lab, plot_spec_Eros, generate_latex_table

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd

def extract_flux(f0, fixscale=False):
    """
    Extract thermal flux from output of TPM.

    Parameters
    ----------
    f : str
        output file of TPM
    fixscale : bool, optional
        whether fix the scale factor (i.e., trust shape model)

    Return
    ------
    df : pandas.DataFrame
        dataframe with extracted fluxes etc.
    """

    epoch_list, jd_list, w_list = [], [], []
    f_obs_list, ferr_obs_list, f_model_list = [], [], []
    with open (f0, "r") as f:
        f = f.readlines()
        for l in f:
            # l is as below:
            # f> 052     2454199.2663881136   000 18.000        1.698843   0.113268    1.441020    0.000000      1.1789
            #    epoch   JD                   n   wavelength    f_obs      ferr_obs    f_model     ferr_model    f_obs/f_model
            if l[0:2] == "f>":
                # Extract info
                l = l.split()
                epoch, jd, n       = l[1], l[2], l[3]
                w, f_obs, ferr_obs = l[4], l[5], l[6]
                f_model            = l[7]
                epoch_list.append(float(epoch))
                jd_list.append(float(jd))
                w_list.append(float(w))
                f_obs_list.append(float(f_obs))
                ferr_obs_list.append(float(ferr_obs))
                f_model_list.append(float(f_model))
            elif l[0:2] == "r>":
                # Extract scale factor
                # r>     100.0   0.0   0.0  0.00  0.12 1.15125549     19.364     21.576  1.000   2259.587    4.18547
                l = l.split()
                # Fix scale factor to 1
                if fixscale:
                    scalefactor = 1
                # Use scale factor in TPM output
                else:
                    scalefactor = float(l[6])

                # This is useful to check intputs of the TPM calculation
                # Haple angle in deg (t.Bar in the output)
                TI = float(l[1])
                Htheta = float(l[2])
                A = float(l[5])

            else:
                continue

        # DataFrame
        df = pd.DataFrame({
            "epoch": epoch_list,
            "jd": jd_list,
            "w": w_list,
            "f_obs": f_obs_list,
            "ferr_obs": ferr_obs_list,
            "f_model": f_model_list,
            })
        df["scalefactor"] = scalefactor
        df["TI"] = TI
        df["Htheta"] = Htheta
        df["A"] = A
        #print(f"Scafe factor = {scalefactor}")
    return df



def plot_all_epochs(df, epoch_list, Nparam, out=None):
    """
    Plots observed/model flux and standardized residuals for all epochs in 5x3 panel layout.
    Also computes and displays chi-squared per epoch and total.
    """

    # Degree of freedom
    Ndata = len(df)
    nu = Ndata - Nparam

    col_obs = "black"
    col_model = "red"
    num_epochs = len(epoch_list)
    rows, cols = 5, 4

    # --- STEP 1: compute max standardized residual for ylim ---
    max_residual = 0
    chi2_list = []
    for epoch in epoch_list:
        df_e = df[df["jd"] == epoch]
        res = (df_e["f_obs"] - df_e["f_model"]) / df_e["ferr_obs"]
        if not res.empty:
            max_residual = max(max_residual, np.abs(res).max())
            chi2 = np.sum(res**2)
            chi2_list.append(chi2)
        else:
            chi2_list.append(0.0)

    ylim_resid = np.ceil(max_residual)
    total_chi2 = np.sum(chi2_list)

    # --- STEP 2: plot ---
    fig = plt.figure(figsize=(15, 16))
    gs = gridspec.GridSpec(rows * 2, cols, figure=fig, hspace=0.5, wspace=0.3)

    for i in range(rows * cols):
        row = (i // cols) * 2
        col = i % cols

        if i < num_epochs:
            epoch = epoch_list[i]
            df_e = df[df["jd"] == epoch]

            ax_main = fig.add_subplot(gs[row, col])
            ax_resid = fig.add_subplot(gs[row + 1, col], sharex=ax_main)

            # Main plot
            ax_main.errorbar(df_e["w"], df_e["f_obs"], df_e["ferr_obs"],
                             fmt="o", ms=5, color=col_obs, label=f"Obs ({epoch})")
            ax_main.scatter(df_e["w"], df_e["f_model"], color=col_model,
                            s=20, facecolor="None", marker="o", label="Model")
            ax_main.set_ylabel("Flux [Jy]", fontsize=10)
            ax_main.tick_params(labelbottom=False)
            ax_main.legend(fontsize=8).set_alpha(1)

            # Chi-squared for this epoch
            chi2 = chi2_list[i]

            # Residuals
            residual = (df_e["f_obs"] - df_e["f_model"]) / df_e["ferr_obs"]
            ax_resid.axhline(0, color='gray', linestyle='--', linewidth=1)
            ax_resid.errorbar(df_e["w"], residual, fmt="o", color="gray", ms=3, label=f"χ²={chi2:.1f}")
            ax_resid.set_xlabel("Wavelength", fontsize=10)
            ax_resid.set_ylabel("Residual", fontsize=9)
            ax_resid.set_ylim(-ylim_resid, ylim_resid)
            ax_resid.legend().set_alpha(1)
        else:
            # 空きパネルに合計chi2表示
            ax_main = fig.add_subplot(gs[row, col])
            ax_main.axis("off")
            if i == num_epochs:
                ax_main.text(0.5, 0.5, f"Total χ² = {total_chi2:.2f} (reduced χ² = {total_chi2/nu:.2f})",
                             fontsize=14, ha="center", va="center", transform=ax_main.transAxes)

            fig.add_subplot(gs[row + 1, col]).axis("off")

    plt.tight_layout()
    if out:
        plt.savefig(out)
        plt.close()


def introduce_sf(df, epochs_spec, df_sf):
    """Introduce scale factors.

    Parameters
    ----------
    df : pandas.DataFrame
        input data with f_model
    epochs_spec : array-like
        epochs of spectroscopy in jd
    df_sf : pandas.DataFrame
        input data with scale factors (sf1, sf2, ...)
    
    Return
    ------
    df : pandas.DataFrame
        dataframe with scale factor introduced
    """
    for idx, e in enumerate(epochs_spec):
        # Extract corresponding scale factor
        sf = df_sf[f"sf{idx+1}"].values[0]
        df.loc[df["jd"]==e, "f_model"] *= sf**2
    return df

if __name__ == "__main__":
    parser = ap(
        description="Plot obs vs. model fluxes.")
    parser.add_argument(
        "res", type=str, 
        help="Results of tpm")
    parser.add_argument(
        "--sf", type=str, default=False,
        help="Scale factors (only for w/NN)")
    parser.add_argument(
        "--NN", action="store_true", default=False,
        help="Whether input file is NN predection")
    parser.add_argument(
        "--Nparam", default=2,
        help="Number of parameters")
    parser.add_argument(
        "--out", type=str, default=None,
        help="Output figure name")
    args = parser.parse_args()

    # Results of TPM
    f_res = args.res
    f_sf = args.sf
    # Number of parameters
    Nparam = args.Nparam
    # Output file
    out = args.out

    if args.NN:
        df = pd.read_csv(f_res, sep=" ")
        # Best fit thermal inertia and thetabar
        TI = args.TI
        Htheta = args.Htheta

        if args.sf:
            # Chi2 and scale factors
            df_sf = pd.read_csv(f_sf, sep=" ")
            # N_spec = N_scalefactor = 13
            epochs = sorted(list(set(df.jd)))
            # Remove akari
            #epochs = [x for x in epochs if x < 2454199]
            # Introuce scale factors for spec (N=13)
            epochs_spec = [x for x in epochs if x < 2454199]
            
            df_selected = df[(df["TI"].round(2) == round(TI, 2)) & (df["Htheta"].round(2) == round(Htheta, 2))]
            # Extract scale factors for a combination of TI and Htheta
            df_sf_selected = df_sf[(df_sf["TI"] == TI) & (df_sf["Htheta"] == Htheta)]
            df_sf_introduced = introduce_sf(df_selected, epochs_spec, df_sf_selected)
            plot_all_epochs(df_sf_introduced, epochs, Nparam, out)

        else:
            df_selected = df[(df["TI"].round(2) == round(TI, 2)) & (df["Htheta"].round(2) == round(Htheta, 2))]
            # TI, Htheta
            Nparam = 2
            
            # Without scalefactors
            plot_all_epochs(df_selected, epochs, Nparam, out)


    else:
        df = extract_flux(f_res)
        epochs = sorted(list(set(df.jd)))
        plot_all_epochs(df, epochs, Nparam, out)
