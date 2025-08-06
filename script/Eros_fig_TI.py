#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Preprocesses of thermal inertia from previous studies.

Hung+2022, PSJ, Table 3
https://content.cld.iop.org/journals/2632-3338/3/3/56/revision2/psjac4d1ft3_mrt.txt

Novakovic+2024, PSJ, Table 11
https://iopscience.iop.org/2632-3338/5/1/11/suppdata/psjad08c0t11_ascii.txt?doi=10.3847/PSJ/ad08c0
"""
import os 
from argparse import ArgumentParser as ap
import pandas as pd
import numpy as np
from astropy.constants import c, au
from astropy.time import Time
import matplotlib.pyplot as plt
import re

from Eros_common import mycolor, mymark


def _parse_float(s):
    try:
        return float(s.strip())
    except:
        return None


def _parse_int(s):
    try:
        return int(s.strip())
    except:
        return None


def read_H22(f):
    """Read Table 11 of Hung+2022, PSJ.

    Parameter
    ---------
    f : str
        downloaded text table

    Return
    ------
    df : pandas.DataFrame
        dataframe with information
    """
    data = []

    with open(f, encoding='utf-8') as file:
        lines = file.readlines()[44:]

        for line in lines:
            if not line.strip():
                continue 
            try:
                row = {
                    "obj": line[0:21].strip().replace(" ", "_").replace("*", ""),
                    "DWISE": _parse_float(line[22:27]),
                    "e_DWISE": _parse_float(line[28:32]),
                    "D_km": _parse_float(line[33:38]),
                    "Derr_km_u": _parse_float(line[39:43]),
                    "Derr_km_l": _parse_float(line[44:48]),
                    "TI": _parse_int(line[49:53]),
                    "TIerr_u": _parse_int(line[54:58]),
                    "TIerr_l": _parse_int(line[59:63]),
                    "A": _parse_float(line[64:69]),
                    "Aerr_u": _parse_float(line[70:75]),
                    "Aerr_l": _parse_float(line[76:81]),
                    "pV": _parse_float(line[82:87]),
                    "pVerr_u": _parse_float(line[88:93]),
                    "pVerr_l": _parse_float(line[94:99]),
                    "theta": _parse_float(line[100:104]),
                    "chi2": _parse_float(line[105:109]),
                    "Model": line[110:112].strip()
                }
                data.append(row)
            except Exception as e:
                print(f"Error parsing line: {line}")
                print(e)

    df = pd.DataFrame(data)
    df = df.fillna(-999)
    return df


def read_N24(f):
    """Read Table 3 of Novakovic+2024, PSJ.

    Parameter
    ---------
    f : str
        downloaded text table

    Return
    ------
    df : pandas.DataFrame
        dataframe with information
    """

    with open(f, encoding='utf-8') as file:
        lines = file.readlines()

    # Header
    header_index = None
    for i, line in enumerate(lines):
        if line.strip().startswith("Asteroid"):
            header_index = i
            break
    if header_index is None:
        raise ValueError("Header line starting with 'Asteroid' not found.")
    
    # Data
    data_lines = []
    for line in lines[header_index + 2:]:
        line = line.strip()
        if not line or line.startswith("Note."):
            break
        data_lines.append(line)

    records = []

    for line in data_lines:
        parts = re.split(r'\t+|\s{2,}', line.strip())
        while len(parts) < 9:
            parts.append(None)

        asteroid_raw = parts[0]
        asteroid = asteroid_raw.replace("(", "").replace(")", "").replace(" ", "_")

        # H
        try:
            H = float(parts[1])
        except:
            H = None

        # D, Derr
        D_km, Derr_km = None, None
        D_match = re.match(r"([\d.]+)\s*\+or-\s*([\d.]+)", parts[2] or "")
        if D_match:
            D_km = float(D_match.group(1))
            Derr_km = float(D_match.group(2))

        # pV, pVerr
        pV, pVerr = None, None
        pV_match = re.match(r"([\d.]+)\s*\+or-\s*([\d.]+)", parts[3] or "")
        if pV_match:
            pV = float(pV_match.group(1))
            pVerr = float(pV_match.group(2))

        # dm, P
        try:
            Delta_m = float(parts[4])
        except:
            Delta_m = None
        try:
            P = float(parts[5])
        except:
            P = None

        # TI_nominal
        TI_nominal, TIerr_l, TIerr_u = None, None, None
        gamma_nom_match = re.match(
            r"\$\{(\d+(?:\.\d+)?)\}_\{\-(\d+(?:\.\d+)?)\}\^\{\+(\d+(?:\.\d+)?)\}\$", parts[6] or "")
        if gamma_nom_match:
            TI_nominal = float(gamma_nom_match.group(1))
            TIerr_l = float(gamma_nom_match.group(2))
            TIerr_u = float(gamma_nom_match.group(3))
        
        # TI_alternative
        TI_alt, TIalt_err_l, TIalt_err_u = None, None, None
        gamma_alt_match = re.match(
            r"\$\{(\d+(?:\.\d+)?)\}_\{\-(\d+(?:\.\d+)?)\}\^\{\+(\d+(?:\.\d+)?)\}\$", parts[7] or "")
        if gamma_alt_match:
            TI_alt = float(gamma_alt_match.group(1))
            TIalt_err_l = float(gamma_alt_match.group(2))
            TIalt_err_u = float(gamma_alt_match.group(3))

        records.append({
            "obj": asteroid,
            "H": H,
            "D_km": D_km,
            "Derr_km": Derr_km,
            "pV": pV,
            "pVerr": pVerr,
            "Delta_m": Delta_m,
            "P": P,
            "TI": TI_nominal,
            "TIerr_u": TIerr_u,
            "TIerr_l": TIerr_l,
            "TI_alt": TI_alt,
            "TIalt_err_u": TIalt_err_u,
            "TIalt_err_l": TIalt_err_l
        })

    df = pd.DataFrame.from_records(records)
    df = df.fillna(-999) 
    return df



if __name__ == "__main__":
    parser = ap(
        description="Preprocess of thermal fluxes of Icarus and Itokawa.")
    parser.add_argument(
        "H22", type=str,
        help="Table 11 of Hung+2022")
    parser.add_argument(
        "N24", type=str,
        help="Table 3 of Novakovic+2024")
    parser.add_argument(
        "--out", type=str, default="TI.pdf",
        help="Output filename")
    args = parser.parse_args()
    
    col_use = ["obj", "TI", "D_km"]
    col_H22, col_N24 = mycolor[0], mycolor[1]
    mark_H22, mark_N24 = mymark[0], mymark[1]

    # Read Hung+2022, Table 11
    df_H22 = read_H22(args.H22)
    df_H22 = df_H22[col_use]

    # Read Novakovic+2024, Table 3
    df_N24 = read_N24(args.N24)
    df_N24 = df_N24[col_use]

    # Plot D km vs. TI
    label_H22 = f"Hung+2022\n  N={len(df_H22)}"
    label_N24 = f"Novakovic+2024\n  N={len(df_N24)}"

    fig = plt.figure(figsize=(12, 12))
    ax1 = fig.add_axes([0.10, 0.10, 0.85, 0.85])
    ax1.set_xlabel("Diameter [km]")
    ax1.set_ylabel("Thermal inertia [tiu]")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.scatter(
        df_H22["D_km"], df_H22["TI"], color=col_H22, marker=mark_H22, 
        label=label_H22)
    ax1.scatter(
        df_N24["D_km"], df_N24["TI"], color=col_N24, marker=mark_N24, 
        label=label_N24)
    # Eros
    ax1.scatter(
        13, 60, color="white", marker="*", facecolor="orange", s=1500,
        label="Eros (This study)")
    
    ax1.legend(loc="lower left", fontsize=20)
    plt.savefig(args.out)
