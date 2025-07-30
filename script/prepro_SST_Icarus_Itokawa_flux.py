#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Preprocesses of thermal fluxes of Icarus and Itokawa from SST.

../data/Itokawa_SST/
        r17760512/ch0/pbcd/
                  SPITZER_S0_17760512_0002_6_E7593789_tune.tbl
                  SPITZER_S0_17760512_0003_6_E7593791_tune.tbl
                  SPITZER_S0_17760512_0004_6_E7593790_tune.tbl
                  SPITZER_S0_17760512_0005_6_E7593792_tune.tbl
                 ch2/pbcd/
                  SPITZER_S2_17760512_0006_6_E7593261_tune.tbl
                  SPITZER_S2_17760512_0007_6_E7593262_tune.tbl
                  SPITZER_S2_17760512_0008_6_E7593263_tune.tbl
                  SPITZER_S2_17760512_0009_6_E7593264_tune.tbl
        r21563648/ch0/pdcb/
                  SPITZER_S0_21563648_0002_6_E7602125_tune.tbl
                  SPITZER_S0_21563648_0002_6_E8038171_bksub.tbl
                  SPITZER_S0_21563648_0003_6_E7602122_tune.tbl
                  SPITZER_S0_21563648_0003_6_E8038178_bksub.tbl
                  SPITZER_S0_21563648_0004_6_E7602124_tune.tbl
                  SPITZER_S0_21563648_0004_6_E8038183_bksub.tbl
                  SPITZER_S0_21563648_0005_6_E7602126_tune.tbl
                  SPITZER_S0_21563648_0005_6_E8038180_bksub.tbl
                  SPITZER_S0_21563648_0006_6_E7602127_tune.tbl
                  SPITZER_S0_21563648_0006_6_E8038181_bksub.tbl
                  SPITZER_S0_21563648_0007_6_E7602121_tune.tbl
                  SPITZER_S0_21563648_0007_6_E8038182_bksub.tbl
                  SPITZER_S0_21563648_0008_6_E7602128_tune.tbl
                  SPITZER_S0_21563648_0008_6_E8038198_bksub.tbl
                  SPITZER_S0_21563648_0009_6_E7602123_tune.tbl
                  SPITZER_S0_21563648_0009_6_E8038188_bksub.tbl
                ch2/pdcb/
                  SPITZER_S2_21563648_0010_6_E7601644_tune.tbl
                  SPITZER_S2_21563648_0010_6_E8038242_bksub.tbl
                  SPITZER_S2_21563648_0011_6_E7601668_tune.tbl
                  SPITZER_S2_21563648_0011_6_E8038236_bksub.tbl
                  SPITZER_S2_21563648_0012_6_E7601649_tune.tbl
                  SPITZER_S2_21563648_0012_6_E8038247_bksub.tbl
                  SPITZER_S2_21563648_0013_6_E7601652_tune.tbl
                  SPITZER_S2_21563648_0013_6_E8038251_bksub.tbl
                  SPITZER_S2_21563648_0014_6_E7601646_tune.tbl
                  SPITZER_S2_21563648_0014_6_E8038235_bksub.tbl
                  SPITZER_S2_21563648_0015_6_E7601645_tune.tbl
                  SPITZER_S2_21563648_0015_6_E8038238_bksub.tbl
                  SPITZER_S2_21563648_0016_6_E7601661_tune.tbl
                  SPITZER_S2_21563648_0016_6_E8038239_bksub.tbl
                  SPITZER_S2_21563648_0017_6_E7601657_tune.tbl
                  SPITZER_S2_21563648_0017_6_E8038240_bksub.tbl
        r21563904/ch2/pdcb
                  SPITZER_S2_21563904_0010_6_E7601651_tune.tbl
                  SPITZER_S2_21563904_0011_6_E7601656_tune.tbl
                  SPITZER_S2_21563904_0012_6_E7601673_tune.tbl
                  SPITZER_S2_21563904_0013_6_E7601663_tune.tbl
                  SPITZER_S2_21563904_0014_6_E7601671_tune.tbl
                  SPITZER_S2_21563904_0015_6_E7601660_tune.tbl
                  SPITZER_S2_21563904_0016_6_E7601676_tune.tbl
                  SPITZER_S2_21563904_0017_6_E7601675_tune.tbl

../data/Icarus_SST/
        r17760256/ch0/pdcb
                  SPITZER_S0_17760256_0002_7_E7491104_tune.tbl  
                  SPITZER_S0_17760256_0004_7_E7491103_tune.tbl
                  SPITZER_S0_17760256_0003_7_E7491107_tune.tbl  
                  SPITZER_S0_17760256_0005_7_E7491106_tune.tbl
                 ch2/pdcb
                  SPITZER_S2_17760256_0006_7_E7490819_tune.tbl  
                  SPITZER_S2_17760256_0008_7_E7490821_tune.tbl
                  SPITZER_S2_17760256_0007_7_E7490820_tune.tbl  
                  SPITZER_S2_17760256_0009_7_E7490822_tune.tbl
        r17761792/ch0/pdcb/
                  SPITZER_S0_17761792_0000_7_E7483961_tune.tbl
                  SPITZER_S0_17761792_0001_7_E7483952_tune.tbl
                  SPITZER_S0_17761792_0002_7_E7483954_tune.tbl
                  SPITZER_S0_17761792_0003_7_E7483960_tune.tbl
                  SPITZER_S0_17761792_0004_7_E7483962_tune.tbl
                  SPITZER_S0_17761792_0005_7_E7483956_tune.tbl
                  SPITZER_S0_17761792_0006_7_E7483964_tune.tbl
                  SPITZER_S0_17761792_0007_7_E7483951_tune.tbl
                  SPITZER_S0_17761792_0008_7_E7483958_tune.tbl
                  SPITZER_S0_17761792_0009_7_E7483959_tune.tbl
                 ch2/pdcb/
                   SPITZER_S2_17761792_0010_7_E7483668_tune.tbl
                   SPITZER_S2_17761792_0011_7_E7483669_tune.tbl
                   SPITZER_S2_17761792_0012_7_E7483667_tune.tbl
                   SPITZER_S2_17761792_0013_7_E7483665_tune.tbl
                   SPITZER_S2_17761792_0014_7_E7483664_tune.tbl
                   SPITZER_S2_17761792_0015_7_E7483660_tune.tbl
                   SPITZER_S2_17761792_0016_7_E7483663_tune.tbl
                   SPITZER_S2_17761792_0017_7_E7483662_tune.tbl
                   SPITZER_S2_17761792_0018_7_E7483661_tune.tbl
                   SPITZER_S2_17761792_0019_7_E7483657_tune.tbl
                   SPITZER_S2_17761792_0020_7_E7483658_tune.tbl
                   SPITZER_S2_17761792_0021_7_E7483656_tune.tbl
        r17762048/ch0/
                   SPITZER_S0_17762048_0000_7_E7483835_tune.tbl
                   SPITZER_S0_17762048_0001_7_E7483843_tune.tbl
                   SPITZER_S0_17762048_0002_7_E7483846_tune.tbl
                   SPITZER_S0_17762048_0003_7_E7483836_tune.tbl
                   SPITZER_S0_17762048_0004_7_E7483847_tune.tbl
                   SPITZER_S0_17762048_0005_7_E7483850_tune.tbl
                   SPITZER_S0_17762048_0006_7_E7483837_tune.tbl
                   SPITZER_S0_17762048_0007_7_E7483839_tune.tbl
                   SPITZER_S0_17762048_0008_7_E7483834_tune.tbl
                   SPITZER_S0_17762048_0009_7_E7483844_tune.tbl
                  ch2/
                   SPITZER_S2_17762048_0010_7_E7483652_tune.tbl
                   SPITZER_S2_17762048_0011_7_E7483653_tune.tbl
                   SPITZER_S2_17762048_0012_7_E7483650_tune.tbl
                   SPITZER_S2_17762048_0013_7_E7483655_tune.tbl
                   SPITZER_S2_17762048_0014_7_E7483649_tune.tbl
                   SPITZER_S2_17762048_0015_7_E7483648_tune.tbl
                   SPITZER_S2_17762048_0016_7_E7483651_tune.tbl
                   SPITZER_S2_17762048_0017_7_E7483654_tune.tbl
                   SPITZER_S2_17762048_0018_7_E7483647_tune.tbl
                   SPITZER_S2_17762048_0019_7_E7483646_tune.tbl
                   SPITZER_S2_17762048_0020_7_E7483643_tune.tbl
                   SPITZER_S2_17762048_0021_7_E7483642_tune.tbl
        r17762304/ch0/
                   SPITZER_S0_17762304_0000_7_E7483852_tune.tbl
                   SPITZER_S0_17762304_0001_7_E7483853_tune.tbl
                   SPITZER_S0_17762304_0002_7_E7483848_tune.tbl
                   SPITZER_S0_17762304_0003_7_E7483856_tune.tbl
                   SPITZER_S0_17762304_0004_7_E7483857_tune.tbl
                   SPITZER_S0_17762304_0005_7_E7483845_tune.tbl
                   SPITZER_S0_17762304_0006_7_E7483849_tune.tbl
                   SPITZER_S0_17762304_0007_7_E7483854_tune.tbl
                   SPITZER_S0_17762304_0008_7_E7483851_tune.tbl
                   SPITZER_S0_17762304_0009_7_E7483855_tune.tbl
                  ch2/
                   SPITZER_S2_17762304_0010_7_E7483640_tune.tbl
                   SPITZER_S2_17762304_0011_7_E7483641_tune.tbl
                   SPITZER_S2_17762304_0012_7_E7483638_tune.tbl
                   SPITZER_S2_17762304_0013_7_E7483637_tune.tbl
                   SPITZER_S2_17762304_0014_7_E7483635_tune.tbl
                   SPITZER_S2_17762304_0015_7_E7483634_tune.tbl
                   SPITZER_S2_17762304_0016_7_E7483633_tune.tbl
                   SPITZER_S2_17762304_0017_7_E7483631_tune.tbl
                   SPITZER_S2_17762304_0018_7_E7483629_tune.tbl
                   SPITZER_S2_17762304_0019_7_E7483628_tune.tbl
                   SPITZER_S2_17762304_0020_7_E7483627_tune.tbl
                   SPITZER_S2_17762304_0021_7_E7483626_tune.tbl
        r4871168/ch0
                   SPITZER_S0_4871168_0002_8_E7362453_tune.tbl
                   SPITZER_S0_4871168_0003_8_E7362454_tune.tbl
                   SPITZER_S0_4871168_0004_8_E7362455_tune.tbl
                   SPITZER_S0_4871168_0005_8_E7362456_tune.tbl
"""
import os 
from argparse import ArgumentParser as ap
import pandas as pd
import numpy as np
from astroquery.jplhorizons import Horizons
from astropy.constants import c, au
from astropy.time import Time

from Eros_common import remove_largevar, remove_edge, SST_ltcor


def extract_DATE_OBS_SST(f0):
    with open(f0, 'r') as f:
        for line in f:
            if 'DATE_OBS' in line:
                start = line.find("'") + 1
                end = line.find("'", start)
                t_utc = line[start:end]
                # Convert to jd
                t_jd = Time(str(t_utc), format='isot', scale='utc').jd
                return t_jd
    return None  


def extract_SST_flux(fdir, out):
    """Extract SST flux and save them.

    Parameters
    ----------
    fdir : str
        Downloaded directory contains SST spectra.
    out : str
        output filename
    """

    # Not used
    key_skip = [
        "\\processing", "\\wavsamp", "\\char", "\\character", 
        "\\COMMENT",  "\\HISTORY", "\\int", "\\float"]
    df_list = []
    
    # The directory structure is as follows:
    #    fdir/
    #      README
    #      r17760256/ch0/pbcd
    #                ch2/pbcd
    #      r17761792/ch0/pbcd
    #                ch2/pbcd
    blocks = [name for name in os.listdir(fdir)
          if os.path.isdir(os.path.join(fdir, name)) and name.startswith('r')]
    blocks = sorted(blocks)
    N_block = len(blocks)
    print(f"Extract SST fluxes from {fdir}")
    print(f"  Number of blocks (starting 'r', e.g., r17760256) N={N_block}")

    for block in blocks:
        print(f"    Block {block}")
        block = os.path.join(fdir, block)
        channels = [name for name in os.listdir(block)
          if os.path.isdir(os.path.join(block, name)) and name.startswith('ch')]
        N_ch = len(channels)
        print(f"      Number of channels (starting 'ch', e.g., ch0) N={N_ch}")

        for channel in channels:
            print(f"        {channel}")
            channel = os.path.join(block, channel, "pbcd")
            tbls = [name for name in os.listdir(channel)
              if os.path.isfile(os.path.join(channel, name)) and name.startswith('S')]
            N_tbl = len(tbls)
            print(f"          Number of files (starting 'SPITZER') N={N_tbl}")

            for f in tbls:
                f0 = os.path.join(channel, f)
                # Extract DATE_OBS independently in jd
                t_jd = extract_DATE_OBS_SST(f0)
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
                df_1spec["jd"] = t_jd
                # Remove bad data with but-flag
                df_1spec = df_1spec[df_1spec["bit-flag"] == "0"]
                df_1spec = df_1spec.astype(float)
                df_list.append(df_1spec)

    df = pd.concat(df_list)
    # Rename
    df = df.rename(columns={"flux_density":"flux", "error":"fluxerr"})
    # Save
    df.to_csv(out, sep=" ", index=False)


def make_ave_SST(df1, df2):
    """Make avareged spectrum for SST (of Eros).

    Parameters
    ----------
    df1 : pandas.DataFrame
        input data1
    df2 : pandas.DataFrame
        input data2

    Return
    ------
    df12 : pandas.DataFrame
        output merged data
    """

    # Most of the wavelengths are the same, but not always all of them.
    w1 = sorted(set(df1.wavelength))
    w2 = sorted(set(df2.wavelength))

    # Final wavelength grids
    w12 = sorted(set(w1 + w2))

    print(f"    Wavelength1  N = {len(w1)}")
    print(f"    Wavelength2  N = {len(w2)}")
    print(f"    -> combined  N = {len(w12)}\n")

    # Set wavelength as index
    df1_ = df1.set_index("wavelength")
    df2_ = df2.set_index("wavelength")

    # Make a new dataframe
    f_mean = []
    ferr_mean = []

    for w in w12:
        f1 = df1_["flux"][w] if w in df1_.index else np.nan
        f2 = df2_["flux"][w] if w in df2_.index else np.nan
        e1 = df1_["fluxerr"][w] if w in df1_.index else np.nan
        e2 = df2_["fluxerr"][w] if w in df2_.index else np.nan

        # Mean
        if not np.isnan(f1) and not np.isnan(f2):
            f = 0.5 * (f1 + f2)
            # Error
            ferr_ave = 0.5 * np.abs(f1 - f2) 
            # Original error
            ferr_ori = (0.5*(e1**2 + e2**2))**0.5
            ferr = (ferr_ave**2 + ferr_ori**2)**0.5
        elif not np.isnan(f1):
            f = f1
            ferr = e1
        elif not np.isnan(f2):
            f = f2
            ferr = e2
        else:
            f = np.nan
            ferr = np.nan
            assert False, "Check the code."

        f_mean.append(f)
        ferr_mean.append(ferr)

    df_12 = pd.DataFrame({
        "wavelength": w12,
        "flux": f_mean,
        "fluxerr": ferr_mean
    })
    return df_12



if __name__ == "__main__":
    parser = ap(
        description="Preprocess of thermal fluxes of Icarus and Itokawa.")
    parser.add_argument(
        "fdir_Icarus", type=str, 
        help="input file directory for Icarus")
    parser.add_argument(
        "fdir_Itokawa", type=str, 
        help="input file directory for Itokawa")
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
    
    out = "SST_flux_Icarus.txt"
    out = os.path.join(args.outdir, out)
    extract_SST_flux(args.fdir_Icarus, out)

    out = "SST_flux_Itokawa.txt"
    out = os.path.join(args.outdir, out)
    extract_SST_flux(args.fdir_Itokawa, out)

    
    assert False, 1
    # # Do light-time correction, save results, and plot
    # Obs time
    # DATE_OBS: starting time
    #   2004-09-30T00:52:07.827
    # SAMPTIME: integtime?
    #   1.0486 s
    # EXPTOT_T: w/deadtime?
    #   14.68 s

    # Relative error
    err_th = 0.05
    # Variation (25%)
    var_th = 0.25

    key_w, key_flux, key_fluxerr = "wavelength", "flux_density", "error"
    df_list_cor = []
    for idx, df in enumerate(df_list):

        # Remove large error
        df = df[df[key_fluxerr]/df[key_flux] < err_th]
        # Remove large variation
        df = remove_largevar(df, key_flux, var_th)
        # Remove edge
        df = remove_edge(df, 3)
        df_list_cor.append(df)
