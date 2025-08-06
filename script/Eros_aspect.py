#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Check aspect data of Eros.
"""
import os 
from argparse import ArgumentParser as ap
import pandas as pd
import numpy as np
from astroquery.jplhorizons import Horizons
from astropy.constants import c, au


if __name__ == "__main__":
    parser = ap(
        description="Check aspect data of Eros.")
    args = parser.parse_args()

    # Harris: 7 (spec)
    jd_H = [2450991.767627034, 2450991.809986014, 2450993.7376662632, 2450993.7432215353, 2450994.7403925494, 2450994.7924731895, 2450994.8341376856]
    # Lim   : 1 (spec)
    jd_L = [2452539.768144042]
    # Wolter: 1 (spec)
    jd_W = [2452545.9060276826]
    # SST   : 2 (merged)
    jd_S = [2453278.535584863, 2453278.5362978536]
    #   (from 2453278.5354158804, 2453278.535753845, 2453278.5361288944, 2453278.5364668127)
    # AKARI : 5 (phot)
    jd_A = [2454199.2591583026, 2454199.397187427, 2454199.8112584595, 2454199.880273481, 2454199.949288874]

    # Checked with JPL HORIZONS ORBITER by eyes.
    # -, -, -, +, +

    for idx, jd_list in enumerate([jd_H, jd_L, jd_W, jd_S, jd_A]):
        print(f"Dataset {idx+1}")

        r_list = []
        alpha_list = []
        for jd in jd_list:
            # Use code 500
            ast = Horizons(location='500', id=433, epochs=jd)
            eph = ast.ephemerides()
            # Heliocentric distance 
            r = eph["r"][0]
            # Phase angle
            alpha = eph["alpha"][0]
            # Solar elongation
            elong = eph['elong'][0]

            print(f"  Epoch = {jd}")
            print(f"    Heliocentric dist. {r:.3f} au")
            print(f"    Phase angle        {alpha:.3f} deg")


            # Get signed phase angle ==========================================
            ## Observer-target
            #vec_obs = ast.vectors()
            ## Sun-target
            #obj_sun = Horizons(
            #    id=433, location='500@10', epochs=jd)
            #vec_sun = obj_sun.vectors()
            #

            ## ターゲットから見た Observer と Sun の位置ベクトル
            ## observer_vec = - (target wrt observer) なのでそのまま
            ## sun_vec = - (target wrt sun)
            #
            #pos_obs = np.vstack([
            #    np.asarray(vec_obs['x']).data,
            #    np.asarray(vec_obs['y']).data,
            #    np.asarray(vec_obs['z']).data
            #]).T
            #
            #pos_sun = np.vstack([
            #    np.asarray(vec_sun['x']).data,
            #    np.asarray(vec_sun['y']).data,
            #    np.asarray(vec_sun['z']).data
            #]).T

            ## pos_obs, pos_sun: shape = (N, 3)
            #
            ## 1. XY平面に投影
            #vec_obs_xy = pos_obs[:, :2]
            #vec_sun_xy = pos_sun[:, :2]
            #
            ## 2. 単位ベクトル
            #u_obs = vec_obs_xy / np.linalg.norm(vec_obs_xy, axis=1)[:, None]
            #u_sun = vec_sun_xy / np.linalg.norm(vec_sun_xy, axis=1)[:, None]


            #
            ## 3. 位相角（0〜180度）
            #dot = np.sum(u_obs * u_sun, axis=1)
            #dot = np.clip(dot, -1, 1)
            #phase_angle = np.arccos(dot) * 180 / np.pi
            #
            ## 4. 方向の符号（XY平面での2次元 cross product）
            ## cross2D = x1*y2 - y1*x2 → 正ならリード（天体が進行方向側）
            #cross2d = u_obs[:, 0] * u_sun[:, 1] - u_obs[:, 1] * u_sun[:, 0]
            #sign = np.sign(cross2d)
            #
            #signed_alpha = phase_angle * sign
            #signed_alpha = signed_alpha[0]
            #print(f"    Signed phase angle {signed_alpha:.3f} deg")
            # Get signed phase angle ==========================================

            r_list.append(r)
            alpha_list.append(alpha)
        r_mean = np.mean(r_list)
        alpha_mean = np.mean(alpha_list)
        print(f"  Mean (r, alpha) = ({r_mean:.3f}, {alpha_mean:.3f})")
        print("")
