#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Compare single and dual component models with F-test.
"""
import os 
from argparse import ArgumentParser as ap
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import f

from Eros_common import mycolor


def plot_f_distribution_with_observed(
    nu1, nu2, F_obs, alpha=0.05, xmax=10, out=None):
    """Compare observed and expected F-ratios.
    
    Parameters
    ----------
    nu1 : float
        degrees of freedom of numerator
    nu2 : float
        degrees of freedom of denominator
    F_obs : float
        observed F-ratio
    alpha : float, optional
        confidence level
    xmax : float, optional
        maximum range of x axis
    out : str, optional
        output filename
    """
    # PDF
    x = np.linspace(0, xmax, 100)
    pdf = f.pdf(x, nu1, nu2)
    
    # Calculate expected F, F_crit, and corresponding p value
    F_crit = f.ppf(1 - alpha, nu1, nu2)
    p_value = 1 - f.cdf(F_obs, nu1, nu2)
    
    # Plot 
    plt.figure(figsize=(9, 5))
    plt.plot(
        x, pdf, color="black", 
        label=r"F-distribusion ($\nu_1, \nu_2$) =" + f"({nu1}, {nu2})")
    
    plt.axvline(
        F_obs, color=mycolor[0], linestyle='--', 
        label=f'Observed F = {F_obs:.3f}')
    plt.axvline(
        F_crit, color=mycolor[1], linestyle=':', 
        label=f'Critical F (α={alpha}) = {F_crit:.3f}')
    
    # Make shadow 
    if F_obs < xmax:
        x_fill = np.linspace(F_obs, xmax, 200)
        plt.fill_between(
            x_fill, f.pdf(x_fill, nu1, nu2), color=mycolor[0], alpha=0.3, 
            label=f'p-value area corresponds to observed F ({p_value:.4f})')
    
    plt.xlabel('F value')
    plt.ylabel('PDF')
    plt.legend()
    plt.grid(True)
    
    if out:
        plt.savefig(out)
    else:
        plt.show()


if __name__ == "__main__":
    parser = ap(
        description="Compare simple/complex models with F-test.")
    parser.add_argument(
        "Ndata", type=int,
        help="Number of data points")
    parser.add_argument(
        "Nparam_s", type=int,
        help="Number of parameters for simple model")
    parser.add_argument(
        "Nparam_c", type=int,
        help="Number of parameters for complex model")
    parser.add_argument(
        "chi2_s", type=float,
        help="Chi-squared statistics for simple model")
    parser.add_argument(
        "chi2_c", type=float,
        help="Chi-squared statistics for complex model")
    parser.add_argument(
        "--out", type=str, default="Ftest.pdf",
        help="Output filename")
    args = parser.parse_args()

    # Confidence level
    alpha = 0.01
    # Maximum x range
    xmax = 6

    # Number of data: N
    # Simgle model  (N_param=p)
    # Complex model (N_param=q)
    N = args.Ndata
    p = args.Nparam_s
    q = args.Nparam_c
    chi2_s = args.chi2_s
    chi2_c = args.chi2_c
    # Output file
    out = args.out

    # Degrees of freedom of numerator
    # q-p 
    nu1 = q - p
    # Degrees of freedom of denominator
    # N-q 
    nu2 = N - q

    F_obs = ((chi2_s - chi2_c)/nu1) / (chi2_c/nu2)
    print(f"Observed F-ratio = {F_obs:.2f}")
    print(f"DoFs: (nu1, nu2) = ({nu1}, {nu2})")
    print(f"Confidence level = {alpha}")
    
    plot_f_distribution_with_observed(nu1, nu2, F_obs, alpha=alpha, xmax=xmax, out=out)

