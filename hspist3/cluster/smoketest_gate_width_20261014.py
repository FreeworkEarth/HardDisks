#!/usr/bin/env python3
"""##CHRIS 2026-10-14: width of the KOA smoke-test c_s gate (261012 sec. 1.10.1). Reads the Mac pi/8 pilot only.

0.05150 is np.std(imp, ddof=1) at cluster/confinement_pilot.py:64 -- the sample STANDARD DEVIATION of the nine
per-mass implied c_s = nu/x, i.e. the scatter, not the standard error of the pilot value.

The pilot value is the through-origin slope c_s = sum(x nu)/sum(x^2) = sum(w c_i)/sum(w) with w = x^2, a weighted mean
of the implied c_i. With a common scatter s, SE = s sqrt(sum w^2)/sum w (DERIVATION). Two independent pilots (Mac, KOA)
with the same scatter differ with sigma_diff = sqrt(2) SE; the gate is set at 2 sigma_diff.
"""
import math, os, sys
import numpy as np
from scipy.stats import norm
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(HS, "validation")); sys.path.insert(0, HS)
import tests_20260913 as T
from paper1_populate_cs_err_20261002 import cell
OUT = os.path.join(HS, "experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/confinement_pilot_20261013/mac_pi8_H10_L10")

def main():
    x, nu = [], []
    for M in T.A1_MASSES:
        c = cell((0.392699, 10.0, M, T.cell_runs(os.path.join(OUT, f"m_{M}"), M))); x.append(T.x_of(M, 10.0)); nu.append(c["nu"])
    x, nu = np.array(x), np.array(nu); imp = nu / x; w = x * x
    cs = float(np.sum(x * nu) / np.sum(w)); s = float(np.std(imp, ddof=1))
    se_slope = s * math.sqrt(np.sum(w * w)) / np.sum(w); se_mean = s / math.sqrt(len(imp)); neff = np.sum(w) ** 2 / np.sum(w * w)
    sd = math.sqrt(2) * se_slope; gate = 2 * sd
    print("### KOA smoke-test gate width, from the Mac pi/8 pilot\n")
    print(f"- pilot c_s (through-origin slope)                  = {cs:.5f}   (reproduces the header's 3.85886)")
    print(f"- s = np.std(imp, ddof=1), confinement_pilot.py:64  = {s:.5f}   -> the 0.05150 is the SCATTER (SD) across the 9 masses")
    print(f"- slope weights w = x^2 (share per mass, M = {T.A1_MASSES[0]} ... {T.A1_MASSES[-1]}): "
          + ", ".join(f"{v:.3f}" for v in w / w.sum()) + f";  effective n = {neff:.2f}")
    print(f"- SE of the pilot c_s (slope)  = s sqrt(sum w^2)/sum w = {se_slope:.5f}   (an unweighted mean would have s/3 = {se_mean:.5f})")
    print(f"- sigma of (KOA - Mac), two independent pilots     = sqrt(2) SE = {sd:.5f}")
    print(f"- **new gate: |c_s(KOA) - c_s(Mac)| <= 2 sigma_diff = {gate:.5f}**")
    print(f"- expected false-fail probability under the null (Gaussian, same scatter on KOA): 2(1 - Phi(2)) = {2*(1-norm.cdf(2)):.4f}")
    print(f"- the old gate 0.05150 sat at {0.05150/sd:.2f} sigma_diff; its false-fail probability was {2*(1-norm.cdf(0.05150/sd)):.4f}")
    print("- caveat (stated, not corrected): s is estimated from 9 single-seed values, so sigma_diff itself is uncertain by about "
          f"1/sqrt(2*8) = {1/math.sqrt(16):.2f} (relative); the per-mass frequencies are quantised by the 200-period spectral bin.")

if __name__ == "__main__":
    main()
