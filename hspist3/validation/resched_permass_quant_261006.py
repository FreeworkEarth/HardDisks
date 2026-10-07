#!/usr/bin/env python3
"""##CHRIS 2026-10-06 (261012 sec. 4.4.9): per-mass old (279282b) vs new (73fc07f replay) divider frequency of the two replayed
method-B cells, with the quantization of the spectral estimator: how many distinct nu values the 25 seeds take, the frequency-bin
width df/nu, and the per-seed spread in bins. A per-seed spread of 1-2 bins makes the per-mass mean lumpy, so Gaussian p-values
of the per-mass z are optimistic. Analysis only. usage (from hspist3/): python3 validation/resched_permass_quant_261006.py"""
import os, numpy as np, pandas as pd
from scipy.stats import chi2
HS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RB = "experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/confinement_B_20261013"
NEW = os.path.join(HS, "experiments_resched_gate_261005")
for cid in ("e0p10_H_H10_L39.25", "epi8_H_H10_L10"):
    zs = []
    print(f"\n{cid}\n| M | old distinct nu (of 25) | new distinct | bin df/nu [%] | per-seed sd / df | z (new - old) |\n|---|---|---|---|---|---|")
    for M in (50, 100, 200, 300, 500, 750, 1000, 1500, 2000):
        o = pd.read_csv(os.path.join(HS, RB, cid, f"m_{M}", "red_nu.csv")).nu; n = pd.read_csv(os.path.join(NEW, RB, cid, f"m_{M}", "red_nu.csv")).nu
        df = float(np.median(np.diff(np.unique(np.r_[o, n]))))
        z = (n.mean() - o.mean()) / np.hypot(o.std(ddof=1) / np.sqrt(len(o)), n.std(ddof=1) / np.sqrt(len(n))); zs.append(z)
        print(f"| {M} | {o.nunique()} | {n.nunique()} | {100 * df / o.mean():.2f} | {o.std(ddof=1) / df:.2f} | {z:+.2f} |")
    zs = np.array(zs); print(f"chi2 of the 9 per-mass z = {np.sum(zs ** 2):.2f} / 9, nominal p = {chi2.sf(np.sum(zs ** 2), 9):.3f}")
