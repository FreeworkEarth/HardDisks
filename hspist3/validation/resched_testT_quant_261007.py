#!/usr/bin/env python3
"""##CHRIS 2026-10-07 (261012 sec. 4.4.12): INFORMATION ONLY, after Test T's verdict (FAIL by the registered rule); nothing here
changes or re-tests that verdict. Per mass of the Test T cell (epi8_H_H10_L10, minimal vs legacy, same binary): the FFT bin of
the per-seed nu (the smallest spacing of its distinct values), the mean difference in bins, the per-seed SD of each policy and
their ratio with the two-sided F-test p, the number of distinct values, and the same-seed correlation (same seed on both
policies; chaotic divergence makes it ~0 if the policies share nothing but the initial state).
usage (from hspist3/): python3 validation/resched_testT_quant_261007.py
"""
import os
import numpy as np, pandas as pd
from scipy.stats import f as F

HS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOC = os.environ.get("HD_RESCHED2_LOC", os.path.join(HS, "experiments_resched_gate2_261007"))
REL = "experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/resched_testT_261007"
MASSES = (50, 100, 200, 300, 500, 750, 1000, 1500, 2000)


def main():
    print("# Test T, information only: quantization and spread per mass (minimal vs legacy, same seeds)\n")
    print("| M | n per policy | FFT bin of nu | mean difference [bins] | SD minimal [bins] | SD legacy [bins] | SD ratio min/leg | "
          "F-test p (two-sided) | distinct values min / leg | same-seed correlation |\n|---|---|---|---|---|---|---|---|---|---|")
    for M in MASSES:
        a, b = (pd.read_csv(os.path.join(LOC, REL, p, "epi8_H_H10_L10", f"m_{M}", "red_nu.csv")).sort_values("run")["nu"].to_numpy(float)
                for p in ("minimal", "legacy"))
        v = np.unique(np.round(np.concatenate([a, b]), 12)); bin_ = float(np.diff(v).min())
        sa, sb = a.std(ddof=1), b.std(ddof=1); n = len(a)
        fr = sa * sa / (sb * sb); p = 2 * min(F.cdf(fr, n - 1, n - 1), F.sf(fr, n - 1, n - 1))
        print(f"| {M} | {n} | {bin_:.4e} | {(a.mean() - b.mean()) / bin_:+.3f} | {sa / bin_:.2f} | {sb / bin_:.2f} | {sa / sb:.3f} | "
              f"{p:.3f} | {len(np.unique(np.round(a, 12)))} / {len(np.unique(np.round(b, 12)))} | {np.corrcoef(a, b)[0, 1]:+.3f} |")


if __name__ == "__main__":
    main()
