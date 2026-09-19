#!/usr/bin/env python3
"""##CHRIS 2026-09-12: is non-convergence of the damped-cosine fit a MISSING measurement
or a SELECTION?

23 % of health-clean traces return no frequency from the time-domain fit, and at the
lightest divider mass it can be 0 of 10. That is only harmless if which traces fail is
independent of the frequency they would have given. The FFT-bin frequency exists for
EVERY trace, converged or not, so it is the instrument for the test:

  per cell   compare the binned-nu distribution of converged vs non-converged traces
             (medians, ratio, two-sample KS p-value)
  pooled     within each cell normalise binned nu by the converged median, then pool
             the non-converged ratios and test them against 1

If the two populations agree, the 23 % is a missing measurement. If they differ, the
surviving mean is biased and the PSD estimator (which needs no convergence) must become
primary immediately.
"""
import csv, math, sys
from collections import defaultdict
import numpy as np
from scipy import stats

WHICH = sys.argv[1]          # "A2" or "A1"
PATH = sys.argv[2]
MINN = int(sys.argv[3]) if len(sys.argv) > 3 else 3


def f(x):
    try: return float(x)
    except (TypeError, ValueError): return float("nan")


rows = list(csv.DictReader(open(PATH)))
cells = defaultdict(lambda: {"conv": [], "non": []})
n_health = 0
for r in rows:
    if WHICH == "A2":
        if any(int(r[k]) for k in ("forced_advance", "clamp_repair", "overlap_repair", "wall_overdue")):
            n_health += 1; continue
        key = (r["eta"], int(r["N"]), int(r["M"]))
        binned = f(r["nu_binned"]); fitted = r["nu_damped"]
    else:
        if r.get("eligible") not in (None, "", "1"):
            n_health += 1; continue
        key = (r["eta"], int(r["M"]))
        binned = f(r["nu"]); fitted = r.get("nu_fitted", "")
    if not (binned > 0):
        continue
    (cells[key]["conv"] if (fitted not in ("", "nan", None)) else cells[key]["non"]).append(binned)

tot_c = sum(len(v["conv"]) for v in cells.values())
tot_n = sum(len(v["non"]) for v in cells.values())
print(f"{WHICH}: health-discarded {n_health}   converged {tot_c}   non-converged {tot_n}"
      f"   ({100*tot_n/max(tot_c+tot_n,1):.1f} % non-converged)")

tested, flagged, pooled = 0, [], []
print(f"\n{'cell':>26} {'nC':>4} {'nX':>4} {'med conv':>12} {'med non':>12} {'ratio':>7} {'KS p':>9}")
for key in sorted(cells):
    c = np.array(cells[key]["conv"]); n = np.array(cells[key]["non"])
    if len(c) >= MINN and len(n) >= MINN:
        mc, mn = float(np.median(c)), float(np.median(n))
        p = float(stats.ks_2samp(c, n).pvalue)
        tested += 1
        lbl = ",".join(str(x) for x in key)
        star = "  <-- differs" if p < 0.01 else ""
        if p < 0.01: flagged.append((lbl, mn/mc, p))
        if p < 0.05 or abs(mn/mc - 1) > 0.10:
            print(f"{lbl:>26} {len(c):>4} {len(n):>4} {mc:>12.6g} {mn:>12.6g} {mn/mc:>7.3f} {p:>9.3g}{star}")
    if len(c) >= MINN and len(n) >= 1:
        mc = float(np.median(c))
        if mc > 0: pooled += list(n / mc)

print(f"\ncells with >= {MINN} in both groups: {tested};  KS p < 0.01 in {len(flagged)}")
print("  (cells whose fits ALL fail have no converged median, so this pooling is blind")
print("   to exactly the worst cells -- the predicted-frequency test below is not.)")

# ---- pooled test against nu_predicted, available for every trace ----
rc, rn = [], []
for r in rows:
    if WHICH == "A2":
        if any(int(r[k]) for k in ("forced_advance", "clamp_repair", "overlap_repair", "wall_overdue")): continue
        binned, fitted, pred = f(r["nu_binned"]), r["nu_damped"], f(r["nu_predicted"])
    else:
        if r.get("eligible") not in (None, "", "1"): continue
        binned, fitted, pred = f(r["nu"]), r.get("nu_fitted", ""), f(r.get("nu_predicted", "nan"))
    if not (binned > 0 and pred > 0): continue
    (rc if fitted not in ("", "nan", None) else rn).append(binned / pred)
if rc and rn:
    a, b = np.array(rc), np.array(rn)
    ks = stats.ks_2samp(a, b)
    mw = stats.mannwhitneyu(a, b, alternative="two-sided")
    print(f"\nbinned nu / nu_predicted, pooled over all cells:")
    print(f"   converged      n={len(a):>5}  median {np.median(a):.4f}  "
          f"IQR [{np.percentile(a,25):.4f}, {np.percentile(a,75):.4f}]")
    print(f"   non-converged  n={len(b):>5}  median {np.median(b):.4f}  "
          f"IQR [{np.percentile(b,25):.4f}, {np.percentile(b,75):.4f}]")
    d = 100*(np.median(b)/np.median(a) - 1)
    print(f"   median shift {d:+.2f} %   KS p = {ks.pvalue:.3g}   Mann-Whitney p = {mw.pvalue:.3g}")
    # restrict to the physical band so gross alias failures do not dominate
    ok_a = a[(a > 1/3.) & (a < 3)]; ok_b = b[(b > 1/3.) & (b < 3)]
    if len(ok_b) > 5:
        d2 = 100*(np.median(ok_b)/np.median(ok_a) - 1)
        ks2 = stats.ks_2samp(ok_a, ok_b)
        print(f"   restricted to ratio in [1/3, 3]: converged n={len(ok_a)} median {np.median(ok_a):.4f}; "
              f"non-converged n={len(ok_b)} median {np.median(ok_b):.4f}")
        print(f"   median shift {d2:+.2f} %   KS p = {ks2.pvalue:.3g}")
        print(f"   fraction of non-converged OUTSIDE the band: {100*(1-len(ok_b)/len(b)):.1f} % "
              f"(converged: {100*(1-len(ok_a)/len(a)):.1f} %)")
if pooled:
    a = np.array(pooled)
    med = float(np.median(a))
    # sign test of the pooled ratios against 1
    k = int((a > 1).sum()); nn = len(a)
    p_sign = float(stats.binomtest(k, nn, 0.5).pvalue)
    w = float(stats.wilcoxon(a - 1.0).pvalue) if nn > 10 else float("nan")
    print(f"pooled non-converged binned-nu / converged-cell-median:")
    print(f"   n = {nn}   median = {med:.4f}   IQR = "
          f"[{np.percentile(a,25):.4f}, {np.percentile(a,75):.4f}]")
    print(f"   above 1: {k}/{nn}   sign-test p = {p_sign:.3g}   Wilcoxon p = {w:.3g}")
    bias = 100 * (med - 1)
    print(f"\nVERDICT: non-converged traces sit {bias:+.2f} % from the converged median.")
    if len(flagged) > 0.05 * max(tested, 1) or abs(bias) > 1.0:
        print("  -> DIFFERENCE BEYOND NOISE. Non-convergence is a SELECTION, not a missing")
        print("     measurement. The PSD estimator must become primary.")
    else:
        print("  -> consistent with a missing measurement, not a selection.")
