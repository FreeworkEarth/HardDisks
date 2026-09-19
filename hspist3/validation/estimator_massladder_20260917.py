#!/usr/bin/env python3
"""##CHRIS 2026-09-17: which frequency estimator does the mass ladder itself prefer?

No equation of state enters. Roman's relation nu(M) = c_s K(M/N) / (2 pi L_eff) demands that ONE slope
put every divider mass on one line through the origin, so the estimator whose per-mass frequencies
scatter less about that line is the better supported one. Two estimators are compared:

  peak : largest FFT bin at f >= nu_pred / 2.5  (the primary rule, section 1 of the method note)
  fit  : per-trajectory resonance fit S(u) = A/((u^2-u0^2)^2 + (g u)^2) + B  (the cross-check)

Reported against the bin-rounding floor 100/(N_cyc sqrt(12)) %, which is what the peak estimator
cannot beat if the maximum locks to the same bin in every seed of a cell.
Analysis only; reads existing traces.
"""
import glob, json, math, os, sys
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tests_20260913 as T
from damping_test_20260915 import fit_one

DILUTE = os.path.join(T.ROOT, "A2_dilute50_20260917")


def rms_about_line(x, y):
    """rms relative residual about the best through-origin line, in percent."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    s = float((x * y).sum() / (x * x).sum())
    return 100.0 * float(np.sqrt(np.mean(((y - s * x) / (s * x)) ** 2)))


def a1v2():
    """9 masses x 25 seeds x 200 periods, from the damping-test cell cache."""
    cells = json.load(open(os.path.join(T.PLOTS, "260915_A1v2_damping_cells.json")))
    pk, ft, win = [], [], {"peak": 0, "fit": 0}
    for leaf in T.a1_leaf_table():
        if leaf["eta"] > 0.69:
            continue
        cs = [c for c in cells if c["eta"] == leaf["eta"]]
        if len(cs) < 5:
            continue
        x = np.array([T.x_of(c["M"], float(leaf["L0"])) for c in cs])
        got = {}
        for key, lab in (("pos", "peak"), ("f0", "fit")):
            y = np.array([c[key] for c in cs]); ok = np.isfinite(y)
            if ok.sum() >= 5:
                got[lab] = rms_about_line(x[ok], y[ok])
        if len(got) == 2:
            pk.append(got["peak"]); ft.append(got["fit"])
            win["fit" if got["fit"] < got["peak"] else "peak"] += 1
    return pk, ft, win


def dilute():
    """5 masses x 10 seeds x 50 periods, per (eta, N) cell, re-analysed from the traces."""
    pk, ft, lock, rows = [], [], [], []
    for tag, eta in (("eta_0p02", 0.02), ("eta_0p05", 0.05)):
        for N in (100, 400, 900, 1600):
            xs, P, F = [], [], []
            for M in (50, 200, 500, 1000, 2000):
                p_, f_, L0 = [], [], None
                for p in sorted(glob.glob(f"{DILUTE}/{tag}/N{N}/m_{M}/r*/wall_x_positions_*run0.csv")):
                    d = pd.read_csv(p, usecols=["Time", "Displacement(σ)", "Predicted_Frequency", "L0"])
                    t = d["Time"].to_numpy(float); x = d["Displacement(σ)"].to_numpy(float)
                    nup = float(d["Predicted_Frequency"].iloc[0]); L0 = float(d["L0"].iloc[0])
                    dt = (t[-1] - t[0]) / (len(t) - 1)
                    Pw, df = T._spectrum(x, dt)
                    k = max(1, int(round(len(x) * dt * nup / 2.5)))
                    p_.append((k + int(np.argmax(Pw[k:]))) * df)
                    u0, g = fit_one(Pw, df, nup)
                    if np.isfinite(u0):
                        f_.append(u0 * nup)
                if len(p_) >= 3:
                    xs.append(T.k_root(M / N) / (2 * math.pi * T.l_eff(L0)))
                    P.append(np.mean(p_)); F.append(np.mean(f_))
                    lock.append(len(set(np.round(p_, 12))) / len(p_))
            a, b = rms_about_line(xs, P), rms_about_line(xs, F)
            pk.append(a); ft.append(b); rows.append((eta, N, a, b))
    return pk, ft, lock, rows


def main():
    print("### Which estimator puts the masses on one line? (no EOS involved)\n")
    pk1, ft1, win = a1v2()
    pk2, ft2, lock, rows = dilute()
    print("| record | N_cyc | rounding floor [%] | median rms residual, largest bin [%] | median rms residual, fit [%] |")
    print("|---|---|---|---|---|")
    print(f"| A1 v2 (24 η ≤ 0.69, 9 masses × 25 seeds) | 200 | {100 / (200 * math.sqrt(12)):.3f} | "
          f"{np.median(pk1):.3f} | {np.median(ft1):.3f} |")
    print(f"| A2 dilute (8 cells, 5 masses × 10 seeds) | 50 | {100 / (50 * math.sqrt(12)):.3f} | "
          f"{np.median(pk2):.3f} | {np.median(ft2):.3f} |")
    print(f"\nA1 v2 per-density winner: fit {win['fit']}, largest bin {win['peak']} of {win['fit'] + win['peak']}")
    print(f"distinct peak bins per 10-seed cell (1.0 = all different, 0.1 = all locked): median {np.median(lock):.2f}\n")
    print("| η | N | rms residual, largest bin [%] | rms residual, fit [%] | winner |")
    print("|---|---|---|---|---|")
    for eta, N, a, b in rows:
        print(f"| {eta:.2f} | {N} | {a:.3f} | {b:.3f} | {'fit' if b < a else 'peak'} |")


if __name__ == "__main__":
    main()


# ##CHRIS 2026-09-17, second pass: the diagnostics ChatGPT asked for on top of plain scatter --
# normalised residual against alpha = M/N (a bias shows up as a trend, not as noise) and
# leave-one-mass-out stability of c_s. Run with `python estimator_massladder_20260917.py extra`.
def extra():
    cells = json.load(open(os.path.join(T.PLOTS, "260915_A1v2_damping_cells.json")))
    print("\n### Residual against α = M/N: a biased estimator leans, an unbiased one scatters")
    print("| estimator | mean residual per α, over 24 densities [%] | Spearman-like trend (first α − last α) [%] |")
    print("|---|---|---|")
    alphas = sorted({c["M"] / 100 for c in cells})
    for key, lab in (("pos", "largest bin"), ("f0", "per-trajectory fit")):
        by = {a: [] for a in alphas}
        for leaf in T.a1_leaf_table():
            if leaf["eta"] > 0.69:
                continue
            cs = [c for c in cells if c["eta"] == leaf["eta"]]
            if len(cs) < 5:
                continue
            x = np.array([T.x_of(c["M"], float(leaf["L0"])) for c in cs])
            y = np.array([c[key] for c in cs]); ok = np.isfinite(y)
            if ok.sum() < 5:
                continue
            s = float((x[ok] * y[ok]).sum() / (x[ok] ** 2).sum())
            for c, xi, yi, o in zip(cs, x, y, ok):
                if o:
                    by[c["M"] / 100].append(100 * (yi - s * xi) / (s * xi))
        m = [float(np.mean(by[a])) for a in alphas]
        print(f"| {lab} | " + ", ".join(f"{v:+.2f}" for v in m) + f" | {m[0] - m[-1]:+.2f} |")
    print(f"(α = " + ", ".join(f"{a:g}" for a in alphas) + ")")

    print("\n### Leave-one-mass-out: how much does c_s move when any single mass is dropped?")
    print("| estimator | median spread of c_s over the 9 leave-one-out fits [%] | worst density [%] |")
    print("|---|---|---|")
    for key, lab in (("pos", "largest bin"), ("f0", "per-trajectory fit")):
        sp = []
        for leaf in T.a1_leaf_table():
            if leaf["eta"] > 0.69:
                continue
            cs = [c for c in cells if c["eta"] == leaf["eta"]]
            if len(cs) < 5:
                continue
            x = np.array([T.x_of(c["M"], float(leaf["L0"])) for c in cs])
            y = np.array([c[key] for c in cs]); ok = np.isfinite(y)
            if ok.sum() < 6:
                continue
            xx, yy = x[ok], y[ok]
            loo = [float((np.delete(xx, i) * np.delete(yy, i)).sum() / (np.delete(xx, i) ** 2).sum())
                   for i in range(len(xx))]
            sp.append(100 * (max(loo) - min(loo)) / np.mean(loo))
        print(f"| {lab} | {np.median(sp):.3f} | {max(sp):.3f} |")


if len(sys.argv) > 1 and sys.argv[1] == "extra":
    extra()
