#!/usr/bin/env python3
"""##CHRIS 2026-09-15: is the dense-eta difference between the velocity-spectrum peak and the position-spectrum
peak a damping shift of the position peak, or contamination of the velocity spectrum?

Per trajectory (A1 v2, first 200 predicted oscillations, mean-subtracted, no transient drop):
  position peak   largest bin of P_x at f >= nu_pred/2.5           (Paper 1 primary)
  velocity peak   largest bin of (2 pi f)^2 P_x at k >= 2
  resonance fit   Whittle fit of S(u) = A/((u^2-u0^2)^2+(g u)^2)+B on [nu_pred/2, 2 nu_pred],
                  u = f/nu_pred: gives f_0 = u0 nu_pred and Gamma = g nu_pred, i.e. Gamma/f_0 = g/u0
For a linearly damped oscillator under flat forcing the position spectrum peaks at
f_0 sqrt(1 - Gamma^2/(2 f_0^2)), i.e. below f_0 by Gamma^2/(4 f_0^2) to leading order, while the velocity
spectrum peaks at f_0 exactly. c_s per eta = through-origin slope over the mass ladder of each estimator."""
import os, sys, math, json
import numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import tests_20260913 as T

TD, X_EDGE, FITW = 200, 2.5, 2.0


def fit_one(P, df, nup, w=FITW):
    """-> (u0, g) in units of nu_pred, or (nan, nan)."""
    from scipy.optimize import minimize
    idx = np.arange(len(P)); f = idx * df
    sel = (idx >= 1) & (f >= nup / w) & (f <= nup * w)
    if sel.sum() < 6:
        return float("nan"), float("nan")
    u = f[sel] / nup; I = P[sel] / P[sel].max()
    B0 = max(float(np.percentile(I, 10)), 1e-6)

    def nll(p):
        lA, u0, lg, lB = p
        if not (1 / w < u0 < w):
            return 1e30
        S = math.exp(lA) / ((u * u - u0 * u0) ** 2 + (math.exp(lg) * u) ** 2) + math.exp(lB)
        if not np.all(np.isfinite(S)) or np.any(S <= 0):
            return 1e30
        return float(np.sum(np.log(S) + I / S))

    best = None
    for u0s in (float(u[int(np.argmax(I))]), 1.0):
        g0 = 0.1
        A0 = max(1.0 - B0, 1e-6) * (g0 * u0s) ** 2
        r = minimize(nll, [math.log(A0), u0s, math.log(g0), math.log(B0)], method="Nelder-Mead",
                     options=dict(maxiter=6000, xatol=1e-7, fatol=1e-9))
        if best is None or r.fun < best.fun:
            best = r
    u0, g = float(best.x[1]), math.exp(float(best.x[2]))
    if not (1 / w * 1.02 < u0 < w / 1.02):
        return float("nan"), float("nan")
    return u0, g


def cell(task):
    eta, L0, M, runs = task
    pos, vel, f0, gf, nfail = [], [], [], [], 0
    for r, p, disc in runs:
        if disc:
            continue
        t, x, nup = T._load(p); dt = (t[-1] - t[0]) / (len(t) - 1)
        n = T._prefix(t, nup, TD); P, df = T._spectrum(x[:n], dt)
        k = int(round(TD / X_EDGE)); pos.append((k + int(np.argmax(P[k:]))) * df)
        f = np.arange(len(P)) * df
        Pv = (2 * np.pi * f) ** 2 * P
        vel.append((2 + int(np.argmax(Pv[2:]))) * df)
        u0, g = fit_one(P, df, nup)
        if np.isfinite(u0):
            f0.append(u0 * nup); gf.append(g / u0)
        else:
            nfail += 1
    return dict(eta=eta, L0=L0, M=M, n=len(pos), nfail=nfail,
                pos=float(np.mean(pos)), vel=float(np.mean(vel)),
                f0=float(np.mean(f0)) if f0 else float("nan"),
                gamma_over_f0=float(np.median(gf)) if gf else float("nan"))


if __name__ == "__main__":
    table = T.a1_leaf_table(); tasks = []
    for leaf in table:
        for M in T.A1_MASSES:
            runs = T.cell_runs(os.path.join(T.DROOT, leaf["leaf"], f"m_{M}"), M)
            if runs:
                tasks.append((leaf["eta"], float(leaf["L0"]), M, runs))
    with Pool(10) as pool:
        cells = pool.map(cell, tasks, chunksize=1)
    json.dump(cells, open(os.path.join(T.PLOTS, "260915_A1v2_damping_cells.json"), "w"))
    rows = []
    for leaf in table:
        eta = leaf["eta"]; cs = [c for c in cells if c["eta"] == eta]; L0 = float(leaf["L0"])
        xs = [T.x_of(c["M"], L0) for c in cs]
        cp = T.slope(xs, [c["pos"] for c in cs])[0]; cv = T.slope(xs, [c["vel"] for c in cs])[0]
        cf, sf, _ = T.slope(xs, [c["f0"] for c in cs])
        gam = float(np.median([c["gamma_over_f0"] for c in cs if np.isfinite(c["gamma_over_f0"])]))
        rows.append(dict(eta=f"{eta:.6f}", L0=leaf["L0"], c_s_position=cp, c_s_velocity=cv, c_s_fit_f0=cf,
                         fit_scatter=sf, dev_vel_pct=100 * (cv - cp) / cp, dev_fit_pct=100 * (cf - cp) / cp,
                         median_gamma_over_f0=gam, predicted_shift_pct=100 * gam ** 2 / 4,
                         n_fit_failed=sum(c["nfail"] for c in cs), n_used=sum(c["n"] for c in cs)))
    T.write_csv(os.path.join(T.PLOTS, "260915_A1v2_damping_test.csv"), rows)
    print("| η | c_s position (X = 2.5) | c_s velocity peak | c_s fitted f_0 | (vel − pos)/pos [%] | (fit − pos)/pos [%] | predicted shift Γ²/4f_0² [%] | median Γ/f_0 | fits failed / used |")
    print("|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        print(f"| {r['eta']} | {r['c_s_position']:.4f} | {r['c_s_velocity']:.4f} | {r['c_s_fit_f0']:.4f} | "
              f"{r['dev_vel_pct']:+.3f} | {r['dev_fit_pct']:+.3f} | {r['predicted_shift_pct']:.3f} | {r['median_gamma_over_f0']:.4f} | "
              f"{r['n_fit_failed']} / {r['n_used']} |")
    for lab, sub in (("η < 0.5", [r for r in rows if float(r["eta"]) < 0.5]),
                     ("η ≥ 0.5", [r for r in rows if float(r["eta"]) >= 0.5]),
                     ("0.5 ≤ η ≤ 0.69", [r for r in rows if 0.5 <= float(r["eta"]) <= 0.69])):
        v = np.array([r["dev_vel_pct"] for r in sub]); fdev = np.array([r["dev_fit_pct"] for r in sub])
        p = np.array([r["predicted_shift_pct"] for r in sub]); g = np.array([r["median_gamma_over_f0"] for r in sub])
        print(f"{lab}: median (vel−pos) {np.median(v):+.3f} %, median (fit−pos) {np.median(fdev):+.3f} %, "
              f"median predicted shift {np.median(p):.3f} %, median Γ/f_0 {np.median(g):.3f}; "
              f"|fit−vel| median {np.median(np.abs(fdev - v)):.3f} %, |fit−pos| median {np.median(np.abs(fdev)):.3f} %")
