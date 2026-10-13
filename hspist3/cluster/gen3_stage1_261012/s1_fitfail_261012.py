#!/usr/bin/env python3
"""##CHRIS 2026-10-13 (plan-author decision 45 item 5 b, addendum 1 item 2; registration 261012 sec. 4.8.5 (group rule), sec. 4.8.6.1-4.8.6.2):
THE nu_d FIT FAILURES OF A FETCHED GROUP, on the Mac, from fetched files only (nothing on KOA is read, nothing runs there, no simulation runs;
nothing is written into the fetched folder):
  per (cell, mass) folder with K or more fit failures (nud_runs.csv's fit_ok 0, written on KOA by reduce_testR_koa.nud_table; K = --min):
    - the cell, M and alpha; expected, present, excluded (the health rule), missing, fit failures, used and lost trajectories, and whether
      the group rule flags it (more than 2 of 25 lost: s1_verdict_261012.per_mass, the function of the INTERIM tables);
    - Q of the folder: Q_k = pi nu_d-bar tau_r of E1 over its used seeds (s1_verdict_261012.per_mass), and Q_meanfit of the KOA reduction
      over all its trajectories (acf_mean_fit.csv);
    - per failed trajectory: run s, exact seed, nu_arg, and WHICH CRITERION FAILED. The registered fit (validation/
      resched_testT_followup_261007.fit_acf: scipy curve_fit, the item-4a model, bounds, maxfev 200000) has one failure criterion: it
      raises. It is repeated here on the fetched ACF with the KOA reduction's inputs (the folder's ACFs cut to the shortest, its dt,
      nu_ref = the folder's argmax mean), and the exception (type and message) is printed. As information: scipy least_squares with the
      fit's model, start, bounds and limit (curve_fit's own call): its status, evaluations, last iterate (A, tau_T, B, tau_r, omega, phi),
      the bounds it ends on, nu_d = omega / 2 pi and Q = omega tau_r / 2 there; and the ACF's sign changes in its first 10 periods of
      nu_ref (resched_testT_followup_261007.crossings);
    - the repetition check: fit_ok of every trajectory of the folder on the Mac against nud_runs.csv, and the largest |nu_d (Mac) /
      nu_d (KOA) - 1| over the trajectories fitted on both;
  per cell with a flagged folder: the masses left for E1 (alpha >= 1, not flagged, its seed-averaged fit done) and whether they are 5 or
  more (sec. 4.8.6.2: "A cell with fewer than 5 masses left has no E1 value").
usage (from hspist3/): python3 cluster/gen3_stage1_261012/s1_fitfail_261012.py <group> [--root R] [--min K]
"""
import argparse, glob, math, os, re, sys, time, warnings
import numpy as np, pandas as pd, scipy
from scipy.optimize import least_squares
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HS, "validation")); sys.path.insert(0, HS)
import s1_verdict_261012 as V                                    # registered(), read_mass(), per_mass(): the INTERIM tables' functions
import resched_testT_followup_261007 as FU                       # fit_acf, _g, crossings: the item-4a fit
warnings.filterwarnings("ignore")
NAMES = ("A", "tau_T", "B", "tau_r", "omega", "phi")


def inputs(d):
    """reduce_testR_koa.nud_table's inputs, read only: {run: ACF cut to the shortest}, dt, nu_ref"""
    red = pd.read_csv(os.path.join(d, "red_nu.csv")).sort_values("run").reset_index(drop=True)
    z = np.load(os.path.join(d, "acf_runs.npz")); nl = min(len(z[f"run{r}"]) for r in red["run"])
    return {int(r): np.asarray(z[f"run{r}"][:nl], dtype=float) for r in red["run"]}, float(red["dt"].iloc[0]), float(red["nu"].mean()), nl


def criterion(msg):
    m = msg.lower()
    if "function evaluations" in m or "maxfev" in m: return "no convergence within maxfev = 200000 evaluations"
    if "not finite" in m or "infs or nans" in m: return "a value of the ACF or of the model at the start is not finite"
    if "infeasible" in m: return "the start lies outside the bounds"
    if "svd" in m: return "the covariance step failed (SVD)"
    return "other"


def lsq(a, dt, nu_ref):
    """curve_fit's call inside fit_acf, by hand: the last iterate even when it does not converge"""
    c = a / a[0]; lag = np.arange(len(c)) * dt; om0 = 2 * math.pi * nu_ref; P = 1 / nu_ref
    lo = [0, 2 * dt, 0, 2 * dt, 0.2 * om0, -math.pi / 2]; hi = [1.5, 1e7, 1.5, 1e7, 5 * om0, math.pi / 2]
    r = least_squares(lambda p: FU._g(lag, *p) - c, [0.3, 2 * P, 0.7, 5 * P, om0, 0.0], jac="2-point", bounds=(lo, hi), method="trf",
                      max_nfev=200000)
    on = [f"{NAMES[i]} {'lower' if r.active_mask[i] < 0 else 'upper'}" for i in range(6) if r.active_mask[i] != 0]
    return r, on


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("group"); ap.add_argument("--root", default=os.path.join(HS, "experiments_gen3_stage1_koa_261012"))
    ap.add_argument("--min", type=int, default=3); a = ap.parse_args()
    base = os.path.join(a.root, a.group, "merged", "stage1", "gen3")
    if not os.path.isdir(base): sys.exit(f"STOP: no fetched folder {base}")
    C = V.registered()
    print(f"# INTERIM -- group {a.group}: the nu_d fit failures (s1_fitfail_261012.py; the fit repeated on the Mac: python {sys.version.split()[0]}, "
          f"scipy {scipy.__version__}, numpy {np.__version__})\n")
    allf = []
    for d in sorted(glob.glob(os.path.join(base, "*", "m_*"))):
        p = os.path.join(d, "nud_runs.csv")
        if not os.path.exists(p): continue
        n = pd.read_csv(p); nf = int((n["fit_ok"] == 0).sum())
        allf.append((nf, os.path.basename(os.path.dirname(d)), int(os.path.basename(d)[2:]), d))
    hist = {}
    for nf, *_ in allf: hist[nf] = hist.get(nf, 0) + 1
    print(f"folders with nud_runs.csv: {len(allf)}; by fit failures: " + ", ".join(f"{hist[k]} with {k}" for k in sorted(hist)) + "\n")
    sel = sorted((x for x in allf if x[0] >= a.min), key=lambda x: (-x[0], V.GEN.GROUPS.index(a.group), x[1], x[2]))
    if not sel: print(f"no folder with {a.min} or more fit failures"); return
    flagged_cells = {}
    for nf, cid, M, d in sel:
        cell = C[cid]; D = V.read_mass(d); r = V.per_mass(cell, M, D, jackknife=False)
        mf = pd.read_csv(os.path.join(d, "acf_mean_fit.csv")).iloc[0] if os.path.exists(os.path.join(d, "acf_mean_fit.csv")) else None
        print(f"## {cid}, M = {M} (alpha {cell['masses'][M]:g}; N {cell['N']}, eta {cell['eta']:.6f}): {nf} fit failures\n")
        print(f"- trajectories: expected {r['expected']}, present {r['present']}, excluded {r['excluded']} (health rule), missing {r['missing']}, "
              f"fit failures {r['fitfail']}, used {r['used']}, lost {r['lost']}: the group rule "
              + ("FLAGS it (more than 2 lost): it is left out of the cell's fits" if r["flagged"] else "does not flag it (2 or fewer lost)"))
        q1 = f"{r['Q']:.2f}" if math.isfinite(r.get("Q", math.nan)) else "-"
        q2 = (f"{mf['Q_meanfit']:.2f} (fit_ok {int(mf['fit_ok'])}, n {int(mf['n'])})" if mf is not None else "no acf_mean_fit.csv")
        print(f"- Q of the folder: Q_k = pi nu_d-bar tau_r over its {r['used']} used seeds = {q1}; the KOA reduction's Q_meanfit over all "
              f"{r['present']} trajectories = {q2}")
        if r["flagged"]: flagged_cells.setdefault(cid, []).append(M)
        acf, dt, nu_ref, nl = inputs(d)
        koa = pd.read_csv(os.path.join(d, "nud_runs.csv")).set_index("run")
        print(f"- the fit's inputs: {len(acf)} ACFs cut to {nl} lags, dt {dt:.6g}, nu_ref {nu_ref:.7g} (the folder's argmax mean)\n")
        print("| run | seed | nu_arg | the registered fit on the Mac | criterion | least_squares (information): status, evaluations | last iterate: "
              "nu_d / nu_ref, tau_r, Q = omega tau_r / 2, A, B, tau_T, phi | bounds it ends on | sign changes (10 periods) |\n|---|---|---|---|---|---|---|---|---|")
        same, dev, t0 = 0, [], time.time()
        for run, x in acf.items():
            try:
                nd = FU.fit_acf(x, dt, nu_ref)[4] / (2 * math.pi); err = None
            except Exception as e:
                nd = math.nan; err = f"{type(e).__name__}: {str(e).strip()}"
            k_ok = int(koa.loc[run, "fit_ok"]); m_ok = int(err is None)
            same += (k_ok == m_ok)
            if k_ok and m_ok: dev.append(abs(nd / float(koa.loc[run, "nu_d"]) - 1))
            if k_ok: continue
            try:
                rr, on = lsq(x, dt, nu_ref); p = rr.x
                diag = (f"{rr.status} ({rr.message.strip()}), {rr.nfev} | {p[4] / (2 * math.pi) / nu_ref:.4f}, {p[3]:.4g}, {p[4] * p[3] / 2:.3g}, "
                        f"{p[0]:.3g}, {p[2]:.3g}, {p[1]:.4g}, {p[5]:+.3f} | {', '.join(on) or 'none'}")
            except Exception as e:
                diag = f"raises {type(e).__name__}: {str(e).strip()} | - | -"
            fit = f"raises {err}" if err else f"CONVERGES on the Mac: nu_d {nd:.7g} (KOA: fit_ok 0)"
            print(f"| {run} | {int(koa.loc[run, 'seed'])} | {float(koa.loc[run, 'nu_arg']):.7g} | {fit} | {criterion(err) if err else '-'} | "
                  f"{diag} | {FU.crossings(x, dt, nu_ref)} |")
        print(f"\nthe repetition check: fit_ok equal on the Mac and on KOA for {same} of {len(acf)} trajectories; largest |nu_d (Mac) / nu_d (KOA) - 1| "
              f"over the {len(dev)} fitted on both: {max(dev) if dev else float('nan'):.2e} ({time.time() - t0:.0f} s)\n")
    for cid, Ms in sorted(flagged_cells.items()):
        cell = C[cid]; left, out = [], []
        for M in sorted(cell["masses"]):
            d = os.path.join(base, cid, f"m_{M}")
            r = V.per_mass(cell, M, V.read_mass(d), jackknife=False) if os.path.isdir(d) else None
            if cell["masses"][M] < 1: continue
            if r and not r["flagged"] and math.isfinite(r.get("c", math.nan)): left.append(M)
            else: out.append(M)
        print(f"## the group rule in cell {cid}: flagged {', '.join(map(str, Ms))}; the masses alpha >= 1 left for E1: {len(left)} of "
              f"{len(left) + len(out)} (left: {', '.join(map(str, left)) or 'none'}" + (f"; out: {', '.join(map(str, out))}" if out else "") + "): "
              + ("5 or more -- the cell keeps its E1 value" if len(left) >= 5 else "FEWER THAN 5 -- the cell has no E1 value") + "\n")
    if not flagged_cells: print("no folder is flagged by the group rule, so every cell keeps all of its masses for E1\n")


if __name__ == "__main__":
    main()
