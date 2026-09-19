#!/usr/bin/env python3
"""##CHRIS 2026-09-14: final A1 v2 figures with the Paper 1 primary estimator.
Per trajectory: mean-subtracted FFT of the first 200 predicted oscillations after release, largest bin at
f >= nu_pred / 2.5 (k_min = 80). Mean over health-clean seeds per (eta, M); c_s = through-origin slope over
the nine masses; error bar = 1 sigma scatter of per-mass c_s. No T_i correction (drift-first, T_i = 1).
Writes 260914_A1v2_final_cs_vs_eta.csv and draws 260914_cs_vs_eta / 260914_cs_idealgas_zoom."""
import os, sys, subprocess
import numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import tests_20260913 as T

TD, X_EDGE = 200, 2.5
LEGEND = ("A1 v2, N = 100, 9 masses × 25 seeds, drift-first seeding, 200 oscillations,\n"
          "largest FFT bin at f ≥ ν_pred/2.5, error bar = 1σ scatter over masses")

def cell(task):
    eta, L0, M, runs = task
    nus, nd = [], 0
    for r, p, disc in runs:
        if disc:
            nd += 1; continue
        t, x, nup = T._load(p); dt = (t[-1] - t[0]) / (len(t) - 1)
        n = T._prefix(t, nup, TD); P, df = T._spectrum(x[:n], dt)
        k = int(round(TD / X_EDGE)); nus.append((k + int(np.argmax(P[k:]))) * df)
    return dict(eta=eta, L0=L0, M=M, nu=float(np.mean(nus)), sd=float(np.std(nus, ddof=1)), n=len(nus), nd=nd)

if __name__ == "__main__":
    table = T.a1_leaf_table(); tasks = []
    for leaf in table:
        for M in T.A1_MASSES:
            runs = T.cell_runs(os.path.join(T.DROOT, leaf["leaf"], f"m_{M}"), M)
            if runs: tasks.append((leaf["eta"], float(leaf["L0"]), M, runs))
    with Pool(10) as pool:
        cells = pool.map(cell, tasks, chunksize=1)
    rows = []
    for leaf in table:
        cs = [c for c in cells if c["eta"] == leaf["eta"]]; L0 = float(leaf["L0"])
        c, s, nm = T.slope([T.x_of(q["M"], L0) for q in cs], [q["nu"] for q in cs])
        rows.append(dict(eta=f"{leaf['eta']:.6f}", L0=leaf["L0"], c_s=f"{c:.5f}", c_s_err="0", c_s_scatter_mass=f"{s:.5f}",
                         n_masses=nm, trajectories_used=sum(q["n"] for q in cs), trajectories_discarded=sum(q["nd"] for q in cs),
                         KR=f"{T.kr_cs(leaf['eta']):.5f}" if leaf["eta"] <= 0.69 else "",
                         dev_KR_pct=f"{100*(c/T.kr_cs(leaf['eta'])-1):+.3f}" if leaf["eta"] <= 0.69 else ""))
    csvp = os.path.join(T.PLOTS, "260914_A1v2_final_cs_vs_eta.csv"); T.write_csv(csvp, rows)
    print("| η | L0 | c_s | ± scatter | KR | dev from KR [%] | used / discarded |")
    print("|---|---|---|---|---|---|---|")
    for r in rows:
        print(f"| {r['eta']} | {r['L0']} | {r['c_s']} | {r['c_s_scatter_mass']} | {r['KR'] or 'n/a'} | {r['dev_KR_pct'] or 'n/a'} | {r['trajectories_used']} / {r['trajectories_discarded']} |")
    res = subprocess.run(["/opt/homebrew/bin/python3", os.path.join(T.PLOTS, "plot_cs_meeting.py"), csvp,
                          os.path.join(T.PLOTS, "routeA_fit_input_manifest_20260909.csv"), os.path.join(T.PLOTS, "260914"), "final", LEGEND],
                         cwd=T.PLOTS, capture_output=True, text=True)
    print(res.stdout[-1500:], res.stderr[-800:])
