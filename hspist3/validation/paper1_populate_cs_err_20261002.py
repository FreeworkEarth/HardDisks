#!/usr/bin/env python3
"""##CHRIS 2026-10-02: populate the c_s_err column of the A1 v2 canonical table. Analysis only.

WHAT WAS ACTUALLY WRONG. `c_s_err` has been the literal string "0" since
final_A1_figures_20260914.py:40 wrote it that way. It is a DEAD COLUMN: no plotter reads it. The
canonical figures draw `c_s_scatter_mass` (paper1_canonical_20260919.py sets S = c_s_scatter_mass
and passes yerr=S), and the legend says so in words -- "error bar = 1 sigma scatter over masses".
So no Paper 1 figure was ever drawn without error bars. What they were drawn with, however, is the
SCATTER of the nine per-mass values, which is not the uncertainty on c_s: it is larger by roughly
sqrt(9) = 3, and it never propagated the per-mass frequency errors from the 25 seeds at all.

WHAT THIS SCRIPT ADDS, without touching a single central value:

  c_s_err        the standard error of the THROUGH-ORIGIN slope, per-mass nu errors propagated:
                     s = sum(x y)/sum(x^2)                     (unweighted -- identical to T.slope,
                                                                which is why c_s cannot move)
                     Var(s) = sum(x^2 sigma_y^2)/(sum x^2)^2
                 with sigma_y the standard error of nu over the seeds of that (eta, M) cell.
  chi2_red       sum((y - s x)^2/sigma_y^2)/(n-1). This is the diagnostic that matters: chi2_red >> 1
                 says the nine masses disagree by more than their own seed errors, i.e. there is a
                 real mass dependence or the cell is not equilibrated. It is exactly what separates
                 "small error bar" from "trustworthy small error bar".
  c_s_err_scaled the error bar the FIGURES should carry: c_s_err * max(1, sqrt(chi2_red)), the
                 standard PDG scale factor. chi2_red comes out 1.3-33, i.e. the nine masses
                 disagree by MORE than their seed errors everywhere, so the bare propagated error
                 is a lower bound and plotting it alone would understate the uncertainty -- the
                 exact failure this exercise exists to prevent.
  c_s_scatter_mass  UNCHANGED, still the sd(nu/x) over masses, so nothing that quoted it moves.

A TRAP FOUND ON THE FIRST ATTEMPT, recorded so nobody re-treads it. The historical 260914 table
was computed when T.x_of used l_eff = L0 - 2r WITHOUT the divider thickness, and
paper1_canonical_20260919.py then multiplied by (L0-2r-t/2)/(L0-2r) to produce the canonical 260919
table. T.x_of NOW includes the thickness itself (T.l_eff = L0 - 2r - 0.5t), so recomputing and then
rescaling applies the factor TWICE. The signature is unmistakable: identical nu, identical mass
scatter, and a c_s shift that grows from -42 ppm at eta = 0.0065 to -1314 ppm at eta = 0.196,
matching (L0-2r-t/2)/(L0-2r) at every density. Recomputing with the current T.x_of and NOT
rescaling reproduces the published canonical c_s to better than 5e-6 everywhere. So this script
writes the canonical 260919 table DIRECTLY and does not touch 260914, which stays as the historical
no-thickness record.

##CHRIS 2026-10-14: BOX-TRUNCATION CORRECTION (methods sec. 14; pre-registered 4db8c9d, verdict REGENERATE 615561c,
GO from the plan author). The binary sizes the box in whole pixels, SIM_WIDTH = (int)(2 L0 * 24) (00ALLINONE.c:323),
so a non-grid L0 runs in a box shorter by delta; the divider oscillates about the truncated centre, each compartment
is L0 - delta/2 on average (the runs' own Center_X confirms -delta/2 to 5.5e-6). This script now writes the
CORRECTED canonical table: eta = eta_true = eta_rec L0/(L0 - delta/2), c_s from the per-mass nu with
L_eff,true = L0 - 2r - t/2 - delta/2, KR and dev_KR_pct at eta_true (for the same 24 cells, chosen by eta_rec <= 0.69),
plus columns eta_rec, delta_sigma, L_eff_true. L0 stays the value the runs were launched with. GATE before writing:
the UNcorrected recomputation must reproduce the pre-correction table (PRE19, the dated copy) to 5e-6 in every cell.
"""
import csv, math, os, sys
import numpy as np
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import tests_20260913 as T
import edmd_acc_guard   # ##CHRIS 2026-10-08 (261012 sec. 4.7.4, decision 2): the loader provenance guard (full name: no alias can be shadowed)

TD, X_EDGE = 200, 2.5           # identical to final_A1_figures_20260914.py
RDISK, WALL_T = 0.5, 0.05
SRC14 = "260914_A1v2_final_cs_vs_eta.csv"
SRC19 = "260919_A1v2_final_cs_vs_eta.csv"
PRE19 = "260919_A1v2_final_cs_vs_eta_pre_boxtrunc_20261014.csv"   # ##CHRIS 2026-10-14: the uncorrected table


def box_delta(L0):
    """##CHRIS 2026-10-14: box shortfall [sigma] of SIM_WIDTH = (int)(2 * L0_UNITS * PIXELS_PER_SIGMA) (00ALLINONE.c:323),
    evaluated in the binary's own float32 arithmetic."""
    w = np.float32(2) * np.float32(L0) * np.float32(24)
    return (float(w) - math.floor(float(w))) / 24.0


def cell(task):
    """Per-(eta, M) frequency, byte-for-byte the estimator of final_A1_figures_20260914.cell()."""
    eta, L0, M, runs = task
    for r, p, disc in runs:
        edmd_acc_guard.guard(p)          # every trajectory of the cell, also the discarded ones
    nus, nd = [], 0
    for r, p, disc in runs:
        if disc:
            nd += 1; continue
        t, x, nup = T._load(p); dt = (t[-1] - t[0]) / (len(t) - 1)
        n = T._prefix(t, nup, TD); P, df = T._spectrum(x[:n], dt)
        k = int(round(TD / X_EDGE)); nus.append((k + int(np.argmax(P[k:]))) * df)
    if not nus:
        return dict(eta=eta, L0=L0, M=M, nu=float("nan"), sd=float("nan"), n=0, nd=nd)
    return dict(eta=eta, L0=L0, M=M, nu=float(np.mean(nus)),
                sd=float(np.std(nus, ddof=1)) if len(nus) > 1 else float("nan"),
                n=len(nus), nd=nd)



def slope_with_errors(x, y, sy):
    """Through-origin slope with propagated + chi2-scaled error. ONE definition, imported by the
    A2 figure scripts too, so the A1 and A2 error bars cannot drift apart.

    Returns (s, err, err_scaled, chi2_red). The slope is UNWEIGHTED, matching T.slope and matching
    what every existing figure already plots as its central value, so adopting this cannot move a
    published number -- only the error bar changes.
    """
    import numpy as _np, math as _math
    x = _np.asarray(x, float); y = _np.asarray(y, float); sy = _np.asarray(sy, float)
    ok = _np.isfinite(x) & _np.isfinite(y) & _np.isfinite(sy) & (sy > 0)
    if ok.sum() < 3:
        return float("nan"), float("nan"), float("nan"), float("nan")
    x, y, sy = x[ok], y[ok], sy[ok]
    sxx = float((x * x).sum())
    s = float((x * y).sum() / sxx)
    err = _math.sqrt(float((x * x * sy * sy).sum()) / sxx ** 2)
    chi2 = float((((y - s * x) / sy) ** 2).sum()) / max(1, len(x) - 1)
    return s, err, err * max(1.0, _math.sqrt(chi2)), chi2

def thickness_factor(L0):
    lo = L0 - 2 * RDISK
    return (lo - 0.5 * WALL_T) / lo


def main():
    table = T.a1_leaf_table()
    tasks = []
    for leaf in table:
        for M in T.A1_MASSES:
            runs = T.cell_runs(os.path.join(T.DROOT, leaf["leaf"], f"m_{M}"), M)
            if runs:
                tasks.append((leaf["eta"], float(leaf["L0"]), M, runs))
    print(f"{len(tasks)} (eta, M) cells; ~{sum(len(t[3]) for t in tasks)} trajectories")
    # Pool(6), not 10: the Md10 long campaign is running on this machine tonight.
    with Pool(6) as pool:
        cells = pool.map(cell, tasks, chunksize=1)

    old = {r["eta"]: r for r in csv.DictReader(open(T.plot_path(SRC14)))}
    pre = {r["eta"]: r for r in csv.DictReader(open(T.plot_path(PRE19)))}    # ##CHRIS 2026-10-14
    rows, moved, gate_fail = [], [], []
    for leaf in table:
        cs = [c for c in cells if c["eta"] == leaf["eta"] and c["n"] > 0]
        if len(cs) < 3:
            continue
        L0 = float(leaf["L0"])
        x = np.array([T.x_of(q["M"], L0) for q in cs])
        y = np.array([q["nu"] for q in cs])
        sy = np.array([(q["sd"] / math.sqrt(q["n"])) if q["n"] > 1 else np.nan for q in cs])

        s, scat, nm = T.slope(x, y)                    # unchanged central value and scatter
        sxx = float((x * x).sum())
        var = float((x * x * sy * sy).sum()) / (sxx ** 2)
        err = math.sqrt(var)
        resid = y - s * x
        chi2 = float(((resid / sy) ** 2).sum()) / max(1, len(x) - 1)

        key = f"{leaf['eta']:.6f}"
        o = old.get(key)
        if o is not None and abs(float(o["c_s"]) - s) > 5e-6 * max(1.0, abs(s)):
            moved.append((key, float(o["c_s"]), s))
        # ##CHRIS 2026-10-14 GATE: the uncorrected recomputation must reproduce the pre-correction canonical table.
        p = pre.get(key)
        if p is None or abs(float(p["c_s"]) - s) > 5e-6 * max(1.0, abs(s)):
            gate_fail.append((key, p["c_s"] if p else None, s))
        # ##CHRIS 2026-10-14 BOX-TRUNCATION CORRECTION (methods sec. 14): same per-mass nu, true length and density.
        delta = box_delta(L0)
        LeT = T.l_eff(L0) - delta / 2
        xT = np.array([T.k_root(q["M"] / (2.0 * T.N_SIDE)) / (2 * math.pi * LeT) for q in cs])
        sT, scatT, _ = T.slope(xT, y)
        sxxT = float((xT * xT).sum()); errT = math.sqrt(float((xT * xT * sy * sy).sum()) / (sxxT ** 2))
        chi2T = float((((y - sT * xT) / sy) ** 2).sum()) / max(1, len(xT) - 1)
        etaT = leaf["eta"] * L0 / (L0 - delta / 2)
        r = dict(o) if o else {}
        r.update(eta=f"{etaT:.6f}", L0=leaf["L0"], c_s=f"{sT:.5f}", c_s_err=f"{errT:.6f}",
                 c_s_scatter_mass=f"{scatT:.5f}", chi2_red=f"{chi2T:.3f}",
                 c_s_err_scaled=f"{errT * max(1.0, math.sqrt(chi2T)):.6f}", n_masses=nm,
                 trajectories_used=sum(q["n"] for q in cs),
                 trajectories_discarded=sum(q["nd"] for q in cs),
                 eta_rec=key, delta_sigma=f"{delta:.6f}", L_eff_true=f"{LeT:.6f}")
        if leaf["eta"] <= 0.69:                      # the KR set is chosen by eta_rec: the same 24 cells as before
            r["KR"] = f"{T.kr_cs(etaT):.5f}"
            r["dev_KR_pct"] = f"{100 * (sT / T.kr_cs(etaT) - 1):+.3f}"
        else:
            r["KR"] = ""; r["dev_KR_pct"] = ""
        rows.append(r)

    fields = ["eta", "L0", "c_s", "c_s_err", "c_s_err_scaled", "c_s_scatter_mass", "chi2_red", "n_masses",
              "trajectories_used", "trajectories_discarded", "KR", "dev_KR_pct",
              "eta_rec", "delta_sigma", "L_eff_true"]          # ##CHRIS 2026-10-14: box-truncation columns
    if gate_fail:                                                   # ##CHRIS 2026-10-14: fail closed
        print("GATE FAILED -- the uncorrected recomputation does not reproduce the pre-correction table; NOTHING WRITTEN:")
        for g in gate_fail: print("   ", g)
        sys.exit(1)
    print(f"gate: the uncorrected recomputation reproduces {PRE19} to 5e-6 in all {len(rows)} cells")
    # Write the CANONICAL table directly. No rescale: see the trap in the module docstring.
    # 260914 is left exactly as it is -- it is the historical no-thickness record.
    out19 = os.path.join(T.PLOTS, SRC19)
    old19 = {r["eta"]: r for r in csv.DictReader(open(T.plot_path(SRC19)))}
    moved19 = []
    for r in rows:
        o = old19.get(r.get("eta_rec", r["eta"]))          # ##CHRIS 2026-10-14: match on the recorded eta
        if o is not None and abs(float(o["c_s"]) - float(r["c_s"])) > 5e-6 * max(1.0, abs(float(r["c_s"]))):
            moved19.append((r["eta"], float(o["c_s"]), float(r["c_s"])))
    with open(out19, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    print(f"wrote {out19}")

    print("\n=== CENTRAL VALUES: did anything move? ===")
    print(f"  canonical 260919: {len(moved19)} of {len(rows)} rows changed beyond 5e-6 relative")
    for k, a, b in moved19[:10]:
        print(f"     eta {k}: {a:.6f} -> {b:.6f}   ({1e6*(b/a-1):+.1f} ppm)")
    if not moved19:
        print("     none -- every published c_s reproduced to better than 5e-6")

    print("\n=== the two error bars, and the diagnostic ===")
    print("|   eta |      c_s | c_s_err | chi2_red | c_s_err_scaled (FOR FIGURES) | scatter_mass (old) |")
    print("|---|---|---|---|---|---|")
    for r in rows:
        e, c = float(r["eta"]), float(r["c_s"])
        er, sc, x2 = float(r["c_s_err"]), float(r["c_s_scatter_mass"]), float(r["chi2_red"])
        es = float(r["c_s_err_scaled"])
        print(f"| {e:.4f} | {c:8.4f} | {er:7.5f} | {x2:8.2f} | {es:28.5f} | {sc:18.5f} |")
    import numpy as _np
    R = _np.array([float(r["c_s_err_scaled"]) / float(r["c_s_scatter_mass"]) for r in rows])
    print(f"\nscaled-error / old-plotted-scatter: median {_np.median(R):.2f}, range {R.min():.2f}-{R.max():.2f}")
    print("i.e. the error bars the canonical figures have been carrying are the right SIZE;")
    print("what was wrong was the justification, not the picture.")


if __name__ == "__main__":
    main()
