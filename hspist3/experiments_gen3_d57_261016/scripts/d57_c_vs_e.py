#!/usr/bin/env python3
"""##CHRIS 2026-10-16 (plan-author decision 57, addendum 3, task 1; notes 261012 sec. 4.7.140): c AGAINST THE RELAXED STRAIN e, over the
8 P2 cells, plus D56-1's eta 0.76 pilot. No run; no model beyond a straight line.
  e per run: the mean over its record's snapshots (p2_frames.csv, phase 1) of the strain of both compartments (p2_lattice.config_strain,
     left and right averaged); per cell: the mean over its 50 runs, its SE (sd / sqrt(n)) and sd; the end strain (the last snapshot) and the
     strain at the hold's end alongside.
  c per cell: E1 over the 5 masses with its SE (Report D57-2, experiments_gen3_p2_261016/report_261016/p2_cells_e1.csv).
  The pilot (experiments_gen3_d56_261015, N 900, eta 0.76, H 30, the default orientation = turned 17 x 27, alpha 5 only, 2 runs): e from its
     snapshot files the same way; c = its E1 of the 2 runs (c_analysis_b/c_cells.csv; no jackknife SE: the error bar is half the difference
     of its two runs' own c, indicative only).
  The slopes d ln c / d e (e as a fraction, so the slope is the % change of c per % of e): a weighted straight line ln c = a + b e with the
     effective variance s_i^2 = (SE_c / c)^2 + b^2 SE_e^2 (iterated); the SE of b from the fit, and scaled by sqrt(chi2/dof) when that
     exceeds 1. eta 0.72: the three boxes (H 30: the two G1 cells, the same box and both turned at the end, combined as the inverse-variance
     mean of c and the mean of e; H 29.833: G2 turned; H 29.708: G2 pin 0, turned at the end) and, as a check, the four cells; eta 0.78:
     the two turned cells (two points: the slope of the line through them, its SE propagated from both points' SE in c and e); and, since
     the three 0.72 boxes do not lie on one line within their SE, the two 0.72 segments (two points each) alongside; and a sensitivity row:
     the three boxes with point 1's e error set to half the two G1 cells' e difference (the two cells share the box but differ in e by
     about 7.7 SE; found by the independent verification of Report D57-3).
Writes c_vs_e/: c_vs_e_cells.csv, c_vs_e_runs.csv, c_vs_e_slopes.csv, c_vs_e.md, fig_d57_c_vs_e.png.
usage (from hspist3/): python3 experiments_gen3_d57_261016/scripts/d57_c_vs_e.py [--out <a new folder>]
"""
import glob, math, os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); E = os.path.join(HS, "experiments_gen3_d57_261016")
sys.path.insert(0, os.path.join(HS, "cluster", "gen3_stage2_261016"))
import p2_lattice as LT
P2 = os.path.join(HS, "experiments_gen3_stage2_koa_261016", "P2", "merged", "stage2", "gen3")
E1F = os.path.join(HS, "experiments_gen3_p2_261016", "report_261016", "p2_cells_e1.csv")
PIL = os.path.join(HS, "experiments_gen3_d56_261015", "data", "pilot", "P2_N0900_eta0.760000", "m_4500")
PILC = os.path.join(HS, "experiments_gen3_d56_261015", "c_analysis_b")
OUT = os.path.abspath(sys.argv[sys.argv.index("--out") + 1]) if "--out" in sys.argv else os.path.join(E, "c_vs_e")
C_PIN = {"0": "#2a6fbb", "turned": "#d9822b"}; MK_G = {"G1": "o", "G2": "s"}


def runs_table():
    rows = []
    for d in sorted(glob.glob(os.path.join(P2, "*", "m_*"))):
        cid = d.split(os.sep)[-2]; M = int(d.split("_")[-1]); R = pd.read_csv(os.path.join(d, "p2_runs.csv")); F = pd.read_csv(os.path.join(d, "p2_frames.csv"))
        F = F[F["phase"] == 1]
        for r in R.to_dict("records"):
            f = F[F["run"] == r["run"]]
            rows.append(dict(cell=cid, eta=r["eta"], G=r["G"], pin=str(r["pin"]), H=r["H"], L0=r["L0"], M=M, run=int(r["run"]), seed=int(r["seed"]),
                             e_rec=float(((f["e_0"] + f["e_1"]) / 2).mean()), e_rec_0=float(f["e_0"].mean()), e_rec_1=float(f["e_1"].mean()), frames=len(f),
                             e_end=(r["end_e_0"] + r["end_e_1"]) / 2, e_hold=(r["hold_e_0"] + r["hold_e_1"]) / 2, line_end_0=r["end_line_0"], line_end_1=r["end_line_1"],
                             psi6=(r["psi6_rec_0"] + r["psi6_rec_1"]) / 2, Z=r["Z"]))
    return pd.DataFrame(rows)


def pilot():
    rows = []
    for p in sorted(glob.glob(os.path.join(PIL, "snapshots_*_run*.bin*"))):
        r = int(p.split("_run")[1].split(".")[0]); S = LT.read_snapshots(p)
        ee = [LT.config_strain(fr["xy"], fr["cen"], S["W"], S["H"]) for fr in S["frames"] if fr["phase"] == 1]
        e0, e1 = np.array([x[0]["e"] for x in ee]), np.array([x[1]["e"] for x in ee])
        rows.append(dict(run=r, frames=len(ee), e_rec=float(np.nanmean((e0 + e1) / 2)), e_rec_0=float(np.nanmean(e0)), e_rec_1=float(np.nanmean(e1)),
                         line_0=ee[-1][0]["line"], line_1=ee[-1][1]["line"]))
    R = pd.DataFrame(rows); cr = pd.read_csv(os.path.join(PILC, "c_runs.csv")); cr = cr[cr["cell"] == 76]; cc = pd.read_csv(os.path.join(PILC, "c_cells.csv"))
    cc = cc[cc["cell"] == 76].iloc[0]; x = float(cc["x"])
    cs = [math.sqrt(a ** 2 + (1 / (2 * math.pi * t)) ** 2) / x for a, t in zip(cr["nu_d"], cr["tau_r"])]
    return R, dict(c=float(cc["c"]), sig=abs(cs[0] - cs[1]) / 2 if len(cs) == 2 else math.nan, c_runs=cs, n=int(cc["n"]),
                   psi6=float(((cr["psi6_rec_0"] + cr["psi6_rec_1"]) / 2).mean()), Z=float(cr["Z"].mean()))


def fit_line(e, se, c, sc):
    """weighted straight line ln c = a + b e, effective variance; returns b, SE_b, SE_b scaled, chi2, dof, a"""
    e, se, y, sy = map(np.asarray, (e, se, np.log(c), np.asarray(sc) / np.asarray(c)))
    if len(e) == 2:
        de = e[1] - e[0]; b = (y[1] - y[0]) / de
        sb = math.sqrt((sy[0] ** 2 + sy[1] ** 2) / de ** 2 + b ** 2 * (se[0] ** 2 + se[1] ** 2) / de ** 2)
        return dict(b=b, sb=sb, sb_scaled=sb, chi2=0.0, dof=0, a=y[0] - b * e[0])
    b = np.polyfit(e, y, 1)[0]
    for _ in range(50):
        w = 1 / (sy ** 2 + b ** 2 * se ** 2); S, Sx, Sy, Sxx, Sxy = w.sum(), (w * e).sum(), (w * y).sum(), (w * e * e).sum(), (w * e * y).sum()
        D = S * Sxx - Sx ** 2; bn = (S * Sxy - Sx * Sy) / D; a = (Sxx * Sy - Sx * Sxy) / D
        if abs(bn - b) < 1e-12: b = bn; break
        b = bn
    sb = math.sqrt(S / D); chi2 = float((w * (y - a - b * e) ** 2).sum()); dof = len(e) - 2
    return dict(b=b, sb=sb, sb_scaled=sb * max(1.0, math.sqrt(chi2 / dof)) if dof > 0 else sb, chi2=chi2, dof=dof, a=a)


def main():
    if os.path.exists(OUT): sys.exit(f"REFUSED: {OUT} exists (a new folder only)")
    R = runs_table(); E1 = pd.read_csv(E1F); cells = []
    for cid, g in R.groupby("cell", sort=False):
        e1 = E1[E1["cell"] == cid].iloc[0]; n = len(g)
        cells.append(dict(cell=cid, eta=g["eta"].iloc[0], G=g["G"].iloc[0], pin=g["pin"].iloc[0], H=g["H"].iloc[0], L0=g["L0"].iloc[0], runs=n,
                          e_rec=g["e_rec"].mean(), e_rec_se=g["e_rec"].std(ddof=1) / math.sqrt(n), e_rec_sd=g["e_rec"].std(ddof=1),
                          e_rec_left=g["e_rec_0"].mean(), e_rec_right=g["e_rec_1"].mean(), e_end=g["e_end"].mean(), e_hold=g["e_hold"].mean(),
                          lines_y_end=int((g[["line_end_0", "line_end_1"]] == 90).to_numpy().sum()), psi6=g["psi6"].mean(), Z=g["Z"].mean(),
                          c=e1["E1_c"], sig_c=e1["E1_sig"]))
    C = pd.DataFrame(cells)
    PR, pc = pilot()
    pil = dict(cell="D56-1 pilot N0900 eta0.76 (H 30, turned 17x27, alpha 5 only)", eta=0.76, G="pilot", pin="turned", H=30.0, L0=15.501279540739111,
               runs=len(PR), e_rec=PR["e_rec"].mean(), e_rec_se=PR["e_rec"].std(ddof=1) / math.sqrt(len(PR)), e_rec_sd=PR["e_rec"].std(ddof=1),
               e_rec_left=PR["e_rec_0"].mean(), e_rec_right=PR["e_rec_1"].mean(), e_end=math.nan, e_hold=math.nan,
               lines_y_end=int(((PR[["line_0", "line_1"]]) == 90).to_numpy().sum()), psi6=pc["psi6"], Z=pc["Z"], c=pc["c"], sig_c=pc["sig"])
    CA = pd.concat([C, pd.DataFrame([pil])], ignore_index=True)
    k72 = C[np.isclose(C["eta"], 0.72)]; g1 = k72[k72["G"] == "G1"]; w = 1 / g1["sig_c"] ** 2
    comb = dict(e=g1["e_rec"].mean(), se=math.sqrt((g1["e_rec_se"] ** 2).sum()) / 2, c=float((w * g1["c"]).sum() / w.sum()), sc=float(w.sum() ** -0.5))
    p72 = sorted([comb] + [dict(e=r["e_rec"], se=r["e_rec_se"], c=r["c"], sc=r["sig_c"]) for r in k72[k72["G"] == "G2"].to_dict("records")],
                 key=lambda q: q["e"])                       # by e: H 30 (combined), 29.833 (G2 turned), 29.708 (G2 pin 0)
    S = []
    for lab, pts in (("eta 0.72, the three boxes (H 30 combined, 29.833, 29.708)", p72),
                     ("eta 0.72, the four cells (check)", [dict(e=r["e_rec"], se=r["e_rec_se"], c=r["c"], sc=r["sig_c"]) for r in k72.to_dict("records")]),
                     ("eta 0.72, the segment H 30 combined -> 29.833", p72[:2]), ("eta 0.72, the segment 29.833 -> 29.708", p72[1:]),
                     ("eta 0.72, the three boxes, point 1's SE_e = half the two G1 cells' e difference (sensitivity)",
                      [dict(p72[0], se=abs(g1["e_rec"].iloc[0] - g1["e_rec"].iloc[1]) / 2)] + p72[1:]),
                     ("eta 0.78, the two turned cells (H 30, 28.667)", [dict(e=r["e_rec"], se=r["e_rec_se"], c=r["c"], sc=r["sig_c"])
                                                                       for r in C[np.isclose(C["eta"], 0.78) & (C["pin"] == "turned")].to_dict("records")])):
        f = fit_line([p["e"] for p in pts], [p["se"] for p in pts], [p["c"] for p in pts], [p["sc"] for p in pts])
        S.append(dict(fit=lab, points=len(pts), e_min=min(p["e"] for p in pts), e_max=max(p["e"] for p in pts), slope=f["b"], slope_se=f["sb"],
                      slope_se_scaled=f["sb_scaled"], chi2=f["chi2"], dof=f["dof"], ln_c_at_e0=f["a"], c_at_e0=math.exp(f["a"])))
    S = pd.DataFrame(S); os.makedirs(OUT)
    R.to_csv(os.path.join(OUT, "c_vs_e_runs.csv"), index=False, float_format="%.8g"); CA.to_csv(os.path.join(OUT, "c_vs_e_cells.csv"), index=False, float_format="%.8g")
    S.to_csv(os.path.join(OUT, "c_vs_e_slopes.csv"), index=False, float_format="%.8g")
    L = ["# c against the relaxed strain e (decision 57 addendum 3, task 1; d57_c_vs_e.py)", "",
         "e: the mean over the record's snapshots of both compartments' strain, then over the runs (SE = sd / sqrt(n)); strains in %. "
         "c: E1 over the 5 masses with its SE (Report D57-2). The pilot: alpha 5 only, 2 runs; its error bar is half the difference of its two runs' c.", "",
         "| eta | G | pin (start) | H [sigma] | runs | e (record) ± SE | sd over runs | e left / right | e at the hold's end | e at the end | lines along y at the end | psi6 (record) | Z | c ± SE |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in CA.sort_values(["eta", "e_rec"]).to_dict("records"):
        L.append(f"| {r['eta']:.2f} | {r['G']} | {r['pin']} | {r['H']:.4f} | {r['runs']} | {100 * r['e_rec']:+.3f} ± {100 * r['e_rec_se']:.3f} | {100 * r['e_rec_sd']:.3f} | "
                 f"{100 * r['e_rec_left']:+.3f} / {100 * r['e_rec_right']:+.3f} | {('%+.3f' % (100 * r['e_hold'])) if np.isfinite(r['e_hold']) else '-'} | "
                 f"{('%+.3f' % (100 * r['e_end'])) if np.isfinite(r['e_end']) else '-'} | {r['lines_y_end']} of {2 * r['runs']} | {r['psi6']:.3f} | {r['Z']:.2f} | "
                 f"{r['c']:.3f} ± {r['sig_c']:.3f}{' (indicative)' if r['G'] == 'pilot' else ''} |")
    L += ["", "The slopes d ln c / d e (the % change of c per % of e; a weighted straight line in ln c, effective variance; the scaled SE multiplies by "
              "sqrt(chi2/dof) when that exceeds 1):", "",
          "| the points | n | e range [%] | d ln c / d e ± SE | SE scaled | chi2 / dof | c at e = 0 (the line) |", "|---|---|---|---|---|---|---|"]
    for r in S.to_dict("records"):
        L.append(f"| {r['fit']} | {r['points']} | {100 * r['e_min']:+.3f} to {100 * r['e_max']:+.3f} | {r['slope']:.3f} ± {r['slope_se']:.3f} | {r['slope_se_scaled']:.3f} | "
                 f"{r['chi2']:.2f} / {r['dof']} | {r['c_at_e0']:.3f} |")
    gd = g1["e_rec"].iloc[0] - g1["e_rec"].iloc[1]; gs = math.sqrt((g1["e_rec_se"] ** 2).sum())
    cd = g1["c"].iloc[0] - g1["c"].iloc[1]; cs_ = math.sqrt((g1["sig_c"] ** 2).sum())
    L += ["", f"- the two G1 cells at 0.72 (the same box, H 30): e differs by {100 * gd:+.3f} ± {100 * gs:.3f} % ({gd / gs:+.1f} SE) while c differs by "
              f"{cd:+.3f} ± {cs_:.3f}; the sensitivity row of the slopes takes half that e difference as point 1's error"]
    L += [f"- the pilot: its two runs' own c {pc['c_runs'][0]:.3f} and {pc['c_runs'][1]:.3f}; record snapshots per run {', '.join(str(v) for v in PR['frames'])}; "
              f"lines along y at its last record snapshot in {int((PR[['line_0', 'line_1']] == 90).to_numpy().sum())} of {2 * len(PR)} compartments"]
    open(os.path.join(OUT, "c_vs_e.md"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
    fig, axs = plt.subplots(1, 3, figsize=(17, 5.4), gridspec_kw=dict(width_ratios=[1.25, 0.6, 1.25]))
    for ax, eta in zip(axs, (0.72, 0.76, 0.78)):
        g = CA[np.isclose(CA["eta"], eta)]
        for r in g.to_dict("records"):
            col = C_PIN[r["pin"]]; mk = MK_G.get(r["G"], "D"); vac = eta == 0.78 and r["pin"] == "0"
            lab = (f"{r['G']}, start pin {r['pin']}" + (", ended turned" if eta == 0.72 and r["pin"] == "0" else "") + (", 12 vacancies" if vac else "")
                   + f" (H {r['H']:.3f})") if r["G"] != "pilot" else "D56-1 pilot (alpha 5, 2 runs; indicative error)"
            ax.errorbar(100 * r["e_rec"], r["c"], xerr=100 * r["e_rec_se"], yerr=r["sig_c"], marker=mk, ms=8, color=col, mfc=col if not vac else "white",
                        mec=col, capsize=3, lw=1.2, ls="none", label=lab)
        for f in S.to_dict("records"):
            if (eta == 0.72 and f["fit"].startswith("eta 0.72, the three")) or (eta == 0.78 and f["fit"].startswith("eta 0.78")):
                xs = np.linspace(min(-0.4, 100 * f["e_min"] - 0.3), 100 * f["e_max"] + 0.3, 50)
                ax.plot(xs, np.exp(f["ln_c_at_e0"] + f["slope"] * xs / 100), color="0.35", lw=1.0, ls="--",
                        label=f"the straight line in ln c: d ln c / d e = {f['slope']:.2f} ± {f['slope_se_scaled']:.2f}")
        ax.axvline(0, color="0.6", lw=0.7, ls=":"); ax.grid(alpha=0.25); ax.set_xlabel("relaxed strain e, the record's mean  [%]", fontsize=9)
        ax.set_ylabel("c  (E1 over the masses)", fontsize=9); ax.set_title(f"eta {eta}" + (": D56-1's pilot only" if eta == 0.76 else ""), fontsize=10)
        ax.legend(fontsize=7, loc="best", framealpha=0.9); ax.tick_params(labelsize=8)
    fig.suptitle("c against the crystal's relaxed strain e (P2, N 900; colour: the start's pin; circles G1, squares G2; open: 12 vacancies)", fontsize=10.5)
    fig.savefig(os.path.join(OUT, "fig_d57_c_vs_e.png"), dpi=120, bbox_inches="tight"); plt.close(fig)


if __name__ == "__main__":
    main()
