#!/usr/bin/env python3
"""##CHRIS 2026-10-11 (plan-author decision 29, task P; 261012 sec. 4.7.59): THE FIGURE AND THE TABLE OF THE PREVIEW SWEEP -- PREVIEW,
INFORMATION; no verdict, never part of stage 1, never in a paper. From two tables only:
  - the preview's E0 table: e0_preview_261011.py --preview's CSV, one row per cell in the canonical table's columns;
  - the canonical A1 v2 table 260919_A1v2_final_cs_vs_eta.csv (engine 279282b; methods sec. 14.2), read unchanged.
FIGURE (the project's plot rules; the shared helpers of plot_speed_of_sound_edmd.py):
  (a) the preview's c_s(eta), blue, with error bars c_s_err_scaled (the error bar the canonical figures carry), at the box-corrected eta;
      KR (red), SPT and Henderson, all with the adiabatic mapping (speed_of_sound_theory_curves, draw_theory_curves); the eta regimes
      shaded (add_eta_regime_shading);
  (b) the canonical A1 v2 values on the same axes;
  (c), (d) below each: the residual to KR, 100 (c_s / KR(eta) - 1), with its error 100 c_s_err_scaled / KR. KR is computed here with
      the current module for both, at each table's own eta and only where eta_rec <= 0.69. The A1 v2 table's own KR column differs
      from it by at most 6e-6 relative (e0_preview_261011.py --gate-a1v2).
  "PREVIEW" stands in the title and in panels (a) and (c); N, radius, L0 and H in the footer; one build generation per panel (gen3
  a643bc7 | gen2 279282b).
TABLE: preview minus A1 v2 per density, wherever the two tables have the same cell (the same L0 and recorded eta): the difference of
c_s with its error sqrt(err_p^2 + err_a^2), from c_s_err_scaled and from c_s_err, and z. A preview cell without such an A1 v2 cell is
listed with the nearest A1 v2 cell, and no difference is formed.
usage (from hspist3/): python3 cluster/gen3_preview_261011/fig_preview_261011.py --preview <E0 csv> --out <folder> [--label TEXT]
  --label   replaces "PREVIEW" in the title (the Mac rehearsal's dry test says what it is)
"""
import argparse, csv, math, os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(HS, "validation")); sys.path.insert(0, HS)
import tests_20260913 as T
import plot_speed_of_sound_edmd as S

X_MAX = 0.65


def resid(df):
    """100 (c_s / KR - 1) and its error, KR from the current module at the row's own eta, where eta_rec <= 0.69"""
    out = []
    for _, r in df.iterrows():
        if float(r["eta_rec"]) <= 0.69:
            kr = T.kr_cs(float(r["eta"])); out.append((100 * (float(r["c_s"]) / kr - 1), 100 * float(r["c_s_err_scaled"]) / kr))
        else: out.append((math.nan, math.nan))
    return np.array(out)


def diff_table(P, A):
    rows = []
    for _, p in P.iterrows():
        m = A[(A["eta_rec"].map(lambda v: f"{float(v):.6f}") == f"{float(p['eta_rec']):.6f}") & (A["L0"].astype(str) == str(p["L0"]))]
        if len(m) == 1:
            a = m.iloc[0]; d = float(p["c_s"]) - float(a["c_s"])
            es = math.hypot(float(p["c_s_err_scaled"]), float(a["c_s_err_scaled"])); e = math.hypot(float(p["c_s_err"]), float(a["c_s_err"]))
            rows.append(dict(cell=p.get("cell", ""), eta=p["eta"], eta_rec=p["eta_rec"], L0=p["L0"], c_s_preview=p["c_s"], c_s_A1v2=a["c_s"],
                             difference=f"{d:+.5f}", rel_pct=f"{100 * d / float(a['c_s']):+.3f}", err_scaled=f"{es:.5f}",
                             z_scaled=f"{d / es:+.2f}", err=f"{e:.5f}", z=f"{d / e:+.2f}", note="the same cell"))
        else:
            k = (A["eta_rec"].astype(float) - float(p["eta_rec"])).abs().idxmin(); a = A.loc[k]
            rows.append(dict(cell=p.get("cell", ""), eta=p["eta"], eta_rec=p["eta_rec"], L0=p["L0"], c_s_preview=p["c_s"], c_s_A1v2="",
                             difference="", rel_pct="", err_scaled="", z_scaled="", err="", z="",
                             note=f"no A1 v2 cell here; the nearest is eta_rec {float(a['eta_rec']):.6f} (L0 {a['L0']})"))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--preview", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--label", default="PREVIEW"); a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    P = pd.read_csv(a.preview, dtype={"L0": str}).sort_values("eta")
    A = pd.read_csv(T.plot_path("260919_A1v2_final_cs_vs_eta.csv"), dtype={"L0": str})
    A = A[A["eta"].astype(float) <= X_MAX].sort_values("eta")
    grid = np.linspace(0.001, X_MAX, 400)
    th = S.speed_of_sound_theory_curves(grid, mapping="adiabatic", kbt=1.0, m=1.0)
    th = {k: th[k] for k in ("ideal", "spt", "henderson", "kolafa_rottner")}
    fig, ax = plt.subplots(2, 2, figsize=(14, 9.6), sharex=True, gridspec_kw=dict(height_ratios=[2.2, 1], hspace=0.12, wspace=0.12,
                                                                                    left=0.06, right=0.99, top=0.865, bottom=0.12))
    for col, (D, lab, who) in enumerate(((P, f"{a.label}: gen3, engine-gen3 a643bc7 on KOA", "preview"),
                                         (A, "canonical A1 v2: gen2, engine 279282b (methods sec. 14.2)", "A1 v2"))):
        top, bot = ax[0, col], ax[1, col]
        S.add_eta_regime_shading(top, x_max=X_MAX, x_min=0.0); S.add_eta_regime_shading(bot, x_max=X_MAX, x_min=0.0)
        S.draw_theory_curves(top, grid, th, theory_label="adiabatic")
        S.draw_measured_series(top, D["eta"].astype(float), D["c_s"].astype(float), D["c_s_err_scaled"].astype(float),
                               label=f"{who}: c_s by E0 (error bar c_s_err_scaled)", color="C0")
        top.set_title(f"({'ab'[col]}) {lab}", loc="left", fontsize=10)
        top.set_ylabel("c_s [reduced units]"); top.legend(loc="center left", fontsize=7.5, framealpha=0.9)
        R = resid(D)
        bot.axhline(0, color="C3", lw=1.6, alpha=0.75)
        bot.errorbar(D["eta"].astype(float), R[:, 0], yerr=R[:, 1], fmt="o", color="C0", ms=4, capsize=3, lw=1)
        bot.set_ylabel("c_s / KR - 1 [%]"); bot.set_xlabel("η (box-corrected, methods sec. 14)")
        bot.set_title(f"({'cd'[col]}) residual to Kolafa–Rottner 2006 (η_rec ≤ 0.69)", loc="left", fontsize=9)
    ax[0, 1].sharey(ax[0, 0]); ax[1, 1].sharey(ax[1, 0]); ax[0, 0].set_xlim(0, X_MAX)
    for b in (ax[1, 0], ax[1, 1]):
        v = np.concatenate([resid(P)[:, 0], resid(A)[:, 0]]); e = np.concatenate([resid(P)[:, 1], resid(A)[:, 1]])
        ok = np.isfinite(v); b.set_ylim(min(-1.0, np.min(v[ok] - e[ok]) - 0.5), max(1.0, np.max(v[ok] + e[ok]) + 0.5))
    for x in (ax[0, 0], ax[1, 0]):
        x.text(0.99, 0.04, a.label, transform=x.transAxes, ha="right", va="bottom", fontsize=13, fontweight="bold", color="#555555", alpha=0.6)
    fig.suptitle(f"{a.label}: c_s(η) at N = 100 by the canonical estimator E0\n(INFORMATION, no verdict; never part of stage 1, "
                 f"never in a paper)", x=0.01, ha="left", fontsize=12, fontweight="bold", y=0.995)
    fig.text(0.01, 0.9, "N = 100 disks (50 per compartment), radius 0.5, H = 10, wall thickness 0.05; L0 per cell from its density (the "
             "preview: " + ", ".join(f"{float(l):g}" for l in P["L0"]) + ");\nthe nine A1 v2 masses, 25 seeds each; the canonical A1 v2 "
             "command (methods sec. 8)", fontsize=8.5, va="bottom")
    fig.text(0.01, 0.01, "E0: per (eta, M) the mean argmax frequency over the seeds; per cell the through-origin slope over the nine masses "
             "against x_M = K(alpha)/(2 pi L_eff,true), with the box-truncation correction of methods sec. 14 (eta = eta_rec L0/(L0 - delta/2), "
             "L_eff,true = L0 - 2r - t/2 - delta/2);\nerror bar c_s_err_scaled = c_s_err max(1, sqrt(chi2_red)) (PDG scale factor). Theory: "
             "the adiabatic mapping c_s^2 = (kT/m)(Z + eta Z' + Z^2) for every EOS; KR compared only for eta <= 0.69. Residuals: KR of the "
             "current module at each table's own eta.", fontsize=7, color="#555555")
    paths = []
    for ext in ("pdf", "png"):
        p = os.path.join(a.out, f"preview_cs_vs_eta.{ext}")
        fig.savefig(p, dpi=170, metadata=({"CreationDate": None} if ext == "pdf" else {"Software": None})); paths.append(p)
    plt.close(fig)
    Dt = diff_table(P, pd.read_csv(T.plot_path("260919_A1v2_final_cs_vs_eta.csv"), dtype={"L0": str}))
    Dt.to_csv(os.path.join(a.out, "preview_minus_A1v2.csv"), index=False)
    print(f"# {a.label}: preview minus A1 v2 per density (c_s; errors from c_s_err_scaled and from c_s_err)\n")
    print("| cell | eta | L0 | preview c_s | A1 v2 c_s | difference | relative [%] | error (scaled) | z (scaled) | error | z | note |\n"
          "|---|---|---|---|---|---|---|---|---|---|---|---|")
    for _, r in Dt.iterrows():
        print(f"| {r['cell']} | {float(r['eta']):.6f} | {r['L0']} | {r['c_s_preview']} | {r['c_s_A1v2']} | {r['difference']} | {r['rel_pct']} | "
              f"{r['err_scaled']} | {r['z_scaled']} | {r['err']} | {r['z']} | {r['note']} |")
    print("\nfigure: " + ", ".join(paths))


if __name__ == "__main__":
    main()
