#!/usr/bin/env python3
"""##CHRIS 2026-10-11 (plan-author decision 35, item 4c; 261012 sec. 4.7.68): THE M1 UPDATE -- figures 1 and 2 with Test R, for the
meeting of Wednesday 2026-10-14. The accepted M1 figures (m1_figures_261011.py, c11b5e48) stay as they are; this script writes two new
ones, from the committed tables only (this folder's tables/: the M1 tables and m1_tables_testR_261011.py's), with the helpers of
m1_figures_261011.py imported unchanged.
  m1_fig1_testG_testR_koa   title "Test G: FAIL at M = 50; replication Test R: ENGINES AGREE". (a) gen2, (b) gen3: the mean nu per mass
                            of both tests and both estimators, with KR, SPT and Henderson as nu = c_s x_M; (c) gen3 - gen2 per mass, both
                            tests and both estimators, 95 % intervals, each test's registered outcome; Test R's NO DIFFERENCE band for
                            nu_d (its 95 % interval within +-0.15 %) shaded where Test R ran
  m1_fig2_deltaM_G2_testR   delta_M per engine: the top row Report G2's (POST HOC, INFORMATION; c_ref = the weighted mean over M >= 300),
                            the bottom row Test R's information rows (INFORMATION; against c_300); the same axes in all four panels
The plot rules as M1's (error bars 95 %; KR red, own data blue; the phase named; N, radius, L0, H on the figure; one build generation per
panel). PROVENANCE as M1's: without --draft it stops unless the tables and this script are committed, and prints their commits.
usage (from hspist3/): python3 experiments_m1_figures_261011/m1_figures_v2_261011.py [--draft --out DIR]
"""
import argparse, hashlib, math, os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import m1_figures_261011 as F                                     # the M1 helpers and style, unchanged
plt = F.plt; TAB, ROOT, git = F.TAB, F.ROOT, F.git
BLUE, DARK, INK, MUTED, Z95 = F.BLUE, "#08306b", F.INK, F.MUTED, F.Z95


def provenance(draft):
    me = os.path.relpath(os.path.abspath(__file__), ROOT); tab = os.path.relpath(TAB, ROOT)
    if draft: return "DRAFT -- tables and script not yet committed"
    dirty = git("status", "--porcelain", "--", me, tab)
    untracked = [f for f in os.listdir(TAB) if f.endswith(".csv") and not git("ls-files", "--", os.path.join(tab, f))]
    if dirty or untracked: sys.exit(f"STOP: commit the tables and this script first (modified or untracked: {dirty or untracked})")
    return f"tables: {tab}/ (main {git('log', '-1', '--format=%h', '--', tab)}); script: {me} (main {git('log', '-1', '--format=%h', '--', me)})"


def fig1(out, prov):
    import tests_20260913 as T
    import plot_speed_of_sound_edmd as S
    G = pd.read_csv(os.path.join(TAB, "fig1_testG_koa_means.csv")); GD = pd.read_csv(os.path.join(TAB, "fig1_testG_koa_difference.csv"))
    R = pd.read_csv(os.path.join(TAB, "fig1_testR_koa_means.csv")); RD = pd.read_csv(os.path.join(TAB, "fig1_testR_koa_difference.csv"))
    R["estimator"] = R["estimator"].replace({"nu_arg": "argmax"}); RD["estimator"] = RD["estimator"].replace({"nu_arg": "argmax"})
    masses = sorted(set(G.M) | set(R.M)); eta = math.pi / 8; L0 = 10.0
    Mg = np.geomspace(masses[0] / 1.4, masses[-1] * 1.4, 300); xg = np.array([T.x_of(m, L0) for m in Mg]); a = np.array([eta])
    cs = {lab: float(S.cs_adiabatic_2d_monatomic(Z(a), dZ(a), a, kbt=1, m=1)[0]) for lab, Z, dZ in
          (("KR", S.Z_kolafa_rottner_2006, S.dZ_kolafa_rottner_2006), ("SPT", S.Z_spt_eos, S.dZ_spt_eos),
           ("Henderson", S.Z_henderson_eos, S.dZ_henderson_eos))}
    fig = plt.figure(figsize=(16, 6.4)); gs = fig.add_gridspec(1, 3, left=0.045, right=0.99, top=0.85, bottom=0.235, wspace=0.2,
                                                               width_ratios=[1, 1, 1.35])
    axa = fig.add_subplot(gs[0]); axb = fig.add_subplot(gs[1], sharex=axa, sharey=axa); axc = fig.add_subplot(gs[2])
    sty = (("Test G", G, "argmax", "o", BLUE, BLUE, 1 / 1.07), ("Test G", G, "nu_d", "s", "white", BLUE, 1 / 1.02),
           ("Test R", R, "argmax", "D", DARK, DARK, 1.02), ("Test R", R, "nu_d", "D", "white", DARK, 1.07))
    for ax, e, lab in ((axa, "gen2", "(a) gen2 (--engine=gen2)"), (axb, "gen3", "(b) gen3 (--engine=gen3)")):
        ax.plot(Mg, cs["KR"] * xg, "-", color=F.KR_C, lw=2.2, alpha=0.75, zorder=2, label=f"Kolafa–Rottner 2006: c_s x_M ({cs['KR']:.3f})")
        ax.plot(Mg, cs["SPT"] * xg, "--", color=F.SPT_C, lw=1.5, alpha=0.75, zorder=3, label=f"SPT: c_s x_M ({cs['SPT']:.3f})")
        ax.plot(Mg, cs["Henderson"] * xg, "-.", color=F.HEN_C, lw=1.3, alpha=0.75, zorder=3, label=f"Henderson a = 0.125 ({cs['Henderson']:.3f})")
        for test, Dm, est, mk, fc, ec, f in sty:
            g = Dm[(Dm.engine == e) & (Dm.estimator == est)].sort_values("M")
            ax.errorbar(g.M * f, g["mean"], yerr=Z95 * g.se, fmt=mk, ms=5, mfc=fc, mec=ec, color=ec, capsize=2.5, lw=1, zorder=5,
                        label=f"{test}: {'argmax ν' if est == 'argmax' else 'ν_d'}")
        ax.set_xscale("log"); ax.set_yscale("log"); F.mass_axis(ax, masses); ax.set_ylabel("mean ν [reduced units]")
        ax.set_title(lab, loc="left"); ax.grid(True, which="major", color="#e6e6e6", lw=0.6); ax.legend(loc="lower left", frameon=False, fontsize=7)
        ax.text(0.98, 0.97, "η = π/8: dense isotropic fluid", transform=ax.transAxes, ha="right", va="top", fontsize=7.5, color=MUTED)
        ax.text(0.98, 0.9, "95 % error bars are smaller\nthan the markers", transform=ax.transAxes, ha="right", va="top", fontsize=7, color=MUTED)
    yt = [0.015, 0.02, 0.03, 0.05, 0.07]; axa.set_yticks(yt); axa.set_yticklabels([f"{v:g}" for v in yt]); axa.minorticks_off()
    ax = axc; ax.axhline(0, color="#888888", lw=0.8, zorder=1); col = {M: i for i, M in enumerate(masses)}
    for M in sorted(set(R.M)):
        ax.fill_between([col[M] - 0.45, col[M] + 0.45], -0.15, 0.15, color="#dddddd", alpha=0.6, lw=0, zorder=0)
    rows = (("Test G", GD.rename(columns={"rel_pct": "d"}), "argmax", "o", BLUE, BLUE, -0.27), ("Test G", GD.rename(columns={"rel_pct": "d"}), "nu_d", "s", "white", BLUE, -0.09),
            ("Test R", RD.rename(columns={"d_pct": "d"}), "argmax", "D", DARK, DARK, 0.09), ("Test R", RD.rename(columns={"d_pct": "d"}), "nu_d", "D", "white", DARK, 0.27))
    for test, Dd, est, mk, fc, ec, dx in rows:
        g = Dd[Dd.estimator == est].sort_values("M")
        ax.errorbar([col[m] + dx for m in g.M], g.d, yerr=[g.d - g.lo95_pct, g.hi95_pct - g.d], fmt=mk, ms=5, mfc=fc, mec=ec, color=ec,
                    capsize=2.5, lw=1, zorder=5, label=f"{test}: {'argmax ν' if est == 'argmax' else 'ν_d'}")
    for _, d in GD[GD.estimator == "argmax"].iterrows():
        o = d.outcome if isinstance(d.outcome, str) and d.outcome else "information"
        hi = GD[GD.M == d.M].hi95_pct.max()
        ax.annotate(f"G: {o}\nz = {d.z:+.2f}", (col[d.M] - 0.18, hi), xytext=(0, 4), textcoords="offset points", ha="center", va="bottom",
                    fontsize=6.6, color=INK, fontweight="bold" if o == "DIFFERENCE" else "normal")
    for M in sorted(set(RD.M)):
        g = RD[RD.M == M]; lo = g.lo95_pct.min()
        ax.annotate(f"R: {g.outcome.iloc[0]}\nz {g.z.iloc[0]:+.2f} / {g.z.iloc[1]:+.2f}", (col[M] + 0.18, lo), xytext=(0, -4),
                    textcoords="offset points", ha="center", va="top", fontsize=6.6, color=DARK)
    F.cat_axis(ax, masses); ax.set_ylabel("gen3 − gen2, relative [%]")
    ax.set_title("(c) gen3 − gen2 per mass, 95 % intervals", loc="left")
    ax.grid(True, axis="y", color="#e6e6e6", lw=0.6); ax.set_ylim(min(GD.lo95_pct.min(), RD.lo95_pct.min()) - 0.3, max(GD.hi95_pct.max(), RD.hi95_pct.max()) + 0.3)
    h, l = ax.get_legend_handles_labels()
    h.append(F.matplotlib.patches.Patch(color="#dddddd", alpha=0.6)); l.append("Test R's NO DIFFERENCE band for ν_d (±0.15 %)")
    ax.legend(h, l, loc="upper right", frameon=False, fontsize=7, ncol=2)
    fig.suptitle("Test G: FAIL at M = 50; replication Test R: ENGINES AGREE", fontsize=12.5, fontweight="bold", x=0.01, ha="left", y=0.975)
    fig.text(0.01, 0.905, "N = 100 disks (50 per compartment), radius 0.5, L0 = 10, H = 10, η = π/8; one KOA binary, engine-gen3 a643bc7 (target koa), "
             "run as gen2 (a) and gen3 (b); the same seeds in both arms. Test G: 400 trajectories per engine and mass (M = 2000: 100); "
             "Test R: 3600", fontsize=8.3, color=INK)
    F.footer(fig, [
        "Cell epi8_H_H10_L10: N = 100, N_s = 50, disk radius 0.5, L0 = 10, H = 10, wall thickness 0.05, eta = pi/8. Test G's rule (sec. 4.7.29): "
        "NO DIFFERENCE if |z| < 2, DIFFERENCE if |z| >= 3, argmax nu, at M = 50, 300, 1500 (M = 2000 information). Test R's rule (decision 28): "
        "per mass DIFFERENCE CONFIRMED if |z| >= 3.5 by either estimator, NO DIFFERENCE if both |z| < 3.5 and nu_d's 95 % interval lies within +-0.15 %.",
        "Error bars 95 %: (a, b) mean +- 1.96 SE; (c) Test G: 100 (d -+ 1.96 s) / mean_gen2 (its verdict); Test R: d -+ 1.96 SE(d), delta method "
        "(its verdict). nu_d: the sec. 4.4.13 item 4a fit per trajectory (Test G: information; Test R: registered). Values reproduce both "
        "verdicts' printed tables (secs. 4.7.54, 4.7.66). Test G and Test R pooled (Report R1): M = 50, nu_arg +0.107 +- 0.058 % (z 1.86), "
        "nu_d -0.006 +- 0.035 %.",
        "Theory lines: nu = c_s(eta) x_M, c_s adiabatic (kT = m = 1) at eta = pi/8; x_M = K(alpha)/(2 pi L_eff), L_eff = 8.975 (tests_20260913.x_of).",
    ], prov)
    return F.save(fig, out, "m1_fig1_testG_testR_koa")


def fig2(out, prov):
    G2 = pd.read_csv(os.path.join(TAB, "fig2_deltaM_G2.csv")); TR = pd.read_csv(os.path.join(TAB, "fig2_deltaM_testR.csv"))
    masses = sorted(set(G2[G2.engine.isin(["gen2", "gen3"])].M.unique()) | set(TR.M.unique()))
    fig = plt.figure(figsize=(12.5, 9.4)); gs = fig.add_gridspec(2, 2, left=0.07, right=0.99, top=0.88, bottom=0.15, wspace=0.08, hspace=0.3)
    axes = [[fig.add_subplot(gs[i, j]) for j in range(2)] for i in range(2)]
    for i in range(2):
        for j in range(2):
            if (i, j) != (0, 0): axes[i][j].sharex(axes[0][0]); axes[i][j].sharey(axes[0][0])
    style = (("nu_arg", "o", BLUE, 1 / 1.08, "argmax ν"), ("nu_d", "s", "white", 1.0, "ν_d"), ("nu_0", "D", "#9ecae1", 1.08, "ν_0"))
    for i, (D, rowlab, tagtxt, ref) in enumerate(((G2, "Report G2 (Test G on KOA)", "POST HOC, INFORMATION", "c_ref: weighted mean over M ≥ 300"),
                                                  (TR, "Test R (its information rows)", "INFORMATION", "against c_300"))):
        for j, e in enumerate(("gen2", "gen3")):
            ax = axes[i][j]; ax.axhline(0, color="#888888", lw=0.8, zorder=1)
            for est, mk, fc, f, elab in style:
                g = D[(D.engine == e) & (D.estimator == est)].sort_values("M")
                if g.empty: continue
                ax.errorbar(g.M * f, g.delta_pct, yerr=Z95 * g.se_pct, fmt=mk, ms=5.5, mfc=fc, mec=BLUE, color=BLUE, capsize=3, lw=1, zorder=5,
                            label=elab + (f" (c_ref = {g.c_ref.iloc[0]:.5f})" if "c_ref" in g and pd.notna(g.c_ref.iloc[0]) else ""))
            if i == 0 and e == "gen2":
                s = D[D.engine == "gen2 (methods sec. 16)"].iloc[0]
                ax.errorbar([50 / 1.3], [s.delta_pct], yerr=[Z95 * s.se_pct], fmt="v", ms=6, mfc="white", mec=INK, color=INK, capsize=3, lw=1,
                            zorder=5, label="methods sec. 16: ν_0 (nine-mass reference)")
            F.mass_axis(ax, masses); ax.set_xlim(50 / 1.7, masses[-1] * 1.45)
            ax.set_title(f"({'abcd'[2 * i + j]}) {rowlab}: {e} ({ref})", loc="left", fontsize=9.5)
            ax.grid(True, which="major", color="#e6e6e6", lw=0.6); ax.legend(loc="lower right", frameon=True, framealpha=0.92, edgecolor="#cccccc", fontsize=7)
            F.tag(ax, tagtxt)
            if j == 0: ax.set_ylabel("δ_M = c_M / c_ref − 1 [%]")
            else: plt.setp(ax.get_yticklabels(), visible=False)
    lo = min(G2.lo95_pct.min(), TR.lo95_pct.min()); hi = max(G2.hi95_pct.max(), TR.hi95_pct.max()); axes[0][0].set_ylim(lo - 0.1, hi + 0.15)
    fig.suptitle("δ_M per engine: Report G2 (POST HOC, INFORMATION) and Test R (INFORMATION)", fontsize=12, fontweight="bold", x=0.01, ha="left", y=0.98)
    fig.text(0.01, 0.915, "N = 100 disks (50 per compartment), radius 0.5, L0 = 10, H = 10, η = π/8 (dense isotropic fluid); one KOA binary, "
             "engine-gen3 a643bc7 (target koa),\nrun as gen2 and gen3; Test G 400 trajectories per engine and mass (M = 2000: 100), Test R 3600",
             fontsize=8.3, color=INK)
    F.footer(fig, [
        "c_M = nu_M / x_M, x_M = K(alpha)/(2 pi L_eff); delta_M = c_M / c_ref - 1. The references differ by row: Report G2's c_ref is the "
        "inverse-variance weighted mean of c_M over M >= 300 (300, 1500, 2000) of the same engine and estimator; Test R's is c_300, its one mass >= 300.",
        "SE: delete-one-group jackknife (groups of 10), the reference's own error included; error bars 95 % (1.96 SE). Values as printed in "
        "Report G2 (main 76f1f115) and in Test R's verdict output (main 7ed5e102, information rows).",
        "Test R at M = 50, nu_d: gen2 -0.389 +- 0.029 %, gen3 -0.428 +- 0.029 %: gen3 - gen2 = -0.039 +- 0.041 % (Report R1).",
    ], prov)
    return F.save(fig, out, "m1_fig2_deltaM_G2_testR")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--draft", action="store_true"); ap.add_argument("--out", default=os.path.join(HERE, "figures"))
    a = ap.parse_args()
    if not a.draft and os.path.abspath(a.out) != os.path.join(HERE, "figures"): sys.exit("STOP: --out only with --draft")
    prov = provenance(a.draft); os.makedirs(a.out, exist_ok=True)
    print(f"# M1 update (decision 35, item 4c): figures 1 and 2 with Test R, printed by experiments_m1_figures_261011/m1_figures_v2_261011.py\n\n{prov}\n")
    print("| figure | file | SHA-256 |\n|---|---|---|")
    for k, f in ((1, fig1), (2, fig2)):
        for p in f(a.out, prov):
            print(f"| {k} | {os.path.basename(p)} | {hashlib.sha256(open(p, 'rb').read()).hexdigest()} |")


if __name__ == "__main__":
    main()
