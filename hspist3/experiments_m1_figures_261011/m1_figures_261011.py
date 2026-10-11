#!/usr/bin/env python3
"""##CHRIS 2026-10-11 (plan-author decision 29, task M1; 261012 sec. 4.7.59): the FOUR FIGURES for the meeting of Wednesday 2026-10-14.
It draws from the committed tables of m1_tables_261011.py (this folder's tables/) and from nothing else. The one exception is the
theory lines of Fig. 1, which are computed from the project's EOS module (plot_speed_of_sound_edmd.py) and x_M (tests_20260913.x_of).
The project's plot rules (decision 29):
  - error bars: 95 % (1.96 SE) wherever an SE exists; Fig. 4 has single timings and says so;
  - KR red, SPT and Henderson next to KR, own data blue (C0, as draw_measured_series);
  - phase: none of the four has an eta axis, so the phase is named in each panel, and Fig. 3's solid cells are tinted with the solid
    colour of add_eta_regime_shading;
  - N / radius / L0 / H on every figure;
  - one build generation per panel: Fig. 1 gen2 | gen3 | their difference; Fig. 2 gen2 | gen3; Fig. 3 gen3 only; Fig. 4 gen3 | the
    window's gen3 / gen2 ratio.
  Fig. 1  Test G on KOA. The title states the registered verdict (FAIL at M = 50) and that the replication (Test R) is running.
  Fig. 2  delta_M per engine (Report G2): POST HOC, INFORMATION.
  Fig. 3  The solid: Z per compartment in the cells of P4-P6, the 0-deg and the turned orientation, with the compressed starts' two
          Z groups: PILOT.
  Fig. 4  The KOA rate table (Gold node) and the gen3 / gen2 window speed-up.
PROVENANCE: without --draft the script stops unless tables/ and this script are committed and unmodified, and it prints their commits
on every figure. --draft writes into --out and stamps DRAFT instead. PDF and PNG are written without a creation date, so a re-print
is byte-identical.
usage (from hspist3/): python3 experiments_m1_figures_261011/m1_figures_261011.py [--draft --out DIR]
"""
import argparse, hashlib, math, os, subprocess, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches

HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE); ROOT = os.path.dirname(HS)
sys.path.insert(0, os.path.join(HS, "validation")); sys.path.insert(0, HS)
TAB = os.path.join(HERE, "tables")
BLUE, KR_C, SPT_C, HEN_C, INK, MUTED = "C0", "C3", "C1", "C2", "#222222", "#555555"
SOLID_TINT = "#ded8f5"                                            # add_eta_regime_shading's "solid" colour
Z95 = 1.96
plt.rcParams.update({"font.size": 9, "axes.titlesize": 9.5, "axes.labelsize": 9, "legend.fontsize": 7.5, "xtick.labelsize": 8,
                     "ytick.labelsize": 8, "axes.edgecolor": "#444444", "axes.linewidth": 0.8})


def git(*a):
    return subprocess.run(["git", "-C", ROOT, *a], check=True, capture_output=True, text=True).stdout.strip()


def provenance(draft):
    me = os.path.relpath(os.path.abspath(__file__), ROOT); tab = os.path.relpath(TAB, ROOT)
    if draft: return "DRAFT -- tables and script not yet committed"
    dirty = git("status", "--porcelain", "--", me, tab)
    untracked = [f for f in os.listdir(TAB) if f.endswith(".csv") and not git("ls-files", "--", os.path.join(tab, f))]
    if dirty or untracked: sys.exit(f"STOP: commit the tables and this script first (modified or untracked: {dirty or untracked})")
    return (f"tables: {tab}/ (main {git('log', '-1', '--format=%h', '--', tab)}); script: {me} (main {git('log', '-1', '--format=%h', '--', me)})")


def footer(fig, lines, prov, y=0.005):
    fig.text(0.01, y, "\n".join(lines + [prov]), ha="left", va="bottom", fontsize=7, color=MUTED, wrap=True)


def tag(ax, text, loc="upper left"):
    x, ha = (0.02, "left") if "left" in loc else (0.98, "right")
    ax.text(x, 0.97 if "upper" in loc else 0.03, text, transform=ax.transAxes, ha=ha, va="top" if "upper" in loc else "bottom",
            fontsize=8, fontweight="bold", color=INK, bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#888888", lw=0.6))


def save(fig, out, stem):
    paths = []
    for ext in ("pdf", "png"):
        p = os.path.join(out, f"{stem}.{ext}")
        fig.savefig(p, dpi=200, metadata=({"CreationDate": None} if ext == "pdf" else {"Software": None}))
        paths.append(p)
    plt.close(fig)
    return paths


def mass_axis(ax, masses):
    ax.set_xscale("log"); ax.set_xticks(masses); ax.minorticks_off()
    labs = ax.set_xticklabels([str(m) for m in masses])
    for m, t in zip(masses, labs):                                  # 1500 and 2000 sit close on a log axis
        if m == 1500: t.set_ha("right")
        if m == 2000: t.set_ha("left")
    ax.set_xlim(masses[0] / 1.45, masses[-1] * 1.45); ax.set_xlabel("divider mass M [disk masses]")


def cat_axis(ax, masses):
    ax.set_xticks(range(len(masses))); ax.set_xticklabels([str(m) for m in masses]); ax.set_xlim(-0.6, len(masses) - 0.4)
    ax.set_xlabel("divider mass M [disk masses] (one column per mass)")


# ---------------------------------------------------------------- Fig. 1
def fig1(out, prov):
    import tests_20260913 as T
    import plot_speed_of_sound_edmd as S
    Mn = pd.read_csv(os.path.join(TAB, "fig1_testG_koa_means.csv")); D = pd.read_csv(os.path.join(TAB, "fig1_testG_koa_difference.csv"))
    masses = sorted(Mn.M.unique()); eta = math.pi / 8; L0 = 10.0
    Mg = np.geomspace(masses[0] / 1.4, masses[-1] * 1.4, 300); xg = np.array([T.x_of(m, L0) for m in Mg]); a = np.array([eta])
    cs = {lab: float(S.cs_adiabatic_2d_monatomic(Z(a), dZ(a), a, kbt=1, m=1)[0]) for lab, Z, dZ in
          (("KR", S.Z_kolafa_rottner_2006, S.dZ_kolafa_rottner_2006), ("SPT", S.Z_spt_eos, S.dZ_spt_eos),
           ("Henderson", S.Z_henderson_eos, S.dZ_henderson_eos))}
    fig = plt.figure(figsize=(15, 6.0)); gs = fig.add_gridspec(1, 3, left=0.05, right=0.99, top=0.85, bottom=0.235, wspace=0.22)
    axa = fig.add_subplot(gs[0]); axb = fig.add_subplot(gs[1], sharex=axa, sharey=axa); axc = fig.add_subplot(gs[2])
    for ax, e, lab in ((axa, "gen2", "(a) gen2 (--engine=gen2)"), (axb, "gen3", "(b) gen3 (--engine=gen3)")):
        ax.plot(Mg, cs["KR"] * xg, "-", color=KR_C, lw=2.2, alpha=0.75, zorder=2, label=f"Kolafa–Rottner 2006: c_s x_M (c_s = {cs['KR']:.3f})")
        ax.plot(Mg, cs["SPT"] * xg, "--", color=SPT_C, lw=1.5, alpha=0.75, zorder=3, label=f"SPT: c_s x_M (c_s = {cs['SPT']:.3f})")
        ax.plot(Mg, cs["Henderson"] * xg, "-.", color=HEN_C, lw=1.3, alpha=0.75, zorder=3,
                label=f"Henderson a = 0.125: c_s x_M (c_s = {cs['Henderson']:.3f})")
        for est, mk, fc, f, elab in (("argmax", "o", BLUE, 1 / 1.05, "argmax ν (the registered estimator)"),
                                     ("nu_d", "s", "white", 1.05, "ν_d (fit per trajectory; information)")):
            g = Mn[(Mn.engine == e) & (Mn.estimator == est)].sort_values("M")
            ax.errorbar(g.M * f, g["mean"], yerr=Z95 * g.se, fmt=mk, ms=5.5, mfc=fc, mec=BLUE, color=BLUE, capsize=3, lw=1, zorder=5,
                        label=elab)
        g = Mn[(Mn.engine == e) & (Mn.M == 2000) & (Mn.estimator == "argmax")].iloc[0]
        ax.annotate(f"M = 2000: information\n({int(g.n)} per engine)", (2000, g["mean"]), xytext=(-6, -26), textcoords="offset points",
                    ha="right", fontsize=7, color=MUTED)
        ax.set_xscale("log"); ax.set_yscale("log"); mass_axis(ax, masses)
        ax.set_ylabel("mean ν [reduced units]"); ax.set_title(lab, loc="left")
        ax.text(0.98, 0.97, "η = π/8: dense isotropic fluid", transform=ax.transAxes, ha="right", va="top", fontsize=7.5, color=MUTED)
        ax.grid(True, which="major", color="#e6e6e6", lw=0.6); ax.legend(loc="lower left", frameon=False)
        ax.text(0.98, 0.89, "95 % error bars are smaller than the markers\n(at most 0.27 % of ν)", transform=ax.transAxes, ha="right",
                va="top", fontsize=7, color=MUTED)
    yt = [0.015, 0.02, 0.03, 0.05, 0.07]; axa.set_yticks(yt); axa.set_yticklabels([f"{v:g}" for v in yt]); axa.minorticks_off()
    plt.setp(axb.get_yticklabels(), visible=True)
    ax = axc; ax.axhline(0, color="#888888", lw=0.8, zorder=1); col = {M: i for i, M in enumerate(masses)}
    for est, mk, fc, dx, elab in (("argmax", "o", BLUE, -0.15, "argmax ν (registered at M = 50, 300, 1500)"),
                                  ("nu_d", "s", "white", 0.15, "ν_d (information)")):
        g = D[D.estimator == est].sort_values("M")
        ax.errorbar([col[m] + dx for m in g.M], g.rel_pct, yerr=[g.rel_pct - g.lo95_pct, g.hi95_pct - g.rel_pct], fmt=mk, ms=5.5, mfc=fc,
                    mec=BLUE, color=BLUE, capsize=3, lw=1, zorder=5, label=elab)
    for _, d in D.iterrows():
        reg = isinstance(d.outcome, str) and d.outcome != ""
        if d.estimator == "argmax":
            txt = f"{d.outcome}\nz = {d.z:+.2f}" if reg else f"information (n = {int(d.n_gen3)})\nz = {d.z:+.2f}"
            ax.annotate(txt, (col[d.M] - 0.15, d.hi95_pct), xytext=(0, 5), textcoords="offset points", ha="center", va="bottom",
                        fontsize=7, color=INK, fontweight="bold" if d.outcome == "DIFFERENCE" else "normal")
        else:
            ax.annotate(f"z = {d.z:+.2f}", (col[d.M] + 0.15, d.lo95_pct), xytext=(0, -4), textcoords="offset points", ha="center",
                        va="top", fontsize=6.8, color=MUTED)
    cat_axis(ax, masses); ax.set_ylabel("gen3 − gen2, relative to gen2's mean [%]")
    ax.set_title("(c) gen3 − gen2 per mass, 95 % intervals", loc="left"); ax.grid(True, which="major", color="#e6e6e6", lw=0.6)
    lo, hi = D.lo95_pct.min(), D.hi95_pct.max(); ax.set_ylim(lo - 0.22, hi + 0.3); ax.legend(loc="upper right", frameon=False)
    fig.suptitle("Test G on KOA: registered verdict FAIL at M = 50 (DIFFERENCE); NO DIFFERENCE at 300 and 1500; "
                 "replication (Test R) running", fontsize=11.5, fontweight="bold", x=0.01, ha="left", y=0.975)
    fig.text(0.01, 0.905, "N = 100 disks (50 per compartment), radius 0.5, L0 = 10, H = 10, η = π/8; one KOA binary, engine-gen3 a643bc7 "
             "(target koa), run as gen2 (a) and gen3 (b); the same seeds in both arms, 400 trajectories per engine and mass", fontsize=8.5, color=INK)
    footer(fig, [
        "Cell epi8_H_H10_L10 (the cell of Test T and of methods sec. 16): N = 100, N_s = 50 per compartment, disk radius r = 0.5, L0 = 10, "
        "H = 10, wall thickness t = 0.05, eta = pi/8 = 0.3927. KOA, 13 nodes (Report G-KOA's inventory), build line '00ALLINONE  git "
        "a643bc7  target koa' in every cell.",
        "argmax nu: the registered estimator (red_nu.csv, reduce_B.py). nu_d: the sec. 4.4.13 item 4a fit of each trajectory's ACF, started "
        "at the mass's pooled argmax mean (the verdict's information procedure); 0 fit failures in 2600. M = 2000: 100 per engine, information.",
        "Error bars 95 %: (a, b) mean +- 1.96 SE, SE = SD/sqrt(n); (c) 100 (d -+ 1.96 s) / mean_gen2 with s = sqrt(SE_gen3^2 + SE_gen2^2), "
        "the verdict's interval. Rule per mass: NO DIFFERENCE if |z| < 2, DIFFERENCE if |z| >= 3 (sec. 4.7.29). Values reproduce Report "
        "G-KOA's printed tables (sec. 4.7.54).",
        f"Theory lines: nu = c_s(eta) x_M, c_s adiabatic (kT = m = 1) from each EOS at eta = pi/8; x_M = K(alpha)/(2 pi L_eff), cot K = alpha K, "
        f"alpha = M/(2 N_s), L_eff = L0 - 2r - t/2 = {T.l_eff(L0):.3f} (tests_20260913.x_of).",
    ], prov)
    return save(fig, out, "m1_fig1_testG_koa")


# ---------------------------------------------------------------- Fig. 2
def fig2(out, prov):
    D = pd.read_csv(os.path.join(TAB, "fig2_deltaM_G2.csv")); masses = sorted(D[D.engine.isin(["gen2", "gen3"])].M.unique())
    fig = plt.figure(figsize=(12.5, 6.0)); gs = fig.add_gridspec(1, 2, left=0.065, right=0.99, top=0.82, bottom=0.21, wspace=0.08)
    ax0 = fig.add_subplot(gs[0]); ax1 = fig.add_subplot(gs[1], sharex=ax0, sharey=ax0)
    style = (("nu_arg", "o", BLUE, 1 / 1.08, "argmax ν"), ("nu_d", "s", "white", 1.0, "ν_d"), ("nu_0", "D", "#9ecae1", 1.08, "ν_0"))
    for ax, e, lab in ((ax0, "gen2", "(a) gen2"), (ax1, "gen3", "(b) gen3")):
        ax.axhline(0, color="#888888", lw=0.8, zorder=1)
        ax.axvspan(240, 2600, color="#f0f0f0", zorder=0, lw=0)
        ax.text(math.sqrt(240 * 2600), 0.88, "reference masses (M ≥ 300)", transform=ax.get_xaxis_transform(), ha="center",
                fontsize=7, color=MUTED)
        for est, mk, fc, f, elab in style:
            g = D[(D.engine == e) & (D.estimator == est)].sort_values("M")
            ax.errorbar(g.M * f, g.delta_pct, yerr=Z95 * g.se_pct, fmt=mk, ms=5.5, mfc=fc, mec=BLUE, color=BLUE, capsize=3, lw=1, zorder=5,
                        label=f"{elab} (c_ref = {g.c_ref.iloc[0]:.5f})")
        if e == "gen2":
            s = D[D.engine == "gen2 (methods sec. 16)"].iloc[0]
            ax.errorbar([50 / 1.3], [s.delta_pct], yerr=[Z95 * s.se_pct], fmt="v", ms=6, mfc="white", mec=INK, color=INK, capsize=3, lw=1,
                        zorder=5, label="methods sec. 16: ν_0, Test T legacy (gen2), nine-mass reference")
        mass_axis(ax, masses); ax.set_xlim(50 / 1.7, masses[-1] * 1.45)
        ax.set_title(lab, loc="left"); ax.grid(True, which="major", color="#e6e6e6", lw=0.6)
        ax.legend(loc="lower center", frameon=True, framealpha=0.92, edgecolor="#cccccc")
        tag(ax, "POST HOC, INFORMATION")
        ax.text(0.98, 0.97, "η = π/8: dense isotropic fluid", transform=ax.transAxes, ha="right", va="top", fontsize=7.5, color=MUTED)
    ax0.set_ylabel("δ_M = c_M / c_ref − 1 [%]"); plt.setp(ax1.get_yticklabels(), visible=False)
    ax0.set_ylim(D.lo95_pct.min() - 0.08, D.hi95_pct.max() + 0.16)
    fig.suptitle("δ_M per engine from Report G2 (Test G on KOA): POST HOC, INFORMATION", fontsize=11.5, fontweight="bold", x=0.01, ha="left", y=0.975)
    fig.text(0.01, 0.875, "N = 100 disks (50 per compartment), radius 0.5, L0 = 10, H = 10, η = π/8; one KOA binary, engine-gen3 a643bc7 "
             "(target koa),\nrun as gen2 (a) and gen3 (b); 400 trajectories per engine and mass (M = 2000: 100)", fontsize=8.5, color=INK)
    footer(fig, [
        "c_M = nu_M / x_M with x_M = K(alpha)/(2 pi L_eff) (tests_20260913.x_of, L_eff = 8.975); delta_M = c_M / c_ref - 1, c_ref the "
        "inverse-variance weighted mean of c_M over M >= 300 (300, 1500, 2000) of the same engine and estimator.",
        "Estimators: argmax nu (registered); nu_d (sec. 4.4.13 item 4a fit per trajectory); nu_0 = sqrt(nu_d^2 + (1/(2 pi tau_r))^2) "
        "(methods sec. 16). SE: delete-one-group jackknife over the trajectories (groups of 10; M = 2000: 5), the reference's own error included.",
        "Error bars 95 % (1.96 SE). Values exactly as Report G2 prints them (sec. 4.7.56, main 76f1f115). Methods sec. 16's -0.424 +- 0.078 % "
        "is against a nine-mass reference, not against G2's M >= 300 reference.",
    ], prov)
    return save(fig, out, "m1_fig2_deltaM_G2")


# ---------------------------------------------------------------- Fig. 3
ORDER = (("P6 turned", "turned\nstart\n(P6)"), ("P6 0-deg", "0° start\n(P6,\naddition)"), ("P4 A", "present\nM4 family\n(P4 A)"),
         ("P4a A [0, 5000)", "P4 A,\nfirst quarter\n[0, 5000)"), ("P4 B u 0.01", "compressed\nu 0.01\n(P4 B)"),
         ("P5 B u 0.001", "compressed\nu 0.001\n(P5 B)"), ("P5 B u 0.0001", "compressed\nu 0.0001\n(P5 B)"))


def fig3(out, prov):
    Z = pd.read_csv(os.path.join(TAB, "fig3_solid_Z_per_compartment.csv")); G = pd.read_csv(os.path.join(TAB, "fig3_solid_groups.csv"))
    cells = list(dict.fromkeys(Z.cell))
    fig = plt.figure(figsize=(16, 7.0)); gs = fig.add_gridspec(1, 3, left=0.045, right=0.995, top=0.82, bottom=0.25, wspace=0.17)
    for k, cell in enumerate(cells):
        ax = fig.add_subplot(gs[k]); C = Z[Z.cell == cell]; c0 = C.iloc[0]
        ax.add_patch(matplotlib.patches.Rectangle((0, 0), 1, 1, transform=ax.transAxes, color=SOLID_TINT, alpha=0.32, lw=0, zorder=0))
        starts = [s for s in ORDER if s[0] in set(C.start)]; xs = {s[0]: i for i, s in enumerate(starts)}
        comp = [xs[s] for s in xs if C[C.start == s].kind.iloc[0] == "compressed"]
        cut = G[(G.cell == cell) & G.row.str.startswith("compressed")]
        if len(cut):
            x0, x1 = min(comp) - 0.45, max(comp) + 0.45
            for _, g in cut.iterrows():
                v = C[(C.kind == "compressed") & (C.group == g.row.split()[1])].Z
                ax.fill_between([x0, x1], v.min() - 0.004, v.max() + 0.004, color="#9e9e9e", alpha=0.22, lw=0, zorder=1)
                up = "higher" in g.row
                ax.text(x1 - 0.02, (v.max() + 0.004) if up else (v.min() - 0.004), f"{g.row.split()[1]} Z group: {int(g.n)} compartments, "
                        f"{g['mean']:.4f} ({g.se:.4f})", ha="right", va="bottom" if up else "top", fontsize=6.8, color=INK)
            cv = cut.cut.iloc[0]
            ax.plot([x0, x1], [cv, cv], "--", color="#555555", lw=0.9, zorder=2)
            ax.text(x0 - 0.03, cv, f"P6's cut {cv:.4f}", ha="right", va="center", fontsize=6.5, color=MUTED)
        t = C[C.start == "P6 turned"].Z.mean()
        ax.axhline(t, color=BLUE, lw=0.8, ls=":", alpha=0.8, zorder=2)
        zr = C.Z.max() - C.Z.min(); ax.set_ylim(C.Z.min() - 0.10 * zr, C.Z.max() + 0.16 * zr)
        for s, _ in starts:
            S_ = C[C.start == s].reset_index(drop=True); x = xs[s] + np.linspace(-0.27, 0.27, len(S_)) - 0.06
            for o, mk, fc in (("0 deg", "o", BLUE), ("turned", "^", BLUE), ("not measured", "o", "white")):
                m = (S_.orientation == o).to_numpy()
                if m.any(): ax.plot(x[m], S_.Z[m], mk, ms=5, mfc=fc, mec=BLUE, mew=0.9, ls="none", zorder=4)
            mu, se = S_.Z.mean(), S_.Z.std(ddof=1) / math.sqrt(len(S_))
            ax.errorbar([xs[s] + 0.33], [mu], yerr=[Z95 * se], fmt="_", ms=9, color=INK, capsize=3, lw=1.1, zorder=5)
        ax.set_xticks(range(len(starts))); ax.set_xticklabels([lab for _, lab in starts], fontsize=7.2)
        ax.set_xlim(-0.6, len(starts) - 0.4); ax.grid(True, axis="y", color="#e6e6e6", lw=0.6, zorder=0)
        ax.set_title(f"N = {c0.N} ({c0.N // 2} per compartment), η = {c0.eta:.3f}\nradius 0.5, L0 = {c0.L0:.3f}, H = {c0.H:g}", loc="left")
        ax.text(0.98, 0.985, "solid (η ≥ 0.720)", transform=ax.transAxes, ha="right", va="top", fontsize=7.5, color=MUTED)
        if k == 0: ax.set_ylabel("Z_c per compartment (virial, Q4 window)")
        tag(ax, "PILOT", "upper left")
        if cell == "N400_eta0.740":
            ax.text(0.98, 0.93, "P6 printed no two-group split at N = 400", transform=ax.transAxes, ha="right", va="top", fontsize=6.8,
                    color=MUTED)
    h = [plt.Line2D([], [], marker="o", ls="none", mfc=BLUE, mec=BLUE, label="measured 0°: close-packed lines along x"),
         plt.Line2D([], [], marker="^", ls="none", mfc=BLUE, mec=BLUE, label="measured ±30°: turned (lines along y)"),
         plt.Line2D([], [], marker="o", ls="none", mfc="white", mec=BLUE, label="orientation not measured (P4/P5: seeds 2-5)"),
         plt.Line2D([], [], marker="_", ms=9, ls="none", color=INK, label="mean ± 1.96 SE over the 12 compartments"),
         plt.Line2D([], [], ls=":", color=BLUE, label="the turned start's mean")]
    fig.legend(handles=h, loc="upper right", bbox_to_anchor=(0.995, 0.935), ncol=5, frameon=False, fontsize=7.5)
    fig.suptitle("The solid: Z per compartment, 0° and turned orientation, in the cells of P4–P6: PILOT (no verdict, no paper use)",
                 fontsize=11.5, fontweight="bold", x=0.01, ha="left", y=0.975)
    footer(fig, [
        "Z_c: the dt-weighted mean of each compartment's virial blocks in Q4 = [1.5e4, 2e4) sigma-time of a 2e4 sigma-time hold with the divider "
        "held (M = 300); for 'P4 A, first quarter' its window [0, 5000), before any compartment turned. Two compartments per run, 6 seeds per start.",
        "Marker = the orientation MEASURED in that window (P6: psi6 per compartment, every seed; P4, P5: snapshots of seeds 0 and 1 only), not "
        "the start's: at N = 1600, eta = 0.740 every compartment of the 0° start reads +-30° in Q4 (it turned).",
        "Grey bands (N = 1600): the compressed starts' compartments in the two Z groups, split at the cut P6 printed (the largest gap); "
        "the label gives each group's mean (SE). Box: L0 = N pi / (8 H eta) exactly (--gen3-exact-box), radius 0.5.",
        "Builds, all gen3: P6 engine-gen3 f72d764; P4, P4a and P5 engine-gen3 06e7248. Values as printed in P6's tables "
        "(experiments_gen3_p6_261010/p6_tables_output.txt, engine-gen3 c540b746; sec. 4.7.43).",
    ], prov)
    return save(fig, out, "m1_fig3_solid_Z_P4_P6")


# ---------------------------------------------------------------- Fig. 4
def fig4(out, prov):
    R = pd.read_csv(os.path.join(TAB, "fig4_koa_rates_gold.csv")); S = pd.read_csv(os.path.join(TAB, "fig4_window_speedup.csv"))
    fig = plt.figure(figsize=(13.5, 6.0)); gs = fig.add_gridspec(1, 2, left=0.06, right=0.99, top=0.84, bottom=0.24, wspace=0.2,
                                                                 width_ratios=[1.25, 1])
    ax = fig.add_subplot(gs[0]); r0 = R.iloc[0]
    shades = {"pi/8": ("#6baed6", "η = π/8 (fluid)"), "0.70": ("#2171b5", "η = 0.70 (fluid limit)"), "0.78": ("#08306b", "η = 0.78 (solid)")}
    for lab, (col, name) in shades.items():
        g = R[R.eta_label == lab].sort_values("N")
        ax.plot(g.N, g.events_per_s_engine / 1e5, "-o", color=col, ms=5, lw=1.4, zorder=4)
        ax.plot(g.N, g.events_per_s_run / 1e5, "--o", color=col, ms=5, mfc="white", lw=1.1, zorder=4)
        ax.annotate(name, (g.N.iloc[-1], g.events_per_s_engine.iloc[-1] / 1e5), xytext=(6, {"pi/8": 0, "0.70": -7, "0.78": 7}[lab]),
                    textcoords="offset points", va="center", fontsize=7.5, color=col)
    ax.set_xscale("log"); ax.set_xticks(sorted(R.N.unique())); ax.set_xticklabels([str(n) for n in sorted(R.N.unique())]); ax.minorticks_off()
    ax.set_xlim(80, 3600); ax.set_ylim(0, None); ax.set_xlabel("N (disks, N/2 per compartment)"); ax.set_ylabel("events per second [10^5]")
    ax.grid(True, color="#e6e6e6", lw=0.6)
    ax.legend(handles=[plt.Line2D([], [], color=INK, marker="o", lw=1.4, label="per second of engine time"),
                       plt.Line2D([], [], color=INK, marker="o", mfc="white", ls="--", lw=1.1, label="per second of the whole run")],
              loc="upper right", frameon=False)
    ax.set_title(f"(a) KOA, Gold node {r0.node} ({r0.cpu}): gen3 rate table, one run per cell", loc="left")
    ax.text(0.02, 0.03, "no error bars: one timing per cell", transform=ax.transAxes, fontsize=7, color=MUTED)
    ax = fig.add_subplot(gs[1]); Ns = sorted(S.N.unique()); w = 0.36
    for j, (ts, col) in enumerate(((1, "#9ecae1"), (10000, BLUE))):
        g = S[S.time_scale == ts].set_index("N").loc[Ns]
        xb = np.arange(len(Ns)) + (j - 0.5) * w
        ax.bar(xb, g.speedup_sigma_time, width=w * 0.94, color=col, label=f"time scale {ts}", zorder=3)
        for x_, v, f3, f2 in zip(xb, g.speedup_sigma_time, g.gen3_frames, g.gen2_frames):
            ax.text(x_, v * 1.08, f"{v:.3g}×\n({f3} / {f2} frames)", ha="center", va="bottom", fontsize=7, color=INK)
    ax.set_yscale("log"); ax.set_ylim(0.5, 2000); ax.axhline(1, color="#888888", lw=0.8, zorder=2)
    ax.set_xticks(range(len(Ns))); ax.set_xticklabels([f"N = {n}\nH = {20 * math.sqrt(n / 400):.2f}, L0 = {25.34 * math.sqrt(n / 400):.2f}" for n in Ns])
    ax.set_ylabel("gen3 / gen2: σ-time advanced per wall second"); ax.grid(True, axis="y", color="#e6e6e6", lw=0.6, zorder=0)
    ax.legend(loc="upper left", frameon=False)
    ax.set_title("(b) the window on the Mac: gen3 / gen2 speed-up (one 30-s run each)", loc="left")
    ax.text(0.98, 0.97, "η = 0.31: dense isotropic fluid", transform=ax.transAxes, ha="right", va="top", fontsize=7.5, color=MUTED)
    fig.suptitle("Speed: the KOA rate table (Gold node) and the gen3 / gen2 window speed-up", fontsize=11.5, fontweight="bold", x=0.01,
                 ha="left", y=0.975)
    fig.text(0.01, 0.9, "(a) H = 10√(N/100), L0 = N π / (8 H η) (exact box), radius 0.5, M = 300; (b) radius 0.5, η = 0.31, "
             "H = 20√(N/400), L0 = 25.34√(N/400), wall mass factor 50", fontsize=8.5, color=INK)
    footer(fig, [
        f"(a) K2-2 part (c), job {r0.job}: engine-gen3 a643bc7 (target koa), --engine=gen3; 500 sigma-time held, 100 released; one run per cell, "
        "one at a time. Solid: events per second of the run record's engine_s; dashed: per second of its run_s (the difference is the record's "
        "driver share). Printed with it (information): an E5-2670 v3 node runs at 0.709 of this node's rate (N 400, eta 0.70).",
        "(b) Decision 19 item 2 (sec. 4.7.50): one binary, engine-gen3 7f6adca (target mac-O3-e0pre), --engine=gen2 and --engine=gen3; SDL's dummy "
        "video driver and the software renderer, --auto-release, 30 wall seconds per run, one run each. Bar label: the ratio, and the frames "
        "drawn by gen3 / gen2.",
        "At time scale 1 the fps ratio is 1.18 and 1.29; sec. 4.7.50's INFERENCE: both engines are limited there by drawing in software. "
        "No error bars: single timings.",
    ], prov)
    return save(fig, out, "m1_fig4_rates_speedup")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--draft", action="store_true"); ap.add_argument("--out", default=os.path.join(HERE, "figures"))
    a = ap.parse_args()
    if not a.draft and os.path.abspath(a.out) != os.path.join(HERE, "figures"): sys.exit("STOP: --out only with --draft")
    prov = provenance(a.draft); os.makedirs(a.out, exist_ok=True)
    print(f"# M1 (decision 29): the four meeting figures, printed by experiments_m1_figures_261011/m1_figures_261011.py\n\n{prov}\n")
    print("| figure | file | SHA-256 |\n|---|---|---|")
    for k, f in enumerate((fig1, fig2, fig3, fig4), 1):
        for p in f(a.out, prov):
            print(f"| {k} | {os.path.basename(p)} | {hashlib.sha256(open(p, 'rb').read()).hexdigest()} |")


if __name__ == "__main__":
    main()
