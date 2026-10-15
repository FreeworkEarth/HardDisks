#!/usr/bin/env python3
"""##CHRIS 2026-10-15 (plan-author decision 54 item 4; registration 261012 sec. 4.8.11, with the list of the plan author's answer in sec.
4.7.111): THE REGISTERED FIGURES of stage 1 Part B, by script, from the fetched folder only (read only). The numbers are the full report's:
s1_verdict_261012.analyse(<root>) over all six groups (not --interim, its jackknife), the registered functions unchanged; nothing here
fits or estimates anything of its own. Per figure a data file (CSV: the numbers drawn) and PNG + PDF, named by date and content.
  fig1  c_s(eta) of E1 for the four sizes with KR (red, drawn to eta 0.69), SPT, Henderson, Roman 2002 (black) and the shaded regions;
        below it Delta = c_s / c_s^KR - 1 [%] with the extrapolated a(eta) of F2 and F1 (sec. 4.8.11 item 1)
  fig2  Delta against u = N^(-1/2) per density with the F2 and F1 lines and a at u = 0 (sec. 4.8.11 item 2; the fits' own variable u)
  fig3  the ideal-gas zoom, eta <= 0.15, blocks B and D, with sqrt(2) (sec. 4.8.11 item 3)
  fig4  Z(eta) for the four sizes with KR, and Delta^Z = Z / Z_KR - 1 [%] with a_Z of F2 (sec. 4.8.11 item 4)
  (sec. 4.8.11 item 5, the identity, waits for Part A: not drawn)
  fig5  delta_0.5(eta, N) [%] (sec. 4.8.11 item 6)
  fig6  Gamma = 2 / tau_r and Q_k per (cell, mass), against eta, one line per alpha, one column per N (sec. 4.8.11 item 6)
  fig7  psi6(eta) per N: the cell's mean +- SE (log scale); beside it eta 0.60-0.69 on a linear scale with every trajectory's
        value as a dot (sec. 4.7.111)
  fig8  Verdict S: the four z per block S cell, with the verdict's 3.5 (sec. 4.7.111)
The conventions of sec. 4.8.11: one build generation (gen3 at f72d764); error bars everywhere; KR red, Roman black, own data blue (the
four sizes in four blues and four marker shapes); SPT and Henderson next to KR; the phase regions shaded (plot_speed_of_sound_edmd's
add_eta_regime_shading); N, the disk radius, L0 and H on the figure. Open symbols: the cells with a stationarity flag (sec. 4.8.7.6),
and in fig6 also the (cell, mass) folders the group rule leaves out (sec. 4.8.5).
usage (from hspist3/): python3 cluster/gen3_stage1_261012/s1_figures_261015.py --data <the fetched folder> --out <a new folder> [--date 261015]
"""
import argparse, csv, math, os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE); sys.path.insert(0, HS)
import s1_verdict_261012 as V                                    # analyse(): the full report's numbers
import plot_speed_of_sound_edmd as SOS                           # the equations of state, Roman 2002, the regions
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker

NS = (100, 400, 900, 1600)
COL = {100: "#6baed6", 400: "#2b7bba", 900: "#0b4f8f", 1600: "#061f4a"}
MK = {100: "o", 400: "s", 900: "^", 1600: "D"}
OFF = {100: -0.0045, 400: -0.0015, 900: 0.0015, 1600: 0.0045}   # eta offsets in figs 5-7 only, for legibility
KRC, ROMC, SPTC, HENC = "#d62728", "black", "#e69f00", "#009e73"
SCRIPT = "cluster/gen3_stage1_261012/s1_figures_261015.py"
GEOM = ("N = 100, 400, 900, 1600 disks, N/2 per side; disk radius 0.5; divider thickness 0.05; H = 10 sqrt(N/100) = 10, 20, 30, 40; "
        "L0 = N pi / (8 H eta) per side (exact box); gen3 engine, one build (f72d764); 25 seeds per (cell, mass), 10 in block S")
nan = math.nan


def cid(blk, N, e, st="default"): return V.GEN.cell_id(blk, N, e, st)
def cs_eos(Zf, dZf, e): return SOS.cs_adiabatic_2d_monatomic(Zf(e), dZf(e), e, kbt=1.0, m=1.0)


def shade(ax, x_max=0.72, x_min=0.0, labels=True):
    """plot_speed_of_sound_edmd's regions; labels=False removes the region names it writes (the shading and the boundaries stay)"""
    n = len(ax.texts); SOS.add_eta_regime_shading(ax, x_max=x_max, x_min=x_min)
    if not labels:
        for t in ax.texts[n:]: t.remove()


def symlog_y(ax, linthresh, linscale, ticks):
    ax.set_yscale("symlog", linthresh=linthresh, linscale=linscale); ax.set_yticks(ticks)
    ax.set_yticklabels([f"{t:g}" for t in ticks]); ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())


def stat_flag(c):
    t, z = c["hold_t"][0], c["rec_z"][0]
    return bool((math.isfinite(t) and abs(t) >= V.T_STAT) or (math.isfinite(z) and abs(z) >= V.T_STAT))


def write_csv(path, rows):
    keys = list(rows[0].keys()) if rows else ["empty"]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)


def footer(fig, data, root, y=0.004):
    fig.text(0.5, y, f"data: {data}\nscript: {SCRIPT} (s1_verdict_261012.analyse on {root}, all six groups, the full report's numbers)\n{GEOM}",
             ha="center", va="bottom", fontsize=7.2, color="0.25")


def save(fig, out, name):
    for ext in ("png", "pdf"): fig.savefig(os.path.join(out, f"{name}.{ext}"), dpi=170 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig); return [f"{name}.png", f"{name}.pdf"]


def pt(ax, x, y, s, N, flagged, ms=6, label=None, zorder=5):
    ax.errorbar([x], [y], yerr=[s], fmt=MK[N], color=COL[N], mfc="white" if flagged else COL[N], mec=COL[N], ms=ms, capsize=2, lw=1,
                zorder=zorder, label=label)


def size_legend(ax, extra=(), loc="upper left", **kw):
    h = [plt.Line2D([], [], marker=MK[N], color=COL[N], ls="-", lw=1, ms=6, label=f"N = {N}") for N in NS]
    h += [plt.Line2D([], [], marker="o", color="0.3", mfc="white", ls="none", ms=6, label="open: stationarity flag (sec. 4.8.7.6)")]
    ax.legend(handles=h + list(extra), loc=loc, fontsize=8, framealpha=0.95, **kw)


# ---------------------------------------------------------------------------------------------------------------- fig1 and fig2
def fig_cs(res, out, root, tag):
    rows, fits = [], []
    for N in NS:
        for e in V.GEN.ETA_B:
            c = res["cell"].get(cid("B", N, e))
            if not c: continue
            E = c["E1"]; kr = V.kr_cs(e)
            rows.append(dict(cell=c["cid"], N=N, eta=e, set="Phi" if e in V.PHI else "approach", E1_cs=E["cs"] if E else nan,
                             E1_sigma=E["sig"] if E else nan, E1_masses=E["n"] if E else 0, E1_chi2=E["chi2"] if E else nan, E1_dof=E["dof"] if E else 0,
                             KR=kr, Delta_pct=100 * (E["cs"] / kr - 1) if E else nan, Delta_sigma_pct=100 * E["sig"] / kr if E else nan,
                             lowQ_flag=int(E["lowQ"]) if E else 0, stationarity_flag=int(stat_flag(c)), hold_t=c["hold_t"][0], record_z=c["rec_z"][0]))
    for e in V.GEN.ETA_B:
        f = res["fitsC"][e]; K = res["K"]
        r = dict(eta=e, set="Phi" if e in V.PHI else "approach", sizes=f["n"])
        for lab in ("F2", "F1", "F1x"):
            F = f.get(lab); nb = 3 if lab == "F2" else 2
            for i in range(nb):
                r[f"{lab}_b{i}_pct"] = 100 * F["beta"][i] if F else nan; r[f"{lab}_se{i}_pct"] = 100 * F["se"][i] if F else nan
            r[f"{lab}_chi2"] = F["chi2"] if F else nan; r[f"{lab}_dof"] = F["dof"] if F else 0; r[f"{lab}_p"] = F["p"] if F else nan
        r["K_z"] = K["z"].get(e, nan) if K else nan
        fits.append(r)
    n1 = f"{tag}_s1_fig1_cs_vs_eta_E1_N100-1600"; n2 = f"{tag}_s1_fig2_Delta_vs_Nm12_per_eta_F1_F2"
    write_csv(os.path.join(out, n1 + ".csv"), rows); write_csv(os.path.join(out, n2 + "_fits.csv"), fits)
    df = pd.DataFrame(rows); ff = pd.DataFrame(fits)
    files = [n1 + ".csv", n2 + "_fits.csv"]
    # fig1
    fig, (ax, ar) = plt.subplots(2, 1, figsize=(13, 11), sharex=True, gridspec_kw=dict(height_ratios=[2.1, 1.5], hspace=0.05))
    shade(ax); shade(ar, labels=False)
    ef = np.linspace(0.002, 0.70, 600); ek = np.linspace(0.002, 0.69, 600)
    ax.plot(ef, cs_eos(SOS.Z_spt_eos, SOS.dZ_spt_eos, ef), "--", color=SPTC, lw=1.5, label="SPT")
    ax.plot(ef, cs_eos(SOS.Z_henderson_eos, SOS.dZ_henderson_eos, ef), "-.", color=HENC, lw=1.5, label="Henderson (a = 0.125)")
    ax.plot(ek, cs_eos(SOS.Z_kolafa_rottner_2006, SOS.dZ_kolafa_rottner_2006, ek), "-", color=KRC, lw=2.4, alpha=0.9,
            label="KR (Kolafa-Rottner 2006, the rho_max 0.90 fit; drawn to eta 0.69)")
    re_, rc, rs = SOS.roman_2002_reference_points(radius=0.5, height=10.0)
    ax.errorbar(re_, rc, yerr=rs, fmt="s", color=ROMC, ms=5, capsize=2, zorder=6, label="Roman et al. 2002 (N = 100)")
    for N in NS:
        d = df[df.N == N].sort_values("eta"); ok = d.E1_cs.notna()
        ax.plot(d.eta[ok], d.E1_cs[ok], "-", color=COL[N], lw=1, zorder=4)
        for r in d[ok].itertuples(): pt(ax, r.eta, r.E1_cs, r.E1_sigma, N, r.stationarity_flag)
    ax.set_ylabel("speed of sound  c_s  [sqrt(kT/m)], E1 +- sigma", fontsize=11); ax.grid(alpha=0.25)
    ax.set_ylim(0, max(df.E1_cs.max() * 1.08, 16)); ax.set_xlim(0, 0.72)
    hh, ll = ax.get_legend_handles_labels()
    size_legend(ax, extra=hh, loc="upper left", bbox_to_anchor=(0.035, 0.88))
    K = res["K"]
    ax.set_title("Stage 1, Part B: speed of sound against packing fraction, E1 (eight masses alpha >= 1), four box sizes\n"
                 + (f"Verdict K on Phi (eta <= 0.60): {'KR-CONSISTENT' if K['consistent'] else 'KR-INCONSISTENT'} (X^2 = {K['X2']:.2f}, {K['dof']} dof, "
                    f"p = {K['p']:.4f}, max |z| = {K['maxz']:.2f}); 0.63-0.69 approach cells, no verdict" if K else "Verdict K: not printed"), fontsize=11)
    # the residual, on a symmetric-log axis: linear within +-2.5 %, logarithmic beyond (the approach cells reach -21 % and +42 %)
    ar.axhline(0, color=KRC, lw=2, alpha=0.9)
    for N in NS:
        d = df[(df.N == N) & df.E1_cs.notna()].sort_values("eta")
        ar.plot(d.eta, d.Delta_pct, "-", color=COL[N], lw=0.8, zorder=4)
        for r in d.itertuples(): pt(ar, r.eta, r.Delta_pct, r.Delta_sigma_pct, N, r.stationarity_flag)
    fp = ff[ff.F2_b0_pct.notna()]
    ar.errorbar(fp.eta - 0.002, fp.F2_b0_pct, yerr=fp.F2_se0_pct, fmt="*", color="k", ms=10, capsize=3, zorder=7,
                label="a(eta), F2 (N -> infinity) +- sigma_a")
    f1 = ff[ff.F1_b0_pct.notna()]
    ar.errorbar(f1.eta + 0.002, f1.F1_b0_pct, yerr=f1.F1_se0_pct, fmt="*", color="0.5", mfc="white", ms=9, capsize=2, zorder=6,
                label="a(eta), F1 (reported)")
    symlog_y(ar, 2.5, 1.8, [-20, -10, -5, -2, -1, 0, 1, 2, 5, 10, 20, 40]); ar.set_ylim(-30, 55)
    ar.axhspan(-2.5, 2.5, color="0.5", alpha=0.06, lw=0, zorder=0)
    ar.text(0.006, 0.5, "linear within +-2.5 %, logarithmic beyond", transform=ar.transAxes, rotation=90, va="center", fontsize=7, color="0.35")
    ar.set_ylabel("Delta = c_s / c_s^KR - 1  [%]", fontsize=11); ar.set_xlabel("packing fraction  eta", fontsize=11)
    ar.grid(alpha=0.25); ar.legend(loc="lower left", fontsize=8, framealpha=0.95, bbox_to_anchor=(0.03, 0.0))
    footer(fig, f"{n1}.csv (cells), {n2}_fits.csv (a of F2 and F1)", root, y=-0.035)
    files += save(fig, out, n1)
    # fig2
    fig, axs = plt.subplots(4, 4, figsize=(17, 16), sharex=True)
    uu = np.linspace(0, 0.105, 200)
    for k, e in enumerate(V.GEN.ETA_B):
        ax = axs.flat[k]; f = res["fitsC"][e]
        ax.axhline(0, color=KRC, lw=1.4)
        for N in NS:
            r = df[(df.N == N) & (df.eta == e)]
            if len(r) and r.E1_cs.notna().iloc[0]:
                r = r.iloc[0]; pt(ax, N ** -0.5, r.Delta_pct, r.Delta_sigma_pct, N, r.stationarity_flag)
        txt = []
        if "F2" in f:
            b = f["F2"]["beta"]; ax.plot(uu, 100 * (b[0] + b[1] * uu + b[2] * uu ** 2), "-", color="k", lw=1.2)
            ax.errorbar([0.0], [100 * b[0]], yerr=[100 * f["F2"]["se"][0]], fmt="*", color="k", ms=10, capsize=3, zorder=7)
            txt.append(f"F2: a = {100 * b[0]:+.3f} +- {100 * f['F2']['se'][0]:.3f} %, z = {b[0] / f['F2']['se'][0]:+.2f}, p = {f['F2']['p']:.3f}")
        if "F1" in f:
            b = f["F1"]["beta"]; ax.plot(uu, 100 * (b[0] + b[1] * uu), "--", color="0.45", lw=1.1)
            ax.errorbar([0.003], [100 * b[0]], yerr=[100 * f["F1"]["se"][0]], fmt="*", color="0.45", mfc="white", ms=9, capsize=2, zorder=6)
            txt.append(f"F1: a = {100 * b[0]:+.3f} +- {100 * f['F1']['se'][0]:.3f} %, p = {f['F1']['p']:.3f}")
        ax.set_title(f"eta = {e:.4g}" + ("  (Phi)" if e in V.PHI else "  (approach cell: no verdict)"), fontsize=10, pad=26)
        ax.text(0.5, 1.012, "\n".join(txt) or "no fit", transform=ax.transAxes, fontsize=7.2, va="bottom", ha="center")
        ax.grid(alpha=0.25); ax.set_xlim(-0.004, 0.106)
        for N in NS: ax.axvline(N ** -0.5, color="0.85", lw=0.6, zorder=0)
        if k >= 12: ax.set_xlabel("u = N^(-1/2)   (0.1, 0.05, 0.033, 0.025 = N 100, 400, 900, 1600)", fontsize=8)
        if k % 4 == 0: ax.set_ylabel("Delta [%]", fontsize=9)
    hl = [plt.Line2D([], [], marker=MK[N], color=COL[N], ls="none", ms=6, label=f"N = {N}") for N in NS]
    hl += [plt.Line2D([], [], color="k", lw=1.2, marker="*", ms=9, label="F2 = a + b u + d u^2 (primary), a +- sigma_a at u = 0"),
           plt.Line2D([], [], color="0.45", lw=1.1, ls="--", marker="*", mfc="white", ms=9, label="F1 = a + b u (reported)"),
           plt.Line2D([], [], color=KRC, lw=1.4, label="KR (Delta = 0)")]
    fig.legend(handles=hl, loc="upper center", ncol=7, fontsize=9, bbox_to_anchor=(0.5, 0.975))
    fig.suptitle("Delta = c_s(N) / c_s^KR - 1 (E1) against u = N^(-1/2) per density, with the finite-size fits of sec. 4.8.7.1\n"
                 + (f"Verdict K on Phi: {'KR-CONSISTENT' if K['consistent'] else 'KR-INCONSISTENT'}, X^2 = {K['X2']:.2f} ({K['dof']} dof), p = {K['p']:.4f}, "
                    f"max |z| = {K['maxz']:.2f}" if K else ""), fontsize=12, y=1.0)
    fig.tight_layout(rect=(0, 0.035, 1, 0.962))
    footer(fig, f"{n1}.csv (the points), {n2}_fits.csv (the fits)", root, y=-0.005)
    files += save(fig, out, n2)
    return files, df


# ---------------------------------------------------------------------------------------------------------------- fig3, the ideal-gas zoom
def fig_zoom(res, out, root, tag):
    rows = []
    for N in NS:
        for blk, etas in (("D", V.GEN.ETA_D), ("B", [e for e in V.GEN.ETA_B if e <= 0.15 + 1e-12])):
            for e in etas:
                c = res["cell"].get(cid(blk, N, e))
                if not c or not c["E1"]: continue
                E = c["E1"]; kr = V.kr_cs(e)
                rows.append(dict(cell=c["cid"], block=blk, N=N, eta=e, E1_cs=E["cs"], E1_sigma=E["sig"], E1_masses=E["n"], KR=kr,
                                 Delta_pct=100 * (E["cs"] / kr - 1), Delta_sigma_pct=100 * E["sig"] / kr,
                                 vs_sqrt2_pct=100 * (E["cs"] / math.sqrt(2) - 1), stationarity_flag=int(stat_flag(c))))
    n = f"{tag}_s1_fig3_cs_ideal_gas_zoom_eta_le_0.15_B_D"; write_csv(os.path.join(out, n + ".csv"), rows); df = pd.DataFrame(rows)
    fig, (ax, ar) = plt.subplots(2, 1, figsize=(11, 9.5), sharex=True, gridspec_kw=dict(height_ratios=[2, 1.3], hspace=0.05))
    shade(ax, x_max=0.16); shade(ar, x_max=0.16, labels=False)
    ef = np.linspace(0.001, 0.16, 400)
    ax.plot(ef, cs_eos(SOS.Z_spt_eos, SOS.dZ_spt_eos, ef), "--", color=SPTC, lw=1.4, label="SPT")
    ax.plot(ef, cs_eos(SOS.Z_henderson_eos, SOS.dZ_henderson_eos, ef), "-.", color=HENC, lw=1.4, label="Henderson (a = 0.125)")
    ax.plot(ef, cs_eos(SOS.Z_kolafa_rottner_2006, SOS.dZ_kolafa_rottner_2006, ef), "-", color=KRC, lw=2.2, alpha=0.9, label="KR")
    ax.axhline(math.sqrt(2), color="0.3", ls=":", lw=1.3, label="ideal gas: c_s = sqrt(2)")
    re_, rc, rs = SOS.roman_2002_reference_points(radius=0.5, height=10.0); m = re_ <= 0.16
    ax.errorbar(re_[m], rc[m], yerr=rs[m], fmt="s", color=ROMC, ms=5, capsize=2, zorder=6, label="Roman et al. 2002 (N = 100)")
    for N in NS:
        d = df[df.N == N].sort_values("eta")
        ax.plot(d.eta, d.E1_cs, "-", color=COL[N], lw=1, zorder=4)
        for r in d.itertuples(): pt(ax, r.eta, r.E1_cs, r.E1_sigma, N, r.stationarity_flag)
    ax.set_ylabel("c_s  [sqrt(kT/m)], E1 +- sigma", fontsize=11); ax.grid(alpha=0.25); ax.set_xlim(0, 0.16)
    hh, ll = ax.get_legend_handles_labels(); size_legend(ax, extra=hh, loc="upper left", bbox_to_anchor=(0.14, 0.86))
    ax.set_title("The ideal-gas zoom, eta <= 0.15: block D (eta 0.02, 0.035, 0.075; N 100 and 400) and block B (eta 0.05, 0.10, 0.15; four sizes)", fontsize=11)
    ar.axhline(0, color=KRC, lw=2)
    for N in NS:
        d = df[df.N == N].sort_values("eta"); ar.plot(d.eta, d.Delta_pct, "-", color=COL[N], lw=0.8)
        for r in d.itertuples(): pt(ar, r.eta, r.Delta_pct, r.Delta_sigma_pct, N, r.stationarity_flag)
    ar.set_ylabel("Delta = c_s / c_s^KR - 1  [%]", fontsize=11); ar.set_xlabel("packing fraction  eta", fontsize=11); ar.grid(alpha=0.25)
    footer(fig, f"{n}.csv", root, y=-0.05)
    return [n + ".csv"] + save(fig, out, n)


# ---------------------------------------------------------------------------------------------------------------- fig4, Z
def fig_z(res, out, root, tag):
    rows, fits = [], []
    for N in NS:
        for e in V.GEN.ETA_B:
            c = res["cell"].get(cid("B", N, e))
            if not c: continue
            zk = V.z_kr(e); Z = c["Z"]
            rows.append(dict(cell=c["cid"], N=N, eta=e, Z=Z[0], Z_se=Z[1], n=Z[2], Z_KR=zk, DeltaZ_pct=100 * (Z[0] / zk - 1),
                             DeltaZ_se_pct=100 * Z[1] / zk, stationarity_flag=int(stat_flag(c))))
    for e in V.GEN.ETA_B:
        F = res["fitsZ"][e].get("F2")
        fits.append(dict(eta=e, set="Phi" if e in V.PHI else "approach", aZ_pct=100 * F["beta"][0] if F else nan, aZ_se_pct=100 * F["se"][0] if F else nan,
                         bZ_pct=100 * F["beta"][1] if F else nan, dZ_pct=100 * F["beta"][2] if F else nan, chi2=F["chi2"] if F else nan,
                         abs_aZ_lt_0p3pct=("yes" if abs(F["beta"][0]) < 0.003 else "no") if F else "-"))
    n = f"{tag}_s1_fig4_Z_vs_eta_and_DeltaZ_N100-1600"
    write_csv(os.path.join(out, n + ".csv"), rows); write_csv(os.path.join(out, n + "_fits.csv"), fits)
    df = pd.DataFrame(rows); ff = pd.DataFrame(fits)
    fig, (ax, ar) = plt.subplots(2, 1, figsize=(13, 10.5), sharex=True, gridspec_kw=dict(height_ratios=[2, 1.5], hspace=0.05))
    shade(ax); shade(ar, labels=False)
    ef = np.linspace(0.002, 0.70, 600); ek = np.linspace(0.002, 0.69, 600)
    ax.plot(ef, SOS.Z_spt_eos(ef), "--", color=SPTC, lw=1.4, label="SPT")
    ax.plot(ef, SOS.Z_henderson_eos(ef), "-.", color=HENC, lw=1.4, label="Henderson (a = 0.125)")
    ax.plot(ek, SOS.Z_kolafa_rottner_2006(ek), "-", color=KRC, lw=2.2, alpha=0.9, label="KR (drawn to eta 0.69)")
    for N in NS:
        d = df[df.N == N].sort_values("eta"); ax.plot(d.eta, d.Z, "-", color=COL[N], lw=1, zorder=4)
        for r in d.itertuples(): pt(ax, r.eta, r.Z, r.Z_se, N, r.stationarity_flag)
    ax.set_ylabel("Z (virial, last quarter of the hold) +- SE", fontsize=11); ax.grid(alpha=0.25); ax.set_xlim(0, 0.72)
    hh, ll = ax.get_legend_handles_labels(); size_legend(ax, extra=hh, loc="upper left", bbox_to_anchor=(0.035, 0.88))
    ax.set_title("The compressibility factor Z against packing fraction, four box sizes, and Delta^Z = Z / Z_KR - 1 with a_Z of F2\n"
                 "(sec. 4.8.7.3: reported, no verdict; the stated expectation |a_Z| < 0.3 % on Phi)", fontsize=11)
    ar.axhline(0, color=KRC, lw=2); ar.axhspan(-0.3, 0.3, color="0.5", alpha=0.18, lw=0, label="|a_Z| < 0.3 % (the stated expectation)")
    for N in NS:
        d = df[df.N == N].sort_values("eta"); ar.plot(d.eta, d.DeltaZ_pct, "-", color=COL[N], lw=0.8)
        for r in d.itertuples(): pt(ar, r.eta, r.DeltaZ_pct, r.DeltaZ_se_pct, N, r.stationarity_flag)
    fp = ff[ff.aZ_pct.notna()]
    ar.errorbar(fp.eta, fp.aZ_pct, yerr=fp.aZ_se_pct, fmt="*", color="k", ms=10, capsize=3, zorder=7, label="a_Z(eta), F2 (N -> infinity)")
    symlog_y(ar, 0.5, 1.0, [-10, -5, -2, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 2, 5, 10]); ar.set_ylim(-20, 15)
    ar.text(0.006, 0.5, "linear within +-0.5 %, logarithmic beyond", transform=ar.transAxes, rotation=90, va="center", fontsize=7, color="0.35")
    ar.set_ylabel("Delta^Z = Z / Z_KR - 1  [%]", fontsize=11); ar.set_xlabel("packing fraction  eta", fontsize=11); ar.grid(alpha=0.25)
    ar.legend(loc="lower left", fontsize=8, framealpha=0.95, bbox_to_anchor=(0.03, 0.0))
    footer(fig, f"{n}.csv (cells), {n}_fits.csv (a_Z, b_Z, d_Z, chi2 of F2)", root, y=-0.04)
    return [n + ".csv", n + "_fits.csv"] + save(fig, out, n)


# ---------------------------------------------------------------------------------------------------------------- fig5, delta_0.5
def fig_d05(res, out, root, tag):
    rows = []
    for N in NS:
        for e in V.GEN.ETA_B:
            c = res["cell"].get(cid("B", N, e))
            if not c: continue
            d = c.get("delta05")
            rows.append(dict(cell=c["cid"], N=N, eta=e, delta05_pct=100 * d[0] if d else nan, delta05_sigma_pct=100 * d[1] if d else nan,
                             stationarity_flag=int(stat_flag(c))))
    n = f"{tag}_s1_fig5_delta05_vs_eta_N100-1600"; write_csv(os.path.join(out, n + ".csv"), rows); df = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(13, 6.5)); shade(ax)
    ax.axhline(0, color="0.3", lw=1.2)
    for N in NS:
        d = df[(df.N == N) & df.delta05_pct.notna()].sort_values("eta")
        ax.plot(d.eta + OFF[N], d.delta05_pct, "-", color=COL[N], lw=0.8)
        for r in d.itertuples(): pt(ax, r.eta + OFF[N], r.delta05_pct, r.delta05_sigma_pct, N, r.stationarity_flag)
    ax.set_xlim(0, 0.72); ax.grid(alpha=0.25); size_legend(ax, loc="upper left", bbox_to_anchor=(0.035, 0.86))
    ax.set_xlabel("packing fraction  eta  (the four sizes drawn offset by -0.0045, -0.0015, +0.0015, +0.0045 for legibility)", fontsize=10)
    ax.set_ylabel("delta_0.5 = c_(alpha = 0.5) / c_s - 1  [%]", fontsize=11)
    ax.set_title("The light-divider term delta_0.5 per cell (sec. 4.8.6.2; Q5): the lightest divider alpha = 0.5 against E1 of the cell", fontsize=11)
    footer(fig, f"{n}.csv", root, y=-0.09)
    return [n + ".csv"] + save(fig, out, n)


# ---------------------------------------------------------------------------------------------------------------- fig6, Gamma and Q
def fig_gq(res, out, root, tag):
    rows = []
    for (c_id, M), r in res["mass"].items():
        c = res["cell"].get(c_id)
        if not c or c["blk"] != "B": continue
        rows.append(dict(cell=c_id, N=c["N"], eta=c["eta"], M=M, alpha=r["alpha"], used=r["used"], group_rule_flag=int(r["flagged"]),
                         nud=r.get("nud", nan), tau_r=r.get("tau_r", nan), Gamma=r.get("Gamma", nan), Q=r.get("Q", nan), nu0=r.get("nu0", nan),
                         stationarity_flag=int(stat_flag(c))))
    n = f"{tag}_s1_fig6_Gamma_and_Q_per_mass_vs_eta_N100-1600"; write_csv(os.path.join(out, n + ".csv"), rows); df = pd.DataFrame(rows)
    al = sorted(df.alpha.unique()); cm = plt.get_cmap("viridis"); acol = {a: cm(i / (len(al) - 1)) for i, a in enumerate(al)}
    fig, axs = plt.subplots(2, 4, figsize=(18, 9.5), sharex=True, sharey="row")
    for j, N in enumerate(NS):
        for i, (q, lab) in enumerate((("Gamma", "Gamma = 2 / tau_r  [1 / sigma-time]"), ("Q", "Q_k = pi nu_d-bar tau_r"))):
            ax = axs[i, j]; shade(ax, labels=(i == 0 and j == 0))
            for a in al:
                d = df[(df.N == N) & (df.alpha == a) & df[q].notna()].sort_values("eta")
                ax.plot(d.eta, d[q], "-", color=acol[a], lw=1)
                for r in d.itertuples():
                    o = r.group_rule_flag or r.stationarity_flag
                    ax.plot([r.eta], [getattr(r, q)], "o", color=acol[a], mfc="white" if o else acol[a], ms=4.5)
            ax.set_yscale("log"); ax.grid(alpha=0.25, which="both"); ax.set_xlim(0, 0.72)
            if q == "Q": ax.axhline(5, color=KRC, ls="--", lw=1.2)
            if j == 0: ax.set_ylabel(lab, fontsize=10)
            if i == 0: ax.set_title(f"N = {N}", fontsize=11)
            if i == 1: ax.set_xlabel("packing fraction  eta", fontsize=10)
    hl = [plt.Line2D([], [], color=acol[a], marker="o", ms=5, label=f"alpha = {a:g}") for a in al]
    hl += [plt.Line2D([], [], color="0.3", marker="o", mfc="white", ls="none", ms=5, label="open: left out by the group rule, or a stationarity flag"),
           plt.Line2D([], [], color=KRC, ls="--", label="Q = 5 (the low-Q flag)")]
    fig.legend(handles=hl, loc="upper center", ncol=6, fontsize=9, bbox_to_anchor=(0.5, 0.99))
    fig.suptitle("Damping per (cell, mass), the inputs for Paper 2 (sec. 4.8.6.4; methods sec. 13): Gamma and Q from the item-4a fit of the "
                 "seed-averaged ACF, block B", fontsize=11.5, y=1.025)
    fig.tight_layout(rect=(0, 0.04, 1, 0.93))
    footer(fig, f"{n}.csv", root, y=-0.02)
    return [n + ".csv"] + save(fig, out, n)


# ---------------------------------------------------------------------------------------------------------------- fig7, psi6
def traj_psi6(root, c, masses, res):
    """the cell's per-trajectory psi6_rec, selected as s1_verdict_261012.per_cell selects them (not flagged; clean runs of red_nu.csv)"""
    v = []
    for M in sorted(masses):
        r = res["mass"].get((c["cid"], M)); d = os.path.join(root, c["g"], "merged", "stage1", "gen3", c["cid"], f"m_{M}")
        p1, pr, ph = (os.path.join(d, f) for f in ("s1_runs.csv", "red_nu.csv", "runs_health.csv"))
        if r is None or r["flagged"] or not os.path.exists(p1) or not os.path.exists(pr): continue
        s1 = pd.read_csv(p1)
        if not len(s1): continue
        ok = set(pd.read_csv(pr)["run"].astype(int))
        if os.path.exists(ph):
            h = pd.read_csv(ph)
            if len(h): ok &= {int(x.run) for x in h.itertuples() if int(x.gen3_clean) == 1}
        v += list(s1[s1["run"].isin(sorted(ok))]["psi6_rec"])
    return v


def fig_psi6(res, out, root, tag):
    C = V.registered(); rows, traj = [], []; worst = 0.0
    for N in NS:
        for e in V.GEN.ETA_B:
            c = res["cell"].get(cid("B", N, e))
            if not c: continue
            v = traj_psi6(root, dict(cid=c["cid"], g=c["g"]), C[c["cid"]]["masses"], res)
            if v: worst = max(worst, abs(float(np.mean(v)) - c["psi6"][0]))
            rows.append(dict(cell=c["cid"], N=N, eta=e, psi6=c["psi6"][0], psi6_se=c["psi6"][1], n=c["psi6"][2], stationarity_flag=int(stat_flag(c))))
            traj += [dict(cell=c["cid"], N=N, eta=e, psi6_rec=x) for x in v]
    n = f"{tag}_s1_fig7_psi6_vs_eta_N100-1600"
    write_csv(os.path.join(out, n + ".csv"), rows); write_csv(os.path.join(out, n + "_trajectories.csv"), traj)
    df = pd.DataFrame(rows); dt = pd.DataFrame(traj)
    fig, (ax, az) = plt.subplots(1, 2, figsize=(17, 7.5), gridspec_kw=dict(width_ratios=[1.6, 1], wspace=0.16))
    shade(ax); shade(az, x_min=0.575, labels=False)
    for N in NS:
        d = df[df.N == N].sort_values("eta"); ax.plot(d.eta + OFF[N], d.psi6, "-", color=COL[N], lw=1.1, zorder=5)
        for r in d.itertuples(): pt(ax, r.eta + OFF[N], r.psi6, r.psi6_se, N, r.stationarity_flag, ms=7, zorder=6)
    ax.set_yscale("log"); ax.set_xlim(0, 0.72); ax.set_ylim(0.015, 1.0); ax.grid(alpha=0.25, which="both")
    tk = [0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0]; ax.set_yticks(tk); ax.set_yticklabels([f"{t:g}" for t in tk])
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel("packing fraction  eta  (sizes offset by -0.0045, -0.0015, +0.0015, +0.0045)", fontsize=10)
    ax.set_ylabel("psi6: the mean of the two compartments' record means", fontsize=10.5)
    size_legend(ax, loc="upper left", bbox_to_anchor=(0.035, 0.86))
    ax.set_title("per cell: the mean +- SE over its trajectories (all masses)", fontsize=10.5)
    rng = np.random.default_rng(0)                                   # horizontal jitter of the dots only (drawing)
    for N in NS:
        t = dt[(dt.N == N) & (dt.eta > 0.575)]
        az.plot(t.eta + OFF[N] + rng.uniform(-0.0012, 0.0012, len(t)), t.psi6_rec, ".", color=COL[N], ms=3, alpha=0.35, zorder=3)
        d = df[(df.N == N) & (df.eta > 0.575)].sort_values("eta"); az.plot(d.eta + OFF[N], d.psi6, "-", color=COL[N], lw=1.1, zorder=5)
        for r in d.itertuples(): pt(az, r.eta + OFF[N], r.psi6, r.psi6_se, N, r.stationarity_flag, ms=7, zorder=6)
    az.set_xlim(0.575, 0.70); az.set_ylim(0, 0.75); az.grid(alpha=0.25)
    az.set_xlabel("packing fraction  eta  (0.60-0.69)", fontsize=10); az.set_ylabel("psi6", fontsize=10.5)
    az.set_title("0.60-0.69, linear: every trajectory's value as a dot (jittered in eta)", fontsize=10.5)
    fig.suptitle("Bond order psi6 against packing fraction, block B (sec. 4.8.6.4)", fontsize=12, y=0.995)
    footer(fig, f"{n}.csv (cells), {n}_trajectories.csv (the dots)", root, y=-0.075)
    return [n + ".csv", n + "_trajectories.csv"] + save(fig, out, n), worst


# ---------------------------------------------------------------------------------------------------------------- fig8, Verdict S
def fig_s(res, out, root, tag):
    rows = []
    for (N, e), v in sorted(res["S"].items()):
        z = v["z"]; rows.append(dict(N=N, eta=e, z_Z=z["Z"], z_psi6=z["psi6"], z_nud_alpha1=z["nud_a1"], z_nud_alpha10=z["nud_a10"],
                                    verdict="no data" if v["independent"] is None else ("START-INDEPENDENT" if v["independent"] else "START-DEPENDENT")))
    n = f"{tag}_s1_fig8_verdictS_four_z_per_cell"; write_csv(os.path.join(out, n + ".csv"), rows); df = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(14, 6.2)); x = np.arange(len(df))
    for k, (q, lab, m, col) in enumerate((("z_Z", "z of Z", "o", "#0072b2"), ("z_psi6", "z of psi6", "s", "#e69f00"),
                                          ("z_nud_alpha1", "z of nu_d-bar, alpha = 1", "^", "#009e73"), ("z_nud_alpha10", "z of nu_d-bar, alpha = 10", "D", "#cc79a7"))):
        ax.plot(x + (k - 1.5) * 0.14, df[q], m, color=col, ms=7, label=lab, ls="none")
    for s in (-1, 1): ax.axhline(s * V.Z_S, color=KRC, ls="--", lw=1.3)
    ax.axhline(0, color="0.4", lw=0.8)
    ax.set_xticks(x); ax.set_xticklabels([f"N {r.N}\n{r.eta:.2f}" for r in df.itertuples()], fontsize=8)
    for b in range(1, len(NS)): ax.axvline(b * 5 - 0.5, color="0.7", lw=0.8)
    ax.set_ylabel("z = (mean_0 - mean_turned) / sqrt(SE_0^2 + SE_turned^2)", fontsize=10); ax.grid(alpha=0.25, axis="y")
    nind = int((df.verdict == "START-INDEPENDENT").sum())
    ax.set_title(f"Verdict S (sec. 4.8.7.5), block S: the lattice start 0 against turned, per cell (N, eta): the four z\n"
                 f"START-INDEPENDENT if all four |z| < {V.Z_S:g} (red dashed): {nind} of {len(df)} cells", fontsize=11)
    ax.legend(loc="upper left", fontsize=8.5, ncol=4, framealpha=0.95)
    lim = max(4.0, float(np.nanmax(np.abs(df[["z_Z", "z_psi6", "z_nud_alpha1", "z_nud_alpha10"]].to_numpy()))) + 0.5); ax.set_ylim(-lim, lim)
    footer(fig, f"{n}.csv", root, y=-0.13)
    return [n + ".csv"] + save(fig, out, n)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--data", required=True); ap.add_argument("--out", required=True); ap.add_argument("--date", default="261015")
    a = ap.parse_args()
    if os.path.exists(a.out): sys.exit(f"REFUSED: {a.out} exists (a new folder only)")
    os.makedirs(a.out)
    res = V.analyse(a.data, interim=False, jackknife=True, quiet=True)
    print(f"# the registered figures (s1_figures_261015.py), from {a.data}: groups {', '.join(res['groups_present'])}; stopped: {res['stopped'] or 'none'}")
    kk = max(abs(V.kr_cs(e) / float(cs_eos(SOS.Z_kolafa_rottner_2006, SOS.dZ_kolafa_rottner_2006, np.array([e]))[0]) - 1) for e in V.GEN.ETA_B)
    print(f"- check: KR's c_s curve (analytic Z') against the verdict's kr_cs (numerical Z') at the 16 densities: largest relative difference {kk:.1e}")
    files = []
    f, _ = fig_cs(res, a.out, a.data, a.date); files += f
    files += fig_zoom(res, a.out, a.data, a.date)
    files += fig_z(res, a.out, a.data, a.date)
    files += fig_d05(res, a.out, a.data, a.date)
    files += fig_gq(res, a.out, a.data, a.date)
    f, worst = fig_psi6(res, a.out, a.data, a.date); files += f
    print(f"- check: the per-trajectory psi6 dots against the verdict's cell means: largest |mean of the dots - psi6 of the cell| {worst:.1e}")
    files += fig_s(res, a.out, a.data, a.date)
    print("- (sec. 4.8.11 item 5, the identity rho_I against 1/N_s and R_C(eta), waits for Part A: not drawn)")
    print(f"- written into {a.out}:"); [print(f"  {x}") for x in files]


if __name__ == "__main__":
    main()
