#!/usr/bin/env python3
"""##CHRIS 2026-10-15 (plan-author decision 56 task 1; decision 55 addendum 1 item 2; notes 261012 sec. 4.7.129, 4.7.133): FIGURE 1 WITH
THE DILUTE BLOCK. c_s(eta) of E1 over 0.02-0.69 for N 100, 400, 900, 1600, blocks B and D (D: eta 0.02, 0.035, 0.075 at N 100 and 400,
the same markers), with the 2D ideal-gas limit c_s = sqrt(2) (gamma = 2, kT = m = 1) dashed, KR as the reference curve from eta -> 0 to
0.69, SPT and Henderson thin, the phase regions shaded, and below it Delta = c_s / c_s^KR - 1 [%] with the extrapolated a(eta) of F2 and
F1 (block B; block D has two sizes, no fit). Linear eta axis. Figure 3 stays as the dilute zoom.
The numbers are the registered full report's (Report D54-1), read from the CSVs written with them (full_261015/figures_261015b/: fig1's
block B cells, fig3's block D cells, fig2's fits); nothing is fitted or estimated here. The CSV written here holds B's and D's rows.
The phase regions (Bernard and Krauth 2011; Engel et al. 2013; the project's constants in plot_speed_of_sound_edmd): fluid below
0.700, liquid-hexatic coexistence 0.700-0.716, hexatic 0.716-0.720, solid above 0.720.
usage (from hspist3/): python3 cluster/gen3_stage1_261012/s1_figures_261015c.py --out <new folder> [--date 261015]
"""
import argparse, csv, math, os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE)); sys.path.insert(0, HS)
import plot_speed_of_sound_edmd as SOS
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker

B = os.path.join(HS, "experiments_gen3_s1_261012", "full_261015", "figures_261015b")
F1 = os.path.join(B, "261015_s1_fig1_cs_vs_eta_E1_N100-1600.csv"); F3 = os.path.join(B, "261015_s1_fig3_cs_ideal_gas_zoom_eta_le_0.15_B_D.csv")
F2 = os.path.join(B, "261015_s1_fig2_Delta_vs_Nm12_per_eta_F1_F2_fits.csv")
SCRIPT = "cluster/gen3_stage1_261012/s1_figures_261015c.py"
NS = (100, 400, 900, 1600)
COL = {100: "#6baed6", 400: "#2b7bba", 900: "#0b4f8f", 1600: "#061f4a"}
MK = {100: "o", 400: "s", 900: "^", 1600: "D"}
KRC, SPTC, HENC, IDC = "#d62728", "#e69f00", "#009e73", "0.25"
GEOM = ("N = 100, 400, 900, 1600 disks, N/2 per side; disk radius 0.5; divider thickness 0.05; H = 10 sqrt(N/100) = 10, 20, 30, 40; "
        "L0 = N pi / (8 H eta) per side (exact box); gen3 engine, one build (f72d764); 25 seeds per (cell, mass)")
REGIONS = [(0.0, SOS.ETA_COEX_LO, "fluid", "#eef4fb"), (SOS.ETA_COEX_LO, SOS.ETA_COEX_HI, "liquid-hexatic\ncoexistence", "#ffd9b8"),
           (SOS.ETA_COEX_HI, SOS.ETA_HEXATIC_SOLID, "hexatic", "#f6c7df"), (SOS.ETA_HEXATIC_SOLID, 0.74, "solid", "#ded8f5")]


def cs_eos(Zf, dZf, e): return SOS.cs_adiabatic_2d_monatomic(Zf(e), dZf(e), e, kbt=1.0, m=1.0)


def shade(ax, labels=True):
    for lo, hi, lab, col in REGIONS:
        ax.axvspan(lo, hi, color=col, alpha=0.55, lw=0, zorder=0)
        if labels:
            if hi - lo > 0.03: ax.text((lo + hi) / 2, 0.975, lab, transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=9, color="0.3")
            else: ax.text((lo + hi) / 2, 0.70 if hi <= SOS.ETA_COEX_HI + 1e-9 else 0.33, lab, transform=ax.get_xaxis_transform(), ha="center",
                          va="center", rotation=90, fontsize=7, color="0.3")
    for x in (SOS.ETA_COEX_LO, SOS.ETA_COEX_HI, SOS.ETA_HEXATIC_SOLID): ax.axvline(x, color="0.45", lw=0.7, ls=":", zorder=1)


def pt(ax, x, y, s, N, flagged, ms=6):
    ax.errorbar([x], [y], yerr=[s], fmt=MK[N], color=COL[N], mfc="white" if flagged else COL[N], mec=COL[N], ms=ms, capsize=2, lw=1, zorder=5)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out", required=True); ap.add_argument("--date", default="261015"); a = ap.parse_args()
    if os.path.exists(a.out) and set(os.listdir(a.out)) - {"RUN_RECORD.txt"}:
        sys.exit(f"REFUSED: {a.out} holds files (a new folder only; its RUN_RECORD.txt alone may be there)")
    b = pd.read_csv(F1); z = pd.read_csv(F3); f = pd.read_csv(F2)
    d = z[z.block == "D"].copy()
    # the B rows of fig3 (eta <= 0.15) must be fig1's
    zb = z[z.block == "B"].merge(b, on=["cell", "N", "eta"], suffixes=("_3", "_1"))
    worst = float(np.max(np.abs(zb.E1_cs_3 - zb.E1_cs_1))) if len(zb) else 0.0
    rows = []
    for r in b.itertuples():
        rows.append(dict(block="B", cell=r.cell, N=r.N, eta=r.eta, E1_cs=r.E1_cs, E1_sigma=r.E1_sigma, E1_masses=r.E1_masses, KR=r.KR,
                         Delta_pct=r.Delta_pct, Delta_sigma_pct=r.Delta_sigma_pct, vs_sqrt2_pct=100 * (r.E1_cs / math.sqrt(2) - 1), stationarity_flag=r.stationarity_flag))
    for r in d.itertuples():
        rows.append(dict(block="D", cell=r.cell, N=r.N, eta=r.eta, E1_cs=r.E1_cs, E1_sigma=r.E1_sigma, E1_masses=r.E1_masses, KR=r.KR,
                         Delta_pct=r.Delta_pct, Delta_sigma_pct=r.Delta_sigma_pct, vs_sqrt2_pct=r.vs_sqrt2_pct, stationarity_flag=r.stationarity_flag))
    df = pd.DataFrame(rows).sort_values(["N", "eta"]); os.makedirs(a.out, exist_ok=True)
    n = f"{a.date}_s1_fig1c_cs_vs_eta_E1_blocks_B_D_N100-1600"
    df.to_csv(os.path.join(a.out, n + ".csv"), index=False)
    fig, (ax, ar) = plt.subplots(2, 1, figsize=(13, 11), sharex=True, gridspec_kw=dict(height_ratios=[2.1, 1.5], hspace=0.05))
    shade(ax); shade(ar, labels=False)
    e_all = np.linspace(0.0005, 0.70, 700); ek = np.linspace(0.0005, 0.69, 700)
    ax.plot(e_all, cs_eos(SOS.Z_spt_eos, SOS.dZ_spt_eos, e_all), "--", color=SPTC, lw=1.0, label="SPT")
    ax.plot(e_all, cs_eos(SOS.Z_henderson_eos, SOS.dZ_henderson_eos, e_all), "-.", color=HENC, lw=1.0, label="Henderson (a = 0.125)")
    ax.plot(ek, cs_eos(SOS.Z_kolafa_rottner_2006, SOS.dZ_kolafa_rottner_2006, ek), "-", color=KRC, lw=2.4, alpha=0.9,
            label="KR (Kolafa-Rottner 2006, the rho_max 0.90 fit), eta -> 0 to 0.69")
    ax.axhline(math.sqrt(2), color=IDC, ls=(0, (6, 3)), lw=1.4, label="ideal gas: c_s = sqrt(2) (2D, gamma = 2, kT = m = 1)")
    for N in NS:
        s = df[df.N == N].sort_values("eta")
        ax.plot(s.eta, s.E1_cs, "-", color=COL[N], lw=1, zorder=4)
        for r in s.itertuples(): pt(ax, r.eta, r.E1_cs, r.E1_sigma, N, r.stationarity_flag)
    ax.set_xlim(0, 0.74); ax.set_ylim(0, max(df.E1_cs.max() * 1.08, 16)); ax.grid(alpha=0.25)
    ax.set_ylabel("speed of sound  c_s  [sqrt(kT/m)], E1 +- sigma", fontsize=11)
    h = [plt.Line2D([], [], marker=MK[N], color=COL[N], ls="-", lw=1, ms=6, label=f"N = {N}" + (" (blocks B and D)" if N in (100, 400) else " (block B)")) for N in NS]
    h += [plt.Line2D([], [], marker="o", color="0.3", mfc="white", ls="none", ms=6, label="open: stationarity flag (none)")]
    hh, ll = ax.get_legend_handles_labels()
    ax.legend(handles=h + hh, loc="upper left", bbox_to_anchor=(0.01, 0.9), fontsize=8.5, framealpha=0.95)
    ax.set_title("Stage 1, Part B: speed of sound against packing fraction, E1, blocks B (eta 0.05-0.69) and D (eta 0.02, 0.035, 0.075)\n"
                 "Verdict K on Phi (eta <= 0.60): KR-INCONSISTENT (X^2 = 32.69, 12 dof, p = 0.0011, max |z| = 3.58 at 0.60); 0.63-0.69 approach cells, no verdict",
                 fontsize=10.5)
    # the inset: eta <= 0.16, where the dilute block sits
    ins = ax.inset_axes([0.475, 0.47, 0.29, 0.40]); shade(ins, labels=False)
    e_in = np.linspace(0.0005, 0.16, 300)
    ins.plot(e_in, cs_eos(SOS.Z_kolafa_rottner_2006, SOS.dZ_kolafa_rottner_2006, e_in), "-", color=KRC, lw=1.8)
    ins.plot(e_in, cs_eos(SOS.Z_spt_eos, SOS.dZ_spt_eos, e_in), "--", color=SPTC, lw=0.8)
    ins.plot(e_in, cs_eos(SOS.Z_henderson_eos, SOS.dZ_henderson_eos, e_in), "-.", color=HENC, lw=0.8)
    ins.axhline(math.sqrt(2), color=IDC, ls=(0, (6, 3)), lw=1.2)
    for N in NS:
        s = df[(df.N == N) & (df.eta <= 0.16)].sort_values("eta"); ins.plot(s.eta, s.E1_cs, "-", color=COL[N], lw=0.8)
        for r in s.itertuples(): pt(ins, r.eta, r.E1_cs, r.E1_sigma, N, r.stationarity_flag, ms=4)
    ins.set_xlim(0, 0.16); ins.set_ylim(1.38, 2.02); ins.tick_params(labelsize=7); ins.grid(alpha=0.25)
    ins.set_title("eta <= 0.16 (figure 3 is the full zoom)", fontsize=7.5)
    # the residual, symmetric-log as figure 1 of figures_261015b
    ar.axhline(0, color=KRC, lw=2, alpha=0.9)
    for N in NS:
        s = df[df.N == N].sort_values("eta"); ar.plot(s.eta, s.Delta_pct, "-", color=COL[N], lw=0.8, zorder=4)
        for r in s.itertuples(): pt(ar, r.eta, r.Delta_pct, r.Delta_sigma_pct, N, r.stationarity_flag)
    fp = f[f.F2_b0_pct.notna()]
    ar.errorbar(fp.eta - 0.002, fp.F2_b0_pct, yerr=fp.F2_se0_pct, fmt="*", color="k", ms=10, capsize=3, zorder=7, label="a(eta), F2 (N -> infinity) +- sigma_a, block B")
    f1 = f[f.F1_b0_pct.notna()]
    ar.errorbar(f1.eta + 0.002, f1.F1_b0_pct, yerr=f1.F1_se0_pct, fmt="*", color="0.5", mfc="white", ms=9, capsize=2, zorder=6, label="a(eta), F1 (reported), block B")
    ar.set_yscale("symlog", linthresh=2.5, linscale=1.8); tk = [-20, -10, -5, -2, -1, 0, 1, 2, 5, 10, 20, 40]
    ar.set_yticks(tk); ar.set_yticklabels([f"{t:g}" for t in tk]); ar.yaxis.set_minor_locator(matplotlib.ticker.NullLocator()); ar.set_ylim(-30, 55)
    ar.axhspan(-2.5, 2.5, color="0.5", alpha=0.06, lw=0, zorder=0)
    ar.text(0.006, 0.5, "linear within +-2.5 %, logarithmic beyond", transform=ar.transAxes, rotation=90, va="center", fontsize=7, color="0.35")
    ar.set_ylabel("Delta = c_s / c_s^KR - 1  [%]", fontsize=11); ar.set_xlabel("packing fraction  eta  (linear)", fontsize=11); ar.grid(alpha=0.25)
    ar.legend(loc="lower left", fontsize=8, framealpha=0.95, bbox_to_anchor=(0.03, 0.0))
    fig.text(0.5, 0.052, f"data: {n}.csv (blocks B and D, from figures_261015b's fig1 and fig3 CSVs), figures_261015b/261015_s1_fig2_Delta_vs_Nm12_per_eta_F1_F2_fits.csv "
             f"(a of F2 and F1)\nscript: {SCRIPT} (the full report's numbers, Report D54-1; nothing fitted here); the phase regions: Bernard and Krauth 2011, "
             f"Engel et al. 2013\n{GEOM}", ha="center", va="top", fontsize=7.2, color="0.25")
    for ext in ("png", "pdf"): fig.savefig(os.path.join(a.out, f"{n}.{ext}"), dpi=170 if ext == "png" else None, bbox_inches="tight")
    print(f"# figure 1 with block D (s1_figures_261015c.py)\n- inputs: {F1}, {F3}, {F2}\n- check: fig3's block B rows against fig1's: largest |c_s difference| {worst:.1e} ({len(zb)} rows)"
          f"\n- rows: {sum(df.block == 'B')} of block B, {sum(df.block == 'D')} of block D\n- written into {a.out}: {n}.csv, {n}.png, {n}.pdf")


if __name__ == "__main__":
    main()
