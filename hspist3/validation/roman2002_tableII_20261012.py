#!/usr/bin/env python3
"""##CHRIS 2026-10-12: Roman 2002 Table II (system-size series) re-mapped against Kolafa-Rottner.

SOURCE, read from the PDF: Roman, Gonzalez, White & Velasco, Am. J. Phys. 70, 847 (2002),
doi 10.1119/1.1482060 (NOT Eur. J. Phys.). Table II, p. 851: "Influence of the size of the system on
the speed of sound of a hard disk gas with packing fraction eta = 0.393":

    N     L_0       c_s
    64    8 sqrt2   3.81
    256   16 sqrt2  3.78
    1024  32 sqrt2  3.76
    4096  64 sqrt2  3.75

Text, p. 851: same aspect ratio (A = L_0) at fixed eta = 0.393; constant K = 1.07687, i.e. piston mass
M = N m; "the extrapolated value of the velocity of sound is c_s = 3.74". NO uncertainties are printed,
so the points are plotted WITHOUT error bars, labelled "as published", and no chi^2 is computed.
N is per compartment: eta = N pi r^2 / (A L_0) with A = L_0 gives pi/8 = 0.392699 for every row
(recomputed below, not copied).

THE L/H AMBIGUITY. Roman grows A and L_0 together (A = L_0 ~ sqrt N), so N^{-1/2}, 1/L_0 and 1/H are
the same variable along this series: a finite-size term ~1/H and one ~1/L cannot be told apart here.

OUR POINTS (blue). N_s = 50, H = L_0 = 10: the N = 50 member of the same square series.
  * canonical: the 9-mass A1v2 value from the 260919 table, with its scaled error;
  * like-for-like: the A1v2 M = 50 cell alone, alpha = M/(2 N_s m) = 0.5 -- Roman's own K.
Both are recomputed with the canonical estimator (paper1_populate_cs_err_20261002.cell) from T.DROOT,
and the 9-mass recomputation must reproduce the table to 5e-6 or our points are marked VOID.
"""
import csv, math, os, sys
import numpy as np
from scipy.optimize import brentq
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import tests_20260913 as T
import plot_speed_of_sound_edmd as sos
from paper1_populate_cs_err_20261002 import cell, slope_with_errors

OUT = os.path.join(os.path.dirname(os.path.dirname(HERE)), "0000_PLAN_OVERALL",
                   "paper1_speedofsound", "experiments", "final")
R = 0.5
TABLE2 = [(64, 8, 3.81), (256, 16, 3.78), (1024, 32, 3.76), (4096, 64, 3.75)]   # (N, L_0/sqrt2, c_s)
ROMAN_EXTRAP, ROMAN_K = 3.74, 1.07687
ETA, L0, NS = 0.392699, 10.0, 50

def kr(eta):
    e = np.array([eta]); Z = sos.Z_kolafa_rottner_2006(e); dZ = sos.dZ_kolafa_rottner_2006(e)
    return float(sos.cs_adiabatic_2d_monatomic(Z, dZ, e, kbt=1.0, m=1.0)[0])

def ours():
    ref = [r for r in csv.DictReader(open(T.plot_path("260919_A1v2_final_cs_vs_eta.csv")))
           if abs(float(r["eta"]) - ETA) < 1e-6][0]
    D = {}
    for M in T.A1_MASSES:
        runs = T.cell_runs(os.path.join(T.DROOT, "eta_0p392699", f"m_{M}"), M)
        D[M] = cell((ETA, L0, M, runs))
    x = np.array([T.x_of(M, L0) for M in T.A1_MASSES]); y = np.array([D[M]["nu"] for M in T.A1_MASSES])
    sy = np.array([D[M]["sd"] / math.sqrt(D[M]["n"]) for M in T.A1_MASSES])
    s, e, es, ch = slope_with_errors(x, y, sy)
    ok = abs(s - float(ref["c_s"])) <= 5e-6 * abs(float(ref["c_s"]))
    m50 = D[50]; x50 = T.x_of(50, L0)
    return dict(ref=ref, s=s, es=es, ok=ok, cs50=m50["nu"] / x50,
                e50=m50["sd"] / math.sqrt(m50["n"]) / x50, n50=m50["n"])

def main():
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    rows = []
    for N, l, c in TABLE2:
        L = l * math.sqrt(2); eta = N * math.pi * R * R / (L * L); rows.append((N, L, eta, c))
    K = kr(math.pi / 8)
    xs = np.array([N ** -0.5 for N, *_ in rows]); ys = np.array([c for *_, c in rows])
    b, a = np.polyfit(xs, ys, 1)
    O = ours(); kk = brentq(lambda k: math.cos(k) / math.sin(k) - 0.5 * k, 1e-6, math.pi - 1e-6)
    print("### Roman 2002 Table II (Am. J. Phys. 70, 847; table on p. 851), re-mapped\n")
    print(f"KR adiabatic c_s at eta = pi/8 = {math.pi/8:.6f}: **{K:.5f}**\n")
    print("| N (per compartment) | L_0 = A | eta recomputed | c_s as published | (c_s - KR)/KR | guide line |")
    print("|---|---|---|---|---|---|")
    for (N, L, eta, c), x in zip(rows, xs):
        print(f"| {N} | {L:.4f} | {eta:.6f} | {c:.2f} | {100*(c/K-1):+.2f} % | {a+b*x:.4f} |")
    print(f"\nGuide: unweighted straight line in N^(-1/2) through the four published values (no errors exist, "
          f"so no chi^2): c_s = {a:.4f} + {b:.4f} N^(-1/2).")
    print(f"Its intercept {a:.4f} vs Roman's own extrapolated 3.74 (as published) and KR {K:.5f} "
          f"({100*(a/K-1):+.2f} %); at N = 50 the guide gives {a+b/math.sqrt(50):.4f}.\n")
    print(f"Our estimator self-check (9-mass recomputation from T.DROOT vs 260919 table {O['ref']['c_s']}): "
          f"{O['s']:.5f} -> {'REPRODUCED to 5e-6' if O['ok'] else '**NOT reproduced -- our points VOID**'}")
    print(f"Our canonical point, N_s = 50, H = L_0 = 10: c_s = {float(O['ref']['c_s']):.5f} ± "
          f"{float(O['ref']['c_s_err_scaled']):.5f} (scaled), {100*(float(O['ref']['c_s'])/K-1):+.2f} % vs KR; "
          f"minus guide at N = 50: {float(O['ref']['c_s'])-(a+b/math.sqrt(50)):+.4f}")
    print(f"Our like-for-like point, M = 50 only (alpha = 0.5, K = {kk:.5f} vs Roman's {ROMAN_K}): c_s = "
          f"{O['cs50']:.5f} ± {O['e50']:.5f} ({O['n50']} seeds), {100*(O['cs50']/K-1):+.2f} % vs KR; "
          f"minus guide at N = 50: {O['cs50']-(a+b/math.sqrt(50)):+.4f}")
    with open(os.path.join(OUT, "261012_roman2002_tableII_vs_KR.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["source", "N_per_compartment", "L0", "eta", "c_s", "c_s_err", "KR", "note"])
        for N, L, eta, c in rows: w.writerow(["Roman2002_TableII", N, f"{L:.6f}", f"{eta:.6f}", c, "", f"{K:.6f}", "as published, no uncertainty"])
        w.writerow(["Roman2002_extrapolated", "inf", "", f"{math.pi/8:.6f}", ROMAN_EXTRAP, "", f"{K:.6f}", "as published"])
        w.writerow(["ours_A1v2_9mass", NS, L0, ETA, O["ref"]["c_s"], O["ref"]["c_s_err_scaled"], f"{K:.6f}", "260919 table, scaled error"])
        w.writerow(["ours_A1v2_M50", NS, L0, ETA, f"{O['cs50']:.6f}", f"{O['e50']:.6f}", f"{K:.6f}", "alpha=0.5 like-for-like, seed SE"])
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    xx = np.linspace(0, 0.15, 50)
    ax.axhline(K, color="red", lw=1.5, label=f"Kolafa–Rottner adiabatic, {K:.4f}")
    ax.plot(xx, a + b * xx, color="black", lw=1, ls=":", label="guide: line through Table II (no fit statistics)")
    ax.plot(xs, ys, "o", color="black", ms=7, label="Román 2002 Table II (as published, no uncertainties)")
    ax.plot([0], [ROMAN_EXTRAP], "D", color="black", mfc="white", ms=7, label="Román's extrapolation 3.74 (as published)")
    x50 = 50 ** -0.5
    if O["ok"]:
        ax.errorbar([x50], [float(O["ref"]["c_s"])], yerr=[float(O["ref"]["c_s_err_scaled"])], fmt="o",
                    color="tab:blue", ms=7, capsize=3, label="this work, N_s = 50, H = L_0 = 10 (9-mass fit)")
        ax.errorbar([x50 + 0.003], [O["cs50"]], yerr=[O["e50"]], fmt="o", color="tab:blue", mfc="white",
                    ms=7, capsize=3, label="this work, M = 50 only (α = 0.5, Román's K); x offset +0.003")
    ax.set_xlabel(r"$N^{-1/2}$  ($N$ disks per compartment; square compartments, $A = L_0 \propto \sqrt{N}$)")
    ax.set_ylabel(r"$c_s$  [$(kT/m)^{1/2}$]")
    ax.set_title(r"Size series at $\eta = \pi/8$: Román 2002 Table II vs Kolafa–Rottner", fontsize=10)
    ax.set_xlim(-0.006, 0.155); ax.legend(fontsize=7.5, frameon=False, loc="upper left"); fig.tight_layout()
    for ext in ("png", "pdf"): fig.savefig(os.path.join(OUT, f"261012_roman2002_tableII_vs_KR.{ext}"), dpi=200)
    print("\nfigure + csv: 0000_PLAN_OVERALL/paper1_speedofsound/experiments/final/261012_roman2002_tableII_vs_KR.{png,pdf,csv}")

if __name__ == "__main__":
    main()
