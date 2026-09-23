#!/usr/bin/env python3
"""##CHRIS 2026-10-02: the melting-region figure for Paper 1 sec:ordering. Analysis only, no runs.

c_s(eta) over 0.66 <= eta <= 0.76 with the CORRECTED error bars (c_s_err_scaled), against the
structural landmarks measured independently by Engel, Anderson, Glotzer, Isobe, Bernard & Krauth,
Phys. Rev. E 87, 042134 (2013):

    eta ~= 0.702, 0.714   the extrema of their Mayer-Wood loop, position independent of N
    eta ~= 0.716          end of liquid-hexatic coexistence ("the region eta >~ 0.716 is thus hexatic")
    eta  = 0.720          positional order jumps; C_q0(r) reaches the KTHNY r^(-1/3) stability limit

Kolafa-Rottner is continued DASHED beyond 0.69 to make visible that it is an extrapolation there --
a fluid EOS fitted below the transition, so a deviation from it in this window measures the
extrapolation and not the gas.

The claim the figure supports: c_s has a local maximum at eta = 0.700 and a local minimum at
eta = 0.715, bracketing the same interval as the loop. Mechanism: between the loop's extrema
Z'(eta) < 0, and c_s^2 = (kT/m)(Z + eta Z' + Z^2), so a loop in Z must appear as a dip in c_s over
the same range.
"""
import csv, math, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import tests_20260913 as T
import plot_speed_of_sound_edmd as sos

BLUE, RED, GREY, ORANGE = "#2a78d6", "#e34948", "#52514e", "#eb6834"
LO, HI = 0.66, 0.765


def main():
    rows = list(csv.DictReader(open(T.plot_path("260919_A1v2_final_cs_vs_eta.csv"))))
    d = [(float(r["eta"]), float(r["c_s"]),
          float(r.get("c_s_err_scaled") or r["c_s_scatter_mass"]),
          float(r.get("c_s_err") or 0.0), float(r.get("chi2_red") or 0.0))
         for r in rows if LO <= float(r["eta"]) <= HI]
    d.sort()
    e = np.array([q[0] for q in d]); c = np.array([q[1] for q in d])
    es = np.array([q[2] for q in d]); er = np.array([q[3] for q in d])

    fig, ax = plt.subplots(figsize=(9.2, 5.8))

    # Engel landmarks
    for x, lab, col, ls in ((0.702, "Mayer–Wood\nextremum", GREY, "--"),
                            (0.714, "Mayer–Wood\nextremum", GREY, "--"),
                            (0.716, "coexistence\nends", ORANGE, "-"),
                            (0.720, "hexatic→solid", ORANGE, ":")):
        ax.axvline(x, color=col, ls=ls, lw=1.6, alpha=0.85, zorder=1)
    ax.axvspan(0.700, 0.716, color=ORANGE, alpha=0.07, zorder=0)
    ax.annotate("liquid–hexatic coexistence\n(Engel et al. 2013)", xy=(0.708, 0.055),
                xycoords=("data", "axes fraction"), ha="center", fontsize=8.5, color=ORANGE)

    # KR, dashed beyond 0.69 to mark it as extrapolation
    g = np.linspace(LO, HI, 400)
    h = 1e-5
    kr = np.array([float(sos.cs_adiabatic_2d_monatomic(
        sos.Z_kolafa_rottner_2006(x),
        (sos.Z_kolafa_rottner_2006(x + h) - sos.Z_kolafa_rottner_2006(x - h)) / (2 * h),
        x, kbt=1, m=1)) for x in g])
    ok = g <= 0.69
    ext = (g > 0.69) & (g <= 0.705)     # stop before the fit's derivative explodes -- see below
    ax.plot(g[ok], kr[ok], "-", color=RED, lw=2.2, zorder=2,
            label="Kolafa–Rottner 2006 (fitted range, $\\eta \\leq 0.69$)")
    ax.plot(g[ext], kr[ext], "--", color=RED, lw=1.8, alpha=0.8, zorder=2,
            label="Kolafa–Rottner, EXTRAPOLATED (no fluid branch here)")
    # Why the dashed line stops at 0.705 rather than running to 0.76: the KR fit does not merely
    # become inaccurate there, it becomes unusable for a SOUND SPEED. c_s depends on Z', and the
    # fitted Z' runs 48.3 (eta=0.67) -> 27.8 (0.69) -> -8.7 (0.700) -> +3386 (0.720). A negative Z'
    # inside the stated fitted range already makes KR's own c_s turn over and fall from eta ~ 0.683
    # -- unphysical. Drawing that divergence would put a spurious factor-4 feature on the axis.
    ax.annotate("KR's fitted $Z'$ changes sign near $\\eta = 0.70$\nand diverges by $0.72$: $c_s$ from it is\n"
                "meaningless here, not merely uncertain",
                xy=(0.6975, 12.3), fontsize=8.2, color=RED, ha="left", va="top")

    ax.errorbar(e, c, yerr=es, fmt="o-", color=BLUE, ms=6, lw=1.3, capsize=3.5, zorder=5,
                label="$N = 100$, 9 masses $\\times$ 25 seeds\nerror = SE of slope $\\times\\max(1,\\sqrt{\\chi^2_{\\rm red}})$")

    # mark the two turning points the claim rests on
    imax = int(np.argmin(np.abs(e - 0.700))); imin = int(np.argmin(np.abs(e - 0.715)))
    for i, lab in ((imax, f"local max\n$\\eta = {e[imax]:.3f}$"), (imin, f"local min\n$\\eta = {e[imin]:.3f}$")):
        ax.annotate(lab, xy=(e[i], c[i]), xytext=(e[i] + (0.012 if i == imax else -0.004), c[i] + 3.2),
                    fontsize=9, color=BLUE, ha="center",
                    arrowprops=dict(arrowstyle="->", color=BLUE, lw=1.3))

    ax.set_xlim(LO, HI); ax.set_ylim(11.5, 39.5)
    ax.set_xlabel(r"packing fraction  $\eta$", fontsize=11.5)
    ax.set_ylabel(r"speed of sound  $c_s$  [$\sqrt{k_BT/m}$]", fontsize=11.5)
    ax.set_title("The melting region: $c_s$ dips across the coexistence interval", fontsize=12.5)
    ax.grid(alpha=0.3); ax.legend(frameon=False, fontsize=8.5, loc="upper left")

    out = os.path.join(T.PLOTS, "261002_p1_melting_region")
    for _ext in ("png", "pdf"):          # not `ext` -- that name is the KR mask above
        fig.savefig(f"{out}.{_ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote 261002_p1_melting_region.png/.pdf")

    print("\nthe dip, in the two error bars:")
    print("| eta | c_s | c_s_err (propagated only) | c_s_err_scaled (plotted) | chi2_red |")
    print("|---|---|---|---|---|")
    for q in d:
        print(f"| {q[0]:.4f} | {q[1]:8.4f} | {q[3]:.5f} | {q[2]:.5f} | {q[4]:6.2f} |")
    dip = c[imax] - c[imin]
    sg = math.sqrt(es[imax] ** 2 + es[imin] ** 2)
    sr = math.sqrt(er[imax] ** 2 + er[imin] ** 2)
    print(f"\ndepth of the dip: {dip:.3f}")
    print(f"  on the PLOTTED (scaled) errors : {dip/sg:.1f} sigma")
    print(f"  on the propagated errors alone : {dip/sr:.1f} sigma")


if __name__ == "__main__":
    main()
