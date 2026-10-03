#!/usr/bin/env python3
"""##CHRIS 2026-10-01: the four Paper 2 figures the draft lists as missing. Analysis only, no runs.

TWO OF THE FOUR ARE NOT WHAT THEY WERE CALLED, and are drawn as what the data actually is:

  * "zeta(tau)" DOES NOT EXIST. zeta is measured on a PINNED piston during an equilibrium hold, and
    the independent variable is the Green-Kubo integration cutoff t_cut, not a push duration tau.
    There is no zeta-versus-tau or zeta-versus-u dataset anywhere. The figure is therefore
    zeta(t_cut) at the two ends of the Level 1 path, which is the measurement that was made, and
    its point is the SOUND-CROSSING zero: the integral dies at L/c_s, ~21.5 and ~18.0 sigma.

  * "the Level 1 path integral" is not force-versus-position. What is integrated is Z against
    ln eta, since int Z dln(eta) = ln(T_f/T_i). Five held-wall points, each with its own error,
    against the Kolafa-Rottner curve.

Numbers for those two come from the committed reports (260917 section 2/6c, 260919 section 1),
which are the only place they are written down; each is tagged with its source in the code below.
The Level 4 v3 bars are quoted from 260926 rather than re-derived: the exact averaging window used
for the published table is not recoverable from any file, and re-deriving gives slightly different
values, so mixing the two would be worse than quoting one.

  1  zeta(t_cut), both path ends      261001_p2_zeta_tcut
  2  Level 1 path integral            261001_p2_level1_path
  3  Level 4 equilibrium ACF + fit    261001_p2_level4_acf        <- the figure that ties the papers
  4  Level 4 v3 bars                  261001_p2_level4_bars

Colours: our data BLUE #2a78d6, Kolafa-Rottner RED #e34948, Roman/theory BLACK.
"""
import glob, math, os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, brentq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, ".."))
import plot_speed_of_sound_edmd as sos

BLUE, RED, GREY, ORANGE, GREEN = "#2a78d6", "#e34948", "#52514e", "#eb6834", "#1baf7a"
REPO = os.path.dirname(os.path.dirname(HERE))
OUTDIR = os.path.join(REPO, "0000_PLAN_OVERALL", "paper2_energytransfer", "experiments", "final")


def save(fig, name):
    os.makedirs(OUTDIR, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUTDIR, f"{name}.{ext}"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {name}.png/.pdf")


# ------------------------------------------------------------------- 1. zeta(t_cut)
# SOURCE: 260917_paper2_level2_REPORT.md section 2 (start of path) and section 6c (both ends).
def fig_zeta():
    tcut = np.array([0.25, 0.50, 1.00, 2.00, 5.00, 10.00, 20.00, 40.00, 60.00])
    z_start = np.array([2.546, 2.486, 2.387, 2.170, 1.777, 1.219, 0.149, -0.053, 0.148])
    # bin-width convergence at five bin sizes, section 6b -> spread is the systematic on each point
    conv = {1.0: [2.377, 1.790, 1.235, 0.149, -0.051], 0.5: [2.378, 1.793, 1.233, 0.142, -0.049],
            0.25: [2.380, 1.780, 1.222, 0.143, -0.049], 0.125: [2.374, 1.785, 1.230, 0.142, -0.063],
            0.0625: [2.368, 1.781, 1.224, 0.136, -0.063]}
    conv_t = np.array([1.0, 5.0, 10.0, 20.0, 40.0])
    conv_sd = np.array([np.std([conv[b][i] for b in conv], ddof=1) for i in range(5)])
    err = np.interp(tcut, conv_t, conv_sd)

    fig, ax = plt.subplots(figsize=(8.4, 5.4))
    ax.errorbar(tcut, z_start, yerr=err, fmt="o-", color=BLUE, capsize=3, ms=6, lw=1.8,
                label=r"start of path, $L = 35.32\,\sigma$, $\eta = 0.1112$")
    ax.axhline(0, color=GREY, lw=1)
    ax.axhline(2.613, color=BLUE, ls=":", lw=1.4, label=r"zero-lag self term $2.613$")
    ax.axhline(2.720, color=RED, ls="--", lw=1.8, label=r"Enskog $\zeta_E = 2.720$ (start)")
    for xv, lab, col in ((19.8, r"predicted $L/c_s = 19.8$", RED), (21.5, r"measured zero $21.5$", BLUE)):
        ax.axvline(xv, color=col, ls="-." if col == RED else "-", lw=1.6, alpha=0.8)
        ax.annotate(lab, xy=(xv, 2.0), rotation=90, fontsize=8, color=col,
                    ha="right", va="top")
    ax.set_xscale("log")
    ax.set_xlabel(r"Green–Kubo cutoff $t_{\rm cut}$  [$\sigma$-time]")
    ax.set_ylabel(r"$\zeta(t_{\rm cut}) = \int_0^{t_{\rm cut}}\langle \delta F(0)\,\delta F(t)\rangle\,{\rm d}t$")
    ax.set_title(r"The friction integral dies at the sound-crossing time, not at a plateau")
    ax.grid(alpha=0.3); ax.legend(frameon=False, fontsize=9, loc="lower left")
    save(fig, "261001_p2_zeta_tcut")


# ------------------------------------------------------------- 2. Level 1 path integral
# SOURCE: 260919_paper2_ramp_linewidth_fast_REPORT.md section 1 (five held-wall path points).
def fig_level1_path():
    L = np.array([38.750000, 37.770833, 36.791667, 35.812500, 34.812500])
    eta_true = np.array([0.101342, 0.103969, 0.106736, 0.109654, 0.112804])
    Zw = np.array([1.2768, 1.2935, 1.3052, 1.3064, 1.3241])
    dZw = np.array([0.0027, 0.0026, 0.0025, 0.0027, 0.0027])
    Zkr = np.array([1.2399, 1.2472, 1.2551, 1.2634, 1.2725])
    lne = np.log(eta_true)
    # The PUBLISHED path integral (260919 report section 1) is 0.139215 -> T_f/T_i = 1.149371,
    # giving W_qs = 7.4685 +- 0.0131. A plain trapezoid over these five points gives 0.139487,
    # 0.2 % higher: the quadrature choice matters at the third digit. The figure quotes the
    # published value, and the trapezoid is printed so the difference is on the record rather
    # than hidden inside a title.
    PUB_INTEGRAL, PUB_RATIO = 0.139215, 1.149371
    integral = PUB_INTEGRAL
    trap = float(np.trapezoid(Zw, lne)) if hasattr(np, "trapezoid") else float(np.trapz(Zw, lne))

    fig, (ax, axr) = plt.subplots(1, 2, figsize=(12.4, 5.2), gridspec_kw={"width_ratios": [1.5, 1]})
    g = np.linspace(eta_true.min() * 0.985, eta_true.max() * 1.015, 200)
    ax.plot(np.log(g), sos.Z_kolafa_rottner_2006(g), "-", color=RED, lw=2,
            label="Kolafa–Rottner 2006 (bulk EOS)")
    ax.errorbar(lne, Zw, yerr=dZw, fmt="o", color=BLUE, capsize=3, ms=7, zorder=3,
                label=r"measured $Z_{\rm wall}$, held wall, 100 seeds/point")
    ax.fill_between(lne, Zw - dZw, Zw + dZw, color=BLUE, alpha=0.15)
    ax.fill_between(lne, Zkr, Zw, color=BLUE, alpha=0.07)
    ax.set_xlabel(r"$\ln \eta$")
    ax.set_ylabel(r"$Z = PV/N k_B T$")
    ax.set_title(r"Level 1: $\int Z\,{\rm d}\ln\eta = %.6f \;\Rightarrow\; T_f/T_i = %.6f$"
                 % (PUB_INTEGRAL, PUB_RATIO))
    ax.grid(alpha=0.3); ax.legend(frameon=False, fontsize=9, loc="upper left")
    ax.text(0.03, 0.05, r"$W_{\rm qs} = 7.4685 \pm 0.0131\,k_BT$", transform=ax.transAxes,
            fontsize=11, bbox=dict(fc="white", ec=GREY, alpha=0.9))

    ratio = 100.0 * (Zw / Zkr - 1.0)
    axr.errorbar(L, ratio, yerr=100.0 * dZw / Zkr, fmt="o-", color=BLUE, capsize=3, ms=6, lw=1.6)
    axr.axhline(0, color=RED, ls="--", lw=1.8, label="bulk KR")
    axr.axhline(ratio.mean(), color=GREY, ls=":", lw=1.6,
                label=f"mean {ratio.mean():+.2f} %")
    axr.set_xlabel(r"box length $L$  [$\sigma$]")
    axr.set_ylabel(r"$Z_{\rm wall}/Z_{\rm KR} - 1$  [%]")
    axr.set_title("The box is over-pressured at every path point")
    axr.grid(alpha=0.3); axr.legend(frameon=False, fontsize=9)
    save(fig, "261001_p2_level1_path")
    print(f"     published integral {PUB_INTEGRAL:.6f} -> T_f/T_i = {PUB_RATIO:.6f}")
    print(f"     plain trapezoid over the 5 plotted points {trap:.6f} ({100*(trap/PUB_INTEGRAL-1):+.2f} % vs published)")
    print(f"     mean Z excess over KR {ratio.mean():+.2f} %")


# -------------------------------------------- 3. Level 4 equilibrium ACF -- the tie between papers
def fig_level4_acf():
    P = os.path.join(REPO, "hspist3", "experiments_energy_transfer",
                     "level4_equilibrium_20260929", "Md10")
    fs = sorted(glob.glob(os.path.join(P, "red_*.csv")))
    if not fs:
        print(f"  MISSING {P}"); return
    NS, dt = 50, 5.0
    D, X = [], []
    for f in fs:
        e = pd.read_csv(f)
        D.append((e["KE_gas_left"].to_numpy(float) - e["KE_gas_right"].to_numpy(float)) / NS)
        X.append(e["W0_x_sigma"].to_numpy(float))
    n = min(len(a) for a in D); lo = int(2000 / dt)
    NL = int(1500 / dt) + 1; lag = np.arange(NL) * dt

    eta = 50 * math.pi * 0.25 / (39.25 * 10)
    Z = float(sos.Z_kolafa_rottner_2006(np.array([eta]))[0])
    dZ = float(sos.dZ_kolafa_rottner_2006(np.array([eta]))[0])
    cs = math.sqrt(Z + eta * dZ + Z * Z)
    K = brentq(lambda k: math.cos(k) / math.sin(k) - 0.1 * k, 1e-9, math.pi - 1e-9)
    LEFF = 37.75
    nu = cs * K / (2 * math.pi * LEFF); OM = 2 * math.pi * nu; PER = 1 / nu
    nu_id = math.sqrt(2.0) * K / (2 * math.pi * LEFF)

    def acf1(x):
        x = x - x.mean(); m = len(x)
        f = np.fft.rfft(x, 2 * m)
        c = np.fft.irfft(f * np.conj(f))[:NL].real
        return c / c[0]

    def model(t, A, tT, B, tr):
        return A * np.exp(-t / tT) + B * np.exp(-t / tr) * np.cos(OM * t)

    fig, (ax, axz) = plt.subplots(1, 2, figsize=(12.8, 5.2), gridspec_kw={"width_ratios": [1.6, 1]})
    res = {}
    for nm, S, colour, mark in ((r"$T_1-T_2$", D, BLUE, "o"), ("divider $x$", X, ORANGE, "s")):
        cs_all = np.array([acf1(a[lo:n]) for a in S])
        c = cs_all.mean(0)
        err = cs_all.std(0, ddof=1) / math.sqrt(len(S))
        p, _ = curve_fit(model, lag, c, p0=[0.5, 500., 0.5, 200.],
                         bounds=([0, 20, 0, 20], [1.5, 1e4, 1.5, 5e3]), maxfev=40000)
        res[nm] = p
        ax.errorbar(lag[::3], c[::3], yerr=err[::3], fmt=mark, color=colour, ms=3.4, lw=0,
                    elinewidth=0.8, alpha=0.65, capsize=0)
        ax.plot(lag, model(lag, *p), "-", color=colour, lw=2.0,
                label=f"{nm}: $\\tau_T = {p[1]:.0f}$, $\\tau_r = {p[3]:.0f}$")
        axz.errorbar(lag[:41], c[:41], yerr=err[:41], fmt=mark, color=colour, ms=4.5, lw=0,
                     elinewidth=0.9, alpha=0.7)
        axz.plot(lag[:41], model(lag[:41], *p), "-", color=colour, lw=2.0)
    for a in (ax, axz):
        a.axhline(0, color=GREY, lw=1)
        a.axvline(PER, color=RED, ls="--", lw=2)
    ax.annotate(f"predicted period {PER:.1f}\n(KR, $\\cot K=\\alpha K$, $\\alpha=0.1$)",
                xy=(PER, 0.72), xytext=(PER * 1.9, 0.80), fontsize=9, color=RED,
                arrowprops=dict(arrowstyle="->", color=RED))
    axz.axvline(1 / nu_id, color="k", ls=":", lw=1.8)
    # The discriminating feature of the whole figure: the second ACF peak sits ON the KR line and
    # clearly OFF the ideal-gas one. Both labels must be inside the axes or the point is lost.
    axz.annotate(f"KR {PER:.1f}", xy=(PER, 0.93), fontsize=9.5, color=RED, ha="right",
                 va="center", fontweight="bold")
    axz.annotate(f"ideal gas {1/nu_id:.1f}", xy=(1 / nu_id, 0.93), fontsize=9.5, color="k",
                 ha="left", va="center")
    ax.set_xlabel(r"lag  [$\sigma$-time]"); ax.set_ylabel("autocorrelation")
    ax.set_title(r"Level 4 equilibrium ACF: a slow isobaric mode $+$ Paper 1's divider resonance")
    ax.grid(alpha=0.3); ax.legend(frameon=False, fontsize=9)
    axz.set_xlim(0, 200); axz.set_ylim(-0.12, 1.06); axz.set_xlabel(r"lag  [$\sigma$-time]")
    axz.set_title("First two periods")
    axz.grid(alpha=0.3)
    save(fig, "261001_p2_level4_acf")
    print(f"     predicted period {PER:.2f} (KR) vs ideal {1/nu_id:.2f};"
          f" fits tau_T={res[r'$T_1-T_2$'][1]:.0f}/{res['divider $x$'][1]:.0f}")


# -------------------------------------------------------------------- 4. Level 4 v3 bars
# SOURCE: 260926_paper2_level4_v3.md section 2. Quoted, not re-derived -- see the module docstring.
def fig_level4_bars():
    M = np.array([10, 50, 200]); seeds = [40, 40, 8]
    far = np.array([0.488, 0.438, 0.276]); dfar = np.array([0.126, 0.127, 0.153])
    setl = np.array([1.867, 2.245, 2.012]); dsetl = np.array([0.281, 0.214, 0.377])
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.6, 4.8))
    xp = np.arange(3)
    a1.bar(xp, far, yerr=dfar, color=BLUE, alpha=0.85, capsize=6, width=0.55, ecolor=GREY)
    a1.axhline(0.5, color="k", ls="--", lw=2, label=r"prediction $1/2$")
    a1.set_xticks(xp); a1.set_xticklabels([f"$M_d = {m}$\n{s} seeds" for m, s in zip(M, seeds)])
    a1.set_ylabel("far gas / total work"); a1.set_ylim(0, 0.75)
    a1.set_title("(i) the work splits in half")
    for i, (v, e) in enumerate(zip(far, dfar)):
        a1.annotate(f"{abs(v-0.5)/e:.1f}$\\sigma$", (i, v + e + 0.02), ha="center", fontsize=9)
    a1.legend(frameon=False); a1.grid(alpha=0.3, axis="y")

    a2.bar(xp, setl, yerr=dsetl, color=BLUE, alpha=0.85, capsize=6, width=0.55, ecolor=GREY)
    a2.axhline(1.965, color="k", ls="--", lw=2, label=r"prediction $\Delta x/2 = 1.965$")
    a2.set_xticks(xp); a2.set_xticklabels([f"$M_d = {m}$\n{s} seeds" for m, s in zip(M, seeds)])
    a2.set_ylabel(r"settled divider displacement  [$\sigma$]"); a2.set_ylim(0, 3.0)
    a2.set_title(r"(iii) the divider takes $\Delta x/2$")
    for i, (v, e) in enumerate(zip(setl, dsetl)):
        a2.annotate(f"{abs(v-1.965)/e:.1f}$\\sigma$", (i, v + e + 0.06), ha="center", fontsize=9)
    a2.legend(frameon=False); a2.grid(alpha=0.3, axis="y")
    save(fig, "261001_p2_level4_bars")


if __name__ == "__main__":
    print("Paper 2 figures ->", OUTDIR)
    fig_zeta()
    fig_level1_path()
    fig_level4_acf()
    fig_level4_bars()
