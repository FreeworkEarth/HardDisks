#!/usr/bin/env python3
"""Meeting version (2026-09-11) of the c_s(eta) figures.
Changes vs plot_cs_final_style.py:
  * Kolafa-Rottner red, Roman et al. black, our data blue.
  * Error bars = 1 sigma scatter of c_s across the 9 divider masses (per-mass estimates from the
    trajectory manifest), added in quadrature with the fit error. The fit error alone (0.1-0.2 %)
    is smaller than the marker and hides the real mass-to-mass systematic.
  * Points are corrected for the seeding temperature T_i (nu -> nu / sqrt(T_i), c_s ~ sqrt(T));
    the uncorrected values are drawn as hollow grey markers in the zoom figure.
  * Ideal-gas limit (eta -> 0, c_s = sqrt 2) shown as the exact point at the origin.
    The first-order low-density expansion c_s = sqrt(2)(1 + 2 eta) was removed
    2026-09-12: it is a one-term truncation that is already 1 % low by eta = 0.05 and
    was being read as a second reference curve next to Kolafa-Rottner.
usage: python3 plot_cs_meeting.py routeA_refit_cs_vs_eta_20260909.csv routeA_fit_input_manifest_20260909.csv 260911_MEETING [raw|corr|both]
  raw  = as measured, blue (meeting default)   corr = T_i-corrected only   both = corrected blue + uncorrected hollow grey
"""
import os, sys, csv, math
import numpy as np
import pandas as pd
from scipy.optimize import brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CSV, MANIFEST, PREFIX = sys.argv[1], sys.argv[2], sys.argv[3]
MODE = sys.argv[4] if len(sys.argv) > 4 else "both"   # raw | corr | both | pre | v2
# ##CHRIS 2026-09-13: mode v2 = A1 v2 data (drift-first seeding, T_i = 1 exactly, no correction); argv[5] names the estimator
V2_LABEL = sys.argv[5] if len(sys.argv) > 5 else "Román 2002 FFT estimator"
N_SIDE = 50  # disks per compartment, route A

KR_COEF = {0: 1.0, 1: 2.0, 2: 1.12801775, 3: 0.00181895291, 4: -0.0526134737, 5: 0.0504960168,
           6: -0.0325537792, 7: 0.0134578632, 8: 0.00140888182, 9: -0.00834273601, 10: 0.00694127367,
           11: -0.00262254723, 12: 0.000355746352, 22: -5.24672938e-9}
def Z_kr(e):  x = e / (1 - e); return sum(c * x ** p for p, c in KR_COEF.items())
def Z_spt(e): return 1.0 / (1 - e) ** 2
def Z_hen(e, a=0.125): return (1 + a * e * e) / (1 - e) ** 2
def cs(Zf, e, h=1e-6):
    e = np.asarray(e, float); Z = Zf(e); dZ = (Zf(e + h) - Zf(e - h)) / (2 * h)
    return np.sqrt(np.maximum(Z + e * dZ + Z * Z, 0))

# Roman et al. (2002), N=100, r=0.5, H=10, keyed by L0
# ##CHRIS 2026-09-13: Table I of Roman et al., Am. J. Phys. 70, 847 (2002), checked against the
# paper. L0 = 25 (eta = 0.157) was transcribed as 2.10; the table gives 2.01 +- 0.02.
ROMAN = {7.5: (5.99, 0.09), 10.0: (3.78, 0.08), 15.0: (2.61, 0.03), 20.0: (2.20, 0.02),
         25.0: (2.01, 0.02), 30.0: (1.89, 0.02), 35.0: (1.81, 0.02)}
rom = sorted(((100 * math.pi * 0.25) / (2 * L * 10.0), v[0], v[1]) for L, v in ROMAN.items())

# ---- per-mass scatter and T_i correction from the trajectory manifest ----
def K_of(M):  # divider mode: cot K = (M / 2 N m) K
    return brentq(lambda K: 1 / np.tan(K) - (M / (2 * N_SIDE)) * K, 1e-6, np.pi / 2 - 1e-9)
man = pd.read_csv(MANIFEST); man = man[man.eligible == 1]
corr = {}
for eta, g in man.groupby("eta"):
    Leff = g.L0.iloc[0] - 1.0
    x = np.array([K_of(M) for M in g.M]) / (2 * np.pi * Leff)
    raw = (g.nu * x).sum() / (x * x).sum()
    nuc = g.nu / np.sqrt(g.T_i_mean)
    cor = (nuc * x).sum() / (x * x).sum()
    pm = g.assign(nuc=nuc).groupby("M").nuc.mean().reset_index()
    xm = np.array([K_of(M) for M in pm.M]) / (2 * np.pi * Leff)
    csM = pm.nuc.values / xm
    corr[round(eta, 4)] = (cor / raw, csM.std(ddof=1) / csM.mean(), len(csM), g.T_i_mean.mean())

sim = []  # (eta, cs_corr, err_total, cs_raw, fit_err)
for r in csv.DictReader(open(CSV)):
    eta, c, fe = float(r["eta"]), float(r["c_s"]), float(r["c_s_err"])
    if MODE in ("pre", "v2", "final"):
        # ##CHRIS 2026-09-12: the input table is ALREADY T_i-corrected upstream (robust
        # estimator). Applying the manifest correction again would double-count it and
        # inflate every deviation by ~1 %. Take c_s and its mass scatter verbatim.
        f, sc = 1.0, float(r.get("c_s_scatter_mass", 0) or 0) / c
    elif round(eta, 4) in corr:
        f, sc, nm, Ti = corr[round(eta, 4)]
    else:  # new rows (drift-first seeding, T_i = 1): no correction, scatter from the CSV if it has it
        f, sc = 1.0, float(r.get("c_s_scatter_mass", 0) or 0) / c
    cc = c * f
    err_raw = math.hypot(fe, sc * c)
    sim.append((eta, cc, math.hypot(fe * f, sc * cc), c, err_raw))
sim.sort()

REGIONS = [(0.000, 0.020, "ideal-gas\nlimit", "#f4f4f4"), (0.020, 0.200, "dilute gas-like\nfluid", "#dceeff"),
           (0.200, 0.500, "dense isotropic\nfluid", "#dff4e3"), (0.500, 0.700, "stiff dense\nisotropic fluid", "#fff0c8"),
           (0.700, 0.716, "fluid-hexatic\ncoexist.", "#ffd9b8"), (0.716, 0.720, "hexatic", "#f6c7df"),
           (0.720, 0.780, "solid", "#ded8f5")]
def shade(ax, x_min, x_max):
    for lo, hi, label, color in REGIONS:
        lo, hi = max(lo, x_min), min(hi, x_max)
        if hi <= lo: continue
        ax.axvspan(lo, hi, color=color, alpha=0.32, lw=0, zorder=0)
        if hi - lo >= 0.03 * (x_max - x_min) / 0.78:
            ax.text((lo + hi) / 2, 0.975, label, transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=8, color="#333333")
        else:
            y = 0.38 if lo >= 0.716 else 0.70
            ax.text((lo + hi) / 2, y, label, transform=ax.get_xaxis_transform(), ha="center", va="center", rotation=90, fontsize=6.5, color="#333333")

KR_RED = "#d62728"
def eos_curves(ax, x_max, lw=1.4):
    e = np.linspace(1e-4, min(x_max, 0.78), 600)
    ax.plot(e, cs(Z_spt, e), "--", color="C1", lw=lw, label="SPT equation of state")
    ax.plot(e, cs(Z_hen, e), "-.", color="C2", lw=lw, label="Henderson (a = 0.125) equation of state")
    ek = e[e <= 0.69]
    ax.plot(ek, cs(Z_kr, ek), "-", color=KR_RED, lw=lw + 1.0, label="Kolafa–Rottner 2006 (main reference, valid to η ≈ 0.69)")
    ax.axhline(math.sqrt(2), color="0.45", lw=1.0, ls=":", label="ideal-gas limit η → 0:  c_s = √(2 k_BT/m) = √2")

ERR_NOTE = "error bar = 1σ scatter of c_s over the 9 divider masses ⊕ fit error"

# ---------------- figure 1: full range ----------------
fig, ax = plt.subplots(figsize=(11, 6.4))
shade(ax, 0, 0.78)
eos_curves(ax, 0.78)
ax.errorbar([p[0] for p in rom], [p[1] for p in rom], yerr=[p[2] for p in rom], fmt="^--", color="black", capsize=3, alpha=0.9, ms=7, lw=1.0,
            label="Román et al. (2002), N = 100 (literature)")
s = [p for p in sim if p[0] <= 0.74]
if MODE in ("v2", "final"):
    # ##CHRIS 2026-09-14: mode "final" takes the whole legend text verbatim from argv[5].
    ax.errorbar([p[0] for p in s], [p[1] for p in s], yerr=[p[2] for p in s], fmt="o-", color="C0", capsize=3, ms=5, lw=1.2, zorder=4,
                label=V2_LABEL if MODE == "final" else
                      "EDMD simulation, A1 v2: N = 100 (50 per compartment), r = 0.5, H = 10, drift-first seeding (T_i = 1)\n"
                      "9 divider masses × 25 seeds per η · c_s from the " + V2_LABEL + "\n"
                      "health-clean · error bar = 1σ scatter of per-mass c_s")
elif MODE == "raw":
    ax.errorbar([p[0] for p in s], [p[3] for p in s], yerr=[p[4] for p in s], fmt="o-", color="C0", capsize=3, ms=5, lw=1.2, zorder=4,
                label="EDMD simulation, route A: N = 100 (50 per compartment), r = 0.5, H = 10,\n"
                      "L₀ = 3.93 σ/η · 9 divider masses × 25 repeats per η · as measured (runs seeded at T_i ≈ 0.98)\n"
                      "health-clean · " + ERR_NOTE)
else:
    ax.errorbar([p[0] for p in s], [p[1] for p in s], yerr=[p[2] for p in s], fmt="o-", color="C0", capsize=3, ms=5, lw=1.2, zorder=4,
                label="EDMD simulation, route A: N = 100 (50 per compartment), r = 0.5, H = 10,\n"
                      "L₀ = 3.93 σ/η · 9 divider masses × 25 repeats per η · T_i-corrected (ν/√T_i)\n"
                      "health-clean · " + ERR_NOTE)
ax.set_xlim(0, 0.78); ax.set_ylim(0, 22)
ax.set_xlabel("Packing fraction  η", fontsize=12); ax.set_ylabel("Speed of sound  c_s  [√(k_BT/m)]", fontsize=12)
ax.set_title("Speed of sound c_s vs packing fraction η — EDMD hard disks vs equations of state", fontsize=13)
ax.grid(True, ls=":", alpha=0.6); ax.legend(loc="upper left", bbox_to_anchor=(0.01, 0.84), fontsize=8.5, framealpha=0.95)
ax.text(0.655, 2.6, "η ≥ 0.65: 6 σ compartment,\nstructure changes during\nmeasurement — not a\nfluid-branch value", fontsize=7.5, color="#333333", ha="right", va="bottom")
fig.text(0.99, 0.005, "adiabatic mapping c_s² = (k_BT/m)[Z + ηZ′ + Z²] · divider mode ν = c_s K/(2π L_eff), L_eff = L₀ − 2r, cot K = (M/2Nm)K · kT = m = σ = 1 · data: " + os.path.basename(CSV) + " + " + os.path.basename(MANIFEST),
         ha="right", va="bottom", fontsize=7, color="0.4")
fig.tight_layout(rect=(0, 0.02, 1, 1))
fig.savefig(PREFIX + "_cs_vs_eta.png", dpi=200); fig.savefig(PREFIX + "_cs_vs_eta.pdf")

# ---------------- figure 2: ideal-gas zoom, two panels ----------------
XMAX = 0.105
fig, (ax, axd) = plt.subplots(2, 1, figsize=(9.5, 7.8), sharex=True, gridspec_kw={"height_ratios": [2.2, 1]})
for a in (ax, axd): shade(a, 0, XMAX)
eos_curves(ax, XMAX, lw=1.4)
e = np.linspace(0, XMAX, 200)
FO = "#7b3fa0"
ax.plot([0], [math.sqrt(2)], "*", ms=14, color="black", zorder=5, clip_on=False, label="exact ideal-gas point (η = 0, c_s = √2 = 1.4142)")
z = [p for p in sim if p[0] <= XMAX]
if MODE in ("v2", "final"):
    ax.errorbar([p[0] for p in z], [p[1] for p in z], yerr=[p[2] for p in z], fmt="o", color="C0", capsize=3, ms=6, lw=1.3, zorder=4,
                label=V2_LABEL if MODE == "final" else
                      "EDMD A1 v2, drift-first seeding (T_i = 1); c_s from the " + V2_LABEL + "\nerror bar = 1σ scatter of per-mass c_s")
elif MODE == "raw":
    ax.errorbar([p[0] for p in z], [p[3] for p in z], yerr=[p[4] for p in z], fmt="o", color="C0", capsize=3, ms=6, lw=1.3, zorder=4,
                label="EDMD route A, as measured (N = 100, L₀ = 200 … 35 σ; runs seeded at T_i ≈ 0.98)\n" + ERR_NOTE)
else:
    if MODE == "both":
        ax.errorbar([p[0] for p in z], [p[3] for p in z], yerr=[p[4] for p in z], fmt="o", mfc="white", mec="0.5", ecolor="0.5", capsize=2, ms=6, lw=1.0, zorder=3,
                    label="EDMD route A, uncorrected (seeding gave T_i ≈ 0.98)")
    ax.errorbar([p[0] for p in z], [p[1] for p in z], yerr=[p[2] for p in z], fmt="o", color="C0", capsize=3, ms=6, lw=1.3, zorder=4,
                label="EDMD route A, T_i-corrected (N = 100, L₀ = 200 … 35 σ)\n" + ERR_NOTE)
ax.set_ylim(1.36, 1.86); ax.set_ylabel("Speed of sound  c_s  [√(k_BT/m)]", fontsize=11)
ax.set_title("Ideal-gas end of the sound-speed validation (η ≤ 0.1)", fontsize=12.5)
ax.grid(True, ls=":", alpha=0.6); ax.legend(loc="upper left", bbox_to_anchor=(0.01, 0.92), fontsize=7.8, framealpha=0.95)
# deviation panel
ez = np.array([p[0] for p in z]); kz = cs(Z_kr, ez)
axd.axhline(0, color=KR_RED, lw=1.6)
axd.axhspan(-0.5, 0.5, color="0.85", alpha=0.6, zorder=0)
axd.plot(e[1:], 100 * (cs(Z_spt, e[1:]) - cs(Z_kr, e[1:])) / cs(Z_kr, e[1:]), "--", color="C1", lw=1.3)
axd.plot(e[1:], 100 * (cs(Z_hen, e[1:]) - cs(Z_kr, e[1:])) / cs(Z_kr, e[1:]), "-.", color="C2", lw=1.3)
if MODE == "raw":
    axd.errorbar(ez, 100 * (np.array([p[3] for p in z]) - kz) / kz, yerr=100 * np.array([p[4] for p in z]) / kz, fmt="o", color="C0", capsize=3, ms=6, lw=1.3, zorder=4)
else:
    if MODE == "both":
        axd.errorbar(ez, 100 * (np.array([p[3] for p in z]) - kz) / kz, yerr=100 * np.array([p[4] for p in z]) / kz, fmt="o", mfc="white", mec="0.5", ecolor="0.5", capsize=2, ms=6, lw=1.0, zorder=3)
    axd.errorbar(ez, 100 * (np.array([p[1] for p in z]) - kz) / kz, yerr=100 * np.array([p[2] for p in z]) / kz, fmt="o", color="C0", capsize=3, ms=6, lw=1.3, zorder=4)
axd.set_ylim(-1.6, 2.0); axd.set_ylabel("deviation from Kolafa–Rottner  [%]", fontsize=10); axd.set_xlabel("Packing fraction  η", fontsize=12)
axd.text(0.002, 0.55, "±0.5 % band", fontsize=8, color="0.35", va="bottom")
if MODE == "both": axd.text(0.103, 1.85, "hollow grey = uncorrected, blue = T_i-corrected", fontsize=7.5, color="0.3", ha="right", va="top")
if MODE == "raw":  axd.text(0.103, -1.45, "as measured; runs were seeded at T_i ≈ 0.98 and c_s ∝ √T, i.e. ≈ 1 % low vs a T = 1 gas", fontsize=7.5, color="0.3", ha="right", va="bottom")
axd.grid(True, ls=":", alpha=0.6)
ax.set_xlim(0, XMAX)
fig.text(0.99, 0.005, "lower panel: (c_s − c_s,KR)/c_s,KR for the EOS (lines) and the simulation (points) · KR line in red" + ("" if MODE in ("raw", "v2", "final") else " · T_i correction: ν → ν/√T_i because c_s ∝ √T for hard disks"),
         ha="right", va="bottom", fontsize=7, color="0.4")
fig.tight_layout(rect=(0, 0.02, 1, 1))
fig.savefig(PREFIX + "_cs_idealgas_zoom.png", dpi=200); fig.savefig(PREFIX + "_cs_idealgas_zoom.pdf")
for p in z: print("eta %.4f  raw %.4f ± %.4f   corr %.4f ± %.4f   dev_KR corr %+.2f %%" % (p[0], p[3], p[4], p[1], p[2], 100 * (p[1] / cs(Z_kr, p[0]) - 1)))
print("wrote", PREFIX + "_cs_vs_eta and _cs_idealgas_zoom")
