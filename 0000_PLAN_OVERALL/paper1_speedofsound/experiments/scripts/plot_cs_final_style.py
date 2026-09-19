#!/usr/bin/env python3
"""c_s(eta) figures in the style of the 2026-08-20 'FINAL speed_of_sound_on_packing_fracture'
plot: regime shading, SPT / Henderson / Kolafa-Rottner EOS through the adiabatic mapping
c_s^2 = (kT/m)[Z + eta Z' + Z^2], the Roman et al. (2002) literature points, and the
route-A refit (2026-09-09, through the origin, eligible trajectories).
  argv[1] = routeA_refit_cs_vs_eta_20260909.csv   argv[2] = output prefix
Produces <prefix>_cs_vs_eta.{png,pdf} and <prefix>_cs_idealgas_zoom.{png,pdf}.
"""
import sys, csv, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CSV, PREFIX = sys.argv[1], sys.argv[2]

KR_COEF = {0: 1.0, 1: 2.0, 2: 1.12801775, 3: 0.00181895291, 4: -0.0526134737, 5: 0.0504960168,
           6: -0.0325537792, 7: 0.0134578632, 8: 0.00140888182, 9: -0.00834273601, 10: 0.00694127367,
           11: -0.00262254723, 12: 0.000355746352, 22: -5.24672938e-9}
def Z_kr(e):  x = e / (1 - e); return sum(c * x ** p for p, c in KR_COEF.items())
def Z_spt(e): return 1.0 / (1 - e) ** 2
def Z_hen(e, a=0.125): return (1 + a * e * e) / (1 - e) ** 2
def cs(Zf, e, h=1e-6):
    e = np.asarray(e, float); Z = Zf(e); dZ = (Zf(e + h) - Zf(e - h)) / (2 * h)
    return np.sqrt(np.maximum(Z + e * dZ + Z * Z, 0))

# Roman et al. (2002) reference values, keyed by L0 (N=100, r=0.5, H=10) as in plot_speed_of_sound_edmd.py
# ##CHRIS 2026-09-13: Roman et al. 2002 Table I gives 2.01 +- 0.02 (was mistranscribed as 2.10)
ROMAN = {7.5: (5.99, 0.09), 10.0: (3.78, 0.08), 15.0: (2.61, 0.03), 20.0: (2.20, 0.02),
         25.0: (2.01, 0.02), 30.0: (1.89, 0.02), 35.0: (1.81, 0.02)}
rom = sorted(((100 * math.pi * 0.25) / (2 * L * 10.0), v[0], v[1]) for L, v in ROMAN.items())

sim = []
for r in csv.DictReader(open(CSV)):
    sim.append((float(r["eta"]), float(r["c_s"]), float(r["c_s_err"])))
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

def eos_curves(ax, x_max, lw=1.6):
    e = np.linspace(1e-4, min(x_max, 0.78), 600)
    ax.plot(e, cs(Z_spt, e), "--", color="C1", lw=lw, label="SPT")
    ax.plot(e, cs(Z_hen, e), "-.", color="C2", lw=lw, label="Henderson (a = 0.125)")
    ek = e[e <= 0.69]
    ax.plot(ek, cs(Z_kr, ek), "-", color="black", lw=lw + 0.8, label="Kolafa–Rottner 2006 (main reference, valid to η ≈ 0.69)")
    ax.axhline(math.sqrt(2), color="0.45", lw=1.0, ls=":", label="ideal gas  c_s = √(2 k_BT/m)")

# ---------------- figure 1: full range ----------------
fig, ax = plt.subplots(figsize=(11, 6.4))
shade(ax, 0, 0.78)
eos_curves(ax, 0.78)
ax.errorbar([p[0] for p in rom], [p[1] for p in rom], yerr=[p[2] for p in rom], fmt="^--", color="C3", capsize=3, alpha=0.85, ms=7, label="Román et al. (2002), N = 100")
s = [p for p in sim if p[0] <= 0.74]
ax.errorbar([p[0] for p in s], [p[1] for p in s], yerr=[p[2] for p in s], fmt="o-", color="C0", capsize=3, ms=5, lw=1.2,
            label="EDMD simulation fit, route A: N = 100 (50 per compartment), r = 0.5, H = 10, L₀ = 3.93 σ/η\n(refit 2026-09-09: through the origin, 9 divider masses × 25 repeats per η, health-clean)")
ax.set_xlim(0, 0.78); ax.set_ylim(0, 22)
ax.set_xlabel("Packing fraction  η", fontsize=12); ax.set_ylabel("Speed of sound  c_s  [√(k_BT/m)]", fontsize=12)
ax.set_title("Speed of sound c_s vs packing fraction η — EDMD hard disks vs equations of state", fontsize=13)
ax.grid(True, ls=":", alpha=0.6); ax.legend(loc="upper left", bbox_to_anchor=(0.01, 0.87), fontsize=9, framealpha=0.95)
ax.text(0.655, 2.6, "η ≥ 0.65: 6 σ compartment,\nstructure changes during\nmeasurement — not a\nfluid-branch value", fontsize=7.5, color="#333333", ha="right", va="bottom")
fig.text(0.99, 0.005, "adiabatic mapping c_s² = (k_BT/m)[Z + ηZ′ + Z²] · ν = c_s K/(2π L_eff), L_eff = L₀ − 2r, cot K = (M/2Nm)K · kT = m = σ = 1 · data: 260909_plots/routeA_refit_cs_vs_eta_20260909.csv",
         ha="right", va="bottom", fontsize=7, color="0.4")
fig.tight_layout(rect=(0, 0.02, 1, 1))
fig.savefig(PREFIX + "_cs_vs_eta.png", dpi=200); fig.savefig(PREFIX + "_cs_vs_eta.pdf")

# ---------------- figure 2: ideal-gas zoom, two panels ----------------
XMAX = 0.105
fig, (ax, axd) = plt.subplots(2, 1, figsize=(9.5, 7.6), sharex=True, gridspec_kw={"height_ratios": [2.2, 1]})
for a in (ax, axd): shade(a, 0, XMAX)
eos_curves(ax, XMAX, lw=1.5)
e = np.linspace(0, XMAX, 200)
ax.plot(e, math.sqrt(2) * (1 + 2 * e), ls=(0, (1, 1.5)), color="0.3", lw=1.2, label="c_s = √2 (1 + 2η)  — exact first order")
ax.plot([0], [math.sqrt(2)], "*", ms=13, color="black", zorder=5, clip_on=False)
z = [p for p in sim if p[0] <= XMAX]
ax.errorbar([p[0] for p in z], [p[1] for p in z], yerr=[p[2] for p in z], fmt="o", color="C0", capsize=3, ms=6, lw=1.2, zorder=4,
            label="EDMD simulation fit, route A (N = 100, L₀ = 200 … 35 σ)")
ax.set_ylim(1.36, 1.86); ax.set_ylabel("Speed of sound  c_s  [√(k_BT/m)]", fontsize=11)
ax.set_title("Ideal-gas end of the sound-speed validation (η ≤ 0.1)", fontsize=12.5)
ax.grid(True, ls=":", alpha=0.6); ax.legend(loc="upper left", bbox_to_anchor=(0.01, 0.90), fontsize=8.5, framealpha=0.95)
# deviation panel
ez = np.array([p[0] for p in z]); kz = cs(Z_kr, ez)
axd.axhline(0, color="black", lw=1.2)
axd.axhspan(-0.5, 0.5, color="0.85", alpha=0.6, zorder=0)
axd.plot(e[1:], 100 * (cs(Z_spt, e[1:]) - cs(Z_kr, e[1:])) / cs(Z_kr, e[1:]), "--", color="C1", lw=1.3)
axd.plot(e[1:], 100 * (cs(Z_hen, e[1:]) - cs(Z_kr, e[1:])) / cs(Z_kr, e[1:]), "-.", color="C2", lw=1.3)
axd.plot(e[1:], 100 * (math.sqrt(2) * (1 + 2 * e[1:]) - cs(Z_kr, e[1:])) / cs(Z_kr, e[1:]), ls=(0, (1, 1.5)), color="0.3", lw=1.2)
axd.errorbar(ez, 100 * (np.array([p[1] for p in z]) - kz) / kz, yerr=100 * np.array([p[2] for p in z]) / kz, fmt="o", color="C0", capsize=3, ms=6, lw=1.2, zorder=4)
axd.set_ylim(-1.6, 1.6); axd.set_ylabel("deviation from KR  [%]", fontsize=10); axd.set_xlabel("Packing fraction  η", fontsize=12)
axd.text(0.002, 0.55, "±0.5 % band", fontsize=8, color="0.35", va="bottom")
axd.grid(True, ls=":", alpha=0.6)
ax.set_xlim(0, XMAX)
fig.text(0.99, 0.005, "same curves and data as the full-range figure · lower panel: (c_s − c_s,KR)/c_s,KR for the EOS (lines) and the simulation (points, error bars = fit error)",
         ha="right", va="bottom", fontsize=7, color="0.4")
fig.tight_layout(rect=(0, 0.02, 1, 1))
fig.savefig(PREFIX + "_cs_idealgas_zoom.png", dpi=200); fig.savefig(PREFIX + "_cs_idealgas_zoom.pdf")
print("wrote", PREFIX + "_cs_vs_eta and _cs_idealgas_zoom; sim points:", len(sim), "roman:", len(rom))
