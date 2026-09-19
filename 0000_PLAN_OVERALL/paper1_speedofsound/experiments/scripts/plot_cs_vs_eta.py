#!/usr/bin/env python3
"""Speed of sound vs packing fraction: EOS prediction, ideal-gas limit, and the
simulation points available to Cowork on 2026-09-09.

Theory: Roman relation c_s^2 = (kT/m) [Z + eta Z' + Z^2] with Z from
Kolafa-Rottner 2006 (valid to eta ~ 0.69). Ideal gas: c_s = sqrt(2 kT/m) in 2D.
Low-density expansion (exact through first order): c_s = sqrt(2) (1 + 2 eta + ...).

Simulation points: the six route-A / route-B densities recorded in
260907_speedsound_recap_and_critical_review_COWORK.md (from CC's fits of the
25-repeat campaigns). The full 33-density route-A sweep lives in
campaign_r25_psi6_20260823/analysis/final_plots/speed_of_sound_summary.csv, which
is one folder too deep for the bridge; pass it as argv[2] (columns: eta, c_s,
c_s_err) to replace the six points.
"""
import sys, csv, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = sys.argv[1]
FULL = sys.argv[2] if len(sys.argv) > 2 else None

BLUE, ORANGE, INK, INK2, MUTED, GRID, SURF = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#8a8984", "#e6e5e1", "#fcfcfb"

KR_COEF = {0: 1.0, 1: 2.0, 2: 1.12801775, 3: 0.00181895291, 4: -0.0526134737,
           5: 0.0504960168, 6: -0.0325537792, 7: 0.0134578632, 8: 0.00140888182,
           9: -0.00834273601, 10: 0.00694127367, 11: -0.00262254723,
           12: 0.000355746352, 22: -5.24672938e-9}
def KR(eta):
    x = eta / (1 - eta)
    return sum(c * x ** p for p, c in KR_COEF.items())
def cs_from_Z(Zfun, eta, h=1e-5):
    Z = Zfun(eta); dZ = (Zfun(eta + h) - Zfun(eta - h)) / (2 * h)
    return np.sqrt(Z + eta * dZ + Z ** 2)

# ---- simulation points on record (route A: r=0.5, L0 per eta, H=10, 50 disks/compartment; route B: L0=20, r swept)
routeA = [(0.026, 1.504), (0.157, 2.001), (0.393, 3.775), (0.550, 6.782), (0.630, 9.581), (0.700, 19.30)]
routeB = [(0.026, 1.504), (0.157, 1.986), (0.393, 3.849), (0.550, 7.075), (0.630, 11.657), (0.700, 19.86)]
full = None
if FULL:
    full = []
    for r in csv.DictReader(open(FULL)):
        try:
            full.append((float(r["eta"]), float(r["c_s"]), float(r.get("c_s_err", "nan"))))
        except (KeyError, ValueError):
            pass

fig, ax = plt.subplots(figsize=(9.2, 5.8))
fig.patch.set_facecolor(SURF); ax.set_facecolor(SURF)
ax.grid(True, color=GRID, lw=0.8); ax.set_axisbelow(True)
for s in ("top", "right"): ax.spines[s].set_visible(False)
for s in ("left", "bottom"): ax.spines[s].set_color(MUTED)
ax.tick_params(colors=INK2, labelsize=9)

eg1 = np.linspace(1e-4, 0.69, 500); eg2 = np.linspace(0.69, 0.76, 60)
ax.plot(eg1, cs_from_Z(KR, eg1), "-", color=INK2, lw=1.8, label="EOS prediction (Kolafa–Rottner via Román, valid to η ≈ 0.69)")
ax.axhline(math.sqrt(2), color=MUTED, lw=1.3, ls=(0, (2, 2)))
ax.text(0.08, math.sqrt(2) * 0.94, "ideal gas  c_s = √(2 k_BT/m) = 1.414", color=MUTED, fontsize=8.5, va="top")
ax.axvspan(0.700, 0.716, color=GRID, alpha=0.8, zorder=0)
ax.text(0.695, 2.05, "liquid–hexatic\ncoexistence\n0.700–0.716 →", color=MUTED, fontsize=7.5, ha="right", va="bottom")
ax.text(0.638, 14.5, "η ≥ 0.65: not a fluid-branch\nmeasurement (see text) →", color=MUTED, fontsize=7.5, ha="right", va="center")

if full:
    ax.errorbar([p[0] for p in full], [p[1] for p in full], yerr=[p[2] for p in full], fmt="o", ms=5, color=BLUE,
                mec=SURF, mew=1, capsize=2, lw=1.2, zorder=3, label="route A: N = 100 (50 per compartment), r = 0.5, H = 10, L₀ = 3.93 σ/η — the compartment shrinks with η\n(refit 2026-09-09, through the origin, eligible trajectories only; error bars = fit error)")
else:
    ax.plot([p[0] for p in routeA], [p[1] for p in routeA], "o", ms=7, color=BLUE, mec=SURF, mew=1, zorder=3,
            label="route A, N = 100 (50 per compartment), r = 0.5, H = 10 — 6 of 33 densities on record")
    ax.plot([p[0] for p in routeB], [p[1] for p in routeB], "o", ms=7, mfc="none", mec=ORANGE, mew=1.8, zorder=3,
            label="route B, N = 100, L0 = 20 fixed, r swept (slit confinement)")

ax.set_xlim(0, 0.76); ax.set_yscale("log"); ax.set_ylim(1.2, 42)
ax.set_yticks([1.414, 2, 3, 5, 10, 20, 40]); ax.set_yticklabels(["√2", "2", "3", "5", "10", "20", "40"])
ax.set_xlabel("packing fraction  η", color=INK2, fontsize=10)
ax.set_ylabel("speed of sound  c_s   [√(k_BT/m)]", color=INK2, fontsize=10)
ax.legend(loc="upper left", fontsize=8.5, frameon=False, labelcolor=INK2)
ax.set_title("Speed of sound of the 2D hard-disk fluid: measurement vs equation of state", loc="left", fontsize=12, color=INK, fontweight="semibold")

# ---- inset: ideal-gas end
ins = ax.inset_axes([0.31, 0.47, 0.32, 0.35])
ins.set_facecolor(SURF); ins.grid(True, color=GRID, lw=0.7); ins.set_axisbelow(True)
for s in ("top", "right"): ins.spines[s].set_visible(False)
ins.tick_params(colors=INK2, labelsize=8)
ez = np.linspace(1e-4, 0.1, 200)
ins.plot(ez, cs_from_Z(KR, ez), "-", color=INK2, lw=1.6)
ins.plot(ez, math.sqrt(2) * (1 + 2 * ez), ls=(0, (1, 1.5)), color=INK2, lw=1.3)
ins.axhline(math.sqrt(2), color=MUTED, lw=1.2, ls=(0, (2, 2)))
ins.plot([0], [math.sqrt(2)], "*", ms=11, color=INK, zorder=4, clip_on=False)
pts = [p for p in (full or routeA) if p[0] <= 0.1]
if pts:
    ins.plot([p[0] for p in pts], [p[1] for p in pts], "o", ms=6, color=BLUE, mec=SURF, mew=1, zorder=3)
ins.set_xlim(0, 0.1); ins.set_ylim(1.38, 1.80)
ins.set_title("ideal-gas end, η ≤ 0.1", fontsize=8.5, color=INK, loc="left")
ins.text(0.045, 1.40, "c_s = √2 (1 + 2η)  first order, exact", fontsize=7.2, color=INK2, va="bottom", ha="left")
ins.text(0.098, cs_from_Z(KR, np.array([0.098]))[0] + 0.01, "KR", fontsize=7.5, color=INK2, ha="right", va="bottom")
ins.set_xlabel("η", color=INK2, fontsize=8); ins.set_ylabel("c_s", color=INK2, fontsize=8)

fig.text(0.98, 0.005,
         "EDMD default core · kT = m = σ = 1 at t = 0 · ν = c_s K/(2π L_eff), L_eff = L₀ − 2r, cot K = (M/2Nm)K · 9 masses × 25 repeats per η · routeA_refit_cs_vs_eta_20260909.csv",
         ha="right", va="bottom", fontsize=7.2, color=MUTED)
fig.tight_layout(rect=(0, 0.03, 1, 1))
fig.savefig(OUT + ".png", dpi=200); fig.savefig(OUT + ".pdf")
print("wrote", OUT, "full sweep:", bool(full))
