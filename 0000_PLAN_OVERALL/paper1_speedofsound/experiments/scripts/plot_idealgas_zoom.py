#!/usr/bin/env python3
"""Ideal-gas zoom of the pressure validation: Z vs eta for eta <= 0.1.

Reads analysis/tables_ABC.md from the 2026-09-07 pressure campaign (per-cell
table B and the 1/sqrt(N) fit table) and draws three panels:
  (a) Z vs eta with the exact anchors: ideal gas Z=1, second virial Z=1+2eta,
      and the Kolafa-Rottner 2006 fit; points per N and the N->inf stars
  (b) deviation from Kolafa-Rottner in percent
  (c) (Z-1)/eta, which must go to exactly 2 (= B2 for hard disks) as eta -> 0
Every number is parsed from the table file; nothing is typed in.
"""
import re, math, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TABLE = sys.argv[1]
OUT = sys.argv[2]

# ---- palette: reference categorical slots 1-3 (validated all-pairs), text/ink tokens
C = {400: "#2a78d6", 900: "#eb6834", 1600: "#1baf7a"}
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e6e5e1"

KR_COEF = {0: 1.0, 1: 2.0, 2: 1.12801775, 3: 0.00181895291, 4: -0.0526134737,
           5: 0.0504960168, 6: -0.0325537792, 7: 0.0134578632, 8: 0.00140888182,
           9: -0.00834273601, 10: 0.00694127367, 11: -0.00262254723,
           12: 0.000355746352, 22: -5.24672938e-9}
def KR(eta):
    x = np.asarray(eta) / (1 - np.asarray(eta))
    return sum(c * x ** p for p, c in KR_COEF.items())

# ---- parse tables_ABC.md
txt = open(TABLE, encoding="utf-8").read()
sectB = txt.split("## B.")[1].split("Fit Z(N)")[0]
cells = []   # (eta, N, nseeds, Z, sd, bsem)
for l in sectB.splitlines():
    if l.startswith("| 0."):
        c = [x.strip() for x in l.strip("|").split("|")]
        cells.append((float(c[0]), int(c[1]), int(c[2]), float(c[3]), float(c[4]), float(c[5])))
fit_sect = txt.split("Z_pair_inf | sigma")[1]
fits = {}
for l in fit_sect.splitlines()[1:]:
    if not l.startswith("|"):
        if fits: break
        continue
    c = [x.strip() for x in l.strip("|").split("|")]
    try:
        fits[float(c[0])] = (float(c[2]), float(c[3]))
    except ValueError:
        pass

ZOOM = 0.105
cells = [c for c in cells if c[0] <= ZOOM]
etas = sorted({c[0] for c in cells})
def err(c):
    e, N, n, Z, sd, bsem = c
    ssem = sd / math.sqrt(n) if n > 1 else 0.0
    return max(bsem, ssem)

fig, axs = plt.subplots(1, 3, figsize=(13.5, 4.6))
fig.patch.set_facecolor("#fcfcfb")
for ax in axs:
    ax.set_facecolor("#fcfcfb")
    ax.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=INK2, labelsize=9)
    ax.set_xlim(0, ZOOM)
    ax.set_xlabel("packing fraction  η", color=INK2, fontsize=10)

eg = np.linspace(1e-4, ZOOM, 400)
# small x-offsets so the three N do not sit on top of each other
off = {400: -0.0022, 900: 0.0, 1600: 0.0022}

# (a) Z vs eta
ax = axs[0]
ax.plot(eg, np.ones_like(eg), ls=(0, (2, 2)), color=MUTED, lw=1.4)
ax.plot(eg, 1 + 2 * eg, ls=(0, (1, 1.5)), color=INK2, lw=1.4)
ax.plot(eg, KR(eg), "-", color=INK2, lw=1.6)
for N in (400, 900, 1600):
    pts = [c for c in cells if c[1] == N]
    ax.errorbar([c[0] + off[N] for c in pts], [c[3] for c in pts], yerr=[err(c) for c in pts],
                fmt="o", ms=5, color=C[N], mfc=C[N], mec="#fcfcfb", mew=1, capsize=2, lw=1.2, zorder=3)
ax.plot([e for e in etas], [fits[e][0] for e in etas], "*", ms=12, color=INK, zorder=4)
ax.set_ylabel("Z = P / (ρ k_B T)", color=INK2, fontsize=10)
ax.set_ylim(0.99, 1.26)
ax.text(0.101, 1.0, "ideal gas  Z = 1", color=MUTED, fontsize=8.5, ha="right", va="bottom")
ax.text(0.062, 1.108, "Z = 1 + 2η  (second virial, exact)", color=INK2, fontsize=8.5, ha="left", va="top", rotation=31)
ax.text(0.099, 1.212, "Kolafa–Rottner 2006", color=INK2, fontsize=8.5, ha="right", va="top")
ax.set_title("(a)  the simulation lands on the exact low-density law", loc="left", fontsize=10.5, color=INK)

# (b) deviation from KR in %
ax = axs[1]
ax.axhline(0, color=INK2, lw=1.2)
for N in (400, 900, 1600):
    pts = [c for c in cells if c[1] == N]
    kr = KR([c[0] for c in pts])
    ax.errorbar([c[0] + off[N] for c in pts], [100 * (c[3] - k) / k for c, k in zip(pts, kr)],
                yerr=[100 * err(c) / k for c, k in zip(pts, kr)],
                fmt="o", ms=5, color=C[N], mfc=C[N], mec="#fcfcfb", mew=1, capsize=2, lw=1.2, zorder=3,
                label=f"N = {N}")
kr_inf = KR(etas)
ax.errorbar(etas, [100 * (fits[e][0] - k) / k for e, k in zip(etas, kr_inf)],
            yerr=[100 * fits[e][1] / k for e, k in zip(etas, kr_inf)],
            fmt="*", ms=12, color=INK, capsize=2, lw=1.0, zorder=4, label="N → ∞ (1/√N fit)")
ax.set_ylabel("Z − Z_KR  [%]", color=INK2, fontsize=10)
ax.set_ylim(-0.35, 0.35)
ax.set_yticks([-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3])
ax.axhspan(-0.2, 0.2, color=GRID, alpha=0.6, zorder=0)
ax.text(0.002, 0.21, "±0.2 % band", color=MUTED, fontsize=8.5, va="bottom")
ax.legend(loc="lower right", fontsize=8.5, frameon=False, labelcolor=INK2)
ax.set_title("(b)  deviation from Kolafa–Rottner", loc="left", fontsize=10.5, color=INK)

# (c) slope test
ax = axs[2]
ax.axhline(2.0, color=INK2, lw=1.2)
ax.plot(eg, (KR(eg) - 1) / eg, "-", color=INK2, lw=1.6)
for N in (400, 900, 1600):
    pts = [c for c in cells if c[1] == N]
    ax.errorbar([c[0] + off[N] for c in pts], [(c[3] - 1) / c[0] for c in pts], yerr=[err(c) / c[0] for c in pts],
                fmt="o", ms=5, color=C[N], mfc=C[N], mec="#fcfcfb", mew=1, capsize=2, lw=1.2, zorder=3)
ax.errorbar(etas, [(fits[e][0] - 1) / e for e in etas], yerr=[fits[e][1] / e for e in etas],
            fmt="*", ms=12, color=INK, capsize=2, lw=1.0, zorder=4)
ax.set_ylabel("(Z − 1) / η", color=INK2, fontsize=10)
ax.set_ylim(1.8, 2.5)
ax.text(0.101, 2.005, "B₂ = 2  (exact, η → 0)", color=INK2, fontsize=8.5, ha="right", va="bottom")
ax.text(0.06, 2.28, "Kolafa–Rottner", color=INK2, fontsize=8.5, ha="left", va="bottom", rotation=17)
ax.set_title("(c)  the slope at η → 0 must be exactly 2", loc="left", fontsize=10.5, color=INK)

fig.suptitle("Hard-disk EDMD pressure validation, ideal-gas end (η ≤ 0.1)", x=0.01, ha="left",
             fontsize=12.5, color=INK, fontweight="semibold")
fig.text(0.99, 0.005,
         "hard-wall square box · r = 0.5 · kT ≈ 1 · N = 400 / 900 / 1600, 5 / 4 / 3 seeds · 400 t.u. equilibration + 30 blocks × 20 t.u. · "
         "pair-virial route · campaign 00_pressure_validation_20260907, analysis/tables_ABC.md",
         ha="right", va="bottom", fontsize=7.5, color=MUTED)
fig.tight_layout(rect=(0, 0.03, 1, 0.95))
fig.savefig(OUT + ".png", dpi=200)
fig.savefig(OUT + ".pdf")
print("wrote", OUT + ".png/.pdf", "cells:", len(cells))
