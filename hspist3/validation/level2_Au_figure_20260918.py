#!/usr/bin/env python3
"""##CHRIS 2026-09-18: A(u) = (W - W_qs)/u^2 from u = 0.005 to 10, with the two predicted plateaus.

Every point is scaled to a common gas travel of 3.93 sigma, because the excess work is proportional
to the travel in the single-hit limit (audit Eq. 7) -- without that the fast runs, which travel
1 or 7 sigma, cannot share an axis with the slow ones.

  acoustic plateau  : A = 25.2, measured (760 trajectories, step protocol)
  single-hit plateau: A = 2 N_s m (3.93 / L) = 10.1, audit Eq. 7, no free parameter
  steady flow       : A = N_s m / 6 = 8.3, audit Eq. 15/18 -- what the RAMP protocol leaves behind
"""
import glob, math, os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tests_20260913 as T
from paper2_geometry_fix_20260918 import w_qs, ET, NS
from paper2_ramp_fast_20260918 import path_points, cell, RAMP, FAST, FAST7, L_I, DIV_FACE

OUT = os.path.join(T.PLOTS, "260918_level2_A_of_u")
DX0 = 3.93

pts = path_points()
WQS, DWQS, _, _ = w_qs(pts, L_I, 74.57 - DIV_FACE)


def gas_travel(d):
    """Gas compression = piston displacement - 0.25: the piston parks outside the box wall, so its
    first quarter sigma sweeps no gas. Measured per cell, never assumed (this is what was wrong in
    the first version of this figure, where the fast sets were scaled by the piston displacement)."""
    tr = sorted(glob.glob(f"{d}/tr_*.csv") + glob.glob(f"{d}/tr_*.csv.gz"))[0]
    t = pd.read_csv(tr, usecols=["PistonR_x_sigma", "PistonR_v"], low_memory=False)
    v = np.abs(t["PistonR_v"].to_numpy(float)); x = t["PistonR_x_sigma"].to_numpy(float)
    mv = np.nonzero(v > 1e-12)[0]
    return abs(x[mv[-1]] - x[mv[0]]) - 0.25


def series(paths, dx=None):
    out = []
    for u, d in paths:
        if not os.path.exists(f"{d}/summary.csv"):
            continue
        m, s, sd, n = cell(d)
        g = dx if dx is not None else gas_travel(d)
        wq, dwq, _, _ = w_qs(pts, L_I, L_I - g)
        A = (m - wq) / u ** 2 * (DX0 / g)
        dA = math.hypot(s, dwq) / u ** 2 * (DX0 / g)
        out.append((u, A, dA, g, m, s, 2 * NS * (g / L_I) * u ** 2))
    return np.array(out)

slow = series([(u, f"{ET}/level2_slope_20260917/u{u:.2f}") for u in
               (0.005, 0.01, 0.02, 0.03, 0.05, 0.10, 0.15, 0.20)], DX0)
step3 = series([(u, f"{RAMP}/step_u{u:.2f}") for u in (0.05, 0.10, 0.20)], DX0)
ramp3 = series([(u, f"{RAMP}/u{u:.2f}") for u in (0.05, 0.10, 0.20)], DX0)
fast1 = series([(u, f"{FAST}/u{u}") for u in (3, 5, 10)])
fast7 = series([(u, f"{FAST7}/u{u}") for u in (3, 5, 10)])
scan = series([(10, f"{ET}/level2_fastdx_20260918/dx{d}") for d in (2, 3, 5)])

fig, (ax, axr) = plt.subplots(1, 2, figsize=(12.4, 5.2), gridspec_kw={"width_ratios": [1.55, 1]})

# ---- left: A(u) across five decades of speed, every point scaled to a common gas travel
ax.axhline(25.2, color="tab:blue", lw=1.4, ls="--", label="acoustic plateau, A = 25.2 (measured, step)")
ax.axhline(2 * NS * DX0 / 38.75, color="tab:red", lw=1.4, ls="-.",
           label=f"single-hit limit, 2N$_s$m Δx/L = {2 * NS * DX0 / 38.75:.1f} (audit Eq. 7)")
ax.axhline(NS / 6, color="tab:green", lw=1.2, ls=":", label="steady flow, N$_s$m/6 = 8.3 (audit Eq. 15/18)")
for d, st in ((slow, dict(m="o", color="k", ms=4.5, mfc="k", label="step, Δx = 3.93 σ (760 runs)")),
              (step3, dict(m="s", color="tab:blue", ms=5, mfc="white", label="step, this batch (60 seeds)")),
              (ramp3, dict(m="D", color="tab:green", ms=5, mfc="tab:green", label="ramp over 2L/c$_s$ (60 seeds)")),
              (fast1, dict(m="^", color="tab:orange", ms=6, mfc="tab:orange", label="fast, gas Δx = 0.75 σ")),
              (scan, dict(m="P", color="tab:purple", ms=6, mfc="tab:purple", label="fast, Δx = 2.1–5.1 σ (new)")),
              (fast7, dict(m="v", color="tab:red", ms=6, mfc="white", label="fast, gas Δx = 7.0 σ"))):
    if len(d):
        ax.errorbar(d[:, 0], d[:, 1], yerr=d[:, 2], elinewidth=1.0, capsize=2.5, lw=0,
                    marker=st.pop("m"), **st)
ax.set_xscale("log"); ax.set_xlabel("piston speed  u  [σ/τ]   (c$_s$ ≈ 1.79)")
ax.set_ylabel("A = (⟨W⟩ − W$_{qs}$)/u²,  scaled to Δx = 3.93 σ  [k$_B$T τ²/σ²]")
ax.set_ylim(-20, 60); ax.set_xlim(3e-3, 20)
ax.axvspan(1.79, 20, color="0.92", zorder=0)
ax.text(3.0, 52, "u > c$_s$", fontsize=9, color="0.4")
ax.grid(True, ls=":", alpha=0.6)
ax.set_title("Excess work per u², from the acoustic to the single-hit regime", fontsize=11)
ax.legend(fontsize=7.5, loc="upper left", framealpha=0.95, ncol=2)

# ---- right: the new travel scan. Eq. 7 has no free parameter, so this ratio should be 1.
allfast = np.vstack([r for r in (fast1, scan, fast7) if len(r)])
for u, col, mk in ((3, "tab:brown", "^"), (5, "tab:cyan", "s"), (10, "tab:purple", "o")):
    k = allfast[:, 0] == u
    if k.sum():
        d = allfast[k][np.argsort(allfast[k][:, 3])]
        axr.errorbar(d[:, 3], d[:, 4] / d[:, 6], yerr=d[:, 5] / d[:, 6], marker=mk, color=col,
                     ms=6, lw=1.1, capsize=3, label=f"u = {u:g}")
axr.axhline(1.0, color="tab:red", lw=1.6, ls="-.")
axr.set_xlabel("gas compression  Δx  [σ]")
axr.set_ylabel("⟨W⟩ / 2N$_s$m(Δx/L)u²")
axr.set_xlim(0, 7.8); axr.set_ylim(0.935, 1.40)
axr.grid(True, ls=":", alpha=0.6)
axr.set_title("The single-hit limit against travel (new, 2026-09-18)", fontsize=11)
axr.legend(fontsize=8.5, loc="upper right", framealpha=0.95)
axr.text(0.25, 0.947, "the excess falls toward Eq. 7 as travel and speed grow; it is not explained",
         fontsize=8, color="0.3", ha="left", va="bottom")
axr.text(5.9, 1.012, "Eq. 7, no free parameter", fontsize=8.5, color="tab:red", ha="right")

fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"{OUT}.{ext}", dpi=200)
print(f"written {OUT}.png/.pdf")
for nm, d in (("slow/step", slow), ("step batch", step3), ("ramp", ramp3),
              ("fast 0.75", fast1), ("scan u=10", scan), ("fast 7.0", fast7)):
    if len(d):
        print(f"  {nm:11s} " + "  ".join(f"u={r[0]:g} Δx={r[3]:.2f}: A={r[1]:6.1f}±{r[2]:.1f} ratio={r[4]/r[6]:.3f}" for r in d))
