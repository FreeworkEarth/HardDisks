#!/usr/bin/env python3
"""##CHRIS 2026-09-16: dilute-end zoom, eta <= 0.15, N = 100 against N = 900 and N = 1600.
Analysis only. Same estimator everywhere (largest FFT bin at f >= nu_pred/2.5); error bars are the
1 sigma scatter of the per-mass c_s. Note A2 has only ONE density below 0.15 (eta = 0.10), so the
larger sizes appear as single points, not curves -- that is the honest state of the data."""
import os, sys, csv, math
from collections import defaultdict
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import tests_20260913 as T
import plot_speed_of_sound_edmd as sos

P = T.PLOTS
OUT = os.path.join(P, sys.argv[2] if len(sys.argv) > 2 else "260916_cs_vs_eta_lowdensity_zoom")
XMAX = 0.15


def cs_of(Z, eta):
    h = 1e-5
    return sos.cs_adiabatic_2d_monatomic(Z(eta), (Z(eta + h) - Z(eta - h)) / (2 * h), eta, kbt=1, m=1)


a1 = [(float(r["eta"]), float(r["c_s"]), float(r["c_s_scatter_mass"]))
      for r in csv.DictReader(open(os.path.join(P, "260914_A1v2_final_cs_vs_eta.csv")))]
a1 = sorted(p for p in a1 if p[0] <= XMAX)
by = defaultdict(list)
for r in csv.DictReader(open(os.path.join(P, (sys.argv[1] if len(sys.argv) > 1 else "260916_A2_cs_per_mass.csv")))):
    by[(float(r["eta"]), int(r["N"]))].append((float(r["nu_mean"]), float(r["c_s_mass"])))
a2 = defaultdict(list)
for (eta, N), v in by.items():
    if eta > XMAX:
        continue
    nu = np.array([q[0] for q in v]); cm = np.array([q[1] for q in v]); x = nu / cm
    a2[N].append((eta, float((x * nu).sum() / (x * x).sum()), float(cm.std(ddof=1)), len(v)))

# ##CHRIS 2026-09-17: same regime shading as 260914_cs_vs_eta, from the one definition in
# plot_speed_of_sound_edmd.add_eta_regime_shading (boundaries 0.02 / 0.20 / 0.50 / 0.70 / 0.716 / 0.72).
# Labels on the top panel only; the deviation panel gets the bands without text.
def shade(axis, xmax, labels=True):
    before = list(axis.texts)
    sos.add_eta_regime_shading(axis, x_max=xmax)
    if not labels:
        for t in list(axis.texts):
            if t not in before:
                t.remove()

fig, (ax, axd) = plt.subplots(2, 1, figsize=(9.5, 8.0), sharex=True, gridspec_kw={"height_ratios": [2.2, 1]})
shade(ax, XMAX); shade(axd, XMAX, labels=False)
# the ideal-gas band is only 0.02 wide, so the shared helper prints it rotated and small; in a figure
# that is entirely about the dilute end it is worth the same top label as its neighbour.
for _t in list(ax.texts):
    if _t.get_rotation() == 90:
        _t.remove()
ax.text(0.010, 0.975, "ideal-gas\nlimit", transform=ax.get_xaxis_transform(),
        ha="center", va="top", fontsize=8, color="#333333")
e = np.linspace(1e-4, XMAX, 400)
ax.plot(e, cs_of(sos.Z_kolafa_rottner_2006, e), "-", color="#e34948", lw=2.4, label="Kolafa–Rottner 2006")
ax.plot(e, cs_of(sos.Z_spt_eos, e), "--", color="#eda100", lw=1.3, label="SPT equation of state")
ax.plot(e, cs_of(sos.Z_henderson_eos, e), "-.", color="#1baf7a", lw=1.3, label="Henderson (a = 0.125)")
ax.plot([0], [math.sqrt(2)], "*", ms=15, color="black", zorder=6, clip_on=False,
        label="exact ideal gas, η → 0:  c_s = √2 = 1.4142")
ax.errorbar([p[0] for p in a1], [p[1] for p in a1], yerr=[p[2] for p in a1], fmt="o-", color="#2a78d6",
            ms=5, lw=1.1, capsize=2.5, zorder=4, label="A1 v2, N = 100 · 10 densities, 9 masses × 25 seeds, 200 periods")
sty = {900: ("#eb6834", True, "A2, N = 900"), 1600: ("#4a3aa7", False, "A2, N = 1600")}
for N in (900, 1600):
    col, filled, lab = sty[N]
    pts = sorted(a2[N])
    ax.errorbar([p[0] for p in pts], [p[1] for p in pts], yerr=[p[2] for p in pts], fmt="o", color=col,
                ms=7, lw=1.5, capsize=3.5, mfc=col if filled else "white", mew=1.6, zorder=5,
                label=f"{lab} · 5 masses × 10–35 seeds, 50 periods (η ≤ 0.05) / 37.5 (η = 0.10)")
ax.set_xlim(0, XMAX); ax.set_ylim(1.38, 1.95)
ax.set_ylabel("Speed of sound  c_s  [√(k_BT/m)]", fontsize=11.5)
ax.set_title("The dilute end: one system size measured densely, two larger sizes at η = 0.10", fontsize=12.5)
ax.grid(True, ls=":", alpha=0.6)
ax.legend(loc="upper left", bbox_to_anchor=(0.0, 0.93), fontsize=8.5, framealpha=0.95)

kr1 = np.array([float(cs_of(sos.Z_kolafa_rottner_2006, np.array([p[0]]))[0]) for p in a1])
axd.axhline(0, color="#e34948", lw=1.6); axd.axhspan(-0.5, 0.5, color="#e1e0d9", alpha=0.7, zorder=0)
axd.errorbar([p[0] for p in a1], 100 * (np.array([p[1] for p in a1]) - kr1) / kr1,
             yerr=100 * np.array([p[2] for p in a1]) / kr1, fmt="o-", color="#2a78d6", ms=5, lw=1.1, capsize=2.5)
for N in (900, 1600):
    col, filled, _ = sty[N]
    pts = sorted(a2[N])
    k = np.array([float(cs_of(sos.Z_kolafa_rottner_2006, np.array([p[0]]))[0]) for p in pts])
    axd.errorbar([p[0] for p in pts], 100 * (np.array([p[1] for p in pts]) - k) / k,
                 yerr=100 * np.array([p[2] for p in pts]) / k, fmt="o", color=col, ms=7, lw=1.5, capsize=3.5,
                 mfc=col if filled else "white", mew=1.6, zorder=5)
axd.set_xlabel("Packing fraction  η", fontsize=11.5)
axd.set_ylabel("deviation from Kolafa–Rottner [%]", fontsize=10)
axd.set_ylim(-1.2, 2.0); axd.grid(True, ls=":", alpha=0.6)
axd.text(0.0015, 0.58, "±0.5 % band", fontsize=8, color="0.35", va="bottom")
axd.annotate("A2 has only one density below 0.15 (η = 0.10),\nso the larger sizes are single points, not curves",
             xy=(0.10, -0.55), xytext=(0.055, -1.05), fontsize=7.8, color="#52514e",
             arrowprops=dict(arrowstyle="->", color="#898781", lw=0.9))
fig.text(0.99, 0.004, "same estimator everywhere: largest FFT bin at f ≥ ν_pred/2.5 · error bars = 1σ scatter of per-mass c_s · "
                      "data: 260914_A1v2_final_cs_vs_eta.csv + " + (sys.argv[1] if len(sys.argv) > 1 else "260916_A2_cs_per_mass.csv") + "",
         ha="right", va="bottom", fontsize=7, color="0.4")
fig.tight_layout(rect=(0, 0.015, 1, 1))
for ext in ("png", "pdf"):
    fig.savefig(f"{OUT}.{ext}", dpi=200)
print("| η | N | c_s | ± scatter | KR | dev [%] |")
print("|---|---|---|---|---|---|")
for p in a1:
    kr = float(cs_of(sos.Z_kolafa_rottner_2006, np.array([p[0]]))[0])
    print(f"| {p[0]:.6f} | 100 | {p[1]:.4f} | {p[2]:.4f} | {kr:.4f} | {100*(p[1]-kr)/kr:+.2f} |")
for N in (900, 1600):
    for eta, c, s, n in sorted(a2[N]):
        kr = float(cs_of(sos.Z_kolafa_rottner_2006, np.array([eta]))[0])
        print(f"| {eta:.6f} | {N} | {c:.4f} | {s:.4f} | {kr:.4f} | {100*(c-kr)/kr:+.2f} |")
print("\nwrote", OUT + ".png/.pdf")
