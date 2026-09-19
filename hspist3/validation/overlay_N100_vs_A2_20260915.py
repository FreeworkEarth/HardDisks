#!/usr/bin/env python3
"""##CHRIS 2026-09-15: c_s(eta) at N = 100 (A1 v2) with the A2 points at N = 900 and N = 1600 overlaid.
Analysis only: reads 260914_A1v2_final_cs_vs_eta.csv and 260916_A2_cs_per_mass.csv. Same estimator everywhere
(largest FFT bin at f >= nu_pred/2.5); error bars are the 1 sigma scatter of the per-mass c_s."""
import os, sys, csv, math
from collections import defaultdict
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import tests_20260913 as T
import plot_speed_of_sound_edmd as sos

P = T.PLOTS
A1 = os.path.join(P, "260914_A1v2_final_cs_vs_eta.csv")
A2 = os.path.join(P, "260916_A2_cs_per_mass.csv")
OUT = os.path.join(P, "260916_cs_vs_eta_N100_vs_A2")


def cs_of(Z, eta):
    h = 1e-5
    dZ = (Z(eta + h) - Z(eta - h)) / (2 * h)
    return sos.cs_adiabatic_2d_monatomic(Z(eta), dZ, eta, kbt=1, m=1)


a1 = [(float(r["eta"]), float(r["c_s"]), float(r["c_s_scatter_mass"])) for r in csv.DictReader(open(A1))]
by = defaultdict(list)
for r in csv.DictReader(open(A2)):
    by[(float(r["eta"]), int(r["N"]))].append((float(r["nu_mean"]), float(r["c_s_mass"])))
a2 = defaultdict(list)
for (eta, N), v in by.items():
    nu = np.array([q[0] for q in v]); csm = np.array([q[1] for q in v]); x = nu / csm
    a2[N].append((eta, float((x * nu).sum() / (x * x).sum()), float(csm.std(ddof=1)), len(v)))
for N in a2:
    a2[N].sort()

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

fig, (ax, axd) = plt.subplots(2, 1, figsize=(10.5, 8.2), sharex=True, gridspec_kw={"height_ratios": [2.3, 1]})
shade(ax, 0.78); shade(axd, 0.78, labels=False)
# the rotated "ideal-gas" label sits at the height this figure's legend occupies; drop it clear of it.
for _t in list(ax.texts):
    if _t.get_rotation() == 90 and _t.get_position()[0] < 0.05:
        _t.set_position((_t.get_position()[0], 0.42))
e = np.linspace(0.001, 0.78, 600); ek = e[e <= 0.69]
ax.plot(e, cs_of(sos.Z_spt_eos, e), "--", color="#eda100", lw=1.4, label="SPT equation of state")
ax.plot(e, cs_of(sos.Z_henderson_eos, e), "-.", color="#1baf7a", lw=1.4, label="Henderson (a = 0.125) equation of state")
ax.plot(ek, cs_of(sos.Z_kolafa_rottner_2006, ek), "-", color="#e34948", lw=2.4, label="Kolafa–Rottner 2006 (valid to η ≈ 0.69)")
ax.axhline(math.sqrt(2), color="0.45", lw=1.0, ls=":", label="ideal-gas limit η → 0:  c_s = √2")
ax.errorbar([p[0] for p in a1], [p[1] for p in a1], yerr=[p[2] for p in a1], fmt="o-", color="#2a78d6", ms=4.5, lw=1.1, capsize=2.5,
            label="A1 v2, N = 100 (50/50), 9 masses × 25 seeds, 200 periods")
sty = {900: ("o", "#eb6834", "A2, N = 900 (5 masses × 10–35 seeds, 37.5 periods)"),
       1600: ("o", "#4a3aa7", "A2, N = 1600 (5–6 masses × 10–35 seeds, 37.5 periods)")}
for N in (900, 1600):
    mk, col, lab = sty[N]
    ax.errorbar([p[0] for p in a2[N]], [p[1] for p in a2[N]], yerr=[p[2] for p in a2[N]], fmt=mk, color=col, ms=5, lw=1.4,
                capsize=3, mfc="white" if N == 1600 else col, mew=1.4, zorder=5, label=lab)
ax.set_ylim(0, 22); ax.set_ylabel("Speed of sound  c_s  [√(k_BT/m)]", fontsize=11.5)
ax.set_title("Speed of sound against packing fraction: one system size (N = 100) and two larger ones", fontsize=12.5)
ax.grid(True, ls=":", alpha=0.6)
ax.legend(loc="upper left", bbox_to_anchor=(0.0, 0.91), fontsize=8.5, framealpha=0.95)
ax.text(0.565, 2.2, "η ≥ 0.65: 6 σ compartment,\nstructure changes during\nmeasurement — not a\nfluid-branch value",
        fontsize=8.5, color="0.3", ha="left", va="bottom")
kr1 = cs_of(sos.Z_kolafa_rottner_2006, np.array([p[0] for p in a1 if p[0] <= 0.69]))
axd.axhline(0, color="#e34948", lw=1.6); axd.axhspan(-0.5, 0.5, color="#e1e0d9", alpha=0.7, zorder=0)
sub = [p for p in a1 if p[0] <= 0.69]
axd.errorbar([p[0] for p in sub], 100 * (np.array([p[1] for p in sub]) - kr1) / kr1,
             yerr=100 * np.array([p[2] for p in sub]) / kr1, fmt="o-", color="#2a78d6", ms=4.5, lw=1.1, capsize=2.5)
for N in (900, 1600):
    mk, col, _ = sty[N]
    pts = [p for p in a2[N] if p[0] <= 0.69]
    k = cs_of(sos.Z_kolafa_rottner_2006, np.array([p[0] for p in pts]))
    axd.errorbar([p[0] for p in pts], 100 * (np.array([p[1] for p in pts]) - k) / k,
                 yerr=100 * np.array([p[2] for p in pts]) / k, fmt=mk, color=col, ms=5, lw=1.4, capsize=3,
                 mfc="white" if N == 1600 else col, mew=1.4, zorder=5)
axd.set_ylim(-6, 8); axd.set_xlim(0, 0.78); axd.set_xlabel("Packing fraction  η", fontsize=11.5)
axd.set_ylabel("deviation from Kolafa–Rottner [%]", fontsize=10); axd.grid(True, ls=":", alpha=0.6)
axd.text(0.005, -5.5, "±0.5 % band", fontsize=8, color="0.35", va="bottom")
fig.text(0.99, 0.004, "same estimator everywhere: largest FFT bin at f ≥ ν_pred/2.5 · error bars = 1σ scatter of per-mass c_s · "
                      "data: 260914_A1v2_final_cs_vs_eta.csv + 260916_A2_cs_per_mass.csv",
         ha="right", va="bottom", fontsize=7, color="0.4")
fig.tight_layout(rect=(0, 0.015, 1, 1))
for ext in ("png", "pdf"):
    fig.savefig(f"{OUT}.{ext}", dpi=200)
print("A1 v2 has no density at exactly 0.10/0.30/0.50/0.60/0.65, so the N = 100 column below is A1 v2 at its OWN")
print("nearest density, and every deviation is against KR at that point's own eta.\n")
print("| A2 η | N = 900 | ± | dev KR [%] | N = 1600 | ± | dev KR [%] | nearest A1 v2 η | N = 100 | ± | dev KR [%] |")
print("|---|---|---|---|---|---|---|---|---|---|---|")
for eta, c9, s9, n9 in a2[900]:
    c16 = next((p for p in a2[1600] if abs(p[0] - eta) < 1e-9), None)
    a1p = min(a1, key=lambda p: abs(p[0] - eta))
    kr = float(cs_of(sos.Z_kolafa_rottner_2006, np.array([eta]))[0])
    kr1p = float(cs_of(sos.Z_kolafa_rottner_2006, np.array([a1p[0]]))[0])
    print(f"| {eta:.2f} | {c9:.4f} | {s9:.4f} | {100*(c9-kr)/kr:+.2f} | {c16[1]:.4f} | {c16[2]:.4f} | {100*(c16[1]-kr)/kr:+.2f} | "
          f"{a1p[0]:.6f} | {a1p[1]:.4f} | {a1p[2]:.4f} | {100*(a1p[1]-kr1p)/kr1p:+.2f} |")
print("\nwrote", OUT + ".png/.pdf")
