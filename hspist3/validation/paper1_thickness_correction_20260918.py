#!/usr/bin/env python3
"""##CHRIS 2026-09-18: Paper 1's c_s(eta) with the divider thickness in L_eff. Analysis only.

The speed-of-sound campaigns all ran with --wall-thickness=0.05 (proved on 2026-09-18 by
reproducing a stored A1 v2 trajectory bit for bit: 6420/6420 sample times identical at t = 0.05,
1/6420 at t = 1.0). The estimator uses L_eff = L0 - 2r, which is the zero-thickness geometry. With a
divider of thickness t centred in the box, the disk centres in one compartment span r to
L0 - t/2 - r, so

    L_eff = L0 - 2r - t/2,      t/2 = 0.025 sigma.

c_s is the slope of nu against x_M = K(alpha) / (2 pi L_eff), so it is exactly proportional to L_eff
and the correction is one multiplication per density -- no trajectory is re-analysed:

    c_s_corrected = c_s * (L0 - 2r - t/2) / (L0 - 2r).

It is one-signed (every c_s moves DOWN) and largest where the box is smallest, i.e. at high eta.
Writes the corrected figure alongside the published one; it does not overwrite anything.
"""
import csv, math, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import tests_20260913 as T
import plot_speed_of_sound_edmd as sos

P = T.PLOTS
OUT = os.path.join(P, "260918_cs_vs_eta_thickness_corrected")
T_WALL, RDISK = 0.05, 0.5
HALF = 0.5 * T_WALL


def factor(L0):
    """How much the divider thickness takes off c_s at this box length."""
    lo = L0 - 2 * RDISK
    return (lo - HALF) / lo


def cs_of(Z, eta):
    h = 1e-5
    return sos.cs_adiabatic_2d_monatomic(Z(eta), (Z(eta + h) - Z(eta - h)) / (2 * h), eta, kbt=1, m=1)


a1 = []
for r in csv.DictReader(open(os.path.join(P, "260914_A1v2_final_cs_vs_eta.csv"))):
    eta, L0 = float(r["eta"]), float(r["L0"])
    c, s = float(r["c_s"]), float(r["c_s_scatter_mass"])
    f = factor(L0)
    a1.append((eta, L0, c, c * f, s, f))

# A2: same correction, per (eta, N) cell, using that cell's own L0
a2 = {}
for r in csv.DictReader(open(os.path.join(P, "260917_A2_cs_per_mass.csv"))):
    key = (float(r["eta"]), int(r["N"]))
    a2.setdefault(key, []).append((float(r["L0"]), float(r["c_s_mass"])))

print("### Paper 1 with L_eff = L0 - 2r - t/2, t = 0.05 sigma\n")
print("| η | L₀ | c_s published | c_s corrected | shift [%] | dev. from KR before | after |")
print("|---|---|---|---|---|---|---|")
for eta, L0, c, cc, s, f in a1:
    if eta in (0.006545, 0.026180, 0.078540, 0.196350, 0.392699, 0.523599, 0.600000, 0.650003, 0.700000) \
            or abs(eta - 0.3) < 0.005:
        kr = float(cs_of(sos.Z_kolafa_rottner_2006, np.array([eta]))[0]) if eta <= 0.69 else float("nan")
        d0 = 100 * (c - kr) / kr if kr == kr else float("nan")
        d1 = 100 * (cc - kr) / kr if kr == kr else float("nan")
        print(f"| {eta:.6f} | {L0:.2f} | {c:.4f} | {cc:.4f} | {100 * (f - 1):+.3f} | "
              f"{d0:+.2f} % | {d1:+.2f} % |" if kr == kr else
              f"| {eta:.6f} | {L0:.2f} | {c:.4f} | {cc:.4f} | {100 * (f - 1):+.3f} | — | — |")

fig, (ax, axd) = plt.subplots(2, 1, figsize=(10.0, 8.4), sharex=True,
                              gridspec_kw={"height_ratios": [2.2, 1]})
sos.add_eta_regime_shading(ax, x_max=0.78)
for t_ in list(ax.texts):
    if t_.get_rotation() == 90 and t_.get_position()[0] < 0.05:
        t_.set_position((t_.get_position()[0], 0.42))
b = list(axd.texts); sos.add_eta_regime_shading(axd, x_max=0.78)
for t_ in list(axd.texts):
    if t_ not in b:
        t_.remove()

e = np.linspace(0.001, 0.78, 600); ek = e[e <= 0.69]
ax.plot(ek, cs_of(sos.Z_kolafa_rottner_2006, ek), "-", color="#e34948", lw=2.4,
        label="Kolafa–Rottner 2006 (valid to η ≈ 0.69)")
ax.axhline(math.sqrt(2), color="0.45", lw=1.0, ls=":", label="ideal-gas limit η → 0:  c_s = √2")
E = np.array([p[0] for p in a1]); C = np.array([p[2] for p in a1])
CC = np.array([p[3] for p in a1]); S = np.array([p[4] for p in a1])
ax.errorbar(E, C, yerr=S, fmt="o", color="0.65", ms=4, lw=0, capsize=2,
            label="published, L$_{eff}$ = L$_0$ − 2r")
ax.errorbar(E, CC, yerr=S, fmt="o-", color="#2a78d6", ms=4.5, lw=1.0, capsize=2.5,
            label="corrected, L$_{eff}$ = L$_0$ − 2r − t/2  (t = 0.05 σ)")
ax.set_ylim(0, 22); ax.set_ylabel("Speed of sound  c_s  [√(k_BT/m)]", fontsize=11.5)
ax.set_title("Paper 1 with the divider thickness carried into L$_{eff}$", fontsize=12.5)
ax.grid(True, ls=":", alpha=0.6)
ax.legend(loc="upper left", bbox_to_anchor=(0.0, 0.91), fontsize=9, framealpha=0.95)

k = E <= 0.69
KR = np.array([float(cs_of(sos.Z_kolafa_rottner_2006, np.array([x]))[0]) for x in E[k]])
axd.axhline(0, color="#e34948", lw=1.6)
axd.axhspan(-0.5, 0.5, color="#e1e0d9", alpha=0.7, zorder=0)
axd.plot(E[k], 100 * (C[k] - KR) / KR, "o", color="0.65", ms=4, label="published")
axd.errorbar(E[k], 100 * (CC[k] - KR) / KR, yerr=100 * S[k] / KR, fmt="o-", color="#2a78d6",
             ms=4.5, lw=1.0, capsize=2.5, label="corrected")
axd.set_ylim(-3, 8); axd.set_xlim(0, 0.78)
axd.set_xlabel("Packing fraction  η", fontsize=11.5)
axd.set_ylabel("deviation from Kolafa–Rottner [%]", fontsize=10)
axd.grid(True, ls=":", alpha=0.6); axd.legend(fontsize=9, loc="upper left", framealpha=0.95)
axd.text(0.005, -2.7, "±0.5 % band", fontsize=8, color="0.35")
fig.text(0.99, 0.004, "correction is exact and multiplicative: c_s ∝ L_eff, so c_s → c_s (L₀ − 1 − 0.025)/(L₀ − 1) · "
                      "no trajectory re-analysed · data: 260914_A1v2_final_cs_vs_eta.csv",
         ha="right", va="bottom", fontsize=7, color="0.4")
fig.tight_layout(rect=(0, 0.015, 1, 1))
for ext in ("png", "pdf"):
    fig.savefig(f"{OUT}.{ext}", dpi=200)
print(f"\nwritten {OUT}.png/.pdf")

print("\n### A2 finite-size points, same correction\n")
print("| η | N | L₀ | c_s published | corrected | shift [%] |")
print("|---|---|---|---|---|---|")
for (eta, N), rows in sorted(a2.items()):
    L0 = rows[0][0]
    xs = np.array([T.k_root(M / N) / (2 * math.pi * (L0 - 1.0)) for M in (50, 200, 500, 1000, 2000)][:len(rows)])
    c = float(np.mean([v for _, v in rows]))
    f = factor(L0)
    if N in (900, 1600) and eta in (0.02, 0.05, 0.10, 0.30, 0.50, 0.60, 0.65):
        print(f"| {eta:.2f} | {N} | {L0:.2f} | {c:.4f} | {c * f:.4f} | {100 * (f - 1):+.3f} |")
