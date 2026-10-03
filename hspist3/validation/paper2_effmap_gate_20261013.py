#!/usr/bin/env python3
"""##CHRIS 2026-10-13: the position-settled gate (261010 sec. 3, pre-registered), replacing the KE_div gate.

    Delta = | <x_bar>_{last P_m} - <x_bar>_{previous P_m} | < 0.01 s_rev(k)

x_bar(t): the 8-seed mean divider trajectory (W0_x_sigma, written %.6f at 00ALLINONE.c:17131-17134).
s_rev: the reversible displacement of sec. 1.3. P_m: the mode period from the exact one-column eigen-equation
    M_s omega^2 = k + k_S K cot K,  K = omega L/c_s,  k_S = N m c_s^2/L^2,
with KR c_s at the settled (eta_f, T_f) of sec. 1.3 (controls: eta_0, T = 1).
Every one of the 36 cells gets PASS or FAIL; the six no-push controls calibrate the gate's noise floor.
Re-renders 261012_p2_effmap(.png/.pdf), 261012_p2_effmap_both and per-k panels 261012_p2_effmap_k<k> with the label.
No new runs; reads the recorded red_*.csv only.
"""
import json, math, os, sys
import numpy as np
from scipy.optimize import brentq
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import paper2_effmap_analysis_20261012 as A
from paper2_effmap_amend_20261012 import reversible, tau_r, Z, dZ, KS, MS, US, N, H, R, L0

def eta(L): return N * math.pi * R * R / (H * L)

def mode_period(M, k, L, T):
    e = eta(L); cs = math.sqrt(T * (Z(e) + e * dZ(e) + Z(e) ** 2)); kS = N * cs * cs / L ** 2
    f = lambda w: M * w * w - k - kS * (w * L / cs) / math.tan(w * L / cs)
    w0 = math.sqrt((k + kS) / (M + N / 3.0)); hi = min(1.6 * w0, 0.999 * math.pi * cs / L)
    return 2 * math.pi / brentq(f, 0.5 * w0, hi)

def gate(ds, P, thr):
    n = min(len(d) for d in ds); t = ds[0]["Time"].to_numpy(float)[:n]
    xb = np.mean([d["W0_x_sigma"].to_numpy(float)[:n] for d in ds], axis=0); te = t[-1]
    last = xb[(t > te - P) & (t <= te)].mean(); prev = xb[(t > te - 2 * P) & (t <= te - P)].mean()
    dlt = abs(last - prev); return dlt, ("PASS" if dlt < thr else "FAIL")

def main():
    Ph, rev = reversible(); TR = {(k, M): tau_r(k, M, rev)["tr"] for k in KS for M in MS}
    G, E = {}, {}
    print("### Position-settled gate, all 36 cells (Delta in sigma; threshold 0.01 s_rev)\n")
    print("| k | M_s | u | P_m | Delta | 0.01 s_rev | label |")
    print("|---|---|---|---|---|---|---|")
    for k in KS:
        thr = 0.01 * rev[k]["s"]
        for M in MS:
            P = mode_period(M, k, rev[k]["Lf"], rev[k]["Tf"]); cs = A.load(f"ctrl_k{k}_M{M}")
            TW = 2 * math.pi / math.sqrt((k + A.KGAS) / M)
            for u in US:
                ds = A.load(f"k{k}_M{M}_u{u}"); dl, lab = gate(ds, P, thr); G[(k, M, u)] = (dl, lab, P)
                E[(k, M, u)] = dict(em=A.eps_window(ds, cs, k, 3 * TW), es=A.eps_window(ds, cs, k, TR[(k, M)]))
                print(f"| {k} | {M} | {u} | {P:.2f} | {dl:.5f} | {thr:.5f} | {lab} |")
    print("\n**PASS count per (k, M_s):** " + "; ".join(f"k = {k}, M_s = {M}: {sum(G[(k, M, u)][1] == 'PASS' for u in US)}/6"
                                                     for k in KS for M in MS)
          + f". **Total {sum(v[1] == 'PASS' for v in G.values())}/36.**")
    print("\n### Calibration: the same gate on the no-push controls (settled by construction)\n")
    print("| k | M_s | P_m (eta_0, T = 1) | Delta | 0.01 s_rev | label |")
    print("|---|---|---|---|---|---|")
    C = {}
    for k in KS:
        for M in MS:
            P = mode_period(M, k, L0, 1.0); dl, lab = gate(A.load(f"ctrl_k{k}_M{M}"), P, 0.01 * rev[k]["s"]); C[(k, M)] = (dl, lab)
            print(f"| {k} | {M} | {P:.2f} | {dl:.5f} | {0.01*rev[k]['s']:.5f} | {lab} |")
    print(f"\ncontrols PASS: {sum(v[1] == 'PASS' for v in C.values())}/6")
    print("\n### DIAGNOSTIC, post hoc, not a gate: the controls' Delta against window length (noise floor of the rule)\n")
    print("| k | M_s | 0.01 s_rev | Delta, 1 P_m | 10 P_m | 30 P_m | tau_r |")
    print("|---|---|---|---|---|---|---|")
    for k in KS:
        for M in MS:
            P = mode_period(M, k, L0, 1.0); cs = A.load(f"ctrl_k{k}_M{M}"); thr = 0.01 * rev[k]["s"]
            row = [gate(cs, w, thr)[0] for w in (P, 10 * P, 30 * P, TR[(k, M)])]
            print(f"| {k} | {M} | {thr:.5f} | " + " | ".join(f"{x:.5f}" for x in row) + " |")
    figures(G, E, rev, sum(v[1] == 'PASS' for v in C.values()))
    json.dump(dict(cells={f"{k}_{M}_{u}": dict(delta=v[0], label=v[1], P_m=v[2]) for (k, M, u), v in G.items()},
                   controls={f"{k}_{M}": dict(delta=v[0], label=v[1]) for (k, M), v in C.items()}),
              open(os.path.join(HERE, "261013_effmap_gate.json"), "w"), indent=1)

def figures(G, E, rev, nctrl):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    def panel(ax, k, which, faint=False):
        for M, off, mk in ((50, 0.94, "o"), (200, 1.06, "s")):
            xs = [u * off for u in US]; ys = [E[(k, M, u)][which]["eps"] for u in US]; es = [E[(k, M, u)][which]["err"] for u in US]
            ax.errorbar(xs, ys, yerr=es, fmt=mk, color="tab:blue", mfc="white" if faint else "tab:blue", alpha=0.45 if faint else 1.0,
                        ms=4 if faint else 5, capsize=2 if faint else 3,
                        label=f"{'ε_mean (last 3 T_w)' if which == 'em' else 'ε_settled (last τ_r)'}, M_s = {M}")
            if not faint:
                bad = [(x, y) for x, y, u in zip(xs, ys, US) if G[(k, M, u)][1] == "FAIL"]
                if bad: ax.plot(*zip(*bad), "x", color="black", ms=10, mew=2, label="position-settled: FAIL" if M == 50 else None)
        n = sum(G[(k, M, u)][1] == "PASS" for M in MS for u in US)
        if faint: return
        ax.axhline(rev[k]["eps"], color="red", lw=1.5, label=f"ε_rev (KR) = {rev[k]['eps']:.4f}")
        ax.axvline(0.0884, color="0.6", ls=":", lw=1, label="d/u = 2L₀/c_s")
        ax.set_xscale("log"); ax.set_xlabel("piston speed u")
        ax.set_title(f"k = {k}   (position-settled: {n}/12 PASS)", fontsize=10)
    for name, layers in (("261012_p2_effmap", [("em", False)]), ("261012_p2_effmap_both", [("es", False), ("em", True)])):
        fig, axs = plt.subplots(1, 3, figsize=(12, 3.9), sharex=True)
        for ax, k in zip(axs, KS):
            for which, faint in layers: panel(ax, k, which, faint)
        axs[0].set_ylabel("ε = ΔE_spring / W_in"); axs[0].legend(fontsize=6.5, frameon=False, loc="lower left")
        fig.suptitle("Efficiency map, geometry C. Settle label: position-settled gate (261010 § 3.2), NOISE-LIMITED: "
                     f"settled no-push controls pass only {nctrl}/6, so × marks noise, not motion. Points offset in u.", fontsize=9)
        fig.tight_layout()
        for ext in ("png", "pdf"): fig.savefig(os.path.join(A.FIG, f"{name}.{ext}"), dpi=200)
    for k in KS:
        fig, ax = plt.subplots(figsize=(4.6, 3.6)); panel(ax, k, "es"); ax.set_ylabel("ε = ΔE_spring / W_in")
        ax.text(0.98, 0.98, f"gate noise-limited: controls {nctrl}/6 PASS", transform=ax.transAxes, ha="right", va="top", fontsize=6.5)
        ax.legend(fontsize=6.5, frameon=False, loc="lower left"); fig.tight_layout()
        for ext in ("png", "pdf"): fig.savefig(os.path.join(A.FIG, f"261012_p2_effmap_k{k}.{ext}"), dpi=200)
    print("\nfigures: 261012_p2_effmap, 261012_p2_effmap_both, 261012_p2_effmap_k{0.25,0.5,1.0} (.png/.pdf) in "
          "0000_PLAN_OVERALL/paper2_energytransfer/experiments/final/")

if __name__ == "__main__":
    main()
