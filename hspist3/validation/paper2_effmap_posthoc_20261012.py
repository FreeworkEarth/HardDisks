#!/usr/bin/env python3
"""##CHRIS 2026-10-12: efficiency map -- POST-HOC summaries, written AFTER the pre-registered analysis ran.

Nothing here changes a pre-registered number. paper2_effmap_analysis_20261012.py (committed before launch)
fitted the prediction-2 slope on epsilon_mean, whose last-3-T_w window leaves the divider's thermal motion
unaveraged (errors 5-20x those of epsilon_settled). This script prints, labelled post hoc:
  (a) epsilon_settled averaged over the quasi-static cells u <= 0.05, per (k, M_s) and pooled, vs epsilon_rev;
  (b) epsilon_settled slope on u in {0.1, 0.2, 0.5};
  (c) the M_s = 50 vs 200 difference per (k, u);
  (d) W_in and s_corr at u <= 0.05 against the reversible reference (where the ratio's excess comes from);
  (e) the no-push controls: Delta KE_gas against -E_dof, per cell and chi^2;
  (f) a figure with both estimators.
"""
import math, os, sys
import numpy as np
from scipy.stats import chi2 as chi2d
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import paper2_effmap_analysis_20261012 as A
from paper2_effmap_amend_20261012 import reversible, tau_r, KS, MS, US, N, L0

def wmean(v, s):
    w = 1 / np.asarray(s) ** 2; m = float(np.sum(w * v) / np.sum(w)); e = 1 / math.sqrt(np.sum(w))
    ch = float(np.sum(w * (np.asarray(v) - m) ** 2)); return m, e, ch, len(v) - 1

def main():
    Ph, rev = reversible(); TR = {(k, M): tau_r(k, M, rev)["tr"] for k in KS for M in MS}
    ctrl = {(k, M): A.load(f"ctrl_k{k}_M{M}") for k in KS for M in MS}
    R = {}
    for k in KS:
        for M in MS:
            for u in US:
                ds = A.load(f"k{k}_M{M}_u{u}"); TW = 2 * math.pi / math.sqrt((k + A.KGAS) / M)
                R[(k, M, u)] = dict(em=A.eps_window(ds, ctrl[(k, M)], k, 3 * TW), es=A.eps_window(ds, ctrl[(k, M)], k, TR[(k, M)]))
    print("### (a) POST HOC: epsilon_settled over the quasi-static cells u <= 0.05 (inverse-variance mean)\n")
    print("| k | M_s | mean eps_settled | chi2/dof across u | eps_rev | ratio to eps_rev | (mean - eps_rev)/sigma |")
    print("|---|---|---|---|---|---|---|")
    for k in KS:
        allv, alls = [], []
        for M in MS + ("both",):
            if M == "both": v, s = allv, alls
            else:
                v = [R[(k, M, u)]["es"]["eps"] for u in US if u <= 0.05]; s = [R[(k, M, u)]["es"]["err"] for u in US if u <= 0.05]
                allv += v; alls += s
            m, e, ch, dof = wmean(v, s); er = rev[k]["eps"]
            print(f"| {k} | {M} | {m:.5f} ± {e:.5f} | {ch:.1f}/{dof} | {er:.4f} | {m/er:.4f} | {(m-er)/e:+.1f} |")
    print("\n### (b) POST HOC: epsilon_settled slope on u in {0.1, 0.2, 0.5} (weighted linear fit)\n")
    print("| k | M_s | slope ± σ | intercept | chi2 (1 dof) | pre-registered eps_mean slope (for reference) |")
    print("|---|---|---|---|---|---|")
    for k in KS:
        for M in MS:
            x = np.array([0.1, 0.2, 0.5]); y = np.array([R[(k, M, u)]["es"]["eps"] for u in x]); s = np.array([R[(k, M, u)]["es"]["err"] for u in x])
            w = 1 / s ** 2; X = np.vstack([x, np.ones_like(x)]).T; C = np.linalg.inv(X.T @ (X * w[:, None])); p = C @ (X.T @ (w * y))
            ch = float(np.sum(w * (y - X @ p) ** 2))
            ym = np.array([R[(k, M, u)]["em"]["eps"] for u in x]); sm = np.array([R[(k, M, u)]["em"]["err"] for u in x])
            wm = 1 / sm ** 2; Cm = np.linalg.inv(X.T @ (X * wm[:, None])); pm = Cm @ (X.T @ (wm * ym))
            print(f"| {k} | {M} | {p[0]:+.4f} ± {math.sqrt(C[0,0]):.4f} | {p[1]:.4f} | {ch:.1f} | {pm[0]:+.4f} ± {math.sqrt(Cm[0,0]):.4f} |")
    print("\n### (c) POST HOC: M_s = 50 vs 200, epsilon_settled, per (k, u)\n")
    zs = []
    print("| k | u | eps(50) | eps(200) | difference / sigma |")
    print("|---|---|---|---|---|")
    for k in KS:
        for u in US:
            a, b = R[(k, 50, u)]["es"], R[(k, 200, u)]["es"]; z = (a["eps"] - b["eps"]) / math.hypot(a["err"], b["err"]); zs.append(z)
            print(f"| {k} | {u} | {a['eps']:.4f} | {b['eps']:.4f} | {z:+.1f} |")
    c2 = float(np.sum(np.square(zs))); print(f"\nchi2 = {c2:.1f} on {len(zs)} dof, p = {chi2d.sf(c2, len(zs)):.3g}")
    print("\n### (d) POST HOC: where the quasi-static excess sits -- W_in and s_corr (last tau_r) at u <= 0.05\n")
    print("| k | M_s | <W_in>/W_rev | s_corr/s_rev | Delta E_spring/E_spring,rev |")
    print("|---|---|---|---|---|")
    for k in KS:
        for M in MS:
            c = [R[(k, M, u)]["es"] for u in US if u <= 0.05]
            W = np.mean([q["W"] for q in c]); s = np.mean([q["s"] for q in c]); dE = np.mean([q["dEs"] for q in c])
            print(f"| {k} | {M} | {W/rev[k]['W']:.4f} | {s/rev[k]['s']:.4f} | {dE/rev[k]['Esp']:.4f} |")
    print("\n### (e) POST HOC: no-push controls, Delta KE_gas (last half) against -E_dof\n")
    print("| k | M_s | Delta KE_gas | -E_dof | z |")
    print("|---|---|---|---|---|")
    zc = []
    for k in KS:
        for M in MS:
            cs = ctrl[(k, M)]
            h = np.array([d["KE_gas_total"].to_numpy(float)[len(d)//2:].mean() - d["KE_gas_total"].iloc[0] for d in cs])
            Tc = np.mean([d["KE_gas_total"].to_numpy(float)[len(d)//2:].mean() / N for d in cs]); pr = -A.e_dof(k, L0, Tc)
            z = (h.mean() - pr) / (h.std(ddof=1) / math.sqrt(len(h))); zc.append(z)
            print(f"| {k} | {M} | {h.mean():+.4f} ± {h.std(ddof=1)/math.sqrt(len(h)):.4f} | {pr:+.4f} | {z:+.1f} |")
    c2 = float(np.sum(np.square(zc))); print(f"\nchi2 = {c2:.1f} on {len(zc)} dof, p = {chi2d.sf(c2, len(zc)):.3g}")
    figure(R, rev)

def figure(R, rev):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.9), sharex=True)
    for ax, k in zip(axs, KS):
        for M, off, mk in ((50, 0.94, "o"), (200, 1.06, "s")):
            ax.errorbar([u * off for u in US], [R[(k, M, u)]["es"]["eps"] for u in US], yerr=[R[(k, M, u)]["es"]["err"] for u in US],
                        fmt=mk, color="tab:blue", ms=5, capsize=3, label=f"ε_settled (last τ_r), M_s = {M}")
            ax.errorbar([u * off * 1.015 for u in US], [R[(k, M, u)]["em"]["eps"] for u in US], yerr=[R[(k, M, u)]["em"]["err"] for u in US],
                        fmt=mk, color="tab:blue", mfc="white", alpha=0.45, ms=4, capsize=2, label=f"ε_mean (last 3 T_w), M_s = {M}")
        ax.axhline(rev[k]["eps"], color="red", lw=1.5, label=f"ε_rev (KR) = {rev[k]['eps']:.4f}")
        ax.axvline(0.0884, color="0.6", ls=":", lw=1, label="d/u = 2L₀/c_s")
        ax.set_xscale("log"); ax.set_title(f"k = {k}", fontsize=10); ax.set_xlabel("piston speed u")
    axs[0].set_ylabel("ε = ΔE_spring / W_in"); axs[0].legend(fontsize=6.5, frameon=False, loc="lower left")
    fig.suptitle("Efficiency map, geometry C: both pre-registered estimators (post-hoc figure; points offset in u for legibility)", fontsize=10)
    fig.tight_layout()
    for ext in ("png", "pdf"): fig.savefig(os.path.join(A.FIG, f"261012_p2_effmap_both.{ext}"), dpi=200)
    print("\nfigure: 0000_PLAN_OVERALL/paper2_energytransfer/experiments/final/261012_p2_effmap_both.{png,pdf}")

if __name__ == "__main__":
    main()
