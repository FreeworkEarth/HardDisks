#!/usr/bin/env python3
"""##CHRIS 2026-10-14: settle gate G2 (261010 sec. 3.4, pre-registered 5190846). No new runs.

x_bar = 8-seed mean divider trajectory (W0_x_sigma). W1 = last tau_r of the record, W0 = the tau_r before it,
tau_r = 2 M_hat/gamma_1 per (k, M_s) (sec. 1.8 table 1). Signed drift d = <x_bar>_W1 - <x_bar>_W0.
sigma_ctrl = std (ddof = 1) of the six signed control drifts; theta(k) = max(0.01 s_rev(k), 3 sigma_ctrl).
VALID only if 6/6 controls have |d| < theta(k); then PASS/FAIL per cell. Re-renders 261012_p2_effmap(_both) and
the per-k panels 261012_p2_effmap_k<k> with the G2 label.
"""
import json, math, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import paper2_effmap_analysis_20261012 as A
from paper2_effmap_amend_20261012 import reversible, tau_r, KS, MS, US

def drift(ds, tr):
    n = min(len(d) for d in ds); t = ds[0]["Time"].to_numpy(float)[:n]
    xb = np.mean([d["W0_x_sigma"].to_numpy(float)[:n] for d in ds], axis=0); te = t[-1]
    w1 = (t > te - tr) & (t <= te); w0 = (t > te - 2 * tr) & (t <= te - tr)
    return float(xb[w1].mean() - xb[w0].mean())

def main():
    Ph, rev = reversible(); TR = {(k, M): tau_r(k, M, rev)["tr"] for k in KS for M in MS}
    ctrl = {(k, M): A.load(f"ctrl_k{k}_M{M}") for k in KS for M in MS}
    dc = {(k, M): drift(ctrl[(k, M)], TR[(k, M)]) for k in KS for M in MS}
    sc = float(np.std(list(dc.values()), ddof=1)); th = {k: max(0.01 * rev[k]["s"], 3 * sc) for k in KS}
    print("### Gate G2 -- controls (signed drift d over the last tau_r vs the tau_r before)\n")
    print("| k | M_s | tau_r | d (control) | theta(k) | label |\n|---|---|---|---|---|---|")
    for (k, M), d in dc.items():
        print(f"| {k} | {M} | {TR[(k, M)]:.0f} | {d:+.5f} | {th[k]:.5f} | {'PASS' if abs(d) < th[k] else 'FAIL'} |")
    nc = sum(abs(d) < th[k] for (k, M), d in dc.items()); valid = nc == 6
    print(f"\nsigma_ctrl = {sc:.5f} sigma (std, ddof = 1, of the six signed control drifts); 3 sigma_ctrl = {3*sc:.5f}")
    print("theta(k) = max(0.01 s_rev, 3 sigma_ctrl): " + "; ".join(f"k = {k}: max({0.01*rev[k]['s']:.5f}, {3*sc:.5f}) = {th[k]:.5f}" for k in KS))
    print(f"controls passing: {nc}/6 -> **G2 {'VALID' if valid else 'NOT VALID -- no cell is labelled'}**\n")
    G, E = {}, {}
    print("### Gate G2 -- the 36 cells\n")
    print("| k | M_s | u | d | theta(k) | label |\n|---|---|---|---|---|---|")
    for k in KS:
        for M in MS:
            cs = ctrl[(k, M)]; TW = 2 * math.pi / math.sqrt((k + A.KGAS) / M)
            for u in US:
                ds = A.load(f"k{k}_M{M}_u{u}"); d = drift(ds, TR[(k, M)])
                lab = ("PASS" if abs(d) < th[k] else "FAIL") if valid else "n/a"; G[(k, M, u)] = (d, lab)
                E[(k, M, u)] = dict(em=A.eps_window(ds, cs, k, 3 * TW), es=A.eps_window(ds, cs, k, TR[(k, M)]))
                print(f"| {k} | {M} | {u} | {d:+.5f} | {th[k]:.5f} | {lab} |")
    print("\n**PASS count per (k, M_s):** " + "; ".join(f"k = {k}, M_s = {M}: {sum(G[(k, M, u)][1] == 'PASS' for u in US)}/6" for k in KS for M in MS)
          + f". **Total {sum(v[1] == 'PASS' for v in G.values())}/36.**")
    figures(G, E, rev, valid)
    json.dump(dict(sigma_ctrl=sc, theta={str(k): v for k, v in th.items()}, valid=valid,
                   controls={f"{k}_{M}": d for (k, M), d in dc.items()},
                   cells={f"{k}_{M}_{u}": dict(d=v[0], label=v[1]) for (k, M, u), v in G.items()}),
              open(os.path.join(HERE, "261014_effmap_gate2.json"), "w"), indent=1)

def figures(G, E, rev, valid):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    def panel(ax, k, which, faint=False):
        for M, off, mk in ((50, 0.94, "o"), (200, 1.06, "s")):
            xs = [u * off for u in US]; ys = [E[(k, M, u)][which]["eps"] for u in US]; es = [E[(k, M, u)][which]["err"] for u in US]
            ax.errorbar(xs, ys, yerr=es, fmt=mk, color="tab:blue", mfc="white" if faint else "tab:blue", alpha=0.45 if faint else 1.0,
                        ms=4 if faint else 5, capsize=2 if faint else 3,
                        label=f"{'ε_mean (last 3 T_w)' if which == 'em' else 'ε_settled (last τ_r)'}, M_s = {M}")
            if not faint:
                bad = [(x, y) for x, y, u in zip(xs, ys, US) if G[(k, M, u)][1] == "FAIL"]
                if bad: ax.plot(*zip(*bad), "x", color="black", ms=10, mew=2, label="G2: FAIL" if M == 50 else None)
        if faint: return
        n = sum(G[(k, M, u)][1] == "PASS" for M in MS for u in US)
        ax.axhline(rev[k]["eps"], color="red", lw=1.5, label=f"ε_rev (KR) = {rev[k]['eps']:.4f}")
        ax.axvline(0.0884, color="0.6", ls=":", lw=1, label="d/u = 2L₀/c_s")
        ax.set_xscale("log"); ax.set_xlabel("piston speed u")
        ax.set_title(f"k = {k}   (settled, G2: {n}/12 PASS)" if valid else f"k = {k}   (G2 not valid: no label)", fontsize=10)
    for name, layers in (("261012_p2_effmap", [("em", False)]), ("261012_p2_effmap_both", [("es", False), ("em", True)])):
        fig, axs = plt.subplots(1, 3, figsize=(12, 3.9), sharex=True)
        for ax, k in zip(axs, KS):
            for which, faint in layers: panel(ax, k, which, faint)
        axs[0].set_ylabel("ε = ΔE_spring / W_in"); axs[0].legend(fontsize=6.5, frameon=False, loc="lower left")
        fig.suptitle("Efficiency map, geometry C. Settle label: gate G2 (261010 § 3.4; " + ("valid: 6/6 no-push controls pass)" if valid else "NOT valid)")
                     + "; points offset in u", fontsize=9.5)
        fig.tight_layout()
        for ext in ("png", "pdf"): fig.savefig(os.path.join(A.FIG, f"{name}.{ext}"), dpi=200)
    for k in KS:
        fig, ax = plt.subplots(figsize=(4.6, 3.6)); panel(ax, k, "es"); ax.set_ylabel("ε = ΔE_spring / W_in")
        ax.legend(fontsize=6.5, frameon=False, loc="lower left"); fig.tight_layout()
        for ext in ("png", "pdf"): fig.savefig(os.path.join(A.FIG, f"261012_p2_effmap_k{k}.{ext}"), dpi=200)
    print("\nfigures: 261012_p2_effmap, 261012_p2_effmap_both, 261012_p2_effmap_k{0.25,0.5,1.0} (.png/.pdf), G2 label")

if __name__ == "__main__":
    main()
