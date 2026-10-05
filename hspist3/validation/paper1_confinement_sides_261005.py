#!/usr/bin/env python3
"""##CHRIS 2026-10-05 (Task BB): which side carries the 1/N_s residual of the identity -- DIAGNOSTIC, no verdict (261012 sec. 3.10).

Per cell, both densities: the dynamic stiffness k_S^dyn (method B, C1) and the static stiffness k_static = k_T + F^2/(N_s kT)
(A-fixed, sec. 3.4), each divided by the bulk KR stiffness
    k_KR = N_s m c_s^KR(eta_true)^2 / L_eff,true^2                     (no wall correction)
and by the same corrected for the measured 1/H wall shift, k_KR,B = k_KR (1 + b/H)^2 with the registered B amplitude b
(sec. 2.2 fit: 0.11797 at eta 0.10, 0.26256 at pi/8; at pi/8 the B form is excluded, so the corrected column there is
indicative only). Each ratio minus one is fitted per density with a free intercept and a 1/N_s slope (weighted, sigma of
the ratio from sigma(k_S^dyn) resp. sigma(k_T)). The side whose slope differs from zero carries the 1/N_s trend.
Note [DERIVATION]: k_S^dyn/k_KR - 1 ~ 2 Delta (the c_s shift of sec. 2.2, heavy masses only), and k_static/k_KR depends on the
same L_eff,true convention; only the DIFFERENCE of the two slopes is convention-free (it is the identity residual's slope).
Figure: paper1_speedofsound/experiments/final/261005_p1_identity_sides.{png,pdf}.
usage (from hspist3/): python3 validation/paper1_confinement_sides_261005.py
"""
import contextlib, io, math, os, sys
import numpy as np
from scipy.stats import chi2 as CHI2
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, HS)
import paper1_confinement_results_261004 as R
import paper1_confinement_afix_261005 as AF

B_AMP = {"0.10": 0.11797, "0.39": 0.26256}          # registered one-amplitude B fits, sec. 2.2 (printed there)


def fit2(y, s, Ns):
    X = np.vstack([np.ones(len(Ns)), 1 / np.array(Ns, float)]).T; w = 1 / np.array(s) ** 2
    W = np.diag(w); cov = np.linalg.inv(X.T @ W @ X); b = cov @ X.T @ W @ np.array(y)
    ch = float((w * (np.array(y) - X @ b) ** 2).sum()); return b, np.sqrt(np.diag(cov)), ch, len(Ns) - 2


def main():
    with contextlib.redirect_stdout(io.StringIO()):
        CS = R.cells()
        for c in CS:
            R.method_B(c); R.method_A(c); R.identity(c); c["af"] = AF.inventory_and_static(c)
    print("## Which side bends (261012 sec. 3.10) -- DIAGNOSTIC, no verdict\n")
    print("| eta | cell | N_s | H | k_KR (bulk) | k_S^dyn/k_KR - 1 [%] | +- | k_static/k_KR - 1 [%] | +- | 1/H factor (1 + b/H)^2 - 1 [%] "
          "| k_S^dyn/k_KR,B - 1 [%] | k_static/k_KR,B - 1 [%] |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for c in CS:
        kKR = c["Ns"] * c["KR"] ** 2 / c["LeT"] ** 2; fB = (1 + B_AMP[c["lab"]] / c["H"]) ** 2
        c.update(kKR=kKR, rd=c["kS"] / kKR - 1, srd=c["s_kS"] / kKR, rs=c["af"]["static"] / kKR - 1, srs=c["af"]["s_kT"] / kKR,
                 rdB=c["kS"] / (kKR * fB) - 1, rsB=c["af"]["static"] / (kKR * fB) - 1, fB=fB)
        print(f"| {c['lab']} | {c['cid']} | {c['Ns']} | {c['H']:g} | {kKR:.6g} | {100 * c['rd']:+.3f} | {100 * c['srd']:.3f} | "
              f"{100 * c['rs']:+.3f} | {100 * c['srs']:.3f} | {100 * (fB - 1):+.3f} | {100 * c['rdB']:+.3f} | {100 * c['rsB']:+.3f} |")
    print("\n### Fits per density: ratio - 1 = a + s/N_s (weighted; free intercept)\n")
    print("| eta | reference | side | intercept a [%] | slope s | +- | s/sigma | chi2/dof |\n|---|---|---|---|---|---|---|---|")
    res = {}
    for lab in ("0.10", "0.39"):
        cc = [c for c in CS if c["lab"] == lab]; Ns = [c["Ns"] for c in cc]
        for ref, kd, ks in (("bulk KR", "rd", "rs"), ("KR x (1 + b/H)^2", "rdB", "rsB")):
            for side, key, sk in (("dynamic", kd, "srd"), ("static", ks, "srs")):
                b, sb, ch, dof = fit2([c[key] for c in cc], [c[sk] for c in cc], Ns)
                res[(lab, ref, side)] = (b, sb)
                print(f"| {lab} | {ref} | {side} | {100 * b[0]:+.3f} +- {100 * sb[0]:.3f} | {b[1]:+.4f} | {sb[1]:.4f} | {b[1] / sb[1]:+.1f} | {ch:.1f}/{dof} |")
    print("\n(slope s in units of 1/N_s: a slope of 0.48 means +0.48/N_s, i.e. +1.9 % at N_s = 25)")
    figure(CS)
    print("\nfigure -> 261005_p1_identity_sides.png/.pdf")


def figure(CS):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 2, figsize=(13, 9), sharex=True)
    for col, lab in enumerate(("0.10", "0.39")):
        cc = [c for c in CS if c["lab"] == lab]; x = np.array([1 / c["Ns"] for c in cc])
        for row, (kd, ks, ttl) in enumerate((("rd", "rs", "against the bulk KR stiffness"),
                                              ("rdB", "rsB", "against KR x (1 + b/H)^2 (measured wall shift removed)"))):
            ax = axs[row, col]
            ax.axhline(0, color=R.RED, lw=2, label="KR 2006 (bulk)" if row == 0 else "KR x (1 + b/H)^2")
            ax.errorbar(x * 0.97, [100 * c[kd] for c in cc], yerr=[100 * c["srd"] for c in cc], fmt="o", color=R.BLUE2, ms=6, capsize=3,
                        label="dynamic  k_S^dyn (method B)")
            ax.errorbar(x * 1.03, [100 * c[ks] for c in cc], yerr=[100 * c["srs"] for c in cc], fmt="s", mfc="white", color=R.BLUE,
                        ms=6, capsize=3, label="static  k_T + F^2/(N_s kT) (A-fixed)")
            ax.set_title(f"eta = {'0.100' if lab == '0.10' else 'pi/8'}: {ttl}", fontsize=9.5); ax.grid(True, ls=":", alpha=0.5)
            if col == 0: ax.set_ylabel("stiffness / reference - 1  [%]")
            if row == 1: ax.set_xlabel("1 / N_s")
            ax.legend(fontsize=7.5)
    fig.suptitle("Which side carries the 1/N_s residual? (diagnostic, 261012 sec. 3.10; L_eff,true = L_0 - 2r - t/2)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(R.OUT, f"261005_p1_identity_sides.{ext}"), dpi=200)


if __name__ == "__main__":
    main()
