#!/usr/bin/env python3
"""##CHRIS 2026-10-04 (Task Z): EXPLORATORY analyses of the confinement campaign (261012 sec. 2.9). No verdicts.
Inputs: the registered per-cell results (paper1_confinement_results_261004: method B Delta, sigma, L_eff,true, k_S^dyn) and
the post-hoc drift-corrected identity residual (paper1_confinement_heldwall_posthoc_261004: rho_I with C4).

Z1 effective acoustic length.
  (a) per density, over the H = 10 cells (anchor + L-scan) and the aspect cells: Delta = -eps/L_eff + b/H + c0, weighted
      least squares (sigma = the registered scaled sigma of Delta); eps > 0 means the acoustic length of the frequency
      formula is LONGER than L_eff = L0 - 2r - t/2 [DERIVATION: nu fixed, c_inferred = 2 pi nu L_eff/K, so a true length
      L_eff + eps gives Delta = Delta_true - eps/L_eff to first order].
  (b) per cell, L_KR = c_s^KR(eta) sqrt(N_s m/k_S^dyn), the length that makes the length-free dynamic stiffness equal the
      bulk KR stiffness N_s m c_KR^2/L^2; its offset from L_eff,true, compared with the (a) model's offset
      eps - (b/H + c0) L_eff (the same Delta model written as a length).
Z2 acoustic height deficit: Delta = (dln c_s/dln eta)(2 delta_H/H) => delta_H = b/(2 s), s = dln c_KR/dln eta (central
  difference of tests_20260913.kr_cs), b = the B amplitude (Delta = b/H). Three b's: the registered one-amplitude B fit over
  all cells (sec. 2.2), the H-scan alone, and Z1(a)'s b. Errors raw and, where chi2/dof > 1, scaled by sqrt(chi2/dof).
Z3 the tension: per density, rho_I(C4) = c/N_s (one parameter, weighted); under hypothesis C the mode shift is
  Delta_C = rho_I/2 = c/(2 N_s), so the L-scan would separate N_s = 25 from N_s = 100 by (c/2)(1/25 - 1/100); compared with
  the measured Delta(N_s = 25) - Delta(N_s = 100) of the L-scan.
usage (from hspist3/): python3 validation/paper1_confinement_exploratory_261004.py
"""
import contextlib, io, math, os, sys
import numpy as np
from scipy.stats import chi2 as CHI2
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, HS)
import tests_20260913 as T
import paper1_confinement_results_261004 as R
import paper1_confinement_heldwall_posthoc_261004 as PH


def wls(X, y, s):
    W = np.diag(1 / s ** 2); cov = np.linalg.inv(X.T @ W @ X); b = cov @ X.T @ W @ y
    ch = float(((y - X @ b) ** 2 / s ** 2).sum()); dof = len(y) - X.shape[1]
    return b, np.sqrt(np.diag(cov)), ch, dof


def slope_s(eta, h=1e-4):
    return (math.log(T.kr_cs(eta * (1 + h))) - math.log(T.kr_cs(eta * (1 - h)))) / (2 * h)


def main():
    with contextlib.redirect_stdout(io.StringIO()):
        CS = R.cells()
        for c in CS:
            R.method_B(c); R.method_A(c); R.identity(c)
        for c in CS:
            d = PH.drift(c); st = d["kTc"] + c["F2term"]
            c["rho_c"] = (c["kS"] - st) / c["kS"]; c["s_rho_c"] = math.sqrt((c["s_kS"] / c["kS"]) ** 2 + (d["s_kTc"] / c["kS"]) ** 2)
        res = R.confinement_verdict(CS)
    print("## EXPLORATORY (261012 sec. 2.9) -- no verdicts\n")

    print("### Z1 (a): Delta = -eps/L_eff + b/H + c0 over the H = 10 cells and the aspect cells\n")
    print("| eta | cells | eps [sigma] | b [sigma] | c0 [%] | chi2 / dof | p | errors scaled by sqrt(chi2/dof) (eps, b, c0) |")
    print("|---|---|---|---|---|---|---|---|")
    Z1 = {}
    for lab in ("0.10", "0.39"):
        cc = [c for c in CS if c["lab"] == lab and ((c["scan"] in ("H", "L") and abs(c["H"] - 10) < 1e-9) or c["scan"] == "aspect")]
        X = np.array([[-1 / c["LeT"], 1 / c["H"], 1.0] for c in cc]); y = np.array([c["D"] for c in cc]); s = np.array([c["sD"] for c in cc])
        b, sb, ch, dof = wls(X, y, s); k = math.sqrt(max(1.0, ch / dof)); Z1[lab] = (b, sb, ch, dof, cc)
        print(f"| {lab} | {len(cc)} ({', '.join(c['cid'].split('_', 1)[1] for c in cc)}) | {b[0]:+.4f} +- {sb[0]:.4f} | {b[1]:+.4f} +- {sb[1]:.4f} | "
              f"{100 * b[2]:+.3f} +- {100 * sb[2]:.3f} | {ch:.2f} / {dof} | {CHI2.sf(ch, dof):.3g} | "
              f"{sb[0] * k:.4f}, {sb[1] * k:.4f}, {100 * sb[2] * k:.3f} % |")

    print("\n### Z1 (b): per cell, the length L_KR that makes k_S^dyn equal the bulk KR stiffness N_s m c_KR^2/L^2\n")
    print("| eta | cell | L_eff,true | L_KR = c_KR sqrt(N_s m/k_S^dyn) | offset L_KR - L_eff,true [sigma] | +- | "
          "Z1(a) model offset eps - (b/H + c0) L_eff | (offset - model)/sigma |")
    print("|---|---|---|---|---|---|---|---|")
    for lab in ("0.10", "0.39"):
        b = Z1[lab][0]
        for c in [c for c in CS if c["lab"] == lab]:
            Lkr = c["KR"] * math.sqrt(c["Ns"] / c["kS"]); off = Lkr - c["LeT"]; s_off = Lkr * 0.5 * c["s_kS"] / c["kS"]
            mod = b[0] - (b[1] / c["H"] + b[2]) * c["LeT"]
            print(f"| {lab} | {c['cid']} | {c['LeT']:.4f} | {Lkr:.4f} | {off:+.4f} | {s_off:.4f} | {mod:+.4f} | {(off - mod) / s_off:+.1f} |")

    print("\n### Z2: acoustic height deficit delta_H = b / (2 dln c_s/dln eta)\n")
    print("| eta | s = dln c_KR/dln eta | b source | b [sigma] | +- raw | chi2/dof of that fit | delta_H [sigma] | +- raw | +- scaled |")
    print("|---|---|---|---|---|---|---|---|---|")
    DH = {}
    for lab in ("0.10", "0.39"):
        anc = [c for c in CS if c["lab"] == lab and c["scan"] == "H" and abs(c["H"] - 10) < 1e-9][0]; s = slope_s(anc["eta"])
        fB = res[lab]["fits"]["B"]
        hs = [c for c in CS if c["lab"] == lab and c["scan"] == "H"]
        bh, sbh, chh, dofh = wls(np.array([[1 / c["H"]] for c in hs]), np.array([c["D"] for c in hs]), np.array([c["sD"] for c in hs]))
        z = Z1[lab]
        for name, bb, sbb, ch, dof in (("registered B fit, all cells (sec. 2.2)", fB["a"], fB["sa"], fB["chi2"], fB["dof"]),
                                       ("H-scan only, b/H", bh[0], sbh[0], chh, dofh),
                                       ("Z1(a) b (with eps, c0)", z[0][1], z[1][1], z[2], z[3])):
            k = math.sqrt(max(1.0, ch / dof)); dh = bb / (2 * s); sdh = sbb / (2 * s)
            DH.setdefault(name, {})[lab] = (dh, sdh, sdh * k)
            print(f"| {lab} | {s:.4f} | {name} | {bb:+.5f} | {sbb:.5f} | {ch:.2f}/{dof} | {dh:+.4f} | {sdh:.4f} | {sdh * k:.4f} |")
    print("\none delta_H for both densities? (difference / combined sigma)\n")
    print("| b source | delta_H(0.10) - delta_H(pi/8) [sigma] | / raw sigma | / scaled sigma |\n|---|---|---|---|")
    for name, v in DH.items():
        a, b = v["0.10"], v["0.39"]
        print(f"| {name} | {a[0] - b[0]:+.4f} | {(a[0] - b[0]) / math.hypot(a[1], b[1]):+.1f} | {(a[0] - b[0]) / math.hypot(a[2], b[2]):+.1f} |")

    print("\n### Z3: the identity residual rho_I(C4) = c/N_s, and what the same c would do to the L-scan under C\n")
    print("| eta | c (rho_I = c/N_s) | +- | chi2/dof | c / A_C (A_C = 2 N_s Delta_C) | predicted Delta(N_s 25) - Delta(N_s 100) "
          "[%] | measured (L-scan) [%] | +- | (measured - predicted)/sigma |")
    print("|---|---|---|---|---|---|---|---|---|")
    for lab in ("0.10", "0.39"):
        cc = [c for c in CS if c["lab"] == lab]
        b, sb, ch, dof = wls(np.array([[1 / c["Ns"]] for c in cc]), np.array([c["rho_c"] for c in cc]), np.array([c["s_rho_c"] for c in cc]))
        anc = [c for c in cc if c["scan"] == "H" and abs(c["H"] - 10) < 1e-9][0]; AC = 2 * anc["Ns"] * anc["DC"]
        pred = 0.5 * b[0] * (1 / 25 - 1 / 100); s_pred = 0.5 * sb[0] * (1 / 25 - 1 / 100)
        L = sorted([c for c in cc if c["scan"] == "L"], key=lambda c: c["Ns"]); meas = L[0]["D"] - L[1]["D"]; s_meas = math.hypot(L[0]["sD"], L[1]["sD"])
        print(f"| {lab} | {b[0]:+.4f} | {sb[0]:.4f} | {ch:.1f}/{dof} | {b[0] / AC:.2f} +- {sb[0] / AC:.2f} | {100 * pred:+.3f} +- {100 * s_pred:.3f} | "
              f"{100 * meas:+.3f} | {100 * s_meas:.3f} | {(meas - pred) / math.hypot(s_meas, s_pred):+.1f} |")
    print("\n(pi/8: the L-scan difference also carries the length effect of Z1, so the pi/8 row is confounded and shown for "
          "completeness only.)")


if __name__ == "__main__":
    main()
