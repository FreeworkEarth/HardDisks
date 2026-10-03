#!/usr/bin/env python3
"""##CHRIS 2026-10-13: over-pressure test for the efficiency map (261010 sec. 3, pre-registered before this ran).

Question: is the quasi-static excess eps_settled/eps_rev(KR) = 1.023 / 1.011 / 1.021 (k = 0.25 / 0.5 / 1.0;
261010 sec. 2.2 (a)) the box's over-pressure? Recompute eps_rev with the box's OWN measured equation of state,
Level 3's held-wall force curve F(L), instead of Kolafa-Rottner, and compare.

F(L) source -- recomputed here from the raw files, not typed (the 260925 table was made inline and its window
cannot be recovered; the comparison with it is printed):
  experiments_energy_transfer/level3_FofL_20260925/c{0,2.5,5,7.5,10}/ev_<seed>.csv  (D0 = divider events)
  experiments_energy_transfer/level3_FofL_20260925/c{...}/tr_<seed>.csv              (KE_gas_total)
  dp is the particle's momentum change, written at edmd.c:1069:
      fprintf(g_edmd_evlog, "%.12g,%s,%.12g,%.12g,%.12g,%.12g,%.12g\\n", S->t / g_edmd_evlog_tscale, kind, u, v0, v1, dE, 1.0 * (v1 - v0));
  KE_gas_total is the sum of segment kinetic energies, written at 00ALLINONE.c:17186:
      fprintf(elog, ",%.12e,%.12e,%.12e,%.12e", ke_tot, ke_left, ke_right, px_gas);
  Event times include the 12000-step hold (200 sigma); trace times start at release.
Window (pre-registered): from t_stop + 2 (2 L_0/c_s) = t_stop + 180 sigma to the record end, t_stop = 0.25/u + d/u
(c0: from 180). Per seed F = sum|dp|/(window length), T = <KE_gas_total>/N; Z_box = F L/(N T), L = 78.5 - d.
EOS model: lambda(eta) = Z_box/Z_KR, weighted linear fit in eta over the five points; Z_F = lambda Z_KR.
"""
import glob, json, math, os, sys
import numpy as np, pandas as pd
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.stats import chi2 as chi2d
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
from paper2_effmap_amend_20261012 import Z, dZ, reversible, tau_r, KS, MS, US, XEQ_RUN, N, H, R, D, L0, ETA0
import paper2_effmap_analysis_20261012 as A
REPO = os.path.dirname(os.path.dirname(HERE))
L3 = os.path.join(os.path.dirname(HERE), "experiments_energy_transfer", "level3_FofL_20260925")
FIG = os.path.join(REPO, "0000_PLAN_OVERALL", "paper2_energytransfer", "experiments", "final")
HOLD, END, U3 = 200.0, 666.65, 0.05
CELLS = [("c0", 0.0), ("c2.5", 1.99), ("c5", 3.98), ("c7.5", 5.97), ("c10", 7.96)]
PUB = {"c0": (1.6087, 0.0064, 1.0000, 1.2628), "c2.5": (1.7218, 0.0078, 1.0347, 1.2732), "c5": (1.8366, 0.0067, 1.0699, 1.2792),
       "c7.5": (1.9715, 0.0059, 1.1070, 1.2917), "c10": (2.1217, 0.0042, 1.1499, 1.3016)}     # 260925 table (F, sF, T, Z_box)
XWIN = 180.0                                                                                # 2 x 2L_0/c_s (261010 sec. 1.8)

def eta(L): return N * math.pi * R * R / (H * L)

def level3(xwin=XWIN):
    out = []
    for c, d in CELLS:
        t0 = (0.25 / U3 + d / U3 if d > 0 else 0.0) + xwin; L = 78.5 - d; F, T = [], []
        for f in sorted(glob.glob(os.path.join(L3, c, "ev_*.csv"))):
            sd = f.split("_")[-1][:-4]
            ev = pd.read_csv(f, usecols=["t_sigma", "kind", "dp"]); dv = ev[ev["kind"] == "D0"]
            te = dv["t_sigma"].to_numpy() - HOLD; m = (te >= t0) & (te <= END)
            F.append(np.abs(dv["dp"].to_numpy()[m]).sum() / (END - t0))
            tr = pd.read_csv(os.path.join(L3, c, f"tr_{sd}.csv"), usecols=["Time", "KE_gas_total"], low_memory=False)
            w = (tr["Time"] >= t0).to_numpy(); T.append(tr["KE_gas_total"].to_numpy()[w].mean() / N)
        F, T = np.array(F), np.array(T); Zs = F * L / (N * T); e = eta(L)
        out.append(dict(c=c, d=d, L=L, eta=e, F=F.mean(), sF=F.std(ddof=1) / math.sqrt(len(F)), T=T.mean(),
                        Z=Zs.mean(), sZ=Zs.std(ddof=1) / math.sqrt(len(Zs)), lam=Zs.mean() / Z(e), slam=Zs.std(ddof=1) / math.sqrt(len(Zs)) / Z(e), n=len(F)))
    return out

def fit_lambda(rows):
    x = np.array([r["eta"] - ETA0 for r in rows]); y = np.array([r["lam"] for r in rows]); s = np.array([r["slam"] for r in rows])
    w = 1 / s ** 2; X = np.vstack([np.ones_like(x), x]).T; C = np.linalg.inv(X.T @ (X * w[:, None])); p = C @ (X.T @ (w * y))
    return p, C, float(np.sum(w * (y - X @ p) ** 2))

def construct(k, a, b, Ti=1.0, variant="matched"):
    lam = lambda e: a + b * (e - ETA0)                       # <- the interpolation of Level 3's F(L) (pre-registered line)
    ZF = lambda e: lam(e) * Z(e)
    Tad = lambda L: Ti * math.exp(-quad(lambda y: ZF(eta(math.exp(y))), math.log(L0), math.log(L))[0])
    F = lambda L: N * Tad(L) * ZF(eta(L)) / L
    xeq = XEQ_RUN[k]; Fs = k * (xeq - 30.5)
    Lf = brentq(lambda L: k * (xeq - (109.0 - D - L)) - F(L), 60, L0); xf = 109.0 - D - Lf
    if variant == "literal":                                 # 261010 sec. 1.3 verbatim: start at x = 30.5, L_0, T = Ti
        Esp = 0.5 * k * (xf - xeq) ** 2 - 0.5 * k * (30.5 - xeq) ** 2
        return Esp / (N * (Tad(Lf) - Ti) + Esp)
    Li = brentq(lambda L: k * (xeq - (109.0 - L)) - F(L), L0 - 3, L0 + 3); xi = 109.0 - Li
    s = xi - xf; dEest = Fs * s + 0.5 * k * s * s                 # the estimator's own spring formula (sec. 1.8 A1)
    W = N * (Tad(Lf) - Tad(Li)) + 0.5 * k * (xf - xeq) ** 2 - 0.5 * k * (xi - xeq) ** 2
    return dEest / W

def measured_plateau():
    Ph, rev = reversible(); TR = {(k, M): tau_r(k, M, rev)["tr"] for k in KS for M in MS}
    out, Tc = {}, {}
    for k in KS:
        v, s, tc = [], [], []
        for M in MS:
            cs = A.load(f"ctrl_k{k}_M{M}")
            tc += [d["KE_gas_total"].to_numpy(float)[len(d) // 2:].mean() / N for d in cs]
            for u in US:
                if u > 0.05: continue
                e = A.eps_window(A.load(f"k{k}_M{M}_u{u}"), cs, k, TR[(k, M)]); v.append(e["eps"]); s.append(e["err"])
        w = 1 / np.array(s) ** 2; out[k] = (float(np.sum(w * np.array(v)) / np.sum(w)), float(1 / math.sqrt(np.sum(w))))
        Tc[k] = float(np.mean(tc))
    return out, Tc, rev

def main():
    rows = level3()
    print("### 1. Level 3 F(L), recomputed from the raw files (window t_stop + 180 -> end), against the 260925 table\n")
    print("| cell | L | eta | F (here) | F (260925) | T (here) | T (260925) | Z_box (here) | Z_box (260925) | lambda = Z_box/Z_KR |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        p = PUB[r["c"]]
        print(f"| {r['c']} | {r['L']:.2f} | {r['eta']:.5f} | {r['F']:.4f} ± {r['sF']:.4f} | {p[0]} ± {p[1]} | {r['T']:.4f} | {p[2]} | "
              f"{r['Z']:.4f} ± {r['sZ']:.4f} | {p[3]} | {r['lam']:.4f} ± {r['slam']:.4f} |")
    p, C, ch = fit_lambda(rows)
    print(f"\nlambda(eta) = {p[0]:.5f} ± {math.sqrt(C[0,0]):.5f} + ({p[1]:+.3f} ± {math.sqrt(C[1,1]):.3f})(eta - {ETA0:.8f}); "
          f"chi2 = {ch:.2f} on 3 dof; pooled constant lambda = {np.average([r['lam'] for r in rows], weights=[1/r['slam']**2 for r in rows]):.5f}")
    print("\nRobustness of the fit to the window start (pre-registered check):\n")
    print("| window start after t_stop | a = lambda(eta_0) | b |")
    print("|---|---|---|")
    for xw in (100.0, 180.0, 300.0):
        pp, CC, _ = fit_lambda(level3(xw)); print(f"| {xw:.0f} | {pp[0]:.5f} ± {math.sqrt(CC[0,0]):.5f} | {pp[1]:+.3f} ± {math.sqrt(CC[1,1]):.3f} |")
    meas, Tc, rev = measured_plateau()
    rng = np.random.default_rng(20261013); draws = rng.multivariate_normal(p, C, size=400)
    print("\n### 2. The test (primary: estimator-matched construction from the box's own pre-push equilibrium)\n")
    print("| k | eps_rev KR (sec. 1.3) | eps_rev KR, matched | **eps_rev F(L)** | eps_settled (u <= 0.05) | ratio_KR | **ratio_F** | (ratio_F - 1)/σ | verdict | eps_rev F/KR |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    res = {}; allpass = True
    for k in KS:
        eKR = rev[k]["eps"]; eKRm = construct(k, 1.0, 0.0); eF = construct(k, p[0], p[1])
        sF = float(np.std([construct(k, a, b) for a, b in draws], ddof=1))
        m, sm = meas[k]; rK = m / eKR; rF = m / eF; srF = rF * math.hypot(sm / m, sF / eF); z = (rF - 1) / srF
        ok = abs(z) <= 2; allpass &= ok
        res[k] = dict(eKR=eKR, eKRm=eKRm, eF=eF, seF=sF, m=m, sm=sm, rK=rK, srK=rK * sm / m, rF=rF, srF=srF, z=z, ok=ok)
        print(f"| {k} | {eKR:.5f} | {eKRm:.5f} | **{eF:.5f} ± {sF:.5f}** | {m:.5f} ± {sm:.5f} | {rK:.4f} ± {rK*sm/m:.4f} | "
              f"**{rF:.4f} ± {srF:.4f}** | {z:+.1f} | {'within 2σ' if ok else 'OUTSIDE 2σ'} | {eF/eKR:.4f} |")
    print(f"\n**VERDICT (pre-registered): {'PASS -- the over-pressure explanation becomes DATA' if allpass else 'FAIL -- stays OPEN (1-2 % unexplained offset)'}**")
    rk = np.array([res[k]["rK"] for k in KS]); sk = np.array([res[k]["srK"] for k in KS]); w = 1 / sk ** 2
    mu = float(np.sum(w * rk) / np.sum(w)); c2 = float(np.sum(w * (rk - mu) ** 2))
    print(f"\nConsistency check (pre-registered): measured ratio_KR against one common value: mean {mu:.4f} ± {1/math.sqrt(np.sum(w)):.4f}, "
          f"chi2 = {c2:.2f} on 2 dof, p = {chi2d.sf(c2, 2):.3g}")
    print("\n### 3. Secondary variants (no verdict)\n")
    print("| k | literal sec. 1.3 construction with F(L): eps_rev | ratio | matched + measured control T_i: T_i | eps_rev | ratio |")
    print("|---|---|---|---|---|---|")
    for k in KS:
        eL = construct(k, p[0], p[1], variant="literal"); eT = construct(k, p[0], p[1], Ti=Tc[k]); m = res[k]["m"]
        print(f"| {k} | {eL:.5f} | {m/eL:.4f} | {Tc[k]:.5f} | {eT:.5f} | {m/eT:.4f} |")
    figure(res)
    json.dump(dict(lambda_fit=dict(a=p[0], b=p[1], cov=C.tolist()), verdict="PASS" if allpass else "FAIL", common=dict(mean=mu, chi2=c2),
                   per_k={str(k): v for k, v in res.items()}), open(os.path.join(HERE, "261013_effmap_overpressure.json"), "w"), indent=1, default=float)

def figure(res):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(5.6, 4.0)); ks = np.array(KS)
    ax.errorbar(ks * 0.97, [res[k]["rK"] for k in KS], yerr=[res[k]["srK"] for k in KS], fmt="o", color="red", capsize=3,
                label=r"$\varepsilon_{\rm settled}/\varepsilon_{\rm rev}$, KR equation of state")
    ax.errorbar(ks * 1.03, [res[k]["rF"] for k in KS], yerr=[res[k]["srF"] for k in KS], fmt="s", color="tab:blue", capsize=3,
                label=r"$\varepsilon_{\rm settled}/\varepsilon_{\rm rev}$, box $F(L)$ (Level 3)")
    ax.axhline(1.0, color="0.3", lw=1); ax.set_xscale("log"); ax.set_xticks(KS); ax.set_xticklabels([str(k) for k in KS]); ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel(r"spring constant $k$ [$kT/\sigma^2$]"); ax.set_ylabel("quasi-static plateau / reversible reference")
    ax.set_title("Over-pressure test (points offset in k for legibility)", fontsize=10); ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    for ext in ("png", "pdf"): fig.savefig(os.path.join(FIG, f"261013_p2_effmap_overpressure.{ext}"), dpi=200)
    print("\nfigure: 0000_PLAN_OVERALL/paper2_energytransfer/experiments/final/261013_p2_effmap_overpressure.{png,pdf}")

if __name__ == "__main__":
    main()
