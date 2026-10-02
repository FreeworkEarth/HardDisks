#!/usr/bin/env python3
"""##CHRIS 2026-10-12: efficiency-map analysis, per 261010 sec. 1 with amendments A1-A3 (sec. 1.8).
Committed BEFORE the map ran. Reads level5_effmap_20261010/<cell>/red_<seed>.csv only.

Per seed, at every sample (A3, the exact ledger; the check):
    W_in(t) = Delta KE_gas(t) + KE_div(t) + Delta E_spring(t)
with W_in = PistonWork, KE_div = M_s v^2/2 from W0_v, E_spring = k (x - x_eq)^2/2 from W0_x_sigma.
(The recorded SpringE is 0 in the first trace row -- a bookkeeping gap in the writer -- so the spring
energy is computed from the recorded x; SpringE is compared with it on rows >= 1 and reported.)
E_qs(L) + X is the KR decomposition of Delta KE_gas, with L(t) from the recorded SegEtas.

Observables (A1). s = -(x - x(0)), the divider displacement away from the gas, as in Level 3 v6.
    epsilon_mean    : window = last 3 T_w, T_w = 2 pi sqrt(M_s/(k + k_gas)) (Level 3's), control-corrected
    epsilon_settled : window = last tau_r (Mansour, ONE gas; 261010 sec. 1.8 table), control-corrected
    Delta E_spring(s) = F_s s + k s^2/2, F_s = k (x_eq - 30.5)   [exact for the spring]
    epsilon = Delta E_spring(s_corr) / <W_in>, ratio of means; jackknife over seeds + control SE.
Gate for epsilon_settled (A1 as implemented, sec. 1.8): coherent energy E_coh = k_eff V_c from the
excess variance of the seed-mean trajectory; PASS / FAIL / UNRESOLVED against 0.01 Delta E_spring.
Reproduction line (sec. 1.7): Level 3 v6's estimator verbatim on k0.5_M200_u0.05 vs 0.8705 +- 0.0049.
"""
import glob, json, math, os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import plot_speed_of_sound_edmd as sos
from paper2_effmap_amend_20261012 import reversible, tau_r, T_ad, KS, MS, US, XEQ_RUN, N, H, R, D, L0, Z, dZ

def e_dof(k, L, Tf):
    """INFERENCE (sec. 1.8 A3): energy the released divider's degree of freedom takes from the gas KE at
    equipartition -- kT/2 kinetic + (k/(k + k_S)) kT/2 in the spring (the rest of its potential share is
    gas compression energy, which stays in KE_gas). k_S = N kT_f (Z + eta Z' + Z^2)/L^2 at the window."""
    e = N * math.pi * R * R / (H * L); kS = N * Tf * (Z(e) + e * dZ(e) + Z(e) ** 2) / L ** 2
    return 0.5 * Tf * (1 + k / (k + kS))
REPO = os.path.dirname(os.path.dirname(HERE))
P = os.environ.get("EFFMAP_DIR") or os.path.join(os.path.dirname(HERE), "experiments_energy_transfer", "level5_effmap_20261010")
FIG = os.path.join(REPO, "0000_PLAN_OVERALL", "paper2_energytransfer", "experiments", "final")
S3, S3E = 0.8705, 0.0049                                  # Level 3 closing number (pooled, sigma)
L3 = 78.5; ETA3 = N * math.pi * R * R / (L3 * H)          # Level 3 v6 constants, verbatim
Z3 = float(sos.Z_kolafa_rottner_2006(np.array([ETA3]))[0]); DZ3 = float(sos.dZ_kolafa_rottner_2006(np.array([ETA3]))[0])
KGAS = N * (Z3 + ETA3 * DZ3 + Z3 * Z3) / L3 ** 2
ZAC = 2.2221                                              # sec. 1.1
COLS = ["Time", "KE_gas_total", "W0_x_sigma", "W0_v", "PistonWork", "PistonR_x_sigma", "SegEtas", "SpringE"]

def load(tag):
    out = []
    for f in sorted(glob.glob(os.path.join(P, tag, "red_*.csv"))):
        d = pd.read_csv(f, low_memory=False)
        miss = [c for c in COLS if c not in d.columns]
        if miss: sys.exit(f"ABORT {f}: columns absent {miss}")
        if abs(d["Time"].iloc[0]) > 1e-9: sys.exit(f"ABORT {f}: trace does not start at release (t0 = {d['Time'].iloc[0]})")
        out.append(d)
    return out

def ledger(d, k, M):
    x = d["W0_x_sigma"].to_numpy(float); v = d["W0_v"].to_numpy(float)
    Es = 0.5 * k * (x - XEQ_RUN[k]) ** 2; KEd = 0.5 * M * v * v
    W = d["PistonWork"].to_numpy(float); KE = d["KE_gas_total"].to_numpy(float)
    res = (W - W[0]) - (KE - KE[0]) - (KEd - KEd[0]) - (Es - Es[0])
    sp = np.abs(d["SpringE"].to_numpy(float)[1:] - Es[1:]).max()
    return float(np.abs(res).max()), float(sp)

def seg_L(d):
    eg = d["SegEtas"].astype(str).str.split(";").str[-1].astype(float).to_numpy()
    return N * math.pi * R * R / (H * eg)

def svec(ds, t0, t1):
    """per-seed mean of s = -(x - x(0)) over t0 <= t <= t1"""
    return np.array([-(d["W0_x_sigma"].to_numpy(float)[(d["Time"] >= t0).to_numpy() & (d["Time"] <= t1).to_numpy()].mean()
                       - d["W0_x_sigma"].iloc[0]) for d in ds])

def eps_window(ds, cs, k, win):
    """epsilon on the last `win` sigma-time of the record, control-corrected; jackknife over seeds."""
    tend = min(d["Time"].iloc[-1] for d in ds); t0 = tend - win
    s = svec(ds, t0, tend); sc = svec(cs, t0, tend)
    W = np.array([d["PistonWork"].iloc[-1] - d["PistonWork"].iloc[0] for d in ds])
    Fs = k * (XEQ_RUN[k] - 30.5); dE = lambda q: Fs * q + 0.5 * k * q * q
    sbar_c = sc.mean(); eps = dE(s.mean() - sbar_c) / W.mean()
    n = len(s); jk = np.array([dE(np.delete(s, i).mean() - sbar_c) / np.delete(W, i).mean() for i in range(n)])
    e_jk = math.sqrt((n - 1) / n * np.sum((jk - jk.mean()) ** 2))
    e_c = (Fs + k * (s.mean() - sbar_c)) * (sc.std(ddof=1) / math.sqrt(len(sc))) / W.mean()
    return dict(eps=eps, err=math.hypot(e_jk, e_c), s=s.mean() - sbar_c, W=W.mean(), dEs=dE(s.mean() - sbar_c),
                t0=t0, tend=tend, sc=sbar_c)

def ecoh(ds, win):
    tend = min(d["Time"].iloc[-1] for d in ds); t0 = tend - win
    X = np.array([d["W0_x_sigma"].to_numpy(float)[(d["Time"] >= t0).to_numpy() & (d["Time"] <= tend).to_numpy()] for d in ds])
    m = min(len(r) for r in X); X = np.array([r[:m] for r in X])
    Tf = np.mean([d["KE_gas_total"].to_numpy(float)[(d["Time"] >= t0).to_numpy()].mean() / N for d in ds])
    def est(A):
        n = len(A); vbar = A.mean(axis=0).var(); vi = A.var(axis=1).mean()
        Vc = (n * vbar - vi) / (n - 1); return Tf * Vc / max(vi - Vc, 1e-12)
    E = est(X); n = len(X); jk = np.array([est(np.delete(X, i, axis=0)) for i in range(n)])
    return E, math.sqrt((n - 1) / n * np.sum((jk - jk.mean()) ** 2)), Tf

def level3_estimator(ds, cs, u, M, k=0.5):
    """paper2_level3_v6_20260924.main(), section 1, verbatim logic (i0 = 0: the piston moves from t = 0)."""
    TW = 2 * math.pi / math.sqrt((k + KGAS) / M)
    n = min(len(d) for d in ds); t = ds[0]["Time"].to_numpy(float)[:n]
    S = np.array([-(d["W0_x_sigma"].to_numpy(float)[:n] - d["W0_x_sigma"].iloc[0]) for d in ds])
    nc = min(len(d) for d in cs); tc = cs[0]["Time"].to_numpy(float)[:nc]
    Sc = np.array([-(d["W0_x_sigma"].to_numpy(float)[:nc] - d["W0_x_sigma"].iloc[0]) for d in cs])
    tau = D / u; t0 = tau + 3 * TW; w = t >= t0
    per = S[:, w].mean(axis=1); wc = (tc >= 0) & (tc <= t[-1] - t0); pc = Sc[:, wc].mean(axis=1)
    sb, sbe = per.mean(), per.std(ddof=1) / math.sqrt(len(per)); sc_, sce = pc.mean(), pc.std(ddof=1) / math.sqrt(len(pc))
    return sb - sc_, math.hypot(sbe, sce), (t[-1] - t0) / TW

def main():
    Ph, rev = reversible(); TR = {(k, M): tau_r(k, M, rev)["tr"] for k in KS for M in MS}
    rows, LED, SPR, out = [], 0.0, 0.0, {}
    ctrl = {(k, M): load(f"ctrl_k{k}_M{M}") for k in KS for M in MS}
    print("### Ledger (A3): max |W_in - (Delta KE_gas + KE_div + Delta E_spring)| per cell, all samples, all seeds\n")
    for k in KS:
        for M in MS:
            cs = ctrl[(k, M)]
            for u in US:
                tag = f"k{k}_M{M}_u{u}"; ds = load(tag)
                if len(ds) < 2 or len(cs) < 2: print(f"  {tag}: {len(ds)} seeds, control {len(cs)} -- skipped"); continue
                lg = [ledger(d, k, M) for d in ds]; LED = max(LED, max(a for a, _ in lg)); SPR = max(SPR, max(b for _, b in lg))
                TW = 2 * math.pi / math.sqrt((k + KGAS) / M)
                em = eps_window(ds, cs, k, 3 * TW); es = eps_window(ds, cs, k, TR[(k, M)])
                E, Ee, Tf = ecoh(ds, TR[(k, M)]); thr = 0.01 * es["dEs"]
                gate = "PASS" if E + 2 * Ee < thr else ("FAIL" if E - 2 * Ee > thr else "UNRESOLVED")
                Lw = np.concatenate([seg_L(d)[(d["Time"] >= es["t0"]).to_numpy()] for d in ds]).mean()
                dKE = np.mean([d["KE_gas_total"].to_numpy(float)[(d["Time"] >= es["t0"]).to_numpy()].mean() - d["KE_gas_total"].iloc[0] for d in ds])
                Eqs = N * (T_ad(Lw) - 1)
                rows.append(dict(k=k, M=M, u=u, n=len(ds), led=max(a for a, _ in lg), em=em, es=es, E=E, Ee=Ee, thr=thr, gate=gate,
                                 Tf=Tf, Lw=Lw, dKE=dKE, Eqs=Eqs, X=dKE - Eqs, Ed=e_dof(k, Lw, Tf)))
    print(f"**max over all cells = {LED:.2e} kT**; recorded SpringE vs k(x - x_eq)^2/2 on rows >= 1: max {SPR:.2e} kT\n")
    print("### epsilon per cell (A1): ratio of means, control-corrected; errors jackknife + control SE\n")
    print("| k | M_s | u | seeds | <W_in> | s_corr (3 T_w) | **epsilon_mean** | s_corr (tau_r) | **epsilon_settled** | E_coh ± σ | 0.01 ΔE_spring | gate | epsilon_rev |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        print(f"| {r['k']} | {r['M']} | {r['u']} | {r['n']} | {r['em']['W']:.4f} | {r['em']['s']:.4f} | **{r['em']['eps']:.4f} ± {r['em']['err']:.4f}** | "
              f"{r['es']['s']:.4f} | **{r['es']['eps']:.4f} ± {r['es']['err']:.4f}** | {r['E']:.4f} ± {r['Ee']:.4f} | {r['thr']:.4f} | {r['gate']} | {rev[r['k']]['eps']:.4f} |")
    print("\n### KR decomposition of Delta KE_gas at settle (A3; a model split, not the check)\n")
    print("| k | M_s | u | L (SegEtas, last tau_r) | T_f = KE/N | Delta KE_gas | E_qs(L) | X = Delta KE_gas - E_qs | E_dof (INFERENCE) | X + E_dof | Z_ac u d |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        print(f"| {r['k']} | {r['M']} | {r['u']} | {r['Lw']:.4f} | {r['Tf']:.5f} | {r['dKE']:.4f} | {r['Eqs']:.4f} | {r['X']:+.4f} | {r['Ed']:.4f} | {r['X']+r['Ed']:+.4f} | {ZAC*r['u']*D:.4f} |")
    print("\n### Prediction 1: epsilon -> epsilon_rev as u -> 0 (u = 0.01 and 0.02 cells)\n")
    print("| k | M_s | u | epsilon_mean | epsilon_rev | (eps - eps_rev)/sigma | epsilon_settled | (eps_s - eps_rev)/sigma |")
    print("|---|---|---|---|---|---|---|---|")
    for r in rows:
        if r["u"] > 0.02: continue
        er = rev[r["k"]]["eps"]
        print(f"| {r['k']} | {r['M']} | {r['u']} | {r['em']['eps']:.4f} | {er:.4f} | {(r['em']['eps']-er)/r['em']['err']:+.1f} | "
              f"{r['es']['eps']:.4f} | {(r['es']['eps']-er)/r['es']['err']:+.1f} |")
    print("\n### Prediction 2 (A2): flat for u <~ 0.09, then roughly linear; slope on u in {0.1, 0.2, 0.5}\n")
    print("| k | M_s | slope d eps_mean/du ± σ | -eps_rev Z_ac d / W_rev (order of magnitude) | mean eps_mean, u <= 0.05 |")
    print("|---|---|---|---|---|")
    for k in KS:
        for M in MS:
            rr = [r for r in rows if r["k"] == k and r["M"] == M]
            hi = [r for r in rr if r["u"] >= 0.1]; lo = [r for r in rr if r["u"] <= 0.05]
            if len(hi) < 3: continue
            x = np.array([r["u"] for r in hi]); y = np.array([r["em"]["eps"] for r in hi]); s = np.array([r["em"]["err"] for r in hi])
            w = 1 / s ** 2; A = np.vstack([x, np.ones_like(x)]).T; C = np.linalg.inv(A.T @ (A * w[:, None])); p = C @ (A.T @ (w * y))
            print(f"| {k} | {M} | {p[0]:+.4f} ± {math.sqrt(C[0,0]):.4f} | {-rev[k]['eps']*ZAC*D/rev[k]['W']:+.4f} | "
                  f"{np.mean([r['em']['eps'] for r in lo]):.4f} |")
    print("\n### Prediction 3: no-push controls (W_in = 0; thermal drift of x; floor on E_spring)\n")
    print("| k | M_s | seeds | max |PistonWork| | <s> whole run | ΔE_spring(<s>) | Delta KE_gas, last half | -E_dof predicted |")
    print("|---|---|---|---|---|---|---|---|")
    for k in KS:
        for M in MS:
            cs = ctrl[(k, M)]
            if not cs: continue
            mw = max(np.abs(d["PistonWork"].to_numpy(float)).max() for d in cs)
            sm = np.mean([-(d["W0_x_sigma"].mean() - d["W0_x_sigma"].iloc[0]) for d in cs]); Fs = k * (XEQ_RUN[k] - 30.5)
            h = [(d['KE_gas_total'].to_numpy(float)[len(d)//2:].mean() - d['KE_gas_total'].iloc[0]) for d in cs]
            Tc = np.mean([d['KE_gas_total'].to_numpy(float)[len(d)//2:].mean() / N for d in cs])
            print(f"| {k} | {M} | {len(cs)} | {mw:.2e} | {sm:+.4f} | {Fs*sm+0.5*k*sm*sm:+.4f} | {np.mean(h):+.4f} ± {np.std(h, ddof=1)/math.sqrt(len(h)):.4f} | {-e_dof(k, L0, Tc):+.4f} |")
    print("\n### Reproduction line (sec. 1.7): Level 3 v6 estimator on k0.5_M200_u0.05\n")
    ds, cs = load("k0.5_M200_u0.05"), ctrl[(0.5, 200)]
    if len(ds) > 1 and len(cs) > 1:
        s, e, nper = level3_estimator(ds, cs, 0.05, 200)
        z = abs(s - S3) / math.hypot(e, S3E)
        print(f"s̄ = **{s:.4f} ± {e:.4f} σ** ({nper:.0f} periods in window) vs Level 3 {S3} ± {S3E}: "
              f"**{z:.2f} σ -> {'PASS' if z < 2 else 'FAIL'}** (rule: within 2σ)")
        out["reproduction"] = dict(s=s, err=e, z=z)
    figure(rows, rev)
    out["ledger_max"] = LED; out["springE_vs_x_max"] = SPR
    out["cells"] = [dict(k=r["k"], M=r["M"], u=r["u"], n=r["n"], eps_mean=r["em"]["eps"], eps_mean_err=r["em"]["err"],
                         eps_settled=r["es"]["eps"], eps_settled_err=r["es"]["err"], W_in=r["em"]["W"], E_coh=r["E"],
                         E_coh_err=r["Ee"], gate=r["gate"], X=r["X"]) for r in rows]
    json.dump(out, open(os.path.join(HERE, "261012_effmap_results.json"), "w"), indent=1)

def figure(rows, rev):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.8), sharex=True)
    for ax, k in zip(axs, KS):
        for M, mfc, off in ((50, "tab:blue", 0.97), (200, "white", 1.03)):
            rr = [r for r in rows if r["k"] == k and r["M"] == M]
            if not rr: continue
            ax.errorbar([r["u"] * off for r in rr], [r["em"]["eps"] for r in rr], yerr=[r["em"]["err"] for r in rr], fmt="o",
                        color="tab:blue", mfc=mfc, ms=6, capsize=3, label=f"ε_mean, M_s = {M}")
        ax.axhline(rev[k]["eps"], color="red", lw=1.5, label=f"ε_rev = {rev[k]['eps']:.4f}")
        ax.axvline(0.0884, color="0.6", ls=":", lw=1, label="d/u = 2L₀/c_s")
        ax.set_xscale("log"); ax.set_title(f"k = {k}", fontsize=10); ax.set_xlabel("piston speed u")
    axs[0].set_ylabel("efficiency ε = ΔE_spring / W_in"); axs[0].legend(fontsize=7, frameon=False)
    fig.suptitle("Paper 2 efficiency map, geometry C (shifted ×0.97/×1.03 in u for legibility)", fontsize=10)
    fig.tight_layout()
    for ext in ("png", "pdf"): fig.savefig(os.path.join(FIG, f"261012_p2_effmap.{ext}"), dpi=200)
    print("\nfigure: 0000_PLAN_OVERALL/paper2_energytransfer/experiments/final/261012_p2_effmap.{png,pdf}")

if __name__ == "__main__":
    main()
