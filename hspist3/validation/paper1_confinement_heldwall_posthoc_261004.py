#!/usr/bin/env python3
"""##CHRIS 2026-10-04 (Task X, POST-HOC -- NOT part of the registered analysis, no verdict): the "held" divider of method A
is not held, and what that does to the registered k_T.

FINDING (DATA). Method A holds the divider for 12000 steps (200 sigma-time), then releases it with mass factor 1e9 for the
5000 sigma-time window of the force measurement (conf_worker.sh mode A: --wall-hold-steps=12000 --wall-mass-factors=
1000000000; reduce_A.py: window = [200, t_last]). The registration's drift estimate (261012 sec. 1.4, Table A) treated only
the random thermal drift of the released divider. At an off-centre position x_j = j dL the net gas force is restoring,
so the divider returns toward the centre deterministically. The full pilot cell (pi/8 anchor, traces fetched) shows it:
the window-mean divider displacement is a fixed fraction f of nominal in every off-centre run.

CONSEQUENCE (DERIVATION, linear order). Each compartment is closed (the divider is its only energy channel), so the slow
return compresses or expands each gas adiabatically. With delta_j = xbar_j - x_j:
    time-averaged force  Fbar_L(j) = F_T(L_0 + x_j) - k_S delta_j      (F_T: the isothermal force at T = 1)
    energy balance       N_s k (Tbar_L(j) - 1) = -Fbar_L(j) delta_j  ->  delta_j = -N_s (Tbar_L - 1)/Fbar_L  (left gas)
                                                                   delta_j = +N_s (Tbar_R - 1)/Fbar_R  (right gas)
so the registered stencil (nominal positions, raw forces) returns k_T,meas = k_T - k_S (1 - f): biased LOW.
For hard disks F = T g(L), so Fbar/Tbar = g(Lbar) to first order: the force normalised to T = 1 at the MEASURED mean
position is the isothermal force there. The corrected estimator therefore uses, per seed, F/T, and the stencil spacing
f dL, with f measured in every cell from the recorded temperatures (T_L, T_R in red_<seed>.csv, already on the Mac):
    1 - f = -sum_j delta_j x_j / sum_j x_j^2 over j = -2, -1, +1, +2 (delta_j = mean of the left- and right-gas values)
    k_T,c = -[Ft(-2) - 8 Ft(-1) + 8 Ft(+1) - Ft(+2)]/(12 f dL),  Ft(j) = mean of (F_L/T_L)(run j) and (F_R/T_R)(run -j)
Cross-check at the pi/8 anchor: f from the temperatures vs f from the pilot's divider trajectories (W0_x_sigma).
Everything below is post-hoc and exploratory: the registered verdicts are those of paper1_confinement_results_261004.py.
usage (from hspist3/): python3 validation/paper1_confinement_heldwall_posthoc_261004.py
"""
import glob, math, os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, HS)
import contextlib, io
import paper1_confinement_results_261004 as R

PILOT = os.path.join(HS, R.REL_A, "pilot_epi8_H_H10_L10")


def pilot_positions():
    """f from the pilot's divider trajectories (the only method-A traces on the Mac)."""
    L0, dL = 10.0, 0.125; out = []
    for lab, j in R.POS:
        if j == 0: continue
        for f in sorted(glob.glob(os.path.join(PILOT, f"x_{lab}", "tr_*.csv"))):
            d = pd.read_csv(f, usecols=["Time", "W0_x_sigma"]); t = d["Time"].to_numpy(); x = d["W0_x_sigma"].to_numpy() - L0
            w = t >= 200; out.append(x[w].mean() / (j * dL))
    return np.array(out)


def drift(c, cell_dir=None):
    """Per position: per-seed F/T (both faces), mean temperatures, delta_j from the energy balance; f by least squares."""
    ad = cell_dir or os.path.join(HS, R.REL_A, c["cid"]); st = {}
    for lab, j in R.POS:
        files = sorted(glob.glob(os.path.join(ad, f"x_{lab}", "red_*.csv"))) if cell_dir else \
            [os.path.join(ad, f"x_{lab}", f"red_{s}.csv") for s in c["seeds"][lab]]
        D = pd.concat([pd.read_csv(f) for f in files], ignore_index=True); n = len(D)
        gL, gR = D["F_L"] / D["T_L"], D["F_R"] / D["T_R"]
        st[j] = dict(n=n, TL=D["T_L"].mean(), TR=D["T_R"].mean(), FL=D["F_L"].mean(), FR=D["F_R"].mean(),
                     gL=gL.mean(), gR=gR.mean(), sgL=gL.std(ddof=1) / math.sqrt(n), sgR=gR.std(ddof=1) / math.sqrt(n))
    dL, Ns = c["dL"], c["Ns"]; num = den = 0.0; dl = {}
    for j in (-2, -1, 1, 2):
        dLft = -Ns * (st[j]["TL"] - 1) / st[j]["FL"]; dRgt = Ns * (st[j]["TR"] - 1) / st[j]["FR"]
        dl[j] = (dLft, dRgt); d = 0.5 * (dLft + dRgt); x = j * dL; num += d * x; den += x * x
    one_f = -num / den; f = 1 - one_f
    fj = {j: 1 + 0.5 * (dl[j][0] + dl[j][1]) / (j * dL) for j in dl}
    s_f = float(np.std(list(fj.values()), ddof=1) / math.sqrt(len(fj)))
    Ft = {j: 0.5 * (st[j]["gL"] + st[-j]["gR"]) for j in (-2, -1, 0, 1, 2)}
    sFt = {j: 0.5 * math.hypot(st[j]["sgL"], st[-j]["sgR"]) for j in (-2, -1, 0, 1, 2)}
    kTc = -(Ft[-2] - 8 * Ft[-1] + 8 * Ft[1] - Ft[2]) / (12 * f * dL)
    s_kTc = math.sqrt(sFt[-2] ** 2 + 64 * sFt[-1] ** 2 + 64 * sFt[1] ** 2 + sFt[2] ** 2) / (12 * f * dL)
    s_kTc = math.hypot(s_kTc, kTc * s_f / f)
    return dict(f=f, s_f=s_f, fj=fj, kTc=kTc, s_kTc=s_kTc, st=st)


def main():
    with contextlib.redirect_stdout(io.StringIO()):
        CS = R.cells()
        for c in CS:
            R.method_B(c); R.method_A(c); R.identity(c)
    print("## POST-HOC (not registered, no verdict): the released 'held' divider of method A\n")
    fp = pilot_positions()
    print(f"pilot traces (pi/8 anchor, 16 off-centre runs, W0_x_sigma over the window [200, end]): "
          f"mean displacement / nominal = f = {fp.mean():.4f} +- {fp.std(ddof=1):.4f} (min {fp.min():.4f}, max {fp.max():.4f})")
    anc = [c for c in CS if c["cid"] == "epi8_H_H10_L10"][0]
    dp = drift(anc, cell_dir=PILOT)
    print(f"same pilot runs, f from the temperatures (energy balance): {dp['f']:.4f} +- {dp['s_f']:.4f}  "
          f"(per position: " + ", ".join(f"{j:+d}dL {v:.4f}" for j, v in sorted(dp['fj'].items())) + ")")
    print("\n### Per cell: f from the temperatures, and the identity with the drift-corrected k_T (exploratory)\n")
    print("| eta | cell | N_s | f (T balance) | sigma_f | k_T registered | k_T,c (F/T, spacing f dL) | sigma | k_T,c/k_T - 1 [%] "
          "| static,c = k_T,c + F^2/(N_s kT) | k_S^dyn | rho_I,c [%] | sigma [%] | rho_I,c/sigma | registered rho_I [%] | 2 Delta_C [%] "
          "| gamma_box,c | bulk gamma |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    rows = []
    for c in CS:
        d = drift(c); static = d["kTc"] + c["F2term"]
        rho = (c["kS"] - static) / c["kS"]; s = math.sqrt((c["s_kS"] / c["kS"]) ** 2 + (d["s_kTc"] / c["kS"]) ** 2)
        g = c["kS"] / d["kTc"]
        rows.append(dict(c=c, d=d, rho=rho, s=s, z=rho / s, g=g))
        print(f"| {c['lab']} | {c['cid']} | {c['Ns']} | {d['f']:.4f} | {d['s_f']:.4f} | {c['kT']:.5f} | {d['kTc']:.5f} | {d['s_kTc']:.5f} | "
              f"{100 * (d['kTc'] / c['kT'] - 1):+.2f} | {static:.5f} | {c['kS']:.5f} | {100 * rho:+.3f} | {100 * s:.3f} | {rho / s:+.2f} | "
              f"{100 * c['rho']:+.3f} | {200 * c['DC']:+.3f} | {g:.4f} | {c['gbulk']:.5f} |")
    for lab in ("0.10", "0.39"):
        rr = [r for r in rows if r["c"]["lab"] == lab]
        n2 = sum(abs(r["z"]) <= 2 for r in rr); chi = sum(r["z"] ** 2 for r in rr)
        w = np.array([1 / r["s"] ** 2 for r in rr]); m = float((w * np.array([r["rho"] for r in rr])).sum() / w.sum())
        print(f"\neta {lab}: drift-corrected rho_I within 2 sigma in {n2} of {len(rr)} cells; sum (rho/sigma)^2 = {chi:.1f} "
              f"({len(rr)} cells); inverse-variance mean rho_I,c = {100 * m:+.3f} +- {100 / math.sqrt(w.sum()):.3f} %")
    # rho_I,c against hypothesis C's fixed signature 2 Delta_C (one free amplitude, information only)
    for lab in ("0.10", "0.39"):
        rr = [r for r in rows if r["c"]["lab"] == lab]
        y = np.array([r["rho"] for r in rr]); s = np.array([r["s"] for r in rr]); f = np.array([2 * r["c"]["DC"] for r in rr])
        w = 1 / s ** 2; a = float((w * f * y).sum() / (w * f * f).sum()); sa = 1 / math.sqrt(float((w * f * f).sum()))
        ch = float((w * (y - a * f) ** 2).sum()); ch0 = float((w * y * y).sum())
        print(f"eta {lab}: rho_I,c = a x 2 Delta_C: a = {a:.2f} +- {sa:.2f}, chi2 {ch:.1f} / {len(y) - 1} dof "
              f"(rho_I,c = 0: chi2 {ch0:.1f} / {len(y)} dof)")
    figure(rows)


def figure(rows):
    """Post-hoc figure: the pilot's divider trajectories, 1 - f against k_S per cell, rho_I before and after."""
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(17, 5.6), gridspec_kw=dict(width_ratios=[1, 1, 1.5]))
    for lab, j in R.POS:
        if j == 0: continue
        for f in sorted(glob.glob(os.path.join(PILOT, f"x_{lab}", "tr_*.csv"))):
            d = pd.read_csv(f, usecols=["Time", "W0_x_sigma"])
            a1.plot(d["Time"], (d["W0_x_sigma"] - 10.0) / (j * 0.125), color=R.BLUE, lw=0.6, alpha=0.6)
    a1.axhline(1, color="k", ls="--", lw=1); a1.axvline(200, color="0.5", ls=":", lw=1)
    a1.set_xlabel("time [sigma-time]  (released at 200)"); a1.set_ylabel("divider displacement / nominal x_j")
    a1.set_title("pi/8 anchor pilot: 16 off-centre runs (DATA)", fontsize=10); a1.grid(True, ls=":", alpha=0.5)
    for lab, col, mk in (("0.10", R.BLUE, "o"), ("0.39", R.BLUE2, "s")):
        rr = [r for r in rows if r["c"]["lab"] == lab]
        a2.errorbar([r["c"]["kS"] for r in rr], [1 - r["d"]["f"] for r in rr], yerr=[r["d"]["s_f"] for r in rr], fmt=mk, color=col,
                    ms=6, capsize=2, label=("eta = 0.100" if lab == "0.10" else "eta = pi/8") + ": 1 - f from the T balance")
    fp = pilot_positions(); anc = [r for r in rows if r["c"]["cid"] == "epi8_H_H10_L10"][0]
    a2.errorbar([anc["c"]["kS"]], [1 - fp.mean()], yerr=[fp.std(ddof=1)], fmt="*", color="k", ms=11, capsize=3,
                label="pilot trajectories (pi/8 anchor)")
    a2.set_xscale("log"); a2.set_yscale("log"); a2.grid(True, which="both", ls=":", alpha=0.5)
    a2.set_xlabel("k_S^dyn of the cell"); a2.set_ylabel("1 - f  (shortfall of the mean displacement)")
    a2.set_title("drift of the released divider, every cell (DATA)", fontsize=10); a2.legend(fontsize=8)
    order = sorted(rows, key=lambda r: (r["c"]["lab"], r["c"]["scan"], r["c"]["Ns"], r["c"]["L0"]))
    for i, r in enumerate(order):
        c = r["c"]; col = R.BLUE if c["lab"] == "0.10" else R.BLUE2; mk = "o" if c["lab"] == "0.10" else "s"
        a3.errorbar(i - 0.15, 100 * c["rho"], yerr=100 * c["s_rho"], fmt=mk, mfc="white", color=col, ms=6, capsize=2, lw=1)
        a3.errorbar(i + 0.15, 100 * r["rho"], yerr=100 * r["s"], fmt=mk, color=col, ms=6, capsize=2, lw=1.3)
        a3.plot(i, 200 * c["DC"], marker="_", color="#8f4fd1", ms=14, mew=2)
    a3.axhline(0, color="k", ls="--", lw=1)
    a3.plot([], [], "o", mfc="white", color=R.BLUE, label="registered k_T (nominal positions)")
    a3.plot([], [], "o", color=R.BLUE, label="drift-corrected k_T (F/T, spacing f dL) -- post-hoc")
    a3.plot([], [], marker="_", color="#8f4fd1", ls="none", ms=14, mew=2, label="hypothesis C: rho_I = 2 Delta_C")
    a3.set_xticks(range(len(order))); a3.set_xticklabels([r["c"]["cid"].replace("e0p10_", "0.10 ").replace("epi8_", "pi/8 ")
                                                          for r in order], rotation=70, ha="right", fontsize=7)
    a3.set_ylabel("rho_I = (k_S^dyn - static)/k_S^dyn  [%]"); a3.grid(True, ls=":", alpha=0.5); a3.legend(fontsize=8, loc="upper left")
    a3.set_title("identity residual: registered vs drift-corrected (POST-HOC, no verdict)", fontsize=10)
    fig.suptitle("POST-HOC: the released 'held' divider of method A returns toward the centre (not part of the registered analysis)",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(R.OUT, f"261004_p1_identity_heldwall_posthoc.{ext}"), dpi=200)
    print("\nfigure -> 261004_p1_identity_heldwall_posthoc.png/.pdf")


if __name__ == "__main__":
    main()
