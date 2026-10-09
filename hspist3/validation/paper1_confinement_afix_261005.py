#!/usr/bin/env python3
"""##CHRIS 2026-10-05: the PRE-REGISTERED analysis of the "A-fixed" campaign (261012 sec. 3), written and committed BEFORE
the A-fixed summaries reached the Mac, applied once.

Registration (261012 sec. 3.4-3.6), implemented here:
  per seed F_L/T_L and F_R/T_R (reduce_AF.py, held window [200, 5200)); Ft(j) = mean of (F_L/T_L)(run j) and
  (F_R/T_R)(run -j); k_T = -[Ft(-2) - 8 Ft(-1) + 8 Ft(+1) - Ft(+2)]/(12 dL) with the nominal dL (f = 1); sigma from the
  seed standard errors. Static side k_T + F(L0)^2/(N_s kT), F(L0) = (F_L + F_R)/2 at x = 0 (raw), kT = the mean
  temperature of the x = 0 seeds. Dynamic side k_S^dyn from method B (C1, unchanged: paper1_confinement_results_261004).
  rho_I = (k_S^dyn - static)/k_S^dyn, sigma(rho_I)^2 = (sigma_kS/k_S)^2 + (sigma_kT/k_S)^2.
  G2 inventory: every seed present; summaries record hold 312000, tail 1200, build 279282b, the cell's geometry; health 0;
     window = 5000, t_last >= 5200.
  G3 per cell: u_wall_max = 0 and W_div = 0 for every seed; |1 - f| < 0.002 with f from the temperatures by the C4 formula
     (paper1_confinement_heldwall_posthoc_261004.drift). A failing cell is flagged and left out of P1-P3; more than two
     flagged cells -> STOP.
  P1: per cell, rho_I(A-fixed) - rho_I(C4) within 2 sigma, sigma^2 = (sigma_kT,AF^2 + sigma_kT,C4^2)/k_S^2; holds if every
      cell does; chi2 reported.
  P2: |rho_I(A-fixed)| <= 2 sigma at every cell (the C1 rule).
  P3: per density, weighted fit rho_I = c/N_s; r = c/A_C, A_C = 2 N_s Delta_C at the anchor (0.634, 0.384); chi2 of the fit
      and of rho_I = 0; the outcomes declared in sec. 3.6 are evaluated (C4's r for comparison: the same fit to rho_I(C4)).
Outputs: tables on stdout; paper1_speedofsound/experiments/final/261005_p1_identity_afix.{png,pdf},
261005_p1_identity_afix_cells.csv.
usage (from hspist3/): python3 validation/paper1_confinement_afix_261005.py   [--test-pilot: mechanics test on the pilot cell
       only (the anchor, seeds 9700-9703; gate-G1 data), prints, writes nothing]
"""
import contextlib, glob, io, math, os, re, sys
import numpy as np, pandas as pd
from scipy.stats import chi2 as CHI2
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, HS)
import paper1_confinement_results_261004 as R
import edmd_acc_guard   # ##CHRIS 2026-10-08 (261012 sec. 4.7.4, decision 2): the loader provenance guard (full name: no alias can be shadowed)
import paper1_confinement_heldwall_posthoc_261004 as PH

REL_AF = "experiments_energy_transfer/paper1_confinement_Afix_261004"
BUILD, HOLD, POST, T0, T1 = "279282b", 312000, 1200, 200.0, 312000 * 0.4 / 24.0
HEALTH = re.compile(r"EDMD-HEALTH")
TEST = "--test-pilot" in sys.argv


def af_dir(c):
    return os.path.join(HS, REL_AF, "pilot_epi8_H_H10_L10" if TEST else c["cid"])


def inventory_and_static(c):
    """G2 per cell, then the registered estimator; returns a dict (or the reason the cell cannot be computed)."""
    d0 = af_dir(c); st = {}; miss = geo = health = uw = wd = win = 0; builds = set()
    for lab, j in R.POS:
        rows = []
        for s in c["seeds"][lab]:
            f = os.path.join(d0, f"x_{lab}", f"red_{s}.csv")
            if not os.path.exists(f) or os.path.getsize(f) == 0:
                miss += 1; continue
            r = pd.read_csv(edmd_acc_guard.guard(f)).iloc[0]; rows.append(r)
            uw += r["u_wall_max"] != 0.0; wd += r["W_div"] != 0.0
            win += abs(r["window"] - (T1 - T0)) > 1e-9 or r["t_last"] < T1
            health += len(HEALTH.findall(open(edmd_acc_guard.guard(os.path.join(d0, f"x_{lab}", f"run_{s}.log")), errors="ignore").read()))
            sm = pd.read_csv(edmd_acc_guard.guard(os.path.join(d0, f"x_{lab}", f"summary_{s}.csv"))).iloc[-1]; builds.add(str(sm["build_git"]))
            eta_rec = c["Ns"] * math.pi * 0.25 / (c["H"] * c["L0"])
            ok = [int(sm["wall_hold_steps"]) == HOLD, int(sm["steps_after_release"]) == POST,
                  abs(sm["L0"] - c["L0"]) < 5e-5, abs(float(sm["height"]) - c["H"]) < 5e-5,
                  int(sm["particles_total"]) == 2 * c["Ns"], abs(sm["wall_thickness_sigma"] - 0.05) < 1e-6,
                  abs(sm["box_width_sigma"] - 2 * c["L0"]) < 1e-5, abs(sm["eta_nominal"] - eta_rec) < 1e-6,
                  abs(float(sm["wall_positions_cli"]) - c["xw"][lab]) < 5e-5]
            geo += not all(ok)
        if not rows: return dict(error="no seeds")
        D = pd.DataFrame(rows); n = len(D); gL, gR = D["F_L"] / D["T_L"], D["F_R"] / D["T_R"]
        st[j] = dict(n=n, gL=gL.mean(), gR=gR.mean(), sgL=gL.std(ddof=1) / math.sqrt(n), sgR=gR.std(ddof=1) / math.sqrt(n),
                     FL=D["F_L"].mean(), FR=D["F_R"].mean(), sFL=D["F_L"].std(ddof=1) / math.sqrt(n),
                     sFR=D["F_R"].std(ddof=1) / math.sqrt(n), T=float(((D["T_L"] + D["T_R"]) / 2).mean()),
                     dT=float(np.max(np.abs(np.r_[D["T_L"] - 1, D["T_R"] - 1]))))
    Ft = {j: 0.5 * (st[j]["gL"] + st[-j]["gR"]) for j in (-2, -1, 0, 1, 2)}
    sFt = {j: 0.5 * math.hypot(st[j]["sgL"], st[-j]["sgR"]) for j in (-2, -1, 0, 1, 2)}
    dL = c["dL"]
    kT = -(Ft[-2] - 8 * Ft[-1] + 8 * Ft[1] - Ft[2]) / (12 * dL)
    s_kT = math.sqrt(sFt[-2] ** 2 + 64 * sFt[-1] ** 2 + 64 * sFt[1] ** 2 + sFt[2] ** 2) / (12 * dL)
    F0 = 0.5 * (st[0]["FL"] + st[0]["FR"]); temp = st[0]["T"]; F2 = F0 * F0 / (c["Ns"] * temp)
    with contextlib.redirect_stdout(io.StringIO()):
        dr = PH.drift(c, cell_dir=d0)
    return dict(error=None, miss=miss, geo=geo, health=health, uw=uw, wd=wd, win=win, builds=builds, st=st, kT=kT, s_kT=s_kT,
                F0=F0, temp=temp, F2=F2, static=kT + F2, f=dr["f"], s_f=dr["s_f"], dTmax=max(st[j]["dT"] for j in st),
                n=sum(st[j]["n"] for j in st))


def fit_c(y, s, Ns):
    x = 1 / np.array(Ns, float); w = 1 / np.array(s) ** 2; y = np.array(y)
    c = float((w * x * y).sum() / (w * x * x).sum()); sc = 1 / math.sqrt(float((w * x * x).sum()))
    return c, sc, float((w * (y - c * x) ** 2).sum()), float((w * y * y).sum())


def main():
    with contextlib.redirect_stdout(io.StringIO()):
        CS = R.cells()
        for c in CS:
            R.method_B(c); R.method_A(c); R.identity(c)
            d = PH.drift(c); c["kTc"], c["s_kTc"] = d["kTc"], d["s_kTc"]
            c["rho_c4"] = (c["kS"] - (d["kTc"] + c["F2term"])) / c["kS"]
            c["s_rho_c4"] = math.hypot(c["s_kS"] / c["kS"], d["s_kTc"] / c["kS"])
    if TEST:
        CS = [c for c in CS if c["cid"] == "epi8_H_H10_L10"]
        for c in CS:
            c["seeds"] = {lab: list(range(9700, 9704)) for lab, _ in R.POS}
        print("MECHANICS TEST on the A-fixed pilot cell (gate-G1 data; the anchor's seeds 9700-9703 only) -- nothing written\n")
    for c in CS:
        c["af"] = inventory_and_static(c)
    print("### G2/G3 -- inventory and drift checks per cell (A-fixed, build 279282b)\n")
    print("| eta | cell | seeds present (exp.) | missing | geometry/flags differ | health | u_wall != 0 | W != 0 | window/t_last bad "
          "| build | 1 - f (T balance) | max abs(T - 1) | G3 |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    flagged = []
    for c in CS:
        a = c["af"]
        if a.get("error"):
            print(f"| {c['lab']} | {c['cid']} | -- | all | | | | | | | | | **{a['error']}** |"); flagged.append(c["cid"]); continue
        g3 = a["uw"] == 0 and a["wd"] == 0 and abs(1 - a["f"]) < 0.002
        g2 = a["miss"] == 0 and a["geo"] == 0 and a["health"] == 0 and a["win"] == 0 and a["builds"] == {BUILD}
        if not (g2 and g3): flagged.append(c["cid"])
        exp = sum(len(v) for v in c["seeds"].values())
        print(f"| {c['lab']} | {c['cid']} | {a['n']} ({exp}) | {a['miss']} | {a['geo']} | {a['health']} | {a['uw']} | {a['wd']} | "
              f"{a['win']} | {','.join(sorted(a['builds']))} | {1 - a['f']:+.2e} | {a['dTmax']:.1e} | {'PASS' if g2 and g3 else '**FLAG**'} |")
    print(f"\nflagged cells: {flagged or 'none'}")
    if len(flagged) > 2:
        print("**STOP: more than two flagged cells (sec. 3.4) -- design failure, no P1-P3.**"); return 1
    use = [c for c in CS if c["cid"] not in flagged]
    for c in use:
        a = c["af"]; c["rho_af"] = (c["kS"] - a["static"]) / c["kS"]; c["s_rho_af"] = math.hypot(c["s_kS"] / c["kS"], a["s_kT"] / c["kS"])
        c["z_af"] = c["rho_af"] / c["s_rho_af"]
        c["s_p1"] = math.hypot(a["s_kT"], c["s_kTc"]) / c["kS"]; c["z_p1"] = (c["rho_af"] - c["rho_c4"]) / c["s_p1"]
        c["g_af"] = c["kS"] / a["kT"]; c["s_g_af"] = c["g_af"] * math.hypot(c["s_kS"] / c["kS"], a["s_kT"] / a["kT"])
    print("\n### Static stiffness with the divider held for the whole record, and the identity (sec. 3.4)\n")
    print("| eta | cell | N_s | k_T (A-fixed) | sigma | k_T (C4) | k_T (registered, released) | F(L_0) | kT | static (A-fixed) | k_S^dyn "
          "| rho_I A-fixed [%] | sigma [%] | rho/sigma | rho_I C4 [%] | P1: (AF - C4)/sigma | 2 Delta_C [%] |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for c in use:
        a = c["af"]
        print(f"| {c['lab']} | {c['cid']} | {c['Ns']} | {a['kT']:.6g} | {a['s_kT']:.2g} | {c['kTc']:.6g} | {c['kT']:.6g} | {a['F0']:.6g} | "
              f"{a['temp']:.8f} | {a['static']:.6g} | {c['kS']:.6g} | {100 * c['rho_af']:+.3f} | {100 * c['s_rho_af']:.3f} | "
              f"{c['z_af']:+.2f} | {100 * c['rho_c4']:+.3f} | {c['z_p1']:+.2f} | {200 * c['DC']:+.3f} |")
    p1_bad = [c["cid"] for c in use if abs(c["z_p1"]) > 2]; chi_p1 = sum(c["z_p1"] ** 2 for c in use)
    p2_bad = [c["cid"] for c in use if abs(c["z_af"]) > 2]; chi_p2 = sum(c["z_af"] ** 2 for c in use)
    print(f"\n**P1 (A-fixed = C4 within 2 sigma at every cell): {'HOLDS' if not p1_bad else 'FAILS'}** -- {len(use) - len(p1_bad)} of "
          f"{len(use)} cells; chi2 = {chi_p1:.1f} / {len(use)} (p = {CHI2.sf(chi_p1, len(use)):.3g})"
          f"{'' if not p1_bad else '; outside: ' + ', '.join(p1_bad)}")
    print(f"**P2 (identity, C1 rule: |rho_I| <= 2 sigma at every cell): {'HOLDS' if not p2_bad else 'FAILS'}** -- "
          f"{len(use) - len(p2_bad)} of {len(use)} cells; chi2 = {chi_p2:.1f} / {len(use)} (p = {CHI2.sf(chi_p2, len(use)):.3g})"
          f"{'' if not p2_bad else '; outside: ' + ', '.join(p2_bad)}")
    print("\n### P3 -- the 1/N_s residual: rho_I = c/N_s per density (sec. 3.6)\n")
    print("| eta | cells | c (A-fixed) | sigma_c | c/sigma_c | A_C = 2 N_s Delta_C | r = c/A_C | chi2 fit / dof | chi2 (rho = 0) / dof "
          "| r of C4 (same fit) | (r - r_C4)/sigma | (r - 1)/sigma_r |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    out3 = {}
    for lab in ("0.10", "0.39"):
        cc = [c for c in use if c["lab"] == lab]
        if len(cc) < 2: continue
        anc = [c for c in CS if c["lab"] == lab and c["scan"] == "H" and abs(c["H"] - 10) < 1e-9][0]; AC = 2 * anc["Ns"] * anc["DC"]
        cf, scf, ch, ch0 = fit_c([c["rho_af"] for c in cc], [c["s_rho_af"] for c in cc], [c["Ns"] for c in cc])
        c4, sc4, _, _ = fit_c([c["rho_c4"] for c in cc], [c["s_rho_c4"] for c in cc], [c["Ns"] for c in cc])
        r, sr, r4, sr4 = cf / AC, scf / AC, c4 / AC, sc4 / AC
        out3[lab] = dict(c=cf, sc=scf, r=r, sr=sr, r4=r4, sr4=sr4)
        print(f"| {lab} | {len(cc)} | {cf:+.4f} | {scf:.4f} | {cf / scf:+.1f} | {AC:.4f} | {r:.2f} +- {sr:.2f} | {ch:.1f} / {len(cc) - 1} | "
              f"{ch0:.1f} / {len(cc)} | {r4:.2f} +- {sr4:.2f} | {(r - r4) / math.hypot(sr, sr4):+.1f} | {(r - 1) / sr:+.1f} |")
    print("\nDeclared outcomes (sec. 3.6), evaluated:")
    if out3 and all(abs(v["c"]) < 2 * v["sc"] for v in out3.values()):
        print("- |c| < 2 sigma_c at both densities: NO 1/N_s residual -- the C4-corrected residual was an artefact of correcting a moving divider.")
    for lab, v in out3.items():
        if v["c"] > 2 * v["sc"]:
            agree = abs(v["r"] - v["r4"]) <= 2 * math.hypot(v["sr"], v["sr4"])
            print(f"- eta {lab}: c > 2 sigma_c; r {'within' if agree else 'NOT within'} 2 sigma of the C4 value"
                  f"{' -> the residual is physics (dynamic stiffness exceeds static, proportional to 1/N_s)' if agree else ''}; "
                  f"r {'within' if abs(v['r'] - 1) <= 2 * v['sr'] else 'not within'} 2 sigma of 1 (C's mechanism at its predicted size).")
        elif v["c"] < -2 * v["sc"]:
            print(f"- eta {lab}: c < -2 sigma_c -- the static stiffness exceeds the dynamic one: new.")
        elif not all(abs(w["c"]) < 2 * w["sc"] for w in out3.values()):
            print(f"- eta {lab}: |c| < 2 sigma_c (no residual at this density).")
    print("\n### gamma_box = k_S^dyn / k_T(A-fixed) (no pass/fail)\n")
    print("| eta | cell | gamma_box | sigma | bulk |\n|---|---|---|---|---|")
    for c in use:
        print(f"| {c['lab']} | {c['cid']} | {c['g_af']:.4f} | {c['s_g_af']:.4f} | {c['gbulk']:.5f} |")
    if TEST:
        return 0
    pd.DataFrame([dict(eta_lab=c["lab"], cell=c["cid"], N_s=c["Ns"], k_T_afix=c["af"]["kT"], sigma_k_T_afix=c["af"]["s_kT"],
                       F_L0=c["af"]["F0"], kT=c["af"]["temp"], static_afix=c["af"]["static"], k_S_dyn=c["kS"], sigma_k_S_dyn=c["s_kS"],
                       rho_I_afix=c["rho_af"], sigma_rho_I_afix=c["s_rho_af"], rho_I_C4=c["rho_c4"], sigma_rho_I_C4=c["s_rho_c4"],
                       rho_I_registered=c["rho"], z_P1=c["z_p1"], two_Delta_C=2 * c["DC"], one_minus_f=1 - c["af"]["f"],
                       gamma_box_afix=c["g_af"], gamma_bulk=c["gbulk"]) for c in use]).to_csv(
        os.path.join(R.OUT, "261005_p1_identity_afix_cells.csv"), index=False)
    figure(use)
    print("\ntables -> 261005_p1_identity_afix_cells.csv; figure -> 261005_p1_identity_afix.png/.pdf")
    return 0


def figure(use):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(14.5, 6.3), gridspec_kw=dict(width_ratios=[1, 1.3]))
    for lab, col, mk, nm in (("0.10", R.BLUE, "o", "eta = 0.100"), ("0.39", R.BLUE2, "s", "eta = pi/8")):
        cc = [c for c in use if c["lab"] == lab]
        a1.errorbar([c["af"]["static"] for c in cc], [c["kS"] for c in cc], xerr=[c["af"]["s_kT"] for c in cc],
                    yerr=[c["s_kS"] for c in cc], fmt=mk, color=col, ms=6, capsize=2.5, label=nm)
    lo = min(c["af"]["static"] for c in use) / 1.6; hi = max(c["af"]["static"] for c in use) * 1.6
    a1.plot([lo, hi], [lo, hi], color="k", ls="--", lw=1.2, label="identity  k_S = k_T + F^2/(N_s kT)")
    a1.set_xscale("log"); a1.set_yscale("log"); a1.set_xlim(lo, hi); a1.set_ylim(lo, hi)
    a1.set_xlabel("static, A-fixed (divider held):  k_T + F^2/(N_s kT)"); a1.set_ylabel("dynamic, method B:  k_S^dyn = M_hat omega_1^2 / 2")
    a1.grid(True, which="both", ls=":", alpha=0.5); a1.legend(fontsize=8.5, loc="upper left")
    a1.set_title("Identity with the divider held for the whole record (261012 sec. 3)", fontsize=10.5)
    order = sorted(use, key=lambda c: (c["lab"], c["scan"], c["Ns"], c["L0"]))
    for i, c in enumerate(order):
        col = R.BLUE if c["lab"] == "0.10" else R.BLUE2; mk = "o" if c["lab"] == "0.10" else "s"
        a2.errorbar(i - 0.2, 100 * c["rho_c4"], yerr=100 * c["s_rho_c4"], fmt=mk, mfc="white", color=col, ms=5, capsize=2, lw=0.9, alpha=0.8)
        a2.errorbar(i + 0.1, 100 * c["rho_af"], yerr=100 * c["s_rho_af"], fmt=mk, color=col, ms=6.5, capsize=3, lw=1.4)
        a2.plot(i, 200 * c["DC"], marker="_", color="#8f4fd1", ms=14, mew=2)
    a2.axhline(0, color="k", ls="--", lw=1)
    a2.plot([], [], "o", color=R.BLUE, label="A-fixed (registered estimator, sec. 3.4)")
    a2.plot([], [], "o", mfc="white", color=R.BLUE, label="C4 (post-hoc drift correction, sec. 2.7)")
    a2.plot([], [], marker="_", color="#8f4fd1", ls="none", ms=14, mew=2, label="hypothesis C: rho_I = 2 Delta_C")
    a2.set_xticks(range(len(order))); a2.set_xticklabels([c["cid"].replace("e0p10_", "0.10 ").replace("epi8_", "pi/8 ") for c in order],
                                                         rotation=70, ha="right", fontsize=7)
    a2.set_ylabel("rho_I = (k_S^dyn - static)/k_S^dyn  [%]"); a2.grid(True, ls=":", alpha=0.5); a2.legend(fontsize=8, loc="upper left")
    a2.set_title("Residual per cell (1 sigma)", fontsize=10.5)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(R.OUT, f"261005_p1_identity_afix.{ext}"), dpi=200)


if __name__ == "__main__":
    sys.exit(main())
