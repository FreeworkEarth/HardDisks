#!/usr/bin/env python3
"""##CHRIS 2026-10-14: box-truncation correction of the canonical A1 v2 table (methods sec. 14, pre-registered 4db8c9d).

SIM_WIDTH = (int)(2 * L0_UNITS * PIXELS_PER_SIGMA) (00ALLINONE.c:323) truncates the box; physics uses
boxW = XW2 - XW1 (15882). delta = box shortfall [sigma]; per compartment the mean shortfall is delta/2 (the divider
starts L_0 from the left wall and oscillates about the truncated centre -- printed below from each run's own
Center_X and Displacement). Correction: L_0,true = L_0 - delta/2, eta_true = eta_rec L_0/L_0,true,
L_eff,true = L_eff,rec - delta/2; c_s re-derived from the per-mass frequencies (canonical cell()) with L_eff,true.
Estimator gate: the recomputation must reproduce the 260919 table to 5e-6 per cell, else that cell is VOID.
Verdict: REGENERATE if any cell's D = (c_s - c_s^KR)/sigma changes by more than 0.5; KEEP otherwise.
"""
import csv, glob, math, os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import tests_20260913 as T
import edmd_acc_guard   # ##CHRIS 2026-10-08 (261012 sec. 4.7.4, decision 2): the loader provenance guard (full name: no alias can be shadowed)
import plot_speed_of_sound_edmd as sos
from paper1_populate_cs_err_20261002 import cell, slope_with_errors
REPO = os.path.dirname(os.path.dirname(HERE))
FIG = os.path.join(REPO, "0000_PLAN_OVERALL", "paper1_speedofsound", "experiments", "final")
XW1, PPS, R = 200.0, 24.0, 0.5
CANON = "260919_A1v2_final_cs_vs_eta.csv"
PRE = "260919_A1v2_final_cs_vs_eta_pre_boxtrunc_20261014.csv"   # the uncorrected table, kept as a dated copy (A2)

def source():
    """The UNcorrected table: the dated copy once it exists (after the canonical one is regenerated), else the canonical."""
    return PRE if os.path.exists(os.path.join(T.PLOTS, PRE)) else CANON

def kr(e):
    a = np.array([e]); return float(sos.cs_adiabatic_2d_monatomic(sos.Z_kolafa_rottner_2006(a), sos.dZ_kolafa_rottner_2006(a), a, kbt=1.0, m=1.0)[0])

def leaf(eta):
    best = None
    for d in os.listdir(T.DROOT):
        v = float(d[4:].replace("p", "."))
        if best is None or abs(v - eta) < abs(best[0] - eta): best = (v, d)
    return os.path.join(T.DROOT, best[1]) if abs(best[0] - eta) < 1e-3 else None

def slope_x(cs, Leff):
    x = np.array([T.k_root(q["M"] / (2.0 * T.N_SIDE)) / (2 * math.pi * Leff) for q in cs])
    y = np.array([q["nu"] for q in cs]); sy = np.array([(q["sd"] / math.sqrt(q["n"])) if q["n"] > 1 else np.nan for q in cs])
    return slope_with_errors(x, y, sy)

def compute():
    rows = list(csv.DictReader(open(os.path.join(T.PLOTS, source()))))
    out = []
    for r in rows:
        eta, L0t = float(r["eta"]), float(r["L0"]); d = leaf(eta)
        cs = []
        for M in T.A1_MASSES:
            c = cell((eta, L0t, M, T.cell_runs(os.path.join(d, f"m_{M}"), M)))
            if c["n"] > 0: cs.append(c)
        tr0 = sorted(glob.glob(os.path.join(d, "m_500", "wall_x_positions_*_run*.csv")))
        h = pd.read_csv(edmd_acc_guard.guard(tr0[0]), nrows=1).iloc[0]
        L0 = float(h["L0"]); N = int(h["Left_Count"]) + int(h["Right_Count"]); erec = float(h["eta"])
        w = np.float32(2) * np.float32(L0) * np.float32(PPS)                    # the binary's float expression
        delta = (float(w) - math.floor(float(w))) / PPS
        cen_off = float(h["Center_X(σ)"]) - (XW1 / PPS + L0)                       # recorded: should be -delta/2
        disp = np.mean([pd.read_csv(edmd_acc_guard.guard(p), usecols=["Displacement(σ)"])["Displacement(σ)"].mean() for p in tr0])
        H = N * math.pi * R * R / (2 * L0 * erec)
        Le = T.l_eff(L0t); s, e, es, ch = slope_x(cs, Le)
        ok = abs(s - float(r["c_s"])) <= 5e-6 * abs(float(r["c_s"]))
        L0T = L0t - delta / 2; eT = eta * L0t / L0T; LeT = Le - delta / 2
        sT, eT_, esT, _ = slope_x(cs, LeT)
        sB, _, esB, _ = slope_x(cs, Le - delta)                                     # plan's L_0 - delta, bound only
        eB = eta * L0t / (L0t - delta)
        sig = float(r["c_s_err_scaled"]); kR, kT, kB = kr(eta), kr(eT), kr(eB)
        if not r["KR"].strip():                       # the table tabulates no KR here: no published deviation
            kR = kT = kB = float("nan")
        Db = (float(r["c_s"]) - kR) / sig; Da = (sT - kT) / (sig * LeT / Le); DB = (sB - kB) / (sig * (Le - delta) / Le)
        out.append(dict(eta=eta, L0=L0t, delta=delta, cen=cen_off, disp=disp, H=H, N=N, eT=eT, Le=Le, LeT=LeT,
                        cs=float(r["c_s"]), sig=sig, s=s, ok=ok, sT=sT, sigT=sig * LeT / Le, kR=kR, kT=kT, tabKR=float(r["KR"]) if r["KR"].strip() else float("nan"),
                        Db=Db, Da=Da, ch=Da - Db, chB=DB - Db, ident=sT / s - LeT / Le, nm=len(cs)))
    return out

def main():
    out = compute()
    print("### Box-truncation correction, every canonical A1 v2 cell (sigma = c_s_err_scaled)\n")
    print("| eta_rec | L_0 | delta | Center_X - (XW1 + L_0) | <Displacement> (m_500) | H from eta_rec | eta_true | L_eff,rec | L_eff,true | "
          "c_s,rec ± σ | gate | c_s,true ± σ | KR(eta_rec) | KR(eta_true) | D before | D after | change | flag | change if L_0 - delta |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for o in out:
        flag = "**> 0.5**" if abs(o["ch"]) > 0.5 else ("no KR in table" if o["kR"] != o["kR"] else "")
        print(f"| {o['eta']:.6f} | {o['L0']} | {o['delta']:.6f} | {o['cen']:+.6f} | {o['disp']:+.4f} | {o['H']:.5f} | {o['eT']:.6f} | "
              f"{o['Le']:.4f} | {o['LeT']:.4f} | {o['cs']:.5f} ± {o['sig']:.5f} | {'ok' if o['ok'] else 'VOID'} | {o['sT']:.5f} ± {o['sigT']:.5f} | "
              f"{o['kR']:.5f} | {o['kT']:.5f} | {o['Db']:+.2f} | {o['Da']:+.2f} | {o['ch']:+.2f} | {flag} | {o['chB']:+.2f} |")
    pi8 = [o for o in out if abs(o["eta"] - 0.392699) < 1e-6][0]
    print(f"\npi/8 anchor: delta = {pi8['delta']:.6f} (must be exactly 0) -> {'OK' if pi8['delta'] == 0 else 'NOT ZERO'}")
    print(f"estimator gate: {sum(o['ok'] for o in out)}/{len(out)} cells reproduce the 260919 c_s to 5e-6; "
          f"KR function reproduces the table's KR column to {max(abs(o['kR'] - o['tabKR']) for o in out if o['kR'] == o['kR']):.1e} "
          f"({sum(o['kR'] == o['kR'] for o in out)} cells with a tabulated KR; the rest carry no published deviation and do not enter the verdict)")
    print(f"identity c_s,true/c_s,rec = L_eff,true/L_eff,rec: max deviation {max(abs(o['ident']) for o in out):.1e}")
    print(f"recorded Center_X offset vs -delta/2: max |difference| {max(abs(o['cen'] + o['delta'] / 2) for o in out):.1e} sigma")
    big = [o for o in out if o["ch"] == o["ch"] and abs(o["ch"]) > 0.5]
    print(f"\ncells with |change| > 0.5: {len(big)} -> **VERDICT: {'REGENERATE' if big else 'KEEP'}**" +
          (f" (eta = {', '.join(f'{o[chr(101)+chr(116)+chr(97)]:.3f}' for o in big)})" if big else ""))
    val = [o for o in out if o["ch"] == o["ch"]]
    print(f"largest |change|: {max(abs(o['ch']) for o in val):.2f} at eta = {max(val, key=lambda o: abs(o['ch']))['eta']:.4f}; "
          f"for eta <= 0.39: {max(abs(o['ch']) for o in val if o['eta'] <= 0.4):.2f}")
    figure(out)

def figure(out):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    out = [o for o in out if o["Db"] == o["Db"]]
    fig, ax = plt.subplots(figsize=(7.2, 4.2)); e = np.array([o["eta"] for o in out])
    ax.errorbar(e, [o["Db"] for o in out], yerr=1.0, fmt="o", color="0.55", mfc="white", ms=5, capsize=2,
                label="before: recorded $\\eta$, $L_{\\rm eff}$ (260919 table)")
    ax.errorbar(e * 1.01, [o["Da"] for o in out], yerr=1.0, fmt="o", color="tab:blue", ms=5, capsize=2,
                label="after: $\\eta_{\\rm true}$, $L_{\\rm eff,true}$ ($\\delta/2$ per compartment)")
    ax.axhline(0, color="red", lw=1.5, label="Kolafa–Rottner")
    ax.set_xscale("log"); ax.set_xlabel(r"$\eta$"); ax.set_ylabel(r"$(c_s - c_s^{\rm KR})/\sigma$  (error bar = 1$\sigma$)")
    ax.set_title("Box-truncation correction of the canonical A1 v2 table (after offset ×1.01 in η)", fontsize=10)
    ax.legend(fontsize=8, frameon=False); fig.tight_layout()
    for ext in ("png", "pdf"): fig.savefig(os.path.join(FIG, f"261014_p1_boxtrunc_shift.{ext}"), dpi=200)
    print("\nfigure: 0000_PLAN_OVERALL/paper1_speedofsound/experiments/final/261014_p1_boxtrunc_shift.{png,pdf}")

def table():
    """methods sec. 14.2: the full per-cell regenerated table, the box-height check, and the A2 magnitude (OPEN)."""
    out = compute()
    print(f"### Regenerated per-cell table (uncorrected input: {source()}; sigma = c_s_err_scaled, rescaled with L_eff)\n")
    print("| eta_rec | L_0 | delta | eta_true | L_eff,rec | L_eff,true | c_s,rec ± σ | c_s,true ± σ | KR(eta_rec) | KR(eta_true) | D before [σ] | D after [σ] |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for o in out:
        fm = lambda v, f: ("n/a" if v != v else format(v, f))
        print(f"| {o['eta']:.6f} | {o['L0']} | {o['delta']:.6f} | {o['eT']:.6f} | {o['Le']:.4f} | {o['LeT']:.4f} | "
              f"{o['cs']:.5f} ± {o['sig']:.5f} | {o['sT']:.5f} ± {o['sigT']:.5f} | {fm(o['kR'], '.5f')} | {fm(o['kT'], '.5f')} | "
              f"{fm(o['Db'], '+.2f')} | {fm(o['Da'], '+.2f')} |")
    print(f"\nestimator gate: {sum(o['ok'] for o in out)}/{len(out)} cells reproduce the uncorrected table to 5e-6")
    print("\n**Box HEIGHT (00ALLINONE.c:324, `SIM_HEIGHT = (int)(HEIGHT_UNITS * PIXELS_PER_SIGMA);`).** Every A1 v2 run was launched "
          "by the harness with `--height=10.0` (tests_20260913.py:78 `H = \"10.0\"`, passed at :285 as `f\"--height={H}\"`):")
    hs = sorted({float(T.H) * PPS for _ in out})
    print(f"H x 24 = {hs} -> integer in all {len(out)} cells, so SIM_HEIGHT is exact and the height is NOT truncated. "
          f"Read back from each run's own eta_rec: H = {min(o['H'] for o in out):.5f} ... {max(o['H'] for o in out):.5f} "
          f"(H x 24 = {24*min(o['H'] for o in out):.3f} ... {24*max(o['H'] for o in out):.3f}; 6-decimal eta print).")
    print("\n**OPEN, not corrected in this batch: A2 (the finite-size ladder) uses non-grid L_0 as well.** Per (eta, N), from the A2 "
          "per-mass tables the draft's zoom and overlay figures read:\n")
    print("| table | eta | N | L_0 | delta [σ] | delta/2 / L_eff (c_s shift) | eta shift |")
    print("|---|---|---|---|---|---|---|")
    worst = 0.0
    for fn in ("260919_A2_cs_per_mass.csv", "260919_A2_cs_per_mass_famB.csv"):
        seen = {}
        for r in csv.DictReader(open(os.path.join(T.PLOTS, fn))):
            seen[(float(r["eta"]), int(r["N"]))] = r["L0"]
        for (e, n), l0s in sorted(seen.items()):
            L0 = float(l0s); w = np.float32(2) * np.float32(L0) * np.float32(PPS); dl = (float(w) - math.floor(float(w))) / PPS
            rel = dl / 2 / T.l_eff(L0); worst = max(worst, rel)
            print(f"| {fn.replace('260919_A2_cs_per_mass', 'A2').replace('.csv', '')} | {e:g} | {n} | {l0s} | {dl:.6f} | {100*rel:.4f} % | {100*(L0/(L0-dl/2)-1):.4f} % |")
    print(f"\nlargest A2 c_s shift: {100*worst:.4f} %. The zoom and N100-vs-A2 figures therefore pair a corrected A1 v2 curve with "
          "uncorrected A2 points; their titles say so.")

if __name__ == "__main__":
    table() if "--table" in sys.argv else main()
