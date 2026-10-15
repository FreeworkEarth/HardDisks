#!/usr/bin/env python3
"""##CHRIS 2026-10-12 (plan-author decision 39 items 5 d and 5 f; registration 261012 sec. 4.8.10, GATE G3): THE ESTIMATOR GATE. The stage 1
reducer, applied to existing data, must reproduce to the printed digits:
  (i)  Test G (Mac), gen3 arm: nu_d-bar = 0.0370095 (M = 300) and 0.0172665 (M = 1500) (notes sec. 4.7.15 addendum; printed by
       validation/testG_verdict_261009.information_nud);
  (ii) Test T legacy, the table of methods sec. 16 (printed by validation/resched_testT_followup_261007.py, item 5): Q per mass (the
       item-4a model, averaged ACF), the nu_0 deficit column, and the nine-mass weighted nu_0 slope 3.81286. "The reducer has a nine-mass
       mode for this gate only": item 5's arithmetic, its weights from the nu_d seed SE. The E1 value (eight masses, the jackknife SE) of
       that cell is printed for the record.
THE STAGE 1 REDUCER is applied literally: the two inputs of each mass folder (red_nu.csv, acf_runs.npz) are COPIED to a new temporary
folder and cluster/gen3_testR_261011/reduce_testR_koa.nud_table (nu_d per trajectory, nu_ref = the folder's own argmax mean: decision 39
item 5 d) and cluster/gen3_stage1_261012/reduce_s1_koa.acf_mean_fit (the seed-averaged ACF's item-4a fit) run there; nothing is written
into the recorded folders. The printed values were computed with nu_ref = the mean of TWO arms' argmax means (Test G: gen3 and gen2;
Test T: minimal and legacy), so the start value of the fits differs; this gate shows whether that changes a printed digit. If it does:
the difference is printed and the gate STOPS (decision 39 item 5 d).
usage (from hspist3/): python3 cluster/gen3_stage1_261012/g3_estimator_gate_261012.py [--tmp DIR]
"""
import contextlib, io, math, os, shutil, sys, tempfile
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE)); GEN3 = os.path.dirname(HS)
MAIN_HS = os.path.join(os.path.dirname(GEN3), "HardDisks", "hspist3")
os.environ["HD_RESCHED2_LOC"] = os.path.join(MAIN_HS, "experiments_resched_gate2_261007")      # Test T's data are in the main worktree
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HS, "cluster", "gen3_testR_261011")); sys.path.insert(0, os.path.join(HS, "validation")); sys.path.insert(0, HS)
import reduce_testR_koa as RT
import reduce_s1_koa as RS
import resched_testT_followup_261007 as FU
import s1_verdict_261012 as V
import tests_20260913 as T
TG = os.path.join(MAIN_HS, "experiments_gen3_gate_261009", "testG", "gen3", "epi8_H_H10_L10")
TARGET_G = {300: "0.0370095", 1500: "0.0172665"}
TARGET_Q = {50: "9.3", 100: "13.5", 200: "20.3", 300: "25.5", 500: "32.2", 750: "40.3", 1000: "43.7", 1500: "54.9", 2000: "64.9"}
TARGET_D0 = {50: "-0.424 +- 0.078", 100: "+0.052 +- 0.160", 200: "+0.169 +- 0.100", 300: "-0.109 +- 0.077", 500: "+0.097 +- 0.073",
             750: "-0.025 +- 0.056", 1000: "+0.047 +- 0.050", 1500: "+0.100 +- 0.049", 2000: "-0.002 +- 0.043"}
TARGET_SLOPE = "3.81286"


def reduce_copy(src, dst):
    os.makedirs(dst); [shutil.copy2(os.path.join(src, f), dst) for f in ("red_nu.csv", "acf_runs.npz")]
    RT.nud_table(dst); RS.acf_mean_fit(dst)
    return pd.read_csv(os.path.join(dst, "nud_runs.csv")), pd.read_csv(os.path.join(dst, "acf_mean_fit.csv")).iloc[0]


def main():
    tmp = sys.argv[sys.argv.index("--tmp") + 1] if "--tmp" in sys.argv else tempfile.mkdtemp(prefix="g3_gate_")
    print("# Gate G3, the estimator gate (sec. 4.8.10), printed by cluster/gen3_stage1_261012/g3_estimator_gate_261012.py\n")
    print(f"the stage 1 reducer on copies under {tmp}; nu_ref = each folder's own argmax mean (decision 39 item 5 d)\n")
    ok_all = []
    print("## (i) Test G (Mac), gen3 arm\n\n| M | trajectories | fit failures | nu_ref (own) | nu_d-bar | printed | reproduced to the printed digits | difference |\n|---|---|---|---|---|---|---|---|")
    for M, tgt in TARGET_G.items():
        nd, mf = reduce_copy(os.path.join(TG, f"m_{M}"), os.path.join(tmp, "testG", f"m_{M}"))
        v = nd.loc[nd["fit_ok"] == 1, "nu_d"].mean(); s = f"{v:.7f}"; ok = s == tgt; ok_all.append(ok)
        print(f"| {M} | {len(nd)} | {int((nd['fit_ok'] == 0).sum())} | {nd['nu_ref'].iloc[0]:.7f} | {v:.9f} | {tgt} | {'yes' if ok else '**NO**'} | "
              f"{v - float(tgt):+.2e} ({100 * (v / float(tgt) - 1):+.5f} %) |")
    # Test T's cell epi8_H_H10_L10 (N 100, H = L0 = 10): the nine A1 v2 masses, N_s = 50, x_M = K(M / 2 N_s) / (2 pi L_eff), L_eff = L0 - 1.025,
    # as resched_testT_followup_261007.load() takes them from paper1_confinement_results' method_B (whose loader reads campaign data
    # that exist only in the main worktree)
    Ms, Ns = list(T.A1_MASSES), 50; x = np.array([T.x_of(M, 10.0) for M in Ms])
    print(f"\n## (ii) Test T legacy, methods sec. 16 (the nine-mass mode: item 5's arithmetic)\n\ncell {FU.CID} (N_s {Ns}); masses {Ms}\n")
    rows = {}
    for M in Ms:
        src = os.path.join(FU.TT.LOC, FU.TT.REL_T, "legacy", FU.CID, f"m_{M}")
        nd, mf = reduce_copy(src, os.path.join(tmp, "testT_legacy", f"m_{M}"))
        v = nd.loc[nd["fit_ok"] == 1, "nu_d"].to_numpy(float)
        nu0 = math.sqrt(v.mean() ** 2 + (1 / (2 * math.pi * mf["tau_r"])) ** 2)
        rows[M] = dict(nd=v, nud=v.mean(), tau_r=mf["tau_r"], Q_meanfit=mf["Q_meanfit"], Q_reg=math.pi * v.mean() * mf["tau_r"], nu0=nu0,
                       fails=int((nd["fit_ok"] == 0).sum()), n=len(nd))
    y = np.array([rows[M]["nu0"] for M in Ms]); sy = np.array([rows[M]["nd"].std(ddof=1) / math.sqrt(len(rows[M]["nd"])) for M in Ms])
    wt = 1 / sy ** 2; sw = float((wt * x * y).sum() / (wt * x * x).sum()); d0 = y / x / sw - 1; s0 = sy / x / sw
    print("| M | trajectories | fit failures | Q (item-4a model, averaged ACF) | printed | equal | deficit nu_0 [%] | printed | equal | Q as E1 defines it (pi nu_d-bar tau_r) |\n"
          "|---|---|---|---|---|---|---|---|---|---|")
    for j, M in enumerate(Ms):
        q = f"{rows[M]['Q_meanfit']:.1f}"; dd = f"{100 * d0[j]:+.3f} +- {100 * s0[j]:.3f}"
        e1, e2 = q == TARGET_Q[M], dd == TARGET_D0[M]; ok_all += [e1, e2]
        print(f"| {M} | {rows[M]['n']} | {rows[M]['fails']} | {q} | {TARGET_Q[M]} | {'yes' if e1 else '**NO**'} | {dd} | {TARGET_D0[M]} | {'yes' if e2 else '**NO**'} | {rows[M]['Q_reg']:.1f} |")
    s = f"{sw:.5f}"; ok = s == TARGET_SLOPE; ok_all.append(ok)
    print(f"\nthe nine-mass weighted nu_0 slope: {sw:.7f} -> {s}; printed {TARGET_SLOPE}: {'yes' if ok else '**NO**'} (difference {sw - float(TARGET_SLOPE):+.2e})")
    cell = dict(cid="testT_legacy", g="-", blk="B", N=100, eta=math.pi / 8, start="default", H=10.0, L0=10.0, Ns=50, Teq=0.0,
                masses={M: M / 100.0 for M in Ms}, nseed=max(r["n"] for r in rows.values()))
    prow = []
    for M in Ms:
        d = os.path.join(tmp, "testT_legacy", f"m_{M}")
        Dm = dict(red=pd.read_csv(os.path.join(d, "red_nu.csv")), nud=pd.read_csv(os.path.join(d, "nud_runs.csv")), health=None,
                  failed=pd.DataFrame(columns=["folder", "run"]), s1=None, snap=None, acf=dict(np.load(os.path.join(d, "acf_runs.npz"))))
        prow.append(V.per_mass(dict(cell, nseed=rows[M]["n"]), M, Dm))      # Test T's seeds per mass differ (400 at M = 50, 100 else)
    E = V.e1(prow)
    print(f"for the record, E1 of that cell (eight masses alpha >= 1, the delete-one jackknife SE): c_s = {E['cs']:.5f} +- {E['sig']:.5f} "
          f"(chi2 {E['chi2']:.1f}, {E['dof']} dof); alpha = 0.5: c = {prow[0]['c']:.5f} +- {prow[0]['sig_c']:.5f}" if E else "for the record, E1: none")
    print(f"\nG3: {'PASS -- every printed digit reproduced' if all(ok_all) else 'NOT REPRODUCED: the differences are printed above; STOP (decision 39 item 5 d)'} "
          f"({sum(ok_all)} of {len(ok_all)} values)")


if __name__ == "__main__":
    main()
