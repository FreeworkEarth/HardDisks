#!/usr/bin/env python3
"""##CHRIS 2026-10-07 (261012 sec. 4.4.10 part E, gate version 2): the numbers of TEST T, computed before any of its data.
Test T: cell epi8_H_H10_L10, method B, the campaign protocol; 100 FRESH seeds per mass (400 at M = 50: amendment of the
second plan-author decision of 2026-10-07, item 2, before any data), the same seeds on both policies
(minimal, --legacy-resched), same binary, same partition, tasks interleaved. Eleven numbers, minimal minus legacy, with the
registered estimators: k_S^dyn, c_s, and the mean nu of each of the nine masses; z with both SEs.
RULE (plan author, fixed 2026-10-07): PASS if all eleven |z| < z* (two-sided Bonferroni, family-wise false-fail 5 %, n = 11)
AND the permutation p of the nine-mass chi2 >= 0.01; otherwise FAIL.
Printed here: z*; the expected sigma_diff of each number with 100 + 100 seeds (from the per-seed spread of the 279282b campaign
and the 73fc07f replay of this cell, scaled by 1/sqrt(100)); the z and the power P(|z| >= z*) of each hypothesis the replay
generated, at its observed size; the fresh seed list (run_seed(BASE_T, 0, mass index, r), r = 0..99) checked against every
seed of the confinement campaign's task lists and the smoke/pilot seeds, and its SHA-256; the cost from the measured KOA times
(minimal: the replay's run logs; legacy: the 279282b campaign's run logs, whose speed the legacy path reproduces, sec. 4.4.9).
usage (from hspist3/): python3 validation/resched_testT_design_261007.py
"""
import contextlib, glob, hashlib, io, math, os, re, sys
import numpy as np
from scipy.stats import norm, ncx2, chi2 as CHI2
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, HS)
import paper1_confinement_results_261004 as R
import tests_20260913 as T

CID, BASE_T, ALPHA_FW, NNUM = "epi8_H_H10_L10", 20261007, 0.05, 11
NSEED = {M: (400 if M == 50 else 100) for M in T.A1_MASSES}   # amended 2026-10-07 (item 2): M = 50 gets 400 per policy
NEW = os.path.join(HS, "experiments_resched_gate_261005")
CONF = os.path.join(HS, "cluster", "confinement_20261013")


def main():
    zs = norm.isf(ALPHA_FW / (2 * NNUM))
    print("# Test T -- design numbers (gate version 2, 261012 sec. 4.4.10 E3), before any of its data\n")
    print(f"z* = Phi^-1(1 - {ALPHA_FW}/(2 x {NNUM})) = {zs:.4f}   (two-sided Bonferroni, family-wise false-fail {ALPHA_FW:.0%}, n = {NNUM})")
    print(f"chi2 part: PASS needs the permutation p of the nine-mass chi2 >= 0.01 (nominal 0.99 quantile of chi2_9 = {CHI2.isf(0.01, 9):.2f})\n")
    with contextlib.redirect_stdout(io.StringIO()):
        c = {x["cid"]: x for x in R.cells()}[CID]
        co = dict(c); R.method_B(co)
        old = R.HS; R.HS = NEW; R.ALLOW_NEW_BUILD = True
        cn = dict(c); R.method_B(cn)
        R.HS = old; R.ALLOW_NEW_BUILD = False
    sd = {}
    for ro, rn in zip(co["B"], cn["B"]):
        sd[ro["M"]] = math.sqrt((np.var(ro["nus"], ddof=1) + np.var(rn["nus"], ddof=1)) / 2)   # pooled per-seed SD
    print("### Expected sigma_diff with the amended seed numbers (n + n per mass), and the power of the rule for the replay's hypotheses\n")
    print("| number | seeds per policy | hypothesis (replay, observed size) | shift | sigma_diff | expected z | P(abs(z) >= z*) |\n|---|---|---|---|---|---|---|")
    hyp = {500: +0.0062, 1500: +0.0057, 50: -0.0105}
    se_t = {M: sd[M] / math.sqrt(NSEED[M]) for M in sd}                 # per-policy SE of the mean nu in Test T
    for ro in co["B"]:
        M = ro["M"]; s = math.sqrt(2.0) * se_t[M]; d = hyp.get(M, 0.0) * ro["nu"]; z = d / s
        p = norm.sf(zs - z) + norm.cdf(-zs - z)
        lab = f"alpha = {M / 100:g} {100 * hyp[M]:+.2f} %" if M in hyp else "none (null)"
        print(f"| mean nu, M = {M} | {NSEED[M]} | {lab} | {d:+.3e} | {s:.3e} | {z:+.2f} | {p:.3f} |")
    # c_s and k_S^dyn: their registered SE formulas with the per-mass Test T SEs (both policies alike), sigma_diff = sqrt(2) SE
    Ms = [ro["M"] for ro in co["B"]]; xs = np.asarray(co["x"], float); nus = np.array([ro["nu"] for ro in co["B"]])
    ses = np.array([se_t[M] for M in Ms])
    s_cs = math.sqrt(2.0) * math.sqrt(float((xs * xs * ses * ses).sum())) / float((xs * xs).sum())
    heavy = np.array([any(abs(M / (2.0 * c["Ns"]) - a) < 1e-9 for a in R.HEAVY) for M in Ms])
    Mh = np.array([M + 2.0 * c["Ns"] / 3.0 for M in Ms]); k = Mh * (2 * math.pi * nus) ** 2 / 2.0; sk = k * 2 * ses / nus
    s_ks = math.sqrt(2.0) / math.sqrt(float((1.0 / sk[heavy] ** 2).sum()))
    for lab, d, s in (("k_S^dyn", 0.0356, s_ks), ("c_s (replay shift, -0.28 %)", -0.0106, s_cs)):
        z = d / s; p = norm.sf(zs - z) + norm.cdf(-zs - z)
        print(f"| {lab} | - | observed in the replay | {d:+.4g} | {s:.4g} | {z:+.2f} | {p:.3f} |")
    lam = sum(((hyp.get(ro["M"], 0.0) * ro["nu"]) / (math.sqrt(2.0) * se_t[ro["M"]])) ** 2 for ro in co["B"])
    print(f"\nchi2 part under the three per-mass hypotheses together: noncentrality {lam:.1f}; P(chi2_9 >= {CHI2.isf(0.01, 9):.2f}) = "
          f"{ncx2.sf(CHI2.isf(0.01, 9), 9, lam):.3f} (nominal quantile; the permutation threshold is printed by the Test T analysis)")
    print("\n### Fresh seeds\n")
    seeds = {(mi, r): T.run_seed(BASE_T, 0, mi, r) for mi, M in enumerate(T.A1_MASSES) for r in range(NSEED[M])}
    used = set()
    for f in glob.glob(os.path.join(CONF, "tasks_*.txt")):
        for l in open(f):
            p = l.split()
            if p[0] == "B": used.add(int(p[4]))
            elif p[0] in ("A", "AF"): used.add(int(p[3]))
    pilot = {T.run_seed(20261013, 0, mi, r) for mi in range(9) for r in range(4)}
    vals = list(seeds.values())
    print(f"seeds: run_seed({BASE_T}, 0, mass index, r), r = 0..n_M - 1 (n = 400 at M = 50, 100 otherwise): {len(vals)} seeds, {len(set(vals))} distinct")
    print(f"overlap with every seed of the campaign task lists ({len(used)} seeds, B/A/AF): {len(set(vals) & used)}; with the smoke/pilot seeds: {len(set(vals) & pilot)}")
    txt = "".join(f"{M} {r} {seeds[(mi, r)]}\n" for mi, M in enumerate(T.A1_MASSES) for r in range(NSEED[M]))
    print(f"SHA-256 of the seed list (lines 'M r seed', mass ascending, r ascending): {hashlib.sha256(txt.encode()).hexdigest()}")
    print("\n### Cost from measured KOA times (seconds per trajectory, mean over the 25 seeds)\n")
    print("| M | seeds per policy | minimal (replay, 73fc07f) [s] | legacy (279282b campaign) [s] | n + n trajectories [core-h] | "
          "array task wall time on 8 cores [min] |\n|---|---|---|---|---|---|")
    RB = os.path.join("experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/confinement_B_20261013", CID)
    tot = 0.0; walls = {}
    for M in T.A1_MASSES:
        tm = [int(x) for x in re.findall(r"##RUN .*?\((\d+) s\)", open(os.path.join(NEW, RB, f"m_{M}", "run.log"), errors="ignore").read())]
        tl = [int(x) for x in re.findall(r"##RUN .*?\((\d+) s\)", open(os.path.join(HS, RB, f"m_{M}", "run.log"), errors="ignore").read())]
        ch = NSEED[M] * (np.mean(tm) + np.mean(tl)) / 3600; tot += ch; wall = ch * 60 / 8; walls[M] = wall
        print(f"| {M} | {NSEED[M]} | {np.mean(tm):.1f} | {np.mean(tl):.1f} | {ch:.2f} | {wall:.1f} |")
    ntr = 2 * sum(NSEED.values())
    print(f"\ntotal: {tot:.1f} core-h for {ntr} trajectories")
    print(f"--time of the array tasks: 2 x the longest measured task ({max(walls.values()):.1f} min, M = {max(walls, key=walls.get)}) "
          f"= {2 * max(walls.values()):.1f} min -> 0:30 for every task; the M = 50 task: {walls[50]:.1f} min measured, 2 x = {2 * walls[50]:.1f} min")


if __name__ == "__main__":
    main()
