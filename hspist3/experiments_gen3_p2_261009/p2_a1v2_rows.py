#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.19; stage F): A1 v2's side of pilot P2's table, at P2's own masses. Written and committed
before P2's table was first printed. Why: sec. 4.7.11's table 2 pools A1 v2's nine masses (50 ... 2000), while P2 runs M = 50, 300,
2000 only, and the record of a trajectory is 200 of its own periods, so a heavy divider's record is longer in sigma-time. A like-for-
like comparison needs A1 v2 at the same masses. Also printed: A1 v2's psi6 at the end of its hold (2000 steps = 33.3 sigma-time) and
at the end of its record, which sec. 4.7.11's table does not print.
Existing data only, nothing is run. The cells, trajectories and estimator are sec. 4.7.11's own (validation/paper1_window_aging_261009:
paper1_window_explore_261008's discovery and filters, its bad-run list, the trace check, and traj2; per mass only the trajectories
with a per-run psi6 record). GATE: the "all nine" row of each eta must equal sec. 4.7.11's printed table 2 (trajectories, mean dnu
(SE), dnu exactly 0, mean dA (SE)) to the printed digits; read from its recorded output.
usage (from hspist3/ of the engine-gen3 worktree): python3 experiments_gen3_p2_261009/p2_a1v2_rows.py [--workers 6] > p2_a1v2_rows_output.txt
"""
import math, os, re, sys
from collections import defaultdict
from multiprocessing import Pool
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import p2_run as P2
sys.path.insert(0, os.path.join(P2.MAIN_HS, "validation")); sys.path.insert(0, P2.MAIN_HS)
import tests_20260913 as T
import paper1_window_explore_261008 as W
import paper1_window_aging_261009 as WA
from paper1_populate_cs_err_20261002 import box_delta

REC = os.path.join(T.PAPER1, "exploratory_261009_followup", "261009_window_aging_amplitude_output.txt")   # sec. 4.7.11's output
HOLD_SIGMA = 2000 * 0.4 / 24.0     # A1 v2: 'Wall hold steps: 2000', 'Fixed dt override: 4.00e-01' (run.log); 24 units = 1 sigma-time


def mse(x):
    x = np.asarray(x, float)
    return (x.mean(), x.std(ddof=1) / math.sqrt(len(x))) if len(x) > 1 else (math.nan, math.nan)


def recorded():
    """sec. 4.7.11 table 2, A1 v2 rows: eta -> (trajectories, 'mean dnu (SE)', 'n0 (pct)', 'mean dA (SE)') as printed"""
    out = {}; on = False
    for l in open(REC, errors="ignore"):
        if l.startswith("## 2."): on = True; continue
        if on and l.startswith("| A1v2_20260914 |"):
            p = [x.strip() for x in l.strip().strip("|").split("|")]
            out[p[1]] = (int(p[5]), p[6], p[7], p[10])
    return out


def main():
    workers = int(sys.argv[sys.argv.index("--workers") + 1]) if "--workers" in sys.argv else 6
    want = {f"{e:.4f}" for e, _, _ in P2.A1V2}
    inv = []
    for cell, ms in W.discover().items():                                   # sec. 4.7.11's (= sec. 4.6's) discovery and filters
        if W.campaign_of(cell) != "A1v2_20260914": continue
        M0 = sorted(ms)[0]; mt = W.meta(cell, sorted(ms[M0])[0][1])
        if mt is None: continue
        L0, en, Ns = mt["L0"], mt["eta"], mt["Ns"]
        et = en * L0 / (L0 - box_delta(L0) / 2)
        if f"{et:.4f}" not in want: continue
        if not (W.date_of(cell)[:10] >= W.NEW_CORE and len(ms) >= 3 and mt["complete"] and Ns > 0): continue
        pr = W.psi6_runs(cell)
        if not pr: continue
        bad, _ = W.bad_runs(cell)
        inv.append(dict(cell=cell, eta=et, ms=ms, bad=bad, pr=pr))
    inv.sort(key=lambda q: q["eta"])
    jobs = [(q["cell"], M, r, p) for q in inv for M, runs in q["ms"].items() for r, p in sorted(runs)
            if (M, r) not in q["bad"] and T.trace_check(p)[0]]
    with Pool(workers) as pool:
        out = pool.map(WA.traj2, [p for _, _, _, p in jobs], chunksize=8)
    res = defaultdict(lambda: defaultdict(list)); refused = 0
    for (cell, M, r, p), o in zip(jobs, out):
        if o == "refused": refused += 1; continue
        if o is not None: o["r"] = r; o["p"] = p; res[cell][M].append(o)
    rec = recorded(); gate_bad = 0
    print("# A1 v2 at P2's masses (261012 sec. 4.7.19), printed by experiments_gen3_p2_261009/p2_a1v2_rows.py -- existing data, "
          "sec. 4.7.11's cells, trajectories and estimator\n")
    if refused: print(f"refused by the loader provenance guard: {refused}\n")
    print(f"| eta_true | M | trajectories | record [sigma-time], mean | mean dnu [%] (SE) | dnu exactly 0 | mean dA [%] (SE) | "
          f"psi6 at the end of the hold ({HOLD_SIGMA:.1f} sigma-time) | psi6 at the end of the record |\n|---|---|---|---|---|---|---|---|---|")
    for q in inv:
        R = res[q["cell"]]; allr = []
        for M in sorted(R):
            if len(R[M]) < 3: continue
            tr = [o for o in R[M] if (M, o["r"]) in q["pr"]]
            if len(tr) < 3: continue
            allr += [(M, o) for o in tr]
        for sel in [m for m in P2.MASSES] + ["all nine"]:
            rows = [o for M, o in allr if sel == "all nine" or M == sel]
            if not rows: print(f"| {q['eta']:.4f} | {sel} | 0 | | | | | | |"); continue
            dn = [100 * (o["nu2"] / o["nu1"] - 1) for o in rows]; da = [100 * (o["a2"] / o["a1"] - 1) for o in rows]
            n0 = sum(1 for o in rows if o["nu2"] == o["nu1"])
            Ms = [M for M, o in allr if sel == "all nine" or M == sel]
            ph = [q["pr"][(M, o["r"])][0] for M, o in allr if sel == "all nine" or M == sel]
            pe = [q["pr"][(M, o["r"])][1] for M, o in allr if sel == "all nine" or M == sel]
            rl = [o["TD"] / float(W.header(o["p"])["Predicted_Frequency"]) for o in rows]   # the analysed prefix: TD predicted periods
            (m1, s1), (m2, s2), (h1, hs), (e1, es), (r1, _) = mse(dn), mse(da), mse(ph), mse(pe), mse(rl)
            gate = ""
            if sel == "all nine":
                key = f"{q['eta']:.4f}"
                mine = (len(rows), f"{m1:+.2f} ({s1:.2f})", f"{n0} ({100.0 * n0 / len(rows):.0f} %)", f"{m2:+.2f} ({s2:.2f})")
                ok = rec.get(key) == mine; gate_bad += (not ok)
                gate = " (= sec. 4.7.11)" if ok else f" (**sec. 4.7.11 prints {rec.get(key)}**)"
            print(f"| {q['eta']:.4f} | {sel} | {len(rows)}{gate} | {r1:.0f} | {m1:+.2f} ({s1:.2f}) | {n0} ({100.0 * n0 / len(rows):.0f} %) | "
                  f"{m2:+.2f} ({s2:.2f}) | {h1:.3f} ({hs:.3f}) | {e1:.3f} ({es:.3f}) |")
    print(f"\ngate: the 'all nine' rows against sec. 4.7.11's printed table 2: "
          f"{'IDENTICAL to the printed digits (4 of 4 eta)' if gate_bad == 0 and len(inv) == 4 else f'**{gate_bad} differ, {len(inv)} of 4 eta found**'}")
    print("(record = the analysed prefix, TD = 200 predicted periods, as P2's records: at the same eta and M the two have the same predicted "
          "length. psi6 = the global |psi6|, the driver's per-run summary.)")


if __name__ == "__main__":
    main()
