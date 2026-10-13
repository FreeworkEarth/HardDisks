#!/usr/bin/env python3
"""##CHRIS 2026-10-12 (plan-author decision 41 item 15 b; information): the 12 gen3 lines of Test G chunk 1 in three run records -- Test G on
KOA (job K2-3, binary a643bc7; main's hspist3/experiments_gen3_gate_koa_261009/testG, fetched 2026-10-11), and this gate's two runs of the
same lines (the gate binary a643bc7 and the campaign binary f72d764; experiments_gen3_s1_gate_koa_261012/s1_gate/testg_gate, testg_campaign).
Per line: clean, the event hash and the event counts of the run record [EDMD3-HEALTH] (s1_gate_testg_compare.record, unchanged). Read-only.
usage (from hspist3/ of engine-gen3): python3 experiments_gen3_s1_261012/gate_koa/check_15b_testG_koa.py <main's experiments_gen3_gate_koa_261009>
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(HS, "cluster", "gen3_s1_gate_261012"))
import s1_gate_testg_compare as C


def main():
    koa = os.path.abspath(sys.argv[1]); G = os.path.join(HS, "experiments_gen3_s1_gate_koa_261012", "s1_gate")
    lines = [l.split() for l in open(os.path.join(G, "testg_chunk1_gen3_lines.txt")) if l.strip()]
    print(f"# Decision 41 item 15 b: the 12 gen3 lines of Test G chunk 1 in Test G on KOA ({os.path.basename(koa)}) and in this gate's two runs\n")
    print("| line | M | r | seed | Test G on KOA (a643bc7): clean, event hash | gate binary (a643bc7) | campaign binary (f72d764) | event counts equal (3 records) | equal |")
    print("|---|---|---|---|---|---|---|---|---|")
    n = 0
    for i, f in enumerate(lines, 1):
        rel, M, r, seed = f[1], f[2], f[3], f[4]
        recs = [C.record(root, rel, r) for root in (koa, os.path.join(G, "testg_gate"), os.path.join(G, "testg_campaign"))]
        ok = all(x is not None for x in recs)
        cnt = ok and all(all(x.get(k) == recs[0].get(k) for k in C.EVK) for x in recs)
        eq = ok and cnt and len({x.get("hash") for x in recs}) == 1 and all(x["clean"] == "1" for x in recs)
        n += eq
        cell = lambda x: f"{x['clean']}, {x.get('hash', '-')}" if x else "no record"
        print(f"| {i} | {M} | {r} | {seed} | {cell(recs[0])} | {cell(recs[1])} | {cell(recs[2])} | {'yes' if cnt else 'no'} | {'EQUAL' if eq else '**DIFFERENT**'} |")
    print(f"\nRESULT (information): event hashes and event counts equal in all three records for {n} of {len(lines)} lines")


if __name__ == "__main__":
    main()
