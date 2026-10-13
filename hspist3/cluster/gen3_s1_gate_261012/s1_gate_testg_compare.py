#!/usr/bin/env python3
"""##CHRIS 2026-10-12 (plan-author decision 39 item 5 c; gate G6 (d)): compares the 12 gen3 lines of Test G chunk 1, run once with the gate
binary (a643bc7) and once with the campaign binary (f72d764) by Test G's worker into two roots: per line, the run record's event hash
and its event counts ([EDMD3-HEALTH] in the cell's run.log, the ##RUN section of that line's run index), and clean. Read-only.
Prints a table and the line "TEST G CHUNK 1 (12 gen3 lines): event hashes equal in <n> of 12" (G6 (d) asks 12 of 12).
usage: python3 s1_gate_testg_compare.py <gate root> <campaign root> <the 12 task lines>
"""
import os, re, sys
EVK = ("ev_pair", "ev_wall", "ev_cross", "ev_div", "ev_piston", "ev_band")


def record(root, rel, r):
    p = os.path.join(root, rel, "run.log")
    if not os.path.exists(p): return None
    txt = open(p, errors="replace").read()
    sec = re.search(rf"^##RUN \S+ \S+ run {r} seed \d+ .*?(?=^##RUN |\Z)", txt, re.M | re.S)
    if not sec: return None
    m = re.findall(r"^\[EDMD3-HEALTH\] .*?: clean=(\S+) (.*)$", sec.group(0), re.M)
    if len(m) != 1: return None
    kv = dict(re.findall(r"(\w+)=(\S+)", m[0][1])); kv["clean"] = m[0][0]
    return kv


def main():
    groot, croot, lines = sys.argv[1:4]
    L = [l.split() for l in open(lines) if l.strip()]
    print("| line | cell | M | r | seed | gate: clean, event hash | campaign: clean, event hash | event counts equal | equal |\n|---|---|---|---|---|---|---|---|---|")
    neq = 0
    for i, f in enumerate(L, 1):
        rel, M, r, seed = f[1], f[2], f[3], f[4]
        a, b = record(groot, rel, r), record(croot, rel, r)
        cnt = a is not None and b is not None and all(a.get(k) == b.get(k) for k in EVK)
        eq = a is not None and b is not None and a.get("hash") == b.get("hash") and a["clean"] == b["clean"] == "1" and cnt
        neq += eq
        print(f"| {i} | {rel} | {M} | {r} | {seed} | {(a['clean'] + ', ' + a.get('hash', '-')) if a else 'no record'} | "
              f"{(b['clean'] + ', ' + b.get('hash', '-')) if b else 'no record'} | {'yes' if cnt else 'no'} | {'EQUAL' if eq else '**DIFFERENT**'} |")
    print(f"\nTEST G CHUNK 1 ({len(L)} gen3 lines): event hashes equal in {neq} of {len(L)}")


if __name__ == "__main__":
    main()
