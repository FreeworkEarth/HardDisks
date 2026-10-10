#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.29; plan-author decision 13, item 2d): PREPARED, NOT RUN. The compare step of the cross-node
determinism test: reads every task folder task_<i>/ of one g3_det_xnode.sbatch array (its node.txt and, per cell, the divider trace,
the psi6(t) file, speed_of_sound_psi6.csv and the run record's event hash), groups the tasks by node, and compares:
  - ACROSS NODES (the test): for every pair of distinct nodes, the first task of each; every compared file must be byte-identical and
    the event hashes equal;
  - within a node (information): every further task against that node's first task.
It prints the node pairs compared, IDENTICAL / DIFFERENT per file and pair, the tasks that did not finish (exit != 0, clean != 1, a
missing file), and "NOT TESTED" if every task ran on one node (then nothing was tested across nodes). Read-only.
usage: python3 cluster/gen3_koa_261009/g3_xnode_compare.py <the array's folder .../gen3_xnode_<array job id>>
"""
import filecmp, glob, os, re, sys
from itertools import combinations


def task_files(t, cell):
    d = os.path.join(t, cell)
    fs = sorted(os.path.basename(p) for p in glob.glob(os.path.join(d, "wall_x_positions_*.csv")) + glob.glob(os.path.join(d, "psi6_t_*.csv")))
    if os.path.exists(os.path.join(d, "speed_of_sound_psi6.csv")): fs.append("speed_of_sound_psi6.csv")
    lg = open(os.path.join(d, "run.log"), errors="ignore").read() if os.path.exists(os.path.join(d, "run.log")) else ""
    h = re.findall(r"^\[EDMD3-HEALTH\] .*? clean=(\S+) .* hash=([0-9a-f]+) ", lg, re.M)
    ex = open(os.path.join(d, "exit.txt")).read().strip() if os.path.exists(os.path.join(d, "exit.txt")) else "no exit.txt"
    return fs, (h[0] if len(h) == 1 else None), ex


def main():
    root = os.path.abspath(sys.argv[1])
    tasks = sorted(glob.glob(os.path.join(root, "task_*")), key=lambda p: int(p.rsplit("_", 1)[1]))
    print(f"cross-node determinism: {len(tasks)} tasks in {root}")
    if not tasks: print("NOT TESTED: no task folder"); return
    node = {t: open(os.path.join(t, "node.txt")).read().strip() if os.path.exists(os.path.join(t, "node.txt")) else "unknown" for t in tasks}
    host = {t: node[t].split(" | ")[0] for t in tasks}
    for t in tasks: print(f"  {os.path.basename(t)}: {node[t]}")
    cells = sorted({os.path.basename(c) for t in tasks for c in glob.glob(os.path.join(t, "*")) if os.path.isdir(c)})
    ok_tasks = []
    for t in tasks:
        bad = []
        for c in cells:
            fs, h, ex = task_files(t, c)
            if ex != "exit 0" or not h or h[0] != "1" or len(fs) != 3: bad.append(f"{c}: {ex}, record {h}, files {len(fs)}")
        if bad: print(f"  {os.path.basename(t)} NOT USED: {'; '.join(bad)}")
        else: ok_tasks.append(t)
    first = {}
    for t in ok_tasks: first.setdefault(host[t], t)
    def compare(a, b, label):
        allsame = True
        for c in cells:
            fa, ha, _ = task_files(a, c); fb, hb, _ = task_files(b, c)
            for f in fa:
                same = f in fb and filecmp.cmp(os.path.join(a, c, f), os.path.join(b, c, f), shallow=False)
                allsame &= same
                print(f"     {c}/{f}: {'IDENTICAL' if same else '**DIFFERENT**'}")
            same = ha[1] == hb[1]; allsame &= same
            print(f"     {c} event hash: {ha[1]} / {hb[1]}: {'IDENTICAL' if same else '**DIFFERENT**'}")
        print(f"   {label}: {'IDENTICAL' if allsame else '**DIFFERENT**'}")
        return allsame
    print(f"\n== across nodes (the test): {len(first)} distinct node(s) among the {len(ok_tasks)} finished tasks")
    if len(first) < 2:
        print("NOT TESTED: every finished task ran on one node, so nothing was compared across nodes")
        across = None
    else:
        across = True
        for (ha, a), (hb, b) in combinations(sorted(first.items()), 2):
            print(f"  pair {ha} ({os.path.basename(a)}) vs {hb} ({os.path.basename(b)}):")
            across &= compare(a, b, f"{ha} vs {hb}")
    print("\n== within a node (information)")
    within = True
    for t in ok_tasks:
        f = first[host[t]]
        if t == f: continue
        print(f"  {os.path.basename(t)} vs {os.path.basename(f)} (both on {host[t]}):")
        within &= compare(t, f, f"{os.path.basename(t)} vs {os.path.basename(f)}")
    print(f"\nRESULT: across nodes {'NOT TESTED' if across is None else ('IDENTICAL' if across else '**DIFFERENT**')} "
          f"({len(first)} node(s)); within nodes {'IDENTICAL' if within else '**DIFFERENT**'}; tasks used {len(ok_tasks)} of {len(tasks)}")


if __name__ == "__main__":
    main()
