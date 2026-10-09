#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.20; stage G): PREPARED, NOT RUN, NOT REGISTERED. The task list of TEST G ON KOA: the design of
stage B (cluster/gen3_gate_261009/gen_testG_261009.py, its blocks and its order, gen2 then gen3 per seed) with a new base seed
BASE_KOA = 20261220 (never used before), and the CHUNKS of the array job: chunk j holds task lines j*K .. j*K + K - 1 (K = 25), one
core each, each chunk writing its own folder (testG_koa.sbatch), at most 64 chunks at once.
It must be REGISTERED (this text, the task list with its SHA-256, the KOA verdict wrapper testG_koa_verdict_261009.py) by the plan
author's go before any KOA trajectory exists, as stage B was.
Output: cluster/gen3_koa_261009/tasks_testG_koa_261009.txt (3000 lines, the Mac's format) and the number of chunks.
usage (from hspist3/): python3 cluster/gen3_koa_261009/gen_testG_koa_261009.py [--extension 300|1500]
"""
import glob, hashlib, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); CL = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(CL, "gen3_gate_261009"))
import gen_testG_261009 as G
BASE_KOA, K = 20261220, 25


def main():
    G.BASE = BASE_KOA                                   # the same blocks with the KOA base (run_seed reads G.BASE at call time)
    ext = int(sys.argv[sys.argv.index("--extension") + 1]) if "--extension" in sys.argv else None
    L = G.lines(ext)
    # the extension list keeps the registered verdict's file name (its inventory globs tasks_testG_ext_M*_261009.txt in V.GG = here)
    out = os.path.join(HERE, "tasks_testG_koa_261009.txt" if ext is None else f"tasks_testG_ext_M{ext}_261009.txt")
    if os.path.exists(out) and open(out).read() != "".join(L): sys.exit(f"STOP: {out} exists with other content -- not overwriting")
    open(out, "w").write("".join(L))
    used = set()                                        # every earlier list, the Mac's Test G included; never this folder's own lists
    for f in glob.glob(os.path.join(CL, "*", "tasks_*.txt")):
        if os.path.dirname(os.path.abspath(f)) == HERE: continue
        for l in open(f):
            p = l.split()
            if p and p[0] == "B": used.add(int(p[4]))
            elif p and p[0] in ("A", "AF"): used.add(int(p[3]))
    used |= G.used_seeds()
    seeds = [int(l.split()[4]) if l.startswith("B ") else int(l.split()[3]) for l in L]
    print(f"{os.path.relpath(out, os.path.dirname(CL))}: {len(L)} lines, SHA-256 {hashlib.sha256(''.join(L).encode()).hexdigest()}; "
          f"chunks of {K}: {-(-len(L) // K)} (array 1-{-(-len(L) // K)}%64); distinct seeds {len(set(seeds))}; overlap with every earlier "
          f"list incl. the Mac's Test G: {len(set(seeds) & used)}")


if __name__ == "__main__":
    main()
