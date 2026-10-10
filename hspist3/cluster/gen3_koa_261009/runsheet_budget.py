#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.29; plan-author decision 13, item 2f): the core-hours and the largest core count of every step of
the KOA runsheet, read from the job scripts' #SBATCH headers (time limit, cores per task, partition) and the runsheet's submission
lines (the array size and its %-limit), so the numbers cannot drift from the scripts. "Requested" = cores x the time limit, summed
over the array tasks: the most a step can use. As information, the Mac's measured core-hours for the same Test G task list (the Mac's
Test G run logs: mean seconds per run per cell, mass and engine, times the KOA list's counts); KOA's per-core speed for gen3 is what the
gate's rate table (step 3, part c) measures.
usage (from hspist3/): python3 cluster/gen3_koa_261009/runsheet_budget.py
"""
import glob, os, re, statistics
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE))


def hdr(name):
    h = {}
    for l in open(os.path.join(HERE, name)):
        m = re.match(r"#SBATCH --([a-z-]+)=(\S+)", l)
        if m: h[m.group(1)] = m.group(2)
    t = [int(x) for x in h["time"].split(":")]
    return h["partition"], int(h.get("cpus-per-task", 1)), t[0] + t[1] / 60 + (t[2] / 3600 if len(t) > 2 else 0)


STEPS = [   # step, what, [(script or None for the interactive shell, tasks, at once)]
    ("1-2", "sandbox shell: clone and the two builds", [(("sandbox", 2, 1.0), 1, 1)]),
    ("3", "the one-node gate job (E0, harness hashes, rate table, long double)", [("g3_gate_koa.sbatch", 1, 1)]),
    ("4", "cross-node determinism: array 1-8 and the compare job", [("g3_det_xnode.sbatch", 8, 8), ("g3_xnode_compare.sbatch", 1, 1)]),
    ("5", "Test G on KOA: array 1-144%64 and the merge job", [("testG_koa.sbatch", 144, 64), ("merge_testG_koa.sbatch", 1, 1)]),
    ("6", "the fetch (Mac) and the verdict (Mac)", []),
]


def mac_testG_core_hours():
    root = os.path.join(HS, "experiments_gen3_gate_261009", "testG")
    if not os.path.isdir(root): return None
    sec = {}
    for eng in ("gen2", "gen3"):
        for d in glob.glob(os.path.join(root, eng, "*", "m_*")):
            p = os.path.join(d, "run.log")
            if os.path.exists(p):
                s = [int(x) for x in re.findall(r"^##RUN .*? \((\d+) s\)", open(p, errors="ignore").read(), re.M)]
                if s: sec[(eng, os.path.basename(os.path.dirname(d)), int(d.rsplit("_", 1)[1]))] = statistics.mean(s)
        s = []
        for f in glob.glob(os.path.join(root, eng, "*", "x_0", "run_*.log")):
            s += [int(x) for x in re.findall(r"^##RUN .*? \((\d+) s\) reduce", open(f, errors="ignore").read(), re.M)]
        if s: sec[(eng, "afix", 0)] = statistics.mean(s)
    n = Counter()
    for l in open(os.path.join(HERE, "tasks_testG_koa_261009.txt")):
        p = l.split()
        if p[0] == "B": n[(p[-2], p[1].split("/")[2], int(p[2]))] += 1
        else: n[(p[-2], "afix", 0)] += 1
    missing = [k for k in n if k not in sec]
    return sum(n[k] * sec[k] for k in n if k in sec) / 3600, missing, max(sec.values())


def main():
    print("| step | what | partition | largest core count at once | requested core-hours (cores x time limit, summed) |\n|---|---|---|---|---|")
    tot = 0.0
    for st, what, jobs in STEPS:
        parts, cores, ch = [], 0, 0.0
        for j, ntask, conc in jobs:
            part, c, h = j if isinstance(j, tuple) else hdr(j)
            parts.append(part); ch += ntask * c * h; cores = max(cores, conc * c)
        tot += ch
        print(f"| {st} | {what} | {', '.join(sorted(set(parts))) or 'none (Mac)'} | {cores} | {ch:.2f} |")
    print(f"| all | | | | {tot:.2f} |")
    m = mac_testG_core_hours()
    if m:
        h, missing, mx = m
        print(f"\nInformation: the Mac's measured core-hours for the same Test G task list (the Mac's Test G run logs; the M = 50 main runs at "
              f"the Mac's M = 50 rate): {h:.2f} core-hours; the slowest Mac run averaged {mx:.0f} s"
              + (f"; cells without a Mac rate: {missing}" if missing else ""))
        # the Mac's gen2 rate in Test G's M = 300 runs (events from [EDMD-CONTACT], seconds from the ##RUN header) against KOA's measured
        # gen2 (7b08827) rate per core at N = 100, pi/8: 1.1e5 events/s (261012 sec. 4.7, "Where we are [DATA, KOA, measured]")
        p = os.path.join(HS, "experiments_gen3_gate_261009", "testG", "gen2", "epi8_H_H10_L10", "m_300", "run.log")
        if os.path.exists(p):
            secs = re.split(r"^##RUN ", open(p, errors="ignore").read(), flags=re.M)[1:]
            ev = [int(re.search(r"\[EDMD-CONTACT\] executed events (\d+)", s).group(1)) for s in secs]
            sc = [int(re.match(r".*? \((\d+) s\)", s).group(1)) for s in secs]
            rate = sum(ev) / sum(sc); f = rate / 1.1e5
            print(f"Information: the Mac ran gen2 at {rate:.3g} events/s per core (Test G, M = 300, {len(ev)} runs); KOA at 1.1e5 (measured, "
                  f"7b08827): a factor {f:.2f}, so the Test G array would take about {h * f:.0f} core-hours on KOA if gen3 slows by the same "
                  f"factor (step 3's rate table measures gen3); its slowest chunk, about {25 * mx * f / 3600:.1f} h, against the 6 h limit")


if __name__ == "__main__":
    main()
