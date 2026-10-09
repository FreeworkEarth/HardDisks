#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.15, stage B of the plan-author programme of sec. 4.7.12): the task list of TEST G, the Mac
part of the gen-3 gate, written and committed before any of its trajectories.
Test G: cell epi8_H_H10_L10 (eta = pi/8, H = L0 = 10, N = 100), the speed-of-sound protocol and command of Test T-prime (the
campaign worker's B command, cluster/confinement_20261013/conf_worker.sh; geometry and stride per mass from the campaign's task
list tasks_B_epi8_H_H10_L10.txt, as gen_testTprime_261007.py) except the seeds and --engine; engines gen2 (the default policy)
and gen3, one binary, one build; masses M = 300 and M = 1500; 400 trajectories per engine per mass.
SEEDS: a new base never used before, BASE = 20261209 (no file of main, engine-gen3 or this session's scratch contains the number
before this script; the derived seeds are checked against every task list by --design); run_seed(BASE, l, mass index, r) with
the stream l per block (below); the SAME seed for both engines; each (block, mass, r) as two consecutive lines, gen2 then gen3,
so the engines run interleaved (the Mac runner keeps the task order).
BLOCKS (stream l; the first is the test, the others information rows without a verdict):
  main     l = 0   epi8_H_H10_L10, M = 300 and 1500, r = 0..399 (the extension, if the rule asks for it: r = 400..799 of the same
                   stream at that mass, --extension M; its seed-list SHA-256 is printed now)
  info_m   l = 0   epi8_H_H10_L10, M = 50 and 2000, r = 0..99 (other mass index, so other seeds)
  dense    l = 1   the gate's dense geometry (audit_runs_261007.case_cmd "dense_M": L0 = 269/48, H = 10, Ns = 50, eta = 0.70073,
                   stride T.d_stride from the KR prediction), M = 300, r = 0..199; the Test G protocol (200 oscillations)
  n400     l = 2   eta = pi/8 at N = 400 (H = L0 = 20, Ns = 200), M = 300, r = 0..99, stride T.d_stride from the KR prediction
  afix     l = 3   the static method: the A-fixed cell epi8_H_H10_L10 at x_0 = 10 (the campaign's AF command: hold 312000 steps,
                   post 1200, trace every 600), the seed run_seed(BASE, 3, 0, r) as --seed, r = 0..199
Output: cluster/gen3_gate_261009/tasks_testG_261009.txt (and with --extension M: tasks_testG_ext_M<M>_261009.txt), lines
  B  <rel>/<engine>/<cid>/m_<M> <M> <r> <seed> <L0> <H> <Ns> <stride> <BASE> <engine> <block>
  AF <rel>/<engine>/epi8_H_H10_L10/x_0 10.000000 <seed> 10.000000 10.000000 50 312000 1200 600 <engine> afix
with rel = testG (the data root is hspist3/experiments_gen3_gate_261009/, untracked).
usage (from hspist3/): python3 cluster/gen3_gate_261009/gen_testG_261009.py [--extension 300|1500]
"""
import glob, hashlib, math, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); CL = os.path.dirname(HERE); HS = os.path.dirname(CL)
sys.path.insert(0, os.path.join(HS, "validation")); sys.path.insert(0, HS)
import tests_20260913 as T

CID, BASE, NSEED, MASSES = "epi8_H_H10_L10", 20261209, 400, (300, 1500)
ENGINES = ("gen2", "gen3")
REL = "testG"
LOC = os.path.join(HS, "experiments_gen3_gate_261009")
L0D = 269 / 48                                                    # the gate's dense cell (audit_runs_261007.case_cmd)
ETA_D = 50 * math.pi * 0.25 / (10.0 * L0D)


def campaign_geo():
    geo = {}
    for l in open(os.path.join(CL, "confinement_20261013", f"tasks_B_{CID}.txt")):
        f = l.split(); geo.setdefault(int(f[2]), (f[5], f[6], f[7], f[8]))      # L0 H Ns stride, per mass
    return geo


def blocks(ext=None):
    """[(block, cid, M, r, seed, L0, H, Ns, stride)] in task order (B tasks), and the AF seeds."""
    geo = campaign_geo(); mi = T.A1_MASSES.index
    out = []
    if ext is None:
        for M in MASSES:
            L0, H, Ns, st = geo[M]
            out += [("main", CID, M, r, T.run_seed(BASE, 0, mi(M), r), L0, H, Ns, st) for r in range(NSEED)]
        for M in (50, 2000):
            L0, H, Ns, st = geo[M]
            out += [("info_m", CID, M, r, T.run_seed(BASE, 0, mi(M), r), L0, H, Ns, st) for r in range(100)]
        st = T.d_stride(T.kr_cs(ETA_D) * T.x_of(300, L0D))
        out += [("dense", "dense_H10_L5.604167", 300, r, T.run_seed(BASE, 1, mi(300), r), f"{L0D:.6f}", "10.000000", "50", str(st))
                for r in range(200)]
        st = T.d_stride(T.kr_cs(math.pi / 8) * T.x_of(300, 20.0))
        out += [("n400", "epi8_N400_H20_L20", 300, r, T.run_seed(BASE, 2, mi(300), r), "20.000000", "20.000000", "200", str(st)) for r in range(100)]
        af = [T.run_seed(BASE, 3, 0, r) for r in range(200)]
    else:
        L0, H, Ns, st = geo[ext]
        out = [("main", CID, ext, r, T.run_seed(BASE, 0, mi(ext), r), L0, H, Ns, st) for r in range(NSEED, 2 * NSEED)]
        af = []
    return out, af


def lines(ext=None):
    B, af = blocks(ext); L = []
    for blk, cid, M, r, s, L0, H, Ns, st in B:
        for e in ENGINES:
            L.append(f"B {REL}/{e}/{cid}/m_{M} {M} {r} {s} {L0} {H} {Ns} {st} {BASE} {e} {blk}\n")
    for s in af:
        for e in ENGINES:
            L.append(f"AF {REL}/{e}/{CID}/x_0 10.000000 {s} 10.000000 10.000000 50 312000 1200 600 {e} afix\n")
    return L


def seed_text(ext=None):
    B, af = blocks(ext)
    return "".join(f"{blk} {cid} {M} {r} {s}\n" for blk, cid, M, r, s, *_ in B) + "".join(f"afix {r} {s}\n" for r, s in enumerate(af))


def used_seeds():
    """Every exact seed of the earlier task lists (campaign B/A/AF, Test T, T-prime) and of today's stage A runs."""
    used = set()
    for f in glob.glob(os.path.join(CL, "*", "tasks_*.txt")):
        if os.path.dirname(os.path.abspath(f)) == HERE: continue          # Test G's own lists
        for l in open(f):
            p = l.split()
            if not p: continue
            if p[0] == "B": used.add(int(p[4]))
            elif p[0] in ("A", "AF"): used.add(int(p[3]))
    used |= {T.run_seed(20261013, 0, m, r) for m in range(9) for r in range(4)}            # smoke / pilot / E2
    used |= {T.run_seed(20261091, 5, 0, r) for r in range(8)} | {9700, 9701, 9702, 9703, 9711}   # stage A (acc5), A-fixed, A4
    return used


def main():
    ext = int(sys.argv[sys.argv.index("--extension") + 1]) if "--extension" in sys.argv else None
    if ext is not None and ext not in MASSES: sys.exit("--extension takes 300 or 1500")
    L = lines(ext)
    out = os.path.join(HERE, "tasks_testG_261009.txt" if ext is None else f"tasks_testG_ext_M{ext}_261009.txt")
    if os.path.exists(out) and open(out).read() != "".join(L): sys.exit(f"STOP: {out} exists with other content -- not overwriting")
    open(out, "w").write("".join(L))
    B, af = blocks(ext); seeds = [b[4] for b in B] + af
    print(f"{os.path.relpath(out, HS)}: {len(L)} lines; seed list SHA-256 {hashlib.sha256(seed_text(ext).encode()).hexdigest()}; "
          f"task list SHA-256 {hashlib.sha256(''.join(L).encode()).hexdigest()}")
    from collections import Counter
    print("lines per block and engine: " + ", ".join(f"{k} {v}" for k, v in sorted(Counter((l.split()[-1], l.split()[-2]) for l in L).items())))
    used = used_seeds()
    print(f"seeds: {len(seeds)}, distinct {len(set(seeds))}; overlap with the earlier task lists and today's stage A runs: "
          f"{len(set(seeds) & used)} of {len(used)}")
    if ext is None:
        for M in MASSES:
            E = "".join(f"main {CID} {M} {r} {T.run_seed(BASE, 0, T.A1_MASSES.index(M), r)}\n" for r in range(NSEED, 2 * NSEED))
            print(f"extension seed list at M = {M} (r = {NSEED}..{2 * NSEED - 1}, lines 'main {CID} M r seed'): SHA-256 "
                  f"{hashlib.sha256(E.encode()).hexdigest()}; overlap with the main list and the earlier seeds: "
                  f"{len({T.run_seed(BASE, 0, T.A1_MASSES.index(M), r) for r in range(NSEED, 2 * NSEED)} & (set(seeds) | used))}")


if __name__ == "__main__":
    main()
