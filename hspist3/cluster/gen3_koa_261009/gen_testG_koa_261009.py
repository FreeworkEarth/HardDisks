#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.20; stage G), REGISTRATION of TEST G ON KOA by plan-author decision 13, item 2g (261012
sec. 4.7.29): this script, the task list it writes (tasks_testG_koa_261009.txt) and that list's SHA-256 (tasks_testG_koa_261009.sha256,
read by testG_koa.sbatch) are committed together, before any KOA trajectory exists; none of them is edited afterwards.
DESIGN: stage B's (cluster/gen3_gate_261009/gen_testG_261009.py, the Mac's Test G: its cell epi8_H_H10_L10, its commands and stride
per mass, its engines, its order -- gen2 then gen3 per seed, the same seed for both -- and its information blocks), with the two
changes of decision 13, item 2g:
  - THREE verdict masses, M = 50, 300 and 1500, 400 trajectories per engine each (block "main", stream l = 0 at each mass index);
  - the information rows as on the Mac WITHOUT the M = 50 row: M = 2000 (100 per engine, block "info_m"), the dense cell (200), N = 400
    at pi/8 (100), the static method (the A-fixed cell at x_0, 200 runs per engine).
SEEDS: BASE_KOA = 20261220 (no file contained the number before stage G; reserved since then for this list), run_seed(BASE_KOA, l, mass
index, r) per block as on the Mac; the derived seeds are checked against every earlier task list, the Mac's Test G included, and
against the seed list of the symmetry test (experiments_gen3_sym_261009/sym_seeds.txt).
EXTENSIONS, one per mass as before (sec. 4.7.15): r = 400..799 of the same stream at that mass, only if the rule asks for it at that
mass; each extension's seed-list SHA-256 is printed now (the registration); its task list (--extension M) and its run are prepared then.
CHUNKS of the array job: chunk j holds task lines (j-1)*25+1 .. j*25, one core each, each chunk writing its own folder (testG_koa.sbatch),
at most 64 chunks at once: 3600 lines = 144 chunks.
Output: cluster/gen3_koa_261009/tasks_testG_koa_261009.txt (the Mac's line format) and tasks_testG_koa_261009.sha256 (one line, the hex
digest); with --extension M: tasks_testG_ext_M<M>_261009.txt.
usage (from hspist3/): python3 cluster/gen3_koa_261009/gen_testG_koa_261009.py [--extension 50|300|1500]
"""
import glob, hashlib, math, os, sys
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__)); CL = os.path.dirname(HERE); HS = os.path.dirname(CL)
sys.path.insert(0, os.path.join(CL, "gen3_gate_261009")); sys.path.insert(0, os.path.join(HS, "validation")); sys.path.insert(0, HS)
import gen_testG_261009 as G
import tests_20260913 as T
BASE_KOA, K, MASSES_KOA, INFO_M = 20261220, 25, (50, 300, 1500), (2000,)
TASKS = os.path.join(HERE, "tasks_testG_koa_261009.txt")


def blocks(ext=None):
    """stage B's blocks (gen_testG_261009.blocks) with the three verdict masses and without the M = 50 information row"""
    geo = G.campaign_geo(); mi = T.A1_MASSES.index; CID = G.CID; NSEED = G.NSEED
    out = []
    if ext is None:
        for M in MASSES_KOA:
            L0, H, Ns, st = geo[M]
            out += [("main", CID, M, r, T.run_seed(BASE_KOA, 0, mi(M), r), L0, H, Ns, st) for r in range(NSEED)]
        for M in INFO_M:
            L0, H, Ns, st = geo[M]
            out += [("info_m", CID, M, r, T.run_seed(BASE_KOA, 0, mi(M), r), L0, H, Ns, st) for r in range(100)]
        st = T.d_stride(T.kr_cs(G.ETA_D) * T.x_of(300, G.L0D))
        out += [("dense", "dense_H10_L5.604167", 300, r, T.run_seed(BASE_KOA, 1, mi(300), r), f"{G.L0D:.6f}", "10.000000", "50", str(st))
                for r in range(200)]
        st = T.d_stride(T.kr_cs(math.pi / 8) * T.x_of(300, 20.0))
        out += [("n400", "epi8_N400_H20_L20", 300, r, T.run_seed(BASE_KOA, 2, mi(300), r), "20.000000", "20.000000", "200", str(st))
                for r in range(100)]
        af = [T.run_seed(BASE_KOA, 3, 0, r) for r in range(200)]
    else:
        L0, H, Ns, st = geo[ext]
        out = [("main", CID, ext, r, T.run_seed(BASE_KOA, 0, mi(ext), r), L0, H, Ns, st) for r in range(NSEED, 2 * NSEED)]
        af = []
    return out, af


def lines(ext=None):
    B, af = blocks(ext); L = []
    for blk, cid, M, r, s, L0, H, Ns, st in B:
        for e in G.ENGINES:
            L.append(f"B {G.REL}/{e}/{cid}/m_{M} {M} {r} {s} {L0} {H} {Ns} {st} {BASE_KOA} {e} {blk}\n")
    for s in af:
        for e in G.ENGINES:
            L.append(f"AF {G.REL}/{e}/{G.CID}/x_0 10.000000 {s} 10.000000 10.000000 50 312000 1200 600 {e} afix\n")
    return L


def ext_seed_text(M):
    return "".join(f"main {G.CID} {M} {r} {T.run_seed(BASE_KOA, 0, T.A1_MASSES.index(M), r)}\n" for r in range(G.NSEED, 2 * G.NSEED))


def main():
    ext = int(sys.argv[sys.argv.index("--extension") + 1]) if "--extension" in sys.argv else None
    if ext is not None and ext not in MASSES_KOA: sys.exit("--extension takes 50, 300 or 1500")
    L = lines(ext); txt = "".join(L); sha = hashlib.sha256(txt.encode()).hexdigest()
    out = TASKS if ext is None else os.path.join(HERE, f"tasks_testG_ext_M{ext}_261009.txt")
    if os.path.exists(out) and open(out).read() != txt: sys.exit(f"STOP: {out} exists with other content -- not overwriting")
    open(out, "w").write(txt)
    if ext is None:
        shp = os.path.join(HERE, "tasks_testG_koa_261009.sha256")
        if os.path.exists(shp) and open(shp).read() != sha + "\n": sys.exit(f"STOP: {shp} exists with other content")
        open(shp, "w").write(sha + "\n")
    used = set()                                        # every earlier list (campaign, Test T, T-prime, the Mac's Test G), never this
    for f in glob.glob(os.path.join(CL, "*", "tasks_*.txt")):   # folder's own lists (G.used_seeds() skips only the Mac's folder)
        if os.path.dirname(os.path.abspath(f)) == HERE: continue
        for l in open(f):
            p = l.split()
            if p and p[0] == "B": used.add(int(p[4]))
            elif p and p[0] in ("A", "AF"): used.add(int(p[3]))
    used |= {T.run_seed(20261013, 0, m, r) for m in range(9) for r in range(4)}             # G.used_seeds()'s stage A and pilot seeds
    used |= {T.run_seed(20261091, 5, 0, r) for r in range(8)} | {9700, 9701, 9702, 9703, 9711}
    sym = os.path.join(HS, "experiments_gen3_sym_261009", "sym_seeds.txt")
    if os.path.exists(sym): used |= {int(l.split()[3]) for l in open(sym) if l.strip()}
    seeds = [int(l.split()[4]) if l.startswith("B ") else int(l.split()[3]) for l in L]
    nch = -(-len(L) // K)
    print(f"{os.path.relpath(out, HS)}: {len(L)} lines, SHA-256 {sha}; chunks of {K}: {nch} (array 1-{nch}%64); distinct seeds "
          f"{len(set(seeds))} (each used by both engines); overlap with every earlier list incl. the Mac's Test G and the symmetry test: "
          f"{len(set(seeds) & used)}")
    print("lines per block, mass and engine: " + ", ".join(f"{b} M={m} {e} {v}" for (b, m, e), v in sorted(
        Counter((l.split()[-1], int(l.split()[2]) if l.startswith("B ") else 0, l.split()[-2]) for l in L).items())))
    if ext is None:
        for M in MASSES_KOA:
            E = ext_seed_text(M); es = {int(l.split()[4]) for l in E.splitlines()}
            print(f"extension seed list at M = {M} (r = {G.NSEED}..{2 * G.NSEED - 1}, lines 'main {G.CID} M r seed'): SHA-256 "
                  f"{hashlib.sha256(E.encode()).hexdigest()}; overlap with the main list and the earlier seeds: {len(es & (set(seeds) | used))}")


if __name__ == "__main__":
    main()
