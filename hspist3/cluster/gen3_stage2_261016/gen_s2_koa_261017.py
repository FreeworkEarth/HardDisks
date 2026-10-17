#!/usr/bin/env python3
"""##CHRIS 2026-10-16 (plan-author decision 58; registration S2, notes 261012 sec. 4.9 and sec. 4.7.142): THE TASK GENERATOR OF STAGE 2
"S2" ON KOA. It writes the task list tasks_s2_koa_261017.txt (2100 lines), its SHA-256 (tasks_s2_koa_261017.sha256), the chunk table
chunks_s2_koa_261017.tsv (one line per chunk) and the design s2_design_261017.txt, and prints the design. As gen_p2_koa_261016.py, with:
CELLS (sec. 4.9): eta 0.70, 0.72, 0.74, 0.76, 0.78, 0.80 x N 900, 1600, pin turned; the isotropy check (900, 0.76) and (900, 0.78), pin 0.
  Per cell THREE box heights, Report D57-3 task 3's rule applied again here to s2_vacancy_table_261015_heights.csv (the binary's lattice fit
  replicated; every height checked through strtof): the arrangement with at most 1 vacancy (pin 0: 0 vacancies) whose heights lie nearest the
  stage 1 height (k 720 for N 900, 960 for N 1600); in it, the heights whose START strain is nearest -0.6 %, 0, +0.6 % (a tie: the k nearer the
  stage 1 height). The 42 heights are checked equal to D57-3's table (experiments_gen3_d57_261016/geometry_cost/geometry.csv) before anything
  is written. Exact box: H = k / 24 (the --height text: whole numbers plain, else 17 significant digits), L0 = N pi / (8 H eta), the divider
  at L0, both compartments alike.
MASSES alpha 1, 2, 5, 10, 20 (M = alpha N). SEEDS: 10 per (cell, H, mass), s = 0..9, by the stage 1 rule: run_seed(BASE, stream, mass index,
  s), BASE = 20261017 (sec. 4.9), one stream per (cell, H) in the order of the task list, the mass index the alpha's index in the nine-mass
  ladder (gen_stage1_koa_261012.ALPHA). Every seed is checked against the earlier task lists, the recorded commands and the extra sets of
  gen_p2_koa_261016.collisions(); seeds found there are LISTED with their origin (BASE 20261017 is also the base of D56-1's lattice starts,
  run_seed(20261017, cell 0..13, start 0..2, 0): start-only checks, 60 held steps, no trajectory) -- the generator does not change the
  registered rule; the plan author decides.
PROTOCOL (stage 1's, as P2): T_eq = 5000 (300000 held steps); 200 predicted periods of the binary's nu_pred (nu_pred_binary); the stride
  d_stride(nu_pred); Delta_psi = max(1, T_rec / 2000); virial blocks of 50; the lattice start with the pin; snapshots every (T_eq + T_rec) / 100;
  the initial dump (p2_worker_koa.sh, unchanged).
TIME (Report D57-3 task 3): events = (1.2004 (Z - 1) + 0.0112) N (T_eq + T_rec) x the pin's factor (1.002 turned, 1.171 pin 0), Z = D56-1's
  planning Z (Liu's global branch at 0.70, the solid branch from 0.72) x 1.268 from 0.72 (the upper); microseconds per event: P2's measured
  medians at N 900 (E5-2670 v3 3.555, E5-4640 v2 3.935), at N 1600 x D56-1's N 1600 / N 900 ratio (3.836 / 3.533, 4.580 / 4.242). The
  array's limit: 40 min (decision 58: 1.5 x the longest line on E5-4640 v2, rounded up to 5 min).
TASK LINE (15 fields, the stage 1 worker's): B <cell folder> <M> <s> <seed> <L0> <H> <N_s> <stride> <BASE> <hold> <Delta_psi> <T_snap> <pin> S2;
  the cell folder stage2/gen3/<cell id>/m_<M>, the cell id S2_N<N>_eta<6 decimals>_k<k, 4 digits>_pin<0|turned> (p2_reduce.py reads eta from
  its third field, the height k as its "G" and the pin from its fifth).
ORDER: N 900 turned by eta, N 900 pin 0 by eta, N 1600 turned by eta; within a cell the heights by target (-0.6, 0, +0.6 %), then alpha, then s.
CHUNKS: one line each (chunk j = line j): one array of 2100 one-core tasks.
usage (from hspist3/): python3 cluster/gen3_stage2_261016/gen_s2_koa_261017.py [--dry-run] [--out <a trial folder>]
"""
import csv, ctypes, ctypes.util, hashlib, math, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); CL = os.path.dirname(HERE); HS = os.path.dirname(CL)
for p in (HERE, os.path.join(CL, "gen3_stage1_261012"), os.path.join(HS, "validation"), HS): sys.path.insert(0, p)
import gen_stage1_koa_261012 as GEN
import tests_20260913 as T
import plot_speed_of_sound_edmd as SOS

BASE = 20261017; ALPHAS = (1, 2, 5, 10, 20); NSEED = 10; TEQ = 5000.0; HOLD = int(round(60 * TEQ))
ETAS = (0.70, 0.72, 0.74, 0.76, 0.78, 0.80); KREF = {900: 720, 1600: 960}; TARGETS = (-0.006, 0.0, 0.006)
CELLS = [(900, eta, "turned") for eta in ETAS] + [(900, eta, "0") for eta in (0.76, 0.78)] + [(1600, eta, "turned") for eta in ETAS]
RATE_A, RATE_B, ZX = 1.2004, 0.0112, 1.268; MEAS = {"0": 1.171, "turned": 1.002}
US = {900: {"v3": 3.555, "v2": 3.935}, 1600: {"v3": 3.555 * 3.836 / 3.533, "v2": 3.935 * 4.580 / 4.242}}; LIMIT_MIN = 40
MIX_US = {900: 2.853, 1600: 3.097}                                  # Report D57-3 task 3: P2's KOA mix, microseconds per event
VAC = os.path.join(HERE, "s2_vacancy_table_261015_heights.csv")
D57 = os.path.join(HS, "experiments_gen3_d57_261016", "geometry_cost", "geometry.csv")
OUTD = os.path.abspath(sys.argv[sys.argv.index("--out") + 1]) if "--out" in sys.argv else HERE
TASKS = os.path.join(OUTD, "tasks_s2_koa_261017.txt"); SHAF = os.path.join(OUTD, "tasks_s2_koa_261017.sha256")
CHUNKF = os.path.join(OUTD, "chunks_s2_koa_261017.tsv"); DESIGN = os.path.join(OUTD, "s2_design_261017.txt")
_libc = ctypes.CDLL(ctypes.util.find_library("c")); _libc.strtof.restype = ctypes.c_float
_libc.strtof.argtypes = [ctypes.c_char_p, ctypes.POINTER(ctypes.c_char_p)]


def height_px(text):
    """the binary's box height in px for the --height text: (int)(strtof(text) * 24.0f), in float arithmetic (gen_p2_koa_261016's)"""
    v = np.float32(_libc.strtof(text.encode(), None)); return int(np.float32(v * np.float32(24.0)))


def htext(k): return f"{k / 24.0:.17g}"


def z_eos(eta):
    a = np.array([eta]); return float(SOS.Z_liu_global(a)[0]) if eta < 0.72 - 1e-12 else float(SOS.Z_liu_solid(a)[0])


def pick(rows, N, eta, pin):
    q = [r for r in rows if int(r["N"]) == N and abs(float(r["eta"]) - eta) < 1e-9 and r["pin"] == pin and r["height_px_ok"] == "1"
         and 0 <= int(r["vacancies"]) <= (1 if pin == "turned" else 0)]
    arrs = {}
    for r in q: arrs.setdefault(r["arrangement"], []).append(r)
    dist = {a: min(abs(int(r["k_px"]) - KREF[N]) for r in g) for a, g in arrs.items()}; best = min(dist.values())
    win = [a for a, d in dist.items() if d == best]
    if len(win) != 1: raise SystemExit(f"STOP: {len(win)} arrangements equally near the stage 1 height at N {N} eta {eta} pin {pin}")
    g = arrs[win[0]]; out = []
    for t in TARGETS:
        r = min(g, key=lambda r: (abs(float(r["e"]) - t), abs(int(r["k_px"]) - KREF[N])))
        out.append(dict(N=N, eta=eta, pin=pin, arrangement=win[0], vacancies=int(r["vacancies"]), target=t, k=int(r["k_px"]), e=float(r["e"]),
                        gap=float(r["gap_sigma"]), dist=best))
    if len({o["k"] for o in out}) != 3: raise SystemExit(f"STOP: two targets share a height at N {N} eta {eta} pin {pin}")
    return out


def geometry():
    rows = list(csv.DictReader(open(VAC))); G = []
    for N, eta, pin in CELLS: G += pick(rows, N, eta, pin)
    ref = {(int(r["N"]), round(float(r["eta"]), 2), r["pin"], round(float(r["target_e"]), 4)): int(r["k"]) for r in csv.DictReader(open(D57))}
    if len(ref) != len(G): raise SystemExit(f"STOP: D57-3's table has {len(ref)} heights, this rule {len(G)}")
    for g in G:
        if ref.get((g["N"], round(g["eta"], 2), g["pin"], round(g["target"], 4))) != g["k"]:
            raise SystemExit(f"STOP: the height differs from D57-3's table: {g}")
    for g in G:
        g["H"] = g["k"] / 24.0; g["htxt"] = htext(g["k"]); g["L0"] = g["N"] * math.pi / (8.0 * g["H"] * g["eta"])
        if height_px(g["htxt"]) != g["k"]: raise SystemExit(f"STOP: --height={g['htxt']} gives {height_px(g['htxt'])} px, not {g['k']}")
        g["cid"] = f"S2_N{g['N']:04d}_eta{g['eta']:.6f}_k{g['k']:04d}_pin{g['pin']}"; g["aspect"] = (g["L0"] - 0.025) / g["H"]
    return G


def plan():
    G = geometry(); rows, lines = [], []
    for stream, g in enumerate(G):
        g["stream"] = stream; zl = z_eos(g["eta"]); zu = zl * (ZX if g["eta"] >= 0.72 - 1e-12 else 1.0); r = RATE_A * (zu - 1.0) + RATE_B
        for a in ALPHAS:
            M = int(round(a * g["N"])); mi = GEN.ALPHA.index(a)
            nup = GEN.nu_pred_binary(g["N"], g["H"], g["L0"], M); steps = math.ceil((200.0 / nup) / GEN.DTF); Trec = steps * GEN.DTF
            stride = T.d_stride(nup); dpsi = max(1.0, Trec / 2000.0); tsnap = (TEQ + Trec) / 100.0
            ev = r * g["N"] * (TEQ + Trec) * MEAS[g["pin"]]; rel = f"stage2/gen3/{g['cid']}/m_{M}"
            seeds = [T.run_seed(BASE, stream, mi, s) for s in range(NSEED)]
            for s, sd in enumerate(seeds):
                lines.append(f"B {rel} {M} {s} {sd} {g['L0']:.17g} {g['htxt']} {g['N'] // 2} {stride} {BASE} {HOLD} {dpsi:.17g} {tsnap:.17g} {g['pin']} S2\n")
            rows.append(dict(cid=g["cid"], N=g["N"], eta=g["eta"], pin=g["pin"], k=g["k"], a=a, M=M, mi=mi, stream=stream, nup=nup, Trec=Trec,
                             stride=stride, dpsi=dpsi, tsnap=tsnap, Zl=zl, Zu=zu, ev=ev, min_v3=ev * US[g["N"]]["v3"] / 6e7,
                             min_v2=ev * US[g["N"]]["v2"] / 6e7, seeds=seeds))
    return G, rows, lines


def collisions(seeds):
    """as gen_p2_koa_261016.collisions(), plus P2's own list (scanned by used_seeds) and the D57-3 Mac run (P2's line 21: in P2's list);
    returns the seeds found with their origin"""
    used, nlists, kinds, cmds = GEN.used_seeds()
    st1 = {int(l.split()[4]) for l in open(GEN.TASKS) if l.startswith("B ")}
    d56p = {T.run_seed(20261015, 0, 4, s) for s in range(2)}
    d56b = {T.run_seed(20261017, ci, si, 0): (ci, si) for ci in range(14) for si in range(3)}
    hit = {}
    for s in seeds:
        if s in d56b: hit[s] = f"D56-1's lattice start (cell {d56b[s][0]}, start {d56b[s][1]}): run_seed(20261017, {d56b[s][0]}, {d56b[s][1]}, 0)"
        elif s in d56p: hit[s] = "D56-1's pilot"
        elif s in st1: hit[s] = "the stage 1 list"
        elif s in used: hit[s] = "an earlier task list or recorded command"
    return hit, len(used | st1 | d56p | set(d56b)), nlists + 1, cmds


def main():
    G, rows, lines = plan()
    L = ["# STAGE 2 \"S2\" (registration sec. 4.9; decision 58): the design", "",
         "## The heights (42: 36 turned, 6 pin 0; Report D57-3 task 3's rule, checked equal to its table)", "",
         "| stream | cell id | N | eta | pin | arrangement | vacancies | target e | H [sigma] (k) | L0 | aspect W/H | start e | gap [sigma] |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for g in G:
        L.append(f"| {g['stream']} | {g['cid']} | {g['N']} | {g['eta']:.2f} | {g['pin']} | {g['arrangement']} | {g['vacancies']} | {100 * g['target']:+.1f} % | "
                 f"{g['htxt']} ({g['k']}) | {g['L0']:.6f} | {g['aspect']:.4f} | {100 * g['e']:+.3f} % | {g['gap']:.4f} |")
    L += ["", "## Per (cell, H, mass): nu_pred, the record, the stride, the cadences, the time (the upper Z)", "",
          "| cell id | alpha | M | nu_pred | T_rec | stride | Delta_psi | T_snap | Z upper | events [M] | min (v3) | min (v2) |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        L.append(f"| {r['cid']} | {r['a']} | {r['M']} | {r['nup']:.6f} | {r['Trec']:.1f} | {r['stride']} | {r['dpsi']:.4f} | {r['tsnap']:.3f} | "
                 f"{r['Zu']:.3f} | {r['ev'] / 1e6:.1f} | {r['min_v3']:.2f} | {r['min_v2']:.2f} |")
    longest = max(r["min_v2"] for r in rows); lim = int(math.ceil(1.5 * longest / 5.0) * 5)
    if lim > LIMIT_MIN: raise SystemExit(f"STOP: the longest line needs {lim} min, above decision 58's {LIMIT_MIN}")
    ch3 = sum(r["min_v3"] for r in rows) * NSEED / 60.0; ch2 = sum(r["min_v2"] for r in rows) * NSEED / 60.0
    chm = sum(r["ev"] * MIX_US[r["N"]] / 6e7 for r in rows) * NSEED / 60.0
    seeds = [s for r in rows for s in r["seeds"]]; hit, nused, nl, ncmd = collisions(seeds)
    L += ["", f"- task lines: {len(lines)} ({len(G)} heights x 5 masses x 10 seeds); distinct seeds: {len(set(seeds))}; seeds found in an earlier "
              f"list, the stage 1 list, D56-1's pilot or lattice starts, or a recorded command ({nused} seeds in {nl} lists, {ncmd} command files): {len(hit)}"]
    for s, why in sorted(hit.items(), key=lambda x: x[1]):
        line = next(l for l in lines if l.split()[4] == str(s)).split()
        L.append(f"  - seed {s} (task {line[1]}, s {line[3]}): {why}")
    L += [f"- predicted core-hours (the upper Z): {ch3:.1f} on E5-2670 v3, {ch2:.1f} on E5-4640 v2; the longest line {longest:.1f} min on v2 "
          f"-> 1.5 x, rounded up to 5 min: {lim} min (the array's limit {LIMIT_MIN} min, decision 58)",
          f"- wall time at 64 cores (v3): about {ch3 * 60 / 64:.0f} min; at P2's KOA mix (Report D57-3: {MIX_US[900]} us per event at N 900, "
          f"{MIX_US[1600]} at N 1600): {chm:.0f} core-hours, about {chm / 64:.1f} h at 64 cores"]
    if len(set(seeds)) != len(seeds): raise SystemExit("\n".join(L) + "\nSTOP: a seed repeats within S2")
    if any(not w.startswith("D56-1's lattice start") for w in hit.values()): raise SystemExit("\n".join(L) + "\nSTOP: a seed was used by a trajectory before")
    txt = "".join(lines); sha = hashlib.sha256(txt.encode()).hexdigest()
    ctxt = "chunk\tfirst\tlast\tlines\tgroup\tcell\tM\tpredicted_s_v2\tlimit_min\tclass_min\n"
    by = {(r["cid"], r["M"]): r for r in rows}
    for j, l in enumerate(lines, 1):
        f = l.split(); cid = f[1].split("/")[2]; r = by[(cid, int(f[2]))]
        ctxt += f"{j}\t{j}\t{j}\t1\tS2\t{cid}\t{r['M']}\t{60 * r['min_v2']:.1f}\t{1.5 * r['min_v2']:.1f}\t{LIMIT_MIN}\n"
    L += [f"- {os.path.relpath(TASKS, HS)}: SHA-256 {sha}; {os.path.relpath(CHUNKF, HS)}: {len(lines)} chunks"]
    print("\n".join(L))
    if "--dry-run" in sys.argv: print("\n--dry-run: nothing written"); return
    for p in (TASKS, SHAF, CHUNKF, DESIGN):
        if os.path.exists(p): raise SystemExit(f"REFUSED: {p} exists; the package is written once")
    open(TASKS, "w").write(txt); open(SHAF, "w").write(sha + "\n"); open(CHUNKF, "w").write(ctxt); open(DESIGN, "w").write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
