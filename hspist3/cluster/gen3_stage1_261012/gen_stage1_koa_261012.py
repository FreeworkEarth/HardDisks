#!/usr/bin/env python3
"""##CHRIS 2026-10-12 (plan-author decision 39 item 5 f; the registration: 261012 sec. 4.8, main 85cb671e; gate G1): THE STAGE 1 TASK
GENERATOR, Part B (blocks B, D and S). This script, the task list it writes (tasks_stage1_koa_261012.txt), its SHA-256
(tasks_stage1_koa_261012.sha256) and the chunk table (chunks_stage1_koa_261012.tsv) are committed together before any stage 1 trajectory
exists; none of them is edited afterwards. Every number below is the registration's (sec. 4.8.2-4.8.3, 4.8.8); where the registration
leaves a choice to CC, it is named "CHOICE" here and in --design.
CELLS (sec. 4.8.2): H = 10 sqrt(N/100), L0 = N pi / (8 H eta), passed with 17 significant digits (--gen3-exact-box); N_s = N/2.
  B: N in {100, 400, 900, 1600} x the 16 densities (0.05 ... 0.60 fluid set Phi; 0.63, 0.66, 0.68, 0.69 approach cells);
  D: eta in {0.02, 0.035, 0.075} x N in {100, 400};   S: eta in {0.60, 0.63, 0.66, 0.68, 0.69} x all four N x two starts (0, turned).
MASSES: alpha = M/N in {0.5, 1, 2, 3, 5, 7.5, 10, 15, 20} (B, D); alpha in {1, 10} (S).   SEEDS: 25 per (cell, mass), s = 0..24 (B, D);
  10 per (cell, start, mass) (S). run_seed(BASE, stream, mass index, s) (the driver's speed_sound_run_seed, tests_20260913.run_seed): one
  stream per (block, N, eta) and in block S per start, numbered in the order of the task list (CHOICE); the mass index is the alpha's
  index in the nine-mass ladder (CHOICE); BASE = 20261311, the registration's proposal (--check-base: as an integer token it may occur
  only in the registration text and this package).
PROTOCOL (sec. 4.8.3): Test G's B command with the changes of sec. 4.8.3.2, run by s1_worker_koa.sh: T_eq = 1000 sigma-time (eta <= 0.60)
  or 5000 (eta >= 0.63) = 60 T_eq hold steps; 200 predicted periods recorded (--target-oscillations=200: the binary's own nu_pred,
  Henderson/Roman, sets the record T_rec = 200 / nu_pred, replicated here exactly); Delta_psi = max(1, T_rec / 2000); virial blocks of 50;
  snapshots every (T_eq + T_rec)/100 on s = 0, 1 of alpha = 5 in B and D, of alpha = 1 in S; the trace stride by the generators' rule of
  32 samples per predicted period, d_stride(nu_KR) with nu_KR = c_s^KR(eta) K(M/N) / (2 pi (L0 - 1.025)) (CHOICE: Test G's x_of took
  N_s = 50 from tests_20260913; here N_s = N/2 of the cell).
TASK LINE (15 fields; s1_worker_koa.sh reads them): B <cell folder> <M> <s> <seed> <L0> <H> <N_s> <stride> <BASE> <hold steps> <Delta_psi>
  <T_snap or 0> <start: default | 0 | turned> <block>; the cell folder is stage1/gen3/<cell id>/m_<M>, cell id <block>_N<NNNN>_eta<6 decimals>
  (S: _start0 or _turned appended).
ORDER (sec. 4.8.8): B at N = 100, B at N = 400, D, B at N = 900, S, B at N = 1600; within a block by eta, then alpha, then s.
CHUNKS (sec. 4.8.8): consecutive lines of one (block, N, eta, mass), predicted time <= 60 min at the old-node planning rate (the 25
  or 10 seeds split into the fewest chunks of near-equal size); a trajectory longer than that is a chunk of its own. Time limit:
  max(30 min, 3 x the predicted chunk time), at most 24 h, rounded UP to the next class 0:30, 1, 2, 3, 4, 6, 8, 12, 16, 24 h
  (CHOICE: one array has one time limit, so the arrays are formed per block and class).
COST (sec. 4.8.8's planning model, here with the project's KR module): events per sigma-time = 1.128 N (Z_KR - 1) + 1.2 N; the
  whole-run event rate of K2-2 (Gold node, sec. 4.7.49 item 5) per N at pi/8 and 0.70, linear in eta between them and constant below
  pi/8; the old-node rate = that / 1.5. The threshold of sec. 4.8.9 item 3: B + D + S <= 5,000 core-hours before the factor 1.5.
  Checked against measured E5 times: the preview (N 100, seven densities, 1575 runs on E5-2670 v3; its own protocol).
usage (from hspist3/): python3 cluster/gen3_stage1_261012/gen_stage1_koa_261012.py [--design] [--check-base] [--dry-run]
  --dry-run   print everything, write nothing
"""
import collections, glob, hashlib, math, os, re, subprocess, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); CL = os.path.dirname(HERE); HS = os.path.dirname(CL); REPO = os.path.dirname(HS)
MAIN_WT = os.path.join(os.path.dirname(REPO), "HardDisks")
sys.path.insert(0, os.path.join(HS, "validation")); sys.path.insert(0, HS)
import tests_20260913 as T
import plot_speed_of_sound_edmd as SOS

BASE = 20261311
ETA_B = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, math.pi / 8, 0.45, 0.50, 0.55, 0.60, 0.63, 0.66, 0.68, 0.69]
ETA_D = [0.02, 0.035, 0.075]; ETA_S = [0.60, 0.63, 0.66, 0.68, 0.69]
ALPHA = [0.5, 1, 2, 3, 5, 7.5, 10, 15, 20]; ALPHA_S = [1, 10]
NS_B, NS_D, NS_S = (100, 400, 900, 1600), (100, 400), (100, 400, 900, 1600)
GROUPS = ["B100", "B400", "D", "B900", "S", "B1600"]                  # the registered order of the arrays (sec. 4.8.8)
DTF = float(np.float32(0.4)) / 24.0                                      # the binary's sample dt: (double) 0.4f / 24
K2_RATE = {100: (7.02e5, 5.79e5), 400: (5.15e5, 4.70e5), 900: (3.74e5, 3.91e5), 1600: (2.61e5, 3.34e5)}   # K2-2 events/s (run), pi/8, 0.70
CLASSES = [30, 60, 120, 180, 240, 360, 480, 720, 960, 1440]             # minutes
TASKS = os.path.join(HERE, "tasks_stage1_koa_261012.txt"); SHAF = os.path.join(HERE, "tasks_stage1_koa_261012.sha256")
CHUNKF = os.path.join(HERE, "chunks_stage1_koa_261012.tsv")


def nu_pred_binary(N, H, L0, M):
    """speed_sound_henderson_frequency_estimate as the B loop calls it under --gen3-exact-box (f72d764 00ALLINONE.c:16019-16039,
    :17051-17065): eta from L0 exactly, L_eff = (float) L0 - 1, K from cot K = alpha K with alpha = M / (2 N_s)"""
    eta = N * math.pi * 0.25 / (2.0 * L0 * H)
    a = 0.125; num = 1.0 + eta + 3.0 * a * eta * eta - a * eta ** 3
    return math.sqrt(2.0 * num / (1.0 - eta) ** 3) * T.k_root(M / N) / (2.0 * math.pi * (float(np.float32(L0)) - 1.0))


def z_kr(eta): return float(SOS.Z_kolafa_rottner_2006(np.array([eta]))[0])


def rate_old(N, eta):
    a, b = K2_RATE[N]; e0, e1 = math.pi / 8, 0.70
    r = a if eta <= e0 else a + (b - a) * (min(eta, e1) - e0) / (e1 - e0)
    return r / 1.5


def events_per_sigma(N, eta): return 1.128 * N * (z_kr(eta) - 1.0) + 1.2 * N


def cells():
    """the registered cells in the registered order: (group, block, N, eta, start, alphas, seeds)"""
    out = []
    for g in GROUPS:
        if g.startswith("B"):
            N = int(g[1:]); out += [(g, "B", N, e, "default", ALPHA, 25) for e in ETA_B]
        elif g == "D":
            out += [(g, "D", N, e, "default", ALPHA, 25) for N in NS_D for e in ETA_D]
        else:
            out += [(g, "S", N, e, st, ALPHA_S, 10) for N in NS_S for e in ETA_S for st in ("0", "turned")]
    return out


def cell_id(blk, N, eta, start):
    return f"{blk}_N{N:04d}_eta{eta:.6f}" + ("" if blk != "S" else ("_start0" if start == "0" else "_turned"))


def plan():
    """per (cell, mass): every quantity of the design; the task lines in order"""
    rows, lines = [], []
    for stream, (g, blk, N, eta, start, alphas, nseed) in enumerate(cells()):
        H = 10.0 * math.sqrt(N / 100.0); L0 = N * math.pi / (8.0 * H * eta); Ns = N // 2; cid = cell_id(blk, N, eta, start)
        Teq = 1000.0 if eta <= 0.60 + 1e-12 else 5000.0; hold = int(round(60 * Teq))
        for a in alphas:
            M = int(round(a * N)); mi = ALPHA.index(a)
            nup = nu_pred_binary(N, H, L0, M); steps = math.ceil((200.0 / nup) / DTF); Trec = steps * DTF
            nukr = T.kr_cs(eta) * T.k_root(M / N) / (2.0 * math.pi * (L0 - 1.025)); stride = T.d_stride(nukr)
            dpsi = max(1.0, Trec / 2000.0); tsnap = (Teq + Trec) / 100.0
            snap_a = 1 if blk == "S" else 5
            ev = events_per_sigma(N, eta) * (Teq + Trec); sec = ev / rate_old(N, eta)
            rel = f"stage1/gen3/{cid}/m_{M}"
            seeds = [T.run_seed(BASE, stream, mi, s) for s in range(nseed)]
            for s, sd in enumerate(seeds):
                ts = tsnap if (a == snap_a and s in (0, 1)) else 0.0
                lines.append((len(rows), f"B {rel} {M} {s} {sd} {L0:.17g} {H:.17g} {Ns} {stride} {BASE} {hold} {dpsi:.17g} "
                                         f"{ts:.17g} {start} {blk}\n"))
            rows.append(dict(g=g, blk=blk, N=N, eta=eta, start=start, cid=cid, a=a, M=M, H=H, L0=L0, Ns=Ns, Teq=Teq, hold=hold, nup=nup,
                             nukr=nukr, steps=steps, Trec=Trec, periods_kr=Trec * nukr, stride=stride, dpsi=dpsi, tsnap=tsnap,
                             nseed=nseed, seeds=seeds, ev=ev, sec=sec, stream=stream, mi=mi, snaps=2 if a == snap_a else 0))
    return rows, lines


def chunk(rows, lines):
    """chunks over consecutive lines of one (cell, mass); returns the chunk table rows"""
    out, i = [], 0
    by = collections.defaultdict(list)
    for k, (ri, l) in enumerate(lines): by[ri].append(k)
    for ri, r in enumerate(rows):
        idx = by[ri]; t = r["sec"]
        n = len(idx) if t > 3600 else max(1, math.ceil(len(idx) * t / 3600.0))
        base, extra = divmod(len(idx), n); p = 0
        for c in range(n):
            sz = base + (1 if c < extra else 0); seg = idx[p:p + sz]; p += sz
            pred = sz * t; lim = min(1440.0, max(30.0, 3 * pred / 60.0)); cls = next(x for x in CLASSES if x >= lim - 1e-9)
            out.append(dict(j=len(out) + 1, first=seg[0] + 1, last=seg[-1] + 1, lines=sz, g=r["g"], cid=r["cid"], M=r["M"],
                            pred_s=pred, limit_min=lim, class_min=cls))
    return out


def used_seeds():
    """every exact seed of every earlier task list (both trees), the extra sets Test R's generator named, K2-2's and the M5 job's,
    and every --speed-sound-exact-seed / --seed of a recorded command (both trees: the pilots P1-P6, the M-series, the window tests)"""
    used, lists, kinds, cmds = set(), set(), collections.Counter(), 0
    for t in (HS, os.path.join(MAIN_WT, "hspist3")):
        for f in glob.glob(os.path.join(t, "cluster", "*", "tasks_*.txt")):
            if os.path.dirname(os.path.abspath(f)) == HERE: continue
            lists.add(os.path.relpath(f, t))
            for l in open(f):
                p = l.split()
                if not p: continue
                kinds[p[0]] += 1
                if p[0] == "B": used.add(int(p[4]))
                elif p[0] in ("A", "AF"): used.add(int(p[3]))
    used |= {T.run_seed(20261013, 0, m, r) for m in range(9) for r in range(4)}
    used |= {T.run_seed(20261091, 5, 0, r) for r in range(8)} | {9700, 9701, 9702, 9703, 9711}
    used |= {T.run_seed(20261219, 0, k, 0) for k in range(12)} | {T.run_seed(20261219, 1, n, r) for n in (100, 400) for r in range(4)}
    used |= {T.run_seed(20261219, 2, 400, r) for r in range(4)}
    for t in (HS, os.path.join(MAIN_WT, "hspist3")):
        sym = os.path.join(t, "experiments_gen3_sym_261009", "sym_seeds.txt")
        if os.path.exists(sym): used |= {int(l.split()[3]) for l in open(sym) if l.strip()}
        for f in glob.glob(os.path.join(t, "**", "command.txt"), recursive=True):
            s = open(f, errors="replace").read(); cmds += 1
            used |= {int(x) for x in re.findall(r"--speed-sound-exact-seed=(\d+)", s)} | {int(x) for x in re.findall(r"--seed=(\d+)", s)}
    return used, len(lists), kinds, cmds


def check_base():
    """BASE as an integer token: both working trees (rg) and the tracked files of three heads (git grep); this package excluded"""
    pat = rf"(^|[^0-9.]){BASE}([^0-9.]|$)"
    rg = subprocess.run(["rg", "-l", "--no-ignore", "--hidden", "-g", "!.git", "-e", pat, REPO, MAIN_WT], capture_output=True, text=True)
    hits = sorted({h for h in rg.stdout.splitlines() if "/cluster/gen3_stage1_261012/" not in h})
    gg = subprocess.run(["git", "-C", REPO, "grep", "-E", "-l", pat, "main", "engine-gen3", "engine-divider-resched", "--"],
                        capture_output=True, text=True)
    ghits = sorted({h for h in gg.stdout.splitlines() if "/cluster/gen3_stage1_261012/" not in h})
    return hits, ghits


def preview_check():
    """the cost model against the preview's measured E5 seconds (N 100; the A1 v2 command: 2000 steps held, 200 predicted periods)"""
    import pandas as pd
    root = os.path.join(MAIN_WT, "hspist3", "experiments_gen3_preview_koa_261011", "merged", "preview", "gen3")
    out = []
    for c in sorted(glob.glob(os.path.join(root, "eta_*"))):
        for m in sorted(glob.glob(os.path.join(c, "m_*")), key=lambda p: int(p.rsplit("_", 1)[1])):
            h = pd.read_csv(os.path.join(m, "runs_health.csv")); red = pd.read_csv(os.path.join(m, "red_nu.csv"))
            log = open(os.path.join(m, "run.log"), errors="replace").read()
            eta = float(re.search(r" eta=(\S+) ", re.search(r"^\[EDMD3-RUN\] .*$", log, re.M).group(0)).group(1))
            nup = float(red["nu_pred"].iloc[0]); Tsim = 2000 * DTF + math.ceil((200.0 / nup) / DTF) * DTF
            pred = events_per_sigma(100, eta) * Tsim / rate_old(100, eta)
            out.append((os.path.basename(c), int(m.rsplit("_", 1)[1]), eta, Tsim, float(h["seconds"].mean()), pred, len(h)))
    return out


def design(rows, lines, ch, used, nlists, kinds, ncmd, base_hits):
    seeds = [int(l.split()[4]) for _, l in lines]
    print("# STAGE 1, PART B -- the design tables of the task generator (gate G1; decision 39 item 5 f; registration 261012 sec. 4.8)\n")
    print(f"base seed {BASE}; {len(lines)} task lines; {len(rows)} (cell, mass) pairs; {len(ch)} chunks\n")
    print("## 1. Trajectories per block (sec. 4.8.2 and its cost table)\n")
    print("| group | cells | (cell, mass) pairs | trajectories |\n|---|---|---|---|")
    for g in GROUPS:
        rr = [r for r in rows if r["g"] == g]
        print(f"| {g} | {len({r['cid'] for r in rr})} | {len(rr)} | {sum(r['nseed'] for r in rr)} |")
    print(f"| **all** | {len({r['cid'] for r in rows})} | {len(rows)} | {sum(r['nseed'] for r in rows)} |")
    print("\n## 2. The cells (per cell; the masses' values in section 3)\n")
    print("| cell | block | N | eta (exact) | H | L0 (17 significant digits) | N_s | T_eq | hold steps | start |\n|---|---|---|---|---|---|---|---|---|---|")
    seen = set()
    for r in rows:
        if r["cid"] in seen: continue
        seen.add(r["cid"])
        print(f"| {r['cid']} | {r['blk']} | {r['N']} | {r['eta']:.17g} | {r['H']:g} | {r['L0']:.17g} | {r['Ns']} | {r['Teq']:g} | {r['hold']} | {r['start']} |")
    print("\n## 3. Per (cell, mass): the binary's nu_pred and record, the KR periods, Delta_psi, the stride, snapshots, events, cost\n")
    print("| cell | M | alpha | nu_pred (binary) | record steps | T_rec [sigma-time] | KR periods in T_rec | Delta_psi | stride | T_snap (snapshot seeds) | events per trajectory | E5 core-min per trajectory |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        snap = f"{r['tsnap']:.6g} (s = 0, 1)" if r["snaps"] else "-"
        print(f"| {r['cid']} | {r['M']} | {r['a']:g} | {r['nup']:.8g} | {r['steps']} | {r['Trec']:.6g} | {r['periods_kr']:.1f} | {r['dpsi']:.6g} | "
              f"{r['stride']} | {snap} | {r['ev']:.4g} | {r['sec'] / 60:.2f} |")
    print("\n## 4. Seeds\n")
    nd = len(set(seeds)); ov = set(seeds) & used
    print(f"- {len(seeds)} task lines, {nd} distinct exact seeds: {'ALL DISTINCT' if nd == len(seeds) else 'NOT DISTINCT'}")
    print(f"- earlier seeds checked: {len(used)} (every task list of cluster/*/ in both working trees: {nlists} lists, line types {dict(kinds)}; "
          f"the extra sets of the Test G and Test R generators; K2-2's and the M5 job's; the symmetry test's; every --speed-sound-exact-seed "
          f"and --seed of {ncmd} recorded command.txt files in both trees); overlap: {len(ov)}" + (" -- ZERO OVERLAP" if not ov else " -- OVERLAP"))
    print(f"- streams: {len(cells())} (one per block, N, eta and, in block S, start; numbered in the task list's order); mass index = the alpha's "
          f"index in the nine-mass ladder")
    if base_hits is not None:
        h, gh = base_hits
        print(f"- the base {BASE} as an integer token, outside this package: {len(h)} files in both working trees, {len(gh)} tracked files of main, "
              f"engine-gen3, engine-divider-resched:")
        for x in sorted(set(os.path.basename(p) for p in h) | set(p.split(':', 1)[1] if ':' in p else p for p in gh))[:20]: print(f"    {x}")
    print("\n## 5. Cost (sec. 4.8.8's model and columns, with the project's KR module: core-h on the old E5 nodes = events / (K2-2's Gold "
          "rate / 1.5); the column x 1.5 is the registration's margin)\n")
    REG = {"B100": 13, "B400": 120, "D": 149, "B900": 526, "S": 66, "B1600": 1680}      # sec. 4.8.8's planning table, core-h
    print("| group | trajectories | core-h | x 1.5 | wall at 64 cores | longest trajectory | events | chunks | registration's planning core-h |\n|---|---|---|---|---|---|---|---|---|")
    tot = 0.0
    for g in GROUPS:
        rr = [r for r in rows if r["g"] == g]; e5 = sum(r["sec"] * r["nseed"] for r in rr) / 3600.0
        lo = max(rr, key=lambda r: r["sec"]); tot += e5
        print(f"| {g} | {sum(r['nseed'] for r in rr)} | {e5:,.0f} | {1.5 * e5:,.0f} | {e5 / 64:.1f} h | {lo['sec'] / 60:.0f} min ({lo['cid']}, alpha {lo['a']:g}) | "
              f"{sum(r['ev'] * r['nseed'] for r in rr):.3g} | {sum(1 for c in ch if c['g'] == g)} | {REG[g]:,} |")
    print(f"| **B + D + S** | {sum(r['nseed'] for r in rows)} | **{tot:,.0f}** | **{1.5 * tot:,.0f}** | {tot / 64:.0f} h | | {sum(r['ev'] * r['nseed'] for r in rows):.3g} | {len(ch)} | {sum(REG.values()):,} |")
    print(f"\nTHRESHOLD (sec. 4.8.9 item 3: the script's total for B + D + S before the factor 1.5; decision 39 item 5 f): {tot:,.0f} core-h "
          f"{'<=' if tot <= 5000 else '>'} 5,000: {'BELOW -- the build goes on' if tot <= 5000 else 'EXCEEDED -- STOP AND REPORT'} "
          f"(with the factor 1.5: {1.5 * tot:,.0f})")
    print("model inputs: events per sigma-time 1.128 N (Z_KR - 1) + 1.2 N; K2-2's whole-run events/s on the Gold node "
          + ", ".join(f"N {n}: {a:.3g} (pi/8), {b:.3g} (0.70)" for n, (a, b) in K2_RATE.items()) + "; linear in eta between pi/8 and 0.70, "
          "constant below pi/8; the old nodes at 1/1.5 of it")
    pc = preview_check()
    rat = [p[4] / p[5] for p in pc]
    print(f"\nthe model against measured E5 times, the preview (N 100, its own protocol: 2000 steps held, 200 predicted periods; {sum(p[6] for p in pc)} runs "
          f"on E5-2670 v3): measured / model per (cell, mass) from {min(rat):.2f} to {max(rat):.2f}, median {sorted(rat)[len(rat) // 2]:.2f}\n")
    print("| preview cell | M | eta (recorded) | simulated sigma-time | measured s (mean) | model s | measured / model |\n|---|---|---|---|---|---|---|")
    for c, M, eta, Ts, ms, pr, n in pc:
        if M in (50, 300, 2000): print(f"| {c} | {M} | {eta:.6f} | {Ts:.0f} | {ms:.1f} | {pr:.1f} | {ms / pr:.2f} |")
    print("\n## 6. Chunks and arrays (sec. 4.8.8)\n")
    print("| group | chunks | lines | predicted chunk time (E5) mean / max | time-limit classes (chunks) | requested core-h |\n|---|---|---|---|---|---|")
    for g in GROUPS:
        cc = [c for c in ch if c["g"] == g]; cls = collections.Counter(c["class_min"] for c in cc)
        print(f"| {g} | {len(cc)} | {sum(c['lines'] for c in cc)} | {sum(c['pred_s'] for c in cc) / len(cc) / 60:.1f} / {max(c['pred_s'] for c in cc) / 60:.1f} min | "
              + ", ".join(f"{m // 60}:{m % 60:02d} ({n})" for m, n in sorted(cls.items())) + f" | {sum(c['class_min'] for c in cc) / 60:,.0f} |")
    print(f"\nchunk 1: lines {ch[0]['first']}-{ch[0]['last']} ({ch[0]['cid']}, M {ch[0]['M']}), predicted {ch[0]['pred_s'] / 60:.1f} min, time limit "
          f"{ch[0]['class_min']} min; the arrays: one per (group, class), indices = the chunk numbers, throttled so that at most 64 run at "
          f"once over all stage 1 arrays (sec. 4.8.8)")
    print(f"largest array index {len(ch)} (MaxArraySize 25001)")
    print("\n## 7. Data on KOA scratch, expected (per trajectory: the trace about 126 bytes per row, psi6(t) about 130 bytes per row, the "
          "virial blocks about 230 bytes per row, the snapshots 16 N bytes each; the merge compresses the traces about 8x)\n")
    print("| group | traces | psi6(t) | virial | snapshots | run folders (about 45 kB) | before the merge | after it |\n|---|---|---|---|---|---|---|---|")
    for g in GROUPS:
        rr = [r for r in rows if r["g"] == g]
        tr = sum(r["nseed"] * (r["steps"] / r["stride"]) * 126 for r in rr) / 1e9
        ps = sum(r["nseed"] * ((r["Teq"] + r["Trec"]) / r["dpsi"]) * 130 for r in rr) / 1e9
        vb = sum(r["nseed"] * ((r["Teq"] + r["Trec"]) / 50) * 230 for r in rr) / 1e9
        sn = sum(r["snaps"] * 101 * 16 * r["N"] for r in rr) / 1e9; rf = sum(r["nseed"] for r in rr) * 45e3 / 1e9
        print(f"| {g} | {tr:.2f} GB | {ps:.2f} GB | {vb:.3f} GB | {sn:.3f} GB | {rf:.2f} GB | {tr + ps + vb + sn + rf:.2f} GB | {tr / 8 + ps / 3 + vb + sn + rf:.2f} GB |")
    print("\n## 8. Choices left to CC by the registration, as made here\n")
    for c in ("streams numbered in the task list's order; the mass index = the alpha's index in the nine-mass ladder",
              "the stride's KR frequency with N_s = N/2 of the cell (Test G's x_of assumed N_s = 50)",
              "the time limit rounded up to a class (one array has one time limit); the arrays formed per group and class",
              "the task line's 15 fields and the folder layout stage1/gen3/<cell id>/m_<M> (Test R's merge layout, the engine level gen3)",
              "the record steps and T_rec of a (cell, mass) replicate the binary's own --target-oscillations arithmetic; Delta_psi and T_snap use that T_rec"):
        print(f"- {c}")


def main():
    rows, lines = plan(); txt = "".join(l for _, l in lines); sha = hashlib.sha256(txt.encode()).hexdigest()
    ch = chunk(rows, lines)
    ctxt = "chunk\tfirst\tlast\tlines\tgroup\tcell\tM\tpredicted_s_E5\tlimit_min\tclass_min\n" + "".join(
        f"{c['j']}\t{c['first']}\t{c['last']}\t{c['lines']}\t{c['g']}\t{c['cid']}\t{c['M']}\t{c['pred_s']:.1f}\t{c['limit_min']:.1f}\t{c['class_min']}\n" for c in ch)
    for p, t in ((TASKS, txt), (SHAF, sha + "\n"), (CHUNKF, ctxt)):
        if os.path.exists(p) and open(p).read() != t: sys.exit(f"STOP: {p} exists with other content -- not overwriting")
    base_hits = check_base() if "--check-base" in sys.argv else None
    used, nlists, kinds, ncmd = used_seeds()
    seeds = [int(l.split()[4]) for _, l in lines]
    if len(set(seeds)) != len(seeds) or set(seeds) & used: sys.exit("STOP: the seeds are not all distinct or overlap an earlier list or command")
    if "--dry-run" not in sys.argv:
        for p, t in ((TASKS, txt), (SHAF, sha + "\n"), (CHUNKF, ctxt)):
            if not os.path.exists(p): open(p, "w").write(t)
    print(f"{os.path.relpath(TASKS, HS)}: {len(lines)} lines, SHA-256 {sha}; {os.path.relpath(CHUNKF, HS)}: {len(ch)} chunks\n")
    if "--design" in sys.argv: design(rows, lines, ch, used, nlists, kinds, ncmd, base_hits)


if __name__ == "__main__":
    main()
