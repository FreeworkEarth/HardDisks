#!/usr/bin/env python3
"""##CHRIS 2026-10-15 (plan-author decision 57; notes 261012 sec. 4.7.135): THE TASK GENERATOR OF THE STAGE 2 PILOT "P2" ON KOA (not
registered; the registration follows its result). It writes the task list tasks_p2_koa_261016.txt (400 lines), its SHA-256
(tasks_p2_koa_261016.sha256) and the chunk table chunks_p2_koa_261016.tsv (one line per chunk), and prints the design.
CELLS (decision 57): N = 900; eta in {0.72, 0.78}; two geometries per cell, two pins per geometry:
  G1  the stage 1 rule: H = 10 sqrt(N/100) = 30, L0 = N pi / (8 H eta);
  G2  strain free: the same pinned arrangement as G1 (rows, n_s x n_r, vacancies: the binary's own choice, p2_lattice.lattice_fit,
      checked against the binary in 42 of 42 starts, Report D56-1 (b)) at the height H that makes its strain e closest to 0, N and eta
      kept, L0 = N pi / (8 H eta) following. H lies on the binary's grid: the box height is (int)(H * 24) px in float arithmetic
      (00ALLINONE.c:326), so H = k/24 with k a whole number, and every H written here is checked through the C library's strtof
      (the binary's parser) and the float product. Searched: k = 600 .. 840 (H 25 .. 35); the smallest |e|, then the k nearest 720.
  pins: --gen3-lattice-orientation 0 and turned (the stage 1 worker's start field: 0 | turned), in both compartments, the divider at L0.
MASSES: alpha = M/N in {1, 2, 5, 10, 20}.   SEEDS: 10 per (cell, G, pin, mass), s = 0..9, by the stage 1 rule (gen_stage1_koa_261012:
  run_seed(BASE, stream, mass index, s), one stream per cell in the order of the task list, the mass index the alpha's index in the
  nine-mass ladder) with BASE = 20261016 (decision 57); every seed checked against the earlier task lists and recorded commands.
PROTOCOL (stage 1, unchanged): T_eq = 5000 sigma-time (300000 held steps); 200 predicted periods of the binary's own nu_pred (Henderson,
  replicated exactly by gen_stage1_koa_261012.nu_pred_binary, which uses the box height k/24 as the binary does); Delta_psi =
  max(1, T_rec / 2000); virial blocks of 50; the lattice start. Snapshots on EVERY run every (T_eq + T_rec) / 100 (stage 1: on its
  snapshot subset only), and the initial state of every run (--gen3-dump-initial; p2_worker_koa.sh). The trace stride: d_stride of the
  binary's nu_pred (the stage 1 rule takes nu from KR, which is not defined beyond eta 0.7069; as in the D56-1 pilot).
TIME (Report D56-1 (a)): events per disk and sigma-time r = 1.2004 (Z - 1) + 0.0112 with the upper Z (Liu's solid branch x 1.268, the
  pilot's measured excess), times the Mac rehearsal's measured ratio per pin (decision 57 item 2; experiments_gen3_p2_261016: the four
  eta 0.78, alpha 5 lines executed 1.171 (G1 pin 0), 1.158 (G2 pin 0), 1.002 (G1 turned) and 0.929 (G2 turned) times the table's
  events; the pin 0 crystals, 12 vacancies per compartment, carry Z 20.7 against the table's 18.1): 1.171 for pin 0, 1.002 for turned,
  at both etas; microseconds per event, the dense cells' medians at N 900: 3.533 (E5-2670 v3), 4.242 (E5-4640 v2); the array's time
  limit = 1.5 x the longest line on the slower model, rounded up to 5 min (decision 57: the table + 50 %).
TASK LINE (15 fields, s1_worker_koa.sh's): B <cell folder> <M> <s> <seed> <L0> <H> <N_s> <stride> <BASE> <hold> <Delta_psi> <T_snap>
  <pin: 0 | turned> P2; the cell folder stage2/gen3/<cell id>/m_<M>, the cell id P2_N0900_eta<6 decimals>_<G1|G2>_pin<0|turned>.
ORDER: by eta, then G, then pin, then alpha, then s. CHUNKS: one line each (chunk j = line j): one array of 400 one-core tasks.
usage (from hspist3/): python3 cluster/gen3_stage2_261016/gen_p2_koa_261016.py [--dry-run] [--out <a trial folder>]
"""
import ctypes, ctypes.util, hashlib, math, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); CL = os.path.dirname(HERE); HS = os.path.dirname(CL)
for p in (HERE, os.path.join(CL, "gen3_stage1_261012"), os.path.join(HS, "validation"), HS): sys.path.insert(0, p)
import p2_lattice as LT
import gen_stage1_koa_261012 as GEN
import tests_20260913 as T
import plot_speed_of_sound_edmd as SOS

BASE = 20261016
N = 900; NS = N // 2; ETAS = (0.72, 0.78); GEOMS = ("G1", "G2"); PINS = ("0", "turned"); ALPHAS = (1, 2, 5, 10, 20); NSEED = 10
TEQ = 5000.0; HOLD = int(round(60 * TEQ)); KRANGE = range(600, 841)
RATE_A, RATE_B, ZX = 1.2004, 0.0112, 1.268; US_V3, US_V2 = 3.533, 4.242; MEAS = {"0": 1.171, "turned": 1.002}
OUTD = os.path.abspath(sys.argv[sys.argv.index("--out") + 1]) if "--out" in sys.argv else HERE      # --out: a trial folder instead
TASKS = os.path.join(OUTD, "tasks_p2_koa_261016.txt"); SHAF = os.path.join(OUTD, "tasks_p2_koa_261016.sha256")
CHUNKF = os.path.join(OUTD, "chunks_p2_koa_261016.tsv")
_libc = ctypes.CDLL(ctypes.util.find_library("c")); _libc.strtof.restype = ctypes.c_float
_libc.strtof.argtypes = [ctypes.c_char_p, ctypes.POINTER(ctypes.c_char_p)]


def height_px(text):
    """the binary's box height in px for the --height text: (int)(strtof(text) * 24.0f), in float arithmetic"""
    v = np.float32(_libc.strtof(text.encode(), None)); return int(np.float32(v * np.float32(24.0)))


def same(a, b): return all(a[k] == b[k] for k in ("orient", "alt", "ns", "nr"))


def geometry(eta, G, pin):
    req = 0 if pin == "0" else 1; R = 0.5 * LT.PX
    def arr(k):
        H = k / 24.0; L0 = N * math.pi / (8.0 * H * eta)
        return H, L0, LT.lattice_fit(NS, R, L0 * LT.PX - 0.5 * LT.TH * LT.PX, float(k), req)
    H1, L01, A1 = arr(720)
    if G == "G1": k, H, L0, A = 720, H1, L01, A1
    else:
        cand = []
        for kk in KRANGE:
            H, L0, b = arr(kk)
            if b is not None and same(b, A1): cand.append((abs(LT.strain(b)), abs(kk - 720), kk, H, L0, b))
        if not cand: raise SystemExit(f"STOP: no height in k = {KRANGE.start}..{KRANGE.stop - 1} keeps G1's arrangement at eta {eta}, pin {pin}")
        _, _, k, H, L0, A = min(cand)
    htxt = "30" if k == 720 else f"{H:.17g}"
    if height_px(htxt) != k: raise SystemExit(f"STOP: --height={htxt} gives {height_px(htxt)} px, not {k}")
    ws = LT.w_star(A, R, L0 * LT.PX - 0.5 * LT.TH * LT.PX, float(k)) / LT.PX
    return dict(eta=eta, G=G, pin=pin, k=k, H=H, htxt=htxt, L0=L0, A=A, e=LT.strain(A), W=L0 - LT.TH / 2, Wstar=ws,
                eta_star=N * math.pi / 4.0 / ((2 * ws + LT.TH) * H), aspect=(L0 - LT.TH / 2) / H)


def plan():
    cells, rows, lines = [], [], []
    for eta in ETAS:
        for G in GEOMS:
            for pin in PINS:
                cells.append(geometry(eta, G, pin))
    for stream, c in enumerate(cells):
        cid = f"P2_N{N:04d}_eta{c['eta']:.6f}_{c['G']}_pin{c['pin']}"; c["cid"] = cid; c["stream"] = stream
        Zu = float(SOS.Z_liu_solid(np.array([c["eta"]]))[0]) * ZX; r = RATE_A * (Zu - 1.0) + RATE_B
        for a in ALPHAS:
            M = int(round(a * N)); mi = GEN.ALPHA.index(a)
            nup = GEN.nu_pred_binary(N, c["H"], c["L0"], M); steps = math.ceil((200.0 / nup) / GEN.DTF); Trec = steps * GEN.DTF
            stride = T.d_stride(nup); dpsi = max(1.0, Trec / 2000.0); tsnap = (TEQ + Trec) / 100.0
            ev = r * N * (TEQ + Trec) * MEAS[c["pin"]]; rel = f"stage2/gen3/{cid}/m_{M}"
            seeds = [T.run_seed(BASE, stream, mi, s) for s in range(NSEED)]
            for s, sd in enumerate(seeds):
                lines.append(f"B {rel} {M} {s} {sd} {c['L0']:.17g} {c['htxt']} {NS} {stride} {BASE} {HOLD} {dpsi:.17g} {tsnap:.17g} {c['pin']} P2\n")
            rows.append(dict(cid=cid, eta=c["eta"], G=c["G"], pin=c["pin"], a=a, M=M, mi=mi, stream=stream, nup=nup, steps=steps, Trec=Trec,
                             stride=stride, dpsi=dpsi, tsnap=tsnap, Zu=Zu, r=r, ev=ev, min_v3=ev * US_V3 / 6e7, min_v2=ev * US_V2 / 6e7, seeds=seeds))
    return cells, rows, lines


def collisions(seeds):
    used, nlists, kinds, cmds = GEN.used_seeds()
    st1 = [int(l.split()[4]) for l in open(GEN.TASKS) if l.startswith("B ")]
    extra = set(st1) | {T.run_seed(20261015, 0, 4, s) for s in range(2)} | {T.run_seed(20261017, ci, si, 0) for ci in range(14) for si in range(3)}
    hit = sorted(set(seeds) & (used | extra))
    return hit, len(used | extra), nlists + 1, cmds


def main():
    cells, rows, lines = plan()
    L = ["# THE STAGE 2 PILOT P2 (decision 57): the design", "",
         "## The cells (N 900): G1 the stage 1 geometry; G2 the same pinned arrangement at the height that makes it strain free", "",
         "| cell | H [sigma] (k px) | L0 | aspect W/H | the pinned arrangement (left = right) | vacancies | gap [sigma] | start strain e | "
         "W* at this H (eta*) |", "|---|---|---|---|---|---|---|---|---|"]
    for c in cells:
        A = c["A"]
        L.append(f"| {c['cid']} | {c['htxt']} ({c['k']}) | {c['L0']:.6f} | {c['aspect']:.4f} | {LT.fmt_arr(A)} | {A['nvac']} | {A['s'] / LT.PX:.4f} | "
                 f"{100 * c['e']:+.3f} % | {c['Wstar']:.4f} ({c['eta_star']:.4f}) |")
    L += ["", "## Per (cell, mass): the binary's nu_pred, the record, the stride, the time (upper Z, Report D56-1 (a), x the rehearsal's ratio per pin)", "",
          "| cell | alpha | M | nu_pred | T_rec | stride | Delta_psi | T_snap | Z upper | events [M] | min (v3) | min (v2) |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        L.append(f"| {r['cid']} | {r['a']} | {r['M']} | {r['nup']:.6f} | {r['Trec']:.1f} | {r['stride']} | {r['dpsi']:.4f} | {r['tsnap']:.3f} | "
                 f"{r['Zu']:.3f} | {r['ev'] / 1e6:.1f} | {r['min_v3']:.2f} | {r['min_v2']:.2f} |")
    longest = max(r["min_v2"] for r in rows); limit = int(math.ceil(1.5 * longest / 5.0) * 5)
    ch_v3 = sum(r["min_v3"] for r in rows) * NSEED / 60.0; ch_v2 = sum(r["min_v2"] for r in rows) * NSEED / 60.0
    seeds = [s for r in rows for s in r["seeds"]]; hit, nused, nl, ncmd = collisions(seeds)
    L += ["", f"- task lines: {len(lines)} (8 cells x 5 masses x 10 seeds); distinct seeds: {len(set(seeds))}; seeds already used by an "
              f"earlier list, the stage 1 list, the D56-1 pilot or its lattice starts, or a recorded command ({nused} seeds in {nl} lists, "
              f"{ncmd} command files): {len(hit)}{' ' + str(hit[:10]) if hit else ''}",
          f"- predicted core-hours (upper Z x the rehearsal's ratio per pin): {ch_v3:.1f} on E5-2670 v3, {ch_v2:.1f} on E5-4640 v2; the longest line {longest:.1f} min on v2 "
          f"-> the array's time limit 1.5 x that, rounded up to 5 min: {limit} min",
          f"- wall time at 64 cores (v3): about {ch_v3 * 60 / 64:.0f} min"]
    if len(set(seeds)) != len(seeds) or hit: raise SystemExit("\n".join(L) + "\nSTOP: a seed repeats or was used before")
    txt = "".join(lines); sha = hashlib.sha256(txt.encode()).hexdigest()
    ctxt = "chunk\tfirst\tlast\tlines\tgroup\tcell\tM\tpredicted_s_v2\tlimit_min\tclass_min\n"
    by = {(r["cid"], r["M"]): r for r in rows}
    for j, l in enumerate(lines, 1):
        f = l.split(); cid = f[1].split("/")[2]; r = by[(cid, int(f[2]))]
        ctxt += f"{j}\t{j}\t{j}\t1\tP2\t{cid}\t{r['M']}\t{60 * r['min_v2']:.1f}\t{1.5 * r['min_v2']:.1f}\t{limit}\n"
    L += [f"- {os.path.relpath(TASKS, HS)}: SHA-256 {sha}; {os.path.relpath(CHUNKF, HS)}: {len(lines)} chunks"]
    print("\n".join(L))
    if "--dry-run" in sys.argv: print("\n--dry-run: nothing written"); return
    for p in (TASKS, SHAF, CHUNKF):
        if os.path.exists(p): raise SystemExit(f"REFUSED: {p} exists; the package is written once")
    open(TASKS, "w").write(txt); open(SHAF, "w").write(sha + "\n"); open(CHUNKF, "w").write(ctxt)
    open(os.path.join(OUTD, "p2_design_261016.txt"), "w").write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
