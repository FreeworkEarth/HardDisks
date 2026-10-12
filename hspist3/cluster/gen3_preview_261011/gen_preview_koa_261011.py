#!/usr/bin/env python3
"""##CHRIS 2026-10-11 (plan-author decision 29, task P; 261012 sec. 4.7.59): REGISTRATION of THE PREVIEW SWEEP ON KOA -- PREVIEW,
INFORMATION; no verdict, never part of stage 1, never in a paper; its own root on scratch (preview_koa) and its own seeds. This script,
the task list it writes (tasks_preview_koa_261011.txt), that list's SHA-256 (tasks_preview_koa_261011.sha256) and the chunk size
(preview_chunk_size.txt) are committed together before any preview trajectory exists, and none of them is edited afterwards. The
binary is Test R's: build_record_a643bc7.txt (cluster/gen3_testR_261011/), read unchanged.
DESIGN (decision 29):
  command     the canonical A1 v2 command (methods sec. 8; tests_20260913.argv_for with jobs_D's arguments) through the B command of
              cluster/gen3_koa_261009/testG_worker_koa.sh (unchanged since 1cdb7ff6), --engine=gen3. Its differences from the A1 v2
              command besides --engine are printed by --design.
  cells       N = 100 (N_s = 50), H = 10, radius 0.5, wall thickness 0.05: CELL_SET below (the order's densities; see CELL_SETS).
  masses      the nine A1 v2 masses tests_20260913.A1_MASSES.
  seeds       25 per (cell, mass), r = 0..24, by jobs_D's rule with a new base: per cell i the base BASE_P + i, and the exact seed
              run_seed(BASE_P + i, 0, mass index in A1_MASSES, r). BASE_P is a number that occurs in no file (--check-base).
  stride      jobs_D's rule: d_stride(nu_pred(nu50, M)), with nu50 the binary's own predicted frequency at M = 50. It is computed by
              f_pred() below, a replica of 00ALLINONE.c's speed_sound_henderson_frequency_estimate (Henderson a = 1/8, the gamma = 2
              mapping, L_eff = L0 - 2r, L0 as float). --design checks the replica against the Predicted_Frequency of every A1 v2 leaf
              and Test G's nu_pred, and checks the rule against the stride of every A1 v2 trace.
  line format Test G's: "B preview/gen3/<cell>/m_<M> <M> <r> <seed> <L0> 10.0 50 <stride> <base> gen3 main".
  ORDER       interleaved, r by r: for r = 0..24, for each cell, for each mass. One r-round is 7 x 9 = 63 lines.
  CHUNKS      K = 63 (preview_chunk_size.txt): every chunk is one r-round, i.e. every (cell, mass) once, so the chunk times are
              even. 1575 lines = 25 chunks, the array 1-25.
--design prints the cells, the strides, the seed checks, the commands and their differences, and the predicted times. The prediction
is A1 v2's measured seconds per trajectory on the Mac (gen2; its run.log headers), times the ratio of Test G's gen3 seconds on KOA to
A1 v2's at pi/8 (M = 50 and 300). Cells that are not A1 v2 leaves are interpolated in log L0.
Writes (only when CELL_SET is set and the file is absent; an existing file with other content stops the script):
tasks_preview_koa_261011.txt, tasks_preview_koa_261011.sha256 (one line, the hex digest), preview_chunk_size.txt (one line, K).
usage (from hspist3/): python3 cluster/gen3_preview_261011/gen_preview_koa_261011.py [--design] [--check-base] [--dry-run] [--cells named|a1v2]
  --cells   design printout for one candidate set (only while CELL_SET is unset); --dry-run   print everything, write nothing
"""
import glob, hashlib, math, os, re, subprocess, sys
from collections import Counter
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); CL = os.path.dirname(HERE); HS = os.path.dirname(CL); REPO = os.path.dirname(HS)
GEN3_WT = os.path.join(os.path.dirname(REPO), "HardDisks_gen3")
sys.path.insert(0, os.path.join(HS, "validation")); sys.path.insert(0, HS)
import tests_20260913 as T
BASE_P, K, NSEED, REL, ENGINE = None, 63, 25, "preview", "gen3"          # BASE_P: chosen by --check-base before the registration
TASKS = os.path.join(HERE, "tasks_preview_koa_261011.txt")
SHAF = os.path.join(HERE, "tasks_preview_koa_261011.sha256")
KF = os.path.join(HERE, "preview_chunk_size.txt")
TESTG_KOA = os.path.join(HS, "experiments_gen3_gate_koa_261009", "testG")
N, NS, R_DISK, H_STR, H = 100, 50, 0.5, "10.0", 10.0
ETA_NAMED = (("0p05", 0.05), ("0p10", 0.10), ("0p20", 0.20), ("0p30", 0.30), ("pi8", math.pi / 8), ("0p50", 0.50), ("0p60", 0.60))
# the two readings of "the canonical A1 v2 command and cells for N = 100 at eta = 0.05, 0.10, 0.20, 0.30, pi/8, 0.50, 0.60":
#   named  the seven densities exactly: L0 = N pi r^2 / (2 H eta), to 6 decimals (A1 v2 has a cell at pi/8 only);
#   a1v2   the A1 v2 leaf nearest each named density, with the leaf's L0 string ("0p60": 0.590001 and 0.609999 are equally near; the
#          first is listed here, and the plan author chooses).
CELL_SETS = {
    "named": [(f"eta{t}_H10_L{N * math.pi * R_DISK ** 2 / (2 * H * e):.6f}", f"{N * math.pi * R_DISK ** 2 / (2 * H * e):.6f}", e, t)
              for t, e in ETA_NAMED],
    "a1v2": [("eta_0p052360", "74.9998", 0.05, "0p05"), ("eta_0p112200", "34.9999", 0.10, "0p10"), ("eta_0p196350", "20.0000", 0.20, "0p20"),
             ("eta_0p261799", "15.0000", 0.30, "0p30"), ("eta_0p392699", "10.0000", math.pi / 8, "pi8"), ("eta_0p523599", "7.5000", 0.50, "0p50"),
             ("eta_0p590", "6.6559", 0.60, "0p60")],
}
CELL_SET = None                                                         # set in the registration commit


def k_bin(alpha):
    """00ALLINONE.c speed_sound_fundamental_k: 120 bisections of cot K = alpha K"""
    lo, hi = 1e-10, 0.5 * math.pi - 1e-10
    for _ in range(120):
        mid = 0.5 * (lo + hi)
        if 1.0 / math.tan(mid) - alpha * mid > 0.0: lo = mid
        else: hi = mid
    return 0.5 * (lo + hi)


def f_pred(L0s, M):
    """the binary's predicted frequency (speed_sound_henderson_frequency_estimate at eta_nominal, L0_UNITS as float)"""
    L0 = float(np.float32(float(L0s))); r = 12.0 / 24.0; Hs = 240.0 / 24.0
    eta = (N * math.pi * r * r) / max(1e-12, 2.0 * L0 * Hs)
    a = 0.125; num = 1.0 + eta + 3.0 * a * eta * eta - a * eta * eta * eta
    cs = math.sqrt(max(0.0, 2.0 * 1.0 / 1.0 * num / ((1 - eta) * (1 - eta) * (1 - eta))))
    return cs * k_bin(M / (2.0 * NS)) / (2.0 * math.pi * (L0 - 2 * r))


def stride(L0s, M):
    """jobs_D's rule: d_stride(nu_pred(nu50, M)), nu50 the binary's f_pred at M = 50"""
    return T.d_stride(T.nu_pred(f_pred(L0s, 50), M))


def cells():
    if CELL_SET is None: sys.exit("STOP: CELL_SET is not set (the registration fixes it)")
    return CELL_SETS[CELL_SET]


def lines(cs):
    mi = T.A1_MASSES.index; L = []
    for r in range(NSEED):
        for i, (cid, L0s, eta, tag) in enumerate(cs):
            for M in T.A1_MASSES:
                L.append(f"B {REL}/{ENGINE}/{cid}/m_{M} {M} {r} {T.run_seed(BASE_P + i, 0, mi(M), r)} {L0s} {H_STR} {NS} {stride(L0s, M)} "
                         f"{BASE_P + i} {ENGINE} main\n")
    return L


def used_seeds():
    """every exact seed of every earlier task list (both trees; this folder excluded), the extra seed sets of the Test G generators,
    the symmetry test's seeds (as Test R's generator), and the seeds of tests_20260913's tests A, B and D (D is A1 v2)"""
    used, nlists, kinds = set(), 0, Counter()
    for t in (HS, os.path.join(GEN3_WT, "hspist3")):
        for f in glob.glob(os.path.join(t, "cluster", "*", "tasks_*.txt")):
            if os.path.basename(os.path.dirname(os.path.abspath(f))) == os.path.basename(HERE): continue
            nlists += 1
            for l in open(f):
                p = l.split()
                if not p: continue
                kinds[p[0]] += 1
                if p[0] == "B": used.add(int(p[4]))
                elif p[0] in ("A", "AF"): used.add(int(p[3]))
    used |= {T.run_seed(20261013, 0, m, r) for m in range(9) for r in range(4)}
    used |= {T.run_seed(20261091, 5, 0, r) for r in range(8)} | {9700, 9701, 9702, 9703, 9711}
    for t in (HS, os.path.join(GEN3_WT, "hspist3")):
        sym = os.path.join(t, "experiments_gen3_sym_261009", "sym_seeds.txt")
        if os.path.exists(sym): used |= {int(l.split()[3]) for l in open(sym) if l.strip()}
    abd = {T.run_seed(T.A_BASE + li, 0, mi, r) for li in range(len(T.L0S)) for mi in range(len(T.ROMAN_MASSES)) for r in range(T.SEEDS)}
    abd |= {T.run_seed(T.B_BASE + li, 0, mi, r) for li in range(len(T.L0S)) for mi in range(len(T.B_MASSES)) for r in range(T.SEEDS)}
    abd |= {T.run_seed(T.D_BASE + ei, 0, mi, r) for ei in range(40) for mi in range(len(T.A1_MASSES)) for r in range(T.SEEDS)}
    return used | abd, nlists, kinds, len(abd)


def check_base(base):
    """the base and base + 1 .. base + 6 as integer tokens in any file: both working trees (rg), the tracked files of three heads (git
    grep); this package folder excluded"""
    hits, ghits = [], []
    for b in range(base, base + 7):
        pat = rf"(^|[^0-9.]){b}([^0-9.]|$)"
        rg = subprocess.run(["rg", "-l", "--no-ignore", "--hidden", "-g", "!.git", "-g", f"!{os.path.relpath(HERE, REPO)}/**", "-e", pat,
                             REPO, GEN3_WT], capture_output=True, text=True)
        hits += [h for h in rg.stdout.splitlines() if "/cluster/gen3_preview_261011/" not in h]
        gg = subprocess.run(["git", "-C", REPO, "grep", "-E", "-l", pat, "main", "engine-gen3", "engine-divider-resched", "--"],
                            capture_output=True, text=True)
        ghits += [h for h in gg.stdout.splitlines() if "/cluster/gen3_preview_261011/" not in h]
    return hits, ghits


def a1v2_checks():
    """the f_pred replica against every A1 v2 leaf's Predicted_Frequency and Test G's nu_pred; jobs_D's stride rule against the stride
    of every A1 v2 trace (the time step between its first two rows / (0.4f / 24))"""
    import csv
    leaves = T.a1_leaf_table(); okf = oks = n = 0; bad = []
    dt = float(np.float32(0.4)) / 24.0
    for l in leaves:
        f50 = sorted(glob.glob(f"{T.A1_DIR}/{l['leaf']}/wall_x_positions_*wallmassfactor_50_run0.csv")
                     + glob.glob(f"{T.LOWETA_DIR}/{l['leaf']}/wall_x_positions_*wallmassfactor_50_run0.csv"))[0]
        okf += f"{f_pred(l['L0'], 50):.12g}" == next(csv.DictReader(open(f50)))["Predicted_Frequency"]
        for M in T.A1_MASSES:
            tr = sorted(glob.glob(os.path.join(T.DROOT, l["leaf"], f"m_{M}", f"wall_x_positions_*_run0.csv")))
            if not tr: continue
            with open(tr[0]) as fh:
                rd = csv.DictReader(fh); a = next(rd); b = next(rd)
            st = round((float(b["Time"]) - float(a["Time"])) / dt); n += 1
            if st == stride(l["L0"], M): oks += 1
            else: bad.append((l["leaf"], M, st, stride(l["L0"], M)))
    import pandas as pd
    tg = [(M, float(pd.read_csv(os.path.join(TESTG_KOA, "gen3", "epi8_H_H10_L10", f"m_{M}", "red_nu.csv"))["nu_pred"].iloc[0]))
          for M in (50, 300, 1500, 2000)]
    okg = sum(f"{f_pred('10.000000', M):.12g}" == f"{v:.12g}" for M, v in tg)
    return okf, len(leaves), oks, n, bad, okg, len(tg)


def seconds_model():
    """per (L0, M): predicted E5 seconds of one gen3 trajectory = A1 v2's mean Mac seconds (gen2, run.log headers) x the ratio of Test
    G's gen3 seconds on KOA to A1 v2's at pi/8 (M = 50 and 300, averaged); cells between A1 v2 leaves interpolated in log L0"""
    leaves = T.a1_leaf_table(); mac = {}
    for l in leaves:
        for M in T.A1_MASSES:
            p = os.path.join(T.DROOT, l["leaf"], f"m_{M}", "run.log")
            s = [int(x) for x in re.findall(r"^##RUN .*\((\d+) s\)", open(p, errors="replace").read(), re.M)] if os.path.exists(p) else []
            if s: mac[(float(l["L0"]), M)] = sum(s) / len(s)
    koa = {}
    for M in (50, 300):
        txt = open(os.path.join(TESTG_KOA, "gen3", "epi8_H_H10_L10", f"m_{M}", "run.log"), errors="replace").read()
        s = [int(x) for x in re.findall(r"^##RUN .* node \S+ \((\d+) s\)", txt, re.M)]; koa[M] = sum(s) / len(s)
    ratio = 0.5 * sum(koa[M] / mac[(10.0, M)] for M in (50, 300))
    Ls = sorted({L for L, M in mac})

    def sec(L0s, M):
        L0 = float(L0s)
        if (L0, M) in mac: return mac[(L0, M)] * ratio
        lo = max([L for L in Ls if L <= L0] or [Ls[0]]); hi = min([L for L in Ls if L >= L0] or [Ls[-1]])
        if lo == hi: return mac[(lo, M)] * ratio
        w = (math.log(L0) - math.log(lo)) / (math.log(hi) - math.log(lo))
        return (mac[(lo, M)] ** (1 - w) * mac[(hi, M)] ** w) * ratio
    return sec, ratio, koa, {M: mac[(10.0, M)] for M in (50, 300)}


def design(cs, L, used, nlists, kinds, nabd, base_hits):
    print(f"# THE PREVIEW SWEEP ON KOA -- the design (decision 29, task P; PREVIEW, INFORMATION): cell set '{CELL_SET or '(candidate)'}'\n")
    okf, nl, oks, ns, bad, okg, ng = a1v2_checks()
    print(f"the binary's predicted frequency, replicated (f_pred): equal to the Predicted_Frequency of {okf} of {nl} A1 v2 leaves "
          f"(279282b) and to Test G's nu_pred at {okg} of {ng} masses (a643bc7), as printed (12 significant digits)")
    print(f"jobs_D's stride rule d_stride(nu_pred(nu50, M)) equals the stride of the A1 v2 trace in {oks} of {ns} cells"
          + ("" if not bad else f"; differing: {bad[:6]}"))
    print("\n| cell | eta named | L0 | eta = N pi r^2/(2 L0 H) | box_delta(L0) [sigma] | nu50 (f_pred) | strides M = " +
          ", ".join(map(str, T.A1_MASSES)) + " | base |\n|---|---|---|---|---|---|---|---|")
    import paper1_populate_cs_err_20261002 as E0
    for i, (cid, L0s, eta, tag) in enumerate(cs):
        L0 = float(L0s)
        print(f"| {cid} | {tag} | {L0s} | {N * math.pi * R_DISK ** 2 / (2 * L0 * H):.6f} | {E0.box_delta(L0):.6f} | {f_pred(L0s, 50):.6g} | "
              + " ".join(str(stride(L0s, M)) for M in T.A1_MASSES) + f" | {(BASE_P + i) if BASE_P else '-'} |")
    if L:
        seeds = [int(l.split()[4]) for l in L]
        print(f"\n{len(L)} task lines ({len(cs)} cells x {len(T.A1_MASSES)} masses x {NSEED} seeds), order r by r; distinct exact seeds "
              f"{len(set(seeds))} of {len(seeds)}; overlap with {nlists} earlier task lists (both trees; line types {dict(kinds)}), "
              f"the Test G extra sets, the symmetry seeds and tests A, B, D ({nabd} seeds; D = A1 v2): {len(set(seeds) & used)}")
    if base_hits is not None:
        h, g = base_hits
        print(f"the bases {BASE_P}..{BASE_P + 6} as integer tokens: {len(h)} files in both working trees, {len(g)} tracked files in main, "
              f"engine-gen3, engine-divider-resched (this folder excluded)" + ("" if (h or g) else " -- they occur in no file"))
    sec, ratio, koa, mac = seconds_model()
    print(f"\n## chunks: K = {K} lines per array task = one r-round (every cell and mass once); {-(-NSEED * len(cs) * 9 // K)} chunks\n")
    print(f"the time model: A1 v2's Mac seconds (gen2) x {ratio:.3f}, the ratio of Test G's gen3 seconds on KOA (M = 50: {koa[50]:.2f} s, "
          f"M = 300: {koa[300]:.2f} s; nodes as they came) to A1 v2's at pi/8 (M = 50: {mac[50]:.1f} s, M = 300: {mac[300]:.1f} s)\n")
    print("| cell | " + " | ".join(f"M = {M} [s]" for M in T.A1_MASSES) + " | one r-round [s] |\n|---|" + "---|" * (len(T.A1_MASSES) + 1))
    tot = 0.0
    for cid, L0s, eta, tag in cs:
        v = [sec(L0s, M) for M in T.A1_MASSES]; tot += sum(v)
        print(f"| {cid} | " + " | ".join(f"{x:.1f}" for x in v) + f" | {sum(v):.0f} |")
    print(f"\npredicted chunk time (one r-round): {tot / 60:.1f} min; the whole sweep {tot * NSEED / 3600:.1f} core-hours; requested "
          f"{-(-NSEED * len(cs) * 9 // K)} chunks x 1 h (preview_koa.sbatch --time=1:00:00)")


def main():
    global CELL_SET
    if "--cells" in sys.argv:
        if CELL_SET is not None: sys.exit("STOP: --cells is for the candidates only; CELL_SET is set")
        CELL_SET = sys.argv[sys.argv.index("--cells") + 1]; cs = CELL_SETS[CELL_SET]; dry = True
    else:
        cs = cells(); dry = "--dry-run" in sys.argv
    if BASE_P is None and not dry: sys.exit("STOP: BASE_P is not set")
    L = lines(cs) if BASE_P is not None else []
    base_hits = check_base(BASE_P) if ("--check-base" in sys.argv and BASE_P is not None) else None
    used, nlists, kinds, nabd = used_seeds()
    if L:
        seeds = [int(l.split()[4]) for l in L]
        if len(set(seeds)) != len(seeds) or set(seeds) & used: sys.exit("STOP: the seeds are not all distinct or overlap an earlier list")
        txt = "".join(L); sha = hashlib.sha256(txt.encode()).hexdigest()
        for p, t in ((TASKS, txt), (SHAF, sha + "\n"), (KF, f"{K}\n")):
            if os.path.exists(p) and open(p).read() != t: sys.exit(f"STOP: {p} exists with other content -- not overwriting")
        if base_hits and (base_hits[0] or base_hits[1]) and not os.path.exists(TASKS):
            sys.exit(f"STOP: a base occurs in {base_hits[0][:5]} {base_hits[1][:5]}")
        if not dry:
            for p, t in ((TASKS, txt), (SHAF, sha + "\n"), (KF, f"{K}\n")):
                if not os.path.exists(p): open(p, "w").write(t)
        print(f"{os.path.relpath(TASKS, HS)}: {len(L)} lines, SHA-256 {sha}; chunk size {K} -> {-(-len(L) // K)} chunks"
              + (" (dry run: nothing written)" if dry else "") + "\n")
    if "--design" in sys.argv: design(cs, L, used, nlists, kinds, nabd, base_hits)


if __name__ == "__main__":
    main()
