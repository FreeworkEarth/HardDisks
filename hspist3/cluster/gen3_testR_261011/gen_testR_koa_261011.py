#!/usr/bin/env python3
"""##CHRIS 2026-10-11 (plan-author decision 28; 261012 sec. 4.7.57): REGISTRATION of TEST R ON KOA, the replication of Test G at the
light masses. This script, the task list it writes (tasks_testR_koa_261011.txt), that list's SHA-256 (tasks_testR_koa_261011.sha256),
the chunk size (testR_chunk_size.txt) and the build record of the KOA binary (build_record_a643bc7.txt) are committed together, before
any Test R trajectory exists; none of them is edited afterwards.
DESIGN (decision 28):
  cell        epi8_H_H10_L10 (N = 100, N_s = 50, H = L0 = 10, pi/8): Test G's main cell; the trace stride by the Test G generator's
              rule (cluster/gen3_gate_261009/gen_testG_261009.campaign_geo: the campaign's stride per mass, 26 / 32 / 51 at M = 50 / 100 /
              300);
  masses      M = 50, 100 and 300;
  seeds       3,600 exact seeds per mass, r = 0..3599: run_seed(BASE_R, 0, mass index in tests_20260913.A1_MASSES, r) (the driver's
              speed_sound_run_seed), the same exact seed in both engines; BASE_R = 20261317, a number that occurred in no file of
              the repository before this registration (--check-base: the integer as a token, not inside a longer number, in every
              file of both working trees and in the tracked files of main, engine-gen3 and engine-divider-resched; the Test R
              package folder itself excluded);
  engines     gen2 and gen3, the B command of testG_worker_koa.sh (unchanged since 1cdb7ff6), --engine=<gen2|gen3>;
  line format Test G's: "B testR/<engine>/epi8_H_H10_L10/m_<M> <M> <r> <seed> 10.000000 10.000000 50 <stride> <BASE_R> <engine> main";
  ORDER       interleaved, r by r: for r = 0..3599, for M = 50, 100, 300, gen2 then gen3 (so every chunk holds the same mix of the six
              arms and the chunk times are even; the merged run.log of every cell is in r order);
  CHUNKS      K = 100 task lines per array task (testR_chunk_size.txt, read by testR_koa.sbatch and the merge job): 21,600 lines = 216
              chunks, the array 1-216, far below KOA's MaxArraySize 25001 (KOA_RUNSHEET_261002.md sec. 14); each chunk writes only its
              own folder chunks/chunk_<j>.
--design prints the tables of the design: lines per mass and engine, the chunks, the predicted chunk time on E5 (Test G's measured
seconds per trajectory on KOA, from the fetched run.log headers: E5-2670 v3 nodes only, worker included; M = 100 interpolated in the
simulated time), the core-hours, the expected data on KOA scratch, the seed checks.
Writes (only if absent; an existing file with other content stops the script): tasks_testR_koa_261011.txt,
tasks_testR_koa_261011.sha256 (one line, the hex digest), testR_chunk_size.txt (one line, K), build_record_a643bc7.txt (below).
THE BUILD RECORD (decision 28: "commit hash and binary SHA-256 read from files written at build time; write them now from the
recorded values"): the KOA build of a643bc7 (runsheet step 2, build_gen3_koa.sh) wrote logs/BUILD_KOA_LAST.hash (the binary's SHA-256)
on KOA and printed it; no file there holds the commit. build_record_a643bc7.txt is written now, on the Mac, from two recorded values:
the 00ALLINONE line of cluster/gen3_koa_261009/k2_binaries_a643bc7.sha256 (the SHA-256 step 2 printed on KOA, checked again on KOA by
K2-1's sha256sum -c: OK) and the full hash of engine-gen3 a643bc7 (the commit the binary's own version line names). testR_koa.sbatch
reads it and checks, without git: the binary's SHA-256 equals the record AND the build-time file logs/BUILD_KOA_LAST.hash on KOA, and
the binary's version line is the recorded one.
usage (from hspist3/): python3 cluster/gen3_testR_261011/gen_testR_koa_261011.py [--design] [--check-base] [--dry-run]
  --dry-run   print everything, write nothing
"""
import glob, hashlib, math, os, re, subprocess, sys
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__)); CL = os.path.dirname(HERE); HS = os.path.dirname(CL); REPO = os.path.dirname(HS)
GEN3_WT = os.path.join(os.path.dirname(REPO), "HardDisks_gen3")
sys.path.insert(0, os.path.join(CL, "gen3_gate_261009")); sys.path.insert(0, os.path.join(HS, "validation")); sys.path.insert(0, HS)
import gen_testG_261009 as G
import tests_20260913 as T
BASE_R, K, MASSES, NSEED = 20261317, 100, (50, 100, 300), 3600
CID, REL, ENGINES = "epi8_H_H10_L10", "testR", ("gen2", "gen3")
TASKS = os.path.join(HERE, "tasks_testR_koa_261011.txt")
SHAF = os.path.join(HERE, "tasks_testR_koa_261011.sha256")
KF = os.path.join(HERE, "testR_chunk_size.txt")
RECF = os.path.join(HERE, "build_record_a643bc7.txt")
COMMIT = "a643bc7edb1b41574c0f9af0c5fea199dc797680"
K2SHA = os.path.join(CL, "gen3_koa_261009", "k2_binaries_a643bc7.sha256")
TESTG_KOA = os.path.join(HS, "experiments_gen3_gate_koa_261009", "testG")


def lines():
    geo = G.campaign_geo(); mi = T.A1_MASSES.index; L = []
    for r in range(NSEED):
        for M in MASSES:
            L0, H, Ns, st = geo[M]; s = T.run_seed(BASE_R, 0, mi(M), r)
            for e in ENGINES: L.append(f"B {REL}/{e}/{CID}/m_{M} {M} {r} {s} {L0} {H} {Ns} {st} {BASE_R} {e} main\n")
    return L


def record_text():
    m = re.search(r"^([0-9a-f]{64})  00ALLINONE$", open(K2SHA).read(), re.M)
    if not m: sys.exit("STOP: no 00ALLINONE line in k2_binaries_a643bc7.sha256")
    full = subprocess.run(["git", "-C", REPO, "rev-parse", COMMIT[:7]], capture_output=True, text=True).stdout.strip()
    if full != COMMIT: sys.exit(f"STOP: a643bc7 resolves to {full}")
    return (f"# the KOA binary ~/harddisks_gen3/hspist3/00ALLINONE (runsheet step 2, build_gen3_koa.sh, not rebuilt since); written "
            f"2026-10-11 on the Mac by gen_testR_koa_261011.py from the recorded values (k2_binaries_a643bc7.sha256; the commit)\n"
            f"commit {COMMIT}\nversion 00ALLINONE  git {COMMIT[:7]}  target koa\nsha256 {m.group(1)}\n")


def used_seeds():
    """every exact seed of every earlier task list (both trees), G.used_seeds()'s extra sets and the symmetry test's seeds"""
    used, nlists, kinds = set(), 0, Counter()
    for t in (HS, os.path.join(GEN3_WT, "hspist3")):
        for f in glob.glob(os.path.join(t, "cluster", "*", "tasks_*.txt")):
            if os.path.dirname(os.path.abspath(f)) == HERE: continue
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
    return used, nlists, kinds


def check_base():
    """the base as an integer token in any file: both working trees (rg), the tracked files of three heads (git grep); this package
    folder excluded"""
    pat = rf"(^|[^0-9.]){BASE_R}([^0-9.]|$)"
    ex = os.path.relpath(HERE, REPO)
    rg = subprocess.run(["rg", "-l", "--no-ignore", "--hidden", "-g", "!.git", "-g", f"!{ex}/**", "-e", pat, REPO, GEN3_WT],
                        capture_output=True, text=True)
    hits = [h for h in rg.stdout.splitlines() if "/cluster/gen3_testR_261011/" not in h]
    gg = subprocess.run(["git", "-C", REPO, "grep", "-E", "-l", pat, "main", "engine-gen3", "engine-divider-resched", "--"],
                        capture_output=True, text=True)
    ghits = [h for h in gg.stdout.splitlines() if "/cluster/gen3_testR_261011/" not in h]
    return hits, ghits


def per_task_seconds():
    """Test G on KOA: mean seconds per trajectory (the ##RUN header, worker included) per engine and mass, E5-2670 v3 nodes only"""
    out = {}
    for e in ENGINES:
        for M in (50, 300):
            txt = open(os.path.join(TESTG_KOA, e, CID, f"m_{M}", "run.log"), errors="replace").read()
            s = [int(x) for x in re.findall(r"^##RUN .* node \S+ \((\d+) s\)", txt, re.M)]
            out[(e, M)] = sum(s) / len(s)
    import pandas as pd                                            # the driver's nu_pred / x_M, from Test G's red_nu.csv (M = 50)
    cpred = float(pd.read_csv(os.path.join(TESTG_KOA, 'gen2', CID, 'm_50', 'red_nu.csv'))['nu_pred'].iloc[0]) / T.x_of(50, 10.0)
    tsim = lambda M: 800.0 + 200.0 / (cpred * T.x_of(M, 10.0))     # the hold (2000 x 0.4) + 200 predicted periods
    for e in ENGINES:   # M = 100: linear in the simulated time between M = 50 and 300
        a, b = out[(e, 50)], out[(e, 300)]
        out[(e, 100)] = a + (b - a) * (tsim(100) - tsim(50)) / (tsim(300) - tsim(50))
    return out, tsim


def design(L, used, nlists, kinds, base_hits):
    seeds = [int(l.split()[4]) for l in L]; per = Counter((int(l.split()[2]), l.split()[10]) for l in L)
    nch = -(-len(L) // K)
    print("# TEST R ON KOA -- the design tables of the registration (decision 28; 261012 sec. 4.7.57)\n")
    print(f"cell {CID} (N = 100, H = L0 = 10, pi/8); masses {', '.join(map(str, MASSES))}; {NSEED} exact seeds per mass (r = 0..{NSEED - 1}); "
          f"engines gen2 and gen3, the same seed in both; base seed {BASE_R}")
    geo = G.campaign_geo()
    print("trace stride by the Test G generator's rule (campaign_geo): " + ", ".join(f"M = {M}: {geo[M][3]}" for M in MASSES) + "\n")
    print("| M | engine | task lines | first seed (r = 0) |\n|---|---|---|---|")
    for M in MASSES:
        for e in ENGINES:
            s0 = T.run_seed(BASE_R, 0, T.A1_MASSES.index(M), 0)
            print(f"| {M} | {e} | {per[(M, e)]} | {s0} |")
    print(f"\ntotal {len(L)} task lines; order r by r (M = 50, 100, 300; gen2 then gen3)")
    print(f"distinct exact seeds: {len(set(seeds))} of {len(seeds) // 2} (each in both engines): {'all distinct' if len(set(seeds)) == len(seeds) // 2 else 'NOT DISTINCT'}")
    ov = len(set(seeds) & used)
    print(f"overlap with every earlier task list ({nlists} lists in cluster/*/ of main and the engine-gen3 worktree; line types "
          f"{dict(kinds)}), the extra seed sets of the Test G generators and the symmetry test's seeds ({len(used)} seeds in all): {ov}"
          + (" -- zero overlap with every earlier task list" if ov == 0 else " -- OVERLAP"))
    if base_hits is not None:
        h, g = base_hits
        print(f"the base seed {BASE_R} as an integer token: {len(h)} files in both working trees, {len(g)} tracked files in main, "
              f"engine-gen3, engine-divider-resched (the Test R package folder excluded)" + ("" if (h or g) else " -- it occurs in no file"))
    T_, tsim = per_task_seconds()
    print(f"\n## chunks: K = {K} task lines per array task; {nch} chunks; array 1-{nch} (KOA MaxArraySize 25001); throttle %64 (128 only "
          f"while sinfo -p shared -s shows >= 10 % of the nodes idle)\n")
    print("| engine | M | simulated time [sigma-time] | seconds per trajectory on E5 (Test G, mean) | source |\n|---|---|---|---|---|")
    for e in ENGINES:
        for M in MASSES:
            print(f"| {e} | {M} | {tsim(M):.0f} | {T_[(e, M)]:.2f} | {'Test G on KOA, 400 trajectories' if M != 100 else 'interpolated in the simulated time'} |")
    ch = [0.0] * nch
    for i, l in enumerate(L): ch[i // K] += T_[(l.split()[10], int(l.split()[2]))]
    tot = sum(ch)
    print(f"\npredicted chunk time on E5: mean {sum(ch) / nch / 60:.1f} min, min {min(ch) / 60:.1f}, max {max(ch) / 60:.1f} min "
          f"(every chunk holds 16-17 trajectories of each of the six arms); the limit asked: <= 15 min")
    print(f"expected core-hours: {tot / 3600:.1f} (E5; Gold nodes are faster, E5/Gold 0.709); requested: {nch} chunks x 1 h = {nch} "
          f"core-hours (testR_koa.sbatch --time=1:00:00)")
    print(f"wall time at %64: about {math.ceil(nch / 64)} waves of {max(ch) / 60:.0f} min, if 64 cores are free at once")
    print("data on KOA scratch, expected: per trajectory the trace (about 0.86 MB until the merge compresses it to about 0.1 MB), gen3's "
          "psi6(t) file (0.15-0.3 MB, compressed by the worker) and the run folder (about 45 kB); about "
          f"{len(L) * 0.86 / 1024 + len(L) / 2 * 0.22 / 1024 + len(L) * 0.045 / 1024:.0f} GB before the merge, about "
          f"{len(L) * 0.10 / 1024 + len(L) / 2 * 0.22 / 1024 + len(L) * 0.045 / 1024:.0f} GB after it")


def main():
    L = lines(); txt = "".join(L); sha = hashlib.sha256(txt.encode()).hexdigest(); rec = record_text()
    for p, t in ((TASKS, txt), (SHAF, sha + "\n"), (KF, f"{K}\n"), (RECF, rec)):
        if os.path.exists(p) and open(p).read() != t: sys.exit(f"STOP: {p} exists with other content -- not overwriting")
    base_hits = check_base() if "--check-base" in sys.argv else None
    if base_hits and (base_hits[0] or base_hits[1]) and not os.path.exists(TASKS):
        sys.exit(f"STOP: the base seed {BASE_R} occurs in {base_hits[0][:5]} {base_hits[1][:5]}")
    used, nlists, kinds = used_seeds()
    seeds = [int(l.split()[4]) for l in L]
    if len(set(seeds)) != len(seeds) // 2 or set(seeds) & used: sys.exit("STOP: the seeds are not all distinct or overlap an earlier list")
    if "--dry-run" not in sys.argv:
        for p, t in ((TASKS, txt), (SHAF, sha + "\n"), (KF, f"{K}\n"), (RECF, rec)):
            if not os.path.exists(p): open(p, "w").write(t)
    print(f"{os.path.relpath(TASKS, HS)}: {len(L)} lines, SHA-256 {sha}; chunk size {K} -> {-(-len(L) // K)} chunks; "
          f"{os.path.relpath(RECF, HS)}: sha256 {rec.split('sha256 ')[1].strip()[:12]}..., version '{rec.split('version ')[1].splitlines()[0]}'\n")
    if "--design" in sys.argv: design(L, used, nlists, kinds, base_hits)


if __name__ == "__main__":
    main()
