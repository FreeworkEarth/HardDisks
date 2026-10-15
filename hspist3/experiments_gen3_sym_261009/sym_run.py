#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.28; plan-author decision 13, item 1): the runner of the SYMMETRY TEST, registered with its
seed list and its verdict script (sym_verdict.py) in one commit before any run; neither is edited afterwards.
Cells (the plan author's item 1): N = 400, H = 20, M4 seeding (--gen3-seeding=lattice, the default jitter), the exact box, the
divider held: eta = pi/8 (L0 = 20, as part 3a's held cell) and eta = 0.70 (L0 = N pi / (8 H eta), P1's cell). The command is P1's
(experiments_gen3_p1_261009/p1_run.py cmd(): M = 300, psi6 every 1 sigma-time, a 1-sigma-time record after the hold) with the hold
set to 5500 sigma-time (330000 steps of 1/60 sigma-time), --seed=20261222 and --gen3-virial-blocks=50 (110 hold blocks; the
verdict discards the first 10 = 500 sigma-time).
SEEDS: a new base never used before, BASE = 20261222 (no file of main, engine-gen3 or this session's scratch contained the number
before this script); run_seed(BASE, 0, k, r) with k = 0 (pi/8), 1 (0.70); r = 0..39 the main set, r = 40..79 the extension
(run only for a cell whose 40 give 2 <= |t| < 3). --write-seeds writes sym_seeds.txt and checks the 160 seeds against every
earlier seed list and today's bases (stage C 20261216, P1 20261217, P2 20261218, the KOA package 20261219/20261220, part 3a
20261221, Test G 20261209).
BINARY: the frozen committed-tree binary of engine-gen3 5d19aa7 (experiments_gen3_p3a_261009/bin_3a/00ALLINONE, SHA-256 in
build_record_5d19aa7.txt); the engine is unchanged since. At most --jobs processes (default 10, programme rule 1: <= 12). Every
run in its own folder <out>/<cell>/seed<r>/ (command.txt, run.log, exit_code.txt); nothing is compressed, moved or deleted; a run
with exit != 0, clean != 1 or no run record is a finding, reported here and never rerun (the verdict's inventory makes it INVALID).
usage (from hspist3/ of the engine-gen3 worktree):
  python3 experiments_gen3_sym_261009/sym_run.py --write-seeds
  python3 experiments_gen3_sym_261009/sym_run.py --bin <00ALLINONE> --expect-bin-sha <sha256> --out <data root> --set main|ext [--cells pi8,070]
"""
import argparse, glob, hashlib, math, os, re, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(HS, "experiments_gen3_p1_261009")); sys.path.insert(0, os.path.join(HS, "experiments_gen3_h_261009"))
sys.path.insert(0, os.path.join(HS, "experiments_gen3_p2_261009")); sys.path.insert(0, os.path.join(HS, "cluster", "gen3_gate_261009"))
import p1_run as P1
T = P1.T
BASE, NMAIN, HOLD_STEPS, VB = 20261222, 40, 330000, 50
SEEDS = os.path.join(HERE, "sym_seeds.txt")


def cells():
    p1 = {(c["N"], round(c["eta"], 3)): c for c in P1.cells()}
    c070 = dict(p1[(400, 0.70)])
    return [dict(k=0, name="pi8", N=400, H=20.0, L0=20.0, eta=math.pi / 8),
            dict(k=1, name="070", N=400, H=c070["H"], L0=c070["L0"], eta=c070["eta"])]


def seed_rows():
    return [(c["name"], r, "main" if r < NMAIN else "ext", T.run_seed(BASE, 0, c["k"], r)) for c in cells() for r in range(2 * NMAIN)]


def seed_text():
    return "".join(f"{n} {r} {s} {seed}\n" for n, r, s, seed in seed_rows())


def cmd(binp, d, c, seed):
    out = []
    for a in P1.cmd(binp, d, c, seed):
        if a.startswith("--wall-hold-steps="): a = f"--wall-hold-steps={HOLD_STEPS}"
        elif a.startswith("--seed="): a = f"--seed={BASE}"
        out.append(a)
    return out + [f"--gen3-virial-blocks={VB}"]


def earlier_seeds():
    """every exact seed of the earlier task lists, and of today's experiment bases"""
    import gen_testG_261009 as G
    used = G.used_seeds()                                              # campaign B/A/AF, Test T, T-prime, stage A
    for f in glob.glob(os.path.join(HS, "cluster", "gen3_gate_261009", "tasks_*.txt")):
        for l in open(f):
            p = l.split()
            if p and p[0] == "B": used.add(int(p[4]))
            elif p and p[0] == "AF": used.add(int(p[3]))
    used |= {s for c in P1.cells() for s in c["seeds"]}                 # P1 (20261217)
    used |= {T.run_seed(b, l, k, r) for b in (20261216, 20261218, 20261219, 20261220, 20261221) for l in range(4)
             for k in range(40) for r in range(800 if b == 20261220 else 50)}   # stage C, P2, the KOA package, part 3a (supersets)
    used |= {T.run_seed(20261219, 1, N, r) for N in (100, 400, 900, 1600) for r in range(8)}   # g3_rate_ld_koa.py: k = N
    return used


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--write-seeds", action="store_true"); ap.add_argument("--bin"); ap.add_argument("--expect-bin-sha")
    ap.add_argument("--out"); ap.add_argument("--set", choices=("main", "ext")); ap.add_argument("--cells", default="070,pi8")
    ap.add_argument("--jobs", type=int, default=10)
    a = ap.parse_args()
    if a.write_seeds:
        txt = seed_text()
        if os.path.exists(SEEDS) and open(SEEDS).read() != txt: sys.exit(f"STOP: {SEEDS} exists with other content")
        open(SEEDS, "w").write(txt)
        seeds = [r[3] for r in seed_rows()]; used = earlier_seeds()
        print(f"{os.path.relpath(SEEDS, HS)}: {len(seeds)} seeds ({NMAIN} main + {NMAIN} extension per cell), distinct {len(set(seeds))}; "
              f"overlap with the earlier seeds ({len(used)}): {len(set(seeds) & used)}; SHA-256 {hashlib.sha256(txt.encode()).hexdigest()}")
        return
    if a.jobs > 12: sys.exit("STOP: at most 12 processes (programme rule 1)")
    if not (a.bin and a.expect_bin_sha and a.out and a.set): sys.exit("--bin, --expect-bin-sha, --out and --set are required")
    if open(SEEDS).read() != seed_text(): sys.exit("STOP: sym_seeds.txt is not this script's seed list")
    binp, O = os.path.abspath(a.bin), os.path.abspath(a.out); os.makedirs(O, exist_ok=True)
    if hashlib.sha256(open(binp, "rb").read()).hexdigest() != a.expect_bin_sha: sys.exit("STOP: binary SHA-256")
    build = subprocess.run([binp, "--version"], capture_output=True, text=True).stdout.splitlines()[0]
    log = open(os.path.join(O, "runner.log"), "a")
    def say(s): print(s, flush=True); log.write(s + "\n"); log.flush()
    say(f"# {time.strftime('%Y-%m-%d %H:%M:%S %Z')} symmetry test, set {a.set}, cells {a.cells}: binary {binp} ({a.expect_bin_sha}; {build}); out {O}; jobs {a.jobs}")
    say(subprocess.run(["df", "-h", O], capture_output=True, text=True).stdout.strip())
    want = a.cells.split(","); C = {c["name"]: c for c in cells()}
    jobs = [(C[n], r, seed, os.path.join(O, n, f"seed{r}")) for w in want for n, r, s, seed in seed_rows() if n == w and s == a.set]
    for c, r, seed, d in jobs:
        if os.path.exists(d): sys.exit(f"STOP: {d} exists -- a run is never repeated into the same folder")
    run, t0 = [], time.time()
    def reap():
        for j in list(run):
            p, c, r, d, ts = j
            if p.poll() is None: continue
            run.remove(j); open(os.path.join(d, "exit_code.txt"), "w").write(f"{p.returncode}\n")
            lg = open(os.path.join(d, "run.log"), errors="ignore").read()
            recs = re.findall(r"^\[EDMD3-HEALTH\] .*?: clean=(\S+) ", lg, re.M); nb = len(re.findall(r"^\[EDMD3\] built #", lg, re.M))
            bad = "" if (recs == ["1"] and nb == 1 and p.returncode == 0) else f"exit {p.returncode}, run record {recs}, builds {nb}"
            say(f"  {c['name']} seed{r}: exit {p.returncode}, {time.time() - ts:.0f} s" + (f"  **FINDING: {bad}**" if bad else ""))
    for c, r, seed, d in jobs:
        while len(run) >= a.jobs: reap(); time.sleep(0.5)
        os.makedirs(d)
        cm = cmd(binp, d, c, seed)
        with open(os.path.join(d, "command.txt"), "w") as fh: fh.write(" ".join(cm) + "\n")
        run.append((subprocess.Popen(cm, cwd=d, stdout=open(os.path.join(d, "run.log"), "w"), stderr=subprocess.STDOUT), c, r, d, time.time()))
    while run: reap(); time.sleep(0.5)
    say(f"# done {time.strftime('%Y-%m-%d %H:%M:%S %Z')} after {time.time() - t0:.0f} s")
    say(subprocess.run(["df", "-h", O], capture_output=True, text=True).stdout.strip())


if __name__ == "__main__":
    main()
