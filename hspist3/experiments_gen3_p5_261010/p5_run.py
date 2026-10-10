#!/usr/bin/env python3
"""##CHRIS 2026-10-10 (261012 sec. 4.7.38; plan-author decision 16, item 3): PILOT P5, the compression speed -- a pilot, no paper use.
Start B only (M4b, --gen3-seeding=compressed: the lattice at eta_0 0.60 in a longer box, T_melt 500 sigma-time, then the pistons at u
sigma per sigma-time until each compartment has the target L0, the KE rescaled), with slower pistons than P4's u = 0.01; then the
divider held 2e4 sigma-time (1,200,000 steps). 6 seeds per cell:
  k = 0  N = 1600 (H = 40), eta 0.740, u = 0.001
  k = 1  N = 1600 (H = 40), eta 0.740, u = 0.0001
  k = 2  N = 1600 (H = 40), eta 0.780, u = 0.001
  k = 3  N = 400  (H = 20), eta 0.740, u = 0.001
Every run: P4's command (p4_run.cmd: P1's command with M = 300, --psi6-every=1, --gen3-virial-blocks=50, --gen3-snapshots=2000 for seeds
0 and 1, a 1-sigma-time record; L0 = N pi / (8 H eta) exactly, --gen3-exact-box) with start B and --gen3-compress-speed=u.
SEEDS: new and independent, BASE = 20261226 (no script, command or seed list holds it; its only whole-word occurrences before this
script are time values, 20261226.16..., in one Szilard trace CSV); run_seed(BASE, 1, k, r) (l = 1 for start B, as in P4), k the cell
index, r = 0..5. --check-seeds prints them against every earlier seed (the symmetry test's list of earlier seeds, its own 160, P3, the
M4b evidence, P4 and the KOA Test G list).
SMOKE (--smoke, before the pilot; not part of P5's data): one small run with the slowest piston, N = 100 (H = 10), eta 0.740,
u = 0.0001, held 10 sigma-time, seed run_seed(BASE, 2, 0, 0): that a 100-times-slower compression runs through and is clean.
BINARY: the frozen committed-tree build of engine-gen3 06e7248 (experiments_gen3_m4b_261009/bin_m4b; SHA-256 checked), as P4.
RULES (decision 16): at most 8 simulation processes -- the cap is the smaller of --jobs (<= 8) and the integer in --cap-file, read
before every launch; before every launch the free space of the data's volume must be >= 4 GiB, else nothing more is launched (the runs
already going finish) and the runner says so. The longest runs first (k = 1, then 2, 0, 3). Every run in its own folder
<out>/<cell>/B/seed<r>/; nothing is compressed, moved or deleted; a run with exit != 0, clean != 1 or no run record is a finding,
reported and never rerun.
usage (from hspist3/ of the engine-gen3 worktree):
  python3 experiments_gen3_p5_261010/p5_run.py --check-seeds
  python3 experiments_gen3_p5_261010/p5_run.py --bin <00ALLINONE> --expect-bin-sha <sha256> --out <data root> [--smoke] [--jobs 8] [--cap-file F]
"""
import argparse, glob, hashlib, math, os, re, shutil, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
for sub in ("experiments_gen3_p4_261010", "experiments_gen3_p1_261009", "experiments_gen3_h_261009", "experiments_gen3_p2_261009",
            os.path.join("cluster", "gen3_gate_261009")):
    sys.path.insert(0, os.path.join(HS, sub))
import p4_run as P4
T = P4.T
BASE, NSEED, HOLD, FREE_MIN = 20261226, 6, 1200000, 4 * 1024 ** 3
ORDER = (1, 2, 0, 3)                                                   # the longest runs first


def cells():
    return [dict(k=0, N=1600, H=40.0, eta=0.740, u=0.001), dict(k=1, N=1600, H=40.0, eta=0.740, u=0.0001),
            dict(k=2, N=1600, H=40.0, eta=0.780, u=0.001), dict(k=3, N=400, H=20.0, eta=0.740, u=0.001)]


def cell(c):
    c = dict(c); c["L0"] = c["N"] * math.pi / (8.0 * c["H"] * c["eta"]); c["hold"] = HOLD
    c["name"] = f"N{c['N']}_eta{c['eta']:.3f}_u{c['u']:g}"; c["p4name"] = f"N{c['N']}_eta{c['eta']:.3f}"; return c


def smoke_cell():
    c = cell(dict(k=0, N=100, H=10.0, eta=0.740, u=0.0001)); c["hold"] = 600; c["name"] = "smoke_N100_eta0.740_u0.0001"; return c


def seed(c, r): return T.run_seed(BASE, 1, c["k"], r)


def cmd(binp, d, c, s, r):
    out = []
    for a in P4.cmd(binp, d, c, "compressed", s, r):
        if a.startswith("--seed="): a = f"--seed={BASE}"
        out.append(a)
    return out + [f"--gen3-compress-speed={c['u']:g}"]


def earlier_seeds():
    """every earlier seed: the symmetry test's list of earlier seeds, its own 160, P3, the M4b evidence, P4, the KOA Test G list"""
    sys.path.insert(0, os.path.join(HS, "experiments_gen3_sym_261009")); sys.path.insert(0, os.path.join(HS, "experiments_gen3_p3_261009"))
    sys.path.insert(0, os.path.join(HS, "experiments_gen3_m4b_261009"))
    import sym_run, p3_run, m4b_evidence
    used = sym_run.earlier_seeds()
    used |= {int(l.split()[3]) for l in open(os.path.join(HS, "experiments_gen3_sym_261009", "sym_seeds.txt")) if l.strip()}
    used |= {p3_run.seed(c, l, r) for c in p3_run.cells() for l in (0, 1) for r in range(p3_run.NSEED)}
    used |= {c["seed"] for c in m4b_evidence.cells()}
    used |= {P4.seed(P4.cell(c), l, r) for c in P4.cells() for l in (0, 1) for r in range(P4.NSEED)}
    for l in open(os.path.join(HS, "cluster", "gen3_koa_261009", "tasks_testG_koa_261009.txt")):
        p = l.split()
        if p and p[0] == "B": used.add(int(p[4]))
    return used


def check_seeds():
    S = [seed(cell(c), r) for c in cells() for r in range(NSEED)]; sm = T.run_seed(BASE, 2, 0, 0)
    used = earlier_seeds()
    print(f"P5 seeds: {len(S)} (4 cells x {NSEED}), distinct {len(set(S))}; the smoke seed {sm} among them: {sm in S}")
    print(f"earlier seeds: {len(used)}; overlap with P5's: {len(set(S) & used)}; the smoke seed among them: {sm in used}")
    for c in map(cell, cells()): print(f"  {c['name']}: " + " ".join(str(seed(c, r)) for r in range(NSEED)))
    return len(set(S)) == len(S) and not (set(S) & used) and sm not in used and sm not in S


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--check-seeds", action="store_true"); ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--bin"); ap.add_argument("--expect-bin-sha"); ap.add_argument("--out"); ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--cap-file")
    a = ap.parse_args()
    if a.check_seeds: sys.exit(0 if check_seeds() else "SEED CHECK FAILED")
    if a.jobs > 8: sys.exit("STOP: at most 8 simulation processes (decision 16)")
    if not (a.bin and a.expect_bin_sha and a.out): sys.exit("--bin, --expect-bin-sha and --out are required")
    binp, O = os.path.abspath(a.bin), os.path.abspath(a.out); os.makedirs(O, exist_ok=True)
    if hashlib.sha256(open(binp, "rb").read()).hexdigest() != a.expect_bin_sha: sys.exit("STOP: binary SHA-256")
    build = subprocess.run([binp, "--version"], capture_output=True, text=True).stdout.splitlines()[0]
    log = open(os.path.join(O, "runner.log"), "a")
    def say(s): print(s, flush=True); log.write(s + "\n"); log.flush()
    def cap():
        try: return max(0, min(a.jobs, 8, int(open(a.cap_file).read().strip()))) if a.cap_file else a.jobs
        except (OSError, ValueError): return a.jobs
    say(f"# {time.strftime('%Y-%m-%d %H:%M:%S %Z')} P5{' SMOKE' if a.smoke else ''}: binary {binp} ({a.expect_bin_sha}; {build}); out {O}; "
        f"jobs {a.jobs}; cap file {a.cap_file}")
    say(subprocess.run(["df", "-h", O], capture_output=True, text=True).stdout.strip())
    say(subprocess.run(["sysctl", "vm.swapusage"], capture_output=True, text=True).stdout.strip())
    if a.smoke:
        c = smoke_cell(); jobs = [(c, 0, T.run_seed(BASE, 2, 0, 0), os.path.join(O, c["name"], "B", "seed0"))]
    else:
        C = {c["k"]: cell(c) for c in cells()}
        jobs = [(C[k], r, seed(C[k], r), os.path.join(O, C[k]["name"], "B", f"seed{r}")) for k in ORDER for r in range(NSEED)]
    for j in jobs:
        if os.path.exists(j[3]): sys.exit(f"STOP: {j[3]} exists -- a run is never repeated into the same folder")
    run, t0, stopped = [], time.time(), None
    def reap():
        for j in list(run):
            p, c, r, d, ts = j
            if p.poll() is None: continue
            run.remove(j); open(os.path.join(d, "exit_code.txt"), "w").write(f"{p.returncode}\n")
            lg = open(os.path.join(d, "run.log"), errors="ignore").read()
            recs = re.findall(r"^\[EDMD3-HEALTH\] .*?: clean=(\S+) ", lg, re.M); nb = len(re.findall(r"^\[EDMD3\] built #", lg, re.M))
            m4b = len(re.findall(r"^\[EDMD3-M4B\] L0=", lg, re.M))
            bad = "" if (recs == ["1"] and nb == 1 and m4b == 1 and p.returncode == 0) else f"exit {p.returncode}, run record {recs}, builds {nb}, M4b lines {m4b}"
            say(f"  {c['name']} B seed{r}: exit {p.returncode}, {time.time() - ts:.0f} s" + (f"  **FINDING: {bad}**" if bad else ""))
    for c, r, s, d in jobs:
        while len(run) >= cap(): reap(); time.sleep(1.0)
        free = shutil.disk_usage(O).free
        if free < FREE_MIN:
            stopped = f"free space {free / 1024 ** 3:.2f} GiB < 4 GiB before {c['name']} seed{r}"
            say(f"**STOPPED LAUNCHING: {stopped}; the runs already going finish; nothing else is launched**"); break
        os.makedirs(d)
        cm = cmd(binp, d, c, s, r)
        with open(os.path.join(d, "command.txt"), "w") as fh: fh.write(" ".join(cm) + "\n")
        run.append((subprocess.Popen(cm, cwd=d, stdout=open(os.path.join(d, "run.log"), "w"), stderr=subprocess.STDOUT), c, r, d, time.time()))
        say(f"  start {c['name']} B seed{r} (running {len(run)}, cap {cap()}, free {free / 1024 ** 3:.1f} GiB)")
    while run: reap(); time.sleep(1.0)
    say(f"# done {time.strftime('%Y-%m-%d %H:%M:%S %Z')} after {time.time() - t0:.0f} s" + (f"; STOPPED LAUNCHING ({stopped})" if stopped else ""))
    say(subprocess.run(["df", "-h", O], capture_output=True, text=True).stdout.strip())
    say(subprocess.run(["sysctl", "vm.swapusage"], capture_output=True, text=True).stdout.strip())


if __name__ == "__main__":
    main()
