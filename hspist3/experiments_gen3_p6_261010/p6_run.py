#!/usr/bin/env python3
"""##CHRIS 2026-10-10 (261012 sec. 4.7.43; plan-author decision 18, item 3): PILOT P6, the turned-lattice start -- a pilot, no paper
use. The M4c start --gen3-seeding=lattice --gen3-lattice-orientation=turned (close-packed lines along y; psi6 orientation +-30 deg),
the divider held 2e4 sigma-time (1,200,000 steps), --psi6-every=1, --gen3-psi6-compartments, --gen3-virial-blocks=50, 6 seeds, in
the three dense cells (P4's construction: L0 = N pi / (8 H eta) exactly, --gen3-exact-box, M = 300, a 1-sigma-time record):
  k = 0  N = 400  (H = 20), eta 0.740
  k = 1  N = 1600 (H = 40), eta 0.740
  k = 2  N = 1600 (H = 40), eta 0.780
A cell infeasible for the turned family stops at the seeding (exit 2, [EDMD3-M4] INFEASIBLE) and is reported and skipped.
ADDITION (not ordered; labelled in every table): the explicit 0-deg start (--gen3-lattice-orientation=0) in the two N = 1600 cells,
6 seeds each. Reason (sec. 4.7.42, found before these runs): at N = 1600, eta 0.74 and 0.78, the present M4 family -- P4's start A --
is the turned crystal (P4's snapshots read +-30 deg at t = 0), so without an explicit 0-deg start the question "which orientation
has the lower Z" has no 0-deg side at N = 1600. At N = 400 the present family is the 0-deg crystal (P4a, sec. 4.7.37).
Snapshots every 2000 sigma-time for seeds 0 and 1 of each start, as P4 and P5.
SEEDS: new and independent, BASE = 20261228 (no script, command or seed list holds it; its only whole-word occurrences are time
values in one Szilard trace CSV); run_seed(BASE, l, k, r), l = 0 the turned start, l = 1 the 0-deg addition, r = 0..5.
--check-seeds prints them against every earlier seed (p5_run.earlier_seeds, P5's seeds, the M4c grid of d18_evidence).
BINARY: the committed-tree build of engine-gen3 f72d764 (experiments_gen3_d18_261010/bin_f72d764; SHA-256 checked).
RULES (decision 18): at most 8 simulation processes (the cap is the smaller of --jobs and the integer in --cap-file, read before every
launch); >= 4 GiB free before every launch, else nothing more is launched; the longest runs first (N = 1600, eta 0.78); every run in
its own folder <out>/<cell>/<start>/seed<r>/; a run with exit != 0 (other than an infeasible seeding), clean != 1 or no run record
is a finding, reported and never rerun.
usage (from hspist3/ of the engine-gen3 worktree):
  python3 experiments_gen3_p6_261010/p6_run.py --check-seeds
  python3 experiments_gen3_p6_261010/p6_run.py --bin <00ALLINONE> --expect-bin-sha <sha256> --out <data root> [--jobs 8] [--cap-file F]
"""
import argparse, hashlib, math, os, re, shutil, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
for sub in ("experiments_gen3_p5_261010", "experiments_gen3_p4_261010", "experiments_gen3_d18_261010", "experiments_gen3_p1_261009",
            "experiments_gen3_h_261009", "experiments_gen3_m4_261009", "experiments_gen3_p2_261009", os.path.join("cluster", "gen3_gate_261009")):
    sys.path.insert(0, os.path.join(HS, sub))
import p4_run as P4
T = P4.T
BASE, NSEED, HOLD, FREE_MIN = 20261228, 6, 1200000, 4 * 1024 ** 3
STARTS = (("turned", "turned", 0), ("zero", "0", 1))              # (label, --gen3-lattice-orientation value, l)


def cells():
    return [dict(k=0, N=400, H=20.0, eta=0.740), dict(k=1, N=1600, H=40.0, eta=0.740), dict(k=2, N=1600, H=40.0, eta=0.780)]


def cell(c):
    c = dict(c); c["L0"] = c["N"] * math.pi / (8.0 * c["H"] * c["eta"]); c["hold"] = HOLD; c["name"] = f"N{c['N']}_eta{c['eta']:.3f}"; return c


def starts_for(c): return [s for s in STARTS if s[0] == "turned" or c["N"] == 1600]


def seed(c, l, r): return T.run_seed(BASE, l, c["k"], r)


def cmd(binp, d, c, orient, s, r):
    out = []
    for a in P4.cmd(binp, d, c, "lattice", s, r):
        if a.startswith("--seed="): a = f"--seed={BASE}"
        out.append(a)
    return out + [f"--gen3-lattice-orientation={orient}", "--gen3-psi6-compartments"]


def check_seeds():
    import p5_run as P5
    import d18_evidence as D18
    S = [seed(cell(c), l, r) for c in cells() for lab, o, l in starts_for(cell(c)) for r in range(NSEED)]
    used = P5.earlier_seeds()
    used |= {P5.seed(P5.cell(c), r) for c in P5.cells() for r in range(P5.NSEED)} | {T.run_seed(P5.BASE, 2, 0, 0)}
    used |= {g[k] for g in D18.grid() for k in ("s1", "s2", "sd", "s0")}
    print(f"P6 seeds: {len(S)}, distinct {len(set(S))}; earlier seeds: {len(used)}; overlap: {len(set(S) & used)}")
    for c in map(cell, cells()):
        for lab, o, l in starts_for(c): print(f"  {c['name']} {lab}: " + " ".join(str(seed(c, l, r)) for r in range(NSEED)))
    return len(set(S)) == len(S) and not (set(S) & used)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--check-seeds", action="store_true"); ap.add_argument("--bin"); ap.add_argument("--expect-bin-sha")
    ap.add_argument("--out"); ap.add_argument("--jobs", type=int, default=8); ap.add_argument("--cap-file")
    a = ap.parse_args()
    if a.check_seeds: sys.exit(0 if check_seeds() else "SEED CHECK FAILED")
    if a.jobs > 8: sys.exit("STOP: at most 8 simulation processes (decision 18)")
    if not (a.bin and a.expect_bin_sha and a.out): sys.exit("--bin, --expect-bin-sha and --out are required")
    binp, O = os.path.abspath(a.bin), os.path.abspath(a.out); os.makedirs(O, exist_ok=True)
    if hashlib.sha256(open(binp, "rb").read()).hexdigest() != a.expect_bin_sha: sys.exit("STOP: binary SHA-256")
    build = subprocess.run([binp, "--version"], capture_output=True, text=True).stdout.splitlines()[0]
    log = open(os.path.join(O, "runner.log"), "a")
    def say(s): print(s, flush=True); log.write(s + "\n"); log.flush()
    def cap():
        try: return max(0, min(a.jobs, 8, int(open(a.cap_file).read().strip()))) if a.cap_file else a.jobs
        except (OSError, ValueError): return a.jobs
    say(f"# {time.strftime('%Y-%m-%d %H:%M:%S %Z')} P6: binary {binp} ({a.expect_bin_sha}; {build}); out {O}; jobs {a.jobs}; cap file {a.cap_file}")
    say(subprocess.run(["df", "-h", O], capture_output=True, text=True).stdout.strip())
    say(subprocess.run(["sysctl", "vm.swapusage"], capture_output=True, text=True).stdout.strip())
    C = {c["k"]: cell(c) for c in cells()}
    jobs = [(C[k], lab, o, r, seed(C[k], l, r), os.path.join(O, C[k]["name"], lab, f"seed{r}")) for k in (2, 1, 0)
            for lab, o, l in starts_for(C[k]) for r in range(NSEED)]
    for j in jobs:
        if os.path.exists(j[5]): sys.exit(f"STOP: {j[5]} exists -- a run is never repeated into the same folder")
    run, t0, stopped = [], time.time(), None
    def reap():
        for j in list(run):
            p, c, lab, r, d, ts = j
            if p.poll() is None: continue
            run.remove(j); open(os.path.join(d, "exit_code.txt"), "w").write(f"{p.returncode}\n")
            lg = open(os.path.join(d, "run.log"), errors="ignore").read()
            if "[EDMD3-M4] INFEASIBLE" in lg and p.returncode == 2:
                say(f"  {c['name']} {lab} seed{r}: INFEASIBLE for this orientation (exit 2) -- the cell is skipped"); continue
            recs = re.findall(r"^\[EDMD3-HEALTH\] .*?: clean=(\S+) ", lg, re.M); nb = len(re.findall(r"^\[EDMD3\] built #", lg, re.M))
            bad = "" if (recs == ["1"] and nb == 1 and p.returncode == 0) else f"exit {p.returncode}, run record {recs}, builds {nb}"
            say(f"  {c['name']} {lab} seed{r}: exit {p.returncode}, {time.time() - ts:.0f} s" + (f"  **FINDING: {bad}**" if bad else ""))
    for c, lab, o, r, s, d in jobs:
        while len(run) >= cap(): reap(); time.sleep(1.0)
        free = shutil.disk_usage(O).free
        if free < FREE_MIN:
            stopped = f"free space {free / 1024 ** 3:.2f} GiB < 4 GiB before {c['name']} {lab} seed{r}"
            say(f"**STOPPED LAUNCHING: {stopped}; the runs already going finish; nothing else is launched**"); break
        os.makedirs(d)
        cm = cmd(binp, d, c, o, s, r)
        with open(os.path.join(d, "command.txt"), "w") as fh: fh.write(" ".join(cm) + "\n")
        run.append((subprocess.Popen(cm, cwd=d, stdout=open(os.path.join(d, "run.log"), "w"), stderr=subprocess.STDOUT), c, lab, r, d, time.time()))
        say(f"  start {c['name']} {lab} seed{r} (running {len(run)}, cap {cap()}, free {free / 1024 ** 3:.1f} GiB)")
    while run: reap(); time.sleep(1.0)
    say(f"# done {time.strftime('%Y-%m-%d %H:%M:%S %Z')} after {time.time() - t0:.0f} s" + (f"; STOPPED LAUNCHING ({stopped})" if stopped else ""))
    say(subprocess.run(["df", "-h", O], capture_output=True, text=True).stdout.strip())
    say(subprocess.run(["sysctl", "vm.swapusage"], capture_output=True, text=True).stdout.strip())


if __name__ == "__main__":
    main()
