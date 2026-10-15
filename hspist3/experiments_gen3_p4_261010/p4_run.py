#!/usr/bin/env python3
"""##CHRIS 2026-10-10 (261012 sec. 4.7.34; plan-author decision 15, item 2): PILOT P4 -- a pilot, no paper use. Two starts, as P3
(experiments_gen3_p3_261009): A = the M4 lattice at the target density (--gen3-seeding=lattice); B = M4b, the compressed fluid
(--gen3-seeding=compressed with its defaults: eta_0 0.60, T_melt 500 sigma-time, u 0.01 sigma per sigma-time). 6 seeds per start and cell.
  a. N = 400 (H = 20) at eta 0.740, the divider held 1e5 sigma-time (6,000,000 steps);
  b. N = 1600 (H = 40) at eta 0.740 and 0.780, the divider held 2e4 sigma-time (1,200,000 steps).
Every run: --psi6-every=1, --gen3-virial-blocks=50; --gen3-snapshots=2000 for seeds 0 and 1 of each start; then a 1-sigma-time record
(the speed-of-sound loop needs one). P1's command otherwise (M = 300); L0 = N pi / (8 H eta) exactly (--gen3-exact-box).
SEEDS: new and independent, BASE = 20261225 (in no file before this script); run_seed(BASE, l, k, r) with l = 0 for start A and l = 1 for
start B, k the cell index, r = 0..5.
BINARY: the frozen committed-tree build of engine-gen3 06e7248 (experiments_gen3_m4b_261009/bin_m4b; SHA-256 checked); the engine and
driver are unchanged since.
RULES (decision 15): at most 8 simulation processes -- the cap is the smaller of --jobs (<= 8) and the integer in --cap-file, read before
every launch, so slots can be freed for another test without passing 8 in all; before every launch the free space of the data's volume
must be >= 4 GiB, else nothing more is launched (the runs already going finish) and the runner says so. The longest runs first (part a).
Every run in its own folder <out>/<cell>/<start>/seed<r>/; nothing is compressed, moved or deleted; a run with exit != 0, clean != 1 or
no run record is a finding, reported and never rerun.
usage (from hspist3/ of the engine-gen3 worktree):
  python3 experiments_gen3_p4_261010/p4_run.py --bin <00ALLINONE> --expect-bin-sha <sha256> --out <data root> [--jobs 8] [--cap-file F]
"""
import argparse, hashlib, math, os, re, shutil, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(HS, "experiments_gen3_p1_261009"))
import p1_run as P1
T = P1.T
BASE, NSEED, SNAP_SEEDS, SNAP_EVERY, FREE_MIN = 20261225, 6, (0, 1), 2000, 4 * 1024 ** 3
STARTS = (("A", "lattice"), ("B", "compressed"))


def cells():
    return [dict(k=0, N=400, H=20.0, eta=0.740, hold=6000000, part="a"),
            dict(k=1, N=1600, H=40.0, eta=0.740, hold=1200000, part="b"),
            dict(k=2, N=1600, H=40.0, eta=0.780, hold=1200000, part="b")]


def cell(c):
    c = dict(c); c["L0"] = c["N"] * math.pi / (8.0 * c["H"] * c["eta"]); c["name"] = f"N{c['N']}_eta{c['eta']:.3f}"; return c


def seed(c, l, r): return T.run_seed(BASE, l, c["k"], r)


def cmd(binp, d, c, start, s, r):
    out = []
    for a in P1.cmd(binp, d, c, s):
        if a.startswith("--wall-hold-steps="): a = f"--wall-hold-steps={c['hold']}"
        elif a.startswith("--seed="): a = f"--seed={BASE}"
        elif a.startswith("--gen3-seeding="): a = f"--gen3-seeding={start}"
        out.append(a)
    return out + ["--gen3-virial-blocks=50"] + ([f"--gen3-snapshots={SNAP_EVERY}"] if r in SNAP_SEEDS else [])


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--bin", required=True); ap.add_argument("--expect-bin-sha", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--jobs", type=int, default=8); ap.add_argument("--cap-file")
    a = ap.parse_args()
    if a.jobs > 8: sys.exit("STOP: at most 8 simulation processes (decision 15)")
    binp, O = os.path.abspath(a.bin), os.path.abspath(a.out); os.makedirs(O, exist_ok=True)
    if hashlib.sha256(open(binp, "rb").read()).hexdigest() != a.expect_bin_sha: sys.exit("STOP: binary SHA-256")
    build = subprocess.run([binp, "--version"], capture_output=True, text=True).stdout.splitlines()[0]
    log = open(os.path.join(O, "runner.log"), "a")
    def say(s): print(s, flush=True); log.write(s + "\n"); log.flush()
    def cap():
        try: return max(0, min(a.jobs, 8, int(open(a.cap_file).read().strip()))) if a.cap_file else a.jobs
        except (OSError, ValueError): return a.jobs
    say(f"# {time.strftime('%Y-%m-%d %H:%M:%S %Z')} P4: binary {binp} ({a.expect_bin_sha}; {build}); out {O}; jobs {a.jobs}; cap file {a.cap_file}")
    say(subprocess.run(["df", "-h", O], capture_output=True, text=True).stdout.strip())
    say(subprocess.run(["sysctl", "vm.swapusage"], capture_output=True, text=True).stdout.strip())
    jobs = [(c, lab, start, r, seed(c, l, r), os.path.join(O, c["name"], lab, f"seed{r}"))
            for c in map(cell, cells()) for l, (lab, start) in enumerate(STARTS) for r in range(NSEED)]
    for j in jobs:
        if os.path.exists(j[5]): sys.exit(f"STOP: {j[5]} exists -- a run is never repeated into the same folder")
    run, t0, stopped = [], time.time(), None
    def reap():
        for j in list(run):
            p, c, lab, r, d, ts = j
            if p.poll() is None: continue
            run.remove(j); open(os.path.join(d, "exit_code.txt"), "w").write(f"{p.returncode}\n")
            lg = open(os.path.join(d, "run.log"), errors="ignore").read()
            recs = re.findall(r"^\[EDMD3-HEALTH\] .*?: clean=(\S+) ", lg, re.M); nb = len(re.findall(r"^\[EDMD3\] built #", lg, re.M))
            bad = "" if (recs == ["1"] and nb == 1 and p.returncode == 0) else f"exit {p.returncode}, run record {recs}, builds {nb}"
            say(f"  {c['name']} {lab} seed{r}: exit {p.returncode}, {time.time() - ts:.0f} s" + (f"  **FINDING: {bad}**" if bad else ""))
    for c, lab, start, r, s, d in jobs:
        while len(run) >= cap(): reap(); time.sleep(1.0)
        free = shutil.disk_usage(O).free
        if free < FREE_MIN:
            stopped = f"free space {free / 1024 ** 3:.2f} GiB < 4 GiB before {c['name']} {lab} seed{r}"
            say(f"**STOPPED LAUNCHING: {stopped}; the runs already going finish; nothing else is launched**"); break
        os.makedirs(d)
        cm = cmd(binp, d, c, start, s, r)
        with open(os.path.join(d, "command.txt"), "w") as fh: fh.write(" ".join(cm) + "\n")
        run.append((subprocess.Popen(cm, cwd=d, stdout=open(os.path.join(d, "run.log"), "w"), stderr=subprocess.STDOUT), c, lab, r, d, time.time()))
        say(f"  start {c['name']} {lab} seed{r} (running {len(run)}, cap {cap()}, free {free / 1024 ** 3:.1f} GiB)")
    while run: reap(); time.sleep(1.0)
    say(f"# done {time.strftime('%Y-%m-%d %H:%M:%S %Z')} after {time.time() - t0:.0f} s" + (f"; STOPPED LAUNCHING ({stopped})" if stopped else ""))
    say(subprocess.run(["df", "-h", O], capture_output=True, text=True).stdout.strip())
    say(subprocess.run(["sysctl", "vm.swapusage"], capture_output=True, text=True).stdout.strip())


if __name__ == "__main__":
    main()
