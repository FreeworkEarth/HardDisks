#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.31; plan-author decision 13, item 5): PILOT P3, TWO STARTS -- a pilot, no paper use.
Starts: A = the M4 lattice at the target density (--gen3-seeding=lattice); B = M4b, the compressed fluid (--gen3-seeding=compressed with
its defaults: eta_0 0.60, T_melt 500 sigma-time, u 0.01 sigma per sigma-time). Cells: N = 400 (H = 20) at eta 0.700, 0.708, 0.716,
0.740; N = 1600 (H = 40) at 0.700, 0.708, 0.716; L0 = N pi / (8 H eta) exactly (--gen3-exact-box). After the start the divider is held
2e4 sigma-time (1,200,000 steps), --psi6-every=1, --gen3-virial-blocks=50 (400 blocks); then a 1-sigma-time record (the
speed-of-sound loop needs one). P1's command otherwise (M = 300). 6 seeds per start and cell.
SEEDS: a new base, BASE = 20261224 (in no file before this script); run_seed(BASE, l, k, r) with l = 0 for start A and l = 1 for start
B (independent seeds for the two starts), k the cell index, r = 0..5.
BINARY: the frozen committed-tree build of engine-gen3 06e7248 (the M4b commit; experiments_gen3_m4b_261009/bin_m4b, SHA-256 checked).
At most 12 processes (programme rule 1), N = 1600 first. Every run in its own folder <out>/<cell>/<start>/seed<r>/; nothing is compressed,
moved or deleted; a run with exit != 0, clean != 1 or no run record is a finding, reported and never rerun.
usage (from hspist3/ of the engine-gen3 worktree):
  python3 experiments_gen3_p3_261009/p3_run.py --bin <00ALLINONE> --expect-bin-sha <sha256> --out <data root> [--jobs 12]
"""
import argparse, hashlib, math, os, re, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(HS, "experiments_gen3_p1_261009"))
import p1_run as P1
T = P1.T
BASE, NSEED, HOLD = 20261224, 6, 1200000
STARTS = (("A", "lattice"), ("B", "compressed"))


def cells():
    out = []
    for N, H, etas in ((1600, 40.0, (0.700, 0.708, 0.716)), (400, 20.0, (0.700, 0.708, 0.716, 0.740))):
        for eta in etas:
            out.append(dict(k=len(out), N=N, H=H, eta=eta, L0=N * math.pi / (8.0 * H * eta), name=f"N{N}_eta{eta:.3f}"))
    return out


def seed(c, l, r): return T.run_seed(BASE, l, c["k"], r)


def cmd(binp, d, c, start, s):
    out = []
    for a in P1.cmd(binp, d, c, s):
        if a.startswith("--wall-hold-steps="): a = f"--wall-hold-steps={HOLD}"
        elif a.startswith("--seed="): a = f"--seed={BASE}"
        elif a.startswith("--gen3-seeding="): a = f"--gen3-seeding={start}"
        out.append(a)
    return out + ["--gen3-virial-blocks=50"]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--bin", required=True); ap.add_argument("--expect-bin-sha", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--jobs", type=int, default=12)
    a = ap.parse_args()
    if a.jobs > 12: sys.exit("STOP: at most 12 processes (programme rule 1)")
    binp, O = os.path.abspath(a.bin), os.path.abspath(a.out); os.makedirs(O, exist_ok=True)
    if hashlib.sha256(open(binp, "rb").read()).hexdigest() != a.expect_bin_sha: sys.exit("STOP: binary SHA-256")
    build = subprocess.run([binp, "--version"], capture_output=True, text=True).stdout.splitlines()[0]
    log = open(os.path.join(O, "runner.log"), "a")
    def say(s): print(s, flush=True); log.write(s + "\n"); log.flush()
    say(f"# {time.strftime('%Y-%m-%d %H:%M:%S %Z')} P3: binary {binp} ({a.expect_bin_sha}; {build}); out {O}; jobs {a.jobs}")
    say(subprocess.run(["df", "-h", O], capture_output=True, text=True).stdout.strip())
    jobs = [(c, lab, start, r, seed(c, l, r), os.path.join(O, c["name"], lab, f"seed{r}"))
            for c in cells() for l, (lab, start) in enumerate(STARTS) for r in range(NSEED)]
    for j in jobs:
        if os.path.exists(j[5]): sys.exit(f"STOP: {j[5]} exists -- a run is never repeated into the same folder")
    run, t0 = [], time.time()
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
        while len(run) >= a.jobs: reap(); time.sleep(0.5)
        os.makedirs(d)
        cm = cmd(binp, d, c, start, s)
        with open(os.path.join(d, "command.txt"), "w") as fh: fh.write(" ".join(cm) + "\n")
        run.append((subprocess.Popen(cm, cwd=d, stdout=open(os.path.join(d, "run.log"), "w"), stderr=subprocess.STDOUT), c, lab, r, d, time.time()))
    while run: reap(); time.sleep(0.5)
    say(f"# done {time.strftime('%Y-%m-%d %H:%M:%S %Z')} after {time.time() - t0:.0f} s")
    say(subprocess.run(["df", "-h", O], capture_output=True, text=True).stdout.strip())


if __name__ == "__main__":
    main()
