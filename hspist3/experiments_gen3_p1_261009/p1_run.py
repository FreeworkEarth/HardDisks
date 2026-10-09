#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.17; stage D of the plan-author programme of sec. 4.7.12): PILOT P1, the structural clock. A
pilot: it sizes the campaign; it is not a result.
gen3, M4 seeding (--gen3-seeding=lattice, the default jitter 0.25), the divider held throughout. Geometry: the gate's dense
construction H = 10 sqrt(N/100), L0 = N pi / (8 H eta) exactly (--gen3-exact-box). N = 100, 400, 1600; eta = 0.60, 0.66, 0.68, 0.70,
0.704, 0.708, 0.712, 0.716, 0.72, 0.74, 0.78; 6 seeds per cell, run_seed(20261217, 0, cell index, r) (base never used before).
Hold 2e4 sigma-time (1,200,000 steps of 1/60 sigma-time), then 1 sigma-time released (the speed-of-sound loop needs a record);
--psi6-every=1. Each run in its own folder; at most 12 processes; order: N = 1600 first (the longest), then 400, then 100, each cell's
seeds together. The psi6(t) file is compressed losslessly after its run (gzip -9; SHA-256 of the uncompressed file in the cell's
.sha256_uncompressed first). A gen3 run with clean=0 or without its run record is a finding (programme rule 10): it is reported in
the runner log and kept, never rerun.
usage (from hspist3/ of the engine-gen3 worktree):
  python3 experiments_gen3_p1_261009/p1_run.py --bin <frozen 00ALLINONE> --expect-bin-sha <sha256> --out <data root> [--jobs 12]
"""
import argparse, hashlib, math, os, re, subprocess, sys, time

MAIN_HS = "/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3"
sys.path.insert(0, os.path.join(MAIN_HS, "validation"))
import tests_20260913 as T
NS = (1600, 400, 100)
ETAS = (0.60, 0.66, 0.68, 0.70, 0.704, 0.708, 0.712, 0.716, 0.72, 0.74, 0.78)
NSEED, BASE, HOLD_STEPS = 6, 20261217, 1200000


def cells():
    out = []
    for N in NS:
        for eta in ETAS:
            k = len(out); H = 10.0 * math.sqrt(N / 100.0)
            out.append(dict(k=k, N=N, eta=eta, H=H, L0=N * math.pi / (8.0 * H * eta),
                            seeds=[T.run_seed(BASE, 0, k, r) for r in range(NSEED)], name=f"N{N}_eta{eta:.3f}"))
    return out


def cmd(binp, d, c, seed):
    ns = c["N"] // 2
    return [binp, "--mode=edmd", "--experiment=speed_of_sound", "--headless", "--kbt1", "--seed-drift-order=drift-first", "--edmd-acc=0",
            f"--particles={c['N']}", f"--particles-boxes={ns},{ns}", f"--height={c['H']:.17g}", "--particle-radius=0.5", "--wall-thickness=0.05",
            "--wall-thickness-vis=0.05", f"--lengths={c['L0']:.17g}", "--wall-masses=300", "--repeats=1", f"--seed={BASE}",
            f"--wall-hold-steps={HOLD_STEPS}", "--fixed-dt=0.4", "--record-sigma-time=1", "--oscillation-min-steps=10", "--speed-sound-log-stride=60",
            f"--speed-sound-run-dir={d}", f"--speed-sound-exact-seed={seed}", "--engine=gen3", "--gen3-seeding=lattice", "--gen3-exact-box",
            "--psi6-every=1"]


def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--bin", required=True); ap.add_argument("--expect-bin-sha", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--jobs", type=int, default=12)
    a = ap.parse_args()
    if a.jobs > 12: sys.exit("STOP: at most 12 processes (programme rule 1)")
    binp, O = os.path.abspath(a.bin), os.path.abspath(a.out); os.makedirs(O, exist_ok=True)
    if sha(binp) != a.expect_bin_sha: sys.exit("STOP: binary SHA-256")
    build = subprocess.run([binp, "--version"], capture_output=True, text=True).stdout.splitlines()[0]
    log = open(os.path.join(O, "runner.log"), "a")
    def say(s): print(s, flush=True); log.write(s + "\n"); log.flush()
    say(f"# {time.strftime('%Y-%m-%d %H:%M:%S %Z')} P1: binary {binp} ({a.expect_bin_sha}; {build}); out {O}; jobs {a.jobs}")
    say(subprocess.run(["df", "-h", O], capture_output=True, text=True).stdout.strip())
    jobs = []
    for c in cells():
        for r, s in enumerate(c["seeds"]):
            jobs.append((c, r, s, os.path.join(O, c["name"], f"seed{r}")))
    run, t0 = [], time.time()
    def finish(c, r, d):
        lg = open(os.path.join(d, "run.log"), errors="ignore").read()
        recs = re.findall(r"^\[EDMD3-HEALTH\] .*?: clean=(\S+) ", lg, re.M); nb = len(re.findall(r"^\[EDMD3\] built #", lg, re.M))
        bad = "INFEASIBLE" in lg and "INFEASIBLE" or ("" if recs == ["1"] and nb == 1 else f"run record {recs}, builds {nb}")
        ps = [f for f in os.listdir(d) if f.startswith("psi6_t_") and f.endswith(".csv")]
        for f in ps:
            with open(os.path.join(os.path.dirname(d), ".sha256_uncompressed"), "a") as fh: fh.write(f"{sha(os.path.join(d, f))}  seed{r}/{f}\n")
            subprocess.run(["gzip", "-9", os.path.join(d, f)], check=True)
        return bad
    def reap():
        for j in list(run):
            p, c, r, d, ts = j
            if p.poll() is not None:
                run.remove(j); bad = finish(c, r, d)
                say(f"  {c['name']} seed{r}: exit {p.returncode}, {time.time() - ts:.0f} s" + (f"  **FINDING: {bad}**" if bad or p.returncode else ""))
    for c, r, s, d in jobs:
        while len(run) >= a.jobs: reap(); time.sleep(0.5)
        os.makedirs(d, exist_ok=True)
        cm = cmd(binp, d, c, s)
        with open(os.path.join(d, "command.txt"), "w") as fh: fh.write(" ".join(cm) + "\n")
        run.append((subprocess.Popen(cm, cwd=d, stdout=open(os.path.join(d, "run.log"), "w"), stderr=subprocess.STDOUT), c, r, d, time.time()))
    while run: reap(); time.sleep(0.5)
    say(f"# done {time.strftime('%Y-%m-%d %H:%M:%S %Z')} after {time.time() - t0:.0f} s")
    say(subprocess.run(["df", "-h", O], capture_output=True, text=True).stdout.strip())


if __name__ == "__main__":
    main()
