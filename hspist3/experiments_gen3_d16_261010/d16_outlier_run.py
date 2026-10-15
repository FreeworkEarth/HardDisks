#!/usr/bin/env python3
"""##CHRIS 2026-10-10 (261012 sec. 4.7.37; plan-author decision 16, item 2a): THE PART-3A OUTLIER, ONE RUN -- information, no verdict.
Part 3a's held cell pi/8, seed 1 (experiments_gen3_p3a_261009, e4/pi8_seed1: Z_L - Z_R = +0.0074 over 5000 sigma-time) once more with
the same command and seed, the hold extended from 5000 to 2e4 sigma-time (--wall-hold-steps 300000 -> 1200000), and the initial state
written by --gen3-dump-initial (an output only: it writes the positions and velocities after the seeding; d16_outlier_tables.py prints
whether the first 100 virial blocks are identical to the original run's).
The command is built as p3a_evidence.py built it (P1's cmd, patched by stageH_evidence.patch, plus --gen3-virial-blocks=50), and this
script checks it against the original's command.txt: only the hold, the run folder and the added dump flag may differ.
BINARY: part 3a's frozen committed-tree build of engine-gen3 5d19aa7 (experiments_gen3_p3a_261009/bin_3a; SHA-256 checked).
One process; the run in its own folder <out>/pi8_seed1_hold2e4/; nothing is compressed, moved or deleted; a run is never repeated
into the same folder.
usage (from hspist3/ of the engine-gen3 worktree):
  python3 experiments_gen3_d16_261010/d16_outlier_run.py --bin <00ALLINONE> --expect-bin-sha <sha256> --out <data root>
"""
import argparse, hashlib, math, os, shutil, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(HS, "experiments_gen3_p1_261009")); sys.path.insert(0, os.path.join(HS, "experiments_gen3_h_261009"))
import p1_run as P1
import stageH_evidence as SH
T = P1.T
BASE, K, R, HOLD = 20261221, 0, 1, 1200000                          # part 3a's base, cell pi8 (k = 0), seed r = 1; 2e4 sigma-time
ORIG = os.path.join(HS, "experiments_gen3_p3a_261009", "evidence_5d19aa7", "e4", "pi8_seed1")
NAME, FREE_MIN = "pi8_seed1_hold2e4", 4 * 1024 ** 3


def cmd(binp, d):
    c = dict(N=400, H=20.0, L0=20.0, eta=math.pi / 8)
    return (SH.patch(P1.cmd(binp, d, c, T.run_seed(BASE, 0, K, R)), {"--wall-hold-steps=": str(HOLD), "--seed=": str(BASE)})
            + ["--gen3-virial-blocks=50", f"--gen3-dump-initial={os.path.join(d, 'initial_state.txt')}"])


def check_against_original(cm, d):
    """the new command against the original's command.txt: the differences, as (original, new) pairs"""
    old = open(os.path.join(ORIG, "command.txt")).read().split()
    ko = {a.split("=", 1)[0]: a for a in old[1:]}; kn = {a.split("=", 1)[0]: a for a in cm[1:]}
    diff = [(ko.get(k), kn.get(k)) for k in sorted(set(ko) | set(kn)) if ko.get(k) != kn.get(k)]
    allowed = {"--wall-hold-steps", "--speed-sound-run-dir", "--gen3-dump-initial"}
    bad = [p for p in diff if (p[0] or p[1]).split("=", 1)[0] not in allowed]
    order = [a.split("=", 1)[0] for a in old[1:]] == [a.split("=", 1)[0] for a in cm[1:] if not a.startswith("--gen3-dump-initial=")]
    return diff, bad, order, old[0] == cm[0]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--bin", required=True); ap.add_argument("--expect-bin-sha", required=True)
    ap.add_argument("--out", required=True); a = ap.parse_args()
    binp, O = os.path.abspath(a.bin), os.path.abspath(a.out); os.makedirs(O, exist_ok=True)
    if hashlib.sha256(open(binp, "rb").read()).hexdigest() != a.expect_bin_sha: sys.exit("STOP: binary SHA-256")
    d = os.path.join(O, NAME)
    if os.path.exists(d): sys.exit(f"STOP: {d} exists -- a run is never repeated into the same folder")
    cm = cmd(binp, d)
    diff, bad, order, same_bin = check_against_original(cm, d)
    print(f"# {time.strftime('%Y-%m-%d %H:%M:%S %Z')} item 2a: seed {T.run_seed(BASE, 0, K, R)}; binary {binp} ({a.expect_bin_sha})")
    print(f"against the original's command.txt: same binary path {same_bin}; same arguments in the same order (the dump flag aside) {order}")
    for o, n in diff: print(f"  original: {o}\n  this run: {n}")
    if bad or not order or not same_bin: sys.exit(f"STOP: the command differs from the original beyond the hold, the folder and the dump: {bad}")
    free = shutil.disk_usage(O).free
    if free < FREE_MIN: sys.exit(f"STOP: free space {free / 1024 ** 3:.2f} GiB < 4 GiB")
    print(subprocess.run(["df", "-h", O], capture_output=True, text=True).stdout.strip())
    os.makedirs(d)
    with open(os.path.join(d, "command.txt"), "w") as fh: fh.write(" ".join(cm) + "\n")
    t0 = time.time()
    rc = subprocess.run(cm, cwd=d, stdout=open(os.path.join(d, "run.log"), "w"), stderr=subprocess.STDOUT).returncode
    open(os.path.join(d, "exit_code.txt"), "w").write(f"{rc}\n")
    print(f"# done {time.strftime('%Y-%m-%d %H:%M:%S %Z')}: exit {rc}, {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
