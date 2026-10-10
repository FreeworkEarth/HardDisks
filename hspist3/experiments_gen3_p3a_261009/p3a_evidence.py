#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.26; plan-author decision 12, part 3a): the evidence of the virial pressure per compartment and
the position snapshots, with the binaries of a committed tree (build_clean.sh plus gen3_checkpoint_test from the same archive).
Acceptance (decision 12): observation does not steer (same event hash on/off); rule 4 identical; in a symmetric held cell the two
compartments' Z agree within their SEs and their weighted mean equals the global Z. With them, because part 3a changes the engine
and the checkpoint: stage H's restart tests again (and part 3b's: two checkpoints of the same state byte-identical).
  h    stage H's evidence runner, unchanged (experiments_gen3_h_261009/stageH_evidence.py): rule 4, the M1/M2 harness outputs, the
       engine checkpoint test (both builds), the driver's restarts on three trajectories, the refusals
  e3   observation does not steer: P1's cell N = 400, eta 0.704 (its seed 0), M = 50, held 2000 sigma-time, a 100-sigma-time record:
       plain; --gen3-virial-blocks=50; --gen3-snapshots=10; both
  e4   the symmetric held cell: N = 400, H = 20 at eta pi/8 (L0 = 20) and eta 0.70 (P1's L0), three seeds each (base 20261221, never
       used before), the divider held 5000 sigma-time, --gen3-virial-blocks=50 (100 blocks), a 1-sigma-time record
  e5   the snapshots' content: P1's cell N = 400, eta 0.704, --gen3-snapshots=10 with --gen3-dump-initial, held 100 sigma-time
  e6   a restart with both observations on: e3's command with both flags; U, and C / R1 + R2 at hold:60000 and record:3000
At most 12 simulation processes at once (stage H's runner runs alone first; then at most 10). Nothing is deleted.
usage (from hspist3/ of the engine-gen3 worktree):
  python3 experiments_gen3_p3a_261009/p3a_evidence.py --bin-dir <clean build dir> --out <evidence dir>
"""
import argparse, math, os, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(HS, "experiments_gen3_p1_261009")); sys.path.insert(0, os.path.join(HS, "experiments_gen3_h_261009"))
import p1_run as P1
import stageH_evidence as SH                       # patch(), Runner (unchanged)
T = P1.T
JOBS = 10
BASE = 20261221


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--bin-dir", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args(); B = os.path.abspath(a.bin_dir); os.makedirs(a.out, exist_ok=True); O = os.path.abspath(a.out)
    binp = os.path.join(B, "00ALLINONE")
    R = SH.Runner(O)
    R.say(f"# part 3a evidence {time.strftime('%Y-%m-%d %H:%M:%S %Z')}: binaries {B}")
    R.say(subprocess.run(["df", "-h", O], capture_output=True, text=True).stdout.strip())
    # ---- h: stage H's runner, alone (it runs up to 12 simulations itself)
    if not os.path.exists(os.path.join(O, "h", "runner.log")) or "# done" not in open(os.path.join(O, "h", "runner.log")).read():
        rc = subprocess.run([sys.executable, os.path.join(HS, "experiments_gen3_h_261009", "stageH_evidence.py"), "--bin-dir", B,
                             "--out", os.path.join(O, "h")], stdout=open(os.path.join(O, "h_runner.txt"), "w"), stderr=subprocess.STDOUT).returncode
        R.say(f"  stage H runner: exit {rc}")
    # ---- the commands
    p1 = {(c["N"], round(c["eta"], 3)): c for c in P1.cells()}
    c704 = p1[(400, 0.704)]
    def cmd704(d, extra=()):
        return SH.patch(P1.cmd(binp, d, c704, c704["seeds"][0]), {"--wall-masses=": "50", "--wall-hold-steps=": "120000",
                        "--record-sigma-time=": "100", "--speed-sound-log-stride=": "6"}) + list(extra)
    VB, SN = "--gen3-virial-blocks=50", "--gen3-snapshots=10"
    runs = []
    # e3
    for tag, ex in (("plain", ()), ("vb", (VB,)), ("snap", (SN,)), ("both", (VB, SN))):
        runs.append((f"e3 {tag}", cmd704(os.path.join(O, "e3", tag), ex), os.path.join(O, "e3", tag)))
    # e4
    cells4 = [("pi8", dict(N=400, H=20.0, L0=20.0, eta=math.pi / 8)), ("070", p1[(400, 0.70)])]
    for k, (name, c) in enumerate(cells4):
        for r in range(3):
            seed = T.run_seed(BASE, 0, k, r)
            d = os.path.join(O, "e4", f"{name}_seed{r}")
            cm = SH.patch(P1.cmd(binp, d, c, seed), {"--wall-hold-steps=": "300000", "--seed=": str(BASE)}) + [VB]
            runs.append((f"e4 {name} seed{r}", cm, d))
    # e5
    d5 = os.path.join(O, "e5")
    runs.append(("e5 snapshots vs the initial dump", SH.patch(P1.cmd(binp, d5, c704, c704["seeds"][0]), {"--wall-hold-steps=": "6000"})
                 + [SN, f"--gen3-dump-initial={os.path.join(d5, 'initial_state.txt')}"], d5))
    # e6: U and the writers (C and R1); R2 after its R1
    base6 = os.path.join(O, "e6")
    runs.append(("e6 U", cmd704(os.path.join(base6, "U"), (VB, SN)), os.path.join(base6, "U")))
    pend = []
    for ck in ("hold:60000", "record:3000"):
        t = ck.replace(":", "")
        runs.append((f"e6 C {ck}", cmd704(os.path.join(base6, f"C_{t}"), (VB, SN, f"--gen3-checkpoint={ck}:{base6}/ck_C_{t}.bin")),
                     os.path.join(base6, f"C_{t}")))
        pend.append((t, f"--gen3-checkpoint={ck}:{base6}/ck_R1_{t}.bin"))
    procs = {}
    for name, cm, d in runs:
        procs[name] = R.start(name, cm, d)
    r1 = {}
    for t, ckarg in pend:
        d = os.path.join(base6, f"R1_{t}")
        r1[t] = R.start(f"e6 R1 {t}", cmd704(d, (VB, SN, ckarg, "--gen3-checkpoint-stop")), d)
    for t, _ in pend:
        R.wait(r1[t])
        d = os.path.join(base6, f"R2_{t}")
        R.start(f"e6 R2 {t}", cmd704(d, (VB, SN, f"--gen3-restart={base6}/ck_R1_{t}.bin")), d)
    R.wait_all()
    R.say(f"# done {time.strftime('%Y-%m-%d %H:%M:%S %Z')}")
    R.say(subprocess.run(["df", "-h", O], capture_output=True, text=True).stdout.strip())


if __name__ == "__main__":
    main()
