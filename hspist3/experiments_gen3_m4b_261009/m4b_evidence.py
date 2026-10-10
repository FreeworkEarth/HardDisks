#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.30; plan-author decision 13, item 4): the evidence of M4b, the compressed-fluid start
(--gen3-seeding=compressed), with the binaries of a committed tree (build_clean.sh <commit>, plus gen3_checkpoint_test from the same
archive). Criteria: m4b_tables.py (committed with this script, before the run).
  r4  rule 4: audit_runs_261007.py on ctrl_min and ctrl_leg, the default build and --engine=gen2 (the tables copy the 7b08827 reference
      outputs in, as stage H's tables do)
  h   the M1/M2 harness outputs (gen3_m1 audit, gen3_m2 audit --quick, gen3_m2 audit)
  c   the four cells of the acceptance: N = 400 (H = 20) and N = 1600 (H = 40) at the target eta 0.70 and 0.716, L0 = N pi / (8 H eta)
      exactly; M4b with its defaults (eta_0 0.60, T_melt 500 sigma-time, u 0.01 sigma per sigma-time); P1's command otherwise (M = 300,
      a 1-sigma-time record) with a 200-sigma-time hold; seed run_seed(20261223, 0, k, 0) (base never used before); HD_CONTACT_AUDIT=1.
      Per cell: "plain" (psi6 every 1 sigma-time), "repeat" (the same command again), "obs" (psi6 every 0.25, --gen3-virial-blocks=50,
      --gen3-snapshots=10)
  k   checkpoint and restart through the compression, on N = 400 at 0.716: C (--gen3-checkpoint=hold:6000, going on), R1 (the same
      checkpoint, then --gen3-checkpoint-stop), R2 (restarted from R1's checkpoint: the preparation runs again, then the checkpoint);
      U is the cell's "plain" run
  x   refusals (each must stop with exit 2): compressed under gen2; a compression flag without compressed; compressed without
      --gen3-exact-box; eta_0 above the target eta
At most 12 simulation processes (programme rule 1: this runner's children, the two audit runners running 2 each). Nothing is deleted.
usage (from hspist3/ of the engine-gen3 worktree):
  python3 experiments_gen3_m4b_261009/m4b_evidence.py --bin-dir <clean build dir> --out <evidence dir>
"""
import argparse, math, os, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(HS, "experiments_gen3_p1_261009")); sys.path.insert(0, os.path.join(HS, "experiments_gen3_h_261009"))
sys.path.insert(0, os.path.join(HS, "experiments_gen3_p2_261009"))
import p1_run as P1
import stageH_evidence as SH                       # patch(), Runner (JOBS = 10 children)
T = P1.T
BASE = 20261223
GATE = SH.GATE


def cells():
    out = []
    for N, H in ((400, 20.0), (1600, 40.0)):
        for eta in (0.70, 0.716):
            k = len(out)
            out.append(dict(k=k, N=N, H=H, eta=eta, L0=N * math.pi / (8.0 * H * eta), name=f"N{N}_eta{eta:.3f}", seed=T.run_seed(BASE, 0, k, 0)))
    return out


def cmd(binp, d, c, extra=()):
    return SH.patch(P1.cmd(binp, d, c, c["seed"]), {"--gen3-seeding=": "compressed", "--wall-hold-steps=": "12000", "--seed=": str(BASE)}) + list(extra)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--bin-dir", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args(); B = os.path.abspath(a.bin_dir); os.makedirs(a.out, exist_ok=True); O = os.path.abspath(a.out)
    binp = os.path.join(B, "00ALLINONE")
    R = SH.Runner(O)
    os.environ["HD_CONTACT_AUDIT"] = "1"
    R.say(f"# M4b evidence {time.strftime('%Y-%m-%d %H:%M:%S %Z')}: binaries {B}")
    R.say(subprocess.run(["df", "-h", O], capture_output=True, text=True).stdout.strip())
    # r4 and h
    with open(os.path.join(B, "bin_engine_gen2.sh"), "w") as fh: fh.write(f'#!/bin/bash\nexec "{binp}" "$@" --engine=gen2\n')
    os.chmod(os.path.join(B, "bin_engine_gen2.sh"), 0o755)
    for tag, b in (("default", binp), ("gen2flag", os.path.join(B, "bin_engine_gen2.sh"))):
        R.start(f"r4 {tag}", [sys.executable, os.path.join(GATE, "audit_runs_261007.py"), "run", "--bin", b, "--out", os.path.join(O, "r4", tag),
                             "--cases", "ctrl_min,ctrl_leg", "--jobs", "2"], os.path.join(O, "r4", f"{tag}_runner"))
    R.start("h", ["/bin/bash", "-c", f'"{B}/gen3_m1" audit > m1_audit_output.txt && "{B}/gen3_m2" audit --quick > m2_audit_quick_output.txt'
                  f' && "{B}/gen3_m2" audit > m2_audit_output.txt'], os.path.join(O, "h"))
    # c: the cells (N = 1600 first, the longer runs)
    C = cells()
    for c in sorted(C, key=lambda c: -c["N"]):
        for tag, extra in (("plain", ()), ("repeat", ()), ("obs", ("--gen3-virial-blocks=50", "--gen3-snapshots=10"))):
            d = os.path.join(O, "c", c["name"], tag)
            cm = cmd(binp, d, c, extra)
            if tag == "obs": cm = SH.patch(cm, {"--psi6-every=": "0.25"})
            R.start(f"c {c['name']} {tag}", cm, d)
    # k: checkpoint and restart through the compression (N = 400, 0.716)
    ck = [c for c in C if c["N"] == 400 and abs(c["eta"] - 0.716) < 1e-9][0]; kb = os.path.join(O, "k")
    R.start("k C", cmd(binp, os.path.join(kb, "C"), ck, (f"--gen3-checkpoint=hold:6000:{kb}/ck_C.bin",)), os.path.join(kb, "C"))
    p = R.start("k R1", cmd(binp, os.path.join(kb, "R1"), ck, (f"--gen3-checkpoint=hold:6000:{kb}/ck_R1.bin", "--gen3-checkpoint-stop")), os.path.join(kb, "R1"))
    R.wait(p)
    R.start("k R2", cmd(binp, os.path.join(kb, "R2"), ck, (f"--gen3-restart={kb}/ck_R1.bin",)), os.path.join(kb, "R2"))
    # x: refusals
    c0 = C[0]; xb = os.path.join(O, "x")
    gen2 = [x for x in cmd(binp, os.path.join(xb, "gen2"), c0) if not x.startswith(("--engine=", "--gen3-", "--psi6-every"))] + ["--engine=gen2", "--gen3-seeding=compressed"]
    for k, (name, cm) in enumerate((
            ("compressed under gen2", gen2),
            ("a compression flag without compressed", SH.patch(cmd(binp, os.path.join(xb, "noflag"), c0), {"--gen3-seeding=": "lattice"}) + ["--gen3-compress-melt=10"]),
            ("compressed without --gen3-exact-box", SH.patch(cmd(binp, os.path.join(xb, "nobox"), c0), {"--gen3-exact-box=": None})),
            ("eta_0 above the target eta", cmd(binp, os.path.join(xb, "eta0"), c0, ("--gen3-compress-eta0=0.75",))))):
        d = os.path.join(xb, f"case{k}"); os.makedirs(d, exist_ok=True); open(os.path.join(xb, f"case{k}.name"), "w").write(name + "\n")
        R.start(f"x {name}", cm, d)
    R.wait_all()
    R.say(f"# done {time.strftime('%Y-%m-%d %H:%M:%S %Z')}")
    R.say(subprocess.run(["df", "-h", O], capture_output=True, text=True).stdout.strip())


if __name__ == "__main__":
    main()
