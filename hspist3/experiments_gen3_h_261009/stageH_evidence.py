#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.21; stage H of the plan-author programme of sec. 4.7.12): the evidence of checkpoint and
restart at an event boundary, with the binaries of a committed tree (build_clean.sh <commit>, plus gen3_checkpoint_test built from
the same archive, default and -DEDMD3_LONG_DOUBLE). Acceptance (sec. 4.7.12): a restarted run is byte-identical to the uninterrupted
one (traces, psi6(t), event hash).
  H0  rule 4: gen2 byte identity with 7b08827 (ctrl_min, ctrl_leg; the default build and --engine=gen2): 00ALLINONE.c changed
  H1  the M1 and M2 harness outputs against the committed ones (the engine gained two functions and nothing else)
  H2  the engine test edmd_core/tests/gen3_checkpoint_test.c (both builds)
  H3  the driver, three trajectories of the day's own commands:
        A  pilot P2's cell eta 0.7113, M = 300, its seed 0, exactly as P2 ran it (N = 100, held 1e4 sigma-time = 600000 steps:
           29 origin shifts in the hold; then 200 periods);
        B  pilot P1's cell N = 400, eta 0.704, its seed 0, with M = 50, a 2000-sigma-time hold and a 200-sigma-time record;
        C  pilot P1's cell N = 1600, eta 0.716, its seed 0, with M = 2000, a 500-sigma-time hold and a 100-sigma-time record;
      U the uninterrupted run; at each checkpoint C (the same command writing the checkpoint and going on) and R (the same command
      stopped right after writing it, R1, then restarted from it, R2); a chain in A (a restarted run writes a later checkpoint and
      stops, R2b; a third process goes on from it, R3)
  H4  the driver's refusals (another seed, mass or initial state; not a checkpoint; another build; a truncated file; gen2; two
      trajectories in one process; a checkpoint step that is never reached; --gen3-checkpoint-stop alone)
At most 12 processes at once (programme rule 1). Every run in its own folder, its stdout in run.log; nothing is deleted.
usage (from hspist3/ of the engine-gen3 worktree):
  python3 experiments_gen3_h_261009/stageH_evidence.py --bin-dir <clean build dir> --out <evidence dir>
"""
import argparse, math, os, shutil, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(HS, "experiments_gen3_p2_261009")); sys.path.insert(0, os.path.join(HS, "experiments_gen3_p1_261009"))
import p2_run as P2
import p1_run as P1
GATE = "/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3/cluster/resched_gate_261005"
P1TAB = os.path.join(HS, "experiments_gen3_p1_261009", "p1_tables_output.txt")
JOBS = 12


def patch(cmd, repl):
    """replace the arguments whose '--name=' prefix is a key of repl (value None drops it)"""
    out = []
    for a in cmd:
        k = a.split("=", 1)[0] + "="
        if k in repl:
            if repl[k] is not None: out.append(k + repl[k])
        else: out.append(a)
    return out


def cells(binp):
    A = [c for c in P2.cells(P2.teq_table(P1TAB)) if c["name"] == "eta0.7113_M300"][0]
    p1 = {(c["N"], round(c["eta"], 3)): c for c in P1.cells()}
    B, C = p1[(400, 0.704)], p1[(1600, 0.716)]
    def cmdA(d): return P2.cmd(binp, d, A, A["seeds"][0])
    def cmdB(d): return patch(P1.cmd(binp, d, B, B["seeds"][0]), {"--wall-masses=": "50", "--wall-hold-steps=": "120000",
                                                                    "--record-sigma-time=": "200", "--speed-sound-log-stride=": "6"})
    def cmdC(d): return patch(P1.cmd(binp, d, C, C["seeds"][0]), {"--wall-masses=": "2000", "--wall-hold-steps=": "30000",
                                                                    "--record-sigma-time=": "100", "--speed-sound-log-stride=": "6"})
    return {"A": (cmdA, ["hold:1", "hold:300000", "record:0", "record:25000"]),
            "B": (cmdB, ["hold:60000", "record:6000"]),
            "C": (cmdC, ["hold:15000", "record:3000"])}


class Runner:
    def __init__(self, out): self.out = out; self.run = []; self.log = open(os.path.join(out, "runner.log"), "a")
    def say(self, s): print(s, flush=True); self.log.write(s + "\n"); self.log.flush()
    def start(self, name, cmd, d):
        while len([p for p in self.run if p[0].poll() is None]) >= JOBS: time.sleep(0.5)
        os.makedirs(d, exist_ok=True)
        open(os.path.join(d, "command.txt"), "w").write(" ".join(cmd) + "\n")
        p = subprocess.Popen(cmd, cwd=d, stdout=open(os.path.join(d, "run.log"), "w"), stderr=subprocess.STDOUT)
        self.run.append((p, name, d, time.time())); self.say(f"  start {name}")
        return p
    def wait(self, p):
        rc = p.wait()
        for q in self.run:
            if q[0] is p: self.say(f"  done {q[1]}: exit {rc}, {time.time() - q[3]:.0f} s"); open(os.path.join(q[2], "exit_code.txt"), "w").write(f"{rc}\n")
        return rc
    def wait_all(self):
        for q in list(self.run):
            if not os.path.exists(os.path.join(q[2], "exit_code.txt")): self.wait(q[0])


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--bin-dir", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args(); B = os.path.abspath(a.bin_dir); os.makedirs(a.out, exist_ok=True); O = os.path.abspath(a.out)
    binp = os.path.join(B, "00ALLINONE")
    R = Runner(O)
    R.say(f"# stage H evidence {time.strftime('%Y-%m-%d %H:%M:%S %Z')}: binaries {B}")
    R.say(subprocess.run(["df", "-h", O], capture_output=True, text=True).stdout.strip())
    # ---- H1 and H2 (4 processes) and H0 (2 at a time, 4 runs)
    os.makedirs(os.path.join(O, "h1"), exist_ok=True)
    h1 = R.start("H1", ["/bin/bash", "-c", f'"{B}/gen3_m1" audit > m1_audit_output.txt && "{B}/gen3_m2" audit --quick > m2_audit_quick_output.txt'
                        f' && "{B}/gen3_m2" audit > m2_audit_output.txt'], os.path.join(O, "h1"))
    h2a = R.start("H2 default", ["/bin/bash", "-c", f'"{B}/gen3_checkpoint_test" > checkpoint_test_output.txt'], os.path.join(O, "h2", "default"))
    h2b = R.start("H2 long double", ["/bin/bash", "-c", f'"{B}/gen3_checkpoint_test_ld" > checkpoint_test_output.txt'], os.path.join(O, "h2", "ld"))
    with open(os.path.join(B, "bin_engine_gen2.sh"), "w") as fh: fh.write(f'#!/bin/bash\nexec "{binp}" "$@" --engine=gen2\n')
    os.chmod(os.path.join(B, "bin_engine_gen2.sh"), 0o755)
    for tag, b in (("default", binp), ("gen2flag", os.path.join(B, "bin_engine_gen2.sh"))):
        R.start(f"H0 {tag}", [sys.executable, os.path.join(GATE, "audit_runs_261007.py"), "run", "--bin", b, "--out", os.path.join(O, "h0", tag),
                              "--cases", "ctrl_min,ctrl_leg", "--jobs", "2"], os.path.join(O, "h0", f"{tag}_runner"))
    # ---- H3: U, C and R1 of every checkpoint at once; R2 after its R1; the chain after A's R2 of hold:300000
    C = cells(binp)
    pend = []
    for cell, (cmd, cks) in C.items():
        base = os.path.join(O, "h3", cell)
        R.start(f"{cell} U", cmd(os.path.join(base, "U")), os.path.join(base, "U"))
        for ck in cks:
            t = ck.replace(":", "")
            R.start(f"{cell} C {ck}", cmd(os.path.join(base, f"C_{t}")) + [f"--gen3-checkpoint={ck}:{base}/ck_C_{t}.bin"], os.path.join(base, f"C_{t}"))
            p = R.start(f"{cell} R1 {ck}", cmd(os.path.join(base, f"R1_{t}")) + [f"--gen3-checkpoint={ck}:{base}/ck_R1_{t}.bin", "--gen3-checkpoint-stop"],
                        os.path.join(base, f"R1_{t}"))
            pend.append((cell, cmd, base, t, p))
    chain = None
    for cell, cmd, base, t, p in pend:
        R.wait(p)
        extra = [f"--gen3-restart={base}/ck_R1_{t}.bin"]
        q = R.start(f"{cell} R2 {t}", cmd(os.path.join(base, f"R2_{t}")) + extra, os.path.join(base, f"R2_{t}"))
        if cell == "A" and t == "hold300000":   # the chain: a restarted run writes a later checkpoint and stops; a third process goes on
            chain = R.start("A R2b (restarted at hold:300000, checkpoint at record:25000, stop)",
                            cmd(os.path.join(base, "R2b_chain")) + extra + [f"--gen3-checkpoint=record:25000:{base}/ck_R2b_chain.bin", "--gen3-checkpoint-stop"],
                            os.path.join(base, "R2b_chain"))
    if chain is not None:
        R.wait(chain)
        base = os.path.join(O, "h3", "A")
        R.start("A R3 (restarted from R2b's checkpoint)", C["A"][0](os.path.join(base, "R3_chain")) + [f"--gen3-restart={base}/ck_R2b_chain.bin"],
                os.path.join(base, "R3_chain"))
    # ---- H4: refusals (each must stop with exit 2 and its message)
    baseA = os.path.join(O, "h3", "A"); ck = os.path.join(baseA, "ck_R1_hold300000.bin")
    R.wait_all()
    hb = os.path.join(O, "h4"); os.makedirs(hb, exist_ok=True)
    raw = open(ck, "rb").read()
    bad_build = bytearray(raw); bad_build[16] = ord("X") if bad_build[16] != ord("X") else ord("Y")
    open(os.path.join(hb, "ck_other_build.bin"), "wb").write(bytes(bad_build))
    open(os.path.join(hb, "ck_truncated.bin"), "wb").write(raw[: len(raw) // 2])
    cmdA, cmdB, cmdCc = C["A"][0], C["B"][0], C["C"][0]
    trace = [f for f in os.listdir(os.path.join(baseA, "U")) if f.startswith("wall_x_positions_")][0]
    gen2 = [x for x in cmdA(os.path.join(hb, "gen2")) if not x.startswith(("--engine=", "--gen3-", "--psi6-every"))] + ["--engine=gen2", f"--gen3-restart={ck}"]
    refusals = [
        ("another seed", patch(cmdA(os.path.join(hb, "seed")), {"--speed-sound-exact-seed=": "12345"}) + [f"--gen3-restart={ck}"]),
        ("another mass", patch(cmdA(os.path.join(hb, "mass")), {"--wall-masses=": "2000"}) + [f"--gen3-restart={ck}"]),
        ("another initial state (jitter 0.2)", cmdA(os.path.join(hb, "jitter")) + ["--gen3-jitter=0.2", f"--gen3-restart={ck}"]),
        ("not a checkpoint (a trace)", cmdA(os.path.join(hb, "notckpt")) + [f"--gen3-restart={os.path.join(baseA, 'U', trace)}"]),
        ("another build (the build field changed)", cmdA(os.path.join(hb, "build")) + [f"--gen3-restart={os.path.join(hb, 'ck_other_build.bin')}"]),
        ("a truncated checkpoint", cmdA(os.path.join(hb, "trunc")) + [f"--gen3-restart={os.path.join(hb, 'ck_truncated.bin')}"]),
        ("--gen3-restart under gen2", gen2),
        ("two trajectories in one process", patch(cmdB(os.path.join(hb, "repeats")), {"--repeats=": "2"}) + [f"--gen3-checkpoint=hold:10:{hb}/never.bin"]),
        ("a checkpoint step never reached (hold:30000 of 30000)", cmdCc(os.path.join(hb, "never")) + [f"--gen3-checkpoint=hold:30000:{hb}/never.bin"]),
        ("--gen3-checkpoint-stop alone", cmdB(os.path.join(hb, "stopalone")) + ["--gen3-checkpoint-stop"]),
    ]
    for k, (name, cmd) in enumerate(refusals):
        d = os.path.join(hb, f"case{k}")
        open(os.path.join(hb, f"case{k}.name"), "w").write(name + "\n")
        R.start(f"H4 {name}", cmd, d)
    R.wait_all()
    R.say(f"# done {time.strftime('%Y-%m-%d %H:%M:%S %Z')}")
    R.say(subprocess.run(["df", "-h", O], capture_output=True, text=True).stdout.strip())


if __name__ == "__main__":
    main()
