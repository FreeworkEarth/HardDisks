#!/usr/bin/env python3
"""##CHRIS 2026-10-10 (plan-author decision 17, "after P5"; 261012 sec. 4.7.40): the acceptance evidence of the PROVISIONAL psi6 per
compartment (00ALLINONE.c: psi6_compartments, --gen3-psi6-compartments, --gen3-psi6-snapshots), with the binaries of a committed
tree (experiments_gen3_m3_261009/build_clean.sh <commit> <dir>). The plan author's acceptance, checked by psi6c_tables.py:
  - the event hash unchanged with the columns on and off;
  - the existing columns byte-identical;
  - rule 4 identical (gen2 byte identity with 7b08827: ctrl_min, ctrl_leg; the default build and --engine=gen2);
  - the M1/M2 harness outputs identical (the committed ones);
  - on P4a's snapshots the new per-compartment psi6 separates the turned and the unturned compartment.
The runs:
  S0  rule 4: resched_gate_261005/audit_runs_261007.py run --cases ctrl_min,ctrl_leg --jobs 2, the new build by default and with
      --engine=gen2 (as stage H's H0)
  S1  the M1/M2 harnesses of the new build: gen3_m1 audit, gen3_m2 audit --quick, gen3_m2 audit (as stage H's H1)
  S2  a short command three times: P4a's start-A seed-1 command (N = 400, eta 0.740, the M4 lattice, its seed) held 2000 sigma-time
      (120000 steps) with --gen3-snapshots=500: "old" with P4's frozen binary (06e7248, bin_m4b), "off" with the new build, "on" with
      the new build and --gen3-psi6-compartments
  S3  P4a's start-A seed-1 run exactly as P4 ran it (1e5 sigma-time, snapshots every 2000; experiments_gen3_p4_261010/data_06e7248/
      N400_eta0.740/A/seed1/command.txt) with the new build and --gen3-psi6-compartments: against P4's original (06e7248, columns off)
  S4  --gen3-psi6-snapshots of the new build on P4a's four snapshot files (A seeds 1 and 0, B seeds 0 and 1) and on S2 on's and
      S3's own snapshot files
  S5  the refusals: --engine=gen2 with the flag; the flag with --gen3-checkpoint; the flag in the energy-transfer experiment
RULES (decision 16/17): at most 8 simulation processes (S0's audits count 2 each); before every launch >= 4 GiB free; every run in its
own folder; nothing is compressed, moved or deleted.
usage (from hspist3/ of the engine-gen3 worktree):
  python3 experiments_gen3_psi6c_261010/psi6c_evidence.py --bin-dir <clean build dir> --expect-sha <00ALLINONE sha256> --out <evidence dir>
"""
import argparse, hashlib, os, shutil, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
GATE = "/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3/cluster/resched_gate_261005"
OLD_BIN = os.path.join(HS, "experiments_gen3_m4b_261009", "bin_m4b", "00ALLINONE")
OLD_SHA = "5d9ceb95e716e50d90a69ac8351e8065081c4160ce645c247d8a7f3136cd19ae"
P4A = os.path.join(HS, "experiments_gen3_p4_261010", "data_06e7248", "N400_eta0.740")
ORIG = os.path.join(P4A, "A", "seed1")
SNAPS = (("A", 1), ("A", 0), ("B", 0), ("B", 1))
CAP, FREE_MIN = 8, 4 * 1024 ** 3


def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()


def p4_cmd(binp, d, repl=None, extra=()):
    """P4a's start-A seed-1 command (its command.txt) with another binary and run folder, arguments replaced by prefix, extras added"""
    old = open(os.path.join(ORIG, "command.txt")).read().split()
    out = [binp]
    for a in old[1:]:
        k = a.split("=", 1)[0] + "="
        if a.startswith("--speed-sound-run-dir="): a = f"--speed-sound-run-dir={d}"
        elif repl and k in repl: a = k + repl[k]
        out.append(a)
    return out + list(extra)


class Runner:
    def __init__(s, out): s.out = out; s.run = []; s.log = open(os.path.join(out, "runner.log"), "a")
    def say(s, m): print(m, flush=True); s.log.write(m + "\n"); s.log.flush()
    def load(s): return sum(w for p, w, *_ in s.run if p.poll() is None)
    def start(s, name, cmd, d, weight=1, stdout="run.log"):
        while s.load() + weight > CAP: s.reap(); time.sleep(0.5)
        free = shutil.disk_usage(s.out).free
        if free < FREE_MIN:
            s.say(f"**STOPPED: free space {free / 1024 ** 3:.2f} GiB < 4 GiB before {name}; nothing more is launched**"); s.wait_all(); sys.exit(3)
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, "command.txt"), "w") as fh: fh.write(" ".join(cmd) + "\n")
        p = subprocess.Popen(cmd, cwd=d, stdout=open(os.path.join(d, stdout), "w"), stderr=subprocess.STDOUT)
        s.run.append((p, weight, name, d, time.time())); s.say(f"  start {name} (load {s.load()}, free {free / 1024 ** 3:.1f} GiB)")
        return p
    def reap(s):
        for j in list(s.run):
            p, w, name, d, ts = j
            if p.poll() is None: continue
            s.run.remove(j); open(os.path.join(d, "exit_code.txt"), "w").write(f"{p.returncode}\n")
            s.say(f"  {name}: exit {p.returncode}, {time.time() - ts:.0f} s")
    def wait_all(s):
        while s.run: s.reap(); time.sleep(0.5)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--bin-dir", required=True); ap.add_argument("--expect-sha", required=True)
    ap.add_argument("--out", required=True); a = ap.parse_args()
    B, O = os.path.abspath(a.bin_dir), os.path.abspath(a.out); new = os.path.join(B, "00ALLINONE")
    if os.path.exists(O) and os.listdir(O): sys.exit(f"STOP: {O} is not empty -- the evidence is never repeated into the same folder")
    os.makedirs(O, exist_ok=True)
    if sha(new) != a.expect_sha: sys.exit("STOP: the new build's SHA-256")
    if sha(OLD_BIN) != OLD_SHA: sys.exit("STOP: P4's frozen binary's SHA-256")
    R = Runner(O)
    R.say(f"# psi6 per compartment, acceptance evidence {time.strftime('%Y-%m-%d %H:%M:%S %Z')}: new build {new} "
          f"({subprocess.run([new, '--version'], capture_output=True, text=True).stdout.splitlines()[0]}); old {OLD_BIN}")
    R.say(subprocess.run(["df", "-h", O], capture_output=True, text=True).stdout.strip())
    R.say(subprocess.run(["sysctl", "vm.swapusage"], capture_output=True, text=True).stdout.strip())
    # S3 first (the longest), then S0, S1, S2, S4, S5
    R.start("S3 P4a A seed 1 replay, columns on", p4_cmd(new, os.path.join(O, "s3"), extra=["--gen3-psi6-compartments"]), os.path.join(O, "s3"))
    with open(os.path.join(B, "bin_engine_gen2.sh"), "w") as fh: fh.write(f'#!/bin/bash\nexec "{new}" "$@" --engine=gen2\n')
    os.chmod(os.path.join(B, "bin_engine_gen2.sh"), 0o755)
    for tag, b in (("default", new), ("gen2flag", os.path.join(B, "bin_engine_gen2.sh"))):
        R.start(f"S0 rule 4 {tag}", [sys.executable, os.path.join(GATE, "audit_runs_261007.py"), "run", "--bin", b, "--out", os.path.join(O, "s0", tag),
                                     "--cases", "ctrl_min,ctrl_leg", "--jobs", "2"], os.path.join(O, "s0", f"{tag}_runner"), weight=2)
    R.start("S1 M1/M2 harnesses", ["/bin/bash", "-c", f'"{B}/gen3_m1" audit > m1_audit_output.txt && "{B}/gen3_m2" audit --quick > m2_audit_quick_output.txt'
                                   f' && "{B}/gen3_m2" audit > m2_audit_output.txt'], os.path.join(O, "s1"))
    short = {"--wall-hold-steps=": "120000", "--gen3-snapshots=": "500"}
    for tag, b, extra in (("old", OLD_BIN, []), ("off", new, []), ("on", new, ["--gen3-psi6-compartments"])):
        R.start(f"S2 {tag}", p4_cmd(b, os.path.join(O, "s2", tag), short, extra), os.path.join(O, "s2", tag))
    for st, r in SNAPS:
        sp = [f for f in os.listdir(os.path.join(P4A, st, f"seed{r}")) if f.startswith("snapshots_")]
        R.start(f"S4 P4a {st} seed {r}", [new, f"--gen3-psi6-snapshots={os.path.join(P4A, st, f'seed{r}', sp[0])}"],
                os.path.join(O, "s4", f"p4a_{st}{r}"), stdout="psi6_compartments.csv")
    gen2 = [x for x in p4_cmd(new, os.path.join(O, "s5", "gen2"), short) if not x.startswith(("--engine=", "--gen3-", "--psi6-every"))]
    refusals = (("gen2", gen2 + ["--engine=gen2", "--gen3-psi6-compartments"]),                # as stage H's gen2 refusal: no other gen3 flag
                ("checkpoint", p4_cmd(new, os.path.join(O, "s5", "checkpoint"), short,
                                      ["--gen3-psi6-compartments", f"--gen3-checkpoint=hold:60000:{os.path.join(O, 's5', 'ck.bin')}"])),
                ("energy transfer", [new, "--mode=edmd", "--experiment=energy_transfer", "--headless", "--engine=gen3", "--gen3-psi6-compartments"]))
    for tag, cm in refusals: R.start(f"S5 refusal: {tag}", cm, os.path.join(O, "s5", tag.replace(" ", "_")))
    R.wait_all()
    for tag in ("on",):                                                   # S4 on S2's own snapshots (after S2)
        d = os.path.join(O, "s2", tag); sp = [f for f in os.listdir(d) if f.startswith("snapshots_")]
        R.start(f"S4 S2 {tag}'s own snapshots", [new, f"--gen3-psi6-snapshots={os.path.join(d, sp[0])}"], os.path.join(O, "s4", f"s2_{tag}"),
                stdout="psi6_compartments.csv")
    d = os.path.join(O, "s3"); sp = [f for f in os.listdir(d) if f.startswith("snapshots_")]
    R.start("S4 S3's own snapshots", [new, f"--gen3-psi6-snapshots={os.path.join(d, sp[0])}"], os.path.join(O, "s4", "s3"), stdout="psi6_compartments.csv")
    R.wait_all()
    R.say(f"# done {time.strftime('%Y-%m-%d %H:%M:%S %Z')}")
    R.say(subprocess.run(["df", "-h", O], capture_output=True, text=True).stdout.strip())
    R.say(subprocess.run(["sysctl", "vm.swapusage"], capture_output=True, text=True).stdout.strip())


if __name__ == "__main__":
    main()
