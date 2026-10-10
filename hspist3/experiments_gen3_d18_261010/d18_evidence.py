#!/usr/bin/env python3
"""##CHRIS 2026-10-10 (plan-author decision 18, items 1 and 2; 261012 sec. 4.7.42): the acceptance evidence of
  item 1 -- checkpoint/restart carries --gen3-psi6-compartments, and
  item 2 -- M4c, the turned orientation in the M4 lattice family (--gen3-lattice-orientation=0|turned),
with the binaries of the committed tree (d18_build.sh). d18_tables.py prints the acceptance from this folder.
Parts (at most 8 simulation processes in all; >= 4 GiB free before every launch; every run in its own folder; nothing deleted):
  H   stage H's evidence runner, unchanged (experiments_gen3_h_261009/stageH_evidence.py), its JOBS set to 5 from outside (so at most
      5 - 2 + 4 = 7 simulations): rule 4 (H0), the M1/M2 harness outputs (H1), the engine checkpoint test (H2), the driver's restarts
      (H3), the refusals (H4); d18_tables.py prints stage H's tables with stageH_tables.py, unchanged
  I   alongside H: P4a's start-A seed-1 run exactly as P4 ran it, with the new build and no new flag (the default): its event hash must
      be P4's bd25077c5f1c9c9f and its snapshots (the initial state at t = 0 among them), psi6(t) and virial blocks P4's
  R   item 1: stage H's cells B and C with --gen3-psi6-compartments (B also with --gen3-virial-blocks=50 --gen3-snapshots=500): the
      uninterrupted run U and, per checkpoint (B: hold:60000, record:6000; C: hold:15000, record:3000), C (writes it and goes on), R1
      (writes it and stops) and R2 (restarts from R1's); the refusals: R1's checkpoint of B at hold:60000 (columns on) restarted
      without the flag, and stage H's own checkpoint of cell B at hold:60000 (columns off) restarted with it
  M   item 2:
      ident  existing commands, the previous committed-tree build 0d7f1ca (experiments_gen3_psi6c_261010/bin_0d7f1ca) against the new
             one: stage C's grid cell N = 400, eta 0.70 as stage C ran it ("long", its seed s1: 400 sigma-time, the initial state
             dumped; also against stage C's own evidence of f39e485: information); the same cell with --gen3-lattice=perpendicular, held
             10 sigma-time; P4's start-B (M4b) N = 400 seed-0 command held 10 sigma-time
      grid   M4's acceptance grid N = 100, 400, 900, 1600 (H = 10 sqrt(N/100)) x eta 0.60, 0.70, 0.716, 0.74, 0.78, 0.85, L0 = N pi /
             (8 H eta) exactly (stage C's command, --seed 20261227), with --gen3-lattice-orientation=turned --gen3-psi6-compartments:
             per cell "long" (seed s1: 400 sigma-time held, 1 released, the initial state dumped), "same" (s1, 1 sigma-time held) and
             "other" (s2, 1 sigma-time held), as stage C; seeds run_seed(20261227, 0, k, 0 and 1)
      load   the same 24 cells held 1 sigma-time with --gen3-psi6-compartments: the default family (run_seed(20261227, 1, k, 0)) and
             --gen3-lattice-orientation=0 (run_seed(20261227, 2, k, 0)): the orientation at the load
  The base 20261227 is in no script, command or seed list (only as time values, 20261227.xx, in one Szilard trace CSV).
usage (from hspist3/ of the engine-gen3 worktree):
  python3 experiments_gen3_d18_261010/d18_evidence.py --bin-dir <bin dir> --expect-sha <00ALLINONE sha256> --out <evidence dir>
"""
import argparse, hashlib, math, os, shutil, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
for sub in ("experiments_gen3_h_261009", "experiments_gen3_m4_261009", "experiments_gen3_p1_261009", "experiments_gen3_p2_261009"):
    sys.path.insert(0, os.path.join(HS, sub))
import stageH_evidence as SH
import stageC_evidence as SC
T = SC.T
OLD = os.path.join(HS, "experiments_gen3_psi6c_261010", "bin_0d7f1ca", "00ALLINONE")
OLD_SHA = "082874142c3c1e6a63ef99e02269123eb90cde09cfd10ef549c3017cbbfc32c1"
P4D = os.path.join(HS, "experiments_gen3_p4_261010", "data_06e7248", "N400_eta0.740")
P4A1, P4B0 = os.path.join(P4D, "A", "seed1"), os.path.join(P4D, "B", "seed0")
NS_LIST, ETAS, BASE = (100, 400, 900, 1600), (0.60, 0.70, 0.716, 0.74, 0.78, 0.85), 20261227
LONG, SHORT = 24000, 60                                       # 400 sigma-time and 1 sigma-time held (60 steps per sigma-time)
RB = ["--gen3-psi6-compartments", "--gen3-virial-blocks=50", "--gen3-snapshots=500"]; RC = ["--gen3-psi6-compartments"]
CAP, FREE_MIN = 8, 4 * 1024 ** 3


def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()


def patch(cmd, repl): return SH.patch(cmd, repl)


def grid():
    out = []
    for N in NS_LIST:
        for eta in ETAS:
            k = len(out); H = 10.0 * math.sqrt(N / 100.0); L0 = N * math.pi / (8.0 * H * eta)
            out.append(dict(k=k, N=N, eta=eta, H=H, L0=L0, name=f"N{N}_eta{eta:.3f}", s1=T.run_seed(BASE, 0, k, 0), s2=T.run_seed(BASE, 0, k, 1),
                            sd=T.run_seed(BASE, 1, k, 0), s0=T.run_seed(BASE, 2, k, 0)))
    return out


def gcmd(binp, d, c, seed, hold, extra):
    return patch(SC.cmd(binp, d, c, seed, hold), {"--seed=": str(BASE)}) + list(extra)


def stagec_cell():
    return [c for c in SC.cells() if c["N"] == 400 and c["lab"] == "0.70"][0]


def p4_cmd(src, binp, d, repl=None):
    old = open(os.path.join(src, "command.txt")).read().split()
    out = [binp]
    for a in old[1:]:
        k = a.split("=", 1)[0] + "="
        if k == "--speed-sound-run-dir=": a = f"--speed-sound-run-dir={d}"
        elif repl and k in repl: a = k + repl[k]
        out.append(a)
    return out


class Runner:
    def __init__(s, out): s.out = out; s.run = []; s.log = open(os.path.join(out, "runner.log"), "a")
    def say(s, m): print(m, flush=True); s.log.write(m + "\n"); s.log.flush()
    def load(s): return sum(w for p, w, *_ in s.run if p.poll() is None)
    def start(s, name, cmd, d, weight=1, cwd=None):
        while s.load() + weight > CAP: s.reap(); time.sleep(0.5)
        free = shutil.disk_usage(s.out).free
        if free < FREE_MIN:
            s.say(f"**STOPPED: free space {free / 1024 ** 3:.2f} GiB < 4 GiB before {name}; nothing more is launched**"); s.wait_all(); sys.exit(3)
        if os.path.exists(d) and os.listdir(d): s.say(f"STOP: {d} is not empty"); s.wait_all(); sys.exit(3)
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, "command.txt"), "w") as fh: fh.write(" ".join(cmd) + "\n")
        p = subprocess.Popen(cmd, cwd=cwd or d, stdout=open(os.path.join(d, "run.log"), "w"), stderr=subprocess.STDOUT)
        s.run.append((p, weight, name, d, time.time())); s.say(f"  start {name} (load {s.load()}, free {free / 1024 ** 3:.1f} GiB)")
        return p
    def reap(s):
        for j in list(s.run):
            p, w, name, d, ts = j
            if p.poll() is None: continue
            s.run.remove(j); open(os.path.join(d, "exit_code.txt"), "w").write(f"{p.returncode}\n")
            s.say(f"  {name}: exit {p.returncode}, {time.time() - ts:.0f} s")
    def wait(s, p):
        while p.poll() is None: s.reap(); time.sleep(0.5)
        s.reap()
    def wait_all(s):
        while s.run: s.reap(); time.sleep(0.5)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--bin-dir", required=True); ap.add_argument("--expect-sha", required=True)
    ap.add_argument("--out", required=True); a = ap.parse_args()
    B, O = os.path.abspath(a.bin_dir), os.path.abspath(a.out); NEW = os.path.join(B, "00ALLINONE")
    if os.path.exists(O) and os.listdir(O): sys.exit(f"STOP: {O} is not empty -- the evidence is never repeated into the same folder")
    os.makedirs(O, exist_ok=True)
    if sha(NEW) != a.expect_sha: sys.exit("STOP: the new build's SHA-256")
    if sha(OLD) != OLD_SHA: sys.exit("STOP: the previous build's SHA-256")
    R = Runner(O)
    R.say(f"# decision 18 items 1-2 evidence {time.strftime('%Y-%m-%d %H:%M:%S %Z')}: new {NEW} "
          f"({subprocess.run([NEW, '--version'], capture_output=True, text=True).stdout.splitlines()[0]}); old {OLD}")
    R.say(subprocess.run(["df", "-h", O], capture_output=True, text=True).stdout.strip())
    R.say(subprocess.run(["sysctl", "vm.swapusage"], capture_output=True, text=True).stdout.strip())
    # ---- phase 1: H (stage H's runner, JOBS = 5: at most 7 simulations) and I (1)
    hcode = (f"import sys; sys.path.insert(0, {os.path.join(HS, 'experiments_gen3_h_261009')!r}); import stageH_evidence as SH; SH.JOBS = 5; "
             f"sys.argv = ['stageH_evidence.py', '--bin-dir', {B!r}, '--out', {os.path.join(O, 'stageH')!r}]; SH.main()")
    ph = R.start("H stage H's evidence runner (JOBS 5)", [sys.executable, "-c", hcode], os.path.join(O, "stageH_runner"), weight=7, cwd=HS)
    R.start("I P4a A seed 1, new build, default", p4_cmd(P4A1, NEW, os.path.join(O, "m4c", "ident", "p4a_A1")), os.path.join(O, "m4c", "ident", "p4a_A1"))
    R.wait(ph)
    # ---- phase 2a: R's U, C, R1; M's runs
    cells = SH.cells(NEW); pend = []
    for cell, extra in (("B", RB), ("C", RC)):
        cmd, cks = cells[cell]; base = os.path.join(O, "ckpt", cell)
        R.start(f"R {cell} U", cmd(os.path.join(base, "U")) + extra, os.path.join(base, "U"))
        for ck in cks:
            t = ck.replace(":", "")
            R.start(f"R {cell} C {ck}", cmd(os.path.join(base, f"C_{t}")) + extra + [f"--gen3-checkpoint={ck}:{base}/ck_C_{t}.bin"], os.path.join(base, f"C_{t}"))
            p = R.start(f"R {cell} R1 {ck}", cmd(os.path.join(base, f"R1_{t}")) + extra + [f"--gen3-checkpoint={ck}:{base}/ck_R1_{t}.bin", "--gen3-checkpoint-stop"],
                        os.path.join(base, f"R1_{t}"))
            pend.append((cell, cmd, extra, base, t, p))
    c = stagec_cell(); md = os.path.join(O, "m4c", "ident")
    for tag, b in (("old", OLD), ("new", NEW)):
        R.start(f"M ident stage C cell N400 eta0.70 long ({tag})", SC.cmd(b, os.path.join(md, f"stageC_{tag}"), c, c["s1"], LONG), os.path.join(md, f"stageC_{tag}"))
        R.start(f"M ident perpendicular ({tag})", SC.cmd(b, os.path.join(md, f"perp_{tag}"), c, c["s1"], 600) + ["--gen3-lattice=perpendicular"],
                os.path.join(md, f"perp_{tag}"))
        R.start(f"M ident M4b P4 B seed 0 ({tag})", p4_cmd(P4B0, b, os.path.join(md, f"m4b_{tag}"), {"--wall-hold-steps=": "600"}), os.path.join(md, f"m4b_{tag}"))
    for g in grid():
        gd = os.path.join(O, "m4c", "grid", g["name"])
        R.start(f"M grid {g['name']} long", gcmd(NEW, os.path.join(gd, "long"), g, g["s1"], LONG, ["--gen3-lattice-orientation=turned", "--gen3-psi6-compartments"]),
                os.path.join(gd, "long"))
        R.start(f"M grid {g['name']} same", gcmd(NEW, os.path.join(gd, "same"), g, g["s1"], SHORT, ["--gen3-lattice-orientation=turned", "--gen3-psi6-compartments"]),
                os.path.join(gd, "same"))
        R.start(f"M grid {g['name']} other", gcmd(NEW, os.path.join(gd, "other"), g, g["s2"], SHORT, ["--gen3-lattice-orientation=turned", "--gen3-psi6-compartments"]),
                os.path.join(gd, "other"))
        ld = os.path.join(O, "m4c", "load", g["name"])
        R.start(f"M load {g['name']} default", gcmd(NEW, os.path.join(ld, "default"), g, g["sd"], SHORT, ["--gen3-psi6-compartments"]), os.path.join(ld, "default"))
        R.start(f"M load {g['name']} 0", gcmd(NEW, os.path.join(ld, "zero"), g, g["s0"], SHORT, ["--gen3-lattice-orientation=0", "--gen3-psi6-compartments"]),
                os.path.join(ld, "zero"))
    for cell, cmd, extra, base, t, p in pend: R.wait(p)
    # ---- phase 2b: R2 (after its R1) and the refusals
    for cell, cmd, extra, base, t, p in pend:
        R.start(f"R {cell} R2 {t}", cmd(os.path.join(base, f"R2_{t}")) + extra + [f"--gen3-restart={base}/ck_R1_{t}.bin"], os.path.join(base, f"R2_{t}"))
    cmdB = cells["B"][0]; rb = os.path.join(O, "ckpt", "refusals")
    R.start("R refusal: columns-on checkpoint, run without the flag",
            cmdB(os.path.join(rb, "on_ckpt_off_run")) + ["--gen3-virial-blocks=50", "--gen3-snapshots=500", f"--gen3-restart={O}/ckpt/B/ck_R1_hold60000.bin"],
            os.path.join(rb, "on_ckpt_off_run"))
    R.start("R refusal: stage H's columns-off checkpoint, run with the flag",
            cmdB(os.path.join(rb, "off_ckpt_on_run")) + ["--gen3-psi6-compartments", f"--gen3-restart={O}/stageH/h3/B/ck_R1_hold60000.bin"],
            os.path.join(rb, "off_ckpt_on_run"))
    R.wait_all()
    R.say(f"# done {time.strftime('%Y-%m-%d %H:%M:%S %Z')}")
    R.say(subprocess.run(["df", "-h", O], capture_output=True, text=True).stdout.strip())
    R.say(subprocess.run(["sysctl", "vm.swapusage"], capture_output=True, text=True).stdout.strip())


if __name__ == "__main__":
    main()
