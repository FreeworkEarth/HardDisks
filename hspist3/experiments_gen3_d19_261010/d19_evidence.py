#!/usr/bin/env python3
"""##CHRIS 2026-10-10 (plan-author decision 19, item 2; 261012 sec. 4.7.50): the acceptance evidence of gen3 in the window, a provisional
DISPLAY feature that no campaign script uses, with the binaries of a committed tree (d19_build.sh). d19_tables.py prints the acceptance
from this folder. Parts (at most 8 simulation processes at once; >= 4 GiB free before every launch; every run in its own folder; nothing
deleted; no window on the screen: every window run uses SDL's dummy video driver and the software renderer, with a link to the
font Roboto-Regular.ttf in its folder, where the window loads it -- the first start of this runner, without it, stopped both window runs
at "Failed to load font" and was stopped at once; its partial folder is kept in the session scratchpad):
  H0  rule 4: gen2 byte identity with 7b08827 after the 00ALLINONE.c change -- audit_runs_261007.py on ctrl_min and ctrl_leg with the
      default build and with --engine=gen2, as stage H's H0 (2 runner entries of 2 simulations each)
  H1  the M1 and M2 harness outputs against the committed ones, as stage H's H1
  P4  headless gen3 unchanged: P4a's start-A seed-1 command exactly as P4 ran it, with the new build: event hash bd25077c5f1c9c9f
      (about 20 min)
  RP  headless gen3 unchanged: stage A's replay cell m2_driven_pi8, variant a (the M2 harness's dump replayed through the driver,
      readers at every target, one stop per interval): the harness's event hash 694fdf127f35d4b4, MATCH
  ST  the display does not steer: the interactive two-compartment scene (N = 400, eta 0.31, a right piston stepping in) with the
      window (dummy driver, no input: --demo) and with --headless, the same seed and protocol, 40 sigma-time: equal event hashes; and
      the same pair under gen2 (information: the gen2 scene's state, no event hash)
  PF  ALONE, after everything else, one run at a time: the scene with the window (dummy driver), started by --auto-release (the
      divider free from the start; no --demo, so the window's own pacing: steps per frame = the time scale x the frame time, the frame
      time capped at 0.1 s), at N = 5000 and 8000, eta 0.31 (the ST geometry scaled: H = 20 sqrt(N/400), L0 = 25.34 sqrt(N/400)),
      gen2 and gen3, at the default time scale 1 and at the largest, 10000 (as fast as the engine goes); 30 wall seconds each
      (--gui-run-seconds=30 --gui-stats): frames per second and sigma-time per wall second from the [GUI-STATS] line
usage (from hspist3/ of the engine-gen3 worktree):
  python3 experiments_gen3_d19_261010/d19_evidence.py --bin-dir <bin dir> --expect-sha <00ALLINONE sha256> --out <evidence dir>
"""
import argparse, hashlib, math, os, shutil, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
GATE = "/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3/cluster/resched_gate_261005"
P4A1 = os.path.join(HS, "experiments_gen3_p4_261010", "data_06e7248", "N400_eta0.740", "A", "seed1")
CAP, FREE_MIN = 8, 4 * 1024 ** 3
DUMMY = {"SDL_VIDEODRIVER": "dummy", "SDL_RENDER_DRIVER": "software"}


def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()


def p4_cmd(src, binp, d):
    old = open(os.path.join(src, "command.txt")).read().split()
    return [binp] + [f"--speed-sound-run-dir={d}" if a.startswith("--speed-sound-run-dir=") else a for a in old[1:]]


def scene(binp, eng, N, d, extra, demo=True):
    """the interactive two-compartment scene at eta 0.31 (the ST cell at N = 400: H 20, L0 25.34), N/2 disks per compartment; demo: a
    capture run (--demo --demo-shot: it starts by itself and paces 64 steps per frame, as the WIP's ST pair); otherwise --auto-release
    (the scene starts by itself, the divider free from the start) and the window's own pacing: the time scale alone sets the steps"""
    f = math.sqrt(N / 400.0); H, L0 = 20.0 * f, 25.34 * f
    return [binp, "--show-simulation", "--mode=edmd", f"--engine={eng}", "--edmd-acc=0", "--seed-drift-order=drift-first", "--kbt1",
            f"--particles={N}", f"--particles-boxes={N // 2},{N // 2}", "--particle-radius=0.5", f"--l0={L0:.4f}", f"--height={H:.4f}",
            "--num-walls=1", f"--wall-positions={L0:.4f}", "--wall-mass-factors=50", "--wall-hold-steps=600", "--fixed-dt=0.4", "--seed=9300"] \
        + (["--demo", f"--demo-shot={os.path.join(d, 'shot.bmp')},700"] if demo else ["--auto-release"]) + extra


class Runner:
    def __init__(s, out): s.out = out; s.run = []; s.log = open(os.path.join(out, "runner.log"), "a")
    def say(s, m): print(m, flush=True); s.log.write(m + "\n"); s.log.flush()
    def load(s): return sum(w for p, w, *_ in s.run if p.poll() is None)
    def start(s, name, cmd, d, weight=1, cwd=None, env=None, font=False):
        while s.load() + weight > CAP: s.reap(); time.sleep(0.5)
        free = shutil.disk_usage(s.out).free
        if free < FREE_MIN:
            s.say(f"**STOPPED: free space {free / 1024 ** 3:.2f} GiB < 4 GiB before {name}; nothing more is launched**"); s.wait_all(); sys.exit(3)
        if os.path.exists(d) and os.listdir(d): s.say(f"STOP: {d} is not empty"); s.wait_all(); sys.exit(3)
        os.makedirs(d, exist_ok=True)
        if font: os.symlink(os.path.join(HS, "Roboto-Regular.ttf"), os.path.join(d, "Roboto-Regular.ttf"))   # the window loads it from its folder
        with open(os.path.join(d, "command.txt"), "w") as fh:
            fh.write((" ".join(f"{k}={v}" for k, v in (env or {}).items()) + " " if env else "") + " ".join(cmd) + "\n")
        p = subprocess.Popen(cmd, cwd=cwd or d, stdout=open(os.path.join(d, "run.log"), "w"), stderr=subprocess.STDOUT,
                             env=dict(os.environ, **(env or {})))
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
    R = Runner(O)
    R.say(f"# decision 19 item 2 evidence {time.strftime('%Y-%m-%d %H:%M:%S %Z')}: {NEW} "
          f"({subprocess.run([NEW, '--version'], capture_output=True, text=True).stdout.splitlines()[0]})")
    for c in (["df", "-h", O], ["sysctl", "vm.swapusage"], ["uptime"]): R.say(subprocess.run(c, capture_output=True, text=True).stdout.strip())
    # ---- phase 1: P4 (the long one), H0 (2 x 2), H1, RP, ST (4)
    R.start("P4 P4a A seed 1, new build", p4_cmd(P4A1, NEW, os.path.join(O, "p4a_A1")), os.path.join(O, "p4a_A1"))
    with open(os.path.join(B, "bin_engine_gen2.sh"), "w") as fh: fh.write(f'#!/bin/bash\nexec "{NEW}" "$@" --engine=gen2\n')
    os.chmod(os.path.join(B, "bin_engine_gen2.sh"), 0o755)
    for tag, b in (("default", NEW), ("gen2flag", os.path.join(B, "bin_engine_gen2.sh"))):
        R.start(f"H0 {tag}", [sys.executable, os.path.join(GATE, "audit_runs_261007.py"), "run", "--bin", b, "--out", os.path.join(O, "h0", tag),
                              "--cases", "ctrl_min,ctrl_leg", "--jobs", "2"], os.path.join(O, "h0", f"{tag}_runner"), weight=2, cwd=HS)
    R.start("H1 M1 and M2 harness outputs", ["/bin/bash", "-c", f'"{B}/gen3_m1" audit > m1_audit_output.txt && "{B}/gen3_m2" audit --quick > m2_audit_quick_output.txt'
                                              f' && "{B}/gen3_m2" audit > m2_audit_output.txt'], os.path.join(O, "h1"))
    rp = os.path.join(O, "replay")
    R.start("RP m2 dump, then the replay of m2_driven_pi8 variant a", ["/bin/bash", "-c",
            f'"{B}/gen3_m2" dump "{rp}" > dump_output.txt && "{NEW}" --engine-replay="{rp}/m2_driven_pi8.g3state" --replay-read-every=1 --replay-substeps=1'],
            rp)
    for eng in ("gen3", "gen2"):
        for mode in ("window", "headless"):
            d = os.path.join(O, "steer", f"{eng}_{mode}")
            extra = ["--piston-right-protocol-mode=step", "--velocity-right-piston-step=0.02", "--max-right-piston-travel=3", "--time=40", "--gui-stats"]
            cmd = scene(NEW, eng, 400, d, extra)
            if mode == "headless": cmd.insert(2, "--headless")
            R.start(f"ST {eng} {mode}", cmd, d, env=DUMMY if mode == "window" else None, font=mode == "window")
    R.wait_all()
    # ---- phase 2: PF, alone, one run at a time
    R.say(subprocess.run(["uptime"], capture_output=True, text=True).stdout.strip())
    for N in (5000, 8000):
        for ts in ("1", "10000"):
            for eng in ("gen2", "gen3"):
                d = os.path.join(O, "perf", f"N{N}_ts{ts}_{eng}")
                R.wait(R.start(f"PF N {N} time scale {ts} {eng}", scene(NEW, eng, N, d, [f"--gui-time-scale={ts}", "--gui-run-seconds=30", "--gui-stats"],
                                                                                   demo=False), d, env=DUMMY, font=True))
    R.say(subprocess.run(["uptime"], capture_output=True, text=True).stdout.strip())
    R.say(f"# done {time.strftime('%Y-%m-%d %H:%M:%S %Z')}")


if __name__ == "__main__":
    main()
