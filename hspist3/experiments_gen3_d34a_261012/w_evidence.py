#!/usr/bin/env python3
"""##CHRIS 2026-10-12 (plan-author decision 34, addition; 261012 sec. 4.7.73): the acceptance evidence of the energy-transfer scene in
the gen3 window (a DISPLAY change: the spring-constant row, the pause at a divider's limit, the UI at the display's scale), with the
binaries of a committed tree (w_build.sh). w_tables.py prints the acceptance from this folder. As decision 19's d19_evidence.py: at most
8 simulation processes at once; >= 4 GiB free before every launch; every run in its own folder; nothing deleted; no window on the
screen -- every window run uses SDL's dummy video driver and the software renderer, with a link to the font Roboto-Regular.ttf in its
folder; screenshots are written as BMP into a scratch folder (--shots, outside the repository) and copied here as PNG.
Parts (decision 34 item 4 as decision 19 item 2, and the addition's item E):
  H0  rule 4: gen2 byte identity with 7b08827 (audit_runs_261007.py on ctrl_min and ctrl_leg, the default build and --engine=gen2)
  H1  the M1 and M2 harness outputs against the committed ones
  P4  headless gen3 unchanged: P4a's start-A seed-1 command, event hash bd25077c5f1c9c9f, its data files byte-identical
  RP  headless gen3 unchanged: stage A's replay cell m2_driven_pi8, variant a, hash 694fdf127f35d4b4, MATCH
  ST  the display does not steer: decision 19's ST scene (N = 400, eta 0.31, the right piston stepping in, 40 sigma-time; --demo
      capture mode) under gen3 headless and with the window at TWO sizes: A 1024x768 points at UI scale 1, and B 4112x2462 at UI
      scale 2 (Chris's display, 2056 x 1231 points of 2 device pixels, emulated): equal event hashes (= decision 19's
      daea5de1378282dc), the output files byte-identical
  SP  item E: the energy-transfer scene with a spring (N = 400 on the right, r = 0.2159 sigma, eta 0.116; spring wall 0, L0 = 25.34
      sigma, M = 50) at k = 11.52 kT/sigma^2 (M2's spring_pi8) and k = 115.2, free from the first step (--auto-release), 200 sigma-time
      after the release (--time=200); with the window (A) and headless: gen3's ledger line (E = KE + spring) conserved to rounding,
      equal event hashes, the output files byte-identical
  SF  item B through the F3 panel: the same scene at k = 11.52 rebuilt by the panel's APPLY (--gui-keys=1:apply-gen3, the k row read
      from the command line as the panel does), window B; its ledger line; a screenshot after step 3000 (the HUD's spring line, the
      F3 panel's k row, the UI at Chris's scale)
  LM  item C: the same scene without a spring (Chris's case), window B: the window pauses before the divider reaches the left outer
      wall ([GUI-LIMIT], the run record), a screenshot of the pause; LO, the same scene with decision 19's binary bin_7f6adca: the
      "Invalid wall_x" exit (information: the behaviour Chris met)
usage (from hspist3/ of the engine-gen3 worktree):
  python3 experiments_gen3_d34a_261012/w_evidence.py --bin-dir <bin dir> --expect-sha <00ALLINONE sha256> --out <evidence dir> --shots <scratch dir>
"""
import argparse, hashlib, math, os, shutil, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
GATE = "/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3/cluster/resched_gate_261005"
P4A1 = os.path.join(HS, "experiments_gen3_p4_261010", "data_06e7248", "N400_eta0.740", "A", "seed1")
OLD_BIN = os.path.join(HS, "experiments_gen3_d19_261010", "bin_7f6adca", "00ALLINONE")
CAP, FREE_MIN = 8, 4 * 1024 ** 3
DUMMY = {"SDL_VIDEODRIVER": "dummy", "SDL_RENDER_DRIVER": "software"}
WIN_A = []                                                         # 1024 x 768 points (the dummy display), UI scale 1
WIN_B = ["--gui-window=4112x2462", "--gui-ui-scale=2"]             # Chris's display: 2056 x 1231 points of 2 device pixels
SPRING_KS = ("11.52", "115.2")


def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()


def p4_cmd(src, binp, d):
    old = open(os.path.join(src, "command.txt")).read().split()
    return [binp] + [f"--speed-sound-run-dir={d}" if a.startswith("--speed-sound-run-dir=") else a for a in old[1:]]


def st_scene(binp, shot):
    """decision 19's ST scene (d19_evidence.py scene(), N = 400, demo capture run): H 20, L0 25.34, the right piston stepping in"""
    return [binp, "--show-simulation", "--mode=edmd", "--engine=gen3", "--edmd-acc=0", "--seed-drift-order=drift-first", "--kbt1",
            "--particles=400", "--particles-boxes=200,200", "--particle-radius=0.5", "--l0=25.3400", "--height=20.0000",
            "--num-walls=1", "--wall-positions=25.3400", "--wall-mass-factors=50", "--wall-hold-steps=600", "--fixed-dt=0.4", "--seed=9300",
            "--demo", f"--demo-shot={shot},700", "--piston-right-protocol-mode=step", "--velocity-right-piston-step=0.02",
            "--max-right-piston-travel=3", "--time=40", "--gui-stats"]


def et_scene(binp, k):
    """the energy-transfer scene: N = 400 in the right compartment (r = 0.2159 sigma, eta 0.116), the divider free from the first
    step; k: the spring (wall 0, L0 = 25.34 sigma) in kT/sigma^2, or None (no spring declared: Chris's case)"""
    spring = [f"--spring-k-sigma={k}", "--spring-wall=0", "--spring-eq=25.34"] if k else []
    return [binp, "--show-simulation", "--mode=edmd", "--engine=gen3", "--edmd-acc=0", "--seed-drift-order=drift-first", "--kbt1",
            "--experiment=energy_transfer", "--particles=400", "--particles-boxes=0,400", "--particle-radius=0.2159", "--l0=25.34",
            "--height=20", "--num-walls=1", "--wall-positions=25.34", "--wall-mass-factors=50"] + spring + \
           ["--auto-release", "--fixed-dt=0.4", "--seed=9300", "--gui-stats"]


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
    def wait_all(s):
        while s.run: s.reap(); time.sleep(0.5)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--bin-dir", required=True); ap.add_argument("--expect-sha", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--shots", required=True); a = ap.parse_args()
    B, O, SH = os.path.abspath(a.bin_dir), os.path.abspath(a.out), os.path.abspath(a.shots)
    NEW = os.path.join(B, "00ALLINONE")
    if os.path.exists(O) and os.listdir(O): sys.exit(f"STOP: {O} is not empty -- the evidence is never repeated into the same folder")
    if os.path.exists(SH) and os.listdir(SH): sys.exit(f"STOP: {SH} is not empty")
    if SH.startswith(HS + os.sep): sys.exit("STOP: the screenshot scratch folder must be outside the repository")
    os.makedirs(O, exist_ok=True); os.makedirs(SH, exist_ok=True)
    if sha(NEW) != a.expect_sha: sys.exit("STOP: the new build's SHA-256")
    R = Runner(O)
    R.say(f"# decision 34 addition evidence {time.strftime('%Y-%m-%d %H:%M:%S %Z')}: {NEW} "
          f"({subprocess.run([NEW, '--version'], capture_output=True, text=True).stdout.splitlines()[0]})")
    for c in (["df", "-h", O], ["sysctl", "vm.swapusage"], ["uptime"]): R.say(subprocess.run(c, capture_output=True, text=True).stdout.strip())
    # ---- the long one first, then the gates, then the window runs
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
    # ST: headless, window A, window B
    for mode, extra in (("headless", None), ("window_A", WIN_A), ("window_B", WIN_B)):
        d = os.path.join(O, "steer", mode)
        cmd = st_scene(NEW, os.path.join(SH, f"ST_{mode}_demo.bmp"))
        if extra is None: cmd.insert(2, "--headless")
        else: cmd += extra
        R.start(f"ST gen3 {mode}", cmd, d, env=None if extra is None else DUMMY, font=extra is not None)
    # SP: the spring scene at two k, window A and headless
    for k in SPRING_KS:
        for mode in ("headless", "window_A"):
            d = os.path.join(O, "spring", f"k{k}_{mode}")
            cmd = et_scene(NEW, k) + ["--time=200"] + (["--headless"] if mode == "headless" else ["--gui-time-scale=2000"])
            R.start(f"SP k={k} {mode}", cmd, d, env=None if mode == "headless" else DUMMY, font=mode != "headless")
    # SF: through the F3 panel's APPLY, window B, a screenshot after step 3000
    d = os.path.join(O, "spring_f3", "k11.52_window_B")
    R.start("SF k=11.52 F3 APPLY window B", et_scene(NEW, "11.52") + ["--time=200", "--gui-time-scale=2000", "--gui-keys=1:apply-gen3"] + WIN_B
            + [f"--gui-shot=3000:{os.path.join(SH, 'SF_k11.52_window_B_step3000.bmp')}"], d, env=DUMMY, font=True)
    # LM: no spring, window B: the pause and its screenshot; LO: decision 19's binary on the same scene
    d = os.path.join(O, "limit", "new_window_B")
    R.start("LM no spring, window B (pause)", et_scene(NEW, None) + ["--gui-time-scale=1000", "--gui-run-seconds=30"] + WIN_B
            + [f"--gui-shot=limit:{os.path.join(SH, 'LM_window_B_pause.bmp')}"], d, env=DUMMY, font=True)
    d = os.path.join(O, "limit", "old_bin_7f6adca_window")
    R.start("LO no spring, bin_7f6adca (information)", et_scene(OLD_BIN, None) + ["--gui-time-scale=1000", "--gui-run-seconds=30"], d,
            env=DUMMY, font=True)
    R.wait_all()
    # the screenshots as PNG (the BMPs stay in the scratch folder)
    try:
        from PIL import Image
        os.makedirs(os.path.join(O, "screenshots"), exist_ok=True)
        for f in sorted(os.listdir(SH)):
            if not f.endswith(".bmp"): continue
            im = Image.open(os.path.join(SH, f)).convert("RGB")
            dst = os.path.join(O, "screenshots", f[:-4] + ".png")
            im.save(dst, optimize=True)
            R.say(f"  screenshot {f} ({im.size[0]} x {im.size[1]}) -> {os.path.relpath(dst, O)} ({os.path.getsize(dst) / 1e6:.2f} MB)")
    except Exception as e:   # the tables script reports a missing PNG
        R.say(f"  screenshots: PNG conversion failed: {e}")
    R.say(subprocess.run(["uptime"], capture_output=True, text=True).stdout.strip())
    R.say(f"# done {time.strftime('%Y-%m-%d %H:%M:%S %Z')}")


if __name__ == "__main__":
    main()
