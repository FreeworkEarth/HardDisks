#!/usr/bin/env python3
"""##CHRIS 2026-10-12 (plan-author decision 34, the rest; decision 38 item 4; 261012 sec. 4.7.67): the acceptance evidence of the gen3
window fitted to the screen (a DISPLAY change: the box fitted to the window, the panels and labels in window points, no exit on a scene
that cannot be built, the density of the filled compartment), with the binaries of a committed tree (s_build.sh). s_tables.py prints the
acceptance. As decision 19's d19_evidence.py and the addition's w_evidence.py: at most 8 simulation processes at once; >= 4 GiB free
before every launch; every run in its own folder; nothing deleted; no window on the screen (SDL's dummy video driver and the software
renderer, a link to the font in every window run's folder); screenshots as BMP into a scratch folder outside the repository, copied here
as PNG.
Parts (decision 34 item 4, as decision 19 item 2):
  H0, H1, P4, RP  rule 4; the M1/M2 harness outputs; headless gen3 unchanged (P4a A seed 1, bd25077c5f1c9c9f; the replay cell
      m2_driven_pi8 variant a, 694fdf127f35d4b4)
  ST  the display does not steer: decision 19's ST scene (N = 400, eta 0.31, the right piston stepping in, 40 sigma-time; --demo
      capture mode) under gen3 headless and with the window at TWO sizes, A 1024x768 points at UI scale 1 and B 4112x2462 device pixels
      at UI scale 2 (Chris's display): equal event hashes (= decision 19's daea5de1378282dc), the output files byte-identical
  SS  screenshots at N = 400, 5000 and 8000 (eta 0.31, decision 19's display pair scaled: H = 20 sqrt(N/400), L0 = 25.34 sqrt(N/400)),
      the divider free from the first step, at Chris's screen size in points (2056 x 1231, UI scale 1, so the files stay small; the
      layout in points is the same at UI scale 2), the frame after step 60; and N = 5000 again at UI scale 2 (4112 x 2462 device pixels)
      and in a smaller window (1440 x 900 points): the box re-fitted
  NX  no exit: (a) the F3 panel's APPLY of a scene the seeder cannot place (spring L0 at 40 sigma, walls 5 sigma thick, r = 0.35 sigma:
      400 disks in a compartment that holds 240) is refused with a note and the scene goes on; (b) a scene gen3 cannot build (r = 0.7
      sigma: the diameter exceeds gen3's 32 px cell) pauses the window with the reason (screenshot); (c) information: (b) with the
      addition's binary bin_8bc9ac9 ends with gen3's STOP (exit 2)
usage (from hspist3/ of the engine-gen3 worktree):
  python3 experiments_gen3_d34s_261012/s_evidence.py --bin-dir <bin dir> --expect-sha <00ALLINONE sha256> --out <evidence dir> --shots <scratch dir>
"""
import argparse, hashlib, math, os, shutil, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
GATE = "/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3/cluster/resched_gate_261005"
P4A1 = os.path.join(HS, "experiments_gen3_p4_261010", "data_06e7248", "N400_eta0.740", "A", "seed1")
W1_BIN = os.path.join(HS, "experiments_gen3_d34a_261012", "bin_8bc9ac9", "00ALLINONE")
CAP, FREE_MIN = 8, 4 * 1024 ** 3
DUMMY = {"SDL_VIDEODRIVER": "dummy", "SDL_RENDER_DRIVER": "software"}
WIN_A = []
WIN_B = ["--gui-window=4112x2462", "--gui-ui-scale=2"]
SCREEN = ["--gui-window=2056x1231"]                      # Chris's display in points, UI scale 1
SMALL = ["--gui-window=1440x900"]


def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()


def p4_cmd(src, binp, d):
    old = open(os.path.join(src, "command.txt")).read().split()
    return [binp] + [f"--speed-sound-run-dir={d}" if a.startswith("--speed-sound-run-dir=") else a for a in old[1:]]


def st_scene(binp, shot):
    return [binp, "--show-simulation", "--mode=edmd", "--engine=gen3", "--edmd-acc=0", "--seed-drift-order=drift-first", "--kbt1",
            "--particles=400", "--particles-boxes=200,200", "--particle-radius=0.5", "--l0=25.3400", "--height=20.0000",
            "--num-walls=1", "--wall-positions=25.3400", "--wall-mass-factors=50", "--wall-hold-steps=600", "--fixed-dt=0.4", "--seed=9300",
            "--demo", f"--demo-shot={shot},700", "--piston-right-protocol-mode=step", "--velocity-right-piston-step=0.02",
            "--max-right-piston-travel=3", "--time=40", "--gui-stats"]


def pair_scene(binp, N):
    f = math.sqrt(N / 400.0); H, L0 = 20.0 * f, 25.34 * f
    return [binp, "--show-simulation", "--mode=edmd", "--engine=gen3", "--edmd-acc=0", "--seed-drift-order=drift-first", "--kbt1",
            f"--particles={N}", f"--particles-boxes={N // 2},{N // 2}", "--particle-radius=0.5", f"--l0={L0:.4f}", f"--height={H:.4f}",
            "--num-walls=1", f"--wall-positions={L0:.4f}", "--wall-mass-factors=50", "--auto-release", "--fixed-dt=0.4", "--seed=9300",
            "--gui-time-scale=200", "--gui-run-seconds=8", "--gui-stats"]


def et_refuse(binp):
    return [binp, "--show-simulation", "--mode=edmd", "--engine=gen3", "--edmd-acc=0", "--seed-drift-order=drift-first", "--kbt1",
            "--experiment=energy_transfer", "--particles=400", "--particles-boxes=0,400", "--particle-radius=0.35", "--wall-thickness=5",
            "--l0=25.34", "--height=20", "--num-walls=1", "--wall-positions=25.34", "--wall-mass-factors=50", "--spring-k-sigma=11.52",
            "--spring-wall=0", "--spring-eq=40", "--auto-release", "--fixed-dt=0.4", "--seed=9300", "--gui-time-scale=200",
            "--gui-keys=1:apply-gen3", "--gui-run-seconds=5", "--gui-stats"]


def big_disks(binp):
    return [binp, "--show-simulation", "--mode=edmd", "--engine=gen3", "--edmd-acc=0", "--seed-drift-order=drift-first", "--kbt1",
            "--particles=100", "--particles-boxes=50,50", "--particle-radius=0.7", "--l0=25.34", "--height=20", "--num-walls=1",
            "--wall-positions=25.34", "--wall-mass-factors=50", "--auto-release", "--fixed-dt=0.4", "--seed=9300", "--gui-time-scale=200",
            "--gui-run-seconds=5", "--gui-stats"]


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
        if font: os.symlink(os.path.join(HS, "Roboto-Regular.ttf"), os.path.join(d, "Roboto-Regular.ttf"))
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
    R.say(f"# decision 34 (the rest) evidence {time.strftime('%Y-%m-%d %H:%M:%S %Z')}: {NEW} "
          f"({subprocess.run([NEW, '--version'], capture_output=True, text=True).stdout.splitlines()[0]})")
    for c in (["df", "-h", O], ["sysctl", "vm.swapusage"], ["uptime"]): R.say(subprocess.run(c, capture_output=True, text=True).stdout.strip())
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
    for mode, extra in (("headless", None), ("window_A", WIN_A), ("window_B", WIN_B)):
        d = os.path.join(O, "steer", mode)
        cmd = st_scene(NEW, os.path.join(SH, f"ST_{mode}_demo.bmp"))
        if extra is None: cmd.insert(2, "--headless")
        else: cmd += extra
        R.start(f"ST gen3 {mode}", cmd, d, env=None if extra is None else DUMMY, font=extra is not None)
    for N, tag, win in ((400, "N400_screen", SCREEN), (5000, "N5000_screen", SCREEN), (8000, "N8000_screen", SCREEN),
                        (5000, "N5000_retina", WIN_B + []), (5000, "N5000_small", SMALL)):
        d = os.path.join(O, "screens", tag)
        R.start(f"SS {tag}", pair_scene(NEW, N) + win + [f"--gui-shot=60:{os.path.join(SH, f'SS_{tag}.bmp')}"], d, env=DUMMY, font=True)
    R.start("NX a: the APPLY of a scene the seeder cannot place", et_refuse(NEW), os.path.join(O, "noexit", "apply_refused"), env=DUMMY, font=True)
    R.start("NX b: a scene gen3 cannot build (r = 0.7 sigma)", big_disks(NEW) + SCREEN + [f"--gui-shot=limit:{os.path.join(SH, 'NX_gen3_refused_pause.bmp')}"],
            os.path.join(O, "noexit", "gen3_refused"), env=DUMMY, font=True)
    R.start("NX c: the same with the addition's binary bin_8bc9ac9 (information)", big_disks(W1_BIN), os.path.join(O, "noexit", "gen3_refused_old_bin_8bc9ac9"),
            env=DUMMY, font=True)
    R.wait_all()
    try:
        from PIL import Image
        os.makedirs(os.path.join(O, "screenshots"), exist_ok=True)
        for f in sorted(os.listdir(SH)):
            if not f.endswith(".bmp"): continue
            im = Image.open(os.path.join(SH, f)).convert("RGB")
            dst = os.path.join(O, "screenshots", f[:-4] + ".png")
            im.save(dst, optimize=True)
            R.say(f"  screenshot {f} ({im.size[0]} x {im.size[1]}) -> {os.path.relpath(dst, O)} ({os.path.getsize(dst) / 1e6:.2f} MB)")
    except Exception as e:
        R.say(f"  screenshots: PNG conversion failed: {e}")
    R.say(subprocess.run(["uptime"], capture_output=True, text=True).stdout.strip())
    R.say(f"# done {time.strftime('%Y-%m-%d %H:%M:%S %Z')}")


if __name__ == "__main__":
    main()
