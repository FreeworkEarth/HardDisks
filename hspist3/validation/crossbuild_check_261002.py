#!/usr/bin/env python3
"""##CHRIS 2026-10-02: Task E5 -- one seed per mode through three release binaries, byte comparison of every data file.

Binaries: the new clean build (hspist3/00ALLINONE), the previous clean build (00ALLINONE_e823187_261002) and the first
box-width build (00ALLINONE_5190846dirty_20261014). All three carry the Box_Width_sigma / box_width_sigma column (source
615561c on), so the files are compared byte for byte, with nothing stripped. The energy-transfer summary is compared field by
field, excluding the fields that name the build or the run (build_git, command, timestamp, trace_path).

Same commands as the box-width determinism gate (methods sec. 14.1): pi/8 cell, N = 100, H = L_0 = 10, M = 50 speed of sound
(10 oscillations, HD_KE_TRACE=1), and a held-divider energy-transfer run (12000 steps, HD_PISTON_EVENTS). Output goes to the
directory given as the first argument (the scratchpad); nothing is written to the repo.

usage: python3 hspist3/validation/crossbuild_check_261002.py <outdir>
"""
import csv, filecmp, glob, os, subprocess, sys
HS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BINS = [("new", "00ALLINONE"), ("e823187", "00ALLINONE_e823187_261002"), ("5190846dirty", "00ALLINONE_5190846dirty_20261014")]
OUT = sys.argv[1]
if os.path.exists(OUT): sys.exit(f"{OUT} exists -- choose a fresh directory")

def sos(b, d):
    return [os.path.join(HS, b), "--mode=edmd", "--experiment=speed_of_sound", "--headless", "--kbt1", "--seed-drift-order=drift-first",
            "--edmd-acc=0", "--particles=100", "--particles-boxes=50,50", "--height=10.0", "--particle-radius=0.5",
            "--wall-thickness=0.05", "--wall-thickness-vis=0.05", "--lengths=10.0000", "--wall-masses=50", "--repeats=1",
            "--seed=20261014", "--wall-hold-steps=2000", "--fixed-dt=0.4", "--target-oscillations=10", "--oscillation-safety=1.0",
            "--oscillation-min-steps=10000", "--oscillation-max-steps=400000000", "--speed-sound-log-stride=8",
            f"--speed-sound-run-dir={d}", "--speed-sound-exact-seed=123456789"], dict(os.environ, HD_KE_TRACE="1")

def et(b, d):
    return [os.path.join(HS, b), "--mode=edmd", "--experiment=energy_transfer", "--headless", "--quiet", "--edmd-acc=0",
            "--seed-drift-order=drift-first", f"--energy-transfer-summary={d}/summary.csv", f"--energy-transfer-trace={d}/tr.csv",
            "--trace-every=60", "--particles=100", "--particles-boxes=50,50", "--particle-radius=0.5", "--l0=10.000000", "--height=10",
            "--num-walls=1", "--wall-positions=10.000000", "--wall-mass-factors=1000000000", "--wall-thickness=0.05",
            "--wall-thickness-vis=0.05", "--eff-output=wall-ke", "--wall-hold-steps=12000", "--steps=12000", "--fixed-dt=0.4",
            "--kbt1", "--seed=9700"], dict(os.environ, HD_PISTON_EVENTS=f"{d}/ev.csv")

for tag, b in BINS:
    v = subprocess.run([os.path.join(HS, b), "--version"], capture_output=True, text=True).stdout.splitlines()
    print(f"{tag:14s} {b:36s} {v[0].strip()} | {v[1].strip()}")
    for mode, fn in (("sos", sos), ("et", et)):
        d = os.path.join(OUT, f"{mode}_{tag}"); os.makedirs(d)
        cmd, env = fn(b, d)
        r = subprocess.run(cmd, env=env, cwd=HS, stdout=open(f"{d}/stdout.log", "w"), stderr=open(f"{d}/stderr.log", "w"))
        if r.returncode != 0: sys.exit(f"{tag} {mode}: exit code {r.returncode} -- STOP")

FILES = {"sos": ["wall_x_positions_*_run0.csv", "speed_of_sound_psi6.csv"], "et": ["tr.csv", "ev.csv"]}
print("\n| mode | file | size [bytes] | new vs e823187 | new vs 5190846dirty | e823187 vs 5190846dirty |")
print("|---|---|---|---|---|---|")
allok = True
for mode, pats in FILES.items():
    for pat in pats:
        p = [glob.glob(os.path.join(OUT, f"{mode}_{t}", pat)) for t, _ in BINS]
        assert all(len(x) == 1 for x in p), (mode, pat, p)
        p = [x[0] for x in p]
        c = [filecmp.cmp(p[i], p[j], shallow=False) for i, j in ((0, 1), (0, 2), (1, 2))]
        allok &= all(c)
        print(f"| {mode} | `{os.path.basename(p[0])}` | {os.path.getsize(p[0])} | " + " | ".join("IDENTICAL" if x else "**DIFFERENT**" for x in c) + " |")
skip = {"build_git", "command", "timestamp", "trace_path"}
rows = [list(csv.DictReader(open(os.path.join(OUT, f"et_{t}", "summary.csv"))))[-1] for t, _ in BINS]
diff = sorted({k for k in rows[0] if k not in skip and len({r.get(k) for r in rows}) > 1})
print(f"\nET summary: {len(rows[0])} fields; fields differing outside {sorted(skip)}: {diff or 'none'}; "
      f"build_git = {[r['build_git'] for r in rows]}; box_width_sigma = {rows[0]['box_width_sigma']}")
hl = sum(open(f).read().count(k) for f in glob.glob(os.path.join(OUT, "*", "std*.log"))
         for k in ("EDMD-HEALTH", "forced_advance", "clamp_repair", "overlap_repair", "wall_overdue"))
print(f"health lines in all logs: {hl}")
allok &= not diff
print("\nCROSS-CHECK:", "IDENTICAL across all three binaries in both modes" if allok else "DIFFERENCES FOUND -- STOP (not investigated, per the task)")
sys.exit(0 if allok else 1)
