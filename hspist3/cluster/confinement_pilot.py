#!/usr/bin/env python3
"""##CHRIS 2026-10-13: the pi/8 pilot of the Paper 1 confinement campaign, method B (261012 sec. 1.4), and the
determinism self-test. ONE script runs and analyses on both machines, so the Mac target and the KOA result go
through identical code.

Pilot cell (pre-registered anchor, 261012 sec. 1.2): eta = pi/8, H = L_0 = 10, N_s = 50, r = 0.5, t = 0.05;
the nine A1v2 masses (alpha = 0.5 ... 20), ONE seed each, 200 oscillations, A1v2 protocol (methods sec. 8).
Exact seed per mass: tests_20260913.run_seed(BASE, 0, mass_index, 0), BASE = 20261013.
c_s: canonical estimator (paper1_populate_cs_err_20261002.cell, TD = 200, X_EDGE = 2.5) per mass, then the
through-origin slope nu = c_s x; sigma = the 1-sigma scatter of the per-mass implied c_s (methods sec. 6 error
bar; the 2026-09-16 mirror-gate rule). With one seed per mass the per-mass seed error does not exist.

usage (from hspist3/):
  python3 cluster/confinement_pilot.py run         --bin ./00ALLINONE --out <dir> [--jobs 9]
  python3 cluster/confinement_pilot.py analyse     --out <dir>
  python3 cluster/confinement_pilot.py determinism --bin ./00ALLINONE --out <dir>
  python3 cluster/confinement_pilot.py det1        --bin ./00ALLINONE --out <dir> --tag A    (one run; KOA: one srun step)
  python3 cluster/confinement_pilot.py detcmp      --out <dir>                               (cmp det_A vs det_B)

##CHRIS 2026-10-02 (Task K2): det1/detcmp split the determinism self-test so that KOA can run the two trajectories as two
separate Slurm steps (koa_smoketest.sh: same node; koa_crossnode_det.sh: two different nodes). Every trajectory now runs
with its own output directory as working directory (cwd=d), so the files the binary writes into its cwd (run_params.json,
00_COMMAND.md; 00ALLINONE.c:1819, 2837) land next to the trace instead of in the git checkout -- a sparse checkout turns
-dirty if a tracked path appears in it. --bin and --out are made absolute first. Trace bytes do not depend on the cwd.
"""
import argparse, glob, math, os, subprocess, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(HS, "validation")); sys.path.insert(0, HS)
import tests_20260913 as T

ETA, L0, H, NS, R, TW = 0.392699, 10.0, 10.0, 50, 0.5, 0.05
BASE, TARGET = 20261013, 200

def cmd(binp, M, mi, out, target=TARGET):
    x = T.x_of(M, L0); nu = T.kr_cs(ETA) * x; seed = T.run_seed(BASE, 0, mi, 0)
    return [binp, "--mode=edmd", "--experiment=speed_of_sound", "--headless", "--kbt1", "--seed-drift-order=drift-first",
            "--edmd-acc=0", f"--particles={2*NS}", f"--particles-boxes={NS},{NS}", f"--height={H}", f"--particle-radius={R}",
            f"--wall-thickness={TW}", f"--wall-thickness-vis={TW}", f"--lengths={L0:.4f}", f"--wall-masses={M}", "--repeats=1",
            f"--seed={BASE}", "--wall-hold-steps=2000", "--fixed-dt=0.4", f"--target-oscillations={target}",
            "--oscillation-safety=1.0", "--oscillation-min-steps=10000", "--oscillation-max-steps=400000000",
            f"--speed-sound-log-stride={T.d_stride(nu)}", f"--speed-sound-run-dir={out}", f"--speed-sound-exact-seed={seed}"]

def run(a):
    a.bin, a.out = os.path.abspath(a.bin), os.path.abspath(a.out)
    procs = []
    for mi, M in enumerate(T.A1_MASSES):
        d = os.path.join(a.out, f"m_{M}"); os.makedirs(d, exist_ok=True)
        c = cmd(a.bin, M, mi, d); open(os.path.join(d, "command.txt"), "w").write("HD_KE_TRACE=1 " + " ".join(c) + "\n")
        while len([p for p in procs if p.poll() is None]) >= a.jobs: time.sleep(0.5)
        procs.append(subprocess.Popen(c, stdout=open(os.path.join(d, "run.log"), "w"), stderr=subprocess.STDOUT,
                                      env=dict(os.environ, HD_KE_TRACE="1"), cwd=d))
    rc = [p.wait() for p in procs]; print("run exit codes:", rc); return 0 if all(r == 0 for r in rc) else 1

def analyse(a):
    from paper1_populate_cs_err_20261002 import cell
    xs, nus, rows, Ttot, etas, L0s, health = [], [], [], 0.0, set(), set(), 0
    for M in T.A1_MASSES:
        d = os.path.join(a.out, f"m_{M}"); runs = T.cell_runs(d, M)
        health += sum(1 for l in open(os.path.join(d, "run.log"), errors="ignore") if any(k in l for k in
                      ("EDMD-HEALTH", "forced_advance", "clamp_repair", "overlap_repair", "wall_overdue")))
        c = cell((ETA, L0, M, runs)); x = T.x_of(M, L0)
        tr = sorted(glob.glob(os.path.join(d, "wall_x_positions_*_run0.csv")))[0]
        import pandas as pd
        h = pd.read_csv(tr, nrows=2); etas.add(round(float(h["eta"].iloc[0]), 6)); L0s.add(float(h["L0"].iloc[0]))
        Ttot += float(h["Planned_Duration"].iloc[0]); xs.append(x); nus.append(c["nu"]); rows.append((M, c["nu"], c["nu"] / x, c["n"]))
    xs, nus = np.array(xs), np.array(nus); cs = float(np.sum(xs * nus) / np.sum(xs * xs)); imp = nus / xs
    print(f"pilot cell: eta (trace) = {sorted(etas)}, L_0 (trace) = {sorted(L0s)}, H = {H}, N_s = {NS}, r = {R}, "
          f"t = {TW} (set by --wall-thickness; not written by speed-of-sound mode), L_eff = L_0 - 2r - t/2 = {T.l_eff(L0, TW):.6f}")
    print(f"T_total (sum of planned durations, 9 trajectories) = {Ttot:.1f} sigma-time; health lines = {health}")
    print("| M | nu | implied c_s | seeds |\n|---|---|---|---|")
    for M, nu, c, n in rows: print(f"| {M} | {nu:.8f} | {c:.5f} | {n} |")
    print(f"\n**c_s = {cs:.5f} +- {np.std(imp, ddof=1):.5f}** (through-origin slope; +- = 1-sigma mass scatter of implied c_s)")
    return 0

def determinism(a):
    a.bin, a.out = os.path.abspath(a.bin), os.path.abspath(a.out)
    outs = []
    for tag in ("A", "B"):
        d = os.path.join(a.out, f"det_{tag}", "m_50"); os.makedirs(d, exist_ok=True)
        subprocess.run(cmd(a.bin, 50, 0, d, target=25), stdout=open(os.path.join(d, "run.log"), "w"), stderr=subprocess.STDOUT,
                       env=dict(os.environ, HD_KE_TRACE="1"), check=True, cwd=d)
        outs.append(sorted(glob.glob(os.path.join(d, "wall_x_positions_*_run0.csv")))[0])
    same = subprocess.run(["cmp", outs[0], outs[1]]).returncode == 0
    print(f"determinism self-test (same binary, same seed, twice): {'IDENTICAL' if same else 'DIFFERENT -- STOP'}")
    return 0 if same else 1

def det1(a):
    """One determinism trajectory (M = 50, 25 oscillations, the same seed as `determinism`) into <out>/det_<tag>/m_50."""
    a.bin, a.out = os.path.abspath(a.bin), os.path.abspath(a.out)
    d = os.path.join(a.out, f"det_{a.tag}", "m_50")
    if glob.glob(os.path.join(d, "wall_x_positions_*_run0.csv")): sys.exit(f"{d} already holds a trace -- not overwriting")
    os.makedirs(d, exist_ok=True)
    import socket; open(os.path.join(d, "host.txt"), "w").write(socket.gethostname() + "\n")
    r = subprocess.run(cmd(a.bin, 50, 0, d, target=25), stdout=open(os.path.join(d, "run.log"), "w"), stderr=subprocess.STDOUT,
                       env=dict(os.environ, HD_KE_TRACE="1"), cwd=d)
    print(f"det1 {a.tag}: host {socket.gethostname()}, exit code {r.returncode}"); return r.returncode

def detcmp(a):
    """cmp the trace and the psi6 file of det_A and det_B; print the hosts they ran on."""
    a.out = os.path.abspath(a.out); ok = True
    hosts = [open(os.path.join(a.out, f"det_{t}", "m_50", "host.txt")).read().strip() for t in "AB"]
    for pat in ("wall_x_positions_*_run0.csv", "speed_of_sound_psi6.csv"):
        p = [sorted(glob.glob(os.path.join(a.out, f"det_{t}", "m_50", pat))) for t in "AB"]
        if not all(len(x) == 1 for x in p): print(f"{pat}: missing in one run -- DIFFERENT -- STOP"); ok = False; continue
        same = subprocess.run(["cmp", p[0][0], p[1][0]]).returncode == 0; ok &= same
        print(f"{os.path.basename(p[0][0])}: {'IDENTICAL' if same else 'DIFFERENT -- STOP'} ({os.path.getsize(p[0][0])} bytes)")
    print(f"determinism self-test (same binary, same seed, run A on {hosts[0]}, run B on {hosts[1]}"
          f"{', different nodes' if hosts[0] != hosts[1] else ', same node'}): {'IDENTICAL' if ok else 'DIFFERENT -- STOP'}")
    return 0 if ok else 1

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("what", choices=["run", "analyse", "determinism", "det1", "detcmp"])
    ap.add_argument("--tag", choices=["A", "B"], default="A")
    ap.add_argument("--bin", default=os.path.join(HS, "00ALLINONE")); ap.add_argument("--out", required=True)
    ap.add_argument("--jobs", type=int, default=9); a = ap.parse_args()
    sys.exit({"run": run, "analyse": analyse, "determinism": determinism, "det1": det1, "detcmp": detcmp}[a.what](a))
