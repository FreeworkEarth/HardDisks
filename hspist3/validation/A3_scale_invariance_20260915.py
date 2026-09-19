#!/usr/bin/env python3
"""##CHRIS 2026-09-15: A3, scale-invariance check. Hard disks have no length scale but the disk, so
"same box, smaller disks" must reproduce "same disks, bigger box" once both are read in disk diameters.

A2 N = 400 at eta = 0.30 is L0 = 26.18, H = 20, r = 0.5 (units of the disk diameter).
A3 runs the same physics with the BOX of the A2 N = 100 cell (L0 = 13.09, H = 10 in code length units) and
smaller disks: N = 400 with r = 0.25 and N = 900 with r = 1/6. Everything measured in disk diameters is then
identical to A2 N = 400 and N = 900 -- including the divider, whose thickness is scaled with the disk
(0.05 -> 0.025 and 0.0166667) so it stays 0.05 diameters. What differs is only the absolute/pixel picture:
a disk is 12 px (N = 400) or 8 px (N = 900) across instead of 24 px. Agreement within the mass scatter means
no absolute-unit artefact (pixel grid, seed pad, wall thickness in pixels); a difference beyond it is reported,
not fixed.

usage: A3_scale_invariance_20260915.py run|analyse"""
import os, sys, math, subprocess, glob
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import tests_20260913 as T

OUT = os.path.join(T.ROOT, "A3_scale_invariance_20260915")
ETA, L0, H = 0.30, "13.090000", "10.000000"
MASSES = (50, 200, 500, 1000, 2000)
SEEDS = 10
CELLS = [dict(N=400, r="0.25", wall="0.025", base=23060915, tag="N400_r0p25"),
         dict(N=900, r="0.166667", wall="0.016667", base=23060916, tag="N900_r0p166667")]
A2ROOT = os.path.join(T.ROOT, "famB_20260911", "eta_0p30")


def argv(c, M, r_idx, run_dir, seed):
    return ["./00ALLINONE", "--mode=edmd", "--experiment=speed_of_sound", "--headless", "--kbt1",
            "--seed-drift-order=drift-first", "--edmd-acc=0",
            f"--particles={c['N']}", f"--particles-boxes={c['N']//2},{c['N']//2}",
            f"--height={H}", f"--particle-radius={c['r']}",
            f"--wall-thickness={c['wall']}", f"--wall-thickness-vis={c['wall']}",
            f"--lengths={L0}", f"--wall-masses={M}", "--repeats=1", f"--seed={c['base']}",
            "--wall-hold-steps=2000", "--fixed-dt=0.4", "--target-oscillations=25", "--oscillation-safety=1.5",
            "--oscillation-min-steps=10000", "--oscillation-max-steps=80000000",
            "--speed-sound-log-stride=auto", f"--speed-sound-run-dir={run_dir}",
            f"--speed-sound-exact-seed={seed}"]


def wait_ac():
    while True:
        out = subprocess.run(["pmset", "-g", "batt"], capture_output=True, text=True).stdout
        if "AC Power" in out:
            return
        print("on battery -- holding", flush=True); __import__("time").sleep(120)


def run():
    jobs = []
    for c in CELLS:
        for mi, M in enumerate(MASSES):
            for r in range(SEEDS):
                d = os.path.join(OUT, c["tag"], f"m_{M}", f"r{r}")
                if glob.glob(os.path.join(d, "wall_x_positions_*run0.csv")):
                    continue
                jobs.append((c, M, mi, r, d))
    print(f"{len(jobs)} runs to execute, 10 at a time", flush=True)
    while jobs:
        batch, jobs = jobs[:10], jobs[10:]
        procs = []
        wait_ac()
        for c, M, mi, r, d in batch:
            os.makedirs(d, exist_ok=True)
            seed = T.run_seed(c["base"], 0, mi, r)
            a = argv(c, M, mi, d, seed)
            open(os.path.join(d, "command.txt"), "w").write("HD_KE_TRACE=1 " + " ".join(a) + "\n")
            procs.append((subprocess.Popen(a, cwd=T.HSP, env=dict(os.environ, HD_KE_TRACE="1"),
                                           stdout=open(os.path.join(d, "stdout.log"), "w"), stderr=subprocess.STDOUT), d))
        for p, d in procs:
            p.wait()
            print("rc", p.returncode, os.path.relpath(d, OUT), flush=True)


def cs_from(cells, n_side_for_alpha):
    xs, ys, rows = [], [], []
    for M, nus, L0v, rad in cells:
        if len(nus) < 3:
            continue
        a = np.array(nus)
        x = T.k_root(M / (2.0 * n_side_for_alpha)) / (2 * math.pi * (L0v - 2 * rad))
        xs.append(x); ys.append(a.mean())
        rows.append((M, len(a), a.mean(), a.std(ddof=1), a.mean() / x))
    xs, ys = np.array(xs), np.array(ys)
    return float((xs * ys).sum() / (xs * xs).sum()), float(np.std(ys / xs, ddof=1)), rows


def analyse():
    import pandas as pd
    print("### A3 (fixed box, smaller disks) against A2 (fixed disks, bigger box), η = 0.30\n")
    for c in CELLS:
        a3 = []
        health = 0
        for M in MASSES:
            nus = []
            for d in sorted(glob.glob(os.path.join(OUT, c["tag"], f"m_{M}", "r*"))):
                log = os.path.join(d, "stdout.log")
                if os.path.exists(log) and "EDMD-HEALTH" in open(log, errors="replace").read():
                    health += 1; continue
                g = glob.glob(os.path.join(d, "wall_x_positions_*run0.csv"))
                if not g:
                    continue
                df = pd.read_csv(g[0], usecols=["Time", "Displacement(σ)", "Predicted_Frequency"])
                t = df["Time"].to_numpy(float); x = df["Displacement(σ)"].to_numpy(float)
                nup = float(df["Predicted_Frequency"].iloc[0]); dt = (t[-1] - t[0]) / (len(t) - 1)
                P, dfq = T._spectrum(x, dt)
                ncyc = len(x) * dt * nup; k = max(1, int(round(ncyc / 2.5)))
                nus.append((k + int(np.argmax(P[k:]))) * dfq)
            a3.append((M, nus, float(L0), float(c["r"])))
        # matching A2 cell in reduced units
        a2cell = os.path.join(A2ROOT, f"N{c['N']}")
        a2 = []
        for M in MASSES:
            nus = []
            for r, p, disc in T.cell_runs(os.path.join(a2cell, f"m_{M}"), M):
                if disc:
                    continue
                df = pd.read_csv(p, usecols=["Time", "Displacement(σ)", "Predicted_Frequency", "L0"])
                t = df["Time"].to_numpy(float); x = df["Displacement(σ)"].to_numpy(float)
                nup = float(df["Predicted_Frequency"].iloc[0]); L0v = float(df["L0"].iloc[0])
                dt = (t[-1] - t[0]) / (len(t) - 1)
                P, dfq = T._spectrum(x, dt)
                ncyc = len(x) * dt * nup; k = max(1, int(round(ncyc / 2.5)))
                nus.append((k + int(np.argmax(P[k:]))) * dfq)
            a2.append((M, nus, L0v, 0.5))
        c3, s3, rows3 = cs_from(a3, c["N"] / 2)
        c2, s2, rows2 = cs_from(a2, c["N"] / 2)
        print(f"**{c['tag']}** (r = {c['r']}, wall = {c['wall']}, box {L0} × {H}) vs A2 N = {c['N']} (r = 0.5)\n")
        print("| M | A3 runs | A3 ν̄ | A3 sd | A2 runs | A2 ν̄ | A2 sd | (A3 − A2)/A2 ν̄ [%] |")
        print("|---|---|---|---|---|---|---|---|")
        for (M, n3, m3, d3, _), (M2, n2, m2, d2, _) in zip(rows3, rows2):
            print(f"| {M} | {n3} | {m3:.8g} | {d3:.3g} | {n2} | {m2:.8g} | {d2:.3g} | {100*(m3-m2)/m2:+.3f} |")
        print(f"\nc_s: A3 {c3:.4f} ± {s3:.4f} (mass scatter) · A2 {c2:.4f} ± {s2:.4f} · difference {100*(c3-c2)/c2:+.3f} %"
              f" · combined scatter {100*math.hypot(s3, s2)/c2:.3f} % · health-flagged A3 runs {health}\n")


if __name__ == "__main__":
    if sys.argv[1] == "run":
        run()
    analyse()
