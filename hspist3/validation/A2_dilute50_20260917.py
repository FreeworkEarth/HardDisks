#!/usr/bin/env python3
"""##CHRIS 2026-09-17: the dilute A2 ladders again with 50-PERIOD records (was 37.5). Same cells, same seeds, same everything else, separate root so the two record lengths never mix in one cell. The first 37.5 periods of every trace must reproduce the 2026-09-16 run byte for byte (same exact seed, same dynamics).

Below eta = 0.15 the size ladder exists at ONE density (0.10), so the larger systems are single points
and "the offset closes with N" rests on that one column. This adds two more dilute densities with the
full ladder N = 100, 400, 900, 1600, so each has its own 1/sqrt(N) extrapolation and the low-density
figure gets three real curves instead of one curve and two points.

Settings are famB's, so the new cells drop straight into the existing A2 analysis: L0 and H scale as
sqrt(N/100), five masses, 10 seeds, drift-first, --edmd-acc=0, 25 target oscillations x 1.5 safety
(37.5 periods), stride auto. Order: the whole eta = 0.05 ladder first (cheap, lands first), then 0.02.

usage: A2_dilute_20260916.py run|analyse"""
import os, sys, math, glob, subprocess, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import tests_20260913 as T

OUT = os.path.join(T.ROOT, "A2_dilute50_20260917")
MASSES = (50, 200, 500, 1000, 2000)
SEEDS = 10
CELLS = [(0.05, N) for N in (100, 400, 900, 1600)] + [(0.02, N) for N in (100, 400, 900, 1600)]
BASE = 24110916


def geom(eta, N):
    f = math.sqrt(N / 100.0)
    return f"{3.926990816987241 / eta * f:.6f}", f"{10 * f:.6f}"


def x_of(M, N, L0):
    return T.k_root(M / float(N)) / (2 * math.pi * T.l_eff(L0))


def base_for(eta, N):
    return BASE + 1000 * CELLS.index((eta, N))


def argv(eta, N, M, run_dir, seed):
    L0, H = geom(eta, N)
    return ["./00ALLINONE", "--mode=edmd", "--experiment=speed_of_sound", "--headless", "--kbt1",
            "--seed-drift-order=drift-first", "--edmd-acc=0",
            f"--particles={N}", f"--particles-boxes={N//2},{N//2}",
            f"--height={H}", "--particle-radius=0.5",
            "--wall-thickness=0.05", "--wall-thickness-vis=0.05",
            f"--lengths={L0}", f"--wall-masses={M}", "--repeats=1",
            f"--seed={base_for(eta, N)}", "--wall-hold-steps=2000", "--fixed-dt=0.4",
            "--target-oscillations=50", "--oscillation-safety=1.0",
            "--oscillation-min-steps=10000", "--oscillation-max-steps=400000000",
            "--speed-sound-log-stride=auto", f"--speed-sound-run-dir={run_dir}",
            f"--speed-sound-exact-seed={seed}"]


def cell_dir(eta, N, M, r):
    return os.path.join(OUT, f"eta_{('%.2f' % eta).replace('.', 'p')}", f"N{N}", f"m_{M}", f"r{r}")


def wait_ac():
    while "AC Power" not in subprocess.run(["pmset", "-g", "batt"], capture_output=True, text=True).stdout:
        print("on battery -- holding", flush=True); time.sleep(120)


def run():
    jobs = []
    for eta, N in CELLS:
        for mi, M in enumerate(MASSES):
            for r in range(SEEDS):
                d = cell_dir(eta, N, M, r)
                g = glob.glob(os.path.join(d, "wall_x_positions_*run0.csv"))
                if g and T.trace_check(g[0])[0]:
                    continue
                for stale in g:
                    os.remove(stale)
                jobs.append((eta, N, M, r, d, T.run_seed(base_for(eta, N), 0, mi, r)))
    print(f"{len(jobs)} runs to execute, 10 at a time, eta = 0.05 ladder first", flush=True)
    t0 = time.time()
    while jobs:
        batch, jobs = jobs[:10], jobs[10:]
        wait_ac()
        procs = []
        for eta, N, M, r, d, seed in batch:
            os.makedirs(d, exist_ok=True)
            a = argv(eta, N, M, d, seed)
            open(os.path.join(d, "command.txt"), "w").write("HD_KE_TRACE=1 " + " ".join(a) + "\n")
            procs.append((subprocess.Popen(a, cwd=T.HSP, env=dict(os.environ, HD_KE_TRACE="1"),
                                           stdout=open(os.path.join(d, "stdout.log"), "w"),
                                           stderr=subprocess.STDOUT), d, time.time()))
        for p, d, ts in procs:
            p.wait()
            print(f"rc {p.returncode}  {(time.time()-ts)/60:6.1f} min  {os.path.relpath(d, OUT)}", flush=True)
        print(f"  batch done, {len(jobs)} runs left, {(time.time()-t0)/3600:.2f} h elapsed", flush=True)


def analyse():
    import pandas as pd
    print("| η | N | L0 | trajectories | health-flagged | c_s | mass scatter [%] | error on mean [%] | KR | dev [%] |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for eta, N in CELLS:
        L0 = float(geom(eta, N)[0]); kr = T.kr_cs(eta)
        xs, ys, per, n_tot, bad = [], [], [], 0, 0
        for M in MASSES:
            nus = []
            for d in sorted(glob.glob(os.path.dirname(cell_dir(eta, N, M, 0)) + "/r*")):
                log = os.path.join(d, "stdout.log")
                if os.path.exists(log) and "EDMD-HEALTH" in open(log, errors="replace").read():
                    bad += 1; continue
                g = glob.glob(os.path.join(d, "wall_x_positions_*run0.csv"))
                if not (g and T.trace_check(g[0])[0]):
                    continue
                df = pd.read_csv(g[0], usecols=["Time", "Displacement(σ)", "Predicted_Frequency"])
                t = df["Time"].to_numpy(float); x = df["Displacement(σ)"].to_numpy(float)
                nup = float(df["Predicted_Frequency"].iloc[0]); dt = (t[-1] - t[0]) / (len(t) - 1)
                P, dfq = T._spectrum(x, dt)
                ncyc = len(x) * dt * nup; k = max(1, int(round(ncyc / 2.5)))
                nus.append((k + int(np.argmax(P[k:]))) * dfq)
            if len(nus) < 3:
                continue
            a = np.array(nus); xv = x_of(M, N, L0)
            xs.append(xv); ys.append(a.mean()); per.append(a.mean() / xv); n_tot += len(a)
        if len(xs) < 3:
            print(f"| {eta:.2f} | {N} | {L0:.2f} | {n_tot} | {bad} | (incomplete) | | | {kr:.4f} | |")
            continue
        xs, ys, per = np.array(xs), np.array(ys), np.array(per)
        cs = float((xs * ys).sum() / (xs * xs).sum()); sc = float(per.std(ddof=1))
        print(f"| {eta:.2f} | {N} | {L0:.2f} | {n_tot} | {bad} | {cs:.4f} | {100*sc/cs:.2f} | "
              f"{100*sc/cs/math.sqrt(len(per)):.2f} | {kr:.4f} | {100*(cs-kr)/kr:+.2f} |")


if __name__ == "__main__":
    if sys.argv[1] == "run":
        run()
    analyse()
