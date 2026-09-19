#!/usr/bin/env python3
"""##CHRIS 2026-09-15: A2 at N = 1600 with 200-period records (option A of 2026-09-15).

The A2 error bar is seed noise, and seed noise is set by the FFT bin width: 37.5-period records give
Δf/ν = 2.7 %, 200-period records give 0.5 %. Error falls as 1/(N_cyc sqrt(n)) at a cost of N_cyc·n, so a
longer record buys precision twice as cheaply as more seeds. This runs the three densities where the error
bar matters for the claim -- eta = 0.50, 0.60, 0.65 -- at N = 1600, everything else identical to famB:
same geometry (L0 = 3.926990816987241/eta * 4, H = 40), same five masses, 10 seeds, drift-first,
--edmd-acc=0. Record 200 periods exactly (--oscillation-safety=1.0), stride ~32 samples per period.

usage: A2_long200_20260915.py run|analyse"""
import os, sys, math, glob, subprocess, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import tests_20260913 as T

OUT = os.path.join(T.ROOT, "A2_long200_20260915")
A2 = os.path.join(T.ROOT, "famB_20260911")
N, H, MASSES, SEEDS, TARGET = 1600, "40.000000", (50, 200, 500, 1000, 2000), 10, 200
# ##CHRIS 2026-09-16: eta = 0.60 dropped on Chris's call after 2 of its 5 mass batches had run.
# It already agrees with KR at 37.5 periods (+0.18 %), so the long records buy nothing there, and
# dropping it takes ~9 h off the queue. Its completed cells (M = 2000) stay on disk and are analysed
# if present. The open question is eta = 0.65, where the 37.5-period result is -1.22 %.
ETAS = (0.50, 0.65)
ETAS_DONE_PARTIAL = (0.60,)
BASE = 23160915


def geom(eta):
    return f"{3.926990816987241 / eta * math.sqrt(N / 100):.6f}"


def x_of(M, L0):
    return T.k_root(M / float(N)) / (2 * math.pi * T.l_eff(L0))


def argv(eta, M, mi, run_dir, seed, stride):
    return ["./00ALLINONE", "--mode=edmd", "--experiment=speed_of_sound", "--headless", "--kbt1",
            "--seed-drift-order=drift-first", "--edmd-acc=0",
            f"--particles={N}", f"--particles-boxes={N//2},{N//2}",
            f"--height={H}", "--particle-radius=0.5",
            "--wall-thickness=0.05", "--wall-thickness-vis=0.05",
            f"--lengths={geom(eta)}", f"--wall-masses={M}", "--repeats=1",
            f"--seed={BASE + int(round(eta * 100))}", "--wall-hold-steps=2000", "--fixed-dt=0.4",
            f"--target-oscillations={TARGET}", "--oscillation-safety=1.0",
            "--oscillation-min-steps=10000", "--oscillation-max-steps=400000000",
            f"--speed-sound-log-stride={stride}", f"--speed-sound-run-dir={run_dir}",
            f"--speed-sound-exact-seed={seed}"]


def wait_ac():
    while "AC Power" not in subprocess.run(["pmset", "-g", "batt"], capture_output=True, text=True).stdout:
        print("on battery -- holding", flush=True); time.sleep(120)


def run():
    jobs = []
    for eta in ETAS:
        L0 = float(geom(eta))
        for mi, M in enumerate(MASSES):
            stride = T.d_stride(T.kr_cs(eta) * x_of(M, L0))
            for r in range(SEEDS):
                d = os.path.join(OUT, f"eta_{('%.2f' % eta).replace('.', 'p')}", f"N{N}", f"m_{M}", f"r{r}")
                # ##CHRIS 2026-09-16: "exists" is not "finished" -- a run killed mid-flight leaves a
                # partial trace, and skipping on existence alone would silently accept it. Require the
                # last sample to reach the planned duration recorded in the trace's own header.
                g = glob.glob(os.path.join(d, "wall_x_positions_*run0.csv"))
                if g and T.trace_check(g[0])[0]:
                    continue
                for stale in g:
                    os.remove(stale)
                jobs.append((eta, M, mi, r, d, stride, T.run_seed(BASE + int(round(eta * 100)), 0, mi, r)))
    # heaviest first: cost ~ 1/nu ~ 1/K(alpha)
    jobs.sort(key=lambda j: -j[1])
    print(f"{len(jobs)} runs to execute, 10 at a time, heaviest mass first", flush=True)
    t0 = time.time()
    while jobs:
        batch, jobs = jobs[:10], jobs[10:]
        wait_ac()
        procs = []
        for eta, M, mi, r, d, stride, seed in batch:
            os.makedirs(d, exist_ok=True)
            a = argv(eta, M, mi, d, seed, stride)
            open(os.path.join(d, "command.txt"), "w").write("HD_KE_TRACE=1 " + " ".join(a) + "\n")
            procs.append((subprocess.Popen(a, cwd=T.HSP, env=dict(os.environ, HD_KE_TRACE="1"),
                                           stdout=open(os.path.join(d, "stdout.log"), "w"), stderr=subprocess.STDOUT), d, time.time()))
        for p, d, ts in procs:
            p.wait()
            print(f"rc {p.returncode}  {(time.time()-ts)/60:5.1f} min  {os.path.relpath(d, OUT)}", flush=True)
        print(f"  batch done, {len(jobs)} left, {(time.time()-t0)/3600:.2f} h elapsed", flush=True)


def cell_nu(paths, cut_periods=None):
    import pandas as pd
    out = []
    for p in paths:
        d = pd.read_csv(p, usecols=["Time", "Displacement(σ)", "Predicted_Frequency"])
        t = d["Time"].to_numpy(float); x = d["Displacement(σ)"].to_numpy(float)
        nup = float(d["Predicted_Frequency"].iloc[0]); dt = (t[-1] - t[0]) / (len(t) - 1)
        if cut_periods:
            n = T._prefix(t, nup, cut_periods)
            if n is None:
                continue
            x = x[:n]
        P, df = T._spectrum(x, dt)
        ncyc = len(x) * dt * nup; k = max(1, int(round(ncyc / 2.5)))
        out.append((k + int(np.argmax(P[k:]))) * df)
    return out


def analyse():
    print("| η | record | seeds/mass | c_s | mass scatter [%] | error on mean [%] | KR | dev [%] | health-flagged |")
    print("|---|---|---|---|---|---|---|---|---|")
    for eta in ETAS:
        L0 = float(geom(eta)); kr = T.kr_cs(eta)
        for lab, root, cut in (("200 periods", OUT, None), ("37.5 periods (famB)", A2, None),
                               ("200 traces cut to 37.5", OUT, 37)):
            xs, ys, per, nseed, bad = [], [], [], [], 0
            for M in MASSES:
                if root is OUT:
                    ps = []
                    for d in sorted(glob.glob(os.path.join(OUT, f"eta_{('%.2f' % eta).replace('.', 'p')}", f"N{N}", f"m_{M}", "r*"))):
                        log = os.path.join(d, "stdout.log")
                        if os.path.exists(log) and "EDMD-HEALTH" in open(log, errors="replace").read():
                            bad += 1; continue
                        ps += glob.glob(os.path.join(d, "wall_x_positions_*run0.csv"))
                else:
                    cell = os.path.join(A2, f"eta_{('%.2f' % eta).replace('.', 'p')}", f"N{N}", f"m_{M}")
                    ps = [p for _, p, disc in T.cell_runs(cell, M) if not disc]
                nus = cell_nu(ps, cut)
                if len(nus) < 3:
                    continue
                a = np.array(nus); x = x_of(M, L0)
                xs.append(x); ys.append(a.mean()); per.append(a.mean() / x); nseed.append(len(a))
            if len(xs) < 3:
                print(f"| {eta:.2f} | {lab} | - | (no data yet) | | | {kr:.4f} | | |")
                continue
            xs, ys, per = np.array(xs), np.array(ys), np.array(per)
            cs = float((xs * ys).sum() / (xs * xs).sum()); sc = float(per.std(ddof=1))
            print(f"| {eta:.2f} | {lab} | {min(nseed)} | {cs:.4f} | {100*sc/cs:.2f} | {100*sc/cs/math.sqrt(len(per)):.2f} | "
                  f"{kr:.4f} | {100*(cs-kr)/kr:+.2f} | {bad} |")


if __name__ == "__main__":
    if sys.argv[1] == "run":
        run()
    analyse()
