#!/usr/bin/env python3
"""##CHRIS 2026-10-03 (Task W3): Round 1 wall times on KOA -> cost per trajectory vs N_s -> --time for the resubmission.

DATA (sacct from Chris's terminal, 2026-10-03; build 279282b target koa): Elapsed per array task and State, below.
Trajectories per cell and parallelism P from the task files and the sbatch (#SBATCH --cpus-per-task = xargs -P).
A cell's n trajectories run in waves of P, so the measured wall per wave is Elapsed / ceil(n/P) [INFERENCE: equal
trajectories per wave; for method B the nine masses differ in run length, so this is a cell average].

Fit (as asked, Task W3): log(wall per wave) vs log(N_s) over the completed H-scan cells H5, H10, H20 (N_s 25, 50, 100;
same L_0) of each array -> exponent p, prediction for H40 (N_s 200). Check: a TIMEOUT cell needs MORE than its Elapsed,
so a prediction below the timeout Elapsed is falsified by the data. Printed per array: the fit exponent, the local
exponent of the last doubling (H10 -> H20), and the lower bound on the H20 -> H40 exponent that the timeout implies.
USED for every prediction: p* = the steepest exponent measured in Round 1 (the largest local exponent over the three
arrays), and p* must exceed every timeout lower bound [the script checks it]. Prediction for a cell at N_s from the
nearest completed cell below it (same array, same L_0 scan): wall per wave x (N_s/N_s,ref)^p* x waves.
--time = 3 x prediction of the WHOLE cell, rounded up to 15 min, capped at 3-00:00:00 (shared); a cell whose 3 x
prediction is below the sbatch default keeps the default. The remaining work of a TIMEOUT cell is less than its whole,
so no count of finished trajectories is needed; cluster/check_cells.sh measures that count.
conf_A_0.39 (no KOA data at pi/8 beyond the pilot): per wave at N_s 50 = pilot Elapsed 41 s / ceil(20/16) = 2 waves,
scaled with p* [INFERENCE; at pi/8, B_0.39 went x8 from N_s 50 to 100].
Throttles: the 64-core cap holds for each set of lines that may run together (printed).
usage: python3 hspist3/cluster/round1_timing_261003.py
"""
import math, os, re
HERE = os.path.dirname(os.path.abspath(__file__)); CONF = os.path.join(HERE, "confinement_20261013")
CAP_H = 72.0
SACCT = {   # array -> [(task, Elapsed, State)]   DATA, Chris's terminal 2026-10-03
    "B_0.10": [(1, "00:06:31", "C"), (2, "00:15:42", "C"), (3, "00:46:07", "C"), (4, "02:47:01", "T"), (5, "00:03:12", "C"),
               (6, "01:22:51", "C"), (7, "00:09:08", "C"), (8, "00:11:44", "C"), (9, "00:15:32", "C"), (10, "00:20:27", "C")],
    "B_0.39": [(1, "00:04:08", "C"), (2, "00:08:14", "C"), (3, "01:05:33", "C"), (4, "02:02:00", "T"), (5, "00:02:30", "C"),
               (6, "01:33:20", "C"), (7, "00:10:19", "C"), (8, "00:13:17", "C"), (9, "00:17:56", "C")],
    "A_0.10": [(1, "00:01:56", "C"), (2, "00:03:03", "C"), (3, "00:06:42", "C"), (4, "00:32:06", "T"), (5, "00:02:04", "C"),
               (6, "00:06:29", "C"), (7, "00:03:46", "C"), (8, "00:03:19", "C"), (9, "00:03:10", "C"), (10, "00:03:07", "C")],
}
PILOT_ELAPSED_S, PILOT_N, PILOT_P = 41.0, 20, 16           # job 14966594, N_s = 50 at pi/8

def sec(t): h, m, s = map(int, t.split(":")); return 3600 * h + 60 * m + s
def hms(h): m = math.ceil(h * 4) * 15; return f"{m // 60:02d}:{m % 60:02d}:00" if m < 1440 else f"{m // 1440}-{(m % 1440) // 60:02d}:{m % 60:02d}:00"
def sbatch(g):
    s = open(os.path.join(CONF, f"conf_{g}.sbatch")).read()
    P = int(re.search(r"--cpus-per-task=(\d+)", s).group(1)); t = re.search(r"--time=(\S+)", s).group(1).split(":")
    return P, int(t[0]) + int(t[1]) / 60
def cells(g): return open(os.path.join(CONF, f"cells_{g}.tsv")).read().split()
def ntraj(g, c): return sum(1 for _ in open(os.path.join(CONF, f"tasks_{g[0]}_{c}.txt")))
def ns(g, c): return int(open(os.path.join(CONF, f"tasks_{g[0]}_{c}.txt")).readline().split()[7 if g[0] == "B" else 6])

def fit(xs, ys):
    lx, ly = [math.log(x) for x in xs], [math.log(y) for y in ys]; mx, my = sum(lx) / len(lx), sum(ly) / len(ly)
    p = sum((a - mx) * (b - my) for a, b in zip(lx, ly)) / sum((a - mx) ** 2 for a in lx); return p, math.exp(my - p * mx)

def main():
    print("### Round 1, measured wall per wave (sacct Elapsed / ceil(n/P))\n")
    print("| array | task | cell | N_s | trajectories | P | waves | Elapsed | state | wall per wave (s) |\n|---|---|---|---|---|---|---|---|---|---|")
    per = {}
    for g, rows in SACCT.items():
        P, _ = sbatch(g); cl = cells(g)
        for task, el, st in rows:
            c = cl[task - 1]; n = ntraj(g, c); w = math.ceil(n / P)
            per[(g, c)] = dict(task=task, Ns=ns(g, c), n=n, P=P, waves=w, el=sec(el), st=st, pw=sec(el) / w)
            print(f"| {g} | {task} | {c} | {ns(g, c)} | {n} | {P} | {w} | {el} | {'TIMEOUT' if st == 'T' else 'COMPLETED'} | "
                  f"{'>' if st == 'T' else ''}{sec(el) / w:.1f} |")
    print("\n### Fit over the H-scan cells H5, H10, H20 (N_s 25, 50, 100) -> H40 (N_s 200)\n")
    print("| array | exponent p (fit) | local exponent H10->H20 | H40 predicted, fit (h) | H40 predicted, local (h) | H40 TIMEOUT Elapsed (h) | fit consistent with TIMEOUT? |\n|---|---|---|---|---|---|---|")
    pred, expo = {}, {}
    for g in SACCT:
        hs = [c for c in cells(g) if "_H_H" in c]; done = [per[(g, c)] for c in hs if per[(g, c)]["st"] == "C"]
        p, a = fit([d["Ns"] for d in done], [d["pw"] for d in done])
        loc = math.log(done[-1]["pw"] / done[-2]["pw"]) / math.log(done[-1]["Ns"] / done[-2]["Ns"])
        h40 = [per[(g, c)] for c in hs if per[(g, c)]["st"] == "T"][0]
        t_fit = a * h40["Ns"] ** p * h40["waves"] / 3600; t_loc = done[-1]["pw"] * (h40["Ns"] / done[-1]["Ns"]) ** loc * h40["waves"] / 3600
        pred[g] = (h40, max(t_fit, t_loc)); expo[g] = (p, loc)
        print(f"| {g} | {p:.2f} | {loc:.2f} | {t_fit:.2f} | {t_loc:.2f} | > {h40['el'] / 3600:.2f} | "
              f"{'yes' if t_fit > h40['el'] / 3600 else '**NO** (fit below the timeout)'} |")
    pstar = max(l for _, l in expo.values())
    lb = {g: math.log(pred[g][0]["pw"] / per[(g, [c for c in cells(g) if "_H_H20" in c][0])]["pw"]) / math.log(2) for g in SACCT}
    print(f"\nlower bound on the H20 -> H40 exponent from the timeouts: " + ", ".join(f"{g} > {v:.2f}" for g, v in lb.items()))
    print(f"USED exponent p* = steepest local exponent measured = {pstar:.2f}; exceeds every timeout bound: "
          f"{'yes' if all(pstar > v for v in lb.values()) else '**NO**'}")
    print("\n### --time for the resubmission (p* scaling, 3 x the whole cell)\n")
    print("| array | task(s) | cell | predicted whole cell (h) | basis | --time |\n|---|---|---|---|---|---|")
    lines, A39 = [], []
    for g in SACCT:
        h40, _ = pred[g]; P, tdef = sbatch(g); c = cells(g)[h40["task"] - 1]
        ref = per[(g, [x for x in cells(g) if "_H_H20" in x][0])]
        tp = ref["pw"] * (h40["Ns"] / ref["Ns"]) ** pstar * h40["waves"] / 3600; tt = min(CAP_H, max(3 * tp, tdef))
        print(f"| {g} | {h40['task']} | {c} | {tp:.2f} | H20 {ref['pw']:.1f} s/wave x 2^{pstar:.2f} x {h40['waves']} waves (> timeout {h40['el'] / 3600:.2f} h: {'yes' if tp > h40['el'] / 3600 else '**NO**'}) | {hms(tt)} |")
        rest = [str(t) for t, _, _ in SACCT[g] if t != h40["task"]]
        print(f"| {g} | {','.join(rest)} | the others (COMPLETED; only missing seeds run) | -- | default | {hms(tdef)} (sbatch) |")
        lines.append((1, f"sbatch --array={h40['task']} --time={hms(tt)} cluster/confinement_20261013/conf_{g}.sbatch", P))
        lines.append((1, f"sbatch --array=1-{h40['task'] - 1},{h40['task'] + 1}-{len(cells(g))}%1 cluster/confinement_20261013/conf_{g}.sbatch", P))
    g = "A_0.39"; P, tdef = sbatch(g); pw50 = PILOT_ELAPSED_S / math.ceil(PILOT_N / PILOT_P)
    groups = {}
    for i, c in enumerate(cells(g), 1):
        n, N = ntraj(g, c), ns(g, c); t = pw50 * max(1.0, N / 50) ** pstar * math.ceil(n / P) / 3600
        tt = hms(min(CAP_H, 3 * t)) if 3 * t > tdef else None
        groups.setdefault(tt, []).append(i)
        print(f"| {g} | {i} | {c} (N_s {N}, {n} traj.) | {t:.2f} | pilot {pw50:.1f} s/wave x (N_s/50)^{pstar:.2f} x {math.ceil(n / P)} waves | "
              f"{tt if tt else hms(tdef) + ' (sbatch default >= 3 x)'} |")
    for tt, idx in groups.items():
        thr = "" if len(idx) == 1 else "%1"
        lines.append((2, f"sbatch --array={','.join(map(str, idx))}{thr}{' --time=' + tt if tt else ''} "
                         f"--dependency=afterok:<A_0.10 task-4 jobid>:<A_0.10 rest jobid> cluster/confinement_20261013/conf_{g}.sbatch", P))
    print("\n### Resubmission lines (from ~/harddisks/hspist3 on login-0102; the skip logic runs only what is missing)\n")
    for r in (1, 2):
        cores = sum(P * (1 if "%" not in l else int(l.split("%")[1].split()[0])) for rr, l, P in lines if rr == r)
        print(f"Set {r} ({'now' if r == 1 else 'Round 2, behind the A_0.10 resubmission'}); at most {cores} cores at once:")
        for rr, l, _ in lines:
            if rr == r: print("    " + l)
    long_b = sum(P for rr, l, P in lines if rr == 1 and "--time" in l and "B_" in l)
    c2 = long_b + sum(P * (1 if "%" not in l else int(l.split("%")[1].split()[0])) for rr, l, P in lines if rr == 2)
    print(f"Set 2 starts while the two long B H40 tasks of set 1 may still run: {long_b} + set 2 = {c2} cores (cap {64}).")


if __name__ == "__main__":
    main()
