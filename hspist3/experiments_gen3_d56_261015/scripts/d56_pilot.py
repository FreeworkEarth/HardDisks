#!/usr/bin/env python3
"""##CHRIS 2026-10-15 (plan-author decision 56, task 2 c; notes 261012 sec. 4.7.133): THE STAGE 2 PILOT CELL ON THE MAC, and three
stage 1 seeds of (N 900, eta 0.69) again, to look at the divider's motion in the solid. Information for the sizing (Report D56-1);
nothing here is registered data.
  P69  the registered stage 1 task lines of B_N0900_eta0.690000, M 4500 (alpha 5, the middle of the nine masses), s = 0, 1, 2, copied
       verbatim from cluster/gen3_stage1_261012/tasks_stage1_koa_261012.txt: on this Mac (arm64) their trajectories need not be bit for
       bit KOA's (x86), but they are trajectories of the same protocol and seeds.
  P76  the pilot cell N 900, eta 0.76, M 4500 (alpha 5), s = 0, 1, built with the stage 1 generator's own rules
       (gen_stage1_koa_261012.plan): H = 10 sqrt(N/100) = 30, L0 = N pi / (8 H eta), T_eq = 5000 sigma-time (300000 hold steps),
       200 predicted periods of the binary's own nu_pred (Henderson), Delta_psi = max(1, T_rec / 2000), snapshots every
       (T_eq + T_rec) / 100 (s = 0, 1 of alpha 5), the default lattice start. Two departures, both forced by eta 0.76: the trace stride
       is d_stride of the binary's nu_pred (the generator's rule takes nu from KR, which is not defined beyond eta 0.7069), and the
       seeds come from a new base, run_seed(20261015, 0, 4, s), which no stage 1 or Part A list holds.
Each line runs through the stage 1 worker itself (cluster/gen3_stage1_261012/s1_worker_koa.sh, unchanged) with HD_BIN = the build of the
merged main (experiments_gen3_d50_261014/bin_0544b7e8-dirty: Part B's task lines EQUAL to the campaign path, Report D55a2) and HD_DATA =
experiments_gen3_d56_261015/data, five processes at once.
usage (from hspist3/): python3 experiments_gen3_d56_261015/scripts/d56_pilot.py [--dry-run]
"""
import math, os, subprocess, sys, time
HS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); E = os.path.join(HS, "experiments_gen3_d56_261015")
sys.path.insert(0, os.path.join(HS, "cluster", "gen3_stage1_261012")); sys.path.insert(0, os.path.join(HS, "validation")); sys.path.insert(0, HS)
import gen_stage1_koa_261012 as GEN
import tests_20260913 as T

BIN = os.path.join(HS, "experiments_gen3_d50_261014", "bin_0544b7e8-dirty", "00ALLINONE")
WORKER = os.path.join(HS, "cluster", "gen3_stage1_261012", "s1_worker_koa.sh")
BASE_P = 20261015


def pilot_lines(N=900, eta=0.76, a=5, seeds=(0, 1)):
    H = 10.0 * math.sqrt(N / 100.0); L0 = N * math.pi / (8.0 * H * eta); Ns = N // 2; M = int(round(a * N)); mi = GEN.ALPHA.index(a)
    Teq = 5000.0; hold = int(round(60 * Teq))
    nup = GEN.nu_pred_binary(N, H, L0, M); steps = math.ceil((200.0 / nup) / GEN.DTF); Trec = steps * GEN.DTF
    stride = T.d_stride(nup); dpsi = max(1.0, Trec / 2000.0); tsnap = (Teq + Trec) / 100.0
    rel = f"pilot/P2_N{N:04d}_eta{eta:.6f}/m_{M}"
    out = []
    for s in seeds:
        sd = T.run_seed(BASE_P, 0, mi, s)
        out.append(f"B {rel} {M} {s} {sd} {L0:.17g} {H:.17g} {Ns} {stride} {BASE_P} {hold} {dpsi:.17g} {tsnap:.17g} default B")
    info = dict(N=N, eta=eta, H=H, L0=L0, M=M, nup=nup, Trec=Trec, steps=steps, stride=stride, dpsi=dpsi, tsnap=tsnap)
    return out, info


def registered_lines(cell="B_N0900_eta0.690000", M=4500, seeds=(0, 1, 2)):
    want = {f"stage1/gen3/{cell}/m_{M}"}; out = []
    for l in open(GEN.TASKS):
        f = l.split()
        if f[1] in want and int(f[3]) in seeds: out.append(l.rstrip("\n"))
    assert len(out) == len(seeds), out
    return out


def main():
    p76, info = pilot_lines(); p69 = registered_lines()
    lines = p69 + p76
    os.makedirs(os.path.join(E, "data"), exist_ok=True)
    tf = os.path.join(E, "tasks_d56.txt")
    if not os.path.exists(tf):
        open(tf, "w").write("\n".join(lines) + "\n")
    elif open(tf).read() != "\n".join(lines) + "\n":
        sys.exit("REFUSED: tasks_d56.txt exists with other lines")
    print(f"pilot cell: {info}")
    print("\n".join(lines))
    if "--dry-run" in sys.argv: return
    env = dict(os.environ, HD_BIN=BIN, HD_DATA=os.path.join(E, "data"))
    procs = []; t0 = time.time()
    for k, l in enumerate(lines):
        lg = open(os.path.join(E, f"worker_{k}.out"), "w")
        procs.append((l, subprocess.Popen(["bash", WORKER] + l.split(), stdout=lg, stderr=subprocess.STDOUT, env=env, cwd=E)))
    for l, p in procs:
        rc = p.wait(); print(f"exit {rc}: {' '.join(l.split()[1:5])}  ({time.time() - t0:.0f} s)", flush=True)


if __name__ == "__main__":
    main()
