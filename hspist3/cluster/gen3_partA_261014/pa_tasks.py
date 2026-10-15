#!/usr/bin/env python3
"""##CHRIS 2026-10-14 (plan-author decision 50 item 3 d; 261012 sec. 4.7.113; registration sec. 4.8.4): PART A TASK LINES, option (ii) --
A DRAFT for amendment A5: the pilot's sizes and seeds, BASE_A, the check subset and the full run's seed counts are the plan author's.
What is registered and used here:
  - the cells: block B's cells of stage 1, with L0, H and N_s exactly as written in cluster/gen3_stage1_261012/tasks_stage1_koa_261012.txt
    (so that the held and the free geometry are the same doubles, gate C2);
  - the stencil x_j = L0 + j dL, j = -2..2, dL = sigma_x / 2 exactly, sigma_x = L_eff / sqrt(2 N_s c_KR^2), L_eff = L0 - 2r - t/2, c_KR^2
    at eta = N_s pi r^2 / (H L0) (the pre-registration script's functions), x_0 = L0 as written;
  - T_eq = 1000 sigma-time (60,000 steps) for eta <= 0.60 and 5000 sigma-time (300,000 steps) for eta >= 0.63 (sec. 4.8.3.1), then the
    window of 5000 sigma-time (300,000 steps), then the released tail of 20 sigma-time (1,200 steps), dt = 1/60 sigma-time;
  - seeds: tests_20260913.run_seed(BASE_A, stream, position index, s), the driver's own function as in stage 1: one stream per cell in
    block B's order, the position index 0..4 for j = -2..2. A seed that stage 1 uses is reported (a collision).
Task line (pa_worker_koa.sh): A <rel position folder> <j> <s> <seed> <L0> <H> <N_s> <x_div> <T_eq steps> <window steps> <tail steps>
<trace every> <log 0|1>; the position folder is partA/gen3/<cell>/x_<m2|m1|0|p1|p2> (the confinement worker's names).
usage (from hspist3/): python3 cluster/gen3_partA_261014/pa_tasks.py --sizes 100,400 --seeds 4 [--etas 0.05,...] [--base N]
       [--log s0j0|none] [--steps TEQ,WIN,TAIL,EVERY (a rehearsal's short durations)] > tasks.txt
"""
import argparse, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(HS, "validation")); sys.path.insert(0, HS)
import importlib.util
spec = importlib.util.spec_from_file_location("prereg", os.path.join(HS, "validation", "paper1_confinement_prereg_20261012.py"))
P = importlib.util.module_from_spec(spec); spec.loader.exec_module(P)
import tests_20260913 as T

S1 = os.path.join(HS, "cluster", "gen3_stage1_261012", "tasks_stage1_koa_261012.txt")
POS = (("m2", -2), ("m1", -1), ("0", 0), ("p1", 1), ("p2", 2))
STEPS_PER_SIGMA = 60                                     # dt = 0.4 / 24 sigma-time


def cells():
    """block B's cells in the stage 1 task list's order: (cell, N, eta label, L0 text, H text, N_s)"""
    seen, out = set(), []
    for l in open(S1):
        f = l.split()
        cell = f[1].split("/")[-2]
        if not cell.startswith("B_N") or cell in seen:
            continue
        seen.add(cell)
        out.append((cell, 2 * int(f[7]), cell.split("_eta")[1], f[5], f[6], int(f[7])))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sizes", default="100,400"); ap.add_argument("--seeds", type=int, default=4)
    ap.add_argument("--etas", default=""); ap.add_argument("--base", type=int, default=20261014)
    ap.add_argument("--log", default="none", choices=("none", "s0j0")); ap.add_argument("--steps", default="")
    a = ap.parse_args()
    sizes = {int(x) for x in a.sizes.split(",")}
    etas = {x for x in a.etas.split(",") if x}
    used = {int(l.split()[4]) for l in open(S1)}
    short = [int(x) for x in a.steps.split(",")] if a.steps else None
    n = coll = 0
    for stream, (cell, N, elab, L0s, Hs, NS) in enumerate(cells()):
        if N not in sizes or (etas and not any(abs(float(elab) - float(e)) < 5e-7 for e in etas)):
            continue
        L0, H = float(L0s), float(Hs)
        eta = P.eta_of(NS, H, L0)
        dL = 0.5 * P.leff(L0) / math.sqrt(2 * NS * P.cs2(eta))
        teq = (1000 if eta <= 0.60 + 1e-9 else 5000) * STEPS_PER_SIGMA
        win, tail, every = 5000 * STEPS_PER_SIGMA, 20 * STEPS_PER_SIGMA, STEPS_PER_SIGMA
        if short:
            teq, win, tail, every = short
        for pi, (lab, j) in enumerate(POS):
            xs = L0s if j == 0 else f"{L0 + j * dL:.17g}"
            for s in range(a.seeds):
                seed = T.run_seed(a.base, stream, pi, s)
                coll += seed in used
                log = int(a.log == "s0j0" and s == 0 and j == 0)
                print(f"A partA/gen3/{cell}/x_{lab} {j} {s} {seed} {L0s} {Hs} {NS} {xs} {teq} {win} {tail} {every} {log}")
                n += 1
    print(f"# {n} tasks; seeds that stage 1 uses: {coll}", file=sys.stderr)


if __name__ == "__main__":
    main()
