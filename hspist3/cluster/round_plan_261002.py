#!/usr/bin/env python3
"""##CHRIS 2026-10-02 (Task R2): round plan for the confinement arrays on KOA -- printed, nothing submitted.

Round 1 = conf_B_0.10, conf_B_0.39, conf_A_0.10 (shared, at the same time); Round 2 = conf_A_0.39 after gate 4.
Per array: tasks (= cells), cores per task (#SBATCH --cpus-per-task), core-hours (pre-registered cost model, the
cells' cpuA / cpuB from paper1_confinement_prereg_20261012.py), wall time at throttle %N, scratch GiB, the sbatch line.

Wall time [INFERENCE]: one cell runs its n trajectories P at a time (xargs -P), so its wall is ceil(n/P) x t_traj with
t_traj = core-h / n; the array's cells start in index order as the %N slots free (list scheduling). The cost model is
the Mac-measured rate; --slow F multiplies every time (KOA Ivy Bridge speed is read off the pilot log, OPEN until then).
Scratch [DERIVATION + DATA]: method B, one 200-period trace per trajectory; the stride rule (tests_20260913.d_stride)
fixes the samples per period, so every trace is the size of the Mac pi/8 pilot traces (DATA: 777157..799567 bytes,
9 files) -> 0.80 MB. Method A, per seed: the event log ev_<seed>.csv logs every wall and divider impact; by the contact
theorem the impact rate per unit wall length is n Z / sqrt(2 pi) (kT = m = 1, n = 4 eta / pi), checked on the E5 log
(pi/8, H 10, L0 10, 400 sigma: 20424 rows over wall length 80 -> 0.638 per sigma per length; theory 0.66) at 69.0
bytes per row; wall length = 2 (2 L0) top/bottom + 2 H sides + 2 H divider faces; record = 5000 + 200 (hold)
sigma-time. Trace tr_<seed>.csv: (300000 + 12000)/600 = 520 rows x 377 bytes (E5 trace: 75449 bytes / 200 rows).
usage: python3 hspist3/cluster/round_plan_261002.py [--slow F]

##CHRIS 2026-10-02 (Task U3): (1) the KOA speed is MEASURED: the method-A pilot (job 14966594, sacct TotalCPU 06:11.764 =
371.764 CPU-s for 20 trajectories of 5000 sigma-time at pi/8, N_s = 50 per side; DATA from Chris's terminal) against the
Mac cost model for the same trajectory (5000 x rate_at(pi/8) ms) -> koa_speed(); it is the default of --slow. It is
measured on method A and applied to method B too [INFERENCE]. (2) Trajectories per cell are read from the task files
(after gate 4 the conf_A_0.39 files carry the gate-4 seeds), and a cell's core-h is the pre-registered cost per
trajectory times that count; --eps0-ratio is gone. (3) The --time rule: --time >= 2 x the longest cell at KOA speed,
rounded up to 15 min, at least 30 min (time_limit_h(); gen_confinement_sbatch.py writes it into the sbatch files).
"""
import math, os, re, sys
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import gate4_pilot_261002 as G
PR = G.PR
CONF = os.path.join(HERE, "confinement_20261013")
SCRATCH = "/mnt/lustre/koa/scratch/charing"
MAX_CORES = 64
B_TRACE_MB = 0.80
A_BYTES_ROW, A_TRACE_MB, A_T = 1409182 / 20424, 520 * 75449 / 200 / 1e6, 5200.0
ROUNDS = [("Round 1", [("B_0.10", 2), ("B_0.39", 2), ("A_0.10", 2)]), ("Round 2", [("A_0.39", 4)])]
KOA_PILOT_CPU_S, KOA_PILOT_TRAJ, KOA_PILOT_NS = 371.764, 20, 50      # ##CHRIS 2026-10-02 (Task U3): sacct, job 14966594

def koa_speed():
    """KOA CPU-s per pilot trajectory / Mac cost model of the same trajectory (5000 sigma-time, pi/8, 2 N_s particles)."""
    import contextlib, io
    with contextlib.redirect_stdout(io.StringIO()): rates = PR.cost_rate()
    mac = PR.T_SEED_A * PR.rate_at(math.pi / 8, rates) * (2 * KOA_PILOT_NS / 100) / 1000
    return KOA_PILOT_CPU_S / KOA_PILOT_TRAJ / mac

def cell_wall_h(n, P, core_h):
    """Wall hours of one cell: n trajectories, P at a time (xargs -P)."""
    return math.ceil(n / P) * core_h / n

def time_limit_h(longest_h):
    return max(0.5, math.ceil(2 * longest_h * 4) / 4)

def arg(name, default):
    return float(sys.argv[sys.argv.index(name) + 1]) if name in sys.argv else default

def sbatch_fields(group):
    s = open(os.path.join(CONF, f"conf_{group}.sbatch")).read()
    get = lambda k: re.search(rf"^#SBATCH --{k}=(\S+)", s, re.M).group(1)
    h, m, *_ = (get("time").split(":") + ["0"])
    return s, get("partition"), int(get("cpus-per-task")), int(h) + int(m) / 60

def a_seed_mb(c):
    n = 4 * c["eta"] / math.pi
    rate = n * PR.Z(c["eta"]) / math.sqrt(2 * math.pi)
    return rate * (4 * c["L0"] + 4 * c["H"]) * A_T * A_BYTES_ROW / 1e6 + A_TRACE_MB

def schedule(walls, N):
    slots = [0.0] * N
    for w in walls:
        i = min(range(N), key=lambda k: slots[k]); slots[i] += w
    return max(slots)

def check_paths(group, s):
    """Static check of where the array writes (sbatch + worker + reducers)."""
    out = []
    out.append(("Slurm stdout/stderr", re.search(r"^#SBATCH --output=(\S+)", s, re.M).group(1) + "  (HOME: ~/harddisks/hspist3/logs)"))
    out.append(("HD_DATA", "$SCRATCH/harddisks/hspist3" if 'HD_DATA="$SCRATCH/harddisks/hspist3"' in s else "**NOT SCRATCH**"))
    out.append(("SCRATCH", SCRATCH if f'SCRATCH="{SCRATCH}"' in s else "**OTHER**"))
    bad = 0
    for cell in open(os.path.join(CONF, f"cells_{group}.tsv")).read().split():
        for l in open(os.path.join(CONF, f"tasks_{group[0]}_{cell}.txt")):
            rel = l.split()[1]; bad += rel.startswith("/") or ".." in rel.split("/")
    out.append(("task output paths (relative to HD_DATA, no '/' or '..')", "all OK" if bad == 0 else f"**{bad} BAD**"))
    return out

def seed_diff(refs=(("plan", "70b2069"), ("gate 4", "303280d"))):
    """##CHRIS 2026-10-02 (Tasks U3, V2): conf_A_0.39 seeds per position in the task files at earlier commits
    (planning seeds 70b2069; gate-4 seeds 303280d) vs now (amendment C3). 'nested' = per position, the shorter seed list
    is the first lines of the longer one, i.e. the files differ only by seeds added or removed at the end."""
    import subprocess
    print("### conf_A_0.39 task files: seeds per position, " + ", ".join(f"{a} (git {r})" for a, r in refs) + " -> now (working tree)\n")
    print("| cell | " + " | ".join(f"{a} ({r})" for a, r in refs) + " | now | lines | nested in each earlier file | seeds now |")
    print("|---|" + "---|" * len(refs) + "---|---|---|---|")
    for cell in open(os.path.join(CONF, "cells_A_0.39.tsv")).read().split():
        f = f"cluster/confinement_20261013/tasks_A_{cell}.txt"
        new = open(os.path.join(HS, f)).read().splitlines(); nn = len(new) // 5; cols, nest = [], True
        for _, r in refs:
            old = subprocess.run(["git", "show", f"{r}:hspist3/{f}"], capture_output=True, text=True, cwd=HS).stdout.splitlines()
            no = len(old) // 5; m = min(no, nn); cols.append(str(no))
            nest &= all(new[k * nn:k * nn + m] == old[k * no:k * no + m] for k in range(5))
        seeds = sorted({int(l.split()[3]) for l in new})
        print(f"| {cell} | {' | '.join(cols)} | {nn} | {len(new)} | {'yes' if nest else '**NO**'} | {seeds[0]}..{seeds[-1]} |")
    print()

def main():
    speed = koa_speed(); slow = arg("--slow", speed)
    seed_diff()
    cells, _ = G.planned()
    by = {}
    for c in cells:
        for meth in ("A", "B"):
            by.setdefault(f"{meth}_{c['eta_lab']}", {})[G.cid(c)] = c
    print(f"KOA speed (measured, pilot 14966594): {KOA_PILOT_CPU_S} CPU-s / {KOA_PILOT_TRAJ} = {KOA_PILOT_CPU_S / KOA_PILOT_TRAJ:.2f} "
          f"CPU-s per trajectory; Mac cost model {KOA_PILOT_CPU_S / KOA_PILOT_TRAJ / speed:.2f} -> factor {speed:.3f}")
    print(f"times below: cost model x {slow:.3f} (--slow); trajectories per cell from the task files\n")
    print("| round | array | partition | tasks | cores/task | traj. | core-h (KOA) | longest cell (h) | --time (h) | rule 2 x longest (h) | throttle | cores at once | wall (h) | scratch GiB |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    lines, checks, tot = [], {}, {}
    for rnd, arrays in ROUNDS:
        for group, N in arrays:
            s, part, P, tlim = sbatch_fields(group)
            order = open(os.path.join(CONF, f"cells_{group}.tsv")).read().split()
            walls, ch, ntr, gib = [], 0.0, 0, 0.0
            for cell in order:
                c = by[group][cell]; n = sum(1 for _ in open(os.path.join(CONF, f"tasks_{group[0]}_{cell}.txt")))
                n_plan = 5 * c["nseed"] if group[0] == "A" else 225                 # trajectories behind cpuA / cpuB
                cpu = (c["cpuA"] if group[0] == "A" else c["cpuB"]) * n / n_plan * slow
                walls.append(cell_wall_h(n, P, cpu)); ch += cpu; ntr += n
                gib += n * (a_seed_mb(c) if group[0] == "A" else B_TRACE_MB) / 1024
            w = schedule(walls, N)
            tot[rnd] = tot.get(rnd, 0) + N * P
            print(f"| {rnd} | conf_{group} | {part} | {len(order)} | {P} | {ntr} | {ch:.1f} | {max(walls):.2f} | {tlim:g} | "
                  f"{time_limit_h(max(walls)):g} | %{N} | {N * P} | {w:.2f} | {gib:.1f} |")
            lines.append((rnd, f"sbatch --array=1-{len(order)}%{N} cluster/confinement_20261013/conf_{group}.sbatch"))
            checks[group] = check_paths(group, s) + [("--time >= 2 x longest cell at KOA speed", f"{'yes' if tlim >= 2 * max(walls) else '**NO**'} "
                                                       f"(margin x{tlim / max(walls):.1f})")]
    print()
    for rnd, n in tot.items():
        print(f"{rnd}: cores at once = {n} (limit {MAX_CORES}) -> {'OK' if n <= MAX_CORES else '**OVER**'}")
    print("\nsbatch lines (from ~/harddisks/hspist3, after `mkdir -p logs`):\n")
    for rnd, l in lines: print(f"    {rnd}:  {l}")
    print("\nWhere each array writes (static check of the sbatch, conf_worker.sh, reduce_A.py, reduce_B.py):\n")
    print("| array | item | result |\n|---|---|---|")
    for g, cc in checks.items():
        for k, v in cc: print(f"| conf_{g} | {k} | {v} |")
    w = open(os.path.join(CONF, "conf_worker.sh")).read()
    resume = ['ls "$cell"/wall_x_positions_L0_*_wallmassfactor_${M}_run${r}.csv >/dev/null 2>&1 && exit 0' in w,
              '[ -s "$d/red_${seed}.csv" ] && exit 0' in w]
    print(f"\nexisting output directory: the arrays do NOT refuse it; conf_worker.sh skips a trajectory whose output exists "
          f"(B: trace run<r>.csv present -> {resume[0]}; A: red_<seed>.csv non-empty -> {resume[1]}), i.e. a resubmitted "
          f"array resumes. mkdir -p creates the directory otherwise.")
    # ##CHRIS 2026-10-02 (Task U4): the resume is now guarded by the build that wrote the directory
    gw = 'guard "$cell"' in w and 'guard "$d"' in w; hs = all("sha256sum --status -c logs/BUILD_KOA_LAST.hash" in sbatch_fields(g)[0]
                                                      for _, arr in ROUNDS for g, _ in arr)
    print(f"build guard in conf_worker.sh (both modes): {gw}; every array checks the binary against logs/BUILD_KOA_LAST.hash: {hs}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
