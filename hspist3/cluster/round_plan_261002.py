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
usage: python3 hspist3/cluster/round_plan_261002.py [--slow F] [--eps0-ratio R]   (R = gate-4 eps0_pilot / eps0_plan)
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

def main():
    slow, ratio = arg("--slow", 1.0), arg("--eps0-ratio", None)
    cells, _ = G.planned()
    by = {}
    for c in cells:
        for meth in ("A", "B"):
            by.setdefault(f"{meth}_{c['eta_lab']}", {})[G.cid(c)] = c
    print(f"cost model x {slow:g} (--slow); gate-4 eps0 ratio: {ratio if ratio else 'not applied (planning seeds)'}\n")
    print("| round | array | partition | tasks | cores/task | traj. | core-h | longest cell (h) | --time (h) | throttle | cores at once | wall (h) | scratch GiB |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    lines, checks, tot = [], {}, {}
    for rnd, arrays in ROUNDS:
        for group, N in arrays:
            s, part, P, tlim = sbatch_fields(group)
            order = open(os.path.join(CONF, f"cells_{group}.tsv")).read().split()
            walls, ch, ntr, gib = [], 0.0, 0, 0.0
            for cell in order:
                c = by[group][cell]; n = sum(1 for _ in open(os.path.join(CONF, f"tasks_{group[0]}_{cell}.txt")))
                cpu = (c["cpuA"] if group[0] == "A" else c["cpuB"]) * slow
                if group == "A_0.39" and ratio:                       # gate 4: Tpos scales with eps0^2
                    ns = math.ceil(c["Tpos"] * ratio ** 2 / PR.T_SEED_A); cpu *= 5 * ns / n; n = 5 * ns
                walls.append(math.ceil(n / P) * cpu / n); ch += cpu; ntr += n
                gib += n * (a_seed_mb(c) if group[0] == "A" else B_TRACE_MB) / 1024
            w = schedule(walls, N)
            tot[rnd] = tot.get(rnd, 0) + N * P
            print(f"| {rnd} | conf_{group} | {part} | {len(order)} | {P} | {ntr} | {ch:.1f} | {max(walls):.2f} | {tlim:g} | "
                  f"%{N} | {N * P} | {w:.2f} | {gib:.1f} |")
            lines.append((rnd, f"sbatch --array=1-{len(order)}%{N} cluster/confinement_20261013/conf_{group}.sbatch"))
            checks[group] = check_paths(group, s) + [("longest cell within --time", f"{'yes' if max(walls) <= tlim else '**NO**'} "
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
    return 0

if __name__ == "__main__":
    sys.exit(main())
