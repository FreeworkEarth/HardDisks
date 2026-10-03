#!/usr/bin/env python3
"""##CHRIS 2026-10-13: generate the KOA sbatch files for the Paper 1 confinement campaign from the PRE-REGISTERED
cell list (261012 sec. 1, tables G, A, B printed by validation/paper1_confinement_prereg_20261012.py; amendments
sec. 1.9). One SLURM array task = one cell; inside it the cell's trajectories run in parallel through xargs.

Writes into hspist3/cluster/confinement_20261013/: tasks_<cell>.txt (one worker call per line), cells_<file>.tsv
(array index -> cell), the sbatch files, cells_summary.txt (the table printed here) and fetch_confinement.sh.
Placeholders, to be filled from the KOA runbook by Chris: __PARTITION__, __ACCOUNT__, __SCRATCH__.
Data layout under $HD_DATA (= __SCRATCH__/harddisks/hspist3) is the Mac layout, relative to hspist3/.

##CHRIS 2026-10-02 (Task K1): placeholders FILLED from the KOA facts (Chris's terminal, 2026-10-03 UTC): partition sandbox
for the method-A pilot and shared for the arrays, account uh, scratch /mnt/lustre/koa/scratch/charing, user charing. Every
job sources cluster/koa_env.sh (compiler/GCC/14.3.0 + ~/envs/hd), prints `gcc --version` and the binary's --version, and
stops unless the binary is the clean koa build of the checkout's HEAD. The cell and task lists are unchanged.

##CHRIS 2026-10-02 (Tasks U3, U4):
  - Gate 4 (261012 sec. 1.7 item 4, sec. 1.8): if cluster/confinement_20261013/gate4_pi8_result.txt (written on PASS by
    cluster/gate4_pilot_261002.py) exists, the pi/8 cells (eta_lab 0.39) use the pilot's eps0. T per position is
    proportional to eps0^2 in the pre-registered rule, so Tpos -> Tpos (eps0_pilot/eps0_plan)^2, seeds per position =
    ceil(Tpos/5000), and the cell's core-h scales with the seed count. eta = 0.10 cells and the pilot are unchanged.
  - --time = round_plan_261002.time_limit_h(longest cell at the measured KOA speed): >= 2 x the longest cell, rounded up
    to 15 min, at least 30 min. The pilot keeps its 1 hour (it ran in 41 s).
  - The arrays check the binary against logs/BUILD_KOA_LAST.hash (written by cluster/build_koa.sh), not against
    `git rev-parse HEAD`, so a later `git pull` cannot stop a running campaign; they export HD_BUILD for the build
    guard in conf_worker.sh. The sbatch comment carries the throttle of the round plan.
    (2026-10-03, Task W1: the "stop if flock is missing" check is removed; the guard uses a mkdir lock.)
"""
import contextlib, io, math, os, sys
from scipy.optimize import brentq
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(HS, "validation")); sys.path.insert(0, HS)
import tests_20260913 as T
import paper1_confinement_prereg_20261012 as PR
sys.path.insert(0, HERE)
import round_plan_261002 as RP                       # ##CHRIS 2026-10-02 (Task U3): KOA speed, --time rule, throttles
OUT = os.path.join(HERE, "confinement_20261013")
REL_B = "experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/confinement_B_20261013"
REL_A = "experiments_energy_transfer/paper1_confinement_A_20261013"
SEEDS_B, T_SEED_A, EVERY_A = 25, 5000.0, 600
SCRATCH = "/mnt/lustre/koa/scratch/charing"           # ##CHRIS 2026-10-02 (Task K1): KOA facts
CPUS = {"B": 8, "A": 16}

def kroot(a): return brentq(lambda k: math.cos(k) / math.sin(k) - a * k, 1e-9, math.pi - 1e-9)
def cid(c): return f"e{'0p10' if c['eta_lab'] == '0.10' else 'pi8'}_{c['scan']}_H{c['H']:g}_L{c['L0']:g}"

def main():
    with contextlib.redirect_stdout(io.StringIO()): C, A = PR.main()
    seen, cells = set(), []
    for c in C:                                   # the pi/8 aspect L0/H = 1 cell IS the anchor: not run twice
        key = (c["eta_lab"], round(c["H"], 6), round(c["L0"], 6))
        if key in seen: continue
        seen.add(key); cells.append(c)
    # ##CHRIS 2026-10-02 (Task U3): gate 4 -- the pi/8 cells take the pilot's eps0 (see the docstring)
    g4 = os.path.join(OUT, "gate4_pi8_result.txt")
    if os.path.exists(g4):
        kv = dict(l.split()[:2] for l in open(g4) if l.strip() and not l.startswith("#"))
        assert kv["verdict"] == "PASS", "gate 4 did not pass"
        # ##CHRIS 2026-10-02 (Task V2): amendment C3 (261012 sec. 1.9) -- the upper 1-sigma bound eps0_pilot + err is used
        eps_use = float(kv["eps0_c3_upper"]) if "eps0_c3_upper" in kv else float(kv["eps0_pilot"])
        r2 = (eps_use / float(kv["eps0_plan"])) ** 2
        for c in cells:
            if c["eta_lab"] != "0.39": continue
            n_old = c["nseed"]; c["Tpos"] *= r2; c["nseed"] = math.ceil(c["Tpos"] / T_SEED_A)
            c["cpuA"] *= c["nseed"] / n_old
        print(f"gate 4 applied to the pi/8 cells: eps0 {float(kv['eps0_plan']):.4f} -> {eps_use:.4f}"
              f"{' (amendment C3: pilot ' + format(float(kv['eps0_pilot']), '.4f') + ' + 1 sigma)' if 'eps0_c3_upper' in kv else ''}\n")
    speed = RP.koa_speed(); thr = {g: N for _, arr in RP.ROUNDS for g, N in arr}; thr["A_pilot"] = 1
    os.makedirs(OUT, exist_ok=True); rows, groups, ntraj = [], {}, {}
    for i, c in enumerate(cells):
        name = cid(c); cs = math.sqrt(PR.cs2(c["eta"])); Le = c["Le"]; base = 20261100 + i
        # method B: nine masses x 25 seeds, exact seeds from the A1v2 seed hash
        lines = []
        for mi, M in enumerate(c["Ms"]):
            nu = cs * kroot(M / (2.0 * c["Ns"])) / (2 * math.pi * Le); stride = T.d_stride(nu)
            for r in range(SEEDS_B):
                lines.append(f"B {REL_B}/{name}/m_{M} {M} {r} {T.run_seed(base, 0, mi, r)} {c['L0']:.6f} {c['H']:.6f} {c['Ns']} {stride} {base}")
        open(os.path.join(OUT, f"tasks_B_{name}.txt"), "w").write("\n".join(lines) + "\n")
        rows.append(("B", name, c, len(lines), c["cpuB"], f"{REL_B}/{name}")); ntraj[("B", name)] = len(lines)
        groups.setdefault(f"B_{c['eta_lab']}", []).append((name, c["cpuB"]))
        # method A: five held positions x nseed seeds of <= 5000 sigma-time
        lines = []
        for j, lab in zip((-2, -1, 0, 1, 2), ("m2", "m1", "0", "p1", "p2")):
            xw = c["L0"] + j * c["dL"]
            for s in range(c["nseed"]):
                lines.append(f"A {REL_A}/{name}/x_{lab} {xw:.6f} {9700 + s} {c['L0']:.6f} {c['H']:.6f} {c['Ns']} "
                             f"{int(60 * T_SEED_A)} {EVERY_A}")
        open(os.path.join(OUT, f"tasks_A_{name}.txt"), "w").write("\n".join(lines) + "\n")
        rows.append(("A", name, c, len(lines), c["cpuA"], f"{REL_A}/{name}")); ntraj[("A", name)] = len(lines)
        groups.setdefault(f"A_{c['eta_lab']}", []).append((name, c["cpuA"]))
    # method A pilot at pi/8 (gate 4): 4 seeds per position at the anchor
    anc = [c for c in cells if c["eta_lab"] == "0.39" and c["scan"] == "H" and abs(c["H"] - 10) < 1e-9][0]
    lines = [f"A {REL_A}/pilot_{cid(anc)}/x_{lab} {anc['L0'] + j * anc['dL']:.6f} {9700 + s} {anc['L0']:.6f} {anc['H']:.6f} "
             f"{anc['Ns']} {int(60 * T_SEED_A)} {EVERY_A}" for j, lab in zip((-2, -1, 0, 1, 2), ("m2", "m1", "0", "p1", "p2")) for s in range(4)]
    open(os.path.join(OUT, f"tasks_A_pilot_{cid(anc)}.txt"), "w").write("\n".join(lines) + "\n")
    groups["A_pilot"] = [(f"pilot_{cid(anc)}", anc["cpuA"] * 20 / (5 * anc["nseed"]))]; ntraj[("A", f"pilot_{cid(anc)}")] = 20
    for g, items in groups.items():
        meth = g[0]; cpus = CPUS[meth]; part = "sandbox" if g == "A_pilot" else "shared"
        # ##CHRIS 2026-10-02 (Task U3): was hrs = max(1, math.ceil(1.5 * max(h for _, h in items) / cpus)) (Mac speed, 1.5 x)
        hrs = 1.0 if g == "A_pilot" else RP.time_limit_h(max(RP.cell_wall_h(ntraj[(meth, n)], cpus, h * speed) for n, h in items))
        tstr = f"{int(hrs)}:{int(round((hrs % 1) * 60)):02d}:00"
        open(os.path.join(OUT, f"cells_{g}.tsv"), "w").write("\n".join(n for n, _ in items) + "\n")
        warn = ("# !! SUBMIT ONLY AFTER the pi/8 method-A pilot (conf_A_pilot.sbatch) has measured the noise coefficient and\n"
                "# !! the seeds per position have been recomputed by the 261012 sec. 1.4 rule (gate 4). Seeds here are the\n"
                "# !! planning values (Table A, eps_0 measured at eta = 0.10 only).\n") if g == "A_0.39" and not os.path.exists(g4) else (
               "# Gate 4 PASSED (2026-10-02): seeds per position recomputed by the 261012 sec. 1.4 rule with the pi/8 pilot's eps0,\n"
               "# taken at its upper 1-sigma bound by amendment C3 (261012 sec. 1.9; cluster/confinement_20261013/gate4_pi8_result.txt).\n"
               "# Submit after the go for Round 2, from the same build as Round 1 (runsheet step 8e).\n") if g == "A_0.39" else ""
        open(os.path.join(OUT, f"conf_{g}.sbatch"), "w").write(f"""#!/bin/bash
# ##CHRIS 2026-10-13: Paper 1 confinement campaign, {'method B (free divider, A1v2 protocol)' if meth == 'B' else 'method A (held divider)'}, group {g}.
# Generated by cluster/gen_confinement_sbatch.py from the pre-registered cell list (261012 sec. 1). NOT LAUNCHED.
# One array task = one cell (cells_{g}.tsv line N); the cell's trajectories run {cpus} at a time.
{warn}# FILLED 2026-10-02 (Task K1): partition {part}, account uh, scratch {SCRATCH} (koa_scratch:
# files are DELETED 90 days after last write -- copy results back before then). Submit from ~/harddisks/hspist3 after
# `mkdir -p logs`; the binary must already be built by cluster/build_koa.sh (this job does not build), which records
# logs/BUILD_KOA_LAST.hash. After any `git pull`, rebuild before a NEW submission (runsheet step 8).
# --time: >= 2 x the longest cell at the measured KOA speed (x{speed:.3f}, pilot 14966594).
#   sbatch --array=1-{len(items)}%{thr[g]} cluster/confinement_20261013/conf_{g}.sbatch
#SBATCH --job-name=conf-{g}
#SBATCH --partition={part}
#SBATCH --account=uh
#SBATCH --time={tstr}
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.out
set -uo pipefail
SCRATCH="{SCRATCH}"
cd "$SLURM_SUBMIT_DIR"                       # hspist3/ (source, scripts, the koa-built binary)
source cluster/koa_env.sh || exit 1
export HD_BIN="$PWD/00ALLINONE" HD_DATA="$SCRATCH/harddisks/hspist3"
gcc --version | head -1; "$HD_BIN" --version | head -2
# ##CHRIS 2026-10-02 (Task U4): checked against the hash recorded at build time, not `git rev-parse HEAD`
[ -s logs/BUILD_KOA_LAST.hash ] || {{ echo "STOP: no logs/BUILD_KOA_LAST.hash -- build with cluster/build_koa.sh first"; exit 1; }}
sha256sum --status -c logs/BUILD_KOA_LAST.hash || {{ echo "STOP: ./00ALLINONE is not the build recorded in logs/BUILD_KOA_LAST.hash"; exit 1; }}
export HD_BUILD="$("$HD_BIN" --version | head -1)"
echo "$HD_BUILD" | grep -Eq -- "git [0-9a-f]+  target koa" || {{ echo "STOP: not a clean koa build: $HD_BUILD"; exit 1; }}
# ##CHRIS 2026-10-03 (Task W1): the `command -v flock` check is gone -- the build guard locks with mkdir now
CELL=$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" cluster/confinement_20261013/cells_{g}.tsv)
TASKS=cluster/confinement_20261013/tasks_{'A_' if meth == 'A' else 'B_'}$CELL.txt
echo "cell $CELL: $(wc -l < "$TASKS") trajectories, {cpus} in parallel"
xargs -P {cpus} -L 1 bash cluster/confinement_20261013/conf_worker.sh < "$TASKS"
{'python3 cluster/confinement_20261013/reduce_B.py "$HD_DATA/' + REL_B + '/$CELL"' if meth == 'B' else 'echo "method A: per-seed red_<seed>.csv written by the worker"'}
echo "cell $CELL done; failures: $(grep -c FAILED logs/conf-{g}_${{SLURM_ARRAY_JOB_ID}}_${{SLURM_ARRAY_TASK_ID}}.out || true)"
""")
    L = ["| method | cell id | N_s | H | L_0 | trajectories (seeds) | est. core-h | output dir (relative to hspist3/) |",
         "|---|---|---|---|---|---|---|---|"]
    for meth, name, c, n, ch, od in rows:
        L.append(f"| {meth} | {name} | {c['Ns']} | {c['H']:g} | {c['L0']:g} | {n} ({'9 masses x 25' if meth == 'B' else '5 positions x ' + str(c['nseed'])}) | {ch:.2f} | `{od}` |")
    tot = {m: sum(r[4] for r in rows if r[0] == m) for m in "AB"}
    L += ["", f"Totals: method A {tot['A']:.1f} core-h, method B {tot['B']:.1f} core-h, all {tot['A'] + tot['B']:.1f} core-h "
          f"(the 261012 cost table counts the pi/8 anchor once; so does this list). Plus the pi/8 method-A pilot "
          f"({len(lines)} trajectories, {groups['A_pilot'][0][1]:.2f} core-h)."]
    fetch = f"""#!/usr/bin/env bash
# ##CHRIS 2026-10-13: copy the confinement campaign back FROM KOA (run on the Mac, from the repo root).
# Summaries only, plus the full pilot cells; full trajectories stay on KOA scratch (deleted after 90 days).
# KOA_USER and SCRATCH filled 2026-10-02 (Task K1). Nothing on either side is deleted.
KOA_USER=charing; SCRATCH="{SCRATCH}"; DTN=$KOA_USER@koa-dtn.its.hawaii.edu; R=$SCRATCH/harddisks/hspist3
SUM=(--prune-empty-dirs --include='*/' --include='red_*.csv' --include='red_nu.csv' --include='acf_runs.npz'
     --include='run.log' --include='run_*.log' --include='summary_*.csv' --include='command*.txt' --exclude='*')
rsync -av "${{SUM[@]}}" "$DTN:$R/{REL_B}/" "hspist3/{REL_B}/"
rsync -av "${{SUM[@]}}" "$DTN:$R/{REL_A}/" "hspist3/{REL_A}/"
# full pilot cells (every file):
rsync -av "$DTN:$R/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/confinement_pilot_20261013/koa_pi8_H10_L10/" \\
          "hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/confinement_pilot_20261013/koa_pi8_H10_L10/"
rsync -av "$DTN:$R/{REL_A}/pilot_{cid(anc)}/" "hspist3/{REL_A}/pilot_{cid(anc)}/"
"""
    open(os.path.join(OUT, "fetch_confinement.sh"), "w").write(fetch)
    txt = "\n".join(L) + "\n"; open(os.path.join(OUT, "cells_summary.txt"), "w").write(txt); print(txt)
    print("rsync back (written to cluster/confinement_20261013/fetch_confinement.sh):\n"); print(fetch)
    print("files:", ", ".join(sorted(f for f in os.listdir(OUT) if f.endswith(".sbatch"))))

if __name__ == "__main__":
    main()
