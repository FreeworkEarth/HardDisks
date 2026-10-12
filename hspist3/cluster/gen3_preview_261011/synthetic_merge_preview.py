#!/usr/bin/env python3
"""##CHRIS 2026-10-11 (plan-author decision 32, task P; 261012 sec. 4.7.63): the merge's expected lines for the runsheet's block H, from a
SYNTHETIC full-size chunk tree. No simulation runs: every file is a few bytes, written as testG_worker_koa.sh leaves a finished run in its
cell folder (the shapes of Test R's unit test test_merge_koa_plan.py, b_run(): the build line, the trace, the run.log section, gen3's
psi6(t) .gz and its .sha256_uncompressed line, the run's own folder .done_runs/run<r>/ with the driver's record folder inside). The
layout and the run indices follow the REGISTERED task list (tasks_preview_koa_261011.txt, chunk size from preview_chunk_size.txt); the
trace names carry the driver's L0 token (int)(L0 * 10). Then Test R's committed merge_koa.py runs as the merge job runs it, with
--no-reduce (the synthetic traces are not real), on a new temporary folder (under TMPDIR, kept). It prints the merge's lines; block H
quotes them, together with the reduction's lines that only the rehearsal and KOA print.
usage (from hspist3/): python3 cluster/gen3_preview_261011/synthetic_merge_preview.py
"""
import os, subprocess, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__)); CL = os.path.dirname(HERE)
MERGE = os.path.join(CL, "gen3_testR_261011", "merge_koa.py")
TASKS = os.path.join(HERE, "tasks_preview_koa_261011.txt"); K = int(open(os.path.join(HERE, "preview_chunk_size.txt")).read())
BUILD = "00ALLINONE  git a643bc7  target koa"


def w(p, text):
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p, "w") as fh: fh.write(text)


def b_run(cell, L0, M, r, seed):
    tok = int(float(L0) * 10)
    w(os.path.join(cell, ".build_git"), BUILD + "\n")
    w(os.path.join(cell, f"wall_x_positions_L0_{tok}_wallmassfactor_{M}_run{r}.csv"), f"t,x\n0,{r}\n")
    with open(os.path.join(cell, "run.log"), "a") as fh: fh.write(f"\n##RUN 2026-10-12 12:00:00 run {r} seed {seed} node cn-x (2 s)\nrun = {r},\n")
    w(os.path.join(cell, f"psi6_t_L0_{tok}_wallmassfactor_{M}_run{r}.csv.gz"), f"gz {r}\n")
    with open(os.path.join(cell, ".sha256_uncompressed"), "a") as fh: fh.write(f"{r:064x}  psi6_t_L0_{tok}_wallmassfactor_{M}_run{r}.csv\n")
    d = os.path.join(cell, ".done_runs", f"run{r}")
    for n in ("00_COMMAND.md", "01_PLOT_COMMANDS.md", "00_COMMAND_AND_PLOT_COMMAND.md", "speed_of_sound_psi6.csv",
              "speed_of_sound_batch_status.json", "stdout.log"):
        w(os.path.join(d, n), f"{n} of run {r}\n")
    w(os.path.join(d, "experiments_speed_of_sound", "EDMD", "mode0_real_units", f"simulation_10_12_26_12_00_{r:02d}", "00_COMMAND.md"), "x\n")


def main():
    L = [l.split() for l in open(TASKS) if l.strip()]
    root = tempfile.mkdtemp(prefix="preview_synthetic_merge_")
    for i, p in enumerate(L):
        j = i // K + 1
        b_run(os.path.join(root, "chunks", f"chunk_{j}", p[1]), p[5], int(p[2]), int(p[3]), int(p[4]))
    print(f"# SYNTHETIC full-size merge of the preview (no simulation): {len(L)} registered lines, chunk size {K}, "
          f"{-(-len(L) // K)} chunks, in {root}\n")
    q = subprocess.run([sys.executable, MERGE, root, "--tasks", TASKS, "--chunk-size", str(K), "--no-reduce"], capture_output=True, text=True)
    print(q.stdout + q.stderr)
    print(f"merge_koa.py exit {q.returncode}")
    sys.exit(q.returncode)


if __name__ == "__main__":
    main()
