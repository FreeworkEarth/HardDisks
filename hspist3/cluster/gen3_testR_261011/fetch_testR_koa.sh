#!/bin/bash
# ##CHRIS 2026-10-11 (plan-author decision 28, task R0; 261012 sec. 4.7.57): THE FETCH OF TEST R ON KOA. Chris runs it on the Mac after
# the merge job has ended; ONE rsync over ssh through the DTN (one password and Duo push). Read-only on KOA: nothing there is changed or
# deleted. Here: a target that exists already is not overwritten (rsync --ignore-existing).
# It copies the root /mnt/lustre/koa/scratch/charing/harddisks_gen3/testR_koa WITHOUT the raw data:
#   merged/                the merge record merged/.merge_complete, and every cell folder merged/testR/<engine>/<cell>/m_<M>/ with its
#                          per-trajectory tables (red_nu.csv, acf_runs.npz, nud_runs.csv, runs_health.csv, failed_runs.csv), run.log,
#                          .build_git, .sha256_uncompressed and any .failed_run<r>_* folder's small files;
#   job_logs/              the chunk logs and the merge job's own log (testR_merge.sbatch copies them there);
#   merge_plan_<time>.tsv  the merge plan(s);
# and NOT: the traces (wall_x_positions_*) and psi6(t) files, compressed or not, the run folders (.done_runs/), the chunk folders
# (chunks/, empty after the merge but for anything the merge left in place).
#   bash hspist3/cluster/gen3_testR_261011/fetch_testR_koa.sh
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; HS="$(dirname "$(dirname "$HERE")")"
DST="$HS/experiments_gen3_testR_koa_261011"; mkdir -p "$DST"
SRC="charing@koa-dtn.its.hawaii.edu:/mnt/lustre/koa/scratch/charing/harddisks_gen3/testR_koa"
rsync -av --ignore-existing --exclude='/chunks/' --exclude='wall_x_positions_*' --exclude='psi6_t_*' --exclude='.done_runs/' "$SRC/" "$DST/"
echo "fetched into $DST: $(du -sh "$DST" | cut -f1); merge record: $(head -1 "$DST/merged/.merge_complete" 2>/dev/null || echo MISSING)"
