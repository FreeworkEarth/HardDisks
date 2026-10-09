#!/bin/bash
# ##CHRIS 2026-10-09 (261012 sec. 4.7.20; stage G): PREPARED, NOT RUN. The merge job of TEST G ON KOA: after every chunk of the array
# has ended (submit with --dependency=afterany:<array job id>, then check the chunk logs for FAILED lines), one process merges the
# chunk folders into merged/, reduces the cells (reduce_B.py) and compresses the traces (merge_testG_koa.py).
#   sbatch --dependency=afterany:<ARRAY_JOB_ID> cluster/gen3_koa_261009/merge_testG_koa.sbatch
#SBATCH --job-name=gen3-testG-merge
#SBATCH --partition=shared
#SBATCH --account=uh
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.out
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:?submit with sbatch from ~/harddisks_gen3/hspist3}"
source cluster/koa_env.sh || exit 1
ROOT_NAME=$(basename "$(dirname "$SLURM_SUBMIT_DIR")")
[ "$ROOT_NAME" = harddisks_gen3 ] || { echo "STOP: run from ~/harddisks_gen3, not $ROOT_NAME"; exit 2; }
ROOT=/mnt/lustre/koa/scratch/charing/harddisks_gen3/testG_koa
[ -d "$ROOT/chunks" ] || { echo "STOP: no chunks under $ROOT"; exit 2; }
[ ! -e "$ROOT/merged" ] || { echo "STOP: $ROOT/merged exists -- a merge ran before; not merging twice"; exit 2; }
echo "chunks: $(ls -d "$ROOT"/chunks/chunk_* | wc -l); FAILED lines in the chunk logs: $(grep -l FAILED logs/gen3-testG_*.out 2>/dev/null | wc -l) logs"
python3 cluster/gen3_koa_261009/merge_testG_koa.py "$ROOT" || { echo "STOP: the merge failed"; exit 1; }
echo "merge done; next: the fetch (runsheet step 6) and the KOA verdict on the Mac"
