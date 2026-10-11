#!/bin/bash
# ##CHRIS 2026-10-10 (plan-author decision 20, task 4; 261012 sec. 4.7.47): RUNSHEET ADDENDUM STEP K2-4, the merge of TEST G ON KOA,
# from the SECOND clone ~/harddisks_gen3_k2 with the merge_testG_koa.py of the fix commit. The registered merge (a643bc7) stops with
# "name collision .../x_0/00_COMMAND.md" at the second chunk that holds A-fixed runs (found by the Mac rehearsal); this one keeps each
# chunk's copy of those two per-run files under <name>_chunk<j>. Everything else is merge_testG_koa.sbatch's: one process, after every
# chunk of the array has ended (--dependency=afterany), into merged/ under the same root; it never merges twice. The chunk logs are in
# ~/harddisks_gen3/hspist3/logs/, where the arrays ran.
#   cd ~/harddisks_gen3_k2/hspist3 && sbatch --dependency=afterany:<ARRAY_JOB_ID> cluster/gen3_koa_261009/merge_testG_k2.sbatch
#SBATCH --job-name=gen3-testG-merge-k2
#SBATCH --partition=shared
#SBATCH --account=uh
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.out
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:?submit with sbatch from ~/harddisks_gen3_k2/hspist3}"
source cluster/koa_env.sh || exit 1
ROOT_NAME=$(basename "$(dirname "$SLURM_SUBMIT_DIR")")
[ "$ROOT_NAME" = harddisks_gen3_k2 ] || { echo "STOP: the K2 jobs run from ~/harddisks_gen3_k2, not $ROOT_NAME"; exit 2; }
ROOT=/mnt/lustre/koa/scratch/charing/harddisks_gen3/testG_koa
[ -d "$ROOT/chunks" ] || { echo "STOP: no chunks under $ROOT"; exit 2; }
[ ! -e "$ROOT/merged" ] || { echo "STOP: $ROOT/merged exists -- a merge ran before; not merging twice"; exit 2; }
G="$HOME/harddisks_gen3/hspist3/logs"
echo "merge of $ROOT by the scripts of clone HEAD $(git rev-parse HEAD)"
echo "chunks: $(ls -d "$ROOT"/chunks/chunk_* | wc -l); FAILED lines in the chunk logs ($G): $(grep -l FAILED "$G"/gen3-testG_*.out 2>/dev/null | wc -l) logs"
python3 cluster/gen3_koa_261009/merge_testG_koa.py "$ROOT" || { echo "STOP: the merge failed"; exit 1; }
echo "merge done; next: the fetch (runsheet step 6) and the KOA verdict on the Mac"
