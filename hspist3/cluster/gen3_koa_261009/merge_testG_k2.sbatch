#!/bin/bash
# ##CHRIS 2026-10-10 (plan-author decision 20, task 4; 261012 sec. 4.7.47): RUNSHEET ADDENDUM STEP K2-4, the merge of TEST G ON KOA,
# from the SECOND clone ~/harddisks_gen3_k2 with the merge_testG_koa.py of the fix commit. The registered merge (a643bc7) stops with
# "name collision .../x_0/00_COMMAND.md" at the second chunk that holds A-fixed runs (found by the Mac rehearsal); this one keeps each
# chunk's copy of those two per-run files under <name>_chunk<j>. Everything else is merge_testG_koa.sbatch's: one process, after every
# chunk of the array has ended, into merged/ under the same root; it never merges twice. The chunk logs are in
# ~/harddisks_gen3/hspist3/logs/, where the arrays ran.
# ##CHRIS 2026-10-10 (plan-author decisions 21 and 22; 261012 sec. 4.7.49): the merge PLANS BEFORE IT MOVES (merge_testG_koa.py), and
# this job does NOT depend on the array's job id: Chris submits it after the array has ended (no --dependency), and before the merge
# the job checks for itself, stopping with exit 2 and NOTHING moved on any failure:
#   - no job named gen3-testG of this user is pending or running (squeue);
#   - the task list of ~/harddisks_gen3 is the registered one (its SHA-256 file, as every chunk checks it); it gives the chunk count;
#   - every chunk j has its folder and a log ~/harddisks_gen3/hspist3/logs/gen3-testG_*_<j>.out with the line "chunk <j> done";
# ##CHRIS 2026-10-10 (plan-author decision 24, point 2b; sec. 4.7.52): that check is testG_chunk_logs_koa.sh (unit test
#     test_testG_chunk_logs_261010.py): EXACTLY ONE log of a chunk has its done line, and every other log of the chunk is listed by
#     name (an attempt that stopped in the start-up checks, with its STOP line; an attempt that started and did not finish);
#   - merged/ does not exist (exit 3; the merge script refuses it as well: a completed merge is never repeated).
# The FAILED lines of the chunk logs are printed (the verdict's inventory counts failed runs; they do not stop the merge).
# exit: 0 merged; 2 a check or the merge plan stopped it, NOTHING moved: once the cause is removed the job may be submitted again;
#       3 refused (merged/ exists); 1 the merge failed during its execution: merged/ has no .merge_complete; inspect by hand.
#   cd ~/harddisks_gen3_k2/hspist3 && sbatch cluster/gen3_koa_261009/merge_testG_k2.sbatch
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
A="$HOME/harddisks_gen3/hspist3"; G="$A/logs"
TASKS="$A/cluster/gen3_koa_261009/tasks_testG_koa_261009.txt"
EXPECT=$(cat "$A/cluster/gen3_koa_261009/tasks_testG_koa_261009.sha256" 2>/dev/null)
[ -n "$EXPECT" ] && [ "$(sha256sum "$TASKS" | cut -d' ' -f1)" = "$EXPECT" ] || { echo "STOP: the task list of ~/harddisks_gen3 is not the registered one"; exit 2; }
NL=$(wc -l < "$TASKS"); NCH=$(( (NL + 24) / 25 ))
[ -d "$ROOT/chunks" ] || { echo "STOP: no chunks under $ROOT"; exit 2; }
[ ! -e "$ROOT/merged" ] || { echo "STOP: $ROOT/merged exists -- a merge ran before; not merging twice"; exit 3; }
Q=$(squeue -h -u "$USER" -n gen3-testG -o '%i %T' 2>&1); qrc=$?
[ "$qrc" -eq 0 ] || { echo "STOP: squeue failed (exit $qrc): $Q"; exit 2; }
[ -z "$Q" ] || { echo "STOP: Test G chunks are still pending or running: $(printf '%s\n' "$Q" | head -3 | tr '\n' ' ')"; exit 2; }
echo "merge of $ROOT by the scripts of clone HEAD $(git rev-parse HEAD)"
echo "the registered list: $NL task lines, $NCH chunks; no gen3-testG job pending or running"
bash cluster/gen3_koa_261009/testG_chunk_logs_koa.sh "$G" "$ROOT/chunks" "$NCH" \
  || { echo "STOP: the chunk set is not complete (listed above); NOTHING was moved"; exit 2; }
NF=$(cat "$G"/gen3-testG_*.out | grep -c FAILED); echo "FAILED lines in the chunk logs ($G): $NF"
[ "$NF" -eq 0 ] || grep -H FAILED "$G"/gen3-testG_*.out | head -100
python3 cluster/gen3_koa_261009/merge_testG_koa.py "$ROOT"; rc=$?
case $rc in
  0) echo "merge done; next: the fetch (runsheet step 6) and the KOA verdict on the Mac" ;;
  2) echo "STOP: the merge plan has problems (listed above); NOTHING was moved; once the cause is removed this job may be submitted again"; exit 2 ;;
  3) echo "STOP: refused (merged/ exists)"; exit 3 ;;
  *) echo "STOP: the merge failed during its execution (exit $rc): merged/ has no .merge_complete; inspect by hand, do not submit again"; exit 1 ;;
esac
