#!/bin/bash
# ##CHRIS 2026-10-12 (plan-author decision 39 item 5 f; registration 261012 sec. 4.8.8, "one merge per block when its array is complete";
# gate G1): THE MERGE OF ONE STAGE 1 GROUP ON KOA, one process, from the scripts clone ~/harddisks_gen3_k2 (updated in place to the stage 1
# package commit), after every chunk of the group has ended. It is cluster/gen3_testR_261011/testR_merge.sbatch (c6285af7) with THESE
# CHANGES, all of them: the group is the first argument (sbatch ... s1_merge.sbatch <group>), so the root is .../stage1/<group>, the job
# name of its chunks gen3-s1-<group>, its chunks the group's rows of the chunk table; the chunk check s1_chunk_logs_koa.sh (per group);
# the merge script merge_s1_koa.py with --chunks and --group; the input manifests s1_inputs.sha256 and s1_merge_inputs.sha256; the scripts
# whose SHA-256 it prints. Every other rule is Test R's:
#   - it does not depend on the arrays' job ids: Chris submits it when squeue lists no gen3-s1-<group> job; it checks that itself;
#   - every check comes before any move, and a failed check stops it with exit 2 and NOTHING moved: the inputs as delivered, the task
#     list as registered, no gen3-s1-<group> job pending or running, the chunk check (every chunk of the group has its folder and EXACTLY
#     ONE log logs/gen3-s1-<group>_*_<j>.out with its line "chunk <j> done");
#   - merged/ must not exist (exit 3): a merge is never repeated;
#   - then merge_s1_koa.py: the plan, and only a clean plan is executed (merged/.merge_complete), then the reduction (reduce_B.py and
#     reduce_s1_koa.py per cell folder) and the compression of the traces; the job logs copied to <root>/job_logs/ (no overwrite).
# exit: 0 merged; 2 a check or the merge plan stopped it, NOTHING moved (it may be submitted again once the cause is removed);
#       3 refused (merged/ exists); 1 the merge failed during its execution: merged/ has no .merge_complete; inspect by hand.
#   cd ~/harddisks_gen3_k2/hspist3 && sbatch --job-name=gen3-s1-merge-<group> cluster/gen3_stage1_261012/s1_merge.sbatch <group>
#SBATCH --job-name=gen3-s1-merge
#SBATCH --partition=shared
#SBATCH --account=uh
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --exclude=cn-05-08-07
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.out
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:?submit with sbatch from ~/harddisks_gen3_k2/hspist3}"
P=cluster/gen3_stage1_261012; GRP=${1:?the group (B100, B400, D, B900, S, B1600) is the first argument}
case "$GRP" in B100|B400|D|B900|S|B1600) ;; *) echo "STOP: '$GRP' is not a stage 1 group; NOTHING was moved"; exit 2 ;; esac
ROOT_NAME=$(basename "$(dirname "$SLURM_SUBMIT_DIR")")
[ "$ROOT_NAME" = harddisks_gen3_k2 ] || { echo "STOP: the stage 1 merge runs from ~/harddisks_gen3_k2, not $ROOT_NAME"; exit 2; }
sha256sum --status -c "$P/s1_inputs.sha256" && sha256sum --status -c "$P/s1_merge_inputs.sha256" \
  || { echo "STOP: the inputs are not the delivered ones (s1_inputs.sha256, s1_merge_inputs.sha256); NOTHING was moved"; exit 2; }
source cluster/koa_env.sh || { echo "STOP: cluster/koa_env.sh failed; NOTHING was moved"; exit 2; }
ROOT=/mnt/lustre/koa/scratch/charing/harddisks_gen3/stage1/$GRP
TASKS="$P/tasks_stage1_koa_261012.txt"; CT="$P/chunks_stage1_koa_261012.tsv"; PFX="gen3-s1-$GRP"
[ "$(sha256sum "$TASKS" | cut -d' ' -f1)" = "$(cat "$P/tasks_stage1_koa_261012.sha256")" ] || { echo "STOP: the task list is not the registered one"; exit 2; }
NCH=$(awk -F'\t' -v g="$GRP" 'NR > 1 && $5 == g' "$CT" | wc -l | tr -d ' ')
[ -d "$ROOT/chunks" ] || { echo "STOP: no chunks under $ROOT"; exit 2; }
[ ! -e "$ROOT/merged" ] || { echo "STOP: $ROOT/merged exists -- a merge ran before; not merging twice"; exit 3; }
Q=$(squeue -h -u "$USER" -n "$PFX" -o '%i %T' 2>&1); qrc=$?
[ "$qrc" -eq 0 ] || { echo "STOP: squeue failed (exit $qrc): $Q"; exit 2; }
[ -z "$Q" ] || { echo "STOP: $GRP chunks are still pending or running: $(printf '%s\n' "$Q" | head -3 | tr '\n' ' ')"; exit 2; }
echo "merge of $ROOT by the scripts of ~/harddisks_gen3_k2 (no git at run time; both input manifests: OK); the scripts it runs:"
sha256sum "$P/merge_s1_koa.py" "$P/s1_chunk_logs_koa.sh" "$P/reduce_s1_koa.py" cluster/confinement_20261013/reduce_B.py | sed 's/^/  /'
echo "group $GRP: $NCH chunks of the chunk table; no $PFX job pending or running"
mkdir -p "$ROOT/job_logs"; cp -n -p logs/"$PFX"_*.out "$ROOT/job_logs/"
echo "chunk logs in $ROOT/job_logs/ (copied, no overwrite): $(ls "$ROOT/job_logs" | grep -c "^${PFX}_")"
bash "$P/s1_chunk_logs_koa.sh" logs "$PFX" "$ROOT/chunks" "$CT" "$GRP" \
  || { echo "STOP: the chunk set is not complete (listed above); NOTHING was moved"; exit 2; }
NF=$(cat logs/"$PFX"_*.out | grep -c FAILED); echo "FAILED lines in the chunk logs (logs/): $NF"
[ "$NF" -eq 0 ] || grep -H FAILED logs/"$PFX"_*.out | head -100
python3 "$P/merge_s1_koa.py" "$ROOT" --tasks "$TASKS" --chunks "$CT" --group "$GRP"; rc=$?
case $rc in
  0) echo "merge done; next: the fetch (fetch_s1_koa.sh $GRP) and the tables on the Mac" ;;
  2) echo "STOP: the merge plan has problems (listed above); NOTHING was moved; once the cause is removed this job may be submitted again"; exit 2 ;;
  3) echo "STOP: refused (merged/ exists)"; exit 3 ;;
  *) echo "STOP: the merge failed during its execution (exit $rc): merged/ has no .merge_complete; inspect by hand, do not submit again"; exit 1 ;;
esac
echo "this job's log is copied to $ROOT/job_logs/ as the last command"
cp -n -p "logs/${SLURM_JOB_NAME:-gen3-s1-merge}_${SLURM_JOB_ID}.out" "$ROOT/job_logs/"
exit 0
