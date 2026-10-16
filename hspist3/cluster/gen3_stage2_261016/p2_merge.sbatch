#!/bin/bash
# ##CHRIS 2026-10-15 (plan-author decision 57; notes 261012 sec. 4.7.135): THE MERGE OF THE STAGE 2 PILOT P2 ON KOA, one process, from
# ~/harddisks_gen3_p2/hspist3, after every chunk has ended. It is cluster/gen3_stage1_261012/s1_merge.sbatch with THESE CHANGES, all of them:
# the one group P2 (no argument), its root .../harddisks_gen3/stage2/P2, its job name gen3-p2, its task list and chunk table, the input
# manifest p2_inputs.sha256, the clone ~/harddisks_gen3_p2, the time limit 24:00:00 (decision 57), and ONE STEP ADDED after the stage 1
# merge (merge_s1_koa.py with its reduction: reduce_B.py and reduce_s1_koa.py per cell folder, unchanged): p2_reduce.py, the P2 tables
# (per run: the own fit, the spectrum share, the drift, psi6 and the orientation per compartment, the strain at the start and per snapshot).
# The chunk check is the stage 1 one (s1_chunk_logs_koa.sh, unchanged). Every rule of s1_merge.sbatch otherwise:
#   - it does not depend on the array's job id: Chris submits it when squeue lists no gen3-p2 job; it checks that itself;
#   - every check comes before any move, and a failed check stops it with exit 2 and NOTHING moved: the inputs as delivered, the task
#     list as delivered, no gen3-p2 job pending or running, the chunk check (every chunk has its folder and EXACTLY ONE log
#     logs/gen3-p2_*_<j>.out with its line "chunk <j> done");
#   - merged/ must not exist (exit 3): a merge is never repeated;
#   - the job logs are copied to <root>/job_logs/ (no overwrite).
# exit: 0 merged and reduced; 2 a check or the merge plan stopped it, NOTHING moved (it may be submitted again once the cause is removed);
#       3 refused (merged/ exists); 1 the merge failed during its execution (merged/ has no .merge_complete) or p2_reduce.py failed.
#   cd ~/harddisks_gen3_p2/hspist3 && sbatch cluster/gen3_stage2_261016/p2_merge.sbatch
#SBATCH --job-name=gen3-p2-merge
#SBATCH --partition=shared
#SBATCH --account=uh
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --exclude=cn-05-08-07
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.out
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:?submit with sbatch from ~/harddisks_gen3_p2/hspist3}"
P=cluster/gen3_stage2_261016; S1=cluster/gen3_stage1_261012; GRP=P2
ROOT_NAME=$(basename "$(dirname "$SLURM_SUBMIT_DIR")")
[ "$ROOT_NAME" = harddisks_gen3_p2 ] || { echo "STOP: the P2 merge runs from ~/harddisks_gen3_p2, not $ROOT_NAME"; exit 2; }
sha256sum --status -c "$P/p2_inputs.sha256" || { echo "STOP: the inputs are not the delivered ones ($P/p2_inputs.sha256); NOTHING was moved"; exit 2; }
source cluster/koa_env.sh || { echo "STOP: cluster/koa_env.sh failed; NOTHING was moved"; exit 2; }
ROOT=/mnt/lustre/koa/scratch/charing/harddisks_gen3/stage2/P2
TASKS="$P/tasks_p2_koa_261016.txt"; CT="$P/chunks_p2_koa_261016.tsv"; PFX="gen3-p2"
[ "$(sha256sum "$TASKS" | cut -d' ' -f1)" = "$(cat "$P/tasks_p2_koa_261016.sha256")" ] || { echo "STOP: the task list is not the delivered one"; exit 2; }
NCH=$(awk -F'\t' -v g="$GRP" 'NR > 1 && $5 == g' "$CT" | wc -l | tr -d ' ')
[ -d "$ROOT/chunks" ] || { echo "STOP: no chunks under $ROOT"; exit 2; }
[ ! -e "$ROOT/merged" ] || { echo "STOP: $ROOT/merged exists -- a merge ran before; not merging twice"; exit 3; }
Q=$(squeue -h -u "$USER" -n "$PFX" -o '%i %T' 2>&1); qrc=$?
[ "$qrc" -eq 0 ] || { echo "STOP: squeue failed (exit $qrc): $Q"; exit 2; }
[ -z "$Q" ] || { echo "STOP: $GRP chunks are still pending or running: $(printf '%s\n' "$Q" | head -3 | tr '\n' ' ')"; exit 2; }
echo "merge of $ROOT by the scripts of ~/harddisks_gen3_p2 (no git at run time; the input manifest: OK); the scripts it runs:"
sha256sum "$S1/merge_s1_koa.py" "$S1/s1_chunk_logs_koa.sh" "$S1/reduce_s1_koa.py" cluster/confinement_20261013/reduce_B.py "$P/p2_reduce.py" "$P/p2_lattice.py" | sed 's/^/  /'
echo "group $GRP: $NCH chunks of the chunk table; no $PFX job pending or running"
mkdir -p "$ROOT/job_logs"; cp -n -p logs/"$PFX"_*.out "$ROOT/job_logs/"
echo "chunk logs in $ROOT/job_logs/ (copied, no overwrite): $(ls "$ROOT/job_logs" | grep -c "^${PFX}_")"
bash "$S1/s1_chunk_logs_koa.sh" logs "$PFX" "$ROOT/chunks" "$CT" "$GRP" \
  || { echo "STOP: the chunk set is not complete (listed above); NOTHING was moved"; exit 2; }
NF=$(cat logs/"$PFX"_*.out | grep -c FAILED); echo "FAILED lines in the chunk logs (logs/): $NF"
[ "$NF" -eq 0 ] || grep -H FAILED logs/"$PFX"_*.out | head -100
python3 "$S1/merge_s1_koa.py" "$ROOT" --tasks "$TASKS" --chunks "$CT" --group "$GRP"; rc=$?
case $rc in
  0) echo "merge done; next: p2_reduce.py" ;;
  2) echo "STOP: the merge plan has problems (listed above); NOTHING was moved; once the cause is removed this job may be submitted again"; exit 2 ;;
  3) echo "STOP: refused (merged/ exists)"; exit 3 ;;
  *) echo "STOP: the merge failed during its execution (exit $rc): merged/ has no .merge_complete; inspect by hand, do not submit again"; exit 1 ;;
esac
python3 "$P/p2_reduce.py" "$ROOT/merged" --tasks "$TASKS"; rc=$?
[ "$rc" -eq 0 ] || { echo "STOP: p2_reduce.py failed (exit $rc); the merge itself is complete (.merge_complete); paste this log"; exit 1; }
echo "P2 merge and reduction done; next: the fetch from the Mac (fetch_p2_koa.sh)"
echo "this job's log is copied to $ROOT/job_logs/ as the last command"
cp -n -p "logs/${SLURM_JOB_NAME:-gen3-p2-merge}_${SLURM_JOB_ID}.out" "$ROOT/job_logs/"
exit 0
