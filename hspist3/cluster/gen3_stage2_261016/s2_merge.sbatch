#!/bin/bash
# ##CHRIS 2026-10-16 (plan-author decision 58; registration S2, notes 261012 sec. 4.9 and sec. 4.7.142): THE MERGE OF STAGE 2 "S2" ON KOA, one
# process, from ~/harddisks_gen3_s2/hspist3. It is p2_merge.sbatch (P2's merge: the stage 1 merge and reducers, unchanged, then p2_reduce.py)
# with THESE CHANGES, all of them: the group S2, its root .../harddisks_gen3/stage2/S2, its job names gen3-s2 / gen3-s2-merge, its task list,
# chunk table and input manifest, the clone ~/harddisks_gen3_s2; and, since decision 58 submits it in the SAME paste as the array, held by
# Slurm (--dependency=afterok:<the array's job id>), the queue check WAITS: up to 10 min (20 checks, 30 s apart) for squeue to list no gen3-s2
# task (a task may linger a few seconds in COMPLETING after the dependency is met), and only then stops. Every other rule is P2's:
#   - every check comes before any move, and a failed check stops it with exit 2 and NOTHING moved: the inputs as delivered, the task
#     list as delivered, no gen3-s2 task pending or running, the chunk check (every chunk has its folder and EXACTLY ONE log
#     logs/gen3-s2_*_<j>.out with its line "chunk <j> done");
#   - merged/ must not exist (exit 3): a merge is never repeated;
#   - the job logs are copied to <root>/job_logs/ (no overwrite).
# If an array task ends with a non-zero exit (a start-up STOP, a time-out), Slurm never starts this job (squeue: "DependencyNeverSatisfied"):
# nothing is merged; the runsheet's read-only checks show it.
# exit: 0 merged and reduced; 2 a check or the merge plan stopped it, NOTHING moved; 3 refused (merged/ exists); 1 the merge failed during its
#       execution (merged/ has no .merge_complete) or p2_reduce.py failed.
#SBATCH --job-name=gen3-s2-merge
#SBATCH --partition=shared
#SBATCH --account=uh
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --exclude=cn-05-08-07
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.out
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:?submit with sbatch from ~/harddisks_gen3_s2/hspist3}"
P=cluster/gen3_stage2_261016; S1=cluster/gen3_stage1_261012; GRP=S2
ROOT_NAME=$(basename "$(dirname "$SLURM_SUBMIT_DIR")")
[ "$ROOT_NAME" = harddisks_gen3_s2 ] || { echo "STOP: the S2 merge runs from ~/harddisks_gen3_s2, not $ROOT_NAME"; exit 2; }
sha256sum --status -c "$P/s2_inputs.sha256" || { echo "STOP: the inputs are not the delivered ones ($P/s2_inputs.sha256); NOTHING was moved"; exit 2; }
source cluster/koa_env.sh || { echo "STOP: cluster/koa_env.sh failed; NOTHING was moved"; exit 2; }
ROOT=/mnt/lustre/koa/scratch/charing/harddisks_gen3/stage2/S2
TASKS="$P/tasks_s2_koa_261017.txt"; CT="$P/chunks_s2_koa_261017.tsv"; PFX="gen3-s2"
[ "$(sha256sum "$TASKS" | cut -d' ' -f1)" = "$(cat "$P/tasks_s2_koa_261017.sha256")" ] || { echo "STOP: the task list is not the delivered one"; exit 2; }
NCH=$(awk -F'\t' -v g="$GRP" 'NR > 1 && $5 == g' "$CT" | wc -l | tr -d ' ')
[ -d "$ROOT/chunks" ] || { echo "STOP: no chunks under $ROOT"; exit 2; }
[ ! -e "$ROOT/merged" ] || { echo "STOP: $ROOT/merged exists -- a merge ran before; not merging twice"; exit 3; }
for w in $(seq 1 20); do   # the dependency is met when every task has ended; a task may still be COMPLETING for a few seconds
  Q=$(squeue -h -u "$USER" -n "$PFX" -o '%i %T' 2>&1); qrc=$?
  [ "$qrc" -eq 0 ] || { echo "STOP: squeue failed (exit $qrc): $Q"; exit 2; }
  [ -z "$Q" ] && break
  [ "$w" -lt 20 ] && { echo "waiting for $(printf '%s\n' "$Q" | wc -l | tr -d ' ') $PFX task(s) to leave the queue ($w of 20, 30 s)"; sleep 30; }
done
[ -z "$Q" ] || { echo "STOP: $GRP chunks are still pending or running after 10 min: $(printf '%s\n' "$Q" | head -3 | tr '\n' ' ')"; exit 2; }
echo "merge of $ROOT by the scripts of ~/harddisks_gen3_s2 (no git at run time; the input manifest: OK); the scripts it runs:"
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
echo "S2 merge and reduction done; next: the fetch from the Mac (fetch_s2_koa.sh)"
echo "this job's log is copied to $ROOT/job_logs/ as the last command"
cp -n -p "logs/${SLURM_JOB_NAME:-gen3-s2-merge}_${SLURM_JOB_ID}.out" "$ROOT/job_logs/"
exit 0
