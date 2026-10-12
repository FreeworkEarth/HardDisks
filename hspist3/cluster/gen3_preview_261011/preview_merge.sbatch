#!/bin/bash
# ##CHRIS 2026-10-11 (plan-author decision 29, task P; 261012 sec. 4.7.59): THE MERGE OF THE PREVIEW SWEEP ON KOA (PREVIEW, INFORMATION),
# one process, from the delivery clone ~/harddisks_gen3_k2 (updated in place to the preview commit), after every chunk has ended.
# THIS FILE IS cluster/gen3_testR_261011/testR_merge.sbatch (b5378f4) WITH THESE CHANGES, and no other:
#   - the job name gen3-preview-merge; the package folder P=cluster/gen3_preview_261011, its two input manifests (preview_inputs.sha256,
#     preview_merge_inputs.sha256), its task list and SHA-256 file, its chunk size file;
#   - the data root /mnt/lustre/koa/scratch/charing/harddisks_gen3/preview_koa;
#   - the job name it waits for in squeue and the chunk check's log prefix: gen3-preview (Test R: gen3-testR);
#   - the chunk check, the merge and the KOA reduction are Test R's committed scripts, run UNCHANGED from cluster/gen3_testR_261011/
#     (chunk_logs_koa.sh, merge_koa.py, and through it reduce_B.py and reduce_testR_koa.py); preview_merge_inputs.sha256 holds them;
#   - the job logs copied to <root>/job_logs/ have the prefix gen3-preview; the words "Test R" in the texts.
# Its rules are Test R's (fd8d8a2's):
#   - it does not depend on the array's job id: Chris submits it when squeue lists no gen3-preview job; it checks that itself;
#   - every check comes before any move, and a failed check stops it with exit 2 and NOTHING moved: the inputs as delivered (both
#     manifests, sha256sum -c), the task list as registered, no gen3-preview job pending or running, the chunk check (every chunk
#     1..<count> has its folder and EXACTLY ONE log logs/gen3-preview_*_<j>.out with its line "chunk <j> done");
#   - merged/ must not exist (exit 3): a merge is never repeated;
#   - then merge_koa.py: the plan, and only a clean plan is executed (merged/.merge_complete), then the reduction (reduce_B.py: E0's
#     per-trajectory frequency into red_nu.csv; reduce_testR_koa.py: the health and nu_d tables) and the compression of the traces.
# exit: 0 merged; 2 a check or the merge plan stopped it, NOTHING moved; 3 refused (merged/ exists); 1 the merge failed during its
# execution: merged/ has no .merge_complete; inspect by hand.
#   cd ~/harddisks_gen3_k2/hspist3 && sbatch cluster/gen3_preview_261011/preview_merge.sbatch
#SBATCH --job-name=gen3-preview-merge
#SBATCH --partition=shared
#SBATCH --account=uh
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --exclude=cn-05-08-07
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.out
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:?submit with sbatch from ~/harddisks_gen3_k2/hspist3}"
P=cluster/gen3_preview_261011; R=cluster/gen3_testR_261011
ROOT_NAME=$(basename "$(dirname "$SLURM_SUBMIT_DIR")")
[ "$ROOT_NAME" = harddisks_gen3_k2 ] || { echo "STOP: the preview merge runs from ~/harddisks_gen3_k2, not $ROOT_NAME"; exit 2; }
sha256sum --status -c "$P/preview_inputs.sha256" && sha256sum --status -c "$P/preview_merge_inputs.sha256" \
  || { echo "STOP: the inputs are not the delivered ones (preview_inputs.sha256, preview_merge_inputs.sha256); NOTHING was moved"; exit 2; }
source cluster/koa_env.sh || { echo "STOP: cluster/koa_env.sh failed; NOTHING was moved"; exit 2; }
ROOT=/mnt/lustre/koa/scratch/charing/harddisks_gen3/preview_koa
TASKS="$P/tasks_preview_koa_261011.txt"
[ "$(sha256sum "$TASKS" | cut -d' ' -f1)" = "$(cat "$P/tasks_preview_koa_261011.sha256")" ] || { echo "STOP: the task list is not the registered one"; exit 2; }
K=$(cat "$P/preview_chunk_size.txt"); case "$K" in ''|*[!0-9]*) echo "STOP: the chunk size '$K' is not a number"; exit 2 ;; esac
NL=$(wc -l < "$TASKS"); NCH=$(( (NL + K - 1) / K ))
[ -d "$ROOT/chunks" ] || { echo "STOP: no chunks under $ROOT"; exit 2; }
[ ! -e "$ROOT/merged" ] || { echo "STOP: $ROOT/merged exists -- a merge ran before; not merging twice"; exit 3; }
Q=$(squeue -h -u "$USER" -n gen3-preview -o '%i %T' 2>&1); qrc=$?
[ "$qrc" -eq 0 ] || { echo "STOP: squeue failed (exit $qrc): $Q"; exit 2; }
[ -z "$Q" ] || { echo "STOP: preview chunks are still pending or running: $(printf '%s\n' "$Q" | head -3 | tr '\n' ' ')"; exit 2; }
echo "merge of $ROOT by the scripts of ~/harddisks_gen3_k2 (no git at run time; both input manifests: OK); the scripts it runs:"
sha256sum "$R/merge_koa.py" "$R/chunk_logs_koa.sh" "$R/reduce_testR_koa.py" cluster/confinement_20261013/reduce_B.py | sed 's/^/  /'
echo "the registered list: $NL task lines, chunk size $K, $NCH chunks; no gen3-preview job pending or running"
mkdir -p "$ROOT/job_logs"; cp -n -p logs/gen3-preview_*.out "$ROOT/job_logs/"   # cp -n: its exit status on a skipped file differs between
echo "chunk logs in $ROOT/job_logs/ (copied, no overwrite): $(ls "$ROOT/job_logs" | grep -c '^gen3-preview_')"   # coreutils versions
bash "$R/chunk_logs_koa.sh" logs gen3-preview "$ROOT/chunks" "$NCH" \
  || { echo "STOP: the chunk set is not complete (listed above); NOTHING was moved"; exit 2; }
NF=$(cat logs/gen3-preview_*.out | grep -c FAILED); echo "FAILED lines in the chunk logs (logs/): $NF"
[ "$NF" -eq 0 ] || grep -H FAILED logs/gen3-preview_*.out | head -100
python3 "$R/merge_koa.py" "$ROOT" --tasks "$TASKS" --chunk-size "$K"; rc=$?
case $rc in
  0) echo "merge done; next: the fetch (fetch_preview_koa.sh) and the E0 reduction on the Mac" ;;
  2) echo "STOP: the merge plan has problems (listed above); NOTHING was moved; once the cause is removed this job may be submitted again"; exit 2 ;;
  3) echo "STOP: refused (merged/ exists)"; exit 3 ;;
  *) echo "STOP: the merge failed during its execution (exit $rc): merged/ has no .merge_complete; inspect by hand, do not submit again"; exit 1 ;;
esac
echo "this job's log is copied to $ROOT/job_logs/ as the last command"
cp -n -p "logs/gen3-preview-merge_${SLURM_JOB_ID}.out" "$ROOT/job_logs/"
exit 0
