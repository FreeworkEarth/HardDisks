#!/bin/bash
# ##CHRIS 2026-10-13 (plan-author decision 47 item 3; the contingency for D's merge, job 15055588): THE REDUCTION COMPLETION OF ONE STAGE 1
# GROUP whose merge job stopped during its reduction (killed at its time limit after merged/.merge_complete). It runs from the clone
# ~/harddisks_gen3_rc (the engine-gen3 commit that adds this folder; the stage 1 package's files there are byte for byte the package's:
# s1_inputs.sha256 and s1_merge_inputs.sha256 pass) and runs complete_s1_koa.py, which runs the package's own reduce_B.py and
# reduce_s1_koa.py, the merge's SHA-256 lines and gzip -9 on the folders the merge did not finish (its description says how). Every
# check comes before any change; a failed check stops it with exit 2 and NOTHING changed:
#   - the group; the clone; the three manifests (s1_inputs.sha256, s1_merge_inputs.sha256, this job's own s1_complete_inputs.sha256: every
#     file it runs); cluster/koa_env.sh;
#   - no gen3-s1-merge-<group> job and no other gen3-s1-complete-<group> job pending or running (a failed squeue stops it);
#   - <root>/merged/.merge_complete exists; <root>/merged/.reduction_completed does not (exit 3: the completion ran before to the end);
#   - the merge log, the second argument, exists and is the log of that group's merge ("merge of <root> by the scripts of ...");
#   - complete_s1_koa.py --plan-only: the plan, read only (exit 2 for a state it does not recognise).
# Then complete_s1_koa.py; the merge log and, as the last command, this job's log are copied to <root>/job_logs/ (no overwrite).
#   cd ~/harddisks_gen3_rc/hspist3 && sbatch --time=24:00:00 --job-name=gen3-s1-complete-<group> cluster/gen3_s1_complete_261013/s1_complete.sbatch <group> <merge log>
# exit: 0 completed (or nothing left to do); 2 a check stopped it, NOTHING changed; 3 refused (it ran before); 1 a step failed (inspect)
#SBATCH --job-name=gen3-s1-complete
#SBATCH --partition=shared
#SBATCH --account=uh
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --exclude=cn-05-08-07
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.out
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:?submit with sbatch from ~/harddisks_gen3_rc/hspist3}"
P=cluster/gen3_s1_complete_261013; S1=cluster/gen3_stage1_261012
GRP=${1:?the group (B100, B400, D, B900, S, B1600) is the first argument}; MLOG=${2:?the merge log is the second argument}
case "$GRP" in B100|B400|D|B900|S|B1600) ;; *) echo "STOP: '$GRP' is not a stage 1 group; NOTHING was changed"; exit 2 ;; esac
ROOT_NAME=$(basename "$(dirname "$SLURM_SUBMIT_DIR")")
[ "$ROOT_NAME" = harddisks_gen3_rc ] || { echo "STOP: the completion runs from ~/harddisks_gen3_rc, not $ROOT_NAME; NOTHING was changed"; exit 2; }
sha256sum --status -c "$S1/s1_inputs.sha256" && sha256sum --status -c "$S1/s1_merge_inputs.sha256" && sha256sum --status -c "$P/s1_complete_inputs.sha256" \
  || { echo "STOP: the inputs are not the delivered ones (s1_inputs.sha256, s1_merge_inputs.sha256, s1_complete_inputs.sha256); NOTHING was changed"; exit 2; }
source cluster/koa_env.sh || { echo "STOP: cluster/koa_env.sh failed; NOTHING was changed"; exit 2; }
ROOT=/mnt/lustre/koa/scratch/charing/harddisks_gen3/stage1/$GRP
Q=$(squeue -h -u "$USER" -n "gen3-s1-merge-$GRP,gen3-s1-complete-$GRP" -o '%i %j %T' 2>&1); qrc=$?
[ "$qrc" -eq 0 ] || { echo "STOP: squeue failed (exit $qrc): $Q; NOTHING was changed"; exit 2; }
Q=$(printf '%s\n' "$Q" | awk -v me="${SLURM_JOB_ID:-none}" 'NF && $1 != me')
[ -z "$Q" ] || { echo "STOP: a merge or another completion of $GRP is pending or running: $(printf '%s\n' "$Q" | head -3 | tr '\n' ' '); NOTHING was changed"; exit 2; }
[ -e "$ROOT/merged/.merge_complete" ] || { echo "STOP: $ROOT/merged/.merge_complete is missing: the merge did not finish its moves; NOTHING was changed"; exit 2; }
[ ! -e "$ROOT/merged/.reduction_completed" ] || { echo "STOP: $ROOT/merged/.reduction_completed exists: the completion ran before to the end"; exit 3; }
[ -f "$MLOG" ] && grep -q "^merge of $ROOT by the scripts of " "$MLOG" \
  || { echo "STOP: $MLOG is not the log of the merge of $ROOT; NOTHING was changed"; exit 2; }
echo "completion of the reduction of $ROOT by the scripts of ~/harddisks_gen3_rc (no git at run time; the three manifests: OK); the merge log $MLOG; the scripts it runs:"
sha256sum "$P/complete_s1_koa.py" "$S1/reduce_s1_koa.py" cluster/confinement_20261013/reduce_B.py cluster/gen3_testR_261011/reduce_testR_koa.py | sed 's/^/  /'
echo "the merge log's last line: $(tail -n 1 "$MLOG")"
python3 "$P/complete_s1_koa.py" "$ROOT" --plan-only; rc=$?
[ "$rc" -eq 0 ] || { echo "STOP: the plan has problems (listed above); NOTHING was changed"; exit 2; }
python3 "$P/complete_s1_koa.py" "$ROOT"; rc=$?
mkdir -p "$ROOT/job_logs"; cp -n -p "$MLOG" "$ROOT/job_logs/"
case $rc in
  0) echo "completion done; next: block H-<group>'s counts from this job's log and the fetch" ;;
  3) echo "STOP: refused (it ran before)"; exit 3 ;;
  2) echo "STOP: refused; NOTHING was changed"; exit 2 ;;
  *) echo "STOP: a step failed (exit $rc): merged/.reduction_completed is not written; inspect by hand, do not submit again"; exit 1 ;;
esac
echo "this job's log is copied to $ROOT/job_logs/ as the last command"
cp -n -p "logs/${SLURM_JOB_NAME:-gen3-s1-complete}_${SLURM_JOB_ID}.out" "$ROOT/job_logs/"
exit 0
