#!/bin/bash
# ##CHRIS 2026-10-16 (plan-author decision 58; registration S2, notes 261012 sec. 4.9 and sec. 4.7.142): STAGE 2 "S2" ON KOA, the array of
# one-core chunks (one task line per chunk, 2100 chunks). It is p2_koa.sbatch (P2's array, decision 57) with THESE CHANGES, all of them:
#   - the clone ~/harddisks_gen3_s2 (one fresh clone: built at 8308fb00 -- the package commit's parent, pushed -- by its cluster/build_koa.sh,
#     then checked out at the package commit, which adds files only; the build inputs are 8308fb00's, equal to P2's 8aaf5fef's);
#   - the inputs s2_inputs.sha256; the task list tasks_s2_koa_261017.txt, its SHA-256 file and the chunk table chunks_s2_koa_261017.tsv;
#   - the build-time record logs/BUILD_RECORD_S2.txt (record_s2_build.sh), naming the build commit 8308fb00;
#   - the group S2; the data root .../harddisks_gen3/stage2/S2/chunks/chunk_<j>;
#   - the time limit 0:40:00 (decision 58: the longest line on E5-4640 v2 at the upper Z, 26.4 min, + 50 %, rounded up to 5 min).
# The worker is P2's p2_worker_koa.sh, unchanged (the stage 1 worker plus the initial dump). The merge (s2_merge.sbatch) is submitted in the
# same paste, held by Slurm until every task of this array has ended with exit 0 (--dependency=afterok:<this array's job id>).
# NO GIT AT RUN TIME (node cn-05-08-07 has none and is excluded). Every start-up check comes BEFORE the line "chunk <j> on ...", so a
# chunk that stops in them ran nothing; every STOP line names its node: "STOP: [<node>] ...".
# Submit from ~/harddisks_gen3_s2/hspist3 with the runsheet's paste B only (the array and the dependent merge in one guarded line).
#SBATCH --job-name=gen3-s2
#SBATCH --partition=shared
#SBATCH --account=uh
#SBATCH --time=0:40:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --exclude=cn-05-08-07
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.out
set -uo pipefail
stop() { echo "STOP: [$(hostname)] $*"; exit 2; }
cd "${SLURM_SUBMIT_DIR:?submit with sbatch from ~/harddisks_gen3_s2/hspist3}"
P=cluster/gen3_stage2_261016
ROOT_NAME=$(basename "$(dirname "$SLURM_SUBMIT_DIR")")
[ "$ROOT_NAME" = harddisks_gen3_s2 ] || stop "S2 runs from ~/harddisks_gen3_s2, not $ROOT_NAME"
sha256sum --status -c "$P/s2_inputs.sha256" || stop "the S2 inputs are not the delivered ones ($P/s2_inputs.sha256)"
TASKS="$P/tasks_s2_koa_261017.txt"; CT="$P/chunks_s2_koa_261017.tsv"
[ "$(sha256sum "$TASKS" | cut -d' ' -f1)" = "$(cat "$P/tasks_s2_koa_261017.sha256")" ] || stop "the task list is not the delivered one"
source cluster/koa_env.sh || stop "cluster/koa_env.sh failed"
BIN="$SLURM_SUBMIT_DIR/00ALLINONE"; BREC="$SLURM_SUBMIT_DIR/logs/BUILD_RECORD_S2.txt"
[ -r "$BREC" ] || stop "no build-time record $BREC"
grep -qx "commit_build 8308fb004fdcb7e9c82197451cc04ab49ecd2f13" "$BREC" || stop "the build-time record does not name the build commit 8308fb00"
want_sha=$(awk '$1 == "sha256" && $3 == "00ALLINONE" {print $2}' "$BREC"); want_ver=$(sed -n 's/^version //p' "$BREC")
[ -n "$want_sha" ] && [ -n "$want_ver" ] || stop "$BREC lacks its 00ALLINONE sha256 line or its version line"
have_sha=$(sha256sum "$BIN" | cut -d' ' -f1)
[ "$have_sha" = "$want_sha" ] || stop "$BIN is not the recorded binary (SHA-256 $have_sha)"
HD_BUILD=$("$BIN" --version | head -1)
[ "$HD_BUILD" = "$want_ver" ] || stop "the binary's version line is '$HD_BUILD', recorded '$want_ver'"
j=${SLURM_ARRAY_TASK_ID:?submit as an array}
row=$(awk -F'\t' -v j="$j" 'NR > 1 && $1 == j' "$CT")
[ -n "$row" ] || stop "chunk $j is not a row of the chunk table"
first=$(printf '%s\n' "$row" | cut -f2); last=$(printf '%s\n' "$row" | cut -f3); grp=$(printf '%s\n' "$row" | cut -f5)
[ "$grp" = S2 ] || stop "chunk $j belongs to group '$grp', not S2"
[ "$first" -ge 1 ] && [ "$last" -ge "$first" ] && [ "$last" -le "$(wc -l < "$TASKS")" ] || stop "chunk $j's lines $first-$last are not within the task list"
export HD_BIN="$BIN" HD_BUILD
export HD_DATA="/mnt/lustre/koa/scratch/charing/harddisks_gen3/stage2/S2/chunks/chunk_$j"; mkdir -p "$HD_DATA"
echo "chunk $j on $(hostname) ($(lscpu | awk -F: '/Model name/{gsub(/^ +/,"",$2); print $2; exit}')): group $grp, lines $first-$last; $HD_BUILD"
sed -n "${first},${last}p" "$TASKS" | while read -r line; do
  set -- $line
  bash "$P/p2_worker_koa.sh" "$@"
done
echo "chunk $j done"
