#!/bin/bash
# ##CHRIS 2026-10-15 (plan-author decision 57; notes 261012 sec. 4.7.135): THE STAGE 2 PILOT P2 ON KOA, the array of one-core chunks
# (one task line per chunk, 400 chunks). It is cluster/gen3_stage1_261012/s1_koa.sbatch with THESE CHANGES, all of them:
#   - the clone: ~/harddisks_gen3_p2 (one fresh clone: built at 8aaf5fef by cluster/build_koa.sh, then checked out at the package
#     commit, whose build inputs are 8aaf5fef's; runsheet block a), so the jobs run from ~/harddisks_gen3_p2/hspist3;
#   - the inputs: p2_inputs.sha256; the task list tasks_p2_koa_261016.txt, its SHA-256 file and the chunk table chunks_p2_koa_261016.tsv;
#   - the binary ~/harddisks_gen3_p2/hspist3/00ALLINONE, checked WITHOUT git against the build-time record logs/BUILD_RECORD_P2.txt that
#     record_p2_build.sh wrote (the build commit 8aaf5fef, the version line, the binary's SHA-256);
#   - the worker p2_worker_koa.sh (the stage 1 worker plus the initial dump); the data root .../harddisks_gen3/stage2/P2/chunks/chunk_<j>;
#   - the time limit 0:25:00 (decision 57: the D56-1 table's longest line on the slower E5-4640 v2 at the upper Z, times the Mac
#     rehearsal's measured ratio for the pin 0 crystals, 1.171 -- 14.1 min -- + 50 %, rounded up to 5 min; gen_p2_koa_261016.py).
# NO GIT AT RUN TIME (node cn-05-08-07 has none and is excluded). Every start-up check comes BEFORE the line "chunk <j> on ...", so a
# chunk that stops in them ran nothing; every STOP line names its node: "STOP: [<node>] ...".
# Submit from ~/harddisks_gen3_p2/hspist3 (runsheet block b):  sbatch --array=1-400%64 cluster/gen3_stage2_261016/p2_koa.sbatch
#SBATCH --job-name=gen3-p2
#SBATCH --partition=shared
#SBATCH --account=uh
#SBATCH --time=0:25:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --exclude=cn-05-08-07
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.out
set -uo pipefail
stop() { echo "STOP: [$(hostname)] $*"; exit 2; }
cd "${SLURM_SUBMIT_DIR:?submit with sbatch from ~/harddisks_gen3_p2/hspist3}"
P=cluster/gen3_stage2_261016
ROOT_NAME=$(basename "$(dirname "$SLURM_SUBMIT_DIR")")
[ "$ROOT_NAME" = harddisks_gen3_p2 ] || stop "P2 runs from ~/harddisks_gen3_p2, not $ROOT_NAME"
sha256sum --status -c "$P/p2_inputs.sha256" || stop "the P2 inputs are not the delivered ones ($P/p2_inputs.sha256)"
TASKS="$P/tasks_p2_koa_261016.txt"; CT="$P/chunks_p2_koa_261016.tsv"
[ "$(sha256sum "$TASKS" | cut -d' ' -f1)" = "$(cat "$P/tasks_p2_koa_261016.sha256")" ] || stop "the task list is not the delivered one"
source cluster/koa_env.sh || stop "cluster/koa_env.sh failed"
BIN="$SLURM_SUBMIT_DIR/00ALLINONE"; BREC="$SLURM_SUBMIT_DIR/logs/BUILD_RECORD_P2.txt"
[ -r "$BREC" ] || stop "no build-time record $BREC"
grep -qx "commit_build 8aaf5fef27c92391fc4fadd5133324740c29ef06" "$BREC" || stop "the build-time record does not name the build commit 8aaf5fef"
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
[ "$grp" = P2 ] || stop "chunk $j belongs to group '$grp', not P2"
[ "$first" -ge 1 ] && [ "$last" -ge "$first" ] && [ "$last" -le "$(wc -l < "$TASKS")" ] || stop "chunk $j's lines $first-$last are not within the task list"
export HD_BIN="$BIN" HD_BUILD
export HD_DATA="/mnt/lustre/koa/scratch/charing/harddisks_gen3/stage2/P2/chunks/chunk_$j"; mkdir -p "$HD_DATA"
echo "chunk $j on $(hostname) ($(lscpu | awk -F: '/Model name/{gsub(/^ +/,"",$2); print $2; exit}')): group $grp, lines $first-$last; $HD_BUILD"
sed -n "${first},${last}p" "$TASKS" | while read -r line; do
  set -- $line
  bash "$P/p2_worker_koa.sh" "$@"
done
echo "chunk $j done"
