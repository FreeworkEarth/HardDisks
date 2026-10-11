#!/bin/bash
# ##CHRIS 2026-10-11 (plan-author decision 28; 261012 sec. 4.7.57): TEST R ON KOA, the array of one-core chunks. PREPARED; typed only after
# the plan author has read Report R0. Chunk j runs the task lines (j-1)*K+1 .. j*K of tasks_testR_koa_261011.txt (K from
# testR_chunk_size.txt: 100, so 21,600 lines = 216 chunks), in order, one at a time, with the B command of
# cluster/gen3_koa_261009/testG_worker_koa.sh (unchanged since 1cdb7ff6), writing ONLY into its own chunk folder
# $DATA/chunks/chunk_<j>. The merge job (testR_merge.sbatch) assembles the cells after every chunk has ended.
# NO GIT AT RUN TIME (decisions 26 and 28: node cn-05-08-07 has none). Every start-up check comes BEFORE the line "chunk <j> on ...", so
# a chunk that stops in them ran nothing (the chunk check lists such an attempt by name; decision 25's rule: it is submitted once more,
# testR_stopped_chunks.sh prints that line). Every STOP line names the node it ran on: "STOP: [<node>] ...".
#   1. the clone: this job runs from ~/harddisks_gen3_r/hspist3 (the delivery clone of Report R0);
#   2. the inputs: testR_inputs.sha256 (sha256sum -c): the task list, its SHA-256 file, the chunk size, the build record, the worker
#      and koa_env.sh are the delivered ones; and the task list's SHA-256 is the registered one (tasks_testR_koa_261011.sha256);
#   3. koa_env.sh (which no longer needs git);
#   4. the binary ~/harddisks_gen3/hspist3/00ALLINONE, without git: its SHA-256 equals build_record_a643bc7.txt's AND the build-time file
#      ~/harddisks_gen3/hspist3/logs/BUILD_KOA_LAST.hash (written by build_gen3_koa.sh on KOA in runsheet step 2), and its version line
#      is the recorded "00ALLINONE  git a643bc7  target koa";
#   5. the chunk size is a number and the chunk index lies within the task list.
# Submit from ~/harddisks_gen3_r/hspist3 (runsheet_testR_261011.txt): chunk 1 alone, then the rest:
#   sbatch --array=1 cluster/gen3_testR_261011/testR_koa.sbatch
#   sbatch --array=2-216%64 cluster/gen3_testR_261011/testR_koa.sbatch
#SBATCH --job-name=gen3-testR
#SBATCH --partition=shared
#SBATCH --account=uh
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --exclude=cn-05-08-07
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.out
set -uo pipefail
stop() { echo "STOP: [$(hostname)] $*"; exit 2; }
cd "${SLURM_SUBMIT_DIR:?submit with sbatch from ~/harddisks_gen3_r/hspist3}"
P=cluster/gen3_testR_261011
ROOT_NAME=$(basename "$(dirname "$SLURM_SUBMIT_DIR")")
[ "$ROOT_NAME" = harddisks_gen3_r ] || stop "Test R runs from ~/harddisks_gen3_r, not $ROOT_NAME"
sha256sum --status -c "$P/testR_inputs.sha256" || stop "the Test R inputs are not the delivered ones ($P/testR_inputs.sha256)"
TASKS="$P/tasks_testR_koa_261011.txt"
[ "$(sha256sum "$TASKS" | cut -d' ' -f1)" = "$(cat "$P/tasks_testR_koa_261011.sha256")" ] || stop "the task list is not the registered one"
source cluster/koa_env.sh || stop "cluster/koa_env.sh failed"
REC="$P/build_record_a643bc7.txt"; A="$HOME/harddisks_gen3/hspist3"; BIN="$A/00ALLINONE"
want_sha=$(awk '$1 == "sha256" {print $2}' "$REC"); want_ver=$(sed -n 's/^version //p' "$REC")
[ -n "$want_sha" ] && [ -n "$want_ver" ] || stop "$REC lacks its sha256 or version line"
have_sha=$(sha256sum "$BIN" | cut -d' ' -f1)
[ "$have_sha" = "$want_sha" ] || stop "$BIN is not the recorded binary (SHA-256 $have_sha)"
grep -qx "$want_sha  00ALLINONE" "$A/logs/BUILD_KOA_LAST.hash" || stop "the build-time file $A/logs/BUILD_KOA_LAST.hash does not hold the recorded SHA-256"
HD_BUILD=$("$BIN" --version | head -1)
[ "$HD_BUILD" = "$want_ver" ] || stop "the binary's version line is '$HD_BUILD', recorded '$want_ver'"
K=$(cat "$P/testR_chunk_size.txt"); case "$K" in ''|*[!0-9]*) stop "the chunk size '$K' is not a number" ;; esac
j=${SLURM_ARRAY_TASK_ID:?submit as an array}
[ "$j" -ge 1 ] && [ $(( (j-1)*K )) -lt "$(wc -l < "$TASKS")" ] || stop "chunk $j is beyond the task list"
export HD_BIN="$BIN" HD_BUILD
export HD_DATA="/mnt/lustre/koa/scratch/charing/harddisks_gen3/testR_koa/chunks/chunk_$j"; mkdir -p "$HD_DATA"
echo "chunk $j on $(hostname) ($(lscpu | awk -F: '/Model name/{gsub(/^ +/,"",$2); print $2; exit}')): lines $(( (j-1)*K+1 ))-$(( j*K )); $HD_BUILD"
sed -n "$(( (j-1)*K+1 )),$(( j*K ))p" "$TASKS" | while read -r line; do
  set -- $line
  bash cluster/gen3_koa_261009/testG_worker_koa.sh "${@:1:11}"
done
echo "chunk $j done"
