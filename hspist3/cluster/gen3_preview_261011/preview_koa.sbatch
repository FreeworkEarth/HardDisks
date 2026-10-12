#!/bin/bash
# ##CHRIS 2026-10-11 (plan-author decision 29, task P; 261012 sec. 4.7.59): THE PREVIEW SWEEP ON KOA -- PREVIEW, INFORMATION; no verdict,
# never part of stage 1, never in a paper; its own root and its own seeds. The array of one-core chunks. PREPARED; typed only after the
# plan author's go on Report P0 and after Test R's array has been submitted. Chunk j runs the task lines (j-1)*K+1 .. j*K of
# tasks_preview_koa_261011.txt (K from preview_chunk_size.txt), in order, one at a time, with the B command of
# cluster/gen3_koa_261009/testG_worker_koa.sh (unchanged since 1cdb7ff6), writing ONLY into its own chunk folder $DATA/chunks/chunk_<j>.
# The merge job (preview_merge.sbatch) assembles the cells after every chunk has ended.
# THIS FILE IS cluster/gen3_testR_261011/testR_koa.sbatch (b5378f4) WITH THESE CHANGES, and no other:
#   - the job name gen3-preview (Test R: gen3-testR), so the chunk logs are logs/gen3-preview_<job>_<task>.out;
#   - the package folder P=cluster/gen3_preview_261011, its input manifest preview_inputs.sha256, its task list
#     tasks_preview_koa_261011.txt with its registered SHA-256 file, its chunk size file preview_chunk_size.txt;
#   - the build record: Test R's committed cluster/gen3_testR_261011/build_record_a643bc7.txt, read unchanged (the same binary, not
#     rebuilt; it is in preview_inputs.sha256);
#   - the data root /mnt/lustre/koa/scratch/charing/harddisks_gen3/preview_koa (Test R: .../testR_koa);
#   - --time (from the design: gen_preview_design_output.txt);
#   - the words "Test R" in the STOP texts and in this header.
# NO GIT AT RUN TIME. Every start-up check comes BEFORE the line "chunk <j> on ...", so a chunk that stops in them ran nothing, and
# every STOP line names its node: "STOP: [<node>] ...". The checks are Test R's:
#   1. the clone: this job runs from ~/harddisks_gen3_k2/hspist3 (the K2 clone, updated in place to the preview commit by Report P0's
#      delivery block);
#   2. the inputs (sha256sum -c preview_inputs.sha256): the task list, its SHA-256 file, the chunk size, the build record, the worker
#      and koa_env.sh are the delivered ones; and the task list's SHA-256 is the registered one;
#   3. koa_env.sh (which needs no git);
#   4. the binary ~/harddisks_gen3/hspist3/00ALLINONE, without git: its SHA-256 equals the build record's AND the build-time file
#      ~/harddisks_gen3/hspist3/logs/BUILD_KOA_LAST.hash, and its version line is the recorded "00ALLINONE  git a643bc7  target koa";
#   5. the chunk size is a number and the chunk index lies within the task list.
# Submit from ~/harddisks_gen3_k2/hspist3 (runsheet_preview_261011.txt): chunk 1 alone, then the rest.
#SBATCH --job-name=gen3-preview
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
cd "${SLURM_SUBMIT_DIR:?submit with sbatch from ~/harddisks_gen3_k2/hspist3}"
P=cluster/gen3_preview_261011
ROOT_NAME=$(basename "$(dirname "$SLURM_SUBMIT_DIR")")
[ "$ROOT_NAME" = harddisks_gen3_k2 ] || stop "the preview runs from ~/harddisks_gen3_k2, not $ROOT_NAME"
sha256sum --status -c "$P/preview_inputs.sha256" || stop "the preview inputs are not the delivered ones ($P/preview_inputs.sha256)"
TASKS="$P/tasks_preview_koa_261011.txt"
[ "$(sha256sum "$TASKS" | cut -d' ' -f1)" = "$(cat "$P/tasks_preview_koa_261011.sha256")" ] || stop "the task list is not the registered one"
source cluster/koa_env.sh || stop "cluster/koa_env.sh failed"
REC="cluster/gen3_testR_261011/build_record_a643bc7.txt"; A="$HOME/harddisks_gen3/hspist3"; BIN="$A/00ALLINONE"
want_sha=$(awk '$1 == "sha256" {print $2}' "$REC"); want_ver=$(sed -n 's/^version //p' "$REC")
[ -n "$want_sha" ] && [ -n "$want_ver" ] || stop "$REC lacks its sha256 or version line"
have_sha=$(sha256sum "$BIN" | cut -d' ' -f1)
[ "$have_sha" = "$want_sha" ] || stop "$BIN is not the recorded binary (SHA-256 $have_sha)"
grep -qx "$want_sha  00ALLINONE" "$A/logs/BUILD_KOA_LAST.hash" || stop "the build-time file $A/logs/BUILD_KOA_LAST.hash does not hold the recorded SHA-256"
HD_BUILD=$("$BIN" --version | head -1)
[ "$HD_BUILD" = "$want_ver" ] || stop "the binary's version line is '$HD_BUILD', recorded '$want_ver'"
K=$(cat "$P/preview_chunk_size.txt"); case "$K" in ''|*[!0-9]*) stop "the chunk size '$K' is not a number" ;; esac
j=${SLURM_ARRAY_TASK_ID:?submit as an array}
[ "$j" -ge 1 ] && [ $(( (j-1)*K )) -lt "$(wc -l < "$TASKS")" ] || stop "chunk $j is beyond the task list"
export HD_BIN="$BIN" HD_BUILD
export HD_DATA="/mnt/lustre/koa/scratch/charing/harddisks_gen3/preview_koa/chunks/chunk_$j"; mkdir -p "$HD_DATA"
echo "chunk $j on $(hostname) ($(lscpu | awk -F: '/Model name/{gsub(/^ +/,"",$2); print $2; exit}')): lines $(( (j-1)*K+1 ))-$(( j*K )); $HD_BUILD"
sed -n "$(( (j-1)*K+1 )),$(( j*K ))p" "$TASKS" | while read -r line; do
  set -- $line
  bash cluster/gen3_koa_261009/testG_worker_koa.sh "${@:1:11}"
done
echo "chunk $j done"
