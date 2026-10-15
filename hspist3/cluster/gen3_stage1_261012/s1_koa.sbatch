#!/bin/bash
# ##CHRIS 2026-10-12 (plan-author decision 39 item 5 f; registration 261012 sec. 4.8.8; gate G1): STAGE 1 ON KOA, Part B, the array of
# one-core chunks. PREPARED; typed only after the plan author's go (gate G7). Chunk j runs the task lines <first>..<last> of the chunk
# table (chunks_stage1_koa_261012.tsv, row j) of tasks_stage1_koa_261012.txt, in order, one at a time, with s1_worker_koa.sh, writing
# ONLY into its own chunk folder $ROOT/<group>/chunks/chunk_<j>. The merge job (s1_merge.sbatch <group>) assembles a group's cells after
# every chunk of the group has ended. It is cluster/gen3_testR_261011/testR_koa.sbatch (decision 28) with THESE CHANGES, all of them:
#   - the task list, its SHA-256 and the CHUNK TABLE (the chunks are not of one size; the array index IS the chunk number j);
#   - the binary: the campaign build of f72d764 in ~/harddisks_gen3_s1 (decision 39 item 5 a), checked WITHOUT git against the committed
#     build record build_record_f72d764.txt (written from the build's printed values after the build) AND the build-time record
#     ~/harddisks_gen3_s1/hspist3/logs/BUILD_RECORD_S1.txt that the build wrote (commit, SHA-256 of every binary);
#   - the worker s1_worker_koa.sh (15-field lines) and the data root per group, .../harddisks_gen3/stage1/<group>/chunks/chunk_<j>;
#   - the time limit comes from the sbatch line (--time=<the chunk class>; the line below is only a default), the job name too
#     (--job-name=gen3-s1-<group>, so that the chunk check and the merge find a group's logs and jobs).
# NO GIT AT RUN TIME (node cn-05-08-07 has none and is excluded). Every start-up check comes BEFORE the line "chunk <j> on ...", so a
# chunk that stops in them ran nothing; every STOP line names its node: "STOP: [<node>] ...".
#   1. the clone: ~/harddisks_gen3_k2/hspist3 (the scripts clone, updated in place to the stage 1 package commit by the runsheet);
#   2. the inputs: s1_inputs.sha256 (sha256sum -c): the task list, its SHA-256 file, the chunk table, the build record, the worker, koa_env.sh;
#      and the task list's SHA-256 is the registered one;
#   3. koa_env.sh;
#   4. the binary: its SHA-256 equals the build record's AND the build-time record's 00ALLINONE line, and its version line is the recorded one;
#   5. the chunk number is a row of the chunk table, and its lines lie within the task list.
# Submit from ~/harddisks_gen3_k2/hspist3 (runsheet_stage1_261012.txt), e.g. chunk 1 alone:
#   sbatch --job-name=gen3-s1-B100 --time=0:30:00 --array=1 cluster/gen3_stage1_261012/s1_koa.sbatch
#SBATCH --job-name=gen3-s1
#SBATCH --partition=shared
#SBATCH --account=uh
#SBATCH --time=0:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --exclude=cn-05-08-07
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.out
set -uo pipefail
stop() { echo "STOP: [$(hostname)] $*"; exit 2; }
cd "${SLURM_SUBMIT_DIR:?submit with sbatch from ~/harddisks_gen3_k2/hspist3}"
P=cluster/gen3_stage1_261012
ROOT_NAME=$(basename "$(dirname "$SLURM_SUBMIT_DIR")")
[ "$ROOT_NAME" = harddisks_gen3_k2 ] || stop "stage 1 runs from ~/harddisks_gen3_k2, not $ROOT_NAME"
sha256sum --status -c "$P/s1_inputs.sha256" || stop "the stage 1 inputs are not the delivered ones ($P/s1_inputs.sha256)"
TASKS="$P/tasks_stage1_koa_261012.txt"; CT="$P/chunks_stage1_koa_261012.tsv"
[ "$(sha256sum "$TASKS" | cut -d' ' -f1)" = "$(cat "$P/tasks_stage1_koa_261012.sha256")" ] || stop "the task list is not the registered one"
source cluster/koa_env.sh || stop "cluster/koa_env.sh failed"
REC="$P/build_record_f72d764.txt"; A="$HOME/harddisks_gen3_s1/hspist3"; BIN="$A/00ALLINONE"; BREC="$A/logs/BUILD_RECORD_S1.txt"
want_sha=$(awk '$1 == "sha256" && $3 == "00ALLINONE" {print $2}' "$REC"); want_ver=$(sed -n 's/^version //p' "$REC")
[ -n "$want_sha" ] && [ -n "$want_ver" ] || stop "$REC lacks its 00ALLINONE sha256 line or its version line"
have_sha=$(sha256sum "$BIN" | cut -d' ' -f1)
[ "$have_sha" = "$want_sha" ] || stop "$BIN is not the recorded binary (SHA-256 $have_sha)"
grep -qx "sha256 $want_sha 00ALLINONE" "$BREC" || stop "the build-time record $BREC does not hold the recorded SHA-256 of 00ALLINONE"
HD_BUILD=$("$BIN" --version | head -1)
[ "$HD_BUILD" = "$want_ver" ] || stop "the binary's version line is '$HD_BUILD', recorded '$want_ver'"
j=${SLURM_ARRAY_TASK_ID:?submit as an array}
row=$(awk -F'\t' -v j="$j" 'NR > 1 && $1 == j' "$CT")
[ -n "$row" ] || stop "chunk $j is not a row of the chunk table"
first=$(printf '%s\n' "$row" | cut -f2); last=$(printf '%s\n' "$row" | cut -f3); grp=$(printf '%s\n' "$row" | cut -f5)
[ "$first" -ge 1 ] && [ "$last" -ge "$first" ] && [ "$last" -le "$(wc -l < "$TASKS")" ] || stop "chunk $j's lines $first-$last are not within the task list"
export HD_BIN="$BIN" HD_BUILD
export HD_DATA="/mnt/lustre/koa/scratch/charing/harddisks_gen3/stage1/$grp/chunks/chunk_$j"; mkdir -p "$HD_DATA"
echo "chunk $j on $(hostname) ($(lscpu | awk -F: '/Model name/{gsub(/^ +/,"",$2); print $2; exit}')): group $grp, lines $first-$last; $HD_BUILD"
sed -n "${first},${last}p" "$TASKS" | while read -r line; do
  set -- $line
  bash "$P/s1_worker_koa.sh" "$@"
done
echo "chunk $j done"
