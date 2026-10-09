#!/bin/bash
# ##CHRIS 2026-10-09 (261012 sec. 4.7.20; stage G): PREPARED, NOT RUN, NOT REGISTERED. TEST G ON KOA, the array of one-core chunks.
# Array task j runs the task lines (j-1)*25+1 .. j*25 of cluster/gen3_koa_261009/tasks_testG_koa_261009.txt in order, one at a time,
# with testG_worker_koa.sh (the Mac worker's commands and health rule; every A-fixed event log kept on scratch), writing ONLY into
# its own chunk folder $DATA/chunks/chunk_<j> (the same relative layout as the Mac): no two nodes write one file. The merge job
# (merge_testG_koa.sbatch, --dependency=afterok) then assembles the cells, reduces them and checks the inventory.
# Run from ~/harddisks_gen3/hspist3 after build_gen3_koa.sh and the registration (runsheet step 5):
#   sbatch --array=1-120%64 cluster/gen3_koa_261009/testG_koa.sbatch
#SBATCH --job-name=gen3-testG
#SBATCH --partition=shared
#SBATCH --account=uh
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.out
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:?submit with sbatch from ~/harddisks_gen3/hspist3}"
source cluster/koa_env.sh || exit 1
ROOT_NAME=$(basename "$(dirname "$SLURM_SUBMIT_DIR")")
[ "$ROOT_NAME" = harddisks_gen3 ] || { echo "STOP: Test G on KOA runs from ~/harddisks_gen3, not $ROOT_NAME"; exit 2; }
head=$(git rev-parse --short HEAD)
sha256sum --status -c logs/BUILD_KOA_LAST.hash && ./00ALLINONE --version | head -1 | grep -q -- "git $head  target koa" \
  || { echo "STOP: ./00ALLINONE is not the recorded clean koa build of HEAD $head"; exit 2; }
TASKS=cluster/gen3_koa_261009/tasks_testG_koa_261009.txt
EXPECT=$(cat cluster/gen3_koa_261009/tasks_testG_koa_261009.sha256 2>/dev/null)   # written at the registration commit
[ -n "$EXPECT" ] && [ "$(sha256sum "$TASKS" | cut -d' ' -f1)" = "$EXPECT" ] || { echo "STOP: the task list is not the registered one"; exit 2; }
j=${SLURM_ARRAY_TASK_ID:?submit as an array}
export HD_BIN="$PWD/00ALLINONE" HD_BUILD="$(./00ALLINONE --version | head -1)"
export HD_DATA="/mnt/lustre/koa/scratch/charing/harddisks_gen3/testG_koa/chunks/chunk_$j"; mkdir -p "$HD_DATA"
echo "chunk $j on $(hostname) ($(lscpu | awk -F: '/Model name/{gsub(/^ +/,"",$2); print $2; exit}')): lines $(( (j-1)*25+1 ))-$(( j*25 )); $HD_BUILD"
sed -n "$(( (j-1)*25+1 )),$(( j*25 ))p" "$TASKS" | while read -r line; do
  set -- $line
  if [ "$1" = AF ]; then KEEP_EV=1 bash cluster/gen3_koa_261009/testG_worker_koa.sh "${@:1:11}"
  else bash cluster/gen3_koa_261009/testG_worker_koa.sh "${@:1:11}"; fi
done
echo "chunk $j done"
