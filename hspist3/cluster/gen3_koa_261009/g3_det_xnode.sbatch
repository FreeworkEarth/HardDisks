#!/bin/bash
# ##CHRIS 2026-10-09 (261012 sec. 4.7.20; stage G; REDESIGNED by plan-author decision 13, item 2d, sec. 4.7.29): PREPARED, NOT RUN.
# Cross-node determinism of gen3 (part (b) of the KOA gate) with ONE-NODE jobs: shared allows one node per job, so this is an ARRAY of
# identical one-core tasks. Every task runs the same two commands with the same seeds in its own folder
# $SCRATCH/harddisks_gen3/gen3_xnode_<array job id>/task_<i>/<cell>/ and records its node (hostname, CPU model). The dependent compare
# job (g3_xnode_compare.sbatch, --dependency=afterany:<array job id>) compares the outputs between tasks that ran on DIFFERENT nodes,
# prints which node pairs were compared, and prints "NOT TESTED" if every task ran on one node. The cells: N = 400 at eta 0.70
# (M = 500, M4 seeding, exact box, 2000 sigma-time held, 200 periods released) and the T-prime cell (N = 100, pi/8, M = 300, the
# T-prime protocol). Compared per cell: the divider trace, the psi6(t) file, speed_of_sound_psi6.csv and the run record's event hash.
# Run from ~/harddisks_gen3/hspist3 after build_gen3_koa.sh (runsheet step 4):
#   sbatch --array=1-8 cluster/gen3_koa_261009/g3_det_xnode.sbatch
#SBATCH --job-name=gen3-xnode
#SBATCH --partition=shared
#SBATCH --account=uh
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.out
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:?submit with sbatch from ~/harddisks_gen3/hspist3}"
source cluster/koa_env.sh || exit 1
ROOT_NAME=$(basename "$(dirname "$SLURM_SUBMIT_DIR")")
[ "$ROOT_NAME" = harddisks_gen3 ] || { echo "STOP: run from ~/harddisks_gen3, not $ROOT_NAME"; exit 2; }
[ -n "${SLURM_ARRAY_TASK_ID:-}" ] || { echo "STOP: submit as an array (sbatch --array=1-8 ...)"; exit 2; }
head=$(git rev-parse --short HEAD)
sha256sum --status -c logs/BUILD_KOA_LAST.hash && ./00ALLINONE --version | head -1 | grep -q -- "git $head  target koa" \
  || { echo "STOP: ./00ALLINONE is not the recorded clean koa build of HEAD $head"; exit 2; }
OUT="/mnt/lustre/koa/scratch/charing/harddisks_gen3/gen3_xnode_${SLURM_ARRAY_JOB_ID}/task_${SLURM_ARRAY_TASK_ID}"
[ ! -e "$OUT" ] || { echo "STOP: $OUT exists -- a task is never repeated into the same folder"; exit 2; }
mkdir -p "$OUT"
NODE="$(hostname) | $(lscpu | awk -F: '/Model name/{gsub(/^ +/,"",$2); print $2; exit}')"
echo "$NODE" > "$OUT/node.txt"
echo "task ${SLURM_ARRAY_TASK_ID} of array ${SLURM_ARRAY_JOB_ID} on $NODE | out $OUT"; ./00ALLINONE --version | head -1
runcell() {   # runcell <cell name> <arguments...>: the command in its own folder
  local cell=$1; shift; mkdir -p "$OUT/$cell"
  ( cd "$OUT/$cell" && HD_CONTACT_AUDIT=1 "$SLURM_SUBMIT_DIR/00ALLINONE" "$@" --speed-sound-run-dir="$OUT/$cell" > run.log 2>&1; echo "exit $?" > exit.txt )
  echo "== $cell: $(cat "$OUT/$cell/exit.txt"); $(grep -o "clean=[01]" "$OUT/$cell/run.log" | tr '\n' ' ')$(grep -o "hash=[0-9a-f]*" "$OUT/$cell/run.log")"
}
L0=$(python3 -c "import math; print('%.17g' % (400 * math.pi / (8 * 20 * 0.70)))")
runcell N400_eta070_M500 --mode=edmd --experiment=speed_of_sound --headless --kbt1 --seed-drift-order=drift-first --edmd-acc=0 \
  --particles=400 --particles-boxes=200,200 --height=20 --particle-radius=0.5 --wall-thickness=0.05 --wall-thickness-vis=0.05 \
  --lengths=$L0 --wall-masses=500 --repeats=1 --seed=20261219 --wall-hold-steps=120000 --fixed-dt=0.4 --target-oscillations=200 \
  --oscillation-safety=1.0 --oscillation-min-steps=10000 --oscillation-max-steps=400000000 --speed-sound-log-stride=60 \
  --speed-sound-exact-seed=4242 --engine=gen3 --gen3-seeding=lattice --gen3-exact-box
runcell Tprime_cell_M300 --mode=edmd --experiment=speed_of_sound --headless --kbt1 --seed-drift-order=drift-first --edmd-acc=0 \
  --particles=100 --particles-boxes=50,50 --height=10.000000 --particle-radius=0.5 --wall-thickness=0.05 --wall-thickness-vis=0.05 \
  --lengths=10.000000 --wall-masses=300 --repeats=1 --seed=20261219 --wall-hold-steps=2000 --fixed-dt=0.4 --target-oscillations=200 \
  --oscillation-safety=1.0 --oscillation-min-steps=10000 --oscillation-max-steps=400000000 --speed-sound-log-stride=51 \
  --speed-sound-exact-seed=4243 --engine=gen3
echo "task ${SLURM_ARRAY_TASK_ID} done"
