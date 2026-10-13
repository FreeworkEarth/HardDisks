#!/bin/bash
# ##CHRIS 2026-10-12 (plan-author decision 39 item 5 c; gate G6 (c), cross-node determinism with the campaign binary): an ARRAY of 8
# identical one-core tasks (shared allows one node per job), each running the same two commands with the same seeds in its own folder
# .../harddisks_gen3/s1_xnode/task_<i>/<cell>/ and recording its node. It is cluster/gen3_koa_261009/g3_det_xnode.sbatch (the a643bc7
# cross-node job, decision 13 item 2d) with THESE CHANGES: it runs from the scripts clone ~/harddisks_gen3_k2 WITHOUT GIT (the a643bc7 job
# read `git rev-parse` in ~/harddisks_gen3); the binary is the campaign build ~/harddisks_gen3_s1/hspist3/00ALLINONE, checked against its
# build-time record; the root is fixed (s1_xnode, no job id: a task never writes into an existing folder); the inputs manifest; the
# --exclude line; and the job name gen3-s1-xnode, shared with the compare job, which waits for the array by --dependency=singleton.
# The two cells and their seeds are the a643bc7 job's, unchanged (N 400 at eta 0.70, M 500; the T-prime cell N 100, pi/8, M 300), so its
# event hashes 0586b1e77639502b and 72111cdaa7d26414 (sec. 4.7.46) are printed next to these (information).
#   cd ~/harddisks_gen3_k2/hspist3 && sbatch --array=1-8 cluster/gen3_s1_gate_261012/s1_xnode.sbatch
#SBATCH --job-name=gen3-s1-xnode
#SBATCH --partition=shared
#SBATCH --account=uh
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --exclude=cn-05-08-07
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.out
set -uo pipefail
stop() { echo "STOP: [$(hostname)] $*"; exit 2; }
cd "${SLURM_SUBMIT_DIR:?submit with sbatch from ~/harddisks_gen3_k2/hspist3}"
P=cluster/gen3_s1_gate_261012
ROOT_NAME=$(basename "$(dirname "$SLURM_SUBMIT_DIR")")
[ "$ROOT_NAME" = harddisks_gen3_k2 ] || stop "run from ~/harddisks_gen3_k2, not $ROOT_NAME"
[ -n "${SLURM_ARRAY_TASK_ID:-}" ] || stop "submit as an array (sbatch --array=1-8 ...)"
sha256sum --status -c "$P/s1_gate_inputs.sha256" || stop "the gate inputs are not the delivered ones"
source cluster/koa_env.sh || stop "cluster/koa_env.sh failed"
S1="$HOME/harddisks_gen3_s1/hspist3"; REC="$S1/logs/BUILD_RECORD_S1.txt"; BIN="$S1/00ALLINONE"
want=$(awk '$1 == "sha256" && $3 == "00ALLINONE" {print $2}' "$REC" 2>/dev/null)
[ -n "$want" ] && [ "$(sha256sum "$BIN" | cut -d' ' -f1)" = "$want" ] || stop "$BIN is not the recorded campaign build"
OUT="/mnt/lustre/koa/scratch/charing/harddisks_gen3/s1_xnode/task_${SLURM_ARRAY_TASK_ID}"
[ ! -e "$OUT" ] || stop "$OUT exists -- a task is never repeated into the same folder"
mkdir -p "$OUT"
NODE="$(hostname) | $(lscpu | awk -F: '/Model name/{gsub(/^ +/,"",$2); print $2; exit}')"
echo "$NODE" > "$OUT/node.txt"
echo "task ${SLURM_ARRAY_TASK_ID} on $NODE | out $OUT"; "$BIN" --version | head -1
runcell() {   # runcell <cell name> <arguments...>: the command in its own folder
  local cell=$1; shift; mkdir -p "$OUT/$cell"
  ( cd "$OUT/$cell" && HD_CONTACT_AUDIT=1 "$BIN" "$@" --speed-sound-run-dir="$OUT/$cell" > run.log 2>&1; echo "exit $?" > exit.txt )
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
echo "(the a643bc7 cross-node run, sec. 4.7.46: N400_eta070_M500 hash=0586b1e77639502b, Tprime_cell_M300 hash=72111cdaa7d26414)"
echo "task ${SLURM_ARRAY_TASK_ID} done"
