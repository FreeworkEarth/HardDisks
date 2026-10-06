#!/bin/bash
# ##CHRIS 2026-10-05 (261012 sec. 4.4, engine gates G-E3/G-E4): three cells of the 279282b confinement campaign replayed with the
# engine-divider-resched build -- the campaign's own task lines, i.e. the same masses, positions, seeds and protocol:
#   array 1 = B  e0p10_H_H10_L39.25   (eta 0.10 anchor, N_s = 50, 9 masses x 25 seeds, free divider)
#   array 2 = B  epi8_H_H10_L10       (pi/8 anchor,    N_s = 50, 9 masses x 25 seeds, free divider)
#   array 3 = AF epi8_H_H10_L10       (pi/8 anchor, divider held for the whole record, 5 positions, 115 seeds)
# Data root $SCRATCH/<clone>/hspist3 = $SCRATCH/harddisks_resched/hspist3: the 279282b tree is never touched, and the worker's
# root guard (conf_worker.sh root_guard) records this root's build generation on the first write.
# Cost on the OLD engine (cells_summary.txt and the A-fixed sacct): about 1.4 + 1.0 + 0.6 core-h; 3 x 8 = 24 cores at most.
#   cd ~/harddisks_resched/hspist3 && mkdir -p logs && sbatch --array=1-3 cluster/resched_gate_261005/replay.sbatch
#SBATCH --job-name=resched-replay
#SBATCH --partition=shared
#SBATCH --account=uh
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.out
set -uo pipefail
SCRATCH="/mnt/lustre/koa/scratch/charing"
cd "${SLURM_SUBMIT_DIR:?submit with sbatch from ~/harddisks_resched/hspist3}"
source cluster/koa_env.sh || exit 1
ROOT_NAME=$(basename "$(dirname "$SLURM_SUBMIT_DIR")")
[ "$ROOT_NAME" != harddisks ] || { echo "STOP: the engine gate runs from the branch clone ~/harddisks_resched, not ~/harddisks"; exit 2; }
export HD_BIN="$PWD/00ALLINONE" HD_DATA="$SCRATCH/$ROOT_NAME/hspist3"
# ##CHRIS 2026-10-05 (after the branch review): the engine's contact audit in every trajectory (print only; [EDMD-CONTACT])
export HD_CONTACT_AUDIT=1
gcc --version | head -1; "$HD_BIN" --version | head -2
[ -s logs/BUILD_KOA_LAST.hash ] || { echo "STOP: no logs/BUILD_KOA_LAST.hash -- build with cluster/build_koa.sh first"; exit 1; }
sha256sum --status -c logs/BUILD_KOA_LAST.hash || { echo "STOP: ./00ALLINONE is not the build recorded in logs/BUILD_KOA_LAST.hash"; exit 1; }
export HD_BUILD="$("$HD_BIN" --version | head -1)"
echo "$HD_BUILD" | grep -Eq -- "git [0-9a-f]+  target koa" || { echo "STOP: not a clean koa build: $HD_BUILD"; exit 1; }
case "${SLURM_ARRAY_TASK_ID:-0}" in
  1) MODE=B;  CELL=e0p10_H_H10_L39.25 ;;
  2) MODE=B;  CELL=epi8_H_H10_L10 ;;
  3) MODE=AF; CELL=epi8_H_H10_L10 ;;
  *) echo "STOP: submit with --array=1-3"; exit 2 ;;
esac
TASKS=cluster/confinement_20261013/tasks_${MODE}_${CELL}.txt
echo "replay $MODE $CELL: $(wc -l < "$TASKS") trajectories (the 279282b campaign's task lines), 8 in parallel, into $HD_DATA"
xargs -P 8 -L 1 bash cluster/confinement_20261013/conf_worker.sh < "$TASKS"
if [ "$MODE" = B ]; then
  python3 cluster/confinement_20261013/reduce_B.py "$HD_DATA/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/confinement_B_20261013/$CELL" \
    || echo "replay B $CELL reduce_B FAILED (exit $?)"   # ##CHRIS 2026-10-05: counted below, so 'failures: 0' cannot hide it
fi
echo "replay $MODE $CELL done; failures: $(grep -c FAILED logs/resched-replay_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out || true)"
