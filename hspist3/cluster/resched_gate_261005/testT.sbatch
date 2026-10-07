#!/bin/bash
# ##CHRIS 2026-10-07 (261012 sec. 4.4.10 E3, gate version 2): TEST T, minimal against legacy on the SAME binary.
# Cell epi8_H_H10_L10, method B, the campaign's protocol; fresh seeds, 400 per policy at M = 50 and 100 for the other masses
# (amended 2026-10-07, second plan-author decision, item 2; cluster/resched_gate_261005/tasks_T_*.txt, written by
# gen_testT_261007.py; seed list SHA-256 0123f894e13741c22d33b9146e3b3896d76175070220572d335e3223e96d45f4). One array task per
# mass: its 2n trajectories (n seeds x minimal and legacy, interleaved line by line) run 8 at a time on one node, so both policies
# share that node; the worker records the node in every ##RUN header; HD_CONTACT_AUDIT=1 in every trajectory.
# Measured cost (validation/resched_testT_design_261007.py): 8.2 core-h in all; the longest task (M = 2000) 10.9 min and the M = 50
# task 9.7 min on 8 cores -> --time 0:30 (2 x the longest, rounded up); %4 -> at most 32 cores.
# Data root $SCRATCH/<clone>/hspist3 of the gate-v2 clone (~/harddisks_resched2): root_guard records its build there.
#   cd ~/harddisks_resched2/hspist3 && sbatch --array=1-9%4 cluster/resched_gate_261005/testT.sbatch
#SBATCH --job-name=resched-testT
#SBATCH --partition=shared
#SBATCH --account=uh
#SBATCH --time=0:30:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.out
set -uo pipefail
SCRATCH="/mnt/lustre/koa/scratch/charing"
cd "${SLURM_SUBMIT_DIR:?submit with sbatch from ~/harddisks_resched2/hspist3}"
source cluster/koa_env.sh || exit 1
ROOT_NAME=$(basename "$(dirname "$SLURM_SUBMIT_DIR")")
case "$ROOT_NAME" in harddisks|harddisks_resched) echo "STOP: Test T runs from the gate-v2 clone (~/harddisks_resched2)"; exit 2 ;; esac
export HD_BIN="$PWD/00ALLINONE" HD_DATA="$SCRATCH/$ROOT_NAME/hspist3" HD_CONTACT_AUDIT=1
gcc --version | head -1; "$HD_BIN" --version | head -2
[ -s logs/BUILD_KOA_LAST.hash ] || { echo "STOP: no logs/BUILD_KOA_LAST.hash -- build with cluster/build_koa.sh first"; exit 1; }
sha256sum --status -c logs/BUILD_KOA_LAST.hash || { echo "STOP: ./00ALLINONE is not the build recorded in logs/BUILD_KOA_LAST.hash"; exit 1; }
export HD_BUILD="$("$HD_BIN" --version | head -1)"
echo "$HD_BUILD" | grep -Eq -- "git [0-9a-f]+  target koa" || { echo "STOP: not a clean koa build: $HD_BUILD"; exit 1; }
TASKS=cluster/resched_gate_261005/tasks_T_epi8_H_H10_L10.txt
echo "task list SHA-256 $(sha256sum "$TASKS" | cut -d' ' -f1) (expected 60a00704e168d95bde3bac28ecd7fd4257428d94989c724d3ca294d5076e6be1)"
MASSES=(50 100 200 300 500 750 1000 1500 2000)
M=${MASSES[$((${SLURM_ARRAY_TASK_ID:-0} - 1))]:-}
[ -n "$M" ] || { echo "STOP: submit with --array=1-9%4"; exit 2; }
echo "Test T, M = $M: $(awk -v M="$M" '$3 == M' "$TASKS" | wc -l) trajectories (minimal and legacy interleaved), 8 in parallel on $(hostname)"
awk -v M="$M" '$3 == M' "$TASKS" | xargs -P 8 -L 1 bash cluster/confinement_20261013/conf_worker.sh
for pol in minimal legacy; do
  python3 cluster/confinement_20261013/reduce_B.py "$HD_DATA/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/resched_testT_261007/$pol/epi8_H_H10_L10" "$M" \
    || echo "testT M=$M $pol reduce_B FAILED (exit $?)"
done
echo "testT M=$M done; failures: $(grep -c FAILED logs/resched-testT_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out || true)"
