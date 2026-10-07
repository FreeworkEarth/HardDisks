#!/bin/bash
# ##CHRIS 2026-10-07 (261012 sec. 4.4.13, gate v3): TEST T-PRIME, minimal against legacy on the binary that ran Test T (7b08827).
# Cell epi8_H_H10_L10, method B, the campaign's protocol; masses M = 300 and M = 1500 only; 400 fresh seeds per policy per mass
# (cluster/resched_gate_261005/tasks_Tprime_*.txt, written by gen_testTprime_261007.py; seed list SHA-256
# 01143aac529b7e9d963f40ef489b33b26ee6d1023256059038c6f1b970090381). One array task per mass: its 800 trajectories (400 seeds x
# minimal and legacy, interleaved line by line) run 8 at a time on one node, so both policies share that node; the worker records
# the node in every ##RUN header; HD_CONTACT_AUDIT=1 in every trajectory. Same partition, cores and interleaving as Test T.
# Measured cost (validation/resched_testTprime_261007.py --design, from Test T's KOA times): 11.4 core-h; the longest task
# (M = 1500) 62.9 min on 8 cores -> --time 2:10:00 (2 x the longest, rounded up to 10 min); %2 -> 16 cores.
# The binary is NOT rebuilt: ./00ALLINONE must be the recorded build (logs/BUILD_KOA_LAST.hash) AND print
# "00ALLINONE  git 7b08827  target koa" -- the binary of Test T -- whatever commit the clone's working tree is at.
#   cd ~/harddisks_resched2/hspist3 && sbatch --array=1-2%2 cluster/resched_gate_261005/testTprime.sbatch
#SBATCH --job-name=resched-testTprime
#SBATCH --partition=shared
#SBATCH --account=uh
#SBATCH --time=2:10:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.out
set -uo pipefail
SCRATCH="/mnt/lustre/koa/scratch/charing"
cd "${SLURM_SUBMIT_DIR:?submit with sbatch from ~/harddisks_resched2/hspist3}"
source cluster/koa_env.sh || exit 1
ROOT_NAME=$(basename "$(dirname "$SLURM_SUBMIT_DIR")")
[ "$ROOT_NAME" = harddisks_resched2 ] || { echo "STOP: Test T-prime runs from the gate-v2 clone (~/harddisks_resched2), not $ROOT_NAME"; exit 2; }
export HD_BIN="$PWD/00ALLINONE" HD_DATA="$SCRATCH/$ROOT_NAME/hspist3" HD_CONTACT_AUDIT=1
gcc --version | head -1; "$HD_BIN" --version | head -2
[ -s logs/BUILD_KOA_LAST.hash ] || { echo "STOP: no logs/BUILD_KOA_LAST.hash"; exit 1; }
sha256sum --status -c logs/BUILD_KOA_LAST.hash || { echo "STOP: ./00ALLINONE is not the build recorded in logs/BUILD_KOA_LAST.hash"; exit 1; }
export HD_BUILD="$("$HD_BIN" --version | head -1)"
[ "$HD_BUILD" = "00ALLINONE  git 7b08827  target koa" ] || { echo "STOP: not the Test T binary (git 7b08827, target koa): $HD_BUILD"; exit 1; }
echo "binary sha256 $(cut -d' ' -f1 logs/BUILD_KOA_LAST.hash) = logs/BUILD_KOA_LAST.hash; $HD_BUILD"
TASKS=cluster/resched_gate_261005/tasks_Tprime_epi8_H_H10_L10.txt
echo "task list SHA-256 $(sha256sum "$TASKS" | cut -d' ' -f1) (expected cdef566b2a74fd0ec9122eee2ed2fef271fa9ee75d9630de378819565831992f)"
[ "$(sha256sum "$TASKS" | cut -d' ' -f1)" = cdef566b2a74fd0ec9122eee2ed2fef271fa9ee75d9630de378819565831992f ] || { echo "STOP: task list differs from the registered one"; exit 2; }
MASSES=(300 1500)
M=${MASSES[$((${SLURM_ARRAY_TASK_ID:-0} - 1))]:-}
[ -n "$M" ] || { echo "STOP: submit with --array=1-2%2"; exit 2; }
echo "Test T-prime, M = $M: $(awk -v M="$M" '$3 == M' "$TASKS" | wc -l) trajectories (minimal and legacy interleaved), 8 in parallel on $(hostname)"
awk -v M="$M" '$3 == M' "$TASKS" | xargs -P 8 -L 1 bash cluster/confinement_20261013/conf_worker.sh
for pol in minimal legacy; do
  python3 cluster/confinement_20261013/reduce_B.py "$HD_DATA/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/resched_testTprime_261007/$pol/epi8_H_H10_L10" "$M" \
    || echo "testTprime M=$M $pol reduce_B FAILED (exit $?)"
done
echo "testTprime M=$M done; failures: $(grep -c FAILED logs/resched-testTprime_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out || true)"
