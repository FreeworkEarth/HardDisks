#!/bin/bash
# ##CHRIS 2026-10-08 HST (2026-10-09 on the plan author's clock; 261012 sec. 4.4.14): TEST T-PRIME IN THE AMENDED JOB SHAPE.
# Registered and unchanged: the binary 7b08827 target koa (hash and version checks as in testTprime.sbatch), the task list
# cdef566b (seeds, cell epi8_H_H10_L10, method B, line order), HD_CONTACT_AUDIT=1, the worker (cluster/confinement_20261013/
# conf_worker.sh, which records the node of every trajectory in its ##RUN header), the verdict script and its rule.
# Amended (job shape only): 32 one-core array tasks of 50 consecutive task-list lines each, run one after another inside the
# task, so the registered interleaving (minimal, legacy, minimal, ... seed by seed) is kept inside every chunk:
#   chunk i = lines 50(i-1)+1 .. 50i; chunks 1-16 are M = 300, chunks 17-32 are M = 1500.
# Each chunk writes into ITS OWN data root, $SCRATCH/harddisks_resched2/tprime_chunks_261009/cNN, never into a directory another
# chunk writes: the worker appends every trajectory's ##RUN section to the cell's run.log under flock, and flock is not known to
# serialise across nodes on KOA's Lustre (in the registered shape one node wrote each cell). reduce_tprime.sbatch, one job
# with --dependency=afterok:<array id>, merges the 32 roots into the registered layout ($SCRATCH/harddisks_resched2/hspist3/...),
# sequentially in one process, and then runs reduce_B per policy and mass.
# --time: 50 trajectories x 40 s (upper end of the measured KOA times per trajectory) x 2 = 66.7 min -> 1:10:00. Measured
#   (Test T's KOA times, 261012 sec. 4.4.13): an M = 1500 chunk = 25 x (30.6 + 44.8) s = 31.4 min, an M = 300 chunk =
#   25 x (10.9 + 16.0) s = 11.2 min. 32 tasks x 1 core = 32 cores, under the 64-core cap (runsheet sec. 14).
#   cd ~/harddisks_resched2/hspist3 && sbatch --array=1-32%32 cluster/resched_gate_261005/testTprime_chunks.sbatch
# DRY RUN (Mac plumbing test only, never on KOA): TPRIME_DRY_TASKS, TPRIME_DRY_BIN, TPRIME_DRY_ROOT (absolute paths) select a
# task list of NON-registered seeds (refused if any seed or base occurs in the registered list), that binary and a data root
# outside KOA; the KOA checks are skipped, everything else is the same code. Plumbing only: its outputs are not analysed.
#SBATCH --job-name=resched-testTprime-chunks
#SBATCH --partition=shared
#SBATCH --account=uh
#SBATCH --time=1:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.out
set -uo pipefail
SCRATCH="/mnt/lustre/koa/scratch/charing"
CHUNK=50; NCHUNK=32
REG=cluster/resched_gate_261005/tasks_Tprime_epi8_H_H10_L10.txt
REL_TP=experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/resched_testTprime_261007
cd "${SLURM_SUBMIT_DIR:?submit with sbatch from ~/harddisks_resched2/hspist3}"
ROOT_NAME=$(basename "$(dirname "$SLURM_SUBMIT_DIR")")
[ "$ROOT_NAME" = harddisks_resched2 ] || { echo "STOP: Test T-prime runs from the gate-v2 clone (~/harddisks_resched2), not $ROOT_NAME"; exit 2; }
if [ -n "${TPRIME_DRY_TASKS:-}" ]; then
  : "${TPRIME_DRY_BIN:?dry run: set TPRIME_DRY_BIN}"; : "${TPRIME_DRY_ROOT:?dry run: set TPRIME_DRY_ROOT}"
  case "$TPRIME_DRY_TASKS$TPRIME_DRY_BIN$TPRIME_DRY_ROOT" in *"/mnt/lustre/"*) echo "STOP: the dry run never runs on KOA"; exit 2 ;; esac
  case "$TPRIME_DRY_ROOT" in /*) ;; *) echo "STOP: TPRIME_DRY_ROOT must be absolute"; exit 2 ;; esac
  awk 'NR == FNR { s[$5] = 1; b[$10] = 1; next } ($5 in s) || ($10 in b) { bad++ } END { exit bad > 0 }' "$REG" "$TPRIME_DRY_TASKS" \
    || { echo "STOP: the dry task list contains a registered seed or base"; exit 2; }
  echo "DRY RUN -- Mac plumbing test, NOT the registered test; outputs are not analysed"
  TASKS=$TPRIME_DRY_TASKS; BASE=$TPRIME_DRY_ROOT
  export HD_BIN="$TPRIME_DRY_BIN" HD_CONTACT_AUDIT=1
  "$HD_BIN" --version | head -2
else
  source cluster/koa_env.sh || exit 1
  BASE="$SCRATCH/$ROOT_NAME"
  export HD_BIN="$PWD/00ALLINONE" HD_CONTACT_AUDIT=1
  gcc --version | head -1; "$HD_BIN" --version | head -2
  [ -s logs/BUILD_KOA_LAST.hash ] || { echo "STOP: no logs/BUILD_KOA_LAST.hash"; exit 1; }
  sha256sum --status -c logs/BUILD_KOA_LAST.hash || { echo "STOP: ./00ALLINONE is not the build recorded in logs/BUILD_KOA_LAST.hash"; exit 1; }
  TASKS=$REG
fi
export HD_BUILD="$("$HD_BIN" --version | head -1)"
if [ -z "${TPRIME_DRY_TASKS:-}" ]; then
  [ "$HD_BUILD" = "00ALLINONE  git 7b08827  target koa" ] || { echo "STOP: not the Test T binary (git 7b08827, target koa): $HD_BUILD"; exit 1; }
  echo "binary sha256 $(cut -d' ' -f1 logs/BUILD_KOA_LAST.hash) = logs/BUILD_KOA_LAST.hash; $HD_BUILD"
  echo "task list SHA-256 $(sha256sum "$TASKS" | cut -d' ' -f1) (expected cdef566b2a74fd0ec9122eee2ed2fef271fa9ee75d9630de378819565831992f)"
  [ "$(sha256sum "$TASKS" | cut -d' ' -f1)" = cdef566b2a74fd0ec9122eee2ed2fef271fa9ee75d9630de378819565831992f ] || { echo "STOP: task list differs from the registered one"; exit 2; }
fi
[ "$(wc -l < "$TASKS")" -eq $((CHUNK * NCHUNK)) ] || { echo "STOP: the task list has $(wc -l < "$TASKS") lines, not $((CHUNK * NCHUNK))"; exit 2; }
# never mix shapes: a registered cell holding T-prime traces must have been filled by reduce_tprime.sbatch (its marker);
# traces without the marker were written by the registered shape (job 15022392) -- then this shape must not run
for c in "$BASE/hspist3/$REL_TP"/*/epi8_H_H10_L10/m_*; do
  [ -d "$c" ] || continue
  if compgen -G "$c/wall_x_positions_L0_*_run*.csv" >/dev/null && [ ! -e "$c/.merged_from_chunks_261009" ]; then
    echo "STOP: $c holds T-prime traces not merged from chunks (the registered shape ran?) -- do not mix shapes"; exit 2
  fi
done
i=${SLURM_ARRAY_TASK_ID:-0}
[ "$i" -ge 1 ] 2>/dev/null && [ "$i" -le "$NCHUNK" ] || { echo "STOP: submit with --array=1-$NCHUNK%$NCHUNK (task id '$i')"; exit 2; }
a=$((CHUNK * (i - 1) + 1)); b=$((CHUNK * i))
Ms=$(sed -n "${a},${b}p" "$TASKS" | awk '{print $3}' | sort -u)
[ "$(printf '%s\n' "$Ms" | wc -l)" -eq 1 ] || { echo "STOP: chunk $i mixes masses: $Ms"; exit 2; }
export HD_DATA="$BASE/tprime_chunks_261009/c$(printf %02d "$i")"
echo "Test T-prime chunk $i of $NCHUNK: task-list lines $a-$b, M = $Ms, $(sed -n "${a},${b}p" "$TASKS" | wc -l) trajectories" \
     "(minimal and legacy interleaved), one core on $(hostname), data root $HD_DATA"
sed -n "${a},${b}p" "$TASKS" | xargs -P 1 -L 1 bash cluster/confinement_20261013/conf_worker.sh
LOG="logs/resched-testTprime-chunks_${SLURM_ARRAY_JOB_ID:-dry}_${i}.out"
echo "testTprime chunk $i (M=$Ms, lines $a-$b) done; failures: $(grep -c FAILED "$LOG" 2>/dev/null || true)"
