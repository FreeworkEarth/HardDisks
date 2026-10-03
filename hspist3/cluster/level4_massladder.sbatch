#!/bin/bash
# ##CHRIS 2026-10-01: Level 4 tau_T mass ladder -- the equilibrium fluctuation route.
#
#   sbatch --array=1-80%80 --partition=shared cluster/level4_massladder.sbatch <M_d>
#
# One array task = one seed. M_d is the FIRST ARGUMENT (10, 50 or 200), not a manifest, because the
# ladder is three separate submissions with very different walltimes.
#
# WHY THIS SHAPE (260929/260930 reports): tau_T is measured from the STATIONARY autocorrelation, so
# there is no decaying step to lose against the noise floor. What limits it is the estimator's BIAS,
# which is common-mode across seeds and therefore immune to adding seeds -- only a longer record
# fixes it. Measured separation between the hard-disk (50.35 M) and ideal (61.84 M) predictions:
#
#     record 8200 sigma / L-over-tau 16 / 20 seeds -> 1.4 sigma
#     record 8200      / 16           / 80 seeds  -> 1.7 sigma   <- seeds alone do not do it
#     record 32800     / 65           / 20 seeds  -> 1.9 sigma
#     record 32800     / 65           / 80 seeds  -> 4.0 sigma   <- THIS
#
# so: 65 tau_T, 80 seeds, and the modelled ACF estimator (OU + damped cosine with omega fixed by
# cot K = alpha K), whose bias slope reaches 0.74 at this record length.
#
# TRACE SIZE. The energy-transfer trace writes ONE ROW PER STEP (~368 bytes) and ignores
# --output-dt; undecimated, the M_d = 200 cell is 14.5 GB PER SEED. --trace-every=N decimates in
# steps (exact, seed-independent; N<=1 is the historical behaviour). N is chosen from the physics --
# ~25 samples per divider-mode period, where the Md10 run's 19 already gave a clean fit -- and then
# checked against a 100 MB ceiling, which it clears by more than 25x:
#
#  M_d  alpha   K       period  tau_T   record    steps   undecimated  N    dt     samp/period  size
#   10  0.10  1.4289     95.2     504    32 728    2.0M      0.7 GB   200   3.33       29       3.6 MB
#   50  0.50  1.0769    126.3    2518   163 638    9.8M      3.6 GB   300   5.00       25      12.0 MB
#  200  2.00  0.6533    208.1   10070   654 550   39.3M     14.5 GB   500   8.33       25      28.9 MB
#
# (period = 2 pi L_eff/(c_s K), the Roman divider mode -- Paper 1's equation, which is also what the
# oscillatory half of this autocorrelation is. Two papers, one eigenvalue problem.)
#
#SBATCH --job-name=hd-l4tau
#SBATCH --partition=shared
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.out
#SBATCH --requeue

set -u
MD="${1:?usage: sbatch --array=1-80 cluster/level4_massladder.sbatch <M_d: 10|50|200>}"
TASK="${SLURM_ARRAY_TASK_ID:?this script must run as an array task}"
ROOT="${HD_ROOT:-$SLURM_SUBMIT_DIR}"
BIN="${HD_BIN:-$ROOT/00ALLINONE}"

# steps and decimation per mass. STEPS = 65 * 50.35 * M_d * 60.0294 steps-per-sigma (fixed-dt=0.4).
case "$MD" in
  10)  STEPS=1964913;  EVERY=200 ;;
  50)  STEPS=9824565;  EVERY=300 ;;
  200) STEPS=39298260; EVERY=500 ;;
  *)   echo "M_d must be 10, 50 or 200 (got '$MD')"; exit 2 ;;
esac

SEED=$((9200 + TASK - 1))
OUT="${HD_SCRATCH:-$HOME/koa_scratch}/level4_massladder/Md${MD}"
mkdir -p "$OUT"

# Restart safety -- required, because kill-shared is preemptable and --requeue is set.
# "Complete" means the reduced file exists and is non-empty; the raw trace is deleted after
# reduction, so its absence is not evidence of failure.
if [ -s "$OUT/red_${SEED}.csv" ]; then
  echo "task $TASK: seed $SEED already reduced, skipping"; exit 0
fi

cmd=("$BIN" --mode=edmd --experiment=energy_transfer --headless --quiet
     --edmd-acc=0 --seed-drift-order=drift-first
     "--energy-transfer-summary=$OUT/summary.csv"
     "--energy-transfer-trace=$OUT/raw_${SEED}.csv"
     "--trace-every=$EVERY"
     --particles=100 --particles-boxes=50,50 --particle-radius=0.5
     --l0=39.25 --height=10 --num-walls=1 --wall-positions=39.25
     "--wall-mass-factors=$MD" --eff-output=wall-ke
     --wall-hold-steps=12000 "--steps=$STEPS" --fixed-dt=0.4 --kbt1 "--seed=$SEED")

t0=$SECONDS
"${cmd[@]}" > "$OUT/run_${SEED}.log" 2>&1
rc=$?
dt=$((SECONDS - t0))

# node provenance -- acceptance across nodes is STATISTICAL, never byte-identical, so every run
# must say which CPU produced it. If a cell ever splits into two clusters, look here first.
{ echo; echo "## node provenance";
  echo "- hostname: $(hostname)";
  echo "- cpu: $(awk -F: '/model name/{print $2; exit}' /proc/cpuinfo | sed 's/^ *//')";
  echo "- slurm: job ${SLURM_ARRAY_JOB_ID:-none} task ${TASK} on ${SLURM_JOB_PARTITION:-none}";
  echo "- binary: $(sha1sum "$BIN" | cut -d' ' -f1)";
  echo "- trace-every: $EVERY   steps: $STEPS   wall-clock: ${dt}s";
} >> "$OUT/00_COMMAND.md"

health=$(grep -cE 'EDMD-HEALTH|forced_advance|clamp_repair|overlap_repair|wall_overdue' \
         "$OUT/run_${SEED}.log" 2>/dev/null) || true
health=${health:-0}
if [ "$rc" -ne 0 ] || [ "$health" -ne 0 ]; then
  echo "task $TASK FAILED rc=$rc health=$health -- raw trace KEPT for diagnosis"; exit 1
fi

# Reduce to the four analysed columns and delete the raw trace. Even decimated, keeping 80 raw
# traces per mass is pointless: the analysis reads Time, the two gas energies and the divider x.
python3 - "$OUT/raw_${SEED}.csv" "$OUT/red_${SEED}.csv" <<'PY'
import sys, pandas as pd
cols = ["Time", "KE_gas_left", "KE_gas_right", "W0_x_sigma"]
pd.read_csv(sys.argv[1], low_memory=False, usecols=cols).to_csv(sys.argv[2], index=False)
PY
rc2=$?
if [ "$rc2" -ne 0 ] || [ ! -s "$OUT/red_${SEED}.csv" ]; then
  echo "task $TASK: reduction FAILED, raw trace kept"; exit 1
fi
rm -f "$OUT/raw_${SEED}.csv"
echo "task $TASK done: M_d=$MD seed=$SEED ${dt}s health=0 -> $(stat -c%s "$OUT/red_${SEED}.csv") bytes"
