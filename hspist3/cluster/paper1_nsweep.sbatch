#!/bin/bash
# ##CHRIS 2026-10-01: Paper 1 c_s(eta) at fixed N -- one array task = one (mass, seed) trajectory.
#
#   sbatch --array=1-50%50 --partition=shared cluster/paper1_nsweep.sbatch <eta> <N>
#
# e.g.  for e in 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.60 0.65; do
#           sbatch --array=1-50%50 cluster/paper1_nsweep.sbatch $e 900
#       done
#
# 50 tasks = 5 masses x 10 seeds. The mass ladder is the measurement: c_s is the through-origin
# slope of nu_bar_M against x_M, so every mass in a cell must complete before that cell means
# anything -- which is why the array is per-cell and not per-campaign.
#
# RELATION TO run_array.sbatch: that one is manifest-driven (make_manifest.py writes one line per
# trajectory) and remains the tool for the big multi-cell campaigns -- a2_long200, n900_sweep,
# mirror. This template is the single-cell form: no manifest, submit and read the numbers. Use
# run_array.sbatch when you want one job ID to cover a whole campaign; use this when adding one
# density, or when the manifest machinery is one moving part too many.
#
# GEOMETRY (famB convention, identical to make_manifest.a2_geometry):
#     L0 = (25 pi / 20) / eta * sqrt(N/100),   H = 10 * sqrt(N/100)
# so the aspect ratio is fixed as N grows and eta = N pi r^2 / (2 L0 H) with r = 0.5 exactly.
# eta is the NOMINAL packing fraction on the physical area -- the convention KR is defined with and
# the one Paper 1's eta_of() uses. L_eff enters the acoustic length only, never the density.
#
# NOTE eta <= 0.69: above that KR is outside its fitted range and the solid branch is a lower bound,
# not c_s. Do not add cells above 0.69 to this campaign without deciding what they mean.
#
#SBATCH --job-name=hd-sos-n
#SBATCH --partition=shared
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.out
#SBATCH --requeue

set -u
ETA="${1:?usage: sbatch --array=1-50 cluster/paper1_nsweep.sbatch <eta> <N>}"
NPART="${2:?usage: sbatch --array=1-50 cluster/paper1_nsweep.sbatch <eta> <N>}"
TASK="${SLURM_ARRAY_TASK_ID:?this script must run as an array task}"
ROOT="${HD_ROOT:-$SLURM_SUBMIT_DIR}"
BIN="${HD_BIN:-$ROOT/00ALLINONE}"

MASSES=(50 200 500 1000 2000)
NSEED=10
mi=$(( (TASK - 1) / NSEED ))
si=$(( (TASK - 1) % NSEED ))
[ "$mi" -lt "${#MASSES[@]}" ] || { echo "task $TASK out of range for ${#MASSES[@]} masses x $NSEED seeds"; exit 2; }
M=${MASSES[$mi]}

read -r L0 H < <(awk -v e="$ETA" -v n="$NPART" 'BEGIN{f=sqrt(n/100.0); printf "%.6f %.6f\n", 3.926990816987241/e*f, 10*f}')
SEED=$(( 24030000 + si ))
TAG="eta$(awk -v e="$ETA" 'BEGIN{printf "%.4f", e}')_N${NPART}"
OUT="${HD_SCRATCH:-$HOME/koa_scratch}/paper1_nsweep/${TAG}/m_${M}/seed_${SEED}"
mkdir -p "$OUT"

# Restart safety: kill-shared is preemptable and --requeue is set, so a task may start more than
# once. A trajectory counts as complete only when its trace REACHES ITS PLANNED DURATION -- the
# planned value is in the trace's own header (field 13 of row 2), so a truncated file is redone
# rather than silently accepted as a short record.
trace="$(ls "$OUT"/wall_x_positions_*run0.csv 2>/dev/null | head -1)"
if [ -n "$trace" ] && awk -F, 'NR==2{p=$13} END{exit !(NR>2 && $1 >= 0.999*p)}' "$trace"; then
  echo "task $TASK: already complete, skipping ($trace)"; exit 0
fi
rm -f "$OUT"/wall_x_positions_*run0.csv

cmd=("$BIN" --mode=edmd --experiment=speed_of_sound --headless --kbt1
     --seed-drift-order=drift-first --edmd-acc=0
     "--particles=$NPART" "--particles-boxes=$((NPART/2)),$((NPART/2))"
     "--height=$H" --particle-radius=0.5
     --wall-thickness=0.05 --wall-thickness-vis=0.05
     "--lengths=$L0" "--wall-masses=$M" --repeats=1
     "--seed=$SEED" --wall-hold-steps=2000 --fixed-dt=0.4
     --target-oscillations=200 --oscillation-safety=1.5
     --oscillation-min-steps=10000 --oscillation-max-steps=40000000
     --speed-sound-log-stride=4 "--speed-sound-run-dir=$OUT")

t0=$SECONDS
"${cmd[@]}" > "$OUT/stdout.log" 2>&1
rc=$?
dt=$((SECONDS - t0))

{ echo; echo "## node provenance";
  echo "- hostname: $(hostname)";
  echo "- cpu: $(awk -F: '/model name/{print $2; exit}' /proc/cpuinfo | sed 's/^ *//')";
  echo "- slurm: job ${SLURM_ARRAY_JOB_ID:-none} task ${TASK} on ${SLURM_JOB_PARTITION:-none}";
  echo "- binary: $(sha1sum "$BIN" | cut -d' ' -f1)";
  echo "- cell: eta=$ETA N=$NPART L0=$L0 H=$H M=$M seed=$SEED   wall-clock: ${dt}s";
} >> "$OUT/00_COMMAND.md"

trace="$(ls "$OUT"/wall_x_positions_*run0.csv 2>/dev/null | head -1)"
ok=0
if [ $rc -eq 0 ] && [ -n "$trace" ] && \
   awk -F, 'NR==2{p=$13} END{exit !(NR>2 && $1 >= 0.999*p)}' "$trace"; then ok=1; fi
health=$(grep -cE 'EDMD-HEALTH|forced_advance|clamp_repair|overlap_repair|wall_overdue' "$OUT/stdout.log" 2>/dev/null) || true
health=${health:-0}

if [ "$ok" -ne 1 ] || [ "$health" -ne 0 ]; then
  echo "task $TASK FAILED rc=$rc ok=$ok health=$health"; exit 1
fi
echo "task $TASK done: $TAG M=$M seed=$SEED ${dt}s health=0"
