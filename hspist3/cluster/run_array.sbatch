#!/bin/bash
# ##CHRIS 2026-09-16: one SLURM array task = one trajectory.
#
# Why an array and not one long job: every run is independent and 0.5-4 h, the scheduler can start
# them as slots appear, a failure costs one trajectory rather than the campaign, and a re-submitted
# array skips whatever already finished. Fill PARTITION and TIME from the KOA facts before first use.
#
# usage: sbatch --array=1-<N>%<CONCURRENCY> cluster/run_array.sbatch <manifest.tsv>
#
#SBATCH --job-name=hd-sos
#SBATCH --partition=PARTITION_TBD
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.out
#SBATCH --requeue

set -u
MANIFEST="${1:?usage: sbatch ... run_array.sbatch manifest.tsv}"
ROOT="${HD_ROOT:-$SLURM_SUBMIT_DIR}"
BIN="${HD_BIN:-$ROOT/00ALLINONE}"
TASK="${SLURM_ARRAY_TASK_ID:?this script must run as an array task}"

# line TASK of the manifest, header skipped
read -r idx campaign outdir particles boxes height radius wall lengths wall_mass \
        target safety stride max_steps edmd_acc seed_base exact_seed \
  < <(awk -v n="$TASK" 'NR == n + 1' "$MANIFEST")
[ -n "${exact_seed:-}" ] || { echo "manifest line $TASK missing or short"; exit 2; }

outdir="$(eval echo "$outdir")"          # expands $SCRATCH in the manifest path
mkdir -p "$outdir"

# Restart safety: a completed trajectory is never redone. "Completed" means the trace exists and its
# last sample reaches the planned duration recorded in its own header.
trace="$(ls "$outdir"/wall_x_positions_*run0.csv 2>/dev/null | head -1)"
if [ -n "$trace" ] && awk -F, 'NR==2{p=$13} END{exit !(NR>2 && $1 >= 0.999*p)}' "$trace"; then
  echo "task $TASK: already complete, skipping ($trace)"
  exit 0
fi
rm -f "$outdir"/wall_x_positions_*run0.csv

cmd=("$BIN" --mode=edmd --experiment=speed_of_sound --headless --kbt1
     --seed-drift-order=drift-first "--edmd-acc=$edmd_acc"
     "--particles=$particles" "--particles-boxes=$boxes"
     "--height=$height" "--particle-radius=$radius"
     "--wall-thickness=$wall" "--wall-thickness-vis=$wall"
     "--lengths=$lengths" "--wall-masses=$wall_mass" --repeats=1
     "--seed=$seed_base" --wall-hold-steps=2000 --fixed-dt=0.4
     "--target-oscillations=$target" "--oscillation-safety=$safety"
     --oscillation-min-steps=10000 "--oscillation-max-steps=$max_steps"
     "--speed-sound-log-stride=$stride" "--speed-sound-run-dir=$outdir"
     "--speed-sound-exact-seed=$exact_seed")

printf 'HD_KE_TRACE=1 %s\n' "${cmd[*]}" > "$outdir/command.txt"
echo "task $TASK: $campaign eta-cell $(basename "$(dirname "$(dirname "$outdir")")") M=$wall_mass seed=$exact_seed"
t0=$SECONDS
HD_KE_TRACE=1 "${cmd[@]}" > "$outdir/stdout.log" 2>&1
rc=$?
dt=$((SECONDS - t0))

trace="$(ls "$outdir"/wall_x_positions_*run0.csv 2>/dev/null | head -1)"
ok=0
if [ $rc -eq 0 ] && [ -n "$trace" ] && \
   awk -F, 'NR==2{p=$13} END{exit !(NR>2 && $1 >= 0.999*p)}' "$trace"; then ok=1; fi
# grep -c prints 0 and exits 1 when there is no match, so `|| echo 0` would append a SECOND zero and
# split the provenance line in two. Take the count as-is and default only when the file is missing.
health=$(grep -c "EDMD-HEALTH" "$outdir/stdout.log" 2>/dev/null) || true
health=${health:-0}

# one provenance line per run; the analysis side reads these, never the SLURM logs
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
  "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$TASK" "$exact_seed" "$rc" "$ok" "$health" "$dt" \
  >> "$outdir/../../../../run_records_${SLURM_ARRAY_JOB_ID:-local}.tsv"

if [ "$ok" -ne 1 ]; then echo "task $TASK FAILED rc=$rc"; exit 1; fi
echo "task $TASK done in ${dt}s, health lines $health"
