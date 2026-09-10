#!/usr/bin/env bash
# ##CHRIS: pressure validation, stage 2 -- dense stable-fluid completion.
#
# The eta <= 0.5 grid is already complete (96/96 accepted, zero discards) and is
# preserved separately; this stage extends the result far enough into the dense
# stable fluid to overlap the speed-of-sound validation, plus a small, explicitly
# EXPLORATORY transition diagnostic.
#
# The change that makes it affordable: calibration is amortised per (eta,N)
# instead of repeated per scientific seed. Running a full disposable calibration
# before every seed cost ~280 simulated time units per trajectory and turned the
# dense tail into a multi-day job.
#
#   for each (eta,N):
#       calibrate on TWO disposable seeds        (one configuration can
#       chunk = 0.8 * min(safe_A, safe_B)         underestimate a bursty rate)
#       run every scientific seed from a FRESH state at that fixed chunk
#
# Caching is a SPEED optimisation only. The guarantee is unchanged: any health
# event anywhere in a production trajectory discards that whole trajectory. A
# production seed is never adapted and continued.
set -uo pipefail
D="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$D/.."
OUT="${1:?usage: run_pressure_campaign2.sh OUTDIR [JOBS]}"; JOBS="${2:-12}"
mkdir -p "$OUT/runs"
CAL="$OUT/chunk_calibration.csv"
[ -f "$CAL" ] || echo "eta,N,safe_chunk_seedA,safe_chunk_seedB,production_chunk" > "$CAL"

# Main validation grid: dense stable fluid, overlapping the sound-speed domain.
MAIN_ETAS=(0.60 0.65 0.67 0.69)
# Exploratory only -- NOT part of the EOS acceptance range.
DIAG_ETAS=(0.702 0.710 0.720)
seeds_for(){ case "$1" in 400) echo 5;; 900) echo 4;; *) echo 3;; esac; }
equil_for(){ awk -v e="$1" 'BEGIN{print (e>0.60)? 400 : 200}'; }
bdt_for(){   awk -v e="$1" 'BEGIN{print (e>0.60)? 20 : 12}'; }
NBLOCKS=30

# Calibration is expensive (up to 14 verification intervals) and there are ~18
# (eta,N) cells. Doing it inline inside the production loop serialises the whole
# campaign behind one probe at a time -- observed 12+ min for a single cell with
# every core idle. So calibrate ALL cells first, in parallel, then produce.
calibrate_cell(){
  local eta="$1" N="$2"
  local a b prod
  # calibration stderr is kept next to the calibration table: the table goes into the
  # paper, and why a candidate failed its verification must stay recoverable.
  a=$(./validation/pressure_validation "$eta" "$N" 911000001 0 0 0 - - calibrate 2>"$(dirname "$CAL")/calib_${eta}_${N}_911000001.err" | tail -1)
  b=$(./validation/pressure_validation "$eta" "$N" 911000002 0 0 0 - - calibrate 2>"$(dirname "$CAL")/calib_${eta}_${N}_911000002.err" | tail -1)
  if [ -z "$a" ] || [ "$a" = "FAILED" ] || [ -z "$b" ] || [ "$b" = "FAILED" ]; then
    echo "$eta,$N,${a:-none},${b:-none},CALIBRATION_FAILED" >> "$CAL"
    printf "  [calib] eta=%-6s N=%-5s FAILED (A=%s B=%s) -- cell will not run\n" "$eta" "$N" "${a:-none}" "${b:-none}"
    return
  fi
  # ##CHRIS 2026-09-08 -- 2a. Production chunk = min(0.6*min(A,B), 320/N) for
  # eta >= 0.6. The core's per-call event budget (EDMD_ADVANCE_MAX_EVENTS =
  # 250000) does not scale with N, so the safe chunk scales ~1/N; the 2026-09-07
  # campaign measured 0.8 / 0.4 / 0.2 at N = 400 / 900 / 1600 -- exactly 320/N.
  # Cells calibrated above that line all failed (17 discards). Applied at the
  # runner level: no core change during a campaign.
  prod=$(awk -v x="$a" -v y="$b" -v n="$N" -v e="$eta" \
         'BEGIN{m=(x<y)?x:y; c=0.6*m; cap=320.0/n; if(e>=0.6 && c>cap) c=cap; printf "%.10g", c}')
  echo "$eta,$N,$a,$b,$prod" >> "$CAL"
  printf "  [calib] eta=%-6s N=%-5s A=%-10s B=%-10s -> %s\n" "$eta" "$N" "$a" "$b" "$prod"
}

calibrate_all(){
  for eta in "${MAIN_ETAS[@]}" "${DIAG_ETAS[@]}"; do
    for N in 400 900 1600; do
      awk -F, -v e="$eta" -v n="$N" 'NR>1 && $1==e && $2==n {found=1} END{exit !found}' "$CAL" && continue
      calibrate_cell "$eta" "$N" &
      while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do sleep 1; done
    done
  done
  wait
}

lookup_chunk(){
  # a CALIBRATION_FAILED row must yield nothing, so the cell is skipped rather
  # than the runner receiving a non-numeric chunk and self-calibrating silently
  awk -F, -v e="$1" -v n="$2" 'NR>1 && $1==e && $2==n && $5!="CALIBRATION_FAILED" {print $5; exit}' "$CAL"
}

run_grid(){   # label eta-array-name
  local label="$1"; shift
  for eta in "$@"; do
    for N in 400 900 1600; do
      if [ "$label" = "diag" ] && [ "$N" = "1600" ]; then continue; fi
      local ns; ns=$(seeds_for $N)
      if [ "$label" = "diag" ]; then ns=3; fi
      local chunk; chunk=$(lookup_chunk "$eta" "$N")
      [ -z "$chunk" ] && { echo "  no calibration for eta=$eta N=$N, skipping" >&2; continue; }
      local eq bd; eq=$(equil_for "$eta"); bd=$(bdt_for "$eta")
      for ((k=0;k<ns;k++)); do
        local seed=$(( 20260907 + 104729*k + 7919*N ))
        local tj="$OUT/runs/${label}_traj_${eta}_${N}_${seed}.csv"
        [ -f "$tj" ] && continue          # resume by (eta,N,seed)
        local bk="$OUT/runs/${label}_blk_${eta}_${N}_${seed}.csv"
        ( ./validation/pressure_validation "$eta" "$N" "$seed" "$NBLOCKS" "$bd" "$eq" \
            "$tj" "$bk" "$chunk" >> "$OUT/run.log" 2>&1 || rm -f "$tj" ) &
        while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do sleep 1; done
      done
    done
  done
}

echo "stage 2: calibrating all (eta,N) cells in parallel"
calibrate_all
echo
echo "stage 2: main validation grid (eta = ${MAIN_ETAS[*]})"
run_grid main "${MAIN_ETAS[@]}"
wait
echo "stage 2: transition DIAGNOSTIC (exploratory, eta = ${DIAG_ETAS[*]})"
run_grid diag "${DIAG_ETAS[@]}"
wait

acc=$(ls "$OUT"/runs/*_traj_*.csv 2>/dev/null | wc -l | tr -d ' ')
dis=$(grep -c DISCARD "$OUT/run.log" 2>/dev/null; true)
echo "accepted trajectories: $acc   discards: ${dis:-0}"
echo "Complete: $OUT"
