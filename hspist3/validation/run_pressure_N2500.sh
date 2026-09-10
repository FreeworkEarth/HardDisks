#!/usr/bin/env bash
# ##CHRIS: TASK 7 (prompt-1 review §7 item 8; go given 2026-09-10): N = 2500 at
# eta = 0.65, 0.67, 0.69, 3 seeds. Equilibrated disposable calibration (seeds
# 911000001/2, stderr kept), production chunk = min(0.6*min(A,B), 0.75*320/N =
# 0.096), 400 + 30x20 protocol, all-or-nothing health contract, resumable by
# (eta, N, seed). Nothing else may start on the machine while it runs.
set -uo pipefail
D="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$D/.."
OUT="${1:?usage: run_pressure_N2500.sh OUTDIR [JOBS]}"; JOBS="${2:-9}"
N=2500; CAP=0.096; ETAS=(0.65 0.67 0.69); NBLOCKS=30; BDT=20; EQ=400
mkdir -p "$OUT/runs"; CAL="$OUT/chunk_calibration_N2500.csv"
[ -f "$CAL" ] || echo "eta,N,safe_chunk_seedA,safe_chunk_seedB,production_chunk,cap_0.75x320overN" > "$CAL"
calibrate_cell(){
  local eta="$1" a b prod
  a=$(./validation/pressure_validation "$eta" "$N" 911000001 0 0 0 - - calibrate 2>"$OUT/calib_${eta}_${N}_911000001.err" | tail -1)
  b=$(./validation/pressure_validation "$eta" "$N" 911000002 0 0 0 - - calibrate 2>"$OUT/calib_${eta}_${N}_911000002.err" | tail -1)
  [ -z "$a" ] && a=FAILED; [ -z "$b" ] && b=FAILED
  if [ "$a" = FAILED ] || [ "$b" = FAILED ]; then echo "$eta,$N,$a,$b,FAILED,$CAP" >> "$CAL"; return; fi
  prod=$(awk -v x="$a" -v y="$b" -v c="$CAP" 'BEGIN{m=(x<y)?x:y; p=0.6*m; if(p>c)p=c; printf "%.10g", p}')
  echo "$eta,$N,$a,$b,$prod,$CAP" >> "$CAL"
}
echo "N=2500 pressure run -> $OUT   (calibrating ${ETAS[*]} in parallel)"
for eta in "${ETAS[@]}"; do
  awk -F, -v e="$eta" 'NR>1 && $1==e {f=1} END{exit !f}' "$CAL" && continue
  calibrate_cell "$eta" &
done; wait
echo "calibration:"; cat "$CAL"
for eta in "${ETAS[@]}"; do
  chunk=$(awk -F, -v e="$eta" 'NR>1 && $1==e {print $5; exit}' "$CAL")
  [ -z "$chunk" ] || [ "$chunk" = FAILED ] && { echo "  eta=$eta: calibration_failed, not launched"; continue; }
  for k in 0 1 2; do
    seed=$(( 20260907 + 104729*k + 7919*N ))
    tj="$OUT/runs/main_traj_${eta}_${N}_${seed}.csv"; [ -f "$tj" ] && continue
    ( ./validation/pressure_validation "$eta" "$N" "$seed" "$NBLOCKS" "$BDT" "$EQ" "$tj" "$OUT/runs/main_blk_${eta}_${N}_${seed}.csv" "$chunk" \
        >> "$OUT/run.log" 2>&1 || rm -f "$tj" ) &
    while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do sleep 2; done
  done
done; wait
echo "accepted: $(ls "$OUT"/runs/main_traj_*.csv 2>/dev/null | wc -l | tr -d ' ')/9   discards: $(grep -c 'DISCARD\|valid=0' "$OUT/run.log"; true)"
echo "Complete: $OUT"
