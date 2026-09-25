#!/usr/bin/env bash
# ##CHRIS 2026-10-09: item 6 -- fill the ladder at the R values the R-collapse had to INTERPOLATE.
# Pre-registered in 261006 section 7.4. Runs in the background; gates nothing.
#
# WHY. Pass 2 agreed with the ladder at R = 1, the one MEASURED ladder node, and sat 15-22 % above
# it at R = 4 and R = 2, where f(R) had to be interpolated between the ladder's nodes (5, 2.5, 1,
# 0.5, 0.25). Refitting the ladder with a smooth convex quadratic made A marginally WORSE, so the
# residual is not an interpolation artefact -- but it is either a real box-size dependence at fixed
# R, or structure in the ladder's own f(R) between its nodes. Measuring the ladder AT R = 4 and 2
# separates them, and then no interpolation enters the comparison at all.
#
# OUTCOME I  : f = 1.12 / 1.40 within 2 sigma -> the collapse holds, the ladder's M = 10 and 20
#              anchors are biased, and "GP reproduced at R = 5 and 2.5" must be softened.
# OUTCOME II : f = 0.98 / 1.18 within 2 sigma -> f has a SECOND variable; N and tau_T/T_mode are
#              the candidates to test next.
# Neither -> report.
#
# Ladder box: N_s = 50, L_c = 38.75, wall 39.25, box 78.5, t = 1.0, eta = 0.10134170.
# tau_GP = 50.349 M. Fractional wall_mass_factors=12.5 verified ACCEPTED by the parser.
#   M     R   alpha    K     period  tau_GP  record(65x1.5)  steps      every  s/per
#   12.5  4.0 0.125  1.3978   97.0    629.4   61 400         3 690 000  240    24.3
#   25    2.0 0.250  1.2646  107.2   1258.7  122 800         7 370 000  300    21.4
# 80 seeds each = 883 M steps, ~3.9 core-hours. Seeds 9700-9779.
set -uo pipefail
cd /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3
P=experiments_energy_transfer/level4_ladderRfill_20261009
MAXJOBS=${MAXJOBS:-9}
run_one() {
  local M=$1 STEPS=$2 EVERY=$3 TAG=$4 sd=$5
  local d=$P/$TAG
  [ -s "$d/red_${sd}.csv" ] && return 0
  nice -n 5 ./00ALLINONE --mode=edmd --experiment=energy_transfer --headless --quiet \
      --edmd-acc=0 --seed-drift-order=drift-first \
      --energy-transfer-summary="$d/summary_${sd}.csv" \
      --energy-transfer-trace="$d/raw_${sd}.csv" --trace-every=$EVERY \
      --particles=100 --particles-boxes=50,50 --particle-radius=0.5 \
      --l0=39.25 --height=10 --num-walls=1 --wall-positions=39.25 \
      --wall-mass-factors=$M --eff-output=wall-ke \
      --wall-hold-steps=12000 --steps=$STEPS --fixed-dt=0.4 --kbt1 --seed=$sd \
      > "$d/run_${sd}.log" 2>&1
  local rc=$? health
  health=$(grep -cE 'EDMD-HEALTH|forced_advance|clamp_repair|overlap_repair|wall_overdue' "$d/run_${sd}.log" 2>/dev/null) || true
  health=${health:-0}
  if [ "$rc" -ne 0 ] || [ "$health" -ne 0 ]; then echo "$TAG seed=$sd FAILED rc=$rc health=$health"; return 1; fi
  python3 - "$d/raw_${sd}.csv" "$d/red_${sd}.csv" <<'PY'
import sys, pandas as pd
want = ["Time","KE_gas_left","KE_gas_right","W0_x_sigma","W0_v","SegCounts","SegEtas"]
e = pd.read_csv(sys.argv[1], low_memory=False)
missing = [c for c in want if c not in e.columns]
if missing: sys.exit("REDUCTION ABORTED: columns absent: %s" % ",".join(missing))
e[want].to_csv(sys.argv[2], index=False)
PY
  if [ -s "$d/red_${sd}.csv" ]; then rm -f "$d/raw_${sd}.csv"; else echo "$TAG seed=$sd reduction FAILED"; return 1; fi
  return 0
}
for t in Md25 Md125; do mkdir -p $P/$t; done
echo "ladderRfill start $(date +%H:%M:%S), max $MAXJOBS concurrent"
for spec in "25 7370000 300 Md25" "12.5 3690000 240 Md125"; do
  set -- $spec; M=$1; ST=$2; EV=$3; TAG=$4
  for i in $(seq 0 79); do
    while [ "$(jobs -rp | wc -l | tr -d ' ')" -ge "$MAXJOBS" ]; do sleep 5; done
    run_one "$M" "$ST" "$EV" "$TAG" $((9700+i)) &
  done
done
wait
echo "ladderRfill done $(date +%H:%M:%S)"
for t in Md125 Md25; do d=$P/$t; echo "  $t: $(find $d -name 'red_*.csv'|wc -l|tr -d ' ')/80; health $(cat $d/run_*.log 2>/dev/null | grep -cE 'EDMD-HEALTH|forced_advance|clamp_repair|overlap_repair|wall_overdue'); aborts $(cat $d/run_*.log 2>/dev/null | grep -c ABORTING)"; done
du -sh $P; df -h . | tail -1
