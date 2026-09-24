#!/usr/bin/env bash
# ##CHRIS 2026-10-04: the Level 4 ladder RERUN at 65 tau_TRUE. Same experiment, same analysis rules,
# only the record length changed -- which is what the pre-registration always intended.
#
# WHY. The 261003 ladder sized its records at 65 tau_PRED = 65 x 50.35 M. tau turned out up to 3.7x
# that, so L/tau ACTUAL collapsed 64 -> 17 across the ladder, the calibration slope fell 0.683 ->
# 0.476 (below its own 0.6 adoption threshold at M = 100 and 200), the three estimators spread to
# 1.86x at M = 200, and b rose monotonically as heavier masses were added -- the estimator
# signature, not an exponent. Records here are 65 x max(modelled, block, S(0)) so that even the
# PESSIMISTIC tau gets 65 lengths; if the linear law is right these are simply longer than needed.
#
#   M   tau_max  record     steps        period  N    samp/period  rows     MB/seed  core-h
#   50    4155     270 075   16 212 446  125.9   350     21.6       46 321    2.4      1.8
#  100   12935     840 775   50 471 238  157.6   450     21.0      112 158    5.9      5.6
#  200   53404   3 471 260  208 377 735  207.6   600     20.8      347 296   18.1     23.2
#
# 30.6 core-hours, ~3.4 h wall clock at 9 concurrent. 2.1 GB reduced.
#
# TWO CONSTRAINTS COLLIDED and the tie was broken deliberately: ">= 20 samples per mode period"
# and "<= 5 MB per seed" are incompatible at M >= 100 (M = 200 would need N >= 2176 for the size
# cap and N <= 622 for the resolution floor). The RESOLUTION floor is kept, because resolving the
# oscillation is what the two-component fit depends on, and 2.1 GB against 27 GB free is not a
# constraint worth honouring over it. M = 200 is 18.1 MB per seed.
#
# eta_phys = 0.101342, from the RECORDED wall_thickness_sigma = 1.0 in every Level 4 summary.csv
# -- not inferred. c_s = 1.74930, tau_T = 50.349 M (hard disk), 61.836 M (ideal).
#
# M = 10 and M = 20 are NOT repeated: L/tau 64 and 61, calibration slope 0.683 and 0.630, both in
# spec by the pre-registered rules.
set -uo pipefail
cd /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3
P=experiments_energy_transfer/level4_ladder_rerun_20261004
MAXJOBS=${MAXJOBS:-9}

run_one() {
  local M=$1 STEPS=$2 EVERY=$3 sd=$4
  local d=$P/Md${M}
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
  if [ "$rc" -ne 0 ] || [ "$health" -ne 0 ]; then
    echo "M=$M seed=$sd FAILED rc=$rc health=$health -- raw trace KEPT"; return 1
  fi
  python3 - "$d/raw_${sd}.csv" "$d/red_${sd}.csv" <<'PY'
import sys, pandas as pd
cols = ["Time","KE_gas_left","KE_gas_right","W0_x_sigma"]
pd.read_csv(sys.argv[1], low_memory=False, usecols=cols).to_csv(sys.argv[2], index=False)
PY
  if [ -s "$d/red_${sd}.csv" ]; then rm -f "$d/raw_${sd}.csv"; else
    echo "M=$M seed=$sd reduction FAILED -- raw trace KEPT"; return 1; fi
  return 0
}

for M in 200 100 50; do mkdir -p $P/Md${M}; done
echo "ladder rerun start $(date +%H:%M:%S), max $MAXJOBS concurrent"
for spec in "200 208377735 600" "100 50471238 450" "50 16212446 350"; do
  set -- $spec; M=$1; STEPS=$2; EVERY=$3
  for i in $(seq 0 79); do
    while [ "$(jobs -rp | wc -l | tr -d ' ')" -ge "$MAXJOBS" ]; do sleep 5; done
    run_one "$M" "$STEPS" "$EVERY" $((9200+i)) &
  done
done
wait
echo "ladder rerun done $(date +%H:%M:%S)"
for M in 50 100 200; do
  d=$P/Md${M}
  echo "  M=$M: $(find $d -name 'red_*.csv' | wc -l | tr -d ' ')/80 reduced; \
health $(cat $d/run_*.log 2>/dev/null | grep -cE 'EDMD-HEALTH|forced_advance|clamp_repair|overlap_repair|wall_overdue'); \
aborts $(cat $d/run_*.log 2>/dev/null | grep -c ABORTING); \
raw left $(find $d -name 'raw_*.csv' 2>/dev/null | wc -l | tr -d ' ')"
done
du -sh $P
