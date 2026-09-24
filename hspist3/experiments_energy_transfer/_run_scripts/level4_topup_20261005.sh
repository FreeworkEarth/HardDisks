#!/usr/bin/env bash
# ##CHRIS 2026-10-05: top-up of the two LIGHT masses, so every mass on the ladder clears L/tau >= 65.
#
# WHY. After the 65 tau_true rerun, M = 50/100/200 sit at L/tau 70/68/87 with calibration slopes
# 0.746/0.741/0.734, while M = 10 and M = 20 -- which were NOT rerun -- are now the WORST sampled:
# L/tau 64 and 58, slopes 0.708 and 0.667. A lower slope means more bias to correct, so if the
# calibration under-corrects at the light end it pulls tau_light down and MANUFACTURES the very
# curvature under test (local b rising 1.18 -> 1.37 -> 1.68 -> 1.70). This run removes that
# possibility. Target is L/tau = 70, not the bare 65, so a slightly larger tau still clears it.
#
#   M   tau now  record now       steps      was        factor  N    rows     MB    core-h
#   10    482     30 735 (L/t 64)  2 144 000  1 964 913  1.09x  200   10 720   0.6    0.24
#   20   1091     63 456 (L/t 58)  4 703 000  3 929 226  1.20x   50   94 060   4.9    0.52
#
# 0.76 core-hours, ~5 minutes at 9 concurrent. Written to a NEW directory so the 261003/261004
# cells stay exactly as analysed; nothing is overwritten.
set -uo pipefail
cd /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3
P=experiments_energy_transfer/level4_topup_20261005
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
    echo "M=$M seed=$sd FAILED rc=$rc health=$health"; return 1
  fi
  python3 - "$d/raw_${sd}.csv" "$d/red_${sd}.csv" <<'PY'
import sys, pandas as pd
cols = ["Time","KE_gas_left","KE_gas_right","W0_x_sigma"]
pd.read_csv(sys.argv[1], low_memory=False, usecols=cols).to_csv(sys.argv[2], index=False)
PY
  if [ -s "$d/red_${sd}.csv" ]; then rm -f "$d/raw_${sd}.csv"; else
    echo "M=$M seed=$sd reduction FAILED"; return 1; fi
  return 0
}

for M in 20 10; do mkdir -p $P/Md${M}; done
echo "topup start $(date +%H:%M:%S), max $MAXJOBS concurrent"
for spec in "20 4703000 50" "10 2144000 200"; do
  set -- $spec; M=$1; STEPS=$2; EVERY=$3
  for i in $(seq 0 79); do
    while [ "$(jobs -rp | wc -l | tr -d ' ')" -ge "$MAXJOBS" ]; do sleep 3; done
    run_one "$M" "$STEPS" "$EVERY" $((9200+i)) &
  done
done
wait
echo "topup done $(date +%H:%M:%S)"
for M in 10 20; do
  d=$P/Md${M}
  echo "  M=$M: $(find $d -name 'red_*.csv'|wc -l|tr -d ' ')/80; health $(cat $d/run_*.log 2>/dev/null | grep -cE 'EDMD-HEALTH|forced_advance|clamp_repair|overlap_repair|wall_overdue'); aborts $(cat $d/run_*.log 2>/dev/null | grep -c ABORTING)"
done
du -sh $P
