#!/usr/bin/env bash
# ##CHRIS 2026-10-05: the R-COLLAPSE test. Is the variable R = N_s m/M, or M itself?
#
# WHY THIS IS THE RIGHT TEST. Cencini et al. 2007 p. 4, Sect. II.B, define their limit explicitly:
#   "We are interested in the limit N, M, L -> infinity in which we keep fixed rho_0 = N/L and the
#    nondimensional mass ratio R = Nm/M."
# Their "timescale proportional to M", and Gruber-Piasecki's infinite-gas derivation, are statements
# AT FIXED R. Our 261003/261005 ladder held N_s = 50 and L fixed and varied M, so R ran 5 -> 0.25;
# nothing in either paper claims linearity along that path. Written as tau_T = M g(R), the ladder is
#   R    = 5.00  2.50  1.00  0.50  0.25
#   g    = 47.5  54.4  76.7 123.3 200.4   (+- 2.0, 2.8, 1.9, 8.6, 8.5)
# and g(R -> infinity) is the Gruber-Piasecki value 50.35 -- which R = 5 reproduces at -1.4 sigma
# and R = 2.5 at +1.4 sigma. The theory is CONFIRMED where the gas is the reservoir and departs,
# smoothly and by a factor 4, once the divider outweighs a compartment of gas.
#
# THE TEST. Double N_s at fixed eta and fixed rho_0 = N_s/L_c, which is the path Cencini's limit is
# defined along, and ask which prediction the new points land on:
#   prediction A -- the variable is R: g(R=2) = 59.2, g(R=1) = 76.7  (interpolated from the ladder)
#   prediction B -- the variable is M: g(M=50) = 76.7, g(M=100) = 123.3
# Separation: 30 % at M = 50 and 61 % at M = 100, against a 5-16 % estimator spread.
#
# GEOMETRY. L_c = 2 x 38.75 = 77.5 EXACTLY -- doubling N_s and L_c together leaves eta identical to
# 8 digits (0.10134170) and is grid-exact (77.5 x 24 = 1860, box 156.0 x 24 = 3744, centre 78.0 x 24
# = 1872). Computing L_c from a rounded eta instead gave 77.4998 and FAILED the grid check.
# rho_0 = N_s/L_c = 1.2903, the same as the ladder. t = 1.0 as recorded everywhere in Level 4.
# tau_GP doubles with L_c to 100.70 M.
#
#   M   R  alpha   K      period  tau_GP  predA  predB  record     steps       N    s/per  MB/seed
#   50  2  0.25  1.2646   217.3    5035   2959   3835    327 269   19 766 000  590   22.1   1.7
#  100  1  0.50  1.0769   255.2   10070   7670  12330    801 450   48 231 000  690   22.2   3.6
#
# Records are 65 x the LARGEST candidate per cell, so every hypothesis gets 65 tau.
set -uo pipefail
cd /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3
P=experiments_energy_transfer/level4_Rcollapse_20261005
MAXJOBS=${MAXJOBS:-9}

run_one() {
  local M=$1 STEPS=$2 EVERY=$3 sd=$4
  local d=$P/Md${M}
  [ -s "$d/red_${sd}.csv" ] && return 0
  nice -n 5 ./00ALLINONE --mode=edmd --experiment=energy_transfer --headless --quiet \
      --edmd-acc=0 --seed-drift-order=drift-first \
      --energy-transfer-summary="$d/summary_${sd}.csv" \
      --energy-transfer-trace="$d/raw_${sd}.csv" --trace-every=$EVERY \
      --particles=200 --particles-boxes=100,100 --particle-radius=0.5 \
      --l0=78.0 --height=10 --num-walls=1 --wall-positions=78.0 \
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

for M in 100 50; do mkdir -p $P/Md${M}; done
echo "Rcollapse start $(date +%H:%M:%S), max $MAXJOBS concurrent"
for spec in "100 48231000 690" "50 19766000 590"; do
  set -- $spec; M=$1; STEPS=$2; EVERY=$3
  for i in $(seq 0 79); do
    while [ "$(jobs -rp | wc -l | tr -d ' ')" -ge "$MAXJOBS" ]; do sleep 5; done
    run_one "$M" "$STEPS" "$EVERY" $((9200+i)) &
  done
done
wait
echo "Rcollapse done $(date +%H:%M:%S)"
for M in 50 100; do
  d=$P/Md${M}
  echo "  M=$M: $(find $d -name 'red_*.csv'|wc -l|tr -d ' ')/80; health $(cat $d/run_*.log 2>/dev/null | grep -cE 'EDMD-HEALTH|forced_advance|clamp_repair|overlap_repair|wall_overdue'); aborts $(cat $d/run_*.log 2>/dev/null | grep -c ABORTING)"
done
du -sh $P
