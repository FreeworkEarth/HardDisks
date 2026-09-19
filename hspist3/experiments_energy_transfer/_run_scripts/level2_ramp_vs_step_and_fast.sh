#!/usr/bin/env bash
# ##CHRIS 2026-09-18: Level 2 ramp-vs-step (item 1) and fast end (item 3). Separate ramp build,
# installed 00ALLINONE untouched. Same box/hold/flags as level2_slope_20260917 so every speed is
# comparable. Stop snapshots on for the flow/compression decomposition.
set -uo pipefail
cd /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3
B=./00ALLINONE_ramp
run () {  # $1 dir  $2 mode  $3 u  $4 steps  $5 nseeds  $6 extra
  mkdir -p "$1"
  for k in $(seq 0 $(($5-1))); do
    sd=$((9200+k)); tj="$1/tr_${sd}.csv"; [ -s "$tj" ] && continue
    ( HD_PISTON_EVENTS="$1/ev_${sd}.csv" HD_STOP_SNAPSHOT="$1/snap_${sd}.csv" nice -n 5 $B \
        --mode=edmd --experiment=energy_transfer --headless --quiet --edmd-acc=0 \
        --seed-drift-order=drift-first --energy-transfer-summary="$1/summary.csv" \
        --energy-transfer-trace="$tj" --particles=100 --particles-boxes=50,50 --particle-radius=0.5 \
        --l0=39.25 --height=10 --num-walls=1 --wall-positions=39.25 --wall-mass-factors=1000000000 \
        --piston-right-protocol-mode=$2 $6 --velocity-right-piston-step=$3 \
        --max-right-piston-travel=3.93 --auto-piston-step --wall-hold-steps=12000 --steps=$4 \
        --fixed-dt=0.4 --energy-measurement --eff-output=wall-ke --kbt1 --seed=$sd >> "$1/run.log" 2>&1 ) &
    while [ "$(jobs -rp | wc -l)" -ge 10 ]; do sleep 0.3; done
  done
}
P=experiments_energy_transfer/level2_ramp_20260918
F=experiments_energy_transfer/level2_fast_20260918
# item 1: ramp, T = 40 sigma-time (clamped to travel/u where that is shorter), 60 seeds each
run "$P/u0.05" ramp 0.05 17000 60 "--piston-ramp-time=40"
run "$P/u0.10" ramp 0.10 15000 60 "--piston-ramp-time=40"
run "$P/u0.20" ramp 0.20 13500 60 "--piston-ramp-time=40"
# and the step control with snapshots, so the decomposition can be done on both
run "$P/step_u0.05" step 0.05 17000 60 ""
run "$P/step_u0.10" step 0.10 15000 60 ""
run "$P/step_u0.20" step 0.20 13500 60 ""
wait
# item 3: fast end, travel 1 sigma (flag is distance-from-right-wall; 3.93 -> 4.18 achieved, so
# 6.68 gives 1.18 and 6.93 gives 0.93; the achieved travel is measured per run, not assumed)
for spec in "3 12500" "5 12500" "10 12500"; do
  set -- $spec
  mkdir -p "$F/u$1"
  for k in $(seq 0 99); do
    sd=$((9200+k)); tj="$F/u$1/tr_${sd}.csv"; [ -s "$tj" ] && continue
    ( HD_PISTON_EVENTS="$F/u$1/ev_${sd}.csv" HD_STOP_SNAPSHOT="$F/u$1/snap_${sd}.csv" nice -n 5 $B \
        --mode=edmd --experiment=energy_transfer --headless --quiet --edmd-acc=0 \
        --seed-drift-order=drift-first --energy-transfer-summary="$F/u$1/summary.csv" \
        --energy-transfer-trace="$tj" --particles=100 --particles-boxes=50,50 --particle-radius=0.5 \
        --l0=39.25 --height=10 --num-walls=1 --wall-positions=39.25 --wall-mass-factors=1000000000 \
        --piston-right-protocol-mode=step --velocity-right-piston-step=$1 \
        --max-right-piston-travel=6.93 --auto-piston-step --wall-hold-steps=12000 --steps=$2 \
        --fixed-dt=0.4 --energy-measurement --eff-output=wall-ke --kbt1 --seed=$sd >> "$F/u$1/run.log" 2>&1 ) &
    while [ "$(jobs -rp | wc -l)" -ge 10 ]; do sleep 0.3; done
  done
done
wait
echo "ramp+step $(ls $P/*/tr_*.csv 2>/dev/null | wc -l | tr -d ' ')/360, fast $(ls $F/*/tr_*.csv 2>/dev/null | wc -l | tr -d ' ')/300"
echo "aborts $(cat $P/*/run.log $F/*/run.log | grep -c ABORTING), health $(cat $P/*/run.log $F/*/run.log | grep -c EDMD-HEALTH)"
