#!/usr/bin/env bash
# ##CHRIS 2026-09-19: Level 0b -- does the master box with both walls held reproduce geometry A?
# Grid-exact (every position an integer number of 1/24 sigma) and t = 1.0, the default
# Levels 0-2 ran with -- which is what makes the gas compartment exactly 38.75.
# Travel 3.93 -- the SAME flag value geometry A used. An earlier attempt used 3.875 ("10 % of
# 38.75") and read 1.08 sigma; the travel, not the box, was the whole difference. Only wall
# positions and the box length must sit on the 1/24 sigma grid; the travel need not.
# 10-sigma compartment behind a held wall, which the gas cannot see. If <W_in> agrees within 1 sigma
# with the 60 existing geometry-A seeds at u = 0.05, Levels 0-2 stand as measured in the master box.
set -uo pipefail
cd /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3
P=experiments_energy_transfer/level0b_masterbox_v3_20260919/u0.05
mkdir -p "$P"
for k in $(seq 0 24); do
  sd=$((9200+k)); [ -s "$P/tr_${sd}.csv" ] && continue
  ( HD_PISTON_EVENTS="$P/ev_${sd}.csv" nice -n 5 ./00ALLINONE_sp --mode=edmd \
      --experiment=energy_transfer --headless --quiet --edmd-acc=0 --seed-drift-order=drift-first \
      --energy-transfer-summary="$P/summary.csv" --energy-transfer-trace="$P/tr_${sd}.csv" \
      --particles=100 --particles-boxes=0,50,50 --particle-radius=0.5 --l0=44.75 --height=10 \
      --num-walls=2 --wall-positions=10.5,50.25 --wall-mass-factors=1000000000,1000000000 \
      --piston-right-protocol-mode=step --velocity-right-piston-step=0.05 \
      --max-right-piston-travel=3.93 --auto-piston-step --wall-hold-steps=12000 --steps=17000 \
      --fixed-dt=0.4 --energy-measurement --eff-output=wall-ke --kbt1 --seed=$sd \
      >> "$P/run.log" 2>&1 ) &
  while [ "$(jobs -rp | wc -l)" -ge 8 ]; do sleep 0.3; done
done
wait
echo "done $(ls $P/tr_*.csv | wc -l | tr -d ' ')/25; aborts $(grep -c ABORTING $P/run.log); health $(grep -c EDMD-HEALTH $P/run.log)"
