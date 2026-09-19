#!/usr/bin/env bash
# ##CHRIS 2026-09-16: Level 1 statistics. Same commands as level0_Wqs_20260911/run.sh (identical geometry,
# hold, steps per speed, event log), 40 NEW seeds 9100-9139 per speed, disjoint from the original 9000-9009.
# Purpose: decide whether <W>(u = 0.02) lying 3.3 sigma below W_qs^finite was a 10-seed fluctuation.
set -uo pipefail
cd /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3
P="experiments_energy_transfer/level1_moreseeds_20260916"
for spec in "0.005 60000" "0.01 36000" "0.02 24000"; do
  set -- $spec; u=$1; st=$2
  mkdir -p "$P/u${u}"
  for k in $(seq 0 39); do
    sd=$(( 9100 + k )); tj="$P/u${u}/tr_${sd}.csv"; [ -s "$tj" ] && continue
    ( HD_PISTON_EVENTS="$P/u${u}/ev_${sd}.csv" nice -n 5 ./00ALLINONE --mode=edmd --experiment=energy_transfer \
        --headless --quiet --edmd-acc=0 --seed-drift-order=drift-first \
        --energy-transfer-summary="$P/u${u}/summary.csv" --energy-transfer-trace="$tj" \
        --particles=100 --particles-boxes=50,50 --particle-radius=0.5 --l0=39.25 --height=10 \
        --num-walls=1 --wall-positions=39.25 --wall-mass-factors=1000000000 \
        --piston-right-protocol-mode=step --velocity-right-piston-step=$u \
        --max-right-piston-travel=3.93 --auto-piston-step \
        --wall-hold-steps=12000 --steps=$st --fixed-dt=0.4 --energy-measurement \
        --eff-output=wall-ke --kbt1 --seed=$sd >> "$P/u${u}/run.log" 2>&1 ) &
    while [ "$(jobs -rp | wc -l)" -ge 10 ]; do sleep 1; done
  done
done
wait
echo "done: $(ls $P/u*/tr_*.csv 2>/dev/null | wc -l | tr -d ' ')/120 traces; summary rows $(cat $P/u*/summary.csv 2>/dev/null | grep -vc '^timestamp'); aborts $(cat $P/u*/run.log | grep -c ABORTING); health lines $(cat $P/u*/run.log | grep -c EDMD-HEALTH)"
