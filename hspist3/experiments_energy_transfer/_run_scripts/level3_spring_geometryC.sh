#!/usr/bin/env bash
# ##CHRIS 2026-09-19: Level 3 pilot, geometry C -- one gas drives a spring-loaded wall.
# wall | spring | free wall (M_s = 200, k = 5) | gas (N_s = 50, eta = 0.1013) | piston.
# Five speeds spanning both limits: T_w = 39.5 sigma-time against pushes of 786 .. 20.
set -uo pipefail
cd /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3
P=experiments_energy_transfer/level3_spring_20260919
for spec in "0.005 60000 25" "0.02 30000 25" "0.05 20000 25" "0.1 18000 25" "0.2 17000 25"; do
  set -- $spec; u=$1; st=$2; n=$3
  mkdir -p "$P/u$u"
  for k in $(seq 0 $((n-1))); do
    sd=$((9200+k)); [ -s "$P/u$u/tr_${sd}.csv" ] && continue
    ( HD_PISTON_EVENTS="$P/u$u/ev_${sd}.csv" nice -n 5 ./00ALLINONE --mode=edmd \
        --experiment=energy_transfer --headless --quiet --edmd-acc=0 \
        --seed-drift-order=drift-first --energy-transfer-summary="$P/u$u/summary.csv" \
        --energy-transfer-trace="$P/u$u/tr_${sd}.csv" --particles=50 --particles-boxes=0,50 \
        --particle-radius=0.5 --l0=39.25 --height=10 --num-walls=1 --wall-positions=39.25 \
        --wall-mass-factors=200 --wall-thickness=0.05 --spring-k=5 --spring-eq=39.25 \
        --energy-measurement --eff-output=spring --piston-right-protocol-mode=step \
        --velocity-right-piston-step=$u --max-right-piston-travel=3.93 --auto-piston-step \
        --wall-hold-steps=12000 --steps=$st --fixed-dt=0.4 --kbt1 --seed=$sd \
        >> "$P/u$u/run.log" 2>&1 ) &
    while [ "$(jobs -rp | wc -l)" -ge 10 ]; do sleep 0.3; done
  done
done
wait
echo "done $(ls $P/*/tr_*.csv | wc -l | tr -d ' ')/125; aborts $(cat $P/*/run.log | grep -c ABORTING); health $(cat $P/*/run.log | grep -c EDMD-HEALTH)"
