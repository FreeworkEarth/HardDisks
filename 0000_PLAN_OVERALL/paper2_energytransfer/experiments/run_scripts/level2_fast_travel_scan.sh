#!/usr/bin/env bash
# ##CHRIS 2026-09-18: fast end against travel. The flag is the GAS compression (the target is
# box_wall - flag); the piston displaces flag + 0.25 because it parks outside the wall.
set -uo pipefail
cd /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3
F=experiments_energy_transfer/level2_fastdx_20260918
for dx in 2 3 5; do
  mkdir -p "$F/dx$dx"
  for k in $(seq 0 99); do
    sd=$((9200+k)); [ -s "$F/dx$dx/tr_${sd}.csv" ] && continue
    ( HD_PISTON_EVENTS="$F/dx$dx/ev_${sd}.csv" nice -n 5 ./00ALLINONE_ramp \
        --mode=edmd --experiment=energy_transfer --headless --quiet --edmd-acc=0 \
        --seed-drift-order=drift-first --energy-transfer-summary="$F/dx$dx/summary.csv" \
        --energy-transfer-trace="$F/dx$dx/tr_${sd}.csv" --particles=100 --particles-boxes=50,50 \
        --particle-radius=0.5 --l0=39.25 --height=10 --num-walls=1 --wall-positions=39.25 \
        --wall-mass-factors=1000000000 --piston-right-protocol-mode=step \
        --velocity-right-piston-step=10 --max-right-piston-travel=$dx --auto-piston-step \
        --wall-hold-steps=12000 --steps=12500 --fixed-dt=0.4 --energy-measurement \
        --eff-output=wall-ke --kbt1 --seed=$sd >> "$F/dx$dx/run.log" 2>&1 ) &
    while [ "$(jobs -rp | wc -l)" -ge 10 ]; do sleep 0.3; done
  done
done
wait
echo "done $(ls $F/*/tr_*.csv | wc -l | tr -d ' ')/300; aborts $(cat $F/*/run.log | grep -c ABORTING); health $(cat $F/*/run.log | grep -c EDMD-HEALTH)"
