#!/usr/bin/env bash
# ##CHRIS 2026-09-19: item 3 -- measure the finite-box wall pressure at PAPER 1's dense geometry
# (eta = 0.5236, L0 = 7.5, N = 100, t = 0.05) instead of extrapolating the eta = 0.10 value.
# Hold only: --steps=1 means the piston never moves, so this is pure equilibrium at that density.
set -uo pipefail
cd /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3
for spec in "7.5 eta0p5236" "20.0 eta0p1963" "10.0 eta0p3927"; do
  set -- $spec; L=$1; tag=$2
  P="experiments_energy_transfer/level1_Zwall_paper1geom_20260919/$tag"
  mkdir -p "$P"
  for k in $(seq 0 39); do
    sd=$((9200+k)); [ -s "$P/ev_${sd}.csv" ] && continue
    ( HD_PISTON_EVENTS="$P/ev_${sd}.csv" nice -n 5 ./00ALLINONE --mode=edmd \
        --experiment=energy_transfer --headless --quiet --edmd-acc=0 \
        --seed-drift-order=drift-first --energy-transfer-summary="$P/summary.csv" \
        --energy-transfer-trace="$P/tr_${sd}.csv" --particles=100 --particles-boxes=50,50 \
        --particle-radius=0.5 --l0=$L --height=10 --num-walls=1 --wall-positions=$L \
        --wall-mass-factors=1000000000 --wall-thickness=0.05 --wall-thickness-vis=0.05 \
        --wall-hold-steps=12000 --steps=1 --fixed-dt=0.4 --energy-measurement \
        --eff-output=wall-ke --kbt1 --seed=$sd >> "$P/run.log" 2>&1 ) &
    while [ "$(jobs -rp | wc -l)" -ge 10 ]; do sleep 0.3; done
  done
done
wait
P=experiments_energy_transfer/level1_Zwall_paper1geom_20260919
echo "done $(ls $P/*/ev_*.csv | wc -l | tr -d ' ')/120; aborts $(cat $P/*/run.log | grep -c ABORTING); health $(cat $P/*/run.log | grep -c EDMD-HEALTH)"
