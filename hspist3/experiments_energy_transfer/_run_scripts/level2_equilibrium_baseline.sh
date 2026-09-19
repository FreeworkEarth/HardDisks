#!/usr/bin/env bash
# ##CHRIS 2026-09-18: equilibrium reference profile at the SAME final compartment. A push so slow
# that the excess work is zero leaves the gas in equilibrium at L_f, so its stop snapshot is the
# baseline the flow/compression decomposition must be measured against (wall layering included).
set -uo pipefail
cd /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3
P=experiments_energy_transfer/level2_ramp_20260918/step_u0.005
mkdir -p $P
for k in $(seq 0 59); do
  sd=$((9200+k)); [ -s "$P/tr_${sd}.csv" ] && continue
  ( HD_PISTON_EVENTS="$P/ev_${sd}.csv" HD_STOP_SNAPSHOT="$P/snap_${sd}.csv" nice -n 5 ./00ALLINONE_ramp \
      --mode=edmd --experiment=energy_transfer --headless --quiet --edmd-acc=0 \
      --seed-drift-order=drift-first --energy-transfer-summary="$P/summary.csv" \
      --energy-transfer-trace="$P/tr_${sd}.csv" --particles=100 --particles-boxes=50,50 \
      --particle-radius=0.5 --l0=39.25 --height=10 --num-walls=1 --wall-positions=39.25 \
      --wall-mass-factors=1000000000 --piston-right-protocol-mode=step --velocity-right-piston-step=0.005 \
      --max-right-piston-travel=3.93 --auto-piston-step --wall-hold-steps=12000 --steps=60000 \
      --fixed-dt=0.4 --energy-measurement --eff-output=wall-ke --kbt1 --seed=$sd >> "$P/run.log" 2>&1 ) &
  while [ "$(jobs -rp | wc -l)" -ge 10 ]; do sleep 0.3; done
done
wait
echo "done $(ls $P/snap_*.csv | wc -l | tr -d ' ')/60; aborts $(grep -c ABORTING $P/run.log); health $(grep -c EDMD-HEALTH $P/run.log)"
