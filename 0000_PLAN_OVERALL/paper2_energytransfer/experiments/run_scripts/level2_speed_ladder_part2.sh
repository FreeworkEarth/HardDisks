#!/usr/bin/env bash
# ##CHRIS 2026-09-17: Level 2, slow end. Test dW_diss/du = zeta*Delta_x with zeta measured independently
# from equilibrium force fluctuations (2.53 +- 0.05 -> predicted slope 9.9 +- 0.2).
# Same geometry, hold and flags as the Level-1 runs (hold 12000 steps), so every speed is comparable:
# u = 0.02 (+50 seeds -> 100), 0.03 (100, new), 0.05 (+90 -> 100), 0.10 (+90 -> 100). Seeds 9200+.
set -uo pipefail
cd /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3
P="experiments_energy_transfer/level2_slope_20260917"
for spec in "0.005 60000 90" "0.01 36000 90" "0.15 14000 100" "0.20 13500 100"; do
  set -- $spec; u=$1; st=$2; n=$3
  mkdir -p "$P/u${u}"
  for k in $(seq 0 $((n-1))); do
    sd=$(( 9200 + k )); tj="$P/u${u}/tr_${sd}.csv"; [ -s "$tj" ] && continue
    ( HD_PISTON_EVENTS="$P/u${u}/ev_${sd}.csv" nice -n 5 ./00ALLINONE --mode=edmd --experiment=energy_transfer \
        --headless --quiet --edmd-acc=0 --seed-drift-order=drift-first \
        --energy-transfer-summary="$P/u${u}/summary.csv" --energy-transfer-trace="$tj" \
        --particles=100 --particles-boxes=50,50 --particle-radius=0.5 --l0=39.25 --height=10 \
        --num-walls=1 --wall-positions=39.25 --wall-mass-factors=1000000000 \
        --piston-right-protocol-mode=step --velocity-right-piston-step=$u \
        --max-right-piston-travel=3.93 --auto-piston-step \
        --wall-hold-steps=12000 --steps=$st --fixed-dt=0.4 --energy-measurement \
        --eff-output=wall-ke --kbt1 --seed=$sd >> "$P/u${u}/run.log" 2>&1 ) &
    while [ "$(jobs -rp | wc -l)" -ge 10 ]; do sleep 0.5; done
  done
done
wait
echo "done: $(ls $P/u*/tr_*.csv 2>/dev/null | wc -l | tr -d ' ')/380; aborts $(cat $P/u*/run.log | grep -c ABORTING); health $(cat $P/u*/run.log | grep -c EDMD-HEALTH)"
