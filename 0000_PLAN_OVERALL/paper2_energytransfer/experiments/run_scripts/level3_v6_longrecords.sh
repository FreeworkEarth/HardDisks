#!/usr/bin/env bash
# ##CHRIS 2026-09-24: v6 item 6 -- re-run the two cells whose records were too short for a valid
# settled window. The settled average needs >= 3 free periods after t > tau + 3 T_w, and
# T_w = 2 pi sqrt(M_s/(k + k_gas)) = 120 sigma-time at M_s = 200. The v3 records ran to 500, which
# left 0.5 period at u = 0.1 and nothing at all at u = 0.05. 75000 steps = 1250 sigma-time.
# Geometry, pre-load and seeds identical to v3; only --steps changes.
set -uo pipefail
cd /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3
P=experiments_energy_transfer/level3_v6_20260924
for u in 0.05 0.1; do
  d=$P/k0.5_M200_u${u}; mkdir -p "$d"
  for i in $(seq 0 39); do
    sd=$((9200+i)); [ -s "$d/tr_${sd}.csv" ] && continue
    ( HD_PISTON_EVENTS="$d/ev_${sd}.csv" nice -n 5 ./00ALLINONE --mode=edmd \
        --experiment=energy_transfer --headless --quiet --edmd-acc=0 \
        --seed-drift-order=drift-first --energy-transfer-summary="$d/summary.csv" \
        --energy-transfer-trace="$d/tr_${sd}.csv" --particles=100 --particles-boxes=0,100 \
        --particle-radius=0.5 --l0=54.75 --height=10 --num-walls=1 --wall-positions=30.5 \
        --wall-mass-factors=200 --spring-k-sigma=0.5 --spring-wall=0 --spring-eq=33.65 \
        --eff-output=spring --piston-right-protocol-mode=step \
        --velocity-right-piston-step=$u --max-right-piston-travel=7.96 --auto-piston-step \
        --wall-hold-steps=12000 --steps=75000 --fixed-dt=0.4 --kbt1 --seed=$sd \
        >> "$d/run.log" 2>&1 ) &
    while [ "$(jobs -rp | wc -l)" -ge 9 ]; do sleep 0.3; done
  done; wait
done
echo "v6 done $(ls $P/*/tr_*.csv | wc -l | tr -d ' ')/80; aborts $(cat $P/*/run.log | grep -c ABORTING); health $(cat $P/*/run.log | grep -cE 'EDMD-HEALTH|forced_advance|clamp_repair|overlap_repair|wall_overdue')"
