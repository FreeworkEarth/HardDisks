#!/usr/bin/env bash
# ##CHRIS 2026-09-26: Level 4 on the Mac. Record lengths set by the B0 prediction, not by guesswork.
#
# B0 (see the plan's Level 4 section): the Brownian-piston friction derived from kinetic theory is
# gamma = 4 n h sqrt(2 m kT/pi), so the VELOCITY relaxes on tau_v = M/gamma = M/4.118 here. The
# TEMPERATURE relaxation does not inherit that: measured/tau_v is 3.7 at M = 50 and 13.8 at M = 200,
# i.e. the ratio itself grows as M, giving tau_T ~ M^2. Predictions used below:
#   M =  10 : the M^2 law has saturated (tau_T would be 1.8 < tau_v = 2.4), predict tau_T ~ tau_v ~ 2.4
#   M =  50 : tau_T ~ 45.5 sigma-time
#   M = 200 : tau_T ~ 729 sigma-time
# Records are >= 5 tau_T and never shorter than 5 x 46 at M = 50, as instructed.
#   M = 10, 50 : 1000 sigma-time (60000 steps), 40 seeds
#   M = 200    : 4000 sigma-time (240000 steps), 8 seeds -- to MEASURE tau(200) rather than bound it
#
# Two-compartment box (the Level 4 pilot's protocol fix): segment 0 IS gas 1, segment 1 IS gas 2,
# so KE_gas_left/right are the two gases. Box 0..78.5, divider centre 39.25, 38.75 sigma each side.
set -uo pipefail
cd /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3
P=experiments_energy_transfer/level4_v3_20260926
for spec in "10 60000 40" "50 60000 40" "200 240000 8"; do
  set -- $spec; md=$1; st=$2; ns=$3
  d=$P/Md${md}; mkdir -p "$d"
  for i in $(seq 0 $((ns-1))); do
    sd=$((9200+i)); [ -s "$d/tr_${sd}.csv" ] && continue
    ( HD_PISTON_EVENTS="$d/ev_${sd}.csv" nice -n 5 ./00ALLINONE --mode=edmd \
        --experiment=energy_transfer --headless --quiet --edmd-acc=0 \
        --seed-drift-order=drift-first --energy-transfer-summary="$d/summary.csv" \
        --energy-transfer-trace="$d/tr_${sd}.csv" \
        --particles=100 --particles-boxes=50,50 --particle-radius=0.5 \
        --l0=39.25 --height=10 --num-walls=1 --wall-positions=39.25 \
        --wall-mass-factors=$md --eff-output=wall-ke \
        --piston-right-protocol-mode=step --velocity-right-piston-step=0.05 \
        --max-right-piston-travel=3.93 --auto-piston-step \
        --wall-hold-steps=12000 --steps=$st --fixed-dt=0.4 --kbt1 --seed=$sd \
        >> "$d/run.log" 2>&1 ) &
    while [ "$(jobs -rp | wc -l)" -ge 9 ]; do sleep 0.3; done
  done; wait
done
echo "level4 v3 done $(ls $P/*/tr_*.csv | wc -l | tr -d ' ')/88; aborts $(cat $P/*/run.log | grep -c ABORTING); health $(cat $P/*/run.log | grep -cE 'EDMD-HEALTH|forced_advance|clamp_repair|overlap_repair|wall_overdue')"
