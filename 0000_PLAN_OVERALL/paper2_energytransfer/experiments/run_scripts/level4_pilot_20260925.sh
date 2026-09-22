#!/usr/bin/env bash
# ##CHRIS 2026-09-25: Level 4 PILOT -- geometry B, transmission through a free divider.
#
# Master box, both gases present, wall_S HELD (it is only the spring compartment's boundary here),
# divider at 70.25 FREE with mass M_d. Piston compresses gas 2 by 10 % (dx = 3.93 sigma on a
# 38.75 sigma compartment, eta 0.101 -> 0.113) at u = 0.05.
#
# Records are set to >= 5 tau_heat(M_d) using the values A4 measured in the spring geometry
# (tau_heat = 48 sigma-time at M = 50, 154 at M = 200), i.e. 250 and 800 sigma-time after the push.
# This is a PILOT: it fixes the protocol, it does not decide a pass. The KOA ladder is
# M_d in {10, 50, 200, 1000, 5000}, 40 seeds, records 5 tau_heat.
set -uo pipefail
cd /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3
P=experiments_energy_transfer/level4_pilot_20260925
for spec in "50 20000" "200 55000"; do
  set -- $spec; md=$1; st=$2
  d=$P/Md${md}; mkdir -p "$d"
  for i in $(seq 0 19); do
    sd=$((9200+i)); [ -s "$d/tr_${sd}.csv" ] && continue
    ( HD_PISTON_EVENTS="$d/ev_${sd}.csv" nice -n 5 ./00ALLINONE --mode=edmd \
        --experiment=energy_transfer --headless --quiet --edmd-acc=0 \
        --seed-drift-order=drift-first --energy-transfer-summary="$d/summary.csv" \
        --energy-transfer-trace="$d/tr_${sd}.csv" \
        --particles=100 --particles-boxes=0,50,50 --particle-radius=0.5 \
        --l0=54.75 --height=10 --num-walls=2 --wall-positions=30.5,70.25 \
        --wall-mass-factors=1000000000,$md --eff-output=wall-ke \
        --piston-right-protocol-mode=step --velocity-right-piston-step=0.05 \
        --max-right-piston-travel=3.93 --auto-piston-step \
        --wall-hold-steps=12000 --steps=$st --fixed-dt=0.4 --kbt1 --seed=$sd \
        >> "$d/run.log" 2>&1 ) &
    while [ "$(jobs -rp | wc -l)" -ge 9 ]; do sleep 0.3; done
  done; wait
done
echo "level4 pilot done $(ls $P/*/tr_*.csv | wc -l | tr -d ' ')/40; aborts $(cat $P/*/run.log | grep -c ABORTING); health $(cat $P/*/run.log | grep -cE 'EDMD-HEALTH|forced_advance|clamp_repair|overlap_repair|wall_overdue')"
