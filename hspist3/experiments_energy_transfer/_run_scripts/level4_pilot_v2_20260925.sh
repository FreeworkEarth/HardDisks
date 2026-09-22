#!/usr/bin/env bash
# ##CHRIS 2026-09-25: Level 4 pilot v2 -- the PROTOCOL FIX found by pilot v1.
#
# v1 ran geometry B inside the master box and could not measure its own central observable. The
# trace's KE_gas_left is segment 0 ONLY (00ALLINONE.c:17137) and KE_gas_right is every other
# segment, so with the master box's empty spring compartment as segment 0 the columns read
# "0" and "both gases combined" -- the work split across the divider is invisible.
#
# Fix, no code change: drop the empty spring compartment for Level 4 and use the TWO-compartment
# box, where segment 0 IS gas 1 and segment 1 IS gas 2. Level 0b already established that an empty
# compartment behind a held wall changes nothing the gas can see, so this costs no continuity.
#   box 0 .. 78.5 (--l0=39.25), divider centre 39.25, thickness 1.0
#   gas 1: 0 .. 38.75  |  gas 2: 39.75 .. 78.5    (38.75 sigma each, 50 disks, eta = 0.1013)
# Piston compresses gas 2 by 3.93 sigma (10 %) at u = 0.05. Records >= 5 tau_heat(M_d).
set -uo pipefail
cd /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3
P=experiments_energy_transfer/level4_pilot_v2_20260925
for spec in "50 20000" "200 55000"; do
  set -- $spec; md=$1; st=$2
  d=$P/Md${md}; mkdir -p "$d"
  for i in $(seq 0 19); do
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
echo "level4 v2 done $(ls $P/*/tr_*.csv | wc -l | tr -d ' ')/40; aborts $(cat $P/*/run.log | grep -c ABORTING); health $(cat $P/*/run.log | grep -cE 'EDMD-HEALTH|forced_advance|clamp_repair|overlap_repair|wall_overdue')"
