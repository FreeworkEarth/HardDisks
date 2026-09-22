#!/usr/bin/env bash
# ##CHRIS 2026-09-25: measure F(L) directly, so the fixed point needs no equation of state.
#
# Geometry C with the spring wall HELD at 1e9. The piston compresses the gas by 0 / 2.5 / 5 / 7.5 /
# 10 % and then stops; after the transient dies the gas is in equilibrium at that length and the
# momentum delivered to the held wall per unit time IS the force. Measuring the gas temperature at
# the same time (KE_gas_total; in 2D U = N kT so T = KE/N) gives Z_box = F L /(N k T) with no
# adiabat and no area convention assumed.
#
# Why: the v6 comparison rested on the bulk EOS evaluated with the nominal area, and the residual
# (+3.6 %) is the finite-size stiffness of the 100-disk box. This measures that stiffness instead of
# inferring it. 20 seeds per compression; the 0 % cell simply never starts the piston.
set -uo pipefail
cd /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3
P=experiments_energy_transfer/level3_FofL_20260925
for spec in "0 0" "1.99 2.5" "3.98 5" "5.97 7.5" "7.96 10"; do
  set -- $spec; dx=$1; pct=$2
  d=$P/c${pct}; mkdir -p "$d"
  EXTRA="--max-right-piston-travel=$dx --auto-piston-step"
  [ "$dx" = "0" ] && EXTRA=""          # 0 %: the piston never starts
  for i in $(seq 0 19); do
    sd=$((9200+i)); [ -s "$d/tr_${sd}.csv" ] && continue
    ( HD_PISTON_EVENTS="$d/ev_${sd}.csv" nice -n 5 ./00ALLINONE --mode=edmd \
        --experiment=energy_transfer --headless --quiet --edmd-acc=0 \
        --seed-drift-order=drift-first --energy-transfer-summary="$d/summary.csv" \
        --energy-transfer-trace="$d/tr_${sd}.csv" \
        --particles=100 --particles-boxes=0,100 --particle-radius=0.5 \
        --l0=54.75 --height=10 --num-walls=1 --wall-positions=30.5 \
        --wall-mass-factors=1000000000 --eff-output=wall-ke \
        --piston-right-protocol-mode=step --velocity-right-piston-step=0.05 $EXTRA \
        --wall-hold-steps=12000 --steps=40000 --fixed-dt=0.4 --kbt1 --seed=$sd \
        >> "$d/run.log" 2>&1 ) &
    while [ "$(jobs -rp | wc -l)" -ge 9 ]; do sleep 0.3; done
  done; wait
done
echo "F(L) done $(ls $P/*/tr_*.csv | wc -l | tr -d ' ')/100; aborts $(cat $P/*/run.log | grep -c ABORTING); health $(cat $P/*/run.log | grep -cE 'EDMD-HEALTH|forced_advance|clamp_repair|overlap_repair|wall_overdue')"
