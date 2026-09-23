#!/usr/bin/env bash
# ##CHRIS 2026-09-28: Level 4 THERMAL -- the adiabatic piston with no piston at all.
#
# Two gases at different temperatures either side of a free divider; nothing else happens. This is
# the Gruber-Piasecki / Cencini stage-2 measurement, and it is the configuration their theory is
# actually about.
#
# NO NEW FLAG WAS NEEDED. --temperature accepts a COMMA-SEPARATED LIST, one value per segment
# (00ALLINONE.c:4579 -> cli_parse_float_list -> apply_segment_temperatures at :6271, applied right
# after the per-segment equalise during initialisation). An earlier report that this could not be
# done was wrong: it read the help text, which documents only the scalar form.
#
# UNITS TRAP, verified: with --kbt1 the runtime temperature is 100 and kB is scaled so kB*T = 1, so
# the LIST is in those raw units. --temperature=1.25,0.75 gives kT = 0.0125/0.0075 -- the right
# ratio and 80x the wrong scale. Use --temperature=125,75 for kT = 1.25/0.75. Measured: kT_left =
# 1.2500, kT_right = 0.7500, ratio 1.6667.
#
# Records >= 5 tau_T with tau_T = 50.4 M (hard-disk isobar): 504 and 2517 sigma-time at M_d = 10
# and 50, so 3000 and 14000 sigma-time. At 0.016643 sigma-time/step that is 181000 and 842000 steps.
set -uo pipefail
cd /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3
P=experiments_energy_transfer/level4_thermal_20260928
for spec in "10 181000" "50 842000"; do
  set -- $spec; md=$1; st=$2
  d=$P/Md${md}; mkdir -p "$d"
  for i in $(seq 0 19); do
    sd=$((9200+i)); [ -s "$d/tr_${sd}.csv" ] && continue
    ( nice -n 5 ./00ALLINONE --mode=edmd --experiment=energy_transfer --headless --quiet \
        --edmd-acc=0 --seed-drift-order=drift-first \
        --energy-transfer-summary="$d/summary.csv" --energy-transfer-trace="$d/tr_${sd}.csv" \
        --particles=100 --particles-boxes=50,50 --particle-radius=0.5 \
        --l0=39.25 --height=10 --num-walls=1 --wall-positions=39.25 \
        --wall-mass-factors=$md --eff-output=wall-ke \
        --wall-hold-steps=12000 --steps=$st --fixed-dt=0.4 --kbt1 \
        --temperature=125,75 --seed=$sd >> "$d/run.log" 2>&1 ) &
    while [ "$(jobs -rp | wc -l)" -ge 9 ]; do sleep 0.3; done
  done; wait
done
echo "level4 thermal done $(ls $P/*/tr_*.csv | wc -l | tr -d ' ')/40; aborts $(cat $P/*/run.log | grep -c ABORTING); health $(cat $P/*/run.log | grep -cE 'EDMD-HEALTH|forced_advance|clamp_repair|overlap_repair|wall_overdue')"
