#!/usr/bin/env bash
# no-fuse check: A1v2 eta=0.1122 cell, 9 masses x 3 seeds, canonical seeds, argv_for() shape.
set -uo pipefail
cd /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3
N=experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/nofuse_check_20261010/eta_0p112200
run_mass() { local M=$1 T=$2 S=$3; shift 3; local d=$N/m_$M; mkdir -p "$d"; local r=0
  for sd in "$@"; do
    [ -s "$d/wall_x_positions_L0_349_wallmassfactor_${M}_run${r}.csv" ] && { r=$((r+1)); continue; }
    local tmp=$d/.run$r; rm -rf "$tmp"; mkdir -p "$tmp"
    ./00ALLINONE --mode=edmd --experiment=speed_of_sound --headless --kbt1 --seed-drift-order=drift-first --edmd-acc=0 \
      --particles=100 --particles-boxes=50,50 --height=10.0 --particle-radius=0.5 --wall-thickness=0.05 --wall-thickness-vis=0.05 \
      --lengths=34.9999 --wall-masses=$M --repeats=1 --seed=$sd --wall-hold-steps=2000 --fixed-dt=0.4 \
      --target-oscillations=$T --oscillation-safety=1.0 --oscillation-min-steps=10000 --oscillation-max-steps=400000000 \
      --speed-sound-log-stride=$S --speed-sound-run-dir="$tmp" --speed-sound-exact-seed=$sd >> "$d/run.log" 2>&1
    f=$(ls "$tmp"/wall_x_positions_*_run0.csv 2>/dev/null | head -1)
    [ -n "$f" ] && cp "$f" "$d/wall_x_positions_L0_349_wallmassfactor_${M}_run${r}.csv"
    r=$((r+1))
  done
}
echo "nofuse check start $(date +%H:%M:%S)"
run_mass 50 200 208 1028397161 3821645941 349093933 &
run_mass 100 200 260 1351329399 252923647 155003720 &
run_mass 200 200 343 2767689153 2570411904 2185563206 &
run_mass 300 200 409 1293110360 3241107528 1615635956 &
run_mass 500 200 517 1839238817 3086486408 71037639 &
run_mass 750 200 627 1793676155 42399776 697831531 &
run_mass 1000 200 720 512098950 2638603432 1519799425 &
run_mass 1500 200 877 2719783013 3277081129 1292647052 &
run_mass 2000 200 1011 3259447140 3952854748 625024150 &
wait
echo "nofuse check done $(date +%H:%M:%S)"
for M in 50 100 200 300 500 750 1000 1500 2000; do echo "  m_$M: $(ls $N/m_$M/wall_x_positions_*run*.csv 2>/dev/null | wc -l | tr -d ' ')/3  health $(grep -cE 'forced_advance=[1-9]|clamp_repair=[1-9]|overlap_repair=[1-9]|wall_overdue=[1-9]' $N/m_$M/run.log 2>/dev/null)"; done
