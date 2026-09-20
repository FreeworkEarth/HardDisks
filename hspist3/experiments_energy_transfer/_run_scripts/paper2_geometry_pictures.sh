#!/usr/bin/env bash
# ##CHRIS 2026-09-19: two pictures per geometry from the ONE master box -- eight PNGs.
# This opens a real SDL window eight times; it needs a screen and it will steal focus.
# Usage:  ./paper2_geometry_pictures.sh [outdir]      (default: ../paper2_pictures)
#
# Per geometry:
#   --render=experiment   the normal black GUI, unchanged -- what you watch a run in
#   --render=paper        the white figure render -- what goes in the paper
# Same seed, same flags, same moment mid-push, so the pair differs only in how it is drawn.
#
# The master box (grid-exact: every wall position an integer number of 1/24 sigma):
#   0 .. 109.50 (--l0=54.75), wall_S centre 30.5, divider centre 70.25, thickness 1.0 (default).
#   empty spring compartment 0 .. 30.0 | gas 1 31.0 .. 69.75 | gas 2 70.75 .. 109.50, 38.75 each.
set -uo pipefail
HS=/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3
OUT=${1:-$HS/experiments_energy_transfer/paper2_pictures}
BIN=$HS/00ALLINONE
"$BIN" --help 2>&1 | grep -q -- '--spring-wall' || BIN=$HS/00ALLINONE_sp
mkdir -p "$OUT"
cd "$HS" || exit 1

# --demo runs 2 steps per frame (the GUI is otherwise ~1 step per 25 frames, far too slow to watch).
COMMON=(--show-simulation --demo --mode=edmd --edmd-acc=0 --seed-drift-order=drift-first
        --particle-radius=0.5 --l0=54.75 --height=10
        --piston-right-protocol-mode=step --velocity-right-piston-step=0.02
        --wall-hold-steps=2000 --steps=40000 --fixed-dt=0.4 --energy-measurement --kbt1 --seed=9200)
GAS2=(--particles=100 --particles-boxes=0,50,50 --num-walls=2 --wall-positions=30.5,70.25)
GAS1=(--particles=100 --particles-boxes=0,100  --num-walls=1 --wall-positions=30.5)
HELD=1000000000

A=("${GAS2[@]}" --wall-mass-factors=$HELD,$HELD --eff-output=wall-ke --max-right-piston-travel=3.93)
B=("${GAS2[@]}" --wall-mass-factors=$HELD,1000 --eff-output=wall-ke --max-right-piston-travel=3.93)
C=("${GAS1[@]}" --wall-mass-factors=200 --spring-k-sigma=0.494 --spring-wall=0 --spring-eq=30.5
   --eff-output=spring --max-right-piston-travel=7.96)
D=("${GAS2[@]}" --wall-mass-factors=200,1000 --spring-k-sigma=0.494 --spring-wall=0 --spring-eq=30.5
   --eff-output=spring --max-right-piston-travel=3.93)

for g in A B C D; do
  eval "flags=(\"\${$g[@]}\")"
  for mode in paper experiment; do
    bmp=$OUT/260920_master_geom${g}_${mode}.bmp
    png=${bmp%.bmp}.png
    rm -f "$bmp" "$png"
    # --demo-shot fires the moment the piston is actually stepping, so every picture is mid-push.
    "$BIN" "${COMMON[@]}" "${flags[@]}" --render=$mode --demo-shot="$bmp",2600 \
        > "$OUT/log_${g}_${mode}.txt" 2>&1 &
    pid=$!
    for _ in $(seq 40); do [ -s "$bmp" ] && break; sleep 0.5; done
    sleep 1; kill "$pid" 2>/dev/null; wait "$pid" 2>/dev/null
    if [ -s "$bmp" ]; then
      sips -s format png "$bmp" --out "$png" > /dev/null 2>&1 && rm -f "$bmp"
      printf "  %s %-11s %s\n" "$g" "$mode" "$png"
    else
      printf "  %s %-11s MISSING -- see %s\n" "$g" "$mode" "$OUT/log_${g}_${mode}.txt"
    fi
  done
done
echo "eight pictures in $OUT"
