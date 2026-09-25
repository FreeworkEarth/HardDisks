#!/usr/bin/env bash
# ##CHRIS 2026-10-10: WATCH any Level-4 experiment live, and/or grab the two pictures.
#
#   ./watch.sh demo            open the GUI and watch the push experiment (B1long)
#   ./watch.sh demo shot       don't watch -- just save the experiment + paper screenshots
#   ./watch.sh equil           the equilibrium (no-piston) run that tau_T comes from
#   ./watch.sh effmap [k]      geometry C, the efficiency map (k = 0.25 | 0.5 | 1.0, default 0.5)
#   ./watch.sh 4b              the Level 4b transmission cell (M = 200, u = 0.2)
#
# Keys in the GUI: P start piston, S screenshot, +/- speed, Q quit.
#
# THE PICTURE RULE, from 2026-10-10: EVERY campaign gets both pictures before it launches.
#   ./watch.sh <experiment> shot
# runs first, the two PNGs go in paper2_pictures/, and the paper-render one is copied into
# 0000_PLAN_OVERALL/paper2_energytransfer/experiments/final/ with the campaign's date prefix.
# A campaign whose geometry nobody has looked at is a campaign nobody can check: the 0.25-sigma
# piston gap was visible in the very first paper render (piston at 78.44 with W_in = 0.000) and
# cost three rounds of inference to find in the numbers.
#
# effmap flags follow 261010 section 1: the spring rest length is NOT free -- it is set per k so
# the system starts in mechanical equilibrium, k (x_eq - 30.5) = P h = 1.57489:
#     k = 0.25 -> x_eq = 36.7996     k = 0.50 -> x_eq = 33.6498     k = 1.00 -> x_eq = 32.0749
# Level 3's 33.65 is the k = 0.5 value; using it at k = 0.25 or 1.0 starts the run with a 50-100 %
# force preload and the measured efficiency is meaningless.
set -uo pipefail
cd /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3
WHAT=${1:-demo}; MODE=${2:-watch}
OUT=experiments_energy_transfer/paper2_pictures; mkdir -p "$OUT"
STAMP=$(date +%y%m%d)

COMMON=(--mode=edmd --experiment=energy_transfer --edmd-acc=0 --seed-drift-order=drift-first
        --particle-radius=0.5 --height=10 --fixed-dt=0.4 --kbt1 --seed=9600)

case "$WHAT" in
  demo)   FLAGS=(--particles=100 --particles-boxes=50,50 --l0=39.25 --num-walls=1
                 --wall-positions=39.25 --wall-mass-factors=200 --eff-output=wall-ke
                 --piston-right-protocol-mode=step --velocity-right-piston-step=0.2
                 --max-right-piston-travel=3.875 --auto-piston-step
                 --wall-hold-steps=12000 --steps=400000); SHOT_AT=12700 ;;
  4b)     FLAGS=(--particles=100 --particles-boxes=50,50 --l0=39.25 --num-walls=1
                 --wall-positions=39.25 --wall-mass-factors=200 --eff-output=wall-ke
                 --piston-right-protocol-mode=step --velocity-right-piston-step=0.2
                 --max-right-piston-travel=3.875 --auto-piston-step
                 --wall-hold-steps=12000 --steps=630000); SHOT_AT=12700 ;;
  equil)  FLAGS=(--particles=200 --particles-boxes=100,100 --l0=78.0 --num-walls=1
                 --wall-positions=78.0 --wall-mass-factors=100 --eff-output=wall-ke
                 --wall-hold-steps=12000 --steps=400000); SHOT_AT=20000 ;;
  effmap) K=${3:-0.5}
          case "$K" in 0.25) XEQ=36.7996 ;; 0.5) XEQ=33.6498 ;; 1.0) XEQ=32.0749 ;;
            *) echo "k must be 0.25, 0.5 or 1.0 (the pre-registered grid)"; exit 1 ;; esac
          echo "  effmap: k = $K, spring rest length x_eq = $XEQ (mechanical equilibrium at t = 0)"
          FLAGS=(--particles=100 --particles-boxes=0,100 --l0=54.75 --num-walls=1
                 --wall-positions=30.5 --wall-mass-factors=200 --spring-k-sigma=$K
                 --spring-wall=0 --spring-eq=$XEQ --eff-output=spring
                 --piston-right-protocol-mode=step --velocity-right-piston-step=0.05
                 --max-right-piston-travel=7.96 --auto-piston-step
                 --wall-hold-steps=12000 --steps=200000); SHOT_AT=14000 ;;
  *) echo "unknown experiment '$WHAT'; try: demo 4b equil effmap"; exit 1 ;;
esac

if [ "$MODE" = "watch" ]; then
  echo "opening the GUI for '$WHAT'. Keys: P piston, S screenshot, +/- speed, Q quit."
  exec ./00ALLINONE "${COMMON[@]}" "${FLAGS[@]}" --show-simulation --demo --render=experiment
fi

for m in experiment paper; do
  bmp="$OUT/${STAMP}_${WHAT}_${m}.bmp"; png="${bmp%.bmp}.png"; rm -f "$bmp" "$png"
  ./00ALLINONE "${COMMON[@]}" "${FLAGS[@]}" --show-simulation --demo --render=$m \
      --demo-shot="$bmp",$SHOT_AT > "$OUT/log_${WHAT}_${m}.txt" 2>&1 &
  pid=$!
  for _ in $(seq 60); do [ -s "$bmp" ] && break; sleep 0.5; done
  sleep 1; kill "$pid" 2>/dev/null; wait "$pid" 2>/dev/null
  if [ -s "$bmp" ]; then
    # sips cannot read the 32-bit BMP SDL writes ("Error 13"); PIL can. Also auto-crop the
    # large empty margin the full-window grab leaves.
    python3 - "$bmp" "$png" <<'PYX'
import sys
from PIL import Image, ImageChops
im = Image.open(sys.argv[1]).convert("RGB")
bg = Image.new("RGB", im.size, im.getpixel((im.size[0]-1, im.size[1]-1)))
bb = ImageChops.difference(im, bg).convert("L").point(lambda v: 255 if v > 8 else 0).getbbox()
if bb:
    pad = 24
    bb = (max(0,bb[0]-pad), max(0,bb[1]-pad), min(im.size[0],bb[2]+pad), min(im.size[1],bb[3]+pad))
    im = im.crop(bb)
im.save(sys.argv[2], "PNG", optimize=True)
print("   cropped to", im.size)
PYX
    rm -f "$bmp"; echo "  $m -> $png"
  else echo "  $m MISSING (no display? see $OUT/log_${WHAT}_${m}.txt)"; fi
done
