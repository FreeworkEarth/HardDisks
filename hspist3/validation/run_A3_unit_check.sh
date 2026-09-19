#!/usr/bin/env bash
# ##CHRIS 2026-09-12: A3 -- reduced-unit consistency check at eta = 0.30.
#
# The box is held FIXED in the code's absolute units at the A2 N=100 geometry
# (L0 = 13.089969, H = 10) and only the particle count and radius change:
#   N = 100, r = 0.5        -> box is 13.09 x 10 particle diameters   (= A2 N=100)
#   N = 400, r = 0.25       -> box is 26.18 x 20 particle diameters   (= A2 N=400)
#   N = 900, r = 1/6        -> box is 39.27 x 30 particle diameters   (= A2 N=900)
# eta = N pi r^2 / (2 L0 H) = 0.30 in all three; verified from the trace header.
#
# Each row is therefore the SAME physical system as the A2 cell with the same N,
# just written with a different particle size. c_s is measured in units of
# sqrt(kT/m), which carries no length unit, so c_s and the per-mass c_s must agree
# with A2 within the fit error. The frequencies must differ by exactly the ratio of
# L_eff in code units (nu_A3 / nu_A2 = 2 at N=400, 3 at N=900). Any departure beyond
# the fit error is a unit bug, not physics.
set -uo pipefail
D="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$D/.."
OUT="${1:?usage: run_A3_unit_check.sh OUTDIR [JOBS]}"
JOBS="${2:-6}"
MASSES=50,200,500,1000,2000
REPEATS=10
BASE=21660913   # disjoint from famB_20260911
L0=13.089969
H=10

if otool -L ./00ALLINONE 2>/dev/null | grep -qi asan; then
  echo "REFUSING: AddressSanitizer build." >&2; exit 1
fi
mkdir -p "$OUT"

i=0
for spec in "100:0.5" "400:0.25" "900:0.16666666666666666"; do
  N="${spec%%:*}"; r="${spec##*:}"
  dir="$OUT/N${N}_r${r}"
  if [ -f "$dir/.done" ]; then i=$((i+1)); continue; fi
  mkdir -p "$dir"
  (
    HD_KE_TRACE=1 ./00ALLINONE --mode=edmd --experiment=speed_of_sound --headless --kbt1 \
      --seed-drift-order=drift-first --edmd-acc=0 \
      --particles=$N --particles-boxes=$((N/2)),$((N/2)) \
      --height=$H --particle-radius=$r \
      --wall-thickness=0.05 --wall-thickness-vis=0.05 \
      --lengths=$L0 --wall-masses=$MASSES --repeats=$REPEATS --seed=$(( BASE + i )) \
      --wall-hold-steps=2000 --fixed-dt=0.4 \
      --target-oscillations=25 --oscillation-safety=1.5 \
      --oscillation-min-steps=10000 --oscillation-max-steps=40000000 \
      --speed-sound-log-stride=auto --speed-sound-run-dir="$dir" > "$dir/run.log" 2>&1
    rc=$?
    inv=$(grep -c "INVALID" "$dir/run.log" 2>/dev/null); inv=${inv:-0}
    hl=$(grep -c "EDMD-HEALTH" "$dir/run.log" 2>/dev/null); hl=${hl:-0}
    tr=$(ls "$dir"/wall_x_positions_*.csv 2>/dev/null | wc -l | tr -d ' ')
    [ "$rc" -eq 0 ] && touch "$dir/.done"
    printf "  done N=%-5s r=%-20s rc=%s traj=%-4s invalid=%-3s health=%-3s [%s]\n" \
           "$N" "$r" "$rc" "$tr" "$inv" "$hl" "$(date +%H:%M)"
  ) &
  i=$((i+1))
  while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do sleep 3; done
done
wait
echo
echo "cells: $(find "$OUT" -name .done | wc -l | tr -d ' ')/3"
echo "trajectories: $(find "$OUT" -name 'wall_x_positions_*.csv' | wc -l | tr -d ' ')/150"
echo "invalid=$(grep -rh INVALID "$OUT" --include=run.log 2>/dev/null | wc -l | tr -d ' ')  health=$(grep -rh EDMD-HEALTH "$OUT" --include=run.log 2>/dev/null | wc -l | tr -d ' ')"
echo "Complete: $OUT"
