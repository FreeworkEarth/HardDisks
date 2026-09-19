#!/usr/bin/env bash
# ##CHRIS 2026-09-12: A4 -- particle-number ladder in a FIXED box.
#   L0 = 20, H = 10, r = 0.5, one held-then-released divider at the middle.
#   N per side = 1, 2, 5, 10, 25, 50, 100, 165  ->  eta = 0.003927 ... 0.647953
#   masses 50 200 500 1000 2000, 10 repeats.
#
# N_side = 1 IS NOT RUNNABLE WITH DRIFT-FIRST SEEDING. remove_drift_segment subtracts
# the compartment centre-of-mass velocity; with one particle per compartment that IS
# its own velocity, so the particle is left exactly at rest. Measured: KE_tot = 0,
# kT_mean = 0. There is no gas and no sound. N_side = 1 is therefore run with the
# DEFAULT seed order (--seed-drift-order=old) and flagged; N_side = 2 is run BOTH ways so
# the size of the seeding difference is measurable rather than assumed.
# More generally drift removal costs each compartment 2 of its 2*N_side momentum
# degrees of freedom, a 1/N_side effect -- precisely the regime this ladder probes.
set -uo pipefail
D="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$D/.."
OUT="${1:?usage: run_A4_particle_ladder.sh OUTDIR [JOBS]}"
JOBS="${2:-10}"
NSIDES=(1 2 5 10 25 50 100 165)
MASSES=50,200,500,1000,2000
REPEATS=10
BASE=21460912   # disjoint from famB_20260911
L0=20
H=10

if otool -L ./00ALLINONE 2>/dev/null | grep -qi asan; then
  echo "REFUSING: AddressSanitizer build." >&2; exit 1
fi
mkdir -p "$OUT"

launch () {  # nside order tag seed
  local ns=$1 order=$2 tag=$3 seed=$4
  local n=$((ns*2))
  local dir="$OUT/nside_${ns}${tag}"
  [ -f "$dir/.done" ] && return 0
  mkdir -p "$dir"
  (
    HD_KE_TRACE=1 ./00ALLINONE --mode=edmd --experiment=speed_of_sound --headless --kbt1 \
      --seed-drift-order="$order" --edmd-acc=0 \
      --particles=$n --particles-boxes=$ns,$ns \
      --height=$H --particle-radius=0.5 \
      --wall-thickness=0.05 --wall-thickness-vis=0.05 \
      --lengths=$L0 --wall-masses=$MASSES --repeats=$REPEATS --seed=$seed \
      --wall-hold-steps=2000 --fixed-dt=0.4 \
      --target-oscillations=25 --oscillation-safety=1.5 \
      --oscillation-min-steps=10000 --oscillation-max-steps=40000000 \
      --speed-sound-log-stride=auto --speed-sound-run-dir="$dir" > "$dir/run.log" 2>&1
    rc=$?
    inv=$(grep -c "INVALID" "$dir/run.log" 2>/dev/null); inv=${inv:-0}
    hl=$(grep -c "EDMD-HEALTH" "$dir/run.log" 2>/dev/null); hl=${hl:-0}
    tr=$(ls "$dir"/wall_x_positions_*.csv 2>/dev/null | wc -l | tr -d ' ')
    [ "$rc" -eq 0 ] && touch "$dir/.done"
    printf "  done nside=%-4s%-12s rc=%s traj=%-4s invalid=%-3s health=%-3s [%s]\n" \
           "$ns" "$tag" "$rc" "$tr" "$inv" "$hl" "$(date +%H:%M)"
  ) &
  while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do sleep 3; done
}

i=0
for ns in "${NSIDES[@]}"; do
  if [ "$ns" -eq 1 ]; then
    launch "$ns" old _oldseed $(( BASE + 1000 + i ))
  else
    launch "$ns" drift-first "" $(( BASE + i ))
    [ "$ns" -eq 2 ] && launch "$ns" old _oldseed $(( BASE + 1000 + i ))
  fi
  i=$((i+1))
done
wait

echo
echo "cells: $(find "$OUT" -name .done | wc -l | tr -d ' ')/9"
echo "trajectories: $(find "$OUT" -name 'wall_x_positions_*.csv' | wc -l | tr -d ' ')/450"
echo "invalid=$(grep -rh INVALID "$OUT" --include=run.log 2>/dev/null | wc -l | tr -d ' ')  health=$(grep -rh EDMD-HEALTH "$OUT" --include=run.log 2>/dev/null | wc -l | tr -d ' ')"
echo "Complete: $OUT"
