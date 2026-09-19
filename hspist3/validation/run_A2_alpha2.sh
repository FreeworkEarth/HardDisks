#!/usr/bin/env bash
# ##CHRIS 2026-09-12: A2 alpha = 2 cells -- close the alpha confound.
# The Roman piston parameter is alpha = M/(2 N_side) = M/N. Roman used N = 100, M = 200,
# i.e. alpha = 2. Our A2 ladder tops out at M = 2000, which is alpha = 1.25 at N = 1600
# and 0.80 at N = 2500, so the large-N cells sample only lighter pistons than Roman and
# c_s-vs-N carries a co-varying change of piston regime. These four cells restore
# alpha = 2 at the two largest sizes:
#     N = 1600 -> M = 3200      N = 2500 -> M = 5000
# eta = 0.10 and 0.30, 25 repeats. Base seed disjoint from famB and from the top-up.
set -uo pipefail
D="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$D/.."
OUT="${1:?usage: run_A2_alpha2.sh OUTDIR [JOBS] [REPEATS]}"
JOBS="${2:-4}"
REPEATS="${3:-25}"
BASE=21860912

if otool -L ./00ALLINONE 2>/dev/null | grep -qi asan; then
  echo "REFUSING: AddressSanitizer build." >&2; exit 1
fi
mkdir -p "$OUT"

i=0
for eta in 0.10 0.30; do
  for pair in "1600:3200" "2500:5000"; do
    N="${pair%%:*}"; M="${pair##*:}"
    fac=$(/opt/homebrew/bin/python3 -c "import math;print(math.sqrt($N/100))")
    L0=$(/opt/homebrew/bin/python3 -c "print(f'{3.926990816987241/$eta*$fac:.6f}')")
    H=$(/opt/homebrew/bin/python3 -c "print(f'{10*$fac:.6f}')")
    dir="$OUT/eta_$(echo "$eta" | tr '.' 'p')/N${N}/m_${M}"
    if [ -f "$dir/.done" ]; then i=$((i+1)); continue; fi
    mkdir -p "$dir"
    (
      HD_KE_TRACE=1 ./00ALLINONE --mode=edmd --experiment=speed_of_sound --headless --kbt1 \
        --seed-drift-order=drift-first --edmd-acc=0 \
        --particles=$N --particles-boxes=$((N/2)),$((N/2)) \
        --height="$H" --particle-radius=0.5 \
        --wall-thickness=0.05 --wall-thickness-vis=0.05 \
        --lengths="$L0" --wall-masses=$M --repeats=$REPEATS --seed=$(( BASE + i )) \
        --wall-hold-steps=2000 --fixed-dt=0.4 \
        --target-oscillations=25 --oscillation-safety=1.5 \
        --oscillation-min-steps=10000 --oscillation-max-steps=80000000 \
        --speed-sound-log-stride=auto --speed-sound-run-dir="$dir" > "$dir/run.log" 2>&1
      rc=$?
      hl=$(grep -c "EDMD-HEALTH" "$dir/run.log" 2>/dev/null); hl=${hl:-0}
      tr=$(ls "$dir"/wall_x_positions_*.csv 2>/dev/null | wc -l | tr -d ' ')
      [ "$rc" -eq 0 ] && touch "$dir/.done"
      printf "  done eta=%-5s N=%-5s M=%-5s alpha=%.2f rc=%s traj=%-4s health=%-3s [%s]\n" \
             "$eta" "$N" "$M" "$(echo "$M/$N" | bc -l)" "$rc" "$tr" "$hl" "$(date +%H:%M)"
    ) &
    i=$((i+1))
    while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do sleep 5; done
  done
done
wait
echo
echo "cells: $(find "$OUT" -name .done | wc -l | tr -d ' ')/4   traj: $(find "$OUT" -name 'wall_x_positions_*.csv' | wc -l | tr -d ' ')"
echo "health=$(grep -rh EDMD-HEALTH "$OUT" --include=run.log 2>/dev/null | wc -l | tr -d ' ')"
echo "Complete: $OUT"
