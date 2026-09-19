#!/usr/bin/env bash
# ##CHRIS: famB -- finite-size c_s at the famA aspect ratio, drift-first seeding.
#   eta = 0.10 0.30 0.50 0.60 0.65 ; N = 100 400 900 1600 ; masses 50 200 500 1000 2000 ;
#   10 repeats ; r = 0.5 fixed ; L0 and H both ~ sqrt(N), so the aspect ratio is fixed
#   at 2*L0/H = 0.7854/eta for each eta, exactly as famA.
#
# T_i is recorded per trajectory via HD_KE_TRACE=1: the audit line printed just
# before each "Running:" line carries KE_left/KE_right after seeding. With
# --seed-drift-order=drift-first it should read exactly 50/50 per compartment,
# i.e. kT = 1.
#
# SEED RULE -- documented deviation. The campaign rule is
# speed_sound_run_seed(base, l, m, r). Running all 5 masses in one invocation would
# reproduce it exactly, but the longest such job is ~40 h and cannot be split. The
# unit of work here is (eta, N, mass) with a distinct base seed per cell, so the
# per-run seeds are deterministic and reproducible, but they are NOT the integers a
# single-invocation run would produce.
set -uo pipefail
D="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$D/.."
OUT="${1:?usage: run_famB_campaign.sh OUTDIR [JOBS]}"
JOBS="${2:-10}"
ETAS=(0.10 0.30 0.50 0.60 0.65)
NS=(100 400 900 1600)
MASSES=(50 200 500 1000 2000)
REPEATS=10
BASE=20260911

if otool -L ./00ALLINONE 2>/dev/null | grep -qi asan; then
  echo "REFUSING: AddressSanitizer build." >&2; exit 1
fi
AV=$(df -m "$(dirname "$OUT")" 2>/dev/null | awk 'NR==2{print $4}')
if [ -n "$AV" ] && [ "$AV" -lt 20000 ]; then
  echo "REFUSING: only ${AV} MB free." >&2; exit 1
fi
mkdir -p "$OUT"

ei=0
for eta in "${ETAS[@]}"; do
  ni=0
  for N in "${NS[@]}"; do
    f=$(/opt/homebrew/bin/python3 -c "import math;print(math.sqrt($N/100))")
    L0=$(/opt/homebrew/bin/python3 -c "print(f'{3.926990816987241/$eta*$f:.6f}')")
    H=$(/opt/homebrew/bin/python3 -c "print(f'{10*$f:.6f}')")
    mi=0
    for M in "${MASSES[@]}"; do
      dir="$OUT/eta_$(echo "$eta" | tr '.' 'p')/N${N}/m_${M}"
      if [ -f "$dir/.done" ]; then mi=$((mi+1)); continue; fi
      mkdir -p "$dir"
      seed=$(( BASE + 100000*ei + 1000*ni + mi ))
      (
        HD_KE_TRACE=1 ./00ALLINONE --mode=edmd --experiment=speed_of_sound --headless --kbt1 \
          --seed-drift-order=drift-first --edmd-acc=0 \
          --particles=$N --particles-boxes=$((N/2)),$((N/2)) \
          --height="$H" --particle-radius=0.5 \
          --wall-thickness=0.05 --wall-thickness-vis=0.05 \
          --lengths="$L0" --wall-masses=$M --repeats=$REPEATS --seed=$seed \
          --wall-hold-steps=2000 --fixed-dt=0.4 \
          --target-oscillations=25 --oscillation-safety=1.5 \
          --oscillation-min-steps=10000 --oscillation-max-steps=10000000 \
          --speed-sound-log-stride=auto --speed-sound-run-dir="$dir" > "$dir/run.log" 2>&1
        inv=$(grep -c "INVALID" "$dir/run.log" 2>/dev/null); inv=${inv:-0}
        hl=$(grep -c "EDMD-HEALTH" "$dir/run.log" 2>/dev/null); hl=${hl:-0}
        touch "$dir/.done"
        printf "  done eta=%-5s N=%-5s M=%-5s invalid=%-3s health=%-3s [%s]\n" \
               "$eta" "$N" "$M" "$inv" "$hl" "$(date +%H:%M)"
      ) &
      mi=$((mi+1))
      while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do sleep 3; done
    done
    ni=$((ni+1))
  done
  ei=$((ei+1))
done
wait

echo
echo "cells: $(find "$OUT" -name .done | wc -l | tr -d ' ')/100"
echo "trajectories: $(find "$OUT" -name 'wall_x_positions_*.csv' | wc -l | tr -d ' ')/1000"
echo "invalid=$(grep -rh INVALID "$OUT" --include=run.log 2>/dev/null | wc -l | tr -d ' ')  health=$(grep -rh EDMD-HEALTH "$OUT" --include=run.log 2>/dev/null | wc -l | tr -d ' ')"
echo "Complete: $OUT"
