#!/usr/bin/env bash
# ##CHRIS 2026-09-12: A2 top-up (fixed-eta N ladder).
#   Part A: eta = 0.10, 0.30 at N = 900, 1600 -- 25 EXTRA repeats per (eta,N,mass) cell.
#   Part B: N = 2500 at eta = 0.10, 0.30 -- 10 repeats, the same 5 masses.
#
# Written to a SEPARATE tree from famB_20260911. The runner names traces
# wall_x_positions_<L0>_wallmassfactor_<M>_run<r>.csv with r starting at 0, so writing
# a top-up into the original leaf would silently overwrite run0..run9. The analysis
# merges the two trees instead; nothing under famB_20260911 is touched.
#
# Seed rule as in famB: one invocation per (eta, N, mass) cell with its own base seed,
# so the per-run seeds are deterministic. BASE is new (20260912) so no top-up run can
# repeat a seed already on disk.
set -uo pipefail
D="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$D/.."
OUT="${1:?usage: run_A2_topup.sh OUTDIR [JOBS] [part]}"
JOBS="${2:-10}"
PART="${3:-all}"          # a | b | all
MASSES=(50 200 500 1000 2000)
BASE=21260912   # disjoint from famB_20260911 (see seed-collision note)

if otool -L ./00ALLINONE 2>/dev/null | grep -qi asan; then
  echo "REFUSING: AddressSanitizer build." >&2; exit 1
fi
AV=$(df -m "$(dirname "$OUT")" 2>/dev/null | awk 'NR==2{print $4}')
if [ -n "$AV" ] && [ "$AV" -lt 30000 ]; then
  echo "REFUSING: only ${AV} MB free." >&2; exit 1
fi
mkdir -p "$OUT"

cell () {   # eta N M repeats seed
  local eta=$1 N=$2 M=$3 rep=$4 seed=$5
  local f L0 H dir
  f=$(/opt/homebrew/bin/python3 -c "import math;print(math.sqrt($N/100))")
  L0=$(/opt/homebrew/bin/python3 -c "print(f'{3.926990816987241/$eta*$f:.6f}')")
  H=$(/opt/homebrew/bin/python3 -c "print(f'{10*$f:.6f}')")
  dir="$OUT/eta_$(echo "$eta" | tr '.' 'p')/N${N}/m_${M}"
  [ -f "$dir/.done" ] && return 0
  mkdir -p "$dir"
  (
    HD_KE_TRACE=1 ./00ALLINONE --mode=edmd --experiment=speed_of_sound --headless --kbt1 \
      --seed-drift-order=drift-first --edmd-acc=0 \
      --particles=$N --particles-boxes=$((N/2)),$((N/2)) \
      --height="$H" --particle-radius=0.5 \
      --wall-thickness=0.05 --wall-thickness-vis=0.05 \
      --lengths="$L0" --wall-masses=$M --repeats=$rep --seed=$seed \
      --wall-hold-steps=2000 --fixed-dt=0.4 \
      --target-oscillations=25 --oscillation-safety=1.5 \
      --oscillation-min-steps=10000 --oscillation-max-steps=40000000 \
      --speed-sound-log-stride=auto --speed-sound-run-dir="$dir" > "$dir/run.log" 2>&1
    rc=$?
    inv=$(grep -c "INVALID" "$dir/run.log" 2>/dev/null); inv=${inv:-0}
    hl=$(grep -c "EDMD-HEALTH" "$dir/run.log" 2>/dev/null); hl=${hl:-0}
    tr=$(ls "$dir"/wall_x_positions_*.csv 2>/dev/null | wc -l | tr -d ' ')
    [ "$rc" -eq 0 ] && touch "$dir/.done"
    printf "  done eta=%-5s N=%-5s M=%-5s rc=%s traj=%-4s invalid=%-3s health=%-3s [%s]\n" \
           "$eta" "$N" "$M" "$rc" "$tr" "$inv" "$hl" "$(date +%H:%M)"
  ) &
  while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do sleep 5; done
}

ei=0
for eta in 0.10 0.30; do
  ni=0
  for N in 900 1600 2500; do
    mi=0
    for M in "${MASSES[@]}"; do
      seed=$(( BASE + 100000*ei + 1000*ni + mi ))
      if [ "$N" -eq 2500 ]; then
        [ "$PART" = a ] || cell "$eta" "$N" "$M" 10 "$seed"
      else
        [ "$PART" = b ] || cell "$eta" "$N" "$M" 25 "$seed"
      fi
      mi=$((mi+1))
    done
    ni=$((ni+1))
  done
  ei=$((ei+1))
done
wait

echo
echo "cells: $(find "$OUT" -name .done | wc -l | tr -d ' ')"
echo "trajectories: $(find "$OUT" -name 'wall_x_positions_*.csv' | wc -l | tr -d ' ')"
echo "invalid=$(grep -rh INVALID "$OUT" --include=run.log 2>/dev/null | wc -l | tr -d ' ')  health=$(grep -rh EDMD-HEALTH "$OUT" --include=run.log 2>/dev/null | wc -l | tr -d ' ')"
echo "Complete: $OUT"
