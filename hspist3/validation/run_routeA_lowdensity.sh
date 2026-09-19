#!/usr/bin/env bash
# ##CHRIS 2026-09-12: route-A extension below eta = 0.02.
#   eta = 0.013090 / 0.009817 / 0.006545  (L0 = 300 / 400 / 600, exact on the pixel grid:
#   2*L0*24 = 14400 / 19200 / 28800 px), N = 100 (50/50), r = 0.5, H = 10,
#   the 9 route-A masses 50..2000, 10 repeats, drift-first seeding, 25 oscillations.
#
# SEED RULE -- reproduces the route-A convention exactly: one invocation per eta with all
# 9 masses, so speed_sound_run_seed(base, l=0, m, r) is the same function route-A used.
# Each eta gets its own base seed; within an eta the (m, r) indices are the campaign's.
#
# STEP BUDGET. nu scales as 1/L_eff, so the route-A cap of 1e7 steps is too small here.
# Measured from pilot traces (Predicted_Frequency at M = 2000, dt = 0.4/24 sigma-time):
#   L0 = 300  nu = 1.71384e-4  -> 1.31e7 steps needed
#   L0 = 400  nu = 1.27586e-4  -> 1.76e7
#   L0 = 600  nu = 8.44287e-5  -> 2.67e7
# The cap below is 4e7, above the largest requirement for the heaviest mass.
#
# HEALTH. Strict contract: a trajectory is accepted only with zero forced_advance,
# zero clamp_repair, zero overlap_repair and zero wall_overdue. This script only
# records; the acceptance test is applied in analyze_routeA_lowdensity.py.
set -uo pipefail
D="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$D/.."
OUT="${1:?usage: run_routeA_lowdensity.sh OUTDIR}"
BIN="${2:-./00ALLINONE}"
L0S=(300 400 600)
MASSES=50,100,200,300,500,750,1000,1500,2000
REPEATS=10
BASE=2026091200

if otool -L "$BIN" 2>/dev/null | grep -qi asan; then
  echo "REFUSING: AddressSanitizer build ($BIN)." >&2; exit 1
fi
AV=$(df -m "$(dirname "$OUT")" 2>/dev/null | awk 'NR==2{print $4}')
if [ -n "$AV" ] && [ "$AV" -lt 20000 ]; then
  echo "REFUSING: only ${AV} MB free." >&2; exit 1
fi
mkdir -p "$OUT"

li=0
for L0 in "${L0S[@]}"; do
  eta=$(/opt/homebrew/bin/python3 -c "print(f'{3.926990816987241/$L0:.6f}')")
  dir="$OUT/eta_$(echo "$eta" | tr '.' 'p')"
  if [ -f "$dir/.done" ]; then li=$((li+1)); continue; fi
  mkdir -p "$dir"
  seed=$(( BASE + li ))
  (
    HD_KE_TRACE=1 "$BIN" --mode=edmd --experiment=speed_of_sound --headless --quiet --kbt1 \
      --seed-drift-order=drift-first --edmd-acc=0 \
      --particles=100 --particles-boxes=50,50 \
      --height=10.0 --particle-radius=0.5 \
      --wall-thickness=0.05 --wall-thickness-vis=0.05 \
      --lengths="$L0" --wall-masses=$MASSES --repeats=$REPEATS --seed=$seed \
      --wall-hold-steps=2000 --fixed-dt=0.4 \
      --target-oscillations=25 --oscillation-safety=1.5 \
      --oscillation-min-steps=10000 --oscillation-max-steps=40000000 \
      --speed-sound-log-stride=auto --speed-sound-run-dir="$dir" > "$dir/run.log" 2>&1
    rc=$?
    inv=$(grep -c "INVALID" "$dir/run.log" 2>/dev/null); inv=${inv:-0}
    hl=$(grep -c "EDMD-HEALTH" "$dir/run.log" 2>/dev/null); hl=${hl:-0}
    tr=$(ls "$dir"/wall_x_positions_*.csv 2>/dev/null | wc -l | tr -d ' ')
    [ "$rc" -eq 0 ] && touch "$dir/.done"
    printf "  done L0=%-4s eta=%-9s rc=%s traj=%-4s invalid=%-3s health=%-3s [%s]\n" \
           "$L0" "$eta" "$rc" "$tr" "$inv" "$hl" "$(date +%H:%M)"
  ) &
  li=$((li+1))
done
wait

echo
echo "leaves: $(find "$OUT" -name .done | wc -l | tr -d ' ')/3"
echo "trajectories: $(find "$OUT" -name 'wall_x_positions_*.csv' | wc -l | tr -d ' ')/270"
echo "invalid=$(grep -rh INVALID "$OUT" --include=run.log 2>/dev/null | wc -l | tr -d ' ')  health=$(grep -rh EDMD-HEALTH "$OUT" --include=run.log 2>/dev/null | wc -l | tr -d ' ')"
echo "Complete: $OUT"
