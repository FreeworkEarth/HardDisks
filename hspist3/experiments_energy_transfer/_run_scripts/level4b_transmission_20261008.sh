#!/usr/bin/env bash
# ##CHRIS 2026-10-08: LEVEL 4b -- acoustic transmission scan. Pre-registered in
# 261008_paper2_level4b.md section 1, committed BEFORE this script ran. Do not edit to change a
# prediction; if the verdict is uncomfortable, the verdict is the result.
#
# QUESTION: of the work the piston does on gas 1, what fraction crosses the divider on the FIRST
# PASS, and is the divider a mass on a string as far as that pass is concerned?
#   |t(w)|^2 = 1/(1 + (w M_d/2Z)^2)   [SOURCE: Morse & Ingard, Theoretical Acoustics, section 4.3]
#   Z = N m c_s/L_c = 50 x 1.749302/38.75 = 2.2572   [identification with the gas impedance: INFERENCE]
#   f_early = |t|^2/2  -- the half is required by the plan's own x<<1 limit of 0.5, and is supplied.
#   x = M_d/(2 Z tau_push) spans 0.029 to 11.4 across the grid.
#
# THREE THINGS THE PLAN AS WRITTEN DID NOT CLOSE, all resolved in section 1.1 of the report:
#   d = 3.875, NOT 3.876 -- 3.876 x 24 = 93.024 fails the pixel grid; 3.875 x 24 = 93 is exact and
#     reproduces B1long's 19.375-sigma push at u = 0.2 exactly.
#   "N = 100 per gas" would give eta = 0.2027, twice the stated 0.10134. N = 100 TOTAL, 50 per side,
#     is the only reading consistent with eta_phys = 0.10134170, with the plan's own Z ~ 2.3, and
#     with reusing B2pilot.
#   f_early = |t|^2/2, see above.
#
# GEOMETRY: l0 = 39.25, L_c = 38.75, box 78.5, t = 1.0 (recorded), eta = 0.10134170, N_s = 50.
# RUN LENGTH: tau_push + L_c/c_s + 5 tau_r(M_d), NOT tau_T. First-pass transmission is the target
#   and tau_T(200) ~ 40 000 is deliberately out of reach. tau_r = 208/586/2090 for M = 10/50/200,
#   from 261006_mode_ladder.json.
# CADENCE: --trace-every=12 (dt = 0.2 sigma) everywhere. That is <= tau_r/20 by a factor 52 at the
#   tightest, AND gives 19 samples across the shortest push (3.875 sigma at u = 1.0) -- which a
#   cadence chosen from tau_r/20 alone would not.
# THE (200, 1.0) CELL IS NOT HERE: it reuses B2pilot (same box, same M_d, same u, same d, 8 seeds).
#   B2pilot's dt = 10 cannot resolve the push, but f_early does not need it -- only W_in from the
#   summary ledger and Delta E2 on the plateau, which carries 627 samples at dt = 10.
#
# COST: 11 cells x 8 seeds = 28.0 M steps = 0.12 core-hours. Essentially free. If the seed errors
# make the chi^2 verdict uninformative, the honest move is to rerun at 32 or 80 seeds (0.5 / 1.2
# core-hours) and say the first pass was under-powered -- NOT to reinterpret the threshold.
set -uo pipefail
cd /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3
P=experiments_energy_transfer/level4b_transmission_20261008
MAXJOBS=${MAXJOBS:-9}

run_one() {
  local M=$1 U=$2 STEPS=$3 tag=$4 sd=$5
  local d=$P/$tag
  [ -s "$d/red_${sd}.csv" ] && return 0
  nice -n 5 ./00ALLINONE --mode=edmd --experiment=energy_transfer --headless --quiet \
      --edmd-acc=0 --seed-drift-order=drift-first \
      --energy-transfer-summary="$d/summary_${sd}.csv" \
      --energy-transfer-trace="$d/raw_${sd}.csv" --trace-every=12 \
      --particles=100 --particles-boxes=50,50 --particle-radius=0.5 \
      --l0=39.25 --height=10 --num-walls=1 --wall-positions=39.25 \
      --wall-mass-factors=$M --eff-output=wall-ke \
      --piston-right-protocol-mode=step --velocity-right-piston-step=$U \
      --max-right-piston-travel=3.875 --auto-piston-step \
      --wall-hold-steps=12000 --steps=$STEPS --fixed-dt=0.4 --kbt1 --seed=$sd \
      > "$d/run_${sd}.log" 2>&1
  local rc=$? health
  health=$(grep -cE 'EDMD-HEALTH|forced_advance|clamp_repair|overlap_repair|wall_overdue' "$d/run_${sd}.log" 2>/dev/null) || true
  health=${health:-0}
  if [ "$rc" -ne 0 ] || [ "$health" -ne 0 ]; then
    echo "$tag seed=$sd FAILED rc=$rc health=$health"; return 1
  fi
  python3 - "$d/raw_${sd}.csv" "$d/red_${sd}.csv" <<'PY'
import sys, pandas as pd
want = ["Time","KE_gas_left","KE_gas_right","W0_x_sigma","W0_v_sigma","PistonWork","W_in"]
e = pd.read_csv(sys.argv[1], low_memory=False)
e[[c for c in want if c in e.columns]].to_csv(sys.argv[2], index=False)
PY
  if [ -s "$d/red_${sd}.csv" ]; then rm -f "$d/raw_${sd}.csv"; else
    echo "$tag seed=$sd reduction FAILED"; return 1; fi
  return 0
}

# "M u steps tag" -- longest first. (200, 1.0) omitted on purpose: B2pilot supplies it.
# (200, 1.0) ADDED POST-HOC, 2026-10-08, after the first 11 cells had run and before any 4b trace
# was opened. It was originally omitted on the grounds that B2pilot supplies it; B2pilot's
# --trace-every=600 (dt = 10) cannot resolve the 3.875-sigma push or the 44.3-sigma return spacing,
# which f_first and the mechanism test both need. 32 seeds (not 8) because it also carries the
# 3-sigma "the divider shares dissipated work" claim. Seeds 9950-9981, distinct from the 9900 block.
SPECS=(
  "200 1.0 629000 M200_u100"
  "200 0.05 633000 M200_u005" "200 0.2 630000 M200_u020" "200 0.5 629000 M200_u050"
  "50 0.05 182000 M50_u005"  "50 0.2 179000 M50_u020"  "50 0.5 178000 M50_u050"  "50 1.0 178000 M50_u100"
  "10 0.05 69000 M10_u005"   "10 0.2 65000 M10_u020"   "10 0.5 65000 M10_u050"   "10 1.0 64000 M10_u100"
)
for s in "${SPECS[@]}"; do set -- $s; mkdir -p $P/$4; done
echo "level4b start $(date +%H:%M:%S), max $MAXJOBS concurrent, 11 cells x 8 seeds"
for s in "${SPECS[@]}"; do
  set -- $s; M=$1; U=$2; ST=$3; TAG=$4
  NSEED=8; BASE=9900
  [ "$TAG" = "M200_u100" ] && { NSEED=32; BASE=9950; }
  for i in $(seq 0 $((NSEED-1))); do
    while [ "$(jobs -rp | wc -l | tr -d ' ')" -ge "$MAXJOBS" ]; do sleep 5; done
    run_one "$M" "$U" "$ST" "$TAG" $((BASE+i)) &
  done
done
wait
echo "level4b done $(date +%H:%M:%S)"
for s in "${SPECS[@]}"; do
  set -- $s; d=$P/$4
  echo "  $4: $(find $d -name 'red_*.csv'|wc -l|tr -d ' '); health $(cat $d/run_*.log 2>/dev/null | grep -cE 'EDMD-HEALTH|forced_advance|clamp_repair|overlap_repair|wall_overdue'); aborts $(cat $d/run_*.log 2>/dev/null | grep -c ABORTING)"
done
du -sh $P; df -h . | tail -1
