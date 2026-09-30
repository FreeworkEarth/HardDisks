#!/usr/bin/env bash
# ##CHRIS 2026-10-10: Level 4b POWER-UP, per 261008 section 1.4 as pre-registered: the three u = 0.2
# cells at 80 seeds instead of 8, because their excess-ledger denominator came out consistent with
# zero and the six-cell T(x) fit was unconstrained. Same geometry, steps and cadence as
# level4b_transmission_20261008b; --trace-every=12 (dt = 0.2). Binary: v1 + -ffp-contract=off +
# --version (branch RED). The 8-seed first pass was on the contraction-ON binary and is reported
# as "first pass, under-powered"; the two are never mixed in one figure or one fit.
# Seeds 9900-9979 (the first 8 coincide with the first pass's seeds, on a different binary).
set -uo pipefail
cd /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3
P=experiments_energy_transfer/level4b_powerup_20261010
MAXJOBS=${MAXJOBS:-9}
run_one() {
  local M=$1 U=$2 STEPS=$3 TAG=$4 sd=$5
  local d=$P/$TAG
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
  if [ "$rc" -ne 0 ] || [ "$health" -ne 0 ]; then echo "$TAG seed=$sd FAILED rc=$rc health=$health"; return 1; fi
  python3 - "$d/raw_${sd}.csv" "$d/red_${sd}.csv" <<'PY'
import sys, pandas as pd
want = ["Time","KE_gas_left","KE_gas_right","W0_x_sigma","W0_v","PistonWork","PistonR_x_sigma","PistonR_v","SegCounts","SegEtas"]
e = pd.read_csv(sys.argv[1], low_memory=False)
missing = [c for c in want if c not in e.columns]
if missing: sys.exit("REDUCTION ABORTED: columns absent: %s" % ",".join(missing))
e[want].to_csv(sys.argv[2], index=False)
PY
  if [ -s "$d/red_${sd}.csv" ]; then rm -f "$d/raw_${sd}.csv"; else echo "$TAG seed=$sd reduction FAILED"; return 1; fi
  return 0
}
# The pre-registered rerun is the three u = 0.2 cells at 80 seeds. The six-cell fit also uses u = 0.5,
# and the rebaseline rule (05215ea) forbids mixing contraction-on and contraction-off data in one
# fit -- so the WHOLE grid is rerun on this binary: the six fit cells (u = 0.2, 0.5) at 80 seeds,
# the six non-fit cells (u = 0.05, 1.0) at 8. "M u steps tag nseeds"
SPECS=("200 0.2 630000 M200_u020 80" "200 0.5 629000 M200_u050 80" "200 0.05 633000 M200_u005 8" "200 1.0 629000 M200_u100 8"
       "50 0.2 179000 M50_u020 80"   "50 0.5 178000 M50_u050 80"   "50 0.05 182000 M50_u005 8"   "50 1.0 178000 M50_u100 8"
       "10 0.2 65000 M10_u020 80"    "10 0.5 65000 M10_u050 80"    "10 0.05 69000 M10_u005 8"    "10 1.0 64000 M10_u100 8")
for s in "${SPECS[@]}"; do set -- $s; mkdir -p $P/$4; done
echo "level4b power-up start $(date +%H:%M:%S), binary: $(./00ALLINONE --version | head -1)"
for s in "${SPECS[@]}"; do
  set -- $s; M=$1; U=$2; ST=$3; TAG=$4; NS=$5
  for i in $(seq 0 $((NS-1))); do
    while [ "$(jobs -rp | wc -l | tr -d ' ')" -ge "$MAXJOBS" ]; do sleep 5; done
    run_one "$M" "$U" "$ST" "$TAG" $((9900+i)) &
  done
done
wait
echo "level4b power-up done $(date +%H:%M:%S)"
for s in "${SPECS[@]}"; do set -- $s; d=$P/$4; echo "  $4: $(find $d -name 'red_*.csv'|wc -l|tr -d ' ')/$5; health $(cat $d/run_*.log 2>/dev/null | grep -cE 'EDMD-HEALTH|forced_advance|clamp_repair|overlap_repair|wall_overdue'); aborts $(cat $d/run_*.log 2>/dev/null | grep -c ABORTING)"; done
du -sh $P; df -h . | tail -1
