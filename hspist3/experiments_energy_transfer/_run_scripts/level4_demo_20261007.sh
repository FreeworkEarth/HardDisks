#!/usr/bin/env bash
# ##CHRIS 2026-10-07: THE DEMONSTRATION RUN -- work in on the right, and what actually arrives on
# the left. Pre-registered in 261007_paper2_level4_Rcollapse2_demo.md section 3 BEFORE any fit.
#
# GEOMETRY (the ladder's own box, so tau_T is already measured here):
#   N_s = 50 per side, L_c = 38.75, wall 39.25, box 78.5, t = 1.0, eta = 0.10134170, M_d = 200.
#   Piston on the RIGHT, step protocol, travel 3.875 sigma = 10 % of L_c.
#   Grid: 3.875 x 24 = 93 EXACT; right wall 78.5 -> 74.625, x24 = 1791 EXACT.
#
# WHY u = 0.2 IS FAST ENOUGH, and it is NOT the tau_v argument in the plan. The plan said the push
# (19.38 sigma) is short against tau_v = 486. tau_v = M/gamma = 200/4.118 = 48.6, not 486 -- gamma
# = 4.118 is the two-sided kinetic friction already used for the recollision parameter. The correct
# and stronger argument is CAUSAL: the compression cannot reach the divider before
#     L_c/c_s = 38.75/1.7493 = 22.15 sigma-time,
# which is AFTER the push has ended at 19.38. The divider is not merely slow to respond, it has not
# yet been told. That is why the right gas takes the whole of W_in and the v3 runs (u = 0.05,
# 77.5 sigma push, 3.5 sound transits) did not.
#
# THE CELLS
#   B1long  u = 0.2, 80 seeds, 7 250 000 steps = 120 833 sigma = 3.0 tau_T(200), every 600 (dt = 10)
#   B1zoom  u = 0.2, 80 seeds,     6 000 steps =     100 sigma, every 1  (dt = 1/60) -- the PUSH
#           itself at full resolution. Same seeds and same flags, so it is a strict PREFIX of
#           B1long: --steps only sets where the run stops. At dt = 10 the 19.38-sigma push is two
#           samples, which is no use for the first panel.
#   B2pilot u = 1.0, 8 seeds, as B1long. See section 3.4: at u = 0.2 the push is near-reversible and
#           stage 1 removes the ENTIRE temperature difference, so stage 2 has nothing to decay.
#           u = 1.0 is Mach 0.57 -- strongly irreversible but NOT a shock. 8 seeds is a pilot to
#           measure the surviving residual before anyone spends 80.
#
# COST ~2.8 core-hours total. Run AFTER the R-collapse campaign: this machine has 12 performance
# cores and campaign A already holds 9 of them, so launching alongside would not finish the pair any
# sooner -- it would only push A out.
set -uo pipefail
cd /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3
P=experiments_energy_transfer/level4_demo_20261007
MAXJOBS=${MAXJOBS:-9}

run_one() {
  local tag=$1 U=$2 STEPS=$3 EVERY=$4 sd=$5
  local d=$P/$tag
  [ -s "$d/red_${sd}.csv" ] && return 0
  nice -n 5 ./00ALLINONE --mode=edmd --experiment=energy_transfer --headless --quiet \
      --edmd-acc=0 --seed-drift-order=drift-first \
      --energy-transfer-summary="$d/summary_${sd}.csv" \
      --energy-transfer-trace="$d/raw_${sd}.csv" --trace-every=$EVERY \
      --particles=100 --particles-boxes=50,50 --particle-radius=0.5 \
      --l0=39.25 --height=10 --num-walls=1 --wall-positions=39.25 \
      --wall-mass-factors=200 --eff-output=wall-ke \
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
cols = ["Time","KE_gas_left","KE_gas_right","W0_x_sigma"]
e = pd.read_csv(sys.argv[1], low_memory=False)
keep = [c for c in cols if c in e.columns]
extra = [c for c in ("PistonWork","W_in","RightPiston_x_sigma") if c in e.columns]
e[keep+extra].to_csv(sys.argv[2], index=False)
PY
  if [ -s "$d/red_${sd}.csv" ]; then rm -f "$d/raw_${sd}.csv"; else
    echo "$tag seed=$sd reduction FAILED"; return 1; fi
  return 0
}

for t in B1long B1zoom B2pilot; do mkdir -p $P/$t; done
echo "demo start $(date +%H:%M:%S), max $MAXJOBS concurrent"
# "tag u steps every nseeds seed0"
for spec in "B1long 0.2 7250000 600 80 9600" "B2pilot 1.0 7250000 600 8 9800" "B1zoom 0.2 6000 1 80 9600"; do
  set -- $spec; tag=$1; U=$2; ST=$3; EV=$4; NS=$5; S0=$6
  for i in $(seq 0 $((NS-1))); do
    while [ "$(jobs -rp | wc -l | tr -d ' ')" -ge "$MAXJOBS" ]; do sleep 5; done
    run_one "$tag" "$U" "$ST" "$EV" $((S0+i)) &
  done
done
wait
echo "demo done $(date +%H:%M:%S)"
for t in B1long B2pilot B1zoom; do
  d=$P/$t
  echo "  $t: $(find $d -name 'red_*.csv'|wc -l|tr -d ' ') done; health $(cat $d/run_*.log 2>/dev/null | grep -cE 'EDMD-HEALTH|forced_advance|clamp_repair|overlap_repair|wall_overdue'); aborts $(cat $d/run_*.log 2>/dev/null | grep -c ABORTING)"
done
du -sh $P; df -h . | tail -1
