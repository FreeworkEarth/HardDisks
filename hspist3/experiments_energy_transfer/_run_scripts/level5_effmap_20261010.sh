#!/usr/bin/env bash
# ##CHRIS 2026-10-10: EFFICIENCY MAP, geometry C, per 261010 section 1 (committed before this ran).
# GO given 2026-10-12 with amendments A1-A3 (261010 sec. 1.8), committed before launch.
# Binary: v1 + -ffp-contract=off + --version (05215ea line). The piston parks 0.25 sigma outside:
# tau_push = piston_stop_t_rel - 0.25/u, compression start = first nonzero PistonWork, d from piston_target.
# Spring rest length per k so the run STARTS in mechanical equilibrium (k (x_eq - 30.5) = P h = 1.57489):
#   k = 0.25 -> 36.7996   k = 0.5 -> 33.65 (Level 3's value; the equilibrium is 33.6498)   k = 1.0 -> 32.0749
# AMENDED 2026-10-12 (261010 sec. 1.8, A1): run length = 0.25/u + d/u + max(5 spring periods, 3 tau_r), tau_r from
# Mansour friction with ONE gas and M_hat = M_s + N m/3 (3270-9430 sigma-time); steps and table printed by
# validation/paper2_effmap_amend_20261012.py (--emit-specs wrote the SPECS block below). --trace-every unchanged
# (dt <= period/40); no-push control per (k, M_s) at the u = 0.01 length. 336 runs, ~388 M steps, ~1.3 core-h.
# The reduction now also keeps PistonR_v, so Level 3 v6's estimator runs verbatim. Seeds 9500-9507.
# (was: d/u + 5 spring periods, ~16 M steps -- too short for epsilon_settled.)
set -uo pipefail
cd /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3
P=experiments_energy_transfer/level5_effmap_20261010
MAXJOBS=${MAXJOBS:-9}
run_one() {
  local TAG=$1 K=$2 XEQ=$3 M=$4 U=$5 STEPS=$6 EVERY=$7 PUSH=$8 sd=$9
  local d=$P/$TAG; mkdir -p "$d"
  [ -s "$d/red_${sd}.csv" ] && return 0
  local pflags="--piston-right-protocol-mode=step --velocity-right-piston-step=$U --max-right-piston-travel=7.96"
  [ "$PUSH" = 1 ] && pflags="$pflags --auto-piston-step"
  nice -n 5 ./00ALLINONE --mode=edmd --experiment=energy_transfer --headless --quiet \
      --edmd-acc=0 --seed-drift-order=drift-first \
      --energy-transfer-summary="$d/summary_${sd}.csv" \
      --energy-transfer-trace="$d/raw_${sd}.csv" --trace-every=$EVERY \
      --particles=100 --particles-boxes=0,100 --particle-radius=0.5 \
      --l0=54.75 --height=10 --num-walls=1 --wall-positions=30.5 \
      --wall-mass-factors=$M --spring-k-sigma=$K --spring-wall=0 --spring-eq=$XEQ --eff-output=spring \
      $pflags --wall-hold-steps=12000 --steps=$STEPS --fixed-dt=0.4 --kbt1 --seed=$sd \
      > "$d/run_${sd}.log" 2>&1
  local rc=$? health
  health=$(grep -cE 'EDMD-HEALTH|forced_advance|clamp_repair|overlap_repair|wall_overdue' "$d/run_${sd}.log" 2>/dev/null) || true
  health=${health:-0}
  if [ "$rc" -ne 0 ] || [ "$health" -ne 0 ]; then echo "$TAG seed=$sd FAILED rc=$rc health=$health"; return 1; fi
  python3 - "$d/raw_${sd}.csv" "$d/red_${sd}.csv" <<'PY'
import sys, pandas as pd
want = ["Time","KE_gas_total","W0_x_sigma","W0_v","PistonWork","PistonR_x_sigma","PistonR_v","SegCounts","SegEtas","SpringE"]
e = pd.read_csv(sys.argv[1], low_memory=False)
missing = [c for c in want if c not in e.columns]
if missing: sys.exit("REDUCTION ABORTED: columns absent: %s" % ",".join(missing))
e[want].to_csv(sys.argv[2], index=False)
PY
  if [ -s "$d/red_${sd}.csv" ]; then rm -f "$d/raw_${sd}.csv"; else echo "$TAG seed=$sd reduction FAILED"; return 1; fi
}
# TAG k x_eq M u steps every push
SPECS=(
  "k0.25_M50_u0.01 0.25 36.7996 50 0.01 656000 133 1"
  "k0.25_M50_u0.02 0.25 36.7996 50 0.02 631000 133 1"
  "k0.25_M50_u0.05 0.25 36.7996 50 0.05 616000 133 1"
  "k0.25_M50_u0.1 0.25 36.7996 50 0.1 611000 133 1"
  "k0.25_M50_u0.2 0.25 36.7996 50 0.2 609000 133 1"
  "k0.25_M50_u0.5 0.25 36.7996 50 0.5 607000 133 1"
  "ctrl_k0.25_M50 0.25 36.7996 50 0.01 656000 133 0"
  "k0.25_M200_u0.01 0.25 36.7996 200 0.01 1746000 266 1"
  "k0.25_M200_u0.02 0.25 36.7996 200 0.02 1722000 266 1"
  "k0.25_M200_u0.05 0.25 36.7996 200 0.05 1707000 266 1"
  "k0.25_M200_u0.1 0.25 36.7996 200 0.1 1702000 266 1"
  "k0.25_M200_u0.2 0.25 36.7996 200 0.2 1700000 266 1"
  "k0.25_M200_u0.5 0.25 36.7996 200 0.5 1698000 266 1"
  "ctrl_k0.25_M200 0.25 36.7996 200 0.01 1746000 266 0"
  "k0.5_M50_u0.01 0.5 33.65 50 0.01 645000 94 1"
  "k0.5_M50_u0.02 0.5 33.65 50 0.02 621000 94 1"
  "k0.5_M50_u0.05 0.5 33.65 50 0.05 606000 94 1"
  "k0.5_M50_u0.1 0.5 33.65 50 0.1 601000 94 1"
  "k0.5_M50_u0.2 0.5 33.65 50 0.2 598000 94 1"
  "k0.5_M50_u0.5 0.5 33.65 50 0.5 597000 94 1"
  "ctrl_k0.5_M50 0.5 33.65 50 0.01 645000 94 0"
  "k0.5_M200_u0.01 0.5 33.65 200 0.01 1717000 188 1"
  "k0.5_M200_u0.02 0.5 33.65 200 0.02 1692000 188 1"
  "k0.5_M200_u0.05 0.5 33.65 200 0.05 1678000 188 1"
  "k0.5_M200_u0.1 0.5 33.65 200 0.1 1673000 188 1"
  "k0.5_M200_u0.2 0.5 33.65 200 0.2 1670000 188 1"
  "k0.5_M200_u0.5 0.5 33.65 200 0.5 1669000 188 1"
  "ctrl_k0.5_M200 0.5 33.65 200 0.01 1717000 188 0"
  "k1.0_M50_u0.01 1.0 32.0749 50 0.01 639000 66 1"
  "k1.0_M50_u0.02 1.0 32.0749 50 0.02 614000 66 1"
  "k1.0_M50_u0.05 1.0 32.0749 50 0.05 600000 66 1"
  "k1.0_M50_u0.1 1.0 32.0749 50 0.1 595000 66 1"
  "k1.0_M50_u0.2 1.0 32.0749 50 0.2 592000 66 1"
  "k1.0_M50_u0.5 1.0 32.0749 50 0.5 591000 66 1"
  "ctrl_k1.0_M50 1.0 32.0749 50 0.01 639000 66 0"
  "k1.0_M200_u0.01 1.0 32.0749 200 0.01 1699000 133 1"
  "k1.0_M200_u0.02 1.0 32.0749 200 0.02 1675000 133 1"
  "k1.0_M200_u0.05 1.0 32.0749 200 0.05 1660000 133 1"
  "k1.0_M200_u0.1 1.0 32.0749 200 0.1 1655000 133 1"
  "k1.0_M200_u0.2 1.0 32.0749 200 0.2 1653000 133 1"
  "k1.0_M200_u0.5 1.0 32.0749 200 0.5 1651000 133 1"
  "ctrl_k1.0_M200 1.0 32.0749 200 0.01 1699000 133 0"
)
echo "effmap start $(date +%H:%M:%S), binary: $(./00ALLINONE --version | head -1)"
for s in "${SPECS[@]}"; do set -- $s
  for i in $(seq 0 7); do
    while [ "$(jobs -rp | wc -l | tr -d ' ')" -ge "$MAXJOBS" ]; do sleep 3; done
    run_one "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" $((9500+i)) &
  done
done
wait
echo "effmap done $(date +%H:%M:%S)"
for s in "${SPECS[@]}"; do set -- $s; d=$P/$1; echo "  $1: $(find $d -name 'red_*.csv'|wc -l|tr -d ' ')/8; health $(cat $d/run_*.log 2>/dev/null | grep -cE 'EDMD-HEALTH|forced_advance|clamp_repair|overlap_repair|wall_overdue')"; done
du -sh $P; df -h . | tail -1
