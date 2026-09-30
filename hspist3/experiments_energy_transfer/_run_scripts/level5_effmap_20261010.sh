#!/usr/bin/env bash
# ##CHRIS 2026-10-10: EFFICIENCY MAP, geometry C, per 261010 section 1 (committed before this ran).
# NOT LAUNCHED until Chris and the plan author have read section 1 and given the go.
# Binary: v1 + -ffp-contract=off + --version (05215ea line). The piston parks 0.25 sigma outside:
# tau_push = piston_stop_t_rel - 0.25/u, compression start = first nonzero PistonWork, d from piston_target.
# Spring rest length per k so the run STARTS in mechanical equilibrium (k (x_eq - 30.5) = P h = 1.57489):
#   k = 0.25 -> 36.7996   k = 0.5 -> 33.65 (Level 3's value; the equilibrium is 33.6498)   k = 1.0 -> 32.0749
# Run length d/u + 5 spring periods; --trace-every so dt <= period/40; no-push control per (k, M_s) at the
# u = 0.01 length. 36 cells x 8 seeds + 6 controls x 8 = 336 runs, ~16 M steps. Seeds 9500-9507.
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
want = ["Time","KE_gas_total","W0_x_sigma","W0_v","PistonWork","PistonR_x_sigma","SegCounts","SegEtas","SpringE"]
e = pd.read_csv(sys.argv[1], low_memory=False)
missing = [c for c in want if c not in e.columns]
if missing: sys.exit("REDUCTION ABORTED: columns absent: %s" % ",".join(missing))
e[want].to_csv(sys.argv[2], index=False)
PY
  if [ -s "$d/red_${sd}.csv" ]; then rm -f "$d/raw_${sd}.csv"; else echo "$TAG seed=$sd reduction FAILED"; return 1; fi
}
# TAG k x_eq M u steps every push
SPECS=(
  "k0.25_M50_u0.01 0.25 36.7996 50 0.01 75000 133 1"
  "k0.25_M50_u0.02 0.25 36.7996 50 0.02 51000 133 1"
  "k0.25_M50_u0.05 0.25 36.7996 50 0.05 37000 133 1"
  "k0.25_M50_u0.1 0.25 36.7996 50 0.1 32000 133 1"
  "k0.25_M50_u0.2 0.25 36.7996 50 0.2 30000 133 1"
  "k0.25_M50_u0.5 0.25 36.7996 50 0.5 28000 133 1"
  "ctrl_k0.25_M50 0.25 36.7996 50 0.01 75000 133 0"
  "k0.25_M200_u0.01 0.25 36.7996 200 0.01 102000 266 1"
  "k0.25_M200_u0.02 0.25 36.7996 200 0.02 78000 266 1"
  "k0.25_M200_u0.05 0.25 36.7996 200 0.05 63000 266 1"
  "k0.25_M200_u0.1 0.25 36.7996 200 0.1 59000 266 1"
  "k0.25_M200_u0.2 0.25 36.7996 200 0.2 56000 266 1"
  "k0.25_M200_u0.5 0.25 36.7996 200 0.5 55000 266 1"
  "ctrl_k0.25_M200 0.25 36.7996 200 0.01 102000 266 0"
  "k0.5_M50_u0.01 0.5 33.65 50 0.01 67000 94 1"
  "k0.5_M50_u0.02 0.5 33.65 50 0.02 43000 94 1"
  "k0.5_M50_u0.05 0.5 33.65 50 0.05 29000 94 1"
  "k0.5_M50_u0.1 0.5 33.65 50 0.1 24000 94 1"
  "k0.5_M50_u0.2 0.5 33.65 50 0.2 22000 94 1"
  "k0.5_M50_u0.5 0.5 33.65 50 0.5 20000 94 1"
  "ctrl_k0.5_M50 0.5 33.65 50 0.01 67000 94 0"
  "k0.5_M200_u0.01 0.5 33.65 200 0.01 86000 188 1"
  "k0.5_M200_u0.02 0.5 33.65 200 0.02 62000 188 1"
  "k0.5_M200_u0.05 0.5 33.65 200 0.05 48000 188 1"
  "k0.5_M200_u0.1 0.5 33.65 200 0.1 43000 188 1"
  "k0.5_M200_u0.2 0.5 33.65 200 0.2 41000 188 1"
  "k0.5_M200_u0.5 0.5 33.65 200 0.5 39000 188 1"
  "ctrl_k0.5_M200 0.5 33.65 200 0.01 86000 188 0"
  "k1.0_M50_u0.01 1.0 32.0749 50 0.01 62000 66 1"
  "k1.0_M50_u0.02 1.0 32.0749 50 0.02 38000 66 1"
  "k1.0_M50_u0.05 1.0 32.0749 50 0.05 23000 66 1"
  "k1.0_M50_u0.1 1.0 32.0749 50 0.1 19000 66 1"
  "k1.0_M50_u0.2 1.0 32.0749 50 0.2 16000 66 1"
  "k1.0_M50_u0.5 1.0 32.0749 50 0.5 15000 66 1"
  "ctrl_k1.0_M50 1.0 32.0749 50 0.01 62000 66 0"
  "k1.0_M200_u0.01 1.0 32.0749 200 0.01 75000 133 1"
  "k1.0_M200_u0.02 1.0 32.0749 200 0.02 51000 133 1"
  "k1.0_M200_u0.05 1.0 32.0749 200 0.05 37000 133 1"
  "k1.0_M200_u0.1 1.0 32.0749 200 0.1 32000 133 1"
  "k1.0_M200_u0.2 1.0 32.0749 200 0.2 30000 133 1"
  "k1.0_M200_u0.5 1.0 32.0749 200 0.5 28000 133 1"
  "ctrl_k1.0_M200 1.0 32.0749 200 0.01 75000 133 0"
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
