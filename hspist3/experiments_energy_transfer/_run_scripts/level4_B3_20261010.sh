#!/usr/bin/env bash
# ##CHRIS 2026-10-10: B3 -- the out-of-equilibrium test of tau_T. Pre-registered in 261007 section
# 3B (committed before this ran; the chain script refuses to start it until the marker /tmp/B3_go
# exists, which is created only after that commit). Geometry and cadence identical to B1long except
# u = 1.0: ladder box, N_s = 50, M_d = 200, d = 3.875, 80 seeds, 3 tau_T = 7 250 000 steps,
# --trace-every=600 (dt = 10). Binary: v1 + -ffp-contract=off (05215ea). Seeds 9800-9879 (the first
# 8 coincide with B2pilot's, on a different binary; not mixed).
set -uo pipefail
cd /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3
P=experiments_energy_transfer/level4_B3_20261010/B3
MAXJOBS=${MAXJOBS:-9}; mkdir -p "$P"
run_one() {
  local sd=$1
  [ -s "$P/red_${sd}.csv" ] && return 0
  nice -n 5 ./00ALLINONE --mode=edmd --experiment=energy_transfer --headless --quiet \
      --edmd-acc=0 --seed-drift-order=drift-first \
      --energy-transfer-summary="$P/summary_${sd}.csv" \
      --energy-transfer-trace="$P/raw_${sd}.csv" --trace-every=600 \
      --particles=100 --particles-boxes=50,50 --particle-radius=0.5 \
      --l0=39.25 --height=10 --num-walls=1 --wall-positions=39.25 \
      --wall-mass-factors=200 --eff-output=wall-ke \
      --piston-right-protocol-mode=step --velocity-right-piston-step=1.0 \
      --max-right-piston-travel=3.875 --auto-piston-step \
      --wall-hold-steps=12000 --steps=7250000 --fixed-dt=0.4 --kbt1 --seed=$sd \
      > "$P/run_${sd}.log" 2>&1
  local rc=$? health
  health=$(grep -cE 'EDMD-HEALTH|forced_advance|clamp_repair|overlap_repair|wall_overdue' "$P/run_${sd}.log" 2>/dev/null) || true
  health=${health:-0}
  if [ "$rc" -ne 0 ] || [ "$health" -ne 0 ]; then echo "B3 seed=$sd FAILED rc=$rc health=$health"; return 1; fi
  python3 - "$P/raw_${sd}.csv" "$P/red_${sd}.csv" <<'PY'
import sys, pandas as pd
want = ["Time","KE_gas_left","KE_gas_right","W0_x_sigma","W0_v","PistonWork","SegCounts","SegEtas"]
e = pd.read_csv(sys.argv[1], low_memory=False)
missing = [c for c in want if c not in e.columns]
if missing: sys.exit("REDUCTION ABORTED: columns absent: %s" % ",".join(missing))
e[want].to_csv(sys.argv[2], index=False)
PY
  if [ -s "$P/red_${sd}.csv" ]; then rm -f "$P/raw_${sd}.csv"; else echo "B3 seed=$sd reduction FAILED"; return 1; fi
}
echo "B3 start $(date +%H:%M:%S), binary: $(./00ALLINONE --version | head -1)"
for i in $(seq 0 79); do
  while [ "$(jobs -rp | wc -l | tr -d ' ')" -ge "$MAXJOBS" ]; do sleep 5; done
  run_one $((9800+i)) &
done
wait
echo "B3 done $(date +%H:%M:%S): $(find $P -name 'red_*.csv'|wc -l|tr -d ' ')/80; health $(cat $P/run_*.log 2>/dev/null | grep -cE 'EDMD-HEALTH|forced_advance|clamp_repair|overlap_repair|wall_overdue'); aborts $(cat $P/run_*.log 2>/dev/null | grep -c ABORTING)"
du -sh $P; df -h . | tail -1
