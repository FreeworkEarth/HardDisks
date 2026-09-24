#!/usr/bin/env bash
# ##CHRIS 2026-10-06: R-COLLAPSE, PASS 2. Pre-registered in 261006_paper2_level4_Rcollapse2.md,
# which was committed BEFORE this script ran. Read §4 of that file for the predictions; they are
# fixed and this script must not be used to change them.
#
# WHAT PASS 1 GOT WRONG, AND WHOSE FAULT IT WAS. Pass 1 phrased "the variable is R" as
# g = tau_T/M carrying over unchanged. tau_GP is proportional to M x L_c, so doubling N_s at fixed
# eta doubles L_c and doubles tau_GP/M (50.349 M -> 100.698 M). Phrasing the hypothesis in g builds
# the L_c scaling out of the thing being tested. That was an analysis error, not a code error. The
# admissible statement is f = tau_T/tau_GP = f(R), each box with its own L_c.
#
# GEOMETRY -- unchanged from pass 1, and grid-exact by construction:
#   L_c = 2 x 38.75 = 77.5 exactly; wall 78.0; box 156.0; t = 1.0 (recorded).
#   77.5 x 24 = 1860, 156.0 x 24 = 3744, 78.0 x 24 = 1872. All integers.
#   eta = 200 pi r^2/(2 x 77.5 x 10) = 0.10134170 -- IDENTICAL to the ladder's to 8 digits.
#   rho_0 = N_s/L_c = 1.29032 -- identical too. That shared rho_0 is Cencini's path.
#   Z = 1.239880, eta Z' = 0.282875, c_s = 1.749302 (Kolafa-Rottner). tau_GP = 100.698 M.
#
# THE CELLS. cot K = alpha K with alpha = M/(2 N_s m), L_eff = L_c - 2r = 76.5:
#   M    R   alpha     K   period  tau_GP  predA f  predB f  record    steps      every s/per
#   100  1.00 0.500 1.0769 255.16  10069.8  1.52    1.22    1 025 245  61 500 000  690  22.2
#   50   2.00 0.250 1.2646 217.28   5034.9  1.18    0.76      471 445  28 300 000  590  22.1
#   25   4.00 0.125 1.3978 196.57   2517.5  0.98    0.59      212 725  12 800 000  580  20.3
#
# RECORD SIZING. 65 x the largest candidate. For M = 50 and 100 that is the tau MEASURED in pass 1
# (7253, 15773), which is above every prediction -- that is why pass 1 failed L/tau >= 60 and why
# sizing from a prediction is not safe here. For M = 25 there is no pass-1 measurement, so the basis
# is 65 x 1.30 x tau_GP. The plan said 1.05; I used 1.30 because pass 1's new-box f came in 22 %
# above the ladder interpolation at R = 2, and a 5 % margin would repeat the exact failure this pass
# removes. R = 25 is the cheapest cell; the margin costs ~1 core-hour. No prediction was touched.
#
# SEEDS. M = 50 and 100 REUSE pass 1's seeds 9200-9279 with pass 1's --trace-every (590, 690), so
# each pass-2 trace is a STRICT SUPERSET of the pass-1 trace for the same seed. The pass-1 vs pass-2
# difference is then purely a record-length effect with ZERO seed noise -- and record length is
# precisely the systematic under scrutiny. Prediction A comes from the N_s = 50 ladder, which shares
# no data with these seeds, so nothing is tuned to them. M = 25 is new and uses 9400-9479.
#
# COST. 8.21e9 steps, 80 seeds x 3 cells = 240 runs, ~36 core-hours, ~4.0 h at 9 concurrent,
# ~670 MB of reduced traces. Volume had 29 GB free at launch (99 % full) -- fits, nothing deleted.
#
# Resumable: a cell/seed with a non-empty red_*.csv is skipped. Health contract is strict.
set -uo pipefail
cd /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3
P=experiments_energy_transfer/level4_Rcollapse2_20261006
MAXJOBS=${MAXJOBS:-9}

run_one() {
  local M=$1 STEPS=$2 EVERY=$3 sd=$4
  local d=$P/Md${M}
  [ -s "$d/red_${sd}.csv" ] && return 0
  nice -n 5 ./00ALLINONE --mode=edmd --experiment=energy_transfer --headless --quiet \
      --edmd-acc=0 --seed-drift-order=drift-first \
      --energy-transfer-summary="$d/summary_${sd}.csv" \
      --energy-transfer-trace="$d/raw_${sd}.csv" --trace-every=$EVERY \
      --particles=200 --particles-boxes=100,100 --particle-radius=0.5 \
      --l0=78.0 --height=10 --num-walls=1 --wall-positions=78.0 \
      --wall-mass-factors=$M --eff-output=wall-ke \
      --wall-hold-steps=12000 --steps=$STEPS --fixed-dt=0.4 --kbt1 --seed=$sd \
      > "$d/run_${sd}.log" 2>&1
  local rc=$? health
  health=$(grep -cE 'EDMD-HEALTH|forced_advance|clamp_repair|overlap_repair|wall_overdue' "$d/run_${sd}.log" 2>/dev/null) || true
  health=${health:-0}
  if [ "$rc" -ne 0 ] || [ "$health" -ne 0 ]; then
    echo "M=$M seed=$sd FAILED rc=$rc health=$health"; return 1
  fi
  python3 - "$d/raw_${sd}.csv" "$d/red_${sd}.csv" <<'PY'
import sys, pandas as pd
cols = ["Time","KE_gas_left","KE_gas_right","W0_x_sigma"]
pd.read_csv(sys.argv[1], low_memory=False, usecols=cols).to_csv(sys.argv[2], index=False)
PY
  if [ -s "$d/red_${sd}.csv" ]; then rm -f "$d/raw_${sd}.csv"; else
    echo "M=$M seed=$sd reduction FAILED"; return 1; fi
  return 0
}

for M in 100 50 25; do mkdir -p $P/Md${M}; done
echo "Rcollapse2 start $(date +%H:%M:%S), max $MAXJOBS concurrent, 240 runs"
# longest mass first; "M steps every seed0"
for spec in "100 61500000 690 9200" "50 28300000 590 9200" "25 12800000 580 9400"; do
  set -- $spec; M=$1; STEPS=$2; EVERY=$3; S0=$4
  for i in $(seq 0 79); do
    while [ "$(jobs -rp | wc -l | tr -d ' ')" -ge "$MAXJOBS" ]; do sleep 5; done
    run_one "$M" "$STEPS" "$EVERY" $((S0+i)) &
  done
done
wait
echo "Rcollapse2 done $(date +%H:%M:%S)"
for M in 25 50 100; do
  d=$P/Md${M}
  echo "  M=$M: $(find $d -name 'red_*.csv'|wc -l|tr -d ' ')/80; health $(cat $d/run_*.log 2>/dev/null | grep -cE 'EDMD-HEALTH|forced_advance|clamp_repair|overlap_repair|wall_overdue'); aborts $(cat $d/run_*.log 2>/dev/null | grep -c ABORTING)"
done
du -sh $P; df -h . | tail -1
