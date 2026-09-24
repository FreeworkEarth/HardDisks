#!/usr/bin/env bash
# ##CHRIS 2026-10-03: the Level 4 mass ladder -- tau_T(M), tau_r(M) and the mode period across
# alpha = 0.2 to 2.0. Same protocol as the M_d = 10 run of 261002; that mass is NOT repeated here,
# its 80 seeds are already on disk and go into the ladder fit as the fifth point.
#
# WHAT THE LADDER TESTS that one mass cannot: the EXPONENT. Gruber-Piasecki's statement is that
# stage 2 is linear in M; 261002 confirmed the PREFACTOR at one mass (tau_T = 491 +- 29 against
# 50.35 M = 503.5, excluding the ideal 61.84 M = 618.4 at 3.5 sigma) but says nothing about b in
# tau_T = a M^b. Five masses do.
#
# GEOMETRY of the divider mode, from cot K = alpha K with alpha = M/(2 N_s m) = M/100, KR c_s =
# 1.744337 and L_eff = 38.75 - 2r = 37.75. Every K below satisfies its equation to <= 1e-12:
#
#   M    alpha   K        period   tau_T    record      steps       N     dt     samp/period  MB
#   20   0.2     1.3138   103.5     1007     65 455     3 929 226    50   0.83      124       4.10
#   50   0.5     1.0769   126.3     2518    163 638     9 823 065   110   1.83       69       4.66
#  100   1.0     0.8603   158.1     5035    327 275    19 646 129   210   3.50       45       4.88
#  200   2.0     0.6533   208.1    10070    654 550    39 292 259   420   7.00       30       4.88
#
# N is the smallest tidy value keeping each reduced file <= 5 MB; it still leaves >= 30 samples per
# divider-mode period everywhere, against the 19 that already gave a clean two-component fit at
# M = 10. Disk: 4 masses x 80 seeds x ~4.6 MB = ~1.5 GB of reduced files, plus one ~35 MB raw
# trace per concurrent job, deleted as soon as its reduction succeeds.
#
# COST, from the measured ~10 s per 2 M steps: 0.44 + 1.09 + 2.18 + 4.37 = 8.1 core-hours.
# At 9 concurrent jobs that is ~1 h wall clock. Longest mass first, so the tail starts early.
#
# RESUMABLE and HEALTH-GATED: a seed whose reduced file exists is skipped; a seed with a non-zero
# health count or exit code keeps its raw trace for diagnosis and is NOT reduced.
set -uo pipefail
cd /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3
P=experiments_energy_transfer/level4_ladder_20261003
MAXJOBS=${MAXJOBS:-9}

run_one() {                      # $1 = M, $2 = steps, $3 = trace-every, $4 = seed
  local M=$1 STEPS=$2 EVERY=$3 sd=$4
  local d=$P/Md${M}
  [ -s "$d/red_${sd}.csv" ] && return 0
  nice -n 5 ./00ALLINONE --mode=edmd --experiment=energy_transfer --headless --quiet \
      --edmd-acc=0 --seed-drift-order=drift-first \
      --energy-transfer-summary="$d/summary_${sd}.csv" \
      --energy-transfer-trace="$d/raw_${sd}.csv" --trace-every=$EVERY \
      --particles=100 --particles-boxes=50,50 --particle-radius=0.5 \
      --l0=39.25 --height=10 --num-walls=1 --wall-positions=39.25 \
      --wall-mass-factors=$M --eff-output=wall-ke \
      --wall-hold-steps=12000 --steps=$STEPS --fixed-dt=0.4 --kbt1 --seed=$sd \
      > "$d/run_${sd}.log" 2>&1
  local rc=$?
  local health
  health=$(grep -cE 'EDMD-HEALTH|forced_advance|clamp_repair|overlap_repair|wall_overdue' "$d/run_${sd}.log" 2>/dev/null) || true
  health=${health:-0}
  if [ "$rc" -ne 0 ] || [ "$health" -ne 0 ]; then
    echo "M=$M seed=$sd FAILED rc=$rc health=$health -- raw trace KEPT"; return 1
  fi
  python3 - "$d/raw_${sd}.csv" "$d/red_${sd}.csv" <<'PY'
import sys, pandas as pd
cols = ["Time","KE_gas_left","KE_gas_right","W0_x_sigma"]
pd.read_csv(sys.argv[1], low_memory=False, usecols=cols).to_csv(sys.argv[2], index=False)
PY
  if [ -s "$d/red_${sd}.csv" ]; then rm -f "$d/raw_${sd}.csv"; else
    echo "M=$M seed=$sd reduction FAILED -- raw trace KEPT"; return 1; fi
  return 0
}

# longest mass first so the critical path starts immediately
for M in 200 100 50 20; do mkdir -p $P/Md${M}; done
echo "ladder start $(date +%H:%M:%S), max $MAXJOBS concurrent"
for spec in "200 39292259 420" "100 19646129 210" "50 9823065 110" "20 3929226 50"; do
  set -- $spec; M=$1; STEPS=$2; EVERY=$3
  for i in $(seq 0 79); do
    while [ "$(jobs -rp | wc -l | tr -d ' ')" -ge "$MAXJOBS" ]; do sleep 5; done
    run_one "$M" "$STEPS" "$EVERY" $((9200+i)) &
  done
done
wait
echo "ladder done $(date +%H:%M:%S)"
for M in 20 50 100 200; do
  d=$P/Md${M}
  echo "  M=$M: $(ls $d/red_*.csv 2>/dev/null | wc -l | tr -d ' ')/80 reduced; \
health $(cat $d/run_*.log 2>/dev/null | grep -cE 'EDMD-HEALTH|forced_advance|clamp_repair|overlap_repair|wall_overdue'); \
aborts $(cat $d/run_*.log 2>/dev/null | grep -c ABORTING); \
raw left $(ls $d/raw_*.csv 2>/dev/null | wc -l | tr -d ' ')"
done
du -sh $P
