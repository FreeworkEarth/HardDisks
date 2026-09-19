#!/usr/bin/env bash
# ##CHRIS 2026-09-13: resume after the battery outage, then finish the weekend queue.
#   0. hold until the machine is on AC power (or the ALLOW_BATTERY flag file exists)
#   1. GATE: a lone exact-seed run must reproduce a completed top-up run byte-for-byte
#   2. resume the A2 top-up part B: only the missing N = 2500 runs
#   3. regenerate 260913_A2_* from famB + top-up
#   4. alpha = 2 cells, 10 repeats, resumable run by run
#   5. regenerate 260913_A2_with_alpha2_*
#   6. A4 particle ladder
# Launched detached from the Claude Code session (new process session) and under caffeinate,
# so neither the session ending nor idle sleep stops it. A hard power loss still would; every
# step skips finished work, so relaunching this script continues where it stopped.
set -uo pipefail
REPO=/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks
ROOT="$REPO/hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN"
STATE="$ROOT/_orchestration_20260913"
TOPUP="$ROOT/A2_topup_20260912"
ALPHA="$ROOT/A2_alpha2_20260912"
A4="$ROOT/A4_particle_ladder_20260912"
PLOTS="$REPO/0000_PLAN_OVERALL/ALL_MARKDOWNS/260909_plots"
PY=/opt/homebrew/bin/python3
FLAG="$STATE/ALLOW_BATTERY"
mkdir -p "$STATE"
cd "$REPO/hspist3" || exit 1

ts() { date '+%F %T %Z'; }
wait_ac() {
  local said=0
  while ! pmset -g batt | grep -q "AC Power" && [ ! -f "$FLAG" ]; do
    if [ "$said" -eq 0 ]; then echo "[$(ts)] on battery -- waiting for AC power (or: touch $FLAG)"; said=1; fi
    sleep 60
  done
  if [ "$said" -eq 1 ]; then echo "[$(ts)] power OK"; fi
  return 0
}

echo "[$(ts)] chain start, pid $$"
wait_ac

if [ ! -f "$STATE/gate.pass" ]; then
  echo "[$(ts)] step 1: exact-seed identity gate"
  if "$PY" validation/resume_runs.py --gate "$TOPUP/eta_0p30/N900/m_50" --gate-run 1 --gate-out "$STATE/gate_exactseed"; then
    touch "$STATE/gate.pass"
  else
    echo "[$(ts)] GATE FAILED -- stopping before any resume; nothing was written to the campaign"
    exit 1
  fi
fi

echo "[$(ts)] step 2: resume A2 top-up part B"
"$PY" validation/resume_runs.py --plan topup-b --root "$TOPUP" --jobs 10 --allow-battery-flag "$FLAG"
echo "[$(ts)] step 2 returned $?"

echo "[$(ts)] step 3: regenerate 260913_A2 (famB + top-up)"
"$PY" validation/build_A2_trace_cache.py "$ROOT/famB_20260911,$TOPUP" "$STATE/A2_cache_famB_topup.csv" 10
"$PY" validation/analyze_A2_cut_sensitivity.py "$STATE/A2_cache_famB_topup.csv" "$PLOTS/260913_A2"
echo "[$(ts)] step 3 returned $?"

wait_ac
echo "[$(ts)] step 4: alpha = 2 cells, 10 repeats"
"$PY" validation/resume_runs.py --plan alpha --root "$ALPHA" --repeats 10 --jobs 10 --allow-battery-flag "$FLAG"
echo "[$(ts)] step 4 returned $?"

echo "[$(ts)] step 5: regenerate with the alpha = 2 cells"
"$PY" validation/build_A2_trace_cache.py "$ROOT/famB_20260911,$TOPUP,$ALPHA" "$STATE/A2_cache_with_alpha2.csv" 10
"$PY" validation/analyze_A2_cut_sensitivity.py "$STATE/A2_cache_with_alpha2.csv" "$PLOTS/260913_A2_with_alpha2"
echo "[$(ts)] step 5 returned $?"

wait_ac
echo "[$(ts)] step 6: A4 particle ladder"
bash validation/run_A4_particle_ladder.sh "$A4" 10
echo "[$(ts)] step 6 returned $?"
echo "[$(ts)] CHAIN DONE"
