#!/usr/bin/env bash
# ##CHRIS 2026-09-22: Level 3 v4 -- decide what the +1.0 to +1.9 kT residual is.
#
# v3 left the parameter-free 1-DOF model low by that much at the quasi-static end (2.8-4.2 sigma).
# Two hypotheses, and these runs separate them:
#   A1  the GAS is the problem -- the pressure at the far wall during the push is not the
#       quasi-static F_ad(L). Measure it: hold the spring wall rigid (M = 1e9), read the momentum
#       delivered to it per collision out of the event log (kind D0, column dp), and drive the same
#       ODE with that measured force instead of F_ad.
#   A3  a slow point, u = 0.02, to see whether eps_meas/eps_ODE falls toward 1 as the push slows
#       (supports A1, a dynamic gas effect) or stays flat (supports A2, gas inertia in the wall
#       model, which is analysis-only and needs no runs).
#   A4  thermalised baseline -- let the wall equilibrate with NO push, so the plateau of
#       E_spring + E_wall,kin can be measured instead of assumed.
#
# ASSUMPTION, stated because it halves A1: with the wall held at M = 1e9 it does not move, so the
# gas sees a rigid boundary and F_meas(t) CANNOT depend on M_s. The two M_s values of the v3 grid
# therefore share one held-wall run per speed; 2 cells, not 4.
#
# Geometry is v3 exactly: master box, geometry C, one gas of 100 disks over 78.5 sigma,
# eta = 0.1001, dx = 7.96 (10 %), spring pre-loaded at 33.65 = wall + F/k, k = 0.5 kT/sigma^2.
set -uo pipefail
cd /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3
BIN=./00ALLINONE
P=experiments_energy_transfer/level3_v4_20260922
SEEDS=${SEEDS:-40}
JOBS=${JOBS:-9}

cell(){   # $1 dir  $2 M_s  $3 u  $4 steps  $5.. extra flags
  local d=$P/$1 ms=$2 u=$3 st=$4; shift 4
  mkdir -p "$d"
  local i sd
  for i in $(seq 0 $((SEEDS - 1))); do
    sd=$((9200 + i)); [ -s "$d/tr_${sd}.csv" ] && continue
    ( HD_PISTON_EVENTS="$d/ev_${sd}.csv" nice -n 5 "$BIN" --mode=edmd \
        --experiment=energy_transfer --headless --quiet --edmd-acc=0 \
        --seed-drift-order=drift-first --energy-transfer-summary="$d/summary.csv" \
        --energy-transfer-trace="$d/tr_${sd}.csv" \
        --particles=100 --particles-boxes=0,100 --particle-radius=0.5 \
        --l0=54.75 --height=10 --num-walls=1 --wall-positions=30.5 \
        --wall-mass-factors=$ms --spring-k-sigma=0.5 --spring-wall=0 --spring-eq=33.65 \
        --eff-output=spring --piston-right-protocol-mode=step \
        --velocity-right-piston-step=$u --max-right-piston-travel=7.96 --auto-piston-step \
        --fixed-dt=0.4 --kbt1 --seed=$sd "$@" >> "$d/run.log" 2>&1 ) &
    while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do sleep 0.3; done
  done
  wait
}

# ---- A1: spring wall HELD, measure the force the gas actually delivers to it -----------------
for u in 0.05 0.1; do
  cell "A1_held_u${u}" 1000000000 "$u" 30000 --wall-hold-steps=12000
done

# ---- A3: the slow point, u = 0.02 (push lasts 398 sigma-time = 8.8 sound traversals) ---------
for ms in 50 200; do
  cell "A3_M${ms}_u0.02" "$ms" 0.02 45000 --wall-hold-steps=12000
done

n=$(ls $P/*/tr_*.csv 2>/dev/null | wc -l | tr -d ' ')
echo "A1+A3 done $n/$((4 * SEEDS)); aborts $(cat $P/*/run.log | grep -c ABORTING); health $(cat $P/*/run.log | grep -cE 'EDMD-HEALTH|forced_advance|clamp_repair|overlap_repair|wall_overdue')"

# ---- A4: thermalised baseline. Wall FREE from t = 0 (--wall-hold-steps=0) and NO push at all
# (--auto-piston-step omitted, so the piston never starts). Watch E_spring + E_wall,kin climb to
# equipartition; the plateau is the baseline v3 had to guess, and the settling time is a first
# tau_heat(M) for Level 4. 2000 sigma-time, because at M_s = 1000 the wall period alone is 268.
a4(){ local ms=$1 d=$P/A4_free_M${ms}; mkdir -p "$d"; local i sd
  for i in $(seq 0 $((SEEDS - 1))); do
    sd=$((9200 + i)); [ -s "$d/tr_${sd}.csv" ] && continue
    ( nice -n 5 "$BIN" --mode=edmd --experiment=energy_transfer --headless --quiet --edmd-acc=0 \
        --seed-drift-order=drift-first --energy-transfer-summary="$d/summary.csv" \
        --energy-transfer-trace="$d/tr_${sd}.csv" --particles=100 --particles-boxes=0,100 \
        --particle-radius=0.5 --l0=54.75 --height=10 --num-walls=1 --wall-positions=30.5 \
        --wall-mass-factors=$ms --spring-k-sigma=0.5 --spring-wall=0 --spring-eq=33.65 \
        --eff-output=spring --wall-hold-steps=0 --steps=120000 --fixed-dt=0.4 --kbt1 \
        --seed=$sd >> "$d/run.log" 2>&1 ) &
    while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do sleep 0.3; done
  done; wait; }
for ms in 2 10 50 200 1000; do a4 "$ms"; done
echo "A4 done $(ls $P/A4_free_M*/tr_*.csv 2>/dev/null | wc -l | tr -d ' ')/$((5 * SEEDS))"
