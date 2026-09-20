#!/usr/bin/env bash
# ##CHRIS 2026-09-20: Level 3 in the master box -- can a spring catch the ACOUSTIC pulse?
#
# The 2026-09-19 pilot asked the spring to catch the quasi-static part. It cannot: gas and spring
# in series give E_max = k_gas dx^2/8 = N kT f^2 (Z + eta Z')/8, second order in the compression
# fraction where W_qs is first order, so 0.098 kT at f = 0.1 against a kT/2 thermal pedestal.
# See validation/level3_design_20260920.py. What IS above the pedestal is what Level 2 measured a
# stepped piston to launch beyond W_qs: A u^2 with A_step = 25.9 +- 4.3, i.e. 1.0 / 6.5 / 25.9 kT
# at u = 0.2 / 0.5 / 1.0. So the spring is asked to catch that instead.
#
# ALL k HERE ARE IN kT/sigma^2 (--spring-k-sigma). The old --spring-k is per PIXEL^2, a factor
# 24^2 = 576 smaller, which is what wrecked the pilot's design: its "k = 5" was 2880 kT/sigma^2,
# so T_w was 1.66 sigma-time, not the 39.5 the pilot reported, and every cell it believed was
# impulsive was in fact deep in the quasi-static regime.
#
# Geometry C in the master box: 30 sigma empty spring compartment | wall_S on the spring (M_s) |
# one gas of 100 disks over 78.5 sigma (eta = 0.1001) | piston. dx = 7.96 sigma = 10 %.
# Box 0 .. 109.50 (--l0=54.75), wall_S centre 30.5. The spring compartment is 30 sigma, not 10,
# because at the soft end the wall's thermal excursion is sqrt(kT/(k+k_gas)) = 4.5 sigma and it
# needs the room; Level 0b at 40 seeds shows the gas does not notice the change (0.83 sigma).
#
# ARM (b), the discriminator: M_s is scanned at FIXED k. The series formula does not contain M_s
# at all, so it predicts a flat line; the impedance formula T_imp = 4 Z_w Z_g/(Z_w+Z_g)^2 with
# Z_w = sqrt(M_s(k+k_gas)) and Z_g = rho c_s H = 2.222 predicts an interior maximum. At k = 0.5
# the match is at M_s = 9.9, so the ladder 2/10/50/200/1000 brackets it: T_imp = 0.86 / 1.00 /
# 0.84 / 0.59 / 0.32. The u ladder is nested inside because the acoustic term grows as u^2 and the
# quasi-static one does not -- a second, independent discriminator at no extra cost.
#
# ARM (a) WAS DROPPED, and why is a result. The anchor cell k = k_gas aborted 2/2 seeds with
# [wall_boundary_contact]: the wall was pushed into the outer boundary. It is not under-roomed, it
# is impossible. A spring anchored at the wall's start position must hold the STANDING gas force
# F = N kT Z / L, so the wall sits a distance F/k from the anchor, and at k = k_gas
#
#     F / k_gas = (N kT Z/L) / (N kT (Z + eta Z')/L^2) = L * Z/(Z + eta Z') = 0.816 L ,
#
# INDEPENDENT of N -- always about eight tenths of the gas length. So the stiffness that maximises
# quasi-static capture is exactly the stiffness at which the spring cannot be statically balanced
# inside its own apparatus. Pre-loading the spring to compensate does not save it either: the
# baseline then couples linearly, and the pedestal fluctuation F*sqrt(kT/(k+k_gas)) becomes 7.1 kT
# at k = k_gas against a 0.195 kT signal. Measured offsets: 64.1 / 15.8 / 3.15 / 0.79 / 0.31 sigma
# at k = 0.0246 / 0.1 / 0.5 / 2 / 5. k = 0.5 sits at 3.15 sigma, comfortably inside 30.
set -uo pipefail
cd /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3
BIN=./00ALLINONE
P=experiments_energy_transfer/level3_master_20260920
SEEDS=${SEEDS:-40}
JOBS=${JOBS:-9}

run_cell(){                       # $1 = cell dir, $2 = k [kT/sigma^2], $3 = M_s, $4 = u, $5 = steps
  local d=$P/$1 k=$2 ms=$3 u=$4 st=$5 sd
  mkdir -p "$d"
  for i in $(seq 0 $((SEEDS - 1))); do
    sd=$((9200 + i)); [ -s "$d/tr_${sd}.csv" ] && continue
    ( HD_PISTON_EVENTS="$d/ev_${sd}.csv" nice -n 5 "$BIN" --mode=edmd \
        --experiment=energy_transfer --headless --quiet --edmd-acc=0 \
        --seed-drift-order=drift-first --energy-transfer-summary="$d/summary.csv" \
        --energy-transfer-trace="$d/tr_${sd}.csv" \
        --particles=100 --particles-boxes=0,100 --particle-radius=0.5 \
        --l0=54.75 --height=10 --num-walls=1 --wall-positions=30.5 \
        --wall-mass-factors=$ms --spring-k-sigma=$k --spring-wall=0 --spring-eq=30.5 \
        --eff-output=spring --piston-right-protocol-mode=step \
        --velocity-right-piston-step=$u --max-right-piston-travel=7.96 --auto-piston-step \
        --wall-hold-steps=12000 --steps=$st --fixed-dt=0.4 --kbt1 --seed=$sd \
        >> "$d/run.log" 2>&1 ) &
    while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do sleep 0.3; done
  done
  wait
}

# arm (b): the M_s ladder at k = 0.5 kT/sigma^2, u nested inside
for ms in 2 10 50 200 1000; do
  for u in 0.2 0.5 1.0; do
    case $u in 0.2) st=30000 ;; 0.5) st=28000 ;; *) st=27000 ;; esac
    run_cell "k0.5_M${ms}_u${u}" 0.5 "$ms" "$u" "$st"
  done
done

# arm (a) anchor deliberately NOT run -- see the header. Kept here, commented, so the decision is
# visible rather than silently absent:
#   run_cell "kgas_M200_u0.2" 0.02457 200 0.2 30000

n=$(ls $P/*/tr_*.csv 2>/dev/null | wc -l | tr -d ' ')
echo "done $n/$(( 5 * 3 * SEEDS )); aborts $(cat $P/*/run.log | grep -c ABORTING); health $(cat $P/*/run.log | grep -cE 'EDMD-HEALTH|forced_advance|clamp_repair|overlap_repair|wall_overdue')"
