#!/usr/bin/env bash
# ##CHRIS 2026-09-29: tau_T from EQUILIBRIUM fluctuations. No temperature step, no piston.
#
# The 0.20 per-seed spread of T_1 - T_2 is not measurement noise but the microcanonical energy-split
# fluctuation: hard disks have no potential energy, so a compartment's kinetic temperature IS its
# energy, f = E_1/E ~ Beta(N,N), sigma_f = 1/(2 sqrt(2N+1)) = 0.0497 and sigma(T_1-T_2) = 4 T
# sigma_f = 0.199 (measured 0.204, 0.179). Relaxation and fluctuation are one OU process, so the
# stationary autocorrelation carries the same tau_T and no step is needed.
#
# RECORD LENGTH MATTERS MORE THAN SEEDS. A synthetic OU control with this estimator shows the
# finite-record bias is -7 % at L/tau = 8.8, -37 % at 5.5 and -67 % at 3.6. Paper 1's records sit
# at 3.6-8.8 and cannot discriminate 50.4M from 61.8M. Here L/tau = 20 at M_d = 10 (tau ~ 504).
#
# DISK: the energy-transfer trace writes ONE ROW PER STEP and ignores --output-dt, so a 600k-step
# run is ~220 MB. A first attempt at 3M steps x 20 seeds reached 21 GB and 2.6 GB free before it
# was killed. Each seed is therefore REDUCED to the four needed columns at 1-in-300 and the raw
# trace deleted immediately; peak disk is one trace.
set -uo pipefail
cd /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3
P=experiments_energy_transfer/level4_equilibrium_20260929
d=$P/Md10; mkdir -p "$d"
for i in $(seq 0 19); do
  sd=$((9200+i))
  [ -s "$d/red_${sd}.csv" ] && continue
  nice -n 5 ./00ALLINONE --mode=edmd --experiment=energy_transfer --headless --quiet \
      --edmd-acc=0 --seed-drift-order=drift-first \
      --energy-transfer-summary="$d/summary.csv" --energy-transfer-trace="$d/raw_${sd}.csv" \
      --particles=100 --particles-boxes=50,50 --particle-radius=0.5 \
      --l0=39.25 --height=10 --num-walls=1 --wall-positions=39.25 \
      --wall-mass-factors=10 --eff-output=wall-ke \
      --wall-hold-steps=12000 --steps=612000 --fixed-dt=0.4 --kbt1 --seed=$sd \
      >> "$d/run.log" 2>&1
  python3 - "$d/raw_${sd}.csv" "$d/red_${sd}.csv" <<'PY'
import sys, pandas as pd
src, dst = sys.argv[1], sys.argv[2]
cols = ["Time","KE_gas_left","KE_gas_right","W0_x_sigma"]
pd.read_csv(src, low_memory=False, usecols=cols).iloc[::300].to_csv(dst, index=False)
PY
  rm -f "$d/raw_${sd}.csv"
  echo "seed $sd reduced -> $(ls -la "$d/red_${sd}.csv" | awk '{print $5}') bytes"
done
echo "equilibrium done $(ls $d/red_*.csv | wc -l | tr -d ' ')/20; aborts $(grep -c ABORTING $d/run.log); health $(grep -cE 'EDMD-HEALTH|forced_advance|clamp_repair|overlap_repair|wall_overdue' $d/run.log); disk $(du -sh $P | cut -f1)"
