#!/usr/bin/env bash
# ##CHRIS 2026-10-02: the M_d = 10 cell of the Level 4 equilibrium route at FULL length.
# Not launched by this session -- written so it can be started the night after the meeting.
#
# WHAT IT SETTLES. tau_T = 50.35 M (hard disks) against 61.84 M (ideal) is a 23 % gap. The 2026-09-29
# run (8200 sigma of usable record, 20 seeds) measured tau_T = 490 +- 190 and could not separate
# them. The binding constraint is the estimator's BIAS, which is common-mode across seeds and
# therefore immune to adding seeds -- only a longer record fixes it. Measured separations:
#
#     8200 sigma  / L-over-tau 16 / 20 seeds -> 1.4 sigma
#     8200        / 16            / 80 seeds -> 1.7 sigma    <- seeds alone do NOT do it
#     16400       / 33            / 20 seeds -> 1.6 sigma
#     32800       / 65            / 20 seeds -> 1.9 sigma
#     32800       / 65            / 80 seeds -> 4.0 sigma    <- THIS (modelled estimator)
#                                            -> 3.3 sigma    (block estimator)
#
# so: 65 tau_T of record and 80 seeds, at ONE mass. The other two masses (50, 200) are the KOA
# arrays in cluster/level4_massladder.sbatch; this cell is a Mac job and needs no cluster.
#
# COST. 1,964,913 steps x 80 seeds = 1.57e8 steps, about 13x the 2026-09-29 run. Overnight.
# Run it with the RELEASE binary: `make release` in hspist3. Plain `make` builds the AddressSanitizer
# debug target, which is ~5x slower and not byte-comparable with a release build.
#
# DISK. The energy-transfer trace writes one row per step (~368 bytes) and ignores --output-dt, so
# undecimated this cell would be 723 MB per seed, 58 GB for the campaign. --trace-every=200 gives
# dt = 3.33 sigma-time, which is 29 samples per divider-mode period (the 2026-09-29 run's 19 already
# gave a clean two-component ACF fit), i.e. 3.6 MB per seed and ~290 MB in total before reduction.
# Each seed is still reduced to the four analysed columns and its raw trace deleted immediately.
#
# RESUMABLE: a seed whose reduced file already exists is skipped, so it can be interrupted.
set -uo pipefail
cd /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3
P=experiments_energy_transfer/level4_equilibrium_KOAlength_20261002
d=$P/Md10; mkdir -p "$d"

STEPS=1964913          # 65 * 50.35 * 10 sigma-time * 60.0294 steps per sigma-time (fixed-dt=0.4)
EVERY=200

for i in $(seq 0 79); do
  sd=$((9200+i))
  [ -s "$d/red_${sd}.csv" ] && continue
  nice -n 5 ./00ALLINONE --mode=edmd --experiment=energy_transfer --headless --quiet \
      --edmd-acc=0 --seed-drift-order=drift-first \
      --energy-transfer-summary="$d/summary.csv" --energy-transfer-trace="$d/raw_${sd}.csv" \
      --trace-every=$EVERY \
      --particles=100 --particles-boxes=50,50 --particle-radius=0.5 \
      --l0=39.25 --height=10 --num-walls=1 --wall-positions=39.25 \
      --wall-mass-factors=10 --eff-output=wall-ke \
      --wall-hold-steps=12000 --steps=$STEPS --fixed-dt=0.4 --kbt1 --seed=$sd \
      >> "$d/run_${sd}.log" 2>&1
  rc=$?
  health=$(grep -cE 'EDMD-HEALTH|forced_advance|clamp_repair|overlap_repair|wall_overdue' "$d/run_${sd}.log" 2>/dev/null) || true
  health=${health:-0}
  if [ "$rc" -ne 0 ] || [ "$health" -ne 0 ]; then
    echo "seed $sd FAILED rc=$rc health=$health -- raw trace KEPT for diagnosis"; continue
  fi
  python3 - "$d/raw_${sd}.csv" "$d/red_${sd}.csv" <<'PY'
import sys, pandas as pd
cols = ["Time","KE_gas_left","KE_gas_right","W0_x_sigma"]
pd.read_csv(sys.argv[1], low_memory=False, usecols=cols).to_csv(sys.argv[2], index=False)
PY
  if [ -s "$d/red_${sd}.csv" ]; then rm -f "$d/raw_${sd}.csv"; fi
  echo "seed $sd done -> $(ls -la "$d/red_${sd}.csv" | awk '{print $5}') bytes"
done

echo "equilibrium KOA-length done $(ls $d/red_*.csv 2>/dev/null | wc -l | tr -d ' ')/80; \
aborts $(cat $d/run_*.log 2>/dev/null | grep -c ABORTING); \
health $(cat $d/run_*.log 2>/dev/null | grep -cE 'EDMD-HEALTH|forced_advance|clamp_repair|overlap_repair|wall_overdue'); \
disk $(du -sh $P | cut -f1)"
echo "analyse with: python3 validation/paper2_level4_fd2_acfmodel_20260930.py   (point P= at this dir)"
