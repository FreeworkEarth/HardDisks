#!/usr/bin/env bash
# ##CHRIS 2026-09-16: Level 1 -- equilibrium wall pressure of the finite box along the compression path.
# Hold-only energy-transfer runs, identical to the W_qs runs except the compartment length:
# L = 37.291667 (1790/48, midpoint) and 35.3125 (1695/48, end of the push): wall positions must sit on the 24 px/sigma grid (2 L x 24 integer) or the validator aborts with initial_wall_position_mismatch. eta = 3.92699/L = 0.105306, 0.111207.
# Per-event log via HD_PISTON_EVENTS; Z_wall from WR and gas-side D0 hits during the 200-sigma hold.
set -uo pipefail
cd /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3
P="experiments_energy_transfer/level1_Zwall_path_20260916"
for L in 37.291667 35.3125; do
  tag="L${L/./p}"; mkdir -p "$P/$tag"
  for k in $(seq 0 19); do
    sd=$(( 9500 + k )); ev="$P/$tag/ev_${sd}.csv"; [ -s "$ev" ] && continue
    ( HD_PISTON_EVENTS="$ev" nice -n 15 ./00ALLINONE --mode=edmd --experiment=energy_transfer \
        --headless --quiet --edmd-acc=0 --seed-drift-order=drift-first \
        --energy-transfer-summary="$P/$tag/summary.csv" --energy-transfer-trace="$P/$tag/tr_${sd}.csv" \
        --particles=100 --particles-boxes=50,50 --particle-radius=0.5 --l0=$L --height=10 \
        --num-walls=1 --wall-positions=$L --wall-mass-factors=1000000000 \
        --piston-right-protocol-mode=step --velocity-right-piston-step=0.005 \
        --max-right-piston-travel=3.93 --auto-piston-step \
        --wall-hold-steps=12000 --steps=10 --fixed-dt=0.4 --energy-measurement \
        --eff-output=wall-ke --kbt1 --seed=$sd >> "$P/$tag/run.log" 2>&1 ) &
    while [ "$(jobs -rp | wc -l)" -ge 3 ]; do sleep 1; done
  done
done
wait
echo "done: $(ls $P/L*/ev_*.csv 2>/dev/null | wc -l | tr -d ' ')/40 event logs; health lines: $(cat $P/L*/run.log | grep -c EDMD-HEALTH)"
