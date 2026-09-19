#!/usr/bin/env bash
# ##CHRIS 2026-09-17: Level 1 -- densify and deepen the equilibrium wall-pressure measurement along the
# compression path. Grid-exact compartment lengths (L = k/48, so 2 L x 24 is an integer):
#   38.270833 (1837/48) and 36.3125 (1743/48): new points, 100 seeds each (9600-9699)
#   37.291667 (1790/48) and 35.3125 (1695/48): 80 more seeds each (9520-9599), on top of 9500-9519
# Hold-only runs, identical to level1_zwall_path.sh otherwise.
set -uo pipefail
cd /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3
P="experiments_energy_transfer/level1_Zwall_path_20260916"
run_point(){ L=$1; s0=$2; n=$3
  tag="L${L/./p}"; mkdir -p "$P/$tag"
  for k in $(seq 0 $((n-1))); do
    sd=$(( s0 + k )); ev="$P/$tag/ev_${sd}.csv"; [ -s "$ev" ] && [ "$(wc -l < "$ev")" -gt 100 ] && continue
    ( HD_PISTON_EVENTS="$ev" nice -n 10 ./00ALLINONE --mode=edmd --experiment=energy_transfer \
        --headless --quiet --edmd-acc=0 --seed-drift-order=drift-first \
        --energy-transfer-summary="$P/$tag/summary.csv" --energy-transfer-trace="$P/$tag/tr_${sd}.csv" \
        --particles=100 --particles-boxes=50,50 --particle-radius=0.5 --l0=$L --height=10 \
        --num-walls=1 --wall-positions=$L --wall-mass-factors=1000000000 \
        --piston-right-protocol-mode=step --velocity-right-piston-step=0.005 \
        --max-right-piston-travel=3.93 --auto-piston-step \
        --wall-hold-steps=12000 --steps=10 --fixed-dt=0.4 --energy-measurement \
        --eff-output=wall-ke --kbt1 --seed=$sd >> "$P/$tag/run.log" 2>&1 ) &
    while [ "$(jobs -rp | wc -l)" -ge 4 ]; do sleep 0.3; done
  done
}
run_point 38.270833 9600 100
run_point 36.3125   9600 100
run_point 37.291667 9520 80
run_point 35.3125   9520 80
wait
for t in L38p270833 L37p291667 L36p3125 L35p3125; do
  echo "$t: event logs with data $(find $P/$t -name 'ev_*.csv' -size +10k | wc -l | tr -d ' '), aborts $(grep -c ABORTING $P/$t/run.log), health $(grep -c EDMD-HEALTH $P/$t/run.log)"
done
