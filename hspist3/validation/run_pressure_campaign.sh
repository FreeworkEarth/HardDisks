#!/usr/bin/env bash
# ##CHRIS: full equilibrium pressure validation campaign.
# Three estimators (pair virial, wall momentum flux x and y), multiple seeds,
# block-resolved uncertainties, across eta and N. See pressure_validation.c.
set -uo pipefail
D="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$D/.."
OUT="${1:?usage: run_pressure_campaign.sh OUTDIR [JOBS]}"; JOBS="${2:-12}"
mkdir -p "$OUT"
TRAJ="$OUT/pressure_trajectories.csv"; BLK="$OUT/pressure_blocks.csv"
echo "eta,N,seed,boxW,boxH,equilibration_time,measurement_time,pair_event_count,wall_L_events,wall_R_events,wall_B_events,wall_T_events,Z_pair,Z_pair_sem,Z_wall_x,Z_wall_x_sem,Z_wall_y,Z_wall_y_sem,T_mean,psi6_global_mean,psi6_local_mean,n_blocks,forced_advance,overlap_repair,clamp_repair,wall_overdue,valid" > "$TRAJ"
echo "eta,N,seed,block_index,t_start,t_end,T,Z_pair,Z_wall_x,Z_wall_y,psi6" > "$BLK"

ETAS=(0.005 0.02 0.05 0.10 0.20 0.30 0.40 0.50 0.60 0.65 0.67 0.69 0.698 0.702 0.710 0.718 0.720)
seeds_for(){ case "$1" in 400) echo 5;; 900) echo 4;; *) echo 3;; esac; }
# Denser states relax more slowly and their blocks are more correlated, so give
# them a longer equilibration and longer blocks rather than assuming one window
# suits every density.
equil_for(){ awk -v e="$1" 'BEGIN{print (e>0.60)? 400 : (e>0.30? 200 : 100)}'; }
bdt_for(){   awk -v e="$1" 'BEGIN{print (e>0.60)? 20 : (e>0.30? 12 : 8)}'; }
NBLOCKS=30

i=0
for eta in "${ETAS[@]}"; do
  for N in 400 900 1600; do
    ns=$(seeds_for $N); eq=$(equil_for $eta); bd=$(bdt_for $eta)
    for ((k=0;k<ns;k++)); do
      seed=$(( 20260907 + 104729*k + 7919*i ))
      tj="$OUT/.part_traj_${eta}_${N}_${k}.csv"; bk="$OUT/.part_blk_${eta}_${N}_${k}.csv"
      ( ./validation/pressure_validation "$eta" "$N" "$seed" "$NBLOCKS" "$bd" "$eq" "$tj" "$bk" \
          >> "$OUT/run.log" 2>&1 ) &
      i=$((i+1))
      while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do sleep 1; done
    done
  done
done
wait
cat "$OUT"/.part_traj_*.csv >> "$TRAJ" 2>/dev/null
cat "$OUT"/.part_blk_*.csv  >> "$BLK"  2>/dev/null
rm -f "$OUT"/.part_traj_*.csv "$OUT"/.part_blk_*.csv
echo "trajectories: $(( $(wc -l < "$TRAJ") - 1 ))   invalid: $(awk -F, 'NR>1 && $26==0' "$TRAJ" | wc -l | tr -d ' ')"
echo "blocks      : $(( $(wc -l < "$BLK") - 1 ))"
echo "Complete: $OUT"
