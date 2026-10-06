#!/usr/bin/env bash
# ##CHRIS 2026-10-05 (261012 sec. 4.4): copy the engine-gate outputs back FROM KOA (run on the Mac, from the repo root).
# Into hspist3/experiments_resched_gate_261005/, a tree of its own: the 279282b data and its registered analyses never see
# these files (validation/resched_gate_261005.py reads them with an explicit flag). Nothing on either side is deleted.
#   bash hspist3/cluster/resched_gate_261005/fetch_resched.sh
KOA_USER=charing; SCRATCH="/mnt/lustre/koa/scratch/charing"; DTN=$KOA_USER@koa-dtn.its.hawaii.edu; R=$SCRATCH/harddisks_resched
LOC=hspist3/experiments_resched_gate_261005
SUM=(--prune-empty-dirs --include='*/' --include='red_*.csv' --include='red_nu.csv' --include='acf_runs.npz'
     --include='run.log' --include='run_*.log' --include='summary_*.csv' --include='command*.txt' --include='.build_git' --exclude='*')
rsync -av "${SUM[@]}" "$DTN:$R/hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/confinement_B_20261013/" \
          "$LOC/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/confinement_B_20261013/"
rsync -av "${SUM[@]}" "$DTN:$R/hspist3/experiments_energy_transfer/paper1_confinement_Afix_261004/" \
          "$LOC/experiments_energy_transfer/paper1_confinement_Afix_261004/"
rsync -av "$DTN:$R/hspist3/.build_generation" "$LOC/"
# G-E2 (every file, about 20 MB) and the profile directories (logs and times; the event logs stay on KOA)
rsync -av "$DTN:$R/resched_gate_261005/" "$LOC/resched_gate_261005/"
rsync -av --exclude='ev.csv' --exclude='tr.csv' --exclude='perf.data' "$DTN:$R/profile_edmd_*" "$LOC/"
