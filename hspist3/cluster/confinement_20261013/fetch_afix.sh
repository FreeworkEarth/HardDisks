#!/usr/bin/env bash
# ##CHRIS 2026-10-04 (Task Y2): copy the A-fixed campaign back FROM KOA (run on the Mac, from the repo root). Summaries of every
# cell, plus the full A-fixed pilot cell (its event logs feed gate G1). Nothing on either side is deleted.
KOA_USER=charing; SCRATCH="/mnt/lustre/koa/scratch/charing"; DTN=$KOA_USER@koa-dtn.its.hawaii.edu; R=$SCRATCH/harddisks/hspist3
case "${1:-all}" in
  pilot) rsync -av "$DTN:$R/experiments_energy_transfer/paper1_confinement_Afix_261004/pilot_epi8_H_H10_L10/" "hspist3/experiments_energy_transfer/paper1_confinement_Afix_261004/pilot_epi8_H_H10_L10/" ;;
  *) rsync -av --prune-empty-dirs --include='*/' --include='red_*.csv' --include='run_*.log' --include='summary_*.csv' \
           --exclude='*' "$DTN:$R/experiments_energy_transfer/paper1_confinement_Afix_261004/" "hspist3/experiments_energy_transfer/paper1_confinement_Afix_261004/" ;;
esac
