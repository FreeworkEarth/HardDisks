#!/usr/bin/env bash
# ##CHRIS 2026-10-07 (261012 sec. 4.4.10, gate version 2): copy the gate-v2 outputs back FROM KOA (run on the Mac, repo root).
# Root $SCRATCH/harddisks_resched2 (the gate-v2 clone), into hspist3/experiments_resched_gate2_261007/, a tree of its own.
# Test T: summaries only (red_nu.csv, acf_runs.npz, run.log, .build_git, failed-run stdout.log); the E0/E2 audit runs: the
# report and the run logs with their audit lines (traces and event logs stay on KOA); the profile: logs and times.
# Nothing on either side is deleted.   bash hspist3/cluster/resched_gate_261005/fetch_resched2.sh
KOA_USER=charing; SCRATCH="/mnt/lustre/koa/scratch/charing"; DTN=$KOA_USER@koa-dtn.its.hawaii.edu; R=$SCRATCH/harddisks_resched2
LOC=hspist3/experiments_resched_gate2_261007
SUM=(--prune-empty-dirs --include='*/' --include='red_nu.csv' --include='acf_runs.npz' --include='run.log' --include='.build_git'
     --include='stdout.log' --exclude='*')
rsync -av "${SUM[@]}" "$DTN:$R/hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/resched_testT_261007/" \
          "$LOC/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/resched_testT_261007/"
rsync -av "$DTN:$R/hspist3/.build_generation" "$LOC/"
rsync -av --prune-empty-dirs --include='*/' --include='report.txt' --include='run.log' --include='version.txt' --include='command.txt' \
          --include='red_9700.csv' --exclude='*' "$DTN:$R/resched_gate_261005/" "$LOC/resched_gate_261005/"
rsync -av --exclude='ev.csv' --exclude='tr.csv' --exclude='perf.data' "$DTN:$R/profile_edmd_*" "$LOC/"
# ##CHRIS 2026-10-07 (261012 sec. 4.4.13, gate v3): Test T-prime (summaries, as Test T) and the ASan job (report, build log, run logs;
# not the sanitizer binary, not the traces). Before those jobs have written anything, these two lines report "No such file".
rsync -av "${SUM[@]}" "$DTN:$R/hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/resched_testTprime_261007/" \
          "$LOC/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/resched_testTprime_261007/"
rsync -av --exclude='00ALLINONE_asan' --exclude='*.csv' --exclude='*.json' "$DTN:$R/asan_261007_*" "$LOC/"
