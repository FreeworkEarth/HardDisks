#!/bin/bash
# ##CHRIS 2026-10-12 (plan-author decision 39 item 5 c): THE FETCH OF THE STAGE 1 GATE (the gate job's root s1_gate and the cross-node root
# s1_xnode), as fetch_preview_koa.sh: ONE rsync over ssh through the DTN per root, in one command line (two passwords and Duo pushes);
# read-only on KOA; here a target that exists is not overwritten. Without the traces and psi6(t) files (the compares ran on KOA); with the
# job logs, the build-time record, the E0 reports, the harness outputs and every run.log; without E0's afix event logs (ev_*.csv, 17 MB
# each, six per build: E0 compares them on KOA). Into hspist3/experiments_gen3_s1_gate_koa_261012/.
#   bash hspist3/cluster/gen3_s1_gate_261012/fetch_s1_gate_koa.sh
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; HS="$(dirname "$(dirname "$HERE")")"
DST="$HS/experiments_gen3_s1_gate_koa_261012"; mkdir -p "$DST"
SRC="charing@koa-dtn.its.hawaii.edu:/mnt/lustre/koa/scratch/charing/harddisks_gen3"
for r in s1_gate s1_xnode; do
  rsync -av --ignore-existing --exclude='wall_x_positions_*' --exclude='psi6_t_*' --exclude='snapshots_*' --exclude='ev_*.csv' --exclude='.done_runs/' "$SRC/$r/" "$DST/$r/"
done
echo "fetched into $DST: $(du -sh "$DST" | cut -f1); build record: $(sed -n 's/^commit //p' "$DST/s1_gate/BUILD_RECORD_S1.txt" 2>/dev/null || echo MISSING)"
