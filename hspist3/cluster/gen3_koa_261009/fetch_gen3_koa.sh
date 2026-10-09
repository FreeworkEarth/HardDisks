#!/bin/bash
# ##CHRIS 2026-10-09 (261012 sec. 4.7.20; stage G): PREPARED, NOT RUN. Fetch the KOA gen3 gate outputs to this Mac -- Chris runs it,
# from the Mac, after the jobs have ended; it asks for the KOA password and Duo push (rsync over ssh through the DTN). It copies:
#   the merged Test G tree WITHOUT the compressed traces and psi6 files (the verdict reads red_nu.csv, acf_runs.npz, run.log, .build_git,
#   .failed_*, and the A-fixed red_/run_ files) -> hspist3/experiments_gen3_gate_koa_261009/testG/
#   the one-node gate job's output folder and the cross-node job's output folder -> hspist3/experiments_gen3_gate_koa_261009/{gate,xnode}/
# Read-only on KOA: nothing there is changed or deleted. Here: a target that exists already is not overwritten (rsync --ignore-existing).
#   bash hspist3/cluster/gen3_koa_261009/fetch_gen3_koa.sh <gate job id> <xnode job id>
set -uo pipefail
GATE=${1:?gate job id}; XN=${2:?xnode job id}
HERE="$(cd "$(dirname "$0")" && pwd)"; HS="$(dirname "$(dirname "$HERE")")"
DST="$HS/experiments_gen3_gate_koa_261009"; mkdir -p "$DST"
SRC="charing@koa-dtn.its.hawaii.edu:/mnt/lustre/koa/scratch/charing/harddisks_gen3"
rsync -av --ignore-existing --exclude='*.csv.gz' --exclude='.done_runs/' "$SRC/testG_koa/merged/testG/" "$DST/testG/"
rsync -av --ignore-existing --exclude='*.csv.gz' "$SRC/gen3_gate_koa_${GATE}/" "$DST/gate/"
rsync -av --ignore-existing --exclude='*.csv.gz' "$SRC/gen3_xnode_${XN}/" "$DST/xnode/"
echo "fetched into $DST: $(du -sh "$DST" | cut -f1)"
