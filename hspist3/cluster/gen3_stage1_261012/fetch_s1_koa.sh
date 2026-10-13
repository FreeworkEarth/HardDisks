#!/bin/bash
# ##CHRIS 2026-10-12 (plan-author decision 39 item 5 f; registration 261012 sec. 4.8.8, "Fetch. Summaries, run records and logs for
# everything; full raw data for the (N = 100, pi/8) cell and for all snapshot files"; gate G1): THE FETCH OF ONE STAGE 1 GROUP. It is
# cluster/gen3_preview_261011/fetch_preview_koa.sh with these changes: the root on KOA (.../harddisks_gen3/stage1/<group>), the folder
# here (hspist3/experiments_gen3_stage1_koa_261012/<group>), and the filter: the (N = 100, pi/8) cell B_N0100_eta0.392699 whole (its
# traces, psi6(t) and virial files too), every snapshot file, and otherwise no trace, psi6(t) or virial file and no run folder. Chris runs
# it on the Mac after the group's merge job has ended: ONE rsync over ssh through the DTN (one password and Duo push). Read-only on KOA.
# Here: a target that exists already is not overwritten (rsync --ignore-existing).
#   bash hspist3/cluster/gen3_stage1_261012/fetch_s1_koa.sh <group: B100 | B400 | D | B900 | S | B1600>
set -uo pipefail
GRP=${1:?the group}; case "$GRP" in B100|B400|D|B900|S|B1600) ;; *) echo "not a stage 1 group: $GRP"; exit 2 ;; esac
HERE="$(cd "$(dirname "$0")" && pwd)"; HS="$(dirname "$(dirname "$HERE")")"
DST="$HS/experiments_gen3_stage1_koa_261012/$GRP"; mkdir -p "$DST"
SRC="charing@koa-dtn.its.hawaii.edu:/mnt/lustre/koa/scratch/charing/harddisks_gen3/stage1/$GRP"
rsync -av --ignore-existing --include='/merged/stage1/gen3/B_N0100_eta0.392699/***' --exclude='/chunks/' --exclude='wall_x_positions_*' \
      --exclude='psi6_t_*' --exclude='virial_blocks_*' --exclude='.done_runs/' "$SRC/" "$DST/"
echo "fetched into $DST: $(du -sh "$DST" | cut -f1); merge record: $(head -1 "$DST/merged/.merge_complete" 2>/dev/null || echo MISSING)"
