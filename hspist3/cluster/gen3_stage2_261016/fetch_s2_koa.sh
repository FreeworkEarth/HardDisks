#!/bin/bash
# ##CHRIS 2026-10-16 (plan-author decision 58; notes 261012 sec. 4.7.142): THE FETCH OF STAGE 2 "S2" (runsheet paste E). It is fetch_p2_koa.sh
# (P2's, with the openrsync-tested filter and the gzipped traces, Reports D57-1 and D57-2) with THESE CHANGES, all of them: the root on KOA
# .../harddisks_gen3/stage2/S2 and the folder here hspist3/experiments_gen3_stage2_koa_261016/S2 (decision 58). Chris runs it on the Mac after
# the merge job has ended: ONE rsync over ssh through the DTN (one password and Duo push). Read-only on KOA; here a file that exists already is
# not overwritten (rsync --ignore-existing). The filter: the merged tree (its tables, run logs, psi6(t) and virial files, the gzipped traces,
# the run folders' small files with the initial states), no snapshot file, no chunk folder, no nested experiments_speed_of_sound/ folder.
#   bash hspist3/cluster/gen3_stage2_261016/fetch_s2_koa.sh
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; HS="$(dirname "$(dirname "$HERE")")"
DST="$HS/experiments_gen3_stage2_koa_261016/S2"; mkdir -p "$DST"
SRC="charing@koa-dtn.its.hawaii.edu:/mnt/lustre/koa/scratch/charing/harddisks_gen3/stage2/S2"
rsync -av --ignore-existing --exclude='/chunks/' --include='wall_x_positions_*.csv.gz' --exclude='wall_x_positions_*' --exclude='snapshots_*' \
      --exclude='experiments_speed_of_sound/' \
      "$SRC/" "$DST/"
echo "fetched into $DST: $(du -sh "$DST" | cut -f1); merge record: $(head -1 "$DST/merged/.merge_complete" 2>/dev/null || echo MISSING); p2_reduce: $(tail -1 "$DST/merged/p2_reduce_summary.txt" 2>/dev/null || echo MISSING)"
