#!/bin/bash
# ##CHRIS 2026-10-15 (plan-author decision 57; notes 261012 sec. 4.7.135): THE FETCH OF THE STAGE 2 PILOT P2 (runsheet block e). It is
# cluster/gen3_stage1_261012/fetch_s1_koa.sh with these changes: the root on KOA (.../harddisks_gen3/stage2/P2), the folder here
# (hspist3/experiments_gen3_stage2_koa_261016/P2, decision 57), and the filter: the merged tree with its tables (the stage 1 reductions,
# p2_runs.csv, p2_frames.csv), the run logs, the psi6(t) and virial files, and the run folders' small files (the initial state, stdout,
# the run parameters; not their nested experiments_speed_of_sound/ folders); no trace and no snapshot file (p2_reduce.py read them on KOA;
# they stay there), no chunk folder. Chris runs it on the Mac after the merge job has ended: ONE rsync over ssh through the DTN (one
# password and Duo push).
# ##CHRIS 2026-10-15 (the Mac rehearsal, Report D57-1): the Mac's only rsync is /usr/bin/rsync = openrsync (protocol 29), which did not
# apply the first filter's '.done_runs/run*/**' exclude; the filter is now four plain excludes, tested with openrsync on the rehearsal's
# merged root (0.4 MB per run, no trace, no snapshot).
# ##CHRIS 2026-10-15 (decision 57 addendum 1, task 2: the P2 report's figure "traces and fits for one seed per cell"): the merged, gzipped
# traces (wall_x_positions_*.csv.gz, about 70 KB per run, about 28 MB for the 400) are now fetched too; an uncompressed trace is still
# left on KOA. Tested with openrsync on the rehearsal's merged root (experiments_gen3_p2_261016/rehearsal_261015c/fetch_filter_test2_openrsync.txt).
# Read-only on KOA. Here: a file that exists already is not overwritten (rsync --ignore-existing).
#   bash hspist3/cluster/gen3_stage2_261016/fetch_p2_koa.sh
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; HS="$(dirname "$(dirname "$HERE")")"
DST="$HS/experiments_gen3_stage2_koa_261016/P2"; mkdir -p "$DST"
SRC="charing@koa-dtn.its.hawaii.edu:/mnt/lustre/koa/scratch/charing/harddisks_gen3/stage2/P2"
rsync -av --ignore-existing --exclude='/chunks/' --include='wall_x_positions_*.csv.gz' --exclude='wall_x_positions_*' --exclude='snapshots_*' \
      --exclude='experiments_speed_of_sound/' \
      "$SRC/" "$DST/"
echo "fetched into $DST: $(du -sh "$DST" | cut -f1); merge record: $(head -1 "$DST/merged/.merge_complete" 2>/dev/null || echo MISSING); p2_reduce: $(tail -1 "$DST/merged/p2_reduce_summary.txt" 2>/dev/null || echo MISSING)"
