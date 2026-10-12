#!/bin/bash
# ##CHRIS 2026-10-12 (plan-author decision 38 item 6; 261012 sec. 4.7.80): THE FETCH OF THE M5 SOLID-CELL JOB, the preview's fetch
# (cluster/gen3_preview_261011/fetch_preview_koa.sh) with these changes: the root on KOA (.../harddisks_gen3/m5_solid_koa), the folder
# here (hspist3/experiments_gen3_m5_solid_koa_261012), and the WHOLE root (the traces and psi6(t) files included: the rule's
# decorrelation time and the argmax nu read them; about 10 MB). Chris runs it on the Mac after the job has ended: ONE rsync over ssh
# through the DTN (one password and Duo push). Read-only on KOA: nothing there is changed or deleted. Here: a target that exists
# already is not overwritten (rsync --ignore-existing).
#   bash hspist3/cluster/gen3_m5_solid_261012/fetch_m5_solid_koa.sh
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; HS="$(dirname "$(dirname "$HERE")")"
DST="$HS/experiments_gen3_m5_solid_koa_261012"; mkdir -p "$DST"
SRC="charing@koa-dtn.its.hawaii.edu:/mnt/lustre/koa/scratch/charing/harddisks_gen3/m5_solid_koa"
rsync -av --ignore-existing "$SRC/" "$DST/"
echo "fetched into $DST: $(du -sh "$DST" | cut -f1); job folders: $(cd "$DST" && ls -d job_* 2>/dev/null | tr '\n' ' ')"
