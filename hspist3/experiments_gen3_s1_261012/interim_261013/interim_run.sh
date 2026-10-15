#!/bin/bash
# ##CHRIS 2026-10-13 (plan-author decision 45 items 5 and 6, addendum 1): the INTERIM report's parts for one group, on the Mac, from the
# fetched folder only (read only). Run from engine-gen3's hspist3/:  bash <this> <group> <min fit failures> [<merge start, UTC>]
# Writes into experiments_gen3_s1_261012/interim_261013/<group>/ (a new folder; never into the fetched folder). For a test on other data:
# S1_FETCHED (the fetched folder), S1_OUT (the output folder), S1_CHUNKS and S1_TASKS (the chunk table and task list of s1_timing).
set -u
G=$1; K=$2; START=${3:-}
P=cluster/gen3_stage1_261012; F=${S1_FETCHED:-experiments_gen3_stage1_koa_261012}; O=${S1_OUT:-experiments_gen3_s1_261012/interim_261013}
CT=${S1_CHUNKS:-$P/chunks_stage1_koa_261012.tsv}; TL=${S1_TASKS:-$P/tasks_stage1_koa_261012.txt}
SP=/private/tmp/claude-501/-Users-chrisharing-Desktop-CCS-complex-coupled-systems-Repo-HardDisks/91cb08ec-0599-4faa-a6af-d5a7ca834255/scratchpad
[ -d "$F/$G/merged/stage1/gen3" ] || { echo "STOP: no fetched folder $F/$G"; exit 2; }
FA=$(cd "$F" && pwd); mkdir -p "$O/$G"
# (a) the verdict's INTERIM tables of this group alone: a view folder holding only a link to the fetched group
V=$SP/d45/view_$(echo "$FA" | shasum | cut -c1-8)_$G; mkdir -p "$V"; [ -e "$V/$G" ] || ln -s "$FA/$G" "$V/$G"
python3 $P/s1_verdict_261012.py --data "$V" --interim > "$O/$G/a_interim_tables.txt" 2>&1; echo "(a) exit $?"
python3 $P/s1_fitfail_261012.py "$G" --root "$FA" --min "$K" > "$O/$G/b_fit_failures.txt" 2>&1; echo "(b) exit $?"
python3 $P/s1_timing_261012.py "$G" --root "$FA" --chunks "$CT" --tasks "$TL" --out "$O/timing" > "$O/$G/c_timing.txt" 2>&1; echo "(c) exit $?"
if [ -n "$START" ]; then python3 $P/s1_merge_time_261012.py "$G" --root "$FA" --start "$START" > "$O/$G/e_merge_time.txt" 2>&1
else python3 $P/s1_merge_time_261012.py "$G" --root "$FA" > "$O/$G/e_merge_time.txt" 2>&1; fi; echo "(e) exit $?"
ls -la "$O/$G"
