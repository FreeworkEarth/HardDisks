#!/bin/bash
# ##CHRIS 2026-10-12 (plan-author decision 39 item 5 f; gate G1): THE STAGE 1 CHUNK CHECK, per group. It is
# cluster/gen3_testR_261011/chunk_logs_koa.sh (c6285af7) with ONE change: the chunks checked are the group's rows of the chunk table
# (Test R: 1 .. <number of chunks>), so the last two arguments are the chunk table and the group. Every rule and printed line is
# Test R's: a chunk j counts as DONE when its folder <chunks>/chunk_<j> exists and EXACTLY ONE of its logs <logs>/<prefix>_*_<j>.out
# has the line "chunk <j> done"; every other log of a chunk is listed by name (stopped in the start-up checks: its STOP line; started
# and not finished: its last line); a chunk without its folder, without a done log or with more than one is listed, and the check fails.
# usage: bash s1_chunk_logs_koa.sh <logs dir> <log prefix: the job name> <chunks dir> <chunk table> <group>     exit 0: every chunk done; 2: not
set -u
G=$1; PFX=$2; C=$3; CT=$4; GRP=$5
JS=$(awk -F'\t' -v g="$GRP" 'NR > 1 && $5 == g {print $1}' "$CT"); NCH=$(printf '%s\n' $JS | sed '/^$/d' | wc -l | tr -d ' ')
[ "$NCH" -gt 0 ] || { echo "NOT DONE: the chunk table has no chunk of group $GRP"; exit 2; }
nlogs=0; ndone=0; nstop=0; nstart=0; bad=""
for j in $JS; do
  nd=0
  for f in "$G"/"$PFX"_*_"$j".out; do
    [ -e "$f" ] || continue
    nlogs=$((nlogs + 1))
    if grep -qx "chunk $j done" "$f"; then nd=$((nd + 1))
    elif grep -q "^chunk $j on " "$f"; then
      nstart=$((nstart + 1)); echo "earlier attempt of chunk $j, started and not finished: $(basename "$f") (its last line: $(tail -n 1 "$f"))"
    else
      nstop=$((nstop + 1)); echo "earlier attempt of chunk $j, stopped in the start-up checks (ran nothing): $(basename "$f"): $(grep -m 1 '^STOP' "$f" || echo 'no STOP line')"
    fi
  done
  if [ ! -d "$C/chunk_$j" ]; then bad="$bad $j(no folder)"
  elif [ "$nd" -eq 0 ]; then bad="$bad $j(no log with its done line)"
  elif [ "$nd" -gt 1 ]; then bad="$bad $j($nd logs with its done line)"
  else ndone=$((ndone + 1)); fi
done
echo "chunk logs: $nlogs; chunks done: $ndone of $NCH (group $GRP; each with its folder and exactly one log with 'chunk <j> done'); earlier attempts:" \
     "$nstop stopped in the start-up checks, $nstart started and not finished"
[ -z "$bad" ] || { echo "NOT DONE:$bad"; exit 2; }
exit 0
