#!/bin/bash
# ##CHRIS 2026-10-11 (plan-author decision 28, task R0; 261012 sec. 4.7.57): THE CHUNK CHECK GENERALISED for Test R from fd8d8a2's
# cluster/gen3_koa_261009/testG_chunk_logs_koa.sh. ONE change: the log names' prefix (the job name) is an argument, where fd8d8a2 has
# gen3-testG fixed; every rule and every printed line is fd8d8a2's. testR_merge.sbatch runs it before the merge plan; unit test
# test_chunk_logs_koa.py.
# fd8d8a2's description follows.
# ##CHRIS 2026-10-10 (plan-author decision 24, point 2b; 261012 sec. 4.7.52): the chunk check of the Test G merge, a script of its own
# so that it has a unit test (test_testG_chunk_logs_261010.py); merge_testG_k2.sbatch runs it before the merge plan.
# A chunk j counts as DONE when its folder <chunks>/chunk_<j> exists and EXACTLY ONE of its logs <logs>/gen3-testG_*_<j>.out has the
# line "chunk <j> done" (fd8d8a2: <logs>/gen3-testG_*_<j>.out; here <logs>/<prefix>_*_<j>.out). Every other log of a chunk is listed by name:
#   - an attempt that stopped in the start-up checks (no line "chunk <j> on ..."): its STOP line is printed (such an attempt ran nothing);
#   - an attempt that started and did not finish (a line "chunk <j> on ...", no done line): listed with its last line.
# A chunk without its folder, without a log with its done line, or with more than one such log is listed, and the check fails: the
# merge refuses an incomplete chunk set.
# usage: bash chunk_logs_koa.sh <logs dir> <log prefix: the job name> <chunks dir> <number of chunks>     exit 0: every chunk done; 2: not
set -u
G=$1; PFX=$2; C=$3; NCH=$4
nlogs=0; ndone=0; nstop=0; nstart=0; bad=""
for j in $(seq 1 "$NCH"); do
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
echo "chunk logs: $nlogs; chunks done: $ndone of $NCH (each with its folder and exactly one log with 'chunk <j> done'); earlier attempts:" \
     "$nstop stopped in the start-up checks, $nstart started and not finished"
[ -z "$bad" ] || { echo "NOT DONE:$bad"; exit 2; }
exit 0
