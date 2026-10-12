#!/bin/bash
# ##CHRIS 2026-10-11 (plan-author decision 29, task P; 261012 sec. 4.7.59): READ ONLY. THIS FILE IS
# cluster/gen3_testR_261011/testR_stopped_chunks.sh (b5378f4) WITH ONE CHANGE: the sbatch file its printed resubmission line names
# (cluster/gen3_preview_261011/preview_koa.sbatch). Its description follows; for the preview read gen3-preview for gen3-testR.
# usage (from ~/harddisks_gen3_k2/hspist3, when squeue lists no gen3-preview job):
#   bash cluster/gen3_preview_261011/preview_stopped_chunks.sh logs gen3-preview <number of chunks>
# ##CHRIS 2026-10-11 (plan-author decision 28, task R0; 261012 sec. 4.7.57): READ ONLY. Test R's rule (decision 28, as decision 25): "a chunk
# that stopped in its start-up checks and ran nothing is submitted once more and listed". This script reads the chunk logs and prints
#   - every chunk without a log that has its line "chunk <j> done", classified: stopped in the start-up checks (each such log: its STOP
#     line, which names its node), started and not finished, or no log at all;
#   - THE RESUBMISSION LINE for the chunks whose every log stopped in the start-up checks, exactly once (one stopped log, no done log):
#       sbatch --array=<their numbers, as ranges> --exclude=cn-05-08-07[,<the nodes of their STOP lines>] cluster/gen3_testR_261011/testR_koa.sbatch
#     or "nothing to resubmit". It submits nothing: Chris pastes the output; the line is typed once, after the plan author's go.
#   - chunks that stopped twice or more (resubmitted once already) and chunks that started and did not finish are NOT in the line: they
#     need the plan author.
# usage (from ~/harddisks_gen3_k2/hspist3, when squeue lists no gen3-testR job):
#   bash cluster/gen3_testR_261011/testR_stopped_chunks.sh logs gen3-testR <number of chunks>
set -u
G=$1; PFX=$2; NCH=$3
stop_once=""; nodes=""; twice=""; started=""; nolog=""
for j in $(seq 1 "$NCH"); do
  nd=0; ns=0; nst=0; nl=0; nn=""
  for f in "$G"/"$PFX"_*_"$j".out; do
    [ -e "$f" ] || continue
    nl=$((nl + 1))
    if grep -qx "chunk $j done" "$f"; then nd=$((nd + 1))
    elif grep -q "^chunk $j on " "$f"; then nst=$((nst + 1))
    else
      ns=$((ns + 1)); l=$(grep -m 1 '^STOP' "$f" || echo 'no STOP line')
      echo "chunk $j: stopped in the start-up checks (ran nothing): $(basename "$f"): $l"
      n=$(printf '%s\n' "$l" | sed -n 's/^STOP: \[\([^]]*\)\].*/\1/p'); [ -n "$n" ] && nn="$nn $n"
    fi
  done
  [ "$nd" -gt 0 ] && continue
  if [ "$nl" -eq 0 ]; then nolog="$nolog $j"
  elif [ "$nst" -gt 0 ]; then started="$started $j"; echo "chunk $j: started and did not finish ($nst log(s)): needs the plan author"
  elif [ "$ns" -eq 1 ]; then stop_once="$stop_once $j"; nodes="$nodes$nn"
  else twice="$twice $j"; echo "chunk $j: stopped in the start-up checks $ns times: not submitted again (once more only): needs the plan author"; fi
done
ranges() { tr ' ' '\n' | sed '/^$/d' | sort -n | awk 'NR==1{a=$1;b=$1;next} $1==b+1{b=$1;next} {printf "%s%s,", a, (b>a?"-" b:""); a=$1; b=$1} END{if (NR) printf "%s%s", a, (b>a?"-" b:"")}'; }
echo "without a done log: stopped once $(echo $stop_once | wc -w | tr -d ' '), stopped twice or more $(echo $twice | wc -w | tr -d ' '), started and not finished $(echo $started | wc -w | tr -d ' '), no log $(echo $nolog | wc -w | tr -d ' ')$( [ -n "$nolog" ] && echo " ($(echo "$nolog" | ranges))")"
if [ -n "$stop_once" ]; then
  ex=$(printf '%s\n' cn-05-08-07 $nodes | sed '/^$/d' | awk '!seen[$0]++' | paste -sd, -)
  echo "THE RESUBMISSION LINE (type it once, after the plan author's go; from ~/harddisks_gen3_k2/hspist3):"
  echo "sbatch --array=$(echo "$stop_once" | ranges) --exclude=$ex cluster/gen3_preview_261011/preview_koa.sbatch"
else
  echo "nothing to resubmit"
fi
exit 0
