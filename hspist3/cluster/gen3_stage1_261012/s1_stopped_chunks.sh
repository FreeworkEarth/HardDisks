#!/bin/bash
# ##CHRIS 2026-10-12 (plan-author decision 39 item 5 f; gate G1): READ ONLY. The stage 1 copy of cluster/gen3_testR_261011/
# testR_stopped_chunks.sh (c6285af7), Test R's rule (decision 28, as decision 25): "a chunk that stopped in its start-up checks and ran
# nothing is submitted once more and listed". TWO changes: (1) the chunks read are the group's rows of the chunk table; (2) the
# resubmission line is printed per time class of those chunks (the chunk table's class_min; one array has one time limit), with the
# group's job name. It prints every chunk of the group without a done log, classified (stopped in the start-up checks, with its STOP line
# and node; started and not finished; no log), and THE RESUBMISSION LINE(S) for the chunks whose every log stopped exactly once:
#   sbatch --job-name=gen3-s1-<group> --time=<class> --array=<their numbers> --exclude=cn-05-08-07[,<their nodes>] cluster/gen3_stage1_261012/s1_koa.sbatch
# or "nothing to resubmit". It submits nothing; chunks that stopped twice or more, and chunks that started and did not finish, need the
# plan author.
# usage (from ~/harddisks_gen3_k2/hspist3, when squeue lists no gen3-s1-<group> job):
#   bash cluster/gen3_stage1_261012/s1_stopped_chunks.sh logs gen3-s1-<group> cluster/gen3_stage1_261012/chunks_stage1_koa_261012.tsv <group>
set -u
G=$1; PFX=$2; CT=$3; GRP=$4
JS=$(awk -F'\t' -v g="$GRP" 'NR > 1 && $5 == g {print $1}' "$CT")
stop_once=""; nodes=""; twice=""; started=""; nolog=""
for j in $JS; do
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
echo "group $GRP without a done log: stopped once $(echo $stop_once | wc -w | tr -d ' '), stopped twice or more $(echo $twice | wc -w | tr -d ' '), started and not finished $(echo $started | wc -w | tr -d ' '), no log $(echo $nolog | wc -w | tr -d ' ')"
if [ -n "$stop_once" ]; then
  ex=$(printf '%s\n' cn-05-08-07 $nodes | sed '/^$/d' | awk '!seen[$0]++' | paste -sd, -)
  echo "THE RESUBMISSION LINE(S) (type each once, after the plan author's go; from ~/harddisks_gen3_k2/hspist3):"
  for cls in $(for j in $stop_once; do awk -F'\t' -v j="$j" 'NR > 1 && $1 == j {print $10}' "$CT"; done | sort -n | uniq); do
    js=$(for j in $stop_once; do awk -F'\t' -v j="$j" -v c="$cls" 'NR > 1 && $1 == j && $10 == c {print $1}' "$CT"; done | tr '\n' ' ')
    echo "sbatch --job-name=gen3-s1-$GRP --time=$((cls / 60)):$(printf '%02d' $((cls % 60))):00 --array=$(echo "$js" | ranges) --exclude=$ex cluster/gen3_stage1_261012/s1_koa.sbatch"
  done
else
  echo "nothing to resubmit"
fi
exit 0
