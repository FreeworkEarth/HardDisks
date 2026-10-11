#!/bin/bash
# decision 24 point 2a, printed by the code: a synthetic KOA-scale case -- 144 chunks; chunk 1 from K2-3 (15043031); chunks 2-144 from
# array 15043056, of which 57-61, 66-70, 110-114, 132-136 stopped in a start-up check (an env line and a STOP line, no folder made);
# those 20 resubmitted as array 15050000 and complete (164 logs). Then (i) K2-4 block A as written in Report K2b (dc158b5's runsheet),
# with HOME pointing to the synthetic home (squeue/sacct are not run: they need Slurm), and (ii) dc158b5's merge-job chunk loop, verbatim.
set -u
B=$1; H=$B/home; G=$H/harddisks_gen3/hspist3/logs; C=$B/lustre/testG_koa/chunks
mkdir -p "$G" "$C"
ENV="env: gcc (GCC) 14.3.0 | Python 3.12.14 at /home/charing/envs/hd/bin/python3 | git version 2.52.0"
FAILED="57 58 59 60 61 66 67 68 69 70 110 111 112 113 114 132 133 134 135 136"
done_log() { printf '%s\nchunk %s on cn-03-38-07 (Intel(R) Xeon(R) CPU E5-2670 v3 @ 2.30GHz): lines %s-%s; 00ALLINONE  git a643bc7  target koa\nchunk %s done\n' "$ENV" "$1" $(( ($1-1)*25+1 )) $(( $1*25 )) "$1"; }
done_log 1 > "$G/gen3-testG_15043031_1.out"; mkdir -p "$C/chunk_1"
for j in $(seq 2 144); do
  case " $FAILED " in
    *" $j "*) printf '%s\nSTOP: ./00ALLINONE is not the recorded clean koa build of HEAD a643bc7\n' "$ENV" > "$G/gen3-testG_15043056_$j.out"
              done_log "$j" > "$G/gen3-testG_15050000_$j.out" ;;
    *) done_log "$j" > "$G/gen3-testG_15043056_$j.out" ;;
  esac
  mkdir -p "$C/chunk_$j"
done
echo "== (i) block A of Report K2b, its four read-only file lines (HOME = the synthetic home):"
export HOME=$H
ls ~/harddisks_gen3/hspist3/logs/gen3-testG_*.out | wc -l
grep -lx 'chunk [0-9]* done' ~/harddisks_gen3/hspist3/logs/gen3-testG_*.out | wc -l
cat ~/harddisks_gen3/hspist3/logs/gen3-testG_*.out | grep -c -e FAILED -e REFUSED -e STOP
grep -h '^chunk [0-9]* on ' ~/harddisks_gen3/hspist3/logs/gen3-testG_*.out | sed 's/: lines.*//; s/^chunk [0-9]* on [^ ]* //' | sort | uniq -c
echo "== (ii) dc158b5's merge-job chunk loop, verbatim (ROOT's chunks = the synthetic ones):"
ROOT=$B/lustre/testG_koa; NCH=144
miss=""; multi=""
for j in $(seq 1 "$NCH"); do
  n=0; ok=0
  for f in "$G"/gen3-testG_*_"$j".out; do [ -e "$f" ] || continue; n=$((n + 1)); grep -qx "chunk $j done" "$f" && ok=1; done
  { [ -d "$ROOT/chunks/chunk_$j" ] && [ "$ok" = 1 ]; } || miss="$miss $j"
  [ "$n" -le 1 ] || multi="$multi $j"
done
[ -z "$miss" ] && echo "(no STOP: every chunk has its folder and a log with 'chunk <j> done')" || echo "STOP: chunks without their folder or without 'chunk <j> done' in ...:$miss"
[ -z "$multi" ] || echo "note: chunks with more than one log:$multi"
NF=$(cat "$G"/gen3-testG_*.out | grep -c FAILED); echo "FAILED lines in the chunk logs: $NF"
echo "== (iii) the new check (testG_chunk_logs_koa.sh) on the same logs:"
bash "$2" "$G" "$C" 144 | tail -3; echo "exit ${PIPESTATUS[0]}"
