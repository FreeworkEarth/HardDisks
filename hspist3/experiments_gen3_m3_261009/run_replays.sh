#!/bin/bash
# ##CHRIS 2026-10-09 (M3, 261012 sec. 4.7.14, acceptances 2 and 3): every harness cell (3 of M1, 13 of M2), dumped by the
# harnesses (gen3_m1 dump, gen3_m2 dump) and replayed THROUGH THE DRIVER (00ALLINONE --engine-replay) in four variants:
#   a  readers at every target                        (acceptance 2: the harness's event hash from the same state)
#   b  no readers                                     (acceptance 3: outputs off)
#   c  readers at every target, 3 stops per interval  (acceptance 3: another output cadence and extra advance stops)
#   d  readers every 7th target, 2 stops per interval (acceptance 3: a third cadence)
# usage (from hspist3/): bash experiments_gen3_m3_261009/run_replays.sh <binary> <harness dir> <state dir> > replay_output.txt
set -u
BIN=$1; HB=$2; D=$3
mkdir -p "$D"
echo "## harness dumps"
( cd "$D" && "$HB/gen3_m1" dump "$D" | grep -E "targets" | sed "s#$D/##" && "$HB/gen3_m2" dump "$D" | grep -E "targets" | sed "s#$D/##" )
echo; echo "## replays through the driver"
ok=0; n=0
for f in "$D"/m1_*.g3state "$D"/m2_*.g3state; do
  for v in a b c d; do
    case $v in a) o="--replay-read-every=1 --replay-substeps=1";; b) o="--replay-read-every=0 --replay-substeps=1";;
               c) o="--replay-read-every=1 --replay-substeps=3";; d) o="--replay-read-every=7 --replay-substeps=2";; esac
    out=$("$BIN" --engine-replay="$f" $o 2>&1); rc=$?
    line=$(echo "$out" | grep "^\[EDMD3-REPLAY\]"); h=$(echo "$out" | grep "^\[EDMD3-HEALTH\]" | sed -E 's/.*: (clean=[01]) .*/\1/')
    echo "$v $(basename "$f" .g3state) rc=$rc $h ${line#\[EDMD3-REPLAY\] }"
    n=$((n+1)); [ $rc -eq 0 ] && echo "$line" | grep -q ": MATCH;" && echo "$h" | grep -q "clean=1" && ok=$((ok+1))
  done
done
echo; echo "replays: $n, MATCH with clean=1: $ok"
