#!/bin/bash
# ##CHRIS 2026-10-08 (261012 sec. 4.7.4, decision 2): compare the "before" run (code without the loader guard) with the "after" run
# (with it) of run_set.sh: exit code, stdout and stderr per script, and the SHA-256 of every file the scripts wrote.
# usage: compare_before_after.sh <before dir> <after dir>
B=$1; A=$2; REPO=/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks
printf "%-14s %-6s %-7s %-7s %s\n" script rc stdout stderr "stdout lines"
for s in $(awk '$1 != "DONE" {print $1}' "$B/times.txt"); do
  o=$(cmp -s "$B/$s.out" "$A/$s.out" && echo same || echo DIFF); e=$(cmp -s "$B/$s.err" "$A/$s.err" && echo same || echo DIFF)
  printf "%-14s %-6s %-7s %-7s %s\n" "$s" "$(cat "$B/$s.rc")/$(cat "$A/$s.rc")" "$o" "$e" "$(wc -l < "$A/$s.out" | tr -d ' ')"
done
norm() { sed "s#  $REPO/#  #" "$1" | grep -v "  hspist3/validation/" | sort -k2; }
echo; echo "files written: before $(norm "$B/written.sha" | wc -l | tr -d ' '), after $(norm "$A/written.sha" | wc -l | tr -d ' ')"
if diff <(norm "$B/written.sha") <(norm "$A/written.sha") > /dev/null; then echo "written files: all hashes identical"
else echo "written files: DIFFERENT"; diff <(norm "$B/written.sha") <(norm "$A/written.sha"); fi
