#!/bin/bash
# ##CHRIS 2026-10-09 (261012 sec. 4.7.18; stage E of the plan-author programme of sec. 4.7.12): the M5 evidence with the binaries of
# the committed tree 78ff48d (build_clean.sh, plus the long-double builds of the same archive; 78ff48d = c685d3b + the run header's time quantum):
#   E2a  the default build stays byte-identical on the M1 and M2 harness outputs (gen3_m1 audit, gen3_m2 audit --quick, gen3_m2 audit)
#   E2b  the long-double builds (-DEDMD3_LONG_DOUBLE) of the harnesses: on this Mac long double is the 8-byte double, so their outputs
#        must equal the default build's byte for byte (a consistency check of the transformation; the numerical check is KOA's)
#   E2c  the long-double driver against the default driver on two replays of harness cells (event hash and events) -- the same reason
#   rule 4 (the Makefile and 00ALLINONE.c changed): gen2 byte identity with 7b08827 (ctrl_min, ctrl_leg; default and --engine=gen2) with 78ff48d's default
# usage (from hspist3/ of the engine-gen3 worktree): bash experiments_gen3_m5_261009/stageE_evidence.sh <binary dir> <out dir>
set -u
B=$(cd "$1" && pwd); mkdir -p "$2"; O=$(cd "$2" && pwd)   # absolute: the runs cd into their folders (the first run of 13:19 used
                                                       # relative paths: E2a, E2b and the --engine=gen2 E0 did not start; kept as *_failed_relpath)
HERE="$(cd "$(dirname "$0")" && pwd)"; HS="$(dirname "$HERE")"
GATE=/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3/cluster/resched_gate_261005
echo "binaries $B; out $O; $(date '+%Y-%m-%d %H:%M:%S %Z')"; df -h "$O" | tail -1
( cd "$O" && mkdir -p e2a e2b e2c )
# E2a and E2b: 2 processes (one per build), each running its three harness modes in order
[ -s "$O/e2a/m2_audit_output.txt" ] || ( cd "$O/e2a" && "$B/gen3_m1" audit > m1_audit_output.txt && "$B/gen3_m2" audit --quick > m2_audit_quick_output.txt && "$B/gen3_m2" audit > m2_audit_output.txt; echo "e2a exit $?" ) &
[ -s "$O/e2b/m2_audit_output.txt" ] || ( cd "$O/e2b" && "$B/gen3_m1_ld" audit > m1_audit_output.txt && "$B/gen3_m2_ld" audit --quick > m2_audit_quick_output.txt && "$B/gen3_m2_ld" audit > m2_audit_output.txt; echo "e2b exit $?" ) &
# rule 4: E0 with c685d3b's default build (2 x 4 runs, 2 at a time)
printf '#!/bin/bash\nexec %s "$@" --engine=gen2\n' "$B/00ALLINONE" > "$B/bin_engine_gen2.sh"; chmod +x "$B/bin_engine_gen2.sh"
for tag in default gen2flag; do
  b="$B/00ALLINONE"; [ "$tag" = gen2flag ] && b="$B/bin_engine_gen2.sh"
  [ -d "$O/e0/$tag" ] || python3 "$GATE/audit_runs_261007.py" run --bin "$b" --out "$O/e0/$tag" --cases ctrl_min,ctrl_leg --jobs 2 > "$O/e0_run_$tag.txt" 2>&1
done
wait
# E2c: the long-double driver against the default driver on two harness cells (the replay states of stage A)
S=$HS/experiments_gen3_m3_261009/evidence_f42befb/replay/states
[ -s "$O/e2c/replays.txt" ] || for f in "$S/m1_fluid.g3state" "$S/m2_spring_pi8.g3state"; do
  for tag in default ld; do
    b="$B/00ALLINONE"; [ "$tag" = ld ] && b="$B/00ALLINONE_ld"
    "$b" --engine-replay="$f" --replay-read-every=1 --replay-substeps=1 2>&1 | grep "^\[EDMD3-REPLAY\]\|^\[EDMD3-HEALTH\]" | sed "s/^/$tag $(basename "$f" .g3state): /" | cut -c1-260
  done
done > "$O/e2c/replays.txt"
echo "done $(date '+%H:%M:%S')"
