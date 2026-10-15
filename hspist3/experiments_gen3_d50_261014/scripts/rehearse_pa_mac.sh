#!/bin/bash
# ##CHRIS 2026-10-14 (plan-author decision 50 item 3 d; 261012 sec. 4.7.113): THE MAC REHEARSAL OF THE PART A WORKER (option ii), with the
# clean build of the final code. Into a NEW data root <out>/data:
#   R1  40 short tasks of pa_tasks.py (N 100 and 400; eta 0.05, pi/8, 0.60, 0.69; the five positions; one seed; T_eq 600, window 2400,
#       tail 60 steps, a trace row every 6 steps; the event log for s = 0, j = 0, the check subset), 6 workers at a time;
#   R2  the same 40 task lines again: every one must be skipped ("done before"), no new file;
#   R3  one task line with a different binary into a folder written by the build: REFUSED by the build guard (exit 3);
#   R4  one task line whose divider lies outside the box: the binary refuses it, the worker reports FAILED and keeps .failed_<seed>_<time>;
#   R5  two task lines of the draft pilot at the registered lengths (T_eq + 5000 sigma-time + 20 sigma-time): N 100 eta 0.05 j = -2
#       (T_eq 1000) and N 400 eta 0.69 j = +1 (T_eq 5000), timed, into their own data root <out>/data_full (their seeds are those of
#       the short tasks: a seed depends on the cell, the position and s, not on the lengths, so in one root they read "done before");
# then every red_<seed>.csv of R1 and R5 is listed (summary.txt).
#   usage (from hspist3/): bash experiments_gen3_d50_261014/scripts/rehearse_pa_mac.sh <new folder>
set -u
OUT=${1:?a new folder}; [ -e "$OUT" ] && { echo "REFUSED: $OUT exists"; exit 1; }
HS=$(pwd); W="$HS/cluster/gen3_partA_261014/pa_worker_koa.sh"
case "$OUT" in /*) ;; *) OUT="$HS/$OUT" ;; esac   # the worker takes absolute data roots only (its own guard refused a relative one)
BIN="$HS/experiments_gen3_d50_261014/bin_8f538a2/00ALLINONE"; OTHER="$HS/experiments_gen3_d48_261013/bin_d2dc177/00ALLINONE"
mkdir -p "$OUT/data"; export HD_DATA="$OUT/data" SDL_VIDEODRIVER=dummy
python3 cluster/gen3_partA_261014/pa_tasks.py --sizes 100,400 --seeds 1 --etas 0.05,0.392699,0.60,0.69 --log s0j0 \
        --steps 600,2400,60,6 > "$OUT/tasks_short.txt" 2> "$OUT/tasks_short.stderr"
python3 cluster/gen3_partA_261014/pa_tasks.py --sizes 100,400 --seeds 1 > "$OUT/tasks_full.txt" 2> "$OUT/tasks_full.stderr"
echo "== R1: $(wc -l < "$OUT/tasks_short.txt" | tr -d ' ') short tasks, 6 at a time"
t0=$(date +%s)
tr '\n' '\0' < "$OUT/tasks_short.txt" | xargs -0 -n 1 -P 6 sh -c 'HD_BIN="$0" bash "$1" $2' "$BIN" "$W" > "$OUT/R1_worker_lines.txt" 2>&1
echo "   $(grep -c ' ok (' "$OUT/R1_worker_lines.txt") ok, $(grep -c FAILED "$OUT/R1_worker_lines.txt") FAILED, $(( $(date +%s) - t0 )) s; red files: $(find "$OUT/data" -name 'red_*.csv' | wc -l | tr -d ' '), event logs: $(find "$OUT/data" -name 'ev_*.csv' | wc -l | tr -d ' ')"
echo "== R2: the same lines again (done before)"
n0=$(find "$OUT/data" -type f | wc -l | tr -d ' ')
tr '\n' '\0' < "$OUT/tasks_short.txt" | xargs -0 -n 1 -P 6 sh -c 'HD_BIN="$0" bash "$1" $2' "$BIN" "$W" > "$OUT/R2_worker_lines.txt" 2>&1
echo "   lines printed: $(wc -l < "$OUT/R2_worker_lines.txt" | tr -d ' '); files before $n0, after $(find "$OUT/data" -type f | wc -l | tr -d ' ')"
echo "== R3: a different binary into a folder written by this build"
l=$(sed -n 3p "$OUT/tasks_short.txt" | awk '{ $5 = $5 + 1; print }')
HD_BIN="$OTHER" bash "$W" $l > "$OUT/R3_worker_lines.txt" 2>&1; echo "   exit $?: $(tail -1 "$OUT/R3_worker_lines.txt")"
echo "== R4: a divider outside the box"
l=$(sed -n 1p "$OUT/tasks_short.txt" | awk '{ $2 = "partA/gen3/R4_outside/x_m2"; $5 = 999; $9 = 2 * $6 + 1; print }')
HD_BIN="$BIN" bash "$W" $l > "$OUT/R4_worker_lines.txt" 2>&1; echo "   exit $?: $(tail -1 "$OUT/R4_worker_lines.txt"); kept: $(ls -d "$OUT"/data/partA/gen3/R4_outside/x_m2/.failed_* 2>/dev/null | wc -l | tr -d ' ') failed folder(s), its log: $(grep -h 'STOP\|ABORT' "$OUT"/data/partA/gen3/R4_outside/x_m2/.failed_*/stdout.log | head -1 | cut -c1-150)"
echo "== R5: two tasks at the registered lengths"
for sel in 'partA/gen3/B_N0100_eta0.050000/x_m2 -2 0 ' 'partA/gen3/B_N0400_eta0.690000/x_p1 1 0 '; do
  l=$(grep -F "A $sel" "$OUT/tasks_full.txt" | head -1); t1=$(date +%s)
  HD_DATA="$OUT/data_full" HD_BIN="$BIN" bash "$W" $l >> "$OUT/R5_worker_lines.txt" 2>&1; echo "   exit $?, $(( $(date +%s) - t1 )) s: $(tail -1 "$OUT/R5_worker_lines.txt")"
done
python3 - "$OUT" > "$OUT/summary.txt" <<'PY'
import glob, os, sys
import pandas as pd
out = sys.argv[1]
fs = sorted(glob.glob(os.path.join(out, "data*", "partA", "gen3", "B_*", "x_*", "red_*.csv")))
if not fs:
    sys.exit("# no red rows")
D = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
cols = ["cell", "j", "seed", "x_div", "T_w", "F_L", "F_R", "n_L", "n_R", "T_L", "T_R", "v_div_t0", "dx_div", "dW_div", "dF_L_rel_log", "dF_R_rel_log", "lattice_left"]
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30); pd.set_option("display.max_colwidth", 60)
print(f"# {len(D)} red rows (R1: T_w = 40 sigma-time; R5: T_w = 5000 sigma-time)\n")
print(D[cols].to_string(index=False))
PY
echo "== the rows: $OUT/summary.txt ($(sed -n 1p "$OUT/summary.txt"))"
