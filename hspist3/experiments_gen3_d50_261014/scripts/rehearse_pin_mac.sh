#!/bin/bash
# ##CHRIS 2026-10-14 (plan-author decision 53 item 2; 261012 sec. 4.7.118): THE WORKER WITH THE ORIENTATION PIN on the Mac (bin_8f538a2;
# no C source changed for the pin: --gen3-lattice-orientation exists since M4c). Into a NEW data root <out>/data:
#   P1  the 15-field task lines of pa_tasks.py for N 100 and N 400 at eta pi/8, 0.05, 0.66, 0.68 (the stencil's five positions, one seed,
#       short lengths: T_eq 600, window 2400, tail 60 steps; the log on s = 0, j = 0), 6 workers at a time;
#   P2  one line with the orientation "default" (no pin; its own folder and seed);
# then every row's pin and lattice per compartment (summary.txt): each cell must show its pinned crystal orientation at all positions.
#   usage (from hspist3/): bash experiments_gen3_d50_261014/scripts/rehearse_pin_mac.sh <new folder>
set -u
OUT=${1:?a new folder}; [ -e "$OUT" ] && { echo "REFUSED: $OUT exists"; exit 1; }
HS=$(pwd); case "$OUT" in /*) ;; *) OUT="$HS/$OUT" ;; esac
W="$HS/cluster/gen3_partA_261014/pa_worker_koa.sh"; BIN="$HS/experiments_gen3_d50_261014/bin_8f538a2/00ALLINONE"
mkdir -p "$OUT/data"; export HD_DATA="$OUT/data" SDL_VIDEODRIVER=dummy
python3 cluster/gen3_partA_261014/pa_tasks.py --sizes 100,400 --seeds 1 --etas 0.392699,0.05,0.66,0.68 --log s0j0 \
        --steps 600,2400,60,6 > "$OUT/tasks_pin.txt" 2> "$OUT/tasks_pin.stderr"
awk 'NR == 1 { sub("/x_m2$", "/x_m2_nopin", $2); $5 = $5 + 7; $15 = "default"; print }' "$OUT/tasks_pin.txt" > "$OUT/tasks_default.txt"
echo "== P1: $(wc -l < "$OUT/tasks_pin.txt" | tr -d ' ') pinned tasks ($(cut -d' ' -f15 "$OUT/tasks_pin.txt" | sort | uniq -c | tr -s ' ' | tr '\n' ';'))"
tr '\n' '\0' < "$OUT/tasks_pin.txt" | xargs -0 -n 1 -P 6 sh -c 'HD_BIN="$0" bash "$1" $2' "$BIN" "$W" > "$OUT/P1_worker_lines.txt" 2>&1
echo "   $(grep -c ' ok (' "$OUT/P1_worker_lines.txt") ok, $(grep -c FAILED "$OUT/P1_worker_lines.txt") FAILED"
echo "== P2: one line without a pin"
HD_BIN="$BIN" bash "$W" $(cat "$OUT/tasks_default.txt") > "$OUT/P2_worker_lines.txt" 2>&1; echo "   exit $?: $(tail -1 "$OUT/P2_worker_lines.txt")"
python3 - "$OUT" > "$OUT/summary.txt" <<'PY'
import glob, os, re, sys
import pandas as pd
out = sys.argv[1]
D = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(out, "data", "partA", "gen3", "B_*", "x_*", "red_*.csv")))], ignore_index=True)
def cr(s):   # the crystal orientation of a run header's compartment, as g3_lat_cryst
    m = re.search(r"orient (\w+) rows (\w+) n_s (\d+) n_r (\d+) a_s (\S+) a_r (\S+) ", s)
    o, ns, nr, a_s, a_r = m.group(1), int(m.group(3)), int(m.group(4)), float(m.group(5)), float(m.group(6))
    along = 1 if ns < 2 else (0 if nr < 2 else int(a_r < 2 * a_s))
    return f"{m.group(2)[:3]} {ns}x{nr} {'0' if along == int(o == 'perpendicular') else 'turned'}"
D["left"] = D["lattice_left"].map(cr); D["right"] = D["lattice_right"].map(cr)
print(D[["cell", "j", "seed", "pin", "pin_left", "pin_right", "left", "right", "F_L", "F_R", "dF_L_rel_log", "dF_R_rel_log"]].sort_values(["cell", "j"]).to_string(index=False))
bad = D[(D["pin"] != "default") & ((D["left"].str.split().str[-1] != D["pin"].astype(str)) | (D["right"].str.split().str[-1] != D["pin"].astype(str)))]
print(f"\n{len(D)} rows; rows whose compartments do not carry the pin: {len(bad)}")
PY
tail -1 "$OUT/summary.txt"
