#!/usr/bin/env bash
# ##CHRIS 2026-10-13: one trajectory of the Paper 1 confinement campaign (261012 sec. 1). Called by the
# generated sbatch files through xargs, one line of tasks_<cell>.txt per call. Writes into the SAME relative
# layout as the Mac, under $HD_DATA (KOA scratch), so copied-back data is read by the Mac scripts unchanged.
#
#   B <rel_cell_dir> <M> <r> <exact_seed> <L0> <H> <Ns> <stride> <base>     free divider, speed-of-sound mode
#   A <rel_pos_dir> <x_wall> <seed> <L0> <H> <Ns> <steps> <every>           held divider, energy-transfer mode
#
# Method B mirrors the A1v2 harness (validation/tests_20260913.py run_one): the binary writes run0 into a fresh
# temporary .run<r>/, the trace is renamed to _run<r>.csv in the cell directory and its stdout appended to the
# cell's run.log under a "##RUN" header. Only freshly written output of this call is renamed; nothing else moves.
set -uo pipefail
: "${HD_BIN:?set HD_BIN to the koa-built 00ALLINONE}"; : "${HD_DATA:?set HD_DATA to the data root (scratch)}"
# ##CHRIS 2026-10-02 (Task K2): the binary now runs with its own output directory as working directory, so the files it
# writes into its cwd (run_params.json, 00_COMMAND.md, experiments_speed_of_sound/...) land on scratch next to the data,
# never in the git checkout (a sparse checkout reports -dirty if a tracked path such as energy_log.csv appears in it).
# HD_BIN and HD_DATA must therefore be absolute; reduce_A.py is called by absolute path.
case "$HD_BIN" in /*) ;; *) echo "HD_BIN must be an absolute path"; exit 2 ;; esac
case "$HD_DATA" in /*) ;; *) echo "HD_DATA must be an absolute path"; exit 2 ;; esac
HERE="$(cd "$(dirname "$0")" && pwd)"
HEALTH='EDMD-HEALTH|forced_advance|clamp_repair|overlap_repair|wall_overdue'
mode=$1; shift
if [ "$mode" = B ]; then
  rel=$1 M=$2 r=$3 seed=$4 L0=$5 H=$6 NS=$7 stride=$8 base=$9
  cell="$HD_DATA/$rel"; mkdir -p "$cell"
  ls "$cell"/wall_x_positions_L0_*_wallmassfactor_${M}_run${r}.csv >/dev/null 2>&1 && exit 0     # done before
  tmp="$cell/.run$r"; mkdir -p "$tmp"
  t0=$SECONDS
  cd "$tmp" || exit 1
  HD_KE_TRACE=1 "$HD_BIN" --mode=edmd --experiment=speed_of_sound --headless --kbt1 --seed-drift-order=drift-first \
     --edmd-acc=0 --particles=$((2*NS)) --particles-boxes=$NS,$NS --height=$H --particle-radius=0.5 \
     --wall-thickness=0.05 --wall-thickness-vis=0.05 --lengths=$L0 --wall-masses=$M --repeats=1 --seed=$base \
     --wall-hold-steps=2000 --fixed-dt=0.4 --target-oscillations=200 --oscillation-safety=1.0 \
     --oscillation-min-steps=10000 --oscillation-max-steps=400000000 --speed-sound-log-stride=$stride \
     --speed-sound-run-dir="$tmp" --speed-sound-exact-seed=$seed > "$tmp/stdout.log" 2>&1
  rc=$?; cd "$cell" || exit 1; h=$(grep -cE "$HEALTH" "$tmp/stdout.log" 2>/dev/null) || true
  tr=$(ls "$tmp"/wall_x_positions_L0_*_wallmassfactor_${M}_run0.csv 2>/dev/null | head -1)
  if [ "$rc" -ne 0 ] || [ "${h:-0}" -ne 0 ] || [ -z "$tr" ]; then
    echo "B $rel M=$M r=$r FAILED rc=$rc health=${h:-0}"; mv "$tmp" "$cell/.failed_run${r}_$(date +%Y%m%d_%H%M%S)"; exit 1; fi
  mv "$tr" "$cell/$(basename "${tr%run0.csv}")run$r.csv"
  ( flock 9; { printf '\n##RUN %s run %s seed %s (%s s)\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$r" "$seed" "$((SECONDS-t0))"
               sed -e "s/run = 0,/run = $r,/" "$tmp/stdout.log"; } >> "$cell/run.log" ) 9>"$cell/.runlog.lock"
  rm -rf "$tmp"
elif [ "$mode" = A ]; then
  rel=$1 xw=$2 seed=$3 L0=$4 H=$5 NS=$6 steps=$7 every=$8
  d="$HD_DATA/$rel"; mkdir -p "$d"
  [ -s "$d/red_${seed}.csv" ] && exit 0
  cd "$d" || exit 1
  HD_PISTON_EVENTS="$d/ev_${seed}.csv" "$HD_BIN" --mode=edmd --experiment=energy_transfer --headless --quiet \
     --edmd-acc=0 --seed-drift-order=drift-first --energy-transfer-summary="$d/summary_${seed}.csv" \
     --energy-transfer-trace="$d/tr_${seed}.csv" --trace-every=$every --particles=$((2*NS)) \
     --particles-boxes=$NS,$NS --particle-radius=0.5 --l0=$L0 --height=$H --num-walls=1 --wall-positions=$xw \
     --wall-mass-factors=1000000000 --wall-thickness=0.05 --wall-thickness-vis=0.05 --eff-output=wall-ke \
     --wall-hold-steps=12000 --steps=$steps --fixed-dt=0.4 --kbt1 --seed=$seed > "$d/run_${seed}.log" 2>&1
  rc=$?; h=$(grep -cE "$HEALTH" "$d/run_${seed}.log" 2>/dev/null) || true
  if [ "$rc" -ne 0 ] || [ "${h:-0}" -ne 0 ]; then echo "A $rel seed=$seed FAILED rc=$rc health=${h:-0}"; exit 1; fi
  python3 "$HERE/reduce_A.py" "$d/ev_${seed}.csv" "$d/tr_${seed}.csv" "$d/red_${seed}.csv" || { echo "A $rel seed=$seed reduction FAILED"; exit 1; }
else
  echo "mode must be A or B"; exit 2
fi
