#!/usr/bin/env bash
# ##CHRIS 2026-10-13: one trajectory of the Paper 1 confinement campaign (261012 sec. 1). Called by the
# generated sbatch files through xargs, one line of tasks_<cell>.txt per call. Writes into the SAME relative
# layout as the Mac, under $HD_DATA (KOA scratch), so copied-back data is read by the Mac scripts unchanged.
#
#   B <rel_cell_dir> <M> <r> <exact_seed> <L0> <H> <Ns> <stride> <base>     free divider, speed-of-sound mode
#   A <rel_pos_dir> <x_wall> <seed> <L0> <H> <Ns> <steps> <every>           held divider, energy-transfer mode
#   AF <rel_pos_dir> <x_wall> <seed> <L0> <H> <Ns> <hold> <post> <every>  divider held for the WHOLE record (sec. 3)
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
# ##CHRIS 2026-10-02 (Task U4): build guard. Every output directory records the build that wrote it in .build_git (the
# first line of `00ALLINONE --version`, e.g. "00ALLINONE  git 70b2069  target koa"); the first worker writes it, under
# flock. A worker whose binary differs from the record REFUSES the directory, and so does a worker that finds outputs
# but no record (written by an unknown build, e.g. before this guard). The sbatch exports HD_BUILD once.
BUILD="${HD_BUILD:-$("$HD_BIN" --version | head -1)}"
# ##CHRIS 2026-10-03 (Task W1): the U4 version locked with `have=$( flock 9 ... ) 9>"$dir/.build_git.lock"`. On an
# assignment the command substitution is expanded BEFORE the redirection opens fd 9, so flock got no file ("flock: 9: Bad
# file descriptor", once per trajectory in Round 1) and never locked. Parallel workers then raced on `printf > .build_git`
# (truncate, then write) against `cat`, and a worker that read the file in between saw another "build" and refused: the
# "FAILED build guard" lines of Round 1. Those workers exited before writing anything; the seeds were skipped, not spoiled.
# Now: a mkdir lock (atomic on every file system, testable on the Mac), up to 60 s of waiting, released by rmdir, also
# on SIGTERM (Slurm's TIMEOUT); .build_git is written to a temporary name and renamed (mv), so it is never seen half
# written. A lock older than the wait is NOT broken: the worker refuses (the seed stays missing and
# cluster/check_cells.sh lists the lock), because breaking a lock cannot be made race-free here.
guard() {   # guard <dir> <glob of finished outputs, relative to dir>
  local dir=$1 lock="$1/.guard.lock" have="" i=0
  until mkdir "$lock" 2>/dev/null; do
    i=$((i + 1)); [ "$i" -le 600 ] || { echo "REFUSED $dir: lock $lock not free within 60 s"; return 1; }
    sleep 0.1
  done
  trap 'rmdir "$lock" 2>/dev/null' EXIT
  trap 'rmdir "$lock" 2>/dev/null; exit 143' TERM INT
  if [ ! -e "$dir/.build_git" ]; then
    if compgen -G "$dir/$2" >/dev/null; then have="(none recorded, outputs present)"
    else printf '%s\n' "$BUILD" > "$dir/.build_git.tmp$$" && mv "$dir/.build_git.tmp$$" "$dir/.build_git"; fi
  fi
  [ -n "$have" ] || have=$(cat "$dir/.build_git")
  rmdir "$lock"; trap - EXIT TERM INT
  [ "$have" = "$BUILD" ] && return 0
  echo "REFUSED $dir: written by '$have', this binary is '$BUILD'"; return 1
}
# ##CHRIS 2026-10-05 (engine gate G-E6, 261012 sec. 4.4): build-GENERATION guard on the whole data root. guard() above keeps
# two builds out of one directory; this keeps them out of one ROOT, so that cells of the 279282b campaign and cells of the
# engine-divider-resched generation can never end up side by side in one tree (and later in one figure). The root records
# its build in $HD_DATA/.build_generation, written by the first worker into an EMPTY root (no .build_git anywhere below).
# A root that holds data but no record is refused: its build is unknown (the 279282b root on KOA is such a root -- write
# its record by hand before using it again, runsheet). A different build is refused unless HD_ALLOW_BUILD_MIX=1 (explicit).
# The sbatch files of the branch derive HD_DATA from the clone's name ($SCRATCH/<clone>/hspist3), so ~/harddisks_resched
# writes to its own root by construction.
root_guard() {
  local f="$HD_DATA/.build_generation" lock="$HD_DATA/.build_generation.lock" have="" i=0
  if [ ! -e "$f" ]; then
    mkdir -p "$HD_DATA"
    until mkdir "$lock" 2>/dev/null; do
      i=$((i + 1)); [ "$i" -le 600 ] || { echo "REFUSED root $HD_DATA: lock $lock not free within 60 s"; return 1; }
      sleep 0.1
    done
    trap 'rmdir "$lock" 2>/dev/null' EXIT
    trap 'rmdir "$lock" 2>/dev/null; exit 143' TERM INT
    if [ ! -e "$f" ]; then
      if [ -n "$(find "$HD_DATA" -maxdepth 9 -name .build_git -print -quit 2>/dev/null)" ]; then
        rmdir "$lock"; trap - EXIT TERM INT
        echo "REFUSED root $HD_DATA: it holds data but no .build_generation record (build unknown); write the record by hand"
        return 1
      fi
      printf '%s\n' "$BUILD" > "$f.tmp$$" && mv "$f.tmp$$" "$f"
    fi
    rmdir "$lock"; trap - EXIT TERM INT
  fi
  have=$(cat "$f")
  [ "$have" = "$BUILD" ] && return 0
  if [ "${HD_ALLOW_BUILD_MIX:-0}" = 1 ]; then
    echo "WARNING root $HD_DATA: generation '$have', this binary '$BUILD' -- mixing ALLOWED by HD_ALLOW_BUILD_MIX=1"; return 0
  fi
  echo "REFUSED root $HD_DATA: build generation '$have', this binary is '$BUILD' (HD_ALLOW_BUILD_MIX=1 overrides)"; return 1
}
root_guard || { echo "$* FAILED root build guard"; exit 3; }
mode=$1; shift
if [ "$mode" = B ]; then
  rel=$1 M=$2 r=$3 seed=$4 L0=$5 H=$6 NS=$7 stride=$8 base=$9
  cell="$HD_DATA/$rel"; mkdir -p "$cell"
  guard "$cell" 'wall_x_positions_L0_*_run*.csv' || { echo "B $rel M=$M r=$r FAILED build guard"; exit 3; }
  ls "$cell"/wall_x_positions_L0_*_wallmassfactor_${M}_run${r}.csv >/dev/null 2>&1 && exit 0     # done before
  tmp="$cell/.run$r"
  # ##CHRIS 2026-10-03 (Task W1): a .run<r> left by a killed run (Round 1 TIMEOUT) is moved aside, neither reused nor
  # deleted (the binary appends to speed_of_sound_psi6.csv in it, 00ALLINONE.c:15355)
  [ -e "$tmp" ] && mv "$tmp" "$cell/.stale_run${r}_$(date +%Y%m%d_%H%M%S)"
  mkdir -p "$tmp"
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
  guard "$d" 'red_*.csv' || { echo "A $rel seed=$seed FAILED build guard"; exit 3; }
  [ -s "$d/red_${seed}.csv" ] && exit 0
  # ##CHRIS 2026-10-03 (Task W1): files of this seed left by a killed run (no non-empty red_: Round 1 TIMEOUT) are moved
  # aside before the rerun, neither overwritten nor deleted -- the binary APPENDS to summary_<seed>.csv
  # (00ALLINONE.c:17287, fopen "a"), which would otherwise carry two rows; ev_ and tr_ are opened "w".
  if compgen -G "$d/*_${seed}.*" >/dev/null; then
    st="$d/.stale_${seed}_$(date +%Y%m%d_%H%M%S)"; mkdir -p "$st"; mv "$d"/*_"${seed}".* "$st"/
  fi
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
elif [ "$mode" = AF ]; then
  # ##CHRIS 2026-10-04 (Task Y, 261012 sec. 3 "A-fixed"): method A with the divider HELD for the whole record. During
  # --wall-hold-steps the divider has mass 0 = immovable (00ALLINONE.c:16890-16892, edmd.c:779, 1174-1183) and every
  # divider collision is still logged (edmd.c:1183). hold = 12000 (equilibration, 200 sigma-time) + 300000 (record,
  # 5000 sigma-time); post = a short released tail only so the trace (written after release only, 00ALLINONE.c:17069)
  # records the temperatures. reduce_AF.py uses the HELD window [200, hold*dt) only.
  rel=$1 xw=$2 seed=$3 L0=$4 H=$5 NS=$6 hold=$7 post=$8 every=$9
  d="$HD_DATA/$rel"; mkdir -p "$d"
  guard "$d" 'red_*.csv' || { echo "AF $rel seed=$seed FAILED build guard"; exit 3; }
  [ -s "$d/red_${seed}.csv" ] && exit 0
  if compgen -G "$d/*_${seed}.*" >/dev/null; then
    st="$d/.stale_${seed}_$(date +%Y%m%d_%H%M%S)"; mkdir -p "$st"; mv "$d"/*_"${seed}".* "$st"/
  fi
  cd "$d" || exit 1
  HD_PISTON_EVENTS="$d/ev_${seed}.csv" "$HD_BIN" --mode=edmd --experiment=energy_transfer --headless --quiet \
     --edmd-acc=0 --seed-drift-order=drift-first --energy-transfer-summary="$d/summary_${seed}.csv" \
     --energy-transfer-trace="$d/tr_${seed}.csv" --trace-every=$every --particles=$((2*NS)) \
     --particles-boxes=$NS,$NS --particle-radius=0.5 --l0=$L0 --height=$H --num-walls=1 --wall-positions=$xw \
     --wall-mass-factors=1000000000 --wall-thickness=0.05 --wall-thickness-vis=0.05 --eff-output=wall-ke \
     --wall-hold-steps=$hold --steps=$post --fixed-dt=0.4 --kbt1 --seed=$seed > "$d/run_${seed}.log" 2>&1
  rc=$?; h=$(grep -cE "$HEALTH" "$d/run_${seed}.log" 2>/dev/null) || true
  if [ "$rc" -ne 0 ] || [ "${h:-0}" -ne 0 ]; then echo "AF $rel seed=$seed FAILED rc=$rc health=${h:-0}"; exit 1; fi
  t1=$(awk -v h="$hold" 'BEGIN{printf "%.9f", h * 0.4 / 24.0}')
  python3 "$HERE/reduce_AF.py" "$d/ev_${seed}.csv" "$d/tr_${seed}.csv" "$d/red_${seed}.csv" 200 "$t1" || { echo "AF $rel seed=$seed reduction FAILED"; exit 1; }
else
  echo "mode must be A, AF or B"; exit 2
fi
