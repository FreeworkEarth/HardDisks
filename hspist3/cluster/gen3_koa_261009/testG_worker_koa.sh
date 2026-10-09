#!/usr/bin/env bash
# ##CHRIS 2026-10-09 (261012 sec. 4.7.20; stage G): PREPARED, NOT RUN. The KOA copy of the registered Mac worker
# cluster/gen3_gate_261009/testG_worker.sh (that file is not edited after its registration). Differences, all for Linux: the
# SHA-256 tool (sha256sum; shasum -a 256 is a Mac tool); the A-fixed event log is always kept (KEEP_EV=1 from testG_koa.sbatch:
# scratch has the room), so the FIFO path is not used. Everything else -- the commands, the health rule, the build guard, the
# ##RUN header, the folders -- is the Mac worker's.
# The Mac worker's own description follows.
# ##CHRIS 2026-10-09 (261012 sec. 4.7.15, Test G): one trajectory of Test G on the Mac. The campaign worker's B and AF commands
# (cluster/confinement_20261013/conf_worker.sh, the commands of Test T-prime and of the A-fixed campaign) unchanged except
# --engine=<gen2|gen3> (and the seeds, from the task list). Called by run_testG_mac.py, one task line per call:
#   B  <rel_cell_dir> <M> <r> <exact_seed> <L0> <H> <Ns> <stride> <base> <engine> <block>
#   AF <rel_pos_dir> <x_wall> <seed> <L0> <H> <Ns> <hold> <post> <every> <engine> afix
# As conf_worker.sh: the binary runs in a fresh .run<r>/ (B) with its own directory as working directory; the trace is renamed to
# _run<r>.csv in the cell directory and the stdout appended to the cell's run.log under a "##RUN ... node <host> (<s> s)" header
# (T-prime's header with the node); the build guard (.build_git, the first line of --version) refuses a directory written by
# another build. HD_CONTACT_AUDIT=1 and HD_KE_TRACE=1 in every trajectory (T-prime).
# HEALTH (programme rule 10): gen2 as conf_worker.sh -- no line matching its HEALTH pattern; gen3 -- exactly one build line, exactly
# one run record [EDMD3-HEALTH] with clean=1, and no other line matching the HEALTH pattern. A failing trajectory is kept in
# .failed_run<r>_<time>/ (B) or .failed_<seed>_<time>/ (AF) and reported (FAILED line); the verdict's inventory counts it; it is
# never rerun by this worker (a rerun would need the task again, which run_testG_mac.py does not do).
# DISK (programme rule 6): gen3's psi6(t) file is compressed losslessly right after its trajectory (gzip -9; it is not read by any
# reduction); the traces after the cell's reduction (run_testG_mac.py). The SHA-256 of every uncompressed file is recorded in the
# cell's .sha256_uncompressed before compression. AF: the event log goes through a FIFO into af_stream_reduce.py (the registered
# reduce_AF.py on the stream), except for r < 10 per engine (KEEP_EV, the first 10 seeds of the task list), whose event log is
# written to disk, reduced by reduce_AF.py from the file, then compressed.
set -uo pipefail
: "${HD_BIN:?set HD_BIN to the frozen 00ALLINONE}"; : "${HD_DATA:?set HD_DATA to the data root}"
case "$HD_BIN" in /*) ;; *) echo "HD_BIN must be an absolute path"; exit 2 ;; esac
case "$HD_DATA" in /*) ;; *) echo "HD_DATA must be an absolute path"; exit 2 ;; esac
HERE="$(cd "$(dirname "$0")" && pwd)"; CF="$(dirname "$HERE")/confinement_20261013"
HEALTH='EDMD-HEALTH|forced_advance|clamp_repair|overlap_repair|wall_overdue'
BUILD="${HD_BUILD:-$("$HD_BIN" --version | head -1)}"
HOST=$(hostname -s)
guard() {   # conf_worker.sh's build guard, unchanged
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
health_ok() {   # health_ok <engine> <stdout log>: 0 if clean by the rule above
  local eng=$1 log=$2 h nb nr nc
  if [ "$eng" = gen2 ]; then
    h=$(grep -cE "$HEALTH" "$log") || true; [ "${h:-0}" -eq 0 ] || return 1
    ! grep -q '^\[EDMD3' "$log"
  else
    h=$(grep -E "$HEALTH" "$log" | grep -cv '^\[EDMD3-HEALTH\]') || true; [ "${h:-0}" -eq 0 ] || return 1
    nb=$(grep -c '^\[EDMD3\] built #' "$log") || true; nr=$(grep -c '^\[EDMD3-HEALTH\]' "$log") || true
    nc=$(grep -c '^\[EDMD3-HEALTH\] .*: clean=1 ' "$log") || true
    [ "$nb" -eq 1 ] && [ "$nr" -eq 1 ] && [ "$nc" -eq 1 ] && grep -q '^\[EDMD3-HEALTH\] .* engine=gen3 ' "$log"
  fi
}
record_sha() { ( cd "$(dirname "$1")" && sha256sum "$(basename "$1")" ) >> "$2"; }   # KOA: sha256sum
with_lock() {   # with_lock <lock dir> <command...>: macOS has no flock; an atomic mkdir lock, as guard(), up to 600 s
  local lock=$1 i=0 rc; shift
  until mkdir "$lock" 2>/dev/null; do i=$((i + 1)); [ "$i" -le 6000 ] || { echo "lock $lock not free within 600 s"; return 1; }; sleep 0.1; done
  "$@"; rc=$?
  rmdir "$lock"; return $rc
}
append_run() {  # append_run <cell> <r> <seed> <t0> <stdout log>: the ##RUN section (conf_worker.sh's, plus the node)
  { printf '\n##RUN %s run %s seed %s node %s (%s s)\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$2" "$3" "$HOST" "$((SECONDS-$4))"
    sed -e "s/run = 0,/run = $2,/" -e "/^\[EDMD3-HEALTH\] L0=/s/ run=0 seed=/ run=$2 seed=/" "$5"; } >> "$1/run.log"
}   # the run index of the cell, also in gen3's run record (its id), as conf_worker.sh does for "run = 0,"; raw stdout in .done_runs/
mode=$1; shift
export HD_CONTACT_AUDIT=1
if [ "$mode" = B ]; then
  rel=$1 M=$2 r=$3 seed=$4 L0=$5 H=$6 NS=$7 stride=$8 base=$9 eng=${10}
  cell="$HD_DATA/$rel"; mkdir -p "$cell"
  guard "$cell" 'wall_x_positions_L0_*_run*.csv*' || { echo "B $rel M=$M r=$r FAILED build guard"; exit 3; }
  ls "$cell"/wall_x_positions_L0_*_wallmassfactor_${M}_run${r}.csv* >/dev/null 2>&1 && exit 0     # done before
  tmp="$cell/.run$r"
  [ -e "$tmp" ] && mv "$tmp" "$cell/.stale_run${r}_$(date +%Y%m%d_%H%M%S)"
  mkdir -p "$tmp"
  t0=$SECONDS
  cd "$tmp" || exit 1
  HD_KE_TRACE=1 "$HD_BIN" --mode=edmd --experiment=speed_of_sound --headless --kbt1 --seed-drift-order=drift-first \
     --edmd-acc=0 --particles=$((2*NS)) --particles-boxes=$NS,$NS --height=$H --particle-radius=0.5 \
     --wall-thickness=0.05 --wall-thickness-vis=0.05 --lengths=$L0 --wall-masses=$M --repeats=1 --seed=$base \
     --wall-hold-steps=2000 --fixed-dt=0.4 --target-oscillations=200 --oscillation-safety=1.0 \
     --oscillation-min-steps=10000 --oscillation-max-steps=400000000 --speed-sound-log-stride=$stride \
     --speed-sound-run-dir="$tmp" --speed-sound-exact-seed=$seed --engine=$eng > "$tmp/stdout.log" 2>&1
  rc=$?; cd "$cell" || exit 1
  tr=$(ls "$tmp"/wall_x_positions_L0_*_wallmassfactor_${M}_run0.csv 2>/dev/null | head -1)
  if [ "$rc" -ne 0 ] || ! health_ok "$eng" "$tmp/stdout.log" || [ -z "$tr" ]; then
    echo "B $rel M=$M r=$r FAILED rc=$rc engine=$eng"; mv "$tmp" "$cell/.failed_run${r}_$(date +%Y%m%d_%H%M%S)"; exit 1; fi
  mv "$tr" "$cell/$(basename "${tr%run0.csv}")run$r.csv"
  ps=$(ls "$tmp"/psi6_t_L0_*_wallmassfactor_${M}_run0.csv 2>/dev/null | head -1)
  if [ -n "$ps" ]; then
    dst="$cell/$(basename "${ps%run0.csv}")run$r.csv"; mv "$ps" "$dst"
    with_lock "$cell/.sha.lockdir" record_sha "$dst" "$cell/.sha256_uncompressed" || { echo "B $rel M=$M r=$r FAILED sha lock"; exit 1; }
    gzip -9 "$dst"
  fi
  with_lock "$cell/.runlog.lockdir" append_run "$cell" "$r" "$seed" "$t0" "$tmp/stdout.log" || { echo "B $rel M=$M r=$r FAILED run.log lock"; exit 1; }
  mkdir -p "$cell/.done_runs"; mv "$tmp" "$cell/.done_runs/run$r"      # the rest of the run folder (no trace left in it), kept
elif [ "$mode" = AF ]; then
  rel=$1 xw=$2 seed=$3 L0=$4 H=$5 NS=$6 hold=$7 post=$8 every=$9 eng=${10}
  d="$HD_DATA/$rel"; mkdir -p "$d"
  guard "$d" 'red_*.csv' || { echo "AF $rel seed=$seed FAILED build guard"; exit 3; }
  [ -s "$d/red_${seed}.csv" ] && exit 0
  if compgen -G "$d/*_${seed}.*" >/dev/null; then
    st="$d/.stale_${seed}_$(date +%Y%m%d_%H%M%S)"; mkdir -p "$st"; mv "$d"/*_"${seed}".* "$st"/
  fi
  t1=$(awk -v h="$hold" 'BEGIN{printf "%.9f", h * 0.4 / 24.0}')
  keep=0; [ -n "${KEEP_EV:-}" ] && keep=1
  cd "$d" || exit 1
  if [ "$keep" = 1 ]; then ev="$d/ev_${seed}.csv"
  else ev="$d/.ev_${seed}.fifo"; mk="$d/.ev_${seed}.exited"; mkfifo "$ev" || exit 1
       python3 "$HERE/af_stream_reduce.py" "$ev" "$mk" "$d/tr_${seed}.csv" "$d/red_${seed}.csv" 200 "$t1" > "$d/stream_${seed}.log" 2>&1 &
       SPID=$!
  fi
  t0=$SECONDS
  HD_PISTON_EVENTS="$ev" "$HD_BIN" --mode=edmd --experiment=energy_transfer --headless --quiet \
     --edmd-acc=0 --seed-drift-order=drift-first --energy-transfer-summary="$d/summary_${seed}.csv" \
     --energy-transfer-trace="$d/tr_${seed}.csv" --trace-every=$every --particles=$((2*NS)) \
     --particles-boxes=$NS,$NS --particle-radius=0.5 --l0=$L0 --height=$H --num-walls=1 --wall-positions=$xw \
     --wall-mass-factors=1000000000 --wall-thickness=0.05 --wall-thickness-vis=0.05 --eff-output=wall-ke \
     --wall-hold-steps=$hold --steps=$post --fixed-dt=0.4 --kbt1 --seed=$seed --engine=$eng > "$d/run_${seed}.log" 2>&1
  rc=$?
  if [ "$keep" = 1 ]; then
    rrc=1; [ "$rc" -eq 0 ] && { python3 "$CF/reduce_AF.py" "$ev" "$d/tr_${seed}.csv" "$d/red_${seed}.csv" 200 "$t1"; rrc=$?; }
    with_lock "$d/.sha.lockdir" record_sha "$ev" "$d/.sha256_uncompressed"; gzip -9 "$ev"
  else
    touch "$mk"
    python3 -c "import os,sys
try: os.close(os.open(sys.argv[1], os.O_WRONLY | os.O_NONBLOCK))   # a reader still blocked in open (no writer came): give it EOF
except OSError: pass" "$ev"
    wait "$SPID"; rrc=$?
    mkdir -p "$d/.fifo_used"; mv "$ev" "$mk" "$d/.fifo_used/"
  fi
  printf '##RUN %s seed %s node %s engine %s (%s s) reduce exit %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$seed" "$HOST" "$eng" "$((SECONDS-t0))" "$rrc" >> "$d/run_${seed}.log"
  if [ "$rc" -ne 0 ] || ! health_ok "$eng" "$d/run_${seed}.log" || [ "$rrc" -ne 0 ] || [ ! -s "$d/red_${seed}.csv" ]; then
    st="$d/.failed_${seed}_$(date +%Y%m%d_%H%M%S)"; mkdir -p "$st"; mv "$d"/*_"${seed}".* "$st"/
    echo "AF $rel seed=$seed FAILED rc=$rc reduce=$rrc engine=$eng"; exit 1
  fi
else
  echo "mode must be B or AF"; exit 2
fi
