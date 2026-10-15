#!/usr/bin/env bash
# ##CHRIS 2026-10-14 (plan-author decision 50 item 3 d; 261012 sec. 4.7.113; registration sec. 4.8.4): THE PART A WORKER, option (ii): one
# held trajectory (one seed at one stencil position of one cell) per call. The command is the confinement worker's AF command
# (cluster/confinement_20261013/conf_worker.sh, mode AF: --experiment=energy_transfer, the divider held with mass 0 for --wall-hold-steps,
# then the released tail of --steps that writes the temperatures) with THESE CHANGES, all of them:
#   (1) the gen3 engine and option (ii): --engine=gen3 --gen3-exact-box --gen3-seeding=lattice (the registered start and box,
#       sec. 4.8.3.1 and 4.8.4) and --gen3-held-readout=<T_eq steps>: the forces from the engine's per-face impulse counters over
#       (T_eq, T_eq + window] ([EDMD3-HELD] lines). No event log, unless the task asks for one (log = 1: HD_PISTON_EVENTS, the check
#       subset);
#   (2) the hold is T_eq + the window (registered: T_eq of sec. 4.8.3.1, the window 5000 sigma-time); the tail is the AF command's
#       --steps (registered: 20 sigma-time);
#   (3) the task line has 14 fields:
#         A <rel position folder> <j> <s> <seed> <L0> <H> <N_s> <x_div> <T_eq steps> <window steps> <tail steps> <trace every> <log 0|1>
#       L0 and x_div are passed as written (the binary reads both as doubles);
#   (4) the health rule of the stage 1 worker (exactly one build line, exactly one run record [EDMD3-HEALTH] with clean=1, no other HEALTH
#       line), plus exactly one [EDMD3-HELD] window line and one [EDMD3-HELD-GEOMETRY] line;
#   (5) pa_reduce.py reduces the run's record lines (and its event log, if any) into red_<seed>.csv, failing closed.
# The build guard, the locks, the "done before" skip and the failure rule are the stage 1 worker's: a failing trajectory is kept in
# .failed_<seed>_<time>/ and reported (FAILED line), never replaced by another seed. The position folder holds red_<seed>.csv,
# run_<seed>.log, tr_<seed>.csv, summary_<seed>.csv (ev_<seed>.csv with a log); the rest of the run folder goes to .done_runs/.
#   usage: HD_BIN=<abs> HD_DATA=<abs> [HD_BUILD=<line>] bash pa_worker_koa.sh A <the 13 fields after A>
set -uo pipefail
: "${HD_BIN:?set HD_BIN to the frozen 00ALLINONE}"; : "${HD_DATA:?set HD_DATA to the data root}"
case "$HD_BIN" in /*) ;; *) echo "HD_BIN must be an absolute path"; exit 2 ;; esac
case "$HD_DATA" in /*) ;; *) echo "HD_DATA must be an absolute path"; exit 2 ;; esac
HERE="$(cd "$(dirname "$0")" && pwd)"
HEALTH='EDMD-HEALTH|forced_advance|clamp_repair|overlap_repair|wall_overdue'
BUILD="${HD_BUILD:-$("$HD_BIN" --version | head -1)}"
guard() {   # conf_worker.sh's build guard, unchanged (as in s1_worker_koa.sh)
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
health_ok() {   # health_ok <stdout log>: the stage 1 worker's gen3 rule, plus the held path's record lines
  local log=$1 h nb nr nc nw ng
  h=$(grep -E "$HEALTH" "$log" | grep -cv '^\[EDMD3-HEALTH\]') || true; [ "${h:-0}" -eq 0 ] || return 1
  nb=$(grep -c '^\[EDMD3\] built #' "$log") || true; nr=$(grep -c '^\[EDMD3-HEALTH\]' "$log") || true
  nc=$(grep -c '^\[EDMD3-HEALTH\] .*: clean=1 ' "$log") || true
  nw=$(grep -c '^\[EDMD3-HELD\] window ' "$log") || true; ng=$(grep -c '^\[EDMD3-HELD-GEOMETRY\] ' "$log") || true
  [ "$nb" -eq 1 ] && [ "$nr" -eq 1 ] && [ "$nc" -eq 1 ] && [ "$nw" -eq 1 ] && [ "$ng" -eq 1 ] && grep -q '^\[EDMD3-HEALTH\] .* engine=gen3 ' "$log"
}
mode=${1:-}; shift || true
export HD_CONTACT_AUDIT=1
[ "$mode" = A ] || { echo "mode must be A"; exit 2; }
[ $# -eq 13 ] || { echo "A lines have 14 fields; got $(( $# + 1 ))"; exit 2; }
rel=$1 j=$2 s=$3 seed=$4 L0=$5 H=$6 NS=$7 xw=$8 teq=$9 win=${10} tail=${11} every=${12} log=${13}
case "$log" in 0|1) ;; *) echo "A $rel: log must be 0 or 1, not '$log'"; exit 2 ;; esac
pos="$HD_DATA/$rel"; mkdir -p "$pos"
guard "$pos" 'red_*.csv' || { echo "A $rel seed=$seed FAILED build guard"; exit 3; }
[ -s "$pos/red_${seed}.csv" ] && exit 0                                   # done before
tmp="$pos/.run_${seed}"
[ -e "$tmp" ] && mv "$tmp" "$pos/.stale_${seed}_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$tmp"
t0=$SECONDS
cd "$tmp" || exit 1
EV=(); [ "$log" = 1 ] && EV=(env HD_PISTON_EVENTS="$tmp/ev.csv")
${EV[@]+"${EV[@]}"} "$HD_BIN" --mode=edmd --experiment=energy_transfer --headless --quiet \
   --edmd-acc=0 --seed-drift-order=drift-first --energy-transfer-summary="$tmp/summary.csv" \
   --energy-transfer-trace="$tmp/tr.csv" --trace-every=$every --particles=$((2*NS)) \
   --particles-boxes=$NS,$NS --particle-radius=0.5 --l0=$L0 --height=$H --num-walls=1 --wall-positions=$xw \
   --wall-mass-factors=1000000000 --wall-thickness=0.05 --wall-thickness-vis=0.05 --eff-output=wall-ke \
   --wall-hold-steps=$((teq + win)) --steps=$tail --fixed-dt=0.4 --kbt1 --seed=$seed \
   --engine=gen3 --gen3-exact-box --gen3-seeding=lattice --gen3-held-readout=$teq > "$tmp/stdout.log" 2>&1
rc=$?; cd "$pos" || exit 1
fail() { echo "A $rel j=$j s=$s seed=$seed FAILED $1"; mv "$tmp" "$pos/.failed_${seed}_$(date +%Y%m%d_%H%M%S)"; exit 1; }
{ [ "$rc" -eq 0 ] && health_ok "$tmp/stdout.log"; } || fail "rc=$rc engine=gen3"
python3 "$HERE/pa_reduce.py" "$tmp" "$tmp/red.csv" "$rel" "$j" "$s" "$seed" "$NS" "$teq" "$win" || fail "reduction"
mv "$tmp/stdout.log" "run_${seed}.log" && mv "$tmp/tr.csv" "tr_${seed}.csv" && mv "$tmp/summary.csv" "summary_${seed}.csv" || fail "keep"
if [ -e "$tmp/ev.csv" ]; then mv "$tmp/ev.csv" "ev_${seed}.csv" || fail "keep the log"; fi
mv "$tmp/red.csv" "red_${seed}.csv" || fail "keep the row"
mkdir -p .done_runs && mv "$tmp" ".done_runs/run_${seed}"
echo "A $rel j=$j s=$s seed=$seed ok ($((SECONDS - t0)) s)"
