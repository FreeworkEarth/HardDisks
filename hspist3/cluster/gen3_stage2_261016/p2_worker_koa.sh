#!/usr/bin/env bash
# ##CHRIS 2026-10-15 (plan-author decision 57; notes 261012 sec. 4.7.135): THE WORKER OF THE STAGE 2 PILOT P2, one trajectory per call.
# It is cluster/gen3_stage1_261012/s1_worker_koa.sh (the stage 1 worker, as stage 1 ran it) with ONE change: the command also writes the
# initial state after the seeding, --gen3-dump-initial=<the run folder>/initial_state.txt (gen3-initial v1: px, C hex floats), which
# stays in the run folder and so moves with it to .done_runs/run<s> (the stage 1 merge moves that folder whole). Every other line --
# the 15-field line, the command, the pin (the start field 0 | turned: --gen3-lattice-orientation), the snapshots (T_snap), the build
# guard, the health rule, the ##RUN section, the locks, the "done before" skip, the kept and compressed files -- is the stage 1 worker's.
# The stage 1 worker's own description follows.
# ##CHRIS 2026-10-12 (plan-author decision 39 item 5 f; registration 261012 sec. 4.8.3; gate G1): THE STAGE 1 WORKER, Part B: one
# trajectory per call. It is cluster/gen3_koa_261009/testG_worker_koa.sh (unchanged since 1cdb7ff6; the worker of Test G, Test R and the
# preview on KOA) with THESE CHANGES, all of them:
#   (1) the B task line has the 15 fields of gen_stage1_koa_261012.py:
#         B <rel_cell_dir> <M> <s> <exact_seed> <L0> <H> <Ns> <stride> <base> <hold steps> <Delta_psi> <T_snap | 0> <default|0|turned> <block>
#       and the engine is gen3 always (no engine field);
#   (2) the command: Test G's B command with the changes of sec. 4.8.3.2 -- --wall-hold-steps=<60 T_eq> (Test G: 2000), --engine=gen3
#       --gen3-seeding=lattice --gen3-exact-box --gen3-psi6-compartments --psi6-every=<Delta_psi> --gen3-virial-blocks=50, and
#       --gen3-snapshots=<T_snap> where T_snap is not 0 (the snapshot subset), --gen3-lattice-orientation=<0|turned> in block S;
#   (3) the run's virial-block file (virial_blocks_*_run0.csv) and snapshot file (snapshots_*_run0.bin) go into the cell folder as the
#       psi6(t) file does: renamed to _run<s>, the SHA-256 of the uncompressed file appended to .sha256_uncompressed, gzip -9;
#   (4) the A-fixed mode (AF) is not here: Part A waits for the plan author's amendment (decision 39 item 5 e).
# Every other line -- the build guard, the health rule (gen3: exactly one build line, exactly one run record [EDMD3-HEALTH] with clean=1,
# no other line matching the HEALTH pattern), the ##RUN section, the locks, the run folder, the "done before" skip -- is the Test G
# worker's. A failing trajectory is kept in .failed_run<s>_<time>/ and reported (FAILED line); it is never replaced by another seed.
#   usage: HD_BIN=<abs> HD_DATA=<abs> [HD_BUILD=<line>] bash s1_worker_koa.sh B <the 14 fields after B>
set -uo pipefail
: "${HD_BIN:?set HD_BIN to the frozen 00ALLINONE}"; : "${HD_DATA:?set HD_DATA to the data root}"
case "$HD_BIN" in /*) ;; *) echo "HD_BIN must be an absolute path"; exit 2 ;; esac
case "$HD_DATA" in /*) ;; *) echo "HD_DATA must be an absolute path"; exit 2 ;; esac
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
health_ok() {   # health_ok <stdout log>: 0 if clean by the gen3 rule above (the Test G worker's gen3 branch)
  local log=$1 h nb nr nc
  h=$(grep -E "$HEALTH" "$log" | grep -cv '^\[EDMD3-HEALTH\]') || true; [ "${h:-0}" -eq 0 ] || return 1
  nb=$(grep -c '^\[EDMD3\] built #' "$log") || true; nr=$(grep -c '^\[EDMD3-HEALTH\]' "$log") || true
  nc=$(grep -c '^\[EDMD3-HEALTH\] .*: clean=1 ' "$log") || true
  [ "$nb" -eq 1 ] && [ "$nr" -eq 1 ] && [ "$nc" -eq 1 ] && grep -q '^\[EDMD3-HEALTH\] .* engine=gen3 ' "$log"
}
record_sha() { ( cd "$(dirname "$1")" && sha256sum "$(basename "$1")" ) >> "$2"; }   # KOA: sha256sum
with_lock() {   # with_lock <lock dir> <command...>: an atomic mkdir lock, as guard(), up to 600 s
  local lock=$1 i=0 rc; shift
  until mkdir "$lock" 2>/dev/null; do i=$((i + 1)); [ "$i" -le 6000 ] || { echo "lock $lock not free within 600 s"; return 1; }; sleep 0.1; done
  "$@"; rc=$?
  rmdir "$lock"; return $rc
}
append_run() {  # append_run <cell> <r> <seed> <t0> <stdout log>: the ##RUN section (conf_worker.sh's, plus the node)
  { printf '\n##RUN %s run %s seed %s node %s (%s s)\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$2" "$3" "$HOST" "$((SECONDS-$4))"
    sed -e "s/run = 0,/run = $2,/" -e "/^\[EDMD3-HEALTH\] L0=/s/ run=0 seed=/ run=$2 seed=/" "$5"; } >> "$1/run.log"
}
keep() {   # keep <file in the run folder> <cell> <r>: rename to _run<r> in the cell, SHA-256 of the uncompressed file, gzip -9 (change 3)
  local f=$1 cell=$2 r=$3 dst
  dst="$cell/$(basename "${f%run0.${f##*.}}")run$r.${f##*.}"; mv "$f" "$dst"
  with_lock "$cell/.sha.lockdir" record_sha "$dst" "$cell/.sha256_uncompressed" || return 1
  gzip -9 "$dst"
}
mode=$1; shift
export HD_CONTACT_AUDIT=1
[ "$mode" = B ] || { echo "mode must be B (Part A waits; decision 39 item 5 e)"; exit 2; }
[ $# -eq 14 ] || { echo "B lines have 15 fields; got $(( $# + 1 ))"; exit 2; }
rel=$1 M=$2 r=$3 seed=$4 L0=$5 H=$6 NS=$7 stride=$8 base=$9 hold=${10} dpsi=${11} tsnap=${12} start=${13} blk=${14}
EXTRA=()
[ "$tsnap" = 0 ] || EXTRA+=("--gen3-snapshots=$tsnap")
case "$start" in default) ;; 0|turned) EXTRA+=("--gen3-lattice-orientation=$start") ;; *) echo "B $rel: start '$start' is not default, 0 or turned"; exit 2 ;; esac
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
   --wall-hold-steps=$hold --fixed-dt=0.4 --target-oscillations=200 --oscillation-safety=1.0 \
   --oscillation-min-steps=10000 --oscillation-max-steps=400000000 --speed-sound-log-stride=$stride \
   --speed-sound-run-dir="$tmp" --speed-sound-exact-seed=$seed --engine=gen3 \
   --gen3-seeding=lattice --gen3-exact-box --gen3-psi6-compartments --psi6-every=$dpsi --gen3-virial-blocks=50 \
   --gen3-dump-initial="$tmp/initial_state.txt" ${EXTRA[@]+"${EXTRA[@]}"} > "$tmp/stdout.log" 2>&1
rc=$?; cd "$cell" || exit 1
tr=$(ls "$tmp"/wall_x_positions_L0_*_wallmassfactor_${M}_run0.csv 2>/dev/null | head -1)
if [ "$rc" -ne 0 ] || ! health_ok "$tmp/stdout.log" || [ -z "$tr" ]; then
  echo "B $rel M=$M r=$r FAILED rc=$rc engine=gen3"; mv "$tmp" "$cell/.failed_run${r}_$(date +%Y%m%d_%H%M%S)"; exit 1; fi
mv "$tr" "$cell/$(basename "${tr%run0.csv}")run$r.csv"
for f in "$tmp"/psi6_t_L0_*_wallmassfactor_${M}_run0.csv "$tmp"/virial_blocks_L0_*_wallmassfactor_${M}_run0.csv \
         "$tmp"/snapshots_L0_*_wallmassfactor_${M}_run0.bin; do
  [ -e "$f" ] || continue
  keep "$f" "$cell" "$r" || { echo "B $rel M=$M r=$r FAILED sha lock"; exit 1; }
done
with_lock "$cell/.runlog.lockdir" append_run "$cell" "$r" "$seed" "$t0" "$tmp/stdout.log" || { echo "B $rel M=$M r=$r FAILED run.log lock"; exit 1; }
mkdir -p "$cell/.done_runs"; mv "$tmp" "$cell/.done_runs/run$r"      # the rest of the run folder (no trace left in it), kept
