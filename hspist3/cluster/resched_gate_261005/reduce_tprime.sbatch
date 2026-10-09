#!/bin/bash
# ##CHRIS 2026-10-08 HST (2026-10-09 on the plan author's clock; 261012 sec. 4.4.14): the dependent job after the T-prime chunk
# array (testTprime_chunks.sbatch). One process, one core:
#   1. MERGE the 32 chunk roots ($SCRATCH/harddisks_resched2/tprime_chunks_261009/c01..c32) into the registered layout
#      ($SCRATCH/harddisks_resched2/hspist3/experiments_speed_of_sound/.../resched_testTprime_261007/<policy>/epi8_H_H10_L10/m_<M>),
#      the layout the registered shape would have written and fetch_resched2.sh and validation/resched_testTprime_261007.py read:
#      traces copied (a name already there must be byte-identical, else STOP); each chunk cell's run.log appended up to its
#      current end -- the marker .merged_into_registered_layout in the chunk cell records the bytes already appended, so a rerun
#      of this job appends only what a rerun chunk added, never a section twice; .failed_run* and .stale_run* directories
#      copied; .build_git checked against this binary; each registered cell it fills gets .merged_from_chunks_261009, so a
#      rerun chunk can tell merged data from data of the registered shape. Nothing is moved or deleted.
#   2. reduce_B per policy and mass (cluster/confinement_20261013/reduce_B.py <cell> <M>), as testTprime.sbatch did.
#   3. Completeness: traces, run.log sections and red_nu.csv rows per policy and mass (400 each expected).
# Same STOP checks as testTprime.sbatch (clone, recorded binary hash, version line, task-list hash).
#   sbatch --dependency=afterok:<array id> cluster/resched_gate_261005/reduce_tprime.sbatch
# DRY RUN (Mac plumbing test only): the TPRIME_DRY_* variables of testTprime_chunks.sbatch; chunk roots that were not run are
# skipped instead of stopping the job, and the expected count is the number of the dry list's lines in the chunks that ran.
#SBATCH --job-name=resched-testTprime-reduce
#SBATCH --partition=shared
#SBATCH --account=uh
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.out
set -uo pipefail
SCRATCH="/mnt/lustre/koa/scratch/charing"
CHUNK=50; NCHUNK=32
REG=cluster/resched_gate_261005/tasks_Tprime_epi8_H_H10_L10.txt
REL_TP=experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/resched_testTprime_261007
cd "${SLURM_SUBMIT_DIR:?submit with sbatch from ~/harddisks_resched2/hspist3}"
ROOT_NAME=$(basename "$(dirname "$SLURM_SUBMIT_DIR")")
[ "$ROOT_NAME" = harddisks_resched2 ] || { echo "STOP: Test T-prime runs from the gate-v2 clone (~/harddisks_resched2), not $ROOT_NAME"; exit 2; }
DRY=0
if [ -n "${TPRIME_DRY_TASKS:-}" ]; then
  : "${TPRIME_DRY_BIN:?dry run: set TPRIME_DRY_BIN}"; : "${TPRIME_DRY_ROOT:?dry run: set TPRIME_DRY_ROOT}"
  case "$TPRIME_DRY_TASKS$TPRIME_DRY_BIN$TPRIME_DRY_ROOT" in *"/mnt/lustre/"*) echo "STOP: the dry run never runs on KOA"; exit 2 ;; esac
  echo "DRY RUN -- Mac plumbing test, NOT the registered test; outputs are not analysed"
  DRY=1; TASKS=$TPRIME_DRY_TASKS; BASE=$TPRIME_DRY_ROOT; HD_BIN="$TPRIME_DRY_BIN"
else
  source cluster/koa_env.sh || exit 1
  BASE="$SCRATCH/$ROOT_NAME"; HD_BIN="$PWD/00ALLINONE"
  [ -s logs/BUILD_KOA_LAST.hash ] || { echo "STOP: no logs/BUILD_KOA_LAST.hash"; exit 1; }
  sha256sum --status -c logs/BUILD_KOA_LAST.hash || { echo "STOP: ./00ALLINONE is not the build recorded in logs/BUILD_KOA_LAST.hash"; exit 1; }
  TASKS=$REG
fi
HD_BUILD="$("$HD_BIN" --version | head -1)"
if [ "$DRY" = 0 ]; then
  [ "$HD_BUILD" = "00ALLINONE  git 7b08827  target koa" ] || { echo "STOP: not the Test T binary (git 7b08827, target koa): $HD_BUILD"; exit 1; }
  [ "$(sha256sum "$TASKS" | cut -d' ' -f1)" = cdef566b2a74fd0ec9122eee2ed2fef271fa9ee75d9630de378819565831992f ] || { echo "STOP: task list differs from the registered one"; exit 2; }
fi
echo "binary: $HD_BUILD; task list $TASKS ($(wc -l < "$TASKS") lines)"
CANON="$BASE/hspist3"; CHROOT="$BASE/tprime_chunks_261009"
# the registered layout's data root carries the build generation that wrote it (Test T, the same binary)
if [ -e "$CANON/.build_generation" ]; then
  [ "$(cat "$CANON/.build_generation")" = "$HD_BUILD" ] || { echo "STOP: $CANON/.build_generation is '$(cat "$CANON/.build_generation")', not '$HD_BUILD'"; exit 2; }
elif [ "$DRY" = 1 ]; then mkdir -p "$CANON" && printf '%s\n' "$HD_BUILD" > "$CANON/.build_generation"
else echo "STOP: $CANON/.build_generation missing (Test T wrote it there)"; exit 2; fi

# ---- 1. merge
fail=0; ran=()
for i in $(seq 1 $NCHUNK); do
  C="$CHROOT/c$(printf %02d "$i")"
  if [ ! -d "$C" ]; then
    if [ "$DRY" = 1 ]; then continue; fi
    echo "STOP: chunk root $C missing (did array task $i run?)"; exit 2
  fi
  ran+=("$i")
  [ "$(cat "$C/.build_generation" 2>/dev/null)" = "$HD_BUILD" ] || { echo "STOP: $C/.build_generation is not '$HD_BUILD'"; exit 2; }
  for src in "$C/$REL_TP"/*/epi8_H_H10_L10/m_*; do
    [ -d "$src" ] || continue
    dst="$CANON/${src#"$C/"}"
    [ "$(cat "$src/.build_git" 2>/dev/null)" = "$HD_BUILD" ] || { echo "STOP: $src/.build_git is not '$HD_BUILD'"; exit 2; }
    mkdir -p "$dst"
    if [ -e "$dst/.build_git" ]; then [ "$(cat "$dst/.build_git")" = "$HD_BUILD" ] || { echo "STOP: $dst/.build_git differs"; exit 2; }
    else cp -p "$src/.build_git" "$dst/.build_git"; fi
    for f in "$src"/wall_x_positions_L0_*_run*.csv; do
      [ -e "$f" ] || continue
      g="$dst/$(basename "$f")"
      if [ -e "$g" ]; then cmp -s "$f" "$g" || { echo "STOP: $g exists and differs from $f"; exit 2; }
      else cp -p "$f" "$g"; fi
    done
    for d in "$src"/.failed_run* "$src"/.stale_run*; do
      [ -d "$d" ] || continue
      [ -e "$dst/$(basename "$d")" ] || cp -Rp "$d" "$dst/"
    done
    printf '%s\n' "$C" >> "$dst/.merged_from_chunks_261009"
    if [ -e "$src/run.log" ]; then
      done_b=$(awk '{print $1}' "$src/.merged_into_registered_layout" 2>/dev/null); done_b=${done_b:-0}
      size_b=$(wc -c < "$src/run.log" | tr -d ' ')
      [ "$size_b" -ge "$done_b" ] || { echo "STOP: $src/run.log is shorter than the part already appended"; exit 2; }
      if [ "$size_b" -gt "$done_b" ]; then
        tail -c +$((done_b + 1)) "$src/run.log" >> "$dst/run.log" || { echo "STOP: appending $src/run.log failed"; exit 2; }
        printf '%s bytes of run.log appended to %s, last at %s\n' "$size_b" "$dst/run.log" "$(date '+%Y-%m-%d %H:%M:%S')" > "$src/.merged_into_registered_layout"
      fi
    fi
  done
done
echo "merged chunk roots: ${ran[*]:-none}"

# ---- 2. reductions, per policy and mass
for pol in minimal legacy; do
  for M in 300 1500; do
    [ -d "$CANON/$REL_TP/$pol/epi8_H_H10_L10/m_$M" ] || continue
    python3 cluster/confinement_20261013/reduce_B.py "$CANON/$REL_TP/$pol/epi8_H_H10_L10" "$M" || { echo "reduce_tprime $pol M=$M reduce_B FAILED (exit $?)"; fail=$((fail + 1)); }
  done
done

# ---- 3. completeness: traces, run.log sections, red_nu.csv rows; expected from the task list (the chunks that ran, for the dry run)
for pol in minimal legacy; do
  for M in 300 1500; do
    exp=0
    for i in ${ran[@]+"${ran[@]}"}; do
      a=$((CHUNK * (i - 1) + 1)); b=$((CHUNK * i))
      exp=$((exp + $(sed -n "${a},${b}p" "$TASKS" | awk -v p="$pol" -v m="$M" '$11 == p && $3 == m' | wc -l)))
    done
    d="$CANON/$REL_TP/$pol/epi8_H_H10_L10/m_$M"
    nt=$(ls "$d"/wall_x_positions_L0_*_run*.csv 2>/dev/null | wc -l)
    ns=$(grep -c '^##RUN' "$d/run.log" 2>/dev/null || true)
    nr=$(( $(wc -l < "$d/red_nu.csv" 2>/dev/null || echo 1) - 1 ))
    ok=yes; { [ "$nt" -eq "$exp" ] && [ "${ns:-0}" -eq "$exp" ] && [ "$nr" -eq "$exp" ]; } || { ok=NO; fail=$((fail + 1)); }
    echo "$pol M=$M: traces $nt, run.log sections ${ns:-0}, red_nu.csv rows $nr, expected $exp -> $ok"
  done
done
[ "$fail" -eq 0 ] || echo "reduce_tprime: $fail check(s) FAILED"
echo "reduce_tprime done; failures: $fail"
