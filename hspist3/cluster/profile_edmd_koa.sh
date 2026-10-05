#!/bin/bash
# ##CHRIS 2026-10-05 (Task DD): engine cost profile on KOA -- ONE held-divider trajectory of the pi/8 H40 cell (N_s = 200)
# and ONE of the pi/8 H10 cell (N_s = 50), 500 sigma-time of record each (after the usual 200 sigma-time equilibration),
# with the installed recorded build (no build, no code change). Under `perf record -g` if perf exists and is allowed,
# else wall-clock only. Prints, per run: wall seconds, the divider and wall events of the event log per second (the
# binary keeps no total event count; disk-disk events are not logged), and the top 15 functions by self time.
#
# What the code says about the scaling (edmd_core/edmd.c, quoted in 261012 / runsheet step 10):
#   every disk-disk event:   resolve_ab(...); grid_build(S); schedule_for(S, e.a); schedule_for(S, e.b);   (:1529-1532)
#                            schedule_for loops over ALL partners: for (int j = 0; j < S->prm.N; ++j) schedule_ab(S, i, j);
#                            (:791-794) and grid_build loops over all N (:318) -> O(N) per event
#   every divider event:     if (e.type==EV_DL || e.type==EV_DR || ...) { reschedule_all_internal(S); }   (:1537-1539)
#                            reschedule_all_internal: for i < N, for j = i+1 .. N: schedule_ab(S, i, j);   (:812-816)
#                            -> O(N^2) per divider event, also for a HELD (mass 0) divider that cannot move
#   every event also drifts all N particles to the event time (:1504-1507) -> O(N)
# usage (from ~/harddisks/hspist3 on login-0101/0102):  sbatch cluster/profile_edmd_koa.sh
#SBATCH --job-name=profile-edmd
#SBATCH --partition=sandbox
#SBATCH --account=uh
#SBATCH --time=0:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.out
set -uo pipefail
cd "$SLURM_SUBMIT_DIR" || exit 1
source cluster/koa_env.sh || exit 1
export HD_BIN="$PWD/00ALLINONE"
[ -s logs/BUILD_KOA_LAST.hash ] && sha256sum --status -c logs/BUILD_KOA_LAST.hash || { echo "STOP: not the recorded build"; exit 1; }
"$HD_BIN" --version | head -1
OUT="/mnt/lustre/koa/scratch/charing/harddisks/profile_edmd_${SLURM_JOB_ID}"
mkdir "$OUT" || { echo "STOP: $OUT exists"; exit 1; }
PERF=$(command -v perf || true)
if [ -n "$PERF" ] && perf stat -e task-clock true >/dev/null 2>&1; then echo "perf: $PERF (usable)"; else echo "perf: not usable -- wall clock only"; PERF=""; fi
lscpu | awk -F: '/Model name/{gsub(/^ +/,"",$2); print "cpu: " $2; exit}'
for spec in "H40 10.000000 40.000000 200" "H10 10.000000 10.000000 50"; do
  set -- $spec; tag=$1 L0=$2 H=$3 NS=$4; d="$OUT/$tag"; mkdir -p "$d"; cd "$d" || exit 1
  CMD=("$HD_BIN" --mode=edmd --experiment=energy_transfer --headless --quiet --edmd-acc=0 --seed-drift-order=drift-first
       --energy-transfer-summary="$d/summary.csv" --energy-transfer-trace="$d/tr.csv" --trace-every=600
       --particles=$((2*NS)) --particles-boxes=$NS,$NS --particle-radius=0.5 --l0=$L0 --height=$H --num-walls=1
       --wall-positions=$L0 --wall-mass-factors=1000000000 --wall-thickness=0.05 --wall-thickness-vis=0.05
       --eff-output=wall-ke --wall-hold-steps=42000 --steps=60 --fixed-dt=0.4 --kbt1 --seed=9700)
  t0=$(date +%s.%N)
  if [ -n "$PERF" ]; then HD_PISTON_EVENTS="$d/ev.csv" perf record -g -o "$d/perf.data" -- "${CMD[@]}" > "$d/run.log" 2>&1
  else HD_PISTON_EVENTS="$d/ev.csv" "${CMD[@]}" > "$d/run.log" 2>&1; fi
  rc=$?; t1=$(date +%s.%N)
  echo; echo "== $tag (N_s = $NS, H = $H, L0 = $L0): exit $rc, wall $(awk -v a=$t0 -v b=$t1 'BEGIN{printf "%.1f", b-a}') s for 700 sigma-time"
  awk -F, -v a=$t0 -v b=$t1 'NR>1{n[$2]++; tot++} END{for(k in n) printf "   %-3s events %8d  (%.0f per wall-second)\n", k, n[k], n[k]/(b-a);
       printf "   wall+divider events %d (%.0f per wall-second)\n", tot, tot/(b-a)}' "$d/ev.csv"
  if [ -n "$PERF" ]; then
    echo "   top 15 functions by self time:"
    perf report -i "$d/perf.data" --stdio --no-children --sort symbol 2>/dev/null | grep -E '^ +[0-9.]+%' | head -15
  fi
  cd "$SLURM_SUBMIT_DIR" || exit 1
done
echo; echo "profile done; output in $OUT"
