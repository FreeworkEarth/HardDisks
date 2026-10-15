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
#
# ##CHRIS 2026-10-05 (engine gate G-E5, 261012 sec. 4.4; branch engine-divider-resched): the same runs ONCE PER RESCHEDULING
# POLICY with the same binary on the same node -- default (minimal: after a divider event only that disk, plus every disk's
# divider event if the divider moved) and --legacy-resched (the old full reschedule) -- and a second kind, "free": held for
# 200 sigma-time, then released with mass alpha 2 N_s, alpha = 5 (M = 500 / 2000), and 500 sigma-time free, so the
# free-divider path of the fix (one O(N) pass over all disks' divider events per divider collision) is timed too.
# Exponent per policy and kind: p = ln(t_H40 / t_H10) / ln 4 (N = 100 -> 400). Factor = t_legacy / t_minimal, same node.
# The data root follows the clone: $SCRATCH/<clone name>/profile_edmd_<job>.
# usage (from ~/<clone>/hspist3 on login-0101/0102):  sbatch cluster/profile_edmd_koa.sh
#SBATCH --job-name=profile-edmd
#SBATCH --partition=sandbox
#SBATCH --account=uh
#SBATCH --time=0:40:00
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
ROOT_NAME=$(basename "$(dirname "$SLURM_SUBMIT_DIR")")
OUT="/mnt/lustre/koa/scratch/charing/$ROOT_NAME/profile_edmd_${SLURM_JOB_ID}"
mkdir -p "$(dirname "$OUT")"; mkdir "$OUT" || { echo "STOP: $OUT exists"; exit 1; }
PERF=$(command -v perf || true)
if [ -n "$PERF" ] && perf stat -e task-clock true >/dev/null 2>&1; then echo "perf: $PERF (usable)"; else echo "perf: not usable -- wall clock only"; PERF=""; fi
lscpu | awk -F: '/Model name/{gsub(/^ +/,"",$2); print "cpu: " $2; exit}'; echo "node: $(hostname)"
printf 'kind\ttag\tN\tpolicy\twall_s\tdiv_events\twall_events\texit\n' > "$OUT/times.tsv"
# ##CHRIS 2026-10-07 (gate v2, E4): a third kind, "dense": a free divider (alpha = 5) near eta = 0.70 (L0 = 269/48 = 5.604167 on
# the 1/48 grid; H = 10 and 40 at the same N_s/H, so the same eta_true 0.70075), 100 sigma-time held + 100 free (shorter: many more
# collisions per sigma-time). Its N = 100 -> 400 gain decides the melting-sweep cost.
for combo in "held H40 10.000000 40.000000 200" "held H10 10.000000 10.000000 50" "free H40 10.000000 40.000000 200" \
             "free H10 10.000000 10.000000 50" "dense D40 5.604167 40.000000 200" "dense D10 5.604167 10.000000 50"; do
  for policy in minimal legacy; do
   set -- $combo; kind=$1 tag=$2 L0=$3 H=$4 NS=$5; d="$OUT/${kind}_${tag}_${policy}"; mkdir -p "$d"; cd "$d" || exit 1
   case "$kind" in held) HOLD=42000 POST=60 MF=1000000000 ;; free) HOLD=12000 POST=30000 MF=$((5 * 2 * NS)) ;; dense) HOLD=6000 POST=6000 MF=$((5 * 2 * NS)) ;; esac
   EXTRA=(); [ "$policy" = legacy ] && EXTRA=(--legacy-resched)
   CMD=("$HD_BIN" --mode=edmd --experiment=energy_transfer --headless --quiet --edmd-acc=0 --seed-drift-order=drift-first
        --energy-transfer-summary="$d/summary.csv" --energy-transfer-trace="$d/tr.csv" --trace-every=600
        --particles=$((2*NS)) --particles-boxes=$NS,$NS --particle-radius=0.5 --l0=$L0 --height=$H --num-walls=1
        --wall-positions=$L0 --wall-mass-factors=$MF --wall-thickness=0.05 --wall-thickness-vis=0.05
        --eff-output=wall-ke --wall-hold-steps=$HOLD --steps=$POST --fixed-dt=0.4 --kbt1 --seed=9700 ${EXTRA[@]+"${EXTRA[@]}"})
   t0=$(date +%s.%N)
   if [ -n "$PERF" ]; then HD_PISTON_EVENTS="$d/ev.csv" perf record -g -o "$d/perf.data" -- "${CMD[@]}" > "$d/run.log" 2>&1
   else HD_PISTON_EVENTS="$d/ev.csv" "${CMD[@]}" > "$d/run.log" 2>&1; fi
   rc=$?; t1=$(date +%s.%N)
   w=$(awk -v a=$t0 -v b=$t1 'BEGIN{printf "%.2f", b-a}')
   nd=$(awk -F, 'NR>1 && $2 ~ /^D/{n++} END{print n+0}' "$d/ev.csv"); nw=$(awk -F, 'NR>1 && $2 ~ /^W/{n++} END{print n+0}' "$d/ev.csv")
   printf '%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\n' "$kind" "$tag" $((2*NS)) "$policy" "$w" "$nd" "$nw" "$rc" >> "$OUT/times.tsv"
   echo; echo "== $kind $tag (N_s = $NS, H = $H, L0 = $L0) $policy: exit $rc, wall $w s; divider events $nd, outer-wall events $nw"
   grep -h "EDMD-RESCHED\|EDMD-HEALTH" "$d/run.log" | sed 's/^/   /'
   if [ -n "$PERF" ]; then
     echo "   top 15 functions by self time:"
     perf report -i "$d/perf.data" --stdio --no-children --sort symbol 2>/dev/null | grep -E '^ +[0-9.]+%' | head -15
   fi
   cd "$SLURM_SUBMIT_DIR" || exit 1
  done
done
echo; echo "| kind | N = 100 minimal [s] | N = 100 legacy [s] | N = 400 minimal [s] | N = 400 legacy [s] | factor N=100 | factor N=400 | p minimal | p legacy |"
echo "|---|---|---|---|---|---|---|---|---|"
awk -F'\t' 'NR>1{t[$1" "$3" "$4]=$5; k[$1]=1; if ($8 != 0) bad[$1]=bad[$1] " " $3 "/" $4 "(exit " $8 ")"}
  END{for (q in k){a=t[q" 100 minimal"]; b=t[q" 100 legacy"]; c=t[q" 400 minimal"]; e=t[q" 400 legacy"];
      printf "| %s | %.2f | %.2f | %.2f | %.2f | %.2f | %.2f | %.2f | %.2f |\n", q, a, b, c, e, b/a, e/c, log(c/a)/log(4), log(e/b)/log(4);
      if (q in bad) printf "  **INVALID row %s: a run exited non-zero:%s -- its time is not a profile**\n", q, bad[q]}}' "$OUT/times.tsv"
echo; echo "profile done; output in $OUT"
