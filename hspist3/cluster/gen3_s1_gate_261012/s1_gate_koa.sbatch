#!/bin/bash
# ##CHRIS 2026-10-12 (plan-author decision 39 item 5 c; registration 261012 sec. 4.8.10, gates G6 (c) and (d)): THE STAGE 1 GATE JOB on one
# node, with the CAMPAIGN binary (the f72d764 build of ~/harddisks_gen3_s1, checked against its build-time record) and the GATE binary
# (the a643bc7 build of ~/harddisks_gen3, untouched). PREPARED; typed only after the plan author has read its report. In order:
#   (c1) E0 (rule 4, the gen2 path byte-identical): cluster/resched_gate_261005/audit_runs_261007.py with the campaign binary against the
#        7b08827 reference of ~/harddisks_resched2: ctrl_min, ctrl_leg, afix, afix_leg; "audit vs plain" and "plain vs ref" must be
#        IDENTICAL; once with the default build and once with every command plus --engine=gen2 -- the a643bc7 gate job's part (a), unchanged
#        but for the binary;
#   (c2) the M1/M2 harness hashes: the campaign build's bin/gen3_m1 audit, gen3_m2 audit --quick, gen3_m2 audit; the SHA-256 of each output
#        against the a643bc7 gate job's KOA outputs (koa_harness_ref_a643bc7.sha256; job 15040861) -- the harness sources and edmd_core/ are
#        byte-identical between the commits, so equal files are expected;
#   (d)  the 12 gen3 lines of Test G chunk 1 (lines 1-25 of cluster/gen3_koa_261009/tasks_testG_koa_261009.txt; the even ones are gen3),
#        each run by Test G's worker once with the gate binary and once with the campaign binary, into two roots; s1_gate_testg_compare.py
#        compares every line's event hash and event counts (the run record [EDMD3-HEALTH]): equal in 12 of 12 is G6 (d).
# NO GIT AT RUN TIME (--exclude=cn-05-08-07). Every start-up check comes BEFORE the line "s1 gate on ...". Its outputs: the fixed root
# /mnt/lustre/koa/scratch/charing/harddisks_gen3/s1_gate (a job never writes into an existing root) and, at the end, a copy of this log there.
#   cd ~/harddisks_gen3_k2/hspist3 && sbatch cluster/gen3_s1_gate_261012/s1_gate_koa.sbatch
#SBATCH --job-name=gen3-s1-gate
#SBATCH --partition=shared
#SBATCH --account=uh
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --exclude=cn-05-08-07
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.out
set -uo pipefail
stop() { echo "STOP: [$(hostname)] $*"; exit 2; }
cd "${SLURM_SUBMIT_DIR:?submit with sbatch from ~/harddisks_gen3_k2/hspist3}"
P=cluster/gen3_s1_gate_261012
ROOT_NAME=$(basename "$(dirname "$SLURM_SUBMIT_DIR")")
[ "$ROOT_NAME" = harddisks_gen3_k2 ] || stop "the stage 1 gate runs from ~/harddisks_gen3_k2, not $ROOT_NAME"
sha256sum --status -c "$P/s1_gate_inputs.sha256" || stop "the gate inputs are not the delivered ones ($P/s1_gate_inputs.sha256)"
source cluster/koa_env.sh || stop "cluster/koa_env.sh failed"
S1="$HOME/harddisks_gen3_s1/hspist3"; REC="$S1/logs/BUILD_RECORD_S1.txt"
[ -f "$REC" ] || stop "no build-time record $REC (record_s1_build.sh)"
for f in 00ALLINONE bin/00ALLINONE_ld bin/gen3_m1 bin/gen3_m2; do
  want=$(awk -v f="$f" '$1 == "sha256" && $3 == f {print $2}' "$REC"); have=$(sha256sum "$S1/$f" | cut -d' ' -f1)
  [ -n "$want" ] && [ "$have" = "$want" ] || stop "$S1/$f is not the recorded campaign build ($have; record: '$want')"
done
V=$("$S1/00ALLINONE" --version | head -1); [ "$V" = "$(sed -n 's/^version //p' "$REC")" ] || stop "the campaign driver's version line is '$V'"
G="$HOME/harddisks_gen3/hspist3"
(cd "$G" && sha256sum --status -c "$SLURM_SUBMIT_DIR/cluster/gen3_koa_261009/k2_binaries_a643bc7.sha256") || stop "the binaries of ~/harddisks_gen3 are not runsheet step 2's a643bc7 builds"
VG=$("$G/00ALLINONE" --version | head -1); [ "$VG" = "00ALLINONE  git a643bc7  target koa" ] || stop "the gate driver's version line is '$VG'"
REF="$HOME/harddisks_resched2/hspist3"
[ -x "$REF/00ALLINONE" ] && (cd "$REF" && sha256sum --status -c logs/BUILD_KOA_LAST.hash) && "$REF/00ALLINONE" --version | head -1 | grep -q -- "git 7b08827  target koa" \
  || stop "the 7b08827 reference binary of ~/harddisks_resched2 is missing or not its recorded build"
OUT="/mnt/lustre/koa/scratch/charing/harddisks_gen3/s1_gate"
[ ! -e "$OUT" ] || stop "$OUT exists -- the gate job never writes into an existing root"
mkdir -p "$OUT" || stop "cannot create $OUT"
cp -p "$REC" "$OUT/BUILD_RECORD_S1.txt"              # the build-time record travels with the fetch (build_record_f72d764.txt is written from it)
echo "s1 gate on $(hostname) ($(lscpu | awk -F: '/Model name/{gsub(/^ +/,"",$2); print $2; exit}')): out $OUT"
echo "campaign: $V ($(sed -n 's/^commit //p' "$REC")); gate: $VG; E0 reference: $("$REF/00ALLINONE" --version | head -1)"
# (c1) E0
W="$OUT/bin_engine_gen2.sh"; printf '#!/bin/bash\nexec %s "$@" --engine=gen2\n' "$S1/00ALLINONE" > "$W"; chmod +x "$W"
for tag in default gen2flag; do
  b="$S1/00ALLINONE"; [ "$tag" = gen2flag ] && b="$W"
  python3 cluster/resched_gate_261005/audit_runs_261007.py run --bin "$b" --ref-bin "$REF/00ALLINONE" --out "$OUT/e0_$tag" --jobs 8 \
    --cases ctrl_min,ctrl_leg,afix,afix_leg || { echo "STOP: an E0 run failed ($tag)"; exit 1; }
  python3 cluster/resched_gate_261005/audit_runs_261007.py report --out "$OUT/e0_$tag" > "$OUT/e0_$tag/report.txt"
  echo "== (c1) E0 $tag:"; grep -E "^\| (ctrl_min|ctrl_leg|afix|afix_leg) \| git" "$OUT/e0_$tag/report.txt" | awk -F'|' '{print $2, "audit vs plain:", $(NF-2), "plain vs ref:", $(NF-1)}'
done
# (c2) the M1/M2 harness outputs (three runs at once) against the a643bc7 gate's KOA outputs
H="$OUT/harness"; mkdir -p "$H"
( cd "$H" && "$S1/bin/gen3_m1" audit > m1_audit_output.txt 2>&1; echo "m1 exit $?" > m1.exit ) &
( cd "$H" && "$S1/bin/gen3_m2" audit --quick > m2_audit_quick_output.txt 2>&1; echo "m2 quick exit $?" > m2q.exit ) &
( cd "$H" && "$S1/bin/gen3_m2" audit > m2_audit_output.txt 2>&1; echo "m2 exit $?" > m2.exit ) &
wait
echo "== (c2) the M1/M2 harness outputs of the campaign build ($(cat "$H"/m1.exit "$H"/m2q.exit "$H"/m2.exit | tr '\n' ';'))"
neq=0
for f in m1_audit_output.txt m2_audit_quick_output.txt m2_audit_output.txt; do
  want=$(awk -v f="$f" '$2 == f {print $1}' "$P/koa_harness_ref_a643bc7.sha256"); k=$(sha256sum "$H/$f" | cut -d' ' -f1)
  if [ "$k" = "$want" ]; then echo "   $f: $k -- EQUAL to the a643bc7 gate's KOA output"; neq=$((neq + 1)); else echo "   $f: $k -- DIFFERS from the a643bc7 gate's $want"; fi
done
echo "   harness outputs equal to the a643bc7 gate's: $neq of 3"
# (d) Test G chunk 1, its 12 gen3 lines, with both binaries
TL=cluster/gen3_koa_261009/tasks_testG_koa_261009.txt
sed -n '1,25p' "$TL" | awk '$11 == "gen3"' > "$OUT/testg_chunk1_gen3_lines.txt"
[ "$(wc -l < "$OUT/testg_chunk1_gen3_lines.txt")" -eq 12 ] || stop "Test G chunk 1 does not hold 12 gen3 lines"
for side in gate campaign; do
  B="$G/00ALLINONE"; [ "$side" = campaign ] && B="$S1/00ALLINONE"
  n=0
  while read -r line; do
    set -- $line
    HD_BIN="$B" HD_DATA="$OUT/testg_$side" HD_BUILD="$("$B" --version | head -1)" bash cluster/gen3_koa_261009/testG_worker_koa.sh "${@:1:11}" < /dev/null &
    n=$((n + 1)); [ $((n % 8)) -eq 0 ] && wait
  done < "$OUT/testg_chunk1_gen3_lines.txt"
  wait
done
echo "== (d) Test G chunk 1, the 12 gen3 lines, the gate binary against the campaign binary:"
python3 "$P/s1_gate_testg_compare.py" "$OUT/testg_gate" "$OUT/testg_campaign" "$OUT/testg_chunk1_gen3_lines.txt"
echo "s1 gate done; paste everything from 's1 gate on' to here"
cp -p "logs/${SLURM_JOB_NAME:-gen3-s1-gate}_${SLURM_JOB_ID}.out" "$OUT/job_log.out" 2>/dev/null || echo "(the job log was not copied into $OUT)"
