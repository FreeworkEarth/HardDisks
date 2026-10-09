#!/bin/bash
# ##CHRIS 2026-10-09 (261012 sec. 4.7.20; stage G): PREPARED, NOT RUN. The KOA part of the gen-3 gate on one node, four parts in order:
#  (a) E0, gen2-path byte identity on KOA: cluster/resched_gate_261005/audit_runs_261007.py with the new binary (./00ALLINONE, the
#      gen3-containing koa build) against the 7b08827 reference binary of ~/harddisks_resched2 (Test T and T-prime's binary, read only,
#      checked against its recorded hash and version line): ctrl_min, ctrl_leg, afix, afix_leg, "plain vs ref" must be IDENTICAL;
#      once with the default build and once with every command plus --engine=gen2 (a wrapper);
#  (c) the KOA rate table of gen3: N = 100, 400, 900, 1600 x eta = pi/8, 0.70, 0.78; M4 seeding, exact box, 500 sigma-time held and
#      100 released, one run per cell, one at a time (events per second from the run records);
#  (d) the long-double spot check: bin/00ALLINONE_ld (koa-ld) against ./00ALLINONE on N = 100 at pi/8 (M = 300) and N = 400 at eta 0.70
#      (M = 500), 4 seeds each: health, the contact gap per time quantum (it must fall by about 2^11 with the long-double quantum),
#      the event counts and the argmax nu per run (information).
# Run from ~/harddisks_gen3/hspist3 after build_gen3_koa.sh (cluster/gen3_koa_261009/runsheet_gen3_koa_261009.txt, step 3):
#   sbatch cluster/gen3_koa_261009/g3_gate_koa.sbatch
#SBATCH --job-name=gen3-gate
#SBATCH --partition=sandbox
#SBATCH --account=uh
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.out
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:?submit with sbatch from ~/harddisks_gen3/hspist3}"
source cluster/koa_env.sh || exit 1
ROOT_NAME=$(basename "$(dirname "$SLURM_SUBMIT_DIR")")
[ "$ROOT_NAME" = harddisks_gen3 ] || { echo "STOP: the gen3 gate runs from ~/harddisks_gen3, not $ROOT_NAME"; exit 2; }
head=$(git rev-parse --short HEAD)
sha256sum --status -c logs/BUILD_KOA_LAST.hash && ./00ALLINONE --version | head -1 | grep -q -- "git $head  target koa" \
  || { echo "STOP: ./00ALLINONE is not the recorded clean koa build of HEAD $head"; exit 2; }
(cd bin && sha256sum --status -c ../logs/BUILD_KOA_LD.hash) && bin/00ALLINONE_ld --version | head -1 | grep -q -- "git $head  target koa-ld" \
  || { echo "STOP: bin/00ALLINONE_ld is not the recorded clean koa-ld build of HEAD $head"; exit 2; }
REF="$HOME/harddisks_resched2/hspist3"
[ -x "$REF/00ALLINONE" ] && (cd "$REF" && sha256sum --status -c logs/BUILD_KOA_LAST.hash) && "$REF/00ALLINONE" --version | head -1 | grep -q -- "git 7b08827  target koa" \
  || { echo "STOP: the 7b08827 reference binary of ~/harddisks_resched2 is missing or not its recorded build"; exit 2; }
OUT="/mnt/lustre/koa/scratch/charing/harddisks_gen3/gen3_gate_koa_${SLURM_JOB_ID}"; mkdir -p "$OUT"
echo "node $(hostname) | $(lscpu | awk -F: '/Model name/{gsub(/^ +/,"",$2); print $2; exit}') | out $OUT"
./00ALLINONE --version | head -1; bin/00ALLINONE_ld --version | head -1; "$REF/00ALLINONE" --version | head -1
# (a) E0
W="$OUT/bin_engine_gen2.sh"; printf '#!/bin/bash\nexec %s "$@" --engine=gen2\n' "$PWD/00ALLINONE" > "$W"; chmod +x "$W"
for tag in default gen2flag; do
  b="$PWD/00ALLINONE"; [ "$tag" = gen2flag ] && b="$W"
  python3 cluster/resched_gate_261005/audit_runs_261007.py run --bin "$b" --ref-bin "$REF/00ALLINONE" --out "$OUT/e0_$tag" --jobs 8 \
    --cases ctrl_min,ctrl_leg,afix,afix_leg || { echo "STOP: an E0 run failed ($tag)"; exit 1; }
  python3 cluster/resched_gate_261005/audit_runs_261007.py report --out "$OUT/e0_$tag" > "$OUT/e0_$tag/report.txt"
  echo "== (a) E0 $tag:"; grep -E "^\| (ctrl_min|ctrl_leg|afix|afix_leg) \| git" "$OUT/e0_$tag/report.txt" | awk -F'|' '{print $2, "audit vs plain:", $(NF-2), "plain vs ref:", $(NF-1)}'
  grep "^E2 (amended)" "$OUT/e0_$tag/report.txt"
done
# (c) rate table and (d) long-double spot check
python3 cluster/gen3_koa_261009/g3_rate_ld_koa.py --bin "$PWD/00ALLINONE" --bin-ld "$PWD/bin/00ALLINONE_ld" --out "$OUT" --jobs 8 \
  || { echo "STOP: the rate / long-double runs failed"; exit 1; }
echo "gen3 gate (KOA, one node) done; paste everything from 'node ' to here"
