#!/bin/bash
# ##CHRIS 2026-10-10 (plan-author decision 20, tasks 4-5; 261012 sec. 4.7.47): RUNSHEET ADDENDUM STEP K2-2 -- parts (c) and (d) of the
# gen3 gate (g3_gate_koa.sbatch), which stopped in job 15040861 on the parsing bug of g3_rate_ld_koa.py that this commit fixes:
#  (c) the KOA rate table of gen3 (12 runs, one at a time) and (d) the long-double spot check (16 runs, at most 8 at a time), with the
#      same script, arguments and seeds as the gate job.
# It runs from the SECOND clone ~/harddisks_gen3_k2 (the scripts of this commit) and uses the binaries of ~/harddisks_gen3 (a643bc7;
# not rebuilt, not touched). Before anything runs: all four of them against the SHA-256 Chris pasted from runsheet step 2
# (k2_binaries_a643bc7.sha256), and the two drivers' version lines against "git a643bc7  target koa" / "koa-ld"; any mismatch: STOP.
# 1 node, 8 cores, 8 GB, at most 3 h on shared, as the gate job. Output: $SCRATCH/harddisks_gen3/gen3_rate_ld_k2_<job id>/.
# Rehearsed end to end on the Mac (rehearse_k2_mac.py) before Chris types it.
#   cd ~/harddisks_gen3_k2/hspist3 && sbatch cluster/gen3_koa_261009/g3_rate_ld_k2.sbatch
#SBATCH --job-name=gen3-rate-ld-k2
#SBATCH --partition=shared
#SBATCH --account=uh
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.out
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:?submit with sbatch from ~/harddisks_gen3_k2/hspist3}"
source cluster/koa_env.sh || exit 1
ROOT_NAME=$(basename "$(dirname "$SLURM_SUBMIT_DIR")")
[ "$ROOT_NAME" = harddisks_gen3_k2 ] || { echo "STOP: the K2 jobs run from ~/harddisks_gen3_k2, not $ROOT_NAME"; exit 2; }
B="$HOME/harddisks_gen3/hspist3"
(cd "$B" && sha256sum --status -c "$SLURM_SUBMIT_DIR/cluster/gen3_koa_261009/k2_binaries_a643bc7.sha256") \
  || { echo "STOP: the binaries of ~/harddisks_gen3 are not the a643bc7 builds of runsheet step 2 (k2_binaries_a643bc7.sha256)"; exit 2; }
[ "$("$B/00ALLINONE" --version | head -1)" = "00ALLINONE  git a643bc7  target koa" ] \
  && [ "$("$B/bin/00ALLINONE_ld" --version | head -1)" = "00ALLINONE  git a643bc7  target koa-ld" ] \
  || { echo "STOP: the drivers of ~/harddisks_gen3 do not print the a643bc7 koa / koa-ld version lines"; exit 2; }
OUT="/mnt/lustre/koa/scratch/charing/harddisks_gen3/gen3_rate_ld_k2_${SLURM_JOB_ID}"
[ ! -e "$OUT" ] || { echo "STOP: $OUT exists -- a job is never repeated into the same folder"; exit 2; }
mkdir -p "$OUT"
echo "node $(hostname) | $(lscpu | awk -F: '/Model name/{gsub(/^ +/,"",$2); print $2; exit}') | out $OUT | scripts: clone HEAD $(git rev-parse HEAD) | binaries: $B"
"$B/00ALLINONE" --version | head -1; "$B/bin/00ALLINONE_ld" --version | head -1
echo "binaries: SHA-256 as recorded in runsheet step 2 (all four)"
python3 cluster/gen3_koa_261009/g3_rate_ld_koa.py --bin "$B/00ALLINONE" --bin-ld "$B/bin/00ALLINONE_ld" --out "$OUT" --jobs 8 \
  || { echo "STOP: the rate / long-double runs failed"; exit 1; }
echo "K2-2 (gate parts (c) and (d)) done; paste everything from 'node ' to here"
