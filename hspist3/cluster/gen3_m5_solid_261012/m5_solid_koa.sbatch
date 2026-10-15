#!/bin/bash
# ##CHRIS 2026-10-12 (plan-author decision 38 item 6; 261012 sec. 4.7.80): THE M5 SOLID-CELL JOB (stage 1 gate G6 b). PREPARED; typed
# only after the plan author has read its report (runsheet_m5_solid_261012.txt). One job, one node, 8 cores: the 8 runs of
# m5_solid_koa.py at once (N 400, eta 0.74, M 500; 4 seeds x the koa and koa-ld builds of a643bc7 in ~/harddisks_gen3, not rebuilt).
# NO GIT AT RUN TIME (decisions 26 and 28; node cn-05-08-07 has none, and it is excluded below). Every start-up check comes BEFORE the
# line "M5 solid cell on ...", so a job that stops in them ran nothing; every STOP line names its node: "STOP: [<node>] ...".
#   1. the clone: ~/harddisks_gen3_k2/hspist3 (the K2 clone, updated in place to this package's commit by the runsheet's block B);
#   2. the inputs: m5_solid_inputs.sha256 (sha256sum -c): the runner and every module of this repository it imports, koa_env.sh and
#      the binaries' hash file k2_binaries_a643bc7.sha256 are the delivered ones;
#   3. koa_env.sh (which no longer needs git);
#   4. the binaries, without git: ~/harddisks_gen3/hspist3/00ALLINONE and bin/00ALLINONE_ld against k2_binaries_a643bc7.sha256 (KOA's
#      runsheet step 2 build) AND the build-time files logs/BUILD_KOA_LAST.hash and logs/BUILD_KOA_LD.hash written by build_gen3_koa.sh,
#      and their version lines "00ALLINONE  git a643bc7  target koa" / "... target koa-ld";
#   5. the output folder .../harddisks_gen3/m5_solid_koa/job_<job id> does not exist (a job never writes into another job's folder).
# Its output: that folder (ld/N400_M500_eta0.74/seed<r>_<build>/: command.txt, run.log, the trace, the psi6(t) file) and, at the end,
# a copy of this job's log as job_log.out. The fetch (fetch_m5_solid_koa.sh, on the Mac) copies .../m5_solid_koa/ whole (about 10 MB).
#   cd ~/harddisks_gen3_k2/hspist3 && sbatch cluster/gen3_m5_solid_261012/m5_solid_koa.sbatch
#SBATCH --job-name=gen3-m5-solid
#SBATCH --partition=shared
#SBATCH --account=uh
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --exclude=cn-05-08-07
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.out
set -uo pipefail
stop() { echo "STOP: [$(hostname)] $*"; exit 2; }
cd "${SLURM_SUBMIT_DIR:?submit with sbatch from ~/harddisks_gen3_k2/hspist3}"
P=cluster/gen3_m5_solid_261012; K=cluster/gen3_koa_261009/k2_binaries_a643bc7.sha256
ROOT_NAME=$(basename "$(dirname "$SLURM_SUBMIT_DIR")")
[ "$ROOT_NAME" = harddisks_gen3_k2 ] || stop "the M5 job runs from ~/harddisks_gen3_k2, not $ROOT_NAME"
sha256sum --status -c "$P/m5_solid_inputs.sha256" || stop "the M5 inputs are not the delivered ones ($P/m5_solid_inputs.sha256)"
source cluster/koa_env.sh || stop "cluster/koa_env.sh failed"
A="$HOME/harddisks_gen3/hspist3"
(cd "$A" && sha256sum --status -c "$SLURM_SUBMIT_DIR/$K") || stop "the binaries of ~/harddisks_gen3 are not runsheet step 2's a643bc7 builds ($K)"
want=$(awk '$2 == "00ALLINONE" {print $1}' "$K"); want_ld=$(awk '$2 == "bin/00ALLINONE_ld" {print $1}' "$K")
[ -n "$want" ] && [ -n "$want_ld" ] || stop "$K lacks the 00ALLINONE or the bin/00ALLINONE_ld line"
grep -qx "$want  00ALLINONE" "$A/logs/BUILD_KOA_LAST.hash" || stop "the build-time file $A/logs/BUILD_KOA_LAST.hash does not hold $want"
grep -qx "$want_ld  00ALLINONE_ld" "$A/logs/BUILD_KOA_LD.hash" || stop "the build-time file $A/logs/BUILD_KOA_LD.hash does not hold $want_ld"
V=$("$A/00ALLINONE" --version | head -1); V_LD=$("$A/bin/00ALLINONE_ld" --version | head -1)
[ "$V" = "00ALLINONE  git a643bc7  target koa" ] || stop "the double build's version line is '$V'"
[ "$V_LD" = "00ALLINONE  git a643bc7  target koa-ld" ] || stop "the long-double build's version line is '$V_LD'"
OUT="/mnt/lustre/koa/scratch/charing/harddisks_gen3/m5_solid_koa/job_${SLURM_JOB_ID:?}"
[ ! -e "$OUT" ] || stop "$OUT exists -- a job never writes into another job's folder"
mkdir -p "$OUT" || stop "cannot create $OUT"
echo "M5 solid cell on $(hostname) ($(lscpu | awk -F: '/Model name/{gsub(/^ +/,"",$2); print $2; exit}')): out $OUT"
echo "binaries: $V ($want); $V_LD ($want_ld)"
python3 "$P/m5_solid_koa.py" --bin "$A/00ALLINONE" --bin-ld "$A/bin/00ALLINONE_ld" --out "$OUT" || stop "the M5 runs failed (the table above)"
echo "M5 solid cell done; paste everything from 'M5 solid cell on' to here"
cp "logs/${SLURM_JOB_NAME:-gen3-m5-solid}_${SLURM_JOB_ID}.out" "$OUT/job_log.out" 2>/dev/null || echo "(the job log was not copied into $OUT)"
