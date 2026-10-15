#!/bin/bash
# ##CHRIS 2026-10-12 (plan-author decision 39 item 5 c; gate G6 (c)): the compare job of the cross-node determinism test with the campaign
# binary. It is cluster/gen3_koa_261009/g3_xnode_compare.sbatch with THESE CHANGES: the scripts clone ~/harddisks_gen3_k2, no git, the
# inputs manifest, --exclude; the root is fixed (.../harddisks_gen3/s1_xnode), so no job id is typed; it shares the job name of the array
# (gen3-s1-xnode) and is submitted with --dependency=singleton, so it starts when every task has ended; it copies the task logs and its own
# log into the root. g3_xnode_compare.py is the a643bc7 job's, unchanged: the node pairs compared, IDENTICAL / DIFFERENT per file, and
# NOT TESTED if every task ran on one node.
#   cd ~/harddisks_gen3_k2/hspist3 && sbatch --dependency=singleton cluster/gen3_s1_gate_261012/s1_xnode_compare.sbatch
#SBATCH --job-name=gen3-s1-xnode
#SBATCH --partition=shared
#SBATCH --account=uh
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --exclude=cn-05-08-07
#SBATCH --output=logs/%x-compare_%j.out
#SBATCH --error=logs/%x-compare_%j.out
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:?submit with sbatch from ~/harddisks_gen3_k2/hspist3}"
P=cluster/gen3_s1_gate_261012
sha256sum --status -c "$P/s1_gate_inputs.sha256" || { echo "STOP: the gate inputs are not the delivered ones"; exit 2; }
source cluster/koa_env.sh || exit 1
R=/mnt/lustre/koa/scratch/charing/harddisks_gen3/s1_xnode
[ -d "$R" ] || { echo "STOP: no $R (the array ran nothing)"; exit 2; }
python3 cluster/gen3_koa_261009/g3_xnode_compare.py "$R"
mkdir -p "$R/job_logs"; cp -n -p logs/gen3-s1-xnode_*.out "$R/job_logs/" 2>/dev/null
echo "cross-node compare done; paste from 'cross-node determinism:' to here"
cp -p "logs/gen3-s1-xnode-compare_${SLURM_JOB_ID}.out" "$R/job_logs/" 2>/dev/null || true
