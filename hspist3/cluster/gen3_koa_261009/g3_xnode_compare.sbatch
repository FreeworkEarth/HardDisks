#!/bin/bash
# ##CHRIS 2026-10-09 (261012 sec. 4.7.29; plan-author decision 13, item 2d): PREPARED, NOT RUN. The compare job of the cross-node
# determinism test: after every task of the g3_det_xnode.sbatch array has ended, one core reads the array's folder and prints the node
# pairs compared, IDENTICAL / DIFFERENT per file, and NOT TESTED if every task ran on one node (g3_xnode_compare.py). Read-only.
# Run from ~/harddisks_gen3/hspist3 right after submitting the array (runsheet step 4), with the array's job id twice:
#   sbatch --dependency=afterany:<array job id> cluster/gen3_koa_261009/g3_xnode_compare.sbatch <array job id>
#SBATCH --job-name=gen3-xnode-compare
#SBATCH --partition=shared
#SBATCH --account=uh
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.out
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:?submit with sbatch from ~/harddisks_gen3/hspist3}"
source cluster/koa_env.sh || exit 1
AJ=${1:?give the array job id as the argument}
case "$AJ" in *[!0-9]*) echo "STOP: the argument must be the array job id (digits only), not '$AJ'"; exit 2;; esac
python3 cluster/gen3_koa_261009/g3_xnode_compare.py "/mnt/lustre/koa/scratch/charing/harddisks_gen3/gen3_xnode_${AJ}"
echo "cross-node compare done; paste from 'cross-node determinism:' to here"
