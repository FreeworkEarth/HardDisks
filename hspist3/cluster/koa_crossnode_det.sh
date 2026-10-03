#!/bin/bash
# ##CHRIS 2026-10-02 (Task K2): cross-node determinism check. The same binary and the same seed on TWO DIFFERENT nodes
# must give byte-identical output -- that is what the portable -march=x86-64-v2 build is for (Makefile CFLAGS_KOA). The
# smoke test runs its two determinism steps inside one single-node job, i.e. on the same node; this job forces two nodes.
# Run it AFTER the smoke test has built and checked the binary (it does not build). Cost: 2 cores for ~1 minute.
#   cd ~/harddisks/hspist3 && sbatch cluster/koa_crossnode_det.sh
#SBATCH --job-name=det-xnode
#SBATCH --partition=sandbox
#SBATCH --account=uh
#SBATCH --time=00:15:00
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.out
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:?submit with sbatch from ~/harddisks/hspist3}"
source cluster/koa_env.sh || exit 1
[ -x ./00ALLINONE ] || { echo "STOP: no ./00ALLINONE -- run the smoke test first (it builds)"; exit 2; }
./00ALLINONE --version | head -2
./00ALLINONE --version | head -1 | grep -q -- "git $(git rev-parse --short HEAD)  target koa" || { echo "STOP: binary is not the clean koa build of HEAD"; exit 2; }
OUT="/mnt/lustre/koa/scratch/charing/harddisks/hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/confinement_pilot_20261013/koa_det_xnode_${SLURM_JOB_ID}"
echo "nodes: $SLURM_JOB_NODELIST"
srun --nodes=1 --ntasks=1 --relative=0 --exact python3 cluster/confinement_pilot.py det1 --tag A --bin ./00ALLINONE --out "$OUT" &
srun --nodes=1 --ntasks=1 --relative=1 --exact python3 cluster/confinement_pilot.py det1 --tag B --bin ./00ALLINONE --out "$OUT" &
wait
python3 cluster/confinement_pilot.py detcmp --out "$OUT"
