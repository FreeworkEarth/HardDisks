#!/bin/bash
# ##CHRIS 2026-10-05 (261012 sec. 4.4, engine gate G-E2): minimal vs legacy rescheduling, same binary, same seed, on the smoke
# cell (speed-of-sound, free divider after a 2000-step hold) and on one A-fixed cell (held divider, 5200 sigma-time).
# Four trajectories (six with the 279282b binary of ~/harddisks, read only) in parallel, then the comparison table (cluster/resched_gate_261005/ge2.py; registered rules in its header).
# Runs ONLY from the engine-branch clone (~/harddisks_resched/hspist3), after the smoke test built and checked the binary.
#   cd ~/harddisks_resched/hspist3 && sbatch cluster/resched_gate_261005/ge2.sbatch
#SBATCH --job-name=resched-ge2
#SBATCH --partition=sandbox
#SBATCH --account=uh
#SBATCH --time=0:30:00
#SBATCH --cpus-per-task=6
#SBATCH --mem=4G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.out
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:?submit with sbatch from ~/harddisks_resched/hspist3}"
source cluster/koa_env.sh || exit 1
ROOT_NAME=$(basename "$(dirname "$SLURM_SUBMIT_DIR")")
[ "$ROOT_NAME" != harddisks ] || { echo "STOP: the engine gate runs from the branch clone ~/harddisks_resched, not ~/harddisks"; exit 2; }
[ -s logs/BUILD_KOA_LAST.hash ] && sha256sum --status -c logs/BUILD_KOA_LAST.hash || { echo "STOP: ./00ALLINONE is not the recorded build -- run the smoke test first"; exit 1; }
./00ALLINONE --version | head -2
./00ALLINONE --version | head -1 | grep -q -- "git $(git rev-parse --short HEAD)  target koa" || { echo "STOP: binary is not the clean koa build of HEAD"; exit 2; }
OUT="/mnt/lustre/koa/scratch/charing/$ROOT_NAME/resched_gate_261005/ge2_${SLURM_JOB_ID}"
mkdir -p "$(dirname "$OUT")"
echo "node $(hostname) | $(lscpu | awk -F: '/Model name/{gsub(/^ +/,"",$2); print $2; exit}') | out $OUT"
# the recorded 279282b binary of ~/harddisks (read only): its legacy-path comparison is reported, not gated
OLD="$HOME/harddisks/hspist3"; OLDARG=()
if [ -x "$OLD/00ALLINONE" ] && (cd "$OLD" && sha256sum --status -c logs/BUILD_KOA_LAST.hash) && "$OLD/00ALLINONE" --version | head -1 | grep -q -- "git 279282b  target koa"; then
  OLDARG=(--bin-old "$OLD/00ALLINONE"); echo "279282b binary: $("$OLD/00ALLINONE" --version | head -1) (recorded hash OK)"
else echo "279282b binary not found or not the recorded build -- the legacy-vs-279282b comparison is skipped"; fi
python3 cluster/resched_gate_261005/ge2.py run --bin "$PWD/00ALLINONE" ${OLDARG[@]+"${OLDARG[@]}"} --out "$OUT" || { echo "STOP: a G-E2 run failed"; exit 1; }
python3 cluster/resched_gate_261005/ge2.py compare --out "$OUT" | tee "$OUT/ge2_compare.txt"
