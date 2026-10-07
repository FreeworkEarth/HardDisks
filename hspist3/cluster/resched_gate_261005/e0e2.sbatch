#!/bin/bash
# ##CHRIS 2026-10-07 (261012 sec. 4.4.10, gate version 2, E0 and E2): on the KOA build of the gate-v2 commit (73fc07f plus the
# --resched-audit switch), run cluster/resched_gate_261005/audit_runs_261007.py:
#   E0  every case also without the audit ("plain") and with the 73fc07f binary of ~/harddisks_resched ("ref", read only, checked
#       against its recorded hash and version line): plain vs ref must be byte-identical (ctrl_min = the G-E2 smoke trajectory
#       on the minimal path, ctrl_leg = the same on the legacy path, afix / afix_leg = the G-E2 A-fixed trajectory on both paths);
#   E2  (judged by the AMENDED rule of the second plan-author decision of 2026-10-07, item 1; audit_runs_261007.py e2_amended)
#       the contact audit at rounding level in every class and run, and the schedule-equivalence audit: free divider lightest
#       (free_M50_long, >= 1e5 divider events) and heaviest mass (free_M2000), held divider (afix, >= 1e4 D0 events), a dense state
#       near eta = 0.70 (dense_M50, dense_M2000); the mode-2 controls ctrl_min/ctrl_leg as information.
# Runs ONLY from the gate-v2 clone (~/harddisks_resched2/hspist3), after its smoke test built and checked the binary.
#   cd ~/harddisks_resched2/hspist3 && sbatch cluster/resched_gate_261005/e0e2.sbatch
#SBATCH --job-name=resched-e0e2
#SBATCH --partition=sandbox
#SBATCH --account=uh
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.out
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:?submit with sbatch from ~/harddisks_resched2/hspist3}"
source cluster/koa_env.sh || exit 1
ROOT_NAME=$(basename "$(dirname "$SLURM_SUBMIT_DIR")")
case "$ROOT_NAME" in harddisks|harddisks_resched) echo "STOP: gate v2 runs from its own clone (~/harddisks_resched2)"; exit 2 ;; esac
[ -s logs/BUILD_KOA_LAST.hash ] && sha256sum --status -c logs/BUILD_KOA_LAST.hash || { echo "STOP: ./00ALLINONE is not the recorded build -- run the smoke test first"; exit 1; }
./00ALLINONE --version | head -1 | grep -q -- "git $(git rev-parse --short HEAD)  target koa" || { echo "STOP: binary is not the clean koa build of HEAD"; exit 2; }
REF="$HOME/harddisks_resched/hspist3"
[ -x "$REF/00ALLINONE" ] && (cd "$REF" && sha256sum --status -c logs/BUILD_KOA_LAST.hash) && "$REF/00ALLINONE" --version | head -1 | grep -q -- "git 73fc07f  target koa" \
  || { echo "STOP: the 73fc07f reference binary of ~/harddisks_resched is missing or not its recorded build"; exit 2; }
./00ALLINONE --version | head -1; "$REF/00ALLINONE" --version | head -1
OUT="/mnt/lustre/koa/scratch/charing/$ROOT_NAME/resched_gate_261005/e0e2_${SLURM_JOB_ID}"
mkdir -p "$(dirname "$OUT")"
echo "node $(hostname) | $(lscpu | awk -F: '/Model name/{gsub(/^ +/,"",$2); print $2; exit}') | out $OUT"
python3 cluster/resched_gate_261005/audit_runs_261007.py run --bin "$PWD/00ALLINONE" --ref-bin "$REF/00ALLINONE" --out "$OUT" --jobs 8 \
  --cases free_M50_long,free_M2000,afix,afix_leg,dense_M50,dense_M2000,ctrl_min,ctrl_leg || { echo "STOP: a run failed"; exit 1; }
python3 cluster/resched_gate_261005/audit_runs_261007.py report --out "$OUT" > "$OUT/report.txt"
# what to paste: the E2 section as amended (2026-10-07, item 1) with the E2 and E0 verdict lines
sed -n '/### E2 as amended/,/^E2 (amended)/p' "$OUT/report.txt"
