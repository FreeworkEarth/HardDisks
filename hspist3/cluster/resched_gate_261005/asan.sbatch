#!/bin/bash
# ##CHRIS 2026-10-07 (261012 sec. 4.4.13, gate v3, third plan-author decision of 2026-10-07, item 3): ASan + UBSan on KOA.
# A SCRATCH build of 7b08827 -- gcc -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -- compiled from this clone's
# sources after checking that every compiled file equals commit 7b08827 (git diff) and that kissfft is the pinned commit, and
# written to $SCRATCH/<clone>/asan_261007_<jobid>/ only: ./00ALLINONE and logs/BUILD_KOA_LAST.hash are never touched (their hash is
# checked before and after). Its version line is "00ALLINONE  git 7b08827  target asan-scratch", so it cannot be taken for the
# recorded binary. Then asan_runs_261007.py runs four trajectories with --resched-audit: the smoke trajectory (M = 50) and the
# first M = 300 trajectory of Test T, each on both policies; LeakSanitizer on. Requirement: zero sanitizer reports, exit 0, audit
# missing = extra = 0 in all four; any report = stop and paste. Measured on the Mac (UBSan -O1 + audit, the same four runs):
# 2 s and 30-33 s per run; KOA with ASan a few minutes -> 1 h on the sandbox partition, 4 cores (the four runs in parallel).
#   cd ~/harddisks_resched2/hspist3 && sbatch cluster/resched_gate_261005/asan.sbatch
#SBATCH --job-name=resched-asan
#SBATCH --partition=sandbox
#SBATCH --account=uh
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.out
set -uo pipefail
SCRATCH="/mnt/lustre/koa/scratch/charing"
cd "${SLURM_SUBMIT_DIR:?submit with sbatch from ~/harddisks_resched2/hspist3}"
source cluster/koa_env.sh || exit 1
ROOT_NAME=$(basename "$(dirname "$SLURM_SUBMIT_DIR")")
[ "$ROOT_NAME" = harddisks_resched2 ] || { echo "STOP: the ASan job runs from the gate-v2 clone (~/harddisks_resched2), not $ROOT_NAME"; exit 2; }
OUT="$SCRATCH/$ROOT_NAME/asan_261007_${SLURM_JOB_ID}"
[ -e "$OUT" ] && { echo "STOP: $OUT exists -- not overwriting"; exit 2; }
mkdir -p "$OUT"
REC_BEFORE=$(sha256sum 00ALLINONE logs/BUILD_KOA_LAST.hash)
SRCS="00ALLINONE.c kissfft/kiss_fft.c kissfft/kiss_fftr.c edmd_core/edmd.c edmd_core/edmd_accelerated.c experiment_validation.c"
CHECK="00ALLINONE.c experiment_validation.c experiment_validation.h edmd_core"
git diff --quiet 7b08827 -- $CHECK || { git diff --stat 7b08827 -- $CHECK; echo "STOP: the sources differ from 7b08827"; exit 2; }
[ "$(git -C kissfft rev-parse HEAD)" = febd4caeed32e33ad8b2e0bb5ea77542c40f18ec ] || { echo "STOP: kissfft is not the pinned commit"; exit 2; }
SAN="-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer"
echo "== scratch build of 7b08827 ($SAN) -> $OUT/00ALLINONE_asan; sources = 7b08827: $CHECK; kissfft febd4ca"
gcc $LDFLAGS -o "$OUT/00ALLINONE_asan" $SRCS $SAN $(pkg-config --cflags sdl2 SDL2_ttf) \
    -DBUILD_GIT='"7b08827"' -DBUILD_TARGET='"asan-scratch"' -DBUILD_CFLAGS="\"$SAN\"" \
    $(pkg-config --libs sdl2 SDL2_ttf) -lGLEW -lGL -lm > "$OUT/build.log" 2>&1 \
  || { tail -20 "$OUT/build.log"; echo "STOP: scratch build failed"; exit 4; }
missing=$(ldd "$OUT/00ALLINONE_asan" | grep -c "not found") || true
[ "${missing:-0}" -eq 0 ] || { ldd "$OUT/00ALLINONE_asan" | grep "not found"; echo "STOP: unresolved shared libraries"; exit 5; }
ldd "$OUT/00ALLINONE_asan" | grep -E "asan|ubsan"
echo "sanitizer build: $("$OUT/00ALLINONE_asan" --version | head -1)"
echo "recorded build:  $(./00ALLINONE --version | head -1)"
echo "node $(hostname) | $(lscpu | awk -F: '/Model name/{gsub(/^ +/,"",$2); print $2; exit}') | out $OUT"
python3 cluster/resched_gate_261005/asan_runs_261007.py run --bin "$OUT/00ALLINONE_asan" --out "$OUT/runs" --jobs 4
python3 cluster/resched_gate_261005/asan_runs_261007.py report --out "$OUT/runs" > "$OUT/report.txt"
[ "$(sha256sum 00ALLINONE logs/BUILD_KOA_LAST.hash)" = "$REC_BEFORE" ] && echo "recorded binary and hash file unchanged" \
  || echo "WARNING: ./00ALLINONE or logs/BUILD_KOA_LAST.hash changed during this job -- stop and paste"
cat "$OUT/report.txt"
