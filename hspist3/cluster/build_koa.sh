#!/bin/bash
# ##CHRIS 2026-09-16: build 00ALLINONE on KOA (UH Manoa) and record exactly what was built.
#
# The simulation never opens a display: initSDL() and TTF_Init() sit behind `if (!cli_headless)` and
# every render function returns early when headless. SDL2/SDL2_ttf/GLEW/GL are therefore needed only
# to LINK.
#
# ##CHRIS 2026-10-02 (Task K2), REWRITTEN for what KOA turned out to be (KOA facts, 2026-10-03 UTC):
#   - /home is mounted noexec on the login node, so this script refuses to run outside a Slurm job (srun or sbatch);
#     the built binary runs fine on compute nodes.
#   - KOA has no SDL2_ttf module; the libraries come from the user conda env ~/envs/hd, set up by cluster/koa_env.sh
#     together with the compiler (module compiler/GCC/14.3.0). The old module-guessing fallback is gone.
#   - Provenance gate: the binary's build_git must equal `git rev-parse --short HEAD` of this checkout, with no -dirty.
#     A -dirty or unknown hash stops the build (a clean clone of the pushed commit is the intended input).
#
# usage (from hspist3/, inside `srun -p sandbox ... --pty /bin/bash` or from koa_smoketest.sh):  bash cluster/build_koa.sh
set -uo pipefail
cd "$(dirname "$0")/.."
[ -n "${SLURM_JOB_ID:-}" ] || { echo "STOP: run inside a Slurm job (srun/sbatch) -- home is noexec on the login node"; exit 2; }
source cluster/koa_env.sh || exit 2

echo "== toolchain";  gcc --version | head -1
echo "== libraries";  pkg-config --exists sdl2 SDL2_ttf glew || { echo "STOP: pkg-config does not find sdl2/SDL2_ttf/glew in ~/envs/hd"; exit 3; }
echo "   sdl2 $(pkg-config --modversion sdl2), SDL2_ttf $(pkg-config --modversion SDL2_ttf), glew $(pkg-config --modversion glew)"

echo "== build (portable ISA -march=x86-64-v2, not -march=native: every KOA node must produce the same bytes)"
make -B koa || { echo "STOP: make koa failed"; exit 4; }

echo "== checks"
missing=$(ldd ./00ALLINONE | grep -c "not found") || true
[ "${missing:-0}" -eq 0 ] || { ldd ./00ALLINONE | grep "not found"; echo "STOP: unresolved shared libraries"; exit 5; }
./00ALLINONE --version | head -2
head=$(git rev-parse --short HEAD)
ver=$(./00ALLINONE --version | head -1)
echo "$ver" | grep -q -- "git $head  target koa" || { echo "STOP: build_git is not the clean HEAD $head: $ver"; exit 6; }
./00ALLINONE --version | grep -q -- "-ffp-contract=off" || { echo "STOP: binary lacks -ffp-contract=off"; exit 6; }

echo "== provenance"
mkdir -p logs
{
  echo "built_utc      $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "host           $(hostname)   job ${SLURM_JOB_ID}"
  echo "cpu            $(lscpu | awk -F: '/Model name/{gsub(/^ +/,"",$2); print $2; exit}')"
  echo "compiler       $(gcc --version | head -1)"
  echo "version        $ver"
  echo "git_commit     $(git rev-parse HEAD)"
  echo "sha256         $(sha256sum 00ALLINONE | cut -d' ' -f1)"
  echo "libs           sdl2 $(pkg-config --modversion sdl2) SDL2_ttf $(pkg-config --modversion SDL2_ttf) glew $(pkg-config --modversion glew) (~/envs/hd)"
} | tee "logs/BUILD_KOA_${SLURM_JOB_ID}.txt"
echo "BUILD OK"
