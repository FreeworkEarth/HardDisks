#!/bin/bash
# ##CHRIS 2026-09-16: build 00ALLINONE on KOA (UH Manoa) and record exactly what was built.
#
# The simulation never opens a display: initSDL() and TTF_Init() sit behind `if (!cli_headless)` and
# every render function returns early when headless. SDL2/SDL2_ttf/GLEW/GL are therefore needed only
# to LINK. First choice is cluster modules; the fallback is a user-space conda-forge environment,
# which needs no admin rights.
#
# usage: bash cluster/build_koa.sh        (run from the hspist3 directory on a KOA login node)
set -euo pipefail
cd "$(dirname "$0")/.."

echo "== toolchain"
command -v gcc >/dev/null || { echo "no gcc found; module load a compiler first"; exit 2; }
gcc --version | head -1

echo "== SDL2 / GLEW"
if pkg-config --exists sdl2 SDL2_ttf glew 2>/dev/null; then
  echo "pkg-config finds sdl2, SDL2_ttf and glew -- using them"
else
  echo "not found. Trying modules..."
  for m in sdl2 SDL2 glew GLEW mesa Mesa; do module load "$m" 2>/dev/null && echo "  loaded $m" || true; done
  if ! pkg-config --exists sdl2 SDL2_ttf glew 2>/dev/null; then
    cat <<'EOF'
  Still not found. User-space fallback (no admin needed):

      module load lang/Anaconda3 2>/dev/null || true
      conda create -y -p $HOME/envs/hd -c conda-forge sdl2 sdl2_ttf glew mesa-libgl-devel-cos7-x86_64 pkg-config gcc_linux-64
      conda activate $HOME/envs/hd
      export PKG_CONFIG_PATH=$HOME/envs/hd/lib/pkgconfig:$PKG_CONFIG_PATH
      export LD_LIBRARY_PATH=$HOME/envs/hd/lib:$LD_LIBRARY_PATH

  then re-run this script.
EOF
    exit 3
  fi
fi

echo "== build (portable ISA, not -march=native: every KOA node must produce the same bytes)"
make -B koa

echo "== provenance"
{
  echo "built_utc      $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "host           $(hostname)"
  echo "compiler       $(gcc --version | head -1)"
  echo "cflags         $(make -n koa | grep -o -- '-O2 .*' | head -1)"
  echo "git_commit     $(git -C .. rev-parse HEAD 2>/dev/null || echo 'not a git checkout')"
  echo "git_dirty      $(git -C .. status --porcelain 2>/dev/null | wc -l | tr -d ' ') modified paths"
  echo "sha256         $(sha256sum 00ALLINONE | cut -d' ' -f1)"
  echo "cpu            $(lscpu | awk -F: '/Model name/{gsub(/^ +/,"",$2); print $2; exit}')"
} | tee cluster/BUILD_KOA.txt

cat <<'EOF'

Next: the acceptance gate. The cluster binary is x86 and the laptop binary is ARM, so traces are NOT
byte-identical and never will be -- hard-disk dynamics is chaotic and the last bit diverges. Run the
mirror campaign and compare statistically before trusting any cluster number:

    python3 cluster/make_manifest.py mirror manifests/mirror.tsv
    mkdir -p logs && sbatch --array=1-150%150 cluster/run_array.sbatch manifests/mirror.tsv
    # when it finishes, rsync back and run, on the laptop:
    python3 validation/cluster_gate_20260916.py   # (written once the mirror data exists)
EOF
