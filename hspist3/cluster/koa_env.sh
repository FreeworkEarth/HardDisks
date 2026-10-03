# ##CHRIS 2026-10-02: the ONE place that sets up the KOA build and run environment (Task K1). It is SOURCED (not run) by
# build_koa.sh, koa_smoketest.sh, koa_crossnode_det.sh and the confinement sbatch files, always INSIDE a Slurm job:
# /home is mounted noexec on the login node, so nothing built in home runs there (KOA facts, 2026-10-03 UTC).
#
# Compiler: the compiler/GCC/14.3.0 module, not the system gcc 11.5 -- a pinned version that does not change when the
# OS is patched (it is the GCC of KOA's 2025b toolchain), so a rebuild months from now uses the same compiler.
# `--version` of the binary does not record the compiler, so every job log prints `gcc --version` (see below).
#
# Libraries: SDL2, SDL2_ttf, GLEW and GL come from the user conda env ~/envs/hd (conda-forge only; exact package list in
# ~/envs/hd_explicit_261003.txt). They are needed to LINK only; headless runs never open a window. The binary is linked
# with an rpath to ~/envs/hd/lib, so it also runs from scripts that do not source this file.
if ! type module >/dev/null 2>&1; then
  for f in /etc/profile.d/lmod.sh /etc/profile.d/z00_lmod.sh /etc/profile.d/modules.sh; do [ -f "$f" ] && . "$f"; done
fi
[ -n "${SLURM_JOB_ID:-}" ] || echo "WARNING: koa_env.sh sourced outside a Slurm job -- build and run only inside srun/sbatch"
module purge
module load compiler/GCC/14.3.0 || { echo "STOP: module load compiler/GCC/14.3.0 failed"; return 1; }
export PKG_CONFIG_PATH=$HOME/envs/hd/lib/pkgconfig:${PKG_CONFIG_PATH:-}
export LD_LIBRARY_PATH=$HOME/envs/hd/lib:${LD_LIBRARY_PATH:-}
export PATH=$HOME/envs/hd/bin:$PATH
export LDFLAGS="-Wl,-rpath,$HOME/envs/hd/lib ${LDFLAGS:-}"     # Makefile:41 uses LDFLAGS +=, so this survives
gcc --version | head -1 | grep -q "14.3" || { echo "STOP: gcc is not 14.3: $(gcc --version | head -1)"; return 1; }
echo "env: $(gcc --version | head -1) | $(python3 --version 2>&1) at $(command -v python3) | $(git --version)"
