#!/bin/bash
# ##CHRIS 2026-10-11 (plan-author decision 28, C1; 261012 sec. 4.7.57): the Mac binary of engine-gen3 a643bc7 for C1, built from a
# COMMITTED tree only, as experiments_gen3_m3_261009/build_clean.sh builds its 00ALLINONE (the same compile line, target mac-O3-e0pre),
# with one change: `git archive` takes only the driver's sources (00ALLINONE.c, edmd_core/, experiment_validation.[ch]) instead of all
# of hspist3 (about 115 MB), so the source copy stays small. kissfft is copied from the worktree's checkout as build_clean.sh does
# (submodule febd4cae plus the project's local 3-line change to kiss_fft_log.h), and recorded.
# usage: build_mac_a643bc7.sh <out dir, absolute>      writes <out>/00ALLINONE, <out>/build_record.txt
set -eu
OUT=$1; case "$OUT" in /*) ;; *) echo "OUT must be an absolute path"; exit 2 ;; esac
REV=a643bc7edb1b41574c0f9af0c5fea199dc797680
WT=/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks_gen3
H=$(git -C "$WT" rev-parse --short=7 "$REV")
[ ! -e "$OUT/00ALLINONE" ] || { echo "STOP: $OUT/00ALLINONE exists"; exit 2; }
mkdir -p "$OUT/src"
git -C "$WT" archive "$REV" hspist3/00ALLINONE.c hspist3/edmd_core hspist3/experiment_validation.c hspist3/experiment_validation.h \
  | tar -x -C "$OUT/src"
K="$OUT/src/hspist3/kissfft"; mkdir -p "$K"
for f in kiss_fft.c kiss_fftr.c kiss_fft.h kiss_fftr.h _kiss_fft_guts.h kiss_fft_log.h; do cp "$WT/hspist3/kissfft/$f" "$K/"; done
R="$OUT/build_record.txt"
{
  echo "# clean build of engine-gen3 $(git -C "$WT" rev-parse "$REV") ($H), $(date '+%Y-%m-%d %H:%M:%S %Z'), $(sysctl -n machdep.cpu.brand_string), $(cc --version | head -1)"
  echo; echo "## sources: git archive $H of hspist3/00ALLINONE.c hspist3/edmd_core hspist3/experiment_validation.c hspist3/experiment_validation.h"
  (cd "$OUT/src/hspist3" && shasum -a 256 00ALLINONE.c edmd_core/*.c edmd_core/*.h experiment_validation.c experiment_validation.h) | sed 's/^/    /'
  echo; echo "## kissfft (submodule; git archive does not carry it)"; echo
  echo "submodule commit $(git -C "$WT/hspist3/kissfft" rev-parse HEAD); local change (the same in main's tree):"
  git -C "$WT/hspist3/kissfft" diff | sed 's/^/    /'
  (cd "$K" && shasum -a 256 *) | sed 's/^/    /'
} > "$R"
cd "$OUT/src/hspist3"
SDLC=$(pkg-config --cflags sdl2 SDL2_ttf); SDLL=$(pkg-config --libs sdl2 SDL2_ttf)
DRV=(00ALLINONE.c kissfft/kiss_fft.c kissfft/kiss_fftr.c edmd_core/edmd.c edmd_core/edmd_accelerated.c edmd_core/edmd_gen3.c experiment_validation.c)
ARGS=(-o "$OUT/00ALLINONE" "${DRV[@]}" -O3 -ffp-contract=off $SDLC "-DBUILD_GIT=\"$H\"" '-DBUILD_TARGET="mac-O3-e0pre"' '-DBUILD_CFLAGS="-O3 -ffp-contract=off"' $SDLL -L/opt/homebrew/lib -lGLEW -framework OpenGL -framework CoreGraphics)
cc "${ARGS[@]}" 2> "$OUT/00ALLINONE.warnings.txt"
{ echo; echo "## build (from hspist3/ of the archive)"; echo; echo "- 00ALLINONE: \`cc ${ARGS[*]}\`"
  echo "  warnings: $(grep -c 'warning:' "$OUT/00ALLINONE.warnings.txt" || true)"
  echo; echo "## SHA-256"; echo; (cd "$OUT" && shasum -a 256 00ALLINONE) | sed 's/^/    /'
  echo; echo "## --version"; echo; "$OUT/00ALLINONE" --version | sed 's/^/    /'; } >> "$R"
echo "built $H: $("$OUT/00ALLINONE" --version | head -1)"
