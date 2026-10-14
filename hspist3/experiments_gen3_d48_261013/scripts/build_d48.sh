#!/bin/bash
# ##CHRIS 2026-10-13 (plan-author decision 48, item 3): the viewer with the cell-width edit, built with the command of the viewer
# Chris runs (experiments_gen3_d34s_261012/bin_406044f: -O3 -ffp-contract=off, target mac-O3-e0pre), into a NEW folder (nothing
# overwritten); then -Wall -Wextra on the driver's sources (syntax only), as the pinned build record did (25 warnings at 406044f).
set -euo pipefail
G=/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks_gen3/hspist3
cd "$G"
GIT=$(git -C "$G" rev-parse --short HEAD)$(git -C "$G" diff --quiet HEAD 2>/dev/null || echo -dirty)
OUT="$G/experiments_gen3_d48_261013/bin_$GIT"
[ -e "$OUT" ] && { echo "REFUSED: $OUT exists"; exit 1; }
mkdir -p "$OUT"
SRC=(00ALLINONE.c kissfft/kiss_fft.c kissfft/kiss_fftr.c edmd_core/edmd.c edmd_core/edmd_accelerated.c edmd_core/edmd_gen3.c experiment_validation.c)
INC=(-I/opt/homebrew/include/SDL2 -I/opt/homebrew/include -D_THREAD_SAFE -I/opt/homebrew/Cellar/harfbuzz/12.2.0_1/include/harfbuzz
     -I/opt/homebrew/Cellar/glib/2.86.3/include/glib-2.0 -I/opt/homebrew/Cellar/glib/2.86.3/lib/glib-2.0/include -I/opt/homebrew/opt/gettext/include
     -I/opt/homebrew/Cellar/pcre2/10.47_1/include -I/opt/homebrew/Cellar/graphite2/1.3.14/include -I/opt/homebrew/opt/freetype/include/freetype2
     -I/opt/homebrew/opt/libpng/include/libpng16)
DEF=(-DBUILD_GIT="\"$GIT\"" -DBUILD_TARGET='"mac-O3-e0pre"' -DBUILD_CFLAGS='"-O3 -ffp-contract=off"')
LIBS=(-L/opt/homebrew/Cellar/sdl2_ttf/2.24.0/lib -lSDL2_ttf -L/opt/homebrew/lib -lSDL2 -L/opt/homebrew/lib -lGLEW -framework OpenGL -framework CoreGraphics)
CMD=(cc -o "$OUT/00ALLINONE" "${SRC[@]}" -O3 -ffp-contract=off "${INC[@]}" "${DEF[@]}" "${LIBS[@]}")
t0=$(date +%s)
"${CMD[@]}"
t1=$(date +%s)
cc -fsyntax-only -Wall -Wextra -O3 -ffp-contract=off "${INC[@]}" "${DEF[@]}" "${SRC[@]}" 2> "$OUT/00ALLINONE.Wall.txt" || true
{
  echo "# the viewer for decision 48 (261012 sec. 4.7.104): $GIT ($(git -C "$G" rev-parse HEAD) plus the uncommitted edit of 00ALLINONE.c), $(date '+%Y-%m-%d %H:%M:%S %Z'), $(sysctl -n machdep.cpu.brand_string), $(cc --version | head -1)"
  echo
  echo "- built from the engine-gen3 worktree (not a clean export: the edit is not committed; the other dirty files of the worktree are not sources of the build)"
  echo "- the command (that of experiments_gen3_d34s_261012/bin_406044f/build_record.txt, with BUILD_GIT and the output changed), ${t1}s-${t0}s = $((t1 - t0)) s:"
  printf '  `'; printf '%q ' "${CMD[@]}"; printf '`\n'
  echo "- -Wall -Wextra on the driver's sources (syntax only): $(grep -c 'warning:' "$OUT/00ALLINONE.Wall.txt") warnings, by file: $(grep 'warning:' "$OUT/00ALLINONE.Wall.txt" | cut -d: -f1 | sort | uniq -c | tr '\n' ';')"
  echo
  echo "## SHA-256 of the sources and of the binary"
  echo
  shasum -a 256 "${SRC[@]}" edmd_core/edmd_gen3.h edmd_core/edmd.h kissfft/_kiss_fft_guts.h kissfft/kiss_fft.h kissfft/kiss_fft_log.h | sed 's/^/    /'
  (cd "$OUT" && shasum -a 256 00ALLINONE | sed 's/^/    /')
  echo
  echo "## the edit (git diff HEAD -- hspist3/00ALLINONE.c)"
  echo
  git -C "$G" diff HEAD -- 00ALLINONE.c | sed 's/^/    /'
} > "$OUT/build_record.txt"
echo "built $OUT/00ALLINONE in $((t1 - t0)) s; warnings: $(grep -c 'warning:' "$OUT/00ALLINONE.Wall.txt")"
