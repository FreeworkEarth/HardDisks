#!/bin/bash
# ##CHRIS 2026-10-09 (261012 sec. 4.7.14; programme rule 5): binaries from a COMMITTED tree only. The sources are `git archive`
# of the given commit (no working-tree file can enter, so the build line cannot be -dirty); kissfft is a submodule whose content
# git archive does not carry: its six source files are copied from the worktree's checkout (febd4cae plus the project's local
# 3-line change to kiss_fft_log.h, the same in every tree), and their state and SHA-256 are recorded.
# usage: build_clean.sh <commit> <out dir>      writes <out>/00ALLINONE, the test binaries, build_record.txt
set -eu
REV=$1; OUT=$2
WT=/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks_gen3
H=$(git -C "$WT" rev-parse --short=7 "$REV")
mkdir -p "$OUT/src"
git -C "$WT" archive "$REV" hspist3 | tar -x -C "$OUT/src"
K="$OUT/src/hspist3/kissfft"; mkdir -p "$K"
for f in kiss_fft.c kiss_fftr.c kiss_fft.h kiss_fftr.h _kiss_fft_guts.h kiss_fft_log.h; do cp "$WT/hspist3/kissfft/$f" "$K/"; done
R="$OUT/build_record.txt"
{
  echo "# clean build of engine-gen3 $(git -C "$WT" rev-parse "$REV") ($H), $(date '+%Y-%m-%d %H:%M:%S %Z'), $(sysctl -n machdep.cpu.brand_string), $(cc --version | head -1)"
  echo; echo "## kissfft (submodule; git archive does not carry it)"; echo
  echo "submodule commit $(git -C "$WT/hspist3/kissfft" rev-parse HEAD); local change (the same in main's tree):"
  git -C "$WT/hspist3/kissfft" diff | sed 's/^/    /'
  (cd "$K" && shasum -a 256 *) | sed 's/^/    /'
} > "$R"
cd "$OUT/src/hspist3"
SDLC=$(pkg-config --cflags sdl2 SDL2_ttf); SDLL=$(pkg-config --libs sdl2 SDL2_ttf)
DRV=(00ALLINONE.c kissfft/kiss_fft.c kissfft/kiss_fftr.c edmd_core/edmd.c edmd_core/edmd_accelerated.c edmd_core/edmd_gen3.c experiment_validation.c)
{ echo; echo "## builds (from hspist3/ of the archive)"; echo; } >> "$R"
build() {   # name, then the compiler arguments
  n=$1; shift
  cc "$@" 2> "$OUT/$n.warnings.txt"
  echo "- $n: \`cc $*\`" >> "$R"
  echo "  warnings: $(grep -c 'warning:' "$OUT/$n.warnings.txt" || true)" >> "$R"
}
build 00ALLINONE -o "$OUT/00ALLINONE" "${DRV[@]}" -O3 -ffp-contract=off $SDLC "-DBUILD_GIT=\"$H\"" '-DBUILD_TARGET="mac-O3-e0pre"' '-DBUILD_CFLAGS="-O3 -ffp-contract=off"' $SDLL -L/opt/homebrew/lib -lGLEW -framework OpenGL -framework CoreGraphics
cc -fsyntax-only -Wall -Wextra "${DRV[@]}" -O3 -ffp-contract=off $SDLC "-DBUILD_GIT=\"$H\"" '-DBUILD_TARGET="mac-O3-e0pre"' '-DBUILD_CFLAGS="-O3 -ffp-contract=off"' 2> "$OUT/00ALLINONE.Wall.txt" || true
echo "  -Wall -Wextra on the driver's sources (syntax only): $(grep -c 'warning:' "$OUT/00ALLINONE.Wall.txt" || true) warnings, by file: $(grep -o '^[a-z_A-Z/0-9.]*\.c' "$OUT/00ALLINONE.Wall.txt" | sort | uniq -c | tr -s ' ' | tr '\n' ';')" >> "$R"
build gen3_m1 -O3 -ffp-contract=off -Wall -Wextra -o "$OUT/gen3_m1" edmd_core/tests/gen3_m1_harness.c edmd_core/edmd_gen3.c edmd_core/edmd.c -lm
build gen3_m2 -std=c11 -O3 -ffp-contract=off -Wall -Wextra -o "$OUT/gen3_m2" edmd_core/tests/gen3_m2_harness.c edmd_core/edmd_gen3.c edmd_core/edmd.c -lm
build gen3_body_rule_test -std=c11 -O2 -ffp-contract=off -Wall -Wextra -o "$OUT/gen3_body_rule_test" edmd_core/tests/gen3_body_rule_test.c -lm
build gen3_band_edge_test -std=c11 -O2 -ffp-contract=off -Wall -Wextra -o "$OUT/gen3_band_edge_test" edmd_core/tests/gen3_band_edge_test.c edmd_core/edmd_gen3.c -lm
build gen3_band_edge_ties -std=c11 -O2 -ffp-contract=off -Wall -Wextra -o "$OUT/gen3_band_edge_ties" edmd_core/tests/gen3_band_edge_ties.c -lm
build edmd_gen3_alone -std=c11 -O3 -ffp-contract=off -Wall -Wextra -c edmd_core/edmd_gen3.c -o "$OUT/edmd_gen3.o"
{ echo; echo "## SHA-256"; echo; (cd "$OUT" && shasum -a 256 00ALLINONE gen3_m1 gen3_m2 gen3_body_rule_test gen3_band_edge_test gen3_band_edge_ties) | sed 's/^/    /'; } >> "$R"
echo "built $H"
