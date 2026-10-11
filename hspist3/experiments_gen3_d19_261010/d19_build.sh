#!/bin/bash
# ##CHRIS 2026-10-10 (plan-author decision 19, item 2; 261012 sec. 4.7.50): the binaries of a COMMITTED tree for the acceptance of gen3 in
# the window (a provisional DISPLAY feature). build_clean.sh builds into a scratch folder (its git-archive source copy stays there: the
# ground rule), and only the binaries, the build record and the warning lists are copied into the experiment's bin folder; the parent
# build's warning list (f72d764, decision 18) is compared name by name: "0 new warnings" is decision 19's acceptance item.
# usage: d19_build.sh <commit> <scratch build dir, absolute> <bin dir, absolute>
set -eu
REV=$1; OUT=$2; BIN=$3
case "$OUT" in /*) ;; *) echo "STOP: the scratch build dir must be absolute"; exit 2;; esac
case "$BIN" in /*) ;; *) echo "STOP: the bin dir must be absolute"; exit 2;; esac
[ -e "$OUT" ] && { echo "STOP: $OUT exists"; exit 2; }
[ -e "$BIN" ] && { echo "STOP: $BIN exists"; exit 2; }
HS=$(cd "$(dirname "$0")/.." && pwd)
bash "$HS/experiments_gen3_m3_261009/build_clean.sh" "$REV" "$OUT"
mkdir -p "$BIN"
for f in 00ALLINONE gen3_m1 gen3_m2 build_record.txt 00ALLINONE.Wall.txt 00ALLINONE.warnings.txt gen3_m1.warnings.txt gen3_m2.warnings.txt; do
  cp -p "$OUT/$f" "$BIN/"
done
# the parent f72d764 (decision 18's build): between it and the WIP b416dfc only 00ALLINONE.c changed among the C sources, so its
# warnings come from the same -Wall -Wextra syntax pass of build_clean.sh with the parent's 00ALLINONE.c in this archive
WT=/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks_gen3
[ -z "$(git -C "$WT" diff --name-only f72d764 "$REV" -- 'hspist3/*.c' 'hspist3/*.h' hspist3/edmd_core hspist3/kissfft | grep -v '^hspist3/00ALLINONE.c$')" ] \
  || { echo "STOP: C sources other than 00ALLINONE.c differ between f72d764 and $REV"; exit 2; }
git -C "$WT" show f72d764:hspist3/00ALLINONE.c > "$OUT/src/hspist3/00ALLINONE_parent_f72d764.c"
( cd "$OUT/src/hspist3" && cc -fsyntax-only -Wall -Wextra 00ALLINONE_parent_f72d764.c kissfft/kiss_fft.c kissfft/kiss_fftr.c edmd_core/edmd.c \
    edmd_core/edmd_accelerated.c edmd_core/edmd_gen3.c experiment_validation.c -O3 -ffp-contract=off $(pkg-config --cflags sdl2 SDL2_ttf) \
    '-DBUILD_GIT="f72d764"' '-DBUILD_TARGET="mac-O3-e0pre"' '-DBUILD_CFLAGS="-O3 -ffp-contract=off"' ) 2> "$BIN/00ALLINONE.Wall_parent_f72d764.txt" || true
echo "the parent's -Wall -Wextra pass: $(grep -c 'warning:' "$BIN/00ALLINONE.Wall_parent_f72d764.txt" || true) warnings; this build's: $(grep -c 'warning:' "$BIN/00ALLINONE.Wall.txt" || true)"
echo "built $REV; binaries, build record and warning lists copied to $BIN; the source copy stays in $OUT/src"
