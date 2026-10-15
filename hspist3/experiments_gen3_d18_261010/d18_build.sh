#!/bin/bash
# ##CHRIS 2026-10-10 (plan-author decision 18; 261012 sec. 4.7.42): the binaries of a COMMITTED tree for items 1-3. build_clean.sh
# builds into a scratch folder (its git-archive source copy stays there: decision 18's ground rule), then stage H's checkpoint test
# from the same archive (default and -DEDMD3_LONG_DOUBLE, as stage H's build record), and only the binaries and the build record are
# copied into the experiment's bin folder.
# usage: d18_build.sh <commit> <scratch build dir, absolute> <bin dir, absolute>
set -eu
REV=$1; OUT=$2; BIN=$3
case "$OUT" in /*) ;; *) echo "STOP: the scratch build dir must be absolute"; exit 2;; esac
case "$BIN" in /*) ;; *) echo "STOP: the bin dir must be absolute"; exit 2;; esac
[ -e "$OUT" ] && { echo "STOP: $OUT exists"; exit 2; }
[ -e "$BIN" ] && { echo "STOP: $BIN exists"; exit 2; }
HS=$(cd "$(dirname "$0")/.." && pwd)
bash "$HS/experiments_gen3_m3_261009/build_clean.sh" "$REV" "$OUT"
cd "$OUT/src/hspist3"
R="$OUT/build_record.txt"
{ echo; echo "## decision 18: stage H's checkpoint test from the same archive (default and -DEDMD3_LONG_DOUBLE)"; echo; } >> "$R"
for v in default ld; do
  if [ "$v" = ld ]; then FLAG=-DEDMD3_LONG_DOUBLE; SUF=_ld; else FLAG=; SUF=; fi
  cc -std=c11 -O3 -ffp-contract=off -Wall -Wextra $FLAG -o "$OUT/gen3_checkpoint_test$SUF" edmd_core/tests/gen3_checkpoint_test.c edmd_core/edmd_gen3.c edmd_core/edmd.c -lm 2> "$OUT/gen3_checkpoint_test$SUF.warnings.txt"
  echo "- gen3_checkpoint_test$SUF: \`cc -std=c11 -O3 -ffp-contract=off -Wall -Wextra $FLAG -o $OUT/gen3_checkpoint_test$SUF edmd_core/tests/gen3_checkpoint_test.c edmd_core/edmd_gen3.c edmd_core/edmd.c -lm\`" >> "$R"
  echo "  warnings: $(grep -c 'warning:' "$OUT/gen3_checkpoint_test$SUF.warnings.txt" || true)" >> "$R"
done
{ echo; (cd "$OUT" && shasum -a 256 gen3_checkpoint_test gen3_checkpoint_test_ld) | sed 's/^/    /'; } >> "$R"
mkdir -p "$BIN"
for f in 00ALLINONE gen3_m1 gen3_m2 gen3_body_rule_test gen3_band_edge_test gen3_band_edge_ties gen3_checkpoint_test gen3_checkpoint_test_ld build_record.txt; do
  cp -p "$OUT/$f" "$BIN/"
done
echo "built $REV; binaries and the build record copied to $BIN; the source copy stays in $OUT/src"
