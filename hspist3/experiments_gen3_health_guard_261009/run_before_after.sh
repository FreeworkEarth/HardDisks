#!/bin/bash
# ##CHRIS 2026-10-09 (261012 sec. 4.7.14, M3 amendment e): the gen3 run-record guard (validation/edmd3_health_guard.py, called
# from edmd_acc_guard.guard()), the decision-2 standard: the paper scripts run in main's working tree once with main's HEAD code
# (before) and once with the guard installed (after); then compare_before_after.sh, the T-prime verdict script in both states,
# and the unit tests. The guard files are staged in the scratch folder G and installed between the two sets (checked by
# content: edmd_acc_guard.py must differ from HEAD by the hook lines only). Only the known outputs of decision 2 are restored;
# anything else is reported.
set -u
SP=/private/tmp/claude-501/-Users-chrisharing-Desktop-CCS-complex-coupled-systems-Repo-HardDisks/91cb08ec-0599-4faa-a6af-d5a7ca834255/scratchpad
G=$SP/g3guard/stage
REPO=/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks
OUT=$SP/g3guard/out; mkdir -p "$OUT"
cd "$REPO" || exit 1
V=hspist3/validation
EXPECTED=hspist3/experiments_loader_guard_261008/written_before.sha256
git status --porcelain --untracked-files=all > "$OUT/status_before.txt"
git diff --name-only --diff-filter=M -z | tar --null -n -T - -cf "$OUT/backup_modified_tracked.tar"    # Chris's modified files, untouched; a copy
[ -z "$(git status --porcelain -- $V/edmd_acc_guard.py)" ] || { echo "STOP: edmd_acc_guard.py is not clean"; exit 1; }
[ ! -e $V/edmd3_health_guard.py ] && [ ! -e $V/test_edmd3_health_guard.py ] || { echo "STOP: a guard file exists already"; exit 1; }
diff $V/edmd_acc_guard.py "$G/edmd_acc_guard.py" > "$OUT/hook_diff.txt"
[ "$(grep -c '^[<>]' "$OUT/hook_diff.txt")" = 3 ] && [ "$(grep -c '^<' "$OUT/hook_diff.txt")" = 0 ] || { echo "STOP: the staged edmd_acc_guard.py is not HEAD plus the 3 hook lines"; exit 1; }
awk '{print $2}' "$EXPECTED" | grep -v '^hspist3/validation/' | while read -r p; do
  if git ls-files --error-unmatch "$p" > /dev/null 2>&1; then [ -z "$(git status --porcelain -- "$p")" ] || echo "STOP"
  else [ ! -e "$p" ] || echo "STOP"; fi
done | grep -q STOP && { echo "STOP: an expected output is modified or present"; exit 1; }
echo "checks: edmd_acc_guard.py clean (HEAD $(git rev-parse --short HEAD)); the staged copy = HEAD + 3 hook lines; the guard files absent; the 41 tracked outputs clean, the 4 untracked absent" > "$OUT/log.txt"
tprime() {
  D=$SP/gguard/$1; touch "$D/.marker_tprime"; sleep 1; t0=$(date +%s)
  (cd hspist3 && MPLBACKEND=Agg SOURCE_DATE_EPOCH=0 PYTHONHASHSEED=0 python3 validation/resched_testTprime_261007.py > "$D/tprime.out" 2> "$D/tprime.err"; echo $? > "$D/tprime.rc")
  echo "tprime $(( $(date +%s) - t0 )) s rc $(cat "$D/tprime.rc")" >> "$D/times.txt"
  find "$REPO" -newer "$D/.marker_tprime" -type f -not -path '*/.git/*' -not -path '*/__pycache__/*' | sort > "$D/tprime_written.txt"
}
# ---- BEFORE: main's HEAD code
echo "before: HEAD code $(date '+%H:%M:%S')" >> "$OUT/log.txt"
bash hspist3/experiments_loader_guard_261008/run_set.sh g3guard_before
tprime g3guard_before
# ---- install the guard
cp "$G/edmd3_health_guard.py" "$G/test_edmd3_health_guard.py" $V/ && cp "$G/edmd_acc_guard.py" $V/edmd_acc_guard.py
cmp -s "$G/edmd_acc_guard.py" $V/edmd_acc_guard.py && cmp -s "$G/edmd3_health_guard.py" $V/edmd3_health_guard.py || { echo "STOP: install"; exit 1; }
echo "installed: edmd3_health_guard.py, test_edmd3_health_guard.py, edmd_acc_guard.py (hook) $(date '+%H:%M:%S')" >> "$OUT/log.txt"
# ---- AFTER: HEAD plus the guard
bash hspist3/experiments_loader_guard_261008/run_set.sh g3guard_after
tprime g3guard_after
# ---- compare and unit tests
bash hspist3/experiments_loader_guard_261008/compare_before_after.sh "$SP/gguard/g3guard_before" "$SP/gguard/g3guard_after" > "$OUT/compare_before_after_output.txt" 2>&1
for s in g3guard_before g3guard_after; do echo "T-prime run $s wrote $(wc -l < "$SP/gguard/$s/tprime_written.txt" | tr -d ' ') files" >> "$OUT/log.txt"; done
(cd hspist3 && python3 -m unittest validation/test_edmd3_health_guard.py validation/test_edmd_acc_guard.py -v) > "$OUT/unit_tests_output.txt" 2>&1
echo "unit tests exit $?" >> "$OUT/log.txt"
# ---- restore the expected outputs: tracked (clean before) -> git checkout; untracked (absent before) -> moved to scratch
mkdir -p "$SP/gguard/g3guard_moved"
awk '{print $2}' "$EXPECTED" | grep -v '^hspist3/validation/' > "$OUT/expected_outputs.txt"
cat "$SP/gguard/g3guard_before/written.txt" "$SP/gguard/g3guard_after/written.txt" "$SP/gguard/g3guard_before/tprime_written.txt" "$SP/gguard/g3guard_after/tprime_written.txt" \
  | sed "s#^$REPO/##" | grep -v '^hspist3/validation/' | sort -u > "$OUT/written_union.txt"
while read -r rel; do
  if ! grep -qxF "$rel" "$OUT/expected_outputs.txt"; then echo "ATTENTION (written, not an expected output; left as is): $rel" >> "$OUT/log.txt"
  elif git ls-files --error-unmatch "$rel" > /dev/null 2>&1; then git checkout -- "$rel" && echo "restored from git: $rel" >> "$OUT/log.txt"
  elif [ -e "$rel" ]; then mkdir -p "$SP/gguard/g3guard_moved/$(dirname "$rel")"; mv "$rel" "$SP/gguard/g3guard_moved/$rel" && echo "moved (new, untracked): $rel" >> "$OUT/log.txt"; fi
done < "$OUT/written_union.txt"
git status --porcelain --untracked-files=all > "$OUT/status_after.txt"
diff "$OUT/status_before.txt" "$OUT/status_after.txt" > "$OUT/status_diff.txt"
echo "working tree: status difference before -> after (expected: edmd_acc_guard.py modified, the two guard files new):" >> "$OUT/log.txt"
cat "$OUT/status_diff.txt" >> "$OUT/log.txt"
echo "DONE $(date '+%H:%M:%S')" >> "$OUT/log.txt"
