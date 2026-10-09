#!/bin/bash
# ##CHRIS 2026-10-09 (261012 sec. 4.7.15): the guard fix (the file token of per-run files is (int)(L0 * 10)) against the paper
# scripts: run_set.sh with the fixed guard, compared with the "after" set of run_before_after.sh (the first version of the guard);
# the outputs restored as there. Expected: identical (the fix touches only gen3 per-run files; the paper data are gen2).
set -u
SP=/private/tmp/claude-501/-Users-chrisharing-Desktop-CCS-complex-coupled-systems-Repo-HardDisks/91cb08ec-0599-4faa-a6af-d5a7ca834255/scratchpad
REPO=/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks
OUT=$SP/g3guard/fixcheck; mkdir -p "$OUT"; cd "$REPO" || exit 1
EXPECTED=hspist3/experiments_loader_guard_261008/written_before.sha256
git status --porcelain --untracked-files=all > "$OUT/status_before.txt"
bash hspist3/experiments_loader_guard_261008/run_set.sh g3guard_fix
D=$SP/gguard/g3guard_fix; touch "$D/.marker_tprime"; sleep 1
(cd hspist3 && MPLBACKEND=Agg SOURCE_DATE_EPOCH=0 PYTHONHASHSEED=0 python3 validation/resched_testTprime_261007.py > "$D/tprime.out" 2> "$D/tprime.err"; echo $? > "$D/tprime.rc")
find "$REPO" -newer "$D/.marker_tprime" -type f -not -path '*/.git/*' -not -path '*/__pycache__/*' | sort > "$D/tprime_written.txt"
bash hspist3/experiments_loader_guard_261008/compare_before_after.sh "$SP/gguard/g3guard_after" "$SP/gguard/g3guard_fix" > "$OUT/compare_output.txt" 2>&1
mkdir -p "$SP/gguard/g3guard_fix_moved"
awk '{print $2}' "$EXPECTED" | grep -v '^hspist3/validation/' > "$OUT/expected_outputs.txt"
cat "$D/written.txt" "$D/tprime_written.txt" | sed "s#^$REPO/##" | grep -v '^hspist3/validation/' | sort -u > "$OUT/written_union.txt"
while read -r rel; do
  if ! grep -qxF "$rel" "$OUT/expected_outputs.txt"; then echo "ATTENTION (written, not an expected output; left as is): $rel" >> "$OUT/log.txt"
  elif git ls-files --error-unmatch "$rel" > /dev/null 2>&1; then git checkout -- "$rel" && echo "restored from git: $rel" >> "$OUT/log.txt"
  elif [ -e "$rel" ]; then mkdir -p "$SP/gguard/g3guard_fix_moved/$(dirname "$rel")"; mv "$rel" "$SP/gguard/g3guard_fix_moved/$rel" && echo "moved (new, untracked): $rel" >> "$OUT/log.txt"; fi
done < "$OUT/written_union.txt"
git status --porcelain --untracked-files=all > "$OUT/status_after.txt"
cmp -s "$OUT/status_before.txt" "$OUT/status_after.txt" && echo "working tree: status identical to before the runs" >> "$OUT/log.txt" || { echo "working tree: STATUS DIFFERS" >> "$OUT/log.txt"; diff "$OUT/status_before.txt" "$OUT/status_after.txt" >> "$OUT/log.txt"; }
echo "DONE $(date '+%H:%M:%S')" >> "$OUT/log.txt"
