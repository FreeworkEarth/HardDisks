#!/bin/bash
# ##CHRIS 2026-10-09 (plan-author decision of 2026-10-09, item 3): the loader provenance guard on engine-divider-resched, the
# decision-2 standard. The paper scripts run in main's working tree (the only tree with the data), once with the branch's code
# (its 14 pre-guard loader files and its reduce_B.py) and once with the branch plus the guard (main's 14 guarded files, the same
# reduce_B.py); then compare_before_after.sh. Added: the T-prime verdict script (it imports tests_20260913) in both states.
# Every swapped file is checked by blob hash. Only the 45 known outputs of decision 2 are restored; anything else is reported.
set -u
SP=/private/tmp/claude-501/-Users-chrisharing-Desktop-CCS-complex-coupled-systems-Repo-HardDisks/91cb08ec-0599-4faa-a6af-d5a7ca834255/scratchpad
REPO=/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks
OUT=$SP/gguard_edr; mkdir -p "$OUT"
cd "$REPO" || exit 1
git status --porcelain --untracked-files=all > "$OUT/status_before.txt"
git diff --name-only --diff-filter=M -z | tar --null -n -T - -cf "$OUT/backup_modified_tracked.tar"    # Chris's modified files, untouched; a copy
L="hspist3/validation/estimator_massladder_20260917.py hspist3/validation/level2_Au_figure_20260918.py
hspist3/validation/paper1_A2_boxtrunc_261002.py hspist3/validation/paper1_boxtrunc_20261014.py
hspist3/validation/paper1_confinement_afix_261005.py hspist3/validation/paper1_confinement_heldwall_posthoc_261004.py
hspist3/validation/paper1_confinement_prereg_20261012.py hspist3/validation/paper1_confinement_results_261004.py
hspist3/validation/paper1_populate_cs_err_20261002.py hspist3/validation/paper2_figures_20261001.py
hspist3/validation/paper2_geometry_fix_20260918.py hspist3/validation/paper2_ramp_fast_20260918.py
hspist3/validation/resched_gate_261005.py hspist3/validation/tests_20260913.py"
RB=hspist3/cluster/confinement_20261013/reduce_B.py
EXPECTED=hspist3/experiments_loader_guard_261008/written_before.sha256
for f in $L $RB; do [ -z "$(git status --porcelain -- "$f")" ] || { echo "STOP: $f is not clean"; exit 1; }; done
for f in $L; do
  [ "$(git rev-parse "HEAD:$f")" = "$(git rev-parse "0b21270:$f")" ] || { echo "STOP: $f changed on main after the guard commit"; exit 1; }
  [ "$(git rev-parse "engine-divider-resched:$f")" = "$(git rev-parse "0b21270~1:$f")" ] || { echo "STOP: the branch's $f is not main's pre-guard version"; exit 1; }
done
awk '{print $2}' "$EXPECTED" | grep -v '^hspist3/validation/' | while read -r p; do
  if git ls-files --error-unmatch "$p" > /dev/null 2>&1; then [ -z "$(git status --porcelain -- "$p")" ] || echo "STOP"
  else [ ! -e "$p" ] || echo "STOP"; fi
done | grep -q STOP && { echo "STOP: an expected output is modified or present"; exit 1; }
echo "checks: the 15 swapped files clean; main's 14 = 0b21270; the branch's 14 = 0b21270~1; the 41 tracked outputs clean, the 4 untracked absent" > "$OUT/log.txt"
trap 'git checkout -- $L $RB' EXIT     # whatever happens, the 15 files go back to main's (they were clean)
tprime() {   # the T-prime verdict script, recorded like run_set.sh's scripts
  D=$SP/gguard/$1; touch "$D/.marker_tprime"; sleep 1; t0=$(date +%s)
  (cd hspist3 && MPLBACKEND=Agg SOURCE_DATE_EPOCH=0 PYTHONHASHSEED=0 python3 validation/resched_testTprime_261007.py > "$D/tprime.out" 2> "$D/tprime.err"; echo $? > "$D/tprime.rc")
  echo "tprime $(( $(date +%s) - t0 )) s rc $(cat "$D/tprime.rc")" >> "$D/times.txt"
  find "$REPO" -newer "$D/.marker_tprime" -type f -not -path '*/.git/*' -not -path '*/__pycache__/*' | sort > "$D/tprime_written.txt"
}
# ---- BEFORE: the branch's code
for f in $L $RB; do git show "engine-divider-resched:$f" > "$f"; done
for f in $L $RB; do [ "$(git hash-object "$f")" = "$(git rev-parse "engine-divider-resched:$f")" ] || { echo "STOP: swap of $f"; exit 1; }; done
echo "before: the branch's 14 loader files and reduce_B.py in place (blob hashes checked) $(date '+%H:%M:%S')" >> "$OUT/log.txt"
bash hspist3/experiments_loader_guard_261008/run_set.sh edr_before
tprime edr_before
# ---- AFTER: the branch plus the guard
git checkout -- $L
for f in $L; do [ "$(git hash-object "$f")" = "$(git rev-parse "HEAD:$f")" ] || { echo "STOP: restore of $f"; exit 1; }; done
[ "$(git hash-object "$RB")" = "$(git rev-parse "engine-divider-resched:$RB")" ] || { echo "STOP: reduce_B.py is not the branch's"; exit 1; }
echo "after: main's 14 guarded files (= the branch plus the guard), the branch's reduce_B.py (blob hashes checked) $(date '+%H:%M:%S')" >> "$OUT/log.txt"
bash hspist3/experiments_loader_guard_261008/run_set.sh edr_after
tprime edr_after
git checkout -- $RB
[ "$(git hash-object "$RB")" = "$(git rev-parse "HEAD:$RB")" ] || { echo "STOP: restore of reduce_B.py"; exit 1; }
echo "restored: reduce_B.py = main's $(date '+%H:%M:%S')" >> "$OUT/log.txt"
# ---- compare and unit tests
bash hspist3/experiments_loader_guard_261008/compare_before_after.sh "$SP/gguard/edr_before" "$SP/gguard/edr_after" > "$OUT/compare_before_after_output.txt" 2>&1
for s in edr_before edr_after; do echo "T-prime run $s wrote $(wc -l < "$SP/gguard/$s/tprime_written.txt" | tr -d ' ') files" >> "$OUT/log.txt"; done
(cd hspist3 && python3 -m unittest validation/test_edmd_acc_guard.py -v) > "$OUT/unit_tests_output.txt" 2>&1
echo "unit tests exit $?" >> "$OUT/log.txt"
# ---- restore the expected outputs: tracked (clean before) -> git checkout; untracked (absent before) -> moved to scratch
mkdir -p "$SP/gguard/edr_moved"
awk '{print $2}' "$EXPECTED" | grep -v '^hspist3/validation/' > "$OUT/expected_outputs.txt"
cat "$SP/gguard/edr_before/written.txt" "$SP/gguard/edr_after/written.txt" "$SP/gguard/edr_before/tprime_written.txt" "$SP/gguard/edr_after/tprime_written.txt" \
  | sed "s#^$REPO/##" | grep -v '^hspist3/validation/' | sort -u > "$OUT/written_union.txt"
while read -r rel; do
  if ! grep -qxF "$rel" "$OUT/expected_outputs.txt"; then echo "ATTENTION (written, not an expected output; left as is): $rel" >> "$OUT/log.txt"
  elif git ls-files --error-unmatch "$rel" > /dev/null 2>&1; then git checkout -- "$rel" && echo "restored from git: $rel" >> "$OUT/log.txt"
  elif [ -e "$rel" ]; then mkdir -p "$SP/gguard/edr_moved/$(dirname "$rel")"; mv "$rel" "$SP/gguard/edr_moved/$rel" && echo "moved (new, untracked): $rel" >> "$OUT/log.txt"; fi
done < "$OUT/written_union.txt"
git status --porcelain --untracked-files=all > "$OUT/status_after.txt"
if cmp -s "$OUT/status_before.txt" "$OUT/status_after.txt"; then echo "working tree: status identical to before the runs" >> "$OUT/log.txt"
else echo "working tree: STATUS DIFFERS from before" >> "$OUT/log.txt"; diff "$OUT/status_before.txt" "$OUT/status_after.txt" >> "$OUT/log.txt"; fi
echo "DONE $(date '+%H:%M:%S')" >> "$OUT/log.txt"
