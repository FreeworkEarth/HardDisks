#!/bin/bash
# ##CHRIS 2026-10-09 (261012 sec. 4.7.10, decision 5): which data files feed the papers' figures and tables. The 20 paper scripts
# of experiments_loader_guard_261008/run_set.sh run exactly as there (python3 validation/<script>.py), with PYTHONPATH set to
# guardhook/, whose usercustomize.py records every path the loaders pass to edmd_acc_guard.guard() (also in Pool workers and
# subprocesses). stdout/stderr are compared with the run of the same code this morning (gguard/edr_after) to show the hook is
# transparent. Only the 45 known outputs are restored afterwards (tracked: git checkout; new: moved to scratch);
# anything else written is reported, not touched. The working tree's status must be the same before and after.
set -u
SP=/private/tmp/claude-501/-Users-chrisharing-Desktop-CCS-complex-coupled-systems-Repo-HardDisks/91cb08ec-0599-4faa-a6af-d5a7ca834255/scratchpad
REPO=/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks
OUT=$SP/d5/trace; mkdir -p "$OUT"
EXPECTED=hspist3/experiments_loader_guard_261008/written_before.sha256
cd "$REPO" || exit 1
git status --porcelain --untracked-files=all > "$OUT/status_before.txt"
awk '{print $2}' "$EXPECTED" | grep -v '^hspist3/validation/' > "$OUT/expected_outputs.txt"
while read -r p; do
  if git ls-files --error-unmatch "$p" > /dev/null 2>&1; then [ -z "$(git status --porcelain -- "$p")" ] || { echo "STOP: $p modified"; exit 1; }
  else [ ! -e "$p" ] || { echo "STOP: $p present"; exit 1; }; fi
done < "$OUT/expected_outputs.txt"
cd hspist3 || exit 1
export MPLBACKEND=Agg SOURCE_DATE_EPOCH=0 PYTHONHASHSEED=0
touch "$OUT/.marker"; sleep 1
run() { name=$1; shift; t0=$(date +%s)
  GUARD_LOG="$OUT/paths_$name.txt" GUARD_VALIDATION_DIR="$PWD/validation" PYTHONPATH="$PWD/experiments_gen2_defect_scan_261009/guardhook" \
    python3 "$@" > "$OUT/$name.out" 2> "$OUT/$name.err"; rc=$?
  echo "$name $(( $(date +%s) - t0 )) s rc $rc paths $(sort -u "$OUT/paths_$name.txt" 2>/dev/null | wc -l | tr -d ' ')" >> "$OUT/times.txt"; }
: > "$OUT/times.txt"
run populate      validation/paper1_populate_cs_err_20261002.py
run figures       validation/paper1_figures_20261001.py
run damping       validation/damping_test_20260915.py
run massladder    validation/estimator_massladder_20260917.py
run a2boxtrunc    validation/paper1_A2_boxtrunc_261002.py
run boxtrunc      validation/paper1_boxtrunc_20261014.py
run boxtrunc_tab  validation/paper1_boxtrunc_20261014.py --table
run conf_results  validation/paper1_confinement_results_261004.py
run conf_afix     validation/paper1_confinement_afix_261005.py
run conf_heldwall validation/paper1_confinement_heldwall_posthoc_261004.py
run conf_prereg   validation/paper1_confinement_prereg_20261012.py
run resched_gate  validation/resched_gate_261005.py
run p2_figures    validation/paper2_figures_20261001.py
run p2_geomfix    validation/paper2_geometry_fix_20260918.py
run p2_rampfast   validation/paper2_ramp_fast_20260918.py
run p2_level2Au   validation/level2_Au_figure_20260918.py
run canonical     validation/paper1_canonical_20260919.py
run melting       validation/paper1_melting_figure_20261002.py
run draft_audit   validation/paper1_draft_audit_20261014.py
run roman         validation/roman2002_remapped_20260922.py
cd "$REPO" || exit 1
find "$REPO" -newer "$OUT/.marker" -type f -not -path '*/.git/*' -not -path '*/__pycache__/*' -not -path '*/.mplconfig/*' | sed "s#^$REPO/##" | sort > "$OUT/written.txt"
# transparency: the same stdout/stderr as this morning's run of the same code (edr_after)
for f in "$OUT"/*.out "$OUT"/*.err; do b=$(basename "$f"); cmp -s "$f" "$SP/gguard/edr_after/$b" && echo "same as this morning: $b" || echo "DIFFERS from this morning: $b"; done > "$OUT/transparency.txt"
# restore: only the expected outputs
mkdir -p "$SP/d5/moved"
while read -r rel; do
  case "$rel" in hspist3/validation/*) continue ;; esac
  if ! grep -qxF "$rel" "$OUT/expected_outputs.txt"; then echo "ATTENTION (written, not an expected output; left as is): $rel"
  elif git ls-files --error-unmatch "$rel" > /dev/null 2>&1; then git checkout -- "$rel" && echo "restored from git: $rel"
  elif [ -e "$rel" ]; then mkdir -p "$SP/d5/moved/$(dirname "$rel")"; mv "$rel" "$SP/d5/moved/$rel" && echo "moved (new, untracked): $rel"; fi
done < "$OUT/written.txt" > "$OUT/restore_log.txt"
git status --porcelain --untracked-files=all > "$OUT/status_after.txt"
cmp -s "$OUT/status_before.txt" "$OUT/status_after.txt" && echo "working tree: status identical to before" >> "$OUT/restore_log.txt" || { echo "working tree: STATUS DIFFERS" >> "$OUT/restore_log.txt"; diff "$OUT/status_before.txt" "$OUT/status_after.txt" >> "$OUT/restore_log.txt"; }
echo "DONE $(date '+%H:%M:%S')" >> "$OUT/times.txt"
