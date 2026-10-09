#!/bin/bash
# ##CHRIS 2026-10-09 (261012 sec. 4.7.8, item 3): cross-checks of run_edr_before_after.sh's output (information; read-only).
# 0. blob identity of the touched files; 1. the 45 written files against decision 2's committed hashes; 2. stdout/stderr against decision 2's "before" run (scratch);
# 3. the written files against the committed blobs; 4. the T-prime script's output against the text recorded in sec. 4.4.15;
# 5. the guard on every directory of the fetched resched data (Test T, T-prime, ASan, gate runs) that holds a run.log.
SP=/private/tmp/claude-501/-Users-chrisharing-Desktop-CCS-complex-coupled-systems-Repo-HardDisks/91cb08ec-0599-4faa-a6af-d5a7ca834255/scratchpad
REPO=/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks
E=$REPO/hspist3/experiments_loader_guard_edr_261009
cd "$REPO" || exit 1
echo "0. blob identity (git rev-parse, 10 digits) of the files the paper scripts run that the guard commit touches or depends on:"
echo "| file (hspist3/...) | engine-divider-resched da89d39 (before) | main 0b21270~1 (before the guard) | main 524f21e (with the guard) |"
echo "|---|---|---|---|"
for p in validation/estimator_massladder_20260917.py validation/level2_Au_figure_20260918.py validation/paper1_A2_boxtrunc_261002.py \
  validation/paper1_boxtrunc_20261014.py validation/paper1_confinement_afix_261005.py validation/paper1_confinement_heldwall_posthoc_261004.py \
  validation/paper1_confinement_prereg_20261012.py validation/paper1_confinement_results_261004.py validation/paper1_populate_cs_err_20261002.py \
  validation/paper2_figures_20261001.py validation/paper2_geometry_fix_20260918.py validation/paper2_ramp_fast_20260918.py \
  validation/resched_gate_261005.py validation/tests_20260913.py validation/edmd_acc_guard.py validation/test_edmd_acc_guard.py \
  validation/provenance_edmd_acc_261009.py cluster/confinement_20261013/reduce_B.py; do
  b() { git rev-parse --short=10 "$1:hspist3/$p" 2>/dev/null || echo "(absent)"; }
  echo "| $p | $(b da89d39) | $(b 0b21270~1) | $(b 524f21e) |"; done
echo "1. the written files (before and after) against decision 2's committed hashes (experiments_loader_guard_261008/written_before.sha256):"
for s in before after; do cmp -s "$E/written_$s.sha256" hspist3/experiments_loader_guard_261008/written_before.sha256 && r=identical || r=DIFFERENT
  echo "   $s: $(wc -l < "$E/written_$s.sha256" | tr -d ' ') files, $r"; done
echo "2. stdout and stderr of the branch's code (edr_before) against decision 2's before run of 2026-10-08 (main's pre-guard code):"
n=0; for f in "$SP"/gguard/before/*.out "$SP"/gguard/before/*.err; do b=$(basename "$f")
  if ! cmp -s "$f" "$SP/gguard/edr_before/$b"; then n=$((n+1)); echo "   differs: $b"; diff "$f" "$SP/gguard/edr_before/$b" | sed 's/^/     /'; fi; done
echo "   $n of $(ls "$SP"/gguard/before/*.out "$SP"/gguard/before/*.err | wc -l | tr -d ' ') files differ"
echo "3. the written files against the committed blobs (git HEAD):"
while read -r h p; do if git cat-file -e HEAD:"$p" 2>/dev/null; then
    [ "$(git show HEAD:"$p" | shasum -a 256 | cut -d' ' -f1)" = "$h" ] && echo "same as committed: .${p##*.}" || echo "differs from committed: .${p##*.}"
  else echo "not tracked: .${p##*.}"; fi; done < "$E/written_before.sha256" | sort | uniq -c | sed 's/^ */   /'
echo "4. the T-prime script's stdout (before = after) against the output recorded in 261012 sec. 4.4.15:"
F=0000_PLAN_OVERALL/ALL_MARKDOWNS/261012_paper1_confinement.md
s=$(grep -n 'Test T-prime, printed by `cd hspist3 && python3 validation/resched_testTprime_261007.py` on the fetched data, verbatim' $F | cut -d: -f1)
awk -v s="$s" 'NR>s && /^```/ {c++; if (c==2) exit; next} NR>s && c==1 {print}' $F | cmp -s - "$SP/gguard/edr_after/tprime.out" && r="identical, byte for byte" || r=DIFFERENT
echo "   $(wc -l < "$SP/gguard/edr_after/tprime.out" | tr -d ' ') lines: $r"
echo "5. the guard on the fetched resched data (hspist3/experiments_resched_gate2_261007; each directory that holds a run.log):"
(cd hspist3 && python3 - <<'PY'
import glob, os, sys
sys.path.insert(0, "validation"); import edmd_acc_guard as G
ds = sorted({os.path.dirname(p) for p in glob.glob("experiments_resched_gate2_261007/**/run.log", recursive=True)})
bad = []
for d in ds:
    try: G.guard(d)
    except G.AcceleratedRunError as e: bad.append(str(e))
groups = {}
for d in ds:
    k = "Test T-prime" if "resched_testTprime_261007" in d else "Test T" if "resched_testT_261007" in d else d.split("/")[1]
    groups[k] = groups.get(k, 0) + 1
print("   directories: " + ", ".join(f"{k} {v}" for k, v in sorted(groups.items())) + f"; accepted {len(ds) - len(bad)}, refused {len(bad)}")
for b in bad: print("   " + b)
PY
)
