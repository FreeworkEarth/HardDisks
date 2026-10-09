#!/bin/bash
# usage: run_set.sh before|after -- runs every paper script that loads raw data (and the derived-table scripts) on the
# current working tree, deterministic plotting; records stdout, stderr, exit codes and every file written (by timestamp).
SP=/private/tmp/claude-501/-Users-chrisharing-Desktop-CCS-complex-coupled-systems-Repo-HardDisks/91cb08ec-0599-4faa-a6af-d5a7ca834255/scratchpad
REPO=/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks
OUT=$SP/gguard/$1; mkdir -p "$OUT"
cd "$REPO/hspist3" || exit 1
export MPLBACKEND=Agg SOURCE_DATE_EPOCH=0 PYTHONHASHSEED=0
touch "$OUT/.marker"; sleep 1
run() { name=$1; shift; t0=$(date +%s); python3 "$@" > "$OUT/$name.out" 2> "$OUT/$name.err"; rc=$?; echo $rc > "$OUT/$name.rc"; echo "$name $(( $(date +%s) - t0 )) s rc $rc" >> "$OUT/times.txt"; }
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
find "$REPO" -newer "$OUT/.marker" -type f -not -path '*/.git/*' -not -path '*/__pycache__/*' -not -path '*/.mplconfig/*' | sort > "$OUT/written.txt"
while read -r f; do shasum -a 256 "$f"; done < "$OUT/written.txt" > "$OUT/written.sha"
echo "DONE $(date '+%H:%M:%S')" >> "$OUT/times.txt"
