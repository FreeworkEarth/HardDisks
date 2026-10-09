#!/bin/bash
# ##CHRIS 2026-10-09 (261012 sec. 4.7.11, decision 6): the outputs that read paper1_kr_sanity_261002.KR, before and after the
# x^57 coefficient of the rho_max = 0.89 fit is corrected. The sanity script in its three reference modes, and the KR range
# table of decision 3 (it reads KR[0.89]). paper1_window_explore_261008 reads only the 0.88 fit; paper1_draft_audit_20261014
# names the script in a string only. usage (from hspist3/): bash experiments_kr_sanity_fix_261009/run_outputs.sh <before|after>
set -u
D=experiments_kr_sanity_fix_261009/$1; mkdir -p "$D"
python3 validation/paper1_kr_sanity_261002.py            > "$D/kr_sanity_default.txt" 2>&1; echo "kr_sanity default rc $?"
python3 validation/paper1_kr_sanity_261002.py --ref 0.88 > "$D/kr_sanity_ref088.txt"  2>&1; echo "kr_sanity --ref 0.88 rc $?"
python3 validation/paper1_kr_sanity_261002.py --ref 0.89 > "$D/kr_sanity_ref089.txt"  2>&1; echo "kr_sanity --ref 0.89 rc $?"
python3 validation/paper1_kr_ranges_261009.py            > "$D/kr_ranges.txt"         2>&1; echo "kr_ranges rc $?"
