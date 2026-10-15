#!/bin/bash
# ##CHRIS 2026-10-14 (decision 50, the report): the whole test suite on the clean build of the final code (bin_8f538a2), in turn
set -u
SP=/private/tmp/claude-501/-Users-chrisharing-Desktop-CCS-complex-coupled-systems-Repo-HardDisks/91cb08ec-0599-4faa-a6af-d5a7ca834255/scratchpad
cd /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks_gen3/hspist3 || exit 1
NEW=experiments_gen3_d50_261014/bin_8f538a2/00ALLINONE; OLD=experiments_gen3_d48_261013/bin_d2dc177/00ALLINONE
run(){ local lab=$1; shift; python3 "$SP/d50/tests_d50.py" "$lab" "$@" > "$SP/d50/$lab.out" 2>&1; echo "$lab exit $?: $(tail -1 "$SP/d50/$lab.out")"; }
run f1_T1 T1 "$NEW" "$OLD"
run f2_T1AF T1AF "$NEW" "$OLD"
run f3_T2 T2 "$NEW" --seeder=lattice
run f4_C2 C2 "$NEW"
run f5_R0 R0 "$NEW"
run f6_R1 R1 "$NEW"
