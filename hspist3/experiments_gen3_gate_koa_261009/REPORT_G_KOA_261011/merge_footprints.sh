#!/bin/bash
# ##CHRIS 2026-10-11 (plan-author decision 26, point 3; 261012 sec. 4.7.54): whether merged/.merge_complete is in the fetched tree, and the
# merge's footprints there. The registered fetch (fetch_gen3_koa.sh) copies merged/testG/ only, so .merge_complete (in merged/) is never
# fetched. merge_testG_koa.py (fd8d8a2) writes, in this order: the moves; .build_git in every merged cell; merged/.merge_complete; then
# reduce_B (red_nu.csv, acf_runs.npz per m_<M>) and the compression (the SHA-256 of each trace appended to .sha256_uncompressed, then gzip).
# usage: bash merge_footprints.sh <the fetched folder experiments_gen3_gate_koa_261009>
cd "$1/testG" || exit 2
echo "merged/.merge_complete in the fetched folder: $(find .. -name '.merge_complete*' | wc -l | tr -d ' ') file(s)"
echo "| cell | ##RUN sections in run.log | trace SHA lines in .sha256_uncompressed | uncompressed traces fetched | red_nu.csv rows | .build_git | KOA mtime of .build_git / red_nu.csv / .sha256_uncompressed |"
echo "|---|---|---|---|---|---|---|"
for c in */*/m_*; do
  printf '| %s | %s | %s | %s | %s | %s | %s / %s / %s |\n' "$c" "$(grep -c '^##RUN' $c/run.log)" "$(grep -c 'wall_x_positions_' $c/.sha256_uncompressed)" \
    "$(ls $c | grep -c '^wall_x_positions_.*\.csv$')" "$(($(wc -l < $c/red_nu.csv) - 1))" "$(cat $c/.build_git)" \
    "$(stat -f '%Sm' -t '%m-%d %H:%M' $c/.build_git)" "$(stat -f '%Sm' -t '%H:%M' $c/red_nu.csv)" "$(stat -f '%Sm' -t '%H:%M' $c/.sha256_uncompressed)"
done
for e in gen2 gen3; do x=$e/epi8_H_H10_L10/x_0
  echo "$x: red_ files $(ls $x | grep -c '^red_'), run_ logs $(ls $x | grep -c '^run_[0-9]'), the merge's per-chunk 00_COMMAND_chunk<j>.md for chunks $(ls $x | sed -n 's/^00_COMMAND_chunk\([0-9]*\)\.md$/\1/p' | sort -n | tr '\n' ' '), run_params_chunk<j>.json $(ls $x | grep -c '^run_params_chunk'), ev SHA lines $(grep -c 'ev_' $x/.sha256_uncompressed), .build_git '$(cat $x/.build_git)'"
done
