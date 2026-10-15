#!/bin/bash
# ##CHRIS 2026-10-14 (plan-author decision 52, addendum 1 item 3; 261012 sec. 4.7.120): THE RETURN-MAP PANEL'S TESTS in the window's gen3
# scene (run_scene_d52.sh: the SDL dummy driver, no window on the screen). Two scenes at N 400 with Chris's usual flags (the viewer
# command of 261012 sec. 4.7.107 with N, L0 and H of the stage 1 cells): eta 0.30 (L0 26.1799...) and eta 0.69 (L0 11.3826...), H 20.
#   H  equal event hashes: each scene with the panel closed, open from the start, toggled during the run, and open with all its keys
#      pressed (1 2 3 4 5); gen3's event hash, time and executed events printed at exactly steps 600, 1200 and 1800 (--gui-keys "hash");
#   F  frames: the panel open from the start, the frame after step 3000 (BMP, converted to PNG by results_d52.py): both scenes with the
#      line of centres, and eta 0.69 with the velocity after the hit (key 1 at step 2).
#   usage: tests_d52.sh <binary> <prefix>      (each run once; a name that exists is refused)
set -u
BIN=${1:?binary}; P=${2:?prefix}
E=/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks_gen3/hspist3/experiments_gen3_d52_261014
SP=/private/tmp/claude-501/-Users-chrisharing-Desktop-CCS-complex-coupled-systems-Repo-HardDisks/91cb08ec-0599-4faa-a6af-d5a7ca834255/scratchpad
R="$E/scripts/run_scene_d52.sh"
BASE=(--show-simulation --mode=edmd --engine=gen3 --edmd-acc=0 --seed-drift-order=drift-first --kbt1 --particles=400 --particles-boxes=200,200
      --particle-radius=0.5 --height=20 --num-walls=1 --wall-mass-factors=50 --wall-hold-steps=600 --fixed-dt=0.4 --seed=9300
      --gui-time-scale=200 --auto-release --gui-window=2000x1200 --gui-run-seconds=120)
CUT="600:hash,1200:hash,1800:hash"
for sc in e030:26.179938779914941 e069:11.382582078223889; do
  n=${sc%%:*}; L=${sc#*:}; GEO=(--l0=$L --wall-positions=$L)
  bash "$R" ${P}_H_${n}_closed  "$BIN" 300 -- "${BASE[@]}" "${GEO[@]}" --gui-keys=$CUT,2000:Q &
  bash "$R" ${P}_H_${n}_open    "$BIN" 300 -- "${BASE[@]}" "${GEO[@]}" --gui-keys=1:I,$CUT,2000:Q &
  bash "$R" ${P}_H_${n}_toggled "$BIN" 300 -- "${BASE[@]}" "${GEO[@]}" --gui-keys=1:I,300:I,500:I,600:hash,700:I,900:I,1200:hash,1300:I,1500:I,1800:hash,2000:Q &
  bash "$R" ${P}_H_${n}_keys    "$BIN" 300 -- "${BASE[@]}" "${GEO[@]}" --gui-keys=1:I,200:1,300:2,400:3,450:3,500:4,550:5,600:hash,700:1,800:2,900:3,1000:4,1100:5,1200:hash,1800:hash,2000:Q &
  wait
done
for f in e030:26.179938779914941:centres e069:11.382582078223889:centres e069:11.382582078223889:velocity; do
  n=${f%%:*}; rest=${f#*:}; L=${rest%%:*}; m=${rest#*:}; GEO=(--l0=$L --wall-positions=$L)
  K="1:I"; [ "$m" = velocity ] && K="1:I,2:1"
  bash "$R" ${P}_F_${n}_${m} "$BIN" 300 -- "${BASE[@]}" "${GEO[@]}" --gui-keys=$K,3200:Q --gui-shot=3000:$SP/d52/runs/${P}_F_${n}_${m}/frame.bmp
done
