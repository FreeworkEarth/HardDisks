#!/bin/bash
# ##CHRIS 2026-10-13 (plan-author decision 48, item 3): the test scenes. Each starts from Chris's command of 20:55 (00_COMMAND.md)
# with the scene's N, radius, wall thickness and spring, builds gen3 at its first step (the startup build), then rebuilds it through
# the F3 panel's APPLY (--gui-keys STEP:apply-gen3: live_scene_sync_from_runtime, the engine row gen3, live_scene_prepare_rebuild,
# as the APPLY button) and ends with Q.  usage: tests_d48.sh <test name> ...  (each name once; a run is never repeated in place)
set -uo pipefail
SP=/private/tmp/claude-501/-Users-chrisharing-Desktop-CCS-complex-coupled-systems-Repo-HardDisks/91cb08ec-0599-4faa-a6af-d5a7ca834255/scratchpad
G=/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks_gen3/hspist3
OLD=$G/experiments_gen3_d34s_261012/bin_406044f/00ALLINONE
NEW=$G/experiments_gen3_d48_261013/bin_4ede56f7-dirty/00ALLINONE
BASE=(--show-simulation --mode=edmd --engine=gen3 --edmd-acc=0 --seed-drift-order=drift-first --kbt1 --l0=89.5904 --height=70.7107
      --num-walls=1 --wall-positions=89.5904 --wall-mass-factors=50 --wall-hold-steps=600 --fixed-dt=0.4 --seed=9300 --gui-time-scale=200
      --auto-release)   # the scene starts by itself (as the decision-34 acceptance runs; otherwise it waits for SPACE)
ET=(--wall-thickness=3.829 --spring-k-sigma=11.52 --spring-wall=0)          # the screenshot's wall thickness and k (the panel default)
S1=(--particles=341 --particles-boxes=0,341 --particle-radius=0.9079419)    # diameter 43.58 px, as in the banner
S2N=(--particles=341 --particles-boxes=0,341 --particle-radius=1.9146)      # radius 1.9146 with the same N (eta ~0.63)
S2E=(--particles=77 --particles-boxes=0,77 --particle-radius=1.9146)        # radius 1.9146 with eta ~0.139 (N 77)
S3=(--particles=341 --particles-boxes=0,341 --particle-radius=0.5)          # the default radius (R = 12 px)
KEYS=(--gui-keys=120:apply-gen3,6000:Q --gui-run-seconds=60)
run(){ bash "$SP/d48/run_scene.sh" "$@"; echo; }
for t in "$@"; do
  case $t in
    t0b_old_scene1)  run $t "$OLD" 90  -- "${BASE[@]}" "${ET[@]}" "${S1[@]}" --gui-keys=120:apply-gen3,6000:Q --gui-run-seconds=15 ;;
    t1_new_scene1)  run $t "$NEW" 240 -- "${BASE[@]}" "${ET[@]}" "${S1[@]}" "${KEYS[@]}" ;;
    t2_new_scene2_N341) run $t "$NEW" 240 -- "${BASE[@]}" "${ET[@]}" "${S2N[@]}" "${KEYS[@]}" ;;
    t2_new_scene2_eta) run $t "$NEW" 240 -- "${BASE[@]}" "${ET[@]}" "${S2E[@]}" "${KEYS[@]}" ;;
    t3_new_scene3_default) run $t "$NEW" 240 -- "${BASE[@]}" "${ET[@]}" "${S3[@]}" "${KEYS[@]}" ;;
    t3_old_scene3_default) run $t "$OLD" 240 -- "${BASE[@]}" "${ET[@]}" "${S3[@]}" "${KEYS[@]}" ;;
    t4_new_scene1_cellpx64) run $t "$NEW" 240 -- "${BASE[@]}" "${ET[@]}" "${S1[@]}" --gen3-cell-px=64 "${KEYS[@]}" ;;
    t5_new_scene1_headless) run $t "$NEW" 90 -- "${BASE[@]}" "${ET[@]}" "${S1[@]}" --headless --time=5 ;;
    t6_new_scene1_shot) run $t "$NEW" 240 -- "${BASE[@]}" "${ET[@]}" "${S1[@]}" --gui-window=2000x1200 \
                          --gui-shot=3000:"$SP/d48/runs/t6_new_scene1_shot/t6_scene1_step3000.bmp" "${KEYS[@]}" ;;
    t7_new_scene1_cellpx32) run $t "$NEW" 90 -- "${BASE[@]}" "${ET[@]}" "${S1[@]}" --gen3-cell-px=32 --gui-window=2000x1200 \
                          --gui-shot=1:"$SP/d48/runs/t7_new_scene1_cellpx32/t7_scene1_refused.bmp" --gui-keys=120:Q --gui-run-seconds=10 ;;
    *) echo "unknown test $t"; exit 1 ;;
  esac
done
