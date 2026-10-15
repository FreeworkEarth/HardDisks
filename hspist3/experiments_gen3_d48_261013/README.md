# Decision 48 (261012 sec. 4.7.104): the viewer's gen3 cell width -- the build and the tests (2026-10-13; Report D48 accepted by decision 48 addendum 1, sec. 4.7.106; committed with the edit, local only until the plan author's word)

- **The commit:** the edit and this folder, committed together on engine-gen3, with the message the plan author gave in decision 48 addendum 1 item 2. The binaries are not committed, as in the earlier evidence folders.
- **The clean rebuild for Chris's viewer** (addendum 1 item 2) was made from that commit, d2dc1773, with experiments_gen3_m3_261009/build_clean.sh: **bin_d2dc177/00ALLINONE**. bin_4ede56f7-dirty/ stays.
  - It reports `git d2dc177  target mac-O3-e0pre`; its 00ALLINONE.c is byte-identical to the tested edit. -Wall -Wextra gives the same 25 warnings; the harnesses give 0.
  - The smoke test (tests/RESULTS_clean.md): scene 1 builds with a 44 px cell, 37,977 events after the APPLY, clean=1. Scene 3 builds with a 32 px cell, clean=1, with the same event hashes as the old viewer (08bf879143d890ed, bda502aee2a9339f).
  - The build's 188 MB source copy and its test binaries stay in the session scratchpad; bin_d2dc177/ holds the viewer and the text records.
- **The edit:** hspist3/00ALLINONE.c, made in the engine-gen3 worktree on HEAD 4ede56f7, two places marked ##CHRIS 2026-10-13 (decision 48):
  - g3_cell_arg, before g3_build: in the window's scene and without --gen3-cell-px, the larger of 32 px and the diameter rounded up to whole pixels; otherwise cli_gen3_cell_px, as before;
  - the status line (g3_gui_engine_line): the cell width next to "engine: gen3 (provisional display)".
  - The diff is in bin_4ede56f7-dirty/build_record.txt and in Report D48 (sec. 4.7.105).
- **bin_4ede56f7-dirty/:** the viewer, built with the command of experiments_gen3_d34s_261012/bin_406044f (the viewer Chris runs): -O3 -ffp-contract=off, target mac-O3-e0pre.
  - build_record.txt: the command, the SHA-256 of the sources and of the binary, the diff.
  - 00ALLINONE.Wall.txt: -Wall -Wextra, syntax only. The same 25 warnings as 406044f; none is new.
- **scripts/:**
  - build_d48.sh: the build;
  - run_scene.sh: one scene, with the SDL dummy driver (no window on the screen), in its own scratch folder with the font linked in, and a hard kill;
  - tests_d48.sh: the scenes;
  - results_d48.py: makes tests/RESULTS.md.
- **tests/:**
  - per run: `<name>.cmd` (time, binary, flags, exit code), `<name>.out` and `<name>.err`;
  - RESULTS.md;
  - two frames (PNG, from --gui-shot's BMP): t6 (scene 1 at step 3000) and t7 (scene 1 with --gen3-cell-px=32, refused).
- **How each scene runs:**
  - It starts from Chris's command of 20:55 (00_COMMAND.md), with the scene's N and radius, --wall-thickness=3.829 and --spring-k-sigma=11.52 --spring-wall=0.
  - --auto-release is added: the scene starts by itself, as in the decision-34 acceptance runs; otherwise it waits for SPACE. The wall is therefore released at once, not after 600 held steps.
  - gen3 is built at the first step (the startup scene, spring off). It is rebuilt through the F3 panel's APPLY at step 120 (--gui-keys 120:apply-gen3, which does what the APPLY button does with the engine row on gen3; spring on, k 11.52).
  - Q ends the scene at step 6000, 98 sigma-time after the APPLY.
- **Two attempts that tested nothing, kept** (nothing is deleted; their names were not reused):
  - T0_old_scene1_r0.9079: zsh passed the flags as one argument ("Unknown option");
  - t0_old_scene1: without --auto-release the scene waited for SPACE; it was killed after 90 s with no gen3 build.
- **The scenes' CSV files** (energy_log.csv and the others) were written into the run folders in the session scratchpad, not here.

## The results (tests/RESULTS.md)

| test | N | radius [σ] (R [px]) | η (APPLY) | cell width of the builds [px] | events before APPLY / after | clean before / after | worst gap [px] | event hash before / after | exit |
|---|---|---|---|---|---|---|---|---|---|
| t0b_old_scene1 | 341 | 0.9079 (21.791) | - (no APPLY) | refused: cell width 32 px < diameter 43.581211090087891 px | 0 (no gen3 state) | 0 (no gen3 state) | - | - | 0 |
| t1_new_scene1 | 341 | 0.9079 (21.791) | 0.139 | 44, 44 | 722 / 37,981 | 1 / 1 | -9.52e-13 | 7daa71fd3902f4bb / bbc410cb626cf101 | 0 |
| t2_new_scene2_eta | 77 | 1.9146 (45.950) | 0.140 | 92, 92 | 84 / 4,162 | 1 / 1 | -7.67e-13 | bcc140c53acb89a2 / 17f5cb2a5c0bd428 | 0 |
| t2_new_scene2_N341 | 341 | 1.9146 (45.950) | 0.613 | 92, 92 | 1,898 / 73,270 | 1 / 1 | -9.95e-13 | 183cdfc3546be0f9 / f770a6f5305e207e | 0 |
| t3_new_scene3_default | 341 | 0.5000 (12.000) | 0.042 | 32, 32 | 878 / 44,308 | 1 / 1 | -9.09e-13 | 08bf879143d890ed / bda502aee2a9339f | 0 |
| t3_old_scene3_default | 341 | 0.5000 (12.000) | 0.042 | 32, 32 | 878 / 44,308 | 1 / 1 | -9.09e-13 | 08bf879143d890ed / bda502aee2a9339f | 0 |
| t4_new_scene1_cellpx64 | 341 | 0.9079 (21.791) | 0.139 | 64, 64 | 503 / 28,382 | 1 / 1 | -9.73e-13 | 5f0e12cf5c74bf2a / 5a32c0fcc434a936 | 0 |
| t5_new_scene1_headless | 341 | 0.9079 (21.791) | - (no APPLY) | refused: cell width 32 px < diameter 43.581211090087891 px | - | - | - | - | 2 |
| t6_new_scene1_shot | 341 | 0.9079 (21.791) | 0.139 | 44, 44 | 718 / 37,997 | 1 / 1 | -9.52e-13 | 834fc942cc4dea62 / d703105942853178 | 0 |
| t7_new_scene1_cellpx32 | 341 | 0.9079 (21.791) | - (no APPLY) | refused: cell width 32 px < diameter 43.581211090087891 px | 0 (no gen3 state) | 0 (no gen3 state) | - | - | 0 |

- "old" is experiments_gen3_d34s_261012/bin_406044f/00ALLINONE; "new" is bin_4ede56f7-dirty/00ALLINONE.
- Every safety counter of every gen3 state is 0 (the .out files' [EDMD3-HEALTH] lines). The worst gap is the most negative surface gap the validator saw, at the rounding level.
- t5 runs with --headless: the scene stops as before (STOP, exit 2), because the new width applies only in the window.
- The t3 pair: the same events from the old and the new viewer at the default radius (equal event hashes at both builds).
