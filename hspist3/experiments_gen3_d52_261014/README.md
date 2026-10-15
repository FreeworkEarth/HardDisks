# The collision return-map panel in the viewer (plan-author decision 52 item 2 and addendum 1; notes 261012 sec. 4.7.116, 4.7.120)

<!-- ##CHRIS 2026-10-14 -->
Display only. A window run with the panel is a demonstration, never data.

## The code

- **edmd_core/edmd_gen3.[ch]: the collision gather hook** (edmd3_set_coll_hook).
  - It is off (NULL) unless the window's panel is open.
  - The hook and its state are kept outside struct EDMD3, so the state's size, the checkpoint format and every event are unchanged.
  - The three collision rules hand it a copy of what they computed: the disks, their compartments, the line of centres (or the direction of the wall or body), and the velocities after.
  - edmd3_destroy clears it.
- **00ALLINONE.c: the panel** (the block before keysSimulation; rm_frame and rm_draw in the window's frame; the help line).
  - **The keys:**
    - `I` opens and closes the panel;
    - `1` switches the angle (the line of centres, or the velocity after the hit);
    - `2` shows or hides the wall / divider / piston points (cyan);
    - `3` cycles the side (both, left, right);
    - `4` / `5` halve or double K (625 .. 1,280,000; default 20,000).
  - **The histogram:** 90 x 90 bins, log colour scale, axes 0 .. 2 pi with ticks at multiples of pi/3, and the number of points shown.
  - **The hook** is set only while the panel is open, only on the window's gen3 state, and again after every rebuild (g3_builds). Under gen2 the key is refused with a note.
- **00ALLINONE.c: the test action `hash` in --gui-keys.** It prints [GUI-HASH] with gen3's event hash, its time and its executed events at exactly that step, so that runs whose frames differ compare exactly.

## The tests (bin_00f4dd0e-dirty-c: the final code, uncommitted at its build)

- `tests/c_RESULTS.md` (scripts/tests_d52.sh, results_d52.py): the window's gen3 scene at N 400 with Chris's usual flags (sec. 4.7.107), eta 0.30 and 0.69.
  - **H:** gen3's event hash, time and events at exactly steps 600, 1200 and 1800 are EQUAL at both densities, for the panel closed, open from the start, toggled during the run (7 notes), and open with all its keys pressed (12 notes).
  - **F:** frames after step 3000 as PNG:
    - c_F_e030_centres (the fluid);
    - c_F_e069_centres (bands at multiples of pi/3);
    - c_F_e069_velocity (the second angle).
    - Each also as the panel cropped (_panel.png).
- `tests/c_T1` (tests_d50.py T1, TESTS_OUTROOT here): Part B's six registered stage 1 task lines EQUAL against bin_d2dc177 (event hash, t_end, every output file).
- **The smoke runs while building:**
  - `tests/s0_smoke` (the first frame, eta 0.30);
  - `tests/s1_smoke_b` (eta 0.69, after the text fixes);
  - `tests/s2_apply_probe` (the record at a scripted APPLY lands on a frame boundary, 1204 steps for 1200, hence the `hash` action).

## The clean build of the commit (a71f394d)

- `bin_a71f394`: build_clean.sh of a71f394d (a git archive export, kissfft's six files with their recorded local change). The binary's build line reads git a71f394, target mac-O3-e0pre; 25 -Wall warnings, as before. This is the binary for Chris's viewer.
- `tests/k_RESULTS.md`: the same tests on it.
  - H EQUAL at both densities, with hashes identical to bin_00f4dd0e-dirty-c's.
  - The frames k_F_*.png.
- `tests/k_T1`: Part B's six task lines EQUAL against bin_d2dc177.

