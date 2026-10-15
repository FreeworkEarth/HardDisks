# Option (ii) for Part A: the held path on gen3 (plan-author decision 50 item 3; notes 261012 sec. 4.7.113, decisions 51-52 sec. 4.7.114-4.7.116)

<!-- ##CHRIS 2026-10-14 -->
The held path is the energy-transfer experiment (the worker's AF command, cluster/confinement_20261013/conf_worker.sh). Option (ii) gives it, under --engine=gen3, the registered geometry and start, and the forces from the engine's per-face impulse counters.

The four steps, one local commit each:
- (a) The exact box and an exact double divider position (--gen3-exact-box). The box and the position are checked at 1e-9 sigma at initialization, and the held geometry is printed as [EDMD3-HELD-GEOMETRY] for gate C2.
- (b) The gen3 lattice start (g3_m4_seed) in the held geometry.
- (c) F_L and F_R from edmd3_divider_impulse / edmd3_divider_events at T_eq and at the end of the hold, in the run record.
- (d) The worker (t0 = T_eq, the hold, the tail) and a reducer of the record lines.

## scripts/

- `build_d50.sh`: a Mac build of the worktree into `bin_(short hash)[-dirty]`, with the command of experiments_gen3_d48_261013/scripts/build_d48.sh (-O3 -ffp-contract=off, target mac-O3-e0pre). The build record holds the sources' SHA-256 and the diff.
- `tests_d50.py`: the tests.
  - T1: registered stage 1 task lines with the worker's command, under OLD and NEW. Equal means the same event hash, the same t_end and the same SHA-256 for every output file.
  - T1AF: the legacy held path (no exact box), with the event log, compared the same way.
  - T2: the exact geometry at the N 400 stencil (16 cells x 5 positions) and at pi/8 for N 100, 900 and 1600.

## Step (a)

- Build: `bin_109b0057-dirty` (00ALLINONE.c of the step on 109b0057). The first line of its build record was rewritten by hand after a quoting error of the build script, which is fixed since; the rest is as written.
- `tests/a1_T1`: 6 of 6 EQUAL against bin_d2dc177 (the clean viewer build; its sources equal 109b0057's).
- `tests/a2_T1AF`: EQUAL (the event hash, t_end, the trace and the event log).
- `tests/a3_T2`, with the legacy seeder: 8 of 83 PASS, the on-grid cells (pi/8 at every N, the five N 400 positions).
  - The 75 off-grid starts stop at the validator's initialization check (particle_wall_penetration). The legacy seeder places the divider as a fraction of the 1/48-grid box, up to 0.47 px from the exact divider.
  - Step (b)'s lattice start is the registered start, and it runs every cell.

Committed: each run's cmd.txt and stdout.log (added with -f, since *.log is ignored) and each test's RESULTS.md. The runs' data files are not committed; RESULTS.md compares their SHA-256.

## Step (b)

- Build: `bin_bcb0e3b1-dirty` (00ALLINONE.c of the step on bcb0e3b1, step a's commit).
- The code change:
  - The M4 start's body takes its geometry (g3_m4_init_core(boxW, c, seed)); g3_m4_init_driver passes the speed-of-sound loop's geometry, computed as before.
  - The held path with --gen3-exact-box --gen3-seeding=lattice seeds it in the box 2 L0 with the divider at the requested position. It calls srand(seed) right before the start, as the speed-of-sound loop does per trajectory.
  - The lattice's doubles replace the loaded positions. --gen3-dump-initial and the run header [EDMD3-RUN] are written as in the speed-of-sound loop.
  - The guard: in the held path, --gen3-exact-box and --gen3-seeding=lattice go together, and --seed is required.
- `tests/b1_T1`: 6 of 6 EQUAL against bin_d2dc177.
- `tests/b2_T1AF`: EQUAL.
- `tests/b3_T2`, with the lattice start: 83 of 83 PASS. Box and divider are within 7.1e-15 sigma, N_L = N_R = N_s, and every run is clean.
- `tests/b4_C2` (gate C2): 6 of 6 EQUAL. The held start at j = 0 and the free start with the same seed give byte-identical initial dumps (the box, the divider, every disk's x, y, vx, vy in %a), at pi/8 for N 100, 400, 900 and 1600 and at eta 0.05 and 0.69 for N 400. The dumps are committed.

## Step (c)

- Build: `bin_dda7f4a6-dirty` (00ALLINONE.c of the step on dda7f4a6, step b's commit).
- The code change: --gen3-held-readout=K, accepted on the held path on the exact box only.
  - It reads gen3's per-face counters of the divider (edmd3_divider_impulse, edmd3_divider_events), its state (edmd3_divider_state), its work (edmd3_work_divider) and the clock (edmd3_time). The reads are pure.
  - Read at step K (T_eq) and at the last step of the hold, before the release.
  - Printed: `[EDMD3-HELD] t0` and `t1` (cumulative values, also in %a), then `[EDMD3-HELD] window` with T_w, F_L = dJ_L/T_w, F_R = -dJ_R/T_w, n_L, n_R, dx_div, v_div at t0 and t1, and dW_div.
- `tests/c1_T1`: 6 of 6 EQUAL against bin_d2dc177.
- `tests/c2_T1AF`: EQUAL.
- `tests/c3_R0`: 5 of 5 EQUAL. The same held runs with and without the readout give the same event hash, t_end and output files.
- `tests/c4_R1`: 5 of 5 PASS. The readout against the event log (reduce_AF.py's formulas on the D0 rows of the readout's window):
  - F_L and F_R agree to at most 2.6e-13 relative (the log prints 12 digits);
  - n_L and n_R agree exactly;
  - the log's u_wall_max and W_div are 0, and the readout's v_div, dx_div and dW_div are 0.
  - The event logs are not committed.
- Note: a step advances (double)0.4f = 0.40000000596 internal units, so 300,000 steps are 5000.0000745 sigma-time. The readout divides by the engine's own t1 - t0.

## Step (d)

- The package: cluster/gen3_partA_261014/ (pa_worker_koa.sh, pa_reduce.py, the draft pa_tasks.py; its README).
- The clean build of the final code: `bin_8f538a2` (build_clean.sh of 8f538a21, a git archive export; the binary's build line reads git 8f538a2). Step (d) changes no C source.
- The whole suite on bin_8f538a2 (`scripts/run_final_tests.sh`):
  - `tests/f1_T1`: 6 of 6 EQUAL; `tests/f2_T1AF`: EQUAL;
  - `tests/f3_T2`: 83 of 83 PASS; `tests/f4_C2`: 6 of 6 EQUAL;
  - `tests/f5_R0`: 5 of 5 EQUAL; `tests/f6_R1`: 5 of 5 PASS.
  - The values are those of the per-step builds.
- `tests/f7_LAT` (information): the lattice start of every block B cell at the five positions. 21 of 64 cells change the arrangement across the stencil.
  - Nearly all of them also change a compartment's crystal orientation.
  - The dense ones are N 400 at eta 0.66 and 0.68, with one compartment turned at j = +-2. The rest are eta <= 0.55.
- `scripts/rehearse_pa_mac.sh`, the worker rehearsed on the Mac with bin_8f538a2:
  - `tests/d3_rehearsal`:
    - R1: 40 short tasks ok (N 100 and 400, four densities, five positions; 8 with the event log, which agrees with the readout to 3e-13).
    - R2: the same 40 lines are skipped as "done before".
    - R3: a different binary is refused by the build guard (exit 3).
    - R4: a divider outside the box is refused by the binary; the worker reports FAILED and keeps the folder.
    - R5: two tasks at the registered lengths, N 100 eta 0.05 j = -2 in 1 s and N 400 eta 0.69 j = +1 in 63 s, with T_w = 5000.0000745 sigma-time.
  - `tests/d1_rehearsal` (the data root was given as a relative path; the worker refused it, so nothing ran).
  - `tests/d2_rehearsal` (R5 ran into R1's root; its tasks share R1's seeds, so they were skipped as "done before"; fixed by giving R5 its own root).
  - Both are kept as they ran.

## The orientation pin (decision 53 item 2, A5 option b)

- No C source changes: --gen3-lattice-orientation (M4c) already restricts the lattice fit to one crystal orientation.
- The package: pa_orient.py and its table pa_orient_261014.tsv; the 15th task field (orientation) in pa_worker_koa.sh, pa_reduce.py (columns pin, pin_left, pin_right; it fails closed unless [EDMD3-M4C] reports the pin in both compartments) and pa_tasks.py.
- `tests/g1_pin` (pa_orient.py with bin_8f538a2):
  - every block B cell's orientation at j = 0 (40 cells `0`, 24 `turned`; the two compartments always agree);
  - the pinned start at the five positions: 64 of 64 cells feasible, no problem;
  - in 61 cells the pinned arrangement equals j = 0's everywhere;
  - in 3 fluid cells the row count still follows the width with the orientation pinned: N 100 eta 0.55, N 400 eta 0.25, N 900 pi/8.
- `tests/g2_C2P`: gate C2 with the pin, 6 of 6 EQUAL. The dumps equal the unpinned C2's: at j = 0 the pin is the default arrangement.
- `tests/g3_pin_rehearsal` (`scripts/rehearse_pin_mac.sh`):
  - 40 pinned 15-field tasks ok (N 100 and 400 at pi/8, 0.05, 0.66, 0.68), and one unpinned line ok;
  - every row's compartments carry the pin;
  - the dense N 400 cells (0.66, 0.68) start with one arrangement at all five positions.

