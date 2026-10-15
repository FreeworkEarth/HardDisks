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
