# Part A on gen3, option (ii): the package (plan-author decision 50 item 3 d; notes 261012 sec. 4.7.113; registration sec. 4.8.4)

<!-- ##CHRIS 2026-10-14 -->
Nothing here runs on KOA before amendment A5 (decision 50 item 2). The binary needs commit 8f538a21 or later: --gen3-exact-box, --gen3-seeding=lattice and --gen3-held-readout on the held path.

- `pa_worker_koa.sh`: one held trajectory per call. It is the confinement worker's AF command plus the gen3 flags of option (ii) and --gen3-held-readout=T_eq, with the hold = T_eq + the window and the released tail. Its task line has 15 fields:

  `A (rel position folder) j s seed L0 H N_s x_div (T_eq steps) (window steps) (tail steps) (trace every) (log 0|1) (orientation default|0|turned)`

  - The orientation is the cell's pin (decision 53 item 2, A5 option b), passed as --gen3-lattice-orientation; `default` passes none.

  - The health rule is the stage 1 worker's, plus exactly one `[EDMD3-HELD] window` and one `[EDMD3-HELD-GEOMETRY]` line.
  - The build guard, the locks, the "done before" skip and the kept failures are the stage 1 worker's.
  - Outputs per position folder: red_(seed).csv, run_(seed).log, tr_(seed).csv, summary_(seed).csv, and ev_(seed).csv with log = 1.
- `pa_reduce.py`: the run's record lines to one row. It holds the geometry, the lattice per compartment, F_L and F_R over (T_eq, T_eq + window] recomputed from the exact counters, n_L, n_R, the held divider's v_div, dx_div and dW_div, the tail's T_L and T_R, and clean, hash and t_end. With a log it adds the log's forces and their difference from the readout. It fails closed.
- `pa_orient.py` and `pa_orient_261014.tsv`: the orientation pin (decision 53 item 2).
  - For every block B cell, the crystal orientation of the default lattice start at j = 0: 40 cells `0`, 24 `turned`.
  - The feasibility check of that pin at all five positions: 64 of 64 feasible. In 61 cells the pinned arrangement equals j = 0's everywhere; in 3 fluid cells the orientation is pinned but the row count follows the width.
  - pa_tasks.py reads the table.
- `pa_tasks.py`: a DRAFT task generator for A5.
  - What is registered: block B's cells with L0, H and N_s as in stage 1's task list; the stencil L0 + j dL with dL = sigma_x / 2 exactly; T_eq per sec. 4.8.3.1; the window of 5000 sigma-time; the tail of 20 sigma-time.
  - What A5 decides: the seed base, the sizes, the seed counts and the event-log subset.
  - `--sizes 100,400 --seeds 4` gives the 640 lines of the epsilon_0 pilot at N 100 and N 400 (no seed shared with stage 1).

The Mac tests and the rehearsal are in experiments_gen3_d50_261014/ (README there).
