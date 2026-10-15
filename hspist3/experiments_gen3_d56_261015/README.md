# experiments_gen3_d56_261015: the Stage 2 and D2 sizing (plan-author decision 56, task 2; Report D56-1)

##CHRIS 2026-10-15. Mac only, no KOA, no engine change. Planning numbers and information; nothing here is registered data. The notes:
0000_PLAN_OVERALL/ALL_MARKDOWNS/261012_paper1_confinement.md sec. 4.7.133 (the decision) and 4.7.134 (Report D56-1). Every command, time,
input hash and output hash is in RUN_RECORD.txt.

| task | script (scripts/) | output |
|---|---|---|
| 2 c, the runs | d56_pilot.py | tasks_d56.txt, data/, run_d56_pilot.out, worker_*.out |
| 2 c, the analysis | d56_c_analysis.py | c_analysis_b/ (c_runs.csv, c_cells.csv, c_snapshots.csv, c_summary.txt, fig_c1, fig_c2, fig_c3); c_analysis/ is the first run, kept and superseded: its orient_rec_0/1 were linear means of an angle that wraps at +-30 deg (the figures and the other tables are byte-identical) |
| 2 b | d56_b_lattice.py | b_lattice/ (REPORT.md, b_starts.csv, runs/) |
| 2 a | d56_a_timing.py | a_timing/ (a_fits.txt, a_stage1_runs.csv, a_stage2_table.csv and .md, fig_a1) |
| 2 d | d56_d_d2.py | d_d2/ (d_summary.txt, d2_table.csv and .md, d_fit_record.csv, d_runs_fit.csv, fig_d1, fig_d2) |

- data/stage1/gen3/B_N0900_eta0.690000/m_4500: the registered stage 1 task lines of (N 900, eta 0.69, M 4500) s = 0, 1, 2, rerun on this
  Mac with the merged main's build (experiments_gen3_d50_261014/bin_0544b7e8-dirty) through the stage 1 worker. Same protocol and seeds as
  KOA's runs, other trajectories (the executed events differ: 87.95 / 87.93 / 87.94 M here, 87.82 / 87.78 / 87.79 M on KOA).
- data/pilot/P2_N0900_eta0.760000/m_4500: the Stage 2 pilot cell, N 900, eta 0.76, M 4500 (alpha 5, the middle of the nine masses),
  s = 0, 1, seeds run_seed(20261015, 0, 4, s), the stage 1 run length (T_eq 5000, 200 periods of the binary's nu_pred).
- The reduced tables in data/ (red_nu.csv, acf_runs.npz, nud_runs.csv, s1_runs.csv, snap_runs.csv, runs_health.csv, failed_runs.csv,
  acf_mean_fit.csv) are written by the registered reducers (reduce_B.py, reduce_s1_koa.py), unchanged. The raw traces
  (wall_x_positions_*.csv) stay on the Mac (the repository ignores them); their SHA-256 are in RUN_RECORD.txt.
