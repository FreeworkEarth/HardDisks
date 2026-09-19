# Paper 1 experiments — what is where

Everything here is a **copy**, made 2026-09-19. The originals are still in
`ALL_MARKDOWNS/260909_plots/`, which is frozen as of that date; new output goes here.

| folder | holds |
|---|---|
| `final/` | the canonical results. `260919_*` is the current answer: `260919_A1v2_final_cs_vs_eta.csv` (35 densities, N = 100) and `260919_A2_cs_per_mass.csv` (finite size), with the three figures `260919_cs_vs_eta`, `260919_cs_vs_eta_lowdensity_zoom`, `260919_cs_vs_eta_N100_vs_A2`. Also the A2 slow mode, the finite-size forms, the thickness-correction comparison figure and the linewidth JSON. |
| `estimator_tests/` | how the estimator was chosen, not results: the 260913 test A/B/C sweeps, every 260914 A1 v2 estimator variant (first 25/50/100 oscillations, k_min scans, per-trajectory fit, window-3 peak, velocity spectrum, wander equipartition) and the 260915 damping test. |
| `archive/` | superseded output from 260909–260913 — refit, FINAL, MEETING, PARTIAL, robust, with_uncorrected, interim, route A. Kept for provenance; do not cite. |
| `scripts/` | copies of the standalone plotting scripts. The originals live in `hspist3/validation/` and `hspist3/` and are what actually run. |

## Which file is canonical

**`final/260919_A1v2_final_cs_vs_eta.csv`.** It carries the divider thickness in the effective length,
L_eff = L0 − 2r − t/2 with t = 0.05 σ, which the 260914 version did not. The two differ by a factor
(L0 − 1 − 0.025)/(L0 − 1): nothing below η = 0.3, −0.4 % at η = 0.52, −0.5 % at 0.65.

## What regenerates what

| output | script (in `hspist3/validation/`) |
|---|---|
| `260919_*` CSVs and the three figures | `paper1_canonical_20260919.py` |
| the A1 v2 CSV from raw traces | `final_A1_figures_20260914.py` (reads 7875 trajectories) |
| A2 finite size | `analyze_A2_X2p5_20260914.py`, `A2_dilute50_20260917.py` |
| estimator comparison | `estimator_massladder_20260917.py`, `damping_test_20260915.py` |
| linewidth / Γ | `paper1_linewidth_20260918.py` |
| thickness comparison figure | `paper1_thickness_correction_20260918.py` |

The estimator and the geometry are defined once, in `tests_20260913.py`: `WALL_T = 0.05` and
`l_eff(L0) = L0 - 2r - t/2`, used by both the launcher and the analysis so they cannot diverge.

