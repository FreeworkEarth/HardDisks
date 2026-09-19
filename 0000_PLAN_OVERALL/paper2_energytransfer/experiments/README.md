# Paper 2 experiments — what is where

Copies made 2026-09-19; originals in `ALL_MARKDOWNS/260909_plots/` (frozen) and
`hspist3/experiments_energy_transfer/`.

| folder | holds |
|---|---|
| `final/` | `260917_level2_work_and_friction` (the running Green–Kubo integral and the work ladder), `260918_level2_A_of_u` (A(u) across five decades of speed plus the travel scan), and the two demo-mode screenshots `260919_demo_geomA/D.png` used as the plan's figure 1. |
| `run_scripts/` | copies of every batch script that produced the Level 0–2 data. Run them from `hspist3/`; they skip cells that already exist, so re-running resumes. |

## The raw data is not here

The trajectories live in `hspist3/experiments_energy_transfer/<campaign>/` and are **not** in git
(6.3 GB, gzipped traces). What is in git per cell is `summary.csv` (one row per run: the work) and
`00_COMMAND.md` (the exact command). Level 2 additionally has `stop_digest.csv`, an 85 kB extract
from which the whole analysis reproduces — verified by re-running it with every trace hidden.

## Levels and their state

| level | result | script |
|---|---|---|
| 0 ledgers | passed, 8e-12 (energy), 1e-13 (momentum) | `paper2_level0_level1_20260916.py` |
| 1 quasi-static work | W(0) − W_qs = +0.0043 ± 0.0158 kT, 0.3 σ, on the corrected geometry | `paper2_geometry_fix_20260918.py` |
| 2 slow end | no linear friction in a closed box; excess is 25.2 u² (step), 11.3 u² (ramp) | `paper2_level2_20260917.py`, `paper2_ramp_fast_20260918.py` |
| 2 fast end | Eq. 7 confirmed once the swept strip is counted properly | `paper2_ramp_fast_20260918.py` |
| 3–6 | not started | — |

