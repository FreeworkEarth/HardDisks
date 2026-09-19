# Run scripts for the Paper 2 energy box

Every script here is the exact one that produced the data quoted in
`0000_PLAN_OVERALL/ALL_MARKDOWNS/260916_paper2_plan_energy_box_stepbystep_COWORK.pdf`
and in the Level 0/1/2 reports. Run them from `hspist3/`:

```sh
cd hspist3
bash experiments_energy_transfer/_run_scripts/<script>.sh
```

They skip any run whose trace already exists, so re-running one is safe and resumes.

| script | what it produces | runs | wall-clock (M-series Mac, 10 parallel) |
|---|---|---|---|
| `level1_Zwall_path_part1.sh` | `level1_Zwall_path_20260916/L*` — held-wall pressure Z at path points | 40 | ~10 min |
| `level1_Zwall_path_part2.sh` | the two extra path points that densify it to five | 200 | ~40 min |
| `level1_more_work_seeds.sh` | `level1_moreseeds_20260916/u*` — slow-push work, 50 seeds/speed | 150 | ~30 min |
| `level2_speed_ladder_part1.sh` | `level2_slope_20260917/u{0.02,0.03,0.05,0.10}` | 330 | ~25 s |
| `level2_speed_ladder_part2.sh` | the same ladder extended to u = 0.005, 0.01, 0.15, 0.20 | 380 | ~5 min |
| `level2_ramp_vs_step_and_fast.sh` | `level2_ramp_20260918/*` (ramp + step control) and `level2_fast_20260918/*` | 660 | ~40 s |
| `level2_fast_1sigma.sh` | `level2_fast1sigma_20260918/*` — the fast end at travel 1 sigma | 300 | ~12 s |
| `level2_equilibrium_baseline.sh` | `level2_ramp_20260918/step_u0.005` — equilibrium profile for the decomposition | 60 | ~3 min |

**Which binary.** Everything except the ramp uses the installed `./00ALLINONE`. The ramp mode
(`--piston-right-protocol-mode=ramp`) and the stop snapshot (`HD_STOP_SNAPSHOT`) live in
`./00ALLINONE_ramp`, built with `make release TG=00ALLINONE_ramp`; in step mode it is byte-identical
to the installed binary, which is the gate documented in `260919_paper2_ramp_linewidth_fast_REPORT.md`.

**Environment variables used.** `HD_PISTON_EVENTS=<path>` writes the per-event wall impulse log
(needed for every ledger and for zeta); `HD_STOP_SNAPSHOT=<path>` writes one particle snapshot at the
instant the piston stops. Both are logging only: without them not a byte of any other output changes.

**Every run directory also records its own command** in `00_COMMAND.md`, written by the binary, so
any single trajectory can be reproduced without reading these scripts at all.
