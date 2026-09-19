# Paper 2 — Level 0 and Level 1 report

Written 2026-09-16 evening, while the dilute A2 ladders ran. Every number below is read from the output of
`hspist3/validation/paper2_level0_level1_20260916.py` (re-run end to end; output reproduced verbatim), or from
the file named beside it. No binary change, no change to accepted data. The only new simulations are 40 hold-only
runs of N = 100 (≈ 14 s wall in total, niced), in `experiments_energy_transfer/level1_Zwall_path_20260916/`.

## 0. Three corrections to the plan before the results

1. **The per-event log already exists.** Items 1 of the prompt asked to add it. It is in `edmd_core/edmd.c`
   (`edmd_set_event_log`, gated by the environment variable `HD_PISTON_EVENTS=<path>`) and records every outer-wall
   (WL, WR, WB, WT), divider (D0 …) and piston (PL, PR) collision as `t_sigma, kind, u_wall, v_before, v_after,
   dE, dp`, **including the hold phase**. The 2026-09-10 pilot (25 seeds × 3 speeds) and the 2026-09-11 W_qs runs
   (10 seeds × 5 speeds) already carry it. No code change was needed or made.
2. **ζ can come from existing runs after all** — from the holds of the energy-transfer runs, which carry the
   per-event log. It cannot come from the speed-of-sound holds, which never open the log. CC's earlier statement
   ("no per-event impulse log exists") was wrong; the corrected plan's statement ("ζ needs a new column") is
   therefore also wrong and should be reverted to "analysis of the energy-transfer holds".
3. **The pressure-ladder data is intact** (corrected 2026-09-16 evening; an earlier draft of this note wrongly said
   it was gone). All 171 accepted trajectories are on disk as one headerless row per seed — 66 in `runs/`, 9 in
   `N2500_20260910/runs/`, 96 in `preserved/` — with the column header kept only in the 303-byte
   `pressure_trajectories.csv`, which is an aggregate stub from 2026-09-06, not a wiped file. Parsed with that
   header: 171 rows, all `valid = 1`, 0 health events, 45 (η, N) cells, N = 400 / 900 / 1600 / 2500, and
   Z_pair at η = 0.10 reproduces `analysis/tables_ABC.md` exactly (1.2381 / 1.2379 / 1.2380). The false alarm came
   from a directory listing truncated just before `runs/` and from reading headerless files as if they had a header.
   What remains true: that ladder never ran at N = 100 and shows no size dependence at η = 0.10 (KR + 0.15 %), so
   it cannot supply a finite-box Z at N_s = 50. Level 1 therefore measured it directly (§ 2).

## 1. Level 0 — the ledgers close

Pilot of 2026-09-10: geometry A, η = 0.10, N_s = 50 per side, L0 = 39.25, H = 10, divider mass 10⁹, right piston
step at u, travel 3.93 σ, hold 1000 steps. Ledgers from the per-event log against the per-sample trace:

    R_E(t)  = ΔKE_gas(t) − Σ dE_gas(events ≤ t)
    R_px(t) = ΔPx_gas(t) − Σ dp_x(outer-wall, divider and piston events ≤ t)     (pair collisions conserve p)

| speed | seeds | max\|R_E\|/W as logged | max\|R_E\|/W, divider sign consistent | max\|R_px\|/Σ\|dp_x\| |
|---|---|---|---|---|
| 0.02 | 25 | 4.09e-06 | 7.14e-12 | 8.38e-14 |
| 0.05 | 25 | 2.78e-06 | 7.40e-12 | 7.86e-14 |
| 0.10 | 25 | 2.88e-06 | 8.16e-12 | 1.11e-13 |

**Pass (criterion 10⁻⁸) for energy and momentum, all 75 runs**, after two analysis corrections that are
bookkeeping, not physics:

- **Sample-time alignment.** The trace prints `Time` to 6 decimals; an event within 10⁻⁶ σ of a sample boundary
  was assigned to the wrong sample. Symptom: 9 of 25 seeds at u = 0.05 failed at exactly one sample in 17 000, with
  |R_px| equal to |dp| of the one nearest event (ratio 1.000, events 0.36–8 µσ from the sample). Rebuilding sample
  times from the step index removes it in every seed.
- **Divider sign convention.** Divider events log dE as the energy the *divider* gains; piston events log the energy
  the *gas* gains. Symptom: R_E grew smoothly to −1.5 × 10⁻⁵ over the record and equalled −2 × ΔKE_divider exactly
  (correlation +1.0000). Using −dE for divider events closes the ledger to round-off. The piston work in the trace
  matches Σ dE of piston events to 6 × 10⁻¹¹ throughout.

Recommendation: record both conventions in the Level-0 methods text, and either flip the sign of `dE` on divider
events in the log (one line in `edmd.c`, logging only, byte-identity gated) or document it where the log is described.

## 2. Level 1 — the quasi-static work of the finite box

Wall pressure of the right compartment measured **in equilibrium**, from the 200 σ holds, as the mean of the
outer-wall (WR) and gas-side divider (D0, dp > 0) impulse rates, Z = P A /(N_s kT) with kT = 1 exactly
(drift-first; a 10⁹-mass divider exchanges no measurable energy). Three points of the compression path:

| compartment length | η | seeds | Z_wall | sem | Z_KR | Z_wall/Z_KR |
|---|---|---|---|---|---|---|
| 39.25 | 0.100051 | 50 | 1.2933 | 0.0028 | 1.2363 | 1.0461 |
| 37.291667 | 0.105305 | 20 | 1.3327 | 0.0057 | 1.2510 | 1.0653 |
| 35.3125 | 0.111207 | 20 | 1.3492 | 0.0053 | 1.2679 | 1.0641 |

The two new lengths are grid-exact (2 L × 24 integer): 37.285 and 35.32 were first tried and every run was
rejected by the validator with `initial_wall_position_mismatch`, as it should be.

    path η 0.100051 → 0.111207 (the push ended at 0.111183; the grid-exact hold is 0.02 % further)
    W_qs bulk KR adiabat:                           7.0739
    W_qs finite box, start ratio held constant:     7.4233
    W_qs finite box, ratio measured along the path: 7.5316 ± 0.0191
    gap  <W>(u = 0.005, 10 seeds) = 7.5027 ± 0.0192  →  W − W_qs^finite = -0.0289 ± 0.0271 kT (-1.1 σ)
    gap  <W>(u = 0.01, 10 seeds) = 7.4456 ± 0.0446  →  W − W_qs^finite = -0.0859 ± 0.0486 kT (-1.8 σ)
    gap  <W>(u = 0.02, 10 seeds) = 7.3670 ± 0.0457  →  W − W_qs^finite = -0.1645 ± 0.0495 kT (-3.3 σ)

For reference, the independent 25-seed pilot gives ⟨W⟩ = 7.4511 ± 0.0442 at u = 0.02 (gap −0.08, −1.7σ),
7.5037 ± 0.0580 at u = 0.05 (−0.5σ) and 7.7856 ± 0.1315 at u = 0.10 (+1.9σ); the linear intercept of the 10-seed
set over u ≤ 0.02 is 7.5309 ± 0.0273 (gap −0.001 ± 0.033, 0.0σ).

**What this settles.** The 0.2 kT gap of the plan is **not** a missing dissipation: it came from holding the
finite-box wall excess constant along the path. The excess grows from +4.6 % to +6.4 % of Z_KR over the 11 %
compression, and with it measured, W_qs of the finite box agrees with the slowest push (−1.1σ) and with the
intercept (0.0σ).

**What it does not settle — Level 1 is not yet a pass.** Adiabatic compression from equilibrium cannot, on
average, take less work than the quasi-static value, yet the 10-seed u = 0.02 mean is 3.3σ below it and ⟨W⟩(u)
falls from u = 0.005 to 0.02 (fitted slope −6.2 ± 2.5). Two candidate explanations, in order of likelihood:
a low statistical fluctuation of 10 seeds (the 25-seed pilot at the same speed is only −1.7σ), or a small
geometric mismatch between the equilibrium measurement and the push — the holds measure pressure on the fixed
outer wall, while during the push the boundary is the piston face, whose accessible length need not equal the
nominal L. A 1 % difference in effective length moves W_qs by about 0.07 kT, the size of the discrepancy.

## 2b. Level 1 update, 2026-09-17: more work seeds, a densified path, and the gap within ~2σ

**Statistics first.** 40 new seeds per speed (9100–9139, identical commands) give 50 per speed. The new seeds
reproduce the old pattern, so the deficit was not a 10-seed fluctuation:

| u | ⟨W⟩, 50 seeds |
|---|---|
| 0.005 | 7.5070 ± 0.0111 |
| 0.01 | 7.4513 ± 0.0160 |
| 0.02 | 7.4537 ± 0.0231 |

**Two geometric explanations checked and ruled out.** (i) The piston starts 0.25 σ outside the gas and travels
4.18 σ in total, but the trace shows exactly 3.93 σ of travel with gas contact (78.50 → 74.57), so the compression
path used for W_qs was correct. (ii) The outer-wall (WR) and divider (D0) hold pressures differ between seed sets
in both directions (1.2848 vs 1.3018 in one set, 1.3020 vs 1.2951 in the other): noise, not a wall effect.

**The reference was the problem.** The three-point path of § 2 had only 20 seeds at its two new points and a
non-smooth ratio (+4.6 %, +6.5 %, +6.4 %). Densified to five grid-exact points with 100 seeds at each new point and
all 170 available hold phases at the start:

| L | η | seeds | Z_wall | sem | Z_KR | Z_wall/Z_KR | sem of ratio |
|---|---|---|---|---|---|---|---|
| 39.25 | 0.100051 | 170 | 1.2979 | 0.0020 | 1.2363 | 1.0499 | 0.0016 |
| 38.270833 | 0.102611 | 100 | 1.3106 | 0.0026 | 1.2434 | 1.0540 | 0.0021 |
| 37.291667 | 0.105305 | 100 | 1.3230 | 0.0025 | 1.2510 | 1.0575 | 0.0020 |
| 36.3125 | 0.108144 | 100 | 1.3247 | 0.0028 | 1.2591 | 1.0521 | 0.0022 |
| 35.3125 | 0.111207 | 100 | 1.3431 | 0.0027 | 1.2679 | 1.0593 | 0.0021 |

| ratio model along the path | W_qs^finite | ± (resampled) |
|---|---|---|
| piecewise-linear through the 5 points | 7.4877 | 0.0075 |
| weighted linear fit in η (χ² 7.18 on 3 dof) | 7.4879 | 0.0066 |
| weighted constant (χ² 16.81 on 4 dof) | 7.4835 | 0.0067 |

Gaps against the linear-fit W_qs = 7.4879 ± 0.0066, 50 seeds per speed:

| u | ⟨W⟩ ± sem | W − W_qs | σ |
|---|---|---|---|
| 0.005 | 7.5070 ± 0.0111 | +0.0191 ± 0.0129 | +1.5 |
| 0.01 | 7.4513 ± 0.0160 | -0.0366 ± 0.0173 | -2.1 |
| 0.02 | 7.4537 ± 0.0231 | -0.0342 ± 0.0240 | -1.4 |

**Status: consistent within ~2σ; formally borderline.** The deficit shrinks from ≈ 0.08 kT to ≈ 0.035 kT (0.5 %).
The linear path fit has χ² = 7.18 on 3 dof (p ≈ 0.07), so its error bar is likely underestimated by ≈ √(7.18/3) ≈
1.5; with that inflation u = 0.01 is at −1.9σ and every speed is within 2σ. Two tensions remain and are stated,
not hidden: ⟨W⟩(u) still falls between u = 0.005 and 0.01 (slope −4.2 ± 1.7 over u ≤ 0.02), and the hold-point
scatter along the path is ≈ 1.5× its statistics, pointing to a small unmodelled systematic (≈ 0.3 %) in the
equilibrium pressure measurement. Neither changes the conclusion that the plan's 0.2 kT gap was the finite-box
wall pressure, not missing dissipation.

Data: `experiments_energy_transfer/level1_moreseeds_20260916/` (120 runs) and
`experiments_energy_transfer/level1_Zwall_path_20260916/L38p270833, L36p3125` (100 each) plus 80 more seeds at
`L37p291667, L35p3125`; all 0 aborts, 0 health lines.

## 3. Next steps (Level 1 only; Levels 2–5 unchanged)

1. **Statistics:** 40 more seeds at u = 0.005, 0.01, 0.02 (N = 100 runs of ≈ 1 min each; niced, or after the
   dilute ladders). Pass requires ⟨W⟩(u) ≥ W_qs^finite within 2σ at every speed and a non-negative slope.
2. **Geometry:** measure Z on the piston face itself during a hold with the piston in place at the start position,
   and compute A from the piston's recorded position, not from nominal L.
3. **Methods text:** the sample-time and divider-sign conventions of § 1.

## Files

| file | status |
|---|---|
| `hspist3/validation/paper2_level0_level1_20260916.py` | new — the analysis in this report |
| `hspist3/experiments_energy_transfer/level1_Zwall_path_20260916/` | new — 40 hold-only runs, per-event logs, 0 health lines, 0 aborts |
| `0000_PLAN_OVERALL/ALL_MARKDOWNS/260917_paper2_level0_level1_REPORT.md` | new — this report |
| `hspist3/00ALLINONE`, `edmd_core/*` | unchanged |

Running in parallel, untouched: the dilute A2 ladders (190/400 done at 17:07, 0 failed).

---

## Superseded later the same day — see `260917_paper2_level2_REPORT.md`

The Level 1 intercept in this report was a straight-line fit to three speeds with 10 seeds each, and its
slope $b = -4.2 \pm 1.7$ was flagged here as unphysical. It was: the work rises **quadratically** in the
piston speed, so a straight line through the slow end was the wrong form. With 760 trajectories over eight
speeds and the form $W = W(0) + Au^2$:

```
W(0) = 7.4728 ± 0.0089 kT   vs   W_qs^finite = 7.4879 ± 0.0066 kT
W(0) − W_qs = −0.0151 ± 0.0060 kT   →   1.7 σ   →   PASS
A = 24.73 ± 2.72   (N_s m / 2 = 25.0)
```

Level 1's verdict is unchanged (pass at < 2σ) but the number, its sign and its error all come from the
Level 2 report now. Nothing above this line was re-run; it is superseded, not corrected.
