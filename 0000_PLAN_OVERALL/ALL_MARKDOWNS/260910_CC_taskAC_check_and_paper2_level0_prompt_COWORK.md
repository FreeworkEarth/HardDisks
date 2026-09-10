# 260910 — Check of CC's TASK A/B/C, the refit figure, the go for N = 2500, and the Paper 2 Level 0 prompt (COWORK)

Written 2026-09-10 from `260909_wall_overdue_and_temperature_resolution_CC.md`, `260909_cs_refit_notes_CC.md`, `260909_plots/routeA_refit_cs_vs_eta_20260909.csv`, `routeA_combined_speed_of_sound_summary.csv`, `validation/refit_sound_speed_manifest.py`, the regenerated `260908_pressure_final_analysis.md`, and CC's pasted session output. Figure regenerated here: `260909_plots/260909_cs_vs_eta_refit.png/.pdf`.

---

## 1. TASK A — verdict (a) accepted, with two corrections to the note

The evidence is what I asked for and it is conclusive: all 675 dilute trajectories reproduced from seed with a print-only trace; 675/675 ledger counts reproduced; zero counter changes during hold or measurement; every one of the 662 events at t = 0, right wall, gap exactly 0, on the two particles the float seeder puts on the face; the count equals the number of those two whose initial v_x points inward. The code path (`wall_time_from_gap` → rc = 2 → `schedule_walls` counts and schedules a bounce at t_col = 0) resolves as an immediate elastic reflection — physically the same as starting that particle with the opposite sign of v_x. The 504 trajectories are good and go back into the fits. This also removes the cloud over the ideal-gas end of the c_s data.

Two things to fix in `260909_wall_overdue_and_temperature_resolution_CC.md` before it becomes provenance:

1. **"The 2000-step hold (≈ 800 σ-time at dt = 0.4/24)" is off by the factor 24.** 2000 × 0.4 = 800 px-time; in σ-time that is 800/24 = **33.3**. The leaf log confirms the unit: `samples=159840 T=2664 sigma-time` gives dt = 0.01667 σ-time = 0.4/24. The argument survives (one reflection at t = 0 is irrelevant either way), but 33 σ-time at η = 0.02 is only ≈ 3 collisions per particle — which is fine only because the initial velocities are already Maxwellian and rescaled to kT = 1; say that instead of "800".
2. **The count distributions at L0 = 200 and L0 = 100 are identical (63/115/47) because the per-run seeds are identical across those η campaigns** (same `--seed` base, same (M, run) index). Harmless for c_s(η), but it means the velocity draws are not independent across densities, and the "expected ¼/½/¼" comparison is one sample, not three. Note it.

And one action item: the float seed pad (`seed_pad_px = fmaxf(1e-4f, 1e-5f·d)` at ≈ 9800 px) must be fixed in the driver before any new c_s run at L0 ≥ 100 — carry the pad in double or inset by a fixed 1e-3 d as the core lattice seeder does. Not for the existing data (nothing changes for them), and not part of any campaign that is running.

Temperature (A.2): accepted. kT = 1 exactly per compartment at t = 0 by construction (hard rescale at lines 5612–5618, then per-segment equalization), conserved through the fixed-divider hold. No durable record, exact by construction, reproducible from seed. Put the K_B vs `kB_effective()` wrinkle in a footnote and move on.

## 2. TASK B — accepted

The χ² column reproduces my numbers exactly (0.67: 7.62, p = 0.006; 0.69: 2.29; 0.65: 0.07). η = 0.30 at χ²₁ = 4.05 (p = 0.044) is one flag among twelve tests at the 4.6 % level — expected about once; CC reported it and did not tune it, which is right. The block-identity table (50/50, 110/110 identical for every seed) is in the document and "16×" is gone. Boundary table replaces the "threshold". Methods paragraph rewritten with the calendar-pop accounting. Calibration stderr now logged. Nothing left here.

## 3. TASK C — the refit is right; what it does and does not say

**Accepted:** 7200 manifest rows, all joins matched, 7186 eligible (14 boundary peaks dropped at 0.55 and 0.76), the 504 flags retained as the t = 0 class, exact cot K = αK, L_eff = L0 − 2r, per-mass means weighted by 1/sem², through-origin and free-intercept both reported. The old summary was a free-intercept fit; CC's free column reproduces it to < 0.7 % except where boundary peaks were dropped. Good work, and the honest part is the intercept column: where origin and free differ by more than their errors, the ν(K) points do not extrapolate to zero and the single-mode relation is being stretched.

**Framing that has to go into the paper, or the figure misleads:** route A is not a density sweep at fixed geometry. N = 100 is fixed, so L0 = 3.93 σ/η: the compartment is 200 σ long at η = 0.02 and 6 σ at η = 0.65. Confinement changes along the axis together with density. The refit vs Kolafa–Rottner (through origin, fit errors):

| range | deviation | significance |
|---|---|---|
| η = 0.020–0.079 (L0 = 200–50) | −0.49, −0.28, −0.21, −0.02, −0.14 % | 4.9 σ at the first point, 1.7–1.9 σ at the others |
| η = 0.11–0.39 (L0 = 35–10) | +0.52 … +1.11 %, monotone | 4.5–12.7 σ |
| η = 0.52–0.55 (L0 = 7.5–7.1) | +2.5, +2.0 % | 8–10 σ |
| η = 0.57–0.63 (L0 = 6.9–6.2) | +0.7, −0.02, −1.4, −0.7 % | but intercepts of −0.001 … −0.003 absorb up to 3.6 % — the model is off here |
| η ≥ 0.65 | +5 % and up | not a fluid-branch measurement (structure, 6 σ box) |

So the N = 100 sweep validates the sound-speed *instrument* at the ≤ 1 % level for η ≤ 0.4 and shows the sign-changing finite-size/confinement systematic beyond that. The *bulk* c_s claim for Paper 1 cannot come from this sweep; it has to come from the fixed-aspect family (N = 100/400/900/1600 at fixed geometry, `finitesize_aspect_20260826`) with the same manifest treatment and a 1/√N (or 1/L) extrapolation with χ² — GPT's handoff item 4, which CC correctly left undone. That is the next c_s task (section 5).

**The dilute −0.5 %.** The first point (η = 0.0196, L0 = 200) sits 4.9 σ below the EOS; the next four are 0.1–0.3 % low. The slope of c_s at η → 0 is right (c_s/√2 − 1 = 0.035 vs 2η = 0.039), so the ideal-gas limit itself is fine; something shaves a few tenths of a percent off the *measured frequency* in the longest compartments. Three candidates, separable without new runs from the per-mass residuals already in the manifest: a damped-frequency shift (ν_d = ν₀√(1 − ζ²), strongest for the lightest divider, i.e. residual grows toward small M); the L_eff convention (L0 − 2r vs L0 − r: 0.25 % at L0 = 200, 0.5 % at 100 — M-independent, but it would grow with η, which the data do not show); finite N per side (M-dependent through α = M/2N). A per-mass residual table at the five dilute η decides it — one script run.

**Figure.** `260909_cs_vs_eta_refit.png/.pdf` regenerated with the refit table: all 32 densities with error bars, the legend now states L0 = 3.93 σ/η, the ≥ 0.65 points are labelled as not fluid-branch, y-axis log to 40 so nothing is clipped. Inset: the ideal-gas end with the five dilute points on the KR curve and the exact first-order law. This is the Paper 1 sound-speed figure candidate, minus the fixed-aspect extrapolation panel that still has to be built.

## 4. Go for N = 2500

The machine is free of c_s runs. Launch TASK 7 now (prompt-1 review section 7 item 8): N = 2500 at η = 0.65, 0.67, 0.69, 3 seeds, equilibrated calibration with stderr logged, chunk = min(0.6·min(A,B), 0.75·320/N = 0.096), same 400 + 30 × 20 protocol and health contract. It runs for hours; everything below is analysis-only and runs beside it.

## 5. Prompt for Claude Code (paste as is)

```
Read 0000_PLAN_OVERALL/ALL_MARKDOWNS/260910_CC_taskAC_check_and_paper2_level0_prompt_COWORK.md
first. Rules unchanged: no core physics changes, no edits to accepted trajectories, no commits,
quote numbers from files. Order matters: launch step 1, then do 2-5 while it runs.

1. N = 2500 pressure run (go given): eta = 0.65, 0.67, 0.69, 3 seeds (seed rule unchanged),
   equilibrated calibration with seeds 911000001/2 and stderr to calib_*.err, production chunk
   = min(0.6*min(A,B), 0.096), 400 + 30x20 protocol, strict health contract, resumable driver.
   Report the calibration rows and the launch state, then leave it alone. When it lands: rerun
   analyze_pressure_campaign.py (4 sizes, chi2 with 2 dof), regenerate the document, show me
   table B before touching the claim text.

2. Fix the two items in 260909_wall_overdue_and_temperature_resolution_CC.md: hold = 2000 x
   0.4/24 = 33.3 sigma-time (not 800); note that the per-run seeds are identical across the
   L0 = 200/150/100 campaigns (identical 63/115/47 counts). Add the float seed-pad fix as an
   open driver item (double-precision pad or fixed 1e-3 d inset), NOT implemented now.

3. Dilute residuals: from routeA_fit_input_manifest_20260909.csv, for eta = 0.0196, 0.0262,
   0.0393, 0.0524, 0.0785, table of (M, mean nu, sem, nu_predicted_from_EOS = c_s_KR * K/(2 pi
   L_eff), residual %) per divider mass. State whether the residual grows toward small M
   (damping), is M-independent (L_eff), or something else. No new runs.

4. Fixed-aspect family A (finitesize_aspect_20260826, famA, N = 100/400/900/1600, eta = 0.4 ...
   0.72 as available): apply the same manifest treatment (eligibility, health classes, exact
   cot K, through-origin and free intercept), then per eta a fit c_s(N) = c_s_inf + a/sqrt(N)
   with chi2 (2 dof), the 900/1600-only intercept, and the comparison to KR (<= 0.69) and Liu.
   This is GPT's handoff item 4 and the actual bulk c_s claim. Report the table; do not write
   claim text yet.

5. PAPER 2, LEVEL 0 (energy ledger) -- analysis only, no new campaign.
   a. Inventory the existing energy-transfer reference runs (directories, 00_COMMAND, N, piston
      speed, work target, wall types, divider mass, spring constant). For each, quote the header
      of every logged per-sample file. Say which of these are logged and which are not: cumulative
      piston work, gas KE (total and per compartment), divider KE, spring energy, heat to thermal
      walls, piston position/velocity, per-collision impulses. "Not logged" is an answer.
   b. Ledger per run and per sample: W_piston(t) - [dKE_gas + dKE_divider + dE_spring + Q_walls](t),
      with the sign convention written down once. Report max |residual| in units of kT and
      relative to W. For elastic walls the residual must be roundoff; state the number you get.
   c. Piston collision rule: quote it from the code; verify v' = 2u - v and dE = 2 m u (u - v) per
      collision (u = piston velocity, v = particle normal velocity before the hit; fix the signs
      from the code). If per-collision data are logged, sum them over one run and compare with
      the logged W. If they are not, specify -- do not launch -- the smallest run that would log
      them (N = 50, one push, per-event log).
   d. Momentum: change of total gas + divider momentum vs the impulse delivered by pistons and
      outer walls, per run.
   e. Time reversal: one short reference run (existing binary, N <= 100, elastic walls), reverse
      all velocities at t1, and report the return error at t = 0 for t1 = 1, 3, 10, 30 sigma-time.
      This bounds what "exact" means for the ledger. If the driver cannot reverse velocities
      without a code change, say so and stop there.
   f. "SpringE_max/W_in": quote its definition and window from the code, and state in two
      sentences why it is not a thermodynamic work.
   Deliver 260910_paper2_level0_energy_ledger_CC.md with the tables. No plots yet.

Report commands and outputs verbatim; git status --short at the end.
```

Sequence for me afterwards: I check step 3 and 4 numbers against the manifest and the famA files, and the Level 0 ledger against the logs, and then we write the Paper 1 claim text for c_s alongside the pressure one. Paper 2 Level 1 (divider fluctuations in a proper equilibrium run) comes after the ledger closes, not before.
