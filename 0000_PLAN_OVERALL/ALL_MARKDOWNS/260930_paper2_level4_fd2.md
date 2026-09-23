# Level 4 fluctuation route, pass 2 — the oscillation is Paper 1's divider mode

2026-09-30. Post-processing only on `level4_equilibrium_20260929/Md10` (20 seeds, 0 aborts, 0
health), plus one gated binary change and one document repair. No new runs.

**The period-95.3 mode is not unexplained. It is Paper 1's divider resonance, and the OPEN flag is
closed at 0.03 %.** The previous report applied the heavy-divider limit of `cot K = αK` to a light
divider — α = 0.1, not α ≫ 1 — and got 44.7 instead of 95.2.

---

## 1. The mode, identified

Román's divider eigenmode satisfies `cot K = αK` with α = M/(2N_s m). Here

> α = 10/(2 × 50 × 1) = **0.1** — a *light* divider.

The gas is not a massless spring; it is a standing wave carrying most of the inertia. Solving the
first branch gives **K = 1.428870**, close to the α → 0 limit K → π/2 = 1.5708. With
Kolafa–Rottner at η = 0.100051 (Z = 1.236285, ηZ′ = 0.278027, c_s² = 3.042711, **c_s = 1.744337**)
and L_eff = 38.75 − 2r = 37.75:

> ν = c_s K /(2π L_eff) = 0.010508 → **T = 95.16 σ-time** [DERIVATION]

| | period [σ-time] |
|---|---|
| **predicted, KR c_s** | **95.16** |
| measured, spectral peak of T₁−T₂ | **95.13 ± 5.15** |
| measured, spectral peak of divider x | **95.22 ± 7.04** |
| measured, ACF fit with ω free (T₁−T₂ / x) | 93.27 ± 0.51 / 94.31 ± 0.57 |
| predicted, ideal gas (c_s = √2) | 117.38 |
| predicted, Paper 1's measured c_s (+1.01 %) | 94.21 |

**0.03 % and 0.06 % on the spectral peaks**, whose errors are the peak FWHM. The ideal-gas variant
is 4.3σ away, so the mode discriminates KR from the ideal gas on its own. The ACF-fit periods land
1–2 % low, between KR and the Paper-1 c_s variant; the three period estimators are mutually
consistent only at the ~2 % level, so no finer claim than "KR, not ideal" is made.

**Why the previous formula was wrong.** ω = √(2k_gas/M) is the α ≫ 1 limit, where cot K ≈ 1/K gives
K ≈ α^(−1/2). At α = 0.1 that would be K = 3.162, which is not even on the branch whose root is
1.4289 — the asymptotic form is not merely inaccurate here, it is inapplicable. Its prediction was
44.71. The "effective mass 45 rather than 10" of the previous report is exactly the standing-wave
inertia Román's solution accounts for, so that observation was the symptom, not a puzzle.

**Two papers, one equation.** Paper 1's mass ladder and Paper 2's Level 4 autocorrelation are the
same eigenvalue problem at opposite ends of α: Paper 1 runs α = 0.5–20, this run α = 0.1.

## 2. The ACF has a known structure, so it can be modelled

With ω fixed by item 1, fit C(t) = A e^(−t/τ_T) + B e^(−t/τ_r) cos(ωt) to the full-resolution
seed-averaged autocorrelation (lags 0–1500 σ-time, Δt = 5.0), jackknifed over seeds:

| observable | A | τ_T | B | τ_r | fit rms |
|---|---|---|---|---|---|
| T₁ − T₂ | 0.472 ± 0.024 | 311 ± 34 | 0.503 ± 0.021 | 195 ± 17 | 0.034 |
| divider x | 0.468 ± 0.023 | 316 ± 38 | 0.511 ± 0.021 | 186 ± 15 | 0.026 |

The two components split the variance almost exactly in half, and the two observables agree on
every parameter.

**Bias calibration of this estimator**, by the same method as before — synthetic series with the
same record, seed count, blocking-free fit and an AR(2) oscillatory component matched to
(A, B, τ_r):

| true τ_T | 250 | 350 | 450 | 504 | 550 | 618 | 700 | 900 |
|---|---|---|---|---|---|---|---|---|
| recovered | 196 | 242 | 286 | **329** | 336 | **377** | 398 | 478 |

Response slope **0.438**, against 0.335 for the block estimator. Inverting, and taking the larger of
the jackknife and the Monte-Carlo spread as the error:

> **τ_T = 482 ± 100** (T₁−T₂) and **488 ± 107** (divider x)
> vs hard disk 503: **−0.5σ, −0.3σ**  ·  vs ideal 618: **−1.9σ, −1.6σ**

**Decision, by the pre-registered rule: the slope is 0.438, below the 0.6 threshold, so the block
estimator remains the headline.** The modelled fit is reported as a corroborating measurement, and
it corroborates well — 482 and 488 against the block estimator's 490 ± 190 and 450 ± 165, from the
same data by an estimator with a different bias. Two independent analyses of one record agreeing at
this level is the strongest statement this run supports. The headline is unchanged:

> **τ_T(M_d = 10) is of order 500 σ-time, consistent with both hypotheses; the 23 % gap is not
> resolved.**

**But the modelled estimator is what KOA should use**, because its slope clears 0.6 once the record
is long enough — and at the recommended size it beats the block estimator outright:

| record [σ] | L/τ | seeds | slope | separation |
|---|---|---|---|---|
| 8 200 (this run) | 16 | 20 | 0.45 | 1.4σ |
| 16 400 | 33 | 20 | 0.56 | 1.6σ |
| 32 800 | 65 | 20 | 0.74 | 1.9σ |
| **32 800** | **65** | **80** | **0.74** | **4.0σ** (block estimator: 3.3σ) |

## 3. τ_r is Paper 1's linewidth at a mass Paper 1 never ran

For a damped-cosine ACF the power spectrum is a Lorentzian of fractional width
Δf/f = 1/(π τ_r ν). Measured here: **Δf/f = 0.155 and 0.163, i.e. Q = 6.4 and 6.1.**

Paper 1's linewidth analysis uses Mansour's piston form, Δf/f = Γ L_y √(2/(M̂N)) with
M̂ = M + mN/3. Extrapolating to M = 10 — the *only* free step, since L_y = 10, N = 100 and η = 0.10
are shared — gives

| | Δf/f | τ_r |
|---|---|---|
| Enskog Γ = 0.331 | 0.0711 | **426** |
| Paper 1's own mean Γ = 0.314 ± 0.038 | 0.0675 | 449 |
| **measured here** | **0.155 / 0.163** | **195 ± 17 / 186 ± 15** |

**The resonance is 2.2× broader than the extrapolation**, i.e. Γ_implied = 0.72 and 0.76 against
0.331. Stated as a measurement, not a discrepancy to be explained away: M_d = 10 is a factor 5
below Paper 1's lightest mass and in the opposite α regime, where the mode is a gas standing wave
rather than a loaded piston, and the Mansour form is derived for a heavy piston. The extrapolation's
own footing is soft — Paper 1's inferred Γ at η = 0.10 scatters from 0.084 to 0.439 across its nine
masses. This is a new number at a new α, and it is the one place where this run's damping and
Paper 1's disagree.

## 4. Trace decimation — built and gated

The energy-transfer trace was the only writer in `00ALLINONE.c` with no cadence control: it ignores
`--output-dt` and emits one row per step. That is 21 GB for one 10 000 σ-time cell, and the
656 000 σ-time record §2 recommends at M_d = 200 would be terabytes. Without this, no KOA-length
record is possible.

**`--trace-every=N`** — decimation in *steps*, not σ-time, so it is exact and seed-independent.
`N ≤ 1` is the historical behaviour. `--trace-every` is collision-free: no other flag begins with
`--trace`, which matters because every branch in the parser uses `strncmp` prefix matching with no
terminator check (the trap that hides `--seeding` behind `--seed` to this day).

The guard wraps the row write only. **`recorded_steps++` is deliberately left outside it**: it is
also the loop's termination counter, so decimating it would shorten the run.

> **GATE PASS.** `binary_gate.sh ./00ALLINONE_base ./00ALLINONE`: energy_transfer 6/6 and
> speed_of_sound 9/9 data files **byte-identical across 3 seeds**, every summary.csv column
> identical. **No `EXPECT_DIFF` was needed** — the default path is untouched, so nothing differs at
> all. Both binaries built `make release` (`-O3 -march=native`); note `make` alone builds the ASan
> debug target and must not be used for science.

Functional test, which the gate deliberately does not do (the gate proves the *default* is
unchanged):

| | rows | bytes |
|---|---|---|
| full | 12 000 | 4 413 881 |
| `--trace-every=300` | 40 | 15 036 |

and **every kept row is byte-identical to the full run's every-300th row** — decimation selects
rows, it does not recompute them.

## 5. Plan hygiene

75 lines were duplicated verbatim (1-indexed 852–926). The duplicate contained not only
"What replaces the impedance criterion" but a **second `\section{Level 4}` heading**, so the
document had Level 4 twice. Deleted. The outline is now
Level 3 → l3ceiling → l3acoustic → Level 4 → `sec:l4tau` → `sec:l4tau-superseded` → Level 5, and
`sec:l3acoustic`, `eq:onedof`, `eq:transmission`, `eq:twobody` are each defined once. The
`-superseded` box is untouched.

## 6. Status of Level 4

| | |
|---|---|
| (i) work splits in half | **holds**, 0.1–1.5σ |
| (iii) divider takes Δx/2 | **holds**, 0.1–1.3σ |
| stage-1 offset, KR vs ideal | **KR favoured 7–8σ** |
| σ(T₁−T₂) = 4T̄/(2√(2N+1)) | **confirmed**, 0.1980 vs 0.199 |
| **the equilibrium oscillation** | **Paper 1's divider mode, 95.13 ± 5.15 vs 95.16 predicted** — closed |
| τ_T at M_d = 10 | **~490**, consistent with 503 and with 618 |
| τ_T ∝ M | **untested** — one mass only |
| τ_r vs Paper 1's linewidth | **measured, 2.2× broader than the extrapolation** — open |
| Eq. (51), Eq. (52) | **untested** |

Route B (stationary heat conduction, per-wall baths with separate ledgers) stays a design in the
plan: it is a core change to the energy accounting and needs the Level 0 gate re-run. After KOA
onboarding, not before.
