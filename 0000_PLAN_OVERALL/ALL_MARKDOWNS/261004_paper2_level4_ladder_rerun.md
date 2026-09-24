# Level 4 ladder rerun at 65 τ_true — pre-registration, η settled, A2 re-derived

2026-10-04. **This header was written and put on disk while the runs were executing**, before any
rerun record was fitted. Sections 3–4 were complete before the runs finished; section 5 is the
result.

---

## 1. PRE-REGISTERED ANALYSIS (fixed before the data)

Identical rules to 261003 — nothing about the estimator, the model or the fit window changes. Only
the record length does.

* ω fixed from cot K = αK, α = M/(2N_s m), with Kolafa–Rottner c_s **at the corrected η = 0.101342**;
* three estimators: modelled two-component ACF fit, block, and S(0) with the **amplitude free**
  (the corner-frequency form — the σ²-fixed version is the one that was buggy);
* bias calibration per mass on synthetic OU + AR(2), grid **centred on the MEASURED τ and spanning
  0.5–2× it**, so it brackets the answer by construction;
* adoption rule: modelled estimator adopted only if its calibration slope **≥ 0.6**;
* self-consistent ACF fit window at 8 τ_measured, iterated;
* report per mass: L/τ_actual, calibration slope, all three τ, their spread;
* exponent from all five masses, from {10, 20, 50}, from {50, 100, 200}, each with all three
  estimators.

**VERDICT RULE, fixed now:** *b is a result only if every mass has L/τ ≥ 60 AND calibration slope
≥ 0.6 AND the light and heavy subsets agree within 2σ. Otherwise report the failing diagnostic and
stop.*

## 2. The runs

M_d = 50, 100, 200, 80 seeds each; M = 10 and M = 20 are **not** repeated (L/τ 64 and 61,
slope 0.683 and 0.630 — both in spec). Record = 65 × **max**(modelled, block, S(0)) from the 261003
table, so even the pessimistic τ gets 65 lengths.

| M | τ_max used | record [σ] | steps/seed | mode period | N | samples/period | rows | MB/seed | core-h |
|---|---|---|---|---|---|---|---|---|---|
| 50 | 4155 | 270 075 | 16 212 446 | 125.9 | 350 | 21.6 | 46 321 | 2.4 | 1.8 |
| 100 | 12 935 | 840 775 | 50 471 238 | 157.6 | 450 | 21.0 | 112 158 | 5.9 | 5.6 |
| 200 | 53 404 | 3 471 260 | 208 377 735 | 207.6 | 600 | 20.8 | 347 296 | 18.1 | 23.2 |

**30.6 core-hours, ≈ 3.4 h wall clock at 9 concurrent; 2.1 GB reduced.**

**Two constraints collided and the tie was broken deliberately.** "≥ 20 samples per mode period"
and "≤ 5 MB per seed" are incompatible at M ≥ 100: M = 200 needs N ≥ 2176 for the size cap and
N ≤ 622 for the resolution floor. **The resolution floor is kept** — resolving the oscillation is
what the two-component fit depends on, and 2.1 GB against 27 GB free is not worth trading for it.
M = 200 is 18.1 MB per seed.

## 3. η — read from the recorded parameters, not inferred

My earlier inference (t = 0.05 from 39.25 − 38.75 = 0.5 = 2r) was **void**: 0.5 is equally t/2 for a
1.0-thick divider, so the arithmetic identifies nothing. From the runs themselves:

> `wall_thickness_sigma = 1` in the recorded `summary.csv` of **every** Level 4 cell — the Md10
> long run, all four 261003 ladder masses, the 260928 thermal set and the 260926 v3 set.

So the divider is **1.0 σ** thick, the gas occupies 0–38.75 and 39.75–78.5, and

> **η_phys = 100πr²/((78.5 − 1.0)·10) = 0.101342.** The plan's long-standing 0.1013 was right;
> my 0.100114 was wrong.

Paper 1's t = 0.05 (the `WALL_T` guard) is a different geometry and is unaffected.

| quantity | at 0.100051 (used in 261003) | at **0.101342** | change |
|---|---|---|---|
| Z | 1.23628 | 1.23988 | +0.291 % |
| c_s | 1.74434 | **1.74930** | +0.285 % |
| 1 + ηZ′/Z | 1.22489 | 1.22815 | +0.266 % |
| τ_T prefactor | 50.483 M | **50.349 M** | −0.265 % |
| mode period, α = 0.1 | 95.16 | **94.89** | −0.284 % |

**The M = 10 verdict is unchanged**, because the plan's 503.5 was already computed at 0.1013:
491 ± 29 is **−0.4σ** from the hard-disk isobar and **−4.4σ** from the ideal gas.

## 4. Paper 1's A2 ladder at the physical η — the surface scaling survives, but shrinks

Paper 1's geometry has t = 0.05 and 2L₀ ∝ √N at fixed η, so the artefact in the KR deviation is
−S(η)·t/(2L₀) with S = d ln c_KR/d ln η, i.e. **∝ 1/√N** — the same form as the surface-scaling
hypothesis it contaminates.

**The prediction is confirmed cell by cell**: at η = 0.65 the measured change is −1.291 % (N = 100)
and −0.312 % (N = 1600) against a predicted −1.269 % and −0.317 %.

| η | N | dev OLD % | dev NEW % | change | artefact predicted |
|---|---|---|---|---|---|
| 0.500 | 100 | +1.709 | +1.098 | −0.611 | −0.600 |
| 0.500 | 1600 | +0.412 | +0.261 | −0.151 | −0.150 |
| 0.600 | 100 | +1.363 | +0.288 | −1.075 | −1.060 |
| 0.600 | 1600 | +0.076 | −0.189 | −0.265 | −0.265 |
| 0.650 | 100 | +2.653 | +1.361 | −1.291 | −1.269 |
| 0.650 | 1600 | −1.324 | −1.636 | −0.312 | −0.317 |

Re-fitting dev = a + b/√N at η ≥ 0.5:

| η | OLD a, b | NEW a, b | pure artefact b | fraction of the slope that was the η axis |
|---|---|---|---|---|
| 0.500 | +0.61, **+11.1** | +0.61, **+5.0** | −6.0 | **54 %** |
| 0.600 | −0.34, **+15.5** | −0.34, **+4.7** | −10.6 | **68 %** |
| 0.650 | −2.17, **+51.8** | −2.16, **+38.7** | −12.7 | **25 %** |

**The N → ∞ intercept is unchanged to 0.01 % at every density** — as it must be, since the artefact
vanishes as N → ∞. So the extrapolated deviations, and any claim resting on them, stand.
**What changes is the slope**: at η = 0.5 and 0.6 **more than half** of the 1/√N "surface scaling"
was the η axis, not surface physics. At η = 0.65 three quarters of it survives.

Paper 1's low-density results, the estimator, the theory and the Román re-mapping are untouched:
the shift is +0.013 % in η at η = 0.02 and −0.001 % in the deviation.

## 5. Ladder result — the verdict rule fails, but the explanation has flipped

**240/240 runs, 0 aborts, 0 health events**, 15:54–19:55 with one clean pause/resume.

### 5.1 The rerun did what it was for

| M | L/τ before → after | cal slope before → after | estimator spread | τ before → after |
|---|---|---|---|---|
| 50 | 40 → **70** | 0.601 → **0.746** | → **1.16×** | 4060 → 3835 |
| 100 | 29 → **68** | 0.574 → **0.741** | → **1.09×** | 11 399 → 12 330 |
| 200 | **17 → 87** | **0.476 → 0.734** | **1.86× → 1.05×** | 37 400 → 40 079 |

Full table, all five masses:

| M | seeds | record | L/τ actual | cal slope | modelled | block | S(0) | spread |
|---|---|---|---|---|---|---|---|---|
| 10 | 80 | 30 735 | 64 | 0.708 | 482 ± 23 | 489 | 448 ± 14 | 1.09× |
| 20 | 80 | 63 456 | **58** | 0.667 | 1091 ± 53 | 1135 | 1003 ± 22 | 1.13× |
| 50 | 80 | 268 079 | 70 | 0.746 | 3835 ± 95 | 4384 | 3792 ± 28 | 1.16× |
| 100 | 80 | 838 786 | 68 | 0.741 | 12 330 ± 855 | 13 091 | 11 989 ± 89 | 1.09× |
| 200 | 80 | 3 469 269 | **87** | 0.734 | 40 079 ± 1708 | 41 460 | 39 596 ± 109 | **1.05×** |

**M = 200 is now the best-sampled point on the ladder** — L/τ = 87, slope 0.734, three estimators
agreeing to 1.05×.

### 5.2 The verdict rule, applied as written

| criterion | result |
|---|---|
| every mass L/τ ≥ 60 | **FAIL** — only at M = 20 (58). The other four are 64, 70, 68, 87. |
| every mass slope ≥ 0.6 | **PASS** — 0.708, 0.667, 0.746, 0.741, 0.734 |
| light and heavy subsets agree within 2σ | **FAIL** — 8.2σ (modelled), 5.2σ (block), 6.9σ (S(0)) |

> **b IS NOT QUOTABLE. Reporting the failing diagnostic and stopping, as the rule requires.**

### 5.3 But last night's explanation is now excluded

The 261003 verdict was "the heavy end is an estimator limit". **That is no longer tenable.** The
heavy end was fixed — L/τ 17 → 87, slope 0.476 → 0.734, spread 1.86× → 1.05× — and

> **b_heavy went UP, not down: 1.602 → 1.692.**

An estimator artefact removed by better sampling does not strengthen. The subset disagreement
survived the treatment designed to cure it, and is now larger and measured with three estimators
that agree with each other at the heavy end to 5 %.

### 5.4 What the data actually are: not a power law

| M range | τ ratio | M ratio | local b |
|---|---|---|---|
| 10 → 20 | 2.26 | 2.0 | **1.18** |
| 20 → 50 | 3.52 | 2.5 | **1.37** |
| 50 → 100 | 3.21 | 2.0 | **1.68** |
| 100 → 200 | 3.25 | 2.0 | **1.70** |

τ/M rises 48.2 → 54.5 → 76.7 → 123.3 → 200.4. **The local logarithmic slope climbs from ≈ 1.2 at
the light end and plateaus at ≈ 1.7 for M ≥ 50.** A single b fitted across the range (1.49) is
therefore a fit to a curved function and is not a meaningful exponent — which is the deeper reason
the subset test fails, and why it should keep failing however well each point is measured.

### 5.5 The one remaining systematic, and how to kill it cheaply

M = 10 and M = 20 were **not** rerun, and they are now the *worst*-sampled masses: L/τ 64 and 58,
calibration slopes 0.708 and 0.667 against 0.73–0.75 at the heavy end. A lower slope means more
bias to correct, so if the calibration under-corrects at the light end, τ_light is pulled down and
manufactures exactly this curvature.

> **Topping up M = 10 and M = 20 to L/τ ≥ 65 costs 0.68 core-hours — about 10 minutes at 9
> concurrent.** M = 10 needs 1.02× its record, M = 20 needs 1.12×.

That closes the L/τ criterion on every mass *and* tests whether the light end is driving the
curvature. **Not launched** — the rule said stop, and this is the next decision, not tonight's.

### 5.6 Status of the claim

**Unchanged and safe:** τ_T at M_d = 10 matches the parameter-free hard-disk theory and excludes
the ideal gas; the mode period tracks cot K = αK across α = 0.1–2.0.

**Open, and sharper than last night:** τ_T(M) is superlinear and *curved*, the local slope reaching
≈ 1.7 by M = 100–200, and this is **no longer attributable to under-sampling at the heavy end**.
Whether it survives the light-end top-up is the next test. If it does, the quantity to report is
not an exponent but the curve τ_T(M)/M, and the adiabatic-piston theory's linear law is wrong in
this box in a way that grows with mass.
