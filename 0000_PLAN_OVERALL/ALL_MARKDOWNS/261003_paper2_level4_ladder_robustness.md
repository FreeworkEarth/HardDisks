# Level 4 ladder — robustness checks. The exponent does NOT survive them.

2026-10-03. Post-processing only on `level4_ladder_20261003` (M = 20, 50, 100, 200; 320 runs,
0 aborts, 0 health) plus `level4_equilibrium_KOAlength_20261002` (M = 10).

> **b = 1.43 is withdrawn: the heavy end is an estimator limit.** That is not the same as saying
> the exponent is 1 — the in-spec masses still give b = 1.32 ± 0.11, ≈ 2.9σ above the linear law.
> **The exponent is OPEN and the rerun at 65 τ_true decides it.** The criterion set before the checks — *"a real exponent is the same in
> both subsets; an estimator effect makes it grow with the heavy subset"* — fires unambiguously.
> **What survives at M = 10 is unchanged**: τ_T = 491 ± 29, hard-disk isobar at −0.6σ, ideal gas
> excluded at 3.5σ.

---

## 0. Two corrections to the inputs, before any check

**η — I was wrong, and the recorded parameter settles it.** I inferred t = 0.05 for the Level 4 box
from 39.25 − 38.75 = 0.5 = 2r. That inference is void: 0.5 is *also* t/2 for a 1.0-thick divider, so
the arithmetic distinguishes nothing. Read from the runs instead, `wall_thickness_sigma = 1`
in the recorded `summary.csv` of **every** Level 4 cell (Md10 long, all four ladder masses, the
260928 thermal set and the v3 set). So the divider is **1.0 σ thick**, the gas occupies 0–38.75 and
39.75–78.5, and

> **η_phys = 100πr²/((78.5 − 1.0)·10) = 0.101342** — the plan's long-standing 0.1013 was right.

Paper 1's t = 0.05 (the `WALL_T` guard) is a different geometry and is unaffected. Corrected
constants: Z = 1.23988, **c_s = 1.74930** (+0.285 %), 1 + ηZ′/Z = 1.22815, **τ_T = 50.349 M**
(−0.265 %), mode period at α = 0.1 **94.89** (−0.284 %). The M = 10 verdict is unchanged, because
the plan's 503.5 was already computed at 0.1013: measured 491 ± 29 is **−0.4σ** from the hard-disk
isobar and **−4.4σ** from the ideal gas.

**The S(0) estimator had a bug on its first run.** Fixing σ² to the *total* variance (0.199²)
is wrong: the OU component carries only A ≈ 0.47–0.56 of it, so τ came out ~2.6× low even at
M = 10 where everything else agreed. Refitting a Lorentzian with the **amplitude free**, so τ comes
from the corner frequency alone and no variance is assumed, removes the assumption entirely. Only
the corrected version is used below.

## 1–2. Three estimators per mass

| M | L/τ actual | cal slope | modelled | block | S(0), amplitude free | τ_r | B/(A+B) | spread |
|---|---|---|---|---|---|---|---|---|
| 10 | 64 | 0.683 | 479 ± 18 | 522 ± 36 | 448 ± 14 | 207 | 0.467 | **1.17×** |
| 20 | 61 | 0.630 | 1044 ± 59 | 1079 ± 49 | 1003 ± 22 | 284 | 0.481 | **1.08×** |
| 50 | 40 | 0.601 | 4060 ± 209 | 4155 ± 114 | 3546 ± 43 | 497 | 0.506 | **1.17×** |
| 100 | 29 | **0.574** | 11399 ± 837 | 12935 ± 706 | 10650 ± 186 | 775 | 0.534 | **1.21×** |
| 200 | 17 | **0.476** | 37400 ± 2685 | 53404 ± 3357 | 28732 ± 329 | 1115 | 0.559 | **1.86×** |

Three things to read off. The **calibration slope falls monotonically** and drops **below the
pre-registered 0.6 adoption threshold at M = 100 and M = 200** — by the stated rule the modelled
estimator is *not adopted* at those two masses. **L/τ actual collapses** from 64 to 17, because the
records were sized at 65 τ_pred and τ is up to 3.7× the prediction. And the **estimator spread
tracks both**: 1.08–1.21× where the slope is in-spec, **1.86× at M = 200** where it is not.

## 3. The exponent from nested subsets — the kill criterion

| subset | modelled | block | S(0) | L/τ range | cal slope range |
|---|---|---|---|---|---|
| **{10, 20}** | **1.125** | **1.047** | **1.163** | 61–64 | 0.630–0.683 |
| {10, 20, 50} | 1.336 | 1.298 | 1.290 | 40–64 | 0.601–0.683 |
| {10, 20, 50, 100} | 1.391 | 1.404 | 1.376 | 29–64 | 0.574–0.683 |
| all five | 1.460 | 1.541 | 1.405 | 17–64 | 0.476–0.683 |
| {50, 100, 200} | 1.602 | 1.842 | 1.509 | 17–40 | 0.476–0.601 |
| {100, 200} | 1.714 | 2.046 | **1.432** | 17–29 | 0.476–0.574 |

**b rises monotonically as heavier masses are added, in every estimator.** Light against heavy is
3.7σ (modelled), 7.4σ (block), 7.0σ (S(0)) apart. That is the signature named in advance.

**But do not read that as "consistent with linear".** {10, 20} gives 1.125 / 1.047 / 1.163, and a
two-point lever over a factor 2 in M carries σ_b ≈ 0.31 from the 15 % estimator spread alone — so
**b = 1.11 ± 0.31 cannot distinguish 1 from 1.4 and settles nothing.** The meaningful in-spec set is
the three masses that pass *both* pre-registered rules (L/τ ≥ 40 and slope ≥ 0.6): {10, 20, 50}
gives **b = 1.32 ± 0.04 statistically**, and taking the 17 % estimator spread at M = 50 as a
systematic adds ±0.10, so **b = 1.32 ± 0.11 — still ≈ 2.9σ above 1**.

**The honest state is therefore: the heavy end is an estimator limit, the exponent is OPEN, and the
rerun at 65 τ_true decides it.** Neither "b = 1.43" nor "nothing departs from Gruber–Piasecki" is
supportable tonight.

## 4. Degeneracy — this one PASSES

τ_r fixed at 0.5×, 1×, 2× the free value:

| M | τ_T raw at 0.5× | at 1× | at 2× | quoted raw error |
|---|---|---|---|---|
| 100 | 8515 | 8508 | 8504 | ±481 |
| 200 | 24530 | 24520 | 24513 | ±1278 |

τ_T moves by **≤ 10** against errors of 481 and 1278. The slow and oscillatory components are
**not degenerate**, so the amplitude-trading mechanism suspected for the heavy end is ruled out.
The failure is the record length, not the two-component model.

## 5. The independent physical route — no usable answer

The 260928 temperature-step data at M_d = 50 (20 seeds, T = 1.25/0.75) fits a log-linear decay with
**r = −0.464**. That is not an exponential; it is the noise the 260928 report already documented.
The number it returns (9081, against 4060 from equilibrium and 2524 predicted) should not be
quoted, and is not. This route needs the ~2300 seeds that report estimated, which is why the
equilibrium route replaced it.

## 6. Verdict, and the fix

**The heavy end is an estimator limit.** All four diagnostics agree: the calibration slope falls
below its own adoption threshold, L/τ collapses to 17, the three estimators diverge to 1.86×, and b
grows monotonically with the mass range included while the degeneracy check comes back clean.

**For the meeting, unchanged and safe:** τ_T at M_d = 10 matches the parameter-free hard-disk
theory and excludes the ideal gas at 3.5σ; the mode period tracks cot K = αK across α = 0.1–2.0,
with the ideal gas excluded at 7.6–24.1σ; the ladder's mass dependence is under robustness checks
and **the heavy end is currently an estimator limit, not a physical result**.

**The fix is another Mac night.** Records at 65 τ_**true** rather than 65 τ_pred:

| M | have | need | factor | steps/seed | core-hours, 80 seeds |
|---|---|---|---|---|---|
| 50 | 161 639 | 263 898 | 1.6× | 15.8 M | 1.8 |
| 100 | 325 277 | 740 962 | 2.3× | 44.5 M | 4.9 |
| 200 | 652 556 | 2 430 974 | 3.7× | 145.9 M | 16.2 |

**22.9 core-hours, ≈ 2.5 h wall clock at 9 concurrent.** The existing M = 10 and M = 20 cells are
already in-spec and are not repeated. If b is still above 1 with every mass at L/τ ≥ 60 and the
calibration slope ≥ 0.6 throughout, it is a result; on tonight's evidence that is not yet known.

**The standing-wave argument is not invoked.** There is a candidate physical mechanism — in the
standing-wave regime the gas at the divider face moves with the divider, so the relative velocity
driving the asymmetric exchange is smaller than V and α-dependent — but a mechanism for an effect
that has not survived its own robustness checks is not worth writing down as anything yet.
