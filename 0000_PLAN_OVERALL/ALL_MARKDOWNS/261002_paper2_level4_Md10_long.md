# Level 4, M_d = 10 at full length — the 23 % gap is resolved

2026-10-02. `experiments_energy_transfer/level4_equilibrium_KOAlength_20261002/Md10`,
**80 runs, 0 aborts, 0 health events**, 1 964 913 steps each, `--trace-every=200`, release binary.
Record after the t > 2000 burn-in: **30 735 σ-time per seed, L/τ ≈ 61**.

> **τ_T(M_d = 10) = 491 ± 29 (T₁−T₂) and 487 ± 29 (divider x).**
> **Hard-disk isobar 503.5: −0.6σ and −0.7σ — consistent.**
> **Ideal gas 618.4: −3.5σ and −3.6σ — excluded.**

The design asked for 3–4σ at this size and delivered 3.5σ. Everything below was registered before
the data existed; the bias calibration was run on synthetics **before** any real record was fitted.

---

## 1. The premise, again

σ(T₁−T₂) = **0.1975** against the Beta(N,N) prediction **0.199** — −0.7 %, on 80 seeds and a record
four times longer than the run that first confirmed it (0.1980 there). The microcanonical
fluctuation model is not in doubt.

## 2. The mode: KR confirmed, ideal excluded, by a period alone

| | predicted |
|---|---|
| Kolafa–Rottner c_s | **95.16** |
| Paper 1's box c_s (+1.01 %) | 94.21 |
| ideal gas | 117.38 |

| measured (spectral peak, FWHM error) | vs KR | vs Paper-1 c_s | vs ideal |
|---|---|---|---|
| T₁−T₂: **92.92 ± 4.53** | 0.5σ | 0.3σ | **5.4σ** |
| divider x: **95.91 ± 4.70** | 0.2σ | 0.4σ | **4.6σ** |

The ideal-gas period is now excluded at 4.6–5.4σ from the oscillation alone, independently of
anything to do with τ_T. KR and the box c_s remain indistinguishable at this precision, as before.

## 3. Calibration first, then the fit

The estimator's bias was measured on synthetic OU + AR(2) series at **this** record length and seed
count, before the real ACF was touched:

| true τ_T | 350 | 450 | **503.5** | 550 | **618.4** | 700 | 850 |
|---|---|---|---|---|---|---|---|
| recovered | 310.4 ± 12.4 | 379.9 ± 15.8 | **431.2 ± 15.0** | 467.4 ± 25.4 | **516.4 ± 23.9** | 573.0 ± 28.6 | 672.8 ± 36.3 |

**Response slope 0.731**, against 0.438 at the old record length and the **0.6 adoption threshold
set in advance — so the modelled estimator is adopted**, as the pre-registered rule required. The
two hypotheses now map to 431.2 and 516.4, separated by 85 against a scatter of ~15–24: the
estimator can tell them apart, which at 8 200 σ-time it could not.

The fit, ω fixed by cot K = αK:

| observable | A | τ_T (recovered) | B | τ_r |
|---|---|---|---|---|
| T₁ − T₂ | 0.518 ± 0.008 | 419 ± 14 | 0.455 ± 0.006 | 207 ± 4 |
| divider x | 0.515 ± 0.008 | 416 ± 15 | 0.459 ± 0.005 | 206 ± 4 |

The two observables agree on every parameter, and the variance still splits almost exactly in half
between the slow isobaric mode and the resonance.

**Cross-check:** the block estimator, which averages the oscillation away instead of modelling it,
returns raw 415 and 410 against the modelled 419 and 416 — agreement to ~1 % between two estimators
with different systematics. (Its own calibration was not run: the pre-registered rule adopted the
modelled estimator once the slope cleared 0.6.)

## 4. The verdict

| | τ_T | vs hard disk 503.5 | vs ideal 618.4 |
|---|---|---|---|
| from T₁ − T₂ | **491 ± 29** | **−0.6σ** | **−3.5σ** |
| from divider x | **487 ± 29** | **−0.7σ** | **−3.6σ** |

**The hard-disk correction to the adiabatic-piston isobar is confirmed, and the ideal-gas value is
excluded at 3.5σ, at one mass.** Carrying Z(η) through the isobar multiplies the relaxation rate by
1 + ηZ′/Z = 1.2281, shortening τ_T from 61.84 M to 50.35 M; that 23 % is the whole content of the
measurement, and it is now resolved in the same direction as stage 1, where the offset favoured KR
at 7–8σ.

**What is still not tested:** the exponent. τ_T ∝ M needs the mass ladder, and this is one mass.
Eq. (51)'s drift velocity and Eq. (52) remain untested.

## 5. τ_r — sharper, and still not explained

τ_r = **207 ± 4** and **206 ± 4**, i.e. Δf/f = 0.1467 and Q = 6.8. The errors are now small enough
that the disagreement is firm rather than indicative:

| | τ_r |
|---|---|
| kinetic friction on the divider, 2M_eff/γ | 22 |
| **measured** | **207 ± 4** |
| Mansour piston form at M = 10 (Enskog Γ) | 426 |
| bulk sound absorption, 2/(Γk²) | 4217 |

Still bracketed by the two limits, and **the mode damps 2.1× faster than Mansour's piston form**
(207 against 426) — which is a *broader* line, not a narrower one: Δf/f = 1/(πτ_rν), so a shorter
τ_r is more damping and more linewidth. It is now at 4 σ-time precision rather than 17. This is the one number in Level 4 that disagrees with an independent prediction, and
the mass ladder measures τ_r(M) at the same time as τ_T(M) for free.

## 6. Cost, corrected

I called this "an overnight job, about 13× the last run". It ran **80 seeds in 13 minutes**, at
~10 s per seed, with 0 health events. The estimate was inferred from step count without ever having
measured this box's rate, and it was wrong by roughly 50×.

The consequence is worth having before the meeting: the remaining masses scale as 5× and 20× the
steps, so **M_d = 50 is ≈ 1.1 h and M_d = 200 is ≈ 4.4 h of single-threaded Mac time for 80 seeds
each**. The whole Level 4 mass ladder — the measurement that tests the exponent — is an overnight
Mac job. That does not make KOA idle (Paper 1's N ≈ 1000 sweep is genuinely large), but Level 4 is
no longer a reason to wait for it.

## 7. Not done, deliberately

The plan's Level 4 verdict and the Paper 2 draft are **unchanged**. The wording of a verdict this
strong is worth settling with the Level 4 scope tomorrow rather than committing tonight. Proposed
wording, for that conversation:

> Level 4 (ii): with the equilibrium fluctuation route at 61 τ_T and 80 seeds, τ_T(M_d = 10) =
> 491 ± 29 σ-time, consistent with the hard-disk isobar (50.35 M → 503.5, −0.6σ) and excluding the
> ideal-gas value (61.84 M → 618.4) at 3.5σ. Two observables and two estimators agree. The exponent
> in M remains untested.
