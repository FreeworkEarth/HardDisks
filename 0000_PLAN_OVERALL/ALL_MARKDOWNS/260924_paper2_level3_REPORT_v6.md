# Level 3 v6 — the settled-state comparison, closed

2026-09-24. Post-processing, plus 80 re-runs of two cells whose records were too short
(`level3_v6_20260924`, **0 aborts, 0 health**). Analysis `validation/paper2_level3_v6_20260924.py`.

**Level 3 passes on ε_settled with the box's own measured inputs: −0.08 ± 0.74 %, 0.1σ.** With
textbook bulk inputs it is +3.58 %, 6.2σ. The difference between those two lines is the whole
result, and every input is a measurement with an error bar — nothing is fitted.

---

## 1. Window validity, and one more invalid cell

The settled average needs the wall to have completed several free periods,
T_w = 2π√(M_s/(k + k_gas)) = **60 σ-time at M_s = 50, 120 at 200**. A window shorter than ~3 periods
is a phase sample, not a settled state.

| u | M_s | T_w | window start | record end | periods | verdict | s̄ corrected |
|---|---|---|---|---|---|---|---|
| 0.02 | 50 | 60 | 578 | 1250 | 11.2 | valid | 0.871 ± 0.008 |
| 0.02 | 200 | 120 | 758 | 1250 | 4.1 | valid | 0.858 ± 0.013 |
| 0.05 | 50 | 60 | 339 | 500 | **2.7** | **INVALID** | — |
| 0.05 | 200 | 120 | 519 | 1250 | 6.1 | valid | 0.882 ± 0.011 |
| 0.1 | 50 | 60 | 259 | 500 | 4.0 | valid | 0.847 ± 0.023 |
| 0.1 | 200 | 120 | 439 | 1250 | 6.8 | valid | 0.872 ± 0.010 |

The two M_s = 200 cells that v5 reported as missing/1.7σ were re-run to 1250 σ-time and are now
valid. A cell v5 *did* count — u = 0.05, M_s = 50 at 2.7 periods — falls just below the bar and is
now excluded; v5's 1.3σ for it should not have been quoted.

**Pooled over the 5 valid cells, inverse variance: s̄ = 0.8705 ± 0.0049 σ.**

## 2. Bulk inputs versus the box's own

| model inputs | s_qs | vs pooled s̄ | σ | ΔE(s_qs) [kT] |
|---|---|---|---|---|
| bulk Kolafa–Rottner (the v3–v5 model) | 0.8404 | +3.58 % | **6.2** | 1.5000 |
| (a) + held-wall standing force, additive | 0.8404 | +3.58 % | 6.2 | 1.5179 |
| (a)+(b) + Paper 1's c_s for this box | 0.8552 | +1.78 % | 3.1 | 1.5479 |
| **(a)+(b)+(c) + Level 1's measured adiabat** | **0.8711** | **−0.08 %** | **0.1** | 1.5801 |

Inputs, all measured:
- **(a)** A1 held wall, F_box = 1.5961 against the EOS value 1.5749, **+1.35 %**.
- **(b)** Paper 1 at η = 0.1001 (interpolated between its bracketing densities),
  c_s(box)/c_s(bulk) = **1.0101 ± 0.0027**, so the gas spring k_gas ∝ c_s² is **+2.02 ± 0.55 %**.
- **(c)** Level 1 measured T_f/T_i = **1.149** against the bulk adiabat's 1.1432, **+0.51 %**.

**Sensitivities, stated explicitly.** ∂ln s/∂ln k_gas = k/(k + k_gas) = **0.910**.
∂ln s/∂ln F for an **additive** offset is **exactly 0** — the fixed point k s = F(L−Δx+s) − F(L)
depends only on the force *difference*, so a constant cancels. A *multiplicative* scaling is a
different object: it scales the slope too and would give ∂ln s = 0.910 × ∂ln F.

**Why (a) is additive, and a correction to v5.** v5 offered the +1.35 % standing force as "a third
to a half" of the residual. That was wrong twice over: a constant offset cancels from the fixed
point, and the excess is not a bulk stiffening anyway. Paper 1 established that the wall pressure
Z_wall is a **surface contact value** — measured at Paper 1's own densities it over-predicts the c_s
offset five-fold, and the claim that it explained the N = 100 offset was **withdrawn on 2026-09-19**.
Scaling the force curve by F_box/F_EOS, and hence its slope, would have repeated that error. It is
applied additively, it leaves s_qs untouched (row 2 of the table proves it), and it enters only the
energy ΔE = F s + ½k s², where F is the actual standing force.

**Result:** box-input model s_qs = 0.8711 ± 0.0041, measurement 0.8705 ± 0.0049, residual
**−0.08 ± 0.74 %, 0.1σ.**

**Do not over-read the 0.1σ.** The total correction is +3.6 %, built from two independent-looking
measurements worth +1.78 % and +1.86 %, and their sum lands on the measurement to better than the
inputs' own precision. That is partly luck. Two caveats:

1. **(b) and (c) are not proven independent.** Both could be manifestations of the same finite-size
   stiffening. They are not *simply* the same, though: a uniform excess δ in Z would need δ = 1.34 %
   to produce (b) and δ = 3.78 % to produce (c) (sensitivity d ln c_s²/d ln Z = 1.502), so they
   cannot both be one uniform Z excess.
2. **(c) is transferred across geometries.** Level 1's T_f/T_i was measured with 50 disks in
   geometry A; it is applied here to 100 disks in geometry C. T_f/T_i is intensive in N at fixed
   compression ratio, so this should hold, but it is an assumption. A second assumption: the adiabat
   is set using the full Δx, whereas the wall moves, so the gas is actually compressed by Δx − s.

## 3. Consistency check on the pre-push offset

| | |
|---|---|
| predicted from the held-wall force, (F_box − F_EOS)/(k + k_gas) = 0.0212/0.5494 | **0.0386 σ** |
| measured, A4 no-push control mean | **0.0503 σ** (M_s = 50: 0.0427; 200: 0.0578) |

Two independent measurements of the same +1.35 % excess standing force, agreeing to the spread
between the two masses. Note this offset does **not** enter the settled residual — the control
removes it, and an additive F offset cancels from the fixed point regardless.

## 4. Verdict

> The parameter-free model with bulk EOS inputs predicts the settled wall displacement to
> **+3.6 ± 0.6 %** (five cells pooled, 6.2σ). With the box's own measured standing force, its
> Paper-1 sound speed and Level 1's measured adiabat as inputs — all measurements, none fitted —
> the residual is **−0.1 ± 0.7 %, 0.1σ**. The transient excess reported in v4 was a peak statistic,
> reproduced by a no-push control.

**Level 3 passes on ε_settled** (criterion: ≤2σ on the box-input model). The peak statistic is
retained as a diagnostic only and is dropped from the pass criteria.

Both papers now say one thing: **the N = 100 box is about one per cent stiffer than bulk in c_s, and
once you know that, its energy transfer is predicted.**

## 5. B2 — decided, not carried again

The canonical Paper 1 traces carry `Time, Wall_X, Displacement, Left_Count, Right_Count, L0, eta,
Center_X, Seed, Target_Oscillations, Predicted_Frequency, Planned_Steps, Planned_Duration` — **particle
counts, no per-compartment kinetic temperatures and no release snapshot**. The regression over
canonical seeds therefore cannot be done on stored data and is **not carried a fourth time**.

What does exist is the dedicated 5-seed test already in the method note §7.2 (L0 = 20, 200 periods,
with KE_L/KE_R logged): correlation of the slow mode with the compartment temperature difference
r = +0.978, amplitude reproduced to 6–7 % by x = (L₀/2)(ΔT/T)·Z/(Z+ηZ′). The mechanism is therefore
directly evidenced; only the regression over the canonical campaign is missing. `log T_1, T_2 at
release` is added to the KOA campaign list.

## 6. Radiation reaction — still not built

Unchanged from v5: it adds damping, reduces overshoot, and the residual is now 0.1σ anyway. Level 4
item, wall–gas coupling Z_g = N m c_s/L = 2.22, Q ≈ ω_w M/Z_g.
