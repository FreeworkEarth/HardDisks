# Level 4 thermal — the adiabatic piston, no piston

2026-09-28. `experiments_energy_transfer/level4_thermal_20260928`, **40 runs, 0 aborts, 0 health
events**. T₁ = 1.25, T₂ = 0.75, N = 50 each, L_c = 38.75, divider free, no piston motion at all.
M_d = 10 (3000 σ-time) and 50 (14000 σ-time), 20 seeds each. Predictions written into the plan
before the runs.

**Stage 1 is observed and roughly as predicted. Stage 2's τ cannot be measured with 20 seeds, and
the report says so rather than quoting the number the fitter returns.**

---

## 0. No new flag was needed, and an earlier report was wrong

`--temperature` already accepts a **comma-separated list**, one value per segment
(`00ALLINONE.c:4579` → `cli_parse_float_list` → `apply_segment_temperatures` at `:6271`, applied
during initialisation right after the per-segment equaliser). The earlier statement that
per-compartment temperatures could not be set was wrong: it was based on the help text, which
documents only the scalar form, and on the signature of `equalize_temperature_per_segment(float)`.
A duplicate implementation was written and then discarded when it failed to compile against the
existing function of the same name; the source was reverted and **no code change was made**.

**Units trap, measured.** With `--kbt1` the runtime temperature is 100 and k_B is scaled so
k_BT = 1, and the list is in those raw units:

| flag | resulting kT |
|---|---|
| `--temperature=1.25,0.75` | 0.0125 / 0.0075 — right ratio, 80× wrong scale |
| `--temperature=125,75` | **1.2500 / 0.7500**, ratio 1.6667 ✓ |

Same family as the `--spring-k` per-pixel trap.

## 1. Two corrections to the predictions

**The stage-1 offset is 4.92 σ ideal, not 5.3.** With the 2D adiabat TL = const and pressure
balance, offset = L_c(r−1)/(r+1), r = √(T₁/T₂) = 1.29099, giving **4.922 σ** for L_c = 38.75. (5.3
would require L_c = 41.7.) Carrying KR — adiabat d ln T = −Z d ln L, pressure balance with Z(η) —
gives **3.989 σ, −18.9 % against the ideal gas**, and that is the value to compare against. The
temperatures at the end of stage 1 are T₁ = 1.108, T₂ = 0.860, so the difference entering stage 2
is 0.249 — **exactly half** the initial 0.5, as anticipated.

**The √T nonlinearity is 0.2 %, not 4 %.** At the start of stage 2 the exact Eq. (51) differs from
its linearisation by 0.20 % (0.80 % at the initial 1.25/0.75). Integrating the exact ODE and fitting
an exponential over 3τ returns 507.5 and 2537.3 σ-time at M_d = 10, 50 against the analytic 507.6
and 2537.9 — under 1 %. A pure exponential fit is safe. (Fitting over 40τ returns a spurious 1.4×;
the tail underflows and its logarithm is noise.)

## 2. (a) Stage 1 — observed

| M_d | R = mN/M | regime | max offset | KR prediction | ideal |
|---|---|---|---|---|---|
| 10 | 5.00 | over-damped | +6.38 σ | 3.99 | 4.92 |
| 50 | 1.00 | critical | +6.97 σ | 3.99 | 4.92 |

The divider swings out and rings, overshooting the predicted equilibrium by 1.6–1.7×, which is what
the damping regime implies. A *settled* stage-1 plateau cannot be quoted: stage 2 runs the offset
back to zero, and at M_d = 50 the ringing has not died before it does. The divider position itself
carries thermal noise of rms √(kT/2k_gas) = 2.2 σ per seed, i.e. 0.49 σ on a 20-seed mean.

## 3. (b) Stage 2 — NOT MEASURABLE at 20 seeds

The fitter returns τ = 670 ± 3 and 5615 ± 9, and the ratio τ(50)/τ(10) = 8.38 ± 0.04 against 5.0
for τ ∝ M. **None of those numbers should be believed.** Three checks say so.

**The curve does not decay.** Seed-averaged T₁ − T₂ at M_d = 10: 0.500, 0.073, 0.329, 0.194, 0.168,
0.035, 0.069, −0.044 at t = 0, 50, 100, 200, 500, 1000, 2000, 3000. It bounces, and goes negative.

**The fit is window-dependent.** At M_d = 10, τ = 777, 865, 846, 1359, 2010 as the window moves
later; at M_d = 50, 5576, 8299, 8699, 10678, 12857. A quantity that doubles with the fit window is
not a measurement.

**Local decay rates change sign.** In successive blocks at M_d = 50: 5459, 5741, **−4553**, 1880,
**−18794**.

**The reason is the noise floor, and it was predictable.** For N_s = 50 in 2D the kinetic
temperature has σ_T/T = 1/√N_s = 0.141, so σ(T₁−T₂) ≈ √2 × 0.141 = 0.200 per seed. **Measured:
0.204 and 0.179.** On 20 seeds that is 0.045 on the mean, and the signal entering stage 2 is only
0.25, so it exceeds the noise for

> t < τ ln(0.25/0.045) = **1.7 τ**

— less than two time constants, which cannot pin an exponential. Tracking 3τ at 3σ needs

> **≈ 2300 seeds**, against the 20 run here.

## 4. (d) Divider kinetic energy

⟨KE_wall⟩ late = 0.5057 ± 0.1459 (M_d = 10) and 0.4841 ± 0.1412 (50), against equipartition 0.5.
**Consistent.** It does *not* discriminate Eq. (52): by that time T₁ ≈ T₂, so ½√(T₁T₂) = 0.4999 and
½kT = 0.5 are the same number. Testing Eq. (52) requires the stage-2 window where T₁ ≠ T₂, and there
the single divider's kinetic energy is itself too noisy on 20 seeds.

## 5. What this run established, and what it did not

**Established:** the apparatus creates a clean temperature step (1.2500 / 0.7500 measured); stage 1
happens, with the divider swinging out and ringing in the regime R predicts; the KR correction to
the stage-1 offset is −18.9 %; equipartition of the divider holds at the end; and the noise model
for this observable is confirmed quantitatively (0.204 measured against 0.200 predicted).

**Not established:** τ_T, its scaling with M, Eq. (51)'s drift velocity, and Eq. (52). All four need
the signal tracked over several τ, which needs ~100× more seeds.

**For KOA:** M_d ∈ {10, 50, 200}, **≈2000 seeds each**, records 5τ from τ_T = 50.4 M. Alternatively
raise N_s: the noise falls as 1/√N_s, so N_s = 500 buys a factor √10 at the cost of changing the
gas. The seed count is the cheaper knob and parallelises perfectly.

**The M² claim stays withdrawn** — not because this run refutes it, but because the literature is
explicit that stage 2 is linear in M and this run cannot test either exponent.
