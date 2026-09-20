# Level 3 v3 — the pre-loaded spring and the one-degree-of-freedom model

2026-09-21. Supersedes v1 (`260920_...`) and v2 (`260921_paper2_level3_REPORT.md`). Runs:
`experiments_energy_transfer/level3_master_preload_20260921`, 760 trajectories, **0 aborts, 0 health
events**. Analysis `validation/paper2_level3_v3_20260921.py`.

**Verdict: Level 3 is NOT passed.** The apparatus and the ceiling are now right, the model is
parameter-free, and it lands within a factor 1.5–2.5 of the data at the quasi-static end — 2.8σ to
4.2σ, not the 2σ required. The residual is a systematic **+1.0 to +1.9 kT excess of measurement over
model**, partly accounted for and partly not (§5).

---

## 1. The ceiling, rewritten: the apparatus is a pre-loaded spring

A spring holding a wall against a gas is not at its natural length: it carries F = N kT Z/L, so the
wall sits d = F/k from the anchor and

> **ΔE(s) = F·s + ½k s² — first order in s.**

The F·s term is work lifting a weight, which is what a Szilard-type engine is judged by. The
pre-load cancels from the *dynamics* (at s = 0 the spring balances the gas) but not from the
*energy*. Limits: k → 0 is a pushrod, the wall follows the piston, ε → 1 trivially and nothing was
captured *from the gas*; k → ∞ is a rigid wall, ε → 0. In between, quasi-static capture is ΔE(s_qs)
with s_qs = Δx·k_gas/(k+k_gas) = 0.715 σ → **1.255 kT**.

The unloaded derivation is retained only as the explanation of why the pilot could not pass: without
the F·s term, ΔE = ½ks² is second order, maximal at k = k_gas with N kT f² c_s²/8 = 0.098 kT against
a kT/2 pedestal, and that optimum cannot be built anyway because F/k_gas = L·Z/(Z+ηZ′) = **0.816 L
independent of N**.

## 2. The gas spring is adiabatic — and it is Paper 1's c_s

There is no heat bath, so compression heats the gas. With U = N kT in 2D, dU = −P dV gives
d ln T/d ln L = −Z, and

> **k_gas^ad = −∂F/∂L|_S = (N kT/L²)(Z + ηZ′ + Z²) = (N kT/L²)·c_s²**

— exactly Paper 1's measured quantity. Here 0.04938 against the isothermal 0.02457, a factor
**2.01**. Using the isothermal stiffness put the model a factor 2.3 below the measured wall
displacement; this was the single largest error in the v2 analysis. Independent check of the adiabat:
it predicts T_f/T_i = 1.1432 for this 10 % compression against **Level 1's measured 1.149**, agreeing
to 0.5 %.

## 3. Baseline — the question answered

`ΔE_coh` **does** subtract the baseline **mean**, per seed, not merely report its rms. Proof: for
M_s = 200, u = 0.2 the raw ensemble-mean peak is 8.315 kT, the baseline mean is 2.560, and the table
reports 5.755 = 8.315 − 2.560. So the suspicion that every ε carried a ~3 kT offset is not borne out
and the v2 numbers stand as ΔE.

Two estimator defects *were* found in the course of checking:

- **Peak-picking bias.** `argmax` over ~18 000 noisy time points is biased high. Measured with the
  same estimator on the causally blind window, where the true value is exactly zero: 0.01–0.72 kT
  depending on cell. Now subtracted, self-calibrating, no model.
- **The blind window was not blind at high u.** At 0.8 L/c_s the "bias" came out at 12 kT at u = 1 —
  a strong push outruns c_s. Tightened to 0.4 L/c_s, which needs a 2.5 c_s front to leak.

## 4. The measurement against the model

Eq.: M s̈ = F_ad(L − x_p + s) − F_ad(L) − k s, x_p = min(ut, Δx), s(0) = ṡ(0) = 0, RK4 at dt = 1e-3,
exact adiabatic F_ad. Nothing is fitted.

| u | τ_push/(L/c_s) | M_s | τ/T_w | ΔE_coh corrected | 1-DOF | diff | σ |
|---|---|---|---|---|---|---|---|
| 0.05 | 3.54 | 50 | 2.66 | 2.819 ± 0.456 | 1.265 | +1.554 | 3.4 |
| 0.05 | 3.54 | 200 | 1.33 | 2.399 ± 0.356 | 1.388 | +1.011 | 2.8 |
| 0.1 | 1.77 | 50 | 1.33 | 2.981 ± 0.520 | 1.388 | +1.593 | 3.1 |
| 0.1 | 1.77 | 200 | 0.66 | 3.601 ± 0.456 | 1.670 | +1.932 | 4.2 |
| 0.2 | 0.88 | 2…1000 | 3.32…0.15 | 4.28…4.35 | 1.23…2.42 | +1.9…+3.7 | 4.1…14.8 |
| 0.5 | 0.35 | 2…1000 | 1.33…0.06 | 5.81…10.66 | 1.39…2.46 | +3.3…+8.4 | 9.2…17.8 |
| 1.0 | 0.18 | 2…1000 | 0.66…0.03 | 9.46…18.84 | 1.67…2.47 | +7.0…+16.4 | 10.0…14.6 |

**Criterion (iii): FAIL.** 0/4 of the genuinely quasi-static cells are within 2σ; the best is 2.8σ.

Note the domain correction: u = 0.2 was chosen for the quasi-static test but
τ_push/(L/c_s) = 0.88 < 1 — the sound does not cross the box during the push, so a single gas
coordinate cannot stand in for the field. Only u ≤ 0.1 qualifies, which is why the extra speeds were
run.

## 5. What is and is not explained

The excess is **systematic and positive** everywhere, and it grows with u exactly as an acoustic
term should — from +1.0 kT at u = 0.05 to +16 kT at u = 1.0, against A u² = 0.06 and 25.9 kT. At the
fast end the deviation is the acoustic contribution and is reported as such, with no pass claimed.

At the quasi-static end, where there should be nothing left, +1.0 to +1.9 kT remains. Identified so
far:

- **Baseline thermal deficit.** The wall is released from rest, so it carries no thermal energy at
  t = 0 and the blind window samples a partly unthermalised state. Measured baselines run 2.98
  (M_s = 2, thermalises fast) down to 2.477 (M_s = 1000) against an equilibrium 2.935 kT, so ΔE_coh
  is biased high by **0 to 0.46 kT**, worst for the heavy walls. Stated rather than corrected:
  correcting would mean assuming an equilibrium the run has not reached.
- That leaves roughly **0.5 to 1.2 kT unexplained**, and I do not have a mechanism for it. The
  obvious next candidate is gas inertia — the 1-DOF model treats the gas as a massless spring, but
  N m = 100 is comparable to M_s, and Paper 1's own linewidth section already carries Mansour's
  M̂ = M + mN/3 for exactly this reason.

## 6. Honest status of the criteria

| | |
|---|---|
| (i) ledger with the spring term at 1e-8 | met in substance — median 4.64e-10 over 600 runs, 593/600 below 1e-8, exceptions at M_s = 2 |
| (ii) coherent peak > 3σ above the pedestal | **PASS** — 6.8σ to 31σ in all cells |
| (iii) parameter-free model within 2σ where the gas is quasi-static | **FAIL** — best 2.8σ, systematic +1.0…+1.9 kT |

Level 3 is closer than it has been, on an apparatus and a ceiling that are now right, with the gas
spring tied to Paper 1's measured c_s. It is not closed.
