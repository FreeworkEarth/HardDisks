# Level 3 v4 — what the residual is

2026-09-22. Runs `experiments_energy_transfer/level3_v4_20260922` (A1 160, A4 200) plus the v3 grid.
**0 aborts, 0 health events throughout.** Analysis `validation/paper2_level3_v4_20260922.py`.

**Headline: a bug in v3 is corrected here, and criterion (iv) selects A1 (gas dynamics) — but not
by the test as written.** The statics are right to 8–10 %; the whole residual is in the transient.

---

## 0. BUG IN v3, corrected — v3's model numbers were 28 % too weak

The adiabatic force table normalised the adiabat at the **left edge of the interpolation grid**
instead of at L:

```python
_lnT = np.concatenate([[0.0], np.cumsum(...)])    # lnT = 0 at _LG[0], NOT at L
```

Every force therefore carried a constant factor exp(−Z̄ ln(L/L_grid0)) ≈ 0.72. Because F₀ carries
it too, the *difference* F_ad(L′) − F₀ that drives the ODE was ~28 % too weak. Effect: the
quasi-static fixed point read **0.625 σ instead of 0.840 σ**, and the v3 model looked low by more
than it was. Fixed by `_lnT -= np.interp(L, _LG, _lnT)`; verified T(L) = 1.000000 and the fixed
point now solves the same equation as the standalone calculation (0.8404 σ → 1.500 kT). **The v3
report's model column and its 2.8–4.2σ distances are superseded by this document.**

## 1. A1 — the force the gas actually delivers

Spring wall held rigid (M = 1e9), momentum per collision read from the event log (`kind D0`,
column `dp`), binned at Δt = 0.5, 40 seeds.

| u | before the signal can arrive | late in the push | after the stop |
|---|---|---|---|
| 0.05 | 1.5961 | 1.6175 | 1.8533 |
| 0.1 | 1.5961 | 1.5860 | 1.7319 |

against the standing force F = N kT Z/L = **1.5749** — so the pre-signal window reproduces the
static force to 1.3 %, an independent check that the held-wall force measurement is sound.

Driving the ODE with F_meas instead of F_ad:

| u | M_s | measured ΔE_coh | ODE v3 (F_ad) | ODE A1 (F_meas) | ODE A2 (M̂) | σ to A1 | σ to A2 | σ to v3 |
|---|---|---|---|---|---|---|---|---|
| 0.05 | 50 | 2.819 ± 0.456 | 1.672 | 4.182 | 1.572 | 3.0 | 2.7 | 2.5 |
| 0.05 | 200 | 2.399 ± 0.356 | 1.859 | 2.160 | 1.805 | **0.7** | 1.7 | 1.5 |
| 0.1 | 50 | 2.981 ± 0.520 | 1.859 | 3.490 | 1.600 | 1.0 | 2.7 | 2.2 |
| 0.1 | 200 | 3.601 ± 0.456 | 2.218 | 2.356 | 2.341 | 2.7 | 2.8 | 3.0 |

A1 does **not** pass the stated <1σ test. But its failure has a direction and a known cause: it
**overshoots at M_s = 50 and not at M_s = 200**. A rigid boundary doubles the pressure of an
incident pulse, so a held wall measures roughly twice the force a *moving* wall would feel — and
the lighter the wall, the more it moves and the worse that overestimate. Feeding a held-wall force
into a moving-wall ODE is therefore biased high, most for light walls, which is exactly the pattern
observed. Correcting it requires the wall's radiation reaction, i.e. the impedance term — which
belongs here as a correction to the *drive*, not, as v2 assumed, as the capture mechanism itself.

## 2. A2 — gas inertia is not the answer

M̂ = M_s + Nm/3 (Mansour–Garcia–Baras Eq. 18) moves the prediction by less than the error bar:
1.672 → 1.572 at M_s = 50, 1.859 → 1.805 at M_s = 200. Worst distance 2.8σ, essentially unchanged
from the v3 baseline's 3.0σ. **Rejected.**

## 3. A3 — the slow point: the ratio does not fall

| u | τ_push/(L/c_s) | M_s | measured ΔE_coh | ODE | ratio |
|---|---|---|---|---|---|
| 0.02 | 8.84 | 50 | 2.743 ± 0.613 | 1.566 | 1.75 |
| 0.02 | 8.84 | 200 | 2.636 ± 0.460 | 1.644 | 1.60 |
| 0.05 | 3.54 | 50 | 2.819 ± 0.456 | 1.672 | 1.69 |
| 0.05 | 3.54 | 200 | 2.399 ± 0.356 | 1.859 | 1.29 |
| 0.1 | 1.77 | 50 | 2.981 ± 0.520 | 1.859 | 1.60 |
| 0.1 | 1.77 | 200 | 3.601 ± 0.456 | 2.218 | 1.62 |

Flat at ≈1.6 from 1.8 to **8.8 sound traversals**. By the criterion as written this favours A2 over
A1 — but A2 cannot produce the magnitude, so the criterion as written does not resolve it. §5 does.

## 4. A4 — the thermalised baseline, and τ_heat(M)

Wall free from t = 0, **no push at all**, 2000 σ-time, 40 seeds. Targets: static pre-load
F²/2k = 2.480 kT, and 1 kT of equipartition above it (½kT in the spring coordinate, ½kT kinetic).

| M_s | ⟨E_spring⟩ plateau | ⟨KE_wall⟩ plateau | sum − pre-load | τ_heat [σ-time] |
|---|---|---|---|---|
| 2 | 2.965 | 0.496 | 0.981 | 1 |
| 10 | 3.019 | 0.510 | 1.049 | 4 |
| 50 | 3.028 | 0.548 | 1.095 | 48 |
| 200 | 2.974 | 0.453 | 0.947 | 154 |
| 1000 | 2.632 | 0.084 | **0.235** | not reached |

**Equipartition confirmed to 5–10 % for M_s ≤ 200.** M_s = 1000 does not thermalise in 2000 σ-time,
so its v3 row rests on an unthermalised wall and should be treated as provisional. τ_heat grows
with mass — 1, 4, 48, 154 σ-time — and is logged here as the first τ_heat(M) number for Level 4.

**How much of the baseline line disappears: almost none.** A4 is a proper no-push control, so it can
be subtracted *time-matched* rather than as an assumed plateau. Doing that:

| u | M_s | blind-window baseline | time-matched control | shift |
|---|---|---|---|---|
| 0.02 | 50 | 3.065 | 2.691 ± 0.623 | −0.374 |
| 0.02 | 200 | 2.677 | 2.717 ± 0.460 | +0.040 |
| 0.05 | 50 | 3.141 | 3.048 ± 0.541 | −0.093 |
| 0.05 | 200 | 2.440 | 2.422 ± 0.316 | −0.018 |
| 0.1 | 50 | 3.303 | 3.427 ± 0.662 | +0.123 |
| 0.1 | 200 | 3.642 | 3.428 ± 0.346 | −0.214 |

Shifts of −0.37 to +0.12 kT with **no consistent sign**. The v3 report's "0–0.46 kT baseline
deficit" is therefore real in size but not systematic, and it is not the residual.

## 5. Criterion (iv): the statics are right, the transient is not

The sharpest diagnostic is to compare the *settled* and the *peak* wall displacement separately:

| u | M_s | settled s, measured | s_qs model | ratio | peak s, measured | peak s, model | ratio |
|---|---|---|---|---|---|---|---|
| 0.02 | 50 | 0.904 | 0.840 | **1.08** | 1.367 | 0.873 | **1.57** |
| 0.02 | 200 | 0.901 | 0.840 | **1.07** | 1.241 | 0.912 | **1.36** |
| 0.05 | 50 | 0.911 | 0.840 | **1.08** | 1.434 | 0.925 | **1.55** |
| 0.05 | 200 | 0.922 | 0.840 | **1.10** | 1.214 | 1.016 | **1.19** |

**The model gets where the wall ends up right to 8–10 %. It gets how far it swings on the way there
wrong by 19–57 %.** The residual is entirely in the transient, not in the equation of state, not in
the pre-load, not in the baseline.

That rules out A2 and rules out a static error, and it points at A1's mechanism: the ODE is driven
by the quasi-static pressure history, and the wall is receiving more impulse than that history
supplies. The extra impulse is the acoustic pulse. Note the model is **undamped**, so it already
gives the largest overshoot a linear 1-DOF system can produce — real gas damping would *lower* the
predicted peak and widen the gap. The excess cannot be explained by anything that removes energy.

**Verdict: (iv) selects A1 (gas dynamics at the far wall), on the statics/transient split rather
than on the F_meas test, which is biased high for light walls by rigid-boundary pressure doubling.
Level 3 remains NOT passed:** criteria (i) and (ii) hold as in v3, and (iii) still fails, now at
1.5–3.0σ rather than 2.8–4.2σ after the §0 bugfix.

## 6. What would close it

The drive needs the radiation reaction of a moving wall, not the held-wall force: F_meas corrected
by the wall's own impedance, or equivalently a two-variable model carrying the first acoustic mode
of the gas alongside the wall coordinate. That is one more ODE, no new runs — the A1 held-wall force
and the A4 control are both already measured.

The plan's §Level 3 verdict is **not** changed by this document beyond the bugfix note; the
apparatus, the pre-loaded ceiling and the adiabatic gas spring all stand.
