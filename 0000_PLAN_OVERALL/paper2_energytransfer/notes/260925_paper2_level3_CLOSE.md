# Level 3 — closed on the observable; the quantitative close follows the measured F(L)

2026-09-25. Post-processing. Answers the direct question put to v6 and revises its verdict.

---

## 1. The question: which L does the fixed-point adiabat use?

**The nominal L = 78.5**, not an effective length. `paper2_level3_v6_20260924.py` line 40:
`N, L, H, R, DX, K = 100, 78.5, 10.0, 0.5, 7.96, 0.5`.

Geometry: wall_S centre 30.5, thickness 1.0, so its gas-side face is at 31.0; the box ends at
109.50; the gas region is 109.50 − 31.0 = **78.5 σ**. Disk *centres* span 31.5 … 109.0 = 77.5 σ.

**And nominal is the correct choice.** Kolafa–Rottner's η is defined on the physical area, and
Paper 1 uses it that way — `tests_20260913.eta_of(L0)` returns `N_TOTAL·πr²/(2·L0·H)` with no
exclusion, while `l_eff(L0) = L0 − 2r − t/2` is defined separately and enters only the *acoustic*
length of the divider mode, never the density. So substituting L_eff into the density would not be
a change of convention; it would be an ad-hoc finite-size correction wearing a convention's clothes.

## 2. L_eff does not reproduce contribution (c)

| | change in s_qs |
|---|---|
| (c) as implemented in v6 — rescale the adiabat to Level 1's T_f/T_i = 1.149 | **+1.86 %** |
| switching the density length 78.5 → 77.5, bulk adiabat throughout | **+3.05 %** |

They are not the same number, so the hypothesis that (c) *is* the effective-length correction in
disguise is **not supported**. (For the record, the natural T_f/T_i moves only 1.1432 → 1.1458,
+0.23 %, on that switch, against Level 1's measured +0.51 %.)

## 3. The reason the attribution is open

Substituting exclusion lengths into the density is inadmissible (§1), but doing it *anyway* shows
how weakly the comparison constrains the answer. Fixed point with bulk KR and no other input:

| length used for the density | L | H | η | s_qs | vs measured 0.8705 ± 0.0049 | σ |
|---|---|---|---|---|---|---|
| **nominal — the admissible one** | 78.5 | 10 | 0.10005 | 0.8404 | +3.59 % | **6.1** |
| exclude r in x | 77.5 | 10 | 0.10134 | 0.8660 | +0.52 % | 0.9 |
| exclude r in x and y | 77.5 | 9 | 0.11260 | 0.9080 | −4.13 % | 7.6 |
| exclude r in y | 78.5 | 9 | 0.11117 | 0.8807 | −1.16 % | 2.1 |

Rows 2–4 are **not** alternative conventions; they are ad-hoc finite-size corrections, and they are
listed only to show their span: **8 %**, against a residual of 3.6 % that v6 attributed to two
measured inputs. An ad-hoc correction of the right size and sign is available for free. That is why
v6's "−0.08 %, 0.1σ" cannot be read as an attribution, and why its "partly luck" caveat was too
weak. The measurement and the model agree to a few per cent; *what* closes the few per cent is not
pinned by this comparison.

## 4. Verdict

> With the bulk equation of state the parameter-free model predicts the settled wall displacement to
> **+3.6 ± 0.6 %**. The residual is the finite-size stiffness of the 100-disk box, the same effect
> Paper 1 measures as +1 % in c_s, but it is not independently pinned here, so no attribution is
> claimed; the transient excess reported earlier was a peak statistic reproduced by a no-push
> control. **Level 3 closes on the observable; the quantitative close follows from the measured
> F(L) set.**

The peak statistic is retained as a diagnostic only.

## 5. What closes it: F(L) measured, not inferred

`_run_scripts/level3_FofL_20260925.sh`, **run 2026-09-25: 100 trajectories, 0 aborts, 0 health**. Geometry C with the spring wall
held at 1e9; the piston compresses by 0 / 2.5 / 5 / 7.5 / 10 % and stops; after the transient the
gas is in equilibrium at that length and the momentum delivered to the held wall per unit time is
the force. Measuring the gas temperature at the same instant (2D: U = N kT, so T = KE_gas/N) gives

> **Z_box(η) = F·L / (N k T)**

with no adiabat, no equation of state and no area convention assumed. 20 seeds per compression.
### Result

| compression | L [σ] | η | F measured | T measured | Z_box | Z_KR | Z_box/Z_KR |
|---|---|---|---|---|---|---|---|
| 0 % | 78.50 | 0.10005 | 1.6087 ± 0.0064 | 1.0000 | 1.2628 ± 0.0050 | 1.2363 | 1.0215 |
| 2.5 % | 76.51 | 0.10265 | 1.7218 ± 0.0078 | 1.0347 | 1.2732 ± 0.0057 | 1.2435 | 1.0238 |
| 5 % | 74.52 | 0.10539 | 1.8366 ± 0.0067 | 1.0699 | 1.2792 ± 0.0046 | 1.2513 | 1.0224 |
| 7.5 % | 72.53 | 0.10829 | 1.9715 ± 0.0059 | 1.1070 | 1.2917 ± 0.0039 | 1.2595 | 1.0256 |
| 10 % | 70.54 | 0.11134 | 2.1217 ± 0.0042 | 1.1499 | 1.3016 ± 0.0026 | 1.2683 | 1.0262 |

**Z_box/Z_KR = 1.0248 ± 0.0014 (+2.48 ± 0.14 %)**, and notably flat across the range.

**An independent reproduction of Level 1 falls out for free.** The measured gas temperature at 10 %
compression is **T = 1.1499**, against Level 1's measured T_f/T_i = 1.149 — from a different
geometry, a different observable and a different analysis.

**Fixed point from the measured adiabat: s_qs = 0.8819**, against the pooled measurement
0.8705 ± 0.0049: **−1.29 %, 2.3σ** (bulk KR gave +3.59 %, 6.1σ).

### And the answer to the attribution question is NO

Does Z_box − Z_KR match Paper 1's +1 % in c_s through c_s² = Z + ηZ′ + Z²? With
d ln c_s²/d ln Z = 1.502, a +2.48 % excess in Z implies **+1.86 % in c_s**. Paper 1 measures
**+1.01 ± 0.27 %** at this η. They disagree by roughly 3σ.

This is the same disagreement Paper 1 already documented and withdrew: the wall force is a
**surface contact value**, and it over-states the bulk stiffness. At Paper 1's higher densities it
over-predicted the c_s offset five-fold; here at η = 0.10 it over-predicts it by about 1.8×. So the
measured-F(L) route removes the equation of state and the area convention, and substitutes a
different systematic in their place.

### What this does and does not settle

The two routes **bracket** the measurement rather than closing on it:

| model for the gas spring | s_qs | vs measured | σ |
|---|---|---|---|
| bulk Kolafa–Rottner (nominal area) | 0.8404 | +3.59 % | 6.1 |
| **measurement 0.8705 ± 0.0049** | — | — | — |
| measured wall force Z_box(η) | 0.8819 | −1.29 % | 2.3 |

The bulk EOS under-stiffens the box; the wall contact value over-stiffens it; the truth sits between
them, and the measurement does too. **Level 3's verdict is unchanged by this set** — it closes on
the observable, and the attribution stays open. What would close it is a *bulk* stiffness for this
box, which is what Paper 1's c_s already is — and the A2 finite-size ladder is the instrument for
that, not the wall.

## 6. Corrections carried from v6

Two of mine and two of the user's, all withdrawn in v6 or here:

- **Mine (v5):** the +1.35 % standing force explains "a third to a half" of the settled residual.
  Wrong — a constant offset cancels from the fixed point, and the excess is a contact value, not a
  bulk stiffening.
- **Mine (v4):** the residual is the acoustic pulse. Wrong — a no-push control reproduces it, and it
  does not scale with u.
- **User's:** ∂ln s/∂ln F = k/(k + k_gas) = 0.89. That is the sensitivity to a *multiplicative*
  scaling, which is the withdrawn Z_wall candidate in disguise; for an additive offset it is exactly 0.
- **User's:** the pre-push offset is (F_box − F_EOS)/k. The restoring stiffness is k + k_gas, giving
  0.0386 σ against a measured 0.0427 (M_s = 50) and 0.0578 (200). The check still stands as an
  independent confirmation that the 1.35 % wall-force excess is real and static.
