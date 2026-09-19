# Paper 2 — ramp against step, Paper 1's linewidth, and the fast end

2026-09-18. Reference: `260918_paper2_theory_audit_COWORK.pdf` §2.4, §2.5, §3, §4.2.
No git commands were run. The new protocol mode was first built separately as `00ALLINONE_ramp`; after
its gates passed it was promoted into the installed `hspist3/00ALLINONE` (sha cc3b9bfe → 1595c7fe, same
output byte for byte) — see appendix A2. Health contract: **0 aborts, 0 `[EDMD-HEALTH]` lines** over all 1080 new runs.

Scripts: `validation/paper2_geometry_fix_20260918.py`, `paper2_ramp_fast_20260918.py`,
`paper1_linewidth_20260918.py`, `level2_Au_figure_20260918.py`.
Figure: `260909_plots/260918_level2_A_of_u.png`.

---

## 0. Housekeeping — the Level 2 traces are compressed

Verified **before** compressing: the whole Level 2 analysis was re-run with every `tr_*.csv` hidden
from `glob`, and every number came back identical — A = 25.21, W(0) = 7.4728 ± 0.0089, ζ self term
2.615, the full momentum table. `stop_digest.csv` (85 kB) carries all of it.

| | before | after |
|---|---|---|
| `level2_slope_20260917/` total | 6.8 GB | **552 MB** |
| raw traces | 6.4 GB | 143 MB (gzip −6) |

Disk free went 17 → 19 GB. Nothing was deleted; `gunzip` restores the originals.

---

## 1. A correction that came out of building the ramp: the compartment is 0.5 σ shorter than we thought

While checking that the ramp reproduced the step's travel I found that the piston moves **4.18 σ**,
not the 3.93 σ the flag requests. Chasing it down produced a geometry correction that touches every
Paper 2 number so far. Three measured facts:

1. `wall_thickness_sigma = 1` in every energy-transfer run (it is in `summary.csv`), and the wall
   flag sets the divider **centre**, so the gas-side face is at 39.25 + 0.5 = **39.75**.
2. The stop snapshots confirm it: over 60 runs the closest disk centre to the divider is 40.2867,
   and a disk of radius 0.5 cannot approach a face at 39.75 closer than 40.25.
3. The gas's right boundary during the hold is the **box wall at 78.50** — the travel flag is
   measured from it (target = 78.50 − 3.93 = 74.57, which is where the piston stops). The piston
   parks at 78.75, so its first 0.25 σ of travel happens outside the gas.

So the compartment is **38.75 σ**, not the nominal L₀ = 39.25, and the gas is compressed by exactly
3.93 σ while the piston displaces 4.18 σ. `paper2_level0_level1_20260916.py` used the nominal length
twice — once to turn the wall impulse into Z (Z = Σ|dp|·L/(T·N_s)) and once to set η — so both move:

| true L | η (was) | η (true) | Z_wall (true) | Z_KR | ratio |
|---|---|---|---|---|---|
| 38.750000 | 0.100051 | 0.101342 | 1.2768 ± 0.0027 | 1.2399 | 1.0298 |
| 37.770833 | 0.102611 | 0.103969 | 1.2935 ± 0.0026 | 1.2472 | 1.0371 |
| 36.791667 | 0.105305 | 0.106736 | 1.3052 ± 0.0025 | 1.2551 | 1.0400 |
| 35.812500 | 0.108144 | 0.109654 | 1.3064 ± 0.0027 | 1.2634 | 1.0340 |
| 34.812500 | 0.111207 | 0.112804 | 1.3241 ± 0.0027 | 1.2725 | 1.0405 |

The finite-box excess over bulk KR drops from ~4.6 % to ~3.0–4.0 %, and

```
∫Z dlnη = 0.139215  →  T_f/T_i = 1.149371
W_qs^finite = 7.4685 ± 0.0131 kT      (was 7.4879 ± 0.0066)
```

**This fixes a sign problem nobody had noticed.** The published Level 1 result was
W(0) − W_qs = −0.0151 ± 0.0110, i.e. the measured work sat *below* the quasi-static bound — which an
adiabatic compression may not do. On the corrected geometry:

```
W(0) − W_qs^finite = +0.0043 ± 0.0158 kT   →   0.3 σ   →   PASS, and on the physical side
```

Everything downstream that used Δx = 3.93 as the piston displacement should use 3.93 for the **gas**
and 4.18 for the **piston**: ζ₀Δx = 10.28 (gas) or 10.93 (piston), against 10.01 as published.

---

## 2. Item 1 — ramp against step: the u² term is mostly protocol-launched

New protocol mode `--piston-right-protocol-mode=ramp --piston-ramp-time=T`: velocity rises linearly
0 → u, holds, then falls linearly to zero, arriving at the **same target with the same travel**. The
deceleration is keyed to the remaining distance (v = √(2ad)), which is the linear-in-time ramp
written without a clock, so the piston lands on the target with v = 0 whatever the step size. T is
clamped at travel/u, where the profile becomes triangular; measured pushes were 123.5 / 81.7 / 41.7 σ
against 123.6 / 81.8 / 41.8 predicted.

**Gates, both passed before any batch.** Step mode byte-identical to the installed binary on 3 seeds
(traces *and* event logs, 9.36 MB and 583 kB each); ramp-mode ledgers on 5 seeds
|R_E|/W = 7.22 × 10⁻¹², |R_p|/Σ|dp| = 2.57 × 10⁻¹⁴, matching the step control (7.02 × 10⁻¹², 2.10 × 10⁻¹⁴).

T = 40 σ-time ≈ 2L/c_s, 60 seeds per cell, W_qs = 7.4685:

| protocol | u | ⟨W⟩ [kT] | W − W_qs [kT] | (W−W_qs)/u² |
|---|---|---|---|---|
| step | 0.05 | 7.4620 ± 0.0471 | −0.0065 ± 0.0488 | −2.62 ± 19.54 |
| step | 0.10 | 7.7332 ± 0.0777 | +0.2646 ± 0.0788 | 26.46 ± 7.88 |
| step | 0.20 | 8.5783 ± 0.2108 | +1.1098 ± 0.2113 | 27.75 ± 5.28 |
| ramp | 0.05 | 7.5247 ± 0.0342 | +0.0561 ± 0.0366 | 22.45 ± 14.64 |
| ramp | 0.10 | 7.5654 ± 0.0477 | +0.0969 ± 0.0495 | 9.69 ± 4.95 |
| ramp | 0.20 | 7.9248 ± 0.1141 | +0.4563 ± 0.1149 | 11.41 ± 2.87 |

```
A_step = 25.91 ± 4.28   (χ²/dof = 1.13)     — reproduces the 760-run campaign value 25.2
A_ramp = 11.30 ± 2.45   (χ²/dof = 0.34)
difference 14.61 ± 4.93  (3.0 σ)
A_ramp against N_s m/6 = 8.3 : +1.2 σ        A_ramp against N_s m/2 = 25 : −5.6 σ
```

**This is the decision the audit asked for.** A_ramp is neither ≈ A_step nor ≈ 0: it lands at 11.3,
which is 5.6 σ below the step value and 1.2 σ from the steady-flow prediction N_s m u²/6. So roughly
**14 of the 25 was energy the step protocol launched** as a pressure pulse, and what a gentle push
leaves behind is consistent with the flow term alone.

**One caveat on how gentle the ramp really is.** Because the piston's first 0.25 σ happens outside
the gas, the gas does not see the ramp start from rest: at contact the piston is already at
½a·t²  = 0.25 σ, i.e. v = 0.025 (u = 0.05), 0.035 (u = 0.10), 0.069 (u = 0.20) — 35–50 % of u. A
ramp that begins at the wall would be gentler still, so **11.3 is an upper bound** on what the
protocol-independent excess is. That is one flag change, not new code.

### The decomposition at stop

From the env-gated stop snapshots (`HD_STOP_SNAPSHOT`), 10 bins across the compressed compartment,
seed-averaged. Two corrections were essential and are worth recording: the naive estimator charges
the flow term (½ n_bins kT/n_seeds = 0.08 kT) and the compression term (c²n_bins/2n_seeds = 0.31 kT)
with pure sampling noise — which is the entire signal at u ≤ 0.1 — and the equilibrium state is not
uniform (the disks layer against the walls), so an equilibrium profile at the same L must be
subtracted. That baseline is 60 runs at u = 0.005, whose excess work is zero. Errors are bootstrap
over seeds.

| protocol | u | excess W−W_qs | flow ½mΣN_b v̄_b² | compression | rest |
|---|---|---|---|---|---|
| step | 0.05 | −0.0065 ± 0.0488 | 0.025 ± 0.080 | −0.044 ± 0.185 | +0.012 |
| step | 0.10 | +0.2646 ± 0.0788 | 0.038 ± 0.075 | −0.065 ± 0.185 | +0.291 |
| step | 0.20 | +1.1098 ± 0.2113 | **0.529 ± 0.157** | −0.256 ± 0.145 | +0.837 |
| ramp | 0.05 | +0.0561 ± 0.0366 | 0.141 ± 0.107 | −0.351 ± 0.134 | +0.267 |
| ramp | 0.10 | +0.0969 ± 0.0495 | 0.140 ± 0.104 | −0.214 ± 0.159 | +0.171 |
| ramp | 0.20 | +0.4563 ± 0.1149 | 0.013 ± 0.068 | +0.051 ± 0.210 | +0.392 |

**Only one row is resolved:** step at u = 0.20, flow = 0.529 ± 0.157 (3.4 σ), against the predicted
N_s m u²/6 = 0.333 (1.2 σ high). Every other entry is consistent with zero at 60 seeds; the
compression column is negative as often as positive, which is what a noise-dominated estimator looks
like. **The decomposition is not yet a measurement** — it needs a few hundred seeds per cell, or a
finer instrument than 10 bins of 50 particles. The conclusion of item 1 rests on A_ramp, not on this
table.

---

## 3. Item 2 — Paper 1's linewidth gives Γ, and Γ agrees with Enskog

FWHM of the seed-averaged position spectrum, 200 periods (one bin = 0.50 % of f), interpolated on
the half-maximum, 25 seeds per mass, A1 v2 data — **no new runs**. Mansour's piston is our divider:
Q = 1/μ with μ = Γ L_y √(2/(M̂N)), M̂ = M + mN/3. The nearest A1 density to 0.10 is η = 0.1122.

| M | α = M/N | M̂ | measured Δf/f [%] | 1/Q with Γ = 0.331 [%] | ratio | Γ implied |
|---|---|---|---|---|---|---|
| 50 | 0.5 | 83.3 | 1.39 | 5.13 | 0.27 | 0.090 |
| 100 | 1 | 133.3 | 4.17 | 4.05 | 1.03 | 0.341 |
| 200 | 2 | 233.3 | 3.42 | 3.06 | 1.11 | 0.369 |
| 300 | 3 | 333.3 | 2.93 | 2.56 | 1.14 | 0.378 |
| 500 | 5 | 533.3 | 1.30 | 2.03 | 0.64 | 0.212 |
| 750 | 7.5 | 783.3 | 1.51 | 1.67 | 0.90 | 0.299 |
| 1000 | 10 | 1033.3 | 1.49 | 1.46 | 1.02 | 0.339 |
| 1500 | 15 | 1533.3 | 1.64 | 1.20 | 1.37 | 0.454 |
| 2000 | 20 | 2033.3 | 1.47 | 1.04 | 1.41 | 0.467 |

```
Γ implied, all nine masses      : 0.327 ± 0.039   (spread 36 %)
Γ implied, α > 2 only           : 0.358 ± 0.039   (spread 27 %)
Enskog, audit Eq. 5 at η = 0.100: 0.331
```

**The value lands on Enskog to 1 %** (0.327 ± 0.039 against 0.331, 0.1 σ). The M = 50 row is the
outlier the audit predicted — α = 0.5 is exactly where K ≪ 1 fails, and it gives Γ = 0.09, a factor
3.7 low; it is outside the model and is reported, not used.

**But the scaling test is only marginally passed.** Γ must be mass-independent, and it wanders from
0.21 (M = 500) to 0.47 (M = 2000), a 27–36 % spread with no clean trend. At 200 periods a 1.5 % line
is three bins wide, so the heavy-mass rows are resolution-limited: that is the likely source, and
longer records would settle it.

Repeating at the next two densities (also asked for):

| η | Γ, all masses | Γ, α > 2 | spread (α > 2) |
|---|---|---|---|
| 0.1122 | 0.327 ± 0.039 | 0.358 ± 0.039 | 27 % |
| 0.2618 | 0.419 ± 0.039 | 0.485 ± 0.022 | 11 % |
| 0.5236 | 0.687 ± 0.053 | 0.753 ± 0.031 | 10 % |

Γ rises with density as Enskog requires, and the spread *shrinks* at higher density (the lines are
wider there, so the 0.5 % bin matters less) — consistent with resolution being the limit at η = 0.11.

**What this buys Paper 2, with no new runs.** ζ_hyd = L_y Γ/X_p = 10 × 0.327/38.75 = **0.0844 ±
0.0101**, so the predicted linear term is

```
dW/du = ζ_hyd Δx = 0.332 ± 0.040      against the measured −0.50 ± 0.67 for u ≤ 0.05  (1.2 σ)
```

Consistent, and still below the resolution of 760 trajectories — exactly as the audit said. The
Level 2 sentence should now read "the Markovian friction is the viscous one, it is 0.33 u, and it is
below our resolution", with Γ sourced from Paper 1 rather than assumed.

---

## 4. Item 3 — the fast end

u = 3, 5, 10, step protocol, 100 seeds each. Two travels: 1.00 σ (as specified) and 7.08 σ (the
first batch, kept because it tests the Δx-linearity of Eq. 7 for free). W_qs is recomputed for each
travel on the corrected geometry.

| Δx | u | ⟨W⟩ [kT] | Eq. 7: 2N_s m(Δx/L)u² | ratio | Var W | Eq. 7 Var | ratio |
|---|---|---|---|---|---|---|---|
| 1.00 σ | 3 | 21.90 ± 1.91 | 23.23 | 0.943 | 363 | 465 | 0.78 |
| 1.00 σ | 5 | 57.15 ± 5.05 | 64.52 | 0.886 | 2546 | 3355 | 0.76 |
| 1.00 σ | 10 | 209.18 ± 18.29 | 258.06 | 0.811 | 33469 | 52129 | 0.64 |
| 7.08 σ | 3 | 183.92 ± 5.28 | 161.42 | 1.139 | 2788 | 3228 | 0.86 |
| 7.08 σ | 5 | 469.97 ± 12.26 | 451.61 | 1.041 | 15023 | 23484 | 0.64 |
| 7.08 σ | 10 | **1802.29 ± 47.80** | **1827.96** | **0.986** | 228441 | 369247 | 0.62 |

**Eq. 7 is confirmed to 1.4 % at u = 10 with Δx = 7.08 σ**, with no free parameter, and the approach
is monotonic in u (1.14 → 1.04 → 0.99) exactly as a single-hit limit should be.

> **Withdrawn — see appendix A3.** The Δx = 1.00 σ row above is wrong: the travel flag is the *gas*
> compression, while the piston displaces 0.25 σ more, so that set compressed 0.75 σ and its Eq. 7 is
> 33 % too large here. Corrected, there is no deficit at small travel but an ~8 % **excess** that
> decays monotonically with travel and is **not explained**. The "strip ambiguity" argument that
> stood here is withdrawn: the centres a moving face can strike lie in a strip exactly Δx wide, so
> there is no ±50 % ambiguity in N_hit.

The variance sits at 0.62–0.86 of Eq. 7 throughout. The prediction treats N_hit as Poisson and the
hits as independent; in a strip that is already 18 % of the compartment the particles are neither,
so a deficit is expected. Worth one line in the paper, not a new campaign.

### The figure

`260918_level2_A_of_u.png` — A(u) = (W − W_qs)/u², every point scaled to a common gas travel of
3.93 σ, from u = 0.005 to u = 10 on a log axis, with the two predicted plateaus drawn:

```
acoustic (step, measured)      A = 25.2      reached at u = 0.1-0.2
single-hit (audit Eq. 7)       A = 10.1      reached at u = 3-10
steady flow (audit Eq. 15/18)  A =  8.3      where the RAMP sits:   9.7, 11.4
```

**Rebuilt 2026-09-18 evening** (appendix A3): the first version scaled the fast sets by the *piston*
displacement, so their points sat ~25 % low. Every cell's gas travel is now measured from its own
trace rather than assumed, and the figure has a second panel carrying the new travel scan —
⟨W⟩/Eq. 7 against gas Δx at u = 3, 5, 10. Corrected values of A (scaled to Δx = 3.93 σ):

| set | u = 3 | u = 5 | u = 10 |
|---|---|---|---|
| gas Δx = 0.75 σ | 12.0 ± 1.1 | 11.7 ± 1.1 | 10.9 ± 1.0 |
| gas Δx = 2.1–5.1 σ (new) | — | — | 11.0 → 10.1 |
| gas Δx = 7.0 σ | 10.6 ± 0.3 | 10.2 ± 0.3 | 9.9 ± 0.3 |

against the parameter-free 10.1 — so the fast plateau is now hit from above and converges onto the
line, instead of sitting below it as the stale figure showed.

The crossover the audit asked for is there, and the ramp points sit on the lower plateau while the
step points sit on the upper one at the *same* speeds — which is the whole result of item 1 in one
picture.

---

## Verdict

| item | result |
|---|---|
| 0 housekeeping | traces 6.8 GB → 552 MB, digest verified first — **done** |
| geometry | compartment 38.75 σ not 39.25; W_qs = 7.4685 ± 0.0131; Level 1 now +0.0043 ± 0.0158 (0.3 σ) and on the physical side — **corrected** |
| 1 ramp vs step | A_step = 25.9 ± 4.3, A_ramp = 11.3 ± 2.5, 3.0 σ apart; the u² term is mostly protocol-launched — **decided** |
| 1b decomposition | only step u = 0.2 flow is resolved (0.53 ± 0.16 vs 0.33 predicted); needs ~300 seeds — **inconclusive** |
| 2 linewidth | Γ = 0.327 ± 0.039 vs Enskog 0.331; ζ_hyd = 0.0844, slope 0.332 ± 0.040 — **holds, value; scaling marginal** |
| 3 fast end | Eq. 7 confirmed to 1.4 % at u = 10, Δx = 7.08 σ; variance 0.62–0.86 of prediction — **passed** |

## Queued, not started

Level 2B, the long or absorbing compartment (L = 4× at the same η, 25 seeds, u = 0.05 and 0.1),
which asks whether the recurrence disappears and ζ climbs from ζ_hyd toward ζ₀. Waiting for KOA.

---

# Appendix, 2026-09-18 afternoon — five follow-ups

## A1. Paper 1's wall thickness: **0.05 σ, proved by reproduction**

The worry was specific: if the speed-of-sound campaigns had also run at the default thickness of 1 σ,
then L_eff = L₀ − 2r would be wrong by 0.5 σ and every c_s would be high by 0.5/(L₀ − 1) — 1.3 % at
η = 0.10, 4 % at 0.30, 8 % at 0.50, shrinking with N exactly like the finite-size offset we have been
attributing to physics.

Inference was not enough, so this was settled by rerunning a stored trajectory:

| campaign | evidence | thickness |
|---|---|---|
| A1 v2 (`A1v2_20260914`) | `tests_20260913.py:233` passes it, **and** the stored trace reproduces | **0.05** |
| A2 famB (`famB_20260911`) | `--wall-thickness=0.05` in the recorded command | 0.05 |
| A2 α = 2 (`A2_alpha2_20260912`) | same | 0.05 |
| A2 top-up (`A2_topup_20260912`) | same | 0.05 |
| A2 dilute 50-period (`A2_dilute50_20260917`) | same | 0.05 |

The reproduction is the strong one. Taking the seed out of the stored A1 v2 trace
(η = 0.5236, M = 1000, run0, seed 1200845163) and re-running it through `--speed-sound-exact-seed`:

```
--wall-thickness=0.05 :  6420/6420 shared sample times identical, max |ΔWall_X| = 0.000e+00
--wall-thickness=1.0  :  1/6420 identical,                        max |ΔWall_X| = 5.357e-01
```

**Paper 1 is safe.** No re-analysis is needed and the finite-size section stands.

**Methods sentence to add:** *divider thickness 0.05 σ; L_eff = L₀ − 2r − t/2.* The disk centres in one
compartment span from r to L₀ − t/2 − r, so the thickness takes another t/2 = 0.025 σ off L_eff. Since
c_s is proportional to L_eff, this lowers every measured c_s by 0.025/(L₀ − 1):

| η | L₀ | shift in c_s |
|---|---|---|
| 0.10 | 39.27 | −0.065 % |
| 0.30 | 13.09 | −0.207 % |
| 0.50 | 7.85 | **−0.365 %** |
| 0.5236 | 7.50 | −0.385 % |
| 0.65 | 6.04 | −0.496 % |

One-signed and below the error bar everywhere except the densest points, where it is comparable to it
— and it moves those points *down*, i.e. slightly toward Kolafa–Rottner.

**The corrected figure exists**: `260909_plots/260918_cs_vs_eta_thickness_corrected.png`, published
and corrected points overlaid with a deviation panel, written by
`validation/paper1_thickness_correction_20260918.py`. The correction is exact and multiplicative —
c_s ∝ L_eff, so c_s → c_s·(L₀ − 1 − 0.025)/(L₀ − 1) — and no trajectory is re-analysed.

| η | L₀ | c_s published | corrected | dev. from KR before | after |
|---|---|---|---|---|---|
| 0.006545 | 600.00 | 1.4430 | 1.4430 | +0.71 % | +0.70 % |
| 0.078540 | 50.00 | 1.6795 | 1.6787 | +0.90 % | +0.85 % |
| 0.196350 | 20.00 | 2.2139 | 2.2110 | +1.56 % | +1.43 % |
| 0.392699 | 10.00 | 3.8194 | 3.8088 | +1.96 % | +1.68 % |
| 0.523599 | 7.50 | 6.1117 | 6.0882 | +3.03 % | +2.63 % |
| 0.650003 | 6.04 | 11.0858 | 11.0308 | +6.20 % | +5.67 % |

**And it does not explain the finite-size offset**, which was the reason to check in the first place.
The correction is larger for smaller boxes, so it shrinks the N = 100 − N = 1600 gap, but only a
little:

| η | gap before | gap after |
|---|---|---|
| 0.10 | +0.84 % | +0.78 % |
| 0.30 | +1.01 % | +0.88 % |
| 0.50 | +2.24 % | +1.93 % |
| 0.60 | +0.94 % | +0.59 % |
| 0.65 | +7.25 % | +6.83 % |

At η = 0.50 it accounts for 0.31 % of a 2.24 % gap — about one part in seven. The finite-size
reading of the A1–A2 offset survives. Had the thickness been 1 σ the shift would have been 20×
larger and would have swallowed the whole effect, which is exactly why this was worth proving rather
than assuming.

## A2. The patched source is promoted

```
old installed binary : cc3b9bfe625d9ad05f2f38596a4bf5858cec33f1  (2026-09-14, backed up as
                       _orchestration_tests_20260913/binaries/00ALLINONE_sha_cc3b9bfe_20260914)
new installed binary : 1595c7fe21e58b1fd6e75758985a358d4c2b8a01  (2026-09-18, from the patched source)
gate, 3 seeds        : traces and event logs IDENTICAL to the old binary (9.36 MB and 583 kB each)
```

`00ALLINONE_ramp` is now a symlink to `00ALLINONE`, so every script keeps working and the ramp mode
and stop snapshot are production features. Nothing was deleted.

## A3. The fast end: my explanation was wrong **and so was my number**

The "strip ambiguity" paragraph is withdrawn, and the deficit it was invented to explain does not
exist. **The travel flag is the gas compression, not the piston displacement** — the piston parks
0.25 σ outside the box wall, so a flag of 0.75 gives a piston displacement of 1.00 and a gas
compression of 0.75. The 2026-09-18 table used 1.00 for that set, making its Eq. 7 prediction 33 %
too large and manufacturing the deficit.

Corrected, with the new travel scan at u = 10 (100 seeds each, 0 aborts, 0 health lines):

| gas Δx [σ] | u | ⟨W⟩ [kT] | Eq. 7 | ratio |
|---|---|---|---|---|
| 0.750 | 3 | 21.90 ± 1.91 | 17.42 | 1.257 |
| 0.750 | 5 | 57.15 ± 5.05 | 48.39 | 1.181 |
| 0.750 | 10 | 209.18 ± 18.29 | 193.55 | 1.081 |
| 2.083 | 10 | 585.14 ± 23.39 | 537.63 | 1.088 |
| 3.083 | 10 | 837.67 ± 29.90 | 795.70 | 1.053 |
| 5.083 | 10 | 1322.92 ± 38.73 | 1311.83 | 1.008 |
| 7.083 | 10 | 1802.29 ± 47.80 | 1827.96 | 0.986 |
| 6.950 | 3 | 183.92 ± 5.28 | 161.42 | 1.139 |
| 7.000 | 5 | 469.97 ± 12.26 | 451.61 | 1.041 |

At fixed u = 10 the ratio falls monotonically with travel, 1.081 → 1.088 → 1.053 → 1.008 → 0.986,
crossing unity near Δx ≈ 5 σ. So there is no deficit at small travel: there is an **excess of about
8 %** that decays as the swept strip lengthens. **This is not explained.** A wall-layer effect is the
obvious candidate — the equilibrium density within the first σ of the piston is enhanced by the
contact value g(η) = 1.20 at this density, so a short push sweeps more particles than N_s Δx/L
predicts, and a long one averages the enhancement away — but nothing here tests that, and it is
written down as a candidate, not a result.

Eq. 7 itself is confirmed: 0.986 at the longest travel and highest speed, with no free parameter.

## A4. Linewidth, with the bin taken out in quadrature

Δf_true = √(Δf² − Δf_bin²), Δf_bin/f = 1/200 by construction:

| M | α | raw Δf/f [%] | deconvolved [%] | 1/Q at Γ = 0.331 [%] | ratio | Γ implied |
|---|---|---|---|---|---|---|
| 50 | 0.5 | 1.39 | 1.30 | 5.13 | 0.25 | 0.084 |
| 100 | 1 | 4.17 | 4.14 | 4.05 | 1.02 | 0.338 |
| 200 | 2 | 3.42 | 3.38 | 3.06 | 1.10 | 0.365 |
| 300 | 3 | 2.93 | 2.88 | 2.56 | 1.12 | 0.372 |
| 500 | 5 | 1.30 | **1.20** | 2.03 | **0.59** | **0.195** |
| 750 | 7.5 | 1.51 | 1.42 | 1.67 | 0.85 | 0.282 |
| 1000 | 10 | 1.49 | 1.40 | 1.46 | 0.96 | 0.319 |
| 1500 | 15 | 1.64 | 1.56 | 1.20 | 1.31 | 0.432 |
| 2000 | 20 | 1.47 | 1.38 | 1.04 | 1.33 | 0.439 |

| η | Γ, all masses | Γ, α > 2 | spread (α > 2) |
|---|---|---|---|
| 0.1122 | **0.314 ± 0.038** | 0.340 ± 0.038 | 28 % |
| 0.2618 | 0.409 ± 0.038 | 0.472 ± 0.021 | 11 % |
| 0.5236 | 0.681 ± 0.053 | 0.744 ± 0.033 | 11 % |

**The claim is now: Γ = 0.31 ± 0.04 at η = 0.112, consistent with Enskog's 0.331 within its 12 %
error (0.4 σ).** Not "to 1 %" — the central value landing on 0.331 was partly luck, and the 36 %
spread across masses is the honest uncertainty on the scaling.

**Deconvolution does not rescue M = 500.** Resolution can only broaden a line, so it explains the
heavy-mass rows reading high (1.38–1.56 % against 1.04–1.20 predicted) but it makes M = 500 *worse*:
1.30 → 1.20 % against 2.03 predicted, ratio 0.59. That row is flagged as unexplained.

Downstream: ζ_hyd = 0.081 ± 0.010, predicted linear term 0.32 ± 0.04 u (was 0.33 ± 0.04).

## A5. Plan document

Added to the SDL section: the **two-wall spring chain** (geometry D, `--num-walls=2
--wall-positions=20,30 --wall-mass-factors=1000,200 --spring-k=5`, three compartments of 50 disks),
verified to start and initialise the spring measurement — because Paper 2's figure 1 is geometry D,
not geometry A. The Level 2 results box now carries A_ramp, Γ and Eq. 7 at 7.08 σ on the corrected
geometry. PDF is 12 pages, 0 errors.

The A(u) figure was rebuilt on the corrected travels and now carries the travel scan as a second
panel; see the amended "The figure" block in section 4.

---

# Appendix B, 2026-09-18 — making the L_eff error impossible, and what else could move N = 100

## B1. The error is now structural, not a habit

The L_eff mistake had one cause: **the launcher and the estimator each carried their own copy of the
geometry, and nothing compared them.** `argv_for()` passed `--wall-thickness=0.05` while `x_of()`
computed `L0 - 2r`. Both were internally consistent; together they were wrong, and no file on disk
said so. Four changes, all in `validation/`:

1. **One constant.** `tests_20260913.py` now defines `WALL_T = 0.05` once. `argv_for()` passes it to
   the binary and `l_eff(L0) = L0 - 2r - t/2` consumes it. The launcher and the estimator cannot
   disagree any more; changing the constant changes the runs and the analysis together.
2. **One definition of L_eff.** `A2_dilute_20260916.py`, `A2_dilute50_20260917.py`,
   `A2_long200_20260915.py` and `estimator_massladder_20260917.py` each had their own `(L0 - 1.0)`;
   all four now call `T.l_eff(L0)`.
3. **A guard that refuses bad data.** `T.assert_wall_thickness(run_dir)` reads the thickness back out
   of a campaign's own command record and raises on a mismatch. It returns `None` — reported, never
   assumed to agree — when a campaign recorded nothing.
4. **The launcher now records its command.** The speed-of-sound harness wrote only a `run.log`, which
   does not carry the argv; that is precisely why the geometry could not be read back from the data
   and the error stayed invisible. It now writes `00_COMMAND.md` beside each cell, including the line
   `Estimator geometry: L_eff = L0 - 2r - t/2 with r = 0.5, t = 0.05`, so the guard has something to
   check on every future campaign.

Verified: `l_eff(7.5) = 6.475` (was 6.5), the guard reads 0.05 back from `famB_20260911`, returns
`None` for `A1v2_20260914` (which predates the change), and all four dependent scripts import and run.

**Consequence to decide on:** any figure regenerated from these scripts now carries the correction.
The three circulated figures are stale until regenerated, and the numbers move by the amounts in
appendix A1 (nothing below η = 0.3, −0.4 % at η = 0.52, −0.5 % at 0.65).

## B2. What else could move N = 100 toward Kolafa–Rottner?

Two further candidates were examined. **One goes the wrong way and one is not an error at all.**

**The gas's own inertia loading the divider.** Mansour's piston carries M̂ = M + mN/3, not the bare
wall mass. If Román's α should be built from M̂, then α rises, K(α) falls, and c_s *rises*:

| M | α bare | α with M̂ | c_s would change by |
|---|---|---|---|
| 50 | 0.500 | 0.833 | **+17.3 %** |
| 200 | 2.000 | 2.333 | +6.9 % |
| 1000 | 10.000 | 10.333 | +1.6 % |
| 2000 | 20.000 | 20.333 | +0.8 % |

One-signed **upward** — it would push N = 100 further above KR, not closer. It is also not a
finite-size effect: α shifts by mN/(2N_s) = 1/3 at every N. And the mass-ladder test of 2026-09-17
is evidence against it: with the bare mass the nine masses already lie on one line to 0.27 %, which a
17 % mass-dependent distortion at M = 50 would destroy. Not applied; recorded as checked.

**The finite box is genuinely stiffer — and Paper 2 measured it.** Level 1 measured the wall pressure
of the *same* geometry (N_s = 50) directly: Z_wall/Z_KR = 1.030 at η = 0.101 rising to 1.040 at
η = 0.113, on 100 seeds per point. Feeding a 3 % enhancement of Z through c_s² = T(Z + ηZ′ + Z²):

| η | predicted from Z × 1.03 | measured A1 offset |
|---|---|---|
| 0.1122 | +2.25 % | +1.17 % |
| 0.2618 | +2.27 % | +1.75 % |
| 0.5236 | +2.37 % | +3.03 % |

Same sign, same order, and it vanishes as N grows — which is exactly what the A2 ladder shows. The
comparison is indicative rather than exact, because Z_wall carries a surface contribution (the disks
layer against the wall) that is not the bulk compressibility entering c_s. But it means **the N = 100
offset is most likely finite-box thermodynamics that we have independently measured, not an
uncorrected analysis error** — and Paper 1 can say so with a Paper 2 number behind it.

**Still unexamined**, for honesty: any aspect-ratio dependence at H = 10 (Wu et al. 2016 report wall
layering in narrow channels). *(The ensemble item that stood here — microcanonical versus canonical T
at N = 100 — was struck on 2026-09-19: drift-first rescales to KE = N_s kT exactly after removing the
drift, so T = 1 by construction and there is no 1/N convention error to chase.)*

---

# Appendix C, 2026-09-19 — canonical figures, the dense-geometry Z, the fast end explained, demo mode

## C1. The thickness correction is canonical, verified against the estimator

`validation/paper1_canonical_20260919.py`. One cell was recomputed from raw traces with the production
estimator (η = 0.5236, M = 1000, 25 seeds, 0 discarded) and compared with the closed-form rescale:

```
c_s from traces, L_eff = L0 - 2r - t/2 : 6.086728
c_s from traces, L_eff = L0 - 2r       : 6.110229
ratio                                   : 0.99615385
(L0 - 1 - 0.025)/(L0 - 1)               : 0.99615385      agreement 1.1e-16
```

So the rescale *is* the code path, and 7875 trajectories did not need re-reading. New canonical files:
`260919_A1v2_final_cs_vs_eta.csv` (deviation column recomputed), `260919_A2_cs_per_mass.csv`,
`260919_A2_cs_per_mass_famB.csv`, and three figures — `260919_cs_vs_eta`,
`260919_cs_vs_eta_lowdensity_zoom`, `260919_cs_vs_eta_N100_vs_A2`. The 260914/260916/260917 versions
are untouched. Both figure scripts now take their CSVs from `HD_A1_CSV` / `HD_A2_CSV`, because the
first attempt silently drew a corrected A2 against a stale hardcoded A1.

## C2. Z of the finite box at Paper 1's own densities — and it does **not** explain the offset

120 hold-only runs (40 per density) in Paper 1's geometry, t = 0.05, N_s = 50, 0 aborts, 0 health lines:

| L₀ | true compartment | η | Z_wall | Z_KR | Z_wall/Z_KR | c_s excess that implies | measured A1 offset |
|---|---|---|---|---|---|---|---|
| 20.0 | 19.975 | 0.1966 | 1.6647 ± 0.0031 | 1.5568 | 1.069 | **+5.22 %** | +1.56 % |
| 10.0 | 9.975 | 0.3937 | 3.1741 ± 0.0048 | 2.7693 | 1.146 | **+11.23 %** | +1.96 % |
| 7.5 | 7.475 | 0.5254 | 5.4838 ± 0.0060 | 4.5651 | 1.201 | **+15.86 %** | +3.03 % |

**This kills the candidate as a quantitative explanation, and the appendix B version of it was too
generous.** At η = 0.10 the wall pressure was 3 % above bulk and predicted +2.25 % against a measured
+1.17 % — close enough to look like agreement. Measured at the densities that matter it is 20 % above
bulk and predicts +15.9 % against +3.0 %, over-predicting by a factor of five, and the discrepancy
grows with η.

The reason is that Z_wall is a *surface* quantity: it is the contact value at the wall, inflated by the
density layer that builds against a hard boundary, while the sound mode samples the bulk. So the
finite box is stiffer at its walls than in its interior, and only the interior sets c_s. **The N = 100
offset remains finite-size in origin but is not quantified by Z_wall.** The instrument that does
measure it correctly is the one Paper 1 already uses — the A2 ladder and its c_∞ extrapolation.

## C3. The fast-end excess is explained: it is the contact layer, counted properly

Appendix A3 left an ~8 % excess at short travel "not explained". It is now explained, and the
explanation was measurable rather than assertable. Eq. 7 assumes the swept strip holds n_bulk·Δx·H
centres. It does not: a disk centre cannot approach the face closer than r, and just inside that
exclusion the density *overshoots*. Counting centres in the strip that a face can actually strike —
width Δx, starting at face − r — over the 60 equilibrium snapshots at the same state:

| Δx [σ] | ⟨N in strip⟩ | n_bulk·Δx·H | ratio | ⟨W⟩/Eq. 7 (u = 10) | Eq. 7 using the **measured** count |
|---|---|---|---|---|---|
| 0.750 | 1.217 | 1.077 | **1.130** | 1.081 | 0.957 |
| 2.083 | 3.200 | 2.991 | 1.070 | 1.088 | 1.017 |
| 3.083 | 4.517 | 4.427 | 1.020 | 1.053 | 1.032 |
| 5.083 | 7.483 | 7.299 | 1.025 | 1.008 | 0.983 |
| 7.083 | 10.483 | 10.171 | 1.031 | 0.986 | 0.957 |

The strip holds 13 % more centres than uniform density at Δx = 0.75 σ and 3 % more at 7 σ — the contact
layer, diluted as the strip lengthens. Feeding the measured count into Eq. 7 removes the trend: the
ratio becomes 0.96–1.03 scattered about unity, instead of falling monotonically 1.08 → 0.99. **Eq. 7
is right; the uniform-density N_hit was the approximation.** A caveat worth keeping: measuring the
strip from the face rather than from face − r gives 0.39 × bulk and the opposite conclusion, so the
r-shift is the whole of the argument.

## C4. Demo mode

New GUI-only flag `--demo` (plus `--demo-shot=PATH[,STEPS]` for unattended capture).

**Gate first:** headless output byte-identical to the installed binary on 3 seeds, traces *and* event
logs. `--demo` changes no physics; it does not exist outside the render loop.

**What it does.** Starts paused, divider held, piston parked. `SPACE` run/pause, `R` release the
divider, `P` start the piston, `+`/`-` double/halve the pace, `S` screenshot, `Q`/`ESC` quit.
`--auto-release-after-hold` and `--auto-piston-step` are ignored so every stage begins on a key. A
protocol summary prints once (geometry, masses, u, travel, hold, spring k) and a HUD line sits at the
top of the window and is echoed to the terminal on every stage change:

```
[DEMO] PISTON step u=0.020 | t=56.8 sigma | piston x=78.35 | divider x=26.17 | W_in=0.390 |
       KE_L=33.92 KE_R=32.87 | E_spring=0.269
```

**Two things were wrong in the plan's SDL section, and this is why the window "flashed and finished".**
First, `--experiment=energy_transfer` is a self-contained driver that **never renders** — with or
without `--headless`. The window belongs to the interactive loop, reached with `--show-simulation`.
Second, the generic GUI advances about one step per 25 frames, so a 12 000-step push would have taken
over an hour; `--demo` sets two steps per frame, which puts the u = 0.02 push at roughly two minutes.
Both commands in the plan are rewritten, verified to run, and each carries a two-line "what you will
see" and a screenshot: `260909_plots/260919_demo_geomA.png` and `_geomD.png`. Plan PDF is 13 pages.

## C5. Struck

The microcanonical-versus-canonical temperature item is removed from appendix B's open list:
drift-first removes the drift and *then* rescales to KE = N_s kT, so T = 1 by construction.
