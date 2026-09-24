# The R-collapse test — R is the variable, but my prediction was dimensionally wrong

2026-10-05. `level4_Rcollapse_20261005`, **160 runs, 0 aborts, 0 health events**, 21:04–23:43.
First use of the N_s = 100 box: N = 200 particles, L_c = 77.5, box 156.0, divider centre 78.0,
recorded t = 1.0, η = 0.10134170 — identical to the ladder's to 8 digits, and ρ₀ = N_s/L_c = 1.2903,
the same. That shared ρ₀ is the point: it is the path Cencini's limit is defined along.

> **The pre-registered verdict rule returns NOT RESOLVED**, on the sampling criterion — L/τ came out
> 45 and 51 against the required 60, because τ is larger than *every* prediction I sized the records
> from.
>
> **But the reason it is larger is that prediction A was mis-stated, and correcting it selects R.**

---

## 1. What was run, and what came out

| M | R | L/τ | cal slope | modelled | block | S(0) | spread | g = τ/M |
|---|---|---|---|---|---|---|---|---|
| 50 | 2 | **45** | 0.605 | 7253 ± 644 | 7251 | 8874 | 1.22× | **145.1 ± 12.9** |
| 100 | 1 | **51** | 0.643 | 15 773 ± 915 | 17 804 | 17 660 | 1.13× | **157.7 ± 9.2** |

Slope ≥ 0.6 **passes** at both masses. L/τ ≥ 60 **fails** at both. Against the predictions as I wrote
them — A: g = 59.2, 76.7; B: g = 76.7, 123.3 — the measurement is **above both**, by 3.8σ to 8.9σ.
A test where the answer beats every hypothesis is usually a sign the hypotheses were built wrong.

## 2. The error: g = τ/M is not the invariant when L_c changes

τ_GP = (4/√(2π))·M·L_c/√(mk T)/(1+ηZ′/Z) is **proportional to M·L_c**. Doubling N_s at fixed η
doubles L_c, so τ_GP per unit mass doubles: 50.349 M → **100.70 M**. Therefore *even if the physics
collapses perfectly onto R*, g = τ/M must double at fixed R. I wrote prediction A as
"g(R) carries over unchanged", which builds the L_c scaling out of the hypothesis it was meant to
test.

**The correct statement of "the variable is R" is τ_T/τ_GP = f(R).**

## 3. In that variable, the datasets collapse

| dataset | N_s | L_c | M | R | τ_T | τ_GP | **f = τ_T/τ_GP** |
|---|---|---|---|---|---|---|---|
| ladder | 50 | 38.75 | 10 | 5.00 | 475 ± 20 | 503 | **0.94 ± 0.04** |
| ladder | 50 | 38.75 | 20 | 2.50 | 1089 ± 57 | 1007 | **1.08 ± 0.06** |
| ladder | 50 | 38.75 | 50 | 1.00 | 3835 ± 95 | 2517 | **1.52 ± 0.04** |
| ladder | 50 | 38.75 | 100 | 0.50 | 12 330 ± 855 | 5035 | **2.45 ± 0.17** |
| ladder | 50 | 38.75 | 200 | 0.25 | 40 079 ± 1708 | 10 070 | **3.98 ± 0.17** |
| R-collapse | **100** | **77.5** | 50 | 2 | 7253 ± 644 | 5035 | **1.44 ± 0.13** |
| R-collapse | **100** | **77.5** | 100 | 1 | 15 773 ± 915 | 10 070 | **1.57 ± 0.09** |

| R | ladder f(R) | R-collapse f(R) | difference |
|---|---|---|---|
| 1 | 1.52 ± 0.04 | 1.57 ± 0.09 | **0.4σ** |
| 2 | 1.18 ± 0.06 (interpolated) | 1.44 ± 0.13 | **1.9σ** |

**Two boxes differing by a factor 2 in N_s, L_c and τ_GP land on the same f(R) curve** — exactly at
R = 1, marginally at R = 2.

## 4. The two hypotheses, restated dimensionally correctly

| hypothesis | M = 50 (R = 2) | σ | M = 100 (R = 1) | σ |
|---|---|---|---|---|
| **A — the variable is R** (f = f(R)) | g = 118.4 | **2.1σ** | g = 153.4 | **0.5σ** |
| **B — the variable is M** (f = h(M)) | g = 153.1 | 0.6σ | g = 246.7 | **9.7σ** |
| A as I originally wrote it | 59.2 | 6.7σ | 76.7 | 8.9σ |
| B as originally written | 76.7 | 5.3σ | 123.3 | 3.8σ |

**B is excluded at 9.7σ at R = 1.** A is the only hypothesis not excluded at either mass, and it
misses "within 2σ at both" by a hair (2.1σ at R = 2).

## 5. Verdict

**Formally: NOT RESOLVED**, because L/τ = 45 and 51 fail the pre-registered sampling criterion, and
because the corrected form of A is a *post-hoc* reformulation. I will not claim a pre-registered
test that I re-stated after seeing the data.

**Substantively: R, not M.** The evidence is that f(R) collapses across a factor-2 change in box
size at R = 1 (0.4σ), and that the M-hypothesis is excluded at 9.7σ there. This is consistent with
Cencini p. 4, Sect. II.B — *"the limit N, M, L → ∞ in which we keep fixed ρ₀ = N/L and the
nondimensional mass ratio R = Nm/M"* — and with the ladder reproducing Gruber–Piasecki at R = 5
(f = 0.94 ± 0.04) and R = 2.5 (f = 1.08 ± 0.06).

## 6. What closes it

Rerun the same two cells at 65 × the **measured** τ, with prediction A stated as f(R) **before** the
data:

| M | τ measured | record needed | have | factor | steps/seed | core-h |
|---|---|---|---|---|---|---|
| 50 | 7253 | 471 432 | 327 280 | 1.4× | 28 M | 6.3 |
| 100 | 15 773 | 1 025 233 | 801 467 | 1.3× | 62 M | 13.7 |

**20 core-hours, ≈ 2.2 h at 9 concurrent.** That fixes the sampling criterion and makes A a genuine
pre-registered prediction rather than a repair. A third point at R = 4 (N_s = 100, M_d = 25) would
test the collapse where the ladder is closest to Gruber–Piasecki and is the cheapest cell of all.

**Not launched.** The rule said report the failing diagnostic and stop.
