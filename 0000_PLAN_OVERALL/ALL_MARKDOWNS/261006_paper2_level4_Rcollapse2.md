# Level 4 — R-collapse, pass 2: PRE-REGISTRATION

2026-10-06. **This section was written and committed before a single record of the new campaign was
fitted.** Pass 1 (`261005_paper2_level4_Rcollapse.md`) returned NOT RESOLVED on two counts: the
sampling criterion L/tau >= 60 failed (45 and 51), and the hypothesis that the data *did* select was
one I re-stated after seeing them. Both are repaired here by fixing the prediction first and the
records second. Nothing below §5 exists yet.

---

## 1. The invariant, stated correctly this time

Pass 1 wrote prediction A as "g = tau_T/M carries over unchanged". That is wrong, and it was my
error, not the code's: the Gruber-Piasecki time is

> tau_GP = (4/sqrt(2 pi)) · M · L_c / sqrt(m k_B T) / (1 + eta Z'/Z)    **proportional to M · L_c**

so doubling N_s at fixed eta doubles L_c and doubles tau_GP per unit mass, from 50.349 M to
**100.698 M**. A hypothesis phrased in g therefore *builds the L_c scaling out of the very thing it
is testing*. The dimensionless statement is the only admissible one:

> **"the variable is R" means f = tau_T / tau_GP is a function of R = N_s m / M alone**,
> with tau_GP evaluated with each box's own L_c.

## 2. Geometry (identical eta, identical rho_0, grid-exact)

| | ladder (pass 0) | this box |
|---|---|---|
| N_s per side | 50 | **100** |
| L_c (free compartment) | 38.75 | **77.5 = 2 x 38.75, exact** |
| wall position / box | 39.25 / 78.5 | **78.0 / 156.0** |
| divider thickness t (recorded) | 1.0 | 1.0 |
| eta = N pi r^2/(2 L_c h) | 0.10134170 | **0.10134170** |
| rho_0 = N_s/L_c | 1.29032 | **1.29032** |
| Z, eta Z', c_s (Kolafa-Rottner) | 1.239880, 0.282875, 1.749302 | same |
| tau_GP | 50.349 M | **100.698 M** |

Grid check: 77.5 x 24 = 1860, 156.0 x 24 = 3744, 78.0 x 24 = 1872 — all integers. (Computing L_c
from a rounded eta instead gives 77.4998, x24 = 1859.9945, and **fails**. Doubling is the correct
construction.) rho_0 held fixed is the point: it is the path Cencini's limit is defined along
[Cencini 2007, p. 4, Sect. II.B].

## 3. The three cells

Roman's divider eigenmode `cot K = alpha K`, alpha = M/(2 N_s m), L_eff = L_c - 2r = 76.5:

| M | R = N_s m/M | alpha | K | mode period | tau_GP |
|---|---|---|---|---|---|
| **25** | **4.00** | 0.125 | 1.3978 | 196.57 | 2517.5 |
| **50** | **2.00** | 0.250 | 1.2646 | 217.28 | 5034.9 |
| **100** | **1.00** | 0.500 | 1.0769 | 255.16 | 10069.8 |

80 seeds each, free divider, no piston, no spring, `--edmd-acc=0`, `--kbt1`, `--fixed-dt=0.4`
(dt_sigma = 1/60 exactly), `--wall-hold-steps=12000`. R = 4 is new and is the cheapest cell; it
tests the collapse where the ladder is *closest* to Gruber-Piasecki, i.e. where A and the theory
agree and B does not.

## 4. THE PREDICTIONS — fixed before the data

**A — the variable is R.** f(R) read off the N_s = 50 ladder by piecewise log-log interpolation in
R (interpolating **f**, not g; the interpolation is quoted so it cannot be re-chosen later):

| R | bracket used | weight | **A: f** | my propagation | pre-registered sigma |
|---|---|---|---|---|---|
| 4 | (2.50, 1.081) - (5.00, 0.943) | 0.6781 | **0.98** | ±0.033 | **±0.05** |
| 2 | (1.00, 1.523) - (2.50, 1.081) | 0.7565 | **1.18** | ±0.047 | **±0.06** |
| 1 | measured, no interpolation | — | **1.52** | ±0.038 | **±0.04** |

**B — the variable is M.** If tau_T(M) is the same function of M it was in the small box, then at
the same M the new box's f is halved, because tau_GP doubled: f_B = f_ladder(M)/2.

| M | source | **B: f** | sigma |
|---|---|---|---|
| 25 | log-log between (M=20, f=1.081) and (M=50, f=1.523), w = 0.2435 | **0.59** | ±0.02 |
| 50 | ladder f(M=50) = 1.523, halved | **0.76** | ±0.02 |
| 100 | ladder f(M=100) = 2.449, halved | **1.22** | ±0.09 |

Separation A/B: **1.68x at R=4, 1.54x at R=2, 1.24x at R=1** — the R=4 cell is the sharpest
discriminator of the three, which is why it was added.

**VERDICT RULE.**
1. **Sampling gates, both required at all three cells:** record/tau >= 60, and calibration response
   slope >= 0.6.
2. **A is selected** iff the measurement is within 2 sigma of A at all three cells **and** B is
   excluded at > 3 sigma at at least two of them. **Symmetrically for B.**
3. Otherwise **NOT RESOLVED**, and the failing diagnostic is reported.
4. sigma is hypot(measurement error, prediction error). Three estimators as always (modelled ACF
   fit, block average, S(0) with free amplitude); the headline is the modelled fit if its
   calibration slope clears 0.6, else the block estimator. **The bias-calibration grid is centred
   on the measured tau, not on any prediction** — the clamping trap of 2026-10-02.

## 5. Record lengths, and one deliberate deviation from the plan as written

| M | sizing basis | tau used | record [sigma] | steps | --trace-every | samples/period | seeds | L/tau if A holds |
|---|---|---|---|---|---|---|---|---|
| 25 | **1.30 x tau_GP** | 3273 | 212 725 | 12 800 000 | 580 | 20.3 | 9400-9479 (new) | 85.7 |
| 50 | measured in pass 1 | 7253 | 471 445 | 28 300 000 | 590 | 22.1 | 9200-9279 (pass 1) | 79.7 |
| 100 | measured in pass 1 | 15 773 | 1 025 245 | 61 500 000 | 690 | 22.2 | 9200-9279 (pass 1) | 66.8 |

**Seed reuse at M = 50 and 100 is deliberate**, with pass 1's `--trace-every` kept unchanged, so each
pass-2 trace is a **strict superset** of the pass-1 trace for the same seed. The pass-1 vs pass-2
difference is then a pure record-length effect with **zero seed noise** — and record length is
exactly the systematic under scrutiny (L/tau 45 and 51). Prediction A is read off the N_s = 50
ladder, which shares no trajectory with these seeds, so nothing is tuned to them. There is no
checkpointing, so reuse costs nothing and buys nothing except the clean diagnostic. M = 25 has no
pass-1 counterpart and uses fresh seeds.

**Deviation, stated because it is mine:** the plan said size R = 4 at 65 x 1.05 x tau_GP. I used
**1.30**, not 1.05. Reason: in pass 1 the new box's f came in **22 % above** the ladder's
interpolated f at R = 2 (1.44 vs 1.18). A 5 % margin at R = 4 would therefore repeat exactly the
failure this pass exists to remove, and R = 4 is the cheapest cell in the campaign — the margin
costs 1 core-hour. A longer record cannot bias tau; it can only satisfy or fail the gate. **No
prediction was touched.** With the 22 % overshoot repeated, R = 4 would still land at L/tau = 70.

Cost: **36 core-hours, 8.21 x 10^9 steps, ~4.0 h wall at 9 concurrent**; ~670 MB of reduced traces.
(The plan's "~2.5 h" was for two cells at the smaller margin.)

**Note on disk: the volume has 29 GB free of 1.8 TB (99 % full).** The campaign needs ~1 GB
including peak raw files, so it fits, but that is thin and nothing was deleted to make room.

## 6. INFERENCE — why R, and not M, is the variable the theory itself points at

*Tagged INFERENCE. No fit of f to this parameter is attempted here; it is recorded now so that it
is on the record before the data, not afterwards.*

Gruber and Piasecki state their own condition for the Boltzmann description to be exact
[SOURCE: Gruber & Piasecki, *Stationary Motion of the Adiabatic Piston*, cond-mat/9812359 /
Physica A **268** (1999) 412, **p. 2, second paragraph of the Introduction**, verbatim]:

> "In the following we study the more physical situation where M >> m. In this case the Boltzmann
> equation describes asymptotically exactly the dynamic of the piston. This is because the massive
> piston has a vanishing probability to interact back with the perturbation it causes in the states
> of the surrounding fluids (no recollisions in the sense of the kinetic theory). These perturbations
> go away to infinity, and the piston at consecutive collisions always "sees" on both sides the
> unperturbed initial equilibrium states."

In a finite box the perturbation does not go away to infinity. It is sound emitted by the divider,
and it returns after 2 L_c/c_s ~ 44 sigma-time in the small box, 89 in this one. The dimensionless
quantity that decides whether the divider still remembers its velocity when its own perturbation
comes back is

> **P = tau_v c_s /(2 L_c)**, with tau_v = M/gamma the velocity-relaxation time of the divider.

Because L_c = N_s/(n h) and the friction coefficient gamma is proportional to n h, with n and h the
same in both boxes,

> **P proportional to M/N_s = 1/R** — a function of R alone, with no residual L_c [DERIVATION].

Numerically, with gamma = 4.118 (the kinetic value used in `260930_paper2_level4_fd2.md`) and
c_s = 1.7493, **P = 0.274/R**, and it is identical in the two boxes at equal R by construction:

| R | 5 | 4 | 2.5 | 2 | 1 | 0.5 | 0.25 |
|---|---|---|---|---|---|---|---|
| **P = 0.274/R** | 0.055 | 0.069 | 0.110 | 0.137 | **0.274** | 0.548 | **1.096** |
| ladder f(R) | 0.94 | (0.98) | 1.08 | (1.18) | 1.52 | 2.45 | 3.98 |

So the theory's **own stated assumption fails along exactly the variable the data collapse on**, and
it fails where the data depart: P << 1 at R = 5 and 2.5, where f = 0.94 and 1.08 reproduce
Gruber-Piasecki; P ~ 1 at R = 0.25, where f = 3.98.

This is the same physics as the standing-wave picture seen from the other side. A divider that still
remembers its velocity after a sound round-trip is a divider locked into the box's own eigenmode —
and a mode is reversible: it stores energy, it does not transport it. That is Level 3's A4 lesson
(coherent momentum returned to the mode rather than lost to the gas) in the geometry of Level 4, and
it is the same alpha = M/(2 N_s m) = 1/(2R) that Paper 1's divider eigenmode already runs on.
**alpha and R are the same variable.** Paper 1's mass ladder and Paper 2's Level 4 are one
eigenvalue problem, and the recollision parameter is what turns that coincidence into a mechanism.

**If A is selected, this is the sentence Paper 2 makes:** not "the linear law is wrong in this box",
but *"the infinite-reservoir theory holds where the gas is the reservoir, and departs, as its own
no-recollision assumption predicts, once the divider outweighs a compartment of gas."*

---

## 7. Results

Campaign `level4_Rcollapse2_20261006`, 00:06:30 - 12:06:11, **240/240 runs, 0 aborts, 0 health
events**, 652 MB reduced. (The machine slept 01:06-09:27 with the lid closed; compute throughput
was ~2650 M steps/h throughout, matching pass 1.)

**Parameter audit, read from the recorded values and not inferred.** All three cells:
`wall_thickness_sigma = 1`, `--l0=78.0`, 200 particles as `100;100`, `radius 0.5`,
`dt_sigma = 0.0166666669`, and per cell `M_d`/`steps_after_release`/`--trace-every`/seed range
exactly 25/12 800 000/580/9400-9479, 50/28 300 000/590/9200-9279, 100/61 500 000/690/9200-9279.
**Every value matches section 5. No mixed values within any cell.** (The summary's
`eta_nominal = 0.100692` measures from the wall *position* and ignores the divider; the physical
eta from the recorded thickness is 0.10134170, as section 2.)

### 7.1 Measurements

| M | R | L/tau | cal slope | modelled | block | S(0) | spread | tau_T | g = tau/M | **f = tau/tau_GP** |
|---|---|---|---|---|---|---|---|---|---|---|
| 25 | 4 | **75** | 0.760 | 2821 ± 85 | 2948 | 2864 | 1.04x | **2821 ± 85** | 112.8 ± 3.4 | **1.120 ± 0.034** |
| 50 | 2 | **67** | 0.729 | 7049 ± 315 | 7271 | 8383 | 1.19x | **7049 ± 315** | 141.0 ± 6.3 | **1.400 ± 0.063** |
| 100 | 1 | **67** | 0.702 | 15 263 ± 589 | 17 490 | 18 636 | 1.22x | **15 263 ± 589** | 152.6 ± 5.9 | **1.516 ± 0.058** |

**Both sampling gates pass this time**, which was the whole point of pass 2: L/tau = 75, 67, 67
against the required 60, and calibration slopes 0.760, 0.729, 0.702 against 0.6, with no clamping.

### 7.2 The pre-registered test, applied exactly as written

| R | M | f measured | A: f(R) | sigma (registered) | sigma (my propagation) | B: f_lad(M)/2 | sigma |
|---|---|---|---|---|---|---|---|
| 4 | 25 | **1.120 ± 0.034** | 0.98 ± 0.05 | **2.3** | 2.8 | 0.59 ± 0.02 | **13.4** |
| 2 | 50 | **1.400 ± 0.063** | 1.18 ± 0.06 | **2.5** | 2.8 | 0.76 ± 0.02 | **9.7** |
| 1 | 100 | **1.516 ± 0.058** | 1.52 ± 0.04 | **0.1** | 0.1 | 1.22 ± 0.09 | **2.8** |

Section 4 promised both readings of the prediction error if they ever disagreed. **They do not
disagree on the verdict**: my propagated errors are *tighter* than the registered ones, so they make
A worse (2.8, 2.8, 0.1), not better.

> ## VERDICT: **NOT RESOLVED**
>
> A requires "within 2 sigma at all three". It is 2.3 and 2.5 sigma out at R = 4 and R = 2.
> **A is not selected.** B requires the same and is 13.4 and 9.7 sigma out. **B is not selected
> either, and is excluded outright.**

**The rule is not ambiguous and has not been rewritten.** What it does not capture is how unequal
the two failures are: chi² over the three points is **11.8 for A (chi²_red 3.9) against 282 for B
(chi²_red 94)** — **A is better by a factor 24** and B is dead. The honest sentence is
*"B is excluded; A's specific ladder-derived values are excluded at 2.3-2.5 sigma at two of three
points, and agree to 0.1 sigma at the third."*

### 7.3 What pass 2 did settle

**(a) The record-length worry is over, and it was never the explanation.** The nested diagnostic
reuses pass 1's seeds and cadence, so truncating pass 2 back to pass 1's length isolates record
length with **zero seed noise**:

| M | pass 1 record | pass 1 tau | pass 2 truncated to it | pass 2 full | full/short |
|---|---|---|---|---|---|
| 50 | 329 426 | 7253 | 7403 ± 273 | 7049 ± 315 | **0.952x** |
| 100 | 803 838 | 15 773 | 15 837 ± 868 | 15 263 ± 589 | **0.964x** |

**Lengthening the record by 1.4x and 1.3x moved tau by only -4.8 % and -3.6 %.** Pass 1's L/tau = 45
and 51 were a real protocol violation but a small numerical one, and pass 1's central values stand:
pass 1 -> pass 2 moved f by -2.8 % (0.3 sigma) at R = 2 and -3.4 % (0.5 sigma) at R = 1.

**(b) The residual is not an interpolation artefact.** A at R = 4 and R = 2 is *interpolated* from
the ladder while A at R = 1 is a *measured* ladder node — and the agreement is perfect exactly where
A is measured. That invites the explanation "the interpolation is wrong", and it is checkable:
refitting the five ladder points with a smooth quadratic in ln f vs ln R (convex, curvature
+0.131) instead of piecewise chords gives

| R | measured | piecewise (registered A) | smooth quadratic | measured/smooth |
|---|---|---|---|---|
| 4 | 1.120 | 0.980 | 0.970 | **1.155** |
| 2 | 1.400 | 1.180 | 1.145 | **1.223** |
| 1 | 1.516 | 1.520 | 1.533 | **0.989** |

**The smooth curve is slightly *lower*, so it makes A marginally worse, not better.** The residual
is real: **the N_s = 100 box sits 15-22 % above the N_s = 50 ladder's f(R) at R = 4 and R = 2, and
on top of it at R = 1.**

### 7.4 What closes it, and it is cheap

The residual sits exactly where f(R) must be *interpolated* from the ladder, and nowhere else. That
is either a genuine box-size dependence at fixed R, or structure in the ladder's own f(R) between
its nodes (R = 5, 2.5, 1, 0.5, 0.25) that no interpolation can know about. **The two are separated
by measuring the ladder at the missing R directly** — N_s = 50 with M_d = 12.5 (R = 4) and
M_d = 25 (R = 2) — after which both boxes have measured points at R = 4, 2 and 1 and **no
interpolation enters the comparison at all**.

Cost at 65 tau: M = 12.5 needs ~2.5 M steps/seed, M = 25 ~5.7 M; 80 seeds each is **652 M steps,
~2.9 core-hours**. **Not launched** — no new campaign without a go.
