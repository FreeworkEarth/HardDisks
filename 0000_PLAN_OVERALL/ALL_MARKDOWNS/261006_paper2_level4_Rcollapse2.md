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

*Empty. To be filled only after the campaign completes and the three estimators run.*
