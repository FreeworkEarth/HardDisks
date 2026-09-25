# Level 4, next batch — the mode across the ladder, and the demonstration run

2026-10-07. Campaign A (R-collapse pass 2) is running; its pre-registration is committed as 3af501c
and its results go in `261006_paper2_level4_Rcollapse2.md`. **Section 3 below is the
pre-registration for the demonstration run and was written before a single record of it existed.**
Section 1 is analysis of traces already on disk — no new runs, no core changes.

---

## 1. The divider eigenmode across alpha = 0.1 – 2.0, in two boxes (item C)

The ladder's `tau_T` fits fixed omega at the Kolafa–Rottner prediction, so the mode period had been
measured at one mass only (M = 10, `260930`). Here **omega is free at every mass in both boxes**,
which turns "the oscillation is Paper 1's divider mode" from one agreement into a seven-point test
of the eigenvalue equation itself. Figure: `261006_p2_mode_ladder.{png,pdf}`.

**The error bar, and why it is not the jackknife.** The first pass of this analysis printed
deviations of "+179 sigma" using the block-jackknife error on the fitted centre — 0.01 to 0.12 %.
That is the mirror image of the mistake `260930` recorded: 260930 quoted a 0.03 % *agreement* off a
5 % measurement, and this would quote a 100-sigma *disagreement* off a 0.02 % statistical error with
a systematic ten times larger sitting unmeasured beside it. The error used is therefore
**max(block jackknife, half the spread between the two independent observables)** — divider position
and temperature difference. That spread is 0.88 % at M = 10 and below 0.1 % for M >= 50, so it
dominates exactly where the jackknife is least believable. The spectral FWHM (0.8–6 %) is reported
as the **line width** and never as the error on the centre.

**eta convention, and a 0.28 % shift.** `261003` predicted periods at eta = 0.100051, the summary's
`eta_nominal`, which measures from the wall *position* and so ignores the 1.0-sigma divider. The
audited value uses the free compartment length and is **0.10134170 in both boxes**, giving
c_s = 1.749302 instead of 1.744337. Every predicted period here is 0.28 % shorter than 261003's.

| box | N_s | L_c | M | alpha | R | period measured | sigma | period (dT) | spectral peak | FWHM | KR | Paper 1 c_s | ideal |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A | 50 | 38.75 | 10 | 0.100 | 5.00 | **94.47** | ±0.41 | 93.64 | 95.84 | ±6.0% | 94.89 | 93.95 | 117.38 |
| A | 50 | 38.75 | 20 | 0.200 | 2.50 | **102.54** | ±0.18 | 102.19 | 102.25 | ±3.5% | 103.20 | 102.17 | 127.66 |
| A | 50 | 38.75 | 50 | 0.500 | 1.00 | **124.63** | ±0.05 | 124.54 | 124.69 | ±2.2% | 125.91 | 124.65 | 155.75 |
| A | 50 | 38.75 | 100 | 1.000 | 0.50 | **155.74** | ±0.03 | 155.71 | 155.15 | ±1.3% | 157.60 | 156.03 | 194.95 |
| A | 50 | 38.75 | 200 | 2.000 | 0.25 | **204.81** | ±0.02 | 204.80 | 204.67 | ±0.8% | 207.56 | 205.48 | 256.74 |
| B | 100 | 77.50 | 50 | 0.250 | 2.00 | **215.50** | ±0.10 | 215.33 | 215.14 | ±1.7% | 217.28 | 215.11 | 268.77 |
| B | 100 | 77.50 | 100 | 0.500 | 1.00 | **252.49** | ±0.04 | 252.44 | 253.27 | ±1.4% | 255.16 | 252.61 | 315.62 |

### 1.1 Inverted: K measured against cot K = alpha K

| box | M | alpha | K from cot K = aK | K measured (KR c_s) | departure | vs Paper 1 c_s | vs ideal gas |
|---|---|---|---|---|---|---|---|
| A | 10 | 0.100 | **1.4289** | 1.4354 ± 0.0063 | **+0.45 ± 0.44 %** | −0.55 % | +24.3 % |
| A | 20 | 0.200 | **1.3138** | 1.3223 ± 0.0023 | **+0.65 ± 0.17 %** | −0.36 % | +24.5 % |
| A | 50 | 0.500 | **1.0769** | 1.0879 ± 0.0004 | **+1.03 ± 0.04 %** | +0.02 % | +25.0 % |
| A | 100 | 1.000 | **0.8603** | 0.8706 ± 0.0002 | **+1.20 ± 0.02 %** | +0.18 % | +25.2 % |
| A | 200 | 2.000 | **0.6533** | 0.6620 ± 0.0000 | **+1.34 ± 0.01 %** | +0.33 % | +25.4 % |
| B | 50 | 0.250 | **1.2646** | 1.2751 ± 0.0006 | **+0.83 ± 0.05 %** | −0.18 % | +24.7 % |
| B | 100 | 0.500 | **1.0769** | 1.0883 ± 0.0002 | **+1.06 ± 0.02 %** | +0.05 % | +25.0 % |

**The one shared alpha = 0.500 is the test that matters.** Box A (N_s = 50, M = 50) gives
**+1.03 ± 0.04 %**; box B (N_s = 100, M = 100) gives **+1.06 ± 0.02 %**. They agree to **0.03 %**
while L_c, N and the period itself all differ by a factor 2 (124.63 vs 252.49). **The departure
from the prediction is a function of alpha and not of box size** — the same statement the R-collapse
is testing for tau_T, holding independently for the mode.

Three further things this table settles:

1. **The departure grows monotonically with alpha**, +0.45 % to +1.34 %, a factor 3. **No single
   sound speed can absorb it**, so it is not a c_s error; it is an alpha-dependent correction to
   `cot K = alpha K` of at most 1.3 %. **OPEN**, and small.
2. **Paper 1's independently measured c_s (+1.01 % on bulk KR) leaves a residual of −0.55 % to
   +0.33 %, crossing zero near alpha = 0.5.** This is the **fourth** independent appearance of that
   +1 % — after Paper 1's c_s, Level 3's box stiffness and the M = 10 mode — and the first time it
   has been seen across a range of alpha and two box sizes at once.
3. **The ideal gas is +24.9 % away on average.** Excluded absolutely, at every alpha, in both boxes.

### 1.2 tau_r: Mansour's piston form becomes *exact* in the heavy limit

| box | N_s | M | M_hat = M + mN/3 | Δf/f Mansour | tau_r Mansour | tau_r measured | sigma | Mansour/measured | tau_r/L_c |
|---|---|---|---|---|---|---|---|---|---|
| A | 50 | 10 | 43.33 | 0.0711 | 425 | **208** | ±3 | **2.04×** | 5.37 |
| A | 50 | 20 | 53.33 | 0.0641 | 513 | **298** | ±6 | **1.72×** | 7.68 |
| A | 50 | 50 | 83.33 | 0.0513 | 782 | **586** | ±9 | **1.33×** | 15.11 |
| A | 50 | 100 | 133.33 | 0.0405 | 1237 | **1085** | ±9 | **1.14×** | 28.01 |
| A | 50 | 200 | 233.33 | 0.0306 | 2156 | **2090** | ±12 | **1.03×** | 53.93 |
| B | 100 | 50 | 116.67 | 0.0306 | 2257 | **1209** | ±23 | **1.87×** | 15.61 |
| B | 100 | 100 | 166.67 | 0.0256 | 3168 | **2091** | ±35 | **1.52×** | 26.98 |

**This replaces the 260930 statement.** That report, working at M = 10 alone, said the mode "damps
2.2× faster than Mansour's piston form" and left it OPEN as a disagreement. With the ladder in hand
the ratio is seen to fall **monotonically from 2.04 at M = 10 to 1.03 at M = 200**: Mansour's form
is a *piston* form, and it assumes the divider carries the inertia. At alpha = 0.1 the gas standing
wave carries it instead and the form is 2× too slow; by alpha = 2.0 the divider does dominate and
the form is **exact to 3 %**. **It is not a standing disagreement — it is a crossover, and the
location of the crossover is the result.**

**tau_r does not collapse on alpha, and it does not collapse on R.** It scales with L_c at fixed M:
tau_r/L_c = 15.11 (box A) vs 15.61 (box B) at M = 50 — 3.2 % apart — and 28.01 vs 26.98 at M = 100
— 3.8 % apart. Mansour's form predicts box ratios 2.89 and 2.56 against the measured 2.06 and 1.93,
so **its box scaling is ~40 % too strong**. Three quantities, three different invariants:
tau_T on R, the mode period on alpha, tau_r on L_c at fixed M. **OPEN.**

---

## 2. Batch status

| item | state |
|---|---|
| A — R-collapse pass 2 | **running**, launched 00:06:30, pre-registration committed 3af501c |
| B — demonstration run | **pre-registered below, chained to launch when A exits** |
| C — mode ladder | **done**, section 1, figure `261006_p2_mode_ladder` |
| D — report, status lines, commit script | this file; `commit_20261007.sh` |
| plan Level 4 section, Paper 2 draft | **untouched**, as instructed |

**Two operational notes.** The machine slept 01:06–09:27 (8 h 21 m, lid closed), so A is ~8 h behind
the clock but exactly on schedule in compute: batches of 9 take ~13 min at M = 100, matching pass 1.
`caffeinate -i -s` is now bound to the driver's lifetime. And the volume has **29 GB free of 1.8 TB
(99 % full)**; this batch needs ~1.7 GB, so it fits, and nothing was deleted.

---

## 3. PRE-REGISTRATION — the demonstration run (item B)

*Written before any record of `level4_demo_20261007` existed.*

### 3.1 Geometry and the causal argument

The ladder's own box, so tau_T is already measured here: N_s = 50 per side, L_c = 38.75, wall 39.25,
box 78.5, t = 1.0, eta = 0.10134170, **M_d = 200**, divider free. Piston on the right, step
protocol, **travel 3.875 sigma = 10 % of L_c** (grid: 3.875 × 24 = 93 exact; right wall
78.5 → 74.625, × 24 = 1791 exact). u = 0.2, so the push lasts **19.38 sigma-time**.

**The plan's reason for "the divider cannot follow" does not survive checking, but the conclusion
does, for a better reason.** The plan said the 20-sigma push is short against tau_v = 486. In fact
tau_v = M/gamma = 200/4.118 = **48.6** — gamma = 4.118 is the two-sided kinetic friction already used
for the recollision parameter, and 486 is ten times it. The correct argument is **causal**: the
compression cannot reach the divider before

> L_c/c_s = 38.75/1.7493 = **22.15 sigma-time**, which is *after* the push ends at 19.38.

The divider is not merely slow to respond — it has not yet been told. That is also why the v3 runs
failed: at u = 0.05 the push lasted 77.5 sigma, three and a half sound transits, and the divider had
every opportunity to follow.

### 3.2 Predictions, fixed before the data

**(i) The work lands entirely on the right.** From the Kolafa–Rottner adiabat
(d ln T = −Z d ln L, integrated over eta = 0.101342 → 0.112602):

| | T_right after push | W_in | initial T_1 − T_2 |
|---|---|---|---|
| **Kolafa–Rottner** | **1.14143** | **7.0715 kT** | **0.14143** |
| ideal gas (Z = 1) | 1.11111 | 5.5556 kT | 0.11111 |

**ΔKE_left/W_in < 1 % at t = 19.38**, by the causal argument. **W_in itself is a 27.3 % KR-vs-ideal
discriminator measured on the piston's own ledger** — a mechanical quantity, not a fluctuating one,
so it is far better determined than the 80-seed temperature noise. This is the sharpest test in
the run.

**(ii) The decay is two-stage, and — stated in advance — the thermal stage is NOT observable here.**
Stage 1 is mechanical, on ~10 tau_r = 20 900 sigma. Stage 2 is thermal, tau_T(200) = 40 079 ± 1708.
They are cleanly separated by a factor 19. But **the amplitude left for stage 2 is zero in the
reversible limit**:

> After the push both gases are still at their **original entropy** (piston compression and divider
> motion are both adiabatic). Equal N and equal entropy at equal pressure fix the same state, so
> mechanical equilibrium is reached at **equal temperature**. The computed equal-pressure position
> is x = 36.812500 against the equal-T-equal-P value 36.812500, and **T_1 − T_2 = 0 to machine
> precision** [DERIVATION, verified numerically against the KR adiabat].

What survives is only the irreversible part, which is localised in the right gas (the compression
wave is launched there and damps there). At Mach = u/c_s = 0.114 the weak-piston excess is
O(Ma²) = 1.3 % of W_in = 0.092 kT, i.e. a residual **T_1 − T_2 ~ 0.0019 against an 80-seed noise
floor of 0.0222**. Reaching 3 sigma on it would need **~104 000 seeds**.

**So prediction (ii) as the plan writes it cannot be tested by this protocol, and that is itself the
result**: compressing one side of a two-compartment box with a free adiabatic divider, reversibly,
produces *no lasting temperature difference at all*. It is the cleanest possible statement of why
the fluctuation route exists, and it belongs in the paper next to the figure. The pre-registered
prediction is therefore: **|T_1 − T_2| falls from 0.1414 to below the 0.0222 noise floor on the
MECHANICAL timescale (~2 × 10⁴ sigma), not on tau_T (4 × 10⁴), and no exponential with
tau = 40 079 is recoverable from the remainder.**

**(ii') ADDED 2026-10-07, before `level4_demo_20261007` had written a single record.** The
reversible result of (ii) is stronger than "the ideal-gas case is degenerate": both gases start in
the same state, each is compressed isentropically, and at mechanical equilibrium they share P and
N, so the same (S, P, N) fixes the same T. **T_1 - T_2 = 0 at mechanical equilibrium for ANY
equation of state**, not only the ideal one [DERIVATION]. That inverts the experiment's purpose:

> **any T_1 - T_2 surviving mechanical equilibrium IS dissipated work**, and the B1long/B2pilot pair
> therefore measures irreversibility directly rather than heat conduction.

The reported quantity is the **dissipated fraction**

> **D = (T_1 - T_2) N k / W_in**, N = 50 per side, with the error taken over seeds,

and the pre-registered expectation is **D ~ Ma² = (u/c_s)²** [INFERENCE, from the weak-piston
excess-work scaling]:

| cell | u | Ma = u/c_s | **predicted D** | predicted T_1 - T_2 | 8- or 80-seed noise on T_1 - T_2 |
|---|---|---|---|---|---|
| **B1long** | 0.2 | 0.1143 | **0.0131** | 0.00185 | 0.0222 (80 seeds) |
| **B2pilot** | 1.0 | 0.5716 | **0.327** | 0.0462 | 0.0704 (8 seeds) |

**CORRECTED 2026-10-08, and the correction strengthens the test.** The retention bound was
first written as "the whole excess stays in gas 1", giving D = excess/W_in. That is a factor two
loose: heat deposited in gas 1 alone still drives the divider, gas 1 expands against gas 2, and for
E = N k T the mechanical equalisation hands a **quarter** of that heat across. **Full retention is
therefore D = (1/2) x excess/W_in**, i.e.

| cell | u | excess = W_in - W_qs | **full-retention T_1 - T_2** | measured | **exclusion** |
|---|---|---|---|---|---|
| B1long | 0.2 | 1.369 | **0.0137** | -0.0028 ± 0.0111 | **1.5 sigma — not excluded** |
| B2pilot | 1.0 | 13.594 | **0.1359** | +0.0088 ± 0.0410 | **3.1 sigma** |

So **"the free divider shares dissipated work" rests on B2pilot alone: eight seeds, 3.1 sigma.**
B1long cannot exclude full retention at all. The 32-seed `M200_u100` cell added to Level 4b on
2026-10-08 carries this claim; until it reports, the statement in section 4.4 is quoted with its
provenance and does not go into a draft.

**The Ma² form is a weak-piston expansion and is not trustworthy at Ma = 0.57**; it is written down
anyway, before the data, so that the size of its failure is a measurement rather than a
rationalisation. What is pre-registered without qualification is the **ordering**: D(u = 1.0) >
D(u = 0.2), and D > 0 at both. B2pilot's 8 seeds give a noise floor of 0.0704 on T_1 - T_2, so
**even the predicted 0.0462 is below 1 sigma there** — B2pilot decides whether an 80-seed
irreversible cell is worth running, and is not itself expected to resolve D.

**(iii) The divider moves LEFT, by exactly half the piston travel, and does not come back.** The
plan says it "drifts toward the hot side by the isobaric amount and returns to centre". Both halves
are wrong. The right gas is *compressed as well as heated*, so its pressure is higher
(0.2078 vs 0.1600 immediately after the push) and the divider is pushed **away** from the hot side.
It settles at

> x = 36.8125, i.e. **1.9375 sigma to the LEFT = exactly half the 3.875-sigma piston travel**,

and stays there, because the box is permanently shorter. (Gruber–Piasecki's drift *toward* the hot
region is a second-order effect at *equal* pressures and is negligible against this.) **The
testable content is the timescale**: the approach to 36.8125 is governed by tau_r = 2090, not by
tau_T = 40 079 — a factor 19, and the direct demonstration that the two stages are separate.

**(iv) The ledger closes.** Total energy conserved to the summary's own accounting;
health contract strictly zero.

### 3.3 Fit rule

Seed-averaged T_1 − T_2, A_inf = 0, over the window where the signal exceeds 3 × the seed noise
(3 × 0.199/sqrt(80) = 0.0667). Report tau with its error against **both** tau_T = 40 079 ± 1708 and
tau_GP = 10 070. If, as predicted in (ii), the window closes before stage 2 begins, report that and
the bound, and do not fit.

### 3.4 Cells, and one addition

| cell | u | seeds | steps | record | --trace-every | dt | purpose |
|---|---|---|---|---|---|---|---|
| **B1long** | 0.2 | 80 | 7 250 000 | 120 833 = 3.0 tau_T | 600 | 10 | the decay |
| **B1zoom** | 0.2 | 80 | 6 000 | 100 sigma | 1 | 1/60 | **the push itself** |
| **B2pilot** | **1.0** | **8** | 7 250 000 | 120 833 | 600 | 10 | **the irreversible residual** |

**B1zoom** exists because at dt = 10 the 19.38-sigma push is two samples, which is useless for the
first panel. Same seeds and flags, only `--steps` differs, so it is a **strict prefix** of B1long.

**B2pilot is an addition to the plan, and the reason is section 3.2(ii).** At u = 0.2 the push is
near-reversible and stage 2 has no amplitude. u = 1.0 is **Mach 0.57 — strongly irreversible but
not a shock** — and should leave a residual roughly (0.57/0.114)² = 25× larger, i.e. of order the
noise floor rather than a tenth of it. **8 seeds, not 80**: this measures the residual and decides
whether an 80-seed irreversible cell is worth running, before anyone spends one. ~0.25 core-hours.

### 3.5 Figure, and the geometry pictures

Three panels, log time axis, error bands throughout: ΔKE_right and ΔKE_left vs t (B1zoom for
t < 100, B1long after); T_1 − T_2 vs t with **both** predicted decays drawn and the noise floor
shaded; divider position vs t with 36.8125 marked and tau_r, tau_T marked as spans.

**The two geometry pictures need a GUI session and cannot be produced headless from here.** The
commands are in `paper2_geometry_pictures.sh`; `--render=experiment` and `--render=paper` both open
an SDL window and capture with the `S` key. Flagged, not silently skipped.

---

## 4. Results — the demonstration run

`level4_demo_20261007`, 12:07:25–12:15:13. **B1long 80/80, B1zoom 80/80, B2pilot 8/8, 0 aborts,
0 health events**, 95 MB. Recorded parameters match section 3.1 in every cell:
`wall_thickness_sigma = 1`, `L0 = 39.25`, `piston_right_travel_sigma = 3.875`,
`piston_target_sigma = 74.625`, `wall_mass_factors = 200`, `steps_after_release = 7 250 000`,
record 120 830 sigma.

### 4.1 (i) The work lands entirely on the right — confirmed

At the recorded end of the push (20.650 sigma, identical across all 80 seeds), from **B1zoom**:

| | measured | as a fraction of W_in | prediction |
|---|---|---|---|
| Delta KE_right | **+8.4438 ± 0.2362** | **+100.0 %** | ~ W_in |
| Delta KE_left | **−0.1838 ± 0.1156** | **−2.18 ± 1.37 %** | < 1 % |
| divider displacement | **+0.0376 ± 0.0653** | — | ~ 0 |

**The causal argument holds.** The divider has not moved (0.6 sigma from zero) and the far gas has
not been touched — Delta KE_left is 1.6 sigma from zero, and its magnitude is within 1 sigma of the
predicted 1 % bound. All of the piston's work is in the pushed gas when the push ends.

**One correction to section 3.1, from the recorded stop time.** The push ends at **20.650**, not
19.375: the piston overshoots its travel by ~0.27 sigma (confirmed against B2pilot's 4.150 vs
3.875). The causal margin against the 22.15-sigma transit is therefore **1.5 sigma-time, not 2.8** —
still clean, and the measurement above confirms it, but with half the headroom claimed.

### 4.2 (ii) The thermal stage is unobservable, exactly as pre-registered

| t | T_1 − T_2 (B1long, 80 seeds) |
|---|---|
| **push end, 20.65** | **+0.17255** |
| 30 (first trace sample after it) | **+0.1507 ± 0.0077** |
| +0.5 tau_r (1070) | **−0.0059 ± 0.0177** |
| +1 tau_r (2120) | −0.0010 ± 0.0175 |
| +2 tau_r (4210) | +0.0004 ± 0.0197 |
| mechanical equilibrium, t > 10 tau_r, time-averaged | **−0.00282 ± 0.01109** |

Noise floor 0.0222; the pre-registered fit window requires |T_1 − T_2| > 3 x noise = 0.0667.
**The signal is below that floor by 0.5 tau_r — about 1000 sigma-time — while tau_T = 40 079.**
The window closes long before the thermal stage begins, so **by the rule in section 3.3 no fit was
performed**. The prediction as written in 3.2(ii) is confirmed, and it was confirmed *faster* than
predicted: the mechanical stage needs ~0.5 tau_r, not the ~10 tau_r allowed for.

At mechanical equilibrium T_1 − T_2 is **−0.0028 ± 0.0111, i.e. 0.25 sigma from zero.**

**A label error corrected 2026-10-08.** The first version of this table called the t = 30 sample
"push end" and then could not reconcile its 0.1507 with Delta KE_right/N = 8.4438/50 = 0.1689. There
is nothing to reconcile: **the push ends at 20.65, where T_1 - T_2 = 0.17255**, which is exactly
Delta KE_right/N + |Delta KE_left|/N = 0.16888 + 0.00368. The drop to 0.15071 by t = 30 is
**1.09 kT leaving the two gases in 9.35 sigma-time** — the pulse reaching the divider at 22.15 and
kicking it, predicted at 0.765 kT (section 1B of `261008`). Two different times were quoted for two
different quantities; that was sloppiness, not a discrepancy.

### 4.3 (iii) The divider moves left by half the piston travel — confirmed, including the sign

| cell | divider displacement at mechanical equilibrium | prediction | sigma |
|---|---|---|---|
| B1long (u = 0.2, 80 seeds) | **−1.9058 ± 0.1447** | −1.9375 | **0.2** |
| B2pilot (u = 1.0, 8 seeds) | **−2.0406 ± 0.4752** | −1.9375 | **0.2** |

**Away from the hot side, by exactly half the 3.875-sigma piston travel, and it does not return.**
The plan's "toward the hot side ... returns to centre" is excluded in both direction and
destination.

### 4.4 (ii') Dissipation, and the failure of the Ma² form

| cell | u | Ma | W_in | **D = (T_1−T_2) N/W_in** | predicted Ma² |
|---|---|---|---|---|---|
| B1long | 0.2 | 0.114 | **8.4403 ± 0.1840** | **−0.0167 ± 0.0657** | 0.0131 |
| B2pilot | 1.0 | 0.572 | **20.6653 ± 2.9402** | **+0.0213 ± 0.0991** | **0.3268** |

At u = 1.0 the Ma² prediction is **3.1 sigma above the measurement**. The weak-piston form fails at
Mach 0.57, as section 3.2 said it would — the size of the failure is now measured rather than
assumed. The pre-registered **ordering** D(1.0) > D(0.2) is nominally satisfied (+0.021 vs −0.017),
**but neither D differs from zero, so the ordering test is uninformative and is reported as such.**

**The physically important number is W_in, not D.** At u = 1.0 the piston does **20.67 kT** against
a Kolafa–Rottner isentropic **7.07** — **13.6 kT, two thirds of the work, is dissipated** — and the
temperature difference surviving mechanical equilibrium is still **+0.009 ± 0.041**. If that
dissipation were deposited in the pushed gas alone it would leave T_1 − T_2 ~ 0.14; the 2-sigma
bound is 0.09. **The free divider shares the dissipated work between the compartments as it
equilibrates — at 3.1 sigma, on eight seeds** (see the corrected bound in section 3.2(ii'); the
2-sigma bound quoted above is superseded by that calculation). That is the strongest form of the section-3.2 result: *no push protocol of this
kind — reversible or violently irreversible — creates a temperature difference to watch decay*,
which is why the fluctuation route is the only route.

### 4.5 (iv) The W_in ledger verdict — and my own pre-registration was wrong about it

> **W_in(u = 0.2) = 8.4403 ± 0.1840 kT.** Kolafa–Rottner isentropic 7.0715: **+19.4 %**.
> Ideal-gas isentropic 5.5556: **+51.9 %**.

W_in is work done *by* the piston, so it is bounded below by the isentropic work and dissipation
only adds. **The test is one-sided, and neither equation of state is excluded.** Section 3.2 called
this "a 27.3 % discriminator ... the sharpest test in the run". **That was wrong, and the error was
mine**: it assumed the dissipation would be O(Ma²) = 1.3 %, and the measured dissipation is
**19.4 %, fifteen times larger and comparable to the 27.3 % separation it was supposed to resolve.**

**The test is recoverable, and Level 4b already contains it.** Dissipation measured at two speeds
gives W_in/W_isen − 1 = 0.194 (u = 0.2) and 1.922 (u = 1.0), scaling as u^1.42. Extrapolated to
**u = 0.05** that is a **2.7 % correction**, against a 27 % KR-vs-ideal separation — a clean
ten-to-one. **The `M200_u005` cell of Level 4b is therefore the real ledger test**, and it is
already on disk.

### 4.7 Recorded values, read rather than inferred

**(a) The piston travels exactly as far as it was told, and runs slow.** `piston_target_sigma`
= 74.625 against a start of 78.5 is a travel of **3.8750** in both cells. The excess is entirely in
the *time*: mean speed 0.18765 (93.8 % of 0.2) and 0.93374 (93.4 % of 1.0). **Not an overshoot and
not a trigger latency — both of which I had inferred from stop times alone, and both wrong.**
The divider agrees independently: 2 x 1.9058 ± 0.1447 = **3.81 ± 0.29**, 0.22 sigma from 3.875 and
1.2 sigma from 4.145. **W_qs depends on travel, not time, so it stays 7.0715 and section 4.5
stands.**

**(a2) The 0.25 sigma, decided by measurement AND code — and my three earlier accounts of it were
all inferences.** I said in turn 1 that it was a fixed 1.275-sigma trigger latency, in turn 2 that it
was a 0.27-sigma travel overshoot, and in turn 3 that `piston_target_sigma` settled it. The first two
were wrong and the third was right for the wrong reason. **Three independent lines now agree, and
one of them is the source code:**

1. **Measurement — `PistonWork(t)` in B1zoom at dt = 1/60.** Across **all 80 seeds the piston does
   exactly zero work for t < 1.25**, and the earliest nonzero increment is at t = 1.26667
   (median 2.108, max 8.917 — the spread is the wait for the first particle collision). 0.25/u at
   u = 0.2 is **1.25**. There is a real gap and the gas is not in it.
2. **Measurement — `SegEtas`.** Total gas length 77.4998 at t = 0 and 73.6235 at the push end:
   a compression of **3.8763**, i.e. the requested 3.875.
3. **Code — `compute_segment_bounds()`, 00ALLINONE.c:3183.** The right-hand boundary of the gas is

   > `float right_plane = fminf((float)XW2, piston_right_x);`

   so it is **the box wall XW2 = 78.5 while the piston is outside it, and the piston only after it
   crosses**. `piston_right_x` is a gas-facing **plane, not a centre** — it enters `right_plane`
   with no half-thickness — and it is parked at 78.75, exactly 0.25 outside. That is why SegEtas
   reads 38.75 at t = 0 and piston-minus-divider later: **the convention never changes, the `fminf`
   does.** (The left boundary is `fmaxf(XW1, piston_left_x + 5.0f)` — a 5-px offset on the left
   piston that the right one does not have. Inactive here, but noted.)

4. **Code — `handle_piston_collisions()`, 00ALLINONE.c:7798.** The work counter is gated on the
   particle reaching the piston and **on nothing else** — there is no `XW2` test anywhere in it:

   > `if (X[i] + Radius[i] >= piston_right_x) {`
   > `    float vrel = Vx[i] - vx_piston_right;`
   > `    if (vrel > 0.f) { ... piston_work_right += dE; ... } }`

   So **any** piston collision is counted, wherever the piston is. That closes the argument by code
   at both ends: had the gap contained gas, work would have been recorded before t = 1.25, and in
   80 seeds it was not. (That is the CCD path; these runs are EDMD, where the counter accumulates
   the engine's own delta, `piston_work_right += (wR - edmd_prev_work_R)` at 00ALLINONE.c:16930 —
   also ungated on XW2.)

**Conclusion, now resting on code and not on inference: the piston crosses a real 0.25-sigma gap at
constant speed doing no work, then compresses the gas by exactly the requested 3.875.**

| quantity | value |
|---|---|
| initial compartment length / total | 38.7499 / **77.4998** |
| **compression** | **3.875** |
| initial eta (both sides) | **0.10134170** |
| final eta, pushed side | 0.1126019 |
| c_s, Z_acoustic at that eta | 1.749302, **2.2572** |
| **W_qs, Kolafa–Rottner** | **7.0715** |
| W_qs, ideal gas | 5.5556 |
| **dissipation, u = 0.2 / u = 1.0** | **19.4 % / 192 %** |
| **tau_push (the compression)** | **d/u = 19.375, 7.75, 3.875** — `piston_stop_t_rel` is longer by 0.25/u, which is the gap |
| 4b ledger targets: intercept / slope Zd | 7.0715 / 8.746 (KR); 5.5556 / 7.071 (ideal) |

**Section 4.5 stands as written, and the "pending geometry" hold on it is released.**

**(a3) LEVEL 3 IS UNAFFECTED — checked the same day, before any 4b analysis.** The concern was that
geometry C's closing number, a −1.3 ± 0.6 % comparison between the settled displacement and an
F(L)-based model built on the piston's travel, would move if the true compression there were also
0.25 longer than the travel parameter. **It is not.** Reproducing
`level3_v6_20260924/k0.5_M200_u0.05` with dt = 1/60 and reading the recorded values:

> box XW2 = 2 x 54.75 = **109.5**; piston parked at **109.7500** = XW2 + 0.25; recorded
> `piston_target_sigma` = **101.539998** = XW2 − 7.96; **XW2 − target = 7.9600 = the recorded travel
> parameter exactly**; `PistonWork` is **zero for all t < 5.000 = 0.25/u**, first nonzero at 5.933;
> and the gas at t = 0 runs to **109.4998**, i.e. to the box wall, not to the piston.

**The target is computed as XW2 − travel in both geometries, and the gap is crossed without doing
work, so the compression equals the recorded travel parameter in Level 3 exactly as in Level 4.
s̄_model vs 0.8705 ± 0.0049 is unchanged and Level 3 stays closed.**

**(b) T is KE/(N k), with no mean-flow subtraction.** Checked directly: (KE_R - KE_L)/N and
((Delta KE_R) - (Delta KE_L))/N agree **to better than 1e-5 at every time**, because every seed
starts at exactly KE = 50.0000 per side. The D and the "Delta T = 0" identity are therefore stated
in the same variable the trace records.

**(c) The per-seed energy ledger closes.** Collisions are elastic and there is no bath, so
W_in = Delta KE_R + Delta KE_L + Delta KE_div must hold **per seed**, not just on average. With
KE_div from the time derivative of the recorded divider position (B1zoom, dt = 1/60, so this is
accurate), at the push end over 80 seeds:

> mean residual **-0.0005**, sd **0.0063**, **max |residual| = 0.0506 kT** — 0.6 % of W_in
> in the worst seed.

**(d) Per-seed W_in, B2pilot, all eight:** 12.635, 14.348, 15.387, 16.559, 19.157, 20.164, 32.795,
34.279. Mean 20.665, sd 8.316 (40 % per seed, as expected for four or five piston collisions at
Mach 0.57), **max |z| = 1.64 — no outlier**. The distribution is visibly skewed, six seeds in
12-21 and two near 33, but eight seeds cannot diagnose that; the 32-seed `M200_u100` cell will.

### 4.8 The B1long time series — and a correction to 4.2

*Post-hoc on B: this test was specified in section 1B of `261008` after B had run.*

**The pre-registered single exponential on [30, 1070] fails, and the reason is the result.**

> A = **-0.0420 ± 0.0056** (negative), tau = **166 ± 33**, **chi²_red = 77.2** over 104 points.

A chi²_red of 77 and a negative amplitude are not a measurement of a sharing time. **The 80-seed
mean T_1 - T_2 is dominated by a coherent oscillation**, not by a decay:

| t | T_1 - T_2 | divider dx |
|---|---|---|
| 20 | **+0.1685 ± 0.0062** | +0.040 ± 0.063 |
| 50 | +0.0536 ± 0.0085 | -1.506 ± 0.102 |
| **100** | **-0.1530 ± 0.0057** | **-4.223 ± 0.093** |
| 150 | +0.0031 ± 0.0096 | -2.221 ± 0.126 |
| **200** | **+0.1645 ± 0.0092** | **-0.113 ± 0.106** |
| 300 | -0.1215 ± 0.0102 | -3.853 ± 0.137 |
| 1000 | -0.0329 ± 0.0155 | -2.887 ± 0.210 |
| 3000 | +0.0335 ± 0.0171 | -2.200 ± 0.255 |

**Is the mode visible in the 80-seed mean? Emphatically yes.** The residual spectrum after removing
the fitted exponential peaks at a period of **173 sigma-time** with **8.6 x 10^4 times the median
residual power** (dt = 10 samples the 208-period mode 21 times over, so this is not aliasing).
**The divider is kicked into a large coherent oscillation — amplitude about 2.1 sigma about its new
equilibrium at -1.94, period ~200 — and T_1 - T_2 swings with it** because the divider
alternately compresses each side. That exchange is **reversible**: it is the mode storing and
returning energy, not heat crossing the divider.

**The envelope, not the instantaneous value, is what decays.** From ~0.165 at t = 200 to ~0.03 at
t = 3000 is an e-folding of roughly **1600 sigma-time**, consistent with tau_r = 2090 — **the mode's
own damping — and not with the 210 predicted in section 1B of `261008` for the excess.** The 210
prediction is neither confirmed nor excluded here: **this observable is the wrong one to test it
with**, because the coherent mode is two orders of magnitude larger than the effect.

> **CORRECTION TO SECTION 4.2.** That section read T_1 - T_2 at four instants (t = 1070, 2120,
> 4210), found -0.0059, -0.0010 and +0.0004, and concluded "the signal is below the noise floor by
> 0.5 tau_r". **Those three times all land near zero crossings of a large oscillation.** The
> envelope at t = 1070 is not small — it is about 0.08, four times the noise floor. The correct
> statement is: **the coherent mode dominates T_1 - T_2 until its envelope damps on tau_r ~ 2090;
> what is below the noise floor is the AVERAGE over mechanical equilibrium (t > 10 tau_r), which
> is -0.0028 ± 0.0111 and averages the oscillation away.** The conclusion of 4.2 — that no
> exponential with tau = 40 079 is recoverable, and that no fit was performed — is unchanged, and
> the corrected reason is stronger: the record is dominated by a reversible mode, not by a thermal
> decay. **Sampling a ringing signal at four points and calling it small was the error.**

### 4.9 The mode, measured properly — period, energy, and the right reference work

*Post-hoc on B.*

**(a) Period from zero crossings, not from a spectrum.** A 2-sigma step leaks badly into a
periodogram, which is why section 4.8's residual spectrum read 173. Taking zero crossings of
x(t) - x_eq(t), with x_eq a running mean over one period, on the 80-seed mean over t = 30-3000:

| running-mean window | crossings | period |
|---|---|---|
| 195 | 32 | 186.1 ± 1.5 |
| **200** | **32** | **186.2 ± 1.4** |
| 210 | 32 | 186.2 ± 1.4 |

**Period = 186.2 ± 1.4**, insensitive to the window.

**The comparison with the ladder needs the temperature, which I left out.** The first version of
this section scaled the ladder period by the geometry and by c_s at the compressed eta **but at
kT = 1**. The gas is not at kT = 1 after the push: the same isentrope integral that gives
E_qs = 3.30 kT per gas heats both sides to **T_ad = 1.066069**, and c_s² = (kT/m)(Z + eta Z' + Z²)
carries a sqrt(T).

| scaling of the ladder's 204.81 ± 0.02 at alpha = 2 | c_s' | predicted period | vs 186.2 ± 1.4 |
|---|---|---|---|
| geometry only, c_s(eta') at kT = 1 *(first version — wrong)* | 1.770036 | 192.0 | 4.1 sigma |
| **geometry + isentropic heating, T_ad = 1.066069** | **1.827573** | **186.0** | **0.2 sigma** |
| + the 1.37 kT excess spread over 100 particles | 1.839279 | 184.8 | 1.0 sigma |

> **That is agreement, not a discrepancy. The 3 % I reported was the missing sqrt(T).**

**And it is a better statement than a period check.** The push changes **L, eta and T** of both
gases, and the mode's period follows all three — geometry, density and temperature — with no free
parameter. The ladder measured this mode at alpha = 2 in an unperturbed box; the same mode in a
compressed, heated box lands on the scaled prediction at 0.2 sigma.

**(b) Mode energy — and the inertia it appears to give is not a measurement.**

> **E_mode = W_qs(3.875) - 2 W_qs(1.9375) = 7.0715 - 6.6069 = 0.4646 kT**

That number stands. **What I withdraw is the M_eff = 217.4 inferred from it.**
E_mode = W_qs(d) - 2 W_qs(d/2) equals (1/4) k_1 d² only in the harmonic limit; the next term is
(1/8) k_2 d³, and with k_2/k_1 ~ (gamma+1)/L ~ 0.08 at d = 3.875 that is a **~15 % correction to
the inferred stiffness, in the direction that puts M_eff below 200**. So "the gas adds 17.4" was a
statement about the anharmonicity of the isentrope, not about inertia, and it is withdrawn.

**The harmonic-limit prediction does exist and is worth recording as untested.** For the standing
wave xi proportional to sin(kx)/sin(K), the gas kinetic energy per side is
(1/2) N_s m V² x [1/2 - sin 2K/(4K)]/sin²K, which at K = 0.6533 is **0.3535 N_s m per side**, i.e.
**M_eff = 235.4** [DERIVATION]; the same integral for the potential energy gives 1.004 x the static
stiffness, so the stiffness is the static one. That is what `cot K = alpha K` contains, and **it
cannot be tested from a 2-sigma swing — the frequency route in (a) already tests it, at 0.2 sigma.**

**(c) Section 4.5 has two reference works, and both belong there.**

| reference | value | what it is |
|---|---|---|
| **W_qs(3.875)** | **7.0715** | reversible work with the divider **held** during the push — the right reference while tau_push = 19.4 << period 186 |
| **2 W_qs(1.9375)** | **6.6069** | reversible **minimum** for the same final state |
| difference | **0.4646** | the mode energy, dissipated later on tau_r |

**Dissipation quoted against both:** W_in = 8.4403 ± 0.1840 is **+19.4 %** on 7.0715 and
**+27.7 %** on 6.6069. The first is the dissipation during the push; the second also counts the
mode energy that tau_r later turns into heat. Section 4.5's 19.4 % is the push figure and stands;
**27.7 % is the total irreversibility of the whole protocol.**

### 4.6 Not done

The two geometry pictures still need a GUI session (section 3.5). Not produced.
