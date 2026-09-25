# Level 4b — acoustic transmission scan: PRE-REGISTRATION

2026-10-08. **Written and committed before `level4b_transmission_20261008` wrote a single record.**
Campaigns A (R-collapse pass 2) and B (demonstration) are running or queued; this is chained behind
them and does not interrupt either.

The question: when the piston does work on gas 1, **what fraction of it crosses the divider on the
first pass**, and is the divider simply a mass on a string as far as that first pass is concerned?
This is the chain's requirement, and it is the mechanism behind `tau_T = M g(R)` seen at one
frequency instead of integrated over all of them.

---

## 1. PRE-REGISTRATION

### 1.1 Geometry, and three places the plan as written does not close

Two-compartment box, `--l0=39.25`, divider free, `wall_thickness_sigma = 1` (recorded, not assumed),
L_c = 38.75, eta_phys = 0.10134170, c_s = 1.749302 (Kolafa–Rottner), transit L_c/c_s = 22.15.

**(1) d = 3.875, not 3.876.** The plan asks for d such that u = 0.2 reproduces B1long's push, and
back-computes 3.876 from the rounded 19.38. **3.876 x 24 = 93.024 — not an integer, so it fails the
pixel grid.** d = **3.875** is exact (x 24 = 93) and gives tau_push = 19.375, which *is* B1long's
push. Used: **d = 3.875 for every cell**, so tau_push = d/u alone distinguishes them.

**(2) "N = 100 per gas" cannot be right, and the plan's own numbers say which reading is.** At
l0 = 39.25, 100 particles per side gives eta = 0.2027, twice the stated 0.10134. The consistent
reading is **N = 100 total, 50 per side** — the ladder box. Two independent checks confirm it: it
is the only reading that gives eta_phys = 0.10134170, and the plan's own Z ~ 2.3 comes out as
**Z = N m c_s/L_c = 50 x 1.749302/38.75 = 2.2572** with N = 50. It is also the only reading under
which B2pilot can be reused, since B2pilot is the ladder box. **N_s = 50 per side is used.**

**(3) f_early = |t|²/2, not |t|².** The plan gives the transmission as |t(w)|² and separately
requires (limit c) that x << 1 give f_early = **0.5**. Those are consistent only with a factor of
one half: at x -> 0 the divider is transparent *and massless*, so the work ends up shared equally,
not transferred entirely. The predicted observable is therefore

> **f_early = (1/2) / (1 + (w_eff M_d/2Z)²) = (1/2)/(1 + pi² x²)**, x = M_d/(2 Z tau_push).

The half is supplied here, before the data, and flagged as supplied.

**w_eff = pi/tau_push is an ansatz** [INFERENCE]: a step push of finite duration has a spectrum, not
a line, and pi/tau_push is a representative choice, not a derived one. If the collapse fails only in
its horizontal scale, this is the first thing to suspect — not the impedance identification.

**(4) tau_push is ~7 % longer than d/u, and the analysis uses the RECORDED value.**

> **LABEL, corrected 2026-10-08. This note is POST-DATA, PRE-ANALYSIS.** An earlier version of it
> said "still before Level 4b writes a record". **That was false.** The chain fired Level 4b at
> 12:16:08, forty-six seconds after campaign B finished, and the campaign was over at 12:16:54 —
> before this note was written. No 4b trace had been *opened*, and the note's content is derivation
> plus B's recorded values, with nothing fitted to a 4b number; but it cannot be called pre-data and
> the claim is withdrawn.

Two revisions, both from **recorded** values rather than inference:

**It is neither a trigger latency nor a travel overshoot — the piston runs slow.** The first version
read B1long alone (stop 20.650 vs nominal 19.375) and called it a fixed +1.275 sigma latency; the
second read B2pilot too (4.150 vs 3.875) and called it a 0.27 sigma travel overshoot. **Both were
inferences where a recorded number exists.** `piston_target_sigma = 74.625` against a start of 78.5
is a **travel of exactly 3.8750 in both cells** — the piston stops where it was told to. What is
long is the *time*:

| u | travel (recorded) | t_stop (recorded) | mean speed | as % of requested |
|---|---|---|---|---|
| 0.2 | **3.8750** | 20.6500 | 0.18765 | **93.8 %** |
| 1.0 | **3.8750** | 4.1500 | 0.93374 | **93.4 %** |

**The piston covers the requested distance at ~93.6 % of the requested speed.** The divider gives an
independent second reading and agrees: with D ~ 0 both gases end in identical states, so the divider
sits at half the travel, and 2 x 1.9058 ± 0.1447 = **3.81 ± 0.29**, which is 0.22 sigma from 3.875
and 1.2 sigma from 4.145.

**Consequence: W_qs is unchanged at 7.0715**, because the quasi-static work depends on the *travel*,
not on the time. The alternative in which travel were 4.145 and W_qs rose to ~7.6, halving the
quoted dissipation, **is excluded by the recorded target**. Section 4.5 of `261007` stands as
written.

**tau_push is taken from each cell's own recorded `piston_stop_t_rel`**, never from d/u. The nominal
table in section 1.2 is left unchanged beside it.

The causal margin in `261007` section 3.1 is correspondingly tighter than claimed there: the push
ends at **20.65** against the 22.15 transit, a margin of **1.5 sigma-time, not 2.8**.

### 1.2 The prediction, cell by cell

|t(w)|² = 1/(1 + (w M_d/2Z)²) is the mass-on-a-string transmission [STANDARD RESULT, citation unverified: Morse & Ingard,
*Theoretical Acoustics* — the section number was written from memory and nobody has opened the
book; `[SOURCE:]` in this project means read from the document, and this was not]. **Identifying Z with the gas impedance N m c_s/L_c is INFERENCE**,
and it is what the verdict rule tests.

| M_d | u | tau_push | w_eff | **x = M_d/(2 Z tau_push)** | \|t\|² | **predicted f_early** | ideal-gas Z | M* = 2 Z tau_push |
|---|---|---|---|---|---|---|---|---|
| 10 | 0.05 | 77.500 | 0.0405 | **0.0286** | 0.9920 | **0.4960** | 0.4939 | 349.9 |
| 10 | 0.2 | 19.375 | 0.1621 | **0.1143** | 0.8857 | **0.4429** | 0.4176 | 87.5 |
| 10 | 0.5 | 7.750 | 0.4054 | **0.2858** | 0.5536 | **0.2768** | 0.2238 | 35.0 |
| 10 | 1.0 | 3.875 | 0.8107 | **0.5717** | 0.2367 | **0.1183** | 0.0842 | 17.5 |
| 50 | 0.05 | 77.500 | 0.0405 | **0.1429** | 0.8322 | **0.4161** | 0.3821 | 349.9 |
| 50 | 0.2 | 19.375 | 0.1621 | **0.5717** | 0.2367 | **0.1183** | 0.0842 | 87.5 |
| 50 | 0.5 | 7.750 | 0.4054 | **1.4291** | 0.0473 | **0.0236** | 0.0157 | 35.0 |
| 50 | 1.0 | 3.875 | 0.8107 | **2.8583** | 0.0123 | **0.0061** | 0.0040 | 17.5 |
| 200 | 0.05 | 77.500 | 0.0405 | **0.5717** | 0.2367 | **0.1183** | 0.0842 | 349.9 |
| 200 | 0.2 | 19.375 | 0.1621 | **2.2866** | 0.0190 | **0.0095** | 0.0063 | 87.5 |
| 200 | 0.5 | 7.750 | 0.4054 | **5.7166** | 0.0031 | **0.0015** | 0.0010 | 35.0 |
| 200 | 1.0 | 3.875 | 0.8107 | **11.4331** | 0.0008 | **0.0004** | 0.0003 | 17.5 |

**x spans 0.029 to 11.4 — a factor 400 — with the crossover x = 1 crossed between cells.**

**(d) The ideal-gas discriminator.** Z_ideal = N m c_s,ideal/L_c = 1.8248, a ratio
c_s,KR/c_s,ideal = **1.2369**. The plan calls this weak and it is, *at the ends*: 0.4960 vs 0.4939 at
the smallest x. **Near the crossover it is not weak** — 0.2768 vs 0.2238 at (10, 0.5), an absolute
gap of 0.053 and a relative gap of 24 %. Reported as a discriminator whose power is concentrated in
three or four cells, not as a uniformly weak one.

### 1.3 Observable and verdict rule

**f_early = Delta E2/W_in** on the plateau **[tau_push + L_c/c_s + 2 tau_r, tau_push + L_c/c_s +
5 tau_r]**, with Delta E2 the *far* gas's energy gain (the piston is on the right, so gas 2 is the
left compartment) and W_in the **recorded** piston work. tau_r from `261006_mode_ladder.json`:
208 (M = 10), 586 (M = 50), 2090 (M = 200). Also reported: f at end of run, and the **divider
kinetic energy on the plateau — if it is not small the window is too early, and that is reported,
not fixed by shifting the window**.

**Verdict.** Collapse if, fitting the one-parameter curve with Z free and c_s fixed over all twelve
points, **chi²_red < 3 AND the fitted Z is within 25 % of 2.2572** (i.e. 1.693–2.822).
**Limits, either failure meaning the identification of Z is wrong and not the run:**
(10, 0.05) must give f_early = 0.5 within 2 sigma; (200, 1.0) must give f_early < 0.2.

### 1.4 Cells and cost

d = 3.875 throughout; 8 seeds per cell; run length tau_push + L_c/c_s + 5 tau_r(M_d), **not tau_T**
— the first-pass transmission is the target and tau_T(200) ~ 40 000 is deliberately out of reach.
`--trace-every=12` (dt = 0.2 sigma) for every cell: that is <= tau_r/20 by a factor 52 at the
tightest (M = 10) and still gives **19 samples across the shortest push** (tau_push = 3.875 at
u = 1.0), which a cadence chosen from tau_r/20 alone would not.

**The (200, 1.0) cell reuses B2pilot and is not rerun.** B2pilot is the same box, M_d = 200,
u = 1.0, d = 3.875, 8 seeds. Its `--trace-every=600` (dt = 10) cannot resolve the 3.875-sigma push
— but **f_early does not need the push resolved**: it needs W_in from the summary ledger and
Delta E2 on the plateau, which at dt = 10 carries 627 samples. The push-resolution panel of the
figure omits that one cell and says so.

**11 new cells, 8 seeds each, 28.0 M steps total = 0.12 core-hours.** This campaign is essentially
free. **Consequence to state in advance:** at 8 seeds the chi² test in 1.3 may be uninformative in
either direction — trivially passed if the seed errors are large, trivially failed if they are
tiny. The seed count is left at the pre-registered 8; **if the errors make chi² uninformative, the
honest move is to rerun at 32 or 80 seeds (0.5 or 1.2 core-hours) and say that the first pass was
under-powered, not to reinterpret the threshold.**

### 1.5 Figure

`261008_p2_level4b_transmission.{png,pdf}`: f_early against x for all twelve cells, **prediction as
a black curve**, our data **blue with error bars over seeds**, the ideal-gas curve dashed grey, the
two pre-registered limit boxes marked, and a second panel showing Delta E2(t)/W_in for the three
M_d at u = 0.2 with the plateau window shaded.

---

### 1A Amendment — the first-pass observable

> **PROVENANCE, stated plainly. I have no record of the message this amendment is said to come
> from.** There is no section 1A anywhere in this repository and no occurrence of `f_first`,
> `tau_ac` or "first kick" in any file. The text below is **reconstructed from the one-line summary
> given in the instruction**. It is **POST-DATA, PRE-ANALYSIS** — written after Level 4b finished
> at 12:16:54 and before any 4b trace was analysed — and its content is derivation plus B's
> recorded values, with nothing fitted to a 4b number. **If the original text differs, the original
> governs.**

**The observable changes from `f_early` to `f_first`, and the factor 1/2 goes away.** Section
1.1(3) supplied a half on the grounds that a transparent massless divider shares work equally. That
is true *after many passes*; on the **first pass** the transmitted fraction of an incident pulse is
the transmission itself. The window is bounded by the pulse's own return:

> the pulse is fully launched by tau_push and reaches the divider at L_c/c_s = 22.15; it reflects
> off the piston and returns **2 L_c/c_s = 44.30** later. So
> **f_first = Delta E_2/W_in averaged on [tau_push + 22.15, 66.45]**.

**tau_push is the COMPRESSION duration, d/u** — not `piston_stop_t_rel`. Section 4.7(a2) of
`261007` settles this from the code: the piston crosses a real 0.25-sigma gap doing exactly zero
work, so `piston_stop_t_rel` exceeds d/u by 0.25/u. Every x below uses d/u.

| M_d | u | tau_push = d/u | x = M_d/(2 Z tau_push) | first-written ansatz \|t\|² |
|---|---|---|---|---|
| 10 | 0.2 | 19.375 | 0.11433 | **0.8857** |
| 10 | 0.5 | 7.750 | 0.28582 | **0.5536** |
| 10 | 1.0 | 3.875 | 0.57165 | **0.2367** |
| 50 | 0.2 | 19.375 | 0.57165 | **0.2367** |
| 50 | 0.5 | 7.750 | 1.42912 | **0.0473** |
| 50 | 1.0 | 3.875 | 2.85824 | **0.0123** |
| 200 | 0.2 | 19.375 | 2.28659 | **0.0190** |
| 200 | 0.5 | 7.750 | 5.71647 | **0.0031** |
| 200 | 1.0 | 3.875 | 11.43295 | **0.0008** |

> **CORRECTED 2026-10-08. There was never an inconsistent limit — the flag was my arithmetic.**
> |t|^2 = 1/(1+pi^2 x^2) is monotone in x, and my (10, 0.2) entry read **0.4707** when the value at
> that x is **0.8857**. I then wrote "see note" against (10, 0.5) because the column was not
> monotone — **I saw the symptom and did not act on it** — and reported the limit as inconsistent
> with the table. One entry in nine was wrong; the other eight were right to four decimals.

**THE u = 0.05 CELLS CANNOT PRODUCE f_first, AND THIS IS STRUCTURAL.** The window needs
tau_push + 22.15 < 66.45, i.e. **tau_push < 44.30**. At u = 0.05 tau_push is 77.5, so the pulse
returns *before the push has finished*: there is no first pass to isolate. **Three of the twelve
cells are excluded from the f_first fit**, reported as excluded rather than as failures. They stay
in the W_in(u) ledger.

**Limits:** (10, 0.2) must give f_first > 0.6; (200, 1.0) must give f_first < 0.1.

**W_in(u) ledger.** W_in = W_qs + Z d u, i.e. **intercept 7.0715 (KR) or 5.5556 (ideal), slope
Z d = 8.746 (KR) or 7.071 (ideal)**. Fitted over u = 0.05, 0.2, 0.5 per M_d; u = 1.0 plotted only.
If u = 0.05 falls off the line at light M_d that is divider feedback, reported and not dropped.

**Ledger identity:** f_plateau = (1 - D)/2, with D from section 4.4 of `261007`.

**Abort gate:** the first-kick divider velocity in the heavy limit is ΔV = 2 Z d/M_d = **0.08747 at
M = 200**, *independent of u* — a slower push lasts proportionally longer, so the pulse momentum
Z d is unchanged. Measure it before any fit; **if ΔV < 0.04 the acoustic kick is absent and no fit
is performed.**

### 1B The work-sharing is one mechanism, not two

> **POST-DATA, PRE-ANALYSIS.** Rewritten 2026-10-08; the earlier two-mechanism version is
> superseded, and the reason is section 1C.

The first version offered "return-leakage" and "mode-parking" as competing mechanisms. **With the
spectral average they were never two mechanisms.** The divider's recoil and re-radiation *is* how
the low-frequency part of the pulse gets through — that is exactly what T(x) integrates. The
question is not *which* mechanism but

> **does the per-encounter transfer follow T(x), and does the excess decay on the timescale T(x)
> implies?**

**The M = 4 N_s m degeneracy, and why M = 200 is the wrong cell for the shape.** The divider's
acoustic radiation time is M/(2Z) and the pulse round trip is 2 L_c/c_s; their ratio is
**M/(4 N_s m)**, exactly 1 at M = 200, N_s = 50:

| M_d | radiation time M/(2Z) | round trip 2L_c/c_s | ratio |
|---|---|---|---|
| 10 | **2.22** | 44.30 | 0.050 |
| **50** | **11.08** | 44.30 | **0.250** |
| 200 | **44.30** | 44.30 | **1.000** |

At M = 200 recoil and return are **indistinguishable in time** — any curve fits both stories.
**The separable cell is (50, 0.2)**: x = 0.5717, radiation time 11.1 against a 44.3 round trip, a
factor four apart. **M = 200 tests the total, not the mechanism.**

**Predictions.**

- **(50, 0.2), the mechanism cell.** x = 0.5717, **T(x) = 0.5278**. The heavy-limit form
  KE_div = (1/2)M(2Zd/M)^2 does **not** apply at x = 0.57 — it is a large-x result. Use T(x):
  **energy into gas 2 after the first encounter = T x E_pulse**, and **KE_div at its peak is
  bounded above by that**. Recoil at 11.1 should be complete before the first return at 44.3.
- **(200, 0.2), the total.** Per-encounter transfer **T = 0.1900**, so the excess should equalise
  with an e-folding time of about **44.3/ln(1/(1-0.190)) = 210 sigma-time**
  [INFERENCE — encounters are not identical, the pulse spreads and the divider recoils].

**That 210 is testable on B1long, which is already on disk**: T_1 - T_2 was below the noise floor by
~1000 sigma-time, which is about five e-folds of 210. The test is section 4.8 of `261007`,
**labelled post-hoc on B**.

### 1C The spectral average replaces the ansatz

> **POST-DATA, PRE-ANALYSIS**, same label and reason. Derivation; nothing fitted to a 4b number.

**w_eff = pi/tau_push was an ansatz, and it underestimates the transmission by an order of magnitude
at large x.** A compression pulse whose velocity is everywhere positive has its energy spectrum
peaked at **w = 0**, not at pi/tau. The medium is linear and lossless, so each component transmits
|t(w)|^2 independently and the first-pass energy fraction is the **spectral average**

> **T = integral |V(w)|^2 |t(w)|^2 dw / integral |V(w)|^2 dw**, V = Fourier transform of the piston
> velocity.

For a rectangular velocity profile of duration tau_push, |V|^2 is proportional to
sinc^2(w tau/2), and with x = M_d/(2 Z tau_push) the integral closes:

> **T(x) = 1 - x (1 - e^(-1/x))**   [DERIVATION: partial fractions
> 1/(s^2(1+4x^2 s^2)) = 1/s^2 - 4x^2/(1+4x^2 s^2); integral sin^2 s/s^2 = pi;
> integral cos 2s/(1+b^2 s^2) = (pi/b) e^(-2/b), b = 2x]

**Limits: T -> 1 as x -> 0, and T -> 1/(2x) as x -> infinity** — a 1/x tail against the ansatz's
1/x^2. **Independent cross-check from momentum:** a heavy divider takes 2Zd from the pulse, so
KE_div = 2 Z^2 d^2/M = (1/x) E_pulse, and it radiates half into gas 2 -> **1/(2x)**. Two
derivations, one answer.

**The rectangular-profile condition is verified, not assumed.** `PistonR_v` in B1zoom is **-0.2
exactly at every sample while moving** (sd 5.6e-17, a single distinct value; a linear fit to x(t)
leaves a maximum residual of 5e-6 sigma). It is a clean step, so the closed form applies and no
numerical spectral average is needed.

| M_d | u | tau_push = d/u | x | ansatz \|t\|² | **T(x), the prediction** |
|---|---|---|---|---|---|
| 10 | 0.2 | 19.375 | 0.11433 | 0.8857 | **0.8857** |
| 10 | 0.5 | 7.750 | 0.28582 | 0.5536 | **0.7228** |
| 10 | 1.0 | 3.875 | 0.57165 | 0.2367 | **0.5278** |
| 50 | 0.2 | 19.375 | 0.57165 | 0.2367 | **0.5278** |
| 50 | 0.5 | 7.750 | 1.42912 | 0.0473 | **0.2808** |
| 50 | 1.0 | 3.875 | 2.85824 | 0.0123 | **0.1562** |
| 200 | 0.2 | 19.375 | 2.28659 | 0.0190 | **0.1900** |
| 200 | 0.5 | 7.750 | 5.71647 | 0.0031 | **0.0826** |
| 200 | 1.0 | 3.875 | 11.43295 | 0.0008 | **0.0425** |

**The verdict rule (Z free, chi²_red < 3, Z within 25 % of 2.2572) is applied to T(x).** The ansatz
curve stays in the figure, **dashed, labelled "first-written ansatz"**, as the record of what was
predicted first. **Limits: 0.8857 against > 0.6, and 0.0425 against < 0.1 — both comfortable.**

### 1D The observable is the EXCESS ledger, not Delta E_2/W_in

> **POST-DATA, PRE-ANALYSIS**, same label and reason as 1A–1C. Written before any 4b trace was
> analysed. This supersedes the definition of f_first in 1A.

**Why 1A's f_first would have killed T(x) for the wrong reason.** The same first-order term that
makes T_1 - T_2 swing in `261007` section 4.8 also dominates Delta E_2 inside the first-pass window.
A divider displaced by X does reversible work P h X on the far gas, and at M = 200 the divider has
already moved about 1.5 sigma by t = 66:

> **P h X = 0.15999 x 10 x 1.5 = 2.4 kT, about 30 % of W_in**, against an acoustic first-pass
> transfer of **T(x) x E_pulse = 0.190 x 1.64 = 0.31 kT** — a factor 8 larger, and of the wrong
> physical kind.

At M = 10 the divider is already past its new equilibrium inside the window, so Delta E_2/W_in
would read 0.5 or more. **The measured "f_first" would have come out ~0.3 at the heavy end and
~0.5 at the light end, T(x) would have been declared dead, and the impedance identification would
have taken the blame for an observable that was never the acoustic one.** This is the third
appearance of the same error class in this batch — a reversible, first-order term mistaken for the
quantity of interest — and it is caught before the fit rather than after.

**The fix is the one that already worked for the plateau: subtract the reversible part exactly.**
Each compartment's length L_i(t) is recorded in `SegEtas`, and the isentropic energy change from
eta_0 = 0.10134170 is the same KR integral used for W_qs:

> **E_qs,i(t) = W_qs(38.75 -> L_i(t))** (signed; negative on expansion)
> **X_i(t) = Delta KE_i(t) - E_qs,i(t)**

X_i is the **non-isentropic** energy in gas i — acoustic plus dissipated heat, and nothing else.
**The divider mode lives entirely inside the E_qs terms and cancels**, which is exactly what makes
it the right observable. Exact per-seed ledger, at every sample:

> **W_in = E_qs,1 + E_qs,2 + X_1 + X_2 + KE_div**

> **f_first = X_2 / (X_1 + X_2 + KE_div)**, reported as the window average on
> [tau_push + 22.15, 66.45] **and** as the value at t = 66.45 with (1/2) KE_div added to X_2
> (the divider is still radiating at that instant; section 1B).

**Prediction.** f_first = **T(x)** in the linear limit, and **below T(x)** by whatever fraction of
the excess has already thermalised in gas 1 before the pulse reaches the divider [INFERENCE —
expected small at u = 0.2, and large at u = 1.0 where the measured excess is 13.6 kT against
Z u d = 8.2 kT, i.e. most of it is not in the coherent pulse at all].

**Verdict rule** (Z free, chi²_red < 3, Z within 25 % of 2.2572) is applied to the **six cells with
u in {0.2, 0.5}**. **u = 1.0 is plotted and not fitted**, and is expected to sit below the curve;
u = 0.05 has no first pass (section 1A). **Limits unchanged: (10, 0.2) > 0.6 and (200, 1.0) < 0.1.**

**The naive Delta E_2/W_in is reported in the same window alongside**, so the record shows what the
first-order term would have done.

## 2. Results

`level4b_transmission_20261008b`, 12 cells, 8 seeds each except (200, 1.0) which has **32**
(added post-hoc 2026-10-08 and labelled as such in SPECS). 0 aborts, 0 health events. Every number
below was printed by `paper2_level4b_final_20261009.py`; none was transcribed by hand.

### 2.1 Validation — both items pass

**(a) The per-seed energy ledger closes.**

> max |W_in - (dKE_1 + dKE_2 + KE_div)| over **every sample, every seed, all 12 cells**
> = **1.03 x 10^-6 kT**

That is float32 rounding on terms that swing over +/- 11 kT, and it is not a tautology: 38 845 of
52 500 samples per seed carry a nonzero residual. `W0_v` is sigma/sigma-time and equipartition
confirms the mass (1/2 M <v^2> = 0.501, 0.567, 0.51 kT at M = 10, 50, 200).

**(c) The first kick matches the heavy-limit prediction at M = 200.**

| u | measured dV (coherent, jackknife) | 2Zd/M (KR) | 2Z_id d/M (ideal) | gate \|dV\| > 0.04 |
|---|---|---|---|---|
| 0.05 | -0.0379 ± 0.0080 | -0.0875 | -0.0707 | **FAIL** (no first-pass window anyway) |
| **0.2** | **-0.0839 ± 0.0113** | **-0.0875** | -0.0707 | **PASS** |
| 0.5 | -0.1099 ± 0.0106 | -0.0875 | -0.0707 | PASS |
| 1.0 | -0.1342 ± 0.0086 | -0.0875 | -0.0707 | PASS |

At u = 0.2 the measurement is **0.3 sigma from the KR value and 1.2 sigma from the ideal-gas one**.
The sign is negative because the pulse pushes the divider away from the piston. The abort gate
passes at every cell that enters the fit.

### 2.2 The work ledger — and why its intercept is not yet an EOS test

| M_d | W_in(0.05) | W_in(0.2) | W_in(0.5) | intercept (3 pt) | slope (3 pt) | intercept (2 pt) |
|---|---|---|---|---|---|---|
| 10 | 6.677 ± 0.134 | 8.165 ± 0.539 | 14.092 ± 1.127 | 5.928 ± 0.187 | 14.44 ± 2.11 | **6.181** |
| 50 | 7.040 ± 0.079 | 7.603 ± 0.550 | 14.092 ± 1.127 | 6.427 ± 0.139 | 11.94 ± 2.09 | **6.853** |
| 200 | 7.172 ± 0.126 | 7.943 ± 0.471 | 14.092 ± 1.127 | 6.549 ± 0.179 | 11.57 ± 2.02 | **6.915** |

Targets: KR **7.0715** / slope 8.747; ideal **5.5556** / slope 7.071.

**The three-point line is contaminated by curvature and its intercept should not be quoted against
KR.** At u = 0.5 the Mach number is 0.29 and W_in = W_qs + Z d u has visibly broken (14.09 against a
linear 11.5); fitting a curve with a straight line drags the intercept down and the slope up.
Dropping u = 0.5 moves M = 200 from 6.549 to **6.915**, within 2 % of 7.0715. **A usable intercept
needs u = 0.01-0.02, which this grid does not have** — that is carried into the efficiency map,
whose grid starts at u = 0.01.

**One clean result here.** W_in at u = 0.05 rises with mass — **6.677, 7.040, 7.172** — because a
light divider follows the push and gas 1 never receives the full compression. That is the divider
feedback section 1A said to report rather than drop, and it is measured.

### 2.3 f_first against T(x): MODEL TEST FAILED

| M | u | x | T(x) | ansatz | **f (ratio of means)** | naive ΔKE₂/W_in | <X₂> | <denominator> | E_pulse |
|---|---|---|---|---|---|---|---|---|---|
| 10 | 0.2 | 0.1143 | 0.8857 | 0.8857 | **-0.211 ± 0.861** | 0.897 | -0.134 ± 0.574 | **+0.634 ± 0.430** | 1.749 |
| 10 | 0.5 | 0.2858 | 0.7228 | 0.5536 | **+0.171 ± 0.287** | 0.598 | +1.100 ± 2.062 | +6.441 ± 1.194 | 4.373 |
| 50 | 0.2 | 0.5717 | 0.5278 | 0.2367 | **-1.125 ± 1.530** | 0.843 | -0.462 ± 0.518 | **+0.411 ± 0.409** | 1.749 |
| 50 | 0.5 | 1.4291 | 0.2807 | 0.0473 | **+0.076 ± 0.049** | 0.651 | +0.483 ± 0.343 | +6.328 ± 0.909 | 4.373 |
| 200 | 0.2 | 2.2866 | 0.1900 | 0.0190 | **-0.010 ± 0.122** | 0.359 | -0.012 ± 0.156 | **+1.246 ± 0.459** | 1.749 |
| 200 | 0.5 | 5.7166 | 0.0826 | 0.0031 | **+0.046 ± 0.022** | 0.347 | +0.331 ± 0.149 | +7.271 ± 1.096 | 4.373 |
| *10* | *1.0* | *0.5717* | *0.5278* | *0.2367* | *+0.248 ± 0.063* | *0.510* | | | *8.747* |
| *50* | *1.0* | *2.8583* | *0.1562* | *0.0123* | *+0.071 ± 0.043* | *0.577* | | | *8.747* |
| *200* | *1.0* | *11.4331* | *0.0425* | *0.0008* | *+0.049 ± 0.008* | *0.354* | | | *8.747* |

*(italic rows = u = 1.0, plotted not fitted. The u = 0.05 cells have no first pass: the window needs
tau_push < 44.30 and tau_push(0.05) = 77.5.)*

> **Six-cell pre-registered fit, Z free: Z = 0.525 ± 0.260 against 2.2572. chi2_red = 0.68.**
> **Z within 25 %: FAIL. chi2_red < 3: PASS.**
> **Limits: (10, 0.2) f = -0.211 against > 0.6 — FAIL. (200, 1.0) f = +0.049 against < 0.1 — PASS.**
>
> ## **MODEL TEST FAILED. Ledger and first kick pass.**

**The chi² passing means nothing here and should not be read as support.** Three of the six fitted
cells are the u = 0.2 row, whose **denominator is consistent with zero** — 0.634 ± 0.430,
0.411 ± 0.409, 1.246 ± 0.459 against an E_pulse of 1.749 — so their f carries errors of 0.86, 1.53
and 0.12 and a chi² fit through them is unconstrained. A negative transmission is unphysical and
says the observable, not the physics, has failed at those cells.

**Why the excess is not resolved at u = 0.2.** X_i = dKE_i - E_qs,i(L_i(t)) removes the *first-order*
reversible term, which is what section 1D was for, and it works: the naive ΔE₂/W_in column reads
0.36-0.90 across the grid with no dependence on x at all — **exactly the flat, meaningless answer
1D predicted the naive observable would give**. But E_qs assumes the gas follows the divider
*quasi-statically*, and at light M the divider swings on a period (95 at M = 10) comparable to the
gas transit time (22), so X_i still contains the non-quasi-static part of the divider's own
compression. That residual is second order and small in absolute terms — but at u = 0.2 the whole
acoustic pulse is only 1.75 kT, so small is not small enough.

**Where the denominator is solid (u = 0.5 and 1.0) f sits systematically below T(x)**, by about half
at moderate x and converging at large x (0.049 ± 0.008 against 0.043 at x = 11.4). That is the
direction section 1D pre-registered as an INFERENCE: at u = 1.0 the total excess is 15.6 kT against
Z u d = 8.75, so **44 % of it was never in the coherent pulse**, and f measures the share of the
total excess while T(x) predicts the share of the pulse. **No further amendment is made and no
further 4b run is done.** Closing T(x) properly needs the excess normalised to the coherent pulse,
which is a definition change and is not made after seeing the data.

### 2.4 D(u, M_d)

| M_d | u = 0.05 | 0.2 | 0.5 | 1.0 |
|---|---|---|---|---|
| 10 | +0.480 ± 0.307 | -0.268 ± 0.368 | +0.243 ± 0.121 | +0.065 ± 0.123 |
| 50 | +0.302 ± 0.347 | -0.448 ± 0.291 | +0.235 ± 0.141 | +0.114 ± 0.089 |
| 200 | +0.075 ± 0.226 | -0.082 ± 0.156 | +0.123 ± 0.111 | **+0.224 ± 0.049** |

**The 32-seed cell: D(200, 1.0) = +0.224 ± 0.049, i.e. 4.5 sigma from zero.** Full retention
predicts 0.347, so the measurement is **2.5 sigma below full retention and 4.5 sigma above zero —
partial sharing, not complete sharing.**

**This is not the same window as B2pilot and must not be merged with it.** B2pilot measured
+0.021 ± 0.099 over t > 10 tau_r on a 120 833-sigma record; this cell's record is 10 476 sigma, so
its plateau (t > 6 286) is about 3 tau_r, where the mode is damped but not gone. **The two numbers
answer different questions and the 8-seed-versus-32-seed comparison the pilot was for is
therefore not yet made.** What section 4.4 of `261007` claims on eight seeds still stands as
written, with its provenance.

### 2.5 Figure

`261008_p2_level4b_transmission.{png,pdf}` — (a) f_first against x with T(x), the first-written
ansatz dashed, and the naive ΔE₂/W_in in red as the trap; (b) the W_in(u) ledger against both
equations of state.
