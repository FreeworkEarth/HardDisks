# Paper 1 confinement campaign — PRE-REGISTRATION (not launched)

Written 2026-10-12 before any run. **Nothing here is launched.** *Amended the same day (§ 1.9, C1–C2); still unlaunched.* The campaign starts after the
Paper 2 efficiency map and only on an explicit go. Every number in the tables of § 1.4 is printed
by `python3 hspist3/validation/paper1_confinement_prereg_20261012.py`, pasted verbatim.

## 1. PRE-REGISTRATION

### 1.1 The question

Paper 1's $N = 100$ box measures $c_s$ above Kolafa–Rottner (KR) by **+1.0 %** at
$\eta \approx 0.10$ and **+1.675 %** at $\eta = \pi/8$. A2 showed the excess closes when the
square box grows (methods § 9), and Román's Table II shows the same at $\eta = \pi/8$ (260913
REPORT, 1a). But both series scale $H$ and $L$ together. The same straight line in $N^{-1/2}$ fits

$$\text{A:}\;\; \Delta = a\Big(\frac{2}{H} + \frac{2}{L_0}\Big), \qquad
\text{B:}\;\; \Delta = \frac{b}{H}\;\;(\text{weak } L), \qquad
\text{C:}\;\; \Delta = \frac{c'}{N_s},$$

where $\Delta \equiv c_s/c_s^{\rm KR}(\eta) - 1$. This campaign moves $H$, $L_0$ and (at fixed
area) the aspect ratio separately, so the three forms predict different numbers.

**Hypothesis C is added here; it is not in the plan.** The divider's thermal excursion makes the
mode slightly anharmonic. For a compartment force $F \propto L^{-q}$ with
$q = c_s^2/(Z\,kT/m) = (Z + \eta Z' + Z^2)/Z$, the restoring force of the divider is a hardening
Duffing force, $\propto x + (q+1)(q+2)x^3/(6L^2)$. A thermally driven Duffing oscillator has
$\langle x^2\rangle = kT/(2k_S)$, and its mean frequency shift
[STANDARD RESULT, citation unverified] gives

$$\Delta_C = \frac{(q+1)(q+2)}{16\,N_s\,q\,Z}.$$

That is +0.63 % at the $\eta \approx 0.10$ anchor and +0.38 % at $\pi/8$. It depends on $N_s$
only, not on $H$ or $L$ separately. It is an INFERENCE (heavy-divider limit, not tested), and it
matters for the design: **in the H-scan $N_s \propto H$, so B and C give identical predictions**
(Table P, last column, 0.0σ). Only the L-scan and the aspect scan separate them. C also predicts a
specific *violation* of the identity in § 1.5 (Table I), because the mode carries the shift and the
local $k_T$ does not.

### 1.2 Box and states

Román geometry, **gas | divider | gas**, $r = 0.5$, divider thickness $t = 0.05$ (Paper 1),
$L_{\rm eff} = L_0 - 2r - t/2$. $N_s$ per side is set by $\eta = N_s\pi\sigma^2/(4HL_0)$.

| anchor | $L_0$ | $H$ | $N_s$ | $\eta$ | why |
|---|---|---|---|---|---|
| $\eta \approx 0.10$ | 39.25 | 10 | 50 | 0.100051 | grid-exact (942/24). The same $L_0, H, N_s$ as Paper 2's Level 4 box (only $t$ differs: 0.05 vs 1.0). A2's $L_0 = 39.2699$ ($\eta = 0.100000$) is not grid-exact; the 0.05 % shift in $\eta$ moves KR's $c_s$ by far less than the error. |
| $\eta = \pi/8$ | 10 | 10 | 50 | 0.392699 | Assumption: the plan's "0.40" is read as Román's $\pi/8$. This is the A1v2 canonical cell, Román's own Table I box, and the $N = 50$ member of his Table II square series, so the anchor has published data at both sources. |

Every cell in Table G is grid-exact (1/24 σ). The aspect cells hold $N_s = 50$ and the anchor area
to the grid, so their $\eta$ moves by at most 0.2 %. KR is evaluated at each cell's exact $\eta$.

### 1.3 Scans, in launch order

1. **H-scan (first).** $H \in \{5, 10, 20, 40\}$ at the anchor $L_0$, $N_s = 5H$, both anchors.
2. **L-scan.** $L_0 \in \{\tfrac12, 1, 2\} \times$ anchor at $H = 10$, $N_s \propto L_0$, both anchors.
3. **Aspect control.** Fixed $N_s = 50$ and fixed area, $L_0/H \in \{1, 2, 4, 8\}$. It costs
   21 core-h for both anchors, and **it is the clean test of C** (C predicts a flat line at fixed
   $N_s$). It is also the strongest A-vs-B lever at $\pi/8$: −4.3σ and −7.0σ at $L_0/H = 4$ and 8.
   **Recommended**, not optional. At $\pi/8$, $L_0/H = 1$ *is* the anchor and is not rerun.

**Where the discriminating power is (Table P).** At $\eta \approx 0.10$ the H-scan separates A from
B by at most 0.7σ, because with $L_0 = 39.25 \gg H$ the $2/L_0$ term is small. **The decision
rests on the $\pi/8$ H-scan (−2.8σ / +2.1σ at $H = 5 / 40$), both L-scans, and the aspect scan.**
This is stated now so that a null result at $\eta \approx 0.10$ is not read as evidence.

### 1.4 Methods per cell

**(A) Held divider → $F(L)$, $k_T$, symmetry.** This runs in energy-transfer mode, the only mode
that writes the divider event log (`HD_PISTON_EVENTS`, `00ALLINONE.c:16651`). `edmd.c:1069` writes
`dp` = the **particle's** momentum change $m(v_{\rm after} - v_{\rm before})$ for every divider
event `D0`. For a held divider, a particle arriving from the left leaves with $dp < 0$ and one from
the right with $dp > 0$. So the **face is the sign of `dp`**, and

$$F_L = -\frac{1}{T}\sum_{dp<0} dp, \qquad F_R = \frac{1}{T}\sum_{dp>0} dp.$$

The reduction fails closed on any $dp = 0$ event. Level 3 measured $F$ this way from `|dp|` and
closed to $Z_{\rm box}/Z_{\rm KR} = 1.0248 \pm 0.0014$.

The divider is held with mass factor $10^9$ (Level 3 convention) at $L_0 + x$,
$x \in \{0, \pm\delta L, \pm 2\delta L\}$. Each run gives $F_L(L_0 + x)$ and $F_R(L_0 - x)$, so every
$L$-point is measured twice, by the two faces of the mirror runs $\pm x$. Then

$$k_T = -\frac{F(L_0 - 2\delta L) - 8F(L_0 - \delta L) + 8F(L_0 + \delta L) - F(L_0 + 2\delta L)}{12\,\delta L},$$

with $T = {\rm KE}/N$ measured on the same seeds (2D: $U = NkT$).

**$\delta L$ rule.** $\delta L = \max(1/24,\ \text{grid-rounded } \sigma_x/2)$, where
$\sigma_x = (kT L_{\rm eff}^2/2N_s m c_s^2)^{1/2}$ is the free divider's thermal rms excursion. The
stencil $\pm 2\delta L$ then spans what the free divider of (B) actually samples, so the identity
compares like with like. A finite box has solvation-force structure in $F(L)$ on the scale σ that
the KR model below does not contain; this rule averages over the same range the mode does.

**Budget, < 1 % on $k_T$.** Truncation is ≤ 0.1 %; the KR model gives ≤ $2\times10^{-5}$ in every
cell (Table A). The rest, ≤ 0.9 %, is noise, and the record per position is chosen to meet it.
**Noise model:** $\sigma_F/F = 0.2710\,(100/N_s)^{1/2}/\sqrt{T}$. The coefficient is measured from
Level 3 c0; the $N_s^{-1/2}$ scaling is ASSUMED, on the argument that wall-force fluctuations
follow ${\rm KE}_x$ fluctuations, not shot noise. **At $\pi/8$ the coefficient is unmeasured**, so a 4-seed pilot at the
anchor fixes it before (A) launches at $\pi/8$, and the record is recomputed by the same rule.
**Seeds are ≤ 5000 σ-time each.** Over the ~$10^6$ σ-time needed per position, a $10^9$ divider
would wander by its thermal amplitude. Per 5000-σ seed the drift is $\le 4\times10^{-3}\,\sigma$
(Table A): ≤ 3 % of $\delta L$ in every cell except $\pi/8$, $H = 40$, where $\delta L$ sits at the
1/24 floor and the drift is 9 % of it. The drift is random in sign, so it adds noise that averages
over the seeds, not bias.

**Symmetry checks**, each within 2σ: $F_L(L_0) = F_R(L_0)$ at $x = 0$; and, for every $x$,
$F_L(L_0 + x)$ from run $+x$ equals $F_R(L_0 + x)$ from run $-x$.

**(B) Free divider → $\nu_1$, $c_s$, $\Gamma$, power.** This runs in speed-of-sound mode exactly as
A1v2 (methods § 8): 25 seeds, 200 oscillations, `drift-first`, `--edmd-acc=0`, `HD_KE_TRACE=1`,
one invocation per (cell, $M$, seed) with `--speed-sound-exact-seed`. The **masses are scaled with
$N_s$**, $M = M_{\rm A1v2} \times N_s/50$, so every cell runs the same set
$\alpha = M/(2N_s m) = 0.5 \ldots 20$ and the same $K$ ladder ("the canonical masses" read as the
canonical $\alpha$ set; stated assumption). $c_s$ comes from the canonical estimator
(`paper1_populate_cs_err_20261002.cell`, TD = 200, X_EDGE = 2.5, through-origin
`slope_with_errors`) with $\nu = c_s K/2\pi L_{\rm eff}$ and $\cot K = \alpha K$. **Per (cell, $M$):**
$\Gamma = 2/\tau_r$ and $P_1 = B$ from the ACF fit, per methods § 13.


### 1.4b Tables (printed by `python3 hspist3/validation/paper1_confinement_prereg_20261012.py`, verbatim)

#### Inputs read from disk

- per-face force noise, Level 3 c0 (20 seeds, 867 sigma-time each, F = 1.6107): sigma_F/F = **0.2710/sqrt(T)** for N = 100 behind the face, H = 10, eta = 0.10005
- CPU cost, A1v2 run.log (N = 100): eta = 0.1122: **0.372 ms per sigma-time**
- CPU cost, A1v2 run.log (N = 100): eta = 0.392699: **2.060 ms per sigma-time**
- anchors (fractional c_s excess over KR at H = 10, N_s = 50): eta ~ 0.10: +1.01 %; eta = 0.392699: +1.675 % (260919 table)
- anchor relative error on c_s (scaled): eta ~ 0.10: 0.292 % (A1v2 eta = 0.1122 used as proxy); eta = 0.392699: 0.294 %

#### Table G -- cell geometry (gas | divider | gas; r = 0.5, t = 0.05; eta = N_s pi r^2 / (H L_0))

| scan | eta anchor | H | L_0 | N_s per side | eta exact | L_eff | L_0/H | grid-exact (1/24) | masses M (alpha = 0.5 ... 20) |
|---|---|---|---|---|---|---|---|---|---|
| H | 0.10 | 5.0000 | 39.2500 | 25 | 0.100051 | 38.2250 | 7.850 | yes | 25 ... 1000 |
| H | 0.10 | 10.0000 | 39.2500 | 50 | 0.100051 | 38.2250 | 3.925 | yes | 50 ... 2000 |
| H | 0.10 | 20.0000 | 39.2500 | 100 | 0.100051 | 38.2250 | 1.962 | yes | 100 ... 4000 |
| H | 0.10 | 40.0000 | 39.2500 | 200 | 0.100051 | 38.2250 | 0.981 | yes | 200 ... 8000 |
| L | 0.10 | 10.0000 | 19.6250 | 25 | 0.100051 | 18.6000 | 1.962 | yes | 25 ... 1000 |
| L | 0.10 | 10.0000 | 78.5000 | 100 | 0.100051 | 77.4750 | 7.850 | yes | 100 ... 4000 |
| aspect | 0.10 | 19.7917 | 19.7917 | 50 | 0.100252 | 18.7667 | 1.000 | yes | 50 ... 2000 |
| aspect | 0.10 | 14.0000 | 28.0000 | 50 | 0.100178 | 26.9750 | 2.000 | yes | 50 ... 2000 |
| aspect | 0.10 | 9.9167 | 39.6250 | 50 | 0.099937 | 38.6000 | 3.996 | yes | 50 ... 2000 |
| aspect | 0.10 | 7.0000 | 56.0417 | 50 | 0.100104 | 55.0167 | 8.006 | yes | 50 ... 2000 |
| H | 0.39 | 5.0000 | 10.0000 | 25 | 0.392699 | 8.9750 | 2.000 | yes | 25 ... 1000 |
| H | 0.39 | 10.0000 | 10.0000 | 50 | 0.392699 | 8.9750 | 1.000 | yes | 50 ... 2000 |
| H | 0.39 | 20.0000 | 10.0000 | 100 | 0.392699 | 8.9750 | 0.500 | yes | 100 ... 4000 |
| H | 0.39 | 40.0000 | 10.0000 | 200 | 0.392699 | 8.9750 | 0.250 | yes | 200 ... 8000 |
| L | 0.39 | 10.0000 | 5.0000 | 25 | 0.392699 | 3.9750 | 0.500 | yes | 25 ... 1000 |
| L | 0.39 | 10.0000 | 20.0000 | 100 | 0.392699 | 18.9750 | 2.000 | yes | 100 ... 4000 |
| aspect | 0.39 | 10.0000 | 10.0000 | 50 | 0.392699 | 8.9750 | 1.000 | yes | 50 ... 2000 |
| aspect | 0.39 | 7.0833 | 14.1250 | 50 | 0.392495 | 13.1000 | 1.994 | yes | 50 ... 2000 |
| aspect | 0.39 | 5.0000 | 20.0000 | 50 | 0.392699 | 18.9750 | 4.000 | yes | 50 ... 2000 |
| aspect | 0.39 | 3.5417 | 28.2917 | 50 | 0.391917 | 27.2667 | 7.988 | yes | 50 ... 2000 |

#### Table P -- predicted fractional c_s excess over KR, by hypothesis

A: Delta = a (2/H + 2/L_0), a fixed by the H = 10 anchor.  B: Delta = Delta_10 (10/H), no L dependence.
C (INFERENCE, heavy-divider estimate, no free parameter): thermal-amplitude anharmonicity,
   Delta_C = (q+1)(q+2) / (16 N_s q Z),  q = c_s^2/(Z kT/m) = (Z + eta Z' + Z^2)/Z.

| scan | eta | H | L_0 | N_s | A [%] | B [%] | C [%] | (A-B)/sigma | (B-C_scaled)/sigma | sigma assumed [%] |
|---|---|---|---|---|---|---|---|---|---|---|
| H | 0.100051 | 5.000 | 39.250 | 25 | +1.815 | +2.020 | +1.269 | -0.7 | +0.0 | 0.292 |
| H | 0.100051 | 10.000 | 39.250 | 50 | +1.010 | +1.010 | +0.634 | +0.0 | +0.0 | 0.292 |
| H | 0.100051 | 20.000 | 39.250 | 100 | +0.608 | +0.505 | +0.317 | +0.4 | +0.0 | 0.292 |
| H | 0.100051 | 40.000 | 39.250 | 200 | +0.406 | +0.252 | +0.159 | +0.5 | +0.0 | 0.292 |
| L | 0.100051 | 10.000 | 19.625 | 25 | +1.215 | +1.010 | +1.269 | +0.7 | -3.5 | 0.292 |
| L | 0.100051 | 10.000 | 78.500 | 100 | +0.907 | +1.010 | +0.317 | -0.4 | +1.7 | 0.292 |
| aspect | 0.100252 | 19.792 | 19.792 | 50 | +0.813 | +0.510 | +0.634 | +1.0 | -1.7 | 0.292 |
| aspect | 0.100178 | 14.000 | 28.000 | 50 | +0.862 | +0.721 | +0.634 | +0.5 | -1.0 | 0.292 |
| aspect | 0.099937 | 9.917 | 39.625 | 50 | +1.015 | +1.018 | +0.634 | -0.0 | +0.0 | 0.292 |
| aspect | 0.100104 | 7.000 | 56.042 | 50 | +1.294 | +1.443 | +0.634 | -0.5 | +1.5 | 0.292 |
| H | 0.392699 | 5.000 | 10.000 | 25 | +2.513 | +3.350 | +0.768 | -2.8 | +0.0 | 0.294 |
| H | 0.392699 | 10.000 | 10.000 | 50 | +1.675 | +1.675 | +0.384 | +0.0 | +0.0 | 0.294 |
| H | 0.392699 | 20.000 | 10.000 | 100 | +1.256 | +0.838 | +0.192 | +1.4 | +0.0 | 0.294 |
| H | 0.392699 | 40.000 | 10.000 | 200 | +1.047 | +0.419 | +0.096 | +2.1 | +0.0 | 0.294 |
| L | 0.392699 | 10.000 | 5.000 | 25 | +2.513 | +1.675 | +0.768 | +2.8 | -5.7 | 0.294 |
| L | 0.392699 | 10.000 | 20.000 | 100 | +1.256 | +1.675 | +0.192 | -1.4 | +2.8 | 0.294 |
| aspect | 0.392699 | 10.000 | 10.000 | 50 | +1.675 | +1.675 | +0.384 | +0.0 | +0.0 | 0.294 |
| aspect | 0.392495 | 7.083 | 14.125 | 50 | +1.775 | +2.365 | +0.384 | -2.0 | +2.3 | 0.294 |
| aspect | 0.392699 | 5.000 | 20.000 | 50 | +2.094 | +3.350 | +0.384 | -4.3 | +5.7 | 0.294 |
| aspect | 0.391917 | 3.542 | 28.292 | 50 | +2.661 | +4.729 | +0.384 | -7.0 | +10.4 | 0.294 |

#### Table K -- bulk thermodynamic ratios at the two anchors (KR)

| eta | Z | eta Z' | c_s^2 | q = c_s^2/Z | gamma = 1 + Z^2/(Z+eta Z') | k_T/k_S (bulk) |
|---|---|---|---|---|---|---|
| 0.100051 | 1.23628 | 0.27803 | 3.04271 | 2.4612 | 2.00930 | 0.49769 |
| 0.392699 | 2.76012 | 3.65471 | 14.03312 | 5.0842 | 2.18760 | 0.45712 |

#### Table I -- the identity test per cell: expected sigma, and what hypothesis C would do to it

rho_I = [N_s m c_s^2/L_eff^2 - k_T - F^2/(N_s kT)] / (N_s m c_s^2/L_eff^2), predicted 0.
sigma(rho_I)^2 = (2 sigma_cs)^2 + ((k_T/k_S) sigma_kT)^2, with sigma_cs = the anchor's relative error and
sigma_kT = 0.9 % (Table A budget, noise part); the F^2 term's own error (< 0.05 %) is neglected.
Under C the mode frequency carries the thermal-amplitude shift and the local k_T does not: rho_I ~ 2 Delta_C.

| scan | eta | H | L_0 | N_s | k_T/k_S (bulk) | sigma(rho_I) [%] | rho_I under C [%] | rho_I(C)/sigma |
|---|---|---|---|---|---|---|---|---|
| H | 0.1001 | 5.000 | 39.250 | 25 | 0.4977 | 0.736 | +2.537 | +3.4 |
| H | 0.1001 | 10.000 | 39.250 | 50 | 0.4977 | 0.736 | +1.269 | +1.7 |
| H | 0.1001 | 20.000 | 39.250 | 100 | 0.4977 | 0.736 | +0.634 | +0.9 |
| H | 0.1001 | 40.000 | 39.250 | 200 | 0.4977 | 0.736 | +0.317 | +0.4 |
| L | 0.1001 | 10.000 | 19.625 | 25 | 0.4977 | 0.736 | +2.537 | +3.4 |
| L | 0.1001 | 10.000 | 78.500 | 100 | 0.4977 | 0.736 | +0.634 | +0.9 |
| aspect | 0.1003 | 19.792 | 19.792 | 50 | 0.4977 | 0.736 | +1.268 | +1.7 |
| aspect | 0.1002 | 14.000 | 28.000 | 50 | 0.4977 | 0.736 | +1.268 | +1.7 |
| aspect | 0.0999 | 9.917 | 39.625 | 50 | 0.4977 | 0.736 | +1.269 | +1.7 |
| aspect | 0.1001 | 7.000 | 56.042 | 50 | 0.4977 | 0.736 | +1.269 | +1.7 |
| H | 0.3927 | 5.000 | 10.000 | 25 | 0.4571 | 0.717 | +1.536 | +2.1 |
| H | 0.3927 | 10.000 | 10.000 | 50 | 0.4571 | 0.717 | +0.768 | +1.1 |
| H | 0.3927 | 20.000 | 10.000 | 100 | 0.4571 | 0.717 | +0.384 | +0.5 |
| H | 0.3927 | 40.000 | 10.000 | 200 | 0.4571 | 0.717 | +0.192 | +0.3 |
| L | 0.3927 | 10.000 | 5.000 | 25 | 0.4571 | 0.717 | +1.536 | +2.1 |
| L | 0.3927 | 10.000 | 20.000 | 100 | 0.4571 | 0.717 | +0.384 | +0.5 |
| aspect | 0.3927 | 10.000 | 10.000 | 50 | 0.4571 | 0.717 | +0.768 | +1.1 |
| aspect | 0.3925 | 7.083 | 14.125 | 50 | 0.4572 | 0.717 | +0.768 | +1.1 |
| aspect | 0.3927 | 5.000 | 20.000 | 50 | 0.4571 | 0.717 | +0.768 | +1.1 |
| aspect | 0.3919 | 3.542 | 28.292 | 50 | 0.4573 | 0.718 | +0.769 | +1.1 |

#### Table A -- held divider (method A): stencil, derivative budget, record, cost

Rule: delta_L = max(1/24, grid-rounded sigma_x/2), sigma_x = (kT L_eff^2/(2 N_s m c_s^2))^(1/2) the free divider's
thermal rms excursion, so the stencil L_0 + {0, +-dL, +-2dL} spans what the free divider samples.
k_T = -[F(-2) - 8F(-1) + 8F(+1) - F(+2)]/(12 dL); noise factor sqrt(130)/12; per L-point the two faces of the
mirror runs (+x, -x) are averaged. Noise model eps(N_s) = eps0 (100/N_s)^(1/2) (ASSUMED; calibrated at eta = 0.39
by a 4-seed pilot before launch). Seeds of <= 5000 sigma-time, held mass 1e+09.

| scan | eta | H | L_0 | N_s | sigma_x | delta_L | 2dL/L_0 | bias (KR model) | T per position needed | seeds/position | divider drift/seed | CPU [core-h] |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| H | 0.1001 | 5.000 | 39.250 | 25 | 3.099 | 1.5417 (37/24) | 0.0786 | -2.0e-05 | 7.07e+05 | 142 | 4.6e-04 | 0.16 |
| H | 0.1001 | 10.000 | 39.250 | 50 | 2.191 | 1.0833 (26/24) | 0.0552 | -4.9e-06 | 7.16e+05 | 144 | 6.5e-04 | 0.33 |
| H | 0.1001 | 20.000 | 39.250 | 100 | 1.550 | 0.7917 (19/24) | 0.0403 | -1.4e-06 | 6.7e+05 | 135 | 9.2e-04 | 0.62 |
| H | 0.1001 | 40.000 | 39.250 | 200 | 1.096 | 0.5417 (13/24) | 0.0276 | -3.1e-07 | 7.16e+05 | 144 | 1.3e-03 | 1.33 |
| L | 0.1001 | 10.000 | 19.625 | 25 | 1.508 | 0.7500 (18/24) | 0.0764 | -1.8e-05 | 7.47e+05 | 150 | 6.5e-04 | 0.17 |
| L | 0.1001 | 10.000 | 78.500 | 100 | 3.141 | 1.5833 (38/24) | 0.0403 | -1.4e-06 | 6.7e+05 | 135 | 6.5e-04 | 0.62 |
| aspect | 0.1003 | 19.792 | 19.792 | 50 | 1.075 | 0.5417 (13/24) | 0.0547 | -4.8e-06 | 7.28e+05 | 146 | 9.1e-04 | 0.34 |
| aspect | 0.1002 | 14.000 | 28.000 | 50 | 1.546 | 0.7917 (19/24) | 0.0565 | -5.4e-06 | 6.82e+05 | 137 | 7.7e-04 | 0.32 |
| aspect | 0.0999 | 9.917 | 39.625 | 50 | 2.213 | 1.1250 (27/24) | 0.0568 | -5.5e-06 | 6.77e+05 | 136 | 6.4e-04 | 0.31 |
| aspect | 0.1001 | 7.000 | 56.042 | 50 | 3.154 | 1.5833 (38/24) | 0.0565 | -5.4e-06 | 6.83e+05 | 137 | 5.4e-04 | 0.32 |
| H | 0.3927 | 5.000 | 10.000 | 25 | 0.339 | 0.1667 (4/24) | 0.0333 | -4.9e-06 | 1.09e+06 | 219 | 1.4e-03 | 1.57 |
| H | 0.3927 | 10.000 | 10.000 | 50 | 0.240 | 0.1250 (3/24) | 0.0250 | -1.6e-06 | 9.7e+05 | 194 | 1.9e-03 | 2.78 |
| H | 0.3927 | 20.000 | 10.000 | 100 | 0.169 | 0.0833 (2/24) | 0.0167 | -3.1e-07 | 1.09e+06 | 219 | 2.7e-03 | 6.27 |
| H | 0.3927 | 40.000 | 10.000 | 200 | 0.120 | 0.0417 (1/24) | 0.0083 | -1.9e-08 | 2.18e+06 | 437 | 3.8e-03 | 25.01 |
| L | 0.3927 | 10.000 | 5.000 | 25 | 0.150 | 0.0833 (2/24) | 0.0333 | -4.9e-06 | 1.09e+06 | 219 | 1.9e-03 | 1.57 |
| L | 0.3927 | 10.000 | 20.000 | 100 | 0.358 | 0.1667 (4/24) | 0.0167 | -3.1e-07 | 1.09e+06 | 219 | 1.9e-03 | 6.27 |
| aspect | 0.3927 | 10.000 | 10.000 | 50 | 0.240 | 0.1250 (3/24) | 0.0250 | -1.6e-06 | 9.7e+05 | 194 | 1.9e-03 | 2.78 |
| aspect | 0.3925 | 7.083 | 14.125 | 50 | 0.350 | 0.1667 (4/24) | 0.0236 | -1.2e-06 | 1.09e+06 | 218 | 1.6e-03 | 3.12 |
| aspect | 0.3927 | 5.000 | 20.000 | 50 | 0.507 | 0.2500 (6/24) | 0.0250 | -1.6e-06 | 9.7e+05 | 194 | 1.4e-03 | 2.78 |
| aspect | 0.3919 | 3.542 | 28.292 | 50 | 0.730 | 0.3750 (9/24) | 0.0265 | -2.0e-06 | 8.66e+05 | 174 | 1.1e-03 | 2.48 |

#### Table B -- free divider (method B): 9 masses x 25 seeds x 200 periods, cost

| scan | eta | H | L_0 | N_s | period range (alpha 0.5 ... 20) | CPU [core-h] | wall [h] at 9 jobs |
|---|---|---|---|---|---|---|---|
| H | 0.1001 | 5.000 | 39.250 | 25 | 127.9 ... 620.9 | 0.71 | 0.08 |
| H | 0.1001 | 10.000 | 39.250 | 50 | 127.9 ... 620.9 | 1.41 | 0.16 |
| H | 0.1001 | 20.000 | 39.250 | 100 | 127.9 ... 620.9 | 2.83 | 0.31 |
| H | 0.1001 | 40.000 | 39.250 | 200 | 127.9 ... 620.9 | 5.65 | 0.63 |
| L | 0.1001 | 10.000 | 19.625 | 25 | 62.2 ... 302.1 | 0.34 | 0.04 |
| L | 0.1001 | 10.000 | 78.500 | 100 | 259.1 ... 1258.4 | 5.73 | 0.64 |
| aspect | 0.1003 | 19.792 | 19.792 | 50 | 62.7 ... 304.7 | 0.69 | 0.08 |
| aspect | 0.1002 | 14.000 | 28.000 | 50 | 90.2 ... 438.0 | 1.00 | 0.11 |
| aspect | 0.0999 | 9.917 | 39.625 | 50 | 129.1 ... 627.1 | 1.43 | 0.16 |
| aspect | 0.1001 | 7.000 | 56.042 | 50 | 184.0 ... 893.5 | 2.03 | 0.23 |
| H | 0.3927 | 5.000 | 10.000 | 25 | 14.0 ... 67.9 | 0.48 | 0.05 |
| H | 0.3927 | 10.000 | 10.000 | 50 | 14.0 ... 67.9 | 0.96 | 0.11 |
| H | 0.3927 | 20.000 | 10.000 | 100 | 14.0 ... 67.9 | 1.91 | 0.21 |
| H | 0.3927 | 40.000 | 10.000 | 200 | 14.0 ... 67.9 | 3.83 | 0.43 |
| L | 0.3927 | 10.000 | 5.000 | 25 | 6.2 ... 30.1 | 0.21 | 0.02 |
| L | 0.3927 | 10.000 | 20.000 | 100 | 29.6 ... 143.5 | 4.04 | 0.45 |
| aspect | 0.3927 | 10.000 | 10.000 | 50 | 14.0 ... 67.9 | 0.96 | 0.11 |
| aspect | 0.3925 | 7.083 | 14.125 | 50 | 20.4 ... 99.1 | 1.40 | 0.16 |
| aspect | 0.3927 | 5.000 | 20.000 | 50 | 29.6 ... 143.5 | 2.02 | 0.22 |
| aspect | 0.3919 | 3.542 | 28.292 | 50 | 42.6 ... 206.7 | 2.90 | 0.32 |

#### Cost per scan (A + B; the H = 10 anchor cell is counted once, in the H-scan)

| scan | eta anchor | cells | CPU A [core-h] | CPU B [core-h] | total [core-h] | wall [h] at 9 jobs |
|---|---|---|---|---|---|---|
| H | 0.10 | 4 | 2.5 | 10.6 | 13.0 | 1.4 |
| H | 0.39 | 4 | 35.6 | 7.2 | 42.8 | 4.8 |
| L | 0.10 | 2 | 0.8 | 6.1 | 6.9 | 0.8 |
| L | 0.39 | 2 | 7.8 | 4.3 | 12.1 | 1.3 |
| aspect | 0.10 | 4 | 1.3 | 5.2 | 6.4 | 0.7 |
| aspect | 0.39 | 3 | 8.4 | 6.3 | 14.7 | 1.6 |

### 1.5 Predictions, written now

**The identity, per compartment** (exact for hard disks, because $F \propto T$ at fixed $L$ and
$C_L = N_s k$):

$$\frac{N_s m\,c_s^2}{L_{\rm eff}^2} = k_S = -\Big(\frac{\partial F}{\partial L}\Big)_T + \frac{F^2}{N_s kT}.$$

Here $c_s$ comes from (B) and $k_T$, $F$, $T$ from (A); $L_{\rm eff}$ is the length in
$\nu = c_sK/2\pi L_{\rm eff}$. **Verdict: agreement within 2σ at every cell.** *(Superseded as the test by § 1.9, C1: the length-free form $2k_S^{\rm dyn} = \hat M\omega_1^2$ is primary.)* The expected
$\sigma(\rho_I)$ is 0.72–0.74 % (Table I). Under C the residual would be $\rho_I \approx 2\Delta_C$,
which is +3.4σ at $N_s = 25$, $\eta \approx 0.10$ and +1.7σ at the anchor. **So a failure
concentrated at small $N_s$ is C's signature, and is read that way.**

**$\gamma_{\rm box}$.** $\gamma_{\rm box} = k_S/k_T$ is compared with the bulk
$1 + Z^2/(Z + \eta Z') = 2.00930$ ($\eta = 0.100051$) and 2.18760 ($\pi/8$). It is also compared
with $1 + F^2/(N_s kT\,k_T)$ from (A) alone, which is the identity restated. This is reported with
its scaling in $H$ and $L$, with no pass/fail: a finite box need not reproduce the bulk ratio.

**Confinement.** Table P gives each hypothesis's prediction with its amplitude fixed by the
$H = 10$ anchor. B gives **+2.0, +1.0, +0.5, +0.25 %** at $H = 5 \ldots 40$ for
$\eta \approx 0.10$, as the plan states.

**Decision rule, per $\eta$, run once on the campaign's own cells.** For each of A, B and C, fit
$\Delta_i$ over all cells at that $\eta$ with **one free amplitude** (so the test is of the
*shape*), weighted by the scaled $\sigma_i$:

- A hypothesis is **excluded** if $p(\chi^2) < 0.01$.
- If exactly one survives, it is the result.
- If several survive, the outcome is "not separated", and $\Delta\chi^2$ is reported.
- If none survive, the two-term forms $b/H + c'/N_s$ and $a(2/H + 2/L_0) + c'/N_s$ are reported as
  exploratory, not as a verdict.

C is additionally reported with its amplitude fixed at $\Delta_C$ (no free parameter).

**Binary rule.** The campaign runs on the post-flag binary (v1 + `-ffp-contract=off` +
`--version`). It **re-measures its own anchors** (the $H = 10$ cells), and every fit uses campaign
cells only. A1v2 (contraction on) and A2 numbers appear only as quoted comparisons, never in a fit
or on a shared axis with campaign data.

### 1.6 Cost and order

Per Table "Cost per scan": H-scan 13.0 + 42.8 core-h (η ≈ 0.10, π/8), L-scan 6.9 + 12.1, aspect
6.4 + 14.7. **≈ 96 core-h in total, ≈ 11 h wall at 9 jobs.** The $\pi/8$ (A) runs dominate,
especially $H = 40$ (25 core-h) where $\delta L$ hits the 1/24 floor. Order: H-scan (π/8 pilot
first, then (B) and (A) both anchors), then L-scan, then aspect.

### 1.7 Gates before launch

1. **The go** (Chris and the plan author), after the Paper 2 map.
2. **Two pictures per cell geometry** (GUI + paper render) via `watch.sh`. That needs a `conf`
   entry in `watch.sh`, which is not written yet. There are 20 cells in 19 distinct geometries,
   plus the four held offsets of one (A) cell as a spot check.
3. **Geometry by code.** Energy-transfer mode must honour `--wall-thickness=0.05`, and the
   `--l0` / `--wall-positions` convention (`--l0` = one compartment, divider at
   $x = L_0 + x_{\rm off}$) must be checked from the summary CSV and the pictures, not assumed.
   Level 3 and Level 4 ran $t = 1.0$.
4. **Pilot at $\pi/8$:** 4 seeds per position at the anchor. This measures the noise coefficient,
   fixes the record, and checks that the steps-to-σ-time conversion gives 5000 σ-time per seed
   (read from the event-log time range, not computed from `--steps`).
5. **Mode-equivalence gate** (§ 1.9, C2): one cell in both modes; $\eta$, $t$, compartment lengths and $L_{\rm eff}$ agree to $10^{-6}$; pictures from both.
6. **Health contract** zero on every run (forced_advance, clamp_repair, overlap_repair,
   wall_overdue). `00_COMMAND.md` per leaf, and a `--version` line in every summary.

### 1.8 What would change this registration

Only the pilot of gate 4: it may change the record length and the cost of (A) at $\pi/8$, by the
rule already written. Nothing else is tuned after data.

---

### 1.9 Amendments C1–C2 (2026-10-12), before any run — the campaign stays unlaunched until the map is analysed and the go is given

#### C1 — the identity in its length-free form is the primary test

$$2\,k_S^{\rm dyn} \equiv \hat M\,\omega_1^2,\qquad \hat M = M + \tfrac{2}{3}N_s m \qquad(\text{heavy masses, }\alpha \ge 5),$$

compared with the static side, which is unchanged:

$$k_S^{\rm dyn} \overset{?}{=} -\Big(\frac{\partial F}{\partial L}\Big)_T + \frac{F^2}{N_s kT}.$$

$k_S^{\rm dyn}$ is computed per heavy mass, and the five values ($\alpha = 5, 7.5, 10, 15, 20$) are combined by inverse-variance weighting. Neither side contains a length. $\omega_1$ and $M$ are measured or set; $F$ and $\partial F/\partial L$ come from method A, where the derivative is with respect to the divider position, so no convention enters.

**$\hat M$ is $M + \tfrac23 N_s m$, not $M + \tfrac13 N_s m$.** The divider drives two gas columns, and each has the linear-profile inertia $N_s m/3$. Mansour's $\hat M = M + mN/3$ (Eq. 18) has $N = 2N_s$ total, and it is the $K \to 0$ limit of the standing-wave mass $M + 2N_s m[\tfrac12 - \sin 2K/4K]/\sin^2K$. Expanding $\cot K = \alpha K$ gives

$$K^2(\alpha + \tfrac13) = 1 - \frac{K^4}{45} + \dots\;\Rightarrow\;\omega^2 = \frac{2k_S}{M + \tfrac23 N_s m}\Big(1 - \frac{K^4}{45}\Big).$$

Table H shows the result. With $\tfrac23$, the heavy form matches the exact standing wave to ≤ 0.08 % for $\alpha \ge 5$. With $\tfrac13$, it would be off by 3.3 % at $\alpha = 5$ and still 0.8 % at $\alpha = 20$, several times the expected σ. The amendment as written in the plan would therefore have built a 1–3 % bias into the test, so $\tfrac23$ is used.

**The check at all $\alpha$** uses the exact standing-wave stiffness, which is also length-free:

$$k_S^{\rm SW} = \frac{N_s m\,\omega_1^2}{K(\alpha)^2},\qquad \cot K = \alpha K.$$

It is reported per mass, with no verdict.

**The $c_s/L$ form is NOT the test.** $N_s m c_s^2/L^2$ needs a length, and the choice between the geometric $L_0$ and the $L_{\rm eff}$ of the frequency formula moves it by $(L_0/L_{\rm eff})^2$. That is **1.0543 at the $\eta \approx 0.10$ anchor and 1.2415 at $\pi/8$** (Table L), far beyond any error bar. The § 1.5 identity in the $c_s/L$ form and Table I are therefore superseded as the test. The comparison of $c_s(\eta)$ with KR, using $L_{\rm eff}$, remains a separate bulk comparison, as in Paper 1.

**Expected σ, recomputed (Table I-ω).** $\sigma(\rho_I)$ is **0.457 % ($\eta \approx 0.10$) and 0.435 % ($\pi/8$)**, against 0.736 % and 0.717 % for the $c_s/L$ form. The heavy masses pin $\omega_1$ to 0.09–0.14 % on $k_S$, so the static side's 0.9 % on $k_T$, weighted by $k_T/k_S \approx 0.5$, now dominates. Under hypothesis C, $\rho_I \approx 2\Delta_C$ is **+2.8σ and +1.8σ** at the two anchors, up from +1.7σ and +1.1σ. The ω-form makes C easier to test, not harder.

**Verdict rule, unchanged:** agreement within 2σ at every cell.

#### C2 — mode-equivalence gate (added to § 1.7)

Before any quantity from method A (energy-transfer mode, held divider) is compared with any quantity from method B (speed-of-sound mode, free divider), one test cell is run in both modes. The recorded $\eta$, divider thickness, both compartment lengths and $L_{\rm eff}$ must agree to $10^{-6}$, read from each mode's own summary or log output, not from the command line. Pictures (GUI + paper render) are taken from both modes. **No cross-mode comparison is made before this gate passes.** If it fails, the difference is reported and the geometry is reconciled before launch. Levels 3 and 4 ran energy-transfer mode with $t = 1.0$, and Paper 1 ran speed-of-sound mode with $t = 0.05$. The convention has never been checked across the two modes.

**Tables printed by `python3 hspist3/validation/paper1_confinement_prereg_20261012.py --c1`** (verbatim):

#### Table L -- the length convention the c_s/L form would depend on

| anchor | eta | L_0 (geometric) | L_eff = L_0 - 2r - t/2 | (L_0/L_eff)^2 |
|---|---|---|---|---|
| 0.10 | 0.100051 | 39.2500 | 38.2250 | 1.0543 |
| 0.39 | 0.392699 | 10.0000 | 8.9750 | 1.2415 |

#### Table H -- heavy-divider form vs the exact standing wave, per alpha (box-independent)

Exact: omega^2 = c_s^2 K^2/L^2 with cot K = alpha K. Heavy form: omega^2 = 2 k_S / M_hat, k_S = N_s m c_s^2/L^2,
M_hat = M + 2 N_s m/3, i.e. omega^2 = (c_s^2/L^2)/(alpha + 1/3). The plan's M + N_s m/3 is shown for comparison.

| alpha | K | heavy/exact omega^2, M_hat = M + 2N_s m/3 | same with M + N_s m/3 | used in primary |
|---|---|---|---|---|
| 0.5 | 1.07687 | 1.03479 | 1.29349 | check only |
| 1 | 0.86033 | 1.01328 | 1.15803 | check only |
| 2 | 0.65327 | 1.00424 | 1.08149 | check only |
| 3 | 0.54716 | 1.00205 | 1.05479 | check only |
| 5 | 0.43284 | 1.00079 | 1.03308 | yes |
| 7.5 | 0.35723 | 1.00037 | 1.02211 | yes |
| 10 | 0.31105 | 1.00021 | 1.01661 | yes |
| 15 | 0.25536 | 1.00010 | 1.01109 | yes |
| 20 | 0.22176 | 1.00005 | 1.00832 | yes |

#### Table I-omega -- expected sigma of the identity residual in the length-free form

Primary (alpha >= 5): k_S^dyn = M_hat omega_1^2 / 2 per mass, inverse-variance mean over the five heavy masses;
sigma(k_S^dyn)/k_S = 2 sigma_nu/nu (M_hat exact). Per-mass sigma_nu/nu = seed SE of the A1v2 cell (canonical
estimator) at the anchor (eta = 0.1122 stands in for 0.10). Static side as in Table I (k_T noise 0.9 %).
Standing-wave check (all alpha): k_S^SW = N_s m omega_1^2 / K(alpha)^2, per mass.

| anchor | per-mass 2 sigma_nu/nu, alpha = 0.5 ... 20 [%] | heavy combined 2 sigma_nu/nu [%] | k_T/k_S | sigma(rho_I) omega-form [%] | sigma(rho_I) c_s/L form (Table I) [%] | rho_I under C [%] | rho_I(C)/sigma |
|---|---|---|---|---|---|---|---|
| 0.10 | 0.81 / 0.69 / 0.42 / 0.36 / 0.20 / 0.24 / 0.20 / 0.24 / 0.18 | 0.093 | 0.4977 | 0.457 | 0.736 | +1.269 | +2.8 |
| 0.39 | 1.01 / 0.80 / 0.68 / 0.46 / 0.39 / 0.39 / 0.33 / 0.31 / 0.24 | 0.142 | 0.4571 | 0.435 | 0.717 | +0.768 | +1.8 |


#### C3 — seeds at π/8 from the upper 1σ bound of the pilot ε₀ (2026-10-02, machine date; before any π/8 held-divider array)

**Amendment.** The seeds per position of conf_A_0.39 are set from the upper 1σ bound of the pilot's noise coefficient, $\epsilon_0 = 0.0820 + 0.0106 = 0.0926$ (§ 1.12, U1), by the **unchanged** § 1.4 rule: $T_{\rm pos} = \big((\sqrt{130}/12)\,\epsilon\,F/K\,/\,(\delta L \cdot 0.009)\big)^2/2$, $\epsilon = \epsilon_0 (100/N_s)^{1/2}$, seeds per position $= \lceil T_{\rm pos}/5000\rceil$.

**Reasons.**
1. The pilot's ε₀ carries a 13 % error, from 30 degrees of freedom (4 seeds × 5 positions × 2 faces).
2. The seed count scales as $\epsilon_0^2$ [DERIVATION, the rule above]. If the coefficient is 1σ low, every π/8 cell needs $(0.0926/0.0820)^2 = 1.28$ times the record the gate-4 seeds give it, so it gets only 78 % of that record, and the pre-registered noise budget (≤ 0.9 % on $k_T$) would not be met.
3. Under-seeding weakens the pre-registered discrimination between the hypotheses (Tables P, I and I-ω).
4. Extra seeds cannot bias the estimate. They are further independent records of the same cell, with the same stencil, record length and analysis.
5. The cost is about +2.5 core-h at KOA speed; the exact figure is printed in § 1.12 (V2).

**§ 1.8 is respected [DATA].** § 1.8 lets only the gate-4 pilot change the record of (A) at π/8, and nothing is tuned after data. No π/8 held-divider array has run. On KOA the only π/8 method-A jobs are the pilot itself (job 14966594, 20 trajectories) and its duplicate submission (job 14966614, which ran nothing; § 1.12 U2). No campaign data exist at π/8, so this choice cannot be informed by results. It changes neither the rule, nor the record per seed (5000 σ-time), nor the stencil, nor the analysis. conf_A_0.10 is unchanged.

---

### 1.10 KOA smoke test, sbatch generation, local gates (2026-10-13; nothing submitted to KOA, nothing launched)

**KOA facts used (SOURCE: the saved runbook pages, `0000_PLAN_OVERALL/ALL_MARKDOWNS/00000_KOA/`).**
- **Partitions.** `sandbox` is for tests (short runs). `shared` allocates by core, with a maximum job time of 3 days. `shared-long` allows 7 days.
- **Storage.** **`koa_scratch` has no per-user quota (800 TiB shared), but files are deleted automatically 90 days after they were last written.** Home is 50 GiB.
- The nodes are a mix of Intel and AMD CPUs from 2014 to now, so the `koa` target stays at `-march=x86-64-v2`.

**What the 90-day purge means for where the data lives (decision for Chris).** Summaries plus the two full pilot cells come back to the Mac, as recommended. **The full trajectories cannot stay on scratch "until there is an external drive"**: unless they are touched or copied, they are deleted 90 days after the run. Either they are copied to permanent storage (KoaStore / lab storage, or the external drive) within 90 days, or losing them is accepted. The summaries are enough for every pre-registered analysis, given the reduction gate below.

#### 3a — smoke test and Mac target

The file is `hspist3/cluster/koa_smoketest.sh`, submitted with sbatch. Its steps:
1. `make -B koa`
2. `--version`, which must show `-ffp-contract=off`
3. the determinism self-test: the same seed twice, and `cmp` must report the files identical
4. the $\pi/8$ pilot: the anchor cell, method B, the nine A1v2 masses with **one** seed each, 200 oscillations, seeds `run_seed(20261013, 0, m, 0)`

One script, `hspist3/cluster/confinement_pilot.py`, runs and analyses on both machines, so both go through the same code.

**The Mac target**, run 2026-10-13 on the release binary (`05215ea`, `-O3 -march=native -ffp-contract=off`). The determinism self-test on the Mac gave **IDENTICAL**. Printed by `python3 cluster/confinement_pilot.py analyse --out <mac_pi8_H10_L10>`:

pilot cell: eta (trace) = [0.392699], L_0 (trace) = [10.0], H = 10.0, N_s = 50, r = 0.5, t = 0.05 (set by --wall-thickness; not written by speed-of-sound mode), L_eff = L_0 - 2r - t/2 = 8.975000
T_total (sum of planned durations, 9 trajectories) = 69944.5 sigma-time; health lines = 0
| M | nu | implied c_s | seeds |
|---|---|---|---|
| 50 | 0.07519225 | 3.93752 | 1 |
| 100 | 0.05870708 | 3.84803 | 1 |
| 200 | 0.04395559 | 3.79433 | 1 |
| 300 | 0.03664223 | 3.77643 | 1 |
| 500 | 0.02912381 | 3.79432 | 1 |
| 750 | 0.02392290 | 3.77643 | 1 |
| 1000 | 0.02092928 | 3.79432 | 1 |
| 1500 | 0.01742544 | 3.84802 | 1 |
| 2000 | 0.01506197 | 3.83012 | 1 |

**c_s = 3.85886 +- 0.05150** (through-origin slope; +- = 1-sigma mass scatter of implied c_s)

**Gates in the script header** (fixed now):
- **determinism:** KOA run twice gives IDENTICAL.
- **geometry:** $\eta$, $L_0$, $H$ and $L_{\rm eff}$ equal to the Mac within $10^{-6}$. $t$ is not written by speed-of-sound mode, so it is checked by the mode-equivalence gate below.
- **statistics:** $|c_s^{\rm KOA} - c_s^{\rm Mac}| \le 0.05150$, the Mac pilot's 1σ mass scatter. This is the 2026-09-16 mirror-gate rule. With one seed per mass, a per-mass seed error does not exist.

**Byte-identity between the Mac and KOA is not required.** The Mac is arm64 (clang, Apple libm) and KOA is x86-64-v2 (gcc, glibc libm). Even with `-ffp-contract=off` on both, transcendental functions (log, cos and exp in the velocity draw) are not correctly rounded, and they differ in the last bit between the two libraries. Chaotic dynamics amplify one ulp within a few hundred collisions, so the two runs are independent realisations and the gate is statistical.

#### 3b — sbatch files, generated from the pre-registered cell list

`hspist3/cluster/gen_confinement_sbatch.py` writes `hspist3/cluster/confinement_20261013/`:
- one array task per cell: `conf_B_0.10`, `conf_B_0.39`, `conf_A_0.10` and `conf_A_0.39`;
- the method-A pilot `conf_A_pilot`;
- the per-cell task lists, the per-trajectory worker, and the two reductions.

**Placeholders** for Chris: `__PARTITION__`, `__ACCOUNT__`, `__SCRATCH__`, plus `__UHID__` in the fetch script. **`conf_A_0.39` is marked "submit only after the method-A pilot"** (gate 4): its seeds per position are the planning values from Table A.

**Data layout.** Under `$HD_DATA = __SCRATCH__/harddisks/hspist3` the paths are relative to `hspist3/`, exactly as on the Mac:
- method B follows the A1v2 run0 → `_run<r>.csv` harness;
- method A writes `x_<position>/ev|tr|summary|red_<seed>`.

**Reduction gate.** Run on the full Mac pilot traces, `reduce_B.py` reproduces `cell()`'s per-mass ν to $8\times10^{-17}$ (one ulp, a CSV round-trip). The same check is to be repeated on the full KOA pilot cell when it comes back.

**Printed by `python3 cluster/gen_confinement_sbatch.py`** (the summary table and the copy-back commands):

| method | cell id | N_s | H | L_0 | trajectories (seeds) | est. core-h | output dir (relative to hspist3/) |
|---|---|---|---|---|---|---|---|
| B | e0p10_H_H5_L39.25 | 25 | 5 | 39.25 | 225 (9 masses x 25) | 0.71 | `experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/confinement_B_20261013/e0p10_H_H5_L39.25` |
| A | e0p10_H_H5_L39.25 | 25 | 5 | 39.25 | 710 (5 positions x 142) | 0.16 | `experiments_energy_transfer/paper1_confinement_A_20261013/e0p10_H_H5_L39.25` |
| B | e0p10_H_H10_L39.25 | 50 | 10 | 39.25 | 225 (9 masses x 25) | 1.41 | `experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/confinement_B_20261013/e0p10_H_H10_L39.25` |
| A | e0p10_H_H10_L39.25 | 50 | 10 | 39.25 | 720 (5 positions x 144) | 0.33 | `experiments_energy_transfer/paper1_confinement_A_20261013/e0p10_H_H10_L39.25` |
| B | e0p10_H_H20_L39.25 | 100 | 20 | 39.25 | 225 (9 masses x 25) | 2.83 | `experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/confinement_B_20261013/e0p10_H_H20_L39.25` |
| A | e0p10_H_H20_L39.25 | 100 | 20 | 39.25 | 675 (5 positions x 135) | 0.62 | `experiments_energy_transfer/paper1_confinement_A_20261013/e0p10_H_H20_L39.25` |
| B | e0p10_H_H40_L39.25 | 200 | 40 | 39.25 | 225 (9 masses x 25) | 5.65 | `experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/confinement_B_20261013/e0p10_H_H40_L39.25` |
| A | e0p10_H_H40_L39.25 | 200 | 40 | 39.25 | 720 (5 positions x 144) | 1.33 | `experiments_energy_transfer/paper1_confinement_A_20261013/e0p10_H_H40_L39.25` |
| B | e0p10_L_H10_L19.625 | 25 | 10 | 19.625 | 225 (9 masses x 25) | 0.34 | `experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/confinement_B_20261013/e0p10_L_H10_L19.625` |
| A | e0p10_L_H10_L19.625 | 25 | 10 | 19.625 | 750 (5 positions x 150) | 0.17 | `experiments_energy_transfer/paper1_confinement_A_20261013/e0p10_L_H10_L19.625` |
| B | e0p10_L_H10_L78.5 | 100 | 10 | 78.5 | 225 (9 masses x 25) | 5.73 | `experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/confinement_B_20261013/e0p10_L_H10_L78.5` |
| A | e0p10_L_H10_L78.5 | 100 | 10 | 78.5 | 675 (5 positions x 135) | 0.62 | `experiments_energy_transfer/paper1_confinement_A_20261013/e0p10_L_H10_L78.5` |
| B | e0p10_aspect_H19.7917_L19.7917 | 50 | 19.7917 | 19.7917 | 225 (9 masses x 25) | 0.69 | `experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/confinement_B_20261013/e0p10_aspect_H19.7917_L19.7917` |
| A | e0p10_aspect_H19.7917_L19.7917 | 50 | 19.7917 | 19.7917 | 730 (5 positions x 146) | 0.34 | `experiments_energy_transfer/paper1_confinement_A_20261013/e0p10_aspect_H19.7917_L19.7917` |
| B | e0p10_aspect_H14_L28 | 50 | 14 | 28 | 225 (9 masses x 25) | 1.00 | `experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/confinement_B_20261013/e0p10_aspect_H14_L28` |
| A | e0p10_aspect_H14_L28 | 50 | 14 | 28 | 685 (5 positions x 137) | 0.32 | `experiments_energy_transfer/paper1_confinement_A_20261013/e0p10_aspect_H14_L28` |
| B | e0p10_aspect_H9.91667_L39.625 | 50 | 9.91667 | 39.625 | 225 (9 masses x 25) | 1.43 | `experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/confinement_B_20261013/e0p10_aspect_H9.91667_L39.625` |
| A | e0p10_aspect_H9.91667_L39.625 | 50 | 9.91667 | 39.625 | 680 (5 positions x 136) | 0.31 | `experiments_energy_transfer/paper1_confinement_A_20261013/e0p10_aspect_H9.91667_L39.625` |
| B | e0p10_aspect_H7_L56.0417 | 50 | 7 | 56.0417 | 225 (9 masses x 25) | 2.03 | `experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/confinement_B_20261013/e0p10_aspect_H7_L56.0417` |
| A | e0p10_aspect_H7_L56.0417 | 50 | 7 | 56.0417 | 685 (5 positions x 137) | 0.32 | `experiments_energy_transfer/paper1_confinement_A_20261013/e0p10_aspect_H7_L56.0417` |
| B | epi8_H_H5_L10 | 25 | 5 | 10 | 225 (9 masses x 25) | 0.48 | `experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/confinement_B_20261013/epi8_H_H5_L10` |
| A | epi8_H_H5_L10 | 25 | 5 | 10 | 1095 (5 positions x 219) | 1.57 | `experiments_energy_transfer/paper1_confinement_A_20261013/epi8_H_H5_L10` |
| B | epi8_H_H10_L10 | 50 | 10 | 10 | 225 (9 masses x 25) | 0.96 | `experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/confinement_B_20261013/epi8_H_H10_L10` |
| A | epi8_H_H10_L10 | 50 | 10 | 10 | 970 (5 positions x 194) | 2.78 | `experiments_energy_transfer/paper1_confinement_A_20261013/epi8_H_H10_L10` |
| B | epi8_H_H20_L10 | 100 | 20 | 10 | 225 (9 masses x 25) | 1.91 | `experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/confinement_B_20261013/epi8_H_H20_L10` |
| A | epi8_H_H20_L10 | 100 | 20 | 10 | 1095 (5 positions x 219) | 6.27 | `experiments_energy_transfer/paper1_confinement_A_20261013/epi8_H_H20_L10` |
| B | epi8_H_H40_L10 | 200 | 40 | 10 | 225 (9 masses x 25) | 3.83 | `experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/confinement_B_20261013/epi8_H_H40_L10` |
| A | epi8_H_H40_L10 | 200 | 40 | 10 | 2185 (5 positions x 437) | 25.01 | `experiments_energy_transfer/paper1_confinement_A_20261013/epi8_H_H40_L10` |
| B | epi8_L_H10_L5 | 25 | 10 | 5 | 225 (9 masses x 25) | 0.21 | `experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/confinement_B_20261013/epi8_L_H10_L5` |
| A | epi8_L_H10_L5 | 25 | 10 | 5 | 1095 (5 positions x 219) | 1.57 | `experiments_energy_transfer/paper1_confinement_A_20261013/epi8_L_H10_L5` |
| B | epi8_L_H10_L20 | 100 | 10 | 20 | 225 (9 masses x 25) | 4.04 | `experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/confinement_B_20261013/epi8_L_H10_L20` |
| A | epi8_L_H10_L20 | 100 | 10 | 20 | 1095 (5 positions x 219) | 6.27 | `experiments_energy_transfer/paper1_confinement_A_20261013/epi8_L_H10_L20` |
| B | epi8_aspect_H7.08333_L14.125 | 50 | 7.08333 | 14.125 | 225 (9 masses x 25) | 1.40 | `experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/confinement_B_20261013/epi8_aspect_H7.08333_L14.125` |
| A | epi8_aspect_H7.08333_L14.125 | 50 | 7.08333 | 14.125 | 1090 (5 positions x 218) | 3.12 | `experiments_energy_transfer/paper1_confinement_A_20261013/epi8_aspect_H7.08333_L14.125` |
| B | epi8_aspect_H5_L20 | 50 | 5 | 20 | 225 (9 masses x 25) | 2.02 | `experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/confinement_B_20261013/epi8_aspect_H5_L20` |
| A | epi8_aspect_H5_L20 | 50 | 5 | 20 | 970 (5 positions x 194) | 2.78 | `experiments_energy_transfer/paper1_confinement_A_20261013/epi8_aspect_H5_L20` |
| B | epi8_aspect_H3.54167_L28.2917 | 50 | 3.54167 | 28.2917 | 225 (9 masses x 25) | 2.90 | `experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/confinement_B_20261013/epi8_aspect_H3.54167_L28.2917` |
| A | epi8_aspect_H3.54167_L28.2917 | 50 | 3.54167 | 28.2917 | 870 (5 positions x 174) | 2.48 | `experiments_energy_transfer/paper1_confinement_A_20261013/epi8_aspect_H3.54167_L28.2917` |

Totals: method A 56.4 core-h, method B 39.6 core-h, all 95.9 core-h (the 261012 cost table counts the pi/8 anchor once; so does this list). Plus the pi/8 method-A pilot (20 trajectories, 0.06 core-h).

rsync back (written to cluster/confinement_20261013/fetch_confinement.sh):

```sh
#!/usr/bin/env bash
# ##CHRIS 2026-10-13: copy the confinement campaign back FROM KOA (run on the Mac, from the repo root).
# Summaries only, plus the full pilot cells; full trajectories stay on KOA scratch (deleted after 90 days).
# Fill KOA_USER and SCRATCH. Nothing on either side is deleted.
KOA_USER=__UHID__; SCRATCH="__SCRATCH__"; DTN=$KOA_USER@koa-dtn.its.hawaii.edu; R=$SCRATCH/harddisks/hspist3
SUM=(--prune-empty-dirs --include='*/' --include='red_*.csv' --include='red_nu.csv' --include='acf_runs.npz'
     --include='run.log' --include='run_*.log' --include='summary_*.csv' --include='command*.txt' --exclude='*')
rsync -av "${SUM[@]}" "$DTN:$R/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/confinement_B_20261013/" "hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/confinement_B_20261013/"
rsync -av "${SUM[@]}" "$DTN:$R/experiments_energy_transfer/paper1_confinement_A_20261013/" "hspist3/experiments_energy_transfer/paper1_confinement_A_20261013/"
# full pilot cells (every file):
rsync -av "$DTN:$R/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/confinement_pilot_20261013/koa_pi8_H10_L10/" \
          "hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/confinement_pilot_20261013/koa_pi8_H10_L10/"
rsync -av "$DTN:$R/experiments_energy_transfer/paper1_confinement_A_20261013/pilot_epi8_H_H10_L10/" "hspist3/experiments_energy_transfer/paper1_confinement_A_20261013/pilot_epi8_H_H10_L10/"
```

#### 3c — local gates on the Mac

**Mode-equivalence gate (C2): PASS.** The test cell is the $\pi/8$ anchor in both modes:
- speed-of-sound: the Mac pilot, $M = 50$;
- energy-transfer: a held divider, 200 σ, seed 9700, run through the campaign worker. That run had health 0 and $F_L = 15.97$, $F_R = 15.88$, $T = 1.000$.

**Disclosure.** The first version of the comparison reported a spurious FAIL of $2.5\times10^{-4}$ on the divider position. It had compared the speed-of-sound trace's first row, which comes one step *after* release, with energy-transfer's held position. The corrected script reads the held position from the run.log line `Initial wall_x` (`00ALLINONE.c:15789`, printed in px to 3 decimals, so 2.1e-5 σ). That tolerance replaces $10^{-6}$ wherever the speed-of-sound print is coarser, and it is marked in the table. The code lines are quoted in the script header.

Printed by `python3 hspist3/validation/paper1_modegate_20261013.py`:

#### Mode-equivalence gate, pi/8 anchor (H = L_0 = 10, N_s = 50, t = 0.05)

| quantity | speed-of-sound | energy-transfer | abs. difference | tolerance | verdict | note |
|---|---|---|---|---|---|---|
| eta (nominal) | 0.392699 | 0.392699 | 8.2e-08 | 1e-06 | PASS | both written, %.6f |
| L_0 | 10.000000 | 10.000000 | 0.0e+00 | 1e-06 | PASS | both written |
| N | 100.000000 | 100.000000 | 0.0e+00 | 0e+00 | PASS | SoS: Left+Right counts; ET: particles_total |
| H | 10.000002 | 10.000000 | 2.1e-06 | 2e-05 | PASS | SoS does not write H; inferred from its eta (6-decimal print -> ~1e-5) |
| divider centre from left wall (held) | 10.000000 | 10.000000 | 0.0e+00 | 2e-05 | PASS | SoS: run.log Initial wall_x (px, 3 dec. -> 2.1e-5); ET: W0_x_sigma |
| t | 0.050000 | 0.050000 | 2.0e-09 | 1e-06 | PASS | SoS does NOT write t (input shown); ET summary |
| left free length | 9.975000 | 9.975000 | 1.0e-09 | 2e-05 | PASS | x - t/2 (inherits the 2.1e-5 of x) |
| right free length | 9.975000 | 9.975000 | 1.0e-09 | 2e-05 | PASS | 2 L_0 - x - t/2 |
| L_eff = free length - 2r | 8.975000 | 8.975000 | 1.0e-09 | 2e-05 | PASS | SoS r not written (input) |

SegEtas cross-check (ET honours t): free length from SegEtas = 9.97501 / 9.97501 vs x - t/2 = 9.97500 (SegEtas printed to 6 decimals -> ~1e-5 sigma); with t ignored it would be 10.00000.

Grid exactness of every campaign geometry (the (int) cast at line 323 truncates 2 L_0 x 24 px):
  20 cells; 2 L_0 x 24 and H x 24 integer in all: YES

OPEN, outside this campaign: the same (int) cast on the canonical A1v2 cells with non-grid L_0 (260919 table):

| eta | L_0 (table) | 2 L_0 x 24 px | box after (int) | box shortened by [sigma] | relative |
|---|---|---|---|---|---|
| 0.019635 | 199.9995 | 9599.9760 | 9599 | 0.0407 | 1.0e-04 |
| 0.026180 | 149.9996 | 7199.9808 | 7199 | 0.0409 | 1.4e-04 |
| 0.039270 | 99.9998 | 4799.9904 | 4799 | 0.0413 | 2.1e-04 |
| 0.052360 | 74.9998 | 3599.9904 | 3599 | 0.0413 | 2.8e-04 |
| 0.078540 | 49.9999 | 2399.9952 | 2399 | 0.0415 | 4.1e-04 |
| 0.112200 | 34.9999 | 1679.9952 | 1679 | 0.0415 | 5.9e-04 |
| 0.130900 | 29.9999 | 1439.9952 | 1439 | 0.0415 | 6.9e-04 |
| 0.157080 | 24.9999 | 1199.9952 | 1199 | 0.0415 | 8.3e-04 |
| 0.549999 | 7.14 | 342.7200 | 342 | 0.0300 | 2.1e-03 |
| 0.569996 | 6.8895 | 330.6960 | 330 | 0.0290 | 2.1e-03 |
| 0.590001 | 6.6559 | 319.4832 | 319 | 0.0201 | 1.5e-03 |
| 0.609999 | 6.4377 | 309.0096 | 309 | 0.0004 | 3.1e-05 |
| 0.630002 | 6.2333 | 299.1984 | 299 | 0.0083 | 6.6e-04 |
| 0.650003 | 6.0415 | 289.9920 | 289 | 0.0413 | 3.4e-03 |
| 0.669998 | 5.8612 | 281.3376 | 281 | 0.0141 | 1.2e-03 |
| 0.679998 | 5.775 | 277.2000 | 277 | 0.0083 | 7.2e-04 |
| 0.689999 | 5.6913 | 273.1824 | 273 | 0.0076 | 6.7e-04 |
| 0.695006 | 5.6503 | 271.2144 | 271 | 0.0089 | 7.9e-04 |
| 0.699998 | 5.61 | 269.2800 | 269 | 0.0117 | 1.0e-03 |
| 0.705000 | 5.5702 | 267.3696 | 267 | 0.0154 | 1.4e-03 |
| 0.709997 | 5.531 | 265.4880 | 265 | 0.0203 | 1.8e-03 |
| 0.714999 | 5.4923 | 263.6304 | 263 | 0.0263 | 2.4e-03 |
| 0.719994 | 5.4542 | 261.8016 | 261 | 0.0334 | 3.1e-03 |
| 0.725005 | 5.4165 | 259.9920 | 259 | 0.0413 | 3.8e-03 |
| 0.730005 | 5.3794 | 258.2112 | 258 | 0.0088 | 8.2e-04 |
| 0.740006 | 5.3067 | 254.7216 | 254 | 0.0301 | 2.8e-03 |
| 0.749998 | 5.236 | 251.3280 | 251 | 0.0137 | 1.3e-03 |
| 0.759999 | 5.1671 | 248.0208 | 248 | 0.0009 | 8.4e-05 |

**GATE: PASS** -- every quantity both modes write agrees to its tolerance. t and r are not written by speed-of-sound mode; they are equal by construction (one global, set in parse_cli_options before either experiment runs), and H is pinned by the equal eta at fixed N, r, L_0.

**OPEN: a systematic in the existing Paper 1 data, found by this gate.**
- `initialize_simulation_dimensions()` sets `SIM_WIDTH = (int)(2 * L0_UNITS * PIXELS_PER_SIGMA)` (`00ALLINONE.c:323`), and the physics box is `prm.boxW = (double)(XW2 - XW1)` (15882 speed-of-sound, 16592 energy-transfer).
- So **any $L_0$ that is not a multiple of 1/48 σ runs in a box shorter than recorded**, while the recorded $\eta$ and the analysis $L_{\rm eff}$ use the untruncated $L_0$.
- In the canonical A1v2 table the shortening is up to 0.042 σ: relative $10^{-4}$ at dilute $\eta$, up to $3.8\times10^{-3}$ at $\eta \approx 0.65$–0.73 (table above).
- The effect on $c_s$ is of the same order, through both $\eta$ and $L_{\rm eff}$. It is not quantified here.
- The $\pi/8$ canonical cell ($L_0 = 10$) and every confinement cell are exact.

**Pictures gate: done for one H-scan cell and one L-scan cell**, $\pi/8$ with $H = 20, L_0 = 10$ and with $H = 10, L_0 = 20$, $N_s = 100$, held divider. They were taken with `watch.sh conf shot H L0 Ns`, a new entry, and copied to `paper1_speedofsound/experiments/final/261013_conf_H20_L10_{paper,experiment}.png` and `261013_conf_H10_L20_{...}.png`. Both show $W_{\rm in} = 0.000$ and $KE_L = KE_R = 100$.

**Found on the way: the GUI capture path starts a step piston on its own.** A capture run auto-starts the piston 200 steps after release, and the shot fires only while it moves (`00ALLINONE.c:20404–20405`). The first `conf` picture therefore showed the right gas being compressed at $u = 1.0$ ($W_{\rm in} = 400$). The `conf` entry now gives the piston $u = 0.01$ and travel 0.25, so it crosses only the 0.25 σ gap and does zero work.
- Headless runs never start a piston without `--auto-piston-step` (the gate run: $T_L = T_R = 1.000$).
- The existing `watch.sh equil` entry has no piston flags and is presumably affected the same way. Its pictures should be re-checked (OPEN).
- `watch.sh` deletes an earlier picture of the same name before shooting. The first, piston-contaminated `conf_H10_L20` render was overwritten that way; it was viewed before it was replaced.

**Not done: pictures from speed-of-sound mode.** That mode has no piston, so the automatic shot (gated on a moving piston) never fires. Capturing it needs a GUI code change, which needs a go. Geometry equality across the modes rests on the numerical gate above.

**Still unlaunched.** Nothing has been submitted. The order stays: the smoke test on KOA (`sandbox`), then the method-A pilot, then the H-scan, the L-scan and aspect, after the go.

#### 1.10.1 Smoke-test gate width (2026-10-14; § 1.10 text left as written)

**What 0.05150 is.** It is the sample **standard deviation** of the nine per-mass implied $c_s$, not a standard error. It is computed at `cluster/confinement_pilot.py:64`:

    print(f"\n**c_s = {cs:.5f} +- {np.std(imp, ddof=1):.5f}** (through-origin slope; +- = 1-sigma mass scatter of implied c_s)")

**The right width (DERIVATION).** The pilot $c_s$ is the through-origin slope, a weighted mean of the implied $c_i$ with weights $w = x^2$:
$$\mathrm{SE} = s\,\frac{\sqrt{\sum w^2}}{\sum w},\qquad \sigma_{\rm diff} = \sqrt2\,\mathrm{SE}\ \text{(two independent pilots)},\qquad \text{gate} = 2\sigma_{\rm diff}.$$

Printed by `python3 cluster/smoketest_gate_width_20261014.py`:

##### KOA smoke-test gate width, from the Mac pi/8 pilot

- pilot c_s (through-origin slope)                  = 3.85886   (reproduces the header's 3.85886)
- s = np.std(imp, ddof=1), confinement_pilot.py:64  = 0.05150   -> the 0.05150 is the SCATTER (SD) across the 9 masses
- slope weights w = x^2 (share per mass, M = 50 ... 2000): 0.368, 0.235, 0.135, 0.095, 0.059, 0.040, 0.031, 0.021, 0.016;  effective n = 4.45
- SE of the pilot c_s (slope)  = s sqrt(sum w^2)/sum w = 0.02441   (an unweighted mean would have s/3 = 0.01717)
- sigma of (KOA - Mac), two independent pilots     = sqrt(2) SE = 0.03452
- **new gate: |c_s(KOA) - c_s(Mac)| <= 2 sigma_diff = 0.06903**
- expected false-fail probability under the null (Gaussian, same scatter on KOA): 2(1 - Phi(2)) = 0.0455
- the old gate 0.05150 sat at 1.49 sigma_diff; its false-fail probability was 0.1357
- caveat (stated, not corrected): s is estimated from 9 single-seed values, so sigma_diff itself is uncertain by about 1/sqrt(2*8) = 0.25 (relative); the per-mass frequencies are quantised by the 200-period spectral bin.

**Consequence.** The slope weights concentrate on the light masses ($M = 50$ alone carries 37 %), so the effective number of masses is 4.45, not 9.
- The SE of the pilot is therefore 0.0244, not $s/3 = 0.0172$, and $\sigma_{\rm diff} = 0.0345$.
- The § 1.10 gate of 0.05150 sat at only $1.49\,\sigma_{\rm diff}$. It would have failed a correct KOA build 13.6 % of the time.
- **The gate in `koa_smoketest.sh` is now $|c_s^{\rm KOA} - c_s^{\rm Mac}| \le 0.06903$** ($2\sigma_{\rm diff}$; false-fail 4.55 % under the null). The header carries a dated amendment block, and the old lines are left in place.

**Caveat, stated.** $s$ comes from nine single-seed values, so $\sigma_{\rm diff}$ is itself uncertain by about 25 %.

### 1.11 KOA build and gates (2026-10-03; written 2026-10-02 HST on the Mac)

The date in the heading is KOA's: `date +%y%m%d` on KOA named the environment list `hd_explicit_261003.txt` while the Mac clock read 2026-10-02 HST, so KOA's shell evidently runs on UTC [INFERENCE]. This section was written on the Mac. Every value below is quoted from the KOA log Chris pasted, unless it is marked otherwise.

**Binary [DATA].** It was built in sandbox job 14966574 on cn-03-33-01, by `cluster/build_koa.sh`. `logs/BUILD_KOA_14966574.txt`:

    make compiler  gcc -> /opt/apps/software/compiler/GCCcore/14.3.0/bin/gcc -> gcc (GCC) 14.3.0   (CC=gcc; Makefile:15 'CC ?= cc')
    version        00ALLINONE  git 70b2069  target koa
    git_commit     70b20698ab21028c3cd301fd01fd7d34b0ab8706
    sha256         f15fb1f107dc0dc12d06ac821e9c471cd9ac43be2193b17fa41a9b8c9a4fe160
    libs           sdl2 2.32.56 SDL2_ttf 2.24.0 glew 2.3.0 (~/envs/hd)
    BUILD OK

`./00ALLINONE --version` prints `CFLAGS: -O2 -march=x86-64-v2 -mtune=generic -ffp-contract=off`.

**Compiler record [DATA].** `Makefile:15` is `CC              ?= cc`. `cluster/koa_env.sh` exports `CC=gcc` (f4c4756). The "make compiler" line above is the first word of `make -n -B koa`, resolved to its path and version. It shows that make invoked the GCC 14.3.0 of module `compiler/GCC/14.3.0`, not the system gcc 11.5.

**Build fixes on KOA.**
1. **`opengl.pc` [DATA].** The first pkg-config check stopped with `Package 'opengl', required by 'glu', not found`. Chris installed the one missing package into `~/envs/hd` inside a sandbox job, without changing anything already installed:

       conda install -y -p $HOME/envs/hd --override-channels -c conda-forge --freeze-installed libopengl-devel

   This added `libopengl-devel-1.7.0 ha4b6fd6_5` (16 KB, conda-forge). The environment list after it was written by `conda list -p $HOME/envs/hd --explicit > $HOME/envs/hd_explicit_$(date +%y%m%d)_opengl.txt`, so on KOA's date it is `~/envs/hd_explicit_261003_opengl.txt` (that exact name not yet confirmed with `ls`); the list before it is `~/envs/hd_explicit_261003.txt` (177 lines).
2. **libm [DATA].** The first link stopped with `undefined reference to symbol 'acos'` / `DSO missing from command line`. Commit 70b2069 adds `SYS_LIBS := -lm` in the non-Darwin branch of the Makefile (`Makefile:52`, `LIBS_BASE := -lGLEW $(GL_LIBS) $(SYS_LIBS)` at `:56`). On the Mac, SYS_LIBS stays empty, and the `make -n` link line is identical before and after (38 of 38 arguments).
3. **No git on the login node [DATA].** On `login-0102`, `git` gives `command not found`. Compute nodes have `/usr/bin/git` 2.52.0. The clone and every `git pull` therefore run inside a sandbox session (runsheet, acb0380).

**Gates [DATA].**

| gate | job | node(s) | result |
|---|---|---|---|
| smoke test: π/8 pilot $c_s$ vs Mac | 14966575 (sandbox) | cn-03-33-01 | $c_s$ = 3.81894 ± 0.02600 vs Mac 3.85886; difference −0.03992, gate ±0.06903 (§ 1.10.1): **PASS** |
| smoke test: η, $L_0$, $L_{\rm eff}$, health | 14966575 | cn-03-33-01 | **PASS** (all four) |
| determinism, same node | 14966575 | cn-03-33-01 | **IDENTICAL** |
| determinism, two nodes | 14966588 (sandbox) | cn-03-33-01 vs cn-03-33-02 | `wall_x_positions_L0_100_wallmassfactor_50_run0.csv: IDENTICAL (104506 bytes)`, `speed_of_sound_psi6.csv: IDENTICAL (353 bytes)` |

**Trace sizes, KOA vs Mac.** Both are the determinism trajectory: M = 50, 25 oscillations, seed 57831576, HD_KE_TRACE = 1.
- Mac [DATA]: `wc -lc` gives `846  104237 $SP/det_e823187/_determinism/det_A/m_50/wall_x_positions_L0_100_wallmassfactor_50_run0.csv`. That is 845 data rows plus the header. The file records `Planned_Steps 21944`, so the row count is fixed by the analytic predicted frequency, not by the trajectory.
- KOA: 104506 bytes. ~~Its row count is OPEN.~~ **Closed 2026-10-02 [DATA]:** `wc -l` on KOA (Chris's terminal) gives **846 lines, 104506 bytes**, the same row count as the Mac.
- Why the sizes differ: the row counts are equal [DATA], so the 269 bytes are digit and sign characters only. That is DATA for the row count; the cause below stays INFERENCE. The same seed on two platforms gives different trajectories. glibc's libm (KOA, gcc 14.3) and Apple's libm (Mac, clang) differ in the last bits of transcendental functions, and the chaotic collision sequence amplifies that difference, so the printed values (for example the number of minus signs; the Mac file has 393) and therefore their character counts differ. If the KOA row count is not 846, the difference is in the number of samples, not only in their characters. That would be a finding, and nothing here explains it.

**Storage decision (Chris and the plan author) [DATA].** There is no lab storage. Raw trajectories and event logs stay on `koa_scratch`, which is purged 90 days after the last write. They can be regenerated from the committed task files and seeds. Per-seed summaries (`red_*.csv`, `red_nu.csv`, `summary_*.csv`, run logs) and the full pilot cells come back to the Mac. Round-plan estimate (`cluster/round_plan_261002.py`, runsheet step 8) [INFERENCE]: conf_A_0.39 needs about 235 GiB of event logs and conf_A_0.10 about 35 GiB, while method B needs under 2 GiB per array. The KOA scratch quota has not been read yet (OPEN).

**Repository decision (Chris and the plan author) [DATA].** The repository stays public: KOA clones it anonymously over https, without a key. The unpublished Paper 3 ideas are listed by Task S, and no file has been changed for that.


### 1.12 Gate 4 result and launch plan (2026-10-02)

**Inputs [DATA, Chris's KOA terminal, 2026-10-03 UTC].**
- **Pilot, job 14966594_1:** COMPLETED 0:0; Elapsed 00:00:41; TotalCPU 06:11.764 (371.8 CPU-s for 20 trajectories, 16 in parallel). The log ends `cell pilot_epi8_H_H10_L10 done; failures: 0`.
- **Duplicate submission, job 14966614_1:** COMPLETED 0:0; Elapsed 00:00:01; TotalCPU 00:01.002. Its log text is identical (373 bytes).
- **`lfs quota`:** 155.7M used, quota 0k, limit 0k, so there is no per-user limit; 355 files.
- **Pilot summaries:** 60 files rsynced to the Mac (`red_*.csv`, `run_*.log`, `summary_*.csv`; 5 positions × seeds 9700–9703). They are committed with this section, so gate 4 can be reproduced from the repository.

#### U1 — Gate 4 (§ 1.4 rule, § 1.7 item 4, § 1.8)

Printed by `python3 cluster/gate4_pilot_261002.py` (verbatim). The 5000 ± 1 % window tolerance is this script's own operational reading of "gives 5000 σ-time"; it is not pre-registered. The ε₀ error, the Bartlett test, the impact-rate line and the "ε₀ + 1σ" column are for information only. The verdict uses the pre-registered rule alone.

##### Reproduction gate: the pre-registered rule with the planning eps0 against the task files at git 70b2069

| cell | seeds/position (rule) | seeds/position (tasks file) | reproduced |
|---|---|---|---|
| epi8_H_H5_L10 | 219 | 219 | yes |
| epi8_H_H10_L10 | 194 | 194 | yes |
| epi8_H_H20_L10 | 219 | 219 | yes |
| epi8_H_H40_L10 | 437 | 437 | yes |
| epi8_L_H10_L5 | 219 | 219 | yes |
| epi8_L_H10_L20 | 219 | 219 | yes |
| epi8_aspect_H7.08333_L14.125 | 218 | 218 | yes |
| epi8_aspect_H5_L20 | 194 | 194 | yes |
| epi8_aspect_H3.54167_L28.2917 | 174 | 174 | yes |

reproduction gate: PASS

##### Pilot: experiments_energy_transfer/paper1_confinement_A_20261013/pilot_epi8_H_H10_L10

seeds found: 20 of 20; missing files: none; health lines: 0
  red_*: 20 files, mtime (UTC) 2026-10-03 07:06:18 .. 2026-10-03 07:06:38
  run_*: 20 files, mtime (UTC) 2026-10-03 07:06:16 .. 2026-10-03 07:06:37
  summary_*: 20 files, mtime (UTC) 2026-10-03 07:06:16 .. 2026-10-03 07:06:37
  run_*.log: 1 distinct content(s); all files: 60, mtime (UTC) 2026-10-03 07:06:16 .. 2026-10-03 07:06:38
window per seed from the event log: 4999.9 .. 5000.0 sigma-time -> conversion check (5000 +- 1 %): PASS

| position | face | F mean | SD over 4 seeds | eps = SD/mean x sqrt(window) |
|---|---|---|---|---|
| x_m2 | L | 16.88662 | 0.02952 | 0.1236 |
| x_m2 | R | 15.07771 | 0.02109 | 0.0989 |
| x_m1 | L | 16.40074 | 0.02635 | 0.1136 |
| x_m1 | R | 15.46030 | 0.02020 | 0.0924 |
| x_0 | L | 15.92653 | 0.01873 | 0.0832 |
| x_0 | R | 15.92009 | 0.01253 | 0.0556 |
| x_p1 | L | 15.48972 | 0.01425 | 0.0650 |
| x_p1 | R | 16.38306 | 0.03963 | 0.1711 |
| x_p2 | L | 15.07554 | 0.03399 | 0.1594 |
| x_p2 | R | 16.89430 | 0.03301 | 0.1381 |

pooled eps at N_s = 50: 0.1160 -> eps0 (pi/8 pilot) = 0.0820 +- 0.0106 (30 degrees of freedom)  (planning value 0.2710, ratio 0.303)
Bartlett test, one relative variance across the 10 (position, face) groups: p = 0.748
impacts per face per sigma-time: 0.628 (eta 0.10005) vs 5.506 (pi/8), ratio 8.76; pure shot noise would scale eps by sqrt(1/ratio) = 0.338; measured eps(pi/8, N_s 50)/eps0(plan) = 0.428

##### Seeds per position for conf_A_0.39, by the pre-registered rule with the pilot's eps0

| cell | T per position (plan) | seeds/position (plan) | T per position (pilot eps0) | seeds/position (new) | seeds/position at eps0 + 1 sigma (info) | core-h plan (Mac model) | core-h new (Mac model) |
|---|---|---|---|---|---|---|---|
| epi8_H_H5_L10 | 1.09e+06 | 219 | 9.99e+04 | 20 | 26 | 1.57 | 0.14 |
| epi8_H_H10_L10 | 9.7e+05 | 194 | 8.88e+04 | 18 | 23 | 2.78 | 0.26 |
| epi8_H_H20_L10 | 1.09e+06 | 219 | 9.99e+04 | 20 | 26 | 6.27 | 0.57 |
| epi8_H_H40_L10 | 2.18e+06 | 437 | 2e+05 | 40 | 51 | 25.01 | 2.29 |
| epi8_L_H10_L5 | 1.09e+06 | 219 | 9.99e+04 | 20 | 26 | 1.57 | 0.14 |
| epi8_L_H10_L20 | 1.09e+06 | 219 | 9.99e+04 | 20 | 26 | 6.27 | 0.57 |
| epi8_aspect_H7.08333_L14.125 | 1.09e+06 | 218 | 9.98e+04 | 20 | 26 | 3.12 | 0.29 |
| epi8_aspect_H5_L20 | 9.7e+05 | 194 | 8.88e+04 | 18 | 23 | 2.78 | 0.26 |
| epi8_aspect_H3.54167_L28.2917 | 8.66e+05 | 174 | 7.93e+04 | 16 | 21 | 2.48 | 0.23 |

conf_A_0.39 core-h: plan 51.8 -> new 4.7 (Mac cost model); at the measured KOA speed x1.804: plan 93.5 -> new 8.6
recorded: cluster/confinement_20261013/gate4_pi8_result.txt

**GATE 4: PASS** -- conf_A_0.39 uses the new seeds per position (its tasks files must be regenerated before Round 2). conf_A_0.10 is NOT affected: its eps0 was measured at eta = 0.10005 (Level 3 c0), and sec. 1.8 allows the pilot to change only (A) at pi/8.

**Reading [INFERENCE].** ε₀ at π/8 is 0.30 of the planning value, so the seeds per position fall by a factor of about 11. That is a large change, so here is why the comparison is like for like:
- **Same definition.** Both numbers are the relative per-face divider force noise × √T, with F = Σ|dp| of the `D0` events over the event-log time range, at H = 10. The planning value comes from `noise_eps0()` (Level 3 c0); the pilot value from `reduce_A.py`.
- **The size of the drop is plausible.** The pre-registration ASSUMED that ε depends on η only through $N_s$. At π/8 the divider is hit 8.8 times more often per face, which on its own (shot noise) would scale ε by 0.34; the measurement gives 0.43.
- **The groups agree.** The ten (position, face) groups share one relative variance (Bartlett p = 0.75).

Over 30 degrees of freedom the ε₀ error is ±13 %. At ε₀ + 1σ, the rule would give about 26 % more seeds; that column is information only.

#### U2 — The duplicate pilot submission (job 14966614)

[DATA] It re-ran nothing:
- every one of the 60 rsynced files has a KOA mtime between 07:06:16 and 07:06:38 UTC (preserved by `rsync -a`), which lies within job 14966594's 41 s, and none is from 07:11;
- the 20 `run_*.log` files have one identical content;
- the `summary_*.csv` timestamps are 07:06.

[DERIVATION] It could not have re-run anything without leaving a trace:
- `conf_worker.sh` mode A exits before the binary starts when `red_<seed>.csv` is non-empty: `[ -s "$d/red_${seed}.csv" ] && exit 0`;
- the binary truncates `run_<seed>.log` (`> "$d/run_${seed}.log"`), so any rerun would have given that file a 07:11 mtime.

So no file was written by both jobs. The event logs and traces (`ev_*`, `tr_*`) stayed on scratch and their mtimes were not inspected, but by the same exit line they were not touched. **The pilot is usable.**

#### U3 — Task files with the gate-4 seeds; round plan at the measured KOA speed

**Seeds.**
- `cluster/gate4_pilot_261002.py` records the result in `cluster/confinement_20261013/gate4_pi8_result.txt` (ε₀ 0.0820079 ± 0.0105872, planning value 0.2709969, PASS).
- `cluster/gen_confinement_sbatch.py` reads that file and scales $T_{\rm pos}$ by $(\epsilon_0^{\rm pilot}/\epsilon_0^{\rm plan})^2$ for the π/8 cells only, which is exact because the pre-registered $T_{\rm pos} \propto \epsilon_0^2$ [DERIVATION]. It then sets seeds per position = ⌈$T_{\rm pos}$/5000⌉ and rewrites the task files.
- A run without the gate-4 file first showed that the generator reproduces every committed task, cell and summary file byte for byte. With the file, only the nine `tasks_A_epi8_*.txt` files and `cells_summary.txt` changed. All B task files, the A_0.10 task files and the pilot's task file are unchanged.

**KOA speed [DATA → DERIVATION].**
- KOA cost: 371.764 CPU-s / 20 = 18.59 CPU-s per held-divider trajectory (5000 σ-time, π/8, $N_s$ = 50 per side, including the 200 σ-time hold and the Python reduction).
- Mac cost model for the same trajectory: 10.30 CPU-s.
- So the factor is **1.804**. It was measured on method A and is applied to method B too [INFERENCE].

**`--time` rule:** `--time` ≥ 2 × the longest cell at KOA speed, rounded up to 15 min, at least 30 min (`round_plan_261002.time_limit_h`). The generator writes it into every array.

Printed by `python3 cluster/round_plan_261002.py` (verbatim):

##### conf_A_0.39 task files: planning seeds (git 70b2069) vs gate-4 seeds (working tree)

| cell | seeds/position at 70b2069 | seeds/position now | lines now = first lines of each position at 70b2069 | seeds now |
|---|---|---|---|---|
| epi8_H_H5_L10 | 219 | 20 | yes | 9700..9719 |
| epi8_H_H10_L10 | 194 | 18 | yes | 9700..9717 |
| epi8_H_H20_L10 | 219 | 20 | yes | 9700..9719 |
| epi8_H_H40_L10 | 437 | 40 | yes | 9700..9739 |
| epi8_L_H10_L5 | 219 | 20 | yes | 9700..9719 |
| epi8_L_H10_L20 | 219 | 20 | yes | 9700..9719 |
| epi8_aspect_H7.08333_L14.125 | 218 | 20 | yes | 9700..9719 |
| epi8_aspect_H5_L20 | 194 | 18 | yes | 9700..9717 |
| epi8_aspect_H3.54167_L28.2917 | 174 | 16 | yes | 9700..9715 |

KOA speed (measured, pilot 14966594): 371.764 CPU-s / 20 = 18.59 CPU-s per trajectory; Mac cost model 10.30 -> factor 1.804
times below: cost model x 1.804 (--slow); trajectories per cell from the task files

| round | array | partition | tasks | cores/task | traj. | core-h (KOA) | longest cell (h) | --time (h) | rule 2 x longest (h) | throttle | cores at once | wall (h) | scratch GiB |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Round 1 | conf_B_0.10 | shared | 10 | 8 | 2250 | 39.4 | 1.33 | 2.75 | 2.75 | %2 | 16 | 2.71 | 1.8 |
| Round 1 | conf_B_0.39 | shared | 9 | 8 | 2025 | 32.0 | 0.94 | 2 | 2 | %2 | 16 | 2.22 | 1.6 |
| Round 1 | conf_A_0.10 | shared | 10 | 16 | 7030 | 8.2 | 0.15 | 0.5 | 0.5 | %2 | 32 | 0.26 | 34.9 |
| Round 2 | conf_A_0.39 | shared | 9 | 16 | 960 | 8.6 | 0.27 | 0.75 | 0.75 | %4 | 64 | 0.27 | 21.6 |

Round 1: cores at once = 64 (limit 64) -> OK
Round 2: cores at once = 64 (limit 64) -> OK

sbatch lines (from ~/harddisks/hspist3, after `mkdir -p logs`):

    Round 1:  sbatch --array=1-10%2 cluster/confinement_20261013/conf_B_0.10.sbatch
    Round 1:  sbatch --array=1-9%2 cluster/confinement_20261013/conf_B_0.39.sbatch
    Round 1:  sbatch --array=1-10%2 cluster/confinement_20261013/conf_A_0.10.sbatch
    Round 2:  sbatch --array=1-9%4 cluster/confinement_20261013/conf_A_0.39.sbatch

**Scratch:** the A_0.39 event logs now come to about 22 GiB (235 GiB planned); with no per-user quota, the § 1.11 storage question is closed. **Cost:** the A_0.39 core-hours drop from 93.5 to 8.6 at KOA speed. Round 1 dominates (about 80 core-hours, about 2.7 h of wall time).

#### U4 — Build-hash guard and the recorded build hash

**Arrays (all five `conf_*.sbatch`, generated).**

Old:

    "$HD_BIN" --version | head -1 | grep -q -- "git $(git rev-parse --short HEAD)  target koa" || { echo "STOP: not the clean koa build of HEAD"; exit 1; }

New:

    [ -s logs/BUILD_KOA_LAST.hash ] || { echo "STOP: no logs/BUILD_KOA_LAST.hash -- build with cluster/build_koa.sh first"; exit 1; }
    sha256sum --status -c logs/BUILD_KOA_LAST.hash || { echo "STOP: ./00ALLINONE is not the build recorded in logs/BUILD_KOA_LAST.hash"; exit 1; }
    export HD_BUILD="$("$HD_BIN" --version | head -1)"
    echo "$HD_BUILD" | grep -Eq -- "git [0-9a-f]+  target koa" || { echo "STOP: not a clean koa build: $HD_BUILD"; exit 1; }
    command -v flock >/dev/null || { echo "STOP: flock not found (conf_worker.sh needs it)"; exit 1; }

**Build (`cluster/build_koa.sh`, new lines before `BUILD OK`).** The build itself still requires build_git = HEAD, with no `-dirty`:

    sha256sum 00ALLINONE > logs/BUILD_KOA_LAST.hash
    echo "recorded       logs/BUILD_KOA_LAST.hash: $(cat logs/BUILD_KOA_LAST.hash)"

**Worker (`conf_worker.sh`, both modes; the pilot path is the same worker).**
- A new function `guard <dir> <glob>` reads `<dir>/.build_git` under `flock`, or creates it if absent. It refuses when:
  - the recorded line differs from this binary's `--version` line, or
  - outputs are present but no record exists.
- Called as:

      guard "$cell" 'wall_x_positions_L0_*_run*.csv' || { echo "B $rel M=$M r=$r FAILED build guard"; exit 3; }
      guard "$d" 'red_*.csv' || { echo "A $rel seed=$seed FAILED build guard"; exit 3; }

  The first line comes before the B "done before" skip, the second before the A skip.
- Tested on the Mac with a stub `flock`, five cases: fresh, resume with the same build, other build, outputs without a record, and the B glob. All behave as specified.
- The locking itself is not tested, because macOS has no `flock`; KOA's presence is checked by the sbatch [OPEN until the first array].
- **Consequence:** the existing KOA pilot directories have outputs but no `.build_git`, so a resubmitted pilot would now be refused. That is intended.
- `cluster/koa_crossnode_det.sh` still compares against HEAD. It is a one-off test, not a campaign, and is unchanged.

**Runsheet (step 8, new rules):**
1. after every `git pull`, rebuild before any NEW submission;
2. never pull or rebuild while array tasks are pending or running;
3. the first submission after this change needs a pull and a rebuild, because the 70b2069 build wrote no hash file.

Step 8e gives the commands.

#### Launch lines (from `~/harddisks/hspist3` on `login-0102`, after the pull and rebuild of step 8e and after the go)

Round 1 (64 cores at once):

    sbatch --array=1-10%2 cluster/confinement_20261013/conf_B_0.10.sbatch
    sbatch --array=1-9%2 cluster/confinement_20261013/conf_B_0.39.sbatch
    sbatch --array=1-10%2 cluster/confinement_20261013/conf_A_0.10.sbatch

Round 2 (64 cores at once; after its go):

    sbatch --array=1-9%4 cluster/confinement_20261013/conf_A_0.39.sbatch

**Free cross-check [DERIVATION].** The anchor cell `epi8_H_H10_L10` reruns the pilot's seeds 9700–9703. The C source is the same as the pilot's, only the build hash differs, so its `red_970[0-3].csv` must equal the pilot's byte for byte.


#### V2 — conf_A_0.39 under amendment C3 (2026-10-02)

The generator (`cluster/gen_confinement_sbatch.py`) reads `eps0_c3_upper 0.09259510508078693` = ε₀ + 1σ from `cluster/confinement_20261013/gate4_pi8_result.txt`, written by `cluster/gate4_pilot_261002.py`. It rewrote the nine conf_A_0.39 task files, `cells_summary.txt`, and the header comment of `conf_A_0.39.sbatch`; no other file changed. Rounding ε₀ to 0.0926 gives the same seeds.

Printed by `python3 cluster/round_plan_261002.py` (verbatim):

##### conf_A_0.39 task files: seeds per position, plan (git 70b2069), gate 4 (git 303280d) -> now (working tree)

| cell | plan (70b2069) | gate 4 (303280d) | now | lines | nested in each earlier file | seeds now |
|---|---|---|---|---|---|---|
| epi8_H_H5_L10 | 219 | 20 | 26 | 130 | yes | 9700..9725 |
| epi8_H_H10_L10 | 194 | 18 | 23 | 115 | yes | 9700..9722 |
| epi8_H_H20_L10 | 219 | 20 | 26 | 130 | yes | 9700..9725 |
| epi8_H_H40_L10 | 437 | 40 | 51 | 255 | yes | 9700..9750 |
| epi8_L_H10_L5 | 219 | 20 | 26 | 130 | yes | 9700..9725 |
| epi8_L_H10_L20 | 219 | 20 | 26 | 130 | yes | 9700..9725 |
| epi8_aspect_H7.08333_L14.125 | 218 | 20 | 26 | 130 | yes | 9700..9725 |
| epi8_aspect_H5_L20 | 194 | 18 | 23 | 115 | yes | 9700..9722 |
| epi8_aspect_H3.54167_L28.2917 | 174 | 16 | 21 | 105 | yes | 9700..9720 |

KOA speed (measured, pilot 14966594): 371.764 CPU-s / 20 = 18.59 CPU-s per trajectory; Mac cost model 10.30 -> factor 1.804
times below: cost model x 1.804 (--slow); trajectories per cell from the task files

| round | array | partition | tasks | cores/task | traj. | core-h (KOA) | longest cell (h) | --time (h) | rule 2 x longest (h) | throttle | cores at once | wall (h) | scratch GiB |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Round 1 | conf_B_0.10 | shared | 10 | 8 | 2250 | 39.4 | 1.33 | 2.75 | 2.75 | %2 | 16 | 2.71 | 1.8 |
| Round 1 | conf_B_0.39 | shared | 9 | 8 | 2025 | 32.0 | 0.94 | 2 | 2 | %2 | 16 | 2.22 | 1.6 |
| Round 1 | conf_A_0.10 | shared | 10 | 16 | 7030 | 8.2 | 0.15 | 0.5 | 0.5 | %2 | 32 | 0.26 | 34.9 |
| Round 2 | conf_A_0.39 | shared | 9 | 16 | 1240 | 11.0 | 0.33 | 0.75 | 0.75 | %4 | 64 | 0.33 | 27.8 |

Printed by `python3 cluster/gate4_pilot_261002.py`: `amendment C3 (eps0 + 1 sigma = 0.0926): conf_A_0.39 core-h 6.1 (Mac cost model), 11.0 at KOA speed (+2.5 over the gate-4 seeds)`. The longest cell is 0.33 h, so `--time` stays at 0:45 (rule ≥ 2 × 0.33 h, rounded up to 15 min). `wc -l cluster/confinement_20261013/tasks_A_epi8_H_H10_L10.txt` gives **115** (23 seeds × 5 positions). The "nested" column shows that every file only adds seeds at the end of each position. The gate-4 and planning seed lists are prefixes of the C3 list.

#### Analysis-plan gate: determinism across jobs (added 2026-10-02, before any π/8 array)

The anchor cell `epi8_H_H10_L10` of conf_A_0.39 runs seeds 9700–9722 at the five positions x_m2 … x_p2 with the same command line as the pilot.
- **Gate:** at every position, the anchor's `red_9700.csv` … `red_9703.csv` must be byte-identical (`cmp`) to the pilot's (`pilot_epi8_H_H10_L10/x_*/red_970[0-3].csv`, committed in df53ba1). That is 20 comparisons.
- **What it tests:** determinism across jobs, nodes and builds of the same C source. The pilot ran at 70b2069; Round 2 runs on the rebuild of step 8e, where only the build hash string differs.
- **If any pair differs:** the difference is reported and **no conf_A_0.39 result is used** until it is explained.
- **When it is checked:** on the Mac, after the summaries of Round 2 are copied back, before any analysis.


### 1.13 Round 1 outcome and repair (2026-10-03)

**Round 1 [DATA, Chris's KOA terminal, 2026-10-03].** Build `279282b target koa`. The arrays were conf_B_0.10 (14967049), conf_B_0.39 (14967050) and conf_A_0.10 (14967051). conf_A_0.39 (14967120) was cancelled by its `afterok` dependency, because A_0.10 task 4 timed out. Scratch use: 3.1 GiB, 41.3 k files (`koa_scratch` limit 800 TiB).

The `done; failures: N` lines (`FAILED build guard`):
- A_0.10: H10 1, L_H10_L19.625 6, aspect_H7_L56.0417 2;
- B_0.10: aspect_H9.91667 2, aspect_H7 2;
- B_0.39: aspect_H5_L20 1, aspect_H3.54167_L28.2917 3;
- all other finished cells 0.

`flock: 9: Bad file descriptor` appears once per trajectory in every array task (166–750 times) and never in the pilot or smoke logs.

#### Cause 1 — the build guard never locked (CC's error in U4)

The U4 guard locked with:

    have=$( flock 9
            ...
            cat "$dir/.build_git" ) 9>"$dir/.build_git.lock"

On an assignment, the command substitution is expanded **before** the redirection opens fd 9. So `flock` had no file and never locked, which is the "Bad file descriptor" line. The Mac test of U4 stubbed `flock` and could not see this.

Without the lock, the 8–16 workers that start together in one directory raced: `printf … > .build_git` (truncate, then write) against `cat`. A worker that read the file in between saw a different "build" and refused (`FAILED build guard`) [DERIVATION from the code; the REFUSED detail lines were not pasted, OPEN]. A refused worker exits **before** anything of its trajectory exists. Its seed was skipped, not spoiled (W2).

#### Cause 2 — the H = 40 cells cost far more than the model

In all three arrays, task 4 is the H = 40 cell ($N_s$ = 200), and it hit `--time`. Printed by `python3 cluster/round1_timing_261003.py` (verbatim):

##### Round 1, measured wall per wave (sacct Elapsed / ceil(n/P))

| array | task | cell | N_s | trajectories | P | waves | Elapsed | state | wall per wave (s) |
|---|---|---|---|---|---|---|---|---|---|
| B_0.10 | 1 | e0p10_H_H5_L39.25 | 25 | 225 | 8 | 29 | 00:06:31 | COMPLETED | 13.5 |
| B_0.10 | 2 | e0p10_H_H10_L39.25 | 50 | 225 | 8 | 29 | 00:15:42 | COMPLETED | 32.5 |
| B_0.10 | 3 | e0p10_H_H20_L39.25 | 100 | 225 | 8 | 29 | 00:46:07 | COMPLETED | 95.4 |
| B_0.10 | 4 | e0p10_H_H40_L39.25 | 200 | 225 | 8 | 29 | 02:47:01 | TIMEOUT | >345.6 |
| B_0.10 | 5 | e0p10_L_H10_L19.625 | 25 | 225 | 8 | 29 | 00:03:12 | COMPLETED | 6.6 |
| B_0.10 | 6 | e0p10_L_H10_L78.5 | 100 | 225 | 8 | 29 | 01:22:51 | COMPLETED | 171.4 |
| B_0.10 | 7 | e0p10_aspect_H19.7917_L19.7917 | 50 | 225 | 8 | 29 | 00:09:08 | COMPLETED | 18.9 |
| B_0.10 | 8 | e0p10_aspect_H14_L28 | 50 | 225 | 8 | 29 | 00:11:44 | COMPLETED | 24.3 |
| B_0.10 | 9 | e0p10_aspect_H9.91667_L39.625 | 50 | 225 | 8 | 29 | 00:15:32 | COMPLETED | 32.1 |
| B_0.10 | 10 | e0p10_aspect_H7_L56.0417 | 50 | 225 | 8 | 29 | 00:20:27 | COMPLETED | 42.3 |
| B_0.39 | 1 | epi8_H_H5_L10 | 25 | 225 | 8 | 29 | 00:04:08 | COMPLETED | 8.6 |
| B_0.39 | 2 | epi8_H_H10_L10 | 50 | 225 | 8 | 29 | 00:08:14 | COMPLETED | 17.0 |
| B_0.39 | 3 | epi8_H_H20_L10 | 100 | 225 | 8 | 29 | 01:05:33 | COMPLETED | 135.6 |
| B_0.39 | 4 | epi8_H_H40_L10 | 200 | 225 | 8 | 29 | 02:02:00 | TIMEOUT | >252.4 |
| B_0.39 | 5 | epi8_L_H10_L5 | 25 | 225 | 8 | 29 | 00:02:30 | COMPLETED | 5.2 |
| B_0.39 | 6 | epi8_L_H10_L20 | 100 | 225 | 8 | 29 | 01:33:20 | COMPLETED | 193.1 |
| B_0.39 | 7 | epi8_aspect_H7.08333_L14.125 | 50 | 225 | 8 | 29 | 00:10:19 | COMPLETED | 21.3 |
| B_0.39 | 8 | epi8_aspect_H5_L20 | 50 | 225 | 8 | 29 | 00:13:17 | COMPLETED | 27.5 |
| B_0.39 | 9 | epi8_aspect_H3.54167_L28.2917 | 50 | 225 | 8 | 29 | 00:17:56 | COMPLETED | 37.1 |
| A_0.10 | 1 | e0p10_H_H5_L39.25 | 25 | 710 | 16 | 45 | 00:01:56 | COMPLETED | 2.6 |
| A_0.10 | 2 | e0p10_H_H10_L39.25 | 50 | 720 | 16 | 45 | 00:03:03 | COMPLETED | 4.1 |
| A_0.10 | 3 | e0p10_H_H20_L39.25 | 100 | 675 | 16 | 43 | 00:06:42 | COMPLETED | 9.3 |
| A_0.10 | 4 | e0p10_H_H40_L39.25 | 200 | 720 | 16 | 45 | 00:32:06 | TIMEOUT | >42.8 |
| A_0.10 | 5 | e0p10_L_H10_L19.625 | 25 | 750 | 16 | 47 | 00:02:04 | COMPLETED | 2.6 |
| A_0.10 | 6 | e0p10_L_H10_L78.5 | 100 | 675 | 16 | 43 | 00:06:29 | COMPLETED | 9.0 |
| A_0.10 | 7 | e0p10_aspect_H19.7917_L19.7917 | 50 | 730 | 16 | 46 | 00:03:46 | COMPLETED | 4.9 |
| A_0.10 | 8 | e0p10_aspect_H14_L28 | 50 | 685 | 16 | 43 | 00:03:19 | COMPLETED | 4.6 |
| A_0.10 | 9 | e0p10_aspect_H9.91667_L39.625 | 50 | 680 | 16 | 43 | 00:03:10 | COMPLETED | 4.4 |
| A_0.10 | 10 | e0p10_aspect_H7_L56.0417 | 50 | 685 | 16 | 43 | 00:03:07 | COMPLETED | 4.3 |

##### Fit over the H-scan cells H5, H10, H20 (N_s 25, 50, 100) -> H40 (N_s 200)

| array | exponent p (fit) | local exponent H10->H20 | H40 predicted, fit (h) | H40 predicted, local (h) | H40 TIMEOUT Elapsed (h) | fit consistent with TIMEOUT? |
|---|---|---|---|---|---|---|
| B_0.10 | 1.41 | 1.55 | 1.98 | 2.26 | > 2.78 | **NO** (fit below the timeout) |
| B_0.39 | 1.99 | 2.99 | 3.45 | 8.70 | > 2.03 | yes |
| A_0.10 | 0.93 | 1.20 | 0.21 | 0.27 | > 0.54 | **NO** (fit below the timeout) |

lower bound on the H20 -> H40 exponent from the timeouts: B_0.10 > 1.86, B_0.39 > 0.90, A_0.10 > 2.19
USED exponent p* = steepest local exponent measured = 2.99; exceeds every timeout bound: yes

##### --time for the resubmission (p* scaling, 3 x the whole cell)

| array | task(s) | cell | predicted whole cell (h) | basis | --time |
|---|---|---|---|---|---|
| B_0.10 | 4 | e0p10_H_H40_L39.25 | 6.12 | H20 95.4 s/wave x 2^2.99 x 29 waves (> timeout 2.78 h: yes) | 18:30:00 |
| B_0.10 | 1,2,3,5,6,7,8,9,10 | the others (COMPLETED; only missing seeds run) | -- | default | 02:45:00 (sbatch) |
| B_0.39 | 4 | epi8_H_H40_L10 | 8.70 | H20 135.6 s/wave x 2^2.99 x 29 waves (> timeout 2.03 h: yes) | 1-02:15:00 |
| B_0.39 | 1,2,3,5,6,7,8,9 | the others (COMPLETED; only missing seeds run) | -- | default | 02:00:00 (sbatch) |
| A_0.10 | 4 | e0p10_H_H40_L39.25 | 0.93 | H20 9.3 s/wave x 2^2.99 x 45 waves (> timeout 0.54 h: yes) | 03:00:00 |
| A_0.10 | 1,2,3,5,6,7,8,9,10 | the others (COMPLETED; only missing seeds run) | -- | default | 00:30:00 (sbatch) |
| A_0.39 | 1 | epi8_H_H5_L10 (N_s 25, 130 traj.) | 0.05 | pilot 20.5 s/wave x (N_s/50)^2.99 x 9 waves | 00:45:00 (sbatch default >= 3 x) |
| A_0.39 | 2 | epi8_H_H10_L10 (N_s 50, 115 traj.) | 0.05 | pilot 20.5 s/wave x (N_s/50)^2.99 x 8 waves | 00:45:00 (sbatch default >= 3 x) |
| A_0.39 | 3 | epi8_H_H20_L10 (N_s 100, 130 traj.) | 0.41 | pilot 20.5 s/wave x (N_s/50)^2.99 x 9 waves | 01:15:00 |
| A_0.39 | 4 | epi8_H_H40_L10 (N_s 200, 255 traj.) | 5.78 | pilot 20.5 s/wave x (N_s/50)^2.99 x 16 waves | 17:30:00 |
| A_0.39 | 5 | epi8_L_H10_L5 (N_s 25, 130 traj.) | 0.05 | pilot 20.5 s/wave x (N_s/50)^2.99 x 9 waves | 00:45:00 (sbatch default >= 3 x) |
| A_0.39 | 6 | epi8_L_H10_L20 (N_s 100, 130 traj.) | 0.41 | pilot 20.5 s/wave x (N_s/50)^2.99 x 9 waves | 01:15:00 |
| A_0.39 | 7 | epi8_aspect_H7.08333_L14.125 (N_s 50, 130 traj.) | 0.05 | pilot 20.5 s/wave x (N_s/50)^2.99 x 9 waves | 00:45:00 (sbatch default >= 3 x) |
| A_0.39 | 8 | epi8_aspect_H5_L20 (N_s 50, 115 traj.) | 0.05 | pilot 20.5 s/wave x (N_s/50)^2.99 x 8 waves | 00:45:00 (sbatch default >= 3 x) |
| A_0.39 | 9 | epi8_aspect_H3.54167_L28.2917 (N_s 50, 105 traj.) | 0.04 | pilot 20.5 s/wave x (N_s/50)^2.99 x 7 waves | 00:45:00 (sbatch default >= 3 x) |

##### Resubmission lines (from ~/harddisks/hspist3 on login-0102; the skip logic runs only what is missing)

Set 1 (now); at most 64 cores at once:
    sbatch --array=4 --time=18:30:00 cluster/confinement_20261013/conf_B_0.10.sbatch
    sbatch --array=1-3,5-10%1 cluster/confinement_20261013/conf_B_0.10.sbatch
    sbatch --array=4 --time=1-02:15:00 cluster/confinement_20261013/conf_B_0.39.sbatch
    sbatch --array=1-3,5-9%1 cluster/confinement_20261013/conf_B_0.39.sbatch
    sbatch --array=4 --time=03:00:00 cluster/confinement_20261013/conf_A_0.10.sbatch
    sbatch --array=1-3,5-10%1 cluster/confinement_20261013/conf_A_0.10.sbatch
Set 2 (Round 2, behind the A_0.10 resubmission); at most 48 cores at once:
    sbatch --array=1,2,5,7,8,9%1 --dependency=afterok:<A_0.10 task-4 jobid>:<A_0.10 rest jobid> cluster/confinement_20261013/conf_A_0.39.sbatch
    sbatch --array=3,6%1 --time=01:15:00 --dependency=afterok:<A_0.10 task-4 jobid>:<A_0.10 rest jobid> cluster/confinement_20261013/conf_A_0.39.sbatch
    sbatch --array=4 --time=17:30:00 --dependency=afterok:<A_0.10 task-4 jobid>:<A_0.10 rest jobid> cluster/confinement_20261013/conf_A_0.39.sbatch
Set 2 starts while the two long B H40 tasks of set 1 may still run: 16 + set 2 = 64 cores (cap 64).


**Reading [DATA → INFERENCE].** The H-scan fit over $N_s$ = 25, 50, 100 is falsified by two of the three timeouts: its H40 prediction is below the time those cells had already run. The cost per trajectory grows steeper than any single power law over that range:
- A_0.10 goes up ≥ 4.6× from $N_s$ 100 to 200 (exponent ≥ 2.19);
- B_0.39 goes up 8× from 50 to 100 (2.99).

The pre-registered cost model (per σ-time ∝ N) and the KOA speed factor (measured at $N_s$ = 50) therefore underestimate the large cells. The cause is not known; something in the code may scale like $N^2$ or worse at N = 400 (OPEN; no profiling, since there are no simulation runs on the Mac). **Used for every `--time`:** the steepest exponent measured, p* = 2.99, which exceeds every timeout bound, × 3, capped at 3 days. conf_A_0.39 gets the same scaling from the pilot ($N_s$ = 50), so its $N_s$ = 100 cells (tasks 3, 6) now get 1:15 instead of 0:45, and its H40 cell 17:30.

#### W1 — the lock fix (`conf_worker.sh`)

The guard now locks with `mkdir`, which is atomic on every file system:

    until mkdir "$lock" 2>/dev/null; do
      i=$((i + 1)); [ "$i" -le 600 ] || { echo "REFUSED $dir: lock $lock not free within 60 s"; return 1; }
      sleep 0.1
    done
    trap 'rmdir "$lock" 2>/dev/null' EXIT
    trap 'rmdir "$lock" 2>/dev/null; exit 143' TERM INT
    if [ ! -e "$dir/.build_git" ]; then
      if compgen -G "$dir/$2" >/dev/null; then have="(none recorded, outputs present)"
      else printf '%s\n' "$BUILD" > "$dir/.build_git.tmp$$" && mv "$dir/.build_git.tmp$$" "$dir/.build_git"; fi
    fi
    [ -n "$have" ] || have=$(cat "$dir/.build_git")
    rmdir "$lock"; trap - EXIT TERM INT

- **Lock and record.** The lock is released on SIGTERM too, which is Slurm's TIMEOUT. `.build_git` is renamed into place, so it is never seen half-written.
- **Stale locks.** A lock that stays held is not broken, because breaking it cannot be made race-free. The worker refuses and the seed stays missing; `check_cells.sh` lists the lock.
- **Partial files.** A refused worker writes nothing.
- **Killed runs.** New lines move the partial files of a killed run aside before the rerun; they are neither reused nor deleted:
  - B: `[ -e "$tmp" ] && mv "$tmp" "$cell/.stale_run${r}_$(date +%Y%m%d_%H%M%S)"`;
  - A: the seed's `ev_/tr_/summary_/run_` files go to `.stale_<seed>_<date>/`. The binary opens `summary_<seed>.csv` in append mode (`00ALLINONE.c:17287`, `fopen(summary_path, "a")`), so a rerun would otherwise leave two rows. `ev_` and `tr_` are opened with `"w"` (`edmd.c:1057`; `00ALLINONE.c:16531`, `FILE *elog = fopen(trace_path, "w")`).
- **sbatch generator.** The `command -v flock` check is removed (`gen_confinement_sbatch.py`; one line per sbatch). The B `run.log` append keeps its `flock`, which has the correct form `( flock 9; … ) 9>file` and works on KOA.

**Mac test** of the real `conf_worker.sh`, with a stub binary and a stub `reduce_A.py`; all as expected:
1. 16 parallel workers, 48 seeds, fresh A directory: 0 FAILED, 48 `red_`, `.build_git` correct, no lock left.
2. Resume of the same seeds: 0 binary calls.
3. 16 new seeds racing into a directory with outputs and a record: 0 FAILED.
4. Another build on that directory: all 16 refused, 0 binary calls, directory unchanged.
5. Outputs without a record: all 16 refused, directory unchanged.
6. A killed seed (ev, summary and run log, no red): moved to `.stale_9790_*`, rerun, the new summary has 1 row.
7. 32 B runs, 16 parallel, with a stale `.run5`: 0 FAILED, 32 traces, `.run5` moved aside.
8. A lock held by someone else: refused after 60 s, nothing written.

#### W2 — data integrity of Round 1 [DERIVATION from `conf_worker.sh`]

**No output can be written twice.** Every output name belongs to one task line, and each line (B: mass M, run r; A: position, seed) appears once in its task file and is run once by `xargs`.
- B runs in its own `.run$r` inside `m_$M` and moves the finished trace to `"$cell/$(basename "${tr%run0.csv}")run$r.csv"`.
- A writes `ev_${seed}.csv`, `tr_${seed}.csv`, `summary_${seed}.csv`, `run_${seed}.log` and `red_${seed}.csv` in its position directory.

The race was in the guard, which runs before any of these. It could skip a trajectory but never write one twice.

What a resume does:
- **(a) A seed skipped by the race:** none of its files exist. The guard passes, because the record matches the unchanged binary; the skip test (`[ -s red_<seed>.csv ]`, or the B trace) fails, so the seed runs.
- **(b) A directory whose `.build_git` holds `00ALLINONE  git 279282b  target koa`:** it matches `HD_BUILD` of the same binary, so the resume proceeds. **A rebuild at a new commit would refuse every such directory**, hence runsheet rule 5: no rebuild while any cell is incomplete.
- **(c) The three TIMEOUT cells:**
  - finished trajectories keep their outputs and are skipped;
  - a trajectory killed in flight left only partial files (B: inside `.run<r>`, since its trace is moved into the cell only on success; A: `ev_/tr_/summary_/run_` without `red_`), which are moved aside and the trajectory is rerun;
  - `reduce_B.py` never ran on the two B H40 cells, and runs at the end of their resubmission.

  On a resubmitted B cell `reduce_B.py` rewrites `red_nu.csv` from the same traces.

**`cluster/check_cells.sh`** (bash plus Python from `~/envs/hd`; read-only; refuses to run on the login node). Per cell it prints:
- trajectories expected (task file) against present (non-empty);
- the missing list as ranges;
- the `.build_git` record(s), and directories without a record;
- zero-size, duplicate, unexpected and malformed (not 2-line) outputs;
- partial seeds, leftovers (`.run<r>`, `.stale_*`, `.failed_*`) and held `.guard.lock` directories;
- a timing line for incomplete cells;
- a verdict: COMPLETE, INCOMPLETE (n missing), PROBLEM or NOT STARTED.

It was tested on a mock tree built from the real task files: complete cells, two race-skipped seeds, a TIMEOUT-like cell with partial seeds and a leftover `.run21`, and a PROBLEM cell with a zero-size `red`, a held lock and a second build record. Every case was reported as built.

#### W3 — resubmission

`--time` is set as above. The lines are printed by `round1_timing_261003.py` (verbatim, above); runsheet step 8f gives the KOA sequence. Set 1 (all six lines together) uses at most 64 cores. Set 2 (conf_A_0.39) waits for both A_0.10 resubmission jobs (`afterok`). When it starts, the two long B H40 tasks may still be running: 16 + 48 = 64 cores.

#### W4 — pull without rebuild

The repair touches no build input. Runsheet rule 4 now allows a pull that touches no `*.c`, `*.h`, `Makefile`, `edmd_core/` or `kissfft` to be followed by a resubmission without a rebuild, because the arrays verify the binary by sha256. Rule 5 forbids a rebuild while any cell is incomplete.


---

## 2. RESULTS (2026-10-04) — the pre-registered analysis of § 1 (with C1–C3), applied once

### Plain summary (for a non-specialist)

Paper 1 found that the speed of sound in its small simulated box comes out about 1 % above the value for an infinitely large gas. This campaign asked why. It changed the box's height and its length separately, and it measured the gas in each box in two independent ways:
- **Swinging divider.** A movable divider between two gas compartments is left free to swing, and its ringing frequency gives the speed of sound.
- **Clamped divider.** The divider is clamped, and the gas's push on it is measured at slightly different positions, which gives the gas's springiness directly.

All 38 measurement series finished. One of 12,545 simulation runs was set aside by the pre-registered health rule; it cannot have changed any result.

What was found:

1. **Dilute gas (η = 0.10).**
   - The 1 % excess comes from the walls that run along the direction of the sound: it grows about in inverse proportion to the box's height, and it barely changes when the box is made longer or shorter.
   - Two of the three pre-registered explanations pass the test: A ("all four walls") and B ("the long walls only"). B fits clearly better, but the pre-registered rule does not separate them.
   - The third, C (a shift caused by the divider's own thermal jiggling), is ruled out as the cause of the excess.
2. **Dense gas (η = π/8).**
   - None of the three explanations fits. The excess again grows as the box gets thinner, but it also depends on the box's length in a way none of them allows: the shortest box gives a sound speed below the infinite-gas value.
   - For the dense gas the question is **not resolved**.
3. **The cross-check between the two methods failed its pre-registered test.**
   - Most of the failure in the dense gas comes from a flaw in the clamped measurement, found afterwards. The "clamped" divider was not perfectly fixed: while the force was being recorded, it slowly slid back towards the middle. That made the gas look less springy than it is, by up to 35 % in the tallest box.
   - How far it slid can be read off the gas temperatures, which every run recorded. Correcting for it, a step that was not pre-registered, brings the two methods into agreement to within about 1.6 %.
   - A small difference remains in both densities. It is largest in the boxes with the fewest particles, which is the pattern pre-registered for explanation C, at roughly 60–75 % of its predicted size.

### 2.1 Data, gates and exclusion (X1, X2)

**Fetched [DATA].** Everything was fetched by `fetch_confinement.sh` on 2026-10-03, build `279282b target koa`:
- Method B: 513 files (19 cells × 9 masses × `red_nu.csv`, `acf_runs.npz`, `run.log`; 22 MB).
- Method A: 24,810 files (8,270 seeds × `red_`, `run_`, `summary_`; 20 MB).
- The two full pilot cells: `koa_pi8_H10_L10`, 133 files; the A pilot with its traces.

Raw trajectories stay on KOA scratch (decision, § 1.11). Every number below is printed by `python3 hspist3/validation/paper1_confinement_results_261004.py`, verbatim, unless marked post-hoc.

**Gates [DATA].**
- **Inventory: PASS.** Every cell is complete, and its geometry as recorded by the binary matches the task file and the registration (A summaries: L₀, H, 2N_s, t = 0.05, box width = 2L₀, η, wall position; B logs: `Initial wall_x` = 200 + 24 L₀ px in all 4,274 used runs).
  - Disclosure: the first run flagged 5,282 A summaries as different. The cause was print precision, not geometry: the summary prints L₀, H and the wall position to 4 decimals (`'19.7917'`), so values on the 1/24 grid differ from the task file by 3.3×10⁻⁵, while my tolerance was 10⁻⁵.
  - The tolerance was set to the print precision (5×10⁻⁵), and nothing else changed. A rerun reproduced every analysis number and both result CSVs byte for byte.
  - The B logs also print a startup banner `L0 (half-length): 20.000 σ` for every cell. It is the global default `L0_UNITS = 20.0f` (`00ALLINONE.c:250`), printed before the experiment loop sets `L0_UNITS = L0` per run (`:15732`); the per-run `Running: L0 = …` and `Initial wall_x` lines are right [DERIVATION].
- **Reduction gate (§ 1.10): PASS.** On the KOA pilot traces, `reduce_B.py` equals the canonical `cell()` to 7.3×10⁻¹⁷.
- **Determinism gate (§ 1.12): PASS.** 20 of 20 `red_970[0-3].csv` of the anchor cell are byte-identical to the pilot's: same source, different build hash, different job and node.
- **Health [DATA].** 0 health lines in every used B `run.log` and A `run_<seed>.log`.
- **Carried over from § 1.10:** the C2 mode-equivalence gate PASSED before launch. The pictures gate was met for two held-divider cells only, with no speed-of-sound pictures.

##### X1 -- inventory of the fetched summaries (expected = task file)

| eta | cell | B masses | B trajectories (exp.) | B files, MB | B nu rows with n missing | B wall_x = 200 + 24 L0 px | A positions x seeds (exp.) | A files, MB | A window min..max | health lines (B run.log / A run logs) | A build | A geometry (L0, H, 2N_s, t, box, eta, x_wall) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.10 | e0p10_H_H5_L39.25 | 9 | 225 (225) | 27, 1.18 | 0 | 225/225 | 5 x 142 (710) | 2130, 1.81 | 4999.4..5000.0 | 0 / 0 | 279282b | all as registered |
| 0.10 | e0p10_H_H10_L39.25 | 9 | 225 (225) | 27, 1.18 | 0 | 225/225 | 5 x 144 (720) | 2160, 1.84 | 4999.2..5000.0 | 0 / 0 | 279282b | all as registered |
| 0.10 | e0p10_H_H20_L39.25 | 9 | 225 (225) | 27, 1.19 | 0 | 225/225 | 5 x 135 (675) | 2025, 1.73 | 4999.5..5000.0 | 0 / 0 | 279282b | all as registered |
| 0.10 | e0p10_H_H40_L39.25 | 9 | 225 (225) | 27, 1.19 | 0 | 225/225 | 5 x 144 (720) | 2160, 1.84 | 4999.6..5000.0 | 0 / 0 | 279282b | all as registered |
| 0.10 | e0p10_L_H10_L19.625 | 9 | 225 (225) | 27, 1.18 | 0 | 225/225 | 5 x 150 (750) | 2250, 1.91 | 4999.3..5000.0 | 0 / 0 | 279282b | all as registered |
| 0.10 | e0p10_L_H10_L78.5 | 9 | 225 (225) | 27, 1.19 | 0 | 225/225 | 5 x 135 (675) | 2025, 1.73 | 4999.7..5000.0 | 0 / 0 | 279282b | all as registered |
| 0.10 | e0p10_aspect_H19.7917_L19.7917 | 9 | 225 (225) | 27, 1.19 | 0 | 225/225 | 5 x 146 (730) | 2190, 1.89 | 4999.3..5000.0 | 0 / 0 | 279282b | all as registered |
| 0.10 | e0p10_aspect_H14_L28 | 9 | 225 (225) | 27, 1.18 | 0 | 225/225 | 5 x 137 (685) | 2055, 1.75 | 4999.4..5000.0 | 0 / 0 | 279282b | all as registered |
| 0.10 | e0p10_aspect_H9.91667_L39.625 | 9 | 225 (225) | 27, 1.19 | 0 | 225/225 | 5 x 136 (680) | 2040, 1.76 | 4999.5..5000.0 | 0 / 0 | 279282b | all as registered |
| 0.10 | e0p10_aspect_H7_L56.0417 | 9 | 225 (225) | 27, 1.18 | 0 | 225/225 | 5 x 137 (685) | 2055, 1.76 | 4999.6..5000.0 | 0 / 0 | 279282b | all as registered |
| 0.39 | epi8_H_H5_L10 | 9 | 225 (225) | 27, 1.20 | 0 | 225/225 | 5 x 26 (130) | 390, 0.33 | 4999.9..5000.0 | 0 / 0 | 279282b | all as registered |
| 0.39 | epi8_H_H10_L10 | 9 | 225 (225) | 27, 1.21 | 0 | 225/225 | 5 x 23 (115) | 345, 0.29 | 4999.9..5000.0 | 0 / 0 | 279282b | all as registered |
| 0.39 | epi8_H_H20_L10 | 9 | 225 (225) | 27, 1.21 | 0 | 225/225 | 5 x 26 (130) | 390, 0.33 | 4999.9..5000.0 | 0 / 0 | 279282b | all as registered |
| 0.39 | epi8_H_H40_L10 | 9 | 224 (225) | 27, 1.21 | 0 | 224/224 | 5 x 51 (255) | 765, 0.65 | 4999.9..5000.0 | 0 / 0 | 279282b | all as registered |
| 0.39 | epi8_L_H10_L5 | 9 | 225 (225) | 27, 1.21 | 0 | 225/225 | 5 x 26 (130) | 390, 0.33 | 4999.8..5000.0 | 0 / 0 | 279282b | all as registered |
| 0.39 | epi8_L_H10_L20 | 9 | 225 (225) | 27, 1.21 | 0 | 225/225 | 5 x 26 (130) | 390, 0.33 | 4999.9..5000.0 | 0 / 0 | 279282b | all as registered |
| 0.39 | epi8_aspect_H7.08333_L14.125 | 9 | 225 (225) | 27, 1.21 | 0 | 225/225 | 5 x 26 (130) | 390, 0.34 | 4999.9..5000.0 | 0 / 0 | 279282b | all as registered |
| 0.39 | epi8_aspect_H5_L20 | 9 | 225 (225) | 27, 1.21 | 0 | 225/225 | 5 x 23 (115) | 345, 0.29 | 4999.9..5000.0 | 0 / 0 | 279282b | all as registered |
| 0.39 | epi8_aspect_H3.54167_L28.2917 | 9 | 225 (225) | 27, 1.21 | 0 | 225/225 | 5 x 21 (105) | 315, 0.27 | 4999.9..5000.0 | 0 / 0 | 279282b | all as registered |

inventory: every cell complete and as registered, except the one excluded B trajectory

##### Reduction gate (sec. 1.10): reduce_B.py vs the canonical cell() on the full KOA pilot cell (pi/8 anchor, 1 run per mass)

| M | nu, reduce_B.py (red_nu.csv) | nu, cell() | abs. difference |
|---|---|---|---|
| 50 | 0.072116205132757294 | 0.072116205132757294 | 0.0e+00 |
| 100 | 0.058980134293145602 | 0.058980134293145671 | 6.9e-17 |
| 200 | 0.044370262967487598 | 0.044370262967487598 | 0.0e+00 |
| 300 | 0.037163207870788702 | 0.037163207870788709 | 6.9e-18 |
| 500 | 0.029535941825679202 | 0.029535941825679275 | 7.3e-17 |
| 750 | 0.024149659505260102 | 0.024149659505260157 | 5.6e-17 |
| 1000 | 0.021028002470520901 | 0.021028002470520973 | 7.3e-17 |
| 1500 | 0.017425438449665001 | 0.017425438449665025 | 2.4e-17 |
| 2000 | 0.0150619720497352 | 0.015061972049735238 | 3.8e-17 |

reduction gate: max abs. difference 7.3e-17 -> PASS

##### Determinism gate (sec. 1.12): anchor cell of conf_A_0.39 vs the pilot, red_970[0-3].csv, cmp

| position | seed | anchor bytes | pilot bytes | cmp |
|---|---|---|---|---|
| x_m2 | 9700 | 145 | 145 | IDENTICAL |
| x_m2 | 9701 | 144 | 144 | IDENTICAL |
| x_m2 | 9702 | 147 | 147 | IDENTICAL |
| x_m2 | 9703 | 147 | 147 | IDENTICAL |
| x_m1 | 9700 | 146 | 146 | IDENTICAL |
| x_m1 | 9701 | 146 | 146 | IDENTICAL |
| x_m1 | 9702 | 143 | 143 | IDENTICAL |
| x_m1 | 9703 | 145 | 145 | IDENTICAL |
| x_0 | 9700 | 147 | 147 | IDENTICAL |
| x_0 | 9701 | 144 | 144 | IDENTICAL |
| x_0 | 9702 | 146 | 146 | IDENTICAL |
| x_0 | 9703 | 147 | 147 | IDENTICAL |
| x_p1 | 9700 | 147 | 147 | IDENTICAL |
| x_p1 | 9701 | 147 | 147 | IDENTICAL |
| x_p1 | 9702 | 147 | 147 | IDENTICAL |
| x_p1 | 9703 | 146 | 146 | IDENTICAL |
| x_p2 | 9700 | 147 | 147 | IDENTICAL |
| x_p2 | 9701 | 147 | 147 | IDENTICAL |
| x_p2 | 9702 | 147 | 147 | IDENTICAL |
| x_p2 | 9703 | 146 | 146 | IDENTICAL |

determinism gate: 20/20 IDENTICAL -> PASS
max |eta_rec - eta_reg| = 4.2e-08, max |dL - dL_reg| = 6.7e-07

geometry vs registration (eta, dL from paper1_confinement_prereg_20261012): all equal; box shortfall delta: max 2.54e-06 sigma (grid-exact boxes)

**The excluded trajectory (X2).** It is `epi8_H_H40_L10`, M = 3000 (α = 7.5), run 5. It ended with rc = 0 and health = 1.
- **The rule [SOURCE, § 1.7 item 6]:** "**Health contract** zero on every run (forced_advance, clamp_repair, overlap_repair, wall_overdue)."
- **How it was applied [SOURCE, `conf_worker.sh` mode B]:** the run was set aside before its trace entered the cell:

      if [ "$rc" -ne 0 ] || [ "${h:-0}" -ne 0 ] || [ -z "$tr" ]; then
        echo "B $rel M=$M r=$r FAILED rc=$rc health=${h:-0}"; mv "$tmp" "$cell/.failed_run${r}_$(date +%Y%m%d_%H%M%S)"; exit 1; fi

  The canonical estimator would have discarded it as well (`tests_20260913.cell_runs`, `r in bad`).
- **The health line itself is OPEN.** It is in `.failed_run5_*/stdout.log` on KOA scratch, which the fetch filter did not copy.
- **What each counter means [SOURCE, `edmd_core/edmd.c`]:**
  - `forced_advance` (`:1457`, `S->forced_advance_count++;`): the event loop exceeded its event or stagnation guard, and time was advanced by force.
  - `wall_clamp_repairs` (`:365`, `S->clamp_repair_count++;`): `grid_build()` had to bounce a particle back into the box.
  - `overlap_repairs` (`:760`, `if(ok == 2) S->overlap_repair_count++;`): an already-overlapping approaching pair had to be rescued.
  - `wall_overdue` (`:589–601`, `if(rc==2) S->wall_overdue_count++;`): a wall collision was already overdue when it was scheduled.
- **Precedent [SOURCE, 260913 STATUS]:** A1 v2 had 2 health lines in 7,875 runs, both `wall_clamp_repairs=1`, both discarded.
- **What the exclusion could have moved [DATA]:**

##### X2 -- the excluded trajectory: how far could it have moved the cell?

cell epi8_H_H40_L10, M = 3000 (alpha = 7.5): 24 of 25 runs used; nu = 0.0236678, sd over seeds = 0.0001899
- a 25th run largest deviation among the 24 used runs (0.0004819) away from the mean moves that mass's nu by 1.93e-05 and c_s by 1.23e-04 = 0.007 of the cell's c_s_err_scaled (0.01721)
- a 25th run 3 sd (0.0005696) away from the mean moves that mass's nu by 2.28e-05 and c_s by 1.46e-04 = 0.008 of the cell's c_s_err_scaled (0.01721)
- the same mass is one of the five heavy masses of the identity: its k_S^dyn = 34.64592; a 3-sd 25th run changes it by 6.67e-02; with its inverse-variance weight 0.168 the cell's combined k_S^dyn moves by 1.12e-02 = 0.24 of its sigma (4.66e-02) and rho_I by 0.033 %

### 2.2 Confinement, method B: the shift Δ = c_s/c_s^KR − 1 and the A/B/C verdicts (§ 1.5)

##### Confinement: per-cell shift from method B (canonical estimator, campaign cells only)

| eta | cell | scan | H | L_0 | N_s | eta_true | c_s | c_s_err | chi2_red (9 masses) | c_s_err_scaled | c_s^KR(eta_true) | Delta = c_s/c_s^KR - 1 [%] | sigma [%] | shape A: 2/H + 2/L_0 | shape B: 1/H | shape C: 1/N_s | Delta_C (fixed) [%] |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.10 | e0p10_H_H5_L39.25 | H | 5 | 39.25 | 25 | 0.100051 | 1.78685 | 0.00345 | 0.88 | 0.00345 | 1.74434 | +2.437 | 0.198 | 0.4510 | 0.2000 | 0.0400 | +1.269 |
| 0.10 | e0p10_H_H10_L39.25 | H | 10 | 39.25 | 50 | 0.100051 | 1.76165 | 0.00186 | 3.19 | 0.00331 | 1.74434 | +0.992 | 0.190 | 0.2510 | 0.1000 | 0.0200 | +0.634 |
| 0.10 | e0p10_H_H20_L39.25 | H | 20 | 39.25 | 100 | 0.100051 | 1.74851 | 0.00260 | 10.14 | 0.00827 | 1.74434 | +0.239 | 0.474 | 0.1510 | 0.0500 | 0.0100 | +0.317 |
| 0.10 | e0p10_H_H40_L39.25 | H | 40 | 39.25 | 200 | 0.100051 | 1.74514 | 0.00203 | 3.08 | 0.00356 | 1.74434 | +0.046 | 0.204 | 0.1010 | 0.0250 | 0.0050 | +0.159 |
| 0.10 | e0p10_L_H10_L19.625 | L | 10 | 19.625 | 25 | 0.100051 | 1.76502 | 0.00352 | 1.02 | 0.00355 | 1.74434 | +1.186 | 0.203 | 0.3019 | 0.1000 | 0.0400 | +1.269 |
| 0.10 | e0p10_L_H10_L78.5 | L | 10 | 78.5 | 100 | 0.100051 | 1.76551 | 0.00160 | 3.64 | 0.00306 | 1.74434 | +1.214 | 0.175 | 0.2255 | 0.1000 | 0.0100 | +0.317 |
| 0.10 | e0p10_aspect_H19.7917_L19.7917 | aspect | 19.7917 | 19.7917 | 50 | 0.100252 | 1.75723 | 0.00385 | 0.84 | 0.00385 | 1.74511 | +0.695 | 0.221 | 0.2021 | 0.0505 | 0.0200 | +0.634 |
| 0.10 | e0p10_aspect_H14_L28 | aspect | 14 | 28 | 50 | 0.100178 | 1.76267 | 0.00319 | 1.90 | 0.00439 | 1.74483 | +1.022 | 0.252 | 0.2143 | 0.0714 | 0.0200 | +0.634 |
| 0.10 | e0p10_aspect_H9.91667_L39.625 | aspect | 9.91667 | 39.625 | 50 | 0.099937 | 1.76276 | 0.00269 | 1.69 | 0.00349 | 1.74390 | +1.081 | 0.200 | 0.2522 | 0.1008 | 0.0200 | +0.634 |
| 0.10 | e0p10_aspect_H7_L56.0417 | aspect | 7 | 56.0417 | 50 | 0.100104 | 1.77716 | 0.00266 | 5.89 | 0.00646 | 1.74454 | +1.870 | 0.370 | 0.3214 | 0.1429 | 0.0200 | +0.634 |
| 0.39 | epi8_H_H5_L10 | H | 5 | 10 | 25 | 0.392699 | 3.90716 | 0.00901 | 6.20 | 0.02243 | 3.74608 | +4.300 | 0.599 | 0.6000 | 0.2000 | 0.0400 | +0.768 |
| 0.39 | epi8_H_H10_L10 | H | 10 | 10 | 50 | 0.392699 | 3.82030 | 0.00671 | 3.41 | 0.01239 | 3.74608 | +1.981 | 0.331 | 0.4000 | 0.1000 | 0.0200 | +0.384 |
| 0.39 | epi8_H_H20_L10 | H | 20 | 10 | 100 | 0.392699 | 3.75811 | 0.00732 | 1.54 | 0.00909 | 3.74608 | +0.321 | 0.243 | 0.3000 | 0.0500 | 0.0100 | +0.192 |
| 0.39 | epi8_H_H40_L10 | H | 40 | 10 | 200 | 0.392699 | 3.71657 | 0.00727 | 5.60 | 0.01721 | 3.74608 | -0.788 | 0.459 | 0.2500 | 0.0250 | 0.0050 | +0.096 |
| 0.39 | epi8_L_H10_L5 | L | 10 | 5 | 25 | 0.392699 | 3.68187 | 0.01305 | 3.70 | 0.02508 | 3.74608 | -1.714 | 0.670 | 0.6000 | 0.1000 | 0.0400 | +0.768 |
| 0.39 | epi8_L_H10_L20 | L | 10 | 20 | 100 | 0.392699 | 3.85055 | 0.00490 | 1.24 | 0.00547 | 3.74608 | +2.789 | 0.146 | 0.3000 | 0.1000 | 0.0100 | +0.192 |
| 0.39 | epi8_aspect_H7.08333_L14.125 | aspect | 7.08333 | 14.125 | 50 | 0.392495 | 3.87945 | 0.00554 | 1.13 | 0.00590 | 3.74367 | +3.627 | 0.158 | 0.4239 | 0.1412 | 0.0200 | +0.384 |
| 0.39 | epi8_aspect_H5_L20 | aspect | 5 | 20 | 50 | 0.392699 | 3.98277 | 0.00539 | 0.91 | 0.00539 | 3.74608 | +6.318 | 0.144 | 0.5000 | 0.2000 | 0.0200 | +0.384 |
| 0.39 | epi8_aspect_H3.54167_L28.2917 | aspect | 3.54167 | 28.2917 | 50 | 0.391917 | 3.99919 | 0.00455 | 0.30 | 0.00455 | 3.73688 | +7.020 | 0.122 | 0.6354 | 0.2824 | 0.0200 | +0.384 |

###### eta 0.10: one-amplitude fits over its 10 cells (sec. 1.5 decision rule; excluded if p < 0.01)

| hypothesis | shape | amplitude | chi2 | dof | p(chi2) | excluded? |
|---|---|---|---|---|---|---|
| A | a (2/H + 2/L_0) | 0.04674 +- 0.00258 | 13.56 | 9 | 0.139 | no |
| B | b/H | 0.11797 +- 0.00642 | 4.45 | 9 | 0.88 | no |
| C | c'/N_s | 0.49190 +- 0.02859 | 46.30 | 9 | 5.3e-07 | YES |
| C_fixed | Delta_C, no free parameter | fixed | 83.65 | 10 | 9.62e-14 | YES |

**Verdict, eta 0.10: not separated: A, B survive; Delta chi2 to the best (B): A +9.11.**

###### eta 0.39: one-amplitude fits over its 9 cells (sec. 1.5 decision rule; excluded if p < 0.01)

| hypothesis | shape | amplitude | chi2 | dof | p(chi2) | excluded? |
|---|---|---|---|---|---|---|
| A | a (2/H + 2/L_0) | 0.10159 +- 0.00134 | 480.97 | 8 | 8.51e-99 | YES |
| B | b/H | 0.26256 +- 0.00336 | 141.94 | 8 | 9.38e-27 | YES |
| C | c'/N_s | 2.54736 +- 0.03521 | 1018.12 | 8 | 1.83e-214 | YES |
| C_fixed | Delta_C, no free parameter | fixed | 5493.24 | 9 | 0 | YES |

**Verdict, eta 0.39: none survives.**

Exploratory two-term forms (registered as exploratory only, not a verdict):

| form | first amplitude | c' | chi2 | dof | p(chi2) |
|---|---|---|---|---|---|
| b/H + c'/N_s | 0.29114 +- 0.00978 | -0.31900 +- 0.10254 | 132.26 | 7 | 2.12e-25 |
| a(2/H + 2/L_0) + c'/N_s | 0.15478 +- 0.00625 | -1.43370 +- 0.16460 | 405.10 | 7 | 1.92e-83 |

**Reading.**
- **η = 0.10.** A (p = 0.14) and B (p = 0.88) survive; C is excluded (p = 5×10⁻⁷, and 10⁻¹³ at its fixed amplitude) [DATA]. **Registered verdict: not separated.** Δχ² = +9.1 for A over B [DATA]. B is preferred by that Δχ² (likelihood ratio ≈ e^4.6), but the registered rule does not declare it [INFERENCE].
  - The η = 0.10 L-scan is flat (+1.19, +0.99, +1.21 % at L₀ = 19.6, 39.25, 78.5) [DATA]. That is B's shape, and it is also why C (which would be 1/N_s, i.e. 1/L₀ here) is excluded [DERIVATION].
  - The campaign's own anchor reads +0.99 ± 0.19 %; Paper 1's A1 v2 value, quoted only and never fitted (binary rule, § 1.5), is +1.01 % [DATA].
- **η = π/8.** A, B and C, and C at its fixed amplitude, are all excluded (p ≤ 10⁻²⁶). The two exploratory two-term forms fail as well (p = 2×10⁻²⁵ and 2×10⁻⁸³) [DATA]. **Registered outcome: not resolved.**
  - The anchor reads +1.98 ± 0.33 %; Paper 1's +1.675 % is quoted only [DATA].
  - What no form contains [DATA]: at H = 10 the L-scan runs opposite to A, with −1.71 ± 0.67 %, +1.98 ± 0.33 % and +2.79 ± 0.15 % at L₀ = 5, 10, 20. The H-scan falls below zero at H = 40 (−0.79 ± 0.46 %). The aspect cells rise to +7.0 % at L₀/H = 8.
  - [INFERENCE] The sign and the 1/L shape of the π/8 L-dependence are those of an offset in the length of the frequency formula, ν = c_s K/(2π L_eff), at high density; the length-free identity (§ 2.3) does not see it.
  - **What would resolve it [OPEN]:**
    - compare the length-free stiffness of every cell (method B, k_S^dyn) with the bulk KR stiffness, and read off the effective length that makes them agree, cell by cell, with no new runs;
    - a finer L-scan at π/8 (L₀ = 5 … 40 at H = 10);
    - the wall-contact density profile from the method-A event logs on KOA scratch.

![confinement shift](../paper1_speedofsound/experiments/final/261004_p1_confinement_shift.png)

`paper1_speedofsound/experiments/final/261004_p1_confinement_shift.png/.pdf`
- Rows: η = 0.10 and π/8. Columns: the H-scan against 1/H, the L-scan against 1/L₀, and the aspect cells against L₀/H.
- Error bars are `c_s_err_scaled`. KR is drawn in red (Δ = 0), the data in blue.
- A, B and C are drawn with the free amplitudes of the fits above; C at its fixed amplitude is the thin line.

### 2.3 Method A and the identity (C1)

##### Method A per cell: F at L_0, k_T (5-point stencil, mirror faces averaged), kT, symmetry checks

| eta | cell | dL | seeds/position | F(L_0) | sigma_F | kT (x = 0 seeds) | k_T | sigma(k_T) | sigma(k_T)/k_T [%] | symmetry z: x=0; x=+1,+2,-1,-2 dL | |z| > 2 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.10 | e0p10_H_H5_L39.25 | 1.5417 | 142 | 0.82215 | 0.00029 | 1.000000 | 0.02638 | 0.00018 | 0.69 | +0.36; +1.37, -0.13, +0.66, +0.02 | 0 |
| 0.10 | e0p10_H_H10_L39.25 | 1.0833 | 144 | 1.63150 | 0.00037 | 1.000000 | 0.05272 | 0.00032 | 0.61 | -1.50; +1.12, -0.29, +0.84, +1.00 | 0 |
| 0.10 | e0p10_H_H20_L39.25 | 0.7917 | 135 | 3.24996 | 0.00052 | 1.000000 | 0.10446 | 0.00065 | 0.62 | -2.46; +0.88, -0.32, -0.25, +0.67 | 1 |
| 0.10 | e0p10_H_H40_L39.25 | 0.5417 | 144 | 6.48580 | 0.00068 | 1.000000 | 0.20819 | 0.00121 | 0.58 | -0.68; -0.38, -0.08, +0.34, -0.84 | 0 |
| 0.10 | e0p10_L_H10_L19.625 | 0.7500 | 150 | 1.67176 | 0.00051 | 1.000000 | 0.10926 | 0.00063 | 0.57 | +0.32; +0.04, -1.23, +1.02, +1.71 | 0 |
| 0.10 | e0p10_L_H10_L78.5 | 1.5833 | 135 | 1.61202 | 0.00031 | 1.000000 | 0.02575 | 0.00017 | 0.67 | +0.85; +1.87, +1.08, +0.01, -1.73 | 0 |
| 0.10 | e0p10_aspect_H19.7917_L19.7917 | 0.5417 | 146 | 3.30693 | 0.00070 | 1.000000 | 0.21204 | 0.00122 | 0.58 | +0.94; -1.26, +0.88, +1.57, +1.29 | 0 |
| 0.10 | e0p10_aspect_H14_L28 | 0.7917 | 137 | 2.30551 | 0.00052 | 1.000000 | 0.10472 | 0.00063 | 0.60 | +0.96; +1.05, -0.40, -0.16, +0.28 | 0 |
| 0.10 | e0p10_aspect_H9.91667_L39.625 | 1.1250 | 136 | 1.61620 | 0.00038 | 1.000000 | 0.05064 | 0.00033 | 0.65 | +0.24; -0.81, +0.10, -0.28, +0.99 | 0 |
| 0.10 | e0p10_aspect_H7_L56.0417 | 1.5833 | 137 | 1.13947 | 0.00027 | 1.000000 | 0.02554 | 0.00016 | 0.64 | -0.66; -0.75, +0.87, +1.87, +2.19 | 1 |
| 0.39 | epi8_H_H5_L10 | 0.1667 | 26 | 8.21301 | 0.00238 | 1.000000 | 1.96591 | 0.01348 | 0.69 | +0.59; +1.52, +0.03, -0.02, -0.39 | 0 |
| 0.39 | epi8_H_H10_L10 | 0.1250 | 23 | 15.92514 | 0.00359 | 1.000000 | 3.62779 | 0.02583 | 0.71 | +0.04; +1.24, +2.73, +1.85, +0.18 | 1 |
| 0.39 | epi8_H_H20_L10 | 0.0833 | 26 | 31.34194 | 0.00391 | 1.000000 | 6.66234 | 0.05131 | 0.77 | -0.56; -1.82, -0.01, +0.99, +0.46 | 0 |
| 0.39 | epi8_H_H40_L10 | 0.0417 | 51 | 62.17875 | 0.00433 | 1.000000 | 11.26896 | 0.10122 | 0.90 | +0.23; -0.63, -2.50, -0.43, +0.25 | 1 |
| 0.39 | epi8_L_H10_L5 | 0.0833 | 26 | 17.56306 | 0.00468 | 1.000000 | 7.45739 | 0.05361 | 0.72 | -0.44; +0.45, -0.68, -0.45, +1.54 | 0 |
| 0.39 | epi8_L_H10_L20 | 0.1667 | 26 | 15.13641 | 0.00236 | 1.000000 | 1.74548 | 0.01324 | 0.76 | +0.78; -0.70, +0.24, -0.06, -0.16 | 0 |
| 0.39 | epi8_aspect_H7.08333_L14.125 | 0.1667 | 26 | 11.09669 | 0.00216 | 1.000000 | 1.86791 | 0.01361 | 0.73 | +0.79; -1.49, -1.62, -1.20, +0.80 | 0 |
| 0.39 | epi8_aspect_H5_L20 | 0.2500 | 23 | 7.84468 | 0.00206 | 1.000000 | 0.93416 | 0.00696 | 0.74 | -0.12; -0.81, -0.81, -0.82, -1.66 | 0 |
| 0.39 | epi8_aspect_H3.54167_L28.2917 | 0.3750 | 21 | 5.54406 | 0.00172 | 1.000000 | 0.45143 | 0.00333 | 0.74 | +0.04; -0.09, -1.94, +1.14, +1.04 | 0 |

symmetry checks beyond 2 sigma: 4 of 95 (registered: each within 2 sigma; expected by chance if all hold: 4.3)

The registered rule "each [symmetry check] within 2 σ" is formally violated by 4 of 95 checks. That is the rate chance alone gives (4.3 expected) [DATA, DERIVATION]. The F² term's own relative error is at most 0.071 % (the registration said < 0.05 %; it is neglected, as registered) [DATA].

##### Identity, length-free form (amendment C1): k_S^dyn (heavy masses) vs k_T + F^2/(N_s kT)

| eta | cell | N_s | k_S^dyn | sigma | chi2_red of the 5 heavy masses | k_T | F^2/(N_s kT) | its own sigma (neglected) | static = k_T + F^2/(N_s kT) | rho_I [%] | sigma(rho_I) [%] | rho_I/sigma | within 2 sigma? | 2 Delta_C (hyp. C) [%] |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.10 | e0p10_H_H5_L39.25 | 25 | 0.054596 | 0.000081 | 1.22 | 0.026383 | 0.027037 | 1.93e-05 | 0.053420 | +2.155 | 0.367 | +5.87 | **NO** | +2.537 |
| 0.10 | e0p10_H_H10_L39.25 | 50 | 0.106551 | 0.000108 | 0.21 | 0.052725 | 0.053236 | 2.41e-05 | 0.105961 | +0.554 | 0.319 | +1.73 | yes | +1.269 |
| 0.10 | e0p10_H_H20_L39.25 | 100 | 0.210768 | 0.000212 | 0.38 | 0.104463 | 0.105622 | 3.40e-05 | 0.210085 | +0.324 | 0.322 | +1.01 | yes | +0.634 |
| 0.10 | e0p10_H_H40_L39.25 | 200 | 0.418797 | 0.000451 | 0.55 | 0.208185 | 0.210328 | 4.40e-05 | 0.418513 | +0.068 | 0.308 | +0.22 | yes | +0.317 |
| 0.10 | e0p10_L_H10_L19.625 | 25 | 0.225581 | 0.000367 | 0.17 | 0.109264 | 0.111792 | 6.83e-05 | 0.221056 | +2.006 | 0.322 | +6.24 | **NO** | +2.537 |
| 0.10 | e0p10_L_H10_L78.5 | 100 | 0.051781 | 0.000037 | 1.20 | 0.025753 | 0.025986 | 9.90e-06 | 0.051739 | +0.081 | 0.340 | +0.24 | yes | +0.634 |
| 0.10 | e0p10_aspect_H19.7917_L19.7917 | 50 | 0.437882 | 0.000646 | 1.20 | 0.212043 | 0.218716 | 9.31e-05 | 0.430758 | +1.627 | 0.316 | +5.15 | **NO** | +1.268 |
| 0.10 | e0p10_aspect_H14_L28 | 50 | 0.212795 | 0.000262 | 1.26 | 0.104717 | 0.106307 | 4.84e-05 | 0.211024 | +0.832 | 0.319 | +2.61 | **NO** | +1.268 |
| 0.10 | e0p10_aspect_H9.91667_L39.625 | 50 | 0.104417 | 0.000092 | 1.12 | 0.050643 | 0.052242 | 2.47e-05 | 0.102885 | +1.467 | 0.326 | +4.51 | **NO** | +1.269 |
| 0.10 | e0p10_aspect_H7_L56.0417 | 50 | 0.051878 | 0.000049 | 2.10 | 0.025541 | 0.025968 | 1.22e-05 | 0.051509 | +0.712 | 0.327 | +2.18 | **NO** | +1.269 |
| 0.39 | epi8_H_H5_L10 | 25 | 4.772302 | 0.005708 | 1.37 | 1.965906 | 2.698141 | 1.56e-03 | 4.664047 | +2.268 | 0.307 | +7.39 | **NO** | +1.536 |
| 0.39 | epi8_H_H10_L10 | 50 | 9.006886 | 0.011494 | 0.92 | 3.627794 | 5.072205 | 2.28e-03 | 8.699998 | +3.407 | 0.314 | +10.85 | **NO** | +0.768 |
| 0.39 | epi8_H_H20_L10 | 100 | 17.503565 | 0.020399 | 1.67 | 6.662338 | 9.823171 | 2.45e-03 | 16.485510 | +5.816 | 0.315 | +18.44 | **NO** | +0.384 |
| 0.39 | epi8_H_H40_L10 | 200 | 34.561988 | 0.046576 | 1.58 | 11.268959 | 19.330988 | 2.69e-03 | 30.599947 | +11.464 | 0.322 | +35.56 | **NO** | +0.192 |
| 0.39 | epi8_L_H10_L5 | 25 | 21.654516 | 0.051033 | 1.18 | 7.457385 | 12.338453 | 6.57e-03 | 19.795838 | +8.583 | 0.342 | +25.11 | **NO** | +1.536 |
| 0.39 | epi8_L_H10_L20 | 100 | 4.127295 | 0.003957 | 0.37 | 1.745481 | 2.291110 | 7.14e-04 | 4.036591 | +2.198 | 0.335 | +6.56 | **NO** | +0.384 |
| 0.39 | epi8_aspect_H7.08333_L14.125 | 50 | 4.386889 | 0.004730 | 1.38 | 1.867908 | 2.462732 | 9.57e-04 | 4.330640 | +1.282 | 0.329 | +3.90 | **NO** | +0.768 |
| 0.39 | epi8_aspect_H5_L20 | 50 | 2.202577 | 0.001810 | 1.58 | 0.934164 | 1.230779 | 6.45e-04 | 2.164943 | +1.709 | 0.326 | +5.24 | **NO** | +0.768 |
| 0.39 | epi8_aspect_H3.54167_L28.2917 | 50 | 1.075638 | 0.000900 | 0.07 | 0.451427 | 0.614731 | 3.82e-04 | 1.066158 | +0.881 | 0.321 | +2.75 | **NO** | +0.769 |

**Identity verdict (C1 rule: agreement within 2 sigma at every cell): FAIL** -- 4 of 19 cells within 2 sigma; outside: e0p10_H_H5_L39.25, e0p10_L_H10_L19.625, e0p10_aspect_H19.7917_L19.7917, e0p10_aspect_H14_L28, e0p10_aspect_H9.91667_L39.625, e0p10_aspect_H7_L56.0417, epi8_H_H5_L10, epi8_H_H10_L10, epi8_H_H20_L10, epi8_H_H40_L10, epi8_L_H10_L5, epi8_L_H10_L20, epi8_aspect_H7.08333_L14.125, epi8_aspect_H5_L20, epi8_aspect_H3.54167_L28.2917.
(information, not the registered rule: sum of (rho_I/sigma)^2 over the 19 cells = 2636.6, p = 0; P(all 19 within 2 sigma | identity exact) = 0.41)

**Registered verdict: FAIL.** 4 of 19 cells are within 2 σ [DATA].
- At π/8 every cell fails, with ρ_I = +0.9 … +11.5 %, growing with H/L₀.
- At η = 0.10 the failures are the small-N_s cells: N_s = 25 gives +2.0 and +2.2 %; the N_s = 50 aspect cells give +0.7 … +1.6 %. Every N_s ≥ 100 cell and the anchor pass.
- § 1.5 registered how a failure concentrated at small N_s is to be read: "C's signature, and is read that way." The π/8 pattern is not that signature; § 2.8 finds its main cause (post-hoc).

##### Standing-wave check at all alpha (C1, no verdict): k_S^SW(alpha)/k_S^dyn, k_S^SW = N_s m omega_1^2/K(alpha)^2

| eta | cell | alpha 0.5 | alpha 1 | alpha 2 | alpha 3 | alpha 5 | alpha 7.5 | alpha 10 | alpha 15 | alpha 20 |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.10 | e0p10_H_H5_L39.25 | 1.0001 | 0.9974 | 1.0032 | 1.0056 | 1.0067 | 0.9950 | 0.9993 | 1.0005 | 1.0013 |
| 0.10 | e0p10_H_H10_L39.25 | 0.9873 | 1.0027 | 1.0050 | 0.9999 | 1.0015 | 0.9999 | 0.9983 | 1.0007 | 1.0007 |
| 0.10 | e0p10_H_H20_L39.25 | 0.9809 | 0.9987 | 0.9998 | 1.0010 | 1.0014 | 0.9995 | 0.9983 | 0.9995 | 1.0014 |
| 0.10 | e0p10_H_H40_L39.25 | 0.9929 | 0.9933 | 0.9980 | 0.9972 | 1.0016 | 1.0028 | 0.9988 | 0.9980 | 0.9996 |
| 0.10 | e0p10_L_H10_L19.625 | 0.9881 | 1.0027 | 1.0125 | 0.9980 | 1.0019 | 1.0035 | 0.9987 | 0.9999 | 0.9995 |
| 0.10 | e0p10_L_H10_L78.5 | 1.0079 | 1.0004 | 1.0008 | 0.9973 | 1.0012 | 1.0024 | 0.9949 | 1.0004 | 0.9996 |
| 0.10 | e0p10_aspect_H19.7917_L19.7917 | 1.0062 | 0.9967 | 0.9963 | 0.9994 | 1.0014 | 1.0026 | 1.0050 | 0.9963 | 0.9978 |
| 0.10 | e0p10_aspect_H14_L28 | 1.0054 | 0.9975 | 1.0066 | 1.0113 | 0.9971 | 1.0050 | 0.9990 | 0.9983 | 0.9990 |
| 0.10 | e0p10_aspect_H9.91667_L39.625 | 0.9993 | 0.9938 | 0.9970 | 1.0053 | 1.0041 | 1.0009 | 1.0013 | 0.9970 | 0.9989 |
| 0.10 | e0p10_aspect_H7_L56.0417 | 1.0057 | 1.0061 | 1.0116 | 1.0030 | 1.0030 | 1.0026 | 1.0006 | 0.9959 | 1.0022 |
| 0.39 | epi8_H_H5_L10 | 0.9861 | 0.9911 | 1.0021 | 0.9948 | 1.0072 | 0.9970 | 0.9995 | 0.9977 | 1.0006 |
| 0.39 | epi8_H_H10_L10 | 1.0125 | 1.0046 | 0.9971 | 1.0061 | 1.0001 | 0.9997 | 1.0023 | 0.9948 | 1.0008 |
| 0.39 | epi8_H_H20_L10 | 0.9966 | 1.0123 | 0.9947 | 1.0023 | 1.0081 | 1.0050 | 0.9989 | 0.9989 | 0.9977 |
| 0.39 | epi8_H_H40_L10 | 0.9819 | 0.9998 | 0.9956 | 0.9914 | 1.0052 | 1.0028 | 1.0044 | 0.9979 | 0.9960 |
| 0.39 | epi8_L_H10_L5 | 0.9699 | 1.0036 | 1.0009 | 1.0110 | 0.9936 | 0.9932 | 1.0036 | 1.0048 | 1.0029 |
| 0.39 | epi8_L_H10_L20 | 0.9960 | 0.9971 | 0.9982 | 1.0012 | 1.0008 | 1.0001 | 0.9982 | 0.9994 | 1.0016 |
| 0.39 | epi8_aspect_H7.08333_L14.125 | 0.9935 | 1.0057 | 1.0039 | 1.0027 | 0.9968 | 0.9979 | 1.0005 | 1.0042 | 0.9994 |
| 0.39 | epi8_aspect_H5_L20 | 0.9979 | 1.0037 | 1.0008 | 0.9990 | 1.0004 | 0.9986 | 1.0030 | 0.9980 | 0.9986 |
| 0.39 | epi8_aspect_H3.54167_L28.2917 | 0.9973 | 1.0016 | 1.0023 | 1.0020 | 1.0005 | 1.0005 | 1.0009 | 1.0002 | 0.9995 |

The dynamic side is internally consistent: across the whole mass ladder, k_S^SW/k_S^dyn lies within 0.991–1.013 for α ≥ 1 (0.970–1.013 at α = 0.5) [DATA]. So the free-divider frequencies follow cot K = αK with one stiffness per cell.

![identity](../paper1_speedofsound/experiments/final/261004_p1_identity.png)

`paper1_speedofsound/experiments/final/261004_p1_identity.png/.pdf`
- Left: k_S^dyn against k_T + F²/(N_s kT) on log axes, with the diagonal.
- Right: ρ_I per cell (1 σ thick, 2 σ thin), with C's registered signature 2Δ_C.

### 2.4 γ_box = k_S/k_T (§ 1.5, no pass/fail)

##### gamma_box = k_S^dyn / k_T (no pass/fail; sec. 1.5)

| eta | cell | scan | H | L_0 | N_s | gamma_box | sigma | bulk 1 + Z^2/(Z + eta Z') | 1 + F^2/(N_s kT k_T) (A alone) |
|---|---|---|---|---|---|---|---|---|---|
| 0.10 | e0p10_H_H5_L39.25 | H | 5 | 39.25 | 25 | 2.0694 | 0.0147 | 2.00930 | 2.0248 |
| 0.10 | e0p10_H_H10_L39.25 | H | 10 | 39.25 | 50 | 2.0209 | 0.0125 | 2.00930 | 2.0097 |
| 0.10 | e0p10_H_H20_L39.25 | H | 20 | 39.25 | 100 | 2.0176 | 0.0126 | 2.00930 | 2.0111 |
| 0.10 | e0p10_H_H40_L39.25 | H | 40 | 39.25 | 200 | 2.0117 | 0.0119 | 2.00930 | 2.0103 |
| 0.10 | e0p10_L_H10_L19.625 | L | 10 | 19.625 | 25 | 2.0646 | 0.0123 | 2.00930 | 2.0231 |
| 0.10 | e0p10_L_H10_L78.5 | L | 10 | 78.5 | 100 | 2.0107 | 0.0135 | 2.00930 | 2.0091 |
| 0.10 | e0p10_aspect_H19.7917_L19.7917 | aspect | 19.7917 | 19.7917 | 50 | 2.0651 | 0.0123 | 2.00934 | 2.0315 |
| 0.10 | e0p10_aspect_H14_L28 | aspect | 14 | 28 | 50 | 2.0321 | 0.0124 | 2.00933 | 2.0152 |
| 0.10 | e0p10_aspect_H9.91667_L39.625 | aspect | 9.91667 | 39.625 | 50 | 2.0618 | 0.0134 | 2.00928 | 2.0316 |
| 0.10 | e0p10_aspect_H7_L56.0417 | aspect | 7 | 56.0417 | 50 | 2.0312 | 0.0131 | 2.00931 | 2.0167 |
| 0.39 | epi8_H_H5_L10 | H | 5 | 10 | 25 | 2.4275 | 0.0169 | 2.18760 | 2.3725 |
| 0.39 | epi8_H_H10_L10 | H | 10 | 10 | 50 | 2.4827 | 0.0180 | 2.18760 | 2.3982 |
| 0.39 | epi8_H_H20_L10 | H | 20 | 10 | 100 | 2.6272 | 0.0205 | 2.18760 | 2.4744 |
| 0.39 | epi8_H_H40_L10 | H | 40 | 10 | 200 | 3.0670 | 0.0279 | 2.18760 | 2.7154 |
| 0.39 | epi8_L_H10_L5 | L | 10 | 5 | 25 | 2.9038 | 0.0220 | 2.18760 | 2.6545 |
| 0.39 | epi8_L_H10_L20 | L | 10 | 20 | 100 | 2.3646 | 0.0181 | 2.18760 | 2.3126 |
| 0.39 | epi8_aspect_H7.08333_L14.125 | aspect | 7.08333 | 14.125 | 50 | 2.3486 | 0.0173 | 2.18736 | 2.3184 |
| 0.39 | epi8_aspect_H5_L20 | aspect | 5 | 20 | 50 | 2.3578 | 0.0177 | 2.18760 | 2.3175 |
| 0.39 | epi8_aspect_H3.54167_L28.2917 | aspect | 3.54167 | 28.2917 | 50 | 2.3828 | 0.0177 | 2.18668 | 2.3618 |

eta 0.10: gamma_box inverse-variance mean 2.0377 +- 0.0040 (chi2 31.2 / 9 dof), range 2.0107 .. 2.0694; bulk 2.00930; A-alone mean 2.0183

eta 0.39: gamma_box inverse-variance mean 2.4926 +- 0.0063 (chi2 1050.1 / 8 dof), range 2.3486 .. 3.0670; bulk 2.18760; A-alone mean 2.4361

- At η = 0.10, γ_box = 2.038 ± 0.004, against the bulk 2.009 [DATA]. It is highest in the small-N_s cells.
- At π/8 the registered γ_box (2.35–3.07) inherits the biased k_T of § 2.8. With the drift correction (post-hoc) it is 2.26–2.41, against the bulk 2.188 [DATA].

### 2.5 Damping (§ 1.4 (B), methods § 13)

##### Damping per (cell, M): Gamma^-1 = tau_r/2 [sigma-time] +- jackknife (methods sec. 13 model, ACF to 20 periods)

| eta | cell | alpha 0.5 | alpha 1 | alpha 2 | alpha 3 | alpha 5 | alpha 7.5 | alpha 10 | alpha 15 | alpha 20 |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.10 | e0p10_H_H5_L39.25 | 244 +- 20 | 422 +- 43 | 772 +- 48 | 1.3e+03 +- 95 | 1.79e+03 +- 2.3e+02 | 3.15e+03 +- 2.5e+02 | 3.86e+03 +- 4.4e+02 | 4.12e+03 +- 3.9e+02 | 6.53e+03 +- 9.9e+02 |
| 0.10 | e0p10_H_H10_L39.25 | 322 +- 36 | 508 +- 42 | 1.07e+03 +- 1e+02 | 1.75e+03 +- 85 | 2.38e+03 +- 1.6e+02 | 3.29e+03 +- 3.3e+02 | 4.8e+03 +- 3.9e+02 | 7.28e+03 +- 4.8e+02 | 9.2e+03 +- 7.7e+02 |
| 0.10 | e0p10_H_H20_L39.25 | 307 +- 18 | 554 +- 55 | 1.12e+03 +- 86 | 1.59e+03 +- 1.4e+02 | 2.67e+03 +- 2.7e+02 | 4.16e+03 +- 4.8e+02 | 5.59e+03 +- 6.2e+02 | 8.78e+03 +- 9.4e+02 | 1.27e+04 +- 1.8e+03 |
| 0.10 | e0p10_H_H40_L39.25 | 361 +- 30 | 646 +- 52 | 1.23e+03 +- 1.3e+02 | 1.69e+03 +- 1.4e+02 | 2.67e+03 +- 1.8e+02 | 4.15e+03 +- 4e+02 | 6.86e+03 +- 7.9e+02 | 8.76e+03 +- 1.1e+03 | 8.55e+03 +- 7.9e+02 |
| 0.10 | e0p10_L_H10_L19.625 | 79.6 +- 1.9 | 148 +- 9.3 | 273 +- 16 | 396 +- 29 | 679 +- 70 | 971 +- 36 | 1.15e+03 +- 86 | 1.88e+03 +- 1.1e+02 | 2.42e+03 +- 1.9e+02 |
| 0.10 | e0p10_L_H10_L78.5 | 1.06e+03 +- 1.3e+02 | 1.86e+03 +- 1.2e+02 | 4.23e+03 +- 5.4e+02 | 5.78e+03 +- 4.9e+02 | 9.87e+03 +- 9.2e+02 | 1.45e+04 +- 2.1e+03 | 1.71e+04 +- 3.3e+03 | 2.45e+04 +- 2.9e+03 | 3.32e+04 +- 3.3e+03 |
| 0.10 | e0p10_aspect_H19.7917_L19.7917 | 97.9 +- 6.1 | 158 +- 12 | 285 +- 15 | 449 +- 28 | 714 +- 73 | 1.15e+03 +- 75 | 1.46e+03 +- 1.1e+02 | 2.26e+03 +- 1.7e+02 | 2.39e+03 +- 2e+02 |
| 0.10 | e0p10_aspect_H14_L28 | 162 +- 8.9 | 279 +- 23 | 621 +- 45 | 787 +- 47 | 1.38e+03 +- 1.3e+02 | 1.99e+03 +- 1.5e+02 | 2.77e+03 +- 2.7e+02 | 3.9e+03 +- 3.3e+02 | 5.87e+03 +- 6.5e+02 |
| 0.10 | e0p10_aspect_H9.91667_L39.625 | 298 +- 13 | 555 +- 42 | 974 +- 89 | 1.59e+03 +- 1e+02 | 2.96e+03 +- 2.2e+02 | 4.17e+03 +- 3.1e+02 | 6.09e+03 +- 5.9e+02 | 6.49e+03 +- 5.8e+02 | 1.03e+04 +- 7.4e+02 |
| 0.10 | e0p10_aspect_H7_L56.0417 | 484 +- 27 | 959 +- 47 | 1.98e+03 +- 2.1e+02 | 2.94e+03 +- 3.2e+02 | 4.69e+03 +- 5.9e+02 | 7.48e+03 +- 7.1e+02 | 9.18e+03 +- 8e+02 | 1.29e+04 +- 9.1e+02 | 1.56e+04 +- 2.1e+03 |
| 0.39 | epi8_H_H5_L10 | 18.1 +- 0.53 | 30 +- 1.6 | 76.2 +- 6.1 | 107 +- 8.6 | 163 +- 11 | 240 +- 16 | 338 +- 21 | 451 +- 44 | 599 +- 40 |
| 0.39 | epi8_H_H10_L10 | 20.8 +- 1 | 41.4 +- 2.9 | 78.1 +- 5.4 | 114 +- 10 | 163 +- 16 | 225 +- 17 | 400 +- 38 | 433 +- 31 | 734 +- 51 |
| 0.39 | epi8_H_H20_L10 | 22.1 +- 1.4 | 39.2 +- 2.1 | 66 +- 7.4 | 103 +- 6.7 | 173 +- 16 | 336 +- 41 | 397 +- 37 | 517 +- 43 | 688 +- 32 |
| 0.39 | epi8_H_H40_L10 | 21.1 +- 0.8 | 39.5 +- 2.1 | 74.3 +- 5.9 | 108 +- 8.7 | 192 +- 14 | 240 +- 16 | 316 +- 20 | 468 +- 40 | 899 +- 89 |
| 0.39 | epi8_L_H10_L5 | 4.84 +- 0.22 | 8.17 +- 0.4 | 15.5 +- 0.93 | 23 +- 1.7 | 39.5 +- 2.7 | 54.7 +- 4.1 | 77.3 +- 4.7 | 113 +- 8.9 | 137 +- 11 |
| 0.39 | epi8_L_H10_L20 | 77.3 +- 5.5 | 174 +- 21 | 334 +- 42 | 431 +- 37 | 751 +- 89 | 904 +- 80 | 1.33e+03 +- 2.2e+02 | 2.19e+03 +- 3.6e+02 | 2.97e+03 +- 4.8e+02 |
| 0.39 | epi8_aspect_H7.08333_L14.125 | 43 +- 3.2 | 67.2 +- 5.2 | 148 +- 14 | 196 +- 17 | 334 +- 28 | 619 +- 45 | 746 +- 57 | 956 +- 89 | 1.28e+03 +- 85 |
| 0.39 | epi8_aspect_H5_L20 | 81.5 +- 7.2 | 150 +- 8.3 | 288 +- 12 | 403 +- 24 | 593 +- 31 | 1.07e+03 +- 97 | 1.47e+03 +- 1.2e+02 | 1.99e+03 +- 2.2e+02 | 2.55e+03 +- 3.9e+02 |
| 0.39 | epi8_aspect_H3.54167_L28.2917 | 113 +- 8 | 260 +- 16 | 532 +- 52 | 744 +- 54 | 1.19e+03 +- 88 | 1.81e+03 +- 1.1e+02 | 2.47e+03 +- 3.3e+02 | 3.11e+03 +- 3.5e+02 | 4.18e+03 +- 6.5e+02 |

ACF fits: 171 of 171 converged

### 2.6 Per-cell summary (X3 (a))

##### X3 (a) -- per-cell table

| eta | cell | H | L_0 | N_s | eta_true | c_s^B +- err (scaled) | k_S^dyn +- | k_T +- | F(L_0) | k_T + F^2/(N_s kT) | gamma_box = k_S/k_T | Gamma^-1 at alpha = 5 (M = 10 N_s) +- |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.10 | e0p10_H_H5_L39.25 | 5 | 39.25 | 25 | 0.100051 | 1.7868 +- 0.0035 | 0.054596 +- 8.1e-05 | 0.026383 +- 0.00018 | 0.82215 | 0.05342 | 2.069 +- 0.015 | 1788 +- 2.3e+02 |
| 0.10 | e0p10_H_H10_L39.25 | 10 | 39.25 | 50 | 0.100051 | 1.7616 +- 0.0033 | 0.10655 +- 0.00011 | 0.052725 +- 0.00032 | 1.6315 | 0.10596 | 2.021 +- 0.013 | 2381 +- 1.6e+02 |
| 0.10 | e0p10_H_H20_L39.25 | 20 | 39.25 | 100 | 0.100051 | 1.7485 +- 0.0083 | 0.21077 +- 0.00021 | 0.10446 +- 0.00065 | 3.25 | 0.21009 | 2.018 +- 0.013 | 2670 +- 2.7e+02 |
| 0.10 | e0p10_H_H40_L39.25 | 40 | 39.25 | 200 | 0.100051 | 1.7451 +- 0.0036 | 0.4188 +- 0.00045 | 0.20819 +- 0.0012 | 6.4858 | 0.41851 | 2.012 +- 0.012 | 2674 +- 1.8e+02 |
| 0.10 | e0p10_L_H10_L19.625 | 10 | 19.625 | 25 | 0.100051 | 1.7650 +- 0.0035 | 0.22558 +- 0.00037 | 0.10926 +- 0.00063 | 1.6718 | 0.22106 | 2.065 +- 0.012 | 678.9 +- 70 |
| 0.10 | e0p10_L_H10_L78.5 | 10 | 78.5 | 100 | 0.100051 | 1.7655 +- 0.0031 | 0.051781 +- 3.7e-05 | 0.025753 +- 0.00017 | 1.612 | 0.051739 | 2.011 +- 0.014 | 9866 +- 9.2e+02 |
| 0.10 | e0p10_aspect_H19.7917_L19.7917 | 19.7917 | 19.7917 | 50 | 0.100252 | 1.7572 +- 0.0039 | 0.43788 +- 0.00065 | 0.21204 +- 0.0012 | 3.3069 | 0.43076 | 2.065 +- 0.012 | 713.8 +- 73 |
| 0.10 | e0p10_aspect_H14_L28 | 14 | 28 | 50 | 0.100178 | 1.7627 +- 0.0044 | 0.2128 +- 0.00026 | 0.10472 +- 0.00063 | 2.3055 | 0.21102 | 2.032 +- 0.012 | 1379 +- 1.3e+02 |
| 0.10 | e0p10_aspect_H9.91667_L39.625 | 9.91667 | 39.625 | 50 | 0.099937 | 1.7628 +- 0.0035 | 0.10442 +- 9.2e-05 | 0.050643 +- 0.00033 | 1.6162 | 0.10289 | 2.062 +- 0.013 | 2963 +- 2.2e+02 |
| 0.10 | e0p10_aspect_H7_L56.0417 | 7 | 56.0417 | 50 | 0.100104 | 1.7772 +- 0.0065 | 0.051878 +- 4.9e-05 | 0.025541 +- 0.00016 | 1.1395 | 0.051509 | 2.031 +- 0.013 | 4693 +- 5.9e+02 |
| 0.39 | epi8_H_H5_L10 | 5 | 10 | 25 | 0.392699 | 3.9072 +- 0.0224 | 4.7723 +- 0.0057 | 1.9659 +- 0.013 | 8.213 | 4.664 | 2.428 +- 0.017 | 163.3 +- 11 |
| 0.39 | epi8_H_H10_L10 | 10 | 10 | 50 | 0.392699 | 3.8203 +- 0.0124 | 9.0069 +- 0.011 | 3.6278 +- 0.026 | 15.925 | 8.7 | 2.483 +- 0.018 | 163.3 +- 16 |
| 0.39 | epi8_H_H20_L10 | 20 | 10 | 100 | 0.392699 | 3.7581 +- 0.0091 | 17.504 +- 0.02 | 6.6623 +- 0.051 | 31.342 | 16.486 | 2.627 +- 0.020 | 173.2 +- 16 |
| 0.39 | epi8_H_H40_L10 | 40 | 10 | 200 | 0.392699 | 3.7166 +- 0.0172 | 34.562 +- 0.047 | 11.269 +- 0.1 | 62.179 | 30.6 | 3.067 +- 0.028 | 192.2 +- 14 |
| 0.39 | epi8_L_H10_L5 | 10 | 5 | 25 | 0.392699 | 3.6819 +- 0.0251 | 21.655 +- 0.051 | 7.4574 +- 0.054 | 17.563 | 19.796 | 2.904 +- 0.022 | 39.47 +- 2.7 |
| 0.39 | epi8_L_H10_L20 | 10 | 20 | 100 | 0.392699 | 3.8506 +- 0.0055 | 4.1273 +- 0.004 | 1.7455 +- 0.013 | 15.136 | 4.0366 | 2.365 +- 0.018 | 751 +- 89 |
| 0.39 | epi8_aspect_H7.08333_L14.125 | 7.08333 | 14.125 | 50 | 0.392495 | 3.8794 +- 0.0059 | 4.3869 +- 0.0047 | 1.8679 +- 0.014 | 11.097 | 4.3306 | 2.349 +- 0.017 | 333.8 +- 28 |
| 0.39 | epi8_aspect_H5_L20 | 5 | 20 | 50 | 0.392699 | 3.9828 +- 0.0054 | 2.2026 +- 0.0018 | 0.93416 +- 0.007 | 7.8447 | 2.1649 | 2.358 +- 0.018 | 592.5 +- 31 |
| 0.39 | epi8_aspect_H3.54167_L28.2917 | 3.54167 | 28.2917 | 50 | 0.391917 | 3.9992 +- 0.0045 | 1.0756 +- 0.0009 | 0.45143 +- 0.0033 | 5.5441 | 1.0662 | 2.383 +- 0.018 | 1186 +- 88 |

tables -> 261004_p1_confinement_cells.csv, 261004_p1_confinement_damping.csv

figure -> 261004_p1_confinement_shift.png/.pdf
figure -> 261004_p1_identity.png/.pdf

gates: inventory PASS, reduction PASS, determinism PASS

Files: `paper1_speedofsound/experiments/final/261004_p1_confinement_cells.csv` (every column above plus η_rec, δ, L_eff,true, χ²_red, Δ, σ_Δ, Δ_C, ρ_I) and `261004_p1_confinement_damping.csv` (per (cell, M): τ_r, Γ, Γ⁻¹, P_1 = B, τ_T, the fitted ω, Q, with jackknife errors).

### 2.7 POST-HOC diagnosis (not registered, no verdict): the "held" divider is released and returns towards the centre

**Finding [DATA].** The full pilot cell (π/8 anchor) has the divider's trajectory. In all 16 off-centre runs, the divider is released at t = 200 with mass 10⁹ (`--wall-hold-steps=12000 --wall-mass-factors=1000000000`) and drifts back towards the centre: by the end of the window it stands at about 0.90 of its nominal offset.
- Its window-mean offset is f = 0.9662 ± 0.0013 of nominal, the same at ±1 and ±2 dL and for every seed.
- The registration's drift estimate (§ 1.4, Table A) covered only the random thermal wander. It missed this deterministic return, which is driven by the net restoring force at an off-centre position.

**Why it biases k_T [DERIVATION].** Each compartment is closed, so the slow return compresses or expands each gas adiabatically. The time-averaged force is then F̄(j) = F_T(L₀ + x_j) − k_S (x̄_j − x_j), and the registered stencil, which uses nominal positions, returns k_T,meas = k_T − k_S(1 − f): too low.

**Measured in every cell, without traces [DATA, DERIVATION].** Energy conservation in each compartment gives N_s(T̄_L − 1) = −F̄_L(x̄_j − x_j), so f follows from the recorded temperatures.
- At the anchor this gives f = 0.9677 ± 0.0003, against 0.9662 ± 0.0013 from the trajectories.
- The temperatures shift exactly as predicted: at x_m2, T_L − 1 = −0.00270 against a predicted −0.00268.
- Across all 19 cells, 1 − f is proportional to k_S, as a harmonic return predicts.

**Correction.** The corrected estimator, which was not registered, normalises each force to T = 1 (F/T, exact for hard disks, F = T g(L)) and uses the measured spacing f·dL. Printed by `python3 hspist3/validation/paper1_confinement_heldwall_posthoc_261004.py`:

pilot traces (pi/8 anchor, 16 off-centre runs, W0_x_sigma over the window [200, end]): mean displacement / nominal = f = 0.9662 +- 0.0013 (min 0.9631, max 0.9689)
same pilot runs, f from the temperatures (energy balance): 0.9677 +- 0.0003  (per position: -2dL 0.9679, -1dL 0.9665, +1dL 0.9680, +2dL 0.9678)

##### Per cell: f from the temperatures, and the identity with the drift-corrected k_T (exploratory)

| eta | cell | N_s | f (T balance) | sigma_f | k_T registered | k_T,c (F/T, spacing f dL) | sigma | k_T,c/k_T - 1 [%] | static,c = k_T,c + F^2/(N_s kT) | k_S^dyn | rho_I,c [%] | sigma [%] | rho_I,c/sigma | registered rho_I [%] | 2 Delta_C [%] | gamma_box,c | bulk gamma |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.10 | e0p10_H_H5_L39.25 | 25 | 0.9998 | 0.0000 | 0.02638 | 0.02639 | 0.00018 | +0.04 | 0.05343 | 0.05460 | +2.133 | 0.367 | +5.81 | +2.155 | +2.537 | 2.0685 | 2.00930 |
| 0.10 | e0p10_H_H10_L39.25 | 50 | 0.9996 | 0.0000 | 0.05272 | 0.05277 | 0.00032 | +0.09 | 0.10601 | 0.10655 | +0.510 | 0.319 | +1.60 | +0.554 | +1.269 | 2.0191 | 2.00930 |
| 0.10 | e0p10_H_H20_L39.25 | 100 | 0.9991 | 0.0000 | 0.10446 | 0.10464 | 0.00065 | +0.17 | 0.21027 | 0.21077 | +0.238 | 0.323 | +0.74 | +0.324 | +0.634 | 2.0141 | 2.00930 |
| 0.10 | e0p10_H_H40_L39.25 | 200 | 0.9983 | 0.0000 | 0.20819 | 0.20891 | 0.00121 | +0.35 | 0.41924 | 0.41880 | -0.106 | 0.309 | -0.34 | +0.068 | +0.317 | 2.0047 | 2.00930 |
| 0.10 | e0p10_L_H10_L19.625 | 25 | 0.9991 | 0.0000 | 0.10926 | 0.10947 | 0.00063 | +0.18 | 0.22126 | 0.22558 | +1.917 | 0.322 | +5.96 | +2.006 | +2.537 | 2.0608 | 2.00930 |
| 0.10 | e0p10_L_H10_L78.5 | 100 | 0.9998 | 0.0000 | 0.02575 | 0.02576 | 0.00017 | +0.04 | 0.05175 | 0.05178 | +0.059 | 0.340 | +0.17 | +0.081 | +0.634 | 2.0098 | 2.00930 |
| 0.10 | e0p10_aspect_H19.7917_L19.7917 | 50 | 0.9982 | 0.0000 | 0.21204 | 0.21281 | 0.00123 | +0.36 | 0.43152 | 0.43788 | +1.452 | 0.317 | +4.58 | +1.627 | +1.268 | 2.0576 | 2.00934 |
| 0.10 | e0p10_aspect_H14_L28 | 50 | 0.9991 | 0.0000 | 0.10472 | 0.10490 | 0.00063 | +0.18 | 0.21121 | 0.21280 | +0.746 | 0.319 | +2.34 | +0.832 | +1.268 | 2.0285 | 2.00933 |
| 0.10 | e0p10_aspect_H9.91667_L39.625 | 50 | 0.9996 | 0.0000 | 0.05064 | 0.05069 | 0.00033 | +0.09 | 0.10293 | 0.10442 | +1.425 | 0.326 | +4.37 | +1.467 | +1.269 | 2.0600 | 2.00928 |
| 0.10 | e0p10_aspect_H7_L56.0417 | 50 | 0.9998 | 0.0000 | 0.02554 | 0.02555 | 0.00016 | +0.04 | 0.05152 | 0.05188 | +0.691 | 0.327 | +2.11 | +0.712 | +1.269 | 2.0303 | 2.00931 |
| 0.39 | epi8_H_H5_L10 | 25 | 0.9830 | 0.0001 | 1.96591 | 2.04600 | 0.01390 | +4.07 | 4.74414 | 4.77230 | +0.590 | 0.315 | +1.87 | +2.268 | +1.536 | 2.3325 | 2.18760 |
| 0.39 | epi8_H_H10_L10 | 50 | 0.9681 | 0.0002 | 3.62779 | 3.91532 | 0.02712 | +7.93 | 8.98752 | 9.00689 | +0.215 | 0.327 | +0.66 | +3.407 | +0.768 | 2.3004 | 2.18760 |
| 0.39 | epi8_H_H20_L10 | 100 | 0.9390 | 0.0003 | 6.66234 | 7.74197 | 0.05715 | +16.20 | 17.56514 | 17.50357 | -0.352 | 0.347 | -1.01 | +5.816 | +0.384 | 2.2609 | 2.18760 |
| 0.39 | epi8_H_H40_L10 | 200 | 0.8847 | 0.0005 | 11.26896 | 15.26485 | 0.12250 | +35.46 | 34.59584 | 34.56199 | -0.098 | 0.379 | -0.26 | +11.464 | +0.192 | 2.2642 | 2.18760 |
| 0.39 | epi8_L_H10_L5 | 25 | 0.9279 | 0.0004 | 7.45739 | 8.97555 | 0.05936 | +20.36 | 21.31400 | 21.65452 | +1.572 | 0.361 | +4.35 | +8.583 | +1.536 | 2.4126 | 2.18760 |
| 0.39 | epi8_L_H10_L20 | 100 | 0.9849 | 0.0001 | 1.74548 | 1.80702 | 0.01354 | +3.53 | 4.09813 | 4.12729 | +0.707 | 0.342 | +2.07 | +2.198 | +0.384 | 2.2840 | 2.18760 |
| 0.39 | epi8_aspect_H7.08333_L14.125 | 50 | 0.9842 | 0.0001 | 1.86791 | 1.93752 | 0.01393 | +3.73 | 4.40026 | 4.38689 | -0.305 | 0.335 | -0.91 | +1.282 | +0.768 | 2.2642 | 2.18736 |
| 0.39 | epi8_aspect_H5_L20 | 50 | 0.9921 | 0.0000 | 0.93416 | 0.95137 | 0.00703 | +1.84 | 2.18215 | 2.20258 | +0.928 | 0.330 | +2.81 | +1.709 | +0.768 | 2.3152 | 2.18760 |
| 0.39 | epi8_aspect_H3.54167_L28.2917 | 50 | 0.9962 | 0.0000 | 0.45143 | 0.45552 | 0.00336 | +0.91 | 1.07025 | 1.07564 | +0.501 | 0.323 | +1.55 | +0.881 | +0.769 | 2.3614 | 2.18668 |

eta 0.10: drift-corrected rho_I within 2 sigma in 4 of 10 cells; sum (rho/sigma)^2 = 122.6 (10 cells); inverse-variance mean rho_I,c = +0.878 +- 0.103 %

eta 0.39: drift-corrected rho_I within 2 sigma in 6 of 9 cells; sum (rho/sigma)^2 = 39.4 (9 cells); inverse-variance mean rho_I,c = +0.421 +- 0.113 %
eta 0.10: rho_I,c = a x 2 Delta_C: a = 0.75 +- 0.07, chi2 11.2 / 9 dof (rho_I,c = 0: chi2 122.6 / 10 dof)
eta 0.39: rho_I,c = a x 2 Delta_C: a = 0.58 +- 0.12, chi2 17.2 / 8 dof (rho_I,c = 0: chi2 39.4 / 9 dof)

**Reading [DATA; INFERENCE where marked].**
- **Size of the bias.** The registered k_T is low by 0.9–35 % at π/8 (35 % at H = 40) and by ≤ 0.4 % at η = 0.10.
- **π/8 after correction.** The residuals fall from +0.9 … +11.5 % to −0.35 … +1.57 %, with 6 of 9 cells within 2 σ. The remaining outliers are L₀ = 5 (+1.57 ± 0.36 %), H5/L20 (+0.93 ± 0.33 %) and L20 (+0.71 ± 0.34 %).
- **η = 0.10.** Nothing changes.
- **Pattern of what remains.** At both densities the residual follows C's registered 1/N_s signature, at 0.75 ± 0.07 (η = 0.10; χ² 11.2/9) and 0.58 ± 0.12 (π/8; χ² 17.2/8) of its fixed amplitude.
- **[INFERENCE] A tension.** Hypothesis C predicts the same mode shift in Δ. In the η = 0.10 L-scan it would separate N_s = 25 from N_s = 100 by about 0.7 % at that amplitude, but they agree (+1.19 and +1.21 %). The source of the residual is therefore OPEN. Candidates: C acting on the heavy masses only, or a static-side effect of order 1/N_s.
- **What this means for the verdict.** The registered verdict above stands. A corrected verdict would need an amendment, C4, which is a decision for the plan author: the k_T estimator with the measured f, applied once to the existing data.

![held divider, post-hoc](../paper1_speedofsound/experiments/final/261004_p1_identity_heldwall_posthoc.png)

`paper1_speedofsound/experiments/final/261004_p1_identity_heldwall_posthoc.png/.pdf`
- Left: the pilot's 16 divider trajectories.
- Middle: 1 − f against k_S for every cell, with the pilot trajectories as a star.
- Right: ρ_I, registered (open markers) and drift-corrected (filled), with 2Δ_C.

### 2.8 Verdicts at a glance

| test (registration) | rule | result |
|---|---|---|
| confinement, η = 0.10 (§ 1.5) | one free amplitude each; excluded if p < 0.01 | **not separated**: A (p 0.14) and B (p 0.88) survive, C excluded; Δχ²(A − B) = +9.1 |
| confinement, π/8 (§ 1.5) | same | **none survives** (A, B, C, C fixed); exploratory two-term forms fail too → **not resolved** |
| identity, length-free (C1) | within 2 σ at every cell | **FAIL**, 4/19 [post-hoc drift-corrected: 10/19; the π/8 failures are mostly the released-divider bias] |
| γ_box (§ 1.5) | no pass/fail | η 0.10: 2.038 ± 0.004 (bulk 2.009); π/8 registered 2.35–3.07, drift-corrected 2.26–2.41 (bulk 2.188) |
| gates | inventory, reduction (§ 1.10), determinism (§ 1.12), health | all PASS; one trajectory excluded by the health rule (effect ≤ 0.008 σ on c_s) |


---

## 3. PRE-REGISTRATION "A-fixed": the identity with a divider that is actually held (2026-10-04, before any A-fixed run)

Written and committed before any A-fixed trajectory exists. Decisions taken by the plan author on 2026-10-04:
- the registered verdicts of § 2 stand;
- C4 (the temperature-based drift correction of § 2.7) is a documented POST-HOC analysis, not the identity result;
- the identity is re-measured, by the design below.

### 3.1 Why

The registered identity test (C1) failed: 4 of 19 cells were within 2 σ (§ 2.3). Afterwards, § 2.7 found the main reason in the data: method A held its divider for the 200 σ-time equilibration only, then released it with mass 10⁹. In the off-centre runs the divider returned towards the centre while the force was being recorded, which biased k_T low, by 0.9–35 % at π/8.

A-fixed measures the same static side again with a divider held for the entire record.

### 3.2 Code facts (Task Y1) [SOURCE, quoted]

- **The hold is on by default:** `bool wall_hold_enabled = true;` (`00ALLINONE.c:294`).
- **While held, the core is given mass 0 and velocity 0:** `if (hold_active) { div_mass[w] = 0.0; div_vx[w] = 0.0; }` (`00ALLINONE.c:16890–16892`).
- **A mass-0 divider has no spring** (`divider_has_spring` requires mass > 0, `edmd.c:253`). Its position update is therefore `S->prm.divider_x[d] += S->prm.divider_vx[d] * dt;` (`edmd.c:779`), which adds exactly 0: the divider does not move.
- **Every collision with it is resolved by the infinite-mass branch and logged:** `if (M <= 0.0){ double v1 = 2.0 * u2 - u1; … edmd_log_event(S, kb, u2, u1, v1, dE); …` (`edmd.c:1174–1183`). Each one appears in the event log as `D0` with dp = v₁ − u₁.
  - **So yes: the forces F_L and F_R are recorded during the hold**, from the event log, exactly as after the release.
- **The run length:** `recorded_steps = 0;` at the release (`00ALLINONE.c:17055`), `if (!wall_is_released) continue;` (`:17069`) and `recorded_steps++;` (`:17205`) with `target_steps = num_steps` (`:16592`). So `--steps` counts released steps only.
  - The trace (KE, hence T_L and T_R) is written after the release only, and its time is `simulation_time - wall_release_time` (`:17086`).
- **[DATA] Pre-check on existing data**, printed by `python3 hspist3/cluster/afix_pilot_check_261004.py`. In all 20 method-A pilot runs, the 200 σ-time hold has about 2,500 divider collisions, each with u_wall = 0 and Σ dE = 0 exactly:

##### Pre-check on existing data: the hold phase (t < 200) of the 20 method-A pilot runs

| position | seed | D0 events, t < 200 | max abs u_wall | sum dE |
|---|---|---|---|---|
| x_m2 | 9700 | 2577 | 0 | 0 |
| x_m2 | 9701 | 2613 | 0 | 0 |
| x_m2 | 9702 | 2589 | 0 | 0 |
| x_m2 | 9703 | 2485 | 0 | 0 |
| x_m1 | 9700 | 2642 | 0 | 0 |
| x_m1 | 9701 | 2527 | 0 | 0 |
| x_m1 | 9702 | 2546 | 0 | 0 |
| x_m1 | 9703 | 2544 | 0 | 0 |
| x_0 | 9700 | 2570 | 0 | 0 |
| x_0 | 9701 | 2550 | 0 | 0 |
| x_0 | 9702 | 2573 | 0 | 0 |
| x_0 | 9703 | 2502 | 0 | 0 |
| x_p1 | 9700 | 2562 | 0 | 0 |
| x_p1 | 9701 | 2601 | 0 | 0 |
| x_p1 | 9702 | 2547 | 0 | 0 |
| x_p1 | 9703 | 2546 | 0 | 0 |
| x_p2 | 9700 | 2538 | 0 | 0 |
| x_p2 | 9701 | 2502 | 0 | 0 |
| x_p2 | 9702 | 2543 | 0 | 0 |
| x_p2 | 9703 | 2511 | 0 | 0 |

pre-check: the divider is immovable during the hold in all 20 runs

**Consequence: no code change.** A-fixed runs on the same binary as method B, `279282b target koa`. The identity figure therefore compares two methods of one build generation, and no byte-identity evidence across builds is needed.

### 3.3 Design

- **Cells, positions, seeds and stencil** are those of method A, by construction. `hspist3/cluster/gen_afix_sbatch_261004.py` turns every line of `tasks_A_<cell>.txt` one-to-one into an `AF` line, and checks that the (position, seed) sets are equal; they are, in all 19 cells (table below).
- **Flags** are those of method A (`conf_worker.sh` mode A), except `--wall-hold-steps=312000` (200 σ-time equilibration + 5000 σ-time record, dt = 1/60 σ-time) and `--steps=1200`. The 20 σ-time released tail exists only so that the trace records the temperatures (mode `AF`).
- **Window:** [200, 5200) σ-time, i.e. the equilibration and the released tail are both excluded (`reduce_AF.py`).
- **Output** goes to `experiments_energy_transfer/paper1_confinement_Afix_261004/` (new). The method-A directories are never written.
- **Seed reuse [INFERENCE].** The first 200 σ-time of every A-fixed run is identical to its method-A run; gate G1(a) checks this. The records then differ (held against released) and decorrelate within a few collision times. P1 treats the two as independent; any residual positive correlation would make P1 conservative.

Printed by `python3 hspist3/cluster/gen_afix_sbatch_261004.py` (verbatim):

##### A-fixed arrays (261012 sec. 3): tasks from the method-A task files, hold 312000 steps, tail 1200, every 600

| group | task | cell | N_s | trajectories | same (position, seed) set as method A | predicted wall (h) | basis | core-h (wall x 16) |
|---|---|---|---|---|---|---|---|---|
| Afix_0.10 | 1 | e0p10_H_H5_L39.25 | 25 | 710 | yes | 0.03 | Round 1 sacct (measured) | 0.5 |
| Afix_0.10 | 2 | e0p10_H_H10_L39.25 | 50 | 720 | yes | 0.05 | Round 1 sacct (measured) | 0.8 |
| Afix_0.10 | 3 | e0p10_H_H20_L39.25 | 100 | 675 | yes | 0.11 | Round 1 sacct (measured) | 1.8 |
| Afix_0.10 | 4 | e0p10_H_H40_L39.25 | 200 | 720 | yes | 0.93 | H20 per wave x 2^2.99 | 14.9 |
| Afix_0.10 | 5 | e0p10_L_H10_L19.625 | 25 | 750 | yes | 0.03 | Round 1 sacct (measured) | 0.6 |
| Afix_0.10 | 6 | e0p10_L_H10_L78.5 | 100 | 675 | yes | 0.11 | Round 1 sacct (measured) | 1.7 |
| Afix_0.10 | 7 | e0p10_aspect_H19.7917_L19.7917 | 50 | 730 | yes | 0.06 | Round 1 sacct (measured) | 1.0 |
| Afix_0.10 | 8 | e0p10_aspect_H14_L28 | 50 | 685 | yes | 0.06 | Round 1 sacct (measured) | 0.9 |
| Afix_0.10 | 9 | e0p10_aspect_H9.91667_L39.625 | 50 | 680 | yes | 0.05 | Round 1 sacct (measured) | 0.8 |
| Afix_0.10 | 10 | e0p10_aspect_H7_L56.0417 | 50 | 685 | yes | 0.05 | Round 1 sacct (measured) | 0.8 |
| Afix_0.10 | -- | sbatch default --time = 3 x the largest non-H40 cell = 0:30:00; H40 override 3:00:00 | | | | | | |
| Afix_0.39 | 1 | epi8_H_H5_L10 | 25 | 130 | yes | 0.05 | pilot per wave x (N_s/50)^2.99 | 0.8 |
| Afix_0.39 | 2 | epi8_H_H10_L10 | 50 | 115 | yes | 0.05 | pilot per wave x (N_s/50)^2.99 | 0.7 |
| Afix_0.39 | 3 | epi8_H_H20_L10 | 100 | 130 | yes | 0.41 | pilot per wave x (N_s/50)^2.99 | 6.5 |
| Afix_0.39 | 4 | epi8_H_H40_L10 | 200 | 255 | yes | 5.78 | pilot per wave x (N_s/50)^2.99 | 92.4 |
| Afix_0.39 | 5 | epi8_L_H10_L5 | 25 | 130 | yes | 0.05 | pilot per wave x (N_s/50)^2.99 | 0.8 |
| Afix_0.39 | 6 | epi8_L_H10_L20 | 100 | 130 | yes | 0.41 | pilot per wave x (N_s/50)^2.99 | 6.5 |
| Afix_0.39 | 7 | epi8_aspect_H7.08333_L14.125 | 50 | 130 | yes | 0.05 | pilot per wave x (N_s/50)^2.99 | 0.8 |
| Afix_0.39 | 8 | epi8_aspect_H5_L20 | 50 | 115 | yes | 0.05 | pilot per wave x (N_s/50)^2.99 | 0.7 |
| Afix_0.39 | 9 | epi8_aspect_H3.54167_L28.2917 | 50 | 105 | yes | 0.04 | pilot per wave x (N_s/50)^2.99 | 0.6 |
| Afix_0.39 | -- | sbatch default --time = 3 x the largest non-H40 cell = 1:15:00; H40 override 17:30:00 | | | | | | |

A-fixed pilot: 20 trajectories (the method-A pilot's tasks, held), sandbox, 1:00:00

total (upper bound, wall x 16 cores): Afix_0.10 23.8 core-h, Afix_0.39 110.0 core-h, together 133.9 core-h; p* = 2.99

submission lines (runsheet step 9; at most 64 cores: 2 x 16 + 2 x 16):

    sbatch --array=1-1 cluster/confinement_20261013/conf_Afix_pilot.sbatch      (first, sandbox; gate G1)
    sbatch --array=4 --time=3:00:00 cluster/confinement_20261013/conf_Afix_0.10.sbatch
    sbatch --array=1,2,3,5,6,7,8,9,10%2 cluster/confinement_20261013/conf_Afix_0.10.sbatch
    sbatch --array=4 --time=17:30:00 cluster/confinement_20261013/conf_Afix_0.39.sbatch
    sbatch --array=1,2,3,5,6,7,8,9%2 cluster/confinement_20261013/conf_Afix_0.39.sbatch

files written: tasks_AF_*.txt (19 cells + pilot), cells_Afix_{0.10,0.39,pilot}.tsv, conf_Afix_{0.10,0.39,pilot}.sbatch, fetch_afix.sh

**Cost note [INFERENCE].** The η = 0.10 cells take their wall times from the Round 1 sacct, except H40. H40 and every π/8 cell use the Round-1 measured-scaling rule (p* = 2.99). The π/8 H40 cell (5.8 h, 92 core-h) dominates the total and is an extrapolation from N_s = 50. The Round-2 sacct of `conf-A_0.39` task 4 would replace it with a measurement. `--time` is 3× the prediction, so an over-estimate costs only queue priority.

### 3.4 Estimator (registered now)

- **Per seed:** F_L/T_L and F_R/T_R, the normalisation of C4. Here T ≡ 1 by construction, since each compartment is closed and the divider does no work.
- **Points:** F̃(j) is the mean of (F_L/T_L)(run j) and (F_R/T_R)(run −j).
- **k_T** = −[F̃(−2) − 8F̃(−1) + 8F̃(+1) − F̃(+2)]/(12 dL), with the nominal spacing dL (f = 1). Its σ comes from the seed standard errors.
- **Static side:** k_T + F(L₀)²/(N_s kT), with F(L₀) = ½(F_L + F_R) at x = 0 and kT the mean temperature of the x = 0 seeds, as in § 2.3.
- **Dynamic side:** k_S^dyn is method B's, unchanged (C1, § 2.3).
- **Residual:** ρ_I = (k_S^dyn − static)/k_S^dyn, with σ(ρ_I)² = (σ_kS/k_S)² + (σ_kT/k_S)².
- **Drift check, per cell.** f from the recorded temperatures by the C4 formula (`paper1_confinement_heldwall_posthoc_261004.drift`) must give |1 − f| < 0.002.
  - A cell that fails is flagged, reported, and left out of P1–P3. More than two flagged cells → stop (design failure).
  - Every seed must also have `u_wall_max` = 0 and `W_div` = 0.
- **The analysis script** will be `hspist3/validation/paper1_confinement_afix_<date>.py`, written before the data are opened and applied once.

### 3.5 Gates

- **G1 (before the arrays): the A-fixed pilot.** It runs the 20 tasks of the method-A pilot, held. `afix_pilot_check_261004.py` checks, per run:
  - (a) the event-log lines before t = 200 are identical to the method-A pilot's;
  - (b) u_wall = 0 for every D0 event before the release;
  - (c) Σ dE = 0;
  - (d) the record is complete;
  - (e) health 0;
  - (f) |1 − f| < 0.002.
  
  If (a) fails, stop: the two flags would be changing the trajectory before 200 σ-time, which must be explained before launch.
- **G2 (after): inventory,** as in § 2.1: every seed, the recorded geometry, the build, health 0.
- **G3 (after, per cell):** the drift check of § 3.4.

### 3.6 Predictions, and what each outcome means (registered now)

**P1: consistency with C4.**
- **Test:** per cell, ρ_I(A-fixed) − ρ_I(C4) must lie within 2 σ, with σ² = (σ_kT,A-fixed² + σ_kT,C4²)/k_S^dyn². k_S^dyn is common to both and cancels. The C4 values are the § 2.7 table (`261004_p1_identity_heldwall_posthoc`).
- **Rule:** P1 holds if every cell is within 2 σ. The χ² over the cells is reported as well.
- **If it holds:** the C4 correction is validated, and the A-fixed values become the paper's identity figure.
- **If it fails:** C4 is not an adequate correction of a moving divider. The A-fixed values supersede it, and the pattern of the differences is reported.

**P2: the identity, by the C1 rule.** |ρ_I(A-fixed)| ≤ 2 σ at every cell.
- **If it holds:** the dynamic and the static stiffness agree in every box; the identity holds in the confined system.
- **If it fails:** there is a real difference between dynamic and static stiffness, characterised by P3.

**P3: the 1/N_s residual.**
- **Fit:** per density, a weighted one-parameter fit ρ_I = c/N_s.
- **Report:** r = c/A_C, where A_C = 2 N_s Δ_C, the amplitude hypothesis C predicts: (q+1)(q+2)/(8qZ) = 0.634 at η = 0.10 and 0.384 at π/8 [DERIVATION, Table P: Δ_C = +0.634 % and +0.384 % at N_s = 50]. Also report r's σ, the χ² of the fit and the χ² of ρ_I = 0.
- **Outcomes, declared now:**
  - |c| < 2 σ_c at both densities: there is no 1/N_s residual. The C4-corrected residual (r = 0.75 ± 0.07 and 0.58 ± 0.12, § 2.7) was an artefact of correcting a moving divider.
  - c > 2 σ_c, with r within 2 σ of the C4 values: the residual is physics; the dynamic stiffness exceeds the static one in proportion to 1/N_s. The paper reports it, and its tension with the flat L-scan of Δ (§ 2.9) is stated as OPEN.
  - r within 2 σ of 1: hypothesis C's mechanism (thermal-amplitude anharmonicity) at its predicted size.
  - c < −2 σ_c: the static stiffness exceeds the dynamic one. This is new, and is reported as such.

### 3.7 What would change this registration

Only gate G1. If any of (a)–(f) fails, the design is revised, and amended in writing, before any array runs. Nothing is tuned after array data.
