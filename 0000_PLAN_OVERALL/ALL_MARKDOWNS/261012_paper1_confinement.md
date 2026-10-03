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
