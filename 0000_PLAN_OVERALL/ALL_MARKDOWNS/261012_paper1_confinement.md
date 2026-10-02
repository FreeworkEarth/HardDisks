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
