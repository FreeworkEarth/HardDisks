# Efficiency map, geometry C — PRE-REGISTRATION

2026-10-10. **Not launched.** — *2026-10-12: go given with amendments A1–A3 (§ 1.8), committed before launch.* Written ledger-first, as agreed: the exact energy accounting comes
before any observable is named, because three observables this week were named first and each one
measured a reversible first-order term instead of the quantity intended. Every number in the tables
below was printed by the script that computed it.

This is Paper 2's core figure and it has not been run.

---

## 1. PRE-REGISTRATION

### 1.1 The exact ledger, written before any observable

Geometry C: one gas of **N = 100** against a spring-loaded divider. `--l0=54.75` (box 109.5),
divider at 30.5 with thickness 1.0, gas from 31.0 to the right piston at 109.5, so
**L_0 = 78.4998**, **eta_0 = 0.10005098**, **c_s = 1.744338**, **Z_ac = N m c_s/L_0 = 2.2221**.
Piston on the right, travel d = 7.96.

At every instant, with L(t) the gas length from the recorded `SegEtas` and x(t) the recorded
divider position:

> **W_in(t) = E_qs,gas(L(t)) + X_gas(t) + KE_div(t) + E_spring(t)**
>
> E_qs,gas(L) = N (T_ad(L) - 1), T_ad from the KR isentrope d ln T = -Z d ln L from eta_0
> E_spring(t) = (1/2) k (x(t) - x_eq)^2 - (1/2) k (x(0) - x_eq)^2
> KE_div(t) = (1/2) M_s v(t)^2 from the recorded `W0_v`
> **X_gas(t) = the remainder** — the non-isentropic energy in the gas, acoustic + dissipated heat

**This is a closed identity and it is the first thing the analysis checks**, per seed, at every
sample, exactly as Level 4b did (which closed to 1.0e-6 kT). *(Wording amended in § 1.8, A3: the check is W_in = ΔKE_gas + KE_div + ΔE_spring; E_qs + X is the KR split of ΔKE_gas.)*

> **At settle, KE_div = 0 and X_gas is thermalised, so**
> **epsilon = E_spring/W_in exactly**, with no modelling.

That is the efficiency, and it is a ledger entry rather than a fitted quantity.

### 1.2 The spring rest length is NOT free, and Level 3 only ever ran one of the three

Mechanical equilibrium at t = 0 requires k (x_eq - 30.5) = P(L_0) h, with **P h = 1.57489**:

| k | **x_eq required** | note |
|---|---|---|
| 0.25 | **36.7996** | |
| 0.50 | **33.6498** | **this is Level 3's 33.65 — it was the k = 0.5 equilibrium** |
| 1.00 | **32.0749** | |

**Running k = 0.25 or k = 1.0 with Level 3's `--spring-eq=33.65` would start the system with a
large preload** (a force imbalance of 0.79 and 1.58 respectively, i.e. 50 % and 100 % of the gas
force) and the measured epsilon would be meaningless. Each k gets its own x_eq. x_eq is a spring
rest length, not a wall position, so the 1/24-sigma grid does not constrain it.

### 1.3 The reversible reference, one equation per cell

The reversible final state is where the spring force balances the KR isentropic pressure force,
with the geometry tying the two together: L = 109 - d - x, so

> **k (x_eq - (109 - d - L)) = P(L) h**, one equation, solved for L_f.

| k | x_eq | L_f | divider moves | E_qs(gas) | E_spring | W_rev | **epsilon_rev** |
|---|---|---|---|---|---|---|---|
| 0.25 | 36.7996 | 72.0387 | +1.4987 | 11.3217 | 2.6411 | 13.9628 | **0.1892** |
| 0.50 | 33.6498 | 71.3803 | +0.8403 | 12.6201 | 1.5000 | 14.1201 | **0.1062** |
| 1.00 | 32.0749 | 70.9880 | +0.4480 | 13.4083 | 0.8060 | 14.2142 | **0.0567** |

**epsilon_rev depends on k and not on M_s or u** — that is the whole point of a reversible
reference, and it is the first thing the data must reproduce as u -> 0.

### 1.4 Predictions, fixed before the data

1. **epsilon(u, k, M_s) -> epsilon_rev(k) as u -> 0**, i.e. 0.1892, 0.1062, 0.0567, with **no M_s
   dependence in the limit**. Any residual M_s dependence at u = 0.01 is divider inertia and is
   reported as such.
2. **1 - epsilon = (E_qs + X_gas)/W_in**, with X_gas the dissipation: **first order in u**, of order
   **Z_ac u d = 17.69 u**, plus the mode energy left in the spring-divider oscillation.
   So epsilon should fall roughly linearly in u with slope set by Z_ac d/W_rev. *(Restricted in § 1.8, A2: only for d/u < 2L_0/c_s = 90, i.e. u > 0.088.)*
3. **The no-push control** (same cell, piston never released) must give W_in = 0 and
   epsilon undefined; it measures the thermal drift of x and hence the floor on E_spring.
4. The **ledger identity of 1.1 closes per seed at every sample**, to float rounding.

### 1.5 Grid, seeds, cost

**u in {0.01, 0.02, 0.05, 0.1, 0.2, 0.5} x k in {0.25, 0.5, 1.0} x M_s in {50, 200}** = 36 cells,
**8 seeds** each, plus a **no-push control per (k, M_s)** = 6 controls x 8 seeds. 336 runs.

Run length **d/u + 5 spring periods**, period = 2 pi sqrt(M_s/k):

| k | M_s | period | 5 periods | run at u = 0.01 | steps | run at u = 0.5 | steps |
|---|---|---|---|---|---|---|---|
| 0.25 | 50 | 88.9 | 444 | 1240 | 75 000 | 460 | 28 000 |
| 0.25 | 200 | 177.7 | 889 | 1685 | 102 000 | 904 | 55 000 |
| 0.50 | 50 | 62.8 | 314 | 1110 | 67 000 | 330 | 20 000 |
| 0.50 | 200 | 125.7 | 628 | 1424 | 86 000 | 644 | 39 000 |
| 1.00 | 50 | 44.4 | 222 | 1018 | 62 000 | 238 | 15 000 |
| 1.00 | 200 | 88.9 | 444 | 1240 | 75 000 | 460 | 28 000 |

> **Total ~16 M steps, ~0.07 core-hours.** The whole map costs less than one cell of the R-collapse.
> *(Superseded by § 1.8, A1: run length max(d/u + 5P, d/u + 3 tau_r), ~388 M steps, ~1.3 core-h.)*

`--trace-every` set so dt <= period/40 in every cell, and the reduction **fails closed** on any
missing column (Time, KE_gas_total, W0_x_sigma, W0_v, PistonWork, PistonR_x_sigma, SegCounts,
SegEtas, SpringE).

### 1.6 The piston gap, as established this week

Geometry C parks the piston **0.25 sigma outside** the gas (verified: `piston_left_x = XW1 - 6.0f`,
`piston_right_x = XW2 + 6.0f`, and `compute_segment_bounds` uses
`right_plane = fminf(XW2, piston_right_x)`). The compression therefore **starts at t = 0.25/u, not
at t = 0**, and `piston_stop_t_rel` exceeds d/u by that amount.

> **The compression start is taken from the first nonzero `PistonWork` per seed**, and **d from
> `piston_target_sigma`**, never from d/u or from the stop time. Confirmed in Level 3's own
> geometry: `XW2 - piston_target = 7.9600` exactly, and PistonWork is zero for all t < 5.000 = 0.25/u.

---

### 1.7 Added 2026-10-10: binary, the v1 gap rule, and the pre-registered Level 3 reproduction

**Binary.** Piston v2 was not adopted (branch RED, `05215ea`). The map runs on **v1 +
`-ffp-contract=off` + `--version`**, and every summary CSV carries `build_git / build_target /
build_cflags`. The **0.25-sigma piston gap therefore stays**: the piston parks at XW2 + 0.25 and
does zero work until t = 0.25/u; the compression equals the recorded travel d = 7.96 exactly
(verified on Level 3's own geometry: `XW2 − piston_target = 7.9600`, PistonWork = 0 for
t < 0.25/u). So, as established: **tau_push = piston_stop_t_rel − 0.25/u**, the compression start is
read from the first nonzero PistonWork per seed, and d from `piston_target_sigma`. No number in
sections 1.1–1.6 changes.

**Pre-registered consistency line.** The cell **(k = 0.5, u = 0.05, M_s = 200)** is Level 3's
geometry C cell (`level3_v6_20260924/k0.5_M200_u0.05`, same L0 = 54.75, spring_eq 33.65 = 33.6498,
travel 7.96). Its settled displacement ratio must reproduce Level 3's closing number
**s̄ = 0.8705 ± 0.0049 within 2 sigma**, computed with Level 3's own estimator on the new run. *(s̄ is a settled displacement in σ, not a ratio; § 1.8.)*
Pass → the map is on the same footing as Level 3. Fail → the map does not proceed to a figure
until the discrepancy is understood; it is reported, not absorbed.

**Never mixed:** Level 3's data (contraction-on binary) and the map (contraction-off) do not appear
in one figure; the reproduction line is a comparison of two numbers, each with its own error.

---

### 1.8 Amendments A1–A3 (2026-10-12), committed before launch — the go was given with these

**What changed, one line each.**
- **A1** adds two efficiency definitions and extends every run to at least $3\tau_r$ after the push, so that $\varepsilon_{\rm settled}$ exists.
- **A2** restricts prediction 2 to $\tau_{\rm push} < 2L_0/c_s$.
- **A3** makes the exact ledger the per-seed check, and labels $E_{\rm qs} + X$ as a model split of $\Delta KE_{\rm gas}$.

#### A1 — two efficiencies, both pre-registered

With $s(t) = -(x(t) - x(0))$, the divider displacement away from the gas (Level 3 v6's sign), the spring energy change for a displacement $s$ is exactly

$$\Delta E_{\rm spring}(s) = F_s\,s + \tfrac12 k s^2,\qquad F_s = k\,(x_{\rm eq} - 30.5).$$

$F_s = 1.57489$ at $k = 0.25$ and $1.0$; at $k = 0.5$ it is $1.5750$, because the runner keeps Level 3's $x_{\rm eq} = 33.65$. The two efficiencies are

$$\varepsilon_{\rm mean} = \frac{\Delta E_{\rm spring}\big(\bar s_{\rm corr}[\text{last } 3T_w]\big)}{\langle W_{\rm in}\rangle},\qquad
\varepsilon_{\rm settled} = \frac{\Delta E_{\rm spring}\big(\bar s_{\rm corr}[\text{last } \tau_r]\big)}{\langle W_{\rm in}\rangle},$$

with $\bar s_{\rm corr} = \bar s - \bar s_{\rm ctrl}$, the control mean taken over the same absolute window of the matching no-push control. $T_w = 2\pi\sqrt{M_s/(k + k_{\rm gas})}$ as in Level 3 v6, and $\tau_r$ is from table 1 below. The ratio is of means, with a jackknife over seeds (numerator and denominator jointly) plus the control's standard error. $W_{\rm in}$ is PistonWork at the end of the record.

**$\varepsilon_{\rm mean}$ is Level 3's estimator in form.** The § 1.7 reproduction line runs Level 3 v6's code path verbatim on k0.5_M200_u0.05 with ctrl_k0.5_M200: the window runs from $d/u + 3T_w$ to the end, with Level 3's control window. **$\bar s$ there is a settled displacement in σ**: Level 3 v6 prints "s̄ = 0.8705 ± 0.0049 sigma", pooled over its six valid cells. That is what "ratio" in § 1.7 refers to.

**Why $\varepsilon_{\rm settled}$ uses the window-mean position, not the mean of the instantaneous $E_{\rm spring}$.** At settle the spring carries a thermal $\tfrac12 kT_f$ on top of its static value: 0.56 kT (table 4), i.e. about +0.04 in $\varepsilon$ at $W_{\rm in} \approx 14$, which is larger than the effects being measured. The position-based form excludes it by construction. $W_{\rm in}$ is unaffected, because the thermal shares come out of $KE_{\rm gas}$.

**The $KE_{\rm div}$ gate cannot be implemented as written.** $\langle KE_{\rm div}\rangle$ has the equipartition floor $kT_f/2 = 0.557$–$0.567$, which is 21–70× the threshold $0.01\,\Delta E_{\rm spring}$ (table 4). It is implemented instead in two parts:

1. **By construction.** Every run continues at least $3\tau_r$ after the push. $\tau_r$ comes from Mansour, and that is an upper bound here: the measured damping exceeded Mansour in all seven 261006 cells (260913 REPORT, 1b). For a sudden step of the equilibrium position, the coherent energy averaged over the last-$\tau_r$ window is then $0.0085$–$0.0101\,\Delta E_{\rm spring}$ (table 4).
2. **Measured.** $E_{\rm coh} = k_{\rm eff}V_c$, with
$$V_c = \frac{n\,{\rm Var}_t\,\bar x - \langle {\rm Var}_t\,x_i\rangle}{n-1},\qquad k_{\rm eff} = \frac{kT_f}{\langle {\rm Var}_t\,x_i\rangle - V_c},$$
from the excess variance of the seed-mean trajectory, with a jackknife σ. Per cell the verdict is **PASS** if $E_{\rm coh} + 2\sigma < 0.01\,\Delta E_{\rm spring}$, **FAIL** if $E_{\rm coh} - 2\sigma > 0.01\,\Delta E_{\rm spring}$, and **UNRESOLVED** otherwise. The 8-seed floor $kT_f/16 \approx 0.07$ exceeds every threshold, so UNRESOLVED is the expected verdict, and then part 1 carries the gate. A FAIL is reported, and that cell's $\varepsilon_{\rm settled}$ is flagged and not used.

**Run length.** $T = 0.25/u + d/u + \max(5P,\ 3\tau_r)$. The 0.25/u gap is now included (the old rule was $d/u + 5P$). Steps are 60 per σ-time after the hold; a 120 000-step run records $t = 0 \to 1999.07$ and starts at release. See tables 1–3.

#### A2 — prediction 2, restricted

Linear-in-$u$ dissipation, $X_{\rm gas} \approx Z_{\rm ac}\,u\,d$, holds only while $\tau_{\rm push} = d/u < 2L_0/c_s = 90.0$, i.e. for $u > 0.088$. Slower pushes leave the gas time to re-equilibrate acoustically, and the excess falls toward the quasi-static value (cf. 4b's cell (200, 0.05)). **Expected shape: $\varepsilon \approx \varepsilon_{\rm rev}$ for $u \lesssim 0.1$, then falling roughly linearly in $u$.** On this grid, $u = 0.01, 0.02, 0.05$ sit on the quasi-static side, 0.1 at the boundary, and 0.2 and 0.5 on the linear side. The slope $d\varepsilon/du$ is fitted on $u \in \{0.1, 0.2, 0.5\}$ and reported beside $-\varepsilon_{\rm rev}Z_{\rm ac}d/W_{\rm rev}$, which is an order-of-magnitude reference only: it ignores the extra spring loading by the heated gas.

#### A3 — the check and the model, separated

The per-seed check, at every sample, uses the recorded PistonWork, KE_gas_total, W0_v and W0_x_sigma:

$$W_{\rm in}(t) = \Delta KE_{\rm gas}(t) + KE_{\rm div}(t) + \Delta E_{\rm spring}(t).$$

The maximum residual is reported; it is limited to about 1e-5 kT by $x$ being printed to 6 decimals. **The recorded `SpringE` is 0 in the first trace row** (from row 1 on it equals $\tfrac12 k(x - x_{\rm eq})^2$), so $E_{\rm spring}$ is computed from $x$. SpringE is compared with it on rows ≥ 1, and the difference is reported.

$E_{\rm qs}(L) + X$ is the KR decomposition of $\Delta KE_{\rm gas}$, with $L(t)$ from `SegEtas`. It is a model split, labelled as such, and not part of the check.

**Added (INFERENCE).** At release, the held divider's degree of freedom thermalises and takes

$$E_{\rm dof} = \tfrac12 kT_f\Big(1 + \frac{k}{k + k_S}\Big) \approx kT_f$$

out of the gas KE: $kT/2$ kinetic, plus the spring's share of the potential $kT/2$ (the gas-compression share stays in $KE_{\rm gas}$). So at settle $X = \Delta KE_{\rm gas} - E_{\rm qs} \approx -kT_f$ even with no dissipation, and **the dissipation is reported as $X + E_{\rm dof}$**. The no-push control tests this directly: its late $\Delta KE_{\rm gas}$ must equal $-E_{\rm dof}$ (prediction 3 table). This was found while smoke-testing the pipeline on throwaway seeds; no number from them is used.

**Analysis.** `hspist3/validation/paper2_effmap_analysis_20261012.py`, committed with this section, reads only `red_*.csv` and fails closed on a missing column or a trace that does not start at release.

**Launch.** After this commit, the pictures come first (`watch.sh effmap shot <k>`, one per k), then `level5_effmap_20261010.sh`.

**Tables printed by `python3 hspist3/validation/paper2_effmap_amend_20261012.py`** (verbatim):

#### 0. Reproduction of the pre-registered sec. 1.3 (reversible reference)

eta_0 = 0.10005098, P h = F(L_0) = 1.57489

| k | x_eq | L_f | divider moves | E_qs(gas) | E_spring | W_rev | epsilon_rev | T_ad(L_f) |
|---|---|---|---|---|---|---|---|---|
| 0.25 | 36.7996 | 72.0387 | +1.4987 | 11.3217 | 2.6411 | 13.9628 | 0.1892 | 1.11322 |
| 0.5 | 33.6498 | 71.3803 | +0.8403 | 12.6201 | 1.5000 | 14.1201 | 0.1062 | 1.12620 |
| 1.0 | 32.0749 | 70.9880 | +0.4480 | 13.4083 | 0.8060 | 14.2142 | 0.0567 | 1.13408 |

#### 1. tau_r per (k, M_s): Mansour friction, ONE gas column, at the settled state

gamma_1 = L_y (eta_s + zeta)/L_f (Eq. 17, one term = half the two-gas rate); M_hat = M_s + N m/3 (Eq. 18);
Enskog at eta_f = eta_0 L_0/L_f and T_f = T_ad(L_f) (KR isentrope); tau_r = 2 M_hat/gamma_1 (amplitude time).

| k | M_s | eta_f | T_f | eta_s + zeta | gamma_1 | M_hat | Gamma_E = gamma_1/M_hat | **tau_r** | 3 tau_r | 5 P (spring) |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.25 | 50 | 0.10902 | 1.1132 | 0.3567 | 0.04951 | 83.33 | 5.941e-04 | **3366** | 10099 | 444 |
| 0.25 | 200 | 0.10902 | 1.1132 | 0.3567 | 0.04951 | 233.33 | 2.122e-04 | **9426** | 28277 | 889 |
| 0.5 | 50 | 0.11003 | 1.1262 | 0.3596 | 0.05038 | 83.33 | 6.046e-04 | **3308** | 9924 | 314 |
| 0.5 | 200 | 0.11003 | 1.1262 | 0.3596 | 0.05038 | 233.33 | 2.159e-04 | **9262** | 27786 | 628 |
| 1.0 | 50 | 0.11064 | 1.1341 | 0.3615 | 0.05092 | 83.33 | 6.110e-04 | **3273** | 9820 | 222 |
| 1.0 | 200 | 0.11064 | 1.1341 | 0.3615 | 0.05092 | 233.33 | 2.182e-04 | **9165** | 27495 | 444 |

Mansour is an UPPER bound on tau_r here: in all 7 cells of the 261006 ladder the measured damping
exceeded Mansour's friction (260913 REPORT, 1b), so 3 tau_r(Mansour) is the conservative length.

#### 2. Run lengths and steps: old (d/u + 5P) vs amended (0.25/u + d/u + max(5P, 3 tau_r))

| cell | P | T old | steps old | T new | steps new | trace-every | rows/trace |
|---|---|---|---|---|---|---|---|
| k0.25_M50_u0.01 | 88.9 | 1240 | 75000 | 10920 | 656000 | 133 | 4932 |
| k0.25_M50_u0.02 | 88.9 | 842 | 51000 | 10510 | 631000 | 133 | 4744 |
| k0.25_M50_u0.05 | 88.9 | 603 | 37000 | 10263 | 616000 | 133 | 4631 |
| k0.25_M50_u0.1 | 88.9 | 524 | 32000 | 10181 | 611000 | 133 | 4593 |
| k0.25_M50_u0.2 | 88.9 | 484 | 30000 | 10140 | 609000 | 133 | 4578 |
| k0.25_M50_u0.5 | 88.9 | 460 | 28000 | 10115 | 607000 | 133 | 4563 |
| ctrl_k0.25_M50 | 88.9 | 1240 | 75000 | 10920 | 656000 | 133 | 4932 |
| k0.25_M200_u0.01 | 177.7 | 1685 | 102000 | 29098 | 1746000 | 266 | 6563 |
| k0.25_M200_u0.02 | 177.7 | 1287 | 78000 | 28688 | 1722000 | 266 | 6473 |
| k0.25_M200_u0.05 | 177.7 | 1048 | 63000 | 28442 | 1707000 | 266 | 6417 |
| k0.25_M200_u0.1 | 177.7 | 968 | 59000 | 28359 | 1702000 | 266 | 6398 |
| k0.25_M200_u0.2 | 177.7 | 928 | 56000 | 28318 | 1700000 | 266 | 6390 |
| k0.25_M200_u0.5 | 177.7 | 904 | 55000 | 28294 | 1698000 | 266 | 6383 |
| ctrl_k0.25_M200 | 177.7 | 1685 | 102000 | 29098 | 1746000 | 266 | 6563 |
| k0.5_M50_u0.01 | 62.8 | 1110 | 67000 | 10745 | 645000 | 94 | 6861 |
| k0.5_M50_u0.02 | 62.8 | 712 | 43000 | 10334 | 621000 | 94 | 6606 |
| k0.5_M50_u0.05 | 62.8 | 473 | 29000 | 10088 | 606000 | 94 | 6446 |
| k0.5_M50_u0.1 | 62.8 | 394 | 24000 | 10006 | 601000 | 94 | 6393 |
| k0.5_M50_u0.2 | 62.8 | 354 | 22000 | 9965 | 598000 | 94 | 6361 |
| k0.5_M50_u0.5 | 62.8 | 330 | 20000 | 9940 | 597000 | 94 | 6351 |
| ctrl_k0.5_M50 | 62.8 | 1110 | 67000 | 10745 | 645000 | 94 | 6861 |
| k0.5_M200_u0.01 | 125.7 | 1424 | 86000 | 28607 | 1717000 | 188 | 9132 |
| k0.5_M200_u0.02 | 125.7 | 1026 | 62000 | 28197 | 1692000 | 188 | 9000 |
| k0.5_M200_u0.05 | 125.7 | 788 | 48000 | 27951 | 1678000 | 188 | 8925 |
| k0.5_M200_u0.1 | 125.7 | 708 | 43000 | 27868 | 1673000 | 188 | 8898 |
| k0.5_M200_u0.2 | 125.7 | 668 | 41000 | 27827 | 1670000 | 188 | 8882 |
| k0.5_M200_u0.5 | 125.7 | 644 | 39000 | 27803 | 1669000 | 188 | 8877 |
| ctrl_k0.5_M200 | 125.7 | 1424 | 86000 | 28607 | 1717000 | 188 | 9132 |
| k1.0_M50_u0.01 | 44.4 | 1018 | 62000 | 10641 | 639000 | 66 | 9681 |
| k1.0_M50_u0.02 | 44.4 | 620 | 38000 | 10230 | 614000 | 66 | 9303 |
| k1.0_M50_u0.05 | 44.4 | 381 | 23000 | 9984 | 600000 | 66 | 9090 |
| k1.0_M50_u0.1 | 44.4 | 302 | 19000 | 9902 | 595000 | 66 | 9015 |
| k1.0_M50_u0.2 | 44.4 | 262 | 16000 | 9861 | 592000 | 66 | 8969 |
| k1.0_M50_u0.5 | 44.4 | 238 | 15000 | 9836 | 591000 | 66 | 8954 |
| ctrl_k1.0_M50 | 44.4 | 1018 | 62000 | 10641 | 639000 | 66 | 9681 |
| k1.0_M200_u0.01 | 88.9 | 1240 | 75000 | 28316 | 1699000 | 133 | 12774 |
| k1.0_M200_u0.02 | 88.9 | 842 | 51000 | 27906 | 1675000 | 133 | 12593 |
| k1.0_M200_u0.05 | 88.9 | 603 | 37000 | 27659 | 1660000 | 133 | 12481 |
| k1.0_M200_u0.1 | 88.9 | 524 | 32000 | 27577 | 1655000 | 133 | 12443 |
| k1.0_M200_u0.2 | 88.9 | 484 | 30000 | 27536 | 1653000 | 133 | 12428 |
| k1.0_M200_u0.5 | 88.9 | 460 | 28000 | 27511 | 1651000 | 133 | 12413 |
| ctrl_k1.0_M200 | 88.9 | 1240 | 75000 | 28316 | 1699000 | 133 | 12774 |

#### 3. Cost

timed here: one throwaway 240 000-step run of 00ALLINONE in this geometry -> **12.02 us per step** (incl. the 12000-step hold)
total steps (8 seeds, 36 cells + 6 controls): old 16.7 M -> **new 388.0 M**
CPU: **1.30 core-h**; wall at 9 jobs ~ **9 min** (longest single run 0.3 min)

#### 4. Floors that decide how the KE_div gate can be implemented

Equipartition gives the divider a thermal <KE_div> = kT_f/2 in its one degree of freedom, and the
instantaneous E_spring a thermal kT_f/2 on top of its static value. The coherent-energy floor of an
8-seed average is kT_f/(2 x 8). Sudden-limit bound on the coherent energy averaged over the last-tau_r
window of a run >= 3 tau_r long: (1 + k_S/k) e^-4 (1 - e^-2)/2 x Delta E_spring, k_S = N m c_s^2/L_f^2.

| k | Delta E_spring (rev) | 0.01 Delta E_spring | kT_f/2 (thermal KE_div) | ratio | kT_f/16 (8-seed floor) | k_S/k | sudden-limit E_coh/Delta E_spring |
|---|---|---|---|---|---|---|---|
| 0.25 | 2.6411 | 0.0264 | 0.5566 | 21x | 0.0696 | 0.2716 | 0.0101 |
| 0.5 | 1.5000 | 0.0150 | 0.5631 | 38x | 0.0704 | 0.1406 | 0.0090 |
| 1.0 | 0.8060 | 0.0081 | 0.5670 | 70x | 0.0709 | 0.0718 | 0.0085 |

A2 boundary: 2 L_0/c_s(eta_0) = 90.0 sigma-time -> tau_push = d/u below it for u > 0.0884

---

## 2. Results

*Empty. Not launched. Chris and the plan author read section 1 first.*
