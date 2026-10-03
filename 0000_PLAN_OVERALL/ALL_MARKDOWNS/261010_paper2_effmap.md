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

**Run.** Launched 2026-10-01 17:16:22 (clock) on the post-flag binary (`05215ea`, release, `-ffp-contract=off`; "-dirty" refers to unrelated files, the core sources are identical to 05215ea), after § 1.8 was committed (b788e82). **336/336 runs, 0 failed, 0 health lines, 0 aborts**, 10.5 min wall, 306 MB. Pictures per k: `261012_effmap_k{0.25,0.5,1.0}_paper.png`.

### 2.1 Pre-registered analysis, run once

**Printed by `python3 hspist3/validation/paper2_effmap_analysis_20261012.py`** (committed before launch), verbatim:

#### Ledger (A3): max |W_in - (Delta KE_gas + KE_div + Delta E_spring)| per cell, all samples, all seeds

**max over all cells = 1.35e-05 kT**; recorded SpringE vs k(x - x_eq)^2/2 on rows >= 1: max 1.59e-05 kT

#### epsilon per cell (A1): ratio of means, control-corrected; errors jackknife + control SE

| k | M_s | u | seeds | <W_in> | s_corr (3 T_w) | **epsilon_mean** | s_corr (tau_r) | **epsilon_settled** | E_coh ± σ | 0.01 ΔE_spring | gate | epsilon_rev |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.25 | 50 | 0.01 | 8 | 14.1501 | 1.4728 | **0.1831 ± 0.0101** | 1.5503 | **0.1938 ± 0.0015** | -0.0030 ± 0.0373 | 0.0274 | UNRESOLVED | 0.1892 |
| 0.25 | 50 | 0.02 | 8 | 14.3209 | 1.6277 | **0.2021 ± 0.0056** | 1.5666 | **0.1937 ± 0.0018** | 0.0094 ± 0.0346 | 0.0277 | UNRESOLVED | 0.1892 |
| 0.25 | 50 | 0.05 | 8 | 14.4435 | 1.5666 | **0.1921 ± 0.0124** | 1.5659 | **0.1920 ± 0.0021** | 0.0095 ± 0.0667 | 0.0277 | UNRESOLVED | 0.1892 |
| 0.25 | 50 | 0.1 | 8 | 15.1256 | 1.5650 | **0.1832 ± 0.0090** | 1.5980 | **0.1875 ± 0.0018** | -0.0282 ± 0.0633 | 0.0284 | UNRESOLVED | 0.1892 |
| 0.25 | 50 | 0.2 | 8 | 17.1641 | 1.7370 | **0.1814 ± 0.0095** | 1.7222 | **0.1796 ± 0.0028** | -0.0047 ± 0.0950 | 0.0308 | UNRESOLVED | 0.1892 |
| 0.25 | 50 | 0.5 | 8 | 25.4512 | 2.1014 | **0.1517 ± 0.0052** | 2.1668 | **0.1571 ± 0.0047** | -0.0042 ± 0.0403 | 0.0400 | UNRESOLVED | 0.1892 |
| 0.25 | 200 | 0.01 | 8 | 14.2585 | 1.5633 | **0.1941 ± 0.0103** | 1.5665 | **0.1945 ± 0.0019** | -0.0330 ± 0.1213 | 0.0277 | UNRESOLVED | 0.1892 |
| 0.25 | 200 | 0.02 | 8 | 14.2703 | 1.5505 | **0.1922 ± 0.0103** | 1.5450 | **0.1914 ± 0.0020** | -0.0786 ± 0.0770 | 0.0273 | UNRESOLVED | 0.1892 |
| 0.25 | 200 | 0.05 | 8 | 14.2929 | 1.6402 | **0.2043 ± 0.0119** | 1.5705 | **0.1946 ± 0.0019** | -0.0691 ± 0.0653 | 0.0278 | UNRESOLVED | 0.1892 |
| 0.25 | 200 | 0.1 | 8 | 15.5385 | 1.5275 | **0.1736 ± 0.0124** | 1.6287 | **0.1864 ± 0.0022** | 0.0102 ± 0.0701 | 0.0290 | UNRESOLVED | 0.1892 |
| 0.25 | 200 | 0.2 | 8 | 17.2524 | 1.6778 | **0.1735 ± 0.0127** | 1.7162 | **0.1780 ± 0.0030** | -0.0382 ± 0.1629 | 0.0307 | UNRESOLVED | 0.1892 |
| 0.25 | 200 | 0.5 | 8 | 25.4191 | 2.2999 | **0.1685 ± 0.0089** | 2.1771 | **0.1582 ± 0.0049** | 0.1906 ± 0.0908 | 0.0402 | UNRESOLVED | 0.1892 |
| 0.5 | 50 | 0.01 | 8 | 14.3414 | 1.0228 | **0.1306 ± 0.0087** | 0.8637 | **0.1079 ± 0.0006** | -0.0161 ± 0.0394 | 0.0155 | UNRESOLVED | 0.1062 |
| 0.5 | 50 | 0.02 | 8 | 14.2787 | 0.9188 | **0.1161 ± 0.0089** | 0.8535 | **0.1069 ± 0.0009** | 0.0385 ± 0.0618 | 0.0153 | UNRESOLVED | 0.1062 |
| 0.5 | 50 | 0.05 | 8 | 14.4228 | 0.8389 | **0.1038 ± 0.0077** | 0.8683 | **0.1079 ± 0.0008** | -0.0089 ± 0.0336 | 0.0156 | UNRESOLVED | 0.1062 |
| 0.5 | 50 | 0.1 | 8 | 14.8496 | 1.0231 | **0.1261 ± 0.0095** | 0.8865 | **0.1073 ± 0.0009** | -0.0310 ± 0.0368 | 0.0159 | UNRESOLVED | 0.1062 |
| 0.5 | 50 | 0.2 | 8 | 17.4588 | 0.9068 | **0.0936 ± 0.0062** | 0.9613 | **0.1000 ± 0.0015** | -0.0362 ± 0.0447 | 0.0175 | UNRESOLVED | 0.1062 |
| 0.5 | 50 | 0.5 | 8 | 25.4245 | 1.2287 | **0.0910 ± 0.0054** | 1.2120 | **0.0895 ± 0.0022** | 0.0037 ± 0.0336 | 0.0228 | UNRESOLVED | 0.1062 |
| 0.5 | 200 | 0.01 | 8 | 14.4762 | 0.8882 | **0.1103 ± 0.0066** | 0.8573 | **0.1060 ± 0.0009** | 0.0185 ± 0.1567 | 0.0153 | UNRESOLVED | 0.1062 |
| 0.5 | 200 | 0.02 | 8 | 14.3214 | 0.8808 | **0.1104 ± 0.0087** | 0.8558 | **0.1069 ± 0.0010** | -0.0531 ± 0.0845 | 0.0153 | UNRESOLVED | 0.1062 |
| 0.5 | 200 | 0.05 | 8 | 14.4494 | 0.8324 | **0.1027 ± 0.0087** | 0.8717 | **0.1082 ± 0.0008** | 0.0965 ± 0.1129 | 0.0156 | UNRESOLVED | 0.1062 |
| 0.5 | 200 | 0.1 | 8 | 14.7210 | 0.9903 | **0.1226 ± 0.0082** | 0.8891 | **0.1085 ± 0.0014** | 0.0117 ± 0.1238 | 0.0160 | UNRESOLVED | 0.1062 |
| 0.5 | 200 | 0.2 | 8 | 17.0360 | 0.8555 | **0.0898 ± 0.0089** | 0.9321 | **0.0989 ± 0.0014** | -0.0762 ± 0.0492 | 0.0169 | UNRESOLVED | 0.1062 |
| 0.5 | 200 | 0.5 | 8 | 25.4213 | 1.0993 | **0.0800 ± 0.0052** | 1.2151 | **0.0898 ± 0.0026** | 0.0748 ± 0.0911 | 0.0228 | UNRESOLVED | 0.1062 |
| 1.0 | 50 | 0.01 | 8 | 14.4137 | 0.4860 | **0.0613 ± 0.0071** | 0.4608 | **0.0577 ± 0.0004** | 0.0055 ± 0.0334 | 0.0083 | UNRESOLVED | 0.0567 |
| 1.0 | 50 | 0.02 | 8 | 14.5058 | 0.5103 | **0.0644 ± 0.0088** | 0.4661 | **0.0581 ± 0.0005** | 0.0243 ± 0.0601 | 0.0084 | UNRESOLVED | 0.0567 |
| 1.0 | 50 | 0.05 | 8 | 14.6093 | 0.3587 | **0.0431 ± 0.0074** | 0.4607 | **0.0569 ± 0.0006** | -0.0058 ± 0.0287 | 0.0083 | UNRESOLVED | 0.0567 |
| 1.0 | 50 | 0.1 | 8 | 15.1579 | 0.5880 | **0.0725 ± 0.0076** | 0.4720 | **0.0564 ± 0.0007** | 0.0045 ± 0.0452 | 0.0085 | UNRESOLVED | 0.0567 |
| 1.0 | 50 | 0.2 | 8 | 18.2446 | 0.6371 | **0.0661 ± 0.0087** | 0.5239 | **0.0527 ± 0.0010** | -0.0165 ± 0.0333 | 0.0096 | UNRESOLVED | 0.0567 |
| 1.0 | 50 | 0.5 | 8 | 25.4237 | 0.7248 | **0.0552 ± 0.0045** | 0.6424 | **0.0479 ± 0.0014** | 0.0212 ± 0.0336 | 0.0122 | UNRESOLVED | 0.0567 |
| 1.0 | 200 | 0.01 | 8 | 14.4324 | 0.4818 | **0.0606 ± 0.0045** | 0.4641 | **0.0581 ± 0.0003** | 0.0965 ± 0.1145 | 0.0084 | UNRESOLVED | 0.0567 |
| 1.0 | 200 | 0.02 | 8 | 14.5782 | 0.4184 | **0.0512 ± 0.0086** | 0.4673 | **0.0580 ± 0.0003** | 0.0062 ± 0.0848 | 0.0085 | UNRESOLVED | 0.0567 |
| 1.0 | 200 | 0.05 | 8 | 14.4441 | 0.4863 | **0.0612 ± 0.0066** | 0.4630 | **0.0579 ± 0.0004** | 0.0131 ± 0.0476 | 0.0084 | UNRESOLVED | 0.0567 |
| 1.0 | 200 | 0.1 | 8 | 15.0841 | 0.4351 | **0.0517 ± 0.0065** | 0.4784 | **0.0575 ± 0.0006** | -0.0488 ± 0.0537 | 0.0087 | UNRESOLVED | 0.0567 |
| 1.0 | 200 | 0.2 | 8 | 17.1097 | 0.5666 | **0.0615 ± 0.0061** | 0.5106 | **0.0546 ± 0.0005** | -0.0223 ± 0.0356 | 0.0093 | UNRESOLVED | 0.0567 |
| 1.0 | 200 | 0.5 | 8 | 25.4456 | 0.6719 | **0.0505 ± 0.0040** | 0.6535 | **0.0488 ± 0.0013** | -0.0228 ± 0.0283 | 0.0124 | UNRESOLVED | 0.0567 |

#### KR decomposition of Delta KE_gas at settle (A3; a model split, not the check)

| k | M_s | u | L (SegEtas, last tau_r) | T_f = KE/N | Delta KE_gas | E_qs(L) | X = Delta KE_gas - E_qs | E_dof (INFERENCE) | X + E_dof | Z_ac u d |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.25 | 50 | 0.01 | 72.1641 | 1.10422 | 10.4220 | 11.0777 | -0.6557 | 0.9874 | +0.3317 | 0.1769 |
| 0.25 | 50 | 0.02 | 72.1743 | 1.10324 | 10.3238 | 11.0580 | -0.7342 | 0.9867 | +0.2525 | 0.3538 |
| 0.25 | 50 | 0.05 | 72.1793 | 1.10436 | 10.4364 | 11.0482 | -0.6118 | 0.9876 | +0.3758 | 0.8844 |
| 0.25 | 50 | 0.1 | 72.2032 | 1.11148 | 11.1475 | 11.0020 | +0.1456 | 0.9934 | +1.1390 | 1.7688 |
| 0.25 | 50 | 0.2 | 72.3296 | 1.12846 | 12.8460 | 10.7574 | +2.0886 | 1.0076 | +3.0962 | 3.5376 |
| 0.25 | 50 | 0.5 | 72.7712 | 1.20340 | 20.3400 | 9.9116 | +10.4284 | 1.0695 | +11.4979 | 8.8440 |
| 0.25 | 200 | 0.01 | 72.1776 | 1.10507 | 10.5073 | 11.0516 | -0.5443 | 0.9882 | +0.4438 | 0.1769 |
| 0.25 | 200 | 0.02 | 72.1629 | 1.10323 | 10.3232 | 11.0802 | -0.7570 | 0.9866 | +0.2296 | 0.3538 |
| 0.25 | 200 | 0.05 | 72.1759 | 1.10595 | 10.5953 | 11.0549 | -0.4597 | 0.9889 | +0.5292 | 0.8844 |
| 0.25 | 200 | 0.1 | 72.2443 | 1.11726 | 11.7262 | 10.9224 | +0.8038 | 0.9983 | +1.8021 | 1.7688 |
| 0.25 | 200 | 0.2 | 72.3342 | 1.13064 | 13.0638 | 10.7485 | +2.3152 | 1.0094 | +3.3246 | 3.5376 |
| 0.25 | 200 | 0.5 | 72.7828 | 1.20347 | 20.3469 | 9.8897 | +10.4573 | 1.0695 | +11.5268 | 8.8440 |
| 0.5 | 50 | 0.01 | 71.4425 | 1.11541 | 11.5414 | 12.4962 | -0.9548 | 1.0474 | +0.0926 | 0.1769 |
| 0.5 | 50 | 0.02 | 71.4341 | 1.11517 | 11.5167 | 12.5129 | -0.9963 | 1.0471 | +0.0509 | 0.3538 |
| 0.5 | 50 | 0.05 | 71.4536 | 1.11672 | 11.6722 | 12.4742 | -0.8020 | 1.0486 | +0.2466 | 0.8844 |
| 0.5 | 50 | 0.1 | 71.4657 | 1.12197 | 12.1974 | 12.4500 | -0.2526 | 1.0532 | +0.8007 | 1.7688 |
| 0.5 | 50 | 0.2 | 71.5521 | 1.14556 | 14.5556 | 12.2784 | +2.2772 | 1.0743 | +3.3515 | 3.5376 |
| 0.5 | 50 | 0.5 | 71.7924 | 1.22046 | 22.0455 | 11.8040 | +10.2415 | 1.1408 | +11.3823 | 8.8440 |
| 0.5 | 200 | 0.01 | 71.4436 | 1.11593 | 11.5928 | 12.4940 | -0.9013 | 1.0478 | +0.1466 | 0.1769 |
| 0.5 | 200 | 0.02 | 71.4423 | 1.11534 | 11.5342 | 12.4967 | -0.9626 | 1.0473 | +0.0847 | 0.3538 |
| 0.5 | 200 | 0.05 | 71.4583 | 1.11957 | 11.9566 | 12.4648 | -0.5082 | 1.0511 | +0.5429 | 0.8844 |
| 0.5 | 200 | 0.1 | 71.4732 | 1.12275 | 12.2750 | 12.4353 | -0.1602 | 1.0539 | +0.8937 | 1.7688 |
| 0.5 | 200 | 0.2 | 71.5198 | 1.13794 | 13.7935 | 12.3427 | +1.4509 | 1.0675 | +2.5183 | 3.5376 |
| 0.5 | 200 | 0.5 | 71.7955 | 1.22253 | 22.2530 | 11.7978 | +10.4552 | 1.1426 | +11.5978 | 8.8440 |
| 1.0 | 50 | 0.01 | 71.0232 | 1.12514 | 12.5140 | 13.3372 | -0.8232 | 1.0878 | +0.2646 | 0.1769 |
| 1.0 | 50 | 0.02 | 71.0257 | 1.12533 | 12.5335 | 13.3321 | -0.7987 | 1.0880 | +0.2893 | 0.3538 |
| 1.0 | 50 | 0.05 | 71.0255 | 1.12619 | 12.6190 | 13.3326 | -0.7136 | 1.0888 | +0.3752 | 0.8844 |
| 1.0 | 50 | 0.1 | 71.0356 | 1.13116 | 13.1161 | 13.3122 | -0.1961 | 1.0934 | +0.8974 | 1.7688 |
| 1.0 | 50 | 0.2 | 71.0849 | 1.16146 | 16.1457 | 13.2126 | +2.9330 | 1.1218 | +4.0549 | 3.5376 |
| 1.0 | 50 | 0.5 | 71.1977 | 1.22976 | 22.9761 | 12.9856 | +9.9904 | 1.1857 | +11.1761 | 8.8440 |
| 1.0 | 200 | 0.01 | 71.0230 | 1.12345 | 12.3450 | 13.3376 | -0.9926 | 1.0862 | +0.0936 | 0.1769 |
| 1.0 | 200 | 0.02 | 71.0252 | 1.12546 | 12.5463 | 13.3332 | -0.7869 | 1.0881 | +0.3012 | 0.3538 |
| 1.0 | 200 | 0.05 | 71.0205 | 1.12493 | 12.4926 | 13.3427 | -0.8501 | 1.0876 | +0.2374 | 0.8844 |
| 1.0 | 200 | 0.1 | 71.0369 | 1.13210 | 13.2102 | 13.3095 | -0.0993 | 1.0943 | +0.9950 | 1.7688 |
| 1.0 | 200 | 0.2 | 71.0685 | 1.15025 | 15.0248 | 13.2458 | +1.7790 | 1.1113 | +2.8904 | 3.5376 |
| 1.0 | 200 | 0.5 | 71.2056 | 1.23129 | 23.1288 | 12.9698 | +10.1590 | 1.1871 | +11.3461 | 8.8440 |

#### Prediction 1: epsilon -> epsilon_rev as u -> 0 (u = 0.01 and 0.02 cells)

| k | M_s | u | epsilon_mean | epsilon_rev | (eps - eps_rev)/sigma | epsilon_settled | (eps_s - eps_rev)/sigma |
|---|---|---|---|---|---|---|---|
| 0.25 | 50 | 0.01 | 0.1831 | 0.1892 | -0.6 | 0.1938 | +3.2 |
| 0.25 | 50 | 0.02 | 0.2021 | 0.1892 | +2.3 | 0.1937 | +2.6 |
| 0.25 | 200 | 0.01 | 0.1941 | 0.1892 | +0.5 | 0.1945 | +2.9 |
| 0.25 | 200 | 0.02 | 0.1922 | 0.1892 | +0.3 | 0.1914 | +1.1 |
| 0.5 | 50 | 0.01 | 0.1306 | 0.1062 | +2.8 | 0.1079 | +2.7 |
| 0.5 | 50 | 0.02 | 0.1161 | 0.1062 | +1.1 | 0.1069 | +0.7 |
| 0.5 | 200 | 0.01 | 0.1103 | 0.1062 | +0.6 | 0.1060 | -0.3 |
| 0.5 | 200 | 0.02 | 0.1104 | 0.1062 | +0.5 | 0.1069 | +0.7 |
| 1.0 | 50 | 0.01 | 0.0613 | 0.0567 | +0.6 | 0.0577 | +2.3 |
| 1.0 | 50 | 0.02 | 0.0644 | 0.0567 | +0.9 | 0.0581 | +3.0 |
| 1.0 | 200 | 0.01 | 0.0606 | 0.0567 | +0.9 | 0.0581 | +5.0 |
| 1.0 | 200 | 0.02 | 0.0512 | 0.0567 | -0.6 | 0.0580 | +4.2 |

#### Prediction 2 (A2): flat for u <~ 0.09, then roughly linear; slope on u in {0.1, 0.2, 0.5}

| k | M_s | slope d eps_mean/du ± σ | -eps_rev Z_ac d / W_rev (order of magnitude) | mean eps_mean, u <= 0.05 |
|---|---|---|---|---|
| 0.25 | 50 | -0.0847 ± 0.0233 | -0.2396 | 0.1924 |
| 0.25 | 200 | -0.0139 ± 0.0350 | -0.2396 | 0.1968 |
| 0.5 | 50 | -0.0489 ± 0.0224 | -0.1331 | 0.1168 |
| 0.5 | 200 | -0.0863 ± 0.0221 | -0.1331 | 0.1078 |
| 1.0 | 50 | -0.0414 ± 0.0201 | -0.0706 | 0.0563 |
| 1.0 | 200 | -0.0144 ± 0.0169 | -0.0706 | 0.0577 |

#### Prediction 3: no-push controls (W_in = 0; thermal drift of x; floor on E_spring)

| k | M_s | seeds | max |PistonWork| | <s> whole run | ΔE_spring(<s>) | Delta KE_gas, last half | -E_dof predicted |
|---|---|---|---|---|---|---|---|
| 0.25 | 50 | 8 | 0.00e+00 | +0.0807 | +0.1279 | -0.9941 ± 0.0512 | -0.9091 |
| 0.25 | 200 | 8 | 0.00e+00 | +0.0723 | +0.1145 | -1.2921 ± 0.1633 | -0.9066 |
| 0.5 | 50 | 8 | 0.00e+00 | +0.0406 | +0.0644 | -0.9701 ± 0.0252 | -0.9462 |
| 0.5 | 200 | 8 | 0.00e+00 | +0.0425 | +0.0674 | -0.8241 ± 0.0568 | -0.9475 |
| 1.0 | 50 | 8 | 0.00e+00 | +0.0211 | +0.0335 | -1.0248 ± 0.0437 | -0.9667 |
| 1.0 | 200 | 8 | 0.00e+00 | +0.0194 | +0.0308 | -1.0090 ± 0.0544 | -0.9668 |

#### Reproduction line (sec. 1.7): Level 3 v6 estimator on k0.5_M200_u0.05

s̄ = **0.8710 ± 0.0078 σ** (229 periods in window) vs Level 3 0.8705 ± 0.0049: **0.06 σ -> PASS** (rule: within 2σ)

figure: 0000_PLAN_OVERALL/paper2_energytransfer/experiments/final/261012_p2_effmap.{png,pdf}

Figure: `261012_p2_effmap.{png,pdf}` ($\varepsilon_{\rm mean}$, the pre-registered figure).

### 2.2 Post-hoc summaries (written after 2.1; no pre-registered number changes)

$\varepsilon_{\rm mean}$'s last-$3T_w$ window leaves the divider's thermal motion unaveraged, so its errors are 5–20× those of $\varepsilon_{\rm settled}$. The pre-registered slope fit used $\varepsilon_{\rm mean}$. The tables below re-express the same cells through $\varepsilon_{\rm settled}$ and are labelled post hoc.

**Printed by `python3 hspist3/validation/paper2_effmap_posthoc_20261012.py`**, verbatim:

#### (a) POST HOC: epsilon_settled over the quasi-static cells u <= 0.05 (inverse-variance mean)

| k | M_s | mean eps_settled | chi2/dof across u | eps_rev | ratio to eps_rev | (mean - eps_rev)/sigma |
|---|---|---|---|---|---|---|
| 0.25 | 50 | 0.19335 ± 0.00099 | 0.6/2 | 0.1892 | 1.0222 | +4.2 |
| 0.25 | 200 | 0.19360 ± 0.00112 | 1.7/2 | 0.1892 | 1.0235 | +4.0 |
| 0.25 | both | 0.19346 ± 0.00074 | 2.3/5 | 0.1892 | 1.0228 | +5.8 |
| 0.5 | 50 | 0.10765 ± 0.00043 | 0.9/2 | 0.1062 | 1.0134 | +3.3 |
| 0.5 | 200 | 0.10707 ± 0.00052 | 3.3/2 | 0.1062 | 1.0079 | +1.6 |
| 0.5 | both | 0.10742 ± 0.00033 | 5.0/5 | 0.1062 | 1.0112 | +3.6 |
| 1.0 | 50 | 0.05765 ± 0.00028 | 2.6/2 | 0.0567 | 1.0168 | +3.4 |
| 1.0 | 200 | 0.05802 ± 0.00019 | 0.2/2 | 0.0567 | 1.0233 | +7.1 |
| 1.0 | both | 0.05791 ± 0.00016 | 4.0/5 | 0.0567 | 1.0212 | +7.8 |

#### (b) POST HOC: epsilon_settled slope on u in {0.1, 0.2, 0.5} (weighted linear fit)

| k | M_s | slope ± σ | intercept | chi2 (1 dof) | pre-registered eps_mean slope (for reference) |
|---|---|---|---|---|---|
| 0.25 | 50 | -0.0761 ± 0.0123 | 0.1950 | 0.0 | -0.0847 ± 0.0233 |
| 0.25 | 200 | -0.0711 ± 0.0133 | 0.1931 | 0.1 | -0.0139 ± 0.0350 |
| 0.5 | 50 | -0.0459 ± 0.0059 | 0.1113 | 2.7 | -0.0489 ± 0.0224 |
| 0.5 | 200 | -0.0481 ± 0.0074 | 0.1114 | 6.9 | -0.0863 ± 0.0221 |
| 1.0 | 50 | -0.0216 ± 0.0039 | 0.0582 | 1.7 | -0.0414 ± 0.0201 |
| 1.0 | 200 | -0.0223 ± 0.0036 | 0.0594 | 0.9 | -0.0144 ± 0.0169 |

#### (c) POST HOC: M_s = 50 vs 200, epsilon_settled, per (k, u)

| k | u | eps(50) | eps(200) | difference / sigma |
|---|---|---|---|---|
| 0.25 | 0.01 | 0.1938 | 0.1945 | -0.3 |
| 0.25 | 0.02 | 0.1937 | 0.1914 | +0.9 |
| 0.25 | 0.05 | 0.1920 | 0.1946 | -0.9 |
| 0.25 | 0.1 | 0.1875 | 0.1864 | +0.4 |
| 0.25 | 0.2 | 0.1796 | 0.1780 | +0.4 |
| 0.25 | 0.5 | 0.1571 | 0.1582 | -0.2 |
| 0.5 | 0.01 | 0.1079 | 0.1060 | +1.8 |
| 0.5 | 0.02 | 0.1069 | 0.1069 | -0.0 |
| 0.5 | 0.05 | 0.1079 | 0.1082 | -0.2 |
| 0.5 | 0.1 | 0.1073 | 0.1085 | -0.8 |
| 0.5 | 0.2 | 0.1000 | 0.0989 | +0.5 |
| 0.5 | 0.5 | 0.0895 | 0.0898 | -0.1 |
| 1.0 | 0.01 | 0.0577 | 0.0581 | -0.8 |
| 1.0 | 0.02 | 0.0581 | 0.0580 | +0.2 |
| 1.0 | 0.05 | 0.0569 | 0.0579 | -1.4 |
| 1.0 | 0.1 | 0.0564 | 0.0575 | -1.3 |
| 1.0 | 0.2 | 0.0527 | 0.0546 | -1.6 |
| 1.0 | 0.5 | 0.0479 | 0.0488 | -0.5 |

chi2 = 13.2 on 18 dof, p = 0.782

#### (d) POST HOC: where the quasi-static excess sits -- W_in and s_corr (last tau_r) at u <= 0.05

| k | M_s | <W_in>/W_rev | s_corr/s_rev | Delta E_spring/E_spring,rev |
|---|---|---|---|---|
| 0.25 | 50 | 1.0245 | 1.0415 | 1.0461 |
| 0.25 | 200 | 1.0223 | 1.0413 | 1.0459 |
| 0.5 | 50 | 1.0161 | 1.0256 | 1.0287 |
| 0.5 | 200 | 1.0209 | 1.0253 | 1.0284 |
| 1.0 | 50 | 1.0208 | 1.0323 | 1.0365 |
| 1.0 | 200 | 1.0190 | 1.0375 | 1.0423 |

#### (e) POST HOC: no-push controls, Delta KE_gas (last half) against -E_dof

| k | M_s | Delta KE_gas | -E_dof | z |
|---|---|---|---|---|
| 0.25 | 50 | -0.9941 ± 0.0512 | -0.9091 | -1.7 |
| 0.25 | 200 | -1.2921 ± 0.1633 | -0.9066 | -2.4 |
| 0.5 | 50 | -0.9701 ± 0.0252 | -0.9462 | -0.9 |
| 0.5 | 200 | -0.8241 ± 0.0568 | -0.9475 | +2.2 |
| 1.0 | 50 | -1.0248 ± 0.0437 | -0.9667 | -1.3 |
| 1.0 | 200 | -1.0090 ± 0.0544 | -0.9668 | -0.8 |

chi2 = 16.3 on 6 dof, p = 0.0121

figure: 0000_PLAN_OVERALL/paper2_energytransfer/experiments/final/261012_p2_effmap_both.{png,pdf}

### 2.3 Reading, item by item

1. **Ledger (A3): closed.** $\max|W_{\rm in} - (\Delta KE_{\rm gas} + KE_{\rm div} + \Delta E_{\rm spring})| = 1.35\times10^{-5}$ kT over every sample of all 288 pushed and 48 control runs, at the precision of the 6-decimal $x$. The `SpringE` row-0 gap is confirmed; from row 1 on, SpringE agrees with $\tfrac12 k(x - x_{\rm eq})^2$ to $1.6\times10^{-5}$.
2. **Reproduction line (§ 1.7): PASS.** $\bar s = 0.8710 \pm 0.0078\,\sigma$ against Level 3's $0.8705 \pm 0.0049$, i.e. 0.06σ. The map stands on Level 3's footing. Level 3 (contraction on) and the map (contraction off) are compared only as two numbers.
3. **Prediction 1, $\varepsilon \to \varepsilon_{\rm rev}$ with no $M_s$ dependence: holds to 1–2 %, and the residual is resolved.**
   - $\varepsilon_{\rm mean}$ agrees with $\varepsilon_{\rm rev}$ within 2σ in 10 of the 12 cells at $u \le 0.02$.
   - $\varepsilon_{\rm settled}$ resolves a small excess. Over $u \le 0.05$ it is $1.0228$, $1.0112$ and $1.0212 \times \varepsilon_{\rm rev}$ for $k = 0.25, 0.5, 1.0$ (+5.8σ, +3.6σ, +7.8σ; post hoc (a)).
   - **There is no $M_s$ dependence at any $u$:** $\chi^2 = 13.2$ on 18 dof, $p = 0.78$ (post hoc (c)).
   - The excess sits in both terms of the ratio: $W_{\rm in}$ is 1.6–2.5 % above $W_{\rm rev}$, and the settled $s$ is 2.5–4.2 % above $s_{\rm rev}$ (post hoc (d)).
   - **INFERENCE, not tested:** this is the box's finite-size over-pressure and over-stiffness, which the KR reversible reference omits. Level 3 measured $Z_{\rm box}/Z_{\rm KR} = 1.0248 \pm 0.0014$ in this same geometry. Recomputing $\varepsilon_{\rm rev}$ with Level 3's measured $F(L)$ would test it; that has not been done.
4. **Prediction 2, A2's shape: confirmed by $\varepsilon_{\rm settled}$.**
   - **Plateau.** $\varepsilon$ is flat for $u \le 0.05$ in all six $(k, M_s)$, with $\chi^2$ across $u$ of 0.2–3.3 on 2 dof. At $u = 0.1$ it is still on the plateau for $k = 0.5$, and just below it for $k = 0.25$ and 1.0.
   - **Decline.** Beyond that it falls roughly linearly. The slopes on $u \in \{0.1, 0.2, 0.5\}$ are $-0.076/-0.071$ ($k = 0.25$), $-0.046/-0.048$ ($k = 0.5$) and $-0.022/-0.022$ ($k = 1.0$), for $M_s = 50/200$. They fall with $k$ roughly in proportion to $\varepsilon_{\rm rev}$, and are about a third of the order-of-magnitude reference $-\varepsilon_{\rm rev}Z_{\rm ac}d/W_{\rm rev}$ (−0.24, −0.13, −0.071). That reference ignores the extra spring loading by the heated gas, which is in the direction observed.
   - **Pre-registered fit.** The slopes fitted on $\varepsilon_{\rm mean}$ (table 2.1) are consistent with these within their larger errors.
   - **Dissipation from the KR split.** $X + E_{\rm dof}$ is close to $Z_{\rm ac}ud$ for $u \ge 0.1$, within about a factor 1.5: $u = 0.5$ gives 11.2–11.6 against 8.84, $u = 0.2$ gives 2.5–4.1 against 3.54, and $u = 0.1$ gives 0.8–1.8 against 1.77. For $u \le 0.05$ it falls below the linear law (0.05–0.54 against 0.18–0.88), which is A2's quasi-static fall-off.
5. **Prediction 3, the controls: $W_{\rm in} = 0$ exactly in all 48 control runs.**
   - **Drift.** The control's drift $\langle s\rangle$ is +0.081/+0.072, +0.041/+0.043 and +0.021/+0.019 σ for $k = 0.25, 0.5, 1.0$. It scales as $1/k$, which is the pre-push offset Level 3 traced to the box's excess standing force.
   - **Spring-energy floor:** 0.03–0.13 kT.
   - **The $E_{\rm dof}$ INFERENCE (§ 1.8, A3) is right in sign and size but not in detail.** The control $\Delta KE_{\rm gas}$ is −0.82 to −1.29 against −0.91 to −0.97 predicted. Per cell, $\chi^2 = 16.3$ on 6 dof ($p = 0.012$), driven by two $M_s = 200$ cells at −2.4σ and +2.2σ with opposite signs, so not by a common offset. It is reported as a magnitude check, not a confirmation.
6. **Gate (A1 as implemented):** all 36 cells are UNRESOLVED, as expected from the 8-seed floor, and none is FAIL. The largest $E_{\rm coh}$ is $0.19 \pm 0.09$ at $(0.25, 200, 0.5)$, below the FAIL line. The construction part ($\ge 3\tau_r$ after the push, with Mansour's $\tau_r$ as an upper bound) carries the gate.

**The map, in one sentence.** The efficiency of storing piston work in the spring is set by $k$ alone, with no measurable $M_s$ dependence. It equals the reversible value (+1–2 %, the box effect) for pushes slower than one acoustic round trip ($d/u > 2L_0/c_s$), and falls roughly linearly in $u$ beyond it.

---

## 3. Over-pressure test — pre-registration (2026-10-13, written and committed before anything below was computed)

**The inference under test** (§ 2.3, item 3; commit 0504d68). The quasi-static plateau sits above the KR reversible reference:

$$\varepsilon_{\rm settled}/\varepsilon_{\rm rev} = 1.023 \pm 0.004,\; 1.011 \pm 0.003,\; 1.021 \pm 0.003 \qquad (k = 0.25,\ 0.5,\ 1.0).$$

This was attributed (INFERENCE) to the box's over-pressure, $Z_{\rm box}/Z_{\rm KR} = 1.0248$ (Level 3). The test replaces the guess by a number for each $k$: recompute $\varepsilon_{\rm rev}$ with the box's own measured equation of state instead of KR.

### Method

**F(L) source (DATA).** The files are `experiments_energy_transfer/level3_FofL_20260925/c{0,2.5,5,7.5,10}/ev_<seed>.csv` (20 seeds per compression) and the matching `tr_<seed>.csv`. These are geometry C with the divider held at mass $10^9$, the gas compressed by $d = 0, 1.99, 3.98, 5.97, 7.96$, so $L = 78.5 - d$. The recorded quantities and the code that writes them:

- `dp` in the event log is the **particle's** momentum change, `edmd.c:1069`:
  `fprintf(g_edmd_evlog, "%.12g,%s,%.12g,%.12g,%.12g,%.12g,%.12g\n", S->t / g_edmd_evlog_tscale, kind, u, v0, v1, dE, 1.0 * (v1 - v0));`
  The divider's events have kind `D0`. Event times include the 12 000-step hold (200 σ); trace times start at release.
- `KE_gas_total` is the sum of the segment kinetic energies, `00ALLINONE.c:17186`:
  `fprintf(elog, ",%.12e,%.12e,%.12e,%.12e", ke_tot, ke_left, ke_right, px_gas);`
  with `ke_tot += segment_ke[sgi]` (lines 17173–17175).

**Estimator.** For each seed,
$$F = \frac{1}{t_1 - t_0}\sum_{D0,\ t_0 \le t \le t_1}|dp|,\qquad T = \frac{\langle KE_{\rm gas}\rangle}{N},\qquad Z_{\rm box} = \frac{F L}{N T},\qquad \lambda = \frac{Z_{\rm box}}{Z_{\rm KR}(\eta)}.$$
The window is fixed by a principle, not a fit: it opens **two acoustic round trips after the piston stops**, $t_0 = t_{\rm stop} + 2\,(2L_0/c_s) = t_{\rm stop} + 180$ (with $t_{\rm stop} = 0.25/u + d/u$; for c0, $t_0 = 180$), and closes at the end of the record, $t_1 = 666.65$. Means are taken over seeds, with standard errors.

**Disclosure.** The 260925 table was computed inline and never committed. While preparing this section I scanned window starts from 0 to 450 σ after $t_{\rm stop}$ against it:
- $T$ reproduces exactly (1.0000 / 1.0347 / 1.0699 / 1.1070 / 1.1499) for every window after the stop.
- **$F$ is not reproduced to 4 decimals by any window.** The best is 0.002 off, and every post-stop window lies within ~0.005 of it, i.e. inside the published errors (0.004–0.008).
- So the published estimator is not recoverable. The table here is recomputed and printed beside the published one. Robustness of the $\lambda$ fit to window starts 100, 180 and 300 is a pre-registered check.

**Interpolation.** $\lambda(\eta)$ is a weighted linear fit over the five points, and $Z_F(\eta) = \lambda(\eta)\,Z_{\rm KR}(\eta)$. In `hspist3/validation/paper2_effmap_overpressure_20261013.py`, the pre-registered lines are:

    lam = lambda e: a + b * (e - ETA0)                       # <- the interpolation of Level 3's F(L) (pre-registered line)
    ZF = lambda e: lam(e) * Z(e)

The map's settled states ($\eta = 0.109$–$0.111$) lie inside Level 3's range ($0.10005$–$0.11134$), so there is no extrapolation.

**Reversible construction.** This is § 1.3 with $Z_F$ in place of $Z_{\rm KR}$, both in the force and in the isentrope:
$$d\ln T = -Z_F\,d\ln L \;(\text{exact for hard disks, whose energy is purely kinetic}),\qquad F(L) = N\,T_{\rm ad}(L)\,Z_F(\eta(L))/L,$$
with $k\,(x_{\rm eq} - (109 - d - L)) = F(L)$ solved for $L$, and $x_{\rm eq}$ as run.

With KR, § 1.3's starting point $x = 30.5$ **is** the gas–spring equilibrium (the $x_{\rm eq}$ were chosen for it). With $F(L)$ it is not: the box pushes ~2.5 % harder, the divider settles at a slightly different $x_i$ before the push, and the measurement subtracts exactly that offset through the no-push control. Two versions are therefore defined.

**The primary version is the estimator-matched construction, started from the box's own pre-push equilibrium $x_i$:**
$$s_F = x_i - x_f,\qquad \Delta E_{\rm est} = F_s s_F + \tfrac12 k s_F^2\;(F_s = k(x_{\rm eq} - 30.5),\ \text{the estimator's own formula, § 1.8}),$$
$$W_{\rm rev,F} = N\,[T_{\rm ad}(L_f) - T_{\rm ad}(L_i)] + \tfrac12 k(x_f - x_{\rm eq})^2 - \tfrac12 k(x_i - x_{\rm eq})^2,\qquad \varepsilon_{\rm rev,F} = \frac{\Delta E_{\rm est}}{W_{\rm rev,F}}.$$
$W_{\rm in}$ counts piston work only, and the relaxation $30.5 \to x_i$ involves none, so $W_{\rm rev,F}$ and $W_{\rm in}$ are on the same footing. With $\lambda \equiv 1$ this construction reproduces § 1.3 exactly; that identity is printed as a check.

**Errors.** $\sigma(\varepsilon_{\rm rev,F})$ comes from 400 draws of the $\lambda$-fit parameters (fit covariance, fixed seed). It is combined in quadrature with $\sigma(\varepsilon_{\rm settled})$.

**Measured plateau.** The inverse-variance mean of $\varepsilon_{\rm settled}$ over $u \le 0.05$ and both $M_s$, exactly as in § 2.2 (a).

### Verdict rules (written before the data)

- **PASS:** all three recomputed ratios $\varepsilon_{\rm settled}/\varepsilon_{\rm rev,F}$ lie within 2σ of 1.00. **The over-pressure explanation becomes DATA.**
- **FAIL:** any ratio lies outside 2σ. **It stays OPEN**, and the paper says "a 1–2 % unexplained offset". At most three candidate causes are listed, as OPEN, with no new runs.

**Consistency check (pre-registered).** $\chi^2$ of the three measured $\varepsilon_{\rm settled}/\varepsilon_{\rm rev,KR}$ against one common value (2 dof). A single factor $\lambda$ predicts $k$-independence only to first order: the spring preload $F_s$ is fixed by $x_{\rm eq}$, not by the gas. So the recomputation also prints $\varepsilon_{\rm rev,F}/\varepsilon_{\rm rev,KR}$ per $k$, which is the actual prediction of "over-pressure only".

**Secondary, no verdict.**
- (i) The literal § 1.3 construction with $Z_F$, started at $x = 30.5$.
- (ii) The primary construction with the gas started at the measured control temperature $T_i = \langle KE_{\rm gas}\rangle/N$ (late half of the no-push controls). This covers the divider's equipartition share, § 1.8 A3.

### The settle gate, replaced (§ 3.2)

The $KE_{\rm div}$ gate measured physics, not settling: the divider always carries $kT/2$. The replacement is **position-settled**:
$$\Delta = \big|\langle \bar x\rangle_{\text{last } P_m} - \langle \bar x\rangle_{\text{previous } P_m}\big| < 0.01\, s_{\rm rev}(k),$$
- $\bar x(t)$ is the 8-seed mean trajectory.
- $s_{\rm rev} = 1.4987, 0.8403, 0.4480$ is the reversible displacement of § 1.3, the plan's "$x_{\rm rev}$".
- $P_m$ is the mode period, from the exact one-column eigen-equation at the settled state:
$$M_s\,\omega^2 = k + k_S\,K\cot K,\qquad K = \frac{\omega L_f}{c_s},\qquad k_S = \frac{N m c_s^2}{L_f^2},$$
with KR $c_s$ at $(\eta_f, T_f)$ of § 1.3. For the controls, $(\eta_0, T = 1)$.

Every one of the 36 cells gets PASS or FAIL. **Calibration (reported, not a reclassification):** the same gate is applied to the six no-push controls, which are settled by construction. Their pass rate measures the gate's noise floor. The map figures are re-rendered with the label, and nothing on them says "UNRESOLVED".

### 3.1 Results — over-pressure test (run once, after the § 3 commit f3c6209)

**Printed by `python3 hspist3/validation/paper2_effmap_overpressure_20261013.py`**, verbatim:

#### 1. Level 3 F(L), recomputed from the raw files (window t_stop + 180 -> end), against the 260925 table

| cell | L | eta | F (here) | F (260925) | T (here) | T (260925) | Z_box (here) | Z_box (260925) | lambda = Z_box/Z_KR |
|---|---|---|---|---|---|---|---|---|---|
| c0 | 78.50 | 0.10005 | 1.6135 ± 0.0051 | 1.6087 ± 0.0064 | 1.0000 | 1.0 | 1.2666 ± 0.0040 | 1.2628 | 1.0245 ± 0.0033 |
| c2.5 | 76.51 | 0.10265 | 1.7208 ± 0.0035 | 1.7218 ± 0.0078 | 1.0347 | 1.0347 | 1.2725 ± 0.0020 | 1.2732 | 1.0232 ± 0.0016 |
| c5 | 74.52 | 0.10539 | 1.8378 ± 0.0053 | 1.8366 ± 0.0067 | 1.0699 | 1.0699 | 1.2800 ± 0.0033 | 1.2792 | 1.0230 ± 0.0026 |
| c7.5 | 72.53 | 0.10829 | 1.9699 ± 0.0070 | 1.9715 ± 0.0059 | 1.1070 | 1.107 | 1.2906 ± 0.0045 | 1.2917 | 1.0247 ± 0.0035 |
| c10 | 70.54 | 0.11134 | 2.1303 ± 0.0058 | 2.1217 ± 0.0042 | 1.1499 | 1.1499 | 1.3069 ± 0.0034 | 1.3016 | 1.0304 ± 0.0026 |

lambda(eta) = 1.02189 ± 0.00181 + (+0.587 ± 0.301)(eta - 0.10005098); chi2 = 2.09 on 3 dof; pooled constant lambda = 1.02471

Robustness of the fit to the window start (pre-registered check):

| window start after t_stop | a = lambda(eta_0) | b |
|---|---|---|
| 100 | 1.02246 ± 0.00206 | +0.187 ± 0.339 |
| 180 | 1.02189 ± 0.00181 | +0.587 ± 0.301 |
| 300 | 1.01990 ± 0.00338 | +0.624 ± 0.477 |

#### 2. The test (primary: estimator-matched construction from the box's own pre-push equilibrium)

| k | eps_rev KR (sec. 1.3) | eps_rev KR, matched | **eps_rev F(L)** | eps_settled (u <= 0.05) | ratio_KR | **ratio_F** | (ratio_F - 1)/σ | verdict | eps_rev F/KR |
|---|---|---|---|---|---|---|---|---|---|
| 0.25 | 0.18915 | 0.18916 | **0.19445 ± 0.00199** | 0.19346 ± 0.00074 | 1.0228 ± 0.0039 | **0.9949 ± 0.0109** | -0.5 | within 2σ | 1.0280 |
| 0.5 | 0.10623 | 0.10624 | **0.10997 ± 0.00127** | 0.10742 ± 0.00033 | 1.0112 ± 0.0031 | **0.9768 ± 0.0117** | -2.0 | within 2σ | 1.0352 |
| 1.0 | 0.05670 | 0.05670 | **0.05897 ± 0.00073** | 0.05791 ± 0.00016 | 1.0212 ± 0.0027 | **0.9820 ± 0.0124** | -1.4 | within 2σ | 1.0400 |

**VERDICT (pre-registered): PASS -- the over-pressure explanation becomes DATA**

Consistency check (pre-registered): measured ratio_KR against one common value: mean 1.0181 ± 0.0018, chi2 = 7.67 on 2 dof, p = 0.0216

#### 3. Secondary variants (no verdict)

| k | literal sec. 1.3 construction with F(L): eps_rev | ratio | matched + measured control T_i: T_i | eps_rev | ratio |
|---|---|---|---|---|---|
| 0.25 | 0.21033 | 0.9198 | 0.98857 | 0.19484 | 0.9929 |
| 0.5 | 0.11874 | 0.9047 | 0.99103 | 0.11003 | 0.9763 |
| 1.0 | 0.06361 | 0.9103 | 0.98983 | 0.05895 | 0.9823 |

Figure: `261013_p2_effmap_overpressure.{png,pdf}`.

**Reading.**

- **The F(L) input (DATA).** Recomputed from the raw files, $T$ reproduces the 260925 table exactly, and $F$ lies within the published errors. The fit is $\lambda(\eta) = 1.02189 \pm 0.00181 + (0.59 \pm 0.30)(\eta - \eta_0)$, with $\chi^2 = 2.09$ on 3 dof. Moving the window start between 100 and 300 σ shifts $\lambda(\eta_0)$ within its error.
- **Verdict by the pre-registered rule: PASS.** The three recomputed ratios are $0.995 \pm 0.011$, $0.977 \pm 0.012$ and $0.982 \pm 0.012$, at $-0.5\sigma$, $-2.0\sigma$ and $-1.4\sigma$. **The over-pressure explanation is tagged DATA**, with this precise content: *within the precision of Level 3's measured $F(L)$, the box's over-pressure accounts for the quasi-static plateau excess.* Three qualifications go with it:
  1. **The pass rests on the reference error.** $\sigma(\varepsilon_{\rm rev,F})$ is 1.0–1.2 %, three to four times the measurement error. The three ratios share the same $\lambda$-fit draws, so they are correlated, not three independent passes.
  2. **The correction overshoots.** All three ratios lie *below* 1. The F(L) reference raises $\varepsilon_{\rm rev}$ by +2.8, +3.5 and +4.0 %, against an observed excess of +2.3, +1.1 and +2.1 %. $k = 0.5$ sits on the 2σ line.
  3. **The pattern in $k$ is not explained (OPEN).** The measured ratios are not consistent with one common value: $\chi^2 = 7.67$ on 2 dof, $p = 0.022$, driven by the dip at $k = 0.5$. Over-pressure alone predicts a ratio rising monotonically with $k$ ($\varepsilon_{\rm rev,F}/\varepsilon_{\rm rev,KR} = 1.028, 1.035, 1.040$), and the data do not follow it. The $k = 0.5$ runs used $x_{\rm eq} = 33.65$ rather than 33.6498, a preload difference of $10^{-4}$, which is too small to matter.
- **Secondary variants.**
  - The literal § 1.3 construction, started at $x = 30.5$, gives 0.90–0.92. This confirms that the pre-push offset, which the control removes, has to be matched.
  - Starting the gas at the measured control temperature $T_i = 0.989$–$0.991$ changes the ratios by at most 0.002.

**For the paper:** the quasi-static efficiency agrees with the reversible value computed from the box's own measured equation of state, within that equation of state's 1 % precision. The residual $k$-dependence at the 1 % level is stated as open.

### 3.2 Results — position-settled gate (pre-registered in § 3; no new runs)

**Printed by `python3 hspist3/validation/paper2_effmap_gate_20261013.py`**, verbatim. The last table is a post-hoc diagnostic, labelled as such.

#### Position-settled gate, all 36 cells (Delta in sigma; threshold 0.01 s_rev)

| k | M_s | u | P_m | Delta | 0.01 s_rev | label |
|---|---|---|---|---|---|---|
| 0.25 | 50 | 0.01 | 112.47 | 0.09614 | 0.01499 | FAIL |
| 0.25 | 50 | 0.02 | 112.47 | 0.03102 | 0.01499 | FAIL |
| 0.25 | 50 | 0.05 | 112.47 | 0.11372 | 0.01499 | FAIL |
| 0.25 | 50 | 0.1 | 112.47 | 0.05232 | 0.01499 | FAIL |
| 0.25 | 50 | 0.2 | 112.47 | 0.12445 | 0.01499 | FAIL |
| 0.25 | 50 | 0.5 | 112.47 | 0.04698 | 0.01499 | FAIL |
| 0.25 | 200 | 0.01 | 172.17 | 0.07130 | 0.01499 | FAIL |
| 0.25 | 200 | 0.02 | 172.17 | 0.08537 | 0.01499 | FAIL |
| 0.25 | 200 | 0.05 | 172.17 | 0.03096 | 0.01499 | FAIL |
| 0.25 | 200 | 0.1 | 172.17 | 0.02271 | 0.01499 | FAIL |
| 0.25 | 200 | 0.2 | 172.17 | 0.01207 | 0.01499 | PASS |
| 0.25 | 200 | 0.5 | 172.17 | 0.05158 | 0.01499 | FAIL |
| 0.5 | 50 | 0.01 | 92.74 | 0.19781 | 0.00840 | FAIL |
| 0.5 | 50 | 0.02 | 92.74 | 0.00580 | 0.00840 | PASS |
| 0.5 | 50 | 0.05 | 92.74 | 0.12182 | 0.00840 | FAIL |
| 0.5 | 50 | 0.1 | 92.74 | 0.03171 | 0.00840 | FAIL |
| 0.5 | 50 | 0.2 | 92.74 | 0.08829 | 0.00840 | FAIL |
| 0.5 | 50 | 0.5 | 92.74 | 0.20781 | 0.00840 | FAIL |
| 0.5 | 200 | 0.01 | 130.02 | 0.06728 | 0.00840 | FAIL |
| 0.5 | 200 | 0.02 | 130.02 | 0.01995 | 0.00840 | FAIL |
| 0.5 | 200 | 0.05 | 130.02 | 0.02276 | 0.00840 | FAIL |
| 0.5 | 200 | 0.1 | 130.02 | 0.07469 | 0.00840 | FAIL |
| 0.5 | 200 | 0.2 | 130.02 | 0.02105 | 0.00840 | FAIL |
| 0.5 | 200 | 0.5 | 130.02 | 0.00610 | 0.00840 | PASS |
| 1.0 | 50 | 0.01 | 82.04 | 0.06280 | 0.00448 | FAIL |
| 1.0 | 50 | 0.02 | 82.04 | 0.05319 | 0.00448 | FAIL |
| 1.0 | 50 | 0.05 | 82.04 | 0.17001 | 0.00448 | FAIL |
| 1.0 | 50 | 0.1 | 82.04 | 0.00270 | 0.00448 | PASS |
| 1.0 | 50 | 0.2 | 82.04 | 0.11607 | 0.00448 | FAIL |
| 1.0 | 50 | 0.5 | 82.04 | 0.01466 | 0.00448 | FAIL |
| 1.0 | 200 | 0.01 | 98.26 | 0.05638 | 0.00448 | FAIL |
| 1.0 | 200 | 0.02 | 98.26 | 0.14157 | 0.00448 | FAIL |
| 1.0 | 200 | 0.05 | 98.26 | 0.02206 | 0.00448 | FAIL |
| 1.0 | 200 | 0.1 | 98.26 | 0.02789 | 0.00448 | FAIL |
| 1.0 | 200 | 0.2 | 98.26 | 0.02834 | 0.00448 | FAIL |
| 1.0 | 200 | 0.5 | 98.26 | 0.00320 | 0.00448 | PASS |

**PASS count per (k, M_s):** k = 0.25, M_s = 50: 0/6; k = 0.25, M_s = 200: 1/6; k = 0.5, M_s = 50: 1/6; k = 0.5, M_s = 200: 1/6; k = 1.0, M_s = 50: 1/6; k = 1.0, M_s = 200: 1/6. **Total 5/36.**

#### Calibration: the same gate on the no-push controls (settled by construction)

| k | M_s | P_m (eta_0, T = 1) | Delta | 0.01 s_rev | label |
|---|---|---|---|---|---|
| 0.25 | 50 | 120.45 | 0.08686 | 0.01499 | FAIL |
| 0.25 | 200 | 178.17 | 0.01015 | 0.01499 | PASS |
| 0.5 | 50 | 103.33 | 0.12421 | 0.00840 | FAIL |
| 0.5 | 200 | 134.23 | 0.07251 | 0.00840 | FAIL |
| 1.0 | 50 | 95.61 | 0.01225 | 0.00448 | FAIL |
| 1.0 | 200 | 104.78 | 0.05309 | 0.00448 | FAIL |

controls PASS: 1/6

#### DIAGNOSTIC, post hoc, not a gate: the controls' Delta against window length (noise floor of the rule)

| k | M_s | 0.01 s_rev | Delta, 1 P_m | 10 P_m | 30 P_m | tau_r |
|---|---|---|---|---|---|---|
| 0.25 | 50 | 0.01499 | 0.08686 | 0.00797 | 0.00719 | 0.00216 |
| 0.25 | 200 | 0.01499 | 0.01015 | 0.04544 | 0.01545 | 0.01168 |
| 0.5 | 50 | 0.00840 | 0.12421 | 0.02248 | 0.00621 | 0.00316 |
| 0.5 | 200 | 0.00840 | 0.07251 | 0.02027 | 0.00167 | 0.00443 |
| 1.0 | 50 | 0.00448 | 0.01225 | 0.00102 | 0.00231 | 0.00142 |
| 1.0 | 200 | 0.00448 | 0.05309 | 0.00711 | 0.00659 | 0.00115 |

Figures re-rendered with the label: `261012_p2_effmap.{png,pdf}`, `261012_p2_effmap_both.{png,pdf}`, and the per-k panels `261012_p2_effmap_k{0.25,0.5,1.0}.{png,pdf}`. FAIL cells carry a black ×. Every figure states the control calibration, and nothing on them says "UNRESOLVED".

**Reading.**

- **The gate as pre-registered (DATA): 5 of 36 cells PASS.**
- **It is noise-limited (DATA).** The six no-push controls are settled by construction, yet only 1 of 6 passes. For a settled system, the one-period mean of the 8-seed trajectory moves by 0.010–0.124 σ, which is 1–28× the thresholds of 0.0045–0.015 σ. **A FAIL from this gate therefore does not indicate an unsettled divider**, and the label cannot serve as the paper figure's settle criterion. A one-period window does not average out the divider's thermal motion.
- **What would work (post-hoc diagnostic, not adopted).** With windows one $\tau_r$ long — the windows $\varepsilon_{\rm settled}$ already uses — the controls' $\Delta$ is 0.0011–0.0117 σ, under threshold in 6 of 6. A "last $\tau_r$ vs previous $\tau_r$" version of the same criterion would be a usable gate. Adopting it after seeing these data would be a post-hoc change, so **it is not applied here; the decision is the plan author's (OPEN).**
- **Unaffected:** the $\varepsilon$ values, their errors, the § 3.1 verdict and every number in § 2. The gate is a label, not an input.
