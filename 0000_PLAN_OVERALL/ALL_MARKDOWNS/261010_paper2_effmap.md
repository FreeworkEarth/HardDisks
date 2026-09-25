# Efficiency map, geometry C — PRE-REGISTRATION

2026-10-10. **Not launched.** Written ledger-first, as agreed: the exact energy accounting comes
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
sample, exactly as Level 4b did (which closed to 1.0e-6 kT).

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
   So epsilon should fall roughly linearly in u with slope set by Z_ac d/W_rev.
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

## 2. Results

*Empty. Not launched. Chris and the plan author read section 1 first.*
