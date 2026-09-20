# Level 3 — can a spring capture the work? Master box, geometry C

2026-09-21. Supersedes `260920_paper2_level3_REPORT.md`, whose numbers are withdrawn: the pilot's
stated design rested on a units error (§1). Runs: `experiments_energy_transfer/level3_master_20260920`
(v1, superseded) and `level3_master_preload_20260921` (v2, the result). Analysis:
`validation/paper2_level3_preload_20260921.py`, design `validation/level3_design_20260920.py`.

**Verdict: Level 3 is NOT passed.** Criteria (i) and (ii) are met; (iii) is not, and the reason is
physical rather than statistical — the impedance formula is not the controlling physics in this
apparatus. Three separate defects were found and fixed along the way, and each is a result in its
own right.

---

## 1. The pilot's design was wrong by a factor of 24 (units)

`--spring-k` is **per pixel², not per σ²**. Confirmed two independent ways on a stored pilot trace:

| check | measured | k = 5 as kT/σ² | k = 5 × 24² = 2880 kT/σ² |
|---|---|---|---|
| equipartition, ⟨SpringE⟩ = 0.526 kT with wall rms 0.01910 σ → implied k | **2885** | 5 | 2880 |
| wall oscillation, FFT peak f = 0.6033 → T_w | **1.657 σ-time** | 39.54 | 1.656 |

So the pilot's "T_w = 39.5 σ-time, τ_push/T_w spans 19.9 … 0.5" was wrong by 24×. The true span was
**474 … 12: every pilot cell was deep in the quasi-static regime and none was ever impulsive.** The
two limits the pilot reported spanning were never spanned. `SpringE` itself is in kT and was always
correct — only k's units were misread.

Fixed by adding `--spring-k-sigma=K` (kT/σ²); `--spring-k` is untouched so every pre-existing
command still means what its recorded command line says. It had to be inserted *above* the
`--spring-k` parse branch, which prefix-matches on 10 characters and silently swallowed it — the
same trap that already made `--seeding` dead code behind `--seed`. Every run now prints k in both
units plus ω_w, T_w and the thermal rms, so the error is visible on line one.

## 2. Arm (a) — the quasi-static optimum — is structurally impossible

The cell k = k_gas aborted 2/2 seeds with `[wall_boundary_contact]`. It is not under-roomed. A
spring anchored at the wall's start must hold the standing gas force F = N kT Z/L, so the wall sits
d = F/k from the anchor, and at k = k_gas

> **F / k_gas = (N kT Z/L) / (N kT (Z + ηZ′)/L²) = L · Z/(Z + ηZ′) = 0.816 L, independent of N.**

The stiffness that maximises quasi-static capture is exactly the stiffness at which the spring
cannot be statically balanced inside its own apparatus, and no compartment size fixes it because the
requirement scales with the gas length. Measured static offsets:

| k [kT/σ²] | 0.0246 (= k_gas) | 0.1 | 0.5 | 2 | 5 |
|---|---|---|---|---|---|
| offset F/k [σ] | **64.1** | 15.8 | 3.15 | 0.79 | 0.31 |
| thermal rms [σ] | 4.51 | 2.83 | 1.38 | 0.70 | 0.45 |
| pedestal rms F·σ_δ [kT] | 7.10 | 4.46 | 2.17 | 1.11 | 0.70 |

Pre-loading to compensate does not rescue it either: the baseline then couples linearly and the
pedestal fluctuation becomes 7.1 kT against a 0.195 kT signal. Arm (a) was dropped and the reason
recorded; the run script keeps the cell commented rather than silently absent.

## 3. v1 ran clean and measured the wrong thing (the start-up transient)

600/600, 0 aborts, 0 health events — and invalid. `--spring-eq` was set to the wall's start, so the
spring began at its natural length while the gas was **already pushing on it** with F = 1.575 kT/σ.
The wall was released into that unbalanced load and rang:

- seed-mean wall swing **7.5 σ** (30.500 → 23.006), peak spring energy **14.0 kT**
- against the **1.04 kT** the piston launches beyond W_qs at u = 0.2

The transient is identical in every seed, so ensemble averaging does not remove it — it *was* the
signal v1 measured, and it is why v1's ε_coh rose monotonically with M_s. Fixed by anchoring the
spring at wall + F/k = 33.65 σ so it balances the standing force at t = 0; the swing falls to 3.32 σ,
which is ~2.4 thermal rms, i.e. mostly thermal.

## 4. The measurement (v2, pre-loaded)

Geometry C in the 30 σ master box: one gas of 100 disks over 78.5 σ (η = 0.1001), Δx = 7.96 σ
(10 % compression), k = 0.5 kT/σ², step protocol, 40 seeds per cell, 600 runs, **0 aborts, 0 health
events**. Primary observable is the **ensemble-mean** E_spring(t); the per-seed maximum is the wrong
estimator and the pilot showed why.

M_s is scanned at fixed k because that is the discriminator: the series formula contains no M_s and
predicts a flat line, the impedance formula T_imp = 4Z_wZ_g/(Z_w+Z_g)², Z_w = √(M_s(k+k_gas)),
Z_g = ρc_sH = 2.222, predicts an interior maximum at M_s = 9.41.

| u | M_s = 2 | 10 | 50 | 200 | 1000 | peak at |
|---|---|---|---|---|---|---|
| 0.2 | 0.300 | 0.304 | 0.337 | **0.348** | 0.263 | 200 |
| 0.5 | 0.419 | 0.441 | **0.462** | 0.374 | 0.250 | 50 |
| 1.0 | 0.336 | 0.380 | **0.456** | 0.336 | 0.237 | 50 |
| T_imp predicts | 0.864 | **1.000** | 0.844 | 0.586 | 0.322 | 9.41 |

(ε_coh = E_coh/⟨W_in⟩. Full table with E_coh, t_coh, significance and pedestals in the script output.)

## 5. Against the criteria

**(i) Ledger with the spring term — met in substance.** Over all 600 runs the residual
|W − (ΔKE_gas + ΔKE_wall + ΔE_spring)|/W has **median 4.64e-10**, 99th percentile 1.08e-08, worst
2.83e-08; **593/600 are below 1e-8** and the 7 exceptions are concentrated at M_s = 2, where the
wall KE term is smallest relative to the total. Reaching this took two precision fixes (spring
diagnostics and wall velocities float → double, both gated as precision-only) and one analysis fix:
the spring is armed one step after the first trace row, so `SpringE[0]` is a spurious exact zero,
and referencing the ledger to it inflated the residual to 2e-1 — an analysis bug, not physics.

**(ii) E_coh above the pedestal — PASS.** The coherent peak stands 6.8σ to 31σ above the pre-signal
level in all 15 cells with u ≥ 0.2. The baseline is taken from the causally blind window: the wall
sits 78.5 σ from the piston, so nothing the piston does can reach it for L/c_s = 45 σ-time.

**(iii) ε_coh maximal at the impedance match — FAIL.** There *is* an interior maximum, which the
series formula cannot produce, so the M_s scan did its job as a discriminator. But it sits at
**M_s ≈ 50 against a predicted 9.41**, a factor 5 in mass and 2.2 in impedance, and it does not move
with u the way an acoustic term should.

## 6. What the data say instead, and the open question

E_coh is **larger than the acoustic excess it was supposed to measure**: 4.9 kT at u = 0.2 against
A u² = 1.04 kT. So the spring is not mainly catching the pulse — it is catching the quasi-static
pressure rise, and it can do so because it is pre-loaded.

That is the open question, and it contradicts the ceiling derivation now in the plan. With the
spring **unloaded**, E = ½kx² and the capture is second order in the wall displacement. With the
spring **pre-loaded**, E = ½k(d + s)² so

> ΔE = F·s + ½k s² — **first order** in s, with F = k d the standing force.

At the quasi-static wall shift this is 0.67 kT against the unloaded formula's 0.035 kT, about 19×.
The "capture is second order in the compression fraction f" argument therefore holds only for an
unloaded spring, which is not the apparatus anyone would build. The plan's §Level 3 ceiling has
**not** been rewritten pending a decision, because the correction changes its conclusion rather than
its arithmetic.

Two consequences if the first-order reading stands:
1. The ε_coh ≈ 0.25–0.46 measured above is a real capture fraction, and a large one — but it is
   quasi-static capture by a pre-loaded spring, not the acoustic capture Level 3 set out to test.
2. Separating the two needs an observable that isolates the pulse: the excess over the *u → 0*
   limit of E_coh at fixed M_s, which subtracts the quasi-static term the way Level 1 subtracts
   W_qs from W. That is one more speed ladder per M_s, not a new apparatus.

## 7. Provenance

Binary sha1 2ec52910 → (spring hygiene, render modes) → (double spring diagnostics) → (double wall
velocities), each step gated with `validation/binary_gate.sh`: data files byte-identical, declared
columns verified as precision-only by rounding the new doubles back to float32 and requiring
agreement within 8 ULPs at the family's scale. Master box Level 0b at the 30 σ compartment:
⟨W_in⟩ = 7.5413 ± 0.0573 against geometry A's 7.4841 ± 0.0381 over 90 seeds = 0.83σ, ledger 1.74e-05.
