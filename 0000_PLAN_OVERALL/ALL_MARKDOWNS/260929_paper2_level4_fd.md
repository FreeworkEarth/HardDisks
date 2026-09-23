# Level 4 — τ_T from equilibrium fluctuations (fluctuation–dissipation)

2026-09-29. `experiments_energy_transfer/level4_equilibrium_20260929/Md10`, **20 runs, 0 aborts, 0
health events**, record 10 195 σ-time each, sampled 1-in-300 to Δt = 5.0 σ-time. Two-compartment
box, N_s = 50 per side, l0 = 39.25, height 10, divider free at x = 39.25 with M_d = 10, held 12 000
steps then released. η = 0.1001. **No temperature step, no piston.**

**τ_T is measured for the first time: 490 ± 190 (from T₁−T₂) and 450 ± 165 (from the divider
position), against 504 predicted for hard disks and 618 for an ideal gas. It is consistent with
both. The 23 % gap is not resolved, and this report says why and what would resolve it.**

The previous report's conclusion — *"needs ≈ 2300 seeds"* — **was wrong**, and wrong in an
expensive direction. The binding constraint is record length, not seed count.

---

## 1. The premise is confirmed

The per-seed spread of T₁ − T₂ is not measurement noise. Hard disks have no potential energy, so a
compartment's kinetic temperature *is* its energy, the microcanonical split is f = E₁/E ~ Beta(N,N),
and σ_f = 1/(2√(2N+1)) = 0.0497, giving σ(T₁−T₂) = 4T̄σ_f = **0.199**.

**Measured over 20 seeds × 1 639 samples: 0.1980.** Three digits. Relaxation and fluctuation are
one process, so the stationary autocorrelation carries τ_T and no temperature step is needed. That
is the whole justification for this run, and it holds.

## 2. The first estimator was wrong, and its failure was diagnostic

A single-exponential fit to the autocorrelation returned **288 ± 16** from T₁−T₂ and **335 ± 10**
from the divider position. Both sat *below* the OU control's prediction for either hypothesis, and
they disagreed with each other by 3.5σ. Two observables of the same process cannot disagree; the
estimator was at fault, not the data.

**The autocorrelation is not a single exponential.** Sampled properly it reads

| lag | 0 | 25 | 45 | 50 | 70 | 90 | 100 | 150 | 200 | 300 | 500 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ACF(T₁−T₂) | 1.000 | 0.350 | −0.009 | 0.015 | 0.400 | 0.686 | 0.625 | 0.118 | 0.388 | 0.238 | 0.072 |
| ACF(x) | 1.000 | 0.393 | 0.004 | 0.010 | 0.349 | 0.657 | 0.625 | 0.097 | 0.399 | 0.249 | 0.086 |

It falls to zero by lag 45 and comes **back up to 0.66** by lag 90. The power spectrum resolves it
into three parts:

| component | weight | what it is |
|---|---|---|
| period > 200 σ-time | **45.1 %** (x), 43.5 % (ΔT) | the slow isobaric mode — the one τ_T describes |
| discrete peak at **period 95.3** | ≈ 30 % | an oscillation, see §3 |
| remainder | ≈ 25 % | fast collisional noise |

A single exponential fitted across all three returns a meaningless average of them, which is
precisely why the two observables disagreed.

**It is not a start-up transient.** The divider is released from rest at t ≈ 200, so a coherent
ring-down was the obvious suspect. It is not: the autocorrelation computed on t > 2000 alone is the
same as on t < 2000 (0.65 vs 0.61 at lag 90), and the 20-seed mean of both observables scatters
about zero throughout with no coherent decay. The oscillation is a stationary equilibrium mode.

## 3. An unexplained factor of 2.1 — flagged, not explained

The adiabatic gas spring predicts a divider oscillation at

> k_gas = N kT c_s²/L² = 0.0988 per side (η = 0.1001, Z = 1.2363, c_s² = 3.0425),
> ω = √(2k_gas/M_d) = 0.1405, **period 44.7 σ-time**.

The measured peak is at **period 95.3** — a factor **2.13** too slow, i.e. an effective stiffness
0.22× the adiabatic two-sided value, or an effective mass 45 rather than 10. An isothermal spring
(dropping Z² from c_s²) gives 63.4 and does not close the gap either.

**No mechanism is claimed.** Hydrodynamic added mass on a piston in a closed tube is the obvious
candidate and would have to supply ≈ 35 m, but nothing here tests that. It is recorded as an open
item because it carries ~30 % of the variance of both observables and any future estimator must
handle it.

## 4. The estimator that works, and its calibration

Block-averaging the series at **100 σ-time** — one block ≈ one oscillation period — removes the
period-95 component by construction. The block autocorrelation is then a clean single exponential:

| lag [σ] | 0 | 100 | 200 | 300 | 400 | 500 | 600 | 700 |
|---|---|---|---|---|---|---|---|---|
| ACF(T₁−T₂) | 1.000 | 0.803 | 0.596 | 0.431 | 0.310 | 0.219 | 0.163 | 0.117 |
| ACF(x) | 1.000 | 0.768 | 0.572 | 0.426 | 0.314 | 0.226 | 0.180 | 0.115 |

successive ratios 0.80, 0.74, 0.72, 0.72, 0.71, 0.74, 0.72 — constant to 4 %.

**The two observables now agree**: τ = **302 ± 65** (T₁−T₂) and **291 ± 55** (x), jackknifed over
seeds, against 288 vs 335 before. Agreement between two independent observables of the same process
is the check that the estimator is no longer mixing modes.

These are *biased low*. The identical procedure — same record length, same blocking, same fit
window, same 20 seeds — run on synthetic Ornstein–Uhlenbeck series of known τ (40 trials each):

| true τ | 250 | 350 | 450 | 504 | 550 | 618 | 700 | 900 |
|---|---|---|---|---|---|---|---|---|
| recovered | 194 | 252 | 292 | **306** | 311 | **349** | 376 | 412 |
| spread (1 trial) | 23 | 29 | 39 | 37 | 36 | 41 | 53 | 48 |

The response slope is **0.335**: the estimator compresses differences threefold at this record
length. Inverting it,

> **τ_T(M_d = 10) = 490 ± 190** from T₁ − T₂
> **τ_T(M_d = 10) = 450 ± 165** from the divider position

| hypothesis | prediction | control gives | measured | |
|---|---|---|---|---|
| hard disk, 50.4 M | 504 | 306 ± 37 | 302 ± 65 | **−0.1σ** |
| ideal, 61.8 M | 618 | 349 ± 41 | 302 ± 65 | **−0.6σ** |

**Consistent with both.** A 23 % gap cannot be resolved with a 35–40 % error bar. No preference
between Kolafa–Rottner and the ideal gas is claimed here — unlike stage 1, where the same run family
separates them at 7–8σ (`260928_paper2_level4_thermal.md`).

What *is* established is that τ_T at M_d = 10 is of order 500 σ-time, which is the first direct
measurement of the adiabatic-piston stage-2 timescale in this apparatus, and is consistent with the
linear-in-M law taken from the literature.

## 5. The KOA design — record length, not seeds

The previous report asked for ≈ 2300 seeds. That was the wrong lesson: it was derived from a
signal-to-noise argument on a *decaying step*, and the equilibrium route has no step to lose. The
controlling quantity is L/τ, because the bias — not the scatter — is what hides the 23 % gap, and
bias is common-mode across seeds so seeds cannot remove it.

Separation between the two hypotheses for one experiment of the stated size:

| record [σ-time] | L/τ | seeds | slope | separation |
|---|---|---|---|---|
| 8 200 (this run) | 16 | 20 | 0.29 | 0.8σ |
| 8 200 | 16 | **80** | 0.33 | 1.7σ |
| **16 400** | 33 | 20 | **0.78** | 1.7σ |
| 32 800 | 65 | 20 | 0.77 | 1.8σ |
| **32 800** | **65** | **80** | **0.85** | **3.3σ** |
| 65 600 | 130 | 20 | 0.93 | 2.5σ |

Quadrupling the record buys more than quadrupling the seeds (0.29 → 0.78 in slope on the first
doubling alone), and the two together reach 3σ.

> **Recommended:** M_d ∈ {10, 50, 200}, record **4× this one** (≈ 65 τ_T at each mass, i.e. 32 800,
> 164 000, 656 000 σ-time), **80 seeds** each. At M_d = 10 that is ≈ 2.0 M steps × 80 seeds, about
> 13× the cost of this run — an overnight job, not a cluster campaign.

The mass ladder is what tests the *exponent* as well as the coefficient; this run has only M_d = 10
and therefore tests neither.

## 6. Route B (stationary heat conduction) — blocked, and the design is recorded

The alternative to a decaying signal is a *stationary* one: hold the two outer walls at T_hot and
T_cold, let the box reach a steady state, and read the conductivity from the steady flux. A
stationary observable has no signal-to-noise decay at all, so it is the natural KOA experiment.

**The current binary cannot do it.** In `edmd.h` the heat bath is a single global — `heatbath_enabled`
(int) and `heatbath_temperature` (double) — applied to every thermalising wall, and `S->heat_bath`
accumulates **one total**. With both walls at the same temperature the net steady flux is ~0 by
construction, and the per-side fluxes are not separable even in principle from one accumulator.

What the experiment needs, stated so it can be built later and **not built now**:

1. `heatbath_temperature` becomes a per-wall array, set from a comma list (the `--temperature`
   list at `00ALLINONE.c:4579` is the pattern to copy, including its per-segment parser).
2. `S->heat_bath` becomes a per-wall accumulator, so flux in at the hot wall and flux out at the
   cold wall are measured separately and their difference is the steady-state check.
3. Trace columns for each accumulator.

That is a core change to the energy ledger and is **not made here**. Level 0's bit-identity gate
and the health contract would both have to be re-run against it.

## 7. Status of Level 4

| | |
|---|---|
| (i) work splits in half | **holds**, 0.1–1.5σ, no mass dependence over ×20 in M_d |
| (iii) divider takes Δx/2 | **holds**, 0.1–1.3σ |
| stage 1 offset, KR vs ideal | **KR favoured 7–8σ** (`260928`) |
| σ(T₁−T₂) = 4T̄/(2√(2N+1)) | **confirmed**, 0.1980 vs 0.199 |
| (ii) τ_T at M_d = 10 | **measured, 490 ± 190** — consistent with 504 and with 618 |
| (ii) τ_T ∝ M | **untested** — one mass only |
| Eq. (51) drift velocity, Eq. (52) | **untested** |

The M² claim remains withdrawn. This run does not test the exponent, but its single measured value
at M_d = 10 sits on the linear law and two orders of magnitude away from what the M² fit implied at
that mass.
