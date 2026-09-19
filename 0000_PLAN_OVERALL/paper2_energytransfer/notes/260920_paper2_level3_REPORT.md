# Paper 2, Level 3 pilot — one gas driving a spring-loaded wall

2026-09-19. Geometry C: rigid wall | spring (k = 5) | free wall (M_s = 200) | gas (N_s = 50,
η = 0.1013) | piston. 125 runs, five speeds × 25 seeds, Δx = 3.93 σ of gas.
**0 aborts, 0 `[EDMD-HEALTH]` lines.** Script `validation/paper2_level3_20260919.py`,
batch `experiments_energy_transfer/_run_scripts/level3_spring_geometryC.sh`,
figure `experiments/final/260919_level3_spring.png`.

**Verdict: the pilot ran and the ledger closes, but the experiment as designed cannot measure what
Level 3 asks. It is noise-limited by construction. Level 3 is NOT passed; the design is the result.**

## Design numbers, computed before the batch

```
eta = 0.10134, Z = 1.2399, c_s = 1.7493, k_gas = N_s kT (Z + eta Z')/L^2 = 0.05071
k/(k + k_gas) = 0.9900          omega_w = 0.15891      T_w = 39.54 sigma-time
2L/c_s = 44.30 sigma-time       pushes: 786, 196, 79, 39, 20 sigma-time
```

The ladder was chosen to cross the wall period: τ_push/T_w runs 19.9 → 0.5, so it spans the
quasi-static limit and the impulsive one.

## 0. Gate — the ledger with the spring term

Closes to the precision of the trace: R_E = W − [ΔKE_gas + ΔKE_wall + E_spring] = **+1.8 × 10⁻⁵**
against W = 7.15, i.e. |R_E|/W = 2.5 × 10⁻⁶. The partition table below sums to ⟨W_in⟩ in every row
to three decimals.

**The 10⁻⁸ criterion is not met and cannot be, this way.** 2.5 × 10⁻⁶ is the six-decimal printing of
the trace columns, not a physics residual — the Level 0 ledgers reach 10⁻¹² only because they use the
event log, which carries no spring term. Either the trace needs more decimals or the spring energy
needs to enter the event log. Recorded as a gap, not waved through.

## 1–2. What came out

| u | τ_push [σ] | τ/T_w | ⟨W_in⟩ [kT] | ⟨E_spring,max⟩ [kT] | t_max after stop [σ] | ε = ⟨E⟩/⟨W⟩ |
|---|---|---|---|---|---|---|
| 0.005 | 786 | 19.9 | 7.214 ± 0.015 | 2.249 ± 0.177 | 74 ± 13 | 0.312 ± 0.025 |
| 0.02 | 196 | 5.0 | 7.284 ± 0.038 | 3.547 ± 0.310 | 122 ± 21 | **0.487 ± 0.043** |
| 0.05 | 79 | 2.0 | 7.323 ± 0.084 | 3.440 ± 0.248 | 153 ± 12 | 0.470 ± 0.034 |
| 0.10 | 39 | 1.0 | 7.334 ± 0.146 | 2.857 ± 0.303 | 164 ± 15 | 0.390 ± 0.041 |
| 0.20 | 20 | 0.5 | 8.283 ± 0.342 | 2.795 ± 0.244 | 165 ± 16 | 0.337 ± 0.029 |

The secondary estimator ⟨E/W⟩ agrees with ε to the third decimal everywhere.

Energy partition at the end of the record — it sums to the work exactly, which is the ledger again:

| u | ΔKE_gas | KE_wall | E_spring | sum | ⟨W_in⟩ |
|---|---|---|---|---|---|
| 0.005 | 6.303 | 0.448 | 0.463 | 7.214 | 7.214 |
| 0.02 | 6.063 | 0.686 | 0.535 | 7.284 | 7.284 |
| 0.05 | 6.195 | 0.832 | 0.296 | 7.323 | 7.323 |
| 0.10 | 6.365 | 0.499 | 0.469 | 7.334 | 7.334 |
| 0.20 | 7.056 | 0.583 | 0.644 | 8.283 | 8.283 |

**About 86 % of the injected work is in the gas at the end, and the spring holds 4–8 %.**

## 3. Why these ε values are not a measurement of capture

The quasi-static prediction is
x_eq = ΔP·H/(k + k_gas) = 0.044 σ and E_qs = ½k x_eq² = **0.0051 kT** — four hundred times smaller
than the ⟨E_spring,max⟩ in the table. That gap is the finding.

The spring is a quadratic mode in contact with a gas at kT = 1, so it carries **½kT = 0.5 kT of
thermal energy whatever the piston does**, with an rms amplitude √(kT/k) = 0.447 σ — ten times the
deterministic displacement. Separating signal from noise:

| u | peak of the seed-**averaged** curve (coherent) | mean of per-seed peaks (what the table reports) |
|---|---|---|
| 0.005 | 0.636 | 2.249 |
| 0.02 | 1.017 | 3.547 |
| 0.05 | 1.315 | 3.440 |
| 0.10 | 0.815 | 2.857 |
| 0.20 | 0.836 | 2.795 |

The per-seed maximum is ~3× the coherent peak, because taking a maximum over a long record of a
fluctuating quantity samples the tail of the thermal distribution. **ε = 0.31–0.49 is therefore an
overestimate; the coherent capture is 0.09–0.18**, and even that sits on a ½kT pedestal.

## 4. What to change before Level 3 is run for real

The deterministic energy a spring can take is E = F²/2k with F = ΔP·H = 0.22 kT/σ, so a *stiffer*
spring captures less, and a softer one needs a displacement the box does not have. At this
compression the capturable energy is at most a few tenths of kT against a ½kT thermal pedestal.
Level 3 needs a bigger drive, not a different spring:

1. **Compress harder.** ΔP ∝ the compression; E ∝ ΔP², so doubling Δx is a factor four.
2. **Denser gas.** At η = 0.5 the pressure is 4× larger for the same N_s.
3. **Report the coherent (seed-averaged) spring energy as primary**, with the per-seed maximum as a
   clearly-labelled secondary — the thermal pedestal must be stated either way.
4. **Subtract the pedestal explicitly**: quote E_spring,max − ½kT, or measure the same box with the
   piston never released and subtract that distribution.

## Open

- The 10⁻⁸ ledger gate needs the spring energy in the event log (or more trace decimals).
- Ring-down time and Mansour's Q for this wall: not fitted, because with the coherent signal at the
  thermal level the decay constant is not identifiable from 25 seeds.
- Both limits (quasi-static and impulsive) remain untested for the same reason.
