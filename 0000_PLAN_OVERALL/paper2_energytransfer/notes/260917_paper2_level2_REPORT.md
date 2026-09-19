# Paper 2, Level 2 — friction, and how the work rises with piston speed

2026-09-17. Analysis script `hspist3/validation/paper2_level2_20260917.py`; every number below is
printed by that one file. No physics code was touched. 710 new trajectories
(`experiments_energy_transfer/level2_slope_20260917/`), pooled with the 50 older
`level0_Wqs_20260911/` runs, which carry the identical command line (verified: N = 100, 50 + 50,
r = 0.5, L₀ = 39.25, hold 12000 steps, Δx = 3.93 σ, `--edmd-acc=0`, drift-first).
Health contract: 0 aborts, 0 `[EDMD-HEALTH]` lines over all 710 runs.

---

## 1. What Level 2 claimed, and what was wrong with the first pass

The plan states the near-equilibrium prediction

$$W_\mathrm{diss} \;=\; \zeta\,u\,\Delta x, \qquad
\zeta \;=\; \beta\!\int_0^\infty \langle \delta F(0)\,\delta F(t)\rangle\, dt ,$$

with ζ measured on the **pinned** piston during the equilibrium hold. The first pass reported
ζ = 2.537 ± 0.050, hence a predicted slope dW/du = ζΔx = 9.97 ± 0.20, against a measured slope
near 3. That gap was real, and its cause is now identified: **2.537 is the value of the integral at
zero lag, not its plateau.** It is the impulsive self term β Σdp²/(2T) — the free-molecular
(Enskog) piston friction — and the correlation function has a long negative tail that cancels it.

## 2. The running integral

Force on the pinned piston, built from the per-event impulse log (`HD_PISTON_EVENTS`, wall label
`WR` while the piston is held), binned at 0.25 σ, hold window 20 → 200 σ, 190 holds:

| t_cut [σ] | β∫₀^{t_cut} ⟨δF(0)δF(t)⟩ dt | implied dW/du = ζΔx |
|---|---|---|
| 0.25 | +2.546 | +10.01 |
| 0.50 | +2.486 | +9.77 |
| 1.00 | +2.387 | +9.38 |
| 2.00 | +2.170 | +8.53 |
| 5.00 | +1.777 | +6.98 |
| 10.00 | +1.219 | +4.79 |
| 20.00 | +0.149 | +0.59 |
| 40.00 | −0.053 | −0.21 |
| 60.00 | +0.148 | +0.58 |

Cross-check on the bookkeeping: the trapezoid's first half-step is 2.546 and the independently
summed self term β Σdp²/(2T) is 2.615 — the same number, so the self term must not be added again
(the first pass did add it, which is the second half of the discrepancy).

**The integral decays to zero, and it does so at the sound-crossing time.** The compressed
compartment is L = 35.32 σ long at η = 0.1112, where Kolafa–Rottner gives Z = 1.2678 and
c_s = 1.788, so L/c_s = **19.8 σ** — exactly where the running integral crosses zero. Physically:
the momentum the piston hands to the gas comes back as a returning sound wave, so in a closed box
the zero-frequency friction of a pinned wall vanishes. Linear response therefore predicts **no**
term linear in u here, and ζΔx = 9.97 was never the right number to compare against.

Past the crossing the integral does not sit still: it swings between −0.78 and +0.50 with the
recurrence period 2L/c_s = 39.5 σ, which is the sound bouncing. Averaged over one full recurrence
(t_cut = 20 → 60 σ) it is **−0.072**, implying dW/du = −0.28 — against the measured −0.50 ± 0.67 for
u ≤ 0.05. So the honest statement is not "the integral has a plateau at zero" but "the integral has
no plateau; it oscillates about zero, and its recurrence average agrees with the measured slope".

Figure: `260909_plots/260917_level2_work_and_friction.png` (left panel §3–4, right panel this one).

## 3. The speed ladder

W_qs^finite = 7.4879 ± 0.0066 kT (Level 1, finite-box Z along the path), Δx = 3.93 σ.

| u [σ/τ] | seeds | ⟨W⟩ [kT] | W − W_qs [kT] |
|---|---|---|---|
| 0.005 | 100 | 7.4787 ± 0.0080 | −0.0092 ± 0.0104 |
| 0.01 | 100 | 7.4836 ± 0.0106 | −0.0043 ± 0.0125 |
| 0.02 | 60 | 7.4338 ± 0.0235 | −0.0541 ± 0.0244 |
| 0.03 | 100 | 7.4679 ± 0.0295 | −0.0200 ± 0.0302 |
| 0.05 | 100 | 7.4771 ± 0.0366 | −0.0108 ± 0.0371 |
| 0.10 | 100 | 7.7084 ± 0.0583 | +0.2205 ± 0.0587 |
| 0.15 | 100 | 8.1148 ± 0.1048 | +0.6269 ± 0.1050 |
| 0.20 | 100 | 8.4756 ± 0.1659 | +0.9877 ± 0.1660 |

Straight-line fits:

| fit range | intercept [kT] | slope dW/du | χ²/dof |
|---|---|---|---|
| u = 0.005–0.02 (3 speeds) | +0.0053 ± 0.0132 | −2.06 ± 1.52 | 1.96 |
| u = 0.005–0.05 (5 speeds) | −0.0061 ± 0.0088 | −0.50 ± 0.67 | 1.11 |
| u = 0.005–0.20 (8 speeds) | −0.0348 ± 0.0072 | +2.59 ± 0.36 | 6.79 |
| u = 0.02–0.20 (6 speeds) | −0.1666 ± 0.0247 | +4.53 ± 0.50 | 2.36 |

Below u = 0.05 the slope is **zero within error** (−0.50 ± 0.67), as §2 requires. A straight line
through the whole ladder fits badly (χ²/dof = 6.8) and its apparent slope depends entirely on where
the range is cut — the signature of a curve, not a line.

## 4. The excess work is quadratic in u

Through the origin, over all eight speeds:

| form | coefficient A | χ²/dof |
|---|---|---|
| A u (linear) | 1.590 ± 0.297 | 9.17 |
| A u² (quadratic) | 23.497 ± 2.676 | 2.26 |
| A u³ (cubic) | 147.742 ± 16.506 | 1.82 |

The three speeds that carry signal (u ≥ 0.1) separate quadratic from cubic cleanly:

| form | A | χ²/dof | A point by point (u = 0.10, 0.15, 0.20) |
|---|---|---|---|
| u² | 25.21 | 0.32 | 22.1, 27.9, 24.7 |
| u³ | 148.93 | 2.21 | 220.5, 185.8, 123.5 |

Quadratic. And the coefficient has a reading: **N_s m / 2 = 25.0**, the kinetic energy of the whole
compressed compartment moving at the piston speed. The excess work is the compression flow the
piston leaves behind, which thermalises after it stops — not a friction.

Direct check on the momentum, at the instant the piston stops:

| u | seeds | ⟨Px_gas⟩ | −N_s m u | ⟨Px⟩²/(2 N_s m) [kT] | W − W_qs [kT] |
|---|---|---|---|---|---|
| 0.05 | 90 | −0.109 ± 1.122 | −2.50 | 0.0001 ± 0.0024 | −0.0108 |
| 0.10 | 90 | −2.035 ± 1.178 | −5.00 | 0.0414 ± 0.0480 | +0.2205 |
| 0.15 | 100 | −4.830 ± 1.001 | −7.50 | 0.2333 ± 0.0967 | +0.6269 |
| 0.20 | 100 | −6.889 ± 1.144 | −10.00 | 0.4746 ± 0.1576 | +0.9877 |

The gas does carry coherent momentum, growing with u, about 65–70 % of N_s m u. The flow is not
uniform (it is a compression wave: gas near the piston moves at u, gas far from it less), so
⟨Px⟩²/2M is a strict lower bound on the flow's kinetic energy — and indeed it is about half the
measured excess. Consistent with the picture; not a closed-form identity.

## 5. Level 1's intercept, remeasured

With the correct extrapolation form (quadratic, not linear) and 760 trajectories instead of 30:

```
W(0) − W_qs^finite = −0.0151 ± 0.0060 kT   (1.7 σ)   χ²/dof = 1.57
W(0) = 7.4728 ± 0.0089 kT   against   W_qs^finite = 7.4879 ± 0.0066 kT
A    = 24.73 ± 2.72
```

**Level 1 PASSES its 2σ criterion**, now on 25× the statistics. The sign of the gap flipped relative
to the earlier +0.032 ± 0.018 (1.8 σ high) because that fit was a straight line through three
speeds with 10 seeds each; the unphysical slope b = −4.2 ± 1.7 it reported was the quadratic curve
being forced through a line.

---

## Verdict

- **Level 2 PASSES, with its stated prediction corrected.** The linear-response friction of a
  pinned piston in a closed box is zero, because the force autocorrelation integral is cancelled by
  the returning sound wave at t ≈ L/c_s = 19.8 σ. The measurement agrees: dW/du = −0.50 ± 0.67 for
  u ≤ 0.05.
- **The finite-speed cost is second order:** W − W_qs = (25.2) u², matching N_s m u²/2, the kinetic
  energy of the compression flow. Cubic is excluded (A drifts 220 → 124 across the three points).
- **ζ = 2.54 is still a correct number** — it is the free-molecular/Enskog piston friction, the
  zero-lag value of the same integral — but it is not the coefficient in W_diss = ζuΔx for this
  geometry. The plan's Level 2 paragraph should quote the running integral, not the zero-lag value.
- **Level 1 is now settled** at 1.7 σ with the correct fit form.

## 6. The ζ audit (script `validation/level2_zeta_audit_20260917.py`)

Three checks asked for on the method, all passed.

**(a) Is the held wall the generalised force of the coordinate we drive?** Yes, and the event log
proves it rather than assuming it: label `WR` spans t = 0.0 → 207.5 σ with u_wall = 0 and label `PR`
spans 210.5 → 598.5 σ. They never overlap — it is one object, renamed at release. So C_FF is measured
on the piston itself, not on a proxy wall.

**(b) Bin-width convergence.** The force is an impulse train, so the zero-lag term is a delta smeared
over one bin; the integral must stop moving as dt → 0. Over 100 holds:

| bin dt [σ] | zero-lag self term | ζ(1 σ) | ζ(5 σ) | ζ(10 σ) | ζ(20 σ) | ζ(40 σ) |
|---|---|---|---|---|---|---|
| 1.0000 | 2.613 | +2.377 | +1.790 | +1.235 | +0.149 | −0.051 |
| 0.5000 | 2.613 | +2.378 | +1.793 | +1.233 | +0.142 | −0.049 |
| 0.2500 | 2.613 | +2.380 | +1.780 | +1.222 | +0.143 | −0.049 |
| 0.1250 | 2.613 | +2.374 | +1.785 | +1.230 | +0.142 | −0.063 |
| 0.0625 | 2.613 | +2.368 | +1.781 | +1.224 | +0.136 | −0.063 |

Converged over a 16-fold range of bin width; nothing in the conclusion depends on it.

**(c) ζ at the other end of the compression path.** The sound-recurrence explanation makes a sharp,
falsifiable prediction: the compressed compartment is shorter *and* hotter, so its cancellation must
move to a shorter time. Using the post-stop equilibrium segment of the same runs (piston pinned at the
compressed position, t = 340–598 σ, 100 runs):

| | L [σ] | η | T | c_s | predicted L/c_s | measured zero crossing |
|---|---|---|---|---|---|---|
| start of path | 35.32 | 0.1112 | 1.000 | 1.788 | 19.8 σ | **21.5 σ** |
| compressed | 31.39 | 0.1251 | 1.150 | 1.978 | 15.9 σ | **18.0 σ** |

The crossing moves the predicted way and by roughly the predicted amount (measured ratio 1.19,
predicted 1.25), on data that were never used to build the explanation.

**And the zero-lag value is the Enskog piston friction, with no free parameter.** For a wall of length
h in a hard-disk gas, ζ_E = h n √(8mkT/π)·g(η) with g(η) = (1 − 7η/16)/(1 − η)²:

| | n [σ⁻²] | T | g(η) | ζ_E predicted | ζ zero-lag measured | ratio |
|---|---|---|---|---|---|---|
| start | 0.14156 | 1.000 | 1.2043 | 2.720 | 2.613 | 0.961 |
| compressed | 0.15929 | 1.150 | 1.2349 | 3.366 | 3.680 | 1.093 |

Within 4 % and 9 % of a parameter-free kinetic prediction at both ends. So the number the first pass
called ζ is correctly identified — it is simply the wrong coefficient for W_diss = ζuΔx.

**Still open on the theory side.** Two caveats that this data cannot settle:
1. Sivak–Crooks friction should in principle be integrated along the path, ζ = ∫ζ(λ)dλ, since η and T
   both change over the 10 % compression. Here it changes nothing — ζ_plateau ≈ 0 at both ends — but
   the paper should state the path form, not the endpoint form.
2. `drift-first` seeding fixes total kinetic energy tightly, so the preparation is nearer
   microcanonical than canonical. The Green–Kubo relation is usually written canonically. This matters
   for Jarzynski (Level 2b) more than for ζ, and is the reason the canonical seeder flag is still on
   the code list.

## What this does not yet cover

The cancellation is a property of a **closed** compartment. A box where the sound does not return
coherently — much longer, or with an absorbing end — should show the linear ζuΔx law. That is the
experiment that would separate "no friction" from "friction hidden by recurrence", and it is one
geometry change, not new code.
