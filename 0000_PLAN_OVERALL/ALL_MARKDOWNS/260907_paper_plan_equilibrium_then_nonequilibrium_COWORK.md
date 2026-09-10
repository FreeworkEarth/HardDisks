# Paper plan: equilibrium validation first, then how the piston experiments get "proved"

Date: 2026-09-07 (Cowork / Claude). Companion to `260907_speedsound_recap_and_critical_review_COWORK.md`.

Read for this: `validation/pressure_validation.c`, `run_pressure_campaign*.sh`, `test_runner_and_smoke.sh`, `EXPERIMENT_CORE_VALIDATION_AUDIT_20260820.md`, `VALIDATED_EXPERIMENT_RUNS_20260819.md`, `TIME_EDMD_DIFFERENCES.md`, `THERMAL_WALL_IMPLEMENTATION.md`, `experiments_energy_transfer/{GOALS,GOALS_ACHIEVE_EXPERIMENTS,EXPERIMENT_RUN_VALIDATION,WORK_TARGET_PISTON_STOP}.md`, `SIMPLE_BOX_WORKFLOW.md`, `ZZZ_PAPER/SIMPLE_GAS_BOX/chatgpt_discuss.md`, `README_ENTROPY_IMPLEMENTATION.md`, the two knowledge-transfer markdowns.

---

## 0. Straight answer to "did we check pressure, Z and c_s?"

**c_s: yes**, from the files on disk (routes A/B, strip ladder, N=1000, fixed-aspect). See the companion review.

**Z / pressure: no.** What I have read is the *design* and last night's debugging, not results:

- `pressure_validation.c` measures Z by two independent routes (pair virial; wall momentum flux, whose x and y components give an isotropy check — not a third method), in a square hard-wall box with no piston, box ∝ √N, N = 400/900/1600, block errors, ψ₆, and the strict health contract. Good design.
- The header comment already records the interesting physics: the wall estimators sit above the pair virial by an amount that shrinks like perimeter/area (~1/√N). That is the same hard-wall boundary effect that drives the c_s finite-size story, so the two halves of the paper share one mechanism.
- `run_pressure_campaign2.sh` says the η ≤ 0.5 grid is complete (96/96 accepted, 0 discards) and "preserved separately", and that the dense grid (0.60, 0.65, 0.67, 0.69 main; 0.702/0.710/0.720 diagnostic) is what stalled on the seeder hang. The seeder regression tests GPT listed exist as `test_lattice_seeder.c`, `lattice_smoke.c`, `test_runner_and_smoke.sh` (built 2026-09-06/07); whether they all pass and whether the dense grid has been relaunched I can't tell from here.
- The campaign writes to an OUTDIR passed on the command line, and I could not find that folder (nothing named `pressure*` in HardDisks/ or hspist3/). **Tell me the OUTDIR of the 96/96 dataset** (it holds `pressure_trajectories.csv`, `pressure_blocks.csv`, `run.log`) and I'll check the Z numbers against Kolafa–Rottner, the three-estimator agreement, the 1/√N wall term, and the block statistics the same way I did c_s.

Also on the c_s side: the `--edmd-acc=1` core FAILED validation (181/450 invalid per the chat export) and is excluded; and the Román reproduction is asserted in chat, not re-checked by me.

---

## 1. Paper 1 — equilibrium validation of the hard-disk piston code

What it has to establish, in one sentence: *this code, with this validity contract, reproduces known equilibrium properties of the 2D hard-disk fluid, and here is exactly where finite hard-wall boxes stop being bulk.*

Sections that are already backed by data or nearly so:

1. **Validity contract and what it caught.** Fail-closed ledger, forced-advance / overlap / clamp / wall-overdue counters, initialization gates, seeds per trajectory. The failure classes found on the way are a result, not an embarrassment: the deterministic seed 2381038820 (missed PP event during hold), the harmonic-wall zero-gap recollision, the event-budget forced advance that produced a plausible Z at η=0.5/N=1600, the inert lattice seed that reported Z=1 as valid, corner-packing, wall-contact seeding, the random seeder hang above η≈0.55. That list is the software-validation contribution.
2. **Z(η)** by two independent routes plus the isotropy check, N = 400/900/1600, wall term extrapolated out, vs Kolafa–Rottner on the stable fluid (η ≤ 0.69), with the 0.70–0.72 diagnostic clearly labelled exploratory. Needs the dense grid.
3. **c_s(η)** by the Román divider-resonance method, fluid branch η ≤ 0.63–0.65 (route A + strip ladder + N=1000), vs EOS-derived adiabatic c_s; route B as a slit-confinement series. Done modulo the fixes in the review.
4. **Finite size and ordering near the transition**: fixed-aspect family, ψ₆(N), the strip lesson, seeding/commensurability. Presented as a finite-size study, not as c_s values. Needs the stationarity audit and the c∞ fits from the review.
5. **Equilibrium fluctuations of the divider itself** — cheap, not yet done, and it is the bridge to Paper 2: for a free divider ½M⟨v_w²⟩ = ½kT with a Maxwellian v_w; for the spring-loaded divider the position variance is kT/k_eff with k_eff = k + 2 N m c_s²/L² (the gas stiffness from the validated EOS); the ring-down / FFT linewidth gives the divider's damping, and the force autocorrelation on the *held* divider or piston gives the Sivak–Crooks friction γ(λ) = β∫⟨δF(0)δF(t)⟩dt (not the velocity autocorrelation — its integral vanishes for a confined object). These use the same runs as c_s (the hold phase and the ring-down) and they certify the divider–gas coupling *before* anything is driven.

What I would *not* put in Paper 1: any thermostat, any driven protocol, TIME mode.

---

## 2. Paper 2 — how you "prove" the piston energy-transfer experiment

There is no EOS to compare a driven process to, so validation has to come from things that are exact regardless of the model details. There are more of those than it looks. Ordered from cheapest/strongest to the actual new physics:

### Level 0 — identities that must hold to machine precision (EDMD)

- **First law per trajectory:** W_piston = ΔKE_gas + ΔKE_walls + ΔE_spring + Q_bath, with Q_bath = 0 for elastic walls. In EDMD every term is exact between events, so the residual should be ~1e-10 relative, not "small". The audit explicitly lists a hard energy-balance gate as *not yet done* because the sign conventions of piston work, spring energy, wall KE and thermostat heat were never audited in one equation. That audit is the first task of Paper 2; until it exists, every "efficiency" number is unaudited.
- **Momentum balance** of the piston impulse against the gas.
- **Time reversal:** reverse all velocities at time τ, integrate back, recover the initial state to roundoff for short τ. Hard-disk chaos makes this degrade exponentially, so it is a short-τ test of event ordering — the class of bug that produced seed 2381038820.

### Level 1 — equilibrium properties of the moving parts (Paper 1, section 5)

Equipartition and Maxwellian divider velocity; spring-divider position variance vs EOS; friction from the ring-down and from the held-object force autocorrelation (Sivak–Crooks). These are the "before driving" baseline, and the friction coefficient is reused in Level 4: for a slow protocol ⟨W_diss⟩ ≈ ∫γ(λ)λ̇²dt, an exact linear-response statement with no fitted parameter.

### Level 2 — exact non-equilibrium results in the ideal-gas limit

Switch particle–particle collisions off (the `simple_prediction_box` already does this) or use η → 0. Then:

- **Work distribution for a piston compressing an ideal gas** is known exactly (the two papers already in `ZZZ_PAPER/SIMPLE_GAS_BOX`: Lua & Grosberg 2005, and the adiabatic-compression work distribution paper). Compare the simulated histogram of W over seeds with the analytic P(W) for a constant-speed piston at several Mach numbers. This tests the piston–particle collision rule, the work accumulator and the sampling, with no free parameter.
- **Jarzynski equality** ⟨e^{−βW}⟩ = e^{−βΔF} with the *initial* canonical ensemble (prepare with a bath, then drive with elastic walls — no thermostat needed during the drive). For the 2D ideal gas ΔF = −N k T ln(A_f/A_i). Use small N (5–20) so the exponential average converges; the identity is exact for every protocol speed, which makes it a strong test of the driven machinery. Also ⟨W⟩ ≥ ΔF and the quasi-static limit W → N k T_i (A_i/A_f − 1) (2D adiabatic ideal gas, γ = 2).
- **Same test with collisions on** at finite η: now ΔF comes from the validated EOS, βΔF/N = ln(η_f/η_i) + ∫_{η_i}^{η_f} (Z−1)/η dη. This is the direct link between the two papers: Paper 1's Z(η) is the reference the Paper 2 free energies are checked against. If the Jarzynski estimate at slow driving does not land on the EOS value, either the code or the sampling is wrong, and you know which side to look at because the ideal-gas version was already checked.

### Level 3 — thermostatted processes (heat, Crooks, the prediction box)

Everything with a heat bath depends on the thermal wall being a real heat bath. As documented in `THERMAL_WALL_IMPLEMENTATION.md` it is not one, on three counts:

1. It samples the *speed* from the 2D Maxwell speed distribution and the *direction* uniformly in the half-plane. A diffuse (Maxwell) wall re-emits the normal velocity component with p(v_n) ∝ v_n e^{−mv_n²/2kT} and the tangential one Gaussian, i.e. flux-weighted with cosine-law angles. The documented kernel re-emits particles with mean energy kT; the correct kernel (and the incoming flux) carries 3kT/2 in 2D. So the wall as documented *cools* the gas, violates detailed balance, and produces an anisotropic distribution near the wall. The doc's claim "includes correct flux weighting" is not right as written.
2. It only fires when |T_gas − T_bath| exceeds a threshold, using the *global* gas temperature. That makes it a global feedback thermostat, not a local heat bath: at T_gas ≈ T_bath the wall is elastic and no heat can flow, which is exactly when a driven process needs heat to flow.
3. The default "adaptive" mode overshoots toward the target. Fine for equilibrating, not for a measurement.

Before any Crooks / heat / prediction-box result: implement the standard diffuse wall (Tehver et al. 1998, "Thermal walls in computer simulations" is the usual reference; or the Andersen-type kernel already in `ANDERSEN_IMPLEMENTATION_00ALLINONE.md` if it is per-collision and unconditional), always on, no threshold, no overshoot, with the heat exchanged per wall collision logged. Validate it in equilibrium (T_gas → T_bath with no feedback; Maxwellian bulk; zero mean heat flux) and out of equilibrium in the free-molecular limit (two walls at T₁ ≠ T₂, no PP collisions: the heat flux is exactly computable), then repeat Level 2's Jarzynski test with the bath on, and add **Crooks**: P_F(W)/P_R(−W) = e^{β(W−ΔF)} for compression vs. expansion. Only after that does a thermostatted number mean anything.

### Level 4 — linear response: the regime where "energy transfer" has a prediction

The piston → gas column → divider(+spring) → second column chain is, for small amplitude and slow driving, a linear acoustic system whose every parameter is now known: c_s(η) and ρ from Paper 1, L, M, k, and the divider friction from Level 1. That gives the transfer function and the energy delivered to the spring for a given protocol with **no free parameters** (the Román transcendental relation is the one-column special case). Measured SpringE(t) vs. this prediction, as a function of protocol duration relative to L/c_s and to the divider's resonance period, is the "proof" of the energy-transfer experiment. Fluctuation–dissipation gives a second check: the excess work at slow driving must equal ∫γ(λ)λ̇²dt with the friction measured from equilibrium force fluctuations (Sivak & Crooks 2012).

Then the deviations at large amplitude and speed are the new physics, and they are reported *as* deviations from a validated linear baseline, not as bare efficiencies.

Two things about the current runs in this light:

- Every validated energy-transfer command drives with `--piston-speed 100` capped by 10σ of travel in 2 time units (so ≈5 thermal speeds, Mach ≈ 2 at η=0.2, higher in the dilute cases) and stops at `--piston-work-target 20000` with N=600, i.e. ≈33 kT per particle. That is a strong shock, the temperature rises ~30×, and nothing in Levels 2–4 converges there (the Jarzynski average is dominated by rare trajectories, linear response does not apply). Those runs are fine as the far end of a speed ladder; they cannot be the only regime. Add Mach 0.05 → 5 at fixed injected work.
- "Efficiency = SpringE_max / W_in" is the maximum over time of an oscillating quantity, taken inside a window (`--eff-window 20`). It is a transient, not a thermodynamic work, and it is sensitive to the window and to when the spring happens to peak. For a paper you want a quantity with a definition that survives: energy delivered to a *load* (a damped divider, or a ratchet, or the spring energy at a defined phase), or better the time-resolved SpringE(t) compared to the Level-4 prediction. The notes in `00_COMMANDS_EXPERIMENTS.md` record 80–90% in TIME mode and "99.9 or more" in EDMD for this observable in comparable one-wall setups; `TIME_EDMD_DIFFERENCES.md` explains why TIME is wrong under strong driving, which is one more reason to stay in EDMD and one more reason not to lean on that observable.

### Level 5 — Thermodynamics of Prediction (the simple prediction box)

The setup in `SIMPLE_BOX_WORKFLOW.md` is the right one: stochastic piston targets (iid vs. correlated), drive + relax, work per step, bath heat, coarse gas features s_t, and the three informations I(s_t; x_t) (memory), I(s_t; x_{t+1}) (prediction), nostalgia = memory − prediction. Still et al. relate the dissipation in the work step to kT × nostalgia; as far as I remember it is an identity given their definition of dissipated work through the nonequilibrium free energy and an instantaneous work step, and how it maps onto a finite-time drive and onto the mechanically measured W − ΔF_eq has to be worked out from the paper (it is in `ZZZ_PAPER`) before the test is designed. To make that a test rather than an illustration:

- It needs Level 3 first (the heat is half of the equation).
- The work step in the simulation is finite-time (the chat already settled: no teleporting), which adds finite-rate dissipation on top of the nostalgia term; so run a protocol-speed ladder and extrapolate to the slow limit, or compare against the inequality only.
- Mutual information from binned coarse features is biased upward with few samples and is a lower bound on the microstate quantity (data processing). So the testable statement is: measured dissipation ≥ kT × (I_mem − I_pred) computed from coarse features, with bias-corrected estimators and shuffle nulls, over thousands of protocol steps. The ideal-gas version (PP off) is where this can be pushed hardest, because the single-particle dynamics between moves is exactly solvable and the informations can be cross-checked semi-analytically.

### Level 6 — Szilard / information engine

Separate project; the same validity contract applies, but it needs protocol-aware crossing rules (the audit says so). Not part of Paper 2.

---

## 3. What this means for the order of work

1. Finish Paper 1's data: dense Z grid after the seeder tests pass; the c_s fixes from the review; the divider-fluctuation section. Write it. It is close.
2. Paper 2, in this order, each step gating the next: energy-balance audit and hard gate → ideal-gas exact work distributions and Jarzynski → finite-η Jarzynski vs EOS → correct thermal wall + its equilibrium and free-molecular validation → Crooks → linear-response prediction of energy transfer vs measurement → speed/amplitude ladder → prediction box.
3. Do not spend more trajectories on the Mach-2, 33 kT/particle energy-transfer grid until the energy-balance gate exists and the slow-driving baseline is in hand; those runs can be re-used as the far end of the ladder later.

---

## 4. Paste-able block for Claude Code (Paper 2 groundwork; no new physics runs before item 3 passes)

```
Non-equilibrium groundwork. No production energy-transfer or prediction-box runs
until items 1-3 are reported.

1. Energy accounting audit. Write down ONE equation
      W_piston = dKE_gas + dKE_walls + dE_spring + Q_bath
   with the sign convention of every accumulator in 00ALLINONE.c (piston work
   accumulators, spring energy, wall KE, thermal-wall heat). Add a per-trajectory
   residual to the energy-transfer trace and summary. For elastic walls in EDMD
   the relative residual must be at roundoff level; report its distribution over
   the existing validated_reference_20260820_all_gates_r1 runs before making it
   a rejection gate.

2. Thermal wall audit. Quote the actual re-emission kernel from the code (speed
   distribution, angular distribution, activation threshold, adaptive overshoot).
   Compare with the diffuse Maxwell wall: normal component p(v_n) ~ v_n exp(-m v_n^2/2kT),
   tangential Gaussian. If the code samples the 2D speed distribution with uniform
   angles, or gates on a global temperature threshold, it is not a heat bath; do
   not fix it silently, report first. Then implement/select a per-collision,
   unconditional diffuse wall, log heat per wall collision, and validate:
   (a) T_gas -> T_bath with no feedback, (b) Maxwellian bulk velocities,
   (c) zero mean heat flux at equilibrium, (d) two walls at T1 != T2 with PP
   collisions off: heat flux vs the free-molecular prediction.

3. Ideal-gas exact tests (PP collisions off, elastic walls during the drive,
   canonical initial ensemble prepared by the validated thermal wall):
   (a) constant-speed piston compression, N = 5, 10, 20, three piston speeds:
       histogram of W over >= 2000 seeds vs the analytic work distribution
       (papers in ZZZ_PAPER/SIMPLE_GAS_BOX);
   (b) Jarzynski: <exp(-W/kT)> vs exp(-dF/kT) with dF = -N kT ln(A_f/A_i);
       report the estimator with bootstrap error and the sample-size dependence;
   (c) slow limit: <W> -> N k T_i (A_i/A_f - 1).

4. Finite-eta Jarzynski at eta = 0.2 and 0.5, N = 50: dF from the EOS,
   beta dF/N = ln(eta_f/eta_i) + integral (Z-1)/eta d eta, using Kolafa-Rottner and
   the campaign's own Z if available. Speed ladder Mach 0.05 .. 2.

5. Divider equilibrium fluctuations from EXISTING speed-of-sound holds and ring-downs:
   <v_w^2> vs kT/M, velocity histogram, position variance of the spring divider vs
   kT/k_eff with k_eff = k + 2*N_side*m*c_s^2/L^2, damping from the ring-down /
   FFT linewidth, and the held-divider force autocorrelation
   gamma = beta * integral <dF(0) dF(t)> dt (Sivak-Crooks). Do NOT use the velocity
   autocorrelation integral for a confined divider (it is zero).

6. Linear-response baseline for energy transfer: derive the small-amplitude
   transfer function of piston -> column -> divider(+spring) from c_s(eta), rho, L,
   M, k; compare SpringE(t) from small-amplitude slow runs (Mach <= 0.1, injected
   work <= 0.1 kT per particle) against it. Only then rerun the existing Mach-2
   grid as the far end of a speed ladder.

Do not use TIME mode for any of this. Do not auto-commit.
```

---

## 5. Questions for GPT / Susanne

- Is the intended Paper 2 quantity the Still et al. nostalgia–dissipation relation (prediction box), the energy-transfer efficiency of the coupled chain, or both? The validation ladder is the same up to Level 3; after that the two need different experiments.
- For the energy-transfer chain: what is the *load*? A spring stores and returns energy; a paper needs a definition of delivered work that does not depend on a window.
- Is there a reason to keep the adaptive/thresholded thermal wall at all, or can it be replaced outright by the diffuse Maxwell wall?
