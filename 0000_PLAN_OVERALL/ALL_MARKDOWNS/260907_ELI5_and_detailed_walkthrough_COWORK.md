# The whole thing, twice: ELI5 first, then the detailed step-by-step version

Date: 2026-09-07 (Cowork / Claude). This folds together the two earlier notes (`260907_speedsound_recap_and_critical_review_COWORK.md`, `260907_paper_plan_equilibrium_then_nonequilibrium_COWORK.md`) and goes deeper. Every section has a plain-language part ("ELI5") and a technical part ("In detail"). Numbers come from files on disk unless marked otherwise.

Contents

- Part A — the machine and the three measurements (what a hard disk is, EDMD, the divider, c_s, Z, ψ₆)
- Part B — what happened, in order, from Aug 19 to Sept 7
- Part C — the eight findings from the review, each with evidence and the fix
- Part D — pressure / Z: what exists, what I could not check
- Part E — Paper 1, section by section, with acceptance criteria
- Part F — Paper 2: the validation ladder for the driven piston, level by level
- Part G — the next two weeks as a checklist
- Part H — glossary and formulas in one place

---

# Part A — the machine and the three measurements

## A1. What is being simulated

**ELI5.** Imagine coins on a frictionless table that never slow down. They fly in straight lines and bounce off each other and off the walls perfectly, like ideal billiard balls, no energy lost, ever. That is a "hard-disk gas". It is the simplest system that behaves like a real fluid: at low density it is a gas, at high density it gets stiff like a liquid, and at very high density the coins lock into a honeycomb and it is a crystal. Because nothing is ever lost, it is the cleanest possible test bed for thermodynamics.

The simulator (`00ALLINONE`) can put a movable wall — "the divider" — across the table, with coins on both sides, and it can also push one end wall in like a piston. The divider can be free, or held in place, or attached to a spring.

**In detail.** N disks of radius r (diameter σ = 2r), mass m, in a box of width L₀ (per compartment) and height H, hard walls. Units: σ = 1 when r = 0.5, m = 1, k_BT = 1 (`--kbt1`), so velocities are in units of √(k_BT/m) and time in σ/√(k_BT/m). Packing fraction η = N_compartment · π r² / (L₀ · H). The `--particles-boxes=50,50` flag puts 50 disks in each of two compartments; "N=100" in the plots is the total.

Two solvers exist. **EDMD** (event-driven molecular dynamics, `--mode=edmd`) computes the exact time of the next collision anywhere in the system, jumps to it, applies the exact elastic collision rule, and repeats. Between events particles fly freely, so the trajectory is exact up to floating point. **TIME** mode (`--mode=time`) steps a fixed dt and repairs overlaps afterwards; under strong driving it invents structure (coherent "sheets" of particles) — `TIME_EDMD_DIFFERENCES.md` explains why. All scientific runs use EDMD with the default core (`--edmd-acc=0`); the cell-list core (`--edmd-acc=1`) failed validation (181 of 450 trajectories invalid — boundary escapes and overlaps) and is excluded.

The divider is a rigid wall of thickness 0.05 and mass M = (wall mass factor) × m, free to move in x after `--wall-hold-steps` steps of being held fixed. A spring-loaded divider follows x_w(t) = x_eq + (x₀ − x_eq) cos ωt + (v₀/ω) sin ωt between collisions (ω² = k/M), and the core finds particle–wall contacts on that curve analytically.

Every trajectory carries a **validity ledger**: counters for forced advances (the core gave up on exact event processing and pushed time forward), overlap repairs, clamp repairs, overdue wall events, plus geometry checks at initialization and per step (no overlaps, no wall penetration, particle stays in its compartment, walls stay ordered, finite state). A trajectory with a violation is written to a failures CSV and never enters an analysis. This is the "fail-closed" contract from `EXPERIMENT_CORE_VALIDATION_AUDIT_20260820.md`.

## A2. The speed of sound, measured Román's way

**ELI5.** Hold the divider still while the gas on both sides settles down. Let go. The gas on each side is a spring made of air: push the divider left and the left gas pushes back harder, the right gas less, so the divider swings back and forth like a mass on a spring. How fast it swings depends on how stiff the gas is, and gas stiffness *is* the speed of sound (sound is just a stiffness wave). So: let go, listen to the note, and the pitch tells you the speed of sound. Do it with several divider masses (heavier divider, lower note) and fit — that is more robust than one note.

**In detail.** Román, White & Velasco (2002) derived the resonance of a divider of mass M between two columns of length L (measured from the divider face to the far wall) of a fluid with sound speed c and total mass N m per column:

    ν = (c_s / 2π L_eff) · K,   with   cot K = (M / 2 N m) · K,

where N is the number of disks *per side* and K = k L is the dimensionless wavenumber of the standing wave in each column. For a heavy divider K is small and this reduces to the Rüchardt result ω² = 2 N m c²/(M L²); for a light divider K → π/2 (a quarter-wave pipe). The relation is exact for a 1D linear acoustic column at any mass ratio, so light dividers are fine *for the formula*; what is not obviously fine is applying a 1D column model to a compartment that is 5 diameters long or partly crystalline (see C4).

The pipeline: hold (2000 steps × dt 0.4 ≈ 800 time units, worth confirming what one hold step is), release, record x_w(t) for ~25 predicted periods, FFT the divider velocity, take the peak frequency, repeat over 9 masses × 25 seeds, fit ν against K/(2π L_eff) per η: the slope is c_s. L_eff = L₀ − 2r (the disk centres cannot reach the walls). Flags like `peak_on_search_boundary` mark cases where the FFT peak sat at the edge of the search window — those points are not trustworthy.

The theoretical reference: for hard disks the internal energy is purely kinetic (U = N k T in 2D), P = ρ k T Z(η), so the adiabatic sound speed follows from the equation of state alone:

    c_s² = (k_B T / m) · [ Z + η Z′(η) + Z² ].

At η → 0, Z → 1 and c_s → √(2 k_B T/m) (γ = 2 for a 2D ideal gas). Z(η) comes from an EOS: Henderson (1 + η²/8)/(1−η)², Kolafa–Rottner 2006 (most accurate on the fluid branch), SPT, or Liu 2021 (global, includes hexatic/solid branches). These differ by ~1–1.5% in c_s at η = 0.6, which is the size of the agreement being claimed, so always say which one you compare to.

## A3. The compressibility factor Z, measured two independent ways (plus an isotropy check)

**ELI5.** Z is "how much harder than an ideal gas does this gas push at the same density and temperature". Z = 1 is ideal. Coins push harder because they take up space: at η = 0.5 they push about four times harder. The pressure runner (`validation/pressure_validation.c`) measures the push two unrelated ways — from how hard the coins hit each other (pair virial) and from how hard they hit the walls (momentum flux) — and checks that the left/right walls and the top/bottom walls feel the same push (isotropy). If the two routes agree in the large-box limit, that is evidence rather than an accounting identity.

**In detail.** Square box of side ∝ √N, no divider, no piston, elastic walls. Estimators:

- pair virial: Z_pair = 1 + (m σ Σ_collisions |Δv_n|) / (2 · KE · t), where Δv_n is the normal relative velocity change per pair collision and KE = N k T is the total kinetic energy; this is the 2D virial theorem for hard disks;
- x-wall momentum flux: Z_wall,x = (I_L + I_R) · W / (2 · KE · t) with I the total impulse delivered to the left/right walls and W the box width;
- y-wall momentum flux: same with top/bottom and H — the x-vs-y comparison is an isotropy check of the wall route, not a third method.

Blocks of length `block_dt` give a standard error; ψ₆ is recorded per block; the health counters are checked at the end and any nonzero count discards the trajectory. The runner's own header notes that the wall estimators sit above the pair virial by an amount that shrinks like perimeter/area (~1/√N): near a hard wall the disks are layered, the local pressure on the wall is not the bulk pressure, and that is a genuine finite-size effect that has to be extrapolated away. The same physics is behind the c_s finite-size story.

## A4. ψ₆, the honeycomb meter

**ELI5.** For each coin, look at its neighbours and ask "are they arranged at 60° steps like a honeycomb?" Score 1 if perfectly, 0 if not. ψ₆ local is the average of those scores (how honeycomb-like each neighbourhood is). ψ₆ global first averages the honeycomb *directions* and then takes the magnitude, so it is only high if the whole box is one honeycomb pointing the same way. Several small honeycombs at different angles give high local, low global.

**In detail.** ψ₆,i = (1/n_i) Σ_{j∈nbrs(i)} exp(6 i θ_ij) with neighbours within 1.4 σ; global = |⟨ψ₆,i⟩_i|, local = ⟨|ψ₆,i|⟩_i. A perfect hexagonal crystal gives 1 and 1. A perfect *square* lattice gives 0 (four neighbours at 90°: the six-fold phases cancel) — that matters in C3. The speed-of-sound runs store ψ₆ at the end of the hold and at the end of the run, plus the mean neighbour count.

## A5. The geometry arithmetic you need for everything below

With `--particles=100 --particles-boxes=50,50 --height=10`:

- route A (r = 0.5, L₀ chosen per η): η = 50 · π · 0.25 / (L₀ · 10) → L₀ = 3.927/η. At η = 0.72, L₀ = 5.45; the compartment is 5.45 σ × 10 σ, about 5 columns × 10 rows of disks.
- route B (L₀ = 20, r chosen per η): r = √(η · 200 / (50 π)) = √(1.273 η). At η = 0.72, r = 0.957, σ = 1.915; the compartment is 20/1.915 = 10.4 σ × 10/1.915 = 5.2 σ, about 10 columns × 5 rows.
- strip ladder (r = 0.5, H = 10, L₀ ∝ N): at η = 0.72, L₀ = 5.45 / 10.9 / 21.8 / 54.5 for N = 100/200/400/1000 — the compartment gets longer, never taller.
- fixed-aspect family A (r = 0.5, L₀ and H ∝ √N): H = 20/30/40 for N = 400/900/1600, L₀ ≈ 11/16.5/22 at η ≈ 0.7.

---

# Part B — what happened, in order

**ELI5.** First the team made the simulator refuse to lie (Aug 19–20): every run now checks itself and throws itself away if anything is off. Then they found and fixed a few real bugs that the self-checks exposed (Aug 20–23). Then they measured the speed of sound two different ways (Aug 23–25), found the two ways disagree at high density, tried making the box bigger in one direction (Aug 25–26), then in both directions (Aug 26), and learned the box shape matters more than the box size. Last night (Sept 6–7) they started the pressure measurement, found the random starting-arrangement code hangs at high density, wrote a lattice starter, found four bugs in it, and stopped to write tests.

**In detail, step by step.**

1. **Aug 19 — validated references.** A 288-trajectory speed-of-sound reference (8 lengths × 12 masses × 3 repeats, 288/288 valid) and a 20-trajectory energy-transfer functional reference. `VALIDATED_EXPERIMENT_RUNS_20260819.md`.
2. **Aug 20 — the core audit.** `EXPERIMENT_CORE_VALIDATION_AUDIT_20260820.md`. Changes that were not cosmetic: the energy-transfer EDMD path had particle–particle collisions *disabled* (zero-initialized parameter struct) — every legacy EDMD energy dataset is suspect; the harmonic spring-wall contact search used a coarse time scan and could skip a collision — replaced by an analytic-extrema + bisection root finder; a zero-gap recollision edge case (particle touching the wall right after a collision, spring later pushes the wall back into it) was caught by the validator at seed 2826082000 and fixed; forced-advance counter exposed and made a hard failure; `-ffast-math` removed from the release build because it lets the compiler assume NaN never happens. Explicitly *not* done: a first-law energy-balance gate.
3. **Aug 21 — high-η initialization.** The long eta-split pilot (r10) failed 120 of 180 high-η trajectories at L₀ = 5.61 and 6.04, all with the same overlap at step 0: the rectangular initializer placed rows microscopically closer than a diameter and nobody checked. Fix: verify spacing, else fall back to a checked hex lattice (Δx = d(1+ε), Δy = (√3/2)Δx). `01_improvements_bugsfxed_dev/26_08_21_HIGH_ETA_INITIALIZATION_HEX_FALLBACK.md`.
4. **Aug 22 — the deterministic seed.** Seed 2381038820 (L₀ = 150, M = 50, repeat 3) failed identically in the r10 and r11 pilots: particles 73/74 overlapping by 0.0926 at t = 2.467, during the hold, i.e. a missed particle–particle event in the core. The per-trajectory seed is a hash of the *indices* (l, m, r) in the sweep, which is why a reduced 4-run reproduction didn't hit it. `20260825_knowledge_transfer_Codex_chat_to_claude.md`.
5. **Aug 23 — core fix and route A.** The later notes refer to a "2026-08-23 tunnelling fix" (`wall_time_from_gap()` treating a zero gap as an overdue wall collision); whether that is also what resolved seed 2381038820 I could not confirm from the files I read — the chat only says "the core is fixed and verified across 14,400 trajectories". Then the full route A campaign: 32 η × 9 masses × 25 repeats = 7200 trajectories, 0 invalid, ψ₆ recorded (`campaign_r25_psi6_20260823`).
6. **Aug 25 — route B and the strip ladder.** Route B: same 32 η by sweeping the radius at fixed L₀ = 20, H = 10; 7200 trajectories, 0 invalid. The grouped analysis collapsed (groups by L₀, which is constant in B, and used r = 0.5 for L_eff) — deleted and redone per η. A vs B diverge above η ≈ 0.55. Timing: N = 1000 costs 391 s per trajectory (O(N²) core), so a ladder N = 100/200/400 at 5 η with H fixed was run instead.
7. **Aug 26, 00:26 — accelerated core check.** `validate_acc_N100_20260826`: route-A geometry at 5 η with `--edmd-acc=1`. The chat export settles it: **269 valid / 181 invalid of 450** (90 boundary escapes, 91 particle overlaps) on settings where the default core gives 450/450 — the cell-list core is unusable and excluded; the overnight N=1000 ran on the default core at ~405 s per trajectory. (My earlier note that it "ran clean 90/90" looked at one η's status file only — corrected.)
8. **Aug 26, 02:39–11:25 — overnight N = 1000.** 4 η × 9 masses × 18 repeats = 648 trajectories at H = 10, one cell with `health=1` (η = 0.700, M = 750) still counted valid, script died after merging with a shell syntax error; analysis produced separately. Plots include a 1/√N extrapolation (c∞ = 10.44 at η = 0.70).
9. **Aug 26, 14:25–20:32+ — fixed-aspect campaign.** Family A (box ∝ √N, aspect fixed, N = 400/900/1600, η = 0.67–0.73) and family C (fixed box, r ∝ N^−1/2, η = 0.710 control); 2540 trajectories, 0 invalid. Conclusion in the figure title: "the strip result does not survive". c_s drops from N = 100 to 400 and then stops moving within errors for η ≥ 0.70; at η = 0.67 it keeps drifting toward the fluid EOS.
10. **Sept 6–7 — pressure campaign.** New runner and strict contract (any health event anywhere → discard). η ≤ 0.5 done: 96/96 accepted. Dense grid stalled because `edmd_init_random_gas()` is unbounded rejection sampling and hangs above η ≈ 0.55. New `edmd_init_lattice_gas()` in the core; four defects found in it via the validity checks (empty event calendar → inert run reported as valid with Z = 1; row-major corner packing; disks placed exactly touching the wall; lattice degeneracy needing jitter). Random vs lattice seeding agree on Z to ±0.25% at η = 0.2–0.5. Seeder regression tests written (`test_lattice_seeder.c`, `lattice_smoke.c`, `test_runner_and_smoke.sh`); dense grid not relaunched as of the files I saw.

---

# Part C — the eight findings, each with evidence

## C1. The chat's geometry was off by two

**ELI5.** The chat said "100 coins in a nearly square room". Actually it is 50 coins in each of two rooms, and neither room is square: room A is a tall narrow corridor (5 coins wide, 10 tall), room B a low wide one (10 wide, 5 tall). Sound travels along the 5-wide direction in A and along the 10-wide direction in B.

**In detail.** See A5. The table in the chat (L₀/σ = 5.45, H/σ = 10) was already per compartment; the prose ("11.2 × 10 diameters") was not. Consequence: route A is short *along* the propagation direction, route B is narrow *across* it. Both are 50-disk boxes.

## C2. Route B is a confinement experiment, not a second validation

**ELI5.** If you measure the stiffness of a crowd in a corridor only five people wide, the walls do half the pushing. Route B's room gets narrower (in coin sizes) as you raise the density, and its answer drifts up exactly as the room narrows — even at densities where the coins are still a liquid. Route A stays on the textbook curve.

**In detail.** Chat values vs EOS-derived adiabatic c_s (Henderson; SHY within 1.5%):

| η | EOS | route A | A/EOS | route B | B/EOS | B compartment height H/σ |
|---|---|---|---|---|---|---|
| 0.026 | 1.491 | 1.504 | 1.009 | 1.504 | 1.009 | 27 |
| 0.157 | 1.985 | 2.001 | 1.008 | 1.986 | 1.000 | 11.2 |
| 0.393 | 3.752 | 3.775 | 1.006 | 3.849 | 1.026 | 7.1 |
| 0.550 | 6.655 | 6.782 | 1.019 | 7.075 | 1.063 | 6.0 |
| 0.630 | 9.658 | 9.581 | 0.992 | 11.657 | 1.207 | 5.6 |
| 0.700 | 14.4 | 19.30 | 1.34 | 19.86 | 1.38 | 5.3 |

Route B's excess over the EOS grows monotonically as H/σ shrinks, and at η = 0.63 its ψ₆ is 0.18 — a liquid. So it is slit confinement (wall layering makes the in-plane response stiffer), not crystallisation. The statement "both routes match Liu to ~2% below 0.55" holds for A only. Fix: plot A + strip ladder as the validation, B as a confinement series with H/σ on the axis; quote residuals vs a named EOS.

## C3. Route B at η ≥ 0.75 is a square lattice, seeded by the initializer

**ELI5.** If you tile a room 5 coins tall with coins, a honeycomb doesn't fit — 5 honeycomb rows only hold 48 coins, and you have 50. A square grid of 10 × 5 does fit, just barely. The starter code lays down the square grid, and at that density nothing can move, so the square grid is what you measure. Room A is the opposite: the square grid doesn't fit (5 × 9 = 45), a honeycomb of 11 rows fits exactly 50. So A starts as a honeycomb and stays one.

**In detail.** `routeB/eta_0p750/speed_of_sound_psi6.csv`: ψ₆ global 0.12, local 0.24, mean neighbours 3.98 (225 runs). Four neighbours and near-zero six-fold order = square packing. The initializer (per the 26_08_21 note: checked rectangular if feasible, else hex) — my arithmetic with usable width = box − σ and ε = 10⁻³: B at η = 0.72: rectangular 10 × 5 fits with spacings 1.049 σ and 1.055 σ; hex rows: 5 rows hold 10+9+10+9+10 = 48 < 50. A at η = 0.72: rectangular 5 × 9 = 45 < 50, so hex: 11 rows, 6·5 + 5·4 = 50 exactly. Verify against `initialize_simulation()` in `00ALLINONE.c` (this is not the new core seeder). Fix: compute ψ₄ for B ≥ 0.72 (expect ≈ 1); drop those points from any c_s(η) figure or label them "seeded square jam".

## C4. At η ≥ 0.70 the structure changes while the measurement runs

**ELI5.** You are measuring the note of a bell while the bell is still being cast. In several boxes the coins are rearranging from disordered to honeycomb (or the reverse) *during* the swing you are timing, and heavier dividers swing longer, so they see more of the rearrangement. The "speed of sound" you fit from that is a number for a moving target.

**In detail.** ψ₆ global at end of hold → end of run, mean over runs:

| case | release | end | note |
|---|---|---|---|
| A 0.63 | 0.33 | 0.33 | stationary |
| A 0.70 | 0.38 | 0.60 | orders during run |
| A 0.72 / 0.75 | 0.91 / 0.97 | 0.90 / 0.97 | crystal from the start |
| B 0.63 | 0.18 | 0.18 | stationary liquid |
| B 0.70 | 0.25 | 0.33 | slowly ordering |
| B 0.72 | 0.26 | 0.69 | crystallises during run; bimodal end (118/225 ≥ 0.8, 36 < 0.4) |
| B 0.75 | 0.12 | 0.12 | square jam |
| strip N=200, 0.72 | 0.53 | 0.23 | global order falls, local stays 0.78: domains |
| strip N=400, 0.70 / 0.72 | 0.63 / 0.72 | 0.64 / 0.73 | stationary |

In B at 0.72 the end-of-run ψ₆ tracks the divider mass (0.50 at M = 50 → 0.79 at M = 300 → 0.61 at M = 2000), and mass is the fit variable. So the structural state enters the slope. The hold (≈800 time units) is not enough at η ≥ 0.70 for 50 disks. Fix: hold → end audit of every campaign including famA/famC (same CSV columns); flag |Δψ₆| > 0.1 or bimodality; then hold-length test 2000 → 20000 steps at 0.70/0.72 for A and famA N = 400, recording ψ₆ during the hold.

## C5. The strip ladder is not a finite-size scaling

**ELI5.** Making the corridor longer while keeping it 10 coins tall is not "a bigger room", it is a longer corridor. The walls are just as close as before.

**In detail.** H = 10 σ fixed, L₀ ∝ N: the 1/√N extrapolation (c∞ = 10.44 at η = 0.70) mixes aspect-ratio change with size. CC dropped it correctly after the fixed-aspect campaign; the figure `finite_size_scaling_with_N1000.png` still looks final on disk. Fix: mark superseded (README or rename).

## C6. "Drops then plateaus" is generous

**ELI5.** With three sizes and error bars about as big as the differences, "it stopped changing" and "it is still drifting slowly" look the same. You fit a line and quote the intercept with an error bar instead of saying "plateau".

**In detail.** Family A: η = 0.67: 13.3 / 12.2 / 11.9 / 11.4 (N = 100/400/900/1600; EOS 11.3) — still falling toward the fluid value; η = 0.71: 15.2 / 12.4 / 13.4 / 13.3; η = 0.73: 18.5 / 14.2 / 16.4 / 14.9, errors ±0.5–1. Fix: fit c(N) = c∞ + a/√N over N = 400–1600 per η, report c∞ ± σ. The 1/√N ansatz is justified by the wall-layer (perimeter/area) mechanism; ψ₆ at 0.67 decays like N^−0.55 ≈ 1/L, consistent with it.

## C7. Above η ≈ 0.70 there is no benchmark on the plot

**ELI5.** The textbook curve is for a liquid. At those densities the coins are (or want to be) a crystal, and a crystal carries sound faster because it resists shear too. Sitting above the liquid curve is not a failure; it is a missing reference.

**In detail.** Bulk hard disks: liquid ≤ 0.700, liquid–hexatic coexistence 0.700–0.716, hexatic to ≈0.720, solid beyond (Bernard–Krauth 2011; Engel et al. 2013). For a solid the longitudinal speed is c_L² = (K + μ)/ρ (2D), K the bulk modulus and μ the shear modulus; the fluid formula has μ = 0, so "Liu solid, bulk mode" is a lower bound (the plot says so). Fix: either build c_L from the solid EOS plus a literature shear modulus for hard-disk crystals, or label the region "no benchmark". Also: in coexistence the divider's oscillation perturbs the phase balance, so c_s in a two-phase box is not a material constant, and ψ₆ = 0.4–0.65 at N = 900–1600 says the simulated systems are ordered partly because hard walls order a hard-disk fluid before the bulk does.

## C8. Provenance items

- `OVERNIGHT.log`: η = 0.700, M = 750, `health=1`, `invalid=0`. Under the strict rule adopted for the pressure campaign this trajectory is discarded. Identify the counter, exclude it, reconcile the two runners' validity rules.
- The overnight script crashed after the merge (`line 135: unexpected EOF`); the `an/` analysis and the FINAL PNGs came from a separate step with no `00_COMMAND.md`. Add one.
- `--edmd-acc=1`: failed validation (181/450 invalid, per the chat export). Excluded from everything scientific until it passes the same gates; a c_s comparison is moot until then.
- The analysis uses cot K = (M/2Nm)K with N per side and L_eff = L₀ − 2r; confirm the exact line in `analyze_speed_of_sound_by_eta.py` and show per-mass residuals for N ≥ 900.

---

# Part D — pressure / Z: what exists and what I could not check

**ELI5.** I read the recipe and the bug reports for the pressure measurement, not the measurement itself. The recipe is good. The low-density half is reported done and clean. The high-density half hasn't started because the starting-arrangement code was fixed only last night. I need the folder the results were written to.

**In detail.** `run_pressure_campaign2.sh` says the η ≤ 0.5 grid is complete (96/96 accepted, 0 discards) and "preserved separately"; the output folder is a command-line argument and nothing named `pressure*` exists at the top of HardDisks/ or hspist3/. With that OUTDIR I would check: Z_pair vs Kolafa–Rottner at each (η, N) with block errors; the three-estimator spread; the wall-term scaling (Z_wall − Z_pair) vs 1/√N and its extrapolation; block-length adequacy (autocorrelation); ψ₆ per block stationarity; the seed-to-seed scatter vs the block error. For the dense grid, the acceptance test is the same plus: lattice vs random seeding agreement where both work (already ±0.25% at 0.2–0.5), and ψ₆ decay during equilibration at η ≥ 0.60.

---

# Part E — Paper 1, section by section

**ELI5.** Paper 1 says: "this program, with this self-checking, reproduces what everyone already knows about coins on a table, and here is exactly where small boxes stop telling you about big ones." Boring on purpose. It is the licence for Paper 2.

**In detail — sections, what goes in, what "done" means.**

1. **Model and validity contract.** EDMD, geometry, units, the ledger, the gates (initialization, per step, completion, batch, external verifier). Done when the list of gates in the paper matches `experiment_validation.c` and `pressure_validation.c` line by line.
2. **What the contract caught** (the software result): the seed-2381038820 missed event; the harmonic-wall coarse scan; the zero-gap recollision; the forced advance that gave Z = 4.1476 vs clean 4.149–4.150 at η = 0.5, N = 1600 (plausible number, thousands of repairs); the inert lattice seed reporting Z = 1 as valid; corner packing; wall-contact seeding; the random-seeder hang. Done when each has a one-line reproduction and the test that now guards it.
3. **Z(η).** Two independent routes (pair virial; wall momentum flux, with x-vs-y as the wall route's isotropy check — GPT's wording, correct), N = 400/900/1600, wall term extrapolated, vs Kolafa–Rottner for η ≤ 0.69; 0.702–0.720 shown as exploratory. Done when: dense grid run with the tested seeder; per (η, N) block errors; wall-term fit; lattice-vs-random check at ≥ 2 dense η.
4. **c_s(η), fluid branch.** Route A + strip ladder N = 100–400 + N = 1000 at 0.63, vs EOS-derived c_s; Román reproduction; route B as confinement. Done when: the 8 fixes in Part C are in; residual table vs a named EOS; peak-on-boundary points removed; health = 1 trajectory excluded.
5. **Finite size and ordering, η = 0.67–0.73.** Fixed-aspect family, c∞ ± σ fits, ψ₆(N), hold → end stationarity audit, the strip lesson, seeding/commensurability. Done when the stationarity audit is clean for every point that is plotted, and the seeding-swap test has been run (it decides how the section is worded).
6. **Equilibrium fluctuations of the divider** (new, cheap, the bridge to Paper 2). From the hold and ring-down phases already on disk: ½M⟨v_w²⟩ = ½k_BT and a Maxwellian v_w histogram (exact, ensemble-independent); the spring divider's position variance vs k_BT/k_eff with k_eff = k + 2 N m c_s²/L² (the gas stiffness; note the isothermal/adiabatic distinction and that the box is microcanonical); and a friction coefficient — from the damping of the divider's oscillation (ring-down decay rate γ/2M, or the linewidth of the FFT peak), cross-checked with the force autocorrelation on the *held* divider, γ = β∫₀^∞⟨δF(0)δF(t)⟩dt (the Sivak–Crooks friction; mind the plateau-value caveat for a finite box). Not the velocity autocorrelation: for a confined divider its time integral is zero. Done when the three numbers are tabulated at 3–4 η with errors.

Not in Paper 1: thermostats, driven protocols, TIME mode, Szilard.

---

# Part F — Paper 2: the validation ladder for the driven piston

**ELI5 of the whole ladder.** Nobody has a textbook curve for "push a piston into coins at Mach 2 and see how much ends up in a spring". So you cannot validate that directly. What you *can* do is climb a ladder where each rung is something exact: (0) energy is never created or lost; (1) a free divider in a calm gas jiggles exactly as much as temperature says; (2) if the coins don't collide with each other (ideal gas) the work done by a piston has a known probability distribution, and a famous identity (Jarzynski) holds for any push speed; (3) walls that are supposed to be radiators must actually be radiators; (4) for gentle slow pushes the whole piston–gas–divider–spring chain is a linear system whose answer you can write down with numbers from Paper 1; (5) only then the real experiment, where deviations from (4) are the new physics, and the information-theory relation (nostalgia ↔ dissipation) can be tested.

## Level 0 — identities that must hold to roundoff

**ELI5.** Count the money before and after. If the books don't balance to the cent, nothing else matters.

**In detail.**

1. Write the single equation W_piston = ΔKE_gas + ΔKE_walls + ΔE_spring + Q_bath, with the sign convention of *every* accumulator in `00ALLINONE.c` (`W_pistonL`, `W_pistonR`, spring energy, wall KE, thermal-wall heat). The audit says these were never put in one equation; do that first.
2. Add the residual per trajectory to trace and summary. In EDMD with elastic walls every term is exact between events; the relative residual should be ~10⁻¹⁰, and its distribution over `validated_reference_20260820_all_gates_r1` should be reported before it becomes a gate. If it is 10⁻³ somewhere, that is a bug in an accumulator, not "numerical error".
3. Momentum: piston impulse = Σ Δp of the gas.
4. Time reversal: at time τ flip every velocity (particles and walls), integrate τ back, compare to the initial state. Exact dynamics is reversible; chaos amplifies roundoff like e^{λτ}, so this works for short τ only, but it directly tests event ordering, which is the class of bug behind seed 2381038820.

## Level 1 — the divider in a calm gas (Paper 1, section 6)

**ELI5.** Before you push anything, check that the divider jiggles the way a warm object should.

**In detail.** Equipartition ½M⟨v_w²⟩ = ½kT; Maxwellian v_w; position variance kT/k_eff; friction from the ring-down / linewidth and from the held-divider force autocorrelation γ(λ) = β∫₀^∞⟨δF(0)δF(t)⟩_λ dt. That last quantity is exactly the "thermodynamic metric" of Sivak & Crooks (PRL 2012), computed at fixed control parameter λ (piston or divider position); it is what Level 4 needs. Do it for the piston too, held at several positions.

## Level 2 — exact results in the ideal-gas limit

**ELI5.** If the coins ignore each other, each coin just bounces between the walls and the piston on its own, and you can work out by hand what a moving piston does to it. Add up 5 or 10 coins and you have an exact prediction for the whole box. And there is a beautiful identity (Jarzynski): however fast and wasteful you push, if you average e^{−W/kT} over many pushes you get exactly the answer for an infinitely slow push. Fast pushes are allowed to be wasteful on average, but not in that exponential average. It is a razor: any error in the piston–coin rule, the work bookkeeping or the sampling breaks it.

**In detail.**

1. Single collision with a piston moving at speed u into the gas: a particle with x-velocity v (toward the piston) leaves with −v − 2u; energy gain 2mu(v + u). Work = Σ over collisions. With particle–particle collisions off (`simple_prediction_box` already does this; or `--no-pp-collisions`) the particles are independent, so P(W) for N particles is the N-fold convolution of the single-particle distribution, which is computable (Lua & Grosberg 2005 do the 1D/3D cases; the adiabatic-compression work-distribution paper in `ZZZ_PAPER/SIMPLE_GAS_BOX` gives the dilute-gas result). Test: histogram of W over ≥ 2000 seeds vs the analytic P(W) at three piston speeds, N = 5, 10, 20.
2. Jarzynski: prepare a canonical ensemble (validated thermal wall, then switch walls to elastic), drive with the piston, compute ⟨e^{−βW}⟩ and compare to e^{−βΔF} with ΔF = −N k T ln(A_f/A_i) for the 2D ideal gas. Report bootstrap errors and the dependence on sample size (the exponential average is dominated by rare low-W trajectories; that is why small N first).
3. Slow limit: ⟨W⟩ → N k T_i (A_i/A_f − 1) (2D adiabatic ideal gas, γ = 2 ⇒ T·A = const); the excess ⟨W⟩ − W_qs should vanish linearly in piston speed.
4. Same with collisions on at η = 0.2 and 0.5, N = 50: ΔF from the EOS, β ΔF/N = ln(η_f/η_i) + ∫_{η_i}^{η_f} (Z − 1)/η dη, with Kolafa–Rottner and with the campaign's own Z. This is the direct use of Paper 1 inside Paper 2.

## Level 3 — the heat bath has to be a heat bath

**ELI5.** A radiator wall must give back, on average, exactly the energy that hits it when the room is at the right temperature — not a bit less, and it must not switch itself off when the room's average temperature looks fine. As documented, the current wall gives back too little (it under-emits by a third), only turns on when the whole room's temperature drifts, and in its default mode deliberately overshoots. That is a room thermostat, not a radiator.

**In detail.** As documented in `THERMAL_WALL_IMPLEMENTATION.md`: speed drawn from the 2D Maxwell *speed* distribution p(v) ∝ v e^{−mv²/2kT}, direction uniform in the half-plane; activated only if |T_gas − T_bath| > 10⁻⁴ T_bath using the global temperature; default mode overshoots toward the target. The diffuse (Maxwell) wall that satisfies detailed balance re-emits the normal component with p(v_n) ∝ v_n e^{−mv_n²/2kT} and the tangential component Gaussian; the emitted mean energy is 3kT/2 in 2D, equal to the flux-weighted mean energy of the particles arriving, whereas the documented kernel emits kT. So, as documented, it cools the gas, violates detailed balance, and gives a wrong angular distribution — and the threshold/adaptive logic masks that by feedback. Verify in the code; if it is as documented, replace with the diffuse Maxwell wall (Tehver et al. 1998 is the standard reference) or a per-collision, unconditional Andersen kernel; log heat per wall collision. Validation: (a) T_gas → T_bath with no feedback, (b) Maxwellian bulk, (c) zero mean heat flux at equilibrium over long runs, (d) two walls at T₁ ≠ T₂ with PP collisions off: the free-molecular heat flux is exactly computable (each particle shuttles between the walls, carrying the difference of the two flux-weighted mean energies per crossing), (e) Level 2's Jarzynski with the bath on during the drive, then Crooks: P_F(W)/P_R(−W) = e^{β(W − ΔF)} for compression vs expansion.

## Level 4 — linear response: where "energy transfer" gets a prediction

**ELI5.** Push the piston gently and slowly and the whole chain — piston, gas, divider, spring, gas — behaves like connected springs and masses. You know every spring (the gas stiffness is c_s from Paper 1), every mass, and the friction (Level 1). So you can write down how much energy the spring should get for a given push, with no fitted numbers. Then push harder and faster and report *how far* you are from that answer; that is the new physics.

**In detail.** For small amplitude and slow driving the gas column is an acoustic transmission line with impedance ρ c_s and length L; the divider is a mass M with spring k and friction γ; the Román transcendental relation is the one-column special case. Derive the transfer function from piston motion to spring energy, predict SpringE(t) for the actual protocols, and compare with small-amplitude runs (Mach ≤ 0.1, injected work ≤ 0.1 kT per particle). Fluctuation–dissipation (Sivak & Crooks 2012): for a slow protocol λ(t) the mean excess work is ⟨W_diss⟩ ≈ ∫ γ(λ) λ̇² dt with γ(λ) the force-autocorrelation friction from Level 1 measured at fixed λ — an exact linear-response statement with no fitted parameter, and the natural first result of the non-equilibrium paper because it ties the driven measurement to an equilibrium one. Current runs: `--piston-speed 100` capped at 10 σ in 2 time units (≈5 thermal speeds, Mach ≈ 2 at η = 0.2), `--piston-work-target 20000` at N = 600 (≈33 kT per particle, temperature up ~30×). Keep them as the far end of a speed ladder Mach 0.05 → 5 at fixed injected work. And replace "SpringE_max / W_in" (the maximum of an oscillating quantity inside a 20-unit window) by a defined delivered work or by the time-resolved SpringE(t) compared with the prediction.

## Level 5 — Thermodynamics of Prediction (the simple prediction box)

**ELI5.** The piston moves to a new random spot, the gas scrambles to adjust, heat leaks to the walls, repeat. The gas "remembers" where the piston was. If the piston's next move is predictable, some of that memory is useful; the rest is nostalgia — remembering the past for no benefit — and Still et al. showed nostalgia costs energy. To test it you need: correct heat accounting (Level 3), a protocol with tunable predictability, and honest information estimates.

**In detail.** x_t = piston target at step t; s_t = coarse gas features after relaxation; I_mem = I(s_t; x_t), I_pred = I(s_t; x_{t+1}), nostalgia = I_mem − I_pred. Still et al. relate the work dissipated in the work step to kT × nostalgia (with their definition of dissipated work through the nonequilibrium free energy and an instantaneous work step; work out from the paper how it maps onto a finite-time drive and onto the mechanically measured W − ΔF_eq before designing the test). Practical requirements: the drive is finite-time (no teleporting — settled), so add a speed ladder and extrapolate; mutual information from binned features is biased upward at small samples and is a lower bound on the microstate value (data processing), so the testable statement is an inequality with bias-corrected estimators and shuffle nulls over thousands of steps; iid vs Markov protocols (`--simple-box-target-protocol=iid|correlated`) give the two limiting cases. The ideal-gas version (PP off) is where the numbers can be cross-checked semi-analytically because the single-particle dynamics between moves is exactly solvable.

## Level 6 — Szilard

Separate; same contract; needs protocol-aware crossing rules (the audit says so). Not in Paper 2.

---

# Part G — the next two weeks as a checklist

Each item names its acceptance criterion. Nothing new is launched until the item above it passes.

**Week 1 — close Paper 1 data**

1. Send me (or open for CC) the pressure OUTDIR with the 96/96 dataset. Criterion: Z_pair vs Kolafa–Rottner residuals tabulated with block errors; wall-term vs 1/√N fitted.
2. Run `test_runner_and_smoke.sh` and the seeder tests; relaunch the dense Z grid. Criterion: all tests pass; 0 discards or every discard explained.
3. Speed-of-sound stationarity audit (ψ₆ hold → end, ψ₄, neighbour count) over every campaign incl. famA/famC. Criterion: a table; flagged cells removed from c_s(η) plots.
4. Refits: famA c∞ ± σ; per-mass residuals for N ≥ 900; exclude the health = 1 trajectory; mark the strip extrapolation superseded; add the missing `00_COMMAND.md`. Criterion: new FINAL figure with a named EOS, residual table, provenance box (N, r, L₀, H, seeds).
5. Seeding-swap and hold-length tests at η = 0.72 (N = 100, minutes). Criterion: ψ₆ at release and at end, ψ₄, c_s ± σ for each of the four cases; this decides the wording of section 5.
6. Divider fluctuation section from existing holds/ring-downs. Criterion: ⟨v_w²⟩/(kT/M), position variance ratio, γ at 3–4 η.

**Week 2 — Paper 2 groundwork (no physics runs before 7 and 8 pass)**

7. Energy accounting audit and residual (Level 0). Criterion: residual distribution over the 20-run reference at roundoff level; then it becomes a gate.
8. Thermal wall audit (Level 3, items a–d). Criterion: T_gas = T_bath without feedback; Maxwellian; zero mean heat flux; free-molecular flux matches.
9. Ideal-gas exact tests (Level 2, items 1–3). Criterion: P(W) histogram vs analytic; Jarzynski within bootstrap error at three speeds; slow limit.
10. Finite-η Jarzynski vs EOS (Level 2, item 4). Criterion: ΔF agreement within error at Mach ≤ 0.2.
11. Linear-response baseline (Level 4). Criterion: SpringE(t) from small-amplitude runs on top of the prediction; friction from FDT equal to Level 1's γ.
12. Only then: the speed/amplitude ladder including the existing Mach-2 grid, and the prediction box.

---

# Part H — glossary and formulas

- **η** packing fraction = disk area / box area, per compartment.
- **σ** disk diameter = 2r; **Z** = P A/(N k T), the compressibility factor.
- **EOS** equation of state, Z(η). Henderson: (1 + η²/8)/(1 − η)². Kolafa–Rottner 2006 for precision. Liu 2021 for the global (fluid/hexatic/solid) curve.
- **c_s** adiabatic sound speed: c_s² = (kT/m)[Z + ηZ′ + Z²]; ideal-gas limit √(2kT/m).
- **Román relation**: ν = c_s K/(2π L_eff), cot K = (M/2Nm) K, N per side, L_eff = L₀ − 2r.
- **ψ₆**: local ⟨|ψ₆,i|⟩, global |⟨ψ₆,i⟩|, neighbours within 1.4 σ; square lattice → 0, hexagonal → 1.
- **Phases (bulk)**: liquid ≤ 0.700; coexistence 0.700–0.716; hexatic to ≈0.720; solid beyond.
- **c_L** in a solid: c_L² = (K + μ)/ρ; fluid has μ = 0.
- **Finite-size ansatz**: c(N) = c∞ + a/√N at fixed aspect ratio (perimeter/area).
- **Virial (2D hard disks)**: Z_pair = 1 + m σ Σ|Δv_n| / (2 KE t).
- **Moving-piston collision**: v → −v − 2u; ΔE = 2mu(v + u).
- **Jarzynski**: ⟨e^{−βW}⟩ = e^{−βΔF}; **Crooks**: P_F(W)/P_R(−W) = e^{β(W − ΔF)}; ideal 2D gas ΔF = −NkT ln(A_f/A_i); hard disks β ΔF/N = ln(η_f/η_i) + ∫(Z − 1)/η dη.
- **Adiabatic 2D ideal gas**: γ = 2, T·A = const, W_qs = N k T_i (A_i/A_f − 1).
- **Diffuse wall**: p(v_n) ∝ v_n e^{−mv_n²/2kT}, v_t Gaussian; emitted mean energy 3kT/2 in 2D.
- **Friction / thermodynamic metric (Sivak–Crooks 2012)**: γ(λ) = β∫₀^∞⟨δF(0)δF(t)⟩_λ dt at fixed control parameter; slow-protocol excess work ⟨W_diss⟩ ≈ ∫γ(λ)λ̇² dt. (For a confined divider the velocity-autocorrelation integral is zero, so it cannot be used for this.)
- **Nostalgia**: I(s_t; x_t) − I(s_t; x_{t+1}).
