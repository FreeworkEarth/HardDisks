# Speed of sound campaign — recap of what was done, and a critical review

Date: 2026-09-07 (Cowork / Claude, "max" pass)
Scope: the route A / route B chat Chris pasted (Aug 25), plus what is actually on disk after it (ladder, overnight N=1000, accelerated-core check, fixed-aspect-ratio campaign, all Aug 25–26), plus the pressure-validation thread from last night (260908_speedsound_validation_low_high_density.md).

Everything below that carries a number was read from the files on disk, not from the chat. Paths are listed at the end. What I could NOT read: anything deeper than 7 folders below the connected HardDisks folder (all `analysis*/final_plots/*.csv`, the route A/B overlay PNG, the fixed-aspect famA/famC per-η results). If you add `.../00_eta_sweep_ROMAN` as a connected folder in the desktop app I can read those too.

---

## 1. What happened, in order (the "last 10 messages", plus the two days after)

**Setup that all of this uses.** `00ALLINONE --mode=edmd --experiment=speed_of_sound`: two compartments separated by a movable divider ("wall"), each compartment L0 long and H high, `--particles-boxes=50,50` so N=100 means **50 disks per compartment**. The divider is held for `--wall-hold-steps=2000` (with `--fixed-dt=0.4`, so presumably ~800 time units of equilibration with the wall fixed — worth confirming what one hold step is), then released; its x(t) is recorded for ~25 predicted oscillations, FFT'd, and the peak frequency is fitted against a mass term for 9 divider masses (M = 50…2000 m). The slope is c_s. This is a Rüchardt-type measurement: the gas column is the spring, the divider is the mass. Validity ledger per trajectory (overlap / clamp / forced-advance counters).

**Step 1 — Route A (campaign_r25_psi6_20260823).** r = 0.5, H = 10, L0 chosen per η (L0 = 200 at η=0.02 down to 5.17 at η=0.76), 32 η values, 9 masses × 25 repeats = 225 trajectories per η, 7200 total, 0 invalid. ψ₆ (global and local, at end of hold and at end of run) was added to the output — that's the "psi6" in the folder name.

**Step 2 — Route B (routeB_radius_N100_L0_20_20260825).** Same code, same seeds, but L0 = 20 and H = 10 fixed and the **radius** is swept (r = 0.158 → 0.984) to set η. Same 32 η, 7200 trajectories, 0 invalid. Purpose: a second way to reach the same η with a different box, to test whether c_s is intensive.

**Step 3 — the broken route-B analysis.** The grouped analysis groups by L0; in route B L0 is 20 for every η, so all 32 η collapsed into one group and the summary had 1 row. It also used the default r = 0.5 for L_eff = L0 − 2r instead of the real r. That analysis (`routeB.../analysis/`) was deleted and redone per η with the correct radius (`analysis_pereta/`, merged into `analysis_combined/`). The file you had opened ("FINAL speed_of_sound_on_packing_fracture.pdf" under `analysis/`) was the broken one.

**Step 4 — A vs B comparison.** Below η ≈ 0.55 the two routes agree to a few percent; above, they diverge (B/A = 1.22 at 0.63, 1.03 at 0.70, 1.38 at 0.73) and ψ₆ differs (A crystallises, B doesn't). CC's reading: at high density you measure the structure the box imposes, not a bulk property; the working split "ready below 0.55 / not ready above" comes from this.

**Step 5 — timing N=1000.** One trajectory at N=1000, η=0.70 takes 391 s with the default O(N²) core (2.2 s at N=100; scaling ~N^2.2). A full curve at N=1000 would be ~780 CPU-hours, so instead:

**Step 6 — the "strip ladder" (ladder_N100/N200/N400_20260825).** N = 100, 200, 400 at η = 0.40, 0.55, 0.63, 0.70, 0.72, 9 masses × 10 repeats. r = 0.5 and **H = 10 fixed**, L0 grown ∝ N (5.45 → 10.9 → 21.8 at η=0.72). Result: flat in N at η ≤ 0.63 (matches the EOS); at 0.70 c_s falls 19.2 → 17.0 → 14.8; at 0.72 it rises 16.7 → 18.9 → 19.8.

**Step 7 — accelerated core check (validate_acc_N100_20260826).** Same as route A at 5 η but with `--edmd-acc=1` (cell lists). Per the chat export it FAILED: 269 valid / 181 invalid of 450 (90 boundary escapes, 91 overlaps) where the default core gives 450/450; excluded. (Earlier text here said "0 invalid" because I had looked at one η's status file only.)

**Step 8 — overnight N=1000 (overnight_N1000_20260826).** N=1000 (500 per side), r=0.5, H=10, L0 = 62.3/58.6/56.1/54.5 at η = 0.63/0.67/0.70/0.72, 9 masses × 18 repeats = 648 trajectories, ~6 h on 12 jobs. Log says invalid=0 everywhere, but one cell reports **health=1** (η=0.700, M=750). The runner script died at the end with a shell syntax error (line 135, unmatched quote) after merging; the analysis under `an/` was produced by a separate step. Result: at 0.70 c_s keeps falling (13.3 at N=1000), at 0.72 it goes back down to 17.3. A 1/√N extrapolation was drawn (c∞ = 10.44 at 0.70, R²=0.99; the other η's have R² ≤ 0.5).

**Step 9 — fixed-aspect-ratio campaign (finitesize_aspect_20260826).** This is the one that changes the picture. Family A: r=0.5, L0 and H both ∝ √N (H = 20/30/40 for N = 400/900/1600, whole system ≈1.1:1), η = 0.67…0.73 (12 values), 5 masses at N = 900/1600 (9 at N=400), 2540 trajectories, 0 invalid. Family C: L0, H fixed, r ∝ N^(-1/2) at η=0.710 only (a route-B-style control with proper aspect). Results (from the PNGs): c_s drops from N=100 to N=400 and then stops moving within error bars for η ≥ 0.70 (12.5–15); at η=0.67 it keeps drifting down (13.3 → 12.2 → 11.9 → 11.4, EOS ≈ 11.3). Families A and C agree at η=0.710 (13.3 vs 13.7 at N=1600). The strip ladder does not agree with either. ψ₆ vs N: η=0.67 melts to a liquid as N grows (0.55 → 0.12); η ≥ 0.70 stays at 0.4–0.65 even at N=1600. CC's own title on that figure: "the strip result does not survive".

**Step 10 — last night's pressure-validation thread (the 260908 markdown).** Different runner (`validation/pressure_validation.c` on `edmd_core` directly, not 00ALLINONE), measuring Z(η) with a strict contract: any forced advance / overlap repair / clamp repair / wall overdue anywhere in a trajectory discards it. 96/96 accepted at η ≤ 0.5. The dense grid stalled: root cause was `edmd_init_random_gas()` — unbounded rejection sampling that hangs above η ≈ 0.55 (the "40-minute η=0.60 runs" never left seeding). Fix: a new `edmd_init_lattice_gas()` in the core (checked rectangular, hex fallback, jitter, wall inset), following the 26_08_21 hex-fallback note; `edmd_init_random_gas()` left untouched. Four bugs found on the way, all in the new seeder or the validity rule (empty event calendar → inert run reported valid with Z=1; corner packing; exact wall contact; lattice degeneracy). Random vs lattice seeding agree on Z to ±0.25% at η = 0.2–0.5. GPT's list of seeder regression tests is queued; dense campaign not relaunched yet. Two things from that thread matter for the speed-of-sound work: the same "plausible number, silent failure" pattern, and GPT's rule that **stationary T and P do not prove equilibration — look at ψ₆ over time**. That rule has not been applied to the speed-of-sound runs, and it bites (section 2d).

---

## 2. Where the chat narrative is wrong or too generous (with evidence)

### 2a. The geometry described in the chat is off by a factor 2

The commands say `--particles=100 --particles-boxes=50,50`. So every "N=100" box is two compartments of 50. Per compartment:

| | route A, η=0.72 | route B, η=0.72 |
|---|---|---|
| command | `--lengths=5.4542 --height=10 --particle-radius=0.5` | `--lengths=20 --height=10 --particle-radius=0.957461` |
| compartment in diameters | **5.45 × 10.0** (≈5 columns × 10 rows) | **10.44 × 5.22** (≈10 columns × 5 rows) |
| disks | 50 | 50 |

Check: 50·π·0.25/(5.4542·10) = 0.720. The chat's "route A near-square box 11.2 × 10 diameters" is not right; A is a 5-diameter-long column along the direction the sound travels, B is a 5-row-high slit across it. Both are 50-particle boxes; neither is close to square. The table in the chat (L0/σ = 5.45, H/σ = 10) was already per-compartment, the prose wasn't.

### 2b. "Two routes agree to 1–3% below 0.55, both match Liu to ~2%" — only route A does

I recomputed the adiabatic sound speed from two fluid EOS (c_s² /(kT/m) = Z + ηZ′ + Z², Henderson and Santos–Haro–Yuste) and compared with the chat's numbers:

| η | EOS (Henderson) | route A | A/EOS | route B | B/EOS | B compartment height H/σ |
|---|---|---|---|---|---|---|
| 0.026 | 1.491 | 1.504 | 1.009 | 1.504 | 1.009 | 27 |
| 0.157 | 1.985 | 2.001 | 1.008 | 1.986 | 1.000 | 11.2 |
| 0.393 | 3.752 | 3.775 | 1.006 | 3.849 | 1.026 | 7.1 |
| 0.550 | 6.655 | 6.782 | 1.019 | 7.075 | 1.063 | 6.0 |
| 0.630 | 9.658 | 9.581 | 0.992 | 11.657 | 1.207 | 5.6 |
| 0.700 | 14.4 | 19.30 | 1.34 | 19.86 | 1.38 | 5.3 |

Route A is within 2% of the EOS all the way to η = 0.63. Route B is 3% high at 0.39, 6% at 0.55, 21% at 0.63 — and the excess grows monotonically as the slit gets narrower in particle diameters. At 0.63 route B's ψ₆ is 0.18 (liquid), so this is **not** crystallisation; it is a liquid confined to a 5–6 row slit, whose response to the divider is stiffer than bulk (wall layering). So route B is not an independent confirmation of the fluid branch; it is a confinement series. That's still useful, but it should be plotted and described as that. The "agreement to 1–3%" holds only where H/σ ≳ 7 (η ≤ 0.4).

Also: the B/A ratio is non-monotonic (1.22 at 0.63, 1.03 at 0.70, 1.38 at 0.73). Without per-point error bars nobody can say which of those are significant. The per-η summaries carry uncertainties; they belong in the table.

### 2c. Route B at η ≥ 0.75 is not "jammed disordered" — it is a square lattice

From `routeB.../eta_0p750/speed_of_sound_psi6.csv`, averaged over 225 runs: ψ₆_global = 0.12, **ψ₆_local = 0.24, mean neighbour count = 3.98**. Four neighbours and near-zero local sixfold order is a square (4-fold) packing, not a glass. It is almost certainly the seeded layout: the initializer (per the 26_08_21 note) places a checked rectangular lattice when one fits and only falls back to hex when it doesn't. In B's compartment at η=0.72–0.76 a 10 × 5 rectangular grid fits (spacing 1.02–1.05 σ), while hexagonal rows cannot: 5 hex rows in that height hold 48 disks, not 50 (my arithmetic from the spacing rules in the note — verify against `initialize_simulation()` in 00ALLINONE.c, which is not the new core seeder). So the box is commensurate with square order and incommensurate with hex, and the seed already is square. Route A's compartment at the same η is the opposite: rectangular doesn't fit (5 × 9 = 45 < 50), hex fits exactly (11 rows, 6·5 + 5·4 = 50), and ψ₆ is 0.91–0.97 at release. The A/B difference at high η is therefore mostly **seeding + commensurability**, not something the fluid chose. Those points (B at ≥ 0.75, plus the peak-on-boundary ones) should not be on any c_s(η) plot except as a labelled artifact. To make this airtight, compute ψ₄ from the stored configurations, it should be ≈ 1 for B at 0.75.

### 2d. At η ≥ 0.70 the structure changes *during* the measurement window

This is the most important finding and it is in files that already exist. ψ₆_global at end of hold (= at release) vs at end of run, mean over all runs:

| case (50 disks/compartment unless noted) | ψ₆ at release | ψ₆ at end | comment |
|---|---|---|---|
| A η=0.63 | 0.33 | 0.33 | stationary |
| A η=0.70 | 0.38 | 0.60 | orders during run |
| A η=0.72 | 0.91 | 0.90 | crystal from the start (hex seed) |
| A η=0.75 | 0.97 | 0.97 | crystal |
| B η=0.63 | 0.18 | 0.18 | stationary liquid |
| B η=0.70 | 0.25 | 0.33 | slowly ordering |
| B η=0.72 | 0.26 | 0.69 | **crystallises during run**; end values bimodal: 118/225 runs end above 0.8, 36 end below 0.4 |
| B η=0.75 | 0.12 | 0.12 | square jam (2c) |
| ladder N=200 η=0.72 | 0.53 | 0.23 | global order *falls*, local stays 0.78 → breaks into domains |
| ladder N=400 η=0.70 / 0.72 | 0.63 / 0.72 | 0.64 / 0.73 | stationary |

So for A at 0.70, B at 0.70/0.72 and the N=200 strip, the divider is oscillating inside a system that is still reorganising. The "c_s" from those cells is a number for a transient state that depends on how long the run was — and run length depends on divider mass (25 predicted cycles: heavier divider → longer run). In B at 0.72 the end-of-run ψ₆ indeed tracks M: 0.50 at M=50, 0.66–0.69 at 100–200, 0.71–0.79 at 300–1500, 0.61 at 2000. A mass-dependent structural state feeds straight into the frequency-vs-mass fit that gives c_s. This is exactly the situation GPT warned about for the pressure campaign, and 800 time units of hold is evidently not enough at η ≥ 0.70 for 50 disks. I could not check the fixed-aspect famA/famC CSVs (depth limit), but they have the same columns; the same hold → end comparison needs to be done there before anything from N = 400–1600 is called converged.

### 2e. The strip ladder is not a finite-size scaling

H = 10σ fixed, L0 ∝ N: the compartment goes from 5 × 10 to 54 × 10 diameters. That changes aspect ratio and confinement, not just size, and the 1/√N extrapolation (c∞ = 10.44 at 0.70) rests on it. CC already dropped it ("the strip result does not survive"), correctly. But the number and the "excess is finite-size, vanishes as the box grows" plot are still on disk in `overnight_N1000_20260826/` and look final. They should be marked superseded so they don't end up in a paper.

### 2f. "Drops then plateaus" is generous

Fixed-aspect family A: at η=0.67 the values are 13.3 / 12.2 / 11.9 / 11.4 for N = 100/400/900/1600 — still falling, towards the fluid EOS (11.3). At 0.71: 15.2 / 12.4 / 13.4 / 13.3. At 0.73: 18.5 / 14.2 / 16.4 / 14.9, error bars ±0.5–1. That is "no longer resolvable beyond the errors with three sizes", not a plateau. The honest statement comes from fitting c(N) = c∞ + a/√N over N = 400–1600 per η and quoting c∞ ± σ. At 0.67 it will land on the EOS; at ≥ 0.70 the uncertainty will be large, which is the true state of knowledge.

### 2g. Above η ≈ 0.70 the reference line is not the right one

Liu's fluid branch (and its homogeneous continuation through coexistence) is what the plots compare to. But for η ≥ 0.70 the equilibrium hard-disk system is in coexistence / hexatic / solid, and ψ₆ = 0.4–0.65 at N = 900–1600 says the simulated systems are ordered. A solid's longitudinal speed is c_L² = (K + μ)/ρ, which is above the fluid value by the shear-modulus term; the "Liu solid (bulk mode)" line on the FINAL plot is explicitly a lower bound. So the N ≥ 900 values sitting above every drawn line is not a discrepancy — it is the absence of a benchmark. Either compute c_L for the ordered phase (solid EOS + a literature shear modulus for the hard-disk crystal) or say plainly that there is nothing to validate against there.

Two more things about the ordered side: (i) hard flat walls order a hard-disk fluid before the bulk does, and ψ₆ at η=0.67 decays like ψ₆ ∝ N^(-0.55) ≈ 1/L — a wall-layer (surface/area) effect. That is the natural reason a 1/√N fit is the right ansatz for the fixed-aspect family, and also the reason the ordered state at 0.70–0.72 for N ≤ 1600 may still be partly wall-induced. (ii) In coexistence the divider's oscillation is itself a perturbation that can shift the phase balance; c_s in a two-phase strip is not a material constant.

### 2h. Small provenance items

- `OVERNIGHT.log`: η=0.700, M=750 has health=1 with invalid=0. Under the rule adopted last night (any health event → discard) that trajectory must be excluded from the N=1000 η=0.70 point, and the speed-of-sound runner's validity rule should be reconciled with the pressure runner's. Check which counter fired.
- The overnight script crashed after the merge; the `an/` analysis and the two FINAL PNGs came from a separate step with no `00_COMMAND.md` next to them. Every other campaign has one. Add it.
- Divider masses 50–2000 m are fixed while N grows to 1000 and 1600 (500 and 800 disks per side), so M/(N_side m) goes down to 0.06: the divider is lighter than the gas it pushes. The analysis uses Román's relation ν = (c_s/2πL₀)·K with cot K = (M/2Nm)·K (per README_ENTROPY_IMPLEMENTATION.md), which is exact for a 1D linear acoustic column at any mass ratio, so a light divider is not a problem *for the formula*. What is untested is whether a 1D column model still describes a compartment that is ordered/layered (η ≥ 0.70) or only 5 diameters long. Cheap check: freq-vs-kterm residuals for N ≥ 900 and per mass, and a refit with only M ≥ 2·N_side·m; if c_s moves, the model, not the statistics, is the problem.
- The accelerated core (`--edmd-acc=1`) failed validation (181/450 invalid) and is excluded; nothing at N ≥ 900 used it. (Corrected after reading the chat export; earlier text here was wrong read.

---

## 3. Where that leaves "paper-ready"

**Fluid branch, η ≤ 0.63 (probably 0.65): yes.** Route A, the strip ladder N = 100–400, and N = 1000 at 0.63 all sit within ~1–3% of EOS-based adiabatic c_s, 0 invalid trajectories, and (per CC, not re-checked here) the pipeline reproduces the Román benchmark. Present route A + ladder as the validation. Present route B as what it is, a slit-confinement series that converges to the bulk value for H/σ ≳ 7, which is a nice bonus figure, not a second validation. Quote the residual vs EOS as an error budget rather than "matches to ~2%".

**Transition region 0.67–0.73: not a c_s(η) number yet, but a defensible finite-size/ordering study once three things are done:** the hold → end ψ₆ stationarity check for every fixed-aspect cell (2d), a fitted c∞ ± σ instead of an eyeballed plateau (2f), and either a proper c_L benchmark or an explicit "no benchmark" (2g). The seeding/commensurability story (2c) is a real methods result and worth a paragraph, but as "the initial lattice and wall commensurability decide the structure in boxes this small", not as "the box imposes the structure" — the seed is doing half the work.

**Route B η ≥ 0.75 and all peak-on-search-boundary points: artifacts.** Drop from c_s plots or label.

On CC's "two papers" question: I'd call it one paper — validation of the method on the fluid branch, plus a finite-size section that is honest about the transition region — unless the seeding-swap test (below) turns the commensurability story into something sharp enough to stand alone.

---

## 4. What I would do next, cheapest first

None of this needs more N=1000 runs at fixed H.

1. **Analysis only, today.** For every campaign including famA/famC: table of ψ₆_hold → ψ₆_end per (η, N, M), flag cells where the change or the run-to-run spread is large. Add ψ₄ (or use the neighbour count already stored) to identify square states. Colour the c_s points by ψ₆_end on the final plots. This decides which existing numbers are usable without running anything.
2. **Seeding-swap test, N=100, minutes.** At η=0.72: route A started from a rectangular/jittered layout, route B started from a jittered hex layout (if the initializer can be forced; otherwise via a stored configuration). If ψ₆ and c_s follow the seed, "box-imposed" becomes "seed-imposed, box-locked". If they don't, the commensurability argument stands on its own.
3. **Hold-length test, cheap.** η = 0.70 and 0.72, route A and fixed-aspect N=400, `--wall-hold-steps` 2000 → 20000 (and record ψ₆ during the hold, not just at its end). If ψ₆ at release changes, all current ≥ 0.70 results are hold-limited.
4. **Refits, no new runs.** c∞ + a/√N over N = 400–1600 for each η in famA; mass-subset robustness (M ≥ 2 N_side m); exclude the health=1 trajectory; mark the strip-ladder extrapolation as superseded.
5. **Benchmark for the ordered side.** c_L from the solid branch (Liu solid EOS plus a literature shear modulus for hard-disk crystals) or an explicit statement that none is used.
6. **Only then**, if the drift at η=0.67–0.70 is still unresolved, one more fixed-aspect size (N=2500) at those two η only, because there the fluid EOS is a valid target and the extra size buys a real statement. With `--edmd-acc=1` once it is validated against the default core at N ≥ 400.
7. Not now: more strip-geometry runs, route-B beyond 0.70, N=1000 at fixed H.

---

## 5. Paste-able instructions for Claude Code

```
Before any new speed-of-sound runs:

1. Structural stationarity audit (analysis only).
   For every campaign under 00_eta_sweep_ROMAN (routeA campaign_r25, routeB, ladder_N*,
   overnight_N1000, finitesize_aspect famA and famC): from speed_of_sound_psi6.csv
   tabulate per (eta, N, M): mean psi6_global_hold, mean psi6_global_end, the
   run-to-run sd of psi6_global_end, mean neighbors_hold/neighbors_end.
   Flag any cell where |end - hold| > 0.1 or where end values are bimodal.
   Report the table. Do not change any analysis formula.

2. Identify square states. From stored configurations (or neighbor counts) compute
   psi4 for routeB eta >= 0.72 and for any flagged cell. Report psi4 next to psi6.

3. Provenance fixes.
   - overnight_N1000: identify the eta=0.700 M=750 trajectory with health=1 (which
     counter), exclude it from the N=1000 eta=0.70 point, re-run that analysis, and
     write the 00_COMMAND.md that produced the an/ analysis and the two FINAL PNGs.
   - Mark finite_size_scaling_with_N1000.* and the c_inf=10.44 extrapolation as
     SUPERSEDED (rename or add a README), per the fixed-aspect result.

4. Refits (no new runs).
   - famA: fit c(N) = c_inf + a/sqrt(N) over N=400,900,1600 per eta; report c_inf +- sigma.
   - Mass-subset robustness for N >= 900: refit using only M >= 2*(N/2)*m and report the
     shift in c_s. Show freq-vs-kterm residuals for N=900 and 1600, per mass.
   - Confirm in analyze_speed_of_sound_by_eta.py that the kterm is K/(2*pi*L_eff) with
     cot K = (M/(2*N_side*m))*K (Román), state which L_eff is used, and cite the line.

5. Two cheap N=100 tests at eta=0.72 (report psi6_hold, psi6_end, psi4, c_s +- sigma):
   a. seeding swap: route A geometry from a jittered rectangular start, route B
      geometry from a jittered hex start (or from a stored hex configuration).
   b. hold length: --wall-hold-steps 2000 vs 20000 for route A and famA N=400 at
      eta = 0.70 and 0.72, recording psi6 during the hold.

6. Reference for eta >= 0.70: either compute c_L = sqrt((K+mu)/rho) for the ordered
   phase from the solid EOS plus a literature shear modulus, with the source cited,
   or remove the fluid-branch reference lines from that region and label it
   "no benchmark".

Do not launch more N=1000 fixed-H runs, more route-B runs above eta=0.70, or any
new geometry until 1–5 are reported. Do not auto-commit.
```

---

## 6. Questions I'd put to GPT and CC

- Do you accept that route B above η ≈ 0.4 is a confinement measurement rather than a second validation? If not, what explains the monotonic B/EOS excess with shrinking H/σ at densities where B's ψ₆ is that of a liquid?
- Given ψ₆ moves from 0.26 to 0.69 during the run in B at 0.72, and depends on divider mass, what is the c_s fitted from those runs a measurement of?
- Which relation is used to turn the divider frequency into c_s, and is it valid when M < N_side m (N = 1000/1600)?
- What is the intended benchmark above η = 0.70?

---

## 7. Sources read (verbatim paths under HardDisks/hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/)

- campaign_r25_psi6_20260823/eta_0p720/00_COMMAND.md, speed_of_sound_batch_status.json, run.log; eta_0p630/0p700/0p720/0p750/speed_of_sound_psi6.csv
- routeB_radius_N100_L0_20_20260825/eta_0p720/00_COMMAND.md, speed_of_sound_batch_status.json, run.log; eta_0p630/0p700/0p720/0p750/speed_of_sound_psi6.csv; eta_0p750/speed_of_sound_batch_status.json, run.log
- ladder_N400_20260825/eta_0p720/00_COMMAND.md, run.log; ladder_N100/N200/N400 eta_0p720 and N400 eta_0p700 speed_of_sound_psi6.csv; ladder_N400_20260825/finite_size_scaling_cs_vs_N.png
- overnight_N1000_20260826/OVERNIGHT.log, finite_size_scaling_with_N1000.png, FINAL speed_of_sound_with_finite_size_correction.png
- validate_acc_N100_20260826/eta_0p700/00_COMMAND.md, speed_of_sound_batch_status.json
- finitesize_aspect_20260826/CAMPAIGN.log, finite_size_fixed_aspect_ratio.png, FINAL cs_vs_eta_fixed_aspect_ratio.png
- 0000_PLAN_OVERALL/ALL_MARKDOWNS/260908_speedsound_validation_low_high_density.md (all 1443 lines)
- Not readable from here (too deep): all analysis*/final_plots/*.csv and *.pdf, routeA_vs_routeB_finite_size.png, cs_and_psi6_vs_eta.png, famA/famC per-eta results, overnight an/ summaries.

EOS numbers in 2b: c_s² /(kT/m) = Z + ηZ′ + Z² with Z from Henderson (1+η²/8)/(1−η)² and from Santos–López de Haro–Yuste; the two differ by ≤1.5% up to η=0.63, which is the size of the effect being claimed, so quote which EOS is used.
