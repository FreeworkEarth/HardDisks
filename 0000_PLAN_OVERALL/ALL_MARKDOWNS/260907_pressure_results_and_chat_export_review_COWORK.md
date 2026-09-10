# Pressure campaign: actual results on disk, plus what the full chat export adds

Date: 2026-09-07 (Cowork / Claude). Sources: the ChatGPT export `harddisk speed sound Usage limit analysis 20260907.md` (26,222 lines, Aug 21 → Sep 7 08:06 HST), and the campaign folder it names: `hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_pressure_validation_20260907/` (`run.log`, `chunk_calibration.csv`, `preserved/`, `runs/`). Numbers below are computed from those files, not quoted from chat.

This supersedes the "Z: not checked" parts of the three earlier notes; the three notes were also corrected in place (accelerated core; "three estimators" wording).

---

## 1. Corrections the export forced on my earlier notes

- **The accelerated core failed validation.** `validate_acc_N100_20260826`: 269 valid / 181 invalid of 450 (90 boundary escapes, 91 particle overlaps) on settings where the default core gives 450/450. I had written "ran clean 90/90" after looking at one η's status file. Wrong; corrected. Nothing scientific used it; the N=1000 overnight ran on the default core at ~405 s per trajectory.
- **"Three independent estimators" is overstated** (GPT caught it on Sep 6, CC accepted): Z_pair (pair virial) and Z_wall (wall momentum flux) are the two independent routes; Z_wall,x vs Z_wall,y is an isotropy check of the wall route. Corrected in my notes.
- **Error bars exist for c_s.** E.g. c_s(η=0.63, N=1000) = 9.622 ± 0.111 (Liu-based 9.488, +1.4%). My "the table has no uncertainties" applied to the chat table, not the analysis.
- **GPT already made the strip-vs-square and "c_eff is not c_s" points on Aug 26**, and CC ran the fixed-aspect campaign in response. My review agrees with that record; it is not new.
- **Engel et al. numbers that matter for the transition section:** orientational correlation length ≈ 50 σ at η = 0.698, positional ≈ 100 σ at 0.718, and finite-size effects exceed statistical errors even at N = 1024². The fixed-aspect boxes are 20–40 σ across, so "plateau at N ≥ 900" near 0.70 is not convergence to bulk; it cannot be, by the literature's own yardstick.
- **Aug 26 findings I had not listed** (all reproducible with a seed/command, all worth the methods section): traces ~290× oversampled (1.3 M rows, 210 MB per low-η trajectory, campaign would have needed 165 GB; the disk was full — 1.8 TB of 1.8 TB); trace filenames encode L₀ as `(int)(L0*10)`, so η = 0.715 and 0.720 silently overwrote each other, which is why no fine scan near the transition existed before; FFT search-window clipping; the tunnelling bug's three-link cause.
- **The pressure OUTDIR** is `…/mode1_normalized_units/00_pressure_validation_20260907` (a sibling of `00_eta_sweep_ROMAN`, which is why I did not find it).

---

## 2. Where the pressure campaign actually stands (from `run.log`, last write ≈ 8 h before this note)

- **Accepted: 144 trajectories.** 96 at η ≤ 0.5 (in `preserved/`) + 48 dense (in `runs/`).
- **Discarded: 18**, all dense, all N ≥ 900: 17 with health events *during equilibration* (forced advances 27–45 per trajectory, followed by 1,776–18,172 overlap repairs) and 1 with health events during production (η = 0.60, N = 900). The contract did exactly what it should — every one of those would have produced a plausible Z.
- **Empty cells as a result:** η = 0.60/N = 1600, 0.65/900, 0.65/1600, 0.67/900, 0.67/1600. So the dense N-scaling exists only at η = 0.60 (N = 400, 900) and η = 0.69 (N = 400, 900, 1600). The diagnostic η = 0.702/0.710/0.720 ran at N = 400 and 900 only (by design).
- CC's last status in the export ("5 accepted, 27 running", Sep 7 08:06) predates all of this; neither CC nor GPT has seen the discards yet.

### Why the discards happened (this is a real finding, not bad luck)

`edmd_core/edmd.c` line 9: `#define EDMD_ADVANCE_MAX_EVENTS 250000L`, and `edmd_advance_to()` counts events per *call* and forces a free-flight advance when the count exceeds it. That budget does not scale with N. So the largest safe chunk scales like 1/N at fixed η, and the campaign's empirical safe values say exactly that: chunk 0.8 works at N = 400, 0.4 at N = 900 (η = 0.69 passed at 0.4, failed at 0.8), 0.2 at N = 1600 (passed at 0.2, failed at 0.4). Rule of thumb from the data: chunk ≤ ~320/N at η ≥ 0.6.

The two-seed calibration did not catch it because it verifies a candidate for 5 time units on a *freshly seeded lattice*, whose event rate is lower than the equilibrated fluid's; the bursts that exceed the budget appear later. [Corrected 260908: I first wrote "the warnings sit at t ≈ 7–45", from the first discards I opened. CC's scan of all warnings says n=345, t from 7.15 to 423.20, median 100.39, 344/345 in equilibration — the bursts are spread over the whole equilibration, so they are a property of the dense fluid at those chunk sizes, not a lattice-start transient. That makes the chunk cap ∝ 1/N the essential fix; calibrating on an equilibrated instance helps but cannot replace it.] GPT's warning that "calibration is a heuristic, production decides" was right; the fix is either (a) calibrate on an equilibrated disposable instance and verify over ≥ 20 chunks, (b) hard-code chunk ∝ 1/N with the table above, or (c) make the budget ∝ N in the core (it exists to catch avalanches, not to cap throughput; 250 k events is only ~156 per particle for N = 1600). (c) is the clean one and is a one-line change with a regression test; (b) is the safe one for a relaunch tonight.

---

## 3. Z results, low density (η ≤ 0.5): this part is done and it is good

Z_pair mean over seeds, vs Kolafa–Rottner 2006 (coefficients from `plot_speed_of_sound_edmd.py`), and the wall − pair gap:

| η | N=400 dev | N=900 dev | N=1600 dev | 1/√N-extrapolated Z∞ | KR2006 | Z∞ vs KR | wall−pair gap (400/900/1600) |
|---|---|---|---|---|---|---|---|
| 0.20 | +0.53% | +0.44% | +0.32% | 1.5723 | 1.5704 | +0.13% | 2.6 / 1.7 / 1.4 % |
| 0.30 | +1.15% | +0.87% | +0.52% | 2.0630 | 2.0633 | −0.01% | 3.2 / 2.1 / 1.6 % |
| 0.40 | +1.56% | +1.04% | +0.76% | 2.8285 | 2.8293 | −0.03% | 3.7 / 2.4 / 1.8 % |
| 0.50 | +1.93% | +1.39% | +1.05% | 4.1146 | 4.1064 | +0.20% | 4.2 / 2.7 / 2.0 % |

(η ≤ 0.1: all deviations ≤ 0.16%; the ideal-gas limit at η = 0.005 is hit to 0.01%.) Seed-to-seed scatter is at or below the block error (e.g. η = 0.5, N = 1600: sd 0.0004, sem 0.0020). Zero health events in all 96.

Two things the table says that the chat had not quite said:

1. **The pair virial in a hard-wall box also carries a positive ~1/√N surface term**, about half the size of the wall route's. So "Z_pair is the bulk route" is not right either; both routes need the extrapolation. Done that way, the N → ∞ pair virial lands on Kolafa–Rottner within 0.2% for η = 0.2–0.5. That is the validation statement for Paper 1, and it is a strong one.
2. The wall − pair gap scales as 1/√N at every η (ratios 1 : 0.66 : 0.49 vs 1 : 0.67 : 0.50), and the two routes extrapolate to the same limit (wall route N → ∞: +0.41 / +0.18 / −0.13 / −0.19 / +0.06 % vs KR at η = 0.1–0.5). That is the "independent mechanical cross-check" GPT wanted, quantified.

---

## 4. Z results, dense grid (η ≥ 0.6): the validated range should stop around 0.65

Accepted cells, Z_pair vs KR (KR is only meaningful to ≈ 0.69; it turns over above 0.695 and is nonsense at 0.71–0.72):

| η | N | seeds | Z_pair | dev vs KR | wall−pair | chunk |
|---|---|---|---|---|---|---|
| 0.60 | 400 | 5 | 6.594 | +2.15% | +4.6% | 0.8 |
| 0.60 | 900 | 3 | 6.544 | +1.38% | +3.0% | 0.8 |
| 0.65 | 400 | 5 | 8.471 | +0.75% | +4.8% | 0.8 |
| 0.67 | 400 | 5 | 9.287 | −0.81% | +4.8% | 0.8 |
| 0.69 | 400 | 5 | 10.085 | −1.10% | +4.9% | 0.8 |
| 0.69 | 900 | 4 | 9.926 | −2.66% | +3.2% | 0.4 |
| 0.69 | 1600 | 3 | 9.884 | −3.07% | +2.4% | 0.2 |
| 0.702 | 400 / 900 | 3 / 3 | 10.549 / 10.325 | — | +5.0 / +3.3% | 0.8 / 0.4 |
| 0.710 | 400 / 900 | 3 / 3 | 10.945 / 10.613 | — | +5.0 / +3.3% | 0.8 / 0.4 |
| 0.720 | 400 / 900 | 3 / 3 | 11.486 / 11.065 | — | +5.0 / +3.3% | 0.8 / 0.4 |

Reading it:

- **η = 0.60 behaves like the fluid branch**: +2.15% → +1.38% is exactly the 1/√N surface term seen at 0.5, heading to ≈ 0 at N → ∞. Fine.
- **From η = 0.65 upward the N = 400 values fall below the surface-term trend** (they should sit ≈ +2% above KR at N = 400; they sit at +0.75, −0.81, −1.10%). At η = 0.69 the deviation *grows* with N (−1.1 → −2.7 → −3.1%), the opposite of a boundary correction. Both routes extrapolate to the same limit: Z∞ ≈ 9.67 (pair) and ≈ 9.64 (wall), i.e. βPσ² ≈ 8.48, about 5% below Kolafa–Rottner (10.20; βPσ² = 8.96) and below what the Bernard–Krauth/Engel EOS gives near 0.69 (their coexistence plateau is βPσ² = 9.185 at 0.700–0.716).
- Within each η = 0.69 trajectory the 30 blocks show **no drift** in Z, T or ψ₆ (Z blocks fluctuate ±1% with no trend; global ψ₆ ≈ 0.1–0.2, local ≈ 0.71, stable). Global ψ₆ falls as 1/√N (0.18 → 0.12 → 0.10) while local stays 0.71: many small six-fold domains, a dense liquid's structure, not a crystal. So the state is *stationary*, and its pressure is 5% low. Candidate explanations, all testable: (a) 400 time units of equilibration from a square-lattice seed is not enough this close to the transition even though the observables look flat (the "mechanically stationary, structurally wandering" case GPT flagged); (b) hard walls change the bulk state at this density (pre-ordering), which would need a periodic-boundary or much larger control; (c) the reference: KR is at the very edge of its fit range, but the Engel plateau value argues the reference is not 5% off.
- The diagnostic points 0.702–0.720 give βPσ² = 9.23 / 9.59 / 10.14 at N = 900, above the bulk coexistence plateau and still rising: a homogeneous state that cannot phase-separate in a 32 σ hard-wall box. Exploratory, as labelled.

**Consequence for Paper 1:** claim the EOS validation for η ≤ 0.60 with N → ∞ extrapolation (≤ 0.2% for 0.2–0.5, ≈ 0 at 0.6), report 0.65–0.69 as the onset of a systematic negative departure that grows with N, and treat ≥ 0.70 as the structural boundary. That is a cleaner and more defensible line than "validated to 0.69", and it matches what the c_s side says (fluid branch solid to ≈ 0.63–0.65, transition region a finite-size study).

---

## 5. Small things worth knowing

- **Temperature is not 1.0 per trajectory in the pressure runner**: T_mean runs 0.913–1.057 between seeds at N = 400 (velocities are drawn from a unit Gaussian and never rescaled; the relative scatter is 1/√N). Harmless for Z (it is T-independent) but the same seeding logic at N = 100 in the speed-of-sound runs would give ±10% in T and ±5% in c_s per trajectory. `--kbt1` presumably rescales; confirm that each speed-of-sound trajectory reports kT = 1.000 after the hold, or that the analysis divides by √T_measured.
- The lattice seeder starts with ψ₆ ≈ 0 (it is 4-fold), so the pressure runs do not start ordered; the question at 0.69 is whether 400 units is enough to *reach* the fluid, not whether they start in a crystal.
- The chunk table CC reported at 08:06 (0.8/0.4/0.2 pattern, "monotone, no discontinuity") was the calibration's opinion; production disagreed at six cells. Keep the table in the methods section together with the discards — that is the honest version of "calibration is a heuristic".

---

## 6. Paste-able for Claude Code

```
Pressure campaign status (from run.log in 00_pressure_validation_20260907):
144 accepted (96 low + 48 dense), 18 discarded. All discards are N >= 900 at
eta = 0.60-0.67: 17 with forced advances during equilibration (chunk 0.8 at
N=900/1600 and 0.4 at N=1600), 1 with health during production (0.60/900).
Empty cells: 0.60/1600, 0.65/900, 0.65/1600, 0.67/900, 0.67/1600.

1. Root cause: EDMD_ADVANCE_MAX_EVENTS = 250000 per call is independent of N,
   so the safe chunk scales ~1/N (empirically 0.8 / 0.4 / 0.2 at N = 400 / 900 /
   1600 for eta >= 0.6). The 5-unit calibration on a fresh lattice does not see
   the equilibrated event rate. Choose ONE fix and document it:
   (a) production chunk = min(1, 320/N) for eta >= 0.6 (safe, no core change), or
   (b) budget proportional to N in the core (cleaner; add a regression test that
       a chunk which passed at N=400 also passes at N=1600 at the same eta).
   Then rerun ONLY the five empty cells with the same seeds. Do not touch the
   accepted trajectories.

2. Report Z_pair AND Z_wall extrapolated in 1/sqrt(N) per eta (both routes carry
   a surface term; the pair route's is ~half the wall route's). For eta <= 0.5
   the pair extrapolation lands on Kolafa-Rottner within 0.2%; state that as the
   validation result. Plot dev-vs-1/sqrt(N) per eta with the fit line.

3. eta = 0.69 anomaly: Z_pair dev vs KR is -1.1 / -2.7 / -3.1 % at N = 400/900/1600,
   growing with N; both routes extrapolate to Z ~ 9.67 (betaP sigma^2 ~ 8.49),
   ~5% below KR and below the Engel coexistence plateau (9.185). Blocks are
   stationary (Z, T, psi6). Before interpreting: equilibration ladder at 0.69,
   N = 900, equil = 400 / 1600 / 6400 with Z and both psi6 measures logged during
   equilibration; and one run seeded from an equilibrated eta = 0.65 fluid
   configuration compressed to 0.69 (if the runner can load a state). If Z
   moves, the current 0.65-0.69 numbers are not equilibrated; if not, the
   departure is real and belongs in the boundary section.

4. Claim range: EOS validated for eta <= 0.60 (N -> inf extrapolated); 0.65-0.69
   reported as the onset of a systematic departure; >= 0.70 exploratory.

5. Confirm kT per trajectory in the speed-of-sound runs (should be exactly 1
   after --kbt1 and the hold), or that c_s is normalised by sqrt(T_measured).
   In the pressure runner T_mean varies 0.91-1.06 between seeds; fine for Z,
   not for c_s if the same happens there.

6. Record: accelerated core (--edmd-acc=1) = 181/450 invalid, excluded; the
   forced-advance discards above are the sixth failure class for the methods
   section (budget exhaustion during equilibration masked by a short calibration).
Do not auto-commit.
```
