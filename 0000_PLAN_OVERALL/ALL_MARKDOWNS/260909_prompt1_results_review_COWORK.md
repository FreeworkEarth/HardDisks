# 260909 — Prompt 1 results: what is established, where the write-up overreaches, what to run next (COWORK)

Written 2026-09-09 by Cowork, from the files on disk, not from CC's report: `00_pressure_validation_20260907/run.log` (67 PROD, 17 DISCARD, 1 `valid=0`, 345 warning lines), `chunk_calibration*.csv`, `analysis/tables_ABC.md`, `analysis/tables_D.md`, `analysis/dev_vs_invsqrtN.png`, `ladder_0p69_N900/{run_ladder.sh, blk_eq*_s0.csv, traj_eq*_s0.csv}`, `edmd_core/edmd.c` lines 1383–1499, and `ALL_MARKDOWNS/260908_pressure_final_analysis.md`. Where I quote CC's numbers I checked them against those files; where I disagree it is with the interpretation, not the arithmetic.

---

## 1. Verdict

- The campaign is complete and clean. 162 accepted trajectories (96 preserved + 66 dense), 18 discards all accounted for, 21/21 reruns clean at the capped chunks, zero new discards, ladder 9/9 with health 0. The runner-level fix did what it was supposed to. Nothing in the core was touched. Good.
- The low-density validation is publication-grade: Z∞ within ±0.17 % of Kolafa–Rottner for η = 0.2–0.5 and +0.13 % at 0.60, the wall route landing on the same limit within 0.25 %, the 1/√N fits with χ² of order 1.
- Above η = 0.65 the write-up overreaches in three places — the "16× equilibration ladder", "systematic departure that grows with N", and "threshold exactly 320/N" — and the event-budget story has an accounting hole that a referee who knows EDMD will find in five minutes. Section 3 and 4.
- Proposed claim range in section 5, one more run (N = 2500) in section 6, and Prompt 1b for CC in section 7.

---

## 2. Checked against disk

| CC's statement | File | Result |
|---|---|---|
| 17 DISCARD + 1 production failure, 345 warnings, all `dominant=AB`, 344/345 in equilibration | `run.log` | confirmed verbatim (`grep -c`); 1 warning has `stagnant≠0`, the rest 0 |
| 162 accepted, 21/21 target cells | `run.log`, `runs/` | 67 PROD lines = 48 original + 1 failed + 18 rerun; 96 + 66 = 162 |
| KR values 6.4553 / 8.4080 / 9.3625 / 10.2007 | recomputed from the `KR2006_COEFFICIENTS` in `plot_speed_of_sound_edmd.py` | 6.4553 / 8.4080 / 9.3624 / 10.1969 (0.04 % at 0.69, irrelevant); KR peaks at 0.700 and falls to 9.42 at 0.71 — "invalid above ≈ 0.695" is right |
| ladder eq=400 reproduces the accepted 0.69/900 runs within −0.4 to −1.0 σ | `run.log`, `traj_eq400_s*.csv` | confirmed (different chunk, 0.3 vs 0.4, so different trajectories; agreement is statistical, as it should be) |
| seed means 9.9026 → 9.9298 → 9.9327 | `tables_D.md` D1 | confirmed, but see 3(a): these are three windows of the same trajectory |
| Z_pair,∞ table and wall route | `tables_ABC.md` | confirmed; χ² per fit added below, which the tables do not report |
| physical collision rate | `traj_eq400_s0.csv`: `pair_event_count` = 5 506 936 over 600 time units, N = 900 | 10.2 pair collisions per particle per unit time; Enskog with kT = m = σ = 1, Γ/2 = 2(Z−1)/√π = 10.4. The simulation's collision rate is textbook. This number is what breaks the event-budget story (section 4) |

---

## 3. Three corrections to the framing

### (a) The ladder is one trajectory per seed, read at three windows

`run_ladder.sh` launches the three rungs with the same seed and the same chunk (0.3). The runner is deterministic, and the accumulator reset at the end of equilibration does not touch the dynamics. I compared the block files: every one of the 20 equilibration + 30 measurement blocks of `blk_eq400_s0.csv` is bit-identical (T, Z_pair, Z_wall_x, Z_wall_y, ψ₆) to the block at the same `t_start` in `blk_eq6400_s0.csv`; likewise all 110 blocks of `blk_eq1600_s0.csv`. So "equil = 400 / 1600 / 6400" are the windows t ∈ [400,1000], [1600,2200], [6400,7000] of one 7000-unit trajectory per seed.

That is a *better* result than "16× longer equilibration does not move Z", and it should be written that way: **at η = 0.69, N = 900, Z_pair is stationary over t = 400–7000 in three independent trajectories (9.90–9.96, seed-to-seed sd ≈ 0.02, one trajectory (k = 1) wandering up by 0.065 ≈ 3 block-σ between the first and second window and staying there)**. What it does *not* show is that a differently equilibrated system (a different route to 0.69, e.g. compression from an equilibrated 0.65 state) lands on the same Z — that test was correctly skipped (no state save/load in the core), and the sentence "robust to a 16× equilibration ladder" should not be used.

### (b) "Onset of a systematic negative departure that grows with N" is not what the data say at 0.67

Per-cell means (`tables_ABC.md`, table B) and the 1/√N fits with the χ² that the tables omit (2 parameters, 3 points, 1 degree of freedom; errors as CC used them, max(block sem, seed sem)):

| η | Z(400) | Z(900) | Z(1600) | Z∞ (3-pt) | vs KR | slope a | χ²₁ | Z∞ from 900/1600 only |
|---|---|---|---|---|---|---|---|---|
| 0.60 | 6.5941 | 6.5475 | 6.5292 | 6.4635 ± 0.0056 | +0.13 % | +2.60 | 0.67 | 6.4743 |
| 0.65 | 8.4710 | 8.4373 | 8.4221 | 8.3727 ± 0.0075 | −0.42 % | +1.96 | 0.07 | 8.3765 |
| 0.67 | 9.2867 | 9.2246 | 9.2340 | 9.1771 ± 0.0132 | −1.98 % | +2.07 | **7.6** | 9.2622 (slope −1.1) |
| 0.69 | 10.0850 | 9.9255 | 9.8835 | 9.6476 ± 0.0286 | −5.42 % | **+8.69** | 2.3 | 9.7575 (slope +5.0) |

- At 0.67 the N = 1600 mean (9.2340 ± 0.0053) sits *above* the N = 900 mean (9.2246 ± 0.0087). The three points are not on a Z∞ + a/√N line (χ²₁ = 7.6, p ≈ 0.006); drop N = 400 and the slope changes sign. The "−1.98 %" is the intercept of a line through inconsistent points. Look at `dev_vs_invsqrtN.png`: the 0.67 points visibly do not follow their line.
- At 0.69 the fit is tolerable (χ²₁ = 2.3) but the slope is four times the value at every other density, and the intercept moves by 1.1 % depending on whether N = 400 is included. With ξ ≈ 30–50 σ and boxes of 21–43 σ, the finite-size dependence has no reason to be a pure perimeter term here.
- At 0.65 the form holds cleanly (χ²₁ = 0.07, a = 1.96 in line with 1.46 at 0.5 and 2.60 at 0.6) and Z∞ is 0.42 ± 0.09 % below KR. That is a 4.7 σ statistical statement. Whether it is physics (walls) or the extrapolation (a 1/N term a 3-point fit cannot see) or the reference (KR's own accuracy at 0.65) is open; a fourth N decides (section 6).

So the honest split is not "0.65–0.69: departure grows with N" but: **surface-term extrapolation valid through 0.65; from 0.67 the finite-size dependence in these hard-wall boxes is not a/√N and no bulk value is extracted; Z(N) is reported (1.4 % below KR at N = 900–1600 for 0.67; 2.7 % and 3.1 % at 0.69).** "Both routes agree" adds no independence here: pair virial and wall momentum flux are tied by momentum balance in the same trajectories; their agreement is a consistency check, not a second measurement.

Also: the seed-to-seed sd is 1.6–2.2 × the block error at 0.67/900, 0.69/900 and 0.69/1600 (0.017 / 0.023 / 0.030 against 0.010 / 0.014 / 0.014) versus ≈ 1 or below at 0.65 — trajectories at ≥ 0.67 differ by more than their internal noise, i.e. they sit in slightly different structural states for hundreds of time units; the k = 1 ladder trajectory, which moved by 0.065 between its first and second window and stayed there, is the direct example. That is the same "mechanically stationary, structurally wandering" behaviour as in the c_s campaign and belongs in the boundary section.

### (c) "The failure pattern confirms the threshold is exactly 320/N" — no

From `run.log` and the calibration tables, in units of N·chunk:

- η = 0.60: N = 1600 at chunk 0.8 (1280) exceeded the budget in only 27 of 500 equilibration calls (5.4 %). The boundary at 0.60 is ≈ 1280/N, four times the cap.
- η = 0.65–0.67: 640–720 exceed in 3.5–9 % of calls; 450 (calibration 0.5 at N = 900) verified clean over 20 calls. Boundary ≈ 600–700/N.
- η = 0.69–0.72: 360 (0.4 at N = 900) ran clean in every accepted trajectory (≈ 10 000 calls at 0.69 alone); 270 (the ladder, 0.3 at N = 900) clean over ≈ 100 000 calls. No failure data above 360 at these densities.

320/N is a conservative cap with a margin between ≈ 1.1 (at 0.69–0.72, unmeasured) and ≈ 4 (at 0.60). Fine as a rule; wrong as a "threshold", and the methods text should say the safe window shrinks ≈ 1/N *at fixed density and steeply with density*.

---

## 4. The event budget: what is actually being counted (this changes the methods paragraph)

`edmd_advance_to()` (edmd.c 1384–1499): `events_processed++` at line 1403 counts **every popped heap entry**, before the entry is validated against the collision counters (lines 1465–1477); stale entries are counted and then skipped. Every pop, valid or not, also moves all N particles to the entry's time (1454–1463). So the 250 000 "events" per call are calendar pops, not collisions.

Put the physical collision rate next to it. At η = 0.69, N = 900: 10.2 collisions per particle per unit time (from `pair_event_count`; Enskog 10.4). At the cap N·chunk = 320 a call contains 320 × (Γ/2) ≈ 2 000 (η = 0.60) to 3 300 (η = 0.69) collisions — **2–4 collisions per particle**. At the observed boundary (N·chunk ≈ 1280 at 0.60, ≈ 650 at 0.65–0.67) a call contains ≈ 5 500–8 000 collisions, so the budget trips at **≈ 30–45 calendar pops per physical collision**. That is the number the methods section has to state, because "2.5 × 10⁵ events per call" will be read as collisions, and a referee who computes N·Γ·chunk/2 gets 3 000 and asks what the other 247 000 are.

The exhaustions are throughput, not avalanches, and the log proves it: all 345 logged exceedances (the first 20 per trajectory — the warning is rate-limited at line 1410, so the totals are the `fa` counters, 27–45) occur in the last 0.0001–0.06 time units of their chunk (median 0.009, i.e. the last 1–8 % of a 0.4 or 0.8 call), never mid-chunk; `stagnant = 0` in 344 of 345. The count crosses 250 000 near the end of a call whose total is a few percent over the budget, in 3–9 % of the calls of a cell that sits at the boundary. Spread over t = 7–423 because the pop rate of the equilibrated fluid is what it is. The word "avalanche" in the warning text is the code's, not a description of what happened.

Consequences:

1. Methods: "the integrator caps the number of calendar entries processed per advance call at 2.5 × 10⁵, including invalidated entries; at η = 0.60–0.69 this corresponds to ≈ 30–45 entries per collision, i.e. 2–4 collisions per particle per call at the window we used". Then the 1/N rule follows and is not mysterious.
2. The clean core fix, after the campaign: count only validated events, or scale the budget with N. One line plus a regression test; it does not change any collision physics. Until then the runner cap stands.
3. `recalibrate6.sh` sends the calibration's stderr to `/dev/null`. The equilibrated calibration returned 0.125 at 0.67/1600 (so 0.25 = N·chunk 400 failed its 20-chunk verification) while 0.69/1600 ran clean at 0.2 (320) in production and 0.65/1600 passed at 0.25. Whether that failure was a budget exhaustion, an overlap repair, or something else is now unrecoverable. Log calibration stderr next time; the calibration table is going into the paper.

---

## 5. Proposed claim range (replaces the three bullets at the top of `260908_pressure_final_analysis.md`)

- **η ≤ 0.65: equation of state validated in the N → ∞ limit.** Z(N) = Z∞ + a/√N over N = 400/900/1600 (χ²₁ ≤ 1 at every η); Z_pair,∞ within ±0.17 % of Kolafa–Rottner 2006 for η = 0.2–0.5, +0.13 ± 0.09 % at 0.60, −0.42 ± 0.09 % at 0.65; the wall-momentum-flux route extrapolates to the same limit within 0.25 % (momentum-balance consistency). The 0.65 value is statistically below KR and awaits the N = 2500 point before it is called either way.
- **η = 0.67–0.69: no bulk value from these boxes.** The finite-size dependence is no longer a perimeter term (χ²₁ = 7.6 at 0.67 with Z(1600) > Z(900); slope 4× larger at 0.69). We report Z(N): −0.8/−1.5/−1.4 % vs KR at N = 400/900/1600 for 0.67 and −1.1/−2.7/−3.1 % for 0.69; stationary over 400–7000 time units at 0.69/900; seed-to-seed scatter up to 2 × the block error. Consistent with correlation lengths comparable to the box (Bernard–Krauth: ξ ≈ 50 σ at 0.698) in hard-wall geometry.
- **η ≥ 0.70: exploratory.** KR is not a reference there (its fit turns over at 0.70). βPσ² at N = 900 (9.23 / 9.59 / 10.14 at 0.702 / 0.710 / 0.720) sits above the coexistence plateau (9.185) and rises: a homogeneous state that cannot phase-separate in a 32 σ hard-wall box; ψ₆ drifts (|t| up to 6) while Z and T are flat; one cell (0.720/900) shows a 3.5 % x–y anisotropy in the wall pressure, which is a physical signal of an anisotropic structure, not noise.

Reference above 0.65: quote Engel et al. 2013 / Bernard–Krauth 2011 liquid-branch pressures next to KR (KR is a fit at the edge of its range there). That is a question for GPT, who knows those tables.

---

## 6. What to run, and what not to

**Run: TASK 7, N = 2500 at η = 0.65, 0.67, 0.69, 3 seeds (seed rule unchanged: 20260907 + 104729·k + 7919·N).** Chunk from the equilibrated calibration × 0.6, capped at 0.75 × 320/N = 0.096 (the extra 0.75 because at 0.69 the cap has no measured margin). Same 400 + 30 × 20 protocol. Cost: the default core is O(N) per pop and pops ∝ N, so ≈ (2500/1600)² ≈ 2.4 × the N = 1600 wall time — about 2–3 h per trajectory, 9 jobs in parallel, an evening. Decides: with four sizes the 1/√N fit has 2 degrees of freedom — does 0.65 stay at −0.4 %; does 0.67 stay non-monotone; does the 0.69 slope stay at 8.7. If the fits with four points are consistent, add the 2500-point Z∞ to the table; if not, section 5 stands as written.

**Do not run an aspect-ratio series at fixed N.** I considered it: at N = 900 a 4:1 box changes perimeter/area by only 25 % relative to the square (2:1 by 6 %), so the surface term shifts by ≈ 0.02 in Z at 0.67 — two sigma, and a 16:1 strip that would double P/A is 8 σ wide and brings the route-B confinement physics back. N is the lever.

**Do not run more equilibration.** Settled by 3(a).

**Later, core, after the campaign (TASK 9, not now):** budget counts validated events only, or scales with N; expose a pop counter so the calibration can measure pops-per-collision directly instead of inferring it.

---

## 7. Prompt 1b for Claude Code (paste as is)

```
Context: same repo and campaign as Prompt 1. Read
0000_PLAN_OVERALL/ALL_MARKDOWNS/260909_prompt1_results_review_COWORK.md first; every
number below is in it with its file. Rules unchanged: no core changes, no edits to
accepted trajectories, no commits, quote numbers from files.

TASK A -- corrections to validation/write_pressure_final_analysis.py and the two analysis
scripts (regenerate the document afterwards):
1. Claim range: replace the three bullets with the text of section 5 of the review note,
   with the numbers parsed from tables_ABC.md as now.
2. Table B: add a chi2 column (weighted, 2 parameters, 3 points -> 1 dof) and a column
   "Zinf from N=900/1600 only". Print p-value for chi2_1 in a footnote. Flag any eta with
   chi2_1 > 4 as "1/sqrt(N) form rejected".
3. D1: add the block-identity check (compare blk_eq400_s*.csv and blk_eq1600_s*.csv with
   blk_eq6400_s*.csv at equal t_start; report identical/different counts per seed) and
   rewrite the D1 text as "one trajectory per seed, three windows; Z stationary over
   t = 400-7000". Remove "16x longer equilibration".
4. TASK 1 text: replace "threshold ... = 320/N" by the boundary table of section 3(c)
   (N*chunk 1280 at 0.60, 600-700 at 0.65-0.67, >= 360 at 0.69-0.72) and call 320/N a
   conservative cap with 1.1-4x margin.
5. Methods paragraph: rewrite per section 4 -- "calendar entries including invalidated
   ones" (edmd.c line 1403 vs 1465-1477), 30-45 entries per collision, 2-4 collisions per
   particle per call at the cap, exhaustion in the last 1-8% of the call (all 345 logged
   cases), "in production" -> "during equilibration". Keep the sixth-failure-class framing;
   drop the word avalanche except when quoting the warning text.
6. Add to section C a one-line note that pair and wall routes are tied by momentum balance
   and their agreement is a consistency check, not an independent measurement.
7. recalibrate6.sh and run_pressure_campaign2.sh: send calibration stderr to a log file
   next to the calibration table (calib_<eta>_<N>_<seed>.err), never to /dev/null.

TASK B -- TASK 7 runs (launch only after TASK A is regenerated and you have shown me the
new claim-range text):
8. Calibrate (equilibrated method, seeds 911000001/2) and run N = 2500 at eta = 0.65,
   0.67, 0.69, 3 seeds each, seed rule unchanged, chunk = min(0.6*min(A,B), 0.75*320/N).
   Same 400 + 30x20 protocol, same health contract, resumable driver. Report the
   calibration rows and the launch state; do not start anything else on the machine.
9. When they land: rerun analyze_pressure_campaign.py, report table B with four sizes and
   the chi2 (2 dof now) per eta, and regenerate the document. If the four-point fit at
   0.65 still gives Zinf below KR by more than 3 sigma, say so; do not adjust the claim
   text yourself -- I want to see the numbers first.

Report actual commands and outputs, git status --short, git diff --stat.
```

---

## 8. For GPT (paste)

```
Pressure campaign done: 162 accepted, 18 discards explained (per-call calendar-pop budget
250k incl. stale pops = 30-45 pops/collision, i.e. 2-4 collisions per particle per call;
throughput, not avalanches -- all exhaustions in the last 1-8% of a chunk). Claim: EOS
validated to eta=0.65 (Zinf within 0.2% of KR to 0.60, -0.42+-0.09% at 0.65, chi2 ok);
at 0.67 the 1/sqrt(N) form fails (chi2_1=7.6, Z(1600)>Z(900)) and at 0.69 the slope
quadruples, so we report Z(N), no bulk value. Question for you: best reference pressures
for the liquid branch at eta=0.65-0.69 (Engel 2013 / Bernard-Krauth tables) and KR2006's
stated accuracy there; and whether a 1/N corner term of the size needed to explain -0.4%
at 0.65 is plausible for hard walls.
```
