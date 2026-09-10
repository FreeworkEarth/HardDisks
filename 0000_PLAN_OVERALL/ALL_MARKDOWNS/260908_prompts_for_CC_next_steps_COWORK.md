# Prompts to continue: three for Claude Code (in order), one for GPT

Date: 2026-09-08 (Cowork / Claude). Use them one at a time; each ends with a report that the next one needs. Start every CC session with the "read first" line so it has the map (the notes hold the map, the disk holds the territory).

---

## Prompt 1 — Pressure: close the dense grid and settle η = 0.69 (Paper 1, Z section)

```
Read first, do not skip:
  0000_PLAN_OVERALL/ALL_MARKDOWNS/260907_pressure_results_and_chat_export_review_COWORK.md
  0000_PLAN_OVERALL/ALL_MARKDOWNS/260907_paper_plan_equilibrium_then_nonequilibrium_COWORK.md (section 0-1)
Campaign folder:
  hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_pressure_validation_20260907/

State (verify against run.log before doing anything): 144 accepted (96 low-eta in
preserved/, 48 dense in runs/), 18 discarded, all at N >= 900 for eta = 0.60-0.67:
17 with forced advances during equilibration (chunk 0.8 at N=900/1600, 0.4 at
N=1600) and 1 with health during production (0.60/900). Empty cells: 0.60/1600,
0.65/900, 0.65/1600, 0.67/900, 0.67/1600.

Cause: EDMD_ADVANCE_MAX_EVENTS = 250000 per edmd_advance_to() call is independent
of N, so the safe chunk scales ~1/N; empirically 0.8 / 0.4 / 0.2 at N = 400 /
900 / 1600 for eta >= 0.6. The 5-unit calibration on a freshly seeded lattice sees
a lower event rate than the equilibrated fluid, so it passed chunks that fail
later in equilibration.

TASK 1 - confirm the cause. From run.log, list every DISCARD with (eta, N, chunk,
fa, orep, crep) and the t values of the avalanche warnings. Confirm that no
accepted trajectory has any nonzero health counter. Report the table verbatim.

TASK 2 - fix at the RUNNER level, not in the core (no collision physics, no core
change during a campaign):
  a. production chunk = min(calibrated_chunk, 320/N) for eta >= 0.6;
  b. calibration verifies the candidate on a disposable instance AFTER 100 time
     units of equilibration on that instance, over >= 20 consecutive chunks;
  c. keep the all-or-nothing production contract unchanged.
Print the resulting chunk table (eta, N, chunk) before launching.

TASK 3 - rerun ONLY the five empty cells with the SAME seeds
(seed = 20260907 + 104729*k + 7919*N). Do not touch accepted trajectories. Any
discard again -> stop and report, do not adapt.

TASK 4 - equilibration ladder at eta = 0.69, N = 900, 3 seeds each, with
equilibration_time = 400 / 1600 / 6400 and the same 30 x 20 measurement blocks.
Log per block DURING equilibration as well: T, Z_pair, Z_wall_x, Z_wall_y,
psi6_global, psi6_local. Report Z_pair(equil) per seed and the block traces.
If the runner can start from a saved configuration, add one run seeded from an
equilibrated eta = 0.65, N = 900 state compressed to 0.69; otherwise skip and say so.

TASK 5 - analysis, four parts, tables verbatim:
  A. numerical validity: accepted/discarded per cell, health counters, chunk table
     incl. the two calibration seeds, T_mean per trajectory (it varies 0.91-1.06;
     state that Z is T-independent).
  B. Z_pair vs Kolafa-Rottner per (eta, N) with block sem and seed sd; fit
     Z(N) = Z_inf + a/sqrt(N) per eta over N = 400/900/1600; report Z_inf +- sigma
     and dev vs KR. Expected: within 0.2% for eta 0.2-0.5, ~0 at 0.60.
  C. wall route: Z_wall_x vs Z_wall_y (isotropy), wall-pair gap vs 1/sqrt(N),
     Z_wall_inf per eta; both routes must agree on Z_inf.
  D. structural boundary: eta = 0.65-0.69 departure (currently -1.1/-2.7/-3.1% at
     0.69 for N = 400/900/1600, both routes -> Z_inf ~ 9.65, KR 10.20, Engel
     coexistence plateau betaP sigma^2 = 9.185), the equilibration ladder result,
     and 0.702-0.720 as exploratory (KR is invalid above ~0.695; do not quote
     deviations vs KR there).
  Plots: dev vs 1/sqrt(N) per eta with fit; Z(eta) with both routes and KR
  clipped at 0.69; stationarity panels for 0.60, 0.69, 0.72.

TASK 6 - write 0000_PLAN_OVERALL/ALL_MARKDOWNS/260908_pressure_final_analysis.md
with the tables, the claim range (EOS validated for eta <= 0.60 with N -> inf
extrapolation; 0.65-0.69 onset of a systematic departure; >= 0.70 exploratory),
and a methods paragraph for the sixth failure class: event-budget exhaustion during
equilibration masked by a short fresh-seed calibration.

Rules: do not change collision physics or the core; do not modify accepted
trajectories; do not auto-commit; quote numbers from files, not from memory.
```

---

## Prompt 2 — Speed of sound: cleanup and the transition section (Paper 1, sections 4–5)

```
Read first:
  0000_PLAN_OVERALL/ALL_MARKDOWNS/260907_speedsound_recap_and_critical_review_COWORK.md
  (sections 2 and 5 are the spec; 2c, 2d, 2f are the load-bearing ones)
Root: hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/

TASK 1 - structural stationarity audit, analysis only. For every campaign
(campaign_r25_psi6_20260823, routeB_radius_N100_L0_20_20260825, ladder_N100/N200/
N400, overnight_N1000_20260826, finitesize_aspect_20260826 famA and famC): from
speed_of_sound_psi6.csv per (eta, N, M): mean psi6_global_hold, mean
psi6_global_end, run-to-run sd of psi6_global_end, mean neighbors_hold/end.
Flag |end - hold| > 0.1, bimodal end values, and neighbours < 4.3 (square-lattice
states). Table verbatim. Compute psi4 from stored configurations for routeB
eta >= 0.72 if configurations exist; otherwise say so.

TASK 2 - per-trajectory temperature. Confirm kT after the hold is 1.000 for every
speed-of-sound trajectory (or that the analysis normalises c_s by sqrt(T_measured)).
Report min/max T over the campaigns. If T varies, c_s must be rescaled before any refit.

TASK 3 - provenance fixes:
  - overnight_N1000: find the eta=0.700 M=750 trajectory with health=1 (which
    counter), exclude it, re-run that eta's analysis, add the missing 00_COMMAND.md
    for the an/ analysis and the two FINAL PNGs;
  - mark finite_size_scaling_with_N1000.* and the c_inf = 10.44 strip extrapolation
    SUPERSEDED (README in the folder);
  - remove peak_on_search_boundary points from every final figure;
  - record that --edmd-acc=1 failed validation (181/450 invalid) and is excluded.

TASK 4 - refits, no new runs:
  - famA: fit c(N) = c_inf + a/sqrt(N) over N = 400/900/1600 per eta; report
    c_inf +- sigma; at eta = 0.67 it should approach the fluid EOS value;
  - per-mass residuals of the freq-vs-kterm fit for N >= 900; refit with only
    M >= 2*N_side*m and report the shift in c_s;
  - confirm in analyze_speed_of_sound_by_eta.py that kterm = K/(2 pi L_eff) with
    cot K = (M/(2 N_side m)) K and L_eff = L0 - 2r; cite the lines.

TASK 5 - final figure: route A + strip ladder + N=1000 as the validation for
eta <= 0.65 vs a NAMED EOS (Kolafa-Rottner, clipped at 0.69) with a residual
table; route B as a confinement series with H/sigma on a top axis; famA c_inf
points with error bars for 0.67-0.73; coexistence band shaded; Liu solid line
labelled "bulk mode, lower bound on c_L"; provenance box (N, r, L0, H, seeds,
hold). Colour points by psi6_end. Drop or grey out every cell flagged in Task 1.

TASK 6 - two cheap N=100 tests at eta = 0.72 (minutes each), report
psi6_hold, psi6_end, psi4, c_s +- sigma:
  a. seeding swap: route A geometry from a jittered rectangular start, route B
     geometry from a jittered hex start (or a stored hex configuration);
  b. hold length: --wall-hold-steps 2000 vs 20000 for route A and famA N=400 at
     eta = 0.70 and 0.72, logging psi6 during the hold.

TASK 7 - write 0000_PLAN_OVERALL/ALL_MARKDOWNS/260908_speedsound_final_analysis.md
with the audit table, the refit table, the claim range (fluid branch validated to
~0.63-0.65; 0.67-0.73 a finite-size/ordering study, no bulk c_s claimed; route B
>= 0.75 a seeded square jam), and what Task 6 decided about the wording.

Rules: no new campaigns beyond Task 6; no core changes; no auto-commit; numbers
from files.
```

---

## Prompt 3 — The divider in a calm gas (Paper 1 section 6, bridge to Paper 2)

```
Read first: 260907_paper_plan_equilibrium_then_nonequilibrium_COWORK.md, section 1
item 5 and Level 1.

From EXISTING speed-of-sound runs (route A campaign_r25_psi6_20260823 and famA),
using the hold phase (divider held) and the ring-down after release:

TASK 1 - equipartition of the released divider: 0.5 M <v_w^2> vs 0.5 kT and the
v_w histogram vs a Maxwellian, per (eta, M), from the trace CSVs. Report the ratio
with error; it must be 1 within error at every eta (this is ensemble-independent).

TASK 2 - spring divider, if runs with --spring-k exist in EDMD with the validity
ledger: position variance vs kT/k_eff, k_eff = k + 2 N_side m c_s^2/L_eff^2 with c_s
from the Paper 1 EOS. State whether the adiabatic or isothermal stiffness applies
in a microcanonical box and quote both. If no such runs exist, run 3 eta x 3 k x
10 seeds at N=100 (cheap) and say so.

TASK 3 - damping: fit the ring-down envelope of x_w(t) after release to
exp(-gamma t / 2M) per (eta, M); compare with the FFT linewidth. Report gamma(eta).

TASK 4 - Sivak-Crooks friction from the HELD divider: during the hold, record the
net force on the divider F(t) (sum of impulses per output interval), compute
gamma(lambda) = beta * integral_0^inf <dF(0) dF(t)> dt with a plateau-value
analysis (plot the running integral; report the plateau and the time at which it
is reached). Do this for 3-4 eta at N=100 and at famA N=400. Do NOT use the
velocity autocorrelation integral (it is zero for a confined divider).

TASK 5 - table: eta, N, M, <v_w^2>/(kT/M), var(x)/(kT/k_eff), gamma_ringdown,
gamma_forceACF, with errors, plus the plots. Write
0000_PLAN_OVERALL/ALL_MARKDOWNS/260908_divider_fluctuations.md.

Rules: analysis of existing runs first; only Task 2's small batch if needed; no
core changes; no auto-commit.
```

---

## Message for GPT (status, so it can review the next CC reports)

```
Status update from the files on disk, 2026-09-08.

Pressure campaign (00_pressure_validation_20260907): 144 accepted, 18 discarded.
All discards are N >= 900 at eta 0.60-0.67, forced advances during equilibration
(fa 27-45, then thousands of overlap repairs) at chunk 0.8 (N=900/1600) and 0.4
(N=1600). Cause: EDMD_ADVANCE_MAX_EVENTS = 250000 per call is N-independent, so the
safe chunk ~ 1/N (0.8/0.4/0.2 at N = 400/900/1600); the 5-unit fresh-lattice
calibration does not see the equilibrated event rate. Fix chosen: runner-level
chunk = min(calibrated, 320/N) plus calibration on an equilibrated disposable
instance; rerun the five empty cells only.

Low density is done: Z_pair extrapolated in 1/sqrt(N) lands on Kolafa-Rottner
within 0.2% for eta 0.2-0.5 (the pair virial carries its own ~1/sqrt(N) surface
term, about half the wall route's); wall route within 0.4%; wall-pair gap scales
as 1/sqrt(N) at every eta; both routes agree on Z_inf.

Open: at eta = 0.69 the deviation from KR is -1.1 / -2.7 / -3.1% for N = 400/900/
1600, growing with N; both routes extrapolate to Z ~ 9.65 (betaP sigma^2 ~ 8.48),
5% below KR and below the Engel plateau 9.185; blocks stationary in Z, T, psi6
(global 0.1-0.2 falling as 1/sqrt(N), local 0.71). Equilibration ladder
400/1600/6400 at 0.69 is being run before interpreting. Proposed claim range: EOS
validated eta <= 0.60; 0.65-0.69 onset of a systematic departure; >= 0.70
exploratory. Also: accelerated core failed validation (181/450 invalid), excluded;
T_mean varies 0.91-1.06 between seeds in the pressure runner (harmless for Z;
being checked for the speed-of-sound runs).

Questions: (1) do you accept the runner-level chunk rule over a core budget change
for this campaign? (2) is 0.65 or 0.60 the right upper edge of the validation
claim given the 0.65/400 point already sits 1.5% below the surface-term trend?
(3) any objection to reporting BOTH routes extrapolated rather than calling the
pair virial the bulk route?
```

---

## How to run these

Prompt 1 today (the reruns are minutes at N=900, the 1600 cells and the 6400-unit ladder are hours — launch, then continue with Prompt 2 while they run). Prompt 3 needs nothing from 1 or 2 and can go in parallel. When the three `260908_*_final_analysis.md` notes exist, Paper 1 has its data; the Paper 2 groundwork prompt (energy-balance audit, thermal wall audit, ideal-gas Jarzynski) is section 4 of the paper-plan note and stays untouched until then.
