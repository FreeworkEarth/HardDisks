# 260909 — In plain words: what the review said, the two plots, why the ideal gas matters, what "nostalgia" is, and where we are on the timeline (COWORK)

Written 2026-09-09. Companion to `260909_prompt1_results_review_COWORK.md` (the technical version) and a reply to GPT's `260909_speedsound_structure_handoff_CODEX.md`. Plots in `260909_plots/`.

---

## 1. The review, without the jargon

You ran a big pressure experiment: put N disks in a box at a given density, let them bounce, measure how hard they push (that is the pressure, written as Z = P/(ρ k_B T), "how many times the ideal-gas pressure"). Then compare with the best known formula for hard disks (Kolafa–Rottner). CC did that, cleanly, for 162 runs. Good.

I checked CC's numbers against the files and then argued about four sentences in CC's write-up. Here they are, in order of how much they matter.

**"16× longer equilibration confirms it."** CC ran the η = 0.69 system with waiting times of 400, 1600 and 6400 before measuring. But all three used the same random seed and the same settings, and the code is deterministic, so they are literally the same movie watched at three later times — I checked: every number in the 400-run appears identically in the 6400-run at the same time stamp. So the right sentence is "the pressure at η = 0.69 does not drift over 7000 time units". Still useful. Just not "three independent equilibrations".

**"Above 0.65 the deviation grows with N."** To get the pressure of an infinitely large system from finite boxes, you plot Z against 1/√N and extend the line to zero. That trick works when the walls only add a thin "skin" effect. At η = 0.60 and 0.65 it works beautifully (the points sit on a line). At η = 0.67 they do not: the N = 1600 point sits *above* the N = 900 point, the line is a bad fit (χ² = 7.6 for one degree of freedom, about 1-in-170 odds), and the extrapolated number changes by 1 % depending on which points you include. At 0.69 the slope of that line is four times larger than everywhere else. Translation: from 0.67 up, the walls are not a skin effect any more — the correlated patches in the fluid are as big as the box — so the extrapolation is not a measurement of the bulk. Report the finite-box values there and say so. "Validated to 0.65" is the honest line, and 0.65 itself is 0.42 % below the formula, which is real but tiny, and we will learn what it means from one more box size (N = 2500).

**"The threshold is exactly 320/N."** The simulator has a safety fuse: if one advance step processes more than 250 000 calendar entries, it forces the disks forward and repairs overlaps — those runs are (rightly) thrown away. CC's rule "step length ≤ 320/N" keeps the fuse from blowing. Fine. But it is a safe cap with a factor 1.1–4 to spare depending on density, not the edge; at η = 0.60 the fuse only blew at 1280/N.

**"250 000 events."** A referee will read "events" as collisions and compute that 900 disks at η = 0.69 make about 3 700 collisions per step, not 250 000. The counter actually counts every entry popped from the event calendar, including stale ones that get thrown away, about 30–45 per real collision. And all 345 fuse trips happened in the last 1–8 % of a step, which means it is a plain throughput limit, not some avalanche. This is a one-paragraph fix in the methods; it just has to be said correctly.

Net: the campaign is done and the low-density result is paper-grade. The high-density story got *more* careful, not worse.

---

## 2. The plots

**`260909_Z_idealgas_zoom.png / .pdf`** — the pressure at the ideal-gas end (η ≤ 0.1), three panels: (a) Z against η with the three exact anchors (Z = 1 ideal gas, Z = 1 + 2η second virial, Kolafa–Rottner), our points for N = 400/900/1600 and the N → ∞ stars sitting on the line; (b) the same as percent deviation, everything inside ±0.2 %; (c) the slope test (Z − 1)/η, which must approach exactly 2 as η → 0 (that "2" is the second virial coefficient of hard disks in these units, a number with no adjustable parameter). It does: 1.98 ± 0.06 at η = 0.005. Every number on that figure is parsed from `analysis/tables_ABC.md`; the provenance line is on the figure.

**`260909_cs_vs_eta.png / .pdf`** — speed of sound against η. The curve is what the equation of state predicts through Román's relation, the dashed line is the ideal gas (c_s = √2 in units of √(kT/m) — that is where γ = 2 for a 2D monatomic gas comes in), the inset zooms into η ≤ 0.1 with the exact first-order law c_s = √2 (1 + 2η). The points are the six route-A/route-B densities I have on record from CC's fits. **This one is provisional**: the full 33-density route-A table lives in `campaign_r25_psi6_20260823/analysis/final_plots/speed_of_sound_summary.csv`, which sits one folder too deep for my bridge to your Mac. Copy it up (Finder drag, or `cp` it into `0000_PLAN_OVERALL/ALL_MARKDOWNS/260909_plots/`) and I regenerate the figure with all 33 points and error bars; the script `plot_cs_vs_eta.py` already takes that file as its second argument. **And there is a catch that GPT found** (section 5): at the three lowest densities (η = 0.0196, 0.026, 0.039) 72–80 % of the route-A trajectories carry a `wall_overdue` health flag. Under our own strict rule those points are not usable until someone explains the flag. So the ideal-gas end of the *sound-speed* data is exactly the part that is currently under a cloud, while the ideal-gas end of the *pressure* data is spotless. That has to be resolved before the c_s figure goes anywhere near Susanne.

---

## 3. Why your professor cares about the ideal-gas limit

Three reasons, all the same reason really: it is the only place where the answer is known exactly, with nothing to fit.

1. **Every formula becomes exact.** PV = NkT, Z = 1, c_s = √(γ kT/m) with γ = 2 in 2D, Maxwell velocities, equipartition. If the simulation is off there, nothing at higher density can be interpreted. And the *first correction* is exact too: Z = 1 + 2η + …, so the slope of Z at small η must be exactly 2 (panel c). That checks the collision counting itself — the virial machinery, the wall-impulse machinery, the time averaging, the temperature definition — with no physics of dense fluids involved.
2. **The walls stop mattering.** The finite-size term a/√N in the pressure fits is 0.004 at η = 0.005 and 2.6 at η = 0.60. At low density a box of 400 disks already *is* the infinite system, so you can validate the instrument without any extrapolation argument.
3. **It is where Paper 2 starts.** The driven-piston work distributions (Lua–Grosberg), the Jarzynski check, the Szilard engine, and the "thermodynamics of prediction" box are all analytically solvable *for an ideal gas*. Your non-equilibrium experiments will be compared to theory first in exactly this regime. A simulator that is exact at η → 0 in equilibrium is the precondition for trusting it there out of equilibrium.

That is also why the `wall_overdue` flags at η = 0.02–0.04 in the sound-speed runs matter more than their count suggests: they sit precisely at the anchor.

---

## 4. "Nostalgia", again

From Still, Sivak, Bell and Crooks, *Thermodynamics of Prediction* (PRL 2012). A system with state s is driven by a signal x. Its state carries information about the *current* signal, I(s_t ; x_t) — call that memory. Part of that information is also useful for the *next* value of the signal, I(s_t ; x_{t+1}) — call that predictive power. The difference,

nostalgia = I(s_t ; x_t) − I(s_t ; x_{t+1}),

is the information the system holds about the past that is useless for the future. The theorem: the work you dissipate per step is at least k_B T times the nostalgia. Useless memory costs energy. In our setting the "system" is the divider (or a small predictor box) and the "signal" is the piston or the gas hitting it; a divider that is heavy and slow remembers what the piston did a while ago (high nostalgia) and dissipates; a divider tuned to track the piston's *next* move dissipates less. That is the physics link between Paper 2's energy transfer and Susanne's information-theoretic program — and it is only measurable once the energy ledger (work in = Δ gas + Δ divider + Δ spring + heat out) closes to round-off, which is Level 0 of the Paper 2 ladder.

---

## 5. GPT's package and its handoff: what I agree with, what I do not

GPT (Codex) built `analysis_paper1_20260908/`: 18 938 sound-speed trajectories screened for structural change (ψ₆ at hold vs end), 4 plots, exclusions, 19 tests, nothing deleted. I read the handoff, the README and the tables. Verdict:

- **The 505 health flags are the real finding.** 504 are `wall_overdue = 1 or 2` at route A's three lowest densities (162, 180, 162 out of 225 each), plus the one clamp repair at N = 1000 / η = 0.70. Earlier status notes (mine included) said "one bad trajectory". That was wrong; the ledgers marked them valid and nobody looked at the counters. What `wall_overdue` means in a dilute gas with L₀ = 100–200 (a divider collision found later than scheduled? a bookkeeping event at release?) decides whether the ideal-gas end of the c_s data is usable. This is now the first task for CC — GPT's item 1 — and I would put it ahead of everything else in the sound-speed work.
- **No measured temperature in the c_s records.** Also true and also matters: c_s scales as √T, and the pressure runner showed seed-to-seed T scatter of ±5 % at N = 400 (velocities drawn, never rescaled). If `--kbt1` rescales at release, fine, but the record has to say so. GPT's item 2. Agree.
- **Structural change during measurement** at route B 0.72 (ψ₆ 0.26 → 0.69) and fixed-aspect N = 1600 at 0.70/0.71: agrees with what I wrote on 260907; those cells are excluded from the fluid-branch claim anyway.
- **What GPT got right in its critique of my plan:** Jarzynski does not need a thermal wall — you can sample canonical initial conditions directly and drive the piston; the finite isolated divider has microcanonical, not canonical, statistics (so "equipartition of the divider" is an approximation at N = 100); the collisionless piston regime (Lua–Grosberg) and the quasistatic isotropic-gas formula are different limits and must not be mixed. All three are bounded corrections to `260907_paper_plan...`; I will fold them in.
- **Where I disagree:** "EOS validated to η ≤ 0.60 is stronger than I would write" — no, with 3 sizes, χ² ≈ 1, ±0.2 % agreement, two routes, that is exactly what the data support (and after the χ² check it extends to 0.65). The sensitivity GPT asks for *is* in the fit table now (intercept from 900/1600 only). And "two-week timeline unreliable" — sure, it was a checklist, not a promise.

---

## 6. Where we are on the walkthrough's checklist

| Step (260907 ELI5 walkthrough) | Status 2026-09-09 |
|---|---|
| 1–2 Pressure campaign, dense reruns, equilibration ladder | **Done.** 162 accepted, claim range fixed (≤ 0.65 validated; 0.67–0.69 finite-box only; ≥ 0.70 exploratory). One optional run left: N = 2500 at 0.65/0.67/0.69 to settle the 0.65 point. |
| 3 Sound-speed stationarity screen | **Done by GPT** (structural screen). Not done: the 504 `wall_overdue` flags, temperature provenance, and the refits on eligible trajectories. |
| 4 Corrected c_s fits and final figures | **Open.** Blocked on step 3's two items. The figure template exists (`plot_cs_vs_eta.py`). |
| 5 Initialization / hold-length comparison | Open; only if the refits need it. |
| 6 Divider equilibrium fluctuations (Level 1 of Paper 2) | Open; needs proper equilibrium runs (not ring-down). |
| 7+ Paper 2 Level 0: energy ledger on the existing driven runs | **Open and unblocked** — can start now, in parallel. |

Paper 1 is roughly two-thirds there: the pressure half is finished; the sound-speed half needs the flag question answered, then one refit pass and the figures. Paper 2 has not started but its first step needs no new simulation. So no, you do not "need a lot" to begin Paper 2 — you need the energy ledger closed on runs you already have; everything else on the Level 0–6 ladder comes after that, and Levels 3, 5, 6 are only prerequisites if Paper 2 claims them.

---

## 7. Prompt for Claude Code (paste as is; combines GPT's handoff with the pressure follow-ups)

```
Read, in this order: 0000_PLAN_OVERALL/ALL_MARKDOWNS/260909_speedsound_structure_handoff_CODEX.md,
260909_prompt1_results_review_COWORK.md, 260909_plain_language_status_plots_and_next_steps_COWORK.md.
Rules unchanged: no core physics changes, no edits to accepted trajectories, no commits, no new
campaigns without my go, quote numbers from files.

TASK A (first, sound speed): the 504 wall_overdue flags.
1. Pick 5 flagged trajectories at eta=0.019635 and 5 unflagged ones from
   campaign_r25_psi6_20260823/eta_0p019635/. From run.log and the wall_x trace, state for each
   flagged one: at what time the overdue event happened (release? during hold? during
   oscillation?), which wall, and by how much it was overdue. Then quote the code path in
   edmd_core/edmd.c that increments wall_overdue and say in one paragraph what the counter means
   physically. Decide: (a) bookkeeping at release -> flag is harmless, document and keep the
   trajectories; (b) a missed/late collision during measurement -> keep the strict rule, the
   three lowest densities are lost and need a rerun with a documented fix. Do not pick (a)
   without the time stamps.
2. Temperature provenance: find where --kbt1 acts (rescaling at release? at start?) and whether
   any durable record has the measured T per trajectory. Report the code line.

TASK B (pressure follow-ups from the review, small):
3. Apply the write-up corrections listed in 260909_prompt1_results_review_COWORK.md section 7
   TASK A (claim-range text, chi2 column, D1 same-trajectory statement, 320/N wording,
   methods paragraph, calibration stderr logging). Regenerate 260908_pressure_final_analysis.md.
4. Copy campaign_r25_psi6_20260823/analysis/final_plots/speed_of_sound_summary.csv and
   combined_speed_of_sound_summary.csv into 0000_PLAN_OVERALL/ALL_MARKDOWNS/260909_plots/
   (copies, so Cowork can reach them) and print their header lines.

TASK C (only after A is answered): the refits GPT's handoff items 3-4 describe, on eligible
trajectories only, with the named EOS and the Roman relation with exact cot K. Then run
0000_PLAN_OVERALL/ALL_MARKDOWNS/260909_plots/plot_cs_vs_eta.py <out> <summary csv> and show me.

Report commands and outputs verbatim; git status --short at the end.
```

Paper 2 Level 0 (energy ledger) is a separate prompt; I will write it once you say the c_s side is moving, so the two do not fight for the machine.
