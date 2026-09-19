# The frequency estimator, its error bar, and the selection cuts

Written 2026-09-12. Every number is measured from the A2 trace cache
(1458 health-clean trajectories: the `famB_20260911` campaign plus the part-A top-up).
This document exists because the 0.03 % selection cut looked like it might be discarding
inconvenient data. The short answer is that it is not selecting on the *value* of the
result, but it **is** selecting on a quantity that turns out not to be a valid error
bar, and that is a real problem with a clean fix. The long answer follows.

---

## 1. What we actually measure

The divider is held fixed while the gas equilibrates, then released. We record its
position along x as a function of time. That is the entire raw observable: one
displacement time series per trajectory, typically ~7800 samples spanning ~30
oscillation periods, about 256 samples per period.

From each series we want one number, the oscillation frequency ν. The speed of sound is
then the slope of ν against x = K/(2π L_eff) across the ladder of divider masses.

## 2. The model we fit

We fit the post-transient series (first 20 % dropped, mean removed) with a damped
cosine:

    x(t) = A · exp(−γ t) · cos(2π ν t + φ) + c

Five free parameters: amplitude A, damping rate γ, frequency ν, phase φ, offset c. It
is seeded by linear prediction (Prony) and refined by non-linear least squares
(`scipy.optimize.curve_fit`).

## 3. Where σ_ν comes from — the covariance matrix, step by step

This is the part that needs to be understood before the cut can be judged.

**Step 1 — the cost function.** Least squares picks the parameter vector θ = (A, γ, ν,
φ, c) that minimises the sum of squared residuals

    S(θ) = Σᵢ [ xᵢ − model(tᵢ; θ) ]²

**Step 2 — curvature at the minimum.** Near the best fit θ̂, expand S to second order:

    S(θ) ≈ S(θ̂) + (θ − θ̂)ᵀ H (θ − θ̂),     H ≈ 2 JᵀJ

where J is the Jacobian, Jᵢⱼ = ∂model(tᵢ)/∂θⱼ. H measures how sharply the cost rises
when a parameter is moved. A sharply curved direction is a well-determined parameter.

**Step 3 — from curvature to an error bar.** If the residuals are independent draws
from a Gaussian of variance σ², standard theory gives

    Cov(θ̂) = σ² (JᵀJ)⁻¹

`curve_fit` does not know σ², so it estimates it from the fit itself,
σ̂² = S(θ̂)/(n − p) with n data points and p = 5 parameters, and returns

    pcov = σ̂² (JᵀJ)⁻¹

**Step 4 — the number we used.** σ_ν = √(pcov[2,2]), the diagonal entry for ν. The cut
`σ_ν/ν < 3×10⁻⁴` is a cut on this.

**The load-bearing assumption is in step 3: the residuals must be independent.** If they
are not, `(JᵀJ)⁻¹` is still the right shape but the scale σ̂² is meaningless, because
n correlated points carry far less information than n independent ones.

## 4. Why the assumption fails here

Measured at η = 0.10, N = 400, per divider mass:

| M | cycles recorded | samples/cycle | γ/ν | **residual RMS / fit amplitude** |
|---|---|---|---|---|
| 200 | 30.0 | 260 | 0.040 | **0.80** |
| 500 | 30.3 | 257 | 0.032 | **0.45** |
| 1000 | 30.3 | 255 | 0.029 | **0.50** |
| 2000 | 30.1 | 256 | 0.049 | **0.55** |

Sampling is not the limitation (256 points per cycle). Record length is not either
(γ/ν ≈ 0.04 means coherence is lost after ~25 cycles and we record 30).

The limitation is the last column: **between 45 % and 80 % of the divider's motion is
not described by the model**. And the reason is physical, not numerical.

A damped cosine describes a *free decay from a known initial displacement*. Our divider
is never displaced. It is released from whatever thermal pressure fluctuation happens to
exist at that instant, and thereafter every single gas collision kicks it again. It is a
**continuously driven stochastic oscillator**, not a decaying one. The residual is not
noise sitting on top of the signal — it is real, physical, driven motion for which the
model has no term.

And driven motion is strongly correlated in time. The residual at one instant looks very
much like the residual a few samples later. That is exactly the condition under which
step 3 breaks.

## 5. How badly it breaks — measured, not argued

The test is direct. Run the same cell many times with different seeds. ν is a physical
eigenfrequency, so the spread of the *fitted* ν across repeats is the true estimation
error. Compare it to what the fit claimed.

Over 72 cells with ≥ 8 repeats, using the median absolute deviation so that outliers do
not dominate:

| quantity | median | 10th pct | 90th pct |
|---|---|---|---|
| (robust scatter of ν) / (reported σ_ν) | **92.4** | 59.0 | 190.3 |
| (raw sd of ν) / (reported σ_ν) | 307.3 | — | — |
| pull (ν − mean)/σ_ν, sd | **246** (should be 1.0) | — | — |

**The reported error bar is too small by roughly a factor of 90 on the bulk of the
traces.** It is not a mildly optimistic error, it is wrong by two orders of magnitude.

This is quantitatively explained by the correlation. For correlated residuals with
correlation time τ, the effective number of independent points is n_eff ≈ T/τ rather
than n, and the error is underestimated by √(n/n_eff). Worked for M = 200, η = 0.10,
N = 400: the residual is driven motion whose correlation time is the mode coherence time
τ = 1/γ. With γ/ν = 0.040 and ν = 0.0038667, τ ≈ 6465 time units against a record of
T ≈ 7800 and Δt ≈ 1. So

    √(n / n_eff) = √(7800 / (7800/6465)) = √6465 ≈ **80**

against a measured median of **92**. The mechanism accounts for the effect.

## 6. Two further problems the same diagnosis exposes

**Non-convergence.** 331 of 1458 health-clean traces (23 %) return no frequency estimate
at all — the fit fails outright. At M = 50, η = 0.10, N = 400 it is **0 of 10**. That is
not an imprecise measurement, it is a missing one, and it must be reported as such.

**Gross outliers.** 133 of the 1127 converged fits (11.8 %) locked onto the wrong
spectral component entirely, with ν/ν_predicted ranging from 1.7×10⁻⁹ to 1.1×10³. These
are caught by the accuracy guard (§8). After that guard, a further ~9 % of traces still
sit more than 5 robust standard deviations from their cell median. Bad fits are common,
not exceptional.

## 7. So was the cut cheating?

Three separate questions, three separate answers.

**Was it selecting on the answer?** No. σ_ν/ν is computed from the fit's own covariance.
It never looks at the value of ν, never compares it to theory, and never compares it to
the other runs. Nothing in the cut can prefer a result that agrees with Kolafa–Rottner.

**Did it change the result?** No, and this was tested before the question was asked.
Re-running with no precision cut at all:

| η | 0.03 % cut | 0.10 % cut | no precision cut | Kolafa–Rottner |
|---|---|---|---|---|
| 0.10 | 1.74358 ± 0.00037 | 1.74408 ± 0.00163 | 1.74212 ± 0.00633 | 1.74414 |
| 0.30 | 2.85727 ± 0.00355 | 2.85282 ± 0.00653 | 2.86002 ± 0.00617 | 2.85150 |
| 0.50 | 5.48292 ± 0.00666 | 5.46425 ± 0.01752 | 5.47359 ± 0.02094 | 5.41622 |
| 0.60 | 8.41288 ± 0.01344 | 8.38413 ± 0.03202 | 8.57709 ± 0.14617 | 8.22261 |

The central values move by less than a percent. Dropping the lightest mass entirely,
which the cut nearly does, shifts c_s by at most **0.07σ** across all twenty cells,
because an inverse-variance-weighted fit gives an imprecise mass almost no weight
whether you keep it or not.

**Was it defensible?** No, and this is the honest part. The cut selected on σ_ν, and
§5 shows σ_ν is not a valid error estimate. Selecting on a quantity that is wrong by a
factor of 90 is selecting on something close to noise. It happened not to bias the
central value, but that was luck, not design, and it cannot be argued from first
principles in a paper. It also did real damage to the *claimed significance*: the
η = 0.60 excess reads 14.2σ with the cut and 2.4σ without. Quoting 14σ off a 45 %
selection is not defensible.

**Verdict: not fraud, but not a method we should publish.** It must be replaced.

## 8. What replaces it

Two filters, one kept and one discarded.

**Keep the accuracy guard.** Reject a fit whose ν differs from the driver's own
predicted frequency by more than a factor of 3. This is not a precision cut and not a
selection on agreement with theory — the predicted frequency comes from the run
parameters, not from the result. Over the traces that pass a precision cut,
ν/ν_predicted spans [0.977, 1.237], so a factor-3 window is 2.4× clear of all real data
and 86× clear of the aliases. Without it, 11.8 % of converged fits are nonsense.

**Discard the σ_ν cut.** Replace it with the error we can actually measure: the
**observed scatter across repeats**. Per divider mass, take the median of ν over
repeats, reject only points beyond 5 robust standard deviations of that median, and
weight by MAD/√n. This uses every trace, never consults σ_ν, and handles the ~9 %
outlier population by robustness rather than by exclusion.

Results, using 85–95 % of health-clean converged traces:

| η | c_s(N→∞) | Kolafa–Rottner | deviation | significance |
|---|---|---|---|---|
| 0.10 | 1.74549 ± 0.00033 | 1.74414 | +0.08 % | 4.1σ |
| 0.30 | 2.86066 ± 0.00470 | 2.85150 | +0.32 % | 1.9σ |
| 0.50 | 5.48279 ± 0.02789 | 5.41622 | +1.23 % | 2.4σ |
| 0.60 | 8.40238 ± 0.02403 | 8.22261 | +2.19 % | 7.5σ |
| 0.65 | 10.26047 ± 0.11847 | 10.43875 | −1.71 % | 1.5σ |

Every η now yields an extrapolation, including η = 0.65 which the 0.03 % cut had starved
to two system sizes. The physics conclusions are unchanged: agreement at the sub-percent
level for η ≤ 0.30, a real excess at 0.50 and 0.60.

## 9. Improving the measurement itself

The cut discussion is about analysis. Three ways to improve the underlying measurement,
in increasing order of cost.

**More repeats (free, already running).** The error on the per-mass mean falls as 1/√n
no matter how noisy a single trace is. This is why η = 0.10 went from ± 0.00576 to
± 0.00033 when the top-up tripled the statistics.

**A matched estimator (analysis only, prototyped, not yet validated).** For a thermally
driven damped oscillator the natural object is the power spectrum, averaged over a
cell's repeats before fitting:

    S(f) = A / [ (f² − f₀²)² + (Γ f)² ]

Prototype at η = 0.10, N = 400:

| M | time-domain ν | averaged-PSD ν | PSD σ_ν/ν |
|---|---|---|---|
| 50 | did not converge | 0.00510162 | 7.1×10⁻⁴ |
| 200 | 0.00585266 (aliased) | 0.00388779 | 1.6×10⁻⁴ |
| 500 | 0.00256684 | 0.00285387 | 1.2×10⁻⁴ |
| 1000 | 0.00213776 | 0.00214205 | 1.3×10⁻⁴ |
| 2000 | 0.00154868 | 0.00155419 | 1.8×10⁻⁴ |

It recovers M = 50 everywhere, uses every trace with no cut, and agrees with the
cut-filtered time-domain answer to 0.2–0.5 % on the heavy masses. It also benefits from
longer runs, which the damped cosine does not, because a driven mode never decays away —
it reaches a steady state whose spectrum keeps improving with record length.
**Caveat: the recovered M = 50 reads ~1.4 % high**, and until that is understood this is
not a drop-in replacement.

**A bigger signal (requires a code change, not done).** The root cause is that signal
and noise are both thermal, so the signal-to-noise per trajectory is order 1 *by
construction*. Displacing the divider by a known small amount before release, staying
within linear response, would give a large reproducible amplitude and collapse the
estimation error. That converts the experiment from an equilibrium-fluctuation
measurement into a driven one, which is a design change, not a tuning knob.

## 10. What goes in Paper 1

1. Primary analysis: **every health-clean trace with a converged fit**, robust
   per-mass median, weighted by measured scatter. No σ_ν cut.
2. The accuracy guard stays, with the alias statistics quoted.
3. Report that 23 % of traces yield no frequency estimate, and that it is concentrated
   at the lightest divider mass.
4. Report σ_ν as diagnostic only, with the explicit statement that it underestimates the
   true error by ~90× because the residuals are driven and correlated. This is a result
   about the method and belongs in the paper, not a footnote to hide.
5. Show the 0.03 % and 0.10 % variants as a robustness table, so a reader can see the
   conclusion does not depend on the choice.

---

## 11. Convergence-bias test (added 2026-09-12)

§6 reported that 23 % of health-clean traces yield no frequency from the time-domain
fit. That is only harmless if *which* traces fail is independent of the frequency they
would have given. If, for instance, realisations that happen to start at a low apparent
frequency fail more often, the surviving mean is biased upward and the 23 % is a
selection, not a missing measurement.

The binned FFT frequency exists for **every** trace, converged or not, so it is the
instrument for the test.

### Distribution test

Pooled over all A2 cells, using ν_binned/ν_predicted so that cell-to-cell differences
divide out (the cells that fail most have *no* converged traces, so normalising by a
converged within-cell median is blind to exactly the worst cases):

| population | n | median | IQR |
|---|---|---|---|
| converged | 1158 | 1.0332 | [0.9999, 1.0999] |
| non-converged | 331 | 0.9999 | [0.9331, 1.0665] |

Raw median shift −3.22 %, KS p = 2.8×10⁻²¹. **That number is misleading on its own.**
Restricting both populations to the physical band ν/ν_predicted ∈ [1/3, 3]:

| population | n | median |
|---|---|---|
| converged | 1135 | 1.0332 |
| non-converged | 249 | 1.0332 |

Median shift **−0.00 %**. The apparent −3.22 % comes entirely from the 24.8 % of
non-converged traces whose *binned* frequency is also garbage (against 2.0 % of
converged). For those, no estimator gives a usable number — the trace itself is
pathological, and the PSD estimator would not rescue them either.

For A1 the non-convergence rate is only 0.7 % (52 of 7186). Pooled against the
converged within-cell median: median ratio 1.0005, i.e. **+0.05 %**, sign-test
p = 0.68, Wilcoxon p = 0.54.

### Operational test

The decisive test is whether including the recoverable non-converged traces changes the
answer. Using the binned frequency for both arms with the identical alias guard, c_s was
computed per (η, N) from converged-only traces and from every trace:

- traces used: 1135 → 1384, i.e. **+249 recovered, 18 % more data**
- median relative shift in c_s: **+0.000 %**
- largest shift anywhere: **0.46 %** at η = 0.60, N = 400, which is **+0.21σ**
- every one of the 20 cells shifts by less than 0.21σ

**Verdict: non-convergence is a missing measurement, not a selection.** The time-domain
estimator stays primary for Paper 1. The PSD estimator remains the identified improvement
for the next paper, on the separate grounds of §9 (it recovers M = 50 and needs no
convergence), not because the present analysis is biased.

## 12. Final estimator and per-figure accounting (2026-09-12)

**Primary estimator.** Every health-clean trace whose fit converged; alias guard
ν/ν_predicted ∈ [1/3, 3]; per divider mass the robust median over repeats with
5-robust-σ rejection; weights from the measured scatter MAD/√n. **No σ_ν cut anywhere.**

**On the alias guard, stated plainly for the methods section.** The predicted frequency
is computed from the run's own parameters via an equation-of-state estimate of c_s, so
this window is centred on theory. It is defensible only because it is far wider than any
real effect: over the traces that pass, ν/ν_predicted spans **0.977 to 1.237**, so the
factor-3 window sits **2.4× clear of the real data** on both sides, while the rejected
traces are wrong by factors of 10⁻⁹ to 10³. It removes gross estimator failures and
cannot express a preference among physically plausible frequencies. Both numbers belong
in the paper.

### Trace accounting

| figure | total | health | non-converged | alias | robust-rejected | used |
|---|---|---|---|---|---|---|
| `260913_A2_cs_vs_N` | 1491 | 1 | 331 | 135 | 102 | **921** |
| `260913_cs_vs_eta` (A1) | 7200 | 14 | 52 | 685 | 444 | **6005** |

### A2 robustness table, c_s(N → ∞)

| η | robust | 0.03 % cut | 0.10 % cut | no cut (binned) | KR |
|---|---|---|---|---|---|
| 0.10 | 1.74540 ± 0.00029 | 1.74348 ± 0.00035 | 1.74408 ± 0.00163 | 1.73844 ± 0.00273 | 1.74414 |
| 0.30 | 2.85623 ± 0.00820 | 2.85817 ± 0.00350 | 2.86416 ± 0.00639 | 2.85606 ± 0.00336 | 2.85150 |
| 0.50 | 5.48279 ± 0.02789 | 5.48292 ± 0.00666 | 5.46425 ± 0.01752 | 5.47674 ± 0.03249 | 5.41622 |
| 0.60 | 8.40238 ± 0.02403 | 8.41288 ± 0.01344 | 8.38413 ± 0.03202 | 8.40717 ± 0.01359 | 8.22261 |
| 0.65 | 10.26047 ± 0.11847 | not fitted | 10.11961 ± 0.17747 | 9.89547 ± 0.03447 | 10.43875 |

The 0.03 % cut is the variant that cannot fit η = 0.65 at all. Against the robust
baseline it differs at η = 0.10 by 4.3σ — the estimator disagreement becoming visible
once the error bars are honest instead of ~90× too small.

The full A1 four-variant comparison over all 32 η is in
`260909_plots/260913_A1_robustness_table.csv`. The variants agree everywhere except
η = 0.710 (14.82 / 16.42 / 15.96 / 15.70), which sits inside the liquid–hexatic
coexistence window where the fluid branch is not defined in any case.

### Files

- `260913_cs_vs_eta`, `260913_cs_idealgas_zoom` — A1, 32 η plus the three low-density
  points, robust estimator
- `260913_A2_cs_vs_N_robust` plus `_cut0p03`, `_cut0p10`, `_nocut` — A2 and its
  robustness set
- `260913_A1_cs_vs_eta.csv`, `260913_A1_robustness_table.csv`,
  `260913_A2_cut_sensitivity.csv`, `260913_A2_cs_per_mass.csv`

One note on the A1 pipeline: the robust table is already T_i-corrected per trajectory,
so `plot_cs_meeting.py` is run in a mode that takes c_s verbatim. Running it in the
older correction mode would apply the per-η scalar correction a second time and inflate
every deviation by about 1 %.

---

## 13. Relation to Román et al. (2002): the FFT estimator (added 2026-09-13)

The experiment reproduces Román, González, White and Velasco, *Am. J. Phys.* **70**, 847
(2002). Their procedure, from §III of the paper:

- N₀ = 100 disks (50 per side), σ = m = k_BT = 1, box height A = 10, half-length
  L₀ = 7.5–35, a zero-width piston of mass M = 20–1000, held at the centre for about 10⁴
  collisions per particle and then released.
- Frequency: "a standard fast Fourier transform (FFT) method to obtain the power
  spectrum of the displacement ξ(t)", taking "the frequency corresponding to the maximum
  peak in the spectrum". One peak frequency per trajectory, **100 trajectories per case,
  averaged, with the standard deviation of the 100 values as the error**.
- c_s: least-squares slope of ⟨ν₁⟩ against K/2π(L₀ − 1), i.e. with the effective length
  L₀ − σ, the same L_eff = L₀ − 2r used here.
- Finite size (their Table II): η = 0.393, fixed aspect ratio A = L₀ and **fixed
  K = 1.07687, i.e. M = Nm scaled with N**; N = 64, 256, 1024, 4096 per side give
  c_s = 3.81, 3.78, 3.76, 3.75, extrapolated 3.74.

What this means for this work:

1. **Román's error bar is the empirical scatter across trajectories, not a fit
   covariance.** The robust estimator of §12 returns to exactly that principle; the σ_ν
   cut of §7 was the departure from it.
2. **The binned "no cut" column of §12 is Román's estimator** (per-trajectory FFT peak,
   averaged over repeats). On A2 it agrees with the robust estimator within errors at
   η = 0.30, 0.50 and 0.60, and differs at η = 0.10 (1.73844 ± 0.00273 against
   1.74540 ± 0.00029, −0.40 %, 2.5σ) and η = 0.65 (9.89547 ± 0.03447 against
   10.26047 ± 0.11847, −3.6 %, 3.0σ).
3. **The FFT peak is quantised to one bin**, Δf = 1/T. A record holding n periods resolves
   ν only to 1/n of itself; the ~30 periods recorded here give the ~3.3 % measured in §6.
   Román reduced this by averaging 100 trajectories; the mode's own linewidth, γ/ν ≈ 0.04,
   sets a comparable floor for any single trajectory.
4. **His finite-size study held K, and therefore α = M/(2N_side m), fixed while N grew.**
   The A2 ladder did not, which is the α confound; the α = 2 cells restore his design at
   N = 1600 and 2500.
5. **The improved version of his method is the averaged-spectrum line-shape fit of §9**:
   the same FFT power spectrum, averaged over a cell's repeats before the peak is located,
   and fitted with the driven damped-oscillator line shape instead of reading off the
   highest bin. It keeps the lineage of the reference method and removes its bin
   quantisation. It is not yet validated: the M = 50 point read 1.4 % high in the
   prototype.

**Correction.** `plot_cs_meeting.py` carried Román's L₀ = 25 value (η = 0.157) as
2.10 ± 0.02. Table I of the paper gives **2.01 ± 0.02**. Corrected 2026-09-13 and
`260913_cs_vs_eta` regenerated.
