# Speed-of-sound estimator for Paper 1 (method note)

Written 2026-09-13, revised 2026-09-14. Supersedes the damped-cosine estimator, the σ_ν cut, the alias guard and the
MAD rejection introduced 2026-09-10 to 2026-09-12, and the unrestricted largest-bin rule of the 2026-09-13 version.
None of those is used anywhere in Paper 1. All numbers below are read from the CSV files named with them.

## 1. Per trajectory

One seed gives one divider position record $x(t_j)$ at fixed sampling interval after release, from the release
on (no transient drop). With $T$ the record length and $\nu_{\mathrm{pred}}$ the predicted mode frequency carried in
the trace header, $N_{\mathrm{cyc}} = T\,\nu_{\mathrm{pred}}$ is the record length in predicted periods.

$$P(f) = \Big|\sum_j \big(x(t_j) - \bar x\big)\, e^{-2\pi i f t_j}\Big|^2,
\qquad \nu_r = \arg\max_{f \,\ge\, \nu_{\mathrm{pred}}/X} P(f), \qquad X = 2.5 .$$

In bins, $k \ge k_{\min} = \mathrm{round}(N_{\mathrm{cyc}}/X)$: 10 at 25 periods, 15 at 37.5 (A2), 80 at 200 (A1 v2).
The frequency of the run is the centre of the largest bin at or above that edge. Bin width $\Delta f = 1/T$, so
$\Delta f/\nu = 1/N_{\mathrm{cyc}}$.

**What the edge is.** A one-sided lower edge at $\nu_{\mathrm{pred}}/2.5$, with no upper edge. It is placed with the
predicted frequency, so the estimator is not theory-free. It sits a factor 2.5 below the predicted frequency, while
the effects the paper tests are 1–3 %, and every edge from $\nu_{\mathrm{pred}}/2.08$ to $\nu_{\mathrm{pred}}/6.25$
gives the identical $c_s$ at 200 periods (table 1.1). Why a lower edge is needed at all is section 7: the divider
has a second, slow mode whose power sits at low frequency.

### 1.1 Sensitivity, A1 v2 (35 η, 9 masses, 25 seeds, 7873 health-clean trajectories)

$c_s$ from each variant against the largest bin inside $[\nu_{\mathrm{pred}}/3,\ 3\nu_{\mathrm{pred}}]$ ($w = 3$).
Entries are max over all 35 η of $|\Delta c_s|/c_s$ in %. $X = N_{\mathrm{cyc}}/k_{\min}$.
Sources: `260909_plots/260914_A1v2_kmin_verification.csv`, `260914_A1v2_kmin_sensitivity.csv`,
`260914_A1v2_kmin_sensitivity_extended.csv`, `260914_A1v2_kmin_sensitivity_X2p5.csv`,
`260914_A1v2_velocity_spectrum_test.csv`, `260914_A1v2_pertraj_fit_and_windowed_peaks.csv`.

| variant | 25 periods: X | max \|Δ\| % | 200 periods: X | max \|Δ\| % |
|---|---|---|---|---|
| literal largest bin, k ≥ 1 (Román 2002 rule) | 25 | 33.339 | 200 | 93.024 |
| k_min = 2 | 12.5 | 7.467 | 100 | 82.977 |
| k_min = 3 | 8.33 | 2.343 | 66.7 | 70.518 |
| k_min = 4 | 6.25 | 1.817 | 50 | 57.788 |
| k_min = 8 | 3.13 | 0.627 | **25 (failed proposal)** | 21.622 |
| k_min = 10 | 2.5 | 0.000 | — | — |
| k_min = 12 | 2.08 | 0.000 | — | — |
| k_min = 16 | 1.56 | 0.000 | 12.5 | 1.281 |
| k_min = 24 | — | — | 8.33 | 1.281 |
| k_min = 32 | — | — | 6.25 | 0.000 |
| k_min = 48 | — | — | 4.17 | 0.000 |
| k_min = 64 | — | — | 3.13 | 0.000 |
| k_min = 80 | — | — | **2.5 (adopted)** | 0.000 |
| k_min = 96 | — | — | 2.08 | 0.000 |
| velocity spectrum $(2\pi f)^2 P(f)$, k ≥ 2 | — | 3.334 | — | 2.226 |

The k_min = 16 and 24 failures at 200 periods are one density (η = 0.609999, −1.281 %). The k_min = 8 failure at 25
periods is one density (η = 0.630002, −0.627 %). On A2 (N = 100 … 2500, 37.5 periods, 23 (η, N) points) the adopted
edge equals $w = 3$ at every point (Δ = +0.000 %, `260914_A2_cs_vs_N` run output).

**Velocity spectrum, tested and rejected as primary.** Weighting by $(2\pi f)^2$ removes the slow mode without any edge
(its largest bin never lies above $3\nu_{\mathrm{pred}}$; worst high-frequency power ratio 0.073), but it moves $c_s$
up by +0.02 to +3.3 %, one-signed and largest at η ≥ 0.52. That is a systematic shift from the spectral shape, not
noise, so it fails the 0.1 % rule and is not used.

## 2. Per cell $(\eta, N, M)$, $R$ seeds

$$\bar\nu = \frac1R\sum_{r=1}^R \nu_r,\qquad
s_\nu = \sqrt{\frac{1}{R-1}\sum_r(\nu_r-\bar\nu)^2},\qquad
\sigma_{\bar\nu} = s_\nu/\sqrt R .$$

Averaging over seeds helps because the thermal driving moves the peak between neighbouring bins from run to run; it
does not beat the bin width by itself.

## 3. Across the divider-mass ladder

$$x_M = \frac{K(M/N)}{2\pi L_{\mathrm{eff}}},\qquad L_{\mathrm{eff}} = L_0 - 2r,\qquad
\cot K = \frac{M}{2N_s m}\,K,$$

with $N_s$ the disks per side. $c_s$ = slope of $\bar\nu$ against $x_M$ through the origin; its error = 1σ scatter
of the per-mass values $\bar\nu_M/x_M$. Masses 50, 100, 200, 300, 500, 750, 1000, 1500, 2000 (A1); 50, 200, 500,
1000, 2000 (A2), with 3200 at $N = 1600$ and 5000 at $N = 2500$ so that $\alpha = M/N$ reaches Román's value 2.
A2 finite size: $c_s = c_\infty + b/\sqrt N$ per η, weighted by the mass scatter.

### 3.1 Which estimator the mass ladder itself prefers (2026-09-17)

The mass ladder is an internal consistency test that never mentions an equation of state: one $c_s$
must put **every** mass on the same line through the origin, so the estimator that leaves less scatter
about that line has more evidence behind it. Metric: rms relative residual about the through-origin fit,
per $(\eta, N)$ cell.

| record | $N_\mathrm{cyc}$ | bin-rounding floor $100/(N_\mathrm{cyc}\sqrt{12})$ [%] | median rms residual, largest bin [%] | median rms residual, per-trajectory fit [%] |
|---|---|---|---|---|
| A1 v2, 24 densities η ≤ 0.69, 9 masses × 25 seeds | 200 | 0.144 | 0.270 | 0.279 |
| A2 dilute, 8 cells, 5 masses × 10 seeds | 50 | 0.577 | 0.499 | 0.453 |

**The size of the scatter does not discriminate.** At 200 periods the two agree to 0.01 % (the fit is
closer to the line at 10 of the 24 densities, the largest bin at 14); at 50 periods the fit is better by
0.05 %. Both sit within a factor two of the bin-rounding floor.

**The structure of the residual does discriminate, and it selects the largest bin.** A biased estimator
does not scatter about the line, it leans. Mean residual per mass, same 24 densities:

| $\alpha = M/N$ | largest bin [%] | per-trajectory fit [%] |
|---|---|---|
| 0.5 | −0.051 ± 0.092 | **+0.349 ± 0.058** |
| 1 | +0.014 ± 0.102 | +0.051 ± 0.023 |
| 2 | +0.066 ± 0.088 | −0.206 ± 0.033 |
| 3 | +0.120 ± 0.057 | −0.289 ± 0.045 |
| 5 | +0.027 ± 0.071 | −0.434 ± 0.083 |
| 7.5 | +0.015 ± 0.080 | −0.488 ± 0.086 |
| 10 | −0.043 ± 0.106 | −0.537 ± 0.111 |
| 15 | −0.167 ± 0.126 | −0.626 ± 0.127 |
| 20 | −0.150 ± 0.116 | **−0.641 ± 0.127** |

Residual against $\ln\alpha$:

| estimator | slope [% per e-fold in $\alpha$] | significance |
|---|---|---|
| largest bin (X = 2.5) | −0.030 ± 0.029 | 1.0 σ — **flat** |
| per-trajectory fit | −0.286 ± 0.020 | 14.4 σ — **leans** |

The fit walks monotonically from +0.35 % at the lightest divider to −0.64 % at the heaviest, a one-percent
swing across the ladder; the largest bin is flat within its errors. The fit is therefore the *more precise*
estimator (leave-one-mass-out spread 0.174 % against 0.284 %) and the *more biased* one — lower variance
bought with a mass-dependent shift. No equation of state entered this judgement.

Conclusion: **the largest-bin rule (X = 2.5, § 1) is primary on its own merit**, not only for
comparability with Román 2002. The per-trajectory fit remains a cross-check, and its α-dependent offset is
reported as an estimator systematic — never split by density or mass regime. Script
`validation/estimator_massladder_20260917.py` (`main`, and `extra` for the two diagnostics above).

## 4. Health contract (the only discard)

Any nonzero EDMD health counter discards that trajectory, never the cell. The binary prints the
`[EDMD-HEALTH]` line only when a counter is nonzero, so a missing line means all counters are zero.

| Counter | Meaning | Action |
|---|---|---|
| `forced_advance` | scheduler could not proceed; time forced forward | discard trajectory |
| `overlap_repair` | two disks overlapped and were separated by code | discard trajectory |
| `wall_overdue` | a wall collision was processed after its due time | discard trajectory |
| `clamp_repair` | wall position forced back into the box | discard trajectory |
| trace missing / truncated | no complete $x(t)$ | discard trajectory |

Everything else is kept. Report per figure: trajectories run, discarded, used.

## 5. Cross-checks (reported, never headline)

- Resonance fit to each trajectory's spectrum in $[\nu_{\mathrm{pred}}/2,\ 2\nu_{\mathrm{pred}}]$, mean over seeds.
- Resonance fit to the seed-averaged spectrum.

Both use a band placed around the resonance. Against $w = 3$ the per-trajectory fit differs by at most 1.408 %
(25 periods) and 1.433 % (200 periods) for η ≤ 0.69 (`260914_A1v2_pertraj_fit_and_windowed_peaks.csv`).
Optional: parabolic interpolation of $\log P$ over the three bins around the maximum.

**Damping shift of the position peak (2026-09-15, `260909_plots/260915_A1v2_damping_test.csv`).** For a linearly
damped oscillator under flat forcing the position spectrum peaks below the eigenfrequency by
$\Gamma^2/4f_0^2$, while the velocity spectrum peaks at $f_0$. The per-trajectory fits give
$\Gamma/f_0 = 0.023$ (median, η < 0.5) and $0.063$ (η ≥ 0.5), so that shift is at most 0.34 % and typically
0.01–0.11 % — an order of magnitude below the +0.87 % median gap between the velocity and position peaks at
η ≥ 0.5. The primary estimator is therefore the position-spectrum maximum, for comparability with Román 2002,
and the bias against the eigenfrequency from damping is bounded by 0.34 % over η ≤ 0.69. What remains of the
velocity–position gap is a spectral-shape effect, not the damping shift (§ 5 of the report addendum of
2026-09-15).

## 6. Campaign status (2026-09-14)

- **A1 v2** complete: 7875 trajectories, 2 discarded (`wall_clamp_repairs=1`, at L0 = 600 and 400), fixed seed pad,
  `--seed-drift-order=drift-first`, 200 target oscillations, 25 seeds, nine masses, 35 densities. $T_i = 1$
  verified: all 7875 audits read `2b after per-segment equalize ... KE_left=50 KE_right=50 ... kT_mean=1`.
  Final figures `260914_cs_vs_eta`, `260914_cs_idealgas_zoom`; table `260914_A1v2_final_cs_vs_eta.csv`.
- **A2** re-analysed with section 1: 1555 trajectories used, 1 discarded. N = 2500 cells of the top-up hold 4–7 of 10
  seeds (part B not resumed). Figure `260914_A2_cs_vs_N`; tables `260914_A2_cs_per_mass.csv`,
  `260914_A2_cs_vs_N_extrapolation.csv`.
- **α = 2 cells** (M = 3200 at N = 1600, M = 5000 at N = 2500, η = 0.10 and 0.30, 10 seeds, A2 settings) running
  since 2026-09-14 00:22 HST via `validation/resume_runs.py --plan alpha`. A4 paused.

## 7. Slow wander of the divider

**What it is.** Besides the acoustic resonance, the divider position has a slow component: a 5-period running mean
of $x$ carries 9–57 % of the variance of $x$ (TEST A, 1000 periods, medians over 25 seeds). It is not a start-up
transient: every run starts with exactly $KE_L = KE_R = 50$, its variance is the same in the second half of the
record as in the first (ratio 0.72–1.25), a 100-fold longer hold does not change it (TEST B), and its correlation
time is 53–102 acoustic periods. It is not an oscillation: it is a relaxation.

**Physics: the adiabatic piston.** The acoustic mode is fast, so each compartment responds adiabatically. Heat
crosses the divider only through divider collisions, which is slow. Thermal fluctuations of the two compartment
temperatures shift the pressure balance, the divider follows, and it relaxes back as heat flows. The fast mode sees
the adiabatic stiffness (the one in $c_s^2 = (k_BT/m)[Z + \eta Z' + Z^2]$); the full, slow equilibrium sees the
isothermal stiffness.

**7.1 Amplitude against equipartition** (`260914_wander_equipartition_A1v2.csv`). The free energy of the two gases
gives $k_{\mathrm{eff}} = 2N_sk_BT[Z + \eta Z']/L^2$ and

$$\langle x^2\rangle_{\mathrm{iso}} = \frac{L^2}{2N_s[Z + \eta Z']},\qquad
\langle x^2\rangle_{\mathrm{adi}} = \frac{L^2}{2N_s[Z + \eta Z' + Z^2]},\qquad N_s = 50,\ Z\ \text{from KR}.$$

rms $x$ about 0 (the symmetric equilibrium), last 800 of 1000 periods, 25 seeds, $L = L_{\mathrm{eff}}$:

| L0 | η | M | rms x | measured / isothermal | measured / adiabatic |
|---|---|---|---|---|---|
| 20.0 | 0.196350 | 20 | 1.2112 | 0.973 | 1.390 |
| 20.0 | 0.196350 | 100 | 1.1976 | 0.962 | 1.374 |
| 20.0 | 0.196350 | 1000 | 1.0009 | 0.804 | 1.148 |
| 7.5 | 0.523599 | 20 | 0.1698 | 1.000 | 1.549 |
| 7.5 | 0.523599 | 100 | 0.1772 | 1.044 | 1.617 |
| 7.5 | 0.523599 | 1000 | 0.1544 | 0.909 | 1.409 |

For light dividers on long records the isothermal form holds to 4 % at both densities; the adiabatic form overshoots
by 37–62 %. Heavy dividers and 200-period records fall below it (A1 v2 M = 1000 median 0.777): the slow mode starts
from $T_L = T_R$ and is slower for heavy dividers, so it has not reached full amplitude. At η ≥ 0.67 the ratios fall to
0.4–0.7 where KR's $Z + \eta Z'$ turns down (41.70 at η = 0.67, 29.37 at η = 0.69), outside KR's range.

**7.2 Direct test with the compartment kinetic energies** (5 seeds, L0 = 20, 200 periods, trace columns
$KE_L, KE_R$; `adiabatic_piston_check_20260913/step3_results.json`). $T_s = KE_s/N_s$; 5-period running means.

| M | mean r(x_slow, ΔT raw) | mean r(x_slow, ΔT slow) | mean sd(x_slow) | mean (L0/2)·sd(ΔT slow)/T · Z/(Z+ηZ′) | ratio |
|---|---|---|---|---|---|
| 100 | +0.531 | +0.978 | 0.5105 | 0.5472 | 0.933 |
| 1000 | +0.302 | +0.910 | 0.1868 | 0.1984 | 0.941 |

The sign is positive: a hotter left compartment pushes the divider towards +x. At M = 100 the lag of the maximum
correlation is ±0.25 periods. The pressure balance $x = (L_0/2)(\Delta T/T)\,Z/(Z + \eta Z')$ reproduces the amplitude
to 6–7 %; the ideal-gas form without $Z/(Z + \eta Z') = 0.6675$ does not, nor does raw $\Delta T$, which carries fast
fluctuations the divider cannot follow. The instantaneous scale is $\mathrm{sd}(T_L - T_R)/T \approx \sqrt{2/N_s} = 0.20$
(measured 0.11–0.17 in these records).

**7.3 Consequence for the estimator.** The slow mode's power is spread over low frequencies with a long tail, not
concentrated in a line: in the seed-averaged spectrum at 200 periods the largest power above bin 8 reaches 12–30 % of
the resonance peak at η ≤ 0.26 and 85–280 % at η ≥ 0.52 for the lightest divider (`260914_A1v2_kmin_sensitivity.csv`).
The longer the record, the more of that power is resolved into separate bins that can outgrow the resonance bin:
the unrestricted largest-bin rule is off by at most 33 % against $w = 3$ at 25 periods and by up to 93 % at 200 periods
(table 1.1). Short records therefore hide much of the problem, but in our data they do not remove it; a lower edge is
needed at every record length we use, and the edge adopted in section 1 is verified at both.

**Why Román 2002 apparently did not need this floor (INFERENCE, not stated in his paper).** His
Ref. 7 is *Numerical Recipes in C*, 2nd ed., chapter 13, whose spectral-estimate routine is the
**Welch segmented periodogram**: the record is cut into overlapping segments, each is windowed and
transformed, and the results are averaged. Segmenting has a specific consequence here — the lowest
frequency a segment can resolve is 1/T_seg, not 1/T, so power below that is not resolved into
separate bins but folded into the DC/first bin, which is exactly the region his rule discards
("the nonzero value of the frequency corresponding to the maximum peak"). A segmented estimator
therefore suppresses the sub-segment low-frequency content that our fixed-record periodogram
resolves, and the slow divider wander of § 7 would be largely invisible to it.

This is **inference and is tagged as such**: his paper does not state that the record was segmented,
does not give a segment length, and does not mention windowing or averaging over segments — only
that 100 independent *trajectories* were averaged. The one supporting detail in the paper is that
his Fig. 3(b) is plotted from f = 0.02 upward with ν₁ ≈ 0.085, i.e. displayed only above
≈ ν₁/4.25, which sits inside our flat range [ν_pred/6.25, ν_pred/2.08]; whether that is an axis
choice or where he searched is not stated. Our records are single fixed-length periodograms, which
resolve the wander into bins that can outgrow the resonance, so the ν_pred/2.5 floor of § 1 is
needed instead.

**Slow-mode origin (2026-09-24, B2 closed).** The canonical A1 v2 traces carry `Time, Wall_X,
Displacement, Left_Count, Right_Count, L0, eta, Center_X, Seed, Target_Oscillations,
Predicted_Frequency, Planned_Steps, Planned_Duration` — particle counts only, with no
per-compartment kinetic temperature and no release snapshot — so the regression of the slow-mode
variance fraction on (ΔT/T)² cannot be done over the canonical seeds and is not pursued further.
The mechanism is nevertheless directly evidenced by the dedicated 5-seed test in § 7.2, which logged
KE_L and KE_R: the slow mode correlates with the compartment temperature difference at r = +0.978,
and its amplitude is reproduced to 6–7 % by x = (L₀/2)(ΔT/T)·Z/(Z+ηZ′), with the expected
instantaneous imbalance √(2/N_s) ≈ 20 %. Origin is therefore the adiabatic-piston relaxation from
the release-time temperature imbalance, stated on that evidence. **`log T_1, T_2 at release` is added
to the KOA campaign list** so the regression can be done properly on the next campaign.
