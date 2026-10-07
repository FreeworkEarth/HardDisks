# Paper 1 — Methods

Draft, 2026-09-12; revised 2026-09-15 (estimator, A1 v2, A2 with the α = 2 cells, the slow divider mode).
Every command, grid and threshold below is transcribed from the `00_COMMAND.md` provenance files written by the
runs themselves, from the campaign drivers in `hspist3/validation/`, or from the CSV files named beside each
number. Nothing here is recalled. Where a number is not yet final it is marked **pending**.

---

## 1. Model and units

Two-dimensional hard disks of diameter σ and mass m in a rectangular box with hard
walls on all four sides, divided into two compartments by a single massive divider
that slides along x. Collisions are elastic and instantaneous; the dynamics is
event-driven (EDMD), so there is no integration timestep in the dynamics itself and no
energy drift by construction.

Reduced units throughout: **k_B T = m = σ = 1**. Lengths are in σ, velocities in
√(k_BT/m), times in σ√(m/k_BT). The packing fraction is

    η = N π r² / (2 L₀ H)

with r = σ/2 the disk radius, 2L₀ the full box width, H its height, and N the total
particle count. The divider sits at x = L₀ with N/2 particles each side.

All scientific runs use the reference EDMD backend, `--edmd-acc=0`. The accelerated
backend has not been validated and is not used for any number in this paper.

## 2. The sound-speed observable

The divider is held fixed while the gas equilibrates, then released. It executes a
damped oscillation in the compartment-pressure difference. For a piston of mass
M·m in a box of effective length L_eff = L₀ − 2r, the fundamental eigenmode is

    ν = c_s K / (2π L_eff),        cot K = α K,    α = M / (2 N_side)

with K the fundamental root, obtained by bisection. Measuring ν for a ladder of
divider masses and fitting ν against x = K/(2π L_eff) through the origin gives c_s as
the slope. The through-origin constraint is physical: ν must vanish as the piston
becomes infinitely heavy.

Theory curves are mapped from an equation of state through the adiabatic relation

    c_s² = (k_BT/m) [ Z + η Z′ + Z² ]

with Z(η) = P/(ρ k_BT). The reference is Kolafa–Rottner 2006, quoted only where it is
valid (η ≤ 0.69). Liu 2021's global equation of state is carried for the dense and
solid branches. The first-order virial truncation c_s ≈ √2(1+2η) was used in early
drafts and has been **removed from all figures**: it is 1.7 % below Kolafa–Rottner by
η = 0.08, several times the effect under test.

The divider carries a **second, slow mode** as well as this resonance; it is a physical
result of the paper and it is what forces the low-frequency edge in the estimator of § 6.
It is described in § 7.

## 3. Seeding, and the temperature it produces

Velocities are drawn Maxwell–Boltzmann, rescaled globally so that KE = N k_B T
exactly, then equalised per compartment. The order of the per-compartment
centre-of-mass drift removal matters and is selectable with `--seed-drift-order`:

- `old` (the historical default): rescale, then remove drift. Removing the drift after
  the rescale takes away ⟨½ N_s m v_cm²⟩ = k_B T per compartment, leaving
  **T_i = 1 − X/50 with X ~ Exp(1)**, i.e. about 2 % low and fluctuating.
- `drift-first`: remove the drift, recompute the compartment kinetic energy, then
  rescale. This gives **T_i = 1 exactly**.

Every production campaign in this paper now uses `drift-first`, so **no temperature
correction is applied to any figure**. The A1 campaign of 2026-08 used `old` and was
rerun as A1 v2 (§ 7). Verification for A1 v2: all **7875** `HD_KE_TRACE=1` audit lines read
`2b after per-segment equalize N=100 KE_tot=100 KE_left=50 KE_right=50 kT_mean=1.000000000`.

**Limit of the method.** Drift removal costs each compartment 2 of its 2·N_side
momentum degrees of freedom, a 1/N_side effect. At N_side = 1 it removes everything, so
one particle per side is unreachable with `drift-first` and is run with `old`, flagged.

## 4. Seeder geometry

Particles are placed on a lattice inset from the walls by a pad of 10⁻³·d, d the
particle diameter. The pad was previously `fmaxf(1e-4, 1e-5·d)` ≈ 2.4×10⁻⁴ px, which
is **below float resolution once the box exceeds ≈ 118 σ**: the outermost column
rounded onto the wall face and every such particle generated a `wall_overdue` event at
t = 0. That accounted for all 662 such events in the 2026-08 A1 campaign. The fixed inset
changes **seed positions only** — radii, box dimensions, particle counts and therefore η
are untouched — and is safe to a box of about 4×10⁵ px.

## 5. Health contract

An accepted production trajectory has **zero** numerical health events from
initialisation through the end of measurement: `forced_advance`, `wall_clamp_repairs`,
`overlap_repairs` and `wall_overdue` must all be 0. There is no tolerance and no
partial credit; a trajectory is used whole or discarded whole. A production seed is
never adapted and continued. The binary prints its `[EDMD-HEALTH]` line only when a
counter is nonzero, so a missing line is a positive statement that all four are zero.

Six distinct failure classes were caught by this layer and not by the observable:
event-budget exhaustion in production; an inert event calendar reporting Z = 1 as
valid; corner-packed seeding; wall-contact seeding; non-terminating random insertion
above η ≈ 0.55; and event-budget exhaustion during equilibration masked by a short
fresh-seed calibration. In each case the reported observable was numerically plausible.

Counts for the campaigns used in this paper: A1 v2 **2 discarded in 7875**
(`wall_clamp_repairs = 1`, at L₀ = 600 M = 100 and L₀ = 400 M = 2000); A2 famB **0 in 1000**;
A2 top-up **1 in 556**; α = 2 cells **0 in 40**; low-η extension **0 in 270**.

## 6. Frequency estimator

Per trajectory the record is the divider displacement from release to the end (no
transient drop), mean removed. With T the record length and ν_pred the predicted mode
frequency in the trace header, N_cyc = T·ν_pred is the record length in predicted periods:

    P(f) = |Σ_j (x(t_j) − x̄) e^(−2πi f t_j)|²,     ν_r = argmax_{f ≥ ν_pred/2.5} P(f)

i.e. the largest bin at or above k_min = round(N_cyc/2.5): 10 bins at 25 periods, 15 at
37.5 (A2), 80 at 200 (A1 v2). Bin width Δf/ν = 1/N_cyc, so 0.5 % at 200 periods and 2.7 %
at 37.5. Per cell the ν_r are averaged over seeds; per (η, N) the masses are fitted through
the origin, and the quoted error bar is the **1σ scatter of the per-mass implied c_s**.

**Why the edge exists.** The slow mode of § 7 puts power at low frequency, and the longer
the record, the more of it is resolved into bins that can outgrow the resonance bin. The
unrestricted largest-bin rule of Román 2002 is off by up to 33 % against a
[ν_pred/3, 3ν_pred] window at 25 periods and by up to 93 % at 200 periods. The edge is
one-sided, has no upper bound, and is placed with ν_pred, so the estimator is **not
theory-free**; what makes it safe is that it is flat against its own parameter.

**Sensitivity** (A1 v2, 35 η, 9 masses, 25 seeds; max over η of |Δc_s|/c_s against the
[ν_pred/3, 3ν_pred] window; `260914_A1v2_kmin_*.csv`):

| edge, as X = N_cyc/k_min | 25 periods | 200 periods |
|---|---|---|
| none (literal Román rule) | 33.339 % | 93.024 % |
| X = 100 … 25 | 7.467 % (X = 12.5) | 82.977 … 21.622 % |
| X = 12.5 | 0.000 % (k_min = 2 is X = 12.5 here) | 1.281 % |
| X = 6.25 | 1.817 % | 0.000 % |
| **X = 2.5 (adopted)** | **0.000 %** | **0.000 %** |
| X = 2.08 | 0.000 % | 0.000 % |

On A2 (23 (η, N) points, 37.5 periods) the adopted edge equals the window value at every
point. Rejected alternatives, reported and not used: the **velocity spectrum**
(2πf)²P(f), which needs no edge but shifts c_s by +0.02 to +3.3 %, one-signed and largest
at η ≥ 0.52; and fixed bin cuts, which are record-length dependent by construction.

**Cross-checks** (`260914_A1v2_pertraj_fit_and_windowed_peaks.csv`): a resonance fit to each
trajectory's spectrum and to the seed-averaged spectrum, both in a band around the
resonance, differ from the window value by at most 1.43 % for η ≤ 0.69.

**Damping bias.** The per-trajectory fits return Γ/f₀ = 0.023 (median, η < 0.5) and 0.063
(η ≥ 0.5), so the position peak lies below the eigenfrequency by Γ²/4f₀² ≤ 0.34 % over
η ≤ 0.69, and by 0.01–0.11 % typically. The measured velocity–position gap at η ≥ 0.5 is
+0.87 % median, which the damping shift does not explain; the fitted f₀ sits +0.47 %
above the position peak, between the two. The primary estimator stays the position-spectrum
maximum, for comparability with Román 2002, with the 0.34 % bound stated.

Cuts are never applied on agreement with another ν estimate, which would bias the slope.

## 7. The slow mode of the divider

Besides the acoustic resonance the divider position carries a slow component: a 5-period
running mean of x holds 9–57 % of the variance of x. It is **not** a start-up transient —
every run begins with exactly KE_L = KE_R = 50, its variance is the same in the second half
of a 1000-period record as in the first (ratio 0.72–1.25), and a 100× longer hold does not
change it — and it is not an oscillation but a relaxation, with correlation time 53–102
acoustic periods.

**Mechanism (adiabatic piston).** The acoustic mode is fast, so each compartment responds
adiabatically; heat crosses the divider only through divider collisions, which is slow.
Temperature fluctuations of the two 50-particle gases shift the pressure balance, the
divider follows, and it relaxes as heat flows. The fast mode sees the adiabatic stiffness
that sets c_s; the slow one sees the **isothermal** stiffness.

**Amplitude** (`260914_wander_equipartition_A1v2.csv`). With ⟨x²⟩ = L²/(2N_s[Z + ηZ′])
isothermal and L²/(2N_s[Z + ηZ′ + Z²]) adiabatic, N_s = 50 and Z from Kolafa–Rottner, the rms
of x about 0 over the last 800 of 1000 periods gives measured/isothermal = 0.962–1.044 for
light dividers at η = 0.196 and 0.524, while measured/adiabatic = 1.37–1.62. Heavy dividers
and 200-period records fall below the isothermal value because the slow mode starts from
T_L = T_R and relaxes more slowly the heavier the divider.

**Direct test** (`adiabatic_piston_check_20260913/`, 5 seeds, L₀ = 20, 200 periods, with the
per-compartment kinetic energies logged in the trace). The correlation of the slow part of
x with the slow part of T_L − T_R is **+0.978** (M = 100) and **+0.910** (M = 1000); the sign is
positive, a hotter left compartment displacing the divider to +x; at M = 100 the lag of the
maximum is ±0.25 periods. The pressure balance x = (L₀/2)(ΔT/T)·Z/(Z + ηZ′) reproduces the
measured amplitude to 6–7 %.

**Size scaling** (`260916_A2_slow_mode_vs_N.png`, A2 traces). In A2 the box grows as √N at
fixed η, so equipartition predicts a slow-mode amplitude that is flat in σ and falls as
N^(−1/2) in units of the box. Measured exponents: −0.03 … +0.11 in σ, and **−0.39 … −0.54**
for sd(x_slow)/L₀ across the five densities, i.e. the slow mode is a finite-size fluctuation
that vanishes as N → ∞.

**Resonance width.** Γ/f₀ measured on A2 does not follow the N^(−1/2) expected of acoustic
attenuation (exponents −0.02 … −0.25). At A2's 37.5-period records the spectral resolution
floor is 1/N_cyc = 0.027 and the measured Γ/f₀ is 0.033–0.13, i.e. at or near the floor over
most of the ladder, so A2 cannot resolve the width. The A1 v2 records (200 periods, floor
0.005, measured 0.017–0.17) can, and are the basis for the damping bound in § 6.

## 8. Campaign A1 v2 — L₀ sweep at fixed N

35 packing fractions from η = 0.006545 to 0.76. N = 100 (50/50), r = 0.5, H = 10; η is set by
the box length alone. Nine divider masses M = 50, 100, 200, 300, 500, 750, 1000, 1500, 2000;
25 seeds; 200 target oscillations; `drift-first`; one invocation per (η, M, seed) with
`--speed-sound-exact-seed`, so each trajectory is individually reproducible.

```sh
HD_KE_TRACE=1 ./00ALLINONE --mode=edmd --experiment=speed_of_sound --headless --kbt1 \
  --seed-drift-order=drift-first --edmd-acc=0 \
  --particles=100 --particles-boxes=50,50 --height=10.0 --particle-radius=0.5 \
  --wall-thickness=0.05 --wall-thickness-vis=0.05 --lengths=20.0000 --wall-masses=500 \
  --repeats=1 --seed=22260914 --wall-hold-steps=2000 --fixed-dt=0.4 \
  --target-oscillations=200 --oscillation-safety=1.0 \
  --oscillation-min-steps=10000 --oscillation-max-steps=400000000 \
  --speed-sound-log-stride=239 --speed-sound-run-dir=.../A1v2_20260914/eta_0p196350/m_500 \
  --speed-sound-exact-seed=417881510
```

**7875 trajectories, 2 discarded** (§ 5), T_i = 1 on all of them. The stride is set per cell to
about 32 samples per predicted period; halving it again changes the peak by 0.0000 %.
Results: `260914_A1v2_final_cs_vs_eta.csv`, figures `260914_cs_vs_eta`, `260914_cs_idealgas_zoom`.
At N = 100 the measured c_s sits +0.42 … +0.90 % above Kolafa–Rottner for η ≤ 0.08,
+1.2 … +2.0 % for 0.11 ≤ η ≤ 0.39, +3.0 % at η = 0.52, and returns to +0.02 % at η = 0.61.
§ 9 shows that this offset is finite size.

## 9. Campaign A2 — fixed η, N ladder

A1 measures c_s(η) at one system size and cannot separate the equation of state from the
finite box. A2 does: η is held fixed and the system is grown with **L₀ and H both scaling as
√N**, so the aspect ratio never changes and N is the only variable.

η = 0.10, 0.30, 0.50, 0.60, 0.65; N = 100, 400, 900, 1600, plus N = 2500 at η = 0.10 and 0.30;
masses M = 50, 200, 500, 1000, 2000, with M = 3200 at N = 1600 and M = 5000 at N = 2500 so that
α = M/N reaches Román's value 2 at the two largest sizes; 10–35 seeds per cell. N_side = N/2,
so α = M/N. Record length 37.5 periods.

```sh
HD_KE_TRACE=1 ./00ALLINONE --mode=edmd --experiment=speed_of_sound --headless --kbt1 \
  --seed-drift-order=drift-first --edmd-acc=0 \
  --particles=1600 --particles-boxes=800,800 --height=40.000000 --particle-radius=0.5 \
  --wall-thickness=0.05 --wall-thickness-vis=0.05 --lengths=157.079633 \
  --wall-masses=3200 --repeats=1 --seed=21860912 --wall-hold-steps=2000 --fixed-dt=0.4 \
  --target-oscillations=25 --oscillation-safety=1.5 \
  --oscillation-min-steps=10000 --oscillation-max-steps=80000000 \
  --speed-sound-log-stride=auto \
  --speed-sound-run-dir=.../A2_alpha2_20260912/eta_0p10/N1600/m_3200 \
  --speed-sound-exact-seed=<per-run>
```

**1595 trajectories used, 1 discarded** at the time of writing; the N = 2500 top-up cells hold
4–7 of their 10 seeds and are being completed (**pending**, § 12).

**Finite size.** c_s(N) is fitted per η with a + b/√N, a + b/N and a + b/√N + c/N, weighted by
the mass scatter (`260916_A2_finite_size_forms.csv`). The three forms agree to 0.8–1.7 % at
η ≤ 0.60 and bracket Kolafa–Rottner:

| η | c_∞ across forms | deviation from KR | best χ² (dof) |
|---|---|---|---|
| 0.10 | 1.7246 … 1.7545 | −1.12 … +0.59 % | 0.094 (2) |
| 0.30 | 2.8413 … 2.8741 | −0.36 … +0.79 % | 0.808 (2) |
| 0.50 | 5.4162 … 5.4575 | −0.00 … +0.76 % | 0.374 (1) |
| 0.60 | 8.2053 … 8.3127 | −0.21 … +1.10 % | 0.421 (1) |
| 0.65 | 9.8894 … 10.3175 | −5.26 … −1.16 % | 0.253 (1) |

The +1 … +3 % offset seen at N = 100 in A1 v2 therefore closes as the box grows: with the
finite-size form as the dominant uncertainty, **no departure from Kolafa–Rottner is claimed
for η ≤ 0.60**. At η = 0.65 all three forms sit below Kolafa–Rottner, by 1.2 % to 5.3 %; the
spread between forms (4.2 %) exceeds the statistical error of any one of them, so the size of
the deficit is not determined by these data, only its sign. The three-parameter form has one
degree of freedom on four sizes and is reported for the spread, not as a preferred fit.

## 10. Pressure campaign

An independent check of the same equation of state through a different observable.
The collisional virial gives Z_pair = 1 + W/(2·KE·dt) with W = Σ m|Δv_n|σ over pair
collisions. The wall momentum flux gives a **second, independent** estimator
Z_wall = (I_L+I_R)·W/(2·dt·KE); the two are never added, which would double-count, and
the x/y components of the wall route serve as an isotropy check.

Grid η = 0.005, 0.02, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.65, 0.67, 0.69 with
an exploratory tail at 0.698, 0.702, 0.710, 0.718, 0.720; N = 400, 900, 1600 (5, 4 and
3 seeds) and a later N = 2500 stage at η = 0.65, 0.67, 0.69. 30 blocks per trajectory.

```sh
./validation/pressure_validation "$eta" "$N" "$seed" "$NBLOCKS" "$bd" "$eq" \
    "$traj_csv" "$blk_csv" "$chunk"
```

The event budget is the limiting resource: the integrator caps calendar entries per
advance call and counts invalidated entries too, so the largest safe window shrinks
roughly as 1/N at fixed density and steeply with density. The runner calibrates the
chunk size on two disposable seeds per (η, N), takes 0.8× the smaller, and then runs
every scientific seed from a fresh state at that fixed chunk.

**Finite-size form uncertainty.** Extrapolating Z(N) with a + b/√N and with a + b/N
gives intercepts that agree only for η ≤ 0.10. At η = 0.65 the two forms give
8.3783 ± 0.0053 (χ² = 1.17) and 8.4067 ± 0.0029 (χ² = 0.42) on 2 degrees of freedom, a
spread of 0.0283 = 4.7× the combined statistical error, so **Z∞(0.65) = 8.3925 ± 0.0053
(stat) ± 0.0141 (form)** and no significant departure is claimed. At η = 0.67 and 0.69,
Z(N) is non-monotone and both forms are rejected by their own χ²; these data do not
support a bulk extrapolation with this model, and Z(N) is reported per size.

---

## 11. Figure captions

**Figure 1 — `260914_cs_vs_eta`.** Speed of sound against packing fraction over the full A1 v2
range, η = 0.0065 to 0.76, against Kolafa–Rottner 2006 (red, drawn only where valid), Liu 2021,
scaled-particle theory and Henderson. Points are the A1 v2 slope estimate with the estimator of
§ 6 on 200-period records; the error bar is the 1σ scatter of c_s across the nine divider masses.
No temperature correction is applied (T_i = 1, measured). Shaded bands mark the fluid,
coexistence, hexatic and solid regimes; no fluid-branch curve is drawn across coexistence, where
Liu's rigidity d(Zη)/dη is negative. Above η ≈ 0.65 the compartment is only ~6σ across and
structure changes during measurement, so those points are not fluid-branch values.

**Figure 2 — `260914_cs_idealgas_zoom`.** The dilute end, η ≤ 0.1, with the exact ideal-gas point
c_s = √2 at the origin. Lower panel: deviation from Kolafa–Rottner in per cent, with a ±0.5 %
band for scale. The offset is +0.42 … +0.90 % and does not close as η → 0 at fixed N = 100;
Figure 3 shows it is a finite-size effect.

**Figure 3 — `260916_A2_cs_vs_N`** (**pending** the N = 2500 top-up). Speed of sound against system
size at fixed η and fixed aspect ratio, one panel per η, against 1/√N. The red line is
Kolafa–Rottner; the orange square is the extrapolated c_∞ with its statistical error; the error
bar on each point is the 1σ scatter over divider masses. The N = 1600 and 2500 points at η = 0.10
and 0.30 include the α = 2 cells. Extrapolated to infinite size, η = 0.10 … 0.60 agree with
Kolafa–Rottner once the finite-size form uncertainty of § 9 is carried.

**Figure 4 — `Z_vs_eta` and `dev_vs_invsqrtN`** (pressure campaign). Compressibility factor against
packing fraction from the collisional virial, with the wall momentum flux as an independent
estimator, against Kolafa–Rottner. The companion panel shows the finite-size deviation against
1/√N per η with both extrapolation forms overlaid.

**Figure 5 — `260916_A2_slow_mode_vs_N`.** The slow divider mode against system size at fixed η and
aspect ratio. Upper row: sd of the 5-period running mean of x in units of L₀, log-log, with an
N^(−1/2) reference; the fitted exponents are −0.39 … −0.54, i.e. the slow mode is an equilibrium
finite-size fluctuation. Lower row: the fitted resonance width Γ/f₀, which is at the 1/N_cyc
resolution floor of these 37.5-period records over most of the ladder and is therefore not
resolved by A2 (§ 7).

---

## 12. Reproducibility

Every campaign leaf carries a `00_COMMAND.md` with the verbatim invocation, the working directory
and a timestamp; resumable cells additionally carry `00_RESUME.md` with the per-run seed, wall time
and health-line count of every run. Per-run seeds are a deterministic hash of the base seed and the
(length, mass, repeat) grid indices, and `--speed-sound-exact-seed` reproduces any single
trajectory. Analysis scripts live in `hspist3/validation/` and are read-only on all trajectory data.

Binary changes are gated by byte identity: the same speed-of-sound invocation through the old and
the new binary must produce identical trajectory CSVs before the new binary is used for science.
The current binary adds optional per-compartment kinetic-energy columns to the trace behind the
environment variable `HD_KE_COLUMNS` (§ 7). Four gates were run before it was installed: the old
binary reproducing an accepted A1 v2 trace; the unmodified source rebuilt; the new source with the
variable unset (all three byte-identical to the accepted trace); and the new source with the
variable set, whose first 13 columns are byte-identical and whose stdout differs only in the
output-folder line. The pre-change binary is kept alongside the new one.

---

## 13. Linewidth convention (fixed 2026-10-12; applies to every spectrum from now on)

**Model.** The divider-position autocorrelation is fitted as

$$C(t) = A\,e^{-t/\tau_T} + B\,e^{-t/\tau_r}\cos\omega_1 t ,$$

where $\tau_r$ is the **amplitude** decay time of the mode. This is the quantity every table up to
2026-10-12 reports, including 261006 §3.

**Linewidth.** $\Gamma$ is the **energy**-decay rate of the mode:

$$\Gamma \equiv \frac{2}{\tau_r}, \qquad E_1(t) \propto e^{-\Gamma t}, \qquad \sigma_\Gamma = \frac{2\,\sigma_{\tau_r}}{\tau_r^2}.$$

It is twice the amplitude rate $1/\tau_r$. Near $\omega_1$ the Fourier transform of the oscillatory
term is a Lorentzian,

$$S(\omega) \simeq \frac{B}{2}\,\frac{\Gamma}{(\omega-\omega_1)^2 + (\Gamma/2)^2},$$

so $\Gamma$ is also the **FWHM in angular frequency**. In ordinary frequency
$\Delta f_{\rm FWHM} = \Gamma/2\pi = 1/(\pi\tau_r)$, and
$Q = \omega_1/\Gamma = \pi\nu_1\tau_r = (\Delta f/f)^{-1}$. The $\Delta f/f$ of
`paper1_linewidth_20260918.py` (FWHM of the seed-averaged position spectrum over $f$) and the
"$\Delta f/f$ Mansour" column of 261006 §3 are therefore $\Gamma/\omega_1$. No conversion factor
is needed.

**Integrated peak power.** With $S(\omega) = \int C(t)\,e^{i\omega t}\,dt$ (two-sided), Parseval
gives $C(0) = \int S\,d\omega/2\pi$. The two Lorentzians at $\pm\omega_1$ carry $B/2$ each, so

$$P_1 \equiv \int_{\rm peak} S(\omega)\,\frac{d\omega}{2\pi} = B \quad [\sigma^2],$$

which is the variance of the divider position carried by the mode. The primary estimator is $B$
from the ACF fit. A periodogram check integrates the one-sided PSD over
$|f-\nu_1| \le 5\,\Delta f_{\rm FWHM}$ and divides by the Lorentzian fraction
$(2/\pi)\arctan 10 = 0.9365$. A periodogram FWHM includes the record-length resolution
broadening of one bin ($1/T_{\rm rec}$), so its raw width is an upper bound on $\Gamma/2\pi$.

**Not converted.** The "$\Gamma/f_0 = 0.023$ / 0.063" of § 6 (damping bias) came from a
per-trajectory resonance fit whose script is no longer on disk. Its width convention cannot be
read from code, so it is left as published and is not compared with any $\Gamma$ defined here.
The slow-mode time $\tau_T$ keeps its ACF meaning; the factor 2 applies to the oscillatory line
only.

---

## 14. Box-truncation correction of the canonical A1 v2 table (`final/260919_A1v2_final_cs_vs_eta.csv`) — pre-registration (2026-10-14)

*(The plan's "§9 of the A1v2 results markdown". No markdown holds the 260919 table itself. A1 v2 is documented in this file's § 8, so the section goes here, numbered 14 because §§ 9–13 exist. Written and committed before anything below was computed.)*

**The finding (261012 § 1.10, DATA).**
- The box width is set in whole pixels, `SIM_WIDTH = (int)(2 * L0_UNITS * PIXELS_PER_SIGMA)` (`00ALLINONE.c:323`, 24 px/σ). The physics uses that truncated box: `prm.boxW = (double)(XW2 - XW1)` (15882 speed-of-sound).
- The recorded $\eta$ (`eta_nominal_const`, 15800–15802) and the analysis length $L_{\rm eff} = L_0 - 2r - t/2$ (`tests_20260913.l_eff`) both use the untruncated $L_0$.
- The box shortfall is
$$\delta = \frac{2L_0 \cdot 24 - \lfloor 2L_0 \cdot 24\rfloor}{24}\ [\sigma].$$

**Per compartment it is $\delta/2$, not $\delta$ (from the runs' own output).** For $\eta = 0.1122$, $L_0 = 34.999901$ (`A1v2/eta_0p112200/m_500`, run 0, first row):
- `Center_X(σ)` = 43.3125, which is $(200 + 200 + 1679)/48$: the centre of the *truncated* box (15797).
- `Wall_X` = 43.333232, i.e. $200/24 + L_0$: the divider *starts* $L_0$ from the left wall.
- `Displacement(σ)` = +0.020732 $= \delta/2$.

So the left compartment starts at $L_0$ and the right at $L_0 - \delta$. With equal $N_s$ the divider oscillates about the truncated centre, where each compartment is $L_0 - \delta/2$. The script prints, per cell, `Center_X` minus $(200/24 + L_0)$, which should equal $-\delta/2$, and the trajectory-mean `Displacement`, which should be about 0 rather than $+\delta/2$. That confirms the mean geometry from the data. The plan's "$L_0 - \delta$" is the shortened *right* compartment at $t = 0$ only. It is printed as an upper-bound variant, not used for the verdict.

**The correction, per cell (primary):**
$$L_{0,\rm true} = L_0 - \frac{\delta}{2},\qquad \eta_{\rm true} = \eta_{\rm rec}\,\frac{L_0}{L_{0,\rm true}},\qquad L_{\rm eff,true} = L_{\rm eff,rec} - \frac{\delta}{2}.$$

**Check against the run.** $\eta_{\rm true}$ is checked against the recorded particle counts (`Left_Count`, `Right_Count`) and the height. The height is not written by speed-of-sound mode. It is read back from the run's own $\eta_{\rm rec} = N\pi r^2/(2L_0H)$ as $H = N\pi r^2/(2L_0\eta_{\rm rec})$ and compared with 10. `SIM_HEIGHT = (int)(H \cdot 24)` is exact for $H = 10$. Then $\eta_{\rm true} = N_s\pi r^2/(H\,L_{0,\rm true})$.

**$c_s$ re-derived from the recorded frequencies.** The canonical estimator (`paper1_populate_cs_err_20261002.py`, lines 117–120):

    L0 = float(leaf["L0"])
    x = np.array([T.x_of(q["M"], L0) for q in cs])
    y = np.array([q["nu"] for q in cs])

with `T.x_of(M, L0) = k_root(M/(2 N_SIDE)) / (2 pi l_eff(L0))` (`tests_20260913.py:182–183`), and $c_s$ the through-origin slope (`slope_with_errors`). The re-derivation keeps the per-mass $\nu$ (from `cell()` on the raw traces) and replaces $l_{\rm eff}(L_0)$ by $L_{\rm eff,true}$.

The mass ratio $\alpha$, and hence $K$, does not change. So $c_{s,\rm true} = c_{s,\rm rec}\,L_{\rm eff,true}/L_{\rm eff,rec}$ exactly, and the error rescales with it. That identity is printed as a check.

**Estimator gate first.** The recomputation from the raw traces must reproduce the 260919 table's $c_s$ to $5\times10^{-6}$ in every cell before any correction is reported. Otherwise the cell is VOID.

**Verdict rule.** With $D = (c_s - c_s^{\rm KR})/\sigma$, $\sigma$ = `c_s_err_scaled`, and $c_s^{\rm KR}$ from the table's own function at the respective $\eta$:

- **REGENERATE:** some canonical cell's $D$ changes by more than $0.5$ (of that cell's own error) between $(c_{s,\rm rec}, \eta_{\rm rec})$ and $(c_{s,\rm true}, \eta_{\rm true})$. Paper 1's table and figure must then be regenerated from the corrected values.
- **KEEP:** otherwise. The correction is documented and the table is unchanged.

The $\pi/8$ anchor ($L_0 = 10$) must come out with $\delta = 0$ exactly.

### 14.1 Results (the plan's "§ 9.1"; computed once, after the § 14 commit 4db8c9d)

**Printed by `python3 hspist3/validation/paper1_boxtrunc_20261014.py`**, verbatim:

#### Box-truncation correction, every canonical A1 v2 cell (sigma = c_s_err_scaled)

| eta_rec | L_0 | delta | Center_X - (XW1 + L_0) | <Displacement> (m_500) | H from eta_rec | eta_true | L_eff,rec | L_eff,true | c_s,rec ± σ | gate | c_s,true ± σ | KR(eta_rec) | KR(eta_true) | D before | D after | change | flag | change if L_0 - delta |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.006545 | 600.0 | 0.000000 | -0.000000 | -0.9402 | 9.99998 | 0.006545 | 598.9750 | 598.9750 | 1.44297 ± 0.00238 | ok | 1.44297 ± 0.00238 | 1.43290 | 1.43290 | +4.24 | +4.24 | +0.00 |  | +0.00 |
| 0.009817 | 400.0 | 0.000000 | -0.000000 | -1.2442 | 10.00049 | 0.009817 | 398.9750 | 398.9750 | 1.45388 ± 0.00434 | ok | 1.45388 ± 0.00434 | 1.44238 | 1.44238 | +2.65 | +2.65 | +0.00 |  | +0.00 |
| 0.013090 | 300.0 | 0.000000 | -0.000000 | -1.0503 | 9.99998 | 0.013090 | 298.9750 | 298.9750 | 1.46003 ± 0.00247 | ok | 1.46003 ± 0.00247 | 1.45195 | 1.45195 | +3.27 | +3.27 | +0.00 |  | +0.00 |
| 0.019635 | 199.9995 | 0.040649 | -0.020329 | +2.9489 | 10.00000 | 0.019637 | 198.9745 | 198.9542 | 1.47744 ± 0.00442 | ok | 1.47729 ± 0.00442 | 1.47138 | 1.47138 | +1.37 | +1.34 | -0.03 |  | -0.07 |
| 0.026180 | 149.9996 | 0.040873 | -0.020436 | -1.4855 | 10.00000 | 0.026184 | 148.9746 | 148.9542 | 1.50376 ± 0.00443 | ok | 1.50356 ± 0.00442 | 1.49118 | 1.49119 | +2.84 | +2.79 | -0.05 |  | -0.10 |
| 0.039270 | 99.9998 | 0.041260 | -0.020635 | +0.2237 | 10.00000 | 0.039278 | 98.9748 | 98.9542 | 1.54272 ± 0.00419 | ok | 1.54240 ± 0.00419 | 1.53196 | 1.53199 | +2.57 | +2.48 | -0.08 |  | -0.16 |
| 0.052360 | 74.9998 | 0.041270 | -0.020635 | +0.0826 | 10.00000 | 0.052374 | 73.9748 | 73.9542 | 1.58740 ± 0.00257 | ok | 1.58696 ± 0.00257 | 1.57439 | 1.57444 | +5.07 | +4.88 | -0.19 |  | -0.38 |
| 0.078540 | 49.9999 | 0.041463 | -0.020734 | -0.3284 | 10.00000 | 0.078573 | 48.9749 | 48.9542 | 1.67866 ± 0.00526 | ok | 1.67795 ± 0.00526 | 1.66454 | 1.66465 | +2.68 | +2.53 | -0.16 |  | -0.31 |
| 0.112200 | 34.9999 | 0.041468 | -0.020734 | +0.0221 | 10.00000 | 0.112267 | 33.9749 | 33.9542 | 1.81155 ± 0.00529 | ok | 1.81044 ± 0.00529 | 1.79189 | 1.79216 | +3.71 | +3.46 | -0.26 |  | -0.51 |
| 0.130900 | 29.9999 | 0.041468 | -0.020734 | +0.0704 | 10.00001 | 0.130991 | 28.9749 | 28.9542 | 1.89683 ± 0.00473 | ok | 1.89547 ± 0.00473 | 1.86885 | 1.86924 | +5.91 | +5.55 | -0.36 |  | -0.73 |
| 0.157080 | 24.9999 | 0.041468 | -0.020734 | +0.1151 | 10.00002 | 0.157210 | 23.9749 | 23.9542 | 2.01616 ± 0.00277 | ok | 2.01442 ± 0.00277 | 1.98495 | 1.98555 | +11.25 | +10.41 | -0.84 | **> 0.5** | -1.68 |
| 0.196350 | 20.0 | 0.000000 | -0.000000 | -0.0053 | 9.99998 | 0.196350 | 18.9750 | 18.9750 | 2.21096 ± 0.00555 | ok | 2.21096 ± 0.00555 | 2.17981 | 2.17981 | +5.62 | +5.62 | -0.00 |  | -0.00 |
| 0.261799 | 15.0 | 0.000000 | -0.000000 | -0.0750 | 10.00001 | 0.261799 | 13.9750 | 13.9750 | 2.61299 ± 0.00766 | ok | 2.61299 ± 0.00766 | 2.57255 | 2.57255 | +5.28 | +5.28 | +0.00 |  | +0.00 |
| 0.392699 | 10.0 | 0.000000 | -0.000000 | -0.0001 | 10.00000 | 0.392699 | 8.9750 | 8.9750 | 3.80884 ± 0.01119 | ok | 3.80884 ± 0.01119 | 3.74608 | 3.74608 | +5.61 | +5.61 | -0.00 |  | -0.00 |
| 0.523599 | 7.5 | 0.000000 | -0.000000 | -0.0039 | 10.00000 | 0.523599 | 6.4750 | 6.4750 | 6.08815 ± 0.01713 | ok | 6.08815 ± 0.01713 | 5.93213 | 5.93213 | +9.11 | +9.11 | +0.00 |  | +0.00 |
| 0.549999 | 7.14 | 0.030000 | -0.015000 | -0.0053 | 9.99999 | 0.551157 | 6.1150 | 6.1000 | 6.73988 ± 0.03585 | ok | 6.72335 ± 0.03576 | 6.60185 | 6.63378 | +3.85 | +2.50 | -1.35 | **> 0.5** | -2.71 |
| 0.569996 | 6.8895 | 0.029001 | -0.014500 | +0.0203 | 10.00001 | 0.571198 | 5.8645 | 5.8500 | 7.21712 ± 0.05177 | ok | 7.19927 ± 0.05164 | 7.18766 | 7.22533 | +0.57 | -0.50 | -1.07 | **> 0.5** | -2.16 |
| 0.590001 | 6.6559 | 0.020134 | -0.010066 | -0.0067 | 10.00001 | 0.590895 | 5.6309 | 5.6208 | 7.92513 ± 0.09137 | ok | 7.91096 ± 0.09121 | 7.85471 | 7.88661 | +0.77 | +0.27 | -0.50 | **> 0.5** | -1.01 |
| 0.609999 | 6.4377 | 0.000399 | -0.000200 | -0.0032 | 10.00000 | 0.610018 | 5.4127 | 5.4125 | 8.57780 ± 0.06569 | ok | 8.57748 ± 0.06569 | 8.61562 | 8.61638 | -0.58 | -0.59 | -0.02 |  | -0.03 |
| 0.630002 | 6.2333 | 0.008268 | -0.004133 | +0.0107 | 10.00000 | 0.630420 | 5.2083 | 5.2042 | 9.51057 ± 0.12953 | ok | 9.50302 ± 0.12943 | 9.48076 | 9.49995 | +0.23 | +0.02 | -0.21 |  | -0.41 |
| 0.650003 | 6.0415 | 0.041334 | -0.020666 | -0.0046 | 9.99999 | 0.652234 | 5.0165 | 4.9958 | 11.03078 ± 0.05298 | ok | 10.98533 ± 0.05276 | 10.43889 | 10.54877 | +11.17 | +8.27 | -2.90 | **> 0.5** | -5.83 |
| 0.669998 | 5.8612 | 0.014066 | -0.007033 | +0.0064 | 10.00000 | 0.670803 | 4.8362 | 4.8292 | 12.99552 ± 0.03927 | ok | 12.97662 ± 0.03922 | 11.37328 | 11.40491 | +41.31 | +40.08 | -1.23 | **> 0.5** | -2.44 |
| 0.679998 | 5.775 | 0.008334 | -0.004166 | +0.0025 | 10.00001 | 0.680489 | 4.7500 | 4.7458 | 14.51554 ± 0.04912 | ok | 14.50280 ± 0.04907 | 11.67118 | 11.67835 | +57.91 | +57.56 | -0.35 |  | -0.69 |
| 0.689999 | 5.6913 | 0.007600 | -0.003800 | +0.0028 | 10.00000 | 0.690460 | 4.6663 | 4.6625 | 16.51931 ± 0.10955 | ok | 16.50586 ± 0.10946 | 11.55102 | 11.52719 | +45.35 | +45.48 | +0.13 |  | +0.28 |
| 0.695006 | 5.6503 | 0.008934 | -0.004466 | -0.0004 | 10.00000 | 0.695556 | 4.6253 | 4.6208 | 17.76362 ± 0.15152 | ok | 17.74646 ± 0.15137 | nan | nan | +nan | +nan | +nan | no KR in table | +nan |
| 0.699998 | 5.61 | 0.011667 | -0.005833 | +0.0016 | 10.00001 | 0.700727 | 4.5850 | 4.5792 | 19.24481 ± 0.22622 | ok | 19.22032 ± 0.22593 | nan | nan | +nan | +nan | +nan | no KR in table | +nan |
| 0.705000 | 5.5702 | 0.015400 | -0.007700 | +0.0283 | 10.00000 | 0.705976 | 4.5452 | 4.5375 | 18.38705 ± 0.67519 | ok | 18.35590 ± 0.67405 | nan | nan | +nan | +nan | +nan | no KR in table | +nan |
| 0.709997 | 5.531 | 0.020334 | -0.010166 | +0.0011 | 9.99999 | 0.711304 | 4.5060 | 4.4958 | 17.11137 ± 0.51285 | ok | 17.07276 ± 0.51169 | nan | nan | +nan | +nan | +nan | no KR in table | +nan |
| 0.714999 | 5.4923 | 0.026267 | -0.013133 | +0.0218 | 10.00000 | 0.716713 | 4.4673 | 4.4542 | 16.44125 ± 0.19129 | ok | 16.39291 ± 0.19073 | nan | nan | +nan | +nan | +nan | no KR in table | +nan |
| 0.719994 | 5.4542 | 0.033399 | -0.016700 | +0.0114 | 10.00000 | 0.722205 | 4.4292 | 4.4125 | 16.75907 ± 0.23220 | ok | 16.69588 ± 0.23133 | nan | nan | +nan | +nan | +nan | no KR in table | +nan |
| 0.725005 | 5.4165 | 0.041334 | -0.020666 | +0.0125 | 10.00000 | 0.727782 | 4.3915 | 4.3708 | 17.86572 ± 0.13509 | ok | 17.78164 ± 0.13446 | nan | nan | +nan | +nan | +nan | no KR in table | +nan |
| 0.730005 | 5.3794 | 0.008799 | -0.004400 | +0.0027 | 10.00000 | 0.730603 | 4.3544 | 4.3500 | 18.58431 ± 0.09869 | ok | 18.56553 ± 0.09859 | nan | nan | +nan | +nan | +nan | no KR in table | +nan |
| 0.740006 | 5.3067 | 0.030067 | -0.015033 | +0.0078 | 10.00000 | 0.742108 | 4.2817 | 4.2667 | 22.76518 ± 0.13520 | ok | 22.68524 ± 0.13472 | nan | nan | +nan | +nan | +nan | no KR in table | +nan |
| 0.749998 | 5.236 | 0.013667 | -0.006833 | +0.0049 | 10.00000 | 0.750978 | 4.2110 | 4.2042 | 27.84595 ± 0.21448 | ok | 27.80076 ± 0.21414 | nan | nan | +nan | +nan | +nan | no KR in table | +nan |
| 0.759999 | 5.1671 | 0.000867 | -0.000433 | +0.0003 | 10.00000 | 0.760063 | 4.1421 | 4.1417 | 37.40656 ± 0.19181 | ok | 37.40264 ± 0.19179 | nan | nan | +nan | +nan | +nan | no KR in table | +nan |

pi/8 anchor: delta = 0.000000 (must be exactly 0) -> OK
estimator gate: 35/35 cells reproduce the 260919 c_s to 5e-6; KR function reproduces the table's KR column to 4.9e-06 (24 cells with a tabulated KR; the rest carry no published deviation and do not enter the verdict)
identity c_s,true/c_s,rec = L_eff,true/L_eff,rec: max deviation 3.3e-16
recorded Center_X offset vs -delta/2: max |difference| 5.5e-06 sigma

cells with |change| > 0.5: 6 -> **VERDICT: REGENERATE** (eta = 0.157, 0.550, 0.570, 0.590, 0.650, 0.670)
largest |change|: 2.90 at eta = 0.6500; for eta <= 0.39: 0.84

Figure: `paper1_speedofsound/experiments/final/261014_p1_boxtrunc_shift.{png,pdf}`.

**Verdict, by the pre-registered rule: REGENERATE (DATA).** Six cells change $D$ by more than 0.5:

| $\eta$ | 0.157 | 0.550 | 0.570 | 0.590 | 0.650 | 0.670 |
|---|---|---|---|---|---|---|
| change in $D$ | −0.84 | −1.35 | −1.07 | −0.50 | −2.90 | −1.23 |

At $\eta = 0.590$ the change sits on the 0.5 line.

**What the correction does.**
- **Direction.** Every change is negative. The correction lowers $c_s$ (a shorter $L_{\rm eff}$) and raises $c_s^{\rm KR}$ (a higher $\eta$), so the deviation from KR shrinks everywhere.
- **Dilute side.** At $\eta \le 0.13$ the changes are at most 0.36σ, and the +1–2 % N = 100 excess of § 8 is essentially untouched. For example, $\eta = 0.1122$ goes from $D = +3.71$ to $+3.46$.
- **Exact cells.** The $\pi/8$ anchor has $\delta = 0$ exactly. So do the cells with $L_0 = 600, 400, 300, 20, 15, 10, 7.5$.
- **Where it matters.** The large shifts are at the dense cells, where both $L_{\rm eff}$ is short and $c_s^{\rm KR}$ is steep in $\eta$.

**Data checks (DATA).**
- **Truncated geometry, confirmed from the runs' own output.** The recorded `Center_X` sits at $-\delta/2$ from $XW_1/24 + L_0$ in every cell, to $5.5\times10^{-6}$ σ, i.e. at print precision.
- **The divider oscillates about the truncated centre.** At dense $\eta$ the trajectory-mean `Displacement` is $0.000$–$0.03$, not $+\delta/2$. At dilute $\eta$ the divider's slow wander (±1–3 σ at $L_0 = 600$) makes this check uninformative.
- **Height and counts.** $H$ read back from $\eta_{\rm rec}$ is $10.0000 \pm 5\times10^{-5}$ in every cell, with $N = 100$.

**Gates.**
- The estimator reproduces all 35 cells to $5\times10^{-6}$.
- The KR function reproduces the table's KR column to $4.9\times10^{-6}$ in the 24 cells that tabulate it.
- The 11 cells with $\eta \ge 0.695$ carry no KR in the table, hence no published deviation. They do not enter the verdict, though their $c_{s,\rm true}$ is printed.
- The rescaling identity holds to $3\times10^{-16}$.

**The plan's "$L_0 - \delta$"** (the right compartment at $t = 0$ only) would roughly double every change; see the last column. It is a bound, not the correction.

**What follows.** Paper 1's canonical table and figure must be regenerated with $\eta_{\rm true}$ and $L_{\rm eff,true}$. This is a recomputation from data that already exist; no rerun is needed, and the script above already computes every $c_{s,\rm true}$. It is **not done in this batch**, because the canonical table and the drafts change only on the plan author's go.

**Recording fix in the binary (output only), and its determinism gate.**

The code that writes or warns (`00ALLINONE.c`):
- **Warning**, line 330 in `initialize_simulation_dimensions()`: `fprintf(stderr, "WARNING: L_0 not on the 1/48 sigma grid: box truncated by %.6f sigma\n", ...)`. The same check is applied to $H$ on the 1/24 grid.
- **Speed-of-sound trace**: the header gets a new last column at 15796, `fputs(",Box_Width_sigma", wall_log);`. Its value is `box_width_sigma_const = (XW2 - XW1) / PIXELS_PER_SIGMA` (15810), written by both row writers at 16103 and 16202.
- **Energy-transfer summary**: the header gains `...,build_cflags,box_width_sigma` (17283), and the row its value (17473).

`--version` before was `git 05215ea-dirty target release`, CFLAGS `-O3 -march=native -ffp-contract=off`. After it was `git 5190846-dirty target release` with the same CFLAGS. The previous binary is kept as `hspist3/00ALLINONE_pre_boxw_20261014` (untracked).

**Gate.** The same seeds went through both binaries in both modes, on the grid-exact $\pi/8$ cell and the non-grid $L_0 = 34.9999$ cell, with outputs written to the scratchpad only:

| file | old vs new | note |
|---|---|---|
| SoS pi8 trace (1252 lines) | IDENTICAL after stripping the last column `Box_Width_sigma` | Box_Width_sigma = ['20.000000'] |
| SoS pi8 psi6 csv | IDENTICAL | |
| SoS nongrid trace (8330 lines) | IDENTICAL after stripping the last column `Box_Width_sigma` | Box_Width_sigma = ['69.958333'] |
| SoS nongrid psi6 csv | IDENTICAL | |
| ET pi8 tr.csv | IDENTICAL | byte compare |
| ET pi8 ev.csv | IDENTICAL | byte compare |
| ET pi8 summary | IDENTICAL on all fields except ['box_width_sigma', 'build_git', 'command', 'timestamp', 'trace_path'] | box_width_sigma = 20.000000 |
| ET nongrid tr.csv | IDENTICAL | byte compare (partial trace of the aborted run) |
| ET nongrid ev.csv | IDENTICAL | byte compare (partial trace of the aborted run) |
| ET nongrid summary | none written by EITHER binary -- both abort identically (same message) | `ABORTING INVALID RUN [initial_wall_position_mismatch]: Wall 0 initialized at 34.9791667 sigma; command requested 34.9999008 sigma.` |
- old sos_pi8: no warning
- old sos_nongrid: no warning
- old et_pi8: no warning
- old et_nongrid: no warning
- new sos_pi8: no warning
- new sos_nongrid: ['WARNING: L_0 not on the 1/48 sigma grid: box truncated by 0.041468 sigma']
- new et_pi8: no warning
- new et_nongrid: ['WARNING: L_0 not on the 1/48 sigma grid: box truncated by 0.041468 sigma']

**DETERMINISM GATE: PASS**

ENERGY-TRANSFER MODE ALREADY REFUSED NON-GRID GEOMETRY: both binaries abort the $L_0 = 34.9999$ run with `initial_wall_position_mismatch` (the divider snaps to the pixel grid, 34.9791667 σ). That is the "grid check" the Paper 2 run scripts relied on. Speed-of-sound mode had no such check, which is why only Paper 1 data are affected.

### 14.2 Regenerated table (2026-10-14; GO from the plan author on the REGENERATE verdict of § 14.1)

The full per-cell correction was printed **before any canonical file was changed**. Columns: recorded and true packing fraction, recorded and true acoustic length, $c_s$ before and after (σ = `c_s_err_scaled`, which rescales with $L_{\rm eff}$), KR at both packing fractions, and the deviation $D = (c_s - c_s^{\rm KR})/\sigma$ before and after. Rows with $\eta_{\rm rec} \ge 0.695$ carry no KR in the canonical table, so $D$ is n/a there.

**Printed by `python3 hspist3/validation/paper1_boxtrunc_20261014.py --table`** (verbatim):

#### Regenerated per-cell table (uncorrected input: 260919_A1v2_final_cs_vs_eta.csv; sigma = c_s_err_scaled, rescaled with L_eff)

| eta_rec | L_0 | delta | eta_true | L_eff,rec | L_eff,true | c_s,rec ± σ | c_s,true ± σ | KR(eta_rec) | KR(eta_true) | D before [σ] | D after [σ] |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.006545 | 600.0 | 0.000000 | 0.006545 | 598.9750 | 598.9750 | 1.44297 ± 0.00238 | 1.44297 ± 0.00238 | 1.43290 | 1.43290 | +4.24 | +4.24 |
| 0.009817 | 400.0 | 0.000000 | 0.009817 | 398.9750 | 398.9750 | 1.45388 ± 0.00434 | 1.45388 ± 0.00434 | 1.44238 | 1.44238 | +2.65 | +2.65 |
| 0.013090 | 300.0 | 0.000000 | 0.013090 | 298.9750 | 298.9750 | 1.46003 ± 0.00247 | 1.46003 ± 0.00247 | 1.45195 | 1.45195 | +3.27 | +3.27 |
| 0.019635 | 199.9995 | 0.040649 | 0.019637 | 198.9745 | 198.9542 | 1.47744 ± 0.00442 | 1.47729 ± 0.00442 | 1.47138 | 1.47138 | +1.37 | +1.34 |
| 0.026180 | 149.9996 | 0.040873 | 0.026184 | 148.9746 | 148.9542 | 1.50376 ± 0.00443 | 1.50356 ± 0.00442 | 1.49118 | 1.49119 | +2.84 | +2.79 |
| 0.039270 | 99.9998 | 0.041260 | 0.039278 | 98.9748 | 98.9542 | 1.54272 ± 0.00419 | 1.54240 ± 0.00419 | 1.53196 | 1.53199 | +2.57 | +2.48 |
| 0.052360 | 74.9998 | 0.041270 | 0.052374 | 73.9748 | 73.9542 | 1.58740 ± 0.00257 | 1.58696 ± 0.00257 | 1.57439 | 1.57444 | +5.07 | +4.88 |
| 0.078540 | 49.9999 | 0.041463 | 0.078573 | 48.9749 | 48.9542 | 1.67866 ± 0.00526 | 1.67795 ± 0.00526 | 1.66454 | 1.66465 | +2.68 | +2.53 |
| 0.112200 | 34.9999 | 0.041468 | 0.112267 | 33.9749 | 33.9542 | 1.81155 ± 0.00529 | 1.81044 ± 0.00529 | 1.79189 | 1.79216 | +3.71 | +3.46 |
| 0.130900 | 29.9999 | 0.041468 | 0.130991 | 28.9749 | 28.9542 | 1.89683 ± 0.00473 | 1.89547 ± 0.00473 | 1.86885 | 1.86924 | +5.91 | +5.55 |
| 0.157080 | 24.9999 | 0.041468 | 0.157210 | 23.9749 | 23.9542 | 2.01616 ± 0.00277 | 2.01442 ± 0.00277 | 1.98495 | 1.98555 | +11.25 | +10.41 |
| 0.196350 | 20.0 | 0.000000 | 0.196350 | 18.9750 | 18.9750 | 2.21096 ± 0.00555 | 2.21096 ± 0.00555 | 2.17981 | 2.17981 | +5.62 | +5.62 |
| 0.261799 | 15.0 | 0.000000 | 0.261799 | 13.9750 | 13.9750 | 2.61299 ± 0.00766 | 2.61299 ± 0.00766 | 2.57255 | 2.57255 | +5.28 | +5.28 |
| 0.392699 | 10.0 | 0.000000 | 0.392699 | 8.9750 | 8.9750 | 3.80884 ± 0.01119 | 3.80884 ± 0.01119 | 3.74608 | 3.74608 | +5.61 | +5.61 |
| 0.523599 | 7.5 | 0.000000 | 0.523599 | 6.4750 | 6.4750 | 6.08815 ± 0.01713 | 6.08815 ± 0.01713 | 5.93213 | 5.93213 | +9.11 | +9.11 |
| 0.549999 | 7.14 | 0.030000 | 0.551157 | 6.1150 | 6.1000 | 6.73988 ± 0.03585 | 6.72335 ± 0.03576 | 6.60185 | 6.63378 | +3.85 | +2.50 |
| 0.569996 | 6.8895 | 0.029001 | 0.571198 | 5.8645 | 5.8500 | 7.21712 ± 0.05177 | 7.19927 ± 0.05164 | 7.18766 | 7.22533 | +0.57 | -0.50 |
| 0.590001 | 6.6559 | 0.020134 | 0.590895 | 5.6309 | 5.6208 | 7.92513 ± 0.09137 | 7.91096 ± 0.09121 | 7.85471 | 7.88661 | +0.77 | +0.27 |
| 0.609999 | 6.4377 | 0.000399 | 0.610018 | 5.4127 | 5.4125 | 8.57780 ± 0.06569 | 8.57748 ± 0.06569 | 8.61562 | 8.61638 | -0.58 | -0.59 |
| 0.630002 | 6.2333 | 0.008268 | 0.630420 | 5.2083 | 5.2042 | 9.51057 ± 0.12953 | 9.50302 ± 0.12943 | 9.48076 | 9.49995 | +0.23 | +0.02 |
| 0.650003 | 6.0415 | 0.041334 | 0.652234 | 5.0165 | 4.9958 | 11.03078 ± 0.05298 | 10.98533 ± 0.05276 | 10.43889 | 10.54877 | +11.17 | +8.27 |
| 0.669998 | 5.8612 | 0.014066 | 0.670803 | 4.8362 | 4.8292 | 12.99552 ± 0.03927 | 12.97662 ± 0.03922 | 11.37328 | 11.40491 | +41.31 | +40.08 |
| 0.679998 | 5.775 | 0.008334 | 0.680489 | 4.7500 | 4.7458 | 14.51554 ± 0.04912 | 14.50280 ± 0.04907 | 11.67118 | 11.67835 | +57.91 | +57.56 |
| 0.689999 | 5.6913 | 0.007600 | 0.690460 | 4.6663 | 4.6625 | 16.51931 ± 0.10955 | 16.50586 ± 0.10946 | 11.55102 | 11.52719 | +45.35 | +45.48 |
| 0.695006 | 5.6503 | 0.008934 | 0.695556 | 4.6253 | 4.6208 | 17.76362 ± 0.15152 | 17.74646 ± 0.15137 | n/a | n/a | n/a | n/a |
| 0.699998 | 5.61 | 0.011667 | 0.700727 | 4.5850 | 4.5792 | 19.24481 ± 0.22622 | 19.22032 ± 0.22593 | n/a | n/a | n/a | n/a |
| 0.705000 | 5.5702 | 0.015400 | 0.705976 | 4.5452 | 4.5375 | 18.38705 ± 0.67519 | 18.35590 ± 0.67405 | n/a | n/a | n/a | n/a |
| 0.709997 | 5.531 | 0.020334 | 0.711304 | 4.5060 | 4.4958 | 17.11137 ± 0.51285 | 17.07276 ± 0.51169 | n/a | n/a | n/a | n/a |
| 0.714999 | 5.4923 | 0.026267 | 0.716713 | 4.4673 | 4.4542 | 16.44125 ± 0.19129 | 16.39291 ± 0.19073 | n/a | n/a | n/a | n/a |
| 0.719994 | 5.4542 | 0.033399 | 0.722205 | 4.4292 | 4.4125 | 16.75907 ± 0.23220 | 16.69588 ± 0.23133 | n/a | n/a | n/a | n/a |
| 0.725005 | 5.4165 | 0.041334 | 0.727782 | 4.3915 | 4.3708 | 17.86572 ± 0.13509 | 17.78164 ± 0.13446 | n/a | n/a | n/a | n/a |
| 0.730005 | 5.3794 | 0.008799 | 0.730603 | 4.3544 | 4.3500 | 18.58431 ± 0.09869 | 18.56553 ± 0.09859 | n/a | n/a | n/a | n/a |
| 0.740006 | 5.3067 | 0.030067 | 0.742108 | 4.2817 | 4.2667 | 22.76518 ± 0.13520 | 22.68524 ± 0.13472 | n/a | n/a | n/a | n/a |
| 0.749998 | 5.236 | 0.013667 | 0.750978 | 4.2110 | 4.2042 | 27.84595 ± 0.21448 | 27.80076 ± 0.21414 | n/a | n/a | n/a | n/a |
| 0.759999 | 5.1671 | 0.000867 | 0.760063 | 4.1421 | 4.1417 | 37.40656 ± 0.19181 | 37.40264 ± 0.19179 | n/a | n/a | n/a | n/a |

estimator gate: 35/35 cells reproduce the uncorrected table to 5e-6

**Box HEIGHT (00ALLINONE.c:324, `SIM_HEIGHT = (int)(HEIGHT_UNITS * PIXELS_PER_SIGMA);`).** Every A1 v2 run was launched by the harness with `--height=10.0` (tests_20260913.py:78 `H = "10.0"`, passed at :285 as `f"--height={H}"`):
H x 24 = [240.0] -> integer in all 35 cells, so SIM_HEIGHT is exact and the height is NOT truncated. Read back from each run's own eta_rec: H = 9.99998 ... 10.00049 (H x 24 = 239.999 ... 240.012; 6-decimal eta print).

**OPEN, not corrected in this batch: A2 (the finite-size ladder) uses non-grid L_0 as well.** Per (eta, N), from the A2 per-mass tables the draft's zoom and overlay figures read:

| table | eta | N | L_0 | delta [σ] | delta/2 / L_eff (c_s shift) | eta shift |
|---|---|---|---|---|---|---|
| A2 | 0.02 | 100 | 196.349548 | 0.032430 | 0.0083 % | 0.0083 % |
| A2 | 0.02 | 400 | 392.699097 | 0.023193 | 0.0030 % | 0.0030 % |
| A2 | 0.02 | 900 | 589.048645 | 0.013997 | 0.0012 % | 0.0012 % |
| A2 | 0.02 | 1600 | 785.398193 | 0.004720 | 0.0003 % | 0.0003 % |
| A2 | 0.05 | 100 | 78.539818 | 0.037964 | 0.0245 % | 0.0242 % |
| A2 | 0.05 | 400 | 157.079636 | 0.034261 | 0.0110 % | 0.0109 % |
| A2 | 0.05 | 900 | 235.619446 | 0.030558 | 0.0065 % | 0.0065 % |
| A2 | 0.05 | 1600 | 314.159271 | 0.026855 | 0.0043 % | 0.0043 % |
| A2 | 0.1 | 100 | 39.269909 | 0.039815 | 0.0521 % | 0.0507 % |
| A2 | 0.1 | 400 | 78.539818 | 0.037964 | 0.0245 % | 0.0242 % |
| A2 | 0.1 | 900 | 117.809723 | 0.036112 | 0.0155 % | 0.0153 % |
| A2 | 0.1 | 1600 | 157.079636 | 0.034261 | 0.0110 % | 0.0109 % |
| A2 | 0.1 | 2500 | 196.349548 | 0.032430 | 0.0083 % | 0.0083 % |
| A2 | 0.3 | 100 | 13.089969 | 0.013270 | 0.0550 % | 0.0507 % |
| A2 | 0.3 | 400 | 26.179939 | 0.026545 | 0.0528 % | 0.0507 % |
| A2 | 0.3 | 900 | 39.269909 | 0.039815 | 0.0521 % | 0.0507 % |
| A2 | 0.3 | 1600 | 52.359879 | 0.011424 | 0.0111 % | 0.0109 % |
| A2 | 0.3 | 2500 | 65.449844 | 0.024689 | 0.0192 % | 0.0189 % |
| A2 | 0.5 | 100 | 7.853982 | 0.041298 | 0.3024 % | 0.2636 % |
| A2 | 0.5 | 400 | 15.707963 | 0.040927 | 0.1394 % | 0.1304 % |
| A2 | 0.5 | 900 | 23.561945 | 0.040558 | 0.0900 % | 0.0861 % |
| A2 | 0.5 | 1600 | 31.415928 | 0.040192 | 0.0661 % | 0.0640 % |
| A2 | 0.6 | 100 | 6.544985 | 0.006636 | 0.0601 % | 0.0507 % |
| A2 | 0.6 | 400 | 13.089969 | 0.013270 | 0.0550 % | 0.0507 % |
| A2 | 0.6 | 900 | 19.634954 | 0.019908 | 0.0535 % | 0.0507 % |
| A2 | 0.6 | 1600 | 26.179939 | 0.026545 | 0.0528 % | 0.0507 % |
| A2 | 0.65 | 100 | 6.041524 | 0.041382 | 0.4125 % | 0.3437 % |
| A2 | 0.65 | 400 | 12.083049 | 0.041097 | 0.1858 % | 0.1704 % |
| A2 | 0.65 | 900 | 18.124573 | 0.040812 | 0.1193 % | 0.1127 % |
| A2 | 0.65 | 1600 | 24.166098 | 0.040527 | 0.0876 % | 0.0839 % |
| A2_famB | 0.02 | 100 | 196.349548 | 0.032430 | 0.0083 % | 0.0083 % |
| A2_famB | 0.02 | 400 | 392.699097 | 0.023193 | 0.0030 % | 0.0030 % |
| A2_famB | 0.02 | 900 | 589.048645 | 0.013997 | 0.0012 % | 0.0012 % |
| A2_famB | 0.02 | 1600 | 785.398193 | 0.004720 | 0.0003 % | 0.0003 % |
| A2_famB | 0.05 | 100 | 78.539818 | 0.037964 | 0.0245 % | 0.0242 % |
| A2_famB | 0.05 | 400 | 157.079636 | 0.034261 | 0.0110 % | 0.0109 % |
| A2_famB | 0.05 | 900 | 235.619446 | 0.030558 | 0.0065 % | 0.0065 % |
| A2_famB | 0.05 | 1600 | 314.159271 | 0.026855 | 0.0043 % | 0.0043 % |
| A2_famB | 0.1 | 100 | 39.269909 | 0.039815 | 0.0521 % | 0.0507 % |
| A2_famB | 0.1 | 400 | 78.539818 | 0.037964 | 0.0245 % | 0.0242 % |
| A2_famB | 0.1 | 900 | 117.809723 | 0.036112 | 0.0155 % | 0.0153 % |
| A2_famB | 0.1 | 1600 | 157.079636 | 0.034261 | 0.0110 % | 0.0109 % |
| A2_famB | 0.1 | 2500 | 196.349548 | 0.032430 | 0.0083 % | 0.0083 % |
| A2_famB | 0.3 | 100 | 13.089969 | 0.013270 | 0.0550 % | 0.0507 % |
| A2_famB | 0.3 | 400 | 26.179939 | 0.026545 | 0.0528 % | 0.0507 % |
| A2_famB | 0.3 | 900 | 39.269909 | 0.039815 | 0.0521 % | 0.0507 % |
| A2_famB | 0.3 | 1600 | 52.359879 | 0.011424 | 0.0111 % | 0.0109 % |
| A2_famB | 0.3 | 2500 | 65.449844 | 0.024689 | 0.0192 % | 0.0189 % |
| A2_famB | 0.5 | 100 | 7.853982 | 0.041298 | 0.3024 % | 0.2636 % |
| A2_famB | 0.5 | 400 | 15.707963 | 0.040927 | 0.1394 % | 0.1304 % |
| A2_famB | 0.5 | 900 | 23.561945 | 0.040558 | 0.0900 % | 0.0861 % |
| A2_famB | 0.5 | 1600 | 31.415928 | 0.040192 | 0.0661 % | 0.0640 % |
| A2_famB | 0.6 | 100 | 6.544985 | 0.006636 | 0.0601 % | 0.0507 % |
| A2_famB | 0.6 | 400 | 13.089969 | 0.013270 | 0.0550 % | 0.0507 % |
| A2_famB | 0.6 | 900 | 19.634954 | 0.019908 | 0.0535 % | 0.0507 % |
| A2_famB | 0.6 | 1600 | 26.179939 | 0.026545 | 0.0528 % | 0.0507 % |
| A2_famB | 0.65 | 100 | 6.041524 | 0.041382 | 0.4125 % | 0.3437 % |
| A2_famB | 0.65 | 400 | 12.083049 | 0.041097 | 0.1858 % | 0.1704 % |
| A2_famB | 0.65 | 900 | 18.124573 | 0.040812 | 0.1193 % | 0.1127 % |
| A2_famB | 0.65 | 1600 | 24.166098 | 0.040527 | 0.0876 % | 0.0839 % |

largest A2 c_s shift: 0.4125 %. The zoom and N100-vs-A2 figures therefore pair a corrected A1 v2 curve with uncorrected A2 points; their titles say so.

#### 14.2.1 Regeneration record and draft audit (2026-10-14)

**Dated copies first (copy, never move).** These were made before anything was regenerated, with `cp -n` and checked with `cmp`. There are twelve, in `paper1_speedofsound/experiments/final/`, each named `<name>_pre_boxtrunc_20261014.<ext>`:
- `260919_A1v2_final_cs_vs_eta.csv`;
- `260919_cs_vs_eta`, `260919_cs_vs_eta_lowdensity_zoom`, `260919_cs_vs_eta_N100_vs_A2` and `261002_p1_melting_region`, each as `.png` and `.pdf`;
- `260922_roman2002_remapped_vs_KR` as `.csv`, `.png` and `.pdf`.

**Regenerated, each by the script that originally made it (DATA).**
- **Canonical table** `260919_A1v2_final_cs_vs_eta.csv`, by `paper1_populate_cs_err_20261002.py`.
  - It first recomputes the uncorrected table from the raw traces. It refuses to write unless that recomputation reproduces the dated copy to $5\times10^{-6}$ in every cell; all 35 passed.
  - It then applies $L_{\rm eff,true}$ and $\eta_{\rm true}$. The `eta` column now holds $\eta_{\rm true}$, and three columns are new: `eta_rec`, `delta_sigma` and `L_eff_true`.
  - `KR` and `dev_KR_pct` are evaluated at $\eta_{\rm true}$ for the same 24 cells as before ($\eta_{\rm rec} \le 0.69$).
- **Main figure, low-density zoom and the N = 100 vs A2 overlay**, by `paper1_canonical_20260919.py`, which calls `lowdensity_zoom_20260917.py` and `overlay_N100_vs_A2_20260915.py`. The A2 inputs it would rebuild (260916/260917) are not on disk, so the A2 tables were left byte-identical to HEAD.
- **Román comparison** by `roman2002_remapped_20260922.py`; **melting region** by `paper1_melting_figure_20261002.py`.
- **Titles and style.** Every regenerated title carries "corrected for box truncation (methods §14)". The two figures that show A2 points add "A2 points not yet corrected". The style is unchanged: data blue, KR red, error bars, no Liu 2021.

**Not regenerated, and why.**
- **The four `261001_p1_*` figures** (`paper1_figures_20261001.py`) are drawn from raw traces, not from this table.
  - The estimator floor and the mass-ladder residuals are relative, per-density quantities. They are invariant when $x$ is rescaled by one factor per density.
  - The worked ladder line and the slow mode are drawn at the recorded geometry of the $\eta_{\rm rec} = 0.1122$ cell. That cell's $c_s$ moves from 1.81155 to 1.81044 (§ 14.1 table), and the draft labels the cell by its recorded $\eta$.
  - **OPEN:** redraw the ladder line at $L_{\rm eff,true}$ if its legend is to match the table.
- **The tracked `writeup/paper1_draft.pdf` was not rebuilt.** The edited tex compiles cleanly into the scratchpad (two passes, no undefined references).

**Downstream readers of the canonical CSV** (`grep -rl --include='*.py' 260919_A1v2_final_cs_vs_eta hspist3`):
- **Unchanged inputs.** `paper1_confinement_prereg_20261012.py` and `roman2002_tableII_20261012.py` read the $\pi/8$ row, where $\delta = 0$, so they are unchanged. `paper1_modegate_20261013.py` reads only `L0`, which is unchanged, and computes the truncation itself.
- **`paper1_nofuse_check_20261010.py`** selects the row with `if abs(float(r["eta"]) - ETA) < 1e-6: REF = r`, where `ETA = 0.112200`. The `eta` column now holds 0.112267, so a rerun would find no row and stop.
  - Its committed result was computed on the pre-correction table and stands as a historical record.
  - Pointing it at the dated copy or at `eta_rec` is a one-line change, left for a go.
- **`paper2_level3_v6_20260924.py:163–172`** interpolates this table's deviation to $\eta = 0.1$.
  - INFERENCE: a rerun would give a slightly smaller box-stiffening ratio, because the deviation at the two bracketing cells fell by 0.05 and 0.08 percentage points (table below).
  - Its committed numbers stand as computed on the pre-correction table.

**Erratum to § 14.1 (DATA, from its own printed table).** Two sentences there are not exact: "Every change is negative" and "the deviation from KR shrinks everywhere".
- At $\eta_{\rm rec} = 0.690$ the change in $D$ is $+0.13$. There $c_s^{\rm KR}(\eta_{\rm true}) = 11.52719$ is below $c_s^{\rm KR}(\eta_{\rm rec}) = 11.55102$, so the higher $\eta$ lowers KR.
- At $\eta_{\rm rec} = 0.610$, $|D|$ grows from 0.58 to 0.59.
- The verdict is unaffected: all six cells that changed by more than 0.5 moved toward KR.

**Two draft numbers that the old table did not reproduce.**
- **Zoom caption (l. 247), "$0.5$–$1.5\,\%$".** The old table gives 0.41 to 1.50. The text was replaced by the new table's 0.40 to 1.40, printed as "$0.4$–$1.4\,\%$".
- **Ladder caption (l. 214), "$1.0\,\%$".** The old table gives 1.10; the corrected table gives 1.02, which matches the text, so it was not edited.

**Draft audit (A3), printed by `python3 hspist3/validation/paper1_draft_audit_20261014.py` before any edit** (verbatim). $D$ here comes from the table's 5-decimal values and can differ from § 14.2 by 0.01.

### Paper 1 draft audit: every number from the canonical A1 v2 table, old -> new (printed BEFORE editing)

| tex line(s) | old text | new text | definition | old table gives | old text reproduced? | new table gives | action |
|---|---|---|---|---|---|---|---|
| 33, 242, 303 | `$+1.04\,\%$` | `$+1.01\,\%$` | mean dev from KR, eta <= 0.4 (Roman script def.; 14 cells) | +1.04 | yes | +1.01 | EDIT |
| 243 | `never worse than $1.68\,\%$` | `never worse than $1.68\,\%$` | max dev, eta <= 0.4 (the pi/8 cell, delta = 0) | 1.68 | yes | 1.68 | UNCHANGED |
| 247 | `sits $0.5$--$1.5\,\%$ above` | `sits $0.4$--$1.4\,\%$ above` | dev range of the N = 100 points, eta <= 0.15 (zoom XMAX) | 0.41 to 1.50 | **NO** | 0.40 to 1.40 | EDIT |
| 214 | `sits $1.0\,\%$ above Kolafa` | `sits $1.0\,\%$ above Kolafa` | dev at the worked-ladder cell, eta_rec = 0.1122 | 1.10 | **NO** | 1.02 | UNCHANGED |
| 204 | `That scatter, $0.27\,\%$` | `That scatter, $0.27\,\%$` | median of c_s_scatter_mass / c_s over eta_rec <= 0.69 | 0.268 | yes | 0.268 | UNCHANGED |
| 310 | `ours sits $+2.6\,\%$` | `ours sits $+2.6\,\%$` | dev at eta = 0.5236 (L_0 = 7.5, delta = 0) | +2.63 | yes | +2.63 | UNCHANGED |
| 343 | `local maximum & $\eta = 0.700$` | `local maximum & $\eta = 0.701$` | eta of the local c_s maximum (melting script) | 0.7000 | yes | 0.7007 | EDIT |
| 344 | `local minimum & $\eta = 0.715$` | `local minimum & $\eta = 0.717$` | eta of the local c_s minimum (melting script) | 0.7150 | yes | 0.7167 | EDIT |
| 349 | `$c_s(0.700) - c_s(0.715) = 2.80$` | `$c_s(0.701) - c_s(0.717) = 2.83$` | depth of the dip (melting script) | 2.804 | yes | 2.827 | EDIT |
| 349 | `$\mathbf{9.5\sigma}$` | `$\mathbf{9.6\sigma}$` | dip / plotted (scaled) errors in quadrature | 9.46 | yes | 9.56 | EDIT |
| 350 | `$18.2\sigma$ on the propagated` | `$18.4\sigma$ on the propagated` | dip / propagated errors in quadrature | 18.22 | yes | 18.43 | EDIT |
| 354 | `$\chi^2_{\mathrm{red}} = 1.8$--$14$` | `$\chi^2_{\mathrm{red}} = 1.8$--$14$` | chi2_red range, 0.695 <= eta <= 0.720 | 1.8-14 (6 cells) | yes | 1.8-14 (5 cells) | UNCHANGED |
| 355 | `$0.40$ inside` | `$0.45$ inside` | mean c_s_scatter_mass, 0.695 <= eta <= 0.720 | 0.404 | yes | 0.451 | EDIT |
| 355 | `against $0.08$ outside` | `against $0.09$ outside` | mean c_s_scatter_mass, rest of the melting-figure range 0.66-0.765 | 0.079 | yes | 0.088 | EDIT |

'a factor five' (l. [355]): inside/outside = 5.13 before, 5.15 after -> unchanged.
'both within our grid spacing of 0.005': |0.7007 - 0.702| = 0.0013, |0.7167 - 0.714| = 0.0027 -> still true.
'0.73' (the nofuse sigma check) in the draft: 0 occurrences -> not quoted, nothing to audit.
abstract 'eta = 0.0065 to 0.76': corrected range 0.0065 to 0.7601 -> unchanged.
l. 101 thickness factor (t/2)/(L_0 - 2r) at eta = 0.65 with L_0,true: 0.498 % (text 0.50 %); at eta = 0.0065 delta = 0 -> unchanged.
Cell labels NOT edited: 'eta = 0.1122' (ladder-line and slow-mode captions, ll. 211, 163) names the cell by its recorded eta, as do the two 261001 figures drawn from its raw traces; its corrected eta is 0.112267. '24 densities with eta <= 0.69' (l. 222) is the canonical KR cut, applied to the recorded eta; the 24th cell's corrected eta is 0.6905.

### Per-cell deviations before and after (KR = the table's own KR column; D in units of c_s_err_scaled)

| eta_rec | eta_true | c_s - KR before | c_s - KR after | dev before [%] | dev after [%] | D before | D after | change in D | abs(D) |
|---|---|---|---|---|---|---|---|---|---|
| 0.006545 | 0.006545 | +0.01007 | +0.01007 | +0.703 | +0.703 | +4.24 | +4.24 | +0.00 | same |
| 0.009817 | 0.009817 | +0.01150 | +0.01150 | +0.797 | +0.797 | +2.65 | +2.65 | +0.00 | same |
| 0.013090 | 0.013090 | +0.00808 | +0.00808 | +0.556 | +0.556 | +3.27 | +3.27 | +0.00 | same |
| 0.019635 | 0.019637 | +0.00606 | +0.00591 | +0.412 | +0.402 | +1.37 | +1.34 | -0.03 | smaller |
| 0.026180 | 0.026184 | +0.01258 | +0.01237 | +0.844 | +0.830 | +2.84 | +2.80 | -0.05 | smaller |
| 0.039270 | 0.039278 | +0.01076 | +0.01041 | +0.702 | +0.680 | +2.57 | +2.48 | -0.08 | smaller |
| 0.052360 | 0.052374 | +0.01301 | +0.01252 | +0.826 | +0.795 | +5.07 | +4.88 | -0.19 | smaller |
| 0.078540 | 0.078573 | +0.01412 | +0.01329 | +0.848 | +0.798 | +2.68 | +2.53 | -0.16 | smaller |
| 0.112200 | 0.112267 | +0.01966 | +0.01828 | +1.097 | +1.020 | +3.71 | +3.45 | -0.26 | smaller |
| 0.130900 | 0.130991 | +0.02798 | +0.02623 | +1.497 | +1.403 | +5.92 | +5.55 | -0.37 | smaller |
| 0.157080 | 0.157210 | +0.03121 | +0.02887 | +1.572 | +1.454 | +11.25 | +10.42 | -0.83 | smaller |
| 0.196350 | 0.196350 | +0.03115 | +0.03115 | +1.429 | +1.429 | +5.61 | +5.61 | +0.00 | same |
| 0.261799 | 0.261799 | +0.04044 | +0.04044 | +1.572 | +1.572 | +5.28 | +5.28 | +0.00 | same |
| 0.392699 | 0.392699 | +0.06276 | +0.06276 | +1.675 | +1.675 | +5.61 | +5.61 | +0.00 | same |
| 0.523599 | 0.523599 | +0.15602 | +0.15602 | +2.630 | +2.630 | +9.11 | +9.11 | +0.00 | same |
| 0.549999 | 0.551157 | +0.13803 | +0.08957 | +2.091 | +1.350 | +3.85 | +2.50 | -1.35 | smaller |
| 0.569996 | 0.571198 | +0.02946 | -0.02605 | +0.410 | -0.361 | +0.57 | -0.50 | -1.07 | smaller |
| 0.590001 | 0.590895 | +0.07042 | +0.02434 | +0.897 | +0.309 | +0.77 | +0.27 | -0.50 | smaller |
| 0.609999 | 0.610018 | -0.03782 | -0.03891 | -0.439 | -0.452 | -0.58 | -0.59 | -0.02 | LARGER |
| 0.630002 | 0.630420 | +0.02981 | +0.00307 | +0.314 | +0.032 | +0.23 | +0.02 | -0.21 | smaller |
| 0.650003 | 0.652234 | +0.59189 | +0.43657 | +5.670 | +4.139 | +11.17 | +8.27 | -2.90 | smaller |
| 0.669998 | 0.670803 | +1.62224 | +1.57171 | +14.264 | +13.781 | +41.31 | +40.08 | -1.23 | smaller |
| 0.679998 | 0.680489 | +2.84436 | +2.82445 | +24.371 | +24.185 | +57.91 | +57.56 | -0.35 | smaller |
| 0.689999 | 0.690460 | +4.96829 | +4.97867 | +43.012 | +43.191 | +45.35 | +45.48 | +0.13 | LARGER |
| 0.695006 | 0.695556 | n/a | n/a | n/a | n/a | n/a | n/a | n/a | no KR in table |
| 0.699998 | 0.700727 | n/a | n/a | n/a | n/a | n/a | n/a | n/a | no KR in table |
| 0.705000 | 0.705976 | n/a | n/a | n/a | n/a | n/a | n/a | n/a | no KR in table |
| 0.709997 | 0.711304 | n/a | n/a | n/a | n/a | n/a | n/a | n/a | no KR in table |
| 0.714999 | 0.716713 | n/a | n/a | n/a | n/a | n/a | n/a | n/a | no KR in table |
| 0.719994 | 0.722205 | n/a | n/a | n/a | n/a | n/a | n/a | n/a | no KR in table |
| 0.725005 | 0.727782 | n/a | n/a | n/a | n/a | n/a | n/a | n/a | no KR in table |
| 0.730005 | 0.730603 | n/a | n/a | n/a | n/a | n/a | n/a | n/a | no KR in table |
| 0.740006 | 0.742108 | n/a | n/a | n/a | n/a | n/a | n/a | n/a | no KR in table |
| 0.749998 | 0.750978 | n/a | n/a | n/a | n/a | n/a | n/a | n/a | no KR in table |
| 0.759999 | 0.760063 | n/a | n/a | n/a | n/a | n/a | n/a | n/a | no KR in table |

cells with |change in D| > 0.5: 6 at eta_rec = 0.157, 0.550, 0.570, 0.590, 0.650, 0.670; all toward KR: True
cells where |D| grows (by > 0.005): eta_rec = 0.610: -0.58 -> -0.59; eta_rec = 0.690: +45.35 -> +45.48

### Regeneration checks

rows changed: 28 of 35; unchanged (delta = 0): L_0 = 600, 400, 300, 20.0000, 15.0000, 10.0000, 7.5000
identity c_s,new = c_s,old * L_eff,true/L_eff,rec: max |difference| 9.0e-06 (table prints 5 decimals)
identity eta_new = eta_old * L_0/(L_0 - delta/2): max |difference| 5.1e-07 (table prints 6 decimals)

### Numbers for the Methods paragraph

delta = 0 at 7 of 35 densities; c_s lowered by at most 0.09 % for eta <= 0.16 (at eta_rec = 0.1571) and by at most 0.47 % overall (at eta_rec = 0.7250); eta raised by at most 0.38 % (at eta_rec = 0.7250); 6 densities changed D by more than 0.5, all toward KR: True; Center_X check quoted from methods sec. 14.1: 5.5\times10^{-6} sigma.

**Applied** with `--apply`, which printed (verbatim):

    edited (3x): $+1.04\,\%$ -> $+1.01\,\%$
    edited (1x) -- old text was NOT reproduced by the old table; new text is the new table value: sits $0.5$--$1.5\,\%$ above -> sits $0.4$--$1.4\,\%$ above
    edited (1x): local maximum & $\eta = 0.700$ -> local maximum & $\eta = 0.701$
    edited (1x): local minimum & $\eta = 0.715$ -> local minimum & $\eta = 0.717$
    edited (1x): $c_s(0.700) - c_s(0.715) = 2.80$ -> $c_s(0.701) - c_s(0.717) = 2.83$
    edited (1x): $\mathbf{9.5\sigma}$ -> $\mathbf{9.6\sigma}$
    edited (1x): $18.2\sigma$ on the propagated -> $18.4\sigma$ on the propagated
    edited (1x): $0.40$ inside -> $0.45$ inside
    edited (1x): against $0.08$ outside -> against $0.09$ outside

    applied: edits + Methods paragraph written to paper1_draft.tex

**The Methods paragraph** ("Integer-pixel box width", five sentences, placed after "Effective length" in § II) takes its numbers from the last block above. The one exception is the Center_X agreement, $5.5\times10^{-6}\,\sigma$, which is quoted from § 14.1. The paragraph says "for every $N = 100$ density" rather than "throughout", and it states that the larger systems of the finite-size section are not yet corrected.

**§ 14 addendum: summary rotation (2026-10-14; no code change).** From source 615561c on (first binary 5190846-dirty, now e823187), the energy-transfer summary header ends in `,build_cflags,box_width_sigma` (`00ALLINONE.c:17283`), and `ensure_csv_header_schema` (called at 17285) renames any shared `summary.csv` whose first line differs to `<path>.legacy_<unix time>` before writing a fresh file, so the five Level 3-era runners in `hspist3/experiments_energy_transfer/_run_scripts/` that append to `--energy-transfer-summary="$d/summary.csv"` (`level3_master_preload_20260921.sh`, `level3_master_20260920.sh`, `level3_v6_longrecords.sh`, `level3_FofL_20260925.sh`, `level4_equilibrium_KOAlength_20261002.sh`) must not be rerun into their existing directories with this or any newer build, and any rerun gets a new output directory. The rename, `00ALLINONE.c:2903–2904`:

    snprintf(backup, sizeof(backup), "%s.legacy_%ld", path, (long)now);
    if (rename(path, backup) != 0) {

### 14.3 A2 size-ladder correction — pre-registration (2026-10-02, machine date; written and committed before anything below is computed)

**What is corrected.** § 14.2 (OPEN table) showed that the A2 ladder uses non-grid $L_0$ too, with $\delta$ up to $0.0414\,\sigma$ and a $c_s$ shift up to $0.41\,\%$ (η = 0.65, N = 100). The geometry is the same as in § 14: the divider starts $L_0$ from the left wall and oscillates about the centre of the truncated box. So the correction is the § 14 one, per (η, N) cell:
$$L_{\rm eff,true} = L_{\rm eff} - \frac{\delta}{2},\qquad \eta_{\rm true} = \eta\,\frac{L_0}{L_0 - \delta/2},\qquad \delta = \frac{2L_0\cdot 24 - \lfloor 2L_0\cdot 24\rfloor}{24},$$
with $\delta$ computed by the binary's float32 expression (`SIM_WIDTH = (int)(2*L0_UNITS*PIXELS_PER_SIGMA)`, `00ALLINONE.c:323`), as in § 14.2.

**The ladder cells (DATA, from the runs' own command files).** $L_0$ is the value recorded in the per-mass tables (`260919_A2_cs_per_mass.csv` and `_famB`; listed with its $\delta$ per cell in the § 14.2 OPEN table). $H = 10\sqrt{N/100}$:

| N | 100 | 400 | 900 | 1600 | 2500 |
|---|---|---|---|---|---|
| H | 10 | 20 | 30 | 40 | 50 |

These are read from `--height=` in the `00_COMMAND.md` files of `A2_dilute_20260916`, `A2_dilute50_20260917`, `A2_topup_20260912`, `A2_alpha2_20260912` and `A2_long200_20260915`. For cells whose directories hold no `00_COMMAND.md`, the script reads $H$ back from the trace's recorded η, $H = N\pi r^2/(2L_0\eta_{\rm rec})$, as in § 14.1.
- **Height cast:** `SIM_HEIGHT = (int)(HEIGHT_UNITS * PIXELS_PER_SIGMA);` (`00ALLINONE.c:324`). $H\cdot 24 = 240, 480, 720, 960, 1200$ are integers, so the height is never truncated. The script checks this per cell.
- **Off-grid cells:** every A2 cell with $\delta > 0$ in § 14.2. In that table, all 60 rows of both families are off-grid (smallest $\delta = 0.0047\,\sigma$).

**What is recomputed, and from what.** No trajectory is re-analysed. The per-mass frequencies $\bar\nu_M$ are data. The correction changes only $x_M = K(\alpha)/(2\pi L_{\rm eff})$, by one factor per (η, N) cell, plus the η at which KR is evaluated.
1. **Estimator gate first.** From `260919_A2_cs_per_mass.csv` ($\bar\nu_M$, $L_0$), the published fit `260917_A2_cs_vs_N_extrapolation.csv` (printed by `analyze_A2_X2p5_20260914.py`) must be reproduced to its printed precision, with that script's own definitions:
   - $x_M$ with $L_{\rm eff} = L_0 - 2r$;
   - per-cell $c_s$ = the unweighted through-origin slope;
   - weight = SD of $\bar\nu_M/x_M$;
   - `sos.weighted_linreg` on $1/\sqrt N$.

   If it is not reproduced, the analysis stops and the verdict is VOID.
2. **Three columns are reported:**
   - **(P)** as published, $L_{\rm eff} = L_0 - 2r$;
   - **(B) before** = the current canonical geometry, $L_{\rm eff} = L_0 - 2r - t/2$, as in the A1 table since 2026-09-18;
   - **(A) after** = $L_0 - 2r - t/2 - \delta/2$ and $\eta_{\rm true}$.

   **The verdict compares B with A only**, i.e. the box truncation alone, as in § 14. P → B is the divider-thickness factor, which the draft's finite-size numbers do not yet carry. That is reported separately and does not enter the verdict.

**Per-point test.**
- For every ladder point, $D = (c_s - c_s^{\rm KR}(\eta))/\sigma$. Here σ is the error bar the draft plots for that point: `slope_with_errors` with $\sigma_\nu = {\rm sd}/\sqrt n$, inflated by $\sqrt{\chi^2_{\rm red}}$, as in `overlay_N100_vs_A2_20260915.py:42` and `lowdensity_zoom_20260917.py:44`. KR is taken at the point's own η (B: η; A: $\eta_{\rm true}$).
- Both per-mass tables are tested.
- Points above η = 0.69 do not exist in A2.

**Fit test.** The draft defines the finite-size fit as $c_s(N) = c_\infty + b/\sqrt N$, weighted by the mass scatter (`analyze_A2_X2p5_20260914.py:120`). The exponent is fixed at 1/2; the draft fits no exponent, so the two fitted parameters are the intercept $c_\infty$ and the coefficient $b$.
- **Why the fit uses deviations.** After the correction, the points of one ladder no longer share one η, because $\eta_{\rm true}$ depends on N. The A fit is therefore made on the deviations $\Delta_N = c_s - c_s^{\rm KR}(\eta_{\rm true,N})$, with the same weights, and reported as $c_\infty = c_s^{\rm KR}(\eta) + \Delta_\infty$.
- In B all points share η, so this is identical to the published form.
- The test is on the draft's source table, `260919_A2_cs_per_mass.csv` (the 260917 data).

**Verdict rule (fixed now):**
- **REGENERATE** if any ladder point's $D$ changes by more than 0.5 (of that point's own σ) between B and A, in either table; **or** the fit's $c_\infty$ or $b$ changes by more than its own 1σ (the B error) at any η.
- **KEEP** otherwise.

**If REGENERATE:**
- dated copies `_pre_boxtrunc_261002` (`cp -n`, `cmp`) of both A2 per-mass tables, `260917_A2_cs_vs_N_extrapolation.csv`, and the two figures that show A2 points;
- the tables are regenerated with $c_s$-like columns at $L_{\rm eff,true}$, plus `eta_true` and `delta_sigma` columns. The figures are regenerated by their original scripts with "corrected for box truncation (methods §14)", and the "A2 points not yet corrected" annotation goes.
- The extrapolation table is regenerated by `paper1_A2_boxtrunc_261002.py` with the original script's fit (same function, same weights) in the A geometry.

**If KEEP:** nothing is regenerated, and the draft says in one sentence that the correction changes no A2 result by more than the stated amount.

#### 14.3.1 Results (2026-10-02; computed once, after the § 14.3 commit 0ddefa9)

**Verdict, by the pre-registered rule: REGENERATE (DATA).**
- **Per-point test.** Four points change $D$ by more than 0.5: η = 0.50 at N = 100 (−0.53) and N = 400 (−0.73), in both per-mass tables.
- **Fit test.** The coefficient $b$ moves by more than its own 1σ at η = 0.50 (−1.10σ) and η = 0.65 (−1.03σ).
- **What does not move.** No intercept $c_\infty$ moves by more than 0.34σ. At η = 0.02 and 0.05, the two values the draft quotes, they move by 0.00σ and 0.02σ.
- **Direction.** Every point's $D$ decreases. At η ≥ 0.6 this pushes some negative $D$ further from zero, e.g. η = 0.65, N = 1600: −1.61 → −2.03.

**Why $b$ moves while $c_\infty$ stays (INFERENCE).** The correction is largest in the smallest box: at N = 100, $(\delta/2)/L_{\rm eff}$ reaches 0.41 %, and $\eta_{\rm true}$ raises KR most where KR is steep. Both act on the small-N end of each ladder ($1/\sqrt N = 0.1$), which tilts the line. The intercept is held by the large-N points, which barely move.

**Estimator gate: PASS.** All seven published rows of `260917_A2_cs_vs_N_extrapolation.csv` are reproduced to their printed digits ($c_\infty$, its error, $b$, its error, $\chi^2$). That table was computed with $L_{\rm eff} = L_0 - 2r$ (geometry P). The divider-thickness step P → B, which the draft's finite-size numbers never carried, moves $c_\infty$ by at most 0.05σ. It does not enter the verdict.

**Deviation from the pre-registration (disclosed).** § 14.3 named five campaign roots as the source of `--height=`, but the N = 100 and 400 ladder cells live in a sixth, `famB_20260911`. The script searches it too; it does not search the `_ABANDONED_seedpad` copy. All 30 (η, N) cells got H from a command file, none by read-back, and H·24 is an integer in every cell.

**Printed by `python3 hspist3/validation/paper1_A2_boxtrunc_261002.py --write`** (verbatim; the analysis part is identical to the first, read-only run apart from the inputs line):

inputs (uncorrected): 260919_A2_cs_per_mass_pre_boxtrunc_261002.csv, 260919_A2_cs_per_mass_famB_pre_boxtrunc_261002.csv, 260917_A2_cs_vs_N_extrapolation_pre_boxtrunc_261002.csv

### Estimator gate: 260917_A2_cs_vs_N_extrapolation.csv reproduced from 260919_A2_cs_per_mass.csv (geometry P, analyze_A2_X2p5_20260914.py definitions)

| eta | N values | c_inf (pub) | c_inf (here) | ± (pub/here) | b (pub/here) | ± (pub/here) | chi2 (pub/here) | reproduced |
|---|---|---|---|---|---|---|---|---|
| 0.02 | 100 400 900 1600 | 1.46140 | 1.46140 | 0.01023 / 0.01023 | 0.18671 / 0.18671 | 0.14755 / 0.14755 | 2.075 / 2.075 | yes |
| 0.05 | 100 400 900 1600 | 1.56988 | 1.56988 | 0.00291 / 0.00291 | 0.08034 / 0.08034 | 0.03084 / 0.03084 | 0.446 / 0.446 | yes |
| 0.10 | 100 400 900 1600 2500 | 1.75361 | 1.75361 | 0.00493 / 0.00493 | -0.01399 / -0.01399 | 0.13503 / 0.13503 | 2.143 / 2.143 | yes |
| 0.30 | 100 400 900 1600 2500 | 2.86354 | 2.86354 | 0.00803 / 0.00803 | 0.46011 / 0.46011 | 0.18528 / 0.18528 | 1.999 / 1.999 | yes |
| 0.50 | 100 400 900 1600 | 5.43030 | 5.43030 | 0.02368 / 0.02368 | 0.99980 / 0.99980 | 0.40674 / 0.40674 | 0.385 / 0.385 | yes |
| 0.60 | 100 400 900 1600 | 8.20534 | 8.20534 | 0.04708 / 0.04708 | 1.36046 / 1.36046 | 1.02568 / 1.02568 | 0.673 / 0.673 | yes |
| 0.65 | 100 400 900 1600 | 10.16323 | 10.16323 | 0.06872 / 0.06872 | 6.43163 / 6.43163 | 1.53622 / 1.53622 | 0.598 / 0.598 | yes |

estimator gate: PASS -- every published row reproduced to its printed digits

### Per ladder point (sigma = plotted error; KR at the point's own eta)

| table | eta_rec | N | L_0 | H (source) | H*24 | delta | eta_true | L_eff,rec (B) | L_eff,true (A) | c_s before ± σ | c_s after ± σ | KR(eta_rec) | KR(eta_true) | D before | D after | change |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A2 | 0.02 | 100 | 196.349548 | 10.00000 (command) | 240.000 | 0.032430 | 0.020002 | 195.324548 | 195.308333 | 1.47943 ± 0.00543 | 1.47930 ± 0.00543 | 1.47247 | 1.47248 | +1.28 | +1.26 | -0.02 |
| A2 | 0.02 | 400 | 392.699097 | 20.00000 (command) | 480.000 | 0.023193 | 0.020001 | 391.674097 | 391.662500 | 1.48936 ± 0.01474 | 1.48931 ± 0.01474 | 1.47247 | 1.47247 | +1.15 | +1.14 | -0.00 |
| A2 | 0.02 | 900 | 589.048645 | 30.00000 (command) | 720.000 | 0.013997 | 0.020000 | 588.023645 | 588.016646 | 1.46132 ± 0.00779 | 1.46131 ± 0.00779 | 1.47247 | 1.47247 | -1.43 | -1.43 | -0.00 |
| A2 | 0.02 | 1600 | 785.398193 | 40.00000 (command) | 960.000 | 0.004720 | 0.020000 | 784.373193 | 784.370833 | 1.47384 ± 0.00765 | 1.47384 ± 0.00765 | 1.47247 | 1.47247 | +0.18 | +0.18 | -0.00 |
| A2 | 0.05 | 100 | 78.539818 | 10.00000 (command) | 240.000 | 0.037964 | 0.050012 | 77.514818 | 77.495836 | 1.57746 ± 0.00563 | 1.57707 ± 0.00563 | 1.56662 | 1.56666 | +1.93 | +1.85 | -0.08 |
| A2 | 0.05 | 400 | 157.079636 | 20.00000 (command) | 480.000 | 0.034261 | 0.050005 | 156.054636 | 156.037505 | 1.57133 ± 0.01034 | 1.57116 ± 0.01034 | 1.56662 | 1.56664 | +0.46 | +0.44 | -0.02 |
| A2 | 0.05 | 900 | 235.619446 | 30.00000 (command) | 720.000 | 0.030558 | 0.050003 | 234.594446 | 234.579167 | 1.57515 ± 0.00475 | 1.57505 ± 0.00475 | 1.56662 | 1.56663 | +1.80 | +1.77 | -0.02 |
| A2 | 0.05 | 1600 | 314.159271 | 40.00000 (command) | 960.000 | 0.026855 | 0.050002 | 313.134271 | 313.120843 | 1.57358 ± 0.00486 | 1.57352 ± 0.00486 | 1.56662 | 1.56663 | +1.43 | +1.42 | -0.02 |
| A2 | 0.10 | 100 | 39.269909 | 10.00000 (command) | 240.000 | 0.039815 | 0.100051 | 38.244909 | 38.225001 | 1.74530 ± 0.01814 | 1.74439 ± 0.01813 | 1.74414 | 1.74434 | +0.06 | +0.00 | -0.06 |
| A2 | 0.10 | 400 | 78.539818 | 20.00000 (command) | 480.000 | 0.037964 | 0.100024 | 77.514818 | 77.495836 | 1.75859 ± 0.01518 | 1.75816 ± 0.01518 | 1.74414 | 1.74424 | +0.95 | +0.92 | -0.03 |
| A2 | 0.10 | 900 | 117.809723 | 30.00000 (command) | 720.000 | 0.036112 | 0.100015 | 116.784723 | 116.766667 | 1.75499 ± 0.00356 | 1.75472 ± 0.00356 | 1.74414 | 1.74420 | +3.05 | +2.96 | -0.09 |
| A2 | 0.10 | 1600 | 157.079636 | 40.00000 (command) | 960.000 | 0.034261 | 0.100011 | 156.054636 | 156.037505 | 1.74845 ± 0.00293 | 1.74826 ± 0.00293 | 1.74414 | 1.74418 | +1.47 | +1.39 | -0.08 |
| A2 | 0.10 | 2500 | 196.349548 | 50.00000 (command) | 1200.000 | 0.032430 | 0.100008 | 195.324548 | 195.308333 | 1.74812 ± 0.00562 | 1.74797 ± 0.00562 | 1.74414 | 1.74417 | +0.71 | +0.68 | -0.03 |
| A2 | 0.30 | 100 | 13.089969 | 10.00000 (command) | 240.000 | 0.013270 | 0.300152 | 12.064969 | 12.058334 | 2.89623 ± 0.01675 | 2.89464 ± 0.01674 | 2.85150 | 2.85270 | +2.67 | +2.51 | -0.17 |
| A2 | 0.30 | 400 | 26.179939 | 20.00000 (command) | 480.000 | 0.026545 | 0.300152 | 25.154939 | 25.141666 | 2.89132 ± 0.01285 | 2.88980 ± 0.01284 | 2.85150 | 2.85270 | +3.10 | +2.89 | -0.21 |
| A2 | 0.30 | 900 | 39.269909 | 30.00000 (command) | 720.000 | 0.039815 | 0.300152 | 38.244909 | 38.225001 | 2.87063 ± 0.00880 | 2.86913 ± 0.00880 | 2.85150 | 2.85270 | +2.17 | +1.87 | -0.30 |
| A2 | 0.30 | 1600 | 52.359879 | 40.00000 (command) | 960.000 | 0.011424 | 0.300033 | 51.334879 | 51.329167 | 2.87156 ± 0.00495 | 2.87124 ± 0.00495 | 2.85150 | 2.85176 | +4.05 | +3.93 | -0.12 |
| A2 | 0.30 | 2500 | 65.449844 | 50.00000 (command) | 1200.000 | 0.024689 | 0.300057 | 64.424844 | 64.412500 | 2.88463 ± 0.00953 | 2.88408 ± 0.00953 | 2.85150 | 2.85195 | +3.48 | +3.37 | -0.10 |
| A2 | 0.50 | 100 | 7.853982 | 10.00000 (command) | 240.000 | 0.041298 | 0.501318 | 6.828982 | 6.808333 | 5.50877 ± 0.08249 | 5.49211 ± 0.08224 | 5.41622 | 5.44323 | +1.12 | +0.59 | -0.53 |
| A2 | 0.50 | 400 | 15.707963 | 20.00000 (command) | 480.000 | 0.040927 | 0.500652 | 14.682963 | 14.662500 | 5.47002 ± 0.02848 | 5.46240 ± 0.02844 | 5.41622 | 5.42956 | +1.89 | +1.15 | -0.73 |
| A2 | 0.50 | 900 | 23.561945 | 30.00000 (command) | 720.000 | 0.040558 | 0.500431 | 22.536945 | 22.516666 | 5.50548 ± 0.09177 | 5.50053 ± 0.09169 | 5.41622 | 5.42502 | +0.97 | +0.82 | -0.15 |
| A2 | 0.50 | 1600 | 31.415928 | 40.00000 (command) | 960.000 | 0.040192 | 0.500320 | 30.390928 | 30.370832 | 5.43854 ± 0.06297 | 5.43494 ± 0.06293 | 5.41622 | 5.42276 | +0.35 | +0.19 | -0.16 |
| A2 | 0.60 | 100 | 6.544985 | 10.00000 (command) | 240.000 | 0.006636 | 0.600304 | 5.519985 | 5.516667 | 8.33467 ± 0.26871 | 8.32966 ± 0.26855 | 8.22261 | 8.23420 | +0.42 | +0.36 | -0.06 |
| A2 | 0.60 | 400 | 13.089969 | 20.00000 (command) | 480.000 | 0.013270 | 0.600304 | 12.064969 | 12.058334 | 8.20991 ± 0.07608 | 8.20540 ± 0.07604 | 8.22261 | 8.23420 | -0.17 | -0.38 | -0.21 |
| A2 | 0.60 | 900 | 19.634954 | 30.00000 (command) | 720.000 | 0.019908 | 0.600304 | 18.609954 | 18.600000 | 8.27056 ± 0.04403 | 8.26613 ± 0.04401 | 8.22261 | 8.23420 | +1.09 | +0.73 | -0.36 |
| A2 | 0.60 | 1600 | 26.179939 | 40.00000 (command) | 960.000 | 0.026545 | 0.600304 | 25.154939 | 25.141666 | 8.22885 ± 0.05385 | 8.22450 ± 0.05382 | 8.22261 | 8.23420 | +0.12 | -0.18 | -0.30 |
| A2 | 0.65 | 100 | 6.041524 | 10.00000 (command) | 240.000 | 0.041382 | 0.652234 | 5.016524 | 4.995833 | 10.71564 ± 0.35051 | 10.67145 ± 0.34907 | 10.43875 | 10.54875 | +0.79 | +0.35 | -0.44 |
| A2 | 0.65 | 400 | 12.083049 | 20.00000 (command) | 480.000 | 0.041097 | 0.651107 | 11.058049 | 11.037500 | 10.60812 ± 0.19953 | 10.58841 ± 0.19916 | 10.43875 | 10.49328 | +0.85 | +0.48 | -0.37 |
| A2 | 0.65 | 900 | 18.124573 | 30.00000 (command) | 720.000 | 0.040812 | 0.650733 | 17.099573 | 17.079167 | 10.35085 ± 0.10208 | 10.33849 ± 0.10196 | 10.43875 | 10.47483 | -0.86 | -1.34 | -0.48 |
| A2 | 0.65 | 1600 | 24.166098 | 40.00000 (command) | 960.000 | 0.040527 | 0.650545 | 23.141098 | 23.120834 | 10.30054 ± 0.08595 | 10.29152 ± 0.08587 | 10.43875 | 10.46561 | -1.61 | -2.03 | -0.42 |
| famB | 0.02 | 100 | 196.349548 | 10.00000 (command) | 240.000 | 0.032430 | 0.020002 | 195.324548 | 195.308333 | 1.48161 ± 0.00755 | 1.48148 ± 0.00755 | 1.47247 | 1.47248 | +1.21 | +1.19 | -0.02 |
| famB | 0.02 | 400 | 392.699097 | 20.00000 (command) | 480.000 | 0.023193 | 0.020001 | 391.674097 | 391.662500 | 1.47770 ± 0.00995 | 1.47766 ± 0.00994 | 1.47247 | 1.47247 | +0.53 | +0.52 | -0.00 |
| famB | 0.02 | 900 | 589.048645 | 30.00000 (command) | 720.000 | 0.013997 | 0.020000 | 588.023645 | 588.016646 | 1.47087 ± 0.01237 | 1.47085 ± 0.01237 | 1.47247 | 1.47247 | -0.13 | -0.13 | -0.00 |
| famB | 0.02 | 1600 | 785.398193 | 40.00000 (command) | 960.000 | 0.004720 | 0.020000 | 784.373193 | 784.370833 | 1.48037 ± 0.00850 | 1.48036 ± 0.00850 | 1.47247 | 1.47247 | +0.93 | +0.93 | -0.00 |
| famB | 0.05 | 100 | 78.539818 | 10.00000 (command) | 240.000 | 0.037964 | 0.050012 | 77.514818 | 77.495836 | 1.58926 ± 0.01206 | 1.58887 ± 0.01206 | 1.56662 | 1.56666 | +1.88 | +1.84 | -0.04 |
| famB | 0.05 | 400 | 157.079636 | 20.00000 (command) | 480.000 | 0.034261 | 0.050005 | 156.054636 | 156.037505 | 1.57788 ± 0.00666 | 1.57771 ± 0.00666 | 1.56662 | 1.56664 | +1.69 | +1.66 | -0.03 |
| famB | 0.05 | 900 | 235.619446 | 30.00000 (command) | 720.000 | 0.030558 | 0.050003 | 234.594446 | 234.579167 | 1.56560 ± 0.00486 | 1.56549 ± 0.00486 | 1.56662 | 1.56663 | -0.21 | -0.23 | -0.02 |
| famB | 0.05 | 1600 | 314.159271 | 40.00000 (command) | 960.000 | 0.026855 | 0.050002 | 313.134271 | 313.120843 | 1.57965 ± 0.00600 | 1.57958 ± 0.00600 | 1.56662 | 1.56663 | +2.17 | +2.16 | -0.01 |
| famB | 0.10 | 100 | 39.269909 | 10.00000 (command) | 240.000 | 0.039815 | 0.100051 | 38.244909 | 38.225001 | 1.74530 ± 0.01814 | 1.74439 ± 0.01813 | 1.74414 | 1.74434 | +0.06 | +0.00 | -0.06 |
| famB | 0.10 | 400 | 78.539818 | 20.00000 (command) | 480.000 | 0.037964 | 0.100024 | 77.514818 | 77.495836 | 1.75859 ± 0.01518 | 1.75816 ± 0.01518 | 1.74414 | 1.74424 | +0.95 | +0.92 | -0.03 |
| famB | 0.10 | 900 | 117.809723 | 30.00000 (command) | 720.000 | 0.036112 | 0.100015 | 116.784723 | 116.766667 | 1.75499 ± 0.00356 | 1.75472 ± 0.00356 | 1.74414 | 1.74420 | +3.05 | +2.96 | -0.09 |
| famB | 0.10 | 1600 | 157.079636 | 40.00000 (command) | 960.000 | 0.034261 | 0.100011 | 156.054636 | 156.037505 | 1.74845 ± 0.00293 | 1.74826 ± 0.00293 | 1.74414 | 1.74418 | +1.47 | +1.39 | -0.08 |
| famB | 0.10 | 2500 | 196.349548 | 50.00000 (command) | 1200.000 | 0.032430 | 0.100008 | 195.324548 | 195.308333 | 1.74812 ± 0.00562 | 1.74797 ± 0.00562 | 1.74414 | 1.74417 | +0.71 | +0.68 | -0.03 |
| famB | 0.30 | 100 | 13.089969 | 10.00000 (command) | 240.000 | 0.013270 | 0.300152 | 12.064969 | 12.058334 | 2.89623 ± 0.01675 | 2.89464 ± 0.01674 | 2.85150 | 2.85270 | +2.67 | +2.51 | -0.17 |
| famB | 0.30 | 400 | 26.179939 | 20.00000 (command) | 480.000 | 0.026545 | 0.300152 | 25.154939 | 25.141666 | 2.89132 ± 0.01285 | 2.88980 ± 0.01284 | 2.85150 | 2.85270 | +3.10 | +2.89 | -0.21 |
| famB | 0.30 | 900 | 39.269909 | 30.00000 (command) | 720.000 | 0.039815 | 0.300152 | 38.244909 | 38.225001 | 2.87063 ± 0.00880 | 2.86913 ± 0.00880 | 2.85150 | 2.85270 | +2.17 | +1.87 | -0.30 |
| famB | 0.30 | 1600 | 52.359879 | 40.00000 (command) | 960.000 | 0.011424 | 0.300033 | 51.334879 | 51.329167 | 2.87156 ± 0.00495 | 2.87124 ± 0.00495 | 2.85150 | 2.85176 | +4.05 | +3.93 | -0.12 |
| famB | 0.30 | 2500 | 65.449844 | 50.00000 (command) | 1200.000 | 0.024689 | 0.300057 | 64.424844 | 64.412500 | 2.88463 ± 0.00953 | 2.88408 ± 0.00953 | 2.85150 | 2.85195 | +3.48 | +3.37 | -0.10 |
| famB | 0.50 | 100 | 7.853982 | 10.00000 (command) | 240.000 | 0.041298 | 0.501318 | 6.828982 | 6.808333 | 5.50877 ± 0.08249 | 5.49211 ± 0.08224 | 5.41622 | 5.44323 | +1.12 | +0.59 | -0.53 |
| famB | 0.50 | 400 | 15.707963 | 20.00000 (command) | 480.000 | 0.040927 | 0.500652 | 14.682963 | 14.662500 | 5.47002 ± 0.02848 | 5.46240 ± 0.02844 | 5.41622 | 5.42956 | +1.89 | +1.15 | -0.73 |
| famB | 0.50 | 900 | 23.561945 | 30.00000 (command) | 720.000 | 0.040558 | 0.500431 | 22.536945 | 22.516666 | 5.50548 ± 0.09177 | 5.50053 ± 0.09169 | 5.41622 | 5.42502 | +0.97 | +0.82 | -0.15 |
| famB | 0.50 | 1600 | 31.415928 | 40.00000 (command) | 960.000 | 0.040192 | 0.500320 | 30.390928 | 30.370832 | 5.43854 ± 0.06297 | 5.43494 ± 0.06293 | 5.41622 | 5.42276 | +0.35 | +0.19 | -0.16 |
| famB | 0.60 | 100 | 6.544985 | 10.00000 (command) | 240.000 | 0.006636 | 0.600304 | 5.519985 | 5.516667 | 8.33467 ± 0.26871 | 8.32966 ± 0.26855 | 8.22261 | 8.23420 | +0.42 | +0.36 | -0.06 |
| famB | 0.60 | 400 | 13.089969 | 20.00000 (command) | 480.000 | 0.013270 | 0.600304 | 12.064969 | 12.058334 | 8.20991 ± 0.07608 | 8.20540 ± 0.07604 | 8.22261 | 8.23420 | -0.17 | -0.38 | -0.21 |
| famB | 0.60 | 900 | 19.634954 | 30.00000 (command) | 720.000 | 0.019908 | 0.600304 | 18.609954 | 18.600000 | 8.27056 ± 0.04403 | 8.26613 ± 0.04401 | 8.22261 | 8.23420 | +1.09 | +0.73 | -0.36 |
| famB | 0.60 | 1600 | 26.179939 | 40.00000 (command) | 960.000 | 0.026545 | 0.600304 | 25.154939 | 25.141666 | 8.22885 ± 0.05385 | 8.22450 ± 0.05382 | 8.22261 | 8.23420 | +0.12 | -0.18 | -0.30 |
| famB | 0.65 | 100 | 6.041524 | 10.00000 (command) | 240.000 | 0.041382 | 0.652234 | 5.016524 | 4.995833 | 10.71564 ± 0.35051 | 10.67145 ± 0.34907 | 10.43875 | 10.54875 | +0.79 | +0.35 | -0.44 |
| famB | 0.65 | 400 | 12.083049 | 20.00000 (command) | 480.000 | 0.041097 | 0.651107 | 11.058049 | 11.037500 | 10.60812 ± 0.19953 | 10.58841 ± 0.19916 | 10.43875 | 10.49328 | +0.85 | +0.48 | -0.37 |
| famB | 0.65 | 900 | 18.124573 | 30.00000 (command) | 720.000 | 0.040812 | 0.650733 | 17.099573 | 17.079167 | 10.35085 ± 0.10208 | 10.33849 ± 0.10196 | 10.43875 | 10.47483 | -0.86 | -1.34 | -0.48 |
| famB | 0.65 | 1600 | 24.166098 | 40.00000 (command) | 960.000 | 0.040527 | 0.650545 | 23.141098 | 23.120834 | 10.30054 ± 0.08595 | 10.29152 ± 0.08587 | 10.43875 | 10.46561 | -1.61 | -2.03 | -0.42 |

height cast: H*24 integer in every cell
largest |change in D|: 0.73; points with |change| > 0.5: 4 [('A2', 0.5, 100, -0.53), ('A2', 0.5, 400, -0.73), ('famB', 0.5, 100, -0.53), ('famB', 0.5, 400, -0.73)]

### Finite-size fit c_s = c_inf + b/sqrt(N), weights = mass scatter, source 260919_A2_cs_per_mass.csv

| eta | N values | P: c_inf ± | P: b ± | B: c_inf ± | B: b ± | A: c_inf ± | A: b ± | (A-B) c_inf / σ_B | (A-B) b / σ_B | KR(eta) | A: (c_inf - KR)/σ | chi2 B / A |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.02 | 100 400 900 1600 | 1.46140 ± 0.01023 | +0.1867 ± 0.1475 | 1.46140 ± 0.01023 | +0.1848 ± 0.1475 | 1.46144 ± 0.01023 | +0.1832 ± 0.1475 | +0.00 | -0.01 | 1.47247 | -1.08 | 2.075 / 2.075 |
| 0.05 | 100 400 900 1600 | 1.56988 ± 0.00291 | +0.0803 ± 0.0308 | 1.56988 ± 0.00291 | +0.0752 ± 0.0308 | 1.56992 ± 0.00291 | +0.0705 ± 0.0308 | +0.02 | -0.15 | 1.56662 | +1.14 | 0.446 / 0.446 |
| 0.10 | 100 400 900 1600 2500 | 1.75361 ± 0.00493 | -0.0140 ± 0.1350 | 1.75362 ± 0.00493 | -0.0255 ± 0.1350 | 1.75368 ± 0.00493 | -0.0372 ± 0.1350 | +0.01 | -0.09 | 1.74414 | +1.93 | 2.144 / 2.145 |
| 0.30 | 100 400 900 1600 2500 | 2.86354 ± 0.00803 | +0.4601 ± 0.1853 | 2.86369 ± 0.00804 | +0.3990 ± 0.1854 | 2.86315 ± 0.00799 | +0.3655 ± 0.1843 | -0.07 | -0.18 | 2.85150 | +1.46 | 2.008 / 1.985 |
| 0.50 | 100 400 900 1600 | 5.43030 ± 0.02368 | +0.9998 ± 0.4067 | 5.43131 ± 0.02369 | +0.7895 ± 0.4065 | 5.43258 ± 0.02370 | +0.3414 ± 0.4065 | +0.05 | -1.10 | 5.41622 | +0.69 | 0.387 / 0.388 |
| 0.60 | 100 400 900 1600 | 8.20534 ± 0.04708 | +1.3605 ± 1.0257 | 8.20717 ± 0.04664 | +0.9717 ± 1.0148 | 8.19146 ± 0.04660 | +0.9631 ± 1.0140 | -0.34 | -0.01 | 8.22261 | -0.67 | 0.664 / 0.663 |
| 0.65 | 100 400 900 1600 | 10.16323 ± 0.06872 | +6.4316 ± 1.5362 | 10.16681 ± 0.06909 | +5.8635 ± 1.5421 | 10.17090 ± 0.06941 | +4.2821 ± 1.5474 | +0.06 | -1.03 | 10.43875 | -3.86 | 0.608 / 0.617 |

fit parameters moving by more than 1 sigma (B): [(0.5, 0.05, -1.1), (0.65, 0.06, -1.03)]

**VERDICT (pre-registered rule, sec. 14.3): REGENERATE** -- 4 point(s) with |change in D| > 0.5 (largest 0.73); 2 fit parameter(s) beyond 1 sigma.
P -> B (divider thickness, not in the verdict): c_inf moves by 0.02: +0.00 σ, 0.05: +0.00 σ, 0.10: +0.00 σ, 0.30: +0.02 σ, 0.50: +0.04 σ, 0.60: +0.04 σ, 0.65: +0.05 σ

### Writing the corrected A2 tables (inputs: the _pre_boxtrunc_261002 copies)

  260919_A2_cs_per_mass.csv: 154 rows; B geometry reproduces the old c_s_mass to 8.3e-07; c_s_mass now at L_eff,true; eta = eta_true; new columns eta_rec, delta_sigma, L_eff_true
  260919_A2_cs_per_mass_famB.csv: 154 rows; B geometry reproduces the old c_s_mass to 8.3e-07; c_s_mass now at L_eff,true; eta = eta_true; new columns eta_rec, delta_sigma, L_eff_true
  260917_A2_cs_vs_N_extrapolation.csv: 7 rows, geometry A (L_eff = L0 - 2r - t/2 - delta/2, KR at eta_true per point, c_inf = KR(eta) + intercept of the deviation fit); same columns as before

**Regeneration record.**
- **Dated copies first.** Seven files were copied with `cp -n` and checked with `cmp`, each named `<name>_pre_boxtrunc_261002.<ext>`:
  - `260919_A2_cs_per_mass.csv` and `260919_A2_cs_per_mass_famB.csv`;
  - `260917_A2_cs_vs_N_extrapolation.csv`;
  - `260919_cs_vs_eta_lowdensity_zoom` and `260919_cs_vs_eta_N100_vs_A2`, each as `.png` and `.pdf`.

  The script reads its uncorrected inputs from these copies, so a re-run never corrects twice.
- **Tables.** `--write` wrote both per-mass tables:
  - `eta` now holds $\eta_{\rm true}$, and `c_s_mass` is at $L_{\rm eff,true}$;
  - the new columns are `eta_rec`, `delta_sigma` and `L_eff_true`;
  - gate: the B geometry reproduces the old `c_s_mass` to 8.3e-7.

  It also wrote the extrapolation table in geometry A, with the same columns as before.
- **Figures.** Regenerated by their original scripts with the environment `paper1_canonical_20260919.py` gives them: `lowdensity_zoom_20260917.py 260919_A2_cs_per_mass.csv 260919_cs_vs_eta_lowdensity_zoom` and `overlay_N100_vs_A2_20260915.py` (with `HD_A2_CSV=260919_A2_cs_per_mass_famB.csv`). The title suffix is now "N = 100 (A1 v2) and A2 corrected for box truncation (methods §14)"; the "A2 points not yet corrected" note is gone.
- **`paper1_canonical_20260919.py`.** Its suffix line is updated. Its old thickness-rescale step (260917/260916 → 260919, sources not on disk) now refuses to overwrite a table that carries `delta_sigma`.
- **One fix in `overlay_N100_vs_A2_20260915.py` (printed table only; the figure is written before it).** The table paired N = 900 and 1600 by exactly equal η. That fails once each point sits at its own $\eta_{\rm true}$.
  - It now pairs by the nearest η within 0.005: `c16 = next((p for p in a2[1600] if abs(p[0] - eta) < 5e-3), None)`.
  - The N = 1600 deviation is now taken at that point's own η (`kr16`), as the script's header always said.
- **Not regenerated (OPEN, not in the draft):**
  - `260917_A2_cs_vs_N.{png,pdf}`, the extrapolation figure;
  - `260916_A2_finite_size_forms.csv`, the bracketing forms; the draft quotes no number from it.

**Found, not fixed (DATA).** Both regenerated figures carry a stale footnote, "error bars = 1σ scatter of per-mass c_s" (`overlay_N100_vs_A2_20260915.py:99`, `lowdensity_zoom_20260917.py:107`). Since 2026-10-02 the plotted bar is the propagated `slope_with_errors` error, inflated by $\sqrt{\chi^2_{\rm red}}$ (`overlay…:42`, `lowdensity…:44`). The zoom script's printed column header "± scatter" (`:113`) is stale in the same way.

**The overlay figure's printed table, after the correction** (verbatim, `overlay_N100_vs_A2_20260915.py`):

A1 v2 has no density at exactly 0.10/0.30/0.50/0.60/0.65, so the N = 100 column below is A1 v2 at its OWN
nearest density, and every deviation is against KR at that point's own eta.
| A2 η | N = 900 | ± | dev KR [%] | N = 1600 | ± | dev KR [%] | nearest A1 v2 η | N = 100 | ± | dev KR [%] |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.02 | 1.4709 | 0.0124 | -0.11 | 1.4804 | 0.0085 | +0.54 | 0.019637 | 1.4773 | 0.0044 | +0.40 |
| 0.05 | 1.5655 | 0.0049 | -0.07 | 1.5796 | 0.0060 | +0.83 | 0.052374 | 1.5870 | 0.0026 | +0.80 |
| 0.10 | 1.7547 | 0.0036 | +0.60 | 1.7483 | 0.0029 | +0.23 | 0.112267 | 1.8104 | 0.0053 | +1.02 |
| 0.30 | 2.8691 | 0.0088 | +0.58 | 2.8712 | 0.0050 | +0.68 | 0.261799 | 2.6130 | 0.0077 | +1.57 |
| 0.50 | 5.5005 | 0.0917 | +1.39 | 5.4349 | 0.0629 | +0.22 | 0.523599 | 6.0881 | 0.0171 | +2.63 |
| 0.60 | 8.2661 | 0.0440 | +0.39 | 8.2245 | 0.0538 | -0.12 | 0.590895 | 7.9110 | 0.0912 | +0.31 |
| 0.65 | 10.3385 | 0.1020 | -1.30 | 10.2915 | 0.0859 | -1.66 | 0.652234 | 10.9853 | 0.0528 | +4.14 |

**F3 draft audit, printed by `python3 hspist3/validation/paper1_draft_audit_20261014.py --a2` before any edit** (verbatim):

sec. 14.3 numbers: max |change in D| = 0.73; max |change in c_inf| = 0.34 sigma; |change in b| > 1 sigma at eta = [0.5, 0.65] (-1.10, -1.03); max 1.10

### Paper 1 draft audit, A2-dependent text (methods sec. 14.3), old -> new (printed BEFORE editing)

| tex line(s) | old text | new text | definition | old table gives | old text reproduced / claim holds? | new table gives | action |
|---|---|---|---|---|---|---|---|
| 272 | `1.4614 \pm 0.0102` | `1.4614 \pm 0.0102` | c_inf at eta = 0.02 (260917 extrapolation) | 1.46140 ± 0.01023 | yes | 1.46144 ± 0.01023 | UNCHANGED |
| 272 | `$1.4725$` | `$1.4725$` | KR at eta = 0.02 | 1.47247 | yes | 1.47247 | UNCHANGED |
| 273 | `($-1.1\sigma$)` | `($-1.1\sigma$)` | (c_inf - KR)/sigma at eta = 0.02 | -1.08 | yes | -1.08 | UNCHANGED |
| 273 | `1.5699 \pm 0.0029` | `1.5699 \pm 0.0029` | c_inf at eta = 0.05 (260917 extrapolation) | 1.56988 ± 0.00291 | yes | 1.56992 ± 0.00291 | UNCHANGED |
| 273 | `$1.5666$` | `$1.5666$` | KR at eta = 0.05 | 1.56662 | yes | 1.56662 | UNCHANGED |
| 273 | `($+1.1\sigma$)` | `($+1.1\sigma$)` | (c_inf - KR)/sigma at eta = 0.05 | +1.12 | yes | +1.14 | UNCHANGED |
| 118 | `The larger systems of \S\ref{sec:finitesize} carry the same shortfall and are not yet corrected.` | `The same rule, pre-registered separately, was applied to the larger systems of \S\ref{sec:finitesize}: the correction moves no point of the size ladder by more than $0.73$ of its error bar and no fitted $c_\infty$ by more than $0.34$ of its error, but it lowers the coefficient of $1/\sqrt N$ by up to $1.10$ of its error (at $\eta = 0.50$ and $0.65$), so the size ladder was regenerated too.` | Methods paragraph: the A2 status sentence | not corrected | yes | corrected (sec. 14.3, REGENERATE) | EDIT |
| 111 | `For every $N = 100$ density we therefore use` | `For every density and system size we therefore use` | Methods paragraph: scope of the correction | N = 100 only | yes | all A1 and A2 cells | EDIT |
| 111 | `long, as the recorded box centre confirms at every density to` | `long, as the recorded box centre confirms at every $N = 100$ density to` | Methods paragraph: the Center_X check covers the A1 runs only (sec. 14.1) | every density | yes | every N = 100 density | EDIT |
| 120 | `% TODO-source: methods sec. 14 (rule, 4db8c9d), 14.1 (Center_X check; 615561c), 14.2; validation/paper1_boxtrunc_20261014.py` | `% TODO-source: methods sec. 14 (rule, 4db8c9d), 14.1 (Center_X check; 615561c), 14.2; validation/paper1_boxtrunc_20261014.py; %   sec. 14.3 (A2 rule 0ddefa9, results 14.3.1); validation/paper1_A2_boxtrunc_261002.py` | source comment | - | yes | + sec. 14.3 | EDIT |
| 264 | `larger systems agree with it within their errors.` | `larger systems agree with it within their errors.` | zoom caption: D (A geometry, plotted sigma) of the N = 900/1600 points shown (eta <= 0.15): 0.02/N900: -1.43, 0.02/N1600: +0.18, 0.05/N900: +1.77, 0.05/N1600: +1.42, 0.10/N900: +2.96, 0.10/N1600: +1.39 | (claim) | **NO** | largest |D| = 2.96 at eta = 0.10, N = 900 | FLAG -- claim not supported at 1 sigma; not edited (wording is the plan author's call) |

Not edited: 'a size ladder to N = 2500 shows the offset is finite size' (abstract) and 'The N = 100 offset closes as the box grows' (sec. finitesize) -- after the correction c_inf sits at -1.08 and +1.14 sigma from KR at eta = 0.02 and 0.05, as before.

**Applied** with `--a2 --apply`: the four EDIT rows only.
- The Methods paragraph's A2 sentence now states the § 14.3 result.
- Its scope reads "every density and system size".
- The Center_X check is scoped to "every $N = 100$ density" (it was run on A1 only, § 14.1).
- Its source comment gains § 14.3.

The draft compiles in two passes with no undefined references (scratchpad build).

**FLAG for the plan author (pre-existing, not edited).** The zoom caption says "larger systems agree with it within their errors." The table above shows that is not true at 1σ for every point: N = 900 at η = 0.10 sits at $D = +2.96$, and was at +3.05 before the correction. The wording is the plan author's call.

#### 14.3.2 Draft dependencies after the corrections: 261001 figures, draft pdf, fresh-clone build (2026-10-02, Task G)

**G1, tracked first.** `paper1_figures_20261001.py` and the four `261001_p1_*` figures (png and pdf) were committed exactly as they were (e10f0d2). Before that, the draft included them while git did not track them.

**G1, redrawn.** The ladder-line and slow-mode figures were redrawn at the corrected geometry of the $\eta_{\rm rec} = 0.1122$ cell. Dated copies were made first (`cp -n`, `cmp`): `261001_p1_massladder_line_pre_boxtrunc_261002` and `261001_p1_slowmode_pre_boxtrunc_261002`, each as `.png` and `.pdf`.
- **Ladder line.** $x_M$ at $L_{\rm eff,true} = 33.9542$, KR at $\eta_{\rm true} = 0.112267$. The title carries "corrected for box truncation (methods §14)". The data label gives both η values.
- **A pre-existing inconsistency was found and fixed (DATA).** The 2026-10-01 figure did not use the canonical estimator:
  - it took the whole record, with floor bin $k = N_{\rm cyc}/2.5$, instead of the first TD = 200 periods with $k = {\rm TD}/X$;
  - its ± was ${\rm sd}_{\rm mass}/\sqrt n$, not `c_s_err_scaled`.
  - So it showed $c_s = 1.8098 \pm 0.0015$ (read from the old figure's legend) against the table's $1.81155 \pm 0.00529$.
  - It now uses `cell()`'s per-trajectory estimator and `slope_with_errors`, both imported. It draws nothing unless the gate below passes.
- **Slow mode.** The title gives $\eta_{\rm true}$ and carries the correction note. The periodogram is now computed on the canonical 200-period record, which gives the same ν, 0.002627.
  - The floor $\nu_{\rm pred}/2.5$ and the $\nu_{\rm pred}$ line are **not** moved. $\nu_{\rm pred}$ is the binary's own prediction at the recorded geometry, and the estimator anchors its floor on it, so this is where the floor really was. The legend says so.

**Printed by `python3 hspist3/validation/paper1_figures_20261001.py ladder slowmode`** (verbatim):

Paper 1 figures -> /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/0000_PLAN_OVERALL/paper1_speedofsound/experiments/final
     gate: recorded geometry 1.81155 (pre-correction table 1.81155); corrected 1.81044 (table 1.81044); error 0.005291 (table c_s_err_scaled 0.005291) -> PASS
  wrote 261001_p1_massladder_line.png/.pdf
     c_s = 1.81044 +- 0.00529   KR = 1.79216   ratio 1.01020
  wrote 261001_p1_slowmode.png/.pdf
     M=1000: dt=11.9999, window=160 samples, nu_pred=0.002601, nu=0.002627

The draft's ladder caption ("sits $1.0\,\%$ above Kolafa–Rottner") still holds: the ratio printed above is 1.01020.

**The two relative-quantity figures stay as they were (INFERENCE from their definitions).**
- **`261001_p1_estimator_floor`** plots each density's deviation of the floor estimator from the windowed reference in percent. Both are computed from the same traces at the same geometry, so the one-factor-per-density correction cancels in the ratio.
- **`261001_p1_massladder_residuals`** plots residuals $\bar\nu_M/(c_s x_M) - 1$ about each density's own ladder line. The correction scales $x_M$ and $c_s$ by inverse factors, so the product $c_s x_M$, and with it every residual, is unchanged.

**G2.** `writeup/paper1_draft.pdf` was rebuilt from the corrected tex in two passes (b5a4c9a):
- 6 pages, 636139 bytes, 0 undefined references or citations;
- 6 warnings, all pre-existing: 2 font-shape, and 4 hyperref warnings about math in the section title at tex line 342, which is dropped from the PDF bookmark only;
- 1 overfull hbox (lines 358–365) and 1 underfull vbox.

**G3, fresh-clone build: PASS.**
- **Setup.** A `--shared` clone of HEAD b5a4c9a in the scratchpad, sparse to `0000_PLAN_OVERALL/paper1_speedofsound/` (1759 tracked files), compiled in two passes.
- **Result.** 6 pages, 636139 bytes, the same as G2. 0 missing files, 0 undefined references, 6 warnings.
- **Inputs.** The draft includes ten graphics files, all tracked:
  - `260919_cs_vs_eta`, `260919_cs_vs_eta_N100_vs_A2`, `260919_cs_vs_eta_lowdensity_zoom`;
  - `260922_apparatus_paper.png`, `260922_roman2002_remapped_vs_KR`;
  - the four `261001_p1_*`;
  - `261002_p1_melting_region`.

  Nothing needed adding.

**FLAG (caption labels, not edited).** The ladder-line and slow-mode captions name the cell "$\eta = 0.1122$", its recorded value. The redrawn figure titles now show $\eta_{\rm true} = 0.1123$; the ladder legend gives both values. Whether the captions should read "$\eta = 0.1123$ (recorded 0.1122)" is the plan author's call.

### 14.4 Date note: series dates against machine dates (2026-10-02, machine date)

**What happened (DATA, table below).** From at least c6e5f62 (committed 2026-09-19) to 9ff3491 (committed 2026-10-02 15:51), the dates written into new file names (`2609xx_…`, `2610xx_…`), `##CHRIS` headers and STATUS timestamps came from a series counter that ran ahead of the machine clock.
- The offset was not constant. Per commit day it was:

  | commit day | 09-19 | 09-20 | 09-22 | 09-23 | 09-25 | 09-30 | 10-01 | 10-02 |
  |---|---|---|---|---|---|---|---|---|
  | offset [days] | +2 | +4 | +3 to +7 | +7 to +12 | +14 to +15 | +10 | +11 | +11 to +12 |

- The "8 days on 2026-09-23" in the plan is the cluster-file case: headers say 2026-10-01, file mtime 2026-09-23. The commits of that day range from +7 to +12.
- STATUS timestamps in that range carry the series date with the machine's time of day.
- The dates inside this file's § 14 headings ("2026-10-14") are series dates too. The real date is 2026-10-02.

**Rules from now on.**
- **Git commit timestamps are authoritative.** Where a date in a file name, a header or a STATUS line disagrees with the commit that added it, the commit date is the real one.
- **Nothing is renamed.** Scripts, drafts and earlier sections refer to the existing names.
- **Machine dates from 4adedde (2026-10-02 18:32) on:** `date +%y%m%d` in new file names, `date` in STATUS lines. These commits show offset +0.
- **The two −1 rows** are files written earlier under the series, with header 2026-10-01, and committed later: the cluster runbook and sbatch drafts (fa300e6) and the 261001 figures (e10f0d2).
- **Commits with no date** in what they added are left blank.

**Printed by `python3 hspist3/validation/date_series_audit_261002.py`** (verbatim; series date = the latest date written into what the commit added):

### Series date written into each commit vs its git commit date (every commit since 2026-09-20, oldest first)

| # | hash | commit date (HST) | series date | offset [days] | latest series date found in | subject |
|---|---|---|---|---|---|---|
| 1 | bf8588c | 2026-09-20 22:20 | 2026-09-24 | +4 | ##CHRIS | Level 3 v6: window validity, box-input model, verdict |
| 2 | 7760690 | 2026-09-22 12:17 | 2026-09-25 | +3 | ##CHRIS | Level 3 closed; measured F(L); Level 4 design and pilot |
| 3 | ba83144 | 2026-09-22 12:42 | 2026-09-23 | +1 | file name | Paper 1 draft v1; status figures; Level 3 wording; tau_heat estimator |
| 4 | 2df7c81 | 2026-09-22 22:03 | 2026-09-28 | +6 | ##CHRIS | Level 4 thermal: adiabatic piston stage 1 measured, stage 2 needs more seeds |
| 5 | 849cd77 | 2026-09-22 23:12 | 2026-09-29 | +7 | ##CHRIS | Level 4: tau_T measured from equilibrium fluctuations; record length is the constraint |
| 6 | 767cd9c | 2026-09-23 10:46 | 2026-09-30 | +7 | ##CHRIS | Level 4: the equilibrium oscillation is Paper 1's divider mode; --trace-every gated |
| 7 | fdcf559 | 2026-09-23 12:17 | 2026-10-02 | +9 | ##CHRIS | Paper 1 error bars populated; melting-region figure; Md10 long run resolves the isobar at  |
| 8 | 5ce3522 | 2026-09-23 20:44 | 2026-10-03 | +10 | ##CHRIS | Level 4 ladder: exponent withdrawn, heavy end is an estimator limit; A2 error bars; Eq 51  |
| 9 | 2777d9c | 2026-09-23 20:44 | 2026-10-04 | +11 | ##CHRIS | Level 4 ladder rerun at 65 tau_true; eta read from recorded t=1.0; A2 surface scaling re-d |
| 10 | 0f299df | 2026-09-23 23:58 | 2026-10-05 | +12 | ##CHRIS | Level 4: ladder as tau_T = M g(R); R-collapse test; Cencini fixed-R limit verified |
| 11 | 3af501c | 2026-09-24 00:05 | 2026-10-06 | +12 | ##CHRIS | Level 4 R-collapse pass 2: PRE-REGISTRATION committed before any data exists |
| 12 | c7b95f6 | 2026-09-25 10:03 | 2026-10-09 | +14 | ##CHRIS | Level 4: R-collapse pass 2 NOT RESOLVED; mode ladder; demo; 4b model test failed, ledger p |
| 13 | 364268e | 2026-09-25 10:03 | 2026-09-26 | +1 | file name | housekeeping: KOA docs, COWORK notes, 260908/260926 notes, knowledge-transfer edit |
| 14 | 50e47b3 | 2026-09-25 10:07 | 2026-10-10 | +15 | ##CHRIS | redo follow-up: files dropped by the explicit-path rewrite |
| 15 | 1f7be49 | 2026-09-25 10:08 | - |  | - | Restore gui_config.c; commit the Makefile/KOA/run-script notes; ignore _commit_scripts/ |
| 16 | 0ab395b | 2026-09-25 10:28 | 2026-10-09 | +14 | file name | Item 6 closed: f(R) not universal at 2.3 sigma; second variable degenerate at fixed eta |
| 17 | 05215ea | 2026-09-30 11:36 | 2026-10-10 | +10 | STATUS | Piston v2 not adopted (branch RED); -ffp-contract=off in release/koa; --version provenance |
| 18 | eb7387e | 2026-09-30 11:38 | 2026-10-10 | +10 | ##CHRIS | 261007 section 4.4 withdrawn; B3 PRE-REGISTERED (section 3B) before any record exists; eff |
| 19 | f6cde70 | 2026-09-30 11:44 | 2026-10-10 | +10 | ##CHRIS | Level 4b closed: full grid on one binary (05215ea); T(x) fails on scale by 3-5x, survives  |
| 20 | 33cf95a | 2026-09-30 11:50 | 2026-10-10 | +10 | STATUS | B3: tau_T confirmed out of equilibrium at 0.85 sigma (29 626 +- 12 246 vs 40 079); a consi |
| 21 | f6bf799 | 2026-09-30 11:53 | 2026-10-10 | +10 | ##CHRIS | Plan: Level 4 as measured, one pass (sec:l4measured); B3 figure |
| 22 | a2ccaac | 2026-09-30 11:54 | 2026-10-10 | +10 | ##CHRIS | Efficiency-map runner per 261010 section 1 -- written and syntax-checked, NOT launched (wa |
| 23 | 84d7a49 | 2026-10-01 15:17 | 2026-10-12 | +11 | ##CHRIS | Paper 1 item 1b PRE-REGISTERED before any fit (linewidth decomposition); 1c linewidth conv |
| 24 | 97938a7 | 2026-10-01 15:20 | 2026-10-12 | +11 | ##CHRIS | Paper 1 items 1a and 1b results |
| 25 | 463ccf9 | 2026-10-01 15:26 | 2026-10-12 | +11 | STATUS | Paper 1 confinement campaign PRE-REGISTERED, not launched (261012 sec. 1); status log for  |
| 26 | b788e82 | 2026-10-01 17:15 | 2026-10-12 | +11 | ##CHRIS | Efficiency map: amendments A1-A3 committed BEFORE launch (261010 sec. 1.8); analysis scrip |
| 27 | 1968d08 | 2026-10-01 17:18 | 2026-10-12 | +11 | ##CHRIS | Confinement amendments C1-C2 (261012 sec. 1.9), not launched; Roman citation corrected in  |
| 28 | 0504d68 | 2026-10-01 17:30 | 2026-10-12 | +11 | STATUS | Efficiency map results (261010 sec. 2): 336/336 runs; ledger 1.35e-5 kT; Level 3 reproduct |
| 29 | f3c6209 | 2026-10-02 14:00 | - |  | - | effmap: over-pressure test pre-registration (§3) |
| 30 | c849b93 | 2026-10-02 14:06 | 2026-10-13 | +11 | ##CHRIS | effmap: over-pressure test results, position-settled gate (§3.1–3.2) |
| 31 | 4f8d854 | 2026-10-02 14:08 | 2026-10-13 | +11 | ##CHRIS | Disk report 2026-10-13 (read-only; nothing deleted or moved): 20 largest campaign director |
| 32 | 930a92d | 2026-10-02 14:31 | 2026-10-13 | +11 | STATUS | paper1 confinement: KOA smoke test, sbatch generation, local gates (§1.10) |
| 33 | 4db8c9d | 2026-10-02 14:52 | - |  | - | paper1: box-truncation correction pre-registration (§9) |
| 34 | 5190846 | 2026-10-02 14:52 | - |  | - | effmap: gate G2 pre-registration (§3.4) |
| 35 | 615561c | 2026-10-02 14:57 | 2026-10-14 | +12 | ##CHRIS | paper1: box-truncation results (§9.1), box_width_sigma column + truncation warning |
| 36 | a5c4e19 | 2026-10-02 14:58 | - |  | - | effmap: over-pressure timeline disclosure and tag split (§3.3) |
| 37 | 1ff8654 | 2026-10-02 15:00 | 2026-10-14 | +12 | ##CHRIS | effmap: gate G2 results (§3.5) -- valid (6/6 controls), 30/36 cells PASS; figures relabell |
| 38 | 06e3df6 | 2026-10-02 15:01 | 2026-10-14 | +12 | ##CHRIS | paper1 confinement: smoke-test gate width (§1.10.1) -- 0.05150 was the mass SD; gate now 2 |
| 39 | 3bf7b55 | 2026-10-02 15:01 | 2026-10-14 | +12 | STATUS | status: box-truncation REGENERATE, effmap §3.3/G2, smoke-test gate width |
| 40 | e823187 | 2026-10-02 15:46 | 2026-10-14 | +12 | ##CHRIS | paper1: canonical table/figure regenerated with box-truncation correction (§14.2), draft n |
| 41 | bbdf135 | 2026-10-02 15:50 | 2026-10-14 | +12 | STATUS | status: clean rebuild e823187 (release, -ffp-contract=off), determinism IDENTICAL, new bui |
| 42 | 9ff3491 | 2026-10-02 15:51 | 2026-10-14 | +12 | STATUS | methods §14 + status: summary.csv header rotation -- Level 3-era runners get new output di |
| 43 | fa300e6 | 2026-10-02 15:52 | 2026-10-01 | -1 | ##CHRIS | cluster: KOA first-hour runbook and first sbatch drafts (recovered) |
| 44 | 4adedde | 2026-10-02 18:32 | - |  | - | build: track experiment_validation.c/.h — provenance gap closed |
| 45 | 38e191a | 2026-10-02 18:33 | 2026-10-02 | +0 | STATUS | status: clean rebuild 4adedde, three-binary cross-check IDENTICAL; E1 inventory and E5 scr |
| 46 | f956b4b | 2026-10-02 18:41 | 2026-10-02 | +0 | ##CHRIS | cluster: KOA placeholders filled, noexec build path, runsheet |
| 47 | 0ddefa9 | 2026-10-02 18:44 | - |  | - | paper1: A2 correction pre-registration (§14.3) |
| 48 | 47dbb5f | 2026-10-02 18:52 | 2026-10-02 | +0 | ##CHRIS | paper1: A2 size-ladder box-truncation correction (§14.3.1, REGENERATE), tables/figures reg |
| 49 | e10f0d2 | 2026-10-02 18:52 | 2026-10-01 | -1 | ##CHRIS | paper1: track 261001 figures and script (draft dependency) |
| 50 | 94f7d6c | 2026-10-02 18:55 | 2026-10-02 | +0 | ##CHRIS | paper1: ladder-line and slow-mode figures redrawn at the corrected geometry, canonical est |
| 51 | b5a4c9a | 2026-10-02 18:55 | - |  | - | paper1: draft pdf rebuilt from the corrected tex (6 pages, 0 undefined references) |
| 52 | 300b743 | 2026-10-02 18:56 | 2026-10-02 | +0 | STATUS | paper1: fresh-clone draft build PASS; figures/pdf record (§14.3.2) and STATUS |
| 53 | 52fee6d | 2026-10-02 18:57 | 2026-10-02 | +0 | STATUS | paper1: nofuse check selects by eta_rec (L_eff,true geometry); re-run 0.73 sigma PASS, unc |

commits: 53; carrying a series date: 45; series date >= 2 days ahead of the commit: 35
first affected: bf8588c (2026-09-20 22:20, series 2026-09-24, +4 d); last affected: 9ff3491 (2026-10-02 15:51, series 2026-10-14, +12 d)
offset by commit day: 2026-09-20: +4 to +4; 2026-09-22: +3 to +7; 2026-09-23: +7 to +12; 2026-09-24: +12 to +12; 2026-09-25: +14 to +15; 2026-09-30: +10 to +10; 2026-10-01: +11 to +11; 2026-10-02: +11 to +12
before 2026-09-20 (commits since 2026-08-01): 3 with a series date >= 2 days ahead; earliest c6e5f62 (2026-09-19 14:56, series 2026-09-21, +2 d)

### 14.5 Draft housekeeping after the review (2026-10-02, Task M; the plan author's decisions)

**Printed by `python3 hspist3/validation/paper1_draft_audit_20261014.py --m` before any edit** (verbatim):

zoom points N = 900/1600 (A geometry, plotted sigma): eta 0.02 N 900: D = -1.43, eta 0.02 N 1600: D = +0.18, eta 0.05 N 900: D = +1.77, eta 0.05 N 1600: D = +1.42, eta 0.10 N 900: D = +2.96, eta 0.10 N 1600: D = +1.39
largest |D| 2.96 at eta = 0.10, N = 900; all others |D| <= 1.77; worked cell eta_true = 0.112267 -> 0.1123

### Task M text changes, old -> new (printed BEFORE editing)

| kind | where | tex line | old | new | why |
|---|---|---|---|---|---|
| TEX | zoom caption (M1) | 267 | `larger systems agree with it within their errors.}` | `larger systems agree with it within their errors, except $N = 900$ at $\eta = 0.10$ ($3\sigma$).}` | D = +2.96; the other five points |D| <= 1.77 |
| TEX | slow-mode caption (M2; first mention in the text) | 182 | `$\eta = 0.1122$, with its five-period running mean` | `$\eta = 0.1123$ (recorded $0.1122$), with its five-period running mean` | eta_true = 0.112267; figure title says 0.1123 |
| TEX | ladder-line caption (M2) | 230 | `One worked mass ladder, $\eta = 0.1122$, nine divider masses` | `One worked mass ladder, $\eta = 0.1123$, nine divider masses` | figure title says 0.1123, legend gives both |
| SCRIPT | lowdensity_zoom_20260917.py:107 footnote (M3) | - | `error bars = 1σ scatter of per-mass c_s` | `error bars = 1σ propagated error of the through-origin slope (σ_ν = sd/√n per mass, × √χ²_red when > 1)` | the plotted bar is slope_with_errors' err_scaled (:44) |
| SCRIPT | overlay_N100_vs_A2_20260915.py:99 footnote (M3) | - | `error bars = 1σ scatter of per-mass c_s` | `error bars = 1σ propagated error of the through-origin slope (σ_ν = sd/√n per mass, × √χ²_red when > 1)` | the plotted bar is slope_with_errors' err_scaled (:42) |
| SCRIPT | lowdensity_zoom_20260917.py:113 printed header (M3) | - | `| η | N | c_s | ± scatter | KR | dev [%] |` | `| η | N | c_s | ± err (propagated) | KR | dev [%] |` | same quantity |
| SCRIPT | lowdensity_zoom_20260917.py:3-4 docstring (M3) | - | `error bars are the 1 sigma scatter of the per-mass c_s.` | `error bars are the 1 sigma propagated error of the through-origin slope (slope_with_errors: sd/sqrt(n) per mass, x sqrt(chi2_red) when > 1).` | same quantity |
| SCRIPT | overlay_N100_vs_A2_20260915.py:4 docstring (M3) | - | `error bars are the 1 sigma scatter of the per-mass c_s.` | `error bars are the 1 sigma propagated error of the through-origin slope (slope_with_errors: sd/sqrt(n) per mass, x sqrt(chi2_red) when > 1).` | same quantity |

`--m --check-scripts` printed:

    lowdensity_zoom_20260917.py:107 footnote (M3): old ABSENT, new present
    overlay_N100_vs_A2_20260915.py:99 footnote (M3): old ABSENT, new present
    lowdensity_zoom_20260917.py:113 printed header (M3): old ABSENT, new present
    lowdensity_zoom_20260917.py:3-4 docstring (M3): old ABSENT, new present
    overlay_N100_vs_A2_20260915.py:4 docstring (M3): old ABSENT, new present

**Applied.**
- The three TEX rows were applied with `--m --apply`.
- The five SCRIPT rows were applied exactly as printed. `--m --check-scripts` reports every old string absent and every new string present.

**Layout change (disclosed; not in the list above).** The corrected footnote no longer fits on one line, so its left end was cut off at the figure edge. In both scripts:
- the footnote is now two lines, with "data: …" on the second;
- `fig.tight_layout(rect=(0, 0.015, 1, 1))` became `rect=(0, 0.03, 1, 1)`.

Both figures were then regenerated with the environment `paper1_canonical_20260919.py` gives them. Their printed tables are identical to the § 14.3.1 run, apart from the zoom table's header.

**M1, wording note (DATA, not edited).** The other five larger-system points in the zoom have |D| ≤ 1.77 (list above). So "agree with it within their errors" holds at 2σ, not at 1σ. "Within $2\sigma$" would be the exact wording; that is the plan author's call.

**M4.** `260917_A2_cs_vs_N.{png,pdf}`: dated copies `_pre_boxtrunc_261002` were made first (`cp -n`, `cmp`). The figure was then redrawn by `paper1_A2_boxtrunc_261002.py --figure`, in geometry A, with the layout and colours of `analyze_A2_X2p5_20260914.py:128–147`.
- **Why not the original script.** It recomputes everything at $L_{\rm eff} = L_0 - 2r$ and would overwrite the corrected extrapolation table.
- **What each point is.** $c_s - c_s^{\rm KR}(\eta_{\rm true}) + c_s^{\rm KR}(\eta)$, i.e. the § 14.3 deviation shifted by the constant $c_s^{\rm KR}(\eta)$ (DERIVATION), so the drawn line $c_\infty + b/\sqrt N$ is the fit itself.
- **Printed:** `figure gate: plotted c_inf, b and errors equal 260917_A2_cs_vs_N_extrapolation.csv to its printed digits -> PASS`. The analysis part of that run is identical, line for line, to the committed § 14.3.1 run.

**M5.**
- **Local build.** `paper1_draft.pdf` rebuilt in two passes: 6 pages, 638006 bytes, 0 undefined references, 6 warnings (all pre-existing, § 14.3.2), 1 overfull hbox.
- **Fresh clone of 15c8b39** (sparse to `paper1_speedofsound/`, 1774 files): 6 pages, 638006 bytes, 0 missing files, 0 undefined references, 6 warnings. **PASS.**

    figure gate: plotted c_inf, b and errors equal 260917_A2_cs_vs_N_extrapolation.csv to its printed digits -> PASS
    wrote 260917_A2_cs_vs_N.png/.pdf (geometry A)

## 15. Kolafa–Rottner module: independent check (2026-10-02, Task L)

**Verdict by the stated rule: FAIL, stopped; the module is not edited.** The verdict is relative to the reference fixed in the script before the module was read. The module itself is an exact implementation of one of the paper's three published equations; which one Paper 1 should use is a decision for the plan author.

**What the paper says.** [SOURCE: Kolafa & Rottner, Mol. Phys. 104, 3435–3441 (2006); PDF in `ZZZ_PAPER/`; every digit below was confirmed against the PDF's text layer.]
- **Form.** Eq. (7), p. 3437: $Z(y) = \sum_i A_i x^i$ with $x = y/(1-y)$, where $y$ is the packing fraction.
- **Three fitted equations,** § 3.2, pp. 3438–3439: $\rho_{\max} = 0.88$ ($s = 0.724$), $0.89$ ($s = 0.966$) and $0.90$ ($s = 0.927$).
  - $\rho = N\sigma^2/A$, so $\eta = (\pi/4)\rho$, and $\rho_{\max} = 0.88 \Leftrightarrow \eta = 0.6912$, $0.90 \Leftrightarrow \eta = 0.7069$.
  - Of the 0.90 equation the paper says: "region $\rho \in [0.89, 0.90]$ of this equation may be affected by finite-size effects".
  - Fig. 3 marks all three versions as the best equations.

**What the module implements (DATA, table below).** The $\rho_{\max} = 0.90$ fit: Z to 4.2e-16, Z′ to 1.5e-16, $c_s$ to 3.4e-16 relative. Its own comment says so: `plot_speed_of_sound_edmd.py:620` reads "Kolafa & Rottner (2006), rho_max=0.90 fit". It draws KR only to `KR2006_PLOT_ETA_MAX = 0.690`.

**The reference I chose.** It was the $\rho_{\max} = 0.88$ fit, because the project documents KR as "valid to η ≈ 0.69", which is exactly that fit's range:
- `writeup/papers/README.md`;
- the figure legends "Kolafa–Rottner 2006 (valid to η ≈ 0.69)", `paper1_canonical_20260919.py:74` and `paper1_thickness_correction_20260918.py:87`.

Against it the module differs by up to 1.4e-6 in Z, 7.1e-5 in Z′ and 1.0e-5 in $c_s$ at η = 0.05–0.65. INFERENCE: that is small against every error bar in Paper 1. The question is which published equation the paper names, and what its "fitted range" sentences mean.

**Draft text that depends on the decision.** None of it is edited.
- **Line 255:** "beyond η ≈ 0.69 … no fluid reference exists".
- **Lines 342–345:** "fitted below the transition … above η ≈ 0.69 measures the extrapolation". The implemented fit extends to η = 0.707.
- **Lines 383–387,** the melting caption: "solid inside its fitted range … its fitted Z′ runs 48 at η = 0.67, 28 at 0.69, −9 at 0.700 and +3386 at 0.720 … The negative Z′ already inside the stated fitted range …". These Z′ values come from the 0.90 fit; with the 0.88 fit, the values beyond its range would differ.
- **Lines 398, 404:** "above 0.69".

**Options (decision for the plan author).**
- **(a) Keep the 0.90 fit,** as implemented and commented. Name it in the draft and README as "the $\rho_{\max} = 0.90$ fit, to η = 0.707, compared with our data only to η = 0.69". Re-run L1 with that version as the stated reference.
- **(b) Switch to the 0.88 fit,** valid to η = 0.691. Every KR value moves by ≤ 1e-5 relative at η ≤ 0.65; the KR columns and figures are regenerated; the melting caption's Z′ values change.

**Printed by `python3 hspist3/validation/paper1_kr_sanity_261002.py`** (verbatim):

module functions: Z_kolafa_rottner_2006(eta: 'np.ndarray') -> 'np.ndarray', dZ_kolafa_rottner_2006(eta: 'np.ndarray') -> 'np.ndarray', cs_adiabatic_2d_monatomic(Z: 'np.ndarray', dZ: 'np.ndarray', eta: 'np.ndarray', *, kbt: 'float', m: 'float') -> 'np.ndarray'
independent reference: rho_max = 0.88 version (eta_max = 0.6912)

| eta | Z paper | Z module | rel diff | Z' paper (analytic) | Z' complex-step | Z' module | rel diff | c_s paper | c_s module | rel diff | tests_20260913.kr_cs | rel diff vs paper |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.05 | 1.10838774198968 | 1.10838774198985 | 1.6e-13 | 2.34761773305187 | 2.34761773305187 | 2.3476177330678 | 6.8e-12 | 1.56661801829139 | 1.56661801829183 | 2.8e-13 | 1.56661801830068 | 5.9e-12 |
| 0.30 | 2.06326085484772 | 2.06326087957605 | 1.2e-08 | 6.03583224566138 | 6.03583224566138 | 6.03583334839163 | 1.8e-07 | 2.8515006371546 | 2.85150071739119 | 2.8e-08 | 2.85150071752567 | 2.8e-08 |
| 0.50 | 4.10635800241687 | 4.10636394336527 | 1.4e-06 | 16.7336733246343 | 16.7336733246343 | 16.7336836142638 | 6.1e-07 | 5.41621368750783 | 5.41621921508607 | 1.0e-06 | 5.41621921569349 | 1.0e-06 |
| 0.65 | 8.40804926299822 | 8.4080435277207 | 6.8e-07 | 45.9481747610299 | 45.9481747610299 | 45.9449284500896 | 7.1e-05 | 10.4388531585933 | 10.4387471941955 | 1.0e-05 | 10.4387471910944 | 1.0e-05 |

largest relative difference module vs paper (Z, Z', c_s; and analytic vs complex-step Z'): 7.1e-05
tests_20260913.kr_cs vs paper: largest relative difference 1.0e-05 (above 1e-10; reported, not part of the verdict -- see how kr_cs forms Z')

for information (NOT the verdict) -- the other two published versions against the module, largest relative difference over the four eta:
| version | Z | Z' (analytic) | c_s |
|---|---|---|---|
| rho_max = 0.89 | 6.8e-07 | 7.2e-06 | 1.5e-06 |
| rho_max = 0.9 | 4.2e-16 | 1.5e-16 | 3.4e-16 |
module comment, plot_speed_of_sound_edmd.py:620: '# Kolafa & Rottner (2006), rho_max=0.90 fit.  Their x is eta/(1-eta).'

**VERDICT: FAIL -- STOP, the module is NOT edited** (criterion: relative difference <= 1e-10)

**L2, printed by `python3 hspist3/validation/untracked_imports_scan_261002.py`** (verbatim):

tracked scripts scanned: 98 (of 98); local modules known: 123
| untracked module | file(s) | imported by (tracked) |
|---|---|---|

untracked modules imported by tracked scripts: 0 -- none

### 15.1 The fit, named and bounded (2026-10-02, Task P)

**Decision (plan author).** Paper 1 uses the Kolafa–Rottner 2006 $\rho_{\max} = 0.90$ fit (Eq. 7, coefficients of § 3.2), which is fitted to $\eta \le 0.7069$. Data are compared with it only for $\eta \le 0.69$, and the module stays as it is.

**Reasoning.** The three published fits differ by at most $10^{-5}$ in $c_s$ up to $\eta = 0.65$ (§ 15 table), so a switch would only churn every KR column. The 0.90 fit reaches furthest toward melting, but anything it gives beyond $\eta = 0.7069$ is an extrapolation, so no $Z'$ from there is quoted any more.

**P1, re-run with the 0.90 fit as reference: PASS (DATA).**
- The largest relative difference between module and paper over Z, Z′ and $c_s$ is 4.6e-16 (criterion 1e-10).
- `--ref 0.88` reproduces the § 15 run exactly (table, worst value and verdict).
- `tests_20260913.kr_cs` differs from the analytic value by up to 3.0e-10, because it takes Z′ by finite difference. That is reported, not part of the verdict.

Printed by `python3 hspist3/validation/paper1_kr_sanity_261002.py` (verbatim):

module functions: Z_kolafa_rottner_2006(eta: 'np.ndarray') -> 'np.ndarray', dZ_kolafa_rottner_2006(eta: 'np.ndarray') -> 'np.ndarray', cs_adiabatic_2d_monatomic(Z: 'np.ndarray', dZ: 'np.ndarray', eta: 'np.ndarray', *, kbt: 'float', m: 'float') -> 'np.ndarray'
independent reference: rho_max = 0.9 version (eta_max = 0.7069)

| eta | Z paper | Z module | rel diff | Z' paper (analytic) | Z' complex-step | Z' module | rel diff | c_s paper | c_s module | rel diff | tests_20260913.kr_cs | rel diff vs paper |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.05 | 1.10838774198985 | 1.10838774198985 | 4.0e-16 | 2.3476177330678 | 2.3476177330678 | 2.3476177330678 | 0.0e+00 | 1.56661801829183 | 1.56661801829183 | 2.8e-16 | 1.56661801830068 | 5.6e-12 |
| 0.30 | 2.06326087957606 | 2.06326087957605 | 2.2e-16 | 6.03583334839163 | 6.03583334839163 | 6.03583334839163 | 1.5e-16 | 2.8515007173912 | 2.85150071739119 | 3.1e-16 | 2.85150071752567 | 4.7e-11 |
| 0.50 | 4.10636394336527 | 4.10636394336527 | 0.0e+00 | 16.7336836142638 | 16.7336836142638 | 16.7336836142638 | 0.0e+00 | 5.41621921508607 | 5.41621921508607 | 0.0e+00 | 5.41621921569349 | 1.1e-10 |
| 0.65 | 8.40804352772069 | 8.4080435277207 | 4.2e-16 | 45.9449284500896 | 45.9449284500896 | 45.9449284500896 | 0.0e+00 | 10.4387471941955 | 10.4387471941955 | 3.4e-16 | 10.4387471910944 | 3.0e-10 |

largest relative difference module vs paper (Z, Z', c_s; and analytic vs complex-step Z'): 4.6e-16
tests_20260913.kr_cs vs paper: largest relative difference 3.0e-10 (above 1e-10; reported, not part of the verdict -- see how kr_cs forms Z')

for information (NOT the verdict) -- the other two published versions against the module, largest relative difference over the four eta:
| version | Z | Z' (analytic) | c_s |
|---|---|---|---|
| rho_max = 0.88 | 1.4e-06 | 7.1e-05 | 1.0e-05 |
| rho_max = 0.89 | 6.8e-07 | 7.2e-06 | 1.5e-06 |
module comment, plot_speed_of_sound_edmd.py:620: '# Kolafa & Rottner (2006), rho_max=0.90 fit.  Their x is eta/(1-eta).'

**VERDICT: PASS** (criterion: relative difference <= 1e-10)

**P2–P4: printed by `python3 hspist3/validation/paper1_draft_audit_20261014.py --p` before any edit** (verbatim). It contains:
- the P3 table, i.e. the four Z′ values of the melting caption with the η at which each is evaluated and whether that η lies inside the fit range;
- the old → new list for every text change.

fit range: eta <= pi*0.90/4 = 0.706858 -> '0.7069' (KR2006_ETA_MAX); comparison cutoff KR2006_PLOT_ETA_MAX = 0.69

### P3: the Z' values quoted in the melting caption

| caption value | eta | Z' analytic (module) | Z' finite difference (h = 1e-5, as the melting script) | reproduced | eta <= fit range? | action |
|---|---|---|---|---|---|---|
| 48 | 0.670 | 48.2615 | 48.2615 | yes | yes | keep, labelled fluid-branch fit |
| 28 | 0.690 | 27.7824 | 27.7824 | yes | yes | keep, labelled fluid-branch fit |
| -9 | 0.700 | -8.6854 | -8.6854 | yes | yes | keep, labelled fluid-branch fit |
| +3386 | 0.720 | 3392.4452 | 3392.4505 | **NO** | **NO** (beyond 0.7069) | REMOVE |

KR's own c_s maximum on 0.60 <= eta <= 0.7069 (analytic Z'): eta = 0.6834 -> caption says 'near eta = 0.683': consistent

### P2 / P3 / P4: every text change, old -> new (printed BEFORE editing)

| kind | where | tex line | old | new | why |
|---|---|---|---|---|---|
| TEX | Methods (sec. Model and apparatus): new paragraph after 'Integer-pixel box width' (P2: fit named once) | 123 | `%   sec. 14.3 (A2 rule 0ddefa9, results 14.3.1); validation/paper1_A2_boxtrunc_261002.py` | `%   sec. 14.3 (A2 rule 0ddefa9, results 14.3.1); validation/paper1_A2_boxtrunc_261002.py  % ##CHRIS 2026-10-02: the reference equation of state, named and bounded (plan author's decision, methods sec. 15.1) \paragraph{Reference equation of state} We compare with the hard-disk equation of state of Kolafa and Rottner~\cite{kolafa2006}: the $\rho_{\max} = 0.90$ fit of their Eq.~(7), with the coefficients of their \S3.2, which is fitted to $\eta \le 0.7069$. $Z'$ is taken from the same fit, and data are compared with it only for $\eta \le 0.69$. % TODO-source: hspist3/plot_speed_of_sound_edmd.py:620-640 (KR2006_COEFFICIENTS, KR2006_ETA_MAX); validation/paper1_kr_sanity_261002.py (PASS)` | fit range 0.7069 = pi*0.90/4; comparison cutoff 0.69 = KR2006_PLOT_ETA_MAX |
| TEX | Fig. csvseta caption (P2) | 255 | `beyond $\eta \approx 0.69$ lie in that region, where no fluid reference exists; they are shown but are not counted as a deviation from Kolafa--Rottner.}` | `beyond $\eta = 0.69$ lie in or near that region; Kolafa--Rottner ($\rho_{\max} = 0.90$ fit, fitted to $\eta \le 0.7069$) is compared with data only for $\eta \le 0.69$, so those points are shown but are not counted as a deviation from it.}` | names the range exactly; 'no fluid reference' was imprecise below 0.7069 |
| TEX | sec:ordering, first sentence (P2) | 344 | `Kolafa--Rottner is a \emph{fluid} equation of state fitted below the transition, so a deviation from it above $\eta \approx 0.69$ measures the extrapolation, not the gas.` | `Kolafa--Rottner is a \emph{fluid} equation of state; the $\rho_{\max} = 0.90$ fit used here is fitted to $\eta \le 0.7069$ and compared with data only for $\eta \le 0.69$, so above that a deviation from it measures the fluid-branch fit, not the gas.` | the fit reaches 0.7069, into the coexistence interval; 'fitted below the transition' was imprecise |
| TEX | melting caption (P3) | 383 | `Kolafa--Rottner is solid inside its fitted range and dashed beyond it, and is cut off at $\eta = 0.705$ deliberately: its fitted $Z'$ runs $48$ at $\eta = 0.67$, $28$ at $0.69$, $-9$ at $0.700$ and $+3386$ at $0.720$, so a $c_s$ built from it there is meaningless rather than merely uncertain. The negative $Z'$ already inside the stated fitted range is also why KR's own $c_s$ turns over near $\eta = 0.683$.}` | `Kolafa--Rottner ($\rho_{\max} = 0.90$ fit, fitted to $\eta \le 0.7069$) is solid where it is compared with data, $\eta \le 0.69$, and dashed from there to $\eta = 0.705$, inside its fit range; the fluid-branch fit gives $Z' = 48$ at $\eta = 0.67$, $28$ at $0.69$ and $-9$ at $0.700$. Beyond its fit range the fluid-branch fit is an extrapolation and its derivative carries no information. The negative $Z'$ already inside the fit range is also why KR's own $c_s$ turns over near $\eta = 0.683$.}` | +3386 is at 0.720 > 0.7069: removed; the three values inside the range stay |
| TEX | zoom caption (P4) | 267 | `larger systems agree with it within their errors, except` | `larger systems agree with it within $2\sigma$, except` | plan author's decision 2 |
| UNCHANGED | sec:ordering title, l. 342 | 342 | `\subsection{Above $\eta \approx 0.69$: the right comparison is not Kolafa--Rottner}` | `(unchanged)` | states where the comparison stops, consistent with the decision |
| UNCHANGED | l. 398 | 398 | `$c_s(\eta)$ above $0.69$ is the melting transition` | `(unchanged)` | about our data, not the KR range |
| UNCHANGED | l. 404 (Outlook) | 404 | `Above $\eta \approx 0.69$ the comparison is with` | `(unchanged)` | where the comparison stops, consistent |
| SCRIPT | paper1_canonical_20260919.py:74 legend, main figure (first KR legend: fit named) | - | `label="Kolafa-Rottner 2006 (valid to eta = 0.69)")` | `label="Kolafa-Rottner 2006, rho_max = 0.90 fit (Eq. 7), fitted to eta <= 0.7069;\ncompared with data for eta <= 0.69")` | the first figure that draws the KR equation of state |
| SCRIPT | overlay_N100_vs_A2_20260915.py:68 legend | - | `label="Kolafa–Rottner 2006 (valid to η ≈ 0.69)")` | `label="Kolafa–Rottner 2006 (fitted to η ≤ 0.7069; compared for η ≤ 0.69)")` | range bounded |
| SCRIPT | paper1_melting_figure_20261002.py:68 legend | - | `label="Kolafa–Rottner 2006 (fitted range, $\\eta \\leq 0.69$)")` | `label="Kolafa–Rottner 2006, compared with data ($\\eta \\leq 0.69$)")` | 0.69 is the comparison cutoff, not the fit range (0.7069) |
| SCRIPT | paper1_melting_figure_20261002.py:70 legend | - | `label="Kolafa–Rottner, EXTRAPOLATED (no fluid branch here)")` | `label="same fit, not compared ($0.69 < \\eta \\leq 0.705$; fit range $\\eta \\leq 0.7069$)")` | 0.69-0.705 is inside the fit range: not an extrapolation of the fit |
| SCRIPT | paper1_melting_figure_20261002.py annotation | - | `"KR's fitted $Z'$ changes sign near $\\eta = 0.70$\nand diverges by $0.72$: $c_s$ from it is\n"                 "meaningless here, not merely uncertain"` | `"KR's fluid-branch fit: $Z'$ changes sign near $\\eta = 0.70$,\ninside its fit range; beyond $\\eta = 0.7069$ the fit is\n"                 "an extrapolation and its derivative carries no information"` | no statement about the fit beyond its range |
| DOC | writeup/papers/README.md, KolafaRottner2006 row | - | `| 10.1080/00268970600880574 | the reference EOS, valid to η ≈ 0.69 |` | `| 10.1080/00268970600967963 | the reference EOS: the ρ_max = 0.90 fit (Eq. 7, coefficients of §3.2), fitted to η ≤ 0.7069; compared with data for η ≤ 0.69 |` | range named; DOI corrected from the paper's own title page [SOURCE: PDF in ZZZ_PAPER/] |

`--p --check-scripts` printed:

    paper1_canonical_20260919.py:74 legend, main figure (first KR legend: fit named): old ABSENT, new present
    overlay_N100_vs_A2_20260915.py:68 legend: old ABSENT, new present
    paper1_melting_figure_20261002.py:68 legend: old ABSENT, new present
    paper1_melting_figure_20261002.py:70 legend: old ABSENT, new present
    paper1_melting_figure_20261002.py annotation: old ABSENT, new present
    writeup/papers/README.md, KolafaRottner2006 row: old ABSENT, new present

**Applied.**
- `--p --apply --apply-scripts` applied the TEX, SCRIPT and DOC rows from the same strings as printed.
- `--p --check-scripts` then reports every old string absent and every new string present.

**P3 in one line (DATA).**
- Of the four Z′ values, three are evaluated inside the fit range and stay, labelled "fluid-branch fit": 48 at η = 0.67, 28 at 0.69, −9 at 0.700.
- One is beyond it and is removed: +3386 at η = 0.720 > 0.7069. The caption's +3386 did not even reproduce; the module gives 3392.4 there.
- "KR's own $c_s$ turns over near η = 0.683" is confirmed (0.6834).

**Disclosed, beyond the printed list.**
- **Melting-figure layout (no data change).** The y axis started at 11.5, which clipped most of the drawn KR segment ($c_s$ 10.08–11.70) and hid the dashed 0.69–0.705 part that the legend and caption describe. The lower limit now follows the drawn segment (9.5, printed by the script). The note about the fit, which sat below the axes over the tick labels in the old figure too, moved inside the axes.
- **`paper1_melting_figure_20261002.py`.** Its docstring and one internal comment carried the old "extrapolation beyond 0.69" wording and the +3386. Each got a dated note; the old text was not erased.
- **README DOI.** The references README row for Kolafa–Rottner gave DOI `10.1080/00268970600880574`. The paper's own title page says `10.1080/00268970600967963` [SOURCE: PDF in ZZZ_PAPER/], so it was corrected in the same row.
- **"The first KR figure legend"** is read as the first figure that draws the KR equation of state, i.e. the main $c_s(\eta)$ figure. The earlier ladder-line figure shows KR only as one value at η = 0.1123.

**P5: figures with changed legends.** Main $c_s(\eta)$, N100-vs-A2 overlay, and melting region.
- Dated copies `_pre_krfit_261002` were made first (`cp -n`, `cmp`; 6 files).
- Regeneration, each by its original script:
  - the main figure via `paper1_canonical_20260919.draw_main()` only, so no other figure was touched;
  - the overlay with the environment `paper1_canonical_20260919.py` gives it; its printed table is identical to the § 14.5 run;
  - the melting figure; its printed dip is unchanged (2.827, 9.6σ plotted, 18.4σ propagated), plus the new line `KR drawn: c_s 10.08 to 11.70 on 0.66 <= eta <= 0.705 -> y axis from 9.5`.
- **Draft pdf rebuilt:** 6 pages, 640839 bytes, 0 undefined references, 6 warnings (all pre-existing), 1 overfull hbox.

## 16. Estimator note: a light-mass term at π/8, and what it does to the unweighted slope (2026-10-07; OPEN; no registered result changes)

Requested by the plan author's third decision of 2026-10-07, item 5 (261012 § 4.4.13). **Nothing below changes a registered estimator or a registered result**; it records what the best plain-fluid data so far say about them, and the choices the plan author has to make.

**Plain summary.**
- **The deficit.** At η = π/8 the lightest divider (M = 50, α = 0.5) gives an implied sound speed 0.56–0.69 % below the other masses, 4.5–5.5 standard errors. This holds on both rescheduling paths alike, so it is not an engine effect.
- **The registered unweighted slope** gives the two lightest masses 60 % of the weight. It therefore sits 0.25–0.30 % below the weighted slope.
- **Not the gas inertia.** The registered model already contains it, exactly: the root of cot K = αK. The decision suggested the effective-mass correction M + 2N_s m/3 as the missing term; it is not missing.
- **Damping explains part of the deficit.** The spectral peak of a damped mode sits below its undamped frequency by 1/(4Q²): 0.29 % at M = 50. After removing it, −0.42 ± 0.08 % remains at M = 50. That remainder is OPEN.

**Data [DATA].**
- **Source and cell.** Test T (261012 §§ 4.4.12–4.4.13), cell epi8_H_H10_L10 (η = π/8, H = L_0 = 10, N_s = 50). There are 400 seeds per policy at M = 50 and 100 at each other mass. The quantity is the implied sound speed c_s,M = ν_M / x_M, with the registered x_M and the argmax ν.
- **M = 50:** 3.79170 ± 0.00478 (minimal path) and 3.78592 ± 0.00481 (legacy path). The weighted single-c_s slopes are 3.81307 ± 0.00115 and 3.81224 ± 0.00116. That gives **−0.56 % (−4.5 SE) and −0.69 % (−5.5 SE)**.
- **M = 100:** 3.80292 ± 0.00797 and 3.80471 ± 0.00794, i.e. −0.27 % and −0.20 % (−1.3 and −1.0 SE).
- **M ≥ 200:** within about 2 SE of the weighted slope, with one exception: Test T legacy at M = 1500 sits at +0.24 % (+3.5 SE) with the argmax estimator, and +2.3 SE with the refined ν_d.
- **Single-c_s fit.** The χ² at the weighted slope (8 dof) is 31.9 for minimal and 46.6 for legacy, both p < 0.001.

**The unweighted slope and the light masses [DERIVATION].**
- **The weights.** The registered c_s = Σ x_M ν_M / Σ x_M² = Σ_M w_M c_s,M, with w_M = x_M² / Σ x². For this cell w = 36.8 % at M = 50, 23.5 % at 100 and 13.5 % at 200, falling to 1.6 % at 2000 (table, column w_M). The two lightest masses carry 60 %.
- **Effect at Test T precision.** Unweighted minus weighted is −0.245 % (minimal) and −0.302 % (legacy).
- **Effect at 25 seeds.** The sign follows the M = 50 fluctuation: the campaign anchor gives +0.248 %, and its M = 50 sat 2.75 SE above Test T's; the replay gives −0.180 % (261012 § 4.4.13 item 4c).

**The gas inertia is already in the registered model [SOURCE].**
- **Decision 3's hypothesis.** It named the gas-inertia correction M + 2N_s m/3 (+67 % at M = 50, +1.7 % at M = 2000) as the likely missing light-mass term.
- **What the registered estimators use** is the exact root:
  - `validation/paper1_confinement_results_261004.py:8`: "x = K(alpha)/(2 pi L_eff,true), cot K = alpha K). alpha = M/(2 N_s m)".
  - `:202`: `rows.append(dict(M=M, alpha=al, K=T.k_root(al), ...`.
  - `:205`: `x = np.array([r["K"] / (2 * math.pi * LeT) for r in rows])`.
  - `validation/tests_20260913.py:144-149`, `k_root`: bisection on `math.cos(mid) / math.sin(mid) - alpha * mid`.
  - The same holds for the canonical A1 v2 table: `tests_20260913.py:182-183` (`x_of` = `k_root(M / (2.0 * N_SIDE)) / (2 * math.pi * l_eff(L0, wall_t))`), used at `paper1_populate_cs_err_20261002.py:136` and `:158`.
- **What cot K = αK contains.** It is the eigenvalue condition of a piston between two gas columns in linear acoustics, and it contains the gas inertia to all orders.
- **Where M + 2N_s m/3 comes from.** It is the α ≫ 1 limit of that root, K ≈ (α + 1/3)^(−1/2). It enters only the k_S^dyn identity, at α ≥ 5 (`paper1_confinement_results_261004.py:300`, `:303`: `Mh = r["M"] + 2.0 * c["Ns"] / 3.0`), where it is within 0.04 % of the exact root (table).
- **At α = 0.5** the effective-mass root would be 1.725 % too high. The registered c_s does not use it.
- So the light-mass deficit is a term **beyond** the continuum gas inertia. [INFERENCE] The project's own eigenvalue equation already covers α = 0.5. The effective-mass form is its heavy-divider limit, and taking it as the missing term would apply a formula outside its regime.

**Table**, printed by `cd hspist3 && python3 validation/resched_testT_followup_261007.py` (the item-5 part), verbatim:

```
## Item 5 -- the gas inertia in the registered model, and the light-mass deficit of the plain-fluid baseline (Test T legacy)

| M | alpha = M/(2 N_s m) | gas inertia (2/3) N_s m / M [%] | K, cot K = alpha K (registered) | K_eff = (alpha + 1/3)^(-1/2) | K_eff/K - 1 [%] | w_M [%] | Q, campaign (sec. 13 model) | Q, Test T legacy (4a model, averaged ACF) | -1/(4Q^2) [%] (Test T Q) | deficit, argmax [%] | deficit, nu_d [%] | deficit, nu_0 [%] |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 50 | 0.5 | 66.7 | 1.07687 | 1.09545 | +1.725 | 36.8 | 9.6 | 9.3 | -0.288 | -0.690 +- 0.126 | -0.549 +- 0.078 | -0.424 +- 0.078 |
| 100 | 1 | 33.3 | 0.86033 | 0.86603 | +0.662 | 23.5 | 15.2 | 13.5 | -0.138 | -0.197 +- 0.208 | +0.002 +- 0.161 | +0.052 +- 0.160 |
| 200 | 2 | 16.7 | 0.65327 | 0.65465 | +0.212 | 13.5 | 21.6 | 20.3 | -0.061 | +0.051 +- 0.123 | +0.157 +- 0.100 | +0.169 +- 0.100 |
| 300 | 3 | 11.1 | 0.54716 | 0.54772 | +0.103 | 9.5 | 26.6 | 25.5 | -0.038 | -0.183 +- 0.124 | -0.109 +- 0.077 | -0.109 +- 0.077 |
| 500 | 5 | 6.7 | 0.43284 | 0.43301 | +0.040 | 5.9 | 30.0 | 32.2 | -0.024 | +0.093 +- 0.106 | +0.104 +- 0.073 | +0.097 +- 0.073 |
| 750 | 7.5 | 4.4 | 0.35723 | 0.35729 | +0.018 | 4.0 | 34.2 | 40.3 | -0.015 | -0.033 +- 0.088 | -0.014 +- 0.056 | -0.025 +- 0.056 |
| 1000 | 10 | 3.3 | 0.31105 | 0.31109 | +0.010 | 3.1 | 52.9 | 43.7 | -0.013 | -0.033 +- 0.079 | +0.059 +- 0.050 | +0.047 +- 0.050 |
| 1500 | 15 | 2.2 | 0.25536 | 0.25538 | +0.005 | 2.1 | 46.9 | 54.9 | -0.008 | +0.239 +- 0.069 | +0.114 +- 0.049 | +0.100 +- 0.049 |
| 2000 | 20 | 1.7 | 0.22176 | 0.22177 | +0.003 | 1.6 | 69.1 | 64.9 | -0.006 | +0.028 +- 0.061 | +0.014 +- 0.043 | -0.002 +- 0.043 |

weighted slopes (Test T legacy): argmax 3.81224, nu_d 3.81215, nu_0 3.81286; deficit = c_s,M / (that slope) - 1, SE from the seeds; nu_0 = sqrt(nu_d^2 + (1/(2 pi tau_r))^2) with tau_r from the fit of the mass's seed-averaged ACF
single-c_s chi2 at the weighted slope, argmax: 46.6 (8 dof, p 1.84e-07)
single-c_s chi2 at the weighted slope, nu_d: 63.2 (8 dof, p 1.08e-10)
single-c_s chi2 at the weighted slope, nu_0: 41.7 (8 dof, p 1.57e-06)
```

**Damping accounts for part of the deficit [DERIVATION; INFERENCE where marked].**
- **What the estimator measures.** The registered ν is the periodogram peak (argmax). The continuum model predicts the undamped frequency ν_0.
- **Two shifts below ν_0** for a mode driven by thermal noise and damped at rate 1/τ_r, with Q = π ν τ_r (methods § 13):
  - the displacement spectrum peaks at ν_0 √(1 − 1/(2Q²)) ≈ ν_0 (1 − 1/(4Q²));
  - the ACF oscillates at ν_d = ν_0 √(1 − 1/(4Q²)) ≈ ν_0 (1 − 1/(8Q²)).
- **Size.** Q is 9.3 at M = 50: from the fit of the seed-averaged ACF of Test T legacy; the campaign's § 13 table gives 9.6. Q rises to 65 at M = 2000. So the argmax peak sits **0.29 % low at M = 50**, 0.14 % at M = 100, 0.06 % at M = 200, and less above.
- **Measured at M = 50 [DATA, table].** Each against its own weighted slope:
  - argmax: −0.690 ± 0.126 %;
  - ν_d (damped cosine with free phase, fitted per trajectory): −0.549 ± 0.078 %;
  - ν_0 = √(ν_d² + (1/(2π τ_r))²): −0.424 ± 0.078 %.
- [INFERENCE] The two steps (0.14 % and 0.13 %) are the predicted 1/(8Q²) = 0.14 % each.
- **M = 100.** The deficit disappears after the correction: ν_0 gives +0.05 ± 0.16 %.
- **M = 50 keeps −0.42 ± 0.08 % (5.4 SE).** The single-c_s χ² with ν_0 is 41.7 (8 dof).

**OPEN.**
- **The −0.42 % at α = 0.5 is unexplained.** It is not the gas inertia (in the model), not the damping (removed by ν_0), and not the rescheduling path (both paths alike).
- [INFERENCE] Candidates, none tested yet:
  - the length that carries the gas mass in α (the continuum assumes the N_s m spread over L_eff, and only small α feels a different inertial length);
  - the light divider's own thermal speed, √(kT/M) = 0.14 of the gas's at M = 50 against 0.02 at M = 2000;
  - a frequency dependence of the effective sound speed at the mode's wavelength, 2π/K ≈ 5.8 L_eff at α = 0.5.
- **Also OPEN:** whether the melting-window "masses disagree" hint (χ²_red 1.8–14, unweighted) contains the same light-mass term. That needs the window cells' per-mass c_s,M against a plain-fluid baseline.

**Consequences. The decision is the plan author's; nothing is changed here [INFERENCE].**
1. **Paper 1 reports both estimators,** the registered unweighted slope and the weighted (minimum-χ²) slope, with their difference as a systematic. At Test T precision it is 0.25–0.30 %, unweighted lower. At 25 seeds per mass its sign is set by the M = 50 fluctuation.
2. **Future campaigns pre-register a better estimator.** The options:
   - the weighted slope;
   - a light-mass term (or no α = 0.5 in the c_s fit);
   - the frequency from the ACF fit instead of the periodogram peak. ν_d has 1.2–1.6× smaller per-seed scatter, needs no bins, and with τ_r gives ν_0, which removes the damping shift.
3. **The melting pre-registration carries a plain-fluid baseline cell** outside the window, with the same masses. The "masses disagree" hint is judged against it.
4. **The argmax peak's damping shift −1/(4Q²) can be computed** from the Q that § 13 already measures, and can be stated per cell.
