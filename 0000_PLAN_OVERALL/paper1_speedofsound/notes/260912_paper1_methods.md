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
