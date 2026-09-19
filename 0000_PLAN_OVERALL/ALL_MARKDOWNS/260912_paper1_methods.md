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
