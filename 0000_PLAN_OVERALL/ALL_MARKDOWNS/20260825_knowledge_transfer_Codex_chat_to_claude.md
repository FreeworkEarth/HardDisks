Yep — for the current r11_c25 run, use this. It will analyze the 384 published valid low-eta traces and skip the one preserved invalid trajectory, just like before.

cd /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks

RUN="hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/simulation_oscillation_design_pilot_r11_c25"

MPLCONFIGDIR=/tmp/mplcfg_speed_r11_20260822 \
/opt/homebrew/bin/python3 hspist3/analyze_speed_of_sound_by_eta.py \
  --dir "$RUN/raw_simulations/01_high_eta" \
  --dir "$RUN/raw_simulations/02_mid_eta" \
  --dir "$RUN/raw_simulations/03_low_eta" \
  --out "$RUN/replot_exploratory_r11_384of385" \
  --window-spec "15:0.002:0.20,50:0.0003:0.08,100:0.00015:0.03,200:0.00008:0.015,inf:0.00004:0.008" \
  --peak-signal vel \
  --min-measured-oscillations 25 \
  --max-relative-bin-width 0.04 \
  --min-peak-snr 5 \
  --min-group-acceptance 0.8 \
  --min-group-repeats 3 \
  --min-fit-r2 0.8 \
  --write-final \
  --roman-ref \
  --theory-cs simple \
  --per-eta-plots

Important: do not add

--require-valid-ledger
--strict-counts

for this exploratory plot, because those will stop analysis at the known failed seed.

The combined output should end up under:

hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/
simulation_oscillation_design_pilot_r11_c25/
replot_exploratory_r11_384of385/

and the main plot should be approximately:

replot_exploratory_r11_384of385/final_plots/FINAL speed_of_sound_on_packing_fracture.pdf

One small naming detail: 384/385 refers specifically to the low-eta batch; the high- and mid-eta folders have their own additional trajectories. So the output name is fine as shorthand for “one known invalid trajectory excluded,” but it is not the total full-campaign trajectory count.

ok this looks pretty decent right ? 


now we have to find out the error right ? and maybe use another EOS fir over 0.7?

Yes — this is a genuinely strong-looking result now. The low/mid-density part is especially convincing: your simulation follows Román where the datasets overlap, approaches the ideal-gas limit smoothly on the left, and then shows the expected rapid stiffening at high packing fraction.

But the plot is also telling us exactly what to investigate next.
1. First priority: the reproducible EDMD failure

Yes, this is now the error to find. The same exact case failed twice:
L0​=150,M=50,seed=2381038820

with the exact same particles 73/74, overlap
0.0925995654

and relative normal velocity
−0.834756336.

That's excellent from a debugging perspective because it means the failure is deterministic. I would make that seed the next regression test and inspect the event history immediately before the overlap: scheduled particle-particle collision time, last event touching particles 73/74, event invalidation/generation counters, positions/velocities before and after that event. We should not call the core fully validated until we understand that one.
2. The high-density part should not just extend Kolafa–Rottner through η>0.7

Your current plot actually handles this reasonably: the red Kolafa–Rottner curve is explicitly labelled fluid only and stops before the transition.

Large-system hard-disk simulations place the first-order liquid–hexatic coexistence region roughly at
0.700≲η≲0.716​

with hexatic around η∼0.718 and the continuous hexatic–solid transition around η∼0.720.

And your high-density points are almost perfectly placed to probe this:
L0​6.546.045.615.455.30​η0.60050.65020.70000.72060.7409​​

So you're effectively sampling:

η≈0.60      η≈0.65      η≈0.700      η≈0.721      η≈0.741
   │            │            │             │             │
dense       pre-transition  coexistence   transition     solid-like
fluid                        onset         region

That's potentially much more interesting than merely extending the original Román curve.
3. Yes, add a newer phase-aware EOS — but carefully

Liu's global hard-disk EOS is specifically designed to cover stable liquid, liquid–hexatic coexistence, hexatic, and then join onto a solid branch. It was published in 2021.

But there is now even newer work: Mier-y-Terán published a 2024 analysis of EOS behavior around the fluid–hexatic transition using 105-disk MD data, and a 2025 follow-up gives analytical EOS/Helmholtz expressions including the hexatic transition.

So for a paper I'd probably use:

    Henderson/SPT: historical Román comparison.

    Kolafa–Rottner: high-accuracy ordinary isotropic-fluid reference.

    Liu 2021: global phase-aware EOS.

    Potentially Mier-y-Terán 2025: modern transition-focused comparison.

But here's the critical part:

    I would not simply convert the solid branch of Liu into the same fluid cs​(η) formula and draw it as though nothing changes at η=0.7.

For a fluid,
cs2​=m1​(∂n∂p​)s​,

and because hard disks have purely kinetic internal energy, for
p=nkB​TZ(η)

you get
cs2​=mkB​T​[Z+ηdηdZ​+Z2].​

That is a thermodynamic fluid sound speed.

Once you have a true solid, shear rigidity exists. The longitudinal collective mode is then governed by elastic constants as well as compressibility. In isotropic elasticity the longitudinal and transverse modes involve Lamé/shear moduli rather than pressure EOS alone.

So beyond the transition, your blue points may still represent a perfectly real longitudinal collective mode — but calling them the same simple “fluid sound speed” becomes questionable.

And that's actually more interesting scientifically.
4. The η=0.72 point is not automatically an error

Your blue point around
η≃0.72,c∼16.5

looks dramatic, but hard disks become extremely incompressible as packing rises. A rapidly increasing longitudinal propagation speed is qualitatively sensible.

What we don't yet know is whether that point represents:

    the same hydrodynamic piston–gas acoustic mode;

    a strongly ordered/hexatic collective mode;

    a solid-like longitudinal elastic mode;

    a finite-size artifact;

    or an FFT/eigenmode-fitting artifact.

Your system is only 100 particles, while the canonical phase-boundary studies used enormous systems precisely because the correlation lengths near melting become very large. So the literature's η=0.700,0.716,0.720 boundaries should be shown as thermodynamic-limit reference regions, not treated as proof that your 100-disk box is actually hexatic/solid at exactly those values.

That's important for a paper.
What I think the final figure could become

I would evolve this toward something like:

                 measured piston mode c(η)
                           ● ●
                       ●
                   ●
               ●
           ●
       ●
    ●
 ●
────────────────────────────────────────────── η

│ ideal/dilute │ isotropic fluid │ coexistence │ hexatic │ solid │
                 ↑
        Henderson / KR / Román

                                 ↑
                         Liu / modern EOS
                         pressure/compressibility
                         comparison

                                                   ↑
                                             elastic-mode
                                             interpretation

In other words, don't force one theory curve across everything. Let the change of physical regime be part of the result.

And I particularly like one thing about your current figure: the blue curve starts separating upward from the ordinary fluid EOS already around η∼0.6−0.65. Recent EOS work specifically discusses orientational-ordering effects appearing in the dense fluid as the transition is approached, so that region deserves closer study rather than immediately assuming your simulation is wrong.

So I would now do this order:

debug seed 2381038820 → verify the high-η FFT/mass fits individually → add Liu/2025 phase-aware EOS as a separate comparison → measure an order parameter such as ψ6​ at the high-η points → only then decide how to label the η>0.7 blue measurements.

If ψ6​, the FFT spectra, and the mass-dependent eigenmode fits all change character right where your blue curve shoots upward, that could become far more interesting than simply reproducing Román.

Yes. I would now separate the theory curves by their physically justified regime instead of letting every EOS run across the whole plot. And yes, we can debug this failed seed without rerunning hundreds of trajectories.
1. EOS strategy for the final figure

For the sound-speed figure, I would currently use this logic:
Packing fraction	Physical regime	What I would plot
η→0	ideal/dilute limit	cs​→2kB​T/m
​
0<η≲0.69−0.70	isotropic hard-disk fluid	SPT, Henderson, Kolafa–Rottner
η∼0.700−0.716	liquid–hexatic coexistence	do not blindly continue ordinary-fluid curves
η∼0.716−0.720	hexatic	phase-aware EOS comparison, cautiously
η≳0.720	solid	don't use ordinary fluid sound-speed formula; need longitudinal elastic-mode treatment

That is consistent with the large-system hard-disk literature: the liquid–hexatic transition is first-order, followed by a narrow hexatic regime and then the continuous transition to solid.

So yes: I would stop/fade SPT and Henderson around η≃0.70 just like Kolafa–Rottner. It doesn't mean the formulas mathematically cease to return numbers. It means those numbers no longer represent the equilibrium phase you are simulating.

Your current figure will become more scientifically honest if it says visually:

                           ordinary isotropic fluid
                     <-------------------------------->

ideal       dilute                 dense                    phase change
 |            |                      |                          |
 v            v                      v                          v
──────────────────────────────────────────────────────────────────── η
0           .1          .3          .5      .65      .700 .716 .720
                                                         │     │
SPT ----------------------------------------------------X
Henderson ----------------------------------------------X
Kolafa-Rottner -----------------------------------------X

                                                        coexist.
                                                           │
                                                           ▼
                                                phase-aware EOS

                                                                 solid
                                                                   │
                                                                   ▼
                                                       elastic treatment

I would not connect the theory curves across the transition.
2. Is Liu the correct EOS after η≈0.70?

Liu 2021 is a sensible reference, but I would not make it the only or automatically “correct” replacement.

Liu specifically constructed a global hard-disk EOS designed to cover stable liquid, liquid–hexatic transition, hexatic phase, and then join onto a solid EOS. So it solves exactly the problem that Henderson/SPT/KR don't solve: describing the pressure across multiple phases.

But since then there has been more work. Mier-y-Terán 2024 compared EOS approximations against equilibrium MD of 105 hard disks near the fluid–hexatic transition, and the 2025 paper constructs analytical pressure and free-energy expressions specifically including the hexatic transition and orientational contribution.

So for the paper I would probably treat the references as:

Román / SPT / Henderson
        │
        └── historical validation

Kolafa–Rottner
        │
        └── high-accuracy isotropic-fluid benchmark

Liu 2021
        │
        └── global phase-aware EOS

Mier-y-Terán 2025
        │
        └── newer transition/hexatic-focused EOS

That's stronger than picking one formula as gospel.
3. One major warning: EOS ≠ automatically sound-speed curve

This is the part we should be very careful with.

For your ordinary isotropic hard-disk fluid we use
p=nkB​TZ(η)

and derive the adiabatic sound speed.

For hard disks with the assumptions we're using, that gave the relation
cs2​=mkB​T​[Z+ηdηdZ​+Z2].

Fine in the stable fluid.

But a first-order coexistence region is fundamentally different. There isn't simply one homogeneous phase whose pressure derivative you differentiate through the whole coexistence interval.

And once you're in the solid, the longitudinal mode isn't determined by the scalar pressure EOS alone. Elastic moduli matter.

Therefore I would do:

MAIN c_s PLOT

η < ~0.70:
    EOS-derived sound-speed curves

η ~0.70–0.716:
    shade coexistence
    do NOT force ordinary-fluid c_s formula through it

η ~0.716–0.720:
    show measured mode
    optionally phase-aware thermodynamic prediction,
    clearly labelled

η > ~0.720:
    measured longitudinal piston mode
    NOT automatically called fluid-EOS sound speed

And separately I would love a pressure/EOS panel.

There Liu 2021 or the 2025 EOS can be shown all the way across the transition because that paper actually gives an EOS. That's the cleaner place for it.
4. Your blue points above 0.7 could become especially interesting

Your high-density grid gives approximately:
L0​=5.61⇒η≃0.7000 L0​=5.45⇒η≃0.7206 L0​=5.30⇒η≃0.7409.

That's almost absurdly nice sampling:

     liquid        coexistence       hexatic       solid
        │                │              │             │
────────┼────────────────┼──────────────┼─────────────┼── η
       .65             .700           .716          .720
                         ●                ●             ●
                    our ~.700       our ~.721       our ~.741

But because you only have 100 disks, those bulk phase boundaries cannot prove that your individual finite system is actually hexatic/solid.

That is why ψ6​ becomes useful later.

If your high-η points simultaneously show:

    EOS departure,

    changing cs​,

    strong ψ6​,

    different structure factor,

    and a clean change in the piston spectral mode,

then it starts becoming a really interesting phase-regime story.
5. Now the bug: YES, the same seed has reproduced exactly

This is now extremely clear.

Twice you got:

L0      = 150
M       = 50
repeat  = 3
seed    = 2381038820

particles = 74,73
overlap   = 0.0925995654 px
tolerance = 0.000024 px
v_rel,n   = -0.834756336

This is basically the best debugging scenario you can ask for.

The probability that independent numerical noise happened to produce the same:

    seed,

    pair,

    overlap to 10 decimal places,

    and velocity

is negligible.

We have a deterministic reproduction.
6. Do we already have the bad trajectory?

Probably partially.

Your validation output previously said:

found 349 published traces
found 1 preserved invalid traces
failure ledger contains 1 failure rows

so the workflow appears to deliberately preserve failed trajectory output rather than delete it.

First thing I would do now — without changing anything — is locate everything associated with the seed.

For the new r11 run:

cd /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks

RUN="hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/simulation_oscillation_design_pilot_r11_c25"

rg -n \
  "2381038820|particle_particle_overlap|0.0925995654|Particles 74 and 73" \
  "$RUN/raw_simulations/03_low_eta"

Then:

find "$RUN/raw_simulations/03_low_eta" -type f | \
  grep -E '2381038820|invalid|failure|ledger'

That should tell us exactly what the workflow preserved.
7. But the normal trajectory probably isn't enough to diagnose the EDMD bug

This distinction matters.

Your ordinary speed-of-sound trace likely stores things such as
t,xw​(t),vw​(t)

and possibly aggregate quantities.

That is enough for FFT.

But the thing we need for this bug is the event-level history of particles 73 and 74.

Specifically, immediately before the failure we want:

current simulation time

particle 73:
    x, y
    vx, vy
    collision generation/version

particle 74:
    x, y
    vx, vy
    collision generation/version

distance(73,74)
surface gap
relative normal velocity

NEXT SCHEDULED EVENTS:
    PP(73,74) ?
    PP(73,something)?
    PP(74,something)?
    wall event?
    spring-wall event?

last ~20 processed EDMD events

The question is:

    Why was the event at which 73 and 74 first touched not processed?

That's the actual bug.
8. What I suspect happened

At physical contact,
∣r74​−r73​∣=2r.

The relative normal velocity is
vn​=(v74​−v73​)⋅n^.

You found:
vn​=−0.8348.

Negative means approaching.

Yet they're already overlapping by:
0.0926.

So conceptually:

what SHOULD happen:

        collision
           X
          / \
---------/---\---------------- time
        gap=0


what happened:

        should have collided
           X
           |
           |   simulator continues
           |------->

particles penetrate

          overlap = 0.0926

The most likely broad category is now:
a particle-particle collision event was missed or invalidated​

rather than anything involving FFT.

Possible subcauses include event queue invalidation after another collision, stale event generation counters, near-simultaneous collisions, or an incorrect quadratic/root selection.

We shouldn't choose which one until we see the event history.
9. Can we rerun ONLY this trajectory?
Yes — and that's exactly what we should do.

But there is a subtlety.

Your command uses:

--seed-base 2026082200

The low-density workflow passes approximately:

base seed = 2026082202

and your C program then deterministically generates the actual per-trajectory seed

2381038820

So simply doing:

--seed=2381038820

may not mean “use this exact trajectory seed.” It depends on how 00ALLINONE interprets --seed.

It may interpret that as a new master seed and derive yet another trajectory seed.

So don't assume that yet.
10. We can nevertheless exploit its determinism immediately

The simplest experiment requiring no new code is to run only:
L0​=150,M=50

with the same low-batch master seed and enough repeats to reach repeat 3.

Something conceptually like:

./hspist3/00ALLINONE \
  --mode=edmd \
  --experiment=speed_of_sound \
  --headless \
  --quiet \
  --kbt1 \
  --particles=100 \
  --particles-boxes=50,50 \
  --height=10.0 \
  --particle-radius=0.5 \
  --wall-thickness=0.05 \
  --wall-thickness-vis=0.05 \
  --lengths=150 \
  --wall-masses=50 \
  --repeats=4 \
  --seed=2026082202 \
  --wall-hold-steps=2000 \
  --fixed-dt=0.4 \
  --edmd-acc=0 \
  --target-oscillations=25 \
  --oscillation-safety=1.5 \
  --oscillation-min-steps=10000 \
  --oscillation-max-steps=10000000

But I would first check what seed repeat 3 prints.

If it prints:

seed = 2381038820

perfect.

Then we have reduced a 350-run experiment to four runs.

If it doesn't, the seed generator depends on the position of L0​,M in the full grid.

Then we should add a tiny exact-seed override.
11. The ideal debugging flag

When Codex/Claude is available, I would ask for something like:

--speed-sound-exact-seed=2381038820

plus:

--edmd-debug-particles=73,74
--edmd-debug-event-history=50

Then the entire problem becomes:

./00ALLINONE \
    ...same physics... \
    --lengths=150 \
    --wall-masses=50 \
    --repeats=1 \
    --speed-sound-exact-seed=2381038820 \
    --edmd-debug-particles=73,74 \
    --edmd-debug-event-history=50

This should not change any dynamics.

It only changes instrumentation.
12. What the debugging run should save

Something like:

event #845372
time = ...

event type = PP
particles = 74,18

before:
  r74 = ...
  v74 = ...

after:
  r74 = ...
  v74 = ...

events invalidated:
  PP 74-73 @ t=...

events rescheduled:
  PP 74-73 -> ???

--------------------------------

event #845373
...

--------------------------------

VALIDATOR FAILURE:

particles 74,73
distance = ...
required = 1.0
overlap = 0.0925995654
v_rel_normal = -0.834756336

Then we can literally point at:

    “Here. Event X changed particle 74, invalidated the 74–73 collision, but the new collision wasn't scheduled.”

or whatever the true cause turns out to be.

That is much better than guessing.
13. And I would debug this BEFORE the 50/100-repeat paper run

The current curve is already excellent enough to design the paper experiment.

So I would freeze the physics now:

CURRENT GOOD CURVE
      │
      ├── preserve raw data
      │
      ├── preserve analysis
      │
      ▼
debug deterministic seed
      │
      ▼
core regression passes
      │
      ▼
FFT robustness tests
      │
      ▼
EOS/phase plotting cleanup
      │
      ▼
c25 → c50 convergence
      │
      ▼
final R=50 or R=100 campaign

That way we do not spend several thousand trajectories on a core for which we already possess a known reproducible failure.

And on the EOS question, the key answer is: yes, terminate the simple fluid EOS curves near the stable-fluid boundary; use Liu/Mier-y-Terán for phase-aware thermodynamic comparison, but don't blindly turn a global pressure EOS into a single fluid sound-speed curve across coexistence and solid. The change of regime should actually become visible in the figure rather than being hidden by one continuous theoretical line.

(.venv) chrisharing@Mac HardDisks % cd /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks

RUN="hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/simulation_oscillation_design_pilot_r11_c25"

rg -n \
  "2381038820|particle_particle_overlap|0.0925995654|Particles 74 and 73" \
  "$RUN/raw_simulations/03_low_eta"
hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/simulation_oscillation_design_pilot_r11_c25/raw_simulations/03_low_eta/speed_of_sound_failures.csv
2:2026-08-22 15:59:52,speed_of_sound,edmd,2381038820,150,50,3,148,2.4666667,"wall_hold","particle_particle_overlap","Particles 74 and 73 overlap by 0.0925995654 px (tolerance 2.4e-05 px, relative normal velocity -0.834756336 px/time).",74,73,-1,-1,0.092599565360004021,2.4000000000000001e-05,"/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/simulation_oscillation_design_pilot_r11_c25/raw_simulations/03_low_eta/invalid_wall_x_positions_L0_1500_wallmassfactor_50_run3.csv","./00ALLINONE --mode=edmd --experiment=speed_of_sound --headless --quiet --kbt1 --particles=100 --particles-boxes=50,50 --height=10.0 --particle-radius=0.5 --wall-thickness=0.05 --wall-thickness-vis=0.05 --lengths=50,75,100,150,200 --wall-masses=10,20,50,100,200,300,500 --repeats=11 --seed=2026082202 --wall-hold-steps=2000 --fixed-dt=0.4 --edmd-acc=0 --target-oscillations=25 --oscillation-safety=1.5 --oscillation-min-steps=10000 --oscillation-max-steps=10000000 --speed-sound-run-dir=experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/simulation_oscillation_design_pilot_r11_c25/raw_simulations/03_low_eta"
(.venv) chrisharing@Mac HardDisks % find "$RUN/raw_simulations/03_low_eta" -type f | \
  grep -E '2381038820|invalid|failure|ledger'
hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/simulation_oscillation_design_pilot_r11_c25/raw_simulations/03_low_eta/invalid_wall_x_positions_L0_1500_wallmassfactor_50_run3.csv
hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/simulation_oscillation_design_pilot_r11_c25/raw_simulations/03_low_eta/speed_of_sound_failures.csv
(.venv) chrisharing@Mac HardDisks %

Yes — this output is much more informative than the previous one. We now know when the bug happens, and that narrows the problem dramatically.

The key row is:

seed = 2381038820
L0 = 150
M = 50
repeat = 3
step = 148
time = 2.4666667
phase = wall_hold

Particles 74 and 73 overlap by 0.0925995654
relative normal velocity = -0.834756336

The biggest discovery: it fails during wall_hold

This is important.

It does not fail after 20,000–50,000 time units during the sound-wave experiment.

It fails at:
t=2.4667​

after only
148 recorded steps​

and specifically during

wall_hold

So the sequence is roughly:

initialize particles
      │
      ▼
hold piston fixed / equilibration
      │
      ├── t = 0
      ├── ...
      ├── t = 2.4667
      │
      X particles 73/74 overlap
      │
      ▼
sound-speed experiment NEVER really starts

This means the failure is almost certainly not caused by:

    the new target-oscillation scheduler;

    low-frequency FFT;

    the 25-cycle duration;

    the spring-wall recollision solver during the actual oscillation experiment.

That is excellent news for interpreting your beautiful sound-speed curve.

It instead points much more narrowly toward the particle-particle EDMD collision machinery during equilibration / wall hold.
And 2.4667 contains another clue

Notice:
60148​=2.4666667.

Exactly.

So whatever validation/output cadence is being used there seems to be approximately
Δtsample​=601​.

That means the validator caught it at its 148th sampled state.

The actual missed collision occurred somewhere before that checkpoint.

And since the overlap is huge,
δ=0.0926,

with inward relative speed
∣vn​∣=0.8348,

a crude estimate of how long ago contact should have happened is
Δtmiss​∼∣vn​∣δ​=0.83480.0926​≈0.111.

So approximately:
tcontact​∼2.4667−0.111≈2.356.

Do not treat 2.356 as an exact collision time, because velocities may have changed during that interval. But it tells us the relevant event history is probably only the last ∼0.1−0.2 time units before failure.

That's fantastic for debugging.
What files did the validator preserve?

You have:

invalid_wall_x_positions_L0_1500_wallmassfactor_50_run3.csv

and

speed_of_sound_failures.csv

Unfortunately, from the filename the first seems to contain the wall trajectory, not the positions and velocities of particles 73 and 74.

Let's inspect it anyway:

RUN="hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/simulation_oscillation_design_pilot_r11_c25"

wc -l \
  "$RUN/raw_simulations/03_low_eta/invalid_wall_x_positions_L0_1500_wallmassfactor_50_run3.csv"

head -n 10 \
  "$RUN/raw_simulations/03_low_eta/invalid_wall_x_positions_L0_1500_wallmassfactor_50_run3.csv"

tail -n 30 \
  "$RUN/raw_simulations/03_low_eta/invalid_wall_x_positions_L0_1500_wallmassfactor_50_run3.csv"

That will tell us exactly what columns were preserved.

My expectation is something like:

time, wall_x, wall_v, ...

If that's all it contains, then we do not yet have enough stored data to identify why particles 73 and 74 missed their collision.

We have enough to reproduce it, though.

And that's the important part.
First thing I would check now: how does the program generate 2381038820?

Before we create any debugging code, find the seed-generation logic.

Run:

cd /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks

rg -n \
  "seed.*repeat|repeat.*seed|run_seed|trajectory_seed|seed_base|rng_seed|speed.*seed" \
  hspist3/00ALLINONE.c hspist3/edmd_core

And specifically:

rg -n "2381038820|2026082202" hspist3

The exact literal probably won't exist in source, but this helps locate relevant logging/code.

Also:

rg -n '"seed"|seed=' hspist3/00ALLINONE.c | head -n 100

What we want to understand is:

master seed 2026082202
        │
        ▼
(L0, M, repeat)
        │
        ▼
trajectory seed 2381038820

If that trajectory seed is a deterministic hash of L0​,M,r, then we can probably reproduce only this one trajectory very easily.
You may already be able to reproduce only four runs

Try a separate debug folder, preserving everything else:

cd /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3

./00ALLINONE \
  --mode=edmd \
  --experiment=speed_of_sound \
  --headless \
  --quiet \
  --kbt1 \
  --particles=100 \
  --particles-boxes=50,50 \
  --height=10.0 \
  --particle-radius=0.5 \
  --wall-thickness=0.05 \
  --wall-thickness-vis=0.05 \
  --lengths=150 \
  --wall-masses=50 \
  --repeats=4 \
  --seed=2026082202 \
  --wall-hold-steps=2000 \
  --fixed-dt=0.4 \
  --edmd-acc=0 \
  --target-oscillations=25 \
  --oscillation-safety=1.5 \
  --oscillation-min-steps=10000 \
  --oscillation-max-steps=10000000 \
  --speed-sound-run-dir=experiments_speed_of_sound/EDMD/mode1_normalized_units/debug_seed_2381038820

Then watch what it prints for run 3.
If it says

run = 3, seed = 2381038820

we've won.

We have reduced:
1460+ trajectories

to:
4.
If it gives a different seed

then the per-run seed depends on the position in the whole L0​×M grid.

That's still easy to fix; we just need to inspect the seed-generation function.
What exactly should we debug?

The critical quantity for particles 73 and 74 is their surface gap:
g73,74​(t)=∣r74​(t)−r73​(t)∣−2r.

Since
r=0.5,

contact occurs at
∣r74​−r73​∣=1.

Therefore:

g > 0       particles separated
g = 0       collision
g < 0       overlap -> invalid

At the failure:
g=−0.0925996.

So their center distance is approximately
d=1−0.0925996=0.9074004​.

They are substantially inside each other.
What EDMD should have done

For two free disks during an event interval, define
r=r74​−r73​,

and
v=v74​−v73​.

A collision occurs when
∣r+vt∣=2r.

Square both sides:
(r+vt)⋅(r+vt)=(2r)2.

Expand:
∣v∣2t2+2(r⋅v)t+(∣r∣2−(2r)2)=0.

That's just a quadratic:
at2+bt+c=0,

where
a=∣v∣2, b=2r⋅v, c=∣r∣2−(2r)2.

EDMD calculates the earliest positive physical root and schedules:

particle 73 ↔ particle 74 collision at t_collision

Something prevented that from happening.
The bug is probably in one of two places

Now that the spring wall is essentially ruled out for this failure, I'd divide the search into:
A. Collision prediction

The PP solver failed to predict the valid root.

73 and 74 approach
       │
       ▼
quadratic collision solver
       │
       X says "no collision"

or
B. Event-calendar invalidation/rescheduling

More likely for a rare event:

PP(73,74) scheduled
       │
       ▼
particle 73 collides with particle 21
       │
       ▼
old PP(73,74) correctly becomes stale
       │
       ▼
should recompute PP(73,74)
       │
       X
not rescheduled / invalidated wrongly

Then particles 73 and 74 eventually meet and pass through each other.

The fact that this happens only for one particular seed out of hundreds makes B particularly plausible, although that's still a hypothesis.
The ideal debug output is tiny

We do NOT need to log the entire million-event simulation.

The failure occurs at
t=2.4667.

So for the exact seed, keep a circular buffer of perhaps the last 50–200 EDMD events.

On failure dump:

FAILURE t=2.4666667

particle 73:
    position
    velocity
    generation/event counter
    last collision time

particle 74:
    position
    velocity
    generation/event counter
    last collision time

pair 73-74:
    center distance
    surface gap
    relative velocity
    relative normal velocity
    predicted collision discriminant
    predicted collision root

LAST 100 EVENTS:

#...
t=2.301...
PP particle ...
...
t=2.35...
...
t=2.46...

Then the bug should become visible.
Something else important: this strengthens the current speed-of-sound result

Counterintuitively, yes.

Before this search we had:

    “One trajectory somewhere in the low-density campaign failed.”

Now we know:

    “One deterministic initial condition triggers a PP-core defect extremely early during the fixed-wall equilibration stage; that trajectory is rejected before entering the scientific FFT dataset.”

The validity system did exactly what it was designed to do:

bad physics
   ↓
validator catches it
   ↓
trace marked invalid
   ↓
not published into valid set
   ↓
not used in exploratory fit

That's far better than silently contaminating the sound-speed curve.

Of course, we still need to fix the underlying EDMD bug before the production campaign.
And I would now slightly change the debugging priority

We no longer need to investigate FFT first. The FFT result looks excellent and this failure clearly occurs before FFT measurement.

So:

1. reproduce 2381038820 in minimal run
              ↓
2. dump event history for particles 73/74
              ↓
3. identify missed PP event
              ↓
4. smallest possible EDMD fix
              ↓
5. permanent regression test:
   seed 2381038820 must pass
              ↓
6. rerun c25 subset and confirm curve unchanged
              ↓
7. c50 FFT convergence tests
              ↓
8. final EOS/phase treatment
              ↓
9. publication statistics

That is now a very clean development path.

And do not delete either the original r10 or r11 failure files. Having the exact same failure from two independent campaigns is excellent provenance for the eventual development log.

tail -n 30 \
  "$RUN/raw_simulations/03_low_eta/invalid_wall_x_positions_L0_1500_wallmassfactor_50_run3.csv"
       1 hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/simulation_oscillation_design_pilot_r11_c25/raw_simulations/03_low_eta/invalid_wall_x_positions_L0_1500_wallmassfactor_50_run3.csv
Time,Wall_X,Displacement(σ),Left_Count,Right_Count,L0,eta,Center_X(σ),Seed,Target_Oscillations,Predicted_Frequency,Planned_Steps,Planned_Duration
Time,Wall_X,Displacement(σ),Left_Count,Right_Count,L0,eta,Center_X(σ),Seed,Target_Oscillations,Predicted_Frequency,Planned_Steps,Planned_Duration
(.venv) chrisharing@Mac HardDisks % cd /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3

./00ALLINONE \
  --mode=edmd \
  --experiment=speed_of_sound \
  --headless \
  --quiet \
  --kbt1 \
  --particles=100 \
  --particles-boxes=50,50 \
  --height=10.0 \
  --particle-radius=0.5 \
  --wall-thickness=0.05 \
  --wall-thickness-vis=0.05 \
  --lengths=150 \
  --wall-masses=50 \
  --repeats=4 \
  --seed=2026082202 \
  --wall-hold-steps=2000 \
  --fixed-dt=0.4 \
  --edmd-acc=0 \
  --target-oscillations=25 \
  --oscillation-safety=1.5 \
  --oscillation-min-steps=10000 \
  --oscillation-max-steps=10000000 \
  --speed-sound-run-dir=experiments_speed_of_sound/EDMD/mode1_normalized_units/debug_seed_2381038820
📁 Output folder: /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/debug_seed_2381038820
Planning each experiment for 25 fundamental oscillations (safety 1.500, clamp 10000..10000000 samples)
🔬 Running: L0 = 150.0, M = 50*m, run = 0, seed = 3828468398
🔍 Initial wall_x = 3800.000, vx_wall = 0.000000
🕒 target=25 cycles  f_pred=0.001714989  samples=1311962  T=21866 sigma-time
✅ Finished valid run 0

🔬 Running: L0 = 150.0, M = 50*m, run = 1, seed = 3010016154
🔍 Initial wall_x = 3800.000, vx_wall = 0.000000
🕒 target=25 cycles  f_pred=0.001714989  samples=1311962  T=21866 sigma-time
✅ Finished valid run 1

🔬 Running: L0 = 150.0, M = 50*m, run = 2, seed = 2086766775
🔍 Initial wall_x = 3800.000, vx_wall = 0.000000
🕒 target=25 cycles  f_pred=0.001714989  samples=1311962  T=21866 sigma-time
✅ Finished valid run 2

🔬 Running: L0 = 150.0, M = 50*m, run = 3, seed = 208027055
🔍 Initial wall_x = 3800.000, vx_wall = 0.000000
🕒 target=25 cycles  f_pred=0.001714989  samples=1311962  T=21866 sigma-time
✅ Finished valid run 3

Speed-of-sound batch complete: 4 valid, 0 invalid, 4 requested.
(.venv) chrisharing@Mac hspist3 % cd /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks

rg -n \
  "seed.*repeat|repeat.*seed|run_seed|trajectory_seed|seed_base|rng_seed|speed.*seed" \
  hspist3/00ALLINONE.c hspist3/edmd_core
hspist3/00ALLINONE.c
14460:                                      unsigned int run_seed,
14467:        "timestamp,experiment,sim_mode,seed,L0,wall_mass_factor,repeat,step,time_sigma,"
14483:            timestamp, experiment_sim_mode_name(), run_seed, (double)L0,
14501:static unsigned int speed_sound_run_seed(unsigned int base_seed,
14717:                const unsigned int run_seed =
14718:                    speed_sound_run_seed(cli_seed, l, m, r);
14719:                srand(run_seed);
14790:                       L0, wall_mass_factor, r, run_seed);
15029:                                    run_seed, cli_speed_sound_target_oscillations,
15120:                                    run_seed, cli_speed_sound_target_oscillations,
15188:                        wall_mass_factor, r, run_seed, failure);
15205:                        wall_mass_factor, r, run_seed,
(.venv) chrisharing@Mac HardDisks % rg -n "2381038820|2026082202" hspist3
hspist3/00_COMMAND.md
7:./00ALLINONE --mode=edmd --experiment=speed_of_sound --headless --quiet --kbt1 --particles=100 --particles-boxes=50,50 --height=10.0 --particle-radius=0.5 --wall-thickness=0.05 --wall-thickness-vis=0.05 --lengths=150 --wall-masses=50 --repeats=4 --seed=2026082202 --wall-hold-steps=2000 --fixed-dt=0.4 --edmd-acc=0 --target-oscillations=25 --oscillation-safety=1.5 --oscillation-min-steps=10000 --oscillation-max-steps=10000000 --speed-sound-run-dir=experiments_speed_of_sound/EDMD/mode1_normalized_units/debug_seed_2381038820

hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/debug_seed_2381038820/00_COMMAND.md
7:./00ALLINONE --mode=edmd --experiment=speed_of_sound --headless --quiet --kbt1 --particles=100 --particles-boxes=50,50 --height=10.0 --particle-radius=0.5 --wall-thickness=0.05 --wall-thickness-vis=0.05 --lengths=150 --wall-masses=50 --repeats=4 --seed=2026082202 --wall-hold-steps=2000 --fixed-dt=0.4 --edmd-acc=0 --target-oscillations=25 --oscillation-safety=1.5 --oscillation-min-steps=10000 --oscillation-max-steps=10000000 --speed-sound-run-dir=experiments_speed_of_sound/EDMD/mode1_normalized_units/debug_seed_2381038820

hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/debug_seed_2381038820/speed_of_sound_batch_status.json
3:  "base_seed": 2026082202,

hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/debug_seed_2381038820/00_COMMAND_AND_PLOT_COMMAND.md
6:./00ALLINONE --mode=edmd --experiment=speed_of_sound --headless --quiet --kbt1 --particles=100 --particles-boxes=50,50 --height=10.0 --particle-radius=0.5 --wall-thickness=0.05 --wall-thickness-vis=0.05 --lengths=150 --wall-masses=50 --repeats=4 --seed=2026082202 --wall-hold-steps=2000 --fixed-dt=0.4 --edmd-acc=0 --target-oscillations=25 --oscillation-safety=1.5 --oscillation-min-steps=10000 --oscillation-max-steps=10000000 --speed-sound-run-dir=experiments_speed_of_sound/EDMD/mode1_normalized_units/debug_seed_2381038820
15:  --dir hspist3//Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/debug_seed_2381038820 \
57:hspist3//Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/debug_seed_2381038820/analysis_by_eta/final_plots/FINAL speed_of_sound_on_packing_fracture.pdf

hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/debug_seed_2381038820/01_PLOT_COMMANDS.md
7:  --dir hspist3//Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/debug_seed_2381038820 \
39:hspist3//Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/debug_seed_2381038820/analysis_by_eta/final_plots/FINAL speed_of_sound_on_packing_fracture.pdf
40:hspist3//Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/debug_seed_2381038820/analysis_by_eta/final_plots/combined_speed_of_sound_summary.csv
41:hspist3//Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/debug_seed_2381038820/analysis_by_eta/final_plots/speed_of_sound_freq_vs_kterm_combined.pdf

hspist3/experiments_speed_of_sound/EDMD/mode0_real_units/simulation_21_08_26_16_56_14/00_COMMAND.md
7:./00ALLINONE --mode=edmd --experiment=speed_of_sound --headless --quiet --kbt1 --particles=100 --particles-boxes=50,50 --height=10.0 --particle-radius=0.5 --wall-thickness=0.05 --wall-thickness-vis=0.05 --lengths=50,75,100,150,200 --wall-masses=10,20,50,100,200,300,500 --repeats=10 --seed=2026082202 --wall-hold-steps=2000 --fixed-dt=0.4 --edmd-acc=0 --target-oscillations=25 --oscillation-safety=1.5 --oscillation-min-steps=10000 --oscillation-max-steps=10000000 --speed-sound-run-dir=experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/simulation_oscillation_design_pilot_r10_c25/raw_simulations/03_low_eta

hspist3/experiments_speed_of_sound/EDMD/mode0_real_units/simulation_21_08_26_16_56_14/00_COMMAND_AND_PLOT_COMMAND.md
6:./00ALLINONE --mode=edmd --experiment=speed_of_sound --headless --quiet --kbt1 --particles=100 --particles-boxes=50,50 --height=10.0 --particle-radius=0.5 --wall-thickness=0.05 --wall-thickness-vis=0.05 --lengths=50,75,100,150,200 --wall-masses=10,20,50,100,200,300,500 --repeats=10 --seed=2026082202 --wall-hold-steps=2000 --fixed-dt=0.4 --edmd-acc=0 --target-oscillations=25 --oscillation-safety=1.5 --oscillation-min-steps=10000 --oscillation-max-steps=10000000 --speed-sound-run-dir=experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/simulation_oscillation_design_pilot_r10_c25/raw_simulations/03_low_eta

hspist3/experiments_speed_of_sound/EDMD/mode0_real_units/simulation_22_08_26_18_36_56/00_COMMAND.md
7:./00ALLINONE --mode=edmd --experiment=speed_of_sound --headless --quiet --kbt1 --particles=100 --particles-boxes=50,50 --height=10.0 --particle-radius=0.5 --wall-thickness=0.05 --wall-thickness-vis=0.05 --lengths=150 --wall-masses=50 --repeats=4 --seed=2026082202 --wall-hold-steps=2000 --fixed-dt=0.4 --edmd-acc=0 --target-oscillations=25 --oscillation-safety=1.5 --oscillation-min-steps=10000 --oscillation-max-steps=10000000 --speed-sound-run-dir=experiments_speed_of_sound/EDMD/mode1_normalized_units/debug_seed_2381038820

hspist3/experiments_speed_of_sound/EDMD/mode0_real_units/simulation_22_08_26_18_36_56/00_COMMAND_AND_PLOT_COMMAND.md
6:./00ALLINONE --mode=edmd --experiment=speed_of_sound --headless --quiet --kbt1 --particles=100 --particles-boxes=50,50 --height=10.0 --particle-radius=0.5 --wall-thickness=0.05 --wall-thickness-vis=0.05 --lengths=150 --wall-masses=50 --repeats=4 --seed=2026082202 --wall-hold-steps=2000 --fixed-dt=0.4 --edmd-acc=0 --target-oscillations=25 --oscillation-safety=1.5 --oscillation-min-steps=10000 --oscillation-max-steps=10000000 --speed-sound-run-dir=experiments_speed_of_sound/EDMD/mode1_normalized_units/debug_seed_2381038820

hspist3/experiments_speed_of_sound/EDMD/mode0_real_units/simulation_22_08_26_15_47_12/00_COMMAND.md
7:./00ALLINONE --mode=edmd --experiment=speed_of_sound --headless --quiet --kbt1 --particles=100 --particles-boxes=50,50 --height=10.0 --particle-radius=0.5 --wall-thickness=0.05 --wall-thickness-vis=0.05 --lengths=50,75,100,150,200 --wall-masses=10,20,50,100,200,300,500 --repeats=11 --seed=2026082202 --wall-hold-steps=2000 --fixed-dt=0.4 --edmd-acc=0 --target-oscillations=25 --oscillation-safety=1.5 --oscillation-min-steps=10000 --oscillation-max-steps=10000000 --speed-sound-run-dir=experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/simulation_oscillation_design_pilot_r11_c25/raw_simulations/03_low_eta

hspist3/experiments_speed_of_sound/EDMD/mode0_real_units/simulation_22_08_26_15_47_12/00_COMMAND_AND_PLOT_COMMAND.md
6:./00ALLINONE --mode=edmd --experiment=speed_of_sound --headless --quiet --kbt1 --particles=100 --particles-boxes=50,50 --height=10.0 --particle-radius=0.5 --wall-thickness=0.05 --wall-thickness-vis=0.05 --lengths=50,75,100,150,200 --wall-masses=10,20,50,100,200,300,500 --repeats=11 --seed=2026082202 --wall-hold-steps=2000 --fixed-dt=0.4 --edmd-acc=0 --target-oscillations=25 --oscillation-safety=1.5 --oscillation-min-steps=10000 --oscillation-max-steps=10000000 --speed-sound-run-dir=experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/simulation_oscillation_design_pilot_r11_c25/raw_simulations/03_low_eta

hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/simulation_oscillation_design_pilot_r10_c25/raw_simulations/03_low_eta/00_COMMAND.md
7:./00ALLINONE --mode=edmd --experiment=speed_of_sound --headless --quiet --kbt1 --particles=100 --particles-boxes=50,50 --height=10.0 --particle-radius=0.5 --wall-thickness=0.05 --wall-thickness-vis=0.05 --lengths=50,75,100,150,200 --wall-masses=10,20,50,100,200,300,500 --repeats=10 --seed=2026082202 --wall-hold-steps=2000 --fixed-dt=0.4 --edmd-acc=0 --target-oscillations=25 --oscillation-safety=1.5 --oscillation-min-steps=10000 --oscillation-max-steps=10000000 --speed-sound-run-dir=experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/simulation_oscillation_design_pilot_r10_c25/raw_simulations/03_low_eta

hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/simulation_oscillation_design_pilot_r10_c25/raw_simulations/03_low_eta/speed_of_sound_failures.csv
2:2026-08-21 17:08:17,speed_of_sound,edmd,2381038820,150,50,3,148,2.4666667,"wall_hold","particle_particle_overlap","Particles 74 and 73 overlap by 0.0925995654 px (tolerance 2.4e-05 px, relative normal velocity -0.834756336 px/time).",74,73,-1,-1,0.092599565360004021,2.4000000000000001e-05,"/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/simulation_oscillation_design_pilot_r10_c25/raw_simulations/03_low_eta/invalid_wall_x_positions_L0_1500_wallmassfactor_50_run3.csv","./00ALLINONE --mode=edmd --experiment=speed_of_sound --headless --quiet --kbt1 --particles=100 --particles-boxes=50,50 --height=10.0 --particle-radius=0.5 --wall-thickness=0.05 --wall-thickness-vis=0.05 --lengths=50,75,100,150,200 --wall-masses=10,20,50,100,200,300,500 --repeats=10 --seed=2026082202 --wall-hold-steps=2000 --fixed-dt=0.4 --edmd-acc=0 --target-oscillations=25 --oscillation-safety=1.5 --oscillation-min-steps=10000 --oscillation-max-steps=10000000 --speed-sound-run-dir=experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/simulation_oscillation_design_pilot_r10_c25/raw_simulations/03_low_eta"

hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/simulation_oscillation_design_pilot_r10_c25/raw_simulations/03_low_eta/00_COMMAND_AND_PLOT_COMMAND.md
6:./00ALLINONE --mode=edmd --experiment=speed_of_sound --headless --quiet --kbt1 --particles=100 --particles-boxes=50,50 --height=10.0 --particle-radius=0.5 --wall-thickness=0.05 --wall-thickness-vis=0.05 --lengths=50,75,100,150,200 --wall-masses=10,20,50,100,200,300,500 --repeats=10 --seed=2026082202 --wall-hold-steps=2000 --fixed-dt=0.4 --edmd-acc=0 --target-oscillations=25 --oscillation-safety=1.5 --oscillation-min-steps=10000 --oscillation-max-steps=10000000 --speed-sound-run-dir=experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/simulation_oscillation_design_pilot_r10_c25/raw_simulations/03_low_eta

hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/simulation_oscillation_design_pilot_r10_c25/raw_simulations/03_low_eta/speed_of_sound_batch_status.json
3:  "base_seed": 2026082202,

hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/simulation_oscillation_design_pilot_r11_c25/raw_simulations/03_low_eta/00_COMMAND.md
7:./00ALLINONE --mode=edmd --experiment=speed_of_sound --headless --quiet --kbt1 --particles=100 --particles-boxes=50,50 --height=10.0 --particle-radius=0.5 --wall-thickness=0.05 --wall-thickness-vis=0.05 --lengths=50,75,100,150,200 --wall-masses=10,20,50,100,200,300,500 --repeats=11 --seed=2026082202 --wall-hold-steps=2000 --fixed-dt=0.4 --edmd-acc=0 --target-oscillations=25 --oscillation-safety=1.5 --oscillation-min-steps=10000 --oscillation-max-steps=10000000 --speed-sound-run-dir=experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/simulation_oscillation_design_pilot_r11_c25/raw_simulations/03_low_eta

hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/simulation_oscillation_design_pilot_r11_c25/raw_simulations/03_low_eta/speed_of_sound_failures.csv
2:2026-08-22 15:59:52,speed_of_sound,edmd,2381038820,150,50,3,148,2.4666667,"wall_hold","particle_particle_overlap","Particles 74 and 73 overlap by 0.0925995654 px (tolerance 2.4e-05 px, relative normal velocity -0.834756336 px/time).",74,73,-1,-1,0.092599565360004021,2.4000000000000001e-05,"/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/simulation_oscillation_design_pilot_r11_c25/raw_simulations/03_low_eta/invalid_wall_x_positions_L0_1500_wallmassfactor_50_run3.csv","./00ALLINONE --mode=edmd --experiment=speed_of_sound --headless --quiet --kbt1 --particles=100 --particles-boxes=50,50 --height=10.0 --particle-radius=0.5 --wall-thickness=0.05 --wall-thickness-vis=0.05 --lengths=50,75,100,150,200 --wall-masses=10,20,50,100,200,300,500 --repeats=11 --seed=2026082202 --wall-hold-steps=2000 --fixed-dt=0.4 --edmd-acc=0 --target-oscillations=25 --oscillation-safety=1.5 --oscillation-min-steps=10000 --oscillation-max-steps=10000000 --speed-sound-run-dir=experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/simulation_oscillation_design_pilot_r11_c25/raw_simulations/03_low_eta"

hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/simulation_oscillation_design_pilot_r11_c25/raw_simulations/03_low_eta/00_COMMAND_AND_PLOT_COMMAND.md
6:./00ALLINONE --mode=edmd --experiment=speed_of_sound --headless --quiet --kbt1 --particles=100 --particles-boxes=50,50 --height=10.0 --particle-radius=0.5 --wall-thickness=0.05 --wall-thickness-vis=0.05 --lengths=50,75,100,150,200 --wall-masses=10,20,50,100,200,300,500 --repeats=11 --seed=2026082202 --wall-hold-steps=2000 --fixed-dt=0.4 --edmd-acc=0 --target-oscillations=25 --oscillation-safety=1.5 --oscillation-min-steps=10000 --oscillation-max-steps=10000000 --speed-sound-run-dir=experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/simulation_oscillation_design_pilot_r11_c25/raw_simulations/03_low_eta

hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/simulation_oscillation_design_pilot_r11_c25/raw_simulations/03_low_eta/speed_of_sound_batch_status.json
3:  "base_seed": 2026082202,
(.venv) chrisharing@Mac HardDisks % rg -n '"seed"|seed=' hspist3/00ALLINONE.c | head -n 100
2846:    printf("  --seed=value                RNG seed (integer) or 'time'\n");
(.venv) chrisharing@Mac HardDisks %

Yes — this explains exactly why your minimal 4-run reproduction did not hit the bad seed.

The important clue is here:

14501:static unsigned int speed_sound_run_seed(unsigned int base_seed,
...
14717:const unsigned int run_seed =
14718:    speed_sound_run_seed(cli_seed, l, m, r);

So the trajectory seed is not simply based on
(L0​,M,repeat)

as physical values. It is generated from the array indices
(l,m,r)

inside the complete sweep.

When you changed

lengths = 50,75,100,150,200
wall masses = 10,20,50,100,200,300,500

to

lengths = 150
wall masses = 50

you changed the indices.

Originally:

L0 = 150  -> l = 3
M  = 50   -> m = 2
repeat 3  -> r = 3

but in the reduced run:

L0 = 150  -> l = 0
M  = 50   -> m = 0
repeat 3  -> r = 3

So naturally:

speed_sound_run_seed(2026082202, 3,2,3)
    = 2381038820

but

speed_sound_run_seed(2026082202, 0,0,3)
    = 208027055

That's almost certainly what we are seeing.

And we learned something else useful: the preserved invalid CSV contains only the header. So no, unfortunately we do not already have particle/event history from the failed trajectory. It failed during wall_hold before any normal trace rows were published.
Let's inspect the seed function now

Run exactly this:

cd /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks

sed -n '14490,14535p' hspist3/00ALLINONE.c

sed -n '14700,14730p' hspist3/00ALLINONE.c

Paste that output here.

The first section should show the complete:

speed_sound_run_seed(...)

function.

Once I see that, we can probably calculate a master --seed for the isolated (l=0,m=0,r=3) run that produces exactly 2381038820, meaning we may be able to reproduce the bug with only four trajectories without changing any code.

If the hash isn't conveniently invertible, the next-best option is a tiny temporary exact-seed/debug override.
Also locate the particle-particle event code

Run this too:

rg -n \
  "particle.*particle|predict.*collision|collision.*time|schedule.*collision|event.*particle|generation|invalidate|quadratic" \
  hspist3/edmd_core hspist3/00ALLINONE.c | head -n 250

We're looking for three separate pieces:

1. PP collision-time prediction
       ↓
2. event insertion / scheduling
       ↓
3. stale-event invalidation + rescheduling

Because our current evidence says the failure is probably somewhere in that chain.
The debugging picture is now actually very clean

We know:

original low sweep
base seed = 2026082202

L0 list:
index 0 -> 50
index 1 -> 75
index 2 -> 100
index 3 -> 150    ← target

M list:
index 0 -> 10
index 1 -> 20
index 2 -> 50     ← target

repeat:
index 3           ← target

             ↓

speed_sound_run_seed(
    2026082202,
    3,
    2,
    3
)

             ↓

2381038820

             ↓

wall_hold
t = 2.4666667

             ↓

particles 73 / 74
overlap = 0.0925995654
v_rel,n = -0.834756336

That is an extremely well-defined bug now.

And importantly, your reduced test producing four valid trajectories does not contradict reproducibility at all. It simply generated four different initial conditions.

Once we reproduce the exact seed in isolation, we don't need to run even one complete sound-speed oscillation. The bug occurs at
t=2.4667

during the initial 2000-step wall-hold stage, while a normal L0​=150,M=50 production trajectory would eventually run to
Tobs​=21866.

So the actual debug run should take almost no time once isolated.

After you paste those two sed outputs, we can determine the cleanest way to recreate exactly seed 2381038820 and then inspect the last particle-particle events before t≈2.47.


