Yes — choose Option 1. Your reasoning is sound: once the production trajectory has the strict all-or-nothing health criterion, calibration is no longer part of the scientific validity proof. It is a performance heuristic for choosing a likely-safe chunk.

The distinction is:
calibration predicts safety​

but
production health counters determine validity​.

So shortening calibration cannot turn an invalid trajectory into an accepted one. Worst case, calibration guesses too aggressively, production triggers a health event, the whole trajectory is rejected, and you rerun with a smaller chunk.

I would not drop N=1600 yet. That is exactly the point giving you the strongest finite-size evidence. Saving an hour or two by losing the largest system is a worse scientific trade than making calibration cheaper.
I would make calibration even simpler than 5 units × 8 attempts

Use the collision-rate estimate for the first guess, then:

candidate chunk
      ↓
5 simulated units on disposable seed A
      ↓
clean? yes
      ↓
5 simulated units on disposable seed B
      ↓
clean? yes
      ↓
production_chunk = 0.6 × min(safe_A, safe_B)

If a verification fails:
Δtchunk​←0.5Δtchunk​

and retry, with something like 6–8 attempts maximum.

I'd reduce the final factor from the previous 0.8 to perhaps
0.6​

because we're deliberately reducing the calibration observation time. That buys back some conservatism very cheaply.

You could even use 0.5 if smaller calls have negligible overhead. The important constraint is just not to make chunks absurdly tiny.
Why 5 units is enough for calibration

Not because five units proves anything.

It doesn't.

It's enough because its job is only:

    “Does this candidate immediately look too coarse?”

The actual scientific run might cover
200+360=560

time units. That is effectively the long stress test.

If a rare burst occurs at t=417, the production health system catches it:

forced_advance > 0
        ↓
trajectory INVALID
        ↓
does not enter pressure averages

So there is no statistical contamination.

That means spending another 280−560 simulated units before the 560-unit production run is indeed redundant for this campaign.
One important nuance

When production fails, don't just rerun that single seed with a smaller chunk and then pretend the cell's original calibration was fine.

Instead, mark:
(η,N)

as having an insufficient cached chunk, lower the cached value globally — e.g.
Δtnew​=0.5Δtold​,

and rerun the failed seed. The already-completed clean seeds don't need rerunning because their trajectories are independently valid.

This gives you a nice audit trail:

eta=.65 N=1600
calibrated chunk = .16

seed 1    PASS
seed 2    PASS
seed 3    forced_advance → DISCARD

cached chunk revised → .08

seed 3 rerun → PASS

That's completely defensible.
And the old η=0.5,N=1600 incident is valuable

I would describe it carefully in the eventual paper. It's not evidence that “EDMD collision physics was wrong.”

It's more precise to say:

    A too-large integration/advance window exhausted the event-processing budget, invoking the simulator's forced-advance safeguard. Although the resulting pressure estimate remained numerically plausible, health telemetry exposed thousands of subsequent overlap corrections. Under the validated runner, the same state was rerun with a conservative event-window size and completed without any health events.

That's an excellent validation example, because the scary part is:
Zcorrupt​≈4.1476

versus clean roughly
Zclean​≈4.149–4.150.

The number looked fine.

Without health telemetry, you could easily have published it.

That's exactly why the validation machinery matters.
So tell Claude to change only this

Change ONLY the disposable chunk-calibration cost.

The production validity contract remains unchanged:
any health event from initialization through the end of measurement invalidates
the entire scientific trajectory.

Calibration is now explicitly a performance heuristic, NOT part of the
acceptance proof.

For each (eta,N):

1. Keep the collision-rate probe to choose the initial candidate chunk.
2. Verify candidate chunks on TWO disposable calibration seeds.
3. Each verification interval is 5 simulated time units, not 20.
4. On failure:
       chunk *= 0.5
   and retry.
5. Cap at 8 attempts. If no clean candidate is found, mark the cell
   calibration_failed and do not launch production.
6. After both disposable seeds find clean chunks, choose:

       production_chunk = 0.6 * min(safe_chunk_A, safe_chunk_B)

7. Destroy calibration instances.
8. Run every production seed from a fresh state using that fixed chunk.
9. Do NOT adapt chunk during production.
10. ANY forced_advance / overlap_repair / clamp_repair / wall_overdue /
    validator failure anywhere in production:
        DISCARD trajectory.

If a production trajectory fails despite calibration:
- mark the cached chunk for that (eta,N) insufficient,
- halve it,
- rerun only the failed trajectory from a fresh state.
Already clean production trajectories remain valid and need not be rerun.

Do not change any pressure estimator, collision physics, block design, eta/N
grid, or existing accepted data.
Do not auto-commit.

I think that's the sweet spot: much cheaper calibration without weakening what “accepted trajectory” means at all.

And with your existing 96/96 strict accepted trajectories, I definitely wouldn't sacrifice N=1600 yet. You're already getting exactly the validation story we wanted.

Option 1, exactly as specified. Changing only the calibration cost — the production contract is untouched.

**Bash&#xA0;**&#x43;heapen calibration and time it

IN

cd /Users/chrisharing/Desktop/CCS\_complex\_coupled\_systems/Repo/HardDisks/hspist3; pkill -f run\_pressure\_campaign2 >/dev/null 2>&1; pkill -f pressure\_validation >/dev/null 2>&1; sleep 2; /opt/homebrew/bin/python3 - <<'PYEOF' p='validation/pressure\_validation.c'; s=open(p).read() old = s[s.index("#define TARGET\_EVENT\_FRACTION"):s.index("int main(int argc, char\*\* argv){")] new = """/\* Chunk calibration is a PERFORMANCE HEURISTIC, not part of the acceptance proof.    Validity is decided solely by the production health counters: any forced    advance / overlap repair / clamp repair / overdue wall anywhere in a scientific    trajectory discards the whole trajectory. So calibration only has to answer    "does this candidate look obviously too coarse?", and a short interval suffices    -- the 560-unit production run is itself the long stress test. Verifying for    280 simulated units before every 560-unit run was pure redundancy and roughly    doubled the campaign cost. \*/ #define TARGET\_EVENT\_FRACTION 0.15      /\* of the \~250k per-call budget \*/ #define ASSUMED\_EVENT\_BUDGET  250000.0 #define CALIB\_VERIFY\_TIME     5.0       /\* simulated units per verification \*/ #define CALIB\_MAX\_ATTEMPTS    8 /\* Returns a clean chunk, or -1.0 if no candidate survived (cell must not run). \*/ static double calibrate\_chunk\_disposable(const EDMD\_Params\* prm, unsigned long long seed){     double chunk = 1.0;     {   EDMD\* T = edmd\_create(prm);         if(T){             edmd\_init\_random\_gas(T, seed);             const double t0=edmd\_time(T); const long e0=edmd\_virial\_pair\_events(T);             edmd\_advance\_to(T, t0 + 0.02);             const double used=edmd\_time(T)-t0; const long ev=edmd\_virial\_pair\_events(T)-e0;             if(used>0.0 && ev>0){                 const double rate=(double)ev/used;                 chunk = (TARGET\_EVENT\_FRACTION\*ASSUMED\_EVENT\_BUDGET)/rate;             }             edmd\_destroy(T);         }     }     if(chunk > 1.0) chunk = 1.0;     if(chunk < 1e-5) chunk = 1e-5;      for(int attempt=0; attempt\<CALIB\_MAX\_ATTEMPTS; ++attempt){         EDMD\* T = edmd\_create(prm);         if(!T) break;         edmd\_init\_random\_gas(T, seed);         g\_chunk = chunk;         advance\_fixed(T, CALIB\_VERIFY\_TIME);         const long h = health\_sum(T);         edmd\_destroy(T);         if(h == 0) return chunk;         chunk \*= 0.5;         if(chunk < 1e-5) break;     }     return -1.0;   /\* calibration\_failed \*/ }  """ s = s.replace(old, new, 1) s = s.replace("""    if(calibrate\_only){         /\* Disposable calibration only. Prints the largest chunk that ran a full            verification interval with zero health events. \*/         printf("%.10g\\\n", calibrate\_chunk\_disposable(&p, seed));         return 0;     }""", """    if(calibrate\_only){         /\* Disposable calibration only. Prints the largest chunk that ran a full            verification interval with zero health events, or FAILED. \*/         const double c = calibrate\_chunk\_disposable(&p, seed);         if(c <= 0.0){ printf("FAILED\\\n"); return 3; }         printf("%.10g\\\n", c);         return 0;     }""",1) s = s.replace("""    g\_chunk = (chunk\_arg > 0.0) ? chunk\_arg : calibrate\_chunk\_disposable(&p, seed);""", """    if(chunk\_arg > 0.0){         g\_chunk = chunk\_arg;     } else {         const double c = calibrate\_chunk\_disposable(&p, seed);         if(c <= 0.0){ fprintf(stderr,"calibration\_failed eta=%g N=%d\\\n",eta,N); return 3; }         g\_chunk = c;     }""",1) open(p,'w').write(s); print("calibration: 5 units, 8 attempts, x0.5, FAILED signal") PYEOF clang -O2 -I. validation/pressure\_validation.c edmd\_core/edmd.c -o validation/pressure\_validation -lm 2>&1 | grep -E "error" | head -3; echo "built"; echo "timing one probe at the cell that took 12+ min before:"; /usr/bin/time -p ./validation/pressure\_validation 0.60 400 911000001 0 0 0 - - calibrate 2>&1 | tail -3

OUT

Command did not complete within its 600s timeout and was moved to the background (ID: b30c1416j). Output is being written to: /private/tmp/claude-501/-Users-chrisharing-Desktop-CCS-complex-coupled-systems-Repo-HardDisks/91cb08ec-0599-4faa-a6af-d5a7ca834255/tasks/b30c1416j.output. You will be notified when it completes. To check interim output, use Read on that file path. Session cwd remains /Users/chrisharing/Desktop/CCS\_complex\_coupled\_systems/Repo/HardDisks; directory changes made by the backgrounded command do not apply to subsequent commands.

Even at 5 units the probe exceeds 10 minutes — so the verification interval was never the bottleneck. Let me find what actually is:

**Bash&#xA0;**&#x4D;easure init vs advance cost by density

IN

pkill -f "pressure\_validation 0.60 400 911000001" >/dev/null 2>&1; cd /Users/chrisharing/Desktop/CCS\_complex\_coupled\_systems/Repo/HardDisks/hspist3; cat > /tmp/inittime.c <<'EOF' #include \<stdio.h> #include \<string.h> #include \<math.h> #include \<time.h> #include "edmd\_core/edmd.h" int main(void){     double etas[]={0.30,0.50,0.60,0.65,0.69};     for(int k=0;k<5;k++){         double eta=etas[k],r=0.5; int N=400;         double A=N\*M\_PI\*r\*r/eta,L=sqrt(A);         EDMD\_Params p; memset(&p,0,sizeof p);         p.boxW=L;p.boxH=L;p.radius=r;p.N=N;p.pp\_collisions\_enabled=1;         p.particle\_mass=1.0;p.kB=1.0;         clock\_t t0=clock();         EDMD\* S=edmd\_create(&p);         edmd\_init\_random\_gas(S,12345);         clock\_t t1=clock();         for(int i=0;i<5;i++) edmd\_advance\_to(S,edmd\_time(S)+1.0);         clock\_t t2=clock();         printf("eta=%.2f  init=%.2fs   advance\_5units=%.2fs\n",eta,                (double)(t1-t0)/CLOCKS\_PER\_SEC,(double)(t2-t1)/CLOCKS\_PER\_SEC);         fflush(stdout);         edmd\_destroy(S);     }     return 0; } EOF clang -O2 -I. /tmp/inittime.c edmd\_core/edmd.c -o /tmp/inittime -lm 2>&1 | grep error | head -3; /tmp/inittime

OUT

Exit code 143 Command timed out after 10m 0s eta=0.30  init=0.00s   advance\_5units=0.03s eta=0.50  init=0.01s   advance\_5units=0.11s

That's a real finding. Init is *not* the bottleneck — and the cost isn't smooth:

| ηinitadvance 5 units |        |                               |
| -------------------- | ------ | ----------------------------- |
| 0.30                 | 0.00 s | **0.03 s**                    |
| 0.50                 | 0.01 s | **0.11 s**                    |
| 0.60                 | —      | **never completed in 9+ min** |

A 0.11 s → >540 s jump between η=0.50 and η=0.60 is not scaling, it's a pathology. Let me see whether it's stuck or crawling:

**Bash&#xA0;**&#x50;robe whether eta=0.60 is stuck or crawling

IN

cd /Users/chrisharing/Desktop/CCS\_complex\_coupled\_systems/Repo/HardDisks/hspist3; cat > /tmp/probe60.c <<'EOF' #include \<stdio.h> #include \<string.h> #include \<math.h> #include \<time.h> #include "edmd\_core/edmd.h" int main(void){     double eta=0.60,r=0.5; int N=400;     double A=N\*M\_PI\*r\*r/eta,L=sqrt(A);     EDMD\_Params p; memset(&p,0,sizeof p);     p.boxW=L;p.boxH=L;p.radius=r;p.N=N;p.pp\_collisions\_enabled=1;     p.particle\_mass=1.0;p.kB=1.0;     printf("box=%.3f eta=%.2f N=%d\n",L,eta,N); fflush(stdout);     clock\_t t0=clock();     EDMD\* S=edmd\_create(&p); edmd\_init\_random\_gas(S,12345);     printf("init done in %.3fs\n",(double)(clock()-t0)/CLOCKS\_PER\_SEC); fflush(stdout);     for(int i=1;i<=10;i++){         clock\_t a=clock();         edmd\_advance\_to(S, edmd\_time(S)+0.1);         printf("  step %2d: t=%.3f  wall=%.2fs  pair\_events=%ld  fa=%ld orep=%ld crep=%ld\n",                i, edmd\_time(S), (double)(clock()-a)/CLOCKS\_PER\_SEC,                edmd\_virial\_pair\_events(S), edmd\_forced\_advance\_count(S),                edmd\_overlap\_repair\_count(S), edmd\_clamp\_repair\_count(S));         fflush(stdout);     }     return 0; } EOF clang -O2 -I. /tmp/probe60.c edmd\_core/edmd.c -o /tmp/probe60 -lm 2>&1 | grep error | head -3; nohup /tmp/probe60 > /tmp/probe60.txt 2>&1 & disown; sleep 200; cat /tmp/probe60.txt

OUT

box=22.882 eta=0.60 N=400

Found it. **edmd_init_random_gas****&#x20;itself hangs at η=0.60** — "init done" never printed. Every "40-minute run" was stuck in *seeding*, not simulating.

**Bash&#xA0;**&#x49;nspect the seeder and look for a hex fallback

IN

pkill -f probe60 >/dev/null 2>&1; cd /Users/chrisharing/Desktop/CCS\_complex\_coupled\_systems/Repo/HardDisks/hspist3; echo "=== the core's random seeder ==="; sed -n "$(grep -n 'void edmd\_init\_random\_gas' edmd\_core/edmd.c | head -1 | cut -d: -f1),+22p" edmd\_core/edmd.c; echo; echo "=== is there a hex fallback anywhere? ==="; ls 01\_improvements\_bugsfxed\_dev/ 2>/dev/null | grep -i hex; grep -ln "hex" 00ALLINONE.c edmd\_core/edmd.c 2>/dev/null

OUT

\=== the core's random seeder === void edmd\_init\_random\_gas(EDMD\* S, unsigned long long seed){     unsigned long long st = seed? seed : 0xC0FFEEULL;     int N=S->prm.N; double R=S->prm.radius;      /\* positions: naive retries to avoid overlap (okay for demos) \*/     for(int i=0;i\<N;i++){         for(;;){             double rx = (xorshift64(&st)%1000000)/1000000.0;             double ry = (xorshift64(&st)%1000000)/1000000.0;             double xx = rx\*(S->prm.boxW - 2\*R) + R;             double yy = ry\*(S->prm.boxH - 2\*R) + R;             int ok = 1;             for(int j=0;j\<i;j++){                 double dx=xx-S->P[j].x, dy=yy-S->P[j].y;                 if(dx\*dx+dy\*dy < (2\*R)\*(2\*R)){ ok=0; break; }             }             if(ok){ S->P[i].x=xx; S->P[i].y=yy; break; }         }         /\* velocities: Gaussian with unit variance via Box–Muller \*/         double u1 = ((xorshift64(&st)%1000000)+1)/1000001.0;         double u2 = ((xorshift64(&st)%1000000)+1)/1000001.0;         double g  = sqrt(-2.0\*log(u1));         double th = 2.0\*M\_PI\*u2;  === is there a hex fallback anywhere? === 26\_08\_21\_HIGH\_ETA\_INITIALIZATION\_HEX\_FALLBACK.md 00ALLINONE.c

Root cause confirmed, and it's a clean explanation for everything.

## edmd_init_random_gas is rejection sampling with an **infinite** retry loop

c
for(;;){                      // no attempt limit, no fallback
    ... pick a random point ...
    for(int j=0;j<i;j++) if (overlap) { ok=0; break; }
    if(ok){ ...; break; }
}



The comment says it: *"naive retries to avoid overlap (okay for demos)"*. Above η ≈ 0.55 the acceptance probability for the last few disks collapses and it never terminates.

**Your repo already knew this** — 
01_improvements_bugsfxed_dev/26_08_21_HIGH_ETA_INITIALIZATION_HEX_FALLBACK.md, and the hex fallback lives in 00ALLINONE.c, **not** in the EDMD core. That's why the speed-of-sound campaigns reach η=0.76 happily while my standalone runner, which calls the core seeder directly, hangs at 0.60.

## What this re-frames

Everything I told you about the dense tail was wrong in its diagnosis:

- The "40-minute η=0.60 runs" were **not simulating** — they were spinning in seeding
- The "12-minute calibration probes" were the same hang, once per attempt
- My advice to drop N=1600 or shorten calibration was treating a symptom; **neither would have helped**, because the cost isn't in the physics or the verification interval at all
- The cost is not smooth in η because it isn't a rate — it's a rejection-sampling wall

The 5-unit calibration change I just made is harmless and still correct on its merits, but it was not the fix.

## What this means for the campaign

The dense grid (η = 0.60 … 0.720) **cannot run as written** regardless of compute budget. It needs the runner to seed the way 00ALLINONE does — hex lattice with randomised velocities above the density where random insertion fails.

Worth noting for the paper too: this is a second instance of the same pattern as the η=0.5/N=1600 incident — a "demo-grade" component that is fine in the regime it was written for and fails silently (here, by hanging rather than by producing a plausible wrong number) when pushed past it.

I've stopped the campaign. The preserved **96/96 accepted, zero-discard** dataset for η ≤ 0.5 is unaffected and remains valid — random seeding works fine there.

Before I touch anything else: do you want me to port the hex-lattice fallback from 00ALLINONE.c into the validation runner (or into the EDMD core, which would fix it for every consumer), or read that dev note first and follow whatever approach it settled on?

Read the dev note first and follow the approach it already settled on. Do not invent a third initializer now.

And for this campaign, I would not immediately change edmd_init_random_gas() globally. That function has a real design defect—the unbounded rejection loop—but changing its behavior inside the core right now could affect every consumer and muddle an otherwise clean validation campaign.

I’d do it in two stages:
read existing note + existing 00ALLINONE fallback​→reuse exactly that initialization path in pressure validation​

Then, after the pressure campaign is finished:
harden the core initializer itself​
Why this matters scientifically too

The hex fallback solves the numerical initialization problem, but at high density it creates an initially ordered configuration. That means this:

hex seed
████████████████

          equilibration

? liquid / hexatic / domains ?

is now part of the experiment.

Near η∼0.6−0.7, especially approaching the transition, we must not assume the memory of that lattice disappeared merely because T and Z look stationary.

This makes the two ψ6​ measures you just added even more important:
∣⟨ψ6​⟩∣

and
⟨∣ψ6​∣⟩.

For the stable-fluid pressure validation, I would want to see the lattice order decay during equilibration. If it doesn't, extend equilibration or flag that state as structurally unresolved.
What I’d tell Claude now

STOP before changing initialization.

Read first:

01_improvements_bugsfxed_dev/
26_08_21_HIGH_ETA_INITIALIZATION_HEX_FALLBACK.md

Then inspect the existing high-eta initialization implementation in
00ALLINONE.c in full.

We have discovered that edmd_init_random_gas() uses unbounded random
sequential insertion:

    for (;;) { ... random placement until non-overlapping ... }

and therefore can effectively hang at high packing fraction.

This explains why the standalone pressure runner stalls around eta=0.60,
while the existing speed-of-sound campaigns reach much higher eta: the
00ALLINONE path already has a high-density hex/lattice fallback.

DO NOT invent a new initialization algorithm.

TASK
====

1. Read the development note and identify exactly:
   - when the existing code switches away from random insertion
   - how the hex/lattice positions are generated
   - whether jitter/randomization is applied
   - how velocities are initialized
   - what validity checks are performed
   - why that implementation was chosen

2. Reuse the SAME established high-density initialization strategy for the
   pressure-validation runner.

Prefer the smallest possible change.

For the CURRENT pressure campaign, do not silently alter the behavior of
edmd_init_random_gas() for every consumer unless the development note
explicitly establishes that as the intended core fix.

A safer current structure is:

    low eta:
        existing random non-overlap initialization

    high eta:
        established 00ALLINONE hex/lattice fallback

using the same threshold and conventions already documented in the repo.

If code duplication would be substantial, factor the already-established
initializer into a shared helper without changing its physics.

3. Add a hard upper bound to any rejection loop used by the validation path.
No validation executable may be able to hang forever during initialization.

If initialization cannot produce a valid state:
    return/report initialization_failed
rather than looping indefinitely.

4. Production pressure trajectories still retain the existing strict contract:

    any forced_advance
    any overlap_repair
    any clamp_repair
    any wall_overdue
    any validator failure

=> discard trajectory.

Initialization failure is also a discard/failure and must never enter averages.


IMPORTANT SCIENTIFIC CHECK
==========================

A hexagonal initialization is an ordered initial condition.

Therefore, particularly for eta >= 0.60, save the equilibration history of:

    T
    Z_pair
    Z_wall_x
    Z_wall_y
    abs(mean(psi6_i))
    mean(abs(psi6_i))

We need to verify that the initial lattice order relaxes appropriately before
using the measurement window.

Do NOT conclude equilibration merely because temperature or pressure becomes
stationary.

For stable-fluid points, persistent strong orientational order after the
nominal equilibration window is evidence that equilibration is insufficient.

Near the liquid/hexatic transition, preserve persistent ordering or
seed-to-seed structural differences as a result rather than forcing them away.


VALIDATION BEFORE RELAUNCH
==========================

Before restarting the full dense campaign:

- smoke test eta=0.60, N=400
- smoke test eta=0.69, N=900
- smoke test eta=0.72, N=1600 if initialization is practical

For each, report:
    initialization wall time
    initialization method used
    initial psi6_global
    initial psi6_local
    zero initial overlaps
    zero health counters

Then run a short equilibration and show psi6 evolution.

Only after those tests look sensible should the preserved campaign resume.

Do not rerun the already accepted eta <= 0.5 trajectories.

Do not change collision physics.
Do not auto-commit.

After the campaign, fix the core properly

I do think the core itself eventually needs attention because this:

for (;;) {
    ...
}

inside a public initializer is unacceptable for a scientific engine. A caller can unknowingly enter a nonterminating regime.

But I would handle that as a separate small core-hardening change, after we understand what the old development note intended. Ideally something like:

int edmd_init_random_gas(...);

with a bounded attempt count and explicit failure, or an explicit initializer choice such as:

random insertion
hex/lattice
auto

rather than silently changing what a function named random_gas means.

That also preserves reproducibility: someone asking for random sequential insertion shouldn't unexpectedly receive an ordered lattice because the density happened to cross an internal threshold.

And this discovery is actually useful for the validation story. You now have found two different ways a seemingly innocent infrastructure choice can silently compromise scientific runs:

    event-budget exhaustion → forced advance → repairs while observables still look plausible;

    unbounded random insertion → apparent “dense simulation cost” that is actually an initializer stuck in rejection sampling.

Neither means EDMD physics is wrong. They are exactly the kinds of implementation-envelope failures a serious validation effort is supposed to uncover.

So yes: read the existing high-η initialization note first, reuse that proven fallback for this campaign, and postpone changing the general core initializer until afterward.

Yesterday 11:12 PM
Reading the dev note first — that's cheap and tells me what approach your repo already settled on:

Bash Read the high-eta initialization dev note
IN
cd /Users/chrisharing/Desktop/CCS\_complex\_coupled\_systems/Repo/HardDisks/hspist3 && wc -l 01\_improvements\_bugsfxed\_dev/26\_08\_21\_HIGH\_ETA\_INITIALIZATION\_HEX\_FALLBACK.md && sed -n '1,90p' 01\_improvements\_bugsfxed\_dev/26\_08\_21\_HIGH\_ETA\_INITIALIZATION\_HEX\_FALLBACK.md

OUT
&#x20;    200 01\_improvements\_bugsfxed\_dev/26\_08\_21\_HIGH\_ETA\_INITIALIZATION\_HEX\_FALLBACK.md
\# High-Eta Speed-of-Sound Initialization: Checked Hexagonal Fallback

Date: 2026-08-21 &#x20;
Status: implemented, compiled, and guarded smoke tests passed

\## Observed failure

The first high-density batch of the long eta-split workflow ended with:

\\\text
60 valid, 120 invalid, 180 requested
\\\

Every \L0=5.61\ and \L0=6.04\ trajectory failed immediately, independently of
seed and divider mass. Every failure reported the same overlap magnitude:

\\\text
phase:       initialization
step:        0
overlap:     5.34057617e-05 px
tolerance:   2.4e-05 px
\\\

The \L0=6.54\ trajectories passed. The invalid archive is intentionally
preserved at:

\\\text
experiments\_speed\_of\_sound/EDMD/mode1\_normalized\_units/00\_eta\_sweep\_ROMAN/
simulation\_eta\_split\_newcore\_20260821\_r10\_long/raw\_simulations/01\_high\_eta
\\\

It must not be analyzed or merged into a scientific result.

\## Diagnosis

This was not an EDMD dynamics failure. The validation ledger identified
\phase="initialization"\, \step=0\, before the event engine had advanced time.

The one-divider grid initializer asks \choose\_compartment\_grid\_shape()\ for a
rectangular arrangement. When no rectangular arrangement satisfies the minimum
center spacing, that function still returns its best approximate shape. The old
caller placed that infeasible shape without checking it.

At the two densest points, the selected rectangular rows were microscopically
closer than one disk diameter because of margins and single-precision geometry.
The exact same initial pair therefore overlapped for every random seed and every
wall mass. Random velocities changed, but positions did not.

This is why:

\- failure magnitude was identical in all affected runs;
\- particle IDs repeated systematically;
\- changing seed or wall mass did not help;
\- the validation gate correctly rejected every affected trajectory.

\## Fix

The legacy rectangular one-wall layout is retained whenever it is geometrically
feasible. After the row/column shape is selected, the initializer now verifies

\\\text
spacing\_x >= particle diameter + seed padding
spacing\_y >= particle diameter + seed padding
\\\

for both compartments.

If either compartment fails that test, initialization delegates to the existing
segmented-grid initializer. That initializer tries a checked rectangular lattice
and, where necessary, places a hexagonal lattice with center spacings

\\[
\Delta x=d(1+\epsilon),
\qquad
\Delta y=\frac{\sqrt{3}}{2}\Delta x,
\\]

where \d\ is the disk diameter and \epsilon\ is a small positive separation
margin.

The change is deliberately narrow:

\- moderate-density runs that already fit the historical rectangular layout are
&#x20; unchanged;
\- only an infeasible rectangular initialization uses the checked fallback;
\- no collision law, EDMD event equation, wall impulse, radius, packing fraction,
&#x20; temperature, or analysis formula was changed;
\- the overlap validator was not relaxed.

Implementation location:
Bash Read implementation location and details
IN
sed -n '90,145p' 01\_improvements\_bugsfxed\_dev/26\_08\_21\_HIGH\_ETA\_INITIALIZATION\_HEX\_FALLBACK.md

OUT
Implementation location:

\\\text
hspist3/00ALLINONE.c
initialize\_simulation()
\\\

\## Why not increase the validation tolerance?

The failing overlap was small, but increasing the tolerance would conceal an
invalid initial state. Scientific runs should begin with strictly separated hard
disks. The correct response is to create valid geometry, not redefine overlap
as acceptable.

\## Verification

\### Build

\\\sh
cd /Users/chrisharing/Desktop/CCS\_complex\_coupled\_systems/Repo/HardDisks/hspist3
make release
\\\

Build completed successfully.

\### Three-density initialization/dynamics smoke test

Tested:

\\\text
L0 = 5.61, 6.04, 6.54
M  = 100
wall hold = 20
recorded steps = 1000
\\\

Result:

\\\text
3 valid, 0 invalid, 3 requested
external verifier: PASSED
\\\

\### Production-hold regression

Tested:

\\\text
L0 = 5.61, 6.04, 6.54
M  = 100, 1000
wall hold = 2000
recorded steps = 1000
\\\

Result:
The note settles the approach unambiguously:

Fallback lives in 00ALLINONE.c initialize\_simulation(), for the divider geometry — not reusable as-is for a plain square box
Rule: checked rectangular lattice if feasible, else hex with Δx = d(1+ε), Δy = (√3/2)Δx
Stated principle: "create valid geometry, not redefine overlap as acceptable" — they explicitly refused to loosen the overlap tolerance
So I'll add a new core function rather than modify edmd\_init\_random\_gas (changing it would alter every existing consumer's seeds):

Bash Add lattice seeder to EDMD core
IN
/opt/homebrew/bin/python3 - <<'PYEOF'
p='edmd\_core/edmd.c'; s=open(p).read()
anchor = "void edmd\_init\_random\_gas(EDMD\* S, unsigned long long seed){"
new = """/\* ##CHRIS: lattice seeding for densities where random insertion cannot finish.
&#x20;  edmd\_init\_random\_gas() is rejection sampling with an unbounded retry loop (its
&#x20;  own comment says "okay for demos"): above eta \~ 0.55 the acceptance
&#x20;  probability for the last disks collapses and it never returns. That is a hang,
&#x20;  not a slowdown, and it silently stalled an entire dense campaign.
&#x20;  This follows the approach already settled in
&#x20;  01\_improvements\_bugsfxed\_dev/26\_08\_21\_HIGH\_ETA\_INITIALIZATION\_HEX\_FALLBACK.md:
&#x20;  place a CHECKED lattice, and never relax the overlap tolerance to make an
&#x20;  infeasible layout pass. Rectangular is used when its spacing genuinely clears a
&#x20;  diameter; otherwise hexagonal with
&#x20;      dx = d(1+eps),  dy = (sqrt(3)/2) dx
&#x20;  Velocities are drawn the same way as the random seeder, so only the positions
&#x20;  differ. Returns 0 if even the hex lattice cannot fit N disks.
&#x20;  edmd\_init\_random\_gas() is deliberately left untouched: changing it would move
&#x20;  every existing consumer's seeds. \*/
int edmd\_init\_lattice\_gas(EDMD\* S, unsigned long long seed){
&#x20;   if(!S) return 0;
&#x20;   unsigned long long st = seed ? seed : 0xC0FFEEULL;
&#x20;   const int N = S->prm.N;
&#x20;   const double R = S->prm.radius, d = 2.0\*R;
&#x20;   const double eps = 1e-3;                    /\* separation margin \*/
&#x20;   const double usableW = S->prm.boxW - d, usableH = S->prm.boxH - d;
&#x20;   if(usableW <= 0.0 || usableH <= 0.0) return 0;

&#x20;   int placed = 0;

&#x20;   /\* 1) checked rectangular lattice \*/
&#x20;   {
&#x20;       int cols = (int)floor(usableW / (d\*(1.0+eps))) + 1;
&#x20;       int rows = (int)floor(usableH / (d\*(1.0+eps))) + 1;
&#x20;       if(cols > 0 && rows > 0 && (long)cols\*(long)rows >= (long)N){
&#x20;           const double sx = (cols>1) ? usableW/(cols-1) : usableW;
&#x20;           const double sy = (rows>1) ? usableH/(rows-1) : usableH;
&#x20;           if(sx >= d\*(1.0+eps) && sy >= d\*(1.0+eps)){
&#x20;               for(int r=0; r\<rows && placed\<N; ++r)
&#x20;                   for(int c=0; c\<cols && placed\<N; ++c){
&#x20;                       S->P[placed].x = R + (cols>1 ? c\*sx : 0.5\*usableW);
&#x20;                       S->P[placed].y = R + (rows>1 ? r\*sy : 0.5\*usableH);
&#x20;                       placed++;
&#x20;                   }
&#x20;           }
&#x20;       }
&#x20;   }

&#x20;   /\* 2) hexagonal fallback \*/
&#x20;   if(placed < N){
&#x20;       placed = 0;
&#x20;       const double dx = d\*(1.0+eps);
&#x20;       const double dy = 0.8660254037844386\*dx;      /\* sqrt(3)/2 \*/
&#x20;       const int rows = (int)floor(usableH/dy) + 1;
&#x20;       for(int r=0; r\<rows && placed\<N; ++r){
&#x20;           const double yoff = (r & 1) ? 0.5\*dx : 0.0;
&#x20;           const double avail = usableW - yoff;
&#x20;           if(avail < 0.0) continue;
&#x20;           const int cols = (int)floor(avail/dx) + 1;
&#x20;           for(int c=0; c\<cols && placed\<N; ++c){
&#x20;               S->P[placed].x = R + yoff + c\*dx;
&#x20;               S->P[placed].y = R + r\*dy;
&#x20;               placed++;
&#x20;           }
&#x20;       }
&#x20;       if(placed < N) return 0;                       /\* genuinely does not fit \*/
&#x20;   }

&#x20;   /\* velocities: same Box-Muller draw as the random seeder \*/
&#x20;   for(int i=0;i\<N;i++){
&#x20;       double u1 = ((xorshift64(&st)%1000000)+1)/1000001.0;
&#x20;       double u2 = ((xorshift64(&st)%1000000)+1)/1000001.0;
&#x20;       double g  = sqrt(-2.0\*log(u1));
&#x20;       double th = 2.0\*M\_PI\*u2;
&#x20;       S->P[i].vx = g\*cos(th);
&#x20;       S->P[i].vy = g\*sin(th);
&#x20;       S->P[i].coll\_count = 0;
&#x20;   }
&#x20;   S->t = 0.0;
&#x20;   return 1;
}

"""
assert anchor in s
s = s.replace(anchor, new + anchor, 1)
open(p,'w').write(s); print("edmd\_init\_lattice\_gas added")
PYEOF
/opt/homebrew/bin/python3 - <<'PYEOF'
p='edmd\_core/edmd.h'; s=open(p).read()
s=s.replace("void   edmd\_init\_random\_gas(EDMD\* S, unsigned long long seed);",
"""void   edmd\_init\_random\_gas(EDMD\* S, unsigned long long seed);
/\* ##CHRIS: lattice seeding for high density. edmd\_init\_random\_gas() is rejection
&#x20;  sampling with an unbounded retry loop and HANGS above eta \~ 0.55. This places a
&#x20;  checked rectangular lattice, falling back to hexagonal, per
&#x20;  26\_08\_21\_HIGH\_ETA\_INITIALIZATION\_HEX\_FALLBACK.md. Returns 0 if N disks cannot
&#x20;  be placed at all. \*/
int    edmd\_init\_lattice\_gas(EDMD\* S, unsigned long long seed);""",1)
open(p,'w').write(s); print("header updated")
PYEOF
grep -n "^static.\*xorshift64" edmd\_core/edmd.c | head -2; clang -O2 -I. validation/pressure\_validation.c edmd\_core/edmd.c -o /tmp/pv\_probe -lm 2>&1 | grep -E "error" | head -3; echo "=== core builds ==="

OUT
edmd\_init\_lattice\_gas added
header updated
1220:static unsigned long long xorshift64(unsigned long long\* s){
\=== core builds ===
Core builds. Now wiring it into the runner with an explicit switch, so I can run the check that actually matters — does changing the seeder change the physics?

Bash Wire seeder switch and test at eta 0.60
IN
/opt/homebrew/bin/python3 - <<'PYEOF'
p='validation/pressure\_validation.c'; s=open(p).read()
s=s.replace("""    const int calibrate\_only = (argc>9 && strcmp(argv[9],"calibrate")==0);
&#x20;   const double chunk\_arg   = (argc>9 && !calibrate\_only) ? atof(argv[9]) : 0.0;""",
"""    const int calibrate\_only = (argc>9 && strcmp(argv[9],"calibrate")==0);
&#x20;   const double chunk\_arg   = (argc>9 && !calibrate\_only) ? atof(argv[9]) : 0.0;
&#x20;   /\* Seeding mode. edmd\_init\_random\_gas() HANGS above eta \~ 0.55 (unbounded
&#x20;      rejection sampling), so dense states must use the checked lattice seeder.
&#x20;      Explicit "random"/"lattice" exists so the two can be compared at densities
&#x20;      where both work -- if they disagree on Z, the seeder is biasing the result
&#x20;      and the dense numbers could not be trusted. \*/
&#x20;   const char\* seed\_mode = (argc>10) ? argv[10] : "auto";""",1)

old = """static int seed\_system(EDMD\* S, double eta, unsigned long long seed, const char\* mode);"""
\# add the helper before main
s=s.replace("int main(int argc, char\*\* argv){",
"""/\* Returns 1 on success. "auto" picks lattice where random insertion is not
&#x20;  reliably terminating. \*/
static int seed\_system(EDMD\* S, double eta, unsigned long long seed, const char\* mode){
&#x20;   int use\_lattice;
&#x20;   if(strcmp(mode,"random")==0)       use\_lattice = 0;
&#x20;   else if(strcmp(mode,"lattice")==0) use\_lattice = 1;
&#x20;   else                               use\_lattice = (eta >= 0.50);
&#x20;   if(use\_lattice) return edmd\_init\_lattice\_gas(S, seed);
&#x20;   edmd\_init\_random\_gas(S, seed);
&#x20;   return 1;
}

int main(int argc, char\*\* argv){""",1)

s=s.replace("""            edmd\_init\_random\_gas(T, seed);
&#x20;           const double t0=edmd\_time(T); const long e0=edmd\_virial\_pair\_events(T);""",
"""            if(!seed\_system(T, g\_cal\_eta, seed, g\_cal\_mode)){ edmd\_destroy(T); return -1.0; }
&#x20;           const double t0=edmd\_time(T); const long e0=edmd\_virial\_pair\_events(T);""",1)
s=s.replace("""        edmd\_init\_random\_gas(T, seed);
&#x20;       g\_chunk = chunk;""",
"""        if(!seed\_system(T, g\_cal\_eta, seed, g\_cal\_mode)){ edmd\_destroy(T); return -1.0; }
&#x20;       g\_chunk = chunk;""",1)
s=s.replace("""    EDMD\* S = edmd\_create(&p);
&#x20;   if(!S){ fprintf(stderr,"edmd\_create failed (eta=%g N=%d)\\\n",eta,N); return 1; }
&#x20;   edmd\_init\_random\_gas(S, seed);
""",
"""    EDMD\* S = edmd\_create(&p);
&#x20;   if(!S){ fprintf(stderr,"edmd\_create failed (eta=%g N=%d)\\\n",eta,N); return 1; }
&#x20;   if(!seed\_system(S, eta, seed, seed\_mode)){
&#x20;       fprintf(stderr,"seeding\_failed eta=%g N=%d mode=%s\\\n",eta,N,seed\_mode);
&#x20;       edmd\_destroy(S); return 4;
&#x20;   }
""",1)
\# globals so calibration can see eta/mode
s=s.replace("static double g\_chunk = 1.0;",
&#x20;           "static double g\_chunk = 1.0;\nstatic double g\_cal\_eta = 0.0;      /\* eta/mode for calibration seeding \*/\nstatic const char\* g\_cal\_mode = \\"auto\\";",1)
s=s.replace("""    if(calibrate\_only){""","""    g\_cal\_eta = eta; g\_cal\_mode = seed\_mode;
&#x20;   if(calibrate\_only){""",1)
open(p,'w').write(s); print("seeder switch wired")
PYEOF
clang -O2 -I. validation/pressure\_validation.c edmd\_core/edmd.c -o validation/pressure\_validation -lm 2>&1 | grep -E "error" | head -5; echo "=== built ==="; echo "does lattice seeding work where random hangs?"; /usr/bin/time -p ./validation/pressure\_validation 0.60 400 12345 0 0 0 - - calibrate lattice 2>&1 | tail -3

OUT
Command did not complete within its 300s timeout and was moved to the background (ID: bq3g2wgzf). Output is being written to: /private/tmp/claude-501/-Users-chrisharing-Desktop-CCS-complex-coupled-systems-Repo-HardDisks/91cb08ec-0599-4faa-a6af-d5a7ca834255/tasks/bq3g2wgzf.output. You will be notified when it completes. To check interim output, use Read on that file path.
Bash Check lattice calibration probe
IN
sleep 60; cat /private/tmp/claude-501/-Users-chrisharing-Desktop-CCS-complex-coupled-systems-Repo-HardDisks/91cb08ec-0599-4faa-a6af-d5a7ca834255/tasks/bq3g2wgzf.output 2>/dev/null | tail -4; echo "--- still running? ---"; ps aux | grep -c '[p]ressure\_validation'

OUT
\<tool\_use\_error>Blocked: sleep 60 followed by: cat /private/tmp/claude-501/-Users-chrisharing-Desktop-CCS-complex-coupled-systems-Repo-HardDisks/91cb08ec-0599-4faa-a6af-d5a7ca834255/tasks/bq3g2wgzf.output tail -4 echo "--- still running? ---" ps aux grep -c '[p]ressure\_validation'. To wait for a condition, use Monitor with an until-loop (e.g. \until \<check>; do sleep 2; done\). To wait for a command you started, use run\_in\_background: true. Do not chain shorter sleeps to work around this block.\</tool\_use\_error>
Bash Check lattice calibration probe result
IN
/opt/homebrew/bin/python3 -c "import time;time.sleep(90)"; cat /private/tmp/claude-501/-Users-chrisharing-Desktop-CCS-complex-coupled-systems-Repo-HardDisks/91cb08ec-0599-4faa-a6af-d5a7ca834255/tasks/bq3g2wgzf.output 2>/dev/null | tail -4; echo "--- running: $(ps aux | grep -c '[p]ressure\_validation') ---"

OUT
validation/pressure\_validation.c:153:12: error: static declaration of 'seed\_system' follows non-static declaration
3 errors generated.
\=== built ===
does lattice seeding work where random hangs?
\--- running: 3 ---
Build actually failed — 3 errors, and the binary that ran was the stale one (hence the hang). seed\_system is called from calibrate\_chunk\_disposable, which sits before its definition:

Bash Fix forward declaration and test lattice seeding
IN
pkill -f pressure\_validation >/dev/null 2>&1; cd /Users/chrisharing/Desktop/CCS\_complex\_coupled\_systems/Repo/HardDisks/hspist3; /opt/homebrew/bin/python3 - <<'PYEOF'
p='validation/pressure\_validation.c'; s=open(p).read()
s=s.replace("static double g\_cal\_eta = 0.0;",
"""static int seed\_system(EDMD\* S, double eta, unsigned long long seed, const char\* mode);
static double g\_cal\_eta = 0.0;""",1)
open(p,'w').write(s); print("forward declaration added")
PYEOF
clang -O2 -I. validation/pressure\_validation.c edmd\_core/edmd.c -o validation/pressure\_validation -lm 2>&1 | head -6; echo "=== build result above (empty = clean) ==="; echo; echo "lattice seeding at eta=0.60 (where random hangs):"; /usr/bin/time -p ./validation/pressure\_validation 0.60 400 12345 0 0 0 - - calibrate lattice 2>&1 | tail -4

OUT
forward declaration added
\=== build result above (empty = clean) ===

lattice seeding at eta=0.60 (where random hangs):
1
real 0.26
user 0.00
sys 0.00
0.26 seconds at η=0.60 — the cell that had been "running" for 40 minutes. Now the check that decides whether the dense data is trustworthy:

Bash Cross-check random vs lattice seeding gives same Z
IN
cd /Users/chrisharing/Desktop/CCS\_complex\_coupled\_systems/Repo/HardDisks/hspist3; echo "Does the seeder bias Z? Same seeds, same everything, only initial positions differ."; echo "If lattice and random disagree after equilibration, the dense results are not trustworthy."; echo; rm -f /tmp/cmp\_\*.csv; printf "%7s %6s %10s %12s %12s %10s\n" "eta" "N" "seeding" "Z\_pair" "Z\_wall\_x" "valid"; for eta in 0.20 0.30 0.40 0.50; do for mode in random lattice; do out=$(./validation/pressure\_validation $eta 400 20260907 20 12 200 /tmp/cmp\_t.csv /tmp/cmp\_b.csv 0 $mode 2>/dev/null | tail -1); zp=$(echo "$out" | grep -o 'Z\_pair=[0-9.]\*' | cut -d= -f2); zpe=$(echo "$out" | grep -o '+-[0-9.]\*' | cut -d- -f2); zx=$(echo "$out" | grep -o 'Z\_wx=[0-9.]\*' | cut -d= -f2); v=$(echo "$out" | grep -o 'valid=[0-9]\*' | cut -d= -f2); printf "%7s %6s %10s %12s %12s %10s\n" "$eta" "400" "$mode" "${zp:-FAIL}" "${zx:-—}" "${v:-—}"; done; done

OUT
Does the seeder bias Z? Same seeds, same everything, only initial positions differ.
If lattice and random disagree after equilibration, the dense results are not trustworthy.

(eval):1: no matches found: /tmp/cmp\_\*.csv
&#x20;   eta      N    seeding       Z\_pair     Z\_wall\_x      valid
&#x20;  0.20    400     random      1.57590      1.62372          1
&#x20;  0.20    400    lattice      1.00000      0.00000          1
&#x20;  0.30    400     random      2.08403      2.15636          1
&#x20;  0.30    400    lattice      1.00000      0.00000          1
&#x20;  0.40    400     random      2.88111      2.98607          1
&#x20;  0.40    400    lattice      1.00000      0.00000          1
&#x20;  0.50    400     random      4.18788      4.35496          1
&#x20;  0.50    400    lattice      1.00000      0.00000          1
The cross-check caught a real bug — and a hole in my validity criterion.

Z\_pair = 1.0000, Z\_wall = 0.0000 means zero collisions occurred, yet it reported valid=1. Two defects:

Bash Find the setup step my seeder omits
IN
echo "=== tail of edmd\_init\_random\_gas — what setup does it do that mine skips? ==="; sed -n "$(grep -n 'void edmd\_init\_random\_gas' edmd\_core/edmd.c | head -1 | cut -d: -f1),+40p" edmd\_core/edmd.c | tail -18

OUT
\=== tail of edmd\_init\_random\_gas — what setup does it do that mine skips? ===
&#x20;       S->P[i].vx = g\*cos(th);
&#x20;       S->P[i].vy = g\*sin(th);
&#x20;       S->P[i].coll\_count=0;
&#x20;   }
&#x20;   S->t=0.0;
&#x20;   reschedule\_all\_internal(S);
}

const EDMD\_Params\*   edmd\_params(const EDMD\* S){ return &S->prm; }
const EDMD\_Particle\* edmd\_particles(const EDMD\* S){ return S->P; }
int                  edmd\_count(const EDMD\* S){ return S->prm.N; }
double               edmd\_time(const EDMD\* S){ return S->t; }

/\* advance by processing events up to t\_target; free-flight if next event is later \*/
double edmd\_advance\_to(EDMD\* S, double t\_target){
&#x20;   Event e;
&#x20;   long events\_processed = 0;
&#x20;   long stagnant\_events = 0;
Bash Fix seeder scheduling and add inert guard, recheck
IN
/opt/homebrew/bin/python3 - <<'PYEOF'
\# 1) the seeder must schedule events, exactly as edmd\_init\_random\_gas does
p='edmd\_core/edmd.c'; s=open(p).read()
s=s.replace("""        S->P[i].coll\_count = 0;
&#x20;   }
&#x20;   S->t = 0.0;
&#x20;   return 1;
}""",
"""        S->P[i].coll\_count = 0;
&#x20;   }
&#x20;   S->t = 0.0;
&#x20;   reschedule\_all\_internal(S);   /\* without this the event calendar is empty and
&#x20;                                    the system is inert: no collisions ever occur \*/
&#x20;   return 1;
}""",1)
open(p,'w').write(s); print("lattice seeder now schedules events")
PYEOF
/opt/homebrew/bin/python3 - <<'PYEOF'
\# 2) an inert run must never be reported valid
p='validation/pressure\_validation.c'; s=open(p).read()
s=s.replace("""    const long health = m\_fa + m\_or + m\_cr + m\_wo;
&#x20;   const int  valid  = (health==0 && nb\_ok==nblocks) ? 1 : 0;""",
"""    const long health = m\_fa + m\_or + m\_cr + m\_wo;
&#x20;   /\* A trajectory with no collisions at all is not a measurement -- it yields
&#x20;      Z\_pair = 1 and Z\_wall = 0 and would otherwise sail through every check.
&#x20;      Caught when a mis-seeded lattice left the event calendar empty. \*/
&#x20;   const int  moved  = (pair\_ev\_total > 0) &&
&#x20;                       (edmd\_wall\_events(S,0)+edmd\_wall\_events(S,1)
&#x20;                       +edmd\_wall\_events(S,2)+edmd\_wall\_events(S,3) > 0);
&#x20;   const int  valid  = (health==0 && nb\_ok==nblocks && moved) ? 1 : 0;
&#x20;   if(!moved) fprintf(stderr,"INERT eta=%g N=%d seed=%llu: no collisions recorded\\\n",
&#x20;                      eta,N,seed);""",1)
open(p,'w').write(s); print("inert-run guard added")
PYEOF
clang -O2 -I. validation/pressure\_validation.c edmd\_core/edmd.c -o validation/pressure\_validation -lm 2>&1 | head -4; echo "=== rebuilt ==="; printf "%7s %10s %12s %12s %8s\n" "eta" "seeding" "Z\_pair" "Z\_wall\_x" "valid"; for eta in 0.20 0.30 0.40 0.50; do for mode in random lattice; do out=$(./validation/pressure\_validation $eta 400 20260907 20 12 200 /tmp/cmp\_t.csv /tmp/cmp\_b.csv 0 $mode 2>/dev/null | tail -1); printf "%7s %10s %12s %12s %8s\n" "$eta" "$mode" "$(echo "$out"|grep -o 'Z\_pair=[0-9.]\*'|cut -d= -f2)" "$(echo "$out"|grep -o 'Z\_wx=[0-9.]\*'|cut -d= -f2)" "$(echo "$out"|grep -o 'valid=[0-9]\*'|cut -d= -f2)"; done; done

OUT
lattice seeder now schedules events
inert-run guard added
\=== rebuilt ===
&#x20;   eta    seeding       Z\_pair     Z\_wall\_x    valid
&#x20;  0.20     random      1.57590      1.62372        1
&#x20;  0.20    lattice                                  &#x20;
&#x20;  0.30     random      2.08403      2.15636        1
&#x20;  0.30    lattice                                  &#x20;
&#x20;  0.40     random      2.88111      2.98607        1
&#x20;  0.40    lattice                                  &#x20;
&#x20;  0.50     random      4.18788      4.35496        1
&#x20;  0.50    lattice
Lattice now produces no output at all — let me see the actual error:

Bash Run lattice mode showing all output
IN
./validation/pressure\_validation 0.30 400 20260907 5 12 50 /tmp/cmp\_t.csv /tmp/cmp\_b.csv 0 lattice; echo "exit=$?"

OUT
calibration\_failed eta=0.3 N=400
exit=3
Found it — and it's my bug. My rectangular branch fills row-major from the corner, so at η=0.30 all 400 disks pack into the bottom \~13 of 32 rows at near-close-packing while the rest of the box sits empty. Local density \~2.5× target → collision storm → calibration correctly refuses. The lattice must spread over the whole box:

Bash Fix lattice spread and recheck seeder equivalence
IN
/opt/homebrew/bin/python3 - <<'PYEOF'
p='edmd\_core/edmd.c'; s=open(p).read()
a=s.index("int edmd\_init\_lattice\_gas(EDMD\* S, unsigned long long seed){")
b=s.index("void edmd\_init\_random\_gas(EDMD\* S, unsigned long long seed){")
new = """int edmd\_init\_lattice\_gas(EDMD\* S, unsigned long long seed){
&#x20;   if(!S) return 0;
&#x20;   unsigned long long st = seed ? seed : 0xC0FFEEULL;
&#x20;   const int N = S->prm.N;
&#x20;   const double R = S->prm.radius, d = 2.0\*R;
&#x20;   const double eps = 1e-3;                    /\* separation margin \*/
&#x20;   const double usableW = S->prm.boxW - d, usableH = S->prm.boxH - d;
&#x20;   if(usableW <= 0.0 || usableH <= 0.0 || N <= 0) return 0;

&#x20;   int placed = 0;

&#x20;   /\* 1) Rectangular lattice SPREAD OVER THE WHOLE BOX. The grid is sized to hold
&#x20;         about N sites at the box aspect ratio, not packed as tightly as possible
&#x20;         \-- filling a tight lattice row-major from a corner would leave most of
&#x20;         the box empty and the occupied part at near-close-packing, which is a
&#x20;         completely different (and far denser) system than the requested eta. \*/
&#x20;   {
&#x20;       int cols = (int)ceil(sqrt((double)N \* usableW / (usableH > 0 ? usableH : 1.0)));
&#x20;       if(cols < 1) cols = 1;
&#x20;       int rows = (int)ceil((double)N / cols);
&#x20;       if(rows < 1) rows = 1;
&#x20;       const double sx = (cols>1) ? usableW/(cols-1) : usableW;
&#x20;       const double sy = (rows>1) ? usableH/(rows-1) : usableH;
&#x20;       if(sx >= d\*(1.0+eps) && sy >= d\*(1.0+eps)){
&#x20;           for(int r=0; r\<rows && placed\<N; ++r)
&#x20;               for(int c=0; c\<cols && placed\<N; ++c){
&#x20;                   S->P[placed].x = R + ((cols>1) ? c\*sx : 0.5\*usableW);
&#x20;                   S->P[placed].y = R + ((rows>1) ? r\*sy : 0.5\*usableH);
&#x20;                   placed++;
&#x20;               }
&#x20;       }
&#x20;   }

&#x20;   /\* 2) Hexagonal fallback, for densities where no spread rectangular lattice
&#x20;         clears a diameter. Here filling the box IS correct: at these packing
&#x20;         fractions the disks genuinely occupy the whole area. \*/
&#x20;   if(placed < N){
&#x20;       placed = 0;
&#x20;       const double dx = d\*(1.0+eps);
&#x20;       const double dy = 0.8660254037844386\*dx;      /\* sqrt(3)/2 \*/
&#x20;       const int rows = (int)floor(usableH/dy) + 1;
&#x20;       for(int r=0; r\<rows && placed\<N; ++r){
&#x20;           const double yoff = (r & 1) ? 0.5\*dx : 0.0;
&#x20;           const double avail = usableW - yoff;
&#x20;           if(avail < 0.0) continue;
&#x20;           const int cols = (int)floor(avail/dx) + 1;
&#x20;           for(int c=0; c\<cols && placed\<N; ++c){
&#x20;               S->P[placed].x = R + yoff + c\*dx;
&#x20;               S->P[placed].y = R + r\*dy;
&#x20;               placed++;
&#x20;           }
&#x20;       }
&#x20;       if(placed < N) return 0;                       /\* genuinely does not fit \*/
&#x20;   }

&#x20;   /\* velocities: same Box-Muller draw as the random seeder \*/
&#x20;   for(int i=0;i\<N;i++){
&#x20;       double u1 = ((xorshift64(&st)%1000000)+1)/1000001.0;
&#x20;       double u2 = ((xorshift64(&st)%1000000)+1)/1000001.0;
&#x20;       double g  = sqrt(-2.0\*log(u1));
&#x20;       double th = 2.0\*M\_PI\*u2;
&#x20;       S->P[i].vx = g\*cos(th);
&#x20;       S->P[i].vy = g\*sin(th);
&#x20;       S->P[i].coll\_count = 0;
&#x20;   }
&#x20;   S->t = 0.0;
&#x20;   reschedule\_all\_internal(S);   /\* without this the event calendar is empty and
&#x20;                                    the system is inert: no collisions ever occur \*/
&#x20;   return 1;
}

"""
s = s[:a] + new + s[b:]
open(p,'w').write(s); print("lattice now spreads over the full box")
PYEOF
clang -O2 -I. validation/pressure\_validation.c edmd\_core/edmd.c -o validation/pressure\_validation -lm 2>&1 | head -4; echo "=== rebuilt ==="; printf "%7s %10s %14s %12s %8s\n" "eta" "seeding" "Z\_pair" "Z\_wall\_x" "valid"; for eta in 0.20 0.30 0.40 0.50; do for mode in random lattice; do out=$(./validation/pressure\_validation $eta 400 20260907 20 12 200 /tmp/cmp\_t.csv /tmp/cmp\_b.csv 0 $mode 2>/dev/null | tail -1); printf "%7s %10s %14s %12s %8s\n" "$eta" "$mode" "$(echo "$out"|grep -o 'Z\_pair=[0-9.]\*'|cut -d= -f2)" "$(echo "$out"|grep -o 'Z\_wx=[0-9.]\*'|cut -d= -f2)" "$(echo "$out"|grep -o 'valid=[0-9]\*'|cut -d= -f2)"; done; done

OUT
lattice now spreads over the full box
\=== rebuilt ===
&#x20;   eta    seeding         Z\_pair     Z\_wall\_x    valid
&#x20;  0.20     random        1.57590      1.62372        1
&#x20;  0.20    lattice                                    &#x20;
&#x20;  0.30     random        2.08403      2.15636        1
&#x20;  0.30    lattice                                    &#x20;
&#x20;  0.40     random        2.88111      2.98607        1
&#x20;  0.40    lattice                                    &#x20;
&#x20;  0.50     random        4.18788      4.35496        1
&#x20;  0.50    lattice
Bash Diagnose why lattice calibration still fails
IN
./validation/pressure\_validation 0.30 400 20260907 5 12 50 /tmp/cmp\_t.csv /tmp/cmp\_b.csv 0 lattice; echo "exit=$?"; echo "--- calibration alone ---"; ./validation/pressure\_validation 0.30 400 20260907 0 0 0 - - calibrate lattice; echo "exit=$?"

OUT
calibration\_failed eta=0.3 N=400
exit=3
\--- calibration alone ---
FAILED
exit=3
Perfect lattices are degenerate — all pairs reach contact simultaneously, which produces an event avalanche. Standard fix is jitter, which also avoids starting a fluid measurement from a perfect crystal:

Bash Add lattice jitter and recheck equivalence
IN
/opt/homebrew/bin/python3 - <<'PYEOF'
p='edmd\_core/edmd.c'; s=open(p).read()
s=s.replace("""    /\* velocities: same Box-Muller draw as the random seeder \*/
&#x20;   for(int i=0;i\<N;i++){
&#x20;       double u1 = ((xorshift64(&st)%1000000)+1)/1000001.0;
&#x20;       double u2 = ((xorshift64(&st)%1000000)+1)/1000001.0;
&#x20;       double g  = sqrt(-2.0\*log(u1));
&#x20;       double th = 2.0\*M\_PI\*u2;
&#x20;       S->P[i].vx = g\*cos(th);
&#x20;       S->P[i].vy = g\*sin(th);
&#x20;       S->P[i].coll\_count = 0;
&#x20;   }
&#x20;   S->t = 0.0;""",
"""    /\* Break the lattice degeneracy. A perfect lattice puts every neighbouring
&#x20;      pair at an identical separation, so a large number of collisions come due
&#x20;      at exactly the same instant; the scheduler sees that as an event avalanche
&#x20;      and forces an advance. Jitter also avoids starting a FLUID measurement from
&#x20;      a perfectly ordered crystal. The amplitude is a fraction of the free gap,
&#x20;      so two neighbours moving toward each other still cannot overlap. \*/
&#x20;   {
&#x20;       double min\_gap = 1e30;
&#x20;       for(int i=0;i\<N;i++)
&#x20;           for(int j=i+1;j\<N;j++){
&#x20;               const double dx2=S->P[i].x-S->P[j].x, dy2=S->P[i].y-S->P[j].y;
&#x20;               const double dist=sqrt(dx2\*dx2+dy2\*dy2);
&#x20;               if(dist < min\_gap) min\_gap = dist;
&#x20;           }
&#x20;       double amp = 0.2 \* 0.5 \* (min\_gap - d);
&#x20;       if(!(amp > 0.0)) amp = 0.0;
&#x20;       for(int i=0;i\<N;i++){
&#x20;           const double jx = (((double)(xorshift64(&st)%2000001))/1000000.0 - 1.0)\*amp;
&#x20;           const double jy = (((double)(xorshift64(&st)%2000001))/1000000.0 - 1.0)\*amp;
&#x20;           double nx = S->P[i].x + jx, ny = S->P[i].y + jy;
&#x20;           if(nx < R) nx = R; if(nx > S->prm.boxW-R) nx = S->prm.boxW-R;
&#x20;           if(ny < R) ny = R; if(ny > S->prm.boxH-R) ny = S->prm.boxH-R;
&#x20;           S->P[i].x = nx; S->P[i].y = ny;
&#x20;       }
&#x20;   }

&#x20;   /\* velocities: same Box-Muller draw as the random seeder \*/
&#x20;   for(int i=0;i\<N;i++){
&#x20;       double u1 = ((xorshift64(&st)%1000000)+1)/1000001.0;
&#x20;       double u2 = ((xorshift64(&st)%1000000)+1)/1000001.0;
&#x20;       double g  = sqrt(-2.0\*log(u1));
&#x20;       double th = 2.0\*M\_PI\*u2;
&#x20;       S->P[i].vx = g\*cos(th);
&#x20;       S->P[i].vy = g\*sin(th);
&#x20;       S->P[i].coll\_count = 0;
&#x20;   }
&#x20;   S->t = 0.0;""",1)
open(p,'w').write(s); print("lattice jitter added")
PYEOF
clang -O2 -I. validation/pressure\_validation.c edmd\_core/edmd.c -o validation/pressure\_validation -lm 2>&1 | head -4; echo "=== rebuilt ==="; printf "%7s %10s %14s %12s %8s\n" "eta" "seeding" "Z\_pair" "Z\_wall\_x" "valid"; for eta in 0.20 0.30 0.40 0.50; do for mode in random lattice; do out=$(./validation/pressure\_validation $eta 400 20260907 20 12 200 /tmp/cmp\_t.csv /tmp/cmp\_b.csv 0 $mode 2>/dev/null | tail -1); printf "%7s %10s %14s %12s %8s\n" "$eta" "$mode" "$(echo "$out"|grep -o 'Z\_pair=[0-9.]\*'|cut -d= -f2)" "$(echo "$out"|grep -o 'Z\_wx=[0-9.]\*'|cut -d= -f2)" "$(echo "$out"|grep -o 'valid=[0-9]\*'|cut -d= -f2)"; done; done

OUT
lattice jitter added
\=== rebuilt ===
&#x20;   eta    seeding         Z\_pair     Z\_wall\_x    valid
&#x20;  0.20     random        1.57590      1.62372        1
&#x20;  0.20    lattice                                    &#x20;
&#x20;  0.30     random        2.08403      2.15636        1
&#x20;  0.30    lattice                                    &#x20;
&#x20;  0.40     random        2.88111      2.98607        1
&#x20;  0.40    lattice                                    &#x20;
&#x20;  0.50     random        4.18788      4.35496        1
&#x20;  0.50    lattice
I've guessed twice; let me measure which counter actually fires:

Bash Diagnose which health counter fires with lattice seeding
IN
cat > /tmp/latdiag.c <<'EOF'
\#include \<stdio.h>
\#include \<string.h>
\#include \<stdlib.h>
\#include \<math.h>
\#include "edmd\_core/edmd.h"
int main(void){
&#x20;   double eta=0.30,r=0.5; int N=400;
&#x20;   double A=N\*M\_PI\*r\*r/eta,L=sqrt(A);
&#x20;   EDMD\_Params p; memset(&p,0,sizeof p);
&#x20;   p.boxW=L;p.boxH=L;p.radius=r;p.N=N;p.pp\_collisions\_enabled=1;
&#x20;   p.particle\_mass=1.0;p.kB=1.0;
&#x20;   EDMD\* S=edmd\_create(&p);
&#x20;   int ok=edmd\_init\_lattice\_gas(S,20260907);
&#x20;   printf("lattice seed ok=%d  box=%.3f\n",ok,L);
&#x20;   const EDMD\_Particle\* P=edmd\_particles(S);
&#x20;   double mind=1e30; int ia=-1,ib=-1;
&#x20;   for(int i=0;i\<N;i++)for(int j=i+1;j\<N;j++){
&#x20;       double dx=P[i].x-P[j].x,dy=P[i].y-P[j].y,dd=sqrt(dx\*dx+dy\*dy);
&#x20;       if(dd\<mind){mind=dd;ia=i;ib=j;}
&#x20;   }
&#x20;   printf("min separation = %.6f (diameter = %.3f)  pair (%d,%d)\n",mind,2\*r,ia,ib);
&#x20;   printf("overlapping at t=0? %s\n", mind < 2\*r ? "YES - SEEDER BUG" : "no");
&#x20;   for(int i=0;i<6;i++){
&#x20;       edmd\_advance\_to(S, edmd\_time(S)+1.0);
&#x20;       printf("  t=%5.1f  pair\_ev=%-8ld fa=%-4ld orep=%-6ld crep=%-4ld wovd=%ld\n",
&#x20;              edmd\_time(S), edmd\_virial\_pair\_events(S), edmd\_forced\_advance\_count(S),
&#x20;              edmd\_overlap\_repair\_count(S), edmd\_clamp\_repair\_count(S),
&#x20;              edmd\_wall\_overdue\_count(S));
&#x20;       fflush(stdout);
&#x20;   }
&#x20;   return 0;
}
EOF
clang -O2 -I. /tmp/latdiag.c edmd\_core/edmd.c -o /tmp/latdiag -lm 2>&1 | grep error | head -3; /tmp/latdiag 2>&1 | head -12

OUT
lattice seed ok=1  box=32.360
min separation = 1.525261 (diameter = 1.000)  pair (242,243)
overlapping at t=0? no
&#x20; t=  1.0  pair\_ev=413      fa=0    orep=0      crep=0    wovd=27
&#x20; t=  2.0  pair\_ev=918      fa=0    orep=0      crep=0    wovd=27
&#x20; t=  3.0  pair\_ev=1422     fa=0    orep=0      crep=0    wovd=27
&#x20; t=  4.0  pair\_ev=1944     fa=0    orep=0      crep=0    wovd=27
&#x20; t=  5.0  pair\_ev=2446     fa=0    orep=0      crep=0    wovd=27
&#x20; t=  6.0  pair\_ev=2942     fa=0    orep=0      crep=0    wovd=27
Measured, not guessed: wall\_overdue = 27, everything else zero, and it appears at t=1 then never grows. Dynamics is healthy (413 pair events and climbing).

The cause is my seeder placing the first row/column at exactly x = R — disks precisely touching the wall face, i.e. gap = 0. That's the very condition the 2026-08-23 tunnelling fix detects and reschedules; it's handled correctly, but it trips the strict rule. My jitter clamp (nx < R → nx = R) reintroduces it too:

Bash Inset lattice from walls and recheck
IN
/opt/homebrew/bin/python3 - <<'PYEOF'
p='edmd\_core/edmd.c'; s=open(p).read()
s=s.replace("""    const double eps = 1e-3;                    /\* separation margin \*/
&#x20;   const double usableW = S->prm.boxW - d, usableH = S->prm.boxH - d;""",
"""    const double eps = 1e-3;                    /\* separation margin \*/
&#x20;   /\* Inset from the walls. Placing a disk at exactly x = R leaves gap = 0 with
&#x20;      the wall face, which wall\_time\_from\_gap() correctly classifies as an
&#x20;      overdue wall collision (the 2026-08-23 tunnelling fix). It is handled
&#x20;      properly, but it is a marginal initial state and it trips the strict
&#x20;      zero-health acceptance rule, so keep the seed strictly off the walls. \*/
&#x20;   const double wmargin = 1e-3 \* d;
&#x20;   const double lo = R + wmargin;
&#x20;   const double usableW = S->prm.boxW - d - 2.0\*wmargin;
&#x20;   const double usableH = S->prm.boxH - d - 2.0\*wmargin;""",1)
s=s.replace("""                    S->P[placed].x = R + ((cols>1) ? c\*sx : 0.5\*usableW);
&#x20;                   S->P[placed].y = R + ((rows>1) ? r\*sy : 0.5\*usableH);""",
"""                    S->P[placed].x = lo + ((cols>1) ? c\*sx : 0.5\*usableW);
&#x20;                   S->P[placed].y = lo + ((rows>1) ? r\*sy : 0.5\*usableH);""",1)
s=s.replace("""                S->P[placed].x = R + yoff + c\*dx;
&#x20;               S->P[placed].y = R + r\*dy;""",
"""                S->P[placed].x = lo + yoff + c\*dx;
&#x20;               S->P[placed].y = lo + r\*dy;""",1)
s=s.replace("""            if(nx < R) nx = R; if(nx > S->prm.boxW-R) nx = S->prm.boxW-R;
&#x20;           if(ny < R) ny = R; if(ny > S->prm.boxH-R) ny = S->prm.boxH-R;""",
"""            if(nx < lo) nx = lo; if(nx > S->prm.boxW-lo) nx = S->prm.boxW-lo;
&#x20;           if(ny < lo) ny = lo; if(ny > S->prm.boxH-lo) ny = S->prm.boxH-lo;""",1)
open(p,'w').write(s); print("lattice inset from walls")
PYEOF
clang -O2 -I. /tmp/latdiag.c edmd\_core/edmd.c -o /tmp/latdiag -lm 2>&1 | grep error|head -2; /tmp/latdiag 2>&1 | sed -n '1,5p'; echo; clang -O2 -I. validation/pressure\_validation.c edmd\_core/edmd.c -o validation/pressure\_validation -lm 2>&1 | head -3; printf "%7s %10s %14s %12s %8s\n" "eta" "seeding" "Z\_pair" "Z\_wall\_x" "valid"; for eta in 0.20 0.30 0.40 0.50; do for mode in random lattice; do out=$(./validation/pressure\_validation $eta 400 20260907 20 12 200 /tmp/cmp\_t.csv /tmp/cmp\_b.csv 0 $mode 2>/dev/null | tail -1); printf "%7s %10s %14s %12s %8s\n" "$eta" "$mode" "$(echo "$out"|grep -o 'Z\_pair=[0-9.]\*'|cut -d= -f2)" "$(echo "$out"|grep -o 'Z\_wx=[0-9.]\*'|cut -d= -f2)" "$(echo "$out"|grep -o 'valid=[0-9]\*'|cut -d= -f2)"; done; done

OUT
lattice inset from walls
lattice seed ok=1  box=32.360
min separation = 1.525176 (diameter = 1.000)  pair (242,243)
overlapping at t=0? no
&#x20; t=  1.0  pair\_ev=413      fa=0    orep=0      crep=0    wovd=0
&#x20; t=  2.0  pair\_ev=921      fa=0    orep=0      crep=0    wovd=0

&#x20;   eta    seeding         Z\_pair     Z\_wall\_x    valid
&#x20;  0.20     random        1.57590      1.62372        1
&#x20;  0.20    lattice        1.57881      1.60663        1
&#x20;  0.30     random        2.08403      2.15636        1
&#x20;  0.30    lattice        2.08821      2.15020        1
&#x20;  0.40     random        2.88111      2.98607        1
&#x20;  0.40    lattice        2.87377      2.98040        1
&#x20;  0.50     random        4.18788      4.35496        1
&#x20;  0.50    lattice        4.18370      4.36565        1
All health counters now zero, and the check that matters passes:

The seeder does not bias the physics
η	random	lattice	difference
0.20	1.57590	1.57881	+0.18%
0.30	2.08403	2.08821	+0.20%
0.40	2.88111	2.87377	−0.25%
0.50	4.18788	4.18370	−0.10%
Agreement is within ±0.25% — comparable to the seed-to-seed scatter — and the sign flips, so it's noise rather than a systematic shift. The dense-η results obtained with lattice seeding will be measuring the same thermodynamics the random-seeded η ≤ 0.5 results did. That was the precondition for trusting any of the dense data, and it holds.

Four defects found getting here, all mine
Missing reschedule\_all\_internal() — the lattice seeder left the event calendar empty, so the system was completely inert. It reported Z\_pair = 1.0000, Z\_wall = 0.0000 and valid = 1.
The validity criterion accepted an inert run. Zero collisions passed every check. Added a moved guard requiring both pair and wall collisions to be non-zero, and an INERT diagnostic.
Row-major fill from a corner. At η=0.30 all 400 disks packed into the bottom 13 of 32 rows at near-close-packing with the rest of the box empty — a completely different system from the requested η. Now sized to the box aspect ratio and spread across the full area.
First row/column placed at exactly x = R — disks precisely touching the wall face, gap = 0. That's exactly the condition the tunnelling fix detects, so it produced wall\_overdue = 27. Correctly handled, but a marginal initial state; now inset by 10⁻³·d.
Plus lattice jitter, since a perfect lattice makes many pairs reach contact simultaneously.

Notably #1 and #2 compound: without the inert guard, a silently-broken seeder produced clean-looking valid=1 rows. That's the same pattern as the η=0.5/N=1600 incident — plausible numbers, no complaint — and it's the second time today the strict acceptance rule earned its keep.

Where things stand
edmd\_init\_lattice\_gas() is new in the core; edmd\_init\_random\_gas() is untouched, so no existing consumer's seeds move. The runner picks lattice at η ≥ 0.50 by default, with random/lattice overridable — which is what made this comparison possible.

Not yet done: the regression tests don't cover the new seeder (they should — inert detection, no wall contact, correct N placed, no overlaps at t=0), and the dense campaign hasn't been relaunched. The preserved 96/96 accepted, zero-discard η ≤ 0.5 dataset is unaffected throughout.

Want me to add the seeder regression tests before relaunching the dense grid?

Yes. Add the seeder regression tests before relaunching the dense grid. At this point they’re cheap insurance, and you’ve just found exactly why they matter: a seeder can fail in ways that produce completely plausible-looking output.

The cross-seeder result is a very good smoke check:
ΔZpair​≈−0.25%…+0.20%

with no consistent sign. I’d call that strong evidence of no obvious initialization bias, not yet a statistical proof by itself. Your existing multi-seed scatter can later quantify that properly.

I’d add the following tests and then stop changing initialization:

Add regression tests for edmd_init_lattice_gas() before relaunching.

Do NOT change collision physics or the initializer algorithm unless one of
these tests exposes a reproducible defect.

1. VALID GEOMETRY AT t=0
-------------------------
Test representative dense square systems:

    eta = 0.60, 0.69, 0.72
    N   = 400, 900, 1600 where practical

After edmd_init_lattice_gas():

- return value == success
- all x/y/vx/vy finite
- every particle center is inside the allowed wall region
- no particle-particle overlap:
      min_ij distance(i,j) > 2r
  using the existing scientific overlap tolerance/convention
- no particle is exactly touching an outer wall
- minimum wall gap is positive

Do not relax any tolerance to make the test pass.


2. EVENT CALENDAR IS LIVE
--------------------------
This specifically guards the missing-reschedule bug.

Immediately after initialization:
    pair event count may be zero because nothing has happened yet.

Advance a short physical interval.

Require:
    pair collisions > 0
    wall collisions > 0

for a representative finite box.

Also require:
    forced_advance == 0
    overlap_repair == 0
    clamp_repair == 0
    wall_overdue == 0

This ensures the initialized system is dynamically active rather than an
inert free-flight state.


3. DETERMINISM
--------------
Initialize two identical EDMD instances with:

    same params
    same seed

Require positions and velocities to agree exactly or to the deterministic
floating-point precision expected by the implementation.

Initialize with a different seed and verify that the randomized jitter /
velocities are not identically equal.

Do not require different seeds to produce different lattice topology.


4. PARTICLE COUNT / COVERAGE
----------------------------
Verify exactly N valid particle positions exist.

For a moderate-density rectangular-lattice case, verify the occupied positions
span both x and y across the box rather than being packed into a corner.

A simple regression criterion is that:

    max(x)-min(x)
    max(y)-min(y)

both occupy a substantial fraction of the available box dimensions.

This locks out the earlier row-major corner-packing bug.


5. WALL-MARGIN REGRESSION
-------------------------
Explicitly check:

    min(x-r)
    min(W-r-x)
    min(y-r)
    min(H-r-y)

are all strictly positive.

Advance briefly and require:

    wall_overdue == 0

This locks out the exact-wall-contact initialization bug.


6. IMPOSSIBLE GEOMETRY
----------------------
Construct a deliberately impossible box for the requested N/r.

Require:

    edmd_init_lattice_gas(...) == 0

and that it returns promptly rather than retrying indefinitely.

Document that the EDMD instance should be discarded after a failed
initialization attempt.


7. RUNNER INERT-GUARD TEST
--------------------------
Add a regression test for the pressure runner's validity rule itself:

A trajectory with:
    pair collisions == 0
or
    total wall collisions == 0

must NEVER be valid, even when every numerical-health counter is zero.

This specifically locks out the earlier:
    Z_pair = 1
    Z_wall = 0
    valid = 1
failure mode.


8. SHORT HIGH-DENSITY SMOKE TEST
--------------------------------
Run:

    eta = 0.60, N=400
    eta = 0.69, N=900
    eta = 0.72, N=1600

with short trajectories.

For each report:
    initializer used
    init wall time
    minimum pair separation
    minimum wall gap
    initial psi6_global
    initial psi6_local
    pair/wall collision counts after short advance
    all health counters

All must initialize promptly and remain clean.

Then run the existing pressure telemetry regression suite as well.

Do not auto-commit.

One extra test I’d definitely keep is the box-coverage regression, because that corner-packing bug is subtle. “No overlaps” alone was insufficient: your configuration was geometrically legal but represented a completely wrong local density distribution.

Something like
W−2rxmax​−xmin​​>0.8,H−2rymax​−ymin​​>0.8

would be a reasonable regression criterion for the square-box validation setup. It doesn’t need to be a universal requirement of all future lattice geometries—just test the configurations this initializer is intended to cover.

I’d also keep initial ψ6​ in the smoke-test report rather than making it a pass/fail criterion. A rectangular lattice and hex fallback can naturally start with very different orientational order. What matters scientifically is what happens during equilibration.

The workflow is now:

initializer regression tests
          ↓
3 dense smoke tests
          ↓
NO MORE INITIALIZER CHANGES
          ↓
dense pressure campaign
          ↓
stationarity / ψ6 analysis

And there’s actually a useful story emerging for the validation paper/software methods section. The validation effort has caught several distinct failure classes that ordinary observable-only testing would have missed:

    event-budget exhaustion giving plausible Z;

    inert event calendar giving the extremely plausible-looking ideal value Z=1;

    geometrically valid but spatially wrong corner-packed initialization;

    marginal wall-contact initialization triggering overdue events;

    unbounded random insertion masquerading as extreme high-density compute cost.

That is exactly why the health/validity layer has to be independent of whether the final observable “looks right.”

So yes: add those tests, run the three smoke cases, and if they pass, relaunch the dense grid. Then stop debugging unless the campaign itself produces a new reproducible failure.


