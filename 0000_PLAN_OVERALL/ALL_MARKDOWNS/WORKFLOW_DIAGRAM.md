# Workflow Diagram: Hard Disk Simulation & Analysis

## Complete Pipeline

```
┌─────────────────────────────────────────────────────────────────────┐
│                    1. SIMULATION (00ALLINONE)                       │
└─────────────────────────────────────────────────────────────────────┘
                                  │
                                  │ --mode=edmd --particles=100 --l0=20
                                  │ --wall-mass-factor=200 --kbt1
                                  │ --auto-release --steps=500000
                                  │
                                  ▼
                    ┌──────────────────────────┐
                    │  Event-Driven MD Engine  │
                    │  - Hard disk collisions  │
                    │  - Wall oscillations     │
                    │  - Auto-release at t=0   │
                    └──────────────────────────┘
                                  │
                    ┌─────────────┴─────────────┐
                    │                           │
                    ▼                           ▼
        ┌───────────────────┐      ┌───────────────────────┐
        │ wall_position.csv │      │ particle_states_*.csv │
        │ - Time            │      │ - All particle data   │
        │ - Wall x-position │      │ - Position, velocity  │
        └───────────────────┘      └───────────────────────┘
                    │
                    │
┌───────────────────┴───────────────────────────────────────────────┐
│                    2. ANALYSIS (wall_x_FFT.py)                    │
└───────────────────────────────────────────────────────────────────┘
                    │
                    ▼
    ┌───────────────────────────────────────┐
    │ Step 1: Remove Flat Initial Region    │
    │ - Find where wall starts moving       │
    │ - Trim t < t_wall_release             │
    └───────────────────────────────────────┘
                    │
                    ▼
    ┌───────────────────────────────────────┐
    │ Step 2: Entropy-Based Detection       │
    │ ┌─────────────────────────────────┐   │
    │ │ Calculate S[p(x)] in windows    │   │
    │ │ - Histogram positions           │   │
    │ │ - S = -Σ p log p                │   │
    │ └─────────────────────────────────┘   │
    │ ┌─────────────────────────────────┐   │
    │ │ Calculate S[p(v)] in windows    │   │
    │ │ - Histogram velocities          │   │
    │ │ - S = -Σ p log p                │   │
    │ └─────────────────────────────────┘   │
    │ ┌─────────────────────────────────┐   │
    │ │ Detect stabilization            │   │
    │ │ - std(S)/mean(S) < 5%?          │   │
    │ │ - Use LATER of position/velocity│   │
    │ └─────────────────────────────────┘   │
    └───────────────────────────────────────┘
                    │
                    ▼
    ┌───────────────────────────────────────┐
    │ Step 3: Remove Transients             │
    │ - Trim first ~10-20% (adaptive!)      │
    │ - Keep only steady-state data         │
    └───────────────────────────────────────┘
                    │
                    ▼
    ┌───────────────────────────────────────┐
    │ Step 4: FFT Analysis                  │
    │ - Remove mean (DC component)          │
    │ - Apply FFT                           │
    │ - Find peak frequency                 │
    └───────────────────────────────────────┘
                    │
                    ▼
    ┌───────────────────────────────────────┐
    │ Step 5: Theory Validation             │
    │ ┌─────────────────────────────────┐   │
    │ │ Calculate theoretical ν         │   │
    │ │ - η = Nπσ²/(4AL₀)               │   │
    │ │ - cs from Henderson equation    │   │
    │ │ - K from cot(K) = (M/2Nm)K      │   │
    │ │ - ν_theory = (cs/2πL₀)K         │   │
    │ └─────────────────────────────────┘   │
    │ ┌─────────────────────────────────┐   │
    │ │ Compare ν_sim vs ν_theory       │   │
    │ │ - Ratio = ν_sim / ν_theory      │   │
    │ │ - Verdict: Excellent/Good/...   │   │
    │ └─────────────────────────────────┘   │
    └───────────────────────────────────────┘
                    │
                    ▼
    ┌───────────────────────────────────────┐
    │ Step 6: Generate Plots                │
    │ ┌─────────────────────────────────┐   │
    │ │ 3-Panel Entropy Plot            │   │
    │ │ - Top: Displacement             │   │
    │ │ - Middle: S[p(x)] over time     │   │
    │ │ - Bottom: S[p(v)] over time     │   │
    │ │ - Red line: steady-state point  │   │
    │ └─────────────────────────────────┘   │
    │ ┌─────────────────────────────────┐   │
    │ │ FFT Power Spectrum              │   │
    │ │ - Power vs frequency            │   │
    │ │ - Mark fundamental peak         │   │
    │ └─────────────────────────────────┘   │
    └───────────────────────────────────────┘
                    │
                    ▼
        ┌───────────────────────┐
        │    OUTPUT FILES       │
        │ ┌─────────────────┐   │
        │ │ PDF Plots       │   │
        │ │ PNG Plots       │   │
        │ │ Terminal Report │   │
        │ └─────────────────┘   │
        └───────────────────────┘
```

---

## Data Flow Detail

### Input Parameters → Simulation

```
User specifies:
  N = 100 particles (50 per side)
  L₀ = 20σ (box half-length)
  M = 200m (wall mass)
  T = 1000 (initial temperature)
  --kbt1 (scale to kBT=1)
  --auto-release (release wall at t=0)

         ↓

Simulation calculates:
  η = 0.1963 (packing fraction)
  Particle collisions (event-driven)
  Wall oscillations (Newton's 2nd law)

         ↓

Outputs:
  wall_position.csv (every 0.01 time units)
  particle_states_*.csv (periodic snapshots)
```

### Raw Data → Cleaned Data

```
Raw wall_position.csv:
┌────────────────────────────────┐
│ Flat region (wall held)        │ ← Removed in Step 1
├────────────────────────────────┤
│ Transient oscillations         │ ← Removed in Step 2
│ (energy equilibrating)         │    (entropy detection)
├────────────────────────────────┤
│ Steady-state oscillations      │ ← Used for FFT
│ (clean sinusoidal)             │
└────────────────────────────────┘
```

### Entropy Evolution

```
S[p(x)] over time:
│
│      ╱─────────────────  ← Plateau (equilibrium)
│    ╱                       std(S)/mean(S) < 5%
│  ╱                         → Steady state detected
│╱
└──────────────────────> time
0    20   40   60   80

S[p(v)] over time:
│
│          ╱───────────────  ← Plateau (thermal equilibrium)
│        ╱                    Takes longer than spatial!
│      ╱                      std(S)/mean(S) < 5%
│    ╱
│  ╱
│╱
└──────────────────────> time
0    20   40   60   80  100

Use the LATER of the two (conservative)
```

---

## Theory Validation Flow

```
┌─────────────────────────────────────────────────────────────┐
│                   THEORETICAL PREDICTION                    │
└─────────────────────────────────────────────────────────────┘
                           │
        ┌──────────────────┼──────────────────┐
        │                  │                  │
        ▼                  ▼                  ▼
┌──────────────┐  ┌──────────────┐  ┌──────────────┐
│ Packing      │  │ Speed of     │  │ Transcendental│
│ Fraction     │  │ Sound        │  │ Equation     │
│              │  │              │  │              │
│ η = Nπσ²     │  │ cs = √(2kBT/m)│  │ cot(K) =     │
│     4AL₀     │  │   × √(EOS)   │  │ (M/2Nm)K     │
│              │  │              │  │              │
│ η = 0.1963   │  │ cs = 2.254   │  │ K = 0.653    │
└──────────────┘  └──────────────┘  └──────────────┘
        │                  │                  │
        └──────────────────┼──────────────────┘
                           ▼
                  ┌─────────────────┐
                  │ ν_theory        │
                  │ = (cs/2πL₀)K    │
                  │ = 0.0117        │
                  └─────────────────┘
                           │
┌──────────────────────────┼──────────────────────────┐
│                          ▼                          │
│              ┌─────────────────────┐                │
│              │ ν_simulation        │                │
│              │ (from FFT peak)     │                │
│              │ = 0.01165           │                │
│              └─────────────────────┘                │
│                          │                          │
│                          ▼                          │
│              ┌─────────────────────┐                │
│              │ Ratio = ν_sim       │                │
│              │         ν_theory    │                │
│              │ = 0.9946            │                │
│              └─────────────────────┘                │
│                          │                          │
│         ┌────────────────┼────────────────┐         │
│         │                │                │         │
│         ▼                ▼                ▼         │
│   ┌─────────┐    ┌──────────┐    ┌───────────┐    │
│   │ 95-105% │    │ 90-110%  │    │ 80-120%   │    │
│   │ ✅      │    │ ✓        │    │ ⚠️        │    │
│   │EXCELLENT│    │ GOOD     │    │ACCEPTABLE │    │
│   └─────────┘    └──────────┘    └───────────┘    │
└───────────────────────────────────────────────────┘
```

---

## Decision Tree: Is My Run Good?

```
START: Run simulation and analysis
  │
  ▼
Check entropy plots
  │
  ├─ Both plateau? ─── NO ──→ ❌ Run too short
  │                            Fix: --steps=1000000
  │
  └─ YES
      │
      ▼
Check transient removal
  │
  ├─ ~10-20%? ────────── NO ──→ ⚠️ Unusual but check plots
  │                             Might be OK if entropies stable
  │
  └─ YES
      │
      ▼
Check FFT peak
  │
  ├─ Clear peak? ────── NO ──→ ❌ Noisy data
  │                            Fix: run longer or check kBT=1
  │
  └─ YES
      │
      ▼
Check theory ratio
  │
  ├─ 0.95-1.05? ────── YES ──→ ✅ EXCELLENT! You're done!
  │
  ├─ 0.90-1.10? ────── YES ──→ ✓ GOOD! Acceptable
  │
  ├─ 0.80-1.20? ────── YES ──→ ⚠️ ACCEPTABLE (finite-size?)
  │
  └─ Outside 0.80-1.20? ────→ ❌ Something wrong!
      │                        Check:
      │                        - kBT=1? (use --kbt1)
      │                        - Long enough? (500 time units)
      │                        - Correct parameters?
      ▼
    DEBUG MODE
```

---

## Entropy Detection Algorithm

```
For each distribution (position and velocity):

1. Divide data into windows (1000 samples each)
   ┌────┬────┬────┬────┬────┬────┬────┬────┐
   │ W1 │ W2 │ W3 │ W4 │ W5 │ W6 │ W7 │ W8 │
   └────┴────┴────┴────┴────┴────┴────┴────┘

2. For each window, calculate entropy
   S = -Σ p(x) log p(x)

   Window 1: S₁ = 2.1
   Window 2: S₂ = 2.4
   Window 3: S₃ = 2.7  ← Rising (equilibrating)
   Window 4: S₄ = 2.9
   Window 5: S₅ = 3.0
   Window 6: S₆ = 3.01 ← Plateau start
   Window 7: S₇ = 3.02
   Window 8: S₈ = 3.01

3. Check stability for last 5 windows
   At window 8:
   - Last 5 entropies: [3.0, 3.01, 3.02, 3.01, 3.01]
   - Mean: 3.01
   - Std: 0.007
   - Relative std: 0.007/3.01 = 0.0023 = 0.23%
   - Is 0.23% < 5%? → YES! ✅ Stable!

4. Report steady-state point
   - Position stable at window 6 (6000 samples)
   - Velocity stable at window 8 (8000 samples)
   - Use LATER: window 8 (conservative)
   - Remove first 8000 samples as transient
```

---

## Timeline: What Happens When

```
Simulation Timeline (500 time units total):

t=0 ──────────┬──────────────────────────────────────────── t=500
              │
              ├─ t=0: Wall released (--auto-release)
              │
              ├─ t=0-20: Flat region removed (Step 1)
              │
              ├─ t=20-30: Position entropy rises
              │           (spatial equilibration)
              │
              ├─ t=30-50: Position entropy plateaus
              │
              ├─ t=20-80: Velocity entropy rises
              │           (thermal equilibration)
              │
              ├─ t=80: Velocity entropy plateaus
              │        → STEADY STATE DETECTED
              │        → Transient removal point
              │
              └─ t=80-500: Clean sinusoidal oscillations
                          Used for FFT analysis
                          ~6 full periods
                          Excellent frequency resolution

Number of samples:
  - Total: 60,000 samples (dt = 0.01)
  - Flat removed: 2,000 samples (3%)
  - Transient removed: 8,000 samples (13%)
  - Used for FFT: 50,000 samples (83%)
```

---

## Physical Interpretation

```
┌─────────────────────────────────────────────────────────┐
│                   PHYSICAL PICTURE                      │
└─────────────────────────────────────────────────────────┘

Initial State (t=0):
  ┌────────────────────────────────┐
  │ ●  ●    ● ●  │ ●  ●    ● ●     │  ← Random positions
  │   ●  ●  ●    │   ●  ●  ●       │     Random velocities
  │ ●    ● ●   ● │ ●    ● ●   ●    │     Wall at x=L₀
  └────────────────────────────────┘
         ↓ Collisions redistribute particles

Early Equilibration (t=20-80):
  ┌────────────────────────────────┐
  │ ● ●  ●   ● ● │●  ●   ● ●  ●    │  ← More uniform spatial
  │  ●  ● ● ●    │ ●  ● ● ●   ●    │     distribution forming
  │ ● ●  ●  ● ●  │● ●  ●  ● ●  ●   │     Velocities thermalizing
  └────────────────────────────────┘
         ↓ More collisions

Steady State (t>80):
  ┌────────────────────────────────┐
  │● ● ● ● ● ● ●│ ● ● ● ● ● ● ●  │  ← Uniform spatial dist.
  │ ● ● ● ● ● ● │● ● ● ● ● ● ●   │     Maxwell-Boltzmann v
  │● ● ● ● ● ● ●│ ● ● ● ● ● ● ●  │     Wall oscillates
  └────────────────────────────────┘
                ↕
          Wall oscillation
          ν = 0.0117 Hz
```

---

## Success Pattern Recognition

### ✅ GOOD RUN (Everything Works!)

```
Entropy plots:
  S[p(x)]  │      ╱─────────  ✓ Clear plateau
  S[p(v)]  │          ╱─────  ✓ Clear plateau (later)

Displacement:
  x(t)    │  ╱╲╱╲╱╲╱╲╱╲╱╲   ✓ Regular oscillations
          │╱

FFT:
  Power   │    ╱\           ✓ Sharp peak
          │___/  \___

Theory:
  Ratio = 0.9985            ✓ 95-105% range

Verdict: ✅ EXCELLENT
```

### ❌ BAD RUN (Too Short)

```
Entropy plots:
  S[p(x)]  │      ╱         ✗ No plateau (still rising!)
  S[p(v)]  │    ╱           ✗ No plateau

Displacement:
  x(t)    │  ╱──╱╲╱─       ✗ Irregular at start
          │╱

FFT:
  Power   │  ╱─\  ╱        ✗ Noisy (transients included)
          │_/   \_/

Theory:
  Ratio = 0.82              ⚠️ Outside 95-105%

Verdict: ⚠️ Run longer! (--steps=1000000)
```

---

## Quick Reference: Key Numbers

```
Parameter          Symbol    Typical Value
─────────────────────────────────────────
Particles (total)  N         100 (50/side)
Box half-length    L₀        20σ
Box height         A         10σ
Wall mass          M         200m
Temperature        kBT       1 (reduced)
Packing fraction   η         0.196
Speed of sound     cs        2.254
Transcendental K   K         0.653
Frequency          ν         0.0117
Period             T         85.3

Time scales:
  Collision time   τ_coll    ~0.5
  Spatial equil.   τ_space   ~10-20
  Thermal equil.   τ_therm   ~50-100
  Oscillation      T_osc     ~85

Data points:
  Total samples    N_total   60,000
  Flat removed     N_flat    2,000 (3%)
  Transient        N_trans   8,000 (13%)
  Used for FFT     N_FFT     50,000 (83%)
```

---

**See Also:**
- [QUICK_START.md](QUICK_START.md) - Commands to run
- [IMPLEMENTATION_COMPLETE.md](IMPLEMENTATION_COMPLETE.md) - Full documentation
- [README_ENTROPY_IMPLEMENTATION.md](README_ENTROPY_IMPLEMENTATION.md) - Overview

---

**Last Updated:** 2025-11-05
