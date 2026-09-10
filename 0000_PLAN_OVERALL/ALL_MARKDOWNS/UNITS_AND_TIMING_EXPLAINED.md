# Complete Guide: Units, Timing, and Speed of Sound Validation

## Table of Contents
1. [Understanding Reduced Units](#reduced-units)
2. [Román et al. (2002) Unit System](#roman-units)
3. [Simulation vs Rendering Speed](#timing)
4. [Adding Timestep Controls](#controls)
5. [Why Your Results Differ from Román](#discrepancy)

---

## 1. Understanding Reduced Units (Explained Like You're 5) {#reduced-units}

### The Big Confusion: What Does T=1000 Mean?

Imagine you're thinking:
> "If I set T=1000, does the gas have temperature 1000 Kelvin?"

**Answer: NO! Not unless you're using physical units!**

Let me show you with concrete examples:

---

## Two Ways to Do the SAME Simulation

### Method 1: Physical Units (Real Argon Gas)

```c
// Real argon gas at room temperature
#define SIGMA 3.4e-10        // 3.4 Angstroms (particle diameter)
#define MASS 6.63e-26        // kg (argon atom)
#define K_B 1.38e-23         // J/K (Boltzmann constant)
#define TEMPERATURE 300.0    // Kelvin (ACTUAL TEMPERATURE!)

// Initialize velocity:
float sigma_v = sqrt(K_B * TEMPERATURE / MASS);
Vx[i] = sigma_v * gaussian_random();
```

**Result:**
- T = 300 means **300 Kelvin** (room temperature)
- Velocities: v_rms = √(kBT/m) = √(1.38×10⁻²³ × 300 / 6.63×10⁻²⁶) = 432 m/s
- This is **real physics** with **real units**

---

### Method 2: Reduced Units (What Román & MD Use)

```c
// Reduced units: choose σ=1, m=1, kBT=1
#define SIGMA 1.0            // Unit length (dimensionless)
#define MASS 1.0             // Unit mass (dimensionless)
#define K_B 1.0              // Dimensionless constant
#define TEMPERATURE 1.0      // Dimensionless number (NOT Kelvin!)

// Initialize velocity:
float sigma_v = sqrt(K_B * TEMPERATURE / MASS);
Vx[i] = sigma_v * gaussian_random();
```

**Result:**
- T = 1.0 is **NOT 1 Kelvin!**
- T = 1.0 is a **dimensionless number** (no units at all!)
- Velocities: v_rms = √(1 × 1 / 1) = 1 (dimensionless)

---

## Are These Two Simulations Equivalent?

**YES!** They describe the EXACT SAME PHYSICS!

**Proof by conversion:**

```
Time unit: τ = σ√(m/ε)   where ε = kBT

Physical units (argon):
τ = (3.4×10⁻¹⁰) × √(6.63×10⁻²⁶ / (1.38×10⁻²³ × 300))
  = 2.2×10⁻¹² seconds = 2.2 picoseconds

Reduced units:
τ = 1 × √(1 / 1) = 1 (dimensionless time unit)
```

**Translation:** 1 reduced time unit = 2.2 picoseconds of real time!

**Speed of sound conversion:**
```
Length unit: σ = 3.4×10⁻¹⁰ m
Time unit: τ = 2.2×10⁻¹² s
Velocity unit: σ/τ = 3.4×10⁻¹⁰ / 2.2×10⁻¹² = 155 m/s

Reduced units: cs = 2.2 (dimensionless)
Physical units: cs = 2.2 × 155 m/s = 341 m/s ✓
```

---

## Why Not Use Physical Units?

**You CAN!** But it causes problems:

### Problem 1: Numerical Precision Loss
```
Position: x = 3.4×10⁻¹⁰ meters
Velocity: v = 432 m/s
Force: F = 1.22×10⁻¹¹ Newtons

Computer precision: ~10⁻¹⁶
Your numbers: range from 10⁻²³ to 10²

YOU LOSE PRECISION! ❌
```

### Problem 2: Equations Get Messy
```
With physical units:
F = (1.38×10⁻²³ × 300 / 3.4×10⁻¹⁰) × (distance)

With reduced units:
F = (1.0 × 1.0 / 1.0) × (distance) = distance

Much cleaner! ✓
```

### Problem 3: Every System Different
```
Argon: σ = 3.4 Å, m = 6.63×10⁻²⁶ kg
Water: σ = 2.8 Å, m = 2.99×10⁻²⁶ kg

With reduced units: ALWAYS σ=1, m=1!
Code works for ANY system! ✓
```

---

### What Does "In Units of kBT" Mean?

When physicists say something is "in units of kBT", they mean that quantity is **divided by kBT** to make it dimensionless.

**Example: Energy**
- Real energy: E = 3.0 × 10⁻²¹ Joules
- In units of kBT (at T=300K): E/(kBT) = 3.0×10⁻²¹ / (1.38×10⁻²³ × 300) = 0.725
- We write: E = 0.725 kBT (or just E = 0.725 in reduced units)

**Key Insight:** The distribution shape is IDENTICAL in both unit systems!

```
Maxwell-Boltzmann: P(v) ∝ exp(-m·v² / 2kBT)

This is dimensionless because:
(m·v²) / (kBT) = energy / energy = dimensionless!
```

### The Standard Reduced Unit System

For molecular dynamics, we choose fundamental units and express everything relative to them:

| Quantity | Real Units | Reduced Units | Conversion |
|----------|------------|---------------|------------|
| **Length** | meters | σ (particle diameter) | x' = x/σ |
| **Mass** | kg | m (particle mass) | M' = M/m |
| **Energy** | Joules | kBT | E' = E/(kBT) |
| **Time** | seconds | τ = √(mσ²/kBT) | t' = t/τ |
| **Velocity** | m/s | √(kBT/m) | v' = v/√(kBT/m) |
| **Temperature** | Kelvin | kBT | T' = T (when kB=1) |

**Key Point**: Once you choose σ=1, m=1, kB=1, T=1, ALL your equations simplify!

---

## 2. Román et al. (2002) Unit System {#roman-units}

### Their Explicit Statement (Page 4, Section III):

> "For simplicity, we have chosen **σ=1** and **m=1**. This choice is equivalent to **scaling distances by σ and masses by m**. **Energies are scaled by kBT**, and **time is scaled by (mσ²/kBT)^(1/2)**."

### Breaking It Down:

#### What They Set to 1:
```
σ = 1  (particle diameter)
m = 1  (particle mass)
kBT = 1  (thermal energy)
```

#### What This Means for Their Simulation:

**Distances:**
- Box half-length: L₀ = 20σ = 20 (in code)
- Cross-section: A = 10σ = 10 (in code)

**Masses:**
- Particle mass: m = 1
- Piston mass: M = 200m = 200 (in code)

**Velocities:**
- Maxwell-Boltzmann: v ~ √(kBT/m) = √(1/1) = 1
- So velocities have magnitudes around 1 in their simulation

**Time:**
- Time unit: τ = √(mσ²/kBT) = √(1×1²/1) = 1
- Simulation runs for t = 3000τ = 3000 time units

**Temperature:**
- They set kBT = 1
- With kB = 1, this means T = 1
- The temperature IS correct with T=1!

### The Speed of Sound Formula (Equation 26):

For a 2D hard disk gas:
```
cs = √(2kBT/m) × √[(1 + η + 3aη² - aη³)/(1-η)³]
```

With **kBT=1, m=1**:
```
cs = √2 × √[(1 + η + 3aη² - aη³)/(1-η)³]
   = 1.414 × √[(1 + η + ...)/(1-η)³]
```

For L₀=20, η=0.196:
```
cs(Henderson) = 1.414 × √[(1 + 0.196 + 0.00145 - 0.00019)/(1-0.196)³]
              = 1.414 × √[1.197/0.519]
              = 1.414 × 1.527
              = 2.16  ✓ (matches their Table I!)
```

**Your Results vs Román:**
```
L₀    η      Your cs    Román cs    Theory cs   Match?
7.5   0.524  7.04       5.99±0.09   5.45        ❌ Too high
10    0.393  4.13       3.78±0.08   3.59        ❌ Too high
15    0.262  2.54       2.61±0.03   2.53        ✅ EXCELLENT!
20    0.196  2.59       2.20±0.02   2.16        ❌ Too high
25    0.157  2.08       2.10±0.02   1.97        ❌ Slightly high
30    0.131  1.24       1.89±0.02   1.86        ❌ Too LOW!
35    0.112  1.45       1.81±0.02   1.79        ❌ Too LOW!
```

---

## 3. Simulation Speed vs Rendering Speed {#timing}

This is the **KEY** to understanding why particles look "slow" or "fast"!

### Three Different "Speeds":

#### A. **Physical Simulation Speed** (What Your Code Calculates)
- How fast simulation time advances
- Controlled by: `FIXED_DT` timestep
- Example: dt = 0.001 means each integration step advances simulation by 0.001 time units
- This is **PHYSICS** and should match Román's τ = √(mσ²/kBT)

#### B. **Rendering Frame Rate** (How Fast You SEE Updates)
- How many times per second SDL draws to screen
- Typically 60 FPS = screen updates every ~16ms
- Controlled by: SDL event loop timing

#### C. **Time Scale Factor** (Interactive Speed Control)
- Your `--timescale` parameter
- Multiplies how many simulation steps happen per frame
- Example: `--timescale=2.0` means simulation runs 2× faster on screen
- **Does NOT change physics!** Just changes visualization speed

### The Problem: Confusion Between These Three

When you say "particles are slow with T=1", you might mean:
1. Velocities in data are small → **Physics issue** (check velocity magnitudes)
2. Particles move slowly ON SCREEN → **Rendering issue** (adjust timescale)

### Your Current Settings (from 00ALLINONE.c):

```c
#define FIXED_DT 0.001f          // Simulation timestep
#define SUBSTEPS 8               // Integration substeps per frame
static float time_scale_runtime = 1.0f;  // Interactive speed multiplier
```

Each SDL frame:
```
Simulation time advanced = FIXED_DT × SUBSTEPS × time_scale_runtime
                         = 0.001 × 8 × 1.0
                         = 0.008 time units per frame
```

At 60 FPS:
```
Simulation speed = 0.008 × 60 = 0.48 time units per second of real time
```

To simulate 3000 time units:
```
Real time needed = 3000 / 0.48 = 6250 seconds = 104 minutes!
```

### Solution: Use Time Scale!

With `--timescale=10.0`:
```
Simulation speed = 0.001 × 8 × 10 × 60 = 4.8 time units/second
Time for 3000 units = 3000 / 4.8 = 625 seconds = 10 minutes
```

**Particles look 10× faster on screen, but physics is UNCHANGED!**

---

## 4. Adding Timestep Controls (O/P Keys) {#controls}

Let me add keyboard controls to adjust simulation speed on-the-fly:

### Implementation Plan:

```c
// In your SDL event handler, add:

case SDLK_o:  // 'O' key - Speed UP
    time_scale_runtime *= 1.5f;
    if (time_scale_runtime > 100.0f) time_scale_runtime = 100.0f;
    printf("⏩ Speed increased: timescale = %.2f×\n", time_scale_runtime);
    break;

case SDLK_p:  // 'P' key - Slow DOWN
    time_scale_runtime /= 1.5f;
    if (time_scale_runtime < 0.01f) time_scale_runtime = 0.01f;
    printf("⏪ Speed decreased: timescale = %.2f×\n", time_scale_runtime);
    break;

case SDLK_0:  // '0' key - RESET to 1×
    time_scale_runtime = 1.0f;
    printf("⏸  Speed reset: timescale = 1.0×\n");
    break;
```

This lets you:
- Press **O** repeatedly to speed up visualization
- Press **P** repeatedly to slow down
- Press **0** to reset to normal speed
- **Physics stays correct!** Only affects how fast you watch it

---

## 5. Why Your Results Differ from Román {#discrepancy}

### The Pattern in Your Data:

Looking at your plot, I see a **systematic issue**:

**High Density (small L₀, large η):**
- L₀=7.5, η=0.524: Your cs=7.04 vs Theory=5.45 (**29% too high**)
- L₀=10, η=0.393: Your cs=4.13 vs Theory=3.59 (**15% too high**)

**Medium Density:**
- L₀=15, η=0.262: Your cs=2.54 vs Theory=2.53 (**PERFECT!** ✅)

**Low Density (large L₀, small η):**
- L₀=30, η=0.131: Your cs=1.24 vs Theory=1.86 (**33% too LOW!**)
- L₀=35, η=0.112: Your cs=1.45 vs Theory=1.79 (**19% too LOW!**)

### Possible Causes:

#### 1. **Finite Size Effects** (Most Likely)
- With N=100 particles, your system is VERY SMALL
- Román notes (Table II) finite size effects are present
- At high density, particles get "caged" → sound travels faster
- At low density, box is "too empty" → sound wavelength doesn't fit well

**Fix**: Run with more particles (N=256, 512, 1024) to approach thermodynamic limit

#### 2. **Temperature Initialization Issue**
- If velocities aren't exactly Maxwell-Boltzmann with T=1, cs will be wrong
- cs ∝ √T, so if T_actual = 1000, then cs_measured = √1000 × cs_theory ≈ 31.6× too high

**Check**: After initialization, measure average kinetic energy:
```c
double ke_avg = 0;
for (int i=0; i<N; i++) {
    ke_avg += 0.5 * m * (Vx[i]*Vx[i] + Vy[i]*Vy[i]);
}
ke_avg /= N;
double T_measured = ke_avg;  // In 2D: <KE> = kBT, with kB=1
printf("Measured temperature: T = %.6f\n", T_measured);
```

If T_measured ≠ 1.0, your velocity scaling is wrong!

#### 3. **Effective Length Correction**
- Román uses L_eff = L₀ - σ in Eq. (18) due to finite particle size
- Your plot uses L₀-1, which is correct
- But maybe the correction should be different for the piston?

#### 4. **Piston Thickness**
- Román assumes "piston of mass M and zero width"
- Your piston has `WALL_THICKNESS` > 0
- This changes the effective volume and thus η

**Check**: What is your wall thickness?
```bash
grep WALL_THICKNESS hspist3/00ALLINONE.c
```

### Recommendation for Validation:

Run a **controlled test** matching Román Fig. 3 EXACTLY:
```bash
./00ALLINONE \
  --mode=edmd \
  --particles=100 \
  --l0=20 \
  --height=10 \
  --wall-mass-factor=200 \
  --particle-radius=0.5 \
  --temperature=1.0 \
  --no-experiments \
  --timescale=10.0
```

Then measure:
1. Average particle velocity magnitude
2. Temperature from kinetic energy
3. Wall oscillation frequency
4. Compare these exact numbers with theory

---

## Summary

1. **Román's units are correct**: σ=1, m=1, kBT=1, T=1
2. **Your code should work with T=1** if velocities are properly initialized
3. **"Slow particles"** is likely a **visual rendering speed issue**, not physics
4. **Add O/P keys** to control visual playback speed
5. **Your L₀=15 match is PERFECT!** This proves your code CAN work correctly
6. **Discrepancies at other L₀** are likely due to **finite size effects** and need more particles

The good news: **Your simulation IS working correctly for L₀=15!** The other cases need investigation of finite-size effects and possibly temperature measurement.

Would you like me to:
1. Add the O/P keyboard controls to your code?
2. Add temperature measurement diagnostics?
3. Create a script to run systematic validation tests?
