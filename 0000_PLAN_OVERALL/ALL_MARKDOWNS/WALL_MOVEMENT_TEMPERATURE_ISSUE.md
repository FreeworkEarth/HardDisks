# CRITICAL ISSUE: Wall Movement Depends on Temperature!

## Your Observation (VERY IMPORTANT!)

You noticed:
> "Wall movement depends on temperature. From wall position FFT, I get speed of sound. But if wall movement depends on temperature, I get different speeds of sound for different T!"

**YOU ARE ABSOLUTELY RIGHT!** This is the KEY issue!

---

## The Problem Explained

### The Román Formula (Equation 4):
```
ν = (cs / 2πL₀) × √(2Nm/M)
```

This frequency depends on:
- cs (speed of sound) ∝ √T
- L₀ (box length)
- N (particles per side)
- M (piston mass)

**Key insight:** cs ∝ √T, so if you change T, you change cs!

---

## What Román Actually Did

Looking at page 4, Section III of the paper:

> "Energies are scaled by kBT, and time is scaled by (mσ²/kBT)^(1/2)"

**They explicitly state:** kBT = 1

### Their Exact Setup:
```
σ = 1  (particle diameter)
m = 1  (particle mass)
kBT = 1  (thermal energy)
```

**This means:**
```
With kB = 1:
T = 1 (reduced units, NOT 1 Kelvin!)
```

### Their Speed of Sound Formula (Equation 26):
```
cs = √(2kBT/m) × √[(1 + η + 3aη² - aη³)/(1-η)³]
```

With kBT = 1, m = 1:
```
cs = √2 × √[(1 + η + ...)/(1-η)³]
```

**This is what you should get!**

---

## Why `--kbt1` with `temperature=1000` DOESN'T Change Physics

When you use `--kbt1 --temperature=1000`:

```c
K_B = 1.0 / temperature = 1/1000 = 0.001
kBT = 0.001 × 1000 = 1.0
```

**In the speed of sound formula:**
```c
cs = √(2kBT/m) = √(2 × 1.0 / 1.0) = √2 = 1.414
```

**The speed of sound is STILL calculated with kBT=1!**

The `temperature=1000` parameter ONLY affects:
1. Visualization speed (rendering)
2. Nothing else!

**So your wall movement frequency should be IDENTICAL** whether you use:
- `temperature=1.0 --kbt1`
- `temperature=1000 --kbt1`

Both give kBT=1, so both give the same cs!

---

## TEST THIS NOW!

Run these two simulations and compare:

### Test 1: T=1 with --kbt1
```bash
./00ALLINONE \
  --mode=edmd \
  --particles=100 \
  --l0=20 \
  --height=10 \
  --wall-mass-factor=200 \
  --temperature=1.0 \
  --kbt1 \
  --steps=300000 \
  --no-experiments
```

### Test 2: T=1000 with --kbt1
```bash
./00ALLINONE \
  --mode=edmd \
  --particles=100 \
  --l0=20 \
  --height=10 \
  --wall-mass-factor=200 \
  --temperature=1000.0 \
  --kbt1 \
  --steps=300000 \
  --no-experiments
```

**Prediction:** Both should give:
- Same wall oscillation frequency ν
- Same speed of sound cs
- Only difference: Test 2 LOOKS faster on screen

**If they give DIFFERENT frequencies, there's a bug in the code!**

---

## Checking Your Current Results

Looking at your plot `speed_of_sound_per_L0_with_eta.pdf`:

```
L₀    η      Your cs    Román cs    Ratio
7.5   0.524  7.04       5.99        1.18
10    0.393  4.13       3.78        1.09
15    0.262  2.54       2.61        0.97  ← PERFECT!
20    0.196  2.59       2.20        1.18
25    0.157  2.08       2.10        0.99
30    0.131  1.24       1.89        0.66  ← TOO LOW!
35    0.112  1.45       1.81        0.80
```

**Question:** What temperature did you use for each of these?

### Hypothesis 1: You Used Different Temperatures
If you ran:
- Small L₀ with T=1 (slow) → Maybe you stopped early?
- Large L₀ with T=1000 (fast) → Got full data?

This could cause artificial differences!

### Hypothesis 2: Finite Size Effects
Román mentions (Table II) finite size effects are real with N=100.

Your L₀=15 match is perfect because:
- η = 0.262 (medium density)
- Box not too small, not too large
- N=100 is "just right" for this size

But at:
- Small L₀ (high η): Particles cage each other → cs too high
- Large L₀ (low η): Too few collisions per box → cs too low

---

## What Temperature Did Román Use?

From the paper (page 4):

> "For simplicity, we have chosen **σ=1** and **m=1**. [...] **Energies are scaled by kBT**, and time is scaled by (mσ²/kBT)^(1/2)."

Then in Figure 3 caption:

> "Plot of the displacement of the piston ξ(t)=L(t)−L₀ for L₀=20, N=N₀/2=50, M=200, A=10, and **kBT=1**."

**Answer: They used kBT = 1 in reduced units!**

With kB=1, this means **T = 1** (reduced units, dimensionless).

**NOT T=1000!**

---

## Converting to Physical Units

### If You Want Real Argon Gas at 300K:

**Step 1:** Define your unit system
```
σ = 3.4×10⁻¹⁰ m  (argon diameter)
m = 6.63×10⁻²⁶ kg  (argon mass)
ε = kBT = 1.38×10⁻²³ × 300 = 4.14×10⁻²¹ J
```

**Step 2:** Calculate conversion factors
```
Length unit: σ = 3.4×10⁻¹⁰ m
Time unit: τ = σ√(m/ε) = 2.2×10⁻¹² s
Velocity unit: σ/τ = 155 m/s
```

**Step 3:** Convert results
```
Your reduced units: cs = 2.54 (at L₀=15)
Physical units: cs = 2.54 × 155 m/s = 394 m/s
```

**But this is for your SIMULATED hard disk gas, not real argon!**

Real argon speed of sound ≈ 323 m/s (close!)

---

## Action Items for You

### 1. Verify kBT=1 in All Simulations
Run this diagnostic for each L₀:

```bash
./00ALLINONE \
  --mode=edmd \
  --particles=100 \
  --l0=15 \
  --height=10 \
  --wall-mass-factor=200 \
  --temperature=1000 \
  --kbt1 \
  --no-experiments \
  2>&1 | grep "k_B.*T.*effective"
```

Should show: `k_B*T (effective): 1.000000`

### 2. Re-run All L₀ with Same Settings
Use IDENTICAL parameters for all:
```bash
for L0 in 7.5 10 15 20 25 30 35; do
  ./00ALLINONE \
    --mode=edmd \
    --particles=100 \
    --l0=$L0 \
    --height=10 \
    --wall-mass-factor=200 \
    --temperature=1000 \
    --kbt1 \
    --steps=300000 \
    --no-experiments
done
```

### 3. Check Wall Position Data Quality
For each L₀, verify:
- Simulation ran full 3000 time units
- Wall oscillates clearly
- FFT has clear peak

---

## Variable Naming Convention

**YES!** Let's rename variables to be clear about units:

### Current (Confusing):
```c
#define TEMPERATURE 1000.0f
float temperature_runtime;
```

### Better (Clear):
```c
#define TEMPERATURE_REDUCED 1000.0f  // Dimensionless, for visual speed only
float T_reduced;  // Actual reduced temperature
float kBT_reduced;  // Should always = 1.0 in reduced units mode
```

### In Comments:
```c
// Initialize Maxwell-Boltzmann velocities
// For reduced units: kBT = 1, so v_rms = √(2/m) = √2
float sigma_v = sqrt(kB_effective() * T_reduced / PARTICLE_MASS);
```

---

## Summary

**Román used:**
- kBT = 1 (reduced units, dimensionless)
- T = 1 (with kB = 1)
- NOT T = 1000!

**Your `--kbt1` flag:**
- Correctly enforces kBT = 1
- The `temperature` parameter becomes just a "visual speed dial"
- Physics is independent of the temperature parameter value (as long as --kbt1 is used)

**Wall movement frequency SHOULD be the same** whether you use T=1 or T=1000, as long as you use `--kbt1`!

**If it's not, there's a bug!**

Let's test this now!
