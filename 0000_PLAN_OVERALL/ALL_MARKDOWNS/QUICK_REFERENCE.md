# Quick Reference: Hard Disk Piston Simulation

## TL;DR - What You Need to Know

### Running Simulations Correctly

**ALWAYS use the `--kbt1` flag for correct physics:**

```bash
./hspist3/00ALLINONE --mode=edmd --particles=100 --l0=20 --height=10 \
  --wall-mass-factor=200 --temperature=1000 --kbt1 --steps=300000
```

Without `--kbt1`, your physics will be wrong (kBT ≠ 1).

### What the Temperature Parameter Does

| Flag | kBT Value | Physics Correct? | Visual Speed |
|------|-----------|------------------|--------------|
| `--temperature=1 --kbt1` | 1.0 | ✅ YES | Slow (hard to see) |
| `--temperature=1000 --kbt1` | 1.0 | ✅ YES | Fast (easy to see) |
| `--temperature=1000` (no flag) | 1000.0 | ❌ NO | Too fast (broken) |

**Bottom line:** With `--kbt1`, the `temperature` parameter is just a "playback speed" for your eyes. The physics is always kBT=1.

## File Organization

```
HardDisks/
├── FINAL_SUMMARY_AND_ANSWERS.md     ← All 8 questions answered
├── UNITS_AND_TIMING_EXPLAINED.md    ← "Explain like 5" units guide
├── WALL_MOVEMENT_TEMPERATURE_ISSUE.md ← Why wall depends on kBT
├── TEMPERATURE_FIX_GUIDE.md          ← Original temperature confusion doc
├── test_wall_frequency_independence.sh ← Proves physics = f(kBT) only
├── QUICK_REFERENCE.md                ← THIS FILE
│
├── hspist3/
│   └── 00ALLINONE.c                  ← Your main simulation
│
└── EDMD_4VALIDATE/
    ├── README.md                     ← Validation suite docs
    ├── validate_roman_params.py      ← Theory calculator
    ├── analyze_edmd_validation.py    ← FFT analyzer
    └── validation_sim/
        ├── validate_roman.c          ← Minimal validation sim
        └── Makefile
```

## Common Tasks

### 1. Check if kBT is Correct

```bash
cd hspist3
./00ALLINONE --mode=edmd --particles=100 --l0=20 --temperature=1000 --kbt1 \
  --steps=10 --no-experiments 2>&1 | grep "k_B.*T"
```

**Should see:**
```
k_B*T (effective): 1.000000
k_B*T=1 mode: ON
✅ Reduced units mode: k_B*T = 1 (velocities correctly scaled)
```

**If you see warning:**
```
⚠️  WARNING: k_B*T = 1000.000 (not 1.0!)
```
→ You forgot the `--kbt1` flag!

### 2. Validate Against Román et al. (2002)

```bash
cd EDMD_4VALIDATE
python3 validate_roman_params.py              # Get theory predictions
cd validation_sim && make && ./validate_roman  # Run simulation
cd .. && python3 analyze_edmd_validation.py    # Analyze results
```

**Expected:** Frequency ratio (sim/theory) ≈ 0.95-1.05

### 3. Test Temperature Independence

```bash
./test_wall_frequency_independence.sh
```

This proves that T=1 and T=1000 give identical physics when using `--kbt1`.

### 4. Run Full L₀ Sweep

```bash
cd hspist3
for L0 in 7.5 10 15 20 25 30 35; do
  ./00ALLINONE --mode=edmd --particles=100 --l0=$L0 --height=10 \
    --wall-mass-factor=200 --temperature=1000 --kbt1 --steps=300000 \
    --no-experiments
done
```

## Understanding Your Results

### Speed of Sound

Your simulation gives cs in **reduced units** (dimensionless).

**To convert to physical units** (e.g., argon at 300K):

```python
# Define your system
sigma = 3.4e-10  # m (argon diameter)
m = 6.63e-26     # kg (argon mass)
kBT = 1.38e-23 * 300  # J (thermal energy at 300K)

# Conversion factor
tau = sigma * sqrt(m / kBT)  # time unit = 2.2e-12 s
velocity_unit = sigma / tau   # = 155 m/s

# Convert
cs_physical = cs_reduced * velocity_unit
# Example: cs_reduced = 2.54 → cs_physical = 394 m/s
```

### What Román Used

From the paper (page 4):
- σ = 1 (diameter)
- m = 1 (mass)
- **kBT = 1** (thermal energy)
- **T = 1** (reduced temperature, NOT 1 Kelvin!)

They did **NOT** use T=1000!

## The 8 Questions (Answered)

See [`FINAL_SUMMARY_AND_ANSWERS.md`](FINAL_SUMMARY_AND_ANSWERS.md) for full details:

1. **Is T=1000 means 1000K?** → NO (dimensionless with --kbt1)
2. **Can I use physical units?** → YES but impractical
3. **Does MB distribution work?** → YES (shape is universal)
4. **Particles 31.6× faster?** → Visual only, physics unchanged
5. **Just multiply by factor?** → YES (cs_physical = cs_reduced × σ/τ)
6. **What T did Román use?** → T=1 (reduced units)
7. **Wall depends on T?** → Depends on kBT, but --kbt1 keeps kBT=1
8. **Rename to reduced units?** → YES (much clearer)

## Keyboard Controls (in GUI mode)

- `-` : Slow down visualization (decrease time_scale_runtime)
- `=` : Speed up visualization (increase time_scale_runtime)
- `0` : Reset to default speed
- `ESC` : Quit

**Note:** These change VISUAL speed only, not physics!

## Warning Signs

### ⚠️ Physics is Broken If:

1. You see: `k_B*T = 1000.000` without `--kbt1` flag
2. Wall frequency changes when you change `temperature` parameter
3. Speed of sound depends on `temperature` parameter
4. Maxwell-Boltzmann distribution is too narrow or too wide

### ✅ Physics is Correct If:

1. You see: `k_B*T (effective): 1.000000` with `--kbt1`
2. Wall frequency is SAME for T=1 and T=1000 (when using --kbt1)
3. Validation frequency ratio ≈ 1.0 (within 5%)
4. Temperature in gas ≈ 1.0 throughout simulation

## Need More Detail?

- **Units confusion?** → Read [`UNITS_AND_TIMING_EXPLAINED.md`](UNITS_AND_TIMING_EXPLAINED.md)
- **Wall movement issue?** → Read [`WALL_MOVEMENT_TEMPERATURE_ISSUE.md`](WALL_MOVEMENT_TEMPERATURE_ISSUE.md)
- **All 8 questions?** → Read [`FINAL_SUMMARY_AND_ANSWERS.md`](FINAL_SUMMARY_AND_ANSWERS.md)
- **Validation failing?** → Read [`EDMD_4VALIDATE/README.md`](EDMD_4VALIDATE/README.md)

## The ONE Thing to Remember

**When using reduced units (σ=1, m=1, kBT=1):**

> Always use `--kbt1` flag. The `temperature` parameter is just for visual speed, not physics.

Without this, everything breaks!

---

Last updated: November 2024
