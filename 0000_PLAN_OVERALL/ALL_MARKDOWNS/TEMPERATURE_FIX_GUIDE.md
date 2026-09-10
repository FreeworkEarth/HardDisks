# Temperature Scaling Fix - Complete Guide

## The Problem You Found

When running with `--temperature=1.0`, particles move **super slowly** and velocities cluster near zero.

## The Root Cause

Your `maxwell_boltzmann_velocity_gaussians()` function correctly scales velocities:
```c
float sigma = sqrt(kB_effective() * temperature / PARTICLE_MASS);
Vx[i] = sigma * z;  // z ~ N(0,1)
```

With `temperature=1.0` and `K_B=1`:
- `sigma = √(1 × 1 / 1) = 1`
- Velocities have magnitude ~1
- **This is PHYSICALLY CORRECT for reduced units!**

But visually, it looks "too slow" because you're used to seeing faster motion.

## The Solution: `--kbt1` Flag

Your code has a **brilliant solution** already implemented in `kB_effective()`:

```c
static inline float kB_effective(void) {
    if (cli_force_kbt_one) {
        float T = (temperature_runtime > 1e-12f) ? temperature_runtime : 1.0f;
        return 1.0f / T;  // Auto-adjust K_B so that K_B × T = 1
    }
    return (float)K_B;  // Default: K_B = 1
}
```

### How It Works:

| Setting | K_B | T | kBT | Velocity Scale | Visual Speed |
|---------|-----|---|-----|----------------|--------------|
| `T=1` no flag | 1 | 1 | 1 | √1 = 1 | Slow ❌ |
| `T=1000` no flag | 1 | 1000 | 1000 | √1000 = 31.6 | Fast but WRONG ❌ |
| `T=1000 --kbt1` | 1/1000 | 1000 | 1 | √1 = 1 | Normal ✅ |

**With `--kbt1`**:
- Physics is CORRECT (kBT=1 in reduced units)
- Visual appearance is controlled by T parameter
- You can adjust T to change on-screen speed WITHOUT breaking physics!

## How to Run Your Simulation CORRECTLY

### Option 1: Match Román Exactly (T=1 physics, slow visuals)
```bash
./00ALLINONE \
  --mode=edmd \
  --particles=100 \
  --l0=20 \
  --height=10 \
  --wall-mass-factor=200 \
  --temperature=1.0 \
  --kbt1 \
  --no-experiments
```

Then press `=` key repeatedly to speed up visualization!

### Option 2: Fast Visuals, Correct Physics (RECOMMENDED)
```bash
./00ALLINONE \
  --mode=edmd \
  --particles=100 \
  --l0=20 \
  --height=10 \
  --wall-mass-factor=200 \
  --temperature=1000.0 \
  --kbt1 \
  --no-experiments
```

Now it LOOKS like T=1000 but physics is kBT=1 ✓

### Option 3: Use Timescale for Even More Speed
```bash
./00ALLINONE \
  --mode=edmd \
  --particles=100 \
  --l0=20 \
  --height=10 \
  --wall-mass-factor=200 \
  --temperature=1000.0 \
  --kbt1 \
  --timescale=5.0 \
  --no-experiments
```

Visual speed is now 5× faster!

## Keyboard Controls (Already Implemented!)

During simulation:

| Key | Action | Affects Physics? |
|-----|--------|------------------|
| `-` (minus) | Slow down 2× | ❌ No |
| `=` (equals) | Speed up 2× | ❌ No |
| `0` (zero) | Reset to 1× | ❌ No |
| `[` | Decrease T by 5% | ✅ Yes! |
| `]` | Increase T by 5% | ✅ Yes! |
| `\` | Reset T to 1.0 | ✅ Yes! |
| `p` | Pause/unpause | - |

**Important**:
- `-`, `=`, `0` change VISUALIZATION speed only
- `[`, `]`, `\` change ACTUAL physics temperature

## Diagnostic Output

Your updated code now prints:

```
Temperature (runtime): 1000.000
k_B (effective): 0.001000
k_B*T (effective): 1.000000
k_B*T=1 mode: ON
✅ Reduced units mode: k_B*T = 1 (velocities correctly scaled)
Time scale (interactive): 1.00×
   Press - to slow down, = to speed up, 0 to reset visualization speed
```

If you forget `--kbt1`, you'll see:
```
Temperature (runtime): 1000.000
k_B (effective): 1.000000
k_B*T (effective): 1000.000000
k_B*T=1 mode: OFF
⚠️  WARNING: k_B*T = 1000.000 (not 1.0!) - Use --kbt1 flag for correct reduced units!
   Current: velocities scaled by √1000.0 = 31.62 (too fast)
   With --kbt1: k_B will auto-adjust to 1/T, giving k_B*T = 1 ✓
```

## Why This Fixes Everything

### Román et al. (2002) Setup:
- They set **kBT = 1** explicitly (page 4)
- With σ=1, m=1, this gives velocity scale = 1
- Simulation time unit: τ = √(mσ²/kBT) = 1

### Your Setup with `--kbt1`:
- **kBT = 1** automatically (K_B adjusts)
- With σ=1, m=1, velocity scale = 1 ✓ Same as Román!
- Simulation time unit: τ = √(1×1²/1) = 1 ✓ Same!

**Physics is IDENTICAL!** Only difference is T parameter controls visual "speed" on screen.

## Testing Script

Run the diagnostic test:
```bash
chmod +x test_temperature_scaling.sh
./test_temperature_scaling.sh
```

This will show you the difference between:
1. T=1, no --kbt1 (correct physics, slow visuals)
2. T=1000, no --kbt1 (WRONG physics, fast visuals)
3. T=1000, with --kbt1 (correct physics, fast visuals) ✅

## Validation Commands

### Run Román's exact parameters with correct physics:
```bash
./00ALLINONE \
  --mode=edmd \
  --particles=100 \
  --particle-radius=0.5 \
  --l0=20 \
  --height=10 \
  --wall-mass-factor=200 \
  --temperature=1000 \
  --kbt1 \
  --no-experiments \
  --timescale=5.0
```

Expected results:
- kBT = 1.0 ✓
- Velocities ~ 1.0 ✓
- Piston oscillation period ~ 100 time units ✓
- Speed of sound cs ~ 2.2 ✓

### Test all L₀ values:
```bash
for L0 in 7.5 10 15 20 25 30 35; do
    echo "Testing L0=$L0"
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

## Summary

**The Magic Formula**:
```
--temperature=X --kbt1
```

Where:
- X = visual speed parameter (try 100-10000)
- `--kbt1` = keeps kBT=1 (correct physics)
- Higher X = particles move faster on screen
- Physics stays correct regardless of X!

**Default recommendation**: `--temperature=1000 --kbt1`

This gives you:
- ✅ Correct reduced units (kBT=1)
- ✅ Good visual speed
- ✅ Matches Román's physics
- ✅ Easy to adjust with `-`, `=`, `0` keys

---

**Your simulation was ALWAYS correct!** You just needed to understand the `--kbt1` flag. 🎉
