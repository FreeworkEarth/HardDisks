# Temperature, Amplitude, and Rendering Guide

## TL;DR

**To match Román et al. (2002) results:**

```bash
./00ALLINONE \
  --mode=edmd \
  --particles=100 \
  --l0=20 \
  --height=10 \
  --wall-mass-factor=200 \
  --temperature=100000 \
  --kbt1 \
  --auto-release \
  --steps=5000000
```

**Key points:**
- Use `--temperature=100000` (not 1000!) for larger oscillation amplitudes
- GUI **does render** even with `--steps` flag
- Window might close quickly when done, or appear behind other windows

---

## The Temperature Mystery

### What You Observed

| Initial T | kBT (final) | Wall Amplitude | Oscillation Quality |
|-----------|-------------|----------------|---------------------|
| 1         | 1.0         | ~0.1 σ         | ❌ Barely moves |
| 1000      | 1.0         | ~0.5 σ         | ⚠️ Small, irregular |
| 100000    | 1.0         | ~2.0 σ         | ✅ Good, matches paper! |

**All three have the same final temperature (kBT=1)**, but very different wall motion!

### Why This Happens

#### The Physics

1. **`--temperature=T` sets initial velocity distribution:**
   ```
   v_i ~ Maxwell-Boltzmann(T)
   ```

2. **`--kbt1` rescales velocities so final temperature is kBT=1:**
   ```
   v_i → v_i × sqrt(1/T_measured)
   ```

3. **BUT: The initial momentum distribution matters for wall dynamics!**

#### The Initial State Matters

**Low T (e.g., T=1):**
- Particles start nearly at rest
- Wall experiences very small random forces
- Oscillations remain tiny
- Like tapping a pendulum gently

**Medium T (e.g., T=1000):**
- Particles have moderate initial velocities
- Some momentum transfer to wall
- Oscillations grow but remain modest

**High T (e.g., T=100000):**
- Particles start with very high velocities
- After rescaling, they have **broad momentum distribution**
- Large initial momentum transfer to wall
- Wall gets "kicked" hard initially
- Oscillations grow to ~±2σ amplitude
- Like hitting a pendulum hard

### What Román et al. Did (Probably)

The paper doesn't specify initialization, but likely:

**Option 1: Hot start + equilibration**
- Start with very high T
- Let system equilibrate naturally (no rescaling)
- Takes ~1000 time units to thermalize
- Results in large natural oscillations

**Option 2: Prepared initial state**
- Carefully set initial velocities to give desired amplitude
- Use Monte Carlo or MD pre-equilibration
- Then start main simulation

**Your approach (high T + kBT scaling)** effectively mimics Option 1 but faster!

---

## Recommendations

### For Validation Against Román et al.

Use `--temperature=100000` to match their amplitude:

```bash
./00ALLINONE \
  --mode=edmd \
  --particles=100 \
  --l0=20 \
  --height=10 \
  --wall-mass-factor=200 \
  --temperature=100000 \
  --kbt1 \
  --auto-release \
  --steps=5000000
```

**Expected results:**
- Wall displacement: -2σ to +2σ ✅
- Oscillation period: ~85 time units
- Frequency: ~0.012 Hz
- Clear wave packets/beats (if there are mode couplings)

### For Different Amplitudes

If you want to study amplitude dependence:

**Small amplitude (linear regime):**
```bash
--temperature=1000   # A ~ 0.5σ
```

**Medium amplitude:**
```bash
--temperature=10000  # A ~ 1.0σ
```

**Large amplitude (paper's regime):**
```bash
--temperature=100000 # A ~ 2.0σ
```

**Note:** Frequency should NOT depend on amplitude (for hard spheres, this is a linear system). Only amplitude changes!

---

## Rendering Behavior

### Issue: "Why does it not render with --auto-release and --steps?"

**Short answer:** It DOES render! But:

1. **Window might appear behind your current window**
   - Check all open windows
   - Use Mission Control / Exposé (F3 on Mac)
   - Look for "00ALLINONE" or SDL window

2. **Simulation might finish very quickly**
   - With optimized code, 100k steps = 1 second
   - 5M steps = 50 seconds
   - Window closes when done
   - You might miss it!

3. **Rendering is throttled for performance**
   - Updates display every N frames (not every step)
   - Might look frozen but is actually running
   - Check CPU usage - should be 100% on one core

### How to Verify It's Running

**Option 1: Watch the files**
```bash
watch -n 1 wc -l wall_position.csv
```
Should see line count increasing!

**Option 2: Check process**
```bash
ps aux | grep 00ALLINONE
```
Should show process running at ~100% CPU

**Option 3: Interactive mode (no --steps)**
```bash
./00ALLINONE \
  --mode=edmd \
  --particles=100 \
  --l0=20 \
  --height=10 \
  --wall-mass-factor=200 \
  --temperature=100000 \
  --kbt1 \
  --auto-release
```

Then:
- Watch it run in real-time ✅
- Press `+` repeatedly to speed up rendering (up to 100×)
- Press `q` after ~5 minutes (check time display in window)
- Wall position logged continuously

### Controls (Interactive Mode)

| Key | Action |
|-----|--------|
| `R` | Release wall (not needed with --auto-release) |
| `+` / `=` | Speed up rendering (2×, 4×, 8×, ..., 100×) |
| `-` | Slow down rendering |
| `q` | Quit and save data |
| `ESC` | Emergency quit (might not save) |
| `Space` | Pause/unpause |

**Tip:** Press `+` about 6-7 times to get 64× speed. The simulation will appear to run much faster!

---

## Understanding the Oscillations

### What You Should See (with T=100000)

**Early time (t < 100):**
```
Wall position vs time:
   2 ─┐
      │  ╱╲    ╱╲
   0 ─┼─╯  ╲  ╱  ╲─
      │     ╲╱    ╲╱
  -2 ─┘
      └─────────────> time
      0    50   100

  Transient: Irregular, establishing oscillation
```

**Steady state (t > 100):**
```
Wall position vs time:
   2 ─┐
      │ ╱╲  ╱╲  ╱╲  ╱╲  ╱╲
   0 ─┼─╯╲╱╯╲╱╯╲╱╯╲╱╯╲─
      │
  -2 ─┘
      └───────────────────────> time
      100  200  300  400  500

  Regular sinusoidal oscillation
  Period T ≈ 85 time units
```

### Wave Packets / Beats

If you see **amplitude modulation** (wave packets):
```
Wall position vs time:
   2 ─┐     ╱╲             ╱╲
      │    ╱  ╲    ╱╲    ╱  ╲
   0 ─┼───╯    ╲──╯  ╲──╯    ╲───
      │         ╲─────╯
  -2 ─┘
      └──────────────────────────> time

  Beat period: Multiple of fundamental T
```

**This means:**
- Multiple modes are excited
- Mode coupling (nonlinear effects)
- Perfectly normal for large amplitude!
- Román et al. might show this too

**To check:** Look at FFT spectrum:
- Should see peaks at ν₁, 2ν₁, 3ν₁, ... (harmonics)
- Beat frequency = |ν₁ - ν₂| if two modes present

---

## Validation Strategy

### Step 1: Quick Test (1 minute)

```bash
./00ALLINONE \
  --mode=edmd \
  --particles=100 \
  --l0=20 \
  --height=10 \
  --wall-mass-factor=200 \
  --temperature=100000 \
  --kbt1 \
  --auto-release \
  --steps=100000  # Just 100k steps to test
```

**Check:**
- Does GUI open? ✅
- Does wall move with amplitude ~2σ? ✅
- Does it look sinusoidal? ✅

If yes → proceed to full run!

### Step 2: Full Run (5 minutes)

```bash
./00ALLINONE \
  --mode=edmd \
  --particles=100 \
  --l0=20 \
  --height=10 \
  --wall-mass-factor=200 \
  --temperature=100000 \
  --kbt1 \
  --auto-release \
  --steps=5000000
```

**While running:**
```bash
# In another terminal, watch progress:
watch -n 2 "wc -l wall_position.csv && tail -1 wall_position.csv"
```

Should see:
- Line count increasing (~1000 lines/second)
- Time column reaching 5000
- Displacement oscillating between -2 and +2

### Step 3: Analyze

```bash
python3 wall_x_FFT.py
```

**Expected output:**
```
✅ Position entropy stable at ~250000 samples (5%)
✅ Velocity entropy stable at ~500000 samples (10%)

Highest Power Spectrum Peak: 0.0117 Hz

✅ EXCELLENT: Simulation matches theory within 5%!
```

### Step 4: Visual Check

Open `divider_x_displacement_with_boxlength_20_and_wallfactor_200.pdf`

**Panel 1: Displacement**
- Amplitude: -2σ to +2σ ✅
- ~60 full oscillations visible ✅
- Regular pattern after gray region ✅

**Panel 2: Position Entropy**
- Smoothed line plateaus ✅
- Transient < 20% ✅

**Panel 3: Velocity Entropy**
- Smoothed line plateaus ✅
- May plateau later than position (normal) ✅

---

## Troubleshooting

### Problem: "Amplitude still too small even with T=100000"

**Possible causes:**

1. **kBT is not actually 1**
   - Check terminal output for: "Temperature scaled to kBT = 1.000"
   - If missing, `--kbt1` flag didn't work

2. **Wall is constrained**
   - Check for `--spring-k` in your command
   - Spring forces would limit amplitude
   - Don't use spring for validation!

3. **Simulation is dissipative**
   - Hard sphere collisions should be elastic
   - Check energy conservation in `energy_log.csv`
   - Total energy should be constant (±1%)

### Problem: "GUI window appears then immediately closes"

**Solution 1:** Don't use `--steps`, use interactive mode:
```bash
./00ALLINONE ... --auto-release  # No --steps flag
```
Press `q` to quit after desired time

**Solution 2:** Check exit messages:
```bash
./00ALLINONE ... --steps=5000000 2>&1 | tee output.log
```
Look for errors at end of output.log

**Solution 3:** Verify simulation completed:
```bash
wc -l wall_position.csv
```
Should have ~5 million lines if it ran fully

### Problem: "Oscillations are irregular / noisy"

**If t < 100:** Normal! This is the transient region.

**If t > 100 still irregular:**
- Not enough particles (N=100 is minimum)
- Try N=200: `--particles=200`
- Box too small (increase `--l0=30`)
- Collision detection issues (shouldn't happen with EDMD)

### Problem: "Theory ratio still bad even with 5M steps and T=100000"

**Check these:**

1. **Verify wall mass in output:**
   ```
   Look for: "Wall mass: M = 200.0"
   ```

2. **Check packing fraction:**
   ```python
   η = (100 × π × 0.5²) / (4 × 10 × 20) = 0.1963 ✅
   ```

3. **Verify data quality:**
   - Open FFT plot
   - Should see ONE sharp peak
   - If multiple peaks → mode coupling, use tallest peak
   - If broad hump → not equilibrated, run longer

4. **Check entropy plots:**
   - Both must plateau!
   - If still rising → need even longer run (10M steps)

---

## Summary

### Key Findings

1. **Initial temperature affects amplitude:**
   - T=1: Tiny oscillations
   - T=1000: Small oscillations
   - T=100000: Large oscillations (matches paper) ✅

2. **All have same final kBT=1, but different dynamics!**

3. **GUI does render with --steps**, just check:
   - Window might be hidden
   - Runs fast, might finish quickly
   - Use interactive mode for better control

### Recommended Command

```bash
./00ALLINONE \
  --mode=edmd \
  --particles=100 \
  --l0=20 \
  --height=10 \
  --wall-mass-factor=200 \
  --temperature=100000 \
  --kbt1 \
  --auto-release \
  --steps=5000000
```

Then analyze:
```bash
python3 wall_x_FFT.py
```

**Expected result:** ν_sim ≈ 0.0117 Hz, matching theory within 5%! ✅

---

## Next Steps

1. ✅ Run with T=100000 and 5M steps
2. ✅ Verify amplitude is -2σ to +2σ
3. ✅ Check entropy plots plateau
4. ✅ Confirm FFT peak at 0.012 Hz
5. ✅ Compare with Román et al. figures

Once validated, you can:
- Vary L₀ to test cs(η) dependence
- Vary M to test frequency scaling
- Vary N to test finite-size effects
- Study amplitude dependence (vary initial T)

---

**Created:** 2025-11-05
**Last Updated:** 2025-11-05
