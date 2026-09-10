# Quick Start Guide: Entropy-Based Hard Disk Simulation

## TL;DR - Run This Now

```bash
cd hspist3

# 1. Run simulation (5 minutes)
./00ALLINONE \
  --mode=edmd \
  --particles=100 \
  --l0=20 \
  --height=10 \
  --wall-mass-factor=200 \
  --temperature=1000 \
  --kbt1 \
  --auto-release \
  --steps=500000

# 2. Analyze (10 seconds)
python3 wall_x_FFT.py

# 3. Check results
# Look for: "✅ EXCELLENT: Simulation matches theory within 5%!"
# Open: divider_x_displacement_with_boxlength_20_and_wallfactor_200.pdf
```

---

## What Just Happened?

### Simulation (`00ALLINONE`)
- Created 100 hard disks in a 2D box with oscillating wall
- Automatically released wall at t=0 (`--auto-release`)
- Ran for 500,000 steps ≈ 500 time units ≈ 6 oscillation periods
- Saved wall position to `wall_position.csv`

### Analysis (`wall_x_FFT.py`)
1. **Removed flat part** (before wall starts moving)
2. **Entropy detection** - Found when S[p(x)] and S[p(v)] stabilize
3. **Removed transients** - Typically 10-20% of data
4. **FFT analysis** - Extracted fundamental frequency ν
5. **Theory check** - Compared with ν = (cs/2πL₀)K prediction
6. **Generated plots:**
   - 3-panel: displacement + position entropy + velocity entropy
   - FFT: power spectrum with peak marked

---

## What To Look For

### ✅ Good Run (Success!)

**Terminal output:**
```
✅ Position entropy stable at 5000 samples (8.3%)
✅ Velocity entropy stable at 8000 samples (13.3%)
🎯 Using transient fraction: 13.3%

Highest Power Spectrum Peak: 0.0117 Hz

✅ EXCELLENT: Simulation matches theory within 5%!
Ratio (sim/theory): 0.9985
```

**In plots:**
- Both entropy curves plateau (flat region after rise)
- Clear oscillations in displacement after gray region
- Sharp FFT peak at ν ≈ 0.0117

### ❌ Bad Run (Needs Fixing)

**Terminal output:**
```
⚠️ Could not detect clear steady state (entropy unstable)

❌ WARNING: Large deviation from theory!
Ratio (sim/theory): 0.094
```

**In plots:**
- Entropy curves still rising (no plateau)
- Irregular displacement pattern
- Noisy FFT spectrum

**Fix:** Run longer!
```bash
./00ALLINONE ... --steps=1000000  # Double the time
```

---

## Key Features (New!)

### 1. Auto-Release (`--auto-release`)
- No need to press 'R' key manually
- Wall starts moving immediately
- Logs position from t=0

### 2. Entropy Detection
- **Position entropy** S[p(x)] → spatial equilibration
- **Velocity entropy** S[p(v)] → thermal equilibration
- **Automatic:** Removes transients when both stabilize
- **Adaptive:** Not fixed 20% - based on actual physics!

### 3. Theory Validation
- Calculates expected ν from Henderson equation
- Compares simulation vs. theory
- Verdict: Excellent/Good/Acceptable/Warning
- Catches common mistakes (wrong kBT, too short, etc.)

### 4. Three-Panel Visualization
- **Top:** Wall displacement (gray = transient)
- **Middle:** Position entropy over time
- **Bottom:** Velocity entropy over time
- Red line marks steady-state transition

---

## Troubleshooting

### "Overpacked: requested total=100 > capacity total=98"
**Fix:** Recompile (capacity bug was fixed)
```bash
cd hspist3
make clean
make
```

### "wall_position.csv is empty"
**Fix:** Add `--auto-release` flag to your command

### "Entropy never stabilizes"
**Fix:** Run longer
```bash
--steps=1000000  # Instead of 500000
```

### "Theory ratio = 0.094 (very low)"
**Fix:** Simulation was too short - need at least 500 time units
```bash
--steps=500000  # Make sure you have this!
```

### "Theory ratio varies: 0.98, 1.02, 0.99..."
**Normal:** Statistical noise with N=100
**Options:**
- Average multiple runs
- Use more particles (`--particles=200`)
- Run longer (`--steps=1000000`)

---

## Understanding the Output

### Entropy Plots (Physical Meaning)

**Position Entropy S[p(x)] = -Σ p(x) log p(x)**
- **Rising:** Particles spreading from initial random positions
- **Plateau:** Equilibrium spatial distribution reached
- **Physics:** Maximum entropy = uniform distribution (for hard disks)

**Velocity Entropy S[p(v)] = -Σ p(v) log p(v)**
- **Rising:** Collisions redistributing kinetic energy
- **Plateau:** Maxwell-Boltzmann distribution reached
- **Physics:** Maximum entropy = thermal equilibrium at temperature T

**Why both?**
- Position equilibrates quickly (few collision times)
- Velocity equilibrates slower (energy exchange takes time)
- Need BOTH stable for true steady state!

### Theory Comparison

**Frequency formula:**
```
ν = (cs / 2πL₀) K
```

Where:
- **cs** = speed of sound (from Henderson equation, depends on η)
- **L₀** = box half-length (20σ in our case)
- **K** = transcendental root satisfying cot(K) = (M/2Nm)K

**Your parameters:**
- η = 0.1963 (packing fraction)
- cs = 2.254 (speed of sound)
- K = 0.653 rad (from transcendental equation)
- **Expected:** ν = 0.0117 reduced units

**Simulation should match within 5%!**

---

## Parameter Sweep (Next Step)

Once validation succeeds, test different box sizes:

```bash
cd hspist3

for L0 in 10 15 20 25 30; do
  echo "=== Running L0=$L0 ==="

  ./00ALLINONE \
    --mode=edmd \
    --particles=100 \
    --l0=$L0 \
    --height=10 \
    --wall-mass-factor=200 \
    --temperature=1000 \
    --kbt1 \
    --auto-release \
    --steps=500000

  python3 wall_x_FFT.py

  echo "=== Completed L0=$L0 ==="
  echo ""
done
```

This generates data for:
- Speed of sound vs. packing fraction (cs vs η)
- Frequency vs. box size (ν vs L₀)
- Entropy equilibration times
- Validation across parameter space

---

## Files Generated

### Every Run
- `wall_position.csv` - Wall x-position vs time
- `experiments_energy_transfer/particle_states_*.csv` - Full particle data

### Every Analysis
- `divider_x_displacement_with_boxlength_*.pdf` - 3-panel entropy plot
- `Power_Freq_Spectrum_with_boxlength_*.pdf` - FFT spectrum
- `divider_x_displacement_with_boxlength_*.png` - PNG version
- `Power_Freq_Spectrum_with_boxlength_*.png` - PNG version

---

## Documentation

**For details, see:**

| Topic | File |
|-------|------|
| Complete implementation summary | [IMPLEMENTATION_COMPLETE.md](IMPLEMENTATION_COMPLETE.md) |
| Latest updates | [FINAL_UPDATES_ENTROPY_AND_THEORY.md](FINAL_UPDATES_ENTROPY_AND_THEORY.md) |
| Entropy method details | [hspist3/ENTROPY_STEADY_STATE_DETECTION.md](hspist3/ENTROPY_STEADY_STATE_DETECTION.md) |
| Auto-release feature | [hspist3/AUTO_RELEASE_WALL_FEATURE.md](hspist3/AUTO_RELEASE_WALL_FEATURE.md) |
| Units and timing | [UNITS_AND_TIMING_EXPLAINED.md](UNITS_AND_TIMING_EXPLAINED.md) |
| Validation checklist | [VALIDATION_CHECKLIST.md](VALIDATION_CHECKLIST.md) |

---

## Success Criteria

Your simulation is **validated** when:

- ✅ Both entropies plateau (not still rising)
- ✅ Transient removal ~10-20% (not 0% or 50%)
- ✅ Clear FFT peak (not noise)
- ✅ Theory ratio 0.95-1.05 ("✅ EXCELLENT" message)
- ✅ Displacement shows regular sinusoidal pattern in steady-state region

If all boxes checked → **Your simulation is working correctly!**

Use the data for:
- Publication figures (entropy plots show equilibration)
- Theory validation (compare cs measurements with Henderson equation)
- Parameter studies (sweep L₀, M, N to test predictions)

---

## Quick Reference

### Command Flags
```bash
--mode=edmd              # Event-driven MD (required)
--particles=100          # Total particles (50 per side)
--l0=20                  # Box half-length in σ
--height=10              # Box height in σ
--wall-mass-factor=200   # M = 200m
--temperature=1000       # Initial temperature
--kbt1                   # Scale to kBT=1 (required!)
--auto-release           # Auto-release wall (new!)
--steps=500000           # Run for 500 time units
```

### Analysis Configuration
In `wall_x_FFT.py` (lines 59-62):
```python
use_entropy_detection = True  # Use entropy method (recommended)
remove_transient_percentage = 0.20  # Fallback if entropy fails
```

---

**Last Updated:** 2025-11-05
**Status:** ✅ All features working and validated
