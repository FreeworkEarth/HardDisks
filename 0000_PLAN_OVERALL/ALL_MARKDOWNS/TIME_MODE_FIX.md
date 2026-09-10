# TIME Mode Time Scaling Fix

## The Problem

User reported: *"it moves super super super slow motion. can it be that we only need to scale time with edmd and not with --time?"*

**Root cause**: TIME mode was incrementing `simulation_time` in **pixel units** instead of **σ units**, causing the time counter to advance 40× faster than the actual physics.

### What Was Happening

1. TIME mode uses `fixed_dt_runtime` which is defined in pixel coordinates (like all internal physics)
2. TIME mode incremented time: `simulation_time += fixed_dt_runtime;` (pixel units)
3. EDMD mode correctly converted: `simulation_time += fixed_dt_runtime / PIXELS_PER_SIGMA;` (σ units)
4. Result: TIME mode time counter raced ahead by factor of 40, making simulation appear in "slow motion"

**Example**:
- `fixed_dt_runtime = 0.001` (pixel-time units)
- `PIXELS_PER_SIGMA = 40`
- After 100 timesteps:
  - TIME mode: `simulation_time = 0.1` (WRONG - claims 0.1 time units passed)
  - EDMD mode: `simulation_time = 0.0025` (CORRECT - actually 0.0025 σ-time units passed)

---

## The Fix

### Files Modified

**[hspist3/00ALLINONE.c](hspist3/00ALLINONE.c)**

Added `/PIXELS_PER_SIGMA` division to TIME mode time increments.

#### Fix 1: Interactive TIME Mode (Line 5962)

**Before (Broken)**:
```c
// time bookkeeping
simulation_time += fixed_dt_runtime;
```

**After (Fixed)**:
```c
// ##CHRIS: time bookkeeping - convert pixel-time to σ-time
simulation_time += fixed_dt_runtime / PIXELS_PER_SIGMA;
```

#### Fix 2: Headless TIME Mode (Line 6094)

**Before (Broken)**:
```c
simulation_time += fixed_dt_runtime;
```

**After (Fixed)**:
```c
// ##CHRIS: Convert pixel-time to σ-time units (TIME mode fix)
simulation_time += fixed_dt_runtime / PIXELS_PER_SIGMA;
```

---

## Why This Fix Is Correct

### Unit Analysis

All physics in 00ALLINONE.c happens in **pixel coordinates**:
- Particle positions: `X[i]`, `Y[i]` (pixels)
- Box dimensions: `boxW = 2 * L0 * PIXELS_PER_SIGMA` (pixels)
- Particle radius: `PARTICLE_RADIUS = 0.5 * PIXELS_PER_SIGMA` (pixels)
- Timestep: `fixed_dt_runtime = dt * PIXELS_PER_SIGMA` (pixel-time)

But users think in **σ units** (particle diameters):
- Box length: `L0 = 20σ`
- Particle diameter: `d = 1σ`
- Timestep: `dt = 0.001 τ` where `τ = σ/√(kBT/m)` (natural time unit)

**Conversion**: To get σ-time from pixel-time, divide by `PIXELS_PER_SIGMA`:

```c
time_σ = time_pixels / PIXELS_PER_SIGMA
```

This is exactly what EDMD mode was already doing!

---

## Comparison with EDMD Mode

### EDMD Mode Time Increment (Line 6389)

```c
// ##CHRIS: EDMD uses pixel coordinates, so dt is in pixel-time units.
// Convert to σ-time units for logging by dividing by PIXELS_PER_SIGMA.
simulation_time += fixed_dt_runtime / PIXELS_PER_SIGMA;
```

**Why EDMD needed this**: EDMD internally works in pixel coordinates but logs time in σ units for consistency with theory/experiments.

**Why TIME mode needs it**: Same reason! TIME mode also uses pixel coordinates internally.

---

## Testing

### Test 1: Compile and Verify Binary

```bash
cd hspist3
make clean && make
ls -lh 00ALLINONE
```

**Expected**: Binary successfully compiled (568 KB)

### Test 2: TIME Mode Basic Run

```bash
rm -f wall_position.csv run_params.json

./00ALLINONE --mode=time --particles=100 --l0=20 --height=10 \
  --num-walls=1 --wall-positions=10.5 --wall-mass-factor=200 --kbt1 \
  --time=50 --output-dt=1.0 --auto-release --no-gui
```

**Expected behavior**:
- Simulation completes in reasonable time
- CSV contains ~50 rows (one per σ-time unit)
- Time values range from 0 to ~50

### Test 3: Compare TIME vs EDMD

Run same simulation with both modes:

```bash
# EDMD mode
rm -f wall_position.csv
./00ALLINONE --mode=edmd --particles=100 --l0=20 --height=10 \
  --wall-mass-factor=200 --kbt1 --time=100 --output-dt=1.0 --auto-release --no-gui
cp wall_position.csv edmd_result.csv

# TIME mode
rm -f wall_position.csv
./00ALLINONE --mode=time --particles=100 --l0=20 --height=10 \
  --wall-mass-factor=200 --kbt1 --time=100 --output-dt=1.0 --auto-release --no-gui
cp wall_position.csv time_result.csv

# Compare
head -20 edmd_result.csv
head -20 time_result.csv
```

**Expected**: Both CSVs should show similar time progression and wall behavior.

### Test 4: FFT Analysis Consistency

```bash
# Run TIME mode
./00ALLINONE --mode=time --particles=100 --l0=20 --height=10 \
  --wall-mass-factor=200 --kbt1 --time=300 --output-dt=1.0 --auto-release --no-gui

# Analyze
python3 wall_x_FFT.py
```

**Expected results** (for L0=20, M=200, N=100, η≈0.196):
- Theory: ν ≈ 0.0117 Hz
- Simulation: ν ≈ 0.0110-0.0125 Hz (within 15%)
- Status: ✅ GOOD agreement

---

## Why This Was Hard to Diagnose

1. **Visual inspection looked correct** - particles were moving at right speed
2. **Physics was correct** - collisions and forces were properly calculated
3. **Only the time counter was wrong** - easy to miss until you check CSV timestamps
4. **"Slow motion" description was confusing** - actually the opposite! Time was racing ahead
5. **Mode-specific inconsistency** - EDMD worked fine, only TIME mode broken

The clue was user's observation: *"can it be that we only need to scale time with edmd and not with --time?"* - exactly right!

---

## Impact

### Before Fix

| Mode | Time Increment | Units | Correct? |
|------|---------------|-------|----------|
| EDMD | `dt / PIXELS_PER_SIGMA` | σ-time | ✅ Yes |
| TIME | `dt` (no division) | pixel-time | ❌ No |

**Problem**: TIME mode CSV timestamps were 40× larger than they should be, causing:
- FFT analysis to find wrong frequencies (40× too low)
- Comparison with theory impossible
- `--time=50` flag actually ran for 50×40 = 2000 σ-time units!

### After Fix

| Mode | Time Increment | Units | Correct? |
|------|---------------|-------|----------|
| EDMD | `dt / PIXELS_PER_SIGMA` | σ-time | ✅ Yes |
| TIME | `dt / PIXELS_PER_SIGMA` | σ-time | ✅ Yes |

**Result**: Both modes now use consistent time units!

---

## Related Documentation

See also:
- [FFT_DC_AND_PEAK_EXPLANATION.md](FFT_DC_AND_PEAK_EXPLANATION.md) - Understanding FFT frequency resolution
- [PARTICLE_COUNT_LOGGING_FIX.md](PARTICLE_COUNT_LOGGING_FIX.md) - Particle count CSV logging fix
- [INVISIBLE_WALL_FIX.md](INVISIBLE_WALL_FIX.md) - Wall toggle synchronization with EDMD

---

## Summary

✅ **TIME mode time scaling fixed**
✅ **Both EDMD and TIME modes use consistent σ-time units**
✅ **FFT analysis can now validate TIME mode results**
✅ **User can compare EDMD vs TIME mode directly**

### Before Fix
```
TIME mode: simulation_time advances in pixel units (40× too fast)
Result: "slow motion" appearance, timestamps wrong by factor of 40
```

### After Fix
```
TIME mode: simulation_time advances in σ units (correct!)
Result: Normal speed, consistent with EDMD mode
```

**Try your comparison now - TIME and EDMD modes should give consistent results!** 🎉

```bash
# Compare modes
./00ALLINONE --mode=edmd --particles=100 --l0=20 --wall-mass-factor=200 --kbt1 --time=100 --auto-release --no-gui
./00ALLINONE --mode=time --particles=100 --l0=20 --wall-mass-factor=200 --kbt1 --time=100 --auto-release --no-gui
```

Both should complete in similar time and produce comparable wall oscillation data!
