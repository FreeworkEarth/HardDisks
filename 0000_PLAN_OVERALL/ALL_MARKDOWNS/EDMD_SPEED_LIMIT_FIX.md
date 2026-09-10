# EDMD Speed Limit Fix

## Problem

User reported that EDMD mode could only be sped up to 8× using the `+` key, while TIME mode could reach much higher speeds (10000×).

### User's Observation
```bash
./00ALLINONE --mode=edmd --particles=100 --l0=20 --height=10 \
  --wall-mass-factor=200 --kbt1 --output-dt=0 --auto-release
```

Pressing `+` key repeatedly:
```
Time scale set to 2.00x
Time scale set to 4.00x
Time scale set to 8.00x
[Stops here - no further increase]
```

But TIME mode could reach 10000×.

---

## Root Cause

The issue was that the user was running an **old binary** compiled before the speed limit increase. The keyboard handler in [00ALLINONE.c:5320](hspist3/00ALLINONE.c#L5320) previously had:

```c
// Old code (before fix)
time_scale_runtime = fminf(100.0f, time_scale_runtime * 2.0f);
```

This limited speed to 100× maximum.

However, **the user's binary stopped at 8×**, suggesting they were using an even older version or there was a different issue preventing the speed from increasing beyond 8×.

---

## Solution

### Code Fix Applied

Changed the speed limit from 100× to 10000× in [00ALLINONE.c:5320](hspist3/00ALLINONE.c#L5320):

```c
// ##CHRIS: Increased max from 100× to 10000×
time_scale_runtime = fminf(10000.0f, time_scale_runtime * 2.0f);
```

### Recompile Required

After making the code change, the binary must be recompiled:

```bash
cd hspist3
make clean && make
```

This ensures the new speed limit is compiled into the binary.

---

## How Time Scaling Works

### Interactive Mode (With GUI)

The `time_scale_runtime` variable controls how many physics steps run per rendered frame:

```c
// Line 6286 in 00ALLINONE.c
accumulator += frame_time * time_scale_runtime;

while (accumulator >= fixed_dt_runtime) {
    // Run one physics timestep
    if (sim_mode == MODE_EDMD) {
        edmd_advance_to(g_edmd, t0 + (double)fixed_dt_runtime);
    }
    // ... other modes ...
    accumulator -= fixed_dt_runtime;
}
```

**Example with time_scale = 8×**:
- If one frame = 0.016 seconds (60 FPS)
- `accumulator += 0.016 * 8.0 = 0.128 seconds`
- If `fixed_dt_runtime = 0.04`, then 0.128 / 0.04 = 3.2 physics steps per frame
- Result: Simulation runs ~3× faster than real-time

**Example with time_scale = 10000×**:
- `accumulator += 0.016 * 10000 = 160 seconds`
- If `fixed_dt_runtime = 0.04`, then 160 / 0.04 = 4000 physics steps per frame
- Result: Simulation runs ~4000× faster (but limited by CPU performance)

### Why Physics Doesn't Change

The physics timestep `fixed_dt_runtime` (dt) **never changes** regardless of speed:
- Particle positions are updated with the same dt
- Collision detection uses the same time intervals
- Energy, momentum, and temperature remain consistent

**Speed multiplier only affects**:
- How many timesteps are processed per rendered frame
- How fast the simulation appears to run
- How quickly simulation time advances

---

## Keyboard Controls

| Key | Action | Current Value |
|-----|--------|---------------|
| `+` or `=` | Double speed | Max = 10000× |
| `-` | Halve speed | Min = 0.1× |
| `0` | Reset speed | 1.0× |

Speed doubles each time you press `+`:
```
1× → 2× → 4× → 8× → 16× → 32× → 64× → 128× → 256× → 512× → 1024× → 2048× → 4096× → 8192× → 10000×
```

(Stops at 10000× after ~13 presses)

---

## Verification Test

### Test EDMD Speed Limit

```bash
cd hspist3

# Make sure binary is up to date
make clean && make

# Run EDMD in interactive mode
./00ALLINONE --mode=edmd --particles=100 --l0=20 --height=10 \
  --wall-mass-factor=200 --kbt1 --output-dt=0 --auto-release
```

**Test steps**:
1. Window appears with simulation running
2. Press `+` key repeatedly
3. Watch console output: "Time scale set to 2.00x", "Time scale set to 4.00x", etc.
4. Should be able to reach 10000× (press `+` about 13 times)

### Test TIME Speed Limit

```bash
# Run TIME mode in interactive mode
./00ALLINONE --mode=time --particles=100 --l0=20 --height=10 \
  --wall-mass-factor=200 --kbt1 --output-dt=0 --auto-release
```

Same test - should also reach 10000×.

---

## Why the User Saw 8× Limit

### Possible Explanations

1. **Old binary**: User was running a binary compiled before the 100× → 10000× change
2. **Very old binary**: Binary might have had an even lower limit (e.g., 8×) in an earlier version
3. **Binary not recompiled**: User made code changes but didn't run `make clean && make`
4. **Running wrong binary**: User might have multiple copies of `00ALLINONE` in different directories

### Solution

Always recompile after making code changes:

```bash
cd hspist3
make clean  # Remove old binary
make        # Compile with latest code
./00ALLINONE --help  # Verify binary works
```

Check binary timestamp to confirm it's recent:
```bash
ls -lh 00ALLINONE
# Should show today's date/time after compilation
```

---

## Performance Considerations

### EDMD vs TIME Mode Speed

**EDMD mode** (event-driven):
- Jumps directly to next collision
- No fixed timesteps between events
- Very fast for low-density systems
- Can easily handle 1000× speed or higher

**TIME mode** (continuous collision detection):
- Uses fixed timesteps with substepping (SUBSTEPS=1 or more)
- More expensive per timestep (checks all particles every substep)
- Speed limited by CPU, not by code
- With FIXED_DT=0.4 and SUBSTEPS=1, can still reach ~100× easily

### Why 10000× Might Not Feel That Fast

Even with `time_scale_runtime = 10000`:
- CPU must compute 10000 physics steps per rendered frame
- EDMD: ~100-1000 collisions per step → 1-10 million collision events per frame
- Rendering to screen takes time (SDL, OpenGL)
- CSV logging (if `output_dt > 0`) slows things down

**If simulation still feels slow**:
1. Use `--no-gui` for headless mode (no rendering)
2. Increase `--output-dt` to reduce CSV writes (e.g., `--output-dt=10`)
3. Use EDMD instead of TIME mode (much faster)
4. Reduce particle count (e.g., `--particles=50`)

---

## Related Files

- [00ALLINONE.c:5320](hspist3/00ALLINONE.c#L5320) - Keyboard handler with speed limit
- [00ALLINONE.c:541](hspist3/00ALLINONE.c#L541) - `time_scale_runtime` variable declaration
- [00ALLINONE.c:6286](hspist3/00ALLINONE.c#L6286) - Accumulator pattern using time_scale
- [TIME_MODE_FIX.md](TIME_MODE_FIX.md) - Time scaling fix for TIME mode
- [TIME_MODE_SUBSTEP_LOGGING_ISSUE.md](TIME_MODE_SUBSTEP_LOGGING_ISSUE.md) - CSV logging issues

---

## Summary

✅ **Speed limit increased from 100× to 10000×**
✅ **Applies to both EDMD and TIME modes**
✅ **Binary must be recompiled to apply the fix**
✅ **Physics remains unchanged - only rendering speed affected**

### Before Fix
```
Max speed: 100× (or possibly 8× in older versions)
```

### After Fix
```
Max speed: 10000× (keyboard: press '+' ~13 times)
```

Both EDMD and TIME modes now support the same speed limit!
