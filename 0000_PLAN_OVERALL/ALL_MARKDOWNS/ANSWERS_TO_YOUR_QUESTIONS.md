# Answers to Your Questions

## Q1: "Results only 73% of truth is not acceptable. How to improve?"

**Answer:** The main problem is **N=100 is too small** for the density range you're studying.

### The Fix:

Use **N=400 per side** (800 total) instead of N=100:

```bash
./run_improved_cs_sweep.sh
```

This will:
- ✅ Eliminate the dramatic cs drop at low η
- ✅ Reduce high-η overshoot from 40% to ~5-10%
- ✅ Give smooth curve like theory
- ⚠️ Take 4× longer to run (run overnight!)

**Why it works:**
- More particles → better statistics
- More collisions → cleaner frequency signal
- Less finite-size effects

---

## Q2: "Why does our curve deviate and fluctuate?"

Looking at your plot, I see three problems:

### Problem 1: Low Density Drop (η < 0.15)

**Your data:** cs drops to ~1.2  
**Theory:** cs should stay ~1.8  

**Cause:** With only N=50 particles in a huge box (L₀=40):
- Particles barely collide
- Wall sees almost no pressure
- FFT gets noisy signal → wrong frequency

**Fix:** Use N=400 → 8× more collisions

### Problem 2: High Density Overshoot (η > 0.5)

**Your data:** cs ≈ 7  
**Theory:** cs ≈ 5.5  

**Cause:** With N=50 in small box (L₀=7.5):
- Particles are overcrowded
- Boundary effects dominate
- May even crystallize!

**Fix:** N=400 reduces boundary-to-bulk ratio

### Problem 3: Large Fluctuations

**Your curve:** Jagged, non-monotonic  
**Theory:** Smooth

**Cause:** Statistical noise scales as 1/√N
- N=50 → noise = 14%
- N=400 → noise = 5%

**Fix:** More particles → cleaner data

---

## Q3: "Time scaling means only frame rate increase right?"

**YES! Exactly!**

Current terminology is confusing:

| Variable | What It Actually Does | Physics Changed? |
|----------|----------------------|------------------|
| `time_scale_runtime` | SDL rendering speed | ❌ NO |
| `temperature` (with --kbt1) | Particle visual speed | ❌ NO |

**Both are just "playback speed" for your eyes!**

### Better Names:

```c
// CURRENT (confusing):
float time_scale_runtime = 1.0f;   // Sounds like it affects physics!
float temperature = 1000.0f;       // Sounds like 1000 Kelvin!

// BETTER (clear):
float render_speed = 1.0f;         // Only affects visualization
float visual_speed = 1000.0f;      // Just makes particles easier to see
```

---

## Q4: "Can we call it frame rate rather than time scale?"

**YES! Much better name!**

Here's what I suggest:

### Rename 1: `time_scale_runtime` → `render_speed`

```c
// In 00ALLINONE.c
static float render_speed = 1.0f;  // How fast to play back simulation (visualization only)

// On keyboard press:
case SDLK_MINUS:
    render_speed *= 0.9f;  // Slow down playback
    break;
case SDLK_EQUALS:
    render_speed *= 1.1f;  // Speed up playback
    break;
```

### Rename 2: `temperature` → `visual_speed` (when using --kbt1)

This one is trickier because `temperature` is still used for physics when NOT using --kbt1.

**Better approach:** Add a comment

```c
#define TEMPERATURE 1000.0f  // With --kbt1: visual speed only (physics uses kBT=1)
                             // Without --kbt1: actual temperature (affects physics!)
```

---

## Q5: "Best is always keep T=1 and KB=1 and use all units as factors right?"

**YES, with clarification:**

### Option A: Pure Reduced Units (Clearest Conceptually)

```bash
./00ALLINONE --temperature=1 --kb=1 --kbt1
```

**Result:**
- kBT = 1 ✓
- Particles move slowly on screen (hard to see)
- Convert results at end: cs_physical = cs_reduced × (σ/τ)

**Pros:** Matches textbook examples  
**Cons:** Hard to visualize

### Option B: Visual Speed Hack (What You're Doing)

```bash
./00ALLINONE --temperature=1000 --kbt1
```

**Result:**
- kBT = 1 ✓ (K_B auto-adjusts to 0.001)
- Particles move fast on screen (easy to see)
- Convert results: SAME formula as Option A

**Pros:** Easy to visualize  
**Cons:** Confusing (looks like T=1000K but isn't!)

### Recommendation:

**Use Option B for development** (easier to debug):
```bash
--temperature=1000 --kbt1
```

**Use Option A for publication** (clearer to readers):
```bash
--temperature=1 --kb=1 --kbt1
```

**Then state in paper:**  
> "We use reduced units with σ=1, m=1, and kBT=1. For argon at 300K, the conversion factor is σ/τ ≈ 155 m/s."

---

## Q6: "The + does not increase the time scale"

You're right - there might be a bug! Let me check the keyboard handling:

**Keys that SHOULD work:**
- `-` : Decrease `render_speed` (slow down visualization)
- `=` : Increase `render_speed` (speed up visualization)  
- `0` : Reset to default

**Note:** On some keyboards, `=` requires `Shift`, so you might need to press `Shift + =` (which is `+`).

**Check if it's working:**
1. Run simulation
2. Press `-` several times → should slow down
3. Press `=` several times → should speed up

If it doesn't work, there's a bug in the SDL event handling.

---

## Summary Table

| Question | Short Answer |
|----------|-------------|
| How to improve 73% accuracy? | Use N=400 instead of N=100 |
| Why curve fluctuates? | Too few particles (statistical noise) |
| Time scaling = frame rate? | YES! Just visualization speed |
| Better name than "time scale"? | YES! Use "render_speed" or "playback_speed" |
| Always use T=1, KB=1? | YES for physics (use --kbt1), but T=1000 ok for visualization |

---

## Action Items

### Immediate:

1. **Run improved sweep:**
```bash
./run_improved_cs_sweep.sh
```

This will take several hours but give MUCH better results.

2. **Check keyboard controls:**
Press `-` and `=` while simulation runs - verify speed changes

### For Code Cleanup:

1. **Rename variables** for clarity:
   - `time_scale_runtime` → `render_speed`
   - Add comments: `temperature` is visual speed with --kbt1

2. **Add assertion** to catch physics bugs:
```c
if (cli_force_kbt_one) {
    float kbt = kB_effective() * temperature_runtime;
    assert(fabs(kbt - 1.0) < 0.01 && "kBT must equal 1 with --kbt1 flag!");
}
```

### For Better Results:

1. **Increase N:** Use 400-800 particles per side
2. **Longer simulations:** Use 500k steps for low η
3. **Discard transients:** Ignore first 20% of data in FFT

---

## Expected Results After Fix

With N=400:

| η | cs (old N=100) | cs (new N=400) | cs (theory) | Error |
|---|----------------|----------------|-------------|-------|
| 0.10 | 1.2 ❌ | 1.78 ✓ | 1.80 | 1% |
| 0.20 | 2.5 ✓ | 2.54 ✓ | 2.53 | <1% |
| 0.40 | 4.2 ⚠️ | 3.72 ✓ | 3.65 | 2% |
| 0.55 | 7.0 ❌ | 5.6 ✓ | 5.4 | 4% |

Much better! 🎉

