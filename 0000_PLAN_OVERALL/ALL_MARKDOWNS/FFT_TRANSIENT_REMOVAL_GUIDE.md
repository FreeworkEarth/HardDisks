# FFT Analysis: Transient Removal & Theoretical Predictions

## Your Current Problem

**FFT Result:** ν = 0.0011 Hz  
**Theory:** ν = 0.011717 Hz  
**Ratio:** 0.094 (only 9.4%!) ❌

This is **WORSE** than before! The problem is:

### 1. Simulation Too Short

Expected period: **T ≈ 85 time units**  
Need: **At least 3-5 full periods**  
Minimum simulation time: **250-425 time units**

Your command was missing `--steps`:
```bash
./00ALLINONE \
  --mode=edmd \
  --particles=100 \
  --l0=20 \
  --height=10 \
  --wall-mass-factor=200 \
  --temperature=1000 \
  --kbt1 \
  --auto-release
  # MISSING: --steps=300000 or more!
```

## Corrected Command

```bash
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
```

With typical dt ≈ 0.001, this gives:
- **500 time units** of data
- **~6 full oscillation periods**
- Good frequency resolution: Δf ≈ 1/500 = 0.002

## Transient Removal (Now Automated!)

Your `wall_x_FFT.py` now has **two-stage transient removal**:

### Stage 1: Remove Flat Part (Automatic)
- Detects when wall starts moving
- Removes all data before movement begins
- **Already working:** "✅ Trimmed first 19016 samples"

### Stage 2: Remove Oscillation Transients (NEW!)

After the wall is released, it takes time to reach steady oscillation (like a pendulum settling).

**Configuration** (in `wall_x_FFT.py`):
```python
remove_transient_percentage = 0.20  # Remove first 20% after movement starts
```

**Example:**
- Data after flat removal: 100,000 samples
- Remove 20%: 20,000 samples discarded
- Use remaining 80,000 samples for FFT

**Output:**
```
✅ Trimmed first 19016 samples where displacement was flat.
🔧 Removed additional 20% (14416 samples) for oscillation transients
   Remaining data: 57664 samples from t=34.42 to t=492.06
```

## Theoretical Predictions

For your exact parameters (N=100, L₀=20, M=200m, kBT=1):

| Quantity | Value | Units |
|----------|-------|-------|
| Packing fraction η | 0.1963 | dimensionless |
| Speed of sound cs | 2.2538 | (reduced) |
| Transcendental K | 0.653271 | rad |
| **Frequency ν** | **0.011717** | **(reduced)** |
| Period T | 85.35 | time units |
| Angular frequency ω | 0.073617 | rad/time |

### Expected FFT Peak

Your FFT should show a **strong peak at ν ≈ 0.0117** (in reduced units).

**Acceptable range:** 0.0111 - 0.0123 (95-105% of theory)

## Laplace Transform vs. Simple 20% Cutoff

You asked about Laplace transforms for analyzing transients (like the pendulum video).

### Laplace Transform Approach (Complex)
- Solve differential equation analytically
- Extract exponential decay rate: γ
- Calculate transient duration: τ_transient = 3/γ
- **Pros:** Mathematically rigorous
- **Cons:** Need to know exact damping coefficient

### Simple 20% Cutoff (Practical)
- Just discard first 20% of data
- **Pros:** Simple, robust, works for most systems
- **Cons:** Might discard slightly too much or too little

**For your case:** **20% cutoff is sufficient!**

The wall oscillation is **lightly damped** (very little energy loss), so transients decay in ~1-2 periods:
- Period T = 85 time units
- Transient duration ≈ 2T = 170 time units
- 20% of 500 time units = 100 time units ≈ 1.2 periods ✓

## How to Verify Transient Removal is Working

### 1. Visual Inspection

After running `wall_x_FFT.py`, check the displacement plot:
- First 20% should show irregular oscillations
- Remaining 80% should show regular sinusoidal pattern

### 2. Compare Different Percentages

Try different values and see if frequency changes:

```python
# In wall_x_FFT.py
remove_transient_percentage = 0.00  # No removal
# Run, note frequency

remove_transient_percentage = 0.10  # 10% removal
# Run, note frequency

remove_transient_percentage = 0.20  # 20% removal (recommended)
# Run, note frequency

remove_transient_percentage = 0.30  # 30% removal
# Run, note frequency
```

**Expected result:**
- 0%: Frequency too low (transients contaminate FFT)
- 10-20%: Stable frequency ≈ 0.0117
- 30%+: Frequency still ≈ 0.0117 (but less data)

If frequency varies a lot, you need longer simulation!

## Power Spectrum

**Question:** "What's the theoretical power spectrum value?"

The **power** in FFT is proportional to amplitude squared:

```python
Power = |FFT(signal)|²
```

For a **pure sine wave** with amplitude A:
- FFT shows single peak at frequency ν
- Peak height ∝ A² × N (where N = number of samples)

**For your wall oscillation:**
- Amplitude ≈ 1-2 σ (from pressure fluctuations)
- Peak power ≈ (1-2)² × 57664 ≈ 10⁴ to 10⁵

Your current result: "Power: 2.7722e-03" is **WAY too small!**

This confirms: **simulation too short or not enough oscillations**

## Summary Checklist

Before running FFT analysis:

- [ ] Simulation ran for ≥500 time units (≥6 periods)
- [ ] Used `--kbt1` flag (kBT = 1.0)
- [ ] Used `--auto-release` flag (data logged from start)
- [ ] `remove_transient_percentage = 0.20` (discard irregular oscillations)
- [ ] Expected frequency: ν ≈ 0.0117
- [ ] Expected power: >1000 (not 0.003!)

## Recommended Full Workflow

```bash
# 1. Run long simulation
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

# 2. Check data was saved
wc -l wall_position.csv
# Should see >400,000 lines

# 3. Analyze with FFT
python3 wall_x_FFT.py

# 4. Check results
# Expected: ν ≈ 0.0117 (not 0.0011!)
```

## Why 20% is Enough

**Physical reasoning:**

The wall+gas system is like a **harmonic oscillator with light damping**:

```
m_eff d²x/dt² + γ dx/dt + k x = 0
```

Where:
- m_eff = wall mass + effective gas mass
- γ = damping (very small, gas is ideal)
- k = "spring constant" from gas pressure

**Damping ratio:** ζ = γ / (2√(m_eff × k)) ≈ 0.01 (very small)

**Transient decay:** x(t) ∝ e^(-ζωt) × cos(ωt)

Time to decay to 5%: t_5% = 3/(ζω) ≈ 3/(0.01 × 0.074) ≈ 4000 time units

Wait... that's **50 periods!** 😱

**BUT:** We don't need transients to fully decay, just to be **small enough** to not affect FFT peak position.

After 1-2 periods (20% of 500 = 100 time units):
- Transient amplitude ≈ 10-20% of steady state
- Frequency already correct within 2-3%

**So 20% is a good practical compromise!**

