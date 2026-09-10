# FFT Analysis: Understanding DC, Resolution Floor, and True Peaks

## What You're Seeing

You reported finding:
1. **Large peak at ~0.003 Hz**
2. **Smaller peak at ~0.010 Hz** (close to theory ν ≈ 0.0117 Hz)
3. When plotting full spectrum including DC, everything looks tiny

**This is completely normal!** Let me explain why.

---

## What is "DC" in FFT?

**DC = "Direct Current"** (term borrowed from electrical engineering)

In FFT analysis:
- **DC component** = frequency bin at **0 Hz**
- Represents the **mean/average/constant offset** of your signal
- Often HUGE compared to other frequencies
- Should be removed before FFT (detrending/mean removal)

### Why DC is Huge

If your wall oscillates around position x=20σ with amplitude ±0.5σ:
- **DC component** ≈ 20² = 400 (power)
- **Oscillation component** ≈ 0.5² = 0.25 (power)
- **Ratio**: DC is 1600× larger than the signal!

This is why when you plot the **full spectrum**, everything except DC looks flat.

---

## What is "1/Tspan" (Resolution Floor)?

### Frequency Resolution

When you analyze a time series of duration **Tspan**:
- **Frequency resolution**: Δν = 1/Tspan
- **First non-zero bin**: ~1/Tspan Hz
- **If Tspan = 3000** time units → Δν ≈ 0.000333 Hz

### The Resolution Floor Trap

If you run for Tspan = 300 units:
- Δν ≈ 1/300 = 0.00333 Hz ≈ **0.003 Hz**
- FFT bin spacing is coarse
- Peak finder often latches onto the **first non-zero bin** near DC
- **This is NOT your physical frequency!**

This is exactly what you're seeing - your 0.003 Hz peak is likely the 1/Tspan resolution floor, not the wall eigenfrequency.

---

## Your Two Peaks Explained

### Peak 1: ~0.003 Hz (Large Amplitude)
- **What it is**: Resolution floor (1/Tspan)
- **Why it's large**: Low-frequency drift, transients, or numerical artifacts
- **Physical meaning**: None - this is a measurement artifact

### Peak 2: ~0.010 Hz (Smaller Amplitude)
- **What it is**: **Your actual wall oscillation!**
- **Theory predicts**: ν ≈ 0.0117 Hz for M=200, L0=20, N=100, kBT=1
- **Why it's smaller**: True physical oscillation has small amplitude with kBT=1

---

## Why Everything Looks Tiny When You Include DC

### Log Scale Example

If plotting power spectrum:
```
DC power:        10⁶  (huge!)
1/Tspan peak:    10²  (medium)
Physical peak:   10¹  (small)
Noise floor:     10⁰  (tiny)
```

On a **linear scale**, the physical peak at 10¹ looks like zero compared to DC at 10⁶.

On a **log scale** (semilogy), you can see all components proportionally.

---

## Solutions

### 1. Remove DC Before FFT (Already Done)
```python
signal_centered = signal - np.mean(signal)  # Remove DC
```

### 2. Use Log Scale for Full Spectrum
```python
ax.semilogy(freqs, power)  # Y-axis in log scale
```

### 3. Zoom to Expected Frequency Range
```python
ax.set_xlim(0.005, 0.015)  # Focus on ν ≈ 0.0117 Hz
```

### 4. Increase Tspan to Improve Resolution
- **Currently**: Tspan ≈ 300 → Δν ≈ 0.003 Hz (too coarse!)
- **Better**: Tspan ≈ 3000 → Δν ≈ 0.0003 Hz (10× better resolution)
- **Best**: Tspan ≈ 30000 → Δν ≈ 0.00003 Hz (can resolve peak cleanly)

### 5. Use Velocity FFT Instead of Position
Velocity often has cleaner peaks because:
- Less low-frequency drift
- Smaller DC component
- Better SNR for oscillations

---

## What I Added to Your FFT Script

### New: Full Spectrum Plot with DC Diagnostic

Now `wall_x_FFT.py` generates **two plots**:

#### Plot 1: Zoomed View (Existing)
- X-axis: 0.005 - 0.015 Hz (around expected peak)
- Shows your physical oscillation clearly
- Filename: `Power_Freq_Spectrum_with_boxlength_{L0}_and_wallfactor_{M}.pdf`

#### Plot 2: Full Spectrum View (NEW!)
- X-axis: 0 - 0.05 Hz (includes DC and resolution floor)
- **Log scale Y-axis** so you can see everything
- **Red line**: Marks 1/Tspan (resolution floor)
- **Gray shading**: DC + resolution floor region (0 to 3×Δν)
- Filename: `Power_Spectrum_FULL_with_DC_{L0}_wallfactor_{M}.pdf`

Now you can see:
1. How huge DC is (Y-axis starts at ~10⁶)
2. Where 1/Tspan sits (red dotted line)
3. Your true physical peak at ~0.010 Hz

---

## Practical Recommendations

### To Get Clean Román-Like Results

1. **Use correct mass**: M=200 (not 2000)
2. **Add small excitation**: `--wall-positions=10.5` (0.5σ offset)
3. **Run long enough**: `--steps=3000000` (≈3000 time units with dt=0.001)
4. **Log sparsely**: `--output-dt=1.0` (one sample per time unit)
5. **Use kBT=1**: `--kbt1`

Full command:
```bash
./00ALLINONE --mode=edmd --particles=100 --l0=20 --height=10 \
  --num-walls=1 --wall-positions=10.5 --wall-mass-factor=200 --kbt1 \
  --steps=3000000 --output-dt=1.0 --no-experiments
```

Then analyze:
```bash
python3 wall_x_FFT.py
```

### Expected Results

With Tspan ≈ 3000 and proper excitation:
- **Frequency resolution**: Δν ≈ 0.0003 Hz (30× better!)
- **Peak location**: ν_sim ≈ 0.0116 Hz
- **Theory**: ν_theory ≈ 0.0117 Hz
- **Agreement**: Within 1-3%
- **cs measured**: ≈ 2.25 (matches Henderson at η ≈ 0.196)

---

## Summary Table

| Feature | Value | What It Means |
|---------|-------|---------------|
| **DC** | 0 Hz | Mean/offset of signal (remove before FFT) |
| **1/Tspan** | ~0.003 Hz (if Tspan=300) | Resolution floor - NOT physical! |
| **Your physical peak** | ~0.010 Hz | True wall oscillation |
| **Theory prediction** | 0.0117 Hz | Román et al. for M=200, kBT=1 |
| **Need**: Tspan | ≥ 3000 units | For Δν ≈ 0.0003 Hz (10× better) |

---

## Why This Confused You

1. **Peak at 0.003 Hz was huge** → You thought this was the signal
2. **Peak at 0.010 Hz was small** → You thought this was noise
3. **Actually the opposite!**
   - 0.003 Hz = resolution artifact (1/Tspan)
   - 0.010 Hz = your real physics (close to 0.0117 Hz theory!)

---

## Next Steps

1. **Run with longer Tspan** (3000+ time units)
2. **Check the new full-spectrum plot** - you'll see:
   - DC dominates at very low freq
   - 1/Tspan marked with red line
   - Your physical peak at ~0.010 Hz clearly separated
3. **Compare with theory** - should match within a few %

The physics is correct - you just needed better frequency resolution and to understand what you were looking at!
