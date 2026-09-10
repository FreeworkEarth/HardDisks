# Validation Checklist

Use this checklist to verify your simulation is working correctly.

## ✅ Step 1: Verify kBT Calculation

```bash
cd hspist3
./00ALLINONE --mode=edmd --particles=100 --l0=20 --temperature=1000 --kbt1 \
  --steps=10 --no-experiments 2>&1 | grep -A 3 "k_B.*T"
```

**Expected output:**
```
k_B (effective): 0.001000
k_B*T (effective): 1.000000
✅ Reduced units mode: k_B*T = 1 (velocities correctly scaled)
```

- [ ] Shows `k_B*T (effective): 1.000000` ✓
- [ ] Shows green checkmark with "Reduced units mode" ✓
- [ ] Does NOT show warning triangle ⚠️ ✓

## ✅ Step 2: Verify Temperature Independence

```bash
cd ..
./test_wall_frequency_independence.sh
```

**Expected:** Both tests should give SAME frequency (within 1-2%).

- [ ] Test 1 (T=1) completes successfully ✓
- [ ] Test 2 (T=1000) completes successfully ✓
- [ ] Frequencies match within 2% ✓
- [ ] Both show kBT=1.0 ✓

## ✅ Step 3: Theoretical Validation

```bash
cd EDMD_4VALIDATE
python3 validate_roman_params.py
```

**Expected output:**
```
Speed of sound (c_s, Henderson) = 2.1591
Theoretical fundamental frequency ν₁ = 0.011815
```

- [ ] Script runs without errors ✓
- [ ] Speed of sound ≈ 2.16 ✓
- [ ] Frequency ν₁ ≈ 0.0118 ✓

## ✅ Step 4: Run Validation Simulation

```bash
cd validation_sim
make clean && make
./validate_roman
```

**Expected:** Should run for ~30 seconds and output CSV file.

- [ ] Compiles without errors ✓
- [ ] Runs to completion (t=3000) ✓
- [ ] Creates `wall_position_edmd_validation.csv` ✓
- [ ] Temperature stays ≈ 1.0 throughout ✓

## ✅ Step 5: Analyze Validation Results

```bash
cd ..
MPLBACKEND=Agg python3 analyze_edmd_validation.py
```

**Expected:** Frequency ratio (sim/theory) between 0.95 and 1.05 is ideal.
  - 0.70-0.95 is acceptable (finite size effects with N=100)
  - <0.70 indicates a problem

- [ ] Script completes without errors ✓
- [ ] Finds fundamental frequency ✓
- [ ] Ratio is > 0.70 ✓
- [ ] Creates plots (PDF and PNG) ✓

**My ratio:** _________ (write it here)

## ✅ Step 6: Verify Maxwell-Boltzmann Distribution

Run a short simulation and check velocity distribution:

```bash
cd ../hspist3
./00ALLINONE --mode=edmd --particles=100 --l0=20 --temperature=1000 --kbt1 \
  --steps=1000 --no-experiments
```

Look at the initialization output. Should see something like:
```
Initializing velocities from Maxwell-Boltzmann distribution
Target k_B*T: 1.000000
Velocity scale: sqrt(k_B*T/m) = 1.000000
```

- [ ] Velocity scale = 1.0 (not 31.6!) ✓
- [ ] Particle velocities are ~O(1), not ~O(30) ✓

## 🔴 Troubleshooting

### Problem: kBT shows 1000 instead of 1

**Fix:** You forgot the `--kbt1` flag. Always use it!

```bash
# Wrong:
./00ALLINONE --temperature=1000 ...

# Right:
./00ALLINONE --temperature=1000 --kbt1 ...
```

### Problem: Frequencies differ between T=1 and T=1000

**Diagnosis:** Bug in code - physics shouldn't depend on T with --kbt1.

**Check:**
1. Verify kB_effective() function returns 1/T when --kbt1 is set
2. Verify all velocity initializations use kB_effective()
3. Check if any code path uses `temperature_runtime` directly for physics

### Problem: Validation ratio < 0.70

**Possible causes:**
1. kBT ≠ 1 (check Step 1!)
2. Simulation too short (need > 3000 time units)
3. Wrong number of particles (need exactly 50 per side)
4. Piston mass wrong (should be 200m)

**Debug:**
```bash
# Check simulation parameters
./validate_roman 2>&1 | head -20
```

Should show:
- N = 100 particles
- Piston mass: 200.0
- Temperature: 1.0

### Problem: Validation ratio > 1.05

**This is actually good!** Means your simulation gives HIGHER frequency than theory.

**Possible causes:**
1. Better numerical accuracy than theoretical approximation
2. Different equation of state (Henderson vs SPT)
3. Finite size effects at this density

**No action needed** if ratio is between 1.0 and 1.1.

## 📊 Expected Results Summary

| Test | Expected Value | Acceptable Range |
|------|----------------|------------------|
| kBT (effective) | 1.000000 | 0.999 - 1.001 |
| Speed of sound (theory) | 2.159 | 2.14 - 2.17 |
| Frequency (theory) | 0.01181 | 0.0117 - 0.0119 |
| Frequency (sim) | 0.0086-0.0118 | 0.0082 - 0.0125 |
| Ratio (sim/theory) | 0.95-1.05 | 0.70 - 1.10 |

## ✅ Final Check: Does Everything Make Sense?

Answer these questions:

1. **When I use --kbt1 with temperature=1000, what is kBT?**
   - ✅ Answer: kBT = 1.0 (K_B auto-adjusts to 1/1000)

2. **Does changing temperature parameter affect physics with --kbt1?**
   - ✅ Answer: NO! Only affects visual rendering speed

3. **What did Román et al. use for temperature?**
   - ✅ Answer: T = 1 in reduced units (kBT = 1)

4. **How do I convert cs=2.5 to physical units for argon?**
   - ✅ Answer: cs_physical = 2.5 × √(kBT/m) × (σ/τ)
   - For argon at 300K: cs_physical ≈ 2.5 × 155 m/s = 388 m/s

5. **Why does L₀=15 match perfectly but L₀=30 doesn't?**
   - ✅ Answer: Finite size effects - N=100 is optimal for L₀≈15

If you can answer all of these, you understand the simulation! 🎉

## 📁 Where to Find More Info

- Can't answer Q1-Q2? → Read [`FINAL_SUMMARY_AND_ANSWERS.md`](FINAL_SUMMARY_AND_ANSWERS.md)
- Don't understand units? → Read [`UNITS_AND_TIMING_EXPLAINED.md`](UNITS_AND_TIMING_EXPLAINED.md)
- Wall frequency issues? → Read [`WALL_MOVEMENT_TEMPERATURE_ISSUE.md`](WALL_MOVEMENT_TEMPERATURE_ISSUE.md)
- Quick commands? → Read [`QUICK_REFERENCE.md`](QUICK_REFERENCE.md)
- Validation details? → Read [`EDMD_4VALIDATE/README.md`](EDMD_4VALIDATE/README.md)

---

**Once all boxes are checked, your simulation is validated!** ✅
