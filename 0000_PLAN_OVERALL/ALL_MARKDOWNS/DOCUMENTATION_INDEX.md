# Documentation Index

Complete guide to all documentation for the Hard Disk Piston Simulation project.

## 🚀 Start Here

| Document | Purpose | Read This If... |
|----------|---------|-----------------|
| **[QUICK_REFERENCE.md](QUICK_REFERENCE.md)** | One-page command reference | You need a quick reminder of commands |
| **[VALIDATION_CHECKLIST.md](VALIDATION_CHECKLIST.md)** | Step-by-step verification | You want to verify everything works |
| **[README.md](README.md)** | Project overview | You're new to this project |

## 📚 Understanding Reduced Units

| Document | Purpose | Read This If... |
|----------|---------|-----------------|
| **[UNITS_AND_TIMING_EXPLAINED.md](UNITS_AND_TIMING_EXPLAINED.md)** | Complete units guide with "Explain Like 5" section | You're confused about σ=1, m=1, kBT=1 |
| **[FINAL_SUMMARY_AND_ANSWERS.md](FINAL_SUMMARY_AND_ANSWERS.md)** | Answers to 8 key questions | You have questions about T=1000 vs T=1K |
| **[WALL_MOVEMENT_TEMPERATURE_ISSUE.md](WALL_MOVEMENT_TEMPERATURE_ISSUE.md)** | Why wall frequency depends on kBT | You're confused about temperature dependence |
| **[TEMPERATURE_FIX_GUIDE.md](TEMPERATURE_FIX_GUIDE.md)** | Original temperature debugging | Historical context (may be outdated) |

## 🔬 Validation Suite

| Document | Purpose | Read This If... |
|----------|---------|-----------------|
| **[EDMD_4VALIDATE/README.md](EDMD_4VALIDATE/README.md)** | Validation suite documentation | You want to validate against Román (2002) |
| **[EDMD_4VALIDATE/validate_roman_params.py](EDMD_4VALIDATE/validate_roman_params.py)** | Theoretical predictions script | You need theory values for comparison |
| **[EDMD_4VALIDATE/analyze_edmd_validation.py](EDMD_4VALIDATE/analyze_edmd_validation.py)** | FFT analysis script | You have wall data and need frequency |

## 🧪 Test Scripts

| Script | Purpose | Run This If... |
|--------|---------|----------------|
| **[test_wall_frequency_independence.sh](test_wall_frequency_independence.sh)** | Tests T=1 vs T=1000 give same physics | You want to verify --kbt1 works correctly |

## 🎯 By Use Case

### "I'm just starting and don't understand anything"

1. Read: [UNITS_AND_TIMING_EXPLAINED.md](UNITS_AND_TIMING_EXPLAINED.md) (especially "Explain Like 5" section)
2. Read: [FINAL_SUMMARY_AND_ANSWERS.md](FINAL_SUMMARY_AND_ANSWERS.md)
3. Follow: [VALIDATION_CHECKLIST.md](VALIDATION_CHECKLIST.md)

### "I want to run simulations the RIGHT way"

1. Read: [QUICK_REFERENCE.md](QUICK_REFERENCE.md) → Section "Running Simulations Correctly"
2. **Key command:**
```bash
./hspist3/00ALLINONE --mode=edmd --particles=100 --l0=20 --height=10 \
  --wall-mass-factor=200 --temperature=1000 --kbt1 --steps=300000
```

### "My results don't match Román et al. (2002)"

1. Run: [VALIDATION_CHECKLIST.md](VALIDATION_CHECKLIST.md) (all steps)
2. Read: [EDMD_4VALIDATE/README.md](EDMD_4VALIDATE/README.md) → Troubleshooting section
3. Check: Did you use `--kbt1` flag? (Most common issue!)

### "I don't understand why T=1000 doesn't mean 1000 Kelvin"

1. Read: [FINAL_SUMMARY_AND_ANSWERS.md](FINAL_SUMMARY_AND_ANSWERS.md) → Question 1
2. Read: [UNITS_AND_TIMING_EXPLAINED.md](UNITS_AND_TIMING_EXPLAINED.md) → "Explain Like You're 5"
3. **TL;DR:** With `--kbt1`, temperature is dimensionless and just controls visual speed

### "Wall frequency changes when I change temperature"

1. Read: [WALL_MOVEMENT_TEMPERATURE_ISSUE.md](WALL_MOVEMENT_TEMPERATURE_ISSUE.md)
2. Run: `./test_wall_frequency_independence.sh`
3. **Diagnosis:** You probably forgot `--kbt1` flag

### "I want to convert results to real units (e.g., argon gas)"

1. Read: [UNITS_AND_TIMING_EXPLAINED.md](UNITS_AND_TIMING_EXPLAINED.md) → "Converting to Physical Units"
2. Read: [FINAL_SUMMARY_AND_ANSWERS.md](FINAL_SUMMARY_AND_ANSWERS.md) → Question 5
3. **Formula:** cs_physical = cs_reduced × √(kBT/m) × (σ/τ)

### "I want to reproduce Román et al. Figure 3"

1. Read: [EDMD_4VALIDATE/README.md](EDMD_4VALIDATE/README.md)
2. Run validation suite:
```bash
cd EDMD_4VALIDATE/validation_sim
make && ./validate_roman
cd .. && python3 analyze_edmd_validation.py
```

## 📖 Document Descriptions

### Core Documentation

#### [QUICK_REFERENCE.md](QUICK_REFERENCE.md)
- One-page cheat sheet
- Common commands
- Keyboard controls
- Warning signs (physics broken vs correct)
- **When to read:** Daily reference

#### [VALIDATION_CHECKLIST.md](VALIDATION_CHECKLIST.md)
- Step-by-step verification procedure
- Expected outputs for each step
- Troubleshooting guide
- **When to read:** After making changes, before publishing results

#### [README.md](README.md)
- Project overview
- Directory structure
- Build instructions
- **When to read:** First time using the project

### Theory & Understanding

#### [FINAL_SUMMARY_AND_ANSWERS.md](FINAL_SUMMARY_AND_ANSWERS.md)
**The "bible" for this project.** Answers 8 critical questions:

1. Is T=1000 means 1000K? (NO!)
2. Can I use physical units? (YES but impractical)
3. Does MB distribution work? (YES, universal)
4. Particles 31.6× faster but physics unchanged? (YES with --kbt1)
5. Just multiply by conversion factor? (YES)
6. What T did Román use? (T=1 reduced)
7. Wall depends on T? (On kBT, but --kbt1 fixes kBT=1)
8. Rename to reduced units? (YES, clearer)

**When to read:** When you have ANY question about units or temperature

#### [UNITS_AND_TIMING_EXPLAINED.md](UNITS_AND_TIMING_EXPLAINED.md)
**Most comprehensive units guide.**

Contains:
- "Explain Like You're 5" section (NEW!)
- Reduced units definition
- Conversion formulas
- Example calculations
- **When to read:** First time learning reduced units

#### [WALL_MOVEMENT_TEMPERATURE_ISSUE.md](WALL_MOVEMENT_TEMPERATURE_ISSUE.md)
**Explains why wall oscillation depends on kBT.**

Key insight: cs ∝ √(kBT), so wall frequency ∝ √(kBT)

With `--kbt1`, kBT is ALWAYS 1 regardless of T parameter!

**When to read:** Confused about temperature dependence

#### [TEMPERATURE_FIX_GUIDE.md](TEMPERATURE_FIX_GUIDE.md)
Original debugging document from when temperature system was broken.

**Status:** Historical reference (superseded by FINAL_SUMMARY)

**When to read:** Understanding project history

### Validation

#### [EDMD_4VALIDATE/README.md](EDMD_4VALIDATE/README.md)
Complete guide to validation suite.

Contains:
- Quick start instructions
- File descriptions
- Expected results table
- Troubleshooting guide
- **When to read:** Before running validation

#### [EDMD_4VALIDATE/validate_roman_params.py](EDMD_4VALIDATE/validate_roman_params.py)
Python script that calculates theoretical predictions.

**Output:**
- Speed of sound from Henderson equation
- Fundamental frequency from transcendental equation
- Period

**When to run:** Before validation to know what to expect

#### [EDMD_4VALIDATE/analyze_edmd_validation.py](EDMD_4VALIDATE/analyze_edmd_validation.py)
Python script that analyzes simulation output.

**Input:** `wall_position_edmd_validation.csv`

**Output:**
- FFT frequency spectrum
- Fundamental frequency
- Comparison with theory
- Plots (PDF/PNG)

**When to run:** After validation simulation completes

### Test Scripts

#### [test_wall_frequency_independence.sh](test_wall_frequency_independence.sh)
Bash script that proves physics is independent of T parameter (with --kbt1).

Runs TWO simulations:
1. T=1 with --kbt1
2. T=1000 with --kbt1

Then compares wall frequencies.

**Expected:** Frequencies should match within 1-2%

**When to run:** To verify --kbt1 flag works correctly

## 🔑 Key Concepts Explained

### What is "Reduced Units"?

Set σ=1, m=1, kBT=1 to make equations dimensionless.

**Why?**
- Cleaner equations (no tiny constants like 10⁻²³)
- Better numerical precision
- Same code works for any system (just convert at end)

**See:** [UNITS_AND_TIMING_EXPLAINED.md](UNITS_AND_TIMING_EXPLAINED.md)

### What does --kbt1 flag do?

Forces K_B = 1/T, so kBT always equals 1.

**Example:**
- `--temperature=1 --kbt1` → K_B=1, kBT=1
- `--temperature=1000 --kbt1` → K_B=0.001, kBT=1

**Result:** Temperature parameter becomes "visual speed dial" only

**See:** [FINAL_SUMMARY_AND_ANSWERS.md](FINAL_SUMMARY_AND_ANSWERS.md) Question 4

### Why does L₀=15 match perfectly but others don't?

**Finite size effects** with N=100 particles.

- Small L₀ (high density): Particle caging → cs too high
- Large L₀ (low density): Too few collisions → cs too low
- L₀=15 (η≈0.26): "Sweet spot"

**Solution:** Use more particles (N=256+)

**See:** [WALL_MOVEMENT_TEMPERATURE_ISSUE.md](WALL_MOVEMENT_TEMPERATURE_ISSUE.md)

## 📊 Results Interpretation

### Speed of Sound

| L₀ | η | Your cs | Theory cs | Status |
|----|---|---------|-----------|--------|
| 15 | 0.262 | 2.54 | 2.53 | ✅ Perfect! |
| 20 | 0.196 | ? | 2.16 | Run validation |

**Goal:** Ratio (your/theory) between 0.95 and 1.05

**See:** [EDMD_4VALIDATE/README.md](EDMD_4VALIDATE/README.md) → Understanding Results

### Frequency

**Theory (Román Fig 3, L₀=20):** ν₁ ≈ 0.0118

**Acceptable range:** 0.0082 - 0.0125 (0.70 - 1.05 ratio)

**If too low (<0.70):** Check kBT value, particle count, piston mass

**See:** [VALIDATION_CHECKLIST.md](VALIDATION_CHECKLIST.md) → Step 5

## 🛠️ Common Commands

See [QUICK_REFERENCE.md](QUICK_REFERENCE.md) for full list.

**Check kBT:**
```bash
./hspist3/00ALLINONE --mode=edmd --l0=20 --temperature=1000 --kbt1 \
  --steps=10 --no-experiments 2>&1 | grep "k_B.*T"
```

**Run validation:**
```bash
cd EDMD_4VALIDATE
python3 validate_roman_params.py
cd validation_sim && make && ./validate_roman
cd .. && python3 analyze_edmd_validation.py
```

**Test independence:**
```bash
./test_wall_frequency_independence.sh
```

## 📅 Reading Order (Recommended)

### For Beginners:
1. [QUICK_REFERENCE.md](QUICK_REFERENCE.md) (5 min)
2. [UNITS_AND_TIMING_EXPLAINED.md](UNITS_AND_TIMING_EXPLAINED.md) → "Explain Like 5" (10 min)
3. [FINAL_SUMMARY_AND_ANSWERS.md](FINAL_SUMMARY_AND_ANSWERS.md) (15 min)
4. [VALIDATION_CHECKLIST.md](VALIDATION_CHECKLIST.md) (follow steps)

### For Researchers:
1. [README.md](README.md) (5 min)
2. [EDMD_4VALIDATE/README.md](EDMD_4VALIDATE/README.md) (10 min)
3. Run validation suite
4. [FINAL_SUMMARY_AND_ANSWERS.md](FINAL_SUMMARY_AND_ANSWERS.md) (reference)

### For Troubleshooting:
1. [VALIDATION_CHECKLIST.md](VALIDATION_CHECKLIST.md) → Troubleshooting
2. [WALL_MOVEMENT_TEMPERATURE_ISSUE.md](WALL_MOVEMENT_TEMPERATURE_ISSUE.md)
3. [QUICK_REFERENCE.md](QUICK_REFERENCE.md) → Warning Signs

## 📝 Document Statistics

| Document | Size | Type | Priority |
|----------|------|------|----------|
| FINAL_SUMMARY_AND_ANSWERS.md | 6.4K | Theory | 🔴 High |
| UNITS_AND_TIMING_EXPLAINED.md | 13K | Theory | 🔴 High |
| QUICK_REFERENCE.md | 5.7K | Reference | 🔴 High |
| VALIDATION_CHECKLIST.md | 7.5K | Procedure | 🟠 Medium |
| WALL_MOVEMENT_TEMPERATURE_ISSUE.md | 6.6K | Theory | 🟠 Medium |
| EDMD_4VALIDATE/README.md | 6.2K | Validation | 🟠 Medium |
| TEMPERATURE_FIX_GUIDE.md | 5.6K | Historical | 🟢 Low |
| README.md | 2.6K | Overview | 🟢 Low |

---

**Last updated:** November 2024

**Questions?** Start with [QUICK_REFERENCE.md](QUICK_REFERENCE.md) or [FINAL_SUMMARY_AND_ANSWERS.md](FINAL_SUMMARY_AND_ANSWERS.md)
