# Packing Fraction (η) Guide

## What is the "REAL" Packing Fraction?

### Standard Literature Definition (Román, Henderson, SPT, etc.)

```
η = (N × π × r²) / A_geometric
```

Where:
- **N** = total number of particles
- **r** = particle radius
- **A_geometric** = geometric box area (2L₀ × H)

**Important:** This does **NOT** account for:
- Wall volume ❌
- Particle exclusion zones near walls ❌
- Any container-related corrections ❌

### Why This Definition?

η should be an **intrinsic property of the fluid**, independent of the container.

For **infinite systems** (no walls):
- η represents the fraction of space occupied by disks
- Theory (Henderson, SPT) is developed for infinite systems
- Using geometric area makes η container-independent

For **confined systems** (with walls):
- Actual local density may be slightly higher due to wall exclusion
- But we still use geometric definition for consistency with theory

---

## Three Definitions Compared

The Python script `wall_x_FFT.py` now calculates and displays all three:

### 1. η_nominal (Standard - Used in Literature)
```
Area = 2L₀ × H
η = (N × π × r²) / (2L₀ × H)
```
**This is what you should use for comparison to literature!**

### 2. η_exclude_wall
```
Area = 2L₀ × H - t × H    (exclude wall volume)
η = (N × π × r²) / (2L₀ × H - t × H)
```
Accounts for wall thickness t, but wall is typically thin (t ~ σ).

### 3. η_exclude_wall_and_exclusion
```
Area = 2L₀ × H - t × H - 2r × H    (exclude wall + dead zones)
η = (N × π × r²) / (2L₀ × H - t × H - 2r × H)
```
Also accounts for the fact that particle **centers** can't get closer than r to wall edges.

### Typical Values for Your System

For N=100, L₀=20, H=10, r=0.5, t=1.0:

| Definition | η | Difference |
|------------|---|------------|
| Nominal | 0.1963 | baseline |
| Exclude wall | 0.2014 | +2.6% |
| Exclude wall+exclusion | 0.2067 | +5.3% |

**The differences are small!** Wall exclusion effects are only ~5% for typical parameters.

---

## How to Maintain Constant η When Changing Particle Radius

### The Problem

You want to:
- Change particle radius r to fit more/fewer particles
- Keep η constant for comparison

### The Solution

From the definition:
```
η = (N × π × r²) / (2L₀ × H)
```

Rearranging:
```
N = η × (2L₀ × H) / (π × r²)
```

### Example: Keep η = 0.196 for Different Radii

**Current setup:**
- L₀ = 20 σ
- H = 10 σ
- η = 0.196

**Calculate N for different r:**

| r (σ) | N (particles) | Command |
|-------|---------------|---------|
| 0.3 | 277 | `--eta=0.196 --l0=20 --height=10` |
| 0.4 | 156 | `--eta=0.196 --l0=20 --height=10` |
| **0.5** | **100** | `--eta=0.196 --l0=20 --height=10` ← Current |
| 0.6 | 69 | `--eta=0.196 --l0=20 --height=10` |
| 0.7 | 51 | `--eta=0.196 --l0=20 --height=10` |

### Command Examples

**Use the `--eta` flag** (automatically calculates radius):

```bash
# η = 0.196 with different particle counts
./00ALLINONE --mode=edmd --particles=100 --eta=0.196 \
  --l0=20 --height=10 --wall-mass-factor=200 --kbt1 --auto-release

./00ALLINONE --mode=edmd --particles=200 --eta=0.196 \
  --l0=20 --height=10 --wall-mass-factor=200 --kbt1 --auto-release
```

**The code automatically calculates the correct radius!**

---

## Román et al. (2002) Parameters

From Table I of the paper:

| L₀ (σ) | N (per side) | η | cs (measured) | cs (theory) |
|--------|--------------|---|---------------|-------------|
| 7.5 | 50 | 0.524 | 3.69 ± 0.03 | 3.63 |
| 10.0 | 50 | 0.393 | 3.04 ± 0.02 | 3.01 |
| 15.0 | 50 | 0.262 | 2.48 ± 0.02 | 2.47 |
| **20.0** | **50** | **0.196** | **2.20 ± 0.02** | **2.16** |
| 25.0 | 50 | 0.157 | 2.04 ± 0.02 | 2.03 |

All use:
- Particle radius: r = 0.5 σ (diameter = 1 σ)
- Height: H = 10 σ
- η calculated using **geometric area** (not excluding wall)

---

## Key Takeaways

1. **Use η_nominal (geometric area) for theory comparison** ✅
   - This is standard in all literature
   - Makes η container-independent
   - Matches Román, Henderson, SPT conventions

2. **Wall exclusion effects are small (~5%)** ✅
   - For thin walls (t ~ σ) in large boxes (L₀ >> σ)
   - Only matters at very high η or very small L₀

3. **To vary N while keeping η constant:** ✅
   - Use `--eta=0.196` flag
   - Code automatically adjusts particle radius
   - OR manually calculate: N = η × (2L₀ × H) / (π × r²)

4. **The 5% deviation from Henderson a=0.125 is NOT due to η definition** ⚠️
   - Wall exclusion only explains ~2.6% (if using Definition 2)
   - The real issue is likely:
     - Short run length (not fully equilibrated)
     - Different initialization than Román (T=1 vs pre-equilibration)
     - Need longer simulation with T=100000 initialization

---

## Recommendations

### For Validation Against Román et al.

```bash
./00ALLINONE --mode=edmd --particles=100 --l0=20 --height=10 \
  --wall-mass-factor=200 --temperature=100000 --kbt1 \
  --auto-release --steps=5000000 --single-test

python3 wall_x_FFT.py
```

This will show:
- All three η definitions
- Predicted cs and ν for each definition
- Comparison to your measured values
- Clear verdict on which definition matches best

### For Future Experiments (Varying η)

```bash
# Low density
./00ALLINONE --mode=edmd --particles=100 --eta=0.10 \
  --l0=20 --height=10 --wall-mass-factor=200 --kbt1 --auto-release

# Medium density (Román's value)
./00ALLINONE --mode=edmd --particles=100 --eta=0.196 \
  --l0=20 --height=10 --wall-mass-factor=200 --kbt1 --auto-release

# High density
./00ALLINONE --mode=edmd --particles=100 --eta=0.40 \
  --l0=20 --height=10 --wall-mass-factor=200 --kbt1 --auto-release
```

The `--eta` flag handles the radius calculation automatically!

---

**Created:** 2025-12-10
**Last Updated:** 2025-12-10
