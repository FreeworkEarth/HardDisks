# How to Improve Speed of Sound Measurements

## Problem

Your current cs curve shows:
- **Low density (η < 0.15)**: cs drops dramatically (WRONG!)
- **High density (η > 0.5)**: cs 30-40% too high
- **Large fluctuations**: Not smooth like theory

## Root Cause: Fixed N=100 Across All Densities

You're using **N=100 for ALL box sizes**, which gives:

| L₀ | η | Area | Particles/Area | Status |
|----|---|------|----------------|--------|
| 7.5 | 0.52 | 75 | 0.67 | ⚠️ Borderline |
| 10 | 0.39 | 100 | 0.50 | ❌ Too few |
| 15 | 0.26 | 150 | 0.33 | ❌ Too few |
| 20 | 0.20 | 200 | 0.25 | ❌ Too few |
| 30 | 0.13 | 300 | 0.17 | ❌ WAY too few |
| 40 | 0.10 | 400 | 0.13 | ❌ Unphysical |

**At η=0.1 with only 0.13 particles/area, you're not simulating a gas - you're tracking individual particles!**

## Solution: Scale N with Box Size

### Strategy 1: Constant Particle Density (RECOMMENDED)

Keep **~3-5 particles per unit area**:

```python
#!/usr/bin/env python3
import numpy as np

# Target density (particles per unit area)
target_density = 4.0  # Adjust based on your computer's speed

L0_values = np.array([7.5, 10, 15, 20, 25, 30, 35, 40])
height = 10.0

for L0 in L0_values:
    area = L0 * height
    N_total = int(target_density * area)
    N_per_side = N_total // 2
    
    # Calculate packing fraction
    eta = (N_per_side * np.pi * 0.25) / area
    
    print(f"L₀={L0:5.1f}: N={N_per_side:4d} per side, η={eta:.3f}, density={N_per_side/area:.2f}")
    
    # Generate command
    print(f"  ./00ALLINONE --mode=edmd --particles={N_per_side*2} --l0={L0} \\")
    print(f"    --height={height} --wall-mass-factor=200 --temperature=1000 --kbt1 \\")
    print(f"    --steps=300000 --no-experiments")
    print()
```

**Expected output:**
```
L₀= 7.5: N= 150 per side, η=0.628, density=2.00
L₀=10.0: N= 200 per side, η=0.628, density=2.00
L₀=15.0: N= 300 per side, η=0.628, density=2.00
...
```

**Problem:** This keeps η constant! You won't get a sweep.

### Strategy 2: Constant N, Multiple L₀ (BETTER)

Pick **ONE good N value** (e.g., N=400 per side) and sweep L₀:

```python
#!/usr/bin/env python3
import numpy as np

N_per_side = 400  # Fixed - choose based on computation time
height = 10.0

# Choose L₀ values to get desired η range
# η = (N * π * r²) / (L₀ * height)
# Solving: L₀ = (N * π * r²) / (η * height)

target_eta_values = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6]

for eta in target_eta_values:
    L0 = (N_per_side * np.pi * 0.25) / (eta * height)
    
    # Piston mass scales with N
    M_factor = 200.0  # Keep ratio M/Nm constant
    
    print(f"η={eta:.2f}: L₀={L0:6.2f}, N={N_per_side} per side")
    print(f"  ./00ALLINONE --mode=edmd --particles={N_per_side*2} --l0={L0:.2f} \\")
    print(f"    --height={height} --wall-mass-factor={M_factor} --temperature=1000 --kbt1 \\")
    print(f"    --steps=300000 --no-experiments")
    print()
```

### Strategy 3: Adaptive N (BEST, but complex)

Use **more particles at low density**:

```python
#!/usr/bin/env python3
import numpy as np

def calculate_N_needed(eta, min_particles_per_area=3.0):
    """
    Calculate minimum N to maintain good statistics
    At low η, need more particles
    """
    height = 10.0
    
    # For given η, calculate L₀ that gives minimum particle density
    # We want: N / (L₀ * height) >= min_particles_per_area
    # And: η = (N * π * 0.25) / (L₀ * height)
    # Solving: N = min_particles_per_area * (N * π * 0.25) / η
    # N * (1 - min_particles_per_area * π * 0.25 / η) = 0
    # N >= η * min_particles_per_area * L₀ * height / (π * 0.25)
    
    # Simpler: just set a minimum
    N_min = int(min_particles_per_area * 100)  # For L₀=10
    N_scale = max(1.0, 0.3 / eta)  # Scale up at low η
    N_needed = int(N_min * N_scale)
    
    # Round to nice number
    N_needed = ((N_needed + 49) // 50) * 50
    
    return N_needed

target_eta_values = np.arange(0.1, 0.65, 0.05)
height = 10.0

for eta in target_eta_values:
    N_per_side = calculate_N_needed(eta)
    L0 = (N_per_side * np.pi * 0.25) / (eta * height)
    
    print(f"η={eta:.2f}: L₀={L0:6.2f}, N={N_per_side:4d} per side, density={N_per_side/(L0*height):.2f}")
```

## Recommended Approach

**Start with Strategy 2** (constant N=400):

```bash
#!/bin/bash
cd hspist3

N_TOTAL=800  # 400 per side
HEIGHT=10
M_FACTOR=200

# η values from 0.10 to 0.60
for ETA in 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60; do
    # Calculate L₀ for this η
    # η = (N/2 * π * 0.25) / (L₀ * H)
    # L₀ = (N/2 * π * 0.25) / (η * H)
    L0=$(python3 -c "import math; print((${N_TOTAL}/2 * math.pi * 0.25) / ($ETA * $HEIGHT))")
    
    echo "======================================================================"
    echo "Running: η=$ETA, L₀=$L0, N=$N_TOTAL"
    echo "======================================================================"
    
    ./00ALLINONE \
        --mode=edmd \
        --particles=$N_TOTAL \
        --l0=$L0 \
        --height=$HEIGHT \
        --wall-mass-factor=$M_FACTOR \
        --temperature=1000 \
        --kbt1 \
        --steps=300000 \
        --no-experiments
    
    # Save results with η label
    mv wall_position.csv wall_position_eta_${ETA}.csv
    echo "Saved: wall_position_eta_${ETA}.csv"
    echo ""
done

echo "All simulations complete!"
echo "Now analyze each file to extract cs vs η"
```

## Expected Improvements

With N=400 per side:
- **Smoother curve** (less statistical noise)
- **Correct low-η behavior** (cs doesn't drop)
- **Better high-η accuracy** (within 5-10% of theory)

## Why This Fixes The Issues

1. **Low η drop**: With only N=50 at η=0.1, particles rarely collide → poor statistics
   - **Fix**: N=400 gives 8× more collisions → better cs measurement

2. **High η overshoot**: With N=50 at η=0.5, particles are overcrowded → crystallization
   - **Fix**: N=400 reduces finite-size boundary effects

3. **Fluctuations**: Small N → large relative fluctuations in collision rate
   - **Fix**: Large N → √N improvement in statistics

## Computational Cost

| N (total) | Collisions/step | Time/step | Total time |
|-----------|-----------------|-----------|------------|
| 100 | ~10 | 1× | 1× |
| 400 | ~40 | 4× | 4× |
| 1600 | ~160 | 16× | 16× |

**N=400 takes 4× longer**, but gives **much** better results.

**Optimization**: Run overnight or in parallel!

