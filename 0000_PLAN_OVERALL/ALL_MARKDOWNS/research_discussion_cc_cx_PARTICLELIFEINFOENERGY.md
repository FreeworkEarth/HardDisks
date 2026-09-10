# CC → CX Research Discussion: PARTICLELIFE INFO ENERGY

**Purpose**: Collaborative research discussion between Claude Code (CC) and Codex (CX)
**Format**: Structured conversation log with automatic responses

---

# Discussion Log

---

## Round 1: Initial PARTICLELIFE Implementation & Testing

**Initiated by**: User & CC
**Date**: 2025-10-31

### User Question:

Implement PARTICLELIFE experiment in hspist3/00ALLINONE.c to visualize correlation matrices as artificial life simulations.

**Requirements**:
1. Different particle types with different colors
2. Interaction matrix defining attraction/repulsion between types (like correlation matrix, values -1 to +1)
3. Cyclic/periodic boundary conditions (wrap-around at edges)
4. Force-based physics (not hard-disk collisions)
5. CLI interface for MODE_TIME simulation
6. Enable different radii per particle type (future)

### CC Response:

**Date**: 2025-10-31

I've successfully implemented the PARTICLELIFE experiment in [00ALLINONE.c](hspist3/00ALLINONE.c). Here's what was built:

#### 1. **Architecture Overview**

The implementation adds a complete particle-life physics engine alongside the existing hard-disk simulation, cleanly separated with `##CHRIS` markers:

- **Data structures** ([00ALLINONE.c:542-576](hspist3/00ALLINONE.c#L542-L576)): Interaction matrix, type arrays, colors
- **Physics engine** ([00ALLINONE.c:3133-3298](hspist3/00ALLINONE.c#L3133-L3298)): Force calculation, cyclic boundaries
- **Integration** ([00ALLINONE.c:5513-5519](hspist3/00ALLINONE.c#L5513-L5519)): Replaces hard-disk physics in MODE_TIME when enabled

#### 2. **Core Physics Implementation**

**Force Formula** (from particle-life-app analysis):
```c
if (dist < β=0.3):
    force = (dist/β) - 1          // Close-range repulsion
else:
    force = a * (1 - |1+β-2*dist|/(1-β))  // Long-range interaction (matrix value 'a')
```

**Cyclic Boundaries**:
- Toroidal topology (wrap-around at edges)
- Minimum image convention for distance calculations
- Handles particles crossing boundaries correctly

**Integration**:
```c
// Apply friction
Vx[i] *= friction;  // Default 0.85
Vy[i] *= friction;

// Update velocity with acceleration
Vx[i] += ax[i] * dt;
Vy[i] += ay[i] * dt;

// Update position
X[i] += Vx[i] * dt;
Y[i] += Vy[i] * dt;

// Wrap around boundaries
pl_apply_cyclic_boundaries(&X[i], &Y[i]);
```

#### 3. **Default Configuration**

**3 particle types** with "rock-paper-scissors" interaction pattern:
```
Interaction Matrix:
   -0.10   0.50  -0.30   (Red: repels self, attracts Green, repels Blue)
   -0.30  -0.10   0.50   (Green: repels Red, repels self, attracts Blue)
    0.50  -0.30  -0.10   (Blue: attracts Red, repels Green, repels self)
```

**Colors**: Red (255,0,0), Green (0,255,0), Blue (0,0,255)

**Parameters**:
- `rmax = 0.02`: Interaction radius
- `friction = 0.85`: Velocity damping
- `β = 0.3`: Close-range repulsion threshold

#### 4. **Usage**

```bash
cd hspist3
./00ALLINONE --experiment=particlelife
```

Automatically:
- Disables walls/dividers
- Enables cyclic boundaries
- Uses MODE_TIME integration (default)
- Renders particles with type-based colors

#### 5. **Key Design Decisions**

**✅ Why MODE_TIME only**:
- EDMD requires exact collision times (incompatible with smooth forces)
- Time-stepping integration perfect for force-based physics
- User requested MODE_TIME specifically

**✅ Cyclic boundaries implementation**:
- Normalized coordinates [0,1] internally
- Converts to/from sigma units
- Minimum image convention for pair distances

**✅ O(N²) pairwise forces**:
- Simple and correct for initial implementation
- Cutoff at `rmax` for efficiency
- Spatial partitioning easy to add later

**✅ Type assignment**:
- Random uniform distribution
- Equal probability for all types
- Easy to extend with CLI parameters

#### 6. **Extensibility Points**

All designed for easy extension:

1. **Custom matrices**: Add `--pl-matrix=file.txt` CLI parameter
2. **More types**: Add `--pl-num-types=N` parameter
3. **Different radii**: `pl_particle_radii[]` array already exists
4. **Position patterns**: Add type-based spatial initialization
5. **Realistic forces**: Modify `pl_compute_force()` for LJ potential, Coulomb, etc.

#### 7. **Current Bug Report**

User reports: **Particles clump in top-left corner when released from grid**

**Hypothesis**: Likely issue in `pl_apply_cyclic_boundaries()` coordinate conversion:
```c
float x_norm = (*x - (double)XW1 / PIXELS_PER_SIGMA) / width_sigma;
```

This mixes pixel coordinates with sigma units. The X, Y arrays are in **pixels**, not sigma units!

**Test Requested**: Simpler matrix (diagonal = 1, off-diagonal = 0) for debugging.

#### 8. **Next Steps**

**Immediate**:
1. Fix coordinate system bug in cyclic boundaries
2. Implement diagonal-only test matrix
3. Verify particles don't clump

**Future**:
- Add CLI for custom matrices
- Add more interaction patterns (chains, snakes, etc.)
- Optimize with spatial partitioning
- Add per-type radii support

---

### CC Response (Bug Fixes):

**Date**: 2025-10-31 (continued)

#### Bug Fix #1: Coordinate System Mismatch ✅

**Root Cause**: Mixed pixel and sigma units in cyclic boundary calculations
- `X[i]`, `Y[i]` arrays store positions in **pixels**
- Original code incorrectly divided by `PIXELS_PER_SIGMA`
- Particles wrapped to wrong positions (top-left corner clumping)

**Fix**: [00ALLINONE.c:3157-3205](hspist3/00ALLINONE.c#L3157-L3205)
```c
// Convert to normalized [0,1] using PIXELS (not sigma)
float x_norm = (float)((*x - (double)XW1) / width_pixels);
// Wrap around
while (x_norm < 0.0f) x_norm += 1.0f;
// Convert back to pixels
*x = (double)XW1 + x_norm * width_pixels;
```

#### Bug Fix #2: Force Formula Normalization ✅

**Root Cause**: Particle-life formula expects normalized distances [0,1], but we had pixel distances
- Original: `beta = 60 pixels` → denominator `(1 - 60) = -59` → broken formula
- Force magnitude calculations were incorrect

**Fix**: [00ALLINONE.c:3207-3247](hspist3/00ALLINONE.c#L3207-L3247)
```c
// Normalize distance by rmax to get [0,1] range
double dist_norm = dist / pl_rmax;
const double beta_normalized = 0.3;  // Fixed from particle-life

// Now formula works correctly
if (dist_norm < beta_normalized) {
    force_mag = (dist_norm / beta_normalized) - 1.0;
} else {
    force_mag = a * (1.0 - |1 + 0.3 - 2*dist_norm| / 0.7);
}

// Scale to pixel world for visible motion
force_mag *= pl_force_multiplier * 1000.0;
```

#### Bug Fix #3: Interaction Parameters ✅

**Changed**:
- `pl_rmax`: 0.02 → **200 pixels** (interaction radius)
- Force scale: 100 → **1000** (10x faster motion)
- Test matrix: Diagonal (types attract themselves, ignore others)

**Result**:
```
Matrix:
  1.0  0.0  0.0   (Red attracts Red only)
  0.0  1.0  0.0   (Green attracts Green only)
  0.0  0.0  1.0   (Blue attracts Blue only)
```

#### Validation: ✅ SUCCESS

**User Report**: "HELL YA THEY CLUSTER !!!"

**Observed Behavior**:
- Particles start in grid (colored by type)
- Press SPACE → particles accelerate
- Each color forms separate clusters
- Clusters move around box
- Cyclic boundaries work correctly (wrap at edges)
- No corner clumping!

---

### Synthesis (Round 1 Complete):

**Status**: ✅ Implementation complete, ✅ All bugs fixed, ✅ User validated

**Key Achievements**:
1. Full particle-life physics engine integrated
2. Cyclic boundaries working correctly
3. Force calculations normalized properly
4. Clustering behavior validated

**Next Direction**: User requests extended features:
1. **Increase force range**: Make interactions fall off more gradually with distance
2. **Rule tensors**: Beyond 2D matrix, add higher-dimensional interaction rules

**Action Items for Round 2**:
1. CC: Research rule tensor architectures for particle interactions
2. CC: Propose force range modifications
3. User: Clarify desired tensor structure

---


---

## Round 2: Rule Tensors & Enhanced Force Models

**Initiated by**: User
**Date**: 2025-10-31

### User Question:

The clustering works great, but:
1. **Increase force over distance**: Current force drops off too quickly - particles only interact within rmax
2. **Rule tensors**: Instead of just 2D matrix `M[i][j]`, can we add higher-dimensional rules?

**Context**: Same initialization produces same behavior - want more complex emergent dynamics.

---

### CC Response:

**Date**: 2025-10-31

#### 1. Force Range Extension Options

**Current Limitation**: Forces cut off sharply at `pl_rmax = 200 pixels`
- Inside rmax: Particle-life formula applies
- Outside rmax: Force = 0 (hard cutoff)

**Option A: Soft Cutoff (Simple)**
Add smooth decay near rmax boundary:

```c
if (dist > pl_rmax) {
    // Soft cutoff: exponential decay beyond rmax
    double excess = (dist - pl_rmax) / pl_rmax;
    force_mag *= exp(-excess * 5.0);  // Decay rate = 5.0
}
```

**Option B: Power Law Falloff (Realistic)**
Replace particle-life formula with physics-inspired potential:

```c
// Lennard-Jones style: attractive + repulsive
double sigma = 10.0;  // Particle diameter
double r_norm = dist / sigma;
force_mag = a * (pow(r_norm, -13) - pow(r_norm, -7));  // 12-6 potential derivative
```

**Option C: Extended Range (Hybrid)**
Keep particle-life formula but extend rmax:

```c
pl_rmax = 400.0f;  // Double the range (50% of box width)
// Or make it adaptive: rmax = 0.5 * min(SIM_WIDTH, SIM_HEIGHT)
```

**Recommendation**: Start with **Option C** (simple) + **Option A** (soft cutoff) for smooth behavior.

---

#### 2. Rule Tensor Architecture

**Current**: 2D interaction matrix `M[i][j]` (type i → type j)

**Tensor Extension Ideas**:

**Option 1: Distance-Dependent Tensors (3D)**
```c
// Tensor: M[i][j][distance_bin]
// Different interaction strength at different distances
float pl_interaction_tensor[MAX_TYPES][MAX_TYPES][N_DISTANCE_BINS];

// Example: Red-Green interaction
// - Close range (0-50 px): Repulsion (-0.5)
// - Medium range (50-100 px): Attraction (+1.0)
// - Far range (100-200 px): Weak attraction (+0.2)
```

**Option 2: Velocity-Dependent Rules (4D)**
```c
// Tensor: M[i][j][v_self][v_other]
// Interaction depends on relative velocities
float pl_interaction_tensor[MAX_TYPES][MAX_TYPES][N_VEL_BINS][N_VEL_BINS];

// Example: Predator-prey dynamics
// - Fast prey + slow predator: Prey escapes (repulsion)
// - Slow prey + fast predator: Predator catches (strong attraction)
```

**Option 3: Angular Rules (5D)**
```c
// Tensor: M[i][j][distance][angle]
// Directional interactions (like vision cones)
float pl_interaction_tensor[MAX_TYPES][MAX_TYPES][N_DIST_BINS][N_ANGLE_BINS];

// Example: Flocking behavior
// - Behind: Weak attraction (follow)
// - Ahead: Strong attraction (catch up)
// - Sideways: Repulsion (maintain spacing)
```

**Option 4: Local Density Rules (Environmental)**
```c
// Tensor: M[i][j][local_density]
// Crowding effects
float pl_interaction_tensor[MAX_TYPES][MAX_TYPES][N_DENSITY_BINS];

// Example: Clustering threshold
// - Low density: Strong attraction (gather)
// - High density: Repulsion (avoid overcrowding)
```

**Recommendation for first implementation**: **Option 1 (Distance bins)** - Most intuitive, easy to visualize, powerful for emergent behavior.

---




#### 3. User's Grand Vision: Connecting the Dots

**Date**: 2025-10-31

User's response reveals a **much larger research program**:

**Phase 1: Particle-Life Foundation** (Current)
- Study simple correlation matrices → emergent behavior
- Distance bins → velocity bins → angular bins → density bins
- Rich exploratory space for rule discovery

**Phase 2: Web Deployment**
- Compile to WebAssembly (Emscripten)
- Upload as interactive web demo
- Share with research community

**Phase 3: Information Theory Integration**
- **Mutual Information (MI)** calculations on particle system
- **Energy** measurements (thermodynamic analogy)
- Model interaction rules as **random variables** with different distributions
- Study information flow in particle interactions

**Phase 4: Connection to IB Package**
- Use particle-life as **visualization tool** for Information Bottleneck
- Multi-dimensional tensor = IB compression space
- Particle types = information channels
- Interaction rules = compression operators
- **Delete similar solutions** (redundancy elimination) → particle clustering?

**Phase 5: Neural Network Analog** (Long-term vision)
- Particle interactions → **Neuroplasticity** (synaptic weight changes)
- Distance/velocity bins → **Spiking neural networks** (temporal dynamics)
- Rule tensors → **Learning rules** (STDP, Hebbian, etc.)
- Energy minimization → **Network training**

---

### CC Response (Strategic Analysis):

**Date**: 2025-10-31

This is **absolutely feasible** and incredibly exciting! Here's why:

#### ✅ Why This Vision Works

**1. Common Mathematical Framework**
All these systems share core principles:

| Particle-Life | Information Theory | Neural Networks |
|---------------|-------------------|-----------------|
| Interaction matrix M[i][j] | Channel capacity C[i][j] | Weight matrix W[i][j] |
| Force calculation | Information transfer | Synaptic transmission |
| Particle types | Message types/symbols | Neuron populations |
| Energy conservation | Entropy/MI bounds | Energy-based learning |
| Clustering | Compression | Feature learning |

**2. Particle-Life as IB Visualizer**

**Brilliant insight**: Map IB concepts to particle dynamics!

```
Information Bottleneck          Particle-Life Analog
-------------------            --------------------
Input X                    →   Particle type i
Compressed T               →   Cluster membership
Output Y                   →   Interaction pattern
MI(T;Y)                    →   Force strength (information preserved)
H(T|X)                     →   Entropy within clusters
Compression ratio          →   Cluster tightness
Lagrangian β               →   Temperature parameter
```

**Visualization idea**: 
- Each IB iteration = particle timestep
- Information loss = particle dissipation (friction)
- Relevant features = stable clusters
- Redundant features = particles that merge

**3. Neural Network Connection**

The progression is natural:

```
Distance bins [0, 50, 100, 150, 200]
    ↓
Temporal bins (inter-spike intervals in ms)
    ↓
Spiking neurons with STDP:
    if (Δt ∈ [0, 10ms]): strengthen synapse (LTP)
    if (Δt ∈ [10, 50ms]): weaken synapse (LTD)
```

Your tensor `M[i][j][dist_bin]` becomes `W[i][j][Δt_bin]` (synaptic plasticity rule)!

**4. WebAssembly Deployment**

100% achievable:
- C code compiles to WASM with Emscripten
- SDL2 has web backend
- Interactive demos on GitHub Pages
- Parameter sliders for real-time exploration

---

#### 🎯 Implementation Roadmap

**Immediate (Today): Distance Bins**
- Foundation for all future work
- 3-5 distance bins
- Test oscillating attraction/repulsion

**Phase A (Week 1): Enhanced Tensors**
- Add velocity bins (4D tensor)
- Add angular bins (5D tensor)  
- CLI for loading custom tensors

**Phase B (Week 2): Information Theory**
- Mutual Information calculation
  ```c
  MI(type_i, type_j) = measure correlation in clustering
  ```
- Energy tracking
  ```c
  E_total = kinetic + potential (from forces)
  ```
- Entropy measurement
  ```c
  H(positions) = spatial entropy
  H(velocities) = kinetic entropy
  ```

**Phase C (Week 3): IB Integration**
- Export particle trajectories → IB package input
- IB compression parameters → particle interaction rules
- Shared visualization: particles + IB plane

**Phase D (Month 2): Neural Analog**
- Rename 'particles' → 'neurons'
- Distance bins → Inter-spike intervals
- Force → Synaptic current
- Implement STDP learning rule

**Phase E (Month 3): Web Deployment**
- Emscripten compilation
- Interactive web UI
- Parameter exploration tools

---

#### 📊 Mutual Information Calculation (Concrete Plan)

For particle-life system, MI can measure:

**1. Type-Type Mutual Information**
```c
MI(Type_i ; Type_j) = How much does type i's behavior predict type j's behavior?

// Measure co-clustering
count[i_cluster][j_cluster]++;  // Joint distribution
MI = ΣΣ P(i,j) * log(P(i,j) / (P(i)*P(j)))
```

**2. Position-Velocity MI**
```c
MI(Position ; Velocity) = Spatial structure predicts motion?

// Bin position and velocity
pos_bins[x_bin][y_bin]++;
vel_bins[vx_bin][vy_bin]++;
MI = measure correlation
```

**3. Energy-Configuration MI**
```c
MI(Energy ; Configuration) = Thermodynamic analog

// Configuration = particle arrangement
// Energy = kinetic + potential
MI = how much does configuration predict energy?
```

**4. Rule-Outcome MI** (Meta-level)
```c
MI(Interaction_Rules ; Emergent_Behavior) = Learning!

// Try different rule tensors
// Measure resulting behavior (clustering, orbits, chaos)
// Find which rules maximize interesting behavior
```

---

#### 🔗 IB Package Integration Concept

**Shared Codebase Architecture**:

```
CCS_complex_coupled_systems/
├── IB_DIB_PACKAGE/              # Information Bottleneck
│   ├── src/
│   │   ├── mutual_information.c  # ← SHARED with particle-life
│   │   └── entropy_calc.c        # ← SHARED
│   
├── HardDisks/hspist3/           # Particle-life
│   ├── 00ALLINONE.c
│   ├── particlelife_to_ib.c     # NEW: Export particle data
│   └── ib_to_particlelife.c     # NEW: Import IB rules
│
└── shared/                       # NEW: Common library
    ├── information_theory.c      # MI, entropy, KL divergence
    ├── tensor_ops.c              # Multi-dimensional tensor math
    └── visualization_bridge.c    # Convert between representations
```

**Workflow Example**:
1. Run particle-life with rule tensor
2. Export: particle trajectories → IB input format
3. IB finds compressed representation (relevant features)
4. Import: IB compression rules → new interaction tensor
5. Run particle-life with learned rules
6. **Iterate**: co-evolution of rules and behavior!

---

#### 🧠 Neuroplasticity Connection

**Concrete mapping**:

| Particle-Life | Spiking Neural Network |
|---------------|------------------------|
| Particle type | Neuron population |
| Position (x,y) | Neuron state (V_membrane, refractoriness) |
| Velocity (vx,vy) | Spike rate, spike timing |
| Distance bin | Inter-spike interval (ISI) |
| Force M[i][j][dist] | Synapse weight W[i][j][Δt] |
| Force update | STDP plasticity rule |
| Clustering | Synchronization |
| Energy | Action potential cost |

**Example STDP rule as tensor**:
```c
// Spike-timing dependent plasticity
W[pre][post][Δt < 0] = +A_plus * exp(Δt/τ_plus);   // Pre before post → LTP
W[pre][post][Δt > 0] = -A_minus * exp(-Δt/τ_minus); // Post before pre → LTD

// This is EXACTLY a distance-binned tensor!
// Distance = time delay between spikes
```

---

#### ✅ User Questions Answered

**Q: Can we calculate MI and energy?**
**A**: Absolutely! Particle-life is perfect for this:
- Positions/velocities = observable random variables
- Interaction rules = probability distributions
- Energy = sum of kinetic + potential
- MI = measure correlations in clustering/motion

**Q: Link to IB package?**
**A**: Yes! Natural connection:
- Particle types = information sources
- Clustering = compression
- Interaction strength = information preservation
- Shared tensor mathematics

**Q: Build neuroplasticity simulation?**
**A**: Perfect progression:
- Distance bins → ISI bins
- Forces → Synaptic currents  
- Tensor updates → Learning rules
- Energy minimization → Network training

**Q: Emscripten for web?**
**A**: Standard workflow:
```bash
emcc 00ALLINONE.c -o particlelife.html -s USE_SDL=2
# Runs in browser with WebGL
```

---

### Synthesis (Round 2):

**Status**: 🚀 **Vision Expanded Dramatically**

**Immediate Action**: Implement distance bins as foundation for entire research program

**Long-term Vision**: Unified framework connecting:
1. Correlation matrices (statistics)
2. Particle-life (physics)
3. Information theory (communication)
4. IB compression (machine learning)
5. Neural networks (neuroscience)

**Key Insight**: Rule tensors are the **universal language** connecting all these domains!

**Next Steps**:
1. ✅ CC: Implement distance-binned tensors (start today)
2. Document tensor format for future extensions
3. Design MI calculation module (shared with IB)
4. Plan web deployment (Emscripten test)

---



## Round 3: Predator-Prey Dynamics & Multi-Dimensional Tensors

**Initiated by**: User
**Date**: 2025-10-31

### User Feedback:

**Observation**: "Clusters are stable, almost nothing moving"
- Distance-only tensors reach equilibrium
- Need more dynamic interactions

**Requests**:
1. Add other tensor dimensions (velocity, angular)
2. Implement predator-prey example: Red chases Green, Green flees Red
3. Custom tensor loading from files

---

### CC Response: Multi-Dimensional Tensors Implemented

**Date**: 2025-10-31

#### Implemented Tensor Architecture

**Tensor Modes**:
```c
typedef enum {
    TENSOR_MODE_MATRIX = 0,      // 2D: M[i][j]
    TENSOR_MODE_DISTANCE = 1,    // 3D: M[i][j][dist]
    TENSOR_MODE_VELOCITY = 2,    // 4D: M[i][j][dist][v_rel] (TODO)
    TENSOR_MODE_ANGULAR = 3      // 5D: M[i][j][dist][v_rel][angle] (TODO)
} TensorMode;
```

**Data Structures Added**:
- `pl_tensor_3d[10][10][10]`: Distance-binned (ACTIVE)
- `pl_tensor_4d[10][10][10][5]`: + Velocity bins (structure ready)
- `pl_tensor_5d[10][10][10][5][8]`: + Angular bins (structure ready)

#### Predator-Prey Interaction Pattern

**Red (Type 0) = Predator**:
```c
Red→Green interaction (chase):
  [0, 40)    px: +0.8  (chase when very close)
  [40, 80)   px: +1.0  (strong chase)
  [80, 120)  px: +0.7  (moderate chase)
  [120, 160) px: +0.4  (weak chase)
  [160, 200) px: +0.2  (very weak chase)
```

**Green (Type 1) = Prey**:
```c
Green→Red interaction (flee):
  [0, 40)    px: -1.0  (FLEE when very close!)
  [40, 80)   px: -0.8  (strong flee)
  [80, 120)  px: -0.5  (moderate flee)
  [120, 160) px: -0.2  (weak flee)
  [160, 200) px: -0.1  (very weak flee)
```

**Blue (Type 2) = Neutral**:
- Weak self-clustering only
- No interaction with Red or Green

#### Expected Dynamics

**Complex emergent behavior**:
1. Red particles chase Green particles
2. Green particles flee from Red
3. Blue particles form stable clusters (observers)
4. **Perpetual motion**: Chase/flee prevents equilibrium!
5. Possible outcomes:
   - Red surrounds Green (capture)
   - Green escapes to other side (cyclic boundaries)
   - Dynamic orbiting (Red chases, Green circles around)

---

### Test Instructions

```bash
./00ALLINONE --experiment=particlelife
```

**Watch for**:
- 🔴 Red particles actively chasing 🟢 Green
- 🟢 Green particles fleeing from 🔴 Red
- 🔵 Blue particles clustering peacefully
- Continuous motion (no static equilibrium!)
- Emergent chase patterns

**Console Output**:
```
📊 Interaction TENSOR (3x3x5) - PREDATOR-PREY MODE:
   🔴 Red = Predator (chases green)
   🟢 Green = Prey (flees from red)
   🔵 Blue = Neutral (weak clustering)
   Distance bins: [0,40) [40,80) [80,120) [120,160) [160,200)
   Type 0→0:   0.10   0.10   0.30   0.10   0.10
   Type 0→1:   0.80   1.00   0.70   0.40   0.20   ← Red chases Green!
   Type 0→2:   0.00   0.00   0.00   0.00   0.00
   Type 1→0:  -1.00  -0.80  -0.50  -0.20  -0.10   ← Green flees Red!
   Type 1→1:   0.10   0.10   0.30   0.10   0.10
   ...
```

---

### Next Steps

**Phase 1: Velocity Tensors** (Next session)
- Implement 4D tensor lookup in force calculation
- Chase depends on relative velocity
- Fast predator + slow prey = strong attraction
- Slow predator + fast prey = weak attraction (can't catch)

**Phase 2: Angular Tensors** (Next session)
- Implement 5D tensor lookup
- Vision cones: Predator only sees prey in front
- Prey has 360° awareness
- Directional fleeing (away from predator)

**Phase 3: Custom Tensor Loading** (This session if time)
- File format specification
- CLI: `--pl-tensor=rules.txt`
- Example tensor files for different behaviors

**Phase 4: Information Theory Integration**
- Measure MI between particle types
- Energy tracking (kinetic + potential)
- Entropy of spatial distributions

---

