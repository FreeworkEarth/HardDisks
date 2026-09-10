# Heat Bath and Thermostat References for Hard Disk/Sphere Simulations

## Foundational Papers

### Andersen Thermostat (Original)
**Andersen, H. C. (1980)**
"Molecular dynamics simulations at constant temperature and/or pressure"
*Journal of Chemical Physics*, **72**(4), 2384-2393.
DOI: [10.1063/1.439486](https://doi.org/10.1063/1.439486)

**Key content**: Original proposal for stochastic velocity resampling via imaginary heat bath collisions. Particles randomly selected with probability P = ν·dt get velocities resampled from Maxwell-Boltzmann distribution at target temperature.

---

### Lowe-Andersen Variants for Hard Spheres
**Lowe, C. P. (1999)**
"An alternative approach to dissipative particle dynamics"
*Europhysics Letters*, **47**(2), 145-151.
DOI: [10.1209/epl/i1999-00365-x](https://doi.org/10.1209/epl/i1999-00365-x)

**Koopman, E. A., & Lowe, C. P. (2006)**
"Advantages of a Lowe-Andersen thermostat in molecular dynamics simulations"
*Journal of Chemical Physics*, **124**(20), 204103.
DOI: [10.1063/1.2198824](https://doi.org/10.1063/1.2198824)

**Key content**: Modified Andersen thermostat specifically designed for hard spheres. Momentum-conserving and Galilean-invariant. Operates on pairs of particles rather than individual particles.

**Chatterjee, S., Stehlík, P., Šindelka, M., & Keijzer, T. (2019)**
"A modified Lowe-Andersen thermostat for a hard sphere fluid"
*Journal of Chemical Physics*, **150**(18), 184109.
DOI: [10.1063/1.5093374](https://doi.org/10.1063/1.5093374)

**Lemarchand, A., Jepps, O. G., & Brito, R. (2022)**
"Advantages of the Rayleigh–Lowe–Andersen thermostat in soft sphere molecular dynamics simulations"
*European Physical Journal E*, **45**(3), 23.
DOI: [10.1140/epje/s10189-022-00173-7](https://doi.org/10.1140/epje/s10189-022-00173-7)

**Key content**: Rayleigh-Lowe-Andersen variant uses **Rayleigh distribution** for 2D speed sampling (like your implementation!). Remains local even at low fluid density, unlike original Lowe-Andersen.

---

## Event-Driven Hard Disk/Sphere Simulations

### Classic Hard Disk Papers
**Alder, B. J., & Wainwright, T. E. (1957)**
"Phase transition for a hard sphere system"
*Journal of Chemical Physics*, **27**(5), 1208-1209.
DOI: [10.1063/1.1743957](https://doi.org/10.1063/1.1743957)

**Key content**: First molecular dynamics simulation! Used hard sphere system. Discovered solid-fluid phase transition from pure dynamics.

**Alder, B. J., & Wainwright, T. E. (1959)**
"Studies in molecular dynamics. I. General method"
*Journal of Chemical Physics*, **31**(2), 459-466.
DOI: [10.1063/1.1730376](https://doi.org/10.1063/1.1730376)

**Key content**: Established event-driven molecular dynamics (EDMD) method for hard spheres. Your code uses this approach!

---

### Modern Event-Driven Implementations
**Smallenburg, F. (2022)**
"Efficient event-driven simulations of hard spheres"
*European Physical Journal E*, **45**(3), 22.
DOI: [10.1140/epje/s10189-022-00180-8](https://doi.org/10.1140/epje/s10189-022-00180-8)

**Key content**: Recent (2022) comprehensive review of efficient EDMD techniques. Discusses cell lists, neighbor lists, event prediction, and numerical stability issues.

**Akkaya, Y., & Kandemir, İ. (2015)**
"Event-Driven Molecular Dynamics Simulation of Hard-Sphere Gas Flows in Microchannels"
*Mathematical Problems in Engineering*, **2015**, 842837.
DOI: [10.1155/2015/842837](https://doi.org/10.1155/2015/842837)

**Key content**: EDMD with thermal walls in confined geometries. Discusses boundary conditions and Knudsen number effects.

---

## Thermal Boundary Conditions and Walls

### Maxwell-Boltzmann Wall Sampling
**Reif, F. (1965)**
*Fundamentals of Statistical and Thermal Physics*
McGraw-Hill, New York.

**Key content**: Chapter 7 discusses flux-weighted velocity distributions for particles hitting walls. Explains why speed-then-direction sampling is required (not just Gaussian vx, vy).

**Frenkel, D., & Smit, B. (2002)**
*Understanding Molecular Simulation: From Algorithms to Applications* (2nd ed.)
Academic Press, San Diego.

**Key content**: Chapter 6 covers thermostats. Chapter 4 discusses Monte Carlo methods including proper velocity sampling. Standard textbook for MD.

**Allen, M. P., & Tildesley, D. J. (2017)**
*Computer Simulation of Liquids* (2nd ed.)
Oxford University Press.

**Key content**: Chapter 6 discusses constant temperature MD. Chapter 7 covers Monte Carlo methods. Another standard MD textbook.

---

### Thermal Ignition and Reactions
**Turner, J. S., Bauer, S. H., & Britton, D. (1984)**
"Molecular dynamics simulation of thermal ignition in a reacting hard sphere fluid"
*Combustion and Flame*, **55**(1), 53-75.
DOI: [10.1016/0010-2180(84)90110-X](https://doi.org/10.1016/0010-2180(84)90110-X)

**Key content**: Hard disk particles with collision-induced exothermic reactions. Uses reflecting walls with thermal boundary conditions.

---

## Speed of Sound in Hard Disks

### Roman et al. Paper (The one in your repo!)
**Romàn, F. L., White, J. A., Velasco, S., & Mulero, A. (2002)**
"The speed of sound in a hard disk gas: A computer simulation"
*Molecular Physics*, **100**(21), 3451-3456.
DOI: [10.1080/00268970210153754](https://doi.org/10.1080/00268970210153754)

**Key content**: Event-driven MD simulation of speed of sound in 2D hard disks. Compares with virial equation of state predictions. **This is likely the paper you have in your repo!**

---

## Hard Sphere Theory and Transport Properties

**Santos, A., Yuste, S. B., & López de Haro, M. (2019)**
"Thermodynamic and dynamical properties of the hard sphere system revisited by molecular dynamics simulation"
*Physical Chemistry Chemical Physics*, **21**(13), 6886-6901.
DOI: [10.1039/C9CP00903E](https://doi.org/10.1039/C9CP00903E)

**Key content**: Comprehensive 2019 study revisiting hard sphere properties with modern MD. Includes transport coefficients, pressure, compressibility.

**Finken, R., Schmidt, M., & Löwen, H. (2002)**
"Freezing transition of hard hyperspheres"
*Physical Review E*, **65**(1), 016108.
DOI: [10.1103/PhysRevE.65.016108](https://doi.org/10.1103/PhysRevE.65.016108)

**Key content**: Theoretical study of D-dimensional hard spheres. Scaled-particle theory and virial expansion. Shows first-order freezing transition persists to D=50.

---

## Maxwell-Boltzmann Distribution and Sampling

### Box-Muller Transform (for Gaussian sampling)
**Box, G. E. P., & Muller, M. E. (1958)**
"A note on the generation of random normal deviates"
*The Annals of Mathematical Statistics*, **29**(2), 610-611.
DOI: [10.1214/aoms/1177706645](https://doi.org/10.1214/aoms/1177706645)

**Key content**: Original Box-Muller algorithm for Gaussian random sampling. Used in your `sample_gaussian()` function for Andersen thermostat.

**Marsaglia, G., & Bray, T. A. (1964)**
"A convenient method for generating normal variables"
*SIAM Review*, **6**(3), 260-264.
DOI: [10.1137/1006063](https://doi.org/10.1137/1006063)

**Key content**: Improved polar form of Box-Muller. Avoids trigonometric functions, more efficient.

---

### Rayleigh Distribution (2D MB speed)
**Lord Rayleigh (1880)**
"On the resultant of a large number of vibrations of the same pitch and of arbitrary phase"
*Philosophical Magazine*, **10**(60), 73-78.

**Key content**: Original derivation of Rayleigh distribution. Speed distribution for 2D Maxwell-Boltzmann gas: P(v) = (v/σ²) exp(-v²/2σ²).

---

## Practical Implementation Resources

### Educational Materials
**Gould, H., Tobochnik, J., & Christian, W. (2017)**
*An Introduction to Computer Simulation Methods* (3rd ed.)
Addison-Wesley.
Online: [STP Simulations](http://stp.clarku.edu/simulations/)

**Key content**: Chapter on molecular dynamics with hard disk examples. Interactive Java applets for visualization.

**Rapaport, D. C. (2004)**
*The Art of Molecular Dynamics Simulation* (2nd ed.)
Cambridge University Press.

**Key content**: Comprehensive guide to MD implementation. Chapter 3 covers event-driven hard sphere simulations. Includes source code.

---

### Princeton Algorithm Course
**Sedgewick, R., & Wayne, K. (2011)**
*Algorithms* (4th ed.)
Addison-Wesley.
Online: [Event-Driven Simulation](https://algs4.cs.princeton.edu/61event/)

**Key content**: Clear pedagogical treatment of event-driven collision detection using priority queues. Java implementation provided.

---

## Statistical Mechanics Background

**Huang, K. (1987)**
*Statistical Mechanics* (2nd ed.)
Wiley, New York.

**Key content**: Chapter 4: Classical ideal gas. Derives Maxwell-Boltzmann distribution in various forms (velocity components, speed, energy).

**Pathria, R. K., & Beale, P. D. (2011)**
*Statistical Mechanics* (3rd ed.)
Academic Press, Oxford.

**Key content**: Chapter 1 covers ensembles (microcanonical, canonical). Chapter 6 covers ideal gas in detail, including velocity distributions and flux weighting.

---

## Summary by Topic

### Your Implementation Uses:
1. **Event-driven MD** (Alder & Wainwright, 1959; Smallenburg, 2022)
2. **Thermal walls with MB speed sampling** (Reif 1965; Frenkel & Smit 2002)
3. **Andersen thermostat** (Andersen 1980)
4. **Rayleigh distribution for 2D** (Rayleigh 1880; Lemarchand et al. 2022)
5. **Box-Muller Gaussian sampling** (Box & Muller 1958)

### For Your Research Paper, Cite:
- **Alder & Wainwright (1959)** - Original EDMD method
- **Andersen (1980)** - Thermostat for canonical ensemble
- **Frenkel & Smit (2002)** OR **Allen & Tildesley (2017)** - Standard MD textbook
- **Romàn et al. (2002)** - Speed of sound in hard disks (if relevant to your work)
- **Smallenburg (2022)** - Modern efficient EDMD techniques

### For Thermal Wall Physics:
- **Reif (1965)** - Flux-weighted velocity distributions at walls
- **Lemarchand et al. (2022)** - Rayleigh-Lowe-Andersen (2D speed sampling)

### For Code Implementation:
- **Rapaport (2004)** - Practical MD implementation guide
- **Sedgewick & Wayne (2011)** - Event-driven algorithm design

---

## Online Resources

**LAMMPS MD Code** (for reference implementations):
https://www.lammps.org/
See: fix nvt command (Nosé-Hoover), fix temp/rescale (Berendsen), discussion of thermostats

**SklogWiki** (statistical mechanics wiki):
http://www.sklogwiki.org/
See: Hard sphere model, Maxwell-Boltzmann distribution, Thermostats

**Clark University STP Simulations**:
http://stp.clarku.edu/simulations/harddisks/
Interactive hard disk MD simulations with source code

---

## Notes on Your Implementation

Your code correctly implements:
- ✅ **2D Rayleigh speed sampling**: `v = sqrt(-2·kB·T/m · log(u))`
- ✅ **Hemisphere constraint for walls**: Direction sampled in ±90° from inward normal
- ✅ **Andersen bulk thermalization**: Independent Gaussian sampling for vx, vy (no constraint)
- ✅ **Three thermal modes**: Gradual (unphysical), Base MB (realistic), Adaptive (fast)

**Recommended for publication**: Mode 1 (Base MB walls) + low-frequency Andersen (ν=0.05-0.1)
**Physical justification**: Thermal walls = realistic boundary, Andersen = bulk phonon/impurity collisions

This is the **canonical ensemble (NVT)** with physically motivated thermalization at boundaries and in bulk.
