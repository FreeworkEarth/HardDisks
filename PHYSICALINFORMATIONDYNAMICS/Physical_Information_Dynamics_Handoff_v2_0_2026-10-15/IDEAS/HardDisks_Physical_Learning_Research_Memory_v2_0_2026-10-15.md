---
title: "From the Román piston to trainable physical networks — HardDisks research memory"
date: 2026-10-12
updated: 2026-10-15
version: "2.0"
type: research-memory-and-hypothesis-ledger
status: living-exploratory-document
primary_research: HardDisks/hspist3
related_research:
  - Román thermal piston and equilibrium/dynamic stiffness
  - hard-disk energy-transfer boxes and spring capture
  - spring-damper mechanical learning networks
  - physical reservoir computing and neuromorphic systems
  - predictive coding and dynamical information bottleneck
  - Physarum-inspired morphological computation
  - empowerment, adaptive networks, physical energy accounting
  - DLA growth and branching, neuronal growth cones, morphogenesis, structural learning
knowledge_workspace: IdeaLabs
other_distinct_projects:
  - NeurAIpil
origin: "Research discussions through 2026-10-12; related 2026-10-08 IdeaLabs research brief"
labels: [hard-disks, roman-piston, energy-transfer, springs, stiffness, physical-neural-network, reservoir-computing, information-bottleneck, empowerment, physarum, neuroplasticity, thermodynamics, research-hypotheses]
---

# From the Román piston to trainable physical networks

**A living research memory: the ideas, questions, connections, scientific restrictions, sources, and possible experiments emerging from our discussion through 12 October 2026.**

> **Guiding question:** Can we build a physically trainable system out of characterized, mechanically driven hard-disk compartments—possibly inspired by axonal transmission, spring networks, and Physarum routing—and show that its physical dynamics or adaptive structure perform a useful predictive task **at a defensible total energetic cost**?

**Additional central question (2026-10-12): Where is intelligence physically encoded—in fast physical dynamics, slowly adapted morphology, local learning rules, or their interaction with the environment? When does building and maintaining additional network complexity pay for itself?** See §10.5, §§17–21 for branch generation and contacts, and §§22–27 for the dynamic-radius structural-disk model; §§28–34 for the five-paper/Fable handoff; §§35–39 for new literature and the Physarum fan/nuclear-trajectories clarification, and **§40 for a deep prior-art/novelty-risk audit connecting Still, empowerment and physical growth costs**.

This is **not** a proof of a new computing principle, a claim that a hard-disk box is a neuron, a finished experimental plan, or permission to interrupt the currently validated HardDisks research. Its job is to make sure the questions do not disappear and to connect them to experiments that might eventually answer them. It complements, rather than replaces, the 2026-10-08 `IdeaLabs_Research_Brief_Physarum_Information_Energy_Morphological_Computing_2026-10-08.md`.

## 0. How to read this notebook

Evidence/status labels:

- **[ESTABLISHED / SOURCE]** Published or standard result, subject to the cited paper's assumptions.
- **[PROJECT STATE / REPORTED]** Work or apparatus configuration described in our project discussions; independently check raw runs and commits before putting a number in a paper.
- **[DERIVATION / ASSUMPTIONS]** Mathematical relation valid under specific hypotheses, not a universal identity.
- **[ANALOGY]** Functional resemblance, **not** shared biological or physical mechanism.
- **[HYPOTHESIS]** An empirically refutable suggestion.
- **[EXPERIMENT PROPOSAL]** Not yet implemented or funded as part of the research program.
- **[VISION]** Longer-term research direction.
- **[OPEN]** An unresolved question or incomplete literature review.

**Do not collapse these categories.** In particular: a visually compelling simulation, an oscillation, a delayed mechanical response, or a changing network does not automatically demonstrate learning, computation, energy dissipation, a refractory period, or energy efficiency.

### The user's questions worth preserving in their original spirit

> “is the energy trasnfer box rather the soma or the axon with myelin and soma shall be something different which activates piston pushes?”

> “and isnt a spring and damper better or is our box what we want with different walls inside maybe these show the dead zone after an action potential.”

> “the roman box could also be a very easy spring system with measured stiffness maybe.”

> “maybe a trainable network from these boxes”

> “I want science here not some random stuff.”

> “these are exactly the questions I have and collect together and see parallel lines to our boxes.”

> “where is the intelligence encoded in nature” — question the user highlighted at approximately 4:04 of *Lens of Time: Slime Lapse* (bioGraphic).

> “i want to figure connections between thermodyn. information IB anad empowerment and adaptive learnable networks and when complexity is worth spending the energy to grow network”

Other enduring questions: Why strengthen some paths and prune others? Could a lossy mechanical channel carry task-useful memory? Can different box geometries implement different temporal filters? When is retaining an extra pathway worth its maintenance cost? Can a Physarum-like system choose and reconfigure such paths, and can a predictive Information Bottleneck specify which memories should survive?

---

## 1. The research program in one view

```text
PAPER 1: characterize microscopic medium and mechanical response
    2D EDMD hard disks / Román thermal piston
    force F(L), k_T, equilibrium k_S, resonant frequency, c_s,
    density, confinement, mode lifetime, structural observables
             |
             v
PAPER 2: controlled physical energy transfer
    driven piston -> hard-disk gas -> movable dividers -> spring wall
    work ledger, wave/pulse transfer, signal fidelity, spring capture
             |
             v
PAPER 3: sensing, predictive compressed state, and actuation
    input histories -> Z_t (Dynamical/Predictive IB) -> controller -> piston
    task-relevant performance versus total energetic cost
             |
             v
PAPER 4 CANDIDATE: coupled boxes, topology, adaptation, memory
    input -> C0 -> Box1 -> C1 -> [Box2 || Box3] -> C2 -> Box4 -> output
    compare fixed, adjustable, split-merge, redundant and grown networks
             |
             v
PAPER 5 CANDIDATE / LONG-TERM VISION
    Physarum-inspired adaptive morphology, empowerment, physical
    reservoir learning, neural-like computation, eventually hardware
```

**Critical dependency:** Paper 1 should supply independently measured mechanical parameters to Paper 2; Paper 2 supplies a validated energetic baseline to Paper 3; later network/computing claims require measured tasks and baselines. The 15-module cinematic network and 4-module split–merge circuit are **proposed concepts**, not completed simulation results.

### Chronological landmarks in the current discussion

- **Before October 2026:** HardDisks/hspist3 2D EDMD spring-wall collision scheduling, Román speed-of-sound design and validation, future piston energy transfer; Information Bottleneck, Physarum and empowerment had already been proposed as distinct research branches.
- **8–10 October 2026:** The earlier IdeaLabs brief consolidated the physics-to-morphology plan, Jones/Tero Physarum models, causal trails versus observational hard-disk traces, and the need for full energetic accounting.
- **11 October 2026:** Film developed as a reverse-complexity descent: real tree/roots, Physarum transport, neurons, engineered mechanical network, and eventually one Román apparatus. The mechanical mapping was refined to **soma/initial segment ≈ controller/trigger** and **selected axonal pathways ≈ transmission modules**, not soma = gas box.
- **11–12 October 2026:** Real HardDisks screenshot established the geometry of multiple gas sections, movable dividers, a terminal spring and an input piston. Film target uses a **left-to-right mirror** for continuity (drive at left, spring at right), whereas the supplied actual simulation screenshot had the spring region at left and drive at the other end. This is a *visual orientation convention*, not a claim the simulated apparatus changed.
- **12 October 2026:** The user brought neuroscience videos, predictive-coding references, spring-network learning literature, physical reservoir computing, dynamic memristors, and neuromorphic scaling into the discussion. An increasingly concrete question emerged: **Is our characterized gas/divider/spring system a useful, trainable physical edge—or could a simpler spring–damper do as well?**

---

## 2. The known physical foundation: Román piston and measured stiffness

### 2.1 The *simple Román system* as an experimentally characterized effective spring

**[ESTABLISHED / PROJECT CONTEXT]** The Román-inspired apparatus is a 2D gas separated by a central movable divider, typically with equal gas compartments on both sides and fixed exterior walls:

```text
fixed outer wall | 2D hard-disk gas | freely movable central divider | 2D hard-disk gas | fixed outer wall
```

For **one** compartment of length `L`, height `H`, `N` identical hard disks of diameter `σ` and mass `m`, at temperature `T`:

\[
\eta=\frac{N\pi\sigma^2}{4LH}, \qquad F(L,T)= \text{mean collision-impulse force on ONE divider face}.
\]

**[DERIVATION / ASSUMPTIONS]** Measure `F(L,T)` while the divider is held fixed during each equilibrated run. The one-sided isothermal stiffness is

\[
\boxed{k_T=-\left(\frac{\partial F}{\partial L}\right)_T.}
\]

For the equilibrium hard-disk system considered, the *fully relaxed adiabatic* stiffness of **one compartment** is

\[
\boxed{k_{S,\mathrm{eq}} = k_T + \frac{F^2}{N k_B T}.}
\]

For a symmetric divider where both sides resist displacement, use the sum of the two sides' restoring stiffness contributions (twice the one-side value in the symmetric limit). **Do not mix one-side stiffness with total-divider stiffness.**

**Core new conceptual insight [ANALOGY, TESTABLE]:** At small displacement and sufficiently slow perturbation, a gas compartment **is an effective spring with a measured, state-dependent stiffness**. There is no need to start by claiming a complicated neural-like response. The first rung is to establish that the fully microscopic 2D gas reproduces a familiar mechanical constitutive element in a specified regime.

For the free-divider acoustic experiment, the Román idealization relates mode frequency `f_*`, mass ratio and apparent sound speed through

\[
\cot K = \frac{M}{2 N m} K, \qquad c_{\rm app}= \frac{2\pi L f_*}{K}.
\]

In an appropriate homogeneous acoustic regime, one-compartment acoustic stiffness is often written

\[
 k_{S,\mathrm{ac}} \simeq \frac{Nm c_s^2}{L^2}.
\]

**Important:** The acoustic expression is not automatically identical to finite-box static equilibrium stiffness for every finite `N`, mass, density, frequency or coexistence state. Long-wavelength assumptions, mode geometry, resolution and relaxation matter. Román mode frequency measures a coupled piston–gas mode, not a standalone intrinsic material modulus.

### 2.2 What about damping, viscosity and a simple spring–damper?

For a *candidate coarse-grained local model* of one wall, around equilibrium,

\[
 M\,\delta\ddot x + \gamma\,\delta\dot x + k_{\rm eff}\,\delta x \approx F_{\rm drive}(t).
\]

This is a **benchmark/approximation**, not an exact replacement for finite event-driven hard-disk dynamics. In a closed elastic system, energy can move out of a selected coherent wall mode into the many-body degrees of freedom or return later. Apparent envelope decay or phase mixing is **not automatically irreversible dissipation into a bath**; effective `γ` depends on the measured observation, time window, coupling and thermodynamic assumptions.

The more general response may be history-dependent:

\[
F_{\rm gas}(t) \approx F_0 - \int_{-\infty}^{t} K(t-s)\,\delta x(s)\,ds,
\]

with an effective memory kernel `K`, if such linear response is established. A frequency-dependent complex stiffness `k*(ω)` may be a useful fitted representation, but any inferred relaxation time needs independent validation.

**Question:** Can we measure when the hard-disk chamber behaves like a simple spring, a spring–damper, a dispersive multi-mode resonator, or a genuinely nonlinear, history-dependent transmission element? Define regimes, not metaphors.

### 2.3 Current validation status and non-overclaiming

**[PROJECT STATE / REPORTED]** The EDMD spring-wall exact harmonic collision scheduling, oscillation-target data acquisition and FFT robustness were being validated; a historical low-density discrepancy was suspected to be partly a frequency-resolution artifact. Near liquid–hexatic densities, a finite-frequency discrepancy is not automatically evidence of structural relaxation. Reproduce EDMD failing seeds and establish estimator invariance before building claims on damping, resonance or a new phase effect. Read the hspist3 development-memory index and current experiments for the authoritative validation state; this Markdown does not supersede those artifacts.

**[SOURCE]** Hard-disk liquid–hexatic transition, 2D ordering and existing EOS are established topics (Bernard–Krauth 2011; Engel et al. 2013). Use local `ψ6` and translational/orientational correlations where relevant, not visual color or center trails alone.

---

## 3. Energy-transfer apparatus: an *edge* of a possible physical network

**[PROJECT STATE / GEOMETRY]** The real visual reference in this conversation shows an empty terminal spring region and multiple gas sections separated by thin movable walls. A faithful schematic of the *original* apparatus orientation:

```text
fixed anchor / spring | spring-backed MOVABLE wall | gas A | MOVABLE divider | gas B | drive piston
```

More dividers can be introduced in a controlled experimental sequence. **Each module is a sealed mechanical system**: particles collide with impermeable dividers but do not pass between separated gas compartments. The moving walls carry forces and displacements. Identical hard disks have equal diameter; an equilibrium 2D gas with no gravity fills each available compartment rather than collecting at its bottom.

The **cinematic mirrored orientation** used for a left-to-right network is:

```text
controller / driven piston -> gas A -> movable divider(s) -> gas B -> movable spring-backed wall / sensor
```

This orientation reversal should be described as a schematic mirror. Neither a 3D hardware shell nor blue CGI particles represent an implemented 3D gas: **the physics engine remains two-dimensional hard disks**.

### 3.1 What the edge transmits

Distinguish at least three things:

1. **Total energy** conserved/repartitioned in an ideal closed elastic model.
2. **Coherent mechanical pulse or wall-mode energy**, which can lose amplitude or fidelity through scattering into other degrees of freedom, dephasing, or reflection, even without an external heat bath.
3. **Task-relevant information** carried by the relation between a coded piston input and measured output wall displacement/velocity (or pressure/force waveform).

The label **“lossy signal transmission”** concerns (2) and/or (3) under a defined input-output task; it is *not* proof that energy mysteriously disappears from the simulation.

With explicit drive and any real environmental damping, use an audited ledger, e.g.

\[
W_{\rm in} = \Delta K_{\rm disks}+\Delta K_{\rm walls}+\Delta U_{\rm spring}+Q_{\rm bath}+R_E,
\]

where bath heat is present **only if modeled**, and `R_E` is a numerical/accounting residual rather than a real reservoir. The spring capture ratio

\[
\epsilon = \frac{\Delta U_{\rm spring}}{W_{\rm in}}
\]

requires a defined measurement endpoint, consistent initial states and `W_in>0`; it does not measure communication reliability by itself.

### 3.2 Design and adaptation parameters

- Geometry: `L`, `H`, gas-compartment lengths, number and positions of dividers.
- Particle medium: `N`, `m`, diameter `σ`, packing fraction `η`, kinetic temperature/energy.
- Movable structure: wall masses, terminal spring constant `k`, offsets/rest lengths.
- Drive: force/displacement protocol, timing, amplitude, frequency, duration, energy budget.
- Measurement: sensor type/precision/rate, readout, filtering, detection threshold.
- Network: which controllers/edges are connected, route duplication, rest/deactivation state, physical switching mechanism and work.

**Not all knobs are easy to change physically.** Changing gas particle number, chamber length, wall mass or topology has a real reconfiguration cost and may require a redesigned apparatus. A clean simulation parameter sweep is not evidence of a feasible low-energy adaptive material.

---

## 4. Neuronal parallels: a useful map with explicit limits

### 4.1 Functional map we selected for the film and research questions

| Biological concept | Candidate engineered analogue | Limitation / test required |
|---|---|---|
| Dendritic synaptic input | Coarse measurements arriving at controller | Biological dendrites themselves perform nonlinear local computations; do not reduce all dendrites to passive wires. |
| Soma + axon initial segment | Sensor/integration/threshold/actuator controller | Real spike initiation typically depends on axon initial segment; piston trigger requires external energy. |
| Axonal propagation | Driven hard-disk transmission module (`edge`) | Axonal AP is actively regenerated via voltage-gated ion channels; a passive gas box is not. |
| Myelination | Tunable transmission reliability, latency or confinement (functional analogy only) | Myelin changes membrane capacitance/resistance and enables saltatory conduction; gas-wall tuning is a different mechanism. |
| Synapse | Output measurement + coupling to receiving controller | Chemical synapses release transmitter across cleft; rods/pipes are not chemical synapses. |
| Synaptic efficacy | Controller coupling gain or adjustable transmission parameters | A thicker axon is **not** a universal indicator of stronger synapse. |
| Structural plasticity | Building/pruning/switching a route, changing stiffness/geometry | Must explicitly count time, sensing/material/actuation cost. |
| Refractory period | Hypothesized state-dependent suppressed response to closely spaced inputs | Not equivalent to ordinary spring oscillation or delay. |

### 4.2 Refractoriness, “dead zones,” active regeneration

**[ESTABLISHED]** Biological refractory behavior arises from channel-state dynamics (especially sodium-channel inactivation, with potassium conductance/recovery effects). Myelinated axons use electrically insulated segments and nodes of Ranvier; a propagating action potential is **actively regenerated**. A passive spring/damper or gas box may attenuate, ring or delay a pulse, but it has no automatic ionic excitation mechanism.

**[HYPOTHESIS]** Could a multiple-divider gas apparatus exhibit a temporary **input-history-dependent period** in which it transmits a second coded pulse less reliably than the first? If yes, is that a recoverable nonlinear state dependence or merely linear superposition, phase cancellation, geometry or depleted spring energy?

A plausible neuron-inspired engineered element might eventually consist of a **passive transmission body + sensor + threshold logic + powered actuator**. That creates an active regeneration function and a measurable cost. It is proposed, not present in the passive hard-disk chambers by default.

### 4.3 Distinguish neuronal morphology pictures from gas dynamics

The supplied *Neurons & Synapses* video evokes the same visual vocabulary—many tiny elements, branching networks, local compartments—but reconstructed neuron processes are **not particles of an ideal gas**, and reveal/zoom animations are not proof of growth. The science connection should come from independently measured *transmission, memory, learning and energy* functions, not superficial visual similarity.

Other neuroscience leads in this conversation: different apical and basal dendrite plasticity rules; predictive coding/local error learning; the distinction between chemical synaptic adaptation and axon conduction; optogenetics as a way to manipulate excitability (not a model of mechanical-gas computation). Maintain specific paper provenance before treating educational videos as primary evidence.

---

## 5. Compare against the obvious simpler baseline: spring–damper networks

**User insight:** *Maybe our Román box is already a spring—with a stiffness we measured. Why not simply build trainable spring networks? Is the gas box useful, or needlessly complicated?*

This is the most important experimental control for a future physical-computing claim.

| Candidate substrate | Physics | Measured or trainable parameters | What it might provide | Limitation |
|---|---|---|---|---|
| A. Ideal Hookean spring | `F=-kx` near equilibrium | `k` | Static stiffness and simple storage | No rich temporal memory on its own. |
| B. Linear spring–mass–damper | `m x¨+γx˙+kx=u(t)` | `m,k,γ` | Resonance, filtering, decaying temporal state | State memory ≠ learning; conventional baseline. |
| C. Román gas + divider | Many-body collisions and moving central divider | `k_T,k_S,f_*,c_s` measured from medium | Microscopic origin of effective stiffness/mode | Must validate finite-size and frequency effects. |
| D. Driven multi-divider hard-disk box | Repeated collisions, gas–wall coupling, spring capture | Geometry, divider masses, density, spring, input | Multiple modes, nonlinear dynamics, time-dependent response | May add noise, inefficiency and control cost, not a guaranteed advantage. |
| E. Sensor + controller + box | Above plus information processing and supplied actuation | Observation/encoding policy, feedback law, drive | Task-specific amplification, gating, prediction | Energy accounting and fairness crucial. |
| F. Reconfigurable network of boxes | Multiple C–box–C edges, switchable structure | Topology, route activity, per-edge parameters | Redundancy, alternate routes, possible physical memory | Adaptation/maintenance costs and feasible hardware unknown. |

A **simple spring is not automatically a useful learning system**, just as a hard-disk chamber is not automatically a neuron. Both need inputs, readable outputs, defined tasks and trainable degrees of freedom to test learning. If a linear baseline does the same task with lower cost, that is a valuable negative result.

### 5.1 Why gas dynamics *could* add something

**[HYPOTHESIS]** A multi-divider chamber may provide richer time-frequency filtering, nonlinearity, microscopic fluctuations, and state-dependent responses than a single linear spring or spring–damper. This might support temporal discrimination, delayed response, nonlinear feature mixing, or physical reservoir computing.

**[ALTERNATIVE]** These degrees of freedom may simply waste input work as inaccessible coherent motion, increase variability and make training difficult. Complexity by itself is not a performance metric.

### 5.2 Does a spring network have negative “weights”?

For a passive linear Hookean spring, the stiffness `k` is nonnegative in its ordinary stable configuration. A neural-network **negative effective weight** (e.g. −1) need not correspond to a literal negative spring constant; it can arise through displacement conventions, topology, mechanical inversions, antagonistic motions, or active/architected components. Do not conclude `k∈[-1,1]` from neural-weight diagrams. What a published device does must be verified from the particular paper/video.

---

## 6. Mechanical networks that already learn: what the uploaded papers actually support

### 6.1 Li & Mao (2024): *Training all-mechanical neural networks for task learning through in situ backpropagation*

**[SOURCE]** Nature Communications 15, 10528 (2024), DOI https://doi.org/10.1038/s41467-024-54849-z. The model uses mechanical nodes joined by **linear springs with adjustable nonnegative stiffness**. In its quasi-static force-response formulation:

\[
D(\mathbf k)\,\mathbf u=\mathbf F,
\qquad D=C^\mathsf{T}\operatorname{diag}(\mathbf k)C.
\]

The paper derives an **adjoint, local-strain-gradient** route to training stiffnesses. In its chosen conventions, the reported local update is built from the product of the forward bond elongation `e_i` and the adjoint bond elongation `e_i^adj` (see its equations 4–7). It experimentally validates gradient access and reports numerical learning for specified behavior/classification tasks; do not casually describe all tasks as fully autonomous experimentally trained hardware.

**What transfers:** Define controllable per-edge parameters and a task loss; ask whether local observables allow learning without differentiating the entire simulated many-body trajectory.

**What does not transfer:** Their exact quasi-static spring-network local gradient identity is **not automatically applicable** to stochastic EDMD, nonlinear collision dynamics, irreversibility, driven controllers or time-dependent tasks. Adjoint sensitivity through events would require separate technical treatment.

**High-priority question:** Can measured transfer functions / coarse-grained differentiable surrogates provide usable gradients for a **network of hard-disk transmission boxes** while preserving physical predictive accuracy?

### 6.2 Wright et al. (2022): *Deep physical neural networks trained with backpropagation*

**[SOURCE]** Nature 601, 549–555 (2022), DOI https://doi.org/10.1038/s41586-021-04223-6. Proposes and demonstrates physics-aware training of layered controllable physical systems, including mechanical, optical and electronic substrates. Crucial distinction: a physical layer need not implement a mathematical neural layer literally for training to work.

**Relevance [PROPOSAL]:** A measured or empirically differentiable box-response model could serve as a physical layer; use a simulator-assisted surrogate or physics-aware training only after verifying prediction error and domain of validity. Digital training and physical deployment must have separately accounted energy costs.

### 6.3 Stepney (2024): *Physical reservoir computing: a tutorial*

**[SOURCE]** Natural Computing 23, 665–685 (2024), DOI https://doi.org/10.1007/s11047-024-09997-y. The tutorial emphasizes the distinction between a dynamical object and a **demonstrated computation**. The physical substrate needs encoded inputs, measurable outputs, a representation/readout and task-independent or task-dependent evaluation; generic nonlinear state evolution is not enough.

**Relevance [PROPOSAL]:** Before attempting full adaptive architecture, drive one box with temporal inputs and measure whether a simple trained readout can perform a held-out task. This is a lower-risk intermediate test than demanding self-learning spring stiffness or growing neural morphology.

### 6.4 Du et al. (2017): *Reservoir computing using dynamic memristors for temporal information processing*

**[SOURCE]** Nature Communications 8, 2204 (2017), DOI https://doi.org/10.1038/s41467-017-02337-y. Exploits physical ionic dynamics and short-term memory in a memristor array as the time-dependent substrate, with trained readout. Useful example of mapping a time-varying input into richer physical state features—**not** evidence our gas box already does so.

### 6.5 Kudithipudi et al. (2025): *Neuromorphic computing at scale*

**[SOURCE]** Nature 637, 801–812 (2025), DOI https://doi.org/10.1038/s41586-024-08253-8. A review of physical/neuromorphic hardware and system-scale challenges: coupling memory and compute, sparse event encoding, local learning, practical energy and scaling. The review supports the relevance of such architectures, not a guarantee of their superiority.

### 6.6 Artem Kirsanov's video on predictive coding / local learning

**[SECONDARY EDUCATIONAL SOURCE]** The supplied screenshot from *The Brain’s Learning Algorithm Isn’t Backpropagation* lists predictive coding, energy formalism, neural connectivity, and local activity/weight update rules, and references Bogacz, Friston, Millidge, Whittington and others. The video motivates **local error signals and predictive control** as a research parallel, but a YouTube teaching model is not proof that biological brains implement one universal predictive-coding training rule or that our controllers do so.

**Open literature question:** Which specific predictive-coding architecture and energy/free-energy expression would map to a physically observed, partially driven, stochastic piston network? Write the state and energy variables first; do not import the vocabulary alone.

---

## 7. The actual trainable-box architecture emerging from the video discussion

### 7.1 Preserve topology: soma-like controllers, axon-like transmission *edges*

The original artistic mapping “one soma becomes one box” was revised. Preferred mapping:

```text
biological view: [soma/initial segment] ==== axonal pathway ==== [receiving soma]
conceptual view: [sensor+controller] -- [hard-disk transmission module] -- [sensor+controller]
```

A whole 15-neuron scene maps into a 3D **layout** of approximately 15 controllers and selected mechanical pathways. The individual hard-disk compartments still have 2D mechanics. Controllers do not grow out of random positions; selected edges replace existing pathways, not invent new connections. The CGI is an architectural visualization and not a validated 15-edge physics simulation.

### 7.2 Minimal four-module split–merge demonstrator

```text
INPUT → [C0] → [BOX 1] → [C1] ─┬─→ [BOX 2 / large, broad] ─┬─→ [C2] → [BOX 4] → OUTPUT
                               └─→ [BOX 3 / slim, short] ──┘
```

**Topology:** exactly four transmission modules and three controller nodes. Box 2 and Box 3 are **parallel**; they can differ in compartment count, divider mass, gas density, length or height, and the network should document any unequal time delays or response properties instead of assuming that “thicker = stronger.” The system flows left-to-right visually. Each box's **drive piston is on the left, spring-backed output wall on the right** in the film; the actual historical simulation orientation can be mirrored.

Controller actions would be: C0 actuates Box 1, C1 senses Box 1's output and triggers two downstream drives, C2 senses both outputs and triggers Box 4. **A controller provides energy from an external source to each driven piston; the network is not a passive through-flow of gas or a free amplifier.** Controller cost and delay are experimentally essential.

A four-box circuit can test timing/phase sensitivity, redundant route usefulness, state-dependent gating, and *possibly* whether a controller can learn input-dependent route selection. It does not demonstrate those capabilities just by having arrows.

### 7.3 Future structure and learning knobs

Possible learning variables—distinct algorithms and physical costs:

- **Controller-only learning:** Change digital/feedback weights and actuation timing; hold all box physics fixed. Easiest first ML baseline, but not “the box mechanically learned.”
- **Readout learning:** Keep the mechanical reservoir fixed, train an inexpensive readout. Physical reservoir computing baseline.
- **Mechanical-parameter learning:** Slowly change spring constants, divider masses or chambers' effective states. Specify realistic actuators and reconfiguration work, distinguish slow learning from fast dynamics.
- **Topology learning:** Activate, deactivate, reconnect, maintain dormant routes, or alter routing—more Physarum-like. Requires switching hardware/algorithm and material or operational costs.
- **Joint predictive/morphological control:** Decide what to preserve about the past and which future actions to keep feasible (IB plus empowerment), possibly adaptive network morphology. Distant and challenging.

---

## 8. Concrete experiments: from spring law to physical learning

### Experiment E0 — *Can we represent the Román box as a spring?*

**[NEXT-STEP CANDIDATE]** Choose a validated density and finite box. Hold divider at several equilibrated displacements; independently estimate `F(L)` and `k_T`. Measure temperature, geometry and collision-force uncertainty. Compute `k_{S,eq}` under the correct one-side/two-side convention. Compare to small slow forcing. **Pass gate:** a reversible linear response within stated tolerances in the selected regime. **Control:** a literal Hookean spring with matched `k`.

### E1 — *Does it behave like a damper, inertial resonator or a system with a memory kernel?*

Drive small-amplitude sinusoids or short pulses over a range of frequencies. Record piston force, divider positions/velocities, particle energy, mode spectra. Fit (a) Hookean response, (b) linear mass–spring–damper, (c) multi-mode passive linear transfer function, then (d) nonlinear/finite-memory model only if justified. **Test on held-out frequencies, amplitudes and times**, not on the fitting records. Quantify predictive residuals and parameter identifiability.

**Important caution:** apparent wall-oscillation damping can arise from exchange into gas degrees of freedom, finite-time phase mixing or multiple modes; a closed reversible model does not acquire a heat bath merely because a fitted `γ` is nonzero.

### E2 — *Two impulses, delay, history dependence, and “refractory” analogy*

Send two same-energy piston pulses separated by `Δt`, across a spectrum of separations, with controlled initial states. A schematic input is

\[
u(t)=A p(t)+A p(t-\Delta t),\]

where `p(t)` is a finite pulse; a delta-function is a mathematical idealization, not an implementable piston displacement. Define output `y(t)` (wall motion, spring capture or task readout).

**First control:** For a linear time-invariant system with identical reset state, the response to the pair must be the sum of shifted single-pulse responses. Test:

\[
R_{\rm NL}(t;\Delta t)=y_{12}(t;\Delta t)-[y_1(t)+y_1(t-\Delta t)] .
\]

Beware ensemble fluctuations and different total injected works. **Nonzero paired-pulse interaction**, after noise and controls, indicates nonlinearity/state dependence; it does not by itself equal a biological refractory period or computational memory. If proposing an analogue of refractoriness, preregister a temporary suppression threshold and recovery-time test.

### E3 — *Can one box perform a temporal task?*

Supply repeatable, amplitude-limited time-series inputs `u(t)` (binary pulse patterns, temporal XOR, delayed classification, NARMA-type temporal benchmarks, or short-horizon prediction). Read a small set of positions, forces, energy observables or pressure proxies. Train **only a simple readout** on a training sequence, then evaluate held-out inputs and randomized seeds. Compare with linear spring–damper, linear multi-mode resonator, shuffled-state baseline, and digital equivalent at matched precision/latency. This is the first serious **physical reservoir computing** candidate.

**Metrics:** accuracy/prediction error, linear memory capacity, nonlinear memory capacity, fading memory, reproducibility, bandwidth, response delay, input work, sensing/readout cost and per-task energy. *Do not call a chaotic scatter plot computation.*

### E4 — *Does the complicated box outperform its best simple mechanical surrogate?*

Place strict caps on number of observed channels, controller computation, input work, runtime, size/energy of reconfiguration. Ask whether the multi-divider gas substrate achieves a Pareto improvement over calibrated spring–damper and multi-mode mechanical controls. If not, report the **negative result**: it reveals what microscopic complexity does and does not buy.

### E5 — *Four-box split–merge circuit with no claimed learning yet*

Implement the simple C0–B1–C1–(B2||B3)–C2–B4 graph in a reproducible coarse model using independently characterized edge impulse responses, then compare with fully coupled EDMD as computationally feasible. Preserve actuation/sensing accounting at each controller. Test parallel-route reliability/latency under changing input spectra, controlled damage and environmental uncertainty.

**A key distinction:** If the two branches are *controllers re-actuating separate physical devices*, they do not passively merge their mechanical energy into a shared chamber. C2 measures responses and uses power to drive Box 4. A separate purely mechanical junction is a different topology/model and needs distinct conservation laws.

### E6 — *Can a small box network learn?*

Start with readout or controller parameters only. Advance to trainable spring stiffness or selected geometry after demonstrating a differentiable and identifiable surrogate—or use gradient-free procedures with transparent optimization budget. For in-situ local learning, explicitly establish what physical measurements and auxiliary drives yield gradients, rather than copying Li–Mao's static identity into EDMD. Compare training sample efficiency, performance, energy and robustness to damaged edges.

### E7 — *Can routes grow, prune, remain dormant, and adapt à la Physarum?*

Introduce dynamic demand/sources, material constraints and interventions. Compare fixed graph, flow-only Physarum/Tero, predictive-memory, empowerment-oriented and hybrid policies. Count edge creation, switching, maintenance and lost capacity. This is an *extension of the transport/computing network*, not a direct claim that live slime mold uses hard-disk gas channels.

---

## 9. The information-theory bridge: learning what physical state is worth preserving

### 9.1 What is the memory variable?

Define accessible history `X_{\le t}` from measured pressures, walls, velocities, energies, local structural observables or external demands. Define compressed internal memory `Z_t`. Define target `Y_{>t}` independently (e.g. upcoming demand, next piston-input class, task reward, future useful spring response). Predictive IB asks what can be discarded from past data while preserving information about the future task:

\[
\mathcal L_{\rm PIB}=I(Z_t;X_{\le t})-\beta I(Z_t;Y_{>t}).
\]

The Markov/causal assumptions, prediction horizon and appropriate information estimator must be explicit. A controller that reads the entire simulation microstate for free is not a convincing physically efficient solution.

### 9.2 Energy cost of memory and prediction

In information thermodynamics, predictive versus nonpredictive memory can relate to thermodynamic inefficiency under a **specific driven-system model** (Still, Sivak, Bell & Crooks 2012). This motivates an information-aware cost function but does not imply the gas simply dissipates energy equal to a mutual information score in every protocol.

Define energy as physical work/heat and controller expenditure; measure useful output under a specified task. Candidate accounting:

\[
E_{\rm total}= E_{\rm drive}+E_{\rm sensing}+E_{\rm compute}+E_{\rm switching}+E_{\rm maintenance}+E_{\rm readout}
\]

with care about double counting, recovered mechanical work, energy at final time and reversible/closed subsystems. `E_transport` may be a decomposition of drive work/irreversible losses rather than an independent additive term if it would double-count energy.

### 9.3 Empowerment is different from IB

Empowerment evaluates capacity for agent action sequences to influence perceivable future states; a standard finite-horizon conditional form is

\[
\mathcal E_h(s_t)=\max_{p(A_{t:t+h-1})} I(A_{t:t+h-1};S_{t+h}\mid s_t).
\]

A route that is underused today may preserve future control after damage or source changes; maintaining it costs resources. **IB asks which history matters for a future target; empowerment asks which actions could still matter.** Neither equals shortest path, generic neural learning or thermodynamic efficiency.

### 9.4 One precise overarching optimization question

**[HYPOTHESIS]** Adaptive physical networks might trade off predictive capacity and future controllability against energy and material costs. A rigorous experiment should favor a **constrained formulation** before arbitrary weighted sums:

\[
\min_{\pi,\,\mathcal G} \;\mathbb E[E_{\rm total}]\quad
\text{subject to}\quad\text{held-out task error}\le\varepsilon,\quad\text{latency}\le\tau,\quad\text{robustness}\ge r.
\]

The control policy `π` and graph/physical design `𝒢` may be fixed or adaptive. Add an IB or empowerment regularizer only when its variable definitions and physical role are concrete.

---

## 10. Physarum: the complementary *morphology-learning* inspiration

### 10.1 Two distinct models, not one generic glowing network

**[SOURCE]** Jeff Jones's agent/trail Physarum model has small moving agents that sense, turn and deposit a diffusion/decay trail. Its occupancy/collision and causal trail update rules matter for results. A fading trace of hard-disk particle centers is only an **observational visualization**, not a chemical trail that steers subsequent hard-disk motion.

**[SOURCE]** Tero et al.'s Physarum-inspired graph model adapts edge conductances with flow, with an often used abstract structure

\[
Q_{ij}=\frac{D_{ij}}{L_{ij}}(p_i-p_j), \qquad \dot D_{ij}=f(|Q_{ij}|)-D_{ij}
\]

(up to parameter/scaling conventions). Reinforcing high-flow routes is a model of adaptive routing, not proof of a minimal total thermodynamic cost. The genuine Tokyo railway Physarum experiment motivates trade-offs among length, transport and fault tolerance; it is not the same as an arbitrarily glowing CGI blob.

The living plasmodium is one multinucleate cell; independent visual simulation agents are not literal Physarum cells.

### 10.2 The bridge to our box network

- Treat a characterized hard-disk transmission module as one **costly physical graph edge** with an empirically measured transfer function, latency, input-work requirement and controllable parameters.
- Treat Physarum adaptation as a *candidate network-routing/reconfiguration policy* that decides which edges to maintain, alter, deactivate or add.
- The mechanical edge is **not a tube of slime mold cytoplasm**; constructing/using an edge must have a realistic physical interface and a credible cost model.
- Introduce demand shifts and hazards; test whether flow-only, predictive-IB and empowerment-inspired policies retain different structures and achieve different Pareto trade-offs.
- Compare fixed/static routes, Tero, a nonadaptive equal-cost graph, randomly changing graphs and fully informed optimized planners before claiming an adaptive advantage.

### 10.3 Parallel visual phenomena are not identical mechanisms

A hard-disk center-trail density map may resemble Jones's pattern images; both can show cellular-looking voids, collective patches or stripes. But true particle trajectories in passive EDMD do not deposit a future-steering chemoattractant trail. In dense hard disks, use `ψ6`, defects, time correlation and packing fraction to study structure. In Physarum models, trace the feedback from signal to agent/flow to remodeling. The mathematical connection (if any) must be derived or tested, not inferred from resemblance.

### 10.4 Specific Physarum references and follow-up leads

- Tero et al. (2010), *Rules for Biologically Inspired Adaptive Network Design*, Science 327, 439–442, DOI https://doi.org/10.1126/science.1177894.
- Jones (2010), *Characteristics of Pattern Formation and Evolution in Approximations of Physarum Transport Networks*, Artificial Life 16, 127–153, DOI https://doi.org/10.1162/artl.2010.16.2.16202.
- Jones (2015), *The Emergence and Dynamical Evolution of Complex Transport Networks from Simple Low-Level Behaviours*, https://arxiv.org/abs/1503.06579.
- Schick, Kramar & Alim (2024), *Dynamic Cost Allocation Allows Network-Forming Forager to Switch Between Search Strategies*, PRX Life 2, 033005, DOI https://doi.org/10.1103/PRXLife.2.033005.
- Sage Jenson visualization: https://cargocollective.com/sagejenson/physarum. Useful visually; explicitly distinguish aesthetic modifications from Jones model fidelity.

### 10.5 Where is intelligence physically encoded? — *Slime Lapse* (new central research thread, 2026-10-12)

**User's prompt, preserved:** “where is the intelligence encoded in nature” (highlighted at ~04:04 of the video **`Lens of Time- Slime Lapse - bioGraphic(1).mp4`**). The user was struck by the time-lapse of an organism that has no centralized neuronal brain, yet rearranges a distributed transport system. **[SOURCE LEAD / OBSERVATION]** The film motivates the question; the images are not a measurement of information, memory, fitness or thermodynamic optimization, and biological claims require checking the underlying research.

**Working research question:** *What physical variables, material structures, and local response rules encode task-useful information in adaptive systems, and under what environmental statistics and energetic costs is investing in additional physical complexity worthwhile?*

The word *encoded* has at least four distinct senses that must not be conflated:

1. **Instantaneous physical state (`x_t`).** In Physarum, concentrations, pressures, oscillations, flow and local contractile activity depend on present and recent inputs; in HardDisks, particle positions/velocities, wall displacements and gas density carry a dynamical state. Not all state information is useful memory.
2. **Slow morphology/parameters (`z_t`).** Tube widths, connectivity, conductance, material investment, and geometric organization can persist after a stimulus changes. In a proposed mechanical network these could be divider masses, spring stiffnesses, lengths, branching topology or controller couplings. Such a state can be a *candidate physical memory*; persistence alone does not establish learning.
3. **The local adaptation mechanism (`\mathcal{R}`).** Biochemical feedback, cytoplasmic pumping and mechanochemical rules in life; controller policies, sensors and actuators in a proposed device. These determine how state and topology change. Biochemical machinery and biological heredity are not equivalent to the proposed EDMD model; passive elastic hard disks do not grow connections.
4. **The environment–task relationship (`e_t`, `Y_future`).** A physical pattern can only be called predictive or useful relative to a defined future target, action and environment. A structure that is useful for one future distribution may be costly or maladaptive for another. Intelligence is not guaranteed to reside in one place or in a single snapshot.

A minimal **hypothesis family**, not a derivation or verified biological objective, is

\[
Z_t=g(	ext{past observations, actions, morphology}),\qquad
	ext{test }I(Z_t;Y_{t+\Delta}\mid	ext{current inputs/constraints})>0,
\]

with task performance independently measured and the **complete physical cost** audited. Crucially, a high mutual-information estimate could arise from correlations or confounders; intervene on state/architecture and compare to size-matched, history-shuffled and adaptation-disabled controls.

**When does growing one more route pay?** For a specified task distribution and evaluation horizon, compare an adaptive system to matched fixed alternatives using

\[
\Delta V_{
m task,energy}\;\stackrel{?}{>}\;
E_{
m construction}+E_{
m reconfiguration}+E_{
m maintenance}+E_{
m sensing}+E_{
m actuation}+E_{
m readout}.
\]

Here `ΔV_task,energy` must be defined as a *value in compatible energy units*, or the better formulation is to **minimize total joules subject to a fixed task-quality/latency/robustness requirement**. Information Bottleneck can test whether morphology retains *future-relevant* information; empowerment can test whether spare connections preserve valuable *future action possibilities*. Neither quantity by itself establishes energetic benefit. More complexity is worthwhile only for a specified task, environment and cost model.

**Proposed disambiguating experiment ladder [OPEN/PROPOSAL]:**

- **A. Is it just mechanics?** Characterize a single Román gas compartment as an effective quasi-static spring, then benchmark dynamic models (spring, spring–damper, multi-mode mechanical response) against held-out wall responses.
- **B. Is there history dependence beyond the best linear mechanical baseline?** Give the same module paired piston impulses with varied delays and controlled equal initial states; test for reproducible nonlinear temporal features. A delayed response is not automatically a refractory period or computation.
- **C. Is any such state computationally useful?** Fix a simple temporal prediction/classification target; build a readout, train only on training data, test on held-out inputs, compare to matched linear filters and spring networks at comparable input/readout budgets.
- **D. Is morphology information or only correlation?** Use a Physarum-inspired adaptive graph with controlled training environments, then freeze, scramble or ablate its network geometry while preserving material budget; test whether its structure carries predictive information and causal task benefit beyond a fixed graph.
- **E. Is complexity worth the energy?** Allow the mechanical 4-box model to tune edge parameters, then conditionally add/remove links under changing demands. Log work input, stored energy, recoverable output, control and training cost, maintenance and reconfiguration; compare quality/energy Pareto frontiers. An EDMD system with ideal elastic collisions redistributes energy; **loss of coherent signal does not itself imply heat lost to an external bath**.

**Important physics/biology boundary:** Living *Physarum* is a single multinucleate cell with coupled biochemical, fluid, and structural dynamics; a time-lapse cannot reveal an internal learning algorithm on its own. A mechanical-neural network of tunable springs can implement physical learning in published work (Li & Mao 2024), but that does not establish adaptive learning in our hard-disk chambers. The attractive common abstraction is **physical state + local dynamics + slower trainable morphology + resource constraints**, to be assessed experimentally.

**High-level candidate contribution:** Discover *when physical memory belongs in instantaneous dynamics, fixed mechanical parameters or grown network topology*, and under which task/uncertainty conditions the extra material and maintenance are justified. A negative result favoring a simpler spring–damper would be informative. The goal is not to assert that neural, Physarum and hard-disk mechanisms are the same; it is to find measured, transferable design criteria.

---

## 11. Possible hypotheses in a falsifiable ledger

| ID | Hypothesis to test (not a result) | Decisive test | What would weaken/refute it? |
|---|---|---|---|
| H0 | A small Román gas chamber is approximable by a measured linear spring in a restricted quasi-static regime. | Independently measure `F(L)`, `k_T`, small-displacement response. | No reproducible linear regime or inconsistent conventions. |
| H1 | Multi-divider gas dynamics add computationally useful temporal features beyond a matched spring–damper. | Held-out temporal task; equal channels, input work, bandwidth and readout; matched baselines. | Complexity worsens or does not improve the Pareto frontier. |
| H2 | Two-pulse nonlinear interactions support a controllable short-term memory timescale. | Compare paired response with linear-superposition null and repeatable recovery curves. | All response follows calibrated linear multi-mode model. |
| H3 | Characterized physical transfer functions compose into predictive multi-box circuit models. | Predict C0–C2–B4 network outputs using independently measured modules. | Large unexplained residuals / nonstationary uncontrolled coupling. |
| H4 | Adaptive controller policies improve useful signal transmission per total energy. | Compare matched open-loop, linear feedback, compressed predictive control. | Sensing/actuation/compute costs erase gain. |
| H5 | Adjustable physical edge parameters permit robust or trainable network responses. | Predeclared learning task, held-out tests, damage recovery, counted training/reconfiguration work. | Only offline training works or cost dominates. |
| H6 | Morphology encodes task-relevant predictive memory. | Held-out `I(Z;Yfuture)` with size, demand, topology controls. | Predictive signal vanishes under interventions/confounds. |
| H7 | Empowerment-favored redundant paths are worth their cost under changing uncertainty. | Matched budget, stochastic demand, damage, future-control test. | Added maintenance energy exceeds useful benefit. |
| H8 | Dense-fluid internal relaxation partly controls frequency response near a liquid–hexatic regime. | Independently measured structure vs mechanical spectra and controlled modes. | FFT/mode/finite-size confounds explain effect. |

**Novelty caution:** Many individual pieces are established research fields. A plausible contribution is **integrating independently characterized microscopic mechanical channels with rigorously benchmarked learning/adaptation, prediction and a complete energetic ledger**. Literature novelty is **unverified** until scoped reviews and replication.

---

## 12. Film as a research-vision tool, not scientific evidence

Working title: **FROM PARTICLES TO PHYSICAL INTELLIGENCE** — *Energy · Information · Structure · Adaptation*.

### Descent: complex living structures to elementary physics

```text
living tree / roots
      ↓
Physarum growing network (real Tokyo experiment + cinematic transition)
      ↓
realistic neuronal morphology / axonal signal / synapse
      ↓
15-neuron network → 3D layout of controller nodes and transmission modules
      ↓
4-module horizontal 2D split–merge architecture
      ↓
zoom into BOX 1 (whole sealed mechanical apparatus)
      ↓
simplify from multi-divider drive/spring chamber to Román two-gas divider
      ↓
ACTUAL EQUATIONS / EDMD FOOTAGE / SCIENCE
```

**Ascent:** characterized equilibrium/dynamic mechanics → energy-transfer experiment → predictive control/IB → coupled adaptive networks → Physarum-like morphology/empowerment → speculative physical computing. The descending four-/fifteen-module systems are cinematic conceptual systems **not validated networks**. The later scientific ascent should explicitly distinguish real EDMD data from proposed simulations and future engineering.

### Correct geometry and non-negotiable rendering rules

1. A **neuron soma transforms to a controller**, not to a gas chamber; selected axons transform into **entire gas-transmission modules** along existing edges.
2. 4-module graph has exactly `C0 -> B1 -> C1 -> (B2 || B3) -> C2 -> B4`, with B2 large and wide, B3 thinner/smaller, and different divider configurations.
3. Box contains **its own spring, gas, dividers and piston**; there are no free-floating external springs between separate boxes.
4. A module is a 2D hard-disk gas shown as an elongated transparent rectangular channel; **identical-size hard disks, distributed through the full accessible height, gravity-free**. No unequal-size blue marbles, liquid level or leaking particles.
5. Solid **impermeable** movable dividers; a perforated grille is wrong unless explicitly modelling permeability (which current apparatus does not).
6. Film left-to-right orientation has driven piston at LEFT, spring-backed output wall RIGHT; original research screenshot showed a mirrored orientation. Be clear which artifact is being depicted.
7. The input-controller–gas module–output-controller chain uses sensed output and re-actuation with external power; do not depict particle mass passing through controllers or amplify mechanical power at no energetic cost.
8. **Exact last frame** of approved clip is the next clip's first frame; don't feed a two-panel zoom *illustration* as a starting frame. If zooming, peripheral objects exit the field of view, not dissolve/explode; preserve sealed box integrity.
9. Use image references of **exact intended final architecture** before video generation to conserve Flow/Veo credits. Generated CGI is mood/architecture reference, not automatically validated apparatus physics.
10. Anatomy: no hollow electrical neon laser inside axon, no direct soma fusion, no claim that thick axons equal strong synapses. Synaptic plasticity develops over slower timescales than millisecond spikes.

**Current cinematic development state as discussed:** tree → Physarum → neuron network, conceptual 15-module mechanical transformation accepted artistically; four-module horizontal target images iteratively corrected; a later close-up zoom attempt failed exact first-frame continuity and produced particle leaks. The next production milestone is a scientifically consistent **four-box exact-frame → one-box zoom**, then the Román transition. Assets exist in the conversation's uploaded files; this Markdown does not itself embed or duplicate large movie files.

---

## 13. Provenance, paper index, and references

### Primary sources specifically uploaded in the current conversation

1. **Li, Shuaifeng & Mao, Xiaoming (2024).** “Training all-mechanical neural networks for task learning through in situ backpropagation.” *Nature Communications* **15**, 10528. https://doi.org/10.1038/s41467-024-54849-z. Local filename: `Training all-mechanical neural networks for task learning through in situ backpropagation.pdf`. **Role:** exact local gradients and stiffness-tunable mechanical networks; note quasi-static spring model and limits.
2. **Wright, Logan G., et al. (2022).** “Deep physical neural networks trained with backpropagation.” *Nature* **601**, 549–555. https://doi.org/10.1038/s41586-021-04223-6. Local filename: `Deep physical neural networks trained with backpropagation.pdf`. **Role:** physics-aware training of controllable physical layers.
3. **Stepney, Susan (2024).** “Physical reservoir computing: a tutorial.” *Natural Computing* **23**, 665–685. https://doi.org/10.1007/s11047-024-09997-y. Local filename: `Physical reservoir computing- a tutorial.pdf`. **Role:** distinction between physical dynamics and validated computation, substrate/input/output/benchmark design.
4. **Du, Chao, et al. (2017).** “Reservoir computing using dynamic memristors for temporal information processing.” *Nature Communications* **8**, 2204. https://doi.org/10.1038/s41467-017-02337-y. Local filename: `Reservoir computing using dynamic memristors for temporal information processing.pdf`. **Role:** experimentally exploited short-term physical memory and time-series tasks.
5. **Kudithipudi, Dhireesha, et al. (2025).** “Neuromorphic computing at scale.” *Nature* **637**, 801–812. https://doi.org/10.1038/s41586-024-08253-8. Local filename: `Neuromorphic computing at scale.pdf`. **Role:** architectures, energy/scale/local learning and engineering constraints.

### Core physics and information anchors

6. **Román et al. (2002).** “The speed of sound in a hard disk gas: A computer simulation.” *American Journal of Physics* **70**, 847. https://doi.org/10.1119/1.1482060. **Role:** two-gas thermal piston and sound-speed determination.
7. **Bernard & Krauth (2011).** “Two-Step Melting in Two Dimensions: First-Order Liquid-Hexatic Transition.” *Physical Review Letters* **107**, 155704. https://doi.org/10.1103/PhysRevLett.107.155704.
8. **Engel et al. (2013).** “Hard-Disk Equation of State: First-Order Liquid-Hexatic Transition in Two Dimensions with Three Simulation Methods.” *Physical Review E* **87**, 042134. https://doi.org/10.1103/PhysRevE.87.042134.
9. **Still, Sivak, Bell & Crooks (2012).** “Thermodynamics of Prediction.” *Physical Review Letters* **109**, 120604. https://doi.org/10.1103/PhysRevLett.109.120604. **Role:** link between nonpredictive memory and dissipation in a specified framework; not a generic inequality for arbitrary boxes.
10. **Sivak & Crooks (2012).** “Thermodynamic Metrics and Optimal Paths.” *Physical Review Letters* **108**, 190602. https://doi.org/10.1103/PhysRevLett.108.190602. **Role:** near-equilibrium thermodynamic friction and finite-time control, with assumptions.
11. **Tero et al. (2010).** “Rules for Biologically Inspired Adaptive Network Design.” *Science* **327**, 439–442. https://doi.org/10.1126/science.1177894.
12. **Jeff Jones (2010/2015)** and **Schick, Kramar & Alim (2024)**—see section 10 and the earlier IdeaLabs brief for complete roles and fidelity cautions.

### Educational video/source leads (need primary-source tracing for research claims)

- `Neurons & Synapses.mp4` — structural/neural-process visual comparison, **not** evidence of particle dynamics.
- `Lens of Time- Slime Lapse - bioGraphic(1).mp4` — Physarum growth/transport visual reference; at approximately **04:04** the user flagged “Where is the intelligence encoded in nature?” as a central motivating research question. The video is not itself evidence that a network is thermodynamically optimal or that its morphology carries predictive information; see §10.5 and check its cited experiments before academic reuse.
- `Why Different Neuron Parts Learn Differently-.mp4` — apical/basal dendrite learning-mechanism discussion; identify the exact cited experiment for a manuscript.
- `The Brain’s Learning Algorithm Isn’t Backpropagation.mp4` — predictive-coding exposition and bibliography (Bogacz, Friston, Millidge, Whittington, etc.); use papers for precise model claims.
- `Optogenetics- How Light Became a Switch for Neurons - Nobel Prize 2026.mp4` — signal/control illustration lead, not a substitute for a reviewed biological source.
- `Screenshot 2026-10-11 at 21.37.14.png` — actual HardDisks multiple-gas/wall/spring simulation visual anchor, stronger authority for module geometry than CGI.
- Saved Flow-generated 4-module and 1-box images and clips — concept visuals only, with known frequent physical errors (unequal disks, gravity, grilles, loose springs, continuity changes).

### Continuity with older research note

- `IdeaLabs_Research_Brief_Physarum_Information_Energy_Morphological_Computing_2026-10-08.md`, **not superseded**. It contains detailed Jones/Tero models, Physarum repositories, particle-center-trail analysis plan, precise energetic hypotheses, SciPy/EDMD correctness concerns, the IdeaLabs/NeurAIpil separation and the earlier hypothesis ledger H1–H6.
- HardDisks/hspist3 repo under `01_improvements_bugsfxed_dev/`, especially its development-memory index and validation/FFT notes, for authoritative actual results. **Do not infer results from this idea document.**

---

## 14. IdeaLabs graph structure and project boundaries

**IdeaLabs** is the research-knowledge capture, source/evidence and idea-graph project. It is *not* the HardDisks simulator, not NeurAIpil, and not DV8. **NeurAIpil** is a separate architecture-deliberation pipeline and can eventually consume a frozen evidence packet to plan an independent Physarum or box-network prototype. **HardDisks/hspist3** remains the source of scientific truth for EDMD and mechanical experiments. Keep a possible **Physarum simulator** as a separate code project rather than merging it into IdeaLabs.

Suggested idea nodes and graph edges:

```text
Román-Static-Stiffness --measures--> Gas-Effective-Spring
Román-Dynamic-Response --tests--> Frequency-Dependent-Mechanical-Response
Gas-Effective-Spring --baseline-for--> HardDisks-MultiDivider-Transmission
SpringDamper-Network --baseline-for--> HardDisks-Physical-Reservoir
HardDisks-Particle-Traces --diagnostic-for--> Microstructure-And-Caging
Two-Pulse-Protocol --tests--> Physical-History-Dependence
LiMao-2024 --prior-work-for--> Mechanical-Local-Learning
Wright-2022 --prior-work-for--> Physics-Aware-Training
Stepney-2024 --defines-benchmark-for--> OneBox-Reservoir-Test
Du-2017 --analogy-to--> Material-Memory-Computing
Paper2-Energy-Ledger --required-for--> Energy-Efficient-Controller-Claim
OneBox-Characterization --supports-model-of--> FourBox-SplitMerge-Circuit
Predictive-IB --candidate-objective-for--> Controller-Memory
Empowerment --candidate-objective-for--> Redundant-Route-Preservation
Physarum-Tero --candidate-policy-for--> Adaptive-Graph-Routing
Physarum-Jones --separate-causal-mechanism-from--> HardDisk-Center-Trails
FourBox-Network --future-prototype-of--> Adaptive-Physical-Computing
IdeaLabs-Research-Scout --provides-provenance-to--> NeurAIpil-Research-Benchmark
```

Important source/evidence semantics for the graph: `supports`, `motivates`, `analogous-to`, `tests`, `contradicts`, `requires`, and `not-established-by` should stay distinct. Never import all the above as one giant validated idea.

### Future research-scout requests

- Identify the exact mechanical-spring physical-network video where visible weights appeared in the interval −1 to +1 and verify whether it is actual spring stiffness or an **effective signed transfer weight**.
- Compare Li–Mao's static in-situ gradient scheme to physics-aware training, reservoir-only readouts, coupled learning, equilibrium propagation, finite-dimensional mechanical waves, and differentiable event-driven simulation.
- Review experiment designs and code for mechanical reservoirs based on springs, tensegrity, soft bodies and acoustic waveguides; retrieve measured energy numbers rather than qualitative efficiency claims.
- Survey dynamical hard-disk transport and kinetic memory effects under impulsive wall forcing, including rigorous finite-size/time controls.
- Verify Physarum morphology-as-information/memory literature, resilience and construction-energy costs; identify whether predictive-information or empowerment hypotheses already exist.
- Define realistic physical actuators to adjust spring stiffness, active connections, geometry, or controller gain in an eventual hardware model.

---

## 15. Next scientifically conservative steps

**Do not replace the user's ongoing Paper-1/Paper-2 validation with a new huge hardware project.** Capture this as a longer-term comparative-computation track.

1. **Complete the Román/spring-wall validation gates** in HardDisks; stabilize sound, stiffness, conservation, reproducibility and estimator robustness.
2. **Write a small independent modeling note:** derive the linear small-oscillation Román effective-spring model and compare to a spring–damper fit, using **held-out** dynamic data. This directly answers the user's “Román = spring?” insight.
3. **Perform the two-pulse/superposition diagnostic** in one stable box; test response-history dependence with energy and state controls.
4. **Prototype a single-box reservoir task** with explicit readout and a simple mechanical baseline. Decide whether any useful advantage exists before adding a 4-box network.
5. **Only then:** simulate the C0–B1–C1–(B2||B3)–C2–B4 circuit with transparent work/sensing/actuation ledger and compare to matched passive/active controls.
6. **Separately:** continue the faithful Physarum Jones/Tero tracks and source discovery as possible *network adaptation policies*, not as alleged gas microphysics.

### What would count as a substantive result?

- A measured region where the multi-particle box is equivalent to a simple spring **and** a reproducible region where a better-specified dynamical model is needed.
- A falsifiable temporal task that the hard-disk box or network performs with clear benchmark gains (or an informative negative result).
- A full physical cost ledger per successful task—not an apparent gain obtained by ignoring powered actuation, training, or readout.
- A demonstrated link between microscopic measurements and macroscopic network performance without hidden refitting.
- An adaptive-route result that survives material, sensor and uncertainty controls.

## 16. One-paragraph research north star

The long-term project is to investigate how **energy, information, physical memory and morphology** can cooperate in adaptive systems. A Román hard-disk piston begins as a rigorously measurable mechanical response—perhaps approximable by a spring. Adding dividers, driving and observing it produces a candidate physical transmission/reservoir element. Putting these elements between sensing/actuating controllers creates a proposed mechanical information network. Spring-learning papers show that physical materials can be trained; predictive Information Bottleneck asks what information a controller or morphology must retain; empowerment asks which future actions an architecture can still influence; Physarum suggests adaptive resource allocation, growth and pruning. The research is to **test whether and when** these connections yield useful, reliable and energetically defensible behavior, **not to presume** that resemblance between neurons, slime mold and bouncing disks is itself a law of learning.

---


---

## 17. New main idea: adaptive morphogenesis — *learning when to grow, maintain or abandon a connection*

**Origin:** Follow-up research conversation on 2026-10-12, with `fogleman/dlaf`, the user's neuronal growth-cone/time-lapse videos, Physarum growth experiments, the uploaded *PhyChip* chapter, and a user-provided screenshot showing two neuronal processes approaching. **Status: [RESEARCH VISION / EXPERIMENT DESIGN], not a demonstrated general theory of biological learning.**

### 17.1 Preserve the user's original questions and intuitions

> "thats also why i wanted a physarum softawre adabpatble to test maybe different strategies."
>
> "maybe even include neurons and neuroplascticity instead sugar peices as food we grow a neural netowrk using these princliples"
>
> "a neuron might be random walker hard disk similar or sthm and strethc at random some connecitons into the network but looses them again until it finally connects for real"
>
> "every growing network actually which trasnports limited things either electiral signals , particles, water, bloood , whatever i think uses this trafde off of indormation, thermodynmaics and empowerment"
>
> "they approach from two sides, fluctuate and somehow sttistically find each other. looks very very similar to the connectors at the hairtials of navi and animals where they can connect to eeach other mutliple little connectors form a connection at some point when they kind of sense a neighbour."

**Interpretation without losing the question:** A network might grow many *candidate* branches or contacts, retain a subset based on usefulness and local cues, and withdraw others. Physarum-like flow feedback, neurite growth-cone exploration, stochastic DLA attachment, and Particle Life collective patterns are **distinct mechanisms** that can be compared using a common representation of candidate connections, their survival, and their costs. The *Na'vi neural-queue connector* reference is an imaginative/visual analogy to multiple filopodial contacts, **not** a claim about real animal neural anatomy or actual Avatar biology.

The user's long-term interest is broader than reproducing Physarum: explain *under which environmental uncertainty and energetic constraints* a physical substrate should invest in **complexity, redundancy, branching, modularity, and future options**. Ask whether generic trade-offs recur across different physical substrates, not whether trees, neurons, Physarum, and vasculature secretly implement identical learning rules.

### 17.2 Core three-stage distinction

```text
 (I) PROPOSE / EXPLORE          (II) TEST / USE                   (III) COMMIT / REVISE
 candidate process growth   ->   encounter + transport      ->    stabilization / reinforcement
 random / guided search          signal, chemical or flow         or weakening / retraction
 local geometric rules           measured functional benefit      subject to material/energy limits
                                                                       |
                                                                       v
                                                             changed possible futures
```

This is a **cross-model abstraction**. It does *not* imply that a DLA particle evaluates utility, that Physarum consciously optimizes empowerment, or that a growth cone runs an explicit information bottleneck algorithm.

Distinguish (a) **morphogenesis/development**: extending a new neurite/tube, (b) **contact selection/synaptogenesis**: forming and stabilizing a particular contact, (c) **functional plasticity**: changing transmission efficacy at an existing synapse, and (d) **topological plasticity**: creating, maintaining, or removing edges. They happen on different timescales and have different biological mechanisms. Adult learning is not reducible to neurons wandering until their cell bodies connect.

### 17.3 A deeper hypothesis that we actually can falsify

**[H-ADAPT-MORPH]** Under finite resource budgets, a locally adaptive growth/retraction policy may produce morphologies that retain more *task-relevant predictive information* or *future controllability* per unit total cost than purely random, current-flow-only, or fixed-network baselines—**for some environmental regimes**.

This needs independent measurement of (1) present performance, (2) future-demand distribution, (3) physical/material cost, (4) memory or structural encoding, and (5) causal value of the structure under intervention. Failure under fair controls is a useful negative result.

---

## 18. Branch creation versus branch adaptation: DLA, Physarum, growth cones, Particle Life

### 18.1 `fogleman/dlaf`: a stochastic branch-creation baseline

Source lead: <https://github.com/fogleman/dlaf>. **[SOURCE LEAD; recheck exact implementation options/license against a pinned commit before coding]**. DLA, *diffusion-limited aggregation*, places mobile stochastic walkers around a seed/aggregate; encounters lead to probabilistic adhesion and branching. If each added element is connected to one parent, the resulting attachment genealogy is a tree. DLA does not by itself require a learned objective, a nutrient transport model, feedback control, or rewiring. **Graph edges representing a DLA attachment history need not have hydraulic or mechanical transmission dynamics.**

Candidate experiment: hold total added material and seed geometry fixed; change attachment/rejection/guidance rules; track branching ratios, depth, loop count (typically zero for a pure single-parent tree), spatial coverage, time to target, and post-growth transport performance. For loops/redundant paths, explicitly add a reconnection/bridging mechanism; do not claim a tree gives Physarum-style loop redundancy automatically.

### 18.2 Physarum: exploration, oscillatory flow, conductance reinforcement, pruning

Keep distinct baselines:

- **Jones-style agents + causal trails:** local sensing, motion, trail deposition, diffusion/decay; growing path morphology emerges. Each simulated agent is **not** a literal cell of the multinucleate plasmodium.
- **Tero-style transport graph:** reinforcement/weakening of conductance according to flow; edges can be less active or removed according to an explicit rule, but the model does not by itself solve arbitrary spatial **birth of new edges**.
- **Actual living Physarum:** a giant multinucleate cell with mechanochemical contractions, cytoplasmic shuttle streaming, sensing, and geometry change. Oscillations may coordinate flow and growth but are **not identical to maximizing empowerment** or evidence of deliberate assessment of future branches.

Published source in uploaded papers: Adamatzky, Jones, Mayne et al., **“PhyChip: Growing Computers with Slime Mould”** (2018 book chapter, DOI: <https://doi.org/10.1007/978-3-319-65826-1_8>). Describes experimental living-substrate computing and hybrid Physarum–electronic interfaces; specifically mentions protoplasmic electrical signals, conductive coatings/nanoparticles, prototype interfaces, and why accurate computer models are needed because experiments are slow. **Evidence of prototypes is not evidence that general-purpose fast, self-growing Physarum computers already work or outperform silicon**. The chapter notes speed, longevity, scalability, and reliability limitations.

**Useful special link:** PhyChip §8.4.3 argues for Jeff Jones's simple-component model as a way to explore emergent adaptation and computation—very close to our desire for a configurable spatial simulator. It also explicitly distinguishes emergent computations performed by the substrate from metaphorical algorithms merely inspired by it.

### 18.3 Neuronal processes: directed exploratory growth rather than hard-disk random walks

In developing neurons, **growth cones** are dynamic structures at the ends of growing axons/dendrites. Fine **filopodia** extend and retract; broad **lamellipodia** and rearrangement of actin/microtubules influence direction and elongation. Adhesion, molecular guidance (attractive and repulsive), extracellular matrix, local activity, and cell identity can bias growth. **This can be stochastic and exploratory, but not purely blind Brownian random walking of whole neuron cell bodies.**

Multiple tentative contacts can be formed, tested and destabilized. **Touching another neurite is not sufficient to assert that a functional synapse formed**; pre- and postsynaptic specialization, appropriate target identity and physiological function require independent evidence. Differentiate developmental neurite guidance, transient filopodial contacts, synapse formation, strengthening/weakening of existing synapses, and neuronal network assembly.

### 18.4 Particle Life and simple-rule pattern formation

The user already has or has experimented with a Particle Life-style mode in their simulation. Attraction/repulsion and interaction rules can cause clustering, filaments, dynamic boundaries or apparently complex morphology. **An aesthetically striking emergent structure is not necessarily a graph with reusable edges, a memory system, or a learning algorithm.** If extracting a graph, state explicit segmentation/threshold/track-persistence rules and compare against a randomized control. Avoid confusing Hamiltonian elastic hard disks with attractor-driven self-propelled or actively forced particles; these have different energy budgets.

### 18.5 Comparison matrix for the proposed research

| Substrate/model | How candidate structures appear | Why a route persists | What is measured for a task | Central scientific caveat |
|---|---|---|---|---|
| DLA | Diffusing walkers attach to aggregate | Persistent deposition/irreversible attachment in simple model | Geometry, material cost, reach, optional post-hoc transport | Not inherently flow-adaptive, predictive or loop-forming |
| Jones agent Physarum | Agent movement deposits and reinforces trail field | Trail-mediated positive feedback plus decay | Resource discovery, spatial patterns, transport proxy | Abstract agents and trail concentration are not hard-disk collisions |
| Tero Physarum graph | Usually starts with graph edges; explicit add-edge needed for growth | Flow-dependent conductance adaptation | Flow/path efficiency, resilience, network volume | Conductance update alone doesn't produce unrestricted branch birth |
| Neuronal growth-cone model | Guided, noisy protrusions of neurites | Adhesion/guidance, identity, contact stabilization and activity | Target specificity, viable contact, survival and function | Genuine neuronal biology is molecularly complex |
| Particle Life | Collective self-organization under chosen local forces | Sustained nonlinear attraction/repulsion | Persistence, patterns, dynamics; graph only if defined | Pattern is not proof of useful memory |
| HardDisks modules | Engineered enclosure and actuators (no spontaneous neurite-like growth) | Explicit controller/architecture design | Transfer response, work, memory, prediction | EDMD is not a causal microscopic model of neural or Physarum growth |

---

## 19. Micro-observation: two growth fronts approaching and fluctuating

**User-provided screenshot, 2026-10-12:** `neurons from two sides.png`, taken from **“FINDING THAT CONNECTION© — neurons connecting to one another in a Petri dish — growth cones”**. On-screen attribution reads **Dr. Lila Landowski / @Rockatscientist**; video/image attribution is preserved for provenance, not a reuse license. The adjacent images and other growth/Physarum films supplied by the user are reference material, not experimental measurements made by this notebook.

![User-provided microscopy screenshot showing thin neuronal processes and filopodia-like protrusions approaching each other](assets/neurons_from_two_sides_2026-10-12.png)

**Directly visible:** narrow neuronal processes and finer protrusions extend into the space between adjacent structures. The image is visually compatible with approaching growth fronts or an encounter, but the still alone does **not** verify causal chemotactic attraction, receptor-mediated recognition, successful fusion, formation of a synapse, or an enduring functional connection. A sustained time-lapse, biological marker(s), and functional assay would be needed to test those claims.

**User's core observation:** local exploratory "little fingers" from *both sides* can appear to fluctuate, approach, and eventually establish some type of stable connection; many apparent attempts fail or retract. The resemblance to fantasy braided tendril connections is a vivid **visual analogy only**. A promising abstract geometry here is a **two-sided encounter process**: two growing tips can each explore, shorten the gap, contact, test compatibility, and either stabilize or withdraw.

### 19.1 How to turn the visual idea into a falsifiable minimal model

Keep the soma/node positions fixed. Represent a neurite tip `i` with position `r_i(t)` and direction `n_i(t)`; represent each filopodium as a finite-length transient exploratory protrusion. Allow direction to change through **both noise and defined local cues**, not assume search is perfectly random:

```text
fixed neuron/node bodies
  -> extend neurites / growth cones from specified start locations
  -> generate ephemeral filopodial candidates toward multiple directions
  -> sample local geometry and optional guidance signal
  -> extend / turn / retract under a resource budget
  -> if a candidate contacts a compatible target, establish provisional edge
  -> stabilize only after a defined contact-duration/activity/adhesion test
  -> else retract; optionally recycle material into new attempts
  -> after edge birth, separately permit functional strengthening, weakening, or loss
```

A minimal *hypothesis-generating* spatial step can be written:

\[
\Delta\mathbf r_i = v_0\,\hat{\mathbf n}_i\,\Delta t
+\chi\,\nabla c(\mathbf r_i,t)\,\Delta t
+\sqrt{2D\Delta t}\,\boldsymbol\xi_i,
\]

where the first term is directed extension, the second is an **optional assumed local guidance field**, and the third is noise. This equation is **an illustrative model, not a validated growth-cone law**. `D` here is effective stochastic exploratory diffusivity, not a hard-disk particle diffusion coefficient by default. Interaction/steric collision mechanics and persistent tip orientation would require further modeling. It is also acceptable to begin with a fully discrete lattice/stochastic model.

**Provisional connection:** make candidate edge `e` exist after contact; mark its state `s_e∈{candidate,provisional,stable,retracting}`. Introduce a stabilization probability depending on cumulative contact time and optional local activity/compatibility, with *all* dependence explicitly selected by the modeler rather than presented as discovered biology.

### 19.2 Key tests for the two-sided encounter hypothesis

- Measure **first-encounter time distribution**, failed-contact fraction, successful provisional contacts, survival time and distance-to-target. Compare *both tips search* versus *one tip searches* with matched total exploratory material/effort; the two-sided version is not automatically more efficient.
- Compare pure isotropic random exploration to persistence, attractive/repulsive gradient following, local contact adhesion, and activity-dependent stabilization. Preserve matched candidate count/total protrusion length.
- Test whether apparent successful contact reflects geometry alone: compare to nulls with shuffled target positions, frozen-tip trajectories and altered density.
- Distinguish formation of *a physical contact* from *a functional synapse*: a circuit test needs measured transfer or a stipulated artificial proxy; don't infer from a still frame.
- Account for exploration/maintenance work, reused material, time-to-connect, reliability, and a delay/penalty for late connection. A naive algorithm whose "neurons" search everywhere for free will not test our thermodynamic conjecture.
- Give each edge its own age/history and future task value; test whether **repeated successful use** predicts survival more than chance after controlling contact duration, distance and thickness.

A main candidate graph process is:

\[
\text{two-sided exploration} \rightarrow \text{provisional contact}
\rightarrow \big\{\text{stabilize, retract}\big\}
\rightarrow \text{functional adaptation}.
\]

This is a precise step beyond merely rendering random branching: it defines observable events and controls.

### 19.3 Relevant time scales

Avoid compressing into one artificial process: subsecond or second-scale membrane voltage and spike propagation are distinct from cytoskeletal tip dynamics (seconds/minutes), neurite elongation/synapse formation (minutes/hours/days), and long-term structural network changes. Physarum contractile oscillations and network reshaping have their own timescales. A simulator should have **separate clocks or rate constants** for fast signal dynamics, exploratory growth, contact stabilization, structural reinforcement, and slow pruning.

---

## 20. Design for a single comparative Adaptive Morphogenesis Lab (proposed separate project)

**Research purpose:** explore *when growing complexity is worth its energetic cost* across models with different physical rules, without falsely equating their microphysics. This extends the previously discussed **Physarum software** concept, not replacing it with a generic neural-network drawing package.

### 20.1 Four layers: keep causal physics and appearance separate

```text
  ENVIRONMENT / TASK GENERATOR
  sources, targets, nutrients/demand, uncertainty, obstacles, shocks, damage
                         |
                         v
  GROWTH / TOPOLOGY ENGINE (swappable, comparable rules)
  DLA | Jones-like Physarum | Tero + birth-rule | neurite growth cones | Particle Life
                         |
                         v
  GRAPH + PHYSICAL EDGE SUBSTRATE
  nodes + candidate/provisional/stable edges + material inventories
  ideal transport | spring-damper | measured HardDisks response surrogate
                         |
                         v
  CONTROLLER / OBSERVER / LEARNING
  fixed rules | flow-only | predictive IB | empowerment | task/energy-aware policy
                         |
                         v
  SCIENTIFIC MEASUREMENT AND RENDERING
  events, costs, conditional information, transfer, held-out tasks
  deterministic simulation state | optional real-time 2D/3D renderer
```

**Important:** Some models do not naturally produce permanent graph edges (Particle Life) and some require explicit edge-birth operators (Tero-style graph). A common *measurement interface* is desirable, but should not force incompatible systems into a false identical causal model.

A faithful **2D Jones baseline** and **2D Tero graph** should exist before an impressive 3D presentation. 3D rendering of a 2D simulator is not equivalent to solving three-dimensional cytoplasmic dynamics, neurite biophysics, or 3D hard-sphere mechanics.

### 20.2 Swappable growth and retention policies

Define policies that do not all collapse into one magical “intelligence” rule:

- **No-growth** fixed graph / pre-engineered optimum for known tasks.
- **Random DLA-like attachment** with geometry/material controls.
- **Current-demand / flow-only** reinforcement and decay.
- **Local candidate-contact rule** inspired by filopodia, with stochastic proposal and stabilization.
- **Prediction-informed** edge survival, learning from a *causal* history of observations to estimate task-relevant future demands; use explicit information-bottleneck representation/estimation.
- **Empowerment-informed** retention of optional future routes under a defined action set, horizon, and sensor model.
- **Joint cost-sensitive planning** under a fixed material and operational budget, including construction/reconfiguration costs.

Crucially, compare **predictable shifting demands**, **unpredictable switching**, **stable demands**, **adversarial/damage events**, and **many-source/many-sink competition**. Pre-register the demand process, train/test split and metric. Future unpredictability changes whether predictive memory is useful; severe surprise/damage can make spare connections valuable.

### 20.3 Energy-cost comparison without fake “thermodynamics”

Let `G_t` be the network, `u_t` external signals/requests, `Y_future` a preregistered task variable, and `Z_t` controller plus physical/morphological state. Track at least:

\[
C_{\mathrm{total}}
=C_{\mathrm{growth}}+C_{\mathrm{maintenance}}+C_{\mathrm{transport}}
+C_{\mathrm{sensing}}+C_{\mathrm{actuation}}+C_{\mathrm{computation}}
+C_{\mathrm{reconfiguration}}.
\]

In early abstract spatial-growth prototypes, these are **declared cost proxies**, *not measured joules* or metabolic ATP. Once a HardDisks edge model is attached, mechanical work/energy must be defined from the actual driven walls, spring response and energy ledger. Comparing physical hardware would require measuring controller electronics, pumps/actuators, fabrication/maintenance, and training budgets too.

One operational edge-creation criterion over horizon `H` is:

\[
\Delta J_e
=\mathbb E[C_{\mathrm{task+ops}}(G_t;H)
-C_{\mathrm{task+ops}}(G_t+e;H)\mid \mathcal I_t]
-C_{\mathrm{build}}(e)-C_{\mathrm{maint}}(e;H)-C_{\mathrm{reconfigure}}(e).
\]

Add edge when `ΔJ_e>0` **if all terms have comparable units and constraints and the inferred benefits are trustworthy**. Task penalties must be converted to the stated value/energy units or the comparison must be expressed as **minimize measured energy subject to a task-quality threshold**. This is a *testable decision rule* and not a claim that Physarum, neurons, or trees consciously compute it.

**IB:** test `I(Z_t;Y_future)` on held-out future demands, with `I(Z_t;X_past)` or state/description-cost constraints. Measure whether morphology itself retains future-relevant information, controlling network size, simple geometry, leakage and easy correlations.

**Empowerment:** for an explicitly defined action horizon, estimate potential influence on reachable/perceivable future state. It may favor unused but recoverable routes even when they waste maintenance energy; test whether that reserve option has value under the chosen future distribution. Do **not** presume that a biological oscillation or protrusion equals empowerment maximization.

### 20.4 Plugging in HardDisks without breaking the validated engine

Start with ideal or low-order edge dynamics. Calibrate an approximate input/output transfer operator from **separate measured one-box simulations** (e.g. impulse response or fitted dynamical state-space), with its valid frequency, density, geometry and amplitude regime explicitly specified. Then optionally replace edges with on-demand or offline validated EDMD experiments. Preserve energy/work and uncertainty provenance, and test held-out driving regimes.

For the minimal **four-module split–merge test**, keep topology:

`INPUT → C0 → BOX 1 → C1 → (BOX 2 || BOX 3) → C2 → BOX 4 → OUTPUT`.

`C0,C1,C2` represent controllers, and **entire HardDisks transmission modules are edges**. Each has a driven input piston, gas compartments separated by impermeable movable walls, and a spring-backed terminal output wall. Models can represent different lengths, widths, particle densities, divider masses and spring stiffnesses. Controllers act through sensors and actuators whose external work is **counted**. Passive spring–damper networks are the first baseline.

The measured effective stiffness and mode dynamics of the Román apparatus are the simplest physical calibration. The growing topology problem can be tested at the graph level **before** a complicated gas is simulated in every branch. Never make the EDMD particles literally imitate actin or a neuron growth cone without a separately validated biological mapping.

### 20.5 Experiment ladder with firm stage gates

**M0 — Geometry baseline.** Reproduce DLA-style branching; record ancestry/tree structure and open-space growth. Compare to simple random attachment and controls at equal material use. Ensure loop count is correctly computed.

**M1 — Physarum baselines.** Independently reproduce faithful Jones-like trails and Tero flow adaptation, with known conservation and model differences. Show growth/decay and adaptive path selection, including response to changed sources/targets.

**M2 — Growth-cone-inspired contact model.** Fixed soma locations, 2D candidate protrusions, contact/retract/stabilize event states, two-sided versus one-sided search, cue ablations, matched resource budgets. Validate encounter distributions and stop conditions; do **not** claim biological synaptogenesis.

**M3 — Shared environment benchmark.** Changing source/target distributions and branch breakage; fairness constraints on material, time, sensors, initial layout. Measure task success, latency, reliability, branch life histories, and full proxy costs.

**M4 — Information objectives.** Test whether predicted future demand, empowerment, or a combination improves performance **on held-out environmental shifts**, controlling for extra policy sophistication and compute cost. Negative results preserved.

**M5 — Substitute physical links.** Same resulting topology under ideal edges, spring–damper, and calibrated HardDisks transmission surrogate; identify if physics adds meaningful task capability and whether its additional cost is worthwhile.

**M6 — Beyond fixed controllers.** Introduce local update rules and topology learning, hybrid Physarum-inspired/electromechanical experiments, physical learning or reservoir tasks *only after* earlier gates succeed.

This should remain a **separate experimental/research implementation**, not an unauthorized rewrite of existing HardDisks work. IdeaLabs tracks sources/hypotheses and can propose a source-scout packet; NeurAIpil may separately deliberate on architecture; DV8 remains a separate project.

---

## 21. Working hypothesis ledger and sources for this new thread

### New falsifiable questions

| ID | Question | A first test / possible disconfirmation |
|---|---|---|
| `H-GROW-1` | Can two-sided exploratory neurite-like growth find compatible targets at lower cost or delay than one-sided exploration? | Match total protrusion length/effort; compare contact/survival distributions; advantage disappears after fair normalization = no effect. |
| `H-GROW-2` | Does stabilizing only activity/flow-relevant provisional contacts outperform indiscriminate stabilization? | Same candidate generation, budget and task sequence; randomized matching and blind held-out tests. |
| `H-GROW-3` | Can morphology encode predictive future-demand information beyond network size and geometry trivialities? | Predict on held-out demand from structural state, condition on confounds and shuffle history; absence of remaining information rejects hypothesis. |
| `H-GROW-4` | When does spare branching improve future controllability enough to offset maintenance cost? | Stable vs volatile tasks, damage, matched budgets, defined finite-horizon empowerment and utility. |
| `H-GROW-5` | Do growth/flow/retention rules yield similar task–cost Pareto regimes across DLA, Physarum and a neurite-inspired model? | Normalize only matched observables; different mechanisms may yield incomparable regimes; no presumption of universality. |
| `H-GROW-6` | Can a simple spring–damper baseline explain the same adaptive transmission task as a complex hard-disk edge? | Match bandwidth/latency/readout/work/physical maintenance; if baseline matches or beats it, complexity is not warranted for that task. |
| `H-GROW-7` | Can a *local* growth/reinforcement rule approximate a global cost-benefit optimum? | Compare local cue/credit-assignment policies to oracle planner with equal sensing/compute budgets; evaluate regret and out-of-distribution shifts. |

### Source and asset register

1. **`fogleman/dlaf`** — <https://github.com/fogleman/dlaf>. DLA growth implementation and candidate rule-intervention baseline. Record commit/license and verify behavior before reuse.
2. **`PhyChip- Growing Computers with Slime Mould.pdf`** — Adamatzky, Jones, Mayne, Whiting, Erokhin, Schumann & Siccardi, *PhyChip: Growing Computers with Slime Mould*, chapter (2018), DOI <https://doi.org/10.1007/978-3-319-65826-1_8>. Living Physarum computational substrates, prototypes, hybrid bioelectric/organic interfaces, and the importance/limits of modeling.
3. **`FINDING THAT CONNECTION© - neurons connecting to one another in a Petri dish - growth cones(1).mp4`** — user-supplied video; experimental microscopy inspiration and candidate encounter analysis. Original video attributions require verification before any public re-use.
4. **`neurons from two sides.png`** — user screenshot of the same growth-cone clip, archived in this notebook's `assets/` folder; still image supports a morphology description but does not alone prove a synapse formed.
5. **`Watch a Neuron Connect and Grow- The Fascinating Journey of Brain Synapses 🌟🧠 #scienceexplained(1).mp4`**, **`Neuron time lapse video(1).mp4`**, **`Neuron time lapse video(1)(1).mp4`** — user video references for neurite growth, culture network emergence and possible synaptic/contact dynamics. Check originals and captions for actual experimental provenance.
6. **`The Physarum Experiments, Study No.024- Interspecies Encounter(1).mp4`** — observational Physarum/organism dynamics, not automatically a quantitative model of learning.
7. **`Neurons & Synapses.mp4`** — biological morphology/visual analogy reference; filmed or reconstructed neuronal structures are not free hard-disk particles.
8. **Jones 2010/2015 and Tero 2010**, already cataloged in the 2026-10-08 IdeaLabs brief. Distinct growth/transport models with different mechanisms.
9. **Li & Mao 2024**, *Training all-mechanical neural networks for task learning through in situ backpropagation*, <https://doi.org/10.1038/s41467-024-54849-z>. Prior work for trainable spring stiffnesses, not neuronal growth or EDMD proof.
10. **Still et al. 2012**, *Thermodynamics of Prediction*, <https://doi.org/10.1103/PhysRevLett.109.120604>. Grounds questions about predictive information and physical cost under defined assumptions; does not establish a universal growth law.

### New IdeaLabs ideas and evidence relations (for future explicit import)

```text
Adaptive-Morphogenesis-Lab --proposed-software-to-test--> Grow-vs-Keep-vs-Prune
DLA-Growth --candidate-generation-baseline-for--> Structural-Exploration
GrowthCone-Filopodia --biological-inspiration-for--> TwoSided-Contact-Proposal
PhyChip --demonstrates-limited-prototypes-for--> Physarum-Hybrid-Computing
Jones-Trail-Agents --alternative-morphology-rule-to--> Neurite-Growth
Tero-Conductance-Adaptation --flow-reinforcement-baseline-for--> Edge-Retention
IB-Predictive-Memory --hypothesis-for--> Candidate-Edge-Retention
Empowerment --hypothesis-for--> Future-Option-Value-Of-Spare-Edges
Energy-Accounting --constraint-on--> Growing-And-Retaining-Complexity
HardDisks-Paper1 --independent-calibration-for--> Physical-Edge-Response
HardDisks-Paper2 --independent-calibration-for--> Useful-Transmission-Work
SpringDamper-Network --mandatory-simpler-control-for--> HardDisks-Physical-Computing
Visual-Resemblance --not-evidence-of--> Identical-Microscopic-Mechanisms
```

**Decision not made:** No existing PhD experiment is paused. No three-dimensional growing-neural/Physarum simulator is declared implemented. The starting candidate is a **narrow 2D comparative growth and pruning benchmark** whose outputs can later be visualized in 3D or connected to calibrated HardDisks transmission models.

**Research north star (preserve wording):** *When is it worth spending energy to grow, retain or remove a physical connection—and can the resulting network morphology store useful predictions or future possibilities?* This gives experimentally distinct roles to DLA (candidate branch generation), neuron-inspired contact rules (transient exploration/stabilization), Physarum (flow-mediated network adaptation), IB (predictive memory), empowerment (potential future influence), and thermodynamics (physical cost and constraints). Their integration is **a proposal to test, not a demonstrated shared biological law**.

---

## 22. New geometric hypothesis: networks from living chains of growing disks — *radii as trainable structural state*

**Date:** 2026-10-12. **Status:** [EXPERIMENT PROPOSAL / HYPOTHESIS]. This section captures the new proposal from the conversation. It is **not** an implemented simulator, a modification of the validated HardDisks EDMD code, or a claim that neurons and *Physarum* grow by hard-disk dynamics.

### 22.1 Preserve the idea in the user's own terms

> “shall we use e.g one of our boxes of hard discs without divider and jsut still circles for it? becuase at some point we could let them wiggle right and we have physics in it already.”
>
> “make start particles larger and let them maybe grow smaller hard discs in directions with one tocuh poiutn of the old hard disc so we have chains of hard discs which can then growht and we can draw a channel by drawing surrounding lines as interpolation between the top and bottom points of two discs next to each other that it seems to be a channel.”
>
> “and like in a soma the further away from soma the smaller the hard disc ? in a physarum we weould jsut use more and more in each drection to show the whole surface covered at first before strhengentig the connection.”
>
> “this could also be used as a nice model to cound costs if we give groeth and radius some energetical constraingt.”

The pivotal refinement, enthusiastically endorsed by the user:

> “Instead, make radius an adjustable state variable: r_i(t). An exploratory tip may initially be thin. A frequently used branch could become thicker. An abandoned branch could shrink. … The system doesn't merely decide where to grow. It can also redistribute limited material among competing pathways.”

**Core proposed mechanism:** A small number of large structural disks (neuron-inspired *somas* or growth seeds) spawn attached smaller structural disks. Their chains and branches explore space under local, potentially stochastic rules. Interpolating a continuous envelope around their centers and radii produces smooth processes/tubes. The radii **remain adjustable with time**, as does topology. One shared geometric substrate can visualize a neuron-inspired tree of processes or a *Physarum*-inspired expanding sheet and selective transport-tube reinforcement, provided their actual growth/adaptation laws remain distinct.

### 22.2 Separate three different “particles” that should never be conflated

| Entity | Interpretation and dynamics | State variables | What it is **not** |
|---|---|---|---|
| **EDMD gas disk** | Free, hard, fixed-radius disk; ballistic flight interrupted by elastic collisions in the validated mechanical gas apparatus. | Position, velocity, fixed radius, mass. | Not a growing actin filament, neurite segment or slime-mould cell. |
| **Morphogenesis structural disk** | Attached geometry element; may grow or retract, change radius, move with a modeled process, and spawn children. | Center `x_i(t)`, radius `r_i(t)`, parent/branch ID, attachment state, material budget. | Not a literal thermal hard disk or proof of energy-conserving microscopic dynamics. |
| **Growth-front exploratory agent** | Optional independent walker / growth cone / Jones-like Physarum sensor that proposes candidate branches or deposits a causal field. | Orientation, heading, cue observations, speed, local state, RNG seed. | Not necessarily a structural tube element or a biological cell. |

Start with a **separate 2D morphogenesis kernel** that can reuse collision/geometry/spatial-indexing concepts from HardDisks without changing the EDMD scheduler, the force/energy ledger, or the physical meaning of gas particles. Structural disks may experience explicitly modeled elastic/overdamped mechanical forces in a later extension; “wiggling” by noise is a *kinematic growth rule* until such forces and dissipation are defined.

### 22.3 Minimal circle-chain geometry

Given a parent structural disk with center `x_i`, radius `r_i`, choose a unit growth direction `d_hat` and child radius `r_{i+1}`. The initial **tangent-contact placement** is

\[
\boxed{\mathbf{x}_{i+1}=\mathbf{x}_{i}+(r_i+r_{i+1})\hat{\mathbf d}.}
\]

This provides a simple collision-free parent–child contact **at placement**, subject to detecting and avoiding collisions with other disk elements. It is not sufficient to guarantee a permanent connection if the disks move or change radius. The implementation must define an explicit link/joint, possible bond stretch, and a procedure for revising geometry when radii change.

One proposed per-element record:

```text
StructuralDisk:
    id, parent_id, branch_id, species_or_growth_mode
    position: (x, y)                 # later optional (x, y, z)
    radius: r(t)                     # dynamic, trainable/geometrical
    material_mass_or_area           # well-defined units
    tip_state: active | paused | stabilized | retracting
    attachment_state: free | candidate | confirmed | deleted
    timestamp_created, age, local_activity, last_contact
```

A minimal abstract tip trajectory may combine orientation persistence, guidance by a cue field `c(x,t)`, and random exploration:

\[
 d\mathbf x_{\mathrm{tip}}=v\hat{\mathbf d}\,dt+\chi\nabla c\,dt+\sqrt{2D}\,d\mathbf W_t.
\]

This is a **phenomenological growth-path rule**. A more anatomically defensible neurite model distinguishes growth-cone advance, extension/retraction of filopodia, contact formation, axon guidance, dendritic elaboration and activity-dependent synapse stabilization. It should not imply that entire neuronal somas randomly walk toward each other.

### 22.4 A touching chain is a skeleton, not a sealed tube

Two tangent disks touch at one mathematical point, so their union does **not** contain an open passage through which fluid can travel. The simulator must distinguish:

1. **Connectivity graph:** explicit parent–child/branch/candidate-contact edges; graph connectivity does not depend on pixel overlap.
2. **Geometric skeleton:** disk centers and a smooth interpolated centerline `c(s)`.
3. **Rendered morphology:** tube/ribbon boundaries generated from a smoothly varying radius `r(s)`; in 2D use `c(s) ± r(s)n_hat(s)` along the curve; in 3D use an actual swept cross-section or signed-distance/implicit surface.
4. **Transport law:** a *separate* physical or graph-level rule for hydraulic, electrical, or mechanical signal transfer along the resulting edge. A visually smooth channel is not evidence that it has been simulated as a continuum fluid or an axon.

At sharp turns and bifurcations, naïve top/bottom interpolation can self-intersect or leave gaps. Prefer a regularized centerline, well-defined junction blending, explicit self-contact policies, and quality tests; rendering should never alter growth dynamics or measurements. If the 2D model uses an occupied-area field for sheet-like growth, treat it as a separate state from its thin tube skeleton.

### 22.5 No hard-coded monotonic taper — `r_i(t)` is a physical memory variable

**Key choice:** Do **not** force every new child disk to be smaller than its parent. A decreasing initial radius may be a useful starting condition for neurite-like taper, but *retention and reinforcement* must be capable of changing radii dynamically.

\[
\boxed{r_i=r_i(t),\qquad r_i(t)\ge r_{\min}>0\quad\text{while the element exists}.}
\]

Illustrative, deliberately **not biologically derived**, feedback law:

\[
\dot r_i=\alpha\,f(a_i,\widehat v_i,\text{local cues})-\beta\,r_i-\lambda_i,
\]

where `a_i` is local activity/transport, `v_hat_i` an optional estimated future-use signal, `α` reinforcement gain, `β` decay and `λ_i` a budget/projection pressure term or local competition. Pure flow-based reinforcement uses only measured activity; predictive or empowerment-guided experiments introduce different terms **as separate policies** and must pay for obtaining those signals. A geometric implementation may instead update the **branch-level** radius `r_e(s,t)` and derive individual disk radii consistently.

Candidate developmental states:

```text
new thin exploratory process
  -> stochastic extension and branching
  -> possible neighbor encounter (candidate contact)
  -> local tests / activity / signals
  -> stabilize, thicken or maintain       [connection retained]
  -> pause, thin or retract               [connection abandoned]
  -> reclaim material at stated fraction [if permitted by the model]
```

The current thickness of a connection now becomes a possible **material memory of its earlier use**, but this observation is not evidence that the memory predicts the future or that the policy is learning. That requires an explicitly measurable task and controlled counterfactuals.

---

## 23. Shared geometry, different biological inspirations: neuron, Physarum and Particle Life

### 23.1 Neuron-inspired mode — soma, neurites, growth-cone contact

- Start with relatively large **fixed-position somas**. Some emanating processes may begin thicker near the soma; axons need not monotonically thin with distance, and synaptic strength must not be equated with axon diameter.
- Create smaller chained structural disks at active neurite tips, with heading persistence, growth-cone-inspired local protrusions, cue gradients and collision avoidance.
- Allow fine terminal filopodia to probe nearby processes: two-sided exploration, provisional contact, retreat, repeat, and possible eventual stabilization.
- Define **axon-to-dendrite or other biologically plausible contact categories** rather than treating every touching disk as a neuron-to-neuron synapse. A stable morphological contact is not proof of neurotransmission.
- Later, introduce separate activity-dependent synaptic weight changes versus slower growth/pruning. Growth-cone attachment and mature synaptic plasticity have different cellular mechanisms and timescales.

**Visual goal:** recreate the *phenomenology* of the uploaded growth-cone videos—fine protrusions meeting, fluctuating and sometimes stabilizing—without claiming to reproduce molecular axon guidance, vesicle release or validated neuroplasticity.

### 23.2 Physarum-inspired mode — dense exploration followed by route concentration

The user imagines using many growing disks initially to create an exploratory **covered area**, followed by changing radii and retracting weak pathways so a sparse transport network remains. This is a compelling two-stage model, but a valid *Physarum* simulator must distinguish the living plasmodium's continuous, multinucleate sheet-like body from independent Jones-model computational agents.

Candidate stages:

```text
seed plasmodium / initial occupied patch
  -> many locally advancing elements or moving front fill accessible territory
  -> connect them to a field/continuum-like occupied-region representation
  -> define transport on an extracted or explicit candidate tube graph
  -> reinforce active routes (e.g. local flow), shrink unused tubes
  -> retain limited loops/redundancy when constraints warrant it
```

For an established graph, a Tero-type reference is `Q_ij=(D_ij/L_ij)(p_i-p_j)` with conductances adapting to flow. This is a **separate phenomenological rule**, not something that follows automatically from geometric disk growth. One model may use radius-dependent conductance; if using an ideal cylindrical laminar hydraulic tube, `G∝r^4/L` is an **assumption-dependent benchmark**, not a law governing neuronal spikes or elastic HardDisks boxes.

### 23.3 Particle Life / DLA modes — pattern formation as controls

- **DLA:** stochastic exploratory walkers attach to a structure and generate dendritic branches. Basic DLA is a branch *creation* baseline, with no necessary flow-dependent pruning or learning.
- **Particle Life:** attraction/repulsion between types produces collective patterns. It is an alternative *local-interaction physics* baseline, not automatically a network with stable, functional edges.
- **Growth-cone mode:** directed, attached tips with transient contacts.
- **Physarum mode:** distributed exploratory frontier, causal trail/flow fields, and tube remodeling.

A single **scene renderer** may visualize all four. A single **scientific model** cannot be assumed to cover them by changing colors: event definitions, conservation laws, transport observables, and attachment/retraction dynamics must be separately specified and compared.

---

## 24. Material- and energy-limited morphology: what does it cost to grow and thicken?

### 24.1 Use one conserved or externally supplied growth resource

Start with **material** before claiming microscopic thermodynamic joules. Every structural addition, enlargement and removal updates a common ledger. Choose and document a consistent geometric approximation:

- **Independent 2D structural disks:** `M_disk = ρ_A Σ_i πr_i²`; use this only when the disks represent distinct nonoverlapping material elements. If disks are only control points for an interpolated ribbon, summing their areas can **double count** material or introduce resolution dependence.
- **Continuous 2D ribbon** around nonoverlapping, sufficiently smooth centerline: approximate area `A ≈ ∫2r(s)ds`, with junction overlap correction. Multiply by an explicit 2D mass-per-area density if useful.
- **3D circular tube:** `V ≈ ∫πr(s)²ds`, correcting junction volumes and caps; material mass `M=ρV`. Double radius at fixed length means approximately four times the volume.
- **2D occupied Physarum sheet:** measure the **union area or explicitly simulated occupancy field**, not the sum of every overlapping growth agent's area.

Growth events spend material; retraction may return a specified recoverable fraction. Somas and anchored seeds may have a distinct one-time initial mass. The budget rule is explicit:

\[
M_{\rm occupied}(t)+M_{\rm free}(t)+M_{\rm irrecoverable}(t)=M_{\rm initial}+M_{\rm supplied}(t),
\]

or use a finite external budget with a documented spending schedule. This is a **designed resource accounting rule**, not a claim about biological conservation at cellular/metabolic scales.

### 24.2 Distinguish growth material, operating energy and information costs

A useful ledger is

\[
C_{\rm total}=C_{\rm construction}+C_{\rm maintenance}+C_{\rm transmission}+C_{\rm sensing}+C_{\rm controller}+C_{\rm reconfiguration}.
\]

At first, `C_*` are **defined model proxies** or material budgets, not automatically joules, ATP or thermodynamic entropy production. Costs should be tracked per time and per branch, with an action and budget audit. If a policy prunes a tube and recovers material, recovering mass is not identical to recovering the energy invested in construction. There may be irreversibility, biosynthesis cost, dissipative work and maintenance overhead.

A simple model cost can begin with `C_build=γ_A ΔA_new` in 2D, where `ΔA_new` is newly occupied area; a 3D version can use `γ_V ΔV_new`. Costs of removal can be nonzero. If one policy has greater computational complexity or needs global information, **charge its additional sensing/inference cost** or report separate budgets so comparisons remain honest.

### 24.3 What makes thickness valuable? Test, don't assume

Thickening may increase transmission capacity or robustness but has material cost. The `r^4` hydraulic law is meaningful for a steady Newtonian fluid in an ideal circular cylindrical tube under low-Reynolds-number laminar conditions. Axonal electrical conduction, distributed chemical oscillations and the gas-divider-piston HardDisks module obey different laws. Define the edge transfer function **per substrate**.

Questions worth testing: one thick path versus several thin paths; thick primary routes plus thin dormant reserve paths; fast regrowth versus paying long-term maintenance; controllability gain from an extra loop; response after branch damage; delayed demand switches. Report quality–cost Pareto frontiers rather than assuming complexity pays for itself.

---

## 25. Minimal simulator R0: two somas, growing disks, tentative contacts

**Goal:** create a scientific instrument with an accessible visual behavior, not a fictitious whole-cell neuronal model.

**R0 state:** two fixed seed disks with radius `R_soma`; finite free material pool; structural disk chains with explicit parent links; growth tips with headings and RNG state; provisional contacts and branch radius history. Initially use a flat 2D domain and deterministic seeds. Set collision/attachment tolerance and maximum branch/tip counts to prevent numerical runaway.

**Per time step:**

1. Update the fixed environmental cue field (or no-cue control) and each active tip's heading/position using a specified growth law.
2. Propose child disks tangent to the parent; reject illegal overlap, out-of-domain growth or unaffordable additions; log failed proposals.
3. Detect **candidate** approach/contact between unlike processes at a distance threshold, possibly allowing several protrusions to compete. Do not instantly create a mature synapse or assume every collision is functional.
4. Update link/contact states: ephemeral → stabilized or retracted according to **one explicitly chosen policy** (random retention first).
5. Update `r_i(t)` and any branch-radius smoothing; apply material budget and contact geometry corrections. Enforce bounds, record reallocated and lost material.
6. Recompute graph connectivity and optional coarse signal/flow task on this structure; sample analysis metrics at a fixed physical/model time cadence. Render separately.

**R0 validity tests:** identity/parent invariants; no unpermitted overlap; parent–child bonds preserved under radius growth; bounded budget residual; deterministic replay by seed; renderer independence; consistent pruning; no orphan active tips; conservative treatment of junction area; time step robustness. Simpler contact mechanics may need constraints and relaxation rather than exact tangent placement after dynamic radius changes.

**Next controlled experiments:**

- **R1:** fixed somas with random one-sided versus two-sided exploratory tips; equal total attempted growth/material/time; measure first contact time, stable-contact probability and retraction count.
- **R2:** introduce local directional cues; ablate gradient sensing and persistence separately; measure whether they help beyond extra exploration effort.
- **R3:** add `r_i(t)` with the *same* flow/activity reinforcement rule across seeds, then compare fixed radius, activity-thickening, and random thickening at matched total material.
- **R4:** introduce time-varying demand/task, dormant routes and damage; compare local activity/Physarum-like flow to a nonadaptive and an oracle benchmark.
- **R5:** separately introduce predictive IB or empowerment-informed decisions. Require well-defined `Z_t` (morphology), future target `Y_{t+H}`, action and perception spaces, estimator uncertainty, and **held-out** environments. Do not infer learning only from network aesthetics.
- **R6:** compare identical grown topologies under ideal connections, springs/dampers and externally calibrated HardDisks edge response surrogates, accounting for sensor/actuator input work and uncertainty. No uncontrolled changes to EDMD research.

**Metrics:** contact latency, successful stable connections, graph loops/redundancy, branch lifetime, radii and material distributions, cumulative material built/pruned/recovered, maintenance cost, response/task quality, prediction accuracy or predictive information on held-out sequences, achievable finite-horizon influence (when defined), controller computation cost, and total task-cost Pareto fronts.

**Critical falsifier:** if stochastic local exploration plus material budget already matches sophisticated policies in matched novel environments, the information-theoretic complexity is not justified for that task. Preserve this as a useful negative result.

---

## 26. Interface to existing HardDisks/hspist3 and future 3D display

**Authoritative boundary:** validated EDMD fixed-radius gas remains a separate model. We do **not** turn actual thermal disks into tethered structural disks or change particle radii in the EDMD engine. Reuse drawing tools, computational geometry and indexing algorithms only if their invariants and assumptions remain valid. The growth-kernel is a separate physics/problem definition.

A future *network edge* can select one of several transport models:

```text
Morphology graph (nodes, processes, radii, branching, lifetimes)
    |  \
    |   \-- geometry/material/flow-only rules (Physarum analogue)
    |---- ideal static conductance or delay
    |---- spring–damper/linear transfer function
    +---- calibrated hard-disk piston/divider/spring input-output response
```

A mechanical HardDisks module may approximate an axon-like lossy temporal pathway **functionally**, with a separate controller triggering an input piston and measuring terminal spring response; it is not an axon and does not automatically regenerate spikes. The complex edge must beat simpler baselines on explicit tasks to warrant its simulation and eventual hardware cost.

Start with 2D scientific computation. A three-dimensional **rendered view** of a 2D network is not 3D physical simulation. A genuine 3D growth kernel needs new collision/contact algorithms, tube geometry, branch intersection rules, boundary effects and cost accounting. Allow attractive volumetric renderings only when the underlying model fidelity is visibly documented.

---

## 27. Team handoff for IdeaLabs and discussions with Astra, Fable and Opus 5.5

**Research concept:** *Trainable Morphogenetic Networks with Resource-Constrained, Radius-Adaptive Connections.* Proposed small standalone simulator, linked as a research item and evidence packet in **IdeaLabs**; not a modification of the existing HardDisks paper code. Use **NeurAIpil** only as a distinct architectural deliberation/evaluation project if desired. DV8 is unrelated to the scientific core.

**Vision:** A network of physical elements explores, attaches, reallocates material among candidate pathways, and learns what to keep under uncertain demand. Distinguish *connection birth* (DLA-inspired), *exploratory contact* (growth-cone-inspired), *reinforcement/pruning* (Physarum-inspired), *future-predictive memory* (IB), *future actionability* (empowerment), and *operating/structural cost* (thermodynamics/material budgets). All are candidate **separable modules**, not one already established biological law.

**Three questions to critique before any implementation:**

1. Can a disk-chain + interpolated-radius envelope provide **stable geometry, valid junctions and a material-conserving radius-update law**, or would a signed-distance/phase-field or graph skeleton plus meshing be scientifically cleaner?
2. What is the **smallest nontrivial held-out environment** where local radius allocation and structural pruning could beat static/random/flow-only baselines under a matched material and task budget? How do we prevent the task definition from baking in an advantage?
3. Where, if anywhere, do **IB and empowerment add testable predictive value** beyond ordinary reinforcement, robust network design and an oracle shortest-path/flow solution? Which information variables are observable, and can energy/sensing overhead be counted?

**Astra/Fable/Opus review prompts:** Ask each independently to (a) find mathematical/biological/engineering counterexamples; (b) propose at least two sharply different minimal models; (c) identify existing research, algorithms, repos and licensing/novelty obstacles; (d) write falsifiable hypotheses and experiment controls; (e) critique budget accounting and model-selection biases; and (f) give a pragmatic R0 implementation plan without interrupting Paper 1/2. Preserve points of disagreement, not only a synthesized consensus.

**Proposed IdeaLabs nodes / relations:**

```text
Radius-Adaptive-Structural-Disks --represents--> Morphological-Growth-Geometry
Structural-Disk-Chain --has --> ParentChild-ContactGraph
Interpolated-Tube-Surface --renders --> Structural-Disk-Chain
Material-Budget-Ledger --constrains --> Birth/Thickening/Pruning
Radius-History-r_i(t) --candidate-physical-memory-for --> Previous-Flow-or-Activity
GrowthCone-Exploration --proposes --> Ephemeral-Contacts
DLA --candidate-growth-baseline-for --> Branch-Birth
Tero-Physarum --flow-retention-baseline-for --> Radius-Reinforcement
IB --tests-predictive-value-of --> Radius-and-Topology-Memory
Empowerment --tests-actionability-of --> Dormant-Redundant-Routes
SpringDamper --baseline-for --> Calibrated-HardDisks-Transmission-Edge
HardDisks-EDMD --must-remain-distinct-from --> Morphogenesis-Growth-Kernel
```

**Do not claim:** validated 3D biophysics, neuron formation/synaptogenesis, Physarum cytoplasm, hard-sphere thermodynamic cost for geometric growth, an IB-optimal natural system, free actuation, or a novel universal law. The new computational physics research begins once one of the R0–R6 tasks produces a reproduced, baseline-controlled effect with clear units and confidence intervals.

**One-sentence north star for reviewers:** *When does a physically growing system benefit from reallocating limited material among exploratory connections, so its radius distribution and topology preserve useful future possibilities and predictions at a lower total task cost?*

---

### Update log

- **2026-10-12 — v1:** Consolidated current conversation and existing Physarum/HardDisks idea brief. Added Román-as-spring, spring–damper comparison, mechanical in-situ learning bibliography, axon/controller functional mapping, four-box split–merge topology, refractory/paired-pulse caution, physical reservoir test ladder, IB/empowerment, energetic accounting, film constraints and cross-project ownership. Status deliberately marked exploratory; none of the later model/network ideas is reported as validated.

- **2026-10-12 — v1.1:** Elevated the *Slime Lapse* question “Where is intelligence encoded in nature?” to a central research question, added four candidate loci of physical encoding, experimental interventions for separating memory/correlation/causation, the complexity-versus-energy decision criterion, and a dedicated video source entry (§10.5). No experimental finding newly claimed.

- **2026-10-12 — v1.2:** Added §§17–21: DLA/Physarum/neuron/Particle-Life branching and adaptation synthesis; exact user exploration questions; neurite-growth-cone/two-sided encounter screenshot and cautious observation; *PhyChip* hybrid wetware research; modular 2D→3D growth simulator proposal; energy/IB/empowerment cost criteria; M0–M6 validation ladder; H-GROW hypothesis ledger, source/provenance map and IdeaLabs graph edges. Archived the provided screenshot under `IDEAS/assets/`. Existing v1.1 content preserved.

- **2026-10-12 — v1.3:** Added §§22–27: user-proposed growing/tangent structural disk chains; channel interpolation and its point-contact limitation; dynamic/trainable radii `r_i(t)` as possible material memory; allocation, thickening, thinning, branch birth/retraction, 2D/3D material accounting and explicit cost-ledger distinctions; neuron/Physarum/Particle Life modes and model-fidelity boundaries; two-soma R0 and R1–R6 controlled experiment ladder; HardDisks EDMD isolation and physical-link comparisons; IdeaLabs graph nodes and critical questions for Astra/Fable/Opus 5.5. Prior sections retained unchanged except metadata/cross-reference.

---

## 28. [2026-10-12 / v1.4] Fable handoff — what we are *actually* trying to learn

**Status:** This is a reproducible research *proposal and memory*, not an executed experiment, established new law, or statement that a living system optimizes IB/empowerment. It complements §§0–27 and the separate October 8 IdeaLabs brief; retain the original sections and original user questions.

**User's request (near-verbatim intent):** Preserve the five-paper HardDisks research staircase originally planned with Fable and add the growing-circle/disc chain concept and comparison with DLA-style point-line graphs. Explore trainable physical networks, Physarum-like resource allocation, neuronal protrusions that reach toward neighboring cells and retract, and *when investing energy in structural complexity is useful*. Keep the real HardDisks research independent and validated, but leave an honest experimental path from characterized boxes to adaptive and perhaps biological/electrical/mechanical network hardware. Share with Fable for critique, not automatic endorsement.

**The question to ask Fable:**

> **Can a physical network learn *where*, *whether*, and *how thick* to grow connections, using local information and a limited material/energy budget, such that the resulting morphology helps predict and control future environmental demands at a lower *total* cost than a simpler fixed network?**

Four separable claims must be tested independently:

1. **Emergence:** simple local rules can generate branching networks. This is already known in DLA, Physarum-inspired models, and other growth models; reimplementing it is not novelty.
2. **Functional adaptation:** altered topology/radius actually improves *held-out task performance* over resource-matched baselines; visually attractive growth is insufficient.
3. **Predictive structure:** morphology contains information about future demand beyond trivial covariates, and causal perturbation shows that retained structure contributes to improved performance.
4. **Energy advantage:** gains persist when creation, upkeep, sensing, actuation, learning, and switching costs are included; growth resources in an abstract simulator are **not** automatically thermodynamic joules.

### 28.1 High-level research map — keep these domains distinct

```text
ESTABLISHED HARD-DISK MECHANICS
  Román equilibrium force/stiffness + wave dynamics [Paper 1]
          |
  driven-piston -> collision-mediated transfer -> spring [Paper 2]
          |
  control / predictive memory / useful work [Paper 3, proposed]
          |
NETWORK ADAPTATION [Paper 4, proposed]
  graph topology + changing coupling/thickness + quantified costs
     |              |                |
     |         growth-rule labs      +-- spring/damper baseline
     |        (DLA / Physarum /         +-- measured hard-disk edge model
     |         neuron-inspired)
     v
LATER HARDWARE DIRECTION [Paper 5 / long-term exploratory]
  3D geometry, distributed sensing and physical learning, hybrid substrates
```

**IdeaLabs** is the knowledge/hypothesis/evidence graph; **HardDisks/hspist3** is the validated physics source; **Adaptive Morphogenesis Lab** would be a separate simulator repository; **NeurAIpil** is a separate software architecture/deliberation system (not the simulation or the science itself). No need to merge these projects.

## 29. Five-paper staircase — mathematical plan and experimental boundaries (for Fable)

This is the **planned progression**, not evidence that all five papers exist or their experiments have been completed. Paper 1 is an ongoing physical characterization, Paper 2 investigates driven transfer, Papers 3–5 are increasingly speculative. The recent morphology concept should be an **extension, not a reason to skip Paper 1 and Paper 2**. Reconcile with Fable's original project documents and the authoritative current repo before deciding experiment counts or submission claims. Earlier uploaded Fable/COWORK PDFs named `260916_paper2_plan_energy_box_stepbystep_COWORK.pdf` and `260918_paper2_theory_audit_COWORK.pdf` could not be programmatically read in this update: do not pretend to reproduce unseen details from those files.

### 29.1 Paper 1 — Román two-compartment hard-disk piston: equilibrium thermodynamics, acoustic response and possible dense-fluid structural dynamics

**Actual minimal physical apparatus:** fixed outer boundaries, gas on *both* sides of a central movable divider, identical-radius hard disks moving ballistically between exact elastic collisions, no gravity. Let one side have length `L`, height `H`, `N` disks of diameter `σ`, mass `m`, temperature `T`, and central-divider mass `M`. Each side's packing fraction is

\[
\eta=\frac{N\pi\sigma^2}{4LH},\qquad n=\frac{N}{LH},\qquad P=nk_BT\,Z(\eta).
\]

`Z(η)` is the 2D hard-disk compressibility factor, not an independently fit 'learning parameter'. For a fixed wall at constant `T`, the *one-sided* force is

\[
F(L,T)=PH=\frac{Nk_BT}{L}Z(\eta).
\]

**Measured force protocol:** clamp the divider at each chosen `L` for an independent equilibrated run; accumulate *one-sided collision impulses* and their time stamps. Do not estimate one-side pressure from net divider force (the symmetric net mean may vanish). Move the clamp *only between separate equilibrium runs*. The measured slope gives

\[
\boxed{k_T=-\Big(\frac{\partial F}{\partial L}\Big)_T
=\frac{Nk_BT}{L^2}\left[Z+\eta\frac{dZ}{d\eta}\right].}
\]

For classical 2D hard disks, equilibrium internal energy `U=Nk_BT` per side, so `C_V=Nk_B`. The equilibrium adiabatic one-sided stiffness is

\[
\boxed{k_{S,\mathrm{eq}}=k_T+\frac{F^2}{Nk_BT}
=\frac{Nk_BT}{L^2}(Z+\eta Z'+Z^2).}
\]

This describes a fully equilibrated adiabatic response under the model assumptions. It is not an automatic identity for the finite-frequency, spatially inhomogeneous resonant mode. For an ideal gas (`Z=1, Z'=0`), `k_T=Nk_BT/L^2`, `k_S=2Nk_BT/L^2` (one 2D compartment), an unusually simple starting point for the future spring–damper comparison.

**Dynamic experiment:** unclamp divider; track `x_wall(t)` at a fixed simulation-time sampling cadence, enough to resolve its natural oscillations. Infer spectral peak `f_*`, linewidth (carefully defined), and cross-check with independent time-domain estimates. Román's idealized symmetric acoustic matching relation, in the stated convention, is

\[
\cot K=\frac{M}{2Nm}K,\qquad \boxed{c_{\mathrm{app}}=\frac{2\pi L f_*}{K}}.
\]

The acoustic one-side quantity under a homogeneous, long-wavelength, nondispersive approximation is

\[
\boxed{k_{S,\mathrm{ac}}=\frac{Nm c_s^2}{L^2}.}
\]

For a thermodynamically consistent equilibrium EOS, the corresponding bulk relation is

\[
\boxed{c_s^2=\frac{k_BT}{m}(Z+\eta Z'+Z^2).}
\]

**Conventions:** If a symmetric central divider is compressed on both sides, the total restoring stiffness is twice the one-side contribution. Record whether `N` refers to one side or total disks. Compare equilibrium and acoustic results with approximation/finite-size/relaxation diagnostics—not by asserting equality.

**Dense-fluid hypothesis and measurement safeguards:** The 2D liquid–hexatic coexistence window deserves testing for finite-frequency structural relaxation. Measure local sixfold order

\[
\psi_6(i,t)=\frac{1}{n_i}\sum_{j\in\mathcal N_i}e^{6\mathrm i\theta_{ij}(t)},
\]

orientational correlations, neighbor exchanges, defect trajectories, and density/flux. An illustrative—not automatically valid—relaxation model is

\[
k^*(\omega)=k_\infty-\frac{k_\infty-k_0}{1+\mathrm i\omega\tau}.
\]

Check actual constitutive sign convention and phase carefully if reporting real/imaginary parts. Finite-size peak shifts, wall loading, acoustic dispersion, estimation artifacts, and incomplete sampling can mimic a relaxation effect. **A linewidth is not directly a viscosity or a unique structural relaxation time.** Run window-function, oscillation-count, predictor, time-domain, mass-variation and equilibration robustness tests first.

**Specific plan recovered from the conversation (draft only):** An October 5 design considered total system sizes `N_total = 100, 400, 900`, seven packing fractions `{0.695, 0.700, 0.704, 0.708, 0.712, 0.716, 0.720}`, roughly nine divider masses, with height scaling `H=10 sqrt(N_total/100)` and `L` derived from η. The `N=100` apparent ~14.7% sound-speed dip and high reduced χ² in the coexistence-window mass fits were **exploratory reported figures**, not established phase physics; an `N^{-1/2}` dip extrapolation was a hypothesis for experiment sizing, **not a prediction to present as validated**. Better oscillator-length targeting is essential (earlier fixed-duration low-density upturn was suspected to involve FFT resolution). See the authoritative repo's current scripts/tables for accurate execution details; don't copy this design as final registration.

**Paper-1 go/no-go:** force/EOS/static/dynamic consistency in ordinary fluid within resolved uncertainties, sound-speed estimator robustness, honest limits near coexistence, energy/momentum/collision-integrity regression tests. A negative dense-fluid result is publishable information if the calibration is strong.

### 29.2 Paper 2 — driven work transfer through compartments with movable dividers and a terminal spring

**Actual apparatus reference** (orientation from the *original simulation screenshot*, not the film's left-to-right schematic):

```text
fixed anchor -- coil spring -- LEFT MOVABLE WALL | hard-disk gas A |
MOVABLE DIVIDER 1 | hard-disk gas B | ... more optional dividers ... |
RIGHT DRIVEN PISTON  <-- prescribed input displacement/force
```

The spring is **inside the terminal empty region**, attached to the left movable wall, not a coil randomly linking two modules. Intermediate movable walls are impermeable and partition independently colliding hard-disk gas populations. All disk radii within a monodisperse physical run are the same. Image/movie representations sometimes mirror the entire apparatus so information flows visually left-to-right; label that as a visualization convention. A pictured 'controller' and multiple module network are not validated features of the existing simulation.

A state-resolved energy ledger for a purely elastic, closed gas + movable-wall + spring system driven from outside is

\[
\boxed{W_{\mathrm{in}}(t)=\Delta K_{\mathrm{disks}}+
\Delta K_{\mathrm{movable\ walls}}+\Delta U_{\mathrm{spring}}+
\Delta U_{\mathrm{other\ conservative}}+Q_{\mathrm{out}}+R_E.}
\]

Include `Q_out` only if a thermal bath, explicitly dissipative damper, or other sink is modeled. `R_E` is a numerical residual and should approach measurement/numerical tolerance. The spring potential is

\[
U_{\mathrm{spring}}=\tfrac12 k_s(x_s-x_{s,0})^2.
\]

For a prescribed piston motion `x_p(t)`, external work requires a consistent signed force/displacement convention:

\[
W_{\mathrm{in}}(t)=\int_0^t F_{\mathrm{external\ on\ piston}}(t')\,\dot x_p(t')\,dt'.
\]

For a specified evaluation time and positive net work input, one illustrative spring-capture measure is

\[
\boxed{\epsilon_{\mathrm{capture}}(t_f)=\frac{U_s(t_f)-U_s(0)}{W_{\mathrm{in}}(t_f)}}\quad(W_{\mathrm{in}}>0).
\]

This is **not** a universal thermodynamic efficiency; the spring may re-release energy, and temporary coherent response can be highly protocol/endpoint dependent. Also report peak spring energy, energy returned to source, particle/wall kinetic energy, and delivered signal at the target.

For small linear acoustic disturbances in a channel, an illustrative characteristic force/velocity impedance is

\[
\mathcal Z=\rho c_s H=\frac{Nm c_s}{L}.
\]

One idealized point-mass divider transmission benchmark is

\[
\mathcal T(\omega)=\left[1+\left(\frac{\omega M_d}{2\mathcal Z}\right)^2\right]^{-1}.
\]

The latter does **not** predict arbitrary finite nonlinear multi-divider boxes. Validate no-divider, one-divider, multiple-divider configurations in a staged ladder. Study dependence on input waveform/amplitude/frequency, density/η, divider masses, geometry, spring constant, and output time window; compare with linear-response and quasi-static baselines.

**Paper-2 go/no-go:** collision integrity, global work-energy conservation, correct mechanical boundary conditions, independent force and wall-position diagnostics, clearly defined task-specific output, and parameter-out-of-sample predictions from Paper 1. Do not use 'energy loss' to mean internal conversion from organized piston motion into unresolved gas degrees of freedom in a purely elastic isolated system.

### 29.3 Paper 3 — predictive information / thermodynamic control of a partially observed box

**Proposed system:** sensors return coarse pressure, wall position/velocity, measured collision rates and/or acoustic indicators; they do *not* grant the controller omniscient particle microstate access. The controller forms a compressed causal memory `Z_t` from observation history `X_{\le t}` and selects future piston actuation `A_t`. A predictive IB-style family of objectives:

\[
\boxed{\mathcal L_{\mathrm{PIB}}=I(Z_t;X_{\mathrm{past}})-\beta I(Z_t;Y_{\mathrm{future}}).}
\]

Define `Y_future` **before** testing (e.g., predictable time-varying load, target response, changing environmental demand) and account for the fact that actions affect futures; naive passive predictive IB can be misleading for feedback control. A task-level control cost might be

\[
J=\mathbb E[\ell(Y_{\mathrm{target}},Y_{\mathrm{out}})]
+\lambda_E E_{\mathrm{actuation}}+\lambda_S E_{\mathrm{sensing}}
+\lambda_M E_{\mathrm{memory/compute}},
\]

where all cost terms have declared units/scales and are compared at matched task quality/latency. Information measures are in bits/nats and require multipliers to combine with joules. Literature link: Still et al. (2012) and related partly observable information-engine models; no result yet that the hard-disk controller outperforms matched open-loop or optimal linear feedback.

**Controls:** passive/no controller; same waveform open-loop; proportional/linear feedback with same sensing/actuation budget; full-state oracle upper bound (clearly labeled unattainable); deliberately non-predictive compressed memory; learned controller. Record fresh-test-sequence performance, energetic cost and uncertainty.

### 29.4 Paper 4 — adaptive physical transport *network*: parameter learning, topology and useful redundancy

**Proposed leap:** integrate several physically characterized modules as network edges; a soma-inspired controller serves as a node receiving sensor output and driving a downstream piston. First benchmark is the **four-box split–merge experiment**:

```text
INPUT -> C0 -> [BOX 1] -> C1 --+--> [BOX 2] --+--> C2 -> [BOX 4] -> OUTPUT
                               |               |
                               +--> [BOX 3] --+
```

Each controller requires a specified source of work; no passive magical amplification. This directed computational graph is a *functional control architecture*, not necessarily one giant sealed mechanical assembly. Boxes 2 and 3 may have unequal length/width/number of dividers and different calibrated response functions. In a directly mechanically coupled variant, clearly specify conservation-compatible couplings instead of conflating a sensor/actuator interface with a spring.

Stage A: learn parameters on a fixed topology (`k_s`, wall masses, geometry, actuation, coupling gain, effective radius/length in graph abstraction). Stage B: permit edge activation/deactivation and latent dormant routes. Stage C: penalize formation/maintenance and test **growth/pruning/reconfiguration** in changing tasks. Stage D: compare to a Physarum-style local flow-rule and to simple greedy shortest-path plus matched material baselines.

The separate DLA–Physarum–neuronal growth lab can generate candidate graph topologies. It should not be represented as actual hard-disk network physics until edge-response and energy models have been calibrated and independently validated.

### 29.5 Paper 5 — future physical/morphological computing substrate, potentially 3D and hybrid

**Vision, not a confirmed paper:** a growing/adaptive network whose topology, connection thickness, and internal dynamics jointly store task-relevant information. Candidate substrates include mechanically adjustable channels, soft/biomaterials, fluid/electromechanical hybrids, live Physarum sensing/conductive hybrids, nanowire/memristive networks. Geometric 3D rendering is *not* a validated 3D hard-sphere EDMD implementation. Any claim of advantages over silicon/neuromorphic chips requires matched workload, fidelity, accuracy, bandwidth, latency, thermal/power overhead, lifetime and total operating energy, not just attractive pattern videos.

**Staircase invariant:** Paper 1 characterizes physics → Paper 2 quantifies transfer → Paper 3 tests constrained prediction/control → Paper 4 tests whether structure/parameters should adapt → Paper 5 only then probes substrate/generalization. Independent simpler controls remain required at every level.

## 30. [NEW] Disk-chain structural model versus DLA-style points/lines — what is needed, what is optional?

**Key decision:** Keep a **graph for topology**, **structural disks/radii for spatial geometry**, and **a replaceable transmission/physics model for functionality**. These are three layers; one should not silently stand in for the others. `fogleman/dlaf` is an instructive fast branching baseline with parent links. We do *not* need circles to count material if every graph edge already has `L_e` and `r_e`; circles are nevertheless useful for explicitly resolving local touch, crowding, branch thickness and contact emergence.

### 30.1 Canonical structural entities and source of truth

Each structural disk is a **material/growth element**, distinct from a **thermal hard disk** of `hspist3`, and also distinct from a **DLA diffusing walker**:

\[
D_i(t)=(\mathbf x_i(t),r_i(t),p_i,\mathrm{type}_i,\mathrm{state}_i),
\]

where `p_i` is a parent or attachment identity. Graph nodes can represent somas, branches and contact junctions; graph edges can represent chains of structural elements. Retain stable IDs, ancestry, birth/death times, contact state, task-use histories and material-source accounting. A disk's *visual radius* need not equal the physically relevant hydraulic or electrical lumen radius unless a specific model defines the mapping.

For an end-to-end tangent placement, one candidate child center is

\[
\boxed{\mathbf x_{i+1}=\mathbf x_i+(r_i+r_{i+1})\hat{\mathbf d}_i.}
\]

Allow persistent tip direction with stochastic exploration, obstacle exclusion and chemotactic guidance, e.g.

\[
\Delta\mathbf x_{\mathrm tip}=v\hat{\mathbf d}\,\Delta t
+\chi\nabla c\,\Delta t+\sqrt{2D\Delta t}\,\boldsymbol\xi,
\]

where `ξ` is standard normal and `D` has units of length²/time. This is a *growth-tip displacement proposal*, **not** the Newtonian collision law for thermally moving hard disks. Step length and angle rules may need normalization and joint treatment; test numerical convergence. Two somas are mostly stationary in the minimal neurite model. Branch tips and fine protrusions fluctuate; don't simulate whole somas as ballistic gas particles.

**What changes after radius updates?** If center distances were set exactly tangent initially, changing `r_i` breaks tangency. Define how junctions relax, overlap, separate or reparameterize their centerlines. Do not let rendering interpolation silently repair invalid computational geometry. A chain of externally tangent discs has *point contacts* and by itself is **not a connected fluid-conducting lumen**; a smooth swept tube, spline plus radius field, 2D union with finite overlap, or level-set/implicit surface must be specified if continuous material is intended.

### 30.2 Graph-only representation is a proper baseline, not an inferior hack

For edge `e` with centerline arc length `L_e` and radius field `r_e(s)`, we can compute costs without material disks:

\[
A_{e,2\mathrm D}\approx\int_0^{L_e} 2r_e(s)\,ds
\quad\text{(ignores caps/overlaps)},\qquad
V_{e,3\mathrm D}\approx\int_0^{L_e}\pi r_e(s)^2\,ds.
\]

For constant radius, `V_e≈πr_e²L_e`. For multiply intersecting branches, aggregate material is **the geometric union volume/area**, not the naive sum. Discretization refinement must not change total cost for an unchanged physical shape. Disk *count* measures numerical resolution unless the model explicitly treats each disk as a real material unit of a fixed physical size. Varying radii changes the geometry and material; subdividing an unchanged link into more numerical disks should not magically consume extra joules.

### 30.3 Disk geometry is most useful for local, mechanistic events

1. **Birth:** a tip adds a candidate daughter disk with available material.
2. **Contact:** two tips, or tip and shaft, fall within a geometric range; provisional contact forms.
3. **Contact state machine:** `exploring -> touching -> provisional -> stabilized`, or `provisional -> retracting -> removed`. A single image does not prove an electrical synapse.
4. **Thickness adaptation:** `r_i(t)` changes using a specified locally measurable signal and a material ledger.
5. **Motion:** optional mechanical springs/bending/viscosity allow fluctuation while preserving topology; motion, creation and actuation must have separate costs.
6. **Pruning:** delete/retract elements while accounting for what material is reused and what energy is lost; reuse cannot erase energy previously expended.

Use fixed-time seeded simulation steps, spatial hashing or neighbor data structures, tests for overlap/collision and continuity, and geometry-vs-renderer invariance. Compare graphs and disk chains at **matched effective geometric/material resolution**, not equal object count.

### 30.4 What's the correct `r_i(t)` law?

Do **not** impose monotonic taper as a universal law. Biological dendrites often taper, many axons have more stable caliber, and Physarum tubes can thicken or regress. Let radii be *dynamic structural state*, but demand a defined causal reinforcement variable:

\[
\dot r_i=\operatorname{Proj}_{[r_{\min},r_{\max}]}
\big[\alpha\,u_i(t)-\beta\,m_i(t)-\gamma(r_i-r_{\mathrm baseline})\big].
\]

Here `u_i` is a **specified** local usefulness indicator (e.g. time-averaged flow, contact activity, prediction-relevant sensory signal), `m_i` is local material scarcity/cost, and the projection/constraint is an algorithmic boundary handling rule, **not** a biological or thermodynamic derivation. Enforce a global material pool via constrained updates or a chemical-potential/resource dual variable rather than clipping all radii independently and then claiming material conservation. The rate coefficients have units chosen to make `dr/dt` length/time.

**If thickness is useful:** test by changing one radius while holding topology and task fixed. Compare gain in flow/latency/reliability to incremental formation/upkeep cost. Two common, non-interchangeable limiting models show that *different physical tasks imply different value of radius*:

\[
G_{\mathrm{hydraulic}}=\frac{\pi r^4}{8\mu L}
\quad\text{(ideal Newtonian Poiseuille tube: assumptions required)},
\]

\[
R_{\mathrm{ohmic}}=\frac{\rho_{\mathrm{electrical}} L}{\pi r^2}
\quad\text{(uniform cylindrical conductor: assumptions required)}.
\]

A neuronal axon uses active membrane electrophysiology and myelin; neither expression is an all-purpose neuron law. Neither automatically models confined hard-disk 2D gas propagation with movable internal walls.

### 30.5 Visually distinct but scientifically compatible modes

- **Neuron-inspired mode:** large soma-like seeds, a few growth cones with persistent stochastic exploration, contact guidance, transient filopodia, local contact-stabilization and synapse-formation *proxy*. The new short video `The Tiny Connections That Make Your Brain Learn 🧠(1).mp4` is a **visual morphology reference**: small dendritic protrusions and variable process thickness visible. Its overlay *“Everything you learn looks like this”* is an evocative simplification, **not** proof from the six-second clip that a particular connection was strengthened by learning. Two-sided approach/provisional contacts reference: `IDEAS/assets/neurons_from_two_sides_2026-10-12.png`.
- **Physarum-inspired mode:** many exploratory local growth elements can approximate broad advancing occupation; reinforce selected tubes and retract other regions under flow-based feedback. For realistic sheet-to-tube development, use an occupancy/density field (or level set) in addition to the sparse graph. Real plasmodia are multinucleate single cells, not separate freely roaming neuronal cells.
- **DLA mode:** random walkers aggregate irreversibly to generate a spanning branched morphology; standard tree-growth does not include functional pruning, predictive memory or closed-loop energy optimization.
- **Particle Life mode:** multi-species attraction/repulsion yields pattern formation. It is a useful emergent-geometry control, not automatically a physically growing biological network.

## 31. Material, energy and information ledger — enough to compare growth strategies fairly

A realistic simulator needs **both** a material/resource budget and energy bookkeeping; neither is equivalent to disk count.

### 31.1 Material accounting

Let `Ω_t` be the occupied 2D or 3D structural region; `ρ_mat` be structural material density. Then

\[
M_t=\rho_{\mathrm{mat}}\,\operatorname{measure}(\Omega_t)
\quad\text{(area in 2D or volume in 3D with appropriate density units)}.
\]

A limited-material experiment imposes `M_t≤M_max`; a sourced material model uses mass flow balance

\[
\boxed{M_{t+\Delta t}-M_t=M_{\mathrm{import}}-M_{\mathrm{export}}-M_{\mathrm{unrecovered}}.}
\]

When a branch retracts, decide how much material returns to the available pool, not just delete the drawing. For a 2D simulator, energy/area is a **model proxy** unless mapped to a thickness and physically calibrated tissue/gel density and chemical potential. Evaluate material *distribution* as a state, not merely total mass.

### 31.2 Cost ledger and control problems

A suitable conceptual total cost over horizon `[0,T]` is

\[
\boxed{E_{\mathrm{total}}=
E_{\mathrm{construct}}+E_{\mathrm{maintain}}+
E_{\mathrm{operate}}+E_{\mathrm{sense}}+
E_{\mathrm{compute/readout}}+E_{\mathrm{actuate}}+
E_{\mathrm{rewire}}.}
\]

Avoid double counting energy terms; any recovered work should be represented explicitly in a signed ledger. Use `C_*` instead of `E_*` for abstract unitless model penalties. Illustrative per-area construction cost `C_construct=γ_A A_new`, per-time upkeep `C_maintain=∫ μ_A A(t)dt`. Branch retraction does not automatically refund construction work. If the lab is entirely simulated, these may be **defined proxy costs**; record units prominently. Real hardware/bath dissipation and operating joules require a separate validated model or measurement.

At matched task performance quality `q≥q_min` and latency `ℓ≤ℓ_max`, compare architectures by minimum cost:

\[
\boxed{\min_{\pi,G_t,r_t}\ \mathbb E[C_{\mathrm{total}}]\quad
\text{subject to } q\ge q_{\min},\ \ell\le\ell_{\max},\ M_t\le M_{\max}.}
\]

A cost-aware growth rule should compare expected marginal future *task value* to total marginal construction/upkeep/reconfiguration cost in consistent units; this is a **candidate policy**, not a proven biological law:

\[
\Delta J_e=\mathbb E[V(G+e)-V(G)]-C_{\mathrm{birth}}(e)-
C_{\mathrm{upkeep}}(e)-C_{\mathrm{switch}}(e).
\]

Allow **not** growing as an explicit action; otherwise the simulator trivially rewards filling space.

### 31.3 Predictive IB and empowerment are different observables or policy priors

Morphological memory can be encoded as `Z_t = encode(G_t,{r_i(t)},selected internal states)` with an explicit, consistent encoding. Pre-register a future target `Y_{t+τ}` (next active sink, spatial reward map, future damage location only if learnable, etc.). The predictive IB trade-off is

\[
\boxed{\min I(Z_t;X_{\mathrm{past}})-\beta I(Z_t;Y_{\mathrm{future}}).}
\]

**Caution:** if `Z_t` depends on actions, the causal conditioning and control problem must be formalized; don't infer biological IB optimization from statistical correlation. Estimate mutual information on held-out trajectories with permutation/size controls; report finite-sample uncertainty and checks against trivial encoding of environmental location.

For empowerment under state `s_t`, action sequence `A_{t:t+k-1}`, future perceptible state `S_{t+k}`:

\[
\boxed{\mathcal E_k(s_t)=\max_{p(a_{t:t+k-1})}
I(A_{t:t+k-1};S_{t+k}\mid s_t).}
\]

This measures an **action-to-future-observation channel capacity**, not directly the expected reward or thermodynamic efficiency. An extra branch might preserve possible future actions, but it also increases costs. Add empowerment-based policies only after the states, actions, horizon, intervention semantics and estimator are explicit. Compare directly with no-empowerment policy under the same material/sensor budgets. *Having many options is not always useful or cheap.*

### 31.4 Causal test — is morphology really carrying useful memory?

Observe a trained network after a history of changing resource demands. Create interventions: reset or shuffle widths `r_i` while holding topology `G` fixed; reset topology while preserving as much local material as possible; replay the same future input sequence; compare to matched random morphology. If disrupting **history-dependent** structure selectively reduces future performance, that is evidence of causal utility beyond appearance. Conversely, if effects disappear with mass, path length, or branch-count controls, do **not** claim information-theoretic memory.

## 32. Minimal reproducible prototype and falsification ladder (proposed)

Start **2D** with 2 fixed soma-like disks and limited material. Make only simple local tip rules; no IB/empowerment loss initially. The *graph* is the topology source of truth, *disk chains* are geometry, and a *renderer* smoothly draws the tissue. Proposed state variables: seed positions; per-element IDs/parent IDs/center/radius; tip orientation and growth state; provisional contacts; resource pool; task inputs; optional per-edge flow. Treat all random numbers as seeded and record units.

| Gate | Minimum demonstration | Necessary control / failure trigger |
|---|---|---|
| **R0 Geometry** | Single soma grows stable touching/tapered or variable-radius chains; renderer draws continuous processes | Refining disk count must not change occupied material/cost; check tangent and junction logic |
| **R1 Encounter** | Two tips stochastically explore, can touch, temporarily couple, stabilize or retract | Compare with Brownian/random walker and no-guidance controls; contact alone ≠ synapse |
| **R2 Budget** | Global pool limits total branch material; growth/thickening competes with other paths; pruning recovers a declared fraction | Verify mass ledger and that removing an edge cannot manufacture energy/material |
| **R3 Function** | A source–sink flow or signal task defined on resulting graph | Fixed graph, DLA tree, geometric shortest paths, and random branching at equal resource budget |
| **R4 Reinforcement** | A local flow/activity rule adjusts `r_i`/edge parameters and sometimes prunes | Ablate reinforcement; show transport/task benefit not just nice veins |
| **R5 Prediction** | Environment switches among learnable demand patterns; IB-inspired local/controlled memory is introduced | Held-out futures, matched sensing and material budget, shuffled-history controls |
| **R6 Controllability** | Preserve or regrow options for changing demands/damage via empowerment-inspired policy | Same task utility, cost, learning horizon and matched baseline; confidence intervals |
| **R7 Physical link** | Replace graph-edge proxy with calibrated spring/damper, then independent HardDisks-based edge response | Validate against reserved EDMD traces; include actuator/sensor energy and latency |

**Proposed first scientific figure:** side-by-side same-seed networks (i) DLA simple branching, (ii) growth cones with local contact, (iii) Physarum-flow reinforcement, (iv) energy-limited radius adaptation; plot `task_score`, `material`, `growth_cost`, `upkeep`, `latency`, `heldout_prediction`, and `dormant_route_resilience`. The image is not a result until experiments run.

**Implementation boundary:** the original EDMD engine keeps fixed particle masses/radii, exact collision scheduling and its energy-conservation validation. *Structural* discs may grow and attach only in a **separate simulation module or repository**. When a future physical branch bends through spring forces, add a separate mechanical model; when it conducts fluid, add a transport model; when it transmits information like an axon, add explicit sensors/actuation/refractory/electrophysiological assumptions or avoid biological equivalence claims.

## 33. Visual sources, hypotheses and the scientific interpretation to preserve

**New video reviewed:** `The Tiny Connections That Make Your Brain Learn 🧠(1).mp4` (~6 s). The visible content includes soma and branched processes with many fine outward protrusions and locally different diameters; a camera/view progression is shown. The overlay *“Everything you learn looks like this”* is motivational rhetoric, **not experimental proof that a specifically observed growth/contact event caused learning**. Preserve the clip as a morphological/animation reference only. Specific neuron growth-cone videos previously supplied likewise show local search/tentative contact; do not infer completed functional synapses or measured energetic selection from microscopy alone.

**Other sources already in notebook:** `FINDING THAT CONNECTION© ... growth cones(1).mp4`; `Neuron time lapse video(1).mp4`; `Lens of Time: Slime Lapse`; `Neurons & Synapses.mp4`; `PhyChip: Growing Computers with Slime Mould.pdf` (Adamatzky et al.); Tero et al. 2010 (`10.1126/science.1177894`); Jones 2010/2015; Li & Mao 2024 (`10.1038/s41467-024-54849-z`); Wright et al. 2022 (`10.1038/s41586-021-04223-6`); Stepney 2024 (`10.1007/s11047-024-09997-y`); Du et al. 2017 (`10.1038/s41467-017-02337-y`); `fogleman/dlaf` https://github.com/fogleman/dlaf. **Source existence and claimed mechanisms do not establish novelty of any combined proposal.**

**My assessment for Fable:** Starting from circles attached to circles is excellent as an *explainable and instrumentable model representation*, not a claim that neurons really grow in discrete disks. The scientifically valuable design choice is to preserve the same topology and allow model-substitution: fast point–line DLA baseline versus radius-bearing skeleton versus explicit touching/overlapping disks. We can then identify when explicit geometry buys scientific fidelity, whether changing thickness captures functional memory, whether any adaptive policy beats flow-only heuristics, and whether a real energy advantage survives comparison to optimal simpler systems. I would not train a huge 3D network yet. I would demand a two-soma, limited-budget contact/retraction benchmark first, then a changing-demand flow task, then parameter learning, and only then expensive HardDisks edge dynamics and 3D growth.

## 34. Direct handoff for Fable, Astra and Opus 5.5 — adversarial review request

> We have a verified EDMD HardDisks/hspist3 program and a staged five-paper plan (Román EOS/stiffness/acoustics → piston-to-spring energy transfer → predictive partially observed control → adaptive physical network → long-term computing substrate). We now propose a *separate* Adaptive Morphogenesis Lab using a graph plus attached structural disks of variable radius `r_i(t)` to generate neuron-like exploratory branch tips and Physarum-like sheet-to-tube network adaptation. DLA and Particle Life are non-learning structure-generation baselines. Our main scientific target is a fair test of when material/energy-expensive growth, reinforcement, alternate paths and physical memory help prediction/task performance/future control. The long-term network may replace ideal graph edges with independently calibrated spring-damper or hard-disk gas transmission modules; nodes can be explicit sensor-controller-actuator components. No experiment yet demonstrates the claimed energetic/predictive superiority. Please audit the mathematics, identify prior art, critique physically implausible assumptions, propose *the smallest decisive test*, and separate near-term Paper-1/2 work from exploratory Papers 3–5.

**Please specifically challenge:**

1. **Conceptual novelty:** Which literature already combines neurite growth, Physarum conductance reinforcement, adaptive radii, reservoir computing, physical learning, predictive routing and cost-based pruning? Cite primary sources; a visual analogy is not novelty.
2. **Cost accounting:** What would constitute honest joules, and what must initially be abstract cost proxies? How do we avoid double counting reuse, operating dissipation and biochemical/actuation costs?
3. **Coarse-graining:** When does a circle-chain geometry approximate a continuous tube? What numerical resolution gives invariant mass and contact behavior?
4. **Biological realism:** What does a cultured-neuron time-lapse actually establish? How much growth-cone guidance, adhesion and synapse stabilization must be modeled before using the word neuroplasticity?
5. **Physics:** What extra signal-processing capacity could hard-disk many-body dynamics offer relative to optimally chosen spring/damper or simple electronic delay/reservoir elements, at matched total cost?
6. **IB and empowerment:** Define target variables, observer/action horizons, causal conditioning, local implementability and estimator bias. Under what nontrivial circumstance would they predict morphology better than flow/usage/rule-based baselines?
7. **Paper-1/2 gating:** Which exact present physical results need further validation before any control or morphology claim can be grounded? The October 2026 size sweep is draft, not an approved final protocol.
8. **Minimal R0 specification:** Recommend model laws, invariants, observables and reproducibility tests for a two-soma disk-chain encounter with evolving `r_i`, plus a fair same-task graph-only control.
9. **Potential disproof:** What result would lead us to abandon the extra disk geometry, IB objective, empowerment objective or hard-disk edge apparatus for a much simpler model?
10. **Roadmap prioritization:** What belongs inside the dissertation versus post-PhD long-horizon vision? Protect rigor and achievable publication timelines.

**Decision requested from reviewers:** Return (a) strongest steelman, (b) three fundamental counterarguments, (c) missing literature, (d) precise minimal experiment, (e) falsification criteria, (f) realistic scope/effort, (g) whether the five-paper staircase requires adjustment. Avoid proposing a monolithic feature-complete 3D simulation as the first deliverable.


---

## 35. [2026-10-12 / v1.5] New Physarum fan image: punctate nuclei are **not** the growing material itself

**Primary visual source provided by user:** `Screenshot 2026-10-12 at 22.26.23.png` (archived to IdeaLabs as `assets/physarum_fan_projected_nuclei_2026-10-12.png`). The screenshot shows panel **C** labeled **Fan**, a magenta fluorescence micrograph at scale **200 μm**. The user noted that bright spots and thin structures appear similar to disks/circles in the proposed morphogenesis engine. Preserve that visual intuition **and** the actual physical identity of the observations.

**Figure provenance verified:** Mathieu Le Verge-Serandour and Karen Alim, *Physarum polycephalum: Smart Network Adaptation*, *Annual Review of Condensed Matter Physics* **15** (2024), 263–289, DOI: https://doi.org/10.1146/annurev-conmatphys-040821-115312 ; preprint https://arxiv.org/abs/2306.09063 and review PDF https://www.bpm.ph.tum.de/wp-content/uploads/2023/08/2306.09063.pdf . **Figure 2C** is identified by the original caption as a growing *fan* with **projected trajectories of SYTO 62–labeled nuclei**, adapted from prior work. It is **not** a photograph showing that the tissue grows by rigid hard-disk objects attaching at contact. The fan is an **advancing thin protoplasmic sheet** influenced by upstream contractions and cytoplasmic flow. Some fans retract; others persist, developing channels/tubes. See further primary microscopy evidence *Spatial transcriptomic and single-nucleus analysis reveals heterogeneity in a gigantic single-celled syncytium*, eLife (2022): https://elifesciences.org/articles/69745/figures (Video 3: nuclei shuttle flow across a fan region). Fluorescent nuclei can be *tracers of flow and intracellular distribution* without being geometric building blocks of the material.

**Important modeling consequence:** keep **three** physically distinct object classes even when they all appear visually as circles:

1. **Structural elements**: attached geometry proxies (structural disks with dynamic radii `r_i(t)` or a continuous area/level-set field) that determine material volume, mechanical boundaries, growth and retraction.
2. **Mobile tracers or nuclei**: moving contents *inside* living Physarum, potentially informative about streaming, pressure-driven transport and history, but not identical to structural disk primitives.
3. **Gas hard disks**: fixed-radius ballistic particles in validated `hspist3` EDMD with elastic collisions, mathematically distinct from both growth elements and nuclei.

This suggests a new test/visualizer: a sheet-like **fan mode** with spatial occupancy `φ(x,t)` and advective transport of marker particles in a velocity field, plus a graph/tube skeleton extracted behind the advancing fan. An initially dense front can become a thinner vein network only through an explicitly modeled redistribution and boundary evolution; rendering a point cloud alone is not sufficient. Track **material conservation/budget** separately from **tracer count**; tagged nuclei do not constitute the material budget. For an actual Physarum mechanism, additional actomyosin contractions, flow feedback and biochemistry would be required. The first simplified simulator is explicitly *inspired*, not a biological reproduction.

**Figure-specific question for Fable:** Does the vivid difference between distributed nuclei trajectories and adaptive tube morphology favor (i) a particle-based **tracer** field plus (ii) a level-set **structural** growth field over purely tangent structural disks for a realistic Physarum mode? Suggest a low-complexity, independently verifiable comparison and clarify which quantities are conserved.

## 36. Four newly submitted PDF sources — role, actual claim and limits

These were provided in the current conversation. Treat **source-derived claims** separately from **our proposed use**. They are not interchangeable evidence for morphology learning.

### 36.1 John Duggan (2018), *Physarum Computing*, 16-slide CS 594 teaching deck

**Local original:** `Physarum+Computing.pdf`; author **John Duggan**, CS 594: Unconventional Computation, dated November 30, 2018. It describes Physarum's response to attractants and repellents; pseudopodial/tube exploration; maze/shortest-path demonstrations; summaries of logical gates, robot controller, graph approximations and primitive memory. On page 6 it lists primary-paper leads: **Nakagaki et al. (2000)** maze, **Tsuda et al. (2004, 2007)** gates/robot controllers, **Shirakawa et al. (2009)** Delaunay-like approximation, **Saigusa et al. (2008)** memory, and several Adamatzky examples. Recheck the deck's own references for exact spelling/bibliography before quoting or publication.

**Use:** a *secondary literature map* and possible scenario/task inventory for Physarum mode: shortest-path and obstacles, multiple food targets, logic/gates, changing demand, and memory. **Limit:** a teaching slide saying 'optimized tube' or listing a solution to an NP-hard problem does not establish globally optimal performance, polynomial-time solution to arbitrary NP-hard instances, energy efficiency, predictive IB, or empowerment optimization. Prefer the cited primary experiments and their conditions.

### 36.2 Institute for the Future (2018), *Immersive Human Networks*, 48-page report

**Local original:** `Immersive2_Human_Networks_Report(1).pdf`; commissioned by Cisco, authored by IFTF's Emerging Media Lab with the organizational-network VR prototype **ONVR** developed with Slanted Theory. It documents interactive multi-user immersive **network analysis**, toggling views between organizational relations and flows, and proposes VR as a means of exploring network structure at different scales. It is an *interface/sensemaking* project, not a validated biological or physical adaptive-growth model.

**Use:** a *visualization and research-workflow* reference for the future Adaptive Morphogenesis Lab and IdeaLabs evidence graph: inspect topology at multiple scales, compare alternative policies, color edges by growth/flow/reliability/energy, explore branch history and interventions, and separate **observed graph** from **model-predicted dynamics**. Do not infer any new mathematical rule of network growth or energetic optimality from it. Check report's third-party CC BY-NC-ND licensing restrictions for any reuses.

### 36.3 Singh et al. (2025), *Transfer learning via distributed brain recordings enables reliable speech decoding*

**Local original:** `Transfer learning via distributed brain recordings enables reliable speech decoding(1).pdf`; *Nature Communications* **16**, 8749 (2025), DOI https://doi.org/10.1038/s41467-025-63825-0 . The authors use **distributed stereo-EEG (sEEG)** recordings for speech-motor tasks, a variable-length sequence/phoneme decoder and **cross-subject transfer learning** to uncover shared latent representations despite differences in coverage and neural activity; grouped pretraining outperforms individual-data-only models in their reported task/cohort. This concerns **decoding task-related activity** and transferable computational representations, *not growth cones, changing axons, biologically driven synapse formation, or morphology-based learning*.

**Potential connection** [INTERPRETATION]: information can be read out from distributed, spatially incomplete neural measurements; a fixed network's activity and a decoder's trained weights may carry useful task information even without adaptive biological topology. This is an important **competing explanation** to test in our models. Compare a fixed physical morphology with a powerful trained *readout* against a growing morphology with equally strong readout, matching total memory/sensing/training budgets. Don't interpret successful speech decoding as proof that neurons grow to implement predictive IB.

### 36.4 Brosler et al. (2026), *Simultaneous speech and gesture decoding for multimodal communication in paralysis*

**Local original:** `simulataneous speech and gesture decoding for multimodal communication in paralysis(1).pdf`; *Nature Neuroscience* (online 14 September 2026), DOI https://doi.org/10.1038/s41593-026-02446-2 . The paper reports parallel decoding of speech and movements/gestures from a **single high-density electrocorticography (ECoG) implant** in three participants with paralysis, including concurrent task attempts and a virtual-avatar interface. Training with both isolated and simultaneous task data helped across contexts. This is about **distributed concurrent neural representations and decoding**, not physical neural branch formation.

**Potential connection** [INTERPRETATION]: a physically adaptive network might face **multiple changing tasks at once**, so an appropriately costed growth rule must evaluate interference, task diversity, shared pathways, and redundancy rather than only one source–sink pair. A careful future benchmark could use two independent output objectives and compare a shared network vs separate pathways. **Not evidence** that a physically growing network performs better or optimizes empowerment.

### 36.5 Relevance ranking for near-term research

- **Core morphogenesis mechanics:** Le Verge-Serandour & Alim review + its primary experiments, Jones/Tero models, the *PhyChip* chapter, and DLA baseline.
- **Useful growth-task inventory:** Duggan's lecture deck (secondary, follow the primary references).
- **Important controls for 'intelligence encoded in structure':** distributed-EEG decoding and multimodal BCI papers; they highlight **readout learning and distributed activity** as alternatives to claiming a memory is structurally encoded.
- **Optional future visualization:** IFTF ONVR, kept **separate** from growth simulation and physical law.

## 37. Updated Adaptive Morphogenesis Lab experiment: measure the value of structural adaptation rather than assuming it

A fair comparison needs **three nested representations**:

```text
TOPOLOGY:       fixed nodes and graph edges, branch parent links, contact state
GEOMETRY:       centerlines + radii r(s,t), optional disks / union of occupied area / continuous fan φ(x,t)
DYNAMICS:       selectable edge model (graph flow, resistance, spring–damper,
                calibrated HardDisks transfer) + separate mobile tracer/nuclei model if used
```

**Required ablation:** compare (A) graph only with dynamic edge radii, (B) same graph plus disk-chain geometric collision/contact representation, and (C) an occupancy/level-set fan model for sheet-to-vein growth. Match available **physical** material area/volume, external driving and task; independently vary numerical disk count and grid resolution to ensure measured cost is not a discretization artifact. If simple graph+radii performs equally well on a given task, retain disks only for contact geometry/visualization; more simulation complexity has to earn its place.

**A clean scientific question, now sharpened by the sources:** Is task performance information encoded in (1) transient internal state (pressure/flow/tracers), (2) persistent structural radii and connectivity, (3) learned readout/controller parameters, or (4) environmental traces? **Causal ablations**: hold environment/task fixed while (i) reset only radius pattern, (ii) reset only readout, (iii) shuffle prior task history, (iv) preserve topology but randomize widths, and (v) change the environment while matching total materials. A morphology-reading decoder should be trained and tested out of sample with leakage controls.

**Benchmark environment:** two growing soma-like seeds in 2D for verifying exploratory contacts, then three or more resources/destinations under a changing demand Markov process. Begin with DLA random growth and a simple local flow-reinforcement policy before adding predictive IB or empowerment. Initially report **dimensionless physical-resource and operation-cost proxies**; only claim joules when a calibrated physical mechanism is in place.

**Experimental cost equation (explicitly *not* a universal biophysical law):**

\[
C_{\mathrm{total}}=C_{\mathrm{grow}}+C_{\mathrm{maintain}}+C_{\mathrm{flow/transport}}+C_{\mathrm{sense}}+C_{\mathrm{controller}}+C_{\mathrm{actuate}}+C_{\mathrm{prune}}.
\]

A positive decision to add branch `e` should reflect the *incremental expected task value over a stated horizon*, conditional on the network's information state, versus incremental total cost. Include **do not grow** as a legitimate action and use a resource-constrained static optimum and flow-only adaptation as baselines. For physical spring/EDMD branches, impose the correct drive/work budget; stored spring energy cannot be counted twice as transferred or created work.

## 38. Handoff update: new research questions for Fable/Astra/Opus 5.5

Review the **full v1.5 notebook** and the short `Fable_Research_Review_Handoff_2026-10-12.md` as a pair; do not limit critique to the screenshots or the aesthetics.

1. Verify **Fig. 2C** of Le Verge-Serandour & Alim: identify the fluorescent objects, the difference between advancing sheet/fan, vein/tube, internal transport, and nuclei trajectories. Would a disk-chain model be scientifically wrong as a *literal Physarum cell description*, while still valuable as a *controlled geometric abstraction*?
2. Can an effective graph-plus-radius model be **identifiable** from measured movies / flow? Which parameters can we infer reliably, and which require perturbations or distinct observation modalities?
3. Recommend the simplest **mass-conserving** discrete growth law and radius dynamics that can simultaneously produce exploration, tentative contacts, loops, thickening, retraction, and *no-growth* outcomes. Compare particle chains and a field/continuum model.
4. Specify which information is stored in physical **morphology** versus transient **activity** versus a trained controller/readout; use causal interventions to attribute improvements.
5. Which of the two uploaded BCI papers, if either, genuinely informs the physical-morphogenesis question versus only providing examples of distributed encoding/decoding? Avoid unsupported mechanistic analogies.
6. For the original five-paper program, audit one-sided versus two-sided Román equations, finite-frequency versus equilibrium stiffness, the driven-piston work ledger, acoustic impedance assumptions, and the physical/energetic limits of active controlled networks.
7. Identify the first experimental result that would justify upgrading the model from **graph+radii** to **disk chains** to **sheet-level Physarum** or **actual hard-disk edge dynamics**. State measurable go/no-go gates.

## 39. Provenance, bibliography supplement and change log

**New user-provided artifacts:** `Screenshot 2026-10-12 at 22.26.23.png` (Physarum fan fluorescence panel); `Physarum+Computing.pdf`; `Immersive2_Human_Networks_Report(1).pdf`; `Transfer learning via distributed brain recordings enables reliable speech decoding(1).pdf`; `simulataneous speech and gesture decoding for multimodal communication in paralysis(1).pdf`. The neuroplasticity and growth-cone clips discussed in §§17–34 remain visual/procedural inspiration with their own evidence limits.

**New external verification references:**

- Le Verge-Serandour & Alim (2024), *Physarum polycephalum: Smart Network Adaptation*, https://doi.org/10.1146/annurev-conmatphys-040821-115312; accessible manuscript https://arxiv.org/abs/2306.09063; Figure 2C and caption are critical to interpreting the screenshot.
- *Spatial transcriptomic and single-nucleus analysis reveals heterogeneity in a gigantic single-celled syncytium*, eLife (2022), https://elifesciences.org/articles/69745/figures; video supporting transport of nuclei through the fan.
- Duggan (2018), *Physarum Computing*, user-uploaded 16-slide class deck; **secondary** background.
- IFTF (2018), *Immersive Human Networks*, user-uploaded report; prototype of network visual analytics, **not** biology.
- Singh et al. (2025), *Transfer learning via distributed brain recordings enables reliable speech decoding*, https://doi.org/10.1038/s41467-025-63825-0 .
- Brosler et al. (2026), *Simultaneous speech and gesture decoding for multimodal communication in paralysis*, https://doi.org/10.1038/s41593-026-02446-2 .

**2026-10-12, v1.4:** Added §§28–34 (Fable handoff and mathematically explicit five-paper plan, graphs versus growing disks, material/energy ledger, concrete tests, visual source and critical reviewer request). **This v1.4 material was prepared earlier but had not been consolidated into the main notebook at the time the user asked whether saving worked.**

**2026-10-12, v1.5:** Added §§35–39, the verified fluorescent-fan interpretation, four uploaded PDFs, their distinct roles and evidence limitations, sharper experimental ablations, expanded Fable reviewer questions and bibliography. Prior original notebook sections 0–27 remain unchanged; the addendum §§28–34 is included verbatim. No new empirical HardDisks simulation results are claimed here.


## 40. [2026-10-12 / v1.6] Deep prior-art audit: **predictive information + empowerment + energetically costly structural adaptation**

**Status:** Literature-based *provisional* originality assessment, not an assertion of undiscovered science. **Search date:** 2026-10-12. Grounding: the four user-uploaded Alim-group papers (*Decision-Making in Light-Trapped Slime Molds*; *Efficient Transport Networks through Edge Pruning by Relative Resistance*; *Flow-Induced Vascular Remodeling on-Chip*; *Self-Organized Homogenization of Flow Networks*), prior uploaded biological/physical papers, and independently verified publisher/preprint/PMC results listed below. This is a structured search of accessible publications, **not a systematic review with guaranteed coverage of all work, unpublished manuscripts, dissertations, or conferences**. Never claim "no one has tried this" on this evidence.

### 40.1 What precisely is the proposed intersection?

The user asks a question narrower than "can physical networks learn?":

> **When does it pay, in genuine construction, maintenance and operating resources, for a locally adapting physical network to grow, thicken, preserve an apparently unused branch, or prune it; and can a policy grounded in Still's predictive information/thermodynamics of prediction and information-theoretic empowerment outperform simpler local transport heuristics under stochastic, changing future demands?**

The three proposed lenses need careful definitions:

1. **Prediction / predictive Information Bottleneck:** a state variable, such as physical morphology `M_t`, retains task-relevant information about future environmental demand `Y_{t+T}` without retaining all historical fluctuations. A prospective test is `I(M_t;Y_{t+T}|O_t)` or a representation-based predictive IB objective; do **not** conflate maximizing raw mutual information with IB compression. Verify the conditioning variable `O_t` and causal/data assumptions before interpretation.
2. **Empowerment / future controllability:** an *action-to-future-observation channel capacity* at a chosen time horizon with an explicit action set, sensor/readout definition, and feasible energetic control budget. It is **not** generic robustness, number of possible graph paths, extent of random branching, or uncertainty itself.
3. **Thermodynamics and resource costs:** distinguish dissipated work of an actual driven physical system (joules), maintenance/actuation expenses, building material or area (e.g. structural length/volume), computational power, and user-defined dimensionless cost proxies. Still's results do **not** imply `energy_to_grow_a_branch = k_B T ln(2) × its bits of structural memory`.

**Operational substrate:** `G_t=(V_t,E_t)`, edge lengths, effective radii `r_e(t)`, continuous geometric growth representation (optional disk chains or sheet/level-set), network physical state (flows/oscillations), controller actions and budget. HardDisks EDMD *thermal particles* remain distinct from structural disks. Its physically characterized edge response is a later replacement for initial ideal-resistor / spring–damper edge dynamics. No unmeasured hard-disk physical learning or energy savings are claimed.

### 40.2 Prior work directly threatens broad novelty — **HIGH confidence that individual pairs and much of the triple exist already**

**A. Susanne Still's thermodynamics of predictive memory (fundamental foundation; not a structural-network growth model).**

- Still, Sivak, Bell & Crooks (2012), **"Thermodynamics of Prediction"**, *Physical Review Letters* **109**, 120604. DOI: https://doi.org/10.1103/PhysRevLett.109.120604. Nonpredictive information maintained by a stochastically driven system is intimately linked to thermodynamic inefficiency/dissipation under the paper's protocol and assumptions. Do not apply the relation as a universal equality to every growth simulator.
- Still (2014), **"Information Bottleneck Approach to Predictive Inference"**, *Entropy* **16**, 968–989. DOI: https://doi.org/10.3390/e16020968. Predictive IB, efficiency and learning dynamics.
- Still (2020), **"Thermodynamic Cost and Benefit of Memory"**, *PRL* **124**, 050601. DOI: https://doi.org/10.1103/PhysRevLett.124.050601. A more general limit for partially observable information engines, including the costs of irrelevant retained information. This is a central theoretical foundation and more relevant to the proposed predictive controller than invoking only 2012.

**B. Empowerment and predictive/control information (not new, even in dynamical physical models).**

- Klyubin, Polani & Nehaniv (2008), **"Keep Your Options Open: An Information-Based Driving Principle for Sensorimotor Systems"**, *PLOS ONE* **3**, e4018. DOI: https://doi.org/10.1371/journal.pone.0004018. Formalizes empowerment as future sensorimotor influence. The paper also uses a bottleneck-like penalty in a context-controller demonstration; therefore, a generic claim that "empowerment has never been combined with an information bottleneck" is **false**.
- Tiomkin, Nemenman, Polani & Tishby (2024), **"Intrinsic Motivation in Dynamical Control Systems"**, *PRX Life* **2**, 033009. DOI: https://doi.org/10.1103/PRXLife.2.033009. Dynamical systems and empowerment-based intrinsic control, including tractable computation; especially pertinent to a mechanical simulator.
- Abel et al. (NeurIPS 2025), **"Plasticity as the Mirror of Empowerment"**: https://proceedings.nips.cc/paper_files/paper/2025/hash/f04957cc30544d62386f402e1da0b001-Abstract-Conference.html. Investigates dual influence directions (being affected by past versus affecting futures). **Crucial terminological guardrail:** their formal *plasticity* is an information-theoretic agent property, **not automatically growth-cone plasticity or edge addition**.

**C. The most important counterexample to a claim of novelty of the 'three elements' together.**

- **Chen & Prokopenko (2025), "Why Collective Behaviours Self-Organize to Criticality: A Primer on Information-Theoretic and Thermodynamic Utility Measures"**, *Royal Society Open Science* **12**, 241655. DOI: https://doi.org/10.1098/rsos.241655; full paper https://pmc.ncbi.nlm.nih.gov/articles/PMC12187418/ . **Read before claiming novelty.** It directly compares **predictive information, empowerment, active inference and thermodynamic efficiency** within an Ising perception–action model. Their objectives favor different parameter regimes; the work proposes a principle of super-efficiency. This is **not** Still's specific thermodynamics-of-prediction relation, nor a geometrically growing, finite-material transport network, but it absolutely invalidates "nobody combines prediction, empowerment and energy in collective behavior."
- Ueltzhöffer (2020), **"On the Thermodynamics of Prediction under Dissipative Adaptation"**, **preprint** https://arxiv.org/abs/2009.04006 . Tries to bridge Still-style prediction-efficiency ideas with self-organizing dissipative structures. Treat as *theoretical discussion / non-peer-reviewed status not verified*, rather than experimental proof of morphology-specific learning.

**D. Physical networks already display structural memory and construction-versus-transport trade-offs.**

- Kramar & Alim (2021), **"Encoding Memory in Tube Diameter Hierarchy of Living Flow Network"**, *PNAS* **118**, e2007815118. DOI: https://doi.org/10.1073/pnas.2007815118 . **Tube diameters are explicitly a persistent memory variable; tubes grow at the expense of others shrinking.** The proposal `r_e(t)` is therefore not original *as a concept*, but is an excellent experimentally grounded choice of state variable.
- **Schick, Kramar & Alim (2024), "Dynamic Cost Allocation Allows Network-Forming Forager to Switch Between Search Strategies"**, *PRX Life* **2**, 033005. DOI: https://doi.org/10.1103/PRXLife.2.033005 . **Extremely close precedent.** Measures/estimates construction versus transport energetic costs across morphologies and foraging/search strategies. Claims that the trade-off between structure-building and transport energy has never been studied in Physarum are incorrect. Important boundary: does not, by its abstract, establish a comparative test of Still-derived predictive bottleneck plus constrained empowerment against flow-only policies with *grow/no-grow* actions.
- Schick et al. (2026), **"Decision-Making in Light-Trapped Slime Molds Involves Active Mechanical Processes"**, *PRX Life* **4**, 023026. DOI: https://doi.org/10.1103/rv7g-d9kx. Observes exploration and changing contraction modes that ultimately reorganize toward an effective escape route. **A mechanical account of a particular behavior**, not proof that Physarum maximizes empowerment or predictive information.
- Le Verge-Serandour, Vasilev, Szeimies & Alim (2026), **"Efficient Transport Networks through Edge Pruning by Relative Resistance"**, **arXiv preprint v1, posted 2026-10-09** https://arxiv.org/abs/2610.13582 . Compare `ρ_e=R_e/R_{alternative}` and their pruning criterion `ρ_e>1`, under their stated resistor/flow assumptions. Local readout through test flows can encode globally determined path redundancy. This is a particularly strong **non-informational baseline**, not a universal energy-optimal rule for nonstationary future conditions, and not yet to be described as peer-reviewed on search date.
- Bouvard et al. (2025), **"Self-Organized Homogenization of Flow Networks"**, *Physical Review X* **15**, 041038. DOI: https://doi.org/10.1103/j5ch-4vkh . Erosion-based changing widths homogenize artificial flow networks: **demonstrated physical hardware response**, useful for radius-physics model design. Performance is defined in their specific flow/perfusion context.
- Mirzapour-Shafiyi et al. (2025), **"Flow-Induced Vascular Remodeling on-Chip: Implications for Anti-VEGF Therapy"**, *Advanced Functional Materials* **35**, 2501416. DOI: https://doi.org/10.1002/adfm.202501416 . Biological flow–chemical interaction changes vascular remodeling; *not* a simple universal reinforcement rule.
- Tero et al. (2010), **"Rules for Biologically Inspired Adaptive Network Design"**, *Science* **327**, 439–442. DOI: https://doi.org/10.1126/science.1177894 . Transport-informed reinforcement/routing baseline; Tokyo-area experiment, not proof of global optimality under all constraints.

**E. Information-guided network pruning and resource economy are already mature topics.**

- **Scholl, Rule & Hennig (2021), "The Information Theory of Developmental Pruning: Optimizing Global Network Architectures Using Local Synaptic Rules"**, *PLOS Computational Biology* **17**, e1009458. DOI: https://doi.org/10.1371/journal.pcbi.1009458 . Locally accessible Fisher-information-based synaptic pruning yields compact, functionally useful connectivity. **Very direct challenge** to the novelty of information-driven structural pruning, albeit a model neural system rather than fluid/EDMD transmission substrate with Still's dissipative cost.
- Bullmore & Sporns (2012), **"The Economy of Brain Network Organization"**, *Nature Reviews Neuroscience* **13**, 336–349. DOI: https://doi.org/10.1038/nrn3214 . Material/metabolic wiring costs versus valuable topological properties. A novel paper cannot simply rediscover network economic trade-offs.
- Du et al. (2019), **"Efficient Network Construction through Structural Plasticity"**, preprint https://arxiv.org/abs/1905.11530 . Neural network growth/pruning and hardware energy/latency benefit; abstract mathematical architecture, not growing physical soft matter.
- Li & Mao (2024), **"Training All-Mechanical Neural Networks for Task Learning through In Situ Backpropagation"**, *Nature Communications* **15**, 10528, DOI: https://doi.org/10.1038/s41467-024-54849-z (uploaded previously): physical parameters/learning in spring-based networks. This shows one can train material couplings, not that hard-disk gas edges are automatically superior.

**Evidence conclusion:** the **broad 'intersection of prediction, empowerment, network morphology and cost' already has direct and very close precedents.** The literature found so far does **not** establish that someone has already performed the *exactly specified, fully budgeted, physical edge-growth/pruning and causal morphology-memory benchmark proposed below*, but this is **an absence of an identified paper, not proof of novelty**. Search coverage and formulation details still need checking with Susanne, Fable, and active domain researchers; pay special attention to dissertations / preprints from these authors.

### 40.3 What might remain publishable? Narrow, falsifiable candidate contribution

**Candidate question N1 (stronger than a broad principle):**

> In a stochastically switching demand environment, under identical finite material, energy, sensing and control budgets, does a *locally implementable* policy that uses compressed predictive representations of demand and/or action-conditioned future controllability produce **measurably superior task performance per unit total cost** to flow-based reinforcement, relative-resistance pruning, Fisher-information importance pruning, greedy cost-aware policies and static optimal networks?

**Candidate question N2 (physical morphology as representation):** Does persistent geometry `M_t=(G_t,r_e(t))` carry **incremental out-of-sample predictive information** about future task requirements, *conditional on current environmental state or sensory readings*, and does selectively removing/shuffling that structural memory causally impair performance? If controller weights alone carry the memory, the morphology-specific claim fails.

**Candidate question N3 (transmission physics and scale transfer):** Does using response functions measured from hard-disk modules instead of fitted spring–damper, ideal resistor or telegraph-edge surrogates change the network's optimal topology or physically achievable information–cost Pareto frontier under the **same drives and work budget**? If simpler surrogates match, publish the reduction/negative result and do not claim an artificial mechanical brain.

**Possibly interesting but independent:** A physical instantiation of the above with active soft-matter/hard-disk/gas chambers. That is **engineering novelty only if demonstrated, energy-measured, and contrasted with a simpler actuator/spring or microfluidic architecture**. A pretty circle-chain animation is not proof of substrate-level learning.

### 40.4 Variables and careful mathematics (proposal; not a derivation of a universal bound)

Let `H_t` be recorded environmental history, `O_t` the present observation, `Y_{t+T}` a specified future task-relevant variable, `M_t` structural morphology (graph + radii), and `A_{t:t+T-1}` feasible controller/build/actuation actions. Network operating state (flow/pressure) and trained controller memory should be **recorded separately** from `M_t`.

A **conditional morphology-predictiveness** test is

\[
P_M(T) := I(M_t;Y_{t+T}\mid O_t).
\]

It is not automatically a *causal* contribution: compare paired interventions that reset, shuffle or transplant `M_t` while controlling environment, current flow and controller state. If `Y` and `M` share environment/task history, mutual information alone may reflect correlation.

An **IB representation policy** can be proposed as

\[
\min_{p(z\mid h)}\ I(Z_t;H_t)-\beta I(Z_t;Y_{t+T}),
\]

with observation/memory information constraints specified carefully, and with actions/feedback included in the predictive target if appropriate. For **physical morphology** the encoder is realized by allowable growth/radius/retraction dynamics, *not an unconstrained digital encoder*. The gap between an optimal abstract IB representation and an attainable physical mechanism is itself measurable.

A **budget-conditioned empowerment** is

\[
\mathcal E_T(m)=\max_{p(a_{t:t+T-1}\mid m):\ \mathbb E[E_{act}]\leq B_{act}}
I(A_{t:t+T-1};O_{t+T}\mid M_t=m).
\]

Empowerment depends strongly on what counts as an action, whether growth is an action, the sensory output, time horizon, constraints, and feasible distributions. Report those or the number has no interpretation. Merely storing unused redundant loops is **not** empowerment unless they enlarge the feasible distinguishable action–outcome channel.

For an actual tube of length `L_e` and radius `r_e`, one *illustrative* material proxy is `V_e≈π r_e^2 L_e`; use actual occupied union area in 2D structural-disk representations to avoid double-counting and avoid counting numerical disks as energy. `R_e∝L_e/r_e^4` assumes stationary laminar Newtonian Poiseuille flow in a rigid circular tube. It is **not** a general resistance law for axons, peristaltic Physarum, or hard-disk transmission modules.

An explicit **decision ledger** over horizon `T` is

\[
C_T=C_{build}+C_{maint}+C_{transport}+C_{sense}+C_{infer}
+C_{actuate}+C_{reconfigure}+C_{dispose},
\]

with separate material quantity `B_M` (mass/area/volume, recycled fraction specified) and system boundary. For genuine hardware, separate injected work, heat, useful output work, and stored/recovered energy, applying the first law and validated EDMD work ledger. In an abstract simulator, `C` is a **chosen cost functional, not measured dissipated joules**.

Rather than adding arbitrarily weighted bits, joules and dimensionless task reward, start with a constrained frontier:

\[
\min_\pi\ \mathbb E[C_T]
\quad\mathrm{s.t.}\quad
\mathbb E[L_{task}]\leq \epsilon,\qquad \mathbb E[M]\leq B_M,
\]

and compare the Pareto set on **task loss, lifetime cost, predictive information and attainable empowerment**. Treat information measures as hypotheses/diagnostics or proposed auxiliary terms, and determine their marginal predictive power beyond trivial network size, flow, current observations and baseline performance. A *scarcity/value* statement about growing edge `e` is conditional: `E[Δ U_task | information state] > Δ C_e`, with compatible units and specific time horizon. No tautological payoff coefficient or free-energy identity is implied.

**Still-specific guardrail.** It is legitimate to ask whether retaining nonpredictive morphological history increases dissipation under a *specified model of driven state transitions* analogous to Still et al. (2012/2020). It is **not legitimate** to assert that those papers already prove a lower bound on biochemical ATP for creating a tube, or that IB-optimal morphology minimizes total organismal energetic cost in general. Distinguish information-engine dissipative-memory overhead from synthetic material purchase and active adaptation.

### 40.5 First convincing benchmark, designed to fail fairly

**Environment:** 2D initially sparse network with 3–8 terminals; stochastic source/sink demand switches with controlled predictability. Two matched environments: (i) Markov transition structure learnable from limited history and (ii) unpredictably resampled demand with the same marginal frequencies. Add rare route failure/damage and a held-out distribution shift. No policy gets privileged access to future demand.

**Network mechanics:** first ideal electrical/hydraulic resistor graph with dynamically adapting radii and conserved build material, optional geometry of branching disk chains at graph contact level; later replace edges with *independently fitted and validated* spring–damper and EDMD transfer kernels. Contain active drive energy in the ledger. Simulate real `no growth` and `prune`, allow loops and material recovery under measured/declared efficiency `ζ<1`.

**Policies (comparative controls, not one combined 'AI' score):** random growth/pruning; greedy current-flow Tero-inspired reinforcement; **relative-resistance pruning** (2026 preprint); locally available Fisher-information importance (Scholl-inspired); fixed robust topology; receding-horizon cost-aware controller with observed present demand; predictive IB representation policy; energy-constrained empowerment policy; factorial combination predictive IB + empowerment; oracle with privileged future demand as **upper bound only** (never as a fair competitor).

**Measurements:** held-out served demand and task error; transported quantity/useful work and latency; actual material area/volume, branch lifetime, memory reset/relearning time, construction and disposal losses, controller work/compute, transmission work; cost per successfully served demand; physical `P_M` and conditional predictiveness with CIs; empowerment estimation bias; topology (loops, degree, resistance, path redundancy). Use matched random seeds, preregistered parameter sweeps, finite-sample bias correction and ablations; vary horizon, environmental predictability, material price, repair price, and damage rate.

**Necessary falsification tests:**

1. **Structure really encodes history:** train/predict future demand from `M_t`; shuffle radii across edges while preserving total material and degree distribution; causal paired morphology transplant/reset; evaluate held-out tasks.
2. **Predictivity really matters:** predictive policy should help in learnable transition regimes, not automatically win in i.i.d. demand. Test leakage (no access to evaluation futures) and compare against a simple adequate Markov predictor.
3. **Empowerment really matters:** quantify distinguishable future outcomes under same action-energy budget; compare damage/novel demands at multiple horizons, not just extra backup edges.
4. **Energy really matters:** sweep `C_build/C_operate` and maintenance/resource prices; exhibit crossover boundaries and a strategy that chooses **not** to grow.
5. **Physically richer edge really matters:** compare identical network tasks with calibrated low-order spring–damper, resistor, and hard-disk module kernels; do not charge one model real joules and another only CPU proxy.
6. **Geometry resolution not a confound:** compare graph + edge radii with disk-chain contact geometry and optional continuum Physarum fan, using the same occupied area; change disk granularity while holding physical morphology constant.

**Potential publication unit:** a reproducible benchmark or short theory/simulation paper with precise, *negative-capable* finding: a finite material/energy versus future-demand-uncertainty diagram of which local/prioritized growth rules are Pareto-competitive. The five-paper dissertation remains the separate main staircase; this is a **candidate side branch** unless approved with Susanne and shown to yield novel results.

### 40.6 Prior-art confidence and action items

**High confidence — not novel in isolation:** Still predictive-memory/dissipation connection; empowerment as future influence; comparing predictive information, empowerment and thermodynamic efficiency; flow-adaptive tube thickness; physical morphology recording nutrient history; construction-versus-transport cost trade-off; informative synaptic pruning; local structural pruning; trainable spring networks. **The user has correctly found an active field, not invented it from scratch.**

**Moderate confidence — promising unverified intersection:** a controlled, explicitly material/energy-constrained, **physically instantiated and morphologically growing** network that measures and compares *Still-derived nonpredictive memory/dissipation, predictive IB, and budgeted empowerment as competing or combined local growth/pruning policies*, with rigorous baselines and changing future tasks. Search has **not identified an exact match**, but cannot exclude one; in particular review physical reservoir/embodied information theory, adaptive network growth, neural efficient wiring, active inference, and unpublished doctoral theses.

**Low confidence / unacceptable claim now:** that Nature uses IB or empowerment universally; that an optimized flow network has learned a causal model; that spring/gas-box or disc-chain growth is analogous to synaptic biochemistry; that complex EDMD edges must compute better than low-order models; that a single formula generalizes all networks; that the above research already supplies a new theorem.

**Researcher handoff (Fable, Astra, Opus, Susanne):** Start with **Chen & Prokopenko 2025, Schick/Kramar/Alim 2024, Scholl/Rule/Hennig 2021, Still 2012+2020, Kramar/Alim 2021, and the October 2026 relative-resistance preprint**. Search these papers' citing/cited-by networks. Decide whether the most defensible novelty is (a) a new bound in a specified physical model, (b) evidence of morphological predictive memory under budgeted growth, (c) a comparative benchmark with honest energy accounting, (d) a minimal physical implementation, or (e) **no novel result at present**. Asking which experiment disproves the hypothesis is mandatory.

### 40.7 Directly relevant source list (verified publication venues/DOIs, as of 2026-10-12)

1. Still et al. 2012 PRL, https://doi.org/10.1103/PhysRevLett.109.120604 .
2. Still 2014 Entropy, https://doi.org/10.3390/e16020968 .
3. Still 2020 PRL, https://doi.org/10.1103/PhysRevLett.124.050601 .
4. Klyubin et al. 2008 PLOS ONE, https://doi.org/10.1371/journal.pone.0004018 .
5. Tiomkin et al. 2024 PRX Life, https://doi.org/10.1103/PRXLife.2.033009 .
6. Chen & Prokopenko 2025 Royal Society Open Science, https://doi.org/10.1098/rsos.241655 .
7. Abel et al. 2025 NeurIPS, https://proceedings.nips.cc/paper_files/paper/2025/hash/f04957cc30544d62386f402e1da0b001-Abstract-Conference.html .
8. Ueltzhöffer 2020 preprint, https://arxiv.org/abs/2009.04006 .
9. Kramar & Alim 2021 PNAS, https://doi.org/10.1073/pnas.2007815118 .
10. Schick, Kramar & Alim 2024 PRX Life, https://doi.org/10.1103/PRXLife.2.033005 .
11. Schick et al. 2026 PRX Life, https://doi.org/10.1103/rv7g-d9kx .
12. Le Verge-Serandour et al. 2026 **preprint**, https://arxiv.org/abs/2610.13582 .
13. Bouvard et al. 2025 PRX, https://doi.org/10.1103/j5ch-4vkh .
14. Mirzapour-Shafiyi et al. 2025 Adv. Funct. Mater., https://doi.org/10.1002/adfm.202501416 .
15. Tero et al. 2010 Science, https://doi.org/10.1126/science.1177894 .
16. Scholl et al. 2021 PLOS Computational Biology, https://doi.org/10.1371/journal.pcbi.1009458 .
17. Bullmore & Sporns 2012 Nature Reviews Neuroscience, https://doi.org/10.1038/nrn3214 .
18. Du et al. 2019 **preprint**, https://arxiv.org/abs/1905.11530 .
19. Li & Mao 2024, *Nature Communications* **15**, 10528, DOI: https://doi.org/10.1038/s41467-024-54849-z (verified publisher record, original uploaded PDF).

**Search terms tried:** `Still thermodynamics of prediction morphology network`, `information bottleneck empowerment thermodynamics growth`, `empowerment morphological development cost network`, `predictive information energy neural pruning`, `Physarum construction transport cost`, `structural plasticity Fisher information`, `thermodynamic prediction dissipative adaptation`, and close variants; scholar/DOI/APS/Royal Society/PLOS/PNAS publisher search. **Limitation:** search is broad and targeted but not comprehensive; all claims of originality remain provisional. Exclude papers with **publication dates after 2026-10-12** from the as-of-date assessment.

**Revision log — v1.6 (2026-10-12):** appended this audited prior-art section to existing §§0–39 *without deleting earlier content*. Theoretical and experimental claims are labeled; major novelty threats highlighted; actionable baseline comparisons specified. Earlier five-paper equations and EDMD validation cautions remain unchanged.


---

## 41. [2026-10-13 / v1.7] Grasso–Bongard / Levin research line: **empowerment, living tissue, and homeostatic growth**

**Provenance:** User provided `https://github.com/caitlingrasso/bio-connectivity`, `https://github.com/caitlingrasso/empowered-nca-II`, and the full article **Grasso & Bongard (2023)**, *Selection for short-term empowerment accelerates the evolution of homeostatic neural cellular automata* (GECCO, DOI [10.1145/3583131.3590469](https://doi.org/10.1145/3583131.3590469); version of record and [arXiv:2305.15220](https://arxiv.org/abs/2305.15220)). These are **strong prior-art constraints** on claims that empowerment-driven morphogenesis or homeostasis is new. Grasso & Bongard's previous *Empowered Neural Cellular Automata* (2022) is [arXiv:2205.06771](https://arxiv.org/abs/2205.06771), code [`empowered-nca`](https://github.com/caitlingrasso/empowered-nca). The 2023 repository's README explicitly notes *the arXiv version contains method corrections and updated data*; resolve reported population sizes/figures against the exact manuscript/repo commit when reproducing. No reproduction of their results has yet been completed by us.

### 41.1 `empowered-nca-II` — what it does and what it does *not* demonstrate

Repository: <https://github.com/caitlingrasso/empowered-nca-II> (directly inspected 2026-10-13). Components include `ca_model.py` (cell update model), `main.py` (evolution launcher), `additional_objectives.py` (ablation objectives), `finetune.py` (transfer to other targets), `visualize_sim.ipynb`, `analysis/`, and configurations/target shapes. It is a source of **comparative baselines and methodological controls**, not a direct physical growth-energy simulator.

- One seed cell on a **discrete 2D grid** grows toward a predefined target shape. Cells have binary live/dead status and a separate local signaling channel. Their neural rules are **evolved** using multiobjective Age-Fitness Pareto Optimization, rather than learned online from changing service demands.
- The investigated auxiliary information objective is a **time-lagged mutual information of emitted cell signals (actions) with later observed neighboring signals (sensors)**, referred to as empowerment. For exact conceptual comparisons, distinguish this *empirical MI under the observed action distribution* from classical state-conditional **channel-capacity empowerment** that maximizes MI over an action policy. They are not automatically numerically equivalent.
- Shorter lag `k`, particularly `k=1` in these tested setups, generally improved target-shape loss, stability, and adaptation to related shapes compared with longer lags and relevant loss-only controls. **It does not imply `k=1` is universally optimal for a physical network.** Transmission length, slow mechanical dynamics, sensing time and target statistics could create different horizon dependence.
- **Very important adverse control:** the paper reports *minimizing individual cells' action entropy* as an auxiliary objective performing **better than short-lag empowerment on the tested morphogenesis task**. Do not build an argument on empowerment necessarily outperforming simple local regularizers. Inspect their `additional_objectives.py` and use those comparisons.
- The model does not explicitly charge for physical tube volume, polymerization, nutrient transport, retraction work, metabolic maintenance, or thermodynamically measured actuator work. It optimizes *desired shape + information/other algorithmic objectives* in an abstract NCA.

### 41.2 `bio-connectivity` — functional network inferred from cellular signals

Repository: <https://github.com/caitlingrasso/bio-connectivity> (README inspected 2026-10-13). Associated work: Blackiston, Dromiack, Grasso et al., *Revealing non-trivial information structures in aneural biological tissues via functional connectivity* (PLOS Computational Biology, 2025; <https://doi.org/10.1371/journal.pcbi.1012149>). Its analysis infers **functional-connectivity networks** from pairwise mutual information of calcium-imaging traces in *Xenopus laevis* organoids before and after perturbation, using JIDT and Brain Connectivity Toolbox. **An MI edge is not proof of an anatomical contact, causal influence, or physical material connection**: common drive and indirect pathways may generate functional dependence. This motivates measuring **three different networks** in our future simulator: (a) actual material geometry/topology, (b) signal-activity-based functional connectivity, and (c) causal influence under intervention. Distinguish morphology encoded memory from learned readout/observer memory.

### 41.3 What genuinely differentiates our planned experiment?

The proposal is **not** a first connection between empowerment and morphogenesis. It asks *which environmental uncertainty, finite resource and physical transport regimes favor one adaptive architecture over another*, when both connectivity and tube/disc radii can change. Information may be retained in a **finite controller state** `Z_t`, in the **physical morphology** `M_t=(G_t,{r_e(t)})`, or jointly. Test the contributions separately.

- **Predictive memory:** Dynamical IB can learn compressed summaries of observed histories, but Still's cited **no-feedback** result cannot be treated as a theorem about a closed-loop, materially growing agent. First use the trained IB representation as a *predictor* with its own sensing, learning, compute and storage ledger.
- **Future control:** time- and energy-constrained empowerment depends on *available physical actions, future observed outputs, and an action policy*. It is not equivalent to maximizing the number of branches or sampled signal MI.
- **Physical cost:** energy/material budgets must be non-circular, explicitly stated, and **resolution invariant** (i.e. adding rendering discs should not artificially add cost). Material is not automatically energy, and user-defined abstract penalties are not calibrated joules.
- **Alternatives:** current-flow reinforcement, relative-resistance pruning (Le Verge-Serandour et al. 2026), random retention, fixed redundancy, simple probabilistic forecasting, Fisher/sensitivity pruning (Scholl et al. 2021), and action-entropy minimization must compete fairly.

**Review takeaway:** Grasso/Bongard, Levin/Blackiston and Alim are *relevant research communities and comparison points*, not proof of our specific proposal's novelty, and not evidence that a HardDisks module already computes or learns. User is interested in the intersection of thermodynamics, structure, and prediction, rather than a purely aesthetic/artificial-life demo.

## 42. [2026-10-13 / v1.7] **Smallest parallel adaptive-network experiment: build/retain/prune under predictive demand**

**Status: proposed feasibility experiment, not performed, not yet validated or approved as a new PhD paper.** Develop as a separate `Adaptive Morphogenesis Lab` proof-of-concept, organized in IdeaLabs. **Do not modify the validated HardDisks EDMD scheduler or accidentally redefine its physical energy ledger.** The simplest scientific model is a graph with physically meaningful edge radii. Add explicit attached structural discs for exploration/encounters only if they yield a testable distinction from graph-only growth.

### 42.1 The tiny two-destination example (derive a decision threshold)

Imagine **one source with existing pathway A** and an underused alternative branch **B**. A and B may be thought of as separate downstream destinations or as alternative response pathways. The environment switches its required destination. At a decision point, it can **retain B** for the next time step, or **prune B** and later rebuild if needed. The following *invented illustrative numbers* are in the **same abstract cost units per decision interval** (not measured joules):

- Maintenance of unused B for one interval: `m = 0.8` units.
- If B is pruned and a demand for B occurs, rebuilding costs `C_rebuild = 3` units and delayed/missed service has **energy-equivalent opportunity penalty** `D = 5` units. The latter is **not energy dissipated by the network**; it is a deliberately defined common utility/cost accounting unit.
- The probability of B being needed next interval, conditional on currently observable history, is `p`.

Assuming retaining B avoids *both* rebuild and delayed-service penalties if B is needed, while all other consequences cancel between the alternatives,

`cost(keep B) = m`, and `expected cost(prune B) = p (C_rebuild + D)`.

Therefore the ideal one-step decision rule under these **particular** assumptions is

\[
\boxed{\text{retain B iff }p(C_{\mathrm{rebuild}}+D)>m,\qquad
p^*=\frac{m}{C_{\mathrm{rebuild}}+D}=\frac{0.8}{8}=0.10.}
\]

At `p = 0.05`: expected prune cost `0.4`, smaller than keeping `0.8`; **prune**. At `p = 0.20`: expected prune cost `1.6`, larger than keeping `0.8`; **retain**. Equality is indifference under the toy assumptions.

**This inequality is elementary expected-value decision theory, not a new theorem about biology or a bound from Still.** Its value is as a **clean test target**: does a constrained adaptive controller approach the appropriate decision threshold, and can a *dynamical IB state* correctly distinguish low/high future probability from the observation history? If an oracle predictor always achieves it, that merely validates the implementation. Novelty requires meaningful comparisons under limited sensing, finite memory, actuation/repair costs, changing tasks and correctly chosen baselines. For multi-step horizons, delayed benefits, topology interaction, and radii changes, the one-step threshold changes; use dynamic programming/MDP or direct simulation rather than extrapolating the formula.

### 42.2 Phase A: abstract scientific graph (recommended first)

- **Topology and world:** one source, destination A and B, initial optional alternative branches; later 4- and 20-node graphs. Graph `G_t=(V,E_t)` with lengths `L_e`, radii `r_e(t)` and physically parameterized conductance `g_e(r_e,L_e)`. Specify whether the transported quantity is viscous fluid, a scalar signal, or a mechanical response; **do not interchange their conductance formulas**.
- **Demand process:** two-state Markov environment with either strong autocorrelation/learnable switching or IID demand, matched to have the **same stationary destination frequencies**; variants introduce shocks, damaged edges and changing transition probabilities. Avoid letting predictors see future states or latent transition parameters unavailable to baseline controllers.
- **Material and energy:** a defined cap `M_t <= M_max`; e.g. `M_t = sum_e π r_e(t)^2 L_e` for nonoverlapping ideal 3D cylindrical tubes, while 2D area/union is appropriate for planar discs. Differentiate material stock from supply, construction/maintenance operating power, energy recovered upon pruning, dissipated work, sensing, learning and controller cost. `E_lifetime = E_transmission+E_construction+E_maintenance+E_actuation+E_sensing+E_control` only if every term really uses compatible units and accounting boundaries. Keep missed-service penalty separate as `C_task` (or explicitly calibrate it to resource/energy equivalent).
- **Candidate policies:** (0) no-growth and fixed sparse/redundant networks; (1) random/provisional; (2) simple recent-flow update, (3) relative-resistance pruning or equivalent topology rule, (4) small Markov or Bayesian predictor, (5) verified Dynamical IB predictor + rational keep/prune decision, (6) action-/energy-constrained empowerment, (7) combined policies *only after* separate performance is known. Include simple local action-entropy regularizer as a strong adversarial baseline inspired by Grasso–Bongard.
- **Metrics:** service accuracy/regret and delays under shocks, material stock, cumulative material moved, total energy proxy and separately calibrated energy, number/rate of growth and pruning decisions, connected components/loop redundancy, structural memory `I(M_t;D_{t+τ}|X_t)` estimated on independent trajectories, prediction score, constrained empowerment with explicit horizon/actions, and repeated-seed confidence intervals.
- **Causal test:** freeze environmental and controller states then swap/erase/reset **just the morphology**; compare future task performance. Otherwise improved prediction may reside only in the controller, not the structure.
- **Scientific success:** a clear, statistically supported **parameter regime** where a policy beats simpler baselines on a *pre-registered cost–performance trade-off*, not just interesting morphology. **Valid negative finding:** no additional value for IB/empowerment after simple forecasters, relative-resistance rules and controller overhead are included.

### 42.3 Phase B: the user's preferred **disk-chain growth** (only if justified)

This remains a worthwhile **mechanistic geometric model**, even if a graph-only model suffices for the first information-theory tests:

1. Two soma-like seeds are large, *stationary structural discs*. Each produces **attached discs** along a persistent, partly stochastic growth direction. Parent `i` and child `j` start touching with `x_j=x_i+(r_i+r_j)d_hat`; the connected disc centers form a geometric skeleton; attached material is not a cloud of elastically colliding EDMD gas particles.
2. Tips extend or retract, branch, encounter other tips, form **provisional** contacts, then either stabilize/thicken or retract. Collisions/contact are geometric, not a proof of true synapse formation; neuron guidance is simplified, not mechanistically recreated. Model growing front/occupancy separately in Physarum mode, which is a continuous multinucleated sheet, not a pile of distinct cells/discs.
3. Evolve each disc radius `r_i(t)` or centerline radius `r(s,t)` as an *explicit state*. Draw smooth **channel boundaries** around the centerline by interpolating radius and normal directions; **touching discs themselves have zero-throughput contact area**, so the tube boundary is a defined geometric construct, not a hydraulic path through point contacts.
4. Charge costs for **union** occupied area or tube volume and changing material, not raw disc count; enforce consistent minimum feature size and validate convergence when doubling the spatial resolution. Model local mobility, bends and stiffness only in a new active/mechanical module with explicit forces and work terms.
5. Use graph-based flow models to reinforce, thicken or prune formed channels. For future physical transmission, fit a spring–damper/response surrogate and later an independently characterized EDMD box. **Many-body EDMD cannot be relabeled as growing biological tissue while retaining its original conserved Hamiltonian dynamics.**

**Decision gate:** demonstrate that disc geometry changes a *measured outcome* (e.g. encounter probability, collision feasibility, branching path length, local material competition, prediction/service) beyond a comparably parameterized graph/radius model before adding it to the scientific benchmark. A 3D renderer is a separate visualization project, not evidence of learning.

### 42.4 Lightweight experiment build plan — bounded, parallel, reproducible

This is a **proposal for the user's side project**, not a request to replace the five-paper PhD roadmap:

| Increment | Scientific deliverable | Explicit non-goal |
|---|---|---|
| 0. Fix the hypotheses | Define action/state/demand observables, common energy/material/task ledger, policies and blind held-out tests | No claim of novelty or physical joules yet |
| 1. Small graph simulator | Seeded A/B demand process, radius-adjustable two-path graph, deterministic tests for conservation/budget and the `p*=0.10` toy threshold | No Physarum/neuronal realism |
| 2. Local adaptation baselines | Random, keep-all, flow reinforcement, simple Markov forecast, relative-resistance and entropy controls, reproducible plots | No imported multiagent/IB complexity |
| 3. Information methods | Use *validated* DIB-derived predictive states; then add carefully specified constrained empowerment; ablations, MI estimator uncertainty, horizon sweep | No assertion Still's lower bound equals growth cost |
| 4. Mechanics/morphology | Optional radius chains, geometric growth cones, structural-memory-reset intervention; spring–damper surrogate | No direct modification to `hspist3` |
| 5. Escalation decision | Short report containing success/failure/negative findings, matched budget curves, confidence intervals, reproduction instructions, and an explicit go/no-go with Susanne/Fable | No full 3D artificial organism by default |

**Practical timebox:** user intends to build **in parallel**. Keep this a narrow, scheduled side stream (illustratively around 10% of research effort over a four-to-six-week feasibility period), with Paper 1 Román validation and Paper 2 energy ledgers the primary commitments, and the nearly completed Dynamical IB package validated on its own terms. Reassess priority based on results and Susanne's guidance. This allocation is an advisory safeguard, not a restriction on future research interest.

## 43. [2026-10-13 / v1.7] The decision we reached: **protect the papers; don't discard the morphogenesis question**

The user's recent decision, preserved close to their wording: **“ok i concentrate on my other paper first”**, then **“ill build this parallel.”** Interpret these together: main PhD energy goes into the established HardDisks papers, while a *small parallel exploratory adaptive-network prototype* can progress without taking over the thesis. The five-paper staircase in §29 remains the long-term plan, **not a claim that five papers have been completed or committed**.

The immediate research priorities, previously discussed, are: **(1)** finish Román static/dynamic stiffness, sound-speed and phase-related validation without overclaiming; **(2)** finish spring-backed piston energy transfer with a complete work ledger and controls; **(3)** complete independent Dynamical IB implementation tests; **(4)** only then decide whether adaptive structural learning earns expansion. Existing numerical/physical results stand or fall independently of the morphogenesis idea. The user is interested in an eventual *Physarum/neuron-like, adaptive, possibly hybrid bioelectrical/mechanical system* and a beautiful 3D simulation, but those are **later directions** contingent on strong simpler results.

**Keep the key principle prominently in handoffs:**

> *The research opportunity is not just to combine three established information-theoretic words. It is to derive or test a prediction that the separate principles or simpler local rules cannot already explain, such as a quantitatively defined threshold at which a constrained network preserves an apparently unneeded branch because it predicts future demand—and to establish the lifetime price of that structural memory.*

**Unproven novelty / falsification:** Even the more specific threshold experiment may be a reproduction of known sequential decision theory or stochastic network optimization. A scientific contribution may instead require a nontrivial physical mechanism, a novel derived bound under defensible assumptions, or experimentally verified predictions with out-of-sample comparisons. Fable should actively seek counterexamples in operations research, stochastic optimal networks, computational mechanics, physical reservoir computing, embodied developmental models, and sensorimotor information theory.

### 43.1 Short handoff to Fable / Astra / Opus 5.5

1. **Do not start a new giant repo.** First inspect prior notebook §§29–32 and 40–43, the prior-art companion, and the explicit `p*=0.10` toy. Retain separation of validated HardDisks physics vs speculative network growth.
2. **Audit six source papers:** Grasso & Bongard (2022, 2023); Salge–Glackin–Polani on classical empowerment; Scholl–Rule–Hennig on developmental pruning; Still's Dynamical IB; Chen–Prokopenko on objective conflict; Schick–Kramar–Alim on cost allocation. Add Kramar–Alim structural memory and latest Alim local pruning. Note the 2023 empowerment-versus-action-entropy reversal and the no-feedback condition in Still.
3. **Propose a one-page, executable minimal benchmark** with clear input process, timeline, two-path topology, actions, appropriate simple controls, matched lifetime accounting and a falsifiable quantitative outcome. Decide whether discs buy us anything for its first phase.
4. **Ask for the strongest contrary paper:** has the same threshold, empowerment objective, predictive memory, and physical topology/resource cost already been treated in stochastic network design? Search proactively for failed-novelty cases.
5. **Do not confuse a nice 3D morphology with proof of learning**, sampled action-sensor MI with classical channel capacity, structural calcium MI with physical contact, abstract penalties with physical joules, or closed-loop actuator costs with Still's open/no-feedback memory-machine bound.

**Sources added for this version:**

- Grasso & Bongard 2022, *Empowered Neural Cellular Automata*, <https://arxiv.org/abs/2205.06771> (code <https://github.com/caitlingrasso/empowered-nca>).
- Grasso & Bongard 2023, *Selection for short-term empowerment accelerates the evolution of homeostatic neural cellular automata*, <https://doi.org/10.1145/3583131.3590469> and <https://arxiv.org/abs/2305.15220>; code <https://github.com/caitlingrasso/empowered-nca-II>. Repo README flags corrected methods/figures in arXiv revision.
- Blackiston, Dromiack, Grasso et al. 2025, *Revealing non-trivial information structures in aneural biological tissues via functional connectivity*, <https://doi.org/10.1371/journal.pcbi.1012149>; code <https://github.com/caitlingrasso/bio-connectivity>.
- Related benchmark: Scholl, Rule & Hennig 2021, <https://doi.org/10.1371/journal.pcbi.1009458>; Schick, Kramar & Alim 2024, <https://doi.org/10.1103/PRXLife.2.033005>; Still 2012, <https://doi.org/10.1103/PhysRevLett.109.120604>.

**Revision log — v1.7 (2026-10-13):** Appended §§41–43 with the last two research discussions, concrete A/B cost-threshold example, disk-chain/graph implementation choice, bounded parallel-build plan, validated current priorities and Fable review instructions; all preexisting sections retained.


---

## 44. [2026-10-13 / v1.8] Corrected unified hypothesis — four **coupled roles**, not competing objectives

**Revision guard:** This addendum supersedes the *framing* of the earlier comparison of IB vs empowerment vs energy where it portrayed these as alternative overall growth policies. Earlier sections remain historically preserved as thinking/proposed baselines, not deleted. The new mechanism combines different causal roles; the previous competitive policies survive as ablation controls.

# Integrated adaptive morphogenesis — coupled-role prior-art audit and falsifiable research plan

**Research memory / handoff addendum | 2026-10-13 | v1.8 | PROPOSED, NOT YET DEMONSTRATED**

**Original clarification from the researcher (preserve intent):**
> I never wanted to say that I'm the first to test empowerment, IB or energy each alone on this problem, but that empowerment maybe is for expansion in slime mold, IB for strengthening, DIB for maintaining [knowledge] about the future, and energetic efficiency as the main build — all together combined responsible for the kind of network growth.

**The precise thesis:** Test whether *different information-theoretic functions* can govern **distinct physical decisions in one growing, learning, resource-limited network**, rather than merely compare empowerment, IB, and energy as interchangeable optimization scores.

- **Empowerment / future action capacity → exploratory candidate formation and future option preservation.** Hypothesis: *may* favor exploration when new morphology increases attainable-and-observable futures. It does **not** imply universal expansion: staying put or conserving an existing channel can maximize empowerment.
- **IB / task-relevant information selection → selective reinforcement or shrinkage.** IB acts on a well-defined *representation and relevance variable*, **not directly on branch radius**. A separately stated actuator rule translates representational relevance or marginal counterfactual value into `r_e(t)` / `r_i(t)` updates.
- **Dynamical IB / economical predictive state → temporally informed retention.** A predictive memory of repeated environmental transitions could prevent premature loss of dormant connections; do **not** double-count IB with DIB. Two variants: (a) separate controller memory; (b) network morphology as the only persistent memory.
- **Physics / thermodynamic and material ledger → hard budgets and lifetime selection.** Growth, maintenance, transmission, sensing, control, and removal have different measurable or assumed costs. Still's information–dissipation results **do not equate mutual information to literal construction energy**.

The mechanisms share an evolving graph, structural geometry, environmental data and finite resources. Their interaction—not any one isolated ingredient—is what we propose to test. The ultimate test is *a quantitative prediction that does not follow from individual ingredients or simpler local rules*.

## A. Evidence and novelty audit — what's been done, and what is NOT established

| Source and identifier | Actual supported contribution | Threat to originality / relevance | What our exact test adds (if it works) |
|---|---|---|---|
| Still, Sivak, Bell & Crooks (2012), **Thermodynamics of Prediction**; [PRL DOI](https://doi.org/10.1103/PhysRevLett.109.120604) | For specified stochastic-driven dynamics, nonpredictive memory relates to thermodynamic inefficiency/dissipation. | Very high for claiming the information–energy connection itself. | Construct an explicit morphology-mediated learner and **derive, not assume**, whether a thermodynamic statement applies to its feedback-driven growth. |
| Still, **Dynamical Information Bottleneck Method** (user's uploaded draft `DynamicIB_2022-2.pdf`, dated December 2021) | Learner update `p(s_{t+1}|s_t,x_t)` compresses while retaining future predictive information, motivated by a dissipation bound; **no feedback** in the analyzed subclass. | High for novelty of time-dependent predictive memory. | Predictive state governs costly *physical retention* decisions; no transfer of bound without new derivation. Verify final publication status separately; do **not** cite the uploaded draft's placeholder DOI as final. |
| Klyubin, Polani & Nehaniv (2008), **Keep Your Options Open**; [PLOS DOI](https://doi.org/10.1371/journal.pone.0004018) | Empowerment influences adaptive sensor/actuator morphology with constraints, formalizing perceived future control. | Very high: **empowerment-driven morphology** is not new. | Distinct mechanical construction, pruning and predictive retention within *transport* networks, with matched tasks and lifetime material cost. |
| Anthony, Polani & Nehaniv (ECAL 2009 proceedings), **Impoverished Empowerment**; [DOI](https://doi.org/10.1007/978-3-642-21314-4_37) | Uses Information Bottleneck-like action-sequence compression under limited sensorimotor bandwidth. | High for claims that empowerment and IB were never combined. | IB selects task-relevant **physical pathways**, distinct from action-sequence compression, within an explicit growth model. |
| Grasso & Bongard (2022), **Empowered Neural Cellular Automata**; [DOI](https://doi.org/10.1145/3520304.3529067), [code](https://github.com/caitlingrasso/empowered-nca) | Empowerment-related time-lagged action–sensor MI as auxiliary evolutionary objective improves target-shape morphogenesis. | Very high for empowerment + growing cells. | No prescribed final silhouette, continuous `r_e(t)`, transport task, finite material budget and predictive retention. Read [updated arXiv methods](https://arxiv.org/abs/2205.06771) when reproducing. |
| Grasso & Bongard (2023), **Selection for short-term empowerment...**; [DOI](https://doi.org/10.1145/3583131.3590469), [code](https://github.com/caitlingrasso/empowered-nca-II) | Short-horizon empowerment-related MI improves shape maintenance and adaptation in tested NCA; minimizing local action entropy **sometimes outperformed their empowerment metric**. | Very high and provides essential ablation. | Include minimum local action entropy, pure shape reward, and equivalent-computation baselines. Do not imply their empirical lagged MI equals full action-channel capacity. |
| Scholl, Rule & Hennig (2021), **Information theory of developmental pruning**; [DOI](https://doi.org/10.1371/journal.pcbi.1009458), [code](https://github.com/carolinscholl/PruningBMs) | Information-importance/Fisher-based local pruning of abstract neural connections outperforms magnitude-only pruning in their model. | High for “information-guided pruning” per se. | Continuous material geometry and physical transmission, plus prediction of dormant branch value under changing demands. |
| Kramar & Alim (2021), **Encoding memory in tube diameter hierarchy...**; [DOI](https://doi.org/10.1073/pnas.2007815118) | Experiments and model show food-location memory imprinted in living Physarum tube diameters. | Extremely high for `r_e(t)` as physical memory. | Quantify *predictive* information vs mere retained past and test cost-adjusted contribution under future changes. |
| Schick, Kramar & Alim (2024), **Dynamic Cost Allocation...**; [DOI](https://doi.org/10.1103/PRXLife.2.033005) | Three Physarum morphodynamic modes; estimates building-related elastic and viscous transport costs, with different search behavior. | Extremely high for “growth vs transport energy trade-off.” | Full accounting for proposed synthetic model; task-specific **value of future options** versus local flow allocation. The 2024 “building cost” is a model-based elastic-energy proxy, **not direct synthesis/metabolic calorimetry**. |
| Schick et al. (2026), **Decision-Making in Light-Trapped Slime Molds...**; [PRX Life DOI](https://doi.org/10.1103/rv7g-d9kx) | Exploration and eventual escape in optical trap correlate with evolving peristaltic contraction modes. | High for distributed exploration/selection. | Distinguish coordination caused by mechanics alone from any additional need for predictive or empowerment-based controllers. |
| Le Verge-Serandour et al. (2026 preprint), **Edge Pruning by Relative Resistance**; [arXiv:2610.13582](https://arxiv.org/abs/2610.13582) | Local-flow-accessible/topology-conditioned relative-resistance pruning yields efficient sparse loop-rich networks. | Essential benchmark for non-information-theoretic pruning. | Test its robustness to *temporally correlated, unseen future demand*, keeping physical budgets matched. Uploaded preprint (October 9) is not a peer-reviewed final publication as of this date. |
| Bouvard et al. (2025), **Self-Organized Homogenization of Flow Networks**; [PRX DOI](https://doi.org/10.1103/j5ch-4vkh) | Hydrogel microfluidic channels physically widen by flow-carried enzyme erosion and homogenize resistance under specified conditions. | High for self-adaptive channel radius by actual matter dynamics. | Separate reaction/flow-driven passive geometry adaptation from information-guided growth; test whether extra sensing pays. |
| Chen & Prokopenko (2025), **Why collective behaviours self-organise to criticality...**; [arXiv:2409.15668](https://arxiv.org/abs/2409.15668) | In an Ising framework compare predictive information, empowerment, active inference and thermodynamic efficiency; preferred regimes differ. | Very high for “compare all objectives”; they can **conflict**, so a combined scalar objective is not justified automatically. | Coupled roles **with changes to physical network architecture** and independently recorded costs; don't assume a single critical optimum. |
| Meng, Piazza, Both, Barzel & Barabási (2026), **Surface optimization governs the local design of physical networks**; [Nature DOI](https://doi.org/10.1038/s41586-025-09784-4) | Measurements across physical networks show length/cylindrical-volume minimization cannot fully predict branching geometry; surface cost provides stronger predictions. | Strong novelty constraint for geometrical material economics. | Account for occupied **area, volume, interfaces, branch junctions**; prove results are robust beyond disk discretization artifacts. |
| Massari et al. (2021), **Experimental Evidence that Empowerment May Drive Exploration...**; [arXiv:2107.07031](https://arxiv.org/abs/2107.07031) | Empowerment-inspired exploration agent uses a variational information bottleneck to predict future sensor state. | Further direct IB–empowerment proximity. | Actual structural construction/retention and physical resource competition. |
| Tero et al. (2010), **Rules for Biologically Inspired Adaptive Network Design**; [Science DOI](https://doi.org/10.1126/science.1177894) | Physarum-like flow reinforcement reproduces a transport-network design trade-off. | Essential classic baseline. | Predictive/controllability-informed growth must beat this or explain its failure domain. |

**Novelty verdict, bounded by this targeted search:** The *individual ingredients and many pairwise combinations* are well established, including empowerment+growth, empowerment+IB, memory in tube diameter, information-based pruning, and building/transport trade-offs. We found **no clearly identified paper performing this exact four-function causal decomposition as a physical, resource-limited, time-varying network with full ablations and the Still-type predictive-dissipation connection derived for its actual dynamics**. This is **an unverified possible gap**, not a proof of priority or a guarantee of publishability. Broader searches (active inference, developmental robotics, plant vascular optimization, resource-aware graph rewiring, causal empowerment, network control) and feedback from Fable/Susanne are still needed.

## B. An implementation-level formulation (not a derived law of nature)

Choose a network state at time `t`:

`B_t = (G_t, {r_e(t)}, {x_tip(t)}, {z_e(t)}, R_t)`.

- `G_t=(V_t,E_t)`: connected skeleton (graph). Its edges **carry finite material geometry**, not free mathematical wires.
- `r_e(t)`: radius/width along each edge (can be sampled at structurally attached disks `r_i(t)`). Graph vertices record branch contacts and tip positions; self-intersections and merges must be explicit.
- `x_tip(t)`: growing frontier, trial protrusion positions, provisional contacts.
- `z_e(t)`: local optional predictive/relevance signals; if persistent, **count their memory and controller costs**.
- `R_t`: available material/energy budget.
- `O_t`: sensory data on resources/flow/pressure/targets. `D_t`: future external demand. `A_t`: available structural or transmission actions.

**Causal sequence at a time step (research proposal):**

1. **Observe and update predictive memory**: `S_{t+1} ~ q(S_{t+1}|S_t,O_t)` via an independently validated DIB estimator; evaluate `q(D_{t+1:t+H}|S_{t+1})` on *held-out* sequences.
2. **Generate exploration candidates** from local stochastic growth and geometry, with an optional causal empowerment contribution. Before any contact, `a_t` might be extend tip, rotate, branch, or stop. Do not equate randomness with empowerment.
3. **Estimate candidate value** using possible reachability/control and task success over a stated horizon. Proposed quantity: `ΔEmp_e = Emp(B_t + e) - Emp(B_t)` where the action-to-sensor channel and observational capability are explicit. True empowerment maximizes over admissible action distributions; sampled lagged MI is only a surrogate.
4. **Select and reinforce relevant pathways** using IB (a compression/relevance trade-off on a specified variable) or a separate counterfactual task criterion. Define, for local edge history `H_e`, relevance target `Y`, representation `Z_e`, `L_IB = I(Z_e;H_e) - β I(Z_e;Y)`. **A positive value of I is not a command to thicken an edge.** A distinct, predeclared learning rule maps predictive contribution to material allocation.
5. **Use predictive memory to evaluate retention**, including a provisionally idle connection that may be needed at `t+H`. DIB may or may not improve retention relative to an optimal Markov predictor or a simple schedule model.
6. **Accept/reject under real resource accounting**. This is *not* merely a fourth additive reward: material conservation and allowable actions are constraints. Record net work/heat/material in appropriate units or explicitly named model proxies.
7. **Transmit through current network**, measure task performance, observe changed flows, repeat. Document closed-loop feedback introduced when structure changes demands or observations.

**Important distinction:** Ordinary spatial IB and DIB are not inherently separate physical processes. Their different roles here are an *engineering decomposition to test*, not a claim that biological tissues implement two algorithms. Spatial `Z_e` and temporal `S_t` must be distinguishable; alternatively compare with a *single* predictive IB memory to see whether the two-component architecture adds value.

### Formal objective/constraints (proposed only)

For comparison, define loss as an actual task measure, e.g. unmet flow or delivery delay, not just a visualization metric:

\[
J_{0:H} = \mathbb E\left[\sum_{t=0}^{H} \ell_{\rm task}(D_t,\widehat D_t, Q_t)
+\lambda_E E_{\rm op}(t)+\lambda_R E_{\rm remod}(t)\right],
\]

subject to `V(B_t) ≤ V_max`, power/actuation bounds, positive radii and valid geometry. For experiments with fully monetized resource costs, use compatible units; information in bits gets **explicit** coefficients and is not itself joules.

A separate lifetime accounting ledger is

\[
E_{\rm life} = E_{\rm construct} + E_{\rm transport} + E_{\rm maintain} + E_{\rm remodel} + E_{\rm sensing} + E_{\rm controller}.
\]

**For idealized cylindrical fluid channels only**, `R_e=8 μ L_e/(π r_e^4)` and `P_{\rm diss,e}=R_e Q_e^2`; volume `V_e=π r_e^2 L_e`. Branch junctions, changing cross section, fluid inertia, and Physarum contractile pumping require other physics. For our toy graph, clearly label any `R_e∝L_e/r_e^α` as *chosen proxy*, not a measured Physarum or hard-disk law. The union volume/area of overlapping disks must avoid double counting, and outcomes must be insensitive to discretization resolution. Growth and repair work are **not** thermodynamically recovered simply because a branch is removed.

**Still linkage, precisely:** In the uploaded DIB draft the learning machine has update `p(S_{t+1}|S_t,X_t)` and a dissipation lower bound involving information captured about current/past inputs minus information predictive of the next observation. The underlying analysis is a **no-feedback learner**. Here the growing network changes its future sensing, motion, pressure and geometry; the dynamics become a feedback loop. Therefore **do not write** `E_construct = k_B T I(...)` or claim Still's existing lower bound is a construction-energy law. A legitimate long-term theoretical project is specifying the physical driving protocol and deriving the correct feedback-aware entropy production or information–work relation for the changing morphology.

## C. Tiny benchmark with a derived, falsifiable retain/prune threshold

**Smallest model:** input with two alternative paths to destinations A/B, adjustable path radii and simple transmission delays. Demands `D_t∈{A,B}` have *equal stationary frequencies* but different transition probabilities/autocorrelation. Fix a two-state Markov model, a renewal/switch cost, controller observation accuracy and a time horizon. Allow keeping, thinning, removing, or rebuilding a currently idle branch. Build a clean numerical test before growth cones or 3D rendering.

**One-step threshold (special assumptions):** Keep a dormant B connection for the next period if its upkeep `h` is less than the *predicted expected avoided cost* of needing to rebuild B. Suppose under a one-period model the conditional probability B becomes needed is `p_t`, rebuilding costs `b`, and delayed service costs `d`, with no other benefits/costs. Then

\[
\textbf{keep B iff }h < p_t(b+d).\qquad (\text{toy criterion, NOT a universal law})
\]

Example `h=3` resource units, `b=8`, `d=12` gives threshold `p_t>3/20=0.15`: if predicted `p_t=0.20`, keeping saves expected `4>3`; if `p_t=0.05`, pruning is preferable `1<3`. **These are model units**, not measured joules. For multiple periods, this is a Bellman/partially observed control problem; sunk costs, finite construction time, continuous radius, discounting, independent pressure sources and limited material invalidate simple threshold reuse. The true test asks whether a **DIB-compressed predictor** estimates relevant probabilities well enough to change decisions after charging for its computation.

**Key experimental manipulation:** Hold stationary demand frequencies and mean total requested transport fixed. Vary *temporal correlations*, volatility, growth/maintenance price, ability to observe changes, and perturbation rate. Predict an empirical boundary between **prune** and **retain as a preparedness investment**. A useful result is a boundary predicted before simulations on held-out environments, with error bars and clear failure domains. Different switching laws may produce multiple boundaries.

**Physical memory test:** In half the runs a digital controller keeps predictive state `S_t`; in the other half **reset controller state after each update** and permit only `(G_t,r_t)` to persist. Measure `I(B_t;D_{t+H}|O_t)` with held-out estimation and shuffle/history-destroying controls. Test counterfactual impairment by shuffling branch radii among matched networks or resetting radius history at identical topology. High MI alone does not imply useful memory: require reduction in held-out task loss or lifetime expenditure.

## D. Causal ablations: distinguish genuine synergy from component coexistence

**Always hold the physical budget fixed across competitive policies.** Ideally run a factorial `2×2×2` experiment toggling:

- E: empowerment-derived candidate score vs random/exactly matched candidate proposal.
- I: IB-based selective reinforcement vs local-flow reinforcement or matched marginal-utility heuristic.
- D: DIB predictive memory vs a memoryless state or simple Markov / Kalman predictor of comparable capacity.

A completely separate unconstrained diagnostic may be informative, but removing the budget from only one policy is an invalid comparison. Compare eight combinations plus **preexisting strong baselines**: Tero-type flow reinforcement, Alim-style relative-resistance pruning, Fisher-information-inspired pruning, shortest-path or Steiner-like material economy, simple min-entropy signaling heuristic (from Grasso), and oracle planner **clearly marked as an upper bound**. Test a single integrated predictive-IB component in place of separate spatial-IB+DIB to avoid a contrived decomposition.

For every algorithm use the same observation channels, permissible actions, random seeds, train/test splits, material budget, and—where possible—matched controller compute/energy. Distinguish intrinsic objective as a **training loss** versus a **within-lifetime action policy**. Hold out new demand sequences and unseen parameters; report mean, variability, uncertainty and compute cost. Include policies where you delete each module and a one-factor-at-a-time sweep. Synergy is a measured **interaction** (joint advantage beyond the effects of individual toggles), not the fact that the four modules can run together.

Success criteria (set before running): (i) full policy improves the *task-performance vs total-lifetime-cost Pareto frontier* on held-out regimes; (ii) advantage persists after ablation of fancy metrics versus fair baselines; (iii) the estimated retain/prune boundary changes predictably with demand autocorrelation and construction cost; (iv) benefit survives geometry-resolution tests; (v) mapping to physical apparatus remains consistent with its energy ledger. A null result showing local rules match information methods also counts as a defensible *scientific finding*, if the tests were well-powered.

## E. Engineering path, independent of ongoing HardDisks papers

1. **R0 (~bounded side experiment):** 2D graph, 2 paths, scalar radius per path, limited material budget, Markov/non-Markov demand, baseline pruning/growth, reproducible logged traces, and numerical cost threshold. No 3D, no biological claims.
2. **R1:** Swap in your nearly finished DIB implementation as an optional predictor, compare to exact Markov and simpler predictors; tune on training sequences, score held-out.
3. **R2:** Add candidate growth, soft disk-chain geometry (`x_i,r_i(t)`), provisional contacts, branching and retraction. Separate structural disks from deterministic elastic thermal EDMD particles. Export graph, swept surfaces, and junction-aware resource costs; compare with DLA point-line baseline.
4. **R3:** Physarum-inspired broad front, contraction/flow and flow-dependent thickening; optionally neuron-like growth cones. *Same comparison interface, different underlying laws.* Validate against actual biological shape/time data before analogical interpretation.
5. **R4 (much later):** Replace idealized transmission with measured response kernels or emulators from **Paper 1 Román physics** and **Paper 2 driven multi-wall energy-transfer**. Only add full EDMD per edge if the improved fidelity changes a measurable conclusion. Do not turn current energy-conserving EDMD into a growing-particle simulation without revisiting its dynamical laws and accounting.

**PhD scope decision remains unchanged:** finish Paper 1, Paper 2, and separately validate DIB first. Time-box the adaptive-network benchmark as a parallel research branch within IdeaLabs; do not make a giant visual simulator a prerequisite for publication. IdeaLabs is the research knowledge/project graph; DV8/NeurAIpil are separate projects/use cases, not the simulator's core architecture.

## F. Assumptions that could invalidate the story

1. **Empowerment may favor no growth.** Can a larger branch increase influence enough to offset loss of power elsewhere? Define sensors/actions and causal perturbations.
2. **IB is a coding principle, not a synapse-growth law.** Explicitly test at least one local mapping from learned relevance to radius; do not infer it from the IB optimum.
3. **DIB may offer no advantage.** Memoryless adaptation or a small explicit Markov predictor could already be optimal under our simple switching environment.
4. **Prediction != physical energy savings.** A predictive model's information bound does not determine biomass synthesis/actuation heat.
5. **Not all structural memory is useful.** Geometry may retain historical correlations yet harm decisions under nonstationarity.
6. **Different functions may fight.** An empowerment proposal may conflict with IB selective reinforcement; a resource constraint may block both. Testing coordination/frustration and phase changes is as important as successes.
7. **Unclear lengthscale laws.** Disk-count is resolution-dependent; radius, union area and 3D smooth-junction/surface cost are physical quantities only under specified geometry.
8. **Closed-loop failure.** Learning state, growth and environment are mutually coupled. A no-feedback predictive-dissipation bound must not be transplanted wholesale.
9. **Active inference / morphological computation prior art.** A deeper search might reveal a near-identical integrated architecture under different terminology. Explicitly ask reviewers to find counterexamples.
10. **Time and opportunity cost.** If the experiment cannot deliver a pre-registered reproducible advantage or illuminating null in a bounded period, archive it and protect the validated HardDisks papers.

## G. Specific questions for Fable / Susanne / Astra / Opus

- Find the closest prior art where **empowerment selects new physical pathways**, **IB or equivalent predictive relevance strengthens connections**, **dynamical predictive state prevents costly pruning**, and **physical construction/maintenance energy constrains all three in one experiment**. Return exact section/equations/DOIs; do not claim novelty without one-to-one comparison.
- Is there a principled derivation from Still's thermodynamics of prediction to **feedback-controlled morphological remodeling**? If not, state the missing terms and what protocol would make such a derivation possible.
- Can a *single* predictive IB state replace our two distinct IB/DIB roles? When does a split architecture offer an identifiably different prediction?
- Can an open-loop empowerment measure meaningfully select physical growth if the system cannot sense the consequences of its new branch? Define a genuine causal channel.
- Derive the optimal retain/prune policy in a small partially observed two-path resource environment; identify exact thresholds and conditions under which simple local policies provably suffice.
- Design the strongest fair baseline: 2026 relative-resistance pruning + 2024 Physarum flow/cost model + Grasso NCA/homeostasis + a minimal Markov/Bellman planner. Require a test where an extra information module passes a **predeclared** improvement criterion.
- Do we need structural disk chains for the first theorem/benchmark? If not, begin graph-only. If yes, show a geometry-dependent result unavailable to graph-radii representation.

## H. Search scope, provenance, and citations

**Search performed (2026-10-13):** targeted English-language web searches pairing empowerment, IB, morphogenesis, energy/material budget, predictive control, structural plasticity, and growth; investigation of direct papers, abstracts, DOIs and code links; grounded reading of user-uploaded PDFs `Empowered Neural Cellular Automata.pdf`, `Selection for short-term empowerment...pdf`, `The information theory of developmental pruning...pdf`, `DynamicIB_2022-2.pdf`, `Why collective behaviours self-organise to criticality...pdf`, `Dynamic Cost Allocation...pdf`, `Decision-Making in Light-Trapped Slime Molds...pdf`, `Efficient Transport Networks through Edge Pruning by Relative Resistance...pdf`, `Self-Organized Homogenization of Flow Networks.pdf`, and `Flow-Induced Vascular Remodeling on-Chip...pdf`. The older full research notebook §§0–43, all associated references and the five-paper trajectory remain preserved.

**Research limits:** search is not exhaustive, DOI/paper coverage may differ by database, citation chains can conceal near-identical models, and no formal systematic-review protocol or code replication was conducted. Specific **absence claims are not established**. All proposed equations above are **our modeling choices**, except standard definitions/equations cited to original works. The separate article sources listed above should be treated as the authoritative primary references. This is a hypothesis and experiment blueprint, not completed research.

**One-sentence north star:** **How does a physical network learn *where to explore*, *what to reinforce*, *what to remember for tomorrow*, and *what material and energy it can afford*—and when does this causal specialization outperform simpler rules?**


**Revision log — v1.8 (2026-10-13):** Added corrected coupled-role hypothesis, newly identified Klyubin/Polani empowerment morphology and IB action-bandwidth prior art, Nature (2026) junction geometry constraints, exact mini-benchmark and causal ablations; preserved all existing five-paper material and previous research discussions.


---

## 45. [2026-10-15 / v1.9] User clarification — **ONE morphology engine + actual EDMD transport inside growing channels**

The preceding discussion needs this explicit implementation-level refinement. Structural **large disks generate a *lumen*** rather than themselves acting as the freely transported collision particles. A connected graph controls topology, smoothly tapered swept 2D channels provide an open interior, and a distinct population of **small EDMD mobile disks** moves through it. Different growth policies control the same structural engine: pulsating Physarum fronts with active peristaltic work, growth-cone neurite exploration/contacts, DLA, and later hypothesis-driven E/IB/DIB. Frozen channels and a fixed Y junction are validation fixtures, not a replacement for free morphogenesis. Reuse measured real physics from HardDisks Papers 1–2 only after their independent validation; do not change their execution protocol.

The user also suggested **delayed pressure-triggered growth**: when an internal spring reaches its mechanical stop, sustained excess pressure may trigger widening or branching instead of additional spring compression. This is a speculative mechanically mediated control law, **not a claim that springs generate growth energy**. Explicitly account for continuing piston/pump work, wall reaction forces, resource supply, material and growth costs. This is analogous to stress sensing rather than proof of biological Physarum mechanics.

**Future—not MVP:** coupling network radius to experimentally characterized hydraulic conductance and pressure response; passive electrical cable and active action-potential propagation; explicit charged mobile ions, selective membranes and maintained ion gradients; state-dependent/memristive conductance; Maxwell-demon-inspired feedback with measurement/reset work. Increasing neutral EDMD particle density is NOT an electrical model. The current engine is a physical model of mechanical and material transport, then later an extensible cross-physics framework after independent model validation.

### 45.1 Canonical engineering handoff (read in full)

**`IDEAS/MorphoTransport_Network_Engine_Blueprint_2026-10-15.md`** is the new technical source of truth for this plan. It specifies exact state variables, swept-disk geometry, physics conventions, open/closed EDMD boundaries, growth policy interfaces, delayed pressure growth, staged validation M0–M7, module layout, falsification gates and copy-ready Fable/Opus instructions. Its equations are explicitly modeling proposals; no simulator or biology claim has been validated by this update.

### 45.2 Division of labor and project boundaries

Fable: one scientific critique of the frozen brief. Opus/Claude Code: new **independent MorphoTransport repo** implementing a minimal, fully tested M1 disk-graph and channel geometry, then M2 EDMD; no full simulator generation in one pass. NeurAIpil A3 Sol/Opus heterogeneous deliberation: optional bounded review experiment using frozen M0 spec, **not a required implementation dependency**. IdeaLabs: store ideas/briefs/provenance, not the simulator engine; DV8 remains separate. HardDisks Paper 1, Paper 2 constant-speed `x_p(t)=x_0+vt`, work/spring-energy ledger and DIB validation retain priority. Provider credits (~$250 reportedly expiring November 4) do not justify architecture churn or unbounded multi-agent calls; verify account-specific restrictions, prioritize measurable tests and preserve capacity for physical papers.

### 45.3 Implementation blueprint — complete embedded handoff


> **User's precise request, 2026-10-15:** Grow a network from connected **large structural disks**, derive a **continuous tapered 2D channel boundary** (a swept-disk union / lumen, NOT a convex hull that fills gaps between branches), and run **smaller collision-resolved hard disks through the interior**. Observe particle/energy/pressure/signal transport; allow slower structural growth, retraction, diameter adaptation and topology changes in feedback. Use a **single engine with pluggable growth policies**, not a different simulator for every organism. In the future investigate physically modeled conductance, membrane ion channels, action-potential propagation, history-dependent/memristive connections, and Maxwell-demon-type sensing/feedback costs. **These are later research options.** Preserve the current HardDisks experiments as-is.

## 0. Hard gate: do not delay or alter the existing research

**Paper 1:** Finish Román two-compartment static/equilibrium stiffness, adiabatic/acoustic response, sound-speed validation, and unresolved finite-size/FFT checks. 

**Paper 2 (IMMEDIATE):** Existing `hspist3`: right-hand **prescribed constant-speed piston** `x_p(t)=x_0+vt` over a declared window; multiple impermeable movable dividers; a spring-backed movable terminal wall on the **left**. Measure input force/work `W_in = ∫ F_p \dot x_p dt`, output spring displacement `x_s(t)` and energy `ΔU_s = 1/2 k[x_s(t)^2 - x_s(0)^2]`, and the full work–energy ledger, with `ΔU_s/W_in` reported only for clear positive input work and fixed endpoints. Mobile gas disks stay within their compartments while transmitting momentum to walls. **No network growth code, ion dynamics or rewritten event scheduler enters this validation branch.** DIB implementation and its no-feedback dissipation derivation remain an independent third effort.

## 1. Why a single engine, not rigid Y pieces or circles alone?

Three **distinct, composable state layers**:

1. **Morphological state**: Graph `G=(V,E)` (parents, growth tips, junctions, loops, provisional contacts), disk/tip elements `D_i=(x_i,y_i,r_i,type_i,age_i)`, and edge centerlines `c_e(s,t)` with continuous local lumen half-width `r_e(s,t)`. **Disk chains drive morphology; centerline graph owns connectivity; smooth swept disks/segments generate collision boundaries.** Tapered slices and parametric Y-junctions are geometry tools/test fixtures, not biological building blocks that constrain topology. Support both growth from large soma and sheet-like exploratory fan.
2. **Transport state**: small mobile particles `p_k=(x_k,v_k,a_k,m_k,tag_k)` inside `Ω(t)`, and optional walls/pistons/springs/reservoirs. Initially **uncharged 2D elastic disks**, collision resolved with appropriate static-wall and particle events. Do not call them literal sugar molecules or blood cells without matching their physics. Open channels, reservoirs and collection outlets are separate from sealed compartments; decide which mode is being run and define boundary exchanges.
3. **Learning/control state**: `state_predictor` (possibly existing DIB), `growth_policy`, `resource_ledger`, `signals`, `measurements`. The four **coupled, not interchangeable** hypothesized roles are empowerment-related proposal of exploratory actions, IB-inspired relevance-based reinforcement, DIB-informed predictive retention of underused branches, and physical budget/thermodynamic constraints. None of these assignments follows automatically from an information-theoretic definition; test against simpler flow rules and through ablations.

**2D collision domain:** `Ω_channel(t) = ∪_{e∈E} ∪_{s∈[0,1]} B(c_e(s,t),r_e(s,t))`, with continuous unions at junctions and loops. For transport disks of radius `a`, the center-accessible region is `Ω_channel ⊖ B_a` (morphological erosion). Avoid building a global convex hull or charging energetic cost by the number of numerical circles. Treat circles as lumen-support discs by *explicit convention*; if instead disks represent solid tissue, calculate the complementary free region accordingly. Store a single unambiguous geometry convention in run metadata.

**Material cost:** occupied wall/tissue area or 3D volume (when later using realistic tubes), appropriate interface/surface costs, creation/disposal/recycling and junction volumes; count overlaps once. A 2D channel's open area is not itself the solid material consumed. At fixed numerical geometry, mesh resolution must not change material cost. A lumen of width `w` must be appreciably wider than `2a` for mobile hard-disk centers to pass; model bottlenecks and occlusion explicitly.

## 2. Swappable growth/actuation policies, shared measurement interface

`growth_policy.propose(observations, morphology, prediction, budget, rng) -> proposed_actions` where possible actions are `extend`, `bifurcate`, `merge_if_physically_valid`, `widen`, `narrow`, `stabilize`, `retract`, `prune`, `oscillate_wall`, `no_op`. Resolve legality, material reservation and topology changes in a **separate executor**, not independently inside policies. Observation availability and sensing costs must be consistent between experiments.

| Growth mode | Candidate mechanism; **hypothesis/approximation** | Crucial checks |
|---|---|---|
| `Neutral` / manual | fixed graph, deterministic taper, editable Y/T junction, straight/curved channels | geometry/EDMD benchmark ground truth |
| `DLA` | random walker attachment to structural tips | mainly tree structures; does not autonomously learn or produce loops |
| `Physarum` | broad exploratory fan, moving tips, local flux-dependent tube growth/pruning, distributed rhythmic peristaltic contractions and delayed mechanical feedback | needs active actuation/work, transport and wave phase; oscillations do **not** imply empowerment |
| `Neurite` | soma, branching axon/dendrite growth cones, filopodia-like transient probes, chemo/activity guidance, provisional contacts and separate stabilization/pruning | touching is not automatically a functional synapse; action potential is not growth |
| `E_I_D` candidate | E proposes; IB relevancy reallocates radius/connection strength; DIB temporal prediction retains/reuses dormant paths; budget approves/denies | matched fixed budget and observables; factorial 2³ ablation, simpler baseline can win |

**Physarum pulse mode:** `a_e(t) = a_{0,e} + δa_e sin(ω_e t + φ_e)` is a simple prescribed actuation *test*, not a demonstrated biological law. Evolving phases, internal wave coordination and resource transport require explicit dynamics. Active contraction costs mechanical work; no free propulsion from passive equilibrium collisions. Allow slow fan expansion on a distinct timescale from contractions. Geometry should provide wide sheets near active leading fronts transitioning to thin tubes. A single centerline approximation may be unsuitable for a fan until skeletonization.

**User's delayed pressure-triggered growth thought:** A spring-backed end wall approaches maximum compression `x_s→x_max`. Beyond its stop the boundary reaction/pressure may rise; use **measured stress/pressure over a dwell interval** to trigger *candidate* growth rather than let a spring pass through the stop. E.g. `q_e(t)= ∫_{t-τ}^t max[P_e(s)-P_*,0]ds`; propose widening/branch creation only if `q_e>Θ`, then allocate material and growth energy. Add hysteresis and cooldown to avoid chattering. A fully compressed spring is not an energy source: the new wall displacement/growth draws energy from continuing external piston work or an explicitly stored chemical/material reservoir. Track all relevant pressure-work, elastic, growth and output terms. This is a **testable mechano-growth hypothesis**, not claimed Physarum physiology.

## 3. Transport and physically meaningful observables

Use EDMD for **elastic collision transport at known particle radius, density and wall geometry**; not automatically viscous cytoplasmic streaming, blood, or ionic conduction. Validate a fixed reflective corridor and an analytic Y branching shape before running growing networks.

- Particle throughput: `J_j(t)≈[N_j^{out}(t+Δt)-N_j^{out}(t)]/Δt` per outlet. Track cumulative arrivals, path choice, time-of-flight/residence-time distribution, lost/trapped particles, collision count, packing fraction and throughput uncertainty.
- Mechanical response: force/impulse flux on walls, local particle momentum, pressure estimates with clearly defined estimator and sampling windows; propagating acoustic/wall modes when driven. For input `U_t` (e.g. randomized piston impulses), estimate `I(U_t;Y_j(t+τ))` with held-out sequences, known causal drive, time-resolved estimators and null controls. MI is **not** velocity or causal proof by itself.
- Efficiency: input piston/pump work, change in internal kinetic + potential + spring energies, boundary particle energy fluxes, geometry motion work, external baths if any, construction/maintenance, controller/sensing expenditure (or explicitly declared proxy). Separate transport throughput, signal fidelity and energetic efficiency; high flux needn't imply high information.
- Closed EDMD compartments can mechanically communicate through impermeable movable dividers **without particle mixing**. The eventual growing Y network uses **open** branch junctions permitting mixing and throughput. Test these as separate protocols.
- Density sweep: high packing can alter collision frequency, mean free path, compressibility, wave propagation and ordering; **not automatically electrical signaling**, and may clog/slow throughput. Probe where simple resistor/Poiseuille graph approximations fail, with appropriate 2D geometry and regime qualifications.

For future geometry updates, start with **frozen wall geometry for each transport window** and explicit pause/validation before updates. Any particle displaced or trapped by remeshing requires a deterministic policy and mass/work accounting. No invalid overlaps, teleportation, or unexplained thermalization. Later implement truly moving-wall event prediction and pressure work only after frozen-boundary EDMD tests pass.

## 4. Future electrical/membrane/memristive modes — PARKED, not EDMD MVP

**Conductance baseline:** A graph-edge transport law `Q_e=ΔP_e/R_e` (hydraulic analog), whose effective `R_e(w,length,ρ,particle_size,forcing)` is **measured from EDMD** or benchmarked with continuum theory. Resistance of 3D cylindrical Poiseuille tubes scales ∝`L/r^4`; slit-like laminar 2D flow per depth scales ∝`L/w^3` under different assumptions. Neither law is universally valid for dense confined 2D hard disks. Fitted `R_e` is **hydraulic**, not automatically electrical resistance.

**Electrical cable mode:** Each edge may later hold membrane potential `V(s,t)`, axial resistance, membrane capacitance and leak/ionic current; a suitable baseline is `c_m ∂_t V = ∂_s(g_a∂_s V) - i_ion(V,gates,[ions]) + i_ext` with geometry-dependent coefficients and current-conserving junctions. A passive cable does not produce actively regenerated action potentials; add excitable gates (e.g. FitzHugh–Nagumo test, then carefully parameterized Hodgkin–Huxley / ion channels), concentration reservoirs and pump work. Thick branches may change conductance and capacitance but do not imply improved neural function.

**Mobile ion-disks:** Require charge `q_i`, field/force/Poisson or electrodiffusion treatment, thermal noise/dissipation if modeled, selective membrane/channel boundary conditions, concentration gradients, charge conservation and electrochemical free-energy ledger. Plain neutral elastic EDMD at very high packing **cannot turn into electricity by density alone**. Force-driven charged dynamics invalidate an unmodified exact ballistic EDMD event scheduler between collisions. First use a separate cable/electrodiffusion module rather than silently modify validated hard disks.

**Memristive links:** A proposed history-dependent conductive state `g_e=g(z_e)`, `\dot z_e=f(V_e,I_e,[ions],z_e)` can represent state-dependent ionic channels or plasticity, with an explicitly checked energy/passivity model, but **hysteresis alone does not prove a physical memristor implementation**. Distinguish biological synaptic plasticity, ion-channel gating, and electronic memristors.

**Maxwell's demon / information engine:** Possible future sensing-feedback experiments on particle transport or channel gates should log measurement, controller, feedback work, information acquired, and reset/erasure costs under an explicit stochastic-thermodynamic model. No spontaneous work generation from information alone; do not transplant Still's no-feedback DIB bound to an interactive growing system without rederivation.

## 5. Exact implementation order and acceptance gates

**M0: Scientific design freeze (few days, bounded).** Independent repository `MorphoTransport` (working name), not a fork of `hspist3`; archive a minimal numerical protocol, data units, geometry convention, correctness gates and baseline sources. Require reviewer to challenge before significant coding. Code should be architected as a reusable simulation engine with policies as plug-ins, not one graph per organism.

**M1: Disk chain -> graph -> smooth collision hull (NO particles).** Spawn/extend/retract/taper/bifurcate/connect disks; validate Y and one cycle; deterministic seeded replay; export polygonal or analytic curved walls; verify no gaps or illegal self-overlaps, stable junction geometry and material conservation. A prepared tapered-slice primitive is an interpolation implementation detail, not the biological growth source of truth. View in simple 2D before 3D.

**M2: Frozen channels + small particle EDMD.** Start corridor, closed box, fixed Y with one inlet/outlet and reflective mode; test exact sphere–sphere and sphere–straight-wall collision times, energy/momentum checks for closed stationary walls, throughput count and no leakage. Only after those pass add curved walls and disk-generated hull. Reuse `hspist3` as **reference and optional adapter** after confirming its assumptions; do not make experimental additions to Paper 2 core.

**M3: Reproducible open-boundary transport.** Piston/reservoir driving with declared injection/removal laws and boundary work; flux per outlet, transit-time distribution, pressure/impulse probes; compare identical geometry with simple graph conductance and measure discrepancies. Provide CSV/Parquet/JSON run manifest with RNG seed, code revision and parameter units, along with rendered trace/video as optional output.

**M4: Slow feedback geometry.** Short EDMD window at frozen geometry -> measure -> call policy -> reserve material/energy -> commit geometry -> validate/reconcile particle placement -> continue. Compare a no-growth baseline; test 1) thickness changes, 2) branch retraction, 3) new branching, 4) pressure+dwell-triggered growth. If moving-wall work remains uncalculated, declare the quasi-static geometry approximation and do not claim full thermodynamic closure.

**M5: Biology-inspired profiles.** Physarum rhythmic actuation and flow-driven adaptation with explicit work and timing; separate Neurite growth-cone mode; DLA control. Same observation/measurement/ledger interfaces. Prioritize scientifically interpretable behaviors; no claim of reproducing living-cell biochemistry.

**M6: Coupled information policies and hard scientific test.** E/IB/DIB roles from research notebook; 2-path predictable vs random demands with matched stationary distributions, build/maintain costs and true held-out trajectories; flow-local, relative-resistance, random, shortest-path, simple Markov, and oracle upper-bound baselines. Score lifetime performance vs total budget, predictive structural MI conditional on observations, resilience, signal fidelity, and complexity. Test independent ablations and interaction. Information per joule requires clear work protocol and separate data-estimation uncertainty.

**M7: LATER:** calibrated mechanical transfer elements from the established Paper 1/2 results, electrical cable mode, ion-selective particle dynamics, ion pumps/membrane potentials, history-dependent conductance, feedback demon and possible 3D. Each is a standalone validation milestone, not simultaneous MVP.

**Predeclared failure tests:** (1) fixed-channel mass and energy checks fail; (2) geometry discretization changes material budget or measured flow; (3) remeshing moves particles without accounting; (4) growth without external energy source; (5) prediction/empowerment gains vanish against matched simple controls; (6) visual similarity to biology substitutes for a measured agreement. Any is reason to stop and debug rather than add another architecture module.

## 6. Suggested software boundary / layout (proposal)

```text
MorphoTransport/
  README.md                     # purpose, scientific limits, Quickstart
  docs/SPEC_M0.md               # sign-off review & invariants
  docs/PHYSICS_LEDGER.md        # dimensions, open/closed boundaries, work
  src/geometry/                # disks, centerline graph, swept channels, junctions
  src/growth/                  # DLA, Physarum, neurite, neutral, E/IB/DIB adapters
  src/transport/               # EDMD adapter/kernel (frozen geometry first)
  src/coupling/                # observations, budget, delayed update executor
  src/measurements/            # flux, pressure, transit time, info and energy
  src/protocols/               # Y tests, piston, demand switching, seeded configs
  src/viewer/                  # thin 2D inspector (do not build 3D first)
  tests/geometry/ tests/edmd/ tests/budget/ tests/replay/
  experiments/README.md       # immutable run manifests; no massive files in Git
```

Practical language choice is **not yet prescribed**: Python for initial geometry, scientific orchestration and plotting is convenient, with compiled EDMD kernel when performance or deterministic event ordering requires it; a small web/SDL renderer is separate. Require a clean interface, not an unnecessary wholesale rewrite of your existing C/SDL2 hard-disk simulator. Use a known geometry library for prototype unions and exact collision boundaries only after numerical validation. Human-readable configurations and deterministic tests matter more than a large UI.

## 7. Who should do what (specific to current project maturity)

- **Fable scientific review FIRST:** One focused adversarial pass over M0–M4, the geometry convention, open-system energy ledger, pressure-triggered growth and novelty scope. Ask for exact counterexamples, impossible assumptions and 3–5 minimum tests. Do not spend weeks on a philosophical panel.
- **Opus 5.5 / Claude Code (fresh independent workspace) for implementation:** Give just this frozen brief, the validated interface requirements, existing HardDisks technical constraints, and minimum acceptance tests. First deliver a tiny M1 prototype plus tests, not the whole spec. Keep code in a **new repo**; ask before touching `HardDisks/hspist3`. Require a test report and exact links/commit IDs for completed work. 
- **NeurAIpil (Sol + Opus) as optional *one-time architecture deliberation test*** after Fable clarifies the spec. Its proven A3 stage uses a frozen content-addressed bundle, heterogeneous proposal/cross-critique/synthesis and **stops at `awaiting_plan_approval`**—do not present it as an autonomous build engine. Give it a frozen M0 spec, compare its output with single reviewer advice and cap provider cost. Building MorphoTransport must **not** depend on upgrading NeurAIpil.
- **IdeaLabs** is the research knowledge graph / record and build-brief hub, **not** the simulator repository or the NeurAIpil architecture engine. Preserve independent project boundaries; attach references, decisions and provenance to the project in IdeaLabs.
- **Available Claude subscription/credit pool until Nov 4:** User says approximately **$250**, but exact terms, expiry and rate limits were not checked. Allocate effort to test-producing coding/reviews with strict bounded prompts; avoid multi-agent loops for every feature or rewriting all existing systems. Reserve meaningful capacity for Paper 1/2 and DIB validation, and give the parallel engine only a capped vertical slice. Assess balances and plan rules in the provider UI before assuming credits will pay for any particular model or automation.

**Decision gate:** If Fable finds basic inconsistencies, resolve before coding. If M1/2 tests fail, fix geometry/EDMD; if the M3 transport experiment adds nothing against graph-based baselines, reduce or postpone the collision-resolved layer. Do not trade Paper 1/2 publication progress for an open-ended simulator platform.

## 8. Two copy-ready handoff briefs

**Fable (critique, not implementation):**

> Critically review `MorphoTransport_Network_Engine_Blueprint_2026-10-15.md` with the IdeaLabs research notebook and coupled-role literature audit. The experiment is a single growable 2D channel network: large disks grow a centerline graph and swept channel walls; separate small EDMD disks transmit material and mechanical signals; slow measured feedback adjusts geometry. Different policy plug-ins approximate Physarum, neurite development, DLA and eventual empowerment/IB/DIB. Examine physical feasibility, solid-versus-lumen geometry, EDMD collision event correctness, material and work ledgers, pressure-triggered growth energy supply, Physarum peristaltic work, strongest minimal benchmark, relevant prior art, and which assertions must be tested before claiming learning. Keep Paper 1/2 untouched. Return (a) blockers, (b) minimal acceptance tests, (c) lowest-risk M1–M3 implementation order. Do not rubber-stamp novelty.

**Opus/Claude Code (implementation, *only* after M0 review):**

> In a NEW `MorphoTransport` repository, implement **M1 ONLY** from the attached blueprint. Do not modify HardDisks/hspist3 or IdeaLabs/NeurAIpil. Deliver a deterministic graph of attached variable-radius structural disks, nonconvex smoothly tapered swept-channel geometry with Y junctions and one loop, and an exportable 2D collision-boundary representation. Add geometry/invariant tests: consistent graph topology, open junction continuity, radius/width bounds, no cracks, seeded replay, geometrical area stable against tessellation, and material-area count without overlap double-counting. Provide a minimal 2D viewer and README containing exact commands. Run tests, record limitations and stop; do not implement transport, biological rules or electrical features until M1 is approved.

## 9. Scientific north star and scope

We are not trying to claim Physarum, neurons, and ion conductors obey one law. We are building a **shared physically explicit substrate** on which multiple growth hypotheses can be compared under consistent geometry, transport and energy constraints. The potentially new scientific result is whether *coupled exploration, information-selective reinforcement, predictive structural retention and energetic selection* produces better/meaningfully different morphology or transmission than strong local rules under matched conditions—and whether explicit collision transport changes the conclusion of a coarse graph model. Every separate effect must be demonstrated, not implied by an attractive animation.

**Status:** Architecture/physics plan only. No Fable/Sol/Opus consultations or software build claimed in this document.


**Revision log v1.9 (2026-10-15):** Added canonical MorphoTransport engineering blueprint, dual-population growing 2D channel architecture, physical transport/pressure/coupling ledger, future membrane ion/cable/memristor/demon layers, delayed pressure growth, and staged collaborator/credit strategy. Preserved all preexisting research content and paper plans.


## 46. [2026-10-15 / v2.0] **Existing HardDisks SDL2/EDMD reuse, two-food growth scene, and scientific naming correction**

**User correction:** Previous v1.9 blueprint captures the **scientific architecture**, but omitted an important implementation shortcut: the HardDisks project reportedly already has an **SDL2 visualization window and an EDMD physics channel**. The first shared screen can show structural growth, coarse-grained geometry and **two sugar/nutrient patches**. It is incorrect to assume a separate renderer or wholly independent EDMD engine must be written. However, the codebase API has not yet been inspected here; code sharing is a *technical decision pending M0 audit*, and arbitrary open/curved collision channels are not assumed to exist. Preserve the validated `hspist3` Paper 1/2 workflow.

**Scientific umbrella (provisional, not a final brand):** **Physical Information Dynamics**. Earlier `MorphoTransport` was an unsatisfactory engineering working title because the research is fundamentally about energy–information–matter–morphology relationships, not just hydraulic transport. Keep the repo/app name undecided pending the implementation audit and user preference.

**First 2D visual experiment:** Start one growing structural-disk seed in a world with **two separate sugar/resource patches**. Grow a connected disk-tip graph, construct a nonconvex smooth 2D channel envelope, and use the *existing SDL2 viewer if integration is feasible* to toggle structural disks, geometry, graph, resource/cue field, and later particle/flux/pressure/signal fields. Initially patches are externally modeled environmental sites (not literal sugar disks); cue fields must have declared rules. Small collision-resolved disks travel inside a frozen grown hull only in the subsequent physically validated transport stage. Then introduce true outlet particle counts and flow; finally couple measured flow to radius modification, Physarum-inspired pulse phases, prospective branch retention, and energy/material costs. No numerical disk count as material/energy proxy.

**Major build correction:** A new versioned **Physical Information Dynamics — Growing Network Engine Blueprint v2.0** supersedes the v1.9 advice to build an independent viewer and repository unconditionally. M0 starts by inspecting actual `HardDisks` SDL2 and EDMD source and identifying components safe to share without altering `hspist3`. New isolated application target sharing code is acceptable; so is a separate repo with clean adapters if actual code coupling makes that safer.

**Scientific roles preserved:** Empowerment-informed exploration; IB-based selective reinforcement with an explicit relevance variable; DIB predictive retention; real physical resource constraints and Still-inspired information-thermodynamic questions that need proper feedback derivations. Long-term transport models include EDMD, physically calibrated edge conductance, and *separately* charge, cable/membrane dynamics, excitable impulses, ionic/memristive states and Maxwell-demon-style controller ledgers. No assertion that dense elastic hard disks directly become electricity.

**Immediate priority stays Paper 1/2 and DIB validation.** Network Demo 0 is a bounded parallel programming task, not an alteration of their scientific protocols. For details, see `Physical_Information_Dynamics_Engine_Blueprint_v2_0_2026-10-15.md` and the updated Fable review brief.

**Revision log v2.0:** Corrected existing SDL2/EDMD reuse, two-resource initial visualization, coarse-grained displays, and provisional information-first naming. Complete historical content maintained.
