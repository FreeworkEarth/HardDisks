---
title: "Physical Information Dynamics — growing-network physics engine blueprint (provisional title)"
date: 2026-10-15
status: "proposal / engineering design; no implementation or new physical result claimed"
project: "Physical Information Dynamics (provisional; repo name undecided)"
workspace: IdeaLabs (research memory only)
related_simulation: HardDisks/hspist3 (keep validated Paper 1/2 protocols untouched; evaluate code reuse)
scientific_supervision: Susanne Still — to discuss, not pre-approved
---

# Physical Information Dynamics — growing physical-network engine (working title)

> **User's precise request, 2026-10-15:** Grow a network from connected **large structural disks**, derive a **continuous tapered 2D channel boundary** (a swept-disk union / lumen, NOT a convex hull that fills gaps between branches), and run **smaller collision-resolved hard disks through the interior**. Observe particle/energy/pressure/signal transport; allow slower structural growth, retraction, diameter adaptation and topology changes in feedback. Use a **single engine with pluggable growth policies**, not a different simulator for every organism. In the future investigate physically modeled conductance, membrane ion channels, action-potential propagation, history-dependent/memristive connections, and Maxwell-demon-type sensing/feedback costs. **These are later research options.** Preserve the current HardDisks experiments as-is.

> **Naming decision (2026-10-15):** `MorphoTransport` was only a temporary engineering label and is rejected by the researcher as too transport-specific. **Physical Information Dynamics** is a *provisional scientific umbrella*, not an agreed permanent project or repository name. Alternative future labels include **Information–Matter Dynamics**, **Adaptive Physical Information**, and **PhysInfo**. The engine is described generically as the **growing physical-network engine** until a name is chosen. The objective is to test candidate first-principles mechanisms, not to announce that those principles are established.

## 0. Hard gate: do not delay or alter the existing research

**Paper 1:** Finish Román two-compartment static/equilibrium stiffness, adiabatic/acoustic response, sound-speed validation, and unresolved finite-size/FFT checks. 

**Paper 2 (IMMEDIATE):** Existing `hspist3`: right-hand **prescribed constant-speed piston** `x_p(t)=x_0+vt` over a declared window; multiple impermeable movable dividers; a spring-backed movable terminal wall on the **left**. Measure input force/work `W_in = ∫ F_p \dot x_p dt`, output spring displacement `x_s(t)` and energy `ΔU_s = 1/2 k[x_s(t)^2 - x_s(0)^2]`, and the full work–energy ledger, with `ΔU_s/W_in` reported only for clear positive input work and fixed endpoints. Mobile gas disks stay within their compartments while transmitting momentum to walls. **No network growth code, ion dynamics or rewritten event scheduler enters this validation branch.** DIB implementation and its no-feedback dissipation derivation remain an independent third effort.

## 0A. CORRECTION v2.0 — reuse the EXISTING HardDisks SDL2 and collision code before writing new tools

**Source of this correction:** The researcher explicitly clarified that their HardDisks codebase already has an SDL2 visualization window and an EDMD hard-disk channel/physics implementation. Their desired first experiment has a visible growth seed and **two sugar/nutrient resource patches**. Prior v1.9 implementation guidance suggested a separate viewer and an independent implementation before any code reconnaissance. That is *too prescriptive*. The correct engineering principle is **reuse-first with isolation of validated Paper 1/2 code**. This is user-reported codebase capability; **no source files or API have yet been inspected**. Do not assert that current `hspist3` supports curved, open, growing channels or that the SDL2 interface is reusable without refactoring.

**Mandatory M0 repo audit before coding:** In the real local `HardDisks` repository, locate the SDL2 window/event loop/drawing functions, particle renderers, coarse-graining or binning overlay (if implemented rather than only visually feasible), existing EDMD event queue and hard-disk/hard-wall collision code, renderer-physics coupling, experiment configuration and regression tests. Produce a factual interface map and note which pieces can be linked without modifications. Check ownership of simulation clocks and rendering frames. Record what exists versus what needs an adapter. Read-only reconnaissance must not modify the validated Paper 1/2 engine.

**Choose source layout ONLY after that audit:** Preferred if the architecture supports it: a **new `growing_network` / `physical_information` executable or module in the HardDisks workspace**, sharing an extracted/linked SDL2 renderer and a well-defined EDMD kernel interface, while leaving `hspist3` unchanged. If existing renderer/event loop and solver are too tightly coupled, build a separate sibling repository that imports isolated components through an adapter or uses a deliberately scoped extraction; do not clone a giant hard-disk solver and claim it was reused. A new git repository is **not** inherently a scientific requirement. Any extraction must preserve exact baselines and tests for Paper 1/2; no broad refactor of `hspist3` until the active physics experiments are secure.

### First SDL2 screen and physics-independent visualization

Initial seeded scene, **growth without transported particles**:

```
              FOOD / SUGAR A                   FOOD / SUGAR B
                 [patch A]                       [patch B]
                   ***                             ***
                  /                                /
    SEED O---o---o---o--- tip                     /
                   \---o---o---tip--- ... -------/

   overlay toggles: [structural disks] [centerline graph] [lumen hull]
                    [nutrient/chemoattractant field] [local material]
                    [later: small EDMD disks, speed, packing, flux, pressure]
```

`A` and `B` are **two environmental food/nutrient patches or controlled resource sources**, not automatically a pair of moving sugar molecules. Begin with two fixed patches and a declared spatial cue field (for example diffusion/decay or specified distance-biased chemoattraction); make their geometry and stimulus schedule deterministic. Later vary which resource is active or valuable to test predictive versus unpredictable demand with identical marginal frequencies. A true chemical-diffusion model must define its PDE/boundary conditions and numerical solver separately.

**SDL2 coarse-graining/inspection modes** (all projections of the *same* underlying state, not independent fake simulations):
- micro: large structural disks + tips, optional small EDMD particles once enabled;
- meso: continuous swept-channel hull and contact/junction status, with center-accessible free space for finite-size transport particles;
- macro: sparse graph, edge widths/radii, age, flow capacity, material distribution and predicted/policy state;
- measured field overlays: binned local counts/packing, `v_x`, speed, pressure estimator, flux, arrival times, controlled-signal response and eventually information observables. Clearly label absent/unmeasured fields rather than invent colorful heatmaps.

**Initial demo acceptance tests:** One seed, two resource patches, at least two seeded growth rules (`unbiased`, `resource-guided`), pause/step/replay, rendering of graph+disk-chain+swept hull, one dynamic radius change, numerical geometry/material-budget invariant, and an SDL2 interactive overlay toggle. Preserve reproducible output in a headless run as well; the GUI must not be the sole source of results. Distinguish visually plausible `Physarum` contraction from a *physically validated* peristaltic mechanism.

### EDMD reuse needs an explicit interface, not a promise

`TransportBackend` should accept a static or quasi-static **collision-boundary representation** and particle states; emit collision events, pressure/flux data and conservation diagnostics. Existing hard-disk particle–particle event prediction can often be reused. **But** the present compartments/piston walls need not support arbitrary circular arcs, junctions, open reservoirs or moving topology. Audit existing event predictions, invalidation, grazing contacts, intersections and acceleration assumptions. Do static corridor and polygonal Y validation before curved hull; no changed geometry while collisions are scheduled unless every affected event is invalidated. Once geometry changes, enforce collision-free remapping and account for boundary work and displaced particle energy; do not silently teleport or rethermalize particles.

### Concrete experimental progression

**Demo 0:** Seed→two sugar patches, structural disk-tip growth, source-cue visualization and graph/hull overlays in reused SDL2, no EDMD yet. **Demo 1:** Same frozen grown hull receives *smaller mobile hard disks* and shows real collisions/flux through a validated static corridor/Y. **Demo 2:** Add inlet reservoir and two outlets and measure outlet event counts/first-passage times, open-boundary energy and matter bookkeeping. **Demo 3:** Slowly adjust channel radii after an EDMD measurement window and compare against frozen/no-growth baseline; mark quasi-static approximations and missing work. **Demo 4:** Introduce Physarum-style pulses, local flow feedback, IB relevance, DIB predictive retention, empowerment-guided exploration and explicit material/energy costs in tested increments and ablations. Only later: ion physics and excitable cable models.

This is a **growth-and-transport multiscale simulation**: render and measure micro (disks), meso (channels), and macro (graph/information) without confusing a coarse-grained *visualization* with a derived equation or validated effective model. Existing SDL2 should accelerate inspection, not predetermine new physical claims.


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

**M0: Scientific design freeze and EXISTING-REPO RECONNAISSANCE (few days, bounded).** Inspect HardDisks SDL2/EDMD and choose an isolated new executable/module or sibling repository *after the audit*, without touching validated hspist3 behavior; archive a minimal numerical protocol, data units, geometry convention, correctness gates and baseline sources. Require reviewer to challenge before significant coding. Code should be architected as a reusable simulation engine with policies as plug-ins, not one graph per organism.

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
growing_physical_network/  # provisional folder name; choose after HardDisks repo audit
  README.md                     # purpose, scientific limits, Quickstart
  docs/SPEC_M0.md               # sign-off review & invariants
  docs/PHYSICS_LEDGER.md        # dimensions, open/closed boundaries, work
  src/geometry/                # disks, centerline graph, swept channels, junctions
  src/growth/                  # DLA, Physarum, neurite, neutral, E/IB/DIB adapters
  src/transport/               # EDMD adapter/kernel (frozen geometry first)
  src/coupling/                # observations, budget, delayed update executor
  src/measurements/            # flux, pressure, transit time, info and energy
  src/protocols/               # Y tests, piston, demand switching, seeded configs
  src/viewer/                  # prefer reuse/adaptation of existing HardDisks SDL2, after API audit
  tests/geometry/ tests/edmd/ tests/budget/ tests/replay/
  experiments/README.md       # immutable run manifests; no massive files in Git
```

Practical language choice is **not yet prescribed**: Python for initial geometry, scientific orchestration and plotting is convenient, with compiled EDMD kernel when performance or deterministic event ordering requires it; an SDL2 renderer is reportedly already present in HardDisks; reuse through a thin interface when proven feasible. Require a clean interface, not an unnecessary wholesale rewrite of your existing C/SDL2 hard-disk simulator. Use a known geometry library for prototype unions and exact collision boundaries only after numerical validation. Human-readable configurations and deterministic tests matter more than a large UI.

## 7. Who should do what (specific to current project maturity)

- **Fable scientific review FIRST:** One focused adversarial pass over M0–M4, the geometry convention, open-system energy ledger, pressure-triggered growth and novelty scope. Ask for exact counterexamples, impossible assumptions and 3–5 minimum tests. Do not spend weeks on a philosophical panel.
- **Opus 5.5 / Claude Code (isolated working branch or sibling workspace) for implementation:** Give just this frozen brief, the validated interface requirements, existing HardDisks technical constraints, and minimum acceptance tests. First deliver a tiny M1 prototype plus tests, not the whole spec. Decide **new isolated application in the HardDisks workspace versus a sibling repo** only after the read-only SDL2/EDMD audit; ask before touching `HardDisks/hspist3`. Require a test report and exact links/commit IDs for completed work. 
- **NeurAIpil (Sol + Opus) as optional *one-time architecture deliberation test*** after Fable clarifies the spec. Its proven A3 stage uses a frozen content-addressed bundle, heterogeneous proposal/cross-critique/synthesis and **stops at `awaiting_plan_approval`**—do not present it as an autonomous build engine. Give it a frozen M0 spec, compare its output with single reviewer advice and cap provider cost. Building the growing network engine must **not** depend on upgrading NeurAIpil.
- **IdeaLabs** is the research knowledge graph / record and build-brief hub, **not** the simulator repository or the NeurAIpil architecture engine. Preserve independent project boundaries; attach references, decisions and provenance to the project in IdeaLabs.
- **Available Claude subscription/credit pool until Nov 4:** User says approximately **$250**, but exact terms, expiry and rate limits were not checked. Allocate effort to test-producing coding/reviews with strict bounded prompts; avoid multi-agent loops for every feature or rewriting all existing systems. Reserve meaningful capacity for Paper 1/2 and DIB validation, and give the parallel engine only a capped vertical slice. Assess balances and plan rules in the provider UI before assuming credits will pay for any particular model or automation.

**Decision gate:** If Fable finds basic inconsistencies, resolve before coding. If M1/2 tests fail, fix geometry/EDMD; if the M3 transport experiment adds nothing against graph-based baselines, reduce or postpone the collision-resolved layer. Do not trade Paper 1/2 publication progress for an open-ended simulator platform.

## 8. Two copy-ready handoff briefs

**Fable (critique, not implementation):**

> Critically review `Physical_Information_Dynamics_Engine_Blueprint_v2_0_2026-10-15.md` with the IdeaLabs research notebook and coupled-role literature audit. The experiment is a single growable 2D channel network: large disks grow a centerline graph and swept channel walls; separate small EDMD disks transmit material and mechanical signals; slow measured feedback adjusts geometry. Different policy plug-ins approximate Physarum, neurite development, DLA and eventual empowerment/IB/DIB. Examine physical feasibility, solid-versus-lumen geometry, EDMD collision event correctness, material and work ledgers, pressure-triggered growth energy supply, Physarum peristaltic work, strongest minimal benchmark, relevant prior art, and which assertions must be tested before claiming learning. Keep Paper 1/2 untouched. Return (a) blockers, (b) minimal acceptance tests, (c) lowest-risk M1–M3 implementation order. Do not rubber-stamp novelty.

**Opus/Claude Code (implementation, *only* after M0 review):**

> [SUPERSEDED by the new v2.0 implementation handoff below] Previously: in a NEW `MorphoTransport` repository, implement **M1 ONLY** from the attached blueprint. Do not modify HardDisks/hspist3 or IdeaLabs/NeurAIpil. Deliver a deterministic graph of attached variable-radius structural disks, nonconvex smoothly tapered swept-channel geometry with Y junctions and one loop, and an exportable 2D collision-boundary representation. Add geometry/invariant tests: consistent graph topology, open junction continuity, radius/width bounds, no cracks, seeded replay, geometrical area stable against tessellation, and material-area count without overlap double-counting. Provide a minimal 2D viewer and README containing exact commands. Run tests, record limitations and stop; do not implement transport, biological rules or electrical features until M1 is approved.

## 9. Scientific north star and scope

We are not trying to claim Physarum, neurons, and ion conductors obey one law. We are building a **shared physically explicit substrate** on which multiple growth hypotheses can be compared under consistent geometry, transport and energy constraints. The potentially new scientific result is whether *coupled exploration, information-selective reinforcement, predictive structural retention and energetic selection* produces better/meaningfully different morphology or transmission than strong local rules under matched conditions—and whether explicit collision transport changes the conclusion of a coarse graph model. Every separate effect must be demonstrated, not implied by an attractive animation.

**Status:** Architecture/physics plan only. No Fable/Sol/Opus consultations or software build claimed in this document.

## Copy-ready revised implementation handoff (v2.0, SUPERCEDES the v1.9 M1 prompt)

> Inspect the existing `HardDisks` source **read-only first**. Identify reusable SDL2 render/window/event loop, coarse-graining/visualization capabilities and exact EDMD physics/collision event APIs, including what does NOT handle arbitrary curved/open channels. Write a short reuse-versus-separate-target technical decision and safeguard the currently validated `hspist3` Paper 1 and 2 behavior. Then implement only **Growth Demo 0** in an isolated executable/branch/sibling module as supported by inspection: one root of connected large variable-radius structural disks, a topology graph, a smooth nonconvex swept-channel hull, **two sugar resource patches**, an unbiased growth policy and a resource-guided policy, fixed-seed replay, material conservation/area invariance, and SDL2 toggles to view disks, hull, graph and an explicitly defined resource field. Reuse the existing renderer if technically viable; do not rebuild a second UI needlessly. No mobile EDMD particles, ion physics, general autonomous agents, or modifications to production piston experiments yet. Stop with tests, screenshot/demo, numerical replay and a list of actual shared components. Once approved, separately integrate frozen-geometry small hard-disk EDMD through the channel, then outlet flux measurement. Treat 'Physical Information Dynamics' as a provisional scientific title, not an imposed repository name.


**Revision log v2.0 (2026-10-15):** User rejected the MorphoTransport project name as too narrow; added provisional *Physical Information Dynamics* scientific title, explicit reuse-first SDL2/EDMD audit, two sugar patches in the first growth demo, coarse-graining/inspection overlays, code-sharing decision gate and corrected copy-ready Opus brief. Older advice for a wholly independent viewer/repo is superseded where it conflicts.
