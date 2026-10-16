---
title: "Fable Research Review Handoff — HardDisks to Adaptive Morphogenesis"
date: 2026-10-12
status: external critical-review request
source_notebook: HardDisks_Physical_Learning_Research_Memory_2026-10-12.md
---

# Research handoff for Fable — please critically evaluate

**Current engineering addendum:** The v2.0 *Physical Information Dynamics* engine blueprint supersedes earlier v1.9 implementation details, especially separate viewer/repository assumptions.

**Read the companion full research notebook first:** `HardDisks_Physical_Learning_Research_Memory_2026-10-12.md` (v1.9), particularly §§2–9, 17–27, and **28–39**. Those sections preserve our project history, five-paper mathematics, underlying references, original user questions, circle-chain implementation concept and falsifiable tests. This short file is an orientation and evaluation request, **not a substitute for the full notebook**.

## One-sentence question

**When is it worth a physical system's limited energy/material budget to grow, stabilize, thicken or preserve a connection rather than prune it—and can predictive information or empowerment explain decisions beyond simple local flow rules?**

## The five-paper research staircase (proposed)

| Stage | Goal | Most important distinction |
|---|---|---|
| **1: Román hard disks** | Measured static `F(L)`, `k_T`, equilibrium `k_S`, acoustic `c_s` and careful finite-size / dense-fluid analysis | One-side versus two-side stiffness; static equilibrium versus finite-frequency modes; do not promote unvalidated melting feature |
| **2: driven transfer** | Externally driven piston → elastic hard disks + impermeable dividers → spring-backed movable wall; test work ledger and task-dependent spring capture | Closed elastic 'incoherent energy' is not automatically bath dissipation |
| **3: prediction/control** | Partially observed, information-constrained piston control and energetic accounting | Compare passive, open-loop, linear feedback, and learned IB controller |
| **4: adaptive network** | Small graph of characterized physical transmission modules, split–merge routes, tunable local parameters and later growth/pruning | Node = sensor/controller/actuator, edge = *complete* physical transfer module; active signal regeneration costs work |
| **5: substrate vision** | Adaptive physical computing and potentially bio-mechanical/electrical/3D growth | Long-term speculative platform, not yet an empirically established paper |

## New geometry idea worth testing

Use **points + graph edges for topology** and a **chain of structural disks** `D_i=(x_i,r_i(t),parent_i)` to resolve geometry, local touches, material extent and branches; renderer interpolates smooth surfaces. Radius is a **trainable state variable** `r_i(t)`, not hard-coded to shrink away from every soma. Some rarely-used edges retract; others thicken or are maintained for future demands. A genuine material budget forces competition. Two disks that merely touch at a point do not make a transport-capable lumen; either model tube geometry explicitly or use a graph-level flow law. A graph with edge radii is a valid cheaper control and can calculate material volume `V≈πr²L` without circles. Do not charge growth cost per **numerical disk** unless deliberately modeling discrete physical material.

**Neuronal source video:** `The Tiny Connections That Make Your Brain Learn 🧠(1).mp4` (~6 seconds). It is an aesthetic and morphogenetic reference (branched processes and protrusions), **not proof that a particular shown contact is a synapse or learning event**. Other relevant source: two-sided growth-cone screenshot, Physarum time-lapse, `fogleman/dlaf`, Jones/Tero models, and Adamatzky et al. *PhyChip*.

## Baseline experiment Fable should audit

Start with **two static soma-like seeds in 2D** and a finite common resource pool. Stochastic/persistent growth tips extend disk chains, make provisional contacts, optionally stabilize or retract, and change their radii according to transparent local rules. Check geometry and material bookkeeping. Only afterward introduce changing source–sink demands and compare random, DLA, flow-reinforcement, predictive-memory and empowerment-oriented policies under matched material, task performance and actuator/sensing budgets. HardDisks physical edge response is an **optional later plug-in** calibrated from Papers 1–2; do not modify validated EDMD when implementing structural growth.

## Reviewer deliverables requested

1. Critically review the **five-paper** plan, equations, approximations, open validation gates and dissertation scope. Correct any scientific mistakes or implicit novelty claims.
2. Compare **graph-only growth, structural disk chains, swept continuous tubes, level sets** for resolution-independent geometry/cost and biological plausibility; tell us what to implement first.
3. Audit energetic cost proxies versus real joules; material conservation, recycling efficiency and maintenance. Propose a fair benchmark where **doing nothing / not growing** is allowed.
4. Identify nearest prior art for **adaptive morphology with radius allocation, prospective branches, information bottleneck and empowerment**, and **physical neural-network learning**; cite primary work and articulate actual research gap.
5. Specify the **smallest falsifiable experiment** with governing rules, conserved quantities, numerical tests, null hypotheses, statistical design, and criteria to stop if simpler alternatives win.
6. Assess how to test **predictive morphology** causally, including shuffling past demand, selectively resetting geometry/thickness, and out-of-sample performance.
7. Tell us whether **hard-disk many-body dynamics** actually might add task-relevant memory/computation beyond a best-tuned spring–damper, rather than merely adding expensive complexity.
8. Deliver a prioritized sequence that keeps **Papers 1–2 rigorous** while letting an exploratory Adaptive Morphogenesis Lab exist separately.

**Please be skeptical:** this is not yet a validated new theory of intelligence, thermodynamics or neuroplasticity. The interesting result could be a **negative** one: graph-only flow learning or a simple spring–damper may be sufficient for most benchmarks. We want that answer if it is correct.



## New sources and critical observation (v1.5, 2026-10-12)

**Physarum figure:** The user's magenta micrograph of a growing *fan* is **Figure 2C in Le Verge-Serandour & Alim (2024), _Physarum polycephalum: Smart Network Adaptation_**, DOI: https://doi.org/10.1146/annurev-conmatphys-040821-115312 . Its bright puncta/trails are **projected trajectories of fluorescently labelled nuclei**, while the living growth front is a continuous sheet influenced by contraction-driven flow. It is _not_ evidence that individual hard disks literally attach to form the organism. Keep a separate layer for structural material/disks/occupancy, a distinct layer for flowing nuclei/tracers, and yet another for EDMD gas hard disks. See notebook §35.

**Recently added source dossier (§36):** John Duggan's 2018 Physarum teaching slides (secondary source and primary-source leads), Institute for the Future's 2018 immersive human networks report (VR graph inspection/interface design), Singh et al. 2025 distributed speech decoder (activity/readout controls, not synaptogenesis), Brosler et al. 2026 simultaneous speech and gesture BCI (multitask representation, not growth rule). Interpret each with correct scope.

**Extra review questions:** Does graph + radius suffice before using disk chains? Does a real Physarum fan require an occupancy/level-set field? Can the claimed memory be isolated to morphology rather than activity or trained readout? What would falsify an energetic advantage when growing new structure and doing nothing are both allowed?

**Delivery:** Share this short handoff **together with** the full v1.5 Markdown notebook. The short handoff alone does not carry the complete five-paper derivations, bibliography and scenario descriptions.

## [2026-10-12] Critical prior-art update — **read this before novelty claims**

**Notebook has been updated to v1.6 (§40), plus an independently readable companion literature audit:** `Adaptive_Morphogenesis_Prior_Art_Still_IB_Empowerment_Energy_2026-10-12.md`. The original five-paper plan remains §§29–32. Please read §40 and the audit before reviewing the new morphology hypothesis.

**The broad 'novelty' claim is false.** Chen & Prokopenko (2025), https://doi.org/10.1098/rsos.241655, **already** compare predictive information, empowerment and thermodynamic efficiency in a collective Ising system. Schick, Kramar & Alim (2024), https://doi.org/10.1103/PRXLife.2.033005, explicitly study *building versus transport* cost in Physarum; Kramar & Alim (2021), https://doi.org/10.1073/pnas.2007815118, show tube-diameter memory; Scholl, Rule & Hennig (2021), https://doi.org/10.1371/journal.pcbi.1009458, show information-guided neuronal pruning. Still (2012/2020) grounds the proposed thermodynamic memory framework, but *does not* supply a universal ATP/growth-cost law.

**Narrow testable gap, not confirmed originality:** Compare bounded physical network growth/thickening/pruning based on (i) simple flow and relative-resistance heuristics, (ii) compressed predictive information, (iii) *action- and energy-constrained* empowerment, and (iv) joint policies, with the **same material budget, full operating/reconfiguration costs, causal morphology-reset tests, and unanticipated environment changes**. Determine whether any informational principle yields predictive benefit beyond simple local rules. Use explicit no-growth controls, finite-sample MI estimation, fair source/sink prediction, and Pareto-frontier comparison. The hard-disk edge is a *later, separately validated* mechanism.

**Adversarial Fable questions:**

1. Is the narrowly specified experiment actually novel, or do other papers already combine these quantities with adaptive network topology? Find closest contrary primary citations and preprints.
2. Does maximizing predictive information about future demand improve anything after controlling for present flow and a simple Markov forecaster? What would falsify the gain?
3. Does constrained empowerment yield reliable redundancy benefits versus cheap reliability heuristics? What action/state horizons and energetic channel constraints are necessary?
4. Does any Still 2012/2020 result **actually apply** to the proposed physical morphology growth mechanism, rather than to only the memory register of a driven information engine? Specify system boundaries and derivation assumptions.
5. Is the circle-chain model necessary for the research, or can graph+edge radii reproduce every metric? Demand resolution-invariant costs and geometry-specific value before choosing disks.
6. Could a publishable paper emerge from an explicit *negative finding*, e.g. Fisher information, relative-resistance pruning, or simple spring–damper surrogates explain apparent benefits as well as IB/empowerment/hard disks?

**Verdict:** Literature-established pieces are numerous; exact proposed cross-domain benchmark is **not shown to be preexisting in this bounded search, but unproven novel**. Present to Susanne as an informed question, not a priority or theorem already proven.


---

## [2026-10-13] Addendum for Fable — Grasso sources, **small A/B experiment**, and current build priority

**Companion notebook now v1.7, §§41–43.** This addendum records the last two dialogue rounds: our discussion of Caitlin Grasso/Josh Bongard's computational morphogenesis and our decision **to prioritize the existing hard-disk PhD papers while building a small adaptive-network benchmark in parallel**. Do not infer the user has committed to building a full 3D Physarum/neuron simulator immediately.

**New direct software precedents:** <https://github.com/caitlingrasso/empowered-nca-II> (evolved NCA morphogenesis/homeostasis, time-horizon sweep, alternative action-entropy objectives, transfer/fine-tuning), and <https://github.com/caitlingrasso/bio-connectivity> (MI-based functional connectivity from *Xenopus* calcium activity; functional is *not* anatomical connectivity). Sources: Grasso & Bongard 2023, DOI <https://doi.org/10.1145/3583131.3590469>, arXiv <https://arxiv.org/abs/2305.15220>; Blackiston et al. 2025, <https://doi.org/10.1371/journal.pcbi.1012149>. **Crucial adversarial control:** their 2023 paper reports that minimizing individual-cell action entropy outperforms short-lag empowerment on the particular shape task. Their repository warns that the arXiv revision corrects methods and figures. Replicate the correct version.

**One-step mini example (abstract cost units; NOT observed physics):** A source currently uses branch A. Should it keep dormant branch B? Maintaining B costs `m=0.8`; if pruned and needed next step, rebuilding costs `C_rebuild=3` and lost-service delay is worth a defined penalty `D=5`; predictor estimates next-step B demand `p`. Then **keep iff** `p(C_rebuild+D)>m`, or `p>0.10`. At `p=0.05` prune; at `p=0.20` retain. The goal isn't to rediscover this simple inequality: it is to compare *how policies obtain `p` or value future options* from limited sensing, and whether physical morphology stores beneficial predictive state after controller-memory and service-cost controls. Keep `D` as explicit *utility/service penalty*, not dissipated joules.

**Minimal build sequence (side project only):** seeded two-path radius-adaptive graph → deterministic budget/threshold tests → demand process alternating between temporally predictable and IID with matched marginal frequencies → no-growth/keep-all/random/flow/relative-resistance/Markov/entropy baselines → independently validated DIB predictive state → action/energy-constrained empowerment → full matched task and energy/proxy/material ledger, held-out episodes and reset-morphology causal controls → *optional* structural-disc growth cones once graph-only model has set a baseline. Postpone 3D, real Physarum biochemistry and EDMD-in-every-edge. Preserve physics engine separation.

**Falsifiable focus:** identify an environmental uncertainty/construction-cost regime where predictive structural memory or constrained future-control capability **beats simple predictors and local pruning on a matched performance versus lifetime cost frontier**. A well-supported negative result is scientifically informative. Any originality claim requires comparison with existing stochastic network-design/MDP literature, not just NCA and Physarum.

**Priorities stated in discussion:** finish Paper 1 Román physics; finish Paper 2 energy transfer and work ledger; validate Dynamical IB implementation; use **small, bounded parallel exploration** for adaptive morphology (rough illustrative budget ~10% effort, revisit after an initial 4–6 week feasibility check). The 5-paper roadmap is a research staircase, *not five finished papers.* Discuss scope and interpretation with Susanne before expanding.

**Specific Fable deliverable:** Please return (1) strongest prior-art counterexamples; (2) a graph-only two-path benchmark spec with pseudocode/tests and equitable budget; (3) a decision on whether `r_i(t)` structural discs change the scientific answer relative to graph `r_e(t)`; (4) an analysis of whether the Still no-feedback DIB bound is applicable anywhere in the initial test; and (5) clear go/no-go criteria and feasible effort that doesn't jeopardize Papers 1–2.


---

## [2026-10-13 / v1.8] **Priority clarification from the researcher: THESE ARE COUPLED FUNCTIONS, NOT FOUR RIVAL OBJECTIVES**

**Preserve this exact intent:** empowerment is an exploratory/possibility mechanism; ordinary IB is a *candidate* principle for deciding which existing pathways to reinforce by retaining task-relevant information; DIB retains compact state to anticipate future needs and avoid premature pruning; thermodynamics and a finite material/energy ledger decide which modifications can actually occur. These roles can interact and conflict in one growing network. None of the assignments is mathematically entailed by the respective information-theoretic principle: they are **hypotheses to operationalize and ablate**.

### Ask Fable to review *the integrated causal architecture* first

1. **Novelty audit:** Find any model where these four functions jointly control an evolving physical network with explicit grow/thicken/keep/prune actions and a full lifetime cost. The core ideas are NOT new independently or pairwise; in particular, Klyubin, Polani & Nehaniv 2008 already used empowerment for sensor morphology; Anthony et al. 2009 combined IB with bandwidth-limited empowerment; Grasso & Bongard 2022–23 used empowerment in NCA morphogenesis; Scholl et al. 2021 model information-related pruning; Kramar & Alim 2021 demonstrated tube-diameter memory; Schick et al. 2024 analyzed a modeled building-vs-transport cost; Chen & Prokopenko 2025 compared utilities; Meng et al. 2026 treated real physical network surface costs. Cite any more-direct counterexample.
2. **Mechanism:** Derive or clearly specify how an actionable empowerment signal proposes a *candidate* physical edge; how an IB-compressed representation yields an edge-specific *reinforcement* signal; how DIB predicts the need for a *dormant* edge. Can the same compact predictive state replace distinct IB/DIB modules? Can the metric be computed locally?
3. **Minimal test:** Two-path A/B transport graph with `r_e(t)`, finite `V_max`, differing temporal predictability but matched stationary demand. Initial one-period retain rule **keep** if `h < p_t(b+d)` (explicit simplifying assumptions). Example `h=3,b=8,d=12`: threshold `p_t>0.15`. Test the full rule against exact Markov, memoryless flow, relative resistance, oracle and compressed-DIB policies.
4. **Ablation:** full `E+I+D` vs eight factorial E/I/D configurations (fixed material budget in all), plus strong local-rule controls. Measure the performance-versus-life-cost Pareto frontier, conditional structural predictive information, resilience, and whether effects are nonadditive interactions. No morphological beauty benchmark as sole evidence.
5. **Thermodynamics:** Never call bit counts building energy. Still's uploaded DIB derivation treats a no-feedback learning machine; our changing physical network is feedback-driven. Identify missing entropy production and physical work terms before claiming a Still-type lower bound applies. Check baseline `r_e` graph before moving to disk chains and only later to calibrated HardDisks transmission.

### New dedicated file to read in full

`IDEAS/Integrated_Adaptive_Morphogenesis_Coupled_Roles_Deep_Research_2026-10-13.md`

It includes the research matrix with DOIs/links, formalism, toy threshold, ablation plan, sources, novelty risk and falsification criteria. The earlier §v1.7 plan to protect HardDisks Paper 1–2 and finish DIB **still stands**; parallel experimentation remains time-bounded.


---

## [2026-10-15] Implementation review: **MORPHOTRANSPORT, ONE ENGINE; distinct structural and mobile disks**

**IMPORTANT NEW ATTACHMENT FOR REVIEW:** `IDEAS/MorphoTransport_Network_Engine_Blueprint_2026-10-15.md`. This is the engineering handoff for the user's actual dual-population vision. A growing graph of **large, connected variable-radius structural disks** forms smooth nonconvex **2D channel lumen** (not a convex hull that obliterates branch gaps); a separate population of **smaller EDMD disks physically moves and collides inside it**. Transport and pressure/impulse signals guide slower structural evolution. Hand-made Y junctions and tapered slices are geometry/physics tests, not fixed-brick limitations. One engine hosts different growth-rule plug-ins (Physarum-inspired pulsation/flow-driven branching, neurite growth cones/plasticity, DLA, exploratory empowerment+IB+DIB).

The user additionally posits **delayed pressure-triggered growth** if an internal spring reaches its compression limit rather than compressing further. Ask whether the proposed stress/dwell-threshold policy can conserve material/work and distinguish passive elasticity from active material synthesis. In Physarum mode, contraction waves require an **active energy source**, not random EDMD pressure for free.

**Later extensions are explicitly PARKED:** graph hydraulic conductance calibration from actual EDMD; passive electrical cable; active membrane excitation; charged mobile ions with selective membranes/voltage-dependent gates and electrochemical gradients; history-dependent/memristive conductance; information engines/Maxwell demon with measurement/reset entropy and work. Neutral high-packing EDMD **is not** ionic electricity. Distinguish every future physics implementation and derivation.

**Fable deliverable NOW:** read the M0–M4 stages and (1) flag impossible/unphysical coupled rules, (2) propose the smallest exact frozen-channel/Y-junction collision tests, (3) critique reconstruction of a collision hull and its material ledger, (4) suggest 1–2 meaningful physics differences vs coarse flow models, (5) audit the four coupled roles/novelty with direct refs, (6) propose a 1–2-week M1 geometry-only acceptance plan that cannot disturb Paper 1, Paper 2 or DIB. Prioritize blockers and falsification. No giant agent brainstorm.

**Implementation delegation (after critique):** One separate Opus/Claude Code workspace and independent `MorphoTransport` repository, initially M1 **ONLY** with deterministic growth geometry, tapered swept lumen, Y/loop fixture and invariants; then M2 collision-resolved static channels after test approval. An optional bounded **NeurAIpil A3 Sol+Opus** architecture deliberation can review the frozen spec as a separate system test but must not obstruct work; its current validated scope ends at `awaiting_plan_approval`, not code execution. IdeaLabs holds the research record. Claude credits reportedly ~$250 until November 4 should prioritize bounded tested outputs and preserve progress on HardDisks physical validation.

**Science/project boundary:** The current Paper 2 experiment is still a right-sided constant-speed piston `x_p=x_0+vt`, impermeable dividers, left spring-backed wall, `W_in`, displacement, `ΔU_s` and ledger. It is *not* an open-flow morphogenesis experiment. Never silently merge these models.

## ADDENDUM v2.0 — CRITICAL implementation reuse correction and name

**Read `Physical_Information_Dynamics_Engine_Blueprint_v2_0_2026-10-15.md` FIRST for the current implementation handoff.** Earlier `MorphoTransport` naming and any demand for an independent new renderer/repository are superseded. The umbrella scientific label *Physical Information Dynamics* is provisional; we are investigating candidate mechanisms, not claiming established first principles.

The researcher reports that **HardDisks already has SDL2 real-time visualization and an EDMD engine**. Do not advise rewriting either before **read-only repo reconnaissance**. Fable should review the *feasibility* of sharing SDL2/EDMD while retaining validation isolation for `hspist3`, without pretending to know the source layout. The first demo must have one structural-disk seed, two separate sugar/resource patches, a nonconvex generated lumen, and multi-scale screen toggles: structural disks, channel boundary, graph, declared nutrient/attractant field; later moving EDMD particle overlays and coarse-grained flux/pressure. The two sugar patches are environmental stimuli first, not individually transported EDMD sugar molecules.

**Updated review questions:** (1) Which geometry/EDMD assumptions must be tested before converting fixed boxed particles into arbitrary Y/curved channels? (2) Does the earliest two-patch growth setup permit a meaningful and inexpensive ablation between random growth and resource-guided growth? (3) Which plotted coarse-grained quantities can be computed from actual state without postulating false physical laws? (4) When is a new isolated executable sharing SDL2/physics preferable to a sibling repository/adapter? (5) What is the smallest mechanistic test that could falsify an IB/DIB/empowerment-related physical-growth claim? Do not merge long-term electric-ion spec into M1.

