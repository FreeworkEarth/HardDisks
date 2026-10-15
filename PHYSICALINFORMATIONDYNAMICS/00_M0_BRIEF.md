# M0 brief — Physical Information Dynamics (PID): read-only reconnaissance of HardDisks

For: the Claude Code session working in `PHYSICALINFORMATIONDYNAMICS/`
From: Chris and the Cowork planning chat, 2026-10-15
Scope: M0 only. One deliverable. Stop when it is committed.

Repository root: `/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks`. All paths below are relative to it.

## 1. What this folder is

`PHYSICALINFORMATIONDYNAMICS/` is the sealed home of a parallel project inside the HardDisks repository: a growing 2D physical network — large structural disks grow a centerline graph, the graph is swept into a continuous channel lumen, later small hard disks are transported inside it by collision dynamics, and later still the growth is driven by pluggable policies (empowerment / IB / dynamical IB under a material and energy budget).

The reference spec is `Physical_Information_Dynamics_Engine_Blueprint_v2_0_2026-10-15.md` next to this file. Read its §0, §0A, §5 (M0–M2) and §6. The rest is long-term context; implement none of it.

The project name is provisional. Spend no effort on naming.

## 2. Rules

These protect two papers in progress. They are not negotiable.

1. **Write only inside `PHYSICALINFORMATIONDYNAMICS/`.** Everything else in the repository is read-only for you, above all `hspist3/`. The Stage 1 and Part A pipelines are registered against SHA-256 manifests of `hspist3/00ALLINONE.c` and `hspist3/edmd_core/*`, and the cluster build is pinned to a commit. One changed byte there and every later build differs from the registered one. If you believe an engine change is needed, write it as a request in the report. Do not make it.
2. **No builds that write into `hspist3/` or any `experiments_*` folder.** The one compile check allowed (item G below) builds into `PHYSICALINFORMATIONDYNAMICS/_m0_scratch/`, and you delete that folder before you finish.
3. **No KOA, no SSH, no cluster, no network.**
4. **Git:** stay on `main`. No branches, no worktrees, no stash, no push. There are about 23 uncommitted local changes on `main` that belong to the HardDisks work — never stage them. Stage by explicit path (`git add PHYSICALINFORMATIONDYNAMICS/<file>`), never `git add -A` or `git add .` from the repo root. Commit-message prefix: `pid:`. Two commits in total for M0 (see §5).
5. **Do not start long runs.** You may run an existing hspist3 binary briefly to confirm what a key or flag does; nothing longer than a minute, nothing that writes into an `experiments_*` folder.
6. **Output of M0 is ONE file:** `PHYSICALINFORMATIONDYNAMICS/M0_AUDIT.md`. No other markdown, no README, no code, no `src/` scaffolding. One file.
7. **Evidence rule:** every statement about the code carries `path:line` and, where it matters, the verbatim signature or line. Where you did not find something, write "not found" — do not infer. Where you are unsure, write "unsure" and why. No superlatives. The report will be checked line by line against the source.

## 3. Why M0 exists

Demo 0 (the next milestone — not yours yet) needs: an SDL2 window that draws a few hundred large disks, a centerline graph, a polyline/arc hull and two resource patches, with overlay toggles; a headless mode writing the same state to disk; deterministic seeded replay. Demo 1 after it needs: small hard disks colliding with each other and with a *static polyline/arc boundary* — not an axis-aligned box — with an energy/momentum check.

The blueprint assumes the HardDisks SDL2 viewer and the EDMD core can be reused for this. Nobody has checked. That is all of M0: find out what is actually there, what can be linked unmodified, what needs an adapter, and what must be new, so that the Demo 0 brief can be written against facts instead of assumptions.

## 4. The audit — answer every item in `M0_AUDIT.md`, same letters and numbers

Open the report with `git rev-parse HEAD`, `git status --short` (verbatim), and the date.

### A. Build and entry point

- A1. How is `hspist3/00ALLINONE` built: which script or Makefile, compiler, flags, SDL2 link line. How a `bin_<shorthash>/` folder (e.g. `hspist3/experiments_gen3_d52_261014/bin_a71f394/`) is produced and what stamps the hash.
- A2. Translation units: is `00ALLINONE.c` one unit that `#include`s `edmd_core/*`, or are they compiled separately? List every file under `hspist3/edmd_core/` with line count and a one-line purpose. List any other source the build uses.
- A3. Where CLI flags are parsed. A table of every flag touching geometry, dynamics, seeding, GUI or output — at least `--mode --engine --edmd-acc --seed-drift-order --kbt1 --particles --particles-boxes --particle-radius --l0 --height --num-walls --wall-positions --wall-mass-factors --wall-hold-steps --fixed-dt --seed --show-simulation --gui-time-scale` plus any others you find. For each: type, default, what it changes, `path:line`.

### B. Scene model

- B1. The data structures for particles, walls/dividers, pistons, springs, heat-bath walls: struct names, fields, `path:line`.
- B2. Which boundary primitives exist: outer box (fixed? axis-aligned?), vertical dividers (movable, mass factor, spring, hold steps), piston (prescribed trajectory — which laws?), anything else. State plainly whether any non-axis-aligned or curved wall exists anywhere in the code. Expected answer: none. Confirm or refute from the collision functions, not from comments.
- B3. Is there a "scene" or "compartment" object that could hold a second, independent domain in the same process, or are geometry and particle arrays file-scope globals? Count the file-scope globals that carry simulation state and list the main ones.

### C. Dynamics core — this decides the project; be exact

- C1. Is gen3 a pure event-driven scheduler (event queue, predicted collision times, invalidation after each collision) or a time-stepped scheme with collision handling inside a step? `--fixed-dt` and `--seed-drift-order=drift-first` suggest the latter or a hybrid. Describe the actual main loop with `path:line`. Same for gen2, briefly, and what `--edmd-acc` switches.
- C2. The collision/event functions: particle–particle, particle–outer wall, particle–divider, divider–divider if any, piston–particle, spring update, heat-bath thermalisation. Signatures verbatim with `path:line`.
- C3. Clock ownership: who advances simulation time; how the GUI frame loop interleaves with physics (`--gui-time-scale`); whether headless and GUI paths run the identical physics code path. This matters because the GUI must never be the sole source of results.
- C4. Callability: is there a function boundary like `step(state, dt)` or `advance_to(t)` that an external program could call with its own state, or is the physics inside `main()`'s loop on globals? Give the exact boundary if one exists.
- C5. Neighbour search: cell lists present? How are cells tied to the rectangular box? A polyline boundary would have to live with this.
- C6. The energy ledger and the event hashes used by the Stage 1 regression: where computed, what they cover. These are the checks that must still pass after anything this project ever does.

### D. Renderer

- D1. SDL2 init, window, render loop: `path:line`. What is drawn (particles, walls, the density / speed / v_x histograms; keys 1–5 and I). Is drawing a function of a state struct, or interleaved with physics?
- D2. Could the renderer draw, unmodified / with an adapter / only by edit — answer per item: (a) a second disk population with its own radius and colour, (b) a polyline/arc set, (c) a scalar field as a heat overlay.
- D3. The histogram/binning code: reusable for a coarse-graining overlay? Where.

### E. Outputs

- E1. What a headless run writes: files, formats, location, what triggers them. Is there a full state dump (positions, velocities at a time) or only summaries?
- E2. The `00_COMMAND.md` convention and any run-manifest convention (seed, commit hash, parameters), so PID can follow the same.

### F. Regression and the do-not-touch set

- F1. The exact file set covered by the Stage 1 / Part A SHA-256 manifests: paths of the manifests and the files they list. This defines the frozen set precisely.
- F2. Existing tests or regression scripts (event-hash comparison, "Rule 4: gen2 unchanged"): where, how to run, how long.

### G. One bounded compile check (optional, at most 30 minutes)

Try to compile the files under `hspist3/edmd_core/` as objects *without* `00ALLINONE.c`, same flags, into `PHYSICALINFORMATIONDYNAMICS/_m0_scratch/`. Report whether it compiles and which undefined symbols appear at link time — those are the globals and functions the core expects from the main file. Delete `_m0_scratch/` afterwards. If `edmd_core` is header-only or is `#include`d into the main unit, say so and skip the check.

### H. Decision — yours to recommend, Chris's to make

For Demo 0 (growth + viewer, no particles) and Demo 1 (frozen hull + small hard disks) separately, grade each option in one paragraph tied to A–G:

- (a) a new C executable in `PHYSICALINFORMATIONDYNAMICS/` that links `hspist3/edmd_core` objects unmodified and reuses the SDL2 renderer through an adapter;
- (b) a new C/SDL2 executable in `PHYSICALINFORMATIONDYNAMICS/` that copies nothing: a small renderer and a small stepped collision loop for polyline boundaries, with hspist3 used only as reference and as a numerical cross-check on the box case;
- (c) a Python prototype (numpy, shapely for the swept union, pygame or matplotlib for the picture) for Demo 0, deciding Demo 1's language later.

Then: the three biggest risks you saw, and every assumption the blueprint makes about hspist3 that turned out false.

## 5. Commits and stop condition

Before starting the audit, commit the two files already in this folder:

```
git add PHYSICALINFORMATIONDYNAMICS/00_M0_BRIEF.md PHYSICALINFORMATIONDYNAMICS/Physical_Information_Dynamics_Engine_Blueprint_v2_0_2026-10-15.md
git commit -m "pid: sealed folder, M0 brief and blueprint v2.0"
```

When the audit is written, end `M0_AUDIT.md` with `git status --short` again — beyond the pre-existing HardDisks changes you saw at the start it must show only `PHYSICALINFORMATIONDYNAMICS/M0_AUDIT.md` — then:

```
git add PHYSICALINFORMATIONDYNAMICS/M0_AUDIT.md
git commit -m "pid: M0 read-only audit of hspist3 SDL2/EDMD for the reuse decision"
```

Do not push. Print the two commit hashes and the audit's line count, and stop. Do not start Demo 0, do not create `src/`, do not write a README, do not propose a name. The Demo 0 brief is written after the audit is reviewed.

If the audit cannot be finished in roughly one working session, write what you have with the open items marked "not done", commit that, and stop. A partial factual report is worth more than a complete guessed one.
