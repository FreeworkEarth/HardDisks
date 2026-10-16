# Demo 0 brief — Physical Information Dynamics (PID): growth and viewer, no particles

For: the Claude Code session working in `PHYSICALINFORMATIONDYNAMICS/`
From: Chris and the Cowork planning chat, 2026-10-15
Builds on: `00_M0_BRIEF.md`, `M0_AUDIT.md` (commit 309fe872), decision of 2026-10-15: option (b) for Demo 0, option (d) for Demo 1.
Scope: Demo 0 only, in three deliveries (0a, 0b, 0c). Stop after each.

Repository root: `/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks`. Paths below are relative to `PHYSICALINFORMATIONDYNAMICS/` unless they start with `hspist3/`.

## 1. What Demo 0 is

A new, small C/SDL2 program that grows a 2D network of *structural disks* from one seed disk in the middle of a rectangular box, draws it, and writes its state to disk so the growth can be replayed and tested without the window. No EDMD, no small particles, no hspist3 code — the audit (M0_AUDIT.md §H) found the renderer cannot be linked and Demo 0 needs no engine.

Three deliveries, each one working program state that Chris looks at before the next starts:

- **0a — unbiased growth.** Seed disk, tips that extend by placing tangent disks with a random heading kick, branching, a material budget, box walls. Window with toggles, headless mode, seeded replay, tests.
- **0b — resource-guided growth.** Two food patches, a cue field, tips biased toward it, branches that reach food thicken, the others thin and vanish, material returned to the pool.
- **0c — the exact hull.** The channel boundary (lumen) swept around the network as exact line segments and circular arcs, with an area test. This is Demo 1's input.

Demo 1 (small hard disks inside the hull, fork of `edmd_gen3.c`) is a separate brief, written after 0c.

## 2. Rules

All of `00_M0_BRIEF.md` §2 still holds: write only inside `PHYSICALINFORMATIONDYNAMICS/`; `hspist3/` and every `experiments_*` folder read-only; no KOA, no SSH, no network; stay on `main`, no branches, worktrees, stash or push; stage by explicit path only (`git add PHYSICALINFORMATIONDYNAMICS/<file>`), never `-A` or `.`; the HardDisks session's uncommitted changes are never staged; commit prefix `pid:`; one session, no subagent fan-out; if `git commit` fails on `index.lock`, wait and retry, never remove the lock.

Added for code:

1. **No libc `rand()`/`srand()`.** One explicit generator (splitmix64 seeding a xoshiro256** or xorshift64*, your choice, written in `pid_rng.c`), seeded only from `--seed`. Every random draw in the program goes through it. This is what makes replay a test instead of a hope.
2. **`-ffp-contract=off` from the first compile** (same reason as `hspist3/Makefile:34-40`), `-std=c11 -O2 -Wall -Wextra`, SDL2 via `pkg-config sdl2`. Zero warnings is the bar.
3. **The window draws; it never owns state.** The GUI and the headless path call the same `step()` function on the same `World` struct. Pause, single-step and reset in the window call that function; nothing else in the GUI touches the world. Every number in a test comes from the headless path.
4. **Outputs go to `runs/<name>/`**, which is git-ignored by a `.gitignore` you create inside `PHYSICALINFORMATIONDYNAMICS/` (also `bin/`, `*.o`). Nothing is ever written to the cwd implicitly, and nothing to `hspist3/`. Committed are: source, Makefile, tests, one screenshot per delivery in `demo0/shots/` (PNG, under 300 KB), and the run manifest of the run that made the screenshot (text).
5. **Doubles are printed with `%.17g`** in state files so a file reads back to the same bits. Replay tests compare files with `cmp`.
6. **Keep it small.** Fewer than ten source files, one Makefile, no dependency beyond SDL2 (SDL2_ttf optional, behind `--font=<path>`; without it the status line goes to the terminal). No Python in the build. Python is allowed in `tests/` only, for the 0c area reference (shapely), and the test must skip cleanly if shapely is missing.
7. **Conventions are written down once**, in `demo0/pid.h` as comments and in every state file header (§4). The blueprint requires one unambiguous geometry convention in run metadata.

## 3. The model (fixed for Demo 0; parameters are flags)

**Units.** World length unit is σ, the diameter of the small transport particles of Demo 1 (so a channel of width 3 σ fits three of them). Areas in σ². Time is the growth tick, an integer. The viewer maps σ to pixels with a scale; the state file never contains pixels.

**Box.** `[0, W] × [0, H]`, default 200 × 120 σ. Hard boundary: a disk whose centre would come within its radius of a wall is not placed.

**Structural disk.** `{id, x, y, r, parent, kind ∈ {seed, body, tip}, born (tick), alive}`. The seed has `parent = -1`.

**Placement is tangent.** A new disk is placed touching its parent at exactly one point: `centre_new = centre_parent + (r_parent + r_new)·(cos θ, sin θ)`. Chris's convention of 2026-10-15: disks touch, they do not overlap. **Non-adjacent disks must never overlap** (`|c_i − c_j| ≥ r_i + r_j` for every pair that is not parent–child); a placement that would violate this is rejected. Consequence: at placement the disks are a non-overlapping set and the material they hold is exactly `Σ π r_i²`.

**Tip.** `{disk id, heading θ, alive}`. The seed starts with `K` tips at evenly spaced headings (default K = 6, plus a random rotation of the whole set so seeds differ).

**Growth tick.** For every live tip, in id order (so replay is deterministic), the policy draws an action and the executor applies it:

- `extend(θ', r')`: θ' = θ + κ·ξ with ξ ~ N(0,1) (Box–Muller from the one generator), r' = max(r_min, r·(1 − τ)); place the disk tangent at θ'; if rejected (box or overlap), redraw θ' up to `n_try` times, then the tip dies.
- `branch(θ1, θ2)`: as extend, plus a second tip on the new disk at heading θ' ± Δ (sign from the generator), with probability `p_branch` per tick.
- `stop`: the tip dies (when r' < r_min, or budget exhausted).

Defaults, to be tuned by you for a readable picture and recorded in the manifest: seed radius R0 = 6, κ = 0.35 rad, τ = 0.03, r_min = 1.0, p_branch = 0.08, Δ = 0.6 rad, n_try = 8, budget M = 3000 σ². Every one is a `--name=value` flag.

**Policy as a function.** `typedef Action (*Policy)(const World*, const Tip*, Rng*, const Params*);` The unbiased policy (0a) and the guided policy (0b) are two functions of this type, chosen by `--rule=unbiased|guided`. The executor (legality, material reservation, placement) is one function shared by all policies. This is the plug-in boundary the blueprint asks for (§2): every later rule — empowerment, IB, DIB — is another function of this type, and only this type.

**Material ledger.** `M_total` (budget), `M_used = Σ π r_i²` over live disks minus the overlap of adjacent pairs (zero at placement; nonzero only after 0b thickening, see below), `M_pool = M_total − M_used`. A placement is charged `π r'²`; it is rejected when `M_pool < π r'²`. The ledger is recomputed from the disks at every state write and must equal the running sum to 1e-9 relative (test).

**0b additions.**

- Two patches `{id, x, y, r_p, s_p, active}`, defaults at (0.2 W, 0.5 H) and (0.8 W, 0.5 H), r_p = 5, s_p = 1, both active.
- Cue field, closed form, no PDE: `c(x) = Σ_p s_p · exp(−|x − x_p| / λ)`, λ = 40 σ. Its gradient is closed form too.
- Guided policy: φ = direction of ∇c at the tip; θ' = θ + κ·ξ + β·wrap(φ − θ), β = 0.3. **With β = 0 the guided rule must reproduce the unbiased rule byte for byte** (test).
- A tip whose new disk intersects a patch is *fed*; the chain of disks from it back to the seed (parent pointers) is a *fed path*; the tip then stops.
- Adaptation, every `T_adapt` ticks (default 10): every disk on a fed path gets `r ← min(r_max, r·(1 + γ))`, γ = 0.02, r_max = 2·R0; every other body disk gets `r ← r·(1 − δ)`, δ = 0.01; a disk whose r falls below r_min/2 is removed, with its downstream subtree and tips. Material: thickening is charged to the pool and refused when the pool is empty or when the larger disk would overlap a non-adjacent disk; thinning and removal return material to the pool. Adjacent (parent–child) disks may come to overlap after thickening; their overlap is the circle–circle lens area, computed exactly and subtracted in `M_used`. That is the only overlap allowed anywhere.
- This rule is the dumb local baseline on purpose. It must *not* know where the patches are except through `c(x)`, and it must not draw the straight line seed→patch. Whether the direct connection emerges is what the picture shows; if it does not emerge, that is a result, not a bug to fix by cheating.

**0c — the hull.** For every parent–child edge (i, j) the *capsule* is the convex hull of the two disks (two external tangent segments and the two outer arcs). The lumen is the union of all capsules over all live edges; its boundary is the hull: a set of closed loops of primitives

- segment `{x1, y1, x2, y2}` with the lumen on the left when walking from 1 to 2;
- arc `{cx, cy, R, a0, a1}` counter-clockwise from a0 to a1, lumen on the left.

Computed exactly in doubles: tangent points from the disk pairs, intersections between every pair of primitives of different capsules, keep the pieces outside every other capsule, chain into loops. O(E²) is fine at E < 2000. The lumen area is computed from the hull by the shoelace formula plus circular-segment terms. Report lumen area, hull perimeter and `M_used` side by side; the budget keeps using `M_used` in Demo 0, which one Demo 3 uses is decided then.

## 4. The state file — write and review this before the growth rule

The contract between Demo 0, Demo 1 and every test. Plain text, one record per line, `#` comments, sections in this order:

```
[header]
format 1
unit sigma
geometry tangent-disks; lumen = union of per-edge capsules (convex hull of two disks)
seed 9300
rule unbiased            # or guided
tick 240
box 200 120
budget 3000
m_used 2871.4...         # %.17g
m_pool 128.5...
params kappa=0.35 tau=0.03 rmin=1 pbranch=0.08 delta=0.6 ntry=8 K=6 R0=6 [beta= lambda= Tadapt= gamma= deltathin= rmax=]
build <git short hash of PHYSICALINFORMATIONDYNAMICS at build>  <cflags>
[disks]
# id x y r parent kind born alive
0 100 60 6 -1 seed 0 1
1 ...
[tips]
# disk heading alive
[patches]          # 0b
# id x y r s active
[hull]             # 0c
# loop <n>
# seg x1 y1 x2 y2
# arc cx cy R a0 a1
[areas]            # 0c
lumen_area ...
hull_perimeter ...
```

`pid_state.c` writes and reads it; a round trip (read → write) must be byte-identical (test). A run writes `state_<tick>.txt` every `--write-every` ticks and `final.txt`, plus `manifest.txt`: argv, seed, all parameters with their values, the build line, the commit hash of the folder's last commit, and the SHA-256 of every state file.

## 5. The window

SDL2, one window, world-to-screen scale from the box size. Draws the state and nothing else:

- `D` toggle structural disks (filled, alpha from the slider; tips in another colour; fed paths in a third colour in 0b);
- `G` toggle the centreline graph (segments between parent and child centres);
- `H` toggle the hull (0c);
- `C` toggle the cue field as a heat overlay on a coarse grid (0b; the audit's `rm_draw` technique — a streaming texture — copied by idea, not by text);
- `P` patches;
- `[` `]` disk transparency down/up (Chris's transparency slider);
- `Space` pause/run, `N` one tick, `R` reset to the same seed, `S` write a BMP into the run folder, `Esc` quit;
- status line (terminal, or on screen with `--font`): tick, live tips, disks, M_used, M_pool.

Command line: `pid_demo0 --seed=S --rule=unbiased|guided --steps=N --out=runs/<name> [--headless] [--write-every=k] [--font=path] [parameters]`. With the window, the run stops at `--steps` and leaves the window open. Headless writes the same files and exits.

## 6. Tests — `make test`, all headless, all green before each delivery

1. **Replay:** two runs with the same seed and flags → every state file `cmp`-identical; a different seed → `final.txt` differs. **GUI = headless:** a window run with `--steps=N` writes the same `final.txt` as a headless run (the window must not change the draw order).
2. **Invariants, checked on every written state:** no non-adjacent overlap; every disk inside the box; every child tangent to its parent at placement (`| |c_i − c_j| − (r_i + r_j) | ≤ 1e-12·(r_i + r_j)`), tested in 0a where no adaptation runs; parent ids form a tree rooted at the seed; dead tips have no further children.
3. **Ledger:** `M_used` recomputed from the disks equals the running sum to 1e-9 relative; `M_used + M_pool = M_total` to 1e-9 after every tick, including 0b thinning and removal.
4. **State round trip** byte-identical.
5. **0b:** `--rule=guided --beta=0` gives `final.txt` `cmp`-identical to `--rule=unbiased` at the same seed. With both patches active and default β, at least one fed path exists by tick 400 on the default seed (if not, tune defaults, record the tuning).
6. **0c hull:** (a) every loop closed (end of each primitive equals start of the next to 1e-9); (b) every disk centre inside the lumen; (c) lumen area from the hull equals the area of a raster of the capsule union at 2, 4, 8, 16 px/σ, with the error decreasing with resolution and the finest within 0.1 %; (d) `tests/hull_reference.py` (shapely, optional): the union of the capsules as polygons with `quad_segs` 16, 64, 256 — areas converge toward the C number from below; (e) reading the edges in reversed order gives the same area to 1e-12 and the same set of primitives up to order.

## 7. Deliveries and stop condition

Each delivery: the tests green (paste `make test` output in the commit message body), one screenshot in `demo0/shots/`, its manifest committed, one commit `pid: demo0a ...` / `pid: demo0b ...` / `pid: demo0c ...`, then **stop** and print the hash and the command Chris runs to see it. Chris looks before the next delivery starts.

Suggested order inside 0a: `pid.h` and `pid_state.c` first (the contract), then `pid_rng.c`, then growth + executor + headless loop + tests, then the window. If the window is the last thing done in 0a, that is right.

Not in Demo 0: moving structural disks, 3D, particles, any information-theoretic rule, a name for the project.

If a delivery cannot be finished in one session, commit what passes, mark what does not in the commit message, and stop.
