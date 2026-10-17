# Demo 1 brief — PID: small hard disks inside the frozen hull (fork of gen3)

For: a Claude Code **cloud** session (headless), then the Mac session for the viewer step
From: Chris and the Cowork planning chat, 2026-10-15/16
Builds on: `M0_AUDIT.md` §H option (d) and §C5/C6 (commit 309fe872); `01_DEMO0_BRIEF.md` §3–4, §8–9; Demo 0 complete at commit cae73936 (0b.1).
Scope: Demo 1 only, three deliveries (1a fork, 1b boundary primitives, 1c viewer). Stop after each.

Repository root on the Mac: `/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks`; in the cloud, the clone's root. Paths below are relative to `PHYSICALINFORMATIONDYNAMICS/` unless they start with `hspist3/`.

## 1. What Demo 1 is

Demo 0 grows a network and writes its lumen as exact segments and arcs (`[hull]` in `final.txt`). Demo 1 puts small hard disks inside that frozen lumen and lets them collide — with each other and with the hull — using the validated event-driven engine of HardDisks, **forked** into this folder and extended with a segment/arc wall primitive. Nothing in `hspist3/` changes.

The audit established (M0_AUDIT.md §C1–C6, §H): gen3 is a pure event-driven scheduler (one binary heap keyed (t, type, a, b), lazy invalidation by collision counters, per-disk clocks, a cell grid with crossing events); its only geometry is the axis-aligned box with vertical slabs and faces, hard-coded as a side index 0–3 in `switch` statements; it links alone and references nothing in the driver; its brute-force schedule audit, contact audit, event hash and ledger live inside `edmd_gen3.c` and are switched on by its API. So the fork is: copy, add one wall class, keep every check.

Deliveries:

- **1a — the fork and the box reference.** `demo1/engine/` holds a copy of `hspist3/edmd_core/edmd_gen3.c`, `edmd_gen3.h` and `edmd.h` at HEAD. Unchanged except the header comment recording source commit and blob hashes. A test builds the M1/M2 harnesses against the fork *and* against the unmodified originals and requires identical output (hashes, final states, health) — the box reference of §H (a). No SDL2 anywhere in 1a/1b.
- **1b — segments and arcs.** New event types for static walls given as the Demo 0 hull; prediction, reflection, ledger, validator, audit, cell lists. Tests on a corridor, a tilted corridor, a Y junction, and a real Demo 0 hull. Headless `pid_demo1` executable: loads a Demo 0 `final.txt`, seeds particles, runs, writes states and measurements.
- **1c — viewer (Mac only).** The Demo 0 window draws the particles inside the hull; `--view=network` plus particles. Done on the Mac after 1a/1b are merged.

Not in Demo 1: inlets and outlets (Demo 2), moving walls or growth while particles run (Demo 3/4), a second particle radius, any policy.

## 2. Rules

`00_M0_BRIEF.md` §2 and `01_DEMO0_BRIEF.md` §2 hold, with these changes for a cloud session:

1. **Branch, not main.** The cloud session works on its own branch (its default `claude/...` name is fine) from current `origin/main`, touches only `PHYSICALINFORMATIONDYNAMICS/demo1/` (and `.gitignore` inside `PHYSICALINFORMATIONDYNAMICS/` if needed), commits with the `pid:` prefix by explicit path, and pushes the branch. It never merges, rebases or pushes `main`. Chris merges on the Mac after reading the delivery; a delivery that touches any path outside `demo1/` is not merged.
2. **`hspist3/` is read-only**, in the clone as on the Mac. The fork reads it once (the copy in 1a) and the box-reference test compiles the originals with `-o` into `demo1/bin/`. No write, no `make` in `hspist3/`.
3. **Network only for `pip`/`apt` of test dependencies** (shapely for the hull reference, nothing else). No KOA, no SSH.
4. **Flags:** `-std=c11 -O2 -Wall -Wextra -ffp-contract=off`, zero warnings under `-Werror`. The box-reference comparison needs both sides built with identical flags on the same machine (M0_AUDIT.md §F2; `hspist3/Makefile:34-39`).
5. **One generator** (`demo0/pid_rng.c`, reused by copy or by relative include) for seeding particles; libc `rand()` never. The fork's heat-bath code paths are refused as in gen3 (`edmd3_create` refuses `heatbath_enabled`).
6. **Outputs to `demo1/runs/<name>/`** (git-ignored). Committed: sources, Makefile, tests, one manifest per delivery, small PNGs only from 1c. No state dumps, no `.o`, no binaries.
7. **Evidence rule as in M0:** every claim in a commit message about gen3's code carries `path:line` in the fork; every number has the command that produced it. Where the fork deliberately diverges from gen3, the divergence is a comment at the site starting `/* PID: `.
8. **Keep the names.** The fork keeps the `edmd3_*` API names. The fork and the unmodified originals are never linked into one binary; the box-reference test runs two executables and compares their output files.
9. **Symbol-level diff.** 1a ends with `diff` of the fork against the originals being empty except the header comment; 1b's commit message lists every function added or changed.

## 3. Units and conventions (fixed)

- Demo 0 lengths are σ (the small particle diameter). gen3 works in px and internal time with **24 px = 1 σ, 24 units = 1 σ-time, unit disk mass, kT = 1** (`hspist3/edmd_core/edmd_gen3.h:22-23`). Demo 1 converts the hull at load: every coordinate × 24. The particle radius is **R = 12 px** (0.5 σ), the hspist3 default. The box is the Demo 0 box × 24 (200 × 120 σ → 4800 × 2880 px). Cell width: the gen3 default `EDMD3_DEFAULT_CELL_PX` = 32 (an integer number of px ≥ 2R, `edmd_gen3.c:1116-1117`).
- The hull stays **inside the box**: the cell grid covers [0, boxW] × [0, boxH] (`edmd_gen3.c:1136`), the four outer walls remain as the outer safety boundary (§H (d) 7), and `edmd3_load` refuses a disk outside the box. A hull that touches the box edge is a Demo 0 problem, not Demo 1's.
- **Lumen convention from Demo 0:** segments have the lumen on the left walking 1 → 2; arcs run counter-clockwise from a0 to a1 with the lumen on the left. The particle-accessible region is the lumen eroded by R: a particle centre must stay at distance ≥ R from every primitive and inside the lumen. Seeding uses the capsule list (disks + parent edges from `[disks]`), which is the same region described differently; the engine uses the primitives.
- **One radius, unit mass**, no species, no gates, no heat bath, no pistons or dividers in Demo 1 (the fork keeps that code; the Demo 1 driver does not use it).
- **Manifest** per run: argv, the Demo 0 `final.txt` path and its SHA-256, px/σ = 24, R, N, T, seed, cell width, the fork's source commit and blob hashes, build line, event hash at the end, `clean`, and the SHA-256 of every state file written.

## 4. The engine change (1b), by the audit's list (M0_AUDIT.md §H (d) 1–8)

1. **Event types.** Append `T_SEG = 6`, `T_ARC = 7` after `T_PAIR = 5` (`edmd_gen3.c:89`), with the primitive index in `b`. A box-only run of the fork then pushes no new events, so its tie order (`ev_less`, 186-191) and event hash (1283) equal unmodified gen3's — that is the 1a regression, and it must still pass after 1b with an empty primitive list.
2. **Static geometry.** The fork gets `edmd3_set_walls(EDMD3*, const Seg*, int nseg, const Arc*, int narc)` (names free), called once after `edmd3_create` and before `edmd3_load`. Primitives are immutable for the life of the state. Per cell, a precomputed list of the primitives that can come within R of that cell (conservative bounding-box test is fine), consulted wherever `wall_candidate` is used today: `predict_all` (498), `exec_cross` (786-787), load (1248).
3. **Prediction with absolute positions** (`x = cx·w + xi`, as `predict_obj` does at 482):
   - segment interior: the signed distance to the line minus R, over the normal velocity — a linear gap/speed rule like `wall_gap` (439-446), valid only while the foot of the perpendicular lies within the segment at contact time;
   - segment endpoints: a stationary partner of radius 0 at the endpoint, contact distance R — `pair_rule` (407-419) with the partner's velocity zero;
   - arcs: a stationary partner at the arc centre with contact distance R_arc − R (particle inside the arc's circle, convex lumen side) or R_arc + R (outside, concave side) — the same quadratic as `pair_rule`, choosing the root that corresponds to approaching the wall, and accepted only if the contact point's angle lies in [a0, a1].
   The first-contact rule for a segment is the minimum over interior and endpoints; a contact that is simultaneously interior and endpoint (the foot exactly at an endpoint) is one event, not two.
4. **Reflection** about the contact normal: `v' = v − 2 (v·n) n` with the unit normal from the contact point, changing `vx` and `vy` together (gen3's `exec_wall` flips one component, 742-745). Energy is unchanged exactly up to rounding; the impulse `m (v' − v)` is booked as external, per axis, into the ledger's `J` as the outer walls do (742, 744), and per primitive into a new array `J_prim[nprim]` (the per-wall records are fixed at 4 today: `J_wall[4]`, `edmd_gen3.h:104`). Bump the disk's collision counter (`cnt`), as every velocity change does (728, 747, 845).
5. **Validator and repair.** `local_check`/`full_check` (567-671) get a primitive class: the distance of a disk centre to every candidate primitive must be ≥ R − tol, with tol derived as for the box (`tol_face` from ulp(boxW), 808) but from the primitive's extent. `reflect_into_box` (609-619) gets a sibling `reflect_off_primitive` that pushes a disk found inside the wall back along the normal and counts it (`overlap_repair`-style counter, new name). `clean` (1538-1543) includes the new counters.
6. **Audits.** `schedule_audit` (995) gets a brute-force block for segments and arcs next to the four-wall block (1050-1056), recomputed from absolute positions, and `contact_audit` (675) a branch for the new types. The harness pattern A/B/C (audits on, off, again; equal hash and state; `gen3_m1_harness.c:3-6`) is reused unchanged on every Demo 1 test geometry.
7. **Cell crossings.** On `exec_cross`, predict against the primitives newly reachable from the new cell (the per-cell list difference), as 786-787 does for the box walls. No bands: bands are x-column constructs for vertical bodies (522-554) and do not apply.
8. **Checkpoints** of the fork are incompatible with gen3's (`sizeof(EDMD3)` is in the header, 1414, 1444); the fork writes its own magic. Fine.

Also: the three file-scope hook statics (`coll_hook_S/fn/ctx`, 704-706) move into the state struct in the fork (one `/* PID: */` comment), so two states in one process are possible later.

**Concave corners** (a Y junction's inner corners, where a disk can touch two primitives at once) are the hard case: they interact with the tie-break key (186-191) and the stagnation guard (`same_limit`, 1157). The "tie" cell of the M1 harness is the template; test 6(d) below is the test.

## 5. The Demo 1 driver (`pid_demo1`, headless)

`pid_demo1 --state=<demo0 final.txt> --N=n --T=<σ-time> --seed=S --out=demo1/runs/<name> [--R=12] [--cell=32] [--write-every=<σ-time>] [--audit-every=k] [--kbt=1]`

- Loads `[disks]`, `[hull]` and `[header]` of the Demo 0 file; refuses a file without `[hull]`.
- Converts to px; builds the eroded-lumen test from the capsules; seeds N particle centres by random sequential insertion inside the eroded lumen with no overlaps (refuses after a bounded number of attempts and says so); velocities Maxwell–Boltzmann at kT, total momentum removed, kinetic energy rescaled to N kT exactly (as the harness does).
- Creates the fork state, sets the walls, loads, runs `edmd3_advance_to` in steps of `--write-every`, and at each stop writes `pstate_<t>.txt` (positions, velocities, `coll_count`, `%.17g`) and appends one line to `measure.csv`: t, KE, event counts by type, `clean`, the event hash, and the impulse per primitive summed since the last stop.
- Writes `density_<t>.txt` on a coarse grid (cell = 2 σ, i.e. 48 px): particle count per cell and the Demo 0 structural material per cell (Σ π r² of structural disks by centre). **This is the "one grid, two fields" measurement** the planning chat wants; it is the only analysis in Demo 1.
- Writes the manifest of §3.

## 6. Tests — `make test` in `demo1/`, all headless, all green before each delivery

1a:
1. **Fork = original, box only.** Build `gen3_m1_harness` and `gen3_m2_harness` twice (against the fork, against the originals in `hspist3/edmd_core/`), run `audit --quick` and `audit` for M1 and `audit --quick` for M2 as separate processes; the outputs must be byte-identical except lines that print a build string. Also build the fork alone with `-c` and check `nm -u` shows only libc/libm.
2. **Source diff.** `diff` fork vs originals is empty except the header comment.

1b (each geometry run with the A/B/C audit pattern, `clean=1` required, at η ≈ 0.05 and η ≈ 0.30 of the lumen area, T = 400 σ-time unless stated):
3. **Box-only with an empty primitive list** reproduces test 1's hashes (the new code is inert when unused).
4. **Straight corridor**, axis-aligned, two segments inside the box: relative KE drift ≤ 1e-12; no centre ever closer than R − 1e-9 px to a primitive at any stop; schedule audit: 0 missing, 0 extra; contact audit max gap ≤ 1e-12 px; **pressure:** mean impulse per σ-time per σ of wall on the two long segments gives Z = P·A/(N kT) within 3 % of Kolafa–Rottner at the run's η (this is Paper 1's pressure check reused; at η = 0.05 KR gives 1.108).
5. **Tilted corridor** (the same corridor rotated by 30° about its centre, still inside the box): identical checks; Z within 3 % of the axis-aligned value. Reflection about a general normal is what this tests.
6. **Y junction** made of segments and arcs with two concave corners: (a)–(c) as in test 4; (d) a constructed start with one disk aimed exactly at a concave corner and one aimed at an arc–segment joint, run with the schedule audit after every event: no missing/extra events, no stagnation, `clean=1`.
7. **A real hull:** `demo0` run `--preset=mat45 --rule=guided --seed=9300 --steps=400` (its `final.txt` is reproducible from the manifest in `demo0/shots/`), N = 2000, T = 200 σ-time: `clean=1`, KE drift ≤ 1e-12, no escape, event rate reported. Every particle centre inside the eroded lumen at every stop (the point-in-capsule test).
8. **Determinism:** two runs of test 7 with the same seed give byte-identical `pstate_*.txt` and `measure.csv`; a different seed differs.
9. **Density files** sum to N per stop and the structural field equals Demo 0's `M_used` to 1e-9 (centres binned, no overlap subtraction — say so in the file header).

1c (Mac): the window run of test 7's configuration writes the same `pstate_*.txt` as headless.

## 7. Deliveries and stop condition

- **1a:** `demo1/engine/` (three files + `README_FORK.md` of ten lines: source commit, blob hashes, what will change), `demo1/Makefile`, `demo1/tests/`, tests 1–2 green. Commit `pid: demo1a fork of gen3 with box reference`, push the branch, print the hash, **stop**. Chris merges and tells the Mac session.
- **1b:** the engine change of §4, `pid_demo1`, tests 3–9 green. The commit message lists every function added or changed in the fork with `path:line`, and the measured Z of tests 4 and 5 with the KR value. One commit (or two: primitives, then driver), push, **stop**.
- **1c (Mac session, after merge):** viewer; one PNG of test 7's configuration at t = 100 σ-time with `--view=network` plus particles. Commit, stop.

If a delivery cannot be finished in one session: commit what passes, mark the rest, push, stop. The box reference (test 1) must never be red at a push; if it is, the delivery is not pushed.

## 8. Amendment after 1a and the cloud pre-check (2026-10-16)

1a is merged on the Mac as `41d2ea07` (cherry-pick of `a79ca2aa`). The pre-check the cloud session wrote is committed beside this brief as `02_DEMO1_BRIEF_PRECHECK.md`. **Its part B is adopted in full**, items B1–B21, and replaces the parts of §3–§6 it corrects. Where it leaves a choice, this section decides. Where this section and §3–§6 differ, this section wins.

**Decisions on the pre-check's open choices**

- **B2/B3 arcs.** Every hull arc is hit from inside its circle; there is no outside case. §4.3's "same quadratic as `pair_rule`, approaching root" is withdrawn: a separate inside rule with the larger root, the unwrapped angle test, arcs with R_arc ≤ R skipped, and the last-primitive guard in `Disk3.last` with its own encoding.
- **B7 hook and checkpoints.** The hook statics stay as they are (§4 "Also" is withdrawn): moving them breaks the box reference that 1a now tests. Checkpoints: `edmd3_checkpoint_write` refuses when primitives are set; with an empty primitive list it is unchanged, so test 3 stays byte-identical.
- **B12 Demo 0 headless.** Adopted as written (stub `view_run`, compile the Demo 0 sources from `demo1/Makefile`, never `make` in `demo0/`).
- **B13/B15/B16 the test-7 hull** is identified by structure (416 disks, 415 edges, 10 loops, 833 segments + 330 arcs, lumen_area 6659.169369 to 1e-6 relative), escape is "inside the lumen and distance to every primitive ≥ R − tol", dead disks are filtered.

**Test changes**

- **Test 2** stays 1a-only, and its "last change to … before the copy commit" line becomes INFO, not PASS/FAIL. On the Mac it failed for `edmd.h` while the content checks (the blob at the copy commit, the byte-identical body) passed: which commit `git log -1 -- path` names depends on clone depth and merge simplification, so it is provenance, not identity. The blob is the identity. In 1b test 2 is replaced by B21 (every changed hunk carries `/* PID: */`).
- **Test 3** = test 1's full comparison (all eight outputs byte-identical, verdict checks, object `cmp`) with an empty primitive list — not only hashes.
- **Test 4 pressure.** Corridor at least 40 σ wide. Z = P·A/(N kT) with A the area between the walls. At η = 0.05 within 3 % of Kolafa–Rottner (1.108). At η = 0.30 the KR deviation is reported, not gated (wall layering shifts the mid-corridor density); gated is straight vs tilted within 2 block standard errors. Gap bounds per B1, relative to the run's tolerance, measured maxima reported.
- **Test 9** compares the binned structural field with Σπr² over alive disks to 1e-9, and prints `M_used` and the lens total beside it (B14).
- **New test 10 — positive controls.** A white-box build (as `gen3_band_edge_ties.c` includes the `.c`) plants three errors and each must be reported: a primitive event deleted from the heap (schedule audit: missing ≥ 1), a disk moved 0.1 R through a segment (validator: count ≥ 1, `clean=0`), a reflection with the normal sign flipped (contact audit or KE drift fails). A detector that cannot fail proves nothing.

**Process**

- One session, no reviewer panels or sub-agent fan-out. Mutation and positive-control tests are the review.
- 1b starts from `origin/main` after Chris pushes it (it then holds `41d2ea07`, this amendment and the pre-check). The old branch is not merged.
- If 1b does not finish in one session, commit the primitives with tests 3, 4, 5, 10 green, push, stop; the driver and tests 6–9 follow.
