# Demo 1 brief: pre-check in a cloud container

Date: 2026-10-16. Checked: `02_DEMO1_BRIEF.md` at commit `6b29aa5b` (origin/main).

Machine: Linux x86-64, 4 vCPUs. The default `cc` is gcc 13.3.0; clang 18.1.3 is also installed. SDL2 is not installed.

This file corrects the brief; it does not replace it. Read the brief first, then this file. Part A is needed for delivery 1a. Part B should be folded into the brief before 1b.

How the check was made:
- Three read-only checkers covered the code citations, a dry run of 1a's tests, and the Demo 0 interface.
- A skeptic then re-ran every finding of each checker. All findings below survived. Severities are the skeptics'.
- Builds and runs went into the session scratchpad only. Nothing in the repository was written.

Line numbers without a path refer to `hspist3/edmd_core/edmd_gen3.c` at `6b29aa5b`. Those sources are identical to `a71f394d`.

## A. Delivery 1a

1. **gcc fails `-Werror` on the unmodified gen3** (blocker).
   - Line 842 (`snprintf(kb, sizeof kb, "D%d", o)` into `char kb[8]`) gives `-Wformat-truncation` at -O0, -O2 and -O3. clang 18 compiles it with 0 warnings.
   - 1a may not edit the file, so the brief's rule 4 (zero warnings under `-Werror`) and rule 9 (fork unchanged) conflict under gcc.
   - **Do:** either use `CC=clang` throughout, or keep gcc and add `-Wno-format-truncation` to the gen3 compile only, identically on the fork side and the original side. Record the choice in the Makefile and the manifest.

2. **`edmd.c` (gen2) must be linked but is not one of the fork files** (blocker).
   - Both harnesses run gen2 as their run D (`tests/gen3_m1_harness.c:7, 152-197`; build lines `:16-17`, `gen3_m2_harness.c:34-35`).
   - Under `-std=c11` on glibc, `edmd.c` does not compile with gcc or clang: `M_PI` is undeclared at 556, 694, 741, 1085, 1108, 1160, 1602 and 1635. With `-D_DEFAULT_SOURCE` clang gives 0 warnings. gcc still gives 4 `-Wmisleading-indentation` warnings (719, 765, 1591, 1592).
   - **Do:** compile the original `hspist3/edmd_core/edmd.c` once into `demo1/bin/edmd.o` with `-D_DEFAULT_SOURCE` (under gcc also `-Wno-misleading-indentation`), and link that one object into both sides.
   - Rule 8 reads correctly as "fork gen3 and original gen3 are never linked together". The unmodified gen2 in both binaries is required. Say so in `README_FORK.md`.

3. **Include trap: a harness compiled in place always uses the ORIGINAL headers** (blocker for 1b; in 1a a silent false green).
   - `#include "../edmd_gen3.h"`, `"../edmd.h"` and `"gen3_state_io.h"` are resolved relative to the harness file before any `-I` or `-iquote`. This was shown with `#error` in a decoy header: the in-place build still compiled.
   - **Do:** stage a mirror under the git-ignored `demo1/bin/`:
     - copy the fork's `edmd.h` and `edmd_gen3.h` to `demo1/bin/hfork/`;
     - copy the unchanged harness files `gen3_m1_harness.c`, `gen3_m2_harness.c` and `gen3_state_io.h` to `demo1/bin/hfork/tests/`;
     - build the fork side from there.
   - Guard: the build fails if `cc -std=c11 -MM demo1/bin/hfork/tests/gen3_m1_harness.c` mentions `hspist3`.

4. **The harness output holds no final states** (trap).
   - Final-state equality is a `memcmp` inside one process, printed only as yes/NO (`gen3_m1_harness.c:199-213`). So the cross-binary comparison covers event hash, counts, health and audits.
   - The hash is FNV-1a over every executed (t, type, a, b) (line 1283), which is strong evidence but not "final states identical".
   - **Do:** state this in the manifest. Add `cmp` of the fork object against the original object: they are byte-identical under gcc and clang at -O2 without `-g`. If final states must be claimed, add a small driver, built against each side, that prints `edmd3_particles()` with `%a`.

5. **Test 1 takes longer than one foreground command may run** (trap).
   - Per side, measured: M1 `audit --quick` about 91 s, M2 `audit --quick` about 55 s, M1 `audit` (full) about 285 s. Both sides in sequence take about 14–15 min, against a 600 s cap on a foreground command.
   - Parallel runs give byte-identical output: audit mode prints no timings and writes no files.
   - **Do:** run the six processes with `make -j4` or in the background, each under `timeout 900`. Check every return code.
   - The harness exits 0 even after `RUN FAILED` or `FATAL`. So grep that no output contains `RUN FAILED`, `FATAL` or a `NO` verdict, and that every "audits do not steer" line reads YES, in addition to `cmp`.

6. **`nm -u` differs by compiler** (trap).
   - gcc emits `__printf_chk`, `__snprintf_chk`, `sincos`, `fmax`/`fmin`; clang emits `printf`, `ceil`, `floor`, `cos`/`sin`.
   - **Do:** link the fork object with a stub `int main(void){return 0;}` and `-lm` only (rc 0 under both), or require `nm -u` of the fork object to equal that of the original object. Never use a fixed list.

7. **Linux hashes are not Mac hashes** (trap).
   - Fluid cell, M1 quick, run A: Linux `24546f3c35d265cf` with 34,860 events; Mac `e06d57c1033cede0` with 34,552 (M0_AUDIT.md:2023).
   - **Do:** compare only fork against original, same compiler, flags and machine. Never compare with committed hspist3 outputs or M0_AUDIT hashes.

8. **The header comment shifts every line number** (nit).
   - A provenance header moves every fork line, so the brief's hspist3 line numbers are off by k in the fork, while rule 7 asks for `path:line` in the fork.
   - **Do:** use a fixed-length header and strip exactly that many lines in test 2, or state k in `README_FORK.md`.
   - Provenance to record: HEAD `6b29aa5b`; last commit touching edmd_core `a71f394d`; blobs `edmd_gen3.c e7d9ce3e`, `edmd_gen3.h 6616f999`, `edmd.h 4ccd5da0`, and `edmd.c b9921b0f` (linked unchanged).

9. **The repository-root `.gitignore` ignores some extensions** (nit).
   - It ignores `*.out`, `*.log`, `*.d` and `*.cmd`. A committed test log with one of those extensions silently drops out of `git add`.
   - **Do:** use `.txt`, and check `git status --short PHYSICALINFORMATIONDYNAMICS/demo1` before committing.

10. **The 1a manifest is unspecified** (nit). Write one containing:
    - compiler and version;
    - every object's exact compile line, including the per-file flags from items 1 and 2;
    - HEAD and the blob hashes;
    - the SHA-256 and wall time of each of the six harness outputs;
    - the object `cmp` result;
    - the `nm -u` and stub-link result.

Measured result of the dry run: with items 1–3 applied, fork-side and original-side outputs are byte-identical for M1 quick, M1 full and M2 quick under gcc and under clang. gcc and clang outputs differ only in the `build:` line. 1a is feasible as briefed once these items are handled.

## B. Delivery 1b: corrections to fold into the brief first

Engine (§4):

1. **Test 4's "contact audit max gap ≤ 1e-12 px" fails for unmodified gen3** (blocker).
   - Over 400 σ-time in the box alone, M1 full gives pair gaps up to 3.95e-12 px and wall gaps up to 1.96e-12 px. 5 of 6 cells exceed 1e-12.
   - **Do:** express the bound relative to the run's derived tolerance (`tol_face`, printed via `edmd3_tolerances`), and report the measured maxima.
   - This also covers test 6's "(a)–(c) as in test 4". Test 4's checks have no letters; read it as "all of test 4 except pressure".

2. **Arcs: the lumen is always inside the circle** (trap).
   - Demo 0 writes arcs counter-clockwise with the lumen on the left (`demo0/pid.h:65-66`), so the "R_arc + R, outside" case never occurs.
   - `pair_rule` cannot be reused for the inside case. It returns "at once" whenever c ≤ 0 (413), and inside the eroded circle c < 0 always holds.
   - **Do:**
     - write a separate inside rule using the larger root (−b + √(b² − vv·c))/vv;
     - fire "at once" only when the disk is outside the eroded circle and moving outward;
     - unwrap the angle test as (φ − a0) mod 2π ≤ a1 − a0 (60 of the 330 arcs in the test-7 hull have a1 > π);
     - skip arcs with R_arc ≤ R. Equality is possible, since Demo 0 removes disks only below r_min/2 = R.

3. **Last-collision guard** (trap).
   - For a disk inside an arc's circle, a later contact on the far side always exists. A general-normal reflection is not exact, unlike the box's sign flip.
   - **Do:** encode the last primitive in `Disk3.last` with a value that cannot collide with -1, j ≥ 0 or -2-o (for example -2-NOBJ-p). Suppress only the "at once" case for that primitive, as `mlast` does (359, 363, 483), and keep the far-side root.

4. **The validator must be signed** (trap). "Distance to every candidate primitive ≥ R − tol" cannot detect a disk that has passed fully through a wall.
   - **Do:** check lumen side plus depth (the box's check is signed, `g[s] = x − R`), or add an inside-lumen test to `full_check`.

5. **Catch-all branches would mis-handle the new types** (trap).
   - Event dispatch ends `else exec_obj(S, &e);` (1288), which would run a T_SEG/T_ARC event as a divider or piston collision.
   - `contact_audit`'s final `else` (691-695) would treat it as an outer wall.
   - The gap readers have fixed sizes (`[2]` at 1361, `[4]` at 1478-1481).
   - **Do:** add explicit branches before both `else`s, and a new reader for primitive gaps.

6. **`wall_candidate` has a fourth use** (trap). The brief lists three; the fourth is the schedule audit at 1061 (missing vs deferred). The audit's wall table `awl` is sized 4·N (1000, 1008).
   - **Do:** model the primitive audit block on 1046-1071, with its own (disk, primitive) table and the same per-cell lists.

7. **Checkpoints copy the whole struct** (trap). `checkpoint_write`/`checkpoint_read` copy the whole `EDMD3` and fix up only the listed pointers (1422-1429, 1453-1470). Every new pointer member (primitive arrays, cell lists, `J_prim`, a moved hook) would be written as an address.
   - **Do:** extend both lists and change the magic, or refuse checkpoints when primitives are set.
   - Moving the hook statics touches 704-709, 730, 749, 847 and 1180, and contradicts the comment at `edmd_gen3.h:176-186`. Dropping the hook in the fork is the simpler option.

8. **Public structs** (trap). Append new fields to `EDMD3_Health`, `EDMD3_Ledger` and `EDMD3_Audit` only at the end, or expose them through new functions. Keep the fork's `edmd.h` byte-identical, because the original `edmd.c` is compiled against it.

9. **Tolerance** (nit, but it bites on real hulls).
   - The box validator uses `tol_wall` (1150; used at 585, 590, 640, 645), not `tol_face` (808 only splits overdue from contact-now). Leave `tol_wall`, `band_margin` and `tol_face` untouched: they feed the M2 hash.
   - Derive the primitive tolerance from coordinate magnitude, never from a primitive's length: the shortest segment is 0.4 px, and real-hull vertex gaps reach 1.7e-12 px.

10. **Ledger** (nit). A general reflection rounds both velocity components. Book `ledSP[0]`, `ledSP[1]` and `ledSE` as `exec_pair` does (726-727), not as `exec_wall` does.
    - For Z: every interior hit on a straight segment pushes along its normal, so the per-segment scalar 2|v·n| is exact. Book endpoint hits separately.

11. **Per-cell lists use closed comparisons** (≤/≥ plus a margin), as `wall_candidate` does (264-271). Constructed starts put disks exactly on cell boundaries.
    - Clamp cell indices to the grid: the mat45 hull comes within 0.0987 σ (2.4 px) of y = 0, which is less than R. Do not refuse such hulls.

Hull and driver (§3, §5, §6):

12. **Building Demo 0 headless** (trap).
    - `demo0/Makefile` needs SDL2. Rule 3 forbids installing it, and running `make` in `demo0/` writes `demo0/bin/.build_id`.
    - **Do:** in `demo1/Makefile`, compile `../demo0/pid_main.c` and `../demo0/{pid_grow,pid_state,pid_rng,pid_sha256,pid_hull}.c` together with a stub in `demo1/` (`#include "pid.h"` / `int view_run(Run* r){ (void)r; return 1; }`). Pass `-DPID_BUILD_ID`, `-DPID_CFLAGS` and `-DPID_SRC_SHA256`. It builds clean under `-Werror` with gcc and clang.
    - Run: `--preset=mat45 --rule=guided --seed=9300 --steps=400 --headless --out=<empty dir>`. `--out` is required, and a non-empty `--out` is refused.
    - Load `final.txt`. `state_000400.txt` has no `[hull]`.

13. **`final.txt` is not byte-reproducible across Mac and Linux** (trap). The last bits differ from tick 50 on, and the file carries a build line.
    - **Do:** identify the hull by structure, not by SHA. Linux matches the Mac commit message:
      - 416 disks, 415 edges, 10 loops;
      - 833 segments + 330 arcs;
      - lumen_area 6659.169369, hull_perimeter 3538.547682, m_used 4451.865832.

14. **Test 9 cannot pass as written** (blocker).
    - `M_used` subtracts the parent–child lens overlaps (`demo0/pid_grow.c:436-446`). Σπr² does not.
    - On the test-7 hull they differ by 2.33 %: 4555.719 against 4451.866 σ², the difference being 103.853 σ² of lenses over 20 pairs.
    - **Do:** compare the binned field with Σπr² over alive disks, to 1e-9 relative, and print `M_used` and the lens total beside it. State the unit (σ²).

15. **Escape test** (blocker).
    - "Inside some capsule by depth ≥ R" wrongly rejects about 0.5 % of legal centres near capsule junctions. It never wrongly accepts one.
    - **Do:** the escape test is "inside the lumen AND distance to every `[hull]` primitive ≥ R − tol". The per-capsule test may be used for seeding only.

16. **Filter dead disks** (trap). `[disks]` lists removed disks too (`alive 0`).
    - **Do:** use alive disks only, and capsules only where both ends are alive. The test-7 file has none dead, so the tests will not catch a missing filter.

17. **Units** (nit). Scale x, y, cx, cy and R by 24; never a0 or a1.

18. **Loops** (nit). The outer loop is not always loop 0. Classify loops by signed area: exactly one positive, the rest holes.

19. **Reflex vertices** (nit).
    - In the test-7 hull, 472 of 503 segment–segment vertices are bit-identical (exact ties between two endpoint events) and 31 differ by up to 7.1e-14 σ.
    - There are no convex corners and no arc–arc corners.
    - **Do:** assert in `set_walls` that every arc end is a tangent joint or a segment endpoint. Run test 6(d) with a bit-identical vertex and with a real-hull vertex that has a gap.

20. **A/B/C driver** (nit). Include primitive events in both of the harness's event counts (`gen3_m1_harness.c:131, 140`).

21. **Test 2 is 1a-only** (nit). In 1b, replace it with "every changed hunk carries a `/* PID: */` marker".

## C. Citation nits (no effect on the work)

- kT = 1 is at `gen3_m1_harness.c:14`, not `edmd_gen3.h:22-23`. (Corrected 2026-10-16 from ":15"; line 15 is the build line.)
- The cell width 32 is defined at `edmd_gen3.h:34`; 1116-1117 are only its checks.
- The flag line is `hspist3/Makefile:40`; 34-39 is the rationale comment.
- The full outer-wall audit block is 1046-1071, not 1050-1056.
- The event hash and the ledger always run; only the audits and the event log are switched by the API.

## D. Confirmed as the brief states

- The enum is at 89, and appending 6/7 keeps every existing tie order.
- The hash (1283) depends only on (t, type, a, b), so an empty primitive list leaves it unchanged.
- The heat bath is refused at 1119.
- Every other cited range (ev_less, wall_gap, pair_rule, predict_obj, predict_all, exec_cross, load, exec_wall, J_wall, the validators, schedule_audit, contact_audit, bands, the hook statics, checkpoint sizeof) is where the brief says.
- Units ×24 with R = 12 px, the box 4800 × 2880 px and the 150 × 90 grid are right.
- The lumen convention matches the writer: 0 of 1163 primitives are on the wrong side.
- `pid_rng.c` needs only `pid.h` and libm, and can be included by relative path.
- Seeding N = 2000 in the test-7 hull is feasible.
