# M0 audit: hspist3 SDL2/EDMD for the PID reuse decision

A read-only reconnaissance of HardDisks, written by the Claude Code session working in `PHYSICALINFORMATIONDYNAMICS/` and following `00_M0_BRIEF.md`. The item letters and numbers are the brief's. Part H recommends; the decision is Chris's.

## State at writing

Date: 2026-10-15 14:04 HST.

`git rev-parse HEAD`:

```
a643c02c4f77dd4b56dc3c52374d54090d98769a
```

`git log --format='%h %s' c37bad5e..HEAD` at writing (every commit since M0 started):

```
a643c02c plan-author decision 57 and a correction to decision 56 recorded (261012 sec. 4.7.135): Report D56-1 accepted; the Stage 2 pilot P2 (N 900, eta 0.72/0.78, G1/G2, two pins, 400 runs): package, Mac rehearsal, KOA lines
8aaf5fef Report D56-1 recorded (261012 sec. 4.7.134): figure 1 with block D; Stage 2 sizing -- E1 survives at eta 0.76 (Q 19), pinned lattice starts feasible, 478-578 core-hours; D2 sizing and block D's fit record
1423d1e5 decision 56 tasks 1-2 (Report D56-1): figure 1 with block D (figures_261015c) and the Stage 2 / D2 sizing (experiments_gen3_d56_261015)
595f8d0b plan-author decision 55 addendum 3 and decision 56 recorded (261012 sec. 4.7.133): consolidation closed; parallel project rule; figure 1 with block D, Stage 2 and D2 sizing
000f190e pid: sealed folder, M0 brief and blueprint v2.0
```

- **HEAD.** HEAD was `c37bad5e` when M0 started. M0's first commit is `000f190e` (the brief and the blueprint). Every other commit listed is the HardDisks session's: plan notes, the stage 1 figure 1c and the d56 sizing data. None of them touches an engine source.
- **Line numbers.** `git diff --stat 000f190e HEAD -- hspist3/00ALLINONE.c hspist3/edmd_core hspist3/experiment_validation.c hspist3/experiment_validation.h hspist3/Makefile` is empty, and so is the same diff against `a71f394d` and against the working tree. The engine sources are therefore identical at HEAD, at `000f190e`, at `a71f394d` and in the working tree, and every engine line number in this report holds for all four. The HardDisks commits do add files elsewhere under `hspist3/` (figures, d56 data, scripts).
- **Binary.** `hspist3/experiments_gen3_d52_261014/bin_a71f394/00ALLINONE` was built from those same sources (A1).

`git status --short` at the start of M0, before the first commit:

- It is 6,938 lines (645 KB): 22 tracked entries and 6,916 untracked (`??`). The SHA-256 of the full output is `1eb4435934f16b6e2435df788707a4a07238c0cdaa94e8e981c47df9b11c0779`.
- **Deviation from the brief, accepted by Chris on 2026-10-15:** the 22 tracked lines are given verbatim, the untracked entries are summarised by folder, and the SHA-256 above covers the full output.
- Of the untracked entries, 4,922 are inside `hspist3/experiments_gen3_d50_261014/`. The only untracked entry of M0 itself was `?? PHYSICALINFORMATIONDYNAMICS/`.

Tracked entries, verbatim:

```
 M .vscode/settings.json
 M 0000_PLAN_OVERALL/paper2_energytransfer/writeup/paper2_draft.pdf
 M 0000_PLAN_OVERALL/paper2_energytransfer/writeup/paper2_draft.tex
 m EDMD_4VALIDATE
 M LAMMPS/log.lammps
 D "Rom\303\241n et al. - 2002 - The speed of sound in a hard disk gas A computer simulation.pdf"
 D "Smallenburg - 2022 - Efficient event-driven simulations of hard spheres.pdf"
 M "hspist3/FINAL speed_of_sound_on_packing_fracture.pdf"
 M hspist3/Power_Freq_Spectrum_with_boxlength_None_and_wallfactor_None.pdf
 M hspist3/Power_Freq_Spectrum_with_boxlength_None_and_wallfactor_None.png
 M hspist3/divider_x_displacement_with_boxlength_None_and_wallfactor_None.pdf
 M hspist3/divider_x_displacement_with_boxlength_None_and_wallfactor_None.png
 M hspist3/energy_log.csv
 M hspist3/experiments_energy_transfer/00_COMMAND.md
 M hspist3/experiments_energy_transfer/energy_transfer_trace.csv
 m hspist3/kissfft
 M hspist3/wall_x_FFT.py
 D log.lammps
 D traj_exact.lammpstrj
 D traj_wall_line.lammpstrj
 D wall_exact.dat
 D wall_line.dat
```

Untracked entries per folder (count, folder; `hspist3/` subfolders listed separately; "(repo root files)" = files directly in the repository root; `"hspist3/` and `"0000_PLAN_OVERALL/` = names git prints quoted because of spaces or non-ASCII characters; the counts sum to 6,916):

```
4922 hspist3/experiments_gen3_d50_261014/
 228 hspist3/experiments_gen3_d52_261014/
 226 hspist3/experiments_gen3_d18_261010/
 199 hspist3/experiments_gen3_p1_261009/
 172 hspist3/experiments_gen3_d55_261015/
 145 hspist3/experiments_gen3_p2_261009/
 128 hspist3/experiments_gen3_preview_koa_261011/
 106 hspist3/experiments_gen3_p3a_261009/
  92 hspist3/experiments_gen3_h_261009/
  84 hspist3/experiments_gen3_p3_261009/
  80 hspist3/experiments_gen3_sym_261009/
  72 hspist3/experiments_gen3_h3b_261009/
  72 hspist3/
  50 hspist3/experiments_gen3_m3_261009/
  48 hspist3/experiments_gen3_p4_261010/
  41 hspist3/experiments_gen3_gate_koa_261009/
  40 hspist3/experiments_gen3_p6_261010/
  34 hspist3/experiments_gen3_p5_261010/
  28 0000_PLAN_OVERALL/
  23 hspist3/experiments_gen3_m4b_261009/
  19 hspist3/experiments_gen3_d19_261010/
  17 hspist3/experiments_gen3_d34a_261012/
  15 hspist3/experiments_gen3_d34s_261012/
  10 hspist3/experiments_gen3_psi6c_261010/
  10 (repo root files)
   9 hspist3/01_improvements_bugsfxed_dev/
   7 hspist3/experiments_gen3_testR_koa_261011/
   4 LAMMPS/
   3 "hspist3/
   2 hspist3/experiments_gen3_d48_261013/
   2 hspist3/experiments_gen3_c1_261011/
   1 particle-life-app-main/
   1 neuraipil/
   1 hspist3/manifests/
   1 hspist3/hspist3/
   1 hspist3/experiments_szilard/
   1 hspist3/experiments_simple_box/
   1 hspist3/experiments_resched_gate_261005/
   1 hspist3/experiments_resched_gate2_261007/
   1 hspist3/experiments_gen3_stage1_koa_261012/
   1 hspist3/experiments_gen3_s1_261012/
   1 hspist3/experiments_gen3_s10_261012/
   1 hspist3/experiments_gen3_m5_261009/
   1 hspist3/experiments_gen3_m4_261009/
   1 hspist3/experiments_gen3_gate_261009/
   1 hspist3/experiments_gen3_d16_261010/
   1 hspist3/experiments_energy_transfer/
   1 hspist3/SINGLE_analysis_results/
   1 hspist3/PRTICLELIFE_EXAMPLES/
   1 hspist3/BUG_SEARCH_FREQUENCIES/
   1 hspist3/02_PRESENTATION_HARD_DISC_SIM/
   1 experiments_energy_transfer/
   1 boxtest_2_incl_MB_heatbath/
   1 boxtest_2/
   1 ZZZ_PAPER/
   1 PHYSICALINFORMATIONDYNAMICS/
   1 IB_Package/
   1 0_5MHT/
   1 "0000_PLAN_OVERALL/
```

## How this report was made

- **Mapping and verification.** Five read-only agents each mapped one area: build, provenance and regression (A1, A2, E2, F1, F2); flags (A3); scene and neighbour search (B1–B3, C5); dynamics core (C1–C4, C6); renderer and outputs (D1–D3, E1).
  - A second agent per area then re-opened every cited line, re-ran each "not found" search with other patterns, and corrected the text. The verifiers checked 1,384 claims, found 108 wrong or unsupported, and made 150 corrections. Per area, checked / wrong: build 214/20, flags 430/23, scene 245/16, core 245/27, view 250/22.
  - Items A–F are the verified text, shortened only where a provenance note repeated the paragraph above.
- **Main session.** It ran item G itself and re-checked against the source the claims that H rests on: gen3's wall representation, the cell grid and origin shift, the location of the schedule audit, the checkpoint header, the `f72d764`..HEAD diff of `edmd_gen3.c`, and the blueprint line numbers. Items G and H and the explicit answers in B2, C5 and C6 are the main session's.
- **Citation conventions.**
  - `path:line` is relative to the repository root.
  - Inside a block that names a file, a bare number is a line of that file.
  - In A3 a bare number is a line of `hspist3/00ALLINONE.c`.
  - `00ALLINONE.c:n` means `hspist3/00ALLINONE.c:n`, and `edmd_gen3.c`, `edmd.c`, `edmd.h`, `edmd_gen3.h` and `tests/` are under `hspist3/edmd_core/`.
- **Binary runs.** The agents' runs were wrapped in a 60 s cap (`perl -e 'alarm 60; exec @ARGV'`), except the build agent's two `--version` calls, which exit at once.
  - The flags agent ran `bin_a71f394` with `--version`, `--help`, `--seeding=grid`, `--timescale=0` and `--edmd-acc=1`. All exit in the parser.
  - The view agent ran `bin_a71f394` with `--version`, `--help`, a 2 s headless gen3 scene and a 1 s window start with `SDL_VIDEODRIVER=dummy`. The window start returned 1 for lack of the font.
  - These ran with cwd `PHYSICALINFORMATIONDYNAMICS/_m0_scratch/<agent>_run/`.
  - The build agent ran `--version` of `bin_a71f394` and `bin_00f4dd0e-dirty` with cwd `/`; `--version` exits before the first file write (A1).
  - The main session ran the gen3 M1 harness it had built (G).
- **Scratch files.** The `_m0_scratch/` files cited in D3 and E1 (e.g. `view_run/energy_log.csv`) were read before that folder was deleted.
- **Not done or not determined** (also marked in the items):
  - no rule-4, replay, T1 or harness batch was rerun; their timings are quoted from records (F2);
  - whether gen2 trajectories depend on the `--fixed-dt` stride or differ between the scene loop and the experiment loop (C1, C3);
  - the unity-build route to the renderer was not compiled (D2);
  - whether a box-only run of a forked gen3 reproduces the unmodified gen3 event hash is deduced from the code, not tested (C5, H);
  - the Part A engine pin (A5 is still to be written, F1);
  - whether gen3's energy-transfer path is the Paper 2 validation branch (H);
  - `render_energy_hud`, `render_live_controls`, the Szilard outputs and `hspist3/glgraph.[ch]` were not itemised;
  - `merge_s1_koa.py`, `make_s1_manifests.py` and `reduce_s1_koa.py` were read by grep only.

## Summary

Each point is detailed, with `path:line`, in the item named in brackets.

- **Engines are libraries; the driver is not.** gen2 (`edmd.c`) and gen3 (`edmd_gen3.c`) are separate translation units with handle APIs (`EDMD*`, `EDMD3*`) and reference no symbol of `00ALLINONE.c`. They compile and link without it: the main session's compile check found only `_main` undefined, and both standalone harnesses built; the M1 harness ran [A2, C4, G].
- **gen3 is purely event-driven.** It has one binary heap keyed (t, type, a, b), lazy invalidation by collision counters and body epochs, per-disk clocks, a cell grid with crossing events, and no time step inside the engine [C1]. In the EDMD modes `--fixed-dt` is the driver's stop interval, not an integration step. A piston protocol that varies in time is sampled once per stop, which makes the coupled driver+engine a hybrid in that one respect [C1]. `--seed-drift-order` only orders the initial velocity normalisation [C1].
- **Geometry is the axis-aligned box only, as the brief expected.** It consists of four outer walls, up to 32 vertical divider slabs and two vertical piston faces. Every wall predictor tests x alone or y alone. There are no curved, sloped or polyline walls, no obstacles and no fixed disks in any engine or in the time-stepped code [B2]. In gen3 the walls are a side index s = 0..3 hard-coded in switch statements, not data [B2].
- **One disk radius, unit mass.** The cell width and capacity derive from that single radius [B1, C5].
- **The renderer cannot be linked or adapted without editing `00ALLINONE.c`.** All drawing lives in the 23,715-line unit that holds `main`. Its primitives are `static`, the draw functions take no state and read file-scope globals, and the frame has no hook. A second disk population, a polyline/arc set and a field overlay are each "only by edit" [D1, D2].
- **The histograms read globals** and are not reusable unmodified [D3]. Keys 1–5 and I are the gen3 return-map panel, not the histograms [D1].
- **Outputs.**
  - A headless scene writes summaries only.
  - Full-state forms: the gen3 initial-state text dump; binary same-build checkpoints; positions-only snapshots (`G3SNAP1`); an environment-gated stop snapshot; the harness replay format [E1].
  - `00_COMMAND.md` is written into the cwd on every run past the parser. It holds argv only, with no build hash and no seed unless the seed is in argv [E1, E2].
- **The frozen set is not what the brief states.** The Stage 1 SHA-256 manifests list 21 files (task list, worker, build record, reducers, validation scripts) and the task-list hash. No engine source is in any manifest. The engine is pinned by commit `f72d764` and by the SHA-256 of the KOA binaries. HEAD's engine differs from `f72d764` in `00ALLINONE.c` (1,978 changed lines) and in `edmd_gen3.c`/`.h` by 36 inserted lines (25 + 11), the display-only collision hook; `edmd.c`, `edmd.h`, the Makefile and `experiment_validation.*` are unchanged. Part A has no manifest yet [F1].
- **Regression tooling exists but is spread out.** Pieces: the harnesses' A/B/C hash equality; rule 4 (gen2 byte identity on ctrl_min, ctrl_leg); the replay through the driver (64 of 64 MATCH); T1 (Part B unchanged); the GUI hash. There is no single runner and no `make test`. Several committed scripts still point at the locked `HardDisks_gen3` worktree [F2].
- **Recommendation (H; the decision is Chris's).**
  - Demo 0: a new small C/SDL2 executable in this folder (option b), with no engine.
  - Demo 1: a fork of `edmd_gen3.c` inside this folder with a segment/arc event type (option d). Its box-only runs are regression-tested against the unmodified gen3 objects linked into a test, which is the cross-check part of (a).
  - The renderer is not reused as code; its techniques are copied as reference.

## A. Build and entry point

### A1

**Makefile targets** (hspist3/Makefile, 128 lines)

| target | lines | CFLAGS | stamped BUILD_TARGET / BUILD_CFLAGS |
|---|---|---|---|
| `all` -> `debug` (plain `make`) | 66, 89-92 | `-g -fsanitize=address -DDEBUG $(INCLUDE)` (30) | `debug` / `-g -fsanitize=address -DDEBUG` (89-90) |
| `release` | 94-97 | `-O3 -march=native -ffp-contract=off $(INCLUDE)` (40) | `release` / `-O3 -march=native -ffp-contract=off` (94-95) |
| `koa` | 99-102 | `-O2 -march=x86-64-v2 -mtune=generic -ffp-contract=off $(INCLUDE)` (80) | `koa` / the same without `$(INCLUDE)` (99-100) |
| `release-ld` | 108-111 | `$(CFLAGS_RELEASE) -DEDMD3_LONG_DOUBLE` (110) | `release-ld` / release + `-DEDMD3_LONG_DOUBLE` (108-109) |
| `koa-ld` | 113-116 | `$(CFLAGS_KOA) -DEDMD3_LONG_DOUBLE` (115) | `koa-ld` / koa + `-DEDMD3_LONG_DOUBLE` (113-114) |
| `clean` / `tar` | 123-124 / 127-128 | `rm -f $(TG)` / `tar czvf $(TG).tgz ...` | - |

- No `-Wall` and no `-std=` in any CFLAGS (30, 40, 80). No test target: `.PHONY: all debug release koa release-ld koa-ld clean tar` (63).
- Compiler: `CC ?= cc` (16). This Mac: `cc --version` -> `Apple clang version 17.0.0 (clang-1700.4.4.1)`, `Target: arm64-apple-darwin25.6.0` (`/usr/bin/cc`). KOA: `compiler gcc (GCC) 14.3.0` (hspist3/cluster/gen3_stage1_261012/build_record_f72d764.txt:9).
- One cc command compiles and links everything (119-120): `$(CC) $(LDFLAGS) -o $@ $(TG).c $(KISSFFT_SRCS) $(EDMD_SRCS) $(VALIDATION_SRCS) $(CFLAGS) $(PROVENANCE) $(LIBS)` (120), with `TG ?= 00ALLINONE` (4), `KISSFFT_SRCS := kissfft/kiss_fft.c kissfft/kiss_fftr.c` (5), `EDMD_SRCS ?= edmd_core/edmd.c edmd_core/edmd_accelerated.c edmd_core/edmd_gen3.c` (12), `VALIDATION_SRCS := experiment_validation.c` (13). Output: `hspist3/00ALLINONE`.
- Prerequisites (119) are the .c files only; no header is a prerequisite. KOA builds force a rebuild with `make -B koa-ld` / `make -B koa` (hspist3/cluster/gen3_s1_gate_261012/build_s1_koa.sh:27, :32; also cluster/build_koa.sh:32).
- Other build files in hspist3/ do not build 00ALLINONE: `CMakeLists.txt` (tracked, 29 lines) globs `src/*.c` (:15, i.e. hspist3/src/hspist3.c) and declares `add_executable(hspist3 ...)` twice (:5, :18); `Makefile.gui` (53 lines) builds `gui_config` / `gui_advanced` from gui_config.c / gui_advanced.c (:7-8, :26-30).

**SDL2 and link line**

- SDL flags from `pkg-config --cflags sdl2 SDL2_ttf` (22) and `--libs` (23); fallbacks `-I/opt/homebrew/include -I/opt/homebrew/include/SDL2 -D_THREAD_SAFE` (26) and `-lSDL2 -lSDL2_ttf` (58).
- `LDFLAGS += -L/opt/homebrew/lib -L/opt/homebrew/opt/mesa/lib` (43). GL: on Darwin `-framework OpenGL -framework CoreGraphics` (47); on Linux `-lGL` (49) and `SYS_LIBS := -lm` (53). `LIBS_BASE := -lGLEW $(GL_LIBS) $(SYS_LIBS)` (57); `LIBS := $(if $(PKG_SDL_LIBS),$(PKG_SDL_LIBS),-lSDL2 -lSDL2_ttf) $(LIBS_BASE)` (58).
- pkg-config on this Mac (`--modversion`: SDL2 2.30.12, SDL2_ttf 2.24.0):
  - cflags: `-I/opt/homebrew/include/SDL2 -I/opt/homebrew/include -D_THREAD_SAFE -I/opt/homebrew/Cellar/harfbuzz/12.2.0_1/include/harfbuzz -I/opt/homebrew/Cellar/glib/2.86.3/include/glib-2.0 -I/opt/homebrew/Cellar/glib/2.86.3/lib/glib-2.0/include -I/opt/homebrew/opt/gettext/include -I/opt/homebrew/Cellar/pcre2/10.47_1/include -I/opt/homebrew/Cellar/graphite2/1.3.14/include -I/opt/homebrew/opt/freetype/include/freetype2 -I/opt/homebrew/opt/libpng/include/libpng16`
  - libs: `-L/opt/homebrew/Cellar/sdl2_ttf/2.24.0/lib -lSDL2_ttf -L/opt/homebrew/lib -lSDL2`
  - The same strings appear in the recorded build line experiments_gen3_d52_261014/bin_a71f394/build_record.txt:30.
- A recorded expansion of Makefile:120 on this Mac (a KOA rehearsal with modified CFLAGS, not a production target): hspist3/experiments_gen3_k2_rehearsal_261010/run6_fd8d8a2_retry/logs/make_koa_home_harddisks_gen3.txt:1 = `cc -L/opt/homebrew/lib -L/opt/homebrew/opt/mesa/lib -o 00ALLINONE 00ALLINONE.c kissfft/kiss_fft.c kissfft/kiss_fftr.c edmd_core/edmd.c edmd_core/edmd_accelerated.c edmd_core/edmd_gen3.c experiment_validation.c -O2 -ffp-contract=off -I/opt/homebrew/include/SDL2 ... -DBUILD_GIT='"a643bc7"' -DBUILD_TARGET='"koa"' -DBUILD_CFLAGS='"-O2 -ffp-contract=off [mac rehearsal]"' -L/opt/homebrew/Cellar/sdl2_ttf/2.24.0/lib -lSDL2_ttf -L/opt/homebrew/lib -lSDL2 -lGLEW -framework OpenGL -framework CoreGraphics`. `make -n` was not run.
- GLEW and OpenGL are linked but no source calls them: `glew|GLEW|glBegin|SDL_GL_|<GL/|OpenGL/` and `\bgl[A-Z][a-zA-Z]+\(` over 00ALLINONE.c, edmd_core/*.[ch], experiment_validation.[ch]: 0 hits. The window has no OpenGL flag (00ALLINONE.c:5786 `SDL_WINDOW_SHOWN | SDL_WINDOW_RESIZABLE | SDL_WINDOW_ALLOW_HIGHDPI`); drawing goes through `_renderer = SDL_CreateRenderer(_window, -1, SDL_RENDERER_ACCELERATED);` (5792), software fallback at 5796. Per comment at cluster/build_koa.sh:4-6 the GUI libraries are "needed only to LINK" on KOA.

**What stamps the hash**

- Makefile:86: `BUILD_GIT := $(shell git -C $(dir $(abspath $(lastword $(MAKEFILE_LIST)))) rev-parse --short HEAD 2>/dev/null || echo unknown)$(shell git -C $(dir $(abspath $(lastword $(MAKEFILE_LIST)))) diff --quiet HEAD 2>/dev/null || echo -dirty)`
- Makefile:87: `PROVENANCE       = -DBUILD_GIT='"$(BUILD_GIT)"' -DBUILD_TARGET='"$(BUILD_TARGET)"' -DBUILD_CFLAGS='"$(strip $(BUILD_CF))"'`
- The -dirty test has no pathspec: it is repo-wide, counts gitlinks, and ignores untracked files (it is `git diff`). On this tree:
  - `git diff --quiet HEAD` rc=1;
  - `git diff --quiet HEAD -- . ':!hspist3/kissfft'` rc=1: 21 other modified or deleted tracked paths (among them hspist3/energy_log.csv and the EDMD_4VALIDATE gitlink; `git diff --name-only HEAD | wc -l` = 22);
  - `git diff --quiet HEAD -- hspist3/kissfft` rc=1 (`Subproject commit febd4caeed32e33ad8b2e0bb5ea77542c40f18ec-dirty`; the nested clone shows ` M kiss_fft_log.h` and `?? .DS_Store`);
  - the engine paths: rc=0.
  - Inference from line 86 and these return codes (nothing was built): `make` in this tree would stamp `<short HEAD>-dirty`, i.e. `595f8d0b-dirty` when this was checked (`--short` gives 8 characters now).
- KOA builds must be clean: build_s1_koa.sh:21 stops on any modified tracked file (`git status --porcelain --untracked-files=no`), :29 and :34 require the `--version` line to contain `git $head  target koa-ld` / `git $head  target koa`. Recorded: build_record_f72d764.txt:3 `version 00ALLINONE  git f72d764  target koa`, :4 `version_ld 00ALLINONE  git f72d764  target koa-ld`.
- Defaults without -D: 00ALLINONE.c:659-667 (`#define BUILD_GIT "unknown"`, the same for BUILD_TARGET and BUILD_CFLAGS).

| where the stamp appears | 00ALLINONE.c | verbatim / fields |
|---|---|---|
| `--version`, then `exit(0)` | 4094-4096 | `printf("00ALLINONE  git %s  target %s\n  CFLAGS: %s\n", BUILD_GIT, BUILD_TARGET, BUILD_CFLAGS);` (4095) |
| gen3 run header `[EDMD3-RUN]` | 17409-17410 | `engine=gen3 build=\"%s\" target=%s ...` |
| gen3 run record `[EDMD3-HEALTH]` | 17437-17451 | `hash=%016llx t_end=%.17g engine=gen3 build=\"%s\" target=%s` (17442), arguments at 17449-17450 |
| gen3 checkpoint identity | 17592-17595 | `c->build`, `c->target` (17592), `c->seed` (17595) |
| energy-transfer summary CSV header | 19641-19642 | `...,command,seed,build_git,build_target,build_cflags,box_width_sigma` |
| its row | 19831 | `fprintf(sf, ",%u,%s,%s,\"%s\",%.6f\n", cli_seed, BUILD_GIT, BUILD_TARGET, BUILD_CFLAGS,` |
| psi6-snapshot reader header | 23463-23464 | `... build %s` |

- `--version` is handled in parse_cli_options (defined at 4071; argument loop 4091; exit at 4096). main calls it at 23494, before init_output_dirs_from_argv0 (23497; it only sets path strings, 3395-3422) and before `write_command_md_to_dir(".", argc, argv);` (23500). I ran no binary. Read-only evidence of the stamps:
  - `strings experiments_gen3_d52_261014/bin_a71f394/00ALLINONE` contains `a71f394`, `mac-O3-e0pre` and `-O3 -ffp-contract=off`; `strings .../bin_00f4dd0e-dirty/00ALLINONE` contains `00f4dd0e-dirty`.
  - Recorded at run time (tracked files): experiments_gen3_d55_261015/tests/g5_T1/B_N0100_eta0.050000_m_50_r0_old/stdout.log:39, the `[EDMD3-RUN] L0=78.5 M=50 run=0 seed=974638985:` line with `engine=gen3 build="a71f394" target=mac-O3-e0pre`; experiments_gen3_d55_261015/e0/default/ctrl_min/plain/version.txt:1 `00ALLINONE  git 0544b7e8-dirty  target mac-O3-e0pre`.

**How `bin_<hash>/` folders are made: by scripts, not by the Makefile**

| script | hash | sources | flags | output |
|---|---|---|---|---|
| hspist3/experiments_gen3_m3_261009/build_clean.sh | `H=$(git -C "$WT" rev-parse --short=7 "$REV")` (10); no dirty test | `git -C "$WT" archive "$REV" hspist3` into `$OUT/src` (11-12); kissfft's 6 files copied from the worktree (13-14); `WT=/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks_gen3` (9) | `-O3 -ffp-contract=off $SDLC "-DBUILD_GIT=\"$H\"" '-DBUILD_TARGET="mac-O3-e0pre"' '-DBUILD_CFLAGS="-O3 -ffp-contract=off"' $SDLL -L/opt/homebrew/lib -lGLEW -framework OpenGL -framework CoreGraphics` (33) | `<out>/00ALLINONE`, the harness binaries (36-41), build_record.txt (15-22, 26, 30-31, 35, 42), *.warnings.txt (29), 00ALLINONE.Wall.txt (34), and the source copy `<out>/src/`; it cds into `$OUT/src/hspist3` (23), so OUT must be absolute; usage at 6 |
| experiments_gen3_d18_261010/d18_build.sh, d19_261010/d19_build.sh, d34s_261012/s_build.sh, d34a_261012/w_build.sh | build_clean.sh's | build_clean.sh into a scratch folder, then `cp -p` of chosen files into a bin folder the caller names | build_clean.sh's | d19_build.sh:6 usage `<commit> <scratch build dir, absolute> <bin dir, absolute>`, :14, :15-18; s_build.sh:15-19; d18_build.sh:14, :25-28 |
| experiments_gen3_d52_261014/scripts/build_d52.sh | `GIT=$(git -C "$G" rev-parse --short HEAD)$(git -C "$G" diff --quiet HEAD -- 00ALLINONE.c edmd_core experiment_validation.c experiment_validation.h kissfft 2>/dev/null \|\| echo -dirty)` (10); optional tag (11-12) | the working tree `G=/Users/.../Repo/HardDisks_gen3/hspist3` (7) | the same flags, SDL flags hard-coded (17-23) | `bin_$GIT[-tag]` (12); refuses an existing folder (14); record (28-44) with the sources' SHA-256 (38) and `git diff HEAD -- 00ALLINONE.c edmd_core` (43) |
| experiments_gen3_d50_261014/scripts/build_d50.sh | the same line 10 | the same G (7) | the same | differs from build_d52.sh (by `diff`) in comments, `OUT="$G/experiments_gen3_d50_261014/bin_$GIT"` with no tag (11), the fixed note ` plus the uncommitted edit of 00ALLINONE.c` (12), and a recorded diff of 00ALLINONE.c only (42) |
| experiments_gen3_d48_261013/scripts/build_d48.sh | repo-wide `git -C "$G" diff --quiet HEAD` (8) | HardDisks_gen3 (6) | as above (12-19) | `bin_$GIT` (9) |

- **bin_a71f394**: experiments_gen3_d52_261014/README.md:41 "`bin_a71f394`: build_clean.sh of a71f394d". build_record.txt:1 `# clean build of engine-gen3 a71f394dac037200df2038d8c6a59d13c192df14 (a71f394), 2026-10-14 22:19:14 HST, Apple M3 Max, Apple clang version 17.0.0 (clang-1700.4.4.1)`; built into a scratchpad (`cc -o /private/tmp/.../scratchpad/d52/build_a71f394/00ALLINONE ...`, :30). Its SHA-256 `5f9f898e...a174ea` (:48) equals `shasum -a 256` of the binary today. The copy from that scratch folder into a bin folder: not found in any script (`grep -rln bin_a71f394` over *.sh and *.py in hspist3/ and 0000_PLAN_OVERALL/: 0 files). The folder holds only 00ALLINONE, 00ALLINONE.Wall.txt, 00ALLINONE.warnings.txt and build_record.txt; the harness binaries hashed at :49-53 are not there.
- **bin_00f4dd0e-dirty**: build_d52.sh. Record :1 `... 00f4dd0e-dirty (00f4dd0ec5955da9c5d2c3a71248e67b2e310f12 plus the uncommitted edits of 00ALLINONE.c and edmd_core) ...`; command at :4-5 (6 s, output path `.../Repo/HardDisks_gen3/hspist3/experiments_gen3_d52_261014/bin_00f4dd0e-dirty/00ALLINONE`); SHA-256 at :10-22; the diff from :24.
- Both binaries were built in the HardDisks_gen3 worktree and moved into this tree with its untracked files by 0000_PLAN_OVERALL/consolidation_d51_261014/d51_move_untracked.py (purpose at :2-9). Its manifest lies outside the repository: Repo/HardDisks_d51_hold/MANIFEST_d51_move.tsv:5632-5633 list `hspist3/experiments_gen3_d52_261014/bin_00f4dd0e-dirty/00ALLINONE main 692952` and `.../bin_a71f394/00ALLINONE main 692952` (tab-separated).
- bin_0544b7e8-dirty (experiments_gen3_d50_261014/) was built in place in this tree at the consolidation by a scratchpad copy of build_d50.sh (0000_PLAN_OVERALL/consolidation_d51_261014/d55_consolidation_log_261015.txt:51-52). Its record :1 says ` plus the uncommitted edit of 00ALLINONE.c`; that is build_d50.sh's fixed text (:12). Per experiments_gen3_d55_261015/README.md:7-8, "-dirty" "is only the kissfft nested clone with the project's local kiss_fft_log.h change".
- Git tracks only the .txt files of these folders. The binaries are untracked, not ignored (`?? hspist3/experiments_gen3_d52_261014/bin_a71f394/00ALLINONE`; `git check-ignore` rc=1). The .gitignore patterns that name 00ALLINONE (.gitignore:64, :65, :84) cover only files directly in hspist3/.
- Mac bin builds differ from `make release`: no `-march=native` (build_clean.sh:33, build_d52.sh:23 against Makefile:40), and `mac-O3-e0pre` is not a Makefile target.
- Today none of the committed build scripts can run as written. build_clean.sh copies kissfft from `$WT/hspist3/kissfft` (:14), and build_d48/d50/d52.sh compile in `HardDisks_gen3/hspist3` (G at :6, :7, :7). That worktree's kissfft folder is empty (`ls -la` shows only `.` and `..`, mtime 2026-10-15 11:10); per comment at d51_move_untracked.py:4-9 the mover takes the nested clone out and recreates the empty gitlink folder. The consolidation used scratchpad copies retargeted to Repo/HardDisks (log :28-31; see F2).

### A2

**Translation units: separate. No source of the 00ALLINONE build includes a .c file.**

- 00ALLINONE.c:17-22, the only quoted includes in the file (lines 1-16 are system and SDL headers):
```
#include "kissfft/kiss_fft.h"
#include "kissfft/kiss_fftr.h"
#include "edmd_core/edmd.h"
#include "edmd_core/edmd_accelerated.h"
#include "edmd_core/edmd_gen3.h"   /* ##CHRIS 2026-10-09 (M3, 261012 sec. 4.7.14): --engine=gen3 */
#include "experiment_validation.h"
```
- `grep -n '#include "[^"]*\.c"'` over 00ALLINONE.c, edmd_core/*.c and experiment_validation.c: 0 hits.
- The seven .c files are separate units on one cc command line (Makefile:120; recorded at experiments_gen3_d52_261014/bin_a71f394/build_record.txt:30).
- The accelerated engine links beside gen2 because it renames its API: edmd_accelerated.c:4-33 (`#define edmd_create edmd_acc_create` ...), then `#include "edmd.h"` (:34).
- Includes inside edmd_core: edmd_gen3.c:66 `#include "edmd_gen3.h"`, which includes `edmd.h` (edmd_gen3.h:28); edmd.c:1 `#include "edmd.h"`.
- Only tests include .c files: tests/gen3_body_rule_test.c:19 `#include "../edmd_gen3.c"`; tests/gen3_band_edge_ties.c:32 `#include "../edmd_gen3.c"` (after defining `EDMD3_AUDIT_HOOK`, :31); tests/gen3_checkpoint_test.c:21-23 `#define main gen3_m2_harness_main` / `#include "gen3_m2_harness.c"` / `#undef main`.
- Builds without 00ALLINONE.c and SDL that are on record, each `warnings: 0` (bin_a71f394/build_record.txt):
  - :33 gen3_m1: `cc -O3 -ffp-contract=off -Wall -Wextra -o ... edmd_core/tests/gen3_m1_harness.c edmd_core/edmd_gen3.c edmd_core/edmd.c -lm`;
  - :39 gen3_band_edge_test: `... edmd_core/tests/gen3_band_edge_test.c edmd_core/edmd_gen3.c -lm` (no edmd.c);
  - :43 edmd_gen3_alone: `cc -std=c11 -O3 -ffp-contract=off -Wall -Wextra -c edmd_core/edmd_gen3.c -o .../edmd_gen3.o`.
  - On KOA: `gcc $HF -o bin/gen3_m1 edmd_core/tests/gen3_m1_harness.c edmd_core/edmd_gen3.c edmd_core/edmd.c -lm` with `HF="-std=gnu11 -O2 -march=x86-64-v2 -mtune=generic -ffp-contract=off"` (cluster/gen3_s1_gate_261012/build_s1_koa.sh:39-40).
  - A standalone build of edmd_accelerated.c: not found (grep for `edmd_accelerated.c` in *.sh, *.txt and *.py under hspist3/, src/ copies excluded: only full-driver build lines).

**Every file under hspist3/edmd_core/** (all 15 tracked; `wc -l`, 7392 in total)

| file | lines | purpose (from its header or code) |
|---|---|---|
| EDMD_ACCELERATED_CORE.md | 139 | Doc of the `edmd_acc_*` backend (:1-8). :8 "Both are compiled into `00ALLINONE` by default, and you choose at runtime via a CLI flag." is stale: `--edmd-acc` bare, empty, or 1/true/yes/on calls `refuse_edmd_acc_backend()` (00ALLINONE.c:5266-5276), which prints STOP and `exit(2)` unless built with `-DHD_ALLOW_ACC_BACKEND` (1360-1367). |
| edmd.c | 2086 | gen2 event-driven engine; no header comment. Code: heap (318-338), grid (comment `grid (uniform spatial hashing)` 353; 355-371), wall, divider and piston collision times (518-954), `resolve_ab` (1051), `resolve_wall` (1286), `EDMD* edmd_create(const EDMD_Params* prm_in){` (1441), `double edmd_advance_to(EDMD* S, double t_target){` (1650). |
| edmd.h | 240 | gen2 API. :4-5 "Minimal, C99, event-driven molecular dynamics for 2D hard disks. Static rectangular box [0..boxW] x [0..boxH]. No pistons here (yet)." The piston part is stale: edmd.c:934 `static int collide_time_piston_L(...)`. |
| edmd_accelerated.c | 1431 | "Accelerated EDMD core: provides the same API as edmd.c, but exported under an `edmd_acc_*` prefix" (:1-3). |
| edmd_accelerated.h | 53 | The `edmd_acc_*` declarations (:4-6). |
| edmd_gen3.c | 1548 | "the generation-3 EDMD core" (:1-2; design comment to :64): cells, crossing events, one heap with lazy invalidation, floating time origin, bodies as bands, ledgers. |
| edmd_gen3.h | 208 | gen3 API: "generation-3 event-driven MD for 2D hard disks, a constant-work-per-event engine" (:5; header comment :4-23). |
| tests/gen3_band_edge_cells.h | 83 | The band-edge stress cells shared by the two band-edge tests (:1-20). |
| tests/gen3_band_edge_test.c | 68 | Band-edge cells with schedule, body and contact audits (:1-12). |
| tests/gen3_band_edge_ties.c | 141 | White box: the numbers behind each schedule-audit finding, classes G/K/V/C/U (:1-23). |
| tests/gen3_body_rule_test.c | 242 | White box: the spring divider's `body_rule` against an independent scan, plus constructed cases (:1-18). |
| tests/gen3_checkpoint_test.c | 188 | Checkpoint and restart at an event boundary, compared bit for bit (:1-20). |
| tests/gen3_m1_harness.c | 374 | M1 harness: gen3 runs A/B/C plus gen2 run D; audit, speed, diverge and dump modes (:1-18, :281). |
| tests/gen3_m2_harness.c | 537 | M2 harness with dividers and pistons: A/B/C/D; audit, speed and dump modes (:1-36, :468). |
| tests/gen3_state_io.h | 54 | Writes a harness cell as a replay file for `00ALLINONE --engine-replay=FILE` (:1-11). |

- No Makefile or build script inside edmd_core/. Build lines for the tests are in the test headers (see F2) and also in experiments_gen3_m3_261009/build_clean.sh:36-41, experiments_gen3_d18_261010/d18_build.sh:20 (gen3_checkpoint_test), cluster/gen3_s1_gate_261012/build_s1_koa.sh:39-41 and the hand-written experiments_gen3_m1_261008/00_COMMAND.md. They differ: gen3_m1's header asks for `-std=c11` (gen3_m1_harness.c:15), build_clean.sh:36 builds it without `-std`, KOA with `-std=gnu11` (build_s1_koa.sh:39).

**Other sources the build uses**

- **kissfft** (gitlink):
  - `git ls-files -s hspist3/kissfft` = `160000 febd4caeed32e33ad8b2e0bb5ea77542c40f18ec 0` (path `kissfft`). No `.gitmodules` (`ls .gitmodules hspist3/.gitmodules`: no such file).
  - The nested clone is at febd4cae with a local 3-line change: kiss_fft_log.h:16-18 `#ifndef DEBUG` / `#define DEBUG 4` / `#endif` (`git -C hspist3/kissfft diff --stat`: 1 file, 3 insertions); plus an untracked `.DS_Store`.
  - Compiled: kiss_fft.c (420 lines), kiss_fftr.c (155) (Makefile:5). Headers: kiss_fft.h (160), kiss_fftr.h (54), _kiss_fft_guts.h (167), kiss_fft_log.h (39); kiss_fft_log.h is included only by _kiss_fft_guts.h:18, not by 00ALLINONE.c.
- **experiment_validation.c** (373 lines) and **.h** (94): the invariant checker; types `ExperimentValidationState` (.h:29), `ExperimentValidationFailure` (:44), `ExperimentValidator` (:54); .h:68 `/* Returns 1 while valid and 0 after the first invariant violation. */`, :69 `int experiment_validator_check(ExperimentValidator *validator,`. Last commit touching the .c: 4adeddeb `build: track experiment_validation.c/.h — provenance gap closed`.
- **Runtime file, not a build input:** `font = TTF_OpenFont("Roboto-Regular.ttf", 18);` (00ALLINONE.c:23520) opens the font from the cwd, only when not headless (23514); the program returns 1 if it is missing (23521-23523). The font is tracked at hspist3/Roboto-Regular.ttf.
- **System libraries:** SDL2, SDL2_ttf, GLEW, OpenGL and CoreGraphics on the Mac (A1); GLEW and OpenGL are linked but not called.

### A3

**Convention: a bare number is a line in hspist3/00ALLINONE.c.**

#### Parser location
- **Function.** `static void parse_cli_options(int argc, char **argv)` spans 4071–5694. It is forward-declared at 690. `main` (23489) calls it once, at 23494, after storing `g_main_argc/argv` (23491–23492).
- **Loop structure.** One loop, `for (int i = 1; i < argc; ++i)` (4091–5315):
  - First a standalone `if` for `--version` (4094).
  - Then one if/else-if chain from `--help` (4098) to the final `else` (5310–5314). That `else` prints "Unknown option '%s'.", calls `print_cli_usage` and exits with EXIT_FAILURE.
  - The chain has 172 `if` heads holding 192 `strcmp(arg, "--x") == 0` (exact) or `strncmp(arg, "--x", n) == 0` (prefix) tests. A 193rd strncmp, at 5182, sits inside the `--gen3-virial-blocks=` branch body and is not a dispatch test.
  - There are 188 distinct literals.
- **No other parser.**
  - There is no getopt and no option table: `grep -n 'getopt\|struct option'` finds nothing.
  - Outside the parser, `argv[` appears in 00ALLINONE.c only at 3218 and 3295 (provenance printing, `fprint_shell_escaped`) and at 23497 (`init_output_dirs_from_argv0(argv[0])`).
  - `grep argv` over hspist3/edmd_core/*.[ch] and hspist3/experiment_validation.[ch] finds nothing.
  - Under hspist3/edmd_core/tests/, three harnesses reference argv: gen3_body_rule_test.c, gen3_m1_harness.c and gen3_m2_harness.c. Their parsing was not examined.
- **Value helper.** `static const char *cli_option_value(const char *arg, int argc, char **argv, int *index)` (3097–3107).
  - It returns the text after `=` if that text is non-empty.
  - Otherwise it consumes and returns the next argv element. With no next element it prints "Option '%s' expects a value." and exits.
  - So `--x=v` and `--x v` both work for branches that use it. `--x=` with nothing after the `=` consumes the next argument.
- **`=`-only branches.** Branches whose literal ends in `=` read `arg + strlen(...)` and accept only the `=` form. These are every `--gen3-*=`, every `--gui-*=`, `--psi6-every=`, `--record-sigma-time=`, `--engine-replay=`, `--replay-*=`, `--validator-every=` and `--gen3-cell-px=` (5124–5258). A space-separated value for these falls through to "Unknown option", because no earlier prefix test matches them.
- **Post-loop derivation, same function (5317–5693), in order:**
  - substeps clamps (5318–5333);
  - `--eff-output` default (5339–5341);
  - preset `switch` (5343–5395);
  - particle-box fractions (5397–5419);
  - wall count (5421–5429);
  - wall positions turned into fractions (5431–5454);
  - `--calibrate-work` (5456–5506);
  - globals overwritten (5508–5536);
  - radius from `--eta` / `--particle-radius` / `--particle-diameter` (5538–5636);
  - `phi_max` auto-raise (5637–5650);
  - `cli_force_no_experiments` clears `enable_speed_of_sound_experiments` (5652–5654);
  - repeats clamp (5661–5663);
  - fractions zeroed unless custom (5665–5667);
  - particle-box fractions mapped again (5669–5693).
- **Usage printer.** `static void print_cli_usage(const char *exe_name)` (3907–4069). It is called from `--help` (4099), from an unknown preset (4454) and from an unknown option (5312).
  - It has 161 `\n` and no conditionals, so it prints 161 lines.
  - 54 of the 188 parser literals are not in its text, by `comm` of parser literals against usage tokens. One of the 54, `--resched-audit=all`, is documented as `--resched-audit[=all]` (4066). That leaves **53 undocumented spellings, or 51 if the three `--engine` spellings count once**:
    - `--version`, `--engine`, `--engine=gen2`, `--engine=gen3`, `--engine-replay`;
    - all 19 `--gen3-*` and all 7 `--gui-*`;
    - `--psi6-every`, `--record-sigma-time`, `--replay-read-every`, `--replay-substeps`, `--validator-every`;
    - `--dt`, `--cyclic-boundaries`, `--pl-dt`, `--szilard-no-pp`, `--szilard-view`;
    - 12 underscore / alternate aliases of the piston and velocity flags, among them `--piston-w-in-target`.
  - Every flag the usage text lists is parsed (the reverse `comm` is empty). Four of the listed flags are dead (see shadowing).
- **Writers of the same variables after parsing (not parsers).**
  - The GUI live-scene panel `live_scene_prepare_rebuild` (15130) rewrites `sim_mode` (15164–15170, 15215), `cli_experiment_preset` and `cli_seeding_mode` (15215–15228). It also rewrites `cli_num_walls`, `cli_requested_walls`, `cli_override_particles`, `cli_wall_positions`, `cli_wall_masses`, `cli_particle_radius_sigma` and `particle_scale_runtime` (15258–15297).
  - Experiment runners reassign CLI variables too: 10598–10618 (`szilard_prepare_interactive_config`), 20665–20669 (`run_simple_prediction_box_experiment`) and 21143–21150 (`run_szilard_engine_experiment`).
- **Non-flag inputs.**
  - Environment variables: HD_CONTACT_AUDIT (1580, 2502, 17320; edmd_core/edmd.c:1750), HD_OVERDUE_TRACE (1823, 1833), HD_KE_TRACE (2458, 2477), HD_KE_COLUMNS (2517), HD_STOP_SNAPSHOT (2526) and HD_PISTON_EVENTS (18929).
  - The font `Roboto-Regular.ttf` is opened from cwd (23520, 23527).
  - There is no configuration-file loader. The only `fopen(...,"r"/"rb")` calls are:
    - `ensure_csv_header_schema` (3437);
    - `pl_load_tensor_file` (8252);
    - `g3_slurp` (17600), called by `g3_ckpt_write` (17623–17626);
    - the `--gen3-restart` read in `g3_restart_apply` (17654);
    - `run_gen3_replay` (23345);
    - `run_gen3_psi6_snapshots` (23453).

#### Prefix-match shadowing and duplicate branches
I simulated first-match dispatch over the 172 chain heads, in source order, with each of the 188 literals tried bare, with `=1` and with a suffix.

**Unreachable branches.** Each is reachable from no spelling, because an earlier branch catches it:
- `--seeding` (4625): `strncmp(arg, "--seed", 6)` at 4343 catches it first. `--seeding=grid` therefore reaches the seed parser, which rejects it with `Invalid seed '%s'. Use an integer or 'time'.` (4353).
- `--timescale` (5116): `strncmp(arg, "--time", 6)` at 4581 catches it first and sets `cli_time_span`, not `time_scale_runtime`. The comment at 1501–1502 says the same.
- `--szilard-unstick-hits` (4240) and `--szilard-unstick-vx-frac` (4246): `--szilard-unstick` (strlen 17) at 4237 catches them first. A value such as `--szilard-unstick-hits=5` sets `cli_szilard_unstick = 1` (4239).
- 4917, 4921 and 4926 (`--particles-boxes`, `--particles-box-left`, `--particles-box-right`) duplicate the tests at 4477, 4499 and 4504. The live branch at 4477 also handles absolute counts (4482–4494); the dead one at 4917 does not.

**Silent value-form traps on switches** (the value lands in another branch):
- `--substeps-auto=1` reaches `--substeps` (prefix 10, 4573) and sets `cli_substeps_override=1`.
- `--szilard-no-pp=1` reaches `--szilard-n` (4262) and sets `cli_szilard_n=1`.
- `--experiments=1` reaches `--experiment` (4374), which reports an unknown preset and exits.

All other `switch=1` forms (for example `--kbt1=1`, `--show-simulation=1`, `--headless=1`, `--pp-toi=1`, `--gen3-exact-box=1`) match no branch and exit through "Unknown option".

The usage text still advertises the dead flags: `--seeding` at 3949, `--timescale` at 4017 and the two unstick flags at 4042–4043.

**Parsed in more than one place:**
- `--particles-boxes`, `--particles-box-left` and `--particles-box-right` (above);
- `--edmd-acc`: 5266 bare, 5269 with a value;
- `--engine`: 5212, 5214, 5249, and `cli_engine = 3` is also set by `--engine-replay=` (5203).

#### --mode
- Parse: `strncmp(arg, "--mode", 6)` at 5283.
- Enum: `typedef enum { MODE_TIME=0, MODE_RK4=1, MODE_EDMD=2, MODE_EDMD_HYBRID=3 } sim_mode_t;` (933).
- Default: `static sim_mode_t sim_mode = MODE_TIME;` (934).

| value | maps to (5285–5291) |
|---|---|
| `time` | MODE_TIME |
| `rk4` | MODE_RK4 |
| `edmd` | MODE_EDMD |
| `edmd-hybrid`, `hybrid`, `event-split`, `event_split` | MODE_EDMD_HYBRID |
| anything else | "Unknown mode", exit (5292) |

**MODE_RK4 is a label only.**
- `grep -n MODE_RK4` finds 6 lines in 00ALLINONE.c: 933 (enum), 2278 (run_params.json name), 5286 (parse), and 11880, 18566 and 19711 (name strings). It finds 0 lines in edmd_core/* and experiment_validation.c.
- Every other `\bsim_mode\b` line tests or assigns MODE_EDMD / MODE_EDMD_HYBRID, assigns MODE_TIME (10669, 15169, 21145), or tests MODE_TIME only at 15171, after the live scene reassigns `sim_mode`. So `rk4` takes the non-EDMD (time-stepped) path wherever the code branches on EDMD (e.g. 18003).
- `grep -in 'rk4\|runge'` finds no integrator. The comments at 651, 9816, 19333, 22832 and 22853 treat "TIME / RK4" as one path.

**Run-time remaps of the requested mode:**
- HYBRID → EDMD in the speed-of-sound run (18003–18006).
- EDMD/HYBRID → TIME in `run_szilard_engine_experiment` (21143–21146).
- HYBRID → TIME in the Szilard interactive setup (10666–10669).
- Live-scene APPLY sets the mode from the panel (15164–15170, 15215).
- EDMD → HYBRID in energy transfer if `num_internal_walls > EDMD_MAX_DIVIDERS` (18858–18861). EDMD_MAX_DIVIDERS is 32 (hspist3/edmd_core/edmd.h:13). This cannot fire from CLI input: `num_internal_walls` is clamped to `MAX_DIVIDER_CAPACITY` 32 (45, 5427–5428), and `--num-walls` accepts 1..5 (4589).

gen3 additionally requires `--mode=edmd`: `why = "gen3 is an event-driven engine (use --mode=edmd)"` (23601).

#### --engine and --edmd-acc
**--engine.**
- Only the exact strings `--engine=gen2` (5212, `cli_engine = 2`) and `--engine=gen3` (5214, `cli_engine = 3`) are accepted.
- Any other `--engine` or `--engine=...` hits `strncmp(arg, "--engine", strlen("--engine")) == 0 && (arg[8] == '=' || arg[8] == '\0')` and exits with "Invalid --engine '%s': use --engine=gen2 (default) or --engine=gen3." (5249–5251). This includes the space form `--engine gen3`.
- Default: `static int cli_engine = 2;` (1431).
- `--engine-replay=FILE` sets `cli_engine = 3` (5203). `main` then returns `run_gen3_replay(...)` (23495) before `init_output_dirs_from_argv0` (23497).

**gen3 checks in `main` (23592–23639).** Each failure goes to `g3_stop(why)` (23638). gen3 is refused:
- with the window outside the interactive scene, since experiments need `--headless` (23599–23600);
- with any mode other than `--mode=edmd` (23601);
- with `--single-test` (23602);
- for anything other than the speed-of-sound experiment, the energy-transfer experiment or the interactive scene (23603–23604);
- for the Szilard and particle-life scenes (23605–23608);
- in the scene with a heat bath or Andersen thermostat (23609–23610). Outside the scene the engine itself refuses a heat bath at build: `FAIL("edmd_gen3: the outer-wall heat bath is not implemented")` (hspist3/edmd_core/edmd_gen3.c:1119);
- with the gen2 schedule options (23611);
- with the gen3-option combinations at 23612–23637.

gen3-only flags given without gen3 exit 2 (23640–23667).

**--edmd-acc.** Default `static int cli_edmd_acc = 0;` (1355).
- Bare `--edmd-acc` (exact, 5266) is refused.
- `=1`, `=true`, `=yes` and `=on` are refused (5274–5276). The empty-value branch (5271–5273) is reachable only with an explicit empty next argument (cli_option_value, 3099–3106).
- `=0`, `=false`, `=no` and `=off` set `cli_edmd_acc = 0` (5277–5278).
- Anything else prints "Invalid edmd-acc" and exits (5280–5281).
- The space form `--edmd-acc 0` is refused, because the exact bare branch at 5266 fires first.

**The guard.** `static void refuse_edmd_acc_backend(void)` (1360–1367) runs `fprintf(stderr, "STOP: --edmd-acc selects the accelerated EDMD backend ..."); exit(2);` under `#ifndef HD_ALLOW_ACC_BACKEND` (1361).
- A second guard sits in `edmd_backend_create`: `if (cli_edmd_acc) refuse_edmd_acc_backend();` (1631). That function holds the only `edmd_acc_create` call (1650).
- `grep HD_ALLOW_ACC_BACKEND` over hspist3/Makefile and edmd_core finds nothing.
- The bin_a71f394 build line, hspist3/experiments_gen3_d52_261014/bin_a71f394/build_record.txt:30, has `-DBUILD_GIT="a71f394" -DBUILD_TARGET="mac-O3-e0pre" -DBUILD_CFLAGS="-O3 -ffp-contract=off"` and no `-DHD_ALLOW_ACC_BACKEND`. `strings` on that binary finds the guard text ("HD_ALLOW_ACC_BACKEND accepts it", 1 hit).
- The guard was introduced in commit 14ba1f22 (author date 2026-10-08), the oldest result of `git log -S'HD_ALLOW_ACC_BACKEND'`. That commit is an ancestor of HEAD (`git merge-base --is-ancestor`).

#### Table: geometry, dynamics, seeding, GUI, output
`(prefix n)` marks a strncmp of n characters. `= only` marks the branches described above. "Default" gives the global initialiser line, plus the macro line when the value comes from a macro. The active macro branch is 200–213, because `#define MOLECULAR_MODE 1` (140) and `#define EXPERIMENT_MODE 0` (141).

| flag | type | default (line) | what it changes: variable, then main reader | parse |
|---|---|---|---|---|
| **Dynamics / engine** | | | | |
| `--mode=V` | enum | MODE_TIME (934) | `sim_mode`; see --mode above | 5283 (prefix 6) |
| `--engine=gen2/gen3` | exact string | 2 (1431) | `cli_engine`. gen3 code is gated by `cli_engine == 3` (e.g. 1615, 1630, 1726). Restrictions at 23592–23639 | 5212, 5214, 5249 |
| `--engine-replay=FILE` | path, `=` only | NULL (1448) | sets `cli_engine = 3`; `main` returns `run_gen3_replay` (23495) | 5202 |
| `--replay-substeps=n` / `--replay-read-every=m` | int 1..1e6 / int 0..1e8 | 1 (1450) / 1 (1449) | `run_gen3_replay` (23402 / 23413) | 5208 / 5204 |
| `--edmd-acc[=v]` | bool | 0 (1355) | refused, see above | 5266, 5269 |
| `--legacy-resched` | switch | 0 (1368) | `edmd_set_legacy_resched` in `edmd_backend_create` (1635); refused under gen3 (23611) | 5260 |
| `--resched-audit[=all]` | switch | 0 (1369) | `edmd_set_resched_audit` (1636). 1 = after divider events, 2 = after all events (printf at 1643–1644). Refused under gen3 (23611) | 5262, 5264 |
| `--fixed-dt=v`, alias `--dt` | float > 0 | `fixed_dt_runtime = FIXED_DT` (929), 0.04f (201). In speed-of-sound without `--kbt1`, dt = FIXED_DT/√T (2334–2336) | `cli_fixed_dt` → `fixed_dt_runtime` in `update_dt_runtime_for_temperature` (2329–2331), called from `main` (23508). Headless adds one `fixed_dt_runtime` per loop iteration (22028) | 4524 (prefix 10 / prefix 4, so any `--dt...` matches) |
| `--substeps=N` | int 1..500 | -1, meaning SUBSTEPS = 1 (671; 207, 9210) | `cli_substeps_override`; time-mode step `update_particles_with_substepping` (9212–9213). Takes precedence over auto (9214) | 4573 (prefix 10) |
| `--substeps-auto`, `-auto-min`, `-auto-max`, `-auto-dx-frac` | switch, int 1..500, int 1..500, float (0,1] | false (674), 1 (675), 200 (676), 0.25 (679) | `update_particles_with_substepping` (9214); each of the three valued forms also sets `cli_substeps_auto` | 4544 exact; 4546, 4555, 4564 |
| `--pp-toi`, `--pp-toi-max-events`, `--pp-toi-neighbor-ring` | switch, int 1..1000, int 0..5 | 0 (1400), 4 (1401), 1 (1402) | time-mode step (9491); both valued forms also set `cli_pp_toi = 1` | 4530, 4532, 4538 |
| `--no-pp-collisions` | switch | false (513) | time-mode step (9490); speed-of-sound EDMD `prm.pp_collisions_enabled` (18029); `g3_m4b_init_driver` (17313). **Not honoured** by the energy-transfer and scene EDMD builds, which set 1 (18904, 21633). gen3 refuses 0: `FAIL("edmd_gen3: pair collisions cannot be switched off")` (edmd_gen3.c:1121) | 4114 |
| `--kbt1` | switch | 0 (1370) | `kB_effective()` returns 1/`temperature_runtime`, so kB·T = 1 (2260–2263). Without it, kB = `K_B` = `1/TEMPERATURE` with TEMPERATURE 100.0f (2264; 205, 203), so kB·T = 1 only while T = 100. Also disables the dt∝1/√T scaling (2334) | 4366 |
| `--temperature=T` or `T1,T2,...` | float > 0, or an unvalidated float list | `temperature_runtime = TEMPERATURE` = 100 (928, 203); list NULL (647) | A single value sets `temperature_runtime` (5531–5532). A list fills `cli_temperature_segments` and leaves T at its default (5097), read by `initialize_simulation` (6920) and `g3_m4_init_core` (17244) | 5092 (prefix 13) |
| `--hb-temp=T` | float > 0 | `heatbath_temperature` 1.0 (940); `heatbath_enabled` 0 (939) | sets the temperature (5536) and `heatbath_enabled = 1` (5115). Read by `handle_boundary_collision` (8426) and by speed-of-sound `prm.heatbath_*` (18022–18023). Under gen3 it is refused in the scene (23609–23610) and at engine build (edmd_gen3.c:1119) | 5106 (prefix 9) |
| `--wall-hold-steps=N` | int ≥ 0 | WALL_HOLD_STEPS = 10000 (297; 208) | `wall_hold_steps`; `update_wall` (8755) and `update_extra_walls` (8834) test `steps_elapsed >= wall_hold_steps`; bounds `--gen3-held-readout` (23621) | 4880 |
| `--auto-release` / `--auto-release-after-hold` | switch | 0 (1380) / 0 (1381) | `simulation_loop` (21968–21969, 21977 / 22132); `run_single_test_headless` (23118) | 4368 / 4370 |
| `--steps=N` | int > 0 | `num_steps = TIME_UNITS_SIMULATION / FIXED_DT` (241; 146: 3000), i.e. 75000 | `num_steps` (5509). Speed-of-sound: `speed_sound_steps_for_target_cycles` returns it when no oscillation target is set (16916, 16919), and a value < 1000 is replaced by 200000/dt (17799–17801). Energy transfer `target_steps` (18850); Szilard (21216). Not read by `simulation_loop` | 5025 (prefix 7) |
| `--time=T` | float > 0 | -1, unlimited (684) | `cli_time_span`, measured from the wall release: `simulation_loop` (22893), `run_single_test_headless` (23300) | 4581 (prefix 6) |
| `--record-sigma-time=T` | double (0,1e9], `=` only | 0 (1445) | `run_speed_of_sound_experiments` (17964) | 5194 |
| `--protocol=step/sigmoidal/linear/sinusoidal/optimal` | enum | PROTOCOL_STEP (645) | `cli_protocol`; `apply_piston_protocol` (7067), `configure_right_piston_protocol_from_cli` (12545) | 4714 |
| `--piston-right-protocol-mode=step/linear/ramp/sigmoidal` (+ `_` alias) | enum | "" (592) | `cli_protocol` and `cli_piston_right_mode`; energy transfer (19780) | 4700–4701 |
| `--piston-ramp-time` (+ `_`) | float > 0 | 0 (605) | `configure_right_piston_protocol_from_cli` (12529) | 4749–4750 |
| `--piston-speed`, `--piston-duration`, `--piston-travel-sigma` | float > 0 | 1.0 (584), 2.0 (586), 0 (588; "20% of box width" per the comment there) | `configure_right_piston_protocol_from_cli` (12526, 12533); `trigger_piston_step_with_protocol` (12610); energy transfer (19770–19778) | 4841, 4850, 4859 |
| `--max-right-piston-travel`, `--velocity-right-piston-step` / `-t0` / `-gradient` / `-max`, `--piston-right-duration`, `--piston-right-sigmoid-steepness` (each + alias) | float | 0 (594, 596, 598, 600, 602, 604, 608) | `configure_right_piston_protocol_from_cli` (12521–12539); energy transfer (19769–19779). Travel, step and max also set the generic travel or speed (4772–4773, 4784–4785, 4817–4818); duration sets `cli_piston_duration` (4829–4830) | 4761, 4774, 4786, 4796, 4806, 4819, 4831 |
| `--piston-work-target` (+ 3 aliases) | double > 0 | 0 (392) | `maybe_stop_right_piston_on_work_target` (2720) | 4868–4871 |
| `--auto-piston-step` | switch | false (472) | energy transfer (19406, 19561); cleared by `demo_begin_once` (2792) | 4698 |
| `--spring-wall=i` | int 0..4 | -1 (523) | `main` aborts if the compartment behind that wall holds particles. With no `--spring-wall`, `spring_constant = 0` (23549–23572) | 4658 |
| `--spring-k-sigma=K` / `--spring-k=k` | double (kT/σ²) / float (per px²) | 0 (526) | `cli_spring_k`; the sigma form stores K/24² (4652); `initialize_energy_measurement` (7184) | 4634 / 4653 (prefix 10) |
| `--spring-eq=x` | float | 0 (528) | `initialize_energy_measurement` (7173); energy transfer (19811) | 4667 (prefix 11) |
| `--gen3-cell-px=W` | double > 0, `=` only | 0 (1432), meaning `EDMD3_DEFAULT_CELL_PX` 32.0 (edmd_gen3.h:34). In the window scene the default is max(32, ⌈2R⌉) (1547) | `g3_cell_arg` (1546), `g3_build` (1552) | 5252 |
| `--gen3-compress-eta0` / `-melt` / `-speed` | double, `=` only; ranges at 5136 | 0.60 (1493), 500 (1494), 0.01 (1495) | `g3_m4b_init_driver` (17283, 17287); needs `--gen3-seeding=compressed` (23630) | 5130–5132 |
| `--cyclic-boundaries=0/1` / `--pl-dt=dt` | int / double | -1 (626) / -1 (627) | particle-life only: `pl_initialize` (8033–8034) / `simulation_loop` (22369) | 4465 / 4471 |
| **Geometry** | | | | |
| `--l0=L` | float > 0 | `L0_UNITS` 20.0f (251); override -1 (452); `cli_l0_d` 0 (1469) | `L0_UNITS` (5512) → `SIM_WIDTH = (int)(2 * L0_UNITS * PIXELS_PER_SIGMA)` (325), with PIXELS_PER_SIGMA 24 (202). Also the wall-position fraction (5446) and eta (5616). `cli_l0_d` gives the exact held-path box (18823, 18865). Speed-of-sound takes L0 from `--lengths` instead (17836, 17841) | 5035 (prefix 4) |
| `--height=H` | float > 0 | `HEIGHT_UNITS` 10.0f (252); override -1 (453) | `HEIGHT_UNITS` (5515) → `SIM_HEIGHT = (int)(HEIGHT_UNITS * PIXELS_PER_SIGMA)` (326) | 5046 (prefix 8) |
| `--lengths=L0,...` | float and double lists | NULL (438; 1463) | `run_speed_of_sound_experiments` (17757, 17838); sets `enable_speed_of_sound_experiments` (4966) | 4962 (prefix 9) |
| `--num-walls=N` | int 1..5 | 1 (643) | `cli_requested_walls` (5422–5424), then `num_internal_walls` (5426). Because the default is 1, this assignment always runs. No CLI value gives zero dividers: `rebuild_internal_walls`, called from `initialize_simulation_dimensions` (362), forces `num_internal_walls >= 1` (3571) | 4586 |
| `--walls=N` | int ≥ 1, clamped to `MAX_DIVIDER_CAPACITY` 32 (45) | `cli_requested_walls` 1 (401) | overwritten at 5423 by `cli_num_walls`, which is ≥ 1 at that point (643, 4589), so the flag has no effect. Preset wall counts set at 5345–5389 are overwritten the same way. Code reading, not run-tested | 4948 (prefix 7) |
| `--wall-positions=x1,...` | float list (σ from the left), at most 5, `atof` without validation | zeros (644; 1470) | applied only if `cli_num_walls > 0 && cli_wall_positions[0] != 0.0f` (5432): `preset_wall_fraction = x / (2 L0)`, clamped to [0,1] (5446–5452), then placed by `rebuild_internal_walls` (3579–3586). A missing entry stays 0 and puts that wall at XW1 (5436, 3581–3585). The double copy sets the gen3 held divider (18879) and is checked at 19060 | 4594 |
| `--wall-mass-factors=m1,...` | float list, at most 5 | {200,...} (640); set flag 0 (641) | `rebuild_internal_walls` repeats the last given value for later walls (3613–3617). Entry [0] also sets `wall_mass_runtime` (5522–5523), which overrides `--wall-mass-factor`. The energy-transfer EDMD builds read `cli_wall_masses[w]` directly, so missing entries get the initialiser 200, not the repeated value (18882–18885, 19199–19200) | 4609 |
| `--wall-mass-factor=f` | float > 0 | `wall_mass_runtime` = 1·200 (292, 247) | `wall_mass_runtime` (5518) | 5056 |
| `--wall-masses=f1,...` | int list | NULL (440) | speed-of-sound mass sweep (17764); enables experiments (4970) | 4967 (prefix 13) |
| `--wall-thickness=t` / `--wall-thickness-vis=t` | float > 0 | 0 (276; "will default to 2R at init" per the comment there) / -1 (287) | `set_wall_thickness_sigma` (5526; 282–284) / `wall_thickness_visual_runtime` (5529), drawn by `draw_wall` (8672) | 5082 (prefix 16) / 5066 (prefix 20) |
| `--particle-radius=r` / `--particle-diameter=d` | float > 0 (σ) | r = 0.5σ: `PARTICLE_RADIUS_UNIT` 0.5 (172) × `particle_scale_runtime` 1.0 (237, 209) | `particle_scale_runtime = r / 0.5` (5631, 5635), which feeds the `PARTICLE_RADIUS` macro in px (239). Precedence: eta, then radius, then diameter (5539, 5627, 5632) | 4888 / 4897 |
| `--eta=η`, alias `--packing-fraction` | float, 0 < η < 0.9069 | unset (466–467) | radius from η: from the limiting segment when per-box counts are given (5551–5612), otherwise global (5614–5624); raises `cli_phi_max` (5641–5650) | 4906–4907 (prefix 5) |
| `--phi-max=φ` | float in (0, 0.9069] | 0.78 (541) | per-segment capacity caps in all four seeders (6086, 6260, 6510, 6722). Overpacking exits 1 only when neither `--eta` nor `--particle-radius` is set (e.g. 6089–6093, 6263–6267) | 4940 |
| `--gen3-exact-box` | switch | 0 (1462; tentative declaration at 323) | speed-of-sound `prm.boxW = 2·L0d·24` (18012) and divider (18018); held path through `g3_held_exact` (18639) | 5155 |
| `--calibrate-work=t,tol,min,max` | 4 floats | unset (530–534) | places the last wall when there are ≥ 2 walls (5459–5506). It runs inside the parser, before `initialize_simulation_dimensions` (first call 23504). So `SIM_WIDTH` is still 0 at 5461 (zero-initialised at 315), and `L0_UNITS` is still the default at 5492 and 5498 (overwritten only at 5512). Code reading, not run-tested | 4672 (prefix 16) |
| **Seeding** | | | | |
| `--seed=N/time` | uint or `time` | `cli_seed` 1 (499); `cli_seed_set` false (498) | `main` calls `srand(cli_seed)` only if set (23501–23503). Speed-of-sound always calls `srand(run_seed)`, derived from `cli_seed` unless `--speed-sound-exact-seed` is given (17830–17834). The held lattice path calls `srand(cli_seed)` (18822). The seed is in the `--dist-log` file name (2230, 2241) and is required on the held exact box (23614) | 4343 (prefix 6) |
| `--seed-drift-order=old/drift-first` | enum | 0 = old (1379) | `equalize_temperature_per_segment`: drift removed before or after the kT rescale (2566–2584) | 4338 (before `--seed`, so reachable) |
| `--seeding=grid/honeycomb/random` | enum | SEEDING_GRID (511) | **unreachable** (see shadowing). `cli_seeding_mode` (switch at 6655) changes only internally (10605, 15220, 15228) | 4625 (dead) |
| `--particles=N` | int > 0 | -1 (470); `particles_active = DEFAULT_PARTICLES_ACTIVE` = 100 (156, 155) | `Nreq` in the seeders (6078, 6253, 6513, 6725), capped at NUM_PARTICLES 120000 (e.g. 6079; 154), then `particles_active` (e.g. 6095) | 4932 (prefix 11) |
| `--particles-boxes=a,b,...` | float list | NULL (536, 538) | any entry > 1 switches to absolute counts (4482–4494), otherwise fractions. Normalised into `preset_segment_fraction` (5409–5419, 5682–5691). Counts are read by the seeders (6119, 6290, 6732) and the spring check (23552) | 4477 (4917 dead) |
| `--particles-box-left=n` / `-right=n` | int | -1 (636–637) | two-box fractions (5399–5408, 5672–5681); eta two-box branch (5583) | 4499 / 4504 (4921 / 4926 dead) |
| `--distribute=area/equal` | enum | false (515) | `initialize_simulation_segmented_grid` (6312), `initialize_simulation` (6765) | 4359 (prefix 12) |
| `--left-empty` | switch | false (512) | seeders (6347, 6772); also set by the energy_transfer preset (4424) | 4112 |
| `--gen3-seeding=lattice/compressed/driver` | enum, `=` only | 0 = driver (1457) | `g3_build` (1584), speed-of-sound (17862), `g3_held_lattice` (18641); gen3 only (23648–23651) | 5124 |
| `--gen3-lattice=auto/parallel/perpendicular`, `--gen3-lattice-orientation=0/turned`, `--gen3-jitter=f` | enum, enum, double in [0, 0.5); `=` only | -1 (1458), -1 (1459), 0.25 (1461) | `g3_m4_seed`: `g3_lattice_fit(..., cli_gen3_orient, ...)` (17188) and `g3_lattice_place(..., cli_gen3_jitter)` (17200); `g3_m4_init_core` passes `cli_gen3_cryst` (17222). Orientation needs lattice seeding (23636) | 5140, 5146, 5151 |
| `--gen3-restart=PATH` | path, `=` only | NULL (1480) | `g3_restart_apply` reads the checkpoint (17654); speed-of-sound only (23623) | 5178 |
| `--speed-sound-exact-seed=N` | uint | unset (490–491) | replaces the per-trajectory seed (17830–17833) | 4120 |
| **GUI** | | | | |
| `--show-simulation` | switch | `cli_force_no_experiments` false (468); `enable_speed_of_sound_experiments` false (165) | sets both: `cli_force_no_experiments = true; enable_speed_of_sound_experiments = false;` (4107–4108). `main` skips every experiment branch (23672–23678) and calls `simulation_loop()` (23681). `--single-test` still wins (23670). This sets `run_scene` (23579–23583), which gates the `--gui-*` flags (23584–23591) | 4105 |
| `--no-experiments` | switch | false (468) | same end state; `enable_speed_of_sound_experiments` is cleared at 5652–5654 | 4103 |
| `--no-gui`, `--headless` | switch | false (471) | no SDL or TTF init (23514–23536). `simulation_loop` auto-starts (21977–21979) and adds one `fixed_dt_runtime` per iteration with no event polling (22026–22028) | 4109 |
| `--gui-time-scale=X` | double 0.01..10000, `=` only | 0.0 = unchanged (1501); `time_scale_runtime` 1.0 (930) | sets `time_scale_runtime` at `simulation_loop` entry (21996). The window path adds `frame_time * time_scale_runtime` to the accumulator (22006), with `frame_time` clamped to 0.1 s (22005). Accepted with `--headless` but unused there (22028). Refused outside the scene (23584–23587) | 5222 |
| `--timescale=X` | n/a | n/a | **unreachable**; parsed as `--time` (4581) | 5116 (dead) |
| `--gui-stats`, `--gui-run-seconds=S`, `--gui-keys=STEP:KEY,...` | switch, double > 0, string | 0 (1499), 0 (1500), NULL (1503) | `simulation_loop` (22991 / 22031); `g3_gui_keys_parse` (call at 21995, read at 21316); `--gui-keys` with `--headless` exits (21317) | 5216, 5218, 5226 |
| `--gui-window=WxH`, `--gui-ui-scale=U`, `--gui-shot=STEP:F.bmp/limit:F.bmp` | 2 ints (320x240..16384, 5231), float 1..4, at most 4 (1506) | 0×0 (1505), 0 = measured (774), none (1507) | `initSDL` (5766); `gui_sync_renderer_canvas_size` (813); `gui_shot_maybe` (21548). Refused when headless or outside the scene (23588–23591) | 5228, 5234, 5238 |
| `--render=experiment/paper/debug` | enum | 0 = experiment (1388) | `simulation_loop` (22910) | 4725 |
| `--demo` / `--demo-shot=PATH[,STEPS]` | switch / string | 0 (1387; tentative declaration at 1062) / "" and 0 (1391–1392) | `demo_update_stage` (2745), `demo_echo_stage_change` (2819), `keysSimulation` (12135, 12152), `simulation_loop` (22088, 22978) | 4734 / 4736 |
| `--colourcode`, `--colorcode=none/temperature/velocity` | enum | COLORCODE_NONE (2272) | `render_particles` (7421, 7433, 7492–7500) | 5293 |
| `--memory-view`, `--memory`, `--szilard-view=species/memory/overlay` | enum | SZ_VIEW_SPECIES (2162) | `render_particles` (7453), `render_energy_hud` (16231) | 5299–5301 |
| `--simple-box-render-every=N`, `--simple-box-gui-delay-ms=N` | int, clamped ≥ 1 / ≥ 0 | 1 (2096), 8 (2095) | `simple_box_gui_effective_render_every` (20545) / `simple_box_render_preview` (20642) | 4328, 4333 |
| **Output** | | | | |
| `--output-dt=T` | float ≥ 0 | 1.0 (654); semantics per the comment at 649–653 | `update_particles_with_substepping` (9819), `simulation_loop` (22514, 22655–22656, 22845–22889), `run_single_test_headless` (23135), `run_params.json` (2310) | 4509 (prefix 11) |
| `--trace-every=N` | int ≥ 1 | 1 (669) | energy-transfer trace decimation (18742, 19477) | 4514 |
| `--energy-transfer-trace=PATH` | path | `<g_energy_transfer_dir>/energy_transfer_trace.csv` (18781). That dir is `<dir of argv[0]>/experiments_energy_transfer` (3419), or `experiments_energy_transfer` relative to cwd when argv[0] has no '/', its dir is empty, or its dir is '.' (3401, 3404–3408, 3413–3416) | a relative path is resolved against the launch cwd (18779; 3424–3433). `00_COMMAND.md` goes into the directory part of the given path (18754–18762). `mkdir_p(g_energy_transfer_dir)` runs regardless (18748) | 4690 |
| `--energy-transfer-summary=PATH` | path | `<g_energy_transfer_dir>/energy_transfer_runs_%dwalls.csv` (18574–18575) | a relative path is resolved against the launch cwd (18572; 3424–3433) | 4682 |
| `--speed-sound-run-dir=PATH` / `--speed-sound-log-stride=N/auto` | path / int ≥ 0 or `auto` | `experiments_speed_of_sound/<EDMD or TIME>/mode1_normalized_units/<timestamp>`, relative to cwd (17720–17721, 17729–17730) / 1 (489). `auto` = 0 gives stride ⌊samples per period / 256⌋ (17952–17960; 488) | mode1 dir (17726–17727). The mode0 dir `experiments_speed_of_sound/.../mode0_real_units/<timestamp>` is always created relative to cwd and receives `00_COMMAND.md` (17724–17725, 17734, 17738). Stride read at 17951 | 4116 / 4132 |
| `--dist-log`, `--dist-dt`, `--dist-bins`, `--dist-max-speed`, `--dist-log-path` | switch, float, int (< 1 → 50, > 2000 → 2000; 4187–4188), float, path | 0 (1403), 0 (1404), 50 (1405), 0 (1406), NULL (1408) | `simulation_loop` with the energy-transfer preset (21903–21915, 22664, 23010); `build_dist_log_path` (2223) | 4179, 4181, 4184, 4189, 4192 |
| `--szilard-log=PATH`, `--szilard-log-dt=T` | path, double clamped to [0, 1e9] | NULL (2082), 0.2 (2083) | `simulation_loop` (21722, 21886) | 4258, 4252 |
| `--simple-box-run-dir`, `--simple-box-run-label`, `--simple-box-log-dt` | path, string, double ≥ 0 | NULL (2097), NULL (2098), 0.25 (2088) | `simple_box_build_run_dir` (20457, 20473); simple-box run (20805) | 4320, 4324, 4281 |
| `--gen3-dump-initial=PATH` | path, `=` only | NULL (1473) | speed-of-sound (18083), energy transfer (18955) | 5157 |
| `--gen3-checkpoint=hold:K:PATH/record:K:PATH`, `--gen3-checkpoint-stop` | spec (`=` only), switch | -1 / -1 / NULL (1476–1478), 0 (1479) | `g3_ckpt_write` (17627, 17645); speed-of-sound only (23623–23625) | 5167, 5176 |
| `--gen3-virial-blocks=T`, `--gen3-snapshots=T` | double (0, 1e9], `=` only | 0 (1484), 0 (1485) | `g3_obs_open` (17484, 17493); speed-of-sound only (23626) | 5180–5181 |
| `--gen3-psi6-compartments`, `--gen3-psi6-snapshots=FILE`, `--psi6-every=T` | switch, path, double | 0 (1488), NULL (1490), 0.25 (1451) | `g3_psi6_comp_header` (16982); `main` returns `run_gen3_psi6_snapshots` (23496); `g3_psi6_series_sample` (17011) | 5186, 5188, 5190 |
| `--gen3-held-readout=K`, `--gen3-audit-every=K`, `--validator-every=K` | int ≥ 1, long ≥ 1, int 1..1e8; `=` only | 0 (1471), 0 (1447), 0 (1443). The validator default resolves to 1 under gen2 and `G3_VALIDATOR_EVERY` 60 under gen3 (1454–1455) | energy transfer (19253); `g3_build` (1581); `validator_every()` read at 17595, 18176, 18227, 18740 | 5159, 5198, 5256 |
| `--edmd-debug-particles=A,B`, `--edmd-debug-event-history=N` | int pair, int ≥ 0 | -1, -1 (495–496), 0 (497); set to 50 when A,B is given and it is ≤ 0 (4166) | speed-of-sound (18033–18034) | 4146, 4167 |
| `--quiet`, `--energy-measurement` | switch | 0 (1397), false (469) | printers, e.g. 5701; `initialize_energy_measurement` (7165) | 4680 (prefix 7), 4364 |
| `--version`, `--help` | switch | n/a | `--version` prints `BUILD_GIT`, `BUILD_TARGET`, `BUILD_CFLAGS` (4095; "unknown" defaults at 659–667) and exits 0. `--help` runs `print_cli_usage` and exits with EXIT_SUCCESS | 4094, 4098 |

#### Remaining flags (outside those categories), with parse line
- `--experiments` 4101.
- `--experiment=NAME[:variant]` 4374. Presets (4396–4456): 2wall_exp_with_1_wall, 2wall_exp_with_2_wall[:1/2], wall_mid, speed_of_sound, energy_transfer (sets `--left-empty` and `--energy-measurement`, 4424–4425), simple_prediction_box / simple_gas_box, multi_wall, szilard_engine, szilard_engine_particle_separation, atp_synthase, protocol_optimization, dynamic_ib, particlelife.
- `--single-test` 4372.
- `--repeats` 4971.
- Oscillations: `--target-oscillations` / `--oscillations` 4982–4983; `--oscillation-safety` 4994; `--oscillation-min-steps` 5005; `--oscillation-max-steps` 5015.
- `--pl-tensor` 4458.
- Szilard: `--szilard-no-pp` 4177; `--szilard-equil-steps` 4220; `--szilard-sep-duration` 4224; `--szilard-species-ratio` 4228; `--szilard-perm-swap` 4234; `--szilard-unstick` 4237; `--szilard-unstick-hits` 4240 (dead); `--szilard-unstick-vx-frac` 4246 (dead); `--szilard-n` 4262.
- Simple box: `--simple-box-cycles` 4266; `--simple-box-equil-steps` 4271; `--simple-box-relax-steps` 4276; `--simple-box-target-protocol` 4286; `--simple-box-target-min` 4295; `--simple-box-target-max` 4298; `--simple-box-initial-target` 4301; `--simple-box-corr-alpha` 4304; `--simple-box-noise` 4310; `--simple-box-period` 4315.
- Efficiency: `--eff-window` 4196; `--eff-stop-after-window` 4203; `--eff-output` 4205.

#### Count
- 188 distinct spellings appear in the chain. These include aliases, the three `--engine` spellings and both `--resched-audit` spellings.
- Merging 22 alias spellings gives **166 distinct flags**: 135 in the table and 31 in the remaining list. I checked the two lists against the 188 spellings: every spelling is covered.
- 4 of the 166 flags are unreachable.
- The 6 HD_* environment variables are not counted.

#### Binary runs
**Verifier note.** These five runs were made by the drafting agent. They were not repeated during verification, because the verifier was not allowed to execute the binary. Read-only checks agree with them:
- the format strings at 4095, 4353 and 4584;
- the usage printer's 161 `\n` with no conditionals (3907–4069);
- build_record.txt:30;
- `strings` on the binary, which contains `00ALLINONE  git %s  target %s`, `Invalid seed '%s'. Use an integer or 'time'.`, `Invalid time '%s'.` and the guard text.

Runs as reported: cwd `PHYSICALINFORMATIONDYNAMICS/_m0_scratch/flags_run`, `SDL_VIDEODRIVER=dummy`, `perl -e 'alarm 60; exec @ARGV' -- <abs path>/hspist3/experiments_gen3_d52_261014/bin_a71f394/00ALLINONE <args>`.
1. `--version`: printed `00ALLINONE  git a71f394  target mac-O3-e0pre` / `CFLAGS: -O3 -ffp-contract=off`; exit 0.
2. `--help`: 161 lines; exit 0. It includes the `--seeding=MODE` line (3949), the `--timescale=value` line (4017), `--edmd-acc=0|1 ... DISABLED ... 0 only` (4064) and `--mode=NAME  time (default) | rk4 | edmd | edmd-hybrid | hybrid | event-split` (4067). It has no `--engine`, `--gen3-` or `--gui-` line.
3. `--seeding=grid`: printed `Invalid seed 'grid'. Use an integer or 'time'.`; exit 1.
4. `--timescale=0`: printed `Invalid time '0'.`; exit 1.
5. `--edmd-acc=1`: printed `STOP: --edmd-acc selects the accelerated EDMD backend (edmd_core/edmd_accelerated.c), which misses collisions ... only a study build compiled with -DHD_ALLOW_ACC_BACKEND accepts it.`; exit 2.

All five exit inside `parse_cli_options`, before `init_output_dirs_from_argv0` (23497) and `write_command_md_to_dir(".", ...)` (23500).

At verification, `ls -la` shows `_m0_scratch/flags_run` empty, created 13:11. In bin_a71f394/, the binary is dated Oct 14 22:19, and 00ALLINONE.Wall.txt, 00ALLINONE.warnings.txt and build_record.txt are dated Oct 15 11:20. Both are before 13:11.

No run got past the parser. A run past the parser writes:
- `00_COMMAND.md` into cwd (23500);
- `run_params.json` into cwd (2321, from `initialize_simulation` 6925);
- in the scene, `energy_log.csv`, `wall_position.csv` and `wall_position_run_params.json` into cwd (21707–21714, 21944–21947).

It also points the energy-transfer, simple-box and Szilard dirs at `<binary dir>/experiments_*` (3418–3420). So the following findings come from code reading only: the `--walls` override, the preset wall-count override, `rk4` = time path, `--calibrate-work`, the `--wall-mass-factors` mismatch, missing `--wall-positions` entries, and `--no-pp-collisions` in energy transfer.

## B. Scene model

### B1

Paths: `edmd.h`, `edmd.c`, `edmd_accelerated.[ch]`, `edmd_gen3.[ch]` are in `hspist3/edmd_core/`; `00ALLINONE.c` is `hspist3/00ALLINONE.c`; `tests/` is `hspist3/edmd_core/tests/`.


**Public engine types (edmd.h), used by the GUI and all three engines**

| type | path:line | fields |
|---|---|---|
| `EDMD_Params` | edmd.h:15-94 | 35 declarators (below) |
| `EDMD_Particle` | edmd.h:96-100 | `double x, y;` `double vx, vy;` `int coll_count;` |
| `EDMD` (opaque) | `typedef struct EDMD EDMD;` edmd.h:102 | two different definitions behind one typedef: edmd.c:128-195 (44 fields), edmd_accelerated.c:92-117 (18 fields) |
| `EDMD3` (opaque) | `typedef struct EDMD3 EDMD3;` edmd_gen3.h:39; defined edmd_gen3.c:134-175 | 81 fields |
| `EDMD3_Health`, `EDMD3_Audit`, `EDMD3_Ledger`, `EDMD3_Tol` | edmd_gen3.h:41-67, 72-93, 101-109, 112-122 | telemetry, schedule audit, momentum/energy ledger, tolerances |
| `EDMD3_Coll` | edmd_gen3.h:184 | `typedef struct { int kind, a, b, comp_a, comp_b; double nx, ny, vax, vay, vbx, vby; } EDMD3_Coll;` |

`EDMD_Params` (D = `EDMD_MAX_DIVIDERS` = 32, edmd.h:13); names and types as in the source, regrouped onto fewer lines:
```
double boxW, boxH; double radius; int N; int max_events_hint; double cell_size;
int divider_count; double divider_x[D]; double divider_thickness[D]; double divider_mass[D];
double divider_vx[D]; double divider_k[D]; double divider_xeq[D];
int has_pistonL; double pistonL_x; double pistonL_vx; double pistonL_mass;
int has_pistonR; double pistonR_x; double pistonR_vx; double pistonR_mass;
int pp_collisions_enabled;
int heatbath_enabled; double heatbath_temperature; int thermal_wall_mode; double mb_overshoot_factor;
double stability_window_percent; double particle_mass; double kB;
const int* species;
int divider_gate_mode[D]; int divider_gate_species[D]; int divider_gate_target_side[D];
double divider_gate_hot_ratio[D]; double divider_gate_speed_ref[D];
```
- One radius for all disks: `double radius;         /* disk radius (all equal) */` (edmd.h:17).
- Mass 0 = infinite for dividers and pistons: per comment edmd.h:27, 36, 41; the code takes the prescribed-velocity branch for `M <= 0.0` (edmd.c:1388-1399, 1414-1427; edmd_gen3.c:830-838).
- No per-disk radius, mass or fixed flag: `EDMD_Particle` (edmd.h:96-100). `EDMD_Params` has one `particle_mass` (edmd.h:52); the gen2 body resolvers use `double m = 1.0;` (edmd.c:1387, 1413).
- The only per-disk attribute is the external `species` pointer (edmd.h:57); gen3 refuses it (edmd_gen3.c:1120).

**Engine-internal structs**

- **gen2, edmd.c.**
  - `EvType` (65): `EV_AB=0, EV_WL, EV_WR, EV_WB, EV_WT, EV_DL, EV_DR, EV_PL, EV_PR`.
  - `Event` (82-87): `double t; int a, b; int ca, cb; EvType type;`. `Cell` (116-119): `int *idx; int count, cap;`. `Heap` (122-125). `EdmdTraceRec` (99-107, debug tracer).
  - `struct EDMD` (128-195), 44 fields: `prm, t, P, work_divider[D], work_pistonL, work_pistonR, heat_bath, forced_advance_count, clamped, clamped_count, clamp_repair_count, overlap_repair_count, wall_overdue_count, virial_accum, virial_t0, virial_pair_events, wall_impulse[4], wall_events[4], wall_thermal_events, grid, gw, gh, cell_size, heap, div_epoch[D], ...`. 19 omitted: heap-compaction control and telemetry (`heap_compact_at`, `heap_compactions`, `past_event_count`, 182-184; `heap_compact_at` steers compaction, 1785-1786), contact audit (186-187), schedule audit (189-194).
  - No divider or piston struct. The bodies are the parallel arrays in the engine's own copy `S->prm`, and the engine writes them: `S->prm.divider_x[d] += S->prm.divider_vx[d] * dt;` (993), `S->prm.divider_vx[d] = v2` (1407), `S->prm.pistonL_vx = v2` (1432).
- **accelerated, edmd_accelerated.c.**
  - `EvType` adds `EV_CC` (54). `Event`, `Cell`, `Heap` as gen2 (72-89).
  - `struct EDMD` (92-117), 18 fields: `prm, t, P, work_divider[D], work_pistonL, work_pistonR, heat_bath, forced_advance_count, clamped, clamped_count, clamp_repair_count, overlap_repair_count, wall_overdue_count, grid, gw, gh, cell_size, heap`.
  - Its public functions are renamed by macros, e.g. `#define edmd_create edmd_acc_create` (edmd_accelerated.c:4-13); API in edmd_accelerated.h:10-51.
- **gen3, edmd_gen3.c.**
  - `Ev3` (91-96): `real t; int type, a, b; int ca, cb; int pad;`. `Heap3` (98). Event types: `enum { T_CROSS = 0, T_BAND = 1, T_WALL = 2, T_DIV = 3, T_PISTON = 4, T_PAIR = 5 };` (89).
  - `Disk3` (100-111): `real xi, zeta; real vx, vy; real tau; int cx, cy; int cnt; int slot; int last; int comp;`. Position relative to the cell's lower-left corner (per comment 101; code 1224), `tau` its time stamp (per comment 103).
  - `Obj3` (117-132), indices: dividers 0..31, `OBJ_PL` = 32, `OBJ_PR` = 33, `NOBJ` = 34 (114-116): `int active, harmonic; real M; real k, xeq, omega; real th, h; real x, v, tau; int epoch, band_gen; int blo, bhi; real t_band; int last; int pad; real Jf[2]; long nf[2]; real work; real J_inf; real J_spring;`.
  - `Harm` (330), spring-contact helper.
  - `struct EDMD3` (134-175), 81 fields: `prm, N, R, d, d2, w, boxW, boxH, gw, gh, ncell, ccap, cell, ccount, D, heap, compact_at, T0, now, ..., obj[NOBJ], objs[NOBJ], nobj, ..., evlog, evlog_tscale`. 57 omitted: validator cadence (144), output buffer `out` (145), health `H` (146), validator tolerances (147), virial incl. per compartment (148-149), wall impulse/events (150), event hash (151), stagnation guard (152), fatal flag/message (153), audits and their scratch (155-162), body tolerances and ledger accumulators (166-168), body contact audit (169), body audit scratch (170).
  - `real` is `double` unless `-DEDMD3_LONG_DOUBLE` (79-87).

**GUI program, 00ALLINONE.c: no particle, wall, piston or spring struct used by the physics; file-scope parallel arrays and scalars.**

- **Particles** (748-755), size `#define NUM_PARTICLES 120000` (154), active count `particles_active` (156): `double X[NUM_PARTICLES]; double Y[NUM_PARTICLES]; double X_old[NUM_PARTICLES]; double V_init[NUM_PARTICLES]; double Vx[NUM_PARTICLES]; double Vy[NUM_PARTICLES]; double Radius[NUM_PARTICLES];`. Every `Radius[...] =` assignment found sets `PARTICLE_RADIUS` (6201, 6419, 6445, 6564, 6583, 6857, 6877, 17227, 17383).
- **Primary divider:** `float wall_x_old;` `float wall_x;` `float vx_wall = 0.0;` `int wall_enabled = 0;` (693-696); `wall_mass_runtime` (292); `wall_hold_enabled`, `wall_is_released`, `wall_hold_steps`, `steps_elapsed` (295-298); `wall_impulse_x_accum` (299).
- **Multi-wall:** `num_internal_walls` (400), `cli_requested_walls` (401), `primary_wall_index`, `leftmost_wall_index`, `float *all_wall_positions`, `float *all_wall_masses`, `extra_wall_count`, `float *extra_wall_positions`, `*extra_wall_old_positions`, `*extra_wall_velocity`, `*extra_wall_masses`, `bool *extra_wall_hold_enabled`, `*extra_wall_is_released`, `double *extra_wall_release_time` (402-413). Per-compartment statistics `segment_*` (414-423).
- **Pistons:** `float piston_left_x = 0;` `float piston_right_x;` `float vx_piston_left = 0;` `float vx_piston_right = 0;` (367-370); step state (371-374); work and hit counters (377-380); protocol parameters `piston_protocol_*` (572-580); `typedef enum {...} PistonProtocol` (562-569).
- **Spring:** `typedef struct {...} EnergyMeasurement` (2634-2654), 18 fields, starting `float spring_constant; float equilibrium_position; double spring_force; double spring_energy; ...`; single instance `static EnergyMeasurement energy_measurement` (2659).
- **Heat bath:** scalars `heatbath_enabled, heatbath_temperature, thermal_wall_mode, mb_overshoot_factor, stability_window_percent, current_sample_temp` (939-944); Andersen thermostat `andersen_enabled, andersen_collision_freq, last_andersen_collision_step` (952-954).
- **Other structs:**
  - `MultiWallConfig` (611-620, 8 fields): defined, never used (grep `MultiWallConfig` finds only 620).
  - `SzilardBranchSnapshot` (2164-2213; static instance 2215) copies globals, incl. `double X[NUM_PARTICLES]`, `Y`, `Vx`, `Vy` and two `int` arrays of `NUM_PARTICLES` (2207-2212).
  - `LiveSlider` (967-976) has `void *target` pointing at globals.
  - Analysis/seeding: `Psi6Result` (16568-16573), `Psi6Comp` (16636), `Psi6Accum` (16681-16684), `G3Lattice` (17076), `G3Obs` (17469-17475), `G3RunRefs` (17566-17570), `G3RunCkpt` (17571-17584), `SimpleBox*` (19866-19915), anonymous `static struct {...} rm` (11956-11964).

**Shared and converted**

- **Shared:** `EDMD_Params` and `EDMD_Particle` are used by the GUI, gen2, accelerated and gen3 (gen3 via `#include "edmd.h"`, edmd_gen3.h:28; accelerated via edmd_accelerated.h:8). Each engine copies the params at create: `S->prm = *prm_in;` (edmd.c:1443; edmd_accelerated.c:1021), `S->prm = *prm;` (edmd_gen3.c:1134).
- **GUI to engine (copy-in), inline loops, no named function:**
  - Params field by field in `simulation_loop` (00ALLINONE.c:21592-21643), e.g. `prm.boxW = (double)(XW2 - XW1);` (21593).
  - Particles are written into the engine's own array through a const-cast of `edmd_particles()` (declared `const`, edmd.h:142): 21655-21667 (`double bx = (double)(X[i] - XW1);` 21658, clamped to [R, boxW-R] at 21660-21661), 22349-22356 (hybrid write-back), 10971-10976 (`szilard_rebuild_edmd_from_current_globals`), 19314-19318 (energy transfer, hybrid).
  - Left piston face offset: `double xL = (double)((piston_left_x + 5.0f) - XW1);` (21646).
- **Engine to GUI (copy-out):**
  - Particles: `X[i] = (double)XW1 + (double)P[i].x;` at 22163, 22287, 22341, 19262, 20388, 23210.
  - Bodies: `piston_left_x  = (float)(XW1 + ep->pistonL_x - 5.0);` (22170), `wall_x = (double)XW1 + ep->divider_x[0];` (22173), `vx_wall = (float)ep->divider_vx[0];` (22174). Targets are `float` (367-370, 693-695): the doubles are narrowed.
- **Dispatch:** `#define edmd_create edmd_backend_create` and 19 more defines (2037-2056) redirect later calls of those 20 names to the dispatch functions at 1629-1860. Other `edmd_*` functions are not redirected and call gen2 directly, e.g. `edmd_set_divider_springs(g_edmd, 2, k, xeq);` (11239).
- **gen3 conversion** (design per comment 1418-1430; the gen2 container stays the setup store):
  - `g3_build` (1549-1594) calls `edmd3_create(p, g3_cell_arg(p), err, sizeof err)` (1552) and `edmd3_load(g_edmd3, edmd_particles(S), t0, err, sizeof err)` (1555); it runs at the first advance (1725-1727).
  - Readers: `g3_particles` (1599-1607), a memcpy from `edmd3_particles`; `g3_params` (1608-1614) refreshes divider and piston x, vx from `edmd3_divider_state` / `edmd3_piston_state`.
  - A write into the handed-out copy is caught by `g3_check_untouched` (1595-1598), called at the next advance (1729) and reschedule (1675).
  - Body changes: `edmd_backend_set_divider_motions` (1698-1712), `g3_piston` (1616-1626). Repositioning a body stops the run (1621, 1694).
- **Inside gen3:** `EDMD_Particle` to `Disk3` in `edmd3_load` (edmd_gen3.c:1218-1230): `A->xi = P[i].x - (real)cx * S->w;` (1224). `EDMD_Params` to `Obj3` in `edmd3_create` (1160-1167). `Disk3` to `EDMD_Particle` in `edmd3_particles` (1303-1311): `S->out[i].x = (real)A->cx * S->w + x;` (1306).

### B2

Paths as in B1.

**Asked explicitly: how gen3 represents walls.**
- There is no wall struct and no wall list. The four outer walls are a side index s = 0..3, "outer wall b: 0 L, 1 R, 2 B, 3 T" (per comment, `edmd_gen3.c:10-11`), i.e. x = 0, x = boxW, y = 0, y = boxH.
- s is hard-coded in `switch (s)` blocks and `for (int s = 0; s < 4; ++s)` loops:
  - the predictor `wall_gap` (439-446);
  - the reach test `wall_candidate` (265-272);
  - the resolver `exec_wall` (737-756);
  - the validator (585, 640);
  - the schedule audit (1050-1056);
  - the ledger copy `L->J_wall[s]` (1533).
- A wall event carries s in its `b` field: `heap_push(S, S->now + dt, T_WALL, i, s, A->cnt, 0);` (453).
- Internal bodies are `Obj3` entries (117-132): divider slabs with half-width `O->h = 0.5 * O->th + S->R;` (1162) and piston faces with `O->h = S->R;` (1166-1167). Each is located by its x position alone (`body_gap`, 327).
- So a polyline or arc wall cannot be expressed as data. It is a new event type (C5, H option d).


**Primitives present**

| primitive | gen2 edmd.c | accelerated edmd_accelerated.c | gen3 edmd_gen3.c | time-stepped / hybrid (00ALLINONE.c) |
|---|---|---|---|---|
| Outer box: 4 fixed walls x=0, x=boxW, y=0, y=boxH | predictors 518-536; resolver 1328-1349 | 292-310; 930-931 | `wall_gap` 439-446; `exec_wall` 737-756 | `handle_boundary_collision` 8414-8524 (XW1, XW2, YW1, YW2) |
| Vertical divider slab (centre x, thickness th) | `collide_time_divider_L/R` 670-766 | 340-417, 418-495 | `Obj3` + `body_rule` 352-402; `O->h = 0.5 * O->th + S->R` (1162) | `ccd_wall_step` 8997ff; `hybrid_faces_correction_pass` 9916-10004 |
| Divider mass (particle units, 0 = infinite) | resolve 1387-1410 | 972-990 | `exec_obj` 824-838 | `wall_mass_runtime` (292); `cli_wall_masses[5]`, default 200 (640) |
| Divider spring k, x_eq | `divider_has_spring` 310-315 (needs k > 0 and mass > 0, 314); analytic orbit `harmonic_advance_1d` 456-470; contact `harmonic_first_contact` 551-640, bracketed at the derivative's zeros (566-582), bisection (597-604) | spring branch 346-409: coarse scan (369-380), fine scan (383-394), 80-step bisection (399-404) | `O->harmonic` (1171); `obj_at` 297-305; closed-form brackets + safeguarded Newton in `body_rule` 362-401 | `update_wall` 8775-8782: `vx_wall += a * dt;` (8779), then `wall_x += vx_wall * dt;` (8782) |
| Hold, then release | engine: mass 0 and v 0 while held; driver `div_mass[w] = 0.0; div_vx[w] = 0.0;` (00ALLINONE.c:19212-19214) | same driver | through `edmd3_set_divider_motion` (00ALLINONE.c:1698-1707) | `update_wall` 8741-8765, release at `steps_elapsed >= wall_hold_steps` (8755); `#define WALL_HOLD_STEPS 10000` (208). Other release sites: 18373 (speed of sound), 19386 (energy transfer), 22132 (GUI EDMD, only with `cli_auto_release_after_hold`) |
| Piston: one vertical face per side; mass 0 = prescribed velocity, mass > 0 = free | `collide_time_piston_L/R` 934-953; resolve 1412-1436 | 510-529; 992-1013 | `OBJ_PL/OBJ_PR`, `O->h = S->R` (1166-1167) | `handle_piston_collisions` 8541-8572; CCD version 8576-8650 |
| Heat-bath outer walls (MB resampling) | `resolve_wall` 1290-1324 | 900-927 | refused: `if (prm->heatbath_enabled) FAIL(...)` (1119) | 8426-8434 and the same block per wall |
| Semipermeable divider gates (Szilard) | edmd.h:59-93; `resolve_wall` 1356-1385 (pass-through repositions `A->x`, 1375, 1381) | `divider_gate_allows_pass` 833ff, used 942 | refused (1120, 1127) | gate in `update_particles_with_substepping` 9440-9443 |
| Periodic boundaries | none | none | none | particle-life (soft-force) mode only: `pl_apply_cyclic_boundaries` 7672-7705, used 7962-7964 |

Notes:
- **Piston mass in practice:** every GUI setup passes piston mass 0: `double mL = 0.0; double mR = 0.0;` (00ALLINONE.c:18921, 21650), `pc.pistonL_mass = 0.0` / `pc.pistonR_mass = 0.0` (17314-17315). Only `run_gen3_replay` reads piston masses from a file (23368-23369). Time-stepped pistons are prescribed-velocity only (8544-8571).
- **Divider count limits:** engine 32 (edmd.h:13; `MAX_DIVIDER_CAPACITY 32`, 00ALLINONE.c:45). `--num-walls` accepts 1-5 (4589); `--walls` is clamped to 32 (4957-4959) but `cli_num_walls` (default 1, 643) overrides it when > 0 (5422-5424). CLI arrays hold 5: `cli_wall_positions[5]` (644), `cli_wall_masses[5]` (640); `LIVE_SCENE_MAX_WALLS 5` (985).
- **Spring placement:** in every GUI setup the spring goes on `primary_wall_index`, set to `leftmost_wall_index` (3604): `divider_k[primary_idx]` / `divider_k[pidx]` at 21642, 18913, 15398, 14520, 23072. `run_gen3_replay` reads `divider_k[d]` per divider from a file (23366-23367).
  - `--spring-wall` (parsed 4658-4666, accepts 0..4, 4663) is used at 523, 4666, 15211, 19808, 23550, 23560: the empty-compartment check (23546-23562) and the record (19808). Without it the spring constant is zeroed (23571). The check reads `cli_particles_box_counts` only when given (23552-23553).
  - Not found: a line that picks the spring's divider by `cli_spring_wall_index` (searched `cli_spring_wall_index` in 00ALLINONE.c). The help text says "Declare which wall carries the spring (0-based, left to right)" (4000). Unsure whether an index other than the leftmost wall was intended.
- **Bodies against each other or the box:** no divider-divider, divider-piston or divider-outer-wall event type (enums edmd.c:65, edmd_accelerated.c:54, edmd_gen3.c:89).
  - gen2 clamps dividers to [th/2+1e-9, boxW-th/2-1e-9] (edmd.c:984-995). A divider face outside (0, boxW) gets no events (675, 726).
  - gen3 counts a body outside the box or two overlapping slabs as a validator finding, not an event (edmd_gen3.c:648-668).
  - The gen3 window path stops a piston in the driver before it would leave the box or reach a divider (`g3_gui_guard_pistons`, 00ALLINONE.c:11902-11931).
- **Stale header comment:** edmd.h:5 says "No pistons here (yet)"; pistons are at edmd.h:32-41.

**Piston trajectory laws.** `apply_piston_protocol`, 00ALLINONE.c:7056-7159. Only the right piston is driven; the protocol sets its velocity.
- Position: time mode `move_right_piston`, `piston_right_x += vx_piston_right * dt;` (7009), snapped to the target when reached (7011-7016). EDMD paths hand the velocity to the engine as a mass-0 face (19239-19244, 22259-22263); the engine moves it: `S->prm.pistonR_x += S->prm.pistonR_vx * dt;` (edmd.c:1734), gen3 `obj_at` (edmd_gen3.c:304).
- t = `simulation_time - piston_protocol_start_time` (7063); σ-time per comment 571 and 8757-8758. Default `cli_protocol = PROTOCOL_STEP` (645); direction default `piston_step_direction = -1.0f` (374).

| law | lines | v(t) |
|---|---|---|
| STEP | 7068-7071 | `piston_step_direction * vmax` |
| SIGMOIDAL | 7073-7082; profile 7038-7053 | normalised logistic from v0 to vmax over `duration`; k = steepness/duration, or 8/duration if steepness <= 0 (7044) |
| LINEAR | 7084-7097 | with gradient: min(vmax, max(0, v0 + g t)); otherwise vmax t/duration, then vmax |
| SINUSOIDAL | 7099-7106 | vmax sin(π t/duration) for t < duration, then vmax. From the code v → 0 as t → duration, then jumps to vmax |
| RAMP | 7108-7127 | min(vmax, a t, sqrt(2 a d)), a = vmax/t_ramp, t_ramp ≤ travel/vmax (7116-7119), d = remaining distance to target (7121) |
| OPTIMAL | 7129-7138 | same call as SIGMOIDAL ("placeholder" per comment 7130) |

- Stop: v = 0 when `piston_right_x` reaches `piston_step_target` (7142-7146). In the EDMD paths the position is read back from the engine (19228-19229, 22242-22243) before `apply_piston_protocol` (19238, 22257) and the engine then advances a full fixed step (19252, 22270); no snap as in 7011-7016. From this order the face can pass the target by up to |v|·fixed_dt (inferred, not run).
- The left piston has no protocol: `vx_piston_left`, `move_left_piston` 6987-7002; EDMD passes it with mass 0 (21645-21650).

**Non-axis-aligned or curved walls: none. Confirmed from the collision functions.** Every wall predictor tests x alone against an x-face or y alone against a y-face:
- **gen2:** outer walls `return wall_time_from_gap((A->x - S->prm.radius) - 0.0, -A->vx, tcol);` (edmd.c:520), `(S->prm.boxW - S->prm.radius) - A->x` (525), `(A->y - S->prm.radius) - 0.0` (530), `(S->prm.boxH - S->prm.radius) - A->y` (535). Dividers `double num = L - (A->x + S->prm.radius);` (715), `double dist = (A->x - S->prm.radius) - Rf;` (761); spring dividers `const double base = (A->x - xeq) + (Rpart + 0.5*th);` (697) and 744. Pistons `double num = xL - (A->x - S->prm.radius);` (940), `double num = xR - (A->x + S->prm.radius);` (950). Resolver `if(type==EV_WL || type==EV_WR){ A->vx = -A->vx; }` and `if(type==EV_WB || type==EV_WT){ A->vy = -A->vy; }` (1330-1331); dividers and pistons change `A->vx` only (1396, 1407, 1423, 1430).
- **accelerated:** same expressions at edmd_accelerated.c:294, 299, 304, 309, 412, 490, 516, 526; spring branch on x (354, 361, 432, 439). Resolver 930-931; `A->vx = v1` only at 980, 989, 1003, 1009.
- **gen3:** `case 0: if (A->vx >= 0.0) return 0; *gap = ((real)A->cx * S->w - S->R) + xl; ...` and cases 1-3 (edmd_gen3.c:441-444). Bodies: `return (real)body_side(o, c, xp) * (c - xp) - O->h;` (327) with `const real xp = (real)A->cx * S->w + xl;` (482), x only. Resolution flips vx or vy (742, 744); bodies set `A->vx = v1;` (845).
- **time-stepped:** outer walls `if (X[i] - Ri < XW1)` (8423), `X[i] + Ri > XW2` (8442), `Y[i] - Ri < YW1` (8461), `Y[i] + Ri > YW2` (8479). Slab `const float wL0 = wall_x0 - 0.5f * WALL_THICKNESS;` (9027), inside test on `xs` only (9048); hybrid `if ((X[i] + Ri) > L && (X[i] - Ri) < R)` (9931). Pistons 8544, 8559, 8588-8589, 8621-8622.
- **The only non-axis directions in collision code:**
  - pair line of centres: edmd.c:1057, edmd_accelerated.c:701, edmd_gen3.c:718;
  - heat-bath random rotation of an axis-aligned normal ("Rotate normal by random angle"): edmd.c:1106-1115, 1158-1167; edmd_accelerated.c:744-753, 796-805; 00ALLINONE.c:7551-7560, 7570-7579. Normals set axis-aligned at edmd.c:1302-1305, edmd_accelerated.c:906-909, 00ALLINONE.c:8429-8486;
  - gen3 hook normals are axis-aligned constants: `WNX[4] = {-1.0, 1.0, 0.0, 0.0}, WNY[4] = {0.0, 0.0, -1.0, 1.0}` (edmd_gen3.c:750), body `dir` (848).

**Keyword search.**
- Engines: `circle|\barcs?\b|obstacle|polygon|polyline|segment|rotat|curved|curve|incline|slope|periodic|cyclic|wrap`, case-insensitive, over edmd.c, edmd_gen3.c, edmd_accelerated.c, edmd.h, edmd_gen3.h, edmd_accelerated.h. Hits: edmd.c:29, 1783 ("periodic flush" of the heap); edmd_gen3.c:52 ("arc", band-expiry comment); heat-bath "Rotate" comments edmd.c:1110, 1162, edmd_accelerated.c:748, 800.
- Engines, `normal|angle|\bnx\b|\bny\b`: the pair and heat-bath lines above, plus edmd.c:1590-1593 (lattice-seeder jitter variables named `nx, ny`, not normals).
- 00ALLINONE.c, `obstacle|polygon|polyline|curved|incline|fixed_disk|pinned|slanted|diagonal|oblique|tilted`: 2079, 4041, 8496 (Szilard "pinned" particles), 6965 (spring drawn as "polyline zig-zag"). `circle|arc|rotate|angle|normal_x`: drawing (`SDL_RenderFillCircle` 7351), heat bath (7529-7579), particle-life angle bins (1294-1332, 7756-7778), return-map angles (11939-12016), CSV rotation (3458-3461). None is boundary geometry.

**Not found:**
- Obstacles (search above).
- Fixed or pinned disks: no mass or fixed field in `EDMD_Particle` (edmd.h:96-100) or `Disk3` (edmd_gen3.c:100-111).
- Internal horizontal walls: dividers carry only `divider_x` (edmd.h:25), `Obj3` only x (edmd_gen3.c:122). grep `HORIZONTAL|horizontal wall|horiz_wall|y_wall|wall_y` (case-insensitive) in 00ALLINONE.c hits only substrings of `primary_wall_*`.

The Andersen thermostat is bulk resampling, not a wall: `andersen_thermostat_step` (7640), active only with the heat bath (7641, 8533).

### B3

Paths as in B1.

**GUI program: no scene or domain object that could hold a second independent domain.** "Compartment" and "scene" exist with other meanings:
- **Compartment** = slice between dividers in the one box: `EDMD3_MAX_COMPARTMENTS (EDMD_MAX_DIVIDERS + 1)` (edmd_gen3.h:37); `Disk3.comp`, the number of dividers left of the disk at load (edmd_gen3.c:109-110, 1227-1228); the GUI's `segment_*` arrays (00ALLINONE.c:414-423).
- **Scene** = the F3 panel's parameter set: `LiveSceneMode` (1008-1013), `live_scene_*` (1020-1033). It rebuilds the single global state (`live_rebuild_edmd_from_current_state` 15341; `live_scene_rebuild_requested = 1` at 11896 (ESC in the gen3 key handler) and 15325).
- `SzilardBranchSnapshot` (2164-2213; instance 2215) saves/restores a branch by copying the globals, incl. `double X[NUM_PARTICLES]` and other per-particle arrays (2207-2212).
- Backend dispatch is single-instance: comment "Only valid for the global g_edmd instance" (1628); `g3_live(S)` tests `g3_owner == S` (1615).

**The engines are state-struct based** (below), so two independent domains can live in one process as two engine handles, subject to the process-wide statics listed. The GUI's globals cannot.

**00ALLINONE.c file-scope variables**

Method: Python brace-depth scan run from stdin (no file written). Comments and string literals blanked, preprocessor lines skipped (so both branches of every `#if` are counted), depth-0 statements ending in `;` kept, minus `typedef`s and prototypes (a `(` before any `=` without `(*`); brace initialisers kept. Reproduced independently in this session: **552 statements, 580 declarators**. Four names twice: `enable_speed_of_sound_experiments` (163/165) and `wall_thickness_runtime` (273/276) are `#if`/`#else` alternatives; `cli_gen3_exact_box` (323/1462) and `cli_demo_mode` (1062/1387) are tentative definitions. **576 distinct file-scope variables.**

Limits:
- Classification is by name prefix only (e.g. FFT buffers `wall_x_samples` 3064, `wall_fft_magnitude` 3070 are not in the `wall_` group because the rule is first prefix at the start of the name).
- Line 231 holds a stray `^` inside a branch not compiled with `#define MOLECULAR_MODE 1` (140; `#if MOLECULAR_MODE` 189, `#else` 214).
- Function-local statics are not included: 37 non-const `static` lines inside functions (grep `^[[:space:]]+static` minus `const`, minus 273/276). Simulation-relevant: `static float last_position` (7305, energy measurement), `static double prev_vxL = 0.0, prev_vxR = 0.0;` (22255, piston-velocity change detection), particle-life buffers `ax`, `ay`, `alloc_size` (7836-7838). The rest are GUI histogram smoothing (13095-14018), return-map buffers (12064-12065), print-once flags (1638, 2457-2774) and an analysis buffer (16990).

Prefix counts (first matching prefix at the start of the name; sums to 576):

| prefix | count |
|---|---|
| `cli_` | 204 |
| `sz_` / `g_sz` | 53 / 4 |
| `g3_` (incl. `g3_gui_`) | 42 |
| `gui_` | 38 |
| `piston` | 22 |
| `live_` | 22 |
| `pl_` | 20 |
| `wall_` | 14 |
| `segment_` | 10 |
| `dist_` | 9 |
| `demo_` | 8 |
| `gs_` | 5 |
| other | 125 |

Core physical state, curated: **104 variables**.

| group | n | variables (00ALLINONE.c line) |
|---|---|---|
| Particles | 8 | `particles_active` 156; `X` 748, `Y` 749, `X_old` 750, `V_init` 752, `Vx` 753, `Vy` 754, `Radius` 755 |
| Box and units | 11 | `L0_UNITS` 251, `HEIGHT_UNITS` 252, `SIM_WIDTH`/`SIM_HEIGHT` 315, `XW1` 316, `XW2`/`YW1`/`YW2` 317, `particle_scale_runtime` 237, `wall_thickness_runtime` 273/276, `wall_mass_runtime` 292 |
| Clock and integrator | 8 | `simulation_time` 257, `steps_elapsed` 298, `num_steps` 241, `sim_mode` 934, `fixed_dt_runtime` 929, `temperature_runtime` 928, `time_scale_runtime` 930, `g_last_substeps_used` 682 |
| Primary divider | 10 | 242, 258, 295-297, 299, 693-696 |
| Multi-wall | 13 | 400, 402-413 |
| Pistons | 21 | 367-374, 377-380, 572-580 |
| Spring | 1 | `energy_measurement` 2659 |
| Heat bath | 9 | 939-944, 952-954 |
| Engine handles | 13 | `g_edmd` 1354, `g_edmd_is_acc` 1356, `g_edmd3` 1433, `g3_owner`/`g3_prm`/`g3_copy`/`g3_seen`/`g3_n`/`g3_stale`/`g3_handed`/`g3_work0`/`g3_work0_p`/`g3_load_time` 1434-1440 |
| Segments | 10 | 414-423 |

The Szilard (57) and particle-life (20) groups also hold dynamic state for their modes.

**Engines: all three are state-struct based, handle passed by pointer.** Engine files include only their own header and libc (edmd.c:1-6, edmd_accelerated.c:34-39, edmd_gen3.c:66-71).

- **gen2.** `EDMD*  edmd_create(const EDMD_Params* prm);` (edmd.h:105). Dynamic state in `struct EDMD` (edmd.c:128-195), bodies inside `S->prm`. Process-wide mutable statics:
  - `static long g_edmd_avalanche_warning_count` (14): rate-limits a warning only (1680-1682).
  - `static int g_edmd_legacy_resched` (35): rescheduling policy for every `EDMD` in the process (read at 1671, 1766, 1784).
  - `static int g_edmd_contact_audit = -1;` (45): set lazily from `HD_CONTACT_AUDIT` (1750).
  - `static int g_edmd_resched_audit` (58).
  - Trace ring `g_trace_a/_b`, `g_trace_cap`, `g_trace[]`, `g_trace_seq`, `g_trace_head/_count/_fired` (109-113).
  - `static FILE*  g_edmd_evlog`, `static double g_edmd_evlog_tscale` (1264-1265): one event log for all states.
  - Heat bath uses libc `rand()` (1082-1083, 1092, 1107, 1159). Seeding uses a caller-local xorshift state (`xorshift64` 1486; `st` at 1510 and 1614).
  - Function-local statics: only `static const char* wn[4]` (847, 1342).
- **accelerated.** `struct EDMD` (edmd_accelerated.c:92-117). One static, `g_edmd_avalanche_warning_count` (47). Heat bath uses `rand()` (720-721, 730, 745, 797). Refused on main: `refuse_edmd_acc_backend` (00ALLINONE.c:1357-1360), called at 1631.
- **gen3.** `EDMD3* edmd3_create(const EDMD_Params* prm, double cell_px, char* err, size_t errlen);` (edmd_gen3.h:126). All state in `struct EDMD3` (edmd_gen3.c:134-175).
  - Const statics: `WALL_NAME` (178), `U_ROUND` (293), `WNX/WNY` local to `exec_wall` (750).
  - Only mutable statics: the collision hook `static EDMD3* coll_hook_S`, `static EDMD3_CollHook coll_hook_fn`, `static void* coll_hook_ctx` (704-706); set at 707-709; cleared when its state is destroyed (1180). Header: "one state at a time" (edmd_gen3.h:186). The hook fires only for `coll_hook_S == S` (730, 749, 847).
- **Several engine states in one process:**
  - tests/gen3_m1_harness.c:263-264 holds `EDMD3* S3` and `EDMD* S2` at once, advances both (269), destroys both (275).
  - 00ALLINONE.c keeps the gen2 container `g_edmd` and `g_edmd3` alive together (1433-1437, 1549-1594).
  - Not found: two `EDMD3` states alive at the same time. Searched `edmd3_create|edmd3_destroy` in tests/*.c: gen3_checkpoint_test.c creates and destroys them one after another (90-115, 131-140); gen3_m1_harness.c:117/148, 339/341; gen3_m2_harness.c:189/233; gen3_band_edge_test.c:26/51; gen3_band_edge_ties.c:119/134 likewise.

## C. Dynamics core

### C1

**Verdict**

*Engines*
- gen3 is event-driven, with a separate clock per disk; a disk is moved only when one of its events runs.
  - One binary heap holds predicted event times.
  - Stale events are dropped when popped, using collision counters and body epochs.
  - There is no time step inside the engine.
- gen2 is event-driven with one global clock.
  - Its heap is ordered by `t` only.
  - Invalidation uses `coll_count` snapshots plus a divider epoch.
  - Every disk is moved to every event time.
- The accelerated backend (`--edmd-acc`) has gen2's structure plus cell-crossing events.
- None of the three integrates with a time step.

*Driver*
- In the EDMD modes `--fixed-dt` is the interval between the driver's `advance_to` stops, not an integration step.
- But the driver changes body velocities only at those stops:
  - hold and release go by step count;
  - the piston protocol is evaluated every step from the driver clock `simulation_time` (`void apply_piston_protocol(float dt)`, hspist3/00ALLINONE.c:7056; `elapsed = simulation_time - piston_protocol_start_time`, 7063; called at 19238 and 22257).
- With a time-varying protocol (sigmoidal, linear, sinusoidal, ramp: 7073-7126) the piston velocity is piecewise constant per `--fixed-dt` interval. In that case the stride changes the trajectory, so the coupled system is a hybrid in that one respect.
- In MODE_TIME and MODE_RK4, `--fixed-dt` is the integration step.
- `--seed-drift-order` only orders the initial velocity normalisation.

**gen3 main loop**, `double edmd3_advance_to(EDMD3* S, double t_abs)` (hspist3/edmd_core/edmd_gen3.c:1266-1299)

| step | code |
|---|---|
| fatal | `if (!S \|\| S->fatal) return ...` (1267) |
| stop test | `const real tt = t_abs - S->T0;` (1269). If the heap is empty or `S->heap.d[0].t > tt`, then `if (tt > S->now) S->now = tt; break;` (1270). No disk moves. |
| pop | `heap_pop` (1271; 221-225) |
| stale drop | `if (!ev_live(S, &e)) { S->H.ev_stale++; continue; }` (1272) |
| guards | Past event: `if (e.t < S->now) { S->H.past_event++; continue; }` (1273). Stagnation: more than `same_limit` = max(5000, 4N) (1157) live events at one identical time sets `fatal` and pushes the event back (1274-1281). |
| clock | `S->now = e.t;` (1282); hash update (1283) |
| execute | `exec_pair` / `exec_wall` / `exec_cross` / `exec_band` / `exec_obj` (1284-1288) |
| housekeeping | heap compaction (1290); optional schedule audit (1291); full validator every `check_interval` = 24 units (1156; 1292-1295); origin shift once `now >= EDMD3_ORIGIN_SHIFT` = 8192 (1296; 873-881; edmd_gen3.h:35) |

**Queue**
- One binary min-heap: `typedef struct { Ev3* d; long n, cap; } Heap3;` (98).
- `Ev3 {t, type, a, b, ca, cb, pad}` (91-96).
- Key (t, type, a, b), compared by `ev_less` (186-191).
- Type order: `enum { T_CROSS = 0, T_BAND = 1, T_WALL = 2, T_DIV = 3, T_PISTON = 4, T_PAIR = 5 };` (89).
- No calendar queue and no per-particle event lists.
- Initial capacity is 32N+1024 (1145). The heap grows by `realloc` and calls `exit(3)` on OOM (194-199).
- Compaction runs when `heap.n > compact_at` (1290). `compact_at` starts at 64N+4096 (1147); after each compaction it is max(2·live, 64N+4096) (242-243).

**Prediction**
- Pairs
  - Quadratic `pair_rule` (407-419).
  - The separation is the local difference plus an integer times w (426-427).
  - Only the 9 neighbouring cells are searched (`predict_cell` 473-477; `predict_all` 499).
  - There is no time cut-off. A pair that just collided and has not changed velocity since is skipped (`mutual_last`, 420; used at 423).
- Outer walls: gap/speed (`wall_gap` 439-446, `predict_wall` 447-454), only from cells that can reach the wall (`wall_candidate` 265-272; used at 498).
- Cell crossings
  - At most one live CROSS event per disk (`cross_rule` 456-467, `predict_cross` 468-472).
  - None when the disk moves toward the box edge from the last cell column or row (`A->cx + 1 < S->gw`, `A->cx > 0`, 458-461; `return -1`, 462).
- Dividers and pistons
  - `body_rule` (352-402): constant velocity in closed form (356-361).
  - A spring divider uses a bracketed, safeguarded Newton root (`harm_root` 335-348).
  - Only disks in the body's band columns are predicted (500, 545-554).
  - A BAND event fires at band expiry (`band_compute` 535-543; `exec_band` 857-863).
- After a velocity change the disk gets `predict_all` (495-501).
- After a crossing (`exec_cross` 757-795):
  - pairs with the 3 cells that just became neighbours (780-785);
  - newly reachable walls (786-787);
  - newly entered bands (788-792);
  - the next crossing (793);
  - then `local_check` (794).

**Invalidation** (`ev_live` 227-234)
- Disks carry the collision counter `cnt` (Disk3, 105); events store it as ca/cb.
  - `cnt` is bumped at every velocity change: pair 728, wall 747, body 845.
  - It is also bumped at the repairs (601, 617).
  - A crossing does not bump it (`exec_cross` 757-795).
- DIV and PISTON events check the body `epoch` (231-232), which is bumped in `obj_changed` (558).
- BAND events check `band_gen` (228; bumped 541).

**Time**
- Each disk carries a time stamp `tau` and a position relative to its cell corner (Disk3, 101-103).
- `advance()` (276-280) runs only for the disks of the executing event (712, 739, 759, 819).
- `edmd3_particles` (1303-1311)
  - It computes positions at `now` without storing them in the disks (1304-1308).
  - It writes them into `S->out` and runs `full_check` (1309). `full_check` changes only Health counters (621-671).
- All disks are synchronised only at the origin shift (867-881) and in `edmd3_sync_all` (1301).
- Absolute time is `T0 + now` (1312).

**gen2** (hspist3/edmd_core/edmd.c)
- Heap ordered by `t` alone, with no tie key:
  - push: `if(H->data[p].t <= H->data[i].t) break;` (334);
  - pop compares `.t` only (345-346).
- `Event {t, a, b, ca, cb, type}` (82-87). Nine kinds: `typedef enum { EV_AB=0, EV_WL, EV_WR, EV_WB, EV_WT, EV_DL, EV_DR, EV_PL, EV_PR } EvType;` (65).
- Invalidation
  - Default (minimal) mode drops a stale event before the position jump: `if (!g_edmd_legacy_resched && !event_live(S, &e))` (1671-1674).
  - `event_live` (785-795) checks `coll_count` and, for DL/DR events, `div_epoch` (792). `div_epoch` is bumped when a collision changes a divider velocity (1408).
  - After the jump, 1740-1749 checks `coll_count` only. This is the only check in `--legacy-resched` mode, where every divider event flushes the heap (1766-1767).
- Time
  - One global clock `S->t`. Event times are absolute `S->t + t` (e.g. 975).
  - Every event moves all N disks, the dividers and the pistons to the event time (1727-1735), which is O(N) per event.
- After an AB event
  - `grid_build` (371-439) is O(N) and also clamps disks into the box and out of divider slabs (383-404).
  - Then `schedule_for` runs for both disks against all N disks (1001-1010, loop `j = 0..N` at 1006).
- The grid is built but never read for neighbour search: `S->grid[` occurs only at 366, 375 and 435.
- Rescheduling: a piston event triggers a full O(N²) reschedule (1766-1767). A divider event gets the minimal reschedule (1768-1776) unless `--legacy-resched` is set.
- Per-call caps
  - 250000 events (9) and 5000 stagnant events (`|e.t - last_event_t| <= 1e-13`, 1677; cap at 12).
  - Exceeding either free-flies all disks to the target, reschedules everything and increments `forced_advance_count` (1679-1707).
  - Collisions in the rest of that interval are skipped.
- If the next event lies beyond the target, everything free-flies to the target and the event is pushed back (1709-1723).

**accelerated** (hspist3/edmd_core/edmd_accelerated.c)
- It exports the gen2 API under `edmd_acc_*` names by `#define` (4-33).
- Same global-clock structure as gen2, plus:
  - an `EV_CC` cell-crossing event (54; `schedule_cell_crossing` 580-622);
  - pair prediction against 9 cells only (`schedule_ab_neighbors_9cell` 624-642; used in `schedule_for` 645-651).
- A CC event bumps `coll_count` (1199), which invalidates all of that disk's events.
- Every divider and piston event triggers a full reschedule (1206-1207).
- It has no stale drop before the jump and no `div_epoch` (validation at 1179-1185).

**How 00ALLINONE.c drives the engines**

| mode | enum (hspist3/00ALLINONE.c:933) | parse | one driver step |
|---|---|---|---|
| time | MODE_TIME=0 | 5285 | `update_particles_with_substepping` (9202-9203). Substeps: `SUBSTEPS` or override (9209-9215). Overlap-based pair handling (9490-9541) or `--pp-toi` (9542+). CCD divider `ccd_wall_step` (9399-9426), pistons `handle_piston_collisions_ccd_substep` (9430), outer walls `handle_boundary_collision` (9431). In the scene, `pl_enabled` replaces this with particle-life physics (22366-22370). |
| rk4 | MODE_RK4=1 | 5286 | No RK4 integrator found: `grep -n -i "rk4\|runge"` gives labels and comments only (651, 933, 2278, 4067, 5286, 5292, 9816, 11880, 18566, 19333, 19711, 22832, 22853). Same branch as time (scene 22365-22378; energy transfer 19332-19338). |
| edmd | MODE_EDMD=2 | 5287 | `edmd_backend_advance_to(t0 + dt)` (macro at 2048): gen3 `edmd3_advance_to` (1726-1739), otherwise gen2 or accelerated (1741) |
| edmd-hybrid / hybrid / event-split / event_split | MODE_EDMD_HYBRID=3 | 5288-5291 | gen2 `edmd_advance_pp_only_to` (pairs only), `hybrid_faces_correction_pass` (9916; it calls `handle_piston_collisions` 9923 and `handle_boundary_collision` 10002), write-back through a cast of the const pointer, then pair reschedule (scene 22329-22357; energy transfer 19302-19321). gen3 refuses it (1744; at start-up 23601). |

**Engine selection and mode overrides**
- `cli_engine` defaults to 2 (1431).
- `--engine=gen3` sets 3 (5214-5215); `--engine-replay=` also sets 3 (5203).
- The F3 panel can switch the engine at a scene rebuild (22038).
- gen3 requires `--mode=edmd` (23601).
- gen3 is built at the first advance from the gen2 container (1726-1727; `g3_build` 1549-1594). After that the gen3 branch returns before the gen2 call (1739), so the container is not advanced again.
- The driver overrides the mode in places:
  - MODE_EDMD with more than 32 dividers falls back to hybrid (scene 21584-21588; energy transfer 18858-18861);
  - speed of sound switches hybrid to EDMD (18003-18006);
  - simple box forces EDMD (20655-20657);
  - the Szilard experiment forces TIME (21143-21145);
  - Szilard interactive in hybrid falls back to TIME (10666-10669).

**`--fixed-dt`**
- Parsed at 4524-4529; `strncmp(arg, "--dt", 4)` also accepts the prefix `--dt`.
- Stored as `static float cli_fixed_dt` (1399) and copied into `static float fixed_dt_runtime` (929) by `update_dt_runtime_for_temperature` (2329-2331).
- Without the flag:
  - `FIXED_DT 0.04f` (201), since MOLECULAR_MODE=1 and EXPERIMENT_MODE=0 (140-141, 189-201);
  - `FIXED_DT / sqrtf(T)` for the speed-of-sound preset without `--kbt1` (2334-2336).
- In the window, the live slider "fixed dt" (1037; range 0.001..0.080) also writes it (`live_slider_set_value` 14475-14490). This is not blocked under gen3; only the HB slider is (14478).
- In the EDMD modes it is used as `edmd_backend_advance_to(g_edmd, t0 + (double)fixed_dt_runtime)`. Call sites: 18170 and 18220 directly; 19252 and 22270 through the macro at 2048.
- It is also the cadence of:
  - the piston protocol (19238, 22257);
  - the readers;
  - the validator (18176, 18227-18229).
- Units: `simulation_time += fixed_dt_runtime / PIXELS_PER_SIGMA;` (22435) with `PIXELS_PER_SIGMA 24.0f` (202). So 24 internal units = 1 σ-time, as also stated at edmd_gen3.h:22-23.
- Campaign value: `--fixed-dt=0.4` (hspist3/cluster/gen3_stage1_261012/s1_worker_koa.sh:85; hspist3/cluster/resched_gate_261005/audit_runs_261007.py:56). That is 1/60 σ-time per the comment at 1454.
- Stop times are `t0 + (double)fixed_dt_runtime` with `t0` read back from the engine (18169-18170, 19251-19252, 22269-22270). Each target is the previous engine time plus (double)0.4f, added in double.
- Recorded evidence: `t_end=1200.0000178813934` after 3000 steps (hspist3/experiments_gen3_d19_261010/evidence_7f6adca/steer/gen3_window/run.log:72; `steps=3000` at hspist3/experiments_gen3_d19_261010/D19_RUN_RECORDS.txt:110).

**Does the stride change the trajectory?**
- gen3
  - A stop only sets `now` (1270); the readers change only counters (1309; 621-671). Per the comment at 37-39, readers and stops do not steer.
  - Body velocity changes from the driver happen at stops, so with a time-varying piston protocol the stride matters (see Verdict).
  - Recorded replay test, all four variants (`--replay-read-every=1/0/1/7` with `--replay-substeps=1/1/3/2`; experiments_gen3_m3_261009/run_replays.sh:17-19): 64 of 64 MATCH with clean=1 (stageA_tables_output.txt:140-147). Build f42befb (stageA_tables_output.txt:14); not rerun here.
- gen2
  - Every stop free-flies all disks (1711-1719).
  - The event caps count per call (1652-1653, 1679).
  - So rounding and forced-advance exposure depend on the stride. How much is unsure; not tested here.

**`--seed-drift-order`**
- Parsed at 4338-4342.
- Used only in `equalize_temperature_per_segment` (2550-2586; the only reads are 2566 and 2584): the compartment drift is removed before (drift-first) or after (old, the default) the rescale.

**`--edmd-acc`**
- Sets `g_edmd_is_acc` (1632), which routes the backend calls to `edmd_acc_*` (e.g. 1650, 1676, 1741).
- Guard on HEAD: `static void refuse_edmd_acc_backend(void)` calls `exit(2)` unless the build defines `-DHD_ALLOW_ACC_BACKEND` (1360-1367).
- It is called for:
  - bare `--edmd-acc` (5266-5268);
  - an empty value (5271-5273);
  - `1/true/yes/on` (5274-5276);
  - again in `edmd_backend_create` (1631).
- `--edmd-acc=0/false/no/off` is accepted (5277-5278).
- `edmd_accelerated.c` is still compiled and linked (`EDMD_SRCS ?= edmd_core/edmd.c edmd_core/edmd_accelerated.c edmd_core/edmd_gen3.c`, hspist3/Makefile:12; link rule 119-120).

### C2

In the tables, "static" means file-static. Opening braces are omitted from the signatures.

**gen3** (hspist3/edmd_core/edmd_gen3.c)

| role | verbatim signature | line | static |
|---|---|---|---|
| pair predict | `static inline int pair_rule(real rx, real ry, real vx, real vy, real d2, real* dt, real* c_out)` | 407 | yes |
| pair schedule | `static void predict_pair(EDMD3* S, int i, int j)` | 422 | yes |
| pair resolve | `static void exec_pair(EDMD3* S, const Ev3* e)` (equal unit mass, normal exchange 719-723) | 710 | yes |
| outer wall reach | `static inline int wall_candidate(const EDMD3* S, int cx, int cy, int s)` | 265 | yes |
| outer wall predict | `static inline int wall_gap(const EDMD3* S, const Disk3* A, real xl, real yl, int s, real* gap, real* speed)` / `static void predict_wall(EDMD3* S, int i, int s)` | 439 / 447 | yes |
| outer wall resolve | `static void exec_wall(EDMD3* S, const Ev3* e)` (specular: `A->vx = -A->vx` or `A->vy = -A->vy`, 742-745) | 737 | yes |
| cell crossing | `static inline int cross_rule(const EDMD3* S, const Disk3* A, real x, real y, real* dt)` / `static void predict_cross(EDMD3* S, int i)` / `static void exec_cross(EDMD3* S, const Ev3* e)` | 456 / 468 / 757 | yes |
| divider and piston predict | `static int body_rule(const Obj3* O, int o, real c, real vo, real xp, real vp, int mlast, real* dt, real* g0)` / `static void predict_obj(EDMD3* S, int i, int o, real xl)` | 352 / 479 | yes |
| spring contact root | `static real harm_root(const Harm* H, real lo, real hi)` | 335 | yes |
| divider and piston resolve | `static void exec_obj(EDMD3* S, const Ev3* e)` (finite mass: 1D elastic, 824-829; mass 0: `v1 = 2.0 * u2 - u1`, 831) | 817 | yes |
| band | `static void band_compute(EDMD3* S, int o)` / `static void exec_band(EDMD3* S, const Ev3* e)` / `static void obj_changed(EDMD3* S, int o, int skip)` | 535 / 857 / 556 | yes |
| spring update | `static inline void obj_at(const Obj3* O, real t, real* x, real* v)` (closed-form harmonic, 300-303) / `static void obj_advance(EDMD3* S, int o)` (books the anchor impulse, 312-317) | 297 / 308 | yes |
| API body change | `int edmd3_set_divider_motion(EDMD3* S, int d, double mass, double vx)` / `int edmd3_set_piston_motion(EDMD3* S, int side, double mass, double vx)`, both through `static int set_motion(EDMD3* S, int o, real mass, real vx)` (1485) | 1505 / 1506 | no |
| repairs (safety nets) | `static void local_check(EDMD3* S, int i)` / `static void reflect_into_box(EDMD3* S, int i)` | 567 / 609 | yes |
| heat bath | none. Refused: `if (prm->heatbath_enabled) FAIL("edmd_gen3: the outer-wall heat bath is not implemented");` | 1119 | – |
| divider–divider, divider–piston | not found. No event type exists (89); overlapping slabs and piston order are only validator findings (`body_findings`, 662-667) | – | – |
| divider–outer wall | not found as an event. The validator counts a divider outside the box as `body_findings` (660). In the window the scene pauses before a divider reaches a wall (`gui_divider_limit_check`, 00ALLINONE.c:21384; message "has no divider-wall event", 21410-21412). | – | – |

**gen2** (hspist3/edmd_core/edmd.c; everything is static unless noted)

| role | verbatim signature | line |
|---|---|---|
| pair predict | `static int collide_time_ab(double xi,double yi,double vxi,double vyi,` / `double xj,double yj,double vxj,double vyj,` / `double R, double* tcol)` (3 lines) | 474-476 |
| pair schedule | `static void schedule_ab(EDMD* S, int i, int j)` | 962 |
| pair resolve | `static void resolve_ab(EDMD* S, int i, int j)` | 1051 |
| outer wall predict | `static int wall_time_from_gap(double gap, double speed, double* tcol)`; `static int collide_time_wall_L(const EDMD* S, const EDMD_Particle* A, double* tcol)` (R 523, B 528, T 533); `static void schedule_walls(EDMD* S, int i)` | 512; 518; 644 |
| divider predict | `static int collide_time_divider_L(const EDMD* S, const EDMD_Particle* A, int d, double* tcol)` / `static int collide_time_divider_R(const EDMD* S, const EDMD_Particle* A, int d, double* tcol)`; `static void schedule_divider_one(EDMD* S, int i, int d)` | 670 / 721; 770 |
| spring contact | `static int harmonic_first_contact(double base, double particle_v,` / `double wall_dx, double wall_v0,` / `double omega, double t_max,` / `int crossing_direction, double* tcol)` (4 lines) | 551-554 |
| spring update | `static inline void harmonic_advance_1d(double k, double m, double x_eq, double dt, double* x_inout, double* v_inout)`; `static void advance_dividers(EDMD* S, double dt)` (clamps the divider at the box, 984-995; this is the only divider–outer wall handling) | 456; 979 |
| piston predict | `static int collide_time_piston_L(const EDMD* S, const EDMD_Particle* A, double* tcol)` / `static int collide_time_piston_R(const EDMD* S, const EDMD_Particle* A, double* tcol)`; `static void schedule_pistons(EDMD* S, int i)` | 934 / 944; 954 |
| wall, divider, piston resolve | `static void resolve_wall(EDMD* S, int i, EvType type, int b)` (walls 1290-1350, dividers 1351-1411, pistons 1412-1436) | 1286 |
| heat bath | `static inline void thermal_wall_bounce_edmd(EDMD_Particle* A,` (1100-1102); `static inline void gradual_damping_bounce_edmd(EDMD_Particle* A,` (1119-1121); `static inline void adaptive_thermal_wall_bounce_edmd(EDMD_Particle* A,` (1139-1141); each continues on the next lines. Called from `resolve_wall` (1291-1323). They draw libc `rand()` (1082, 1083, 1092, 1107, 1159). | 1100-1141 |
| semipermeable gate | `static inline int divider_gate_allows_pass(const EDMD* S, int d, EvType type, int species, int i)` | 1195 |
| repair (public) | `void edmd_divider_resolve_overlaps(EDMD* S)` (not static) | 1942 |
| divider–divider | not found: EvType (65) has no such kind, and `advance_dividers` limits only at the box (984-995) | – |

**accelerated** (hspist3/edmd_core/edmd_accelerated.c; all static; first lines verbatim)
- `static int collide_time_ab(double xi,double yi,double vxi,double vyi,` (257)
- `static void schedule_walls(EDMD* S, int i){` (314)
- `static int collide_time_divider_L(const EDMD* S, const EDMD_Particle* A, int d, double* tcol){` (340); `_R` (418)
- `static int collide_time_piston_L(const EDMD* S, const EDMD_Particle* A, double* tcol){` (510); `_R` (520)
- `static void schedule_cell_crossing(EDMD* S, int i){` (580)
- `static void schedule_ab_neighbors_9cell(EDMD* S, int i){` (624)
- `static void resolve_ab(EDMD* S, int i, int j){` (695)
- `static void resolve_wall(EDMD* S, int i, EvType type, int b){` (890)
- heat bath: 738, 757, 777

**00ALLINONE.c, time-stepped modes** (hspist3/00ALLINONE.c)

| role | verbatim signature | line | static |
|---|---|---|---|
| step driver | `void update_particles_with_substepping(float dt, int* out_left, int* out_right,` / `double current_sim_time, bool is_wall_released, double wall_released_at)` | 9202-9203 | no |
| pair–pair | inline in the function above: float overlap test, normal velocity exchange, push-apart (9507-9535); or the `--pp-toi` branch (9542+). No separate function. | 9490-9541 | – |
| outer walls (and heat bath) | `void handle_boundary_collision(int i)` (heat-bath calls 8426-8488; called at 9431) | 8414 | no |
| divider | `static void ccd_wall_step(int i,           // particle index` (8 lines with inline comments, ending `bool force_locked)`; called at 9399-9426) | 8997-9004 | yes |
| piston | `void handle_piston_collisions(int i, float sub_dt)` (called from `hybrid_faces_correction_pass`, 9923); `static inline void handle_piston_collisions_ccd_substep(int i, float sub_dt,` (3 lines; called at 9430) | 8541; 8576-8578 | no; yes |
| piston motion | `void move_left_piston(float dt)` / `void move_right_piston(float dt)` / `void apply_piston_protocol(float dt)` | 6987 / 7005 / 7056 | no |
| divider motion | `void update_wall(float dt, float L0)` / `static void update_extra_walls(float dt, float L0)` | 8737 / 8819 | no / yes |
| spring | inline explicit update of the wall velocity: `vx_wall += accel * sub_dt;` | 9319-9327 | – |
| heat bath | `static inline float sample_MB_speed_2D(float temperature)`; `static inline void adaptive_thermal_wall_bounce(double *vx, double *vy,` (3 lines); `static inline void thermal_wall_bounce(double *vx, double *vy,` (3 lines); `static inline void gradual_damping_bounce(double *vx, double *vy,` (2 lines); `static inline void andersen_thermostat_step(double dt)` (gen3 refuses the heat bath and Andersen in the scene, 23609-23610) | 7519; 7528-7530; 7564-7566; 7583-7584; 7640 | yes |
| hybrid faces | `static void hybrid_faces_correction_pass(float dt)` | 9916 | yes |
| divider–divider | not found. Searched `wall-wall`, `wall_wall`, `divider-divider`, `divider_divider`, `walls collide`, `wall pair`, `adjacent wall`, `wall_sep` in 00ALLINONE.c and edmd_core/*.c. `update_wall` and `update_extra_walls` (8737-8864) clamp only at the box limits (8794-8812, 8847-8855). | – | – |

### C3

**Who owns the clock**
- The engine owns simulation time:
  - gen3: `double edmd3_time(const EDMD3* S){ return S->T0 + S->now; }` (edmd_gen3.c:1312);
  - gen2: `edmd_time` returns `S->t` (edmd.c:1647).
- The driver reads the engine time and asks for `t0 + (double)fixed_dt_runtime` (hspist3/00ALLINONE.c:18169-18170, 18219-18220, 19251-19252, 22269-22270).
- The engine processes every event up to that target (C1). The exception is gen2's per-call caps, which free-fly the rest of the interval (edmd.c:1679-1707).
- A second clock lives in the driver: `double simulation_time` (257).
  - It is incremented by `fixed_dt_runtime / PIXELS_PER_SIGMA` each step (scene 22435, energy transfer 19351).
  - The piston protocol reads it (7059-7063).
  - The speed-of-sound loop takes its times from the engine instead (18186, 18206).

**The GUI frame loop** (`void simulation_loop()`, 21567; `while (running)` 21999-22987)
1. Window: `accumulator += frame_time * time_scale_runtime;`. `frame_time` is in seconds and capped at 0.1 (22002-22006); `accumulator` is a float (21987). Then SDL events and keys are handled (22010-22025).
2. `while (accumulator >= fixed_dt_runtime)` runs whole physics steps (22090-22899). The accumulator is decremented after each pass, also when paused (22898). Each EDMD step ends by copying all particles from the engine into the driver arrays (22285-22291).
3. One render per frame, window only (22903-22986): `SDL_RenderPresent` (22971), then `SDL_Delay(5)` (22985).

Time scale and step size:
- The GUI decides how many steps of size `fixed_dt_runtime` run per wall second.
- In the window the step size can also be changed:
  - the live slider "fixed dt" (1037, range 0.001..0.080) writes `fixed_dt_runtime` (14475-14490), also under gen3;
  - the keys `[`, `]`, `\` call `update_dt_runtime_for_temperature` (12236-12256). That changes dt only for the speed-of-sound preset without `--fixed-dt`/`--kbt1` (2325-2338), and gen3 refuses these keys (11892-11893; filter at 12126).
- `--gui-time-scale=X`
  - Parsed at 5222-5225 (range 0.01..10000).
  - Sets `time_scale_runtime` (a float, 930) once at the start of the scene (21996).
  - Refused outside the scene (23584-23586).
- Since the accumulator compares wall seconds × X with `fixed_dt` in internal units, X is internal units per wall second (X/24 σ-time per second). That holds only while frames take under 0.1 s (22005) and the compute keeps up.
- Changing the scale in the window:
  - keys: `-` halves with floor 0.1 (12366), `=`/`+` doubles to at most 10000 (12373), `0` resets to 1 (12379);
  - the live slider "time scale" covers 0.10..256 (1036).
- `--timescale` (5116-5123) is shadowed by the earlier `--time` prefix test (4581), per comment at 1501-1502.
- `--demo` adds `(demo_steps_per_frame - 1)` extra steps per frame (22088-22089).
- Headless scene: `accumulator += fixed_dt_runtime;` per loop iteration (22027-22028), i.e. one step per iteration, no rendering.

**Scene in the window vs the headless scene.** Both run `simulation_loop`, and the step code 22090-22899 is shared. Differences found:
- events, keys and live sliders: window only (22010-22025);
- `gui_divider_limit_check` returns 0 when headless (21384) and can pause the window (22094);
- `g3_gui_keys_pump`: window only (22438);
- `rm_frame()`, the display-only collision hook: window only (22908; hook edmd_gen3.c:700-709);
- failure handling (`g3_soft_window = !cli_headless;`, 21997): a refused gen3 build or a stopped engine pauses the window; headless calls `exit(2)` (1527-1536, 1553, 1561-1562, 1736-1737);
- `if (!running && cli_engine == 3) break;` (22091). Under gen2 a frame finishes its remaining steps after a stop.
- `--gui-run-seconds` ends the scene by wall-clock time (22031), so its step count is not reproducible.

Recorded evidence (build 7f6adca, not rerun):
- gen3 window and headless gave the same hash `daea5de1378282dc`, with 4 of 4 output files equal (hspist3/experiments_gen3_d19_261010/d19_tables_output.txt:33).
- gen2: wall_position.csv in the window had 30 more rows (d19_tables_output.txt:34). `[GUI-STATS]` shows 3030 steps in the window against 3000 headless (D19_RUN_RECORDS.txt:97, 100).

**Experiments vs the scene.** These are different functions.
- main dispatches single-test, simple-box, energy-transfer, Szilard and speed-of-sound runs instead of `simulation_loop` (23670-23682).
- The experiment loops do not render: `SDL_RenderPresent` occurs only at 20641 (inside `simple_box_render_preview`, 20553) and at 22971.
- Under gen3 the experiments must run headless: `"with the window gen3 runs the interactive scene only (--show-simulation); the experiments need --headless"` (23599-23600).

Per-step EDMD calls, scene vs energy-transfer experiment:

| call | scene `simulation_loop` | energy transfer `run_energy_transfer_experiment` (18747) |
|---|---|---|
| `edmd_config_dividers` | every step, gen2 only (22251) | not called in the step |
| `edmd_set_divider_motions` | every step (22252); under gen3 it reaches the engine only on a real change (1699-1708) | only if changed (19222-19236) |
| piston protocol | `apply_piston_protocol` (22257) plus `g3_gui_guard_pistons` (22258; 11902-11930), which can set a piston velocity to 0 under gen3 (11928-11929) | `apply_piston_protocol` (19238) |
| `edmd_config_pistons` | if the velocity differs from a function-static previous value (22255-22265) | if it differs from the engine's (19239-19246) |
| `edmd_reschedule_all` | **every step** (22267). gen2: heap cleared and all pairs rescheduled (edmd.c:1883, 1021-1035). gen3: no-op apart from the reader-write check (1675). | only on a configuration change (19247-19249) |
| advance | 22270 | 19252 |
| particle copy (reader) | every step (22285-22291); under gen3 each read runs `full_check` (1671, 1601, edmd_gen3.c:1309) | under gen3 only on the steps that need it (`et_particles_due`, 19258-19267) |
| `edmd_divider_resolve_overlaps` | after every step, gen2 only and not with Szilard gates (22275) | none |
| Szilard repairs | called (22278-22283); no-ops outside the Szilard preset (1891-1893, 1979-1980, 10849-10852) | none |
| Andersen thermostat | if enabled, writes velocities into the engine array (22385-22396); gen3 refuses it in the scene (23609-23610) | not in the EDMD step |

What this means:
- Under gen3 both paths reach the same engine API (`edmd3_set_*_motion`, `edmd3_advance_to`). The driver protocol differs (piston guard, read cadence). The read cadence changes only the `full_checks` counter, plus any findings those extra checks make.
- Under gen2 the scene adds a full reschedule and an overlap repair every step. Whether that changes trajectories bitwise against the experiment loop is unsure: no scene-vs-experiment test was found (d19 compared the window and headless scenes only).
- The scene also writes files itself (`wall_log`, 22862-22886; files closed at 23002-23006).

### C4

**Boundary**
- gen3 has a library boundary on state the caller owns, through the opaque handle `typedef struct EDMD3 EDMD3;` (hspist3/edmd_core/edmd_gen3.h:39). The core call is `double edmd3_advance_to(EDMD3* S, double t_abs);` (133).
- gen2 has `double edmd_advance_to(EDMD* S, double t_target);` (edmd.h:146) on `typedef struct EDMD EDMD;` (edmd.h:102).
- The driver's own physics is not callable as a library. The time, hybrid-face and experiment loops are functions of 00ALLINONE.c that work on file-scope globals:
  - `double X[NUM_PARTICLES];`, `Y`, `X_old`, `Vx`, `Vy`, `Radius` (hspist3/00ALLINONE.c:748-755);
  - `static EDMD* g_edmd` (1354) and `static EDMD3* g_edmd3` (1433).
- The driver's dispatch layer is global too: `// EDMD backend dispatch (default vs accelerated). Only valid for the global g_edmd instance.` (1628).

**Public API of edmd_gen3.h** (hspist3/edmd_core/edmd_gen3.h; trailing comments omitted)

Types: Health 41-67, Audit 72-93, Ledger 101-109, Tol 112-122, Coll 184.

```
126 EDMD3* edmd3_create(const EDMD_Params* prm, double cell_px, char* err, size_t errlen);
127 void   edmd3_destroy(EDMD3* S);
130 int    edmd3_load(EDMD3* S, const EDMD_Particle* P, double t_abs, char* err, size_t errlen);
133 double edmd3_advance_to(EDMD3* S, double t_abs);
134 void   edmd3_sync_all(EDMD3* S);
135 const EDMD_Particle* edmd3_particles(EDMD3* S);
136 double edmd3_time(const EDMD3* S);
137 int    edmd3_count(const EDMD3* S);
138 double edmd3_cell_px(const EDMD3* S);
139 int    edmd3_fatal(const EDMD3* S, const char** msg);
143 int    edmd3_set_divider_motion(EDMD3* S, int d, double mass, double vx);
144 int    edmd3_set_piston_motion(EDMD3* S, int side, double mass, double vx);
145 int    edmd3_divider_state(const EDMD3* S, int d, double* x, double* vx);
146 int    edmd3_piston_state(const EDMD3* S, int side, double* x, double* vx);
147 double edmd3_divider_impulse(const EDMD3* S, int d, int face);
148 long   edmd3_divider_events(const EDMD3* S, int d, int face);
149 double edmd3_piston_impulse(const EDMD3* S, int side);
150 long   edmd3_piston_events(const EDMD3* S, int side);
151 double edmd3_work_divider(const EDMD3* S, int d);
152 double edmd3_work_piston(const EDMD3* S, int side);
153 void   edmd3_ledger(EDMD3* S, EDMD3_Ledger* L);
154 int    edmd3_health_clean(const EDMD3* S);
155 void   edmd3_tolerances(const EDMD3* S, EDMD3_Tol* T);
158 const EDMD3_Health* edmd3_health(const EDMD3* S);
159 uint64_t edmd3_event_hash(const EDMD3* S);
160 double edmd3_kinetic_energy(EDMD3* S);
161 void   edmd3_reset_virial(EDMD3* S);
162 double edmd3_compressibility_Z(EDMD3* S);
167 int    edmd3_virial_compartments(const EDMD3* S, int max, double* W, double* KE, long* npair, int* ndisk, double* t_window);
168 long   edmd3_virial_pair_events(const EDMD3* S);
169 double edmd3_virial_sum(const EDMD3* S);
170 double edmd3_wall_impulse(const EDMD3* S, int wall);
171 long   edmd3_wall_events(const EDMD3* S, int wall);
172 void   edmd3_set_check_interval(EDMD3* S, double units);
175 void   edmd3_set_event_log(EDMD3* S, FILE* f, double time_scale);
185 typedef void (*EDMD3_CollHook)(void* ctx, const EDMD3_Coll* c);
186 void   edmd3_set_coll_hook(EDMD3* S, EDMD3_CollHook fn, void* ctx);
189 void   edmd3_set_contact_audit(EDMD3* S, int on);
190 long   edmd3_contact_audit_stats(const EDMD3* S, double max_gap_px[2]);
192 long   edmd3_tie_stats(const EDMD3* S, long* n_live, long* n_now);
197 int    edmd3_checkpoint_write(const EDMD3* S, FILE* f);
198 int    edmd3_checkpoint_read(EDMD3* S, FILE* f, char* err, size_t errlen);
199 long   edmd3_contact_audit_stats4(const EDMD3* S, double max_gap_px[4]);
200 void   edmd3_set_schedule_audit(EDMD3* S, long every);
201 void   edmd3_set_schedule_audit_bodies(EDMD3* S, int on);
202 void   edmd3_schedule_audit_now(EDMD3* S);
203 const EDMD3_Audit* edmd3_schedule_audit_stats(const EDMD3* S);
```

What the gen3 API fixes:
- Input is `EDMD_Params` (edmd.h:15-94) and `EDMD_Particle {x, y, vx, vy, coll_count}` (edmd.h:96-100).
- Geometry is only `boxW`, `boxH`, one `radius` (`/* disk radius (all equal) */`, edmd.h:17), up to `EDMD_MAX_DIVIDERS 32` vertical dividers (edmd.h:13, 23-31) and 2 pistons (33-41). There is no field for any other boundary.
- `edmd3_create` refuses (1114-1130):
  - N <= 0, a radius <= 0, a cell narrower than a diameter, a non-integer cell width, a box smaller than a disk;
  - the heat bath, species, switched-off pairs;
  - a divider of thickness <= 0, gates, bad body values.
- `edmd3_load` refuses an overlapping input or a disk outside the box (1238-1244).
- Units are px and internal time; disk mass is unit (edmd_gen3.h:22-23).

**Public API of edmd.h** (hspist3/edmd_core/edmd.h)
```
105 EDMD*  edmd_create(const EDMD_Params* prm);
106 void   edmd_destroy(EDMD* S);
107 void   edmd_init_random_gas(EDMD* S, unsigned long long seed);
113 int    edmd_init_lattice_gas(EDMD* S, unsigned long long seed);
119 void   edmd_set_event_log(const char* path, double time_scale);
120 void   edmd_close_event_log(void);
125 void   edmd_set_legacy_resched(int on);
126 int    edmd_legacy_resched(void);
128 long   edmd_past_event_count(const EDMD* S);
131 long   edmd_contact_audit_stats(const EDMD* S, double max_gap_px[4]);
135 void   edmd_set_resched_audit(int mode);
136 int    edmd_resched_audit_mode(void);
137 void   edmd_resched_audit_stats(const EDMD* S, long out[7], double* maxdt, double* maxrel);
140 double                edmd_time(const EDMD* S);
141 const EDMD_Params*    edmd_params(const EDMD* S);
142 const EDMD_Particle*  edmd_particles(const EDMD* S);
143 int                   edmd_count(const EDMD* S);
146 double edmd_advance_to(EDMD* S, double t_target);
152 double edmd_advance_pp_only_to(EDMD* S, double t_target);
155 void   edmd_reschedule_all(EDMD* S);
159 void   edmd_reschedule_all_pp_only(EDMD* S);
162 void   edmd_config_dividers(EDMD* S, int count, const double* cx, const double* thickness);
165 void   edmd_set_divider_motions(EDMD* S, int count, const double* mass, const double* vx);
168 void   edmd_set_divider_springs(EDMD* S, int count, const double* k, const double* xeq);
171 void   edmd_config_divider(EDMD* S, int enabled, double cx, double thickness);
172 void   edmd_set_divider_motion(EDMD* S, double mass, double vx);
173 void   edmd_set_divider_spring(EDMD* S, double k, double xeq);
176 void   edmd_config_pistons(EDMD* S,  (3 lines, 176-178: int hasL, double xL, double vxL, double mL, / int hasR, double xR, double vxR, double mR);)
181 double edmd_total_kinetic_energy(const EDMD* S, double mass);
184 double edmd_gas_temperature(const EDMD* S, double mass, double kB);
188 double edmd_work_divider(const EDMD* S);
190 double edmd_work_divider_i(const EDMD* S, int divider_index);
191 double edmd_work_pistonL(const EDMD* S);
192 double edmd_work_pistonR(const EDMD* S);
195 double edmd_heat_bath(const EDMD* S);
196 void   edmd_reset_work(EDMD* S);
202 long   edmd_forced_advance_count(const EDMD* S);
210 void   edmd_reset_virial(EDMD* S);
211 double edmd_compressibility_Z(const EDMD* S);
212 double edmd_virial_accum(const EDMD* S);
218 double edmd_wall_impulse(const EDMD* S, int wall);
219 long   edmd_wall_events(const EDMD* S, int wall);
220 double edmd_wall_Z_x(const EDMD* S);
221 double edmd_wall_Z_y(const EDMD* S);
222 long   edmd_virial_pair_events(const EDMD* S);
223 double edmd_virial_window(const EDMD* S);
225 long   edmd_clamp_repair_count(const EDMD* S);
226 long   edmd_overlap_repair_count(const EDMD* S);
227 long   edmd_wall_overdue_count(const EDMD* S);
232 void   edmd_debug_set_watch(int a, int b, int history);
235 void   edmd_divider_resolve_overlaps(EDMD* S);
```
- gen2 has no load call. `edmd_particles` returns a const pointer to the engine's own array (edmd.c:1645), and callers cast it and write into it: the driver (00ALLINONE.c:21655-21667, 22349-22356) and the M1 harness (gen3_m1_harness.c:158), followed by `edmd_reschedule_all`.
- The accelerated header (edmd_accelerated.h:10-51) exports the same calls as `edmd_acc_*`, without the lattice seeder, the springs, the event log or the audits.

**Standalone harnesses** (hspist3/edmd_core/tests/)

| file | #include lines | recorded build (comment) |
|---|---|---|
| gen3_m1_harness.c | `#include "../edmd.h"`, `#include "../edmd_gen3.h"`, `#include "gen3_state_io.h"` plus libc (20-28) | `cc -std=c11 -O3 -ffp-contract=off -Wall -Wextra -o <out> edmd_core/tests/gen3_m1_harness.c` / `edmd_core/edmd_gen3.c edmd_core/edmd.c -lm` (15-16) |
| gen3_m2_harness.c | the same three (44-46); `#define _POSIX_C_SOURCE 200809L` (37) | the same with gen3_m2_harness.c (34-35) |
| gen3_checkpoint_test.c | `#define main gen3_m2_harness_main` (21), `#include "gen3_m2_harness.c"` (22) | `cc -std=c11 -O3 -ffp-contract=off -Wall -Wextra -o <out> edmd_core/tests/gen3_checkpoint_test.c` / `edmd_core/edmd_gen3.c edmd_core/edmd.c -lm` (18-19) |
| gen3_band_edge_test.c | `#include "../edmd.h"`, `#include "../edmd_gen3.h"`, `#include "gen3_band_edge_cells.h"` (18-20) | `cc -std=c11 -O2 -ffp-contract=off -Wall -Wextra -o <out> edmd_core/tests/gen3_band_edge_test.c edmd_core/edmd_gen3.c -lm` (13) |
| gen3_band_edge_ties.c | `#define EDMD3_AUDIT_HOOK(...)` (31), `#include "../edmd_gen3.c"` (32): white box | `cc -std=c11 -O2 -ffp-contract=off -Wall -Wextra -o <out> edmd_core/tests/gen3_band_edge_ties.c -lm` (24) |
| gen3_body_rule_test.c | `#include "../edmd_gen3.c"` (19) | `cc -std=c11 -O2 -ffp-contract=off -Wall -Wextra -o <out> edmd_core/tests/gen3_body_rule_test.c -lm` (7) |

What the harnesses show:
- The M1 harness drives gen3 through `edmd3_create`, `edmd3_load`, `edmd3_advance_to`, `edmd3_health`, `edmd3_event_hash` and `edmd3_particles` (117-148).
- It also drives gen2 standalone (`edmd_create`, cast `edmd_particles`, `edmd_reschedule_all`, `edmd_advance_to`; 157-195). So gen2 too can be driven without 00ALLINONE.c.
- The band-edge test uses only the gen3 API (26-51) and links `edmd_gen3.c` alone.
- Recorded builds of a71f394 (hspist3/experiments_gen3_d52_261014/bin_a71f394/build_record.txt), each with "warnings: 0":
  - `gen3_band_edge_test: cc -std=c11 -O2 ... edmd_core/tests/gen3_band_edge_test.c edmd_core/edmd_gen3.c -lm` (39-40);
  - `edmd_gen3_alone: cc -std=c11 -O3 -ffp-contract=off -Wall -Wextra -c edmd_core/edmd_gen3.c` (43-44);
  - `gen3_m1: cc -O3 -ffp-contract=off ... gen3_m1_harness.c edmd_core/edmd_gen3.c edmd_core/edmd.c -lm` without `-std` (33-34).
- So gen3 links without edmd.c and without 00ALLINONE.c, per the record. The main session rebuilt the objects and both harnesses and linked them in item G.

**References from the engines to 00ALLINONE.c: none found**
- `#include` lines:
  - edmd.c: `edmd.h` plus libc only (1-6);
  - edmd_gen3.c: `edmd_gen3.h` (which includes `edmd.h`, edmd_gen3.h:28) plus libc (66-71);
  - edmd_accelerated.c: `edmd.h` plus libc (34-39).
- Searched `extern`, `cli_`, `g_edmd\b`, `XW1`, `PIXELS_PER_SIGMA` in all three engine files. The only hit is comment text at edmd.c:1257.
- Every `identifier(` in edmd_gen3.c and edmd.c is file-local, a macro, or libc/libm.
  - gen3 uses `abs`, `acos`, `calloc`, `exit`, `fprintf`, `fread`, `free`, `fwrite`, `isfinite`, `isnan`, `malloc`, `memcmp`, `memcpy`, `memset`, `printf`, `qsort`, `realloc`, `snprintf`, plus math through `R_()`.
  - edmd.c also uses `rand`, `getenv`, `fopen`, `fclose`, `fflush`.
- gen3 calls no `edmd_*` function (grep `\bedmd_[a-z_]*(` in edmd_gen3.c: none). It uses edmd.h only for types and `EDMD_MAX_DIVIDERS`.

**File-scope mutable state a library caller inherits**
- gen3
  - `static EDMD3* coll_hook_S`, `coll_hook_fn`, `coll_hook_ctx` (edmd_gen3.c:704-706; "one state at a time", edmd_gen3.h:186).
  - `exit(3)` on heap OOM (197).
  - Audit `printf` to stdout (893).
  - No other mutable statics: the rest are `const` (178, 293, 750).
- gen2
  - `g_edmd_avalanche_warning_count` (14), `g_edmd_legacy_resched` (35), `g_edmd_resched_audit` (58).
  - `g_edmd_contact_audit`, read once from `getenv("HD_CONTACT_AUDIT")` (45, 1750).
  - The trace ring (109-113) and `g_edmd_evlog` (1264-1265).
  - libc `rand()` in the heat bath (1082-1159). It shares its seed with the driver's `srand` (00ALLINONE.c:17834, 18822, 23501-23503).

### C5

Paths as in B1.

**Asked explicitly: how the cell grid is tied to the box (gen3).**
- **Grid.** `S->gw = (int)R_(ceil)(S->boxW / w); S->gh = (int)R_(ceil)(S->boxH / w);` (`edmd_gen3.c:1136`). The origin is the box corner (0, 0):
  - at load, `int cx = (int)R_(floor)(P[i].x / S->w), cy = (int)R_(floor)(P[i].y / S->w);` (1220), clamped to the grid (1221-1222);
  - the position is stored relative to the cell corner, `A->xi = P[i].x - (real)cx * S->w; A->zeta = P[i].y - (real)cy * S->w;` (1224).
- **Cell width w.** An integer number of px and at least one diameter (1116-1117). Default `EDMD3_DEFAULT_CELL_PX 32.0` (`edmd_gen3.h:34`); the window scene uses max(32, ⌈2R⌉) (`00ALLINONE.c:1545-1548`). The cell capacity is derived from the one diameter (1140).
- **No rebuild.** The grid is never rebuilt. Crossing events move disks between cells (757-795) and never leave the grid (458-461).
- **The origin shift is temporal, not spatial.** Every `EDMD3_ORIGIN_SHIFT 8192.0` internal units (= 341.3 σ-time, `edmd_gen3.h:35`), `origin_shift` synchronises all disks. It then subtracts the shift from `now`, the disk time stamps, the heap times and the body times (`edmd_gen3.c:873-881`). Positions and cell indices are untouched.
- **What a static segment list per cell would have to join:**
  - the three call sites of `wall_candidate`: `predict_all` (498), `exec_cross` (786-787) and load (1248);
  - the validator's wall checks (585, 640);
  - the schedule audit's wall block (1050-1056);
  - the ledger (1533).
  Details are in the numbered list at the end of this item.


| engine | cell list | cell size and count | tie to the box | used for neighbours? | rebuilt when |
|---|---|---|---|---|---|
| gen2 edmd.c | `Cell* grid; int gw, gh;` `double cell_size;` (176-177) | `S->cell_size = (S->prm.cell_size>0.0)? S->prm.cell_size : fmax(2.5*S->prm.radius, 1.0*S->prm.radius);` (1462); `S->gw = (int)fmax(1.0, floor(S->prm.boxW / S->cell_size));` (1463), gh likewise (1464). Driver passes `prm.cell_size = 0.0; // auto` (00ALLINONE.c:21597) | origin at the box corner: `int cx = (int)floor(S->P[i].x / S->cell_size);` then `cx = clampi(cx, 0, S->gw-1);` (431-434) | **No.** The grid is only written (373, 375, 435-437) and freed (366); no other `grid[`/`idx[` use (grep `grid\[|->grid|idx\[|->idx` in edmd.c). `schedule_for` loops `for (int j = 0; j < S->prm.N; ++j)` (1006-1009); full reschedule takes all pairs (1027-1031). The comment at 1000 says "9-cell stencil" | `grid_build` after every executed event (1756, 1771, 1778) and in `reschedule_all_internal` (1023), O(N) each. It also clamps disks into the box and out of slabs (384-404), counted (425-426). The GUI's EDMD loop calls `edmd_reschedule_all` every fixed step (00ALLINONE.c:22267) = `reschedule_all_internal` (edmd.c:1883), all pairs |
| accelerated | same struct (113-114) | same formulas (1035-1037) | `particle_cell_xy` 572-578 (floor and clamp) | yes: `schedule_ab_neighbors_9cell` (624-642), plus cell-crossing events `EV_CC` (580-622) | `grid_build` after every event (1190, 1200, 1209); full reschedule after any divider or piston event (1206-1207). Backend refused on main (00ALLINONE.c:1357-1360) |
| gen3 | `int *cell, *ccount;` fixed slots, "ncell * ccap disk indices; occupancy" (edmd_gen3.c:139); capacity `{ const real q = w / S->d + 1.0; S->ccap = (int)(q * q) + 4; }` (1140); overflow is fatal (600, 778, 1229) | w must be an integer number of px and at least the diameter: `if (!(w >= 2.0 * prm->radius)) FAIL` (1116), `if (w != R_(floor)(w) \|\| w > 1048576.0) FAIL` (1117). Default `EDMD3_DEFAULT_CELL_PX 32.0` (edmd_gen3.h:34); the window scene uses `fmax(EDMD3_DEFAULT_CELL_PX, ceil(2.0 * p->radius))` (00ALLINONE.c:1545-1548). `S->gw = (int)R_(ceil)(S->boxW / w); S->gh = (int)R_(ceil)(S->boxH / w);` (1136) | at load `int cx = (int)R_(floor)(P[i].x / S->w), cy = (int)R_(floor)(P[i].y / S->w);` (1220), clamp (1221-1222), local coordinates (1224); `cell_of` = `cy * S->gw + cx` (249). Outer walls predicted only for cells that can reach them, `wall_candidate` (265-272), e.g. `case 1: return (real)(cx + 1) * S->w >= S->boxW - S->R;` (268) | yes: 9 cells at a velocity change (`predict_all` 499) and the 3 newly adjacent cells at a crossing (780-785) | **never rebuilt as a whole.** Membership changes through crossing events that carry the new cell (`exec_cross` 757-795: `A->xi -= S->w; A->cx++;` 764). Crossings never leave the grid (`A->cx + 1 < S->gw`, 458-461). Safety net: `local_check` re-files a disk found outside its cell (592-603), counted in `cell_repair` (594). Dividers and pistons are not in the cells; they use bands of cell columns (`band_columns` 522-534, `band_compute` 535-543), renewed at band expiry (`exec_band` 857-863; constant v expires after one cell width, `ds = S->w / R_(fabs)(v)`, 518) or at a body velocity change (`obj_changed` 556-561) |
| time-stepped (00ALLINONE.c, `update_particles_with_substepping` 9202) | `cell_head/cell_next/particle_cell`, malloc'd and freed per call (9283-9285, 9852-9854) | `const float cell_size = fmaxf(4.0f * PARTICLE_RADIUS, 1.0f);` (9274); `int grid_cols = (int)ceilf(domain_width  / cell_size);` with `const float domain_width  = (float)(XW2 - XW1);` (9275-9278) | `float x_local = (float)(X[i] - XW1);`, floor and clamp (9477-9483) | yes: 3x3 (9497-9504), or a ring of `cli_pp_toi_neighbor_ring` (9543-9558) | every substep (loop 9311; rebuild 9473-9487) |

**What a polyline boundary would have to live with.** Each point follows from the code cited.

1. **The hull must lie inside [0, boxW]×[0, boxH].** In gen3 the grid covers that rectangle (1136); the four outer walls are always predicted for candidate cells (498, 787, 1248); crossings never leave the grid (458-461); `edmd3_load` refuses disks outside the box (full check 1239-1244). A polyline adds walls; it cannot replace the box.
2. **Each cell needs the list of segments it can reach.** gen3 limits wall work by a closed-form test on the cell index (265-272). A static polyline would need a precomputed per-cell list, consulted where `predict_all` (495-501) and `exec_cross` (786-787) call `wall_candidate`.
   - It needs a new event type. The enum is `T_CROSS = 0 ... T_PAIR = 5` (89); ties are ordered by (t, type, a, b) (186-191); the event hash includes the type (1283). A type appended after `T_PAIR`, with the segment index in `b`, is never pushed in a box-only run, so that run's event order and hash are unchanged (deduced from 186-191, 1283).
   - Any field added to `struct EDMD3` changes the checkpoint format: the header stores `sizeof(EDMD3)` (1414), the whole struct is written (1430), and the reader refuses a size mismatch (1444).
3. **Bands do not carry over to segments.** They work on x-columns only (522-534; `predict_column` walks all rows of one column, 545-554), so they fit vertical bodies only. Static segments need no bands.
4. **A segment predictor must rebuild absolute positions.** Disk positions are cell-local and lazy (`Disk3.xi, zeta, tau` 101-103; `advance` 276-280); the predictor must use x = cx·w + xi as `predict_obj` does (482). Reflection about a general normal changes vx and vy; `exec_wall` flips one component (742-744). The momentum ledger is per axis (`ledJ[2]`, 167); per-wall records are fixed at 4 (`J_wall[4]` edmd_gen3.h:104; `wall_impulse[4]` edmd_gen3.c:150).
5. **Tolerances and the validator know only the box and the bodies.** `tol_face` uses ulp(boxW) (808); `local_check`, `reflect_into_box`, `full_check` (567-671) check the four walls and the bodies, and `reflect_into_box` repairs against the box only (609-619). A polyline needs its own checks and repair rule.
6. **One radius sets the cell width and capacity.** w ≥ 2R (1116); one `d2` for every pair (429, 1134); `ccap` from the one diameter (1140). With two radii the 9-cell stencil is complete only if w is at least the largest diameter, and `ccap` would have to use the smallest diameter.
7. **gen2 has no effective cell list.** A polyline there costs O(segments) per disk per scheduling, like `schedule_walls` (644-665), and the box clamp in `grid_build` (384-387) still runs after every event.
8. **The time-stepped mode has no wall events.** Its grid is sized from the rectangle (9275-9278); boundaries are position tests and clamping (8423-8494).

### C6

**Asked explicitly: where the brute-force schedule audit lives. In the engine.**
- **Functions.** `static void schedule_audit(EDMD3* S){` (`edmd_gen3.c:995`) and `static void schedule_audit_bodies(EDMD3* S){` (937) are file-static functions of `edmd_gen3.c`.
  - Walls and pairs are recomputed "brute force from ABSOLUTE positions (an independent computation)" (comment at 1046). The wall block (1050-1056) repeats the four gap formulas. Pairs use the shared `pair_rule` (comment at 406).
  - The contact audit is `static void contact_audit(EDMD3* S, const Ev3* e){` (675).
- **When they run:** from `edmd3_advance_to` every `audit_every` events (1291); at band events when the body audit is on (1287); after an API body change (1502).
- **Public switches:** `edmd3_set_schedule_audit`, `edmd3_set_schedule_audit_bodies`, `edmd3_schedule_audit_now`, `edmd3_schedule_audit_stats` (`edmd_gen3.c:1362-1365`; `edmd_gen3.h:200-203`) and `edmd3_set_contact_audit` (`edmd_gen3.h:189`).
- **The harnesses only switch them on and print the counts:** `gen3_m1_harness.c:120, 131, 137-138`; `gen3_m2_harness.c:192`; `gen3_band_edge_test.c:29, 36-37`; `gen3_checkpoint_test.c:52, 57`. What the harnesses add is the run pattern A/B/C: audits on, audits off, and again, with equal hash and final state required (`gen3_m1_harness.c:3-6`).
- **Consequence for a segment primitive in a fork (H, option d).** The brute-force template is in the engine itself. The fork adds a segment class beside the wall block and a segment branch in `contact_audit`; the harness pattern is reused unchanged.

**Event hash (gen3)**
- Function: 64-bit FNV-1a.
  - `static inline void fnv(uint64_t* h, const void* p, size_t n)` with prime `1099511628211ULL` (edmd_gen3.c:286-289).
  - Offset basis `S->hash = 1469598103934665603ULL;` (1158).
- Where: in `edmd3_advance_to`, after the stale, past and stagnation checks and before execution: `{ const double th = (double)e.t; fnv(&S->hash, &th, sizeof th); } ... fnv(&S->hash, &e.type, sizeof e.type); fnv(&S->hash, &e.a, sizeof e.a); fnv(&S->hash, &e.b, sizeof e.b);` (1283).
- Covers every executed live event: PAIR, WALL, DIV, PISTON, and also the bookkeeping CROSS and BAND events. For each it takes the double of the origin-relative time `e.t` and the bytes of `int type, a, b`.
- Not covered: stale events (1272), past events (1273), and the event that tripped the stagnation guard (pushed back at 1278, before 1283).
- It hashes no velocities or positions; they enter only through later event times.
- It depends on:
  - the cell width (CROSS times and their a, b);
  - the band logic (BAND events);
  - the origin shift (relative `e.t`, 873-881);
  - `int` size and byte order.
- Read through `uint64_t edmd3_event_hash(const EDMD3* S)` (edmd_gen3.h:159; .c:1317).
- A separate hash in the driver, `g3_init_hash` (00ALLINONE.c:1565-1574; declared 1481): FNV-1a of the gen3 start state.
  - Covered: `boxW`, `boxH`, `radius`, `N`, `divider_count`, the divider x/thickness/mass/vx arrays, `t0`, and every x, y, vx, vy.
  - Not covered: pistons and springs.
  - Used as the restart identity check (17614, 17669).

**Energy and momentum ledger (gen3)**
- Struct `EDMD3_Ledger {P[2], P0[2], J[2], scale_P[2], E, E0, W, scale_E, J_wall[4], J_div[], J_spring[], J_piston[2], J_api}` (edmd_gen3.h:101-109); definition in the comment at 95-100.
- P and E cover all bodies of finite mass: unit-mass disks, free and spring dividers, finite-mass pistons; E includes the spring potential (`mech_now`, edmd_gen3.c:1189-1210).
- J is the impulse from outside:
  - outer walls (742, 744);
  - mass-0 dividers and pistons (833-834);
  - spring anchors (314);
  - API changes (1496).
- W is the work of driven mass-0 bodies (834) and of API changes (1496).
- `void edmd3_ledger(EDMD3* S, EDMD3_Ledger* L)` (1523-1537) assembles these from `mech_now` and the accumulators.
- Rounding scale: `U_ROUND` = 2^-53 (293) times the sum of |terms| (1531-1532). Its terms are added in `exec_pair` (726-727), `exec_wall` (742-745), `exec_obj` (828-836), `obj_advance` (315-316) and `set_motion` (1497).
- The accumulators `ledJ`, `ledJw`, `ledSP`, `ledSE` are read only inside their own accumulation lines and in `edmd3_ledger` (1529-1536). `ledW` is also read by `energy_bound_update` (812).
- That read feeds `E_bound` → `v_ref` → `c_tol`, `tol_face` (`tol_update` 798-809). These tolerances only decide which Health counter an at-once event bumps:
  - `overlap_repair` vs `contact_now` (432-433);
  - `wall_overdue` vs `wall_contact_now` (451);
  - `obj_overlap_repair` vs `obj_contact_now` (488-489).
- Some of those counters enter `clean` (1540-1542). No event time or velocity depends on the tolerances.

**Printed lines** (hspist3/00ALLINONE.c)

| line | where | content |
|---|---|---|
| `[EDMD3-LEDGER]` | `g3_run_extras` 17018-17026, printed just before every run record (17436) | momentum x/y and energy residuals `(P-P0)-J`, `(E-E0)-W`, scale, ratio, E0, E, W |
| `[EDMD3-GAP]` | 17027-17032 | contact-audit maxima when there are audited events; the audit is on iff `HD_CONTACT_AUDIT` is set (`g3_build`, 1580) |
| `[EDMD3-AUDIT]` | 17033-17042 | schedule-audit counts when audits ran (`--gen3-audit-every`, 1581) |
| `[EDMD3-HEALTH] <id>: clean=… hash=%016llx t_end=…` | `g3_run_record` 17433-17452 | every Health counter, `validator_every`, the hash, `t_end`, `engine=gen3`, build, target, timings. `clean = edmd3_health_clean` (edmd_gen3.c:1538-1543): all safety nets, validator and body counts are 0 and the run is not fatal. |
| `[EDMD3-HEALTH]` call sites | 18462 (speed of sound), 19621 (energy transfer), 23441 (replay), 22037 / 23000 / 21417 / 1736 (scene) | – |
| `[GUI-STATS] … events=%ld …` | 22991-22999 | `events` = `g3_gui_events()` = `ev_pair+ev_wall+ev_cross+ev_div+ev_piston+ev_band` of the current gen3 state (11864-11868); 0 without a gen3 state. Printed under gen3 or with `--gui-stats`. |
| `[GUI-HASH]` | 21344-21348 | hash, time and events at a `--gui-keys` "hash" step |
| `[EDMD3-REPLAY]` | 23431-23438 | driver-path hash and event count against the harness: MATCH only if both are equal (23434) |
| `[EDMD-ENERGY]` (gen2 and gen3) | `edmd_energy_audit` 2476-2488 | gas KE and divider KE; only with `HD_KE_TRACE` (2477) |

gen2 has no ledger and no hash in the engine. It has work accumulators (`resolve_wall`, edmd.c:1395-1433) and `heat_bath` (1320).

**How Stage 1 and the gates use them**

Stage 1 worker, per trajectory (hspist3/cluster/gen3_stage1_261012/s1_worker_koa.sh):
- `health_ok` (42-48) requires:
  - exactly 1 `[EDMD3] built #` line;
  - exactly 1 `[EDMD3-HEALTH]` line, with `clean=1` and `engine=gen3`;
  - no line matching `HEALTH='EDMD-HEALTH|forced_advance|clamp_repair|overlap_repair|wall_overdue'` (22) other than `[EDMD3-HEALTH]` lines.
- A failure moves the run to `.failed_run*` (92-93).
- The worker sets `HD_CONTACT_AUDIT=1` (67) and `HD_KE_TRACE=1` (82), so Stage 1 logs also carry `[EDMD3-GAP]` and `[EDMD-ENERGY]` lines.
- The log is appended to the cell's run.log (`append_run`, 56-59; called 100).
- The worker does not compare hashes.

Loader guard (hspist3/validation/edmd3_health_guard.py):
- It refuses gen3 data whose record has `clean` other than 1 (120-121), has more build lines than records (122-123), has no record at all (127), or has no matching per-run record (138).
- Regex `RECORD = re.compile(r"\[EDMD3-HEALTH\] (.*?): clean=(\S+)")` (47); rule text in the docstring at 26-34.
- No call of it found in cluster/gen3_stage1_261012/*.py or cluster/confinement_20261013/reduce_B.py (searched `edmd3_health_guard`, `edmd_acc_guard`).

Gate G6(d) (hspist3/cluster/gen3_s1_gate_261012/s1_gate_testg_compare.py:12-37):
- Per Test G gen3 line, it compares `hash`, `clean` and the counts `EVK = ("ev_pair", "ev_wall", "ev_cross", "ev_div", "ev_piston", "ev_band")` (9) between the gate binary a643bc7 and the campaign binary f72d764 (run at s1_gate_koa.sbatch:77-92).
- The verdict requires 12 of 12 (s1_gate_verdict_261012.py:50-58).

Gate G6(c):
- E0 / rule 4, gen2 byte identity (s1_gate_koa.sbatch:55-63; verdict 29-38)
  - `audit_runs_261007.py` compares `wall_x_positions_*_run0.csv`, `speed_of_sound_psi6.csv`, `ev_9700.csv`, `tr_9700.csv`, `red_9700.csv` with `cmp -s` (hspist3/cluster/resched_gate_261005/audit_runs_261007.py:112-120).
  - Cases ctrl_min, ctrl_leg, afix, afix_leg, against the 7b08827 reference, with and without `--engine=gen2`.
- Harness outputs (s1_gate_koa.sbatch:64-76; verdict 39-44)
  - The SHA-256 of the M1/M2 harness outputs is checked against `koa_harness_ref_a643bc7.sha256` (lines 2-4).
  - The M1 output contains the per-cell hash tables (gen3_m1_harness.c:206-210).
- Cross-node determinism (verdict 45-49)
  - A separate job: s1_xnode.sbatch, then s1_xnode_compare.sbatch:25 runs `cluster/gen3_koa_261009/g3_xnode_compare.py`.
  - That script requires the files byte-identical (`filecmp.cmp`, 50) and the event hashes equal (53-54) across nodes, with each task's record clean=1 (21, 40).
- PASS rule: verdict 55-57.

Other readers:
- Part A: `pa_reduce.py` refuses `clean` other than 1 and records the hash (hspist3/cluster/gen3_partA_261014/pa_reduce.py:115-118, 133).
- `[EDMD3-LEDGER]` is not read by any Stage 1 or gate script: searched `EDMD3-LEDGER` and `LEDGER` in cluster/gen3_stage1_261012, cluster/gen3_s1_gate_261012 and validation (only unrelated `LEDGER` print headers in validation/paper2_level4b_*.py).
- It is parsed in `ledger()`, hspist3/cluster/gen3_koa_261009/koa_k2_m5_rates_261010.py:34-35, which cluster/gen3_m5_solid_261012/m5_solid_verdict_261012.py uses (17, 37). Also in experiments_gen3_m3_261009/stageA_tables.py:95.
- Not read in full: merge_s1_koa.py, make_s1_manifests.py, reduce_s1_koa.py.

**What must keep passing**
- hash and count equality for identical inputs and binary (G6(d), cross-node);
- `clean=1` in every record (worker, guard, Part A);
- gen2 output byte identity (E0);
- the harness-output SHA-256.

What a source change does:
- The registration pins the commit and the binaries' SHA-256: `commit f72d764…` and `sha256 … 00ALLINONE` etc. (hspist3/cluster/gen3_stage1_261012/build_record_f72d764.txt:2-8). That record is itself in the Stage 1 input manifest (s1_inputs.sha256:4), and the gate job refuses other binaries (s1_gate_koa.sbatch:38-41).
- A byte change to `edmd_core/*` or 00ALLINONE.c leaves the registered binaries and their outputs as they are. It produces a build that is not the registered one, and these checks would have to be rerun for it.
- Compiling or copying the files unmodified into a separate build under PHYSICALINFORMATIONDYNAMICS/ touches none of this.

## D. Renderer

### D1

Shorthand: `00ALLINONE.c:n` = `hspist3/00ALLINONE.c:n`.

**Lifecycle**

| Step | path:line | Verbatim / note |
|---|---|---|
| globals | 00ALLINONE.c:763-764 | `SDL_Window* _window = NULL;` `SDL_Renderer* _renderer = NULL;` |
| init function | 00ALLINONE.c:5759-5816 | `void initSDL() {` |
| SDL_Init | 00ALLINONE.c:5760 | `if (SDL_Init(SDL_INIT_VIDEO) < 0) {` |
| SDL_CreateWindow | 00ALLINONE.c:5780-5787 | title `"Hard Spheres Simulation"`, flags `SDL_WINDOW_SHOWN \| SDL_WINDOW_RESIZABLE \| SDL_WINDOW_ALLOW_HIGHDPI`. Size: `--gui-window` (5766) or 98 % x 95 % of the display's usable area (5767-5777; again after creation, 5804-5814). |
| SDL_CreateRenderer | 00ALLINONE.c:5792, 5796 | `SDL_RENDERER_ACCELERATED`, then `SDL_RENDERER_SOFTWARE` with the message `"[GUI] no accelerated renderer (%s); using the software renderer\n"` (5795). |
| GL | not found | grep for `SDL_GL_`, `SDL_WINDOW_OPENGL`, `glBegin`, `glDraw`, `GL/gl`, `OpenGL`, `GLEW` in 00ALLINONE.c: no match. The link line still links GL: `LIBS_BASE       := -lGLEW $(GL_LIBS) $(SYS_LIBS)` (hspist3/Makefile:57; the pinned build's line in hspist3/experiments_gen3_d52_261014/bin_a71f394/build_record.txt:30 ends `-lGLEW -framework OpenGL -framework CoreGraphics`). |
| call site | 00ALLINONE.c:23514-23531 | `if (!cli_headless) {` `initSDL();` `TTF_Init();`, then `font = TTF_OpenFont("Roboto-Regular.ttf", 18);` (23520); `return 1;` if it fails (23521-23524), i.e. before the run dispatch at 23670-23682. The font path is relative to the cwd, also at 23527 and 14366-14367. |
| time scale | 00ALLINONE.c:21996 | `if (cli_gui_time_scale > 0.0) time_scale_runtime = (float)cli_gui_time_scale;` |
| event loop | 00ALLINONE.c:21999, 22010 | `while (running) {` ... `while (SDL_PollEvent(&e)) {`. Per event: SDL_QUIT (22011); resize (22012-22019); `if (live_controls_handle_event(&e)) {` `continue;` (22020-22022), so F2/F3 and panel mouse events reach no key handler; then `keysSimulation(e);` `keysPiston(e);` (22023-22024): keysPiston receives every event that reached keysSimulation, also when keysSimulation returned early. |
| physics sub-steps | 00ALLINONE.c:22090-22899 | `while (accumulator >= fixed_dt_runtime) {`. Window: `accumulator += frame_time * time_scale_runtime;` with frame_time capped at 0.1 s (22002-22006). Headless: `accumulator += fixed_dt_runtime;` (22028), i.e. about one step per pass. |
| render phase | 00ALLINONE.c:22901-22986 | `if (!cli_headless) {` (22903) |
| present | 00ALLINONE.c:22971 | `SDL_RenderPresent(_renderer);`; `gui_shot_maybe(gui_steps);` just before (22970); `SDL_Delay(5);` after (22985). |
| second loop | 00ALLINONE.c:20481-20535, 20553-20649 | simple-box preview: own `SDL_PollEvent` (20483) and `SDL_RenderPresent` (20641); it calls the same scene functions (20561-20567). grep finds no other `SDL_PollEvent` or `SDL_RenderPresent` site. |
| teardown | 00ALLINONE.c:23693-23702 | inside `if (!cli_headless) {` |

**What is drawn per frame, in order (00ALLINONE.c:22906-22971)**

0. Before drawing: `rm_frame();` (22908) attaches the return-map collision hook to the current gen3 state while that panel is open (12000-12013).
1. `--render=paper` (parsed at 4725) calls `render_paper_scene();` (22913; defined 2851) instead of items 2-11 below (the `else` branch runs to 22948). It colours disks by compartment and draws them with `SDL_RenderFillCircle` (2866-2874), without the gs_* transform.
2. `draw_clear_screen();` (22915; 14288): black.
3. `gui_begin_scene_view();` (22916; 1220): `gui_scene_fit()` (1191-1211) and `gs_set(gui_scene_s, gui_scene_bx, gui_scene_s, gui_scene_by, gui_ui_scale);` (1224).
4. `draw_coordinate_system(_renderer);` (22917; 14305): red axes from world point (20,20) (14307-14324).
5. `draw_simulation_boundary();` (22918; 14293): box outline `{XW1, YW1, SIM_WIDTH, SIM_HEIGHT}` (14300), red when `heatbath_enabled`, else yellow (14295-14299).
6. `render_particles();` (22919; 7409): every disk via `gs_disk((float)X[i], (float)Y[i], rad_w);` with `const float rad_w = (float)PARTICLE_RADIUS;` (7446) for all disks. Colour, in code order: Szilard preset, species/memory/overlay (7452-7473, ring via `rad_w + 2.0f` at 7466); MB speed palette when `gui_mb_speed_view && !pl_enabled` (7474-7479; `gui_mb_speed_view` default 1 at 881, toggled in the F2 panel at 15568); particle-life type palette (7482-7491); `colorcode_mode` none/temperature/velocity (7492-7507; `--colourcode`, 5293-5298). With the default MB view on, `--colourcode` has no visible effect, because the MB branch `continue`s first (7474-7479).
7. `render_pistons();` (22920; 6947): two green 5-px-wide rects `{ (int)piston_left_x , 0, 5, YW2 }` / right (6948-6952).
8. `if (wall_enabled) draw_wall();` (22921; 8664): the primary divider as a rect of height SIM_HEIGHT (8675), white (8689), red when `wall_x` is outside [XW1, XW2] (8670, 8691); Szilard colours (8677-8687, 8695-8718); extra dividers as rects (8720-8725) in whatever colour was set last (no colour call in that loop).
9. `render_left_spring();` (22923; 6956): zig-zag, only `if (!energy_measurement.spring_enabled \|\| num_internal_walls <= 0) return;` fails (6957).
10. Panels below the box between `gui_begin_panel_view();` (22926; 1234) and `gui_end_panel_view();` (22934), chosen by `gui_dist_mode` (default 3 at 856; key V):
    - 1: `render_velocity_histograms();` (12763): global vx, vy, speed bars plus MB curves (12841-12881).
    - 2: `render_velocity_histograms_by_segment();` (13413), per-compartment vx/vy.
    - 3: `render_segment_distributions_stacked();` (14210): global x-density (13963), per-compartment x-density only with key X (13421), per-compartment speed (13656), per-compartment vx/vy (13036); layout by panel height (14219-14276); below 70 px it falls back to mode 2 (14211-14214).
11. `gui_end_scene_view();` (22935), then UI pass `gui_begin_ui_view();` (22936; 793): `gui_draw_compartment_labels` (21425), `render_energy_hud` (16100), `render_live_controls` (15757, F2/F3 panels), `gui_div_limit_banner` (21518), `rm_draw` (12051, return map when open).
12. `--demo` (and not paper): HUD bar (22952-22969).

Defined but never called (grep finds only the definition line): `draw_maxwell_boltzmann_curve` (12735), `render_fft_histogram_bottom` (16490), `render_density_histograms_by_segment` (13649).

**Keys** (window only; `--gui-keys` requires the window, 21317). Keys 1-5 and I are the gen3 return-map panel (code 12020-12042; description per comment at 11937-11952), not the histograms; 1-5 act only while it is open (`if (!rm.open) return 0;`, 12033). keysSimulation checks the return-map keys first (12125), then, under gen3, `g3_gui_key` (12126).

| Key | Effect | path:line |
|---|---|---|
| window close | `running = 0` | 00ALLINONE.c:22011 |
| F2 | toggle the slider panel; its rows are mouse-driven (incl. the MB-view row, 15565-15575) | 00ALLINONE.c:15622-15628, 15526-15617 |
| F3 | toggle the scene-reset panel; mouse handling. APPLY rebuilds the scene (per comments at 11829, 12129; flag set via `live_scene_prepare_rebuild`, used at 21354) | 00ALLINONE.c:15631-15640, 15647-15744 |
| I | open/close the return map; under gen2 only the note "I (return map) needs gen3" | 00ALLINONE.c:12021-12032 |
| 1 | angle: line of centres / velocity after; clears the map | 00ALLINONE.c:12035 |
| 2 | wall/divider/piston points cyan or hidden | 00ALLINONE.c:12036 |
| 3 | side: both, left, right | 00ALLINONE.c:12037 |
| 4 / 5 | K halved (min 625) / doubled (max 1,280,000) | 00ALLINONE.c:12038-12039, 11954-11955 |
| gen3: W, `[` `]` `\`, B, K, `,` `.` | refused with an on-screen note | 00ALLINONE.c:11891-11895 |
| gen3: ESC / Q | `live_scene_rebuild_requested = 1;` / `g3_gui_quit = 1;` (loop ends at 22030) | 00ALLINONE.c:11896-11897 |
| `--demo`: SPACE, P, S, + / - | pause toggle; piston step (after release); BMP (`demo_screenshot.bmp` in the cwd unless `--demo-shot`); pace | 00ALLINONE.c:12135-12168 |
| SPACE | `simulation_started = 1;` (starts; initial 0 at 921; auto-started headless or with `--auto-release`, 21977-21980) | 00ALLINONE.c:12172-12174 |
| Q (gen2) | `exit(0);` | 00ALLINONE.c:12176-12178 |
| W | `wall_enabled` toggle; `edmd_config_dividers` + `edmd_reschedule_all` | 00ALLINONE.c:12180-12207 |
| R | release held dividers | 00ALLINONE.c:12209-12224 |
| P | pause toggle (P/SPACE ignored while the divider-limit banner is up, 12127-12132) | 00ALLINONE.c:12226-12229 |
| ESC (gen2) | `reset_simulation_state();` | 00ALLINONE.c:12231-12234 |
| `[` `]` `\` | rescale all velocities: T x0.95 (min 0.01) / x1.05 / T = 1.0 | 00ALLINONE.c:12237-12256 |
| G H J O N | Szilard actions (Szilard preset only) | 00ALLINONE.c:12257-12288 |
| B | heat bath toggle; target = current measured T; written into the EDMD params | 00ALLINONE.c:12291-12309 |
| K | Andersen thermostat toggle | 00ALLINONE.c:12312-12319 |
| `,` `.` | Andersen frequency -/+0.01 (0.001..2.0), or Szilard duration x1.5 / /1.5 | 00ALLINONE.c:12324-12362 |
| - = + 0 | time scale /2 (min 0.1) / x2 (max 10000) / 1 | 00ALLINONE.c:12365-12383 |
| V | `gui_dist_mode = (gui_dist_mode + 1) % 4;` | 00ALLINONE.c:12385-12397 |
| Z L F C U | exact bins / marker lines / freeze vmax / fit overlay / uniform baseline | 00ALLINONE.c:12399-12445 |
| Y | Szilard view cycle, otherwise global x-density curve overlay | 00ALLINONE.c:12447-12466 |
| X | per-compartment x-density panel | 00ALLINONE.c:12468-12473 |
| `;` `'` | Szilard duration | 00ALLINONE.c:12475-12499 |
| PgUp / PgDn | window height | 00ALLINONE.c:12502-12507 |
| A / D, S | left piston velocity -/+25 per press (`const float acceleration = 25;` 12514); S: `vx_piston_left = 0;` | 00ALLINONE.c:12627-12643 |
| Left / Right, Up | right piston velocity -/+25; Up stops; all three clear `piston_step_active` | 00ALLINONE.c:12646-12665 |
| T | piston step protocol (Szilard: separation) | 00ALLINONE.c:12666-12672 |
| E / M | print the energy measurement (Szilard: erase memory) / Szilard measure or particle-life tensor mode | 00ALLINONE.c:12674-12693 |
| key-up A D Left Right | velocity moved 25 towards 0 | 00ALLINONE.c:12696-12711 |
| `--gui-keys=STEP:KEY,...` | pushes key-down/up events at scene steps; also `apply-gen3`, `apply-gen2`, `hash` | 00ALLINONE.c:21297-21368 (parsed 5226) |

Simple-box preview keys: ESC/Q quit, SPACE/P pause, T step, A autorun, G +20 steps, H help, - / = + / 0 time scale (max 256) (00ALLINONE.c:20486-20531).

**Is drawing a function of a state struct? No. Is it interleaved with the physics? The draw calls are not; data collection and input are, in places.**

- No draw function takes a state struct. Signatures (verbatim, line):

| Signature | Line |
|---|---|
| `void render_particles() {` | 7409 |
| `void render_pistons() {` | 6947 |
| `static void render_left_spring() {` | 6956 |
| `void draw_wall() {` | 8664 |
| `void draw_clear_screen() {` | 14288 |
| `void draw_simulation_boundary() {` | 14293 |
| `void draw_coordinate_system(SDL_Renderer *renderer) {` | 14305 |
| `void render_velocity_histograms() {` | 12763 |
| `static void render_velocity_histograms_by_segment(void) {` | 13413 |
| `static void render_segment_distributions_stacked(void) {` | 14210 |
| `static void render_velocity_histograms_by_segment_panel(int y0, int h) {` | 13036 |
| `static void render_density_histograms_by_segment_panel(int y0, int h) {` | 13421 |
| `static void render_speed_histograms_by_segment_panel(int y0, int h) {` | 13656 |
| `static void render_density_histogram_global_panel(int y0, int h) {` | 13963 |
| `static inline void render_energy_hud(SDL_Renderer *r, TTF_Font *f) {` | 16100 |
| `static void render_live_controls(SDL_Renderer *renderer, TTF_Font *font) {` | 15757 |
| `static void rm_draw(SDL_Renderer* r, TTF_Font* f){` | 12051 |
| `static void gui_draw_compartment_labels(int left_particles, int right_particles) {` | 21425 |
| `static void render_paper_scene(void) {` | 2851 |

- They read file-scope globals: particle arrays `double X[NUM_PARTICLES];` ... `double Radius[NUM_PARTICLES];` (748-755), `particles_active` (156), `PARTICLE_RADIUS` (macro, 239, depends on `particle_scale_runtime`, 237); bodies `wall_x` (694), `extra_wall_positions`, `piston_left_x` / `piston_right_x` (e.g. 6949-6951, 8720-8721); box `XW1`, `YW1`, `SIM_WIDTH`, `SIM_HEIGHT` (1198, 1205-1206, 14300); display flags `gui_*`, `segment_*`, `_renderer`.
- Renderer argument: in three regions (6939-6986, 8656-8733, 12727-14337) `SDL_RenderDrawLine`, `SDL_RenderFillRect`, `SDL_RenderDrawRect`, `SDL_RenderDrawPoint`, `draw_text`, `draw_text_clipped`, `TTF_SizeUTF8` are macros onto gs_* that drop the renderer, e.g. `#define SDL_RenderDrawLine(r_, ...) gs_line(__VA_ARGS__)` (12728); the gs_* primitives draw on `_renderer` (gs_fill_dev 1076-1083). Inside those regions `SDL_SetRenderDrawColor` still uses the argument (`SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255);` 14307), and `draw_thick_line_vertical` has `(void)renderer;` (12753). The UI-pass functions outside the regions use their argument (e.g. `hud_text_wrapped(r, f, ...` 16333; `SDL_SetRenderDrawColor(r, 10, 12, 22, 238);` 12059).
- No draw call sits in the step loop: awk over 22090-22899 for `render_`, `draw_`, `SDL_Render`, `SDL_SetRender`, `gs_*(` finds none. After each EDMD step the engine state is copied into the globals for drawing: 22284-22291 (EDMD, comment "copy back for rendering/diagnostics"), 22338-22345 (hybrid). Under gen3 `edmd_particles` is `edmd_backend_particles` (2041), which returns `g3_particles()` (1671; 1599-1607); `g3_check_untouched()` (1595-1598, called at 1729) stops the run if a reader wrote into that copy.
- Exceptions to the separation: (a) the return-map panel's data is gathered inside gen3's advance by the collision hook `rm_hook` (11989-11999), attached by `edmd3_set_coll_hook(g_edmd3, rm_hook, NULL);` (12011) while the panel is open; per comment at 11950-11952 this leaves event hashes unchanged (not re-tested here). (b) Under gen2, key handlers change physics state directly: velocity rescale (12237-12256), divider reconfiguration (12195-12201), heat-bath params (12298-12299).
- The draw functions are not pure readers: the panel functions call `sync_all_wall_positions();` (13042, 13427, 13662, 13968), which writes `all_wall_positions[]` (3724-3736); they keep static EMA arrays across frames (13095-13108, 13479-13492, 14014-14025); `gui_begin_panel_view` overwrites the globals `HIST_HEIGHT`/`HIST_WIDTH` (1236-1238, restored 1244); the UI pass overwrites `MAX_X`/`MAX_Y` (`gui_ui_coords_begin`, 781-787, restored 788-791).
- Headless runs the same `simulation_loop` step code with the render block skipped (22903); differences in the loop: the accumulator rule (22001-22029) and `g3_soft_window = !cli_headless;` (21997: a refused gen3 build pauses the window instead of stopping).

Runtime: not re-run by this verifier. The drafting agent's window run (E1, run 4) reported the software-renderer message (format at 5795) and `Failed to load font: ...` (format at 23522), rc 1; its 00_COMMAND.md (389 B) is in PHYSICALINFORMATIONDYNAMICS/_m0_scratch/view_run and shows that command.

### D2

**Drawing primitives in hspist3/00ALLINONE.c**

| Primitive | path:line | Radius / colour / notes |
|---|---|---|
| `static void gs_disk(float cx, float cy, float r_world) {` | 00ALLINONE.c:1109 | Filled disk, one `SDL_RenderFillRect` per scanline (1119-1127); a square below 1.5 device px (1114-1117); minimum radius `gui_particle_draw_min_px` (F2 slider row at 1041; applied at 1113). Radius is a float argument; colour is the current draw colour set by the caller (e.g. 7459-7461). |
| `void SDL_RenderFillCircle(SDL_Renderer* renderer, int centerX, int centerY, int radius_pixels) {` | 00ALLINONE.c:7351 | point by point (7352-7359); called only from `render_paper_scene` (2874); outside the gs transform. |
| `static int gs_line(int x1, int y1, int x2, int y2) {` | 00ALLINONE.c:1090 | one segment, integer world-px endpoints; axis-aligned segments become filled rects (1093-1094), others `SDL_RenderDrawLine` repeated t times (1095-1097). |
| `gs_fill_rect`, `gs_draw_rect`, `gs_point`, `gs_text` | 00ALLINONE.c:1084, 1100, 1108, 1135 | rect, rect outline (4 x gs_line), point, text (`TTF_RenderUTF8_Blended`, 1139). gs_* block: 1062-1172. |
| `void draw_text(SDL_Renderer *renderer, TTF_Font *font, const char *text, int x, int y, SDL_Color color) {` | 00ALLINONE.c:14343 | text via `SDL_CreateTextureFromSurface` + `SDL_RenderCopy` (14353-14356) |
| streaming texture | 00ALLINONE.c:12091, 12097 | `SDL_CreateTexture(r, SDL_PIXELFORMAT_ARGB8888, SDL_TEXTUREACCESS_STREAMING, RM_BINS, RM_BINS)`, `SDL_UpdateTexture` + `SDL_RenderCopy`, nearest scaling (12093). The only `SDL_CreateTexture` call. |
| alpha blending | e.g. 00ALLINONE.c:12764, 13041, 22958 | `SDL_SetRenderDrawBlendMode(..., SDL_BLENDMODE_BLEND)` (10 sites) |
| colour maps | 00ALLINONE.c:7374, 7364, 12043 | `static inline void gui_speed_heat_rgb(float a01, Uint8 *out_r, Uint8 *out_g, Uint8 *out_b) {`, `static inline void set_draw_color_blue_purple_red(float a01) {`, `static void rm_colour(double t, Uint8* r, Uint8* g, Uint8* b){` |

- **Not used:** `SDL_RenderGeometry`, `SDL_RenderDrawLines`, `SDL_RenderDrawPoints`, `SDL_RenderFillRects`, the float `*F` variants, GL. Method: `grep -o "SDL_Render[A-Za-z]*"` over 00ALLINONE.c gives only RenderClear 2, RenderCopy 3, RenderDrawLine 55, RenderDrawPoint 8, RenderDrawRect 25, RenderFillCircle 3 (the program's own), RenderFillRect 60, RenderPresent 2, RenderReadPixels 2, RenderSetScale 3. The installed SDL2 headers (/opt/homebrew/include/SDL2, version 2.30.12) declare `SDL_RenderGeometry` (SDL_render.h:1644), so it is available but unused.
- **No arc or point-list routine:** `grep -n -i` for `\barc\b`, `polyline`, `_arc(`, `bezier`, `circle(` matches only the comment `// draw as polyline zig-zag` (6965) and the FillCircle lines. Connected segments exist only as inline loops: the spring zig-zag (6966-6975), `paper_draw_coil` (2838-2849) and computed curve overlays inside panel functions (e.g. the density curve 14142-14158). None takes a point list.
- **Colours are literals per draw function** (6948, 8689, 14296-14298, 14307); there is no colour or population parameter.
- **Radius:** `const float rad_w = (float)PARTICLE_RADIUS;` (7446) for every disk. `double Radius[NUM_PARTICLES];` (755) is not read by `render_particles` and is set to `PARTICLE_RADIUS` at every seeding site (e.g. 6201, 6564). Particle-life has `static float pl_particle_radii[MAX_PARTICLE_TYPES];` (1260), set to 1.0 at 8052; grep finds no reader. The engine has one radius: `double radius;         /* disk radius (all equal) */` (hspist3/edmd_core/edmd.h:17); `EDMD_Particle` holds x, y, vx, vy, coll_count only (edmd.h:96-100).
- **Walls drawn:** axis-aligned rects only: dividers 8675, 8723; box 14300; pistons 6949-6951.
- **Outside the build:** hspist3/glgraph.[ch] (glgraph.c 1731 lines) is an older SDL2 primitive library with argument-driven `circ`, `circf`, `ellipse`, `rect`, `bresh` (glgraph.h:94-109). It is not on the 00ALLINONE build line (hspist3/Makefile:119-120) and not included by 00ALLINONE.c (includes 1-22). Not audited further.

**Linking constraint, common to (a)-(c).**

- Every primitive and draw function is in hspist3/00ALLINONE.c, the translation unit with `int main(int argc, char* argv[]) {` (23489). gs_* and most panel functions are `static` (1090, 1109, 13036, 13421, ...). The scene functions (`render_particles` 7409, `draw_wall` 8664, `render_pistons` 6947) have external linkage but no header declares them (00ALLINONE.c:1-22 includes only system, SDL, kissfft, edmd_core and experiment_validation headers). The build compiles everything on one cc line (hspist3/Makefile:119-120).
- A separate PID translation unit therefore cannot link this renderer without also linking `main` and the simulation. The only no-edit route is a unity build (`#include` the 23,715-line file with `-Dmain=...`), which brings in every global and `simulation_loop`. Unsure whether it compiles; not tried.
- The frame is a fixed call sequence with no hook for extra drawing (22915-22942).

| Item | Grade | Evidence |
|---|---|---|
| (a) second disk population, own radius and colour | **only by edit** | `render_particles()` loops `for (int i = 0; i < n; i++) {` (7443) over `int n = particles_active > 0 ? particles_active : 0;` (7410) at one radius (7446); colour by global mode (7417-7507). Nearest case: Szilard two colours plus a ring at `rad_w + 2.0f` (7466), same radius. `gs_disk` itself takes radius and draws in the current colour, but it is static (1109) and nothing in the frame calls it for other data. An in-TU (unity-build) adapter could call `render_particles()` twice, swapping `X`/`Y`/`particles_active`, `particle_scale_runtime` (237, which also sets the physics macro PARTICLE_RADIUS, 239) and the particle-life palette (`pl_enabled`, `pl_particle_types`, `pl_type_colors` 1256-1266); it would write physics globals between steps; untested. |
| (b) polyline / arc set | **only by edit** | `gs_line` (1090) is the only line primitive and is static; no routine takes a point list (see above); arcs would have to be split into segments by new code. Endpoints are `int` world px, so geometry is rounded to 1/PIXELS_PER_SIGMA sigma before the transform (`#define PIXELS_PER_SIGMA 24.0f` at 202, active branch for MOLECULAR_MODE 1 / EXPERIMENT_MODE 0, 140-141). |
| (c) scalar field as heat overlay | **only by edit** | No field overlay on the scene exists. The technique exists in `rm_draw`: pixel buffer (12064-12089), streaming texture (12090-12097), colour map (12043), alpha blending available (12058). `rm_draw` reads the static `rm` struct (11956-11964) and draws into a fixed UI rect (12053-12057, 12096). A box-aligned overlay would also need the scene transform `gui_scene_s`, `gui_scene_bx`, `gui_scene_by` (1065, 1204-1206). |

### D3

All code in hspist3/00ALLINONE.c.

| Function | path:line | Bins | Geometry | Inputs | Notes |
|---|---|---|---|---|---|
| `void compute_velocity_histogram() {` | 00ALLINONE.c:10063 | speed \|v\| | 1D, `NUM_BINS 400` (736) over [0, `MAX_VELOCITY 10000.0`] (737): bin width 25 | globals Vx, Vy, particles_active; writes `int velocity_histogram[NUM_BINS]` (757) | feeds `float S_v = entropy();` in `log_energy` (10206-10207). All 4589 data rows of the scratch energy_log.csv (PHYSICALINFORMATIONDYNAMICS/_m0_scratch/view_run) have `S_v: 0.00000000e+00`, i.e. one occupied bin. |
| `void compute_position_histogram() {` | 00ALLINONE.c:10134 | x, y | 2D, 400 x 400, `(int)(X[i] / MAX_X * NUM_BINS)` / `(int)(Y[i] / MAX_Y * NUM_BINS)` (10141-10142): over the window extents MAX_X x MAX_Y (set at 356-357 or 5766-5777), not the box | globals; writes `int position_histogram[NUM_BINS * NUM_BINS]` (758) | Reset loop clears only the first 400 of 160,000 entries: `for (int i = 0; i < NUM_BINS; i++) position_histogram[i] = 0;` (10136), so counts accumulate across calls. Its entropy S_p (10210) is not written (10215-10217). |
| `void render_velocity_histograms() {` | 00ALLINONE.c:12763 | vx, vy, speed | 1D, NUM_BINS: vx, vy over [-MAX_VELOCITY, MAX_VELOCITY] (12782-12783), speed over [0, MAX_VELOCITY] (12785) | globals; local arrays (12773-12775) | panel mode 1 |
| `static void render_velocity_histograms_by_segment_panel(int y0, int h) {` | 00ALLINONE.c:13036 | vx, vy per compartment | 1D; bins `gui_bins_for_velocity` (13023): `floorf(6.0f * sqrtf((float)n))` clamped to 24..256, capped at the compartment's device-px width (13028-13031); range +-vmax_disp from a per-compartment EMA (13143-13165); binned with edge clamping (13180-13186) | globals; static EMA (13095-13108) | |
| `static void render_density_histograms_by_segment_panel(int y0, int h) {` | 00ALLINONE.c:13421 | x per compartment | 1D in x: `float u = (X[i] - xl) / (seg_w > 0.0f ? seg_w : 1.0f);` (13520); bins `gui_bins_for_xdensity` (13008): width / (2 PARTICLE_RADIUS), clamped 16..256 (13009-13018) | globals; static EMA (13479-13492) | shown only with key X |
| `static void render_speed_histograms_by_segment_panel(int y0, int h) {` | 00ALLINONE.c:13656 | \|v\| per compartment | 1D over [0, vmax_disp] (13787-13791) | globals | |
| `static void render_density_histogram_global_panel(int y0, int h) {` | 00ALLINONE.c:13963 | x | 1D between max(XW1, left piston + 5) and min(XW2, right piston) (13973-13977); up to 512 bins (13970, 13979); binning 14041-14053 | globals; static EMA (14014-14025) | |
| `static void log_distributions(double t_abs, double t_rel, int released) {` | 00ALLINONE.c:11762 | x, y, speed | three separate 1D histograms (not 2D) over [0, W] and [0, H] sigma (box: 21913-21914) and [0, vmax] (11770-11796), values clamped into range (11777-11780, 11790-11792); bins `--dist-bins` (default 50, max 2000; 21904-21905) | globals `dist_hist_*`, `dist_box_*`, `dist_vmax_runtime` | the only one that writes a file (11798-11813; E1). Active only with `--dist-log` and the energy-transfer preset in the scene (21903). |
| `static void recompute_segment_stats_counts_and_temperature(void) {` | 00ALLINONE.c:2341 | count, KE, T, vx/vy/speed mean and variance per compartment | compartments = x-intervals between sorted divider positions (`segment_index_for_position`, 3760-3781) | globals | called in the step loop (22187, 22322, 22359); read by the temperature colour mode (7436, 7503), the HUD (16373), compartment labels (21451) and logs (18408, 22465). The panels do their own binning. |
| `static int psi6_compartments(const double* x, const double* y, int N, double sigma, const double* div_x, int ndiv,` / `double cutoff_sigma, Psi6Comp* out, int max, double* global_abs){` | 00ALLINONE.c:16637-16638 | psi6 per compartment | not a binning; all pairs O(N^2) (16645-16653); compartment = number of `div_x` left of the disk (16655) | **arrays and bounds as arguments** | the only analysis function of this kind that takes its inputs as arguments; limited by `EDMD3_MAX_COMPARTMENTS` (16640) |
| return map: `rm.hist` (11961), `rm_point` (11976) | 00ALLINONE.c:11953-11999 | (theta_n, theta_n+1) per disk | 2D, `RM_BINS 90` x 90 (11953), split by side and class | fed by gen3's collision hook (11989-11999) | display only |

**Reusable for a coarse-graining overlay? Not unmodified.**

- Every binning function except `psi6_compartments` reads the globals X, Y, Vx, Vy, particles_active and box/wall/piston globals; none takes arrays or bounds as arguments.
- The bin-count helpers read globals too: `PARTICLE_RADIUS` (13009) and `gs_ax` (13015, 13029).
- The display versions keep static EMA state across frames and call `sync_all_wall_positions()` (13042, 13427, 13662, 13968).
- No function bins in 2D over the box. The only 2D spatial histogram (10134-10147) covers the window extents and does not reset (10136).
- The fixed velocity range MAX_VELOCITY = 10000.0 (737) with 400 bins gives bin width 25 for speed and 50 for vx/vy in `compute_velocity_histogram` / `render_velocity_histograms`.
- All of them are in the translation unit that holds `main` (23489); see D2.
- Size of the binning cores: 11774-11796 (23 lines) and 14041-14053 (13 lines).

## E. Outputs

### E1

All lines in hspist3/00ALLINONE.c unless another path is given. No global output-directory flag: grep for `"--out`, `"--output-dir`, `"--outdir`, `chdir(` finds only `--output-dt` (4509).

**What writes what**

| Trigger | Files | Location | Format |
|---|---|---|---|
| Every run past option parsing. Exceptions: `--version` (4094-4097) and `--help` (4098-4100) exit in the parser; `--engine-replay` (23495) and `--gen3-psi6-snapshots` (23496) return before 23500. | `00_COMMAND.md` | **cwd, unconditional**: `write_command_md_to_dir(".", argc, argv);` (23500), overwritten (`FILE *f = fopen(out_path, "w");` 3192). Further copies: SoS mode0 and mode1 folders (17738-17739), energy-transfer cmd dir (18762), simple-box dir and run dir (20653, 20705), Szilard dirs (3365, 21140). | `# Command`, `- Timestamp: ...`, `- CWD: \`...\``, then argv shell-escaped in a sh block (3212-3220). argv only: no environment variables, no build hash. |
| `run_params.json` (writer 2275-2319) | `run_params.json` | **cwd**: scene (`write_run_params_json();` 21944, plus `wall_position_run_params.json` 21947), `--single-test` (23112), gen3 lattice seeding (17245, 17393), and `initialize_simulation` (6925) only on its single-wall grid path: it returns earlier for honeycomb (6666), random (6669), multi-wall grid (6671) and the hex fallback (6836). Callers of `initialize_simulation`: 3008, 17864, 18824, 20675, 21177, 21575, 23021. Szilard scene and simple box write into their run dir (21939-21942, 20715). | JSON keys mode, particles, L0, height, num_walls, wall_mass_factor, wall_thickness, wall_thickness_vis, temperature, optional temperature_segments, kbt1, kB_effective, kBT, dt, output_dt, time_span, units, radius_sigma, pixels_per_sigma (2282-2317). **No seed, no build hash, no engine.** |
| Interactive scene `void simulation_loop() {` (21567): `--show-simulation` (4105-4108), or the default when no experiment is selected (23680-23681); window or `--headless` (4109) | `energy_log.csv`, `wall_position.csv` | **cwd** (`char energy_log_path[1024] = "energy_log.csv";` 21707, 21708, opened 21713-21714); Szilard preset: `g_szilard_run_dir` (21709-21711) | see "Main formats" |
| scene + `--experiment=energy_transfer` | `energy_transfer_wall_positions.csv` | cwd (21892) | header 21894-21899 |
| scene + energy_transfer + `--dist-log` (4179) | `energy_transfer_distributions.csv` (2224) or `--dist-log-path` (4192), which gets `_seed<seed>` inserted unless it contains `{seed}` or `{ts}` (2225-2233) | cwd by default | header `TimeAbs,TimeAfterRelease,Released,Bins,XMin,XMax,YMin,YMax,VMax,HistX,HistY,HistV` (21927); each histogram a quoted `;`-list (11798-11813); cadence `--dist-dt` (4181) |
| scene + Szilard | `szilard_info.csv` or `--szilard-log` (4258) | Szilard run dir when relative (21722-21727) | header 21731-21754 |
| `--single-test` (4372; dispatch 23670-23671) | `energy_log.csv`, `wall_position.csv`, `run_params.json` | cwd (23100-23101, 23112) | headers 23108-23109 |
| Speed of sound: `--experiments` (4101), `--experiment=speed_of_sound` (4417-4419) and the other presets without their own runner (e.g. 4396-4416), dispatch 23678-23679 | per trajectory `wall_x_positions_L0_<(int)(L0*10)>_wallmassfactor_<M>_run<r>.csv`, written as `.partial_...` and renamed (18508), or `invalid_...` on failure (18493-18494); `speed_of_sound_psi6.csv` (16716), `speed_of_sound_failures.csv` (16821), `speed_of_sound_batch_status.json` (16939); `00_COMMAND.md` and `01_PLOT_COMMANDS.md` in mode0 and mode1 (17738-17741) | `mode0` = `experiments_speed_of_sound/<EDMD\|TIME>/mode0_real_units/<timestamp>` relative to the cwd, **always created** (17724-17725, 17734). `mode1` = `--speed-sound-run-dir` (4116; relative paths resolved against the launch cwd, 3424-3433) or `experiments_speed_of_sound/<E>/mode1_normalized_units/<timestamp>` (17729-17730). Trajectory files go to mode1, because `#define MOLECULAR_MODE 1` (140) selects `mode_folder = mode1_folder` (17743-17747). | `Time,Wall_X,Displacement(σ),Left_Count,Right_Count,L0,eta,Center_X(σ),Seed,Target_Oscillations,Predicted_Frequency,Planned_Steps,Planned_Duration`, then `,KE_L,KE_R` if `HD_KE_COLUMNS` is set (17907; 2515-2519), then `,Box_Width_sigma` (17904-17909) |
| same + gen3 | `psi6_t_L0_...csv` (18147; header begins `t_sigma,phase,psi6_global,psi6_local,mean_neighbors` 18151, compartment columns with `--gen3-psi6-compartments`, 18152). `--gen3-virial-blocks=T` (5180): `virial_blocks_...csv` (17485; header 17488-17490). `--gen3-snapshots=T` (5181): `snapshots_...bin` (17494). | the mode folder | G3SNAP1: `"G3SNAP1\n"`, int32 N, int32 dividers, 4 doubles (box W, H, radius in sigma, interval) (17497-17501); per snapshot: double t, int32 phase, int32 nd, nd divider centres, N x (double x, double y) in sigma (17533-17539). **Positions only.** Reader: `--gen3-psi6-snapshots=PATH` prints psi6 CSV to stdout (23452-23484). |
| Energy transfer (`--experiment=energy_transfer`, 4420-4425, without `--show-simulation`; dispatch 23674-23675) | trace `energy_transfer_trace.csv` or `--energy-transfer-trace` (4690; 18777-18782); summary `energy_transfer_runs_<n>walls.csv` or `--energy-transfer-summary` (4682; 18570-18576); `00_COMMAND.md` next to the trace (18754-18762); piston/divider event log if `HD_PISTON_EVENTS` names a file (18929-18932; header `t_sigma,kind,u_wall,v_before,v_after,dE,dp` at 1811 for gen3 and hspist3/edmd_core/edmd.c:1273 for gen2) | `g_energy_transfer_dir` = `<dirname argv[0]>/experiments_energy_transfer` (3419), or cwd-relative `experiments_energy_transfer` when argv[0] has no `/` (3401, 3404-3406); **created unconditionally** (`mkdir_p(g_energy_transfer_dir);` 18748) | trace header begins `Time,Wall_X,SpringF,SpringE,EnergyTransferred,Left_Count,Right_Count,PistonWork,PistonWorkDeltaMax,SpringE_Max,` (18802) |
| Simple box (preset 4426-4429) | `simple_box_trace.csv`, `simple_box_state.csv`, `simple_box_step.csv` (20722-20728), `01_PLOT_COMMANDS.md` (20440), `run_params.json` (20715) | run dir `<g_simple_box_dir>/<label_>simulation_<stamp>` (20453-20478; `g_simple_box_dir` per argv[0], 3400/3418) or `--simple-box-run-dir` (4320); `g_simple_box_dir` itself is created and gets `00_COMMAND.md` unconditionally (20652-20653) | headers 20737-20749; "state" rows are summaries (20741-20744), not particles |
| Szilard (4433-4438) | various | `<dirname argv[0]>/experiments_szilard` (3420; created 21139) | not itemised |
| Window only: `--gui-shot=STEP:FILE.bmp` or `limit:FILE.bmp` (5238-5248), `--demo-shot` (4736), key S in `--demo` (12151-12153, default `demo_screenshot.bmp` in the cwd) | BMP | given path | `SDL_RenderReadPixels` + `SDL_SaveBMP` (21540, 2989-2990) |
| `HD_STOP_SNAPSHOT=path` (environment) | one CSV | given path | see state dumps below |
| stdout only | `[EDMD3-RUN]` with build and target (17409-17420); `[EDMD3-HEALTH]` with event hash, build, target (17437-17451); `[GUI-STATS]` (22994-22998) | stdout | one line each |

**Main formats (verbatim)**

- `energy_log.csv`: header `KE_particles,KE_wall,KE_total,T,S_v,F` (21955; 23109); rows every 100 steps (22481-22484) as `"KE_particles: %.8e, KE_wall: %.8e, KE_total: %.8e, T: %.8e, S_v: %.8e, F: %.8e\n"` (10216). Rows are `key: value` pairs that do not match the header; no time or step column. The scratch file PHYSICALINFORMATIONDYNAMICS/_m0_scratch/view_run/energy_log.csv shows exactly this (header plus 4589 `KE_particles: ...` rows).
- `wall_position.csv`: with `bool log_packing_fraction = 1;` (157) the header is `Time, Wall_X, Displacement(σ), Left_Count, Right_Count, Packing_UNITLESS, Packing_ACTUAL_BOUND, Packing_ACTUAL_WITH_WALLBUFFER` (16536), otherwise the first five (21953). Rows only after the divider is released (22827): per EDMD step if `--output-dt` <= 1 (22845-22847), else every `--output-dt` (22850); TIME/RK4 per substep through `g_substep_log` (21932; comment 22832-22836). Row format `"%.3f, %.3f, %.3f, %d, %d, %.6f, %.6f, %.6f\n"` (22880) or `"%.3f, %.3f, %.3f, %d, %d\n"` (22862): time and Wall_X with 3 decimals.

**Full state dump?** Partial, gated forms only:

- `--gen3-dump-initial=PATH` (5157-5158): text `gen3-initial v1`, N, N_L, boxW, boxH, radius, divider, then N lines `%a %a %a %a` (x, y, vx, vy) (17395-17401). Initial state only (written at 18083-18084, 18955-18956); accepted only with the speed-of-sound experiment or the energy-transfer held path (23612-23618).
- `--gen3-checkpoint=hold|record:K:PATH` (5167-5175): binary file with the `G3RunCkpt` struct (17571-17584; magic `G3_RUN_CKPT_MAGIC "g3-run-ckpt-v3"`, 17585), the bytes of the run's output files, then `edmd3_checkpoint_write(g_edmd3, f)` (17630-17633), which writes the EDMD3 struct copy, the disk array, cell arrays and the event heap (hspist3/edmd_core/edmd_gen3.c:1410-1437). Restart refuses another build (17659; hspist3/edmd_core/edmd_gen3.h:193-196) and is wired to the speed-of-sound experiment only (23623-23624).
- `HD_STOP_SNAPSHOT=path` (environment, 2526): once per run, when a piston step reaches its target (7142-7147), CSV header `# t_sigma=%.9f piston_x_sigma=%.9f divider_x_sigma=%.9f N=%d` and `i,x_sigma,y_sigma,vx,vy` (2535-2539), then one row per disk, positions divided by PIXELS_PER_SIGMA, velocities as stored (2540-2545); read through `edmd_backend_particles` (2531), so gen2 or gen3.
- Test-harness format `gen3-state v1`: N, box, radius, cell_px, t0, dividers, pistons, particles `%a %a %a %a`, advance targets, body changes, expected hash and event count (hspist3/edmd_core/tests/gen3_state_io.h:5-11, 34-53); written by gen3_m1_harness.c:318 and gen3_m2_harness.c:325; read by `--engine-replay` (00ALLINONE.c:23344, 23353), which prints `[EDMD3-REPLAY]` to stdout (23435) and writes no file (no `fopen` with write mode in 23344-23450).
- `--gen3-snapshots`: positions only (above).
- Otherwise summaries. The interactive scene writes no per-particle state (no per-particle writer among its `fopen` sites 21713-21925).

**Files written into the cwd (correction to the task premise)**

- Unconditional: only `00_COMMAND.md` (23500).
- `energy_log.csv` and `wall_position.csv`: only `simulation_loop` (21713-21714) and `--single-test` (23100-23101); these are the only `fopen` sites for those names (full `fopen(` list of 00ALLINONE.c checked).
- `run_params.json`: the scene path always (21944), other paths as listed above.
- Tracked in git: hspist3/energy_log.csv, hspist3/wall_position.csv, hspist3/Roboto-Regular.ttf (`git ls-files`); hspist3/00_COMMAND.md and hspist3/run_params.json are not tracked; 1965 tracked `00_COMMAND.md` and 962 tracked `run_params.json` exist elsewhere, mostly under experiments_* folders.

**Runs**

The verifier did not run the binary (assignment: read-only, no binary runs). The drafting agent reported these runs, all with cwd PHYSICALINFORMATIONDYNAMICS/_m0_scratch/view_run and the form `SDL_VIDEODRIVER=dummy perl -e 'alarm 60; exec @ARGV' -- B ...`, B = absolute path of hspist3/experiments_gen3_d52_261014/bin_a71f394/00ALLINONE:

1. `--version`: `00ALLINONE  git a71f394  target mac-O3-e0pre`, `CFLAGS: -O3 -ffp-contract=off`, rc 0. Consistent with the format at 4095 and `-DBUILD_GIT="a71f394" -DBUILD_TARGET="mac-O3-e0pre" -DBUILD_CFLAGS="-O3 -ffp-contract=off"` in bin_a71f394/build_record.txt:30.
2. `--help` (three invocations): 161 lines, none matching `--gen3`, `--engine`, `--gui-`. Consistent with the code: `print_cli_usage` (3907-4069) emits 161 newlines and contains none of those flags.
3. `--show-simulation --headless --mode=edmd --engine=gen3 --seed=1 --gui-run-seconds=2 --gui-stats`, rc 0; reported `[GUI-STATS] engine=gen3 window=0 N=100 frames=0 steps=458908 ...`. Files read by the verifier in view_run: `energy_log.csv` 4590 lines (header + 4589 = 458908/100 rows, key: value form, S_v 0 in every row), `run_params.json` 421 B without seed (pixels_per_sigma 24, L0 20, height 10), `wall_position_run_params.json` byte-identical to it (`cmp`), `wall_position.csv` header only (no release).
4. `--show-simulation --mode=edmd --engine=gen3 --seed=1 --gui-run-seconds=1` (window): reported software-renderer fallback and `Failed to load font: ...`, rc 1. The verifier read view_run/00_COMMAND.md (389 B): it records this command, consistent with only 00_COMMAND.md being rewritten (font failure returns at 23523, before any run_params write).

hspist3/experiments_gen3_d52_261014/bin_a71f394/ holds four files (00ALLINONE, 00ALLINONE.Wall.txt, 00ALLINONE.warnings.txt, build_record.txt), listed by the verifier. `git status` was not re-checked by the verifier.

### E2

**Who writes 00_COMMAND.md**

- `static void write_command_md_to_dir(const char *dir_path, int argc, char **argv) {` (00ALLINONE.c:3186). It builds `<dir>/00_COMMAND.md` (3190) and opens it with `fopen(out_path, "w")` (3192), so it overwrites.
- Content (3212-3220), nothing else: `# Command`; `- Timestamp: %s` (localtime, `"%Y-%m-%d %H:%M:%S %Z"`, 3204); ``- CWD: `%s` `` (getcwd, 3208); a `sh` code block with argv, each argument shell-escaped (fprint_shell_escaped, 3166-3184).
- It carries no BUILD_GIT, BUILD_TARGET or BUILD_CFLAGS, and a seed only when the seed is in argv.

| call site | directory |
|---|---|
| 23500 `write_command_md_to_dir(".", argc, argv);` | the cwd, on every run that gets past 23494-23496: `--version` (4094-4096) and `--help` (4098-4100) exit inside parse_cli_options; `--engine-replay` (23495) and `--gen3-psi6-snapshots` (23496) return first |
| 17738-17739 | speed of sound. The mode0 folder `experiments_speed_of_sound/<EDMD\|TIME>/mode0_real_units/simulation_<dd_mm_yy_HH_MM_SS>` is relative to the cwd (17717-17725) and is always created (mkdir_p at 17734), also when `--speed-sound-run-dir` is given. The mode1 folder is `--speed-sound-run-dir` resolved against the launch cwd (17726-17727), else `.../mode1_normalized_units/<stamp>` (17729-17730) |
| 18762 | energy transfer: the directory of `--energy-transfer-trace`, else `g_energy_transfer_dir` (18755-18761); `g_energy_transfer_dir` is created first in either case (18748) |
| 20653, 20705 | simple box: `g_simple_box_dir` (created at 20652) and the run folder (20703-20705) |
| 3365, 21140 | Szilard: the run folder (3363-3364) and `g_szilard_dir` (created at 21139) |

- Companion files, in the speed-of-sound folders only (17740-17741): `01_PLOT_COMMANDS.md` (3228) and `00_COMMAND_AND_PLOT_COMMAND.md` (3276). Under gen3 the latter adds the trace-column notes and names the `[EDMD3-HEALTH]` run record (3281-3289).
- `run_params.json`: writer `write_run_params_json_to` (2275-2319). `write_run_params_json()` writes it into the cwd (2320-2322) and is called from initialisation paths (6925, 17245, 17393, 23112, 21944). Simple box and Szilard write it into their run folders (20714-20715, 21939-21942). Fields (2282-2317): mode, particles, L0, height, num_walls, wall_mass_factor, wall_thickness, wall_thickness_vis, temperature, optional temperature_segments (2298-2305), kbt1, kB_effective, kBT, dt, output_dt, time_span, units, radius_sigma, pixels_per_sigma. No seed, no build.
- `speed_of_sound_batch_status.json` (16939-16957): complete, `base_seed` (= cli_seed, 16957), requested/valid/invalid/in_progress runs, failure_csv. It is written to a `.tmp` file first (16940-16941).
- Where outputs land:
  - `experiments_simple_box`, `experiments_energy_transfer` and `experiments_szilard`: the path strings are set relative to argv[0]'s directory (init_output_dirs_from_argv0, 3395-3422; 3418-3420). They are relative to the cwd when argv[0] has no slash or its directory is "." (3404-3408, 3413-3416). The function only sets strings; the folders are created by the modes above (18748, 20652, 21139).
  - CLI output paths are resolved against the launch cwd (resolve_cli_output_path, 3424-3433; g_launch_cwd from getcwd at 3396).
- Seed and build in data outputs:
  - the energy-transfer summary CSV `seed,build_git,build_target,build_cflags` (19641-19642, 19831);
  - stdout `[EDMD3-RUN]` and `[EDMD3-HEALTH]` with `build=` and `hash=` (17409-17410, 17437-17451); the run id carries the seed, e.g. `[EDMD3-RUN] L0=78.5 M=50 run=0 seed=974638985:` (experiments_gen3_d55_261015/tests/g5_T1/B_N0100_eta0.050000_m_50_r0_old/stdout.log:39);
  - the gen3 checkpoint's build, target and seed (17592-17595).

**Run-manifest conventions of the experiment folders.** These are written by scripts, not by the binary.

- **Per binary:** `bin_<hash>/build_record.txt` with the commit, date, CPU, compiler, the full cc line, warning counts, and SHA-256 of sources and binaries (build_clean.sh:15-42; build_d52.sh:28-44, sources at :38).
- **Per run:**
  - `command.txt`, `version.txt` (the binary's `--version`) and `run.log`, from cluster/resched_gate_261005/audit_runs_261007.py `run()` (:95, :92, :100). Each run gets its own folder as cwd (:100); an existing folder is refused (:89).
  - `cmd.txt` and `stdout.log`, from experiments_gen3_d50_261014/scripts/tests_d50.py (:69, :71). experiments_gen3_d50_261014/README.md:29: "Committed: each run's cmd.txt and stdout.log (added with -f, since *.log is ignored) and each test's RESULTS.md."
  - `<name>.cmd` with the date, the binary with the first 16 hex digits of its SHA-256, the flags (experiments_gen3_d52_261014/scripts/run_scene_d52.sh:16), then `# exit code N after S s` (:21).
  - On KOA, per cell: `.sha256_uncompressed`, the SHA-256 of each uncompressed file before `gzip -9` (per comment at cluster/gen3_stage1_261012/s1_worker_koa.sh:12).
- **Per folder:**
  - `*SHA256_FILES.txt` as `sha256  size  relpath`, e.g. experiments_gen3_d34s_261012/D34S_SHA256_FILES.txt:1 `f3bd70b0...      675672  bin_406044f/00ALLINONE` (311 lines);
  - `SHA256_UNCOMPRESSED.txt` (experiments_gen3_m3_261009/evidence_f42befb/SHA256_UNCOMPRESSED.txt:1 `# SHA-256 of each raw file before its lossless compression (gzip -9), ...`);
  - `RUN_RECORDS.txt`, the grep of the run-record lines, because the logs are ignored (evidence_f42befb/RUN_RECORDS.txt:1-2);
  - `README.md` at the top of 5 of the 44 `experiments_gen3_*` folders on disk: d48, d50, d52 and d55 tracked; d56 untracked, written 2026-10-15 13:35 by the parallel HardDisks work. Also tracked: experiments_gen3_s1_261012/full_261015/README.md.
  - Older folders have a hand-written `00_COMMAND.md` (experiments_gen3_m1_261008/, experiments_gen3_m2_261009/), the same file name the binary writes. The repository tracks 1965 files named 00_COMMAND.md.
- **Seeds:** `def run_seed(base, l, m, r):` (hspist3/validation/tests_20260913.py:135; docstring :136 "Exact replica of speed_sound_run_seed() in 00ALLINONE.c."). The Stage 1 task lines carry the exact seed and the base: cluster/gen3_stage1_261012/tasks_stage1_koa_261012.txt:1 `B stage1/gen3/B_N0100_eta0.050000/m_50 50 0 974638985 ... 20261311 ...`. The field list is at s1_worker_koa.sh:6.
- **Cluster build guard:** a `.build_git` file per cell folder holds the `--version` line, and a different binary is refused (cluster/gen3_partA_261014/pa_worker_koa.sh:31-47; :32 "conf_worker.sh's build guard, unchanged").

**Examples**

1. experiments_gen3_d52_261014/tests/c_T1/B_N0100_eta0.050000_m_500_r0_new/. cmd.txt and stdout.log are tracked; the rest is untracked.
   - 00_COMMAND.md:3 `- Timestamp: 2026-10-14 22:17:52 HST`; :4 the CWD, which is that folder in HardDisks_gen3.
   - :7 is the argv, with `--seed=20261311 ... --speed-sound-exact-seed=2749684526 --engine=gen3 --gen3-seeding=lattice ...`.
   - The folder also holds run_params.json (no seed field), speed_of_sound_batch_status.json (`"base_seed": 20261311`) and an `experiments_speed_of_sound/` subfolder (the mode0 folder, 17734).
2. experiments_gen3_d55_261015/e0/default/ctrl_min/plain/, a rule-4 run. command.txt and version.txt are tracked.
   - The binary also wrote 00_COMMAND.md (:3 `- Timestamp: 2026-10-15 11:25:44 HST`), 00_COMMAND_AND_PLOT_COMMAND.md, 01_PLOT_COMMANDS.md, run_params.json and speed_of_sound_batch_status.json.
   - It created `experiments_speed_of_sound/EDMD/mode0_real_units/simulation_15_10_26_11_25_44/`, although command.txt carries `--speed-sound-run-dir=` with that folder.
3. The tracked hspist3/experiments_energy_transfer/00_COMMAND.md is modified in the working tree (` M`).
   - `git diff` changes :3 from `2026-09-18 13:28:33 HST` to `2026-10-10 12:10:05 HST` and changes the argv at :7 (`./00ALLINONE --mode=edmd --experiment=energy_transfer ...`); the CWD (:4) is hspist3.
   - hspist3/energy_log.csv (1000 lines changed) and experiments_energy_transfer/energy_transfer_trace.csv changed with it.
   - This is what a run of the binary from hspist3/ does.

## F. Regression and the do-not-touch set

### F1

**Stage 1 manifests** (hspist3/cluster/gen3_stage1_261012/). Paths are relative to hspist3/: make_s1_manifests.py:2-3 says "sha256sum -c format, paths relative to hspist3/", and HS is set at :12.

| manifest | entries | files |
|---|---|---|
| s1_inputs.sha256 | 6 | `cluster/gen3_stage1_261012/tasks_stage1_koa_261012.txt`, `.../tasks_stage1_koa_261012.sha256`, `.../chunks_stage1_koa_261012.tsv`, `.../build_record_f72d764.txt`, `.../s1_worker_koa.sh`, `cluster/koa_env.sh` (:1-6) |
| s1_merge_inputs.sha256 | 15 | Stage 1 scripts (3): `s1_chunk_logs_koa.sh`, `merge_s1_koa.py`, `reduce_s1_koa.py`. Reducers (2): `cluster/confinement_20261013/reduce_B.py`, `cluster/gen3_testR_261011/reduce_testR_koa.py`. Top level (1): `plot_speed_of_sound_edmd.py`. validation/ (9): `edmd3_health_guard.py`, `edmd_acc_guard.py`, `paper1_confinement_prereg_20261012.py`, `paper1_confinement_results_261004.py`, `paper1_populate_cs_err_20261002.py`, `provenance_edmd_acc_261009.py`, `resched_testT_261007.py`, `resched_testT_followup_261007.py`, `tests_20260913.py` (:1-15) |
| tasks_stage1_koa_261012.sha256 | 1 bare hash, no file name | `e3e3e945...0b24` of the 16,550-line task list. It cannot be checked with `-c`; s1_koa.sbatch:41 and s1_merge.sbatch:39 compare it. |

- Maker: make_s1_manifests.py:35-38. It writes list `a`, and list `b` extended by an import probe of reduce_s1_koa.py, tests_20260913 and paper1_populate_cs_err_20261002 (:19-30). It stops without build_record_f72d764.txt (:34). Usage (:9): `python3 cluster/gen3_stage1_261012/make_s1_manifests.py`, from hspist3/.
- The job scripts that do the checking (s1_koa.sbatch, s1_merge.sbatch) are themselves in no manifest.

**Engine sources in these manifests: none.**

- No entry for 00ALLINONE.c, edmd_core/*, the Makefile, experiment_validation.* or kissfft.
- Search: all 64 `*.sha256` and `*SHA256*` files under hspist3/cluster/ and hspist3/experiments_gen3_* (src/ copies excluded) for `00ALLINONE\.c|edmd_core/|Makefile|experiment_validation|kiss_fft`: 0 files match.
- SHA-256 lists of the engine sources exist only in Mac `bin_*/build_record.txt` files (10 found, e.g. experiments_gen3_d52_261014/bin_00f4dd0e-dirty/build_record.txt:10-21).

**How the Stage 1 engine is pinned instead: by commit and by binary SHA-256**

- build_record_f72d764.txt:
  - :2 `commit f72d764a8a7e837644c26bd8abc175d9e2cac268`;
  - :3 `version 00ALLINONE  git f72d764  target koa`; :4 `version_ld 00ALLINONE  git f72d764  target koa-ld`;
  - :5 `sha256 0bb1929e92424242ecbcd57d1dec183ae08ff055fa4c29c21357dc5058348d58 00ALLINONE`;
  - :6-8 the hashes of `bin/00ALLINONE_ld`, `bin/gen3_m1` and `bin/gen3_m2`;
  - :9 `compiler gcc (GCC) 14.3.0`; :10 `built 2026-10-13 07:03:35 UTC on cn-03-33-01`.
- Checked at every job start (s1_koa.sbatch):
  - :39 `sha256sum --status -c "$P/s1_inputs.sha256"`;
  - :41 the task-list hash;
  - :43-50 the binary's SHA-256 against the record and against `~/harddisks_gen3_s1/hspist3/logs/BUILD_RECORD_S1.txt`, and the `--version` line.
- The merge job checks both manifests (s1_merge.sbatch:34-35).
- The build itself (cluster/gen3_s1_gate_261012/build_s1_koa.sh):
  - :19 HEAD must be `f72d764a8a7e837644c26bd8abc175d9e2cac268`;
  - :20 it refuses when a binary or the record already exists;
  - :21 no modified tracked files;
  - :22 kissfft at `febd4caeed32e33ad8b2e0bb5ea77542c40f18ec`;
  - :27 and :32 run `make -B koa-ld` and `make -B koa`.

**f72d764 against HEAD (engine paths)**

- `git diff --stat f72d764 HEAD -- hspist3/00ALLINONE.c hspist3/edmd_core hspist3/experiment_validation.c hspist3/experiment_validation.h hspist3/Makefile hspist3/kissfft`: 00ALLINONE.c 1978 lines, edmd_gen3.c 25, edmd_gen3.h 11 (1643 insertions, 371 deletions). edmd.c, edmd.h, edmd_accelerated.*, the Makefile, experiment_validation.* and the kissfft gitlink are unchanged.
- f72d764 is an ancestor of HEAD.
- The Stage 1 binary was therefore built from the f72d764 tree, not from the working tree.

**Part A** (hspist3/cluster/gen3_partA_261014/; tracked files: README.md, pa_orient.py, pa_orient_261014.tsv, pa_reduce.py, pa_tasks.py, pa_worker_koa.sh)

- Manifest and maker: **not found.**
  - There is no `*.sha256` file in the folder.
  - `grep -n "sha256\|SHA256\|shasum\|hashlib\|manifest"` over its files: 0 hits.
- The Part A engine pin is not registered yet:
  - Registration §4.8.5 "Engine, health and exclusion" (0000_PLAN_OVERALL/ALL_MARKDOWNS/261012_paper1_confinement.md:18186), "One build" at :18188: "gen3 at one campaign commit whose build inputs (every file compiled into the binary) are identical to the accepted `f72d764`".
  - cluster/gen3_partA_261014/README.md:4: "Nothing here runs on KOA before amendment A5 (decision 50 item 2). The binary needs commit 8f538a21 or later". `git diff --stat f72d764 8f538a21` on the engine paths changes only 00ALLINONE.c (1780 lines).
  - Notes :29338: "The orientation pin of Report D53-2 has since been committed locally on engine-gen3 as a102d011. That is the commit A5's KOA build names (decision 53 item 2)." a102d011 is an ancestor of HEAD. HEAD's engine differs from it in 00ALLINONE.c (202 lines), edmd_gen3.c (25) and edmd_gen3.h (11).
  - A5 is still to be written by the plan author (notes :31860).
- The only Part A guard today is the per-folder `.build_git` check (pa_worker_koa.sh:31-47).

**Check against the working tree** (from hspist3/, read-only)

- `shasum -a 256 -c cluster/gen3_stage1_261012/s1_inputs.sha256`: 6 OK, 0 FAILED.
- `shasum -a 256 -c cluster/gen3_stage1_261012/s1_merge_inputs.sha256`: 15 OK, 0 FAILED.
- `shasum -a 256 cluster/gen3_stage1_261012/tasks_stage1_koa_261012.txt` = `e3e3e945ca1e6015a505faf1632efd1081918294f0b35863a592298b2f1c0b24`, equal to the .sha256 file.
- The consolidation recorded the same count of 21 OK (0000_PLAN_OVERALL/consolidation_d51_261014/d55_consolidation_log_261015.txt:48-49).

**The frozen set, by what the files actually check**

- (a) the 21 manifest entries and the task-list hash;
- (b) the Stage 1 KOA binaries by SHA-256, i.e. the f72d764 tree with kissfft febd4ca and no local change (build_s1_koa.sh:19-22);
- (c) Part A: a102d011 is named, A5 is open, and no manifest exists;
- (d) the engine sources at HEAD only by policy (brief rule 1). No manifest covers them. Mac regressions compare against bin_a71f394, built from a71f394d, whose engine sources equal HEAD's (A1; F2 T1).

### F2

**edmd_core/tests/** (build lines copied from each file; every one builds from hspist3/; output to a caller-chosen `<out>`)

| file | what it tests | build (verbatim) | run | last recorded time |
|---|---|---|---|---|
| gen3_m1_harness.c | A: gen3 with schedule and contact audits; B: no audits, equal hash; C: repeat; D: gen2; speed; diverge (:2-13) | `cc -std=c11 -O3 -ffp-contract=off -Wall -Wextra -o <out> edmd_core/tests/gen3_m1_harness.c edmd_core/edmd_gen3.c edmd_core/edmd.c -lm` (:15-16) | `<out> audit\|speed\|diverge [--quick]` (:17), `dump <dir>` (:281); separate processes (:17-18) | `gen3_m1 audit` (144.7 s), experiments_gen3_m2_261009/00_COMMAND.md:37 (sources of ddae96c, :35). speed and diverge: not found |
| gen3_m2_harness.c | M2 cells with dividers and pistons, runs A/B/C/D, ledgers (:1-32) | `cc -std=c11 -O3 -ffp-contract=off -Wall -Wextra -o <out> edmd_core/tests/gen3_m2_harness.c edmd_core/edmd_gen3.c edmd_core/edmd.c -lm` (:34-35) | `<out> audit\|speed [--quick]` (:36), `dump <dir>` (:468). Speed mode writes gen2's `gen3_m2_speed_evlog.csv` into its cwd (00_COMMAND.md:33) | `gen3_m2 audit --quick` (30.3 s) and `gen3_m2 audit` (227.4 s), 00_COMMAND.md:38; `gen3_m2 speed` (32.0 s), :39 |
| gen3_body_rule_test.c | spring-divider contact rule against an independent scan, plus constructed cases (:1-18) | `cc -std=c11 -O2 -ffp-contract=off -Wall -Wextra -o <out> edmd_core/tests/gen3_body_rule_test.c -lm` (:7) | `<out> [n]`, default 20000 (:8) | "`gen3_body_rule_test 20000` (about 10 s)", 00_COMMAND.md:37 |
| gen3_band_edge_test.c | band-edge stress cells under all audits (:1-12) | `cc -std=c11 -O2 -ffp-contract=off -Wall -Wextra -o <out> edmd_core/tests/gen3_band_edge_test.c edmd_core/edmd_gen3.c -lm` (:13) | `<out>`, no arguments (`int main(void){` :55) | no individual time found; part of the harness batch below |
| gen3_band_edge_ties.c | classification of the audit findings (:1-23) | `cc -std=c11 -O2 -ffp-contract=off -Wall -Wextra -o <out> edmd_core/tests/gen3_band_edge_ties.c -lm` (:24) | `<out>` (`int main(void){` :108) | no individual time found |
| gen3_checkpoint_test.c | checkpoint and restart U, C, R, R', R2 bit for bit; 2 refusals (:1-17) | `cc -std=c11 -O3 -ffp-contract=off -Wall -Wextra -o <out> edmd_core/tests/gen3_checkpoint_test.c edmd_core/edmd_gen3.c edmd_core/edmd.c -lm` (:18-19) | `<out>`, exit 0 iff all equal and both refusals refuse (:20) | `done H2 default: exit 0, 394 s` (experiments_gen3_h_261009/evidence_ab80304_runner.txt:48). This is an upper bound: the runner takes the time when it reaps the process (stageH_evidence.py:70-73) and reaps in start order (:75-77), with up to 10 children at once (:31 `JOBS = 10`). Verdict: evidence_ab80304/h2/default/checkpoint_test_output.txt:101 `VERDICT: PASS -- 329 comparisons, 0 different; refusals both refused`. Not built by build_clean.sh (36-41); built in stage H (build_record_ab80304.txt:57) and by d18_build.sh:20 |

- The harness batch (experiments_gen3_m3_261009/stageA_evidence.py:100-106) runs m1 audit, m2 audit --quick, m2 audit, body_rule 20000, band_edge_test and band_edge_ties in sequence.
  - Run: `python3 experiments_gen3_m3_261009/stageA_evidence.py --bin-dir <frozen binaries> --out <evidence dir> --only harness` (usage :19-20, "from hspist3/ of the engine-gen3 worktree").
  - Time: `done harness: rc 0, 408 s` (evidence_f42befb_runner.txt:40), with at most 12 processes at once (:3, :31 `MAXP = 12`).
  - The outputs are compared byte for byte with the committed experiments_gen3_m1_261008/m1_audit_output.txt and experiments_gen3_m2_261009/m2_audit*_output.txt: IDENTICAL (stageA_tables_output.txt:48-50).
- There is no `make test` (Makefile:63) and no single committed runner; each stage has its own *_evidence.py.

**Rule 4: gen2 byte identity**

- Source text, numbered: 0000_PLAN_OVERALL/ALL_MARKDOWNS/261012_paper1_confinement.md:8098-8101, under "GROUND RULES" (:8089): "4. gen2 stays frozen: after every stage that touches 00ALLINONE.c, the Makefile or edmd.c, the default build and --engine=gen2 must be byte-identical to 7b08827 on ctrl_min and ctrl_leg. If not: stop all code stages and report. Do not explain it away." The summary bullet is at :8043.
- Runner: hspist3/cluster/resched_gate_261005/audit_runs_261007.py.
  - Usage (:32-34), from hspist3/: `python3 cluster/resched_gate_261005/audit_runs_261007.py run --bin <new> [--ref-bin <ref>] --out <dir> [--jobs 8] [--cases a,b]`, then `... report --out <dir>`.
  - Variants per case: audit, plain, and ref with `--ref-bin` (:87). Each run uses its folder as cwd (:100). An existing folder is refused (:89).
  - Compared files: `wall_x_positions_*_run0.csv` and `speed_of_sound_psi6.csv`, plus `ev_/tr_/red_9700.csv` when present (outputs(), :112-114), each by `cmp -s` (same(), :117-120: `return "IDENTICAL" if all(subprocess.run(["cmp", "-s", x, y]).returncode == 0 for x, y in zip(f1, f2)) else "**DIFFERENT**"`).
- The exact commands of the last run (experiments_gen3_d55_261015/e0/REPORT_rule4.txt:2; consolidation log :88):
  - `python3 cluster/resched_gate_261005/audit_runs_261007.py run --bin <the new build> --out hspist3/experiments_gen3_d55_261015/e0/default --cases ctrl_min,ctrl_leg --jobs 4`;
  - the same with `--bin experiments_gen3_d55_261015/bin_engine_gen2.sh`, which is `exec .../experiments_gen3_d50_261014/bin_0544b7e8-dirty/00ALLINONE "$@" --engine=gen2` (bin_engine_gen2.sh:3);
  - the 7b08827 references were copied into `<case>/ref` from experiments_gen3_m3_261009/evidence_f42befb/e0/<tag>/<case>/ref (d55 README.md:12-14); then the report was run.
  - The reference's version.txt: `00ALLINONE  git 7b08827  target mac-O3-gen3count`.
- Result: the plain-vs-ref column is IDENTICAL for ctrl_min and ctrl_leg in both builds (REPORT_rule4.txt:9-10 default, :61-62 gen2flag).
  - The summary line `E0 (plain vs ref IDENTICAL in ctrl_min, ctrl_leg): FAIL (or not run)` (REPORT_rule4.txt:21, :73) comes from `e0_ok = len(e0) == 4 and ...` over ctrl_min, ctrl_leg, afix and afix_leg (audit_runs_261007.py:151-152). Only 2 cases ran.
- Time:
  - m3: `done e0 default: rc 0, 7 s` and `done e0 gen2flag: rc 0, 7 s` (evidence_f42befb_runner.txt:7-8); `wall time 6 s` (evidence_f42befb/e0/run_default.txt:5).
  - d55: no recorded time. File mtimes give command.txt at 11:25:44 and the last run.log at 11:25:52 (about 8 s; mtimes, not a record).
- Reference storage:
  - In the m3 references, `ref/wall_x_positions_L0_100_wallmassfactor_50_run0.csv` and `ref/run.log` are git-ignored (`!!`; .gitignore:61 `hspist3/experiments_*/**/wall_x_positions_*.csv`, :70 `*.log`); `ref/speed_of_sound_psi6.csv` is tracked.
  - 28 copies of `ctrl_min/ref/wall_x_positions_L0_100_wallmassfactor_50_run0.csv` exist under experiments_gen3_*.
  - Their SHA-256 `0b10c0ec...a590` is in tracked lists, e.g. experiments_gen3_d34s_261012/D34S_SHA256_FILES.txt:88 (path `evidence_406044f/h0/default/ctrl_min/ref/...`).
  - Checked today: m3 ref = m3 plain = d55 plain = d55 ref = `0b10c0ec...a590`.

**Event-hash comparisons**

- The event hash exists for gen3 only:
  - edmd_gen3.h:159 `uint64_t edmd3_event_hash(const EDMD3* S);           /* FNV-1a over (t, type, a, b) of every executed event */`;
  - initialised at edmd_gen3.c:1158, updated at :1283, returned at :1317.
  - gen2 has none: the only "hash" in edmd.c is the comment `grid (uniform spatial hashing)` (:353).
- (i) **Replay through the driver:** experiments_gen3_m3_261009/run_replays.sh.
  - Usage (:8): `bash experiments_gen3_m3_261009/run_replays.sh <binary> <harness dir> <state dir> > replay_output.txt`, from hspist3/. 16 cells x 4 reader/stop variants (:2-7). The harness dir must hold gen3_m1 and gen3_m2 (:13).
  - The `--engine-replay` runs return at 00ALLINONE.c:23495, before the cwd write at 23500.
  - Match test in the driver: 00ALLINONE.c:23433-23434 (`const int match = !strcmp(got, hexs) && ev == expect_events;`).
  - Recorded: `replays: 64, MATCH with clean=1: 64` (evidence_f42befb/replay/replay_output.txt:85); `done replay: rc 0, 120 s` (evidence_f42befb_runner.txt:38).
- (ii) **T1, Part B unchanged:** experiments_gen3_d50_261014/scripts/tests_d50.py (:4-6).
  - Usage (:30): `python3 tests_d50.py <label> <T1|T1AF|T2|C2|C2P|R0|R1|LAT> <NEW binary> [OLD binary] [--seeder=legacy|lattice]`, from hspist3/.
  - "Equal" means the run record's event hash, t_end and the SHA-256 of every output file except stdout.log, cmd.txt, summary.csv, *.failures.csv and *.md (outputs(), :100-103).
  - The last invocation (consolidation log :72): `TESTS_OUTROOT=hspist3/experiments_gen3_d55_261015/tests python3 $SP/d50/tests_d50.py g5_T1 T1 <the new build> experiments_gen3_d52_261014/bin_a71f394/00ALLINONE`.
  - `$SP` is not defined in the log. A copy at /private/tmp/claude-501/-Users-chrisharing-Desktop-CCS-complex-coupled-systems-Repo-HardDisks/91cb08ec-0599-4faa-a6af-d5a7ca834255/scratchpad/d50/tests_d50.py (mtime 2026-10-15 11:22; the consolidation ran from 11:20:32, log :1) differs from the committed file only at :34 (`HS = ".../Repo/HardDisks/hspist3"`). That it is the `$SP` copy is an inference from path and time.
  - Result: EQUAL in 6 of 6 on event hash, t_end and output SHA-256 (experiments_gen3_d55_261015/tests/g5_T1/RESULTS.md:8-13). Seconds old/new: 2/3, 5/5, 2/3, 14/14, 4/4, 9/9 (same lines).
- (iii) **GUI hash:** experiments_gen3_d52_261014/scripts/tests_d52.sh:5-6 (`--gui-keys` action `hash`). Usage (:9): `tests_d52.sh <binary> <prefix>`.
  - Print site: 00ALLINONE.c:21345 `[GUI-HASH] step %ld: hash=%016llx t=%.17g events=%ld`.
  - Result: `verdict: EQUAL (the four runs at the three steps)` (tests/c_RESULTS.md:12, :23).
  - Time: `# exit code 0 after 5 s` (tests/c_H_e030_closed.cmd:4); 4 s (c_H_e069_closed.cmd:4).
- (iv) **Inside the harnesses:** runs A, B and C must give equal hashes (gen3_m1_harness.c:3-6, gen3_m2_harness.c:4-8).

**Hard-coded paths in the committed test tooling**

- tests_d50.py:34 `HS = "/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks_gen3/hspist3"`, with the default OUTROOT (:35), TASKS (:36) and the prereg module (:39) under it.
- tests_d52.sh:12-13 and run_scene_d52.sh:9-10 (HardDisks_gen3 and a session scratchpad). The scene runs from a scratch folder with the font symlinked from HardDisks_gen3 (run_scene_d52.sh:15, :17).
- build_clean.sh:9 and the kissfft copy at :14 (A1: that kissfft folder is empty now).
- experiments_gen3_m3_261009/bin_stageA/bin_engine_gen2.sh:2.
- stageA_evidence.py:30 (E0REF under /private/tmp), and :94-96, which writes bin_engine_gen2.sh into the binary folder.
- `git worktree list`: `/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks_gen3  a75e4cd4 [engine-gen3] locked`.

## G. One bounded compile check

Run by the main session, 2026-10-15 13:06 HST. The three object compiles took 1.53 s in total (timings below), and the bounded harness run took 60 s. Everything was written into `PHYSICALINFORMATIONDYNAMICS/_m0_scratch/g/`, compiled from `hspist3/` as the cwd with `-o` paths into that folder. The folder was deleted afterwards (see "Writes of this session"). The output below is copied verbatim from the files the script wrote there; long lines are unchanged.

**Answer.** `edmd_core` is neither header-only nor `#include`d into the main unit (A2). `edmd.c`, `edmd_accelerated.c` and `edmd_gen3.c` compile as separate objects with the release flags and zero warnings under `-Wall -Wextra`. Linked together with an empty `main`, they need nothing beyond libc/libm. Without any `main`, the only undefined symbol is `_main`. **The core expects no global or function from `00ALLINONE.c`.**

**Compiler and flags.** These are the `release` flags of `hspist3/Makefile:40` with `$(INCLUDE)` from pkg-config (`Makefile:22`). The provenance `-D`s of `Makefile:87` are left out: they are used only in `00ALLINONE.c`.

```
cc: Apple clang version 17.0.0 (clang-1700.4.4.1)
flags: -O3 -march=native -ffp-contract=off -I/opt/homebrew/include/SDL2 -I/opt/homebrew/include -D_THREAD_SAFE -I/opt/homebrew/Cellar/harfbuzz/12.2.0_1/include/harfbuzz -I/opt/homebrew/Cellar/glib/2.86.3/include/glib-2.0 -I/opt/homebrew/Cellar/glib/2.86.3/lib/glib-2.0/include -I/opt/homebrew/opt/gettext/include -I/opt/homebrew/Cellar/pcre2/10.47_1/include -I/opt/homebrew/Cellar/graphite2/1.3.14/include -I/opt/homebrew/opt/freetype/include/freetype2 -I/opt/homebrew/opt/libpng/include/libpng16
```

Each object was compiled as `cc -c <flags> -Wall -Wextra -o _m0_scratch/g/<name>.o edmd_core/<name>.c`:

```
compile edmd_core/edmd.c rc=0 time=.64s warnings=0
compile edmd_core/edmd_accelerated.c rc=0 time=.32s warnings=0
compile edmd_core/edmd_gen3.c rc=0 time=.57s warnings=0
```

**Undefined symbols per object (`nm -u`), verbatim:**

```
[edmd.o]
___sincos_stret ___stack_chk_fail ___stack_chk_guard ___stderrp _asin _atan2 _calloc _cos _fclose _fflush _fmod _fopen _fprintf _free _fwrite _getenv _hypot _log _malloc _memcpy _printf _rand _realloc _snprintf 
[edmd_accelerated.o]
___sincos_stret ___stack_chk_fail ___stack_chk_guard ___stderrp _calloc _cos _fprintf _free _hypot _log _malloc _memcpy _rand _realloc 
[edmd_gen3.o]
___sincos_stret ___stack_chk_fail ___stack_chk_guard ___stderrp _acos _atan2 _bzero _calloc _exit _fmod _fprintf _fread _free _fwrite _hypot _ldexp _malloc _memcpy _memset _memset_pattern16 _nextafter _printf _qsort _realloc _snprintf 
```

All of these are libc, libm or compiler runtime symbols. Three of them matter for reuse:
- `_rand` in `edmd.o` and `edmd_accelerated.o`: the heat bath draws from libc `rand()`, at `edmd.c:1082-1083, 1092, 1107, 1159` and `edmd_accelerated.c:720-721, 730, 745, 797`.
- `_getenv` in `edmd.o`: `HD_CONTACT_AUDIT`, read once per process (`edmd.c:1750`).
- `_exit` in `edmd_gen3.o`: `exit(3)` on heap out-of-memory (`edmd_gen3.c:197`).

Defined global symbols (`nm -gU`): `edmd.o` 51, `edmd_accelerated.o` 30, `edmd_gen3.o` 46. The accelerated object does not collide with gen2 at link time because it renames its API to `edmd_acc_*` by `#define` (`edmd_accelerated.c:4-33`).

**Links:**

```
--- link the three objects with a stub main (int main(void){return 0;})
stub link rc=0
--- link without any main (shows only what the objects need from outside)
no-main link rc=1
Undefined symbols for architecture arm64:
  "_main", referenced from:
      <initial-undefines>
ld: symbol(s) not found for architecture arm64
clang: error: linker command failed with exit code 1 (use -v to see invocation)
```

**The two standalone harnesses, built by the lines in their own headers** (`edmd_core/tests/gen3_m1_harness.c:15-16`, `gen3_m2_harness.c:34-35`):

```
cc -std=c11 -O3 -ffp-contract=off -Wall -Wextra -o _m0_scratch/g/gen3_m1_harness edmd_core/tests/gen3_m1_harness.c edmd_core/edmd_gen3.c edmd_core/edmd.c -lm
cc -std=c11 -O3 -ffp-contract=off -Wall -Wextra -o _m0_scratch/g/gen3_m2_harness edmd_core/tests/gen3_m2_harness.c edmd_core/edmd_gen3.c edmd_core/edmd.c -lm
```
```
m1 harness build rc=0 warnings=0
m2 harness build rc=0 warnings=0
```

**One bounded run.** Command: `perl -e 'alarm 60; exec @ARGV' -- ./gen3_m1_harness audit --quick > m1_audit_quick.txt 2>&1`, with cwd `_m0_scratch/g/`.
- Result: `rc=142 elapsed=60.0s`. The process was **stopped by the 60 s cap**; it did not finish.
- The output file had 139 lines and exactly 8,192 bytes, a multiple of a typical stdio buffer size. Stdout was redirected to a file, so it was block-buffered, and the process may have run further than what was flushed.
- The flushed output holds three complete audit cells (fluid, dense, lattice; `audit_cell` calls at `gen3_m1_harness.c:326-328`). It then holds the fourth cell, tie (`:329`), up to its `[EDMD-HEALTH gen2]` line; it ends in the middle of the tie cell's KE table header, `| engine | relative KE `.
- The solid cell (`:330`) and the later parts of the audit mode were not reached in the flushed output.

The run shows that gen3 and gen2 are driven to completion of whole cells without `00ALLINONE.c`. Verbatim excerpt (first cell header and table, then the four cells' lines):

```
# gen3 M1 harness: audit (quick)

build: Apple LLVM 17.0.0 (clang-1700.4.4.1), double 8 bytes, long double 8 bytes; cell width 32 px, origin shift 8192 units

### Cell fluid: N = 400, eta = 0.3927, box 960.0000 x 480.0000 px (40.0000 x 20.0000 sigma), T = 40 sigma-time, cell width 32 px

| run | event hash | physical events | crossings | stale | final state equal to A |
|---|---|---|---|---|---|
| A gen3, audits on | e06d57c1033cede0 | 34552 | 17773 | 146625 | - |
| B gen3, audits off | e06d57c1033cede0 | 34552 | 17773 | 146625 | yes |
| C gen3, audits off, again | e06d57c1033cede0 | 34552 | 17773 | 146625 | yes |

audits do not steer (A = B): YES; deterministic (B = C): YES
```

| cell (line in the output) | event hash A = B = C | A = B, B = C | schedule audit: missing / extra (pairs, walls, crossings) | gen3 relative KE drift |
|---|---|---|---|---|
| fluid, eta 0.3927 (7) | `e06d57c1033cede0` | YES, YES | 0 / 0 in all three classes | 1.48e-16 |
| dense, eta 0.7000 (42) | `a6bfa7169d42deea` | YES, YES | 0 / 0 | -4.12e-16 |
| lattice, eta 0.7000 (77) | `5c5f0c1958f9f1fc` | YES, YES | 0 / 0 | -2.18e-15 |
| tie, eta 0.4007 (112) | `2a1a7a8c4769809b` | YES, YES | 0 / 0 | not flushed |

Each complete cell also printed `[EDMD3-HEALTH] overlap_repair=0 wall_overdue=0 past_event=0 clamp_repair=0 cell_repair=0 grid_escape=0 stagnation=0 contact_now=0 ...`, with `local_findings=0 full_findings=0`. `heap_max` was 3721, 21522, 21565 and 1840 for the four cells.

The harness also prints a virial Z per cell, for example fluid gen3 2.79117 against gen2 2.81409. Those are 40 σ-time runs of N = 400 in a hard-wall box; the harness itself labels the comparison "information, not a test". They are not a measurement and are not used here.

**What G adds to C4.** The core's link boundary is the header pair `edmd.h` + `edmd_gen3.h`. A PID executable can compile `hspist3/edmd_core/edmd_gen3.c` (and `edmd.c` if gen2 is wanted) unmodified into its own build folder with the same flags. Doing so writes nothing into `hspist3/`, provided every `-o` and every object path points into `PHYSICALINFORMATIONDYNAMICS/`.

## H. Decision: recommended here, Chris decides

Option (d) is added at the planning chat's request (2026-10-15): a fork of `edmd_gen3.c` inside `PHYSICALINFORMATIONDYNAMICS/` with a new segment–disk primitive.

Facts the grades rest on:
- (i) The engines link alone and reference nothing in the driver [A2, C4, G].
- (ii) The engines' only geometry is the rectangle, vertical slabs and vertical faces, hard-coded per side in `static` functions; `EDMD_Params` has no boundary field [B2, C4].
- (iii) The renderer lives in the unit with `main`, with `static` primitives, global state and no frame hook [D1, D2].
- (iv) The binning reads globals [D3].
- (v) gen3's checks (event hash, ledger, Health, brute-force schedule audit, contact audit) are inside `edmd_gen3.c` and are switched on by its API [C6].
- (vi) The registered Stage 1 engine is pinned by commit `f72d764` and by binary SHA-256, not by source manifests [F1].

### Demo 0: growth and viewer, no particles

- **(a) Link `edmd_core` unmodified and reuse the SDL2 renderer through an adapter.** Not feasible as described. Demo 0 needs no engine, so the linking half has nothing to do. The adapter half has no seam to attach to: no draw function takes a state argument (`void render_particles() {`, `00ALLINONE.c:7409`), the primitives are `static` (`gs_line` 1090, `gs_disk` 1109), and no header declares them (D2). All three Demo 0 drawing needs (a second disk population, a polyline/arc set, a field overlay) are "only by edit" (D2 table). The one route that edits nothing is an untested unity build (`#include` of the 23,715-line file with `-Dmain=...`). It would bring the 576 file-scope variables (B3) and `simulation_loop` into the PID process.
- **(b) A new C/SDL2 executable in this folder that copies nothing.** Fits Demo 0.
  - Toolchain: the same as hspist3's (Apple clang 17, SDL2 2.30.12 / SDL2_ttf 2.24.0 via pkg-config; A1), with `-ffp-contract=off` from the start (`Makefile:34-40`).
  - What hspist3 offers as reference only, copied by idea and not by text:
    - the scanline disk fill (`gs_disk`, 1109-1127);
    - the streaming-texture heat map with a colour map and alpha blending (`rm_draw`, 12064-12097; `gui_speed_heat_rgb` 7374);
    - the fixed-step accumulator loop (21999-22090).
  - A filled polygon or a polyline hull can use `SDL_RenderGeometry`. The installed headers declare it (`/opt/homebrew/include/SDL2/SDL_render.h:1644`); hspist3 never calls it (D2).
  - Headless mode and seeded replay are new code in any option. The lesson of C3: hspist3's window and headless scenes share their step code (`00ALLINONE.c:22090-22899`) but differ in clock handling, guards and failure paths. A new program can call one step function from both paths from the first commit.
  - Cost: small. The drawing needs are disks, lines, filled polygons and one texture.
- **(c) A Python prototype (numpy, shapely, pygame or matplotlib); Demo 1's language decided later.** Feasible for Demo 0 and the quickest way to try growth rules.
  - shapely's buffer and union return polygonal approximations of the swept disks, not exact arcs.
  - Demo 1 needs an exact boundary, as segments and arcs, for exact disk–boundary collision times. With (c) the hull would therefore be built twice: once in Python for Demo 0 and again in C for Demo 1. A conversion would also come in between: a geometry discretisation the blueprint names as failure test (2), §5 line 139: "geometry discretization changes material budget or measured flow".
  - It also adds a second toolchain to version-pin (numpy and shapely versions, on top of the scipy-version sensitivity already recorded in this project's notes).
- **(d) A fork of `edmd_gen3.c` with a segment primitive.** Premature for Demo 0: there are no particles, so 1,548 + 208 lines would be copied before any caller exists.

**Recommendation for Demo 0: (b).** Write the hull as exact geometry from the start: a list of line segments and circular arcs with their centres, radii and angle ranges. The Demo 1 engine then reads exactly what Demo 0 draws, and the headless state file is that list plus the structural disks and the graph. If iterating on growth rules turns out to dominate, (c) is the cheaper fallback for that one purpose, provided its hull exports the same exact segment/arc list.

### Demo 1: frozen hull and small hard disks

- **(a) Link gen3 unmodified and reuse the renderer.** Fails the Demo 1 requirement. The engine accepts only the box with vertical bodies: there is no boundary field in `EDMD_Params` (`edmd.h:15-94`) and no callback in the API (`edmd_gen3.h:126-203`), and the walls are a hard-coded side index (B2). It cannot hold a polyline or arc boundary without an edit, and an edit in `hspist3/` is forbidden. What (a) does give is the **box reference**. G shows the unmodified objects link into a PID test with only `-lm`, so a PID engine can be compared with them on the box case: event hash, event counts, final state and ledger, at the same flags and platform (F2; byte identity needs `-ffp-contract=off` on both sides, `Makefile:34-39`).
- **(b) A new small stepped collision loop for polyline boundaries.** Feasible but the weakest physics.
  - A time-stepped hard-disk loop detects contacts after they happen and needs overlap pushes and clamps. hspist3's own stepped path shows this: overlap-based pair handling with push-apart (`00ALLINONE.c:9507-9535`), clamping at the walls (8423-8494), substeps (9209-9215).
  - Demo 1's energy/momentum check would then measure the integrator's errors as much as the geometry's.
  - The box cross-check against gen3 could only be statistical, never hash-level, because a stepped loop and an event-driven one do not execute the same collisions.
- **(c) Python.** Not a candidate for Demo 1's collision core at the event rates seen in G (21,000 to 186,847 physical events per 40 σ-time cell at N = 400). Language deferred, per the option's wording.
- **(d) A fork of `edmd_gen3.c` inside this folder, adding segment and arc event types.** The strongest candidate for Demo 1. From C5 and C6, the fork must:
  1. Append new event types after `T_PAIR = 5` (`edmd_gen3.c:89`). A box-only run of the fork then pushes no new events, so its tie order (186-191) and event hash (1283) should equal those of unmodified gen3. This is the first regression: fork vs the unmodified objects of (a), on the M1/M2 harness cells, equal hashes and final states. That equality is deduced from the code (C5 point 2), not yet tested.
  2. Give each cell a static list of the segments it can reach, consulted where `wall_candidate` is used today: `predict_all` (498), `exec_cross` (786-787), load (1248).
  3. Predict with absolute positions (`x = cx·w + xi`, as `predict_obj` does at 482) and reflect about a general normal, changing vx and vy together (`exec_wall` flips one component, 742-745).
  4. Add segment terms to the ledger. Momentum is per axis (`ledJ[2]`) and per-wall records are fixed at 4 (`J_wall[4]`, `edmd_gen3.h:104`). For a closed static boundary, energy is conserved exactly by elastic reflection and momentum is not: the wall impulse is booked, as for the outer walls (742, 744).
  5. Add checks and a repair rule to the validator (`local_check`, `full_check`, 567-671; `reflect_into_box` knows only the box, 609-619), and derive its own tolerances. Today's derive from ulp(boxW) (808), and `clean` depends on them (1540-1542).
  6. Add a brute-force segment class to `schedule_audit` (a block next to the four-wall block at 1050-1056) and a branch to `contact_audit` (675). The harness pattern A/B/C (audits on, audits off, again; equal hash and state, `gen3_m1_harness.c:3-6`) is reusable unchanged.
  7. Keep the outer box: the cell grid covers [0, boxW]×[0, boxH] (1136), and `edmd3_load` refuses disks outside it (1238-1244). The hull lies inside the box, with the box walls as an outer safety boundary, unless the fork adds an option to switch them off.
  8. Expect the fork's checkpoints to be incompatible with gen3's, because `sizeof(EDMD3)` is in the checkpoint header (1414, 1444).

  Concave corners at Y junctions, where a disk can touch two segments at the same time, are the hard case. The existing tie-break key and the "tie" stress cell (G) are the template for testing them.

  Fork source: HEAD's `edmd_gen3.c`/`.h`, which are `f72d764`'s plus the display-only hook (25 + 11 inserted lines; F1). The fork's header records the source blob hash and commit. The hook is the one process-global mutable state in gen3 (`edmd_gen3.c:704-706`); a fork can drop it or make it per state.

  One radius suffices for Demo 1: the large structural disks are frozen into the hull, and only the small disks move. A second mobile population would need per-disk radius and mass (B1, C5 point 6).

**Recommendation for Demo 1: (d), with the unmodified gen3 objects of (a) linked into a PID test as the box reference.** (b)'s stepped loop is not recommended for Demo 1. Upstreaming the primitive into `hspist3/edmd_core` is a Paper 1/2 decision for later and is not part of this recommendation.

### The three largest risks

1. **Correctness of the new boundary events (Demo 1).**
   - Grazing contacts at segment joints and at concave corners, and simultaneous contacts with two segments, interact with the tie-break key (`edmd_gen3.c:186-191`), the stagnation guard (`same_limit`, 1157; 1274-1281) and tolerances that are derived for the box only (806-808).
   - The validator's repair is box-only (609-619). A boundary bug would surface as Health counters (`clean=0`), not as a crash, so `clean=1` and the brute-force audit class must be in place before any transport number is read.
2. **Provenance and determinism of PID builds.**
   - The Makefile's `-dirty` test is repo-wide (`Makefile:86`). The main tree carries 22 modified or deleted tracked paths today, two of them gitlinks, so any PID build stamped the same way is permanently `-dirty`. A PID stamp needs a pathspec limited to its own inputs, as `build_d52.sh:10` does for hspist3.
   - Byte-level cross-checks need identical flags and platform (`Makefile:34-39, 71-73`).
   - Binaries written next to a PID executable must not write into the tree. hspist3 writes `00_COMMAND.md` into its cwd and `experiments_*` folders next to `argv[0]` (E1), so any hspist3 run for a cross-check goes in a scratch cwd with the binary called from its `bin_*` folder.
3. **Viewer scope.** The renderer is not reused (D2), so Demo 0's viewer is new code, and the risk is that it grows. hspist3's viewer shows how: key handlers change physics state (`00ALLINONE.c:12237-12256`, 12180-12207), and draw functions write globals and keep state between frames (13042, 13095-13108, 1236-1238). PID's viewer should only read state and own nothing; every result comes from the headless path.

### Blueprint (and brief) assumptions about hspist3 that turned out false

Blueprint, `Physical_Information_Dynamics_Engine_Blueprint_v2_0_2026-10-15.md`:
1. **§0A** (line 25), "an EDMD hard-disk channel/physics implementation": only an axis-aligned rectangle with vertical slabs and faces exists (B2). The blueprint itself marks the statement as user-reported and not inspected.
2. **§0A** (line 29), "sharing an extracted/linked SDL2 renderer ... while leaving `hspist3` unchanged": linking is not possible without editing `00ALLINONE.c` (D2). Only copying, or an untested unity build, remains.
3. **§0A** (line 64), "graph/hull overlays in reused SDL2": the frame has no overlay path, and there are no arc, point-list or scene-aligned field primitives (D1, D2).
4. **§6** (line 159), "reuse through a thin interface when proven feasible": not feasible without edits. The draw entry points take no state and read globals (D1).
5. **§0A** (line 60), "Existing hard-disk particle–particle event prediction can often be reused": only by copying the source or using the whole gen3 state. `pair_rule` is `static inline` (`edmd_gen3.c:407-419`) and is reached only through the box-sized cell grid (C2, C5).
6. **§5 M2** (line 127), "Reuse `hspist3` as reference and optional adapter after confirming its assumptions": holds for the box case (reference via G). An adapter for the polyline/curved case is impossible without an engine change (C4).
7. **§0A** (line 27), "coarse-graining or binning overlay": binning exists only as global-reading 1D histogram panels, plus a 2D position histogram over the window extents that does not reset (`00ALLINONE.c:10136`). There is no scene-aligned field overlay (D3).
8. **§0** (line 21), "a spring-backed movable terminal wall on the **left**": the spring sits on the leftmost internal divider (`00ALLINONE.c:3604`, 21642), not on the outer left wall, which is fixed (`edmd_gen3.c:441`) (B2). It acts as a terminal wall because the compartment behind it is empty.
9. **Unsure.** §0 (line 21), "No ... rewritten event scheduler enters this validation branch": gen3, a rewritten scheduler, is wired into `run_energy_transfer_experiment` (`00ALLINONE.c:23598, 23612-23622`). Whether that path is the Paper 2 validation branch was not determined here.

Brief (`00_M0_BRIEF.md`), for the record:
1. **Rule 1** (line 21), "registered against SHA-256 manifests of `hspist3/00ALLINONE.c` and `hspist3/edmd_core/*`": the manifests hold no engine source. The engine is pinned by commit `f72d764` and by binary SHA-256 (F1). Rule 1 still protects the registered state, because a source edit on main produces an unregistered build, and E0/G6 would have to be rerun.
2. **Item C1**, "`--fixed-dt` and `--seed-drift-order=drift-first` suggest stepping": refuted (C1).
3. **Item D1**, "keys 1–5 and I" as histogram keys: they belong to the gen3 return-map panel. The histogram panels are cycled with V (D1).
4. **Item A1** suggests the Makefile produces `bin_<hash>/`. Scripts produce it (`build_clean.sh`, `build_d52.sh`, ...; A1), and they hard-code the locked `HardDisks_gen3` worktree.

### Engine-change requests

**None is needed for the recommended path.** Option (d) puts every engine change into the fork inside this folder, and (a)'s box reference uses the unmodified objects.

The verifiers listed requests that would help a later reuse. They are recorded here as requests only, for Chris and the HardDisks work to judge. Each changes a build input or a registered binary's bytes:
- a `-dirty` test limited to the build inputs (`Makefile:86`);
- headers as Makefile prerequisites (119);
- a target that builds `edmd_gen3.c`/`edmd.c` objects into a caller-chosen folder;
- an output-root flag and a font-path flag for `00ALLINONE` (E1; `00ALLINONE.c:23500, 23520, 3395-3421`);
- seed, build and engine in `run_params.json` and `00_COMMAND.md` (2282-2317, 3212-3220);
- a renderer split into its own unit with a state struct (D2);
- a per-state collision hook in gen3 (`edmd_gen3.c:704-706`);
- a library-safe error instead of `exit(3)` (197);
- fixes to the dead flags `--seeding`, `--timescale`, `--szilard-unstick-*` and to the duplicate branches at `00ALLINONE.c:4917-4930` (A3).

## Writes of this session and the stop condition

**What this session wrote:**
- Commit `000f190e`, "pid: sealed folder, M0 brief and blueprint v2.0": the two files, staged by explicit path, with `git diff --cached --name-only` empty beforehand.
- `PHYSICALINFORMATIONDYNAMICS/_m0_scratch/`, deleted by the main session with `rm -r` before the status below. The brief (§2 rule 2, §4 G) and Chris's instruction of 2026-10-15 both call for this deletion. It held `g/` (item G) and two folders the binary runs wrote into: `flags_run/` (empty) and `view_run/` (`00_COMMAND.md`, `energy_log.csv` 637,909 B, `run_params.json`, `wall_position.csv`, `wall_position_run_params.json`). Listing at deletion (size, date, path):

```
160 Oct 15 13:16 _m0_scratch
64 Oct 15 13:11 _m0_scratch/flags_run
544 Oct 15 13:07 _m0_scratch/g
43976 Oct 15 13:06 _m0_scratch/g/edmd.o
0 Oct 15 13:06 _m0_scratch/g/edmd.warnings.txt
27160 Oct 15 13:06 _m0_scratch/g/edmd_accelerated.o
0 Oct 15 13:06 _m0_scratch/g/edmd_accelerated.warnings.txt
57312 Oct 15 13:06 _m0_scratch/g/edmd_gen3.o
0 Oct 15 13:06 _m0_scratch/g/edmd_gen3.warnings.txt
141584 Oct 15 13:06 _m0_scratch/g/gen3_m1_harness
158416 Oct 15 13:06 _m0_scratch/g/gen3_m2_harness
8192 Oct 15 13:07 _m0_scratch/g/m1_audit_quick.txt
0 Oct 15 13:06 _m0_scratch/g/m1_harness_build.txt
0 Oct 15 13:06 _m0_scratch/g/m2_harness_build.txt
223 Oct 15 13:06 _m0_scratch/g/nomain_link.txt
143096 Oct 15 13:06 _m0_scratch/g/stub_link
0 Oct 15 13:06 _m0_scratch/g/stub_link.txt
26 Oct 15 13:06 _m0_scratch/g/stub_main.c
224 Oct 15 13:17 _m0_scratch/view_run
389 Oct 15 13:17 _m0_scratch/view_run/00_COMMAND.md
637909 Oct 15 13:17 _m0_scratch/view_run/energy_log.csv
421 Oct 15 13:17 _m0_scratch/view_run/run_params.json
128 Oct 15 13:17 _m0_scratch/view_run/wall_position.csv
421 Oct 15 13:17 _m0_scratch/view_run/wall_position_run_params.json
```

- This file.

The audit commit and a third commit (the handoff bundle, by explicit path, without `.DS_Store`) follow this file. Their hashes are printed by the session and are not knowable inside the file.

**Nothing outside `PHYSICALINFORMATIONDYNAMICS/`.** Evidence:
1. Every compile used `-o` into `_m0_scratch/g/`. No `.o` or harness binary anywhere in `hspist3/` is newer than the brief (`find hspist3 ... -newer PHYSICALINFORMATIONDYNAMICS/00_M0_BRIEF.md`).
2. The agents' transcripts show binary runs only with cwd `_m0_scratch/<agent>_run/`, or with cwd `/` for `--version`; no `make`; and Python only as `python3 -I -` reading files.
3. `hspist3/experiments_gen3_d52_261014/bin_a71f394/` still holds its four files with their old dates. The scene runs called that binary by absolute path, and some modes create `experiments_*` folders next to `argv[0]` (E1); the scene path did not.

**Relaxed stop condition (Chris, 2026-10-15).** Two sessions write to one tree, so the test is "no new entry outside `PHYSICALINFORMATIONDYNAMICS/` that this session wrote". The entries outside the folder that are new since the start of M0 come from the parallel HardDisks work (decision 56):
- `hspist3/cluster/gen3_stage1_261012/s1_figures_261015c.py` (13:09), now committed in `1423d1e5`.
- `hspist3/experiments_gen3_s1_261012/full_261015/figures_261015c/` (13:09), committed in `1423d1e5`.
- `hspist3/experiments_gen3_d56_261015/` (from 13:04), committed in `1423d1e5`; no untracked entry of it remains in the final status below.
- `hspist3/cluster/gen3_stage1_261012/__pycache__/gen_stage1_koa_261012.cpython-313.pyc` (13:13, git-ignored). `import gen_stage1_koa_261012` appears in that session's committed `experiments_gen3_d56_261015/scripts/d56_pilot.py:22` and `d56_d_d2.py:30`.
- `hspist3/experiments_gen3_d50_261014/bin_0544b7e8-dirty/experiments_energy_transfer/` (13:22:29, empty). This is the energy-transfer folder created next to `argv[0]` (`00ALLINONE.c:3419, 18748`). That session's `d56_b_lattice.py` runs the held path with `--experiment=energy_transfer` (its line 5) and names `bin_0544b7e8-dirty`. No agent of this session executed that binary: its last binary run ended at 13:17, and the transcripts show only `ls`/`sed` on that folder.
- `hspist3/kissfft/.git` has a newer mtime (13:47). Unsure what rewrote it; git commands that inspect the nested clone (`git status` in the superproject, `git -C hspist3/kissfft diff`) can refresh its index. Its status is unchanged: ` M kiss_fft_log.h`, `?? .DS_Store`.
- Outside `hspist3/`: `neuraipil/` and the four HardDisks commits listed at the top.

The attribution rests on the git record, timestamps and transcripts. A direct answer from the HardDisks session was not received at writing.

## `git status --short` at the end

Taken after `_m0_scratch/` was deleted and before the audit commit. It is 6939 lines.

Tracked entries, verbatim (22 lines; identical to the 22 at the start):

```
 M .vscode/settings.json
 M 0000_PLAN_OVERALL/paper2_energytransfer/writeup/paper2_draft.pdf
 M 0000_PLAN_OVERALL/paper2_energytransfer/writeup/paper2_draft.tex
 m EDMD_4VALIDATE
 M LAMMPS/log.lammps
 D "Rom\303\241n et al. - 2002 - The speed of sound in a hard disk gas A computer simulation.pdf"
 D "Smallenburg - 2022 - Efficient event-driven simulations of hard spheres.pdf"
 M "hspist3/FINAL speed_of_sound_on_packing_fracture.pdf"
 M hspist3/Power_Freq_Spectrum_with_boxlength_None_and_wallfactor_None.pdf
 M hspist3/Power_Freq_Spectrum_with_boxlength_None_and_wallfactor_None.png
 M hspist3/divider_x_displacement_with_boxlength_None_and_wallfactor_None.pdf
 M hspist3/divider_x_displacement_with_boxlength_None_and_wallfactor_None.png
 M hspist3/energy_log.csv
 M hspist3/experiments_energy_transfer/00_COMMAND.md
 M hspist3/experiments_energy_transfer/energy_transfer_trace.csv
 m hspist3/kissfft
 M hspist3/wall_x_FFT.py
 D log.lammps
 D traj_exact.lammpstrj
 D traj_wall_line.lammpstrj
 D wall_exact.dat
 D wall_line.dat
```

Untracked entries:
- 6917 now, against 6916 at the start.
- 2 are new: 2 inside `PHYSICALINFORMATIONDYNAMICS/` and 0 outside.
- Gone since the start: 1 inside the folder, the start's `?? PHYSICALINFORMATIONDYNAMICS/` (after commit `000f190e` git lists the folder's untracked content per entry instead), and 0 outside it.

New entries inside `PHYSICALINFORMATIONDYNAMICS/`, verbatim:

```
?? PHYSICALINFORMATIONDYNAMICS/M0_AUDIT.md
?? PHYSICALINFORMATIONDYNAMICS/Physical_Information_Dynamics_Handoff_v2_0_2026-10-15/
```

`Physical_Information_Dynamics_Handoff_v2_0_2026-10-15/` was added by Chris at 13:25. Its blueprint copy is byte-identical to the top-level one (`cmp` exit 0; both SHA-256 `ca05004ce5d627a59552ae67825e809d8df2c7cfe1b249baf3660dca0564164d`). It goes into the third commit, without its `.DS_Store`.

New entries outside the folder, per folder; none of them written by this session (see above):

```
(none)
```

Entries gone since the start, per folder:

```
    1 PHYSICALINFORMATIONDYNAMICS/
```
