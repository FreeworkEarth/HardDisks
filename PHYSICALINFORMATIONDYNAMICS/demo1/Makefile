# PID Demo 1 -- the one Makefile (brief 02_DEMO1_BRIEF.md sec. 2, 6 and 8; pre-check 02_DEMO1_BRIEF_PRECHECK.md parts A, B).
# Run make in this folder. It writes only into bin/ and runs/ here (both git-ignored); hspist3/edmd_core/ and ../demo0/ are
# only read (no make there).
#   make                 the box-reference binaries for both sides (fork and unmodified original), the final-state driver,
#                        the wall tests and controls, bin/pid_demo1 (the driver, sec. 5), Demo 0 built headless, and where
#                        SDL2 is installed bin/pid_view1 (the window, sec. 9), with SDL2_ttf too bin/pid_studio (grow the network
#                        and run the disks from one window with panels)
#   make test            deliveries 1a-1c and the studio, tests 2-12 (tests/run_tests.sh): a full rebuild, then about 9 min on 4 cores;
#                        outputs in runs/box_reference_<cc>_<UTC time>/. One make test at a time per checkout: bin/ is shared.
#   make shots           1c: the two PNGs of sec. 9 into shots/ (tests/make_shots.sh; SDL's dummy driver, no display needed)
#   make test CC=clang   the same with another compiler (one run always compares the two sides under ONE compiler and flag set)
CC       ?= cc
CSTD     := -std=c11 -O2 -Wall -Wextra -ffp-contract=off -Werror
HSP      := ../../hspist3/edmd_core
# 1 under clang (Apple clang included), 0 under gcc: only clang's preprocessor replaces __clang__ by 1.
IS_CLANG := $(shell echo __clang__ | $(CC) -E -P - 2>/dev/null | grep -c '^1$$')

# Flags for code that 1a may not edit, the same on both sides (pre-check A1, A2):
#   gen2, hspist3/edmd_core/edmd.c (read-only; linked unchanged for the harnesses' run D): under -std=c11 glibc declares M_PI
#   only with a feature macro (M_PI is used on 9 lines from edmd.c:556; gcc reports 8), hence -D_DEFAULT_SOURCE (macOS ignores
#   it); gcc also warns at edmd.c:719, 765, 1591, 1592 (-Wmisleading-indentation).
#   gen3, fork and original alike: gcc warns at line 842 of the original (fork line 846: snprintf "D%d" into char kb[8],
#   -Wformat-truncation); the line cannot change in 1a (test 2), so gcc's warning is off for that file on both sides.
GEN2_X   := -D_DEFAULT_SOURCE $(if $(filter 1,$(IS_CLANG)),,-Wno-misleading-indentation)
GEN3_X   := $(if $(filter 1,$(IS_CLANG)),,-Wno-format-truncation)

# The include trap (pre-check A3): the harnesses include "../edmd_gen3.h", "../edmd.h" and "gen3_state_io.h", which the
# compiler looks up next to the harness file before any -I, so a harness compiled in place always sees hspist3's headers.
# The fork side therefore compiles unchanged COPIES of the harness files in bin/hfork/tests/, whose parent bin/hfork/ holds
# copies of the fork's headers; each fork-side link stops if a copy differs from its source or if its dependency list names
# a header of the other side.
MIR      := bin/hfork
HARN     := gen3_m1_harness gen3_m2_harness
FORK_H   := engine/edmd_gen3.h engine/edmd.h
ORIG_H   := $(HSP)/edmd_gen3.h $(HSP)/edmd.h
BINS     := $(addprefix bin/orig/,$(HARN) box_final_state) $(addprefix bin/fork/,$(HARN) box_final_state)
# 1b, fork only: the wall tests on hand-built geometry (tests 4-6) and the positive controls (test 10, white-box)
BINS     += bin/fork/prim_tests bin/fork/prim_controls
# 1b: the driver and Demo 0 built headless (test 7)
BINS     += bin/pid_demo1 bin/demo0/pid_demo0_headless
# 1c: the window (brief sec. 9), only where SDL2 is installed (pkg-config sdl2); everything else builds without it
SDL_OK   := $(shell pkg-config --exists sdl2 2>/dev/null && echo 1)
ifeq ($(SDL_OK),1)
SDL_CF   := $(shell pkg-config --cflags sdl2)
SDL_LB   := $(shell pkg-config --libs sdl2)
BINS     += bin/pid_view1
endif
# the studio (pid_studio.c: grow the network and run the disks from one window with panels) also needs SDL2_ttf
TTF_OK   := $(shell pkg-config --exists sdl2 SDL2_ttf 2>/dev/null && echo 1)
ifeq ($(TTF_OK),1)
TTF_CF   := $(shell pkg-config --cflags sdl2 SDL2_ttf)
TTF_LB   := $(shell pkg-config --libs sdl2 SDL2_ttf)
BINS     += bin/pid_studio
endif
# Rebuild when the compiler, the flags or this Makefile change (make test CC=clang after make test): every object and
# program depends on the stamp (rewritten only when its content changes) and on the Makefile. Quotes are dropped from the
# compiler's version line so that the shell line below cannot break.
STAMP    := bin/.compiler
CC_LINE  := $(subst ',,$(CC) | $(shell $(CC) --version 2>/dev/null | head -1) | $(CSTD) | $(GEN2_X) | $(GEN3_X))
$(shell mkdir -p bin; echo '$(CC_LINE)' | cmp -s - $(STAMP) || echo '$(CC_LINE)' > $(STAMP))
DEPS     := $(STAMP) Makefile

.PHONY: all test shots
all: $(BINS)

bin/gen2.o: $(HSP)/edmd.c $(HSP)/edmd.h $(DEPS)
	@mkdir -p $(@D)
	$(CC) $(CSTD) $(GEN2_X) -c -o $@ $<
bin/orig/edmd_gen3.o: $(HSP)/edmd_gen3.c $(ORIG_H) $(DEPS)
	@mkdir -p $(@D)
	$(CC) $(CSTD) $(GEN3_X) -c -o $@ $<
bin/fork/edmd_gen3.o: engine/edmd_gen3.c $(FORK_H) $(DEPS)
	@mkdir -p $(@D)
	$(CC) $(CSTD) $(GEN3_X) -c -o $@ $<

$(MIR)/edmd_gen3.h: engine/edmd_gen3.h
	@mkdir -p $(@D) && cp $< $@
$(MIR)/edmd.h: engine/edmd.h
	@mkdir -p $(@D) && cp $< $@
$(MIR)/tests/%: $(HSP)/tests/%
	@mkdir -p $(@D) && cp $< $@

$(addprefix bin/orig/,$(HARN)): bin/orig/%: $(HSP)/tests/%.c $(HSP)/tests/gen3_state_io.h $(ORIG_H) bin/orig/edmd_gen3.o bin/gen2.o $(DEPS)
	$(CC) $(CSTD) -MM $< > $@.deps
	@if grep -q -e engine/ -e hfork $@.deps || ! grep -q 'edmd_core/tests/\.\./edmd_gen3\.h' $@.deps; then \
	  echo "include check: $@ is not built against hspist3's headers (see $@.deps)"; exit 1; fi
	$(CC) $(CSTD) -o $@ $< bin/orig/edmd_gen3.o bin/gen2.o -lm
$(addprefix bin/fork/,$(HARN)): bin/fork/%: $(MIR)/tests/%.c $(MIR)/tests/gen3_state_io.h $(MIR)/edmd_gen3.h $(MIR)/edmd.h bin/fork/edmd_gen3.o bin/gen2.o $(DEPS)
	@for f in edmd_gen3.h edmd.h; do cmp -s engine/$$f $(MIR)/$$f || { echo "stale copy: $(MIR)/$$f differs from engine/$$f"; exit 1; }; done; \
	 for f in $*.c gen3_state_io.h; do cmp -s $(HSP)/tests/$$f $(MIR)/tests/$$f || { echo "stale copy: $(MIR)/tests/$$f differs from its source"; exit 1; }; done
	$(CC) $(CSTD) -MM $< > $@.deps
	@if grep -q hspist3 $@.deps || ! grep -q 'hfork/tests/\.\./edmd_gen3\.h' $@.deps; then \
	  echo "include check: $@ is not built against the fork's headers (see $@.deps)"; exit 1; fi
	$(CC) $(CSTD) -o $@ $< bin/fork/edmd_gen3.o bin/gen2.o -lm

# The final-state driver (tests/box_final_state.c, pre-check A4): one source, built once per side; its "edmd_gen3.h" is found
# through -I only (tests/ holds no such file), and the same dependency check applies.
bin/orig/box_final_state: tests/box_final_state.c $(ORIG_H) bin/orig/edmd_gen3.o $(DEPS)
	@mkdir -p $(@D)
	$(CC) $(CSTD) -I$(HSP) -MM $< > $@.deps
	@if grep -q -e engine/ -e hfork $@.deps || ! grep -q 'edmd_core/edmd_gen3\.h' $@.deps; then \
	  echo "include check: $@ is not built against hspist3's headers (see $@.deps)"; exit 1; fi
	$(CC) $(CSTD) -I$(HSP) -o $@ $< bin/orig/edmd_gen3.o -lm
bin/fork/box_final_state: tests/box_final_state.c $(FORK_H) bin/fork/edmd_gen3.o $(DEPS)
	@mkdir -p $(@D)
	$(CC) $(CSTD) -Iengine -MM $< > $@.deps
	@if grep -q hspist3 $@.deps || ! grep -q 'engine/edmd_gen3\.h' $@.deps; then \
	  echo "include check: $@ is not built against the fork's headers (see $@.deps)"; exit 1; fi
	$(CC) $(CSTD) -Iengine -o $@ $< bin/fork/edmd_gen3.o -lm

# 1b: Demo 0's generator, the one generator for seeding (brief sec. 2 rule 5), compiled unchanged from ../demo0 (read only)
D0       := ../demo0
bin/demo0/pid_rng.o: $(D0)/pid_rng.c $(D0)/pid.h $(DEPS)
	@mkdir -p $(@D)
	$(CC) $(CSTD) -c -o $@ $<

# 1b: tests/prim_tests.c against the fork (the original has no walls); tests/prim_controls.c includes engine/edmd_gen3.c with
# -DPID_WHITEBOX (the planted-error switch exists only in that build), as hspist3's gen3_band_edge_ties.c includes the original
bin/fork/prim_tests: tests/prim_tests.c $(FORK_H) $(D0)/pid.h bin/fork/edmd_gen3.o bin/demo0/pid_rng.o $(DEPS)
	@mkdir -p $(@D)
	$(CC) $(CSTD) -Iengine -I$(D0) -MM $< > $@.deps
	@if grep -q hspist3 $@.deps || ! grep -q 'engine/edmd_gen3\.h' $@.deps; then \
	  echo "include check: $@ is not built against the fork's headers (see $@.deps)"; exit 1; fi
	$(CC) $(CSTD) -Iengine -I$(D0) -o $@ $< bin/fork/edmd_gen3.o bin/demo0/pid_rng.o -lm
bin/fork/prim_controls: tests/prim_controls.c engine/edmd_gen3.c $(FORK_H) $(D0)/pid.h bin/demo0/pid_rng.o $(DEPS)
	@mkdir -p $(@D)
	$(CC) $(CSTD) $(GEN3_X) -I$(D0) -DPID_WHITEBOX -o $@ $< bin/demo0/pid_rng.o -lm

# 1b: the driver (pid_demo1.c, brief sec. 5) and Demo 0 built headless for test 7 (pre-check B12): Demo 0's unchanged sources,
# compiled here from ../demo0 (read only; no make there, no SDL2), with tests/demo0_view_stub.c in place of the window and the
# provenance macros demo0/Makefile passes (its build id, flags, folder commit and source hash, computed the same way).
D0_CORE  := $(addprefix $(D0)/,pid_grow.c pid_state.c pid_rng.c pid_sha256.c pid_hull.c)
D0_LIST  := Makefile pid.h pid_rng.c pid_sha256.c pid_state.c pid_grow.c pid_main.c pid_view.c pid_check.c pid_hull.c
SHA256   := $(if $(shell command -v shasum),shasum -a 256,sha256sum)
D0_SHORT := $(shell git -C .. log -1 --format=%h -- . 2>/dev/null)
D0_FULL  := $(shell git -C .. log -1 --format=%H -- . 2>/dev/null)
D0_DIRTY := $(shell git -C .. status --porcelain -- . 2>/dev/null | head -1)
D0_ID    := $(if $(D0_SHORT),$(D0_SHORT),unknown)$(if $(D0_DIRTY),-dirty,)
D0_SHA   := $(shell cd $(D0) && cat $(D0_LIST) | $(SHA256) | cut -c1-64)
D0_DEFS  := -DPID_BUILD_ID='"$(D0_ID)"' -DPID_CFLAGS='"$(CSTD)"' -DPID_FOLDER_COMMIT='"$(D0_FULL)"' -DPID_SRC_SHA256='"$(D0_SHA)"'
bin/demo0/pid_demo0_headless: $(D0)/pid_main.c $(D0_CORE) $(D0)/pid.h tests/demo0_view_stub.c $(DEPS)
	@mkdir -p $(@D)
	$(CC) $(CSTD) $(D0_DEFS) -I$(D0) -o $@ $(D0)/pid_main.c $(D0_CORE) tests/demo0_view_stub.c -lm
# pid_demo1: the fork object, Demo 0's reader, hull, generator and sha256 (unchanged), and the provenance of the manifest:
# demo1's build id, the flags, the fork's copy commit and source blobs (from the engine headers) and the sources' sha256.
D1_SHORT := $(shell git log -1 --format=%h -- . 2>/dev/null)
D1_DIRTY := $(shell git status --porcelain -- . 2>/dev/null | head -1)
D1_ID    := $(if $(D1_SHORT),$(D1_SHORT),unknown)$(if $(D1_DIRTY),-dirty,)
D1_SHA   := $(shell cat engine/edmd_gen3.c engine/edmd_gen3.h engine/edmd.h pid_demo1.h pid_demo1.c pid_draw1.h pid_draw1.c pid_view1.c pid_studio.c | $(SHA256) | cut -c1-64)
FORK_CP  := $(shell sed -n 's/^   copied at origin\/main \([0-9a-f]*\);.*/\1/p' engine/edmd_gen3.c)
FORK_BL  := $(foreach f,edmd_gen3.c edmd_gen3.h edmd.h,$(f) $(shell sed -n 's/^   source blob \([0-9a-f]*\)\..*/\1/p' engine/$(f)))
D1_DEFS  := -DPID1_BUILD_ID='"$(D1_ID)"' -DPID1_CFLAGS='"$(CC) $(CSTD)"' -DPID1_SRC_SHA256='"$(D1_SHA)"' \
            -DPID1_FORK='"hspist3/edmd_core at origin/main $(FORK_CP); source blobs $(FORK_BL)"'
bin/pid_demo1: pid_demo1.c pid_demo1.h $(FORK_H) $(D0)/pid.h bin/fork/edmd_gen3.o $(D0_CORE) $(DEPS)
	@mkdir -p $(@D)
	$(CC) $(CSTD) -Iengine -I$(D0) -MM $< > $@.deps
	@if grep -q hspist3 $@.deps || ! grep -q 'engine/edmd_gen3\.h' $@.deps; then \
	  echo "include check: $@ is not built against the fork's headers (see $@.deps)"; exit 1; fi
	$(CC) $(CSTD) $(D1_DEFS) -Iengine -I$(D0) -o $@ $< bin/fork/edmd_gen3.o $(D0_CORE) -lm
# 1c: the window, the same driver with -DPID1_VIEW and pid_view1.c (SDL2); the same include check on both sources
bin/pid_view1: pid_demo1.c pid_view1.c pid_draw1.c pid_draw1.h pid_demo1.h $(FORK_H) $(D0)/pid.h bin/fork/edmd_gen3.o $(D0_CORE) $(DEPS)
	@mkdir -p $(@D)
	$(CC) $(CSTD) -DPID1_VIEW $(SDL_CF) -Iengine -I$(D0) -MM pid_demo1.c pid_view1.c pid_draw1.c > $@.deps
	@if grep -q hspist3 $@.deps || [ "$$(grep -c 'engine/edmd_gen3\.h' $@.deps)" -lt 3 ]; then \
	  echo "include check: $@ is not built against the fork's headers (see $@.deps)"; exit 1; fi
	$(CC) $(CSTD) $(D1_DEFS) -DPID1_VIEW $(SDL_CF) -Iengine -I$(D0) -o $@ pid_demo1.c pid_view1.c pid_draw1.c bin/fork/edmd_gen3.o $(D0_CORE) $(SDL_LB) -lm

# the studio: the driver's functions (pid_demo1.c without its main), the drawing, Demo 0's growth, SDL2 and SDL2_ttf
bin/pid_studio: pid_studio.c pid_demo1.c pid_draw1.c pid_draw1.h pid_demo1.h $(FORK_H) $(D0)/pid.h bin/fork/edmd_gen3.o $(D0_CORE) $(DEPS)
	@mkdir -p $(@D)
	$(CC) $(CSTD) -DPID1_NO_MAIN $(TTF_CF) -Iengine -I$(D0) -MM pid_studio.c pid_demo1.c pid_draw1.c > $@.deps
	@if grep -q hspist3 $@.deps || [ "$$(grep -c 'engine/edmd_gen3\.h' $@.deps)" -lt 3 ]; then \
	  echo "include check: $@ is not built against the fork's headers (see $@.deps)"; exit 1; fi
	$(CC) $(CSTD) $(D1_DEFS) -DPID1_NO_MAIN $(TTF_CF) -Iengine -I$(D0) -o $@ pid_studio.c pid_demo1.c pid_draw1.c bin/fork/edmd_gen3.o $(D0_CORE) $(TTF_LB) -lm

# make test rebuilds everything first (a few seconds), so the tested programs are always built from the files as they are;
# every command of that build, as make runs it, goes into the run's manifest (pre-check A10).
test:
	@$(MAKE) -s --no-print-directory -B all CC="$(CC)"
	@$(MAKE) -s --no-print-directory -n -B all CC="$(CC)" > bin/build_lines.txt
	@CC="$(CC)" IS_CLANG="$(IS_CLANG)" CSTD="$(CSTD)" sh tests/run_tests.sh
shots: bin/pid_view1 bin/demo0/pid_demo0_headless $(if $(filter 1,$(TTF_OK)),bin/pid_studio,)
	@CC="$(CC)" IS_CLANG="$(IS_CLANG)" sh tests/make_shots.sh
