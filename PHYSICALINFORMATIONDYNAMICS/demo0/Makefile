# PID Demo 0 -- the one Makefile (brief 01_DEMO0_BRIEF.md sec. 2). Run make in this folder.
#   make          bin/pid_demo0 (window and headless) and bin/pid_check (test tool)
#   make test     zero-warning check, then tests/run_tests.sh (headless; the GUI = headless test uses SDL's dummy driver)
CC       ?= cc
# -ffp-contract=off from the first compile: no fused multiply-add, so builds do not re-round the geometry by inlining
# (the reason is hspist3/Makefile:34-40).
CSTD     := -std=c11 -O2 -Wall -Wextra -ffp-contract=off
SDL_CF   := $(shell pkg-config --cflags sdl2)
SDL_LB   := $(shell pkg-config --libs sdl2)
# SDL2_ttf is optional: used only for --font=<path>; without it the status line goes to the terminal.
HAVE_TTF := $(shell pkg-config --exists SDL2_ttf && echo 1)
ifeq ($(HAVE_TTF),1)
  TTF_CF := $(shell pkg-config --cflags SDL2_ttf) -DPID_HAVE_TTF
  SDL_LB := $(shell pkg-config --libs sdl2 SDL2_ttf)
endif

# Provenance. The dirty test covers PHYSICALINFORMATIONDYNAMICS/ only (not the whole repository, where the HardDisks work
# always has uncommitted changes). SRC_SHA256 is the sha256 of the files below, concatenated in this order; recompute it
# from any commit with: cat $(SRC_LIST) | shasum -a 256
MK_DIR    := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
PID_DIR   := $(abspath $(MK_DIR)..)
SRC_LIST  := Makefile pid.h pid_rng.c pid_sha256.c pid_state.c pid_grow.c pid_main.c pid_view.c pid_check.c pid_hull.c
GIT_SHORT := $(shell git -C $(PID_DIR) log -1 --format=%h -- . 2>/dev/null)
GIT_FULL  := $(shell git -C $(PID_DIR) log -1 --format=%H -- . 2>/dev/null)
GIT_DIRTY := $(shell git -C $(PID_DIR) status --porcelain -- . 2>/dev/null | head -1)
BUILD_ID  := $(if $(GIT_SHORT),$(GIT_SHORT),unknown)$(if $(GIT_DIRTY),-dirty,)
SRC_SHA   := $(shell cd $(MK_DIR) && cat $(SRC_LIST) | shasum -a 256 | cut -c1-64)
DEFS      := -DPID_BUILD_ID='"$(BUILD_ID)"' -DPID_CFLAGS='"$(CSTD)"' -DPID_FOLDER_COMMIT='"$(GIT_FULL)"' -DPID_SRC_SHA256='"$(SRC_SHA)"'

# Rebuild when the build id or the sources' hash changes (a commit changes the id without touching a source).
STAMP := bin/.build_id
$(shell mkdir -p $(MK_DIR)bin; echo '$(BUILD_ID) $(SRC_SHA)' | cmp -s - $(MK_DIR)$(STAMP) || echo '$(BUILD_ID) $(SRC_SHA)' > $(MK_DIR)$(STAMP))

CORE := pid_grow.c pid_state.c pid_rng.c pid_sha256.c pid_hull.c
HDRS := pid.h

.PHONY: all test clean
all: bin/pid_demo0 bin/pid_check

bin/pid_demo0: pid_main.c pid_view.c $(CORE) $(HDRS) Makefile $(STAMP)
	$(CC) $(CSTD) $(SDL_CF) $(TTF_CF) $(DEFS) -o $@ pid_main.c pid_view.c $(CORE) $(SDL_LB) -lm

bin/pid_check: pid_check.c $(CORE) $(HDRS) Makefile $(STAMP)
	$(CC) $(CSTD) $(DEFS) -o $@ pid_check.c $(CORE) -lm

# Zero warnings is the bar: every source compiled once more with -Werror (objects to /dev/null), then the tests.
test: all
	@for f in pid_main.c pid_view.c pid_check.c $(CORE); do \
	  $(CC) $(CSTD) -Werror $(SDL_CF) $(TTF_CF) $(DEFS) -c -o /dev/null $$f || exit 1; done; \
	  echo "PASS build: 8 sources, 0 warnings under $(CSTD) -Werror"
	@sh tests/run_tests.sh

clean:
	rm -f bin/pid_demo0 bin/pid_check
