# -----------------------------
# Project
# -----------------------------
TG              ?= 00ALLINONE
KISSFFT_SRCS    := kissfft/kiss_fft.c kissfft/kiss_fftr.c
# Optional EDMD core (event-driven hard disks)
# Default: build both EDMD backends (default + accelerated) so runtime can choose.
# You can still override, e.g.:
#   make EDMD_SRCS=edmd_core/edmd.c
#   make EDMD_SRCS=edmd_core/edmd_accelerated.c
# ##CHRIS 2026-10-09 (M3, 261012 sec. 4.7.14): the generation-3 engine is linked into every build; --engine=gen3 selects it
EDMD_SRCS       ?= edmd_core/edmd.c edmd_core/edmd_accelerated.c edmd_core/edmd_gen3.c
VALIDATION_SRCS := experiment_validation.c

# Compiler
CC              ?= cc

# Detect platform (Darwin = macOS, else assume Linux)
PLATFORM        := $(shell uname)

# Try to pull SDL2 flags via pkg-config (falls back to Homebrew paths if not found)
PKG_SDL_CFLAGS  := $(shell pkg-config --cflags sdl2 SDL2_ttf 2>/dev/null)
PKG_SDL_LIBS    := $(shell pkg-config --libs   sdl2 SDL2_ttf 2>/dev/null)

# Include paths (Homebrew defaults kept as fallback)
INCLUDE_FALLBACK := -I/opt/homebrew/include -I/opt/homebrew/include/SDL2 -D_THREAD_SAFE
INCLUDE         := $(if $(PKG_SDL_CFLAGS),$(PKG_SDL_CFLAGS),$(INCLUDE_FALLBACK))

# Build types
CFLAGS_DEBUG    := -g -fsanitize=address -DDEBUG $(INCLUDE)
# Scientific runs depend on IEEE NaN/Inf checks and ordered collision times.
# -ffast-math is intentionally excluded because it lets the compiler assume
# non-finite values never occur and can invalidate experiment safety guards.
# ##CHRIS 2026-10-10: -ffp-contract=off. Without it clang fuses a*b+c into one FMA wherever it
# likes, and WHERE it likes depends on inlining -- so adding code anywhere in the program silently
# re-rounds the collision arithmetic by one ulp and two builds of nominally identical physics stop
# being byte-identical. That is what made gate (a) red for a day. Byte-identity comparisons across
# builds are meaningless without this flag. Costs a few per cent; buys a gate that means something,
# and matters more on KOA where nodes may fuse differently.
CFLAGS_RELEASE  := -O3 -march=native -ffp-contract=off $(INCLUDE)

# Linker search paths (keep your Homebrew/Mesa ones)
LDFLAGS        += -L/opt/homebrew/lib -L/opt/homebrew/opt/mesa/lib

# OpenGL / Core frameworks per platform
ifeq ($(PLATFORM),Darwin)
  GL_LIBS      := -framework OpenGL -framework CoreGraphics
else
  GL_LIBS      := -lGL
  # ##CHRIS 2026-10-02: on Linux the maths library is separate and must be named; the first KOA link failed with
  # "undefined reference to symbol 'acos@@GLIBC_2.2.5' ... libm.so.6: DSO missing from command line". On macOS
  # libm is part of libSystem, so SYS_LIBS stays empty there and the Mac link line is unchanged.
  SYS_LIBS     := -lm
endif

# Libraries: prefer pkg-config libs if present; add GLEW and GL last
LIBS_BASE       := -lGLEW $(GL_LIBS) $(SYS_LIBS)
LIBS            := $(if $(PKG_SDL_LIBS),$(PKG_SDL_LIBS),-lSDL2 -lSDL2_ttf) $(LIBS_BASE)

# -----------------------------
# Targets
# -----------------------------
.PHONY: all debug release koa clean tar

# Default target
all: debug

# ##CHRIS 2026-09-16: cluster build (KOA, UH Manoa). Same sources and same physics as `release`;
# the only difference is the instruction set. -march=native would bake in whatever login node built
# it and could then differ from the compute nodes, so a fixed, portable ISA is used instead: every
# KOA node produces the same bytes. Cluster traces are NEVER byte-identical to the macOS/ARM ones
# (different architecture, different FMA contraction) and must not be mixed with them inside one
# cell or figure; they are gated statistically -- see cluster/README.md.
#
# ##CHRIS 2026-10-01: x86-64-v3 CORRECTED TO v2. KOA's shared partition is heterogeneous and
# includes Ivy Bridge nodes; v3 requires AVX2/FMA (Haswell and later), so a v3 binary dies with
# SIGILL the moment Slurm lands the job on an older node -- intermittently, and only on some array
# tasks, which is the worst way to find out. v2 (SSE4.2/POPCNT, Nehalem and later) runs on every
# node in the partition. Do NOT raise this without pinning --constraint to a node feature.
CFLAGS_KOA      := -O2 -march=x86-64-v2 -mtune=generic -ffp-contract=off $(INCLUDE)

# ##CHRIS 2026-10-10: build provenance, injected at compile time. The binary is no longer tracked
# in git (it is a build artefact that every `make release` dirties), so the record of WHICH binary
# produced a result travels with the result instead: --version prints these, and every
# energy-transfer summary CSV carries them as build_git / build_target / build_cflags.
BUILD_GIT       := $(shell git -C $(dir $(abspath $(lastword $(MAKEFILE_LIST)))) rev-parse --short HEAD 2>/dev/null || echo unknown)$(shell git -C $(dir $(abspath $(lastword $(MAKEFILE_LIST)))) diff --quiet HEAD 2>/dev/null || echo -dirty)
PROVENANCE       = -DBUILD_GIT='"$(BUILD_GIT)"' -DBUILD_TARGET='"$(BUILD_TARGET)"' -DBUILD_CFLAGS='"$(strip $(BUILD_CF))"'

debug: BUILD_TARGET=debug
debug: BUILD_CF=-g -fsanitize=address -DDEBUG
debug: CFLAGS=$(CFLAGS_DEBUG)
debug: $(TG)

release: BUILD_TARGET=release
release: BUILD_CF=-O3 -march=native -ffp-contract=off
release: CFLAGS=$(CFLAGS_RELEASE)
release: $(TG)

koa: BUILD_TARGET=koa
koa: BUILD_CF=-O2 -march=x86-64-v2 -mtune=generic -ffp-contract=off
koa: CFLAGS=$(CFLAGS_KOA)
koa: $(TG)

# Link rule (libs at the end)
$(TG): $(TG).c $(KISSFFT_SRCS) $(EDMD_SRCS) $(VALIDATION_SRCS)
	$(CC) $(LDFLAGS) -o $@ $(TG).c $(KISSFFT_SRCS) $(EDMD_SRCS) $(VALIDATION_SRCS) $(CFLAGS) $(PROVENANCE) $(LIBS)

# Clean up the build
clean:
	rm -f $(TG)

# Create a tarball for the project (don’t include the tarball itself)
tar:
	tar czvf $(TG).tgz --exclude='$(TG).tgz' *
