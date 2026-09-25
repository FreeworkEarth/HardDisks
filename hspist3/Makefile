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
EDMD_SRCS       ?= edmd_core/edmd.c edmd_core/edmd_accelerated.c
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
CFLAGS_RELEASE  := -O3 -march=native $(INCLUDE)

# Linker search paths (keep your Homebrew/Mesa ones)
LDFLAGS        += -L/opt/homebrew/lib -L/opt/homebrew/opt/mesa/lib

# OpenGL / Core frameworks per platform
ifeq ($(PLATFORM),Darwin)
  GL_LIBS      := -framework OpenGL -framework CoreGraphics
else
  GL_LIBS      := -lGL
endif

# Libraries: prefer pkg-config libs if present; add GLEW and GL last
LIBS_BASE       := -lGLEW $(GL_LIBS)
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
CFLAGS_KOA      := -O2 -march=x86-64-v2 -mtune=generic $(INCLUDE)

debug: CFLAGS=$(CFLAGS_DEBUG)
debug: $(TG)

release: CFLAGS=$(CFLAGS_RELEASE)
release: $(TG)

koa: CFLAGS=$(CFLAGS_KOA)
koa: $(TG)

# Link rule (libs at the end)
$(TG): $(TG).c $(KISSFFT_SRCS) $(EDMD_SRCS) $(VALIDATION_SRCS)
	$(CC) $(LDFLAGS) -o $@ $(TG).c $(KISSFFT_SRCS) $(EDMD_SRCS) $(VALIDATION_SRCS) $(CFLAGS) $(LIBS)

# Clean up the build
clean:
	rm -f $(TG)

# Create a tarball for the project (don’t include the tarball itself)
tar:
	tar czvf $(TG).tgz --exclude='$(TG).tgz' *
