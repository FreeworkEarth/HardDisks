# -----------------------------
# Project
# -----------------------------
TG              ?= 00ALLINONE
KISSFFT_SRCS    := kissfft/kiss_fft.c kissfft/kiss_fftr.c

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
CFLAGS_RELEASE  := -O3 -ffast-math -march=native $(INCLUDE)

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
.PHONY: all debug release clean tar

# Default target
all: debug

debug: CFLAGS=$(CFLAGS_DEBUG)
debug: $(TG)

release: CFLAGS=$(CFLAGS_RELEASE)
release: $(TG)

# Link rule (libs at the end)
$(TG): $(TG).c $(KISSFFT_SRCS)
	$(CC) $(LDFLAGS) -o $@ $(TG).c $(KISSFFT_SRCS) $(CFLAGS) $(LIBS)

# Clean up the build
clean:
	rm -f $(TG)

# Create a tarball for the project (don’t include the tarball itself)
tar:
	tar czvf $(TG).tgz --exclude='$(TG).tgz' *