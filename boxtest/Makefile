# fill TG in with program you want to compile
# for lldb use: make DBG=-g
# for webassembly: make em

TG = boxtest

CC      ?= cc
DBG     ?= -O3
LTO     ?= -flto

# Try pkg-config first
PKG_SDL_CFLAGS := $(shell pkg-config --cflags sdl2 SDL2_ttf 2>/dev/null)
PKG_SDL_LIBS   := $(shell pkg-config --libs   sdl2 SDL2_ttf 2>/dev/null)

# Fallback (Homebrew default prefix)
INCLUDE_FALLBACK := -I/opt/homebrew/include -I/opt/homebrew/include/SDL2 -D_THREAD_SAFE
LIBS_FALLBACK    := -L/opt/homebrew/lib -lSDL2 -lSDL2_ttf

CFLAGS  = $(DBG) $(LTO) $(if $(PKG_SDL_CFLAGS),$(PKG_SDL_CFLAGS),$(INCLUDE_FALLBACK))
LIBS    = $(if $(PKG_SDL_LIBS),$(PKG_SDL_LIBS),$(LIBS_FALLBACK)) -lm
LDFLAGS = $(LTO)

OBJ = $(TG).o macsdl2.o SDL_FontCache.o

.PHONY: all clean tar em gzip

all: $(TG)

$(TG): $(OBJ)
	$(CC) $(LDFLAGS) -o $@ $(OBJ) $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *.o *.bc *.js *.wasm *.data $(TG)

tar:
	tar czvf $(TG).tgz *

gzip:
	gzip -k index.wasm

# -------- Emscripten (optional) --------
EM_CFLAGS  = -c -s WASM=1 -O3 -flto
EM_LDFLAGS = -s USE_SDL=2 -s USE_SDL_TTF=2 --no-heap-copy -s TOTAL_MEMORY=268435456 -s ALLOW_MEMORY_GROWTH=1

em:
	emcc $(TG).c $(EM_CFLAGS) -o $(TG).bc
	emcc macsdl2.c $(EM_CFLAGS) -o macsdl2.bc
	emcc SDL_FontCache.c $(EM_CFLAGS) -o SDL_FontCache.bc
	emcc $(TG).bc macsdl2.bc SDL_FontCache.bc $(EM_LDFLAGS) -o index.js --preload-file Crimson-Roman.ttf
