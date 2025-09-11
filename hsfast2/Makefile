# fill TG in with program you want to compile
# for lldb use  % make DBG=-g
# for webassembly, % make em

#TG = circle
#TG = flipo
#TG = hsfast
#TG = hsfast1
#TG = gaussplt
TG = hsfast2

CC ?= cc
DBG ?= -O3

# SDL2 include and libs (auto-detect or fallback)
PKG_SDL_CFLAGS  := $(shell pkg-config --cflags sdl2 SDL2_ttf 2>/dev/null)
PKG_SDL_LIBS    := $(shell pkg-config --libs   sdl2 SDL2_ttf 2>/dev/null)
INCLUDE_FALLBACK := -I/opt/homebrew/include -I/opt/homebrew/include/SDL2 -D_THREAD_SAFE
INCLUDE         := $(if $(PKG_SDL_CFLAGS),$(PKG_SDL_CFLAGS),$(INCLUDE_FALLBACK))
# Add -L fallback so linking works without pkg-config
LIBS            := $(if $(PKG_SDL_LIBS),$(PKG_SDL_LIBS),-L/opt/homebrew/lib -lSDL2 -lSDL2_ttf) -lm

CFLAGS = $(DBG) $(INCLUDE)

OBJ = $(TG).o macsdl2.o gaussrand.o SDL_FontCache.o

$(TG): $(OBJ)
	$(CC) -o $(TG) $(OBJ) $(LIBS)

HEADERS = macsdl2.h
$(OBJ): $(HEADERS)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *.o *.bc *.js *.wasm $(TG)

tar:
	tar czvf $(TG).tgz *

# ---------------- Emscripten (optional) ----------------
IDIR = /Users/rob/emscript/emsdk/fastcomp/enscripten/system/include
EM_ENV=LLVM=/Users/rob/emscript/emsdk NODE_JS=node EMSCRIPTEN_ROOT=/Users/rob/emscript/emsdk
EM_CFLAGS = -s WASM=1 -O3 -I$(IDIR)
EM_LDFLAGS = -s USE_SDL=2 -s USE_FREETYPE=1 -s USE_SDL_TTF=2 --no-heap-copy -s TOTAL_MEMORY=100663296

em:
	emcc $(TG).c $(EM_CFLAGS) $(EM_LDFLAGS) -o $(TG).bc
	emcc macsdl2.c $(EM_CFLAGS) $(EM_LDFLAGS) -o macsdl2.bc
	emcc SDL_FontCache.c $(EM_CFLAGS) $(EM_LDFLAGS) -o SDL_FontCache.bc
	emcc gaussrand.c $(EM_CFLAGS) $(EM_LDFLAGS) -o gaussrand.bc
	emcc $(TG).bc macsdl2.bc SDL_FontCache.bc gaussrand.bc $(EM_CFLAGS) $(EM_LDFLAGS) -o index.js --preload-file Crimson-Roman.ttf