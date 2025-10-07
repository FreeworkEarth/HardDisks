# fill TG in with program you want to compile
# for lldb use  % make DBG=-g
# for webassembly, % make em
# original compiler around 10% slower than new one

# Target program
TG = boxtest

# Compiler and flags
CC = gcc                    
DBG = -O3                   # Use -g for debug, or -O3 for optimized builds
CFLAGS = $(DBG) -flto -I/opt/homebrew/include -D_THREAD_SAFE
LDFLAGS = -flto -L/opt/homebrew/lib
LIBS = -lm -lSDL2 -lSDL2_ttf

# Object files to compile
OBJ = $(TG).o macsdl2.o SDL_FontCache.o

# Default target
all: $(TG)

$(TG): $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS) $(LIBS)

# Clean up build artifacts
clean:
	rm -f *.o *.bc *.js *.wasm *.data $(TG)

# Archive the project
tar:
	tar czvf $(TG).tgz *

# Optional gzip target (for WebAssembly builds)
gzip:
	gzip -k index.wasm

# WebAssembly build with Emscripten
EM_CFLAGS = -c -s WASM=1 -O3 -flto
EM_LDFLAGS = -s USE_SDL=2 -s USE_SDL_TTF=2 --no-heap-copy \
             -s TOTAL_MEMORY=294912000 -s ALLOW_MEMORY_GROWTH=1

em:
	emcc $(TG).c $(EM_CFLAGS) -o $(TG).bc
	emcc macsdl2.c $(EM_CFLAGS) -o macsdl2.bc
	emcc SDL_FontCache.c $(EM_CFLAGS) -o SDL_FontCache.bc
	emcc $(TG).bc macsdl2.bc SDL_FontCache.bc $(EM_LDFLAGS) -o index.js \
		--preload-file Crimson-Roman.ttf
