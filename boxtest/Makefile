# fill TG in with program you want to compile
# for lldb use  % make DBG=-g
# for webassembly, % make em
# original compiler around 10% slower than new one

#TG = randtet
TG = boxtest

DBG = -O3
CFLAGS = $(DBG) -flto
LDFLAGS = -flto
LIBS= -lm -lSDL2 -lSDL2_ttf
CC = cc
#CC = /Library/Developer/CommandLineTools_slim/usr/bin/cc

OBJ = $(TG).o macsdl2.o SDL_FontCache.o
#OBJ = $(TG).o macsdl2.o

$(TG):	$(OBJ)
	$(CC) -o $(TG) $(OBJ) $(LIBS)

clean:
	rm -f *.o *.bc *.js *.wasm *.data $(TG)

tar:
	tar czvf $(TG).tgz *

gzip:
	gzip -k index.wasm

EM_CFLAGS = -c -s WASM=1 -O3 -flto
EM_LDFLAGS= -s USE_SDL=2 -s USE_SDL_TTF=2 --no-heap-copy -s TOTAL_MEMORY=294912000 -s ALLOW_MEMORY_GROWTH=1

em:
	emcc $(TG).c $(EM_CFLAGS) -o $(TG).bc
	emcc macsdl2.c $(EM_CFLAGS) -o macsdl2.bc
	emcc SDL_FontCache.c $(EM_CFLAGS) -o SDL_FontCache.bc
	emcc $(TG).bc macsdl2.bc SDL_FontCache.bc $(EM_LDFLAGS) -o index.js --preload-file Crimson-Roman.ttf
