# Program to compile
TG = 00ALLINONE

# Debug vs Release
DEBUG_FLAGS = -g -fsanitize=address -DDEBUG
RELEASE_FLAGS = -O3 -ffast-math -march=native
KISSFFT_SRCS = kissfft/kiss_fft.c kissfft/kiss_fftr.c

# Common include flags
INCLUDE = -I/opt/homebrew/include -I/opt/homebrew/include/SDL2 -D_THREAD_SAFE

# Choose one of these manually (uncomment one)
CFLAGS_DEBUG = -g -fsanitize=address -DDEBUG $(INCLUDE)
CFLAGS_RELEASE = -O3 -ffast-math -march=native $(INCLUDE)

debug:
	cc -o $(TG) $(TG).c $(KISSFFT_SRCS) $(CFLAGS_DEBUG) $(LIBS)

release:
	cc -o $(TG) $(TG).c $(KISSFFT_SRCS) $(CFLAGS_RELEASE) $(LIBS)

LIBS = -lGL -L/opt/homebrew/lib -lSDL2 -L/opt/homebrew/opt/mesa/lib -lSDL2_ttf -lGLEW -framework CoreGraphics
LDFLAGS = -L/opt/homebrew/opt/mesa/lib

# Compile the main program
$(TG): $(TG).c
	cc -o $(TG) $(TG).c $(KISSFFT_SRCS) $(CFLAGS) $(LIBS)

# Clean up the build
clean:
	rm -f $(TG)

# Create a tarball for the project
tar:
	tar czvf $(TG).tgz *