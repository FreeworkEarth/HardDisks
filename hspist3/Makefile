# Program to compile
TG = 00ALLINONE

# Debug flags or optimization flags
DBG = -O3
CFLAGS = $(DBG) -fsanitize=address -g -Wmost -Wno-deprecated-declarations -I/opt/homebrew/include -I/opt/homebrew/include/SDL2 -D_THREAD_SAFE
LIBS = -lGL -L/opt/homebrew/lib -lSDL2 -L/opt/homebrew/opt/mesa/lib -lSDL2_ttf -lGLEW -framework CoreGraphics -fsanitize=address
LDFLAGS = -L/opt/homebrew/opt/mesa/lib -fsanitize=address

# Compile the main program
$(TG): $(TG).c
	cc -o $(TG) $(TG).c $(CFLAGS) $(LIBS)

# Clean up the build
clean:
	rm -f $(TG)

# Create a tarball for the project
tar:
	tar czvf $(TG).tgz *