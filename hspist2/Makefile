# makefile
# fill TG in with program you want to compile
# for gdb use  % make DBG=-g

# Program to compile
TG = hspist2

# Debug flags or optimization flags
DBG = -O3
CFLAGS = $(DBG) -Wmost -Wno-deprecated-declarations \
         -I/opt/homebrew/Cellar/freeglut/3.6.0/include

# Linking with FreeGLUT and OpenGL for ARM architecture
LIBS = -L/opt/homebrew/Cellar/freeglut/3.6.0/lib -lGL -lGLU -lfreeglut

# Object files
OBJ = $(TG).o gaussrand.o glgraph.o kbd.o vWindow.o fades.o

# Compile the main program
$(TG): $(OBJ)
	cc -o $(TG) $(OBJ) $(LIBS)

# Include the headers in dependencies
HEADERS = $(TG).h rstuff.h
$(OBJ): $(HEADERS)

# Clean up the build
clean:
	rm -f *.o $(TG)

# Create a tarball for the project
tar:
	tar czvf $(TG).tgz *