# makefile
# fill TG in with program you want to compile
# for gdb use  % make DBG=-g

#TG= lorentzrevn
#TG= hspist_mb
#TG= hspist_nsinai_j
#TG= hspist_nsinai_jld
#TG= hspist_diff_ld
#TG= hspist_dual
#TG= hspist_diff_all
TG= hspist2

DBG = -O3
CFLAGS = $(DBG) -Wmost -Wno-deprecated-declarations -m32 -isysroot /Developer/SDKs/MacOSX10.6.sdk
LIBS= -framework OpenGL -framework GLUT -framework Carbon -m32 -isysroot /Developer/SDKs/MacOSX10.6.sdk

OBJ = $(TG).o gaussrand.o glgraph.o kbd.o vWindow.o fades.o

$(TG):	$(OBJ)
	cc -o $(TG) $(OBJ) $(LIBS)

HEADERS = $(TG).h rstuff.h
$(OBJ): $(HEADERS)

clean:
	rm *.o $(TG)

tar:
	tar czvf $(TG).tgz *
