// boxtest.c - example code, disks bouncing in a box
#include <stdio.h>
#include "macsdl2.h"
#include "myrand4a.h"
#if __EMSCRIPTEN__
#include <emscripten/emscripten.h>
#endif
#include <SDL2/SDL.h>

#define REAL float		// WebAssembly likes floats
#define SIZEX 1024		// size of display screen
#define SIZEY 512
#define XFAC (REAL)SIZEX/SIZEY
#define NDISKS 10000	// maximum number disks...
#define DELT .001		// time step

REAL X[NDISKS],Y[NDISKS],VX[NDISKS],VY[NDISKS];

static int quit = 0;
int ncount = 0;			// total number of updates
REAL xscale,yscale;		// display scaling parameters
unsigned int _seed = 1;	// for random number generator...  
int framecount = 0;
int ndisks = 100;		// number of disks displayed
REAL radius = .05;		// radius of disks
int slowctl = 10;		// number of updates per display
REAL bx1= -1,bx2 = 1;	// box size
REAL by1= -1,by2 = 1;

void initdisks(int nd,float rad)
{
for(int i=0;i<123;i++) myrand();	// scramble random number generator
for(int i=0;i<ndisks;i++)			// initial positions
	{
	float x,y;
	x = f1_rand();		// random between -1 and 1
	x *= (1 - radius);	// keep entire disk inside box
	X[i] = x;
	y = f1_rand();
	y *= (1 - radius);
	Y[i] = y;
	}
for(int i=0;i<ndisks;i++)			// initial velocities
	{
//	float theta;
	float theta = ang_rand();
	VX[i] = 1.0 * cos(theta);
	VY[i] = 1.0 * sin(theta);
	}
}

void update()
{
register REAL dt=DELT;
for(int i=0;i<ndisks;i++)
	{
	X[i] += dt*VX[i];
	Y[i] += dt*VY[i];
	if(X[i] < bx1+radius) VX[i] = -VX[i];	// perform bounces...
	if(X[i] > bx2-radius) VX[i] = -VX[i];
	if(Y[i] < by1+radius) VY[i] = -VY[i];
	if(Y[i] > by2-radius) VY[i] = -VY[i];
	}
++ncount;
}

void display()
{
color(WHITE);			// sets color
clear();				// clears Plane1[] with color
color(BLACK);
// draw surrounding box
move(bx1,by1); pendn(); move(bx2,by1); move(bx2,by2); move(bx1,by2); move(bx1,by1); penup();
// draw disks
for(int i=0;i<ndisks;i++)
	{
	color(RED);
	circf(X[i],Y[i],radius);
	color(BLACK);
	circ(X[i],Y[i],radius);
	}

drawtexture();			// copies Plane1[] to texture
drawnum(ncount,820,0);	// draws on texture
iPrintf(10,0,"%d disks in the box",ndisks);
iPrintf(10,100,"n - change ndisks");
iPrintf(10,150,"s - change speed");
iPrintf(10,200,"x - exit");
drawbuffer();			// copies texture to display
++framecount;
}

void main_tick() {
SDL_Event event;

while (SDL_PollEvent(&event)) {
	switch (event.type) {
		case SDL_KEYDOWN: {
			switch (event.key.keysym.sym) {
				case SDLK_ESCAPE:
				case SDLK_x: {
					quit = 1;
					break; }
				case SDLK_n: {
					ndisks *= 10;
					if(ndisks > 10000) ndisks = 10;
					initdisks(ndisks,radius);
					}
				case SDLK_o: {
					float dt;
					static int oldncount = 0, oldticks = 0;
					int newticks = SDL_GetTicks();
					dt = .001 * (newticks - oldticks);
					oldticks = newticks;
					printf("elapsed time %.2f\t",dt);
					printf("frames per second %.2f\n",framecount/dt);
					framecount = 0;
					printf("iterations per second %.2f\n",(ncount - oldncount)/dt);
					oldncount = ncount;
					break; }
				case SDLK_s: {
					if(SDL_GetModState() & KMOD_SHIFT)
						{
						if(SDL_GetModState() & KMOD_ALT) slowctl /= 2;
						else slowctl /= 10;
						if(slowctl == 0) slowctl = 1;
						}
					else
						{
						if(SDL_GetModState() & KMOD_ALT) slowctl *= 2;
						else slowctl *= 10;
						if(slowctl > 10000) slowctl = 1;
						}
					printf("slowctl %d\n",slowctl);
					break; }
				}
			}
		}
	}
for(int i = 0;i<slowctl;i++)
	update();
display();
}

void main_loop()
{
#if __EMSCRIPTEN__
	emscripten_set_main_loop(main_tick, -1, 1);
#else
    while (0 == quit){main_tick();}
#endif
}

int main(int argc, char *argv[])
{
if(argc == 2 && argv[1] != NULL) ndisks = atoi(argv[1]);
prefsize(SIZEX,SIZEY);
initgraph("hi there!");
yscale = 1.2;
xscale = yscale * XFAC;
scale(-xscale,xscale,-yscale,yscale);
initdisks(ndisks,radius);
labelsize(30);
labelcolor(0,0,255,255);
initlabel();
main_loop();
exitgraph();
return 0;
}
