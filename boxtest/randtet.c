//	randtet.c - fills square with random dots
//	looks like it can generate 20 * 10^6 random numbers, and draw the dots, at 10 Hz
//	on the 2.9 Mhz macbook pro.  that's 200,000,000 dots per second!
//	10.7 Hz macbook pro
//	 7.6 Hz safari (localhost)
//	 9.4 Hz firefox (local)
//	 6.9 Hz chrome (local)
//	 7.0 Hz brave (local)
#include <stdio.h>
#include "macsdl2.h"
#include "myrand4a.h"
#if __EMSCRIPTEN__
#include <emscripten/emscripten.h>
#endif
#include <SDL2/SDL.h>

#define REAL float
#define SIZEX 1024
#define SIZEY 512
#define XFAC (REAL)SIZEX/SIZEY

static int quit = 0;
int ncount = 0;
REAL xscale,yscale;		// display scaling parameters
unsigned int _seed = 1;	// for random number generator...
int framecount = 0;
int nits = 100000;		// dots per display frame
int slowctl = 1;		// slow display flag...

void display()
{
drawtexture();			// copies Plane1[] to texture
drawnum(ncount,820,0);	// draws on texture
iPrintf(10,0,"%d dots per frame",nits);
drawbuffer();			// copies texture to display
++framecount;
}

void update()
{
color(WHITE);			// sets color
clear();				// clears Plane1[] with color
color(BLACK);
for(int i=0;i<nits;i++)
	{
	REAL x = f1_rand();
	REAL y = f1_rand();
	dotat(x,y);
//	dotat(0,0);			// even this is slower in browser?
	}
++ncount;
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
if(argc == 2 && argv[1] != NULL) nits = atoi(argv[1]);
prefsize(SIZEX,SIZEY);
initgraph("hi there!");
yscale = 1.2;
xscale = yscale * XFAC;
scale(-xscale,xscale,-yscale,yscale);
labelsize(30);
labelcolor(0,0,255,255);
initlabel();
main_loop();
exitgraph();
return 0;
}
