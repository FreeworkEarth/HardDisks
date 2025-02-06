#include "globals.h"
#include <stdio.h>
#include <stdlib.h>		// for exit()
#include <math.h>
#include <string.h>
#include <SDL2/SDL.h>
#include "glgraph.h"
#include "rstuff.h"
#include "myrand4.h"

// SDL Window & Renderer
SDL_Renderer* _renderer = NULL;
SDL_Window* _window = NULL;



// Particle properties
REAL X[NSPHERES], Y[NSPHERES], VX[NSPHERES], VY[NSPHERES];
REAL Radius[NSPHERES], LastBumpTime[NSPHERES], BumpInterval[NSPHERES];

// Histograms
unsigned long Histo[NENTS], Histo1[NENTS], Histo2[NENTS], Histowall[NENTS], HistoW[NENTS];

// Cells and neighborhoods
long Cells[NCELLS * CSIZE];
long DiskNbrs[NSPHERES * CSIZE];

// Piston properties
REAL xpist1 = XW1, xpist2 = XW2, vpist1 = 0, vpist2 = 0, dxpist = .001;

// Simulation parameters
REAL Ktot = 0, KE0 = 0, KE1 = 0, KE2 = 0;
long bumpnum = 0, wallbump = 0, ncount = 0, nruns = 0;
REAL total_time = 0;

// Histogram controls
long histctl = 0;
REAL vmax = 3.0;

// Simulation modes
int dispctl = 0, radctl = 0, pistctl = 0, modectl = 1;

// Collision tracking
struct bump bumps[MAXBUMPS];