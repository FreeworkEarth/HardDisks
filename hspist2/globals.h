#ifndef GLOBALS_H
#define GLOBALS_H

#include <SDL2/SDL.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>    // for exit()
#include <string.h>
#include <SDL2/SDL.h>
#include "rstuff.h"
#include "myrand4.h"

// SDL related globals
extern SDL_Renderer* _renderer;
extern SDL_Window* _window;

// Random seed and simulation parameters
extern unsigned long int _seed;
extern long nspheres;

// Particle arrays
extern REAL X[], Y[], VX[], VY[], Radius[];  // Positions, velocities, and radii of particles
extern REAL DM2[], LastBumpTime[], BumpInterval[];  // More particle properties
extern long bumpnum, bumptot, wallbump;  // Bump counters
extern unsigned long Histo[], Histo1[], Histo2[], Histowall[], HistoW[];  // Histogram arrays

// Grid and cell parameters
extern long Cells[], DiskNbrs[];  // Cell-related variables
extern REAL xw1, xw2, yw1, yw2;  // Simulation grid boundaries
extern long nrad1, nrad2;  // Radius-related parameters
extern long rowl, nrows, ncells;  // Grid row/column counts
extern long slowctl, keyslowctl, Slowctl;  // Slow control settings
extern int dispctl, Dispctl, radctl;  // Display control settings

// Time and velocity parameters
extern REAL total_time;  // Total simulation time
extern REAL xpist1, vpist1, xpist2, vpist2;  // Piston positions and velocities
extern REAL dxpist;  // Piston displacement
extern REAL Ktot, KE0, KE1, KE2;  // Kinetic energy values
extern long histctl;  // Histogram control
extern REAL vmax;  // Maximum velocity

// Particle interaction parameters
extern REAL Phend;  // End of phase value
extern long trailctl;  // Trail control
extern char *mypath;  // Path for file saving/loading
extern char pathName[1024];  // Path name string
extern long datacnt;  // Data count
extern REAL binwidth;  // Bin width for histograms
extern long nbins;  // Number of bins for histograms
extern int forward;  // Forward simulation control flag
extern long ndone;  // Number of steps completed
extern int pistctl;  // Piston control

// Piston velocity and delay
extern REAL llo, lhi, rlo, rhi;  // Piston limits
extern REAL lvel, rvel;  // Piston velocities
extern REAL pdelay;  // Piston delay
extern REAL dt1, dt2;  // Time steps for pistons

// Mode control
extern int modectl;  // Mode control flag
extern int ppctl;  // Additional control flag
extern REAL dx;  // Displacement for updates

// Structure for particle bump interactions
struct bump {
    long disknum;
    long iden;    // positive for second disk, negative for walls, force centers...
    REAL fcenx;   // Force center coordinates...
    REAL fceny;
};

extern struct bump bumps[MAXBUMPS];  // Array for bump interactions

// Function declarations for core simulation tasks
void update_particles(); // Function to update particle positions
void handle_collisions(); // Function to detect and resolve collisions

#endif // GLOBALS_H