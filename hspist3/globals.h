#ifndef GLOBALS_H
#define GLOBALS_H

// SDL Window & Renderer
#include <SDL2/SDL.h>

// Simulation Constants
#define SEED 123456
#define NSPHERES 1000
#define NBEGIN 500
#define NENTS 100
#define NCELLS 100
#define CSIZE 10
#define MAXBUMPS 500
#define SIZEX 800
#define SIZEY 600
#define XW1 0.0
#define XW2 10.0
#define YW1 0.0
#define YW2 10.0

// Particle and simulation types
typedef double REAL;  // Change this to 'float' or 'double' as necessary

// External declarations
extern SDL_Renderer* _renderer;
extern SDL_Window* _window;

// Particle properties
extern REAL X[NSPHERES], Y[NSPHERES], VX[NSPHERES], VY[NSPHERES], DM2[NSPHERES];
extern REAL Radius[NSPHERES], LastBumpTime[NSPHERES], BumpInterval[NSPHERES];

// Histograms
extern unsigned long Histo[NENTS], Histo1[NENTS], Histo2[NENTS], Histowall[NENTS], HistoW[NENTS];

// Cells and neighborhoods
extern long Cells[NCELLS * CSIZE];
extern long DiskNbrs[NSPHERES * CSIZE];

// Output display
extern unsigned long Vout[SIZEX * SIZEY];

// Grid and cell parameters
extern REAL xw1, xw2, yw1, yw2;
extern REAL rad1, rad2;  // Disk radii
extern long nrad1, nrad2;  // Number of disks with RAD1, RAD2...
extern long nspheres;  // Number of disks, must be smaller than NSPHERES

// Simulation counters and flags
extern long bumpnum, bumptot, wallbump;
extern long ncount, nruns;  // Iteration and run counters
extern REAL total_time;  // Total simulation time

// Piston properties
extern REAL xpist1, vpist1, xpist2, vpist2, dxpist;

// Kinetic energy and control flags
extern REAL Ktot, KE0, KE1, KE2;  // Total and kinetic energies
extern long histctl;  // Histogram control
extern REAL vmax;  // Maximum for velocity histogram

// Other simulation control flags
extern int dispctl, Dispctl, radctl;  // Display control flags
extern long slowctl, keyslowctl, Slowctl;  // Slow display flags
extern long trailctl;  // Trail control flag

// File and path variables
extern char *mypath;
extern char pathName[1024];  // Path name constructed with loadpath()

// Histogram settings
extern long datacnt;  // Number of times histogram totals are summed
extern REAL binwidth;  // Histogram bin width
extern long nbins;  // Number of histogram bins

// Simulation mode flags
extern int forward;  // Flag for trajectory direction
extern long ndone;  // Number of ncounts per run
extern int pistctl;  // Piston protocol control flag
extern REAL llo, lhi, rlo, rhi;  // Left, right piston starting and ending positions
extern REAL lvel, rvel;  // Left, right piston velocities
extern REAL pdelay;  // Time delay for second (right) piston
extern REAL dt1, dt2;  // Time for first and second piston motion
extern int modectl;  // Demo mode flag
extern int ppctl;  // Program protocol control flag
extern REAL dx;  // Distance of piston stroke

// Collision tracking structure
struct bump {
    long disknum;
    long iden;  // Positive for second disk, negative for walls, force centers...
    REAL fcenx;  // Force center coordinates...
    REAL fceny;
};

extern struct bump bumps[MAXBUMPS];

#endif  // GLOBALS_H