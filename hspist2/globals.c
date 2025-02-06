#include "globals.h"
#include "rstuff.h"

// SDL Window & Renderer
SDL_Renderer* _renderer = NULL;
SDL_Window* _window = NULL;

// Random seed and simulation parameters
unsigned long int _seed = SEED;  // Defined in globals.h or another header
long nspheres = NBEGIN;  // Number of spheres (must be smaller than NSPHERES)

// Particle properties
REAL X[NSPHERES], Y[NSPHERES], VX[NSPHERES], VY[NSPHERES], DM2[NSPHERES];
REAL Radius[NSPHERES], LastBumpTime[NSPHERES], BumpInterval[NSPHERES];

// Histograms
unsigned long Histo[NENTS], Histo1[NENTS], Histo2[NENTS], Histowall[NENTS], HistoW[NENTS];

// Cells and neighborhoods
long Cells[NCELLS * CSIZE];
long DiskNbrs[NSPHERES * CSIZE];

// Output display
unsigned long Vout[SIZEX * SIZEY];  // Output display image

// Grid and cell parameters
REAL xw1 = XW1, xw2 = XW2, yw1 = YW1, yw2 = YW2;
REAL rad1, rad2;  // Disk radii
long nrad1, nrad2;  // Number of disks with RAD1, RAD2...
long rowl, nrows, ncells;  // Neighborhood indices

// Simulation counters and flags
long bumpnum = 0, bumptot = 0, wallbump = 0;
long ncount = 0, nruns = 0;  // Iteration and run counters
REAL total_time = 0;  // Total simulation time

// Piston properties
REAL xpist1 = XW1, vpist1 = 0, xpist2 = XW2, vpist2 = 0, dxpist = .001;

// Kinetic energy and control flags
REAL Ktot = 0, KE0 = 0, KE1 = 0, KE2 = 0;  // Total and kinetic energies
long histctl = 0;  // Histogram control
REAL vmax = 3.0;  // Maximum for velocity histogram

// Other simulation control flags
int dispctl = 0, Dispctl = 0, radctl = 0;  // Display control flags
long slowctl = 1, keyslowctl = 1, Slowctl = 1;  // Slow display flags
long trailctl = 0;  // Trail control flag

// File and path variables
char *mypath = NULL;
char pathName[1024];  // Path name constructed with loadpath()

// Histogram settings
long datacnt = 0;  // Number of times histogram totals are summed
REAL binwidth = 0.1;  // Histogram bin width
long nbins = 50;  // Number of histogram bins

// Simulation mode flags
int forward = 1;  // Flag for trajectory direction
long ndone = 0;  // Number of ncounts per run
int pistctl = 0;  // Piston protocol control flag
REAL llo, lhi, rlo, rhi;  // Left, right piston starting and ending positions
REAL lvel = 0.1, rvel = 0.1;  // Left, right piston velocities
REAL pdelay = 0.5;  // Time delay for second (right) piston
REAL dt1 = 0.01, dt2 = 0.01;  // Time for first and second piston motion
int modectl = 1;  // Demo mode flag
int ppctl = 0;  // Program protocol control flag
REAL dx = 0.01;  // Distance of piston stroke

// Collision tracking structure
struct bump bumps[MAXBUMPS];  // Array of bump interactions