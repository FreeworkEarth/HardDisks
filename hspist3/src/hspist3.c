#include <GLFW/glfw3.h>  // OpenGL support
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>
#include <string>
#include <iostream>

// Define constants with constexpr for better type safety
constexpr int NSPHERES = 1000;   // Number of spheres in the simulation
constexpr int NENTS = 100;       // Number of histogram entries
constexpr int NCELLS = 50;       // Number of cells in the grid
constexpr int CSIZE = 10;        // Cell size for the grid
constexpr int MAXBUMPS = 100;    // Maximum number of bumps allowed in the simulation

// Modernize global variables
unsigned long _seed = SEED;  // Random seed for simulation

// Arrays for position, velocity, and properties of the spheres
std::vector<REAL> X(NSPHERES), Y(NSPHERES), VX(NSPHERES), VY(NSPHERES), DM2(NSPHERES);
std::vector<REAL> Radius(NSPHERES), LastBumpTime(NSPHERES), BumpInterval(NSPHERES);

// Histograms for various simulation outputs
std::vector<unsigned long> Histo(NENTS), Histo1(NENTS), Histo2(NENTS), Histowall(NENTS), HistoW(NENTS);

// Grid and neighbor information
std::vector<long> Cells(NCELLS * CSIZE);  // Grid cells storing sphere info
std::vector<long> DiskNbrs(NSPHERES * CSIZE);  // Neighbors of spheres in grid
std::vector<unsigned long> Vout(SIZEX * SIZEY);  // Output display image

// Simulation state variables and settings
REAL xw1 = XW1, xw2 = XW2, yw1 = YW1, yw2 = YW2;  // Boundaries for the simulation space
REAL rad1, rad2;  // Disk radii
long nrad1, nrad2;  // Number of disks with rad1 and rad2
long nspheres = NBEGIN;  // Number of spheres in the simulation
long seed = SEED;  // Seed for random number generation
long bumpnum = 0;  // Number of bumps encountered
long bumptot = 0;  // Total number of bumps
long wallbump = 0;  // Number of bumps against the wall
long rowl, nrows, ncells;  // Indices for neighborhood calculations
long slowctl = 1;  // Slow display flag
long keyslowctl = 1;  // Flag for slow display in keyboard input
long Slowctl = 1;  // Screen refresh slow display flag
int dispctl = 0;  // Display control
int Dispctl = 0;  // Secondary display control
int radctl = 0;  // Flag for adjusting disk radii
long ncount = 0;  // Iteration counter
long nruns = 0;  // Number of simulation runs
REAL total_time = 0;  // Total simulation time
REAL xpist1 = XW1;  // Position of the left piston
REAL vpist1 = 0;  // Velocity of the left piston
REAL xpist2 = XW2;  // Position of the right piston
REAL vpist2 = 0;  // Velocity of the right piston
REAL dxpist = 0.001;  // Step size for piston motion
REAL Ktot;  // Total kinetic energy of the system
REAL KE0, KE1, KE2;  // Kinetic energy at different points in the simulation
long histctl = 0;  // Histogram display flag
REAL vmax = 3.0;  // Maximum velocity for histogram
REAL Phend;  // Henderson approximation for hard disk pressure
long trailctl = 0;  // Trail control flag
char *mypath;  // Program path
char pathName[1024];  // Path name constructed with loadpath()
long datacnt = 0;  // Number of times histogram totals are summed
REAL binwidth;  // Width of histogram bins
long nbins;  // Number of bins used for histograms
int forward = 1;  // Flag for trajectory direction
long ndone;  // Number of iterations completed per run
int pistctl = 0;  // Piston protocol control
REAL llo, lhi, rlo, rhi;  // Starting and ending positions of the left and right pistons
REAL lvel, rvel;  // Velocities of the left and right pistons
REAL pdelay;  // Time delay for second (right) piston motion
REAL dt1, dt2;  // Time intervals for piston motion
int modectl = 1;  // Demo mode flag
int ppctl = 0;  // Program protocol control
REAL dx;  // Distance of piston stroke

// Structure to store bump information
struct bump {
    long disknum;  // Disk number
    long iden;  // Identifier for the bump (positive for second disk, negative for walls)
    REAL fcenx;  // Force center x-coordinate
    REAL fceny;  // Force center y-coordinate
};
std::vector<bump> bumps(MAXBUMPS);  // List of bumps