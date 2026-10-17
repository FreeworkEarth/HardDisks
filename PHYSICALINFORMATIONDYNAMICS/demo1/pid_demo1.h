/* pid_demo1.h -- the run of pid_demo1 (pid_demo1.c), shared by the headless path and the window (pid_view1.c, brief sec. 9).
   The window draws and never owns state: it changes the run only through run1_advance, the function the headless path uses,
   so a window run writes the same pstate_*.txt, density_*.txt and measure.csv as a headless one (make test, test 11). */
#ifndef PID_DEMO1_H
#define PID_DEMO1_H
#include <stdint.h>
#include "edmd_gen3.h"
#include "pid.h"

#define PID1_PXS  24.0     /* px per sigma */
#define PID1_TU   24.0     /* internal time units per sigma-time */
#define PID1_GRID 2.0      /* density grid cell [sigma] */

typedef struct { double x1, y1, r1, x2, y2, r2, ux, uy, h; int whole; } Pid1Cap;   /* whole: 1 = disk 1 alone, 2 = disk 2 alone */

/* the window's options (bin/pid_view1 only); view: 0 disks, 1 network (alpha 0, hull and graph), 2 both; overlay: 0 off,
   1 particle density, 2 particle kinetic energy, 3 structural material; show (if given, it replaces what view shows): the
   PID1_SHOW_* bits */
enum { PID1_SHOW_DISKS = 1, PID1_SHOW_HULL = 2, PID1_SHOW_GRAPH = 4, PID1_SHOW_PARTICLES = 8, PID1_SHOW_FOOD = 16 };
typedef struct { double alpha, speed; int view, overlay, exit_at_end, nshot, show, show_given; double shot_at[16]; } Pid1View;

typedef struct {
    /* options */
    const char* state; const char* out;
    long N; double T, R, cell, write_every, kbt; uint64_t seed; long audit_every; int have_seed;
    int window; Pid1View v;
    int argc; char** argv;
    /* the Demo 0 world (sigma) and the walls (px) */
    World w; char fsha[65];
    int nseg, narc; EDMD3_Seg* seg; EDMD3_Arc* arc;
    Pid1Cap* cap; int ncap;
    int alive; double disk_area, lens;
    /* the state */
    EDMD3* S; double ke0, eta; long trials;
    /* the density grid: per cell the structural material (fixed), and the particle count and kinetic energy (run1_bin) */
    int gx, gy; double* mat; long* cnt; double* kec;
    /* the stops */
    long nstop, kstop; int done, fatal; double* Jprev; long esc_total; double min_dist; double cpu0;
} Run1;

int  run1_setup(Run1* r, int argc, char** argv, int window);   /* 0: go on; else the exit code (the reason on stderr) */
int  run1_advance(Run1* r, double t_sigma);   /* advance to min(t_sigma, T), writing every stop on the way; 1 while T is ahead */
double run1_time(const Run1* r);              /* sigma-time of the state */
void run1_bin(Run1* r, const EDMD_Particle* P);   /* cnt and kec of the density grid from the particles P */
void run1_note_file(const char* name);        /* a file written into --out (listed with its sha256 in manifest.txt) */
int  run1_report(Run1* r);                    /* summary, checks, manifest.txt; the exit code (the state stays readable) */
void run1_free(Run1* r);                      /* frees what run1_setup made, also after a failed setup */
int  run1_finish(Run1* r);                    /* run1_report, then run1_free */
int  view1_run(Run1* r);                      /* pid_view1.c: the window */

#endif
