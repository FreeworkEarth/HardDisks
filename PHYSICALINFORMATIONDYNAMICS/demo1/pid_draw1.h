/* pid_draw1.h -- the drawing of Demo 1's windows (pid_view1.c, pid_studio.c): the Demo 0 network, the hull, the particles and
   the coarse-grid overlay, in sigma through one viewport. It only reads: the World, the Run1 (pid_demo1.h) and the particles. */
#ifndef PID_DRAW1_H
#define PID_DRAW1_H
#include <SDL.h>
#include "pid_demo1.h"

typedef struct { double ox, oy, scale, H; } Pid1Vp;   /* sigma -> screen: x' = ox + x scale, y' = oy + (H - y) scale (y up) */

/* what a scene shows: PID1_SHOW_* bits; overlay 0 off, 1 particle density, 2 particle kinetic energy, 3 structural material;
   alpha 0..255 the disks' fill; growth 1: Demo 0's colours of a growing network (live tips, the new front, dead ends) */
typedef struct { int show, overlay, alpha, growth; } Pid1Scene;

void pd_to_screen(const Pid1Vp* v, double x, double y, float* sx, float* sy);
void pd_fill_disk(SDL_Renderer* R, float cx, float cy, float rad);
void pd_ring(SDL_Renderer* R, float cx, float cy, float rad);
void pd_ramp(double s, Uint8* r, Uint8* g, Uint8* b);          /* dark blue .. yellow, s in [0, 1] */
void pd_speed_colour(double s, Uint8* r, Uint8* g, Uint8* b);  /* the particles: blue, white, red */
/* the box, then (bottom to top) the overlay (needs run), the structural disks, the food patches, the graph, the particles (needs
   run), the hull (w->hull if valid); tex: the overlay's texture, created on first use and kept by the caller */
void pd_scene(SDL_Renderer* R, const Pid1Vp* v, const World* w, Run1* run, const Pid1Scene* sc, SDL_Texture** tex);

#endif
