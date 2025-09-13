#ifndef EDMD_H
#define EDMD_H

/* Minimal, C99, event-driven molecular dynamics for 2D hard disks.
   Static rectangular box [0..boxW] x [0..boxH]. No pistons here (yet). */

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double boxW, boxH;     /* container size */
    double radius;         /* disk radius (all equal) */
    int    N;              /* number of particles */

    int    max_events_hint;/* initial heap capacity hint (optional) */
    double cell_size;      /* 0 => auto (~2.5*radius) */
} EDMD_Params;

typedef struct {
    double x, y;           /* position */
    double vx, vy;         /* velocity */
    int    coll_count;     /* increments on each collision (event invalidation) */
} EDMD_Particle;

typedef struct EDMD EDMD;

/* lifecycle */
EDMD*  edmd_create(const EDMD_Params* prm);
void   edmd_destroy(EDMD* S);
void   edmd_init_random_gas(EDMD* S, unsigned long long seed);

/* access */
double                edmd_time(const EDMD* S);
const EDMD_Params*    edmd_params(const EDMD* S);
const EDMD_Particle*  edmd_particles(const EDMD* S);
int                   edmd_count(const EDMD* S);

/* advance to absolute time t_target (process events chronologically) */
double edmd_advance_to(EDMD* S, double t_target);

/* optional: rebuild grid + reschedule all (after bulk manual changes) */
void   edmd_reschedule_all(EDMD* S);

/* analytics helpers */
double edmd_total_kinetic_energy(const EDMD* S, double mass);

#ifdef __cplusplus
}
#endif
#endif /* EDMD_H */