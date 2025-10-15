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

    /* Optional static vertical divider (slab) inside the box */
    int    has_divider;    /* 1 to enable, 0 to disable */
    double divider_x;      /* center x position (0..boxW) */
    double divider_thickness; /* slab thickness (>=0) */
    double divider_mass;   /* wall mass in particle-mass units (>=0), 0 => infinite mass */
    double divider_vx;     /* divider velocity (x), piecewise-constant between events */
    /* Optional pistons as moving faces */
    int    has_pistonL;    /* left piston (right face) enabled */
    double pistonL_x;      /* face x position */
    double pistonL_vx;     /* face velocity */
    double pistonL_mass;   /* mass in particle units (0 => infinite) */

    int    has_pistonR;    /* right piston (left face) enabled */
    double pistonR_x;      /* face x position */
    double pistonR_vx;     /* face velocity */
    double pistonR_mass;   /* mass in particle units (0 => infinite) */
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

/* Advance considering only particle–particle (AB) collisions up to t_target.
   Ignores/does not schedule any boundary events (container walls, divider,
   pistons). Particles free-fly between AB events. Useful for hybrid schemes
   where faces are handled by a separate CCD pass. */
double edmd_advance_pp_only_to(EDMD* S, double t_target);

/* optional: rebuild grid + reschedule all (after bulk manual changes) */
void   edmd_reschedule_all(EDMD* S);

/* Clear and reschedule ONLY AB events for the current state. Does not build
   wall events and does not clamp/push particles versus boundaries. */
void   edmd_reschedule_all_pp_only(EDMD* S);

/* configure/update a static divider slab (call then reschedule) */
void   edmd_config_divider(EDMD* S, int enabled, double cx, double thickness);

/* set divider dynamic properties (mass, velocity). mass=0 => infinite mass */
void   edmd_set_divider_motion(EDMD* S, double mass, double vx);

/* configure pistons (faces) */
void   edmd_config_pistons(EDMD* S,
                           int hasL, double xL, double vxL, double mL,
                           int hasR, double xR, double vxR, double mR);

/* analytics helpers */
double edmd_total_kinetic_energy(const EDMD* S, double mass);

#ifdef __cplusplus
}
#endif
#endif /* EDMD_H */
