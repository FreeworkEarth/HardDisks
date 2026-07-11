#ifndef EDMD_H
#define EDMD_H

/* Minimal, C99, event-driven molecular dynamics for 2D hard disks.
   Static rectangular box [0..boxW] x [0..boxH]. No pistons here (yet). */

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define EDMD_MAX_DIVIDERS 32

typedef struct {
    double boxW, boxH;     /* container size */
    double radius;         /* disk radius (all equal) */
    int    N;              /* number of particles */

    int    max_events_hint;/* initial heap capacity hint (optional) */
    double cell_size;      /* 0 => auto (~2.5*radius) */

    /* Optional vertical divider slabs (internal walls) inside the box */
    int    divider_count;  /* number of active dividers (0..EDMD_MAX_DIVIDERS) */
    double divider_x[EDMD_MAX_DIVIDERS];        /* center x positions (0..boxW) */
    double divider_thickness[EDMD_MAX_DIVIDERS];/* slab thickness per divider (>=0) */
    double divider_mass[EDMD_MAX_DIVIDERS];     /* wall mass in particle-mass units (>=0), 0 => infinite mass */
    double divider_vx[EDMD_MAX_DIVIDERS];       /* divider velocities (x), piecewise-constant between events */
    /* Optional harmonic springs for divider motion (event-driven; requires TOI root-finding). */
    double divider_k[EDMD_MAX_DIVIDERS];        /* spring constant (k). 0 => no spring (constant-velocity divider) */
    double divider_xeq[EDMD_MAX_DIVIDERS];      /* equilibrium position (x_eq) for spring (same units as divider_x) */
    /* Optional pistons as moving faces */
    int    has_pistonL;    /* left piston (right face) enabled */
    double pistonL_x;      /* face x position */
    double pistonL_vx;     /* face velocity */
    double pistonL_mass;   /* mass in particle units (0 => infinite) */

    int    has_pistonR;    /* right piston (left face) enabled */
    double pistonR_x;      /* face x position */
    double pistonR_vx;     /* face velocity */
    double pistonR_mass;   /* mass in particle units (0 => infinite) */

    /* ##CHRIS: Disable particle-particle collisions (0 = off, 1 = on, default 1) */
    int    pp_collisions_enabled;

    /* ##CHRIS: Heat bath parameters for outer walls */
    int    heatbath_enabled;        /* 1 to enable heat bath on outer walls */
    double heatbath_temperature;    /* Target temperature for heat bath */
    int    thermal_wall_mode;       /* 0=gradual, 1=base MB, 2=adaptive MB (default 2) */
    double mb_overshoot_factor;     /* Overshoot factor for adaptive mode (default 2.0) */
    double stability_window_percent;/* Stability window (default 0.025 = 2.5%) */
    double particle_mass;           /* Particle mass (for MB sampling, default 1.0) */
    double kB;                      /* Boltzmann constant (for MB sampling, default 1.0) */

    /* Optional per-particle species labels (external pointer, length N).
       Used to implement Szilard-style semipermeable dividers in EDMD. */
    const int* species;             /* NULL => all species 0 */

    /* Optional semipermeable divider behavior (per divider, independent).
       The divider is treated as a hard slab by default. When a gate is enabled,
       collisions for the "allowed" direction are treated as pass-through (no impulse).

       divider_gate_mode:
         0 = disabled (hard divider)
         1 = single-species one-way gate (see divider_gate_species/target_side)
         2 = dual-species gate: species 0 targets LEFT, species 1 targets RIGHT
         3 = single-species one-way *block* gate:
               - divider_gate_species[d] is the BLOCKED species
               - divider_gate_target_side[d] is the side the blocked species is allowed to move toward
               - all other species pass through in both directions
         4 = single-species one-way gate with temperature-triggered leak for the opposite species:
               - divider_gate_species[d] is the normally allowed species
               - divider_gate_target_side[d] is the allowed direction target side
               - the opposite species is blocked unless its source-side temperature exceeds the
                 allowed-species temperature on the target side by divider_gate_hot_ratio[d]

       Optional hot-particle leak:
         - if divider_gate_speed_ref[d] > 0 and divider_gate_hot_ratio[d] > 0, then
           blocked particles may pass when their speed exceeds
           divider_gate_hot_ratio[d] * divider_gate_speed_ref[d]
         - for mode 1 this applies to the opposite species moving toward target side
         - for mode 3 this applies to the blocked species moving away from target side

       divider_gate_target_side (mode 1): 0=LEFT, 1=RIGHT
       The gate blocks motion *away* from the target side, i.e.
         target LEFT  => block left->right (EV_DL), allow right->left (EV_DR)
         target RIGHT => block right->left (EV_DR), allow left->right (EV_DL)
    */
    int divider_gate_mode[EDMD_MAX_DIVIDERS];
    int divider_gate_species[EDMD_MAX_DIVIDERS];
    int divider_gate_target_side[EDMD_MAX_DIVIDERS];
    double divider_gate_hot_ratio[EDMD_MAX_DIVIDERS];
    double divider_gate_speed_ref[EDMD_MAX_DIVIDERS];
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

/* configure/update divider slabs (call then reschedule) */
void   edmd_config_dividers(EDMD* S, int count, const double* cx, const double* thickness);

/* set divider dynamic properties (mass, velocity). mass=0 => infinite mass */
void   edmd_set_divider_motions(EDMD* S, int count, const double* mass, const double* vx);

/* set harmonic spring properties (k, x_eq). k<=0 disables the spring for that divider. */
void   edmd_set_divider_springs(EDMD* S, int count, const double* k, const double* xeq);

/* legacy single-divider helpers (index 0) */
void   edmd_config_divider(EDMD* S, int enabled, double cx, double thickness);
void   edmd_set_divider_motion(EDMD* S, double mass, double vx);
void   edmd_set_divider_spring(EDMD* S, double k, double xeq);

/* configure pistons (faces) */
void   edmd_config_pistons(EDMD* S,
                           int hasL, double xL, double vxL, double mL,
                           int hasR, double xR, double vxR, double mR);

/* analytics helpers */
double edmd_total_kinetic_energy(const EDMD* S, double mass);

/* ##CHRIS: Heat bath helper - compute gas temperature from kinetic energy */
double edmd_gas_temperature(const EDMD* S, double mass, double kB);

/* Work accumulation for moving boundaries (cumulative since create/reset) */
/* Total work across all dividers (sum). */
double edmd_work_divider(const EDMD* S);
/* Work for a specific divider index (0-based). */
double edmd_work_divider_i(const EDMD* S, int divider_index);
double edmd_work_pistonL(const EDMD* S);
double edmd_work_pistonR(const EDMD* S);
/* Heat exchanged with the heat bath (outer thermal walls), cumulative since create/reset.
   Positive means energy injected into the gas by the bath. */
double edmd_heat_bath(const EDMD* S);
void   edmd_reset_work(EDMD* S);

/* Resolve any particles currently overlapping divider slabs (push + reflect). */
void   edmd_divider_resolve_overlaps(EDMD* S);

#ifdef __cplusplus
}
#endif
#endif /* EDMD_H */
