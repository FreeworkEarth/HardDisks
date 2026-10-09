#ifndef EDMD_GEN3_H
#define EDMD_GEN3_H

/* ##CHRIS 2026-10-08 (261012 sec. 4.7 design, sec. 4.7.1 amendments; milestone M1, the core engine):
   generation-3 event-driven MD for 2D hard disks, a constant-work-per-event engine.
     - square cells of width w >= the diameter (default 32 px), disk-local coordinates (position relative to the
       lower-left corner of the disk's cell);
     - cell-crossing events that carry the new cell (never recomputed from a rounded position);
     - a time stamp per disk (lazy advancing); sync_all() for every reader of all positions;
     - predictions only among the 9 neighbouring cells; outer-wall events only for disks in cells that can reach a wall;
     - one binary heap with lazy invalidation, ordered by (t, type, a, b): the only arbiter of simultaneous events;
     - a floating time origin: every EDMD3_ORIGIN_SHIFT internal units all disks are synchronised and the shift is
       subtracted from every stored time (exact: all of them are >= the shift then);
     - the safety nets and the validator of sec. 4.7.1 (h) with counters, and two read-only audits (contact distance at
       every executed event; the live heap against a brute-force all-pairs prediction from the synchronised state).
   M1 has no dividers, pistons, heat bath or semipermeable gates: edmd3_create() refuses such parameters (M2 adds the
   divider and piston bands). Units as in edmd.c: px and internal time units (24 px = 1 sigma; 24 units = 1 sigma-time). */

#include <stddef.h>
#include <stdint.h>
#include "edmd.h"

#ifdef __cplusplus
extern "C" {
#endif

#define EDMD3_DEFAULT_CELL_PX   32.0      /* 1.33 diameters at R = 12 px (sec. 4.7.1, amendment a) */
#define EDMD3_ORIGIN_SHIFT      8192.0    /* 2^13 internal units = 341.3 sigma-time */

typedef struct EDMD3 EDMD3;

typedef struct {
    long ev_pair, ev_wall, ev_cross;     /* executed events by kind */
    long ev_stale;                       /* invalidated events dropped at pop */
    /* safety nets: each fires 0 times in a correct run */
    long overlap_repair;   /* a pair predicted while overlapping beyond rounding and approaching: scheduled at once */
    long wall_overdue;     /* a wall predicted at or past its face while approaching: scheduled at once */
    long past_event;       /* a live event popped with t < now: skipped */
    long clamp_repair;     /* a disk found outside the box by the post-event check: put back and reflected */
    long cell_repair;      /* a disk found outside its cell beyond tolerance: re-filed */
    long grid_escape;      /* a crossing event that would leave the grid: dropped */
    long stagnation;       /* more than the guard's number of live events at one time: the run is stopped (fatal) */
    long contact_now;      /* pairs at contact within rounding and approaching: scheduled at once (not a repair) */
    /* validator (sec. 4.7.1, h) */
    long local_checks, local_findings;   /* after every event: the event's disks against their 9 cells and the walls */
    long full_checks, full_findings;     /* every check interval and at every sync_all: every pair in 9 cells, every wall */
    double local_worst, full_worst;      /* most negative surface gap seen [px] (0 if none) */
    double cross_residual_max;           /* max |local coordinate - cell boundary| when a crossing re-files a disk [px] */
    long origin_shifts, syncs, heap_compactions, heap_max;
} EDMD3_Health;

typedef struct {
    long audits;                                               /* audited states */
    long pair_cmp, pair_missing, pair_extra, pair_dt, pair_dt_rel, pair_deferred;
    long wall_cmp, wall_missing, wall_extra, wall_dt, wall_dt_rel, wall_deferred;
    long cross_cmp, cross_missing, cross_extra, cross_dt, cross_dt_rel, cross_dup;   /* cross_dup: a second live crossing of one disk */
    long dup_disagree, cell_inconsistent;
    /* a true event not yet scheduled ("deferred") must come after the next crossing of one of its disks (the crossing
       that will make it eligible); these count deferred events EARLIER than that crossing: 0 in a correct engine */
    long pair_deferred_early, wall_deferred_early;
    double max_dt, max_rel;     /* max |t_bruteforce - t_heap|, and that over the horizon, over matched events */
} EDMD3_Audit;

/* lifecycle; cell_px <= 0 selects EDMD3_DEFAULT_CELL_PX. NULL on error, with the reason in err. */
EDMD3* edmd3_create(const EDMD_Params* prm, double cell_px, char* err, size_t errlen);
void   edmd3_destroy(EDMD3* S);
/* load positions and velocities (absolute px) at absolute time t_abs and schedule everything. 0 on error (an overlap
   or a disk outside the box beyond the validator's tolerance), with the reason in err. */
int    edmd3_load(EDMD3* S, const EDMD_Particle* P, double t_abs, char* err, size_t errlen);

/* process every event with time <= t_abs; disks stay lazy (call edmd3_particles or edmd3_sync_all to read them) */
double edmd3_advance_to(EDMD3* S, double t_abs);
void   edmd3_sync_all(EDMD3* S);                      /* advance every disk to the current time (O(N)), full check */
const EDMD_Particle* edmd3_particles(EDMD3* S);       /* synchronised absolute copies (calls edmd3_sync_all) */
double edmd3_time(const EDMD3* S);
int    edmd3_count(const EDMD3* S);
double edmd3_cell_px(const EDMD3* S);
int    edmd3_fatal(const EDMD3* S, const char** msg);   /* nonzero once the stagnation guard stopped the run */

/* telemetry */
const EDMD3_Health* edmd3_health(const EDMD3* S);
uint64_t edmd3_event_hash(const EDMD3* S);           /* FNV-1a over (t, type, a, b) of every executed event */
double edmd3_kinetic_energy(EDMD3* S);               /* sum v^2 / 2, unit mass */
void   edmd3_reset_virial(EDMD3* S);
double edmd3_compressibility_Z(EDMD3* S);            /* 1 + W / (2 KE t) over the window, as edmd_compressibility_Z */
long   edmd3_virial_pair_events(const EDMD3* S);
double edmd3_wall_impulse(const EDMD3* S, int wall); /* L, R, B, T */
long   edmd3_wall_events(const EDMD3* S, int wall);
void   edmd3_set_check_interval(EDMD3* S, double units); /* full validator cadence (default 24 = 1 sigma-time); <= 0 off */

/* read-only audits: neither changes any stored number, so outputs are identical with and without them */
void   edmd3_set_contact_audit(EDMD3* S, int on);
long   edmd3_contact_audit_stats(const EDMD3* S, double max_gap_px[2]);   /* [0] pairs, [1] outer walls */
void   edmd3_set_schedule_audit(EDMD3* S, long every);                   /* after every k-th executed event; 0 off */
void   edmd3_schedule_audit_now(EDMD3* S);
const EDMD3_Audit* edmd3_schedule_audit_stats(const EDMD3* S);

#ifdef __cplusplus
}
#endif
#endif /* EDMD_GEN3_H */
