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
   ##CHRIS 2026-10-08 (261012 sec. 4.7.4, sec. 4.7.6; milestone M2): dividers and pistons as BANDS (sec. 4.7, item 2e).
   A divider (a vertical slab: held = mass 0 and velocity 0, driven = mass 0 and velocity != 0, free = mass > 0, spring =
   mass > 0 and k > 0) or a piston (one vertical face; mass 0 = prescribed velocity) carries a band of cell columns
   outside which no disk can touch it before the band's expiry; only band disks carry its events. The rounding
   tolerances are derived from the time resolution (amendment b); a momentum and an energy ledger with their rounding
   scales are kept (amendment d). Refused: the heat bath, species and semipermeable gates, switched-off pair collisions,
   a divider of thickness <= 0 (edmd.c ignores such a divider). Units as in edmd.c: px and internal time units
   (24 px = 1 sigma; 24 units = 1 sigma-time); unit disk mass. */

#include <stddef.h>
#include <stdio.h>
#include <stdint.h>
#include "edmd.h"

#ifdef __cplusplus
extern "C" {
#endif

#define EDMD3_DEFAULT_CELL_PX   32.0      /* 1.33 diameters at R = 12 px (sec. 4.7.1, amendment a) */
#define EDMD3_ORIGIN_SHIFT      8192.0    /* 2^13 internal units = 341.3 sigma-time */
#define EDMD3_TOL_K             4.0       /* the factor K of the rounding tolerances (the derived bound is 2.5; sec. 4.7.6) */

typedef struct EDMD3 EDMD3;

typedef struct {
    long ev_pair, ev_wall, ev_cross;     /* executed events by kind */
    long ev_stale;                       /* invalidated events dropped at pop */
    /* safety nets: each fires 0 times in a correct run */
    long overlap_repair;   /* a pair predicted while overlapping beyond rounding (c < -c_tol) and approaching: scheduled at once */
    long wall_overdue;     /* a wall predicted past its face beyond rounding (gap < -tol_face) while approaching: scheduled at once */
    long past_event;       /* a live event popped with t < now: skipped */
    long clamp_repair;     /* a disk found outside the box by the post-event check: put back and reflected */
    long cell_repair;      /* a disk found outside its cell beyond tolerance: re-filed */
    long grid_escape;      /* a crossing event that would leave the grid: dropped */
    long stagnation;       /* more than the guard's number of live events at one time: the run is stopped (fatal) */
    long contact_now;      /* pairs at contact within rounding (-c_tol <= c <= 0) and approaching: scheduled at once (not a repair) */
    /* validator (sec. 4.7.1, h) */
    long local_checks, local_findings;   /* after every event: the event's disks against their 9 cells, the walls and the bodies */
    long full_checks, full_findings;     /* every check interval and at every sync_all: every pair in 9 cells, every wall and body */
    double local_worst, full_worst;      /* most negative surface gap seen [px] (0 if none) */
    double cross_residual_max;           /* max |local coordinate - cell boundary| when a crossing re-files a disk [px] */
    long origin_shifts, syncs, heap_compactions, heap_max;
    /* M2 */
    long ev_div, ev_piston, ev_band;     /* executed divider, piston and band events */
    long obj_overlap_repair;             /* a divider or piston face predicted with gap < -tol_face while approaching: at once (safety net) */
    long obj_contact_now;                /* the same within rounding (-tol_face <= gap <= 0): at once (not a repair) */
    long wall_contact_now;               /* an outer wall at contact within rounding while approaching: at once (not a repair) */
    long body_findings;                  /* a divider or piston outside the box or out of order (validator; must be 0) */
    double contact_c_min;                /* most negative c = |r|^2 - d^2 among contact_now pairs [px^2] (0 if none) */
    double obj_contact_gap_min;          /* most negative face gap among obj_contact_now events [px] (0 if none) */
} EDMD3_Health;

/* schedule-audit classes (amendment e: per-class maxima over ALL matched events, no floor) */
enum { EDMD3_CLS_PAIR = 0, EDMD3_CLS_WALL = 1, EDMD3_CLS_CROSS = 2, EDMD3_CLS_DIV = 3, EDMD3_CLS_PISTON = 4, EDMD3_NCLS = 5 };

typedef struct {
    long audits;                                               /* audited states */
    long pair_cmp, pair_missing, pair_extra, pair_dt, pair_dt_rel, pair_deferred;
    long wall_cmp, wall_missing, wall_extra, wall_dt, wall_dt_rel, wall_deferred;
    long cross_cmp, cross_missing, cross_extra, cross_dt, cross_dt_rel, cross_dup;   /* cross_dup: a second live crossing of one disk */
    long dup_disagree, cell_inconsistent;
    /* a true event not yet scheduled ("deferred") must come after the next crossing of one of its disks (the crossing
       that will make it eligible); these count deferred events EARLIER than that crossing: 0 in a correct engine */
    long pair_deferred_early, wall_deferred_early;
    double max_dt, max_rel;     /* max |t_bruteforce - t_heap|, and that over the horizon (|dt| > 1e-9 only), over matched events */
    /* M2: the divider and piston classes. Deferred: the disk's column is outside the body's band; such an event must
       come after the disk's next crossing or the band's expiry, whichever is first. Bands: a moving body without its
       live BAND event, a second live one, a band that misses a column its body can reach before the expiry. */
    long div_cmp, div_missing, div_extra, div_dt, div_dt_rel, div_deferred, div_deferred_early;
    long pis_cmp, pis_missing, pis_extra, pis_dt, pis_dt_rel, pis_deferred, pis_deferred_early;
    long band_missing, band_extra, band_short;
    /* amendment e: per class, over ALL matched events (no 1e-9 floor): max |dt| and max |dt| / horizon (horizon =
       max(t_bruteforce, t_heap) - now >= |dt|, so the ratio is at most 1; it is 1 for an event due now in one of the two) */
    double cls_max_dt[EDMD3_NCLS], cls_max_rel[EDMD3_NCLS];
    double cls_rel_hz[EDMD3_NCLS];   /* the horizon of the event that gave cls_max_rel */
    double cls_dt_hz[EDMD3_NCLS];    /* the horizon of the event that gave cls_max_dt */
} EDMD3_Audit;

/* amendment d: the momentum and energy ledgers. P = total momentum of every body of finite mass (the disks, free and
   spring dividers, pistons of finite mass); J = the impulses delivered to them from outside: the outer walls, dividers
   and pistons of mass 0 (held or driven), the spring anchors, and velocity or mass changes made through the API. Then
   P - P0 = J up to rounding; scale_P = u x (the sum of |result| over the rounded operations of the dynamics and the
   ledger), u = 2^-53: a first-order bound of the rounding (sec. 4.7.6). The same for the mechanical energy E (kinetic +
   spring) and the work W of driven bodies and API changes. Per axis: [0] x, [1] y. */
typedef struct {
    double P[2], P0[2], J[2], scale_P[2];
    double E, E0, W, scale_E;
    double J_wall[4];                    /* signed impulse of each outer wall on the disks: L, R (x), B, T (y) */
    double J_div[EDMD_MAX_DIVIDERS];     /* signed x impulse of each divider of mass 0 on the disks (0 while free) */
    double J_spring[EDMD_MAX_DIVIDERS];  /* signed x impulse of each spring anchor on its divider */
    double J_piston[2];                  /* signed x impulse of each piston of mass 0 on the disks (L, R) */
    double J_api;                        /* x momentum given to bodies of finite mass by API velocity or mass changes */
} EDMD3_Ledger;

/* the tolerances in force (amendment b, M2 acceptance 5): derived, not tuned */
typedef struct {
    double u_time;      /* ulp(EDMD3_ORIGIN_SHIFT) = 2^-39 internal units: the spacing of event times in [2^13, 2^14) */
    double E_bound;     /* mechanical energy (kinetic + spring) at load plus the work of driven bodies and API changes since */
    double m_min;       /* min(1, the lightest body of finite mass) */
    double v_ref;       /* sqrt(2 E_bound (1 + 1/m_min)) + the fastest body of mass 0: no relative speed exceeds it [px/unit] */
    double K;           /* EDMD3_TOL_K */
    double c_tol;       /* K * 2 d * v_ref * u_time [px^2]: a pair contact within rounding vs an overlap repair */
    double tol_face;    /* K * v_ref * u_time + 8 ulp(box width) [px]: the same for outer walls, divider faces and pistons */
    double band_margin; /* closed-interval margin of the band test [px] (= the validator's wall tolerance) */
    double tol_pair, tol_wall, tol_cell;   /* the validator's tolerances (experiment_validation.c) */
} EDMD3_Tol;

/* lifecycle; cell_px <= 0 selects EDMD3_DEFAULT_CELL_PX. NULL on error, with the reason in err.
   Dividers and pistons are read from prm as edmd.c reads them; their time stamps are set by edmd3_load. */
EDMD3* edmd3_create(const EDMD_Params* prm, double cell_px, char* err, size_t errlen);
void   edmd3_destroy(EDMD3* S);
/* load positions and velocities (absolute px) at absolute time t_abs and schedule everything. 0 on error (an overlap
   or a disk outside the box or inside a body beyond the validator's tolerance), with the reason in err. */
int    edmd3_load(EDMD3* S, const EDMD_Particle* P, double t_abs, char* err, size_t errlen);

/* process every event with time <= t_abs; disks stay lazy (call edmd3_particles or edmd3_sync_all to read them) */
double edmd3_advance_to(EDMD3* S, double t_abs);
void   edmd3_sync_all(EDMD3* S);                      /* advance every disk to the current time (O(N)), full check */
const EDMD_Particle* edmd3_particles(EDMD3* S);       /* synchronised absolute copies (computed, not stored) */
double edmd3_time(const EDMD3* S);
int    edmd3_count(const EDMD3* S);
double edmd3_cell_px(const EDMD3* S);
int    edmd3_fatal(const EDMD3* S, const char** msg);   /* nonzero once the stagnation guard stopped the run */

/* M2: bodies. A change acts at the current time (hold, release, a step of a piston protocol); mass 0 = infinite. The
   readers compute the state at the current time without storing it. 0 on a bad index or an inactive body. */
int    edmd3_set_divider_motion(EDMD3* S, int d, double mass, double vx);
int    edmd3_set_piston_motion(EDMD3* S, int side, double mass, double vx);   /* side 0 = left piston, 1 = right */
int    edmd3_divider_state(const EDMD3* S, int d, double* x, double* vx);
int    edmd3_piston_state(const EDMD3* S, int side, double* x, double* vx);
double edmd3_divider_impulse(const EDMD3* S, int d, int face);   /* x impulse received from disks on the left (face 0) or right (1) */
long   edmd3_divider_events(const EDMD3* S, int d, int face);
double edmd3_piston_impulse(const EDMD3* S, int side);          /* x impulse the piston received from disks */
long   edmd3_piston_events(const EDMD3* S, int side);
double edmd3_work_divider(const EDMD3* S, int d);    /* as edmd_work_divider_i: mass 0: energy given to disks; else the divider's KE gain */
double edmd3_work_piston(const EDMD3* S, int side);  /* as edmd_work_pistonL/R */
void   edmd3_ledger(EDMD3* S, EDMD3_Ledger* L);       /* reads only */
int    edmd3_health_clean(const EDMD3* S);            /* the run flag: 1 iff every safety net, validator and body count is 0 */
void   edmd3_tolerances(const EDMD3* S, EDMD3_Tol* T);

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
/* ##CHRIS 2026-10-09 (M3): edmd.c's gated event log (t_sigma, kind, u_wall, v_before, v_after, dE, dp) per wall, divider-face and
   piston collision; f NULL = off (the default); time_scale converts internal time to sigma-time (the driver passes 24) */
void   edmd3_set_event_log(EDMD3* S, FILE* f, double time_scale);

/* read-only audits: neither changes any stored number, so outputs are identical with and without them */
void   edmd3_set_contact_audit(EDMD3* S, int on);
long   edmd3_contact_audit_stats(const EDMD3* S, double max_gap_px[2]);   /* [0] pairs, [1] outer walls */
long   edmd3_contact_audit_stats4(const EDMD3* S, double max_gap_px[4]);  /* [0] pairs, [1] outer walls, [2] divider faces, [3] pistons */
void   edmd3_set_schedule_audit(EDMD3* S, long every);                   /* after every k-th executed event; 0 off */
void   edmd3_set_schedule_audit_bodies(EDMD3* S, int on);              /* M2: also after every BAND event and API body change */
void   edmd3_schedule_audit_now(EDMD3* S);
const EDMD3_Audit* edmd3_schedule_audit_stats(const EDMD3* S);

#ifdef __cplusplus
}
#endif
#endif /* EDMD_GEN3_H */
