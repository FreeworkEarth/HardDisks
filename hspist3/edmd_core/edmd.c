#include "edmd.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <stdio.h>

#ifndef EDMD_ADVANCE_MAX_EVENTS
#define EDMD_ADVANCE_MAX_EVENTS 250000L
#endif
#ifndef EDMD_ADVANCE_MAX_STAGNANT_EVENTS
#define EDMD_ADVANCE_MAX_STAGNANT_EVENTS 5000L
#endif
static long g_edmd_avalanche_warning_count = 0;

/* ------------------------ internal types ------------------------ */

/* event kind: AB = particle-particle, WL/WR/WB/WT = walls (left/right/bottom/top) */
typedef enum { EV_AB=0, EV_WL, EV_WR, EV_WB, EV_WT, EV_DL, EV_DR, EV_PL, EV_PR } EvType;
static const char* ev_name(EvType t){
    switch(t){
        case EV_AB: return "AB";
        case EV_WL: return "WL";
        case EV_WR: return "WR";
        case EV_WB: return "WB";
        case EV_WT: return "WT";
        case EV_DL: return "DL";
        case EV_DR: return "DR";
        case EV_PL: return "PL";
        case EV_PR: return "PR";
        default: return "?";
    }
}

/* scheduled event at absolute time t */
typedef struct {
    double t;
    int    a, b;        /* if b==-1 => wall event for particle a */
    int    ca, cb;      /* coll_count snapshots (invalidate stale events fast) */
    EvType type;
} Event;

/* ##CHRIS: ------------------- event-history tracer -------------------
   Debug-only instrumentation for hunting missed particle-particle collisions.
   Disabled unless edmd_debug_set_watch() is called, and then it only costs one
   distance evaluation per processed event.

   It keeps a ring buffer of the most recent events and, after every position
   advance, measures the surface gap of the watched pair. The first time that gap
   goes negative it dumps the ring, so you see the exact event sequence that
   produced the overlap instead of guessing from the post-mortem state. */
#define EDMD_TRACE_CAP 4096
typedef struct {
    long   seq;
    double t;
    EvType type;
    int    a, b, ca, cb;
    int    ca_now, cb_now;   /* live coll_count at pop time (stale iff != ca/cb) */
    int    accepted;         /* 1 = resolved, 0 = discarded as stale */
    double gap;              /* watched-pair surface gap after the position jump */
} EdmdTraceRec;

static int          g_trace_a = -1, g_trace_b = -1;
static int          g_trace_cap = 0;
static EdmdTraceRec g_trace[EDMD_TRACE_CAP];
static long         g_trace_seq = 0;
static int          g_trace_head = 0, g_trace_count = 0, g_trace_fired = 0;

/* uniform grid cell (list of particle indices) */
typedef struct {
    int *idx;
    int  count, cap;
} Cell;

/* min-heap of events (priority by smallest t) */
typedef struct {
    Event* data;
    int    n, cap;
} Heap;

/* main state */
struct EDMD {
    EDMD_Params prm;
    double t;
    EDMD_Particle* P;
    /* Work accumulation on moving boundaries (kinetic energy change) */
    double work_divider[EDMD_MAX_DIVIDERS];
    double work_pistonL;
    double work_pistonR;
    /* Heat exchange at thermal (heat bath) outer walls (kinetic energy change) */
    double heat_bath;
    long forced_advance_count;

    /* ##CHRIS: grid_build()'s safety clamp is a real (if rare) boundary bounce: it can
       move a particle and reverse its velocity. That used to happen silently, leaving the
       particle with a stale event set and no scheduled events for its new direction --
       which is how two particles could end up on an unscheduled collision course and pass
       straight through each other. grid_build now records every particle it mutates here
       and bumps its coll_count; callers must reschedule these before returning. */
    int*  clamped;            /* indices mutated by the most recent grid_build() */
    int   clamped_count;
    long  clamp_repair_count; /* total clamps applied (telemetry) */
    long  overlap_repair_count; /* PP collisions scheduled for already-overlapping pairs */
    long  wall_overdue_count;   /* wall collisions that were already due when scheduled */
    /* ##CHRIS: equilibrium pressure by the collisional virial. For hard disks the
       only contribution to the virial is the impulse exchanged at contact:
           W = sum_collisions  m * |dv_n| * sigma
       and, since 2D kinetic energy is KE = N kB T (kB=1, m=1),
           Z = P/(rho kB T) = 1 + W / (2 * KE * t).
       Accumulated here rather than post-hoc because it needs the pre-collision
       normal velocity, which is gone once resolve_ab() has run. Pure telemetry:
       nothing here feeds back into the dynamics. */
    double virial_accum;        /* sum of m*|dv_n|*sigma over PAIR collisions only */
    double virial_t0;           /* sim time when accumulation last reset */
    long   virial_pair_events;  /* number of pair collisions counted */

    /* ##CHRIS: SECOND, INDEPENDENT pressure estimator -- direct momentum flux on
       the four stationary outer walls. These impulses are deliberately NOT added
       into virial_accum: doing so would double-count the pressure, since the pair
       virial above is already the complete hard-disk collisional virial. They are
       accumulated per wall so that P_x and P_y can be formed separately, which
       makes isotropy (Z_x == Z_y) an independent check rather than an assumption.
       Only the elastic reflection path contributes: a heat-bath wall resamples the
       velocity and is not a momentum-conserving reflector, so its "impulse" is not
       a mechanical pressure. */
    double wall_impulse[4];     /* |dp| on L,R,B,T outer walls */
    long   wall_events[4];      /* collision counts, same order */
    long   wall_thermal_events; /* heat-bath bounces skipped (must be 0 for EOS) */

    Cell* grid; int gw, gh;    /* grid width/height */
    double cell_size;

    Heap heap;
};

static inline int edmd_trace_active(const EDMD* S){
    return g_trace_cap > 0 && !g_trace_fired &&
           g_trace_a >= 0 && g_trace_b >= 0 &&
           g_trace_a < S->prm.N && g_trace_b < S->prm.N;
}

/* surface gap of the watched pair: <0 means overlapping */
static double edmd_trace_gap(const EDMD* S){
    const double dx = S->P[g_trace_b].x - S->P[g_trace_a].x;
    const double dy = S->P[g_trace_b].y - S->P[g_trace_a].y;
    return sqrt(dx*dx + dy*dy) - 2.0*S->prm.radius;
}

static void edmd_trace_dump(const EDMD* S, const char* where, double gap){
    const EDMD_Particle* A = &S->P[g_trace_a];
    const EDMD_Particle* B = &S->P[g_trace_b];
    const double dx = B->x - A->x, dy = B->y - A->y;
    const double dist = sqrt(dx*dx + dy*dy);
    const double vn = (dist > 0.0)
        ? ((B->vx - A->vx)*dx + (B->vy - A->vy)*dy) / dist : NAN;

    fprintf(stderr,
        "\n[EDMD-TRACE] watched pair (%d,%d) went overlapping at t=%.17g via %s\n"
        "  gap=%.17g  dist=%.17g  sigma=%.17g  v_rel_normal=%.17g\n"
        "  p%d: x=%.17g y=%.17g vx=%.17g vy=%.17g coll_count=%d\n"
        "  p%d: x=%.17g y=%.17g vx=%.17g vy=%.17g coll_count=%d\n"
        "  last %d events (oldest first); '*' marks an event touching the pair:\n",
        g_trace_a, g_trace_b, S->t, where,
        gap, dist, 2.0*S->prm.radius, vn,
        g_trace_a, A->x, A->y, A->vx, A->vy, A->coll_count,
        g_trace_b, B->x, B->y, B->vx, B->vy, B->coll_count,
        g_trace_count);

    for (int k = 0; k < g_trace_count; ++k) {
        const int idx = (g_trace_head - g_trace_count + k + EDMD_TRACE_CAP) % EDMD_TRACE_CAP;
        const EdmdTraceRec* r = &g_trace[idx];
        const int touches = (r->a == g_trace_a || r->a == g_trace_b ||
                             r->b == g_trace_a || r->b == g_trace_b);
        fprintf(stderr,
            "   %c #%-8ld t=%-22.17g %-2s a=%-4d b=%-4d ca=%d/%d cb=%d/%d %-8s gap=%.17g\n",
            touches ? '*' : ' ', r->seq, r->t, ev_name(r->type), r->a, r->b,
            r->ca, r->ca_now, r->cb, r->cb_now,
            r->accepted ? "RESOLVED" : "stale",
            r->gap);
    }
    fprintf(stderr, "[EDMD-TRACE] end of history\n\n");
    fflush(stderr);
}

static void edmd_trace_record(const EDMD* S, const Event* e, int accepted){
    if (!edmd_trace_active(S)) return;
    EdmdTraceRec* r = &g_trace[g_trace_head];
    r->seq      = g_trace_seq++;
    r->t        = e->t;
    r->type     = e->type;
    r->a        = e->a;
    r->b        = e->b;
    r->ca       = e->ca;
    r->cb       = e->cb;
    r->ca_now   = (e->a >= 0 && e->a < S->prm.N) ? S->P[e->a].coll_count : -1;
    r->cb_now   = (e->type == EV_AB && e->b >= 0 && e->b < S->prm.N) ? S->P[e->b].coll_count : -1;
    r->accepted = accepted;
    r->gap      = edmd_trace_gap(S);
    g_trace_head = (g_trace_head + 1) % EDMD_TRACE_CAP;
    if (g_trace_count < g_trace_cap) g_trace_count++;
}

/* Call after any position advance. Fires (once) when the pair first overlaps. */
static void edmd_trace_check(const EDMD* S, const char* where){
    if (!edmd_trace_active(S)) return;
    const double gap = edmd_trace_gap(S);
    /* same scale as the experiment validator: 1e-6 of a diameter */
    const double tol = fmax(1e-7, 1e-6 * 2.0 * S->prm.radius);
    if (gap < -tol) {
        g_trace_fired = 1;           /* set before dump: dump must not re-enter */
        edmd_trace_dump(S, where, gap);
    }
}

/* ##CHRIS: log every scheduling decision made for the watched pair, so a collision that
   is never scheduled is distinguishable from one that is scheduled and then lost. */
static void edmd_trace_schedule_attempt(const EDMD* S, int i, int j, int ok, double t_abs){
    const double rx = S->P[j].x - S->P[i].x, ry = S->P[j].y - S->P[i].y;
    const double vx = S->P[j].vx - S->P[i].vx, vy = S->P[j].vy - S->P[i].vy;
    const double b  = rx*vx + ry*vy;
    const double rr = rx*rx + ry*ry, vv = vx*vx + vy*vy;
    const double sig = 2.0*S->prm.radius;
    const double c  = rr - sig*sig;
    fprintf(stderr,
        "[EDMD-TRACE] schedule_ab(%d,%d) at t=%.17g -> %s%s%.17g  "
        "b=%.17g c=%.17g vv=%.17g disc=%.17g gap=%.17g\n",
        i, j, S->t, ok ? "SCHEDULED t_abs=" : "no collision",
        ok ? "" : " ", ok ? t_abs : 0.0,
        b, c, vv, b*b - vv*c, sqrt(rr) - sig);
    fflush(stderr);
}

void edmd_debug_set_watch(int a, int b, int history){
    g_trace_a = a;
    g_trace_b = b;
    if (history < 0) history = 0;
    if (history > EDMD_TRACE_CAP) history = EDMD_TRACE_CAP;
    g_trace_cap   = history;
    g_trace_head  = 0;
    g_trace_count = 0;
    g_trace_seq   = 0;
    g_trace_fired = 0;
}

/* ------------------------ small helpers ------------------------ */

static inline int clampi(int v, int lo, int hi) { return v<lo?lo:(v>hi?hi:v); }
static inline int clamp_dividers(int v) { return v < 0 ? 0 : (v > EDMD_MAX_DIVIDERS ? EDMD_MAX_DIVIDERS : v); }
static inline int divider_has_spring(const EDMD* S, int d){
    if (!S) return 0;
    if (S->prm.divider_count <= 0) return 0;
    if (d < 0 || d >= clamp_dividers(S->prm.divider_count)) return 0;
    return (S->prm.divider_k[d] > 0.0) && (S->prm.divider_mass[d] > 0.0);
}

/* heap utils */
static void heap_init(Heap* H, int cap){
    if(cap < 64) cap = 64;
    H->data = (Event*)malloc((size_t)cap*sizeof(Event));
    H->n = 0; H->cap = cap;
}
static void heap_free(Heap* H){ free(H->data); H->data=NULL; H->n=H->cap=0; }
static void heap_swap(Event* a, Event* b){ Event t=*a; *a=*b; *b=t; }
static void heap_push(Heap* H, Event e){
    if(H->n>=H->cap){
        H->cap = H->cap + H->cap/2 + 64;
        H->data = (Event*)realloc(H->data, (size_t)H->cap*sizeof(Event));
    }
    int i = H->n++;
    H->data[i] = e;
    while(i>0){
        int p=(i-1)>>1;
        if(H->data[p].t <= H->data[i].t) break;
        heap_swap(&H->data[p], &H->data[i]); i=p;
    }
}
static int heap_pop(Heap* H, Event* out){
    if(H->n==0) return 0;
    *out = H->data[0];
    H->data[0] = H->data[--H->n];
    int i=0;
    for(;;){
        int l=2*i+1, r=l+1, s=i;
        if(l<H->n && H->data[l].t < H->data[s].t) s=l;
        if(r<H->n && H->data[r].t < H->data[s].t) s=r;
        if(s==i) break;
        heap_swap(&H->data[i], &H->data[s]); i=s;
    }
    return 1;
}

/* ------------------------ grid (uniform spatial hashing) ------------------------ */

static void cell_reserve(Cell* c, int need){
    if(need <= c->cap) return;
    int ncap = c->cap? c->cap*2 : 32;
    while(ncap < need) ncap*=2;
    c->idx = (int*)realloc(c->idx, (size_t)ncap*sizeof(int));
    c->cap = ncap;
}

static void grid_free(EDMD* S){
    if(!S->grid) return;
    int n = S->gw * S->gh;
    for(int i=0;i<n;i++) free(S->grid[i].idx);
    free(S->grid); S->grid=NULL;
}

/* rebuild whole grid from current positions (simple & robust) */
static void grid_build(EDMD* S){
    if(!S->grid){
        S->grid = (Cell*)calloc((size_t)S->gw*S->gh, sizeof(Cell));
    } else {
        for(int i=0;i<S->gw*S->gh;i++) S->grid[i].count=0;
    }
    double R=S->prm.radius; double eps=1e-9;
    S->clamped_count = 0;   /* ##CHRIS */
    for(int i=0;i<S->prm.N;i++){
        /* ##CHRIS: snapshot so a clamp can be detected, bookkept and traced. */
        const double dbg_vx0 = S->P[i].vx, dbg_vy0 = S->P[i].vy;
        const double dbg_x0  = S->P[i].x,  dbg_y0  = S->P[i].y;
        /* safety clamp to ensure inside box */
        if (S->P[i].x < R) { S->P[i].x = R + eps; if (S->P[i].vx < 0) S->P[i].vx = -S->P[i].vx; }
        if (S->P[i].x > S->prm.boxW - R) { S->P[i].x = S->prm.boxW - R - eps; if (S->P[i].vx > 0) S->P[i].vx = -S->P[i].vx; }
        if (S->P[i].y < R) { S->P[i].y = R + eps; if (S->P[i].vy < 0) S->P[i].vy = -S->P[i].vy; }
        if (S->P[i].y > S->prm.boxH - R) { S->P[i].y = S->prm.boxH - R - eps; if (S->P[i].vy > 0) S->P[i].vy = -S->P[i].vy; }
        /* if inside divider slabs, push to nearest face and reflect in divider frame */
        if (S->prm.divider_count > 0) {
            int dcount = clamp_dividers(S->prm.divider_count);
            for (int d = 0; d < dcount; ++d) {
                if (S->prm.divider_thickness[d] <= 0.0) continue;
                double L = S->prm.divider_x[d] - 0.5*S->prm.divider_thickness[d];
                double Rf = S->prm.divider_x[d] + 0.5*S->prm.divider_thickness[d];
                if (L > 0.0 && Rf < S->prm.boxW) {
                    double gapL = (S->P[i].x - R) - L;
                    double gapR = Rf - (S->P[i].x + R);
                    if (gapL > -eps && gapR > -eps) {
                        if (gapL < gapR) { S->P[i].x = L - eps + R; S->P[i].vx = 2.0*S->prm.divider_vx[d] - S->P[i].vx; }
                        else { S->P[i].x = Rf + eps - R; S->P[i].vx = 2.0*S->prm.divider_vx[d] - S->P[i].vx; }
                    }
                }
            }
        }
        /* ##CHRIS: a clamp above is a real boundary bounce. Bookkeep it like one:
           bump coll_count so this particle's stale events are invalidated, and record
           the index so the caller reschedules it. Without this the particle keeps flying
           with an event set computed for its OLD velocity - the defect behind the
           particle-pair tunnelling seen in the speed-of-sound wall_hold phase. */
        /* Only a VELOCITY change needs bookkeeping: that is a genuine bounce, and it is what
           invalidates the particle's pending events. A pure sub-nanometre position nudge
           (the common case) leaves the event set valid to within 1e-9 and is left alone, so
           this fix does not perturb trajectories that were never broken. */
        if (S->P[i].vx != dbg_vx0 || S->P[i].vy != dbg_vy0) {
            if (edmd_trace_active(S) && (i == g_trace_a || i == g_trace_b)) {
                fprintf(stderr,
                    "[EDMD-TRACE] grid_build CLAMPED p%d at t=%.17g (coll_count %d -> %d, will reschedule)\n"
                    "             pos (%.17g,%.17g) -> (%.17g,%.17g)\n"
                    "             vel (%.17g,%.17g) -> (%.17g,%.17g)\n",
                    i, S->t, S->P[i].coll_count, S->P[i].coll_count + 1,
                    dbg_x0, dbg_y0, S->P[i].x, S->P[i].y,
                    dbg_vx0, dbg_vy0, S->P[i].vx, S->P[i].vy);
                fflush(stderr);
            }
            S->P[i].coll_count++;
            S->clamp_repair_count++;
            if (S->clamped && S->clamped_count < S->prm.N) {
                S->clamped[S->clamped_count++] = i;
            }
        }
        int cx = (int)floor(S->P[i].x / S->cell_size);
        int cy = (int)floor(S->P[i].y / S->cell_size);
        cx = clampi(cx, 0, S->gw-1);
        cy = clampi(cy, 0, S->gh-1);
        Cell* c = &S->grid[cy*S->gw + cx];
        cell_reserve(c, c->count+1);
        c->idx[c->count++] = i;
    }
}

/* ##CHRIS: reschedule everything grid_build() had to clamp. Must be called after any
   grid_build() that is NOT immediately followed by a full reschedule. Safe to call when
   nothing was clamped (the common case) - it is then a no-op. */
static void schedule_for(EDMD* S, int i);
static void reschedule_clamped(EDMD* S){
    const int n = S->clamped_count;
    if (n <= 0) return;
    S->clamped_count = 0;   /* clear first: schedule_for() does not clamp, but be safe */
    for (int k = 0; k < n; ++k) {
        schedule_for(S, S->clamped[k]);
    }
}

/* ------------------------ collision-time solvers ------------------------ */

static inline void harmonic_advance_1d(double k, double m, double x_eq, double dt, double* x_inout, double* v_inout){
    if (!(k > 0.0) || !(m > 0.0)) return;
    if (!(dt > 0.0)) return;
    const double w = sqrt(k / m);
    if (!(w > 0.0)) return;
    const double c = cos(w * dt);
    const double s = sin(w * dt);
    const double x0 = *x_inout;
    const double v0 = *v_inout;
    const double dx = x0 - x_eq;
    const double x1 = x_eq + dx * c + (v0 / w) * s;
    const double v1 = -dx * w * s + v0 * c;
    *x_inout = x1;
    *v_inout = v1;
}

/* disk–disk, equal radii R: solve ||r + t v|| = 2R with r·v < 0 (approaching).
   quadratic; take smallest positive root. */
static int collide_time_ab(double xi,double yi,double vxi,double vyi,
                           double xj,double yj,double vxj,double vyj,
                           double R, double* tcol)
{
    double rx=xj-xi, ry=yj-yi;
    double vx=vxj-vxi, vy=vyj-vyi;
    double b = rx*vx + ry*vy;              /* r·v */
    if(b>=0.0) return 0;                    /* separating */
    double rr=rx*rx+ry*ry, vv=vx*vx+vy*vy;
    double sig = 2.0*R;
    double c = rr - sig*sig;
    if(vv<=0.0) return 0;
    /* ##CHRIS: SAFETY NET. c < 0 means the pair is already overlapping, and b < 0 means it
       is still approaching - the collision is overdue, so schedule it immediately.
       Without this branch the earlier root is negative and the `t<=1e-12` guard below
       discards it, so an overlapping approaching pair becomes PERMANENTLY invisible to the
       scheduler: c only grows more negative until they separate on the far side. That is
       what turned a single missed event into two particles passing through each other.
       This should never fire once the scheduler is correct, so every occurrence is counted
       and reported by edmd_overlap_repair_count() rather than being silently absorbed. */
    if(c<0.0){ *tcol = 0.0; return 2; }
    double disc = b*b - vv*c;
    if(disc<=0.0) return 0;
    double t = (-b - sqrt(disc)) / vv;     /* earlier root */
    if(t<=1e-12) return 0;
    *tcol = t; return 1;
}

/* left wall x=0: hit when x - R = 0; particle must move left (vx<0). */
/* ##CHRIS: `gap` is the distance still available before the face is reached and `speed` is
   the (positive) closing speed. gap <= 0 means the particle is already at or past the face
   while still moving outward, i.e. the wall collision is OVERDUE - schedule it now (t=0)
   rather than discarding it. The old `t<=1e-12` guard rejected exactly this case, so a
   particle seeded precisely on a face (gap == 0) got no wall event at all, drifted out of
   the box, and was then silently repaired by grid_build()'s clamp. That was the first link
   in the particle-pair tunnelling bug (seed 2381038820).
   The guard still applies for gap > 0, and cannot mis-fire on a just-resolved collision:
   resolve_wall() reverses the velocity, so the velocity-sign test in each caller rejects it. */
static int wall_time_from_gap(double gap, double speed, double* tcol){
    if (gap <= 0.0) { *tcol = 0.0; return 2; }   /* overdue */
    const double t = gap / speed;
    if (t <= 1e-12) return 0;
    *tcol = t; return 1;
}
static int collide_time_wall_L(const EDMD* S, const EDMD_Particle* A, double* tcol){
    if(A->vx >= 0.0) return 0;
    return wall_time_from_gap((A->x - S->prm.radius) - 0.0, -A->vx, tcol);
}
/* right wall x=boxW: hit when x + R = boxW; particle must move right (vx>0). */
static int collide_time_wall_R(const EDMD* S, const EDMD_Particle* A, double* tcol){
    if(A->vx <= 0.0) return 0;
    return wall_time_from_gap((S->prm.boxW - S->prm.radius) - A->x, A->vx, tcol);
}
/* bottom wall y=0: hit when y - R = 0; particle must move down (vy<0). */
static int collide_time_wall_B(const EDMD* S, const EDMD_Particle* A, double* tcol){
    if(A->vy >= 0.0) return 0;
    return wall_time_from_gap((A->y - S->prm.radius) - 0.0, -A->vy, tcol);
}
/* top wall y=boxH: hit when y + R = boxH; particle must move up (vy>0). */
static int collide_time_wall_T(const EDMD* S, const EDMD_Particle* A, double* tcol){
    if(A->vy <= 0.0) return 0;
    return wall_time_from_gap((S->prm.boxH - S->prm.radius) - A->y, A->vy, tcol);
}

static double harmonic_contact_value(double t, double base, double particle_v,
                                     double wall_dx, double wall_v0,
                                     double omega){
    return base + particle_v * t
         - wall_dx * cos(omega * t)
         - (wall_v0 / omega) * sin(omega * t);
}

/* Find the first directed zero of
     f(t) = base + particle_v*t - wall_dx*cos(w*t) - wall_v0/w*sin(w*t).
   The derivative is a sinusoid plus a constant. Its analytically known zeros
   split time into monotonic intervals, so checking those interval endpoints
   cannot skip a fast wall crossing the way a fixed coarse scan can. */
static int harmonic_first_contact(double base, double particle_v,
                                  double wall_dx, double wall_v0,
                                  double omega, double t_max,
                                  int crossing_direction, double* tcol){
    if (!(omega > 0.0) || !(t_max > 0.0) || !tcol) return 0;
    const double two_pi = 2.0 * M_PI;
    const double period = two_pi / omega;
    const double derivative_amplitude = hypot(wall_dx * omega, wall_v0);
    const double root_tol = 1e-12;
    double previous_t = 0.0;
    double previous_f = harmonic_contact_value(
        previous_t, base, particle_v, wall_dx, wall_v0, omega);

    double offsets[2] = {NAN, NAN};
    int offset_count = 0;
    if (derivative_amplitude > 1e-15 &&
        fabs(particle_v) < derivative_amplitude) {
        const double q = fmax(-1.0, fmin(1.0,
            -particle_v / derivative_amplitude));
        const double phase = atan2(wall_v0, wall_dx * omega);
        const double alpha = asin(q);
        double theta_a = fmod(phase + alpha, two_pi);
        double theta_b = fmod(phase + M_PI - alpha, two_pi);
        if (theta_a < 0.0) theta_a += two_pi;
        if (theta_b < 0.0) theta_b += two_pi;
        offsets[0] = theta_a / omega;
        offsets[1] = theta_b / omega;
        if (offsets[1] < offsets[0]) {
            const double tmp = offsets[0]; offsets[0] = offsets[1]; offsets[1] = tmp;
        }
        offset_count = (fabs(offsets[1] - offsets[0]) <= 1e-14) ? 1 : 2;
    }

    const long max_cycles = (long)ceil(t_max / period) + 1L;
    for (long cycle = 0; cycle <= max_cycles; ++cycle) {
        for (int oi = 0; oi < offset_count; ++oi) {
            double current_t = (double)cycle * period + offsets[oi];
            if (current_t <= previous_t + 1e-14) continue;
            if (current_t > t_max) current_t = t_max;
            const double current_f = harmonic_contact_value(
                current_t, base, particle_v, wall_dx, wall_v0, omega);
            const int crosses = (crossing_direction > 0)
                ? (previous_f <= 0.0 && current_f >= 0.0)
                : (previous_f >= 0.0 && current_f <= 0.0);
            if (crosses) {
                double lo = previous_t, hi = current_t;
                for (int it = 0; it < 80; ++it) {
                    const double mid = 0.5 * (lo + hi);
                    const double fm = harmonic_contact_value(
                        mid, base, particle_v, wall_dx, wall_v0, omega);
                    if ((crossing_direction > 0 && fm >= 0.0) ||
                        (crossing_direction < 0 && fm <= 0.0)) hi = mid;
                    else lo = mid;
                }
                const double root = 0.5 * (lo + hi);
                if (root > root_tol) { *tcol = root; return 1; }
                /* f(0)==0 is normal immediately after a wall collision.  If
                   the particle initially separates, the harmonic wall can
                   still catch it later.  Skip only this zero-time root and
                   continue through the remaining monotonic intervals. */
            }
            previous_t = current_t;
            previous_f = current_f;
            if (current_t >= t_max - 1e-15) return 0;
        }
        if (offset_count == 0) break;
    }

    if (previous_t < t_max) {
        const double final_f = harmonic_contact_value(
            t_max, base, particle_v, wall_dx, wall_v0, omega);
        const int crosses = (crossing_direction > 0)
            ? (previous_f <= 0.0 && final_f >= 0.0)
            : (previous_f >= 0.0 && final_f <= 0.0);
        if (crosses) {
            double lo = previous_t, hi = t_max;
            for (int it = 0; it < 80; ++it) {
                const double mid = 0.5 * (lo + hi);
                const double fm = harmonic_contact_value(
                    mid, base, particle_v, wall_dx, wall_v0, omega);
                if ((crossing_direction > 0 && fm >= 0.0) ||
                    (crossing_direction < 0 && fm <= 0.0)) hi = mid;
                else lo = mid;
            }
            const double root = 0.5 * (lo + hi);
            if (root > root_tol) { *tcol = root; return 1; }
        }
    }
    return 0;
}

/* ------------------------ scheduling ------------------------ */

static void schedule_walls(EDMD* S, int i){
    /* ##CHRIS: rc==2 means the wall collision was already overdue when scheduled. The old
       code discarded exactly those, which is how a particle escaped the box. Counted so the
       exposure is measurable rather than invisible. */
    double t; int rc;
    if((rc = collide_time_wall_L(S, &S->P[i], &t))) {
        if(rc==2) S->wall_overdue_count++;
        heap_push(&S->heap, (Event){ S->t+t, i,-1, S->P[i].coll_count,0, EV_WL });
    }
    if((rc = collide_time_wall_R(S, &S->P[i], &t))) {
        if(rc==2) S->wall_overdue_count++;
        heap_push(&S->heap, (Event){ S->t+t, i,-1, S->P[i].coll_count,0, EV_WR });
    }
    if((rc = collide_time_wall_B(S, &S->P[i], &t))) {
        if(rc==2) S->wall_overdue_count++;
        heap_push(&S->heap, (Event){ S->t+t, i,-1, S->P[i].coll_count,0, EV_WB });
    }
    if((rc = collide_time_wall_T(S, &S->P[i], &t))) {
        if(rc==2) S->wall_overdue_count++;
        heap_push(&S->heap, (Event){ S->t+t, i,-1, S->P[i].coll_count,0, EV_WT });
    }
}

/* divider faces (vertical slab): faces at x = cx - th/2 and x = cx + th/2 */
/* For a vertical slab at [L,R], hits occur from the left side onto L when vx>0 and x+R < L,
   and from the right side onto R when vx<0 and x-R > R. */
static int collide_time_divider_L(const EDMD* S, const EDMD_Particle* A, int d, double* tcol){
    if(S->prm.divider_count <= 0) return 0;
    if(d < 0 || d >= clamp_dividers(S->prm.divider_count)) return 0;
    if(S->prm.divider_thickness[d] <= 0.0) return 0;
    double L = S->prm.divider_x[d] - 0.5*S->prm.divider_thickness[d];
    if(!(L > 0.0 && L < S->prm.boxW)) return 0;
    if (divider_has_spring(S, d)) {
        /* Find earliest t>=0 such that (x+R) = (xw(t) - th/2). Wall center xw(t) is harmonic. */
        const double th = S->prm.divider_thickness[d];
        const double x0 = S->prm.divider_x[d];
        const double v0 = S->prm.divider_vx[d];
        const double k  = S->prm.divider_k[d];
        const double mW = S->prm.divider_mass[d];
        const double xeq = S->prm.divider_xeq[d];
        const double Rpart = S->prm.radius;
        const double gap = (x0 - 0.5*th) - (A->x + Rpart);
        /* A zero gap is the post-collision contact state.  It must remain
           schedulable because an oscillating wall may catch the particle
           again after they initially separate.  Only reject a particle that
           is already materially through this face. */
        if (gap < -1e-9) return 0;

        const double w = sqrt(k / mW);
        if (!(w > 0.0)) return 0;
        const double T = 2.0 * M_PI / w;
        const double amp = hypot((x0 - xeq), (v0 / w));
        /* f(t) = x_p(t) - (xw(t) - th/2 - R) */
        const double base = (A->x - xeq) + (Rpart + 0.5*th);
        const double v = A->vx;

        /* If v<=0, any crossing must happen early (wall moving left into particle). */
        double t_max;
        if (v > 0.0) {
            t_max = (fabs(base) + amp + 1.0) / v + 2.0 * T;
        } else {
            if (base + amp < 0.0) return 0;
            t_max = fmax(2.0 * T, 0.5 * T);
        }
        if (!(t_max > 0.0)) return 0;

        return harmonic_first_contact(base, v, x0 - xeq, v0, w,
                                      t_max, +1, tcol);
    }
    /* left face, approached from left side: (x+R) + vx t = L + vW t */
    double rel = A->vx - S->prm.divider_vx[d];
    double num = L - (A->x + S->prm.radius);
    if(num <= 1e-12) return 0;      /* not strictly left of face */
    if(rel <= 0.0) return 0;        /* not approaching in wall frame */
    double t = num / rel;
    if(t<=1e-12) return 0; *tcol=t; return 1;
}
static int collide_time_divider_R(const EDMD* S, const EDMD_Particle* A, int d, double* tcol){
    if(S->prm.divider_count <= 0) return 0;
    if(d < 0 || d >= clamp_dividers(S->prm.divider_count)) return 0;
    if(S->prm.divider_thickness[d] <= 0.0) return 0;
    double Rf = S->prm.divider_x[d] + 0.5*S->prm.divider_thickness[d];
    if(!(Rf > 0.0 && Rf < S->prm.boxW)) return 0;
    if (divider_has_spring(S, d)) {
        /* Earliest t>=0 such that (x-R) = (xw(t) + th/2). */
        const double th = S->prm.divider_thickness[d];
        const double x0 = S->prm.divider_x[d];
        const double v0 = S->prm.divider_vx[d];
        const double k  = S->prm.divider_k[d];
        const double mW = S->prm.divider_mass[d];
        const double xeq = S->prm.divider_xeq[d];
        const double Rpart = S->prm.radius;
        const double gap = (A->x - Rpart) - (x0 + 0.5*th);
        if (gap < -1e-9) return 0;

        const double w = sqrt(k / mW);
        if (!(w > 0.0)) return 0;
        const double T = 2.0 * M_PI / w;
        const double amp = hypot((x0 - xeq), (v0 / w));
        /* f(t) = x_p(t) - (xw(t) + th/2 + R) */
        const double base = (A->x - xeq) - (Rpart + 0.5*th);
        const double v = A->vx;

        double t_max;
        if (v < 0.0) {
            t_max = (fabs(base) + amp + 1.0) / (-v) + 2.0 * T;
        } else {
            if (base - amp > 0.0) return 0;
            t_max = fmax(2.0 * T, 0.5 * T);
        }
        if (!(t_max > 0.0)) return 0;

        return harmonic_first_contact(base, v, x0 - xeq, v0, w,
                                      t_max, -1, tcol);
    }
    /* right face, approached from right side: (x-R) + vx t = Rf + vW t */
    double rel = A->vx - S->prm.divider_vx[d];
    double dist = (A->x - S->prm.radius) - Rf; /* >0 if strictly to the right */
    if(dist <= 1e-12) return 0;      /* not strictly right of face */
    if(rel >= 0.0) return 0;         /* must move left relative to wall */
    double t = dist / (-rel);        /* since (vx-vW) t = -dist */
    if(t<=1e-12) return 0; *tcol=t; return 1;
}

static void schedule_divider(EDMD* S, int i){
    if(S->prm.divider_count <= 0) return;
    int dcount = clamp_dividers(S->prm.divider_count);
    for (int d = 0; d < dcount; ++d) {
        double t;
        if(collide_time_divider_L(S, &S->P[i], d, &t))
            heap_push(&S->heap, (Event){ S->t+t, i, d, S->P[i].coll_count,0, EV_DL });
        if(collide_time_divider_R(S, &S->P[i], d, &t))
            heap_push(&S->heap, (Event){ S->t+t, i, d, S->P[i].coll_count,0, EV_DR });
    }
}

/* Pistons: left piston right face at xL, right piston left face at xR */
static int collide_time_piston_L(const EDMD* S, const EDMD_Particle* A, double* tcol){
    if(!S->prm.has_pistonL) return 0;
    double xL = S->prm.pistonL_x;
    /* (x - R) + vx t = xL + vW t  => t = (xL - (x-R)) / (vx - vW), need rel<0 and (x-R) > xL (numerator<0) */
    double rel = A->vx - S->prm.pistonL_vx;
    if(rel >= 0.0) return 0;
    double num = xL - (A->x - S->prm.radius);
    if(num >= -1e-12) return 0; /* require particle to the right of face */
    double t = num / rel; if(t<=1e-12) return 0; *tcol=t; return 1;
}
static int collide_time_piston_R(const EDMD* S, const EDMD_Particle* A, double* tcol){
    if(!S->prm.has_pistonR) return 0;
    double xR = S->prm.pistonR_x;
    /* (x + R) + vx t = xR + vW t  => t = (xR - (x+R)) / (vx - vW), need rel>0 and numerator>0 */
    double rel = A->vx - S->prm.pistonR_vx;
    if(rel <= 0.0) return 0;
    double num = xR - (A->x + S->prm.radius);
    if(num <= 1e-12) return 0;
    double t = num / rel; if(t<=1e-12) return 0; *tcol=t; return 1;
}
static void schedule_pistons(EDMD* S, int i){
    double t;
    if (collide_time_piston_L(S, &S->P[i], &t))
        heap_push(&S->heap, (Event){ S->t+t, i,-1, S->P[i].coll_count,0, EV_PL });
    if (collide_time_piston_R(S, &S->P[i], &t))
        heap_push(&S->heap, (Event){ S->t+t, i,-1, S->P[i].coll_count,0, EV_PR });
}

static void schedule_ab(EDMD* S, int i, int j){
    if (!S->prm.pp_collisions_enabled) return; /* ##CHRIS: skip PP if disabled */
    double t;
    const int ok = collide_time_ab(S->P[i].x,S->P[i].y,S->P[i].vx,S->P[i].vy,
                                   S->P[j].x,S->P[j].y,S->P[j].vx,S->P[j].vy,
                                   S->prm.radius, &t);
    /* ##CHRIS: trace every scheduling decision for the watched pair. */
    if (edmd_trace_active(S) &&
        ((i == g_trace_a && j == g_trace_b) || (i == g_trace_b && j == g_trace_a))) {
        edmd_trace_schedule_attempt(S, i, j, ok, ok ? S->t + t : NAN);
    }
    if(ok){
        if(ok == 2) S->overlap_repair_count++;   /* ##CHRIS: overdue (already-overlapping) pair */
        heap_push(&S->heap, (Event){ S->t+t, i,j, S->P[i].coll_count,S->P[j].coll_count, EV_AB });
    }
}

static void advance_dividers(EDMD* S, double dt){
    if (S->prm.divider_count <= 0) return;
    int dcount = clamp_dividers(S->prm.divider_count);
    for (int d = 0; d < dcount; ++d) {
        const double th = (S->prm.divider_thickness[d] > 0.0) ? S->prm.divider_thickness[d] : 0.0;
        const double left_limit  = 0.5 * th + 1e-9;
        const double right_limit = S->prm.boxW - 0.5 * th - 1e-9;
        if (divider_has_spring(S, d)) {
            harmonic_advance_1d(S->prm.divider_k[d], S->prm.divider_mass[d], S->prm.divider_xeq[d], dt,
                                &S->prm.divider_x[d], &S->prm.divider_vx[d]);
            /* Safety clamp: prevent the slab from leaving the box due to numerical/root errors. */
            if (S->prm.divider_x[d] < left_limit) { S->prm.divider_x[d] = left_limit; if (S->prm.divider_vx[d] < 0.0) S->prm.divider_vx[d] = 0.0; }
            if (S->prm.divider_x[d] > right_limit) { S->prm.divider_x[d] = right_limit; if (S->prm.divider_vx[d] > 0.0) S->prm.divider_vx[d] = 0.0; }
        } else {
            S->prm.divider_x[d] += S->prm.divider_vx[d] * dt;
            if (S->prm.divider_x[d] < left_limit) { S->prm.divider_x[d] = left_limit; }
            if (S->prm.divider_x[d] > right_limit) { S->prm.divider_x[d] = right_limit; }
        }
    }
}

/* schedule for one particle i: walls + neighbors from 9-cell stencil */
static void schedule_for(EDMD* S, int i){
    /* Robust: schedule walls and AB with all others. For N up to a few thousands this is fine. */
    schedule_walls(S, i);
    schedule_divider(S, i);
    schedule_pistons(S, i);
    for (int j = 0; j < S->prm.N; ++j) {
        if (j == i) continue;
        schedule_ab(S, i, j);
    }
}

/* schedule only AB partners for particle i (PP-only mode) */
static void schedule_for_pponly(EDMD* S, int i){
    for (int j = 0; j < S->prm.N; ++j) {
        if (j == i) continue;
        schedule_ab(S, i, j);
    }
}

/* clear heap, rebuild grid, schedule all */
static void reschedule_all_internal(EDMD* S){
    S->heap.n = 0;
    grid_build(S);
    /* schedule walls for all */
    for (int i = 0; i < S->prm.N; ++i){ schedule_walls(S, i); schedule_divider(S,i); schedule_pistons(S, i);} 
    /* schedule AB for all pairs (i<j) */
    for (int i = 0; i < S->prm.N; ++i) {
        for (int j = i+1; j < S->prm.N; ++j) {
            schedule_ab(S, i, j);
        }
    }
    /* ##CHRIS: every particle was just scheduled, so anything grid_build clamped is
       already covered - drop the pending list rather than scheduling it twice. */
    S->clamped_count = 0;
}

/* reschedule AB only (no walls/pistons/divider; do not clamp/push) */
static void reschedule_all_internal_pponly(EDMD* S){
    S->heap.n = 0;
    /* schedule AB for all pairs (i<j) */
    for (int i = 0; i < S->prm.N; ++i) {
        for (int j = i+1; j < S->prm.N; ++j) {
            schedule_ab(S, i, j);
        }
    }
}

/* ------------------------ resolution ------------------------ */

/* equal-mass elastic collision for AB: exchange normal component */
static void resolve_ab(EDMD* S, int i, int j){
    (void)S;
    EDMD_Particle *A=&S->P[i], *B=&S->P[j];
    double dx=B->x-A->x, dy=B->y-A->y;
    double dist = sqrt(dx*dx+dy*dy);
    if(dist<=0.0){ dx = S->prm.radius; dy=0.0; dist=S->prm.radius; }
    double nx=dx/dist, ny=dy/dist;

    double dvx=B->vx - A->vx, dvy=B->vy - A->vy;
    double dvn=dvx*nx + dvy*ny;

    /* ##CHRIS: virial BEFORE the velocities are overwritten. dvn < 0 for an
       approaching pair, so -dvn is the closing speed; unit mass, and dist is the
       contact separation (= sigma up to the scheduler's tolerance). */
    S->virial_accum += (-dvn) * dist;
    S->virial_pair_events++;

    A->vx += dvn*nx; A->vy += dvn*ny;
    B->vx -= dvn*nx; B->vy -= dvn*ny;

    A->coll_count++; B->coll_count++;
}

//============================================================================
// ##CHRIS: BEGIN HEAT BATH FUNCTIONS FOR EDMD MODE
//============================================================================

/* Sample from Gaussian distribution (Box-Muller) */
static inline double sample_gaussian_edmd(double mean, double sigma) {
    double u1, u2, z;
    do {
        u1 = (double)rand() / (double)RAND_MAX;
        u2 = (double)rand() / (double)RAND_MAX;
        if (u1 < 1e-12) u1 = 1e-12; /* avoid log(0) */
        z = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
    } while (isnan(z));
    return mean + sigma * z;
}

/* Sample speed from 2D Maxwell-Boltzmann distribution */
static inline double sample_MB_speed_2D(double temperature, double mass, double kB) {
    double rand1 = (double)rand() / (double)RAND_MAX;
    if (rand1 < 1e-12) rand1 = 1e-12; /* avoid log(0) */
    double speed = sqrt(-2.0 * kB * temperature / mass * log(rand1));
    return speed;
}

/* Adaptive thermal wall bounce for EDMD */
/* ##CHRIS: Base thermal wall bounce - Mode 1 (standard MB, no overshoot) */
static inline void thermal_wall_bounce_edmd(EDMD_Particle* A,
                                              double normal_x, double normal_y,
                                              const EDMD_Params* prm) {
    /* Sample speed from MB distribution at heat bath temperature */
    double speed = sample_MB_speed_2D(prm->heatbath_temperature, prm->particle_mass, prm->kB);

    /* Random angle in hemisphere (±90° from normal) */
    double rand_val = (double)rand() / (double)RAND_MAX;
    double angle = (rand_val - 0.5) * M_PI;

    /* Rotate normal by random angle */
    double cos_a = cos(angle);
    double sin_a = sin(angle);

    A->vx = speed * (normal_x * cos_a - normal_y * sin_a);
    A->vy = speed * (normal_x * sin_a + normal_y * cos_a);
}

/* ##CHRIS: Gradual damping bounce - Mode 0 (friction-like) */
static inline void gradual_damping_bounce_edmd(EDMD_Particle* A,
                                                 double gas_temp,
                                                 const EDMD_Params* prm) {
    double dV = 0.01;
    double dV_sign = (prm->heatbath_temperature - gas_temp) / fabs(prm->heatbath_temperature - gas_temp);

    /* Add small random perturbations */
    double sigma = 0.1 * sqrt(A->vx * A->vx + A->vy * A->vy);
    A->vx += sample_gaussian_edmd(0.0, sigma);
    A->vy += sample_gaussian_edmd(0.0, sigma);

    /* Gradually adjust velocity */
    if (fabs(A->vx) > fabs(A->vy)) {
        A->vx += dV * dV_sign;
    } else {
        A->vy += dV * dV_sign;
    }
}

/* ##CHRIS: Adaptive thermal wall bounce - Mode 2 (overshoots for faster convergence) */
static inline void adaptive_thermal_wall_bounce_edmd(EDMD_Particle* A,
                                                       double normal_x, double normal_y,
                                                       double gas_temp, const EDMD_Params* prm) {
    double sample_temp;
    double temp_error = fabs(gas_temp - prm->heatbath_temperature) / prm->heatbath_temperature;

    if (temp_error < prm->stability_window_percent) {
        sample_temp = prm->heatbath_temperature; /* Within stability window */
    } else {
        if (gas_temp < prm->heatbath_temperature) {
            sample_temp = prm->heatbath_temperature * prm->mb_overshoot_factor; /* Too cold */
        } else {
            sample_temp = prm->heatbath_temperature / prm->mb_overshoot_factor; /* Too hot */
        }
    }

    /* Sample speed from MB distribution */
    double speed = sample_MB_speed_2D(sample_temp, prm->particle_mass, prm->kB);

    /* Random angle in hemisphere (±90° from normal) */
    double rand_val = (double)rand() / (double)RAND_MAX;
    double angle = (rand_val - 0.5) * M_PI;

    /* Rotate normal by random angle */
    double cos_a = cos(angle);
    double sin_a = sin(angle);

    A->vx = speed * (normal_x * cos_a - normal_y * sin_a);
    A->vy = speed * (normal_y * cos_a + normal_x * sin_a);
}

//============================================================================
// ##CHRIS: END HEAT BATH FUNCTIONS
//============================================================================

/* reflect on static walls */
/* ##CHRIS: MODIFIED - Added heat bath support for outer walls (WL, WR, WB, WT) */
static double divider_species_side_temperature(const EDMD* S, int d, int species, int side){
    if (!S) return 0.0;
    if (d < 0 || d >= clamp_dividers(S->prm.divider_count)) return 0.0;
    const double cx = S->prm.divider_x[d];
    double ke = 0.0;
    int count = 0;
    for (int i = 0; i < S->prm.N; ++i) {
        const int sp = (S->prm.species ? S->prm.species[i] : 0);
        if (sp != species) continue;
        const int pside = (S->P[i].x < cx) ? 0 : 1;
        if (pside != side) continue;
        ke += 0.5 * S->prm.particle_mass * (S->P[i].vx * S->P[i].vx + S->P[i].vy * S->P[i].vy);
        count++;
    }
    if (count <= 0) return 0.0;
    if (!(S->prm.kB > 0.0)) return 0.0;
    return ke / ((double)count * S->prm.kB);
}

static inline int divider_gate_allows_pass(const EDMD* S, int d, EvType type, int species, int i){
    if (!S) return 0;
    if (d < 0 || d >= EDMD_MAX_DIVIDERS) return 0;
    const int mode = S->prm.divider_gate_mode[d];
    if (mode == 0) return 0;
    const double speed = (i >= 0 && i < S->prm.N)
        ? sqrt(S->P[i].vx * S->P[i].vx + S->P[i].vy * S->P[i].vy)
        : 0.0;
    const double speed_ref = (S->prm.divider_gate_speed_ref[d] > 0.0) ? S->prm.divider_gate_speed_ref[d] : 0.0;
    const double speed_mult = (S->prm.divider_gate_hot_ratio[d] > 0.0) ? S->prm.divider_gate_hot_ratio[d] : 0.0;
    const int hot_enough = (speed_ref > 0.0 && speed_mult > 0.0 && speed >= speed_mult * speed_ref);

    int target_side = -1; /* 0=LEFT, 1=RIGHT */
    if (mode == 1) {
        target_side = (S->prm.divider_gate_target_side[d] != 0) ? 1 : 0;
        if (!S->prm.species) return 0;
        if (species != S->prm.divider_gate_species[d]) {
            const int moving_toward_target = (target_side == 0) ? (type == EV_DR) : (type == EV_DL);
            return (moving_toward_target && hot_enough) ? 1 : 0;
        }
    } else if (mode == 2) {
        if (species == 0) target_side = 0;
        else if (species == 1) target_side = 1;
        else return 0;
    } else if (mode == 3) {
        /* Blocked-species one-way block gate: all other species pass freely. */
        const int blocked = S->prm.divider_gate_species[d];
        if (species != blocked) return 1; /* transparent */
        target_side = (S->prm.divider_gate_target_side[d] != 0) ? 1 : 0;
        {
            const int moving_toward_target = (target_side == 0) ? (type == EV_DR) : (type == EV_DL);
            if (!moving_toward_target && hot_enough) return 1;
        }
    } else if (mode == 4) {
        const int allowed = S->prm.divider_gate_species[d];
        target_side = (S->prm.divider_gate_target_side[d] != 0) ? 1 : 0;
        if (species != allowed) {
            const int moving_toward_target = (target_side == 0) ? (type == EV_DR) : (type == EV_DL);
            if (!moving_toward_target) return 0;
            const int source_side = (target_side == 0) ? 1 : 0;
            const int target_side_idx = target_side;
            const double t_hot = divider_species_side_temperature(S, d, species, source_side);
            const double t_ref = divider_species_side_temperature(S, d, allowed, target_side_idx);
            const double ratio = (S->prm.divider_gate_hot_ratio[d] > 0.0) ? S->prm.divider_gate_hot_ratio[d] : 0.0;
            if (!(ratio > 0.0)) return 0;
            if (t_ref <= 1e-12) return (t_hot > 0.0);
            return (t_hot >= ratio * t_ref);
        }
    } else {
        return 0;
    }

    /* Allow pass-through only in the direction *toward* the target side. */
    if (target_side == 0) return (type == EV_DR); /* right->left */
    return (type == EV_DL); /* left->right */
}

/* ##CHRIS: Paper 2 Level 0 -- gated per-event log of piston and divider
   collisions. PRINT ONLY: nothing here reads back into the dynamics, and every
   value is captured from the same locals the collision rule already computed.
   Enabled by the driver via edmd_set_event_log(); disabled (NULL) by default, so
   no existing run changes. time_scale converts the core's internal time to
   sigma-time (the driver passes PIXELS_PER_SIGMA); the core itself stays free of
   driver constants.
   Columns: t_sigma, kind (PL/PR/D<idx>), u_wall, v_before, v_after, dE, dp.
   dE is the quantity the core books as work (particle KE change for a
   prescribed-velocity wall; the WALL's KE change for a finite-mass wall -- these
   are different quantities, see the note in the Level-0 document). dp is always
   the particle's momentum change m(v_after - v_before). */
static FILE*  g_edmd_evlog = NULL;
static double g_edmd_evlog_tscale = 1.0;

void edmd_set_event_log(const char* path, double time_scale){
    if (g_edmd_evlog) { fclose(g_edmd_evlog); g_edmd_evlog = NULL; }
    g_edmd_evlog_tscale = (time_scale > 0.0) ? time_scale : 1.0;
    if (!path || !*path) return;
    g_edmd_evlog = fopen(path, "w");
    if (g_edmd_evlog) {
        fprintf(g_edmd_evlog, "t_sigma,kind,u_wall,v_before,v_after,dE,dp\n");
    }
}
void edmd_close_event_log(void){
    if (g_edmd_evlog) { fflush(g_edmd_evlog); fclose(g_edmd_evlog); g_edmd_evlog = NULL; }
}
static inline void edmd_log_event(const EDMD* S, const char* kind,
                                  double u, double v0, double v1, double dE){
    if (!g_edmd_evlog) return;
    fprintf(g_edmd_evlog, "%.12g,%s,%.12g,%.12g,%.12g,%.12g,%.12g\n",
            S->t / g_edmd_evlog_tscale, kind, u, v0, v1, dE, 1.0 * (v1 - v0));
}

static void resolve_wall(EDMD* S, int i, EvType type, int b){
    EDMD_Particle* A=&S->P[i];

    /* ##CHRIS: Handle outer wall collisions with optional heat bath */
    if(type==EV_WL || type==EV_WR || type==EV_WB || type==EV_WT){
        if(S->prm.heatbath_enabled){
            /* Get current gas temperature */
            double gas_temp = edmd_gas_temperature(S, S->prm.particle_mass, S->prm.kB);
            double temp_diff = fabs(gas_temp - S->prm.heatbath_temperature);

            if(temp_diff > 0.0001 * S->prm.heatbath_temperature){
                const double E0 = 0.5 * (A->vx*A->vx + A->vy*A->vy);
                /* Heat bath active - select mode based on thermal_wall_mode */
                double normal_x = 0.0, normal_y = 0.0;

                /* Determine wall normal */
                if(type==EV_WL) { normal_x = 1.0; normal_y = 0.0; }       /* Left wall, normal right */
                else if(type==EV_WR) { normal_x = -1.0; normal_y = 0.0; } /* Right wall, normal left */
                else if(type==EV_WB) { normal_x = 0.0; normal_y = 1.0; }  /* Bottom wall, normal up */
                else if(type==EV_WT) { normal_x = 0.0; normal_y = -1.0; } /* Top wall, normal down */

                /* Apply thermal wall mode */
                if(S->prm.thermal_wall_mode == 1) {
                    /* Mode 1: Base Maxwell-Boltzmann (most realistic) */
                    thermal_wall_bounce_edmd(A, normal_x, normal_y, &S->prm);
                } else if(S->prm.thermal_wall_mode == 2) {
                    /* Mode 2: Adaptive MB (fastest, overshoots) */
                    adaptive_thermal_wall_bounce_edmd(A, normal_x, normal_y, gas_temp, &S->prm);
                } else {
                    /* Mode 0: Gradual damping (friction-like) */
                    gradual_damping_bounce_edmd(A, gas_temp, &S->prm);
                }

                const double E1 = 0.5 * (A->vx*A->vx + A->vy*A->vy);
                S->heat_bath += (E1 - E0);
                A->coll_count++;
                return;
            }
        }
        /* No heat bath or at equilibrium - normal elastic reflection */
        /* ##CHRIS: record the ACTUAL momentum change rather than assuming 2|v_n|,
           so the estimator stays correct if the reflection rule ever changes. */
        {
            const double vx0 = A->vx, vy0 = A->vy;
            if(type==EV_WL || type==EV_WR){ A->vx = -A->vx; }
            if(type==EV_WB || type==EV_WT){ A->vy = -A->vy; }
            const double m = (S->prm.particle_mass > 0.0) ? S->prm.particle_mass : 1.0;
            const int w = (type==EV_WL) ? 0 : (type==EV_WR) ? 1 : (type==EV_WB) ? 2 : 3;
            const double dp = (w < 2) ? fabs(A->vx - vx0) : fabs(A->vy - vy0);
            S->wall_impulse[w] += m * dp;
            S->wall_events[w]++;
            /* ##CHRIS: outer-wall events in the gated log, so the momentum balance
               can close. The wall is stationary (u = 0) and the reflection is
               elastic (dE = 0); v_before/v_after are the NORMAL component, so
               WL/WR carry x-momentum and WB/WT carry y-momentum. */
            {
                static const char* wn[4] = {"WL","WR","WB","WT"};
                const double n0 = (w < 2) ? vx0 : vy0;
                const double n1 = (w < 2) ? A->vx : A->vy;
                edmd_log_event(S, wn[w], 0.0, n0, n1, 0.0);
            }
        }
        A->coll_count++;
        return;
    }
    if(type==EV_DL || type==EV_DR){
        int d = b;
        int dcount = clamp_dividers(S->prm.divider_count);
        if (d < 0 || d >= dcount) { A->coll_count++; return; }

        /* Szilard-style semipermeable divider: pass-through (no impulse) for allowed direction/species. */
        const int sp = (S->prm.species ? S->prm.species[i] : 0);
        if (divider_gate_allows_pass(S, d, type, sp, i)) {
            const double th = S->prm.divider_thickness[d];
            const double cx = S->prm.divider_x[d];
            const double R  = S->prm.radius;
            const double eps = 1e-9;

            double x_min = R;
            double x_max = S->prm.boxW - R;
            if (S->prm.has_pistonL) x_min = fmax(x_min, S->prm.pistonL_x + R + eps);
            if (S->prm.has_pistonR) x_max = fmin(x_max, S->prm.pistonR_x - R - eps);
            if (x_max < x_min) { double mid = 0.5 * (x_min + x_max); x_min = mid; x_max = mid; }

            if (type == EV_DL) {
                /* left->right: place just outside the right face */
                double xf = cx + 0.5 * th + R + eps;
                if (xf < x_min) xf = x_min;
                if (xf > x_max) xf = x_max;
                A->x = xf;
            } else {
                /* right->left: place just outside the left face */
                double xf = cx - 0.5 * th - R - eps;
                if (xf < x_min) xf = x_min;
                if (xf > x_max) xf = x_max;
                A->x = xf;
            }
            A->coll_count++;
            return;
        }

        double m = 1.0;
        double M = (S->prm.divider_mass[d] > 0.0 ? S->prm.divider_mass[d] : 0.0);
        double u1 = A->vx;
        double u2 = S->prm.divider_vx[d];
        if (M <= 0.0){
            /* Infinite-mass divider with prescribed velocity: track work via particle KE change */
            double v1 = 2.0 * u2 - u1;
            double dE = 0.5 * m * (v1*v1 - u1*u1);
            S->work_divider[d] += dE;
            A->vx = v1;
            { char kb[8]; snprintf(kb, sizeof kb, "D%d", d); edmd_log_event(S, kb, u2, u1, v1, dE); }
            A->coll_count++;
            return;
        } else {
            double v1 = ((m - M)*u1 + 2.0*M*u2) / (m + M);
            double v2 = ((M - m)*u2 + 2.0*m*u1) / (m + M);
            /* Work: change in divider KE */
            double dE = 0.5 * M * (v2*v2 - u2*u2);
            S->work_divider[d] += dE;
            { char kb[8]; snprintf(kb, sizeof kb, "D%d", d); edmd_log_event(S, kb, u2, u1, v1, dE); }
            A->vx = v1; S->prm.divider_vx[d] = v2; A->coll_count++; return;
        }
    }
    if(type==EV_PL || type==EV_PR){
        double m = 1.0;
        double M = (type==EV_PL ? S->prm.pistonL_mass : S->prm.pistonR_mass);
        double u1 = A->vx;
        double u2 = (type==EV_PL ? S->prm.pistonL_vx : S->prm.pistonR_vx);
        if (M <= 0.0){
            /* Infinite-mass piston with prescribed velocity: track work via particle KE change */
            double v1 = 2.0 * u2 - u1;
            double dE = 0.5 * m * (v1*v1 - u1*u1);
            if (type==EV_PL) { S->work_pistonL += dE; }
            else             { S->work_pistonR += dE; }
            A->vx = v1;
            edmd_log_event(S, (type==EV_PL) ? "PL" : "PR", u2, u1, v1, dE);
            A->coll_count++;
            return;
        }
        double v1 = ((m - M)*u1 + 2.0*M*u2) / (m + M);
        double v2 = ((M - m)*u2 + 2.0*m*u1) / (m + M);
        A->vx = v1;
        double dE = 0.5 * M * (v2*v2 - u2*u2);
        if (type==EV_PL) { S->prm.pistonL_vx = v2; S->work_pistonL += dE; }
        else             { S->prm.pistonR_vx = v2; S->work_pistonR += dE; }
        edmd_log_event(S, (type==EV_PL) ? "PL" : "PR", u2, u1, v1, dE);
        A->coll_count++; return;
    }
}

/* ------------------------ public API ------------------------ */

EDMD* edmd_create(const EDMD_Params* prm_in){
    EDMD* S = (EDMD*)calloc(1, sizeof(EDMD));
    S->prm = *prm_in;
    S->prm.divider_count = clamp_dividers(S->prm.divider_count);
    S->t = 0.0;
    for (int d = 0; d < EDMD_MAX_DIVIDERS; ++d) S->work_divider[d] = 0.0;
    S->work_pistonL = 0.0;
    S->work_pistonR = 0.0;
    S->heat_bath = 0.0;
    S->P = (EDMD_Particle*)calloc((size_t)S->prm.N, sizeof(EDMD_Particle));
    /* ##CHRIS: scratch list of particles clamped by grid_build(), so they can be rescheduled */
    S->clamped = (int*)calloc((size_t)(S->prm.N > 0 ? S->prm.N : 1), sizeof(int));
    S->clamped_count = 0;
    S->clamp_repair_count = 0;
    S->overlap_repair_count = 0;
    S->virial_accum = 0.0;
    S->virial_t0 = 0.0;
    S->virial_pair_events = 0;
    for(int w=0;w<4;w++){ S->wall_impulse[w]=0.0; S->wall_events[w]=0; }
    S->wall_thermal_events = 0;
    S->wall_overdue_count = 0;
    S->cell_size = (S->prm.cell_size>0.0)? S->prm.cell_size : fmax(2.5*S->prm.radius, 1.0*S->prm.radius);
    S->gw = (int)fmax(1.0, floor(S->prm.boxW / S->cell_size));
    S->gh = (int)fmax(1.0, floor(S->prm.boxH / S->cell_size));

    int heap_cap = (S->prm.max_events_hint>0)? S->prm.max_events_hint : (S->prm.N*16);
    heap_init(&S->heap, heap_cap);

    /* allocate grid arrays and build empty structure */
    reschedule_all_internal(S);
    return S;
}

void edmd_destroy(EDMD* S){
    if(!S) return;
    heap_free(&S->heap);
    grid_free(S);
    free(S->P);
    free(S->clamped);   /* ##CHRIS */
    free(S);
}

/* tiny xorshift64 rng */
static unsigned long long xorshift64(unsigned long long* s){
    unsigned long long x = *s;
    x ^= x<<13; x ^= x>>7; x ^= x<<17;
    *s = x; return x;
}

/* random non-overlapping init + Gaussian velocities (Box-Muller) */
/* ##CHRIS: lattice seeding for densities where random insertion cannot finish.
   edmd_init_random_gas() is rejection sampling with an unbounded retry loop (its
   own comment says "okay for demos"): above eta ~ 0.55 the acceptance
   probability for the last disks collapses and it never returns. That is a hang,
   not a slowdown, and it silently stalled an entire dense campaign.
   This follows the approach already settled in
   01_improvements_bugsfxed_dev/26_08_21_HIGH_ETA_INITIALIZATION_HEX_FALLBACK.md:
   place a CHECKED lattice, and never relax the overlap tolerance to make an
   infeasible layout pass. Rectangular is used when its spacing genuinely clears a
   diameter; otherwise hexagonal with
       dx = d(1+eps),  dy = (sqrt(3)/2) dx
   Velocities are drawn the same way as the random seeder, so only the positions
   differ. Returns 0 if even the hex lattice cannot fit N disks.
   edmd_init_random_gas() is deliberately left untouched: changing it would move
   every existing consumer's seeds. */
int edmd_init_lattice_gas(EDMD* S, unsigned long long seed){
    if(!S) return 0;
    unsigned long long st = seed ? seed : 0xC0FFEEULL;
    const int N = S->prm.N;
    const double R = S->prm.radius, d = 2.0*R;
    const double eps = 1e-3;                    /* separation margin */
    /* Inset from the walls. Placing a disk at exactly x = R leaves gap = 0 with
       the wall face, which wall_time_from_gap() correctly classifies as an
       overdue wall collision (the 2026-08-23 tunnelling fix). It is handled
       properly, but it is a marginal initial state and it trips the strict
       zero-health acceptance rule, so keep the seed strictly off the walls. */
    const double wmargin = 1e-3 * d;
    const double lo = R + wmargin;
    const double usableW = S->prm.boxW - d - 2.0*wmargin;
    const double usableH = S->prm.boxH - d - 2.0*wmargin;
    if(usableW <= 0.0 || usableH <= 0.0 || N <= 0) return 0;

    int placed = 0;

    /* 1) Rectangular lattice SPREAD OVER THE WHOLE BOX. The grid is sized to hold
          about N sites at the box aspect ratio, not packed as tightly as possible
          -- filling a tight lattice row-major from a corner would leave most of
          the box empty and the occupied part at near-close-packing, which is a
          completely different (and far denser) system than the requested eta. */
    {
        int cols = (int)ceil(sqrt((double)N * usableW / (usableH > 0 ? usableH : 1.0)));
        if(cols < 1) cols = 1;
        int rows = (int)ceil((double)N / cols);
        if(rows < 1) rows = 1;
        const double sx = (cols>1) ? usableW/(cols-1) : usableW;
        const double sy = (rows>1) ? usableH/(rows-1) : usableH;
        if(sx >= d*(1.0+eps) && sy >= d*(1.0+eps)){
            for(int r=0; r<rows && placed<N; ++r)
                for(int c=0; c<cols && placed<N; ++c){
                    S->P[placed].x = lo + ((cols>1) ? c*sx : 0.5*usableW);
                    S->P[placed].y = lo + ((rows>1) ? r*sy : 0.5*usableH);
                    placed++;
                }
        }
    }

    /* 2) Hexagonal fallback, for densities where no spread rectangular lattice
          clears a diameter. Here filling the box IS correct: at these packing
          fractions the disks genuinely occupy the whole area. */
    if(placed < N){
        placed = 0;
        const double dx = d*(1.0+eps);
        const double dy = 0.8660254037844386*dx;      /* sqrt(3)/2 */
        const int rows = (int)floor(usableH/dy) + 1;
        for(int r=0; r<rows && placed<N; ++r){
            const double yoff = (r & 1) ? 0.5*dx : 0.0;
            const double avail = usableW - yoff;
            if(avail < 0.0) continue;
            const int cols = (int)floor(avail/dx) + 1;
            for(int c=0; c<cols && placed<N; ++c){
                S->P[placed].x = lo + yoff + c*dx;
                S->P[placed].y = lo + r*dy;
                placed++;
            }
        }
        if(placed < N) return 0;                       /* genuinely does not fit */
    }

    /* Break the lattice degeneracy. A perfect lattice puts every neighbouring
       pair at an identical separation, so a large number of collisions come due
       at exactly the same instant; the scheduler sees that as an event avalanche
       and forces an advance. Jitter also avoids starting a FLUID measurement from
       a perfectly ordered crystal. The amplitude is a fraction of the free gap,
       so two neighbours moving toward each other still cannot overlap. */
    {
        double min_gap = 1e30;
        for(int i=0;i<N;i++)
            for(int j=i+1;j<N;j++){
                const double dx2=S->P[i].x-S->P[j].x, dy2=S->P[i].y-S->P[j].y;
                const double dist=sqrt(dx2*dx2+dy2*dy2);
                if(dist < min_gap) min_gap = dist;
            }
        double amp = 0.2 * 0.5 * (min_gap - d);
        if(!(amp > 0.0)) amp = 0.0;
        for(int i=0;i<N;i++){
            const double jx = (((double)(xorshift64(&st)%2000001))/1000000.0 - 1.0)*amp;
            const double jy = (((double)(xorshift64(&st)%2000001))/1000000.0 - 1.0)*amp;
            double nx = S->P[i].x + jx, ny = S->P[i].y + jy;
            if(nx < lo) nx = lo; if(nx > S->prm.boxW-lo) nx = S->prm.boxW-lo;
            if(ny < lo) ny = lo; if(ny > S->prm.boxH-lo) ny = S->prm.boxH-lo;
            S->P[i].x = nx; S->P[i].y = ny;
        }
    }

    /* velocities: same Box-Muller draw as the random seeder */
    for(int i=0;i<N;i++){
        double u1 = ((xorshift64(&st)%1000000)+1)/1000001.0;
        double u2 = ((xorshift64(&st)%1000000)+1)/1000001.0;
        double g  = sqrt(-2.0*log(u1));
        double th = 2.0*M_PI*u2;
        S->P[i].vx = g*cos(th);
        S->P[i].vy = g*sin(th);
        S->P[i].coll_count = 0;
    }
    S->t = 0.0;
    reschedule_all_internal(S);   /* without this the event calendar is empty and
                                     the system is inert: no collisions ever occur */
    return 1;
}

void edmd_init_random_gas(EDMD* S, unsigned long long seed){
    unsigned long long st = seed? seed : 0xC0FFEEULL;
    int N=S->prm.N; double R=S->prm.radius;

    /* positions: naive retries to avoid overlap (okay for demos) */
    for(int i=0;i<N;i++){
        for(;;){
            double rx = (xorshift64(&st)%1000000)/1000000.0;
            double ry = (xorshift64(&st)%1000000)/1000000.0;
            double xx = rx*(S->prm.boxW - 2*R) + R;
            double yy = ry*(S->prm.boxH - 2*R) + R;
            int ok = 1;
            for(int j=0;j<i;j++){
                double dx=xx-S->P[j].x, dy=yy-S->P[j].y;
                if(dx*dx+dy*dy < (2*R)*(2*R)){ ok=0; break; }
            }
            if(ok){ S->P[i].x=xx; S->P[i].y=yy; break; }
        }
        /* velocities: Gaussian with unit variance via Box–Muller */
        double u1 = ((xorshift64(&st)%1000000)+1)/1000001.0;
        double u2 = ((xorshift64(&st)%1000000)+1)/1000001.0;
        double g  = sqrt(-2.0*log(u1));
        double th = 2.0*M_PI*u2;
        S->P[i].vx = g*cos(th);
        S->P[i].vy = g*sin(th);
        S->P[i].coll_count=0;
    }
    S->t=0.0;
    reschedule_all_internal(S);
}

const EDMD_Params*   edmd_params(const EDMD* S){ return &S->prm; }
const EDMD_Particle* edmd_particles(const EDMD* S){ return S->P; }
int                  edmd_count(const EDMD* S){ return S->prm.N; }
double               edmd_time(const EDMD* S){ return S->t; }

/* advance by processing events up to t_target; free-flight if next event is later */
double edmd_advance_to(EDMD* S, double t_target){
    Event e;
    long events_processed = 0;
    long stagnant_events = 0;
    long type_counts[9] = {0};
    double last_event_t = S->t;
    while(S->t < t_target){
        if(!heap_pop(&S->heap, &e)){
            /* no events: free-flight to target */
            double dt = t_target - S->t;
            for(int i=0;i<S->prm.N;i++){
                S->P[i].x += S->P[i].vx*dt;
                S->P[i].y += S->P[i].vy*dt;
            }
            advance_dividers(S, dt);
            S->t = t_target;
            edmd_trace_check(S, "free-flight (heap empty)"); /* ##CHRIS */
            break;
        }
        events_processed++;
        if (e.type >= EV_AB && e.type <= EV_PR) type_counts[(int)e.type]++;
        if (fabs(e.t - last_event_t) <= 1e-13) stagnant_events++;
        else { stagnant_events = 0; last_event_t = e.t; }
        if (events_processed > EDMD_ADVANCE_MAX_EVENTS || stagnant_events > EDMD_ADVANCE_MAX_STAGNANT_EVENTS) {
            g_edmd_avalanche_warning_count++;
            S->forced_advance_count++;
            if (g_edmd_avalanche_warning_count <= 20 || (g_edmd_avalanche_warning_count % 1000) == 0) {
                int dominant = 0;
                for (int k = 1; k <= (int)EV_PR; ++k) {
                    if (type_counts[k] > type_counts[dominant]) dominant = k;
                }
                fprintf(stderr,
                        "[EDMD] warning: event avalanche/stagnation at t=%.17g target=%.17g "
                        "(events=%ld stagnant=%ld dominant=%s:%ld); forcing advance and rescheduling%s.\n",
                        S->t, t_target, events_processed, stagnant_events,
                        ev_name((EvType)dominant), type_counts[dominant],
                        g_edmd_avalanche_warning_count == 20 ? " (further warnings rate-limited)" : "");
            }
            double dt = t_target - S->t;
            if (dt > 0.0) {
                for(int i=0;i<S->prm.N;i++){
                    S->P[i].x += S->P[i].vx*dt;
                    S->P[i].y += S->P[i].vy*dt;
                }
                advance_dividers(S, dt);
                if (S->prm.has_pistonL) S->prm.pistonL_x += S->prm.pistonL_vx * dt;
                if (S->prm.has_pistonR) S->prm.pistonR_x += S->prm.pistonR_vx * dt;
                S->t = t_target;
            }
            edmd_trace_check(S, "forced advance (avalanche/stagnation)"); /* ##CHRIS */
            reschedule_all_internal(S);
            break;
        }
        if(e.t > t_target){
            /* event lies in future: free-flight to t_target and requeue the event */
            double dt = t_target - S->t;
            for(int i=0;i<S->prm.N;i++){
                S->P[i].x += S->P[i].vx*dt;
                S->P[i].y += S->P[i].vy*dt;
            }
            advance_dividers(S, dt);
            if (S->prm.has_pistonL) S->prm.pistonL_x += S->prm.pistonL_vx * dt;
            if (S->prm.has_pistonR) S->prm.pistonR_x += S->prm.pistonR_vx * dt;
            S->t = t_target;
            edmd_trace_check(S, "free-flight (next event beyond target)"); /* ##CHRIS */
            heap_push(&S->heap, e);
            break;
        }
        if(e.t < S->t) continue; /* guard (shouldn’t happen often) */

        /* jump all particles to event time */
        double dt = e.t - S->t;
        for(int i=0;i<S->prm.N;i++){
            S->P[i].x += S->P[i].vx*dt;
            S->P[i].y += S->P[i].vy*dt;
        }
        advance_dividers(S, dt);
        if (S->prm.has_pistonL) S->prm.pistonL_x += S->prm.pistonL_vx * dt;
        if (S->prm.has_pistonR) S->prm.pistonR_x += S->prm.pistonR_vx * dt;
        S->t = e.t;

        /* validate using coll_count snapshots
           ##CHRIS: expressed as a flag rather than three `continue`s so the event
           tracer can log rejected (stale) events too - behaviour is unchanged. */
        int ev_ok = 1;
        if(e.a<0 || e.a>=S->prm.N) ev_ok = 0;
        else if(S->P[e.a].coll_count != e.ca) ev_ok = 0;
        else if(e.type==EV_AB){
            if(e.b<0 || e.b>=S->prm.N) ev_ok = 0;
            else if(S->P[e.b].coll_count != e.cb) ev_ok = 0;
        }
        edmd_trace_record(S, &e, ev_ok);      /* ##CHRIS */
        edmd_trace_check(S, "position jump to event time"); /* ##CHRIS */
        if(!ev_ok) continue;

        /* resolve */
        if(e.type==EV_AB) {
            resolve_ab(S, e.a, e.b);
            grid_build(S);
            schedule_for(S, e.a);
            schedule_for(S, e.b);
            reschedule_clamped(S);   /* ##CHRIS: cover particles grid_build had to bounce */
        } else {
            resolve_wall(S, e.a, e.type, e.b);
            /* moving boundary velocities may have changed (divider/pistons): rebuild and reschedule all */
            if (e.type==EV_DL || e.type==EV_DR || e.type==EV_PL || e.type==EV_PR) {
                reschedule_all_internal(S);
            } else {
                grid_build(S);
                schedule_for(S, e.a);
                reschedule_clamped(S);   /* ##CHRIS */
            }
        }
    }
    return S->t;
}

long edmd_forced_advance_count(const EDMD* S){
    return S ? S->forced_advance_count : 0;
}

/* ##CHRIS: health telemetry. Both should stay 0 in a correct run.
   clamp_repair_count   > 0 : grid_build() had to bounce a particle back into the box.
   overlap_repair_count > 0 : an already-overlapping approaching pair had to be rescued. */

/* ##CHRIS: pressure telemetry. edmd_reset_virial() starts a measurement window
   (call it after equilibration); edmd_compressibility_Z() closes it. Returns NaN
   if the window is empty so a caller can never mistake "no data" for Z=1. */
void edmd_reset_virial(EDMD* S){
    if(!S) return;
    S->virial_accum = 0.0;
    S->virial_pair_events = 0;
    for(int w=0;w<4;w++){ S->wall_impulse[w]=0.0; S->wall_events[w]=0; }
    S->wall_thermal_events = 0;
    S->virial_t0 = S->t;
}
double edmd_virial_accum(const EDMD* S){ return S ? S->virial_accum : 0.0; }
double edmd_wall_impulse(const EDMD* S, int wall){
    return (S && wall>=0 && wall<4) ? S->wall_impulse[wall] : 0.0;
}
long   edmd_wall_events(const EDMD* S, int wall){
    return (S && wall>=0 && wall<4) ? S->wall_events[wall] : 0;
}

/* ##CHRIS: wall-momentum-flux pressure, INDEPENDENT of the pair virial.
     P_x = (I_L + I_R) / (2 H dt),   Z_x = P_x / (rho kB T)
   With rho = N/(W H) and 2D kinetic energy KE = N kB T, rho kB T = KE/(W H), so
     Z_x = (I_L + I_R) * W / (2 dt KE).
   Agreement between Z_x, Z_y and Z_pair is the actual validation; none of the
   three is derived from either of the others. */
static double edmd_ke_total(const EDMD* S){
    double ke=0.0;
    const double m = (S->prm.particle_mass > 0.0) ? S->prm.particle_mass : 1.0;
    for(int i=0;i<S->prm.N;i++) ke += 0.5*m*(S->P[i].vx*S->P[i].vx + S->P[i].vy*S->P[i].vy);
    return ke;
}
double edmd_wall_Z_x(const EDMD* S){
    if(!S) return NAN;
    const double dt = S->t - S->virial_t0; if(!(dt>0.0)) return NAN;
    const double ke = edmd_ke_total(S);    if(!(ke>0.0)) return NAN;
    return (S->wall_impulse[0] + S->wall_impulse[1]) * S->prm.boxW / (2.0 * dt * ke);
}
double edmd_wall_Z_y(const EDMD* S){
    if(!S) return NAN;
    const double dt = S->t - S->virial_t0; if(!(dt>0.0)) return NAN;
    const double ke = edmd_ke_total(S);    if(!(ke>0.0)) return NAN;
    return (S->wall_impulse[2] + S->wall_impulse[3]) * S->prm.boxH / (2.0 * dt * ke);
}
long   edmd_virial_pair_events(const EDMD* S){ return S ? S->virial_pair_events : 0; }
double edmd_virial_window(const EDMD* S){ return S ? (S->t - S->virial_t0) : 0.0; }

double edmd_compressibility_Z(const EDMD* S){
    if(!S) return NAN;
    const double dt = S->t - S->virial_t0;
    if(!(dt > 0.0)) return NAN;
    /* KE = N kB T in 2D with kB = 1, m = 1, so 2*N*kB*T = 2*KE. */
    double ke = 0.0;
    for(int i=0;i<S->prm.N;i++) ke += 0.5*(S->P[i].vx*S->P[i].vx + S->P[i].vy*S->P[i].vy);
    if(!(ke > 0.0)) return NAN;
    return 1.0 + S->virial_accum / (2.0 * ke * dt);
}

long edmd_clamp_repair_count(const EDMD* S){
    return S ? S->clamp_repair_count : 0;
}
long edmd_overlap_repair_count(const EDMD* S){
    return S ? S->overlap_repair_count : 0;
}
long edmd_wall_overdue_count(const EDMD* S){
    return S ? S->wall_overdue_count : 0;
}

void edmd_reschedule_all(EDMD* S){ reschedule_all_internal(S); }

double edmd_advance_pp_only_to(EDMD* S, double t_target){
    Event e;
    while (S->t < t_target) {
        if (!heap_pop(&S->heap, &e)) {
            /* no AB events scheduled: free flight to target */
            double dt = t_target - S->t;
            for (int i = 0; i < S->prm.N; ++i) {
                S->P[i].x += S->P[i].vx * dt;
                S->P[i].y += S->P[i].vy * dt;
            }
            S->t = t_target;
            break;
        }
        if (e.t > t_target) {
            /* free-flight to target and requeue the event */
            double dt = t_target - S->t;
            for (int i = 0; i < S->prm.N; ++i) {
                S->P[i].x += S->P[i].vx * dt;
                S->P[i].y += S->P[i].vy * dt;
            }
            S->t = t_target;
            heap_push(&S->heap, e);
            break;
        }
        if (e.t < S->t) continue;

        /* Jump to event time */
        double dt = e.t - S->t;
        for (int i = 0; i < S->prm.N; ++i) {
            S->P[i].x += S->P[i].vx * dt;
            S->P[i].y += S->P[i].vy * dt;
        }
        S->t = e.t;

        /* Skip any non-AB events defensively (shouldn't exist in PP-only mode) */
        if (e.type != EV_AB) {
            continue;
        }

        if (e.a < 0 || e.a >= S->prm.N) continue;
        if (e.b < 0 || e.b >= S->prm.N) continue;
        if (S->P[e.a].coll_count != e.ca) continue;
        if (S->P[e.b].coll_count != e.cb) continue;

        /* Resolve AB and reschedule for the two participants (AB-only) */
        resolve_ab(S, e.a, e.b);
        schedule_for_pponly(S, e.a);
        schedule_for_pponly(S, e.b);
    }
    return S->t;
}

void edmd_reschedule_all_pp_only(EDMD* S){
    reschedule_all_internal_pponly(S);
}

/* Push particles out of divider slab if overlaps exist; reflect vx in divider frame. */
void edmd_divider_resolve_overlaps(EDMD* S){
    if (!S || S->prm.divider_count <= 0) return;
    int dcount = clamp_dividers(S->prm.divider_count);
    int changed = 0;
    for (int d = 0; d < dcount; ++d) {
        if (S->prm.divider_thickness[d] <= 0.0) continue;
        double L = S->prm.divider_x[d] - 0.5 * S->prm.divider_thickness[d];
        double Rf = S->prm.divider_x[d] + 0.5 * S->prm.divider_thickness[d];
        double vx_div = S->prm.divider_vx[d];
        for (int i = 0; i < S->prm.N; ++i) {
            EDMD_Particle* A = &S->P[i];
            double R = S->prm.radius;
            double left_face  = L - R;
            double right_face = Rf + R;
            if (A->x >= left_face && A->x <= right_face) {
                /* Choose nearest face */
                double gapL = fabs(A->x - left_face);
                double gapR = fabs(right_face - A->x);
                if (gapL < gapR) {
                    A->x = left_face - 1e-9;
                } else {
                    A->x = right_face + 1e-9;
                }
                /* Reflect in divider frame */
                A->vx = 2.0 * vx_div - A->vx;
                A->coll_count++;
                changed = 1;
            }
        }
    }
    if (changed) reschedule_all_internal(S);
}

void edmd_config_dividers(EDMD* S, int count, const double* cx, const double* thickness){
    if (!S) return;
    int n = clamp_dividers(count);
    S->prm.divider_count = n;
    for (int i = 0; i < n; ++i) {
        S->prm.divider_x[i] = cx ? cx[i] : 0.0;
        S->prm.divider_thickness[i] = thickness ? thickness[i] : 0.0;
    }
    for (int i = n; i < EDMD_MAX_DIVIDERS; ++i) {
        S->prm.divider_x[i] = 0.0;
        S->prm.divider_thickness[i] = 0.0;
    }
}

void edmd_set_divider_motions(EDMD* S, int count, const double* mass, const double* vx){
    if (!S) return;
    int n = clamp_dividers(count);
    if (n > S->prm.divider_count) n = S->prm.divider_count;
    for (int i = 0; i < n; ++i) {
        if (mass) S->prm.divider_mass[i] = mass[i];
        if (vx)   S->prm.divider_vx[i] = vx[i];
    }
}

void edmd_set_divider_springs(EDMD* S, int count, const double* k, const double* xeq){
    if (!S) return;
    int n = clamp_dividers(count);
    if (n > S->prm.divider_count) n = S->prm.divider_count;
    for (int i = 0; i < n; ++i) {
        if (k)   S->prm.divider_k[i] = k[i];
        if (xeq) S->prm.divider_xeq[i] = xeq[i];
    }
}

void edmd_config_divider(EDMD* S, int enabled, double cx, double thickness){
    if (!S) return;
    if (enabled) {
        edmd_config_dividers(S, 1, &cx, &thickness);
    } else {
        edmd_config_dividers(S, 0, NULL, NULL);
    }
}

void edmd_set_divider_motion(EDMD* S, double mass, double vx){
    if (!S) return;
    double m = mass, v = vx;
    edmd_set_divider_motions(S, 1, &m, &v);
}

void edmd_set_divider_spring(EDMD* S, double k, double xeq){
    if (!S) return;
    S->prm.divider_k[0] = k;
    S->prm.divider_xeq[0] = xeq;
}

void edmd_config_pistons(EDMD* S,
                         int hasL, double xL, double vxL, double mL,
                         int hasR, double xR, double vxR, double mR)
{
    S->prm.has_pistonL = hasL ? 1 : 0;
    S->prm.pistonL_x = xL;
    S->prm.pistonL_vx = vxL;
    S->prm.pistonL_mass = mL;

    S->prm.has_pistonR = hasR ? 1 : 0;
    S->prm.pistonR_x = xR;
    S->prm.pistonR_vx = vxR;
    S->prm.pistonR_mass = mR;
}

/* total kinetic energy (mass may be 1.0) */
double edmd_total_kinetic_energy(const EDMD* S, double mass){
    double sum=0.0;
    for(int i=0;i<S->prm.N;i++){
        double v2 = S->P[i].vx*S->P[i].vx + S->P[i].vy*S->P[i].vy;
        sum += 0.5 * mass * v2;
    }
    return sum;
}

/* ##CHRIS: Compute gas temperature from kinetic energy (for 2D ideal gas) */
double edmd_gas_temperature(const EDMD* S, double mass, double kB){
    double ke_total = edmd_total_kinetic_energy(S, mass);
    int N = (S->prm.N > 0) ? S->prm.N : 1;
    /* For 2D: T = KE_total / (N * kB) */
    return ke_total / (N * kB);
}

double edmd_work_divider(const EDMD* S){
    if (!S) return 0.0;
    const int dcount = clamp_dividers(S->prm.divider_count);
    double sum = 0.0;
    for (int d = 0; d < dcount; ++d) sum += S->work_divider[d];
    return sum;
}

double edmd_work_divider_i(const EDMD* S, int divider_index){
    if (!S) return 0.0;
    const int dcount = clamp_dividers(S->prm.divider_count);
    if (divider_index < 0 || divider_index >= dcount) return 0.0;
    return S->work_divider[divider_index];
}
double edmd_work_pistonL(const EDMD* S){ return S ? S->work_pistonL : 0.0; }
double edmd_work_pistonR(const EDMD* S){ return S ? S->work_pistonR : 0.0; }
double edmd_heat_bath(const EDMD* S){ return S ? S->heat_bath : 0.0; }
void   edmd_reset_work(EDMD* S){
    if (!S) return;
    for (int d = 0; d < EDMD_MAX_DIVIDERS; ++d) S->work_divider[d] = 0.0;
    S->work_pistonL = 0.0;
    S->work_pistonR = 0.0;
    S->heat_bath = 0.0;
}
