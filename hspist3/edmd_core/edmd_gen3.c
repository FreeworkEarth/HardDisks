/* ##CHRIS 2026-10-08 (261012 sec. 4.7 design, sec. 4.7.1 amendments, sec. 4.7.2 log; milestone M1): the generation-3
   EDMD core. One block comment here; the pieces it replaces in edmd.c (7b08827) are named in sec. 4.7, item 2.

   STATE. Disk i lives in square cell (cx, cy) of width w (>= the diameter; an integer number of px, so cx*w is exact).
   Its position is stored relative to the cell's lower-left corner, (xi, zeta), valid at its own time stamp tau; the
   position at time t is corner + (xi, zeta) + v (t - tau). Times are relative to a floating origin T0 (absolute time =
   T0 + now). Only the disks of an executing event are advanced ("lazy").

   EVENTS. One binary heap, key (t, type, a, b) -- the only arbiter of simultaneous events (sec. 4.7.1, c). Types, in
   tie-break order: CROSS (disk a leaves its cell in direction b: 0 +x, 1 -x, 2 +y, 3 -y), WALL (disk a, outer wall b:
   0 L, 1 R, 2 B, 3 T), PAIR (disks a < b). The values 1, 3, 4 are kept free for the divider band, divider and piston
   events of M2, so adding them will not reorder these. Lazy invalidation: every event carries the collision counters
   (cnt) of its disks at prediction time; a velocity change bumps the counter. A crossing changes no velocity and bumps
   nothing, so a disk's other events survive its crossings.

   WHY NOTHING IS MISSED. Cells are >= one diameter wide, so disks in cells that are not neighbours (index difference >= 2
   in x or y) are >= w >= d apart and cannot touch. Before two disks touch they must become neighbours, which takes a
   crossing of one of them; every disk always has its next crossing in the heap, and at a crossing the disk is predicted
   against the disks of the three cells that just became neighbours (from the current straight lines). A velocity change
   re-predicts the disk against all nine cells. A prediction made while two disks were neighbours stays valid when they
   drift apart (straight lines do not care about cells) until a velocity change invalidates it. A contact at exactly the
   time of a crossing runs after it (CROSS < PAIR). Walls: a disk is predicted against a wall only in a cell from which
   it can reach that wall without leaving the cell, which it must enter by a crossing first.

   NUMBERS. Pair separations are local differences plus an exact multiple of w, so the rounding floor is ulp(w) instead
   of ulp(box) (sec. 4.7, item 3). A crossing re-expresses the coordinate by -w (exact by Sterbenz) or +w (rounded to
   ulp(w)/2), it does not snap it to the boundary (sec. 4.7.2). The origin moves by EDMD3_ORIGIN_SHIFT = 2^13 units once
   now passes it: every disk is synchronised first, then 2^13 is subtracted from now, every tau and every heap time,
   all of which are >= 2^13 then, so the subtraction is exact and the heap order and ties are unchanged.

   TOLERANCES. A pair predicted while overlapping (c = |r|^2 - d^2 < 0) and approaching is scheduled at once: within
   rounding (c >= -C_TOL) it is a contact, beyond it an overlap repair (counted). A wall predicted at or past its face
   while approaching is overdue (counted) and scheduled at once. There is no "t <= 1e-12" cut-off: a disk pair that has
   just collided is not re-predicted while neither has changed velocity since ("mutual last partners" -- receding
   straight lines cannot meet), and a reflected disk moves away from its wall by the sign test.

   READERS NEVER STEER. edmd3_particles() computes positions at the current time without storing them; the validator
   and both audits read the same way. Only the origin shift synchronises (stores), at fixed times. So outputs, audits
   and checks cannot change a trajectory (the audit gate compares hashes with and without them).

   ##CHRIS 2026-10-08 (261012 sec. 4.7.4, sec. 4.7.6; milestone M2): BODIES. A divider d (slab of thickness th centred at
   c) or a piston (one face at c) is a body with a state (c, v) at its own time stamp, a mass M (0 = infinite: held if
   v = 0, else driven), for a divider optionally a spring (k > 0 and M > 0: harmonic about x_eq, as edmd.c), a velocity
   epoch (its DIV and PISTON events carry it; any velocity change bumps it) and a BAND: the cell columns [blo, bhi] that
   meet the contact positions of a disk centre, [lo - h, hi + h] (divider, h = th/2 + R; one side for a piston), where
   [lo, hi] bounds the body's position from now to the band's expiry t_band (the closed-interval test of wall_candidate
   plus a margin). WHY NOTHING IS MISSED: a disk filed in a column outside the band is, until it crosses, inside that
   column (to the cell tolerance), so it cannot touch the body before t_band; it gets the body's event at the crossing
   into a band column (exec_cross) or at the BAND event at t_band that recomputes the band and predicts the new columns
   (BAND < DIV, PISTON at equal times). A velocity change of the body (a collision with a body of finite mass, or an API
   change) bumps the epoch, recomputes the band and re-predicts its disks only: O(sqrt N), not O(N). Expiry: v = 0 (held)
   and a spring at rest: never; constant v: after it moved one cell width; spring: after the arc can have moved one cell
   width, |v| dt + omega^2 A dt^2 / 2 = w. Contact rule (body_rule): with g the face gap, a body of constant velocity is
   met at g / closing speed; a spring divider at the first downward zero of g(t), bracketed between the zeros of g'
   (known in closed form) and solved by safeguarded Newton; g <= 0 while approaching is scheduled at once (a contact
   within rounding if g >= -tol_face, else an overlap repair), except right after a collision of the same disk and body
   ("mutual last", as for pairs), so there is no time cut-off.
   TOLERANCES (amendment b). An executed contact is wrong by at most the time rounding: the event time (ulp(t)/2 x the
   closing speed) and the position of each disk (and body) evaluated at the event and at the prediction (ulp(t)/2 x its
   speed), with t < 2^14 (the origin shift keeps now < 2^13 + one event): a gap error <= 2.5 v_ref u_time, u_time =
   ulp(2^13), v_ref = the largest possible relative speed (energy bound). Hence c_tol = K 2 d v_ref u_time and
   tol_face = K v_ref u_time + 8 ulp(box) with K = EDMD3_TOL_K = 4 (sec. 4.7.6 derives the 2.5 and prints the values).
   LEDGERS (amendment d): momentum and energy of the bodies of finite mass against the impulses and work from outside,
   with the first-order rounding bound accumulated alongside; pure accumulation, never read by the dynamics. */

#include "edmd_gen3.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <float.h>

/* ##CHRIS 2026-10-09 (stage E2, M5; 261012 sec. 4.7.18): THE LONG-DOUBLE BUILD OPTION. With -DEDMD3_LONG_DOUBLE the engine's internal
   arithmetic (positions, velocities, times, predictions, bodies) is long double -- on x86-64 a 64-bit mantissa, so the time quantum
   ulp(2^13) becomes 2^-50 and every rounding tolerance derived from it shrinks by 2^11; on this Mac (arm64) long double is the
   8-byte double, so the numerical check is for KOA. The public API (edmd_gen3.h, EDMD_Particle, EDMD_Params) stays double. The
   default build compiles the same arithmetic as before: real = double, R_(f) = f, and the event hash reads the double value of the
   event time in both builds (a long double has padding bytes). */
#ifdef EDMD3_LONG_DOUBLE
typedef long double real;
#define R_(f) f##l
#define R_MANT_DIG LDBL_MANT_DIG
#else
typedef double real;
#define R_(f) f
#define R_MANT_DIG DBL_MANT_DIG
#endif

enum { T_CROSS = 0, T_BAND = 1, T_WALL = 2, T_DIV = 3, T_PISTON = 4, T_PAIR = 5 };

typedef struct {
    real t;              /* origin-relative event time */
    int    type, a, b;     /* the tie-break key with t */
    int    ca, cb;         /* collision counters at prediction (cb: partner's, PAIR only; else 0) */
    int    pad;            /* always 0: no uninitialised bytes anywhere (sec. 4.7.1, c) */
} Ev3;

typedef struct { Ev3* d; long n, cap; } Heap3;

typedef struct {
    real xi, zeta;       /* position relative to the lower-left corner of the cell [px] */
    real vx, vy;         /* velocity [px per internal unit] */
    real tau;            /* time of (xi, zeta), origin-relative */
    int    cx, cy;         /* cell */
    int    cnt;            /* collision counter */
    int    slot;           /* index in the cell's list */
    int    last;           /* partner of the last velocity change if that was a pair collision, -2 - body if a body
                              collision (M2), else -1 */
    int    pad;
} Disk3;

/* ##CHRIS 2026-10-08 (M2): bodies. Index: dividers 0 .. EDMD_MAX_DIVIDERS-1, then the left and the right piston. */
#define NOBJ    (EDMD_MAX_DIVIDERS + 2)
#define OBJ_PL  EDMD_MAX_DIVIDERS
#define OBJ_PR  (EDMD_MAX_DIVIDERS + 1)
typedef struct {
    int    active, harmonic;   /* harmonic: a divider with k > 0 and M > 0 (edmd.c's divider_has_spring) */
    real M;                  /* mass in disk masses; 0 = infinite (held if v = 0, driven otherwise) */
    real k, xeq, omega;      /* spring constant [kT/px^2], its rest position, sqrt(k/M) */
    real th, h;              /* slab thickness (0 for a piston); h = th/2 + R, the distance face-to-contact centre */
    real x, v, tau;          /* centre (divider) or face (piston) position and velocity at time stamp tau */
    int    epoch, band_gen;    /* velocity epoch (DIV/PISTON events), band generation (BAND events) */
    int    blo, bhi;           /* the band: cell columns [blo, bhi] (empty if blo > bhi) */
    real t_band;             /* the band's expiry (origin-relative), INFINITY if it never expires */
    int    last;               /* disk of the last collision with this body, -1 after an API change */
    int    pad;
    real Jf[2]; long nf[2];  /* x impulse received from disks on face 0 (left side) / 1 (right side), event counts */
    real work;               /* as edmd.c's work_divider / work_pistonL/R */
    real J_inf;              /* ledger: x impulse given to the disks while of mass 0 */
    real J_spring;           /* ledger: x impulse of the spring anchor on the body */
} Obj3;

struct EDMD3 {
    EDMD_Params prm;
    int    N;
    real R, d, d2, w, boxW, boxH;
    int    gw, gh, ncell, ccap;
    int   *cell, *ccount;            /* ncell * ccap disk indices; occupancy */
    Disk3 *D;
    Heap3  heap;
    long   compact_at;
    real T0, now;
    real check_interval, next_check;
    EDMD_Particle* out;
    EDMD3_Health H;
    real tol_pair, tol_wall, tol_cell, c_tol;
    real virial_accum, virial_t0_abs; long virial_pair_events;
    real wall_impulse[4]; long wall_events[4];
    uint64_t hash;
    real same_t; long same_n, same_limit;
    int    fatal; char fatal_msg[256];
    /* audits */
    int    contact_audit; real contact_max[2]; long contact_events;
    long   audit_every, audit_count; EDMD3_Audit A;
    int    audit_bodies;             /* M2: also audit right after every BAND event and every API body change */
    real *alx, *aly;               /* audit scratch: local position of every disk at now */
    real *acr, *awl;               /* live crossing time per disk (and direction in acd), live wall time per (disk, side) */
    int    *acd;
    long   ahcap; uint64_t* ahk; real* aht; unsigned char* ahm;   /* pair table: key, time, matched */
    long   aprinted;
    /* M2: bodies, tolerances, ledgers, the body audit's scratch */
    Obj3   obj[NOBJ];
    int    objs[NOBJ], nobj;         /* the active bodies, in index order */
    real band_margin, tol_face, u_time, E_bound, v_ref, m_min;
    real ledJ[2], ledJw[4], ledW, J_api, ledSP[2], ledSE;   /* ledger accumulators (amendment d) */
    real P0[2], E0, P0s[2], E0s;  /* at load, and the rounding scale of those sums */
    real contact_max_obj[2];       /* contact audit: [0] divider faces, [1] pistons */
    real *aot; int *aob;           /* body audit: live heap time and event code per (disk, active body) */
    /* ##CHRIS 2026-10-09 (M3, 261012 sec. 4.7.14): edmd.c's gated event log (edmd_set_event_log), same columns and values,
       one row per outer-wall, divider-face and piston collision. PRINT ONLY: written from the values the collision rule
       already computed; NULL (off) unless the driver sets it. */
    FILE*  evlog; real evlog_tscale;
};

/* M3: one row of edmd.c's edmd_log_event (t_sigma, kind, u_wall, v_before, v_after, dE, dp), the same format */
static const char* const WALL_NAME[4] = {"WL", "WR", "WB", "WT"};
static void evlog_row(const EDMD3* S, const char* kind, real u, real v0, real v1, real dE){
    fprintf(S->evlog, "%.12g,%s,%.12g,%.12g,%.12g,%.12g,%.12g\n", (double)((S->T0 + S->now) / S->evlog_tscale), kind, (double)u, (double)v0,
            (double)v1, (double)dE, (double)(1.0 * (v1 - v0)));   /* ##CHRIS (E2): printed as double in both builds */
}

/* ------------------------------------------------------------------ heap, ordered by (t, type, a, b) */

static inline int ev_less(const Ev3* x, const Ev3* y){
    if (x->t != y->t) return x->t < y->t;
    if (x->type != y->type) return x->type < y->type;
    if (x->a != y->a) return x->a < y->a;
    return x->b < y->b;
}
static void heap_push(EDMD3* S, real t, int type, int a, int b, int ca, int cb){
    Heap3* H = &S->heap;
    if (H->n >= H->cap) {
        long nc = H->cap + H->cap / 2 + 1024;
        Ev3* nd = (Ev3*)realloc(H->d, (size_t)nc * sizeof(Ev3));
        if (!nd) { fprintf(stderr, "edmd_gen3: out of memory (heap)\n"); exit(3); }
        H->d = nd; H->cap = nc;
    }
    Ev3 e; e.t = t; e.type = type; e.a = a; e.b = b; e.ca = ca; e.cb = cb; e.pad = 0;
    long i = H->n++;
    while (i > 0) {
        long p = (i - 1) >> 1;
        if (!ev_less(&e, &H->d[p])) break;
        H->d[i] = H->d[p]; i = p;
    }
    H->d[i] = e;
    if (H->n > S->H.heap_max) S->H.heap_max = H->n;
}
static void heap_sift_down(Heap3* H, long i){
    const Ev3 e = H->d[i];
    for (;;) {
        long l = 2 * i + 1;
        if (l >= H->n) break;
        long r = l + 1, s = (r < H->n && ev_less(&H->d[r], &H->d[l])) ? r : l;
        if (!ev_less(&H->d[s], &e)) break;
        H->d[i] = H->d[s]; i = s;
    }
    H->d[i] = e;
}
static void heap_pop(Heap3* H, Ev3* out){
    *out = H->d[0];
    H->d[0] = H->d[--H->n];
    if (H->n > 0) heap_sift_down(H, 0);
}

static inline int ev_live(const EDMD3* S, const Ev3* e){
    if (e->type == T_BAND) return S->obj[e->a].band_gen == e->ca;          /* M2: a = body, ca = band generation */
    if (S->D[e->a].cnt != e->ca) return 0;
    if (e->type == T_PAIR && S->D[e->b].cnt != e->cb) return 0;
    if (e->type == T_DIV) return S->obj[e->b >> 1].epoch == e->cb;           /* M2: b = 2 d + face, cb = epoch */
    if (e->type == T_PISTON) return S->obj[OBJ_PL + e->b].epoch == e->cb;   /* M2: b = side */
    return 1;
}

/* drop stale events and re-heapify (Floyd); the pop order is the key order, so this changes nothing observable */
static void heap_compact(EDMD3* S){
    Heap3* H = &S->heap; long n = 0;
    for (long k = 0; k < H->n; ++k) if (ev_live(S, &H->d[k])) H->d[n++] = H->d[k];
    H->n = n;
    for (long i = n / 2 - 1; i >= 0; --i) heap_sift_down(H, i);
    const long floor_n = 64L * S->N + 4096;
    S->compact_at = (2 * n > floor_n) ? 2 * n : floor_n;
    S->H.heap_compactions++;
}

/* ------------------------------------------------------------------ cells */

static inline int cell_of(const EDMD3* S, int cx, int cy){ return cy * S->gw + cx; }
static int cell_insert(EDMD3* S, int i){
    Disk3* A = &S->D[i]; const int c = cell_of(S, A->cx, A->cy);
    if (S->ccount[c] >= S->ccap) return 0;
    A->slot = S->ccount[c];
    S->cell[(long)c * S->ccap + S->ccount[c]++] = i;
    return 1;
}
static void cell_remove(EDMD3* S, int i){
    Disk3* A = &S->D[i]; const int c = cell_of(S, A->cx, A->cy);
    int* L = &S->cell[(long)c * S->ccap];
    const int last = L[--S->ccount[c]];
    L[A->slot] = last; S->D[last].slot = A->slot;
    A->slot = -1;
}
/* can a disk in cell column cx (row cy) reach outer wall s without leaving the cell? (closed cell intervals) */
static inline int wall_candidate(const EDMD3* S, int cx, int cy, int s){
    switch (s) {
        case 0: return (real)cx * S->w <= S->R;
        case 1: return (real)(cx + 1) * S->w >= S->boxW - S->R;
        case 2: return (real)cy * S->w <= S->R;
        default: return (real)(cy + 1) * S->w >= S->boxH - S->R;
    }
}

/* ------------------------------------------------------------------ motion */

static inline void advance(EDMD3* S, int i){
    Disk3* A = &S->D[i];
    const real dt = S->now - A->tau;
    if (dt != 0.0) { A->xi += A->vx * dt; A->zeta += A->vy * dt; A->tau = S->now; }
}
/* local position of disk j at now, not stored */
static inline void local_now(const EDMD3* S, int j, real* x, real* y){
    const Disk3* B = &S->D[j]; const real dt = S->now - B->tau;
    *x = B->xi + B->vx * dt; *y = B->zeta + B->vy * dt;
}
static inline void fnv(uint64_t* h, const void* p, size_t n){
    const unsigned char* c = (const unsigned char*)p;
    for (size_t k = 0; k < n; ++k) { *h ^= c[k]; *h *= 1099511628211ULL; }
}

/* ------------------------------------------------------------------ bodies (M2): motion, bands, the contact rule */

static const real U_ROUND = 1.1102230246251565e-16;   /* 2^-53, the unit roundoff of the ledger scales */
#define G3_PI 3.14159265358979323846     /* M_PI is not in strict C11 */

/* position and velocity of body O at origin-relative time t (computed, not stored); exact at its own time stamp */
static inline void obj_at(const Obj3* O, real t, real* x, real* v){
    const real dt = t - O->tau;
    if (dt == 0.0) { *x = O->x; *v = O->v; return; }
    if (O->harmonic) {             /* as edmd.c's harmonic_advance_1d */
        const real c = R_(cos)(O->omega * dt), s = R_(sin)(O->omega * dt), dx = O->x - O->xeq;
        *x = O->xeq + dx * c + (O->v / O->omega) * s;
        *v = -dx * O->omega * s + O->v * c;
    } else { *x = O->x + O->v * dt; *v = O->v; }
}
static inline real obj_amp(const Obj3* O, real x, real v){ return R_(hypot)(x - O->xeq, v / O->omega); }
/* store body o's state at now; a spring's anchor impulse goes into the ledger */
static void obj_advance(EDMD3* S, int o){
    Obj3* O = &S->obj[o];
    if (O->tau == S->now) return;
    real x, v; obj_at(O, S->now, &x, &v);
    if (O->harmonic) {
        const real A = obj_amp(O, O->x, O->v), dp = O->M * (v - O->v);
        O->J_spring += dp; S->ledJ[0] += dp;
        S->ledSP[0] += 8.0 * O->M * O->omega * A + R_(fabs)(dp) + R_(fabs)(S->ledJ[0]);
        S->ledSE += 8.0 * O->k * A * A;
    }
    O->x = x; O->v = v; O->tau = S->now;
}
static inline int in_band(const Obj3* O, int cx){ return O->blo <= cx && cx <= O->bhi; }
/* side of a disk centre at absolute x relative to body o: +1 = the body is to its right (sigma of the gap) */
static inline int body_side(int o, real c, real xp){
    if (o < EDMD_MAX_DIVIDERS) return xp < c ? 1 : -1;
    return o == OBJ_PL ? -1 : 1;
}
/* the face gap g = sigma (c - xp) - h of a disk centre xp and body o at position c */
static inline real body_gap(const Obj3* O, int o, real c, real xp){ return (real)body_side(o, c, xp) * (c - xp) - O->h; }

/* spring divider: g(t) = G + sg (D cos wt + (V/w) sin wt - vp t) and its derivative */
typedef struct { real G, sg, D, V, w, vp; } Harm;
static inline real harm_g(const Harm* H, real t){ return H->G + H->sg * (H->D * R_(cos)(H->w * t) + (H->V / H->w) * R_(sin)(H->w * t) - H->vp * t); }
static inline real harm_gp(const Harm* H, real t){ return H->sg * (-H->D * H->w * R_(sin)(H->w * t) + H->V * R_(cos)(H->w * t) - H->vp); }
/* the zero of g in [lo, hi], where g decreases from > 0 to <= 0: safeguarded Newton (a step that leaves the bracket
   is replaced by bisection), to a step of 4 eps |t| or adjacent bracket ends (then the end with g <= 0) */
static real harm_root(const Harm* H, real lo, real hi){
    real t = lo + 0.5 * (hi - lo);
    for (int it = 0; it < 400; ++it) {
        const real g = harm_g(H, t);
        if (g <= 0.0) hi = t; else lo = t;
        const real gp = harm_gp(H, t);
        real tn = (gp < 0.0) ? t - g / gp : lo + 0.5 * (hi - lo);
        if (!(tn > lo && tn < hi)) tn = lo + 0.5 * (hi - lo);
        if (!(tn > lo && tn < hi)) return hi;          /* lo and hi are adjacent doubles */
        if (R_(fabs)(tn - t) <= 4.0 * 2.220446049250313e-16 * R_(fabs)(t)) return tn;
        t = tn;
    }
    return hi;
}
/* the contact rule shared by the engine and the schedule audit. Disk centre xp, velocity vp, at now; body o at (c, vo).
   0: no contact; 1: contact after *dt; 2: at once (g <= 0 and approaching, not right after a collision of this disk and
   this body: mlast). *g0 = the gap now. */
static int body_rule(const Obj3* O, int o, real c, real vo, real xp, real vp, int mlast, real* dt, real* g0){
    const int sgi = body_side(o, c, xp);
    const real sg = (real)sgi;
    *g0 = sg * (c - xp) - O->h;
    if (!O->harmonic) {
        const real s = sg * (vp - vo);                   /* closing speed: g' = -s */
        if (!(s > 0.0)) return 0;
        if (*g0 <= 0.0) { if (mlast) return 0; *dt = 0.0; return 2; }
        *dt = *g0 / s; return 1;
    }
    Harm H; H.sg = sg; H.D = c - O->xeq; H.V = vo; H.w = O->omega; H.vp = vp; H.G = sg * (O->xeq - xp) - O->h;
    if (*g0 <= 0.0 && sg * (vo - vp) < 0.0 && !mlast) { *dt = 0.0; return 2; }
    const real T = 2.0 * G3_PI / H.w, Q = R_(hypot)(H.D * H.w, H.V), A = Q / H.w;
    const real drift = sg * vp;            /* g(t + T) = g(t) - drift T exactly: drift > 0 iff the disk approaches the face */
    if (!(Q > R_(fabs)(vp))) {                   /* g' = sg Q cos(w t + phi) - drift has the sign of -drift everywhere: monotone */
        if (!(drift > 0.0)) return 0;
        if (!(harm_g(&H, 0.0) > 0.0)) return 0;
        real tb = R_(fmax)(0.0, (sg * (O->xeq - xp) + A - O->h) / drift) + T;   /* the disk is past the face's extreme by then */
        while (harm_g(&H, tb) > 0.0) tb *= 2.0;                                  /* rounding guard */
        *dt = harm_root(&H, 0.0, tb); return 1;
    }
    /* g oscillates: in each period it falls from a local max (w t + phi = sg acos(vp/Q)) to a local min (w t + phi =
       -sg acos(vp/Q)), phi = atan2(D w, V), and every extremum is drift T lower one period later */
    const real phi = R_(atan2)(H.D * H.w, H.V), al = acos(vp / Q), tp = 2.0 * G3_PI;
    real tmx = R_(fmod)(-phi + sg * al, tp), tmn = R_(fmod)(-phi - sg * al, tp);
    if (tmx < 0.0) tmx += tp;
    if (tmn < 0.0) tmn += tp;
    tmx /= H.w; tmn /= H.w;                                    /* the first local max and min in [0, T) */
    const real ddec = tmn > tmx ? tmn - tmx : tmn - tmx + T; /* the length of a falling stretch */
    /* the first two periods, monotone piece by piece */
    real ta = 0.0, ga = harm_g(&H, 0.0);
    for (int n = 0; n < 2; ++n)
        for (int k = 0; k < 2; ++k) {
            const real tb = (real)n * T + (k == 0 ? R_(fmin)(tmx, tmn) : R_(fmax)(tmx, tmn));
            if (tb <= ta) continue;
            const real gb = harm_g(&H, tb);
            if (ga > 0.0 && gb <= 0.0) { *dt = harm_root(&H, ta, tb); return 1; }
            ta = tb; ga = gb;
        }
    if (!(drift > 0.0)) return 0;            /* the minima do not fall: none after the first period touches zero */
    /* the minima fall by drift T per period: jump to the first one at or below zero (O(1) for any slow approach) */
    const real t1 = tmn + T, m1 = harm_g(&H, t1);           /* the local min in [T, 2T), scanned above */
    if (!(m1 > 0.0)) return 0;                                 /* g <= 0 throughout: an overlapping state (the validator's) */
    real kk = R_(ceil)(m1 / (drift * T));
    if (kk < 1.0) kk = 1.0;
    while (kk > 1.0 && harm_g(&H, t1 + (kk - 1.0) * T) <= 0.0) kk -= 1.0;   /* rounding guards */
    while (harm_g(&H, t1 + kk * T) > 0.0) kk += 1.0;
    const real hi = t1 + kk * T, lo = hi - ddec;            /* the falling stretch that ends at that minimum */
    if (!(harm_g(&H, lo) > 0.0)) return 0;
    *dt = harm_root(&H, lo, hi); return 1;
}

/* ------------------------------------------------------------------ predictions (disk i is synchronised: tau_i == now) */

/* the pair rule shared with the brute-force audit: 0 none, 1 future contact at dt, 2 at once (c < 0, approaching) */
static inline int pair_rule(real rx, real ry, real vx, real vy, real d2, real* dt, real* c_out){
    const real b = rx * vx + ry * vy;
    if (b >= 0.0) return 0;
    const real vv = vx * vx + vy * vy;
    const real c = rx * rx + ry * ry - d2;
    *c_out = c;
    if (c <= 0.0) { *dt = 0.0; return 2; }
    const real disc = b * b - vv * c;
    if (disc <= 0.0) return 0;
    *dt = (-b - R_(sqrt)(disc)) / vv;
    if (*dt < 0.0) *dt = 0.0;
    return 1;
}
static inline int mutual_last(const EDMD3* S, int i, int j){ return S->D[i].last == j && S->D[j].last == i; }

static void predict_pair(EDMD3* S, int i, int j){
    if (mutual_last(S, i, j)) return;
    const Disk3 *A = &S->D[i], *B = &S->D[j];
    real bx, by; local_now(S, j, &bx, &by);
    const real rx = (bx - A->xi) + (real)(B->cx - A->cx) * S->w;
    const real ry = (by - A->zeta) + (real)(B->cy - A->cy) * S->w;
    real dt = 0.0, c = 0.0;
    const int rc = pair_rule(rx, ry, B->vx - A->vx, B->vy - A->vy, S->d2, &dt, &c);
    if (!rc) return;
    if (rc == 2) {   /* ##CHRIS 2026-10-08 (M2, amendment b): c_tol from the time resolution (tol_update) */
        if (c < -S->c_tol) S->H.overlap_repair++;
        else { S->H.contact_now++; if (c < S->H.contact_c_min) S->H.contact_c_min = c; }
    }
    const int a = i < j ? i : j, b = i < j ? j : i;
    heap_push(S, S->now + dt, T_PAIR, a, b, S->D[a].cnt, S->D[b].cnt);
}
/* wall gap of disk i (synchronised) to wall s, and the closing speed; 0 if moving away or parallel */
static inline int wall_gap(const EDMD3* S, const Disk3* A, real xl, real yl, int s, real* gap, real* speed){
    switch (s) {
        case 0: if (A->vx >= 0.0) return 0; *gap = ((real)A->cx * S->w - S->R) + xl; *speed = -A->vx; return 1;
        case 1: if (A->vx <= 0.0) return 0; *gap = (S->boxW - S->R - (real)A->cx * S->w) - xl; *speed = A->vx; return 1;
        case 2: if (A->vy >= 0.0) return 0; *gap = ((real)A->cy * S->w - S->R) + yl; *speed = -A->vy; return 1;
        default: if (A->vy <= 0.0) return 0; *gap = (S->boxH - S->R - (real)A->cy * S->w) - yl; *speed = A->vy; return 1;
    }
}
static void predict_wall(EDMD3* S, int i, int s){
    const Disk3* A = &S->D[i]; real gap, speed;
    if (!wall_gap(S, A, A->xi, A->zeta, s, &gap, &speed)) return;
    real dt = 0.0;
    if (gap <= 0.0) { if (gap < -S->tol_face) S->H.wall_overdue++; else S->H.wall_contact_now++; }   /* ##CHRIS M2: split by tol_face */
    else dt = gap / speed;
    heap_push(S, S->now + dt, T_WALL, i, s, A->cnt, 0);
}
/* next crossing from local coordinates (x, y) and the velocity: direction and time from now; -1 if none */
static inline int cross_rule(const EDMD3* S, const Disk3* A, real x, real y, real* dt){
    real tx = INFINITY, ty = INFINITY; int dx = -1, dy = -1;
    if (A->vx > 0.0 && A->cx + 1 < S->gw) { tx = (S->w - x) / A->vx; dx = 0; }
    else if (A->vx < 0.0 && A->cx > 0)    { tx = x / (-A->vx); dx = 1; }
    if (A->vy > 0.0 && A->cy + 1 < S->gh) { ty = (S->w - y) / A->vy; dy = 2; }
    else if (A->vy < 0.0 && A->cy > 0)    { ty = y / (-A->vy); dy = 3; }
    if (dx < 0 && dy < 0) return -1;
    int dir; real t;
    if (dy < 0 || (dx >= 0 && tx <= ty)) { dir = dx; t = tx; } else { dir = dy; t = ty; }
    *dt = t > 0.0 ? t : 0.0;
    return dir;
}
static void predict_cross(EDMD3* S, int i){
    const Disk3* A = &S->D[i]; real dt;
    const int dir = cross_rule(S, A, A->xi, A->zeta, &dt);
    if (dir >= 0) heap_push(S, S->now + dt, T_CROSS, i, dir, A->cnt, 0);
}
static void predict_cell(EDMD3* S, int i, int cx, int cy){
    if (cx < 0 || cy < 0 || cx >= S->gw || cy >= S->gh) return;
    const int c = cell_of(S, cx, cy); const int* L = &S->cell[(long)c * S->ccap];
    for (int k = 0; k < S->ccount[c]; ++k) if (L[k] != i) predict_pair(S, i, L[k]);
}
/* M2: disk i, local x coordinate xl at now, against body o */
static void predict_obj(EDMD3* S, int i, int o, real xl){
    const Disk3* A = &S->D[i]; const Obj3* O = &S->obj[o];
    real c, vo; obj_at(O, S->now, &c, &vo);
    const real xp = (real)A->cx * S->w + xl;
    const int mlast = A->last == -2 - o && O->last == i;
    real dt = 0.0, g0 = 0.0;
    const int rc = body_rule(O, o, c, vo, xp, A->vx, mlast, &dt, &g0);
    if (!rc) return;
    if (rc == 2) {
        if (g0 < -S->tol_face) S->H.obj_overlap_repair++;
        else { S->H.obj_contact_now++; if (g0 < S->H.obj_contact_gap_min) S->H.obj_contact_gap_min = g0; }
    }
    if (o < EDMD_MAX_DIVIDERS) heap_push(S, S->now + dt, T_DIV, i, 2 * o + (body_side(o, c, xp) > 0 ? 0 : 1), A->cnt, O->epoch);
    else heap_push(S, S->now + dt, T_PISTON, i, o - OBJ_PL, A->cnt, O->epoch);
}
/* after a velocity change of i (or at load): its crossing, its walls, all nine cells; M2: the bodies whose band holds it */
static void predict_all(EDMD3* S, int i){
    const Disk3* A = &S->D[i];
    predict_cross(S, i);
    for (int s = 0; s < 4; ++s) if (wall_candidate(S, A->cx, A->cy, s)) predict_wall(S, i, s);
    for (int dy = -1; dy <= 1; ++dy) for (int dx = -1; dx <= 1; ++dx) predict_cell(S, i, A->cx + dx, A->cy + dy);
    for (int q = 0; q < S->nobj; ++q) if (in_band(&S->obj[S->objs[q]], A->cx)) predict_obj(S, i, S->objs[q], A->xi);
}

/* M2: the band of body o from its state at now: the position bound [lo, hi] until the expiry, the contact positions of a
   disk centre, the columns that meet them (closed test plus the margin); pushes the BAND event if the band expires */
static void band_bounds(const EDMD3* S, int o, real t0, real* lo, real* hi, real* T){
    const Obj3* O = &S->obj[o];
    real x, v; obj_at(O, t0, &x, &v);
    *lo = x; *hi = x; *T = INFINITY;
    if (O->harmonic) {
        const real A = obj_amp(O, x, v);
        if (A > 0.0) {
            const real acc = O->omega * O->omega * A;   /* |x''| <= omega^2 A: |x(t0+s) - x| <= |v| s + acc s^2 / 2 */
            const real ds = 2.0 * S->w / (R_(fabs)(v) + R_(sqrt)(v * v + 2.0 * acc * S->w));
            *lo = R_(fmin)(x, R_(fmax)(x - S->w, O->xeq - A)); *hi = R_(fmax)(x, R_(fmin)(x + S->w, O->xeq + A));
            *T = t0 + ds;
        }
    } else if (v != 0.0) {
        const real ds = S->w / R_(fabs)(v), x1 = x + v * ds;
        *lo = R_(fmin)(x, x1); *hi = R_(fmax)(x, x1); *T = t0 + ds;
    }
}
static void band_columns(const EDMD3* S, int o, real lo, real hi, int* blo, int* bhi){
    const Obj3* O = &S->obj[o];
    real X0, X1;
    if (o < EDMD_MAX_DIVIDERS) { X0 = lo - O->h; X1 = hi + O->h; }
    else if (o == OBJ_PL) { X0 = lo + O->h; X1 = hi + O->h; }
    else { X0 = lo - O->h; X1 = hi - O->h; }
    /* column c meets [X0, X1] (closed, margin m) iff (c+1) w >= X0 - m and c w <= X1 + m */
    real a = R_(ceil)((X0 - S->band_margin) / S->w - 1.0), b = R_(floor)((X1 + S->band_margin) / S->w);
    if (a < 0.0) a = 0.0;
    if (b > (real)(S->gw - 1)) b = (real)(S->gw - 1);
    if (a > b) { *blo = 1; *bhi = 0; return; }       /* the band misses the grid: empty */
    *blo = (int)a; *bhi = (int)b;
}
static void band_compute(EDMD3* S, int o){
    Obj3* O = &S->obj[o];
    real lo, hi, T;
    band_bounds(S, o, S->now, &lo, &hi, &T);
    band_columns(S, o, lo, hi, &O->blo, &O->bhi);
    O->t_band = T;
    O->band_gen++;
    if (isfinite(T)) heap_push(S, T, T_BAND, o, 0, O->band_gen, 0);
}
/* body o against every disk filed in column cx, except disk skip */
static void predict_column(EDMD3* S, int o, int cx, int skip){
    for (int cy = 0; cy < S->gh; ++cy) {
        const int c = cell_of(S, cx, cy); const int* L = &S->cell[(long)c * S->ccap];
        for (int k = 0; k < S->ccount[c]; ++k) {
            const int j = L[k]; if (j == skip) continue;
            real xl, yl; local_now(S, j, &xl, &yl);
            predict_obj(S, j, o, xl);
        }
    }
}
/* body o changed velocity or mass at now (it is synchronised): new epoch, new band, its band disks re-predicted */
static void obj_changed(EDMD3* S, int o, int skip){
    Obj3* O = &S->obj[o];
    O->epoch++;
    band_compute(S, o);
    for (int cx = O->blo; cx <= O->bhi; ++cx) predict_column(S, o, cx, skip);
}

/* ------------------------------------------------------------------ validator (sec. 4.7.1, h); read-only except the box repair */

static void reflect_into_box(EDMD3* S, int i);
/* the event's disk i (synchronised) against every disk of its nine cells and the four walls */
static void local_check(EDMD3* S, int i){
    Disk3* A = &S->D[i]; real worst = 0.0; int found = 0;
    for (int dy = -1; dy <= 1; ++dy) for (int dx = -1; dx <= 1; ++dx) {
        const int cx = A->cx + dx, cy = A->cy + dy;
        if (cx < 0 || cy < 0 || cx >= S->gw || cy >= S->gh) continue;
        const int c = cell_of(S, cx, cy); const int* L = &S->cell[(long)c * S->ccap];
        for (int k = 0; k < S->ccount[c]; ++k) {
            const int j = L[k]; if (j == i) continue;
            real bx, by; local_now(S, j, &bx, &by);
            const real rx = (bx - A->xi) + (real)dx * S->w, ry = (by - A->zeta) + (real)dy * S->w;
            const real gap = R_(sqrt)(rx * rx + ry * ry) - S->d;
            if (gap < worst) worst = gap;
            if (gap < -S->tol_pair) found = 1;
        }
    }
    const real x = (real)A->cx * S->w + A->xi, y = (real)A->cy * S->w + A->zeta;
    const real g[4] = { x - S->R, (S->boxW - S->R) - x, y - S->R, (S->boxH - S->R) - y };
    int out = 0;
    for (int s = 0; s < 4; ++s) { if (g[s] < worst) worst = g[s]; if (g[s] < -S->tol_wall) out = 1; }
    for (int q = 0; q < S->nobj; ++q) {       /* M2: the bodies (a finding only; no repair) */
        const int o = S->objs[q]; real c, vo; obj_at(&S->obj[o], S->now, &c, &vo);
        const real gb = body_gap(&S->obj[o], o, c, x);
        if (gb < worst) worst = gb;
        if (gb < -S->tol_wall) found = 1;
    }
    if (A->xi < -S->tol_cell || A->xi > S->w + S->tol_cell || A->zeta < -S->tol_cell || A->zeta > S->w + S->tol_cell) {
        /* outside its own cell: a bookkeeping failure; re-file it from its absolute position */
        S->H.cell_repair++;
        cell_remove(S, i);
        int ncx = (int)R_(floor)(x / S->w), ncy = (int)R_(floor)(y / S->w);
        ncx = ncx < 0 ? 0 : (ncx >= S->gw ? S->gw - 1 : ncx);
        ncy = ncy < 0 ? 0 : (ncy >= S->gh ? S->gh - 1 : ncy);
        A->xi = x - (real)ncx * S->w; A->zeta = y - (real)ncy * S->w; A->cx = ncx; A->cy = ncy;
        if (!cell_insert(S, i)) { S->fatal = 1; snprintf(S->fatal_msg, sizeof S->fatal_msg, "cell overflow at re-filing disk %d", i); }
        A->cnt++; A->last = -1;
        predict_all(S, i);
    }
    S->H.local_checks++;
    if (found || out) S->H.local_findings++;
    if (worst < S->H.local_worst) S->H.local_worst = worst;
    if (out) reflect_into_box(S, i);
}
static void reflect_into_box(EDMD3* S, int i){
    Disk3* A = &S->D[i];
    real x = (real)A->cx * S->w + A->xi, y = (real)A->cy * S->w + A->zeta;
    if (x < S->R)           { x = S->R;           if (A->vx < 0.0) A->vx = -A->vx; }
    if (x > S->boxW - S->R) { x = S->boxW - S->R; if (A->vx > 0.0) A->vx = -A->vx; }
    if (y < S->R)           { y = S->R;           if (A->vy < 0.0) A->vy = -A->vy; }
    if (y > S->boxH - S->R) { y = S->boxH - S->R; if (A->vy > 0.0) A->vy = -A->vy; }
    A->xi = x - (real)A->cx * S->w; A->zeta = y - (real)A->cy * S->w;
    S->H.clamp_repair++; A->cnt++; A->last = -1;
    predict_all(S, i);
}
/* every pair within nine cells and every wall, positions computed at now (not stored) */
static void full_check(EDMD3* S){
    real worst = 0.0; long found = 0;
    for (int i = 0; i < S->N; ++i) {
        const Disk3* A = &S->D[i]; real ax, ay; local_now(S, i, &ax, &ay);
        for (int dy = -1; dy <= 1; ++dy) for (int dx = -1; dx <= 1; ++dx) {
            const int cx = A->cx + dx, cy = A->cy + dy;
            if (cx < 0 || cy < 0 || cx >= S->gw || cy >= S->gh) continue;
            const int c = cell_of(S, cx, cy); const int* L = &S->cell[(long)c * S->ccap];
            for (int k = 0; k < S->ccount[c]; ++k) {
                const int j = L[k]; if (j <= i) continue;
                real bx, by; local_now(S, j, &bx, &by);
                const real rx = (bx - ax) + (real)dx * S->w, ry = (by - ay) + (real)dy * S->w;
                const real gap = R_(sqrt)(rx * rx + ry * ry) - S->d;
                if (gap < worst) worst = gap;
                if (gap < -S->tol_pair) found++;
            }
        }
        const real x = (real)A->cx * S->w + ax, y = (real)A->cy * S->w + ay;
        const real g[4] = { x - S->R, (S->boxW - S->R) - x, y - S->R, (S->boxH - S->R) - y };
        for (int s = 0; s < 4; ++s) { if (g[s] < worst) worst = g[s]; if (g[s] < -S->tol_wall) found++; }
        for (int q = 0; q < S->nobj; ++q) {   /* M2: the bodies */
            const int o = S->objs[q]; real c, vo; obj_at(&S->obj[o], S->now, &c, &vo);
            const real gb = body_gap(&S->obj[o], o, c, x);
            if (gb < worst) worst = gb;
            if (gb < -S->tol_wall) found++;
        }
    }
    /* M2: every body inside the box, no two slabs overlapping, every slab between the pistons, the pistons in order */
    if (S->nobj) {
        real l[NOBJ], r[NOBJ];
        for (int q = 0; q < S->nobj; ++q) {
            const int o = S->objs[q]; const Obj3* O = &S->obj[o]; real c, vo; obj_at(O, S->now, &c, &vo);
            l[q] = c - 0.5 * O->th; r[q] = c + 0.5 * O->th;
            /* ##CHRIS 2026-10-09 (M3, 261012 sec. 4.7.14, decision log): a piston PARKED behind its own wall is valid -- the
               driver's energy-transfer setup parks them at x = -1 px and boxW + 6 px, as edmd.c allows: no disk can reach
               such a face before its wall, which reflects it first. The left piston only must not be right of the box, the
               right one not left of it; dividers stay inside the box; the order checks below are unchanged. */
            if (o == OBJ_PL) { if (r[q] > S->boxW + S->tol_wall) S->H.body_findings++; }
            else if (o == OBJ_PR) { if (l[q] < -S->tol_wall) S->H.body_findings++; }
            else if (l[q] < -S->tol_wall || r[q] > S->boxW + S->tol_wall) S->H.body_findings++;
        }
        for (int q = 0; q < S->nobj; ++q) for (int p = q + 1; p < S->nobj; ++p) {
            const int a = S->objs[q], b = S->objs[p];
            if (a < EDMD_MAX_DIVIDERS && b < EDMD_MAX_DIVIDERS) { if (l[q] < r[p] - S->tol_wall && l[p] < r[q] - S->tol_wall) S->H.body_findings++; }
            else if (b == OBJ_PL) { if (l[q] < l[p] - S->tol_wall) S->H.body_findings++; }    /* a divider left of the left piston */
            else if (b == OBJ_PR) { if (r[q] > r[p] + S->tol_wall) S->H.body_findings++; }    /* right of the right piston (or PL > PR) */
        }
    }
    S->H.full_checks++; S->H.full_findings += found;
    if (worst < S->H.full_worst) S->H.full_worst = worst;
}

/* ------------------------------------------------------------------ event execution */

static void contact_audit(EDMD3* S, const Ev3* e){
    const Disk3* A = &S->D[e->a]; real gap; int k;
    if (e->type == T_PAIR) {
        const Disk3* B = &S->D[e->b];
        const real rx = (B->xi - A->xi) + (real)(B->cx - A->cx) * S->w, ry = (B->zeta - A->zeta) + (real)(B->cy - A->cy) * S->w;
        gap = R_(sqrt)(rx * rx + ry * ry) - S->d; k = 0;
    } else if (e->type == T_DIV || e->type == T_PISTON) {     /* M2: the face gap of the body's side of the event */
        const int o = e->type == T_DIV ? e->b >> 1 : OBJ_PL + e->b; const Obj3* O = &S->obj[o];
        real c, vo; obj_at(O, S->now, &c, &vo);
        const real xp = (real)A->cx * S->w + A->xi;
        const real sg = o < EDMD_MAX_DIVIDERS ? ((e->b & 1) ? -1.0 : 1.0) : (o == OBJ_PL ? -1.0 : 1.0);
        gap = sg * (c - xp) - O->h;
        k = e->type == T_DIV ? 0 : 1;
        if (R_(fabs)(gap) > S->contact_max_obj[k]) S->contact_max_obj[k] = R_(fabs)(gap);
        S->contact_events++;
        return;
    } else {
        const real x = (real)A->cx * S->w + A->xi, y = (real)A->cy * S->w + A->zeta;
        switch (e->b) { case 0: gap = x - S->R; break; case 1: gap = (S->boxW - S->R) - x; break;
                        case 2: gap = y - S->R; break; default: gap = (S->boxH - S->R) - y; }
        k = 1;
    }
    if (R_(fabs)(gap) > S->contact_max[k]) S->contact_max[k] = R_(fabs)(gap);
    S->contact_events++;
}
static void exec_pair(EDMD3* S, const Ev3* e){
    const int i = e->a, j = e->b;
    advance(S, i); advance(S, j);
    if (S->contact_audit) contact_audit(S, e);
    Disk3 *A = &S->D[i], *B = &S->D[j];
    real dx = (B->xi - A->xi) + (real)(B->cx - A->cx) * S->w, dy = (B->zeta - A->zeta) + (real)(B->cy - A->cy) * S->w;
    real dist = R_(sqrt)(dx * dx + dy * dy);
    if (dist <= 0.0) { dx = S->R; dy = 0.0; dist = S->R; }
    const real nx = dx / dist, ny = dy / dist;
    const real dvn = (B->vx - A->vx) * nx + (B->vy - A->vy) * ny;
    S->virial_accum += (-dvn) * dist; S->virial_pair_events++;     /* as resolve_ab in edmd.c */
    A->vx += dvn * nx; A->vy += dvn * ny;
    B->vx -= dvn * nx; B->vy -= dvn * ny;
    /* M2 ledger scales (read-only here): the same product dvn*n goes to both disks, so momentum changes only by the two
       roundings per axis; the energy by those and by |n|^2 - 1 (a few u, times dvn^2 <= 2 (vA^2 + vB^2)) */
    S->ledSP[0] += R_(fabs)(A->vx) + R_(fabs)(B->vx); S->ledSP[1] += R_(fabs)(A->vy) + R_(fabs)(B->vy);
    S->ledSE += 4.0 * (A->vx * A->vx + A->vy * A->vy + B->vx * B->vx + B->vy * B->vy);
    A->cnt++; B->cnt++; A->last = j; B->last = i;
    S->H.ev_pair++;
    predict_all(S, i); predict_all(S, j);
    local_check(S, i); local_check(S, j);    /* after the predictions: a repair bumps cnt and so retires them */
}
static void exec_wall(EDMD3* S, const Ev3* e){
    const int i = e->a, s = e->b;
    advance(S, i);
    if (S->contact_audit) contact_audit(S, e);
    Disk3* A = &S->D[i];
    if (s < 2) { const real v0 = A->vx; A->vx = -A->vx; S->wall_impulse[s] += R_(fabs)(A->vx - v0); S->ledJw[s] += A->vx - v0; S->ledJ[0] += A->vx - v0; S->ledSP[0] += R_(fabs)(S->ledJ[0]);
                 if (S->evlog) evlog_row(S, WALL_NAME[s], 0.0, v0, A->vx, 0.0); }
    else       { const real v0 = A->vy; A->vy = -A->vy; S->wall_impulse[s] += R_(fabs)(A->vy - v0); S->ledJw[s] += A->vy - v0; S->ledJ[1] += A->vy - v0; S->ledSP[1] += R_(fabs)(S->ledJ[1]);
                 if (S->evlog) evlog_row(S, WALL_NAME[s], 0.0, v0, A->vy, 0.0); }
    S->wall_events[s]++;
    A->cnt++; A->last = -1;
    S->H.ev_wall++;
    predict_all(S, i);
    local_check(S, i);
}
static void exec_cross(EDMD3* S, const Ev3* e){
    const int i = e->a; Disk3* A = &S->D[i];
    advance(S, i);
    const int ocx = A->cx, ocy = A->cy;
    real res;
    cell_remove(S, i);
    switch (e->b) {
        case 0:  res = A->xi - S->w;   A->xi -= S->w;   A->cx++; break;
        case 1:  res = A->xi;          A->xi += S->w;   A->cx--; break;
        case 2:  res = A->zeta - S->w; A->zeta -= S->w; A->cy++; break;
        default: res = A->zeta;        A->zeta += S->w; A->cy--; break;
    }
    if (R_(fabs)(res) > S->H.cross_residual_max) S->H.cross_residual_max = R_(fabs)(res);
    if (A->cx < 0 || A->cy < 0 || A->cx >= S->gw || A->cy >= S->gh) {     /* never: cross_rule does not leave the grid */
        S->H.grid_escape++;
        A->cx = ocx; A->cy = ocy;
        if (e->b == 0) A->xi += S->w; else if (e->b == 1) A->xi -= S->w; else if (e->b == 2) A->zeta += S->w; else A->zeta -= S->w;
        cell_insert(S, i);
        predict_cross(S, i);
        return;
    }
    if (!cell_insert(S, i)) { S->fatal = 1; snprintf(S->fatal_msg, sizeof S->fatal_msg, "cell overflow at crossing of disk %d", i); return; }
    S->H.ev_cross++;
    switch (e->b) {   /* the three cells that just became neighbours */
        case 0:  for (int k = -1; k <= 1; ++k) predict_cell(S, i, A->cx + 1, A->cy + k); break;
        case 1:  for (int k = -1; k <= 1; ++k) predict_cell(S, i, A->cx - 1, A->cy + k); break;
        case 2:  for (int k = -1; k <= 1; ++k) predict_cell(S, i, A->cx + k, A->cy + 1); break;
        default: for (int k = -1; k <= 1; ++k) predict_cell(S, i, A->cx + k, A->cy - 1); break;
    }
    for (int s = 0; s < 4; ++s)
        if (wall_candidate(S, A->cx, A->cy, s) && !wall_candidate(S, ocx, ocy, s)) predict_wall(S, i, s);
    if (A->cx != ocx)                         /* M2: entering a body's band */
        for (int q = 0; q < S->nobj; ++q) {
            const Obj3* O = &S->obj[S->objs[q]];
            if (in_band(O, A->cx) && !in_band(O, ocx)) predict_obj(S, i, S->objs[q], A->xi);
        }
    predict_cross(S, i);
    local_check(S, i);
}

/* M2: tolerances from the time resolution (amendment b); called whenever the energy bound or a body changes */
static void tol_update(EDMD3* S){
    real vdrv = 0.0, mmin = 1.0;
    for (int q = 0; q < S->nobj; ++q) {
        const Obj3* O = &S->obj[S->objs[q]];
        if (O->M > 0.0) { if (O->M < mmin) mmin = O->M; }
        else if (R_(fabs)(O->v) > vdrv) vdrv = R_(fabs)(O->v);
    }
    S->m_min = mmin;
    S->v_ref = R_(sqrt)(2.0 * S->E_bound * (1.0 + 1.0 / mmin)) + vdrv;
    S->c_tol = EDMD3_TOL_K * 2.0 * S->d * S->v_ref * S->u_time;
    S->tol_face = EDMD3_TOL_K * S->v_ref * S->u_time + 8.0 * (R_(nextafter)(S->boxW, INFINITY) - S->boxW);
}
/* E_bound = the largest mechanical energy so far (E0 plus the work done on the system), so v_ref bounds every speed */
static void energy_bound_update(EDMD3* S){
    const real e = S->E0 + S->ledW;
    if (e > S->E_bound) { S->E_bound = e; tol_update(S); }
}

/* M2: a disk meets a divider face or a piston (edmd.c's resolve_wall for EV_DL/DR/PL/PR, unit disk mass) */
static void exec_obj(EDMD3* S, const Ev3* e){
    const int i = e->a, o = e->type == T_DIV ? e->b >> 1 : OBJ_PL + e->b, f = e->type == T_DIV ? (e->b & 1) : 0;
    advance(S, i); obj_advance(S, o);
    if (S->contact_audit) contact_audit(S, e);
    Disk3* A = &S->D[i]; Obj3* O = &S->obj[o];
    const real u1 = A->vx, u2 = O->v;
    real v1, v2;
    if (O->M > 0.0) {
        v1 = ((1.0 - O->M) * u1 + 2.0 * O->M * u2) / (1.0 + O->M);
        v2 = ((O->M - 1.0) * u2 + 2.0 * u1) / (1.0 + O->M);
        O->work += 0.5 * O->M * (v2 * v2 - u2 * u2);
        S->ledSP[0] += 4.0 * (R_(fabs)(u1) + R_(fabs)(v1)) + 4.0 * O->M * (R_(fabs)(u2) + R_(fabs)(v2));
        S->ledSE += 4.0 * (u1 * u1 + v1 * v1) + 4.0 * O->M * (u2 * u2 + v2 * v2);
    } else {
        v1 = 2.0 * u2 - u1; v2 = u2;
        const real dE = 0.5 * (v1 * v1 - u1 * u1), imp = v1 - u1;
        O->work += dE; O->J_inf += imp;
        S->ledJ[0] += imp; S->ledW += dE;
        S->ledSP[0] += R_(fabs)(v1) + R_(fabs)(imp) + R_(fabs)(S->ledJ[0]);
        S->ledSE += v1 * v1 + u1 * u1 + R_(fabs)(S->ledW);
        if (u2 != 0.0) energy_bound_update(S);
    }
    O->Jf[f] -= v1 - u1; O->nf[f]++;
    if (S->evlog) {                                      /* M3: edmd.c's row; dE as edmd.c books it (the body's KE change if M > 0) */
        char kb[8];
        if (e->type == T_DIV) snprintf(kb, sizeof kb, "D%d", o); else snprintf(kb, sizeof kb, "%s", e->b == 0 ? "PL" : "PR");
        evlog_row(S, kb, u2, u1, v1, O->M > 0.0 ? 0.5 * O->M * (v2 * v2 - u2 * u2) : 0.5 * (v1 * v1 - u1 * u1));
    }
    A->vx = v1; A->cnt++; A->last = -2 - o; O->last = i;
    if (e->type == T_DIV) S->H.ev_div++; else S->H.ev_piston++;
    if (v2 != u2) { O->v = v2; obj_changed(S, o, i); }
    predict_all(S, i);
    local_check(S, i);
}
/* M2: a band expires: the next band; the disks of its new columns are predicted */
static void exec_band(EDMD3* S, const Ev3* e){
    const int o = e->a; Obj3* O = &S->obj[o];
    const int oblo = O->blo, obhi = O->bhi;
    band_compute(S, o);
    for (int cx = O->blo; cx <= O->bhi; ++cx) if (cx < oblo || cx > obhi) predict_column(S, o, cx, -1);
    S->H.ev_band++;
}

/* ------------------------------------------------------------------ origin, synchronisation */

static void sync_all_internal(EDMD3* S){
    for (int i = 0; i < S->N; ++i) advance(S, i);
    for (int q = 0; q < S->nobj; ++q) obj_advance(S, S->objs[q]);     /* M2 */
    S->H.syncs++;
    full_check(S);
}
static void origin_shift(EDMD3* S){
    sync_all_internal(S);
    const real s = EDMD3_ORIGIN_SHIFT;
    S->now -= s; S->T0 += s; S->next_check -= s;
    for (int i = 0; i < S->N; ++i) S->D[i].tau -= s;
    for (long k = 0; k < S->heap.n; ++k) S->heap.d[k].t -= s;
    for (int q = 0; q < S->nobj; ++q) { Obj3* O = &S->obj[S->objs[q]]; O->tau -= s; O->t_band -= s; }   /* M2: >= s too; INFINITY stays */
    S->H.origin_shifts++;
}

/* ------------------------------------------------------------------ schedule audit (read-only) */

static void audit_print(EDMD3* S, const char* kind, const char* what, int a, int b, real tbf, real theap){
#ifdef EDMD3_AUDIT_HOOK
    /* ##CHRIS 2026-10-09 (M3, amendment a): a white-box test that includes this file sees every finding at the moment it is
       found (tests/gen3_band_edge_ties.c); compiled out otherwise */
    EDMD3_AUDIT_HOOK(S, kind, what, a, b, tbf, theap);
#endif
    if (S->aprinted >= 50) return;
    S->aprinted++;
    printf("[EDMD3-AUDIT] %s %s a=%d b=%d now=%.17g t_bruteforce=%.17g t_heap=%.17g\n", kind, what, a, b, (double)S->now, (double)tbf, (double)theap);
}
static inline uint64_t akey(int a, int b){ return ((uint64_t)(uint32_t)a << 32) | (uint32_t)b; }
static long atable_slot(const EDMD3* S, uint64_t k){
    uint64_t h = k * 0x9E3779B97F4A7C15ULL; long m = S->ahcap - 1, p = (long)(h >> 20) & m;
    while (S->ahk[p] != UINT64_MAX && S->ahk[p] != k) p = (p + 1) & m;
    return p;
}
static void audit_cmp(EDMD3* S, int cls, real tbf, real th, long* cmp, long* dtc, long* dtrel, const char* what, int a, int b){
    (*cmp)++;
    const real d = R_(fabs)(tbf - th);
    if (d > S->A.max_dt) S->A.max_dt = d;
    {   /* amendment e: every matched event, no floor; horizon = max(t_bruteforce, t_heap) - now (>= d, so finite) */
        const real hz = R_(fmax)(tbf, th) - S->now, rel = hz > 0.0 ? d / hz : 0.0;
        if (d > S->A.cls_max_dt[cls]) { S->A.cls_max_dt[cls] = d; S->A.cls_dt_hz[cls] = hz; }
        if (rel > S->A.cls_max_rel[cls]) { S->A.cls_max_rel[cls] = rel; S->A.cls_rel_hz[cls] = hz; }
    }
    if (d > 1e-9) {
        const real hz = tbf - S->now, rel = hz > 0.0 ? d / hz : INFINITY;
        (*dtc)++;
        if (rel > S->A.max_rel) S->A.max_rel = rel;
        if (rel > 1e-10) { (*dtrel)++; audit_print(S, "dt_rel", what, a, b, tbf, th); }
    }
}
/* M2: the needed position bound of body o from now to its band's expiry (closed form; turning points of a spring) */
static void body_reach(const EDMD3* S, int o, real* lo, real* hi){
    const Obj3* O = &S->obj[o];
    real x, v; obj_at(O, S->now, &x, &v);
    *lo = x; *hi = x;
    if (!isfinite(O->t_band)) return;                  /* a static band: held, or a spring at rest */
    real x1, v1; obj_at(O, O->t_band, &x1, &v1);
    *lo = R_(fmin)(*lo, x1); *hi = R_(fmax)(*hi, x1);
    if (O->harmonic) {                                  /* v(s) = 0 at omega s = atan2(V, D omega) + n pi */
        const real D = x - O->xeq, s0 = R_(atan2)(v, D * O->omega), span = O->t_band - S->now;
        for (int n = -1; n <= 64; ++n) {
            const real s = (s0 + n * G3_PI) / O->omega;
            if (s <= 0.0) continue;
            if (s >= span) break;
            real xt, vt; obj_at(O, S->now + s, &xt, &vt);
            *lo = R_(fmin)(*lo, xt); *hi = R_(fmax)(*hi, xt);
        }
    }
}
/* M2: the live DIV and PISTON events per (disk, body) against body_rule from the absolute state; the BAND events */
static void schedule_audit_bodies(EDMD3* S){
    const int N = S->N, nb = S->nobj; EDMD3_Audit* A = &S->A;
    if (!nb) return;
    if (!S->aot) { S->aot = (real*)malloc((size_t)N * NOBJ * sizeof(real)); S->aob = (int*)malloc((size_t)N * NOBJ * sizeof(int)); }
    for (long k = 0; k < (long)N * nb; ++k) { S->aot[k] = NAN; S->aob[k] = -1; }
    int qof[NOBJ]; long nband[NOBJ]; real tband[NOBJ];
    for (int q = 0; q < nb; ++q) { qof[S->objs[q]] = q; nband[q] = 0; tband[q] = NAN; }
    for (long q = 0; q < S->heap.n; ++q) {
        const Ev3* e = &S->heap.d[q];
        if (e->type != T_BAND && e->type != T_DIV && e->type != T_PISTON) continue;
        if (!ev_live(S, e)) continue;
        if (e->type == T_BAND) { const int k = qof[e->a]; nband[k]++; tband[k] = e->t; continue; }
        const int o = e->type == T_DIV ? e->b >> 1 : OBJ_PL + e->b;
        const long s = (long)e->a * nb + qof[o];
        if (isnan(S->aot[s])) { S->aot[s] = e->t; S->aob[s] = e->b; }
        else {
            if (S->aob[s] != e->b || R_(fabs)(S->aot[s] - e->t) > 1e-9) A->dup_disagree++;
            if (e->t < S->aot[s]) { S->aot[s] = e->t; S->aob[s] = e->b; }
        }
    }
    for (int q = 0; q < nb; ++q) {
        const int o = S->objs[q]; const Obj3* O = &S->obj[o];
        const int isdiv = o < EDMD_MAX_DIVIDERS;
        const char* nm = isdiv ? "DIV" : "PISTON";
        if (isfinite(O->t_band)) {
            if (nband[q] == 0 || tband[q] != O->t_band) { A->band_missing++; audit_print(S, "missing", "BAND", o, 0, O->t_band, tband[q]); }
            if (nband[q] > 1) { A->band_extra++; audit_print(S, "extra", "BAND", o, 0, O->t_band, tband[q]); }
        } else if (nband[q] > 0) { A->band_extra++; audit_print(S, "extra", "BAND", o, 0, NAN, tband[q]); }
        {
            real lo, hi; int nlo, nhi; body_reach(S, o, &lo, &hi); band_columns(S, o, lo, hi, &nlo, &nhi);
            if (nlo <= nhi && (nlo < O->blo || nhi > O->bhi)) { A->band_short++; audit_print(S, "short", "BAND", o, nlo * 10000 + nhi, lo, hi); }
        }
        real c, vo; obj_at(O, S->now, &c, &vo);
        long *cmp = isdiv ? &A->div_cmp : &A->pis_cmp, *dtc = isdiv ? &A->div_dt : &A->pis_dt, *dtr = isdiv ? &A->div_dt_rel : &A->pis_dt_rel;
        long *mis = isdiv ? &A->div_missing : &A->pis_missing, *ext = isdiv ? &A->div_extra : &A->pis_extra;
        long *dfr = isdiv ? &A->div_deferred : &A->pis_deferred, *dfe = isdiv ? &A->div_deferred_early : &A->pis_deferred_early;
        for (int i = 0; i < N; ++i) {
            const Disk3* D = &S->D[i];
            const real xp = (real)D->cx * S->w + S->alx[i];
            const int mlast = D->last == -2 - o && O->last == i;
            real dt = 0.0, g0 = 0.0;
            const int rc = body_rule(O, o, c, vo, xp, D->vx, mlast, &dt, &g0);
            const int bx = isdiv ? 2 * o + (body_side(o, c, xp) > 0 ? 0 : 1) : o - OBJ_PL;
            const long s = (long)i * nb + q; const real th = S->aot[s]; const int hb = S->aob[s];
            if (rc) {
                const real tbf = S->now + dt;
                if (!isnan(th) && hb == bx) audit_cmp(S, isdiv ? EDMD3_CLS_DIV : EDMD3_CLS_PISTON, tbf, th, cmp, dtc, dtr, nm, i, bx);
                else if (!isnan(th)) { (*mis)++; (*ext)++; audit_print(S, "wrong-face", nm, i, bx, tbf, th); }
                else if (in_band(O, D->cx)) { (*mis)++; audit_print(S, "missing", nm, i, bx, tbf, NAN); }
                else {
                    (*dfr)++;
                    const real tc = R_(fmin)(isnan(S->acr[i]) ? INFINITY : S->acr[i], O->t_band);
                    if (tbf < tc - 1e-9 * R_(fmax)(1.0, tc - S->now)) { (*dfe)++; audit_print(S, "deferred-early", nm, i, bx, tbf, tc); }
                }
            } else if (!isnan(th)) { (*ext)++; audit_print(S, "extra", nm, i, hb, NAN, th); }
        }
    }
}
static void schedule_audit(EDMD3* S){
    const int N = S->N; EDMD3_Audit* A = &S->A;
    if (!S->alx) {
        S->alx = (real*)malloc((size_t)N * sizeof(real)); S->aly = (real*)malloc((size_t)N * sizeof(real));
        S->acr = (real*)malloc((size_t)N * sizeof(real)); S->acd = (int*)malloc((size_t)N * sizeof(int));
        S->awl = (real*)malloc((size_t)4 * N * sizeof(real));
    }
    for (int i = 0; i < N; ++i) {
        local_now(S, i, &S->alx[i], &S->aly[i]);
        if (S->alx[i] < -S->tol_cell || S->alx[i] > S->w + S->tol_cell || S->aly[i] < -S->tol_cell || S->aly[i] > S->w + S->tol_cell)
            A->cell_inconsistent++;
        S->acr[i] = NAN; S->acd[i] = -1;
    }
    for (long k = 0; k < 4L * N; ++k) S->awl[k] = NAN;
    long npair = 0;
    for (long q = 0; q < S->heap.n; ++q) if (S->heap.d[q].type == T_PAIR && ev_live(S, &S->heap.d[q])) npair++;
    long cap = 1024; while (cap < 4 * npair + 16) cap <<= 1;
    if (cap > S->ahcap) {
        free(S->ahk); free(S->aht); free(S->ahm);
        S->ahk = (uint64_t*)malloc((size_t)cap * sizeof(uint64_t)); S->aht = (real*)malloc((size_t)cap * sizeof(real));
        S->ahm = (unsigned char*)malloc((size_t)cap); S->ahcap = cap;
    }
    for (long p = 0; p < S->ahcap; ++p) S->ahk[p] = UINT64_MAX;
    /* the live heap */
    for (long q = 0; q < S->heap.n; ++q) {
        const Ev3* e = &S->heap.d[q];
        if (!ev_live(S, e)) continue;
        real* slot = NULL;
        if (e->type == T_CROSS) {   /* a disk has at most ONE live crossing: a second one would be executed twice */
            slot = &S->acr[e->a];
            if (isnan(*slot)) S->acd[e->a] = e->b; else { A->cross_dup++; audit_print(S, "duplicate", "CROSS", e->a, e->b, NAN, e->t); }
        }
        else if (e->type == T_WALL) slot = &S->awl[4L * e->a + e->b];
        else if (e->type == T_PAIR) {
            const long p = atable_slot(S, akey(e->a, e->b));
            if (S->ahk[p] == UINT64_MAX) { S->ahk[p] = akey(e->a, e->b); S->aht[p] = NAN; S->ahm[p] = 0; }
            slot = &S->aht[p];
        }
        if (!slot) continue;
        if (isnan(*slot)) *slot = e->t;
        else { if (R_(fabs)(*slot - e->t) > 1e-9) A->dup_disagree++; if (e->t < *slot) *slot = e->t; }
    }
    /* crossings, from the local positions at now */
    for (int i = 0; i < N; ++i) {
        const Disk3* D = &S->D[i]; real dt;
        const int dir = cross_rule(S, D, S->alx[i], S->aly[i], &dt);
        if (dir < 0) { if (!isnan(S->acr[i])) { A->cross_extra++; audit_print(S, "extra", "CROSS", i, S->acd[i], NAN, S->acr[i]); } continue; }
        if (isnan(S->acr[i])) { A->cross_missing++; audit_print(S, "missing", "CROSS", i, dir, S->now + dt, NAN); continue; }
        if (S->acd[i] != dir) { A->cross_missing++; audit_print(S, "wrong-direction", "CROSS", i, dir, S->now + dt, S->acr[i]); }
        audit_cmp(S, EDMD3_CLS_CROSS, S->now + dt, S->acr[i], &A->cross_cmp, &A->cross_dt, &A->cross_dt_rel, "CROSS", i, dir);
    }
    /* walls and pairs, brute force from ABSOLUTE positions (an independent computation) */
    for (int i = 0; i < N; ++i) {
        const Disk3* D = &S->D[i];
        const real x = (real)D->cx * S->w + S->alx[i], y = (real)D->cy * S->w + S->aly[i];
        for (int s = 0; s < 4; ++s) {
            real gap, speed;
            switch (s) {
                case 0: if (D->vx >= 0.0) continue; gap = x - S->R; speed = -D->vx; break;
                case 1: if (D->vx <= 0.0) continue; gap = (S->boxW - S->R) - x; speed = D->vx; break;
                case 2: if (D->vy >= 0.0) continue; gap = y - S->R; speed = -D->vy; break;
                default: if (D->vy <= 0.0) continue; gap = (S->boxH - S->R) - y; speed = D->vy; break;
            }
            const real tbf = S->now + (gap <= 0.0 ? 0.0 : gap / speed);
            const real th = S->awl[4L * i + s];
            if (isnan(th)) {
                if (wall_candidate(S, D->cx, D->cy, s)) { A->wall_missing++; audit_print(S, "missing", "WALL", i, s, tbf, NAN); }
                else {
                    A->wall_deferred++;
                    const real tc = isnan(S->acr[i]) ? INFINITY : S->acr[i];
                    if (tbf < tc - 1e-9 * R_(fmax)(1.0, tc - S->now)) { A->wall_deferred_early++; audit_print(S, "deferred-early", "WALL", i, s, tbf, tc); }
                }
            } else audit_cmp(S, EDMD3_CLS_WALL, tbf, th, &A->wall_cmp, &A->wall_dt, &A->wall_dt_rel, "WALL", i, s);
            S->awl[4L * i + s] = NAN;     /* consumed: what is left afterwards is extra */
        }
    }
    for (long k = 0; k < 4L * N; ++k) if (!isnan(S->awl[k])) { A->wall_extra++; audit_print(S, "extra", "WALL", (int)(k / 4), (int)(k % 4), NAN, S->awl[k]); }
    for (int i = 0; i < N; ++i) {
        const Disk3* Di = &S->D[i];
        const real xi = (real)Di->cx * S->w + S->alx[i], yi = (real)Di->cy * S->w + S->aly[i];
        for (int j = i + 1; j < N; ++j) {
            if (mutual_last(S, i, j)) continue;
            const Disk3* Dj = &S->D[j];
            const real rx = ((real)Dj->cx * S->w + S->alx[j]) - xi, ry = ((real)Dj->cy * S->w + S->aly[j]) - yi;
            real dt = 0.0, c = 0.0;
            const int rc = pair_rule(rx, ry, Dj->vx - Di->vx, Dj->vy - Di->vy, S->d2, &dt, &c);
            if (!rc) continue;              /* no brute-force event; a heap entry left unmatched is counted as extra below */
            const long p = atable_slot(S, akey(i, j));
            const int inheap = S->ahk[p] != UINT64_MAX;
            const real tbf = S->now + dt;
            if (!inheap) {
                const int nb = abs(Di->cx - Dj->cx) <= 1 && abs(Di->cy - Dj->cy) <= 1;
                if (nb) { A->pair_missing++; audit_print(S, "missing", "PAIR", i, j, tbf, NAN); }
                else {
                    A->pair_deferred++;
                    const real ci = isnan(S->acr[i]) ? INFINITY : S->acr[i], cj = isnan(S->acr[j]) ? INFINITY : S->acr[j];
                    const real tc = ci < cj ? ci : cj;
                    if (tbf < tc - 1e-9 * R_(fmax)(1.0, tc - S->now)) { A->pair_deferred_early++; audit_print(S, "deferred-early", "PAIR", i, j, tbf, tc); }
                }
                continue;
            }
            S->ahm[p] = 1;
            audit_cmp(S, EDMD3_CLS_PAIR, tbf, S->aht[p], &A->pair_cmp, &A->pair_dt, &A->pair_dt_rel, "PAIR", i, j);
        }
    }
    for (long p = 0; p < S->ahcap; ++p)
        if (S->ahk[p] != UINT64_MAX && !S->ahm[p]) {
            A->pair_extra++;
            audit_print(S, "extra", "PAIR", (int)(S->ahk[p] >> 32), (int)(S->ahk[p] & 0xffffffffu), NAN, S->aht[p]);
        }
    schedule_audit_bodies(S);     /* M2 */
    A->audits++;
}

/* ------------------------------------------------------------------ public API */

EDMD3* edmd3_create(const EDMD_Params* prm, double cell_px, char* err, size_t errlen){
    const real w = cell_px > 0.0 ? cell_px : EDMD3_DEFAULT_CELL_PX;
#define FAIL(...) do { if (err && errlen) snprintf(err, errlen, __VA_ARGS__); return NULL; } while (0)
    if (!prm || prm->N <= 0) FAIL("edmd_gen3: no disks");
    if (!(prm->radius > 0.0)) FAIL("edmd_gen3: radius must be > 0");
    if (!(w >= 2.0 * prm->radius)) FAIL("edmd_gen3: cell width %.17g px < diameter %.17g px (sec. 4.7.1, a)", (double)w, 2.0 * prm->radius);
    if (w != R_(floor)(w) || w > 1048576.0) FAIL("edmd_gen3: cell width %.17g px must be an integer number of px (cx*w exact)", (double)w);
    if (!(prm->boxW >= 2.0 * prm->radius) || !(prm->boxH >= 2.0 * prm->radius)) FAIL("edmd_gen3: box smaller than a disk");
    if (prm->heatbath_enabled) FAIL("edmd_gen3: the outer-wall heat bath is not implemented");
    if (prm->species) FAIL("edmd_gen3: species and semipermeable gates are not implemented");
    if (!prm->pp_collisions_enabled) FAIL("edmd_gen3: pair collisions cannot be switched off");
    /* M2: dividers and pistons, read as edmd.c reads them */
    const int ndiv = prm->divider_count < 0 ? 0 : (prm->divider_count > EDMD_MAX_DIVIDERS ? EDMD_MAX_DIVIDERS : prm->divider_count);
    for (int d = 0; d < ndiv; ++d) {
        if (!(prm->divider_thickness[d] > 0.0)) FAIL("edmd_gen3: divider %d has thickness %.17g; edmd.c ignores such a divider, gen3 refuses it", d, prm->divider_thickness[d]);
        if (!(prm->divider_mass[d] >= 0.0) || !isfinite(prm->divider_vx[d]) || !isfinite(prm->divider_x[d])) FAIL("edmd_gen3: divider %d: bad mass, position or velocity", d);
        if (prm->divider_gate_mode[d] != 0) FAIL("edmd_gen3: divider %d: semipermeable gates are not implemented", d);
    }
    if (prm->has_pistonL && !(prm->pistonL_mass >= 0.0 && isfinite(prm->pistonL_x) && isfinite(prm->pistonL_vx))) FAIL("edmd_gen3: left piston: bad mass, position or velocity");
    if (prm->has_pistonR && !(prm->pistonR_mass >= 0.0 && isfinite(prm->pistonR_x) && isfinite(prm->pistonR_vx))) FAIL("edmd_gen3: right piston: bad mass, position or velocity");
#undef FAIL
    EDMD3* S = (EDMD3*)calloc(1, sizeof(EDMD3));
    if (!S) return NULL;
    S->prm = *prm; S->N = prm->N; S->R = prm->radius; S->d = 2.0 * prm->radius; S->d2 = S->d * S->d;
    S->w = w; S->boxW = prm->boxW; S->boxH = prm->boxH;
    S->gw = (int)R_(ceil)(S->boxW / w); S->gh = (int)R_(ceil)(S->boxH / w);
    if (S->gw < 1) S->gw = 1;
    if (S->gh < 1) S->gh = 1;
    S->ncell = S->gw * S->gh;
    { const real q = w / S->d + 1.0; S->ccap = (int)(q * q) + 4; }   /* > the most non-overlapping disks a cell can hold */
    S->cell = (int*)calloc((size_t)S->ncell * (size_t)S->ccap, sizeof(int));
    S->ccount = (int*)calloc((size_t)S->ncell, sizeof(int));
    S->D = (Disk3*)calloc((size_t)S->N, sizeof(Disk3));
    S->out = (EDMD_Particle*)calloc((size_t)S->N, sizeof(EDMD_Particle));
    S->heap.cap = 32L * S->N + 1024; S->heap.d = (Ev3*)calloc((size_t)S->heap.cap, sizeof(Ev3));
    if (!S->cell || !S->ccount || !S->D || !S->out || !S->heap.d) { edmd3_destroy(S); if (err && errlen) snprintf(err, errlen, "edmd_gen3: out of memory"); return NULL; }
    S->compact_at = 64L * S->N + 4096;
    /* the experiment validator's tolerances (experiment_validation.c), so a finding means the same thing */
    S->tol_pair = R_(fmax)(1e-7, 1e-6 * S->d);
    S->tol_wall = R_(fmax)(1e-6, 1e-6 * R_(fmax)(1.0, S->R));
    S->tol_cell = 1e-9;
    /* ##CHRIS 2026-10-08 (M2, amendment b): c_tol and tol_face come from the time resolution and the energy bound, set
       at edmd3_load by tol_update (M1 had c_tol = 64 ulp(d^2) = 8.2e-12 px^2, below the time-rounding scale) */
    S->u_time = R_(ldexp)((real)1, 13 - (R_MANT_DIG - 1));   /* ##CHRIS (E2): ulp(2^13): 2^-39 (double), 2^-50 (x86 long double) */                         /* ulp(2^13) */
    S->band_margin = S->tol_wall;
    S->check_interval = 24.0;                            /* 1 sigma-time */
    S->same_limit = 5000L > 4L * S->N ? 5000L : 4L * S->N;
    S->hash = 1469598103934665603ULL;
    /* M2: the bodies */
    for (int d = 0; d < ndiv; ++d) {
        Obj3* O = &S->obj[d];
        O->active = 1; O->th = prm->divider_thickness[d]; O->h = 0.5 * O->th + S->R;
        O->x = prm->divider_x[d]; O->v = prm->divider_vx[d]; O->M = prm->divider_mass[d];
        O->k = prm->divider_k[d] > 0.0 ? prm->divider_k[d] : 0.0; O->xeq = prm->divider_xeq[d];
    }
    if (prm->has_pistonL) { Obj3* O = &S->obj[OBJ_PL]; O->active = 1; O->h = S->R; O->x = prm->pistonL_x; O->v = prm->pistonL_vx; O->M = prm->pistonL_mass; }
    if (prm->has_pistonR) { Obj3* O = &S->obj[OBJ_PR]; O->active = 1; O->h = S->R; O->x = prm->pistonR_x; O->v = prm->pistonR_vx; O->M = prm->pistonR_mass; }
    for (int o = 0; o < NOBJ; ++o) {
        Obj3* O = &S->obj[o];
        O->last = -1; O->blo = 1; O->bhi = 0; O->t_band = INFINITY;
        O->harmonic = o < EDMD_MAX_DIVIDERS && O->k > 0.0 && O->M > 0.0;
        O->omega = O->harmonic ? R_(sqrt)(O->k / O->M) : 0.0;
        if (O->active) S->objs[S->nobj++] = o;
    }
    return S;
}

void edmd3_destroy(EDMD3* S){
    if (!S) return;
    free(S->cell); free(S->ccount); free(S->D); free(S->out); free(S->heap.d);
    free(S->alx); free(S->aly); free(S->acr); free(S->acd); free(S->awl); free(S->ahk); free(S->aht); free(S->ahm);
    free(S->aot); free(S->aob);
    free(S);
}

/* M2: mechanical energy and momentum of the bodies of finite mass now (computed, not stored), with the rounding scale of
   the sums (sum of |partial sums| and of the squares rounded); the pending spring impulse since each spring's stamp */
static void mech_now(const EDMD3* S, real P[2], real sP[2], real* E, real* sE, real* Jpend, real Jspr[NOBJ]){
    P[0] = P[1] = sP[0] = sP[1] = 0.0; *E = *sE = 0.0; *Jpend = 0.0;
    for (int i = 0; i < S->N; ++i) {
        const Disk3* A = &S->D[i];
        P[0] += A->vx; P[1] += A->vy; sP[0] += R_(fabs)(P[0]); sP[1] += R_(fabs)(P[1]);
        *E += 0.5 * (A->vx * A->vx + A->vy * A->vy); *sE += R_(fabs)(*E) + A->vx * A->vx + A->vy * A->vy;
    }
    for (int q = 0; q < S->nobj; ++q) {
        const int o = S->objs[q]; const Obj3* O = &S->obj[o];
        if (Jspr) Jspr[o] = O->J_spring;
        if (!(O->M > 0.0)) continue;
        real x, v; obj_at(O, S->now, &x, &v);
        P[0] += O->M * v; sP[0] += R_(fabs)(P[0]) + 2.0 * O->M * R_(fabs)(v);
        *E += 0.5 * O->M * v * v; *sE += R_(fabs)(*E) + 2.0 * O->M * v * v;
        if (O->harmonic) {
            const real A = obj_amp(O, x, v), dp = O->M * (v - O->v);
            *E += 0.5 * O->k * (x - O->xeq) * (x - O->xeq); *sE += R_(fabs)(*E) + 8.0 * O->k * A * A;
            *Jpend += dp; if (Jspr) Jspr[o] += dp;
            sP[0] += 8.0 * O->M * O->omega * A;
        }
    }
}

int edmd3_load(EDMD3* S, const EDMD_Particle* P, double t_abs, char* err, size_t errlen){
    if (!S || !P) return 0;
    S->T0 = EDMD3_ORIGIN_SHIFT * R_(floor)(t_abs / EDMD3_ORIGIN_SHIFT);
    S->now = t_abs - S->T0;
    S->heap.n = 0;
    memset(S->ccount, 0, (size_t)S->ncell * sizeof(int));
    for (int i = 0; i < S->N; ++i) {
        Disk3* A = &S->D[i];
        int cx = (int)R_(floor)(P[i].x / S->w), cy = (int)R_(floor)(P[i].y / S->w);
        cx = cx < 0 ? 0 : (cx >= S->gw ? S->gw - 1 : cx);
        cy = cy < 0 ? 0 : (cy >= S->gh ? S->gh - 1 : cy);
        A->cx = cx; A->cy = cy;
        A->xi = P[i].x - (real)cx * S->w; A->zeta = P[i].y - (real)cy * S->w;
        A->vx = P[i].vx; A->vy = P[i].vy; A->tau = S->now;
        A->cnt = 0; A->last = -1; A->pad = 0;
        if (!cell_insert(S, i)) { if (err && errlen) snprintf(err, errlen, "edmd_gen3: cell (%d,%d) overflows at disk %d (overlapping input?)", cx, cy, i); return 0; }
    }
    /* M2: the bodies' time stamps, the ledgers' starting values, the tolerances (before any prediction uses them) */
    for (int q = 0; q < S->nobj; ++q) S->obj[S->objs[q]].tau = S->now;
    {
        real Jp; mech_now(S, S->P0, S->P0s, &S->E0, &S->E0s, &Jp, NULL);
        S->E_bound = S->E0;
        tol_update(S);
    }
    const long f0 = S->H.full_findings, b0 = S->H.body_findings;
    full_check(S);
    if (S->H.full_findings != f0 || S->H.body_findings != b0) {
        if (err && errlen) snprintf(err, errlen, "edmd_gen3: the input overlaps, leaves the box or meets a body (worst surface gap %.3g px; body findings %ld)",
                                    S->H.full_worst, S->H.body_findings - b0);
        return 0;
    }
    for (int i = 0; i < S->N; ++i) {
        const Disk3* A = &S->D[i];
        predict_cross(S, i);
        for (int s = 0; s < 4; ++s) if (wall_candidate(S, A->cx, A->cy, s)) predict_wall(S, i, s);
        for (int dy = -1; dy <= 1; ++dy) for (int dx = -1; dx <= 1; ++dx) {
            const int cx = A->cx + dx, cy = A->cy + dy;
            if (cx < 0 || cy < 0 || cx >= S->gw || cy >= S->gh) continue;
            const int c = cell_of(S, cx, cy); const int* L = &S->cell[(long)c * S->ccap];
            for (int k = 0; k < S->ccount[c]; ++k) if (L[k] > i) predict_pair(S, i, L[k]);
        }
    }
    for (int q = 0; q < S->nobj; ++q) {          /* M2: each body's band and its disks */
        const int o = S->objs[q];
        band_compute(S, o);
        for (int cx = S->obj[o].blo; cx <= S->obj[o].bhi; ++cx) predict_column(S, o, cx, -1);
    }
    S->next_check = S->now + S->check_interval;
    S->virial_t0_abs = t_abs;
    return 1;
}

double edmd3_advance_to(EDMD3* S, double t_abs){
    if (!S || S->fatal) return S ? S->T0 + S->now : 0.0;
    for (;;) {
        const real tt = t_abs - S->T0;
        if (S->heap.n == 0 || S->heap.d[0].t > tt) { if (tt > S->now) S->now = tt; break; }
        Ev3 e; heap_pop(&S->heap, &e);
        if (!ev_live(S, &e)) { S->H.ev_stale++; continue; }
        if (e.t < S->now) { S->H.past_event++; continue; }
        if (e.t == S->same_t) {
            if (++S->same_n > S->same_limit) {
                S->H.stagnation++; S->fatal = 1;
                snprintf(S->fatal_msg, sizeof S->fatal_msg, "stagnation guard: more than %ld live events at t = %.17g", S->same_limit, (double)(S->T0 + e.t));
                heap_push(S, e.t, e.type, e.a, e.b, e.ca, e.cb);
                break;
            }
        } else { S->same_t = e.t; S->same_n = 1; }
        S->now = e.t;
        { const double th = (double)e.t; fnv(&S->hash, &th, sizeof th); }   /* ##CHRIS (E2): the double value in both builds */ fnv(&S->hash, &e.type, sizeof e.type); fnv(&S->hash, &e.a, sizeof e.a); fnv(&S->hash, &e.b, sizeof e.b);
        if (e.type == T_PAIR) exec_pair(S, &e);
        else if (e.type == T_WALL) exec_wall(S, &e);
        else if (e.type == T_CROSS) exec_cross(S, &e);
        else if (e.type == T_BAND) { exec_band(S, &e); if (S->audit_bodies) schedule_audit(S); }   /* M2 */
        else exec_obj(S, &e);                                /* M2: T_DIV, T_PISTON */
        if (S->fatal) break;
        if (S->heap.n > S->compact_at) heap_compact(S);
        if (S->audit_every > 0 && ++S->audit_count >= S->audit_every) { S->audit_count = 0; schedule_audit(S); }
        if (S->check_interval > 0.0 && S->now >= S->next_check) {
            full_check(S);
            while (S->next_check <= S->now) S->next_check += S->check_interval;
        }
        if (S->now >= EDMD3_ORIGIN_SHIFT) { origin_shift(S); S->same_t -= EDMD3_ORIGIN_SHIFT; }
    }
    return S->T0 + S->now;
}

void edmd3_sync_all(EDMD3* S){ if (S) sync_all_internal(S); }

const EDMD_Particle* edmd3_particles(EDMD3* S){
    for (int i = 0; i < S->N; ++i) {
        const Disk3* A = &S->D[i]; real x, y; local_now(S, i, &x, &y);
        S->out[i].x = (real)A->cx * S->w + x; S->out[i].y = (real)A->cy * S->w + y;
        S->out[i].vx = A->vx; S->out[i].vy = A->vy; S->out[i].coll_count = A->cnt;
    }
    full_check(S);
    return S->out;
}
double edmd3_time(const EDMD3* S){ return S->T0 + S->now; }
int    edmd3_count(const EDMD3* S){ return S->N; }
double edmd3_cell_px(const EDMD3* S){ return S->w; }
int    edmd3_fatal(const EDMD3* S, const char** msg){ if (msg) *msg = S->fatal_msg; return S->fatal; }
const EDMD3_Health* edmd3_health(const EDMD3* S){ return &S->H; }
uint64_t edmd3_event_hash(const EDMD3* S){ return S->hash; }
double edmd3_kinetic_energy(EDMD3* S){
    real ke = 0.0;
    for (int i = 0; i < S->N; ++i) ke += 0.5 * (S->D[i].vx * S->D[i].vx + S->D[i].vy * S->D[i].vy);
    return ke;
}
void   edmd3_reset_virial(EDMD3* S){ S->virial_accum = 0.0; S->virial_pair_events = 0; S->virial_t0_abs = S->T0 + S->now; }
double edmd3_compressibility_Z(EDMD3* S){
    const real t = (S->T0 + S->now) - S->virial_t0_abs, ke = edmd3_kinetic_energy(S);
    if (!(t > 0.0) || !(ke > 0.0)) return NAN;
    return 1.0 + S->virial_accum / (2.0 * ke * t);
}
long   edmd3_virial_pair_events(const EDMD3* S){ return S->virial_pair_events; }
double edmd3_wall_impulse(const EDMD3* S, int wall){ return (wall >= 0 && wall < 4) ? S->wall_impulse[wall] : NAN; }
long   edmd3_wall_events(const EDMD3* S, int wall){ return (wall >= 0 && wall < 4) ? S->wall_events[wall] : 0; }
void   edmd3_set_check_interval(EDMD3* S, double units){ S->check_interval = units; S->next_check = S->now + (units > 0.0 ? units : 0.0); }
void   edmd3_set_event_log(EDMD3* S, FILE* f, double time_scale){ if (!S) return; S->evlog = f; S->evlog_tscale = time_scale > 0.0 ? time_scale : 1.0; }
void   edmd3_set_contact_audit(EDMD3* S, int on){ S->contact_audit = on ? 1 : 0; }
long   edmd3_contact_audit_stats(const EDMD3* S, double max_gap_px[2]){ max_gap_px[0] = S->contact_max[0]; max_gap_px[1] = S->contact_max[1]; return S->contact_events; }
void   edmd3_set_schedule_audit(EDMD3* S, long every){ S->audit_every = every; S->audit_count = 0; }
void   edmd3_set_schedule_audit_bodies(EDMD3* S, int on){ S->audit_bodies = on ? 1 : 0; }
void   edmd3_schedule_audit_now(EDMD3* S){ schedule_audit(S); }
const EDMD3_Audit* edmd3_schedule_audit_stats(const EDMD3* S){ return &S->A; }

/* ##CHRIS 2026-10-09 (stage C, M4; 261012 sec. 4.7.16): read-only. The live events now: how many, how many share their time
   EXACTLY with another live event (the tie count of an initial state), how many are due at once (t == now). */
static int cmp_double_asc(const void* a, const void* b){ const real x = *(const real*)a, y = *(const real*)b; return (x > y) - (x < y); }
long edmd3_tie_stats(const EDMD3* S, long* n_live, long* n_now){
    long n = 0, now = 0, tied = 0;
    real* t = (real*)malloc((size_t)(S->heap.n > 0 ? S->heap.n : 1) * sizeof(real));
    if (!t) { if (n_live) *n_live = -1; if (n_now) *n_now = -1; return -1; }
    for (long q = 0; q < S->heap.n; ++q) {
        const Ev3* e = &S->heap.d[q];
        if (!ev_live(S, e)) continue;
        t[n++] = e->t;
        if (e->t == S->now) now++;
    }
    qsort(t, (size_t)n, sizeof(real), cmp_double_asc);
    for (long k = 0; k < n; ++k) if ((k > 0 && t[k] == t[k - 1]) || (k + 1 < n && t[k] == t[k + 1])) tied++;
    free(t);
    if (n_live) *n_live = n;
    if (n_now) *n_now = now;
    return tied;
}

/* ##CHRIS 2026-10-09 (stage H; 261012 sec. 4.7.12 stage H, sec. 4.7 item 6): CHECKPOINT AND RESTART AT AN EVENT BOUNDARY.
   Between two calls of edmd3_advance_to the engine is at an event boundary: every event with t <= now has run, none is half
   done. edmd3_checkpoint_write stores the WHOLE dynamic state there, byte for byte: the state struct (every counter, time,
   tolerance, body, ledger, audit total, the stagnation guard and the event-hash state), every disk (local position, velocity,
   time stamp, cell, collision counter, slot, last partner), the cell lists in their order, and the heap ARRAY as it stands
   (live and stale events, in heap order). edmd3_checkpoint_read puts it into a state built from the same parameters
   (edmd3_create + edmd3_load of the same run, as a restarting driver builds it); only that state's own buffers, its event-log
   file and its copy of the parameters stay. The restarted engine is then the writer's bytes, so every later event, time,
   hash and counter is the writer's. Why the heap is stored and not rebuilt (the design note said "rebuilt in canonical
   order"): a live event of the writer was predicted at its own time from the states then (lazy invalidation, origin shifts
   in between), so a rebuild from the stored states can differ in the last bits and the trajectories part. The file is for
   the same build: a header with the struct sizes and the geometry, and a mismatch is refused. On a failed read the state
   may be partly overwritten: the caller destroys it. */
typedef struct {
    char   magic[16];                         /* "edmd3-ckpt-v1" */
    long   s_state, s_disk, s_ev, s_real;     /* sizeof EDMD3, Disk3, Ev3, real of the writer */
    long   N, ncell, ccap, heap_n;
    double w, boxW, boxH;
} Ckpt3Head;
#define CKPT3_MAGIC "edmd3-ckpt-v1"

int edmd3_checkpoint_write(const EDMD3* S, FILE* f){
    if (!S || !f) return 0;
    Ckpt3Head h; memset(&h, 0, sizeof h);
    memcpy(h.magic, CKPT3_MAGIC, sizeof CKPT3_MAGIC);
    h.s_state = (long)sizeof(EDMD3); h.s_disk = (long)sizeof(Disk3); h.s_ev = (long)sizeof(Ev3); h.s_real = (long)sizeof(real);
    h.N = S->N; h.ncell = S->ncell; h.ccap = S->ccap; h.heap_n = S->heap.n;
    h.w = (double)S->w; h.boxW = (double)S->boxW; h.boxH = (double)S->boxH;
    const size_t nc = (size_t)S->ncell * (size_t)S->ccap;
    return fwrite(&h, sizeof h, 1, f) == 1 && fwrite(S, sizeof *S, 1, f) == 1
        && fwrite(S->D, sizeof(Disk3), (size_t)S->N, f) == (size_t)S->N
        && fwrite(S->ccount, sizeof(int), (size_t)S->ncell, f) == (size_t)S->ncell
        && fwrite(S->cell, sizeof(int), nc, f) == nc
        && (S->heap.n == 0 || fwrite(S->heap.d, sizeof(Ev3), (size_t)S->heap.n, f) == (size_t)S->heap.n);
}

int edmd3_checkpoint_read(EDMD3* S, FILE* f, char* err, size_t errlen){
#define FAIL(...) do { if (err && errlen) snprintf(err, errlen, __VA_ARGS__); return 0; } while (0)
    if (!S || !f) FAIL("edmd_gen3: checkpoint: no state or no file");
    Ckpt3Head h;
    if (fread(&h, sizeof h, 1, f) != 1 || memcmp(h.magic, CKPT3_MAGIC, sizeof CKPT3_MAGIC) != 0) FAIL("edmd_gen3: not a gen3 checkpoint (v1)");
    if (h.s_state != (long)sizeof(EDMD3) || h.s_disk != (long)sizeof(Disk3) || h.s_ev != (long)sizeof(Ev3) || h.s_real != (long)sizeof(real))
        FAIL("edmd_gen3: checkpoint of another build (struct sizes %ld %ld %ld %ld, here %ld %ld %ld %ld)", h.s_state, h.s_disk, h.s_ev,
             h.s_real, (long)sizeof(EDMD3), (long)sizeof(Disk3), (long)sizeof(Ev3), (long)sizeof(real));
    if (h.N != S->N || h.ncell != S->ncell || h.ccap != S->ccap || h.w != (double)S->w || h.boxW != (double)S->boxW || h.boxH != (double)S->boxH)
        FAIL("edmd_gen3: checkpoint of another geometry (N %ld, %ld cells of capacity %ld, w %.17g, box %.17g x %.17g px)", h.N, h.ncell,
             h.ccap, h.w, h.boxW, h.boxH);
    EDMD3* T = (EDMD3*)malloc(sizeof *T);
    if (!T) FAIL("edmd_gen3: out of memory (checkpoint)");
    if (fread(T, sizeof *T, 1, f) != 1 || T->heap.n != h.heap_n || T->heap.n < 0) { free(T); FAIL("edmd_gen3: checkpoint truncated or inconsistent (state)"); }
    /* this state's own buffers, event log and parameters stay; everything else is the writer's */
    T->prm = S->prm;
    T->cell = S->cell; T->ccount = S->ccount; T->D = S->D; T->out = S->out;
    T->heap.d = S->heap.d; T->heap.cap = S->heap.cap;
    T->alx = S->alx; T->aly = S->aly; T->acr = S->acr; T->awl = S->awl; T->acd = S->acd;
    T->ahcap = S->ahcap; T->ahk = S->ahk; T->aht = S->aht; T->ahm = S->ahm; T->aot = S->aot; T->aob = S->aob;
    T->evlog = S->evlog; T->evlog_tscale = S->evlog_tscale;
    if (T->heap.cap < T->heap.n) {
        Ev3* nd = (Ev3*)realloc(T->heap.d, (size_t)(T->heap.n + 1024) * sizeof(Ev3));
        if (!nd) { free(T); FAIL("edmd_gen3: out of memory (checkpoint heap)"); }
        T->heap.d = nd; T->heap.cap = T->heap.n + 1024; S->heap.d = nd; S->heap.cap = T->heap.cap;
    }
    const size_t nc = (size_t)T->ncell * (size_t)T->ccap;
    if (fread(T->D, sizeof(Disk3), (size_t)T->N, f) != (size_t)T->N || fread(T->ccount, sizeof(int), (size_t)T->ncell, f) != (size_t)T->ncell
        || fread(T->cell, sizeof(int), nc, f) != nc || (T->heap.n > 0 && fread(T->heap.d, sizeof(Ev3), (size_t)T->heap.n, f) != (size_t)T->heap.n)) {
        free(T); FAIL("edmd_gen3: checkpoint truncated (disks, cells or heap)");
    }
    *S = *T;
    free(T);
    return 1;
#undef FAIL
}

/* ------------------------------------------------------------------ M2: bodies, ledgers, tolerances (public) */

long edmd3_contact_audit_stats4(const EDMD3* S, double max_gap_px[4]){
    max_gap_px[0] = S->contact_max[0]; max_gap_px[1] = S->contact_max[1];
    max_gap_px[2] = S->contact_max_obj[0]; max_gap_px[3] = S->contact_max_obj[1];
    return S->contact_events;
}
/* a velocity or mass change of body o at now; the change of its momentum and energy (as a body of finite mass) is an
   external input in the ledgers */
static int set_motion(EDMD3* S, int o, real mass, real vx){
    if (!S || o < 0 || o >= NOBJ || !S->obj[o].active || !(mass >= 0.0) || !isfinite(vx)) return 0;
    Obj3* O = &S->obj[o];
    obj_advance(S, o);
    const real sp0 = O->harmonic ? 0.5 * O->k * (O->x - O->xeq) * (O->x - O->xeq) : 0.0;
    const real p0 = O->M > 0.0 ? O->M * O->v : 0.0, e0 = (O->M > 0.0 ? 0.5 * O->M * O->v * O->v : 0.0) + sp0;
    O->M = mass; O->v = vx;
    O->harmonic = o < EDMD_MAX_DIVIDERS && O->k > 0.0 && O->M > 0.0;
    O->omega = O->harmonic ? R_(sqrt)(O->k / O->M) : 0.0;
    const real sp1 = O->harmonic ? 0.5 * O->k * (O->x - O->xeq) * (O->x - O->xeq) : 0.0;
    const real p1 = O->M > 0.0 ? O->M * O->v : 0.0, e1 = (O->M > 0.0 ? 0.5 * O->M * O->v * O->v : 0.0) + sp1;
    S->ledJ[0] += p1 - p0; S->J_api += p1 - p0; S->ledW += e1 - e0;
    S->ledSP[0] += R_(fabs)(p0) + R_(fabs)(p1) + R_(fabs)(S->ledJ[0]); S->ledSE += R_(fabs)(e0) + R_(fabs)(e1) + R_(fabs)(S->ledW);
    O->last = -1;
    tol_update(S);
    energy_bound_update(S);
    obj_changed(S, o, -1);
    if (S->audit_bodies) schedule_audit(S);
    return 1;
}
int edmd3_set_divider_motion(EDMD3* S, int d, double mass, double vx){ return (d >= 0 && d < EDMD_MAX_DIVIDERS) ? set_motion(S, d, mass, vx) : 0; }
int edmd3_set_piston_motion(EDMD3* S, int side, double mass, double vx){ return (side == 0 || side == 1) ? set_motion(S, OBJ_PL + side, mass, vx) : 0; }
static int body_state(const EDMD3* S, int o, double* x, double* vx){   /* ##CHRIS (E2): the API's doubles */
    if (!S || o < 0 || o >= NOBJ || !S->obj[o].active) return 0;
    real a, b; obj_at(&S->obj[o], S->now, &a, &b);
    if (x) *x = a;
    if (vx) *vx = b;
    return 1;
}
int edmd3_divider_state(const EDMD3* S, int d, double* x, double* vx){ return (d >= 0 && d < EDMD_MAX_DIVIDERS) ? body_state(S, d, x, vx) : 0; }
int edmd3_piston_state(const EDMD3* S, int side, double* x, double* vx){ return (side == 0 || side == 1) ? body_state(S, OBJ_PL + side, x, vx) : 0; }
double edmd3_divider_impulse(const EDMD3* S, int d, int face){ return (d >= 0 && d < EDMD_MAX_DIVIDERS && (face == 0 || face == 1)) ? S->obj[d].Jf[face] : NAN; }
long   edmd3_divider_events(const EDMD3* S, int d, int face){ return (d >= 0 && d < EDMD_MAX_DIVIDERS && (face == 0 || face == 1)) ? S->obj[d].nf[face] : 0; }
double edmd3_piston_impulse(const EDMD3* S, int side){ return (side == 0 || side == 1) ? S->obj[OBJ_PL + side].Jf[0] : NAN; }
long   edmd3_piston_events(const EDMD3* S, int side){ return (side == 0 || side == 1) ? S->obj[OBJ_PL + side].nf[0] : 0; }
double edmd3_work_divider(const EDMD3* S, int d){ return (d >= 0 && d < EDMD_MAX_DIVIDERS) ? S->obj[d].work : NAN; }
double edmd3_work_piston(const EDMD3* S, int side){ return (side == 0 || side == 1) ? S->obj[OBJ_PL + side].work : NAN; }

void edmd3_ledger(EDMD3* S, EDMD3_Ledger* L){
    memset(L, 0, sizeof *L);
    real sP[2], sE, Jpend, Jspr[NOBJ], Pm[2], Em;   /* ##CHRIS (E2): through reals into the API's doubles */
    mech_now(S, Pm, sP, &Em, &sE, &Jpend, Jspr);
    L->P[0] = Pm[0]; L->P[1] = Pm[1]; L->E = Em;
    L->P0[0] = S->P0[0]; L->P0[1] = S->P0[1]; L->E0 = S->E0;
    L->J[0] = S->ledJ[0] + Jpend; L->J[1] = S->ledJ[1];
    L->W = S->ledW;
    for (int a = 0; a < 2; ++a) L->scale_P[a] = U_ROUND * (S->ledSP[a] + sP[a] + S->P0s[a] + R_(fabs)(Jpend));
    L->scale_E = U_ROUND * (S->ledSE + sE + S->E0s);
    for (int s = 0; s < 4; ++s) L->J_wall[s] = S->ledJw[s];
    for (int d = 0; d < EDMD_MAX_DIVIDERS; ++d) { L->J_div[d] = S->obj[d].J_inf; L->J_spring[d] = S->obj[d].active ? Jspr[d] : 0.0; }
    L->J_piston[0] = S->obj[OBJ_PL].J_inf; L->J_piston[1] = S->obj[OBJ_PR].J_inf;
    L->J_api = S->J_api;
}
int edmd3_health_clean(const EDMD3* S){
    const EDMD3_Health* h = &S->H;
    return !S->fatal && h->overlap_repair == 0 && h->wall_overdue == 0 && h->past_event == 0 && h->clamp_repair == 0 &&
           h->cell_repair == 0 && h->grid_escape == 0 && h->stagnation == 0 && h->local_findings == 0 && h->full_findings == 0 &&
           h->obj_overlap_repair == 0 && h->body_findings == 0;
}
void edmd3_tolerances(const EDMD3* S, EDMD3_Tol* T){
    T->u_time = S->u_time; T->E_bound = S->E_bound; T->m_min = S->m_min; T->v_ref = S->v_ref; T->K = EDMD3_TOL_K;
    T->c_tol = S->c_tol; T->tol_face = S->tol_face; T->band_margin = S->band_margin;
    T->tol_pair = S->tol_pair; T->tol_wall = S->tol_wall; T->tol_cell = S->tol_cell;
}
