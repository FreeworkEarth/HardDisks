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
   and checks cannot change a trajectory (the audit gate compares hashes with and without them). */

#include "edmd_gen3.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

enum { T_CROSS = 0, T_BAND = 1, T_WALL = 2, T_DIV = 3, T_PISTON = 4, T_PAIR = 5 };

typedef struct {
    double t;              /* origin-relative event time */
    int    type, a, b;     /* the tie-break key with t */
    int    ca, cb;         /* collision counters at prediction (cb: partner's, PAIR only; else 0) */
    int    pad;            /* always 0: no uninitialised bytes anywhere (sec. 4.7.1, c) */
} Ev3;

typedef struct { Ev3* d; long n, cap; } Heap3;

typedef struct {
    double xi, zeta;       /* position relative to the lower-left corner of the cell [px] */
    double vx, vy;         /* velocity [px per internal unit] */
    double tau;            /* time of (xi, zeta), origin-relative */
    int    cx, cy;         /* cell */
    int    cnt;            /* collision counter */
    int    slot;           /* index in the cell's list */
    int    last;           /* partner of the last velocity change if that was a pair collision, else -1 */
    int    pad;
} Disk3;

struct EDMD3 {
    EDMD_Params prm;
    int    N;
    double R, d, d2, w, boxW, boxH;
    int    gw, gh, ncell, ccap;
    int   *cell, *ccount;            /* ncell * ccap disk indices; occupancy */
    Disk3 *D;
    Heap3  heap;
    long   compact_at;
    double T0, now;
    double check_interval, next_check;
    EDMD_Particle* out;
    EDMD3_Health H;
    double tol_pair, tol_wall, tol_cell, c_tol;
    double virial_accum, virial_t0_abs; long virial_pair_events;
    double wall_impulse[4]; long wall_events[4];
    uint64_t hash;
    double same_t; long same_n, same_limit;
    int    fatal; char fatal_msg[256];
    /* audits */
    int    contact_audit; double contact_max[2]; long contact_events;
    long   audit_every, audit_count; EDMD3_Audit A;
    double *alx, *aly;               /* audit scratch: local position of every disk at now */
    double *acr, *awl;               /* live crossing time per disk (and direction in acd), live wall time per (disk, side) */
    int    *acd;
    long   ahcap; uint64_t* ahk; double* aht; unsigned char* ahm;   /* pair table: key, time, matched */
    long   aprinted;
};

/* ------------------------------------------------------------------ heap, ordered by (t, type, a, b) */

static inline int ev_less(const Ev3* x, const Ev3* y){
    if (x->t != y->t) return x->t < y->t;
    if (x->type != y->type) return x->type < y->type;
    if (x->a != y->a) return x->a < y->a;
    return x->b < y->b;
}
static void heap_push(EDMD3* S, double t, int type, int a, int b, int ca, int cb){
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
    if (S->D[e->a].cnt != e->ca) return 0;
    if (e->type == T_PAIR && S->D[e->b].cnt != e->cb) return 0;
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
        case 0: return (double)cx * S->w <= S->R;
        case 1: return (double)(cx + 1) * S->w >= S->boxW - S->R;
        case 2: return (double)cy * S->w <= S->R;
        default: return (double)(cy + 1) * S->w >= S->boxH - S->R;
    }
}

/* ------------------------------------------------------------------ motion */

static inline void advance(EDMD3* S, int i){
    Disk3* A = &S->D[i];
    const double dt = S->now - A->tau;
    if (dt != 0.0) { A->xi += A->vx * dt; A->zeta += A->vy * dt; A->tau = S->now; }
}
/* local position of disk j at now, not stored */
static inline void local_now(const EDMD3* S, int j, double* x, double* y){
    const Disk3* B = &S->D[j]; const double dt = S->now - B->tau;
    *x = B->xi + B->vx * dt; *y = B->zeta + B->vy * dt;
}
static inline void fnv(uint64_t* h, const void* p, size_t n){
    const unsigned char* c = (const unsigned char*)p;
    for (size_t k = 0; k < n; ++k) { *h ^= c[k]; *h *= 1099511628211ULL; }
}

/* ------------------------------------------------------------------ predictions (disk i is synchronised: tau_i == now) */

/* the pair rule shared with the brute-force audit: 0 none, 1 future contact at dt, 2 at once (c < 0, approaching) */
static inline int pair_rule(double rx, double ry, double vx, double vy, double d2, double* dt, double* c_out){
    const double b = rx * vx + ry * vy;
    if (b >= 0.0) return 0;
    const double vv = vx * vx + vy * vy;
    const double c = rx * rx + ry * ry - d2;
    *c_out = c;
    if (c <= 0.0) { *dt = 0.0; return 2; }
    const double disc = b * b - vv * c;
    if (disc <= 0.0) return 0;
    *dt = (-b - sqrt(disc)) / vv;
    if (*dt < 0.0) *dt = 0.0;
    return 1;
}
static inline int mutual_last(const EDMD3* S, int i, int j){ return S->D[i].last == j && S->D[j].last == i; }

static void predict_pair(EDMD3* S, int i, int j){
    if (mutual_last(S, i, j)) return;
    const Disk3 *A = &S->D[i], *B = &S->D[j];
    double bx, by; local_now(S, j, &bx, &by);
    const double rx = (bx - A->xi) + (double)(B->cx - A->cx) * S->w;
    const double ry = (by - A->zeta) + (double)(B->cy - A->cy) * S->w;
    double dt = 0.0, c = 0.0;
    const int rc = pair_rule(rx, ry, B->vx - A->vx, B->vy - A->vy, S->d2, &dt, &c);
    if (!rc) return;
    if (rc == 2) { if (c < -S->c_tol) S->H.overlap_repair++; else S->H.contact_now++; }
    const int a = i < j ? i : j, b = i < j ? j : i;
    heap_push(S, S->now + dt, T_PAIR, a, b, S->D[a].cnt, S->D[b].cnt);
}
/* wall gap of disk i (synchronised) to wall s, and the closing speed; 0 if moving away or parallel */
static inline int wall_gap(const EDMD3* S, const Disk3* A, double xl, double yl, int s, double* gap, double* speed){
    switch (s) {
        case 0: if (A->vx >= 0.0) return 0; *gap = ((double)A->cx * S->w - S->R) + xl; *speed = -A->vx; return 1;
        case 1: if (A->vx <= 0.0) return 0; *gap = (S->boxW - S->R - (double)A->cx * S->w) - xl; *speed = A->vx; return 1;
        case 2: if (A->vy >= 0.0) return 0; *gap = ((double)A->cy * S->w - S->R) + yl; *speed = -A->vy; return 1;
        default: if (A->vy <= 0.0) return 0; *gap = (S->boxH - S->R - (double)A->cy * S->w) - yl; *speed = A->vy; return 1;
    }
}
static void predict_wall(EDMD3* S, int i, int s){
    const Disk3* A = &S->D[i]; double gap, speed;
    if (!wall_gap(S, A, A->xi, A->zeta, s, &gap, &speed)) return;
    double dt = 0.0;
    if (gap <= 0.0) S->H.wall_overdue++; else dt = gap / speed;
    heap_push(S, S->now + dt, T_WALL, i, s, A->cnt, 0);
}
/* next crossing from local coordinates (x, y) and the velocity: direction and time from now; -1 if none */
static inline int cross_rule(const EDMD3* S, const Disk3* A, double x, double y, double* dt){
    double tx = INFINITY, ty = INFINITY; int dx = -1, dy = -1;
    if (A->vx > 0.0 && A->cx + 1 < S->gw) { tx = (S->w - x) / A->vx; dx = 0; }
    else if (A->vx < 0.0 && A->cx > 0)    { tx = x / (-A->vx); dx = 1; }
    if (A->vy > 0.0 && A->cy + 1 < S->gh) { ty = (S->w - y) / A->vy; dy = 2; }
    else if (A->vy < 0.0 && A->cy > 0)    { ty = y / (-A->vy); dy = 3; }
    if (dx < 0 && dy < 0) return -1;
    int dir; double t;
    if (dy < 0 || (dx >= 0 && tx <= ty)) { dir = dx; t = tx; } else { dir = dy; t = ty; }
    *dt = t > 0.0 ? t : 0.0;
    return dir;
}
static void predict_cross(EDMD3* S, int i){
    const Disk3* A = &S->D[i]; double dt;
    const int dir = cross_rule(S, A, A->xi, A->zeta, &dt);
    if (dir >= 0) heap_push(S, S->now + dt, T_CROSS, i, dir, A->cnt, 0);
}
static void predict_cell(EDMD3* S, int i, int cx, int cy){
    if (cx < 0 || cy < 0 || cx >= S->gw || cy >= S->gh) return;
    const int c = cell_of(S, cx, cy); const int* L = &S->cell[(long)c * S->ccap];
    for (int k = 0; k < S->ccount[c]; ++k) if (L[k] != i) predict_pair(S, i, L[k]);
}
/* after a velocity change of i (or at load): its crossing, its walls, all nine cells */
static void predict_all(EDMD3* S, int i){
    const Disk3* A = &S->D[i];
    predict_cross(S, i);
    for (int s = 0; s < 4; ++s) if (wall_candidate(S, A->cx, A->cy, s)) predict_wall(S, i, s);
    for (int dy = -1; dy <= 1; ++dy) for (int dx = -1; dx <= 1; ++dx) predict_cell(S, i, A->cx + dx, A->cy + dy);
}

/* ------------------------------------------------------------------ validator (sec. 4.7.1, h); read-only except the box repair */

static void reflect_into_box(EDMD3* S, int i);
/* the event's disk i (synchronised) against every disk of its nine cells and the four walls */
static void local_check(EDMD3* S, int i){
    Disk3* A = &S->D[i]; double worst = 0.0; int found = 0;
    for (int dy = -1; dy <= 1; ++dy) for (int dx = -1; dx <= 1; ++dx) {
        const int cx = A->cx + dx, cy = A->cy + dy;
        if (cx < 0 || cy < 0 || cx >= S->gw || cy >= S->gh) continue;
        const int c = cell_of(S, cx, cy); const int* L = &S->cell[(long)c * S->ccap];
        for (int k = 0; k < S->ccount[c]; ++k) {
            const int j = L[k]; if (j == i) continue;
            double bx, by; local_now(S, j, &bx, &by);
            const double rx = (bx - A->xi) + (double)dx * S->w, ry = (by - A->zeta) + (double)dy * S->w;
            const double gap = sqrt(rx * rx + ry * ry) - S->d;
            if (gap < worst) worst = gap;
            if (gap < -S->tol_pair) found = 1;
        }
    }
    const double x = (double)A->cx * S->w + A->xi, y = (double)A->cy * S->w + A->zeta;
    const double g[4] = { x - S->R, (S->boxW - S->R) - x, y - S->R, (S->boxH - S->R) - y };
    int out = 0;
    for (int s = 0; s < 4; ++s) { if (g[s] < worst) worst = g[s]; if (g[s] < -S->tol_wall) out = 1; }
    if (A->xi < -S->tol_cell || A->xi > S->w + S->tol_cell || A->zeta < -S->tol_cell || A->zeta > S->w + S->tol_cell) {
        /* outside its own cell: a bookkeeping failure; re-file it from its absolute position */
        S->H.cell_repair++;
        cell_remove(S, i);
        int ncx = (int)floor(x / S->w), ncy = (int)floor(y / S->w);
        ncx = ncx < 0 ? 0 : (ncx >= S->gw ? S->gw - 1 : ncx);
        ncy = ncy < 0 ? 0 : (ncy >= S->gh ? S->gh - 1 : ncy);
        A->xi = x - (double)ncx * S->w; A->zeta = y - (double)ncy * S->w; A->cx = ncx; A->cy = ncy;
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
    double x = (double)A->cx * S->w + A->xi, y = (double)A->cy * S->w + A->zeta;
    if (x < S->R)           { x = S->R;           if (A->vx < 0.0) A->vx = -A->vx; }
    if (x > S->boxW - S->R) { x = S->boxW - S->R; if (A->vx > 0.0) A->vx = -A->vx; }
    if (y < S->R)           { y = S->R;           if (A->vy < 0.0) A->vy = -A->vy; }
    if (y > S->boxH - S->R) { y = S->boxH - S->R; if (A->vy > 0.0) A->vy = -A->vy; }
    A->xi = x - (double)A->cx * S->w; A->zeta = y - (double)A->cy * S->w;
    S->H.clamp_repair++; A->cnt++; A->last = -1;
    predict_all(S, i);
}
/* every pair within nine cells and every wall, positions computed at now (not stored) */
static void full_check(EDMD3* S){
    double worst = 0.0; long found = 0;
    for (int i = 0; i < S->N; ++i) {
        const Disk3* A = &S->D[i]; double ax, ay; local_now(S, i, &ax, &ay);
        for (int dy = -1; dy <= 1; ++dy) for (int dx = -1; dx <= 1; ++dx) {
            const int cx = A->cx + dx, cy = A->cy + dy;
            if (cx < 0 || cy < 0 || cx >= S->gw || cy >= S->gh) continue;
            const int c = cell_of(S, cx, cy); const int* L = &S->cell[(long)c * S->ccap];
            for (int k = 0; k < S->ccount[c]; ++k) {
                const int j = L[k]; if (j <= i) continue;
                double bx, by; local_now(S, j, &bx, &by);
                const double rx = (bx - ax) + (double)dx * S->w, ry = (by - ay) + (double)dy * S->w;
                const double gap = sqrt(rx * rx + ry * ry) - S->d;
                if (gap < worst) worst = gap;
                if (gap < -S->tol_pair) found++;
            }
        }
        const double x = (double)A->cx * S->w + ax, y = (double)A->cy * S->w + ay;
        const double g[4] = { x - S->R, (S->boxW - S->R) - x, y - S->R, (S->boxH - S->R) - y };
        for (int s = 0; s < 4; ++s) { if (g[s] < worst) worst = g[s]; if (g[s] < -S->tol_wall) found++; }
    }
    S->H.full_checks++; S->H.full_findings += found;
    if (worst < S->H.full_worst) S->H.full_worst = worst;
}

/* ------------------------------------------------------------------ event execution */

static void contact_audit(EDMD3* S, const Ev3* e){
    const Disk3* A = &S->D[e->a]; double gap; int k;
    if (e->type == T_PAIR) {
        const Disk3* B = &S->D[e->b];
        const double rx = (B->xi - A->xi) + (double)(B->cx - A->cx) * S->w, ry = (B->zeta - A->zeta) + (double)(B->cy - A->cy) * S->w;
        gap = sqrt(rx * rx + ry * ry) - S->d; k = 0;
    } else {
        const double x = (double)A->cx * S->w + A->xi, y = (double)A->cy * S->w + A->zeta;
        switch (e->b) { case 0: gap = x - S->R; break; case 1: gap = (S->boxW - S->R) - x; break;
                        case 2: gap = y - S->R; break; default: gap = (S->boxH - S->R) - y; }
        k = 1;
    }
    if (fabs(gap) > S->contact_max[k]) S->contact_max[k] = fabs(gap);
    S->contact_events++;
}
static void exec_pair(EDMD3* S, const Ev3* e){
    const int i = e->a, j = e->b;
    advance(S, i); advance(S, j);
    if (S->contact_audit) contact_audit(S, e);
    Disk3 *A = &S->D[i], *B = &S->D[j];
    double dx = (B->xi - A->xi) + (double)(B->cx - A->cx) * S->w, dy = (B->zeta - A->zeta) + (double)(B->cy - A->cy) * S->w;
    double dist = sqrt(dx * dx + dy * dy);
    if (dist <= 0.0) { dx = S->R; dy = 0.0; dist = S->R; }
    const double nx = dx / dist, ny = dy / dist;
    const double dvn = (B->vx - A->vx) * nx + (B->vy - A->vy) * ny;
    S->virial_accum += (-dvn) * dist; S->virial_pair_events++;     /* as resolve_ab in edmd.c */
    A->vx += dvn * nx; A->vy += dvn * ny;
    B->vx -= dvn * nx; B->vy -= dvn * ny;
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
    if (s < 2) { const double v0 = A->vx; A->vx = -A->vx; S->wall_impulse[s] += fabs(A->vx - v0); }
    else       { const double v0 = A->vy; A->vy = -A->vy; S->wall_impulse[s] += fabs(A->vy - v0); }
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
    double res;
    cell_remove(S, i);
    switch (e->b) {
        case 0:  res = A->xi - S->w;   A->xi -= S->w;   A->cx++; break;
        case 1:  res = A->xi;          A->xi += S->w;   A->cx--; break;
        case 2:  res = A->zeta - S->w; A->zeta -= S->w; A->cy++; break;
        default: res = A->zeta;        A->zeta += S->w; A->cy--; break;
    }
    if (fabs(res) > S->H.cross_residual_max) S->H.cross_residual_max = fabs(res);
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
    predict_cross(S, i);
    local_check(S, i);
}

/* ------------------------------------------------------------------ origin, synchronisation */

static void sync_all_internal(EDMD3* S){
    for (int i = 0; i < S->N; ++i) advance(S, i);
    S->H.syncs++;
    full_check(S);
}
static void origin_shift(EDMD3* S){
    sync_all_internal(S);
    const double s = EDMD3_ORIGIN_SHIFT;
    S->now -= s; S->T0 += s; S->next_check -= s;
    for (int i = 0; i < S->N; ++i) S->D[i].tau -= s;
    for (long k = 0; k < S->heap.n; ++k) S->heap.d[k].t -= s;
    S->H.origin_shifts++;
}

/* ------------------------------------------------------------------ schedule audit (read-only) */

static void audit_print(EDMD3* S, const char* kind, const char* what, int a, int b, double tbf, double theap){
    if (S->aprinted >= 50) return;
    S->aprinted++;
    printf("[EDMD3-AUDIT] %s %s a=%d b=%d now=%.17g t_bruteforce=%.17g t_heap=%.17g\n", kind, what, a, b, S->now, tbf, theap);
}
static inline uint64_t akey(int a, int b){ return ((uint64_t)(uint32_t)a << 32) | (uint32_t)b; }
static long atable_slot(const EDMD3* S, uint64_t k){
    uint64_t h = k * 0x9E3779B97F4A7C15ULL; long m = S->ahcap - 1, p = (long)(h >> 20) & m;
    while (S->ahk[p] != UINT64_MAX && S->ahk[p] != k) p = (p + 1) & m;
    return p;
}
static void audit_cmp(EDMD3* S, double tbf, double th, long* cmp, long* dtc, long* dtrel, const char* what, int a, int b){
    (*cmp)++;
    const double d = fabs(tbf - th);
    if (d > S->A.max_dt) S->A.max_dt = d;
    if (d > 1e-9) {
        const double hz = tbf - S->now, rel = hz > 0.0 ? d / hz : INFINITY;
        (*dtc)++;
        if (rel > S->A.max_rel) S->A.max_rel = rel;
        if (rel > 1e-10) { (*dtrel)++; audit_print(S, "dt_rel", what, a, b, tbf, th); }
    }
}
static void schedule_audit(EDMD3* S){
    const int N = S->N; EDMD3_Audit* A = &S->A;
    if (!S->alx) {
        S->alx = (double*)malloc((size_t)N * sizeof(double)); S->aly = (double*)malloc((size_t)N * sizeof(double));
        S->acr = (double*)malloc((size_t)N * sizeof(double)); S->acd = (int*)malloc((size_t)N * sizeof(int));
        S->awl = (double*)malloc((size_t)4 * N * sizeof(double));
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
        S->ahk = (uint64_t*)malloc((size_t)cap * sizeof(uint64_t)); S->aht = (double*)malloc((size_t)cap * sizeof(double));
        S->ahm = (unsigned char*)malloc((size_t)cap); S->ahcap = cap;
    }
    for (long p = 0; p < S->ahcap; ++p) S->ahk[p] = UINT64_MAX;
    /* the live heap */
    for (long q = 0; q < S->heap.n; ++q) {
        const Ev3* e = &S->heap.d[q];
        if (!ev_live(S, e)) continue;
        double* slot = NULL;
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
        else { if (fabs(*slot - e->t) > 1e-9) A->dup_disagree++; if (e->t < *slot) *slot = e->t; }
    }
    /* crossings, from the local positions at now */
    for (int i = 0; i < N; ++i) {
        const Disk3* D = &S->D[i]; double dt;
        const int dir = cross_rule(S, D, S->alx[i], S->aly[i], &dt);
        if (dir < 0) { if (!isnan(S->acr[i])) { A->cross_extra++; audit_print(S, "extra", "CROSS", i, S->acd[i], NAN, S->acr[i]); } continue; }
        if (isnan(S->acr[i])) { A->cross_missing++; audit_print(S, "missing", "CROSS", i, dir, S->now + dt, NAN); continue; }
        if (S->acd[i] != dir) { A->cross_missing++; audit_print(S, "wrong-direction", "CROSS", i, dir, S->now + dt, S->acr[i]); }
        audit_cmp(S, S->now + dt, S->acr[i], &A->cross_cmp, &A->cross_dt, &A->cross_dt_rel, "CROSS", i, dir);
    }
    /* walls and pairs, brute force from ABSOLUTE positions (an independent computation) */
    for (int i = 0; i < N; ++i) {
        const Disk3* D = &S->D[i];
        const double x = (double)D->cx * S->w + S->alx[i], y = (double)D->cy * S->w + S->aly[i];
        for (int s = 0; s < 4; ++s) {
            double gap, speed;
            switch (s) {
                case 0: if (D->vx >= 0.0) continue; gap = x - S->R; speed = -D->vx; break;
                case 1: if (D->vx <= 0.0) continue; gap = (S->boxW - S->R) - x; speed = D->vx; break;
                case 2: if (D->vy >= 0.0) continue; gap = y - S->R; speed = -D->vy; break;
                default: if (D->vy <= 0.0) continue; gap = (S->boxH - S->R) - y; speed = D->vy; break;
            }
            const double tbf = S->now + (gap <= 0.0 ? 0.0 : gap / speed);
            const double th = S->awl[4L * i + s];
            if (isnan(th)) {
                if (wall_candidate(S, D->cx, D->cy, s)) { A->wall_missing++; audit_print(S, "missing", "WALL", i, s, tbf, NAN); }
                else {
                    A->wall_deferred++;
                    const double tc = isnan(S->acr[i]) ? INFINITY : S->acr[i];
                    if (tbf < tc - 1e-9 * fmax(1.0, tc - S->now)) { A->wall_deferred_early++; audit_print(S, "deferred-early", "WALL", i, s, tbf, tc); }
                }
            } else audit_cmp(S, tbf, th, &A->wall_cmp, &A->wall_dt, &A->wall_dt_rel, "WALL", i, s);
            S->awl[4L * i + s] = NAN;     /* consumed: what is left afterwards is extra */
        }
    }
    for (long k = 0; k < 4L * N; ++k) if (!isnan(S->awl[k])) { A->wall_extra++; audit_print(S, "extra", "WALL", (int)(k / 4), (int)(k % 4), NAN, S->awl[k]); }
    for (int i = 0; i < N; ++i) {
        const Disk3* Di = &S->D[i];
        const double xi = (double)Di->cx * S->w + S->alx[i], yi = (double)Di->cy * S->w + S->aly[i];
        for (int j = i + 1; j < N; ++j) {
            if (mutual_last(S, i, j)) continue;
            const Disk3* Dj = &S->D[j];
            const double rx = ((double)Dj->cx * S->w + S->alx[j]) - xi, ry = ((double)Dj->cy * S->w + S->aly[j]) - yi;
            double dt = 0.0, c = 0.0;
            const int rc = pair_rule(rx, ry, Dj->vx - Di->vx, Dj->vy - Di->vy, S->d2, &dt, &c);
            if (!rc) continue;              /* no brute-force event; a heap entry left unmatched is counted as extra below */
            const long p = atable_slot(S, akey(i, j));
            const int inheap = S->ahk[p] != UINT64_MAX;
            const double tbf = S->now + dt;
            if (!inheap) {
                const int nb = abs(Di->cx - Dj->cx) <= 1 && abs(Di->cy - Dj->cy) <= 1;
                if (nb) { A->pair_missing++; audit_print(S, "missing", "PAIR", i, j, tbf, NAN); }
                else {
                    A->pair_deferred++;
                    const double ci = isnan(S->acr[i]) ? INFINITY : S->acr[i], cj = isnan(S->acr[j]) ? INFINITY : S->acr[j];
                    const double tc = ci < cj ? ci : cj;
                    if (tbf < tc - 1e-9 * fmax(1.0, tc - S->now)) { A->pair_deferred_early++; audit_print(S, "deferred-early", "PAIR", i, j, tbf, tc); }
                }
                continue;
            }
            S->ahm[p] = 1;
            audit_cmp(S, tbf, S->aht[p], &A->pair_cmp, &A->pair_dt, &A->pair_dt_rel, "PAIR", i, j);
        }
    }
    for (long p = 0; p < S->ahcap; ++p)
        if (S->ahk[p] != UINT64_MAX && !S->ahm[p]) {
            A->pair_extra++;
            audit_print(S, "extra", "PAIR", (int)(S->ahk[p] >> 32), (int)(S->ahk[p] & 0xffffffffu), NAN, S->aht[p]);
        }
    A->audits++;
}

/* ------------------------------------------------------------------ public API */

EDMD3* edmd3_create(const EDMD_Params* prm, double cell_px, char* err, size_t errlen){
    const double w = cell_px > 0.0 ? cell_px : EDMD3_DEFAULT_CELL_PX;
#define FAIL(...) do { if (err && errlen) snprintf(err, errlen, __VA_ARGS__); return NULL; } while (0)
    if (!prm || prm->N <= 0) FAIL("edmd_gen3: no disks");
    if (!(prm->radius > 0.0)) FAIL("edmd_gen3: radius must be > 0");
    if (!(w >= 2.0 * prm->radius)) FAIL("edmd_gen3: cell width %.17g px < diameter %.17g px (sec. 4.7.1, a)", w, 2.0 * prm->radius);
    if (w != floor(w) || w > 1048576.0) FAIL("edmd_gen3: cell width %.17g px must be an integer number of px (cx*w exact)", w);
    if (!(prm->boxW >= 2.0 * prm->radius) || !(prm->boxH >= 2.0 * prm->radius)) FAIL("edmd_gen3: box smaller than a disk");
    if (prm->divider_count > 0) FAIL("edmd_gen3 M1: dividers are not implemented yet (M2)");
    if (prm->has_pistonL || prm->has_pistonR) FAIL("edmd_gen3 M1: pistons are not implemented yet (M2)");
    if (prm->heatbath_enabled) FAIL("edmd_gen3: the outer-wall heat bath is not implemented");
    if (prm->species) FAIL("edmd_gen3: species and semipermeable gates are not implemented");
    if (!prm->pp_collisions_enabled) FAIL("edmd_gen3: pair collisions cannot be switched off");
#undef FAIL
    EDMD3* S = (EDMD3*)calloc(1, sizeof(EDMD3));
    if (!S) return NULL;
    S->prm = *prm; S->N = prm->N; S->R = prm->radius; S->d = 2.0 * prm->radius; S->d2 = S->d * S->d;
    S->w = w; S->boxW = prm->boxW; S->boxH = prm->boxH;
    S->gw = (int)ceil(S->boxW / w); S->gh = (int)ceil(S->boxH / w);
    if (S->gw < 1) S->gw = 1;
    if (S->gh < 1) S->gh = 1;
    S->ncell = S->gw * S->gh;
    { const double q = w / S->d + 1.0; S->ccap = (int)(q * q) + 4; }   /* > the most non-overlapping disks a cell can hold */
    S->cell = (int*)calloc((size_t)S->ncell * (size_t)S->ccap, sizeof(int));
    S->ccount = (int*)calloc((size_t)S->ncell, sizeof(int));
    S->D = (Disk3*)calloc((size_t)S->N, sizeof(Disk3));
    S->out = (EDMD_Particle*)calloc((size_t)S->N, sizeof(EDMD_Particle));
    S->heap.cap = 32L * S->N + 1024; S->heap.d = (Ev3*)calloc((size_t)S->heap.cap, sizeof(Ev3));
    if (!S->cell || !S->ccount || !S->D || !S->out || !S->heap.d) { edmd3_destroy(S); if (err && errlen) snprintf(err, errlen, "edmd_gen3: out of memory"); return NULL; }
    S->compact_at = 64L * S->N + 4096;
    /* the experiment validator's tolerances (experiment_validation.c), so a finding means the same thing */
    S->tol_pair = fmax(1e-7, 1e-6 * S->d);
    S->tol_wall = fmax(1e-6, 1e-6 * fmax(1.0, S->R));
    S->tol_cell = 1e-9;
    S->c_tol = 64.0 * S->d2 * 2.220446049250313e-16;   /* c = |r|^2 - d^2 within 64 ulp(d^2): a contact, not an overlap */
    S->check_interval = 24.0;                            /* 1 sigma-time */
    S->same_limit = 5000L > 4L * S->N ? 5000L : 4L * S->N;
    S->hash = 1469598103934665603ULL;
    return S;
}

void edmd3_destroy(EDMD3* S){
    if (!S) return;
    free(S->cell); free(S->ccount); free(S->D); free(S->out); free(S->heap.d);
    free(S->alx); free(S->aly); free(S->acr); free(S->acd); free(S->awl); free(S->ahk); free(S->aht); free(S->ahm);
    free(S);
}

int edmd3_load(EDMD3* S, const EDMD_Particle* P, double t_abs, char* err, size_t errlen){
    if (!S || !P) return 0;
    S->T0 = EDMD3_ORIGIN_SHIFT * floor(t_abs / EDMD3_ORIGIN_SHIFT);
    S->now = t_abs - S->T0;
    S->heap.n = 0;
    memset(S->ccount, 0, (size_t)S->ncell * sizeof(int));
    for (int i = 0; i < S->N; ++i) {
        Disk3* A = &S->D[i];
        int cx = (int)floor(P[i].x / S->w), cy = (int)floor(P[i].y / S->w);
        cx = cx < 0 ? 0 : (cx >= S->gw ? S->gw - 1 : cx);
        cy = cy < 0 ? 0 : (cy >= S->gh ? S->gh - 1 : cy);
        A->cx = cx; A->cy = cy;
        A->xi = P[i].x - (double)cx * S->w; A->zeta = P[i].y - (double)cy * S->w;
        A->vx = P[i].vx; A->vy = P[i].vy; A->tau = S->now;
        A->cnt = 0; A->last = -1; A->pad = 0;
        if (!cell_insert(S, i)) { if (err && errlen) snprintf(err, errlen, "edmd_gen3: cell (%d,%d) overflows at disk %d (overlapping input?)", cx, cy, i); return 0; }
    }
    const long f0 = S->H.full_findings;
    full_check(S);
    if (S->H.full_findings != f0) {
        if (err && errlen) snprintf(err, errlen, "edmd_gen3: the input overlaps or leaves the box (worst surface gap %.3g px)", S->H.full_worst);
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
    S->next_check = S->now + S->check_interval;
    S->virial_t0_abs = t_abs;
    return 1;
}

double edmd3_advance_to(EDMD3* S, double t_abs){
    if (!S || S->fatal) return S ? S->T0 + S->now : 0.0;
    for (;;) {
        const double tt = t_abs - S->T0;
        if (S->heap.n == 0 || S->heap.d[0].t > tt) { if (tt > S->now) S->now = tt; break; }
        Ev3 e; heap_pop(&S->heap, &e);
        if (!ev_live(S, &e)) { S->H.ev_stale++; continue; }
        if (e.t < S->now) { S->H.past_event++; continue; }
        if (e.t == S->same_t) {
            if (++S->same_n > S->same_limit) {
                S->H.stagnation++; S->fatal = 1;
                snprintf(S->fatal_msg, sizeof S->fatal_msg, "stagnation guard: more than %ld live events at t = %.17g", S->same_limit, S->T0 + e.t);
                heap_push(S, e.t, e.type, e.a, e.b, e.ca, e.cb);
                break;
            }
        } else { S->same_t = e.t; S->same_n = 1; }
        S->now = e.t;
        fnv(&S->hash, &e.t, sizeof e.t); fnv(&S->hash, &e.type, sizeof e.type); fnv(&S->hash, &e.a, sizeof e.a); fnv(&S->hash, &e.b, sizeof e.b);
        if (e.type == T_PAIR) exec_pair(S, &e);
        else if (e.type == T_WALL) exec_wall(S, &e);
        else exec_cross(S, &e);
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
        const Disk3* A = &S->D[i]; double x, y; local_now(S, i, &x, &y);
        S->out[i].x = (double)A->cx * S->w + x; S->out[i].y = (double)A->cy * S->w + y;
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
    double ke = 0.0;
    for (int i = 0; i < S->N; ++i) ke += 0.5 * (S->D[i].vx * S->D[i].vx + S->D[i].vy * S->D[i].vy);
    return ke;
}
void   edmd3_reset_virial(EDMD3* S){ S->virial_accum = 0.0; S->virial_pair_events = 0; S->virial_t0_abs = S->T0 + S->now; }
double edmd3_compressibility_Z(EDMD3* S){
    const double t = (S->T0 + S->now) - S->virial_t0_abs, ke = edmd3_kinetic_energy(S);
    if (!(t > 0.0) || !(ke > 0.0)) return NAN;
    return 1.0 + S->virial_accum / (2.0 * ke * t);
}
long   edmd3_virial_pair_events(const EDMD3* S){ return S->virial_pair_events; }
double edmd3_wall_impulse(const EDMD3* S, int wall){ return (wall >= 0 && wall < 4) ? S->wall_impulse[wall] : NAN; }
long   edmd3_wall_events(const EDMD3* S, int wall){ return (wall >= 0 && wall < 4) ? S->wall_events[wall] : 0; }
void   edmd3_set_check_interval(EDMD3* S, double units){ S->check_interval = units; S->next_check = S->now + (units > 0.0 ? units : 0.0); }
void   edmd3_set_contact_audit(EDMD3* S, int on){ S->contact_audit = on ? 1 : 0; }
long   edmd3_contact_audit_stats(const EDMD3* S, double max_gap_px[2]){ max_gap_px[0] = S->contact_max[0]; max_gap_px[1] = S->contact_max[1]; return S->contact_events; }
void   edmd3_set_schedule_audit(EDMD3* S, long every){ S->audit_every = every; S->audit_count = 0; }
void   edmd3_schedule_audit_now(EDMD3* S){ schedule_audit(S); }
const EDMD3_Audit* edmd3_schedule_audit_stats(const EDMD3* S){ return &S->A; }
