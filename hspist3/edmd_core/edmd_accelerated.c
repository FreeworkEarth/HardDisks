// Accelerated EDMD core: provides the same API as edmd.c, but exported under an
// `edmd_acc_*` prefix so it can be linked alongside the default core and
// selected at runtime by the main program.
#define edmd_create edmd_acc_create
#define edmd_destroy edmd_acc_destroy
#define edmd_init_random_gas edmd_acc_init_random_gas
#define edmd_time edmd_acc_time
#define edmd_params edmd_acc_params
#define edmd_particles edmd_acc_particles
#define edmd_count edmd_acc_count
#define edmd_advance_to edmd_acc_advance_to
#define edmd_advance_pp_only_to edmd_acc_advance_pp_only_to
#define edmd_reschedule_all edmd_acc_reschedule_all
#define edmd_reschedule_all_pp_only edmd_acc_reschedule_all_pp_only
#define edmd_config_dividers edmd_acc_config_dividers
#define edmd_set_divider_motions edmd_acc_set_divider_motions
#define edmd_config_divider edmd_acc_config_divider
#define edmd_set_divider_motion edmd_acc_set_divider_motion
#define edmd_config_pistons edmd_acc_config_pistons
#define edmd_total_kinetic_energy edmd_acc_total_kinetic_energy
#define edmd_gas_temperature edmd_acc_gas_temperature
#define edmd_work_divider edmd_acc_work_divider
#define edmd_work_divider_i edmd_acc_work_divider_i
#define edmd_work_pistonL edmd_acc_work_pistonL
#define edmd_work_pistonR edmd_acc_work_pistonR
#define edmd_heat_bath edmd_acc_heat_bath
#define edmd_reset_work edmd_acc_reset_work
#define edmd_divider_resolve_overlaps edmd_acc_divider_resolve_overlaps
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

/* event kind: AB = particle-particle, WL/WR/WB/WT = walls (left/right/bottom/top)
 * DL/DR = divider faces, PL/PR = pistons, CC = grid-cell crossing (for neighbor scheduling correctness)
 */
typedef enum { EV_AB=0, EV_WL, EV_WR, EV_WB, EV_WT, EV_DL, EV_DR, EV_PL, EV_PR, EV_CC } EvType;
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
        case EV_CC: return "CC";
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

    Cell* grid; int gw, gh;    /* grid width/height */
    double cell_size;

    Heap heap;
};

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
    for(int i=0;i<S->prm.N;i++){
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
        int cx = (int)floor(S->P[i].x / S->cell_size);
        int cy = (int)floor(S->P[i].y / S->cell_size);
        cx = clampi(cx, 0, S->gw-1);
        cy = clampi(cy, 0, S->gh-1);
        Cell* c = &S->grid[cy*S->gw + cx];
        cell_reserve(c, c->count+1);
        c->idx[c->count++] = i;
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
    double disc = b*b - vv*c;
    if(disc<=0.0 || vv<=0.0) return 0;
    double t = (-b - sqrt(disc)) / vv;     /* earlier root */
    if(t<=1e-12) return 0;
    *tcol = t; return 1;
}

/* left wall x=0: hit when x - R = 0; particle must move left (vx<0). */
static int collide_time_wall_L(const EDMD* S, const EDMD_Particle* A, double* tcol){
    (void)S;
    if(A->vx >= 0.0) return 0;
    double dist = (A->x - S->prm.radius) - 0.0;
    double t = -dist / A->vx;
    if(t<=1e-12) return 0;
    *tcol = t; return 1;
}
/* right wall x=boxW: hit when x + R = boxW; particle must move right (vx>0). */
static int collide_time_wall_R(const EDMD* S, const EDMD_Particle* A, double* tcol){
    if(A->vx <= 0.0) return 0;
    double dist = (S->prm.boxW - S->prm.radius) - A->x;
    double t = dist / A->vx;
    if(t<=1e-12) return 0;
    *tcol = t; return 1;
}
/* bottom wall y=0: hit when y - R = 0; particle must move down (vy<0). */
static int collide_time_wall_B(const EDMD* S, const EDMD_Particle* A, double* tcol){
    (void)S;
    if(A->vy >= 0.0) return 0;
    double dist = (A->y - S->prm.radius) - 0.0;
    double t = -dist / A->vy;
    if(t<=1e-12) return 0;
    *tcol = t; return 1;
}
/* top wall y=boxH: hit when y + R = boxH; particle must move up (vy>0). */
static int collide_time_wall_T(const EDMD* S, const EDMD_Particle* A, double* tcol){
    if(A->vy <= 0.0) return 0;
    double dist = (S->prm.boxH - S->prm.radius) - A->y;
    double t = dist / A->vy;
    if(t<=1e-12) return 0;
    *tcol = t; return 1;
}

/* ------------------------ scheduling ------------------------ */

static void schedule_walls(EDMD* S, int i){
    double t;
    if(collide_time_wall_L(S, &S->P[i], &t))
        heap_push(&S->heap, (Event){ S->t+t, i,-1, S->P[i].coll_count,0, EV_WL });
    if(collide_time_wall_R(S, &S->P[i], &t))
        heap_push(&S->heap, (Event){ S->t+t, i,-1, S->P[i].coll_count,0, EV_WR });
    if(collide_time_wall_B(S, &S->P[i], &t))
        heap_push(&S->heap, (Event){ S->t+t, i,-1, S->P[i].coll_count,0, EV_WB });
    if(collide_time_wall_T(S, &S->P[i], &t))
        heap_push(&S->heap, (Event){ S->t+t, i,-1, S->P[i].coll_count,0, EV_WT });
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
        const double th = S->prm.divider_thickness[d];
        const double x0 = S->prm.divider_x[d];
        const double v0 = S->prm.divider_vx[d];
        const double k  = S->prm.divider_k[d];
        const double mW = S->prm.divider_mass[d];
        const double xeq = S->prm.divider_xeq[d];
        const double Rpart = S->prm.radius;
        const double gap = (x0 - 0.5*th) - (A->x + Rpart);
        if (gap <= 1e-12) return 0;

        const double w = sqrt(k / mW);
        if (!(w > 0.0)) return 0;
        const double T = 2.0 * M_PI / w;
        const double amp = hypot((x0 - xeq), (v0 / w));
        const double base = (A->x - xeq) + (Rpart + 0.5*th);
        const double v = A->vx;

        double t_max;
        if (v > 0.0) t_max = (fabs(base) + amp + 1.0) / v + 2.0 * T;
        else { if (base + amp < 0.0) return 0; t_max = fmax(2.0 * T, 0.5 * T); }
        if (!(t_max > 0.0)) return 0;

        const double coarse_step = fmin(T / 16.0, fmax(1e-6, 0.25 * gap / (fabs(v) + fabs(v0) + 1e-9)));
        double t0s = 0.0;
        double f0s = base + v * t0s - (x0 - xeq) * cos(w * t0s) - (v0 / w) * sin(w * t0s);
        double ta = NAN, tb = NAN;
        for (int it = 0; it < 20000; ++it) {
            double t1s = t0s + coarse_step;
            if (t1s > t_max) t1s = t_max;
            double f1s = base + v * t1s - (x0 - xeq) * cos(w * t1s) - (v0 / w) * sin(w * t1s);
            if ((f0s <= 0.0 && f1s >= 0.0) || (f0s >= 0.0 && f1s <= 0.0)) { ta = t0s; tb = t1s; break; }
            t0s = t1s; f0s = f1s;
            if (t0s >= t_max - 1e-15) break;
        }
        if (!isfinite(ta) || !isfinite(tb) || !(tb > ta)) return 0;

        const double fine_step = fmin(T / 256.0, fmax(1e-7, (tb - ta) / 64.0));
        double t_prev = ta;
        double f_prev = base + v * t_prev - (x0 - xeq) * cos(w * t_prev) - (v0 / w) * sin(w * t_prev);
        double t_lo = ta, t_hi = tb;
        for (int it = 0; it < 4096; ++it) {
            double t_cur = t_prev + fine_step;
            if (t_cur > tb) t_cur = tb;
            double f_cur = base + v * t_cur - (x0 - xeq) * cos(w * t_cur) - (v0 / w) * sin(w * t_cur);
            if ((f_prev <= 0.0 && f_cur >= 0.0) || (f_prev >= 0.0 && f_cur <= 0.0)) { t_lo = t_prev; t_hi = t_cur; break; }
            t_prev = t_cur; f_prev = f_cur;
            if (t_prev >= tb - 1e-15) break;
        }

        double flo = base + v * t_lo - (x0 - xeq) * cos(w * t_lo) - (v0 / w) * sin(w * t_lo);
        double fhi = base + v * t_hi - (x0 - xeq) * cos(w * t_hi) - (v0 / w) * sin(w * t_hi);
        if (!((flo <= 0.0 && fhi >= 0.0) || (flo >= 0.0 && fhi <= 0.0))) return 0;
        for (int it = 0; it < 80; ++it) {
            double tm = 0.5 * (t_lo + t_hi);
            double fm = base + v * tm - (x0 - xeq) * cos(w * tm) - (v0 / w) * sin(w * tm);
            if ((flo <= 0.0 && fm >= 0.0) || (flo >= 0.0 && fm <= 0.0)) { t_hi = tm; fhi = fm; }
            else { t_lo = tm; flo = fm; }
        }
        double t = 0.5 * (t_lo + t_hi);
        if (t <= 1e-12) return 0;
        *tcol = t;
        return 1;
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
        const double th = S->prm.divider_thickness[d];
        const double x0 = S->prm.divider_x[d];
        const double v0 = S->prm.divider_vx[d];
        const double k  = S->prm.divider_k[d];
        const double mW = S->prm.divider_mass[d];
        const double xeq = S->prm.divider_xeq[d];
        const double Rpart = S->prm.radius;
        const double gap = (A->x - Rpart) - (x0 + 0.5*th);
        if (gap <= 1e-12) return 0;

        const double w = sqrt(k / mW);
        if (!(w > 0.0)) return 0;
        const double T = 2.0 * M_PI / w;
        const double amp = hypot((x0 - xeq), (v0 / w));
        const double base = (A->x - xeq) - (Rpart + 0.5*th);
        const double v = A->vx;

        double t_max;
        if (v < 0.0) t_max = (fabs(base) + amp + 1.0) / (-v) + 2.0 * T;
        else { if (base - amp > 0.0) return 0; t_max = fmax(2.0 * T, 0.5 * T); }
        if (!(t_max > 0.0)) return 0;

        const double coarse_step = fmin(T / 16.0, fmax(1e-6, 0.25 * gap / (fabs(v) + fabs(v0) + 1e-9)));
        double t0s = 0.0;
        double f0s = base + v * t0s - (x0 - xeq) * cos(w * t0s) - (v0 / w) * sin(w * t0s);
        double ta = NAN, tb = NAN;
        for (int it = 0; it < 20000; ++it) {
            double t1s = t0s + coarse_step;
            if (t1s > t_max) t1s = t_max;
            double f1s = base + v * t1s - (x0 - xeq) * cos(w * t1s) - (v0 / w) * sin(w * t1s);
            if ((f0s <= 0.0 && f1s >= 0.0) || (f0s >= 0.0 && f1s <= 0.0)) { ta = t0s; tb = t1s; break; }
            t0s = t1s; f0s = f1s;
            if (t0s >= t_max - 1e-15) break;
        }
        if (!isfinite(ta) || !isfinite(tb) || !(tb > ta)) return 0;

        const double fine_step = fmin(T / 256.0, fmax(1e-7, (tb - ta) / 64.0));
        double t_prev = ta;
        double f_prev = base + v * t_prev - (x0 - xeq) * cos(w * t_prev) - (v0 / w) * sin(w * t_prev);
        double t_lo = ta, t_hi = tb;
        for (int it = 0; it < 4096; ++it) {
            double t_cur = t_prev + fine_step;
            if (t_cur > tb) t_cur = tb;
            double f_cur = base + v * t_cur - (x0 - xeq) * cos(w * t_cur) - (v0 / w) * sin(w * t_cur);
            if ((f_prev <= 0.0 && f_cur >= 0.0) || (f_prev >= 0.0 && f_cur <= 0.0)) { t_lo = t_prev; t_hi = t_cur; break; }
            t_prev = t_cur; f_prev = f_cur;
            if (t_prev >= tb - 1e-15) break;
        }

        double flo = base + v * t_lo - (x0 - xeq) * cos(w * t_lo) - (v0 / w) * sin(w * t_lo);
        double fhi = base + v * t_hi - (x0 - xeq) * cos(w * t_hi) - (v0 / w) * sin(w * t_hi);
        if (!((flo <= 0.0 && fhi >= 0.0) || (flo >= 0.0 && fhi <= 0.0))) return 0;
        for (int it = 0; it < 80; ++it) {
            double tm = 0.5 * (t_lo + t_hi);
            double fm = base + v * tm - (x0 - xeq) * cos(w * tm) - (v0 / w) * sin(w * tm);
            if ((flo <= 0.0 && fm >= 0.0) || (flo >= 0.0 && fm <= 0.0)) { t_hi = tm; fhi = fm; }
            else { t_lo = tm; flo = fm; }
        }
        double t = 0.5 * (t_lo + t_hi);
        if (t <= 1e-12) return 0;
        *tcol = t;
        return 1;
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
    double t;
    if(collide_time_ab(S->P[i].x,S->P[i].y,S->P[i].vx,S->P[i].vy,
                       S->P[j].x,S->P[j].y,S->P[j].vx,S->P[j].vy,
                       S->prm.radius, &t))
    {
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
            if (S->prm.divider_x[d] < left_limit) { S->prm.divider_x[d] = left_limit; if (S->prm.divider_vx[d] < 0.0) S->prm.divider_vx[d] = 0.0; }
            if (S->prm.divider_x[d] > right_limit) { S->prm.divider_x[d] = right_limit; if (S->prm.divider_vx[d] > 0.0) S->prm.divider_vx[d] = 0.0; }
        } else {
            S->prm.divider_x[d] += S->prm.divider_vx[d] * dt;
            if (S->prm.divider_x[d] < left_limit) { S->prm.divider_x[d] = left_limit; }
            if (S->prm.divider_x[d] > right_limit) { S->prm.divider_x[d] = right_limit; }
        }
    }
}

static inline void particle_cell_xy(const EDMD* S, int i, int* out_cx, int* out_cy){
    int cx = (int)floor(S->P[i].x / S->cell_size);
    int cy = (int)floor(S->P[i].y / S->cell_size);
    cx = clampi(cx, 0, S->gw - 1);
    cy = clampi(cy, 0, S->gh - 1);
    *out_cx = cx; *out_cy = cy;
}

static void schedule_cell_crossing(EDMD* S, int i){
    if (S->cell_size <= 0.0) return;
    const double cs = S->cell_size;
    const double x = S->P[i].x, y = S->P[i].y;
    const double vx = S->P[i].vx, vy = S->P[i].vy;
    const double eps = fmax(1e-9, 1e-9 * cs);

    /* Use a tiny motion-direction bias for the scheduling cell. Without this,
     * a particle that lands numerically on a grid boundary can repeatedly
     * schedule another almost-zero-time CC event on the same boundary.
     */
    double x_probe = x;
    double y_probe = y;
    if (vx > 0.0) x_probe += eps;
    else if (vx < 0.0) x_probe -= eps;
    if (vy > 0.0) y_probe += eps;
    else if (vy < 0.0) y_probe -= eps;
    x_probe = fmin(fmax(x_probe, 0.0), S->prm.boxW);
    y_probe = fmin(fmax(y_probe, 0.0), S->prm.boxH);

    int cx = clampi((int)floor(x_probe / cs), 0, S->gw - 1);
    int cy = clampi((int)floor(y_probe / cs), 0, S->gh - 1);

    double tx = DBL_MAX, ty = DBL_MAX;
    if (vx > 0.0) {
        const double xb = (double)(cx + 1) * cs;
        tx = (xb - x) / vx;
    } else if (vx < 0.0) {
        const double xb = (double)(cx) * cs;
        tx = (xb - x) / vx; /* vx<0 => positive if xb<x */
    }
    if (vy > 0.0) {
        const double yb = (double)(cy + 1) * cs;
        ty = (yb - y) / vy;
    } else if (vy < 0.0) {
        const double yb = (double)(cy) * cs;
        ty = (yb - y) / vy;
    }

    double t = (tx < ty) ? tx : ty;
    if (!(t > eps) || !isfinite(t)) return;
    heap_push(&S->heap, (Event){ S->t + t, i, -1, S->P[i].coll_count, 0, EV_CC });
}

static void schedule_ab_neighbors_9cell(EDMD* S, int i){
    if (!S->grid) return;
    int cx, cy;
    particle_cell_xy(S, i, &cx, &cy);
    for (int dy = -1; dy <= 1; ++dy) {
        int yy = cy + dy;
        if (yy < 0 || yy >= S->gh) continue;
        for (int dx = -1; dx <= 1; ++dx) {
            int xx = cx + dx;
            if (xx < 0 || xx >= S->gw) continue;
            Cell* c = &S->grid[yy * S->gw + xx];
            for (int k = 0; k < c->count; ++k) {
                int j = c->idx[k];
                if (j == i) continue;
                schedule_ab(S, i, j);
            }
        }
    }
}

/* schedule for one particle i: walls + divider/pistons + cell crossing + AB (local neighbors) */
static void schedule_for(EDMD* S, int i){
    schedule_walls(S, i);
    schedule_divider(S, i);
    schedule_pistons(S, i);
    schedule_cell_crossing(S, i);
    schedule_ab_neighbors_9cell(S, i);
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
    for (int i = 0; i < S->prm.N; ++i) {
        schedule_for(S, i);
    }
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
        if(type==EV_WL || type==EV_WR){ A->vx = -A->vx; }
        if(type==EV_WB || type==EV_WT){ A->vy = -A->vy; }
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
            A->coll_count++;
            return;
        } else {
            double v1 = ((m - M)*u1 + 2.0*M*u2) / (m + M);
            double v2 = ((M - m)*u2 + 2.0*m*u1) / (m + M);
            /* Work: change in divider KE */
            double dE = 0.5 * M * (v2*v2 - u2*u2);
            S->work_divider[d] += dE;
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
            A->coll_count++;
            return;
        }
        double v1 = ((m - M)*u1 + 2.0*M*u2) / (m + M);
        double v2 = ((M - m)*u2 + 2.0*m*u1) / (m + M);
        A->vx = v1;
        double dE = 0.5 * M * (v2*v2 - u2*u2);
        if (type==EV_PL) { S->prm.pistonL_vx = v2; S->work_pistonL += dE; }
        else             { S->prm.pistonR_vx = v2; S->work_pistonR += dE; }
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
    free(S);
}

/* tiny xorshift64 rng */
static unsigned long long xorshift64(unsigned long long* s){
    unsigned long long x = *s;
    x ^= x<<13; x ^= x>>7; x ^= x<<17;
    *s = x; return x;
}

/* random non-overlapping init + Gaussian velocities (Box-Muller) */
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
    long type_counts[10] = {0};
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
            break;
        }
        events_processed++;
        if (e.type >= EV_AB && e.type <= EV_CC) type_counts[(int)e.type]++;
        if (fabs(e.t - last_event_t) <= 1e-13) stagnant_events++;
        else { stagnant_events = 0; last_event_t = e.t; }
        if (events_processed > EDMD_ADVANCE_MAX_EVENTS || stagnant_events > EDMD_ADVANCE_MAX_STAGNANT_EVENTS) {
            g_edmd_avalanche_warning_count++;
            if (g_edmd_avalanche_warning_count <= 20 || (g_edmd_avalanche_warning_count % 1000) == 0) {
                int dominant = 0;
                for (int k = 1; k <= (int)EV_CC; ++k) {
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

        /* validate using coll_count snapshots */
        if(e.a<0 || e.a>=S->prm.N) continue;
        if(S->P[e.a].coll_count != e.ca) continue;
        if(e.type==EV_AB){
            if(e.b<0 || e.b>=S->prm.N) continue;
            if(S->P[e.b].coll_count != e.cb) continue;
        }

        /* resolve */
        if(e.type==EV_AB) {
            resolve_ab(S, e.a, e.b);
            grid_build(S);
            schedule_for(S, e.a);
            schedule_for(S, e.b);
        } else if (e.type==EV_CC) {
            /* Cell-crossing: no collision, but it is a scheduling epoch.
             * Increment coll_count so duplicate stale CC events for this particle
             * do not remain valid and avalanche at the same timestamp.
             */
            S->P[e.a].coll_count++;
            grid_build(S);
            schedule_for(S, e.a);
        } else {
            resolve_wall(S, e.a, e.type, e.b);
            /* moving boundary velocities may have changed (divider/pistons): rebuild and reschedule all */
            if (e.type==EV_DL || e.type==EV_DR || e.type==EV_PL || e.type==EV_PR) {
                reschedule_all_internal(S);
            } else {
                grid_build(S);
                schedule_for(S, e.a);
            }
        }
    }
    return S->t;
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
