#include "edmd.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* ------------------------ internal types ------------------------ */

/* event kind: AB = particle-particle, WL/WR/WB/WT = walls (left/right/bottom/top) */
typedef enum { EV_AB=0, EV_WL, EV_WR, EV_WB, EV_WT, EV_DL, EV_DR, EV_PL, EV_PR } EvType;

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

    Cell* grid; int gw, gh;    /* grid width/height */
    double cell_size;

    Heap heap;
};

/* ------------------------ small helpers ------------------------ */

static inline int clampi(int v, int lo, int hi) { return v<lo?lo:(v>hi?hi:v); }

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
        /* if inside divider slab, push to nearest face and reflect in divider frame */
        if (S->prm.has_divider && S->prm.divider_thickness > 0.0) {
            double L = S->prm.divider_x - 0.5*S->prm.divider_thickness;
            double Rf = S->prm.divider_x + 0.5*S->prm.divider_thickness;
            if (L > 0.0 && Rf < S->prm.boxW) {
                double gapL = (S->P[i].x - R) - L;
                double gapR = Rf - (S->P[i].x + R);
                if (gapL > -eps && gapR > -eps) {
                    if (gapL < gapR) { S->P[i].x = L - eps + R; S->P[i].vx = 2.0*S->prm.divider_vx - S->P[i].vx; }
                    else { S->P[i].x = Rf + eps - R; S->P[i].vx = 2.0*S->prm.divider_vx - S->P[i].vx; }
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
static int collide_time_divider_L(const EDMD* S, const EDMD_Particle* A, double* tcol){
    if(!S->prm.has_divider) return 0;
    double L = S->prm.divider_x - 0.5*S->prm.divider_thickness;
    if(!(L > 0.0 && L < S->prm.boxW)) return 0;
    /* left face, approached from left side: (x+R) + vx t = L + vW t */
    double rel = A->vx - S->prm.divider_vx;
    double num = L - (A->x + S->prm.radius);
    if(num <= 1e-12) return 0;      /* not strictly left of face */
    if(rel <= 0.0) return 0;        /* not approaching in wall frame */
    double t = num / rel;
    if(t<=1e-12) return 0; *tcol=t; return 1;
}
static int collide_time_divider_R(const EDMD* S, const EDMD_Particle* A, double* tcol){
    if(!S->prm.has_divider) return 0;
    double Rf = S->prm.divider_x + 0.5*S->prm.divider_thickness;
    if(!(Rf > 0.0 && Rf < S->prm.boxW)) return 0;
    /* right face, approached from right side: (x-R) + vx t = Rf + vW t */
    double rel = A->vx - S->prm.divider_vx;
    double dist = (A->x - S->prm.radius) - Rf; /* >0 if strictly to the right */
    if(dist <= 1e-12) return 0;      /* not strictly right of face */
    if(rel >= 0.0) return 0;         /* must move left relative to wall */
    double t = dist / (-rel);        /* since (vx-vW) t = -dist */
    if(t<=1e-12) return 0; *tcol=t; return 1;
}

static void schedule_divider(EDMD* S, int i){
    if(!S->prm.has_divider || S->prm.divider_thickness <= 0.0) return;
    double t;
    if(collide_time_divider_L(S, &S->P[i], &t))
        heap_push(&S->heap, (Event){ S->t+t, i,-1, S->P[i].coll_count,0, EV_DL });
    if(collide_time_divider_R(S, &S->P[i], &t))
        heap_push(&S->heap, (Event){ S->t+t, i,-1, S->P[i].coll_count,0, EV_DR });
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

/* reflect on static walls */
static void resolve_wall(EDMD* S, int i, EvType type){
    EDMD_Particle* A=&S->P[i];
    if(type==EV_WL || type==EV_WR){ A->vx = -A->vx; A->coll_count++; return; }
    if(type==EV_WB || type==EV_WT){ A->vy = -A->vy; A->coll_count++; return; }
    if(type==EV_DL || type==EV_DR){
        double m = 1.0;
        double M = (S->prm.divider_mass > 0.0 ? S->prm.divider_mass : 0.0);
        double u1 = A->vx;
        double u2 = S->prm.divider_vx;
        if (M <= 0.0){
            A->vx = 2*u2 - u1; A->coll_count++; return;
        } else {
            double v1 = ((m - M)*u1 + 2.0*M*u2) / (m + M);
            double v2 = ((M - m)*u2 + 2.0*m*u1) / (m + M);
            A->vx = v1; S->prm.divider_vx = v2; A->coll_count++; return;
        }
    }
    if(type==EV_PL || type==EV_PR){
        double m = 1.0;
        double M = (type==EV_PL ? S->prm.pistonL_mass : S->prm.pistonR_mass);
        double u1 = A->vx;
        double u2 = (type==EV_PL ? S->prm.pistonL_vx : S->prm.pistonR_vx);
        if (M <= 0.0){ A->vx = 2*u2 - u1; A->coll_count++; return; }
        double v1 = ((m - M)*u1 + 2.0*M*u2) / (m + M);
        double v2 = ((M - m)*u2 + 2.0*m*u1) / (m + M);
        A->vx = v1;
        if (type==EV_PL) S->prm.pistonL_vx = v2; else S->prm.pistonR_vx = v2;
        A->coll_count++; return;
    }
}

/* ------------------------ public API ------------------------ */

EDMD* edmd_create(const EDMD_Params* prm_in){
    EDMD* S = (EDMD*)calloc(1, sizeof(EDMD));
    S->prm = *prm_in;
    S->t = 0.0;
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
    while(S->t < t_target){
        if(!heap_pop(&S->heap, &e)){
            /* no events: free-flight to target */
            double dt = t_target - S->t;
            for(int i=0;i<S->prm.N;i++){
                S->P[i].x += S->P[i].vx*dt;
                S->P[i].y += S->P[i].vy*dt;
            }
            if (S->prm.has_divider) S->prm.divider_x += S->prm.divider_vx * dt;
            S->t = t_target;
            break;
        }
        if(e.t > t_target){
            /* event lies in future: free-flight to t_target and requeue the event */
            double dt = t_target - S->t;
            for(int i=0;i<S->prm.N;i++){
                S->P[i].x += S->P[i].vx*dt;
                S->P[i].y += S->P[i].vy*dt;
            }
            if (S->prm.has_divider) S->prm.divider_x += S->prm.divider_vx * dt;
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
        if (S->prm.has_divider) S->prm.divider_x += S->prm.divider_vx * dt;
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
        } else {
            resolve_wall(S, e.a, e.type);
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

void edmd_config_divider(EDMD* S, int enabled, double cx, double thickness){
    S->prm.has_divider = enabled ? 1 : 0;
    S->prm.divider_x = cx;
    S->prm.divider_thickness = thickness;
}

void edmd_set_divider_motion(EDMD* S, double mass, double vx){
    S->prm.divider_mass = mass;
    S->prm.divider_vx = vx;
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
