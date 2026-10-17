/* PID Demo 1b: the engine tests of the static walls on hand-built geometry (02_DEMO1_BRIEF.md sec. 6 tests 4-6 as amended
   by sec. 8; pre-check 02_DEMO1_BRIEF_PRECHECK.md part B). Fork only: the original gen3 has no walls. Every geometry runs the
   harnesses' A/B/C pattern (A both audits on, B audits off, C B again: A = B and B = C in event hash, health and every disk),
   with a stop every sigma-time at which every disk centre is checked against every primitive by brute force
   (edmd3_prim_gap: signed, so a centre outside the lumen is caught too).
     corridor ETA TILT T SEED   a closed rectangle of 100 x 40 sigma (four segments) centred in the 200 x 120 sigma box,
                                turned by TILT degrees; N disks at packing ETA of the rectangle's area; T sigma-time; the
                                pressure on the four walls in 10 blocks after T/10 (test 4: TILT 0; test 5: TILT 30)
     y ETA T SEED               a Y of three arms (half-width 5, length 30 sigma) with semicircular caps: 6 segments, 3 arcs,
                                3 reflex vertices, 6 tangent joints (test 6 (a)-(c))
     yaim SPLIT T SEED          the Y with disks aimed exactly at a reflex vertex and at an arc-segment joint, the schedule audit
                                after EVERY event; SPLIT 1 moves one copy of a vertex by 1.7e-12 px, as in real hulls (test 6 (d))
     empty T SEED               a box run with edmd3_set_walls(empty) against one without the call: equal hash, health, disks
   Output: "key value ..." lines; the verdict lines start with "check " and end with "OK" or "FAIL". Exit 0 iff no FAIL.
   Units as gen3: 24 px = 1 sigma, 24 units = 1 sigma-time, unit disk mass, kT = 1; R = 12 px. Seeding: Demo 0's one generator
   (demo0/pid_rng.c, linked unchanged; brief sec. 2 rule 5). */
#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include "edmd_gen3.h"
#include "pid.h"          /* demo0/: the generator only (rng_seed, rng_uniform, rng_normal) */

#define PX   24.0
#define SIGT 24.0
static const double R = 12.0, PI = 3.14159265358979323846;
static const double BOXW = 4800.0, BOXH = 2880.0;   /* the Demo 0 box, 200 x 120 sigma (02_DEMO1_BRIEF.md sec. 3) */
static int nfail = 0;

static Rng rng;

static void check(int ok, const char* what){ printf("check %s %s\n", what, ok ? "OK" : "FAIL"); if (!ok) nfail++; }

/* Kolafa-Rottner 2006 hard-disk Z(eta), the coefficients of hspist3/plot_speed_of_sound_edmd.py:623-639 */
static double z_kr(double eta){
    static const int pw[15] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 22, 57};
    static const double cf[15] = {1.0, 2.0, 1.12801775, 0.00181895291, -0.0526134737, 0.0504960168, -0.0325537792, 0.0134578632,
                                  0.00140888182, -0.00834273601, 0.00694127367, -0.00262254723, 0.000355746352, -5.24672938e-9, 5.88054639e-23};
    const double x = eta / (1.0 - eta); double z = 0.0;
    for (int k = 0; k < 15; ++k) z += cf[k] * pow(x, pw[k]);
    return z;
}

/* ------------------------------------------------------------------ geometry */
typedef struct {
    int nseg, narc; EDMD3_Seg seg[16]; EDMD3_Arc arc[8];
    double area;                     /* the area inside the walls [px^2] (Green's theorem) */
    double area_acc, perim_acc;      /* corridor: the centre-accessible rectangle (sides less 2R) */
    double cx, cy, ux, uy, Lx, Ly;   /* corridor: centre, unit axis, sides */
    double vx_[3], vy_[3];           /* y: the three reflex vertices */
    double jx, jy, jnx, jny;         /* y: one arc-segment joint and the normal into the lumen there */
} Geo;
static double loop_area(const Geo* G){
    double a = 0.0;
    for (int k = 0; k < G->nseg; ++k) a += 0.5 * (G->seg[k].x1 * G->seg[k].y2 - G->seg[k].x2 * G->seg[k].y1);
    for (int k = 0; k < G->narc; ++k) {
        const EDMD3_Arc* A = &G->arc[k];
        a += 0.5 * (A->r * (A->cx * (sin(A->a1) - sin(A->a0)) - A->cy * (cos(A->a1) - cos(A->a0))) + A->r * A->r * (A->a1 - A->a0));
    }
    return a;
}
/* a closed rectangle Lx x Ly [px] centred at (cx, cy), axis at angle th; counter-clockwise, so the lumen is on the left */
static void make_corridor(Geo* G, double Lx, double Ly, double th){
    memset(G, 0, sizeof *G);
    G->cx = 0.5 * BOXW; G->cy = 0.5 * BOXH; G->ux = cos(th); G->uy = sin(th); G->Lx = Lx; G->Ly = Ly;
    const double a = 0.5 * Lx, b = 0.5 * Ly, vx = -G->uy, vy = G->ux;
    double X[4], Y[4];
    const double sa[4] = {-1, 1, 1, -1}, sb[4] = {-1, -1, 1, 1};
    for (int k = 0; k < 4; ++k) { X[k] = G->cx + sa[k] * a * G->ux + sb[k] * b * vx; Y[k] = G->cy + sa[k] * a * G->uy + sb[k] * b * vy; }
    for (int k = 0; k < 4; ++k) { EDMD3_Seg s = {X[k], Y[k], X[(k + 1) % 4], Y[(k + 1) % 4]}; G->seg[G->nseg++] = s; }
    G->area = loop_area(G);
    G->area_acc = (Lx - 2 * R) * (Ly - 2 * R); G->perim_acc = 2.0 * ((Lx - 2 * R) + (Ly - 2 * R));
}
/* a Y: arms at 30, 150 and 270 degrees from the box centre, half-width h, cap centres at distance L, semicircular caps */
static void make_y(Geo* G, double h, double L, int split){
    memset(G, 0, sizeof *G);
    const double c0 = 0.5 * BOXW, c1 = 0.5 * BOXH, th[3] = {PI / 6.0, 5.0 * PI / 6.0, 3.0 * PI / 2.0};
    double Vx[3], Vy[3];                       /* V[k]: between arm k and arm k+1 (counter-clockwise) */
    for (int k = 0; k < 3; ++k) { Vx[k] = c0 + (2.0 * h / sqrt(3.0)) * cos(th[k] + PI / 3.0); Vy[k] = c1 + (2.0 * h / sqrt(3.0)) * sin(th[k] + PI / 3.0); }
    for (int k = 0; k < 3; ++k) {
        const double Cx = c0 + L * cos(th[k]), Cy = c1 + L * sin(th[k]), a0 = th[k] - PI / 2.0, a1 = th[k] + PI / 2.0;
        const double Erx = Cx + h * cos(a0), Ery = Cy + h * sin(a0), Elx = Cx + h * cos(a1), Ely = Cy + h * sin(a1);
        const int kp = (k + 2) % 3;            /* the vertex before arm k */
        EDMD3_Seg right = {Vx[kp], Vy[kp], Erx, Ery}; G->seg[G->nseg++] = right;
        EDMD3_Arc cap = {Cx, Cy, h, a0, a1}; G->arc[G->narc++] = cap;
        EDMD3_Seg left = {Elx, Ely, Vx[k], Vy[k]};
        if (split && k == 0) left.x2 += 1.7e-12;   /* the two copies of vertex 0 differ by 1.7e-12 px */
        G->seg[G->nseg++] = left;
        if (k == 0) { G->jx = Erx; G->jy = Ery; G->jnx = -sin(th[k]); G->jny = cos(th[k]); }   /* the joint at the start of arm 0's cap */
    }
    for (int k = 0; k < 3; ++k) { G->vx_[k] = Vx[k]; G->vy_[k] = Vy[k]; }
    G->area = loop_area(G);
}

/* ------------------------------------------------------------------ states */
static EDMD_Params box_params(int N){
    EDMD_Params p; memset(&p, 0, sizeof p);
    p.boxW = BOXW; p.boxH = BOXH; p.radius = R; p.N = N; p.pp_collisions_enabled = 1; p.particle_mass = 1.0; p.kB = 1.0;
    return p;
}
static EDMD3* make_state(const Geo* G, int N, int walls){
    char err[512]; EDMD_Params p = box_params(N);
    EDMD3* S = edmd3_create(&p, 0.0, err, sizeof err);
    if (!S) { fprintf(stderr, "create: %s\n", err); exit(2); }
    if (walls && !edmd3_set_walls(S, G->seg, G->nseg, G->arc, G->narc, err, sizeof err)) { fprintf(stderr, "walls: %s\n", err); exit(2); }
    return S;
}
/* random sequential insertion of N centres inside the walls (signed gap >= 1e-6 px to every wall, no overlap), MB velocities */
static EDMD_Particle* seed(const Geo* G, int N, uint64_t s){
    EDMD3* W = make_state(G, 1, 1);           /* a one-disk state used only for edmd3_prim_gap */
    EDMD_Particle* P = (EDMD_Particle*)calloc((size_t)N, sizeof(EDMD_Particle));
    double x0 = INFINITY, x1 = -INFINITY, y0 = INFINITY, y1 = -INFINITY;
    for (int k = 0; k < G->nseg; ++k) { x0 = fmin(x0, fmin(G->seg[k].x1, G->seg[k].x2)); x1 = fmax(x1, fmax(G->seg[k].x1, G->seg[k].x2));
                                         y0 = fmin(y0, fmin(G->seg[k].y1, G->seg[k].y2)); y1 = fmax(y1, fmax(G->seg[k].y1, G->seg[k].y2)); }
    for (int k = 0; k < G->narc; ++k) { x0 = fmin(x0, G->arc[k].cx - G->arc[k].r); x1 = fmax(x1, G->arc[k].cx + G->arc[k].r);
                                         y0 = fmin(y0, G->arc[k].cy - G->arc[k].r); y1 = fmax(y1, G->arc[k].cy + G->arc[k].r); }
    rng_seed(&rng, s);
    for (int i = 0; i < N; ++i) {
        for (long tries = 0;; ++tries) {
            if (tries > 50000000L) { fprintf(stderr, "seeding failed at disk %d of %d\n", i, N); exit(2); }
            const double x = x0 + rng_uniform(&rng) * (x1 - x0), y = y0 + rng_uniform(&rng) * (y1 - y0);
            if (!(edmd3_prim_gap(W, x, y) >= 1e-6)) continue;
            int ok = 1;
            for (int j = 0; j < i && ok; ++j) { const double dx = x - P[j].x, dy = y - P[j].y; if (dx * dx + dy * dy < 4.0 * R * R + 1e-6) ok = 0; }
            if (ok) { P[i].x = x; P[i].y = y; break; }
        }
    }
    double px = 0.0, py = 0.0, ke = 0.0;
    for (int i = 0; i < N; ++i) {
        P[i].vx = rng_normal(&rng); P[i].vy = rng_normal(&rng); px += P[i].vx; py += P[i].vy;
    }
    px /= N; py /= N;
    for (int i = 0; i < N; ++i) { P[i].vx -= px; P[i].vy -= py; ke += 0.5 * (P[i].vx * P[i].vx + P[i].vy * P[i].vy); }
    const double sc = sqrt((double)N / ke);
    for (int i = 0; i < N; ++i) { P[i].vx *= sc; P[i].vy *= sc; }
    edmd3_destroy(W);
    return P;
}

/* ------------------------------------------------------------------ one run */
typedef struct {
    uint64_t hash; EDMD3_Health H; EDMD3_Audit A; EDMD3_Tol tol; double cmax[2]; long cev;
    double ke0, ke1; int clean, fatal; double min_gap; long outside, stops;
    EDMD_Particle* fin; int N;
    double Jb[10], Jend; int nb;    /* the walls' normal impulse per pressure block; the endpoint impulse over the run */
    double tb;                      /* block length [units] */
} GRun;
/* audits: 0 off; else the schedule audit after every event for the first `full` physical events, then every `every` */
static void run(const Geo* G, const EDMD_Particle* P0, int N, double T, int audits, long full, long every, int blocks, GRun* r){
    memset(r, 0, sizeof *r); r->N = N; r->min_gap = INFINITY;
    EDMD3* S = make_state(G, N, 1);
    char err[512];
    if (!edmd3_load(S, P0, 0.0, err, sizeof err)) { fprintf(stderr, "load: %s\n", err); exit(2); }
    if (audits) { edmd3_set_contact_audit(S, 1); edmd3_set_schedule_audit(S, full > 0 ? 1 : every); }
    int nprim = edmd3_wall_primitives(S, NULL, NULL);
    r->ke0 = edmd3_kinetic_energy(S);
    const double teq = 0.1 * T;
    r->nb = blocks; r->tb = blocks ? (T - teq) / blocks : 0.0;
    double Jlast = 0.0; int ib = 0;
    const long nstop = (long)(T / SIGT + 0.5);
    for (long k = 1; k <= nstop; ++k) {
        const double t = k * SIGT;
        edmd3_advance_to(S, t);
        const char* msg; if (edmd3_fatal(S, &msg)) { fprintf(stderr, "fatal: %s\n", msg); r->fatal = 1; break; }
        const EDMD3_Health* H = edmd3_health(S);
        if (audits && full > 0 && H->ev_pair + H->ev_wall + H->ev_seg + H->ev_arc >= full) { edmd3_set_schedule_audit(S, every); full = 0; }
        const EDMD_Particle* P = edmd3_particles(S);
        for (int i = 0; i < N; ++i) {
            const double g = edmd3_prim_gap(S, P[i].x, P[i].y);
            if (g < r->min_gap) r->min_gap = g;
            if (g < -1e-9) r->outside++;
        }
        r->stops++;
        if (blocks) {                         /* the walls' normal impulse on their lines and arcs; the blocks after T/10 */
            double J = 0.0;
            for (int p = 0; p < nprim; ++p) J += edmd3_prim_impulse(S, p, 0, NULL);
            if (t <= teq + 1e-9) Jlast = J;
            else if (ib < blocks && t >= teq + (ib + 1) * r->tb - 1e-9) { r->Jb[ib++] = J - Jlast; Jlast = J; }
        }
    }
    for (int p = 0; p < nprim; ++p) r->Jend += edmd3_prim_impulse(S, p, 1, NULL);
    if (audits) edmd3_schedule_audit_now(S);
    r->hash = edmd3_event_hash(S); r->H = *edmd3_health(S); r->A = *edmd3_schedule_audit_stats(S);
    r->cev = edmd3_contact_audit_prims(S, r->cmax); edmd3_tolerances(S, &r->tol);
    r->ke1 = edmd3_kinetic_energy(S); r->clean = edmd3_health_clean(S);
    r->fin = (EDMD_Particle*)malloc((size_t)N * sizeof(EDMD_Particle));
    memcpy(r->fin, edmd3_particles(S), (size_t)N * sizeof(EDMD_Particle));
    edmd3_destroy(S);
}
static int same_disks(const EDMD_Particle* a, const EDMD_Particle* b, int N){
    for (int i = 0; i < N; ++i)
        if (memcmp(&a[i].x, &b[i].x, sizeof(double)) || memcmp(&a[i].y, &b[i].y, sizeof(double)) || memcmp(&a[i].vx, &b[i].vx, sizeof(double))
            || memcmp(&a[i].vy, &b[i].vy, sizeof(double)) || a[i].coll_count != b[i].coll_count) return 0;
    return 1;
}
/* A, B, C and the checks every geometry shares (test 4's list without the pressure; B1: gaps against tol_face) */
static void abc(const char* tag, const Geo* G, const EDMD_Particle* P0, int N, double T, long full, long every, int blocks, GRun* B){
    GRun A, C;
    run(G, P0, N, T, 1, full, every, blocks, &A);
    run(G, P0, N, T, 0, 0, 0, blocks, B);
    run(G, P0, N, T, 0, 0, 0, blocks, &C);
    const EDMD3_Audit* a = &A.A;
    printf("%s N %d T %.0f sigma-time: hash A %016llx B %016llx C %016llx; events seg %ld arc %ld pair %ld cross %ld\n", tag, N, T / SIGT,
           (unsigned long long)A.hash, (unsigned long long)B->hash, (unsigned long long)C.hash, B->H.ev_seg, B->H.ev_arc, B->H.ev_pair, B->H.ev_cross);
    char w[256];
    snprintf(w, sizeof w, "%s audits do not steer (A = B: hash, health, disks)", tag);
    check(A.hash == B->hash && !memcmp(&A.H, &B->H, sizeof A.H) && same_disks(A.fin, B->fin, N), w);
    snprintf(w, sizeof w, "%s deterministic (B = C: hash, health, disks)", tag);
    check(B->hash == C.hash && !memcmp(&B->H, &C.H, sizeof C.H) && same_disks(B->fin, C.fin, N), w);
    printf("%s run flags A %d B %d C %d; fatal %d %d %d; repairs: wall overlap %ld, findings %ld, wall repair %ld, receding %ld; contact_now %ld\n", tag,
           A.clean, B->clean, C.clean, A.fatal, B->fatal, C.fatal, B->H.prim_overlap_repair, B->H.prim_findings, B->H.prim_repair, B->H.prim_receding, B->H.prim_contact_now);
    snprintf(w, sizeof w, "%s clean = 1 in A, B and C", tag); check(A.clean && B->clean && C.clean && !A.fatal && !B->fatal && !C.fatal, w);
    const double drift = (B->ke1 - B->ke0) / B->ke0;
    printf("%s relative KE drift %.3e (B)\n", tag, drift);
    snprintf(w, sizeof w, "%s |KE drift| <= 1e-12", tag); check(fabs(drift) <= 1e-12, w);
    const double mg = fmin(A.min_gap, fmin(B->min_gap, C.min_gap));
    printf("%s every stop (every sigma-time, %ld stops x 3 runs): smallest signed centre-to-wall gap %.3e px; centres outside the lumen %ld\n", tag,
           B->stops, mg, A.outside + B->outside + C.outside);
    snprintf(w, sizeof w, "%s no centre closer than R - 1e-9 px to a wall, none outside", tag); check(mg >= -1e-9 && A.outside + B->outside + C.outside == 0, w);
    printf("%s schedule audit (A): %ld audits; walls compared %ld, missing %ld, extra %ld, |dt| > 1e-9 %ld, of them > 1e-10 x horizon %ld, deferred %ld, "
           "deferred early %ld, boundary %ld, max |dt| %.3e; all classes: pair missing %ld extra %ld, wall missing %ld extra %ld, cross missing %ld extra %ld dup %ld; "
           "duplicate disagreements %ld\n", tag, a->audits, a->prim_cmp, a->prim_missing, a->prim_extra, a->prim_dt, a->prim_dt_rel, a->prim_deferred,
           a->prim_deferred_early, a->prim_boundary, a->prim_max_dt, a->pair_missing, a->pair_extra, a->wall_missing, a->wall_extra, a->cross_missing,
           a->cross_extra, a->cross_dup, a->dup_disagree);
    snprintf(w, sizeof w, "%s schedule audit: 0 missing, 0 extra, 0 early, 0 off by > 1e-10 x horizon (walls and every class)", tag);
    check(a->audits > 0 && a->prim_cmp > 0 && !a->prim_missing && !a->prim_extra && !a->prim_dt_rel && !a->prim_deferred_early && !a->pair_missing && !a->pair_extra
          && !a->wall_missing && !a->wall_extra && !a->cross_missing && !a->cross_extra && !a->cross_dup && !a->dup_disagree && !a->pair_dt_rel && !a->cross_dt_rel
          && !a->pair_deferred_early, w);
    printf("%s contact audit (A, %ld events): max |gap| segments %.3e px, arcs %.3e px; tol_face %.3e px (pre-check B1: gaps against it, not 1e-12)\n",
           tag, A.cev, A.cmax[0], A.cmax[1], A.tol.tol_face);
    snprintf(w, sizeof w, "%s contact-audit gaps <= tol_face", tag); check(A.cmax[0] <= A.tol.tol_face && A.cmax[1] <= A.tol.tol_face, w);
    free(A.fin); free(C.fin);
}

/* ------------------------------------------------------------------ the modes */
static double zblocks(const Geo* G, const GRun* r, int N, double kT, double* se, double* zfull){
    double z[10], m = 0.0, v = 0.0, mf = 0.0;
    for (int b = 0; b < r->nb; ++b) {
        const double P = r->Jb[b] / (r->tb * G->perim_acc);           /* force per length on the accessible contact lines */
        z[b] = P * G->area_acc / (N * kT); m += z[b];
        mf += (r->Jb[b] / (r->tb * 2.0 * (G->Lx + G->Ly))) * G->area / (N * kT);
    }
    m /= r->nb; mf /= r->nb;
    for (int b = 0; b < r->nb; ++b) v += (z[b] - m) * (z[b] - m);
    *se = sqrt(v / (r->nb - 1) / r->nb); *zfull = mf;
    return m;
}
/* returns bit 0: Z at the accessible area within 3 % of KR there; bit 1: Z with the full area within 3 % of KR (sec. 8 as written) */
static int mode_corridor(double eta, double tilt, double Ts, uint64_t sd, double* zout, double* seout){
    Geo G; make_corridor(&G, 100.0 * PX, 40.0 * PX, tilt * PI / 180.0);
    const int N = (int)lround(eta * G.area / (PI * R * R));
    EDMD_Particle* P0 = seed(&G, N, sd);
    char tag[96]; snprintf(tag, sizeof tag, "corridor eta %.2f tilt %.0f:", eta, tilt);
    printf("%s rectangle 100 x 40 sigma at %.0f degrees, area %.6f sigma^2 (Green), N = %d, eta = %.5f, eta of the accessible area = %.5f\n", tag, tilt,
           G.area / (PX * PX), N, N * PI * R * R / G.area, N * PI * R * R / G.area_acc);
    GRun B;
    const long full = N > 1000 ? 200 : 2000, every = N > 1000 ? 2000 : 500;
    abc(tag, &G, P0, N, Ts * SIGT, full, every, 10, &B);
    const double kT = B.ke0 / N;
    double se, zf; const double z = zblocks(&G, &B, N, kT, &se, &zf);
    const double eacc = N * PI * R * R / G.area_acc, efull = N * PI * R * R / G.area;
    printf("%s pressure on the walls (B, 10 blocks after T/10): Z = P A_acc / (N kT) = %.5f +- %.5f (block SE); KR(eta_acc = %.5f) = %.5f, deviation %+.2f %%; "
           "with the full area and wall length (brief sec. 8 as written): %.5f vs KR(%.5f) = %.5f, %+.2f %%; endpoint impulse %.3g\n", tag, z, se, eacc,
           z_kr(eacc), 100.0 * (z / z_kr(eacc) - 1.0), zf, efull, z_kr(efull), 100.0 * (zf / z_kr(efull) - 1.0), B.Jend);
    *zout = z; *seout = se;
    free(P0); free(B.fin);
    return (100.0 * fabs(z / z_kr(eacc) - 1.0) <= 3.0) | (100.0 * fabs(zf / z_kr(efull) - 1.0) <= 3.0) << 1;
}
static void mode_y(double eta, double Ts, uint64_t sd){
    Geo G; make_y(&G, 5.0 * PX, 30.0 * PX, 0);
    const int N = (int)lround(eta * G.area / (PI * R * R));
    EDMD_Particle* P0 = seed(&G, N, sd);
    char tag[64]; snprintf(tag, sizeof tag, "y eta %.2f:", eta);
    printf("%s 6 segments, 3 arcs, area %.6f sigma^2 (Green), N = %d\n", tag, G.area / (PX * PX), N);
    GRun B; abc(tag, &G, P0, N, Ts * SIGT, 2000, 500, 0, &B);
    free(P0); free(B.fin);
}
static void mode_yaim(int split, double Ts, uint64_t sd){
    Geo G; make_y(&G, 5.0 * PX, 30.0 * PX, split);
    const int Nr = 12, N = Nr + 4;
    EDMD_Particle* P0 = seed(&G, N, sd);      /* 16 random disks; the first four are then placed by hand */
    /* disk 0 and 1: towards reflex vertices 0 and 1 along their bisector, from 50 px out; disk 2: perpendicular onto the joint;
       disk 3: towards vertex 2 */
    const double c0 = 0.5 * BOXW, c1 = 0.5 * BOXH;
    for (int k = 0; k < 3; ++k) {
        const int i = k == 2 ? 3 : k;
        const double wx = G.vx_[k] - c0, wy = G.vy_[k] - c1, wn = hypot(wx, wy);
        P0[i].x = G.vx_[k] - (R + 50.0) * wx / wn; P0[i].y = G.vy_[k] - (R + 50.0) * wy / wn; P0[i].vx = wx / wn; P0[i].vy = wy / wn;
    }
    P0[2].x = G.jx + (R + 50.0) * G.jnx; P0[2].y = G.jy + (R + 50.0) * G.jny; P0[2].vx = -G.jnx; P0[2].vy = -G.jny;
    /* move any random disk that overlaps the four, inside the cap of arm 1 */
    for (int i = 4; i < N; ++i)
        for (int j = 0; j < 4; ++j) if (hypot(P0[i].x - P0[j].x, P0[i].y - P0[j].y) < 2 * R + 1.0) { P0[i].x = 0.0; break; }
    int m = 4;
    for (int i = 4; i < N; ++i) if (P0[i].x != 0.0) P0[m++] = P0[i];
    char tag[64]; snprintf(tag, sizeof tag, "yaim split %d:", split);
    printf("%s 4 disks aimed at reflex vertices 0, 1, 2 and at an arc-segment joint, %d random; the schedule audit after every event\n", tag, m - 4);
    GRun B; abc(tag, &G, P0, m, Ts * SIGT, 1L << 40, 1, 0, &B);
    check(B.H.stagnation == 0, "yaim no stagnation");
    free(P0); free(B.fin);
}
static void mode_empty(double Ts, uint64_t sd){
    Geo G; make_corridor(&G, 100.0 * PX, 40.0 * PX, 0.0);
    const int N = 400;
    EDMD_Particle* P0 = seed(&G, N, sd);      /* seeded inside the corridor, but run in the plain box */
    uint64_t h[2]; EDMD3_Health H[2]; EDMD_Particle* F[2];
    for (int w = 0; w < 2; ++w) {
        char err[512]; EDMD_Params p = box_params(N);
        EDMD3* S = edmd3_create(&p, 0.0, err, sizeof err);
        if (w && !edmd3_set_walls(S, NULL, 0, NULL, 0, err, sizeof err)) { fprintf(stderr, "%s\n", err); exit(2); }
        if (!edmd3_load(S, P0, 0.0, err, sizeof err)) { fprintf(stderr, "%s\n", err); exit(2); }
        edmd3_advance_to(S, Ts * SIGT);
        h[w] = edmd3_event_hash(S); H[w] = *edmd3_health(S);
        F[w] = (EDMD_Particle*)malloc((size_t)N * sizeof(EDMD_Particle)); memcpy(F[w], edmd3_particles(S), (size_t)N * sizeof(EDMD_Particle));
        edmd3_destroy(S);
    }
    printf("empty: box run of %d disks, %.0f sigma-time: hash without the call %016llx, with an empty edmd3_set_walls %016llx\n", N, Ts,
           (unsigned long long)h[0], (unsigned long long)h[1]);
    check(h[0] == h[1] && !memcmp(&H[0], &H[1], sizeof H[0]) && same_disks(F[0], F[1], N), "empty: an empty wall list changes nothing (hash, health, disks)");
    free(P0); free(F[0]); free(F[1]);
}

int main(int argc, char** argv){
    if (argc < 2) { fprintf(stderr, "usage: %s corridor ETA TILT T SEED | y ETA T SEED | yaim SPLIT T SEED | empty T SEED | pair ETA T SEED\n", argv[0]); return 2; }
    const char* m = argv[1];
    if (!strcmp(m, "corridor") && argc == 6) {
        double z, se; const int ok = mode_corridor(atof(argv[2]), atof(argv[3]), atof(argv[4]), strtoull(argv[5], NULL, 0), &z, &se);
        printf("corridor result Z %.6f SE %.6f\n", z, se);
        if (atof(argv[2]) < 0.1) { check(ok & 1, "corridor: Z within 3 % of KR at the accessible packing");
                                   check(ok >> 1, "corridor: Z with the full area within 3 % of KR (sec. 8 as written)"); }
    } else if (!strcmp(m, "pair") && argc == 5) {   /* tests 4 and 5 together: straight and tilted at one packing, Z compared */
        const double eta = atof(argv[2]), T = atof(argv[3]); const uint64_t sd = strtoull(argv[4], NULL, 0);
        double z0, s0, z1, s1;
        const int ok0 = mode_corridor(eta, 0.0, T, sd, &z0, &s0), ok1 = mode_corridor(eta, 30.0, T, sd, &z1, &s1);
        const double dz = z1 - z0, sc = sqrt(s0 * s0 + s1 * s1);
        printf("pair eta %.2f: Z straight %.5f +- %.5f, tilted %.5f +- %.5f, difference %+.5f = %.2f combined SE\n", eta, z0, s0, z1, s1, dz, fabs(dz) / sc);
        check(fabs(dz) <= 2.0 * sc, "pair: tilted Z within 2 block SE of straight");
        if (eta < 0.1) {
            check(ok0 & 1, "pair: straight Z within 3 % of KR at the accessible packing"); check(ok1 & 1, "pair: tilted Z within 3 % of KR at the accessible packing");
            check(ok0 >> 1, "pair: straight Z with the full area within 3 % of KR (sec. 8 as written)");
            check(ok1 >> 1, "pair: tilted Z with the full area within 3 % of KR (sec. 8 as written)");
        }
        else printf("pair eta %.2f: the KR deviations are reported, not gated (02_DEMO1_BRIEF.md sec. 8)\n", eta);
    } else if (!strcmp(m, "y") && argc == 5) mode_y(atof(argv[2]), atof(argv[3]), strtoull(argv[4], NULL, 0));
    else if (!strcmp(m, "yaim") && argc == 5) mode_yaim(atoi(argv[2]), atof(argv[3]), strtoull(argv[4], NULL, 0));
    else if (!strcmp(m, "empty") && argc == 4) mode_empty(atof(argv[2]), strtoull(argv[3], NULL, 0));
    else { fprintf(stderr, "bad arguments\n"); return 2; }
    printf("failed checks: %d\n", nfail);
    if (fflush(stdout) != 0 || ferror(stdout)) return 1;
    return nfail ? 1 : 0;
}
