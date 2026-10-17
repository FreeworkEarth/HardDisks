/* PID Demo 1b, test 10 (02_DEMO1_BRIEF.md sec. 8): positive controls. A detector that cannot fail proves nothing, so four
   errors are planted in a white-box build of the fork (this file includes engine/edmd_gen3.c, compiled with -DPID_WHITEBOX, as
   hspist3/edmd_core/tests/gen3_band_edge_ties.c includes the original), and each must be reported by its own detector, while
   the same state without the error reports nothing:
     (i)   a live wall event deleted from the heap      -> the schedule audit counts it missing (prim_missing >= 1)
     (ii)  a disk moved 0.1 R through a segment         -> the validator counts it (prim_findings >= 1) and the run flag is 0
     (iii) the reflection's sign flipped (v + 2(v.n)n)  -> the kinetic energy drifts and the run flag is 0
     (iv)  a live wall event moved by 1e-3 time units   -> the schedule audit counts it off (prim_dt_rel >= 1); its baseline
           is a disk grazing the wall of the corridor turned by 30 degrees (gap 1e-3 px, v_n 1e-4 px per unit), whose heap
           and brute-force times differ by far more than 1e-10 of the horizon but by v_n |dt| ~ 1e-13 px along the normal:
           the same contact to rounding, counted apart (prim_dt_cond >= 1) and not as an error (prim_dt_rel = 0)
   The geometry is the corridor of tests 4 and 5 (a closed 100 x 40 sigma rectangle in the 200 x 120 sigma box), 300 disks.
   Seeding: Demo 0's one generator (demo0/pid_rng.c, linked unchanged; brief sec. 2 rule 5).
   Output: "control ..." lines and "check ... OK|FAIL"; exit 0 iff no FAIL. */
#define PID_WHITEBOX 1
#include "../engine/edmd_gen3.c"
#include "pid.h"          /* demo0/: the generator only */

static const double PXS = 24.0, RR = 12.0, BW = 4800.0, BH = 2880.0, PIV = 3.14159265358979323846;
static int nfail = 0;
static void check(int ok, const char* what){ printf("check %s %s\n", what, ok ? "OK" : "FAIL"); if (!ok) nfail++; }
static Rng rng;
static double ur(void){ return rng_uniform(&rng); }

static EDMD3_Seg SEG[4];
static void corridor(void){
    const double cx = 0.5 * BW, cy = 0.5 * BH, a = 50.0 * PXS, b = 20.0 * PXS;
    const double X[4] = {cx - a, cx + a, cx + a, cx - a}, Y[4] = {cy - b, cy - b, cy + b, cy + b};
    for (int k = 0; k < 4; ++k) { SEG[k].x1 = X[k]; SEG[k].y1 = Y[k]; SEG[k].x2 = X[(k + 1) % 4]; SEG[k].y2 = Y[(k + 1) % 4]; }
}
/* a state of N disks seeded inside the corridor (RSA, MB velocities), walls set and loaded */
static EDMD3* state(int N){
    EDMD_Params p; memset(&p, 0, sizeof p);
    p.boxW = BW; p.boxH = BH; p.radius = RR; p.N = N; p.pp_collisions_enabled = 1; p.particle_mass = 1.0; p.kB = 1.0;
    char err[512];
    EDMD3* S = edmd3_create(&p, 0.0, err, sizeof err);
    if (!S || !edmd3_set_walls(S, SEG, 4, NULL, 0, err, sizeof err)) { fprintf(stderr, "%s\n", err); exit(2); }
    EDMD_Particle* P = (EDMD_Particle*)calloc((size_t)N, sizeof(EDMD_Particle));
    rng_seed(&rng, 0x5EED10ULL);
    for (int i = 0; i < N; ++i)
        for (;;) {
            const double x = SEG[0].x1 + RR + 1e-3 + ur() * (100.0 * PXS - 2 * RR - 2e-3), y = SEG[0].y1 + RR + 1e-3 + ur() * (40.0 * PXS - 2 * RR - 2e-3);
            int ok = 1;
            for (int j = 0; j < i && ok; ++j) if ((x - P[j].x) * (x - P[j].x) + (y - P[j].y) * (y - P[j].y) < 4 * RR * RR + 1e-6) ok = 0;
            if (ok) { P[i].x = x; P[i].y = y; break; }
        }
    for (int i = 0; i < N; ++i) { P[i].vx = rng_normal(&rng); P[i].vy = rng_normal(&rng); }
    if (!edmd3_load(S, P, 0.0, err, sizeof err)) { fprintf(stderr, "%s\n", err); exit(2); }
    free(P);
    return S;
}
/* (iv): the corridor turned by 30 degrees about its centre (counter-clockwise), N disks seeded in its frame; disk 0 grazes the
   bottom wall (segment 0) at its middle: gap 1e-3 px, closing at v_n = 1e-4 px per unit, tangential speed 1 */
static EDMD3* state_tilted_graze(int N){
    const double th = PIV / 6.0, ux = cos(th), uy = sin(th), cx = 0.5 * BW, cy = 0.5 * BH, a = 50.0 * PXS, b = 20.0 * PXS;
    const double sa[4] = {-1, 1, 1, -1}, sb[4] = {-1, -1, 1, 1};
    double X[4], Y[4];
    for (int k = 0; k < 4; ++k) { X[k] = cx + sa[k] * a * ux - sb[k] * b * uy; Y[k] = cy + sa[k] * a * uy + sb[k] * b * ux; }
    for (int k = 0; k < 4; ++k) { SEG[k].x1 = X[k]; SEG[k].y1 = Y[k]; SEG[k].x2 = X[(k + 1) % 4]; SEG[k].y2 = Y[(k + 1) % 4]; }
    EDMD_Params p; memset(&p, 0, sizeof p);
    p.boxW = BW; p.boxH = BH; p.radius = RR; p.N = N; p.pp_collisions_enabled = 1; p.particle_mass = 1.0; p.kB = 1.0;
    char err[512];
    EDMD3* S = edmd3_create(&p, 0.0, err, sizeof err);
    if (!S || !edmd3_set_walls(S, SEG, 4, NULL, 0, err, sizeof err)) { fprintf(stderr, "%s\n", err); exit(2); }
    EDMD_Particle* P = (EDMD_Particle*)calloc((size_t)N, sizeof(EDMD_Particle));
    const double g = 1e-3, vn = 1e-4, nx = -uy, ny = ux;              /* the bottom wall's normal into the lumen */
    P[0].x = cx - (b - RR - g) * nx; P[0].y = cy - (b - RR - g) * ny;
    P[0].vx = ux - vn * nx; P[0].vy = uy - vn * ny;
    rng_seed(&rng, 0x5EED11ULL);
    for (int i = 1; i < N; ++i)
        for (;;) {
            const double s = (2.0 * ur() - 1.0) * (a - RR - 1e-3), t = (2.0 * ur() - 1.0) * (b - RR - 1e-3);
            const double x = cx + s * ux - t * uy, y = cy + s * uy + t * ux;
            int ok = 1;
            for (int j = 0; j < i && ok; ++j) if ((x - P[j].x) * (x - P[j].x) + (y - P[j].y) * (y - P[j].y) < 4 * RR * RR + 1e-6) ok = 0;
            if (ok) { P[i].x = x; P[i].y = y; break; }
        }
    for (int i = 1; i < N; ++i) { P[i].vx = rng_normal(&rng); P[i].vy = rng_normal(&rng); }
    if (!edmd3_load(S, P, 0.0, err, sizeof err)) { fprintf(stderr, "%s\n", err); exit(2); }
    free(P);
    return S;
}
static int in_list(const EDMD3* S, int c, int p){ for (int k = S->pcs[c]; k < S->pcs[c + 1]; ++k) if (S->pcl[k] == p) return 1; return 0; }

int main(void){
    corridor();
    /* (i) a deleted wall event */
    {
        EDMD3* S = state(300);
        edmd3_advance_to(S, 5.0 * PXS);
        schedule_audit(S);
        const long m0 = S->A.prim_missing, x0 = S->A.prim_extra;
        long q = -1;
        for (long k = 0; k < S->heap.n && q < 0; ++k) {
            const Ev3* e = &S->heap.d[k];
            if (e->type == T_SEG && ev_live(S, e) && e->t > S->now && in_list(S, cell_of(S, S->D[e->a].cx, S->D[e->a].cy), e->b)) q = k;
        }
        if (q < 0) { printf("control (i): no deletable wall event found\n"); check(0, "control (i) setup"); }
        else {
            const Ev3 del = S->heap.d[q];
            S->heap.d[q] = S->heap.d[--S->heap.n];
            for (long i = S->heap.n / 2 - 1; i >= 0; --i) heap_sift_down(&S->heap, i);
            schedule_audit(S);
            printf("control (i): the live T_SEG event of disk %d and segment %d at t = %.6g deleted; schedule audit missing %ld -> %ld (extra %ld -> %ld)\n",
                   del.a, del.b, (double)(S->T0 + del.t), m0, S->A.prim_missing - m0, x0, S->A.prim_extra - x0);
            check(m0 == 0 && x0 == 0, "control (i) baseline: the audit finds nothing missing or extra");
            check(S->A.prim_missing - m0 >= 1, "control (i) the deleted event is reported missing");
        }
        edmd3_destroy(S);
    }
    /* (ii) a disk 0.1 R through a segment */
    {
        EDMD3* S = state(300);
        edmd3_advance_to(S, 5.0 * PXS);
        sync_all_internal(S);                 /* every disk at now (this runs the validator once: the baseline) */
        const long f0 = S->H.prim_findings; const int c0 = edmd3_health_clean(S);
        const int i = 0; Disk3* A = &S->D[i];
        const double x = 0.5 * (SEG[0].x1 + SEG[0].x2), y = SEG[0].y1 - 0.1 * RR;   /* below the bottom wall's line by 0.1 R */
        cell_remove(S, i);
        A->cx = (int)floor(x / S->w); A->cy = (int)floor(y / S->w); A->xi = x - A->cx * S->w; A->zeta = y - A->cy * S->w; A->tau = S->now;
        cell_insert(S, i);
        full_check(S);
        printf("control (ii): disk %d put 0.1 R beyond segment 0 (centre at %.3f, %.3f px); validator findings %ld -> %ld, run flag %d -> %d\n",
               i, x, y, f0, S->H.prim_findings, c0, edmd3_health_clean(S));
        check(f0 == 0 && c0 == 1, "control (ii) baseline: no finding, run flag 1");
        check(S->H.prim_findings - f0 >= 1 && edmd3_health_clean(S) == 0, "control (ii) the validator reports it and the run flag is 0");
        edmd3_destroy(S);
    }
    /* (iii) the reflection's sign flipped */
    for (int plant = 0; plant < 2; ++plant) {
        pid_plant_reflect_sign = plant;
        EDMD3* S = state(300);
        edmd3_set_contact_audit(S, 1);
        const double ke0 = edmd3_kinetic_energy(S);
        edmd3_advance_to(S, 20.0 * PXS);
        const char* msg; const int fat = edmd3_fatal(S, &msg);
        const double drift = (edmd3_kinetic_energy(S) - ke0) / ke0;
        double cm[2]; edmd3_contact_audit_prims(S, cm);
        printf("control (iii) %s: KE drift %.3e, run flag %d, fatal %d, wall events %ld, wall overlap repairs %ld, findings %ld, contact-audit max gap %.3e px\n",
               plant ? "planted" : "baseline", drift, edmd3_health_clean(S), fat, S->H.ev_seg, S->H.prim_overlap_repair, S->H.prim_findings, cm[0]);
        if (!plant) check(fabs(drift) <= 1e-12 && edmd3_health_clean(S) == 1 && S->H.ev_seg > 0, "control (iii) baseline: KE conserved, run flag 1");
        else check(fabs(drift) > 1e-6 && edmd3_health_clean(S) == 0, "control (iii) the flipped reflection is reported: KE drift and run flag 0");
        edmd3_destroy(S);
    }
    pid_plant_reflect_sign = 0;
    /* (iv) the time test of the schedule audit: a grazing disk is not an error, a moved event is */
    {
        EDMD3* S = state_tilted_graze(300);
        edmd3_set_schedule_audit(S, 1);
        edmd3_advance_to(S, 12.0);            /* the grazing disk reaches the wall at t = 10 */
        const long r0 = S->A.prim_dt_rel, c0 = S->A.prim_dt_cond;
        printf("control (iv) baseline: corridor turned by 30 degrees, disk 0 grazing segment 0 (gap 1e-3 px, v_n 1e-4): %ld audits, wall time "
               "disagreements > 1e-10 x horizon: %ld grazing (the same contact to rounding), %ld errors; max |dt| %.3e, max rel %.3e\n",
               S->A.audits, c0, r0, (double)S->A.prim_max_dt, (double)S->A.prim_max_rel);
        check(r0 == 0 && c0 >= 1, "control (iv) baseline: a grazing approach agrees to rounding (counted apart, not an error)");
        edmd3_set_schedule_audit(S, 0);
        long q = -1;
        for (long k = 0; k < S->heap.n && q < 0; ++k) {
            const Ev3* e = &S->heap.d[k];
            if (e->type != T_SEG || !ev_live(S, e) || e->t <= S->now || e->cb != 0) continue;
            const real* G = &S->pseg[(long)e->b * PS_N];
            if (R_(fabs)(S->D[e->a].vx * G[PS_NX] + S->D[e->a].vy * G[PS_NY]) > 0.1) q = k;
        }
        if (q < 0) { printf("control (iv): no live line event with v_n > 0.1 found\n"); check(0, "control (iv) setup"); }
        else {
            const Ev3 mv = S->heap.d[q];
            S->heap.d[q].t += 1e-3;
            for (long i = S->heap.n / 2 - 1; i >= 0; --i) heap_sift_down(&S->heap, i);
            schedule_audit(S);
            printf("control (iv): the live T_SEG event of disk %d and segment %d at t = %.6g moved by 1e-3; schedule audit errors %ld -> %ld\n",
                   mv.a, mv.b, (double)(S->T0 + mv.t), r0, S->A.prim_dt_rel);
            check(S->A.prim_dt_rel - r0 >= 1, "control (iv) the moved wall event is reported (prim_dt_rel)");
        }
        edmd3_destroy(S);
    }
    printf("failed checks: %d\n", nfail);
    return nfail ? 1 : 0;
}
