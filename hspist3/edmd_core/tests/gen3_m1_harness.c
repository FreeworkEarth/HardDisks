/* ##CHRIS 2026-10-08 (261012 sec. 4.7.1, item 5; milestone M1): the M1 harness for the generation-3 core.
   For each cell it builds one initial state and runs, from that same state:
     A  gen3 with both audits: the schedule audit (live heap against brute force from the synchronised state) after EVERY
        event for the first --full events, then after every --every-th event; the contact audit at every event;
     B  gen3 without audits            -> event hash and final state equal to A: the audits do not steer;
     C  gen3 again without audits      -> equal to B: deterministic;
     D  gen2 (edmd.c, the code of the binary's legacy path) with its contact audit (HD_CONTACT_AUDIT);
   and prints the health counters, audit counts, energy drift, the virial pressure in time blocks, and events per second.
   Cells: fluid (pi/8, random sequential insertion), dense (0.70, jittered triangular lattice), lattice (0.70, the same
   lattice WITHOUT jitter), tie (a square lattice on the cell corners with equal speeds: every collision of a row happens
   at exactly the same time, the tie-break stress test of sec. 4.7.1, d). Then the speed table: events per second of
   gen2 and gen3 at N = 400 and 1600, pi/8 and 0.70, from the same states. Then the divergence of gen2 and gen3 from one
   state (rounding amplified by chaos).
   Units as in the engine: px, internal time (24 px = 1 sigma, 24 units = 1 sigma-time); unit mass, kT = 1.
   build (from hspist3/): cc -std=c11 -O3 -ffp-contract=off -Wall -Wextra -o <out> edmd_core/tests/gen3_m1_harness.c
                          edmd_core/edmd_gen3.c edmd_core/edmd.c -lm
   usage: <out> audit|speed|diverge [--quick]   (separate processes: edmd.c reads HD_CONTACT_AUDIT once per process, so
          the gen2 speed is measured in a process whose gen2 never had the contact audit on) */
#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include "../edmd.h"
#include "../edmd_gen3.h"

#define PX 24.0            /* px per sigma */
#define SIGT 24.0          /* internal units per sigma-time */
static const double R = 12.0;
static const double PI = 3.14159265358979323846;

static uint64_t rng_state;
static uint64_t xs64(void){ uint64_t x = rng_state; x ^= x << 13; x ^= x >> 7; x ^= x << 17; return rng_state = x; }
static double urand(void){ return ((double)(xs64() >> 11) + 0.5) / 9007199254740992.0; }
static double now_s(void){ struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + 1e-9 * ts.tv_nsec; }

typedef struct { const char* name; int N; double eta, W, H; EDMD_Particle* P; } Cell;

static void gauss_velocities(EDMD_Particle* P, int N){
    for (int i = 0; i < N; ++i) {
        const double u1 = urand(), u2 = urand(), g = sqrt(-2.0 * log(u1)), th = 2.0 * PI * u2;
        P[i].vx = g * cos(th); P[i].vy = g * sin(th); P[i].coll_count = 0;
    }
}
static double box_eta(int N, double W, double H){ return N * PI * R * R / (W * H); }

/* random sequential insertion (fluid only) */
static int make_rsa(Cell* c, uint64_t seed){
    rng_state = seed;
    c->P = (EDMD_Particle*)calloc((size_t)c->N, sizeof(EDMD_Particle));
    for (int i = 0; i < c->N; ++i) {
        for (long tries = 0;; ++tries) {
            if (tries > 10000000L) return 0;
            const double x = R + urand() * (c->W - 2 * R), y = R + urand() * (c->H - 2 * R);
            int ok = 1;
            for (int j = 0; j < i && ok; ++j) { const double dx = x - c->P[j].x, dy = y - c->P[j].y; if (dx * dx + dy * dy < 4 * R * R) ok = 0; }
            if (ok) { c->P[i].x = x; c->P[i].y = y; break; }
        }
    }
    gauss_velocities(c->P, c->N);
    return 1;
}
/* triangular lattice, nr rows x nc columns, spread over the box with a margin m from the walls; jitter amp in units of
   half the smallest free gap (0 = none) */
static int make_hex(Cell* c, int nc, int nr, double jitter, uint64_t seed){
    rng_state = seed;
    c->N = nc * nr;
    c->P = (EDMD_Particle*)calloc((size_t)c->N, sizeof(EDMD_Particle));
    const double m = 0.5;
    const double ax = (c->W - 2 * R - 2 * m) / (nc - 0.5), ay = (c->H - 2 * R - 2 * m) / (nr - 1);
    const double dmin = fmin(ax, sqrt(0.25 * ax * ax + ay * ay));
    if (!(dmin > 2 * R)) return 0;
    const double amp = jitter * 0.5 * (dmin - 2 * R);
    int k = 0;
    for (int r = 0; r < nr; ++r) for (int q = 0; q < nc; ++q, ++k) {
        double x = R + m + q * ax + ((r & 1) ? 0.5 * ax : 0.0), y = R + m + r * ay;
        if (amp > 0.0) { x += (2 * urand() - 1) * amp; y += (2 * urand() - 1) * amp; }
        c->P[k].x = x; c->P[k].y = y;
    }
    gauss_velocities(c->P, c->N);
    return 1;
}
/* square lattice with every site on a cell corner (multiples of 32 px) and speed 1 along x, alternating in each row:
   neighbours (0,1), (2,3), ... meet head-on at exactly t = 4, and so on: exact ties every few units */
static void make_tie(Cell* c, int n){
    c->N = n * n; c->W = 32.0 * (n + 1); c->H = 32.0 * (n + 1);
    c->P = (EDMD_Particle*)calloc((size_t)c->N, sizeof(EDMD_Particle));
    for (int r = 0, k = 0; r < n; ++r) for (int q = 0; q < n; ++q, ++k) {
        c->P[k].x = 32.0 * (q + 1); c->P[k].y = 32.0 * (r + 1);
        c->P[k].vx = (q & 1) ? -1.0 : 1.0; c->P[k].vy = 0.0; c->P[k].coll_count = 0;
    }
    c->eta = box_eta(c->N, c->W, c->H);
}

static EDMD_Params params(const Cell* c){
    EDMD_Params p; memset(&p, 0, sizeof p);
    p.boxW = c->W; p.boxH = c->H; p.radius = R; p.N = c->N; p.pp_collisions_enabled = 1;
    p.particle_mass = 1.0; p.kB = 1.0;
    return p;
}

static double g_cell_px = 0.0;     /* 0 = the engine default (32 px) */

typedef struct {
    uint64_t hash; double secs; long ev_phys, ev_cross; double ke0, ke1, Zmean, Zse; int nblk;
    EDMD3_Health H; EDMD3_Audit A; double cmax[2]; long cev; int fatal; char msg[256];
    EDMD_Particle* fin;
} R3;

/* one gen3 run from c->P to T (internal units) in steps of 1 sigma-time; virial in blocks of Tblk after T/10 */
static int run3(const Cell* c, double T, int audits, long full, long every, R3* r){
    char err[256]; memset(r, 0, sizeof *r);
    EDMD_Params p = params(c);
    EDMD3* S = edmd3_create(&p, g_cell_px, err, sizeof err);
    if (!S) { fprintf(stderr, "%s\n", err); return 0; }
    if (!edmd3_load(S, c->P, 0.0, err, sizeof err)) { fprintf(stderr, "%s: %s\n", c->name, err); edmd3_destroy(S); return 0; }
    if (audits) { edmd3_set_contact_audit(S, 1); edmd3_set_schedule_audit(S, 1); }
    r->ke0 = edmd3_kinetic_energy(S);
    const double teq = 0.1 * T; const int nb = 10; double zs[10]; int nz = 0; const double tb = (T - teq) / nb;
    double t = 0.0, tnext_blk = teq + tb; int vir_on = 0;
    long nev_switch = full;
    const double t0 = now_s();
    while (t < T) {
        t = fmin(T, t + SIGT);
        edmd3_advance_to(S, t);
        const EDMD3_Health* H = edmd3_health(S);
        if (audits && nev_switch > 0 && H->ev_pair + H->ev_wall + H->ev_cross >= nev_switch) { edmd3_set_schedule_audit(S, every); nev_switch = 0; }
        if (!vir_on && t >= teq) { edmd3_reset_virial(S); vir_on = 1; }
        if (vir_on && t >= tnext_blk - 1e-9 && nz < nb) { zs[nz++] = edmd3_compressibility_Z(S); edmd3_reset_virial(S); tnext_blk += tb; }
        const char* msg; if (edmd3_fatal(S, &msg)) { r->fatal = 1; snprintf(r->msg, sizeof r->msg, "%s", msg); break; }
    }
    r->secs = now_s() - t0;
    if (audits) edmd3_schedule_audit_now(S);
    r->H = *edmd3_health(S); r->A = *edmd3_schedule_audit_stats(S);
    r->cev = edmd3_contact_audit_stats(S, r->cmax);
    r->ev_phys = r->H.ev_pair + r->H.ev_wall; r->ev_cross = r->H.ev_cross;
    r->hash = edmd3_event_hash(S); r->ke1 = edmd3_kinetic_energy(S);
    const EDMD_Particle* P = edmd3_particles(S);
    r->fin = (EDMD_Particle*)malloc((size_t)c->N * sizeof(EDMD_Particle)); memcpy(r->fin, P, (size_t)c->N * sizeof(EDMD_Particle));
    r->H = *edmd3_health(S);
    double m = 0, v = 0; for (int k = 0; k < nz; ++k) m += zs[k]; m /= (nz ? nz : 1);
    for (int k = 0; k < nz; ++k) v += (zs[k] - m) * (zs[k] - m);
    r->Zmean = m; r->Zse = nz > 1 ? sqrt(v / (nz - 1) / nz) : NAN; r->nblk = nz;
    edmd3_destroy(S);
    return 1;
}

typedef struct { double secs; long ev_phys; double ke0, ke1, Zmean, Zse; long forced, clamp, ovl, wod, past; double cmax[4]; long cev; } R2;
static int run2(const Cell* c, double T, int contact, R2* r){
    memset(r, 0, sizeof *r);
    if (contact) setenv("HD_CONTACT_AUDIT", "1", 1); else unsetenv("HD_CONTACT_AUDIT");
    EDMD_Params p = params(c);
    EDMD* S = edmd_create(&p);
    EDMD_Particle* P = (EDMD_Particle*)edmd_particles(S);
    memcpy(P, c->P, (size_t)c->N * sizeof(EDMD_Particle));
    edmd_reschedule_all(S);
    r->ke0 = edmd_total_kinetic_energy(S, 1.0);
    const double teq = 0.1 * T; const int nb = 10; double zs[10]; int nz = 0; const double tb = (T - teq) / nb;
    double t = 0.0, tnext_blk = teq + tb; int vir_on = 0;
    const double t0 = now_s();
    while (t < T) {
        t = fmin(T, t + SIGT);
        edmd_advance_to(S, t);
        if (!vir_on && t >= teq) { edmd_reset_virial(S); vir_on = 1; }
        if (vir_on && t >= tnext_blk - 1e-9 && nz < nb) { zs[nz++] = edmd_compressibility_Z(S); edmd_reset_virial(S); tnext_blk += tb; }
    }
    r->secs = now_s() - t0;
    r->ke1 = edmd_total_kinetic_energy(S, 1.0);
    r->forced = edmd_forced_advance_count(S); r->clamp = edmd_clamp_repair_count(S); r->ovl = edmd_overlap_repair_count(S);
    r->wod = edmd_wall_overdue_count(S); r->past = edmd_past_event_count(S);
    r->cev = edmd_contact_audit_stats(S, r->cmax);
    double m = 0, v = 0; for (int k = 0; k < nz; ++k) m += zs[k]; m /= (nz ? nz : 1);
    for (int k = 0; k < nz; ++k) v += (zs[k] - m) * (zs[k] - m);
    r->Zmean = m; r->Zse = nz > 1 ? sqrt(v / (nz - 1) / nz) : NAN;
    edmd_destroy(S);
    unsetenv("HD_CONTACT_AUDIT");
    return 1;
}
/* gen2 physical events are not counted inside edmd.c: pairs from the virial counter (reset per block), walls from the
   wall counters -- so run2 is repeated without virial resets for the count */
static long count2(const Cell* c, double T, double* secs){
    EDMD_Params p = params(c);
    EDMD* S = edmd_create(&p);
    EDMD_Particle* P = (EDMD_Particle*)edmd_particles(S);
    memcpy(P, c->P, (size_t)c->N * sizeof(EDMD_Particle));
    edmd_reschedule_all(S); edmd_reset_virial(S);
    const double t0 = now_s();
    for (double t = 0.0; t < T;) { t = fmin(T, t + SIGT); edmd_advance_to(S, t); }
    *secs = now_s() - t0;
    long n = edmd_virial_pair_events(S); for (int w = 0; w < 4; ++w) n += edmd_wall_events(S, w);
    edmd_destroy(S);
    return n;
}

static int same_state(const EDMD_Particle* a, const EDMD_Particle* b, int N){ return memcmp(a, b, (size_t)N * sizeof(EDMD_Particle)) == 0; }

static void audit_cell(Cell* c, double T, long full, long every){
    R3 A, B, C; R2 D;
    printf("\n### Cell %s: N = %d, eta = %.4f, box %.4f x %.4f px (%.4f x %.4f sigma), T = %.0f sigma-time, cell width %.0f px\n\n",
           c->name, c->N, c->eta, c->W, c->H, c->W / PX, c->H / PX, T / SIGT, g_cell_px > 0.0 ? g_cell_px : EDMD3_DEFAULT_CELL_PX);
    if (!run3(c, T, 1, full, every, &A) || !run3(c, T, 0, 0, 0, &B) || !run3(c, T, 0, 0, 0, &C)) { printf("RUN FAILED\n"); return; }
    run2(c, T, 1, &D);
    if (A.fatal) printf("FATAL (A): %s\n", A.msg);
    printf("| run | event hash | physical events | crossings | stale | final state equal to A |\n|---|---|---|---|---|---|\n");
    printf("| A gen3, audits on | %016llx | %ld | %ld | %ld | - |\n", (unsigned long long)A.hash, A.ev_phys, A.ev_cross, A.H.ev_stale);
    printf("| B gen3, audits off | %016llx | %ld | %ld | %ld | %s |\n", (unsigned long long)B.hash, B.ev_phys, B.ev_cross, B.H.ev_stale, same_state(A.fin, B.fin, c->N) ? "yes" : "NO");
    printf("| C gen3, audits off, again | %016llx | %ld | %ld | %ld | %s |\n", (unsigned long long)C.hash, C.ev_phys, C.ev_cross, C.H.ev_stale, same_state(A.fin, C.fin, c->N) ? "yes" : "NO");
    printf("\naudits do not steer (A = B): %s; deterministic (B = C): %s\n", (A.hash == B.hash && same_state(A.fin, B.fin, c->N)) ? "YES" : "NO",
           (B.hash == C.hash && same_state(B.fin, C.fin, c->N)) ? "YES" : "NO");
    const EDMD3_Audit* a = &A.A;
    printf("\nschedule audit (gen3, run A): %ld audited states (every event for the first %ld events, then every %ld-th, and at the end)\n\n", a->audits, full, every);
    printf("| class | compared | missing | extra | abs dt > 1e-9 | of them > 1e-10 x horizon | deferred (not yet eligible) |\n|---|---|---|---|---|---|---|\n");
    printf("| pairs | %ld | %ld | %ld | %ld | %ld | %ld |\n", a->pair_cmp, a->pair_missing, a->pair_extra, a->pair_dt, a->pair_dt_rel, a->pair_deferred);
    printf("| outer walls | %ld | %ld | %ld | %ld | %ld | %ld |\n", a->wall_cmp, a->wall_missing, a->wall_extra, a->wall_dt, a->wall_dt_rel, a->wall_deferred);
    printf("| crossings | %ld | %ld | %ld | %ld | %ld | - |\n", a->cross_cmp, a->cross_missing, a->cross_extra, a->cross_dt, a->cross_dt_rel);
    printf("\nsecond live crossing of one disk (must be 0): %ld\n", a->cross_dup);
    printf("\nmax |dt| over matched events %.3g, max |dt|/horizon (|dt| > 1e-9) %.3g; duplicate disagreements %ld; disks outside their cell %ld\n",
           a->max_dt, a->max_rel, a->dup_disagree, a->cell_inconsistent);
    printf("deferred events earlier than the next crossing of one of their disks (must be 0): pairs %ld, walls %ld\n",
           a->pair_deferred_early, a->wall_deferred_early);
    printf("\ncontact audit: gen3 %ld events, max |gap| pairs %.3g px, walls %.3g px; gen2 %ld events, max |gap| pairs %.3g px, walls %.3g px\n",
           A.cev, A.cmax[0], A.cmax[1], D.cev, D.cmax[0], D.cmax[1]);
    const EDMD3_Health* h = &A.H;
    printf("\n[EDMD3-HEALTH] overlap_repair=%ld wall_overdue=%ld past_event=%ld clamp_repair=%ld cell_repair=%ld grid_escape=%ld stagnation=%ld "
           "contact_now=%ld local_checks=%ld local_findings=%ld full_checks=%ld full_findings=%ld local_worst=%.3g full_worst=%.3g "
           "cross_residual_max=%.3g origin_shifts=%ld syncs=%ld heap_compactions=%ld heap_max=%ld\n",
           h->overlap_repair, h->wall_overdue, h->past_event, h->clamp_repair, h->cell_repair, h->grid_escape, h->stagnation,
           h->contact_now, h->local_checks, h->local_findings, h->full_checks, h->full_findings, h->local_worst, h->full_worst,
           h->cross_residual_max, h->origin_shifts, h->syncs, h->heap_compactions, h->heap_max);
    printf("[EDMD-HEALTH gen2] forced_advance=%ld clamp_repair=%ld overlap_repair=%ld wall_overdue=%ld past_event=%ld\n",
           D.forced, D.clamp, D.ovl, D.wod, D.past);
    printf("\n| engine | relative KE drift over the run | Z (virial), mean of %d blocks after T/10 | SE |\n|---|---|---|---|\n", A.nblk);
    printf("| gen3 (B) | %.3g | %.5f | %.5f |\n| gen2 | %.3g | %.5f | %.5f |\n", (B.ke1 - B.ke0) / B.ke0, B.Zmean, B.Zse, (D.ke1 - D.ke0) / D.ke0, D.Zmean, D.Zse);
    {   /* one trajectory per engine: block SEs ignore correlations slower than a block (structure at eta >= 0.70) */
        const double se = sqrt(B.Zse * B.Zse + D.Zse * D.Zse), dz = B.Zmean - D.Zmean;
        printf("\nZ gen3 - gen2 = %.5f, z = %.2f (block SEs of one trajectory each; information, not a test)\n", dz, se > 0.0 ? dz / se : 0.0);
    }
    free(A.fin); free(B.fin); free(C.fin);
}

/* gen2 over T, gen3 over 10 T (more events for its clock), both from the same state, no audits; each timed REP times
   (the event counts are identical every time: deterministic); median and range of the rates */
#define REP 3
static int cmp_d(const void* a, const void* b){ const double x = *(const double*)a, y = *(const double*)b; return (x > y) - (x < y); }
static void speed_row(Cell* c, double T){
    double r2[REP], r3[REP]; long n2 = 0; R3 g3; memset(&g3, 0, sizeof g3);
    for (int k = 0; k < REP; ++k) { double s2; n2 = count2(c, T, &s2); r2[k] = n2 / s2; }
    for (int k = 0; k < REP; ++k) { if (k) free(g3.fin); run3(c, 10.0 * T, 0, 0, 0, &g3); r3[k] = g3.ev_phys / g3.secs; }
    qsort(r2, REP, sizeof(double), cmp_d); qsort(r3, REP, sizeof(double), cmp_d);
    const double ph = (double)(g3.ev_phys ? g3.ev_phys : 1);
    printf("| %d | %.4f | %.0f | %ld | %.3g (%.3g-%.3g) | %.0f | %ld | %.2f | %.2f | %.3g (%.3g-%.3g) | %.1f |\n", c->N, c->eta, T / SIGT, n2,
           r2[REP / 2], r2[0], r2[REP - 1], 10.0 * T / SIGT, g3.ev_phys, (double)g3.ev_cross / ph, (double)g3.H.ev_stale / ph,
           r3[REP / 2], r3[0], r3[REP - 1], r3[REP / 2] / r2[REP / 2]);
    free(g3.fin);
}

static void divergence(Cell* c, double T){
    char err[256]; EDMD_Params p = params(c);
    EDMD3* S3 = edmd3_create(&p, 0.0, err, sizeof err); edmd3_load(S3, c->P, 0.0, err, sizeof err);
    EDMD* S2 = edmd_create(&p); EDMD_Particle* P2 = (EDMD_Particle*)edmd_particles(S2);
    memcpy(P2, c->P, (size_t)c->N * sizeof(EDMD_Particle)); edmd_reschedule_all(S2);
    const double thr[4] = {1e-12, 1e-9, 1e-6, 1.0}; double first[4] = {NAN, NAN, NAN, NAN};
    for (double t = 0.0; t < T;) {
        t += 0.25 * SIGT;
        edmd3_advance_to(S3, t); edmd_advance_to(S2, t);
        const EDMD_Particle *a = edmd3_particles(S3), *b = edmd_particles(S2); double mx = 0;
        for (int i = 0; i < c->N; ++i) { const double d = fmax(fabs(a[i].x - b[i].x), fabs(a[i].y - b[i].y)); if (d > mx) mx = d; }
        for (int k = 0; k < 4; ++k) if (isnan(first[k]) && mx > thr[k]) first[k] = t / SIGT;
    }
    printf("| %s | %d | %.4f | %.4g | %.4g | %.4g | %.4g |\n", c->name, c->N, c->eta, first[0], first[1], first[2], first[3]);
    edmd3_destroy(S3); edmd_destroy(S2);
}

int main(int argc, char** argv){
    const char* mode = argc > 1 ? argv[1] : "";
    const int quick = argc > 2 && strcmp(argv[2], "--quick") == 0;
    if (strcmp(mode, "audit") && strcmp(mode, "speed") && strcmp(mode, "diverge")) { fprintf(stderr, "usage: %s audit|speed|diverge [--quick]\n", argv[0]); return 2; }
    const double Ta = quick ? 40 * SIGT : 400 * SIGT;      /* the audit cells cross one origin shift (341 sigma-time) unless --quick */
    const long full = quick ? 2000 : 10000, every = 500;
    printf("# gen3 M1 harness: %s%s\n\nbuild: %s, double %zu bytes, long double %zu bytes; cell width %.0f px, origin shift %.0f units\n",
           mode, quick ? " (quick)" : "", __VERSION__, sizeof(double), sizeof(long double), EDMD3_DEFAULT_CELL_PX, EDMD3_ORIGIN_SHIFT);
    /* the three M1 cells and the tie stress, N = 400 */
    Cell fluid = { "fluid", 400, PI / 8, 0, 480.0, NULL };
    fluid.W = fluid.N * PI * R * R / (fluid.eta * fluid.H);
    if (!make_rsa(&fluid, 0x6E3A1ULL)) { printf("RSA failed\n"); return 1; }
    Cell dense = { "dense", 400, 0.70, 0, 480.0, NULL };
    dense.W = dense.N * PI * R * R / (dense.eta * dense.H);
    if (!make_hex(&dense, 20, 20, 0.2, 0xD3A5EULL)) { printf("lattice failed\n"); return 1; }
    dense.eta = box_eta(dense.N, dense.W, dense.H);
    Cell lat = { "lattice", 400, 0.70, dense.W, 480.0, NULL };
    if (!make_hex(&lat, 20, 20, 0.0, 0x1A771CEULL)) { printf("lattice failed\n"); return 1; }
    lat.eta = box_eta(lat.N, lat.W, lat.H);
    Cell tie = { "tie", 0, 0, 0, 0, NULL }; make_tie(&tie, 20);
    /* eta 0.85 needs a box commensurate with the lattice: nearest-neighbour distance a with
       (19.5 a + 2R + 1)((19 sqrt(3)/2) a + 2R + 1) = N pi R^2 / 0.85, i.e. the margins of make_hex */
    Cell solid = { "solid", 400, 0.85, 0, 0, NULL };
    {
        const double A2 = 19.5 * 19.0 * sqrt(3.0) / 2.0, B1 = (19.5 + 19.0 * sqrt(3.0) / 2.0) * (2 * R + 1), C0 = (2 * R + 1) * (2 * R + 1) - solid.N * PI * R * R / solid.eta;
        const double a = (-B1 + sqrt(B1 * B1 - 4 * A2 * C0)) / (2 * A2);
        solid.W = 19.5 * a + 2 * R + 1; solid.H = 19.0 * sqrt(3.0) / 2.0 * a + 2 * R + 1;
    }
    if (!make_hex(&solid, 20, 20, 0.2, 0x501DULL)) { printf("lattice failed\n"); return 1; }
    solid.eta = box_eta(solid.N, solid.W, solid.H);
    if (!strcmp(mode, "audit")) {
        printf("\n## Audit cells\n");
        audit_cell(&fluid, Ta, full, every);
        audit_cell(&dense, Ta, full, every);
        audit_cell(&lat, Ta, full, every);
        audit_cell(&tie, quick ? 40 * SIGT : 100 * SIGT, full, every);
        audit_cell(&solid, quick ? 40 * SIGT : 100 * SIGT, full, every);   /* eta >= 0.85, as the gate requires */
        /* the cell width is a runtime parameter (sec. 4.7.1, a): one audited cell at 48 px, and the start-up invariant */
        g_cell_px = 48.0;
        audit_cell(&dense, quick ? 40 * SIGT : 100 * SIGT, full, every);
        g_cell_px = 0.0;
        printf("\n### Start-up invariant of the cell width (w >= d = 24 px, integer px)\n\n| cell width [px] | accepted | message |\n|---|---|---|\n");
        const double ws[4] = {20.0, 23.0, 24.0, 32.5};
        for (int k = 0; k < 4; ++k) {
            char err[256] = ""; EDMD_Params p = params(&fluid);
            EDMD3* S = edmd3_create(&p, ws[k], err, sizeof err);
            printf("| %.1f | %s | %s |\n", ws[k], S ? "yes" : "no", S ? "-" : err);
            edmd3_destroy(S);
        }
    }
    if (!strcmp(mode, "speed")) {
        printf("\n## Events per second on this Mac, same initial state for both engines (gen2 = edmd.c, minimal policy; no dividers; no audits)\n\n");
        printf("events = physical events (pair collisions + outer-wall bounces); gen3 also executes crossings and pops stale events, given per physical event\n\n");
        printf("each rate timed %d times: median (min-max); the event counts are the same every time\n\n", REP);
        printf("| N | eta | gen2 T [sigma-time] | gen2 events | gen2 events/s | gen3 T [sigma-time] | gen3 events | gen3 crossings per event | gen3 stale pops per event | gen3 events/s | gen3 / gen2 (medians) |\n"
               "|---|---|---|---|---|---|---|---|---|---|---|\n");
        const int Ns[2] = {400, 1600}; const double etas[2] = {PI / 8, 0.70};
        for (int a = 0; a < 2; ++a) for (int b = 0; b < 2; ++b) {
            Cell c = { "speed", Ns[a], etas[b], 0, 480.0 * sqrt(Ns[a] / 400.0), NULL };
            c.W = c.N * PI * R * R / (c.eta * c.H);
            if (b == 0) { if (!make_rsa(&c, 0x5EED0ULL + a)) continue; }
            else {
                const int n = (int)lround(sqrt((double)c.N));
                if (!make_hex(&c, n, n, 0.2, 0x5EED1ULL + a)) continue;
                c.eta = box_eta(c.N, c.W, c.H);
            }
            const double T = (quick ? 10.0 : 50.0) * SIGT * (b == 0 ? 4.0 : 1.0) * (a == 0 ? 1.0 : 0.25);
            speed_row(&c, T);
            free(c.P);
        }
    }
    if (!strcmp(mode, "diverge")) {
        printf("\n## Divergence of gen2 and gen3 from one state: first sigma-time at which the largest coordinate difference exceeds\n\n");
        printf("| cell | N | eta | 1e-12 px | 1e-9 px | 1e-6 px | 1 px |\n|---|---|---|---|---|---|---|\n");
        divergence(&fluid, 20 * SIGT);
        divergence(&dense, 20 * SIGT);
        divergence(&lat, 20 * SIGT);
    }
    free(fluid.P); free(dense.P); free(lat.P); free(tie.P); free(solid.P);
    return 0;
}
