/* PID Demo 1a: the final-state check of the box reference. The brief (02_DEMO1_BRIEF.md sec. 1, "1a") asks for identical
   "hashes, final states, health" from the fork and the unmodified original; the M1/M2 harnesses print hashes and health but
   no final state (pre-check 02_DEMO1_BRIEF_PRECHECK.md A4). This driver prints the final state, and it also runs the three
   engine paths no harness output covers: the collision hook, checkpoint and restart, and the event log.

   The Makefile builds it twice, against the fork (-Iengine) and against the original (-I hspist3/edmd_core), never both in
   one binary (brief sec. 2 rule 8); tests/run_tests.sh compares the two outputs byte for byte and checks each one's verdicts.
   Two cells from fixed seeds, both with the contact audit on, the schedule audit after every 500th event, and the hook:
     box      N = 400 at eta = pi/8 in 960 x 480 px, the four outer walls only, 400 sigma-time (one origin shift, at 341.3);
     divided  an M2-type cell (cf. the free-divider and piston cells of hspist3/edmd_core/tests/gen3_m2_harness.c:21-31): two
              480 x 480 px compartments of 200 disks at eta = pi/8 either side of a divider of thickness 1.2 px, held until
              40 sigma-time and then free with mass 500; pistons of mass 0 at both outer walls, the right one moving in at
              0.05 sigma per sigma-time over [40, 120) sigma-time; 160 sigma-time; also the schedule audit after every band
              event and body change, and the event log.
   Each cell writes a checkpoint at T/2 and goes on to T; a new state, created and loaded the same way, reads the checkpoint
   and runs from T/2 to T; its event hash, health, ledger and disks must equal the uninterrupted run's, bit for bit.
   Output: per cell the event hash, event counts, the run flag, the schedule-audit and contact-audit results, the hook's
   counts and digest, the restart verdict, the event log's size and digest, the bodies, the impulse, virial and tie
   accessors, the kinetic energy, the health, audit, ledger and tolerance records as raw bytes (only as far as gen3 at
   6b29aa5b defines them: a fork that appends fields, as 1b will, prints the appended bytes on a separate "tail" line), and
   every disk's final position, velocity and collision count. Every double in %a (exact bits); the last line is "# end",
   written only when every cell ran to the end.
   Seeding uses the harnesses' generator (xorshift64, gen3_m2_harness.c:55-57), not demo0/pid_rng.c: this is a test of the
   engine against its own reference; brief sec. 2 rule 5 applies to pid_demo1 (1b). Units as gen3: 24 px = 1 sigma,
   24 units = 1 sigma-time, unit disk mass, kT = 1 (hspist3/edmd_core/tests/gen3_m1_harness.c:14). */
#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <math.h>
#include <stdint.h>
#include "edmd_gen3.h"

#define SIGT 24.0          /* internal units per sigma-time */
static const double R  = 12.0;
static const double PI = 3.14159265358979323846;
static const double TH = 1.2;   /* divider thickness [px], as gen3_m2_harness.c:53 */

/* the public records as far as gen3 at 6b29aa5b defines them (edmd_gen3.h, fork lines 47-129); 1b appends at the end only
   (pre-check B8) and keeps EDMD3_NCLS = 5, so these prefixes stay comparable between a 1b fork and the original */
#define HEALTH_ORIG (offsetof(EDMD3_Health, obj_contact_gap_min) + sizeof(double))
#define AUDIT_ORIG  (offsetof(EDMD3_Audit, cls_dt_hz) + 5 * sizeof(double))
#define LEDGER_ORIG (offsetof(EDMD3_Ledger, J_api) + sizeof(double))
#define TOL_ORIG    (offsetof(EDMD3_Tol, tol_cell) + sizeof(double))

static uint64_t rng_state;   /* xorshift64, as the harnesses (gen3_m2_harness.c:55-57) */
static uint64_t xs64(void){ uint64_t x = rng_state; x ^= x << 13; x ^= x >> 7; x ^= x << 17; return rng_state = x; }
static double urand(void){ return ((double)(xs64() >> 11) + 0.5) / 9007199254740992.0; }

/* FNV-1a 64 over bytes */
static uint64_t fnv(uint64_t h, const void* p, size_t n){
    const unsigned char* b = (const unsigned char*)p;
    for (size_t k = 0; k < n; ++k) { h ^= b[k]; h *= 1099511628211ULL; }
    return h;
}
#define FNV0 1469598103934665603ULL

/* random sequential insertion of disks [a, b) into [x0, x1] x [0, H], centres at least R from every face */
static int rsa(EDMD_Particle* P, int a, int b, double x0, double x1, double H){
    for (int i = a; i < b; ++i) {
        for (long tries = 0;; ++tries) {
            if (tries > 10000000L) return 0;
            const double x = x0 + R + urand() * (x1 - x0 - 2 * R), y = R + urand() * (H - 2 * R);
            int ok = 1;
            for (int j = a; j < i && ok; ++j) { const double dx = x - P[j].x, dy = y - P[j].y; if (dx * dx + dy * dy < 4 * R * R) ok = 0; }
            if (ok) { P[i].x = x; P[i].y = y; P[i].coll_count = 0; break; }
        }
    }
    return 1;
}
/* Maxwell-Boltzmann velocities at kT = 1 (unit mass, px per unit = sigma per sigma-time), total momentum removed, kinetic
   energy rescaled to N kT */
static void mb_velocities(EDMD_Particle* P, int N){
    double px = 0.0, py = 0.0, ke = 0.0;
    for (int i = 0; i < N; ++i) {
        const double u1 = urand(), u2 = urand(), g = sqrt(-2.0 * log(u1)), th = 2.0 * PI * u2;
        P[i].vx = g * cos(th); P[i].vy = g * sin(th); px += P[i].vx; py += P[i].vy;
    }
    px /= N; py /= N;
    for (int i = 0; i < N; ++i) { P[i].vx -= px; P[i].vy -= py; ke += 0.5 * (P[i].vx * P[i].vx + P[i].vy * P[i].vy); }
    const double s = sqrt((double)N / ke);
    for (int i = 0; i < N; ++i) { P[i].vx *= s; P[i].vy *= s; }
}

/* the collision hook: counts per kind and a digest of every field of every call, in event order */
typedef struct { long n[4]; uint64_t h; } HookAcc;
static void hook_fn(void* ctx, const EDMD3_Coll* c){
    HookAcc* A = (HookAcc*)ctx;
    if (c->kind >= 0 && c->kind < 4) A->n[c->kind]++;
    const int iv[5] = {c->kind, c->a, c->b, c->comp_a, c->comp_b};
    const double dv[6] = {c->nx, c->ny, c->vax, c->vay, c->vbx, c->vby};
    A->h = fnv(fnv(A->h, iv, sizeof iv), dv, sizeof dv);
}

/* a record as raw bytes: the prefix gen3 at 6b29aa5b defines, then any appended bytes on a "tail" line */
static void put_record(const char* name, const void* p, size_t n, size_t orig){
    const unsigned char* b = (const unsigned char*)p;
    const size_t m = n < orig ? n : orig;
    printf("%s %zu bytes ", name, m);
    for (size_t k = 0; k < m; ++k) printf("%02x", b[k]);
    printf("\n");
    if (n > orig) { printf("tail %s %zu bytes ", name, n - orig); for (size_t k = orig; k < n; ++k) printf("%02x", b[k]); printf("\n"); }
}

static EDMD3* make_state(const char* name, const EDMD_Params* p, const EDMD_Particle* P0, int divided){
    char err[256];
    EDMD3* S = edmd3_create(p, 0.0, err, sizeof err);
    if (!S) { fprintf(stderr, "%s: %s\n", name, err); return NULL; }
    if (!edmd3_load(S, P0, 0.0, err, sizeof err)) { fprintf(stderr, "%s: %s\n", name, err); edmd3_destroy(S); return NULL; }
    edmd3_set_contact_audit(S, 1); edmd3_set_schedule_audit(S, 500);
    if (divided) edmd3_set_schedule_audit_bodies(S, 1);
    return S;
}
/* advance from the current time to t_end, making the divided cell's body changes on the way: at 40 sigma-time the divider is
   released (mass 500) and the right piston starts in at 0.05 px/unit; at 120 sigma-time the piston stops. 0 on any error. */
static int advance_protocol(EDMD3* S, double t_end, int divided, const char* name){
    static const double tch[2] = {40.0 * SIGT, 120.0 * SIGT};
    const char* msg;
    for (int k = 0; divided && k < 2; ++k) {
        if (tch[k] <= edmd3_time(S) || tch[k] > t_end) continue;
        edmd3_advance_to(S, tch[k]);
        if (edmd3_fatal(S, &msg)) { fprintf(stderr, "%s: fatal at %.17g: %s\n", name, edmd3_time(S), msg); return 0; }
        const int ok = k == 0 ? edmd3_set_divider_motion(S, 0, 500.0, 0.0) && edmd3_set_piston_motion(S, 1, 0.0, -0.05)
                              : edmd3_set_piston_motion(S, 1, 0.0, 0.0);
        if (!ok) { fprintf(stderr, "%s: a body change was refused\n", name); return 0; }
    }
    edmd3_advance_to(S, t_end);
    if (edmd3_fatal(S, &msg)) { fprintf(stderr, "%s: fatal at %.17g: %s\n", name, edmd3_time(S), msg); return 0; }
    return 1;
}

static int same_disks(const EDMD_Particle* a, const EDMD_Particle* b, int N){
    for (int i = 0; i < N; ++i)
        if (memcmp(&a[i].x, &b[i].x, sizeof(double)) || memcmp(&a[i].y, &b[i].y, sizeof(double)) || memcmp(&a[i].vx, &b[i].vx, sizeof(double))
            || memcmp(&a[i].vy, &b[i].vy, sizeof(double)) || a[i].coll_count != b[i].coll_count) return 0;
    return 1;
}

/* one cell: the uninterrupted run with a checkpoint at T/2, the restart from it, then the output; 0 on any error */
static int run_cell(const char* name, const EDMD_Params* p, const EDMD_Particle* P0, double T, int divided){
    const int N = p->N;
    EDMD3* S = make_state(name, p, P0, divided);
    if (!S) return 0;
    HookAcc hook; memset(&hook, 0, sizeof hook); hook.h = FNV0;
    edmd3_set_coll_hook(S, hook_fn, &hook);
    FILE* elog = NULL;
    if (divided) { elog = tmpfile(); if (!elog) { fprintf(stderr, "%s: no temporary file\n", name); edmd3_destroy(S); return 0; } edmd3_set_event_log(S, elog, SIGT); }
    FILE* ck = tmpfile();
    int ok = ck != NULL && advance_protocol(S, 0.5 * T, divided, name);
    if (ok && (!edmd3_checkpoint_write(S, ck) || fflush(ck) != 0)) { fprintf(stderr, "%s: checkpoint write failed\n", name); ok = 0; }
    ok = ok && advance_protocol(S, T, divided, name);
    edmd3_set_coll_hook(S, NULL, NULL);
    EDMD3* S2 = NULL;
    if (ok) {   /* the restart: a new state, created and loaded the same way, the checkpoint read into it, on to T */
        char err[256];
        S2 = make_state(name, p, P0, divided);
        rewind(ck);
        if (!S2 || !edmd3_checkpoint_read(S2, ck, err, sizeof err)) { fprintf(stderr, "%s: checkpoint read failed: %s\n", name, S2 ? err : "-"); ok = 0; }
        ok = ok && advance_protocol(S2, T, divided, name);
    }
    if (ck) fclose(ck);
    if (!ok) { if (S2) edmd3_destroy(S2); if (elog) fclose(elog); edmd3_destroy(S); return 0; }
    const uint64_t hash = edmd3_event_hash(S), hash2 = edmd3_event_hash(S2);
    EDMD_Particle* P = (EDMD_Particle*)malloc((size_t)N * sizeof(EDMD_Particle));
    EDMD_Particle* P2 = (EDMD_Particle*)malloc((size_t)N * sizeof(EDMD_Particle));
    if (!P || !P2) { fprintf(stderr, "%s: out of memory\n", name); return 0; }
    memcpy(P, edmd3_particles(S), (size_t)N * sizeof(EDMD_Particle));
    memcpy(P2, edmd3_particles(S2), (size_t)N * sizeof(EDMD_Particle));
    EDMD3_Ledger L, L2;
    edmd3_ledger(S, &L); edmd3_ledger(S2, &L2);
    const int r_health = memcmp(edmd3_health(S), edmd3_health(S2), sizeof(EDMD3_Health)) == 0;
    const int r_ledger = memcmp(&L, &L2, sizeof L) == 0, r_disks = same_disks(P, P2, N);
    edmd3_destroy(S2);
    long elog_bytes = -1; uint64_t elog_h = FNV0;
    if (elog) {
        fflush(elog); rewind(elog); elog_bytes = 0;
        unsigned char buf[4096]; size_t got;
        while ((got = fread(buf, 1, sizeof buf, elog)) > 0) { elog_h = fnv(elog_h, buf, got); elog_bytes += (long)got; }
        fclose(elog);
    }
    edmd3_schedule_audit_now(S);
    const EDMD3_Health* H = edmd3_health(S);
    const EDMD3_Audit* A = edmd3_schedule_audit_stats(S);
    printf("## cell %s: N = %d, box %a x %a px, T = %a units\n", name, N, p->boxW, p->boxH, T);
    printf("event_hash %016llx\n", (unsigned long long)hash);
    printf("events pair %ld wall %ld cross %ld div %ld piston %ld band %ld stale %ld\n",
           H->ev_pair, H->ev_wall, H->ev_cross, H->ev_div, H->ev_piston, H->ev_band, H->ev_stale);
    printf("run_flag %d time %a origin_shifts %ld\n", edmd3_health_clean(S), edmd3_time(S), H->origin_shifts);
    printf("schedule_audit audits %ld missing %ld %ld %ld %ld %ld extra %ld %ld %ld %ld %ld deferred_early %ld %ld %ld %ld "
           "cross_dup %ld dup_disagree %ld cell_inconsistent %ld bands %ld %ld %ld\n",
           A->audits, A->pair_missing, A->wall_missing, A->cross_missing, A->div_missing, A->pis_missing,
           A->pair_extra, A->wall_extra, A->cross_extra, A->div_extra, A->pis_extra,
           A->pair_deferred_early, A->wall_deferred_early, A->div_deferred_early, A->pis_deferred_early,
           A->cross_dup, A->dup_disagree, A->cell_inconsistent, A->band_missing, A->band_extra, A->band_short);
    double cmax[4];
    const long cev = edmd3_contact_audit_stats4(S, cmax);
    printf("contact_audit events %ld max_gap pairs %a walls %a dividers %a pistons %a\n", cev, cmax[0], cmax[1], cmax[2], cmax[3]);
    printf("hook pair %ld wall %ld divider %ld piston %ld digest %016llx\n", hook.n[0], hook.n[1], hook.n[2], hook.n[3], (unsigned long long)hook.h);
    printf("restart_from_checkpoint at %a units: hash %016llx equal hash %d health %d ledger %d disks %d\n", 0.5 * T,
           (unsigned long long)hash2, hash2 == hash, r_health, r_ledger, r_disks);
    if (divided) {
        double x, v, xl, vl, xr, vr;
        edmd3_divider_state(S, 0, &x, &v); edmd3_piston_state(S, 0, &xl, &vl); edmd3_piston_state(S, 1, &xr, &vr);
        printf("event_log bytes %ld digest %016llx\n", elog_bytes, (unsigned long long)elog_h);
        printf("divider x %a vx %a impulse_left %a impulse_right %a work %a\n", x, v, edmd3_divider_impulse(S, 0, 0),
               edmd3_divider_impulse(S, 0, 1), edmd3_work_divider(S, 0));
        printf("pistons left x %a vx %a right x %a vx %a work_right %a\n", xl, vl, xr, vr, edmd3_work_piston(S, 1));
    }
    long n_live = 0, n_now = 0;
    const long n_ties = edmd3_tie_stats(S, &n_live, &n_now);
    printf("accessors wall_impulse %a %a %a %a wall_events %ld %ld %ld %ld piston_impulse %a %a piston_events %ld %ld "
           "virial_sum %a virial_pair_events %ld ties %ld live %ld now %ld\n",
           edmd3_wall_impulse(S, 0), edmd3_wall_impulse(S, 1), edmd3_wall_impulse(S, 2), edmd3_wall_impulse(S, 3),
           edmd3_wall_events(S, 0), edmd3_wall_events(S, 1), edmd3_wall_events(S, 2), edmd3_wall_events(S, 3),
           edmd3_piston_impulse(S, 0), edmd3_piston_impulse(S, 1), edmd3_piston_events(S, 0), edmd3_piston_events(S, 1),
           edmd3_virial_sum(S), edmd3_virial_pair_events(S), n_ties, n_live, n_now);
    printf("kinetic_energy %a\n", edmd3_kinetic_energy(S));
    EDMD3_Tol Tol;
    edmd3_tolerances(S, &Tol);
    put_record("health", H, sizeof *H, HEALTH_ORIG);
    put_record("audit", A, sizeof *A, AUDIT_ORIG);
    put_record("ledger", &L, sizeof L, LEDGER_ORIG);
    put_record("tolerances", &Tol, sizeof Tol, TOL_ORIG);
    for (int i = 0; i < N; ++i) printf("disk %d %a %a %a %a %d\n", i, P[i].x, P[i].y, P[i].vx, P[i].vy, P[i].coll_count);
    free(P); free(P2);
    edmd3_destroy(S);
    return 1;
}

int main(void){
    printf("# PID Demo 1a final-state check (tests/box_final_state.c): two cells, every double in %%a\n");
    int ok = 1;
    {   /* box: 400 disks at eta = pi/8 */
        EDMD_Params p; memset(&p, 0, sizeof p);
        p.boxW = 960.0; p.boxH = 480.0; p.radius = R; p.N = 400; p.pp_collisions_enabled = 1; p.particle_mass = 1.0; p.kB = 1.0;
        EDMD_Particle* P = (EDMD_Particle*)calloc((size_t)p.N, sizeof(EDMD_Particle));
        rng_state = 0x5EED1A0001ULL;
        if (!P || !rsa(P, 0, p.N, 0.0, p.boxW, p.boxH)) { fprintf(stderr, "box: seeding failed\n"); return 1; }
        mb_velocities(P, p.N);
        ok = ok && run_cell("box", &p, P, 400.0 * SIGT, 0);
        free(P);
    }
    if (ok) {   /* divided: 200 + 200 disks either side of a divider; pistons of mass 0 at both outer walls */
        const double L0 = 480.0, H = 480.0;
        EDMD_Params p; memset(&p, 0, sizeof p);
        p.boxW = 2.0 * L0 + TH; p.boxH = H; p.radius = R; p.N = 400; p.pp_collisions_enabled = 1; p.particle_mass = 1.0; p.kB = 1.0;
        p.divider_count = 1; p.divider_x[0] = L0 + 0.5 * TH; p.divider_thickness[0] = TH; p.divider_mass[0] = 0.0;
        p.has_pistonL = 1; p.pistonL_x = 0.0;    p.pistonL_vx = 0.0; p.pistonL_mass = 0.0;
        p.has_pistonR = 1; p.pistonR_x = p.boxW; p.pistonR_vx = 0.0; p.pistonR_mass = 0.0;
        EDMD_Particle* P = (EDMD_Particle*)calloc((size_t)p.N, sizeof(EDMD_Particle));
        rng_state = 0x5EED1A0002ULL;
        if (!P || !rsa(P, 0, 200, 0.0, L0, H) || !rsa(P, 200, 400, L0 + TH, p.boxW, H)) { fprintf(stderr, "divided: seeding failed\n"); return 1; }
        mb_velocities(P, p.N);
        ok = ok && run_cell("divided", &p, P, 160.0 * SIGT, 1);
        free(P);
    }
    if (ok) printf("# end\n");
    if (fflush(stdout) != 0 || ferror(stdout)) return 1;
    return ok ? 0 : 1;
}
