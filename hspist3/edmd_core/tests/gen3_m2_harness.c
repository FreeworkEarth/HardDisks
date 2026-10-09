/* ##CHRIS 2026-10-09 (261012 sec. 4.7.4 M2 acceptance, sec. 4.7.6; milestone M2): the M2 harness for the generation-3
   engine with dividers and pistons. For each cell it builds one initial state and runs, from that same state and with
   the same protocol (a hold-then-release of the divider, a piston push):
     A  gen3 with both audits: the schedule audit (live heap against brute force, with the divider, piston and band
        classes) after EVERY event for the first --full events, then after every --every-th event, after every band
        expiry and every API change of a body, and at the end; the contact audit at every event;
     B  gen3 without audits            -> event hash, final disks and bodies equal to A: the audits do not steer;
     C  gen3 again without audits      -> equal to B: same-seed bit identity;
     D  gen2 (edmd.c as at 7b08827, the binary's gen2 engine, default minimal rescheduling policy) with its contact audit
        (HD_CONTACT_AUDIT) and, in this harness, an overlap check of every sampled state (all pairs, walls, divider faces);
   and prints the schedule audit (amendment e: per class max |dt| and max |dt|/horizon over ALL matched events), the
   contact gaps of both engines per class, the health counters and the run flag, the tolerances in force with the
   measured contact errors beside them (amendment b), the momentum and energy ledgers with their rounding scales
   (amendment d), and, as information only, Z (pair virial) and the divider's position and period of both engines and
   the static method's face forces. Cells (N = 400, 400 sigma-time unless named):
     cradle_exact, cradle_round (amendment c; N = 300, 100 sigma-time): rows of three-disk cradles, A -> B <- C with B at rest, so the contacts
        A-B and B-C are simultaneous and share disk B; on dyadic positions (exact ties, c = 0 exactly when B is
        re-predicted against C) and on non-dyadic ones (ties to rounding);
     cradle_round_late (amendments b, c): the non-dyadic cradles with gen3 loaded at t = 8100 units, so the contacts fall
        at origin-relative times next to 2^13, where the time quantum is largest (ulp 9.1e-13 below, 1.8e-12 above);
     free_pi8_M50, free_pi8_M500, free_070_M50, free_070_M500: the production geometry of the profile cells (N_s = 200
        per side, H = 40 sigma, compartments L0 = 10 and 5.604167 sigma, thickness 0.05 sigma), divider held for 40
        sigma-time, then free;
     heavy_070_M4e7: the same at M = 4e7;  held_pi8, held_070: never released (the static method);
     driven_pi8: the divider, held 40 sigma-time, then of mass 0 driven at +-0.2 px/unit, reversed every 10 sigma-time (band
        expiries and velocity changes through the API, as the driver's extra_wall_velocity path);
     spring_pi8: Paper 2 geometry C type: one gas (N = 400, pi/8) against a pre-loaded spring divider (M = 50,
        k = 0.02 kT/px^2), an empty compartment behind it; held 40 sigma-time, then free;
     piston_push: Paper 2 geometry B type: two gases of 200 at eta 0.1013 (H = 20 sigma, 77.5 sigma each), divider held
        40 sigma-time then free (M = 100), pistons of infinite mass at both outer walls (as the driver configures them);
        the right one moves in at 0.05 sigma per sigma-time for 155 sigma-time (7.75 sigma), then stops.
   speed mode: events per second of gen2 and gen3 with a divider at N = 400 (pi/8 and 0.70, held and free M = 500).
   Units as in the engine: px, internal time (24 px = 1 sigma, 24 units = 1 sigma-time); unit disk mass, kT = 1.
   build (from hspist3/): cc -std=c11 -O3 -ffp-contract=off -Wall -Wextra -o <out> edmd_core/tests/gen3_m2_harness.c
                          edmd_core/edmd_gen3.c edmd_core/edmd.c -lm
   usage: <out> audit|speed [--quick]   (separate processes: edmd.c reads HD_CONTACT_AUDIT once per process) */
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
#define DTS 6.0            /* sampling step: 0.25 sigma-time */
static const double R = 12.0;
static const double PI = 3.14159265358979323846;
static const double TH = 1.2;   /* divider thickness: 0.05 sigma (production: --wall-thickness=0.05) */

static uint64_t rng_state;
static uint64_t xs64(void){ uint64_t x = rng_state; x ^= x << 13; x ^= x >> 7; x ^= x << 17; return rng_state = x; }
static double urand(void){ return ((double)(xs64() >> 11) + 0.5) / 9007199254740992.0; }
static double now_s(void){ struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + 1e-9 * ts.tv_nsec; }

typedef struct {
    const char* name; const char* what;
    int N; double W, H; EDMD_Particle* P;
    int ndiv; double dx, dM, dk, dxeq;      /* divider centre, its mass once free, spring constant and rest position */
    double t_rel;                            /* release time (internal units); -1: never released (held throughout) */
    int pistons; double t_p0, t_p1, u_p;     /* pistons of mass 0 at both outer walls; right one at -u_p in [t_p0, t_p1) */
    double T;                                /* run length (internal units) */
    int NL;                                  /* disks left of the divider */
    double t0;                               /* gen3 only: the absolute time of the load (gen2 starts at 0); protocol times relative */
    double v_drive, t_switch;                /* driven divider: after t_rel, mass 0 and velocity +-v_drive, sign flipped every t_switch */
} Cell;

static void gauss_velocities(EDMD_Particle* P, int a, int b){
    for (int i = a; i < b; ++i) {
        const double u1 = urand(), u2 = urand(), g = sqrt(-2.0 * log(u1)), th = 2.0 * PI * u2;
        P[i].vx = g * cos(th); P[i].vy = g * sin(th); P[i].coll_count = 0;
    }
}
/* random sequential insertion of disks [a, b) into [x0, x1] x [0, H] (centres R from every face) */
static int rsa(EDMD_Particle* P, int a, int b, double x0, double x1, double H){
    for (int i = a; i < b; ++i) {
        for (long tries = 0;; ++tries) {
            if (tries > 10000000L) return 0;
            const double x = x0 + R + urand() * (x1 - x0 - 2 * R), y = R + urand() * (H - 2 * R);
            int ok = 1;
            for (int j = a; j < i && ok; ++j) { const double dx = x - P[j].x, dy = y - P[j].y; if (dx * dx + dy * dy < 4 * R * R) ok = 0; }
            if (ok) { P[i].x = x; P[i].y = y; break; }
        }
    }
    gauss_velocities(P, a, b);
    return 1;
}
/* triangular lattice of nc x nr disks [a, a + nc nr) in [x0, x1] x [0, H], margin m from the faces, jitter as in M1 */
static int hexfill(EDMD_Particle* P, int a, int nc, int nr, double x0, double x1, double H, double jitter){
    const double m = 0.5;
    const double ax = (x1 - x0 - 2 * R - 2 * m) / (nc - 0.5), ay = (H - 2 * R - 2 * m) / (nr - 1);
    const double dmin = fmin(ax, sqrt(0.25 * ax * ax + ay * ay));
    if (!(dmin > 2 * R)) return 0;
    const double amp = jitter * 0.5 * (dmin - 2 * R);
    int k = a;
    for (int r = 0; r < nr; ++r) for (int q = 0; q < nc; ++q, ++k) {
        double x = x0 + R + m + q * ax + ((r & 1) ? 0.5 * ax : 0.0), y = R + m + r * ay;
        if (amp > 0.0) { x += (2 * urand() - 1) * amp; y += (2 * urand() - 1) * amp; }
        P[k].x = x; P[k].y = y;
    }
    gauss_velocities(P, a, a + nc * nr);
    return 1;
}

/* the production geometry: N_s disks per side, compartments L0 wide (free width), H high, divider of thickness TH */
static int make_divided(Cell* c, double L0px, double Hpx, int dense, uint64_t seed){
    rng_state = seed;
    c->N = 400; c->NL = 200; c->H = Hpx; c->W = 2.0 * L0px + TH; c->ndiv = 1; c->dx = L0px + 0.5 * TH;
    c->P = (EDMD_Particle*)calloc((size_t)c->N, sizeof(EDMD_Particle));
    if (!dense) return rsa(c->P, 0, 200, 0.0, L0px, Hpx) && rsa(c->P, 200, 400, L0px + TH, c->W, Hpx);
    return hexfill(c->P, 0, 5, 40, 0.0, L0px, Hpx, 0.2) && hexfill(c->P, 200, 5, 40, L0px + TH, c->W, Hpx, 0.2);
}
/* cradles: 10 x 10 triples A -> B <- C along x; exact: integer positions, speeds 1, gaps 4..20 px; round: offsets and
   gaps not dyadic, speeds 0.7 */
static void make_cradle(Cell* c, int exact){
    c->N = 300; c->NL = 0; c->ndiv = 0; c->W = 1152.0; c->H = 448.0;
    c->P = (EDMD_Particle*)calloc((size_t)c->N, sizeof(EDMD_Particle));
    for (int r = 0, k = 0; r < 10; ++r) for (int q = 0; q < 10; ++q) {
        const double g = exact ? 4.0 * (1 + (q + r) % 5) : (4.0 * (1 + (q + r) % 5)) / 3.0 + 0.1;
        const double X = exact ? 64.0 + 112.0 * q : 64.3 + 112.1 * q + 0.01 * r, Y = exact ? 32.0 + 40.0 * r : 32.17 + 40.0 * r;
        const double u = exact ? 1.0 : 0.7;
        EDMD_Particle* A = &c->P[k++]; EDMD_Particle* B = &c->P[k++]; EDMD_Particle* C = &c->P[k++];
        A->x = X - 2 * R - g; A->y = Y; A->vx = u;  A->vy = 0.0;
        B->x = X;             B->y = Y; B->vx = 0.0; B->vy = 0.0;
        C->x = X + 2 * R + g; C->y = Y; C->vx = -u; C->vy = 0.0;
        A->coll_count = B->coll_count = C->coll_count = 0;
    }
}

static EDMD_Params params(const Cell* c){
    EDMD_Params p; memset(&p, 0, sizeof p);
    p.boxW = c->W; p.boxH = c->H; p.radius = R; p.N = c->N; p.pp_collisions_enabled = 1;
    p.particle_mass = 1.0; p.kB = 1.0;
    if (c->ndiv) {
        p.divider_count = 1; p.divider_x[0] = c->dx; p.divider_thickness[0] = TH;
        p.divider_mass[0] = c->t_rel > 0.0 || c->t_rel < -0.5 ? 0.0 : c->dM;    /* held at the start unless free from t = 0 */
        p.divider_vx[0] = 0.0; p.divider_k[0] = c->dk; p.divider_xeq[0] = c->dxeq;
    }
    if (c->pistons) {
        p.has_pistonL = 1; p.pistonL_x = 0.0; p.pistonL_vx = 0.0; p.pistonL_mass = 0.0;
        p.has_pistonR = 1; p.pistonR_x = c->W; p.pistonR_vx = 0.0; p.pistonR_mass = 0.0;
    }
    return p;
}
static int cmp_d(const void* a, const void* b){ const double x = *(const double*)a, y = *(const double*)b; return (x > y) - (x < y); }

/* the protocol's change points on the sampling grid */
typedef struct { double t; int kind; double val; } Change;   /* kind 0 release, 1 piston start, 2 piston stop, 3 divider driven at val */
#define MAXCH 256
static int changes(const Cell* c, Change* ch){
    int n = 0;
    if (c->ndiv && c->t_rel > 0.0) {
        if (c->v_drive > 0.0) { double s = 1.0; for (double t = c->t_rel; t < c->T - 1e-9 && n < MAXCH - 4; t += c->t_switch, s = -s) { ch[n].t = t; ch[n].kind = 3; ch[n++].val = s * c->v_drive; } }
        else { ch[n].t = c->t_rel; ch[n].kind = 0; ch[n++].val = 0.0; }
    }
    if (c->pistons && c->u_p > 0.0) { ch[n].t = c->t_p0; ch[n].kind = 1; ch[n++].val = -c->u_p; ch[n].t = c->t_p1; ch[n].kind = 2; ch[n++].val = 0.0; }
    for (int i = 1; i < n; ++i) for (int j = i; j > 0 && ch[j].t < ch[j - 1].t; --j) { const Change x = ch[j]; ch[j] = ch[j - 1]; ch[j - 1] = x; }
    return n;
}

typedef struct {
    uint64_t hash; double secs; EDMD3_Health H; EDMD3_Audit A; double cmax[4]; long cev; int fatal; char msg[256];
    EDMD_Particle* fin; double body[4];
    EDMD3_Ledger L; EDMD3_Tol tol; int clean;
    double Zm, Zse; int nz;
    double* xs; int nxs;
    double JL, JR, tJ; long nL, nR;    /* divider face impulses over [T/10, T] (the static method) */
    double Wp;                         /* work of the right piston */
    long ev_phys;
} R3;
typedef struct { double secs; double Zm, Zse; int nz; long forced, clamp, ovl, wod, past; double cmax[4]; long cev; double* xs; int nxs; double Wp;
                 double worst; long nbad, nsampled; } R2;

static void zblocks(const double* zs, int nz, double* m, double* se){
    double a = 0, v = 0; for (int k = 0; k < nz; ++k) a += zs[k]; a /= (nz ? nz : 1);
    for (int k = 0; k < nz; ++k) v += (zs[k] - a) * (zs[k] - a);
    *m = a; *se = nz > 1 ? sqrt(v / (nz - 1) / nz) : NAN;
}

static double g_cell_px = 0.0;
static int g_check2 = 0;          /* audit mode: check every sampled gen2 state for overlaps */
static int run3(const Cell* c, int audits, long full, long every, R3* r){
    char err[256]; memset(r, 0, sizeof *r);
    EDMD_Params p = params(c);
    EDMD3* S = edmd3_create(&p, g_cell_px, err, sizeof err);
    if (!S) { fprintf(stderr, "%s: %s\n", c->name, err); return 0; }
    if (!edmd3_load(S, c->P, c->t0, err, sizeof err)) { fprintf(stderr, "%s: %s\n", c->name, err); edmd3_destroy(S); return 0; }
    if (audits) { edmd3_set_contact_audit(S, 1); edmd3_set_schedule_audit(S, 1); edmd3_set_schedule_audit_bodies(S, 1); }
    Change ch[MAXCH]; const int nch = changes(c, ch); int ich = 0;
    const double T = c->T, teq = 0.1 * T; const int nb = 10; double zs[10]; int nz = 0; const double tb = (T - teq) / nb;
    double tnext_blk = teq + tb; int vir_on = 0, jon = 0; double J0L = 0, J0R = 0; long n0L = 0, n0R = 0;
    const int nsamp = (int)(T / DTS) + 2; r->xs = (double*)malloc((size_t)nsamp * sizeof(double));
    long nev_switch = full;
    const double t0 = now_s();
    for (int k = 1; k * DTS <= T + 1e-9; ++k) {
        const double t = k * DTS;
        edmd3_advance_to(S, c->t0 + t);
        while (ich < nch && fabs(ch[ich].t - t) < 1e-9) {
            if (ch[ich].kind == 0) edmd3_set_divider_motion(S, 0, c->dM, 0.0);
            else if (ch[ich].kind == 3) edmd3_set_divider_motion(S, 0, 0.0, ch[ich].val);
            else edmd3_set_piston_motion(S, 1, 0.0, ch[ich].val);
            ++ich;
        }
        if (c->ndiv && (c->t_rel < 0.0 || t >= c->t_rel)) { double x; edmd3_divider_state(S, 0, &x, NULL); r->xs[r->nxs++] = x; }
        const EDMD3_Health* H = edmd3_health(S);
        if (audits && nev_switch > 0 && H->ev_pair + H->ev_wall + H->ev_cross + H->ev_div + H->ev_piston >= nev_switch) { edmd3_set_schedule_audit(S, every); nev_switch = 0; }
        if (!vir_on && t >= teq - 1e-9) { edmd3_reset_virial(S); vir_on = 1; }
        if (c->ndiv && !jon && t >= teq - 1e-9) { J0L = edmd3_divider_impulse(S, 0, 0); J0R = edmd3_divider_impulse(S, 0, 1); n0L = edmd3_divider_events(S, 0, 0); n0R = edmd3_divider_events(S, 0, 1); jon = 1; r->tJ = t; }
        if (vir_on && t >= tnext_blk - 1e-9 && nz < nb) { zs[nz++] = edmd3_compressibility_Z(S); edmd3_reset_virial(S); tnext_blk += tb; }
        const char* msg; if (edmd3_fatal(S, &msg)) { r->fatal = 1; snprintf(r->msg, sizeof r->msg, "%s", msg); break; }
    }
    r->secs = now_s() - t0;
    if (audits) edmd3_schedule_audit_now(S);
    r->H = *edmd3_health(S); r->A = *edmd3_schedule_audit_stats(S);
    r->cev = edmd3_contact_audit_stats4(S, r->cmax);
    r->ev_phys = r->H.ev_pair + r->H.ev_wall + r->H.ev_div + r->H.ev_piston;
    r->hash = edmd3_event_hash(S);
    edmd3_ledger(S, &r->L); edmd3_tolerances(S, &r->tol);
    if (c->ndiv) { edmd3_divider_state(S, 0, &r->body[0], &r->body[1]);
                   r->JL = edmd3_divider_impulse(S, 0, 0) - J0L; r->JR = edmd3_divider_impulse(S, 0, 1) - J0R;
                   r->nL = edmd3_divider_events(S, 0, 0) - n0L; r->nR = edmd3_divider_events(S, 0, 1) - n0R; r->tJ = T - r->tJ; }
    if (c->pistons) { edmd3_piston_state(S, 1, &r->body[2], &r->body[3]); r->Wp = edmd3_work_piston(S, 1); }
    const EDMD_Particle* P = edmd3_particles(S);
    r->fin = (EDMD_Particle*)malloc((size_t)c->N * sizeof(EDMD_Particle)); memcpy(r->fin, P, (size_t)c->N * sizeof(EDMD_Particle));
    r->H = *edmd3_health(S);
    r->clean = edmd3_health_clean(S);
    zblocks(zs, nz, &r->Zm, &r->Zse); r->nz = nz;
    edmd3_destroy(S);
    return 1;
}

static void gen2_change(EDMD* S, const Cell* c, const Change* h){
    if (h->kind == 0) { double m = c->dM, v = 0.0; edmd_set_divider_motions(S, 1, &m, &v); }
    else if (h->kind == 3) { double m = 0.0, v = h->val; edmd_set_divider_motions(S, 1, &m, &v); }
    else {
        const EDMD_Params* q = edmd_params(S);
        edmd_config_pistons(S, 1, q->pistonL_x, q->pistonL_vx, q->pistonL_mass, 1, q->pistonR_x, h->val, q->pistonR_mass);
    }
    edmd_reschedule_all(S);      /* as the driver does after a configuration change (00ALLINONE.c:17011) */
}
static int run2(const Cell* c, R2* r){
    memset(r, 0, sizeof *r);
    EDMD_Params p = params(c);
    EDMD* S = edmd_create(&p);
    EDMD_Particle* P = (EDMD_Particle*)edmd_particles(S);
    memcpy(P, c->P, (size_t)c->N * sizeof(EDMD_Particle));
    edmd_reschedule_all(S);
    Change ch[MAXCH]; const int nch = changes(c, ch); int ich = 0;
    const double T = c->T, teq = 0.1 * T; const int nb = 10; double zs[10]; int nz = 0; const double tb = (T - teq) / nb;
    double tnext_blk = teq + tb; int vir_on = 0;
    const int nsamp = (int)(T / DTS) + 2; r->xs = (double*)malloc((size_t)nsamp * sizeof(double));
    const double t0 = now_s();
    for (int k = 1; k * DTS <= T + 1e-9; ++k) {
        const double t = k * DTS;
        edmd_advance_to(S, t);
        while (ich < nch && fabs(ch[ich].t - t) < 1e-9) gen2_change(S, c, &ch[ich++]);
        if (c->ndiv && (c->t_rel < 0.0 || t >= c->t_rel)) r->xs[r->nxs++] = edmd_params(S)->divider_x[0];
        if (g_check2) {   /* gen2 has no validator here: every pair, outer wall and divider face of the sampled state (audit mode only; not timed) */
            const EDMD_Particle* Q = edmd_particles(S); const EDMD_Params* q = edmd_params(S); double worst = 0.0;
            for (int i = 0; i < c->N; ++i) {
                for (int j = i + 1; j < c->N; ++j) { const double dx = Q[j].x - Q[i].x, dy = Q[j].y - Q[i].y; const double g = sqrt(dx * dx + dy * dy) - 2 * R; if (g < worst) worst = g; }
                const double gw[4] = { Q[i].x - R, (q->boxW - R) - Q[i].x, Q[i].y - R, (q->boxH - R) - Q[i].y };
                for (int s = 0; s < 4; ++s) if (gw[s] < worst) worst = gw[s];
                if (c->ndiv) { const double cx = q->divider_x[0], g = (Q[i].x < cx ? cx - Q[i].x : Q[i].x - cx) - (0.5 * TH + R); if (g < worst) worst = g; }
            }
            if (worst < r->worst) r->worst = worst;
            if (worst < -2.4e-5) r->nbad++;          /* the validator's pair tolerance at R = 12 px */
            r->nsampled++;
        }
        if (!vir_on && t >= teq - 1e-9) { edmd_reset_virial(S); vir_on = 1; }
        if (vir_on && t >= tnext_blk - 1e-9 && nz < nb) { zs[nz++] = edmd_compressibility_Z(S); edmd_reset_virial(S); tnext_blk += tb; }
    }
    r->secs = now_s() - t0;
    r->forced = edmd_forced_advance_count(S); r->clamp = edmd_clamp_repair_count(S); r->ovl = edmd_overlap_repair_count(S);
    r->wod = edmd_wall_overdue_count(S); r->past = edmd_past_event_count(S);
    r->cev = edmd_contact_audit_stats(S, r->cmax);
    r->Wp = edmd_work_pistonR(S);
    zblocks(zs, nz, &r->Zm, &r->Zse); r->nz = nz;
    edmd_destroy(S);
    return 1;
}

/* the divider's oscillation (information): mean and SD of the sampled position; the period of the largest peak of the
   power spectrum of the linearly detrended, Hann-windowed record, searched from 3 periods per record up (a slow drift after
   the release would otherwise win at 1 period); jpk = the number of periods in the record (0 if the body did not move) */
static void period(const double* x, int n, double* mean, double* sd, double* p_dft, int* jpk){
    double m = 0, v = 0; for (int k = 0; k < n; ++k) m += x[k]; m /= (n ? n : 1);
    for (int k = 0; k < n; ++k) v += (x[k] - m) * (x[k] - m);
    *mean = m; *sd = n > 1 ? sqrt(v / (n - 1)) : NAN; *jpk = 0; *p_dft = NAN;
    double lo = INFINITY, hi = -INFINITY; for (int k = 0; k < n; ++k) { if (x[k] < lo) lo = x[k]; if (x[k] > hi) hi = x[k]; }
    if (!(n > 8) || !(hi > lo)) return;                 /* identical samples (a held body): no period */
    double sk = 0, skk = 0, sx = 0, skx = 0;
    for (int k = 0; k < n; ++k) { sk += k; skk += (double)k * k; sx += x[k]; skx += k * x[k]; }
    const double b = (n * skx - sk * sx) / (n * skk - sk * sk), a = (sx - b * sk) / n;
    double best = -1; int jb = 0;
    for (int j = 3; j <= n / 2; ++j) {
        double re = 0, im = 0;
        for (int k = 0; k < n; ++k) {
            const double r = (x[k] - a - b * k) * 0.5 * (1.0 - cos(2.0 * PI * k / (n - 1))), ph = 2.0 * PI * j * k / n;
            re += r * cos(ph); im -= r * sin(ph);
        }
        const double pw = re * re + im * im;
        if (pw > best) { best = pw; jb = j; }
    }
    *jpk = jb; *p_dft = jb > 0 ? (double)n * DTS / jb / SIGT : NAN;
}
static int same_state(const EDMD_Particle* a, const EDMD_Particle* b, int N){ return memcmp(a, b, (size_t)N * sizeof(EDMD_Particle)) == 0; }

static void audit_cell(Cell* c, long full, long every){
    R3 A, B, C; R2 D;
    printf("\n### Cell %s: %s\n\nN = %d, box %.4f x %.4f px (%.4f x %.4f sigma), T = %.0f sigma-time, cell width %.0f px", c->name, c->what,
           c->N, c->W, c->H, c->W / PX, c->H / PX, c->T / SIGT, g_cell_px > 0.0 ? g_cell_px : EDMD3_DEFAULT_CELL_PX);
    if (c->ndiv) {
        printf("; divider at %.4f px, thickness %.2f px, ", c->dx, TH);
        if (c->t_rel < 0.0) printf("held (mass 0, velocity 0) throughout");
        else if (c->v_drive > 0.0) printf("held (mass 0, velocity 0) until %g sigma-time, then driven (mass 0) at +-%g px/unit, reversed every %g sigma-time", c->t_rel / SIGT, c->v_drive, c->t_switch / SIGT);
        else printf("held (mass 0, velocity 0) until %g sigma-time, then free with mass %g", c->t_rel / SIGT, c->dM);
        if (c->dk > 0.0) printf(", spring k = %g kT/px^2 about x_eq = %.4f px", c->dk, c->dxeq);
    }
    if (c->pistons) printf("; pistons of mass 0 at x = 0 and x = %.4f px, the right one at -%g px/unit from %g to %g sigma-time", c->W, c->u_p, c->t_p0 / SIGT, c->t_p1 / SIGT);
    if (c->t0 > 0.0) printf("; gen3 loaded at t = %g units (origin-relative times up to and past 2^13 = 8192; gen2 starts at 0)", c->t0);
    printf("\n\n");
    fflush(stdout);
    if (!run3(c, 1, full, every, &A) || !run3(c, 0, 0, 0, &B) || !run3(c, 0, 0, 0, &C)) { printf("RUN FAILED\n"); return; }
    run2(c, &D);
    if (A.fatal) printf("FATAL (A): %s\n", A.msg);
    const int sAB = same_state(A.fin, B.fin, c->N) && !memcmp(A.body, B.body, sizeof A.body);
    const int sBC = same_state(B.fin, C.fin, c->N) && !memcmp(B.body, C.body, sizeof B.body);
    printf("| run | event hash | pair | wall | divider | piston | band | crossings | stale | disks and bodies equal to A |\n|---|---|---|---|---|---|---|---|---|---|\n");
    const R3* rr[3] = {&A, &B, &C}; const char* rn[3] = {"A gen3, audits on", "B gen3, audits off", "C gen3, audits off, again"};
    for (int k = 0; k < 3; ++k)
        printf("| %s | %016llx | %ld | %ld | %ld | %ld | %ld | %ld | %ld | %s |\n", rn[k], (unsigned long long)rr[k]->hash, rr[k]->H.ev_pair, rr[k]->H.ev_wall,
               rr[k]->H.ev_div, rr[k]->H.ev_piston, rr[k]->H.ev_band, rr[k]->H.ev_cross, rr[k]->H.ev_stale,
               k == 0 ? "-" : (same_state(A.fin, rr[k]->fin, c->N) && !memcmp(A.body, rr[k]->body, sizeof A.body) ? "yes" : "NO"));
    printf("\naudits do not steer (A = B): %s; same-seed bit identity (B = C): %s\n", (A.hash == B.hash && sAB) ? "YES" : "NO", (B.hash == C.hash && sBC) ? "YES" : "NO");
    const EDMD3_Audit* a = &A.A;
    printf("\nschedule audit (gen3, run A): %ld audited states (every event for the first %ld events, then every %ld-th, after every band expiry and API body change, and at the end)\n\n", a->audits, full, every);
    printf("| class | compared | missing | extra | abs dt > 1e-9 | of them > 1e-10 x horizon | deferred (not yet eligible) | deferred earlier than eligible | max abs dt, all matched | max abs dt / horizon, all matched (horizon = max(t_bruteforce, t_heap) - now) |\n|---|---|---|---|---|---|---|---|---|---|\n");
    printf("| pairs | %ld | %ld | %ld | %ld | %ld | %ld | %ld | %.3g (at horizon %.3g) | %.3g (at horizon %.3g) |\n", a->pair_cmp, a->pair_missing, a->pair_extra, a->pair_dt, a->pair_dt_rel, a->pair_deferred, a->pair_deferred_early, a->cls_max_dt[EDMD3_CLS_PAIR], a->cls_dt_hz[EDMD3_CLS_PAIR], a->cls_max_rel[EDMD3_CLS_PAIR], a->cls_rel_hz[EDMD3_CLS_PAIR]);
    printf("| outer walls | %ld | %ld | %ld | %ld | %ld | %ld | %ld | %.3g (at horizon %.3g) | %.3g (at horizon %.3g) |\n", a->wall_cmp, a->wall_missing, a->wall_extra, a->wall_dt, a->wall_dt_rel, a->wall_deferred, a->wall_deferred_early, a->cls_max_dt[EDMD3_CLS_WALL], a->cls_dt_hz[EDMD3_CLS_WALL], a->cls_max_rel[EDMD3_CLS_WALL], a->cls_rel_hz[EDMD3_CLS_WALL]);
    printf("| crossings | %ld | %ld | %ld | %ld | %ld | - | - | %.3g (at horizon %.3g) | %.3g (at horizon %.3g) |\n", a->cross_cmp, a->cross_missing, a->cross_extra, a->cross_dt, a->cross_dt_rel, a->cls_max_dt[EDMD3_CLS_CROSS], a->cls_dt_hz[EDMD3_CLS_CROSS], a->cls_max_rel[EDMD3_CLS_CROSS], a->cls_rel_hz[EDMD3_CLS_CROSS]);
    printf("| divider faces | %ld | %ld | %ld | %ld | %ld | %ld | %ld | %.3g (at horizon %.3g) | %.3g (at horizon %.3g) |\n", a->div_cmp, a->div_missing, a->div_extra, a->div_dt, a->div_dt_rel, a->div_deferred, a->div_deferred_early, a->cls_max_dt[EDMD3_CLS_DIV], a->cls_dt_hz[EDMD3_CLS_DIV], a->cls_max_rel[EDMD3_CLS_DIV], a->cls_rel_hz[EDMD3_CLS_DIV]);
    printf("| pistons | %ld | %ld | %ld | %ld | %ld | %ld | %ld | %.3g (at horizon %.3g) | %.3g (at horizon %.3g) |\n", a->pis_cmp, a->pis_missing, a->pis_extra, a->pis_dt, a->pis_dt_rel, a->pis_deferred, a->pis_deferred_early, a->cls_max_dt[EDMD3_CLS_PISTON], a->cls_dt_hz[EDMD3_CLS_PISTON], a->cls_max_rel[EDMD3_CLS_PISTON], a->cls_rel_hz[EDMD3_CLS_PISTON]);
    printf("\nbands: missing %ld, extra %ld, short %ld; second live crossing of one disk %ld; duplicate disagreements %ld; disks outside their cell %ld (all must be 0)\n",
           a->band_missing, a->band_extra, a->band_short, a->cross_dup, a->dup_disagree, a->cell_inconsistent);
    printf("\ncontact audit, max |gap| at executed events [px]:\n\n| engine | events | pairs | outer walls | divider faces | pistons |\n|---|---|---|---|---|---|\n");
    printf("| gen3 (A) | %ld | %.3g | %.3g | %.3g | %.3g |\n| gen2 | %ld | %.3g | %.3g | %.3g | %.3g |\n", A.cev, A.cmax[0], A.cmax[1], A.cmax[2], A.cmax[3], D.cev, D.cmax[0], D.cmax[1], D.cmax[2], D.cmax[3]);
    const EDMD3_Health* h = &A.H;
    printf("\n[EDMD3-HEALTH] overlap_repair=%ld wall_overdue=%ld obj_overlap_repair=%ld past_event=%ld clamp_repair=%ld cell_repair=%ld grid_escape=%ld stagnation=%ld "
           "local_findings=%ld full_findings=%ld body_findings=%ld | contact_now=%ld wall_contact_now=%ld obj_contact_now=%ld contact_c_min=%.3g obj_contact_gap_min=%.3g | "
           "local_checks=%ld full_checks=%ld local_worst=%.3g full_worst=%.3g cross_residual_max=%.3g origin_shifts=%ld syncs=%ld heap_compactions=%ld heap_max=%ld\n",
           h->overlap_repair, h->wall_overdue, h->obj_overlap_repair, h->past_event, h->clamp_repair, h->cell_repair, h->grid_escape, h->stagnation,
           h->local_findings, h->full_findings, h->body_findings, h->contact_now, h->wall_contact_now, h->obj_contact_now, h->contact_c_min, h->obj_contact_gap_min,
           h->local_checks, h->full_checks, h->local_worst, h->full_worst, h->cross_residual_max, h->origin_shifts, h->syncs, h->heap_compactions, h->heap_max);
    printf("run flag (edmd3_health_clean): A %d, B %d, C %d\n", A.clean, B.clean, C.clean);
    printf("[EDMD-HEALTH gen2] forced_advance=%ld clamp_repair=%ld overlap_repair=%ld wall_overdue=%ld past_event=%ld | overlap check of %ld sampled states "
           "(every 0.25 sigma-time, all pairs, walls, divider faces): worst surface gap %.3g px, states with a gap below -2.4e-5 px: %ld\n",
           D.forced, D.clamp, D.ovl, D.wod, D.past, D.nsampled, D.worst, D.nbad);
    {   /* amendment b: the tolerances in force beside the measured contact errors */
        const EDMD3_Tol* t = &A.tol;
        printf("\ntolerances (amendment b; end of run A): u_time = %.4g units, E_bound = %.6g kT, m_min = %g, v_ref = %.6g px/unit, K = %g -> c_tol = %.4g px^2 "
               "(= 2 d x %.4g px), tol_face = %.4g px; band margin %.3g px\n", t->u_time, t->E_bound, t->m_min, t->v_ref, t->K, t->c_tol, t->c_tol / (4.0 * R), t->tol_face, t->band_margin);
        printf("measured: pair contact max |gap| %.3g px -> |c| ~ 2 d |gap| = %.3g px^2 = %.3g c_tol; most negative c of a contact_now %.3g px^2 (%.3g c_tol); "
               "face contacts max |gap| %.3g px (%.3g tol_face) [walls %.3g, divider %.3g, pistons %.3g]\n",
               A.cmax[0], 4.0 * R * A.cmax[0], 4.0 * R * A.cmax[0] / t->c_tol, h->contact_c_min, -h->contact_c_min / t->c_tol,
               fmax(A.cmax[1], fmax(A.cmax[2], A.cmax[3])), fmax(A.cmax[1], fmax(A.cmax[2], A.cmax[3])) / t->tol_face, A.cmax[1], A.cmax[2], A.cmax[3]);
        printf("validator scale: tol_pair = %.3g px (c = %.3g px^2 = %.3g c_tol), tol_wall = %.3g px (%.3g tol_face)\n",
               t->tol_pair, 4.0 * R * t->tol_pair, 4.0 * R * t->tol_pair / t->c_tol, t->tol_wall, t->tol_wall / t->tol_face);
    }
    {   /* amendment d: the ledgers (run B) */
        const EDMD3_Ledger* L = &B.L;
        printf("\nledgers (amendment d; run B, end): momentum of the bodies of finite mass, P - P0, against the impulses from outside, J\n\n"
               "| axis | P - P0 | J | residual | scale (u x sum of rounded terms) | residual / scale |\n|---|---|---|---|---|---|\n");
        for (int ax = 0; ax < 2; ++ax) {
            const double dP = L->P[ax] - L->P0[ax], res = dP - L->J[ax];
            printf("| %s | %.10g | %.10g | %.3g | %.3g | %.3g |\n", ax ? "y" : "x", dP, L->J[ax], res, L->scale_P[ax], fabs(res) / L->scale_P[ax]);
        }
        const double dE = L->E - L->E0, resE = dE - L->W;
        printf("| energy | %.10g (E - E0) | %.10g (W) | %.3g | %.3g | %.3g |\n", dE, L->W, resE, L->scale_E, fabs(resE) / L->scale_E);
        printf("\nimpulses from outside (x unless named): walls L %.6g, R %.6g, B (y) %.6g, T (y) %.6g; divider of mass 0 %.6g; spring anchor %.6g; pistons L %.6g, R %.6g; API changes %.6g\n",
               L->J_wall[0], L->J_wall[1], L->J_wall[2], L->J_wall[3], L->J_div[0], L->J_spring[0], L->J_piston[0], L->J_piston[1], L->J_api);
    }
    {   /* observables, information only */
        printf("\nobservables (information, not a test):\n\n| engine | Z (pair virial), mean of %d blocks after T/10 | SE |", B.nz);
        if (c->ndiv) printf(" divider mean x [px] | SD [px] | period [sigma-time] (periods in the record; peak of the detrended, windowed spectrum, >= 3) |");
        if (c->pistons) printf(" work of the right piston [kT] |");
        printf("\n|---|---|---|%s%s\n", c->ndiv ? "---|---|---|" : "", c->pistons ? "---|" : "");
        double m, sd, pd; int jb;
        printf("| gen3 (B) | %.5f | %.5f |", B.Zm, B.Zse);
        if (c->ndiv) { period(B.xs, B.nxs, &m, &sd, &pd, &jb); if (jb) printf(" %.4f | %.4f | %.4g (%d) |", m, sd, pd, jb); else printf(" %.4f | %.4f | - (does not move) |", m, sd); }
        if (c->pistons) printf(" %.6g |", B.Wp);
        printf("\n| gen2 | %.5f | %.5f |", D.Zm, D.Zse);
        if (c->ndiv) { period(D.xs, D.nxs, &m, &sd, &pd, &jb); if (jb) printf(" %.4f | %.4f | %.4g (%d) |", m, sd, pd, jb); else printf(" %.4f | %.4f | - (does not move) |", m, sd); }
        if (c->pistons) printf(" %.6g |", D.Wp);
        printf("\n");
        if (c->ndiv && c->NL > 0 && c->NL < c->N && !(c->v_drive > 0.0)) {     /* held or free; not driven */
            /* the static method: the mean x force of each gas on the divider over [T/10, T], and Z = F L / (N_s kT) with the
               compartment's own kT from its disks at the end (a held divider keeps the compartments' energies separate) */
            double keL = 0, keR = 0; int nL = 0, nR = 0;
            for (int i = 0; i < c->N; ++i) { const double ke = 0.5 * (B.fin[i].vx * B.fin[i].vx + B.fin[i].vy * B.fin[i].vy);
                                             if (B.fin[i].x < B.body[0]) { keL += ke; nL++; } else { keR += ke; nR++; } }
            const double FL = B.JL / B.tJ, FR = -B.JR / B.tJ, LL = B.body[0] - 0.5 * TH, LR = c->W - (B.body[0] + 0.5 * TH);
            printf("\nstatic method (gen3 run B, %ld + %ld divider events over %.0f sigma-time; information): F_left = %.6g kT/px (Z = F L / (N kT) = %.5f, kT_left = %.5f, N = %d), "
                   "F_right = %.6g kT/px (Z = %.5f, kT_right = %.5f, N = %d)\n", B.nL, B.nR, B.tJ / SIGT, FL, FL * LL / (nL * (keL / nL)), keL / nL, nL,
                   FR, nR ? FR * LR / (nR * (keR / nR)) : NAN, nR ? keR / nR : NAN, nR);
        }
    }
    free(A.fin); free(B.fin); free(C.fin); free(A.xs); free(B.xs); free(C.xs); free(D.xs);
}

/* events per second, gen2 against gen3 with a divider, same states, no audits; gen2's divider and wall events counted
   from its gated event log in a separate (untimed) run, its pair events from the virial counter */
#define REP 3
static long count2(const Cell* c, const char* logpath){
    edmd_set_event_log(logpath, SIGT);
    EDMD_Params p = params(c); EDMD* S = edmd_create(&p);
    EDMD_Particle* P = (EDMD_Particle*)edmd_particles(S); memcpy(P, c->P, (size_t)c->N * sizeof(EDMD_Particle));
    edmd_reschedule_all(S); edmd_reset_virial(S);
    Change ch[MAXCH]; const int nch = changes(c, ch); int ich = 0;
    for (int k = 1; k * DTS <= c->T + 1e-9; ++k) { edmd_advance_to(S, k * DTS); while (ich < nch && fabs(ch[ich].t - k * DTS) < 1e-9) gen2_change(S, c, &ch[ich++]); }
    long n = edmd_virial_pair_events(S);
    edmd_destroy(S); edmd_close_event_log();
    FILE* f = fopen(logpath, "r"); char line[512]; long nl = 0;
    if (f) { while (fgets(line, sizeof line, f)) if (line[0] >= '0' && line[0] <= '9') nl++; fclose(f); }
    return n + nl;
}
static void speed_row(Cell* c, const char* logpath){
    double r2[REP], r3[REP]; R3 g3; R2 g2;
    const long n2 = count2(c, logpath);
    for (int k = 0; k < REP; ++k) { run2(c, &g2); r2[k] = n2 / g2.secs; free(g2.xs); }
    for (int k = 0; k < REP; ++k) { run3(c, 0, 0, 0, &g3); r3[k] = g3.ev_phys / g3.secs; if (k < REP - 1) { free(g3.fin); free(g3.xs); } }
    qsort(r2, REP, sizeof(double), cmp_d); qsort(r3, REP, sizeof(double), cmp_d);
    printf("| %s | %.0f | %ld | %.3g (%.3g-%.3g) | %ld | %ld | %.2f | %.2f | %.3g (%.3g-%.3g) | %.1f |\n", c->name, c->T / SIGT, n2, r2[REP / 2], r2[0], r2[REP - 1],
           g3.ev_phys, g3.H.ev_div, (double)g3.H.ev_cross / g3.ev_phys, (double)g3.H.ev_stale / g3.ev_phys, r3[REP / 2], r3[0], r3[REP - 1], r3[REP / 2] / r2[REP / 2]);
    free(g3.fin); free(g3.xs);
}

static Cell divided(const char* name, const char* what, double L0px, int dense, uint64_t seed, double M, double t_rel, double T){
    Cell c; memset(&c, 0, sizeof c); c.name = name; c.what = what;
    if (!make_divided(&c, L0px, 40.0 * PX, dense, seed)) { fprintf(stderr, "%s: seeding failed\n", name); exit(1); }
    c.dM = M; c.t_rel = t_rel; c.T = T;
    return c;
}

int main(int argc, char** argv){
    const char* mode = argc > 1 ? argv[1] : "";
    const int quick = argc > 2 && strcmp(argv[2], "--quick") == 0;
    if (strcmp(mode, "audit") && strcmp(mode, "speed")) { fprintf(stderr, "usage: %s audit|speed [--quick]\n", argv[0]); return 2; }
    const double T = quick ? 40 * SIGT : 400 * SIGT, trel = quick ? 4 * SIGT : 40 * SIGT;
    const long full = quick ? 2000 : 10000, every = 500;
    printf("# gen3 M2 harness: %s%s\n\nbuild: %s, double %zu bytes; cell width %.0f px, origin shift %.0f units, tolerance factor K = %g\n",
           mode, quick ? " (quick)" : "", __VERSION__, sizeof(double), EDMD3_DEFAULT_CELL_PX, EDMD3_ORIGIN_SHIFT, EDMD3_TOL_K);
    const double L8 = 10.0 * PX, L70 = 5.604167 * PX;     /* the profile cells' compartments (cluster/profile_edmd_koa.sh) */
    if (!strcmp(mode, "audit")) {
        setenv("HD_CONTACT_AUDIT", "1", 1);               /* gen2's contact audit, read once per process */
        g_check2 = 1;
        printf("\n## Audit cells\n");
        Cell ce; memset(&ce, 0, sizeof ce); ce.name = "cradle_exact"; ce.what = "100 three-disk cradles on dyadic positions, speeds 1 (amendment c)";
        make_cradle(&ce, 1); ce.T = quick ? 20 * SIGT : 100 * SIGT; audit_cell(&ce, full, every); free(ce.P);
        Cell cr; memset(&cr, 0, sizeof cr); cr.name = "cradle_round"; cr.what = "100 three-disk cradles on non-dyadic positions, speeds 0.7 (amendment c)";
        make_cradle(&cr, 0); cr.T = quick ? 20 * SIGT : 100 * SIGT; audit_cell(&cr, full, every); free(cr.P);
        Cell cl; memset(&cl, 0, sizeof cl); cl.name = "cradle_round_late"; cl.what = "the same, gen3 loaded at t = 8100 units: the contacts at the full time quantum near 2^13 (amendments b, c)";
        make_cradle(&cl, 0); cl.T = quick ? 20 * SIGT : 100 * SIGT; cl.t0 = 8100.0; audit_cell(&cl, full, every); free(cl.P);
        struct { const char* n; const char* w; double L; int dense; uint64_t seed; double M, trel; } dv[7] = {
            {"free_pi8_M50",   "free divider, pi/8, M = 50",            L8,  0, 0xA11CE1ULL, 50.0,  trel},
            {"free_pi8_M500",  "free divider, pi/8, M = 500",           L8,  0, 0xA11CE2ULL, 500.0, trel},
            {"free_070_M50",   "free divider, eta 0.70, M = 50",        L70, 1, 0xA11CE3ULL, 50.0,  trel},
            {"free_070_M500",  "free divider, eta 0.70, M = 500",       L70, 1, 0xA11CE4ULL, 500.0, trel},
            {"heavy_070_M4e7", "heavy free divider, eta 0.70, M = 4e7", L70, 1, 0xA11CE5ULL, 4e7,   trel},
            {"held_pi8",       "held divider (static method), pi/8",    L8,  0, 0xA11CE6ULL, 0.0,  -1.0},
            {"held_070",       "held divider (static method), eta 0.70", L70, 1, 0xA11CE7ULL, 0.0, -1.0}};
        for (int k = 0; k < 7; ++k) { Cell c = divided(dv[k].n, dv[k].w, dv[k].L, dv[k].dense, dv[k].seed, dv[k].M, dv[k].trel, T); audit_cell(&c, full, every); free(c.P); }
        {   /* band expiries and API changes: a divider of mass 0 driven back and forth (the driver's extra_wall_velocity path) */
            Cell c = divided("driven_pi8", "divider driven back and forth, pi/8 (band expiries and API velocity changes)", L8, 0, 0xA11CEAULL, 0.0, trel, T);
            c.v_drive = 0.2; c.t_switch = 10.0 * SIGT;
            audit_cell(&c, full, every); free(c.P);
        }
        {   /* Paper 2 geometry C type: one gas against a pre-loaded spring divider, an empty compartment behind it */
            Cell c; memset(&c, 0, sizeof c); c.name = "spring_pi8"; c.what = "one gas (N = 400, pi/8) against a pre-loaded spring divider (Paper 2 geometry C type)";
            rng_state = 0xA11CE8ULL;
            const double Lg = 2.0 * L8; c.N = 400; c.NL = 400; c.H = 40.0 * PX; c.ndiv = 1; c.dx = Lg + 0.5 * TH; c.W = Lg + TH + L8;
            c.P = (EDMD_Particle*)calloc(400, sizeof(EDMD_Particle));
            if (!rsa(c.P, 0, 400, 0.0, Lg, c.H)) { fprintf(stderr, "spring: seeding failed\n"); return 1; }
            /* pre-load: k (x - x_eq) = the gas force N kT Z / L, Z = 2.76 (KR at pi/8): x_eq = x - N Z / (k L) */
            c.dk = 0.02; c.dxeq = c.dx - 400.0 * 2.76 / (c.dk * Lg); c.dM = 50.0; c.t_rel = trel; c.T = T;
            audit_cell(&c, full, every); free(c.P);
        }
        {   /* Paper 2 geometry B type: the piston push */
            Cell c; memset(&c, 0, sizeof c); c.name = "piston_push"; c.what = "two gases of 200 at eta 0.1013, free divider M = 100, right piston pushes 7.75 sigma (Paper 2 geometry B type)";
            const double Lp = 77.5 * PX;
            rng_state = 0xA11CE9ULL;
            c.N = 400; c.NL = 200; c.H = 20.0 * PX; c.W = 2.0 * Lp + TH; c.ndiv = 1; c.dx = Lp + 0.5 * TH;
            c.P = (EDMD_Particle*)calloc(400, sizeof(EDMD_Particle));
            if (!rsa(c.P, 0, 200, 0.0, Lp, c.H) || !rsa(c.P, 200, 400, Lp + TH, c.W, c.H)) { fprintf(stderr, "piston: seeding failed\n"); return 1; }
            c.dM = 100.0; c.t_rel = trel; c.T = T;
            c.pistons = 1; c.u_p = 0.05; c.t_p0 = trel; c.t_p1 = trel + (quick ? 15.5 : 155.0) * SIGT;
            audit_cell(&c, full, every); free(c.P);
        }
    }
    if (!strcmp(mode, "speed")) {
        const char* logpath = "gen3_m2_speed_evlog.csv";   /* gen2's event log for the counts, in the working directory (overwritten per cell, kept) */
        printf("\n## Events per second on this Mac with a divider at N = 400, same initial state and protocol for both engines (gen2 = edmd.c, minimal policy; no audits)\n\n");
        printf("events = physical events (pair + outer wall + divider + piston); gen2's divider and wall events counted from its event log in a separate untimed run; each rate timed %d times: median (min-max)\n\n", REP);
        printf("| cell | T [sigma-time] | gen2 events | gen2 events/s | gen3 events | of them divider | gen3 crossings per event | gen3 stale pops per event | gen3 events/s | gen3 / gen2 (medians) |\n|---|---|---|---|---|---|---|---|---|---|\n");
        const double Ts = (quick ? 10.0 : 100.0) * SIGT;
        Cell c1 = divided("held_pi8", "", L8, 0, 0x5EED8ULL, 0.0, -1.0, 2.0 * Ts);  speed_row(&c1, logpath); free(c1.P);
        Cell c2 = divided("free_pi8_M500", "", L8, 0, 0x5EED9ULL, 500.0, DTS * round(0.25 * Ts / DTS), 2.0 * Ts); speed_row(&c2, logpath); free(c2.P);
        Cell c3 = divided("held_070", "", L70, 1, 0x5EEDAULL, 0.0, -1.0, 0.5 * Ts); speed_row(&c3, logpath); free(c3.P);
        Cell c4 = divided("free_070_M500", "", L70, 1, 0x5EEDBULL, 500.0, DTS * round(0.0625 * Ts / DTS), 0.5 * Ts); speed_row(&c4, logpath); free(c4.P);
    }
    return 0;
}
