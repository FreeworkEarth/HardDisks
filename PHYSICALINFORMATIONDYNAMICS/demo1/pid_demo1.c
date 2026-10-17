/* pid_demo1.c -- PID Demo 1, the headless driver (02_DEMO1_BRIEF.md sec. 5 as amended by sec. 8 and pre-check part B):
   small hard disks inside the frozen lumen of a Demo 0 network, run by the fork of gen3 with its static walls.
     pid_demo1 --state=<demo0 final.txt> --N=n --T=<sigma-time> --seed=S --out=<dir> [--R=12] [--cell=32]
               [--write-every=<sigma-time>] [--audit-every=k] [--kbt=1]
   1. Reads the Demo 0 file with Demo 0's own reader (demo0/pid_state.c); refuses a file without [hull]. Alive disks only, and
      capsules only where child and parent are both alive (B16). Loops are classified by signed area: exactly one outer
      boundary (positive), the rest holes (B18).
   2. Converts to px: x, y, r, cx, cy, R x 24; never the angles (B17). The walls are the [hull] primitives as they are; arcs
      with r <= R stay in the list (the engine never schedules them: B2).
   3. Seeds N centres by random sequential insertion (Demo 0's generator, demo0/pid_rng.c, seeded with --seed): a trial centre
      is kept when it lies at depth >= R + 1e-6 px inside some capsule (the convex hull of a disk and its parent: the lumen
      described by the disks, sec. 3) and overlaps no earlier centre; at most 2000 N trials, then a refusal. Velocities
      Maxwell-Boltzmann at kT (unit mass), total momentum removed, kinetic energy rescaled to N kT exactly.
   4. Runs the fork in stops of --write-every up to T. At t = 0 and at every stop: pstate_<t>.txt (every disk: x, y, vx, vy
      in px and px per unit time, coll_count; %.17g), density_<t>.txt (2 sigma grid: particle count and structural material
      per cell), one line of measure.csv (t, KE, events by type, the run flag, the event hash, the normal impulse delivered to
      each primitive since the previous stop), and the escape test: every centre inside the lumen (Demo 0's hull_contains on
      the [hull] in sigma) and at distance >= R - 1e-9 px from every [hull] primitive (brute force, this file; B15).
   5. Prints "key value" lines and "check ... OK|FAIL" lines (exit 0 iff none fails) and writes manifest.txt (sec. 3).
   Nothing is written outside --out, which must be new or empty. Units as gen3 (sec. 3): 24 px = 1 sigma, 24 units = 1
   sigma-time, unit mass, kT = 1 by default. */
#define _POSIX_C_SOURCE 200809L
#include <errno.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <time.h>
#include "edmd_gen3.h"
#include "pid.h"

#ifndef PID1_BUILD_ID
#define PID1_BUILD_ID "unknown"
#endif
#ifndef PID1_CFLAGS
#define PID1_CFLAGS "unknown"
#endif
#ifndef PID1_SRC_SHA256
#define PID1_SRC_SHA256 "unknown"
#endif
#ifndef PID1_FORK
#define PID1_FORK "unknown"
#endif

#define PXS  24.0          /* px per sigma */
#define TU   24.0          /* internal time units per sigma-time */
#define GRID 2.0           /* density grid cell [sigma] */
#define ESC_TOL 1e-9       /* escape test: centre at distance >= R - ESC_TOL px from every primitive */

/* Demo 0's sources call pid_build_id() (pid_grow.c: world_init); this program never grows a world */
const char* pid_build_id(void){ return "pid_demo1 " PID1_BUILD_ID; }

static int nfail = 0;
static void check(int ok, const char* what){ printf("check %s %s\n", what, ok ? "OK" : "FAIL"); if (!ok) nfail++; }

/* ------------------------------------------------------------------ command line */
typedef struct {
    const char* state; const char* out;
    long N; double T, R, cell, write_every, kbt; uint64_t seed; long audit_every;
    int have_seed;
} Opt;
static int num_d(const char* v, double* x){ char* e; errno = 0; *x = strtod(v, &e); return errno || e == v || *e || !isfinite(*x); }
static int num_l(const char* v, long* x){ char* e; errno = 0; *x = strtol(v, &e, 10); return errno || e == v || *e; }
static int num_u64(const char* v, uint64_t* x){ char* e; errno = 0; *x = strtoull(v, &e, 0); return errno || e == v || *e || v[0] == '-'; }
static void usage(FILE* f){
    fprintf(f, "usage: pid_demo1 --state=<demo0 final.txt> --N=n --T=<sigma-time> --seed=S --out=<dir> [--R=12] [--cell=32]\n"
               "                 [--write-every=<sigma-time>] [--audit-every=k] [--kbt=1]\n"
               "  T and --write-every (default 10) are whole sigma-times; --audit-every=k runs gen3's schedule audit after every\n"
               "  k-th event and its contact audit throughout (0, the default: off; the audits are read-only).\n");
}

/* ------------------------------------------------------------------ geometry in px */
typedef struct { double x1, y1, r1, x2, y2, r2, ux, uy, h; int whole; } Cap;   /* whole: 1 = disk 1 alone, 2 = disk 2 alone */
/* signed distance to the convex hull of two disks (the uneven capsule, in the frame of its axis) */
static double cap_sdf(const Cap* c, double px, double py){
    if (c->whole == 1) return hypot(px - c->x1, py - c->y1) - c->r1;
    if (c->whole == 2) return hypot(px - c->x2, py - c->y2) - c->r2;
    const double dx = px - c->x1, dy = py - c->y1;
    const double along = dx * c->ux + dy * c->uy, lat = fabs(dx * c->uy - dy * c->ux);
    const double b = (c->r1 - c->r2) / c->h, a = sqrt(1.0 - b * b), k = -b * lat + a * along;
    if (k < 0.0) return hypot(lat, along) - c->r1;
    if (k > a * c->h) return hypot(lat, along - c->h) - c->r2;
    return a * lat + b * along - c->r1;
}
/* the unsigned distance from (x, y) to segment k or arc k of the walls [px] */
static double seg_dist(const EDMD3_Seg* s, double x, double y){
    const double dx = s->x2 - s->x1, dy = s->y2 - s->y1, L2 = dx * dx + dy * dy;
    double t = ((x - s->x1) * dx + (y - s->y1) * dy) / L2;
    t = t < 0.0 ? 0.0 : t > 1.0 ? 1.0 : t;
    return hypot(x - (s->x1 + t * dx), y - (s->y1 + t * dy));
}
static double arc_dist(const EDMD3_Arc* a, double x, double y){
    const double qx = x - a->cx, qy = y - a->cy;
    double phi = atan2(qy, qx) - a->a0;
    phi -= 2.0 * PID_PI * floor(phi / (2.0 * PID_PI));
    if (phi <= a->a1 - a->a0) return fabs(hypot(qx, qy) - a->r);
    const double d0 = hypot(x - (a->cx + a->r * cos(a->a0)), y - (a->cy + a->r * sin(a->a0)));
    const double d1 = hypot(x - (a->cx + a->r * cos(a->a1)), y - (a->cy + a->r * sin(a->a1)));
    return d0 < d1 ? d0 : d1;
}

/* ------------------------------------------------------------------ files */
static int mkdir_p(const char* path){
    char buf[1024]; snprintf(buf, sizeof buf, "%s", path);
    for (char* p = buf + 1; *p; ++p) if (*p == '/') { *p = '\0'; if (mkdir(buf, 0755) && errno != EEXIST) return -1; *p = '/'; }
    return mkdir(buf, 0755) && errno != EEXIST ? -1 : 0;
}
static int dir_empty(const char* path){
    DIR* d = opendir(path); if (!d) return 0;
    struct dirent* e; int n = 0;
    while ((e = readdir(d))) if (strcmp(e->d_name, ".") && strcmp(e->d_name, "..")) ++n;
    closedir(d);
    return n == 0;
}
static char (*files)[64] = NULL; static int nfiles = 0;
static void note_file(const char* name){
    for (int i = 0; i < nfiles; ++i) if (!strcmp(files[i], name)) return;
    files = realloc(files, (size_t)(nfiles + 1) * sizeof *files);
    if (!files) { fprintf(stderr, "pid_demo1: out of memory\n"); exit(3); }
    snprintf(files[nfiles++], 64, "%s", name);
}
static FILE* open_out(const Opt* o, const char* name, const char* mode){
    char path[1200]; snprintf(path, sizeof path, "%s/%s", o->out, name);
    FILE* f = fopen(path, mode);
    if (!f) { fprintf(stderr, "pid_demo1: cannot write %s\n", path); exit(2); }
    note_file(name);
    return f;
}
static void close_out(FILE* f, const char* name){
    if (ferror(f) | fclose(f)) { fprintf(stderr, "pid_demo1: write error on %s\n", name); exit(2); }
}

/* ------------------------------------------------------------------ the run */
typedef struct {
    int N; double R; int nseg, narc; EDMD3_Seg* seg; EDMD3_Arc* arc;
    const World* w; const Hull* H;
    int gx, gy;                       /* density grid */
    double* mat;                      /* structural material per density cell [sigma^2] */
    double mat_sum, disk_area, lens; int alive;
    double* Jprev;                    /* per primitive: normal impulse up to the previous stop */
    long esc_total; double min_dist;  /* escape test over every stop; smallest centre-to-primitive distance - R [px] */
} Ctx;

static void write_stop(const Opt* o, Ctx* c, EDMD3* S, double t_sig, const char* ttag){
    const EDMD_Particle* P = edmd3_particles(S);
    const EDMD3_Health* H = edmd3_health(S);
    const uint64_t hash = edmd3_event_hash(S);
    const int clean = edmd3_health_clean(S);
    char name[96];
    /* pstate */
    snprintf(name, sizeof name, "pstate_%s.txt", ttag);
    FILE* f = open_out(o, name, "w");
    fprintf(f, "# pid_demo1 particle state: t %s sigma-time (%.17g internal units); N %d; R %.17g px; event hash %016" PRIx64 "; clean %d\n",
            ttag, t_sig * TU, c->N, c->R, hash, clean);
    fprintf(f, "# units: px (24 px = 1 sigma) and px per internal time unit (24 units = 1 sigma-time); unit mass\n# i x y vx vy coll_count\n");
    for (int i = 0; i < c->N; ++i)
        fprintf(f, "%d %.17g %.17g %.17g %.17g %ld\n", i, P[i].x, P[i].y, P[i].vx, P[i].vy, (long)P[i].coll_count);
    close_out(f, name);
    /* density */
    long* cnt = calloc((size_t)c->gx * c->gy, sizeof(long));
    if (!cnt) { fprintf(stderr, "pid_demo1: out of memory\n"); exit(3); }
    for (int i = 0; i < c->N; ++i) {
        int ix = (int)floor(P[i].x / (GRID * PXS)), iy = (int)floor(P[i].y / (GRID * PXS));
        ix = ix < 0 ? 0 : ix >= c->gx ? c->gx - 1 : ix; iy = iy < 0 ? 0 : iy >= c->gy ? c->gy - 1 : iy;
        cnt[(long)iy * c->gx + ix]++;
    }
    snprintf(name, sizeof name, "density_%s.txt", ttag);
    f = open_out(o, name, "w");
    fprintf(f, "# pid_demo1 density: one grid, two fields (02_DEMO1_BRIEF.md sec. 5); t %s sigma-time\n", ttag);
    fprintf(f, "# grid %d x %d cells of %.0f sigma (%.0f px), cell (ix, iy) = [ix %.0f, (ix + 1) %.0f) x [iy %.0f, (iy + 1) %.0f) sigma\n",
            c->gx, c->gy, GRID, GRID * PXS, GRID, GRID, GRID, GRID);
    fprintf(f, "# particles: the number of particle centres in the cell (sums to N = %d)\n", c->N);
    fprintf(f, "# material: the sum of pi r^2 [sigma^2] of the alive structural disks whose centre lies in the cell: centres binned,\n"
               "#   no overlap subtraction; it sums to %.17g sigma^2 = sum of pi r^2 over the %d alive disks. Demo 0's M_used %.17g\n"
               "#   subtracts the parent-child lens areas (%.17g sigma^2 in all), this field does not (pre-check B14).\n",
            c->disk_area, c->alive, c->w->m_used, c->lens);
    fprintf(f, "# ix iy particles material\n");
    for (int iy = 0; iy < c->gy; ++iy)
        for (int ix = 0; ix < c->gx; ++ix)
            fprintf(f, "%d %d %ld %.17g\n", ix, iy, cnt[(long)iy * c->gx + ix], c->mat[(long)iy * c->gx + ix]);
    close_out(f, name);
    free(cnt);
    /* measure.csv */
    f = open_out(o, "measure.csv", "a");
    const int np = c->nseg + c->narc;
    fprintf(f, "%s,%.17g,%ld,%ld,%ld,%ld,%ld,%ld,%d,%016" PRIx64, ttag, edmd3_kinetic_energy(S), H->ev_pair, H->ev_seg, H->ev_arc, H->ev_wall,
            H->ev_cross, H->ev_stale, clean, hash);
    for (int p = 0; p < np; ++p) {
        const double J = edmd3_prim_impulse(S, p, 0, NULL) + (p < c->nseg ? edmd3_prim_impulse(S, p, 1, NULL) : 0.0);
        fprintf(f, ",%.17g", J - c->Jprev[p]); c->Jprev[p] = J;
    }
    fputc('\n', f);
    close_out(f, "measure.csv");
    /* the escape test (B15): inside the lumen and at distance >= R - tol from every [hull] primitive */
    long esc = 0;
    for (int i = 0; i < c->N; ++i) {
        double d = INFINITY;
        for (int k = 0; k < c->nseg; ++k) { const double e = seg_dist(&c->seg[k], P[i].x, P[i].y); if (e < d) d = e; }
        for (int k = 0; k < c->narc; ++k) { const double e = arc_dist(&c->arc[k], P[i].x, P[i].y); if (e < d) d = e; }
        if (d - c->R < c->min_dist) c->min_dist = d - c->R;
        if (!hull_contains(c->H, P[i].x / PXS, P[i].y / PXS) || !(d >= c->R - ESC_TOL)) esc++;
    }
    c->esc_total += esc;
    printf("stop t %s: events pair %ld seg %ld arc %ld box %ld cross %ld; clean %d; hash %016" PRIx64 "; escaped %ld\n", ttag, H->ev_pair,
           H->ev_seg, H->ev_arc, H->ev_wall, H->ev_cross, clean, hash, esc);
}

int main(int argc, char** argv){
    Opt o; memset(&o, 0, sizeof o);
    o.R = 12.0; o.cell = 32.0; o.write_every = 10.0; o.kbt = 1.0; o.N = -1; o.T = -1.0;
    for (int i = 1; i < argc; ++i) {
        const char* a = argv[i]; const char* eq = strchr(a, '=');
        if (strncmp(a, "--", 2) || !eq) { if (!strcmp(a, "--help")) { usage(stdout); return 0; } fprintf(stderr, "pid_demo1: bad argument '%s'\n", a); usage(stderr); return 2; }
        const size_t kl = (size_t)(eq - a); const char* v = eq + 1; int bad = 0;
        #define KEY(k) (kl == strlen(k) && !strncmp(a, k, kl))
        if      (KEY("--state"))       { o.state = v; bad = !v[0]; }
        else if (KEY("--out"))         { o.out = v; bad = !v[0]; }
        else if (KEY("--N"))           bad = num_l(v, &o.N) || o.N < 1;
        else if (KEY("--T"))           bad = num_d(v, &o.T) || o.T <= 0.0 || o.T != floor(o.T);
        else if (KEY("--seed"))        { bad = num_u64(v, &o.seed); o.have_seed = 1; }
        else if (KEY("--R"))           bad = num_d(v, &o.R) || o.R <= 0.0;
        else if (KEY("--cell"))        bad = num_d(v, &o.cell) || o.cell <= 0.0;
        else if (KEY("--write-every")) bad = num_d(v, &o.write_every) || o.write_every <= 0.0 || o.write_every != floor(o.write_every);
        else if (KEY("--audit-every")) bad = num_l(v, &o.audit_every) || o.audit_every < 0;
        else if (KEY("--kbt"))         bad = num_d(v, &o.kbt) || o.kbt <= 0.0;
        else { fprintf(stderr, "pid_demo1: unknown option '%s'\n", a); usage(stderr); return 2; }
        #undef KEY
        if (bad) { fprintf(stderr, "pid_demo1: bad value in '%s'\n", a); return 2; }
    }
    if (!o.state || !o.out || o.N < 1 || o.T <= 0.0 || !o.have_seed) { fprintf(stderr, "pid_demo1: --state, --N, --T, --seed and --out are required\n"); usage(stderr); return 2; }
    if (mkdir_p(o.out)) { fprintf(stderr, "pid_demo1: cannot create %s\n", o.out); return 2; }
    if (!dir_empty(o.out)) { fprintf(stderr, "pid_demo1: %s is not empty; choose a new --out\n", o.out); return 2; }

    /* 1. the Demo 0 file */
    World w; double mpool_file = 0.0; char err[512];
    if (state_read(&w, o.state, &mpool_file, err, sizeof err)) { fprintf(stderr, "pid_demo1: %s\n", err); return 2; }
    const Hull* H = &w.hull;
    if (!H->valid || H->np == 0) { fprintf(stderr, "pid_demo1: %s has no [hull] (run Demo 0 with --hull=final, the default; state_*.txt files have none)\n", o.state); return 2; }
    int sum = 0; for (int l = 0; l < H->nloop; ++l) sum += H->loop[l];
    if (sum != H->np) { fprintf(stderr, "pid_demo1: [hull] loops hold %d primitives, the file %d\n", sum, H->np); return 2; }
    char fsha[65]; if (sha256_file(o.state, fsha)) { fprintf(stderr, "pid_demo1: cannot read %s\n", o.state); return 2; }
    Ctx c; memset(&c, 0, sizeof c);
    c.N = (int)o.N; c.R = o.R; c.w = &w; c.H = H; c.min_dist = INFINITY;
    /* alive disks, capsules (B16) */
    int ncap = 0; Cap* cap = calloc((size_t)w.nd + 1, sizeof(Cap));
    for (int i = 0; i < w.nd; ++i) {
        const Disk* d = &w.d[i];
        if (!d->alive) continue;
        c.alive++; c.disk_area += PID_PI * d->r * d->r;
        if (d->parent < 0 || !w.d[d->parent].alive) continue;
        const Disk* q = &w.d[d->parent];
        c.lens += pid_lens_area(hypot(d->x - q->x, d->y - q->y), d->r, q->r);   /* zero for tangent pairs (B14) */
        Cap k; memset(&k, 0, sizeof k);
        k.x1 = q->x * PXS; k.y1 = q->y * PXS; k.r1 = q->r * PXS; k.x2 = d->x * PXS; k.y2 = d->y * PXS; k.r2 = d->r * PXS;
        k.h = hypot(k.x2 - k.x1, k.y2 - k.y1);
        if (k.h + k.r2 <= k.r1) k.whole = 1; else if (k.h + k.r1 <= k.r2) k.whole = 2;
        else { k.ux = (k.x2 - k.x1) / k.h; k.uy = (k.y2 - k.y1) / k.h; }
        cap[ncap++] = k;
    }
    /* walls in px; loops by signed area (B18) */
    c.seg = calloc((size_t)H->np, sizeof(EDMD3_Seg)); c.arc = calloc((size_t)H->np, sizeof(EDMD3_Arc));
    int nouter = 0, nhole = 0, k0 = 0;
    for (int l = 0; l < H->nloop; ++l) {
        double area = 0.0;
        for (int j = 0; j < H->loop[l]; ++j) {
            const HPrim* q = &H->p[k0 + j];
            if (q->type == 0) {
                EDMD3_Seg s = {q->v[0] * PXS, q->v[1] * PXS, q->v[2] * PXS, q->v[3] * PXS}; c.seg[c.nseg++] = s;
                area += 0.5 * (q->v[0] * q->v[3] - q->v[2] * q->v[1]);
            } else {
                EDMD3_Arc s = {q->v[0] * PXS, q->v[1] * PXS, q->v[2] * PXS, q->v[3], q->v[4]}; c.arc[c.narc++] = s;
                area += 0.5 * (q->v[2] * (q->v[0] * (sin(q->v[4]) - sin(q->v[3])) - q->v[1] * (cos(q->v[4]) - cos(q->v[3]))) + q->v[2] * q->v[2] * (q->v[4] - q->v[3]));
            }
        }
        if (area > 0.0) nouter++; else nhole++;
        k0 += H->loop[l];
    }
    if (nouter != 1) { fprintf(stderr, "pid_demo1: %d loops of positive area in [hull] (one outer boundary expected)\n", nouter); return 2; }
    int small_arcs = 0; for (int k = 0; k < c.narc; ++k) small_arcs += c.arc[k].r <= o.R;
    printf("demo0 %s: sha256 %s; tick %ld; %d disks, %d alive, %d capsules; [hull] %d loops (1 outer, %d holes), %d segments, %d arcs (%d with r <= R, never hit);"
           " lumen_area %.9f sigma^2; m_used %.9f sigma^2\n", o.state, fsha, w.tick, w.nd, c.alive, ncap, H->nloop, nhole, c.nseg, c.narc, small_arcs, H->area, w.m_used);

    /* 2. the engine state */
    EDMD_Params p; memset(&p, 0, sizeof p);
    p.boxW = w.p.W * PXS; p.boxH = w.p.H * PXS; p.radius = o.R; p.N = c.N; p.pp_collisions_enabled = 1; p.particle_mass = 1.0; p.kB = 1.0;
    EDMD3* S = edmd3_create(&p, o.cell, err, sizeof err);
    if (!S) { fprintf(stderr, "pid_demo1: %s\n", err); return 2; }
    if (!edmd3_set_walls(S, c.seg, c.nseg, c.arc, c.narc, err, sizeof err)) { fprintf(stderr, "pid_demo1: walls: %s\n", err); return 2; }

    /* 3. seeding */
    double bx0 = INFINITY, bx1 = -INFINITY, by0 = INFINITY, by1 = -INFINITY;
    for (int k = 0; k < ncap; ++k) {
        bx0 = fmin(bx0, fmin(cap[k].x1 - cap[k].r1, cap[k].x2 - cap[k].r2)); bx1 = fmax(bx1, fmax(cap[k].x1 + cap[k].r1, cap[k].x2 + cap[k].r2));
        by0 = fmin(by0, fmin(cap[k].y1 - cap[k].r1, cap[k].y2 - cap[k].r2)); by1 = fmax(by1, fmax(cap[k].y1 + cap[k].r1, cap[k].y2 + cap[k].r2));
    }
    Rng g; rng_seed(&g, o.seed);
    EDMD_Particle* P0 = calloc((size_t)c.N, sizeof(EDMD_Particle));
    const long max_trials = 2000L * c.N; long trials = 0;
    const double depth = o.R + 1e-6, d2 = 4.0 * o.R * o.R;
    for (int i = 0; i < c.N; ++i) {
        for (;;) {
            if (++trials > max_trials) {
                fprintf(stderr, "pid_demo1: seeding refused: %d of %d centres placed after %ld trials (the lumen eroded by R is too full)\n", i, c.N, max_trials);
                return 2;
            }
            const double x = bx0 + rng_uniform(&g) * (bx1 - bx0), y = by0 + rng_uniform(&g) * (by1 - by0);
            int in = 0;
            for (int k = 0; k < ncap && !in; ++k) in = cap_sdf(&cap[k], x, y) <= -depth;
            if (!in) continue;
            int ok = 1;
            for (int j = 0; j < i && ok; ++j) { const double dx = x - P0[j].x, dy = y - P0[j].y; if (dx * dx + dy * dy <= d2) ok = 0; }
            if (ok) { P0[i].x = x; P0[i].y = y; break; }
        }
    }
    const double sv = sqrt(o.kbt);
    double mx = 0.0, my = 0.0, ke = 0.0;
    for (int i = 0; i < c.N; ++i) { P0[i].vx = sv * rng_normal(&g); P0[i].vy = sv * rng_normal(&g); mx += P0[i].vx; my += P0[i].vy; }
    mx /= c.N; my /= c.N;
    for (int i = 0; i < c.N; ++i) { P0[i].vx -= mx; P0[i].vy -= my; ke += 0.5 * (P0[i].vx * P0[i].vx + P0[i].vy * P0[i].vy); }
    const double sc = sqrt(c.N * o.kbt / ke);
    for (int i = 0; i < c.N; ++i) { P0[i].vx *= sc; P0[i].vy *= sc; }
    const double eta = c.N * PID_PI * (o.R / PXS) * (o.R / PXS) / H->area;
    printf("seeding: %d centres after %ld trials (Demo 0's generator, seed %" PRIu64 "); eta = N pi (R/sigma)^2 / lumen_area = %.5f\n", c.N, trials, o.seed, eta);
    if (!edmd3_load(S, P0, 0.0, err, sizeof err)) { fprintf(stderr, "pid_demo1: load: %s\n", err); return 2; }
    if (o.audit_every > 0) { edmd3_set_schedule_audit(S, o.audit_every); edmd3_set_contact_audit(S, 1); }

    /* density grid and the structural field */
    c.gx = (int)ceil(w.p.W / GRID); c.gy = (int)ceil(w.p.H / GRID);
    c.mat = calloc((size_t)c.gx * c.gy, sizeof(double));
    for (int i = 0; i < w.nd; ++i) {
        const Disk* d = &w.d[i];
        if (!d->alive) continue;
        int ix = (int)floor(d->x / GRID), iy = (int)floor(d->y / GRID);
        ix = ix < 0 ? 0 : ix >= c.gx ? c.gx - 1 : ix; iy = iy < 0 ? 0 : iy >= c.gy ? c.gy - 1 : iy;
        c.mat[(long)iy * c.gx + ix] += PID_PI * d->r * d->r;
    }
    c.Jprev = calloc((size_t)(c.nseg + c.narc), sizeof(double));
    {
        FILE* f = open_out(&o, "measure.csv", "w");
        fprintf(f, "t_sigma,KE,ev_pair,ev_seg,ev_arc,ev_box,ev_cross,ev_stale,clean,event_hash");
        for (int k = 0; k < c.nseg + c.narc; ++k) fprintf(f, ",J_%s%d", k < c.nseg ? "s" : "a", k < c.nseg ? k : k - c.nseg);
        fputc('\n', f);
        close_out(f, "measure.csv");
    }

    /* 4. the run */
    const double ke0 = edmd3_kinetic_energy(S);
    const long nstop = (long)llround(ceil(o.T / o.write_every));
    char tt[32];
    snprintf(tt, sizeof tt, "%06ld", 0L);
    write_stop(&o, &c, S, 0.0, tt);
    const clock_t c0 = clock();
    int fatal = 0;
    for (long k = 1; k <= nstop; ++k) {
        const double ts = k == nstop ? o.T : k * o.write_every;
        edmd3_advance_to(S, ts * TU);
        const char* msg; if (edmd3_fatal(S, &msg)) { printf("fatal at t %.3f sigma-time: %s\n", edmd3_time(S) / TU, msg); fatal = 1; }
        snprintf(tt, sizeof tt, "%06ld", (long)ts);
        write_stop(&o, &c, S, ts, tt);
        if (fatal) break;
    }
    const double cpu = (double)(clock() - c0) / CLOCKS_PER_SEC;
    const EDMD3_Health* Hh = edmd3_health(S);
    const long nev = Hh->ev_pair + Hh->ev_seg + Hh->ev_arc + Hh->ev_wall + Hh->ev_cross;
    const double drift = (edmd3_kinetic_energy(S) - ke0) / ke0;
    const int clean = edmd3_health_clean(S);
    const uint64_t hash = edmd3_event_hash(S);
    printf("events: %ld (pair %ld, segment %ld, arc %ld, box wall %ld, cell crossing %ld; stale dropped %ld) in %.0f sigma-time: %.1f per sigma-time, "
           "%.2f per sigma-time and particle (collisions only: %.3f); cpu %.2f s (stops and checks included)\n", nev, Hh->ev_pair, Hh->ev_seg, Hh->ev_arc,
           Hh->ev_wall, Hh->ev_cross, Hh->ev_stale, o.T, nev / o.T, nev / o.T / c.N, (Hh->ev_pair * 2.0 + Hh->ev_seg + Hh->ev_arc) / o.T / c.N, cpu);
    printf("health: clean %d; wall overlap repairs %ld, findings %ld, wall repairs %ld, receding %ld, contact-now %ld; pair overlap repairs %ld\n",
           clean, Hh->prim_overlap_repair, Hh->prim_findings, Hh->prim_repair, Hh->prim_receding, Hh->prim_contact_now, Hh->overlap_repair);
    printf("relative KE drift %.3e; escape test (inside the lumen and >= R - %.0e px from every primitive, %ld stops x %d centres): %ld escaped; "
           "smallest centre-to-primitive distance - R %.6e px\n", drift, ESC_TOL, nstop + 1, c.N, c.esc_total, c.min_dist);
    if (o.audit_every > 0) {
        edmd3_schedule_audit_now(S);
        const EDMD3_Audit* A = edmd3_schedule_audit_stats(S);
        double cm[2], cb[2]; const long cev = edmd3_contact_audit_prims(S, cm); edmd3_contact_audit_stats(S, cb);
        EDMD3_Tol tol; edmd3_tolerances(S, &tol);
        printf("schedule audit: %ld audits; walls compared %ld, missing %ld, extra %ld, early %ld, off by > 1e-10 x horizon %ld (grazing, the same contact "
               "to rounding: %ld), boundary %ld; pairs missing %ld extra %ld; box walls missing %ld extra %ld; crossings missing %ld extra %ld\n", A->audits,
               A->prim_cmp, A->prim_missing, A->prim_extra, A->prim_deferred_early, A->prim_dt_rel, A->prim_dt_cond, A->prim_boundary, A->pair_missing,
               A->pair_extra, A->wall_missing, A->wall_extra,
               A->cross_missing, A->cross_extra);
        printf("contact audit: %ld wall events; max |gap| segments %.3e px, arcs %.3e px, pairs %.3e px; tol_face %.3e px\n", cev, cm[0], cm[1], cb[0], tol.tol_face);
        check(A->audits > 0 && !A->prim_missing && !A->prim_extra && !A->prim_deferred_early && !A->prim_dt_rel && !A->pair_missing && !A->pair_extra
              && !A->wall_missing && !A->wall_extra && !A->cross_missing && !A->cross_extra && !A->dup_disagree, "schedule audit: nothing missing, extra, early or off");
        check(cm[0] <= tol.tol_face && cm[1] <= tol.tol_face, "contact-audit wall gaps <= tol_face");
    }
    check(clean == 1 && !fatal, "clean = 1");
    check(fabs(drift) <= 1e-12, "|KE drift| <= 1e-12");
    check(c.esc_total == 0, "no escape at any stop");

    /* 5. the manifest (sec. 3) */
    FILE* f = open_out(&o, "manifest.txt", "w");
    fprintf(f, "# pid_demo1 run manifest (02_DEMO1_BRIEF.md sec. 3)\nargv");
    for (int i = 0; i < argc; ++i) fprintf(f, " %s", argv[i]);
    fprintf(f, "\ndemo0_state %s\ndemo0_state_sha256 %s\npx_per_sigma %.0f\nR_px %.17g\nN %d\nT_sigma %.17g\nwrite_every_sigma %.17g\nseed %" PRIu64 "\n",
            o.state, fsha, PXS, o.R, c.N, o.T, o.write_every, o.seed);
    fprintf(f, "kT %.17g\ncell_px %.17g\naudit_every %ld\nfork %s\nbuild %s %s\nsrc_sha256 %s   # demo1/engine/edmd_gen3.c, .h, edmd.h, demo1/pid_demo1.c, concatenated\n",
            o.kbt, edmd3_cell_px(S), o.audit_every, PID1_FORK, PID1_BUILD_ID, PID1_CFLAGS, PID1_SRC_SHA256);
    fprintf(f, "event_hash %016" PRIx64 "\nclean %d\nke_drift %.17g\nescaped %ld\nchecks_failed %d\n", hash, clean, drift, c.esc_total, nfail);
    fprintf(f, "# sha256 of every file written before this one\n");
    for (int i = 0; i < nfiles; ++i) {
        if (!strcmp(files[i], "manifest.txt")) continue;
        char fp[1200], hex[65]; snprintf(fp, sizeof fp, "%s/%s", o.out, files[i]);
        if (sha256_file(fp, hex)) snprintf(hex, sizeof hex, "unreadable");
        fprintf(f, "%s  %s\n", hex, files[i]);
    }
    close_out(f, "manifest.txt");
    edmd3_destroy(S); world_free(&w);
    free(cap); free(c.seg); free(c.arc); free(c.mat); free(c.Jprev); free(P0); free(files);
    printf("failed checks: %d\n", nfail);
    if (fflush(stdout) != 0 || ferror(stdout)) return 1;
    return nfail ? 1 : 0;
}
