/* pid.h -- Physical Information Dynamics (provisional name), Demo 0: growth and viewer, no particles.
   Brief: ../01_DEMO0_BRIEF.md. This header is the one place where the conventions are written down; every state file
   repeats the geometry line in its header.

   CONVENTIONS
   - Length unit sigma = the diameter of Demo 1's small transport disks. Areas in sigma^2. Time is the growth tick, an integer;
     the seed is born at tick 0 and the disks placed while tick t is processed are born at t.
   - Box [0, W] x [0, H], y up, origin at the lower-left corner. A disk is placed only if r <= x <= W - r and r <= y <= H - r.
   - Headings in radians, counter-clockwise from +x, kept in (-pi, pi].
   - Structural disks touch, they do not overlap. A child is placed tangent to its parent:
       c_child = c_parent + (r_parent + r_child) (cos theta, sin theta).
     Non-adjacent disks (every pair that is not parent-child) never overlap. Touching is legal: the test is
     |c_i - c_j| < (r_i + r_j) - 1e-9 (r_i + r_j) for an overlap (pid_disks_overlap), so hex-packed neighbours whose rounded
     distance falls an ulp short of r_i + r_j are not counted. Parent ids are smaller than child ids; the seed is disk 0 with
     parent -1.
   - Start (0a.1): --start-rings=n puts n hex rings of tangent disks of radius R0 around the seed (ring k: 6k disks on the
     triangular lattice of spacing 2 R0, rotated as a whole by one random angle; ids ring by ring, inside a ring by angle).
     A ring-1 disk's parent is the seed; a ring-k disk's parent is the lowest-id ring-(k-1) disk it touches. Every disk of
     the outer ring gets one tip pointing radially outward. All start disks are born at tick 0. n = 0 is the 0a start:
     K tips on the seed at evenly spaced headings, rotated by the same one random angle.
   - Material: M_used = sum of pi r^2 over live disks minus the lens area of every live parent-child pair (zero at placement);
     M_pool = M_total - M_used. A placement is charged pi r'^2 and refused when M_pool < pi r'^2.
   - Tips: {disk, heading, alive}. A tip moves to the disk it places; a dead tip stays on its last disk for good.
     Slots: every disk holds as many tip slots as tips ever sat on it, and grown children (born > 0) + tips_on(d) equals
     that number: for a disk born after tick 0, 1, or 2 after a branch; for the start disks, the tips they started with
     (K for the seed when n = 0; 1 for an outer-ring disk; 0 otherwise). This is how "dead tips have no further children"
     is checked.
   - Branching (0a.1): a tip on a disk of radius r branches with probability pbranch (r / R0)^pbranch_exp; exponent 0 is
     the 0a rule (pow(x, 0) is exactly 1, so the draws and decisions are the same).
   - Randomness: one generator, xoshiro256** seeded through splitmix64 from --seed (pid_rng.c). No libc rand()/srand().
     Normals by Box-Muller, two uniforms per normal, no cached second value.
   - State files print every double with %.17g, so a file reads back to the same bits (pid_state.c).
   - Built with -ffp-contract=off (Makefile): no fused multiply-add, so a build does not re-round the geometry by inlining.

   Demo 0a reading of the brief (sec. 3) where it leaves a choice, recorded here once:
   - r' = max(r_min, r_parent (1 - tau)) as written; the clamp means the "r' < r_min" stop cannot fire, so a tip stops on the
     budget or after its placement attempts.
   - "redraw theta' up to n_try times": 1 + n_try placement attempts per tip and tick; each attempt calls the policy again.
   - The tips present when a tick starts act in it; a tip created by a branch acts from the next tick.
   - A disk's kind is "tip" from its birth until it gets its first child, then "body"; the seed stays "seed". */
#ifndef PID_H
#define PID_H

#include <stdint.h>
#include <stdio.h>
#include <stddef.h>

#define PID_PI 3.14159265358979323846
#define PID_FORMAT 1
#define PID_GEOMETRY "tangent-disks; lumen = union of per-edge capsules (convex hull of two disks)"

/* ---------------------------------------------------------------- generator (pid_rng.c) */
typedef struct { uint64_t s[4]; } Rng;
void     rng_seed(Rng* g, uint64_t seed);
uint64_t rng_next(Rng* g);
double   rng_uniform(Rng* g);          /* [0, 1), 53 bits */
double   rng_normal(Rng* g);           /* N(0, 1), Box-Muller, consumes exactly two draws */

/* ---------------------------------------------------------------- world */
enum { KIND_SEED = 0, KIND_BODY = 1, KIND_TIP = 2 };

typedef struct { int id; double x, y, r; int parent; int kind; long born; int alive; } Disk;
typedef struct { int disk; double heading; int alive; } Tip;

typedef struct {
    double W, H, budget;               /* box and material budget M_total */
    double kappa, tau, rmin, pbranch, delta;
    int    ntry, K;
    double R0;
    int    rings;                      /* 0a.1: hex rings at the start (0: the 0a start with K tips) */
    double pbranch_exp;                /* 0a.1: p_branch(r) = pbranch (r / R0)^pbranch_exp */
} Params;

typedef struct World {
    uint64_t seed;
    char     rule[16];
    long     tick;                     /* ticks completed */
    Params   p;
    double   m_used;                   /* running ledger: charged at every placement */
    Disk*    d; int nd, cap_d;
    Tip*     t; int nt, cap_t;
    Rng      rng;
    char     build[256];               /* "<folder short hash>[-dirty] <cflags>" of the build that wrote the state */
} World;

typedef enum { ACT_STOP = 0, ACT_EXTEND = 1, ACT_BRANCH = 2 } ActType;
typedef struct { ActType type; double theta, r, theta2; } Action;

/* The plug-in boundary (brief sec. 3): every growth rule is a function of this type and only this type. */
typedef Action (*Policy)(const World*, const Tip*, Rng*, const Params*);

/* ---------------------------------------------------------------- growth (pid_grow.c) */
void   params_default(Params* p);
int    params_check(const Params* p, char* err, size_t n);
void   world_init(World* w, const Params* p, uint64_t seed, const char* rule);
void   world_free(World* w);
int    world_add_disk(World* w, Disk d);           /* returns the new id */
int    world_add_tip(World* w, Tip t);
Policy policy_by_name(const char* name);
Action policy_unbiased(const World* w, const Tip* t, Rng* g, const Params* p);
int    executor_apply(World* w, int tip, Policy pol); /* 1 if a disk was placed */
void   step(World* w, Policy pol);                 /* one growth tick; the only function that changes a World */
double pid_wrap(double a);
int    pid_disks_overlap(const Disk* a, const Disk* b);
int    pid_disk_in_box(const Disk* d, const Params* p);
double pid_lens_area(double dist, double r1, double r2);
double ledger_recompute(const World* w);
double world_pool(const World* w);
int    world_live_tips(const World* w);
int    pid_start_ring(int id);                     /* ring of a start disk by id: 0 for the seed, k for ring k */
/* Invariants of brief sec. 6 items 2-3 on one state. tangency: 1 in 0a. Returns the number of failures; the first is in err. */
int    world_check(const World* w, int tangency, double m_pool_file, char* err, size_t n);

/* ---------------------------------------------------------------- state file (pid_state.c) */
int  state_write(const World* w, const char* path);
int  state_write_fp(const World* w, FILE* f);
int  state_read(World* w, const char* path, double* m_pool_file, char* err, size_t n);
void params_line(const Params* p, char* buf, size_t n);  /* the "params ..." values, %.17g */

/* ---------------------------------------------------------------- sha256 (pid_sha256.c) */
int  sha256_file(const char* path, char hex[65]);

/* ---------------------------------------------------------------- run (pid_main.c) and window (pid_view.c) */
typedef struct {
    World  w;
    Policy policy;
    Params p; uint64_t seed; char rule[16];
    char   out[1024];
    long   steps, write_every;
    int    headless, exit_at_end, shot_at_end, tpf, tick_ms;
    const char* font;
    const char* preset;          /* 0a.2: --preset=NAME, applied before the explicit flags */
    int    argc; char** argv;
    char (*files)[64]; int nfiles, cap_files;    /* state files written, in first-write order */
    int    final_written;
} Run;

int  run_tick(Run* r);           /* one step() plus the writes it is due; 0 when --steps is reached */
void run_reset(Run* r);          /* back to tick 0 with the same seed; writes state_000000.txt again */
int  run_finish(Run* r, int completed);           /* manifest.txt */
void run_status(const Run* r, FILE* f);
const char* pid_build_id(void);
int  view_run(Run* r);           /* the window; draws, never owns state */

#endif
