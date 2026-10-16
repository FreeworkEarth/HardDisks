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
   - A disk's kind is "tip" from its birth until it gets its first child, then "body"; the seed stays "seed".

   Demo 0b (brief sec. 3 "0b additions" as amended by the decisions relayed on 2026-10-15 as sec. 8; that section was not on
   disk when this was written, so the readings below are the reference until it is):
   - Patches {id, x, y, r, s, active}, absolute sigma, from --patch=x,y,r,s (repeatable) or a preset. A patch whose disk
     contains the seed centre is the home patch (A): it neither emits cue nor feeds a tip, because the network starts on
     it (otherwise every start disk would be "fed" at tick 0 and the cue would pull every tip back to the seed).
   - Cue c(x) = sum over active non-home patches of s exp(-|x - x_p| / lambda), closed form, and its gradient.
   - Guided rule: theta' = theta + kappa xi + beta wrap(phi - theta), phi the direction of grad c at the tip's disk (no bias
     where the gradient is zero). Same draws in the same order as the unbiased rule, so beta = 0 reproduces it exactly.
   - Fed: a disk that intersects a non-home active patch (|c - x_p| < r + r_p). The fed paths are those disks and all their
     ancestors; a tip whose new disk is fed stops (and makes no branch tip). Fed status is derived from the disks, never
     stored.
   - Adaptation every Tadapt ticks (0 = off, the 0a/0a.1 behaviour): fed disks r <- min(rmax, r (1 + gamma)), charged to
     the pool with the change of the lens areas with parent and children, refused if the pool cannot pay, if the larger
     disk leaves the box, or if it would overlap a non-adjacent disk; then every non-fed body disk r <- r (1 - deltathin),
     material returned; then every disk below rmin / 2 is removed with its subtree and their tips. Parent-child pairs may
     overlap after thickening; their lens is subtracted in M_used (pid_lens_area). Order: ids ascending.
   - Executor (sec. 8): with budget_wait = 1 a tip the pool cannot afford waits (stays alive, nothing placed); a tip blocked
     by geometry still dies after 1 + ntry attempts. With tip_order = random the live tips act in a fresh Fisher-Yates
     permutation each tick drawn from the one generator. Both are off by default (0a/0a.1 behaviour, regression intact).

   Demo 0c (brief sec. 3 "0c", sec. 6 test 6, sec. 9; pid_hull.c):
   - Lumen = union of the capsules of the live edges (child and parent both alive); a capsule is the convex hull of its two
     disks, or the larger disk when one contains the other. The seed alone has no edge: empty lumen, area 0.
   - Hull = closed loops of primitives with the lumen on the left: "seg x1 y1 x2 y2", "arc cx cy R a0 a1" (counter-clockwise,
     a0 in (-pi, pi], a1 > a0). Outer boundaries run counter-clockwise, holes clockwise. Written as [hull] ("loop n" and n
     primitives per loop) and [areas] (lumen_area, hull_perimeter) into final.txt by default (--hull=final), into every
     state file (--hull=all) or nowhere (--hull=off). Area and perimeter are sums over the primitives (Green's theorem), so
     they do not depend on the chaining. Edges are processed in ascending child id, so any input order gives the same bits.
   - Test 6 readings: (c) "error decreasing with resolution" is checked as a falling log-log trend of the raster error over
     2, 4, 8, 16 px per sigma, 16 px better than 2 px, and the finest within 0.1 %, for lumens of at least 1000 sigma^2
     (pixel-centre sampling noise is a few pixels at every resolution, comparable with 0.1 % of a small lumen); a single
     capsule is held to the closed form r1^2 (pi - a) + r2^2 a + (r1 + r2) d sin a, cos a = (r1 - r2)/d, to 1e-12.
     (d) needs shapely and skips when it is missing (the case on this Mac). (e) gives identical primitives and area.
   - View (sec. 9): --alpha and [ ] apply to every fill (disks of every kind and patches); at alpha 0 only outlines remain.
     H draws the hull; --view=network is alpha 0 with the hull and the graph; --view=both keeps the disks and adds the hull.

   Demo 0b.1 (brief sec. 9, last paragraph; preset mat45r). Every new behaviour is a flag that is off by default, so 0b
   (preset mat45) stays the frozen baseline bit for bit:
   - --thin-all=1: thinning applies to every non-fed disk except the seed, dead ends and rim disks included (sec. 9 corrects
     the sec. 3 wording "every other body disk").
   - --regrow-tips=n: when a placed disk is the first to touch an active non-home patch, that disk gets n new tips with
     headings spread evenly over a fan of --regrow-fan radians centred on the arrival heading (the direction from its parent
     to it), i.e. away from where the network came from. They act from the next tick through the normal executor, paid from
     the pool, waiting when it is empty. The first disk to touch a patch is its lowest-id alive touching disk (fed disks
     are never thinned, so it stays); it holds 1 + n tip slots.
   - With regrow on, a disk feeds (and its tip stops) only when its lineage enters the patch for the first time: it touches
     the patch and none of its ancestors does. The regrown front starts inside the food and must grow out of it, and a front
     branch that curls back into the food is still front; strands of the old mat that arrive from outside feed as in 0b.
   - --cue-off-reached=1: a reached patch (some alive disk touches it) keeps feeding but stops emitting cue, so the new
     front is not pulled back into the food it has engulfed. The reached flags are refreshed at the start of every tick
     (and when a state is read), so every tip in a tick sees the same field.
   - The new front = the disks descending from a regrow disk that are not on a fed path. */
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

#define PID_MAX_PATCHES 16
typedef struct { double x, y, r, s; int active; } Patch;

typedef struct {
    double W, H, budget;               /* box and material budget M_total */
    double kappa, tau, rmin, pbranch, delta;
    int    ntry, K;
    double R0;
    int    rings;                      /* 0a.1: hex rings at the start (0: the 0a start with K tips) */
    double pbranch_exp;                /* 0a.1: p_branch(r) = pbranch (r / R0)^pbranch_exp */
    double beta, lambda;               /* 0b: guided bias and cue length */
    int    tadapt;                     /* 0b: adaptation every tadapt ticks, 0 = off */
    double gamma, deltathin, rmax;     /* 0b: thickening, thinning, largest radius */
    int    budget_wait, tip_order;     /* 0b (sec. 8): 1 = budget-blocked tips wait; 0 = id order, 1 = random order */
    int    thin_all, regrow_tips, cue_off_reached;   /* 0b.1: off (0) by default */
    double regrow_fan;                 /* 0b.1: radians */
    int    npatch; Patch patch[PID_MAX_PATCHES];
} Params;

/* 0c: the hull, a list of closed loops of primitives with the lumen on the left (pid_hull.c) */
typedef struct { int type; double v[5]; } HPrim;   /* type 0 seg: x1 y1 x2 y2; type 1 arc: cx cy R a0 a1 (ccw, a1 > a0) */
typedef struct { int valid; HPrim* p; int np, cap; int* loop; int nloop, cap_loop; int open; double area, perim; } Hull;

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
    Hull     hull;                     /* 0c: written as [hull] and [areas] when valid; derived from the disks */
    unsigned char reached[PID_MAX_PATCHES];   /* 0b.1: derived at the start of each tick (world_refresh_reached) */
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
Action policy_guided(const World* w, const Tip* t, Rng* g, const Params* p);
int    pid_patch_home(const World* w, int k);      /* 1 if patch k contains the seed centre */
int    pid_disk_in_patch(const Disk* d, const Patch* q);
void   pid_cue(const World* w, double x, double y, double* c, double* gx, double* gy);
int    world_fed(const World* w, unsigned char* fed);   /* fed[i] = 1 on a fed path; returns how many */
int    world_alive_disks(const World* w);
void   world_refresh_reached(World* w);                 /* 0b.1 */
int    world_regrow_disk(const World* w, int k);        /* 0b.1: the first disk to touch patch k, or -1 */
int    world_front(const World* w, unsigned char* front);   /* 0b.1: front[i] = 1 on the new front; returns how many */
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

/* ---------------------------------------------------------------- hull (pid_hull.c) */
int    hull_edges(const World* w, int* edge_child);            /* live edges as child ids, ascending; returns the count */
int    hull_compute(const World* w, const int* edge_child, int ne, Hull* H);   /* any edge order; processed ascending */
int    hull_compute_world(const World* w, Hull* H);
void   hull_free(Hull* H);
void   hull_add(Hull* H, int type, double v0, double v1, double v2, double v3, double v4);
void   hull_add_loop(Hull* H, int n);
void   hull_prim_ends(const HPrim* q, double* xs, double* ys, double* xe, double* ye);
double hull_max_gap(const Hull* H);                            /* test 6 (a) */
int    hull_contains(const Hull* H, double x, double y);       /* test 6 (b) */
double hull_raster_area(const World* w, double px_per_sigma);  /* test 6 (c) */

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
    int    hull_mode;                  /* 0c: 0 off, 1 in final.txt (default), 2 in every state file */
    double alpha;                      /* 0c: initial disk fill opacity 0..1 (window) */
    int    view_network;               /* 0c: --view: 0 disks (default), 1 network (alpha 0, hull and graph), 2 both */
    long   shot_ticks[16]; int nshot_ticks;   /* 0b.1: --shot-ticks=a,b,...: a BMP at each of these ticks (window) */
} Run;

int  run_tick(Run* r);           /* one step() plus the writes it is due; 0 when --steps is reached */
void run_reset(Run* r);          /* back to tick 0 with the same seed; writes state_000000.txt again */
int  run_finish(Run* r, int completed);           /* manifest.txt */
void run_status(const Run* r, FILE* f);
const char* pid_build_id(void);
int  view_run(Run* r);           /* the window; draws, never owns state */

#endif
