/* pid_grow.c -- the growth model of brief sec. 3 for Demo 0a: seed, tips, tangent placement, branching, the material
   ledger and the invariant checks. step() is the only function that changes a World after world_init(); the window and
   the headless loop both call it (pid_main.c run_tick). */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "pid.h"

void params_default(Params* p){
    /* The brief's defaults (sec. 3) except R0, tuned for a readable picture (Demo 0a, 2026-10-15). Every value is a
       --name=value flag and is written into every state file and the manifest.
       Tuning, 400 ticks, headless: with R0 = 6 the budget of 3000 sigma^2 is spent by tick 7-11 on 31-34 large disks
       (seeds 9300-9302), so no network forms. With R0 = 3 (tau 0.03 as briefed) seeds 9300-9305 grow 251-577 disks and
       stop placing by tick 58-161; tau 0.02 was tried and gives fewer, thicker arms. */
    p->W = 200.0; p->H = 120.0; p->budget = 3000.0;
    p->kappa = 0.35; p->tau = 0.03; p->rmin = 1.0; p->pbranch = 0.08; p->delta = 0.6;
    p->ntry = 8; p->K = 6; p->R0 = 3.0;
}

int params_check(const Params* p, char* err, size_t n){
    if (!(p->W > 0 && p->H > 0))                 { snprintf(err, n, "--W and --H must be > 0"); return -1; }
    if (!(p->budget > 0))                        { snprintf(err, n, "--budget must be > 0"); return -1; }
    if (!(p->tau >= 0 && p->tau < 1))            { snprintf(err, n, "--tau must be in [0, 1)"); return -1; }
    if (!(p->rmin > 0))                          { snprintf(err, n, "--rmin must be > 0"); return -1; }
    if (!(p->pbranch >= 0 && p->pbranch <= 1))   { snprintf(err, n, "--pbranch must be in [0, 1]"); return -1; }
    if (!(p->kappa >= 0) || !(p->delta >= 0))    { snprintf(err, n, "--kappa and --delta must be >= 0"); return -1; }
    if (p->ntry < 0 || p->K < 1)                 { snprintf(err, n, "--ntry must be >= 0 and --K >= 1"); return -1; }
    if (!(p->R0 > 0 && 2 * p->R0 <= p->W && 2 * p->R0 <= p->H)) { snprintf(err, n, "--R0 must be > 0 and fit in the box"); return -1; }
    if (!(PID_PI * p->R0 * p->R0 <= p->budget))  { snprintf(err, n, "--budget must hold the seed disk (pi R0^2)"); return -1; }
    return 0;
}

double pid_wrap(double a){
    while (a > PID_PI) a -= 2.0 * PID_PI;
    while (a <= -PID_PI) a += 2.0 * PID_PI;
    return a;
}

int world_add_disk(World* w, Disk d){
    if (w->nd == w->cap_d) {
        w->cap_d = w->cap_d ? 2 * w->cap_d : 256;
        Disk* nd = (Disk*)realloc(w->d, (size_t)w->cap_d * sizeof(Disk));
        if (!nd) { fprintf(stderr, "pid: out of memory (disks)\n"); exit(3); }
        w->d = nd;
    }
    d.id = w->nd; w->d[w->nd] = d;
    return w->nd++;
}
int world_add_tip(World* w, Tip t){
    if (w->nt == w->cap_t) {
        w->cap_t = w->cap_t ? 2 * w->cap_t : 64;
        Tip* nt = (Tip*)realloc(w->t, (size_t)w->cap_t * sizeof(Tip));
        if (!nt) { fprintf(stderr, "pid: out of memory (tips)\n"); exit(3); }
        w->t = nt;
    }
    w->t[w->nt] = t;
    return w->nt++;
}

void world_init(World* w, const Params* p, uint64_t seed, const char* rule){
    memset(w, 0, sizeof *w);
    w->seed = seed; w->p = *p; w->tick = 0;
    snprintf(w->rule, sizeof w->rule, "%s", rule);
    snprintf(w->build, sizeof w->build, "%s", pid_build_id());
    rng_seed(&w->rng, seed);
    const Disk s = {0, 0.5 * p->W, 0.5 * p->H, p->R0, -1, KIND_SEED, 0, 1};
    world_add_disk(w, s);
    w->m_used = PID_PI * p->R0 * p->R0;
    const double rot = 2.0 * PID_PI * rng_uniform(&w->rng);   /* one draw: rotates the whole set of K headings */
    for (int k = 0; k < p->K; ++k) {
        const Tip t = {0, pid_wrap(rot + 2.0 * PID_PI * (double)k / (double)p->K), 1};
        world_add_tip(w, t);
    }
}
void world_free(World* w){ free(w->d); free(w->t); w->d = NULL; w->t = NULL; w->nd = w->nt = w->cap_d = w->cap_t = 0; }

/* ---------------------------------------------------------------- geometry */
int pid_disks_overlap(const Disk* a, const Disk* b){
    return hypot(a->x - b->x, a->y - b->y) < a->r + b->r;
}
int pid_disk_in_box(const Disk* d, const Params* p){
    return d->x >= d->r && d->x <= p->W - d->r && d->y >= d->r && d->y <= p->H - d->r;
}
/* Area of the intersection of two circles at centre distance dist (exact; 0 when they at most touch). */
double pid_lens_area(double dist, double r1, double r2){
    if (dist >= r1 + r2) return 0.0;
    const double rs = r1 < r2 ? r1 : r2, rl = r1 < r2 ? r2 : r1;
    if (dist <= rl - rs) return PID_PI * rs * rs;
    /* clamp: at a tangency the rounded cosines can exceed 1 by an ulp, and acos would return NaN */
    double c1 = (dist * dist + r1 * r1 - r2 * r2) / (2.0 * dist * r1), c2 = (dist * dist + r2 * r2 - r1 * r1) / (2.0 * dist * r2);
    c1 = c1 > 1.0 ? 1.0 : (c1 < -1.0 ? -1.0 : c1);
    c2 = c2 > 1.0 ? 1.0 : (c2 < -1.0 ? -1.0 : c2);
    const double a1 = acos(c1), a2 = acos(c2);
    const double k = (-dist + r1 + r2) * (dist + r1 - r2) * (dist - r1 + r2) * (dist + r1 + r2);
    return r1 * r1 * a1 + r2 * r2 * a2 - 0.5 * sqrt(k > 0 ? k : 0);
}

/* ---------------------------------------------------------------- the unbiased rule (0a) */
Action policy_unbiased(const World* w, const Tip* t, Rng* g, const Params* p){
    Action a = {ACT_EXTEND, 0.0, 0.0, 0.0};
    double r = w->d[t->disk].r * (1.0 - p->tau);
    if (r < p->rmin) r = p->rmin;                      /* r' = max(r_min, r (1 - tau)) */
    const double xi = rng_normal(g);                   /* two draws */
    a.theta = pid_wrap(t->heading + p->kappa * xi);
    a.r = r;
    if (rng_uniform(g) < p->pbranch) {                 /* one draw */
        const double s = rng_uniform(g) < 0.5 ? 1.0 : -1.0;   /* one more draw, only when branching */
        a.type = ACT_BRANCH;
        a.theta2 = pid_wrap(a.theta + s * p->delta);
    }
    return a;
}
Policy policy_by_name(const char* name){
    if (!strcmp(name, "unbiased")) return policy_unbiased;
    return NULL;                                       /* "guided" arrives with Demo 0b */
}

/* ---------------------------------------------------------------- executor: legality, material, placement */
static int placement_legal(const World* w, const Disk* nd, int parent){
    if (!pid_disk_in_box(nd, &w->p)) return 0;
    for (int j = 0; j < w->nd; ++j) {
        if (j == parent || !w->d[j].alive) continue;
        if (pid_disks_overlap(nd, &w->d[j])) return 0;
    }
    return 1;
}

int executor_apply(World* w, int ti, Policy pol){
    for (int attempt = 0; attempt <= w->p.ntry; ++attempt) {
        const Action a = pol(w, &w->t[ti], &w->rng, &w->p);
        if (a.type == ACT_STOP) { w->t[ti].alive = 0; return 0; }
        const double area = PID_PI * a.r * a.r;
        if (world_pool(w) < area) { w->t[ti].alive = 0; return 0; }          /* budget exhausted */
        const int pid = w->t[ti].disk;
        const Disk* par = &w->d[pid];
        Disk nd = {0, par->x + (par->r + a.r) * cos(a.theta), par->y + (par->r + a.r) * sin(a.theta), a.r, pid, KIND_TIP, w->tick, 1};
        if (!placement_legal(w, &nd, pid)) continue;                         /* redraw */
        const int id = world_add_disk(w, nd);
        if (w->d[pid].kind == KIND_TIP) w->d[pid].kind = KIND_BODY;
        w->m_used += area;
        w->t[ti].disk = id; w->t[ti].heading = a.theta;
        if (a.type == ACT_BRANCH) { const Tip b = {id, a.theta2, 1}; world_add_tip(w, b); }
        return 1;
    }
    w->t[ti].alive = 0;                                                      /* 1 + ntry attempts failed */
    return 0;
}

void step(World* w, Policy pol){
    w->tick += 1;                                      /* the tick being processed; disks placed now are born at it */
    const int n = w->nt;                               /* tips created in this tick act from the next one */
    for (int i = 0; i < n; ++i) if (w->t[i].alive) executor_apply(w, i, pol);
}

/* ---------------------------------------------------------------- ledger and checks */
double world_pool(const World* w){ return w->p.budget - w->m_used; }

double ledger_recompute(const World* w){
    double s = 0.0;
    for (int i = 0; i < w->nd; ++i) if (w->d[i].alive) s += PID_PI * w->d[i].r * w->d[i].r;
    for (int i = 0; i < w->nd; ++i) {
        const Disk* d = &w->d[i];
        if (!d->alive || d->parent < 0 || !w->d[d->parent].alive) continue;
        const Disk* q = &w->d[d->parent];
        s -= pid_lens_area(hypot(d->x - q->x, d->y - q->y), d->r, q->r);
    }
    return s;
}

int world_live_tips(const World* w){ int k = 0; for (int i = 0; i < w->nt; ++i) k += w->t[i].alive != 0; return k; }

#define CHECK(cond, ...) do { if (!(cond)) { if (!fails) snprintf(err, n, __VA_ARGS__); ++fails; } } while (0)

int world_check(const World* w, int tangency, double m_pool_file, char* err, size_t n){
    int fails = 0; err[0] = '\0';
    CHECK(w->nd >= 1 && w->d[0].parent == -1 && w->d[0].kind == KIND_SEED, "disk 0 is not the seed");
    int* children = (int*)calloc((size_t)(w->nd > 0 ? w->nd : 1), sizeof(int));
    int* tips_on  = (int*)calloc((size_t)(w->nd > 0 ? w->nd : 1), sizeof(int));
    if (!children || !tips_on) { fprintf(stderr, "pid: out of memory (check)\n"); exit(3); }
    for (int i = 0; i < w->nd; ++i) {
        const Disk* d = &w->d[i];
        CHECK(d->id == i, "disk %d has id %d", i, d->id);
        CHECK(pid_disk_in_box(d, &w->p), "disk %d is not inside the box", i);
        if (i == 0) continue;
        CHECK(d->kind != KIND_SEED, "disk %d is a second seed", i);
        CHECK(d->parent >= 0 && d->parent < i, "disk %d: parent %d is not an earlier disk (tree rooted at the seed)", i, d->parent);
        if (d->parent < 0 || d->parent >= i) continue;
        const Disk* q = &w->d[d->parent];
        children[d->parent]++;
        CHECK(d->born > q->born, "disk %d born at %ld, not after its parent %d (%ld)", i, d->born, d->parent, q->born);
        if (tangency) {
            const double dist = hypot(d->x - q->x, d->y - q->y), rs = d->r + q->r;
            CHECK(fabs(dist - rs) <= 1e-12 * rs, "disk %d is not tangent to its parent: |dist - (ri + rj)| = %.3g", i, fabs(dist - rs));
        }
    }
    for (int i = 0; i < w->nd; ++i)                    /* non-adjacent pairs never overlap */
        for (int j = i + 1; j < w->nd; ++j) {
            const Disk* a = &w->d[i]; const Disk* b = &w->d[j];
            if (!a->alive || !b->alive || a->parent == j || b->parent == i) continue;
            CHECK(!pid_disks_overlap(a, b), "disks %d and %d overlap and are not parent and child", i, j);
        }
    for (int k = 0; k < w->nt; ++k) {
        CHECK(w->t[k].disk >= 0 && w->t[k].disk < w->nd, "tip %d sits on no disk", k);
        if (w->t[k].disk >= 0 && w->t[k].disk < w->nd) tips_on[w->t[k].disk]++;
    }
    for (int i = 0; i < w->nd; ++i) {                  /* tip slots: dead tips have no further children */
        const int slots = children[i] + tips_on[i];
        if (i == 0) CHECK(slots == w->p.K, "seed: children %d + tips %d != K %d", children[i], tips_on[i], w->p.K);
        else        CHECK(slots == 1 || slots == 2, "disk %d: children %d + tips %d is not 1 or 2", i, children[i], tips_on[i]);
        if (i > 0) CHECK((w->d[i].kind == KIND_BODY) == (children[i] > 0), "disk %d: kind does not match its children", i);
    }
    const double led = ledger_recompute(w);
    CHECK(fabs(led - w->m_used) <= 1e-9 * fabs(led), "ledger: recomputed M_used %.17g != running %.17g", led, w->m_used);
    CHECK(fabs(w->m_used + m_pool_file - w->p.budget) <= 1e-9 * w->p.budget, "ledger: M_used + M_pool != M_total");
    free(children); free(tips_on);
    return fails;
}
