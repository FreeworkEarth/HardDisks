/* pid_hull.c -- Demo 0c: the exact hull (brief sec. 3 "0c", sec. 6 test 6, sec. 9).
   The lumen is the union of the capsules of all live parent-child edges; a capsule is the convex hull of the two disks.
   Its boundary consists of pieces of the disks' circles and of the capsules' external tangent segments, so:
   1. candidates: one full circle per disk that belongs to a live edge (once, however many capsules share the disk: the
      siblings at a branch point and the six capsules at the hex seed share circles, and arcs are never duplicated), and
      the two tangent segments of every capsule (none when one disk contains the other: that capsule is the larger disk);
   2. every candidate is cut at its intersections with every other candidate and at its own tangent points;
   3. a piece is kept when its midpoint is not strictly inside any capsule (exact signed distance of the hull of two disks,
      tolerance 1e-9 sigma), oriented with the lumen on its left (arcs counter-clockwise, segments as below);
   4. the kept pieces are chained into closed loops (outer boundaries counter-clockwise, holes clockwise);
   5. area = 1/2 of the sum of the Green integrals of x dy - y dx over every piece (segments: the shoelace term; arcs: the
      closed form), perimeter = the sum of the lengths. Both are sums over the pieces, so they do not depend on the chaining.
   Edges are processed in ascending child id whatever order they are given in, so the result does not depend on it. */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "pid.h"

#define HULL_EPS 1e-9      /* inside tolerance (sigma), and the endpoint match of the chaining */

typedef struct {
    int a, b, whole;       /* disk ids; whole: 0 proper capsule, 1 the capsule is disk a, 2 it is disk b */
    double x1, y1, r1, x2, y2, r2, ux, uy, h;
    double bx0, by0, bx1, by1;
} Cap;

typedef struct {
    int kind;              /* 0 circle of disk `disk`, 1 segment of capsule `cap` */
    int disk, cap;
    double cx, cy, R;      /* circle */
    double x1, y1, x2, y2; /* segment, lumen on the left when walking from 1 to 2 */
    double bx0, by0, bx1, by1;
    double* cut; int ncut, capcut;   /* angles (circle) or parameters t in (0, 1) (segment) */
} Cand;

static void add_cut(Cand* c, double v){
    if (c->ncut == c->capcut) {
        c->capcut = c->capcut ? 2 * c->capcut : 8;
        c->cut = (double*)realloc(c->cut, (size_t)c->capcut * sizeof(double));
        if (!c->cut) { fprintf(stderr, "pid: out of memory (hull)\n"); exit(3); }
    }
    c->cut[c->ncut++] = v;
}

/* signed distance to the convex hull of two disks (Inigo Quilez's uneven capsule, in the frame of the axis) */
static double cap_sdf(const Cap* c, double px, double py){
    if (c->whole == 1) return hypot(px - c->x1, py - c->y1) - c->r1;
    if (c->whole == 2) return hypot(px - c->x2, py - c->y2) - c->r2;
    const double dx = px - c->x1, dy = py - c->y1;
    const double along = dx * c->ux + dy * c->uy, lat = fabs(dx * c->uy - dy * c->ux);
    const double b = (c->r1 - c->r2) / c->h, a = sqrt(1.0 - b * b);
    const double k = -b * lat + a * along;
    if (k < 0.0) return hypot(lat, along) - c->r1;
    if (k > a * c->h) return hypot(lat, along - c->h) - c->r2;
    return a * lat + b * along - c->r1;
}

static double wrap2pi(double a){ while (a < 0.0) a += 2.0 * PID_PI; while (a >= 2.0 * PID_PI) a -= 2.0 * PID_PI; return a; }
static int cmp_d(const void* x, const void* y){ const double a = *(const double*)x, b = *(const double*)y; return a < b ? -1 : a > b; }
static int cmp_i(const void* x, const void* y){ return *(const int*)x - *(const int*)y; }
typedef struct { double x; int k; } XKey;
static int cmp_xk(const void* p, const void* q){
    const XKey* a = (const XKey*)p; const XKey* b = (const XKey*)q;
    return a->x < b->x ? -1 : a->x > b->x ? 1 : a->k - b->k;
}
static int bb_overlap(double a0, double b0, double a1, double b1, double c0, double d0, double c1, double d1){
    return !(a1 < c0 || c1 < a0 || b1 < d0 || d1 < b0);
}

void hull_add(Hull* H, int type, double v0, double v1, double v2, double v3, double v4){
    if (H->np == H->cap) {
        H->cap = H->cap ? 2 * H->cap : 256;
        H->p = (HPrim*)realloc(H->p, (size_t)H->cap * sizeof(HPrim));
        if (!H->p) { fprintf(stderr, "pid: out of memory (hull)\n"); exit(3); }
    }
    HPrim* q = &H->p[H->np++]; q->type = type; q->v[0] = v0; q->v[1] = v1; q->v[2] = v2; q->v[3] = v3; q->v[4] = v4;
}
void hull_add_loop(Hull* H, int n){
    if (H->nloop == H->cap_loop) {
        H->cap_loop = H->cap_loop ? 2 * H->cap_loop : 32;
        H->loop = (int*)realloc(H->loop, (size_t)H->cap_loop * sizeof(int));
        if (!H->loop) { fprintf(stderr, "pid: out of memory (hull)\n"); exit(3); }
    }
    H->loop[H->nloop++] = n;
}
void hull_free(Hull* H){ free(H->p); free(H->loop); memset(H, 0, sizeof *H); }

void hull_prim_ends(const HPrim* q, double* xs, double* ys, double* xe, double* ye){
    if (q->type == 0) { *xs = q->v[0]; *ys = q->v[1]; *xe = q->v[2]; *ye = q->v[3]; return; }
    *xs = q->v[0] + q->v[2] * cos(q->v[3]); *ys = q->v[1] + q->v[2] * sin(q->v[3]);
    *xe = q->v[0] + q->v[2] * cos(q->v[4]); *ye = q->v[1] + q->v[2] * sin(q->v[4]);
}
static double prim_green(const HPrim* q){    /* 1/2 of the integral of x dy - y dx along the piece */
    if (q->type == 0) return 0.5 * (q->v[0] * q->v[3] - q->v[2] * q->v[1]);
    const double cx = q->v[0], cy = q->v[1], R = q->v[2], a0 = q->v[3], a1 = q->v[4];
    return 0.5 * (R * cx * (sin(a1) - sin(a0)) - R * cy * (cos(a1) - cos(a0)) + R * R * (a1 - a0));
}
static double prim_len(const HPrim* q){
    if (q->type == 0) return hypot(q->v[2] - q->v[0], q->v[3] - q->v[1]);
    return q->v[2] * (q->v[4] - q->v[3]);
}

/* ---------------------------------------------------------------- intersections */
static void cut_circle_circle(Cand* A, Cand* B){
    const double dx = B->cx - A->cx, dy = B->cy - A->cy, d = hypot(dx, dy);
    if (!(d < A->R + B->R) || !(d > fabs(A->R - B->R)) || d == 0.0) return;   /* apart, touching or nested */
    const double a = (A->R * A->R - B->R * B->R + d * d) / (2.0 * d), hh2 = A->R * A->R - a * a;
    if (!(hh2 > 0.0)) return;
    const double hh = sqrt(hh2), mx = A->cx + a * dx / d, my = A->cy + a * dy / d;
    for (int s = -1; s <= 1; s += 2) {
        const double px = mx - s * hh * dy / d, py = my + s * hh * dx / d;
        add_cut(A, wrap2pi(atan2(py - A->cy, px - A->cx)));
        add_cut(B, wrap2pi(atan2(py - B->cy, px - B->cx)));
    }
}
static void cut_seg_circle(Cand* S, Cand* C){
    const double vx = S->x2 - S->x1, vy = S->y2 - S->y1, fx = S->x1 - C->cx, fy = S->y1 - C->cy;
    const double a = vx * vx + vy * vy, b = 2.0 * (fx * vx + fy * vy), c = fx * fx + fy * fy - C->R * C->R;
    const double disc = b * b - 4.0 * a * c;
    if (!(disc > 0.0) || a == 0.0) return;
    const double sq = sqrt(disc);
    for (int s = -1; s <= 1; s += 2) {
        const double t = (-b + s * sq) / (2.0 * a);
        if (!(t > 0.0 && t < 1.0)) continue;
        const double px = S->x1 + t * vx, py = S->y1 + t * vy;
        add_cut(S, t);
        add_cut(C, wrap2pi(atan2(py - C->cy, px - C->cx)));
    }
}
static void cut_seg_seg(Cand* P, Cand* Q){
    const double rx = P->x2 - P->x1, ry = P->y2 - P->y1, sx = Q->x2 - Q->x1, sy = Q->y2 - Q->y1;
    const double den = rx * sy - ry * sx;
    if (den == 0.0) return;                                  /* parallel or collinear: no proper crossing */
    const double qx = Q->x1 - P->x1, qy = Q->y1 - P->y1;
    const double t = (qx * sy - qy * sx) / den, u = (qx * ry - qy * rx) / den;
    if (!(t > 0.0 && t < 1.0 && u > 0.0 && u < 1.0)) return;
    add_cut(P, t); add_cut(Q, u);
}

/* ---------------------------------------------------------------- the hull */
static int inside_any(const Cap* caps, int nc, double px, double py){
    for (int k = 0; k < nc; ++k) {
        const Cap* c = &caps[k];
        if (px < c->bx0 || px > c->bx1 || py < c->by0 || py > c->by1) continue;
        if (cap_sdf(c, px, py) < -HULL_EPS) return 1;
    }
    return 0;
}

int hull_compute(const World* w, const int* edge_child, int ne, Hull* H){
    memset(H, 0, sizeof *H);
    int* kids = (int*)malloc((size_t)(ne > 0 ? ne : 1) * sizeof(int));
    if (!kids) { fprintf(stderr, "pid: out of memory (hull)\n"); exit(3); }
    memcpy(kids, edge_child, (size_t)ne * sizeof(int));
    qsort(kids, (size_t)ne, sizeof(int), cmp_i);             /* canonical order: ascending child id */
    Cap* caps = (Cap*)calloc((size_t)(ne > 0 ? ne : 1), sizeof(Cap));
    int* circ_of = (int*)malloc((size_t)(w->nd > 0 ? w->nd : 1) * sizeof(int));
    Cand* cand = (Cand*)calloc((size_t)(3 * ne + w->nd + 1), sizeof(Cand));
    if (!caps || !circ_of || !cand) { fprintf(stderr, "pid: out of memory (hull)\n"); exit(3); }
    for (int i = 0; i < w->nd; ++i) circ_of[i] = -1;
    int nc = 0, nk = 0;
    for (int e = 0; e < ne; ++e) {                           /* capsules, and the circles of their disks */
        const int a = kids[e], b = w->d[a].parent;
        const Disk* A = &w->d[a]; const Disk* B = &w->d[b];
        Cap* c = &caps[nc++];
        c->a = a; c->b = b; c->x1 = A->x; c->y1 = A->y; c->r1 = A->r; c->x2 = B->x; c->y2 = B->y; c->r2 = B->r;
        const double dx = B->x - A->x, dy = B->y - A->y; c->h = hypot(dx, dy);
        if (c->h <= fabs(A->r - B->r)) c->whole = A->r >= B->r ? 1 : 2;
        else { c->whole = 0; c->ux = dx / c->h; c->uy = dy / c->h; }
        c->bx0 = fmin(A->x - A->r, B->x - B->r); c->bx1 = fmax(A->x + A->r, B->x + B->r);
        c->by0 = fmin(A->y - A->r, B->y - B->r); c->by1 = fmax(A->y + A->r, B->y + B->r);
    }
    for (int k = 0; k < nc; ++k) { circ_of[caps[k].a] = 1; circ_of[caps[k].b] = 1; }   /* disks that belong to an edge */
    for (int i = 0; i < w->nd; ++i) {                        /* their circles, numbered in ascending disk id */
        if (circ_of[i] < 0) continue;
        Cand* C = &cand[nk]; circ_of[i] = nk++;
        C->kind = 0; C->disk = i; C->cx = w->d[i].x; C->cy = w->d[i].y; C->R = w->d[i].r;
        C->bx0 = C->cx - C->R; C->bx1 = C->cx + C->R; C->by0 = C->cy - C->R; C->by1 = C->cy + C->R;
    }
    for (int k = 0; k < nc; ++k) {                           /* tangent segments, and the tangent points on the circles */
        const Cap* c = &caps[k]; if (c->whole) continue;
        const double phi = atan2(c->uy, c->ux), al = acos((c->r1 - c->r2) / c->h);
        const double tp = phi + al, tm = phi - al;
        const double p1px = c->x1 + c->r1 * cos(tp), p1py = c->y1 + c->r1 * sin(tp), p1mx = c->x1 + c->r1 * cos(tm), p1my = c->y1 + c->r1 * sin(tm);
        const double p2px = c->x2 + c->r2 * cos(tp), p2py = c->y2 + c->r2 * sin(tp), p2mx = c->x2 + c->r2 * cos(tm), p2my = c->y2 + c->r2 * sin(tm);
        const double seg[2][4] = {{p1mx, p1my, p2mx, p2my}, {p2px, p2py, p1px, p1py}};   /* lower 1->2, upper 2->1 */
        for (int s = 0; s < 2; ++s) {
            Cand* S = &cand[nk++];
            S->kind = 1; S->cap = k; S->disk = -1;
            S->x1 = seg[s][0]; S->y1 = seg[s][1]; S->x2 = seg[s][2]; S->y2 = seg[s][3];
            S->bx0 = fmin(S->x1, S->x2); S->bx1 = fmax(S->x1, S->x2); S->by0 = fmin(S->y1, S->y2); S->by1 = fmax(S->y1, S->y2);
        }
        add_cut(&cand[circ_of[c->a]], wrap2pi(tp)); add_cut(&cand[circ_of[c->a]], wrap2pi(tm));
        add_cut(&cand[circ_of[c->b]], wrap2pi(tp)); add_cut(&cand[circ_of[c->b]], wrap2pi(tm));
    }
    for (int i = 0; i < nk; ++i)                             /* every pair of candidates, pruned by bounding boxes */
        for (int j = i + 1; j < nk; ++j) {
            Cand* P = &cand[i]; Cand* Q = &cand[j];
            if (!bb_overlap(P->bx0, P->by0, P->bx1, P->by1, Q->bx0, Q->by0, Q->bx1, Q->by1)) continue;
            if (P->kind == 0 && Q->kind == 0) cut_circle_circle(P, Q);
            else if (P->kind == 1 && Q->kind == 1) { if (P->cap != Q->cap) cut_seg_seg(P, Q); }
            else {
                Cand* S = P->kind == 1 ? P : Q; Cand* C = P->kind == 1 ? Q : P;
                if (C->disk == caps[S->cap].a || C->disk == caps[S->cap].b) continue;   /* its own circles: the tangent points */
                cut_seg_circle(S, C);
            }
        }
    /* pieces: keep those whose midpoint is not strictly inside any capsule */
    Hull pieces; memset(&pieces, 0, sizeof pieces);
    for (int i = 0; i < nk; ++i) {
        Cand* C = &cand[i];
        if (C->kind == 0) {
            if (C->ncut == 0) { if (!inside_any(caps, nc, C->cx + C->R, C->cy)) hull_add(&pieces, 1, C->cx, C->cy, C->R, -PID_PI, PID_PI); continue; }
            qsort(C->cut, (size_t)C->ncut, sizeof(double), cmp_d);
            for (int k = 0; k < C->ncut; ++k) {
                const double a0 = C->cut[k], a1 = k + 1 < C->ncut ? C->cut[k + 1] : C->cut[0] + 2.0 * PID_PI;
                if (!(a1 - a0 > 1e-15)) continue;
                const double am = 0.5 * (a0 + a1);
                if (inside_any(caps, nc, C->cx + C->R * cos(am), C->cy + C->R * sin(am))) continue;
                const double s0 = a0 > PID_PI ? a0 - 2.0 * PID_PI : a0;   /* a0 in (-pi, pi], a1 = a0 + span */
                hull_add(&pieces, 1, C->cx, C->cy, C->R, s0, s0 + (a1 - a0));
            }
        } else {
            add_cut(C, 0.0); add_cut(C, 1.0);
            qsort(C->cut, (size_t)C->ncut, sizeof(double), cmp_d);
            const double vx = C->x2 - C->x1, vy = C->y2 - C->y1;
            for (int k = 0; k + 1 < C->ncut; ++k) {
                const double t0 = C->cut[k], t1 = C->cut[k + 1];
                if (!(t1 - t0 > 1e-15)) continue;
                const double tm = 0.5 * (t0 + t1);
                if (inside_any(caps, nc, C->x1 + tm * vx, C->y1 + tm * vy)) continue;
                const double xa = t0 == 0.0 ? C->x1 : C->x1 + t0 * vx, ya = t0 == 0.0 ? C->y1 : C->y1 + t0 * vy;
                const double xb = t1 == 1.0 ? C->x2 : C->x1 + t1 * vx, yb = t1 == 1.0 ? C->y2 : C->y1 + t1 * vy;
                hull_add(&pieces, 0, xa, ya, xb, yb, 0.0);
            }
        }
    }
    H->area = 0.0; H->perim = 0.0;
    for (int k = 0; k < pieces.np; ++k) { H->area += prim_green(&pieces.p[k]); H->perim += prim_len(&pieces.p[k]); }
    /* chain into closed loops: from the end of a piece, the unused piece whose start is nearest (within HULL_EPS) */
    const int np = pieces.np;
    double* sx = (double*)malloc((size_t)(np > 0 ? np : 1) * 4 * sizeof(double));
    char* used = (char*)calloc((size_t)(np > 0 ? np : 1), 1);
    int* order = (int*)malloc((size_t)(np > 0 ? np : 1) * sizeof(int));
    if (!sx || !used || !order) { fprintf(stderr, "pid: out of memory (hull)\n"); exit(3); }
    XKey* xk = (XKey*)malloc((size_t)(np > 0 ? np : 1) * sizeof(XKey));
    if (!xk) { fprintf(stderr, "pid: out of memory (hull)\n"); exit(3); }
    for (int k = 0; k < np; ++k) { hull_prim_ends(&pieces.p[k], &sx[4 * k], &sx[4 * k + 1], &sx[4 * k + 2], &sx[4 * k + 3]); xk[k].x = sx[4 * k]; xk[k].k = k; }
    qsort(xk, (size_t)np, sizeof(XKey), cmp_xk);             /* order by start x for the search */
    for (int k = 0; k < np; ++k) order[k] = xk[k].k;
    free(xk);
    for (int s0 = 0; s0 < np; ++s0) {
        if (used[s0]) continue;
        int n = 0, cur = s0; used[cur] = 1; hull_add(H, pieces.p[cur].type, pieces.p[cur].v[0], pieces.p[cur].v[1], pieces.p[cur].v[2], pieces.p[cur].v[3], pieces.p[cur].v[4]); ++n;
        for (;;) {
            const double ex = sx[4 * cur + 2], ey = sx[4 * cur + 3];
            if (hypot(ex - sx[4 * s0], ey - sx[4 * s0 + 1]) <= HULL_EPS) break;   /* back at the start: the loop is closed */
            int lo = 0, hi = np;                             /* first in `order` with start x >= ex - eps */
            while (lo < hi) { const int m = (lo + hi) / 2; if (sx[4 * order[m]] < ex - HULL_EPS) lo = m + 1; else hi = m; }
            int best = -1; double bd = HULL_EPS;
            for (int m = lo; m < np && sx[4 * order[m]] <= ex + HULL_EPS; ++m) {
                const int q = order[m]; if (used[q]) continue;
                const double dd = hypot(sx[4 * q] - ex, sx[4 * q + 1] - ey);
                if (dd <= bd) { bd = dd; best = q; }
            }
            if (best < 0) { H->open++; break; }              /* should not happen: counted, and test (a) fails */
            cur = best; used[cur] = 1; ++n;
            hull_add(H, pieces.p[cur].type, pieces.p[cur].v[0], pieces.p[cur].v[1], pieces.p[cur].v[2], pieces.p[cur].v[3], pieces.p[cur].v[4]);
        }
        hull_add_loop(H, n);
    }
    free(sx); free(used); free(order); hull_free(&pieces);
    for (int i = 0; i < nk; ++i) free(cand[i].cut);
    free(cand); free(caps); free(circ_of); free(kids);
    H->valid = 1;
    return 0;
}

int hull_edges(const World* w, int* edge_child){             /* the live edges, as child ids, ascending */
    int ne = 0;
    for (int i = 1; i < w->nd; ++i) if (w->d[i].alive && w->d[i].parent >= 0 && w->d[w->d[i].parent].alive) edge_child[ne++] = i;
    return ne;
}
int hull_compute_world(const World* w, Hull* H){
    int* e = (int*)malloc((size_t)(w->nd > 0 ? w->nd : 1) * sizeof(int));
    if (!e) { fprintf(stderr, "pid: out of memory (hull)\n"); exit(3); }
    const int ne = hull_edges(w, e);
    hull_compute(w, e, ne, H);
    free(e);
    return 0;
}

/* ---------------------------------------------------------------- checks for test 6 */
double hull_max_gap(const Hull* H){                          /* (a): end of each primitive to the start of the next */
    double g = 0.0; int k0 = 0;
    for (int l = 0; l < H->nloop; ++l) {
        for (int k = 0; k < H->loop[l]; ++k) {
            const int a = k0 + k, b = k0 + (k + 1) % H->loop[l];
            double xs, ys, xe, ye, xs2, ys2, xe2, ye2;
            hull_prim_ends(&H->p[a], &xs, &ys, &xe, &ye); hull_prim_ends(&H->p[b], &xs2, &ys2, &xe2, &ye2);
            const double d = hypot(xe - xs2, ye - ys2); if (d > g) g = d;
        }
        k0 += H->loop[l];
    }
    return g + (H->open ? 1e300 : 0.0);
}
/* (b): crossing number of a ray from (x, y) in a fixed direction; half-open rules at the endpoints */
int hull_contains(const Hull* H, double x, double y){
    const double dx = cos(0.6180339887498949), dy = sin(0.6180339887498949);
    int cross = 0;
    for (int k = 0; k < H->np; ++k) {
        const HPrim* q = &H->p[k];
        if (q->type == 0) {
            const double sx = q->v[2] - q->v[0], sy = q->v[3] - q->v[1], den = dx * sy - dy * sx;
            if (den == 0.0) continue;
            const double qx = q->v[0] - x, qy = q->v[1] - y;
            const double t = (qx * sy - qy * sx) / den, u = (qx * dy - qy * dx) / den;
            if (t > 0.0 && u >= 0.0 && u < 1.0) cross++;
        } else {
            const double fx = x - q->v[0], fy = y - q->v[1], R = q->v[2];
            const double b = 2.0 * (fx * dx + fy * dy), c = fx * fx + fy * fy - R * R, disc = b * b - 4.0 * c;
            if (!(disc > 0.0)) continue;
            for (int s = -1; s <= 1; s += 2) {
                const double t = (-b + s * sqrt(disc)) / 2.0;
                if (!(t > 0.0)) continue;
                double a = atan2(y + t * dy - q->v[1], x + t * dx - q->v[0]);
                while (a < q->v[3]) a += 2.0 * PID_PI;
                while (a >= q->v[3] + 2.0 * PID_PI) a -= 2.0 * PID_PI;
                if (a < q->v[4]) cross++;
            }
        }
    }
    return cross & 1;
}
/* (c): area of the capsule union on a raster of `res` pixels per sigma (pixel centres) */
double hull_raster_area(const World* w, double res){
    const int nx = (int)ceil(w->p.W * res), ny = (int)ceil(w->p.H * res);
    unsigned char* m = (unsigned char*)calloc((size_t)nx * (size_t)ny, 1);
    int* e = (int*)malloc((size_t)(w->nd > 0 ? w->nd : 1) * sizeof(int));
    if (!m || !e) { fprintf(stderr, "pid: out of memory (raster)\n"); exit(3); }
    const int ne = hull_edges(w, e);
    long count = 0;
    for (int k = 0; k < ne; ++k) {
        const Disk* A = &w->d[e[k]]; const Disk* B = &w->d[A->parent];
        Cap c; memset(&c, 0, sizeof c);
        c.x1 = A->x; c.y1 = A->y; c.r1 = A->r; c.x2 = B->x; c.y2 = B->y; c.r2 = B->r;
        const double ddx = B->x - A->x, ddy = B->y - A->y; c.h = hypot(ddx, ddy);
        if (c.h <= fabs(A->r - B->r)) c.whole = A->r >= B->r ? 1 : 2; else { c.ux = ddx / c.h; c.uy = ddy / c.h; }
        int i0 = (int)floor(fmin(A->x - A->r, B->x - B->r) * res), i1 = (int)ceil(fmax(A->x + A->r, B->x + B->r) * res);
        int j0 = (int)floor(fmin(A->y - A->r, B->y - B->r) * res), j1 = (int)ceil(fmax(A->y + A->r, B->y + B->r) * res);
        if (i0 < 0) i0 = 0;
        if (j0 < 0) j0 = 0;
        if (i1 > nx) i1 = nx;
        if (j1 > ny) j1 = ny;
        for (int j = j0; j < j1; ++j)
            for (int i = i0; i < i1; ++i) {
                unsigned char* px = &m[(size_t)j * (size_t)nx + (size_t)i];
                if (*px) continue;
                if (cap_sdf(&c, (i + 0.5) / res, (j + 0.5) / res) < 0.0) { *px = 1; ++count; }
            }
    }
    free(m); free(e);
    return (double)count / (res * res);
}
