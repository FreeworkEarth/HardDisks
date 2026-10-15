/* ##CHRIS 2026-10-09 (M3, 261012 sec. 4.7.14, amendment a): the numbers behind every schedule-audit finding of the band-edge
   cells (gen3_band_edge_cells.h), taken at the moment the audit finds it (white box: this file includes edmd_gen3.c and sets
   its EDMD3_AUDIT_HOOK). The runs are those of gen3_band_edge_test.c: same cells, same advance targets, same body changes,
   same audits. Per distinct finding (kind, class, a, b), the first occurrence is shown and the occurrences are counted.
   Classes, each decided from the doubles the engine and the audit actually hold (no tolerance enters a class):
     G  exact graze: the relative velocity has no component along one axis (dv_y = 0) and the separation along that axis is
        exactly the contact distance (|r_y| = d), or x and y swapped. In exact arithmetic the discriminant b^2 - |dv|^2 c is
        then dv_x^2 (d^2 - r_y^2) = 0: the disks touch with zero normal speed (zero impulse). Its sign in floating point is
        rounding: the engine (local coordinates) and the audit (absolute coordinates) may disagree on whether the touch
        happens, and when (the root of a zero discriminant is ill-conditioned).
     K  near graze, the discriminant zero within rounding: no exact tie in the doubles, but the evaluation that finds the event
        (or both, for a time disagreement) has 0 < disc <= 4 eps b^2 (eps = 2^-52), so its sign is rounding; or, for a time
        disagreement, |t_heap - t_bruteforce| <= 4 eps b^2 / (2 sqrt(disc) |dv|^2), the shift of the root by an error of 4 eps b^2
        in disc (the root of a near-zero discriminant is ill-conditioned). (Class K was added after the first run of this file,
        in which one finding, the spring cell's pair 12/18 time disagreement, was U: r_x = d - 3.3e-12 px, dv_x = -2.2e-14
        px/unit, disc = 1.5e-11 against b^2 = 81.)
     V  relative velocity at rounding level: |dv| <= 1e-12 px/unit, two velocities equal in exact arithmetic; a contact, if
        any, lies >= 1e12 units ahead.
     C  corner tie: the x and y crossing times of the disk agree to within 8 ulp of the larger one: the disk reaches a cell
        corner; which crossing comes first is rounding (both orders lead to the same cell).
     U  unclassified: none of the above (a finding that is not a tie; counted, must be 0).
   Also printed for pairs: whether the two disks are on opposite sides of the divider (the brute-force contact is then
   fictitious: the audit's pair rule is ballistic and does not know the divider).
   build (from hspist3/): cc -std=c11 -O2 -ffp-contract=off -Wall -Wextra -o <out> edmd_core/tests/gen3_band_edge_ties.c -lm */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
struct EDMD3;
static void be_hook(struct EDMD3* S, const char* kind, const char* what, int a, int b, double tbf, double th);
#define EDMD3_AUDIT_HOOK(S, kind, what, a, b, tbf, th) be_hook(S, kind, what, a, b, tbf, th)
#include "../edmd_gen3.c"
#include "gen3_band_edge_cells.h"

typedef struct { char kind[24], what[8]; int a, b; long n; } Seen;
static Seen seen[256]; static int nseen = 0;
static const BECell* cur = NULL;
static long n_occ = 0, n_dist = 0, n_unclass = 0, n_cls[4] = {0, 0, 0, 0};   /* G, K, V, C */

static double heap_time(const EDMD3* S, int type, int a, int b, int* dir){
    double t = NAN;
    for (long q = 0; q < S->heap.n; ++q) {
        const Ev3* e = &S->heap.d[q];
        if (e->type != type || !ev_live(S, e)) continue;
        if (type == T_PAIR && !((e->a == a && e->b == b) || (e->a == b && e->b == a))) continue;
        if (type == T_CROSS && e->a != a) continue;
        if (isnan(t) || e->t < t) { t = e->t; if (dir) *dir = e->b; }
    }
    return t;
}

static void be_hook(struct EDMD3* S, const char* kind, const char* what, int a, int b, double tbf, double th){
    n_occ++;
    for (int k = 0; k < nseen; ++k)
        if (seen[k].a == a && seen[k].b == b && !strcmp(seen[k].kind, kind) && !strcmp(seen[k].what, what)) { seen[k].n++; return; }
    if (nseen < 256) { Seen* z = &seen[nseen++]; snprintf(z->kind, sizeof z->kind, "%s", kind); snprintf(z->what, sizeof z->what, "%s", what); z->a = a; z->b = b; z->n = 1; }
    n_dist++;
    char q[700] = ""; char cls = 'U';
    if (!strcmp(what, "PAIR")) {
        const Disk3 *A = &S->D[a], *B = &S->D[b];
        double ax, ay, bx, by; local_now(S, a, &ax, &ay); local_now(S, b, &bx, &by);
        const double rxl = (bx - ax) + (double)(B->cx - A->cx) * S->w, ryl = (by - ay) + (double)(B->cy - A->cy) * S->w;
        const double rxa = ((double)B->cx * S->w + bx) - ((double)A->cx * S->w + ax), rya = ((double)B->cy * S->w + by) - ((double)A->cy * S->w + ay);
        const double dvx = B->vx - A->vx, dvy = B->vy - A->vy, vv = dvx * dvx + dvy * dvy;
        const double bl = rxl * dvx + ryl * dvy, cl = rxl * rxl + ryl * ryl - S->d2, ba = rxa * dvx + rya * dvy, ca = rxa * rxa + rya * rya - S->d2;
        const double discl = bl * bl - vv * cl, disca = ba * ba - vv * ca;
        /* the contact axis: the separation along it is exactly d (local or absolute); dv along it decides G or N */
        const int ax_y = fabs(ryl) == S->d || fabs(rya) == S->d, ax_x = fabs(rxl) == S->d || fabs(rxa) == S->d;
        const int G = (ax_y && dvy == 0.0 && dvx != 0.0) || (ax_x && dvx == 0.0 && dvy != 0.0);
        const double eps4 = 4.0 * 0x1p-52;
        int K = 0; double kbound = NAN;
        if (!strcmp(kind, "dt_rel")) {
            const double dm = fmin(discl, disca), b2 = fmax(bl * bl, ba * ba);
            if (dm > 0.0) { kbound = eps4 * b2 / (2.0 * sqrt(dm) * vv); K = fabs(tbf - th) <= kbound; }
        } else {   /* missing: the audit (absolute) finds the event; extra: the engine (local) does */
            const double dsc = !strcmp(kind, "missing") ? disca : discl, bb = !strcmp(kind, "missing") ? ba * ba : bl * bl;
            kbound = eps4 * bb; K = dsc > 0.0 && dsc <= kbound;
        }
        const int V = sqrt(vv) <= 1e-12;
        cls = G ? 'G' : (K ? 'K' : (V ? 'V' : 'U'));
        const double xa = (double)A->cx * S->w + ax, xb = (double)B->cx * S->w + bx;
        double xdiv = NAN, vdiv; if (S->nobj) obj_at(&S->obj[S->objs[0]], S->now, &xdiv, &vdiv);
        const int opp = S->nobj && ((xa < xdiv) != (xb < xdiv));
        snprintf(q, sizeof q, "dv = (%.6g, %.6g); r_x local %.17g / absolute %.17g; r_y local %.17g / absolute %.17g; disc local %.3g / absolute %.3g; "
                 "heap %.17g, brute force %.17g; K bound %.3g; opposite sides of the divider: %s",
                 dvx, dvy, rxl, rxa, ryl, rya, discl, disca, heap_time(S, T_PAIR, a, b, NULL), tbf, kbound, opp ? "yes" : "no");
    } else if (!strcmp(what, "CROSS")) {
        const Disk3* D = &S->D[a]; double x, y; local_now(S, a, &x, &y);
        double tx = INFINITY, ty = INFINITY;
        if (D->vx > 0.0) tx = (S->w - x) / D->vx; else if (D->vx < 0.0) tx = x / (-D->vx);
        if (D->vy > 0.0) ty = (S->w - y) / D->vy; else if (D->vy < 0.0) ty = y / (-D->vy);
        /* the engine's prediction basis: the disk's stamp */
        double tx0 = INFINITY, ty0 = INFINITY;
        if (D->vx > 0.0) tx0 = (S->w - D->xi) / D->vx; else if (D->vx < 0.0) tx0 = D->xi / (-D->vx);
        if (D->vy > 0.0) ty0 = (S->w - D->zeta) / D->vy; else if (D->vy < 0.0) ty0 = D->zeta / (-D->vy);
        const double big = fmax(S->now + tx, S->now + ty);
        const int C = isfinite(tx) && isfinite(ty) && fabs((S->now + tx) - (S->now + ty)) <= 8.0 * (nextafter(big, INFINITY) - big);
        cls = C ? 'C' : 'U';
        int hdir = -1; const double ht = heap_time(S, T_CROSS, a, 0, &hdir);
        snprintf(q, sizeof q, "cell (%d,%d), local (%.17g, %.17g), v (%.6g, %.6g); x and y crossing at %.17g / %.17g (from the stamp at "
                 "%.17g: %.17g / %.17g); heap: direction %d at %.17g; audit: direction %d",
                 D->cx, D->cy, x, y, D->vx, D->vy, S->now + tx, S->now + ty, D->tau, D->tau + tx0, D->tau + ty0, hdir, ht, b);
    } else snprintf(q, sizeof q, "brute force %.17g, heap %.17g", tbf, th);
    if (cls == 'U') n_unclass++; else n_cls[cls == 'G' ? 0 : (cls == 'K' ? 1 : (cls == 'V' ? 2 : 3))]++;
    printf("| %s (%s) | %s %s a=%d b=%d | %.17g | %s | %c |\n", cur->name, cur->generic ? "generic y" : "exact", kind, what, a, b, S->now, q, cls);
}

int main(void){
    BECell cells[BE_NCELL]; be_cells(cells);
    printf("# gen3 band-edge cells: the numbers behind each schedule-audit finding (amendment a), printed by "
           "edmd_core/tests/gen3_band_edge_ties.c\n\n"
           "| cell | finding (first occurrence) | at t | quantities | class |\n|---|---|---|---|---|\n");
    long occ_by_cell[BE_NCELL];
    for (int k = 0; k < BE_NCELL; ++k) {
        const BECell* c = &cells[k]; cur = c; nseen = 0; const long occ0 = n_occ;
        EDMD_Particle P[64]; const int N = be_place(c, P);
        EDMD_Params p; be_params(c, N, &p);
        char err[256];
        EDMD3* S = edmd3_create(&p, BE_W, err, sizeof err);
        if (!S || !edmd3_load(S, P, 0.0, err, sizeof err)) { printf("%s: %s\n", c->name, err); return 2; }
        edmd3_set_contact_audit(S, 1); edmd3_set_schedule_audit(S, 1); edmd3_set_schedule_audit_bodies(S, 1);
        S->aprinted = 50;                                   /* the engine's own audit lines are printed by gen3_band_edge_test */
        double v = c->v;
        for (double t = 1.0; t <= c->T + 1e-9; t += 1.0) {
            edmd3_advance_to(S, t);
            if (c->flip > 0.0 && fmod(t, c->flip) == 0.0) { v = -v; edmd3_set_divider_motion(S, 0, 0.0, v); }
            const char* msg; if (edmd3_fatal(S, &msg)) { printf("%s: fatal: %s\n", c->name, msg); break; }
        }
        edmd3_schedule_audit_now(S);
        occ_by_cell[k] = n_occ - occ0;
        for (int z = 0; z < nseen; ++z)
            printf("| %s (%s) | occurrences of %s %s a=%d b=%d: %ld |  |  |  |\n", c->name, c->generic ? "generic y" : "exact",
                   seen[z].kind, seen[z].what, seen[z].a, seen[z].b, seen[z].n);
        edmd3_destroy(S);
    }
    printf("\nfindings (audit lines, all cells): %ld occurrences, %ld distinct; by cell:", n_occ, n_dist);
    for (int k = 0; k < BE_NCELL; ++k) printf(" %s (%s) %ld;", cells[k].name, cells[k].generic ? "generic y" : "exact", occ_by_cell[k]);
    printf("\ndistinct findings by class: G exact graze %ld, K near graze within rounding %ld, V rounding-level relative velocity %ld, C corner tie %ld, "
           "U unclassified %ld\n", n_cls[0], n_cls[1], n_cls[2], n_cls[3], n_unclass);
    return n_unclass ? 1 : 0;
}
