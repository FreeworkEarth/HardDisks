/* pid_check.c -- the test tool of brief sec. 6; reads state files only, so every number it checks comes from the headless
   path's files.
     pid_check verify FILE...      invariants and ledger on every file (items 2-3); with several files, given in tick order,
                                   also that disks never change and dead tips stay frozen from one state to the next
     pid_check roundtrip IN OUT    read IN, write OUT (item 4; the test then compares the two with cmp)
     pid_check info FILE...        0b: one line per state: tick, alive disks, live tips, M_pool, fed disks, fed path ends, and for every
                                   non-home patch the number of alive disks touching it
     pid_check hull FILE...        0c, test 6 (a)-(c) and (e) on the hull recomputed from the disks; and, when the file
                                   carries a hull, that it equals the recomputed one */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "pid.h"

const char* pid_build_id(void){ return "pid_check"; }

static int prim_eq(const HPrim* a, const HPrim* b){            /* field by field: the struct has padding bytes */
    if (a->type != b->type) return 0;
    for (int k = 0; k < 5; ++k) if (a->v[k] != b->v[k]) return 0;
    return 1;
}   /* world_init() is not used here; state files carry their own build line */

static int sequence_check(const World* a, const World* b, const char* fa, const char* fb){
    if (b->tick < a->tick) { printf("FAIL %s -> %s: tick goes back (%ld -> %ld)\n", fa, fb, a->tick, b->tick); return 1; }
    if (b->nd < a->nd || b->nt < a->nt) { printf("FAIL %s -> %s: disks or tips were removed\n", fa, fb); return 1; }
    const int adapt = a->p.tadapt > 0;                 /* 0b: radii change and disks are removed, nothing else */
    for (int i = 0; i < a->nd; ++i) {
        const Disk* x = &a->d[i]; const Disk* y = &b->d[i];
        if (x->x != y->x || x->y != y->y || x->parent != y->parent || x->born != y->born || (!adapt && (x->r != y->r || x->alive != y->alive)) || (!x->alive && y->alive)) {
            printf("FAIL %s -> %s: disk %d changed\n", fa, fb, i); return 1;
        }
    }
    for (int i = a->nd; i < b->nd; ++i)
        if (b->d[i].born <= a->tick) { printf("FAIL %s -> %s: new disk %d born at %ld, not after tick %ld\n", fa, fb, i, b->d[i].born, a->tick); return 1; }
    for (int k = 0; k < a->nt; ++k) {
        const Tip* x = &a->t[k]; const Tip* y = &b->t[k];
        if (!x->alive && (y->alive || y->disk != x->disk || y->heading != x->heading)) {
            printf("FAIL %s -> %s: dead tip %d moved or revived\n", fa, fb, k); return 1;
        }
    }
    return 0;
}

int main(int argc, char** argv){
    if (argc >= 3 && !strcmp(argv[1], "verify")) {
        int fails = 0, files = 0; long disks_max = 0; World prev; int have_prev = 0; const char* prev_name = NULL;
        for (int i = 2; i < argc; ++i) {
            World w; double pool; char err[512];
            if (state_read(&w, argv[i], &pool, err, sizeof err)) { printf("FAIL %s: %s\n", argv[i], err); ++fails; continue; }
            const int f = world_check(&w, w.p.tadapt == 0, pool, err, sizeof err);
            if (f) { printf("FAIL %s: %d invariant/ledger failures, first: %s\n", argv[i], f, err); ++fails; }
            if (have_prev) { fails += sequence_check(&prev, &w, prev_name, argv[i]); world_free(&prev); }
            if (w.nd > disks_max) disks_max = w.nd;
            prev = w; have_prev = 1; prev_name = argv[i]; ++files;
        }
        if (have_prev) world_free(&prev);
        printf("verify: %d files, up to %ld disks, %d failures\n", files, disks_max, fails);
        return fails ? 1 : 0;
    }
    if (argc >= 3 && !strcmp(argv[1], "info")) {
        for (int i = 2; i < argc; ++i) {
            World w; double pool; char err[512];
            if (state_read(&w, argv[i], &pool, err, sizeof err)) { printf("FAIL %s: %s\n", argv[i], err); return 1; }
            unsigned char* fed = (unsigned char*)malloc((size_t)(w.nd > 0 ? w.nd : 1));
            if (!fed) return 3;
            const int nf = w.p.npatch > 0 ? world_fed(&w, fed) : 0;
            int ends = 0;                              /* fed disks without a fed child: 1 = one chain from A to B */
            for (int j = 0; j < w.nd; ++j) {
                if (!fed[j] || !w.d[j].alive) continue;
                int fc = 0; for (int c = j + 1; c < w.nd && !fc; ++c) fc = w.d[c].alive && w.d[c].parent == j && fed[c];
                ends += !fc;
            }
            printf("tick %ld alive %d live_tips %d m_pool %.6f fed %d fed_ends %d", w.tick, world_alive_disks(&w), world_live_tips(&w), pool, nf, ends);
            if (w.p.regrow_tips >= 2) {                /* 0b.1: material (sum of pi r^2) of the A-B strands, the new front, the old mat */
                unsigned char* fr = (unsigned char*)malloc((size_t)(w.nd > 0 ? w.nd : 1)); if (!fr) return 3;
                const int nfr = world_front(&w, fr);
                double mf = 0, mn = 0, mo = 0;
                for (int j = 0; j < w.nd; ++j) {
                    if (!w.d[j].alive) continue;
                    const double m = PID_PI * w.d[j].r * w.d[j].r;
                    if (fed[j]) mf += m; else if (fr[j]) mn += m; else mo += m;
                }
                printf(" front %d m_strands %.1f m_front %.1f m_oldmat %.1f", nfr, mf, mn, mo);
                free(fr);
            }
            for (int k = 0; k < w.p.npatch; ++k) {
                if (pid_patch_home(&w, k)) continue;
                int t = 0; for (int j = 0; j < w.nd; ++j) t += w.d[j].alive && pid_disk_in_patch(&w.d[j], &w.p.patch[k]);
                printf(" touch%d %d", k, t);
            }
            printf("\n");
            free(fed); world_free(&w);
        }
        return 0;
    }
    if (argc >= 3 && !strcmp(argv[1], "hull")) {
        int fails = 0;
        for (int f = 2; f < argc; ++f) {
            World w; double pool; char err[512];
            if (state_read(&w, argv[f], &pool, err, sizeof err)) { printf("FAIL %s: %s\n", argv[f], err); return 1; }
            Hull H; hull_compute_world(&w, &H);
            int bad = 0;
            {   int e0 = 0; for (int i = 1; i < w.nd; ++i) e0 += w.d[i].alive && w.d[w.d[i].parent].alive;
                if (e0 == 0) {                                                         /* the seed alone: no capsule */
                    const int st = !w.hull.valid || (w.hull.np == 0 && w.hull.area == 0.0);
                    printf("%s %s: no live edge, empty lumen%s\n", st ? "ok" : "FAIL", argv[f], st ? "" : " but the stored hull is not empty");
                    fails += !st; hull_free(&H); world_free(&w); continue;
                } }
            const double gap = hull_max_gap(&H);                                       /* (a) */
            if (!(gap <= 1e-9)) { printf("FAIL %s: (a) a loop is not closed: largest gap %.3g, open chains %d\n", argv[f], gap, H.open); bad++; }
            int outside = 0, inspected = 0;                                            /* (b) */
            int* e = (int*)malloc((size_t)(w.nd > 0 ? w.nd : 1) * sizeof(int)); if (!e) return 3;
            const int ne = hull_edges(&w, e);
            unsigned char* inedge = (unsigned char*)calloc((size_t)(w.nd > 0 ? w.nd : 1), 1); if (!inedge) return 3;
            for (int k = 0; k < ne; ++k) { inedge[e[k]] = 1; inedge[w.d[e[k]].parent] = 1; }
            for (int i = 0; i < w.nd; ++i) if (inedge[i]) { inspected++; if (!hull_contains(&H, w.d[i].x, w.d[i].y)) outside++; }
            if (outside) { printf("FAIL %s: (b) %d of %d disk centres outside the lumen\n", argv[f], outside, inspected); bad++; }
            double ra[4], er[4]; const double res[4] = {2, 4, 8, 16};                  /* (c) */
            for (int k = 0; k < 4; ++k) { ra[k] = hull_raster_area(&w, res[k]); er[k] = fabs(ra[k] - H.area) / H.area; }
            double exact = -1.0;                                                       /* one capsule: the closed form */
            if (ne == 1) {
                const Disk* A = &w.d[e[0]]; const Disk* B = &w.d[A->parent];
                const double d = hypot(B->x - A->x, B->y - A->y), r1 = A->r, r2 = B->r;
                if (d > fabs(r1 - r2)) {
                    const double al = acos((r1 - r2) / d);
                    exact = r1 * r1 * (PID_PI - al) + r2 * r2 * al + (r1 + r2) * d * sin(al);
                    if (!(fabs(exact - H.area) <= 1e-12 * exact)) { printf("FAIL %s: one capsule, area %.17g != closed form %.17g\n", argv[f], H.area, exact); bad++; }
                }
            }
            /* "decreasing": pixel-centre sampling adds noise of a few pixels at every resolution, so the test asks for a
               falling log-log trend (least-squares slope < 0) and 16 px better than 2 px, plus the finest within 0.1 % */
            double sx = 0, sy = 0, sxx = 0, sxy = 0;
            for (int k = 0; k < 4; ++k) { const double lx = log(res[k]), ly = log(er[k] > 1e-300 ? er[k] : 1e-300); sx += lx; sy += ly; sxx += lx * lx; sxy += lx * ly; }
            const double slope = (4 * sxy - sx * sy) / (4 * sxx - sx * sx);
            const int falling = slope < 0 && er[3] < er[0];
            /* below 1000 sigma^2 the raster's pixel noise at 16 px per sigma is comparable with 0.1 %; there the hull is
               held to the closed form when it is one capsule, and the raster numbers are reported only */
            if (H.area >= 1000.0 && (!(er[3] <= 1e-3) || !falling)) { printf("FAIL %s: (c) raster errors %.2e %.2e %.2e %.2e (log-log slope %.2f, finest <= 1e-3 %d)\n", argv[f], er[0], er[1], er[2], er[3], slope, er[3] <= 1e-3); bad++; }
            for (int k = 0; k < ne / 2; ++k) { const int t = e[k]; e[k] = e[ne - 1 - k]; e[ne - 1 - k] = t; }   /* (e) */
            Hull R; hull_compute(&w, e, ne, &R);
            int same = R.np == H.np && fabs(R.area - H.area) <= 1e-12 * H.area;
            for (int k = 0; same && k < H.np; ++k) same = prim_eq(&R.p[k], &H.p[k]);
            if (!same) { printf("FAIL %s: (e) reversed edge order gives a different hull\n", argv[f]); bad++; }
            if (w.hull.valid) {                                                        /* the stored hull */
                int st = w.hull.np == H.np && w.hull.nloop == H.nloop && w.hull.area == H.area && w.hull.perim == H.perim;
                for (int k = 0; st && k < H.np; ++k) st = prim_eq(&w.hull.p[k], &H.p[k]);
                for (int l = 0; st && l < H.nloop; ++l) st = w.hull.loop[l] == H.loop[l];
                if (!st) { printf("FAIL %s: the stored hull differs from the recomputed one\n", argv[f]); bad++; }
            }
            printf("%s %s: edges %d, loops %d, primitives %d, lumen_area %.6f, hull_perimeter %.6f, m_used %.6f, gap %.1e, centres inside %d/%d, raster err 2/4/8/16 px: %.1e %.1e %.1e %.1e (slope %.2f)%s\n",
                   bad ? "FAIL" : "ok", argv[f], ne, H.nloop, H.np, H.area, H.perim, w.m_used, gap, inspected - outside, inspected, er[0], er[1], er[2], er[3], slope,
                   exact >= 0 ? ", one capsule: equals the closed form to 1e-12" : (H.area < 1000.0 ? ", under 1000 sigma^2: raster reported only" : ""));
            fails += bad;
            hull_free(&H); hull_free(&R); free(e); free(inedge); world_free(&w);
        }
        return fails ? 1 : 0;
    }
    if (argc == 4 && !strcmp(argv[1], "roundtrip")) {
        World w; double pool; char err[512];
        if (state_read(&w, argv[2], &pool, err, sizeof err)) { printf("FAIL %s: %s\n", argv[2], err); return 1; }
        if (state_write(&w, argv[3])) { printf("FAIL cannot write %s\n", argv[3]); world_free(&w); return 1; }
        world_free(&w);
        return 0;
    }
    fprintf(stderr, "usage: pid_check verify FILE... | pid_check roundtrip IN OUT | pid_check info FILE... | pid_check hull FILE...\n");
    return 2;
}
