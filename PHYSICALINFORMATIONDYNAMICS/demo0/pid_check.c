/* pid_check.c -- the test tool of brief sec. 6; reads state files only, so every number it checks comes from the headless
   path's files.
     pid_check verify FILE...      invariants and ledger on every file (items 2-3); with several files, given in tick order,
                                   also that disks never change and dead tips stay frozen from one state to the next
     pid_check roundtrip IN OUT    read IN, write OUT (item 4; the test then compares the two with cmp)
     pid_check info FILE...        0b: one line per state: tick, alive disks, live tips, M_pool, fed disks, fed path ends, and for every
                                   non-home patch the number of alive disks touching it */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "pid.h"

const char* pid_build_id(void){ return "pid_check"; }   /* world_init() is not used here; state files carry their own build line */

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
    if (argc == 4 && !strcmp(argv[1], "roundtrip")) {
        World w; double pool; char err[512];
        if (state_read(&w, argv[2], &pool, err, sizeof err)) { printf("FAIL %s: %s\n", argv[2], err); return 1; }
        if (state_write(&w, argv[3])) { printf("FAIL cannot write %s\n", argv[3]); world_free(&w); return 1; }
        world_free(&w);
        return 0;
    }
    fprintf(stderr, "usage: pid_check verify FILE... | pid_check roundtrip IN OUT | pid_check info FILE...\n");
    return 2;
}
