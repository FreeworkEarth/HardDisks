/* ##CHRIS 2026-10-09 (M3, 261012 sec. 4.7.14, acceptance 2): a harness cell as a replay file for the driver
   (00ALLINONE --engine-replay=FILE). Everything is written exactly: doubles as C99 hexadecimal floats (%a), so the driver
   builds the same EDMD_Params and particles, advances to the same targets and applies the same body changes as the harness.
   The harness records its own advance targets and changes while it runs (g3rec_*), so the file holds what was actually done.
   Format (one item per line):
     gen3-state v1 / name <s> / N <n> / boxW, boxH, radius, cell_px, t0 <%a> / divider_count <n>
     divider <d> <x> <th> <mass> <vx> <k> <xeq> / pistonL <has> <x> <vx> <mass> / pistonR <has> <x> <vx> <mass>
     particles <N>, then N lines <x> <y> <vx> <vy>
     targets <n>, then n lines <t_abs>
     changes <m>, then m lines <after target index> <D|PR> <mass> <vx>     (D: divider 0; PR: the right piston)
     expect_hash <16 hex digits> / expect_events <executed events, all kinds> / end */
#ifndef GEN3_STATE_IO_H
#define GEN3_STATE_IO_H
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "../edmd.h"

typedef struct { long after; char what[4]; double mass, vx; } G3Change;
typedef struct { double* t; long nt, capt; G3Change* ch; long nch, capch; } G3Rec;
static G3Rec* g3rec = NULL;                       /* non-NULL while a harness run records */

static inline void g3rec_target(double t){
    if (!g3rec) return;
    if (g3rec->nt == g3rec->capt) { g3rec->capt = g3rec->capt ? 2 * g3rec->capt : 1024; g3rec->t = (double*)realloc(g3rec->t, (size_t)g3rec->capt * sizeof(double)); }
    g3rec->t[g3rec->nt++] = t;
}
static inline void g3rec_change(const char* what, double mass, double vx){
    if (!g3rec) return;
    if (g3rec->nch == g3rec->capch) { g3rec->capch = g3rec->capch ? 2 * g3rec->capch : 64; g3rec->ch = (G3Change*)realloc(g3rec->ch, (size_t)g3rec->capch * sizeof(G3Change)); }
    G3Change* c = &g3rec->ch[g3rec->nch++];
    c->after = g3rec->nt - 1; snprintf(c->what, sizeof c->what, "%s", what); c->mass = mass; c->vx = vx;
}
static inline int g3state_write(const char* path, const char* name, const EDMD_Params* p, const EDMD_Particle* P, double cell_px, double t0,
                         const G3Rec* r, uint64_t hash, long events){
    FILE* f = fopen(path, "w");
    if (!f) return 0;
    fprintf(f, "gen3-state v1\nname %s\nN %d\nboxW %a\nboxH %a\nradius %a\ncell_px %a\nt0 %a\ndivider_count %d\n",
            name, p->N, p->boxW, p->boxH, p->radius, cell_px, t0, p->divider_count);
    for (int d = 0; d < p->divider_count && d < EDMD_MAX_DIVIDERS; ++d)
        fprintf(f, "divider %d %a %a %a %a %a %a\n", d, p->divider_x[d], p->divider_thickness[d], p->divider_mass[d], p->divider_vx[d],
                p->divider_k[d], p->divider_xeq[d]);
    fprintf(f, "pistonL %d %a %a %a\npistonR %d %a %a %a\n", p->has_pistonL, p->pistonL_x, p->pistonL_vx, p->pistonL_mass,
            p->has_pistonR, p->pistonR_x, p->pistonR_vx, p->pistonR_mass);
    fprintf(f, "particles %d\n", p->N);
    for (int i = 0; i < p->N; ++i) fprintf(f, "%a %a %a %a\n", P[i].x, P[i].y, P[i].vx, P[i].vy);
    fprintf(f, "targets %ld\n", r->nt);
    for (long k = 0; k < r->nt; ++k) fprintf(f, "%a\n", r->t[k]);
    fprintf(f, "changes %ld\n", r->nch);
    for (long k = 0; k < r->nch; ++k) fprintf(f, "%ld %s %a %a\n", r->ch[k].after, r->ch[k].what, r->ch[k].mass, r->ch[k].vx);
    fprintf(f, "expect_hash %016llx\nexpect_events %ld\nend\n", (unsigned long long)hash, events);
    return fclose(f) == 0;
}
#endif
