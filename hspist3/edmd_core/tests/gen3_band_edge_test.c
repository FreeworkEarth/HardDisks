/* ##CHRIS 2026-10-09 (M3, 261012 sec. 4.7.14, amendment a): band-edge stress cells (the cells: gen3_band_edge_cells.h).
   Every cell runs with the schedule audit after EVERY event, the body audit after every BAND event and API change, and the
   contact audit. Printed per cell: events, audited states, missing / extra / deferred early (all classes), bands missing /
   extra / short, contacts within rounding, the health line and clean.
   Each cell runs twice: "exact" (y positions and velocities dyadic too) and "generic y" (x exactly the same, the x ties kept,
   y and v_y shifted by non-dyadic amounts). The exact variant also makes pairs graze EXACTLY (vertical offset 2R, no
   vertical relative velocity: the discriminant is 0 in exact arithmetic, the impulse of such a touch is 0), gives pairs
   velocities equal in exact arithmetic (a relative velocity at rounding level, a contact predicted ~1e14 units ahead or not
   at all), and lets disks reach cell CORNERS exactly (an x and a y crossing at one time). Rounding then decides, differently in
   the engine's local coordinates and in the audit's absolute ones, whether such a pair meets and which crossing comes first.
   gen3_band_edge_ties.c prints the numbers behind each such finding. The band classes and clean must hold in all 8 cells;
   the generic-y cells must have no finding in any class.
   build (from hspist3/): cc -std=c11 -O2 -ffp-contract=off -Wall -Wextra -o <out> edmd_core/tests/gen3_band_edge_test.c edmd_core/edmd_gen3.c -lm */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../edmd.h"
#include "../edmd_gen3.h"
#include "gen3_band_edge_cells.h"

static int run(const BECell* c){
    EDMD_Particle P[64]; const int N = be_place(c, P);
    EDMD_Params p; be_params(c, N, &p);
    char err[256];
    EDMD3* S = edmd3_create(&p, BE_W, err, sizeof err);
    if (!S) { printf("%s: create: %s\n", c->name, err); return 3; }
    if (!edmd3_load(S, P, 0.0, err, sizeof err)) { printf("%s: load: %s\n", c->name, err); edmd3_destroy(S); return 3; }
    edmd3_set_contact_audit(S, 1); edmd3_set_schedule_audit(S, 1); edmd3_set_schedule_audit_bodies(S, 1);
    double v = c->v;
    for (double t = 1.0; t <= c->T + 1e-9; t += 1.0) {
        edmd3_advance_to(S, t);
        if (c->flip > 0.0 && fmod(t, c->flip) == 0.0) { v = -v; edmd3_set_divider_motion(S, 0, 0.0, v); }   /* at the band expiry */
        const char* msg; if (edmd3_fatal(S, &msg)) { printf("%s: fatal: %s\n", c->name, msg); break; }
    }
    edmd3_schedule_audit_now(S);
    const EDMD3_Health* H = edmd3_health(S); const EDMD3_Audit* A = edmd3_schedule_audit_stats(S);
    double g[4]; edmd3_contact_audit_stats4(S, g);
    const long miss = A->pair_missing + A->wall_missing + A->cross_missing + A->div_missing + A->pis_missing;
    const long extra = A->pair_extra + A->wall_extra + A->cross_extra + A->div_extra + A->pis_extra;
    const long early = A->pair_deferred_early + A->wall_deferred_early + A->div_deferred_early + A->pis_deferred_early;
    const long drel = A->pair_dt_rel + A->wall_dt_rel + A->cross_dt_rel + A->div_dt_rel + A->pis_dt_rel + A->cross_dup;
    printf("| %s (%s) | %d | %.0f | %ld / %ld / %ld / %ld / %ld | %ld | %ld | %ld | %ld | %ld | %ld / %ld / %ld | %ld | %.2g | %d |\n",
           c->name, c->generic ? "generic y" : "exact", N, c->T, H->ev_pair, H->ev_wall, H->ev_div, H->ev_band, H->ev_cross, A->audits,
           miss, extra, early, drel, A->band_missing, A->band_extra, A->band_short, H->obj_contact_now + H->contact_now, g[2], edmd3_health_clean(S));
    printf("|  | health: overlap_repair %ld, obj_overlap_repair %ld, wall_overdue %ld, past_event %ld, clamp %ld, cell %ld, grid %ld, "
           "local/full findings %ld/%ld, body findings %ld |  |  |  |  |  |  |  |  |  |  |  |  |\n",
           H->overlap_repair, H->obj_overlap_repair, H->wall_overdue, H->past_event, H->clamp_repair, H->cell_repair, H->grid_escape,
           H->local_findings, H->full_findings, H->body_findings);
    const int bad = ((miss || extra || early || drel) ? 1 : 0) | ((A->band_missing || A->band_extra || A->band_short || !edmd3_health_clean(S)) ? 2 : 0);
    edmd3_destroy(S);
    return bad;
}

int main(void){
    BECell cells[BE_NCELL]; be_cells(cells);
    printf("# gen3 band-edge stress cells (amendment a): dyadic geometry, the band reach ending exactly on a column boundary\n\n"
           "w = %g px, R = %g px, divider thickness %g px, h = %g px; disks on the boundary, just outside and just inside the reach; "
           "schedule audit after every event, body audit after every BAND event and API change, contact audit\n\n", BE_W, BE_R, BE_TH, BE_HH);
    printf("| cell | N | T [units] | events pair / wall / divider / band / crossing | audited states | missing (all classes) | "
           "extra | deferred early | time disagreements (rel > 1e-10) and duplicates | bands missing / extra / short | contacts within rounding | "
           "max divider contact gap [px] | clean |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|\n");
    int bad_band = 0, bad_generic = 0;
    for (int k = 0; k < BE_NCELL; ++k) { const int b = run(&cells[k]); if (b & 2) bad_band++; if (cells[k].generic && b) bad_generic++; }
    printf("\ncells with a band finding (bands missing, extra or short), clean = 0 or a load failure: %d of %d; generic-y cells with any "
           "finding in any class: %d of %d\n", bad_band, BE_NCELL, bad_generic, BE_NCELL / 2);
    return (bad_band || bad_generic) ? 1 : 0;
}
