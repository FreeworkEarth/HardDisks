/* ##CHRIS 2026-10-09 (M3, 261012 sec. 4.7.14, amendment a): the band-edge stress cells, shared by gen3_band_edge_test.c (the
   audited runs) and gen3_band_edge_ties.c (the classification of the exact-tie findings). A body's band holds the cell columns
   that meet its contact positions over the band's life, [lo - h, hi + h], by a CLOSED test with the margin m = tol_wall. The
   geometry is dyadic, so the reach ends EXACTLY on a column boundary: cell width w = 32 px, R = 12 px, divider thickness 1 px,
   h = th/2 + R = 12.5 px.
     held_right   a held divider (mass 0, at rest) whose right contact edge x_d + h is a column boundary (288 px)
     held_left    a held divider whose left contact edge x_d - h is a column boundary (256 px)
     driven       a divider of mass 0 at +-0.25 px/unit: each band ends after w/|v| = 128 units exactly at reach x + 32, and
                  x + 32 + h is a column boundary; the velocity flips at every expiry (an API change tied with the BAND event)
     spring       a spring divider (M = 50, period 200 units) started at its turning point; xeq + A + h is a column boundary
   Disks (16 rows of 32 px, box 512 x 512 px):
     rows 0-7, ON the boundary and just OUTSIDE the reach (+0, +1e-12, +1 ulp, +1e-9, +1e-6 px), in the boundary column and in
       two columns further out (40 px apart, the third one moving the other way);
     rows 8-15, just INSIDE the reach, in the boundary column. Where the body sits at its reach at t = 0 (held cells, the
       spring at its turning point) "inside" means touching within rounding (1 ulp .. 2e-12 px of overlap, below tol_face);
       for the driven body, which reaches the boundary only at t = 128, the insets go to 1e-6 px and 0.5 px;
     ties: several rows with the same x and velocity (simultaneous contacts), velocities towards, away and at rest;
     the far side of the divider: 8 disks at rest and moving.
   "exact": every y and velocity dyadic too. "generic y": x exactly as in "exact" (the x ties kept), y and v_y shifted by
   non-dyadic amounts. */
#ifndef GEN3_BAND_EDGE_CELLS_H
#define GEN3_BAND_EDGE_CELLS_H
#include <string.h>
#include "../edmd.h"

#define BE_W 32.0
#define BE_R 12.0
#define BE_TH 1.0
#define BE_HH (0.5 * BE_TH + BE_R)
#define BE_ULP 4.440892098500626e-14      /* rounds to 1 ulp at 256 and 288 px */
#define BE_NCELL 8

typedef struct { const char* name; double xd, mass, v, k, xeq; double edge; int side; double T, flip; int generic; int reach_at_t0; } BECell;

static inline int be_place(const BECell* c, EDMD_Particle* P){
    static const double off[8] = {0.0, 1e-12, BE_ULP, 1e-9, 0.0, 0.0, 1e-6, 0.0};       /* outside the reach */
    static const double vel[8] = {-0.25, -0.25, -0.25, 0.0, -0.25, 0.25, -1e-6, -1.0};   /* towards the body < 0 */
    static const double ins_touch[8] = {1e-12, BE_ULP, 1e-12, BE_ULP, 2e-12, 1e-12, BE_ULP, 2e-12};   /* inside, within rounding */
    static const double ins_free[8]  = {1e-12, BE_ULP, 1e-9, 1e-6, 0.5, 1e-12, 1e-6, BE_ULP};         /* inside, body not there yet */
    static const double vin[8] = {-0.25, -0.25, 0.25, 0.0, -0.25, 0.25, 0.0, -1.0};
    int n = 0;
    for (int col = 0; col < 3; ++col)
        for (int row = 0; row < 8; ++row) {
            EDMD_Particle* p = &P[n++]; memset(p, 0, sizeof *p);
            p->x = c->edge + c->side * (off[row] + 40.0 * col);
            p->y = 16.0 + 32.0 * row + (c->generic ? 0.3141592653589793 * (1 + row % 3) + 0.0271828 * col : 0.0);
            p->vx = c->side * vel[row] * (col == 2 ? -1.0 : 1.0);
            p->vy = ((row % 3 == 0) ? 0.0 : 0.0625 * (row - 4)) + (c->generic ? 0.001414213562 * (row + 1) + 0.000577 * col : 0.0);
        }
    for (int r = 0; r < 8; ++r) {
        EDMD_Particle* p = &P[n++]; memset(p, 0, sizeof *p);
        const int row = 8 + r;
        p->x = c->edge - c->side * (c->reach_at_t0 ? ins_touch[r] : ins_free[r]);
        p->y = 16.0 + 32.0 * row + (c->generic ? 0.2718281828459045 * (1 + r % 3) : 0.0);
        p->vx = c->side * vin[r];
        p->vy = ((r % 3 == 0) ? 0.0 : 0.0625 * (r - 4)) + (c->generic ? 0.001732050808 * (r + 1) : 0.0);
    }
    for (int row = 0; row < 16; row += 2) {
        EDMD_Particle* p = &P[n++]; memset(p, 0, sizeof *p);
        p->x = c->xd - c->side * (BE_HH + 20.0 + 3.0 * (row % 8));
        p->y = 16.0 + 32.0 * row + (c->generic ? 0.1732050808 * (row + 1) : 0.0);
        p->vx = -c->side * 0.125 * ((row % 8) - 3); p->vy = c->generic ? 0.000707 * (row + 2) : 0.0;
    }
    return n;
}

static inline void be_params(const BECell* c, int N, EDMD_Params* p){
    memset(p, 0, sizeof *p);
    p->boxW = 16.0 * BE_W; p->boxH = 16.0 * BE_W; p->radius = BE_R; p->N = N; p->pp_collisions_enabled = 1; p->particle_mass = 1.0; p->kB = 1.0;
    p->divider_count = 1; p->divider_x[0] = c->xd; p->divider_thickness[0] = BE_TH; p->divider_mass[0] = c->mass; p->divider_vx[0] = c->v;
    p->divider_k[0] = c->k; p->divider_xeq[0] = c->xeq;
}

static inline void be_cells(BECell cells[BE_NCELL]){
    const double w2 = 2.0 * 3.14159265358979323846 / 200.0, M = 50.0, A = 16.0;
    const BECell base[4] = {
        {"held_right", 288.0 - BE_HH, 0.0, 0.0, 0.0, 0.0, 288.0, +1, 480.0, 0.0, 0, 1},
        {"held_left",  256.0 + BE_HH, 0.0, 0.0, 0.0, 0.0, 256.0, -1, 480.0, 0.0, 0, 1},
        {"driven",     288.0 - BE_HH - BE_W, 0.0, 0.25, 0.0, 0.0, 288.0, +1, 640.0, 128.0, 0, 0},
        {"spring",     288.0 - BE_HH, M, 0.0, M * w2 * w2, 288.0 - BE_HH - A, 288.0, +1, 800.0, 0.0, 0, 1}};
    for (int k = 0; k < BE_NCELL; ++k) { cells[k] = base[k % 4]; cells[k].generic = k >= 4; }
}
#endif
