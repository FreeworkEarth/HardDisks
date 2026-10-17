/* pid_draw1.c -- the drawing of Demo 1's windows (pid_draw1.h). Ideas (scanline disk fill with alpha, the hull as polylines, a
   streaming texture for a field, Demo 0's colours of a growing network) are taken from demo0/pid_view.c, which is not linked.
   Nothing here changes the World, the run or the engine state: the particles are read through edmd3_particles (computed copies). */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "pid_draw1.h"

void pd_to_screen(const Pid1Vp* v, double x, double y, float* sx, float* sy){
    *sx = (float)(v->ox + x * v->scale);
    *sy = (float)(v->oy + (v->H - y) * v->scale);
}
void pd_fill_disk(SDL_Renderer* R, float cx, float cy, float rad){
    const int y0 = (int)floorf(cy - rad), y1 = (int)ceilf(cy + rad);
    for (int y = y0; y <= y1; ++y) {
        const float dy = ((float)y + 0.5f) - cy, h = rad * rad - dy * dy;
        if (h <= 0.0f) continue;
        const float hw = sqrtf(h);
        const SDL_FRect row = {cx - hw, (float)y, 2.0f * hw, 1.0f};
        SDL_RenderFillRectF(R, &row);
    }
}
void pd_ring(SDL_Renderer* R, float cx, float cy, float rad){
    for (int k = 0; k < 96; ++k) {
        const float a0 = 6.2831853f * (float)k / 96.0f, a1 = 6.2831853f * (float)(k + 1) / 96.0f;
        SDL_RenderDrawLineF(R, cx + rad * cosf(a0), cy + rad * sinf(a0), cx + rad * cosf(a1), cy + rad * sinf(a1));
    }
}
/* a five-stop ramp (dark blue, blue, teal, yellow-green, yellow) of s in [0, 1] */
void pd_ramp(double s, Uint8* r, Uint8* g, Uint8* b){
    static const double C[5][3] = {{20, 30, 90}, {40, 90, 200}, {30, 170, 160}, {170, 220, 70}, {255, 235, 60}};
    s = s < 0.0 ? 0.0 : s > 1.0 ? 1.0 : s;
    const double u = s * 4.0; int k = (int)u; if (k > 3) k = 3;
    const double f = u - k;
    *r = (Uint8)(C[k][0] + f * (C[k + 1][0] - C[k][0])); *g = (Uint8)(C[k][1] + f * (C[k + 1][1] - C[k][1])); *b = (Uint8)(C[k][2] + f * (C[k + 1][2] - C[k][2]));
}
/* the particles' speed scale: bright blue (at rest), white (the most probable speed region), red (2.5 sqrt(kT) and above),
   every colour visible on the dark background */
void pd_speed_colour(double s, Uint8* r, Uint8* g, Uint8* b){
    s = s < 0.0 ? 0.0 : s > 1.0 ? 1.0 : s;
    if (s < 0.5) { const double f = s / 0.5; *r = (Uint8)(70 + f * 165); *g = (Uint8)(130 + f * 105); *b = 255; }
    else { const double f = (s - 0.5) / 0.5; *r = 255; *g = (Uint8)(235 - f * 170); *b = (Uint8)(255 - f * 205); }
}

/* the coarse grid (the density files' 2 sigma cells) as a texture; the top colour is 3 x the mean over the non-empty cells of
   the frame (clipped above), so that a field near equilibrium shows its spread and a pulse stands out (not the largest cell,
   which a single fluctuation sets) */
static void overlay(SDL_Renderer* R, const Pid1Vp* v, Run1* r, int which, SDL_Texture** tex){
    if (!*tex) {
        *tex = SDL_CreateTexture(R, SDL_PIXELFORMAT_ARGB8888, SDL_TEXTUREACCESS_STREAMING, r->gx, r->gy);
        if (!*tex) return;
        SDL_SetTextureBlendMode(*tex, SDL_BLENDMODE_BLEND);
    }
    run1_bin(r, edmd3_particles(r->S));
    const long nc = (long)r->gx * r->gy;
    double vsum = 0.0; long nz = 0;
    for (long c = 0; c < nc; ++c) {
        const double x = which == 1 ? (double)r->cnt[c] : which == 2 ? r->kec[c] : r->mat[c];
        if (x > 0.0) { vsum += x; ++nz; }
    }
    const double vmax = nz ? 3.0 * vsum / nz : 0.0;
    void* px; int pitch;
    if (SDL_LockTexture(*tex, NULL, &px, &pitch)) return;
    for (int iy = 0; iy < r->gy; ++iy) {
        Uint32* row = (Uint32*)((char*)px + (size_t)(r->gy - 1 - iy) * (size_t)pitch);   /* world y up */
        for (int ix = 0; ix < r->gx; ++ix) {
            const long c = (long)iy * r->gx + ix;
            const double x = which == 1 ? (double)r->cnt[c] : which == 2 ? r->kec[c] : r->mat[c];
            if (!(x > 0.0) || !(vmax > 0.0)) { row[ix] = 0; continue; }
            Uint8 cr, cg, cb; pd_ramp(x / vmax, &cr, &cg, &cb);
            row[ix] = (210u << 24) | ((Uint32)cr << 16) | ((Uint32)cg << 8) | cb;
        }
    }
    SDL_UnlockTexture(*tex);
    float x0, y0, x1, y1; pd_to_screen(v, 0.0, 0.0, &x0, &y0); pd_to_screen(v, r->gx * PID1_GRID, r->gy * PID1_GRID, &x1, &y1);
    const SDL_FRect dst = {x0, y1, x1 - x0, y0 - y1};
    SDL_RenderCopyF(R, *tex, NULL, &dst);
}

static void hull(SDL_Renderer* R, const Pid1Vp* v, const Hull* H){
    SDL_SetRenderDrawColor(R, 120, 230, 255, 255);
    for (int k = 0; k < H->np; ++k) {
        const HPrim* q = &H->p[k];
        if (q->type == 0) {
            float ax, ay, bx, by; pd_to_screen(v, q->v[0], q->v[1], &ax, &ay); pd_to_screen(v, q->v[2], q->v[3], &bx, &by);
            SDL_RenderDrawLineF(R, ax, ay, bx, by);
        } else {
            const int n = 1 + (int)((q->v[4] - q->v[3]) / 0.05);
            float px = 0, py = 0;
            for (int j = 0; j <= n; ++j) {
                const double a = q->v[3] + (q->v[4] - q->v[3]) * j / n;
                float sx, sy; pd_to_screen(v, q->v[0] + q->v[2] * cos(a), q->v[1] + q->v[2] * sin(a), &sx, &sy);
                if (j) SDL_RenderDrawLineF(R, px, py, sx, sy);
                px = sx; py = sy;
            }
        }
    }
}

static unsigned char* g_buf = NULL; static int g_cap = 0;   /* fed, live tip and front flags per disk (scratch) */

void pd_scene(SDL_Renderer* R, const Pid1Vp* v, const World* w, Run1* run, const Pid1Scene* sc, SDL_Texture** tex){
    float x0, y0, x1, y1; pd_to_screen(v, 0.0, 0.0, &x0, &y0); pd_to_screen(v, w->p.W, w->p.H, &x1, &y1);
    const SDL_FRect box = {x0, y1, x1 - x0, y0 - y1};
    SDL_SetRenderDrawColor(R, 34, 38, 48, 255); SDL_RenderFillRectF(R, &box);
    SDL_SetRenderDrawColor(R, 120, 124, 140, 255); SDL_RenderDrawRectF(R, &box);
    if (sc->overlay && run) overlay(R, v, run, sc->overlay, tex);
    if ((sc->show & PID1_SHOW_DISKS) && sc->alpha > 0 && w->nd > 0) {
        if (3 * w->nd > g_cap) { g_cap = 6 * w->nd; g_buf = (unsigned char*)realloc(g_buf, (size_t)g_cap); if (!g_buf) { g_cap = 0; return; } }
        unsigned char *fed = g_buf, *live = g_buf + w->nd, *front = g_buf + 2 * w->nd;
        if (w->p.npatch > 0) world_fed(w, fed); else memset(fed, 0, (size_t)w->nd);
        if (sc->growth) {
            memset(live, 0, (size_t)w->nd);
            for (int k = 0; k < w->nt; ++k) if (w->t[k].alive) live[w->t[k].disk] = 1;
            world_front(w, front);
        }
        const Uint8 a = (Uint8)sc->alpha;
        for (int i = 0; i < w->nd; ++i) {
            const Disk* d = &w->d[i]; if (!d->alive) continue;
            if (d->kind == KIND_SEED) SDL_SetRenderDrawColor(R, 225, 95, 80, a);
            else if (fed[i])          SDL_SetRenderDrawColor(R, 250, 235, 90, a);                   /* on a fed path (Demo 0b) */
            else if (sc->growth && front[i] && !live[i]) SDL_SetRenderDrawColor(R, 70, 210, 170, a);   /* the new front (0b.1) */
            else if (sc->growth && live[i])             SDL_SetRenderDrawColor(R, 250, 150, 60, a);   /* a live tip */
            else if (sc->growth && d->kind == KIND_TIP) SDL_SetRenderDrawColor(R, 150, 152, 165, a);  /* a dead end */
            else                      SDL_SetRenderDrawColor(R, 80, 150, 235, a);
            float sx, sy; pd_to_screen(v, d->x, d->y, &sx, &sy);
            pd_fill_disk(R, sx, sy, (float)(d->r * v->scale));
        }
    }
    if (sc->show & PID1_SHOW_FOOD)
        for (int k = 0; k < w->p.npatch; ++k) {
            const Patch* q = &w->p.patch[k]; float sx, sy; pd_to_screen(v, q->x, q->y, &sx, &sy);
            if (pid_patch_home(w, k)) SDL_SetRenderDrawColor(R, 120, 240, 160, 255); else SDL_SetRenderDrawColor(R, 255, 120, 230, 255);
            pd_ring(R, sx, sy, (float)(q->r * v->scale));
        }
    if (sc->show & PID1_SHOW_GRAPH) {
        SDL_SetRenderDrawColor(R, 240, 240, 240, 200);
        for (int i = 1; i < w->nd; ++i) {
            const Disk* d = &w->d[i]; if (!d->alive || d->parent < 0 || !w->d[d->parent].alive) continue;
            const Disk* q = &w->d[d->parent];
            float ax, ay, bx, by; pd_to_screen(v, d->x, d->y, &ax, &ay); pd_to_screen(v, q->x, q->y, &bx, &by);
            SDL_RenderDrawLineF(R, ax, ay, bx, by);
        }
    }
    if ((sc->show & PID1_SHOW_PARTICLES) && run) {
        const EDMD_Particle* P = edmd3_particles(run->S);
        const float rad = (float)(run->R / PID1_PXS * v->scale);
        const double vs = 2.5 * sqrt(run->kbt);                    /* the speed of the top colour [px per unit] */
        for (int i = 0; i < run->N; ++i) {
            Uint8 cr, cg, cb; pd_speed_colour(hypot(P[i].vx, P[i].vy) / vs, &cr, &cg, &cb);
            SDL_SetRenderDrawColor(R, cr, cg, cb, 255);
            float sx, sy; pd_to_screen(v, P[i].x / PID1_PXS, P[i].y / PID1_PXS, &sx, &sy);
            pd_fill_disk(R, sx, sy, rad);
        }
    }
    if ((sc->show & PID1_SHOW_HULL) && w->hull.valid) hull(R, v, &w->hull);
}
