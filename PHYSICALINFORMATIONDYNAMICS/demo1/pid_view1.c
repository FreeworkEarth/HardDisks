/* pid_view1.c -- the window of Demo 1 (02_DEMO1_BRIEF.md sec. 9): the frozen Demo 0 network and the hard disks running in
   its lumen. It draws and never owns state: the run moves only through run1_advance (pid_demo1.c), the function the headless
   path uses, so a window run writes the same pstate_*.txt, density_*.txt and measure.csv (make test, test 11). Drawing ideas
   (scanline disk fill with alpha, the hull as polylines, a streaming texture for a field) are taken from demo0/pid_view.c,
   which is not linked. Layers, bottom to top: the coarse-grid overlay (C), the structural disks (D, fill alpha [ ]), the
   food patches (F), the graph (G), the particles (P, coloured by speed), the hull (H). */
#define SDL_MAIN_HANDLED
#include <SDL.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pid_demo1.h"

typedef struct { int show, overlay, alpha, running; double scale; int margin, win_w, win_h; } V1;

static void to_screen(const V1* v, const World* w, double x, double y, float* sx, float* sy){   /* sigma, y up */
    *sx = (float)(v->margin + x * v->scale);
    *sy = (float)(v->margin + (w->p.H - y) * v->scale);
}
static void fill_disk(SDL_Renderer* R, float cx, float cy, float rad){
    const int y0 = (int)floorf(cy - rad), y1 = (int)ceilf(cy + rad);
    for (int y = y0; y <= y1; ++y) {
        const float dy = ((float)y + 0.5f) - cy, h = rad * rad - dy * dy;
        if (h <= 0.0f) continue;
        const float hw = sqrtf(h);
        const SDL_FRect row = {cx - hw, (float)y, 2.0f * hw, 1.0f};
        SDL_RenderFillRectF(R, &row);
    }
}
static void ring(SDL_Renderer* R, float cx, float cy, float rad){
    for (int k = 0; k < 96; ++k) {
        const float a0 = 6.2831853f * (float)k / 96.0f, a1 = 6.2831853f * (float)(k + 1) / 96.0f;
        SDL_RenderDrawLineF(R, cx + rad * cosf(a0), cy + rad * sinf(a0), cx + rad * cosf(a1), cy + rad * sinf(a1));
    }
}
/* a five-stop ramp (dark blue, blue, teal, yellow-green, yellow) of s in [0, 1] */
static void ramp(double s, Uint8* r, Uint8* g, Uint8* b){
    static const double C[5][3] = {{20, 30, 90}, {40, 90, 200}, {30, 170, 160}, {170, 220, 70}, {255, 235, 60}};
    s = s < 0.0 ? 0.0 : s > 1.0 ? 1.0 : s;
    const double u = s * 4.0; int k = (int)u; if (k > 3) k = 3;
    const double f = u - k;
    *r = (Uint8)(C[k][0] + f * (C[k + 1][0] - C[k][0])); *g = (Uint8)(C[k][1] + f * (C[k + 1][1] - C[k][1])); *b = (Uint8)(C[k][2] + f * (C[k + 1][2] - C[k][2]));
}

/* the particles' speed scale: bright blue (at rest), white (the most probable speed region), red (2.5 sqrt(kT) and above),
   every colour visible on the dark background */
static void speed_colour(double s, Uint8* r, Uint8* g, Uint8* b){
    s = s < 0.0 ? 0.0 : s > 1.0 ? 1.0 : s;
    if (s < 0.5) { const double f = s / 0.5; *r = (Uint8)(70 + f * 165); *g = (Uint8)(130 + f * 105); *b = 255; }
    else { const double f = (s - 0.5) / 0.5; *r = 255; *g = (Uint8)(235 - f * 170); *b = (Uint8)(255 - f * 205); }
}

/* the coarse grid (the density files' 2 sigma cells) as a texture; the top colour is 3 x the mean over the non-empty cells of
   the frame (clipped above), so that a field near equilibrium shows its spread and a pulse stands out (not the largest cell,
   which a single fluctuation sets) */
static void draw_overlay(SDL_Renderer* R, const V1* v, Run1* r, SDL_Texture** tex){
    if (!*tex) {
        *tex = SDL_CreateTexture(R, SDL_PIXELFORMAT_ARGB8888, SDL_TEXTUREACCESS_STREAMING, r->gx, r->gy);
        if (!*tex) return;
        SDL_SetTextureBlendMode(*tex, SDL_BLENDMODE_BLEND);
    }
    run1_bin(r, edmd3_particles(r->S));
    const long nc = (long)r->gx * r->gy;
    double vsum = 0.0; long nz = 0;
    for (long c = 0; c < nc; ++c) {
        const double x = v->overlay == 1 ? (double)r->cnt[c] : v->overlay == 2 ? r->kec[c] : r->mat[c];
        if (x > 0.0) { vsum += x; ++nz; }
    }
    const double vmax = nz ? 3.0 * vsum / nz : 0.0;
    void* px; int pitch;
    if (SDL_LockTexture(*tex, NULL, &px, &pitch)) return;
    for (int iy = 0; iy < r->gy; ++iy) {
        Uint32* row = (Uint32*)((char*)px + (size_t)(r->gy - 1 - iy) * (size_t)pitch);   /* world y up */
        for (int ix = 0; ix < r->gx; ++ix) {
            const long c = (long)iy * r->gx + ix;
            const double x = v->overlay == 1 ? (double)r->cnt[c] : v->overlay == 2 ? r->kec[c] : r->mat[c];
            if (!(x > 0.0) || !(vmax > 0.0)) { row[ix] = 0; continue; }
            Uint8 cr, cg, cb; ramp(x / vmax, &cr, &cg, &cb);
            row[ix] = (210u << 24) | ((Uint32)cr << 16) | ((Uint32)cg << 8) | cb;
        }
    }
    SDL_UnlockTexture(*tex);
    float x0, y0, x1, y1; to_screen(v, &r->w, 0.0, 0.0, &x0, &y0); to_screen(v, &r->w, r->gx * PID1_GRID, r->gy * PID1_GRID, &x1, &y1);
    const SDL_FRect dst = {x0, y1, x1 - x0, y0 - y1};
    SDL_RenderCopyF(R, *tex, NULL, &dst);
}

static void draw_hull(SDL_Renderer* R, const V1* v, const World* w){
    const Hull* H = &w->hull;
    SDL_SetRenderDrawColor(R, 120, 230, 255, 255);
    for (int k = 0; k < H->np; ++k) {
        const HPrim* q = &H->p[k];
        if (q->type == 0) {
            float ax, ay, bx, by; to_screen(v, w, q->v[0], q->v[1], &ax, &ay); to_screen(v, w, q->v[2], q->v[3], &bx, &by);
            SDL_RenderDrawLineF(R, ax, ay, bx, by);
        } else {
            const int n = 1 + (int)((q->v[4] - q->v[3]) / 0.05);
            float px = 0, py = 0;
            for (int j = 0; j <= n; ++j) {
                const double a = q->v[3] + (q->v[4] - q->v[3]) * j / n;
                float sx, sy; to_screen(v, w, q->v[0] + q->v[2] * cos(a), q->v[1] + q->v[2] * sin(a), &sx, &sy);
                if (j) SDL_RenderDrawLineF(R, px, py, sx, sy);
                px = sx; py = sy;
            }
        }
    }
}

static void draw(SDL_Renderer* R, const V1* v, Run1* r, unsigned char* fed, SDL_Texture** tex){
    const World* w = &r->w;
    SDL_SetRenderDrawColor(R, 14, 16, 22, 255);
    SDL_RenderClear(R);
    float x0, y0, x1, y1; to_screen(v, w, 0.0, 0.0, &x0, &y0); to_screen(v, w, w->p.W, w->p.H, &x1, &y1);
    const SDL_FRect box = {x0, y1, x1 - x0, y0 - y1};
    SDL_SetRenderDrawColor(R, 34, 38, 48, 255); SDL_RenderFillRectF(R, &box);
    SDL_SetRenderDrawColor(R, 120, 124, 140, 255); SDL_RenderDrawRectF(R, &box);
    if (v->overlay) draw_overlay(R, v, r, tex);
    if ((v->show & PID1_SHOW_DISKS) && v->alpha > 0) {
        if (w->p.npatch > 0) world_fed(w, fed); else memset(fed, 0, (size_t)w->nd);
        for (int i = 0; i < w->nd; ++i) {
            const Disk* d = &w->d[i]; if (!d->alive) continue;
            if (d->kind == KIND_SEED) SDL_SetRenderDrawColor(R, 225, 95, 80, (Uint8)v->alpha);
            else if (fed[i])          SDL_SetRenderDrawColor(R, 250, 235, 90, (Uint8)v->alpha);   /* on a fed path (Demo 0b) */
            else                      SDL_SetRenderDrawColor(R, 80, 150, 235, (Uint8)v->alpha);
            float sx, sy; to_screen(v, w, d->x, d->y, &sx, &sy);
            fill_disk(R, sx, sy, (float)(d->r * v->scale));
        }
    }
    if (v->show & PID1_SHOW_FOOD)
        for (int k = 0; k < w->p.npatch; ++k) {
            const Patch* q = &w->p.patch[k]; float sx, sy; to_screen(v, w, q->x, q->y, &sx, &sy);
            if (pid_patch_home(w, k)) SDL_SetRenderDrawColor(R, 120, 240, 160, 255); else SDL_SetRenderDrawColor(R, 255, 120, 230, 255);
            ring(R, sx, sy, (float)(q->r * v->scale));
        }
    if (v->show & PID1_SHOW_GRAPH) {
        SDL_SetRenderDrawColor(R, 240, 240, 240, 200);
        for (int i = 1; i < w->nd; ++i) {
            const Disk* d = &w->d[i]; if (!d->alive || d->parent < 0 || !w->d[d->parent].alive) continue;
            const Disk* q = &w->d[d->parent];
            float ax, ay, bx, by; to_screen(v, w, d->x, d->y, &ax, &ay); to_screen(v, w, q->x, q->y, &bx, &by);
            SDL_RenderDrawLineF(R, ax, ay, bx, by);
        }
    }
    if (v->show & PID1_SHOW_PARTICLES) {
        const EDMD_Particle* P = edmd3_particles(r->S);
        const float rad = (float)(r->R / PID1_PXS * v->scale);
        const double vs = 2.5 * sqrt(r->kbt);                      /* the speed of the top colour [px per unit] */
        for (int i = 0; i < r->N; ++i) {
            Uint8 cr, cg, cb; speed_colour(hypot(P[i].vx, P[i].vy) / vs, &cr, &cg, &cb);
            SDL_SetRenderDrawColor(R, cr, cg, cb, 255);
            float sx, sy; to_screen(v, w, P[i].x / PID1_PXS, P[i].y / PID1_PXS, &sx, &sy);
            fill_disk(R, sx, sy, rad);
        }
    }
    if (v->show & PID1_SHOW_HULL) draw_hull(R, v, w);
}

static void save_shot(SDL_Renderer* R, Run1* r, double t){
    int w = 0, h = 0; SDL_GetRendererOutputSize(R, &w, &h);
    SDL_Surface* s = SDL_CreateRGBSurfaceWithFormat(0, w, h, 24, SDL_PIXELFORMAT_RGB24);   /* a 24-bit BI_RGB BMP (tests/bmp2png.py) */
    if (!s) { fprintf(stderr, "pid_view1: shot: %s\n", SDL_GetError()); return; }
    char name[96], path[1200]; int k = 0;
    snprintf(name, sizeof name, "shot_%06.0f.bmp", t);
    for (;;) {
        snprintf(path, sizeof path, "%s/%s", r->out, name);
        FILE* f = fopen(path, "r"); if (!f) break;
        fclose(f); snprintf(name, sizeof name, "shot_%06.0f_%d.bmp", t, ++k);
    }
    if (SDL_RenderReadPixels(R, NULL, SDL_PIXELFORMAT_RGB24, s->pixels, s->pitch) || SDL_SaveBMP(s, path))
        fprintf(stderr, "pid_view1: shot: %s\n", SDL_GetError());
    else { run1_note_file(name); printf("shot %s at t %.3f sigma-time (%dx%d)\n", path, t, w, h); }
    SDL_FreeSurface(s);
}

int view1_run(Run1* r){
    SDL_SetMainReady();
    if (SDL_Init(SDL_INIT_VIDEO) != 0) { fprintf(stderr, "pid_view1: SDL_Init: %s\n", SDL_GetError()); return 2; }
    V1 v; memset(&v, 0, sizeof v);
    v.margin = 24; v.running = 1; v.overlay = r->v.overlay;
    v.alpha = (int)lround(r->v.alpha * 255.0);
    v.show = PID1_SHOW_DISKS | PID1_SHOW_PARTICLES | PID1_SHOW_FOOD;
    if (r->v.view == 1) { v.alpha = 0; v.show = PID1_SHOW_HULL | PID1_SHOW_GRAPH | PID1_SHOW_PARTICLES | PID1_SHOW_FOOD; }   /* --view=network */
    if (r->v.view == 2) v.show |= PID1_SHOW_HULL;                                                                            /* --view=both */
    if (r->v.show_given) v.show = r->v.show;
    const World* w = &r->w;
    v.scale = fmin((1400.0 - 2.0 * v.margin) / w->p.W, (860.0 - 2.0 * v.margin) / w->p.H);
    v.win_w = (int)ceil(w->p.W * v.scale) + 2 * v.margin;
    v.win_h = (int)ceil(w->p.H * v.scale) + 2 * v.margin;
    SDL_Window* win = SDL_CreateWindow("PID Demo 1 -- hard disks in the frozen lumen", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, v.win_w, v.win_h, SDL_WINDOW_SHOWN);
    if (!win) { fprintf(stderr, "pid_view1: SDL_CreateWindow: %s\n", SDL_GetError()); SDL_Quit(); return 2; }
    int vsync = 1;
    SDL_Renderer* R = SDL_CreateRenderer(win, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    if (!R) { vsync = 0; R = SDL_CreateRenderer(win, -1, SDL_RENDERER_SOFTWARE); }
    if (!R) { fprintf(stderr, "pid_view1: SDL_CreateRenderer: %s\n", SDL_GetError()); SDL_DestroyWindow(win); SDL_Quit(); return 2; }
    SDL_SetRenderDrawBlendMode(R, SDL_BLENDMODE_BLEND);
    unsigned char* fed = (unsigned char*)calloc((size_t)w->nd + 1, 1);
    SDL_Texture* tex = NULL;
    int nshot = r->v.nshot; double shots[16]; memcpy(shots, r->v.shot_at, sizeof shots);
    int quit = 0, announced = 0;
    static const char* OV[4] = {"off", "particle density", "particle kinetic energy", "structural material"};
    printf("window %dx%d, %.2f px per sigma; Space run/pause, Esc quit, S shot, D H G P F layers, [ ] alpha, C overlay\n", v.win_w, v.win_h, v.scale);
    while (!quit) {
        SDL_Event e;
        int shot = 0;
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT) quit = 1;
            else if (e.type == SDL_KEYDOWN && !e.key.repeat) {
                switch (e.key.keysym.sym) {
                case SDLK_ESCAPE: case SDLK_q: quit = 1; break;
                case SDLK_SPACE: v.running = !v.running; printf(v.running ? "run\n" : "pause at t %.3f sigma-time\n", run1_time(r)); break;
                case SDLK_d: v.show ^= PID1_SHOW_DISKS; break;
                case SDLK_h: v.show ^= PID1_SHOW_HULL; break;
                case SDLK_g: v.show ^= PID1_SHOW_GRAPH; break;
                case SDLK_p: v.show ^= PID1_SHOW_PARTICLES; break;
                case SDLK_f: v.show ^= PID1_SHOW_FOOD; break;
                case SDLK_c: v.overlay = (v.overlay + 1) % 4; printf("overlay: %s\n", OV[v.overlay]); break;
                case SDLK_LEFTBRACKET:  v.alpha = v.alpha > 16 ? v.alpha - 16 : 0; printf("fill alpha %.2f\n", v.alpha / 255.0); break;
                case SDLK_RIGHTBRACKET: v.alpha = v.alpha < 239 ? v.alpha + 16 : 255; printf("fill alpha %.2f\n", v.alpha / 255.0); break;
                case SDLK_s: shot = 1; break;
                default: break;
                }
            }
        }
        /* a shot due at the current time (t = 0, or a time already reached) is taken on this frame */
        const double now = run1_time(r);
        for (int k = 0; k < nshot; ++k) if (shots[k] <= now + 1e-9) { shot = 1; shots[k] = shots[--nshot]; --k; }
        if (v.running && !r->done && !shot) {
            /* the next frame time: speed sigma-time ahead, or the next shot if that comes first (the frame then shows exactly it) */
            double tt = now + r->v.speed; int at_shot = -1;
            for (int k = 0; k < nshot; ++k) if (shots[k] <= tt) { tt = shots[k]; at_shot = k; }
            run1_advance(r, tt);
            if (at_shot >= 0 && run1_time(r) >= shots[at_shot] - 1e-9) { shot = 1; shots[at_shot] = shots[--nshot]; }
        }
        if (r->done && !announced) { printf("reached T = %.0f sigma-time%s\n", r->T, r->v.exit_at_end ? "" : "; the window stays open (Esc quits)"); announced = 1; }
        if (r->done && !shot)                  /* shots after the end of the run are never due */
            for (int k = 0; k < nshot; ++k) if (shots[k] > run1_time(r) + 1e-9) { printf("no shot at t %.3f: the run ended at %.3f\n", shots[k], run1_time(r)); shots[k] = shots[--nshot]; --k; }
        draw(R, &v, r, fed, &tex);
        if (shot) save_shot(R, r, run1_time(r));
        SDL_RenderPresent(R);
        if (r->done && r->v.exit_at_end && nshot == 0) quit = 1;
        if (!vsync && !r->v.exit_at_end) SDL_Delay(15);
    }
    free(fed); if (tex) SDL_DestroyTexture(tex);
    SDL_DestroyRenderer(R); SDL_DestroyWindow(win); SDL_Quit();
    return 0;
}
