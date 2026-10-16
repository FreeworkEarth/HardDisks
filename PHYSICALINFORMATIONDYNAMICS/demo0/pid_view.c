/* pid_view.c -- the window of brief sec. 5. It draws the World and nothing else. It changes the world only through
   run_tick() (Space runs, N steps one tick) and run_reset() (R); both are the functions the headless path uses, so a window
   run writes the same files as a headless one. Disks are filled by scanlines with alpha (an idea, not code, taken from
   hspist3's gs_disk; M0_AUDIT.md sec. D2). */
#define SDL_MAIN_HANDLED
#include <SDL.h>
#ifdef PID_HAVE_TTF
#include <SDL_ttf.h>
#endif
#include <math.h>
#include <string.h>
#include "pid.h"

typedef struct {
    int show_disks, show_graph, show_hull, show_cue, show_patches;
    int alpha, running;
    double scale; int margin, bar, win_w, win_h;
} View;

static void to_screen(const View* v, const Params* p, double x, double y, float* sx, float* sy){
    *sx = (float)(v->margin + x * v->scale);
    *sy = (float)(v->margin + v->bar + (p->H - y) * v->scale);    /* world y up, screen y down */
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

/* 0b cue overlay: c(x) on a 1-sigma grid into a streaming texture (the idea of hspist3's rm_draw, M0_AUDIT.md sec. D2) */
static void draw_cue(SDL_Renderer* R, const View* v, const World* w, SDL_Texture** tex, int* tw, int* th){
    const Params* p = &w->p;
    const int nx = (int)ceil(p->W), ny = (int)ceil(p->H);
    if (!*tex || *tw != nx || *th != ny) {
        if (*tex) SDL_DestroyTexture(*tex);
        *tex = SDL_CreateTexture(R, SDL_PIXELFORMAT_ARGB8888, SDL_TEXTUREACCESS_STREAMING, nx, ny);
        if (!*tex) return;
        SDL_SetTextureBlendMode(*tex, SDL_BLENDMODE_BLEND); *tw = nx; *th = ny;
    }
    double cmax = 0.0;
    for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) { double c, gx, gy; pid_cue(w, i + 0.5, j + 0.5, &c, &gx, &gy); if (c > cmax) cmax = c; }
    void* px; int pitch;
    if (SDL_LockTexture(*tex, NULL, &px, &pitch)) return;
    for (int j = 0; j < ny; ++j) {
        Uint32* row = (Uint32*)((char*)px + (size_t)(ny - 1 - j) * (size_t)pitch);   /* world y up */
        for (int i = 0; i < nx; ++i) {
            double c, gx, gy; pid_cue(w, i + 0.5, j + 0.5, &c, &gx, &gy);
            const double a = cmax > 0 ? c / cmax : 0.0;
            const Uint32 A = (Uint32)(170.0 * a), Rr = (Uint32)(120 + 135 * a), G = (Uint32)(40 + 120 * a), B = (Uint32)(90 * (1 - a));
            row[i] = (A << 24) | (Rr << 16) | (G << 8) | B;
        }
    }
    SDL_UnlockTexture(*tex);
    float x0, y0, x1, y1; to_screen(v, p, 0.0, 0.0, &x0, &y0); to_screen(v, p, p->W, p->H, &x1, &y1);
    const SDL_FRect dst = {x0, y1, x1 - x0, y0 - y1};
    SDL_RenderCopyF(R, *tex, NULL, &dst);
}

static void draw(SDL_Renderer* R, const View* v, const Run* run, unsigned char* live_tip, unsigned char* fed, SDL_Texture** cue, int* cw, int* ch){
    const World* w = &run->w; const Params* p = &w->p;
    SDL_SetRenderDrawColor(R, 14, 16, 22, 255);
    SDL_RenderClear(R);
    float x0, y0, x1, y1;
    to_screen(v, p, 0.0, 0.0, &x0, &y0); to_screen(v, p, p->W, p->H, &x1, &y1);
    const SDL_FRect box = {x0, y1, x1 - x0, y0 - y1};
    SDL_SetRenderDrawColor(R, 34, 38, 48, 255); SDL_RenderFillRectF(R, &box);
    SDL_SetRenderDrawColor(R, 120, 124, 140, 255); SDL_RenderDrawRectF(R, &box);
    if (v->show_cue && p->npatch > 0) draw_cue(R, v, w, cue, cw, ch);
    if (v->show_patches)
        for (int k = 0; k < p->npatch; ++k) {
            const Patch* q = &p->patch[k]; float sx, sy; to_screen(v, p, q->x, q->y, &sx, &sy);
            if (pid_patch_home(w, k)) SDL_SetRenderDrawColor(R, 90, 210, 130, q->active ? 110 : 40);   /* A, home: green */
            else                      SDL_SetRenderDrawColor(R, 235, 90, 210, q->active ? 130 : 40);   /* food: magenta */
            fill_disk(R, sx, sy, (float)(q->r * v->scale));
        }
    if (p->npatch > 0) world_fed(w, fed); else memset(fed, 0, (size_t)w->nd);
    memset(live_tip, 0, (size_t)w->nd);
    for (int k = 0; k < w->nt; ++k) if (w->t[k].alive) live_tip[w->t[k].disk] = 1;
    if (v->show_disks) {
        for (int i = 0; i < w->nd; ++i) {
            const Disk* d = &w->d[i]; if (!d->alive) continue;
            if (d->kind == KIND_SEED)      SDL_SetRenderDrawColor(R, 225, 95, 80, (Uint8)v->alpha);
            else if (fed[i])               SDL_SetRenderDrawColor(R, 250, 235, 90, (Uint8)v->alpha);    /* 0b: fed path */
            else if (live_tip[i])          SDL_SetRenderDrawColor(R, 250, 150, 60, (Uint8)v->alpha);
            else if (d->kind == KIND_TIP)  SDL_SetRenderDrawColor(R, 150, 152, 165, (Uint8)v->alpha);   /* dead end */
            else                           SDL_SetRenderDrawColor(R, 80, 150, 235, (Uint8)v->alpha);
            float sx, sy; to_screen(v, p, d->x, d->y, &sx, &sy);
            fill_disk(R, sx, sy, (float)(d->r * v->scale));
        }
    }
    if (v->show_graph) {
        SDL_SetRenderDrawColor(R, 240, 240, 240, 230);
        for (int i = 1; i < w->nd; ++i) {
            const Disk* d = &w->d[i]; if (!d->alive || d->parent < 0) continue;
            const Disk* q = &w->d[d->parent];
            float ax, ay, bx, by; to_screen(v, p, d->x, d->y, &ax, &ay); to_screen(v, p, q->x, q->y, &bx, &by);
            SDL_RenderDrawLineF(R, ax, ay, bx, by);
        }
    }
    if (v->show_patches)                               /* patch outlines on top, so they stay visible under the mat */
        for (int k = 0; k < p->npatch; ++k) {
            const Patch* q = &p->patch[k]; float sx, sy; to_screen(v, p, q->x, q->y, &sx, &sy);
            if (pid_patch_home(w, k)) SDL_SetRenderDrawColor(R, 120, 240, 160, 255); else SDL_SetRenderDrawColor(R, 255, 120, 230, 255);
            ring(R, sx, sy, (float)(q->r * v->scale));
        }
    /* the hull (0c) has nothing to draw yet */
}

static void save_shot(SDL_Renderer* R, const Run* run){
    int w = 0, h = 0; SDL_GetRendererOutputSize(R, &w, &h);
    SDL_Surface* s = SDL_CreateRGBSurfaceWithFormat(0, w, h, 24, SDL_PIXELFORMAT_RGB24);   /* 24-bit BI_RGB BMP: readable by sips */
    if (!s) { fprintf(stderr, "pid: shot: %s\n", SDL_GetError()); return; }
    char path[1200]; int k = 0;
    snprintf(path, sizeof path, "%s/shot_%06ld.bmp", run->out, run->w.tick);
    for (FILE* f; (f = fopen(path, "r")) != NULL; ) { fclose(f); snprintf(path, sizeof path, "%s/shot_%06ld_%d.bmp", run->out, run->w.tick, ++k); }
    if (SDL_RenderReadPixels(R, NULL, SDL_PIXELFORMAT_RGB24, s->pixels, s->pitch) || SDL_SaveBMP(s, path))
        fprintf(stderr, "pid: shot: %s\n", SDL_GetError());
    else printf("shot %s (%dx%d)\n", path, w, h);
    SDL_FreeSurface(s);
}

int view_run(Run* run){
    SDL_SetMainReady();
    if (SDL_Init(SDL_INIT_VIDEO) != 0) { fprintf(stderr, "pid: SDL_Init: %s\n", SDL_GetError()); return 2; }
    View v = {1, 1, 0, 0, 1, 200, 1, 0.0, 24, 0, 0, 0};
    const Params* p = &run->p;
    v.scale = fmin((1400.0 - 2.0 * v.margin) / p->W, (860.0 - 2.0 * v.margin) / p->H);
#ifdef PID_HAVE_TTF
    TTF_Font* font = NULL;
    if (run->font) {
        if (TTF_Init() == 0) font = TTF_OpenFont(run->font, 16);
        if (!font) fprintf(stderr, "pid: font %s not loaded (%s); status line in the terminal\n", run->font, TTF_GetError());
        else v.bar = 26;
    }
#else
    if (run->font) fprintf(stderr, "pid: built without SDL2_ttf; --font ignored, status line in the terminal\n");
#endif
    v.win_w = (int)ceil(p->W * v.scale) + 2 * v.margin;
    v.win_h = (int)ceil(p->H * v.scale) + 2 * v.margin + v.bar;
    SDL_Window* win = SDL_CreateWindow("PID Demo 0b -- growth toward food", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, v.win_w, v.win_h, SDL_WINDOW_SHOWN);
    if (!win) { fprintf(stderr, "pid: SDL_CreateWindow: %s\n", SDL_GetError()); SDL_Quit(); return 2; }
    int vsync = 1;
    SDL_Renderer* R = SDL_CreateRenderer(win, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    if (!R) { vsync = 0; R = SDL_CreateRenderer(win, -1, SDL_RENDERER_SOFTWARE); }
    if (!R) { fprintf(stderr, "pid: SDL_CreateRenderer: %s\n", SDL_GetError()); SDL_DestroyWindow(win); SDL_Quit(); return 2; }
    SDL_SetRenderDrawBlendMode(R, SDL_BLENDMODE_BLEND);

    unsigned char* live_tip = NULL; unsigned char* fed = NULL; int cap_live = 0;
    SDL_Texture* cue = NULL; int cue_w = 0, cue_h = 0;
    int quit = 0, shot = 0, announced_end = 0; long last_print = -1; Uint32 last_tick_ms = 0;
    run_status(run, stdout);
    while (!quit) {
        SDL_Event e;
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT) quit = 1;
            else if (e.type == SDL_KEYDOWN && !e.key.repeat) {
                switch (e.key.keysym.sym) {
                case SDLK_ESCAPE: quit = 1; break;
                case SDLK_d: v.show_disks = !v.show_disks; break;
                case SDLK_g: v.show_graph = !v.show_graph; break;
                case SDLK_h: v.show_hull = !v.show_hull; printf("hull: arrives with Demo 0c\n"); break;
                case SDLK_c: v.show_cue = !v.show_cue; if (!run->p.npatch) printf("cue field: no patches in this run\n"); break;
                case SDLK_p: v.show_patches = !v.show_patches; if (!run->p.npatch) printf("patches: none in this run\n"); break;
                case SDLK_LEFTBRACKET:  v.alpha = v.alpha > 31 ? v.alpha - 16 : 15; printf("disk alpha %d\n", v.alpha); break;
                case SDLK_RIGHTBRACKET: v.alpha = v.alpha < 239 ? v.alpha + 16 : 255; printf("disk alpha %d\n", v.alpha); break;
                case SDLK_SPACE: v.running = !v.running; printf(v.running ? "run\n" : "pause\n"); run_status(run, stdout); break;
                case SDLK_n: v.running = 0; if (!run_tick(run)) printf("at --steps=%ld\n", run->steps); run_status(run, stdout); break;
                case SDLK_r: run_reset(run); announced_end = 0; printf("reset to tick 0, seed %llu\n", (unsigned long long)run->seed); run_status(run, stdout); break;
                case SDLK_s: shot = 1; break;
                default: break;
                }
            }
        }
        const Uint32 now_ms = SDL_GetTicks();   /* pacing only: how often ticks run, never what they do */
        if (v.running && (run->tick_ms == 0 || now_ms - last_tick_ms >= (Uint32)run->tick_ms)) {
            last_tick_ms = now_ms;
            for (int k = 0; k < run->tpf; ++k) if (!run_tick(run)) break;
        }
        if (run->w.tick != last_print && (run->w.tick % 10 == 0 || run->w.tick == run->steps)) { run_status(run, stdout); last_print = run->w.tick; }
        const int at_end = run->w.tick >= run->steps;
        if (at_end && !announced_end) { printf("reached --steps=%ld; window stays open (Esc quits)\n", run->steps); announced_end = 1; if (run->shot_at_end) shot = 1; }
        if (run->w.nd > cap_live) {
            cap_live = run->w.nd * 2;
            live_tip = (unsigned char*)realloc(live_tip, (size_t)cap_live); fed = (unsigned char*)realloc(fed, (size_t)cap_live);
            if (!live_tip || !fed) { fprintf(stderr, "pid: out of memory\n"); return 3; }
        }
        draw(R, &v, run, live_tip, fed, &cue, &cue_w, &cue_h);
#ifdef PID_HAVE_TTF
        if (font) {
            char s[256]; snprintf(s, sizeof s, "tick %ld   live tips %d   disks %d   M_used %.2f   M_pool %.2f%s", run->w.tick, world_live_tips(&run->w),
                                  world_alive_disks(&run->w), run->w.m_used, world_pool(&run->w), v.running ? "" : "   (paused)");
            const SDL_Color c = {230, 230, 235, 255};
            SDL_Surface* ts = TTF_RenderUTF8_Blended(font, s, c);
            if (ts) { SDL_Texture* tx = SDL_CreateTextureFromSurface(R, ts); const SDL_Rect dst = {v.margin, 4, ts->w, ts->h}; SDL_RenderCopy(R, tx, NULL, &dst); SDL_DestroyTexture(tx); SDL_FreeSurface(ts); }
        }
#endif
        if (shot) { save_shot(R, run); shot = 0; }
        SDL_RenderPresent(R);
        if (at_end && run->exit_at_end) quit = 1;
        if (!vsync) SDL_Delay(run->exit_at_end ? 0 : 15);
    }
    free(live_tip); free(fed); if (cue) SDL_DestroyTexture(cue);
#ifdef PID_HAVE_TTF
    if (font) TTF_CloseFont(font);
    if (run->font) TTF_Quit();
#endif
    SDL_DestroyRenderer(R); SDL_DestroyWindow(win); SDL_Quit();
    return 0;
}
