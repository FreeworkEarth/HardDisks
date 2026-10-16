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

static void draw(SDL_Renderer* R, const View* v, const Run* run, unsigned char* live_tip){
    const World* w = &run->w; const Params* p = &w->p;
    SDL_SetRenderDrawColor(R, 14, 16, 22, 255);
    SDL_RenderClear(R);
    float x0, y0, x1, y1;
    to_screen(v, p, 0.0, 0.0, &x0, &y0); to_screen(v, p, p->W, p->H, &x1, &y1);
    const SDL_FRect box = {x0, y1, x1 - x0, y0 - y1};
    SDL_SetRenderDrawColor(R, 34, 38, 48, 255); SDL_RenderFillRectF(R, &box);
    SDL_SetRenderDrawColor(R, 120, 124, 140, 255); SDL_RenderDrawRectF(R, &box);
    memset(live_tip, 0, (size_t)w->nd);
    for (int k = 0; k < w->nt; ++k) if (w->t[k].alive) live_tip[w->t[k].disk] = 1;
    if (v->show_disks) {
        for (int i = 0; i < w->nd; ++i) {
            const Disk* d = &w->d[i]; if (!d->alive) continue;
            if (d->kind == KIND_SEED)      SDL_SetRenderDrawColor(R, 225, 95, 80, (Uint8)v->alpha);
            else if (live_tip[i])          SDL_SetRenderDrawColor(R, 250, 190, 70, (Uint8)v->alpha);
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
    /* patches (0b), cue field (0b) and hull (0c) have nothing to draw in Demo 0a */
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
    SDL_Window* win = SDL_CreateWindow("PID Demo 0a.1 -- unbiased growth", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, v.win_w, v.win_h, SDL_WINDOW_SHOWN);
    if (!win) { fprintf(stderr, "pid: SDL_CreateWindow: %s\n", SDL_GetError()); SDL_Quit(); return 2; }
    int vsync = 1;
    SDL_Renderer* R = SDL_CreateRenderer(win, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    if (!R) { vsync = 0; R = SDL_CreateRenderer(win, -1, SDL_RENDERER_SOFTWARE); }
    if (!R) { fprintf(stderr, "pid: SDL_CreateRenderer: %s\n", SDL_GetError()); SDL_DestroyWindow(win); SDL_Quit(); return 2; }
    SDL_SetRenderDrawBlendMode(R, SDL_BLENDMODE_BLEND);

    unsigned char* live_tip = NULL; int cap_live = 0;
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
                case SDLK_c: v.show_cue = !v.show_cue; printf("cue field: arrives with Demo 0b\n"); break;
                case SDLK_p: v.show_patches = !v.show_patches; printf("patches: arrive with Demo 0b\n"); break;
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
        if (run->w.nd > cap_live) { cap_live = run->w.nd * 2; live_tip = (unsigned char*)realloc(live_tip, (size_t)cap_live); if (!live_tip) { fprintf(stderr, "pid: out of memory\n"); return 3; } }
        draw(R, &v, run, live_tip);
#ifdef PID_HAVE_TTF
        if (font) {
            char s[256]; snprintf(s, sizeof s, "tick %ld   live tips %d   disks %d   M_used %.2f   M_pool %.2f%s", run->w.tick, world_live_tips(&run->w),
                                  run->w.nd, run->w.m_used, world_pool(&run->w), v.running ? "" : "   (paused)");
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
    free(live_tip);
#ifdef PID_HAVE_TTF
    if (font) TTF_CloseFont(font);
    if (run->font) TTF_Quit();
#endif
    SDL_DestroyRenderer(R); SDL_DestroyWindow(win); SDL_Quit();
    return 0;
}
