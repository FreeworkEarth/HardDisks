/* pid_view1.c -- the window of Demo 1 (02_DEMO1_BRIEF.md sec. 9): the frozen Demo 0 network and the hard disks running in
   its lumen. It draws and never owns state: the run moves only through run1_advance (pid_demo1.c), the function the headless
   path uses, so a window run writes the same pstate_*.txt, density_*.txt and measure.csv (make test, test 11). The drawing is
   pid_draw1.c (shared with the studio, pid_studio.c). Layers, bottom to top: the coarse-grid overlay (C), the structural disks (D, fill alpha [ ]), the
   food patches (F), the graph (G), the particles (P, coloured by speed), the hull (H). */
#define SDL_MAIN_HANDLED
#include <SDL.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pid_draw1.h"

typedef struct { Pid1Vp vp; Pid1Scene sc; int running, margin, win_w, win_h; } V1;

static void draw(SDL_Renderer* R, const V1* v, Run1* r, SDL_Texture** tex){
    SDL_SetRenderDrawColor(R, 14, 16, 22, 255);
    SDL_RenderClear(R);
    pd_scene(R, &v->vp, &r->w, r, &v->sc, tex);
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
    v.margin = 24; v.running = 1; v.sc.overlay = r->v.overlay;
    v.sc.alpha = (int)lround(r->v.alpha * 255.0);
    v.sc.show = PID1_SHOW_DISKS | PID1_SHOW_PARTICLES | PID1_SHOW_FOOD;
    if (r->v.view == 1) { v.sc.alpha = 0; v.sc.show = PID1_SHOW_HULL | PID1_SHOW_GRAPH | PID1_SHOW_PARTICLES | PID1_SHOW_FOOD; }   /* --view=network */
    if (r->v.view == 2) v.sc.show |= PID1_SHOW_HULL;                                                                               /* --view=both */
    if (r->v.show_given) v.sc.show = r->v.show;
    const World* w = &r->w;
    v.vp.scale = fmin((1400.0 - 2.0 * v.margin) / w->p.W, (860.0 - 2.0 * v.margin) / w->p.H);
    v.vp.ox = v.margin; v.vp.oy = v.margin; v.vp.H = w->p.H;
    v.win_w = (int)ceil(w->p.W * v.vp.scale) + 2 * v.margin;
    v.win_h = (int)ceil(w->p.H * v.vp.scale) + 2 * v.margin;
    SDL_Window* win = SDL_CreateWindow("PID Demo 1 -- hard disks in the frozen lumen", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, v.win_w, v.win_h, SDL_WINDOW_SHOWN);
    if (!win) { fprintf(stderr, "pid_view1: SDL_CreateWindow: %s\n", SDL_GetError()); SDL_Quit(); return 2; }
    int vsync = 1;
    SDL_Renderer* R = SDL_CreateRenderer(win, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    if (!R) { vsync = 0; R = SDL_CreateRenderer(win, -1, SDL_RENDERER_SOFTWARE); }
    if (!R) { fprintf(stderr, "pid_view1: SDL_CreateRenderer: %s\n", SDL_GetError()); SDL_DestroyWindow(win); SDL_Quit(); return 2; }
    SDL_SetRenderDrawBlendMode(R, SDL_BLENDMODE_BLEND);
    SDL_Texture* tex = NULL;
    int nshot = r->v.nshot; double shots[16]; memcpy(shots, r->v.shot_at, sizeof shots);
    int quit = 0, announced = 0;
    static const char* OV[4] = {"off", "particle density", "particle kinetic energy", "structural material"};
    printf("window %dx%d, %.2f px per sigma; Space run/pause, Esc quit, S shot, D H G P F layers, [ ] alpha, C overlay\n", v.win_w, v.win_h, v.vp.scale);
    while (!quit) {
        SDL_Event e;
        int shot = 0;
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT) quit = 1;
            else if (e.type == SDL_KEYDOWN && !e.key.repeat) {
                switch (e.key.keysym.sym) {
                case SDLK_ESCAPE: case SDLK_q: quit = 1; break;
                case SDLK_SPACE: v.running = !v.running; printf(v.running ? "run\n" : "pause at t %.3f sigma-time\n", run1_time(r)); break;
                case SDLK_d: v.sc.show ^= PID1_SHOW_DISKS; break;
                case SDLK_h: v.sc.show ^= PID1_SHOW_HULL; break;
                case SDLK_g: v.sc.show ^= PID1_SHOW_GRAPH; break;
                case SDLK_p: v.sc.show ^= PID1_SHOW_PARTICLES; break;
                case SDLK_f: v.sc.show ^= PID1_SHOW_FOOD; break;
                case SDLK_c: v.sc.overlay = (v.sc.overlay + 1) % 4; printf("overlay: %s\n", OV[v.sc.overlay]); break;
                case SDLK_LEFTBRACKET:  v.sc.alpha = v.sc.alpha > 16 ? v.sc.alpha - 16 : 0; printf("fill alpha %.2f\n", v.sc.alpha / 255.0); break;
                case SDLK_RIGHTBRACKET: v.sc.alpha = v.sc.alpha < 239 ? v.sc.alpha + 16 : 255; printf("fill alpha %.2f\n", v.sc.alpha / 255.0); break;
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
        draw(R, &v, r, &tex);
        if (shot) save_shot(R, r, run1_time(r));
        SDL_RenderPresent(R);
        if (r->done && r->v.exit_at_end && nshot == 0) quit = 1;
        if (!vsync && !r->v.exit_at_end) SDL_Delay(15);
    }
    if (tex) SDL_DestroyTexture(tex);
    SDL_DestroyRenderer(R); SDL_DestroyWindow(win); SDL_Quit();
    return 0;
}
