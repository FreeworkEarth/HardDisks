/* pid_studio.c -- the PID studio: grow a Demo 0 network and run Demo 1's hard disks in it, from one window with panels in the
   manner of HardDisks' F2/F3 panels (hspist3/00ALLINONE.c: sliders, APPLY rebuilds). Three panels:
     1 GROW (Demo 0)        preset, rule, seed, ticks and the growth parameters; APPLY grows a new network from tick 0, shown tick
                            by tick (Demo 0's own code: demo0/pid_grow.c step(), linked unchanged)
     2 PARTICLES (Demo 1)   the network at a chosen tick, packing fraction (or N), seed, T, kT, the stop interval; APPLY computes the
                            hull of that tick, writes it as a Demo 0 state file and starts pid_demo1's run on it (run1_setup with
                            the same arguments pid_demo1 would get), drawn as it runs
     3 VIEW                 speed, fill alpha, overlay and layers; they act at once and never change a run
   The studio owns no physics: the growth is demo0's step(), the hull demo0's hull_compute_world, the particles pid_demo1.c's run
   (run1_advance, as headless). So a network grown here equals Demo 0 headless with the same flags (make test, test 12: every
   preset), and a run started here writes the same pstate_*.txt, density_*.txt and measure.csv as pid_demo1 (test 12). Every
   APPLY prints the equivalent headless commands. Files go to --out (default runs/studio_<UTC time>/): net_<n>_<tick>/final.txt,
   run_<n>/ (pid_demo1's outputs) and shot_*.bmp (S, or --shot-*; tests/bmp2png.py converts them).
   Text: SDL2_ttf with HardDisks' font (../../hspist3/Roboto-Regular.ttf, read only) or a system font; --font=<path> chooses. */
#define _POSIX_C_SOURCE 200809L
#define SDL_MAIN_HANDLED
#include <SDL.h>
#include <SDL_ttf.h>
#include <errno.h>
#include <inttypes.h>
#include <stdarg.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include "pid_draw1.h"

/* ------------------------------------------------------------------ Demo 0's presets, copied from demo0/pid_main.c:67-86 (the
   function is static there); test 12 grows every preset here and in Demo 0 headless and compares the state files */
static const char* const PRESETS[] = {"0a", "0a1", "mat", "mat45", "mat45r"};
enum { NPRESET = 5 };
static void apply_preset(int k, Params* p){
    params_default(p);
    if (k >= 2) { p->pbranch = 1.0; p->pbranch_exp = 0.0; p->delta = 1.57; p->ntry = 48; p->budget = 9000.0; }                 /* mat */
    if (k == 0) { p->rings = 0; p->pbranch_exp = 0.0; p->pbranch = 0.08; p->ntry = 8; p->budget = 3000.0; }                    /* 0a */
    if (k >= 3) {                                                                                                               /* mat45 */
        p->budget = 4500.0; p->tadapt = 10; p->budget_wait = 1; p->tip_order = 1;
        p->npatch = 2;
        p->patch[0] = (Patch){0.5 * p->W, 0.5 * p->H, 5.0, 1.0, 1};
        p->patch[1] = (Patch){0.8 * p->W, 0.5 * p->H, 5.0, 1.0, 1};
    }
    if (k == 4) { p->thin_all = 1; p->regrow_tips = 6; p->regrow_fan = PID_PI; p->cue_off_reached = 1; }                         /* mat45r */
}
static const char* const RULES[] = {"unbiased", "guided"};
static const char* const OVERLAYS[] = {"off", "density", "kinetic energy", "material"};

/* ------------------------------------------------------------------ the studio's state */
typedef struct {
    /* 1 GROW */
    int preset, rule; double gseed, steps, tpf;
    double budget, kappa, tau, pbranch, beta, lambda, tadapt, gamma, deltathin;
    World w; int have_w, growing; double w_ms;    /* the network being grown; the time of its last step batch */
    long grow_n;                                  /* APPLY count of the growth */
    Params gp; uint64_t gs; char grule[16];       /* the parameters of the network in w */
    long hull_tick;                               /* the tick w.hull was computed for (-1: none) */
    /* 2 PARTICLES */
    double ptick, eta, pseed, T, wevery, kbt; long N_cli;
    Run1 run; int have_run, reported; long run_n;
    char argbuf[10][1400]; char* argv[10];
    /* 3 VIEW */
    double speed, alpha; int overlay, show, running;
    /* files, messages */
    char root[1024]; char msg[512]; Uint32 msg_ms;
    int auto_mode, grow_only, exit_at_end; long shot_tick; double shot_at; int shot_tick_done, shot_at_done;
} Studio;

static void say(Studio* st, const char* fmt, ...) __attribute__((format(printf, 2, 3)));
static void say(Studio* st, const char* fmt, ...){
    va_list ap; va_start(ap, fmt); vsnprintf(st->msg, sizeof st->msg, fmt, ap); va_end(ap);
    st->msg_ms = SDL_GetTicks(); printf("studio: %s\n", st->msg); fflush(stdout);
}
static int mkdir_p(const char* path){
    char buf[1024]; snprintf(buf, sizeof buf, "%s", path);
    for (char* p = buf + 1; *p; ++p) if (*p == '/') { *p = '\0'; if (mkdir(buf, 0755) && errno != EEXIST) return -1; *p = '/'; }
    return mkdir(buf, 0755) && errno != EEXIST ? -1 : 0;
}

/* ------------------------------------------------------------------ growth */
static void grow_params(const Studio* st, Params* p){
    apply_preset(st->preset, p);
    p->budget = st->budget; p->kappa = st->kappa; p->tau = st->tau; p->pbranch = st->pbranch; p->beta = st->beta; p->lambda = st->lambda;
    p->tadapt = (int)lround(st->tadapt); p->gamma = st->gamma; p->deltathin = st->deltathin;
}
static void load_preset_sliders(Studio* st){
    Params p; apply_preset(st->preset, &p);
    st->budget = p.budget; st->kappa = p.kappa; st->tau = p.tau; st->pbranch = p.pbranch; st->beta = p.beta; st->lambda = p.lambda;
    st->tadapt = p.tadapt; st->gamma = p.gamma; st->deltathin = p.deltathin;
}
/* the Demo 0 command that grows the same network (every value the panel sets, %.17g) */
static void demo0_command(const Studio* st, long ticks, char* buf, size_t n){
    snprintf(buf, n, "bin/demo0/pid_demo0_headless --preset=%s --rule=%s --seed=%" PRIu64 " --steps=%ld --budget=%.17g --kappa=%.17g --tau=%.17g"
             " --pbranch=%.17g --beta=%.17g --lambda=%.17g --tadapt=%d --gamma=%.17g --deltathin=%.17g --headless --out=<new folder>",
             PRESETS[st->preset], st->grule, st->gs, ticks, st->gp.budget, st->gp.kappa, st->gp.tau, st->gp.pbranch, st->gp.beta,
             st->gp.lambda, st->gp.tadapt, st->gp.gamma, st->gp.deltathin);
}
static int grow_start(Studio* st){
    Params p; grow_params(st, &p);
    char err[256];
    if (params_check(&p, err, sizeof err)) { say(st, "GROW refused: %s", err); return 0; }
    if (st->have_w) world_free(&st->w);
    st->gp = p; st->gs = (uint64_t)llround(st->gseed); snprintf(st->grule, sizeof st->grule, "%s", RULES[st->rule]);
    world_init(&st->w, &st->gp, st->gs, st->grule);
    st->have_w = 1; st->growing = 1; st->hull_tick = -1; st->grow_n++;
    char cmd[1024]; demo0_command(st, (long)llround(st->steps), cmd, sizeof cmd);
    say(st, "GROW %ld: preset %s, rule %s, seed %" PRIu64 ", %ld ticks", st->grow_n, PRESETS[st->preset], st->grule, st->gs, (long)llround(st->steps));
    printf("studio: the same network headless: %s\n", cmd);
    return 1;
}
/* grow w forward to tick k (k >= w.tick), or regrow it from tick 0 when k is behind */
static void grow_to(Studio* st, long k){
    if (!st->have_w) return;
    if (st->w.tick > k) { world_free(&st->w); world_init(&st->w, &st->gp, st->gs, st->grule); st->hull_tick = -1; }
    const Policy pol = policy_by_name(st->grule);
    while (st->w.tick < k) step(&st->w, pol);
}
static void grow_hull(Studio* st){
    if (st->hull_tick == st->w.tick) return;
    hull_free(&st->w.hull); hull_compute_world(&st->w, &st->w.hull); st->hull_tick = st->w.tick;
}
/* the network at its current tick into <root>/net_<n>_<tick>/final.txt, with its hull, as Demo 0 writes final.txt */
static int write_net(Studio* st, char* path, size_t n){
    char dir[1100]; snprintf(dir, sizeof dir, "%s/net_%ld_%06ld", st->root, st->grow_n, st->w.tick);
    if (mkdir_p(dir)) { say(st, "cannot create %s", dir); return 0; }
    snprintf(path, n, "%s/final.txt", dir);
    grow_hull(st);
    if (state_write(&st->w, path)) { say(st, "cannot write %s", path); return 0; }
    char err[256];
    const int fails = world_check(&st->w, st->w.p.tadapt == 0, world_pool(&st->w), err, sizeof err);   /* as Demo 0 at every write */
    if (fails) say(st, "Demo 0's invariant check: %d failures at tick %ld: %s", fails, st->w.tick, err);
    return 1;
}

/* ------------------------------------------------------------------ particles */
static void run_end(Studio* st){
    if (!st->have_run) return;
    if (!st->reported) { run1_report(&st->run); st->reported = 1; }
    run1_free(&st->run); st->have_run = 0;
}
static int particles_start(Studio* st){
    if (!st->have_w) { say(st, "PARTICLES: grow a network first (panel 1, APPLY)"); return 0; }
    run_end(st);
    const long k = (long)llround(st->ptick);
    if (k < 1) { say(st, "PARTICLES: tick 0 has no lumen (only the seed)"); return 0; }
    st->growing = 0;
    grow_to(st, k);
    grow_hull(st);
    if (!st->w.hull.valid || st->w.hull.np == 0 || !(st->w.hull.area > 0.0)) { say(st, "PARTICLES: the network at tick %ld has no lumen", k); return 0; }
    char net[1200]; if (!write_net(st, net, sizeof net)) return 0;
    const double area1 = PID_PI * 0.25;                                   /* one particle, R = 0.5 sigma */
    long N = st->N_cli > 0 ? st->N_cli : lround(st->eta * st->w.hull.area / area1);
    if (st->N_cli > 0) { st->eta = N * area1 / st->w.hull.area; st->N_cli = 0; }   /* the command line's N, once */
    if (N < 1) N = 1;
    st->run_n++;
    snprintf(st->argbuf[0], sizeof st->argbuf[0], "pid_studio");
    snprintf(st->argbuf[1], sizeof st->argbuf[1], "--state=%s", net);
    snprintf(st->argbuf[2], sizeof st->argbuf[2], "--N=%ld", N);
    snprintf(st->argbuf[3], sizeof st->argbuf[3], "--T=%ld", (long)llround(st->T));
    snprintf(st->argbuf[4], sizeof st->argbuf[4], "--seed=%" PRIu64, (uint64_t)llround(st->pseed));
    snprintf(st->argbuf[5], sizeof st->argbuf[5], "--write-every=%ld", (long)llround(st->wevery));
    snprintf(st->argbuf[6], sizeof st->argbuf[6], "--kbt=%.17g", st->kbt);
    snprintf(st->argbuf[7], sizeof st->argbuf[7], "--out=%s/run_%ld", st->root, st->run_n);
    for (int i = 0; i < 8; ++i) st->argv[i] = st->argbuf[i];
    st->argv[8] = NULL;
    say(st, "PARTICLES %ld: tick %ld, lumen %.1f sigma^2, N %ld (eta %.3f), seed %s, T %ld", st->run_n, k, st->w.hull.area, N,
        N * area1 / st->w.hull.area, st->argbuf[4] + 7, (long)llround(st->T));
    printf("studio: the same run headless: bin/pid_demo1");
    for (int i = 1; i < 8; ++i) printf(" %s", st->argv[i]);
    printf("\n");
    const int rc = run1_setup(&st->run, 8, st->argv, 1);
    if (rc) { run1_free(&st->run); say(st, "PARTICLES refused (see the terminal): eta %.3f may be too dense for this lumen", N * area1 / st->w.hull.area); return 0; }
    st->have_run = 1; st->reported = 0; st->running = 1;
    return 1;
}

/* ------------------------------------------------------------------ widgets */
enum { WK_SLIDER, WK_CYCLE, WK_TOGGLES, WK_BUTTONS };
enum { A_NONE, A_GROW_PAUSE, A_PART_STOP, A_RUN_PAUSE, A_SHOT, A_QUIT };
typedef struct {
    int kind; const char* label;
    double* val; double lo, hi; int logs, integer; int live;   /* slider */
    int* ival; const char* const* names; int n;                /* cycle; toggles (bits 1 << k); buttons (actions act[k]) */
    int act[6];
    SDL_Rect r, track, sub[6];
} Widget;
typedef struct { const char* title; Widget w[20]; int n; int open, dirty; SDL_Rect title_r, apply_r, sec_r; int has_apply, sec_act; } Panel;

static TTF_Font *font, *font_b;
static float ui = 1.0f;
static int text(SDL_Renderer* R, TTF_Font* f, int x, int y, SDL_Color c, const char* s){
    if (!s[0]) return 0;
    SDL_Surface* ts = TTF_RenderUTF8_Blended(f, s, c);
    if (!ts) return 0;
    SDL_Texture* tx = SDL_CreateTextureFromSurface(R, ts);
    const SDL_Rect d = {x, y, ts->w, ts->h}; const int w = ts->w;
    if (tx) { SDL_RenderCopy(R, tx, NULL, &d); SDL_DestroyTexture(tx); }
    SDL_FreeSurface(ts);
    return w;
}
static int text_w(TTF_Font* f, const char* s){ int w = 0, h = 0; TTF_SizeUTF8(f, s, &w, &h); return w; }
static void fill(SDL_Renderer* R, SDL_Rect r, int cr, int cg, int cb, int a){ SDL_SetRenderDrawColor(R, (Uint8)cr, (Uint8)cg, (Uint8)cb, (Uint8)a); SDL_RenderFillRect(R, &r); }
static void frame(SDL_Renderer* R, SDL_Rect r, int cr, int cg, int cb){ SDL_SetRenderDrawColor(R, (Uint8)cr, (Uint8)cg, (Uint8)cb, 255); SDL_RenderDrawRect(R, &r); }
static int U(double x){ return (int)lround(x * ui); }

static double slider_frac(const Widget* w){
    const double v = *w->val;
    double t = w->logs ? log(v / w->lo) / log(w->hi / w->lo) : (v - w->lo) / (w->hi - w->lo);
    return t < 0.0 ? 0.0 : t > 1.0 ? 1.0 : t;
}
static void slider_set(Widget* w, double t){
    t = t < 0.0 ? 0.0 : t > 1.0 ? 1.0 : t;
    double v = w->logs ? w->lo * pow(w->hi / w->lo, t) : w->lo + t * (w->hi - w->lo);
    if (w->integer) v = round(v);
    *w->val = v < w->lo ? w->lo : v > w->hi ? w->hi : v;
}
static void slider_nudge(Widget* w, int dir){
    double v = *w->val;
    if (w->integer) v += dir;
    else if (w->logs) v *= pow(1.02, dir);
    else v += dir * (w->hi - w->lo) / 200.0;
    *w->val = v < w->lo ? w->lo : v > w->hi ? w->hi : v;
}
static void fmt_val(const Widget* w, char* b, size_t n){
    if (w->integer) snprintf(b, n, "%.0f", *w->val); else snprintf(b, n, "%.4g", *w->val);
}

/* lays a panel out from (x, y) with width pw; returns the y below it. One row per widget: a slider is its label (left), its track
   and its value (right); the APPLY row also holds the panel's second button */
static int panel_layout(Panel* P, int x, int y, int pw){
    const int pad = U(10), rh = U(26);
    P->title_r = (SDL_Rect){x, y, pw, U(28)}; y += U(28);
    if (!P->open) return y + U(6);
    y += U(4);
    const int iw = pw - 2 * pad, lw = (int)(0.47 * iw), vw = U(62);
    for (int i = 0; i < P->n; ++i) {
        Widget* w = &P->w[i];
        w->r = (SDL_Rect){x + pad, y, iw, rh};
        if (w->kind == WK_SLIDER) w->track = (SDL_Rect){x + pad + lw, y + U(9), iw - lw - vw - U(8), U(8)};
        else if (w->kind == WK_CYCLE) w->sub[0] = (SDL_Rect){x + pad + lw, y + U(1), iw - lw, U(23)};
        else {
            const int n = w->n, bw = (iw - (n - 1) * U(4)) / (n ? n : 1);
            for (int k = 0; k < n; ++k) w->sub[k] = (SDL_Rect){x + pad + k * (bw + U(4)), y + U(1), bw, U(23)};
        }
        y += rh;
    }
    if (P->has_apply) {
        const int sw = P->sec_act ? (iw - U(6)) / 2 : 0;
        P->apply_r = (SDL_Rect){x + pad, y + U(4), iw - sw - (sw ? U(6) : 0), U(28)};
        P->sec_r = (SDL_Rect){x + pad + iw - sw, y + U(4), sw, U(28)};
        y += U(36);
    }
    return y + U(6);
}
static void button(SDL_Renderer* R, SDL_Rect r, int on, const char* label, TTF_Font* f){
    const SDL_Color tc = {225, 228, 236, 255};
    fill(R, r, on ? 60 : 36, on ? 110 : 40, on ? 90 : 54, 255); frame(R, r, 90, 96, 116);
    text(R, f, r.x + (r.w - text_w(f, label)) / 2, r.y + (r.h - TTF_FontHeight(f)) / 2, tc, label);
}
static const char* sec_label(const Studio* st, int act){
    if (act == A_GROW_PAUSE) return st->growing ? "pause growth" : "resume growth";
    if (act == A_PART_STOP) return "stop run";
    return "";
}
static void panel_draw(SDL_Renderer* R, const Panel* P, const Studio* st){
    const SDL_Color tc = {225, 228, 236, 255}, dim = {150, 155, 170, 255}, val = {255, 214, 120, 255};
    fill(R, P->title_r, 44, 50, 66, 255);
    char t[96]; snprintf(t, sizeof t, "%s %s", P->open ? "-" : "+", P->title);
    text(R, font_b, P->title_r.x + U(8), P->title_r.y + U(5), tc, t);
    if (!P->open) return;
    for (int i = 0; i < P->n; ++i) {
        const Widget* w = &P->w[i];
        const int ty = w->r.y + (w->r.h - TTF_FontHeight(font)) / 2;
        if (w->kind == WK_SLIDER) {
            char vb[48]; fmt_val(w, vb, sizeof vb);
            text(R, font, w->r.x, ty, w->live ? tc : dim, w->label);
            text(R, font, w->r.x + w->r.w - text_w(font, vb), ty, val, vb);
            fill(R, w->track, 30, 33, 42, 255); frame(R, w->track, 80, 86, 104);
            SDL_Rect f = w->track; f.w = (int)lround(w->track.w * slider_frac(w)); fill(R, f, 70, 120, 200, 255);
            const SDL_Rect h = {w->track.x + f.w - U(4), w->track.y - U(5), U(8), w->track.h + U(10)}; fill(R, h, 220, 225, 235, 255);
        } else if (w->kind == WK_CYCLE) {
            text(R, font, w->r.x, ty, w->live ? tc : dim, w->label);
            fill(R, w->sub[0], 36, 40, 54, 255); frame(R, w->sub[0], 90, 96, 116);
            char b[96]; snprintf(b, sizeof b, "< %s >", w->names[*w->ival]);
            text(R, font, w->sub[0].x + (w->sub[0].w - text_w(font, b)) / 2, w->sub[0].y + (w->sub[0].h - TTF_FontHeight(font)) / 2, val, b);
        } else {
            for (int k = 0; k < w->n; ++k) {
                int on = 0; const char* lab = w->names[k];
                if (w->kind == WK_TOGGLES) on = (*w->ival >> k) & 1;
                if (w->kind == WK_BUTTONS && w->act[k] == A_RUN_PAUSE) { on = st->running; lab = st->running ? "pause" : "run"; }
                button(R, w->sub[k], on, lab, font);
            }
        }
    }
    if (P->has_apply) {
        if (P->dirty) fill(R, P->apply_r, 220, 130, 50, 255); else fill(R, P->apply_r, 60, 70, 90, 255);
        frame(R, P->apply_r, 140, 146, 166);
        const char* a = P->dirty ? "APPLY (changed)" : "APPLY";
        text(R, font_b, P->apply_r.x + (P->apply_r.w - text_w(font_b, a)) / 2, P->apply_r.y + (P->apply_r.h - TTF_FontHeight(font_b)) / 2, tc, a);
        if (P->sec_act) button(R, P->sec_r, P->sec_act == A_GROW_PAUSE && st->growing, sec_label(st, P->sec_act), font);
    }
}
static int inside(SDL_Rect r, int x, int y){ return x >= r.x && x < r.x + r.w && y >= r.y && y < r.y + r.h; }

/* ------------------------------------------------------------------ screenshots */
static void save_shot(SDL_Renderer* R, Studio* st, const char* tag){
    int w = 0, h = 0; SDL_GetRendererOutputSize(R, &w, &h);
    SDL_Surface* s = SDL_CreateRGBSurfaceWithFormat(0, w, h, 24, SDL_PIXELFORMAT_RGB24);
    if (!s) return;
    char path[1300]; int k = 0;
    snprintf(path, sizeof path, "%s/shot_%s.bmp", st->root, tag);
    for (FILE* f; (f = fopen(path, "r")) != NULL; ) { fclose(f); snprintf(path, sizeof path, "%s/shot_%s_%d.bmp", st->root, tag, ++k); }
    if (SDL_RenderReadPixels(R, NULL, SDL_PIXELFORMAT_RGB24, s->pixels, s->pitch) || SDL_SaveBMP(s, path)) say(st, "shot failed: %s", SDL_GetError());
    else say(st, "shot %s (%dx%d)", path, w, h);
    SDL_FreeSurface(s);
}

/* ------------------------------------------------------------------ command line */
static int num_d(const char* v, double* x){ char* e; errno = 0; *x = strtod(v, &e); return errno || e == v || *e || !isfinite(*x); }
static void usage(FILE* f){
    fprintf(f, "usage: pid_studio [--out=<new folder>] [--font=<ttf>] [growth: --preset=0a|0a1|mat|mat45|mat45r --rule=unbiased|guided\n"
               "       --grow-seed=S --steps=K --tpf=<ticks per frame> --budget= --kappa= --tau= --pbranch= --beta= --lambda= --tadapt=\n"
               "       --gamma= --deltathin=] [particles: --tick=k --eta=x | --N=n --seed=S --T=t --write-every=w --kbt=x]\n"
               "       [view: --speed=<sigma-time per frame> --alpha=0..1 --overlay=off|density|ke|material --show=disks,hull,graph,particles,food]\n"
               "       [--auto (APPLY growth, grow to --steps, APPLY particles at --tick) [--grow-only] --exit-at-end --shot-tick=k --shot-at=t]\n"
               "Every value starts a panel's slider; the panels' APPLY buttons start growth and particles. Keys: Space run/pause,\n"
               "S screenshot, D H G P F layers, C overlay, [ ] fill alpha, Esc quit.\n");
}
static int parse(Studio* st, int argc, char** argv, const char** out, const char** fontp){
    for (int i = 1; i < argc; ++i)                        /* the preset first, so that the explicit values override it (as Demo 0) */
        if (!strncmp(argv[i], "--preset=", 9)) {
            int k = -1; for (int j = 0; j < NPRESET; ++j) if (!strcmp(argv[i] + 9, PRESETS[j])) k = j;
            if (k < 0) { fprintf(stderr, "pid_studio: unknown preset '%s' (0a, 0a1, mat, mat45, mat45r)\n", argv[i] + 9); return 2; }
            st->preset = k; load_preset_sliders(st);
        }
    for (int i = 1; i < argc; ++i) {
        const char* a = argv[i]; const char* eq = strchr(a, '=');
        if (!strcmp(a, "--help")) { usage(stdout); exit(0); }
        if (!strcmp(a, "--auto")) { st->auto_mode = 1; continue; }
        if (!strcmp(a, "--grow-only")) { st->grow_only = 1; continue; }
        if (!strcmp(a, "--exit-at-end")) { st->exit_at_end = 1; continue; }
        if (strncmp(a, "--", 2) || !eq) { fprintf(stderr, "pid_studio: bad argument '%s'\n", a); usage(stderr); return 2; }
        const size_t kl = (size_t)(eq - a); const char* v = eq + 1; int bad = 0; double x = 0.0;
        #define KEY(k) (kl == strlen(k) && !strncmp(a, k, kl))
        if      (KEY("--out"))  *out = v;
        else if (KEY("--font")) *fontp = v;
        else if (KEY("--preset")) { }   /* applied above */
        else if (KEY("--rule"))   { if (!strcmp(v, "unbiased")) st->rule = 0; else if (!strcmp(v, "guided")) st->rule = 1; else bad = 1; }
        else if (KEY("--overlay")) { if (!strcmp(v, "off")) st->overlay = 0; else if (!strcmp(v, "density")) st->overlay = 1; else if (!strcmp(v, "ke")) st->overlay = 2;
                                     else if (!strcmp(v, "material")) st->overlay = 3; else bad = 1; }
        else if (KEY("--show"))   { char b[256]; snprintf(b, sizeof b, "%s", v); st->show = 0;
                                    for (char* t = strtok(b, ","); t && !bad; t = strtok(NULL, ",")) {
                                        if (!strcmp(t, "disks")) st->show |= PID1_SHOW_DISKS; else if (!strcmp(t, "hull")) st->show |= PID1_SHOW_HULL;
                                        else if (!strcmp(t, "graph")) st->show |= PID1_SHOW_GRAPH; else if (!strcmp(t, "particles")) st->show |= PID1_SHOW_PARTICLES;
                                        else if (!strcmp(t, "food")) st->show |= PID1_SHOW_FOOD; else if (strcmp(t, "none")) bad = 1; } }
        else if (KEY("--seed") || KEY("--grow-seed") || KEY("--N")) {
            char* e; errno = 0; const unsigned long long u = strtoull(v, &e, 0); bad = errno || e == v || *e || v[0] == '-' || u > 9007199254740992ULL;
            if (!bad) { if (KEY("--seed")) st->pseed = (double)u; else if (KEY("--grow-seed")) st->gseed = (double)u; else st->N_cli = (long)u; } }
        else if (!(bad = num_d(v, &x))) {
            if      (KEY("--steps")) st->steps = x;   else if (KEY("--tpf")) st->tpf = x;       else if (KEY("--budget")) st->budget = x;
            else if (KEY("--kappa")) st->kappa = x;   else if (KEY("--tau")) st->tau = x;       else if (KEY("--pbranch")) st->pbranch = x;
            else if (KEY("--beta")) st->beta = x;     else if (KEY("--lambda")) st->lambda = x; else if (KEY("--tadapt")) st->tadapt = x;
            else if (KEY("--gamma")) st->gamma = x;   else if (KEY("--deltathin")) st->deltathin = x;
            else if (KEY("--tick")) st->ptick = x;    else if (KEY("--eta")) st->eta = x;       else if (KEY("--T")) st->T = x;
            else if (KEY("--write-every")) st->wevery = x; else if (KEY("--kbt")) st->kbt = x;
            else if (KEY("--speed")) st->speed = x;   else if (KEY("--alpha")) st->alpha = x;
            else if (KEY("--shot-tick")) st->shot_tick = (long)llround(x); else if (KEY("--shot-at")) st->shot_at = x;
            else { fprintf(stderr, "pid_studio: unknown option '%s'\n", a); usage(stderr); return 2; }
        }
        #undef KEY
        if (bad) { fprintf(stderr, "pid_studio: bad value in '%s'\n", a); return 2; }
    }
    return 0;
}

/* ------------------------------------------------------------------ main */
int main(int argc, char** argv){
    static Studio st;
    st.preset = 3; st.rule = 1; st.gseed = 9300; st.steps = 400; st.tpf = 2; load_preset_sliders(&st);
    st.ptick = -1; st.eta = 0.236; st.pseed = 1; st.T = 200; st.wevery = 10; st.kbt = 1.0;
    st.speed = 0.1; st.alpha = 200.0 / 255.0; st.overlay = 0; st.show = PID1_SHOW_DISKS | PID1_SHOW_PARTICLES | PID1_SHOW_FOOD; st.running = 1;
    st.shot_tick = -1; st.shot_at = -1.0;
    const char* out = NULL; const char* fontp = NULL;
    if (parse(&st, argc, argv, &out, &fontp)) return 2;
    if (st.ptick < 0) st.ptick = st.steps;
    if (out) snprintf(st.root, sizeof st.root, "%s", out);
    else { char ts[32]; const time_t now = time(NULL); strftime(ts, sizeof ts, "%Y%m%dT%H%M%SZ", gmtime(&now)); snprintf(st.root, sizeof st.root, "runs/studio_%s", ts); }
    if (mkdir_p(st.root)) { fprintf(stderr, "pid_studio: cannot create %s\n", st.root); return 2; }

    SDL_SetMainReady();
    if (SDL_Init(SDL_INIT_VIDEO) != 0) { fprintf(stderr, "pid_studio: SDL_Init: %s\n", SDL_GetError()); return 2; }
    if (TTF_Init() != 0) { fprintf(stderr, "pid_studio: TTF_Init: %s\n", TTF_GetError()); return 2; }
    SDL_Window* win = SDL_CreateWindow("PID studio -- grow the network, run the disks", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 1400, 840,
                                       SDL_WINDOW_SHOWN | SDL_WINDOW_RESIZABLE | SDL_WINDOW_ALLOW_HIGHDPI);
    if (!win) { fprintf(stderr, "pid_studio: SDL_CreateWindow: %s\n", SDL_GetError()); return 2; }
    int vsync = 1;
    SDL_Renderer* R = SDL_CreateRenderer(win, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    if (!R) { vsync = 0; R = SDL_CreateRenderer(win, -1, SDL_RENDERER_SOFTWARE); }
    if (!R) { fprintf(stderr, "pid_studio: SDL_CreateRenderer: %s\n", SDL_GetError()); return 2; }
    SDL_SetRenderDrawBlendMode(R, SDL_BLENDMODE_BLEND);
    { int ww, wh, ow, oh; SDL_GetWindowSize(win, &ww, &wh); SDL_GetRendererOutputSize(R, &ow, &oh); ui = ww > 0 ? (float)ow / ww : 1.0f; }
    const char* fonts[] = {fontp, getenv("PID_FONT"), "../../hspist3/Roboto-Regular.ttf", "/System/Library/Fonts/Supplemental/Arial.ttf",
                           "/Library/Fonts/Arial.ttf", "/System/Library/Fonts/Helvetica.ttc", "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"};
    for (size_t i = 0; i < sizeof fonts / sizeof fonts[0] && !font; ++i) if (fonts[i]) { font = TTF_OpenFont(fonts[i], U(14)); font_b = TTF_OpenFont(fonts[i], U(15)); }
    if (!font || !font_b) { fprintf(stderr, "pid_studio: no font found (pass --font=<file.ttf>)\n"); return 2; }

    /* the panels */
    static int view_show_bits;
    Panel P[3]; memset(P, 0, sizeof P);
    static const char* const tog_names[] = {"disks", "hull", "graph", "parts", "food"};
    static const char* const view_btn[] = {"run", "screenshot"};
    P[0] = (Panel){.title = "1  GROW the network (Demo 0)", .open = 1, .has_apply = 1};
    {
        Panel* p = &P[0]; int i = 0;
        p->w[i++] = (Widget){.kind = WK_CYCLE, .label = "preset", .ival = &st.preset, .names = PRESETS, .n = NPRESET};
        p->w[i++] = (Widget){.kind = WK_CYCLE, .label = "rule", .ival = &st.rule, .names = RULES, .n = 2};
        p->w[i++] = (Widget){.kind = WK_SLIDER, .label = "seed", .val = &st.gseed, .lo = 1, .hi = 99999, .logs = 1, .integer = 1};
        p->w[i++] = (Widget){.kind = WK_SLIDER, .label = "ticks", .val = &st.steps, .lo = 1, .hi = 4000, .logs = 1, .integer = 1};
        p->w[i++] = (Widget){.kind = WK_SLIDER, .label = "budget [sigma^2]", .val = &st.budget, .lo = 300, .hi = 30000, .logs = 1};
        p->w[i++] = (Widget){.kind = WK_SLIDER, .label = "kappa (heading noise)", .val = &st.kappa, .lo = 0, .hi = 1.5};
        p->w[i++] = (Widget){.kind = WK_SLIDER, .label = "tau (radius shrink)", .val = &st.tau, .lo = 0, .hi = 0.2};
        p->w[i++] = (Widget){.kind = WK_SLIDER, .label = "pbranch", .val = &st.pbranch, .lo = 0, .hi = 1};
        p->w[i++] = (Widget){.kind = WK_SLIDER, .label = "beta (guidance)", .val = &st.beta, .lo = 0, .hi = 1};
        p->w[i++] = (Widget){.kind = WK_SLIDER, .label = "lambda (cue) [sigma]", .val = &st.lambda, .lo = 2, .hi = 400, .logs = 1};
        p->w[i++] = (Widget){.kind = WK_SLIDER, .label = "Tadapt [ticks] (0 off)", .val = &st.tadapt, .lo = 0, .hi = 50, .integer = 1};
        p->w[i++] = (Widget){.kind = WK_SLIDER, .label = "gamma (thickening)", .val = &st.gamma, .lo = 0, .hi = 0.1};
        p->w[i++] = (Widget){.kind = WK_SLIDER, .label = "deltathin (thinning)", .val = &st.deltathin, .lo = 0, .hi = 0.05};
        p->w[i++] = (Widget){.kind = WK_SLIDER, .label = "ticks per frame (live)", .val = &st.tpf, .lo = 1, .hi = 200, .logs = 1, .integer = 1, .live = 1};
        p->n = i; p->sec_act = A_GROW_PAUSE;
    }
    P[1] = (Panel){.title = "2  PARTICLES in its lumen (Demo 1)", .open = 1, .has_apply = 1};
    {
        Panel* p = &P[1]; int i = 0;
        p->w[i++] = (Widget){.kind = WK_SLIDER, .label = "network at tick", .val = &st.ptick, .lo = 1, .hi = 4000, .logs = 1, .integer = 1};
        p->w[i++] = (Widget){.kind = WK_SLIDER, .label = "packing fraction eta", .val = &st.eta, .lo = 0.005, .hi = 0.45, .logs = 1};
        p->w[i++] = (Widget){.kind = WK_SLIDER, .label = "seed", .val = &st.pseed, .lo = 1, .hi = 9999999, .logs = 1, .integer = 1};
        p->w[i++] = (Widget){.kind = WK_SLIDER, .label = "T [sigma-time]", .val = &st.T, .lo = 10, .hi = 10000, .logs = 1, .integer = 1};
        p->w[i++] = (Widget){.kind = WK_SLIDER, .label = "files every [sigma-time]", .val = &st.wevery, .lo = 1, .hi = 500, .logs = 1, .integer = 1};
        p->w[i++] = (Widget){.kind = WK_SLIDER, .label = "kT", .val = &st.kbt, .lo = 0.1, .hi = 10, .logs = 1};
        p->n = i; p->sec_act = A_PART_STOP;
    }
    P[2] = (Panel){.title = "3  VIEW (acts at once)", .open = 1};
    view_show_bits = 0;
    {
        Panel* p = &P[2]; int i = 0;
        p->w[i++] = (Widget){.kind = WK_TOGGLES, .ival = &view_show_bits, .names = tog_names, .n = 5, .live = 1};
        p->w[i++] = (Widget){.kind = WK_CYCLE, .label = "overlay (particles)", .ival = &st.overlay, .names = OVERLAYS, .n = 4, .live = 1};
        p->w[i++] = (Widget){.kind = WK_SLIDER, .label = "speed [sigma-time/frame]", .val = &st.speed, .lo = 0.002, .hi = 5, .logs = 1, .live = 1};
        p->w[i++] = (Widget){.kind = WK_SLIDER, .label = "disk fill alpha", .val = &st.alpha, .lo = 0, .hi = 1, .live = 1};
        p->w[i++] = (Widget){.kind = WK_BUTTONS, .names = view_btn, .n = 2, .act = {A_RUN_PAUSE, A_SHOT}, .live = 1};
        p->n = i;
    }
    /* the toggles' bits in panel order (disks, hull, graph, particles, food) are the PID1_SHOW_* bits 1, 2, 4, 8, 16 */
    view_show_bits = st.show;

    printf("studio: files go to %s\n", st.root);
    grow_start(&st);                                        /* the window opens with the panel's network growing */
    int quit = 0, drag_p = -1, drag_w = -1, scroll = 0;
    Widget* drag = NULL;
    SDL_Texture* tex = NULL;
    while (!quit) {
        int ow, oh, ww, wh; SDL_GetRendererOutputSize(R, &ow, &oh); SDL_GetWindowSize(win, &ww, &wh);
        const float mx_scale = ww > 0 ? (float)ow / ww : 1.0f;
        const int PW = U(390), px0 = ow - PW;
        /* layout */
        int y = U(8) - scroll;
        for (int k = 0; k < 3; ++k) y = panel_layout(&P[k], px0 + U(8), y, PW - U(16));
        const int panels_h = y + scroll;
        int shot = 0; const char* shot_tag = "manual"; char tagbuf[64];
        SDL_Event e;
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT) quit = 1;
            else if (e.type == SDL_KEYDOWN && !e.key.repeat) {
                switch (e.key.keysym.sym) {
                case SDLK_ESCAPE: quit = 1; break;
                case SDLK_SPACE: if (st.have_run) st.running = !st.running; else if (st.have_w) st.growing = !st.growing; break;
                case SDLK_d: view_show_bits ^= PID1_SHOW_DISKS; break;
                case SDLK_h: view_show_bits ^= PID1_SHOW_HULL; break;
                case SDLK_g: view_show_bits ^= PID1_SHOW_GRAPH; break;
                case SDLK_p: view_show_bits ^= PID1_SHOW_PARTICLES; break;
                case SDLK_f: view_show_bits ^= PID1_SHOW_FOOD; break;
                case SDLK_c: st.overlay = (st.overlay + 1) % 4; break;
                case SDLK_LEFTBRACKET:  st.alpha = st.alpha > 0.0625 ? st.alpha - 0.0625 : 0.0; break;
                case SDLK_RIGHTBRACKET: st.alpha = st.alpha < 0.9375 ? st.alpha + 0.0625 : 1.0; break;
                case SDLK_s: shot = 1; break;
                default: break;
                }
            } else if (e.type == SDL_MOUSEBUTTONDOWN) {
                const int mx = (int)(e.button.x * mx_scale), my = (int)(e.button.y * mx_scale);
                for (int k = 0; k < 3; ++k) {
                    Panel* p = &P[k];
                    if (inside(p->title_r, mx, my)) { p->open = !p->open; break; }
                    if (!p->open) continue;
                    if (p->has_apply && inside(p->apply_r, mx, my)) {
                        if (k == 0) { run_end(&st); if (grow_start(&st)) p->dirty = 0; }
                        else if (particles_start(&st)) { p->dirty = 0; P[0].open = 0; }   /* the growth panel folds away during a run */
                        break;
                    }
                    if (p->has_apply && p->sec_act && inside(p->sec_r, mx, my)) {
                        if (p->sec_act == A_GROW_PAUSE && st.have_w && !st.have_run) st.growing = !st.growing;
                        if (p->sec_act == A_PART_STOP && st.have_run) { run_end(&st); say(&st, "run stopped (its manifest is written)"); }
                        break;
                    }
                    for (int i = 0; i < p->n; ++i) {
                        Widget* w = &p->w[i];
                        if (w->kind == WK_SLIDER && inside((SDL_Rect){w->track.x - U(6), w->track.y - U(8), w->track.w + U(12), w->track.h + U(16)}, mx, my)) {
                            drag = w; drag_p = k; drag_w = i; slider_set(w, (double)(mx - w->track.x) / w->track.w); if (!w->live) p->dirty = 1;
                        } else if (w->kind == WK_CYCLE && inside(w->sub[0], mx, my)) {
                            const int dir = (e.button.button == SDL_BUTTON_RIGHT || mx < w->sub[0].x + w->sub[0].w / 3) ? -1 : 1;
                            *w->ival = (*w->ival + dir + w->n) % w->n;
                            if (w->ival == &st.preset) load_preset_sliders(&st);
                            if (!w->live) p->dirty = 1;
                        } else if (w->kind == WK_TOGGLES || w->kind == WK_BUTTONS) {
                            for (int j = 0; j < w->n; ++j) {
                                if (!inside(w->sub[j], mx, my)) continue;
                                if (w->kind == WK_TOGGLES) { *w->ival ^= 1 << j; continue; }
                                const int act = w->act[j];
                                if (act == A_GROW_PAUSE && st.have_w) st.growing = !st.growing;
                                if (act == A_PART_STOP && st.have_run) { run_end(&st); say(&st, "run stopped (its manifest is written)"); }
                                if (act == A_RUN_PAUSE) st.running = !st.running;
                                if (act == A_SHOT) shot = 1;
                            }
                        }
                    }
                }
            } else if (e.type == SDL_MOUSEBUTTONUP) { drag = NULL; drag_p = drag_w = -1; }
            else if (e.type == SDL_MOUSEMOTION && drag) {
                const int mx = (int)(e.motion.x * mx_scale);
                slider_set(drag, (double)(mx - drag->track.x) / drag->track.w);
                if (!drag->live && drag_p >= 0) P[drag_p].dirty = 1;
            } else if (e.type == SDL_MOUSEWHEEL) {
                int mx, my; SDL_GetMouseState(&mx, &my); mx = (int)(mx * mx_scale); my = (int)(my * mx_scale);
                int used = 0;
                if (mx >= px0)
                    for (int k = 0; k < 3 && !used; ++k) for (int i = 0; i < P[k].n && !used; ++i) {
                        Widget* w = &P[k].w[i];
                        if (P[k].open && w->kind == WK_SLIDER && inside(w->r, mx, my)) { slider_nudge(w, e.wheel.y > 0 ? 1 : -1); if (!w->live) P[k].dirty = 1; used = 1; }
                    }
                if (!used && mx >= px0) { scroll -= e.wheel.y * U(40); if (scroll < 0) scroll = 0; if (scroll > panels_h - oh + U(40) && panels_h > oh) scroll = panels_h - oh + U(40); if (panels_h <= oh) scroll = 0; }
            }
        }
        (void)drag_w;
        /* the growth: tpf ticks per frame up to --steps (stopping at a --shot-tick) */
        if (st.have_w && st.growing && !st.have_run) {
            const long K = (long)llround(st.steps);
            long n = (long)llround(st.tpf);
            if (st.shot_tick >= 0 && !st.shot_tick_done && st.w.tick < st.shot_tick && st.w.tick + n > st.shot_tick) n = st.shot_tick - st.w.tick;
            const Policy pol = policy_by_name(st.grule);
            for (long k = 0; k < n && st.w.tick < K; ++k) step(&st.w, pol);
            if (st.w.tick >= K) { st.growing = 0; say(&st, "grown: tick %ld, %d disks alive", st.w.tick, world_alive_disks(&st.w)); }
            if (st.shot_tick >= 0 && !st.shot_tick_done && st.w.tick == st.shot_tick) { shot = 1; snprintf(tagbuf, sizeof tagbuf, "grow_tick%06ld", st.w.tick); shot_tag = tagbuf; st.shot_tick_done = 1; }
        }
        /* --auto: after the growth, APPLY the particles once (or, with --grow-only, write the network and stop) */
        if (st.auto_mode == 1 && st.have_w && !st.growing && st.w.tick >= (long)llround(st.steps)) {
            st.auto_mode = 2;
            if (st.grow_only) { char path[1200]; st.ptick = st.steps; grow_hull(&st); write_net(&st, path, sizeof path); printf("studio: network written to %s\n", path); if (st.exit_at_end) quit = 1; }
            else if (particles_start(&st)) P[0].open = 0;
            else if (st.exit_at_end) quit = 1;
        }
        /* the particles */
        if (st.have_run && st.running && !st.run.done) {
            const double now = run1_time(&st.run);
            double tt = now + st.speed; int at = 0;
            if (st.shot_at >= 0.0 && !st.shot_at_done && st.shot_at > now + 1e-9 && st.shot_at <= tt) { tt = st.shot_at; at = 1; }
            run1_advance(&st.run, tt);
            if (at && run1_time(&st.run) >= st.shot_at - 1e-9) { shot = 1; snprintf(tagbuf, sizeof tagbuf, "run%ld_t%06.0f", st.run_n, st.shot_at); shot_tag = tagbuf; st.shot_at_done = 1; }
        }
        if (st.have_run && st.run.done && !st.reported) {
            const int rc = run1_report(&st.run); st.reported = 1;
            say(&st, "run %ld reached T = %.0f: %s (summary and checks in the terminal; files in %s)", st.run_n, st.run.T, rc ? "a check FAILED" : "all checks OK", st.run.out);
        }
        if (st.exit_at_end && st.auto_mode == 2 && (!st.have_run || st.reported) && (st.shot_at < 0.0 || st.shot_at_done || !st.have_run || st.run.done)) quit = 1;

        /* draw: the scene left of the panels */
        SDL_SetRenderDrawColor(R, 14, 16, 22, 255); SDL_RenderClear(R);
        const int hud_h = U(56), m = U(16);
        const double aw = px0 - 2.0 * m, ah = oh - hud_h - 2.0 * m;
        const World* w = st.have_run ? &st.run.w : &st.w;
        Pid1Vp vp; vp.H = 120.0; vp.scale = fmin(aw / 200.0, ah / 120.0); vp.ox = m + (aw - 200.0 * vp.scale) / 2.0; vp.oy = hud_h + m;
        if (st.have_w || st.have_run) {
            if (!st.have_run && (view_show_bits & PID1_SHOW_HULL) && st.w.nd > 1) grow_hull(&st);
            Pid1Scene sc = {view_show_bits, st.have_run ? st.overlay : 0, (int)lround(st.alpha * 255.0), !st.have_run};
            vp.H = w->p.H; vp.scale = fmin(aw / w->p.W, ah / w->p.H); vp.ox = m + (aw - w->p.W * vp.scale) / 2.0;
            pd_scene(R, &vp, w, st.have_run ? &st.run : NULL, &sc, &tex);
        }
        /* HUD */
        {
            const SDL_Color tc = {225, 228, 236, 255}, hc = {255, 214, 120, 255};
            char l1[512], l2[512];
            if (st.have_run) {
                const EDMD3_Health* H = edmd3_health(st.run.S);
                const double drift = (edmd3_kinetic_energy(st.run.S) - st.run.ke0) / st.run.ke0;
                snprintf(l1, sizeof l1, "PARTICLES %ld   t %.2f / %.0f sigma-time   N %ld  eta %.3f   network tick %ld   %s", st.run_n, run1_time(&st.run), st.run.T, st.run.N,
                         st.run.eta, st.run.w.tick, st.run.done ? "(reached T)" : st.running ? "" : "(paused)");
                snprintf(l2, sizeof l2, "events: pair %ld  wall %ld  arc %ld   KE drift %.1e   clean %d   escaped %ld   overlay %s   files %s", H->ev_pair, H->ev_seg,
                         H->ev_arc, drift, edmd3_health_clean(st.run.S), st.run.esc_total, OVERLAYS[st.overlay], st.run.out);
            } else if (st.have_w) {
                snprintf(l1, sizeof l1, "GROW %ld   tick %ld / %.0f   preset %s  rule %s  seed %" PRIu64 "   %s", st.grow_n, st.w.tick, st.steps, PRESETS[st.preset],
                         st.grule, st.gs, st.growing ? "" : "(stopped)");
                snprintf(l2, sizeof l2, "disks %d alive   live tips %d   M_used %.1f   M_pool %.1f [sigma^2]   %s", world_alive_disks(&st.w), world_live_tips(&st.w),
                         st.w.m_used, world_pool(&st.w), st.w.hull.valid && st.hull_tick == st.w.tick ? "hull shown" : "");
            } else { snprintf(l1, sizeof l1, "no network yet: panel 1, APPLY"); l2[0] = '\0'; }
            text(R, font_b, m, U(8), hc, l1);
            text(R, font, m, U(30), tc, l2);
            if (st.msg[0] && SDL_GetTicks() - st.msg_ms < 8000) {
                const int tw = text_w(font, st.msg);
                fill(R, (SDL_Rect){m - U(4), oh - U(30), tw + U(8), U(24)}, 40, 44, 58, 230);
                text(R, font, m, oh - U(27), (SDL_Color){255, 236, 200, 255}, st.msg);
            }
        }
        /* the panels */
        fill(R, (SDL_Rect){px0, 0, PW, oh}, 22, 24, 31, 255);
        for (int k = 0; k < 3; ++k) panel_draw(R, &P[k], &st);
        if (shot) save_shot(R, &st, shot_tag);
        SDL_RenderPresent(R);
        if (!vsync && !st.exit_at_end) SDL_Delay(15);
    }
    run_end(&st);
    if (st.have_w) world_free(&st.w);
    if (tex) SDL_DestroyTexture(tex);
    TTF_CloseFont(font); TTF_CloseFont(font_b); TTF_Quit();
    SDL_DestroyRenderer(R); SDL_DestroyWindow(win); SDL_Quit();
    return 0;
}
