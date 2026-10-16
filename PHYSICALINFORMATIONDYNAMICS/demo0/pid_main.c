/* pid_main.c -- command line, the run loop shared by the headless path and the window, and the run manifest.
   pid_demo0 --seed=S --rule=unbiased --steps=N --out=runs/<name> [--headless] [--write-every=k] [--font=path] [parameters]
   Nothing is written anywhere but --out. Flags are matched exactly; an unknown flag is an error. */
#include <dirent.h>
#include <errno.h>
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include "pid.h"

#ifndef PID_BUILD_ID
#define PID_BUILD_ID "unknown"
#endif
#ifndef PID_CFLAGS
#define PID_CFLAGS "unknown"
#endif
#ifndef PID_FOLDER_COMMIT
#define PID_FOLDER_COMMIT "unknown"
#endif
#ifndef PID_SRC_SHA256
#define PID_SRC_SHA256 "unknown"
#endif

const char* pid_build_id(void){ return PID_BUILD_ID " " PID_CFLAGS; }

static void usage(FILE* f){
    fprintf(f,
        "usage: pid_demo0 --seed=S --steps=N --out=DIR [options] [parameters]\n"
        "  --rule=unbiased|guided growth rule; guided bends each step toward the cue of the food patches (0b)\n"
        "  --headless             no window; same files, then exit\n"
        "  --write-every=k        state_<tick>.txt every k ticks (0: only final.txt); default 50\n"
        "  --font=PATH            status line on screen (SDL2_ttf build); otherwise in the terminal\n"
        "  --tpf=n                window: ticks per frame, default 1\n"
        "  --tick-ms=m            window: at most one batch of --tpf ticks every m ms, default 50 (0: every frame)\n"
        "  --exit-at-end          window: close when --steps is reached (tests)\n"
        "  --shot-at-end          window: write shot_<tick>.bmp into DIR when --steps is reached\n"
        "  --overwrite            allow a non-empty DIR\n"
        "  --preset=NAME          a saved parameter set, applied first; explicit flags override it:\n"
        "                           mat  the dense branching mat (Chris, 2026-10-15): --pbranch=1 --pbranch-exp=0 --delta=1.57\n"
        "                                --ntry=48 --budget=9000, the rest as the defaults\n"
        "                           0a1  the defaults (Demo 0a.1)\n"
        "                           0a   Demo 0a: --start-rings=0 --pbranch-exp=0 --pbranch=0.08 --ntry=8 --budget=3000\n"
        "                           mat45 Demo 0b (sec. 8): mat at --budget=4500, patch A under the seed, patch B at (0.8W, 0.5H),\n"
        "                                --tadapt=10 --budget-wait=1 --tip-order=random; run it with --rule=guided\n"
        "  --version, --help\n"
        "parameters (defaults in brackets): --W [200] --H [120] --budget [4500] --kappa [0.35] --tau [0.03] --rmin [1]\n"
        "  --pbranch [0.8] --delta [0.6] --ntry [24] --K [6] --R0 [3]\n"
        "  --start-rings [1]      hex rings of R0 disks around the seed; their outer ring carries the tips (0: K tips on the seed)\n"
        "  --pbranch-exp [6]      p_branch(r) = pbranch (r / R0)^pbranch-exp for a tip on a disk of radius r (0: constant)\n"
        "  the 0a defaults: --start-rings=0 --pbranch-exp=0 --pbranch=0.08 --ntry=8 --budget=3000\n"
        "  0b: --patch=x,y,r,s (repeatable, sigma; replaces the preset's patches) --beta [0.3] --lambda [40] --tadapt [0, off]\n"
        "  --gamma [0.02] --deltathin [0.01] --rmax [6] --budget-wait=0|1 [0] --tip-order=id|random [id]\n"
        "window keys: D disks, G graph, H hull (0c), C cue (0b), P patches (0b), [ ] disk transparency, Space pause/run,\n"
        "  N one tick, R reset to the same seed, S BMP into DIR, Esc quit\n");
}

/* ---------------------------------------------------------------- presets (0a.2): saved parameter sets */
static int apply_preset(const char* name, Params* p){
    params_default(p);
    if (!strcmp(name, "0a1")) return 0;
    if (!strcmp(name, "0a"))  { p->rings = 0; p->pbranch_exp = 0.0; p->pbranch = 0.08; p->ntry = 8; p->budget = 3000.0; return 0; }
    if (!strcmp(name, "mat")) { p->pbranch = 1.0; p->pbranch_exp = 0.0; p->delta = 1.57; p->ntry = 48; p->budget = 9000.0; return 0; }
    if (!strcmp(name, "mat45")) {   /* Demo 0b (sec. 8 as relayed 2026-10-15) */
        apply_preset("mat", p); p->budget = 4500.0;
        p->tadapt = 10; p->budget_wait = 1; p->tip_order = 1;
        p->npatch = 2;
        p->patch[0] = (Patch){0.5 * p->W, 0.5 * p->H, 5.0, 1.0, 1};      /* A, under the seed */
        p->patch[1] = (Patch){0.8 * p->W, 0.5 * p->H, 5.0, 1.0, 1};      /* B */
        return 0;
    }
    return -1;
}

/* ---------------------------------------------------------------- strict number parsing */
static int num_d(const char* v, double* out){ char* e; errno = 0; const double x = strtod(v, &e); if (e == v || *e || errno) return -1; *out = x; return 0; }
static int num_l(const char* v, long* out){ char* e; errno = 0; const long x = strtol(v, &e, 10); if (e == v || *e || errno) return -1; *out = x; return 0; }
static int num_u64(const char* v, uint64_t* out){
    if (v[0] == '-' || v[0] == '\0') return -1;
    char* e; errno = 0; const unsigned long long x = strtoull(v, &e, 10); if (*e || errno) return -1; *out = (uint64_t)x; return 0;
}

/* ---------------------------------------------------------------- run folder */
static int mkdir_p(const char* path){
    char buf[1024]; snprintf(buf, sizeof buf, "%s", path);
    for (char* p = buf + 1; *p; ++p) if (*p == '/') { *p = '\0'; if (mkdir(buf, 0755) && errno != EEXIST) return -1; *p = '/'; }
    if (mkdir(buf, 0755) && errno != EEXIST) return -1;
    return 0;
}
static int dir_empty(const char* path){
    DIR* d = opendir(path); if (!d) return 1;
    struct dirent* e; int empty = 1;
    while ((e = readdir(d))) if (strcmp(e->d_name, ".") && strcmp(e->d_name, "..")) { empty = 0; break; }
    closedir(d); return empty;
}

static void note_file(Run* r, const char* name){
    for (int i = 0; i < r->nfiles; ++i) if (!strcmp(r->files[i], name)) return;
    if (r->nfiles == r->cap_files) {
        r->cap_files = r->cap_files ? 2 * r->cap_files : 64;
        r->files = realloc(r->files, (size_t)r->cap_files * sizeof *r->files);
        if (!r->files) { fprintf(stderr, "pid: out of memory (file list)\n"); exit(3); }
    }
    snprintf(r->files[r->nfiles++], 64, "%s", name);
}

/* Every written state is checked against the invariants of brief sec. 6 items 2-3 before the run goes on. */
static void write_state(Run* r, const char* name){
    char path[1200], err[512];
    snprintf(path, sizeof path, "%s/%s", r->out, name);
    if (state_write(&r->w, path)) { fprintf(stderr, "pid: cannot write %s\n", path); exit(2); }
    note_file(r, name);
    const int fails = world_check(&r->w, r->w.p.tadapt == 0, world_pool(&r->w), err, sizeof err);   /* tangency until adaptation runs */
    if (fails) { fprintf(stderr, "pid: INVARIANT FAILED at tick %ld (%d failures): %s\n", r->w.tick, fails, err); exit(4); }
}

/* ---------------------------------------------------------------- the run loop (headless and window) */
int run_tick(Run* r){
    if (r->w.tick >= r->steps) return 0;
    step(&r->w, r->policy);
    if (r->write_every > 0 && r->w.tick % r->write_every == 0) {
        char name[64]; snprintf(name, sizeof name, "state_%06ld.txt", r->w.tick); write_state(r, name);
    }
    if (r->w.tick == r->steps) { write_state(r, "final.txt"); r->final_written = 1; run_finish(r, 1); }
    return 1;
}
static void run_start(Run* r){
    world_init(&r->w, &r->p, r->seed, r->rule);
    if (r->write_every > 0) write_state(r, "state_000000.txt");
    if (r->steps == 0) { write_state(r, "final.txt"); r->final_written = 1; run_finish(r, 1); }
}
void run_reset(Run* r){ world_free(&r->w); run_start(r); }

void run_status(const Run* r, FILE* f){
    fprintf(f, "tick %ld | live tips %d | disks %d | M_used %.3f | M_pool %.3f\n",
            r->w.tick, world_live_tips(&r->w), world_alive_disks(&r->w), r->w.m_used, world_pool(&r->w));
}

static void fprint_arg(FILE* f, const char* a){
    int plain = a[0] != '\0';
    for (const char* c = a; *c; ++c) if (!strchr("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_./=:,+-", *c)) plain = 0;
    if (plain) { fputs(a, f); return; }
    fputc('\'', f);
    for (const char* c = a; *c; ++c) { if (*c == '\'') fputs("'\\''", f); else fputc(*c, f); }
    fputc('\'', f);
}

int run_finish(Run* r, int completed){
    char path[1200]; snprintf(path, sizeof path, "%s/manifest.txt", r->out);
    FILE* f = fopen(path, "w"); if (!f) { fprintf(stderr, "pid: cannot write %s\n", path); return -1; }
    char pl[1024]; params_line(&r->p, pl, sizeof pl);
    fprintf(f, "# pid_demo0 run manifest (brief sec. 4)\n");
    fprintf(f, "argv"); for (int i = 0; i < r->argc; ++i) { fputc(' ', f); fprint_arg(f, r->argv[i]); } fputc('\n', f);
    fprintf(f, "seed %" PRIu64 "\n", r->seed);
    fprintf(f, "rule %s\n", r->rule);
    fprintf(f, "steps %ld\n", r->steps);
    fprintf(f, "write_every %ld\n", r->write_every);
    fprintf(f, "box %.17g %.17g\n", r->p.W, r->p.H);
    fprintf(f, "budget %.17g\n", r->p.budget);
    fprintf(f, "params %s\n", pl);
    fprintf(f, "build %s\n", pid_build_id());
    fprintf(f, "folder_commit %s   # last commit of PHYSICALINFORMATIONDYNAMICS/ when the binary was built\n", PID_FOLDER_COMMIT);
    fprintf(f, "src_sha256 %s   # sha256 of the concatenated sources in the Makefile's SRC_LIST order\n", PID_SRC_SHA256);
    fprintf(f, "preset %s\n", r->preset ? r->preset : "none");
    fprintf(f, "mode %s\n", r->headless ? "headless" : "window");
    fprintf(f, "completed %d   # 1: final.txt written at tick = steps\n", completed);
    fprintf(f, "tick %ld\n", r->w.tick);
    fprintf(f, "# sha256 of every state file (check inside the run folder: shasum -a 256 -c, after removing these comment lines)\n");
    for (int i = 0; i < r->nfiles; ++i) {
        char fp[1200], hex[65]; snprintf(fp, sizeof fp, "%s/%s", r->out, r->files[i]);
        if (sha256_file(fp, hex)) snprintf(hex, sizeof hex, "unreadable");
        fprintf(f, "%s  %s\n", hex, r->files[i]);
    }
    return fclose(f) == 0 ? 0 : -1;
}

/* ---------------------------------------------------------------- main */
int main(int argc, char** argv){
    Run r; memset(&r, 0, sizeof r);
    params_default(&r.p);
    snprintf(r.rule, sizeof r.rule, "unbiased");
    r.write_every = 50; r.tpf = 1; r.tick_ms = 50; r.steps = -1; r.argc = argc; r.argv = argv;
    int have_seed = 0, overwrite = 0, patches_given = 0;
    for (int i = 1; i < argc; ++i)                      /* the preset first, so that explicit flags override it */
        if (!strncmp(argv[i], "--preset=", 9)) {
            if (r.preset) { fprintf(stderr, "pid: --preset given twice\n"); return 2; }
            r.preset = argv[i] + 9;
            if (apply_preset(r.preset, &r.p)) { fprintf(stderr, "pid: unknown preset '%s' (mat45, mat, 0a1, 0a)\n", r.preset); return 2; }
        }
    for (int i = 1; i < argc; ++i) {
        const char* a = argv[i];
        const char* eq = strchr(a, '=');
        char key[64]; const size_t kl = eq ? (size_t)(eq - a) : strlen(a);
        if (kl >= sizeof key || strncmp(a, "--", 2)) { fprintf(stderr, "pid: bad argument '%s'\n", a); usage(stderr); return 2; }
        memcpy(key, a, kl); key[kl] = '\0';
        const char* v = eq ? eq + 1 : NULL;
        long lv = 0; int bad = 0;
        #define NEEDV() do { if (!v) { fprintf(stderr, "pid: %s needs =value\n", key); return 2; } } while (0)
        #define NOV()   do { if (v)  { fprintf(stderr, "pid: %s takes no value\n", key); return 2; } } while (0)
        if      (!strcmp(key, "--help"))        { NOV(); usage(stdout); return 0; }
        else if (!strcmp(key, "--version"))     { NOV(); printf("pid_demo0 build %s\nfolder_commit %s\nsrc_sha256 %s\n", pid_build_id(), PID_FOLDER_COMMIT, PID_SRC_SHA256); return 0; }
        else if (!strcmp(key, "--seed"))        { NEEDV(); bad = num_u64(v, &r.seed); have_seed = 1; }
        else if (!strcmp(key, "--rule"))        { NEEDV(); snprintf(r.rule, sizeof r.rule, "%s", v); }
        else if (!strcmp(key, "--steps"))       { NEEDV(); bad = num_l(v, &r.steps) || r.steps < 0; }
        else if (!strcmp(key, "--out"))         { NEEDV(); snprintf(r.out, sizeof r.out, "%s", v); bad = v[0] == '\0'; }
        else if (!strcmp(key, "--preset"))      { NEEDV(); }  /* applied above */
        else if (!strcmp(key, "--headless"))    { NOV(); r.headless = 1; }
        else if (!strcmp(key, "--write-every")) { NEEDV(); bad = num_l(v, &r.write_every) || r.write_every < 0; }
        else if (!strcmp(key, "--font"))        { NEEDV(); r.font = v; }
        else if (!strcmp(key, "--tpf"))         { NEEDV(); bad = num_l(v, &lv) || lv < 1 || lv > 1000000; r.tpf = (int)lv; }
        else if (!strcmp(key, "--tick-ms"))     { NEEDV(); bad = num_l(v, &lv) || lv < 0 || lv > 60000; r.tick_ms = (int)lv; }
        else if (!strcmp(key, "--exit-at-end")) { NOV(); r.exit_at_end = 1; }
        else if (!strcmp(key, "--shot-at-end")) { NOV(); r.shot_at_end = 1; }
        else if (!strcmp(key, "--overwrite"))   { NOV(); overwrite = 1; }
        else if (!strcmp(key, "--W"))           { NEEDV(); bad = num_d(v, &r.p.W); }
        else if (!strcmp(key, "--H"))           { NEEDV(); bad = num_d(v, &r.p.H); }
        else if (!strcmp(key, "--budget"))      { NEEDV(); bad = num_d(v, &r.p.budget); }
        else if (!strcmp(key, "--kappa"))       { NEEDV(); bad = num_d(v, &r.p.kappa); }
        else if (!strcmp(key, "--tau"))         { NEEDV(); bad = num_d(v, &r.p.tau); }
        else if (!strcmp(key, "--rmin"))        { NEEDV(); bad = num_d(v, &r.p.rmin); }
        else if (!strcmp(key, "--pbranch"))     { NEEDV(); bad = num_d(v, &r.p.pbranch); }
        else if (!strcmp(key, "--delta"))       { NEEDV(); bad = num_d(v, &r.p.delta); }
        else if (!strcmp(key, "--R0"))          { NEEDV(); bad = num_d(v, &r.p.R0); }
        else if (!strcmp(key, "--ntry"))        { NEEDV(); bad = num_l(v, &lv) || lv < 0 || lv > 1000000; r.p.ntry = (int)lv; }
        else if (!strcmp(key, "--K"))           { NEEDV(); bad = num_l(v, &lv) || lv < 1 || lv > 1000000; r.p.K = (int)lv; }
        else if (!strcmp(key, "--start-rings")) { NEEDV(); bad = num_l(v, &lv) || lv < 0 || lv > 50; r.p.rings = (int)lv; }
        else if (!strcmp(key, "--beta"))        { NEEDV(); bad = num_d(v, &r.p.beta); }
        else if (!strcmp(key, "--lambda"))      { NEEDV(); bad = num_d(v, &r.p.lambda); }
        else if (!strcmp(key, "--tadapt"))      { NEEDV(); bad = num_l(v, &lv) || lv < 0 || lv > 1000000; r.p.tadapt = (int)lv; }
        else if (!strcmp(key, "--gamma"))       { NEEDV(); bad = num_d(v, &r.p.gamma); }
        else if (!strcmp(key, "--deltathin"))   { NEEDV(); bad = num_d(v, &r.p.deltathin); }
        else if (!strcmp(key, "--rmax"))        { NEEDV(); bad = num_d(v, &r.p.rmax); }
        else if (!strcmp(key, "--budget-wait")) { NEEDV(); bad = num_l(v, &lv) || (lv != 0 && lv != 1); r.p.budget_wait = (int)lv; }
        else if (!strcmp(key, "--tip-order"))   { NEEDV(); if (!strcmp(v, "id")) r.p.tip_order = 0; else if (!strcmp(v, "random")) r.p.tip_order = 1; else bad = 1; }
        else if (!strcmp(key, "--patch"))       { NEEDV(); Patch q = {0, 0, 0, 0, 1}; char tail;
                                                  if (!patches_given) { r.p.npatch = 0; patches_given = 1; }
                                                  bad = r.p.npatch >= PID_MAX_PATCHES || sscanf(v, "%lf,%lf,%lf,%lf%c", &q.x, &q.y, &q.r, &q.s, &tail) != 4;
                                                  if (!bad) r.p.patch[r.p.npatch++] = q; }
        else if (!strcmp(key, "--pbranch-exp")) { NEEDV(); bad = num_d(v, &r.p.pbranch_exp); }
        else { fprintf(stderr, "pid: unknown option '%s'\n", key); usage(stderr); return 2; }
        if (bad) { fprintf(stderr, "pid: bad value in '%s'\n", a); return 2; }
    }
    if (!have_seed || r.steps < 0 || !r.out[0]) { fprintf(stderr, "pid: --seed, --steps and --out are required\n"); usage(stderr); return 2; }
    r.policy = policy_by_name(r.rule);
    if (!r.policy) { fprintf(stderr, "pid: --rule=%s is not a rule (unbiased, guided)\n", r.rule); return 2; }
    char err[256];
    if (params_check(&r.p, err, sizeof err)) { fprintf(stderr, "pid: %s\n", err); return 2; }
    if (mkdir_p(r.out)) { fprintf(stderr, "pid: cannot create %s\n", r.out); return 2; }
    if (!overwrite && !dir_empty(r.out)) { fprintf(stderr, "pid: %s is not empty; choose a new --out or pass --overwrite\n", r.out); return 2; }

    run_start(&r);
    int rc = 0;
    if (r.headless) { while (run_tick(&r)) {} run_status(&r, stdout); }
    else            { rc = view_run(&r); if (!r.final_written) run_finish(&r, 0); }
    world_free(&r.w); free(r.files);
    return rc;
}
