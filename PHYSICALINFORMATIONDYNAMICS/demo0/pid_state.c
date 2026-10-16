/* pid_state.c -- the state file of brief sec. 4: plain text, one record per line, '#' comments, sections in order.
   Every double is printed with %.17g, so read -> write is byte-identical (tests/run_tests.sh, round trip).
   Demo 0a writes [header], [disks], [tips]; [patches] (0b), [hull] and [areas] (0c) come later. */
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include "pid.h"

static const char* KIND_NAME[3] = {"seed", "body", "tip"};

void params_line(const Params* p, char* buf, size_t n){
    snprintf(buf, n, "kappa=%.17g tau=%.17g rmin=%.17g pbranch=%.17g delta=%.17g ntry=%d K=%d R0=%.17g rings=%d pbranch_exp=%.17g"
             " beta=%.17g lambda=%.17g Tadapt=%d gamma=%.17g deltathin=%.17g rmax=%.17g budget_wait=%d tip_order=%s",
             p->kappa, p->tau, p->rmin, p->pbranch, p->delta, p->ntry, p->K, p->R0, p->rings, p->pbranch_exp,
             p->beta, p->lambda, p->tadapt, p->gamma, p->deltathin, p->rmax, p->budget_wait, p->tip_order ? "random" : "id");
}

int state_write_fp(const World* w, FILE* f){
    char pl[1024]; params_line(&w->p, pl, sizeof pl);
    fprintf(f, "[header]\n");
    fprintf(f, "format %d\n", PID_FORMAT);
    fprintf(f, "unit sigma\n");
    fprintf(f, "geometry %s\n", PID_GEOMETRY);
    fprintf(f, "seed %" PRIu64 "\n", w->seed);
    fprintf(f, "rule %s\n", w->rule);
    fprintf(f, "tick %ld\n", w->tick);
    fprintf(f, "box %.17g %.17g\n", w->p.W, w->p.H);
    fprintf(f, "budget %.17g\n", w->p.budget);
    fprintf(f, "m_used %.17g\n", w->m_used);
    fprintf(f, "m_pool %.17g\n", world_pool(w));
    fprintf(f, "params %s\n", pl);
    fprintf(f, "build %s\n", w->build);
    fprintf(f, "[disks]\n# id x y r parent kind born alive\n");
    for (int i = 0; i < w->nd; ++i) {
        const Disk* d = &w->d[i];
        fprintf(f, "%d %.17g %.17g %.17g %d %s %ld %d\n", d->id, d->x, d->y, d->r, d->parent, KIND_NAME[d->kind], d->born, d->alive);
    }
    fprintf(f, "[tips]\n# disk heading alive\n");
    for (int i = 0; i < w->nt; ++i) fprintf(f, "%d %.17g %d\n", w->t[i].disk, w->t[i].heading, w->t[i].alive);
    if (w->p.npatch > 0) {                             /* 0b; absent when there are no patches (0a, 0a.1) */
        fprintf(f, "[patches]\n# id x y r s active\n");
        for (int k = 0; k < w->p.npatch; ++k) {
            const Patch* q = &w->p.patch[k];
            fprintf(f, "%d %.17g %.17g %.17g %.17g %d\n", k, q->x, q->y, q->r, q->s, q->active);
        }
    }
    return ferror(f) ? -1 : 0;
}

int state_write(const World* w, const char* path){
    FILE* f = fopen(path, "w"); if (!f) return -1;
    const int rc = state_write_fp(w, f);
    return (fclose(f) == 0 && rc == 0) ? 0 : -1;
}

/* ---------------------------------------------------------------- reader */
#define FAIL(...) do { snprintf(err, n, __VA_ARGS__); fclose(f); world_free(w); return -1; } while (0)

static int parse_params(Params* p, const char* s){
    char buf[1024]; snprintf(buf, sizeof buf, "%s", s);
    int seen = 0;
    for (char* tok = strtok(buf, " "); tok; tok = strtok(NULL, " ")) {
        char* eq = strchr(tok, '='); if (!eq) return -1;
        *eq = '\0'; const char* k = tok; const char* v = eq + 1; char* end;
        if      (!strcmp(k, "kappa"))   { p->kappa   = strtod(v, &end); }
        else if (!strcmp(k, "tau"))     { p->tau     = strtod(v, &end); }
        else if (!strcmp(k, "rmin"))    { p->rmin    = strtod(v, &end); }
        else if (!strcmp(k, "pbranch")) { p->pbranch = strtod(v, &end); }
        else if (!strcmp(k, "delta"))   { p->delta   = strtod(v, &end); }
        else if (!strcmp(k, "R0"))      { p->R0      = strtod(v, &end); }
        else if (!strcmp(k, "ntry"))    { p->ntry    = (int)strtol(v, &end, 10); }
        else if (!strcmp(k, "K"))       { p->K       = (int)strtol(v, &end, 10); }
        else if (!strcmp(k, "rings"))       { p->rings = (int)strtol(v, &end, 10); --seen; }        /* 0a.1 keys; absent in 0a files */
        else if (!strcmp(k, "pbranch_exp")) { p->pbranch_exp = strtod(v, &end); --seen; }
        else if (!strcmp(k, "beta"))        { p->beta = strtod(v, &end); --seen; }             /* 0b keys; absent before 0b */
        else if (!strcmp(k, "lambda"))      { p->lambda = strtod(v, &end); --seen; }
        else if (!strcmp(k, "Tadapt"))      { p->tadapt = (int)strtol(v, &end, 10); --seen; }
        else if (!strcmp(k, "gamma"))       { p->gamma = strtod(v, &end); --seen; }
        else if (!strcmp(k, "deltathin"))   { p->deltathin = strtod(v, &end); --seen; }
        else if (!strcmp(k, "rmax"))        { p->rmax = strtod(v, &end); --seen; }
        else if (!strcmp(k, "budget_wait")) { p->budget_wait = (int)strtol(v, &end, 10); --seen; }
        else if (!strcmp(k, "tip_order"))   { if (!strcmp(v, "id")) p->tip_order = 0; else if (!strcmp(v, "random")) p->tip_order = 1; else return -1; end = (char*)v + strlen(v); --seen; }
        else return -1;
        if (*end != '\0' || end == v) return -1;
        ++seen;
    }
    return seen == 8 ? 0 : -1;
}

int state_read(World* w, const char* path, double* m_pool_file, char* err, size_t n){
    memset(w, 0, sizeof *w);
    FILE* f = fopen(path, "r");
    if (!f) { snprintf(err, n, "cannot open %s", path); return -1; }
    char line[4096]; int sec = 0, have = 0; long lineno = 0;
    while (fgets(line, sizeof line, f)) {
        ++lineno;
        size_t L = strlen(line);
        if (L == 0 || line[L - 1] != '\n') FAIL("%s:%ld: line too long or not terminated", path, lineno);
        line[--L] = '\0';
        if (line[0] == '#' || L == 0) continue;
        if (line[0] == '[') {
            if      (!strcmp(line, "[header]")) sec = 1;
            else if (!strcmp(line, "[disks]"))  sec = 2;
            else if (!strcmp(line, "[tips]"))   sec = 3;
            else if (!strcmp(line, "[patches]")) sec = 4;
            else FAIL("%s:%ld: unknown section %s", path, lineno, line);
            continue;
        }
        if (sec == 1) {
            char key[32]; int off = 0;
            if (sscanf(line, "%31s %n", key, &off) != 1) FAIL("%s:%ld: bad header line", path, lineno);
            const char* v = line + off; char* end;
            if      (!strcmp(key, "format"))   { if (strtol(v, &end, 10) != PID_FORMAT || *end) FAIL("%s: format is not %d", path, PID_FORMAT); }
            else if (!strcmp(key, "unit"))     { if (strcmp(v, "sigma")) FAIL("%s: unit is not sigma", path); }
            else if (!strcmp(key, "geometry")) { if (strcmp(v, PID_GEOMETRY)) FAIL("%s: unknown geometry convention", path); }
            else if (!strcmp(key, "seed"))     { w->seed = strtoull(v, &end, 10); if (*end) FAIL("%s: bad seed", path); }
            else if (!strcmp(key, "rule"))     { snprintf(w->rule, sizeof w->rule, "%s", v); }
            else if (!strcmp(key, "tick"))     { w->tick = strtol(v, &end, 10); if (*end) FAIL("%s: bad tick", path); }
            else if (!strcmp(key, "box"))      { if (sscanf(v, "%lf %lf", &w->p.W, &w->p.H) != 2) FAIL("%s: bad box", path); }
            else if (!strcmp(key, "budget"))   { w->p.budget = strtod(v, &end); if (*end) FAIL("%s: bad budget", path); }
            else if (!strcmp(key, "m_used"))   { w->m_used = strtod(v, &end); if (*end) FAIL("%s: bad m_used", path); }
            else if (!strcmp(key, "m_pool"))   { *m_pool_file = strtod(v, &end); if (*end) FAIL("%s: bad m_pool", path); }
            else if (!strcmp(key, "params"))   { if (parse_params(&w->p, v)) FAIL("%s: bad params line", path); }
            else if (!strcmp(key, "build"))    { snprintf(w->build, sizeof w->build, "%s", v); }
            else FAIL("%s:%ld: unknown header key %s", path, lineno, key);
            have |= 1;
        } else if (sec == 2) {
            Disk d; char kind[16];
            if (sscanf(line, "%d %lf %lf %lf %d %15s %ld %d", &d.id, &d.x, &d.y, &d.r, &d.parent, kind, &d.born, &d.alive) != 8)
                FAIL("%s:%ld: bad disk line", path, lineno);
            if      (!strcmp(kind, "seed")) d.kind = KIND_SEED;
            else if (!strcmp(kind, "body")) d.kind = KIND_BODY;
            else if (!strcmp(kind, "tip"))  d.kind = KIND_TIP;
            else FAIL("%s:%ld: bad disk kind %s", path, lineno, kind);
            if (d.id != w->nd) FAIL("%s:%ld: disk ids must be 0, 1, 2, ... in order", path, lineno);
            world_add_disk(w, d);
        } else if (sec == 3) {
            Tip t;
            if (sscanf(line, "%d %lf %d", &t.disk, &t.heading, &t.alive) != 3) FAIL("%s:%ld: bad tip line", path, lineno);
            world_add_tip(w, t);
        } else if (sec == 4) {
            int id; Patch q;
            if (w->p.npatch >= PID_MAX_PATCHES || sscanf(line, "%d %lf %lf %lf %lf %d", &id, &q.x, &q.y, &q.r, &q.s, &q.active) != 6 || id != w->p.npatch)
                FAIL("%s:%ld: bad patch line", path, lineno);
            w->p.patch[w->p.npatch++] = q;
        } else FAIL("%s:%ld: record outside a section", path, lineno);
    }
    if (ferror(f)) FAIL("%s: read error", path);
    fclose(f);
    if (!have) { snprintf(err, n, "%s: no header", path); world_free(w); return -1; }
    return 0;
}
