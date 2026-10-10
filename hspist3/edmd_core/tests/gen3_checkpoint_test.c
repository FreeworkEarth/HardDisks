/* ##CHRIS 2026-10-09 (stage H; 261012 sec. 4.7.12 stage H, sec. 4.7 item 6): the engine-level test of checkpoint and restart at an
   event boundary (edmd3_checkpoint_write / edmd3_checkpoint_read). The cells and protocols are the M2 harness's own (its source is
   included unchanged; its main is renamed and not called): the cradles (exact, round, round loaded at 8100 units), the free,
   heavy and held dividers at pi/8 and 0.70, the driven divider, the spring divider and the piston push, at the harness's full
   lengths. For each cell:
     U  the uninterrupted run (no checkpoint);
     C  the same run writing a checkpoint after sample k and going on    -> must equal U: writing does not steer;
     R  a new state, created and loaded from the same cell as a restarting driver builds it, the checkpoint read into it, then
        on from sample k                                                  -> must equal U;
     R' the same, but the new state first runs to another sample (k/2) before the read: the read must overwrite everything;
     R2 a chain: R's state writes a second checkpoint at a later sample k2, a new state reads it and goes on -> must equal U;
   at every k of: sample 1, a third of the run, the last sample before the first origin shift and the first after it (where
   the run reaches 2^13 units), two thirds, and the second-last sample. "Equal" is: the event hash, every health counter
   (the whole EDMD3_Health), the momentum and energy ledgers (the whole EDMD3_Ledger), the body states and impulses, and every
   disk's synchronised position and velocity, all compared bit for bit. Two cells are also run with both audits on (contact
   audit at every event, schedule audit every 500 events), where the audit totals must be equal too. Two refusals: a
   checkpoint read into a state of another cell, and a truncated checkpoint.
   build (from hspist3/): cc -std=c11 -O3 -ffp-contract=off -Wall -Wextra -o <out> edmd_core/tests/gen3_checkpoint_test.c
                          edmd_core/edmd_gen3.c edmd_core/edmd.c -lm
   usage: <out>   (prints the table; exit 0 iff every comparison is equal and both refusals refuse) */
#define main gen3_m2_harness_main
#include "gen3_m2_harness.c"
#undef main

typedef struct {
    uint64_t hash; EDMD3_Health H; EDMD3_Ledger L; EDMD3_Audit A; double body[8]; double cmax[4]; long cev;
    EDMD_Particle* P; int N; long events;
} Snap;

static int nsamples(const Cell* c){ int n = 0; while ((n + 1) * DTS <= c->T + 1e-9) ++n; return n; }

/* the harness's protocol from sample k0 (exclusive) to k1 (inclusive): advance on the 0.25-sigma-time grid, the changes after
   the advance to their time (run3's loop without its measurements) */
static void drive(EDMD3* S, const Cell* c, int k0, int k1){
    Change ch[MAXCH]; const int nch = changes(c, ch);
    int ich = 0; while (ich < nch && ch[ich].t < k0 * DTS + 1e-9) ++ich;     /* the changes up to and including sample k0 are done */
    for (int k = k0 + 1; k <= k1; ++k) {
        const double t = k * DTS;
        edmd3_advance_to(S, c->t0 + t);
        while (ich < nch && fabs(ch[ich].t - t) < 1e-9) {
            if (ch[ich].kind == 0) edmd3_set_divider_motion(S, 0, c->dM, 0.0);
            else if (ch[ich].kind == 3) edmd3_set_divider_motion(S, 0, 0.0, ch[ich].val);
            else edmd3_set_piston_motion(S, 1, 0.0, ch[ich].val);
            ++ich;
        }
    }
}
static EDMD3* fresh(const Cell* c, int audits){
    char err[256]; EDMD_Params p = params(c);
    EDMD3* S = edmd3_create(&p, g_cell_px, err, sizeof err);
    if (!S || !edmd3_load(S, c->P, c->t0, err, sizeof err)) { fprintf(stderr, "%s: %s\n", c->name, err); exit(2); }
    if (audits) { edmd3_set_contact_audit(S, 1); edmd3_set_schedule_audit(S, 500); }
    return S;
}
static void snap(EDMD3* S, const Cell* c, Snap* s){
    memset(s, 0, sizeof *s);
    s->hash = edmd3_event_hash(S); s->H = *edmd3_health(S); s->A = *edmd3_schedule_audit_stats(S);
    s->cev = edmd3_contact_audit_stats4(S, s->cmax);
    edmd3_ledger(S, &s->L);
    if (c->ndiv) { edmd3_divider_state(S, 0, &s->body[0], &s->body[1]); s->body[2] = edmd3_divider_impulse(S, 0, 0); s->body[3] = edmd3_divider_impulse(S, 0, 1);
                   s->body[4] = edmd3_work_divider(S, 0); }
    if (c->pistons) { edmd3_piston_state(S, 1, &s->body[5], &s->body[6]); s->body[7] = edmd3_work_piston(S, 1); }
    s->events = s->H.ev_pair + s->H.ev_wall + s->H.ev_cross + s->H.ev_div + s->H.ev_piston + s->H.ev_band;
    s->N = c->N; s->P = (EDMD_Particle*)malloc((size_t)c->N * sizeof(EDMD_Particle));
    const EDMD_Particle* P = edmd3_particles(S);           /* after the counters: this reader runs the full check */
    memcpy(s->P, P, (size_t)c->N * sizeof(EDMD_Particle));
}
static int same(const Snap* a, const Snap* b){
    if (a->hash != b->hash || memcmp(&a->H, &b->H, sizeof a->H) || memcmp(&a->L, &b->L, sizeof a->L) || memcmp(&a->A, &b->A, sizeof a->A)
        || memcmp(a->body, b->body, sizeof a->body) || memcmp(a->cmax, b->cmax, sizeof a->cmax) || a->cev != b->cev) return 0;
    for (int i = 0; i < a->N; ++i)
        if (memcmp(&a->P[i].x, &b->P[i].x, sizeof(double)) || memcmp(&a->P[i].y, &b->P[i].y, sizeof(double)) || memcmp(&a->P[i].vx, &b->P[i].vx, sizeof(double))
            || memcmp(&a->P[i].vy, &b->P[i].vy, sizeof(double)) || a->P[i].coll_count != b->P[i].coll_count) return 0;
    return 1;
}
static FILE* ckpt_at(EDMD3* S, long* bytes){
    FILE* f = tmpfile();
    if (!f || !edmd3_checkpoint_write(S, f) || fflush(f)) { fprintf(stderr, "checkpoint write failed\n"); exit(2); }
    *bytes = ftell(f); rewind(f);
    return f;
}
static void restore(EDMD3* S, FILE* f, const char* who){
    char err[256]; rewind(f);
    if (!edmd3_checkpoint_read(S, f, err, sizeof err)) { fprintf(stderr, "%s: %s\n", who, err); exit(2); }
}

static int g_bad = 0, g_cmp = 0;
static void test_cell(Cell* c, int audits){
    const int K = nsamples(c);
    EDMD3* S = fresh(c, audits); drive(S, c, 0, K); Snap U; snap(S, c, &U); const long shifts = U.H.origin_shifts; edmd3_destroy(S);
    /* the samples: 1, K/3, around the first origin shift, 2K/3, K - 1 */
    int ks[8], nk = 0; ks[nk++] = 1; ks[nk++] = K / 3;
    { const double ts = EDMD3_ORIGIN_SHIFT * floor(c->t0 / EDMD3_ORIGIN_SHIFT) + EDMD3_ORIGIN_SHIFT;   /* the first shift, absolute */
      const int kb = (int)floor((ts - c->t0) / DTS - 1e-9);                                             /* last sample before it */
      if (kb >= 1 && kb + 1 < K) { ks[nk++] = kb; ks[nk++] = kb + 1; } }
    ks[nk++] = 2 * K / 3; ks[nk++] = K - 1;
    for (int q = 0; q < nk; ++q) {
        const int k = ks[q];
        int dup = 0; for (int p = 0; p < q; ++p) if (ks[p] == k) dup = 1;
        if (dup || k < 1 || k >= K) continue;
        long bytes = 0, bytes2 = 0;
        /* C: write at k and go on */
        EDMD3* C = fresh(c, audits); drive(C, c, 0, k);
        const long sh_k = edmd3_health(C)->origin_shifts; const double t_k = edmd3_time(C);
        FILE* f = ckpt_at(C, &bytes); drive(C, c, k, K); Snap sc; snap(C, c, &sc); edmd3_destroy(C);
        /* R: a fresh state reads it */
        EDMD3* R = fresh(c, audits); restore(R, f, "R"); drive(R, c, k, K); Snap sr; snap(R, c, &sr); edmd3_destroy(R);
        /* R': a state that ran elsewhere first reads it */
        EDMD3* Q = fresh(c, audits); drive(Q, c, 0, k / 2 > 0 ? k / 2 : 1); restore(Q, f, "R'"); drive(Q, c, k, K); Snap sq; snap(Q, c, &sq); edmd3_destroy(Q);
        /* R2: the chain, a second checkpoint at k2 from a restored state */
        const int k2 = k + (K - k) / 2;
        int chain = -1; Snap s2; memset(&s2, 0, sizeof s2);
        if (k2 > k && k2 < K) {
            EDMD3* A = fresh(c, audits); restore(A, f, "R2a"); drive(A, c, k, k2); FILE* g = ckpt_at(A, &bytes2); edmd3_destroy(A);
            EDMD3* B = fresh(c, audits); restore(B, g, "R2b"); drive(B, c, k2, K); snap(B, c, &s2); edmd3_destroy(B); fclose(g);
            chain = same(&s2, &U);
        }
        fclose(f);
        const int eC = same(&sc, &U), eR = same(&sr, &U), eQ = same(&sq, &U);
        g_cmp += 3 + (chain >= 0); g_bad += !eC + !eR + !eQ + (chain == 0);
        printf("| %s%s | %ld | %ld | %d (t = %.6g units, %ld shifts) | %s | %s | %s | %s | %ld |\n", c->name, audits ? " (audits on)" : "", U.events, shifts,
               k, t_k, sh_k, eC ? "equal" : "**DIFFERENT**", eR ? "equal" : "**DIFFERENT**", eQ ? "equal" : "**DIFFERENT**",
               chain < 0 ? "-" : (chain ? "equal" : "**DIFFERENT**"), bytes);
        free(sc.P); free(sr.P); free(sq.P); free(s2.P);
    }
    free(U.P);
}

static int test_refusals(Cell* a, Cell* b){
    char err[256]; int ok = 1; long bytes;
    EDMD3* S = fresh(a, 0); drive(S, a, 0, 10); FILE* f = ckpt_at(S, &bytes); edmd3_destroy(S);
    EDMD3* T = fresh(b, 0); rewind(f);
    const int r1 = edmd3_checkpoint_read(T, f, err, sizeof err); edmd3_destroy(T);
    printf("| read into a state of another cell (%s into %s) | %s | %s |\n", a->name, b->name, r1 ? "**ACCEPTED**" : "refused", r1 ? "" : err);
    ok &= !r1;
    FILE* h = tmpfile(); rewind(f);
    char* buf = (char*)malloc((size_t)bytes); const size_t got = fread(buf, 1, (size_t)bytes, f);
    fwrite(buf, 1, got / 2, h); rewind(h); free(buf);
    EDMD3* U = fresh(a, 0);
    const int r2 = edmd3_checkpoint_read(U, h, err, sizeof err); edmd3_destroy(U);
    printf("| a truncated checkpoint (%s, half of %ld bytes) | %s | %s |\n", a->name, bytes, r2 ? "**ACCEPTED**" : "refused", r2 ? "" : err);
    ok &= !r2;
    fclose(f); fclose(h);
    return ok;
}

int main(void){
    printf("# gen3 checkpoint and restart at an event boundary (261012 sec. 4.7.12 stage H), printed by edmd_core/tests/gen3_checkpoint_test.c\n\n");
    printf("build: %s, double %zu bytes; the M2 harness's cells and protocols at its full lengths\n\n", __VERSION__, sizeof(double));
    printf("| cell | events (all kinds) | origin shifts | checkpoint after sample k | C (wrote, went on) = U | R (fresh state read it) = U | "
           "R' (state that ran elsewhere read it) = U | R2 (chain: second checkpoint at k + (K-k)/2) = U | checkpoint bytes |\n"
           "|---|---|---|---|---|---|---|---|---|\n");
    const double T = 400 * SIGT, trel = 40 * SIGT;
    const double L8 = 10.0 * PX, L70 = 5.604167 * PX;
    Cell ce; memset(&ce, 0, sizeof ce); ce.name = "cradle_exact"; make_cradle(&ce, 1); ce.T = 100 * SIGT; test_cell(&ce, 0);
    Cell cr; memset(&cr, 0, sizeof cr); cr.name = "cradle_round"; make_cradle(&cr, 0); cr.T = 100 * SIGT; test_cell(&cr, 0);
    Cell cl; memset(&cl, 0, sizeof cl); cl.name = "cradle_round_late"; make_cradle(&cl, 0); cl.T = 100 * SIGT; cl.t0 = 8100.0; test_cell(&cl, 0);
    struct { const char* n; double L; int dense; uint64_t seed; double M, trel; } dv[7] = {
        {"free_pi8_M50", L8, 0, 0xA11CE1ULL, 50.0, trel}, {"free_pi8_M500", L8, 0, 0xA11CE2ULL, 500.0, trel},
        {"free_070_M50", L70, 1, 0xA11CE3ULL, 50.0, trel}, {"free_070_M500", L70, 1, 0xA11CE4ULL, 500.0, trel},
        {"heavy_070_M4e7", L70, 1, 0xA11CE5ULL, 4e7, trel}, {"held_pi8", L8, 0, 0xA11CE6ULL, 0.0, -1.0}, {"held_070", L70, 1, 0xA11CE7ULL, 0.0, -1.0}};
    Cell keep_a, keep_b; memset(&keep_a, 0, sizeof keep_a); memset(&keep_b, 0, sizeof keep_b);
    for (int k = 0; k < 7; ++k) {
        Cell c = divided(dv[k].n, "", dv[k].L, dv[k].dense, dv[k].seed, dv[k].M, dv[k].trel, T);
        test_cell(&c, 0);
        if (k == 3) test_cell(&c, 1);
        if (k == 0) keep_a = c; else if (k == 2) keep_b = c; else free(c.P);
    }
    { Cell c = divided("driven_pi8", "", L8, 0, 0xA11CEAULL, 0.0, trel, T); c.v_drive = 0.2; c.t_switch = 10.0 * SIGT; test_cell(&c, 0); free(c.P); }
    {   Cell c; memset(&c, 0, sizeof c); c.name = "spring_pi8"; rng_state = 0xA11CE8ULL;
        const double Lg = 2.0 * L8; c.N = 400; c.NL = 400; c.H = 40.0 * PX; c.ndiv = 1; c.dx = Lg + 0.5 * TH; c.W = Lg + TH + L8;
        c.P = (EDMD_Particle*)calloc(400, sizeof(EDMD_Particle));
        if (!rsa(c.P, 0, 400, 0.0, Lg, c.H)) { fprintf(stderr, "spring: seeding failed\n"); return 2; }
        c.dk = 0.02; c.dxeq = c.dx - 400.0 * 2.76 / (c.dk * Lg); c.dM = 50.0; c.t_rel = trel; c.T = T;
        test_cell(&c, 0); free(c.P); }
    {   Cell c; memset(&c, 0, sizeof c); c.name = "piston_push"; const double Lp = 77.5 * PX; rng_state = 0xA11CE9ULL;
        c.N = 400; c.NL = 200; c.H = 20.0 * PX; c.W = 2.0 * Lp + TH; c.ndiv = 1; c.dx = Lp + 0.5 * TH;
        c.P = (EDMD_Particle*)calloc(400, sizeof(EDMD_Particle));
        if (!rsa(c.P, 0, 200, 0.0, Lp, c.H) || !rsa(c.P, 200, 400, Lp + TH, c.W, c.H)) { fprintf(stderr, "piston: seeding failed\n"); return 2; }
        c.dM = 100.0; c.t_rel = trel; c.T = T; c.pistons = 1; c.u_p = 0.05; c.t_p0 = trel; c.t_p1 = trel + 155.0 * SIGT;
        test_cell(&c, 0); test_cell(&c, 1); free(c.P); }
    printf("\n## Refusals\n\n| case | result | message |\n|---|---|---|\n");
    const int refused = test_refusals(&keep_a, &keep_b);
    free(keep_a.P); free(keep_b.P);
    printf("\nVERDICT: %s -- %d comparisons, %d different; refusals %s\n", (g_bad == 0 && refused) ? "PASS" : "**FAIL**", g_cmp, g_bad,
           refused ? "both refused" : "**NOT REFUSED**");
    return (g_bad == 0 && refused) ? 0 : 1;
}
