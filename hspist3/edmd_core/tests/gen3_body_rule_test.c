/* ##CHRIS 2026-10-09 (261012 sec. 4.7.6; milestone M2): white-box test of the spring divider's contact rule. The schedule
   audit compares the heap with body_rule itself, so it cannot test body_rule's root search; this file includes
   edmd_gen3.c and compares body_rule, for random disks and spring dividers (seeded), with an independent search: a scan
   of the face gap g(t) in steps of T/2048 up to a horizon, the first fall from > 0 to <= 0 refined by bisection. Slow
   approaches whose contact lies beyond the scan's horizon (the O(1) jump of body_rule) are checked by the contact itself:
   g(dt) = 0 to rounding, every local minimum of the scanned periods > 0, and the minimum one period before the contact > 0.
   build (from hspist3/): cc -std=c11 -O2 -ffp-contract=off -Wall -Wextra -o <out> edmd_core/tests/gen3_body_rule_test.c -lm
   usage: <out> [n]   (n random cases, default 20000) */
/* ##CHRIS 2026-10-09 (M3, 261012 sec. 4.7.14, amendment b): CONSTRUCTED cases after the random ones (whose output line is
   unchanged), counted per category, each with an answer known by construction:
   near-tangent    the gap g has a double root: its first local minimum is set to -delta (touching: a contact just before
                   the minimum) or +delta (missing: none there), delta = 1e-12 x the position scale; with drift 0 (the first
                   two periods) and with the disk approaching (drift > 0: the k-th minimum is the tangent one, the jump path)
   turning point   a transversal contact exactly when the divider turns (velocity 0) at its extreme towards the disk, in
                   the first periods or the k-th, the disk slower than the face (oscillating g) or faster (monotone g)
   slow approach   drift u with u T << the gap: the contact in the falling stretch before the k-th minimum, k up to 1e9
   Shifting the disk's start position moves g by a constant and leaves its shape (minima times, drift) unchanged, so a
   minimum can be put at any value. */
#include "../edmd_gen3.c"

static uint64_t rs = 0x5EEDB0D1ULL;
static double worst_ulps = 0.0;
static double math_ulp(double x){ return nextafter(fabs(x), INFINITY) - fabs(x); }
static double ur(void){ rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return ((double)(rs >> 11) + 0.5) / 9007199254740992.0; }

/* independent: first t in (0, tlim] with g falling from > 0 to <= 0, scanning in steps h; 0 if none; at once rule as body_rule */
static int scan(const Harm* H, double g0, double gp0, int mlast, double h, double tlim, double* tc){
    if (g0 <= 0.0 && gp0 < 0.0 && !mlast) { *tc = 0.0; return 2; }
    double ta = 0.0, ga = harm_g(H, 0.0);
    for (long k = 1; ; ++k) {
        double tb = k * h; if (tb > tlim) tb = tlim;
        const double gb = harm_g(H, tb);
        if (ga > 0.0 && gb <= 0.0) {
            double lo = ta, hi = tb;
            for (int it = 0; it < 200; ++it) { const double m = lo + 0.5 * (hi - lo); if (m <= lo || m >= hi) break; if (harm_g(H, m) <= 0.0) hi = m; else lo = m; }
            *tc = hi; return 1;
        }
        if (tb >= tlim) return 0;
        ta = tb; ga = gb;
    }
}

/* ---------------------------------------------------------------- M3, amendment b: constructed cases */
typedef struct { const char* name; long n, bad; double worst; } Cat;
static Cat cats[] = {
    {"near-tangent touching, drift 0 (first two periods)", 0, 0, 0},
    {"near-tangent missing, drift 0", 0, 0, 0},
    {"near-tangent touching at the k-th minimum (jump path)", 0, 0, 0},
    {"near-tangent missing at the k-th minimum (contact one period later)", 0, 0, 0},
    {"turning point, oscillating g, first period", 0, 0, 0},
    {"turning point, oscillating g, k-th period (jump path)", 0, 0, 0},
    {"turning point, monotone g (disk faster than the face)", 0, 0, 0},
    {"slow approach, k = 10 .. 1e9 periods (jump path)", 0, 0, 0}};
enum { C_TT, C_TM, C_JT, C_JM, C_TP1, C_TPK, C_TPM, C_SLOW, NCAT };
static void cat_fail(int c, long q, const char* what, double a, double b){
    cats[c].bad++;
    if (cats[c].bad <= 5) printf("constructed %s (case %ld): %s %.17g %.17g\n", cats[c].name, q, what, a, b);
}
/* a random spring divider at t = 0 (as the random cases) */
static void rand_spring(Obj3* O, double* A){
    memset(O, 0, sizeof *O);
    O->active = 1; O->harmonic = 1; O->th = 1.2; O->h = 0.5 * O->th + 12.0;
    O->M = 50.0 * pow(10.0, 3.0 * ur()); O->k = 0.002 * pow(10.0, 2.0 * ur()); O->omega = sqrt(O->k / O->M);
    O->xeq = 300.0 + 100.0 * (ur() - 0.5); O->tau = 0.0;
    *A = 5.0 + 55.0 * ur(); const double ph = 2.0 * G3_PI * ur();
    O->x = O->xeq + *A * cos(ph); O->v = -*A * O->omega * sin(ph);
}
static Harm harm_of(const Obj3* O, double xp, double vp, int side){
    Harm H; H.sg = side; H.D = O->x - O->xeq; H.V = O->v; H.w = O->omega; H.vp = vp; H.G = side * (O->xeq - xp) - O->h; return H;
}
/* the first local max and min of g in [0, T), as body_rule finds them (oscillating g: Q > |vp|) */
static void extrema(const Harm* H, double* tmx, double* tmn){
    const double Q = hypot(H->D * H->w, H->V), phi = atan2(H->D * H->w, H->V), al = acos(H->vp / Q), tp = 2.0 * G3_PI;
    double a = fmod(-phi + H->sg * al, tp), b = fmod(-phi - H->sg * al, tp);
    if (a < 0.0) a += tp;
    if (b < 0.0) b += tp;
    *tmx = a / H->w; *tmn = b / H->w;
}
/* the root of g in [lo, hi] (g(lo) > 0 >= g(hi)) by plain bisection to adjacent doubles: an independent refinement */
static double bisect(const Harm* H, double lo, double hi){
    for (int it = 0; it < 2000; ++it) { const double m = lo + 0.5 * (hi - lo); if (m <= lo || m >= hi) break; if (harm_g(H, m) <= 0.0) hi = m; else lo = m; }
    return hi;
}
/* the rounding scale of g at t: the size of its terms times a few eps */
static double gscale(const Harm* H, double t){
    const double A = hypot(H->D, H->V / H->w);
    return 64.0 * 2.220446049250313e-16 * (fabs(H->G) + A * (1.0 + H->w * fabs(t)) + fabs(H->vp * t) + 1.0);
}
static long constructed(void){
    rs = 0x5EEDC0DEULL;                                        /* its own stream */
    const long per = 2000;
    for (int c = 0; c < NCAT; ++c) {
        for (long q = 0; q < per; ++q) {
            Obj3 O; double A; rand_spring(&O, &A);
            const int side = ur() < 0.5 ? 1 : -1;
            const double T = 2.0 * G3_PI / O.omega;
            double dt = 0.0, g0 = 0.0;
            if (c <= C_JM) {                                  /* near-tangent */
                const int jump = c >= C_JT;
                const double u = jump ? 1e-4 * A * O.omega * (0.1 + ur()) : 0.0;      /* drift, below Q */
                const long kk = jump ? (long)(2 + (long)(99.0 * ur())) : 0;           /* the tangent minimum: k periods after the first */
                double xp = side > 0 ? O.xeq - A - O.h - 50.0 : O.xeq + A + O.h + 50.0;
                const double vp = side * u;
                Harm H = harm_of(&O, xp, vp, side);
                double tmx, tmn; extrema(&H, &tmx, &tmn);
                if (tmn < 0.1 * T || tmn > 0.9 * T) { --q; continue; }                 /* the minimum away from t = 0 */
                const double tm = tmn + (double)kk * T;                                /* the tangent minimum */
                const double delta = 1e-12 * fabs(O.xeq);
                const double target = (c == C_TT || c == C_JT) ? -delta : delta;
                const double m = harm_g(&H, tm);
                xp += side * (m - target); H = harm_of(&O, xp, vp, side);              /* g shifted by target - m */
                const double mt = harm_g(&H, tm);
                if (fabs(mt - target) > 0.25 * delta) { --q; continue; }                /* rounding of the shift too large: skip */
                if (!(harm_g(&H, 0.0) > 100.0 * delta)) { --q; continue; }
                cats[c].n++;
                const int rc = body_rule(&O, 0, O.x, O.v, xp, vp, 0, &dt, &g0);
                if (c == C_TM && rc != 0) { cat_fail(c, q, "a contact where none is (dt, min)", dt, mt); continue; }
                if (c == C_TM) continue;
                if (rc != 1) { cat_fail(c, q, "no contact (rc, min)", rc, mt); continue; }
                /* the expected falling stretch: the one ending at the tangent minimum (touching) or one period later (missing) */
                const double tend = c == C_JM ? tm + T : tm;
                const double tbeg = tend - (tmn > tmx ? tmn - tmx : tmn - tmx + T);
                if (!(dt > tbeg - 1e-9 * T && dt <= tend + 1e-9 * T)) { cat_fail(c, q, "contact outside the expected falling stretch (dt, tend)", dt, tend); continue; }
                const double tb = bisect(&H, tbeg, tend);
                const double err = fabs(dt - tb) / fmax(1.0, tb);
                if (err > cats[c].worst) cats[c].worst = err;
                if (err > 1e-6 && fabs(harm_g(&H, dt)) > gscale(&H, dt)) cat_fail(c, q, "contact time (engine, bisection)", dt, tb);
            } else if (c <= C_TPM) {                          /* turning point */
                const int mono = c == C_TPM;
                const double Q = A * O.omega;
                const double u = mono ? Q * (2.0 + 8.0 * ur()) : Q * (0.05 + 0.85 * ur());
                const double vp = side * u;
                /* the turning points: O.x(t) = xeq + A cos(w t + psi); towards the disk (side +1: the disk is left): x = xeq - A */
                const double psi = atan2(-O.v / O.omega, O.x - O.xeq);                 /* x - xeq = A cos psi, v = -A w sin psi */
                const double want = side > 0 ? G3_PI : 0.0;                            /* w t + psi = want (mod 2 pi) */
                double ts = fmod(want - psi, 2.0 * G3_PI); if (ts < 0.0) ts += 2.0 * G3_PI;
                ts /= O.omega;
                const long kk = c == C_TPK ? (long)(2 + (long)(999.0 * ur())) : 0;
                if (ts < 0.05 * T) ts += T;
                ts += (double)kk * T;
                /* the disk at contact at ts: centre = turning position - side h; at t = 0 it was vp ts earlier */
                const double xturn = O.xeq + A * cos(O.omega * ts + psi);
                const double xp = xturn - side * O.h - vp * ts;
                Harm H = harm_of(&O, xp, vp, side);
                if (!(harm_g(&H, 0.0) > 0.0)) { --q; continue; }
                cats[c].n++;
                const int rc = body_rule(&O, 0, O.x, O.v, xp, vp, 0, &dt, &g0);
                if (rc != 1) { cat_fail(c, q, "no contact (rc, ts)", rc, ts); continue; }
                const double err = fabs(dt - ts) / fmax(1.0, ts);
                if (err > cats[c].worst) cats[c].worst = err;
                /* g(ts) = 0 up to the rounding of the construction; the engine's contact is the root of the computed g */
                if (fabs(harm_g(&H, dt)) > gscale(&H, dt) || fabs(dt - ts) * u > 4.0 * gscale(&H, ts))
                    cat_fail(c, q, "contact time (engine, constructed)", dt, ts);
            } else {                                          /* slow approach, k periods */
                static const double ks[4] = {10.0, 1e3, 1e6, 1e9};
                const double k = ks[q % 4];
                double xp = side > 0 ? O.xeq - A - O.h - 50.0 : O.xeq + A + O.h + 50.0;
                Harm H0 = harm_of(&O, xp, 0.0, side);
                double tmx, tmn; extrema(&H0, &tmx, &tmn);
                const double m1 = harm_g(&H0, tmn + T);                               /* with vp = 0: every minimum = m1 */
                const double u = m1 / ((k - 0.5) * T);                                 /* the minimum in period k is -u T / 2 */
                const double vp = side * u;
                Harm H = harm_of(&O, xp, vp, side);
                double tmx2, tmn2; extrema(&H, &tmx2, &tmn2);
                cats[c].n++;
                const int rc = body_rule(&O, 0, O.x, O.v, xp, vp, 0, &dt, &g0);
                if (rc != 1) { cat_fail(c, q, "no contact (rc, k)", rc, k); continue; }
                /* the expected stretch, from the drifting g itself: its minimum in [T, 2T) is m2 = g(t1), and every minimum is
                   u T lower than the one before (g(t + T) = g(t) - u T exactly), so the first one at or below zero is the j-th
                   after t1 with j = ceil(m2 / (u T)); the contact is in the falling stretch that ends there */
                const double t1 = tmn2 + T, m2 = harm_g(&H, t1);
                double j = ceil(m2 / (u * T)); if (j < 0.0) j = 0.0;
                while (j > 0.0 && harm_g(&H, t1 + (j - 1.0) * T) <= 0.0) j -= 1.0;   /* rounding guards, both ways */
                while (harm_g(&H, t1 + j * T) > 0.0) j += 1.0;
                const double tend = t1 + j * T;
                const double rel = fabs(dt - tend) / T;
                if (!(rel < 1.0)) { cat_fail(c, q, "contact not in the k-th period (dt, expected minimum)", dt, tend); continue; }
                if (rel > cats[c].worst) cats[c].worst = rel;
                if (fabs(harm_g(&H, dt)) > gscale(&H, dt)) cat_fail(c, q, "g(dt) above its rounding scale (g, scale)", harm_g(&H, dt), gscale(&H, dt));
                if (dt > T && harm_g(&H, dt - T) <= 0.0) cat_fail(c, q, "g(dt - T) <= 0: an earlier contact (dt, k)", dt, k);
            }
        }
    }
    long bad = 0;
    printf("\nconstructed cases (amendment b), %ld per category; worst = max relative contact-time difference to the independent answer "
           "(for the slow approach: max |dt - expected minimum| / T, which is below 1 inside the expected stretch):\n\n"
           "| category | cases | failures | worst |\n|---|---|---|---|\n", per);
    for (int c = 0; c < NCAT; ++c) { printf("| %s | %ld | %ld | %.3g |\n", cats[c].name, cats[c].n, cats[c].bad, cats[c].worst); bad += cats[c].bad; }
    printf("\nconstructed failures: %ld\n", bad);
    return bad;
}

int main(int argc, char** argv){
    const long n = argc > 1 ? atol(argv[1]) : 20000;
    long n1 = 0, n0 = 0, n2 = 0, nslow = 0, bad = 0, nscan = 0; double worst = 0.0, worst_g = 0.0;
    for (long q = 0; q < n; ++q) {
        Obj3 O; memset(&O, 0, sizeof O);
        O.active = 1; O.harmonic = 1; O.th = 1.2; O.h = 0.5 * O.th + 12.0;
        O.M = 50.0 * pow(10.0, 3.0 * ur()); O.k = 0.002 * pow(10.0, 2.0 * ur()); O.omega = sqrt(O.k / O.M);
        O.xeq = 300.0 + 100.0 * (ur() - 0.5); O.tau = 0.0;
        const double A = 60.0 * ur(), ph = 2.0 * G3_PI * ur();
        O.x = O.xeq + A * cos(ph); O.v = -A * O.omega * sin(ph);
        /* a disk on a random side, at a random distance beyond the face's reach or inside it */
        const int side = ur() < 0.5 ? 1 : -1;     /* +1: disk left of the divider */
        const double cmax = O.xeq + A, cmin = O.xeq - A;
        const double xp = side > 0 ? cmin - O.h - 200.0 * ur() + 10.0 * ur() : cmax + O.h + 200.0 * ur() - 10.0 * ur();
        if ((side > 0 && xp >= O.x - O.h) || (side < 0 && xp <= O.x + O.h)) { --q; continue; }   /* must start separated */
        /* velocity: thermal, or slow, or very slow approach */
        const double r = ur(); double vp = 2.0 * (ur() - 0.5);
        if (r < 0.25) vp = side * 1e-3 * ur(); else if (r < 0.35) vp = side * 1e-7 * (0.5 + ur());
        double dt = 0.0, g0 = 0.0;
        const int rc = body_rule(&O, 0, O.x, O.v, xp, vp, 0, &dt, &g0);
        Harm H; H.sg = side; H.D = O.x - O.xeq; H.V = O.v; H.w = O.omega; H.vp = vp; H.G = side * (O.xeq - xp) - O.h;
        const double T = 2.0 * G3_PI / O.omega, h = T / 2048.0;
        const double gp0 = side * (O.v - vp);
        if (rc == 2) n2++; else if (rc == 1) n1++; else n0++;
        double tlim = 64.0 * T;
        if (rc == 1 && dt + T > tlim) {           /* beyond the scan: the contact itself and the minima before it */
            nslow++;
            const double gc = harm_g(&H, dt);
            if (fabs(gc) > 1e-9 * (1.0 + A)) { bad++; printf("slow: g(dt) = %.3g at dt = %.17g (q = %ld)\n", gc, dt, q); }
            if (harm_g(&H, dt - T) <= 0.0 && dt > T) { bad++; printf("slow: g(dt - T) <= 0 (q = %ld)\n", q); }
            if (fabs(gc) > worst_g) worst_g = fabs(gc);
            double tc; const int rs2 = scan(&H, g0, gp0, 0, h, tlim, &tc);   /* within 64 periods it must find the same contact */
            if (rs2 && fabs(tc - dt) > 1e-9 * fmax(1.0, dt)) { bad++; printf("slow: the scan finds a contact at %.17g, the engine at %.17g (q = %ld)\n", tc, dt, q); }
            {   /* the contact error in time: |g(dt)| over |g'(dt)| against the spacing of doubles at dt */
                const double gpc = harm_gp(&H, dt), et = fabs(gc) / fmax(fabs(gpc), 1e-300);
                if (et / math_ulp(dt) > worst_ulps) worst_ulps = et / math_ulp(dt);
            }
            continue;
        }
        nscan++;
        double tc = 0.0; const int rs2 = scan(&H, g0, gp0, 0, h, tlim, &tc);
        if (rs2 != rc) { bad++; printf("mismatch rc engine %d scan %d: dt %.17g tc %.17g (q = %ld, vp %.3g, A %.3g)\n", rc, rs2, dt, tc, q, vp, A); continue; }
        if (rc) { const double d = fabs(dt - tc) / fmax(1.0, tc); if (d > worst) worst = d; if (d > 1e-9) { bad++; printf("dt engine %.17g scan %.17g (q = %ld)\n", dt, tc, q); } }
    }
    printf("gen3 body_rule (spring divider) against an independent scan: %ld cases: contact %ld (of them beyond 64 periods: %ld), "
           "none %ld, at once %ld; compared by the scan %ld; max |dt_engine - dt_scan| / max(1, dt) = %.3g; slow cases max |g(dt)| = %.3g px, max |g(dt)| / |g'(dt)| = %.3g ulp(dt); "
           "failures %ld\n", n, n1, nslow, n0, n2, nscan, worst, worst_g, worst_ulps, bad);
    const long cbad = constructed();   /* ##CHRIS 2026-10-09 (M3, amendment b) */
    return (bad || cbad) ? 1 : 0;
}
