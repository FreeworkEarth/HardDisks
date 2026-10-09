/* ##CHRIS 2026-10-09 (261012 sec. 4.7.6; milestone M2): white-box test of the spring divider's contact rule. The schedule
   audit compares the heap with body_rule itself, so it cannot test body_rule's root search; this file includes
   edmd_gen3.c and compares body_rule, for random disks and spring dividers (seeded), with an independent search: a scan
   of the face gap g(t) in steps of T/2048 up to a horizon, the first fall from > 0 to <= 0 refined by bisection. Slow
   approaches whose contact lies beyond the scan's horizon (the O(1) jump of body_rule) are checked by the contact itself:
   g(dt) = 0 to rounding, every local minimum of the scanned periods > 0, and the minimum one period before the contact > 0.
   build (from hspist3/): cc -std=c11 -O2 -ffp-contract=off -Wall -Wextra -o <out> edmd_core/tests/gen3_body_rule_test.c -lm
   usage: <out> [n]   (n random cases, default 20000) */
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
    return bad ? 1 : 0;
}
