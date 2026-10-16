/* pid_rng.c -- the one generator of Demo 0: xoshiro256** (Blackman and Vigna) seeded through splitmix64.
   Every random draw of the program goes through these functions; --seed is the only seed. */
#include <math.h>
#include "pid.h"

static uint64_t splitmix64(uint64_t* x){
    uint64_t z = (*x += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}
static inline uint64_t rotl(uint64_t x, int k){ return (x << k) | (x >> (64 - k)); }

void rng_seed(Rng* g, uint64_t seed){
    uint64_t x = seed;
    for (int i = 0; i < 4; ++i) g->s[i] = splitmix64(&x);
}
uint64_t rng_next(Rng* g){
    uint64_t* s = g->s;
    const uint64_t result = rotl(s[1] * 5, 7) * 9;
    const uint64_t t = s[1] << 17;
    s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3];
    s[2] ^= t;
    s[3] = rotl(s[3], 45);
    return result;
}
double rng_uniform(Rng* g){ return (double)(rng_next(g) >> 11) * 0x1.0p-53; }
double rng_normal(Rng* g){
    const double u1 = ((double)(rng_next(g) >> 11) + 1.0) * 0x1.0p-53;   /* (0, 1]: log(u1) is finite */
    const double u2 = rng_uniform(g);
    return sqrt(-2.0 * log(u1)) * cos(2.0 * PID_PI * u2);
}
