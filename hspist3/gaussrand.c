// returns a pair of gaussian distributed random numbers
// box-muller algorithm, see http://www.taygeta.com/random/gaussian.html
// one issue, use the exernal definition of _seed in myrand4.h for a continued random stream
// runs at about 5.5 mits on the macbook, returns two values per iteration.
// rob shaw may 2011

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gaussrand.h"


void gaussrand(double* g1, double* g2) {
    static int hasSpare = 0;
    static double spare;

    if (hasSpare) {
        hasSpare = 0;
        *g1 = spare;
        *g2 = (rand() / ((double) RAND_MAX)) * 2.0 - 1.0;  // Generate a new random value for g2
    } else {
        double u, v, s;
        do {
            u = (rand() / ((double) RAND_MAX)) * 2.0 - 1.0;
            v = (rand() / ((double) RAND_MAX)) * 2.0 - 1.0;
            s = u * u + v * v;
        } while (s >= 1.0 || s == 0.0);

        s = sqrt(-2.0 * log(s) / s);
        spare = v * s;
        *g1 = u * s;
        *g2 = v * s;
        hasSpare = 1;
    }
}