// returns a pair of gaussian distributed random numbers
// box-muller algorithm, see http://www.taygeta.com/random/gaussian.html
// one issue, use the exernal definition of _seed in myrand4.h for a continued random stream
// runs at about 5.5 mits on the macbook, returns two values per iteration.
// rob shaw may 2011

#include <stdio.h>
#include <math.h>
#include "myrand4.h"

void gaussrand(double *g1, double *g2)
{
register double x1,x2,w;

do
	{	
	x1 = f1_rand();
	x2 = f1_rand();
	w = x1*x1 + x2*x2;
	}
while(w > 1.0);

w = sqrt( (-2.0 * log( w ) ) / w );
*g1 = w * x1;
*g2 = w * x2;
}
