// fades.c - hopefully screen fade control...
// cc -c fades.c -framework ApplicationServices
#include <stdio.h>
#include <unistd.h>
#include <ApplicationServices/ApplicationServices.h>

#define kMyFadeTime    0.3
#define kMyFadeSteps    40
#define kMyFadeInterval    ( kMyFadeTime/(double) kMyFadeSteps)

CGGammaValue    redMin, redMax, redGamma,
                greenMin, greenMax, greenGamma,
                blueMin, blueMax, blueGamma;
int step;
double fadeValue;

void fadeout()
{
CGGetDisplayTransferByFormula (kCGDirectMainDisplay,
            &redMin, &redMax, &redGamma,
            &greenMin, &greenMax, &greenGamma,
            &blueMin, &blueMax, &blueGamma );

for ( step = 0; step <= kMyFadeSteps; ++step )
	{
	fadeValue = 1.0 - step/(double)kMyFadeSteps;
//	printf("%f fadevalue\n",fadeValue);
    CGSetDisplayTransferByFormula (kCGDirectMainDisplay,
            redMin, fadeValue*redMax, redGamma,
            greenMin, fadeValue*greenMax, greenGamma,
            blueMin, fadeValue*blueMax, blueGamma );

    usleep( (useconds_t)(1000000.0 * kMyFadeInterval) );
	}
}

void fadein()
{
for ( step = 0; step <= kMyFadeSteps; ++step )
	{
	fadeValue = step / (double)kMyFadeSteps;
    CGSetDisplayTransferByFormula( kCGDirectMainDisplay,
            redMin, fadeValue*redMax, redGamma,
            greenMin, fadeValue*greenMax, greenGamma,
            blueMin, fadeValue*blueMax, blueGamma );
    usleep( (useconds_t)(1000000.0 * kMyFadeInterval) );
	}
CGDisplayRestoreColorSyncSettings();
}

void restorecolor()
{
CGDisplayRestoreColorSyncSettings();
}
