#include "globals.h"
#include "piston.h"

// Move piston based on simulation parameters
void movePiston() {
    xpist1 += vpist1 * DELT;
    xpist2 += vpist2 * DELT;
}