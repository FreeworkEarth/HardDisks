#include "globals.h"
#include "physicssimulation.h"
#include <math.h>
#include "constants.h"
#include <math.h>


REAL dt1 = 0.1, dt2 = 0.2; // Initial values; these can be updated later
REAL pdelay = 0.0; // Initialize pdelay here

// Function to update dt1 and dt2 from user input
void setPistonTimes(REAL user_dt1, REAL user_dt2) {
    dt1 = user_dt1;
    dt2 = user_dt2;
    pdelay = dt1 + 0.00; // Recalculate pdelay based on new dt1 value
}

// Function to update particle positions
void update_particles() {
    for (int i = 0; i < nspheres; i++) {
        X[i] += VX[i] * (i % 2 == 0 ? dt1 : dt2);  // Update X position based on velocity and time step
        Y[i] += VY[i] * (i % 2 == 0 ? dt1 : dt2);  // Update Y position based on velocity and time step

        // Apply boundary conditions (reflective walls)
        if (X[i] < 0) {
            X[i] = 0;  // Reflect particle off left boundary
            VX[i] = -VX[i];  // Reverse velocity
        } else if (X[i] > MAX_X) {
            X[i] = MAX_X;  // Reflect particle off right boundary
            VX[i] = -VX[i];  // Reverse velocity
        }

        if (Y[i] < 0) {
            Y[i] = 0;  // Reflect particle off bottom boundary
            VY[i] = -VY[i];  // Reverse velocity
        } else if (Y[i] > MAX_Y) {
            Y[i] = MAX_Y;  // Reflect particle off top boundary
            VY[i] = -VY[i];  // Reverse velocity
        }
    }
}

// Function to detect and resolve collisions
void handle_collisions() {
    for (int i = 0; i < nspheres; i++) {
        for (int j = i + 1; j < nspheres; j++) {
            REAL dx = X[j] - X[i];
            REAL dy = Y[j] - Y[i];
            REAL dist = sqrt(dx * dx + dy * dy);
            if (dist < Radius[i] + Radius[j]) {
                // Simple elastic collision (swap velocities)
                REAL temp = VX[i];
                VX[i] = VX[j];
                VX[j] = temp;

                // You could also handle the Y velocities, implement energy conservation, etc.
                temp = VY[i];
                VY[i] = VY[j];
                VY[j] = temp;
            }
        }
    }
}

void initialize_simulation() {
    xpist1 = XW1;
    xpist2 = rlo;
    initgauss(nrad1, nrad2, RAD1, 0);
    reseteng();
    
    // Initial pdelay calculation
    pdelay = dt1 + 0.00;
    lvel = 1;
    rvel = 1;
    KE0 = kengy();
    KE1 = KE2 = 0;

    // Configure different modes
    switch (modectl) {
    case 1:
        lvel = 0.1;
        rvel = 0.1;
        dt1 = (lhi - llo) / lvel;
        pdelay = dt1 + 0.30;
        break;
    case 2:
        lvel = 2;
        rvel = 0.1;
        dt1 = (lhi - llo) / lvel;
        pdelay = dt1 + 0.30;
        break;
    case 3:
        lvel = 2;
        rvel = 2;
        dt1 = (lhi - llo) / lvel;
        pdelay = dt1 + 1;
        break;
    case 4:
        lvel = 2;
        rvel = 2;
        dt1 = (lhi - llo) / lvel;
        pdelay = dt1 + 0.08;
        break;
    case 5:
        lvel = 2;
        rvel = 2;
        dt1 = (lhi - llo) / lvel;
        pdelay = 0;
        break;
    }

    keyslowctl = 2;
}