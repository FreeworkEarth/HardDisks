#include "io_log.h"

#include <math.h>

void log_energy(FILE *fptr) {
    double ke_particles = kinetic_energy();
    double ke_wall = 0.5 * WALL_MASS * vx_wall * vx_wall;
    double ke_total = ke_particles + ke_wall;
    double T = temperature();

    double S_kinetic = entropy_kinetic();
    double S_position = entropy_position();
    double S_total = S_kinetic + S_position;

    compute_velocity_histogram();
    float S_v = entropy();

    compute_position_histogram();
    float S_p = position_entropy();
    (void)S_p;

    float F = (float)total_free_energy(T, S_total);

    fprintf(fptr,
            "KE_particles: %.8e, KE_wall: %.8e, KE_total: %.8e, T: %.8e, S_v: %.8e, F: %.8e\n",
            ke_particles, ke_wall, ke_total, T, S_v, F);
}

void log_packing_fractions(FILE *log, const char *label) {
    float particle_area_unitless = (float)(M_PI * PARTICLE_RADIUS_UNIT * PARTICLE_RADIUS_UNIT);
    float particle_area_actual = (float)(M_PI * PARTICLE_RADIUS * PARTICLE_RADIUS);

    float box_area_unitless = 2.0f * L0_UNITS * HEIGHT_UNITS;

    float L0_scaled = L0_UNITS * PIXELS_PER_SIGMA;
    float height_scaled = HEIGHT_UNITS * PIXELS_PER_SIGMA;
    float box_area_scaled = 2.0f * L0_scaled * height_scaled;
    float wall_buffer_area = 2.0f * WALL_THICKNESS * height_scaled;

    float packing_unitless = NUM_PARTICLES * particle_area_unitless / box_area_unitless;
    float packing_actual_bound = NUM_PARTICLES * particle_area_actual / box_area_scaled;
    float packing_actual_excl_wall = NUM_PARTICLES * particle_area_actual / (box_area_scaled - wall_buffer_area);

    fprintf(log,
            "Time, Wall_X, Displacement(σ), Left_Count, Right_Count, Packing_UNITLESS, Packing_ACTUAL_BOUND, Packing_ACTUAL_WITH_WALLBUFFER\n");
    (void)label;
    (void)packing_unitless;
    (void)packing_actual_bound;
    (void)packing_actual_excl_wall;
}
