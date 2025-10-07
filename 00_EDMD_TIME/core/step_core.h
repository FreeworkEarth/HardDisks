#ifndef CORE_STEP_CORE_H
#define CORE_STEP_CORE_H

#include "state.h"

#ifdef __cplusplus
extern "C" {
#endif

void remove_drift(void);
double kinetic_energy(void);
float compute_measured_temperature_from_ke(void);
void maxwell_boltzmann_2D(float temperature, float *vx, float *vy);
float maxwell_boltzmann_velocity_gaussians(float temperature);
float sample_gaussian_velocity(float temperature);
void maxwell_boltzmann_velocity_ARRAY_ALL(float temperature);
float maxwell_boltzmann_speed3D(float temperature);

void initialize_simulation_random(void);
void initialize_simulation_honeycomb(void);
void initialize_simulation(void);

void update_particles_with_substepping(float dt, int *out_left, int *out_right);

double kinetic_energy_expected(void);
double average_kinetic_energy_per_particle(void);
double kinetic_energy_total_system(void);
double temperature(void);

void compute_velocity_histogram(void);
void compute_position_histogram(void);

double entropy_kinetic(void);
double entropy_position(void);
float entropy(void);
float position_entropy(void);

double total_free_energy(double T, double S);

#ifdef __cplusplus
}
#endif

#endif /* CORE_STEP_CORE_H */
