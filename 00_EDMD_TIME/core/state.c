#include "state.h"

#include <math.h>
#include <stdio.h>

bool log_packing_fraction = true;

#if EXPERIMENT_MODE
bool enable_speed_of_sound_experiments = true;
#else
bool enable_speed_of_sound_experiments = false;
#endif

int num_steps = (int)(TIME_UNITS_SIMULATION / FIXED_DT);
bool wall_position_is_managed_externally = false;

float L0_UNITS = 20.0f;
float HEIGHT_UNITS = 10.0f;

double simulation_time = 0.0;
double wall_release_time = -1.0;

#if EXPERIMENT_MODE
float wall_thickness_runtime = 1.0f;
#else
float wall_thickness_runtime = 2.0f * PARTICLE_RADIUS;
#endif

void set_wall_thickness_sigma(float thickness_in_sigma) {
    wall_thickness_runtime = thickness_in_sigma * PIXELS_PER_SIGMA;
}

float wall_mass_runtime = PARTICLE_MASS * WALL_MASS_FACTOR;

bool wall_hold_enabled = true;
bool wall_is_released = false;
int wall_hold_steps = WALL_HOLD_STEPS;
int steps_elapsed = 0;
float wall_impulse_x_accum = 0.f;

double TAU_time_scale_factor_to_molecular = 0.0;

void initialize_time_scale(void) {
#if MOLECULAR_MODE
    TAU_time_scale_factor_to_molecular = 1.0;
#else
    TAU_time_scale_factor_to_molecular = sqrtf((PARTICLE_MASS * powf(PARTICLE_RADIUS * METERS_PER_PIXEL, 2)) /
                                              (K_B * TEMPERATURE));
    printf("✅ [Time Scale] Calculated TAU = %.15e seconds\n", TAU_time_scale_factor_to_molecular);
#endif
}

int SIM_WIDTH = 0, SIM_HEIGHT = 0;
int XW1 = 50;
int XW2 = 0, YW1 = 0, YW2 = 0;
int MAX_X = 0, MAX_Y = 0;
int HIST_WIDTH = 0, HIST_HEIGHT = 0;

void initialize_simulation_dimensions(void) {
    SIM_WIDTH = (int)(2 * L0_UNITS * PIXELS_PER_SIGMA);
    SIM_HEIGHT = (int)(HEIGHT_UNITS * PIXELS_PER_SIGMA);

    XW2 = XW1 + SIM_WIDTH;
    YW2 = YW1 + SIM_HEIGHT;

    MAX_X = XW2 + 100;
    MAX_Y = YW2 + 100;

    HIST_WIDTH = MAX_X - SIM_WIDTH;
    HIST_HEIGHT = MAX_Y - SIM_HEIGHT;
}

float piston_left_x = 0.0f;
float piston_right_x = 0.0f;
float vx_piston_left = 0.0f;
float vx_piston_right = 0.0f;

double piston_work_left = 0.0;
double piston_work_right = 0.0;
long piston_hits_left = 0;
long piston_hits_right = 0;
long missed_collision_events = 0;
double worst_penetration_observed = 0.0;

float wall_x_old = 0.0f;
float wall_x = 0.0f;
float vx_wall = 0.0f;
int wall_enabled = 0;

void initialize_simulation_parameters(void) {
    piston_right_x = (float)MAX_X;
    wall_x = XW1 + (L0_UNITS * PIXELS_PER_SIGMA);
    wall_x_old = wall_x;
    vx_wall = 0.0f;

    piston_left_x = 0.0f;
    vx_piston_left = 0.0f;
    vx_piston_right = 0.0f;

    piston_work_left = 0.0;
    piston_work_right = 0.0;
    piston_hits_left = 0;
    piston_hits_right = 0;
    missed_collision_events = 0;
    worst_penetration_observed = 0.0;
}

int left_count = 0;
int right_count = 0;

double X[NUM_PARTICLES] = {0};
double Y[NUM_PARTICLES] = {0};
double X_old[NUM_PARTICLES] = {0};
double V_init[NUM_PARTICLES] = {0};
double Vx[NUM_PARTICLES] = {0};
double Vy[NUM_PARTICLES] = {0};
double Radius[NUM_PARTICLES] = {0};

int velocity_histogram[NUM_BINS] = {0};
int position_histogram[NUM_BINS * NUM_BINS] = {0};
int energy_histogram[NUM_BINS] = {0};

SDL_Window* _window = NULL;
SDL_Renderer* _renderer = NULL;
TTF_Font* font = NULL;

int simulation_started = 0;
bool paused = false;

bool live_fft = false;
float wall_x_samples[FFT_SIZE] = {0};
float energy_samples[FFT_SIZE] = {0};
kiss_fft_scalar in_wall[FFT_SIZE] = {0};
kiss_fft_scalar in_energy[FFT_SIZE] = {0};
kiss_fft_cpx out_wall[FFT_SIZE/2 + 1] = {0};
kiss_fft_cpx out_energy[FFT_SIZE/2 + 1] = {0};
float wall_fft_magnitude[FFT_SIZE/2 + 1] = {0};
float energy_fft_magnitude[FFT_SIZE/2 + 1] = {0};
int sample_index = 0;
kiss_fftr_cfg fft_cfg = NULL;

void print_simulation_info(void) {
    printf("\n========== Simulation Info ==========\n");
    printf("Mode: %s\n", MOLECULAR_MODE ? "MOLECULAR (Normalized units)" : "MACROSCOPIC (Real units)");
    printf("Experiment Mode: %s\n", EXPERIMENT_MODE ? "ON (fine dt/substeps)" : "OFF (default dt)");
    printf("Fixed Time Step (dt): %.2e\n", (double)FIXED_DT);
    printf("TAU scaling factor: %.4f\n", TAU_time_scale_factor_to_molecular);
    printf("Wall hold mode: %s\n", wall_hold_enabled ? "Enabled" : "Disabled");
    printf("Wall hold steps: %d\n", wall_hold_steps);
    printf("Target total simulation time: %.1f\n", 3000.0f);
    printf("Calculated steps to reach target: %d\n",
           (int)(3000.0f / (FIXED_DT / (MOLECULAR_MODE ? 1.0f : TAU_time_scale_factor_to_molecular))));
    printf("======================================\n\n");
}
