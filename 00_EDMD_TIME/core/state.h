#ifndef CORE_STATE_H
#define CORE_STATE_H

#include <stdbool.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_ttf.h>
#include "kissfft/kiss_fft.h"
#include "kissfft/kiss_fftr.h"

#ifdef __cplusplus
extern "C" {
#endif

//////////// SIMULATION SETTINGS ////////////
// Global simulation flags
#define MOLECULAR_MODE 1       // 1 = reduced units; 0 = real units
#define EXPERIMENT_MODE 0      // 1 = tuned dt/substep values for better experiments
// --- Simulation time ---
#define TIME_UNITS_SIMULATION 500
#define SUBSTEPPING 1 // BOOL
#define INIT_SCALING 0 // BOOL
#define INIT_SCALING_WALL_RELEASE_STABLE 0 // BOOL

// PARTICLE PARAMETERS
#define NUM_PARTICLES 100      // total number of particles
extern bool log_packing_fraction; // 1 = log packing fraction, 0 = don't log

#if EXPERIMENT_MODE
extern bool enable_speed_of_sound_experiments;
#else
extern bool enable_speed_of_sound_experiments;
#endif

// --- Physical constants ---
#if MOLECULAR_MODE
    #define PARTICLE_RADIUS_UNIT 0.5f
    #define PARTICLE_MASS 1.0f
    #define K_B 1.0f
    #define METERS_PER_PIXEL 1.0f
#else
    #define PARTICLE_RADIUS_UNIT 1e-9f
    #define PARTICLE_MASS 1.67e-27
    #define K_B 1.38e-23f
    #define TEMPERATURE 300.0f
    #define METERS_PER_PIXEL 1e-9f
#endif

// --- Time step and substeps ---
#if MOLECULAR_MODE
    #if EXPERIMENT_MODE
        #define FIXED_DT 0.001f
        #define PIXELS_PER_SIGMA 1.0f
        #define TEMPERATURE 1.0f
        #define SUBSTEPS 20
        #define WALL_HOLD_STEPS 10000
        #define PARTICLE_SCALE_PARAM 1.0f
    #else
        #define FIXED_DT 0.0016f
        #define PIXELS_PER_SIGMA 40.0f
        #define TEMPERATURE 10000.0f
        #define SUBSTEPS 10
        #define WALL_HOLD_STEPS 10000
        #define PARTICLE_SCALE_PARAM 0.01f
    #endif
#else
    #if EXPERIMENT_MODE
        #define FIXED_DT 2e-4f
        #define SUBSTEPS 20
        #define WALL_HOLD_STEPS 10000
        #define PIXELS_PER_SIGMA 1.0f
        #define PARTICLE_SCALE_PARAM 1.0f
    #else
        #define FIXED_DT 1.6e-3f
        #define SUBSTEPS 10
        #define WALL_HOLD_STEPS 5000
        #define PIXELS_PER_SIGMA 20.0f
        #define PARTICLE_SCALE_PARAM 1.0f
    #endif
#endif

#define DIAMETER (2 * PARTICLE_RADIUS_UNIT)
#define PARTICLE_RADIUS (PARTICLE_RADIUS_UNIT * PIXELS_PER_SIGMA)

extern int num_steps;
extern bool wall_position_is_managed_externally;

// WALL PARAMETERS
#define WALL_MASS_FACTOR 50
extern float L0_UNITS;
extern float HEIGHT_UNITS;
#define SIGMA_UNIT 1.0f

extern double simulation_time;
extern double wall_release_time;

extern float wall_thickness_runtime;
#define WALL_THICKNESS (wall_thickness_runtime)
void set_wall_thickness_sigma(float thickness_in_sigma);

extern float wall_mass_runtime;
#define WALL_MASS (wall_mass_runtime)

extern bool wall_hold_enabled;
extern bool wall_is_released;
extern int wall_hold_steps;
extern int steps_elapsed;
extern float wall_impulse_x_accum;

extern double TAU_time_scale_factor_to_molecular;
void initialize_time_scale(void);

// --- Simulation dimensions ---
extern int SIM_WIDTH, SIM_HEIGHT;
extern int XW1, XW2, YW1, YW2;
extern int MAX_X, MAX_Y;
extern int HIST_WIDTH, HIST_HEIGHT;
void initialize_simulation_dimensions(void);

// Piston and wall state
extern float piston_left_x;
extern float piston_right_x;
extern float vx_piston_left;
extern float vx_piston_right;
extern double piston_work_left;
extern double piston_work_right;
extern long piston_hits_left;
extern long piston_hits_right;
extern long missed_collision_events;
extern double worst_penetration_observed;
extern float wall_x_old;
extern float wall_x;
extern float vx_wall;
extern int wall_enabled;
void initialize_simulation_parameters(void);

extern int left_count;
extern int right_count;

// SIMULATION PARAMETERS
#define MIN_SPATIAL_WINDOW ((MAX_X/10 * MAX_Y/10)/ NUM_PARTICLES)
#define NUM_BINS 400
#define MAX_VELOCITY 10000.0
#define SEED 42
#define SCALE_Vxy 1
#define Scale_Energy (1e27/1000)

// Particle states
extern double X[NUM_PARTICLES];
extern double Y[NUM_PARTICLES];
extern double X_old[NUM_PARTICLES];
extern double V_init[NUM_PARTICLES];
extern double Vx[NUM_PARTICLES];
extern double Vy[NUM_PARTICLES];
extern double Radius[NUM_PARTICLES];

extern int velocity_histogram[NUM_BINS];
extern int position_histogram[NUM_BINS * NUM_BINS];
extern int energy_histogram[NUM_BINS];

// SDL globals
extern SDL_Window* _window;
extern SDL_Renderer* _renderer;
extern TTF_Font* font;

extern int simulation_started;
extern bool paused;

// FFT
#define FFT_SIZE 256
extern bool live_fft;
extern float wall_x_samples[FFT_SIZE];
extern float energy_samples[FFT_SIZE];
extern kiss_fft_scalar in_wall[FFT_SIZE];
extern kiss_fft_scalar in_energy[FFT_SIZE];
extern kiss_fft_cpx out_wall[FFT_SIZE/2 + 1];
extern kiss_fft_cpx out_energy[FFT_SIZE/2 + 1];
extern float wall_fft_magnitude[FFT_SIZE/2 + 1];
extern float energy_fft_magnitude[FFT_SIZE/2 + 1];
extern int sample_index;
extern kiss_fftr_cfg fft_cfg;

void print_simulation_info(void);

#ifdef __cplusplus
}
#endif

#endif /* CORE_STATE_H */
