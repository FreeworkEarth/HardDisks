#include <SDL2/SDL.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <errno.h>
#include <stdlib.h>
#include <SDL2/SDL_ttf.h>
#include "kissfft/kiss_fft.h"
#include "kissfft/kiss_fftr.h"

#define MAX_DIVIDER_CAPACITY 32

static void rebuild_internal_walls(void);
static void free_internal_wall_arrays(void);
static void reset_simulation_state(void);
void initialize_simulation(void);

/*1.	Particle Simulation:
	•	Position and Velocity: Each particle has a position (X, Y) and velocity (VX, VY), which updates at each simulation step.
	•	Collisions: Particles can collide with each other, and elastic collisions should swap their velocities. Boundary conditions should also apply to reflect particles off the walls.
	•	Piston Interaction: The simulation includes two pistons that control the movement of the particles, influencing their velocity and position in a manner that simulates pressure or external forces.
	2.	Control Parameters:
	•	Time Steps (dt1, dt2): The simulation uses different time steps for particles based on certain conditions, potentially updating dynamically.
	•	Piston Movement: The pistons move at specified velocities (lvel, rvel), and the time step (dt1, dt2) changes based on these velocities.
	•	Simulation Modes: The simulation can run in different modes, which affect how the pistons and particles behave. For instance, different piston velocities, initial conditions, and time delays.
	3.	Energy and Statistics:
	•	Kinetic Energy: The simulation tracks kinetic energy, with energy calculated before and after piston movement.
	•	Histograms: The simulation tracks and updates various histograms for velocities, collision events, etc.
	4.	Boundary and Interaction Control:
	•	Reflective Boundaries: When particles hit the boundaries of the simulation box, they reflect off the walls, reversing their velocity in the corresponding direction.
	•	Collision Detection: Particles detect and resolve collisions based on their radius and distance between centers. These collisions swap their velocities for now (ideal elastic collisions).
	5.	Visualization:
	•	Display Control: A graphical display is available to visualize the particles’ movement and interactions. This is controlled by variables like dispctl, dispctl2, and histctl.

The simulation will run iteratively, with particle positions updated based on their velocities, collisions resolved, and piston interactions applied, while tracking energy and histogram data.

In summary, the simulation simulates a system of particles with piston-driven motion, boundary reflections, and collisions between particles. 
It includes control parameters for adjusting the simulation and tracking energy and statistical properties. Visualization and interaction modes can be adjusted as needed.
*/

/*  
✅ Yes, it is standard practice in molecular dynamics (MD) to set:
	•	k_B = 1
	•	m = 1
	•	\sigma = 1 (length unit — usually particle diameter)
	•	T = 1 (or whatever you want)

This is called reduced units or natural units.
It simplifies all equations.

Why?
	•	No need to carry around ugly constants (no 1.38 \times 10^{-23} all the time)
	•	Simulations become easier, more stable
	•	Scaling back to real-world units is easy afterward if needed.

🔵 Every major MD engine (like LAMMPS, GROMACS) supports “reduced units” mode.
HERE: MOLECULAR MODE
*/


//////////// SIMULATION SETTINGS ////////////
// Global simulation flags
#define MOLECULAR_MODE 1       // 1 = reduced units; 0 = real units
#define EXPERIMENT_MODE 0      // 1 = tuned dt/substep values for better experiments
// --- Simulation time ---
/* So your TIME_UNITS_SIMULATION = 3000 works in BOTH modes!
	•	In molecular mode, it means 3000 units of normalized time.
	•	In real mode, it means 3000 seconds (if you scale τ properly). */
#define TIME_UNITS_SIMULATION 500
#define SUBSTEPPING 10 // BOOL
#define INIT_SCALING 0 // BOOL
#define INIT_SCALING_WALL_RELEASE_STABLE 0 // BOOL

// PARTICLE PARAMETERS
#define NUM_PARTICLES 5000      // total number of particles
bool log_packing_fraction = 1; // 1 = log packing fraction, 0 = don't log




#if EXPERIMENT_MODE
    bool enable_speed_of_sound_experiments = true;
#else
    bool enable_speed_of_sound_experiments = false;
#endif



// --- Physical constants ---
#if MOLECULAR_MODE
    #define PARTICLE_RADIUS_UNIT 1.0f      // unit length               normalized radius = 1 (works in both modes!)
    #define PARTICLE_MASS 1.0 // Mass (m) 1 unit mass
    #define K_B 1
    #define METERS_PER_PIXEL 1.0f // pixel size for mapping to meters (important for real mode)



#else
    #define PARTICLE_RADIUS_UNIT 1e-9 // 1 nm radius
    #define PARTICLE_MASS 1.67e-27     // kg (proton mass)
    #define K_B 1.38e-23f                // J/K
    #define TEMPERATURE 300.0f           // Kelvin
    #define METERS_PER_PIXEL 1e-9 // pixel size for mapping to meters (important for real mode)

#endif

// --- Time step and substeps ---
#if MOLECULAR_MODE
    #if EXPERIMENT_MODE
        #define FIXED_DT 0.001f  // in normalized time units They don’t state the exact dt, but normalized simulations use dt = 10^{-3} often.
        #define PIXELS_PER_SIGMA 1.0f      // rendering scale
        #define TEMPERATURE 1.0f
        #define SUBSTEPS 20 // default, can be changed externally
        #define WALL_HOLD_STEPS 10000
        #define PARTICLE_SCALE_PARAM 1.0f // scaling factor for particle size to make more space for particles



    #else
        #define FIXED_DT 0.0016f  // in normalized time units They don’t state the exact dt, but normalized simulations use dt = 10^{-3} often.
        #define PIXELS_PER_SIGMA 40.0f      // rendering scaleq
        #define TEMPERATURE 300.0f            //SET to 0 if want to ssimulte traversing wave
        #define SUBSTEPS 2 // default, can be changed externally
        #define WALL_HOLD_STEPS 10000
        #define PARTICLE_SCALE_PARAM .10f // scaling factor for particle size


    #endif
#else
    #if EXPERIMENT_MODE
        #define FIXED_DT 2e-4f
        #define SUBSTEPS 20 // default, can be changed externallyq
        #define WALL_HOLD_STEPS 10000
        #define PIXELS_PER_SIGMA 1.0f      // rendering scale
        #define PARTICLE_SCALE_PARAM 1.0f // scaling factor for particle size


    #else
        #define FIXED_DT 1.6e-3f
        #define SUBSTEPS 10 // default, can be changed externally
        #define WALL_HOLD_STEPS 5000
        #define PIXELS_PER_SIGMA 20.0f      // rendering scale
        #define PARTICLE_SCALE_PARAM 1.0f // scaling factor for particle size


^
    #endif
#endif


#define DIAMETER (2 * PARTICLE_RADIUS_UNIT)
#define PARTICLE_RADIUS (PARTICLE_RADIUS_UNIT * PIXELS_PER_SIGMA * PARTICLE_SCALE_PARAM) // in meters

int num_steps = TIME_UNITS_SIMULATION / FIXED_DT;  // = 3,000,000 with 3000 units and 0.001 dt
bool wall_position_is_managed_externally = false;



// WALL PARAMETERS
#define WALL_MASS_FACTOR 200                         // choose factor relative to PARTICLE_MASS
/// SIMULATION  PARAMETERS
// Set these in main or before simulation starts
float L0_UNITS = 20.0f;       // 20 × radius = 10 × diameter Box half lenght paper 7.5 until 35
float HEIGHT_UNITS = 10.0f;   // 10 × radius = 5 × diameter
#define SIGMA_UNIT 1.0f             // radius = σ = 1

// --- GLOBAL TIMERS & WALL STATE ---
// Remove any earlier duplicate of wall_release_time and simulation_time.
double simulation_time = 0.0;         // advances only inside main loop
double wall_release_time = -1.0;      // < 0 means "not yet released"
// ...existing code...


// /////////   WALL THICKNESSS
/* 
#if EXPERIMENT_MODE
    #define WALL_THICKNESS 1.0f
#else
    #define WALL_THICKNESS fmaxf(PARTICLE_RADIUS*2, 0.10f)  // always enough to prevent tunneling //WALL_THICKNESS ((PARTICLE_RADIUS * 3 / 1.0f) < 1.0f  ? 1.0f : (PARTICLE_RADIUS * 3 / 1.0f))
#endif
*/

// Near your other runtime globals:
#if EXPERIMENT_MODE
    static float wall_thickness_runtime = 1.0f;

#else
    static float wall_thickness_runtime = 2.0f * PARTICLE_RADIUS; // default
#endif

// Helper accessors:
#define WALL_THICKNESS (wall_thickness_runtime)
// You can set this anywhere (UI key, config load, etc.):
void set_wall_thickness_sigma(float thickness_in_sigma) {
    wall_thickness_runtime = thickness_in_sigma * PIXELS_PER_SIGMA; // if your positions are in pixels
}



//////v   WEIGHTED WALL MASS
float wall_mass_runtime = PARTICLE_MASS * WALL_MASS_FACTOR;
#define WALL_MASS (wall_mass_runtime)

bool wall_hold_enabled = true;         // Enable or disable hold wall mode
bool wall_is_released = false;          // Has the user released the wall?
int wall_hold_steps = WALL_HOLD_STEPS;             // Hold the wall for this many steps
int steps_elapsed = 0;    
static float wall_impulse_x_accum = 0.f;




double TAU_time_scale_factor_to_molecular = 0.0; // Scaling factor for time in molecular mode
void initialize_time_scale() {
#if MOLECULAR_MODE
    TAU_time_scale_factor_to_molecular = 1.0; // normalized units
#else
    TAU_time_scale_factor_to_molecular = sqrtf((PARTICLE_MASS * powf(PARTICLE_RADIUS * METERS_PER_PIXEL, 2)) / (K_B * TEMPERATURE));
    printf("✅ [Time Scale] Calculated TAU = %.15e seconds\n", TAU_time_scale_factor_to_molecular);
#endif
}

// --- Simulation dimensions ---
int SIM_WIDTH, SIM_HEIGHT;
int XW1 = 50;
int XW2, YW1 = 0, YW2;
int MAX_X, MAX_Y;
int HIST_WIDTH, HIST_HEIGHT;


void initialize_simulation_dimensions() {
    SIM_WIDTH = (int)(2 * L0_UNITS * PIXELS_PER_SIGMA);
    SIM_HEIGHT = (int)(HEIGHT_UNITS * PIXELS_PER_SIGMA);

    XW2 = XW1 + SIM_WIDTH;
    YW2 = YW1 + SIM_HEIGHT;

    MAX_X = XW2 + 100;  // Add padding
    MAX_Y = YW2 + 100;

    HIST_WIDTH = MAX_X - SIM_WIDTH;
    HIST_HEIGHT = MAX_Y - SIM_HEIGHT;

    rebuild_internal_walls();
}


// Global variables for pistons
float piston_left_x = 0;
float piston_right_x;
float vx_piston_left = 0; // Slow movement for piston velocity
float vx_piston_right = 0; // Slow movement for piston velocity
static bool  piston_step_active = false;
static float piston_step_target = 0.0f;
static float piston_step_speed = 0.0f;

// piston diagnostic
static double piston_work_left  = 0.0;
static double piston_work_right = 0.0;
static long   piston_hits_left  = 0;
static long   piston_hits_right = 0;

// Diagnostics for potential collision misses (large overlaps after integration)
static long   missed_collision_events      = 0;
static double worst_penetration_observed   = 0.0;

// Divider wall configuration
static int    num_internal_walls = 1;
static int    cli_requested_walls = 1;
static int    primary_wall_index = 0;
static float *all_wall_positions = NULL;            // length = num_internal_walls
static int    extra_wall_count = 0;
static float *extra_wall_positions = NULL;          // excludes primary wall
static float *extra_wall_old_positions = NULL;
static float *extra_wall_velocity = NULL;           // mostly zero (static walls)
static bool  *extra_wall_hold_enabled = NULL;
static bool  *extra_wall_is_released = NULL;
static double *extra_wall_release_time = NULL;
static int    segment_count = 0;
static int   *segment_counts = NULL;                // length = num_internal_walls + 1
static double *segment_ke = NULL;
static double *segment_temperature = NULL;

static double primary_wall_release_time = -1.0;
static double leftmost_wall_release_time = -1.0;
static float  leftmost_wall_baseline = 0.0f;
static float  leftmost_wall_max_delta = 0.0f;
static float  leftmost_wall_current_delta = 0.0f;
static double leftmost_wall_time_of_max = 0.0;
static double leftmost_wall_average_speed = 0.0;
static float  leftmost_wall_velocity_at_max = 0.0f;
static bool   leftmost_wall_detected = false;
static bool   leftmost_wall_peak_recorded = false;
static bool   leftmost_wall_logged = false;

// -------- CLI configuration --------
static float *cli_experiment_lengths = NULL;
static size_t cli_experiment_lengths_count = 0;
static int   *cli_experiment_wall_masses = NULL;
static size_t cli_experiment_wall_masses_count = 0;
static int    cli_experiment_repeats = 1;
static int    cli_override_num_steps = -1;
static float  cli_override_L0_units = -1.0f;
static float  cli_override_height_units = -1.0f;
static float  cli_override_wall_mass_factor = -1.0f;
static float  cli_override_wall_thickness_sigma = -1.0f;
static bool   cli_force_no_experiments = false;

typedef enum {
    EXPERIMENT_PRESET_NONE = 0,
    EXPERIMENT_PRESET_OFFCENTER_SINGLE,
    EXPERIMENT_PRESET_TWO_WALLS
} ExperimentPreset;

static ExperimentPreset cli_experiment_preset = EXPERIMENT_PRESET_NONE;

static bool  preset_custom_positions = false;
static int   preset_wall_count = 0;
static float preset_wall_fraction[MAX_DIVIDER_CAPACITY];
static bool  preset_custom_fractions = false;
static float preset_segment_fraction[MAX_DIVIDER_CAPACITY + 1];

static void parse_cli_options(int argc, char **argv);
static void print_cli_usage(const char *exe_name);
// Wall parameters
float wall_x_old;
float wall_x;
float vx_wall = 0.0;  // Velocity of the wall
int wall_enabled = 0;  // 0 = Wall is disabled, 1 = Wall is enabled

/*
 * initialize_simulation_parameters
 * ---------------------------------
 * Resets the piston and divider-wall state so the interactive loop can
 * assume a clean configuration: divider centered at L0, stationary wall,
 * pistons parked at the container edges, and diagnostic counters cleared.
 */
void initialize_simulation_parameters() {
    piston_left_x = (float)XW1 - 6.0f;   // park just outside left wall
    piston_right_x = (float)XW2 + 6.0f;  // park just outside right wall
    vx_piston_left = 0.0f;
    vx_piston_right = 0.0f;
    missed_collision_events = 0;
    worst_penetration_observed = 0.0;
    wall_is_released = false;
    wall_hold_enabled = true;
    wall_impulse_x_accum = 0.0f;
    if (all_wall_positions) {
        all_wall_positions[primary_wall_index] = wall_x;
    }
}

// Global:
 // VARIABLES TO AVOID TUNNELING die to time step and stationary particles
int left_count = 0;
int right_count = 0;



// SIMULATION PARAMETERS
#define MIN_SPATIAL_WINDOW (MAX_X/10 * MAX_Y/10)/ NUM_PARTICLES// Prevent too small window
#define NUM_BINS 400 // NUM bins for histogram
#define MAX_VELOCITY 10000.0
#define SEED 42  // Seed for random number generation (ensures repeatable randomness)
#define SCALE_Vxy  1 // Mass of particle (e.g., hydrogen atom) in kg
#define Scale_Energy 1e27/1000






// Particle states
double X[NUM_PARTICLES];
double Y[NUM_PARTICLES];
double X_old[NUM_PARTICLES];    // previous X positions

double V_init[NUM_PARTICLES];   // initial velocity
double Vx[NUM_PARTICLES];
double Vy[NUM_PARTICLES];
double Radius[NUM_PARTICLES];

int velocity_histogram[NUM_BINS] = {0};
int position_histogram[NUM_BINS * NUM_BINS] = {0};
int energy_histogram[NUM_BINS] = {0};


// Global Variables
SDL_Window* _window = NULL;
SDL_Renderer* _renderer = NULL;
TTF_Font *font = NULL;

// Global variables for simulation control
int simulation_started = 0;  // Flag to track simulation start
// Global flag to pause the simulation
bool paused = false;





static void reset_simulation_state(void) {
    rebuild_internal_walls();
    initialize_simulation_parameters();
    initialize_simulation();
    simulation_time = 0.0;
    wall_release_time = -1.0;
    steps_elapsed = 0;
    simulation_started = 0;
    paused = false;
    piston_work_left = piston_work_right = 0.0;
    piston_hits_left = piston_hits_right = 0;
    missed_collision_events = 0;
    worst_penetration_observed = 0.0;
    piston_step_active = false;
    piston_step_speed = 0.0f;
    vx_piston_right = 0.0f;
    vx_piston_left = 0.0f;
    for (int j = 0; j < extra_wall_count; ++j) {
        if (extra_wall_positions) {
            extra_wall_old_positions[j] = extra_wall_positions[j];
        }
        if (extra_wall_velocity) {
            extra_wall_velocity[j] = 0.0f;
        }
        if (extra_wall_hold_enabled) {
            extra_wall_hold_enabled[j] = true;
        }
        if (extra_wall_is_released) {
            extra_wall_is_released[j] = false;
        }
        if (extra_wall_release_time) {
            extra_wall_release_time[j] = -1.0;
        }
    }
    primary_wall_release_time = -1.0;
    leftmost_wall_release_time = -1.0;
    leftmost_wall_detected = false;
    leftmost_wall_peak_recorded = false;
    leftmost_wall_logged = false;
    leftmost_wall_max_delta = 0.0f;
    leftmost_wall_average_speed = 0.0;
    leftmost_wall_velocity_at_max = 0.0f;
    if (segment_counts && segment_count > 0) {
        memset(segment_counts, 0, (size_t)segment_count * sizeof(int));
    }
}

//FFT
#define FFT_SIZE 256
bool live_fft = false;  // Flag to enable/disable live FFT
float wall_x_samples[FFT_SIZE];
float energy_samples[FFT_SIZE];
kiss_fft_scalar in_wall[FFT_SIZE];
kiss_fft_scalar in_energy[FFT_SIZE];
kiss_fft_cpx out_wall[FFT_SIZE/2 + 1];
kiss_fft_cpx out_energy[FFT_SIZE/2 + 1];
float wall_fft_magnitude[FFT_SIZE/2 + 1];
float energy_fft_magnitude[FFT_SIZE/2 + 1];
int sample_index = 0;

kiss_fftr_cfg fft_cfg;  // <== global config





// -------------------- CLI helpers --------------------
static void *cli_checked_realloc(void *ptr, size_t size) {
    void *mem = realloc(ptr, size);
    if (!mem) {
        fprintf(stderr, "Out of memory while parsing CLI options.\n");
        exit(EXIT_FAILURE);
    }
    return mem;
}

static char *cli_strdup(const char *src) {
    size_t len = strlen(src) + 1;
    char *dst = (char *)cli_checked_realloc(NULL, len);
    memcpy(dst, src, len);
    return dst;
}

static const char *cli_option_value(const char *arg, int argc, char **argv, int *index) {
    const char *eq = strchr(arg, '=');
    if (eq && *(eq + 1) != '\0') {
        return eq + 1;
    }
    if (*index + 1 >= argc) {
        fprintf(stderr, "Option '%s' expects a value.\n", arg);
        exit(EXIT_FAILURE);
    }
    return argv[++(*index)];
}

static void cli_parse_float_list(const char *value, float **buffer, size_t *count) {
    char *copy = cli_strdup(value);
    char *token = strtok(copy, ",");
    float *list = NULL;
    size_t items = 0;

    while (token) {
        errno = 0;
        char *endptr = NULL;
        float v = strtof(token, &endptr);
        if (errno != 0 || endptr == token || *endptr != '\0') {
            fprintf(stderr, "Invalid float value '%s' in list '%s'.\n", token, value);
            free(copy);
            exit(EXIT_FAILURE);
        }
        list = (float *)cli_checked_realloc(list, (items + 1) * sizeof(float));
        list[items++] = v;
        token = strtok(NULL, ",");
    }

    if (items == 0) {
        fprintf(stderr, "Empty list provided for option value '%s'.\n", value);
        free(copy);
        exit(EXIT_FAILURE);
    }

    free(*buffer);
    *buffer = list;
    *count = items;
    free(copy);
}

static void cli_parse_int_list(const char *value, int **buffer, size_t *count) {
    char *copy = cli_strdup(value);
    char *token = strtok(copy, ",");
    int *list = NULL;
    size_t items = 0;

    while (token) {
        errno = 0;
        char *endptr = NULL;
        long v = strtol(token, &endptr, 10);
        if (errno != 0 || endptr == token || *endptr != '\0') {
            fprintf(stderr, "Invalid integer value '%s' in list '%s'.\n", token, value);
            free(copy);
            exit(EXIT_FAILURE);
        }
        list = (int *)cli_checked_realloc(list, (items + 1) * sizeof(int));
        list[items++] = (int)v;
        token = strtok(NULL, ",");
    }

    if (items == 0) {
        fprintf(stderr, "Empty list provided for option value '%s'.\n", value);
        free(copy);
        exit(EXIT_FAILURE);
    }

    free(*buffer);
    *buffer = list;
    *count = items;
    free(copy);
}

static void free_internal_wall_arrays(void) {
    free(all_wall_positions); all_wall_positions = NULL;
    free(extra_wall_positions); extra_wall_positions = NULL;
    free(extra_wall_old_positions); extra_wall_old_positions = NULL;
    free(extra_wall_velocity); extra_wall_velocity = NULL;
    free(extra_wall_hold_enabled); extra_wall_hold_enabled = NULL;
    free(extra_wall_is_released); extra_wall_is_released = NULL;
    free(extra_wall_release_time); extra_wall_release_time = NULL;
    free(segment_counts); segment_counts = NULL;
    free(segment_ke); segment_ke = NULL;
    free(segment_temperature); segment_temperature = NULL;
    extra_wall_count = 0;
    segment_count = 0;
}

static void rebuild_internal_walls(void) {
    if (num_internal_walls < 1) num_internal_walls = 1;

    all_wall_positions = (float *)cli_checked_realloc(all_wall_positions,
                                                     (size_t)num_internal_walls * sizeof(float));
    float width_px = (float)SIM_WIDTH;

    if (preset_custom_positions && preset_wall_count == num_internal_walls) {
        for (int w = 0; w < num_internal_walls; ++w) {
            float frac = preset_wall_fraction[w];
            if (frac < 0.0f || frac > 1.0f) {
                frac = (float)(w + 1) / (float)(num_internal_walls + 1);
            }
            all_wall_positions[w] = XW1 + frac * width_px;
        }
    } else {
        float segment_unit_px = width_px / (float)(num_internal_walls + 1);
        for (int w = 0; w < num_internal_walls; ++w) {
            all_wall_positions[w] = XW1 + (w + 1) * segment_unit_px;
        }
    }

    float center_x = XW1 + width_px * 0.5f;
    primary_wall_index = 0;
    float best_dist = fabsf(all_wall_positions[0] - center_x);
    for (int w = 1; w < num_internal_walls; ++w) {
        float d = fabsf(all_wall_positions[w] - center_x);
        if (d < best_dist) {
            best_dist = d;
            primary_wall_index = w;
        }
    }

    wall_x = all_wall_positions[primary_wall_index];
    wall_x_old = wall_x;
    vx_wall = 0.0f;
    wall_impulse_x_accum = 0.0f;
    wall_is_released = false;
    wall_hold_enabled = true;
    primary_wall_release_time = -1.0;

    extra_wall_count = num_internal_walls - 1;
    if (extra_wall_count > 0) {
        extra_wall_positions = (float *)cli_checked_realloc(extra_wall_positions,
                                                           (size_t)extra_wall_count * sizeof(float));
        extra_wall_old_positions = (float *)cli_checked_realloc(extra_wall_old_positions,
                                                                (size_t)extra_wall_count * sizeof(float));
        extra_wall_velocity = (float *)cli_checked_realloc(extra_wall_velocity,
                                                           (size_t)extra_wall_count * sizeof(float));
        extra_wall_hold_enabled = (bool *)cli_checked_realloc(extra_wall_hold_enabled,
                                                             (size_t)extra_wall_count * sizeof(bool));
        extra_wall_is_released = (bool *)cli_checked_realloc(extra_wall_is_released,
                                                             (size_t)extra_wall_count * sizeof(bool));
        extra_wall_release_time = (double *)cli_checked_realloc(extra_wall_release_time,
                                                                (size_t)extra_wall_count * sizeof(double));
        int idx = 0;
        for (int w = 0; w < num_internal_walls; ++w) {
            if (w == primary_wall_index) continue;
            float pos = all_wall_positions[w];
            extra_wall_positions[idx] = pos;
            extra_wall_old_positions[idx] = pos;
            extra_wall_velocity[idx] = 0.0f;
            extra_wall_hold_enabled[idx] = true;
            extra_wall_is_released[idx] = false;
            extra_wall_release_time[idx] = -1.0;
            ++idx;
        }
    } else {
        free(extra_wall_positions); extra_wall_positions = NULL;
        free(extra_wall_old_positions); extra_wall_old_positions = NULL;
        free(extra_wall_velocity); extra_wall_velocity = NULL;
        free(extra_wall_hold_enabled); extra_wall_hold_enabled = NULL;
        free(extra_wall_is_released); extra_wall_is_released = NULL;
        free(extra_wall_release_time); extra_wall_release_time = NULL;
    }

    segment_count = num_internal_walls + 1;
    segment_counts = (int *)cli_checked_realloc(segment_counts, (size_t)segment_count * sizeof(int));
    segment_ke = (double *)cli_checked_realloc(segment_ke, (size_t)segment_count * sizeof(double));
    segment_temperature = (double *)cli_checked_realloc(segment_temperature, (size_t)segment_count * sizeof(double));
    memset(segment_counts, 0, (size_t)segment_count * sizeof(int));
    memset(segment_ke, 0, (size_t)segment_count * sizeof(double));
    memset(segment_temperature, 0, (size_t)segment_count * sizeof(double));

    leftmost_wall_baseline = (num_internal_walls > 0) ? all_wall_positions[0] : wall_x;
    for (int w = 1; w < num_internal_walls; ++w) {
        if (all_wall_positions[w] < leftmost_wall_baseline) {
            leftmost_wall_baseline = all_wall_positions[w];
        }
    }
    if (num_internal_walls <= 0) {
        leftmost_wall_baseline = wall_x;
    }
    leftmost_wall_max_delta = 0.0f;
    leftmost_wall_time_of_max = 0.0;
    leftmost_wall_average_speed = 0.0;
    leftmost_wall_velocity_at_max = 0.0f;
    leftmost_wall_release_time = -1.0;
    leftmost_wall_detected = false;
    leftmost_wall_peak_recorded = false;
    leftmost_wall_logged = false;
}

static void get_sorted_wall_positions(float *out_positions) {
    if (num_internal_walls <= 0) return;
    for (int i = 0; i < num_internal_walls; ++i) {
        out_positions[i] = all_wall_positions[i];
    }
    // primary wall might have moved; ensure array uses latest value
    out_positions[primary_wall_index] = wall_x;
    // simple insertion sort (wall count small)
    for (int i = 1; i < num_internal_walls; ++i) {
        float key = out_positions[i];
        int j = i - 1;
        while (j >= 0 && out_positions[j] > key) {
            out_positions[j + 1] = out_positions[j];
            --j;
        }
        out_positions[j + 1] = key;
    }
}

static void sync_all_wall_positions(void) {
    if (!all_wall_positions || num_internal_walls <= 0) return;
    int idx = 0;
    for (int w = 0; w < num_internal_walls; ++w) {
        if (w == primary_wall_index) {
            all_wall_positions[w] = wall_x;
        } else {
            if (idx < extra_wall_count && extra_wall_positions) {
                all_wall_positions[w] = extra_wall_positions[idx++];
            }
        }
    }
}

static void reset_leftmost_metrics(double release_time) {
    sync_all_wall_positions();
    float min_pos = wall_x;
    if (extra_wall_count > 0 && extra_wall_positions) {
        for (int j = 0; j < extra_wall_count; ++j) {
            if (extra_wall_positions[j] < min_pos) {
                min_pos = extra_wall_positions[j];
            }
        }
    }
    leftmost_wall_baseline = min_pos;
    leftmost_wall_max_delta = 0.0f;
    leftmost_wall_current_delta = 0.0f;
    leftmost_wall_time_of_max = release_time;
    leftmost_wall_average_speed = 0.0;
    leftmost_wall_velocity_at_max = 0.0f;
    leftmost_wall_release_time = release_time;
    leftmost_wall_detected = false;
    leftmost_wall_peak_recorded = false;
    leftmost_wall_logged = false;
}

static int segment_index_for_position(double px) {
    if (segment_count <= 0) return 0;
    if (num_internal_walls <= 0) {
        return (px < XW1) ? 0 : (px > XW2 ? 0 : 0);
    }
    float tmp_positions[MAX_DIVIDER_CAPACITY];
    if (num_internal_walls > (int)(sizeof(tmp_positions)/sizeof(tmp_positions[0]))) {
        // fallback without sorting beyond capacity
        return 0;
    }
    get_sorted_wall_positions(tmp_positions);
    int seg = 0;
    for (int w = 0; w < num_internal_walls; ++w) {
        float boundary = tmp_positions[w];
        if (px < boundary) {
            return seg;
        }
        seg++;
    }
    if (seg >= segment_count) seg = segment_count - 1;
    return seg;
}

static void compute_segment_bounds(float *left_bounds, float *right_bounds) {
    int segs = num_internal_walls + 1;
    if (segs <= 0) return;
    float tmp[MAX_DIVIDER_CAPACITY];
    if (num_internal_walls > (int)(sizeof(tmp)/sizeof(tmp[0]))) {
        for (int i = 0; i < segs; ++i) {
            left_bounds[i] = (float)XW1;
            right_bounds[i] = (float)XW2;
        }
        return;
    }
    get_sorted_wall_positions(tmp);
    float prev = (float)XW1;
    for (int s = 0; s < num_internal_walls; ++s) {
        float boundary = tmp[s];
        left_bounds[s] = prev;
        right_bounds[s] = boundary - WALL_THICKNESS * 0.5f;
        prev = boundary + WALL_THICKNESS * 0.5f;
    }
    left_bounds[segs - 1] = prev;
    right_bounds[segs - 1] = (float)XW2;
}

static void print_cli_usage(const char *exe_name) {
    printf("Usage: %s [options]\n", exe_name);
    printf("Options:\n");
    printf("  --experiments               Run the speed-of-sound experiment batch\n");
    printf("  --no-experiments            Force interactive SDL mode even if experiments are enabled\n");
    printf("  --walls=N                   Number of internal divider walls (default 1)\n");
    printf("  --experiment=NAME           Use a canned experiment preset (2wall_exp_with_1_wall, 2wall_exp_with_2_wall)\n");
    printf("  --lengths=L0,...            Override experiment half-lengths (sigma units)\n");
    printf("  --wall-masses=f1,...        Override experiment wall-mass factors (multiples of particle mass)\n");
    printf("  --repeats=N                 Number of repeats per experiment configuration\n");
    printf("  --steps=N                   Steps to simulate after wall release (default derived from TIME_UNITS_SIMULATION)\n");
    printf("  --l0=value                  Set initial half-length for interactive mode\n");
    printf("  --height=value              Set box height (in sigma units) for interactive mode\n");
    printf("  --wall-mass-factor=value    Set wall mass factor for interactive mode (M_wall = factor * m_particle)\n");
    printf("  --wall-thickness=value      Set wall thickness in sigma units\n");
    printf("  --help                      Show this help message and exit\n");

}

static void parse_cli_options(int argc, char **argv) {
    for (int i = 1; i < argc; ++i) {
        const char *arg = argv[i];
        

        if (strcmp(arg, "--help") == 0) {
            print_cli_usage(argv[0]);
            exit(EXIT_SUCCESS);
        } else if (strcmp(arg, "--experiments") == 0) {
            enable_speed_of_sound_experiments = true;
            cli_force_no_experiments = false;
        } else if (strcmp(arg, "--no-experiments") == 0) {
            cli_force_no_experiments = true;
        } else if (strncmp(arg, "--experiment", 12) == 0) {
            const char *value_raw = cli_option_value(arg, argc, argv, &i);
            const char *variant_str = NULL;
            char preset_name[128];
            size_t name_len = strcspn(value_raw, ":");
            if (name_len >= sizeof(preset_name)) {
                name_len = sizeof(preset_name) - 1;
            }
            memcpy(preset_name, value_raw, name_len);
            preset_name[name_len] = '\0';
            if (value_raw[name_len] == ':') {
                variant_str = value_raw + name_len + 1;
            }

            if ((!variant_str || *variant_str == '\0') && (i + 1) < argc) {
                const char *candidate = argv[i + 1];
                if (candidate[0] != '-') {
                    variant_str = candidate;
                    ++i;
                }
            }

            if (strcmp(preset_name, "2wall_exp_with_1_wall") == 0) {
                cli_experiment_preset = EXPERIMENT_PRESET_OFFCENTER_SINGLE;
                enable_speed_of_sound_experiments = true;
                cli_force_no_experiments = false;
            } else if (strcmp(preset_name, "2wall_exp_with_2_wall") == 0) {
                preset_two_wall_active_walls = 2;
                if (variant_str && *variant_str) {
                    errno = 0;
                    char *endptr = NULL;
                    long v = strtol(variant_str, &endptr, 10);
                    if (errno != 0 || endptr == variant_str || *endptr != '\0' || v < 1 || v > 2) {
                        fprintf(stderr, "Invalid variant '%s' for preset '%s'. Use 1 or 2.\n",
                                variant_str, preset_name);
                        exit(EXIT_FAILURE);
                    }
                    preset_two_wall_active_walls = (int)v;
                }
                cli_experiment_preset = EXPERIMENT_PRESET_TWO_WALLS;
                enable_speed_of_sound_experiments = true;
                cli_force_no_experiments = false;
            } else {
                fprintf(stderr, "Unknown experiment preset '%s'.\n", value_raw);
                print_cli_usage(argv[0]);
                exit(EXIT_FAILURE);
            }

        } else if (strncmp(arg, "--walls", 7) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0;
            char *endptr = NULL;
            long v = strtol(value, &endptr, 10);
            if (errno != 0 || endptr == value || *endptr != '\0' || v < 1) {
                fprintf(stderr, "Invalid wall count '%s'.\n", value);
                exit(EXIT_FAILURE);
            }
            if (v > MAX_DIVIDER_CAPACITY) {
                printf("⚠️ Requested wall count %ld exceeds capacity %d. Clamping.\n", v, MAX_DIVIDER_CAPACITY);
                v = MAX_DIVIDER_CAPACITY;
            }
            cli_requested_walls = (int)v;
        } else if (strncmp(arg, "--lengths", 9) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            cli_parse_float_list(value, &cli_experiment_lengths, &cli_experiment_lengths_count);
            enable_speed_of_sound_experiments = true;
        } else if (strncmp(arg, "--wall-masses", 13) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            cli_parse_int_list(value, &cli_experiment_wall_masses, &cli_experiment_wall_masses_count);
            enable_speed_of_sound_experiments = true;
        } else if (strncmp(arg, "--repeats", 9) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0;
            char *endptr = NULL;
            long v = strtol(value, &endptr, 10);
            if (errno != 0 || endptr == value || *endptr != '\0' || v <= 0) {
                fprintf(stderr, "Invalid repeats value '%s'.\n", value);
                exit(EXIT_FAILURE);
            }
            cli_experiment_repeats = (int)v;
            enable_speed_of_sound_experiments = true;
        } else if (strncmp(arg, "--steps", 7) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0;
            char *endptr = NULL;
            long v = strtol(value, &endptr, 10);
            if (errno != 0 || endptr == value || *endptr != '\0' || v <= 0) {
                fprintf(stderr, "Invalid steps value '%s'.\n", value);
                exit(EXIT_FAILURE);
            }
            cli_override_num_steps = (int)v;
            enable_speed_of_sound_experiments = true;
        } else if (strncmp(arg, "--l0", 4) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0;
            char *endptr = NULL;
            float v = strtof(value, &endptr);
            if (errno != 0 || endptr == value || *endptr != '\0' || v <= 0.0f) {
                fprintf(stderr, "Invalid L0 value '%s'.\n", value);
                exit(EXIT_FAILURE);
            }
            cli_override_L0_units = v;
        } else if (strncmp(arg, "--height", 8) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0;
            char *endptr = NULL;
            float v = strtof(value, &endptr);
            if (errno != 0 || endptr == value || *endptr != '\0' || v <= 0.0f) {
                fprintf(stderr, "Invalid height value '%s'.\n", value);
                exit(EXIT_FAILURE);
            }
            cli_override_height_units = v;
        } else if (strncmp(arg, "--wall-mass-factor", 19) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0;
            char *endptr = NULL;
            float v = strtof(value, &endptr);
            if (errno != 0 || endptr == value || *endptr != '\0' || v <= 0.0f) {
                fprintf(stderr, "Invalid wall mass factor '%s'.\n", value);
                exit(EXIT_FAILURE);
            }
            cli_override_wall_mass_factor = v;
        } else if (strncmp(arg, "--wall-thickness", 16) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0;
            char *endptr = NULL;
            float v = strtof(value, &endptr);
            if (errno != 0 || endptr == value || *endptr != '\0' || v <= 0.0f) {
                fprintf(stderr, "Invalid wall thickness '%s'.\n", value);
                exit(EXIT_FAILURE);
            }
            cli_override_wall_thickness_sigma = v;
        } else {
            fprintf(stderr, "Unknown option '%s'.\n", arg);
            print_cli_usage(argv[0]);
            exit(EXIT_FAILURE);
        }
    }

    switch (cli_experiment_preset) {
        case EXPERIMENT_PRESET_OFFCENTER_SINGLE:
            cli_requested_walls = 1;
            preset_custom_positions = true;
            preset_wall_count = 1;
            preset_wall_fraction[0] = 1.0f / 3.0f;
            preset_custom_fractions = true;
            preset_segment_fraction[0] = 1.0f / 3.0f;
            preset_segment_fraction[1] = 2.0f / 3.0f;
            break;
        case EXPERIMENT_PRESET_TWO_WALLS:
            preset_custom_positions = true;
            preset_custom_fractions = true;
            memset(preset_wall_fraction, 0, sizeof(preset_wall_fraction));
            memset(preset_segment_fraction, 0, sizeof(preset_segment_fraction));
            if (preset_two_wall_active_walls <= 1) {
                cli_requested_walls = 1;
                preset_wall_count = 1;
                preset_wall_fraction[0] = 1.0f / 3.0f;
                preset_segment_fraction[0] = 1.0f / 3.0f;
                preset_segment_fraction[1] = 2.0f / 3.0f;
            } else {
                cli_requested_walls = 2;
                preset_wall_count = 2;
                preset_wall_fraction[0] = 1.0f / 3.0f;
                preset_wall_fraction[1] = 2.0f / 3.0f;
                preset_segment_fraction[0] = 1.0f / 3.0f;
                preset_segment_fraction[1] = 1.0f / 3.0f;
                preset_segment_fraction[2] = 1.0f / 3.0f;
            }
            break;
        case EXPERIMENT_PRESET_NONE:
        default:
            break;
    }

    num_internal_walls = cli_requested_walls;
    if (num_internal_walls > MAX_DIVIDER_CAPACITY) {
        num_internal_walls = MAX_DIVIDER_CAPACITY;
    }

    if (cli_override_num_steps > 0) {
        num_steps = cli_override_num_steps;
    }
    if (cli_override_L0_units > 0.0f) {
        L0_UNITS = cli_override_L0_units;
    }
    if (cli_override_height_units > 0.0f) {
        HEIGHT_UNITS = cli_override_height_units;
    }
    if (cli_override_wall_mass_factor > 0.0f) {
        wall_mass_runtime = PARTICLE_MASS * cli_override_wall_mass_factor;
    }
    if (cli_override_wall_thickness_sigma > 0.0f) {
        set_wall_thickness_sigma(cli_override_wall_thickness_sigma);
    }
    if (cli_force_no_experiments) {
        enable_speed_of_sound_experiments = false;
    }
    if (cli_experiment_repeats < 1) {
        cli_experiment_repeats = 1;
    }

    if (!preset_custom_fractions) {
        memset(preset_segment_fraction, 0, sizeof(preset_segment_fraction));
    }
}




///////////  SIMULALTION INFO ////////////
void print_simulation_info() {
    printf("\n========== Simulation Info ==========\n");
    printf("Mode: %s\n", MOLECULAR_MODE ? "MOLECULAR (Normalized units)" : "MACROSCOPIC (Real units)");
    printf("Experiment Mode: %s\n", EXPERIMENT_MODE ? "ON (fine dt/substeps)" : "OFF (default dt)");
    printf("Fixed Time Step (dt): %.2e\n", FIXED_DT);
    printf("TAU scaling factor: %.4f\n", TAU_time_scale_factor_to_molecular);
    printf("Wall hold mode: %s\n", wall_hold_enabled ? "Enabled" : "Disabled");
    printf("Wall hold steps: %d\n", wall_hold_steps);
    printf("Internal walls: %d\n", num_internal_walls);
    printf("Target total simulation time: %.1f\n", 3000.0f);
    printf("Calculated steps to reach target: %d\n", (int)(3000.0f / (FIXED_DT / (MOLECULAR_MODE ? 1.0f : TAU_time_scale_factor_to_molecular))));
    printf("======================================\n\n");
}

//////////////// SDL INITIALIZATION //////////////////////
// Function to Initialize SDL
void initSDL() {
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        printf("SDL could not initialize! SDL_Error: %s\n", SDL_GetError());
        exit(1);
    }
    _window = SDL_CreateWindow("Hard Spheres Simulation", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, MAX_X, MAX_Y, SDL_WINDOW_SHOWN);
    if (!_window) {
        printf("Window could not be created! SDL_Error: %s\n", SDL_GetError());
        exit(1);
    }
    _renderer = SDL_CreateRenderer(_window, -1, SDL_RENDERER_ACCELERATED);
    if (!_renderer) {
        printf("Renderer could not be created! SDL_Error: %s\n", SDL_GetError());
        exit(1);
    }
}





////// INITIALIZATION FUNCTIONS //////
// Each component of velocity (v_x, v_y) in thermal equilibrium is Gaussian:
/* Function to generate Maxwell-Boltzmann distributed velocity 2D (v_x ,y, z are normal distributed separately!!!!!)


// f(v_x) = sqrt(m / (2π k_B T)) * exp(-m v_x^2 / (2 k_B T))
// Mean = 0, StdDev = sqrt(k_B T / m)

// This explains why histograms of v_x and v_y are symmetric and Gaussian.

// Speed v = sqrt(v_x^2 + v_y^2) in 2D follows a Rayleigh distribution:

// f(v) = (m v / k_B T) * exp(-m v^2 / (2 k_B T))
// Not symmetric: peaks at v > 0, long tail

// In 3D, speed follows the Maxwell-Boltzmann distribution:

// f(v) ∝ v^2 * exp(-m v^2 / (2 k_B T))

// Box-Muller transform can generate Gaussian (normal) random numbers for v_x, v_y

*/

////////////////// no drift due to statistical sampling from distribution == only thermal , kinetic energy
void remove_drift() {
    float vx_sum = 0.0f, vy_sum = 0.0f;

    for (int i = 0; i < NUM_PARTICLES; i++) {
        vx_sum += Vx[i];
        vy_sum += Vy[i];
    }

    float vx_mean = vx_sum / NUM_PARTICLES;
    float vy_mean = vy_sum / NUM_PARTICLES;


    for (int i = 0; i < NUM_PARTICLES; i++) {
        Vx[i] -= vx_mean;
        Vy[i] -= vy_mean;
    }

    printf("✅ Drift removed: mean Vx = %.4f, mean Vy = %.4f\n", vx_mean, vy_mean);
}







///////////////////////////   ENERGY CALCULATIONS   ///////////////////////
double kinetic_energy() {
    double ke_total = 0.0;
    for (int i = 0; i < NUM_PARTICLES; i++) {
        double particle_ke = 0.5 * PARTICLE_MASS * (Vx[i] * Vx[i] + Vy[i] * Vy[i]);  // In Joules
        ke_total += particle_ke;
    }
    return ke_total;

}

float compute_measured_temperature_from_ke() {
    // USES RELATION:
    // T = \frac{2}{d \cdot N k_B} \sum \frac{1}{2}mv^2 = \frac{E_{\text{kin}}}{N k_B} \quad \text{(for 2D)}
    double ke_total = kinetic_energy();  // reuse existing function
    return (float)(ke_total / (NUM_PARTICLES * K_B));
}



// Returns a 2D Maxwell-Boltzmann distributed velocity vector (vx, vy)
void maxwell_boltzmann_2D(float temperature, float *vx, float *vy) {
    float sigma = sqrt(K_B * temperature / PARTICLE_MASS);
    float u1, u2, z0, z1;

    do {
        u1 = (float)rand() / RAND_MAX;
        u2 = (float)rand() / RAND_MAX;
    } while (u1 <= 1e-10f); // prevent log(0)

    z0 = sqrtf(-2.0f * logf(u1)) * cosf(2.0f * M_PI * u2);
    z1 = sqrtf(-2.0f * logf(u1)) * sinf(2.0f * M_PI * u2);

    *vx = sigma * z0;
    *vy = sigma * z1;
}




// Function to generate Maxwell-Boltzmann distributed velocity in 1D
float maxwell_boltzmann_velocity_gaussians(float temperature) {
    float sigma = sqrt(K_B * temperature / PARTICLE_MASS);  // Standard deviation (velocity scale)
    printf("sigma = %f\n", sigma);
    float u1, u2, z;

    // Box-Muller transform: Generate Gaussian-distributed value
    /*
    The Box-Muller transform converts two uniformly distributed random numbers into two independent standard normally distributed numbers. It works as follows:
	1.	Draw two uniform random numbers  u_1, u_2  from  U(0,1) .
	2.	Transform them using: z = \sqrt{-2 \ln u_1} \cos(2\pi u_2)
    */
    do {
        u1 = (float)rand() / RAND_MAX;        // Generate two independent random numbers in the range (0,1)
        u2 = (float)rand() / RAND_MAX;
        z = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);  // Generate a normal-distributed variable Box-Muller transform equations: Z = \sqrt{-2 \ln u_1} \cos(2\pi u_2) where  Z \sim N(0,1)  (a standard normal variable).
        printf("u1 = %f, u2 = %f, z = %f\n", u1, u2, z);
    } while (isnan(z));  // Ensure numerical stability
    return sigma * z ;  // To scale this to our required normal distribution  N(0, \sigma^2) , we multiply by  \sigma : v_x = \sigma Z = \sqrt{\frac{k_B T}{m}} \cdot Z Scale the Gaussian value by sigma to match Maxwell-Boltzmann distribution
    /*
    •	Box-Muller Transform:
        Marsaglia, G., & Bray, T. A. (1964). “A convenient method for generating normal variables.” SIAM Review, 6(3), 260–264.
	•	Maxwell-Boltzmann Distribution:
        Reif, F. (1965). Fundamentals of Statistical and Thermal Physics. McGraw-Hill.
    */
}






// NEW: Generate one Gaussian velocity sample
/*
*/
float sample_gaussian_velocity(float temperature) {
    float sigma = sqrt(K_B * temperature / PARTICLE_MASS);
    float u1, u2, z;

    do {
        u1 = (float)rand() / RAND_MAX;
        u2 = (float)rand() / RAND_MAX;
        z = sqrtf(-2.0f * logf(u1)) * cosf(2.0f * M_PI * u2);
    } while (isnan(z));

    return sigma * z;
}




void maxwell_boltzmann_velocity_ARRAY_ALL(float temperature) {
    float sigma = sqrt(K_B * temperature / PARTICLE_MASS);

    for (int i = 0; i < NUM_PARTICLES; i += 2) {
        float u1 = (float)rand() / RAND_MAX;
        float u2 = (float)rand() / RAND_MAX;

        float z0 = sqrt(-2.0f * log(u1)) * cos(2.0f * M_PI * u2);
        float z1 = sqrt(-2.0f * log(u1)) * sin(2.0f * M_PI * u2);

        Vx[i] = sigma * z0;
        Vy[i] = sigma * z1;

        // If odd number of particles:
        if (i + 1 < NUM_PARTICLES) {
            float u3 = (float)rand() / RAND_MAX;
            float u4 = (float)rand() / RAND_MAX;
            float z2 = sqrt(-2.0f * log(u3)) * cos(2.0f * M_PI * u4);
            float z3 = sqrt(-2.0f * log(u3)) * sin(2.0f * M_PI * u4);
            Vx[i + 1] = sigma * z2;
            Vy[i + 1] = sigma * z3;
        }
    }
}






  
//////////////////   in 3D for v complete: ////////////////////////////
// f(v) = \left(\frac{m}{2\pi k_B T}\right)^{3/2} 4\pi v^2 \exp\left(-\frac{m v^2}{2 k_B T}\right)
// Function to generate speed v sampled from the Maxwell-Boltzmann distribution
float maxwell_boltzmann_speed3D(float temperature) {
    float sigma = sqrt(K_B * temperature / PARTICLE_MASS);  // Compute sigma
    printf("sigma = %f\n", sigma);

    float u1, u2, u3, z1, z2, z3;

    // Generate three independent standard normal variables using Box-Muller
    do {
        u1 = (float)rand() / RAND_MAX;
        u2 = (float)rand() / RAND_MAX;
        u3 = (float)rand() / RAND_MAX;

        z1 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
        z2 = sqrt(-2.0 * log(u1)) * sin(2.0 * M_PI * u2);
        z3 = sqrt(-2.0 * log(u3)) * cos(2.0 * M_PI * u3);
        
        printf("u1 = %f, u2 = %f, u3 = %f, z1 = %f, z2 = %f, z3 = %f\n", u1, u2, u3, z1, z2, z3);

    } while (isnan(z1) || isnan(z2) || isnan(z3));  // Ensure numerical stability

    // Compute speed from the three normal components
    float v = sigma * sqrt(z1 * z1 + z2 * z2 + z3 * z3);
        // Generate random direction on a sphere
        float r1 = (float)rand() / RAND_MAX;
        float r2 = (float)rand() / RAND_MAX;
        float theta = acos(2.0 * r1 - 1.0);  // Correctly distributed polar angle
        float phi = 2.0 * M_PI * r2;         // Uniform azimuthal angle

        // Convert to Cartesian components
        float vx = v * sin(theta) * cos(phi);
        float vy = v * sin(theta) * sin(phi);
        float vz = v * cos(theta);

    return v, vx, vy, vz;  // Return the speed
}




void initialize_simulation_random() {
    int segments = num_internal_walls + 1;
    if (segments <= 0) segments = 1;

    float seg_left[MAX_DIVIDER_CAPACITY + 1];
    float seg_right[MAX_DIVIDER_CAPACITY + 1];

    if (segments > MAX_DIVIDER_CAPACITY + 1) {
        segments = MAX_DIVIDER_CAPACITY + 1;
    }

    if (num_internal_walls > 0) {
        compute_segment_bounds(seg_left, seg_right);
    } else {
        seg_left[0] = (float)XW1;
        seg_right[0] = (float)XW2;
    }

    int base = NUM_PARTICLES / segments;
    int remainder = NUM_PARTICLES % segments;
    int target_counts[MAX_DIVIDER_CAPACITY + 1];

    if (preset_custom_fractions && preset_wall_count + 1 == segments) {
        double accum = 0.0;
        int assigned = 0;
        for (int seg = 0; seg < segments - 1; ++seg) {
            double frac = preset_segment_fraction[seg];
            if (frac < 0.0) frac = 0.0;
            accum += frac;
            int c = (int)llround(frac * NUM_PARTICLES);
            if (c < 0) c = 0;
            if (c > NUM_PARTICLES) c = NUM_PARTICLES;
            target_counts[seg] = c;
            assigned += c;
        }
        int last = NUM_PARTICLES - assigned;
        if (last < 0) last = 0;
        target_counts[segments - 1] = last;
        for (int seg = 0; seg < segments; ++seg) {
            if (target_counts[seg] < 0) target_counts[seg] = 0;
        }
    } else {
        for (int seg = 0; seg < segments; ++seg) {
            target_counts[seg] = base + (seg < remainder ? 1 : 0);
        }
    }

    int index = 0;
    for (int seg = 0; seg < segments; ++seg) {
        int target = target_counts[seg];
        float left_bound = seg_left[seg] + PARTICLE_RADIUS;
        float right_bound = seg_right[seg] - PARTICLE_RADIUS;
        if (right_bound <= left_bound) {
            right_bound = left_bound + PARTICLE_RADIUS;
        }

        for (int c = 0; c < target && index < NUM_PARTICLES; ++c, ++index) {
            int valid_position = 0;
            int attempts = 0;

            while (!valid_position) {
                float rx = (float)rand() / (float)RAND_MAX;
                X[index] = left_bound + rx * (right_bound - left_bound);
                float ry = (float)rand() / (float)RAND_MAX;
                Y[index] = YW1 + PARTICLE_RADIUS + ry * ((YW2 - YW1) - 2.0f * PARTICLE_RADIUS);

                valid_position = 1;
                for (int j = 0; j < index; ++j) {
                    float dx = X[index] - X[j];
                    float dy = Y[index] - Y[j];
                    float dist_sq = dx * dx + dy * dy;
                    if (dist_sq < (4.0f * PARTICLE_RADIUS * PARTICLE_RADIUS)) {
                        valid_position = 0;
                        break;
                    }
                }

                attempts++;
                if (attempts > 20000) {
                    printf("❌ Could not place particle %d in segment %d after many tries!\n", index, seg);
                    exit(1);
                }
            }

            Vx[index] = maxwell_boltzmann_velocity_gaussians(TEMPERATURE);
            Vy[index] = maxwell_boltzmann_velocity_gaussians(TEMPERATURE);
            Radius[index] = PARTICLE_RADIUS;
        }
    }

    while (index < NUM_PARTICLES) {
        Vx[index] = maxwell_boltzmann_velocity_gaussians(TEMPERATURE);
        Vy[index] = maxwell_boltzmann_velocity_gaussians(TEMPERATURE);
        Radius[index] = PARTICLE_RADIUS;
        ++index;
    }

    if (all_wall_positions) {
        all_wall_positions[primary_wall_index] = wall_x;
    }
}





/// place particles on grid honeycomb, works
void initialize_simulation_honeycomb() {
    if (num_internal_walls > 1) {
        initialize_simulation_random();
        return;
    }
    float wall_buffer = WALL_THICKNESS + PARTICLE_RADIUS * 0.20f;  // space around wall
    float small_gap = PARTICLE_RADIUS * 0.01f;  // just a tiny gap

    float left_width = wall_x - wall_buffer - XW1;
    float right_width = XW2 - wall_x - wall_buffer;
    float height = YW2 - YW1;

    int particles_left = NUM_PARTICLES / 2;
    if (preset_custom_fractions && segment_count >= 2) {
        particles_left = (int)llround(preset_segment_fraction[0] * NUM_PARTICLES);
        if (particles_left < 0) particles_left = 0;
        if (particles_left > NUM_PARTICLES) particles_left = NUM_PARTICLES;
    }
    int particles_right = NUM_PARTICLES - particles_left;

    // 🧠💡 Plan: determine how many rows and columns we need
    int num_rows_left = (int)sqrtf(particles_left);
    int num_cols_left = (particles_left + num_rows_left - 1) / num_rows_left;  // ceil division

    int num_rows_right = (int)sqrtf(particles_right);
    int num_cols_right = (particles_right + num_rows_right - 1) / num_rows_right;

    float horizontal_spacing_left = left_width / num_cols_left;
    float vertical_spacing_left = height / num_rows_left;

    float horizontal_spacing_right = right_width / num_cols_right;
    float vertical_spacing_right = height / num_rows_right;

    int idx = 0;

    // --- Place Left Particles ---
    for (int row = 0; row < num_rows_left; row++) {
        for (int col = 0; col < num_cols_left; col++) {
            if (idx >= particles_left) break;

            float x = XW1 + (col + 0.5f) * horizontal_spacing_left;
            float y = YW1 + (row + 0.5f) * vertical_spacing_left;

            X[idx] = x;
            Y[idx] = y;
            Vx[idx] = maxwell_boltzmann_velocity_gaussians(TEMPERATURE);
            Vy[idx] = maxwell_boltzmann_velocity_gaussians(TEMPERATURE);
            Radius[idx] = PARTICLE_RADIUS;

            printf("✅ Particle %d (LEFT) placed at (%.2f, %.2f)\n", idx, x, y);
            idx++;
        }
    }

    // --- Place Right Particles ---
    for (int row = 0; row < num_rows_right; row++) {
        for (int col = 0; col < num_cols_right; col++) {
            if (idx >= NUM_PARTICLES) break;

            float x = wall_x + wall_buffer + (col + 0.5f) * horizontal_spacing_right;
            float y = YW1 + (row + 0.5f) * vertical_spacing_right;

            X[idx] = x;
            Y[idx] = y;
            Vx[idx] = maxwell_boltzmann_velocity_gaussians(TEMPERATURE);
            Vy[idx] = maxwell_boltzmann_velocity_gaussians(TEMPERATURE);
            Radius[idx] = PARTICLE_RADIUS;

            printf("✅ Particle %d (RIGHT) placed at (%.2f, %.2f)\n", idx, x, y);
            idx++;
        }
    }

    printf("✅ Finished placing all %d particles evenly across box!\n", NUM_PARTICLES);
    printf("🔍 Initial particle velocities (first 10):\n");
    for (int i = 0; i < 10; i++) {
        printf("Particle %2d: Vx = %.6f, Vy = %.6f\n", i, Vx[i], Vy[i]);    
    }
}





void initialize_simulation(void) {
    if (num_internal_walls > 1) {
        initialize_simulation_random();
        return;
    }
    float wall_buffer = WALL_THICKNESS + DIAMETER * 0.25f; // space from wall
    float margin = DIAMETER * 1.25f;

    float left_width = wall_x - wall_buffer - XW1;
    float right_width = XW2 - wall_x - wall_buffer;
    float height = YW2 - YW1;

    int particles_left = NUM_PARTICLES / 2;
    if (preset_custom_fractions && segment_count >= 2) {
        particles_left = (int)llround(preset_segment_fraction[0] * NUM_PARTICLES);
        if (particles_left < 0) particles_left = 0;
        if (particles_left > NUM_PARTICLES) particles_left = NUM_PARTICLES;
    }
    int particles_right = NUM_PARTICLES - particles_left;

    int num_rows_left = (int)sqrt(particles_left);
    int num_cols_left = (particles_left + num_rows_left - 1) / num_rows_left;

    int num_rows_right = (int)sqrt(particles_right);
    int num_cols_right = (particles_right + num_rows_right - 1) / num_rows_right;

    float spacing_left_x = left_width / num_cols_left;
    float spacing_left_y = height / num_rows_left;

    float spacing_right_x = right_width / num_cols_right;
    float spacing_right_y = height / num_rows_right;

    int idx = 0;

    // Left half
    for (int row = 0; row < num_rows_left; row++) {
        for (int col = 0; col < num_cols_left; col++) {
            if (idx >= particles_left) break;
            X[idx] = XW1 + margin + col * spacing_left_x;
            Y[idx] = YW1 + margin + row * spacing_left_y;

            float vx, vy;
            maxwell_boltzmann_2D(TEMPERATURE, &vx, &vy);
            Vx[idx] = vx;
            Vy[idx] = vy;
            V_init[idx] = sqrtf(vx * vx + vy * vy);  // Optional: store the speed if needed
            Radius[idx] = PARTICLE_RADIUS;
            idx++;
        }
    }

    // Right half
    for (int row = 0; row < num_rows_right; row++) {
        for (int col = 0; col < num_cols_right; col++) {
            if (idx >= NUM_PARTICLES) break;
            X[idx] = XW2 - margin - col * spacing_right_x;
            Y[idx] = YW1 + margin + row * spacing_right_y;
            float vx, vy;
            maxwell_boltzmann_2D(TEMPERATURE, &vx, &vy);
            Vx[idx] = vx;
            Vy[idx] = vy;
            V_init[idx] = sqrtf(vx * vx + vy * vy);  // Optional: store the speed if needed
            Radius[idx] = PARTICLE_RADIUS;
            idx++;
        }
    }

    printf("✅ Particles placed with buffer. Wall_x = %.3f\n", wall_x);
    fflush(stdout);
    // After initializing all Vx, Vy:
    double actual_ke = kinetic_energy();  // sum of 0.5*m*v^2
    double target_ke = NUM_PARTICLES * K_B * TEMPERATURE;
    double scale = sqrt(target_ke / actual_ke);  // to adjust v_i → v_i * scale
    printf("✅ Scaling factor for velocities: %.4f\n", scale);

    for (int i = 0; i < NUM_PARTICLES; i++) {
        Vx[i] *= scale;
        Vy[i] *= scale;
}



}






/////////////////// PISTON FUNCTIONS /////////////////////////

// Function to render pistons
void render_pistons() {
    SDL_SetRenderDrawColor(_renderer, 0, 255, 0, 150);
    SDL_Rect left_piston = { (int)piston_left_x , 0, 5, YW2 };
    SDL_RenderFillRect(_renderer, &left_piston);
    SDL_Rect right_piston = { (int)piston_right_x , 0, 5, YW2 };
    SDL_RenderFillRect(_renderer, &right_piston);
}



void move_left_piston(float dt) {
    float left_box_edge = (float)XW1;
    float piston_width = 20.0f;
    float right_limit = piston_right_x - piston_width; // Ensure pistons don't overlap

    piston_left_x += vx_piston_left * dt;

    if (piston_left_x < left_box_edge) {
        piston_left_x = left_box_edge;
        vx_piston_left = 0;
    }
    if (piston_left_x > right_limit) {
        piston_left_x = right_limit;
        vx_piston_left = 0;
    }
}


void move_right_piston(float dt) {
    float right_box_edge = (float)XW2;
    float left_limit = piston_left_x + 20.0f; // Ensure pistons don't overlap

    piston_right_x += vx_piston_right * dt;

    if (piston_step_active && piston_right_x <= piston_step_target) {
        piston_right_x = piston_step_target;
        vx_piston_right = 0.0f;
        piston_step_active = false;
    }

    if (piston_right_x > right_box_edge) {
        piston_right_x = right_box_edge;
        vx_piston_right = 0;
    }
    if (piston_right_x < left_limit) {
        piston_right_x = left_limit;
        vx_piston_right = 0;
    }
}

// Function to update piston positions based on velocity and time step
void update_pistons(float dt) {
    move_left_piston(dt);
    move_right_piston(dt);
}














//////////////////// PARTICLE FUNCTIONS ///////////////////////////

// Function to render a filled circle
void SDL_RenderFillCircle(SDL_Renderer* renderer, int centerX, int centerY, int radius_pixels) {
    for (int w = 0; w < radius_pixels * 2; w++) {
        for (int h = 0; h < radius_pixels * 2; h++) {
            int dx = radius_pixels - w;
            int dy = radius_pixels - h;
            if ((dx * dx + dy * dy) <= (radius_pixels * radius_pixels)) {
                SDL_RenderDrawPoint(renderer, centerX + dx, centerY + dy);
            }
        }
    }
}

// Function to render particles
void render_particles() {
    SDL_SetRenderDrawColor(_renderer, 255, 255, 255, 255);  // White color for particles

    for (int i = 0; i < NUM_PARTICLES; i++) {
        int pixel_x = (int)(X[i]);
        int pixel_y = (int)(Y[i]);
        int radius_pixels = (int)(PARTICLE_RADIUS);// * PIXELS_PER_SIGMA);  // 👈 SCALE properly to screen

        SDL_RenderFillCircle(_renderer, pixel_x, pixel_y, radius_pixels);
    }
}


// Function to check and handle boundary collisions
// newtons secon law
void handle_boundary_collision(int i) {
    if (X[i] - PARTICLE_RADIUS < XW1) {
        X[i] = XW1 + PARTICLE_RADIUS;
        Vx[i] = -Vx[i];
    }
    if (X[i] + PARTICLE_RADIUS > XW2) {
        X[i] = XW2 - PARTICLE_RADIUS;
        Vx[i] = -Vx[i];
    }
    if (Y[i] - PARTICLE_RADIUS < YW1) {
        Y[i] = YW1 + PARTICLE_RADIUS;
        Vy[i] = -Vy[i];
    }
    if (Y[i] + PARTICLE_RADIUS > YW2) {
        Y[i] = YW2 - PARTICLE_RADIUS;
        Vy[i] = -Vy[i];
    }
}





void handle_piston_collisions(int i, float sub_dt) {
    const float eps = 1e-6f;

   if (X[i] - Radius[i] <= piston_left_x) {
    float vrel = Vx[i] - vx_piston_left;
    if (vrel < 0.f) {
        double u = Vx[i];                       // before
        double v = 2.f*vx_piston_left - Vx[i];  // after
        double dE = 0.5*PARTICLE_MASS*(v*v - u*u); // change in kinetic energy
        piston_work_left += dE;                 // energy into gas (can be -)
        piston_hits_left++;

        Vx[i] = (float)v;
    }
    X[i] = piston_left_x + Radius[i] + eps; // prevent sticking
    }

    // Right piston
    if (X[i] + Radius[i] >= piston_right_x) {
        float vrel = Vx[i] - vx_piston_right;
        if (vrel > 0.f) {
            double u = Vx[i];
            double v = 2.f*vx_piston_right - Vx[i];
            double dE = 0.5*PARTICLE_MASS*(v*v - u*u);
            piston_work_right += dE;
            piston_hits_right++;

            Vx[i] = (float)v;
        }
        X[i] = piston_right_x - Radius[i] - eps; // prevent sticking
    }
}



//////// WALL FUNCTIONS /////////

void draw_wall() {
    if (!wall_enabled) return;  // Don't draw if the wall is disabled
       if (isnan(wall_x) || wall_x < XW1 || wall_x > XW2) {
        printf("🚨 Invalid wall_x = %f\n", wall_x);
            exit(1);
        }
    SDL_SetRenderDrawColor(_renderer, 255, 255, 255, 255);  // White color
    SDL_Rect wall_rect = { (int)(wall_x - (WALL_THICKNESS/2)), 0, (int)WALL_THICKNESS, SIM_HEIGHT };  // Thin vertical wall


    SDL_RenderFillRect(_renderer, &wall_rect);

    for (int i = 0; i < extra_wall_count; ++i) {
        float wx = extra_wall_positions[i];
        SDL_Rect extra_rect = { (int)(wx - (WALL_THICKNESS/2)), 0, (int)WALL_THICKNESS, SIM_HEIGHT };
        SDL_RenderFillRect(_renderer, &extra_rect);
    }
}


// Function to update wall position and velocity
void update_wall(float dt, float L0) {
    //steps_elapsed++;  // Increment the step counter

    // Inside update_wall()
    if (wall_hold_enabled && !wall_is_released) {
        vx_wall = 0.0f;

        if (!wall_position_is_managed_externally) {
            bool constrain_to_L0 = (num_internal_walls <= 1) && !preset_custom_positions;
            if (constrain_to_L0) {
                wall_x = XW1 + (L0 * PIXELS_PER_SIGMA);  // Only reset for single-wall setups
            }
        }

        if (all_wall_positions) {
            all_wall_positions[primary_wall_index] = wall_x;
        }

        if (steps_elapsed >= wall_hold_steps) {
            wall_is_released = true;
            wall_release_time = steps_elapsed * dt;
            primary_wall_release_time = simulation_time;
            printf("🔔 Wall released at step %d (t = %.3f)\n", steps_elapsed, wall_release_time);
        }
        return;
    }

    // === Normal Wall Motion Phase ===
    wall_x_old = wall_x;  // Store previous wall position
    wall_x += vx_wall * dt;  // Move wall based on current velocity

    if (steps_elapsed % 100 == 0) {
    printf("📊 Step %d | wall_x = %.6f | vx_wall = %.8f\n", steps_elapsed, wall_x, vx_wall);
    printf("wall_enabled = %d, wall_released = %d, left = %d, right = %d\n",
        wall_enabled, wall_is_released, left_count, right_count);
    }

    // === Enforce Wall Boundaries ===
    if (wall_x < XW1 + 0.1 * L0) {
        printf("🛑 Wall hit left boundary, clamping and resetting vx_wall\n");
        wall_x = XW1 + 0.1 * L0;
        vx_wall = -vx_wall;                 // elastic bounce
    }
    if (wall_x > XW2 - 0.1 * L0) {
        printf("🛑 Wall hit right boundary, clamping and resetting vx_wall\n");
        wall_x = XW2 - 0.1 * L0;
        vx_wall = -vx_wall;                 // elastic bounce
    }

}

static void update_extra_walls(float dt, float L0) {
    if (extra_wall_count <= 0 || !extra_wall_positions) {
        sync_all_wall_positions();
        return;
    }

    float margin = 0.1f * L0 * PIXELS_PER_SIGMA;
    float left_limit_base = (float)XW1 + margin;
    float right_limit_base = (float)XW2 - margin;

    for (int j = 0; j < extra_wall_count; ++j) {
        if (extra_wall_hold_enabled && extra_wall_hold_enabled[j] && !extra_wall_is_released[j]) {
            extra_wall_velocity[j] = 0.0f;
            extra_wall_positions[j] = extra_wall_old_positions[j];
            if (steps_elapsed >= wall_hold_steps) {
                extra_wall_is_released[j] = true;
                extra_wall_hold_enabled[j] = false;
                if (extra_wall_release_time) extra_wall_release_time[j] = simulation_time;
            }
            continue;
        }

        extra_wall_old_positions[j] = extra_wall_positions[j];
        extra_wall_positions[j] += extra_wall_velocity[j] * dt;

        float left_limit = left_limit_base;
        float right_limit = right_limit_base;

        if (extra_wall_positions[j] < left_limit) {
            extra_wall_positions[j] = left_limit;
            extra_wall_velocity[j] = -extra_wall_velocity[j];
        }
        if (extra_wall_positions[j] > right_limit) {
            extra_wall_positions[j] = right_limit;
            extra_wall_velocity[j] = -extra_wall_velocity[j];
        }
    }

    sync_all_wall_positions();
}

static void update_wall_metrics(void) {
    if (num_internal_walls < 1) {
        leftmost_wall_current_delta = leftmost_wall_baseline - wall_x;
        if (primary_wall_release_time >= 0.0 && leftmost_wall_release_time < 0.0) {
            leftmost_wall_release_time = primary_wall_release_time;
        }
        if (leftmost_wall_release_time >= 0.0 && !leftmost_wall_peak_recorded) {
            float delta = leftmost_wall_baseline - wall_x;
            if (delta > leftmost_wall_max_delta) {
                leftmost_wall_max_delta = delta;
                leftmost_wall_time_of_max = simulation_time;
                leftmost_wall_velocity_at_max = vx_wall;
                if (leftmost_wall_time_of_max > leftmost_wall_release_time) {
                    leftmost_wall_average_speed = delta / (leftmost_wall_time_of_max - leftmost_wall_release_time);
                }
                leftmost_wall_detected = true;
            } else if (leftmost_wall_detected && delta < leftmost_wall_max_delta * 0.95f) {
                leftmost_wall_peak_recorded = true;
            }
        }
        return;
    }

    sync_all_wall_positions();

    float min_pos = wall_x;
    float min_vel = vx_wall;
    double release_time = primary_wall_release_time;

    if (extra_wall_count > 0 && extra_wall_positions) {
        for (int j = 0; j < extra_wall_count; ++j) {
            float pos = extra_wall_positions[j];
            if (pos < min_pos) {
                min_pos = pos;
                min_vel = extra_wall_velocity ? extra_wall_velocity[j] : 0.0f;
                if (extra_wall_release_time) release_time = extra_wall_release_time[j];
            }
        }
    }

    leftmost_wall_current_delta = leftmost_wall_baseline - min_pos;

    if (release_time >= 0.0 && leftmost_wall_release_time < 0.0) {
        leftmost_wall_release_time = release_time;
    }

    if (leftmost_wall_release_time >= 0.0 && !leftmost_wall_peak_recorded) {
        float delta = leftmost_wall_current_delta;
        if (delta > leftmost_wall_max_delta) {
            leftmost_wall_max_delta = delta;
            leftmost_wall_time_of_max = simulation_time;
            leftmost_wall_velocity_at_max = min_vel;
            if (leftmost_wall_time_of_max > leftmost_wall_release_time) {
                leftmost_wall_average_speed = delta / (leftmost_wall_time_of_max - leftmost_wall_release_time);
            }
            leftmost_wall_detected = true;
        } else if (leftmost_wall_detected && delta < leftmost_wall_max_delta * 0.95f) {
            leftmost_wall_peak_recorded = true;
        }
    }
}


// =====================================================================
// Continuous Collision Detection (CCD) for a moving vertical slab wall.
/*  *
 * Discrete collision detection (DCD) only checks positions at the ends of
 * fixed time steps (Δt). If a fast object (or a thin/moving barrier) passes
 * through another object *between* those checks, DCD can miss the contact and
 * you see tunneling or “teleporting.”
 *
 * CCD fixes this by solving, within each step, for the exact time-of-impact
 * (TOI) assuming linear motion during Δt. We:
 *   1) Predict motion over Δt.
 *   2) Solve analytically for the first t_hit ∈ [0, Δt] when the shapes first touch.
 *   3) Advance to t_hit, resolve the collision (e.g., elastic reflection in the
 *      colliding frame), then advance the *remaining* time (Δt - t_hit).
 *   4) Repeat if there could be a second hit in the same Δt (rare).
 *
 * Why we use it here:
 * - The center wall and pistons can move and/or be thin.
 * - Particles can be fast relative to Δt.
 * - A piston can push a particle “inside” the wall slab within one substep.
 *   CCD prevents pops and tunneling by resolving the hit exactly when it happens.
 */
// Handles BOTH faces (left & right) of the wall with thickness,
// computes precise time-of-impact within dt, reflects in wall frame,
// advances remaining time, repeats (up to 2 hits per substep).
// =====================================================================


// =====================================================================
// Robust CCD vs moving slab (thickness) using Minkowski sum (inflate by R).
// - Works even when particle and wall move in the SAME direction (grazing).
// - Uses relative motion in wall frame: x_p(t)=xs+vxs*t, slab L(t)=L0+vwx*t, R(t)=R0+vwx*t
// - Inflate slab by +R: L'(t)=L(t)-R, R'(t)=R(t)+R, test point xs against [L',R']
// - Finite-mass elastic exchange along x updates BOTH vxs and vx_wall
// - Handles start-inside deterministically and adds small eps push-outs
// =====================================================================








// =====================================================================
// CCD against a finite-thickness, horizontally moving slab (runtime thickness).
// Robust, energy-conserving elastic response with optional finite-mass wall.
// =====================================================================
/*
In ccd_wall_step(i, dt, wall_x0, wall_vx_live, prev_side, wall_impulse_x_accum, force_locked):
	•	i — particle index.
	•	dt — duration of this (sub)step.
	•	wall_x0 — wall center x position frozen at start of this substep.
	•	wall_vx_live — pointer to the current wall velocity, updated instantly after every impulse so later collisions "see" the new speed.
	•	prev_side — where particle was relative to the un-inflated wall last step: -1 left, 0 inside thickness, +1 right (used as a hint for “start-inside” push-out direction).
	•	wall_impulse_x_accum — pointer to an accumulator (scalar) where we add the negative of the particle’s impulse along x (so later you do vx_wall += J_total / WALL_MASS once per fixed step).
	•	force_locked — if true the wall behaves as an infinite-mass mirror regardless of global hold/release flags.
Local names you see:
	•	X[i], Y[i] — particle position.
	•	Vx[i], Vy[i] — particle velocity.
	•	Radius[i] — particle radius.
	•	WALL_THICKNESS (or runtime var; see §3 below) — slab thickness.
	•	wL0 = wall_x0 - 0.5*WALL_THICKNESS, wR0 = wall_x0 + 0.5*WALL_THICKNESS — left/right slab faces at substep start.
	•	Inflated faces (Minkowski sum): L0p = wL0 - R, R0p = wR0 + R so we can treat the particle as a point.
	•	vrel = vxs - vwx — relative velocity in wall frame; collisions are solved in this frame.
	•	inside — whether the point is inside [L0p, R0p] (within a tiny tolerance).
	•	locked_wall — if the wall is held (pre-release). We treat it as infinite mass, i.e., mirror reflection in wall frame.
	•	eps_pos — minimal push-out distance (prevents numerical re-penetration).
	•	eps_t — tiny time epsilon (guards against division and end-point mistakes).
	•	v_hys — relative-velocity hysteresis; if |vrel| < v_hys we skip collision math and just free-fly (avoids chatter).

*/

static void ccd_wall_step(int i,           // particle index
                          float dt,        // substep time
                          const float wall_x0,   // wall center x at substep start (frozen pose)
                          float *wall_vx_live,   // live wall velocity (updated on-the-fly)
                          const int   prev_side, // -1 left, 0 inside, +1 right (hint)
                          float *wall_impulse_x_accum, // accumulates -J_particle for diagnostics
                          bool force_locked)
{
    // Bail if wall is off: advance particle by ballistic motion and return
    if (!wall_enabled) {                   // if wall disabled
        X[i] += Vx[i] * dt;                //   x ← x + vx·dt
        Y[i] += Vy[i] * dt;                //   y ← y + vy·dt
        return;                            //   done
    }

    // Nothing to do if dt is non-positive
    if (dt <= 0.f) return;                 // guard

    // Determine whether the wall is 'locked' (infinite mass behavior) this step
    const bool locked_wall = force_locked || (wall_hold_enabled && !wall_is_released); // held before release or explicitly locked

    // Snapshot particle state at the beginning of the substep
    float xs  = X[i];                      // particle x
    float ys  = Y[i];                      // particle y
    float vxs = Vx[i];                     // particle vx
    float vys = Vy[i];                     // particle vy
    const float R = Radius[i];             // particle radius

    // Compute slab faces at t=0 of this substep (frozen pose)
    const float wL0 = wall_x0 - 0.5f * WALL_THICKNESS; // slab left face
    const float wR0 = wall_x0 + 0.5f * WALL_THICKNESS; // slab right face
    float vwx = wall_vx_live ? *wall_vx_live : vx_wall; // live wall speed (updates after impulses)

    // Inflate faces by particle radius -> treat particle as a point
    const float L0p = wL0 - R;                // inflated left face
    const float R0p = wR0 + R;                // inflated right face

    // Numerical tolerances (tight → less artificial dissipation)
    const float eps_pos = fmaxf(1e-8f, 1e-6f * fmaxf(R, 1.f)); // push-out distance
    const float eps_t   = 1e-12f;                               // small time epsilon
    const float v_hys   = 1e-10f;                               // grazing hysteresis

    // Remaining time to advance in this substep (we may split on collision)
    float t_remain = dt;                    // start with full dt

    // Allow a bounded number of collision/advance cycles in this substep
    for (int iter = 0; iter < 100 && t_remain > eps_t; ++iter) { // up to 100 events
        float vrel = vxs - vwx;             // relative vx (particle in wall frame)

        // Inside test: is point strictly between inflated faces (with tiny slack)?
        const bool inside = (xs > L0p - eps_pos) && (xs < R0p + eps_pos); // inside?

        // If relative motion is too small (grazing), just free-flight to end
        if (fabsf(vrel) < v_hys) {          // nearly zero relative speed
            if (inside) {                    // if numerically inside
                int face = (prev_side != 0)  // pick face using previous side if known
                         ? prev_side
                         : (fabsf(xs - L0p) <= fabsf(R0p - xs) ? -1 : +1);

                const float m = (float)PARTICLE_MASS;
                const float M = (float)WALL_MASS;
                const float u1 = vxs;
                const float u2 = vwx;

                float v1;
                if (locked_wall) {
                    v1 = 2.f * u2 - u1;
                } else {
                    v1 = ((m - M) * u1 + 2.f * M * u2) / (m + M);
                    float Jp = m * (v1 - u1);
                    if (wall_impulse_x_accum) {
                        *wall_impulse_x_accum -= Jp;
                    }
                    if (wall_vx_live) {
                        *wall_vx_live -= Jp / WALL_MASS;
                        vwx = *wall_vx_live;
                    }
                }
                vxs = v1;

                xs = (face == -1)
                   ? (L0p - eps_pos)
                   : (R0p + eps_pos);
                xs += vxs * (1e-12f);
            }
            xs += vxs * t_remain;            // advance x by remaining time
            ys += vys * t_remain;            // advance y by remaining time
            t_remain = 0.f;                  // consumed all time
            break;                           // exit loop
        }

        // If we start inside the slab, resolve deterministically & push out
        if (inside) {                        // if currently inside
            int face = (prev_side != 0)      // decide nearest face (use hint first)
                     ? prev_side
                     : (fabsf(xs - L0p) <= fabsf(R0p - xs) ? -1 : +1);

            const float m = (float)PARTICLE_MASS; // particle mass
            const float M = (float)WALL_MASS;     // wall mass (finite or large)
            const float u1 = vxs;                 // particle pre-collision vx
            const float u2 = vwx;                 // wall vx (frozen pose during substep)

            float v1;                             // particle post-collision vx
            if (locked_wall) {                    // infinite-mass mirror mode
                v1 = 2.f * u2 - u1;               // reflect in wall frame
            } else {                              // finite-mass elastic 1D along x
                v1 = ((m - M)*u1 + 2.f*M*u2) / (m + M);
                float Jp = m * (v1 - u1);         // particle impulse (signed)
                if (wall_impulse_x_accum) {       // accumulate opposite impulse to wall
                    *wall_impulse_x_accum -= Jp;  // wall gets -Jp
                }
                if (wall_vx_live) {               // update wall velocity instantly
                    *wall_vx_live -= Jp / WALL_MASS;
                    vwx = *wall_vx_live;          // keep local copy in sync for later hits
                }
            }
            vxs = v1;                             // commit new vx

            xs = (face == -1)                     // place just outside face
               ? (L0p - eps_pos)
               : (R0p + eps_pos);
            xs += vxs * (1e-12f);                 // tiny drift to avoid re-penetration
            continue;                             // still have time → iterate again
        }

        // Compute earliest time-of-impact with either face within remaining time
        float thit = INFINITY;                    // initialize TOI
        int   face = 0;                           // -1=left, +1=right, 0=no hit

        if (vrel > 0.f) {                         // moving right in wall frame
            float tL = (L0p - xs) / vrel;         // time to left face
            if (tL >= -eps_t && tL <= t_remain + eps_t) { thit = tL; face = -1; }
        }
        if (vrel < 0.f) {                         // moving left in wall frame
            float tR = (R0p - xs) / vrel;         // time to right face
            if (tR >= -eps_t && tR <= t_remain + eps_t && tR < thit) { thit = tR; face = +1; }
        }

        // If no valid hit within remaining time, free-flight and finish
        if (face == 0 || !isfinite(thit) || thit > t_remain - eps_t) {
            xs += vxs * t_remain;                 // advance to end
            ys += vys * t_remain;
            t_remain = 0.f;                       // done
            break;
        }

        // Advance exactly to the collision time
        xs += vxs * thit;                         // x to impact
        ys += vys * thit;                         // y to impact

        // Elastic collision response along x at impact
        {
            const float m = (float)PARTICLE_MASS; // masses again
            const float M = (float)WALL_MASS;
            const float u1 = vxs;                 // pre-collision particle vx
            const float u2 = vwx;                 // wall vx at impact

            float v1;                             // post-collision particle vx
            if (locked_wall) {                    // infinite-mass wall behavior
                v1 = 2.f * u2 - u1;               // reflect
            } else {                              // finite-mass elastic formula
                v1 = ((m - M)*u1 + 2.f*M*u2) / (m + M);
                float Jp = m * (v1 - u1);         // particle impulse
                if (wall_impulse_x_accum) {       // accumulate wall impulse for diagnostics
                    *wall_impulse_x_accum -= Jp;  // wall accumulates opposite
                }
                if (wall_vx_live) {               // integrate wall velocity immediately
                    *wall_vx_live -= Jp / WALL_MASS;
                    vwx = *wall_vx_live;
                }
            }
            vxs = v1;                             // commit new vx
        }

        // Snap just outside the hit face to avoid re-penetration and add tiny drift
        xs = (face == -1) ? (L0p - eps_pos) : (R0p + eps_pos); // outside
        xs += vxs * (1e-12f);                   // micro drift along motion

        // Consume the time used up to the hit and continue with remainder
        t_remain -= thit;                        // reduce remaining time
    }

    // Write back the updated particle state at end of substep
    X[i]  = xs;                                  // commit x
    Y[i]  = ys;                                  // commit y
    Vx[i] = vxs;                                 // commit vx
    Vy[i] = vys;                                 // commit vy
}





/*
 * update_particles_with_substepping
 * ----------------------------------
 * Advances all particles through a fixed time step using optional sub-stepping.
 *  - Each sub-step runs CCD against the moving divider wall, piston collisions,
 *    and boundary reflections.
 *  - A uniform grid broad phase limits pairwise collision checks to neighbouring
 *    cells, keeping the cost reasonable for thousands of disks.
 *  - Any unusually deep overlaps (beyond a small tolerance) are recorded as
 *    potential "missed" collisions so we can monitor numerical stability.
 */
void update_particles_with_substepping(float dt, int* out_left, int* out_right) {
    if (dt <= 0.0f) return;

    float sub_dt = dt;

    #if SUBSTEPPING
        sub_dt = dt / SUBSTEPS;
    #endif

    int Side_old[NUM_PARTICLES];
    left_count = right_count = 0;

    const float missed_threshold = PARTICLE_RADIUS * 0.05f; // ~5% radius tolerance

    // --- Broad-phase acceleration structure (uniform grid) ---
    const float cell_size = fmaxf(4.0f * PARTICLE_RADIUS, 1.0f); // wide enough so neighbours capture overlaps
    const float domain_width  = (float)(XW2 - XW1);
    const float domain_height = (float)(YW2 - YW1);
    int grid_cols = (int)ceilf(domain_width  / cell_size);
    int grid_rows = (int)ceilf(domain_height / cell_size);
    if (grid_cols < 1) grid_cols = 1;
    if (grid_rows < 1) grid_rows = 1;
    int max_cells = grid_cols * grid_rows;

    int *cell_head = (int*)malloc((size_t)max_cells * sizeof(int));
    int *cell_next = (int*)malloc((size_t)NUM_PARTICLES * sizeof(int));
    int *particle_cell = (int*)malloc((size_t)NUM_PARTICLES * sizeof(int));
    if (!cell_head || !cell_next || !particle_cell) {
        printf("❌ Failed to allocate collision grid.\n");
        exit(1);
    }

    if (extra_wall_count > 0) {
        for (int j = 0; j < extra_wall_count; ++j) {
            extra_wall_old_positions[j] = extra_wall_positions[j];
        }
    }

    for (int step = 0; step < (int)SUBSTEPS; ++step) {
        for (int i = 0; i < NUM_PARTICLES; i++) {
            X_old[i] = X[i];

            if (wall_enabled) {
                if (X_old[i] < wall_x - WALL_THICKNESS / 2)      Side_old[i] = -1;
                else if (X_old[i] > wall_x + WALL_THICKNESS / 2) Side_old[i] =  1;
                else                                             Side_old[i] =  0;
            } else {
                Side_old[i] = 0;
            }

            const float w_x0_sub = wall_x_old;   // frozen pose for this fixed step
            const int prev_side = Side_old[i];   // -1 left, 0 inside, +1 right

            // CCD handles ballistic advance + wall interaction (now updating wall velocity instantly).
            ccd_wall_step(i, sub_dt, w_x0_sub, &vx_wall, prev_side, &wall_impulse_x_accum, false);

            if (wall_enabled && extra_wall_count > 0) {
                for (int w = 0; w < extra_wall_count; ++w) {
                    float wx = extra_wall_old_positions[w];
                    float wall_left = wx - WALL_THICKNESS / 2.0f;
                    float wall_right = wx + WALL_THICKNESS / 2.0f;
                    int extra_prev_side = 0;
                    if (X_old[i] < wall_left) extra_prev_side = -1;
                    else if (X_old[i] > wall_right) extra_prev_side = 1;
                    ccd_wall_step(i, sub_dt, extra_wall_old_positions[w], &extra_wall_velocity[w], extra_prev_side, NULL, true);
                }
            }

            handle_piston_collisions(i, sub_dt);
            handle_boundary_collision(i);
        }

        // --- rebuild uniform grid ---
        for (int c = 0; c < max_cells; ++c) cell_head[c] = -1;

        for (int i = 0; i < NUM_PARTICLES; ++i) {
            float x_local = (float)(X[i] - XW1);
            float y_local = (float)(Y[i] - YW1);
            int cx = (int)floorf(x_local / cell_size);
            int cy = (int)floorf(y_local / cell_size);
            if (cx < 0) cx = 0; else if (cx >= grid_cols) cx = grid_cols - 1;
            if (cy < 0) cy = 0; else if (cy >= grid_rows) cy = grid_rows - 1;
            int cell_index = cy * grid_cols + cx;
            particle_cell[i] = cell_index;
            cell_next[i] = cell_head[cell_index];
            cell_head[cell_index] = i;
        }

        // --- narrow-phase: only visit local cell neighbourhoods ---
        for (int i = 0; i < NUM_PARTICLES; ++i) {
            int cell_index = particle_cell[i];
            int cx = cell_index % grid_cols;
            int cy = cell_index / grid_cols;

            for (int dy = -1; dy <= 1; ++dy) {
                int ny = cy + dy;
                if (ny < 0 || ny >= grid_rows) continue;
                for (int dx = -1; dx <= 1; ++dx) {
                    int nx = cx + dx;
                    if (nx < 0 || nx >= grid_cols) continue;
                    int neighbor_idx = ny * grid_cols + nx;
                    for (int j = cell_head[neighbor_idx]; j != -1; j = cell_next[j]) {
                        if (j <= i) continue; // ensure each pair only once

                        float dxp = (float)(X[i] - X[j]);
                        float dyp = (float)(Y[i] - Y[j]);
                        float dist2 = dxp * dxp + dyp * dyp;
                        float min_d = 2.f * PARTICLE_RADIUS;

                        if (dist2 < min_d * min_d && dist2 > 0.f) {
                            float dist = sqrtf(dist2);
                            float nxp = dxp / dist;
                            float nyp = dyp / dist;

                            float relx = (float)(Vx[i] - Vx[j]);
                            float rely = (float)(Vy[i] - Vy[j]);
                            float vn = relx * nxp + rely * nyp;
                            if (vn < 0.f) {
                                Vx[i] -= vn * nxp;
                                Vy[i] -= vn * nyp;
                                Vx[j] += vn * nxp;
                                Vy[j] += vn * nyp;

                                float overlap = min_d - dist;
                                if (overlap > missed_threshold) {
                                    missed_collision_events++;
                                    if (overlap > worst_penetration_observed) {
                                        worst_penetration_observed = overlap;
                                    }
                                }
                                float push = 0.5f * overlap;
                                X[i] += nxp * push;  Y[i] += nyp * push;
                                X[j] -= nxp * push;  Y[j] -= nyp * push;
                            }
                        }
                    }
                }
            }
        }
    }

    free(cell_head);
    free(cell_next);
    free(particle_cell);

    if (segment_counts && segment_count > 0) {
        memset(segment_counts, 0, (size_t)segment_count * sizeof(int));
        for (int i = 0; i < NUM_PARTICLES; ++i) {
            int seg = segment_index_for_position(X[i]);
            if (seg < 0) seg = 0;
            if (seg >= segment_count) seg = segment_count - 1;
            segment_counts[seg]++;
        }
        left_count = segment_counts[0];
        right_count = segment_counts[segment_count - 1];
        if (segment_ke && segment_temperature) {
            memset(segment_ke, 0, (size_t)segment_count * sizeof(double));
            memset(segment_temperature, 0, (size_t)segment_count * sizeof(double));
            for (int i = 0; i < NUM_PARTICLES; ++i) {
                int seg = segment_index_for_position(X[i]);
                if (seg < 0) seg = 0;
                if (seg >= segment_count) seg = segment_count - 1;
                double ke = 0.5 * PARTICLE_MASS * (Vx[i] * Vx[i] + Vy[i] * Vy[i]);
                segment_ke[seg] += ke;
            }
            for (int seg = 0; seg < segment_count; ++seg) {
                if (segment_counts[seg] > 0) {
                    segment_temperature[seg] = segment_ke[seg] / (segment_counts[seg] * K_B);
                } else {
                    segment_temperature[seg] = 0.0;
                }
            }
        }
    } else {
        if (wall_enabled) {
            left_count = right_count = 0;
            for (int i = 0; i < NUM_PARTICLES; i++) {
                if (X[i] < wall_x - WALL_THICKNESS / 2) left_count++;
                else if (X[i] > wall_x + WALL_THICKNESS / 2) right_count++;
            }
        }
        if (segment_ke && segment_temperature) {
            memset(segment_ke, 0, (size_t)segment_count * sizeof(double));
            memset(segment_temperature, 0, (size_t)segment_count * sizeof(double));
        }
    }

    if (out_left)  *out_left  = left_count;
    if (out_right) *out_right = right_count;
}



////////////////////////////   ENERGY CALCULATION   ///////////////////////



// Expected kinetic energy for ideal gas: (3/2) N k_B T

// --> only applies if no interactions of particles
double kinetic_energy_expected() {
    return (3.0 / 2.0) * NUM_PARTICLES * K_B * TEMPERATURE;  // In Joules
}

double average_kinetic_energy_per_particle() {
    return kinetic_energy() / NUM_PARTICLES;
}

double kinetic_energy_total_system() {
    double ke_particles = kinetic_energy(); // existing function (particles only)
    double ke_wall = 0.5 * WALL_MASS * vx_wall * vx_wall;
    return ke_particles + ke_wall;
}




// /////// TMEPERATURE with dimension d ///////
//T = \frac{2}{d N k_B} \sum_{i=1}^{N} \frac{1}{2} m v_i^2
double temperature() {
    return (2.0 / 2.0) * kinetic_energy() / (NUM_PARTICLES * K_B);  // In Kelvin
}





///////////////////////////   ENTROPY CALCULATION   ///////////////////////
/// VELOCITY
////////  single velocity histogram update ///////

void compute_velocity_histogram() {
    // Reset histogram
    for (int i = 0; i < NUM_BINS; i++) velocity_histogram[i] = 0;

    // Compute speeds and bin them
    for (int i = 0; i < NUM_PARTICLES; i++) {
        // v = \sqrt{V_x^2 + V_y^2}

        float speed = sqrt(Vx[i] * Vx[i] + Vy[i] * Vy[i]);
        int bin = (int)(speed / MAX_VELOCITY * NUM_BINS);
        if (bin >= 0 && bin < NUM_BINS) velocity_histogram[bin]++;
    }
}





///////////////////////////   ENTROPY CALCULATION   ///////////////////////

///////////// THermodynamic Entropy S = -k_B ln (THETA) //////////////
// Function to compute kinetic entropy
double entropy_kinetic() {
    double S_kinetic = 0.0;
    
    for (int i = 0; i < NUM_PARTICLES; i++) {
        double speed_square = Vx[i] * Vx[i] + Vy[i] * Vy[i];
        //double speed = sqrt(speed_square);
        
        // Using the Maxwell-Boltzmann distribution entropy formula
        if (speed_square > 0) {
            double f_v = exp(-speed_square / (2 * K_B * TEMPERATURE / PARTICLE_MASS));
            S_kinetic -= f_v * log(f_v);
        }
    }
    
    return S_kinetic * K_B;
}
// Function to compute position entropy
double entropy_position() {
    // The number of accessible positions is proportional to the area of the box
    return K_B * log((XW2 - XW1) * (YW2 - YW1));
}


///////////// SHANNON with histograms  Entropy S = sum pi * ln(pi)  //////////////

double entropy() {
    double S = 0.0;
    int total_counts = 0;
    
    for (int i = 0; i < NUM_BINS; i++) {
        total_counts += velocity_histogram[i];
    }
    
    if (total_counts == 0) return 0.0;  // Avoid division by zero
    
    for (int i = 0; i < NUM_BINS; i++) {
        if (velocity_histogram[i] > 0) {
            double p_i = (double)velocity_histogram[i] / total_counts;
            S -= p_i * log(p_i);
        }
    }
    return S * K_B;  // Entropy in J/K
}



/// POSITION
void compute_position_histogram() {
    // Reset histogram
    for (int i = 0; i < NUM_BINS; i++) position_histogram[i] = 0;

    // Compute position bins
    for (int i = 0; i < NUM_PARTICLES; i++) {
        int bin_x = (int)(X[i] / MAX_X * NUM_BINS);  // X-coordinate bin
        int bin_y = (int)(Y[i] / MAX_Y * NUM_BINS);  // Y-coordinate bin
        int bin = bin_x + bin_y * NUM_BINS;  // Convert to a 1D index

        if (bin >= 0 && bin < NUM_BINS * NUM_BINS) position_histogram[bin]++;
    }
}

double position_entropy() {
    double S = 0.0;
    int total_counts = 0;
    
    for (int i = 0; i < NUM_BINS * NUM_BINS; i++) {
        total_counts += position_histogram[i];
    }
    
    if (total_counts == 0) return 0.0;  // Avoid division by zero
    
    for (int i = 0; i < NUM_BINS * NUM_BINS; i++) {
        if (position_histogram[i] > 0) {
            double p_i = (double)position_histogram[i] / total_counts;
            S -= p_i * log(p_i);
        }
    }
    return S * K_B;  // Entropy in J/K
}






///////////////////////////   FREE ENERGY CALCULATION   ///////////////////////
double total_free_energy(double T, double S) {
    double ke = kinetic_energy();  // Internal energy (J)
    double F = ke - T * S;  // Free energy (J)
    return F;
}




//////LOGGING ////////
// Function to log energy and temperature to a file


void log_energy(FILE *fptr) {
    double ke_particles = kinetic_energy();
    double ke_wall = 0.5 * WALL_MASS * vx_wall * vx_wall;
    double ke_total = ke_particles + ke_wall;
    double T = temperature(); // still based on particle KE (as you have)


    // Compute kinetic and position entropy
    double S_kinetic = entropy_kinetic();
    double S_position = entropy_position();
    double S_total = S_kinetic + S_position;
    
    /*
    printf("Kinetic Entropy: %e J/K\n", S_kinetic);
    printf("Position Entropy: %e J/K\n", S_position);
    printf("Total Entropy: %e J/K\n", S_total);
    */


    compute_velocity_histogram();
    float S_v = entropy();

    compute_position_histogram();
    float S_p = position_entropy();


    float F = total_free_energy(T, S_total);

    fprintf(fptr,
      "KE_particles: %.8e, KE_wall: %.8e, KE_total: %.8e, T: %.8e, S_v: %.8e, F: %.8e\n",
      ke_particles, ke_wall, ke_total, T, S_v, F);    // Print the calculated values to the console
    //printf("Kinetic Energy: %.6e J, Kinetic Energy expected: %.6e , AVG Kinetic Energy per particle: %.6e, Temperature: %.6e K, Entropy velocity: %.6e J/K, Entropy position: %.6e J/K, Total Free Energy: %.6e J\n", ke,ke_expected, avg_ke_pp, T, S_v, S_p, F);
}








////////// KEYBOARD INPUT HANDLING //////////
// Keyboard input handling
void keysSimulation(SDL_Event e) {
    if (e.type == SDL_KEYDOWN) {
        switch (e.key.keysym.sym) {

            case SDLK_SPACE:
                simulation_started = 1;
                break;
            
            case SDLK_q:
                exit(0);
                break;

            case SDLK_w:  // Toggle wall with 'w' or 'W'
                wall_enabled = !wall_enabled;
                printf("Wall %s\n", wall_enabled ? "Enabled" : "Disabled");
                break;

            case SDLK_r:  // 🔥 Release wall manually
            if (wall_hold_enabled && !wall_is_released) {
                wall_is_released = true;
                wall_hold_enabled = false;
                primary_wall_release_time = simulation_time;
                if (extra_wall_count > 0) {
                    for (int j = 0; j < extra_wall_count; ++j) {
                        if (extra_wall_hold_enabled) extra_wall_hold_enabled[j] = false;
                        if (extra_wall_is_released) extra_wall_is_released[j] = true;
                        if (extra_wall_release_time) extra_wall_release_time[j] = simulation_time;
                    }
                }
                leftmost_wall_release_time = simulation_time;
                printf("🔔 Wall manually released by key 'r'.\n");}
                break;

            case SDLK_p:  // Pause/unpause simulation
            paused = !paused;
            printf("Simulation %s\n", paused ? "Paused" : "Running");
            break;

            case SDLK_ESCAPE:
                printf("🔄 Resetting simulation.\n");
                reset_simulation_state();
                break;
        }
    }
}


/// Constants for acceleration INIT piston
const float acceleration = 25;  // Increased acceleration
const float max_velocity = 1000000.0;   // Increased max speed
const float deceleration = 25;


// Function to handle piston movement based on user input
void keysPiston(SDL_Event e) {
    if (e.type == SDL_KEYDOWN) {
        switch (e.key.keysym.sym) {

        // Left piston movement
        case SDLK_a:
            if (piston_left_x > XW1) {
                vx_piston_left -= acceleration;
                if (vx_piston_left < -max_velocity) vx_piston_left = -max_velocity;
            }
            break;
        case SDLK_d:
            if (piston_left_x < piston_right_x - 20) {
                vx_piston_left += acceleration;
                if (vx_piston_left > max_velocity) vx_piston_left = max_velocity;
            }
            break;
        case SDLK_s: // Stop left piston immediately
            vx_piston_left = 0;
            break;

        // Right piston movement
        case SDLK_LEFT:
            if (piston_right_x > piston_left_x + 20) {
                vx_piston_right -= acceleration;
                if (vx_piston_right < -max_velocity) vx_piston_right = -max_velocity;
            }
            piston_step_active = false;
            break;
        case SDLK_RIGHT:
            if (piston_right_x < XW2) {
                vx_piston_right += acceleration;
                if (vx_piston_right > max_velocity) vx_piston_right = max_velocity;
            }
            piston_step_active = false;
            break;
        case SDLK_UP: // Stop right piston immediately
            vx_piston_right = 0;
            piston_step_active = false;
            break;
        case SDLK_t: { // trigger piston step
            float travel = (float)SIM_WIDTH * 0.20f;
            float target = (float)XW2 - travel;
            float min_target = piston_left_x + 25.0f;
            if (target < min_target) target = min_target;
            piston_step_target = target;
            piston_step_speed = -fabsf((float)SIM_WIDTH) * 0.5f;
            vx_piston_right = piston_step_speed;
            piston_step_active = true;
            simulation_started = 1;
            printf("➡️  Piston step triggered: target %.2f px, speed %.2f\n", piston_step_target, piston_step_speed);
        }
            break;
        }
    }
    else if (e.type == SDL_KEYUP) {
        switch (e.key.keysym.sym) {
        case SDLK_a:
        case SDLK_d:
            if (vx_piston_left > 0) vx_piston_left -= deceleration;
            if (vx_piston_left < 0) vx_piston_left += deceleration;
            if (fabs(vx_piston_left) < deceleration) vx_piston_left = 0;
            break;

        case SDLK_LEFT:
        case SDLK_RIGHT:
            if (vx_piston_right > 0) vx_piston_right -= deceleration;
            if (vx_piston_right < 0) vx_piston_right += deceleration;
            if (fabs(vx_piston_right) < deceleration) vx_piston_right = 0;
            break;
        }
    }
}










/////////  HISTO AND SIMULATION RENDERING /////////
// Function to handle particle- piston collisions

void draw_maxwell_boltzmann_curve(float temperature, int max_count) {
    float sigma = sqrt(K_B * temperature / PARTICLE_MASS);
    float norm_factor = NUM_PARTICLES * (SIM_WIDTH / NUM_BINS);

    for (int i = 0; i < NUM_BINS; i++) {
        float v = MAX_VELOCITY * i / NUM_BINS;
        float f_v = (v * v) * expf(-PARTICLE_MASS * v * v / (2.0f * K_B * temperature));
        float f_scaled = f_v * norm_factor;

        int height = (int)((f_scaled * HIST_HEIGHT) / max_count);

        SDL_SetRenderDrawColor(_renderer, 255, 255, 0, 255); // Yellow curve
        SDL_RenderDrawPoint(_renderer, XW1 + i * (SIM_WIDTH / NUM_BINS), YW1 + SIM_HEIGHT + (HIST_HEIGHT - height));
    }
}


void draw_thick_line_vertical(SDL_Renderer *renderer, int x, int y1, int y2, int thickness) {
    if (y1 > y2) {
        int temp = y1;
        y1 = y2;
        y2 = temp;
    }
    SDL_Rect rect = {x - thickness / 2, y1, thickness, y2 - y1 + 1};
    SDL_RenderFillRect(renderer, &rect);
}

void render_velocity_histograms() {
    SDL_SetRenderDrawBlendMode(_renderer, SDL_BLENDMODE_BLEND);

    float measured_T = compute_measured_temperature_from_ke();
    float sigma = sqrtf(K_B * measured_T / PARTICLE_MASS);

    float center_x = XW1 + L0_UNITS * PIXELS_PER_SIGMA;
    int bin_width = SIM_WIDTH / NUM_BINS;

    // === Histogram bins ===
    int binsX[NUM_BINS] = {0};
    int binsY[NUM_BINS] = {0};
    int binsSpeed[NUM_BINS] = {0};

    int max_binX = 1, max_binY = 1, max_binSpeed = 1;

    for (int i = 0; i < NUM_PARTICLES; i++) {
        // vx, vy, v
        int binX = (int)((Vx[i] + MAX_VELOCITY) / (2 * MAX_VELOCITY) * NUM_BINS);
        int binY = (int)((Vy[i] + MAX_VELOCITY) / (2 * MAX_VELOCITY) * NUM_BINS);
        float speed = sqrtf(Vx[i] * Vx[i] + Vy[i] * Vy[i]);
        int binV = (int)(speed / MAX_VELOCITY * NUM_BINS);

        if (binX >= 0 && binX < NUM_BINS) binsX[binX]++;
        if (binY >= 0 && binY < NUM_BINS) binsY[binY]++;
        if (binV >= 0 && binV < NUM_BINS) binsSpeed[binV]++;
    }

    for (int i = 0; i < NUM_BINS; i++) {
        if (binsX[i] > max_binX) max_binX = binsX[i];
        if (binsY[i] > max_binY) max_binY = binsY[i];
        if (binsSpeed[i] > max_binSpeed) max_binSpeed = binsSpeed[i];
    }

    // === Vx Histogram (white, bottom, centered at XW1 + L0) ===
    for (int i = 0; i < NUM_BINS; i++) {
        int height = (binsX[i] * HIST_HEIGHT) / max_binX;
        SDL_SetRenderDrawColor(_renderer, 255, 255, 255, 150);
        SDL_Rect rect = {
            (int)(center_x + (i - NUM_BINS / 2) * bin_width),
            YW1 + SIM_HEIGHT + (HIST_HEIGHT - height),
            bin_width,
            height
        };
        SDL_RenderFillRect(_renderer, &rect);
    }

    // === Speed Histogram (yellow, bottom, same center) ===
    for (int i = 0; i < NUM_BINS; i++) {
        int height = (binsSpeed[i] * HIST_HEIGHT) / max_binSpeed;
        SDL_SetRenderDrawColor(_renderer, 255, 255, 0, 255);
        SDL_Rect rect = {
            (int)(center_x + (i - NUM_BINS / 2) * bin_width),
            YW1 + SIM_HEIGHT + (HIST_HEIGHT - height),
            bin_width,
            height
        }; 
        SDL_RenderFillRect(_renderer, &rect);
    }

    // === Vy Histogram (light blue, right) ===
    for (int i = 0; i < NUM_BINS; i++) {
        int width = (binsY[i] * HIST_WIDTH) / max_binY;
        SDL_SetRenderDrawColor(_renderer, 100, 200, 255, 180);
        SDL_Rect rect = {
            XW1 + SIM_WIDTH - HIST_WIDTH,
            YW1 + i * (SIM_HEIGHT / NUM_BINS),
            width,
            SIM_HEIGHT / NUM_BINS
        };
        SDL_RenderFillRect(_renderer, &rect);
    }

    // === Overlay Maxwell-Boltzmann curve for v (yellow: live, red: init) ===
    float normalization = NUM_PARTICLES * bin_width * 50.0f;  // 250 is a scaling factor for visibility



    for (int i = 1; i < NUM_BINS; i++) {
        float v1 = ((float)(i - 1) / NUM_BINS) * MAX_VELOCITY;
        float v2 = ((float)i / NUM_BINS) * MAX_VELOCITY;

        float f1 = (PARTICLE_MASS * v1 / (K_B * measured_T)) * expf(-PARTICLE_MASS * v1 * v1 / (2 * K_B * measured_T));
        float f2 = (PARTICLE_MASS * v2 / (K_B * measured_T)) * expf(-PARTICLE_MASS * v2 * v2 / (2 * K_B * measured_T));

        int x1 = center_x + (i - 1 - NUM_BINS / 2) * bin_width;
        int x2 = center_x + (i - NUM_BINS / 2) * bin_width;

        int y1 = YW1 + SIM_HEIGHT + HIST_HEIGHT - (int)(f1 * normalization);
        int y2 = YW1 + SIM_HEIGHT + HIST_HEIGHT - (int)(f2 * normalization);

        SDL_SetRenderDrawColor(_renderer, 255, 255, 0, 255);
        draw_thick_line_vertical(_renderer, x1, y1, y2, 2); // e.g., thickness = 2
    }

    // --- Initial MB curve (red)
    float sigma0 = sqrtf(K_B * TEMPERATURE / PARTICLE_MASS);
    for (int i = 1; i < NUM_BINS; i++) {
        float v1 = ((float)(i - 1) / NUM_BINS) * MAX_VELOCITY;
        float v2 = ((float)i / NUM_BINS) * MAX_VELOCITY;

        float f1 = (PARTICLE_MASS * v1 / (K_B * TEMPERATURE)) * expf(-PARTICLE_MASS * v1 * v1 / (2 * K_B * TEMPERATURE));
        float f2 = (PARTICLE_MASS * v2 / (K_B * TEMPERATURE)) * expf(-PARTICLE_MASS * v2 * v2 / (2 * K_B * TEMPERATURE));

        int x1 = center_x + (i - 1 - NUM_BINS / 2) * bin_width;
        int x2 = center_x + (i - NUM_BINS / 2) * bin_width;

        int y1 = YW1 + SIM_HEIGHT + HIST_HEIGHT - (int)(f1 * normalization);
        int y2 = YW1 + SIM_HEIGHT + HIST_HEIGHT - (int)(f2 * normalization);

        SDL_SetRenderDrawColor(_renderer, 255, 0, 0, 120);
        draw_thick_line_vertical(_renderer, x1, y1, y2, 2); // e.g., thickness = 2
    }
}










void draw_clear_screen() {
    SDL_SetRenderDrawColor(_renderer, 0, 0, 0, 255); // Black color
    SDL_RenderClear(_renderer);
}

void draw_simulation_boundary() {
    // Draw a yellow border around the simulation area
    SDL_SetRenderDrawColor(_renderer, 255, 255, 0, 255); // Yellow color
    SDL_Rect sim_border = {XW1, YW1, SIM_WIDTH, SIM_HEIGHT};
    SDL_RenderDrawRect(_renderer, &sim_border);
}


void draw_coordinate_system(SDL_Renderer *renderer) {
    // Set color to red
    SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255);

    // Define start point (origin)
    int origin_x = 20;
    int origin_y = 20;
    
    // Define arrow length
    int axis_length = 50;

    // Draw X-axis (right)
    SDL_RenderDrawLine(renderer, origin_x, origin_y, origin_x + axis_length, origin_y);
    SDL_RenderDrawLine(renderer, origin_x + axis_length, origin_y, origin_x + axis_length - 5, origin_y - 5); // Arrow tip
    SDL_RenderDrawLine(renderer, origin_x + axis_length, origin_y, origin_x + axis_length - 5, origin_y + 5);

    // Draw Y-axis (down)
    SDL_RenderDrawLine(renderer, origin_x, origin_y, origin_x, origin_y + axis_length);
    SDL_RenderDrawLine(renderer, origin_x, origin_y + axis_length, origin_x - 5, origin_y + axis_length - 5); // Arrow tip
    SDL_RenderDrawLine(renderer, origin_x, origin_y + axis_length, origin_x + 5, origin_y + axis_length - 5);

    // (Optional) Add text labels: "X" and "Y" using SDL_ttf or simple rectangle markers
}



// draw text
void draw_text(SDL_Renderer *renderer, TTF_Font *font, const char *text, int x, int y, SDL_Color color) {
    SDL_Surface *text_surface = TTF_RenderText_Blended(font, text, color);
    SDL_Texture *text_texture = SDL_CreateTextureFromSurface(renderer, text_surface);
    SDL_Rect dst = { x, y, text_surface->w, text_surface->h };
    SDL_RenderCopy(renderer, text_texture, NULL, &dst);
    SDL_FreeSurface(text_surface);
    SDL_DestroyTexture(text_texture);
}


// ---------- HUD: energy & temperature ----------
static inline void render_energy_hud(SDL_Renderer *r, TTF_Font *f) {
    // compute energies
    const double KEp = kinetic_energy();
    const double KEw = 0.5 * WALL_MASS * vx_wall * vx_wall;
    const double KEt = KEp + KEw;

    // temperatures
    const double Tp  = KEp / (NUM_PARTICLES * K_B);          // 2D gas
    const double Tw  = KEw / (0.5 * K_B);                    // 1 DOF wall
    const double Tt  = KEt / ((NUM_PARTICLES + 0.5) * K_B);  // effective total

    // format
    char line0[1280], line1[128], line2[128], line3[128], line4[128], line5[128];
    if (num_internal_walls > 0) {
        double delta_max_sigma = leftmost_wall_max_delta / PIXELS_PER_SIGMA;
        snprintf(line0, sizeof(line0), "T_total: %.3f    maxΔx_left_wall: %.4f σ", Tt, delta_max_sigma);
    } else {
        snprintf(line0, sizeof(line0), "T_total: %.3f", Tt);
    }
    snprintf(line1, sizeof(line1), "KE_tot: %.6e", KEt);
    snprintf(line2, sizeof(line2), "KE_particles: %.6e", KEp);
    snprintf(line3, sizeof(line3), "KE_wall: %.6e", KEw);
    snprintf(line4, sizeof(line4), "T_particles: %.3f", Tp);
    snprintf(line5, sizeof(line5), "T_wall: %.3f", Tw);

    // draw starting where T_measured used to be
    SDL_Color white = {255,255,255,255};
    int x = XW1 + 1000;  // right side
    int y = YW1 + HEIGHT_UNITS * PIXELS_PER_SIGMA;

    draw_text(r, f, line0, x, y, white);  y += 20;
    draw_text(r, f, line1, x, y, white);  y += 20;
    draw_text(r, f, line2, x, y, white);  y += 20;
    draw_text(r, f, line3, x, y, white);  y += 20;
    draw_text(r, f, line4, x, y, white);  y += 20;
    draw_text(r, f, line5, x, y, white);
}

// FFT BARS
void render_fft_histogram_bottom(SDL_Renderer* renderer, float* fft_data, int fft_len) {
    // Normalize FFT data
    float max_val = 1.0f;
    for (int i = 0; i < fft_len; ++i) {
        if (fft_data[i] > max_val) max_val = fft_data[i];
    }

    // FFT bars drawn like Vx histogram, but different color (e.g. red)
    for (int i = 0; i < fft_len; ++i) {
        float normalized = fft_data[i] / max_val;
        int height = (int)(normalized * HIST_HEIGHT);

        SDL_SetRenderDrawColor(renderer, 255, 0, 0, 200); // Red with some transparency
        SDL_Rect rect = {
            XW1 + i * (SIM_WIDTH / fft_len),
            YW1 + SIM_HEIGHT + (HIST_HEIGHT - height),
            SIM_WIDTH / fft_len,
            height
        };
        SDL_RenderFillRect(renderer, &rect);
    }
}





//// CALCULATE PACKING FRACTION
void log_packing_fractions(FILE *log, const char *label) {
    float particle_area_unitless = M_PI * PARTICLE_RADIUS_UNIT * PARTICLE_RADIUS_UNIT;
    float particle_area_actual = M_PI * PARTICLE_RADIUS * PARTICLE_RADIUS;

    float box_area_unitless = 2.0f * L0_UNITS * HEIGHT_UNITS;

    float L0_scaled = L0_UNITS * PIXELS_PER_SIGMA;
    float height_scaled = HEIGHT_UNITS * PIXELS_PER_SIGMA;
    float box_area_scaled = 2.0f * L0_scaled * height_scaled;
    float wall_buffer_area = 2.0f * WALL_THICKNESS * height_scaled;

    float packing_unitless = NUM_PARTICLES * particle_area_unitless / box_area_unitless;
    float packing_actual_bound = NUM_PARTICLES * particle_area_actual / box_area_scaled;
    float packing_actual_excl_wall = NUM_PARTICLES * particle_area_actual / (box_area_scaled - wall_buffer_area);

    //fprintf(log, "# Packing_UNITLESS (%s) = %.6f\n", label, packing_unitless);
    //fprintf(log, "# Packing_ACTUAL_BOUND (%s) = %.6f\n", label, packing_actual_bound);
    //fprintf(log, "# Packing_ACTUAL_WITH_WALLBUFFER (%s) = %.6f\n", label, packing_actual_excl_wall);
    fprintf(log, "Time, Wall_X, Displacement(σ), Left_Count, Right_Count, Packing_UNITLESS, Packing_ACTUAL_BOUND, Packing_ACTUAL_WITH_WALLBUFFER\n");
}








/*
 * run_speed_of_sound_experiments
 * -------------------------------
 * Automates the divider-release study for multiple box lengths and wall masses.
 *  - Rebuilds output folders, sweeps through the configured parameter grid,
 *    and re-initialises the simulation for each combination.
 *  - After release, records wall displacement (and optional packing fractions)
 *    so we can post-process the oscillation frequency / effective sound speed.
 *  - Reuses the live missed-collision diagnostics to highlight any numerical
 *    issues that might appear during long unattended batches.
 */
void run_speed_of_sound_experiments() {
    system("mkdir -p experiments_speed_of_sound/mode0_real_units");
    system("mkdir -p experiments_speed_of_sound/mode1_normalized_units");
    system("rm -f experiments_speed_of_sound/mode0_real_units/*.csv");
    system("rm -f experiments_speed_of_sound/mode1_normalized_units/*.csv");

    wall_enabled = 1;                        // ensure wall exists for CCD
    wall_position_is_managed_externally = 1; // you control center via L0

    static const float default_lengths[] = {7.5f, 10.0f, 15.0f};
    static const int   default_wall_mass_factors[] = {20, 50, 100, 200, 500, 1000};

    const float *lengths = cli_experiment_lengths_count
        ? cli_experiment_lengths
        : default_lengths;
    size_t lengths_count = cli_experiment_lengths_count
        ? cli_experiment_lengths_count
        : sizeof(default_lengths) / sizeof(default_lengths[0]);

    const int *wall_mass_factors = cli_experiment_wall_masses_count
        ? cli_experiment_wall_masses
        : default_wall_mass_factors;
    size_t wall_mass_count = cli_experiment_wall_masses_count
        ? cli_experiment_wall_masses_count
        : sizeof(default_wall_mass_factors) / sizeof(default_wall_mass_factors[0]);

    int num_repeats = cli_experiment_repeats > 0 ? cli_experiment_repeats : 1;

    printf("Running each experiment for %d steps after wall release\n", num_steps);
    if (num_steps > 10000000) {
        printf("⚠️ WARNING: Simulation time is very long! Consider reducing the number of steps.\n");
    } else if (num_steps < 1000) {
        printf("⚠️ WARNING: Simulation time is very short, SET EMERGENCY STEPS!\n");
        num_steps = (int)(200000 / FIXED_DT);
    }

    for (size_t l = 0; l < lengths_count; ++l) {
        for (size_t m = 0; m < wall_mass_count; ++m) {
            for (int r = 0; r < num_repeats; ++r) {

                float L0 = lengths[l];
                int wall_mass_factor = wall_mass_factors[m];
                wall_mass_runtime = PARTICLE_MASS * wall_mass_factor;
                L0_UNITS = L0;
                initialize_simulation_dimensions();

                // center wall at L0 (in pixels)
                wall_x     = XW1 + (L0 * PIXELS_PER_SIGMA);
                wall_x_old = wall_x;   // ← important for the first CCD step
                vx_wall    = 0.0f;

                wall_hold_enabled = true;
                wall_is_released  = false;
                steps_elapsed     = 0;
                simulation_time   = 0.0;

                initialize_simulation();

                #if MOLECULAR_MODE
                    const char* mode_folder = "experiments_speed_of_sound/mode1_normalized_units/";
                #else
                    const char* mode_folder = "experiments_speed_of_sound/mode0_real_units/";
                #endif

                char filename[512];
                int L0_int = (int)(L0 * 10);
                sprintf(filename, "%swall_x_positions_L0_%d_wallmassfactor_%d_run%d.csv",
                        mode_folder, L0_int, wall_mass_factor, r);
                FILE *wall_log = fopen(filename, "w");
                if (!wall_log) { printf("❌ Could not open log file.\n"); continue; }

                if (log_packing_fraction) {
                    log_packing_fractions(wall_log, "main");
                } else {
                    fprintf(wall_log, "Time, Wall_X, Displacement(σ), Left_Count, Right_Count\n");
                }

                printf("🔬 Running: L0 = %.1f, M = %d*m, run = %d\n", L0, wall_mass_factor, r);
                printf("🔍 Initial wall_x = %.3f, vx_wall = %.6f\n", wall_x, vx_wall);

                int left_particles = 0, right_particles = 0;
                int recorded_steps = 0;
                missed_collision_events = 0;
                worst_penetration_observed = 0.0;

                while (recorded_steps < num_steps) {
                    // --- START OF STEP ---
                    wall_x_old = wall_x;   // ← CCD uses wall pose from start of step
                    steps_elapsed++;       // ← tick BEFORE update_wall so hold/release can trigger

                    // physics
                    wall_impulse_x_accum = 0.f; // (re)start accumulator for this step
                    update_particles_with_substepping(FIXED_DT, &left_particles, &right_particles);

                    // wall velocity already updated inside ccd_wall_step; accumulator is purely diagnostic now.
                    // (optional) clear accumulator here if you prefer step-scoped use
                    // wall_impulse_x_accum = 0.f;

                    update_pistons(FIXED_DT);
                    update_wall(FIXED_DT, L0_UNITS);
                    update_extra_walls(FIXED_DT, L0_UNITS);
                    sync_all_wall_positions();
                    update_wall_metrics();

                    // time bookkeeping
                    simulation_time += FIXED_DT;

                    // release check using true sim time
                    if (steps_elapsed == wall_hold_steps && !wall_is_released) {
                        wall_release_time = simulation_time;
                        wall_is_released  = true;
                        recorded_steps    = 0;
                        printf("🔔 Wall released at step %d (t = %.3f)\n", steps_elapsed, wall_release_time);
                    }

                    if (!wall_is_released) continue;

                    float time_after_release = (float)(simulation_time - wall_release_time);
                    float wall_x_sigma = wall_x / PIXELS_PER_SIGMA;
                    float disp = wall_x_sigma - (XW1 + XW2) / (2.0f * PIXELS_PER_SIGMA);

                    if (!log_packing_fraction) {
                        fprintf(wall_log, "%.3f, %.3f, %.3f, %d, %d\n",
                                time_after_release, wall_x_sigma, disp, left_particles, right_particles);
                    } else {
                        float particle_area_unitless = (float)(M_PI * PARTICLE_RADIUS_UNIT * PARTICLE_RADIUS_UNIT);
                        float particle_area_actual   = (float)(M_PI * PARTICLE_RADIUS * PARTICLE_RADIUS);
                        float box_area_unitless      = 2.0f * L0_UNITS * HEIGHT_UNITS;
                        float box_area_scaled        = 2.0f * (L0_UNITS * PIXELS_PER_SIGMA) * (HEIGHT_UNITS * PIXELS_PER_SIGMA);
                        float wall_area              = 0.5f * WALL_THICKNESS * (HEIGHT_UNITS * PIXELS_PER_SIGMA);
                        float wall_buffer_area       = 2.0f * wall_area;

                        float packing_unitless       = NUM_PARTICLES * particle_area_unitless / box_area_unitless;
                        float packing_actual_bound   = NUM_PARTICLES * particle_area_actual / box_area_scaled;
                        float packing_actual_excl    = NUM_PARTICLES * particle_area_actual / (box_area_scaled - wall_buffer_area);

                        fprintf(wall_log, "%.3f, %.3f, %.3f, %d, %d, %.6f, %.6f, %.6f\n",
                                time_after_release, wall_x_sigma, disp,
                                left_particles, right_particles,
                                packing_unitless, packing_actual_bound, packing_actual_excl);
                    }

                    recorded_steps++;

                    if ((recorded_steps % 1000) == 0 && missed_collision_events > 0) {
                        printf("[CCD][exp] step=%d  potential_misses=%ld  worst_penetration=%.6g\n",
                               recorded_steps, missed_collision_events, worst_penetration_observed);
                    }
                    if ((recorded_steps % 1000) == 0 && segment_counts && segment_count > 0) {
                        printf("[SEG][exp] step=%d", recorded_steps);
                        for (int seg = 0; seg < segment_count; ++seg) {
                            double ke = segment_ke ? segment_ke[seg] : 0.0;
                            double temp = segment_temperature ? segment_temperature[seg] : 0.0;
                            printf("  S%d=%d KE=%.3e T=%.3f", seg, segment_counts[seg], ke, temp);
                        }
                        printf("\n");
                    }
                    if (leftmost_wall_peak_recorded && !leftmost_wall_logged) {
                        printf("[WALL][exp] first-peak Δx=%.4f px at t=%.4f s, avg_v=%.4f, v_inst=%.4f\n",
                               leftmost_wall_max_delta,
                               leftmost_wall_time_of_max,
                               leftmost_wall_average_speed,
                               leftmost_wall_velocity_at_max);
                        leftmost_wall_logged = true;
                    }

                    // bail-out if wall didn't get any motion for a long time
                    if (fabs(vx_wall) < 1e-15 && recorded_steps > wall_hold_steps * 10) {
                        printf("⚠️ Wall not moving — exiting early.\n");
                        break;
                    }
                }

                fclose(wall_log);
                printf("✅ Finished run %d\n\n", r);
            }
        }
    }
}









///////////////////// VISUAL SIMULATION FUNCTIONS //////////////////////
/*
 * simulation_loop
 * ----------------
 * Drives the interactive SDL visualisation:
 *  - waits for the wall release command, updates physics at a fixed dt, and
 *    renders particles, histograms, and HUD overlays each frame.
 *  - Logs energy and wall displacement traces for later analysis.
 *  - Surfaces collision diagnostics, piston work, and other statistics on a
 *    fixed cadence so long runs are easy to monitor from the terminal.
 */
void simulation_loop() {
    initialize_simulation();

    int running = 1;
    if (live_fft) {
        fft_cfg = kiss_fftr_alloc(FFT_SIZE, 0, NULL, NULL);
    }

    SDL_Event e;

    FILE *logFile = fopen("energy_log.csv", "w");
    FILE *wall_log = fopen("wall_position.csv", "w");
    if (!logFile || !wall_log) { printf("❌ Failed to open log files.\n"); exit(1); }

    if (log_packing_fraction) {
        log_packing_fractions(wall_log, "main");
    } else {
        fprintf(wall_log, "Time, Wall_X, Displacement(σ), Left_Count, Right_Count\n");
    }
    fprintf(logFile, "KE_particles,KE_wall,KE_total,T,S_v,F\n");

    // initial sim state
    steps_elapsed        = 0;
    wall_hold_enabled    = true;
    wall_is_released     = false;
    wall_enabled         = 1;               // ensure wall is drawn/active
    wall_x_old           = wall_x;          // first CCD step sees correct "old" pose
    missed_collision_events = 0;
    worst_penetration_observed = 0.0;

    Uint32 last_time = SDL_GetTicks();
    float accumulator = 0.0f;

    while (running) {
        Uint32 current_time = SDL_GetTicks();
        float frame_time = (current_time - last_time) / 1000.0f;
        last_time = current_time;
        if (frame_time > 0.1f) frame_time = 0.1f;
        accumulator += frame_time;

        // events
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT) running = 0;
            keysSimulation(e);
            keysPiston(e);
        }

        int left_particles = 0, right_particles = 0;

        while (accumulator >= FIXED_DT) {
            if (simulation_started && !paused) {
                // --- START OF STEP ---
                wall_x_old = wall_x;
                steps_elapsed++;

                // reset accumulator at the *start* of the step
                wall_impulse_x_accum = 0.f;

                // particles (this fills wall_impulse_x_accum via ccd_wall_step)
                update_particles_with_substepping(FIXED_DT, &left_particles, &right_particles);

                // wall velocity is already updated inside ccd_wall_step (impulses applied immediately).
                // wall_impulse_x_accum is kept for diagnostics only.

                update_pistons(FIXED_DT);
                update_wall(FIXED_DT, L0_UNITS);
                update_extra_walls(FIXED_DT, L0_UNITS);
                sync_all_wall_positions();
                update_wall_metrics();

                simulation_time += FIXED_DT;

                // Diagnostics (add piston work, hits)
                if ((steps_elapsed % 1000) == 0) {
                    double KEp = kinetic_energy();
                    double KEw = 0.5 * WALL_MASS * vx_wall * vx_wall;
                    double Wp  = piston_work_left + piston_work_right;
                    printf("[E] step=%d  KEp=%.6e  KEw=%.6e  KE=%.6e  Pwork=%.6e  (KE-Pwork)=%.6e  hits(L,R)=(%ld,%ld)\n",
                        steps_elapsed, KEp, KEw, KEp+KEw, Wp, (KEp+KEw)-Wp, piston_hits_left, piston_hits_right);
                    if (missed_collision_events > 0) {
                        printf("[CCD] step=%d  potential_misses=%ld  worst_penetration=%.6g\n",
                               steps_elapsed, missed_collision_events, worst_penetration_observed);
                    }
                    if (segment_counts && segment_count > 0) {
                        printf("[SEG] step=%d", steps_elapsed);
                        for (int seg = 0; seg < segment_count; ++seg) {
                            double ke = segment_ke ? segment_ke[seg] : 0.0;
                            double temp = segment_temperature ? segment_temperature[seg] : 0.0;
                            printf("  S%d=%d KE=%.3e T=%.3f", seg, segment_counts[seg], ke, temp);
                        }
                        printf("\n");
                    }
                    if (leftmost_wall_peak_recorded && !leftmost_wall_logged) {
                        printf("[WALL] first-peak Δx=%.4f px at t=%.4f s, avg_v=%.4f, v_inst=%.4f\n",
                               leftmost_wall_max_delta,
                               leftmost_wall_time_of_max,
                               leftmost_wall_average_speed,
                               leftmost_wall_velocity_at_max);
                        leftmost_wall_logged = true;
                    }
                }

                // inside the fixed-step block, after simulation_time += FIXED_DT
                if ((steps_elapsed % 100) == 0) {
                    log_energy(logFile);
                    fflush(logFile);
                }

                if ((steps_elapsed % 500) == 0) {
                double ke_p = kinetic_energy();
                double ke_w = 0.5 * WALL_MASS * vx_wall * vx_wall;
                printf("[ENERGY] step=%d  KE_p=%.6e  KE_w=%.6e  KE_tot=%.6e  vx_wall=%.6g\n",
                    steps_elapsed, ke_p, ke_w, ke_p + ke_w, vx_wall);
}

                if ((steps_elapsed % 1000) == 0) {
                    printf("KE=%.9f  vx_wall=%.6g  J_accum=%.6g\n",
                        kinetic_energy(), vx_wall, wall_impulse_x_accum);
                }

                // establish release time once
                if (wall_is_released && wall_release_time < 0.0) {
                    wall_release_time = simulation_time;
                    printf("✅ Wall released at t = %.3f\n", wall_release_time);
                }

                // logging after release
                if (wall_is_released) {
                    double time_after_release = simulation_time - wall_release_time;
                    double wall_x_sigma = wall_x / PIXELS_PER_SIGMA;
                    double disp = wall_x_sigma - (XW1 + XW2) / (2.0 * PIXELS_PER_SIGMA);

                    if (!log_packing_fraction) {
                        fprintf(wall_log, "%.3f, %.3f, %.3f, %d, %d\n",
                                time_after_release, wall_x_sigma, disp,
                                left_particles, right_particles);
                    } else {
                        float particle_area_unitless = (float)(M_PI * PARTICLE_RADIUS_UNIT * PARTICLE_RADIUS_UNIT);
                        float particle_area_actual   = (float)(M_PI * PARTICLE_RADIUS * PARTICLE_RADIUS);
                        float box_area_unitless      = 2.0f * L0_UNITS * HEIGHT_UNITS;

                        float L0_scaled      = L0_UNITS * PIXELS_PER_SIGMA;
                        float height_scaled  = HEIGHT_UNITS * PIXELS_PER_SIGMA;
                        float box_area_scaled= 2.0f * L0_scaled * height_scaled;
                        float wall_area      = 0.5f * WALL_THICKNESS * height_scaled;
                        float wall_buffer_area = 2.0f * wall_area;

                        float packing_unitless     = NUM_PARTICLES * particle_area_unitless / box_area_unitless;
                        float packing_actual_bound = NUM_PARTICLES * particle_area_actual / box_area_scaled;
                        float packing_actual_excl  = NUM_PARTICLES * particle_area_actual / (box_area_scaled - wall_buffer_area);

                        fprintf(wall_log, "%.3f, %.3f, %.3f, %d, %d, %.6f, %.6f, %.6f\n",
                                time_after_release, wall_x_sigma, disp,
                                left_particles, right_particles,
                                packing_unitless, packing_actual_bound, packing_actual_excl);
                    }
                }
            }
            accumulator -= FIXED_DT;
        }

        // === RENDERING ===
        draw_clear_screen();
        draw_coordinate_system(_renderer);
        draw_simulation_boundary();
        render_particles();
        render_pistons();
        if (wall_enabled) draw_wall();


        // <-- add this line
        render_energy_hud(_renderer, font);

        SDL_Color red = {255, 0, 0, 255};
        SDL_Color blue = {0, 0, 255, 255};
        SDL_Color yellow = {255, 255, 0, 255};

        char buffer[64];
        sprintf(buffer, "Left: %d", left_particles);
        draw_text(_renderer, font, buffer, XW1 + 5, YW1 + 10, red);
        sprintf(buffer, "Right: %d", right_particles);
        draw_text(_renderer, font, buffer, XW2 - 100, YW1 + 10, blue);

        if (wall_hold_enabled && !wall_is_released) {
            sprintf(buffer, "Wall held: %d steps left", wall_hold_steps - steps_elapsed);
            draw_text(_renderer, font, buffer, XW1 + 5, YW1 + 50, yellow);
        }

        render_velocity_histograms();
        float T_measured = compute_measured_temperature_from_ke();
        char temp_label[128];
        sprintf(temp_label, "T_measured: %.2f", T_measured);
        //draw_text(_renderer, font, temp_label, XW1 + 5, YW1 + HEIGHT_UNITS*PIXELS_PER_SIGMA, yellow);

        SDL_RenderPresent(_renderer);
        SDL_Delay(5);
    }

    fclose(logFile);
    fclose(wall_log);
}










/////////////// MAIN FUNCTION /////////////
// Main Function
int main(int argc, char* argv[]) {

    parse_cli_options(argc, argv);
    initialize_simulation_dimensions();     // sets SIM_WIDTH, SIM_HEIGHT, XW2 etc.
    initialize_time_scale();
    initialize_simulation_parameters();     // sets wall and piston regarding simulation dimensions
    print_simulation_info();  // Print simulation parameters

    initSDL();
    TTF_Init();
    //kiss_fftr_cfg fft_cfg = kiss_fftr_alloc(FFT_SIZE, 0, NULL, NULL);
    //initialize_simulation();
    
    font = TTF_OpenFont("Roboto-Regular.ttf", 18);
    if (!font) {
        printf("Failed to load font: %s\n", TTF_GetError());
        return 1;
    }

  
   // printf("wall_x = %f, wall_enabled = %d\n", wall_x, wall_enabled);
    if (enable_speed_of_sound_experiments) {
        run_speed_of_sound_experiments();
    } else {
        simulation_loop();  // your SDL interactive version
    }
    

    // Print the first few velocities for testing
    for (int i = 0; i < 10; i++) {
        printf("Particle %d: Vx = %f, Vy = %f\n", i, Vx[i], Vy[i]);
        
    }
    SDL_DestroyRenderer(_renderer);
    SDL_DestroyWindow(_window);
    //free(fft_cfg);
    TTF_CloseFont(font);
    TTF_Quit();
    SDL_Quit();
    free_internal_wall_arrays();
    free(cli_experiment_lengths); cli_experiment_lengths = NULL;
    free(cli_experiment_wall_masses); cli_experiment_wall_masses = NULL;
    return 0;
}












/////////////////////////////////////


void simulation_loop_old() {
    

    int running = 1;
    if (live_fft == true) {
        fft_cfg = kiss_fftr_alloc(FFT_SIZE, 0, NULL, NULL); // forward FFT
    }

    SDL_Event e;

    FILE *logFile = fopen("energy_log.csv", "w");
    if (!logFile) {
        printf("Error: Unable to open log file!\n");
        return;
    }
    fprintf(logFile, "Kinetic Energy, Temperature, Entropy, Total Free Energy\n");


    FILE *wall_log = fopen("wall_position.csv", "w");
    if (!wall_log) {
        printf("Error: Unable to open wall position log file!\n");
        return;
    }

    float particle_area_unitless = M_PI * PARTICLE_RADIUS_UNIT * PARTICLE_RADIUS_UNIT;
    float particle_area_actual = M_PI * PARTICLE_RADIUS * PARTICLE_RADIUS;
    float box_area = 2.0f * L0_UNITS * HEIGHT_UNITS;

    float packing_unit = NUM_PARTICLES * particle_area_unitless / box_area;
    float packing_actual = NUM_PARTICLES * particle_area_actual / box_area;

    fprintf(wall_log, "# PackingFraction_UNIT = %.6f\n", packing_unit);
    fprintf(wall_log, "# PackingFraction_ACTUAL = %.6f\n", packing_actual);
    fprintf(wall_log, "Time, Wall_X, Displacement(σ), Left_Count, Right_Count, Packing_UNITLESS, Packing_ACTUAL_BOUND, Packing_ACTUAL_WITH_WALL\n");
    //fprintf(wall_log, "Time, Wall_X, Displacement(σ), Left_Count, Right_Count\n");

    wall_hold_enabled = true;  
    wall_is_released = false;
    steps_elapsed = 0;

    Uint32 last_time = SDL_GetTicks();
    float accumulator = 0.0f;

    double time_after_release = 0.0;     // ← This is what we log

    while (running) {
        Uint32 current_time = SDL_GetTicks();
        float frame_time = (current_time - last_time) / 1000.0f;
        last_time = current_time;

        if (frame_time > 0.1f) frame_time = 0.1f;  
        accumulator += frame_time;

        // Poll events
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT) {
                running = 0;
            }
            keysSimulation(e);
            keysPiston(e);
        }




        int left_particles = 0;
        int right_particles = 0;

        // PHYSICS UPDATES
        while (accumulator >= FIXED_DT) {
            if (simulation_started == 1 && !paused) {
                update_particles_with_substepping(FIXED_DT, &left_particles, &right_particles);
                update_pistons(FIXED_DT);
                update_wall(FIXED_DT, L0_UNITS);
                //check_wall_collisions(FIXED_DT);

                simulation_time += FIXED_DT;  // ✅ advance the true simulation time

                // LIVE FFT
                if (live_fft == true) {
                    in_wall[sample_index] = wall_x;
                    in_energy[sample_index] = kinetic_energy();
                    sample_index = (sample_index + 1) % FFT_SIZE;

                    if (sample_index == 0) {
                        kiss_fftr(fft_cfg, in_wall, out_wall);
                        kiss_fftr(fft_cfg, in_energy, out_energy);

                        for (int i = 0; i < FFT_SIZE/2 + 1; i++) {
                            wall_fft_magnitude[i] = sqrtf(out_wall[i].r * out_wall[i].r + out_wall[i].i * out_wall[i].i);
                            energy_fft_magnitude[i] = sqrtf(out_energy[i].r * out_energy[i].r + out_energy[i].i * out_energy[i].i);
                        }
                    }
                }

                // Check if wall was just released
                if (wall_is_released && wall_release_time < 0.0) {
                    wall_release_time = simulation_time;
                    printf("✅ Wall released at t = %.6f s\n", wall_release_time);
                }

                // LOGGING ONLY AFTER WALL RELEASE
                if (wall_is_released) {
                    time_after_release = simulation_time - wall_release_time;

                    // Center of the box in σ units
                    float wall_x_sigma = wall_x / PIXELS_PER_SIGMA;
                    float center_sigma = (XW1 + XW2) / 2.0f / PIXELS_PER_SIGMA;
                    // Displacement in σ units
                    float displacement_centered = wall_x_sigma - center_sigma;
                    // Displacement in pixels
                    float displacement_pixels = wall_x - (L0_UNITS * PIXELS_PER_SIGMA);
                    float displacement_normalized = displacement_pixels / PIXELS_PER_SIGMA;

                    float wall_x_sigma_centered = (wall_x - XW1) / PIXELS_PER_SIGMA;
                    float wall_displacement = wall_x_sigma_centered - L0_UNITS;  // ⬅️ Δx from center at L0

                    fprintf(wall_log, "%.6f, %.10e, %.10e, %d, %d, %.6f, %.6f\n",
                        time_after_release,
                        wall_x / PIXELS_PER_SIGMA,
                        displacement_normalized,
                        left_particles,
                        right_particles,
                        packing_unit,
                        packing_actual);

                // Log energy always if needed
                log_energy(logFile);
            }
            accumulator -= FIXED_DT;
        }



        // RENDERING
        SDL_SetRenderDrawBlendMode(_renderer, SDL_BLENDMODE_BLEND);
        draw_clear_screen();
        draw_coordinate_system(_renderer);
        draw_simulation_boundary();
        render_velocity_histograms();
        render_particles();
        render_pistons();
        render_fft_histogram_bottom(_renderer, wall_fft_magnitude, NUM_BINS);

        if (wall_enabled) {
            draw_wall();
        }

        SDL_Color red = {255, 0, 0, 255};
        SDL_Color blue = {0, 0, 255, 255};
        SDL_Color yellow = {255, 255, 0, 255};

        char buffer[64];
        sprintf(buffer, "Left: %d", left_count);
        draw_text(_renderer, font, buffer, XW1 + 5, YW1 + 10, red);

        sprintf(buffer, "Right: %d", right_count);
        draw_text(_renderer, font, buffer, XW2 - 100, YW1 + 10, blue);

        if (wall_hold_enabled && !wall_is_released) {
            sprintf(buffer, "Wall held: %d steps left", wall_hold_steps - steps_elapsed);
            draw_text(_renderer, font, buffer, XW1 + 5, YW1 + 50, yellow);
        }



        SDL_RenderPresent(_renderer);
        SDL_Delay(10);  // throttle FPS
    }

    fclose(logFile);
    fclose(wall_log);
    }
}
