#include <SDL2/SDL.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <errno.h>
#include <stdlib.h>
#include <float.h>
#include <SDL2/SDL_ttf.h>
#include "kissfft/kiss_fft.h"
#include "kissfft/kiss_fftr.h"
#include "edmd_core/edmd.h"



/* 

00ALLINONE.c
   2023 Chris H. 
   All-in-one file for hard-disk piston simulation.
   Compile with:
   gcc -o hspist3 00ALLINONE.c -lm -lSDL2 -lSDL2_ttf

   TODO: Implement piston velocity control and dynamic time stepping
   implement SPRING instead of particles on left side

   BEST METHOD TO GET OUT WORK OF LEFT MOST WALL
    1. Track position of leftmost wall over time
    2. Calculate work done by integrating force over displacement
    3. Use piston velocity and pressure data for accuracy

*/


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
#define TIME_UNITS_SIMULATION 3000
#define SUBSTEPPING 1 // BOOL
#define INIT_SCALING 0 // BOOL
#define INIT_SCALING_WALL_RELEASE_STABLE 0 // BOOL

// PARTICLE PARAMETERS
// MAX buffer; runtime active count is configurable via --particles
#define NUM_PARTICLES 1000      // MAX_PARTICLES (buffer size)
static int particles_active = NUM_PARTICLES;  // active particles used by seeding/HUD/thermo
bool log_packing_fraction = 1; // 1 = log packing fraction, 0 = don't log




#if EXPERIMENT_MODE
    bool enable_speed_of_sound_experiments = true;
#else
    bool enable_speed_of_sound_experiments = false;
#endif



// --- Physical constants ---
#if MOLECULAR_MODE
    #define PARTICLE_RADIUS_UNIT .50f      // unit length               normalized radius = 1 (works in both modes!)
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
        #define FIXED_DT 0.001f  // in normalized time units They don’t state the exact dt, but normalized simulations use dt = 10^{-3} often.
        #define PIXELS_PER_SIGMA 40.0f      // rendering scaleq
        #define TEMPERATURE 1.0f            // default target temperature (reduced units)
        #define SUBSTEPS 8   // default, can be changed externally
        #define WALL_HOLD_STEPS 10000
        #define PARTICLE_SCALE_PARAM 1.0f // scaling factor for particle size



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


// Runtime-controllable particle size scaling (defaults to compile-time PARAM)
static float particle_scale_runtime = PARTICLE_SCALE_PARAM;
#define DIAMETER (2 * PARTICLE_RADIUS_UNIT)
#define PARTICLE_RADIUS (PARTICLE_RADIUS_UNIT * PIXELS_PER_SIGMA * particle_scale_runtime)

int num_steps = TIME_UNITS_SIMULATION / FIXED_DT;  // = 3,000,000 with 3000 units and 0.001 dt
bool wall_position_is_managed_externally = false;



// WALL PARAMETERS
#define WALL_MASS_FACTOR 200                         // choose factor relative to PARTICLE_MASS

/// SIMULATION  PARAMETERS
// Set these in main or before simulation starts
float L0_UNITS = 20.0f;       // 20 × radius = 10 × diameter Box half lenght paper 7.5 until 35
float HEIGHT_UNITS = 10.0f;   // 10 × radius = 5 × diameterx
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
    static float wall_thickness_runtime = 1.0f;   // fixed thin wall in experiment mode
#else
    // Note: depends on runtime particle radius; set default later in initialize_simulation_parameters()
    static float wall_thickness_runtime = 0.0f;   // will default to 2R at init if not overridden
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
    TAU_time_scale_factor_to_molecular = sqrtf((PARTICLE_MASS * powf(PARTICLE_RADIUS * METERS_PER_PIXEL, 2)) / (kBT_effective()));
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

// Track energy transferred by pistons since a push started (interactive mode)
static bool   piston_push_tracking     = false;
static double piston_work_baseline     = 0.0;   // baseline total piston work at push start
static double piston_work_delta_max    = 0.0;   // max (piston_work_total - baseline) since start

static inline void start_piston_push_tracking(void);

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
static float  cli_override_temperature = -1.0f;  // --temperature=VALUE (reduced units)
// Particle size CLI
static int    cli_particle_radius_set = 0;     // 1 if user set --particle-radius (sigma units)
static float  cli_particle_radius_sigma = 0.0f;
static int    cli_particle_diameter_set = 0;   // 1 if user set --particle-diameter (sigma units)
static float  cli_particle_diameter_sigma = 0.0f;
static bool   cli_force_no_experiments = false;
static bool   cli_enable_energy_measurement = false;
static int    cli_override_particles = -1;   // --particles=N
// Seeding mode CLI
typedef enum { SEEDING_GRID=0, SEEDING_HONEYCOMB=1, SEEDING_RANDOM=2 } seeding_mode_t;
static seeding_mode_t cli_seeding_mode = SEEDING_GRID; // default grid for reproducibility
static bool cli_left_empty = false;          // seed zero particles in leftmost segment
static bool cli_no_pp_collisions = false;    // disable particle-particle collisions
static bool cli_distribute_area = false;     // distribute particles by segment area
// Spring CLI overrides
static int   cli_spring_k_set = 0;          // 1 if user set --spring-k
static float cli_spring_k = 0.0f;
static int   cli_spring_eq_set = 0;         // 1 if user set --spring-eq (sigma units)
static float cli_spring_eq_sigma = 0.0f;     // in sigma units from left boundary
// Calibrate work: approximate placement of rightmost wall to match target
static int   cli_calibrate_work_set = 0;
static float cli_calib_target = 0.0f;        // target work (code units)
static float cli_calib_tol = 0.0f;           // not used in approx, reserved
static float cli_calib_min_sigma = 0.0f;     // min pos (sigma units)
static float cli_calib_max_sigma = 0.0f;     // max pos (sigma units)
// Particles per-segment specification
static float *cli_particles_box_fractions = NULL; // existing fractions list
static size_t cli_particles_box_fractions_count = 0;
static int   *cli_particles_box_counts = NULL;     // absolute counts list
static size_t cli_particles_box_counts_count = 0;
static bool   preset_custom_absolute_counts = false;
static float  cli_phi_max = 0.78f;          // densest packing cap

typedef enum {
    EXPERIMENT_PRESET_NONE = 0,
    EXPERIMENT_PRESET_OFFCENTER_SINGLE,
    EXPERIMENT_PRESET_TWO_WALLS,
    EXPERIMENT_PRESET_WALL_MID,
    // NEW RESEARCH EXPERIMENTS
    EXPERIMENT_PRESET_SPEED_OF_SOUND,
    EXPERIMENT_PRESET_ENERGY_TRANSFER,
    EXPERIMENT_PRESET_MULTI_WALL_SYSTEM,
    EXPERIMENT_PRESET_SZILARD_ENGINE,
    EXPERIMENT_PRESET_ATP_SYNTHASE_ANALOG,
    EXPERIMENT_PRESET_PROTOCOL_OPTIMIZATION,
    EXPERIMENT_PRESET_DYNAMIC_IB
} ExperimentPreset;

// Protocol types for different piston movements
typedef enum {
    PROTOCOL_STEP = 0,           // Instant acceleration (unphysical)
    PROTOCOL_SIGMOIDAL,         // Smooth sigmoidal acceleration
    PROTOCOL_LINEAR,            // Linear acceleration
    PROTOCOL_SINUSOIDAL,        // Sinusoidal protocol
    PROTOCOL_OPTIMAL              // AI-optimized protocol
} PistonProtocol;

// Multi-wall configuration
typedef struct {
    int num_walls;                              // 1-5 walls
    float wall_positions[5];                    // x positions (left to right)
    float wall_masses[5];                       // mass factors (left to right)
    PistonProtocol protocol;                    // Movement protocol
    float protocol_duration;                    // Time for protocol
    float max_velocity;                        // Max piston velocity
    bool energy_tracking;                      // Track energy transfer
    bool information_tracking;                // Track information entropy
} MultiWallConfig;

static ExperimentPreset cli_experiment_preset = EXPERIMENT_PRESET_NONE;

static bool  preset_custom_positions = false;
static int   preset_wall_count = 0;
static float preset_wall_fraction[MAX_DIVIDER_CAPACITY];
static bool  preset_custom_fractions = false;
static float preset_segment_fraction[MAX_DIVIDER_CAPACITY + 1];
static int   preset_two_wall_active_walls = 2;  // controls variant for two-wall preset
// CLI overrides for segment particle allocation
static int    cli_particles_box_left = -1;  // for 2-box shorthand
static int    cli_particles_box_right = -1;

// Multi-wall mass parameters (left to right, up to 5 walls)
static float cli_wall_masses[5] = {200.0f, 200.0f, 200.0f, 200.0f, 200.0f}; // default mass factors
static int   cli_num_walls = 1;  // default number of walls
static float cli_wall_positions[5] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f}; // wall positions
static PistonProtocol cli_protocol = PROTOCOL_SIGMOIDAL; // default protocol
// Optional: per-segment temperature overrides
static float *cli_temperature_segments = NULL;        // from --temperature list
static size_t cli_temperature_segments_count = 0;

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
    // Set a sane default wall thickness if not explicitly provided (e.g., via CLI)
    // Default to 2×particle radius in pixels to avoid tunneling
#if !EXPERIMENT_MODE
    if (wall_thickness_runtime <= 0.0f) {
        wall_thickness_runtime = 2.0f * PARTICLE_RADIUS;
    }
#endif
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

// Forward declarations for helpers referenced in temperature tools
static int segment_index_for_position(double px);
double kinetic_energy(void);

// Runtime temperature control (defaults to compile-time TEMPERATURE)
static float temperature_runtime = TEMPERATURE;
static float fixed_dt_runtime = FIXED_DT; // dt actually used by the integrator
static float time_scale_runtime = 1.0f;   // scales real-time → sim-time in interactive mode
// Integrator / simulation mode
typedef enum { MODE_TIME=0, MODE_RK4=1, MODE_EDMD=2, MODE_EDMD_HYBRID=3 } sim_mode_t;
static sim_mode_t sim_mode = MODE_TIME;
// EDMD state (optional)
static EDMD* g_edmd = NULL;
static int   cli_force_kbt_one = 0; // --kbt1: force K_B*T == 1 (reduced units)

static inline float kB_effective(void) {
    // If user forces K_B*T=1, adjust k_B at runtime so k_B*T stays 1
    if (cli_force_kbt_one) {
        float T = (temperature_runtime > 1e-12f) ? temperature_runtime : 1.0f;
        return 1.0f / T;
    }
    return (float)K_B;
}
static inline float kBT_effective(void) {
    return kB_effective() * temperature_runtime;
}

// Particle color-coding mode
typedef enum { COLORCODE_NONE=0, COLORCODE_TEMPERATURE=1, COLORCODE_VELOCITY=2 } ColorCodeMode;
static ColorCodeMode colorcode_mode = COLORCODE_NONE; // set via --colourcode / --colorcode

// Recompute runtime dt to preserve v·dt/σ when temperature changes (speed_of_sound experiment)
static inline void update_dt_runtime_for_temperature(void) {
    // Default: keep compile-time dt
    fixed_dt_runtime = FIXED_DT;
    // Scale dt with 1/sqrt(T) only for speed_of_sound preset and when not forcing kB*T=1
    if (cli_experiment_preset == EXPERIMENT_PRESET_SPEED_OF_SOUND && !cli_force_kbt_one) {
        float T = fmaxf(temperature_runtime, 1e-12f);
        fixed_dt_runtime = FIXED_DT / sqrtf(T);
    }
}

/* Recompute per-segment counts, KE, and temperature from current arrays. */
static void recompute_segment_stats_counts_and_temperature(void) {
    if (!(segment_count > 0 && segment_count <= MAX_DIVIDER_CAPACITY + 1)) return;
    int n = particles_active > 0 ? particles_active : 0;
    if (!segment_counts) return;
    memset(segment_counts, 0, (size_t)segment_count * sizeof(int));
    if (segment_ke) memset(segment_ke, 0, (size_t)segment_count * sizeof(double));
    if (segment_temperature) memset(segment_temperature, 0, (size_t)segment_count * sizeof(double));
    for (int i = 0; i < n; ++i) {
        int seg = segment_index_for_position(X[i]);
        if (seg < 0) seg = 0; if (seg >= segment_count) seg = segment_count - 1;
        segment_counts[seg]++;
        if (segment_ke) {
            double ke = 0.5 * PARTICLE_MASS * (Vx[i] * Vx[i] + Vy[i] * Vy[i]);
            segment_ke[seg] += ke;
        }
    }
    if (segment_temperature && segment_ke) {
        for (int seg = 0; seg < segment_count; ++seg) {
            if (segment_counts[seg] > 0) {
                segment_temperature[seg] = segment_ke[seg] / (segment_counts[seg] * K_B);
            } else {
                segment_temperature[seg] = 0.0;
            }
        }
    }
}
static inline void rescale_all_velocities_to_temperature(float T_new) {
    int n = particles_active > 0 ? particles_active : 0;
    if (n <= 0) return;
    double ke_now = kinetic_energy();
    if (ke_now <= 0.0) return;
    double ke_target = (double)n * (double)kB_effective() * (double)T_new;
    double scale = sqrt(ke_target / ke_now);
    for (int i = 0; i < n; ++i) { Vx[i] *= scale; Vy[i] *= scale; }
    temperature_runtime = T_new;
}

static inline void remove_drift_segment(int seg_idx) {
    if (!segment_count || seg_idx < 0 || seg_idx >= segment_count) return;
    int n = particles_active > 0 ? particles_active : 0;
    double sx = 0.0, sy = 0.0; int c = 0;
    for (int i = 0; i < n; ++i) {
        int s = segment_index_for_position(X[i]);
        if (s == seg_idx) { sx += Vx[i]; sy += Vy[i]; c++; }
    }
    if (c <= 0) return;
    double mx = sx / c, my = sy / c;
    for (int i = 0; i < n; ++i) {
        int s = segment_index_for_position(X[i]);
        if (s == seg_idx) { Vx[i] -= mx; Vy[i] -= my; }
    }
}

static inline void equalize_temperature_per_segment(float T_target) {
    if (segment_count <= 0) return;
    int n = particles_active > 0 ? particles_active : 0;
    for (int seg = 0; seg < segment_count; ++seg) {
        // compute seg KE and count
        double ke_seg = 0.0; int c = 0;
        for (int i = 0; i < n; ++i) {
            int s = segment_index_for_position(X[i]);
            if (s == seg) { ke_seg += 0.5 * PARTICLE_MASS * (Vx[i]*Vx[i] + Vy[i]*Vy[i]); c++; }
        }
        if (c <= 0 || ke_seg <= 0.0) continue;
        double T_seg = ke_seg / (c * K_B);
        if (T_seg > 0.0) {
            double scale = sqrt((double)T_target / T_seg);
            for (int i = 0; i < n; ++i) {
                int s = segment_index_for_position(X[i]);
                if (s == seg) { Vx[i] *= scale; Vy[i] *= scale; }
            }
        }
        remove_drift_segment(seg);
    }
}

// Apply per-segment temperatures (left→right). If temps_count < segment_count,
// remaining segments keep current temperature. Uses kB_effective at runtime.
static void apply_segment_temperatures(const float* temps, size_t temps_count) {
    if (!temps || temps_count == 0) return;
    if (segment_count <= 0) return;
    int n = particles_active > 0 ? particles_active : 0;
    if (n <= 0) return;

    // Accumulate KE and counts per segment
    double* ke = (double*)malloc((size_t)segment_count * sizeof(double));
    int* cnt   = (int*)malloc((size_t)segment_count * sizeof(int));
    if (!ke || !cnt) { if (ke) free(ke); if (cnt) free(cnt); return; }
    for (int s = 0; s < segment_count; ++s) { ke[s] = 0.0; cnt[s] = 0; }
    for (int i = 0; i < n; ++i) {
        int s = segment_index_for_position(X[i]);
        if (s < 0) s = 0; if (s >= segment_count) s = segment_count - 1;
        double v2 = (double)Vx[i]*(double)Vx[i] + (double)Vy[i]*(double)Vy[i];
        ke[s] += 0.5 * (double)PARTICLE_MASS * v2;
        cnt[s] += 1;
    }

    // Compute scale per segment and apply
    for (int s = 0; s < segment_count; ++s) {
        float T_target = (s < (int)temps_count) ? temps[s] : temperature_runtime;
        if (cnt[s] > 0 && ke[s] > 0.0 && T_target > 0.0f) {
            double kB = (double)kB_effective();
            double ke_target = (double)cnt[s] * kB * (double)T_target;
            double scale = sqrt(ke_target / ke[s]);
            for (int i = 0; i < n; ++i) {
                int seg_i = segment_index_for_position(X[i]);
                if (seg_i < 0) seg_i = 0; if (seg_i >= segment_count) seg_i = segment_count - 1;
                if (seg_i == s) { Vx[i] = (float)((double)Vx[i] * scale); Vy[i] = (float)((double)Vy[i] * scale); }
            }
            remove_drift_segment(s);
        }
    }
    free(ke); free(cnt);
}

// Energy measurement system
typedef struct {
    float spring_constant;           // Spring constant for leftmost wall
    float equilibrium_position;      // Equilibrium position of spring
    float spring_force;             // Current spring force
    float spring_energy;            // Potential energy in spring
    float spring_energy_max;        // Peak potential energy observed
    float energy_transferred;       // Total energy transferred to leftmost wall
    float max_displacement;         // Maximum displacement from equilibrium
    bool spring_enabled;            // Enable spring system for energy measurement
} EnergyMeasurement;

// Forward declaration for energy measurement initializer
void initialize_energy_measurement(void);

static EnergyMeasurement energy_measurement = {
    .spring_constant = 1000.0f,      // Strong spring for good measurement
    .equilibrium_position = 0.0f,    // Will be set to initial wall position
    .spring_force = 0.0f,
    .spring_energy = 0.0f,
    .spring_energy_max = 0.0f,
    .energy_transferred = 0.0f,
    .max_displacement = 0.0f,
    .spring_enabled = false
};

// Define after energy_measurement is declared
static inline void start_piston_push_tracking(void) {
    piston_push_tracking  = true;
    piston_work_baseline  = piston_work_left + piston_work_right;
    piston_work_delta_max = 0.0;
    // Only initialize energy measurement if explicitly requested or in energy_transfer preset
    if (!energy_measurement.spring_enabled && num_internal_walls > 0) {
        if (cli_enable_energy_measurement || cli_experiment_preset == EXPERIMENT_PRESET_ENERGY_TRANSFER) {
            initialize_energy_measurement();
        }
    }
}
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
    piston_push_tracking = false;
    piston_work_baseline = 0.0;
    piston_work_delta_max = 0.0;
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
    printf("  --experiment=NAME           Use a preset:\n");
    printf("                               speed_of_sound     - Standard speed of sound experiment\n");
    printf("                               energy_transfer    - Energy transfer efficiency (right→left)\n");
    printf("\nInteractive control\n");
    printf("  --show-simulation           Open SDL window (alias of --no-experiments)\n");
    printf("  --left-empty                Seed zero particles in the leftmost segment\n");
    printf("  --no-pp-collisions          Disable particle–particle collisions (ideal gas limit)\n");
    printf("                               multi_wall         - Multi-wall system analysis\n");
    printf("                               szilard_engine     - Szilard engine entropy experiments\n");
    printf("                               atp_synthase       - ATP synthase analog experiments\n");
    printf("                               protocol_optimization - Find optimal piston protocols\n");
    printf("                               dynamic_ib         - Dynamic information bottleneck\n");
    printf("                               wall_mid           - Single wall in middle\n");
    printf("                               2wall_exp_with_1_wall, 2wall_exp_with_2_wall\n");
    printf("  --lengths=L0,...            Override experiment half-lengths (sigma units)\n");
    printf("  --wall-masses=f1,...        Override experiment wall-mass factors (multiples of particle mass)\n");
    printf("  --num-walls=N               Number of walls (1-5, default 1)\n");
    printf("  --wall-positions=x1,x2,...  Wall positions from left to right (sigma units)\n");
    printf("  --wall-mass-factors=m1,m2,... Wall mass factors from left to right (default 200 each)\n");
    printf("  --seeding=MODE              Seeding mode: grid (default), honeycomb, random\n");
    printf("  --distribute=area           Distribute particles ∝ segment area (width×height)\n");
    printf("  --protocol=NAME             Piston protocol: step, sigmoidal, linear, sinusoidal, optimal\n");
    printf("  --spring-k=value            Spring constant for energy measurement (code units)\n");
    printf("  --spring-eq=value          Spring equilibrium position (sigma units from left)\n");
    printf("  --calibrate-work=target,tol,min,max  Approximate place rightmost wall to match target piston work\n");
    printf("  --repeats=N                 Number of repeats per experiment configuration\n");
    printf("  --steps=N                   Steps to simulate after wall release (default derived from TIME_UNITS_SIMULATION)\n");
    printf("  --l0=value                  Set initial half-length for interactive mode\n");
    printf("  --height=value              Set box height (in sigma units) for interactive mode\n");
    printf("  --wall-mass-factor=value    Set wall mass factor for interactive mode (M_wall = factor * m_particle)\n");
    printf("  --wall-thickness=value      Set wall thickness in sigma units\n");
    printf("  --particle-radius=value     Set particle radius in σ units (default 0.5)\n");
    printf("  --particle-diameter=value   Set particle diameter in σ units\n");
    printf("  --particles=N               Total active particles (<= MAX buffer)\n");
    printf("  --phi-max=value             Max packing fraction cap for seeding (default 0.78)\n");
    printf("  --particles-boxes=a,b,...   Target counts per segment (left→right)\n");
    printf("  --particles-box-left=n      Shorthand: left count for 2 boxes\n");
    printf("  --particles-box-right=n     Shorthand: right count for 2 boxes\n");
    printf("  --energy-measurement        Enable spring-based energy measurement\n");
    printf("  --colourcode=mode           Particle coloring: none (default), temperature, velocity\n");
    printf("  --colorcode=mode            Alias of --colourcode\n");
    printf("  --timescale=value           Speed multiplier for interactive run (default 1.0)\n");
    printf("  --kbt1                      Force k_B*T = 1 (reduced units); adjusts k_B at runtime\n");
    printf("  --mode=NAME                 time (default) | rk4 | edmd | edmd-hybrid\n");
    printf("  --help                      Show this help message and exit\n");
}

static void parse_cli_options(int argc, char **argv) {
    preset_two_wall_active_walls = 2;  // default to both internal walls active for preset
    for (int i = 1; i < argc; ++i) {
        const char *arg = argv[i];

        if (strcmp(arg, "--help") == 0) {
            print_cli_usage(argv[0]);
            exit(EXIT_SUCCESS);
        } else if (strcmp(arg, "--experiments") == 0) {
            enable_speed_of_sound_experiments = true;
        } else if (strcmp(arg, "--no-experiments") == 0) {
            cli_force_no_experiments = true;
        } else if (strcmp(arg, "--show-simulation") == 0) {
            // Alias: explicitly run interactive SDL and do not start batch experiments
            cli_force_no_experiments = true;
            enable_speed_of_sound_experiments = false;
        } else if (strcmp(arg, "--left-empty") == 0) {
            cli_left_empty = true;
        } else if (strcmp(arg, "--no-pp-collisions") == 0) {
            cli_no_pp_collisions = true;
        } else if (strncmp(arg, "--distribute", 12) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            if (strcmp(value, "area") == 0) cli_distribute_area = 1; else if (strcmp(value, "equal") == 0) cli_distribute_area = 0; else {
                fprintf(stderr, "Unknown distribute mode '%s'. Use: area or equal\n", value); exit(EXIT_FAILURE);
            }
        } else if (strcmp(arg, "--energy-measurement") == 0) {
            cli_enable_energy_measurement = true;
        } else if (strcmp(arg, "--kbt1") == 0) {
            cli_force_kbt_one = 1;
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
            } else if (strcmp(preset_name, "wall_mid") == 0) {
                cli_experiment_preset = EXPERIMENT_PRESET_WALL_MID;
                enable_speed_of_sound_experiments = true;
            } else if (strcmp(preset_name, "speed_of_sound") == 0) {
                cli_experiment_preset = EXPERIMENT_PRESET_SPEED_OF_SOUND;
                enable_speed_of_sound_experiments = true;
        } else if (strcmp(preset_name, "energy_transfer") == 0) {
            cli_experiment_preset = EXPERIMENT_PRESET_ENERGY_TRANSFER;
            enable_speed_of_sound_experiments = true;
            // Default behavior for energy transfer: left box empty and enable spring energy measurement
            cli_left_empty = true;
            cli_enable_energy_measurement = true;
        } else if (strcmp(preset_name, "multi_wall") == 0) {
                cli_experiment_preset = EXPERIMENT_PRESET_MULTI_WALL_SYSTEM;
                enable_speed_of_sound_experiments = true;
            } else if (strcmp(preset_name, "szilard_engine") == 0) {
                cli_experiment_preset = EXPERIMENT_PRESET_SZILARD_ENGINE;
                enable_speed_of_sound_experiments = true;
            } else if (strcmp(preset_name, "atp_synthase") == 0) {
                cli_experiment_preset = EXPERIMENT_PRESET_ATP_SYNTHASE_ANALOG;
                enable_speed_of_sound_experiments = true;
            } else if (strcmp(preset_name, "protocol_optimization") == 0) {
                cli_experiment_preset = EXPERIMENT_PRESET_PROTOCOL_OPTIMIZATION;
                enable_speed_of_sound_experiments = true;
            } else if (strcmp(preset_name, "dynamic_ib") == 0) {
                cli_experiment_preset = EXPERIMENT_PRESET_DYNAMIC_IB;
                enable_speed_of_sound_experiments = true;
            } else {
                fprintf(stderr, "Unknown experiment preset '%s'.\n", value_raw);
                print_cli_usage(argv[0]);
                exit(EXIT_FAILURE);
            }

        } else if (strncmp(arg, "--particles-boxes", strlen("--particles-boxes")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            // Parse as floats first
            cli_parse_float_list(value, &cli_particles_box_fractions, &cli_particles_box_fractions_count);
            // Detect absolute counts if any entry > 1.0
            bool any_gt_one = false;
            for (size_t k = 0; k < cli_particles_box_fractions_count; ++k) {
                if (cli_particles_box_fractions[k] > 1.0f + 1e-6f) { any_gt_one = true; break; }
            }
            if (any_gt_one) {
                // Convert to integer counts
                size_t n = cli_particles_box_fractions_count;
                free(cli_particles_box_counts);
                cli_particles_box_counts = (int*)malloc(n * sizeof(int));
                cli_particles_box_counts_count = n;
                for (size_t k = 0; k < n; ++k) cli_particles_box_counts[k] = (int)llround(cli_particles_box_fractions[k]);
                preset_custom_absolute_counts = true;
                preset_custom_fractions = false;
            } else {
                preset_custom_fractions = true;
                preset_custom_absolute_counts = false;
            }
        } else if (strncmp(arg, "--particles-box-left", strlen("--particles-box-left")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0; char *endptr = NULL; long v = strtol(value, &endptr, 10);
            if (errno != 0 || endptr == value) { fprintf(stderr, "Invalid left box count '%s'.\n", value); exit(EXIT_FAILURE);} 
            cli_particles_box_left = (int)v; preset_custom_fractions = true;
        } else if (strncmp(arg, "--particles-box-right", strlen("--particles-box-right")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0; char *endptr = NULL; long v = strtol(value, &endptr, 10);
            if (errno != 0 || endptr == value) { fprintf(stderr, "Invalid right box count '%s'.\n", value); exit(EXIT_FAILURE);} 
            cli_particles_box_right = (int)v; preset_custom_fractions = true;
        } else if (strncmp(arg, "--num-walls", strlen("--num-walls")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0; char *endptr = NULL; long v = strtol(value, &endptr, 10);
            if (errno != 0 || endptr == value || v < 1 || v > 5) { 
                fprintf(stderr, "Invalid number of walls '%s'. Use 1-5.\n", value); 
                exit(EXIT_FAILURE);
            } 
            cli_num_walls = (int)v;
        } else if (strncmp(arg, "--wall-positions", strlen("--wall-positions")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            // Parse comma-separated positions into cli_wall_positions array
            char *value_copy = strdup(value);
            char *token = strtok(value_copy, ",");
            int pos_count = 0;
            while (token && pos_count < 5) {
                // Skip whitespace
                while (*token == ' ') token++;
                cli_wall_positions[pos_count] = atof(token);
                token = strtok(NULL, ",");
                pos_count++;
            }
            free(value_copy);
        } else if (strncmp(arg, "--wall-mass-factors", strlen("--wall-mass-factors")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            // Parse comma-separated masses into cli_wall_masses array
            char *value_copy = strdup(value);
            char *token = strtok(value_copy, ",");
            int mass_count = 0;
            while (token && mass_count < 5) {
                // Skip whitespace
                while (*token == ' ') token++;
                cli_wall_masses[mass_count] = atof(token);
                token = strtok(NULL, ",");
                mass_count++;
            }
            free(value_copy);
        } else if (strncmp(arg, "--seeding", strlen("--seeding")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            if (strcmp(value, "grid") == 0) cli_seeding_mode = SEEDING_GRID;
            else if (strcmp(value, "honeycomb") == 0) cli_seeding_mode = SEEDING_HONEYCOMB;
            else if (strcmp(value, "random") == 0) cli_seeding_mode = SEEDING_RANDOM;
            else {
                fprintf(stderr, "Unknown seeding mode '%s'. Use: grid, honeycomb, random\n", value);
                exit(EXIT_FAILURE);
            }
        } else if (strncmp(arg, "--spring-k", 10) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0; char *endptr = NULL; float v = strtof(value, &endptr);
            if (errno != 0 || endptr == value) { fprintf(stderr, "Invalid spring-k '%s'.\n", value); exit(EXIT_FAILURE);} 
            cli_spring_k_set = 1; cli_spring_k = v;
        } else if (strncmp(arg, "--spring-eq", 11) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0; char *endptr = NULL; float v = strtof(value, &endptr);
            if (errno != 0 || endptr == value) { fprintf(stderr, "Invalid spring-eq '%s'.\n", value); exit(EXIT_FAILURE);} 
            cli_spring_eq_set = 1; cli_spring_eq_sigma = v;
        } else if (strncmp(arg, "--calibrate-work", 16) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            // Parse target,tol,min,max
            float list[4]; size_t n = 0; char *copy = cli_strdup(value); char *tok = strtok(copy, ",");
            while (tok && n < 4) { list[n++] = strtof(tok, NULL); tok = strtok(NULL, ","); }
            free(copy);
            if (n != 4) { fprintf(stderr, "--calibrate-work expects 4 values: target,tol,min,max\n"); exit(EXIT_FAILURE);} 
            cli_calibrate_work_set = 1; cli_calib_target = list[0]; cli_calib_tol = list[1]; cli_calib_min_sigma = list[2]; cli_calib_max_sigma = list[3];
        } else if (strncmp(arg, "--protocol", strlen("--protocol")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            if (strcmp(value, "step") == 0) cli_protocol = PROTOCOL_STEP;
            else if (strcmp(value, "sigmoidal") == 0) cli_protocol = PROTOCOL_SIGMOIDAL;
            else if (strcmp(value, "linear") == 0) cli_protocol = PROTOCOL_LINEAR;
            else if (strcmp(value, "sinusoidal") == 0) cli_protocol = PROTOCOL_SINUSOIDAL;
            else if (strcmp(value, "optimal") == 0) cli_protocol = PROTOCOL_OPTIMAL;
            else {
                fprintf(stderr, "Unknown protocol '%s'. Use: step, sigmoidal, linear, sinusoidal, optimal\n", value);
                exit(EXIT_FAILURE);
            }
        } else if (strncmp(arg, "--particle-radius", strlen("--particle-radius")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0; char *endptr = NULL; float v = strtof(value, &endptr);
            if (errno != 0 || endptr == value || *endptr != '\0' || v <= 0.0f) {
                fprintf(stderr, "Invalid particle radius '%s' (σ units).\n", value);
                exit(EXIT_FAILURE);
            }
            cli_particle_radius_set = 1;
            cli_particle_radius_sigma = v;
        } else if (strncmp(arg, "--particle-diameter", strlen("--particle-diameter")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0; char *endptr = NULL; float v = strtof(value, &endptr);
            if (errno != 0 || endptr == value || *endptr != '\0' || v <= 0.0f) {
                fprintf(stderr, "Invalid particle diameter '%s' (σ units).\n", value);
                exit(EXIT_FAILURE);
            }
            cli_particle_diameter_set = 1;
            cli_particle_diameter_sigma = v;
        } else if (strncmp(arg, "--particles-boxes", strlen("--particles-boxes")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            cli_parse_float_list(value, &cli_particles_box_fractions, &cli_particles_box_fractions_count);
            preset_custom_fractions = true;
        } else if (strncmp(arg, "--particles-box-left", strlen("--particles-box-left")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0; char *endptr = NULL; long v = strtol(value, &endptr, 10);
            if (errno != 0 || endptr == value) { fprintf(stderr, "Invalid left box count '%s'.\n", value); exit(EXIT_FAILURE);} 
            cli_particles_box_left = (int)v; preset_custom_fractions = true;
        } else if (strncmp(arg, "--particles-box-right", strlen("--particles-box-right")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0; char *endptr = NULL; long v = strtol(value, &endptr, 10);
            if (errno != 0 || endptr == value) { fprintf(stderr, "Invalid right box count '%s'.\n", value); exit(EXIT_FAILURE);} 
            cli_particles_box_right = (int)v; preset_custom_fractions = true;

        } else if (strncmp(arg, "--particles", strlen("--particles")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0; char *endptr = NULL; long v = strtol(value, &endptr, 10);
            if (errno != 0 || endptr == value || *endptr != '\0' || v <= 0) {
                fprintf(stderr, "Invalid particles value '%s'.\n", value);
                exit(EXIT_FAILURE);
            }
            cli_override_particles = (int)v;
        } else if (strncmp(arg, "--phi-max", strlen("--phi-max")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0; char *endptr = NULL; float v = strtof(value, &endptr);
            if (errno != 0 || endptr == value || *endptr != '\0' || v <= 0.0f || v > 0.9069f) {
                fprintf(stderr, "Invalid phi-max '%s' (0<phi<=0.9069).\n", value);
                exit(EXIT_FAILURE);
            }
            cli_phi_max = v;
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
        } else if (strncmp(arg, "--wall-mass-factor", strlen("--wall-mass-factor")) == 0) {
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
        } else if (strncmp(arg, "--temperature", 13) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            if (strchr(value, ',') != NULL) {
                // parse per-segment list
                cli_parse_float_list(value, &cli_temperature_segments, &cli_temperature_segments_count);
                cli_override_temperature = -1.0f; // ignore single override in favor of list
            } else {
                errno = 0; char *endptr = NULL; float v = strtof(value, &endptr);
                if (errno != 0 || endptr == value || *endptr != '\0' || v <= 0.0f) {
                    fprintf(stderr, "Invalid temperature '%s'.\n", value);
                    exit(EXIT_FAILURE);
                }
                cli_override_temperature = v;
            }
        } else if (strncmp(arg, "--timescale", 11) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0; char *endptr = NULL; float v = strtof(value, &endptr);
            if (errno != 0 || endptr == value || *endptr != '\0' || v <= 0.0f) {
                fprintf(stderr, "Invalid timescale '%s'.\n", value);
                exit(EXIT_FAILURE);
            }
            time_scale_runtime = v;
        } else if (strncmp(arg, "--mode", 6) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            if (strcmp(value, "time") == 0) sim_mode = MODE_TIME;
            else if (strcmp(value, "rk4") == 0) sim_mode = MODE_RK4;
            else if (strcmp(value, "edmd") == 0) sim_mode = MODE_EDMD;
            else if (strcmp(value, "edmd-hybrid") == 0) sim_mode = MODE_EDMD_HYBRID;
            else { fprintf(stderr, "Unknown mode '%s'. Use: time, rk4, edmd\n", value); exit(EXIT_FAILURE);} 
        } else if (strncmp(arg, "--colourcode", 12) == 0 || strncmp(arg, "--colorcode", 11) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            if (strcmp(value, "none") == 0) colorcode_mode = COLORCODE_NONE;
            else if (strcmp(value, "temperature") == 0) colorcode_mode = COLORCODE_TEMPERATURE;
            else if (strcmp(value, "velocity") == 0) colorcode_mode = COLORCODE_VELOCITY;
            else { fprintf(stderr, "Unknown colourcode '%s'. Use: none | temperature | velocity\n", value); exit(EXIT_FAILURE);} 
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
                preset_wall_fraction[0] = 1.0f / 3.0f; // divider at 1/3 from left
                preset_segment_fraction[0] = 1.0f / 3.0f; // left share
                preset_segment_fraction[1] = 2.0f / 3.0f; // right share
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
        case EXPERIMENT_PRESET_WALL_MID:
            cli_requested_walls = 1;
            preset_custom_positions = true;
            preset_wall_count = 1;
            preset_wall_fraction[0] = 0.5f;
            preset_custom_fractions = true;
            preset_segment_fraction[0] = 0.5f;
            preset_segment_fraction[1] = 0.5f;
            break;
        case EXPERIMENT_PRESET_NONE:
        default:
            break;
    }

    // Apply explicit particle-box allocations if provided
    int segs = num_internal_walls + 1; if (segs < 1) segs = 1;
    if (cli_particles_box_left >= 0 || cli_particles_box_right >= 0) {
        if (segs >= 2) {
            double left = cli_particles_box_left >= 0 ? (double)cli_particles_box_left : 0.0;
            double right = cli_particles_box_right >= 0 ? (double)cli_particles_box_right : 0.0;
            if (left <= 0.0 && right <= 0.0) { left = right = 1.0; }
            double sum = left + right; if (sum <= 0.0) sum = 1.0;
            preset_segment_fraction[0] = (float)(left / sum);
            preset_segment_fraction[1] = (float)(right / sum);
            preset_custom_fractions = true;
        }
    } else if (cli_particles_box_fractions_count > 0) {
        double sum = 0.0; for (size_t k = 0; k < cli_particles_box_fractions_count; ++k) sum += cli_particles_box_fractions[k];
        if (sum <= 0.0) sum = 1.0;
        for (int s = 0; s < segs; ++s) {
            double v = 0.0; if (s < (int)cli_particles_box_fractions_count) v = cli_particles_box_fractions[s] / sum;
            preset_segment_fraction[s] = (float)v;
        }
        double acc = 0.0; for (int s = 0; s < segs; ++s) acc += preset_segment_fraction[s];
        if (acc < 1.0 && segs > 0) preset_segment_fraction[segs - 1] += (float)(1.0 - acc);
        preset_custom_fractions = true;
    }

    // Apply CLI wall count override
    if (cli_num_walls > 0) {
        cli_requested_walls = cli_num_walls;
    }
    
    num_internal_walls = cli_requested_walls;
    if (num_internal_walls > MAX_DIVIDER_CAPACITY) {
        num_internal_walls = MAX_DIVIDER_CAPACITY;
    }
    
    // Apply custom wall positions if provided (values in sigma units from left)
    if (cli_num_walls > 0 && cli_wall_positions[0] != 0.0f) {
        preset_custom_positions = true;
        preset_wall_count = cli_num_walls;
        for (int w = 0; w < cli_num_walls && w < MAX_DIVIDER_CAPACITY; w++) {
            float pos_sigma = cli_wall_positions[w];
            float span_sigma = 2.0f * L0_UNITS; // total width in sigma units
            float frac = pos_sigma / fmaxf(1e-6f, span_sigma);
            // clamp to [0,1]
            if (frac < 0.0f) frac = 0.0f;
            if (frac > 1.0f) frac = 1.0f;
            preset_wall_fraction[w] = frac;
        }
    }

    // Optional: approximate calibration to match target piston work by placing the rightmost wall
    // Assumptions: 2 internal walls → 3 segments; we place the last internal wall so that
    // initial rightmost width w satisfies W_target ≈ N_right * kT * (Δx/w). Here Δx ≈ 0.2 * box_width.
    if (cli_calibrate_work_set && cli_requested_walls >= 2) {
        // Estimate stroke in sigma units (uses same 20% travel used by 't' hotkey)
        float box_width_px = (float)SIM_WIDTH;
        float stroke_px = 0.20f * box_width_px;
        float stroke_sigma = stroke_px / PIXELS_PER_SIGMA;

        // Determine how many particles intended for rightmost box
        int segments = cli_requested_walls + 1;
        int Nreq = (cli_override_particles > 0) ? cli_override_particles : NUM_PARTICLES;
        int N_right = 0;
        if (preset_custom_absolute_counts && (int)cli_particles_box_counts_count == segments) {
            N_right = cli_particles_box_counts[segments - 1];
        } else if (preset_custom_fractions) {
            float frac_sum = 0.0f; for (int s = 0; s < segments; ++s) frac_sum += preset_segment_fraction[s];
            float frac_right = (segments - 1 < MAX_DIVIDER_CAPACITY + 1) ? preset_segment_fraction[segments - 1] : 0.0f;
            if (frac_sum <= 0.0f) frac_sum = 1.0f;
            N_right = (int)llround((frac_right / frac_sum) * Nreq);
        } else if (cli_left_empty && segments == 3) {
            // If left-empty and no explicit fractions, split equally over 2 and 3
            N_right = Nreq / 2;
        } else {
            // Fallback: assume equal share among segments
            N_right = Nreq / segments;
        }

        if (N_right < 1) N_right = 1;
        float kBT = kBT_effective();
        float target = cli_calib_target;
        if (target <= 0.0f) target = 1.0f; // guard

        // w_needed = N_right * kT * Δx / W_target
        float w_needed_sigma = (N_right * kBT * stroke_sigma) / target;
        // Clamp resulting wall position within [min,max]
        float pos_sigma = (2.0f * L0_UNITS) - w_needed_sigma;
        if (cli_calib_min_sigma < cli_calib_max_sigma) {
            if (pos_sigma < cli_calib_min_sigma) pos_sigma = cli_calib_min_sigma;
            if (pos_sigma > cli_calib_max_sigma) pos_sigma = cli_calib_max_sigma;
        }
        // Program as a fraction for the last wall
        float frac = pos_sigma / fmaxf(1e-6f, 2.0f * L0_UNITS);
        if (frac < 0.0f) frac = 0.0f; if (frac > 1.0f) frac = 1.0f;
        preset_custom_positions = true;
        preset_wall_count = cli_requested_walls;
        preset_wall_fraction[cli_requested_walls - 1] = frac;
        // Note: we leave earlier walls as user-defined or defaults
        printf("[CAL] Approx calibrated last wall to %.3f σ (frac=%.4f) for target W=%.3g with N_right=%d\n",
               pos_sigma, frac, (double)target, N_right);
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
    
    // Apply multi-wall mass factors (use first wall mass if multiple walls)
    if (cli_num_walls > 0 && cli_wall_masses[0] > 0.0f) {
        wall_mass_runtime = PARTICLE_MASS * cli_wall_masses[0];
    }
    if (cli_override_wall_thickness_sigma > 0.0f) {
        set_wall_thickness_sigma(cli_override_wall_thickness_sigma);
    }
    if (cli_override_temperature > 0.0f) {
        temperature_runtime = cli_override_temperature;
    }
    // Apply particle size overrides (sigma units)
    if (cli_particle_radius_set) {
        // PARTICLE_RADIUS = PARTICLE_RADIUS_UNIT * PIXELS_PER_SIGMA * particle_scale_runtime
        // We want: PARTICLE_RADIUS = cli_particle_radius_sigma * PIXELS_PER_SIGMA
        // => particle_scale_runtime = cli_particle_radius_sigma / PARTICLE_RADIUS_UNIT
        particle_scale_runtime = cli_particle_radius_sigma / PARTICLE_RADIUS_UNIT;
    } else if (cli_particle_diameter_set) {
        // diameter provided → radius = diameter/2
        float r_sigma = 0.5f * cli_particle_diameter_sigma;
        particle_scale_runtime = r_sigma / PARTICLE_RADIUS_UNIT;
    }
    if (cli_force_no_experiments) {
        enable_speed_of_sound_experiments = false;
    }
    
    // Enable energy measurement if requested
    if (cli_enable_energy_measurement) {
        energy_measurement.spring_enabled = true;
        printf("🔬 Energy measurement enabled\n");
    }
    if (cli_experiment_repeats < 1) {
        cli_experiment_repeats = 1;
    }

    if (!preset_custom_fractions) {
        memset(preset_segment_fraction, 0, sizeof(preset_segment_fraction));
    }

    // Map CLI particle-boxes overrides to segment fractions
    {
        int segs = (cli_requested_walls > 0 ? cli_requested_walls + 1 : 1);
        if (cli_particles_box_left >= 0 || cli_particles_box_right >= 0) {
            if (segs >= 2) {
                double left = cli_particles_box_left >= 0 ? (double)cli_particles_box_left : 0.0;
                double right = cli_particles_box_right >= 0 ? (double)cli_particles_box_right : 0.0;
                if (left <= 0.0 && right <= 0.0) { left = right = 1.0; }
                double sum = left + right; if (sum <= 0.0) sum = 1.0;
                preset_segment_fraction[0] = (float)(left / sum);
                preset_segment_fraction[1] = (float)(right / sum);
                preset_custom_fractions = true;
            }
        } else if (cli_particles_box_fractions_count > 0) {
            double sum = 0.0; for (size_t k = 0; k < cli_particles_box_fractions_count; ++k) sum += cli_particles_box_fractions[k];
            if (sum <= 0.0) sum = 1.0;
            for (int s = 0; s < segs; ++s) {
                double v = 0.0; if (s < (int)cli_particles_box_fractions_count) v = cli_particles_box_fractions[s] / sum;
                preset_segment_fraction[s] = (float)v;
            }
            double acc = 0.0; for (int s = 0; s < segs; ++s) acc += preset_segment_fraction[s];
            if (acc < 1.0 && segs > 0) preset_segment_fraction[segs - 1] += (float)(1.0 - acc);
            preset_custom_fractions = true;
        }
    }
}




///////////  SIMULALTION INFO ////////////
void print_simulation_info() {
    printf("\n========== Simulation Info ==========\n");
    printf("Mode: %s\n", MOLECULAR_MODE ? "MOLECULAR (Normalized units)" : "MACROSCOPIC (Real units)");
    printf("Experiment Mode: %s\n", EXPERIMENT_MODE ? "ON (fine dt/substeps)" : "OFF (default dt)");
    printf("Fixed Time Step (dt): %.2e\n", FIXED_DT);
    printf("Runtime Time Step (dt_runtime): %.2e\n", fixed_dt_runtime);
    printf("TAU scaling factor: %.4f\n", TAU_time_scale_factor_to_molecular);
    // Dimensionless tau for frequency normalization (tau = sqrt(m*sigma^2/(kB*T))) in reduced units
    const double sigma_unit = (double)DIAMETER; // in reduced units, DIAMETER = 1 when PARTICLE_RADIUS_UNIT=0.5
    const double tau_dimless = sqrt((double)PARTICLE_MASS * sigma_unit * sigma_unit / ((double)kBT_effective()));
    printf("Tau (dimensionless, for f* = f·tau): %.6f\n", tau_dimless);
    // Runtime thermodynamic parameters
    printf("Temperature (runtime): %.3f\n", (double)temperature_runtime);
    printf("k_B (effective): %.6f\n", (double)kB_effective());
    printf("k_B*T (effective): %.6f\n", (double)kBT_effective());
    printf("k_B*T=1 mode: %s\n", cli_force_kbt_one ? "ON" : "OFF");
    printf("Time scale (interactive): %.2f×\n", (double)time_scale_runtime);
    // Effective length (half-length minus one diameter)
    const double L0_sigma = (double)L0_UNITS;
    const double L_eff = L0_sigma - 1.0; // σ = 1 in reduced units
    printf("L0 (half-length): %.3f σ,  L_eff = L0 − 1 = %.3f σ\n", L0_sigma, L_eff);
    // Particle size summary
    const double r_sigma = (double)(PARTICLE_RADIUS) / (double)PIXELS_PER_SIGMA;
    const double d_sigma = 2.0 * r_sigma;
    printf("Particle radius: %.3f σ (diameter: %.3f σ)\n", r_sigma, d_sigma);
    printf("Wall hold mode: %s\n", wall_hold_enabled ? "Enabled" : "Disabled");
    printf("Wall hold steps: %d\n", wall_hold_steps);
    printf("Internal walls: %d\n", num_internal_walls);
    printf("Target total simulation time: %.1f\n", 3000.0f);
    printf("Calculated steps to reach target: %d\n", (int)(3000.0f / (fixed_dt_runtime / (MOLECULAR_MODE ? 1.0f : TAU_time_scale_factor_to_molecular))));
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

    for (int i = 0; i < particles_active; i++) {
        vx_sum += Vx[i];
        vy_sum += Vy[i];
    }

    float vx_mean = vx_sum / (particles_active > 0 ? particles_active : 1);
    float vy_mean = vy_sum / (particles_active > 0 ? particles_active : 1);


    for (int i = 0; i < particles_active; i++) {
        Vx[i] -= vx_mean;
        Vy[i] -= vy_mean;
    }

    printf("✅ Drift removed: mean Vx = %.4f, mean Vy = %.4f\n", vx_mean, vy_mean);
}







///////////////////////////   ENERGY CALCULATIONS   ///////////////////////
double kinetic_energy() {
    double ke_total = 0.0;
    int n = particles_active > 0 ? particles_active : 0;
    for (int i = 0; i < n; i++) {
        double particle_ke = 0.5 * PARTICLE_MASS * (Vx[i] * Vx[i] + Vy[i] * Vy[i]);  // In Joules
        ke_total += particle_ke;
    }
    return ke_total;
}

float compute_measured_temperature_from_ke() {
    // USES RELATION:
    // T = E_kin / (N k_B) in 2D (with our conventions)
    double ke_total = kinetic_energy();
    int n = particles_active > 0 ? particles_active : 1;
    return (float)(ke_total / (n * kB_effective()));
}



// Returns a 2D Maxwell-Boltzmann distributed velocity vector (vx, vy)
void maxwell_boltzmann_2D(float temperature, float *vx, float *vy) {
    float sigma = sqrt(kB_effective() * temperature / PARTICLE_MASS);
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
    float sigma = sqrt(kB_effective() * temperature / PARTICLE_MASS);  // Standard deviation (velocity scale)
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
    float sigma = sqrt(kB_effective() * temperature / PARTICLE_MASS);
    float u1, u2, z;

    do {
        u1 = (float)rand() / RAND_MAX;
        u2 = (float)rand() / RAND_MAX;
        z = sqrtf(-2.0f * logf(u1)) * cosf(2.0f * M_PI * u2);
    } while (isnan(z));

    return sigma * z;
}




void maxwell_boltzmann_velocity_ARRAY_ALL(float temperature) {
    float sigma = sqrt(kB_effective() * temperature / PARTICLE_MASS);

    int n = particles_active > 0 ? particles_active : 0;
    for (int i = 0; i < n; i += 2) {
        float u1 = (float)rand() / RAND_MAX;
        float u2 = (float)rand() / RAND_MAX;

        float z0 = sqrt(-2.0f * log(u1)) * cos(2.0f * M_PI * u2);
        float z1 = sqrt(-2.0f * log(u1)) * sin(2.0f * M_PI * u2);

        Vx[i] = sigma * z0;
        Vy[i] = sigma * z1;

        // If odd number of particles:
        if (i + 1 < n) {
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
    float sigma = sqrt(kB_effective() * temperature / PARTICLE_MASS);  // Compute sigma

    float u1, u2, u3, z1, z2, z3;

    // Generate three independent standard normal variables using Box-Muller
    do {
        u1 = (float)rand() / RAND_MAX;
        u2 = (float)rand() / RAND_MAX;
        u3 = (float)rand() / RAND_MAX;

        z1 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
        z2 = sqrt(-2.0 * log(u1)) * sin(2.0 * M_PI * u2);
        z3 = sqrt(-2.0 * log(u3)) * cos(2.0 * M_PI * u3);
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

    return v;  // Return the speed
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

    int Nreq = (cli_override_particles > 0) ? cli_override_particles : NUM_PARTICLES;
    if (Nreq < 0) Nreq = 0; if (Nreq > NUM_PARTICLES) Nreq = NUM_PARTICLES;
    // compute per-segment area caps (approx)
    float height_eff = (YW2 - YW1) - 2.0f * PARTICLE_RADIUS;
    float area_p = (float)(M_PI * PARTICLE_RADIUS * PARTICLE_RADIUS);
    int caps[MAX_DIVIDER_CAPACITY + 1]; int cap_sum = 0;
    for (int s = 0; s < segments; ++s) {
        float w = fmaxf(0.0f, (seg_right[s] - seg_left[s]) - 2.0f * PARTICLE_RADIUS);
        int cap = (int)floorf(cli_phi_max * (w * height_eff) / area_p);
        if (cap < 0) cap = 0; caps[s] = cap; cap_sum += cap;
    }
    if (Nreq > cap_sum) {
        fprintf(stderr, "❌ Overpacked: requested N=%d exceeds capacity=%d at phi_max=%.3f.\n",
                Nreq, cap_sum, (double)cli_phi_max);
        fprintf(stderr, "   Tip: reduce --particles or increase box size (L0/height) or decrease wall count.\n");
        exit(1);
    }
    particles_active = Nreq;

    int base = (segments > 0) ? (Nreq / segments) : 0;
    int remainder = (segments > 0) ? (Nreq % segments) : 0;
    int target_counts[MAX_DIVIDER_CAPACITY + 1];
    for (int i = 0; i < segments; ++i) target_counts[i] = 0;

    // Absolute two-box override if user provided explicit counts
    bool applied_absolute_two_box = false;
    if (num_internal_walls == 1 && (cli_particles_box_left >= 0 || cli_particles_box_right >= 0)) {
        int left = cli_particles_box_left >= 0 ? cli_particles_box_left : (Nreq - cli_particles_box_right);
        if (left < 0) left = 0; if (left > Nreq) left = Nreq;
        int right = Nreq - left;
        if (left > caps[0]) { right += (left - caps[0]); left = caps[0]; }
        if (right > caps[1]) { left += (right - caps[1]); right = caps[1]; if (left > caps[0]) left = caps[0]; }
        target_counts[0] = left;
        target_counts[1] = right;
        applied_absolute_two_box = true;
    }

    if (!applied_absolute_two_box && preset_custom_fractions && preset_wall_count + 1 == segments) {
        double accum = 0.0;
        int assigned = 0;
        for (int seg = 0; seg < segments - 1; ++seg) {
            double frac = preset_segment_fraction[seg];
            if (frac < 0.0) frac = 0.0;
            accum += frac;
            int c = (int)llround(frac * Nreq);
            if (c < 0) c = 0;
            if (c > caps[seg]) c = caps[seg];
            target_counts[seg] = c;
            assigned += c;
        }
        int last = Nreq - assigned;
        if (last < 0) last = 0;
        target_counts[segments - 1] = (last > caps[segments - 1]) ? caps[segments - 1] : last;
        for (int seg = 0; seg < segments; ++seg) {
            if (target_counts[seg] < 0) target_counts[seg] = 0;
        }
    } else if (!applied_absolute_two_box) {
        int sum = 0;
        for (int seg = 0; seg < segments; ++seg) {
            int c = base + (seg < remainder ? 1 : 0);
            if (c > caps[seg]) c = caps[seg];
            target_counts[seg] = c; sum += c;
        }
        int leftover = Nreq - sum;
        for (int seg = 0; seg < segments && leftover > 0; ++seg) {
            int room = caps[seg] - target_counts[seg];
            int add = room < leftover ? room : leftover;
            target_counts[seg] += add; leftover -= add;
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

        for (int c = 0; c < target && index < particles_active; ++c, ++index) {
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

            Vx[index] = maxwell_boltzmann_velocity_gaussians(temperature_runtime);
            Vy[index] = maxwell_boltzmann_velocity_gaussians(temperature_runtime);
            Radius[index] = PARTICLE_RADIUS;
        }
    }

    // ensure we seeded exactly particles_active; target_counts logic should guarantee this

    if (all_wall_positions) {
        all_wall_positions[primary_wall_index] = wall_x;
    }

    // Populate initial segment_counts so HUD shows correct numbers before SPACE
    if (segment_counts && segment_count > 0) {
        memset(segment_counts, 0, (size_t)segment_count * sizeof(int));
        for (int i = 0; i < particles_active; ++i) {
            int seg = segment_index_for_position(X[i]);
            if (seg < 0) seg = 0; if (seg >= segment_count) seg = segment_count - 1;
            segment_counts[seg]++;
        }
    }
    // Equalize temperature across segments at initialization
    equalize_temperature_per_segment(temperature_runtime);
}




// Deterministic grid seeding for arbitrary number of segments (num_internal_walls+1)
// Particles are placed on a near-square rectangular lattice per segment, centered with a 1R margin,
// so the piston interacts with a consistent slice each run.
void initialize_simulation_segmented_grid() {
    int segments = num_internal_walls + 1;
    if (segments <= 0) segments = 1;

    float seg_left[MAX_DIVIDER_CAPACITY + 1];
    float seg_right[MAX_DIVIDER_CAPACITY + 1];
    if (segments > MAX_DIVIDER_CAPACITY + 1) segments = MAX_DIVIDER_CAPACITY + 1;
    if (num_internal_walls > 0) {
        compute_segment_bounds(seg_left, seg_right);
    } else {
        seg_left[0] = (float)XW1; seg_right[0] = (float)XW2;
    }

    // Capacity caps and requested N
    int Nreq = (cli_override_particles > 0) ? cli_override_particles : NUM_PARTICLES;
    if (Nreq < 0) Nreq = 0; if (Nreq > NUM_PARTICLES) Nreq = NUM_PARTICLES;
    float height_eff = (YW2 - YW1) - 2.0f * PARTICLE_RADIUS;
    float area_p = (float)(M_PI * PARTICLE_RADIUS * PARTICLE_RADIUS);
    int caps[MAX_DIVIDER_CAPACITY + 1]; int cap_sum = 0;
    for (int s = 0; s < segments; ++s) {
        float w = fmaxf(0.0f, (seg_right[s] - seg_left[s]) - 2.0f * PARTICLE_RADIUS);
        int cap = (int)floorf(cli_phi_max * (w * height_eff) / area_p);
        if (cap < 0) cap = 0; caps[s] = cap; cap_sum += cap;
    }
    if (Nreq > cap_sum) {
        fprintf(stderr, "❌ Overpacked: requested N=%d exceeds capacity=%d at phi_max=%.3f.\n",
                Nreq, cap_sum, (double)cli_phi_max);
        fprintf(stderr, "   Tip: reduce --particles or increase box size (L0/height) or decrease wall count.\n");
        exit(1);
    }
    particles_active = Nreq;

    // Distribute target counts per segment
    int target_counts[MAX_DIVIDER_CAPACITY + 1];
    for (int i = 0; i < segments; ++i) target_counts[i] = 0;

    bool applied_absolute_two_box = false;
    bool applied_absolute_multi   = false;
    if (num_internal_walls == 1 && (cli_particles_box_left >= 0 || cli_particles_box_right >= 0)) {
        int left = cli_particles_box_left >= 0 ? cli_particles_box_left : (Nreq - cli_particles_box_right);
        if (left < 0) left = 0; if (left > Nreq) left = Nreq;
        int right = Nreq - left;
        if (left > caps[0]) { right += (left - caps[0]); left = caps[0]; }
        if (right > caps[1]) { left += (right - caps[1]); right = caps[1]; if (left > caps[0]) left = caps[0]; }
        target_counts[0] = left; target_counts[1] = right;
        applied_absolute_two_box = true;
    }
    // Absolute counts override for multi-wall
    if (!applied_absolute_two_box && preset_custom_absolute_counts && (int)cli_particles_box_counts_count == segments) {
        int sum = 0;
        for (int s = 0; s < segments; ++s) {
            int c = cli_particles_box_counts[s]; if (c < 0) c = 0;
            if (c > caps[s]) {
                fprintf(stderr, "❌ Overpacked in segment %d: requested %d > cap %d at phi_max=%.3f.\n", s, c, caps[s], (double)cli_phi_max);
                exit(1);
            }
            target_counts[s] = c; sum += c;
        }
        particles_active = sum;
        applied_absolute_multi = true;
    }

    if (!applied_absolute_two_box && !applied_absolute_multi && preset_custom_fractions && preset_wall_count + 1 == segments) {
        int assigned = 0;
        for (int s = 0; s < segments - 1; ++s) {
            double frac = preset_segment_fraction[s];
            if (frac < 0.0) frac = 0.0;
            int c = (int)llround(frac * Nreq);
            if (c < 0) c = 0; if (c > caps[s]) c = caps[s];
            target_counts[s] = c; assigned += c;
        }
        int last = Nreq - assigned; if (last < 0) last = 0;
        target_counts[segments - 1] = (last > caps[segments - 1]) ? caps[segments - 1] : last;
    } else if (!applied_absolute_two_box && !applied_absolute_multi && cli_distribute_area) {
        // area-proportional distribution across segments
        float weights[MAX_DIVIDER_CAPACITY + 1]; float wsum = 0.0f;
        for (int s = 0; s < segments; ++s) {
            float w = fmaxf(0.0f, (seg_right[s] - seg_left[s]) - 2.0f * PARTICLE_RADIUS);
            float a = w * height_eff; if (a < 0.0f) a = 0.0f; weights[s] = a; wsum += a;
        }
        if (wsum <= 0.0f) wsum = 1.0f;
        int sum = 0;
        for (int s = 0; s < segments; ++s) {
            int c = (int)llround((weights[s] / wsum) * Nreq); if (c > caps[s]) c = caps[s]; if (c < 0) c = 0;
            target_counts[s] = c; sum += c;
        }
        int leftover = Nreq - sum;
        for (int s = 0; s < segments && leftover > 0; ++s) {
            int room = caps[s] - target_counts[s]; int add = room < leftover ? room : leftover; target_counts[s] += add; leftover -= add;
        }
    } else if (!applied_absolute_two_box && !applied_absolute_multi) {
        int base = (segments > 0) ? (Nreq / segments) : 0;
        int remainder = (segments > 0) ? (Nreq % segments) : 0;
        int sum = 0;
        for (int s = 0; s < segments; ++s) {
            int c = base + (s < remainder ? 1 : 0);
            if (c > caps[s]) c = caps[s];
            target_counts[s] = c; sum += c;
        }
        int leftover = Nreq - sum;
        for (int s = 0; s < segments && leftover > 0; ++s) {
            int room = caps[s] - target_counts[s];
            int add = room < leftover ? room : leftover;
            target_counts[s] += add; leftover -= add;
        }
    }

    // If left-empty is requested, force segment 0 to zero and re-distribute
    if (!applied_absolute_multi && cli_left_empty) {
        // Start with equal shares across segments 1..end within caps
        int sum = 0;
        for (int s = 0; s < segments; ++s) { sum += target_counts[s]; }
        int leftover = Nreq - (sum - target_counts[0]);
        target_counts[0] = 0;
        for (int s = 1; s < segments && leftover > 0; ++s) {
            int room = caps[s] - target_counts[s];
            int add = room < leftover ? room : leftover;
            target_counts[s] += add; leftover -= add;
        }
        // If still leftover, we were over capacity excluding left box
        if (leftover > 0) {
            fprintf(stderr, "❌ Overpacked (left-empty): requested N exceeds capacity of non-left segments.\n");
            exit(1);
        }
    }

    // Place on grid per segment
    int index = 0;
    const float margin = PARTICLE_RADIUS; // 1R margin from each side
    for (int s = 0; s < segments; ++s) {
        int target = target_counts[s];
        if (target <= 0) continue;
        float left_bound  = seg_left[s]  + margin;
        float right_bound = seg_right[s] - margin;
        float top = YW1 + margin;
        float bottom = YW2 - margin;
        float w = fmaxf(0.0f, right_bound - left_bound);
        float h = fmaxf(0.0f, bottom - top);
        if (w <= 0.0f || h <= 0.0f) continue;

        // Choose rows/cols to be near-square, then enforce spacing >= 2R if possible
        int rows = (int)floorf(sqrtf((float)target));
        if (rows < 1) rows = 1;
        int cols = (target + rows - 1) / rows;
        float dx = (cols > 0) ? (w / cols) : w;
        float dy = (rows > 0) ? (h / rows) : h;
        float min_gap = 2.0f * PARTICLE_RADIUS;
        // Increase rows until dy >= 2R (or rows == target)
        while (rows < target && dy < min_gap) {
            rows++;
            cols = (target + rows - 1) / rows;
            dx = (cols > 0) ? (w / cols) : w;
            dy = (rows > 0) ? (h / rows) : h;
        }
        // If still too tight horizontally, increase cols similarly by reducing rows minimally
        while (cols < target && dx < min_gap) {
            cols++;
            rows = (target + cols - 1) / cols;
            dx = (cols > 0) ? (w / cols) : w;
            dy = (rows > 0) ? (h / rows) : h;
        }

        int placed = 0;
        for (int r = 0; r < rows && placed < target; ++r) {
            for (int c = 0; c < cols && placed < target; ++c) {
                if (index >= particles_active) break;
                float cx = left_bound + (c + 0.5f) * dx;
                float cy = top        + (r + 0.5f) * dy;
                X[index] = cx;
                Y[index] = cy;
                Vx[index] = maxwell_boltzmann_velocity_gaussians(temperature_runtime);
                Vy[index] = maxwell_boltzmann_velocity_gaussians(temperature_runtime);
                Radius[index] = PARTICLE_RADIUS;
                index++; placed++;
            }
        }
    }

    // Velocity rescale to match target temperature exactly
    double actual_ke = kinetic_energy();
    double target_ke = particles_active * K_B * temperature_runtime;
    if (actual_ke > 0.0) {
        double scale = sqrt(target_ke / actual_ke);
        for (int i = 0; i < particles_active; ++i) {
            Vx[i] *= scale;
            Vy[i] *= scale;
        }
    }

    // Initialize segment_counts for HUD
    if (segment_counts && segment_count > 0) {
        memset(segment_counts, 0, (size_t)segment_count * sizeof(int));
        for (int i = 0; i < particles_active; ++i) {
            int seg = segment_index_for_position(X[i]);
            if (seg < 0) seg = 0; if (seg >= segment_count) seg = segment_count - 1;
            segment_counts[seg]++;
        }
    }
    // Equalize temperature across segments at initialization
    equalize_temperature_per_segment(temperature_runtime);
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

    // compute per-segment caps (approx area-based)
    float height_eff = (YW2 - YW1) - 2.0f * PARTICLE_RADIUS;
    float left_eff = fmaxf(0.0f, left_width - 2.0f * PARTICLE_RADIUS);
    float right_eff = fmaxf(0.0f, right_width - 2.0f * PARTICLE_RADIUS);
    float area_p = (float)(M_PI * PARTICLE_RADIUS * PARTICLE_RADIUS);
    int cap_left = (int)floorf(cli_phi_max * (left_eff * height_eff) / area_p);
    int cap_right = (int)floorf(cli_phi_max * (right_eff * height_eff) / area_p);

    int Nreq = (cli_override_particles > 0) ? cli_override_particles : NUM_PARTICLES;
    Nreq = (int)fminf((float)Nreq, (float)(cap_left + cap_right));
    if (Nreq < 0) Nreq = 0;
    particles_active = Nreq;

    int particles_left = Nreq / 2;
    if (cli_particles_box_left >= 0 || cli_particles_box_right >= 0) {
        int left = cli_particles_box_left >= 0 ? cli_particles_box_left : (Nreq - cli_particles_box_right);
        if (left < 0) left = 0; if (left > Nreq) left = Nreq;
        particles_left = left;
    } else if (preset_custom_fractions) {
        int desired_left = (int)llround(preset_segment_fraction[0] * Nreq);
        if (desired_left >= 0 && desired_left <= NUM_PARTICLES) particles_left = desired_left;
    }
    // clamp by caps and re-balance
    if (particles_left > cap_left) particles_left = cap_left;
    int particles_right = Nreq - particles_left;
    if (particles_right > cap_right) {
        particles_right = cap_right;
        particles_left = Nreq - particles_right;
        if (particles_left > cap_left) particles_left = cap_left;
        particles_active = particles_left + particles_right;
    }

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
            Vx[idx] = maxwell_boltzmann_velocity_gaussians(temperature_runtime);
            Vy[idx] = maxwell_boltzmann_velocity_gaussians(temperature_runtime);
            Radius[idx] = PARTICLE_RADIUS;

            printf("✅ Particle %d (LEFT) placed at (%.2f, %.2f)\n", idx, x, y);
            idx++;
        }
    }

    // --- Place Right Particles ---
    for (int row = 0; row < num_rows_right; row++) {
        for (int col = 0; col < num_cols_right; col++) {
            if (idx >= particles_left + particles_right) break;

            float x = wall_x + wall_buffer + (col + 0.5f) * horizontal_spacing_right;
            float y = YW1 + (row + 0.5f) * vertical_spacing_right;

            X[idx] = x;
            Y[idx] = y;
            Vx[idx] = maxwell_boltzmann_velocity_gaussians(temperature_runtime);
            Vy[idx] = maxwell_boltzmann_velocity_gaussians(temperature_runtime);
            Radius[idx] = PARTICLE_RADIUS;

            printf("✅ Particle %d (RIGHT) placed at (%.2f, %.2f)\n", idx, x, y);
            idx++;
        }
    }

    printf("✅ Finished placing all %d particles evenly across box!\n", particles_active);
    // Diagnostics: capacity and targets
    printf("[SEED] MAX(buffer)=%d requested=%d active=%d phi_max=%.3f\n",
           NUM_PARTICLES,
           (cli_override_particles > 0 ? cli_override_particles : NUM_PARTICLES),
           particles_active, cli_phi_max);
    printf("[SEED] caps: L=%d R=%d  targets: L=%d R=%d\n", cap_left, cap_right, particles_left, particles_right);

    // Populate initial segment_counts so HUD shows correct numbers before SPACE
    if (segment_counts && segment_count > 0) {
        memset(segment_counts, 0, (size_t)segment_count * sizeof(int));
        for (int i = 0; i < particles_active; ++i) {
            int seg = segment_index_for_position(X[i]);
            if (seg < 0) seg = 0; if (seg >= segment_count) seg = segment_count - 1;
            segment_counts[seg]++;
        }
    }
    printf("🔍 Initial particle velocities (first 10):\n");
    for (int i = 0; i < 10; i++) {
        printf("Particle %2d: Vx = %.6f, Vy = %.6f\n", i, Vx[i], Vy[i]);    
    }
    // Equalize temperature across segments at initialization
    equalize_temperature_per_segment(temperature_runtime);
}





void initialize_simulation(void) {
    // Choose seeding by CLI, with sensible defaults:
    // - grid for reproducibility and consistent piston contact
    // - honeycomb as alternative lattice
    // - random for stochastic tests
    switch (cli_seeding_mode) {
        case SEEDING_GRID:
            if (num_internal_walls > 1) {
                initialize_simulation_segmented_grid();
            } else {
                // existing two-box grid initializer
                ; // fall through to single-wall grid logic below
            }
            break;
        case SEEDING_HONEYCOMB:
            initialize_simulation_honeycomb();
            return;
        case SEEDING_RANDOM:
            initialize_simulation_random();
            return;
    }
    if (num_internal_walls > 1) return; // already handled segmented grid
    float wall_buffer = WALL_THICKNESS + DIAMETER * 0.25f; // space from wall
    float margin = DIAMETER * 1.25f;

    float left_width = wall_x - wall_buffer - XW1;
    float right_width = XW2 - wall_x - wall_buffer;
    float height = YW2 - YW1;

    // Compute per-box capacity (area-based)
    float height_eff = (YW2 - YW1) - 2.0f * PARTICLE_RADIUS;
    float left_eff = fmaxf(0.0f, left_width - 2.0f * PARTICLE_RADIUS);
    float right_eff = fmaxf(0.0f, right_width - 2.0f * PARTICLE_RADIUS);
    float area_p = (float)(M_PI * PARTICLE_RADIUS * PARTICLE_RADIUS);
    int cap_left = (int)floorf(cli_phi_max * (left_eff * height_eff) / area_p);
    int cap_right = (int)floorf(cli_phi_max * (right_eff * height_eff) / area_p);

    int requested_total = (cli_override_particles > 0) ? cli_override_particles : NUM_PARTICLES;
    int req_left = (cli_particles_box_left >= 0) ? cli_particles_box_left : -1;
    int req_right = (cli_particles_box_right >= 0) ? cli_particles_box_right : -1;

    int particles_left = 0, particles_right = 0;
    // Absolute counts for two segments via --particles-boxes counts (length==2)
    if (preset_custom_absolute_counts && cli_particles_box_counts_count == 2) {
        int left = cli_particles_box_counts[0]; if (left < 0) left = 0; if (left > cap_left) {
            printf("❌ Overpacked (left): requested %d > cap %d.\n", left, cap_left); exit(1);
        }
        int right = cli_particles_box_counts[1]; if (right < 0) right = 0; if (right > cap_right) {
            printf("❌ Overpacked (right): requested %d > cap %d.\n", right, cap_right); exit(1);
        }
        particles_left = left; particles_right = right; particles_active = left + right;
    } else
    if (req_left >= 0 && req_right >= 0) {
        particles_active = req_left + req_right;
        if (particles_active > cap_left + cap_right || req_left > cap_left || req_right > cap_right) {
            printf("❌ Overpacked: requested total=%d (L=%d,R=%d) > capacity total=%d (L_cap=%d,R_cap=%d).\n",
                   particles_active, req_left, req_right, cap_left + cap_right, cap_left, cap_right);
            exit(1);
        }
        particles_left = req_left;
        particles_right = req_right;
    } else {
        if (requested_total > cap_left + cap_right) {
            printf("❌ Overpacked: requested total=%d > capacity total=%d (L_cap=%d,R_cap=%d).\n",
                   requested_total, cap_left + cap_right, cap_left, cap_right);
            exit(1);
        }
        particles_active = requested_total;
        if (cli_distribute_area && !preset_custom_fractions && !preset_custom_absolute_counts) {
            // area-proportional split across left/right
            float wL = fmaxf(0.0f, left_eff);
            float wR = fmaxf(0.0f, right_eff);
            float sumw = wL + wR; if (sumw <= 0.0f) sumw = 1.0f;
            particles_left = (int)llround((wL / sumw) * particles_active);
            if (particles_left > cap_left) particles_left = cap_left;
        } else if (cli_left_empty) {
            particles_left = 0;
            particles_right = particles_active;
            if (particles_right > cap_right) {
                printf("❌ Overpacked (left-empty): requested %d exceeds right capacity %d.\n", particles_right, cap_right);
                exit(1);
            }
        } else if (preset_custom_fractions) {
            particles_left = (int)llround(preset_segment_fraction[0] * particles_active);
            if (particles_left < 0) particles_left = 0; if (particles_left > particles_active) particles_left = particles_active;
        } else {
            particles_left = particles_active / 2;
        }
        particles_right = particles_active - particles_left;
        if (particles_left > cap_left || particles_right > cap_right) {
            printf("❌ Per-box overpack: targets L=%d,R=%d exceed caps L_cap=%d,R_cap=%d.\n",
                   particles_left, particles_right, cap_left, cap_right);
            exit(1);
        }
    }

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
            maxwell_boltzmann_2D(temperature_runtime, &vx, &vy);
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
            if (idx >= particles_left + particles_right) break;
            X[idx] = XW2 - margin - col * spacing_right_x;
            Y[idx] = YW1 + margin + row * spacing_right_y;
            float vx, vy;
            maxwell_boltzmann_2D(temperature_runtime, &vx, &vy);
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
    double target_ke = particles_active * K_B * temperature_runtime;
    double scale = sqrt(target_ke / actual_ke);  // to adjust v_i → v_i * scale
    printf("✅ Scaling factor for velocities: %.4f\n", scale);

    for (int i = 0; i < particles_active; i++) {
        Vx[i] *= scale;
        Vy[i] *= scale;
}

    // Diagnostics and HUD seed counts
    printf("[SEED] MAX(buffer)=%d requested=%d active=%d phi_max=%.3f\n",
           NUM_PARTICLES, requested_total, particles_active, cli_phi_max);
    printf("[SEED] caps: L=%d R=%d  targets: L=%d R=%d\n", cap_left, cap_right, particles_left, particles_right);

    if (segment_counts && segment_count > 0) {
        memset(segment_counts, 0, (size_t)segment_count * sizeof(int));
        for (int i = 0; i < particles_active; ++i) {
            int seg = segment_index_for_position(X[i]);
            if (seg < 0) seg = 0; if (seg >= segment_count) seg = segment_count - 1;
            segment_counts[seg]++;
        }
    }
    // Equalize initial temperature per segment to a global target, then apply
    // per-segment overrides if provided.
    equalize_temperature_per_segment(temperature_runtime);
    if (cli_temperature_segments && cli_temperature_segments_count > 0) {
        apply_segment_temperatures(cli_temperature_segments, cli_temperature_segments_count);
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

// Simple spring visualization between left boundary and leftmost internal wall
static void render_left_spring() {
    if (!energy_measurement.spring_enabled || num_internal_walls <= 0) return;
    int x1 = XW1 + 4;
    int x2 = (int)(all_wall_positions[0] - WALL_THICKNESS * 0.5f);
    if (x2 <= x1 + 4) return;
    int y = YW1 + SIM_HEIGHT / 2;
    const int coils = 14;
    const int amp = (int)fminf(20.0f, SIM_HEIGHT / 6.0f);
    SDL_SetRenderDrawColor(_renderer, 200, 200, 255, 220);
    // draw as polyline zig-zag
    int segs = coils * 2;
    for (int k = 0; k < segs; ++k) {
        float t0 = (float)k / segs;
        float t1 = (float)(k + 1) / segs;
        int px0 = x1 + (int)((x2 - x1) * t0);
        int px1 = x1 + (int)((x2 - x1) * t1);
        int py0 = y + ((k % 2) ? amp : -amp);
        int py1 = y + (((k + 1) % 2) ? amp : -amp);
        SDL_RenderDrawLine(_renderer, px0, py0, px1, py1);
    }
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

// Sigmoid function for smooth piston acceleration
float sigmoid_velocity(float t, float max_velocity, float duration) {
    if (t <= 0) return 0.0f;
    if (t >= duration) return max_velocity;
    
    // Sigmoid: v(t) = v_max * (1 / (1 + exp(-k*(t - duration/2))))
    float k = 8.0f / duration;  // Steepness parameter
    float sigmoid = 1.0f / (1.0f + expf(-k * (t - duration / 2.0f)));
    return max_velocity * sigmoid;
}

// Protocol-based piston movement
void apply_piston_protocol(float dt) {
    static float protocol_start_time = -1.0f;
    static float protocol_duration = 2.0f;  // 2 seconds for protocol
    static float max_protocol_velocity = 500.0f;
    
    if (piston_step_active) {
        if (protocol_start_time < 0) {
            protocol_start_time = simulation_time;
        }
        
        float elapsed = simulation_time - protocol_start_time;
        
        switch (cli_protocol) {
            case PROTOCOL_STEP:
                // Instant acceleration (current behavior)
                vx_piston_right = -max_protocol_velocity;
                break;
                
            case PROTOCOL_SIGMOIDAL:
                // Smooth sigmoid acceleration
                vx_piston_right = -sigmoid_velocity(elapsed, max_protocol_velocity, protocol_duration);
                break;
                
            case PROTOCOL_LINEAR:
                // Linear acceleration
                if (elapsed < protocol_duration) {
                    vx_piston_right = -max_protocol_velocity * (elapsed / protocol_duration);
                } else {
                    vx_piston_right = -max_protocol_velocity;
                }
                break;
                
            case PROTOCOL_SINUSOIDAL:
                // Sinusoidal acceleration
                if (elapsed < protocol_duration) {
                    vx_piston_right = -max_protocol_velocity * sinf(M_PI * elapsed / protocol_duration);
                } else {
                    vx_piston_right = -max_protocol_velocity;
                }
                break;
                
            case PROTOCOL_OPTIMAL:
                // AI-optimized protocol (placeholder for future implementation)
                vx_piston_right = -sigmoid_velocity(elapsed, max_protocol_velocity, protocol_duration);
                break;
        }
        
        // Stop when target reached
        if (piston_right_x <= piston_step_target) {
            vx_piston_right = 0;
            piston_step_active = false;
            protocol_start_time = -1.0f;
        }
    }
}

// Energy measurement functions
void initialize_energy_measurement() {
    // Hard gate: only allow spring if explicitly requested, or via energy_transfer preset
    if (!(cli_enable_energy_measurement || cli_experiment_preset == EXPERIMENT_PRESET_ENERGY_TRANSFER)) {
        return;
    }
    if (num_internal_walls > 0) {
        // Default equilibrium: leftmost internal wall position (pixels)
        float eq = 0.0f;
        if (cli_spring_eq_set) {
            // User override in σ units from left boundary
            eq = XW1 + cli_spring_eq_sigma * PIXELS_PER_SIGMA;
        } else {
            // Find leftmost among all current internal wall positions (including primary moved wall)
            float min_pos = wall_x;
            if (all_wall_positions) {
                // Ensure we use the latest positions
                // Note: primary wall is tracked by wall_x; others in extra_wall_positions
                min_pos = wall_x;
                for (int w = 0; w < num_internal_walls; ++w) {
                    float wx = all_wall_positions[w];
                    if (wx < min_pos) min_pos = wx;
                }
            }
            eq = min_pos;
        }
        if (cli_spring_k_set) {
            energy_measurement.spring_constant = cli_spring_k;
        }
        energy_measurement.equilibrium_position = eq;
        energy_measurement.spring_enabled = true;
        energy_measurement.spring_force = 0.0f;
        energy_measurement.spring_energy = 0.0f;
        energy_measurement.spring_energy_max = 0.0f;
        energy_measurement.energy_transferred = 0.0f;
        energy_measurement.max_displacement = 0.0f;
        printf("🔬 Energy measurement initialized: equilibrium at %.2f\n", energy_measurement.equilibrium_position);
    }
}

void update_energy_measurement(float dt) {
    if (!energy_measurement.spring_enabled || num_internal_walls <= 0) return;
    
    // Get current position of leftmost wall
    float current_position = all_wall_positions[0];
    float displacement = current_position - energy_measurement.equilibrium_position;
    
    // Update maximum displacement
    if (fabsf(displacement) > fabsf(energy_measurement.max_displacement)) {
        energy_measurement.max_displacement = displacement;
    }
    
    // Calculate spring force and energy
    energy_measurement.spring_force = -energy_measurement.spring_constant * displacement;
    energy_measurement.spring_energy = 0.5f * energy_measurement.spring_constant * displacement * displacement;
    if (energy_measurement.spring_energy > energy_measurement.spring_energy_max) {
        energy_measurement.spring_energy_max = energy_measurement.spring_energy;
    }
    
    // Calculate energy transferred (work done on the spring)
    static float last_position = 0.0f;
    if (last_position != 0.0f) {
        float work = energy_measurement.spring_force * (current_position - last_position);
        energy_measurement.energy_transferred += work;
    }
    last_position = current_position;
}

void print_energy_measurement() {
    if (!energy_measurement.spring_enabled) return;
    
    printf("🔬 Energy Measurement:\n");
    printf("   Spring Force: %.2f\n", energy_measurement.spring_force);
    printf("   Spring Energy: %.2f\n", energy_measurement.spring_energy);
    printf("   Spring Energy Max: %.2f\n", energy_measurement.spring_energy_max);
    printf("   Energy Transferred: %.2f\n", energy_measurement.energy_transferred);
    printf("   Max Displacement: %.2f\n", energy_measurement.max_displacement);
    printf("   Current Position: %.2f (equilibrium: %.2f)\n", 
           all_wall_positions[0], energy_measurement.equilibrium_position);
}

// Function to update piston positions based on velocity and time step
void update_pistons(float dt) {
    // Apply protocol-based movement if active
    apply_piston_protocol(dt);
    
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

// Helper: blue(0) -> purple(0.5) -> red(1) gradient
static inline void set_draw_color_blue_purple_red(float a01) {
    if (a01 < 0.0f) a01 = 0.0f; if (a01 > 1.0f) a01 = 1.0f;
    Uint8 r = (Uint8)(255.0f * a01);
    Uint8 g = 0;
    Uint8 b = (Uint8)(255.0f * (1.0f - a01));
    SDL_SetRenderDrawColor(_renderer, r, g, b, 255);
}

// Function to render particles with optional color coding
void render_particles() {
    int n = particles_active > 0 ? particles_active : 0;
    if (n <= 0) return;

    // Defaults: white color
    SDL_SetRenderDrawColor(_renderer, 255, 255, 255, 255);

    float vmin = FLT_MAX, vmax = 0.0f;
    if (colorcode_mode == COLORCODE_VELOCITY) {
        for (int i = 0; i < n; ++i) {
            float vx = (float)Vx[i], vy = (float)Vy[i];
            float s = sqrtf(vx*vx + vy*vy);
            if (s > vmax) vmax = s;
            if (s < vmin) vmin = s;
        }
        if (!(vmax > vmin)) { vmin = 0.0f; vmax = 1.0f; }
    }

    // For temperature, compute min/max across segments (non-empty)
    double Tmin = 0.0, Tmax = 1.0; int haveT = 0;
    if (colorcode_mode == COLORCODE_TEMPERATURE && segment_count > 0 && segment_temperature) {
        for (int s = 0; s < segment_count; ++s) {
            if (!segment_counts || segment_counts[s] <= 0) continue;
            double T = segment_temperature[s];
            if (!haveT) { Tmin = Tmax = T; haveT = 1; }
            else { if (T < Tmin) Tmin = T; if (T > Tmax) Tmax = T; }
        }
        if (!haveT || fabs(Tmax - Tmin) < 1e-20) { Tmin = 0.0; Tmax = 1.0; }
    }

    for (int i = 0; i < n; i++) {
        if (colorcode_mode == COLORCODE_NONE) {
            SDL_SetRenderDrawColor(_renderer, 255, 255, 255, 255);
        } else if (colorcode_mode == COLORCODE_VELOCITY) {
            float vx = (float)Vx[i], vy = (float)Vy[i];
            float s = sqrtf(vx*vx + vy*vy);
            float denom = (vmax - vmin);
            float a = denom > 1e-20f ? (s - vmin) / denom : 0.0f;
            set_draw_color_blue_purple_red(a);
        } else if (colorcode_mode == COLORCODE_TEMPERATURE) {
            int seg = segment_index_for_position(X[i]);
            if (seg < 0) seg = 0; if (seg >= segment_count) seg = segment_count - 1;
            double T = (segment_temperature && seg >= 0) ? segment_temperature[seg] : 0.0;
            double a = (T - Tmin) / (Tmax - Tmin);
            if (a < 0.0) a = 0.0; if (a > 1.0) a = 1.0;
            set_draw_color_blue_purple_red((float)a);
        }

        int pixel_x = (int)(X[i]);
        int pixel_y = (int)(Y[i]);
        int radius_pixels = (int)(PARTICLE_RADIUS);
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

    // Apply spring force to the primary (leftmost) internal wall if enabled
    // Treat the spring as connected between the left boundary and the wall with
    // rest position at energy_measurement.equilibrium_position (pixels).
    if (energy_measurement.spring_enabled && num_internal_walls > 0) {
        float disp = wall_x - energy_measurement.equilibrium_position; // pixels
        float Fspring = -energy_measurement.spring_constant * disp;     // arbitrary pixel-units
        float a = Fspring / (float)WALL_MASS;                           // dv/dt
        vx_wall += a * dt;
    }

    wall_x += vx_wall * dt;  // Move wall based on current velocity

    if (steps_elapsed % 100 == 0) {
    printf("📊 Step %d | wall_x = %.6f | vx_wall = %.8f\n", steps_elapsed, wall_x, vx_wall);
    printf("wall_enabled = %d, wall_released = %d, left = %d, right = %d\n",
        wall_enabled, wall_is_released, left_count, right_count);
    }

    // === Enforce Wall Boundaries ===
    // In spring mode, prevent interpenetration without adding non-physical bounce.
    float half_thick = WALL_THICKNESS * 0.5f;
    float left_limit_px  = (float)XW1 + half_thick + 1.0f;
    float right_limit_px = (float)XW2 - half_thick - 1.0f;
    if (energy_measurement.spring_enabled && num_internal_walls > 0) {
        if (wall_x < left_limit_px) {
            wall_x = left_limit_px;
            if (vx_wall < 0.0f) vx_wall = 0.0f;
        }
        if (wall_x > right_limit_px) {
            wall_x = right_limit_px;
            if (vx_wall > 0.0f) vx_wall = 0.0f;
        }
    } else {
        // Legacy elastic bounce when no spring is active
        if (wall_x < left_limit_px) {
            wall_x = left_limit_px;
            vx_wall = -vx_wall;
        }
        if (wall_x > right_limit_px) {
            wall_x = right_limit_px;
            vx_wall = -vx_wall;
        }
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

    if (release_time >= 0.0 && leftmost_wall_release_time < 0.0) {
        leftmost_wall_release_time = release_time;
    }

    if (leftmost_wall_release_time >= 0.0 && !leftmost_wall_peak_recorded) {
        float delta = leftmost_wall_baseline - min_pos;
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
    const float eps_pos = fmaxf(1e-7f, 1e-6f * fmaxf(R, 1.f)); // push-out distance
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

    int n_active = particles_active > 0 ? particles_active : 0;
    for (int step = 0; step < (int)SUBSTEPS; ++step) {
        for (int i = 0; i < n_active; i++) {
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
                    bool hold_active = false;
                    if (extra_wall_hold_enabled) {
                        hold_active = extra_wall_hold_enabled[w];
                        if (extra_wall_is_released) {
                            hold_active = hold_active && !extra_wall_is_released[w];
                        }
                    }
                    ccd_wall_step(i, sub_dt, extra_wall_old_positions[w], &extra_wall_velocity[w], extra_prev_side, NULL, hold_active);
                }
            }

            handle_piston_collisions(i, sub_dt);
            handle_boundary_collision(i);
        }

        // --- rebuild uniform grid ---
        for (int c = 0; c < max_cells; ++c) cell_head[c] = -1;

        for (int i = 0; i < n_active; ++i) {
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
        if (!cli_no_pp_collisions) for (int i = 0; i < n_active; ++i) {
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
        for (int i = 0; i < n_active; ++i) {
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
            for (int i = 0; i < n_active; ++i) {
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
            for (int i = 0; i < n_active; i++) {
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


/*
 * Hybrid face correction pass (for EDMD-hybrid mode)
 * --------------------------------------------------
 * After EDMD advances particle–particle (AB) collisions to t+dt ignoring
 * boundaries, correct any particle that ended beyond container walls,
 * inside/intersecting the divider slab, or overlapping pistons. For the
 * divider, apply finite-mass 1D elastic exchange along x to update both the
 * particle Vx and the wall velocity `vx_wall` (unless the wall is locked).
 * For container walls/pistons/extra walls, treat as infinite-mass mirrors.
 *
 * This pass does not advance time; it only adjusts positions/velocities at
 * the end of the step to eliminate tunneling and reflect impulses.
 */
static void hybrid_faces_correction_pass(float dt) {
    (void)dt; /* reserved for future use (e.g., diagnostics) */
    const float eps = 1e-6f;
    int n_active = particles_active > 0 ? particles_active : 0;

    for (int i = 0; i < n_active; ++i) {
        /* Pistons (infinite mass) */
        handle_piston_collisions(i, dt);

        /* Primary divider slab: finite/infinite mass depending on state */
        if (wall_enabled && WALL_THICKNESS > 0.0f) {
            float L = wall_x - 0.5f * WALL_THICKNESS;
            float R = wall_x + 0.5f * WALL_THICKNESS;
            float Ri = Radius[i] > 0.0f ? Radius[i] : PARTICLE_RADIUS;
            /* inside or overlapping slab interval? */
            if ((X[i] + Ri) > L && (X[i] - Ri) < R) {
                /* choose nearest face to push out */
                float dL = (X[i] + Ri) - L;
                float dR = R - (X[i] - Ri);
                int hit_left_face = (dL <= dR) ? 1 : 0;
                double u1 = Vx[i];
                double u2 = vx_wall;
                int locked_wall = (wall_hold_enabled && !wall_is_released);
                double v1, v2;
                if (locked_wall || WALL_MASS <= 0.0f) {
                    /* infinite-mass mirror in wall frame */
                    v1 = 2.0 * u2 - u1;
                    v2 = u2;
                } else {
                    double m = (double)PARTICLE_MASS;
                    double M = (double)WALL_MASS;
                    v1 = ((m - M) * u1 + 2.0 * M * u2) / (m + M);
                    v2 = ((M - m) * u2 + 2.0 * m * u1) / (m + M);
                }
                /* impulse accounting (opposite on wall) */
                double Jp = (double)PARTICLE_MASS * (v1 - u1);
                wall_impulse_x_accum -= (float)Jp;
                Vx[i] = (float)v1;
                vx_wall = (float)v2;

                if (hit_left_face) {
                    X[i] = L - Ri - eps;
                } else {
                    X[i] = R + Ri + eps;
                }
            }
        }

        /* Extra divider walls (treated as infinite mass mirrors) */
        if (extra_wall_count > 0 && extra_wall_positions) {
            for (int w = 0; w < extra_wall_count; ++w) {
                float wx = extra_wall_positions[w];
                float Lw = wx - 0.5f * WALL_THICKNESS;
                float Rw = wx + 0.5f * WALL_THICKNESS;
                float Ri = Radius[i] > 0.0f ? Radius[i] : PARTICLE_RADIUS;
                if ((X[i] + Ri) > Lw && (X[i] - Ri) < Rw) {
                    float dL = (X[i] + Ri) - Lw;
                    float dR = Rw - (X[i] - Ri);
                    int hit_left = (dL <= dR) ? 1 : 0;
                    /* infinite mass mirror with wall's instantaneous vx (if provided) */
                    double u1 = Vx[i];
                    double u2 = (extra_wall_velocity ? extra_wall_velocity[w] : 0.0f);
                    double v1 = 2.0 * u2 - u1;
                    Vx[i] = (float)v1;
                    if (hit_left) X[i] = Lw - Ri - eps; else X[i] = Rw + Ri + eps;
                }
            }
        }

        /* Container boundaries (static infinite mass) */
        handle_boundary_collision(i);
    }
}



////////////////////////////   ENERGY CALCULATION   ///////////////////////



// Expected kinetic energy for ideal gas: (3/2) N k_B T

// --> only applies if no interactions of particles
double kinetic_energy_expected() {
    return (3.0 / 2.0) * particles_active * K_B * TEMPERATURE;  // In Joules
}

double average_kinetic_energy_per_particle() {
    int n = particles_active > 0 ? particles_active : 1;
    return kinetic_energy() / n;
}

double kinetic_energy_total_system() {
    double ke_particles = kinetic_energy(); // existing function (particles only)
    double ke_wall = 0.5 * WALL_MASS * vx_wall * vx_wall;
    return ke_particles + ke_wall;
}




// /////// TMEPERATURE with dimension d ///////
//T = \frac{2}{d N k_B} \sum_{i=1}^{N} \frac{1}{2} m v_i^2
double temperature() {
    int n = particles_active > 0 ? particles_active : 1;
    return (2.0 / 2.0) * kinetic_energy() / (n * kB_effective());  // In reduced units
}





///////////////////////////   ENTROPY CALCULATION   ///////////////////////
/// VELOCITY
////////  single velocity histogram update ///////

void compute_velocity_histogram() {
    // Reset histogram
    for (int i = 0; i < NUM_BINS; i++) velocity_histogram[i] = 0;

    // Compute speeds and bin them
    int n = particles_active > 0 ? particles_active : 0;
    for (int i = 0; i < n; i++) {
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
    
    int n_ent = particles_active > 0 ? particles_active : 0;
    for (int i = 0; i < n_ent; i++) {
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
    int n = particles_active > 0 ? particles_active : 0;
    for (int i = 0; i < n; i++) {
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

            // Temperature controls (GUI): '[' to decrease, ']' to increase, '\\' to reset to 1.0
            case SDLK_LEFTBRACKET: { // decrease by 5%
                float T_new = fmaxf(0.01f, temperature_runtime * 0.95f);
                rescale_all_velocities_to_temperature(T_new);
                equalize_temperature_per_segment(temperature_runtime);
                update_dt_runtime_for_temperature();
                printf("🌡️ Temperature set to %.4f (−5%%)\n", temperature_runtime);
            } break;
            case SDLK_RIGHTBRACKET: { // increase by 5%
                float T_new = temperature_runtime * 1.05f;
                rescale_all_velocities_to_temperature(T_new);
                equalize_temperature_per_segment(temperature_runtime);
                update_dt_runtime_for_temperature();
                printf("🌡️ Temperature set to %.4f (+5%%)\n", temperature_runtime);
            } break;
            case SDLK_BACKSLASH: { // reset to 1.0
                rescale_all_velocities_to_temperature(1.0f);
                equalize_temperature_per_segment(temperature_runtime);
                update_dt_runtime_for_temperature();
                printf("🌡️ Temperature reset to 1.0000\n");
            } break;

            // Time scale controls: '-' halves, '=' doubles, '0' resets
            case SDLK_MINUS:
                time_scale_runtime = fmaxf(0.1f, time_scale_runtime * 0.5f);
                printf("Time scale set to %.2fx\n", (double)time_scale_runtime);
                break;
            case SDLK_EQUALS:
                time_scale_runtime = fminf(100.0f, time_scale_runtime * 2.0f);
                printf("Time scale set to %.2fx\n", (double)time_scale_runtime);
                break;
            case SDLK_0:
                time_scale_runtime = 1.0f;
                printf("Time scale reset to 1.00x\n");
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
                if (!piston_push_tracking) start_piston_push_tracking();
            }
            break;
        case SDLK_d:
            if (piston_left_x < piston_right_x - 20) {
                vx_piston_left += acceleration;
                if (vx_piston_left > max_velocity) vx_piston_left = max_velocity;
                if (!piston_push_tracking) start_piston_push_tracking();
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
            if (!piston_push_tracking) start_piston_push_tracking();
            break;
        case SDLK_RIGHT:
            if (piston_right_x < XW2) {
                vx_piston_right += acceleration;
                if (vx_piston_right > max_velocity) vx_piston_right = max_velocity;
            }
            piston_step_active = false;
            if (!piston_push_tracking) start_piston_push_tracking();
            break;
        case SDLK_UP: // Stop right piston immediately
            vx_piston_right = 0;
            piston_step_active = false;
            break;
            case SDLK_t: { // trigger piston step with protocol
            float travel = (float)SIM_WIDTH * 0.20f;
            float target = (float)XW2 - travel;
            float min_target = piston_left_x + 25.0f;
            if (target < min_target) target = min_target;
            piston_step_target = target;
            piston_step_speed = -fabsf((float)SIM_WIDTH) * 0.5f;
            vx_piston_right = piston_step_speed;
            piston_step_active = true;
            simulation_started = 1;
            
            // Initialize energy measurement only if requested or in energy_transfer experiment
            if (cli_enable_energy_measurement || cli_experiment_preset == EXPERIMENT_PRESET_ENERGY_TRANSFER) {
                initialize_energy_measurement();
            }
            // Start piston push energy tracking
            start_piston_push_tracking();
            
            printf("➡️  Piston step triggered: target %.2f px, protocol: %s\n", 
                   piston_step_target, 
                   cli_protocol == PROTOCOL_SIGMOIDAL ? "SIGMOIDAL" : 
                   cli_protocol == PROTOCOL_LINEAR ? "LINEAR" :
                   cli_protocol == PROTOCOL_SINUSOIDAL ? "SINUSOIDAL" : "STEP");
        }
            break;
            
            case SDLK_e: // Print energy measurement
                print_energy_measurement();
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
    float sigma = sqrt(kB_effective() * temperature / PARTICLE_MASS);
    float norm_factor = NUM_PARTICLES * (SIM_WIDTH / NUM_BINS);

    for (int i = 0; i < NUM_BINS; i++) {
        float v = MAX_VELOCITY * i / NUM_BINS;
        float f_v = (v * v) * expf(-PARTICLE_MASS * v * v / (2.0f * kB_effective() * temperature));
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
    float sigma = sqrtf(kB_effective() * measured_T / PARTICLE_MASS);

    float center_x = XW1 + L0_UNITS * PIXELS_PER_SIGMA;
    int bin_width = SIM_WIDTH / NUM_BINS;

    // === Histogram bins ===
    int binsX[NUM_BINS] = {0};
    int binsY[NUM_BINS] = {0};
    int binsSpeed[NUM_BINS] = {0};

    int max_binX = 1, max_binY = 1, max_binSpeed = 1;

    int n = particles_active > 0 ? particles_active : 0;
    for (int i = 0; i < n; i++) {
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

        float f1 = (PARTICLE_MASS * v1 / (kB_effective() * measured_T)) * expf(-PARTICLE_MASS * v1 * v1 / (2 * kB_effective() * measured_T));
        float f2 = (PARTICLE_MASS * v2 / (kB_effective() * measured_T)) * expf(-PARTICLE_MASS * v2 * v2 / (2 * kB_effective() * measured_T));

        int x1 = center_x + (i - 1 - NUM_BINS / 2) * bin_width;
        int x2 = center_x + (i - NUM_BINS / 2) * bin_width;

        int y1 = YW1 + SIM_HEIGHT + HIST_HEIGHT - (int)(f1 * normalization);
        int y2 = YW1 + SIM_HEIGHT + HIST_HEIGHT - (int)(f2 * normalization);

        SDL_SetRenderDrawColor(_renderer, 255, 255, 0, 255);
        draw_thick_line_vertical(_renderer, x1, y1, y2, 2); // e.g., thickness = 2
    }

    // --- Initial MB curve (red)
    float sigma0 = sqrtf(kB_effective() * temperature_runtime / PARTICLE_MASS);
    for (int i = 1; i < NUM_BINS; i++) {
        float v1 = ((float)(i - 1) / NUM_BINS) * MAX_VELOCITY;
        float v2 = ((float)i / NUM_BINS) * MAX_VELOCITY;

        float f1 = (PARTICLE_MASS * v1 / (kB_effective() * temperature_runtime)) * expf(-PARTICLE_MASS * v1 * v1 / (2 * kB_effective() * temperature_runtime));
        float f2 = (PARTICLE_MASS * v2 / (kB_effective() * temperature_runtime)) * expf(-PARTICLE_MASS * v2 * v2 / (2 * kB_effective() * temperature_runtime));

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

    // temperatures (use runtime active N)
    int n_active = (particles_active > 0) ? particles_active : 1;
    const double Tp  = KEp / (n_active * K_B);               // particles only
    const double Tw  = KEw / (0.5 * K_B);                    // 1 DOF wall
    const double Tt  = KEt / ((n_active + 0.5) * K_B);       // total incl. wall

    // format
    char line0[1280], line1[128], line2[128], line3[128], line4[128], line5[128];
    snprintf(line0, sizeof(line0), "T_particles: %.3f   T_total(incl wall): %.3f", Tp, Tt);
    snprintf(line1, sizeof(line1), "KE_tot: %.6e", KEt);
    snprintf(line2, sizeof(line2), "KE_particles: %.6e", KEp);
    snprintf(line3, sizeof(line3), "KE_wall: %.6e", KEw);
    snprintf(line4, sizeof(line4), "");
    snprintf(line5, sizeof(line5), "T_wall: %.3f", Tw);
    // Runtime clocks HUD: simulation time, dt, substeps, interactive time scale
    char line_clock[256];
    snprintf(line_clock, sizeof(line_clock), "time: %.3f  dt: %.2e  sub: %d  x%.2f",
             (double)simulation_time, (double)fixed_dt_runtime, SUBSTEPS, (double)time_scale_runtime);
    // Runtime thermo HUD: T, kB_eff, kBT, kBT1 mode
    char line_thermo[256];
    snprintf(line_thermo, sizeof(line_thermo), "T: %.3f  kB: %.6f  kBT: %.6f  kBT1:%s",
             (double)temperature_runtime, (double)kB_effective(), (double)kBT_effective(),
             cli_force_kbt_one ? "ON" : "OFF");
    char line6[256];
    if (piston_push_tracking) {
        double Wp_now = piston_work_left + piston_work_right;
        double dW     = Wp_now - piston_work_baseline;
        snprintf(line6, sizeof(line6), "PistonWork Δ: %.6e (max: %.6e)", dW, piston_work_delta_max);
    } else {
        snprintf(line6, sizeof(line6), "PistonWork Δ: (waiting for push)");
    }

    // draw HUD under the simulation window
    SDL_Color white = {255,255,255,255};
    int x = XW1 + 5;  // left margin
    int y = YW1 + SIM_HEIGHT + 5;  // just below the simulation area

    draw_text(r, f, line0, x, y, white);  y += 20;
    draw_text(r, f, line_clock, x, y, white); y += 20;
    draw_text(r, f, line_thermo, x, y, white); y += 20;
    if (num_internal_walls > 0) {
        double delta_sigma = leftmost_wall_current_delta / PIXELS_PER_SIGMA;
        char line_delta[128];
        snprintf(line_delta, sizeof(line_delta), "Δx_left_wall: %.4f σ", delta_sigma);
        draw_text(r, f, line_delta, x, y, white);  y += 20;
    }
    draw_text(r, f, line1, x, y, white);  y += 20;
    draw_text(r, f, line2, x, y, white);  y += 20;
    draw_text(r, f, line3, x, y, white);  y += 20;
    if (line4[0] != '\0') { draw_text(r, f, line4, x, y, white);  y += 20; }

    // Right column under sim: detailed energy + piston work + wall T
    int x2 = XW1 + SIM_WIDTH/2 + 10;
    int y2 = YW1 + SIM_HEIGHT + 5;
    // Wall temperature and piston work
    draw_text(r, f, line5, x2, y2, white);  y2 += 20;
    draw_text(r, f, line6, x2, y2, white);  y2 += 20;
    // Spring energy summary
    char line7[256];
    if (energy_measurement.spring_enabled) {
        snprintf(line7, sizeof(line7), "SpringE: %.6e (max: %.6e)",
                 (double)energy_measurement.spring_energy,
                 (double)energy_measurement.spring_energy_max);
    } else {
        snprintf(line7, sizeof(line7), "SpringE: (disabled)");
    }
    draw_text(r, f, line7, x2, y2, white); y2 += 20;

    // Detailed energy measurement HUD lines (if enabled)
    if (energy_measurement.spring_enabled) {
        char em0[256], em1[256], em2[256], em3[256], em4[256], em5[256];
        snprintf(em0, sizeof(em0), "SpringF: %.3e", (double)energy_measurement.spring_force);
        snprintf(em1, sizeof(em1), "SpringE_now: %.6e", (double)energy_measurement.spring_energy);
        snprintf(em2, sizeof(em2), "SpringE_max: %.6e", (double)energy_measurement.spring_energy_max);
        snprintf(em3, sizeof(em3), "E_transferred: %.6e", (double)energy_measurement.energy_transferred);
        snprintf(em4, sizeof(em4), "MaxDisp: %.3f px", (double)energy_measurement.max_displacement);
        snprintf(em5, sizeof(em5), "X: %.2f (X_eq: %.2f)",
                 (double)all_wall_positions[0], (double)energy_measurement.equilibrium_position);
        draw_text(r, f, em0, x2, y2, white); y2 += 20;
        draw_text(r, f, em1, x2, y2, white); y2 += 20;
        draw_text(r, f, em2, x2, y2, white); y2 += 20;
        draw_text(r, f, em3, x2, y2, white); y2 += 20;
        draw_text(r, f, em4, x2, y2, white); y2 += 20;
        draw_text(r, f, em5, x2, y2, white); y2 += 20;
    }

    y += 10;
    if (segment_count > 0 && segment_count <= MAX_DIVIDER_CAPACITY + 1 &&
        segment_counts && segment_ke) {
        for (int seg = 0; seg < segment_count; ++seg) {
            int count = segment_counts ? segment_counts[seg] : 0;
            double seg_ke = segment_ke ? segment_ke[seg] : 0.0;
            double seg_temp = 0.0;
            if (segment_temperature) {
                seg_temp = segment_temperature[seg];
            } else if (count > 0) {
                seg_temp = seg_ke / (count * K_B);
            }

            char seg_line[160];
            snprintf(seg_line, sizeof(seg_line),
                     "Box %d: N=%d KE=%.6e T=%.3f",
                     seg + 1, count, seg_ke, seg_temp);
            draw_text(r, f, seg_line, x, y, white);
            y += 20;
        }
    }
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
        num_steps = (int)(200000 / fmaxf(fixed_dt_runtime, 1e-9f));
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
                snprintf(filename, sizeof(filename), "%swall_x_positions_L0_%d_wallmassfactor_%d_run%d.csv",
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
                    update_particles_with_substepping(fixed_dt_runtime, &left_particles, &right_particles);

                    // wall velocity already updated inside ccd_wall_step; accumulator is purely diagnostic now.
                    // (optional) clear accumulator here if you prefer step-scoped use
                    // wall_impulse_x_accum = 0.f;

                    if (sim_mode != MODE_EDMD) update_pistons(fixed_dt_runtime);
                    if (sim_mode != MODE_EDMD) {
                        update_wall(fixed_dt_runtime, L0_UNITS);
                        update_extra_walls(fixed_dt_runtime, L0_UNITS);
                    }
                    sync_all_wall_positions();
                    update_wall_metrics();

                    // time bookkeeping
                    simulation_time += fixed_dt_runtime;

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

/*
 * run_energy_transfer_experiment
 * -------------------------------
 * Runs a single or repeated interactive-like batch with the spring-based
 * energy measurement enabled. Logs wall position and spring metrics over time
 * into experiments_energy_transfer. Honors multi-wall setups; primary wall is
 * the one coupled to the spring.
 */
static void run_energy_transfer_experiment(void) {
    system("mkdir -p experiments_energy_transfer");

    if (!energy_measurement.spring_enabled) {
        energy_measurement.spring_enabled = true;
        initialize_energy_measurement();
    }

    wall_enabled = 1;
    wall_position_is_managed_externally = 0;

    int left_particles = 0, right_particles = 0;

    char filename[512];
    snprintf(filename, sizeof(filename), "experiments_energy_transfer/energy_transfer_trace.csv");
    FILE *elog = fopen(filename, "w");
    if (!elog) { printf("❌ Could not open energy transfer log file.\n"); return; }
    fprintf(elog, "Time,Wall_X,SpringF,SpringE,EnergyTransferred,Left_Count,Right_Count\n");

    // Reset sim
    initialize_simulation_dimensions();
    initialize_simulation_parameters();
    initialize_simulation();
    steps_elapsed = 0;
    wall_is_released = false;
    simulation_time = 0.0f;
    wall_x_old = wall_x;
    missed_collision_events = 0; worst_penetration_observed = 0.0;

    int recorded_steps = 0;
    int target_steps = num_steps > 0 ? num_steps : (int)(3000.0f / fmaxf(fixed_dt_runtime, 1e-9f));

    while (recorded_steps < target_steps) {
        wall_x_old = wall_x;
        steps_elapsed++;

        wall_impulse_x_accum = 0.f;
        update_particles_with_substepping(fixed_dt_runtime, &left_particles, &right_particles);

        if (sim_mode != MODE_EDMD) update_pistons(fixed_dt_runtime);
        if (sim_mode != MODE_EDMD) {
            update_wall(fixed_dt_runtime, L0_UNITS);
            update_extra_walls(fixed_dt_runtime, L0_UNITS);
        }
        sync_all_wall_positions();
        update_wall_metrics();
        update_energy_measurement(fixed_dt_runtime);

        simulation_time += fixed_dt_runtime;

        if (steps_elapsed == wall_hold_steps && !wall_is_released) {
            wall_release_time = simulation_time;
            wall_is_released  = true;
            recorded_steps    = 0;
            printf("🔔 Wall released at step %d (t = %.3f)\n", steps_elapsed, wall_release_time);
        }

        if (!wall_is_released) continue;

        double t = simulation_time - wall_release_time;
        fprintf(elog, "%.6f, %.6f, %.6e, %.6e, %.6e, %d, %d\n",
                t,
                (double)wall_x,
                (double)energy_measurement.spring_force,
                (double)energy_measurement.spring_energy,
                (double)energy_measurement.energy_transferred,
                left_particles, right_particles);
        recorded_steps++;
    }

    fclose(elog);
    printf("✅ Energy transfer experiment done. Output → experiments_energy_transfer/energy_transfer_trace.csv\n");
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
    // If EDMD or hybrid mode selected, initialize EDMD state from current particle arrays
    if (sim_mode == MODE_EDMD || sim_mode == MODE_EDMD_HYBRID) {
        if (g_edmd) { edmd_destroy(g_edmd); g_edmd = NULL; }
        EDMD_Params prm = {0};
        prm.boxW = (double)(XW2 - XW1);
        prm.boxH = (double)(YW2 - YW1);
        prm.radius = (double)PARTICLE_RADIUS;
        prm.N = particles_active;
        prm.cell_size = 0.0; // auto
        /* divider/pistons parameters are unused in PP-only mode but safe to fill */
        prm.has_divider = (num_internal_walls > 0 ? 1 : 0);
        prm.divider_x = (double)(wall_x - XW1);
        prm.divider_thickness = (double)WALL_THICKNESS;
        prm.divider_mass = (double)(WALL_MASS / PARTICLE_MASS);
        prm.divider_vx = (double)vx_wall;
        int hasL = 1; int hasR = 1;
        double xL = (double)((piston_left_x + 5.0f) - XW1);
        double xR = (double)(piston_right_x - XW1);
        double vxL = (double)vx_piston_left;
        double vxR = (double)vx_piston_right;
        double mL = 0.0; double mR = 0.0;
        g_edmd = edmd_create(&prm);
        if (!g_edmd) { fprintf(stderr, "EDMD create failed.\n"); exit(1);} 
        edmd_config_pistons(g_edmd, hasL, xL, vxL, mL, hasR, xR, vxR, mR);
        // load current positions/velocities into EDMD state
        EDMD_Particle* P = (EDMD_Particle*)edmd_particles(g_edmd);
        for (int i=0;i<particles_active;i++) {
            double R = (double)PARTICLE_RADIUS;
            double bx = (double)(X[i] - XW1);
            double by = (double)(Y[i] - YW1);
            if (bx < R) bx = R; if (bx > prm.boxW - R) bx = prm.boxW - R;
            if (by < R) by = R; if (by > prm.boxH - R) by = prm.boxH - R;
            P[i].x = bx;
            P[i].y = by;
            P[i].vx = (double)Vx[i];
            P[i].vy = (double)Vy[i];
            P[i].coll_count = 0;
        }
        if (sim_mode == MODE_EDMD) {
            edmd_reschedule_all(g_edmd);
            printf("[EDMD] initialized: N=%d box=(%.1f,%.1f) R=%.3f\n", prm.N, prm.boxW, prm.boxH, prm.radius);
        } else {
            edmd_reschedule_all_pp_only(g_edmd);
            printf("[EDMD-HYBRID] initialized (PP-only): N=%d box=(%.1f,%.1f) R=%.3f\n", prm.N, prm.boxW, prm.boxH, prm.radius);
        }
    }

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
        accumulator += frame_time * time_scale_runtime;

        // events
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT) running = 0;
            keysSimulation(e);
            keysPiston(e);
        }

        int left_particles = 0, right_particles = 0;

        while (accumulator >= fixed_dt_runtime) {
            if (simulation_started && !paused) {
                // --- START OF STEP ---
                wall_x_old = wall_x;
                steps_elapsed++;

                // reset accumulator at the *start* of the step
                wall_impulse_x_accum = 0.f;

                if (sim_mode == MODE_EDMD) {
                    // advance EDMD to next time
                    // Keep pistons in sync with interactive controls; if velocities changed, update
                    static double prev_vxL = 0.0, prev_vxR = 0.0;
                    // Apply protocol-driven velocities (do not move positions here; EDMD owns positions)
                    apply_piston_protocol(fixed_dt_runtime);
                    if (vx_piston_left != prev_vxL || vx_piston_right != prev_vxR) {
                        const EDMD_Params* epc = edmd_params(g_edmd);
                        edmd_config_pistons(g_edmd,
                            1, epc->pistonL_x, (double)vx_piston_left, epc->pistonL_mass,
                            1, epc->pistonR_x, (double)vx_piston_right, epc->pistonR_mass);
                        edmd_reschedule_all(g_edmd);
                        prev_vxL = vx_piston_left; prev_vxR = vx_piston_right;
                    }
                    double t0 = edmd_time(g_edmd);
                    edmd_advance_to(g_edmd, t0 + (double)fixed_dt_runtime);
                    // copy back for rendering/diagnostics
                    const EDMD_Particle* P = edmd_particles(g_edmd);
                    for (int i=0;i<particles_active;i++) {
                        X[i] = (double)XW1 + (double)P[i].x;
                        Y[i] = (double)YW1 + (double)P[i].y;
                        Vx[i] = (double)P[i].vx;
                        Vy[i] = (double)P[i].vy;
                    }
                    const EDMD_Params* ep = edmd_params(g_edmd);
                    wall_x = (double)XW1 + ep->divider_x;
                    vx_wall = ep->divider_vx;
                    // Sync piston face positions for visualization (assume left piston width=5)
                    piston_left_x  = (float)(XW1 + ep->pistonL_x - 5.0);
                    piston_right_x = (float)(XW1 + ep->pistonR_x);
                    /* rebuild per-segment stats */
                    recompute_segment_stats_counts_and_temperature();
                } else if (sim_mode == MODE_EDMD_HYBRID) {
                    /* PP-only EDMD, followed by CCD-style face correction */
                    /* Keep pistons synced (EDMD ignores faces, but velocities affect UI) */
                    {
                        /* optional: if velocities changed, no EDMD reschedule needed in PP-only */
                    }
                    /* Advance PP-only to t+dt */
                    double t0 = edmd_time(g_edmd);
                    edmd_advance_pp_only_to(g_edmd, t0 + (double)fixed_dt_runtime);
                    /* Copy particle state back for rendering + face correction */
                    const EDMD_Particle* P = edmd_particles(g_edmd);
                    for (int i=0;i<particles_active;i++) {
                        X[i] = (double)XW1 + (double)P[i].x;
                        Y[i] = (double)YW1 + (double)P[i].y;
                        Vx[i] = (double)P[i].vx;
                        Vy[i] = (double)P[i].vy;
                    }
                    /* Correct against faces (container/divider/pistons) and update wall vx */
                    hybrid_faces_correction_pass(fixed_dt_runtime);
                    /* Write corrected state back into EDMD and reschedule AB for next step */
                    EDMD_Particle* Pw = (EDMD_Particle*)edmd_particles(g_edmd);
                    for (int i=0;i<particles_active;i++) {
                        Pw[i].x = (double)(X[i] - XW1);
                        Pw[i].y = (double)(Y[i] - YW1);
                        Pw[i].vx = (double)Vx[i];
                        Pw[i].vy = (double)Vy[i];
                        /* coll_count left as-is to keep event invalidation; full reschedule will rebuild */
                    }
                    edmd_reschedule_all_pp_only(g_edmd);
                    /* HUD: counts + KE + T */
                    recompute_segment_stats_counts_and_temperature();
                } else {
                    // time-driven (Euler ballistic with CCD)
                    update_particles_with_substepping(fixed_dt_runtime, &left_particles, &right_particles);
                }

                // wall velocity is already updated inside ccd_wall_step (impulses applied immediately).
                // wall_impulse_x_accum is kept for diagnostics only.

                if (sim_mode != MODE_EDMD) update_pistons(fixed_dt_runtime);
                if (sim_mode != MODE_EDMD) {
                    update_wall(fixed_dt_runtime, L0_UNITS);
                    update_extra_walls(fixed_dt_runtime, L0_UNITS);
                }
                sync_all_wall_positions();
                update_wall_metrics();
                update_energy_measurement(fixed_dt_runtime);

                simulation_time += fixed_dt_runtime;

                // Update piston push energy tracking (interactive mode)
                if (piston_push_tracking) {
                    double Wp = piston_work_left + piston_work_right;
                    double delta = Wp - piston_work_baseline;
                    if (delta > piston_work_delta_max) {
                        piston_work_delta_max = delta;
                    }
                }

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

                // inside the fixed-step block, after simulation_time += fixed_dt_runtime
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
            accumulator -= fixed_dt_runtime;
        }

        // === RENDERING ===
        draw_clear_screen();
        draw_coordinate_system(_renderer);
        draw_simulation_boundary();
        render_particles();
        render_pistons();
        if (wall_enabled) draw_wall();
        // Draw spring visualization on leftmost segment
        render_left_spring();


        // <-- add this line
        render_energy_hud(_renderer, font);

        SDL_Color segment_palette[] = {
            {255, 120, 0, 255},
            {0, 180, 255, 255},
            {170, 255, 0, 255},
            {255, 80, 180, 255},
            {240, 240, 240, 255}
        };
        const int palette_size = (int)(sizeof(segment_palette) / sizeof(segment_palette[0]));
        SDL_Color yellow = {255, 255, 0, 255};

        if (segment_counts && segment_count > 0 && segment_count <= MAX_DIVIDER_CAPACITY + 1) {
            float left_bounds[MAX_DIVIDER_CAPACITY + 1];
            float right_bounds[MAX_DIVIDER_CAPACITY + 1];
            compute_segment_bounds(left_bounds, right_bounds);

            for (int seg = 0; seg < segment_count; ++seg) {
                int count = segment_counts[seg];
                char seg_buffer[64];
                snprintf(seg_buffer, sizeof(seg_buffer), "Box %d: %d", seg + 1, count);

                float seg_width = right_bounds[seg] - left_bounds[seg];
                float seg_x = left_bounds[seg] + fmaxf(6.0f, seg_width * 0.05f);
                int text_x = (int)seg_x;
                int text_y = YW1 + 10;

                SDL_Color color = segment_palette[seg % palette_size];
                draw_text(_renderer, font, seg_buffer, text_x, text_y, color);

                double ke_val = (segment_ke) ? segment_ke[seg] : 0.0;
                double temp_val = 0.0;
                if (segment_temperature) {
                    temp_val = segment_temperature[seg];
                } else if (count > 0) {
                    temp_val = ke_val / (count * K_B);
                }
                char ke_buffer[80];
                snprintf(ke_buffer, sizeof(ke_buffer), "KE=%.3e T=%.2f", ke_val, temp_val);
                draw_text(_renderer, font, ke_buffer, text_x, text_y + 18, color);
            }
        } else {
            SDL_Color white = {255, 255, 255, 255};
            char fallback[64];
            snprintf(fallback, sizeof(fallback), "Left: %d", left_particles);
            draw_text(_renderer, font, fallback, XW1 + 5, YW1 + 10, white);
            snprintf(fallback, sizeof(fallback), "Right: %d", right_particles);
            draw_text(_renderer, font, fallback, XW2 - 120, YW1 + 10, white);
        }

        char buffer[64];

        if (wall_hold_enabled && !wall_is_released) {
            snprintf(buffer, sizeof(buffer), "Wall held: %d steps left", wall_hold_steps - steps_elapsed);
            draw_text(_renderer, font, buffer, XW1 + 5, YW1 + 50, yellow);
        }

        render_velocity_histograms();
        float T_measured = compute_measured_temperature_from_ke();
        char temp_label[128];
        snprintf(temp_label, sizeof(temp_label), "T_measured: %.2f", T_measured);
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
    // Set dt runtime based on preset and temperature
    update_dt_runtime_for_temperature();
    // Initialize energy measurement only if explicitly enabled or via preset
    if (energy_measurement.spring_enabled) {
        initialize_energy_measurement();
    }
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
    if (!cli_force_no_experiments && cli_experiment_preset == EXPERIMENT_PRESET_ENERGY_TRANSFER) {
        run_energy_transfer_experiment();
    } else if (enable_speed_of_sound_experiments && !cli_force_no_experiments) {
        run_speed_of_sound_experiments();
    } else {
        simulation_loop();  // SDL interactive
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

#include "edmd_core/edmd.h"
