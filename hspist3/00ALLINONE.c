#include <SDL2/SDL.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <errno.h>
#include <stdlib.h>
#include <float.h>
#include <limits.h>
#include <unistd.h>
#include <ctype.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <SDL2/SDL_ttf.h>
#include "kissfft/kiss_fft.h"
#include "kissfft/kiss_fftr.h"
#include "edmd_core/edmd.h"
#include "edmd_core/edmd_accelerated.h"



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

// Persist argv for per-run metadata (e.g., 00_COMMAND.md in output folders).
static int g_main_argc = 0;
static char **g_main_argv = NULL;

static void rebuild_internal_walls(void);
static void free_internal_wall_arrays(void);
static void reset_simulation_state(void);
void initialize_simulation(void);
static void szilard_prepare_interactive_config(void);
static void szilard_setup_interactive_state(void);
static void szilard_measure_memory_from_species(void);
static int szilard_particle_gate_label(int i, int use_memory);
static void szilard_speed_stats(double *mean_speed, double *std_speed);
static int szilard_particle_hot_enough_for_gate_leak(int i);
static double szilard_memory_gate_information_bits_at_split(float split_x);
static double szilard_partial_uncertainty_fraction_now(void);
static void szilard_update_partial_uncertainty_peak(void);
static void szilard_memory_side_counts_at_split(float split_x,
                                                int *mL0, int *mL1, int *mLU,
                                                int *mR0, int *mR1, int *mRU);
static double szilard_mutual_information_3state_nats(int nL0, int nL1, int nLU,
                                                     int nR0, int nR1, int nRU);
static void szilard_apply_species_flips(double dt_sigma);
static void szilard_apply_memory_erase_min_heat(void);
static void szilard_begin_separation(void);
static void szilard_begin_recombine(void);
static void szilard_begin_recombine_wrong(void);
static void szilard_begin_recombine_partial(void);
static void szilard_manual_erase_memory(void);
static void szilard_update_interactive_poststep(void);
static void szilard_sync_edmd_from_globals(void);
static void szilard_clamp_particles_left_half(void);
static void szilard_save_branch_snapshot(void);
static void szilard_restore_branch_snapshot(void);
static void szilard_rebuild_edmd_from_current_globals(void);
static void szilard_apply_interactive_velocity_current_phase(void);
static int szilard_particle_moving_gate_target_side(int i, int *target_side_out);
static int szilard_particle_can_cross_fixed_gate_right_to_left(int i);
static int szilard_repair_edmd_particles_for_pistons_and_pockets(void);
static const char* szilard_particle_view_mode_label(int mode);
static double szilard_reference_kbt(void);
static double szilard_available_free_energy_from_nats(double info_nats, int total_count);
static double szilard_current_visible_free_energy_from_counts(int nL0, int nL1, int nR0, int nR1);
static void szilard_reset_visible_load_state(double fx_start);
static void szilard_update_visible_load_state(double fx_now);
static void run_simple_prediction_box_experiment(void);
static void simple_box_layout_metrics(double *mi_bits, double *H_bits, double *H_norm);

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
// Bump high so we can run large-N tests (targeting 100k scale).
#define NUM_PARTICLES 120000     // MAX_PARTICLES (buffer size)
#define DEFAULT_PARTICLES_ACTIVE 100
static int particles_active = DEFAULT_PARTICLES_ACTIVE;  // active particles used by seeding/HUD/thermo
bool log_packing_fraction = 1; // 1 = log packing fraction, 0 = don't log




#if EXPERIMENT_MODE
    bool enable_speed_of_sound_experiments = true;
#else
    bool enable_speed_of_sound_experiments = false;
#endif



// --- Physical constants ---
#if MOLECULAR_MODE
    #define PARTICLE_RADIUS_UNIT 0.5f      // unit length               normalized radius = 1 (works in both modes!)
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
        #define FIXED_DT 0.04f  // in normalized time units They don't state the exact dt, but normalized simulations use dt = 10^{-3} often.
        #define PIXELS_PER_SIGMA 24.0f      // rendering scale (slightly larger for interactive SDL readability)
        #define TEMPERATURE 100.0f            // default target temperature (reduced units)
        #define K_B 1/TEMPERATURE 

        #define SUBSTEPS 1   // default, can be changed externally
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
        #define PIXELS_PER_SIGMA 24.0f      // rendering scale (slightly larger for interactive SDL readability)
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

// Optional: visual thickness (render-only). If <0, falls back to physics thickness.
static float wall_thickness_visual_runtime = -1.0f; // pixels



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
int XW1 = 200;   // left margin for HUD/padding
int XW2, YW1 = 0, YW2;
int MAX_X, MAX_Y;
int HIST_WIDTH, HIST_HEIGHT;

#define MIN_WINDOW_WIDTH 1300  // Minimum width for consistent UI (matches L0=20)

void initialize_simulation_dimensions() {
    SIM_WIDTH = (int)(2 * L0_UNITS * PIXELS_PER_SIGMA);
    SIM_HEIGHT = (int)(HEIGHT_UNITS * PIXELS_PER_SIGMA);

    XW2 = XW1 + SIM_WIDTH;
    YW2 = YW1 + SIM_HEIGHT;

    // Expand padding on all sides for HUD/plots; left margin handled via XW1
    // Extra space for HUD (right) and distribution panels (bottom).
    // This is pure rendering/UI padding; it does not change physics (box size stays SIM_WIDTH×SIM_HEIGHT).
    int pad_x = 560;
    // Default bottom space for live distributions and stacked F2/F3 control panels.
    // This is UI-only padding; the physical simulation box dimensions above are unchanged.
    int pad_y = 640;

    // Calculate base window width
    int calculated_width = XW2 + pad_x;

    // Enforce minimum width for consistent UI (especially for small L0 values)
    MAX_X = (calculated_width < MIN_WINDOW_WIDTH) ? MIN_WINDOW_WIDTH : calculated_width;
    MAX_Y = YW2 + pad_y;

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
static float piston_step_direction = -1.0f;

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

// Optional: stop the right piston automatically once a target input work is reached.
// Work is measured as ΔW = (W_pistonL + W_pistonR) - baseline at push start (same units as HUD).
static int    cli_piston_work_target_set = 0;   // --piston-work-target=...
static double cli_piston_work_target = 0.0;     // target ΔW (>=0)

// Diagnostics for potential collision misses (large overlaps after integration)
static long   missed_collision_events      = 0;
static double worst_penetration_observed   = 0.0;
static bool   time_mode_walls_integrated_this_step = false;

// Divider wall configuration
static int    num_internal_walls = 1;
static int    cli_requested_walls = 1;
static int    primary_wall_index = 0;
static int    leftmost_wall_index = 0;
static float *all_wall_positions = NULL;            // length = num_internal_walls
static float *all_wall_masses = NULL;               // length = num_internal_walls (absolute mass units)
static int    extra_wall_count = 0;
static float *extra_wall_positions = NULL;          // excludes primary wall
static float *extra_wall_old_positions = NULL;
static float *extra_wall_velocity = NULL;           // mostly zero (static walls)
static float *extra_wall_masses = NULL;             // excludes primary wall (absolute mass units)
static bool  *extra_wall_hold_enabled = NULL;
static bool  *extra_wall_is_released = NULL;
static double *extra_wall_release_time = NULL;
static int    segment_count = 0;
static int   *segment_counts = NULL;                // length = num_internal_walls + 1
static double *segment_ke = NULL;
static double *segment_temperature = NULL;
static double *segment_vx_mean = NULL;
static double *segment_vx_var = NULL;
static double *segment_vy_mean = NULL;
static double *segment_vy_var = NULL;
static double *segment_speed_mean = NULL;
static double *segment_speed_var = NULL;

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
// Render-only wall thickness override (in sigma units); physics uses --wall-thickness
static float  cli_override_wall_thickness_vis_sigma = -1.0f;
static float  cli_override_temperature = -1.0f;  // --temperature=VALUE (reduced units)
static float  cli_override_hb_temp = -1.0f;     // ##CHRIS: --hb-temp=VALUE (heat bath temperature)
// Particle size CLI
static int    cli_particle_radius_set = 0;     // 1 if user set --particle-radius (sigma units)
static float  cli_particle_radius_sigma = 0.0f;
static int    cli_particle_diameter_set = 0;   // 1 if user set --particle-diameter (sigma units)
static float  cli_particle_diameter_sigma = 0.0f;
// ##CHRIS: Packing fraction CLI - auto-calculate particle radius for desired eta
static bool   cli_packing_fraction_set = false;
static float  cli_packing_fraction = -1.0f;    // --eta=VALUE or --packing-fraction=VALUE
static bool   cli_force_no_experiments = false;
static bool   cli_enable_energy_measurement = false;
static int    cli_override_particles = -1;   // --particles=N
static bool   cli_headless = false;          // --no-gui or --headless: run without SDL window
static bool   cli_auto_piston_step = false;  // auto-trigger piston step (same as pressing 't')
static char  *cli_command_line = NULL;       // argv joined for logging
static char  *cli_speed_sound_run_dir = NULL; // --speed-sound-run-dir=PATH
static bool   cli_seed_set = false;
static unsigned int cli_seed = 1;
// Energy transfer summary output (optional override)
static bool   cli_energy_transfer_summary_set = false;
static char   cli_energy_transfer_summary_path[512] = {0};
static bool   cli_energy_transfer_trace_set = false;
static char   cli_energy_transfer_trace_path[512] = {0};
// Output directories (resolved from argv[0])
static char g_program_dir[512] = ".";
static char g_launch_cwd[1024] = ".";
static char g_energy_transfer_dir[512] = "experiments_energy_transfer";
// Seeding mode CLI
typedef enum { SEEDING_GRID=0, SEEDING_HONEYCOMB=1, SEEDING_RANDOM=2 } seeding_mode_t;
static seeding_mode_t cli_seeding_mode = SEEDING_GRID; // default grid for reproducibility
static bool cli_left_empty = false;          // seed zero particles in leftmost segment
static bool cli_no_pp_collisions = false;    // disable particle-particle collisions
static bool cli_szilard_no_pp_during_g = false; /* ##CHRIS: disable PP during Szilard G phase only */
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
    EXPERIMENT_PRESET_SIMPLE_BOX_PREDICTION,
    EXPERIMENT_PRESET_ATP_SYNTHASE_ANALOG,
    EXPERIMENT_PRESET_PROTOCOL_OPTIMIZATION,
    EXPERIMENT_PRESET_DYNAMIC_IB,
    // ##CHRIS: Particle-life artificial life simulation
    EXPERIMENT_PRESET_PARTICLELIFE
} ExperimentPreset;

// Protocol types for different piston movements
typedef enum {
    PROTOCOL_STEP = 0,           // Instant acceleration (unphysical)
    PROTOCOL_SIGMOIDAL,         // Smooth sigmoidal acceleration
    PROTOCOL_LINEAR,            // Linear acceleration
    PROTOCOL_SINUSOIDAL,        // Sinusoidal protocol
    PROTOCOL_OPTIMAL              // AI-optimized protocol
} PistonProtocol;

// Piston protocol parameters (σ-time units; velocities are numerically σ/σ-time)
static float piston_protocol_start_time = -1.0f;
static float piston_protocol_duration = 2.0f;        // how long to ramp to max speed
static float piston_protocol_max_speed = 1.0f;       // magnitude of max speed (σ/σ-time)
static float piston_protocol_v0 = 0.0f;              // initial speed for linear/sigmoidal (σ/σ-time)
static float piston_protocol_linear_gradient = 0.0f; // σ / (σ-time)^2
static bool  piston_protocol_linear_use_gradient = false;
static float piston_protocol_sigmoid_steepness = 0.0f; // dimensionless; if <=0 uses default

// CLI overrides for piston protocol parameters
static int   cli_piston_speed_set = 0;
static float cli_piston_speed = 1.0f;                // magnitude (positive)
static int   cli_piston_duration_set = 0;
static float cli_piston_duration = 2.0f;             // σ-time units
static int   cli_piston_travel_set = 0;
static float cli_piston_travel_sigma = 0.0f;         // σ units; if unset we'll use 20% of box width

// Right piston protocol (more explicit names for energy experiments)
static bool  cli_piston_right_mode_set = false;
static char  cli_piston_right_mode[32] = {0}; // "step"|"linear"|"sigmoidal"
static bool  cli_piston_right_travel_set = false;
static float cli_piston_right_travel_sigma = 0.0f;
static bool  cli_piston_right_step_speed_set = false;
static float cli_piston_right_step_speed = 0.0f;
static bool  cli_piston_right_v0_set = false;
static float cli_piston_right_v0 = 0.0f;
static bool  cli_piston_right_gradient_set = false;
static float cli_piston_right_gradient = 0.0f; // σ / (σ-time)^2
static bool  cli_piston_right_vmax_set = false;
static float cli_piston_right_vmax = 0.0f;
static bool  cli_piston_right_duration_set = false;
static float cli_piston_right_duration = 0.0f; // σ-time units
static bool  cli_piston_right_sigmoid_steepness_set = false;
static float cli_piston_right_sigmoid_steepness = 0.0f; // dimensionless "k*duration" style

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

// ##CHRIS: Particle-life tensor file loading
static char* cli_pl_tensor_file = NULL;  // Path to custom tensor file
static int   cli_pl_cyclic_boundaries = -1;  // -1 = use tensor file, 0/1 = CLI override
static double cli_pl_dt = -1.0;  // Particle-life timestep override (-1 = use default)

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
static float cli_wall_masses[5] = {200.0f, 200.0f, 200.0f, 200.0f, 200.0f}; // mass factors (multiples of particle mass)
static int   cli_wall_masses_set = 0;   // becomes 1 only if --wall-mass-factors provided
static int   cli_wall_masses_count = 0; // how many factors were provided
static int   cli_num_walls = 1;  // default number of walls
static float cli_wall_positions[5] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f}; // wall positions
static PistonProtocol cli_protocol = PROTOCOL_STEP; // default protocol (constant speed once triggered)
// Optional: per-segment temperature overrides
static float *cli_temperature_segments = NULL;        // from --temperature list
static size_t cli_temperature_segments_count = 0;
// CSV logging cadence control
//  - cli_output_dt <= 1.0 :
//        * TIME/RK4 modes → per-substep logging (highest resolution)
//        * EDMD/HYBRID    → per-step logging (highest resolution)
//  - cli_output_dt  > 1.0 : log every output_dt σ-time units (all modes)
static float cli_output_dt = 1.0f;
// Override collision substepping at runtime (<=0 => use compile-time SUBSTEPS).
static int cli_substeps_override = -1;
// Adaptive collision substepping (auto): choose substeps per fixed step based on vmax*dt.
// This keeps runs fast when velocities are low, and increases resolution only during shocks/collapse.
static bool  cli_substeps_auto = false;
static int   cli_substeps_auto_min = 1;
static int   cli_substeps_auto_max = 200;
// Max *relative* displacement per substep as a fraction of particle diameter.
// Heuristic: require (2*vmax*dt/substeps) <= dx_frac * (2R).
static float cli_substeps_auto_dx_frac = 0.25f;
// Last substep count actually used by the CCD time-mode integrator (incl. adaptive mode).
// This is for HUD + diagnostics logs; it is not the compile-time SUBSTEPS constant.
static int g_last_substeps_used = SUBSTEPS;
// Stop interactive run automatically after this many sim units post-release
static float cli_time_span = -1.0f; // <0 => unlimited; >0 => stop when t_rel >= cli_time_span

// ##CHRIS: Global pointers for substep-level CSV logging
static FILE *g_substep_log = NULL;
static double g_substep_base_time = 0.0;

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

static void gui_recompute_hist_extents(void) {
    HIST_WIDTH = MAX_X - SIM_WIDTH;
    HIST_HEIGHT = MAX_Y - SIM_HEIGHT;
    if (HIST_WIDTH < 1) HIST_WIDTH = 1;
    if (HIST_HEIGHT < 1) HIST_HEIGHT = 1;
}

static void gui_sync_renderer_canvas_size(void) {
    int w = 0, h = 0;
    if (_renderer && SDL_GetRendererOutputSize(_renderer, &w, &h) == 0 && w > 0 && h > 0) {
        MAX_X = w;
        MAX_Y = h;
        gui_recompute_hist_extents();
        return;
    }
    if (_window) {
        SDL_GetWindowSize(_window, &w, &h);
        if (w > 0 && h > 0) {
            MAX_X = w;
            MAX_Y = h;
            gui_recompute_hist_extents();
        }
    }
}

static void gui_mouse_to_canvas_xy(int raw_x, int raw_y, int *canvas_x, int *canvas_y) {
    double sx = 1.0;
    double sy = 1.0;
    if (_window && _renderer) {
        int ww = 0, wh = 0, rw = 0, rh = 0;
        SDL_GetWindowSize(_window, &ww, &wh);
        if (SDL_GetRendererOutputSize(_renderer, &rw, &rh) == 0 && ww > 0 && wh > 0 && rw > 0 && rh > 0) {
            sx = (double)rw / (double)ww;
            sy = (double)rh / (double)wh;
        }
    }
    if (canvas_x) *canvas_x = (int)lround((double)raw_x * sx);
    if (canvas_y) *canvas_y = (int)lround((double)raw_y * sy);
}

// GUI toggles (interactive only)
// 0 = off, 1 = global velocity hists, 2 = per-segment velocity hists (aligned to compartments)
// GUI distribution panel mode:
// 0 = off
// 1 = velocity (global)
// 2 = velocity components per-segment
// 3 = stacked per-segment: x-density + velocity components (+ speed)
static int gui_dist_mode = 3;
TTF_Font *font = NULL;
static TTF_Font *font_small = NULL; // UI-only: used for in-panel histogram stats

// For smoother GUI rendering with fixed-step solvers: time remainder (accumulator) after last physics step.
// Used ONLY for visualization (interpolating boundaries/histogram windows); does not affect physics.
static float gui_render_dt_rem = 0.0f;

// GUI window sizing (pure UI; never changes physics).
// Some displays (e.g. external monitors) have limited vertical resolution; allow shrinking the bottom plots.
// Use PageUp/PageDown at runtime.
static int gui_window_min_extra_y = 240;   // minimum pixels reserved below SIM_HEIGHT
static int gui_window_height_step = 80;    // pixels per PgUp/PgDn

// GUI toggle: freeze velocity display x-range (vmax) so heating/cooling is visually obvious.
// When enabled, velocity panels use a fixed vmax for all segments instead of per-segment auto-scaling.
static int gui_vel_range_freeze = 0;
static float gui_vel_range_vmax_fixed = 0.0f;

// GUI toggle: overlay a fitted distribution curve on velocity histograms (vx/vy Gaussian, |v| Rayleigh).
// Toggle with 'c'. Purely visual; does not affect physics.
static int gui_vel_fit_overlay = 0;

// GUI toggle: Maxwell-Boltzmann style speed view.
// Purely visual: particles and speed histogram bars are colored by |v|.
static int gui_mb_speed_view = 1;
static float gui_mb_speed_vmax_ref = 0.0f;

// GUI toggle: draw a "uniform baseline" overlay for x-density panels (expected flat density if well-mixed).
// Toggle with 'u'. Purely visual; does not affect physics.
static int gui_density_uniform_overlay = 1;

// GUI toggle: overlay a smoothed "fit" curve on the global x-density panel (like velocity fit curves).
// Toggle with 'g'. Purely visual; does not affect physics.
static int gui_density_curve_overlay = 1;

// GUI toggle: show per-compartment x-density panel in stacked distributions mode.
// Toggle with 'x'. Default off because it can be visually misleading during compression.
static int gui_show_x_density_by_segment = 0;

// GUI toggle: "EXACT" histogram mode (no smoothing, no fit/baseline overlays).
// Goal: each particle contributes to exactly one bin; visuals match raw bin counts.
// Toggle with 'z'. Purely visual; does not affect physics.
static int gui_hist_exact = 0;

// GUI toggle: show/hide vertical marker lines in distribution panels
// (wall boundaries, piston boundary, anchor line, vx=0 line, mean/median markers).
// Toggle with 'l'. Purely visual; does not affect physics.
static int gui_dist_markers = 1;

// Forward declarations (used by histogram renderers before the function bodies appear)
void draw_text(SDL_Renderer *renderer, TTF_Font *font, const char *text, int x, int y, SDL_Color color);
static void draw_text_clipped(SDL_Renderer *renderer, TTF_Font *font, const char *text, int x, int y, SDL_Color color, int max_w);
static float gui_compute_vmax_fixed(void);
static void gui_adjust_window_height(int delta);
static float gui_mb_speed_vmax_for_render(void);
static void gui_mb_capture_speed_scale_now(void);
static void live_apply_spring_enabled(int enabled);
static void live_rebuild_edmd_from_current_state(void);
static int live_controls_handle_event(const SDL_Event *e);
static void render_live_controls(SDL_Renderer *renderer, TTF_Font *font);
static void live_controls_layout(int *x, int *y, int *w, int *h);
static void live_scene_layout(int *x, int *y, int *w, int *h);

// Global variables for simulation control
int simulation_started = 0;  // Flag to track simulation start

// Forward declarations for helpers referenced in temperature tools
static int segment_index_for_position(double px);
double kinetic_energy(void);

// Runtime temperature control (defaults to compile-time TEMPERATURE)
static float temperature_runtime = TEMPERATURE;
static float fixed_dt_runtime = FIXED_DT; // dt actually used by the integrator
static float time_scale_runtime = 1.0f;   // scales real-time → sim-time in interactive mode
static float gui_particle_draw_min_px = 1.0f; // render-only lower bound so tiny particles remain visible
// Integrator / simulation mode
typedef enum { MODE_TIME=0, MODE_RK4=1, MODE_EDMD=2, MODE_EDMD_HYBRID=3 } sim_mode_t;
static sim_mode_t sim_mode = MODE_TIME;

//============================================================================
// ##CHRIS: BEGIN HEAT BATH VARIABLES AND CONFIGURATION
//============================================================================
static int heatbath_enabled = 0;              // Toggle with 'b' key
static float heatbath_temperature = 1.0f;     // Set to gas temperature when activated
static int thermal_wall_mode = 2;             // 0=gradual damping, 1=base MB, 2=adaptive MB
static float mb_overshoot_factor = 2.0f;      // Overshoot for faster convergence
static float stability_window_percent = 0.025f; // 2.5% stability window
static float current_sample_temp = 0.0f;      // For diagnostic display
//============================================================================
// ##CHRIS: END HEAT BATH VARIABLES
//============================================================================

//============================================================================
// ##CHRIS: BEGIN ANDERSEN THERMOSTAT VARIABLES
//============================================================================
static int andersen_enabled = 0;                     // Toggle with 'a' key
static double andersen_collision_freq = 0.1;        // Collision frequency (collisions per particle per time unit)
static int* last_andersen_collision_step = NULL;    // Track when each particle was last resampled (allocated dynamically)
//============================================================================
// ##CHRIS: END ANDERSEN THERMOSTAT VARIABLES
//============================================================================

//============================================================================
// Live SDL control surface
//============================================================================
typedef enum {
    LIVE_PARAM_FLOAT = 0,
    LIVE_PARAM_DOUBLE = 1
} LiveParamType;

typedef struct {
    const char *label;
    const char *unit;
    LiveParamType type;
    void *target;
    double min_value;
    double max_value;
    double default_value;
    int logarithmic;
} LiveSlider;

static int live_controls_visible = 1;
static int live_scene_visible = 1;
static int live_controls_drag_index = -1;
static int live_scene_drag_index = -1;
static int live_scene_drag_visible_row = -1;
static int live_scene_rebuild_requested = 0;
#define LIVE_SCENE_MAX_WALLS 5
#define LIVE_SCENE_ROW_COUNT (5 + LIVE_SCENE_MAX_WALLS)
typedef enum {
    LIVE_SCENE_MODE_ENERGY_TRANSFER = 0,
    LIVE_SCENE_MODE_SZILARD = 1,
    LIVE_SCENE_MODE_SPEED_OF_SOUND = 2,
    LIVE_SCENE_MODE_COUNT = 3
} LiveSceneMode;
typedef enum {
    LIVE_SCENE_ENGINE_TIME_TOI = 0,
    LIVE_SCENE_ENGINE_HYBRID = 1,
    LIVE_SCENE_ENGINE_EDMD = 2
} LiveSceneEngine;
static LiveSceneMode live_scene_mode = LIVE_SCENE_MODE_ENERGY_TRANSFER;
static int live_scene_mode_dropdown_open = 0;
static LiveSceneEngine live_scene_engine = LIVE_SCENE_ENGINE_EDMD;
static int live_scene_wall_count = 2;
static int live_scene_particle_count = DEFAULT_PARTICLES_ACTIVE;
static float live_scene_eta = 0.20f;
static float live_scene_radius_sigma = 0.10f;
static float live_scene_spring_l0_sigma = 20.0f;
static float live_scene_wall_mass[LIVE_SCENE_MAX_WALLS] = {200.0f, 2000.0f, 200.0f, 200.0f, 200.0f};

static LiveSlider live_sliders[] = {
    { "time scale", "x", LIVE_PARAM_FLOAT,  &time_scale_runtime,       0.10,   256.0, 1.0,      1 },
    { "fixed dt",   "",  LIVE_PARAM_FLOAT,  &fixed_dt_runtime,         0.001,  0.080, FIXED_DT,  1 },
    { "piston vmax","",  LIVE_PARAM_FLOAT,  &piston_protocol_max_speed,0.10,   200.0, 1.0,      1 },
    { "piston ramp","t", LIVE_PARAM_FLOAT,  &piston_protocol_duration, 0.10,   20.0,  2.0,      1 },
    { "HB target T","",  LIVE_PARAM_FLOAT,  &heatbath_temperature,     0.010,  20.0,  1.0,      1 },
    { "draw min r", "px",LIVE_PARAM_FLOAT,  &gui_particle_draw_min_px, 0.0,    8.0,   1.0,      0 }
};

static const int live_slider_count = (int)(sizeof(live_sliders) / sizeof(live_sliders[0]));

static float gui_scene_view_scale = 1.0f;
static int gui_scene_saved_max_x = 0;
static int gui_scene_saved_max_y = 0;
static int gui_scene_saved_hist_w = 0;
static int gui_scene_saved_hist_h = 0;
static int gui_scene_view_active = 0;

static float gui_renderer_device_scale(void) {
    if (!_window || !_renderer) return 1.0f;
    int ww = 0, wh = 0, rw = 0, rh = 0;
    SDL_GetWindowSize(_window, &ww, &wh);
    if (SDL_GetRendererOutputSize(_renderer, &rw, &rh) != 0) return 1.0f;
    if (ww <= 0 || wh <= 0 || rw <= 0 || rh <= 0) return 1.0f;
    const float sx = (float)rw / (float)ww;
    const float sy = (float)rh / (float)wh;
    const float s = fminf(sx, sy);
    return (s > 1.0f) ? s : 1.0f;
}

static int gui_controls_column_left(void) {
    int left = MAX_X;
    if (live_controls_visible) {
        int x = 0, y = 0, w = 0, h = 0;
        live_controls_layout(&x, &y, &w, &h);
        if (w > 0 && x < left) left = x;
    }
    if (live_scene_visible) {
        int x = 0, y = 0, w = 0, h = 0;
        live_scene_layout(&x, &y, &w, &h);
        if (w > 0 && x < left) left = x;
    }
    return left;
}

static float gui_compute_scene_view_scale(void) {
    const int controls_left = gui_controls_column_left();
    const float hud_min_w = 520.0f;
    const float right_limit = (float)controls_left - hud_min_w - 24.0f;
    const float max_scale_x = right_limit / fmaxf(1.0f, (float)XW2 + 10.0f);
    const float min_plot_h = 300.0f;
    const float max_scale_y = (float)MAX_Y / fmaxf(1.0f, (float)SIM_HEIGHT + min_plot_h);
    float s = fminf(max_scale_x, max_scale_y);
    const float device_s = gui_renderer_device_scale();
    if (s < 1.0f) s = 1.0f;
    if (s < device_s && max_scale_x >= device_s && max_scale_y >= device_s) s = device_s;
    if (s > 2.6f) s = 2.6f;
    return s;
}

static int gui_scene_screen_x(float world_x) {
    return (int)lroundf(world_x * gui_scene_view_scale);
}

static void gui_begin_scene_view(void) {
    if (!_renderer || gui_scene_view_active) return;
    gui_scene_saved_max_x = MAX_X;
    gui_scene_saved_max_y = MAX_Y;
    gui_scene_saved_hist_w = HIST_WIDTH;
    gui_scene_saved_hist_h = HIST_HEIGHT;
    gui_scene_view_scale = gui_compute_scene_view_scale();
    if (gui_scene_view_scale < 1.0f) gui_scene_view_scale = 1.0f;

    if (gui_scene_view_scale > 1.001f) {
        MAX_X = (int)floorf((float)gui_scene_saved_max_x / gui_scene_view_scale);
        MAX_Y = (int)floorf((float)gui_scene_saved_max_y / gui_scene_view_scale);
        if (MAX_X < XW2 + 20) MAX_X = XW2 + 20;
        if (MAX_Y < YW2 + 80) MAX_Y = YW2 + 80;
        gui_recompute_hist_extents();
        SDL_RenderSetScale(_renderer, gui_scene_view_scale, gui_scene_view_scale);
    }
    gui_scene_view_active = 1;
}

static void gui_end_scene_view(void) {
    if (!_renderer || !gui_scene_view_active) return;
    SDL_RenderSetScale(_renderer, 1.0f, 1.0f);
    MAX_X = gui_scene_saved_max_x;
    MAX_Y = gui_scene_saved_max_y;
    HIST_WIDTH = gui_scene_saved_hist_w;
    HIST_HEIGHT = gui_scene_saved_hist_h;
    gui_scene_view_active = 0;
}



//============================================================================
// ##CHRIS: BEGIN PARTICLE-LIFE EXPERIMENT DATA STRUCTURES
//============================================================================
#define MAX_PARTICLE_TYPES 10

static int   pl_enabled = 0;                    // Is particle-life mode active?
static int   pl_num_types = 3;                  // Number of particle types
static int*  pl_particle_types = NULL;          // Array: type of each particle (0 to num_types-1)
static float pl_interaction_matrix[MAX_PARTICLE_TYPES][MAX_PARTICLE_TYPES]; // Interaction matrix [-1,1]
static float pl_particle_radii[MAX_PARTICLE_TYPES]; // Radius for each type (in sigma units)
static float pl_rmax = 200.0f;                  // Maximum interaction distance in PIXELS (0.02 was for normalized [0,1] world!)
static float pl_friction = 0.85f;               // Velocity damping (0=full damp, 1=no damp)
static float pl_force_multiplier = 1.0f;        // Global force scaling

// Type-to-color mapping (RGB, 0-255)
static int pl_type_colors[MAX_PARTICLE_TYPES][3] = {
    {255, 0, 0},      // Type 0: Red
    {0, 255, 0},      // Type 1: Green
    {0, 0, 255},      // Type 2: Blue
    {255, 255, 0},    // Type 3: Yellow
    {255, 0, 255},    // Type 4: Magenta
    {0, 255, 255},    // Type 5: Cyan
    {255, 128, 0},    // Type 6: Orange
    {128, 0, 255},    // Type 7: Purple
    {255, 255, 255},  // Type 8: White
    {128, 128, 128}   // Type 9: Gray
};

// Cyclic boundary condition flag (wrap particles around edges)
static int pl_cyclic_boundaries = 1;  // 1 = wrap, 0 = clamp

//============================================================================
// ##CHRIS: MULTI-DIMENSIONAL RULE TENSORS
//============================================================================
#define MAX_DISTANCE_BINS 10
#define MAX_VELOCITY_BINS 5
#define MAX_ANGLE_BINS 8

// Tensor mode selection
typedef enum {
    TENSOR_MODE_MATRIX = 0,      // Simple 2D matrix M[i][j]
    TENSOR_MODE_DISTANCE = 1,    // 3D: M[i][j][dist]
    TENSOR_MODE_VELOCITY = 2,    // 4D: M[i][j][dist][v_relative]
    TENSOR_MODE_ANGULAR = 3      // 5D: M[i][j][dist][v_rel][angle]
} TensorMode;

static TensorMode pl_tensor_mode = TENSOR_MODE_DISTANCE;

// Distance bins
static int   pl_num_distance_bins = 5;
static float pl_distance_bin_edges[MAX_DISTANCE_BINS + 1];
static float pl_tensor_3d[MAX_PARTICLE_TYPES][MAX_PARTICLE_TYPES][MAX_DISTANCE_BINS];

// Velocity bins (relative velocity magnitude)
static int   pl_num_velocity_bins = 3;  // Slow, medium, fast
static float pl_velocity_bin_edges[MAX_VELOCITY_BINS + 1];
static float pl_tensor_4d[MAX_PARTICLE_TYPES][MAX_PARTICLE_TYPES][MAX_DISTANCE_BINS][MAX_VELOCITY_BINS];

// Angular bins (relative angle between velocity and position vector)
static int   pl_num_angle_bins = 4;  // Front, side, back, other-side (4 quadrants)
static float pl_angle_bin_edges[MAX_ANGLE_BINS + 1];
static float pl_tensor_5d[MAX_PARTICLE_TYPES][MAX_PARTICLE_TYPES][MAX_DISTANCE_BINS][MAX_VELOCITY_BINS][MAX_ANGLE_BINS];

// Helpers
static inline int pl_get_distance_bin(double dist) {
    for (int b = 0; b < pl_num_distance_bins; b++) {
        if (dist < pl_distance_bin_edges[b + 1]) return b;
    }
    return pl_num_distance_bins - 1;
}

static inline int pl_get_velocity_bin(double v_rel) {
    for (int b = 0; b < pl_num_velocity_bins; b++) {
        if (v_rel < pl_velocity_bin_edges[b + 1]) return b;
    }
    return pl_num_velocity_bins - 1;
}

static inline int pl_get_angle_bin(double angle) {
    // Angle in radians [-π, π]
    // Normalize to [0, 2π]
    double a = angle;
    while (a < 0) a += 2.0 * M_PI;
    while (a >= 2.0 * M_PI) a -= 2.0 * M_PI;

    for (int b = 0; b < pl_num_angle_bins; b++) {
        if (a < pl_angle_bin_edges[b + 1]) return b;
    }
    return pl_num_angle_bins - 1;
}
//============================================================================
// ##CHRIS: END MULTI-DIMENSIONAL TENSORS
//============================================================================

//============================================================================
// ##CHRIS: END PARTICLE-LIFE DATA STRUCTURES
//============================================================================





// EDMD state (optional)
static EDMD* g_edmd = NULL;
static int   cli_edmd_acc = 0;          // --edmd-acc: use accelerated EDMD backend (if linked)
static int   g_edmd_is_acc = 0;         // backend used for current g_edmd instance
static int   cli_force_kbt_one = 0; // --kbt1: force K_B*T == 1 (reduced units)
static int   cli_auto_release_wall = 0; // --auto-release: automatically release wall at start
static int   cli_auto_release_after_hold = 0; // --auto-release-after-hold: release automatically once hold-steps reached (GUI + headless)
static int   cli_single_test_mode = 0;  // --single-test: headless mode with main dir output
static int   cli_quiet = 0;             // --quiet: suppress verbose per-particle logs
static int   cli_fixed_dt_set = 0;      // --fixed-dt: override fixed timestep at runtime
static float cli_fixed_dt = 0.0f;
static int   cli_pp_toi = 0;            // --pp-toi: particle-particle time-of-impact correction
static int   cli_pp_toi_max_events = 4; // --pp-toi-max-events: iterations per substep
static int   cli_pp_toi_neighbor_ring = 1; // --pp-toi-neighbor-ring: grid neighbor radius
static int   cli_dist_log = 0;          // --dist-log: enable distribution logging (energy_transfer)
static float cli_dist_dt = 0.0f;        // --dist-dt: cadence in sigma-time (0=every step)
static int   cli_dist_bins = 50;        // --dist-bins: histogram bins
static float cli_dist_vmax = 0.0f;      // --dist-max-speed: vmax for speed histogram (sigma-time units)
static float dist_vmax_runtime = 0.0f;
static char *cli_dist_log_path = NULL;  // --dist-log-path=PATH (supports {seed},{ts})
static FILE *dist_log = NULL;
static int   dist_bins_runtime = 0;
static int  *dist_hist_x = NULL;
static int  *dist_hist_y = NULL;
static int  *dist_hist_v = NULL;
static double dist_box_w_sigma = 0.0;
static double dist_box_h_sigma = 0.0;
static double dist_next_log_abs_time = 0.0;

// EDMD backend dispatch (default vs accelerated). Only valid for the global g_edmd instance.
static inline EDMD* edmd_backend_create(const EDMD_Params* prm) {
    g_edmd_is_acc = cli_edmd_acc ? 1 : 0;
    return g_edmd_is_acc ? edmd_acc_create(prm) : edmd_create(prm);
}
static inline void edmd_backend_destroy(EDMD* S) {
    if (!S) return;
    if (g_edmd_is_acc) edmd_acc_destroy(S);
    else               edmd_destroy(S);
}
static inline double edmd_backend_time(const EDMD* S) {
    return g_edmd_is_acc ? edmd_acc_time(S) : edmd_time(S);
}
static inline const EDMD_Params* edmd_backend_params(const EDMD* S) {
    return g_edmd_is_acc ? edmd_acc_params(S) : edmd_params(S);
}
static inline EDMD_Params* edmd_backend_params_mut(const EDMD* S) {
    return (EDMD_Params*)(g_edmd_is_acc ? edmd_acc_params(S) : edmd_params(S));
}
static inline const EDMD_Particle* edmd_backend_particles(const EDMD* S) {
    return g_edmd_is_acc ? edmd_acc_particles(S) : edmd_particles(S);
}
static inline void edmd_backend_reschedule_all(EDMD* S) {
    if (g_edmd_is_acc) edmd_acc_reschedule_all(S);
    else               edmd_reschedule_all(S);
}
static inline void edmd_backend_reschedule_all_pp_only(EDMD* S) {
    if (g_edmd_is_acc) edmd_acc_reschedule_all_pp_only(S);
    else               edmd_reschedule_all_pp_only(S);
}
static inline void edmd_backend_reset_work(EDMD* S) {
    if (g_edmd_is_acc) edmd_acc_reset_work(S);
    else               edmd_reset_work(S);
}
static inline void edmd_backend_config_dividers(EDMD* S, int count, const double* cx, const double* thickness) {
    if (g_edmd_is_acc) edmd_acc_config_dividers(S, count, cx, thickness);
    else               edmd_config_dividers(S, count, cx, thickness);
}
static inline void edmd_backend_set_divider_motions(EDMD* S, int count, const double* mass, const double* vx) {
    if (g_edmd_is_acc) edmd_acc_set_divider_motions(S, count, mass, vx);
    else               edmd_set_divider_motions(S, count, mass, vx);
}
static inline void edmd_backend_config_pistons(EDMD* S,
                                              int hasL, double xL, double vxL, double mL,
                                              int hasR, double xR, double vxR, double mR)
{
    if (g_edmd_is_acc) edmd_acc_config_pistons(S, hasL, xL, vxL, mL, hasR, xR, vxR, mR);
    else               edmd_config_pistons(S, hasL, xL, vxL, mL, hasR, xR, vxR, mR);
}
static inline double edmd_backend_advance_to(EDMD* S, double t_target) {
    return g_edmd_is_acc ? edmd_acc_advance_to(S, t_target) : edmd_advance_to(S, t_target);
}
static inline double edmd_backend_advance_pp_only_to(EDMD* S, double t_target) {
    return g_edmd_is_acc ? edmd_acc_advance_pp_only_to(S, t_target) : edmd_advance_pp_only_to(S, t_target);
}
static inline void edmd_backend_divider_resolve_overlaps(EDMD* S) {
    if (g_edmd_is_acc) edmd_acc_divider_resolve_overlaps(S);
    else               edmd_divider_resolve_overlaps(S);
}
static inline double edmd_backend_work_divider(const EDMD* S) {
    return g_edmd_is_acc ? edmd_acc_work_divider(S) : edmd_work_divider(S);
}
static inline double edmd_backend_work_divider_i(const EDMD* S, int divider_index) {
    return g_edmd_is_acc ? edmd_acc_work_divider_i(S, divider_index) : edmd_work_divider_i(S, divider_index);
}
static inline double edmd_backend_work_pistonL(const EDMD* S) {
    return g_edmd_is_acc ? edmd_acc_work_pistonL(S) : edmd_work_pistonL(S);
}
static inline double edmd_backend_work_pistonR(const EDMD* S) {
    return g_edmd_is_acc ? edmd_acc_work_pistonR(S) : edmd_work_pistonR(S);
}
static inline double edmd_backend_heat_bath(const EDMD* S) {
    return g_edmd_is_acc ? edmd_acc_heat_bath(S) : edmd_heat_bath(S);
}

// --- EDMD diagnostics/repairs (Szilard-focused) ---
// EDMD assumes an overlap-free initial state. If we introduce overlaps by manual
// "push-out" corrections (e.g. semipermeable divider slabs), EDMD may not fix
// them automatically because there is no well-defined future TOI for already-
// overlapping pairs. For Szilard runs (small N), we do a tiny relaxation pass.
static double edmd_compute_max_pp_overlap(const EDMD* S) {
    if (!S) return 0.0;
    const EDMD_Params* ep = edmd_params(S);
    const EDMD_Particle* P = edmd_particles(S);
    if (!ep || !P || ep->N <= 1) return 0.0;
    const double min_d = 2.0 * ep->radius;
    const double min_d2 = min_d * min_d;
    double max_ov = 0.0;
    for (int i = 0; i < ep->N; ++i) {
        for (int j = i + 1; j < ep->N; ++j) {
            const double dx = P[i].x - P[j].x;
            const double dy = P[i].y - P[j].y;
            const double d2 = dx * dx + dy * dy;
            if (d2 < min_d2) {
                const double d = sqrt(fmax(0.0, d2));
                const double ov = min_d - d;
                if (ov > max_ov) max_ov = ov;
            }
        }
    }
    return max_ov;
}

static void edmd_relax_pp_overlaps_szilard(EDMD* S) {
    if (!S) return;
    if (cli_experiment_preset != EXPERIMENT_PRESET_SZILARD_ENGINE) return;
    if (sim_mode != MODE_EDMD) return;

    const EDMD_Params* ep = edmd_params(S);
    if (!ep || ep->N <= 1) return;
    // Keep this tiny and safe: Szilard runs are small; avoid O(N^2) on large-N experiments.
    if (ep->N > 256) return;

    const double max_ov0 = edmd_compute_max_pp_overlap(S);
    if (!(max_ov0 > 1e-7)) return;
    if (!cli_quiet) {
        printf("[EDMD] Szilard: detected PP overlap (max=%.3g px). Relaxing...\n", max_ov0);
    }

    EDMD_Particle* P = (EDMD_Particle*)edmd_particles(S);
    if (!P) return;

    const double boxW = ep->boxW;
    const double boxH = ep->boxH;
    const double R = ep->radius;
    const double min_d = 2.0 * R;
    const double min_d2 = min_d * min_d;
    const double eps = 1e-6; // (pixels) small separation buffer
    const int gate_d = 0;
    const int enforce_sep_gate =
        (ep->species &&
         ep->divider_count > gate_d &&
         ep->divider_thickness[gate_d] > 0.0 &&
         (ep->divider_gate_mode[gate_d] == 1 || ep->divider_gate_mode[gate_d] == 2));
    const double gate_x = enforce_sep_gate ? ep->divider_x[gate_d] : 0.0;
    const double gate_margin = 0.5 * (enforce_sep_gate ? ep->divider_thickness[gate_d] : 0.0) + R + eps;
    const double gate_left_max = gate_x - gate_margin;
    const double gate_right_min = gate_x + gate_margin;

    for (int it = 0; it < 16; ++it) {
        int any = 0;
        for (int i = 0; i < ep->N; ++i) {
            for (int j = i + 1; j < ep->N; ++j) {
                double dx = P[i].x - P[j].x;
                double dy = P[i].y - P[j].y;
                double d2 = dx * dx + dy * dy;
                if (!(d2 < min_d2)) continue;

                double d = sqrt(fmax(0.0, d2));
                // If nearly identical centers, choose an arbitrary separation axis.
                if (d < 1e-12) { dx = 1.0; dy = 0.0; d = 1.0; }
                const double inv_d = 1.0 / d;
                const double nx = dx * inv_d;
                const double ny = dy * inv_d;
                const double ov = (min_d - d);
                const double push = 0.5 * (ov + eps);

                P[i].x += nx * push; P[i].y += ny * push;
                P[j].x -= nx * push; P[j].y -= ny * push;

                // Clamp to box (EDMD outer walls are hard at [R, box-R])
                if (P[i].x < R) P[i].x = R; if (P[i].x > boxW - R) P[i].x = boxW - R;
                if (P[i].y < R) P[i].y = R; if (P[i].y > boxH - R) P[i].y = boxH - R;
                if (P[j].x < R) P[j].x = R; if (P[j].x > boxW - R) P[j].x = boxW - R;
                if (P[j].y < R) P[j].y = R; if (P[j].y > boxH - R) P[j].y = boxH - R;

                if (enforce_sep_gate) {
                    const int spi = (ep->species[i] != 0) ? 1 : 0;
                    const int spj = (ep->species[j] != 0) ? 1 : 0;
                    if (spi == 0) {
                        if (P[i].x > gate_left_max) P[i].x = gate_left_max;
                    } else {
                        if (P[i].x < gate_right_min) P[i].x = gate_right_min;
                    }
                    if (spj == 0) {
                        if (P[j].x > gate_left_max) P[j].x = gate_left_max;
                    } else {
                        if (P[j].x < gate_right_min) P[j].x = gate_right_min;
                    }
                }

                P[i].coll_count++;
                P[j].coll_count++;
                any = 1;
            }
        }
        if (!any) break;
    }
}

static void edmd_enforce_szilard_gate_sides(EDMD* S) {
    if (!S) return;
    if (cli_experiment_preset != EXPERIMENT_PRESET_SZILARD_ENGINE) return;
    if (sim_mode != MODE_EDMD) return;

    const EDMD_Params* ep = edmd_params(S);
    EDMD_Particle* P = (EDMD_Particle*)edmd_particles(S);
    if (!ep || !P || !ep->species || ep->divider_count <= 0) return;

    const int d = 0;
    if (ep->divider_thickness[d] <= 0.0) return;
    if (!(ep->divider_gate_mode[d] == 1 || ep->divider_gate_mode[d] == 2)) return;

    const double gate_x = ep->divider_x[d];
    const double R = ep->radius;
    const double eps = 1e-6;
    const double margin = 0.5 * ep->divider_thickness[d] + R + eps;
    const double left_ok = gate_x - margin;
    const double right_ok = gate_x + margin;
    const double repair_band = 8.0 * R;
    int changed = 0;

    for (int i = 0; i < ep->N; ++i) {
        const int sp = (ep->species[i] != 0) ? 1 : 0;
        int target_side = -1;
        if (ep->divider_gate_mode[d] == 1) {
            if (sp == ep->divider_gate_species[d]) {
                target_side = ep->divider_gate_target_side[d] ? 1 : 0;
            } else {
                target_side = 1 - (ep->divider_gate_target_side[d] ? 1 : 0);
            }
        } else if (ep->divider_gate_mode[d] == 2) {
            if (sp == 0) target_side = 0;
            else if (sp == 1) target_side = 1;
        }
        if (target_side < 0) continue;

        if (target_side == 0) {
            if (P[i].x > left_ok && P[i].x < gate_x + repair_band) {
                P[i].x = left_ok;
                P[i].coll_count++;
                changed = 1;
            }
        } else {
            if (P[i].x < right_ok && P[i].x > gate_x - repair_band) {
                P[i].x = right_ok;
                P[i].coll_count++;
                changed = 1;
            }
        }
    }

    if (changed) {
        edmd_relax_pp_overlaps_szilard(S);
        edmd_reschedule_all(S);
    }
}

// Convenience: keep the rest of the code using the canonical `edmd_*` names.
// These macros map calls to the selected backend (default vs accelerated).
#define edmd_create                  edmd_backend_create
#define edmd_destroy                 edmd_backend_destroy
#define edmd_time                    edmd_backend_time
#define edmd_params                  edmd_backend_params
#define edmd_particles               edmd_backend_particles
#define edmd_reschedule_all          edmd_backend_reschedule_all
#define edmd_reschedule_all_pp_only  edmd_backend_reschedule_all_pp_only
#define edmd_reset_work              edmd_backend_reset_work
#define edmd_config_dividers         edmd_backend_config_dividers
#define edmd_set_divider_motions     edmd_backend_set_divider_motions
#define edmd_config_pistons          edmd_backend_config_pistons
#define edmd_advance_to              edmd_backend_advance_to
#define edmd_advance_pp_only_to      edmd_backend_advance_pp_only_to
#define edmd_divider_resolve_overlaps edmd_backend_divider_resolve_overlaps
#define edmd_work_divider            edmd_backend_work_divider
#define edmd_work_divider_i          edmd_backend_work_divider_i
#define edmd_work_pistonL            edmd_backend_work_pistonL
#define edmd_work_pistonR            edmd_backend_work_pistonR
#define edmd_heat_bath               edmd_backend_heat_bath

// Windowed "efficiency/gain" metric for energy_transfer: measure peak spring response in a fixed
// window after piston motion ends, baseline-subtracted at piston stop.
static int    cli_eff_window_set = 0;           // --eff-window=T
static float  cli_eff_window_sigma = 20.0f;     // sigma-time window length after piston stop
static int    cli_eff_stop_after_window = 0;    // --eff-stop-after-window: end energy_transfer run after window completes

typedef enum {
    EFF_OUTPUT_SPRING = 0,   // spring potential energy (legacy; applies a spring force in TIME/EDMD via impulses)
    EFF_OUTPUT_WALL_KE = 1   // kinetic energy of the primary (leftmost) internal wall (pure hard-body EDMD compatible)
} EffOutputMode;
static int          cli_eff_output_set = 0;     // --eff-output=...
static EffOutputMode cli_eff_output_mode = EFF_OUTPUT_SPRING;

// Szilard engine experiment (particle separation)
static int    cli_szilard_equil_steps = 10000;
/* Szilard interactive: separation/recombine duration (sigma-time).
   Make the default a bit slower to reduce jamming/teleport-like artifacts during G/H. */
static double cli_szilard_sep_duration = 6.0;   // sigma-time
static int    cli_szilard_perm_swap = 1;
static float  cli_szilard_species_ratio = 0.5f; // fraction of species 0
static int    cli_szilard_n = -1;              // optional override
static int    cli_szilard_unstick = 1;          // Szilard: nudge pinned particles off outer walls if they pin with Vx≈0
static int    cli_szilard_unstick_hits = 200;   // consecutive boundary hits before applying nudge
static float  cli_szilard_unstick_vx_frac = 0.05f; // fraction of speed to redirect into x
static char  *cli_szilard_log_path = NULL;      // --szilard-log=PATH
static double cli_szilard_log_dt = 0.2;         // --szilard-log-dt=T (sigma-time), 0=every step
// Simple prediction box experiment
static int    cli_simple_box_cycles = 40;
static int    cli_simple_box_equil_steps = 4000;
static int    cli_simple_box_relax_steps = 800;
static double cli_simple_box_log_dt = 0.25;
static double cli_simple_box_target_min_sigma = -1.0;
static double cli_simple_box_target_max_sigma = -1.0;
static double cli_simple_box_initial_target_sigma = -1.0;
static double cli_simple_box_corr_alpha = 0.85;
static double cli_simple_box_noise_sigma = 2.0;
static int    cli_simple_box_period = 8;
static int    cli_simple_box_gui_delay_ms = 8;
static int    cli_simple_box_render_every = 1;
static char  *cli_simple_box_run_dir = NULL;    // --simple-box-run-dir=PATH
static char  *cli_simple_box_run_label = NULL;  // --simple-box-run-label=NAME
static char   cli_simple_box_target_protocol[32] = "correlated";
static int   *sz_species = NULL;
static int   *sz_memory = NULL;
static int   *sz_gate_labels = NULL;
static int    sz_perm_active = 0;
static int    sz_perm_species = 0;
static int    sz_perm_use_memory = 0;
static int    sz_perm_wall_index = 0;
static int    sz_perm_target_side = 0; // 0=left, 1=right
static int    sz_perm_dual_active = 0; // 1=both species gated at L0 (species0->left, species1->right)
static int    sz_fixed_gate_active = 0;
static int    sz_fixed_gate_species = 0;   // species blocked by fixed L0 gate
static int    sz_fixed_gate_block_dir = 0; // 0=L->R blocked, 1=R->L blocked
static int    sz_fixed_gate_use_memory = 0;
static int    sz_interactive_enabled = 0;
static int    sz_interactive_phase = 0;
static float  sz_interactive_home_px = 0.0f;
static float  sz_interactive_target_px = 0.0f;
static float  sz_interactive_vx_px = 0.0f;
static float  sz_interactive_piston_home_px = 0.0f;
static float  sz_interactive_piston_target_px = 0.0f;
static char   g_simple_box_dir[512] = "experiments_simple_box";
static char   g_szilard_dir[512] = "experiments_szilard";
static char   g_szilard_run_dir[1024] = {0}; /* per-run timestamped output dir (Szilard interactive) */
static int    g_szilard_run_active = 0;
static int    sz_unstick_hit_count[NUM_PARTICLES] = {0};
static int    sz_branch_restore_requested = 0;
static long long sz_log_index = 0;
static int    sz_log_branch_id = 0;
static char   sz_log_branch_label[32] = "Eq";
static int    sz_log_reset_requested = 0;
static double sz_sep_perfect_since = -1.0;
static double sz_sep_perfect_hold_time = 10.0;
static int    sz_xflip_enabled = 0;
static double sz_xflip_rate = 0.05;          /* sigma-time^-1; Stage-1 baseline: unbiased x-relaxation over ~10-50 sigma-time */
static double sz_xflip_tail_sigma = 1.5;     /* only upper-tail particles can flip x */
static double sz_gate_hot_speed_mult = 3.0;  /* hot particles above this * mean speed can leak gates */
static double sz_xeq_ready_since = -1.0;
static double sz_xeq_ready_hold_time = 1.0;  /* sigma-time */
static double sz_xeq_info_threshold_bits = 0.05;
static int    sz_branch_ready_min_xy_mismatch = 5; /* allow a partially decorrelated branch before full Xeq */
static int    sz_branch_ready_reason = 0; /* 0=unset, 1=full Xeq info threshold, 2=xy mismatch threshold */
static double sz_memory_erase_heat_min = 0.0;
static int    sz_partial_obs_enabled = 0;
static double sz_partial_obs_ref_imem_gate_bits = 1.0;
static double sz_partial_obs_peak_uncertainty_frac = 0.0;
static int    sz_pending_recombine_mode = 0; /* 0=none, 1=H, 2=J, 3=O */
static int    sz_stage1_baseline_mode = 1;   /* Stage-1 clean baseline: visible x equilibrates, hidden y stays fixed */
static int    sz_load_enabled = 1;           /* Stage-2 bookkeeping-only spring load during Xeq */
static double sz_load_eta = 1.0;             /* visible free-energy drop -> load efficiency */
static double sz_load_k = 1.0;               /* abstract spring constant */
static double sz_load_fx_start = 0.0;        /* visible free energy at start of Xeq */
static double sz_load_fx_drop = 0.0;         /* max(0, F_x(start_Xeq) - F_x(now)) */
static double sz_load_z = 0.0;               /* spring extension for bookkeeping only */
static double sz_load_energy = 0.0;          /* 0.5*k*z^2 */
static double sz_load_work = 0.0;            /* extracted work into the abstract load */
static double sz_load_qdiss = 0.0;           /* visible free-energy drop not captured by the load */

typedef enum {
    SZ_VIEW_SPECIES = 0,
    SZ_VIEW_MEMORY = 1,
    SZ_VIEW_OVERLAY = 2
} szilard_view_mode_t;
static szilard_view_mode_t sz_view_mode = SZ_VIEW_SPECIES;

typedef struct {
    int valid;
    int particles_active;
    int steps_elapsed;
    int simulation_started;
    int sz_interactive_phase;
    int wall_enabled;
    int wall_hold_enabled;
    int wall_is_released;
    int wall_hold_steps;
    int sz_perm_active;
    int sz_perm_species;
    int sz_perm_use_memory;
    int sz_perm_wall_index;
    int sz_perm_target_side;
    int sz_perm_dual_active;
    int sz_fixed_gate_active;
    int sz_fixed_gate_species;
    int sz_fixed_gate_block_dir;
    int sz_fixed_gate_use_memory;
    int sz_xflip_enabled;
    int sz_partial_obs_enabled;
    int sz_branch_ready_reason;
    int sz_load_enabled;
    double simulation_time;
    double wall_release_time;
    double sz_xeq_ready_since;
    double sz_memory_erase_heat_min;
    double sz_load_eta;
    double sz_load_k;
    double sz_load_fx_start;
    double sz_load_fx_drop;
    double sz_load_z;
    double sz_load_energy;
    double sz_load_work;
    double sz_load_qdiss;
    float wall_x;
    float wall_x_old;
    float vx_wall;
    float piston_left_x;
    float piston_right_x;
    float vx_piston_left;
    float vx_piston_right;
    double X[NUM_PARTICLES];
    double Y[NUM_PARTICLES];
    double Vx[NUM_PARTICLES];
    double Vy[NUM_PARTICLES];
    int sz_species[NUM_PARTICLES];
    int sz_memory[NUM_PARTICLES];
} SzilardBranchSnapshot;

static SzilardBranchSnapshot g_sz_branch_snapshot = {0};

static bool dist_path_has_placeholder(const char *s) {
    if (!s) return false;
    return (strstr(s, "{seed}") != NULL) || (strstr(s, "{ts}") != NULL);
}

static void build_dist_log_path(char *out, size_t out_sz) {
    const char *base = (cli_dist_log_path && cli_dist_log_path[0]) ? cli_dist_log_path
                                                                    : "energy_transfer_distributions.csv";
    if (cli_dist_log_path && !dist_path_has_placeholder(base)) {
        const char *dot = strrchr(base, '.');
        if (dot && strcmp(dot, ".csv") == 0) {
            int prefix_len = (int)(dot - base);
            if (prefix_len < 0) prefix_len = 0;
            snprintf(out, out_sz, "%.*s_seed%u%s", prefix_len, base, cli_seed, dot);
        } else {
            snprintf(out, out_sz, "%s_seed%u", base, cli_seed);
        }
        return;
    }

    // Replace placeholders {seed} and {ts}
    char buf[1024];
    size_t bi = 0;
    const char *p = base;
    const unsigned int seed_val = cli_seed;
    const unsigned int ts_val = (unsigned int)time(NULL);
    while (*p && bi + 1 < sizeof(buf)) {
        if (strncmp(p, "{seed}", 6) == 0) {
            bi += (size_t)snprintf(buf + bi, sizeof(buf) - bi, "%u", seed_val);
            p += 6;
        } else if (strncmp(p, "{ts}", 4) == 0) {
            bi += (size_t)snprintf(buf + bi, sizeof(buf) - bi, "%u", ts_val);
            p += 4;
        } else {
            buf[bi++] = *p++;
        }
    }
    buf[bi] = '\0';
    snprintf(out, out_sz, "%s", buf);
}

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

// Write current run parameters for Python post-processing (wall_x_FFT.py)
static void write_run_params_json_to(const char *path) {
    FILE *f = fopen(path, "w");
    if (!f) return;
    const char* mode_str = (sim_mode==MODE_TIME?"time":(sim_mode==MODE_RK4?"rk4":(sim_mode==MODE_EDMD?"edmd":"edmd-hybrid")));
    // Derive wall mass factor in particle units
    double wall_mass_factor_out = (double)WALL_MASS / (double)PARTICLE_MASS;
    // If per-segment temperatures were provided, write as CSV list; else single
    fprintf(f, "{\n");
    fprintf(f, "  \"mode\": \"%s\",\n", mode_str);
    fprintf(f, "  \"particles\": %d,\n", (cli_override_particles > 0 ? cli_override_particles : particles_active));
    fprintf(f, "  \"L0\": %.6f,\n", (double)L0_UNITS);
    fprintf(f, "  \"height\": %.6f,\n", (double)HEIGHT_UNITS);
    fprintf(f, "  \"num_walls\": %d,\n", (int)num_internal_walls);
    fprintf(f, "  \"wall_mass_factor\": %.6f,\n", wall_mass_factor_out);
    // wall thicknesses in sigma units
    double wall_thick_sigma = (double)(WALL_THICKNESS / PIXELS_PER_SIGMA);
    double wall_thick_vis_sigma = (wall_thickness_visual_runtime > 0.0f)
                                  ? (double)(wall_thickness_visual_runtime / PIXELS_PER_SIGMA)
                                  : wall_thick_sigma; // falls back to physics thickness
    fprintf(f, "  \"wall_thickness\": %.6f,\n", wall_thick_sigma);
    fprintf(f, "  \"wall_thickness_vis\": %.6f,\n", wall_thick_vis_sigma);
    fprintf(f, "  \"temperature\": %.9g,\n", (double)temperature_runtime);
    // optional temperature segments array
    if (cli_temperature_segments && cli_temperature_segments_count > 0) {
        fprintf(f, "  \"temperature_segments\": [");
        for (size_t i=0;i<cli_temperature_segments_count;i++){
            if(i) fprintf(f, ", ");
            fprintf(f, "%.9g", (double)cli_temperature_segments[i]);
        }
        fprintf(f, "],\n");
    }
    fprintf(f, "  \"kbt1\": %d,\n", (int)cli_force_kbt_one);
    fprintf(f, "  \"kB_effective\": %.9g,\n", (double)kB_effective());
    fprintf(f, "  \"kBT\": %.9g,\n", (double)kBT_effective());
    fprintf(f, "  \"dt\": %.9g,\n", (double)fixed_dt_runtime);
    fprintf(f, "  \"output_dt\": %.9g,\n", (double)cli_output_dt);
    fprintf(f, "  \"time_span\": %.9g,\n", (double)cli_time_span);
    fprintf(f, "  \"units\": \"%s\",\n", (MOLECULAR_MODE?"reduced":"real"));
    // Export actual particle radius in sigma units (accounts for --eta, --particle-radius, etc.)
    float actual_radius_sigma = (float)(PARTICLE_RADIUS / PIXELS_PER_SIGMA);
    fprintf(f, "  \"radius_sigma\": %.6f,\n", (double)actual_radius_sigma);
    fprintf(f, "  \"pixels_per_sigma\": %.6f\n", (double)PIXELS_PER_SIGMA);
    fprintf(f, "}\n");
    fclose(f);
}
static void write_run_params_json(void) {
    write_run_params_json_to("run_params.json");
}

// Recompute runtime dt to preserve v·dt/σ when temperature changes (speed_of_sound experiment)
static inline void update_dt_runtime_for_temperature(void) {
    // Default: keep compile-time dt
    fixed_dt_runtime = FIXED_DT;
    // If user explicitly overrides dt, keep it fixed regardless of temperature
    if (cli_fixed_dt_set && cli_fixed_dt > 0.0f) {
        fixed_dt_runtime = cli_fixed_dt;
        return;
    }
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

    // Velocity distribution parameters per segment (mean/variance); used for HUD diagnostics.
    double mean_vx[MAX_DIVIDER_CAPACITY + 1];
    double mean_vy[MAX_DIVIDER_CAPACITY + 1];
    double mean_sp[MAX_DIVIDER_CAPACITY + 1];
    double m2_vx[MAX_DIVIDER_CAPACITY + 1];
    double m2_vy[MAX_DIVIDER_CAPACITY + 1];
    double m2_sp[MAX_DIVIDER_CAPACITY + 1];
    for (int s = 0; s < segment_count; ++s) {
        mean_vx[s] = mean_vy[s] = mean_sp[s] = 0.0;
        m2_vx[s] = m2_vy[s] = m2_sp[s] = 0.0;
    }

    if (segment_vx_mean) memset(segment_vx_mean, 0, (size_t)segment_count * sizeof(double));
    if (segment_vx_var) memset(segment_vx_var, 0, (size_t)segment_count * sizeof(double));
    if (segment_vy_mean) memset(segment_vy_mean, 0, (size_t)segment_count * sizeof(double));
    if (segment_vy_var) memset(segment_vy_var, 0, (size_t)segment_count * sizeof(double));
    if (segment_speed_mean) memset(segment_speed_mean, 0, (size_t)segment_count * sizeof(double));
    if (segment_speed_var) memset(segment_speed_var, 0, (size_t)segment_count * sizeof(double));

    for (int i = 0; i < n; ++i) {
        int seg = segment_index_for_position(X[i]);
        if (seg < 0) seg = 0; if (seg >= segment_count) seg = segment_count - 1;
        segment_counts[seg]++;
        if (segment_ke) {
            double ke = 0.5 * PARTICLE_MASS * (Vx[i] * Vx[i] + Vy[i] * Vy[i]);
            segment_ke[seg] += ke;
        }

        // Welford updates for vx/vy/speed
        const int nn = segment_counts[seg];
        const double vx = (double)Vx[i];
        const double vy = (double)Vy[i];
        const double sp = sqrt(vx * vx + vy * vy);

        double dx = vx - mean_vx[seg];
        mean_vx[seg] += dx / (double)nn;
        m2_vx[seg] += dx * (vx - mean_vx[seg]);

        double dy = vy - mean_vy[seg];
        mean_vy[seg] += dy / (double)nn;
        m2_vy[seg] += dy * (vy - mean_vy[seg]);

        double ds = sp - mean_sp[seg];
        mean_sp[seg] += ds / (double)nn;
        m2_sp[seg] += ds * (sp - mean_sp[seg]);
    }
    if (segment_temperature && segment_ke) {
        for (int seg = 0; seg < segment_count; ++seg) {
            if (segment_counts[seg] > 0) {
                // Use kB_effective so reported T is consistent when --kbt1 is enabled
                segment_temperature[seg] = segment_ke[seg] / (segment_counts[seg] * kB_effective());
            } else {
                segment_temperature[seg] = 0.0;
            }
        }
    }

    if (segment_vx_mean && segment_vx_var && segment_vy_mean && segment_vy_var && segment_speed_mean && segment_speed_var) {
        for (int seg = 0; seg < segment_count; ++seg) {
            const int nn = segment_counts[seg];
            segment_vx_mean[seg] = mean_vx[seg];
            segment_vy_mean[seg] = mean_vy[seg];
            segment_speed_mean[seg] = mean_sp[seg];
            if (nn >= 2) {
                segment_vx_var[seg] = m2_vx[seg] / (double)(nn - 1);
                segment_vy_var[seg] = m2_vy[seg] / (double)(nn - 1);
                segment_speed_var[seg] = m2_sp[seg] / (double)(nn - 1);
            } else {
                segment_vx_var[seg] = 0.0;
                segment_vy_var[seg] = 0.0;
                segment_speed_var[seg] = 0.0;
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
        // Use kB_effective so that with --kbt1 the product kB*T remains invariant
        double T_seg = ke_seg / (c * kB_effective());
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
    bool  eff_enabled;              // Enable efficiency/gain measurement (spring or wall-KE)
    bool spring_enabled;            // Enable spring system for energy measurement
    EffOutputMode eff_output_mode;  // What SpringE_* fields represent (spring or wall KE)
    // Windowed "gain" metrics after piston stops (baseline-subtracted peak in [t_stop, t_stop+Tw])
    bool  eff_window_pending_start; // set when piston reaches target; start window on next energy update
    bool  eff_window_started;
    bool  eff_window_active;
    bool  eff_window_done;
    float eff_piston_stop_time;     // absolute sigma-time when piston motion ended
    float spring_energy_at_stop;    // SpringE at piston stop (baseline)
    float spring_energy_peak_window;// Peak SpringE observed during window
    float spring_energy_peak_window_time; // absolute sigma-time when SpringE_peak_window occurred
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
    .eff_enabled = false,
    .spring_enabled = false,
    .eff_output_mode = EFF_OUTPUT_SPRING,
    .eff_window_pending_start = false,
    .eff_window_started = false,
    .eff_window_active = false,
    .eff_window_done = false,
    .eff_piston_stop_time = 0.0f,
    .spring_energy_at_stop = 0.0f,
    .spring_energy_peak_window = 0.0f,
    .spring_energy_peak_window_time = 0.0f
};

// Define after energy_measurement is declared
static inline void start_piston_push_tracking(void) {
    piston_push_tracking  = true;
    piston_work_baseline  = piston_work_left + piston_work_right;
    piston_work_delta_max = 0.0;
    // Reset windowed gain metric (energy_transfer "efficiency" proxy).
    energy_measurement.eff_window_pending_start = false;
    energy_measurement.eff_window_started = false;
    energy_measurement.eff_window_active = false;
    energy_measurement.eff_window_done = false;
    energy_measurement.eff_piston_stop_time = 0.0f;
    energy_measurement.spring_energy_at_stop = 0.0f;
    energy_measurement.spring_energy_peak_window = 0.0f;
    energy_measurement.spring_energy_peak_window_time = 0.0f;
    // Only initialize energy measurement if explicitly requested or in energy_transfer preset
    if (!energy_measurement.eff_enabled && num_internal_walls > 0) {
        if (cli_enable_energy_measurement || cli_experiment_preset == EXPERIMENT_PRESET_ENERGY_TRANSFER) {
            initialize_energy_measurement();
        }
    }
}

static inline void stop_right_piston_due_to_work_target(void) {
    vx_piston_right = 0.0f;
    piston_step_active = false;
    piston_step_direction = -1.0f;
    piston_protocol_start_time = -1.0f;
    // Start windowed response measurement on the next energy update.
    if (energy_measurement.eff_enabled) {
        energy_measurement.eff_window_pending_start = true;
        energy_measurement.eff_window_started = false;
        energy_measurement.eff_window_active = false;
        energy_measurement.eff_window_done = false;
    }
}

static inline void maybe_stop_right_piston_on_work_target(void) {
    if (!cli_piston_work_target_set) return;
    if (!piston_step_active) return;
    if (!piston_push_tracking) return;
    if (!(cli_piston_work_target > 0.0)) return;
    // Use max observed input work since push start to avoid transient sign flips.
    if (piston_work_delta_max + 1e-12 >= cli_piston_work_target) {
        stop_right_piston_due_to_work_target();
        if (!cli_quiet) {
            printf("⏹️  Piston stopped on work target: W_in=%.6g (target=%.6g)\n",
                   piston_work_delta_max, cli_piston_work_target);
        }
    }
}
// Global flag to pause the simulation
bool paused = false;





static void reset_simulation_state(void) {
    if (cli_experiment_preset == EXPERIMENT_PRESET_SZILARD_ENGINE) {
        szilard_prepare_interactive_config();
    }
    rebuild_internal_walls();
    initialize_simulation_parameters();
    initialize_simulation();
    if (cli_experiment_preset == EXPERIMENT_PRESET_SZILARD_ENGINE) {
        szilard_setup_interactive_state();
    }
    memset(sz_unstick_hit_count, 0, sizeof(sz_unstick_hit_count));
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
    piston_step_direction = -1.0f;
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

static void csv_put_escaped(FILE *f, const char *s) {
    fputc('"', f);
    if (s) {
        for (const char *p = s; *p; ++p) {
            if (*p == '"') fputc('"', f);
            fputc(*p, f);
        }
    }
    fputc('"', f);
}

static void mkdir_p(const char *path) {
    if (!path || !*path) return;
    char tmp[512];
    snprintf(tmp, sizeof(tmp), "%s", path);
    size_t len = strlen(tmp);
    if (len == 0) return;
    if (tmp[len - 1] == '/') tmp[len - 1] = '\0';
    for (char *p = tmp + 1; *p; ++p) {
        if (*p == '/') {
            *p = '\0';
            if (mkdir(tmp, 0755) != 0 && errno != EEXIST) {
                *p = '/';
                return;
            }
            *p = '/';
        }
    }
    (void)mkdir(tmp, 0755);
}

static void mkdir_p_for_file(const char *file_path) {
    if (!file_path || !*file_path) return;
    char tmp[512];
    snprintf(tmp, sizeof(tmp), "%s", file_path);
    char *slash = strrchr(tmp, '/');
    if (!slash) return;
    *slash = '\0';
    if (tmp[0] == '\0') return;
    mkdir_p(tmp);
}

static bool shell_arg_needs_quotes(const char *s) {
    if (!s || !*s) return true;
    for (const unsigned char *p = (const unsigned char *)s; *p; ++p) {
        unsigned char c = *p;
        if (isalnum(c)) continue;
        switch (c) {
            case '_': case '-': case '.': case '/': case ':': case '=': case ',': case '+': case '@':
                continue;
            default:
                return true;
        }
    }
    return false;
}

static void fprint_shell_escaped(FILE *f, const char *s) {
    if (!s || !*s) {
        fputs("''", f);
        return;
    }
    if (!shell_arg_needs_quotes(s)) {
        fputs(s, f);
        return;
    }
    fputc('\'', f);
    for (const char *p = s; *p; ++p) {
        if (*p == '\'') {
            fputs("'\"'\"'", f);
        } else {
            fputc(*p, f);
        }
    }
    fputc('\'', f);
}

static void write_command_md_to_dir(const char *dir_path, int argc, char **argv) {
    if (!dir_path || !*dir_path) return;

    char out_path[1024];
    snprintf(out_path, sizeof(out_path), "%s/00_COMMAND.md", dir_path);

    FILE *f = fopen(out_path, "w");
    if (!f) {
        if (!cli_quiet) {
            fprintf(stderr, "⚠️ Could not write %s: %s\n", out_path, strerror(errno));
        }
        return;
    }

    time_t now = time(NULL);
    struct tm *tm_now = localtime(&now);
    char ts[64] = {0};
    if (tm_now) {
        strftime(ts, sizeof(ts), "%Y-%m-%d %H:%M:%S %Z", tm_now);
    }

    char cwd[1024] = {0};
    if (!getcwd(cwd, sizeof(cwd))) {
        snprintf(cwd, sizeof(cwd), "(getcwd failed: %s)", strerror(errno));
    }

    fprintf(f, "# Command\n\n");
    fprintf(f, "- Timestamp: %s\n", ts[0] ? ts : "(unknown)");
    fprintf(f, "- CWD: `%s`\n\n", cwd);
    fprintf(f, "```sh\n");
    for (int i = 0; i < argc; ++i) {
        if (i) fputc(' ', f);
        fprint_shell_escaped(f, argv[i]);
    }
    fprintf(f, "\n```\n");
    fclose(f);
}

static void write_speed_of_sound_plot_commands_to_dir(const char *dir_path) {
    if (!dir_path || !*dir_path) return;

    char plot_path[1200];
    snprintf(plot_path, sizeof(plot_path), "%s/01_PLOT_COMMANDS.md", dir_path);
    FILE *f = fopen(plot_path, "w");
    if (!f) {
        if (!cli_quiet) {
            fprintf(stderr, "⚠️ Could not write %s: %s\n", plot_path, strerror(errno));
        }
        return;
    }

    fprintf(f, "# Plot commands\n\n");
    fprintf(f, "Run from the repo root:\n\n");
    fprintf(f, "```sh\n");
    fprintf(f, "python3 hspist3/analyze_speed_of_sound_by_eta.py \\\n");
    fprintf(f, "  --dir hspist3/%s \\\n", dir_path);
    fprintf(f, "  --write-final \\\n");
    fprintf(f, "  --roman-ref \\\n");
    fprintf(f, "  --theory-cs simple\n");
    fprintf(f, "```\n\n");
    fprintf(f, "Combine the newest two separated eta simulation folders:\n\n");
    fprintf(f, "```sh\n");
    fprintf(f, "python3 hspist3/analyze_speed_of_sound_by_eta.py \\\n");
    fprintf(f, "  --separated-eta-simulations=true \\\n");
    fprintf(f, "  --latest-count 2 \\\n");
    fprintf(f, "  --write-final \\\n");
    fprintf(f, "  --roman-ref \\\n");
    fprintf(f, "  --theory-cs simple\n");
    fprintf(f, "```\n\n");
    fprintf(f, "Combined output:\n\n");
    fprintf(f, "```text\n");
    fprintf(f, "hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/combined_latest_2_analysis/final_plots/FINAL speed_of_sound_on_packing_fracture.pdf\n");
    fprintf(f, "```\n\n");
    fprintf(f, "Run the full high/mid/low eta split workflow and plot automatically:\n\n");
    fprintf(f, "```sh\n");
    fprintf(f, "python3 hspist3/run_speed_of_sound_eta_split_workflow.py\n");
    fprintf(f, "```\n\n");
    fprintf(f, "Main output:\n\n");
    fprintf(f, "```text\n");
    fprintf(f, "hspist3/%s/analysis_by_eta/final_plots/FINAL speed_of_sound_on_packing_fracture.pdf\n", dir_path);
    fprintf(f, "hspist3/%s/analysis_by_eta/final_plots/combined_speed_of_sound_summary.csv\n", dir_path);
    fprintf(f, "hspist3/%s/analysis_by_eta/final_plots/speed_of_sound_freq_vs_kterm_combined.pdf\n", dir_path);
    fprintf(f, "```\n\n");
    fprintf(f, "Notes:\n");
    fprintf(f, "- The analyzer stages each `L0`/eta group separately and uses eta-specific FFT windows.\n");
    fprintf(f, "- This is preferred over one global FFT window for all eta values.\n");
    fprintf(f, "- If the raw CSV `eta` column is historically wrong, the analyzer recomputes eta from `N`, `r`, `L0`, and `H`.\n");
    fclose(f);

    char combined_path[1200];
    snprintf(combined_path, sizeof(combined_path), "%s/00_COMMAND_AND_PLOT_COMMAND.md", dir_path);
    FILE *g = fopen(combined_path, "w");
    if (!g) return;

    fprintf(g, "# Command and Plot Commands\n\n");
    fprintf(g, "## Simulation Command\n\n");
    fprintf(g, "```sh\n");
    for (int i = 0; i < g_main_argc; ++i) {
        if (i) fputc(' ', g);
        fprint_shell_escaped(g, g_main_argv[i]);
    }
    fprintf(g, "\n```\n\n");
    fprintf(g, "## Plot Command\n\n");
    fprintf(g, "Run from the repo root:\n\n");
    fprintf(g, "```sh\n");
    fprintf(g, "python3 hspist3/analyze_speed_of_sound_by_eta.py \\\n");
    fprintf(g, "  --dir hspist3/%s \\\n", dir_path);
    fprintf(g, "  --write-final \\\n");
    fprintf(g, "  --roman-ref \\\n");
    fprintf(g, "  --theory-cs simple\n");
    fprintf(g, "```\n\n");
    fprintf(g, "## Combine Newest Two Eta Runs\n\n");
    fprintf(g, "Use this after running separate dense/high-eta and low-eta simulation commands:\n\n");
    fprintf(g, "```sh\n");
    fprintf(g, "python3 hspist3/analyze_speed_of_sound_by_eta.py \\\n");
    fprintf(g, "  --separated-eta-simulations=true \\\n");
    fprintf(g, "  --latest-count 2 \\\n");
    fprintf(g, "  --write-final \\\n");
    fprintf(g, "  --roman-ref \\\n");
    fprintf(g, "  --theory-cs simple\n");
    fprintf(g, "```\n\n");
    fprintf(g, "Combined output:\n\n");
    fprintf(g, "```text\n");
    fprintf(g, "hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/combined_latest_2_analysis/final_plots/FINAL speed_of_sound_on_packing_fracture.pdf\n");
    fprintf(g, "```\n\n");
    fprintf(g, "## One-Command High/Mid/Low Workflow\n\n");
    fprintf(g, "This runs three eta-specific simulation batches and then creates one combined plot:\n\n");
    fprintf(g, "```sh\n");
    fprintf(g, "python3 hspist3/run_speed_of_sound_eta_split_workflow.py\n");
    fprintf(g, "```\n\n");
    fprintf(g, "Output is written to a timestamped archive folder:\n\n");
    fprintf(g, "```text\n");
    fprintf(g, "hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/simulation_eta_split_<timestamp>/final_plots/FINAL speed_of_sound_on_packing_fracture.pdf\n");
    fprintf(g, "```\n\n");
    fprintf(g, "## Main Output\n\n");
    fprintf(g, "```text\n");
    fprintf(g, "hspist3/%s/analysis_by_eta/final_plots/FINAL speed_of_sound_on_packing_fracture.pdf\n", dir_path);
    fprintf(g, "```\n");
    fclose(g);
}

static void write_latest_txt(const char *base_dir, const char *run_dir) {
    if (!base_dir || !*base_dir || !run_dir || !*run_dir) return;
    char out_path[1024];
    snprintf(out_path, sizeof(out_path), "%s/LATEST.txt", base_dir);
    FILE *f = fopen(out_path, "w");
    if (!f) return;
    fprintf(f, "%s\n", run_dir);
    fclose(f);
}

static void szilard_setup_timestamped_run_dir(void) {
    g_szilard_run_active = 0;
    g_szilard_run_dir[0] = '\0';
    if (cli_experiment_preset != EXPERIMENT_PRESET_SZILARD_ENGINE) return;

    time_t now = time(NULL);
    struct tm *t = localtime(&now);
    char stamp[64];
    if (t) {
        /* Match speed-of-sound naming style for consistency. */
        snprintf(stamp, sizeof(stamp), "simulation_%02d_%02d_%02d_%02d_%02d",
                 t->tm_mday, t->tm_mon + 1, t->tm_year % 100, t->tm_hour, t->tm_min);
    } else {
        snprintf(stamp, sizeof(stamp), "simulation_%u", (unsigned int)now);
    }

    snprintf(g_szilard_run_dir, sizeof(g_szilard_run_dir), "%s/%s", g_szilard_dir, stamp);
    mkdir_p(g_szilard_run_dir);
    write_command_md_to_dir(g_szilard_run_dir, g_main_argc, g_main_argv);
    write_latest_txt(g_szilard_dir, g_szilard_run_dir);
    g_szilard_run_active = 1;

    printf("📁 Szilard run dir: %s\n", g_szilard_run_dir);

    /* Save plot commands (manual; no auto plotting). */
    {
        char plot_md[1100];
        snprintf(plot_md, sizeof(plot_md), "%s/01_PLOT_COMMANDS.md", g_szilard_run_dir);
        FILE *f = fopen(plot_md, "w");
        if (f) {
            fprintf(f, "# Plot commands (manual)\n\n");
            fprintf(f, "Run these from the `hspist3/` directory (repo-relative paths):\n\n");
            fprintf(f, "```sh\n");
            fprintf(f, "cd hspist3\n");
            fprintf(f, "# Szilard info/work/heat log (auto-written by default):\n");
            fprintf(f, "python3 plot_szilard_info_thermo.py --csv %s/szilard_info.csv --out %s/plots\n",
                    g_szilard_run_dir, g_szilard_run_dir);
            fprintf(f, "# (Optional) Inspect wall motion trace:\n");
            fprintf(f, "python3 wall_x_FFT.py --csv %s/wall_position.csv\n", g_szilard_run_dir);
            fprintf(f, "```\n\n");
            fprintf(f, "Notes:\n");
            fprintf(f, "- `szilard_info.csv` includes branch labels so G/H/J/restores from one session stay distinguishable.\n");
            fprintf(f, "- Logs always written here: `energy_log.csv`, `wall_position.csv`, `run_params.json`.\n");
            fclose(f);
        }
    }
}

static void init_output_dirs_from_argv0(const char *argv0) {
    if (!getcwd(g_launch_cwd, sizeof(g_launch_cwd))) {
        snprintf(g_launch_cwd, sizeof(g_launch_cwd), ".");
    }
    snprintf(g_program_dir, sizeof(g_program_dir), ".");
    snprintf(g_simple_box_dir, sizeof(g_simple_box_dir), "experiments_simple_box");
    snprintf(g_energy_transfer_dir, sizeof(g_energy_transfer_dir), "experiments_energy_transfer");
    snprintf(g_szilard_dir, sizeof(g_szilard_dir), "experiments_szilard");

    if (!argv0 || !*argv0) return;
    const char *slash = strrchr(argv0, '/');
    if (!slash) return;
    size_t dlen = (size_t)(slash - argv0);
    if (dlen == 0) return;
    if (dlen >= sizeof(g_program_dir)) dlen = sizeof(g_program_dir) - 1;
    memcpy(g_program_dir, argv0, dlen);
    g_program_dir[dlen] = '\0';

    if (strcmp(g_program_dir, ".") == 0) {
        snprintf(g_simple_box_dir, sizeof(g_simple_box_dir), "experiments_simple_box");
        snprintf(g_energy_transfer_dir, sizeof(g_energy_transfer_dir), "experiments_energy_transfer");
        snprintf(g_szilard_dir, sizeof(g_szilard_dir), "experiments_szilard");
    } else {
        snprintf(g_simple_box_dir, sizeof(g_simple_box_dir), "%s/experiments_simple_box", g_program_dir);
        snprintf(g_energy_transfer_dir, sizeof(g_energy_transfer_dir), "%s/experiments_energy_transfer", g_program_dir);
        snprintf(g_szilard_dir, sizeof(g_szilard_dir), "%s/experiments_szilard", g_program_dir);
    }
}

static void resolve_cli_output_path(char *dst, size_t dst_sz, const char *cli_path) {
    if (!dst || dst_sz == 0) return;
    dst[0] = '\0';
    if (!cli_path || !cli_path[0]) return;
    if (cli_path[0] == '/') {
        snprintf(dst, dst_sz, "%s", cli_path);
    } else {
        snprintf(dst, dst_sz, "%s/%s", g_launch_cwd, cli_path);
    }
}

static void ensure_csv_header_schema(const char *path, const char *expected_header) {
    if (!path || !expected_header) return;
    FILE *f = fopen(path, "r");
    if (!f) return;
    char line[4096];
    if (!fgets(line, sizeof(line), f)) {
        fclose(f);
        return;
    }
    fclose(f);

    line[strcspn(line, "\r\n")] = '\0';
    char expected[4096];
    snprintf(expected, sizeof(expected), "%s", expected_header);
    expected[strcspn(expected, "\r\n")] = '\0';

    if (strcmp(line, expected) == 0) return;

    char backup[1024];
    time_t now = time(NULL);
    snprintf(backup, sizeof(backup), "%s.legacy_%ld", path, (long)now);
    if (rename(path, backup) != 0) {
        if (!cli_quiet) {
            printf("⚠️ Could not rotate CSV with mismatched header: %s\n", path);
        }
    } else if (!cli_quiet) {
        printf("ℹ️ Rotated old CSV with mismatched header → %s\n", backup);
    }
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
    free(all_wall_masses); all_wall_masses = NULL;
    free(extra_wall_positions); extra_wall_positions = NULL;
    free(extra_wall_old_positions); extra_wall_old_positions = NULL;
    free(extra_wall_velocity); extra_wall_velocity = NULL;
    free(extra_wall_masses); extra_wall_masses = NULL;
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
    all_wall_masses = (float *)cli_checked_realloc(all_wall_masses,
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

    // Default: use leftmost wall as the primary (spring anchor + dynamic wall)
    leftmost_wall_index = 0;
    float min_pos = all_wall_positions[0];
    for (int w = 1; w < num_internal_walls; ++w) {
        float wx = all_wall_positions[w];
        if (wx < min_pos) {
            min_pos = wx;
            leftmost_wall_index = w;
        }
    }
    primary_wall_index = leftmost_wall_index;

    // Assign per-wall masses (absolute mass units) from CLI wall-mass-factors when present.
    // If fewer masses than walls were provided, repeat the last provided mass.
    // Otherwise, default to the current wall_mass_runtime (single-mass mode).
    double base_mf = (double)wall_mass_runtime / (double)PARTICLE_MASS;
    if (!isfinite(base_mf) || base_mf <= 0.0) base_mf = (double)WALL_MASS_FACTOR;
    for (int w = 0; w < num_internal_walls; ++w) {
        double mf = base_mf;
        if (cli_wall_masses_set && cli_wall_masses_count > 0) {
            int idx = (w < cli_wall_masses_count) ? w : (cli_wall_masses_count - 1);
            if (idx < 0) idx = 0;
            if (idx < (int)(sizeof(cli_wall_masses) / sizeof(cli_wall_masses[0])) && cli_wall_masses[idx] > 0.0f) {
                mf = (double)cli_wall_masses[idx];
            }
        }
        if (!isfinite(mf) || mf <= 0.0) mf = base_mf;
        all_wall_masses[w] = (float)((double)PARTICLE_MASS * mf);
    }

    wall_x = all_wall_positions[primary_wall_index];
    wall_x_old = wall_x;
    vx_wall = 0.0f;
    wall_impulse_x_accum = 0.0f;
    wall_is_released = false;
    wall_hold_enabled = true;
    primary_wall_release_time = -1.0;
    wall_mass_runtime = all_wall_masses[primary_wall_index];

    extra_wall_count = num_internal_walls - 1;
    if (extra_wall_count > 0) {
        extra_wall_positions = (float *)cli_checked_realloc(extra_wall_positions,
                                                           (size_t)extra_wall_count * sizeof(float));
        extra_wall_old_positions = (float *)cli_checked_realloc(extra_wall_old_positions,
                                                                (size_t)extra_wall_count * sizeof(float));
        extra_wall_velocity = (float *)cli_checked_realloc(extra_wall_velocity,
                                                           (size_t)extra_wall_count * sizeof(float));
        extra_wall_masses = (float *)cli_checked_realloc(extra_wall_masses,
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
            extra_wall_masses[idx] = all_wall_masses[w];
            extra_wall_hold_enabled[idx] = true;
            extra_wall_is_released[idx] = false;
            extra_wall_release_time[idx] = -1.0;
            ++idx;
        }
    } else {
        free(extra_wall_positions); extra_wall_positions = NULL;
        free(extra_wall_old_positions); extra_wall_old_positions = NULL;
        free(extra_wall_velocity); extra_wall_velocity = NULL;
        free(extra_wall_masses); extra_wall_masses = NULL;
        free(extra_wall_hold_enabled); extra_wall_hold_enabled = NULL;
        free(extra_wall_is_released); extra_wall_is_released = NULL;
        free(extra_wall_release_time); extra_wall_release_time = NULL;
    }

    segment_count = num_internal_walls + 1;
    segment_counts = (int *)cli_checked_realloc(segment_counts, (size_t)segment_count * sizeof(int));
    segment_ke = (double *)cli_checked_realloc(segment_ke, (size_t)segment_count * sizeof(double));
    segment_temperature = (double *)cli_checked_realloc(segment_temperature, (size_t)segment_count * sizeof(double));
    segment_vx_mean = (double *)cli_checked_realloc(segment_vx_mean, (size_t)segment_count * sizeof(double));
    segment_vx_var = (double *)cli_checked_realloc(segment_vx_var, (size_t)segment_count * sizeof(double));
    segment_vy_mean = (double *)cli_checked_realloc(segment_vy_mean, (size_t)segment_count * sizeof(double));
    segment_vy_var = (double *)cli_checked_realloc(segment_vy_var, (size_t)segment_count * sizeof(double));
    segment_speed_mean = (double *)cli_checked_realloc(segment_speed_mean, (size_t)segment_count * sizeof(double));
    segment_speed_var = (double *)cli_checked_realloc(segment_speed_var, (size_t)segment_count * sizeof(double));
    memset(segment_counts, 0, (size_t)segment_count * sizeof(int));
    memset(segment_ke, 0, (size_t)segment_count * sizeof(double));
    memset(segment_temperature, 0, (size_t)segment_count * sizeof(double));
    memset(segment_vx_mean, 0, (size_t)segment_count * sizeof(double));
    memset(segment_vx_var, 0, (size_t)segment_count * sizeof(double));
    memset(segment_vy_mean, 0, (size_t)segment_count * sizeof(double));
    memset(segment_vy_var, 0, (size_t)segment_count * sizeof(double));
    memset(segment_speed_mean, 0, (size_t)segment_count * sizeof(double));
    memset(segment_speed_var, 0, (size_t)segment_count * sizeof(double));

    leftmost_wall_baseline = (num_internal_walls > 0)
                             ? all_wall_positions[leftmost_wall_index]
                             : wall_x;
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
    // Use piston planes as dynamic outer boundaries when pistons are active.
    // Left piston plane is at piston_left_x + 5 (matches CCD/EDMD convention).
    // Right piston plane is at piston_right_x.
    float left_plane = fmaxf((float)XW1, piston_left_x + 5.0f);
    float right_plane = fminf((float)XW2, piston_right_x);
    if (right_plane < left_plane + 1.0f) right_plane = left_plane + 1.0f;

    float prev = left_plane;
    for (int s = 0; s < num_internal_walls; ++s) {
        float boundary = tmp[s];
        left_bounds[s] = prev;
        right_bounds[s] = boundary - WALL_THICKNESS * 0.5f;
        prev = boundary + WALL_THICKNESS * 0.5f;
    }
    left_bounds[segs - 1] = prev;
    right_bounds[segs - 1] = right_plane;
}

// GUI-only variant: compute segment bounds using boundary positions extrapolated by gui_render_dt_rem.
// This makes moving walls/pistons look smoother in the live histogram panels when fixed_dt is not tiny.
// IMPORTANT: visualization only; physics uses compute_segment_bounds().
static void get_sorted_wall_positions_display(float *out_positions) {
    if (num_internal_walls <= 0) return;
    if (!out_positions) return;

    // Build instantaneous+extrapolated wall positions (px, global coords).
    int idx = 0; // index into extra_* arrays (excludes primary wall)
    for (int w = 0; w < num_internal_walls; ++w) {
        float xw = wall_x;
        float vxw = vx_wall;
        if (w != primary_wall_index) {
            if (extra_wall_count > 0 && extra_wall_positions && extra_wall_velocity && idx < extra_wall_count) {
                xw = extra_wall_positions[idx];
                vxw = extra_wall_velocity[idx];
                idx++;
            }
        }
        out_positions[w] = xw + vxw * gui_render_dt_rem;
    }

    // insertion sort (wall count small)
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

static void compute_segment_bounds_display(float *left_bounds, float *right_bounds) {
    int segs = num_internal_walls + 1;
    if (segs <= 0) return;
    float tmp[MAX_DIVIDER_CAPACITY];
    if (num_internal_walls > (int)(sizeof(tmp) / sizeof(tmp[0]))) {
        for (int i = 0; i < segs; ++i) {
            left_bounds[i] = (float)XW1;
            right_bounds[i] = (float)XW2;
        }
        return;
    }

    get_sorted_wall_positions_display(tmp);

    // Extrapolate piston planes for the same dt remainder.
    const float pistonL_x_disp = piston_left_x + vx_piston_left * gui_render_dt_rem;
    const float pistonR_x_disp = piston_right_x + vx_piston_right * gui_render_dt_rem;
    float left_plane = fmaxf((float)XW1, pistonL_x_disp + 5.0f);
    float right_plane = fminf((float)XW2, pistonR_x_disp);
    if (right_plane < left_plane + 1.0f) right_plane = left_plane + 1.0f;

    float prev = left_plane;
    for (int s = 0; s < num_internal_walls; ++s) {
        float boundary = tmp[s];
        left_bounds[s] = prev;
        right_bounds[s] = boundary - WALL_THICKNESS * 0.5f;
        prev = boundary + WALL_THICKNESS * 0.5f;
    }
    left_bounds[segs - 1] = prev;
    right_bounds[segs - 1] = right_plane;
}

// Build semicolon-separated instantaneous packing fractions per segment.
// Uses current segment particle counts and geometric segment widths.
static void build_segment_eta_csv(char *buf, size_t bufsz) {
    if (!buf || bufsz == 0) return;
    buf[0] = '\0';
    if (!(segment_counts && segment_count > 0)) return;
    if (segment_count > MAX_DIVIDER_CAPACITY + 1) return;

    float seg_left[MAX_DIVIDER_CAPACITY + 1];
    float seg_right[MAX_DIVIDER_CAPACITY + 1];
    compute_segment_bounds(seg_left, seg_right);

    const double Hpx = (double)(YW2 - YW1);
    const double rpx = (double)PARTICLE_RADIUS;
    if (!(Hpx > 0.0) || !(rpx > 0.0)) return;

    size_t off = 0;
    for (int s = 0; s < segment_count && off + 32 < bufsz; ++s) {
        const double wpx = (double)seg_right[s] - (double)seg_left[s];
        const double eta_seg = (wpx > 1e-12)
            ? ((double)segment_counts[s] * M_PI * rpx * rpx) / (wpx * Hpx)
            : NAN;
        int n = snprintf(buf + off, bufsz - off, "%s%.6g", (s == 0) ? "" : ";", eta_seg);
        if (n < 0) break;
        off += (size_t)n;
    }
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
    printf("                               simple_prediction_box - Single gas box with stochastic piston protocol\n");
    printf("                               simple_gas_box    - Alias for simple_prediction_box\n");
    printf("                               particlelife       - Particle-life artificial life simulation\n");
    printf("\nInteractive control\n");
    printf("  --show-simulation           Open SDL window (alias of --no-experiments)\n");
    printf("  --no-gui, --headless        Run without SDL window (faster, still writes CSV)\n");
    printf("  --left-empty                Seed zero particles in the leftmost segment\n");
    printf("  --no-pp-collisions          Disable particle–particle collisions (ideal gas limit)\n");
    printf("  --seed=value                RNG seed (integer) or 'time'\n");
    printf("                               multi_wall         - Multi-wall system analysis\n");
    printf("                               szilard_engine     - Szilard engine entropy experiments\n");
    printf("                               atp_synthase       - ATP synthase analog experiments\n");
    printf("                               protocol_optimization - Find optimal piston protocols\n");
    printf("                               dynamic_ib         - Dynamic information bottleneck\n");
    printf("                               wall_mid           - Single wall in middle\n");
    printf("                               2wall_exp_with_1_wall, 2wall_exp_with_2_wall\n");
    printf("  --lengths=L0,...            Override experiment half-lengths (sigma units)\n");
    printf("  --wall-masses=f1,...        Override experiment wall-mass factors (multiples of particle mass)\n");
    printf("  --speed-sound-run-dir=PATH  Speed-of-sound: write mode1 CSVs into this exact directory\n");
    printf("  --num-walls=N               Number of walls (1-5, default 1)\n");
    printf("  --wall-positions=x1,x2,...  Wall positions from left to right (sigma units)\n");
    printf("  --wall-mass-factors=m1,m2,... Wall mass factors from left to right (default 200 each)\n");
    printf("  --seeding=MODE              Seeding mode: grid (default), honeycomb, random\n");
    printf("  --distribute=area           Distribute particles ∝ segment area (width×height)\n");
    printf("  --protocol=NAME             Piston protocol: step, sigmoidal, linear, sinusoidal, optimal\n");
    printf("  --piston-speed=value        Piston max speed for protocol (σ per σ-time)\n");
    printf("  --piston-duration=value     Ramp duration for non-step protocols (σ-time units)\n");
    printf("  --piston-travel-sigma=value Right piston travel distance for step trigger (σ)\n");
    printf("  --energy-transfer-summary=PATH  Append energy-transfer run summary rows to PATH (CSV)\n");
    printf("  --energy-transfer-trace=PATH    Write energy-transfer trace CSV to PATH (default: experiments_energy_transfer/energy_transfer_trace.csv)\n");
    printf("  --piston-right-protocol-mode=MODE   step|linear|sigmoidal (alias for --protocol)\n");
    printf("  --max-right-piston-travel=value     Right piston travel (σ) (alias for --piston-travel-sigma)\n");
    printf("  --velocity-right-piston-step=value  Step protocol constant speed (σ/σ-time)\n");
    printf("  --velocity-right-piston-t0=value    Initial speed for linear/sigmoidal (σ/σ-time)\n");
    printf("  --velocity-right-piston-gradient=value Linear accel (σ/(σ-time)^2)\n");
    printf("  --velocity-right-piston-max=value   Max speed cap (σ/σ-time)\n");
    printf("  --piston-right-duration=value       Duration for linear/sigmoidal (σ-time)\n");
    printf("  --piston-right-sigmoid-steepness=value Sigmoid steepness (dimensionless)\n");
    printf("  --piston-work-target=value      Stop piston when ΔW_in reaches target (same units as PistonWork Δ)\n");
    printf("  --spring-k=value            Spring constant for energy measurement (code units)\n");
    printf("  --spring-eq=value          Spring equilibrium position (sigma units from left)\n");
    printf("  --calibrate-work=target,tol,min,max  Approximate place rightmost wall to match target piston work\n");
    printf("  --repeats=N                 Number of repeats per experiment configuration\n");
    printf("  --steps=N                   Steps to simulate after wall release (default derived from TIME_UNITS_SIMULATION)\n");
    printf("  --single-test               Run headless (no SDL window) but output to main directory for quick testing\n");
    printf("  --l0=value                  Set initial half-length for interactive mode\n");
    printf("  --height=value              Set box height (in sigma units) for interactive mode\n");
    printf("  --wall-mass-factor=value    Set wall mass factor for interactive mode (M_wall = factor * m_particle)\n");
        printf("  --wall-thickness=value      Set wall thickness in sigma units (physics)\n");
        printf("  --wall-thickness-vis=value  Set visual wall thickness in sigma units (render-only)\n");
        printf("  --wall-hold-steps=value    Steps to hold wall before release (default %d)\n", WALL_HOLD_STEPS);
        printf("  --auto-release             Start with wall released (no hold)\n");
        printf("  --auto-piston-step         Auto-trigger piston step (like pressing 't')\n");
        printf("  --auto-release-after-hold  Auto-release wall when hold steps reached (GUI)\n");
        printf("  --quiet                    Suppress per-particle seeding logs (keeps progress updates)\n");
    printf("  --particle-radius=value     Set particle radius in σ units (default 0.5)\n");
    printf("  --particle-diameter=value   Set particle diameter in σ units\n");
    printf("  --eta=value                 Auto-calculate radius for desired packing fraction (0 < eta < 0.9069)\n");
    printf("  --packing-fraction=value    Alias for --eta\n");
    printf("  --particles=N               Total active particles (<= MAX buffer)\n");
    printf("  --phi-max=value             Max packing fraction cap for seeding (default 0.78)\n");
    printf("  --particles-boxes=a,b,...   Target counts per segment (left→right)\n");
    printf("  --particles-box-left=n      Shorthand: left count for 2 boxes\n");
    printf("  --particles-box-right=n     Shorthand: right count for 2 boxes\n");
    printf("  --pl-tensor=FILE            Load custom particle-life interaction tensor from file\n");
    printf("  --energy-measurement        Enable spring-based energy measurement\n");
    printf("  --colourcode=mode           Particle coloring: none (default), temperature, velocity\n");
    printf("  --colorcode=mode            Alias of --colourcode\n");
    printf("  --memory-view=mode          Szilard particle view: species (default), memory, overlay\n");
    printf("  --memory=mode               Alias for --memory-view\n");
    printf("  --temperature=value         Set initial gas temperature (reduced units)\n");
    printf("  --hb-temp=value             Set heat bath temperature and auto-enable it (reduced units)\n");
    printf("  --timescale=value           Speed multiplier for interactive run (default 1.0)\n");
    printf("  --output-dt=value           CSV logging interval (sim units), 0=every step\n");
    printf("  --fixed-dt=value            Override fixed timestep (σ-time units)\n");
    printf("  --substeps=N                Collision substeps per fixed step (reduces tunneling at high speeds)\n");
    printf("  --substeps-auto             Adaptive collision substeps (auto-increase during shocks/collapse)\n");
    printf("  --substeps-auto-min=N       Min substeps when --substeps-auto is enabled (default %d)\n", cli_substeps_auto_min);
    printf("  --substeps-auto-max=N       Max substeps when --substeps-auto is enabled (default %d)\n", cli_substeps_auto_max);
    printf("  --substeps-auto-dx-frac=f   Max relative displacement per substep as fraction of diameter (default %.2f)\n", (double)cli_substeps_auto_dx_frac);
    printf("  --pp-toi                    Use particle–particle time-of-impact correction (slower, more accurate)\n");
    printf("  --pp-toi-max-events=N       Max PP-TOI correction iterations per substep (default %d)\n", cli_pp_toi_max_events);
    printf("  --pp-toi-neighbor-ring=N    Grid neighbor radius for PP-TOI (default %d)\n", cli_pp_toi_neighbor_ring);
    printf("  --dist-log                  Log x/y/speed distributions (energy_transfer only)\n");
    printf("  --dist-dt=value             Distribution log cadence (sigma-time), 0=every step\n");
    printf("  --dist-bins=N               Histogram bins for distributions (default %d)\n", cli_dist_bins);
    printf("  --dist-max-speed=value      Speed histogram vmax (sigma-time), 0=auto\n");
    printf("  --dist-log-path=PATH        Distribution CSV path (supports {seed} and {ts})\n");
    printf("  --eff-window=T              Energy-transfer: window length Tw after piston stops for Gain_win (sigma-time, default %.3g)\n", (double)cli_eff_window_sigma);
    printf("  --eff-stop-after-window     Energy-transfer: stop run after Tw completes (faster EDMD)\n");
    printf("  --eff-output=MODE           Energy-transfer: output metric for W_out (spring | wall-ke; default spring)\n");
    printf("  --szilard-equil-steps=N     Szilard: equilibration steps before measurement (default %d)\n", cli_szilard_equil_steps);
    printf("  --szilard-sep-duration=T    Szilard: separation/recombine duration (sigma-time, default %.3f)\n", cli_szilard_sep_duration);
    printf("  --szilard-species-ratio=r   Szilard: fraction of species-0 (default %.2f)\n", (double)cli_szilard_species_ratio);
    printf("  --szilard-perm-swap=0|1     Szilard: swap permeability on recombine (default %d)\n", cli_szilard_perm_swap);
    printf("  --szilard-unstick=0|1       Szilard: nudge pinned particles off outer walls (default %d)\n", cli_szilard_unstick);
    printf("  --szilard-unstick-hits=N    Szilard: consecutive boundary hits before nudge (default %d)\n", cli_szilard_unstick_hits);
    printf("  --szilard-unstick-vx-frac=f Szilard: redirect fraction of speed into x (default %.2f)\n", (double)cli_szilard_unstick_vx_frac);
    printf("  --szilard-log=PATH          Szilard: write info/work/heat log CSV to PATH\n");
    printf("  --szilard-log-dt=T          Szilard: log cadence (sigma-time), 0=every step (default %.3f)\n", cli_szilard_log_dt);
    printf("  --szilard-n=N               Szilard: override particle count (alias for --particles)\n");
    printf("  --simple-box-cycles=N       Simple box: number of drive/relax protocol steps (default %d)\n", cli_simple_box_cycles);
    printf("  --simple-box-equil-steps=N  Simple box: equilibration steps before logging state 0 (default %d)\n", cli_simple_box_equil_steps);
    printf("  --simple-box-relax-steps=N  Simple box: relaxation steps after each drive (default %d)\n", cli_simple_box_relax_steps);
    printf("  --simple-box-log-dt=T       Simple box: trace cadence in sigma-time, 0=every step (default %.3f)\n", cli_simple_box_log_dt);
    printf("  --simple-box-target-protocol=MODE  Simple box: iid | correlated | periodic (default %s)\n", cli_simple_box_target_protocol);
    printf("  --simple-box-target-min=value  Simple box: minimum right-piston position x_t in sigma from left wall\n");
    printf("  --simple-box-target-max=value  Simple box: maximum right-piston position x_t in sigma from left wall\n");
    printf("  --simple-box-initial-target=value  Simple box: initial equilibration piston position in sigma from left wall\n");
    printf("  --simple-box-corr-alpha=value   Simple box: persistence of correlated protocol in [0,1] (default %.2f)\n", cli_simple_box_corr_alpha);
    printf("  --simple-box-noise=value        Simple box: stochastic target noise scale in sigma (default %.3f)\n", cli_simple_box_noise_sigma);
    printf("  --simple-box-period=N           Simple box: period for periodic protocol (default %d)\n", cli_simple_box_period);
    printf("  --simple-box-run-dir=PATH       Simple box: write this run into a specific directory\n");
    printf("  --simple-box-run-label=NAME     Simple box: prefix timestamped run directory with NAME\n");
    printf("  --simple-box-render-every=N     Simple box GUI: render every N simulation steps (default %d)\n", cli_simple_box_render_every);
    printf("  --simple-box-gui-delay-ms=N     Simple box GUI: delay after each rendered frame (default %d ms)\n", cli_simple_box_gui_delay_ms);
    printf("  --time=value                Stop after this many sim units post-release (interactive)\n");
    printf("  --kbt1                      Force k_B*T = 1 (reduced units); adjusts k_B at runtime\n");
    printf("  --edmd-acc=0|1              Use accelerated EDMD backend (requires it built-in; default 0)\n");
    printf("  --mode=NAME                 time (default) | rk4 | edmd | edmd-hybrid | hybrid | event-split\n");
    printf("  --help                      Show this help message and exit\n");
}

static void parse_cli_options(int argc, char **argv) {
    preset_two_wall_active_walls = 2;  // default to both internal walls active for preset

    // Capture full command line for downstream logs.
    if (cli_command_line) {
        free(cli_command_line);
        cli_command_line = NULL;
    }
    {
        size_t total = 0;
        for (int a = 0; a < argc; ++a) total += strlen(argv[a]) + 1;
        cli_command_line = (char *)malloc(total + 1);
        if (!cli_command_line) { fprintf(stderr, "OOM building command line\n"); exit(EXIT_FAILURE); }
        cli_command_line[0] = '\0';
        for (int a = 0; a < argc; ++a) {
            strcat(cli_command_line, argv[a]);
            if (a + 1 < argc) strcat(cli_command_line, " ");
        }
    }

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
        } else if (strcmp(arg, "--no-gui") == 0 || strcmp(arg, "--headless") == 0) {
            // ##CHRIS: Run without SDL window (faster, still writes CSV files)
            cli_headless = true;
        } else if (strcmp(arg, "--left-empty") == 0) {
            cli_left_empty = true;
        } else if (strcmp(arg, "--no-pp-collisions") == 0) {
            cli_no_pp_collisions = true;
        } else if (strncmp(arg, "--speed-sound-run-dir", strlen("--speed-sound-run-dir")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            if (cli_speed_sound_run_dir) { free(cli_speed_sound_run_dir); cli_speed_sound_run_dir = NULL; }
            if (value && value[0]) cli_speed_sound_run_dir = strdup(value);
        } else if (strcmp(arg, "--szilard-no-pp") == 0) {
            cli_szilard_no_pp_during_g = true; /* ##CHRIS: disable PP during G phase only */
        } else if (strcmp(arg, "--dist-log") == 0) {
            cli_dist_log = 1;
        } else if (strncmp(arg, "--dist-dt", 9) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            cli_dist_dt = (float)atof(value);
        } else if (strncmp(arg, "--dist-bins", 11) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            cli_dist_bins = (int)atol(value);
            if (cli_dist_bins < 1) cli_dist_bins = 50;
            if (cli_dist_bins > 2000) cli_dist_bins = 2000;
        } else if (strncmp(arg, "--dist-max-speed", 16) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            cli_dist_vmax = (float)atof(value);
        } else if (strncmp(arg, "--dist-log-path", 15) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            if (cli_dist_log_path) { free(cli_dist_log_path); cli_dist_log_path = NULL; }
            if (value && value[0]) cli_dist_log_path = strdup(value);
        } else if (strncmp(arg, "--eff-window", strlen("--eff-window")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            double v = atof(value);
            if (v < 0.0) v = 0.0;
            if (v > 1e9) v = 1e9;
            cli_eff_window_sigma = (float)v;
            cli_eff_window_set = 1;
        } else if (strcmp(arg, "--eff-stop-after-window") == 0) {
            cli_eff_stop_after_window = 1;
        } else if (strncmp(arg, "--eff-output", strlen("--eff-output")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            if (!value || !value[0]) {
                fprintf(stderr, "Invalid --eff-output (empty). Use: spring | wall-ke\n");
                exit(EXIT_FAILURE);
            }
            if (strcmp(value, "spring") == 0) {
                cli_eff_output_mode = EFF_OUTPUT_SPRING;
            } else if (strcmp(value, "wall-ke") == 0 || strcmp(value, "wallke") == 0 || strcmp(value, "wall_ke") == 0) {
                cli_eff_output_mode = EFF_OUTPUT_WALL_KE;
            } else {
                fprintf(stderr, "Unknown --eff-output '%s'. Use: spring | wall-ke\n", value);
                exit(EXIT_FAILURE);
            }
            cli_eff_output_set = 1;
        } else if (strncmp(arg, "--szilard-equil-steps", strlen("--szilard-equil-steps")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            cli_szilard_equil_steps = (int)atol(value);
            if (cli_szilard_equil_steps < 0) cli_szilard_equil_steps = 0;
        } else if (strncmp(arg, "--szilard-sep-duration", strlen("--szilard-sep-duration")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            cli_szilard_sep_duration = atof(value);
            if (cli_szilard_sep_duration <= 0.0) cli_szilard_sep_duration = 1.0;
        } else if (strncmp(arg, "--szilard-species-ratio", strlen("--szilard-species-ratio")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            double v = atof(value);
            if (v < 0.0) v = 0.0;
            if (v > 1.0) v = 1.0;
            cli_szilard_species_ratio = (float)v;
        } else if (strncmp(arg, "--szilard-perm-swap", strlen("--szilard-perm-swap")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            cli_szilard_perm_swap = (int)atoi(value) ? 1 : 0;
        } else if (strncmp(arg, "--szilard-unstick", strlen("--szilard-unstick")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            if (value && value[0]) cli_szilard_unstick = (int)atoi(value) ? 1 : 0;
        } else if (strncmp(arg, "--szilard-unstick-hits", strlen("--szilard-unstick-hits")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            int v = (int)atol(value);
            if (v < 1) v = 1;
            if (v > 1000000) v = 1000000;
            cli_szilard_unstick_hits = v;
        } else if (strncmp(arg, "--szilard-unstick-vx-frac", strlen("--szilard-unstick-vx-frac")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            double v = atof(value);
            if (v < 0.0) v = 0.0;
            if (v > 0.95) v = 0.95;
            cli_szilard_unstick_vx_frac = (float)v;
        } else if (strncmp(arg, "--szilard-log-dt", strlen("--szilard-log-dt")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            double v = atof(value);
            if (v < 0.0) v = 0.0;
            if (v > 1e9) v = 1e9;
            cli_szilard_log_dt = v;
        } else if (strncmp(arg, "--szilard-log", strlen("--szilard-log")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            if (cli_szilard_log_path) { free(cli_szilard_log_path); cli_szilard_log_path = NULL; }
            if (value && value[0]) cli_szilard_log_path = strdup(value);
        } else if (strncmp(arg, "--szilard-n", strlen("--szilard-n")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            int v = (int)atol(value);
            if (v > 0) cli_szilard_n = v;
        } else if (strncmp(arg, "--simple-box-cycles", strlen("--simple-box-cycles")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            int v = (int)atol(value);
            if (v < 1) v = 1;
            cli_simple_box_cycles = v;
        } else if (strncmp(arg, "--simple-box-equil-steps", strlen("--simple-box-equil-steps")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            int v = (int)atol(value);
            if (v < 0) v = 0;
            cli_simple_box_equil_steps = v;
        } else if (strncmp(arg, "--simple-box-relax-steps", strlen("--simple-box-relax-steps")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            int v = (int)atol(value);
            if (v < 0) v = 0;
            cli_simple_box_relax_steps = v;
        } else if (strncmp(arg, "--simple-box-log-dt", strlen("--simple-box-log-dt")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            double v = atof(value);
            if (v < 0.0) v = 0.0;
            cli_simple_box_log_dt = v;
        } else if (strncmp(arg, "--simple-box-target-protocol", strlen("--simple-box-target-protocol")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            if (strcmp(value, "iid") != 0 &&
                strcmp(value, "correlated") != 0 &&
                strcmp(value, "periodic") != 0) {
                fprintf(stderr, "Unknown --simple-box-target-protocol '%s'. Use: iid | correlated | periodic\n", value);
                exit(EXIT_FAILURE);
            }
            snprintf(cli_simple_box_target_protocol, sizeof(cli_simple_box_target_protocol), "%s", value);
        } else if (strncmp(arg, "--simple-box-target-min", strlen("--simple-box-target-min")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            cli_simple_box_target_min_sigma = atof(value);
        } else if (strncmp(arg, "--simple-box-target-max", strlen("--simple-box-target-max")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            cli_simple_box_target_max_sigma = atof(value);
        } else if (strncmp(arg, "--simple-box-initial-target", strlen("--simple-box-initial-target")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            cli_simple_box_initial_target_sigma = atof(value);
        } else if (strncmp(arg, "--simple-box-corr-alpha", strlen("--simple-box-corr-alpha")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            double v = atof(value);
            if (v < 0.0) v = 0.0;
            if (v > 1.0) v = 1.0;
            cli_simple_box_corr_alpha = v;
        } else if (strncmp(arg, "--simple-box-noise", strlen("--simple-box-noise")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            double v = atof(value);
            if (v < 0.0) v = 0.0;
            cli_simple_box_noise_sigma = v;
        } else if (strncmp(arg, "--simple-box-period", strlen("--simple-box-period")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            int v = (int)atol(value);
            if (v < 2) v = 2;
            cli_simple_box_period = v;
        } else if (strncmp(arg, "--simple-box-run-dir", strlen("--simple-box-run-dir")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            if (cli_simple_box_run_dir) { free(cli_simple_box_run_dir); cli_simple_box_run_dir = NULL; }
            if (value && value[0]) cli_simple_box_run_dir = strdup(value);
        } else if (strncmp(arg, "--simple-box-run-label", strlen("--simple-box-run-label")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            if (cli_simple_box_run_label) { free(cli_simple_box_run_label); cli_simple_box_run_label = NULL; }
            if (value && value[0]) cli_simple_box_run_label = strdup(value);
        } else if (strncmp(arg, "--simple-box-render-every", strlen("--simple-box-render-every")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            int v = (int)atol(value);
            if (v < 1) v = 1;
            cli_simple_box_render_every = v;
        } else if (strncmp(arg, "--simple-box-gui-delay-ms", strlen("--simple-box-gui-delay-ms")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            int v = (int)atol(value);
            if (v < 0) v = 0;
            cli_simple_box_gui_delay_ms = v;
        } else if (strncmp(arg, "--seed", 6) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            if (strcmp(value, "time") == 0) {
                cli_seed_set = true;
                cli_seed = (unsigned int)time(NULL);
            } else {
                errno = 0;
                char *endptr = NULL;
                unsigned long v = strtoul(value, &endptr, 10);
                if (errno != 0 || endptr == value || *endptr != '\0') {
                    fprintf(stderr, "Invalid seed '%s'. Use an integer or 'time'.\n", value);
                    exit(EXIT_FAILURE);
                }
                cli_seed_set = true;
                cli_seed = (unsigned int)v;
            }
        } else if (strncmp(arg, "--distribute", 12) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            if (strcmp(value, "area") == 0) cli_distribute_area = 1; else if (strcmp(value, "equal") == 0) cli_distribute_area = 0; else {
                fprintf(stderr, "Unknown distribute mode '%s'. Use: area or equal\n", value); exit(EXIT_FAILURE);
            }
        } else if (strcmp(arg, "--energy-measurement") == 0) {
            cli_enable_energy_measurement = true;
        } else if (strcmp(arg, "--kbt1") == 0) {
            cli_force_kbt_one = 1;
        } else if (strcmp(arg, "--auto-release") == 0) {
            cli_auto_release_wall = 1;
        } else if (strcmp(arg, "--auto-release-after-hold") == 0) {
            cli_auto_release_after_hold = 1;
        } else if (strcmp(arg, "--single-test") == 0) {
            cli_single_test_mode = 1;
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
        } else if (strcmp(preset_name, "simple_prediction_box") == 0 ||
                   strcmp(preset_name, "simple_gas_box") == 0) {
            cli_experiment_preset = EXPERIMENT_PRESET_SIMPLE_BOX_PREDICTION;
            enable_speed_of_sound_experiments = true;
        } else if (strcmp(preset_name, "multi_wall") == 0) {
                cli_experiment_preset = EXPERIMENT_PRESET_MULTI_WALL_SYSTEM;
                enable_speed_of_sound_experiments = true;
            } else if (strcmp(preset_name, "szilard_engine") == 0) {
                cli_experiment_preset = EXPERIMENT_PRESET_SZILARD_ENGINE;
                enable_speed_of_sound_experiments = true;
            } else if (strcmp(preset_name, "szilard_engine_particle_separation") == 0) {
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
            } else if (strcmp(preset_name, "particlelife") == 0) {
                // ##CHRIS: Particle-life artificial life simulation
                cli_experiment_preset = EXPERIMENT_PRESET_PARTICLELIFE;
                printf("🦠 Particle-life experiment selected\n");
            } else {
                fprintf(stderr, "Unknown experiment preset '%s'.\n", value_raw);
                print_cli_usage(argv[0]);
                exit(EXIT_FAILURE);
            }

        } else if (strncmp(arg, "--pl-tensor", 11) == 0) {
            // ##CHRIS: Custom particle-life tensor file
            const char *value = cli_option_value(arg, argc, argv, &i);
            if (cli_pl_tensor_file) free(cli_pl_tensor_file);
            cli_pl_tensor_file = strdup(value);
            printf("📁 Particle-life tensor file: %s\n", cli_pl_tensor_file);

        } else if (strncmp(arg, "--cyclic-boundaries", 19) == 0) {
            // ##CHRIS: Override cyclic boundaries setting from tensor file
            const char *value = cli_option_value(arg, argc, argv, &i);
            cli_pl_cyclic_boundaries = atoi(value);
            printf("🔄 Cyclic boundaries CLI override: %s\n", cli_pl_cyclic_boundaries ? "ENABLED (wrapping)" : "DISABLED (reflective walls)");

        } else if (strncmp(arg, "--pl-dt", 7) == 0) {
            // ##CHRIS: Particle-life timestep override
            const char *value = cli_option_value(arg, argc, argv, &i);
            cli_pl_dt = atof(value);
            printf("⏱️  Particle-life timestep override: dt = %.6f\n", cli_pl_dt);

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
        } else if (strncmp(arg, "--output-dt", 11) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0; char *endptr = NULL; float v = strtof(value, &endptr);
            if (errno != 0 || endptr == value || v < 0.0f) { fprintf(stderr, "Invalid output-dt '%s'.\n", value); exit(EXIT_FAILURE);} 
            cli_output_dt = v;
        } else if (strncmp(arg, "--fixed-dt", 10) == 0 || strncmp(arg, "--dt", 4) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0; char *endptr = NULL; float v = strtof(value, &endptr);
            if (errno != 0 || endptr == value || v <= 0.0f) { fprintf(stderr, "Invalid fixed-dt '%s'.\n", value); exit(EXIT_FAILURE);} 
            cli_fixed_dt_set = 1;
            cli_fixed_dt = v;
        } else if (strcmp(arg, "--pp-toi") == 0) {
            cli_pp_toi = 1;
        } else if (strncmp(arg, "--pp-toi-max-events", strlen("--pp-toi-max-events")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0; char *endptr = NULL; long v = strtol(value, &endptr, 10);
            if (errno != 0 || endptr == value || v < 1 || v > 1000) { fprintf(stderr, "Invalid pp-toi-max-events '%s'.\n", value); exit(EXIT_FAILURE);} 
            cli_pp_toi_max_events = (int)v;
            cli_pp_toi = 1;
        } else if (strncmp(arg, "--pp-toi-neighbor-ring", strlen("--pp-toi-neighbor-ring")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0; char *endptr = NULL; long v = strtol(value, &endptr, 10);
            if (errno != 0 || endptr == value || v < 0 || v > 5) { fprintf(stderr, "Invalid pp-toi-neighbor-ring '%s'.\n", value); exit(EXIT_FAILURE);} 
            cli_pp_toi_neighbor_ring = (int)v;
            cli_pp_toi = 1;
        } else if (strcmp(arg, "--substeps-auto") == 0) {
            cli_substeps_auto = true;
        } else if (strncmp(arg, "--substeps-auto-min", strlen("--substeps-auto-min")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0; char *endptr = NULL; long v = strtol(value, &endptr, 10);
            if (errno != 0 || endptr == value || *endptr != '\0' || v < 1 || v > 500) {
                fprintf(stderr, "Invalid substeps-auto-min '%s' (must be 1..500).\n", value);
                exit(EXIT_FAILURE);
            }
            cli_substeps_auto_min = (int)v;
            cli_substeps_auto = true;
        } else if (strncmp(arg, "--substeps-auto-max", strlen("--substeps-auto-max")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0; char *endptr = NULL; long v = strtol(value, &endptr, 10);
            if (errno != 0 || endptr == value || *endptr != '\0' || v < 1 || v > 500) {
                fprintf(stderr, "Invalid substeps-auto-max '%s' (must be 1..500).\n", value);
                exit(EXIT_FAILURE);
            }
            cli_substeps_auto_max = (int)v;
            cli_substeps_auto = true;
        } else if (strncmp(arg, "--substeps-auto-dx-frac", strlen("--substeps-auto-dx-frac")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0; char *endptr = NULL; float v = strtof(value, &endptr);
            if (errno != 0 || endptr == value || !isfinite(v) || v <= 0.0f || v > 1.0f) {
                fprintf(stderr, "Invalid substeps-auto-dx-frac '%s' (must be 0<frac<=1).\n", value);
                exit(EXIT_FAILURE);
            }
            cli_substeps_auto_dx_frac = v;
            cli_substeps_auto = true;
        } else if (strncmp(arg, "--substeps", 10) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0; char *endptr = NULL; long v = strtol(value, &endptr, 10);
            if (errno != 0 || endptr == value || *endptr != '\0' || v < 1 || v > 500) {
                fprintf(stderr, "Invalid substeps '%s' (must be 1..500).\n", value);
                exit(EXIT_FAILURE);
            }
            cli_substeps_override = (int)v;
        } else if (strncmp(arg, "--time", 6) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0; char *endptr = NULL; float v = strtof(value, &endptr);
            if (errno != 0 || endptr == value || v <= 0.0f) { fprintf(stderr, "Invalid time '%s'.\n", value); exit(EXIT_FAILURE);} 
            cli_time_span = v;
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
            cli_wall_masses_set = 1;
            cli_wall_masses_count = mass_count;
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
        } else if (strncmp(arg, "--quiet", strlen("--quiet")) == 0) {
            cli_quiet = 1;
        } else if (strncmp(arg, "--energy-transfer-summary", strlen("--energy-transfer-summary")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            if (!value || value[0] == '\0') {
                fprintf(stderr, "Invalid --energy-transfer-summary (empty path)\n");
                exit(EXIT_FAILURE);
            }
            cli_energy_transfer_summary_set = true;
            snprintf(cli_energy_transfer_summary_path, sizeof(cli_energy_transfer_summary_path), "%s", value);
        } else if (strncmp(arg, "--energy-transfer-trace", strlen("--energy-transfer-trace")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            if (!value || value[0] == '\0') {
                fprintf(stderr, "Invalid --energy-transfer-trace (empty path)\n");
                exit(EXIT_FAILURE);
            }
            cli_energy_transfer_trace_set = true;
            snprintf(cli_energy_transfer_trace_path, sizeof(cli_energy_transfer_trace_path), "%s", value);
        } else if (strncmp(arg, "--auto-piston-step", strlen("--auto-piston-step")) == 0) {
            cli_auto_piston_step = true;
        } else if (strncmp(arg, "--piston-right-protocol-mode", strlen("--piston-right-protocol-mode")) == 0 ||
                   strncmp(arg, "--piston_right_protocol_mode", strlen("--piston_right_protocol_mode")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            if (!value) { fprintf(stderr, "Missing value for --piston-right-protocol-mode\n"); exit(EXIT_FAILURE); }
            cli_piston_right_mode_set = true;
            snprintf(cli_piston_right_mode, sizeof(cli_piston_right_mode), "%s", value);
            if (strcmp(value, "step") == 0) cli_protocol = PROTOCOL_STEP;
            else if (strcmp(value, "sigmoidal") == 0) cli_protocol = PROTOCOL_SIGMOIDAL;
            else if (strcmp(value, "linear") == 0) cli_protocol = PROTOCOL_LINEAR;
            else {
                fprintf(stderr, "Unknown piston-right-protocol-mode '%s'. Use: step, linear, sigmoidal\n", value);
                exit(EXIT_FAILURE);
            }
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
        } else if (strncmp(arg, "--max-right-piston-travel", strlen("--max-right-piston-travel")) == 0 ||
                   strncmp(arg, "--max_right_piston_travel", strlen("--max_right_piston_travel")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0; char *endptr = NULL; float v = strtof(value, &endptr);
            if (errno != 0 || endptr == value || *endptr != '\0' || v <= 0.0f) {
                fprintf(stderr, "Invalid --max-right-piston-travel '%s' (must be > 0, in σ).\n", value);
                exit(EXIT_FAILURE);
            }
            cli_piston_right_travel_set = true;
            cli_piston_right_travel_sigma = v;
            // Also feed the generic parameter so existing code paths work.
            cli_piston_travel_set = 1;
            cli_piston_travel_sigma = v;
        } else if (strncmp(arg, "--velocity-right-piston-step", strlen("--velocity-right-piston-step")) == 0 ||
                   strncmp(arg, "--velocity_right_piston_step_constant", strlen("--velocity_right_piston_step_constant")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0; char *endptr = NULL; float v = strtof(value, &endptr);
            if (errno != 0 || endptr == value || *endptr != '\0' || v <= 0.0f) {
                fprintf(stderr, "Invalid --velocity-right-piston-step '%s' (must be > 0, in σ/σ-time).\n", value);
                exit(EXIT_FAILURE);
            }
            cli_piston_right_step_speed_set = true;
            cli_piston_right_step_speed = v;
            cli_piston_speed_set = 1;
            cli_piston_speed = v;
        } else if (strncmp(arg, "--velocity-right-piston-t0", strlen("--velocity-right-piston-t0")) == 0 ||
                   strncmp(arg, "--velocity_right_piston_t0", strlen("--velocity_right_piston_t0")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0; char *endptr = NULL; float v = strtof(value, &endptr);
            if (errno != 0 || endptr == value || *endptr != '\0' || v < 0.0f) {
                fprintf(stderr, "Invalid --velocity-right-piston-t0 '%s' (must be >= 0, in σ/σ-time).\n", value);
                exit(EXIT_FAILURE);
            }
            cli_piston_right_v0_set = true;
            cli_piston_right_v0 = v;
        } else if (strncmp(arg, "--velocity-right-piston-gradient", strlen("--velocity-right-piston-gradient")) == 0 ||
                   strncmp(arg, "--velocity_right_piston_gradient", strlen("--velocity_right_piston_gradient")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0; char *endptr = NULL; float v = strtof(value, &endptr);
            if (errno != 0 || endptr == value || *endptr != '\0') {
                fprintf(stderr, "Invalid --velocity-right-piston-gradient '%s' (in σ/(σ-time)^2).\n", value);
                exit(EXIT_FAILURE);
            }
            cli_piston_right_gradient_set = true;
            cli_piston_right_gradient = v;
        } else if (strncmp(arg, "--velocity-right-piston-max", strlen("--velocity-right-piston-max")) == 0 ||
                   strncmp(arg, "--velocity_right_piston_max", strlen("--velocity_right_piston_max")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0; char *endptr = NULL; float v = strtof(value, &endptr);
            if (errno != 0 || endptr == value || *endptr != '\0' || v <= 0.0f) {
                fprintf(stderr, "Invalid --velocity-right-piston-max '%s' (must be > 0, in σ/σ-time).\n", value);
                exit(EXIT_FAILURE);
            }
            cli_piston_right_vmax_set = true;
            cli_piston_right_vmax = v;
            // For existing code paths, treat this as generic max speed.
            cli_piston_speed_set = 1;
            cli_piston_speed = v;
        } else if (strncmp(arg, "--piston-right-duration", strlen("--piston-right-duration")) == 0 ||
                   strncmp(arg, "--piston_right_duration", strlen("--piston_right_duration")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0; char *endptr = NULL; float v = strtof(value, &endptr);
            if (errno != 0 || endptr == value || *endptr != '\0' || v <= 0.0f) {
                fprintf(stderr, "Invalid --piston-right-duration '%s' (must be > 0, in σ-time).\n", value);
                exit(EXIT_FAILURE);
            }
            cli_piston_right_duration_set = true;
            cli_piston_right_duration = v;
            cli_piston_duration_set = 1;
            cli_piston_duration = v;
        } else if (strncmp(arg, "--piston-right-sigmoid-steepness", strlen("--piston-right-sigmoid-steepness")) == 0 ||
                   strncmp(arg, "--piston_right_sigmoid_steepness", strlen("--piston_right_sigmoid_steepness")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0; char *endptr = NULL; float v = strtof(value, &endptr);
            if (errno != 0 || endptr == value || *endptr != '\0' || v <= 0.0f) {
                fprintf(stderr, "Invalid --piston-right-sigmoid-steepness '%s' (must be > 0).\n", value);
                exit(EXIT_FAILURE);
            }
            cli_piston_right_sigmoid_steepness_set = true;
            cli_piston_right_sigmoid_steepness = v;
        } else if (strncmp(arg, "--piston-speed", strlen("--piston-speed")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0; char *endptr = NULL; float v = strtof(value, &endptr);
            if (errno != 0 || endptr == value || *endptr != '\0' || v <= 0.0f) {
                fprintf(stderr, "Invalid --piston-speed '%s' (must be > 0, in σ per σ-time).\n", value);
                exit(EXIT_FAILURE);
            }
            cli_piston_speed_set = 1;
            cli_piston_speed = v;
        } else if (strncmp(arg, "--piston-duration", strlen("--piston-duration")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0; char *endptr = NULL; float v = strtof(value, &endptr);
            if (errno != 0 || endptr == value || *endptr != '\0' || v <= 0.0f) {
                fprintf(stderr, "Invalid --piston-duration '%s' (must be > 0, in σ-time units).\n", value);
                exit(EXIT_FAILURE);
            }
            cli_piston_duration_set = 1;
            cli_piston_duration = v;
        } else if (strncmp(arg, "--piston-travel-sigma", strlen("--piston-travel-sigma")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0; char *endptr = NULL; float v = strtof(value, &endptr);
            if (errno != 0 || endptr == value || *endptr != '\0' || v <= 0.0f) {
                fprintf(stderr, "Invalid --piston-travel-sigma '%s' (must be > 0, in σ).\n", value);
                exit(EXIT_FAILURE);
            }
            cli_piston_travel_set = 1;
            cli_piston_travel_sigma = v;
        } else if (strncmp(arg, "--piston-work-target", strlen("--piston-work-target")) == 0 ||
                   strncmp(arg, "--piston_work_target", strlen("--piston_work_target")) == 0 ||
                   strncmp(arg, "--piston-w-in-target", strlen("--piston-w-in-target")) == 0 ||
                   strncmp(arg, "--piston_w_in_target", strlen("--piston_w_in_target")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0; char *endptr = NULL; double v = strtod(value, &endptr);
            if (errno != 0 || endptr == value || *endptr != '\0' || !(v > 0.0)) {
                fprintf(stderr, "Invalid --piston-work-target '%s' (must be > 0, in work/energy units).\n", value);
                exit(EXIT_FAILURE);
            }
            cli_piston_work_target_set = 1;
            cli_piston_work_target = v;
        } else if (strncmp(arg, "--wall-hold-steps", strlen("--wall-hold-steps")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0; char *endptr = NULL; long v = strtol(value, &endptr, 10);
            if (errno != 0 || endptr == value || *endptr != '\0' || v < 0) {
                fprintf(stderr, "Invalid --wall-hold-steps '%s' (must be >= 0)\n", value);
                exit(EXIT_FAILURE);
            }
            wall_hold_steps = (int)v;
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
        } else if (strncmp(arg, "--eta", strlen("--eta")) == 0 ||
                   strncmp(arg, "--packing-fraction", strlen("--packing-fraction")) == 0) {
            // ##CHRIS: Auto-calculate particle radius from desired packing fraction
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0; char *endptr = NULL; float v = strtof(value, &endptr);
            if (errno != 0 || endptr == value || *endptr != '\0' || v <= 0.0f || v >= 0.9069f) {
                fprintf(stderr, "Invalid packing fraction '%s' (must be 0 < eta < 0.9069).\n", value);
                exit(EXIT_FAILURE);
            }
            cli_packing_fraction_set = true;
            cli_packing_fraction = v;
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
        } else if (strncmp(arg, "--wall-thickness-vis", 21) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0;
            char *endptr = NULL;
            float v = strtof(value, &endptr);
            if (errno != 0 || endptr == value || *endptr != '\0' || v <= 0.0f) {
                fprintf(stderr, "Invalid visual wall thickness '%s'.\n", value);
                exit(EXIT_FAILURE);
            }
            cli_override_wall_thickness_vis_sigma = v;
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
        } else if (strncmp(arg, "--hb-temp", 9) == 0) {
            // ##CHRIS: Heat bath temperature override
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0; char *endptr = NULL; float v = strtof(value, &endptr);
            if (errno != 0 || endptr == value || *endptr != '\0' || v <= 0.0f) {
                fprintf(stderr, "Invalid heat bath temperature '%s'.\n", value);
                exit(EXIT_FAILURE);
            }
            cli_override_hb_temp = v;
            heatbath_enabled = 1;  // Auto-enable heat bath when temperature specified
        } else if (strncmp(arg, "--timescale", 11) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            errno = 0; char *endptr = NULL; float v = strtof(value, &endptr);
            if (errno != 0 || endptr == value || *endptr != '\0' || v <= 0.0f) {
                fprintf(stderr, "Invalid timescale '%s'.\n", value);
                exit(EXIT_FAILURE);
            }
            time_scale_runtime = v;
        } else if (strcmp(arg, "--edmd-acc") == 0) {
            cli_edmd_acc = 1;
        } else if (strncmp(arg, "--edmd-acc", strlen("--edmd-acc")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            if (!value || !value[0]) {
                cli_edmd_acc = 1;
            } else if (strcmp(value, "1") == 0 || strcmp(value, "true") == 0 || strcmp(value, "yes") == 0 || strcmp(value, "on") == 0) {
                cli_edmd_acc = 1;
            } else if (strcmp(value, "0") == 0 || strcmp(value, "false") == 0 || strcmp(value, "no") == 0 || strcmp(value, "off") == 0) {
                cli_edmd_acc = 0;
            } else {
                fprintf(stderr, "Invalid edmd-acc '%s'. Use 0|1 (or true|false).\n", value);
                exit(EXIT_FAILURE);
            }
        } else if (strncmp(arg, "--mode", 6) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            if (strcmp(value, "time") == 0) sim_mode = MODE_TIME;
            else if (strcmp(value, "rk4") == 0) sim_mode = MODE_RK4;
            else if (strcmp(value, "edmd") == 0) sim_mode = MODE_EDMD;
            else if (strcmp(value, "edmd-hybrid") == 0 ||
                     strcmp(value, "hybrid") == 0 ||
                     strcmp(value, "event-split") == 0 ||
                     strcmp(value, "event_split") == 0) sim_mode = MODE_EDMD_HYBRID;
            else { fprintf(stderr, "Unknown mode '%s'. Use: time, rk4, edmd, edmd-hybrid\n", value); exit(EXIT_FAILURE);} 
        } else if (strncmp(arg, "--colourcode", 12) == 0 || strncmp(arg, "--colorcode", 11) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            if (strcmp(value, "none") == 0) colorcode_mode = COLORCODE_NONE;
            else if (strcmp(value, "temperature") == 0) colorcode_mode = COLORCODE_TEMPERATURE;
            else if (strcmp(value, "velocity") == 0) colorcode_mode = COLORCODE_VELOCITY;
            else { fprintf(stderr, "Unknown colourcode '%s'. Use: none | temperature | velocity\n", value); exit(EXIT_FAILURE);} 
        } else if (strncmp(arg, "--memory-view", strlen("--memory-view")) == 0 ||
                   strncmp(arg, "--memory", strlen("--memory")) == 0 ||
                   strncmp(arg, "--szilard-view", strlen("--szilard-view")) == 0) {
            const char *value = cli_option_value(arg, argc, argv, &i);
            if (strcmp(value, "species") == 0) sz_view_mode = SZ_VIEW_SPECIES;
            else if (strcmp(value, "memory") == 0) sz_view_mode = SZ_VIEW_MEMORY;
            else if (strcmp(value, "overlay") == 0) sz_view_mode = SZ_VIEW_OVERLAY;
            else {
                fprintf(stderr, "Unknown memory-view '%s'. Use: species | memory | overlay\n", value);
                exit(EXIT_FAILURE);
            }
        } else {
            fprintf(stderr, "Unknown option '%s'.\n", arg);
            print_cli_usage(argv[0]);
            exit(EXIT_FAILURE);
        }
    }

    // Finalize adaptive substeps parameters (after parsing all options).
    if (cli_substeps_auto) {
        if (cli_substeps_auto_min < 1) cli_substeps_auto_min = 1;
        if (cli_substeps_auto_max < 1) cli_substeps_auto_max = 1;
        if (cli_substeps_auto_max < cli_substeps_auto_min) {
            if (!cli_quiet) {
                fprintf(stderr, "⚠️ substeps-auto-max (%d) < substeps-auto-min (%d); clamping max=min.\n",
                        cli_substeps_auto_max, cli_substeps_auto_min);
            }
            cli_substeps_auto_max = cli_substeps_auto_min;
        }
        if (!isfinite(cli_substeps_auto_dx_frac) || cli_substeps_auto_dx_frac <= 0.0f) {
            cli_substeps_auto_dx_frac = 0.25f;
        } else if (cli_substeps_auto_dx_frac > 1.0f) {
            cli_substeps_auto_dx_frac = 1.0f;
        }
    }

    // Default energy-transfer output mode:
    // Prefer spring by default so energy_transfer behaves consistently across TIME and EDMD.
    // (EDMD implements a harmonic divider via divider_k/divider_xeq; no timestep impulses required.)
    // Users can still choose --eff-output=wall-ke for a purely hard-body load proxy.
    if (!cli_eff_output_set) {
        cli_eff_output_mode = EFF_OUTPUT_SPRING;
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
        case EXPERIMENT_PRESET_PARTICLELIFE:
            // ##CHRIS: Particle-life needs no walls/dividers
            cli_requested_walls = 0;
            num_internal_walls = 0;
            break;
        case EXPERIMENT_PRESET_SIMPLE_BOX_PREDICTION:
            cli_requested_walls = 0;
            num_internal_walls = 0;
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
        int Nreq = (cli_override_particles > 0) ? cli_override_particles : particles_active;
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
    
    // Apply multi-wall mass factors ONLY if explicitly provided via CLI
    if (cli_wall_masses_set && cli_num_walls > 0 && cli_wall_masses[0] > 0.0f) {
        wall_mass_runtime = PARTICLE_MASS * cli_wall_masses[0];
    }
    if (cli_override_wall_thickness_sigma > 0.0f) {
        set_wall_thickness_sigma(cli_override_wall_thickness_sigma);
    }
    if (cli_override_wall_thickness_vis_sigma > 0.0f) {
        wall_thickness_visual_runtime = cli_override_wall_thickness_vis_sigma * PIXELS_PER_SIGMA;
    }
    if (cli_override_temperature > 0.0f) {
        temperature_runtime = cli_override_temperature;
    }
    // ##CHRIS: Apply heat bath temperature override
    if (cli_override_hb_temp > 0.0f) {
        heatbath_temperature = cli_override_hb_temp;
    }
    // Apply particle size overrides (sigma units)
    if (cli_packing_fraction_set) {
        // ##CHRIS: Auto-calculate particle radius from desired packing fraction.
        // Default: global eta based on total area and total N.
        // Improved: if per-compartment counts are provided, honor each segment area
        // and pick the limiting (smallest) radius so all segments satisfy eta.

        int segments = num_internal_walls + 1;
        bool used_segment_eta = false;
        float r_sigma = 0.0f;

        // Build per-segment counts when we have explicit counts.
        // Priority: exact list (--particles-boxes) or 2-box shorthand (--particles-box-left/right).
        if (segments > 0 && cli_particles_box_counts_count == (size_t)segments) {
            // widths from wall positions (sigma). If none provided, assume equal widths.
            float widths[MAX_DIVIDER_CAPACITY + 1];
            // Left boundary at 0, right at 2*L0_UNITS
            float last = 0.0f;
            for (int s = 0; s < segments; ++s) {
                float right = (s < num_internal_walls) ? cli_wall_positions[s] : (2.0f * L0_UNITS);
                widths[s] = fmaxf(0.0f, right - last);
                last = right;
            }
            float height = HEIGHT_UNITS;

            float r_min = 1e9f;
            for (int s = 0; s < segments; ++s) {
                int cnt = cli_particles_box_counts[s];
                if (cnt <= 0) continue; // skip empty
                float area = widths[s] * height;
                if (area <= 0.0f) continue;
                float r_sq = (cli_packing_fraction * area) / (cnt * (float)M_PI);
                if (r_sq > 0.0f) {
                    float rseg = sqrtf(r_sq);
                    if (rseg < r_min) r_min = rseg;
                }
            }
            if (r_min < 1e8f) {
                r_sigma = r_min;
                used_segment_eta = true;
                if (!cli_quiet) {
                    printf("✅ [Packing Fraction] eta=%.4f → r=%.4fσ (limiting segment) using per-box counts.\n",
                           cli_packing_fraction, r_sigma);
                }
            }
        } else if (segments == 2 && (cli_particles_box_left >= 0 || cli_particles_box_right >= 0)) {
            // Two-box shorthand
            int left = cli_particles_box_left >= 0 ? cli_particles_box_left : 0;
            int right = cli_particles_box_right >= 0 ? cli_particles_box_right : 0;
            int counts[2] = {left, right};
            float widths[2];
            widths[0] = (num_internal_walls > 0) ? cli_wall_positions[0] : L0_UNITS;
            widths[1] = (num_internal_walls > 0) ? (2.0f * L0_UNITS - cli_wall_positions[0]) : L0_UNITS;
            float height = HEIGHT_UNITS;
            float r_min = 1e9f;
            for (int s = 0; s < 2; ++s) {
                int cnt = counts[s];
                if (cnt <= 0) continue;
                float area = widths[s] * height;
                if (area <= 0.0f) continue;
                float r_sq = (cli_packing_fraction * area) / (cnt * (float)M_PI);
                if (r_sq > 0.0f) {
                    float rseg = sqrtf(r_sq);
                    if (rseg < r_min) r_min = rseg;
                }
            }
            if (r_min < 1e8f) {
                r_sigma = r_min;
                used_segment_eta = true;
                if (!cli_quiet) {
                    printf("✅ [Packing Fraction] eta=%.4f → r=%.4fσ (limiting segment) using two-box counts.\n",
                           cli_packing_fraction, r_sigma);
                }
            }
        }

        if (!used_segment_eta) {
            // Fallback to global computation
            float box_area = 2.0f * L0_UNITS * HEIGHT_UNITS;  // Total box area in σ² units
            int N = (cli_override_particles > 0) ? cli_override_particles : particles_active;
            float r_squared = (cli_packing_fraction * box_area) / (N * M_PI);
            r_sigma = sqrtf(r_squared);  // Particle radius in σ units
            if (!cli_quiet) {
                printf("✅ [Packing Fraction] eta=%.4f → r=%.4fσ (d=%.4fσ) for N=%d, Area=%.1fσ²\n",
                       cli_packing_fraction, r_sigma, 2.0f * r_sigma, N, box_area);
            }
        }

        particle_scale_runtime = r_sigma / PARTICLE_RADIUS_UNIT;
    } else if (cli_particle_radius_set) {
        // PARTICLE_RADIUS = PARTICLE_RADIUS_UNIT * PIXELS_PER_SIGMA * particle_scale_runtime
        // We want: PARTICLE_RADIUS = cli_particle_radius_sigma * PIXELS_PER_SIGMA
        // => particle_scale_runtime = cli_particle_radius_sigma / PARTICLE_RADIUS_UNIT
        particle_scale_runtime = cli_particle_radius_sigma / PARTICLE_RADIUS_UNIT;
    } else if (cli_particle_diameter_set) {
        // diameter provided → radius = diameter/2
        float r_sigma = 0.5f * cli_particle_diameter_sigma;
        particle_scale_runtime = r_sigma / PARTICLE_RADIUS_UNIT;
    }
    // Auto-raise phi_max when eta is explicitly set so capacity checks don't
    // reject valid configurations. The capacity check uses effective area
    // (reduced by particle radius margins) while eta uses raw segment area,
    // so we need generous headroom to avoid false overpacking errors.
    if (cli_packing_fraction_set) {
        float needed_phi_max = fminf(0.9069f, cli_packing_fraction * 1.25f + 0.05f);
        if (needed_phi_max > cli_phi_max) {
            if (!cli_quiet) {
                printf("[Packing] Auto-raising phi_max from %.3f to %.3f to accommodate eta=%.4f\n",
                       (double)cli_phi_max, (double)needed_phi_max, (double)cli_packing_fraction);
            }
            cli_phi_max = needed_phi_max;
        }
    }

    if (cli_force_no_experiments) {
        enable_speed_of_sound_experiments = false;
    }

    // Energy/efficiency measurement is initialized later (initialize_energy_measurement()).
    // Do not force-enable the spring here: in EDMD we often want --eff-output=wall-ke (no spring forces).
    if (cli_enable_energy_measurement && !cli_quiet) {
        printf("🔬 Energy/efficiency measurement requested\n");
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
    if (cli_quiet) return;
    printf("\n========== Simulation Info ==========\n");
    printf("Mode: %s\n", MOLECULAR_MODE ? "MOLECULAR (Normalized units)" : "MACROSCOPIC (Real units)");
    printf("Experiment Mode: %s\n", EXPERIMENT_MODE ? "ON (fine dt/substeps)" : "OFF (default dt)");
    printf("Fixed Time Step (dt): %.2e\n", FIXED_DT);
    printf("Runtime Time Step (dt_runtime): %.2e\n", fixed_dt_runtime);
    if (cli_fixed_dt_set) {
        printf("Fixed dt override: %.2e\n", cli_fixed_dt);
    }
    if (cli_substeps_override > 0) {
        printf("Collision substeps: %d (override)\n", cli_substeps_override);
    } else if (cli_substeps_auto) {
        printf("Collision substeps: auto (min=%d, max=%d, dx<=%.2f*diam)\n",
               cli_substeps_auto_min, cli_substeps_auto_max, (double)cli_substeps_auto_dx_frac);
    } else {
        printf("Collision substeps: %d\n", SUBSTEPS);
    }
    if (cli_pp_toi) {
        printf("PP-TOI: enabled (max_events=%d, neighbor_ring=%d)\n",
               cli_pp_toi_max_events, cli_pp_toi_neighbor_ring);
    }
    if (sim_mode == MODE_EDMD || sim_mode == MODE_EDMD_HYBRID) {
        printf("EDMD backend: %s\n", cli_edmd_acc ? "accelerated" : "default");
    }
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
    printf("Wall Mass Factor:%f\n", wall_mass_runtime);

    // Steps/time target (time mode): default 3000 steps unless overridden by --steps
    int max_steps = (cli_override_num_steps > 0) ? cli_override_num_steps : 3000;
    double target_time_sigma = (fixed_dt_runtime / (MOLECULAR_MODE ? 1.0f : TAU_time_scale_factor_to_molecular)) * max_steps;
    printf("Target total simulation time: %.3f\n", target_time_sigma);
    printf("Calculated steps to reach target: %d\n", max_steps);
    printf("======================================\n\n");
}

//////////////// SDL INITIALIZATION //////////////////////
// Function to Initialize SDL
void initSDL() {
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        printf("SDL could not initialize! SDL_Error: %s\n", SDL_GetError());
        exit(1);
    }
    // Pre-clamp against the primary display so the initial window isn't absurdly large.
    // (We do a second clamp after window creation using the window's actual display index.)
    {
        SDL_DisplayMode dm;
        if (SDL_GetCurrentDisplayMode(0, &dm) == 0) {
            const int max_w = (int)(0.98 * (double)dm.w);
            const int max_h = (int)(0.95 * (double)dm.h);

            // If we have room, expand width to use the screen for HUD/plots.
            // Do NOT auto-expand height to full-screen: users often want vertical space for other apps.
            if (gui_dist_mode != 0) {
                if (MAX_X < max_w) MAX_X = max_w;
            }

            // Clamp to fit display.
            if (MAX_X > max_w) MAX_X = max_w;
            if (MAX_Y > max_h) MAX_Y = max_h;

            // Ensure some bottom area remains for plots.
            const int min_h = YW2 + gui_window_min_extra_y;
            if (MAX_Y < min_h) MAX_Y = min_h;

            gui_recompute_hist_extents();
        }
    }
    _window = SDL_CreateWindow(
        "Hard Spheres Simulation",
        SDL_WINDOWPOS_UNDEFINED,
        SDL_WINDOWPOS_UNDEFINED,
        MAX_X,
        MAX_Y,
        SDL_WINDOW_SHOWN | SDL_WINDOW_RESIZABLE | SDL_WINDOW_ALLOW_HIGHDPI
    );
    if (!_window) {
        printf("Window could not be created! SDL_Error: %s\n", SDL_GetError());
        exit(1);
    }
    _renderer = SDL_CreateRenderer(_window, -1, SDL_RENDERER_ACCELERATED);
    if (!_renderer) {
        printf("Renderer could not be created! SDL_Error: %s\n", SDL_GetError());
        exit(1);
    }

    // Post-clamp to the actual display the window ended up on (important for external monitors).
    {
        int display_index = SDL_GetWindowDisplayIndex(_window);
        if (display_index < 0) display_index = 0;
        SDL_DisplayMode dm;
        if (SDL_GetCurrentDisplayMode(display_index, &dm) == 0) {
            const int max_w = (int)(0.98 * (double)dm.w);
            const int max_h = (int)(0.95 * (double)dm.h);
            if (gui_dist_mode != 0 && MAX_X < max_w) MAX_X = max_w;
            if (MAX_X > max_w) MAX_X = max_w;
            if (MAX_Y > max_h) MAX_Y = max_h;
            const int min_h = YW2 + gui_window_min_extra_y;
            if (MAX_Y < min_h) MAX_Y = min_h;
            SDL_SetWindowSize(_window, MAX_X, MAX_Y);
            gui_sync_renderer_canvas_size();
        }
    }
    gui_sync_renderer_canvas_size();
}

static void gui_adjust_window_height(int delta) {
    if (!_window) return;

    int w = 0, h = 0;
    SDL_GetWindowSize(_window, &w, &h);
    int desired_h = h + delta;

    int display_index = SDL_GetWindowDisplayIndex(_window);
    if (display_index < 0) display_index = 0;

    int max_h = desired_h;
    SDL_DisplayMode dm;
    if (SDL_GetCurrentDisplayMode(display_index, &dm) == 0) {
        max_h = (int)(0.95 * (double)dm.h);
    }

    const int min_h = YW2 + gui_window_min_extra_y;
    if (desired_h < min_h) desired_h = min_h;
    if (desired_h > max_h) desired_h = max_h;

    SDL_SetWindowSize(_window, w, desired_h);
    gui_sync_renderer_canvas_size();

    if (!cli_quiet) {
        printf("📐 Window height: %d (HIST_HEIGHT=%d) (PgUp/PgDn)\n", MAX_Y, HIST_HEIGHT);
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

    int Nreq = (cli_override_particles > 0) ? cli_override_particles : particles_active;
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
    if (!cli_packing_fraction_set && !cli_particle_radius_set && Nreq > cap_sum) {
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

    // Absolute counts for segmented seeding
    if (!applied_absolute_two_box && preset_custom_absolute_counts && (int)cli_particles_box_counts_count == segments) {
        int assigned = 0;
        for (int seg = 0; seg < segments; ++seg) {
            int c = cli_particles_box_counts[seg];
            if (c < 0) c = 0;
            if (c > caps[seg]) c = caps[seg];
            target_counts[seg] = c;
            assigned += c;
        }
        // If user requested more than capacity, clamp; if less, leave empty.
        applied_absolute_two_box = true;
        (void)assigned;
    } else if (!applied_absolute_two_box && preset_custom_fractions && preset_wall_count + 1 == segments) {
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
    int Nreq = (cli_override_particles > 0) ? cli_override_particles : particles_active;
    if (Nreq < 0) Nreq = 0; if (Nreq > NUM_PARTICLES) Nreq = NUM_PARTICLES;
    float height_eff = (YW2 - YW1) - 2.0f * PARTICLE_RADIUS;
    float area_p = (float)(M_PI * PARTICLE_RADIUS * PARTICLE_RADIUS);
    int caps[MAX_DIVIDER_CAPACITY + 1]; int cap_sum = 0;
    for (int s = 0; s < segments; ++s) {
        float w = fmaxf(0.0f, (seg_right[s] - seg_left[s]) - 2.0f * PARTICLE_RADIUS);
        int cap = (int)floorf(cli_phi_max * (w * height_eff) / area_p);
        if (cap < 0) cap = 0; caps[s] = cap; cap_sum += cap;
    }
    if (!cli_packing_fraction_set && !cli_particle_radius_set && Nreq > cap_sum) {
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
            if (!cli_packing_fraction_set && !cli_particle_radius_set && c > caps[s]) {
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
        if (!cli_packing_fraction_set && !cli_particle_radius_set && leftover > 0) {
            fprintf(stderr, "❌ Overpacked (left-empty): requested N exceeds capacity of non-left segments.\n");
            exit(1);
        }
    }

    // Place on an automatic lattice per segment:
    // low/moderate density -> rectangular grid; dense compartments -> hexagonal close lattice.
    int index = 0;
    const float diameter_px = 2.0f * (float)PARTICLE_RADIUS;
    const float seed_pad_px = fmaxf(0.05f, 0.005f * diameter_px);
    const float margin = (float)PARTICLE_RADIUS + seed_pad_px; // keep centers away from boundaries
    const float pistonL_plane = piston_left_x + 5.0f;          // matches EDMD pistonL_x convention
    const float pistonR_plane = piston_right_x;                // piston plane at x = piston_right_x
    for (int s = 0; s < segments; ++s) {
        int target = target_counts[s];
        if (target <= 0) continue;
        float left_bound  = seg_left[s]  + margin;
        float right_bound = seg_right[s] - margin;
        // Also respect piston planes if they sit inside this segment.
        if (pistonL_plane > seg_left[s] && pistonL_plane < seg_right[s]) {
            left_bound = fmaxf(left_bound, pistonL_plane + margin);
        }
        if (pistonR_plane > seg_left[s] && pistonR_plane < seg_right[s]) {
            right_bound = fminf(right_bound, pistonR_plane - margin);
        }
        float top = YW1 + margin;
        float bottom = YW2 - margin;
        float w = fmaxf(0.0f, right_bound - left_bound);
        float h = fmaxf(0.0f, bottom - top);
        if (w <= 0.0f || h <= 0.0f) continue;

        const float min_gap = 2.0f * (float)PARTICLE_RADIUS * 1.0005f;
        const double r_sigma = (double)PARTICLE_RADIUS / (double)PIXELS_PER_SIGMA;
        const double seg_area_sigma = fmax(1e-12,
                                           ((double)(seg_right[s] - seg_left[s]) / (double)PIXELS_PER_SIGMA) *
                                           (double)HEIGHT_UNITS);
        const double local_eta = (double)target * M_PI * r_sigma * r_sigma / seg_area_sigma;

        int cols = (int)ceil(sqrt((double)target * (double)w / fmax(1.0, (double)h)));
        if (cols < 1) cols = 1;
        int rows = (target + cols - 1) / cols;
        if (rows < 1) rows = 1;
        float dx = w / (float)cols;
        float dy = h / (float)rows;
        const bool rectangular_ok = (dx >= min_gap && dy >= min_gap && local_eta < 0.30);

        int placed = 0;
        if (rectangular_ok) {
            for (int r = 0; r < rows && placed < target; ++r) {
                for (int c = 0; c < cols && placed < target; ++c) {
                    if (index >= particles_active) break;
                    const float cx = left_bound + (c + 0.5f) * dx;
                    const float cy = top + (r + 0.5f) * dy;
                    X[index] = cx;
                    Y[index] = cy;
                    Vx[index] = maxwell_boltzmann_velocity_gaussians(temperature_runtime);
                    Vy[index] = maxwell_boltzmann_velocity_gaussians(temperature_runtime);
                    Radius[index] = PARTICLE_RADIUS;
                    index++;
                    placed++;
                }
            }
        } else {
            const float hex_dx = min_gap;
            const float hex_dy = min_gap * 0.8660254038f;
            int max_rows = (hex_dy > 0.0f) ? ((int)floorf(h / hex_dy) + 1) : 1;
            if (max_rows < 1) max_rows = 1;
            for (int r = 0; r < max_rows && placed < target; ++r) {
                const float y = top + (float)r * hex_dy;
                if (y > bottom + 1e-3f) break;
                const float offset = (r & 1) ? 0.5f * hex_dx : 0.0f;
                int row_cols = (hex_dx > 0.0f) ? ((int)floorf((w - offset) / hex_dx) + 1) : 1;
                if (row_cols < 1) continue;
                const float row_width = (float)(row_cols - 1) * hex_dx;
                float x0 = left_bound + 0.5f * fmaxf(0.0f, w - row_width);
                if (r & 1) x0 += 0.25f * hex_dx;
                while (x0 + row_width > right_bound + 1e-3f && x0 > left_bound) x0 -= 0.25f * hex_dx;
                if (x0 < left_bound) x0 = left_bound;
                for (int c = 0; c < row_cols && placed < target; ++c) {
                    const float cx = x0 + (float)c * hex_dx;
                    if (cx < left_bound - 1e-3f || cx > right_bound + 1e-3f) continue;
                    if (index >= particles_active) break;
                    X[index] = cx;
                    Y[index] = y;
                    Vx[index] = maxwell_boltzmann_velocity_gaussians(temperature_runtime);
                    Vy[index] = maxwell_boltzmann_velocity_gaussians(temperature_runtime);
                    Radius[index] = PARTICLE_RADIUS;
                    index++;
                    placed++;
                }
            }
        }

        if (placed < target) {
            fprintf(stderr,
                    "❌ Could only place %d/%d particles in segment %d. "
                    "Try lower eta/radius, fewer particles, fewer walls, or a larger box.\n",
                    placed, target, s);
            exit(1);
        }
    }

    // ##CHRIS: Velocity rescale to match target temperature exactly
    // CRITICAL: Must use kB_effective() not K_B to respect --kbt1 flag!
    double actual_ke = kinetic_energy();
    double target_ke = particles_active * kB_effective() * temperature_runtime;
    if (actual_ke > 0.0) {
        double scale = sqrt(target_ke / actual_ke);
        for (int i = 0; i < particles_active; ++i) {
            Vx[i] *= scale;
            Vy[i] *= scale;
        }
        if (!cli_quiet) {
            printf("✅ Velocity rescaling: actual_ke=%.3f, target_ke=%.3f (kBT=%.3f), scale=%.4f\n",
                   actual_ke, target_ke, kB_effective() * temperature_runtime, scale);
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

    int Nreq = (cli_override_particles > 0) ? cli_override_particles : particles_active;
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

            if (!cli_quiet) printf("✅ Particle %d (LEFT) placed at (%.2f, %.2f)\n", idx, x, y);
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

            if (!cli_quiet) printf("✅ Particle %d (RIGHT) placed at (%.2f, %.2f)\n", idx, x, y);
            idx++;
        }
    }

    if (!cli_quiet) {
        printf("✅ Finished placing all %d particles evenly across box!\n", particles_active);
        // Diagnostics: capacity and targets
        printf("[SEED] MAX(buffer)=%d requested=%d active=%d phi_max=%.3f\n",
               NUM_PARTICLES,
               (cli_override_particles > 0 ? cli_override_particles : particles_active),
               particles_active, cli_phi_max);
        printf("[SEED] caps: L=%d R=%d  targets: L=%d R=%d\n", cap_left, cap_right, particles_left, particles_right);
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
    if (!cli_quiet) {
        printf("🔍 Initial particle velocities (first 10):\n");
        for (int i = 0; i < 10; i++) {
            printf("Particle %2d: Vx = %.6f, Vy = %.6f\n", i, Vx[i], Vy[i]);    
        }
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
    // NOTE: positions X/Y are in pixels. DIAMETER is in σ-units (based on PARTICLE_RADIUS_UNIT),
    // so we must derive margins in pixels from PARTICLE_RADIUS instead of using DIAMETER directly.
    const float diameter_px = 2.0f * (float)PARTICLE_RADIUS;
    const float seed_pad_px = fmaxf(1.0f, 0.10f * diameter_px);   // tiny extra padding beyond geometric R
    const float margin_px = 1.25f * diameter_px + seed_pad_px;    // keep centers away from boundaries/walls

    const float wall_half_thick_px = 0.5f * (float)WALL_THICKNESS;
    const float pistonL_plane = piston_left_x + 5.0f; // matches EDMD pistonL_x convention
    const float pistonR_plane = piston_right_x;       // piston plane at x = piston_right_x
    float left_x_min = (float)XW1 + margin_px;
    const float left_x_max = wall_x - wall_half_thick_px - margin_px;
    const float right_x_min = wall_x + wall_half_thick_px + margin_px;
    float right_x_max = (float)XW2 - margin_px;
    // If piston planes sit inside the box, keep seed points away from them too.
    if (pistonL_plane > (float)XW1 && pistonL_plane < wall_x) {
        left_x_min = fmaxf(left_x_min, pistonL_plane + margin_px);
    }
    if (pistonR_plane > wall_x && pistonR_plane < (float)XW2) {
        right_x_max = fminf(right_x_max, pistonR_plane - margin_px);
    }

    const float y_min = (float)YW1 + margin_px;
    const float y_max = (float)YW2 - margin_px;

    const float left_width = fmaxf(0.0f, left_x_max - left_x_min);
    const float right_width = fmaxf(0.0f, right_x_max - right_x_min);
    const float height = fmaxf(0.0f, y_max - y_min);

    // Compute per-box capacity (area-based) - FIXED to use L0_UNITS directly
    // The old calculation was too conservative due to pixel/sigma confusion
    float height_eff = HEIGHT_UNITS - DIAMETER;  // Leave 1 diameter margin (top+bottom)
    float left_eff = L0_UNITS - DIAMETER;        // Leave 1 diameter margin (wall+boundary)
    float right_eff = L0_UNITS - DIAMETER;       // Same for right side
    // ##CHRIS: Use actual particle radius (accounts for --eta, --particle-radius, etc.)
    float radius_sigma = PARTICLE_RADIUS / PIXELS_PER_SIGMA;  // Convert to σ units
    float area_p = (float)(M_PI * radius_sigma * radius_sigma);
    int cap_left = (int)floorf(cli_phi_max * (left_eff * height_eff) / area_p);
    int cap_right = (int)floorf(cli_phi_max * (right_eff * height_eff) / area_p);

    int requested_total = (cli_override_particles > 0) ? cli_override_particles : particles_active;
    int req_left = (cli_particles_box_left >= 0) ? cli_particles_box_left : -1;
    int req_right = (cli_particles_box_right >= 0) ? cli_particles_box_right : -1;

    int particles_left = 0, particles_right = 0;
    // Absolute counts for two segments via --particles-boxes counts (length==2)
    if (preset_custom_absolute_counts && cli_particles_box_counts_count == 2) {
        int left = cli_particles_box_counts[0]; if (left < 0) left = 0;
        int right = cli_particles_box_counts[1]; if (right < 0) right = 0;
        // When eta is explicitly set, the particle radius was already computed to
        // achieve the target packing fraction per segment. The capacity formula
        // uses reduced effective area (boundary margins) which is more conservative
        // than the eta formula, so skip the cap check to avoid false rejections.
        if (!cli_packing_fraction_set && !cli_particle_radius_set) {
            if (left > cap_left) {
                printf("❌ Overpacked (left): requested %d > cap %d.\n", left, cap_left); exit(1);
            }
            if (right > cap_right) {
                printf("❌ Overpacked (right): requested %d > cap %d.\n", right, cap_right); exit(1);
            }
        }
        particles_left = left; particles_right = right; particles_active = left + right;
    } else
    if (req_left >= 0 && req_right >= 0) {
        particles_active = req_left + req_right;
        if (!cli_packing_fraction_set && !cli_particle_radius_set &&
            (particles_active > cap_left + cap_right || req_left > cap_left || req_right > cap_right)) {
            printf("❌ Overpacked: requested total=%d (L=%d,R=%d) > capacity total=%d (L_cap=%d,R_cap=%d).\n",
                   particles_active, req_left, req_right, cap_left + cap_right, cap_left, cap_right);
            exit(1);
        }
        particles_left = req_left;
        particles_right = req_right;
    } else {
        if (!cli_packing_fraction_set && !cli_particle_radius_set && requested_total > cap_left + cap_right) {
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
            if (!cli_packing_fraction_set && !cli_particle_radius_set && particles_left > cap_left) particles_left = cap_left;
        } else if (cli_left_empty) {
            particles_left = 0;
            particles_right = particles_active;
            if (!cli_packing_fraction_set && !cli_particle_radius_set && particles_right > cap_right) {
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
        if (!cli_packing_fraction_set && !cli_particle_radius_set && (particles_left > cap_left || particles_right > cap_right)) {
            printf("❌ Per-box overpack: targets L=%d,R=%d exceed caps L_cap=%d,R_cap=%d.\n",
                   particles_left, particles_right, cap_left, cap_right);
            exit(1);
        }
    }

    int num_rows_left = (particles_left > 0) ? (int)sqrt((double)particles_left) : 1;
    if (num_rows_left < 1) num_rows_left = 1;
    int num_cols_left = (particles_left > 0) ? ((particles_left + num_rows_left - 1) / num_rows_left) : 1;
    if (num_cols_left < 1) num_cols_left = 1;

    int num_rows_right = (particles_right > 0) ? (int)sqrt((double)particles_right) : 1;
    if (num_rows_right < 1) num_rows_right = 1;
    int num_cols_right = (particles_right > 0) ? ((particles_right + num_rows_right - 1) / num_rows_right) : 1;
    if (num_cols_right < 1) num_cols_right = 1;

    const float spacing_left_x = (num_cols_left > 0) ? (left_width / (float)num_cols_left) : 0.0f;
    const float spacing_left_y = (num_rows_left > 0) ? (height / (float)num_rows_left) : 0.0f;

    const float spacing_right_x = (num_cols_right > 0) ? (right_width / (float)num_cols_right) : 0.0f;
    const float spacing_right_y = (num_rows_right > 0) ? (height / (float)num_rows_right) : 0.0f;

    int idx = 0;

    // Left half
    for (int row = 0; row < num_rows_left; row++) {
        for (int col = 0; col < num_cols_left; col++) {
            if (idx >= particles_left) break;
            // Center within each grid cell to avoid starting exactly on a wall.
            X[idx] = left_x_min + (col + 0.5f) * spacing_left_x;
            Y[idx] = y_min + (row + 0.5f) * spacing_left_y;

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
            X[idx] = right_x_min + (col + 0.5f) * spacing_right_x;
            Y[idx] = y_min + (row + 0.5f) * spacing_right_y;
            float vx, vy;
            maxwell_boltzmann_2D(temperature_runtime, &vx, &vy);
            Vx[idx] = vx;
            Vy[idx] = vy;
            V_init[idx] = sqrtf(vx * vx + vy * vy);  // Optional: store the speed if needed
            Radius[idx] = PARTICLE_RADIUS;
            idx++;
        }
    }

    if (!cli_quiet) {
        printf("✅ Particles placed with buffer. Wall_x = %.3f\n", wall_x);
        fflush(stdout);
    }
    // After initializing all Vx, Vy:
    double actual_ke = kinetic_energy();  // sum of 0.5*m*v^2
    double target_ke = particles_active * K_B * temperature_runtime;
    double scale = sqrt(target_ke / actual_ke);  // to adjust v_i → v_i * scale
    if (!cli_quiet) {
        printf("✅ Scaling factor for velocities: %.4f\n", scale);
    }

    for (int i = 0; i < particles_active; i++) {
        Vx[i] *= scale;
        Vy[i] *= scale;
}

    // Diagnostics and HUD seed counts
    if (!cli_quiet) {
        printf("[SEED] MAX(buffer)=%d requested=%d active=%d phi_max=%.3f\n",
               NUM_PARTICLES, requested_total, particles_active, cli_phi_max);
        printf("[SEED] caps: L=%d R=%d  targets: L=%d R=%d\n", cap_left, cap_right, particles_left, particles_right);
    }

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

    // Persist run parameters early so analysis tools read the exact settings used
    write_run_params_json();



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

    if (piston_step_active &&
        ((piston_step_direction < 0.0f && piston_right_x <= piston_step_target) ||
         (piston_step_direction > 0.0f && piston_right_x >= piston_step_target))) {
        piston_right_x = piston_step_target;
        vx_piston_right = 0.0f;
        piston_step_active = false;
        piston_protocol_start_time = -1.0f;
        // Start windowed spring response measurement on the next energy update.
        if (energy_measurement.spring_enabled) {
            energy_measurement.eff_window_pending_start = true;
            energy_measurement.eff_window_started = false;
            energy_measurement.eff_window_active = false;
            energy_measurement.eff_window_done = false;
        }
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

// Normalized logistic "sigmoid" velocity profile: v(t) transitions from v0 → vmax over [0,duration].
static float sigmoid_velocity_profile(float t, float v0, float vmax, float duration, float steepness) {
    if (duration <= 1e-9f) return vmax;
    if (t <= 0.0f) return v0;
    if (t >= duration) return vmax;

    // Logistic: s(t)=1/(1+exp(-k*(t-duration/2))). Normalize so s(0)=0, s(duration)=1.
    float k = (steepness > 0.0f) ? (steepness / duration) : (8.0f / duration);
    float s0 = 1.0f / (1.0f + expf(-k * (0.0f - duration * 0.5f)));
    float s1 = 1.0f / (1.0f + expf(-k * (duration - duration * 0.5f)));
    float st = 1.0f / (1.0f + expf(-k * (t - duration * 0.5f)));
    float denom = (s1 - s0);
    float s = (fabsf(denom) > 1e-12f) ? ((st - s0) / denom) : 1.0f;
    if (s < 0.0f) s = 0.0f;
    if (s > 1.0f) s = 1.0f;
    return v0 + (vmax - v0) * s;
}

// Protocol-based piston movement
void apply_piston_protocol(float dt) {
    if (piston_step_active) {
        if (piston_protocol_start_time < 0.0f) {
            piston_protocol_start_time = simulation_time;
        }
        
        float elapsed = simulation_time - piston_protocol_start_time;
        float duration = (piston_protocol_duration > 1e-9f) ? piston_protocol_duration : 1e-9f;
        float vmax = fabsf(piston_protocol_max_speed);
        
        switch (cli_protocol) {
            case PROTOCOL_STEP:
                // Constant max speed immediately
                vx_piston_right = piston_step_direction * vmax;
                break;
                
            case PROTOCOL_SIGMOIDAL:
                // Smooth sigmoid acceleration from v0 -> vmax
                vx_piston_right = piston_step_direction * sigmoid_velocity_profile(
                    elapsed,
                    fmaxf(0.0f, piston_protocol_v0),
                    vmax,
                    duration,
                    piston_protocol_sigmoid_steepness
                );
                break;
                
            case PROTOCOL_LINEAR:
                // Linear acceleration: either explicit v(t)=v0+g*t, or ramp 0→vmax over duration.
                if (piston_protocol_linear_use_gradient) {
                    float v = fmaxf(0.0f, piston_protocol_v0 + piston_protocol_linear_gradient * elapsed);
                    if (v > vmax) v = vmax;
                    vx_piston_right = piston_step_direction * v;
                } else {
                    if (elapsed < duration) {
                        vx_piston_right = piston_step_direction * vmax * (elapsed / duration);
                    } else {
                        vx_piston_right = piston_step_direction * vmax;
                    }
                }
                break;
                
            case PROTOCOL_SINUSOIDAL:
                // Sinusoidal acceleration
                if (elapsed < duration) {
                    vx_piston_right = piston_step_direction * vmax * sinf((float)M_PI * elapsed / duration);
                } else {
                    vx_piston_right = piston_step_direction * vmax;
                }
                break;
                
            case PROTOCOL_OPTIMAL:
                // AI-optimized protocol (placeholder for future implementation)
                vx_piston_right = piston_step_direction * sigmoid_velocity_profile(
                    elapsed,
                    fmaxf(0.0f, piston_protocol_v0),
                    vmax,
                    duration,
                    piston_protocol_sigmoid_steepness
                );
                break;
        }
        
        // Stop when target reached
        if ((piston_step_direction < 0.0f && piston_right_x <= piston_step_target) ||
            (piston_step_direction > 0.0f && piston_right_x >= piston_step_target)) {
            vx_piston_right = 0;
            piston_step_active = false;
            piston_protocol_start_time = -1.0f;
            // Start windowed spring response measurement on the next energy update.
            if (energy_measurement.eff_enabled) {
                energy_measurement.eff_window_pending_start = true;
                energy_measurement.eff_window_started = false;
                energy_measurement.eff_window_active = false;
                energy_measurement.eff_window_done = false;
            }
        }
    } else {
        piston_protocol_start_time = -1.0f;
    }
}

// Energy measurement functions
void initialize_energy_measurement() {
    // Hard gate: only allow energy-transfer efficiency metrics if explicitly requested,
    // or via the energy_transfer preset.
    if (!(cli_enable_energy_measurement || cli_experiment_preset == EXPERIMENT_PRESET_ENERGY_TRANSFER)) {
        return;
    }
    if (num_internal_walls > 0) {
        // Default equilibrium: current leftmost wall position (pixels) so spring starts unloaded
        float eq = wall_x;
        if (cli_spring_eq_set) {
            // User override in σ units from left boundary
            eq = XW1 + cli_spring_eq_sigma * PIXELS_PER_SIGMA;
        } else if (all_wall_positions) {
            // If multi-wall, anchor to leftmost wall to avoid pre-load
            float min_pos = all_wall_positions[leftmost_wall_index];
            for (int w = 0; w < num_internal_walls; ++w) {
                float wx = all_wall_positions[w];
                if (wx < min_pos) min_pos = wx;
            }
            eq = min_pos;
        }
        if (cli_spring_k_set) {
            energy_measurement.spring_constant = cli_spring_k;
        }
        energy_measurement.equilibrium_position = eq;
        energy_measurement.eff_enabled = true;
        energy_measurement.eff_output_mode = cli_eff_output_mode;
        // Only enable/actuate a spring when the output mode is spring.
        // In wall-KE mode we want a "pure hard-body" load metric compatible with EDMD.
        energy_measurement.spring_enabled = (cli_eff_output_mode == EFF_OUTPUT_SPRING);
        energy_measurement.spring_force = 0.0f;
        energy_measurement.spring_energy = 0.0f;
        energy_measurement.spring_energy_max = 0.0f;
        energy_measurement.energy_transferred = 0.0f;
        energy_measurement.max_displacement = 0.0f;
        energy_measurement.eff_window_pending_start = false;
        energy_measurement.eff_window_started = false;
        energy_measurement.eff_window_active = false;
        energy_measurement.eff_window_done = false;
        energy_measurement.eff_piston_stop_time = 0.0f;
        energy_measurement.spring_energy_at_stop = 0.0f;
        energy_measurement.spring_energy_peak_window = 0.0f;
        energy_measurement.spring_energy_peak_window_time = 0.0f;
        if (!cli_quiet) {
            const char *mode = (cli_eff_output_mode == EFF_OUTPUT_WALL_KE) ? "wall-ke" : "spring";
            printf("🔬 Energy measurement initialized: mode=%s, equilibrium at %.2f\n",
                   mode, energy_measurement.equilibrium_position);
        }
    }
}

void update_energy_measurement(float dt) {
    if (!energy_measurement.eff_enabled || num_internal_walls <= 0) return;
    
    // Get current position of leftmost wall
    float current_position = all_wall_positions[leftmost_wall_index];
    float displacement = current_position - energy_measurement.equilibrium_position;
    
    // Update maximum displacement
    if (fabsf(displacement) > fabsf(energy_measurement.max_displacement)) {
        energy_measurement.max_displacement = displacement;
    }
    
    // Compute the current "output energy" according to the selected mode.
    // NOTE: We store the chosen output into spring_force/spring_energy fields for backward
    // compatibility with existing CSV schemas and plotting scripts.
    //
    // - spring: SpringE = 0.5*k*x^2, SpringF = -k*x
    // - wall-ke: SpringE = 0.5*M_wall*vx_wall^2, SpringF = 0
    if (energy_measurement.eff_output_mode == EFF_OUTPUT_SPRING) {
        energy_measurement.spring_force = -energy_measurement.spring_constant * displacement;
        energy_measurement.spring_energy = 0.5f * energy_measurement.spring_constant * displacement * displacement;
    } else {
        energy_measurement.spring_force = 0.0f;
        // vx_wall is the velocity of the primary (leftmost) internal wall
        energy_measurement.spring_energy = 0.5f * (float)wall_mass_runtime * (float)(vx_wall * vx_wall);
    }
    if (energy_measurement.spring_energy > energy_measurement.spring_energy_max) {
        energy_measurement.spring_energy_max = energy_measurement.spring_energy;
    }

    // Windowed "gain": start when piston stops, then track peak within a fixed sigma-time window.
    // Use t_now (state after this dt) to align with the updated wall position.
    {
        const float dt_sigma = dt / (float)PIXELS_PER_SIGMA;
        const float t_now = (float)simulation_time + dt_sigma;

        if (energy_measurement.eff_window_pending_start && !energy_measurement.eff_window_started) {
            energy_measurement.eff_window_pending_start = false;
            energy_measurement.eff_window_started = true;
            energy_measurement.eff_window_active = true;
            energy_measurement.eff_window_done = false;
            energy_measurement.eff_piston_stop_time = t_now;
            energy_measurement.spring_energy_at_stop = energy_measurement.spring_energy;
            energy_measurement.spring_energy_peak_window = energy_measurement.spring_energy;
            energy_measurement.spring_energy_peak_window_time = t_now;
        }

        if (energy_measurement.eff_window_active) {
            const float Tw = (cli_eff_window_sigma > 0.0f) ? cli_eff_window_sigma : 0.0f;
            const float dt_since = t_now - energy_measurement.eff_piston_stop_time;
            if (dt_since <= Tw) {
                if (energy_measurement.spring_energy > energy_measurement.spring_energy_peak_window) {
                    energy_measurement.spring_energy_peak_window = energy_measurement.spring_energy;
                    energy_measurement.spring_energy_peak_window_time = t_now;
                }
            } else {
                energy_measurement.eff_window_active = false;
                energy_measurement.eff_window_done = true;
            }
        }
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
           all_wall_positions[leftmost_wall_index], energy_measurement.equilibrium_position);
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

// MB-style speed palette: dark blue -> light blue/cyan -> green -> yellow -> orange/red.
// This is render-only and maps normalized speed |v|/vmax to color.
static inline void gui_speed_heat_rgb(float a01, Uint8 *out_r, Uint8 *out_g, Uint8 *out_b) {
    if (a01 < 0.0f) a01 = 0.0f;
    if (a01 > 1.0f) a01 = 1.0f;

    const float stops[] = {0.00f, 0.18f, 0.36f, 0.55f, 0.74f, 0.88f, 1.00f};
    const Uint8 colors[][3] = {
        {  0,  45, 210},  // slow: dark blue
        {  0, 135, 255},  // blue
        {  0, 235, 245},  // cyan
        { 45, 235,  80},  // green
        {245, 235,  25},  // yellow
        {255, 125,  20},  // orange
        {235,  20,  30}   // fast: red
    };

    int k = 0;
    const int nstops = (int)(sizeof(stops) / sizeof(stops[0]));
    while (k + 1 < nstops && a01 > stops[k + 1]) k++;
    if (k + 1 >= nstops) k = nstops - 2;
    float t = (a01 - stops[k]) / fmaxf(1e-6f, stops[k + 1] - stops[k]);
    if (t < 0.0f) t = 0.0f;
    if (t > 1.0f) t = 1.0f;

    *out_r = (Uint8)lroundf((1.0f - t) * colors[k][0] + t * colors[k + 1][0]);
    *out_g = (Uint8)lroundf((1.0f - t) * colors[k][1] + t * colors[k + 1][1]);
    *out_b = (Uint8)lroundf((1.0f - t) * colors[k][2] + t * colors[k + 1][2]);
}

static inline void set_draw_color_speed_heat(float a01, Uint8 alpha) {
    Uint8 r = 0, g = 0, b = 0;
    gui_speed_heat_rgb(a01, &r, &g, &b);
    SDL_SetRenderDrawColor(_renderer, r, g, b, alpha);
}

// Function to render particles with optional color coding
void render_particles() {
    int n = particles_active > 0 ? particles_active : 0;
    if (n <= 0) return;

    // Defaults: white color
    SDL_SetRenderDrawColor(_renderer, 255, 255, 255, 255);

    float vmin = FLT_MAX, vmax = 0.0f;
    if (gui_mb_speed_view) {
        vmax = gui_mb_speed_vmax_for_render();
        if (!(vmax > 1e-6f)) vmax = 1.0f;
        vmin = 0.0f;
    } else if (colorcode_mode == COLORCODE_VELOCITY) {
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
        // Use rounding (not truncation) to reduce 1–2px visual "overlaps" at contact.
        int pixel_x = (int)lroundf((float)X[i]);
        int pixel_y = (int)lroundf((float)Y[i]);
        int radius_pixels = (int)lroundf(fmaxf((float)PARTICLE_RADIUS, gui_particle_draw_min_px));

        // Szilard mode:
        //  - species view: full fill encodes x
        //  - memory view: full fill encodes y
        //  - overlay view: fill encodes x, outer ring encodes y
        if (cli_experiment_preset == EXPERIMENT_PRESET_SZILARD_ENGINE && sz_species) {
            if (sz_view_mode == SZ_VIEW_MEMORY) {
                if (sz_memory && sz_memory[i] == 0) SDL_SetRenderDrawColor(_renderer, 80, 220, 120, 255);
                else if (sz_memory && sz_memory[i] == 1) SDL_SetRenderDrawColor(_renderer, 255, 220, 80, 255);
                else SDL_SetRenderDrawColor(_renderer, 215, 215, 215, 255);
                SDL_RenderFillCircle(_renderer, pixel_x, pixel_y, radius_pixels);
            } else if (sz_view_mode == SZ_VIEW_SPECIES) {
                if (sz_species[i] == 0) SDL_SetRenderDrawColor(_renderer, 255, 80, 80, 255);
                else                    SDL_SetRenderDrawColor(_renderer, 80, 120, 255, 255);
                SDL_RenderFillCircle(_renderer, pixel_x, pixel_y, radius_pixels);
            } else {
                if (sz_memory && (sz_memory[i] == 0 || sz_memory[i] == 1)) {
                    if (sz_memory[i] == 0) SDL_SetRenderDrawColor(_renderer, 80, 220, 120, 255);
                    else                   SDL_SetRenderDrawColor(_renderer, 255, 220, 80, 255);
                    SDL_RenderFillCircle(_renderer, pixel_x, pixel_y, radius_pixels + 2);
                }
                if (sz_species[i] == 0) SDL_SetRenderDrawColor(_renderer, 255, 80, 80, 255);
                else                    SDL_SetRenderDrawColor(_renderer, 80, 120, 255, 255);
                SDL_RenderFillCircle(_renderer, pixel_x, pixel_y, radius_pixels);
            }
            continue;
        }
        if (gui_mb_speed_view && !pl_enabled) {
            float vx = (float)Vx[i], vy = (float)Vy[i];
            float speed = sqrtf(vx * vx + vy * vy);
            set_draw_color_speed_heat(speed / vmax, 255);
            SDL_RenderFillCircle(_renderer, pixel_x, pixel_y, radius_pixels);
            continue;
        }
        // ##CHRIS: Particle-life mode uses type-based colors
        if (pl_enabled && pl_particle_types) {
            int type = pl_particle_types[i];
            if (type >= 0 && type < pl_num_types) {
                int r = pl_type_colors[type][0];
                int g = pl_type_colors[type][1];
                int b = pl_type_colors[type][2];
                SDL_SetRenderDrawColor(_renderer, r, g, b, 255);
            } else {
                SDL_SetRenderDrawColor(_renderer, 255, 255, 255, 255);  // fallback white
            }
        } else if (colorcode_mode == COLORCODE_NONE) {
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

        SDL_RenderFillCircle(_renderer, pixel_x, pixel_y, radius_pixels);
    }
}


//============================================================================
// ##CHRIS: BEGIN HEAT BATH FUNCTIONS
//============================================================================

// Sample speed from 2D Maxwell-Boltzmann distribution
static inline float sample_MB_speed_2D(float temperature) {
    // For 2D: v = sqrt(-2*kB*T/m * log(rand))
    float rand1 = (float)rand() / (float)RAND_MAX;
    if (rand1 < 1e-9f) rand1 = 1e-9f; // Avoid log(0)
    float speed = sqrtf(-2.0f * kB_effective() * temperature / PARTICLE_MASS * logf(rand1));
    return speed;
}

// Adaptive thermal wall bounce - overshoots to converge faster, then stabilizes
static inline void adaptive_thermal_wall_bounce(double *vx, double *vy,
                                                  float normal_x, float normal_y,
                                                  float gas_temp, float hb_temp) {
    float sample_temp;
    float temp_error = fabsf(gas_temp - hb_temp) / hb_temp;

    if (temp_error < stability_window_percent) {
        // Within stability window - sample from correct MB distribution
        sample_temp = hb_temp;
    } else {
        // Outside stability window - use overshoot
        if (gas_temp < hb_temp) {
            sample_temp = hb_temp * mb_overshoot_factor;  // Too cold - overshoot hot
        } else {
            sample_temp = hb_temp / mb_overshoot_factor;  // Too hot - undershoot cool
        }
    }

    current_sample_temp = sample_temp;

    // Sample speed from MB distribution
    float speed = sample_MB_speed_2D(sample_temp);

    // Random angle in hemisphere (±90° from normal)
    float rand_val = (float)rand() / (float)RAND_MAX;
    float angle = (rand_val - 0.5f) * M_PI;

    // Rotate normal by random angle
    float cos_a = cosf(angle);
    float sin_a = sinf(angle);

    *vx = (double)(speed * (normal_x * cos_a - normal_y * sin_a));
    *vy = (double)(speed * (normal_y * cos_a + normal_x * sin_a));
}

// Base thermal wall bounce - standard Maxwell-Boltzmann (no overshoot)
static inline void thermal_wall_bounce(double *vx, double *vy,
                                        float normal_x, float normal_y,
                                        float hb_temp) {
    // Sample speed from MB distribution at heat bath temperature
    float speed = sample_MB_speed_2D(hb_temp);

    // Random angle in hemisphere (±90° from normal)
    float rand_val = (float)rand() / (float)RAND_MAX;
    float angle = (rand_val - 0.5f) * M_PI;

    // Rotate normal by random angle
    float cos_a = cosf(angle);
    float sin_a = sinf(angle);

    *vx = (double)(speed * (normal_x * cos_a - normal_y * sin_a));
    *vy = (double)(speed * (normal_x * sin_a + normal_y * cos_a));
}

// Gradual damping approach - simple friction-like thermalization
static inline void gradual_damping_bounce(double *vx, double *vy,
                                           float gas_temp, float hb_temp) {
    float dV = 0.01f;
    float dV_sign = (hb_temp - gas_temp) / fabsf(hb_temp - gas_temp);

    // Add small perturbation (roughbounce effect)
    float original_mag = sqrtf((*vx) * (*vx) + (*vy) * (*vy));
    float percent_error = 0.1f;

    // Simple gaussian sampling using Box-Muller
    float u1 = (float)rand() / (float)RAND_MAX;
    float u2 = (float)rand() / (float)RAND_MAX;
    if (u1 < 1e-9f) u1 = 1e-9f;
    float gaussian = sqrtf(-2.0f * logf(u1)) * cosf(2.0f * M_PI * u2);

    *vx += (double)(gaussian * percent_error * original_mag);
    *vy += (double)(gaussian * percent_error * original_mag);

    float new_mag = sqrtf((*vx) * (*vx) + (*vy) * (*vy));
    float scale = original_mag / new_mag;
    *vx *= (double)scale;
    *vy *= (double)scale;

    // Adjust velocity toward heat bath temperature
    if (fabsf(*vx) > fabsf(*vy)) {
        *vx += (double)(dV * dV_sign);
    } else {
        *vy += (double)(dV * dV_sign);
    }
}

//============================================================================
// ##CHRIS: END HEAT BATH FUNCTIONS
//============================================================================

//============================================================================
// ##CHRIS: BEGIN ANDERSEN THERMOSTAT FUNCTION
//============================================================================

// Helper function: Sample from Gaussian(0, sigma) using Box-Muller
static inline double sample_gaussian_andersen(double sigma) {
    double u1, u2, z;
    do {
        u1 = (double)rand() / RAND_MAX;
        u2 = (double)rand() / RAND_MAX;
        z = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
    } while (isnan(z));
    return sigma * z;
}

// Andersen thermostat: stochastic velocity resampling
// Each particle has probability P = nu*dt of "colliding" with heat bath
// If collision occurs, velocity is resampled from MB distribution at heatbath_temperature
static inline void andersen_thermostat_step(double dt) {
    if (!andersen_enabled || !heatbath_enabled) return;  // Only works if heat bath is active

    double collision_prob = andersen_collision_freq * dt;  // P = nu * dt

    // For each particle, check if it "collides" with the heat bath
    for (int i = 0; i < particles_active; i++) {
        double rand_val = (double)rand() / RAND_MAX;
        if (rand_val < collision_prob) {
            // Resample velocity from Maxwell-Boltzmann distribution
            // For bulk resampling (no directional constraint), we CAN use independent Gaussians
            double sigma = sqrt(K_B * heatbath_temperature / PARTICLE_MASS);
            Vx[i] = sample_gaussian_andersen(sigma);
            Vy[i] = sample_gaussian_andersen(sigma);

            // Track this collision for diagnostics
            if (last_andersen_collision_step != NULL) {
                last_andersen_collision_step[i] = steps_elapsed;
            }
        }
    }
}

//============================================================================
// ##CHRIS: END ANDERSEN THERMOSTAT FUNCTION
//============================================================================

//============================================================================
// ##CHRIS: BEGIN PARTICLE-LIFE PHYSICS FUNCTIONS
//============================================================================

// Apply cyclic/periodic boundary conditions (wrap-around at edges) or reflective walls
static inline void pl_apply_cyclic_boundaries(double *x, double *y) {
    // X and Y are in PIXELS, not sigma units!
    float width_pixels = (float)SIM_WIDTH;
    float height_pixels = (float)SIM_HEIGHT;

    if (pl_cyclic_boundaries) {
        // CYCLIC: Wrap around at edges
        // Convert to normalized coordinates [0,1]
        float x_norm = (float)((*x - (double)XW1) / width_pixels);
        float y_norm = (float)((*y - (double)YW1) / height_pixels);

        // Wrap around
        while (x_norm < 0.0f) x_norm += 1.0f;
        while (x_norm >= 1.0f) x_norm -= 1.0f;
        while (y_norm < 0.0f) y_norm += 1.0f;
        while (y_norm >= 1.0f) y_norm -= 1.0f;

        // Convert back to pixels
        *x = (double)XW1 + x_norm * width_pixels;
        *y = (double)YW1 + y_norm * height_pixels;
    } else {
        // REFLECTIVE: Clamp to boundaries (bounce handled by velocity reflection in update loop)
        double x_min = (double)XW1 + PARTICLE_RADIUS;
        double x_max = (double)(XW1 + SIM_WIDTH) - PARTICLE_RADIUS;
        double y_min = (double)YW1 + PARTICLE_RADIUS;
        double y_max = (double)(YW1 + SIM_HEIGHT) - PARTICLE_RADIUS;

        // Clamp position to stay inside boundaries
        if (*x < x_min) *x = x_min;
        if (*x > x_max) *x = x_max;
        if (*y < y_min) *y = y_min;
        if (*y > y_max) *y = y_max;
    }
}

// Compute minimum distance considering periodic boundaries
static inline void pl_minimum_image_vector(double x1, double y1, double x2, double y2,
                                             double *dx_out, double *dy_out) {
    if (!pl_cyclic_boundaries) {
        *dx_out = x2 - x1;
        *dy_out = y2 - y1;
        return;
    }

    // X and Y are in PIXELS
    float width_pixels = (float)SIM_WIDTH;
    float height_pixels = (float)SIM_HEIGHT;

    double dx = x2 - x1;
    double dy = y2 - y1;

    // Apply minimum image convention
    if (dx > 0.5 * width_pixels) dx -= width_pixels;
    else if (dx < -0.5 * width_pixels) dx += width_pixels;

    if (dy > 0.5 * height_pixels) dy -= height_pixels;
    else if (dy < -0.5 * height_pixels) dy += height_pixels;

    *dx_out = dx;
    *dy_out = dy;
}

// Particle-life force calculation (from particle-life-app Main.java)
// Returns acceleration vector based on interaction matrix OR tensor
// ##CHRIS: Extended to support velocity and angular tensors
static inline void pl_compute_force(int type_i, int type_j, double dx, double dy,
                                      double vx_i, double vy_i, double vx_j, double vy_j,
                                      double *ax_out, double *ay_out) {
    // Particle-life formula expects normalized distance [0,1]
    // We normalize by rmax to get distances in [0,1] range
    const double beta_normalized = 0.3;  // Fixed from original particle-life

    double dist = sqrt(dx*dx + dy*dy);

    if (dist < 1e-12) {
        *ax_out = 0.0;
        *ay_out = 0.0;
        return;
    }

    // ##CHRIS: Get interaction strength from tensor (mode-dependent) or matrix
    float a;

    if (pl_tensor_mode == TENSOR_MODE_ANGULAR) {
        // 5D tensor: M[i][j][dist][v_rel][angle]
        int dist_bin = pl_get_distance_bin(dist);

        // Compute relative velocity magnitude
        double dvx = vx_j - vx_i;
        double dvy = vy_j - vy_i;
        double v_rel = sqrt(dvx*dvx + dvy*dvy);
        int vel_bin = pl_get_velocity_bin(v_rel);

        // Compute angle between velocity vector of i and position vector from i to j
        // This gives "vision cone" - particle i only strongly interacts with j if j is in i's "field of view"
        double angle;
        if (fabs(vx_i) < 1e-6 && fabs(vy_i) < 1e-6) {
            angle = 0.0;  // Stationary particle sees all directions equally
        } else {
            // Angle between velocity of i and direction to j
            double dot = vx_i * dx + vy_i * dy;
            double v_mag = sqrt(vx_i*vx_i + vy_i*vy_i);
            double cos_angle = dot / (v_mag * dist);
            cos_angle = fmax(-1.0, fmin(1.0, cos_angle));  // Clamp to [-1,1]
            angle = acos(cos_angle);  // Result in [0, pi]
        }
        int angle_bin = pl_get_angle_bin(angle);

        a = pl_tensor_5d[type_i][type_j][dist_bin][vel_bin][angle_bin];

    } else if (pl_tensor_mode == TENSOR_MODE_VELOCITY) {
        // 4D tensor: M[i][j][dist][v_rel]
        int dist_bin = pl_get_distance_bin(dist);

        // Compute relative velocity magnitude
        double dvx = vx_j - vx_i;
        double dvy = vy_j - vy_i;
        double v_rel = sqrt(dvx*dvx + dvy*dvy);
        int vel_bin = pl_get_velocity_bin(v_rel);

        a = pl_tensor_4d[type_i][type_j][dist_bin][vel_bin];

    } else if (pl_tensor_mode == TENSOR_MODE_DISTANCE) {
        // 3D tensor: M[i][j][dist]
        int bin = pl_get_distance_bin(dist);
        a = pl_tensor_3d[type_i][type_j][bin];

    } else {
        // 2D matrix: M[i][j]
        a = pl_interaction_matrix[type_i][type_j];  // Fallback to simple matrix
    }

    // Normalize distance by rmax (so dist is in [0,1] when within interaction range)
    double dist_norm = dist / pl_rmax;

    // Normalize direction
    double dir_x = dx / dist;
    double dir_y = dy / dist;

    // Compute force magnitude (from particle-life formula using normalized distance)
    double force_mag;
    if (dist_norm < beta_normalized) {
        // Close range: repulsion (even if a > 0)
        force_mag = (dist_norm / beta_normalized) - 1.0;  // Ranges from -1 (at 0) to 0 (at beta)
    } else {
        // Long range: interaction depends on matrix value a
        force_mag = a * (1.0 - fabs(1.0 + beta_normalized - 2.0 * dist_norm) / (1.0 - beta_normalized));
    }

    // Apply global force multiplier and scale to pixels (multiply by some constant for visible motion)
    force_mag *= pl_force_multiplier * 1000.0;  // Scale factor for pixel-based world (increase for faster motion)

    // Return acceleration (force * direction)
    *ax_out = dir_x * force_mag;
    *ay_out = dir_y * force_mag;
}

// Update particle-life physics for all particles
static void pl_update_physics(double dt) {
    if (!pl_enabled || !pl_particle_types) return;

    int n = particles_active;

    // Allocate acceleration arrays
    static double *ax = NULL;
    static double *ay = NULL;
    static int alloc_size = 0;

    if (n > alloc_size) {
        ax = (double*)realloc(ax, n * sizeof(double));
        ay = (double*)realloc(ay, n * sizeof(double));
        alloc_size = n;
    }

    // Zero accelerations
    for (int i = 0; i < n; i++) {
        ax[i] = 0.0;
        ay[i] = 0.0;
    }

    // Compute pairwise forces (O(N^2) - can optimize with spatial partitioning later)
    for (int i = 0; i < n; i++) {
        int type_i = pl_particle_types[i];

        for (int j = i + 1; j < n; j++) {
            int type_j = pl_particle_types[j];

            // Compute minimum image vector (handles periodic boundaries)
            double dx, dy;
            pl_minimum_image_vector(X[i], Y[i], X[j], Y[j], &dx, &dy);

            double dist = sqrt(dx*dx + dy*dy);

            // Skip if beyond interaction radius
            if (dist > pl_rmax) continue;

            // ##CHRIS: Extremely weak repulsion - barely noticeable "safety net"
            // Lets tensor forces fully dominate to allow stable emergent structures
            double hard_core_dist = 2.0 * PARTICLE_RADIUS;
            if (dist < hard_core_dist * 0.95 && dist > 1e-10) {
                // Extremely weak repulsive force - only activates when particles actually overlap
                double sigma_over_r = hard_core_dist / dist;
                double repulsion_strength = 0.1;  // Extremely weak - barely affects dynamics
                double force_magnitude = repulsion_strength * pow(sigma_over_r, 6);  // r^-6 softer repulsion

                // Very low cap - structures can stay together
                if (force_magnitude > 2.0) force_magnitude = 2.0;

                double nx = dx / dist;
                double ny = dy / dist;

                // Apply repulsion equally to both particles
                ax[i] -= force_magnitude * nx;
                ay[i] -= force_magnitude * ny;
                ax[j] += force_magnitude * nx;
                ay[j] += force_magnitude * ny;
            }

            // Compute force i->j (particle i interacting with particle j)
            double ax_ij, ay_ij;
            pl_compute_force(type_i, type_j, dx, dy, Vx[i], Vy[i], Vx[j], Vy[j], &ax_ij, &ay_ij);

            // Compute force j->i (particle j interacting with particle i)
            double ax_ji, ay_ji;
            pl_compute_force(type_j, type_i, -dx, -dy, Vx[j], Vy[j], Vx[i], Vy[i], &ax_ji, &ay_ji);

            // Accumulate accelerations
            ax[i] += ax_ij;
            ay[i] += ay_ij;
            ax[j] += ax_ji;
            ay[j] += ay_ji;

            // ##CHRIS: Check for hard-disk collision and resolve immediately
            double collision_dist = 2.0 * PARTICLE_RADIUS;  // Sum of radii (assuming equal size)
            if (dist < collision_dist && dist > 1e-10) {
                // Elastic collision response (exchange velocity components along collision axis)
                double nx = dx / dist;  // Collision normal (normalized)
                double ny = dy / dist;

                // Relative velocity
                double dvx = Vx[i] - Vx[j];
                double dvy = Vy[i] - Vy[j];
                double dvn = dvx * nx + dvy * ny;  // Relative velocity along normal

                // Resolve collision (both approaching and overlapping cases)
                if (dvn > 0) {
                    // Elastic collision: exchange normal velocity components (equal mass)
                    Vx[i] -= dvn * nx;
                    Vy[i] -= dvn * ny;
                    Vx[j] += dvn * nx;
                    Vy[j] += dvn * ny;
                }

                // Always separate overlapping particles with extra buffer
                double overlap = collision_dist - dist;
                double sep = overlap * 0.5 + 0.3;  // Larger buffer to prevent persistent overlaps
                X[i] += sep * nx;
                Y[i] += sep * ny;
                X[j] -= sep * nx;
                Y[j] -= sep * ny;
            }
        }
    }

    // Apply friction and integrate (Euler method)
    for (int i = 0; i < n; i++) {
        // Apply friction to velocity (exponential decay)
        Vx[i] *= pl_friction;
        Vy[i] *= pl_friction;

        // Update velocity with acceleration
        Vx[i] += ax[i] * dt;
        Vy[i] += ay[i] * dt;

        // ##CHRIS: Add thermal noise to maintain constant temperature
        // Stronger noise to balance friction energy loss
        double noise_strength = 1.5 * (1.0 - pl_friction);
        double angle = ((double)rand() / RAND_MAX) * 2.0 * M_PI;
        Vx[i] += noise_strength * cos(angle);
        Vy[i] += noise_strength * sin(angle);

        // Store old position for collision detection
        double old_x = X[i];
        double old_y = Y[i];

        // Update position
        X[i] += Vx[i] * dt;
        Y[i] += Vy[i] * dt;

        // ##CHRIS: Handle boundaries - either cyclic wrapping or elastic wall reflection
        if (pl_cyclic_boundaries) {
            // CYCLIC: wrap around at edges
            pl_apply_cyclic_boundaries(&X[i], &Y[i]);
        } else {
            // REFLECTIVE: elastic wall bouncing with proper collision physics
            double x_min = (double)XW1 + PARTICLE_RADIUS;
            double x_max = (double)(XW1 + SIM_WIDTH) - PARTICLE_RADIUS;
            double y_min = (double)YW1 + PARTICLE_RADIUS;
            double y_max = (double)(YW1 + SIM_HEIGHT) - PARTICLE_RADIUS;

            // Handle X boundaries (left/right walls)
            if (X[i] < x_min) {
                // Hit left wall - reflect position and velocity
                double overshoot = x_min - X[i];
                X[i] = x_min + overshoot;  // Bounce back by overshoot distance
                Vx[i] = -Vx[i];  // Reverse X velocity (elastic collision)
            } else if (X[i] > x_max) {
                // Hit right wall - reflect position and velocity
                double overshoot = X[i] - x_max;
                X[i] = x_max - overshoot;  // Bounce back by overshoot distance
                Vx[i] = -Vx[i];  // Reverse X velocity (elastic collision)
            }

            // Handle Y boundaries (top/bottom walls)
            if (Y[i] < y_min) {
                // Hit bottom wall - reflect position and velocity
                double overshoot = y_min - Y[i];
                Y[i] = y_min + overshoot;  // Bounce back by overshoot distance
                Vy[i] = -Vy[i];  // Reverse Y velocity (elastic collision)
            } else if (Y[i] > y_max) {
                // Hit top wall - reflect position and velocity
                double overshoot = Y[i] - y_max;
                Y[i] = y_max - overshoot;  // Bounce back by overshoot distance
                Vy[i] = -Vy[i];  // Reverse Y velocity (elastic collision)
            }
        }
    }
}

//============================================================================
// ##CHRIS: END PARTICLE-LIFE PHYSICS FUNCTIONS
//============================================================================

//============================================================================
// ##CHRIS: BEGIN PARTICLE-LIFE INITIALIZATION AND UTILITIES
//============================================================================

// Forward declaration
static bool pl_load_tensor_file(const char *filename);

// Initialize particle-life experiment
static void pl_initialize(void) {
    printf("🦠 Initializing PARTICLE-LIFE experiment...\n");

    pl_enabled = 1;
    int n = particles_active;

    // ##CHRIS: Load custom tensor file FIRST (sets pl_num_types before assigning particle types)
    extern char* cli_pl_tensor_file;  // Declared in CLI section
    extern int cli_pl_cyclic_boundaries;  // CLI override for boundary mode
    bool custom_tensor_loaded = false;
    if (cli_pl_tensor_file) {
        if (pl_load_tensor_file(cli_pl_tensor_file)) {
            custom_tensor_loaded = true;
            printf("✅ Custom tensor loaded from %s - pl_num_types=%d\n", cli_pl_tensor_file, pl_num_types);
        } else {
            fprintf(stderr, "❌ Failed to load tensor file, using defaults\n");
        }
    }

    // Apply CLI override for cyclic boundaries (if specified)
    if (cli_pl_cyclic_boundaries >= 0) {
        pl_cyclic_boundaries = cli_pl_cyclic_boundaries;
        printf("🔄 Boundary mode overridden by CLI: %s\n", pl_cyclic_boundaries ? "CYCLIC (wrapping)" : "REFLECTIVE (walls)");
    } else {
        printf("🔄 Boundary mode from tensor: %s\n", pl_cyclic_boundaries ? "CYCLIC (wrapping)" : "REFLECTIVE (walls)");
    }

    // Allocate type array
    if (pl_particle_types) free(pl_particle_types);
    pl_particle_types = (int*)malloc(n * sizeof(int));

    // Assign types randomly (uniform distribution) - NOW uses correct pl_num_types from tensor file
    for (int i = 0; i < n; i++) {
        pl_particle_types[i] = rand() % pl_num_types;
    }
    printf("🎨 Assigned %d particles to %d types\n", n, pl_num_types);

    // Initialize all radii to 1.0 (same size for now)
    for (int i = 0; i < MAX_PARTICLE_TYPES; i++) {
        pl_particle_radii[i] = 1.0f;
    }

    // ##CHRIS: Initialize distance bins (5 bins: close/near/medium/far/very-far)
    // Only initialize defaults if no custom tensor was loaded
    if (!custom_tensor_loaded) {
    pl_distance_bin_edges[0] = 0.0f;
    pl_distance_bin_edges[1] = 40.0f;    // Bin 0: Very close [0, 40)
    pl_distance_bin_edges[2] = 80.0f;    // Bin 1: Close [40, 80)
    pl_distance_bin_edges[3] = 120.0f;   // Bin 2: Medium [80, 120)
    pl_distance_bin_edges[4] = 160.0f;   // Bin 3: Far [120, 160)
    pl_distance_bin_edges[5] = 200.0f;   // Bin 4: Very far [160, 200)

    // ##CHRIS: Initialize velocity bins (3 bins: slow/medium/fast relative velocity)
    pl_num_velocity_bins = 3;
    pl_velocity_bin_edges[0] = 0.0f;     // Bin 0: Slow [0, 50) px/s
    pl_velocity_bin_edges[1] = 50.0f;    // Bin 1: Medium [50, 150) px/s
    pl_velocity_bin_edges[2] = 150.0f;   // Bin 2: Fast [150+) px/s
    pl_velocity_bin_edges[3] = 1000.0f;  // Upper bound

    // ##CHRIS: Initialize angle bins (8 bins: 45° sectors for vision cones)
    pl_num_angle_bins = 8;
    for (int i = 0; i <= pl_num_angle_bins; i++) {
        pl_angle_bin_edges[i] = (float)i * M_PI / (float)pl_num_angle_bins;
    }
    // Bins: [0°,22.5°), [22.5°,45°), [45°,67.5°), [67.5°,90°), [90°,112.5°), [112.5°,135°), [135°,157.5°), [157.5°,180°)

    // Initialize interaction TENSOR with PREDATOR-PREY pattern
    // Red (type 0) = Predator, Green (type 1) = Prey, Blue (type 2) = Neutral
    for (int i = 0; i < pl_num_types; i++) {
        for (int j = 0; j < pl_num_types; j++) {
            for (int b = 0; b < pl_num_distance_bins; b++) {
                // Default: weak self-attraction
                if (i == j) {
                    pl_tensor_3d[i][j][b] = (b == 2) ? 0.3f : 0.1f;
                } else {
                    pl_tensor_3d[i][j][b] = 0.0f;
                }
            }
        }
    }

    // RED (0) chases GREEN (1): Strong attraction at all distances
    pl_tensor_3d[0][1][0] = +0.8f;   // Very close: Chase!
    pl_tensor_3d[0][1][1] = +1.0f;   // Close: Strong chase
    pl_tensor_3d[0][1][2] = +0.7f;   // Medium: Moderate chase
    pl_tensor_3d[0][1][3] = +0.4f;   // Far: Weak chase
    pl_tensor_3d[0][1][4] = +0.2f;   // Very far: Very weak chase

    // GREEN (1) flees from RED (0): Repulsion
    pl_tensor_3d[1][0][0] = -1.0f;   // Very close: FLEE!
    pl_tensor_3d[1][0][1] = -0.8f;   // Close: Strong flee
    pl_tensor_3d[1][0][2] = -0.5f;   // Medium: Moderate flee
    pl_tensor_3d[1][0][3] = -0.2f;   // Far: Weak flee
    pl_tensor_3d[1][0][4] = -0.1f;   // Very far: Very weak flee

    // BLUE (2) is neutral - just weak self-clustering
    // (already set above in default loop)

    // ##CHRIS: Initialize 4D tensor (velocity-dependent forces)
    // Predator-prey with velocity dependence: Fast predator chases harder, fast prey flees harder
    for (int i = 0; i < pl_num_types; i++) {
        for (int j = 0; j < pl_num_types; j++) {
            for (int d = 0; d < pl_num_distance_bins; d++) {
                for (int v = 0; v < pl_num_velocity_bins; v++) {
                    // Copy from 3D tensor as base
                    pl_tensor_4d[i][j][d][v] = pl_tensor_3d[i][j][d];
                }
            }
        }
    }

    // RED chasing GREEN: Stronger chase when prey is slow (easier to catch)
    for (int d = 0; d < pl_num_distance_bins; d++) {
        pl_tensor_4d[0][1][d][0] = pl_tensor_3d[0][1][d] * 1.5f;  // Slow prey: 1.5x chase
        pl_tensor_4d[0][1][d][1] = pl_tensor_3d[0][1][d] * 1.0f;  // Medium prey: normal
        pl_tensor_4d[0][1][d][2] = pl_tensor_3d[0][1][d] * 0.6f;  // Fast prey: weaker chase (can't catch)
    }

    // GREEN fleeing RED: Stronger flee when predator is fast (panic!)
    for (int d = 0; d < pl_num_distance_bins; d++) {
        pl_tensor_4d[1][0][d][0] = pl_tensor_3d[1][0][d] * 0.8f;  // Slow predator: moderate flee
        pl_tensor_4d[1][0][d][1] = pl_tensor_3d[1][0][d] * 1.2f;  // Medium predator: stronger flee
        pl_tensor_4d[1][0][d][2] = pl_tensor_3d[1][0][d] * 1.5f;  // Fast predator: PANIC flee!
    }

    // ##CHRIS: Initialize 5D tensor (angular vision cones)
    // Predator only chases prey in front, prey flees in all directions
    for (int i = 0; i < pl_num_types; i++) {
        for (int j = 0; j < pl_num_types; j++) {
            for (int d = 0; d < pl_num_distance_bins; d++) {
                for (int v = 0; v < pl_num_velocity_bins; v++) {
                    for (int a = 0; a < pl_num_angle_bins; a++) {
                        // Copy from 4D tensor as base
                        pl_tensor_5d[i][j][d][v][a] = pl_tensor_4d[i][j][d][v];
                    }
                }
            }
        }
    }

    // RED chasing GREEN: Only chase prey in forward vision cone (bins 0-2: 0°-67.5°)
    for (int d = 0; d < pl_num_distance_bins; d++) {
        for (int v = 0; v < pl_num_velocity_bins; v++) {
            // Forward vision (0°-67.5°): Full chase strength
            pl_tensor_5d[0][1][d][v][0] = pl_tensor_4d[0][1][d][v] * 1.2f;  // Dead ahead: strong
            pl_tensor_5d[0][1][d][v][1] = pl_tensor_4d[0][1][d][v] * 1.1f;  // Slight angle
            pl_tensor_5d[0][1][d][v][2] = pl_tensor_4d[0][1][d][v] * 0.9f;  // Peripheral

            // Side vision (67.5°-112.5°): Weak chase
            pl_tensor_5d[0][1][d][v][3] = pl_tensor_4d[0][1][d][v] * 0.3f;
            pl_tensor_5d[0][1][d][v][4] = pl_tensor_4d[0][1][d][v] * 0.3f;

            // Rear vision (112.5°-180°): No chase (can't see behind)
            pl_tensor_5d[0][1][d][v][5] = 0.0f;
            pl_tensor_5d[0][1][d][v][6] = 0.0f;
            pl_tensor_5d[0][1][d][v][7] = 0.0f;
        }
    }

    // GREEN fleeing RED: Flee in all directions (omnidirectional fear), stronger away from predator
    for (int d = 0; d < pl_num_distance_bins; d++) {
        for (int v = 0; v < pl_num_velocity_bins; v++) {
            // When prey is moving away from predator (bins 5-7: 112.5°-180°), flee harder
            pl_tensor_5d[1][0][d][v][5] = pl_tensor_4d[1][0][d][v] * 1.3f;
            pl_tensor_5d[1][0][d][v][6] = pl_tensor_4d[1][0][d][v] * 1.4f;
            pl_tensor_5d[1][0][d][v][7] = pl_tensor_4d[1][0][d][v] * 1.5f;  // Running directly away: maximum flee
        }
    }

    // Also set fallback matrix (for when tensors disabled)
    for (int i = 0; i < pl_num_types; i++) {
        for (int j = 0; j < pl_num_types; j++) {
            pl_interaction_matrix[i][j] = (i == j) ? 1.0f : 0.0f;
        }
    }
    } // End of default initialization (!custom_tensor_loaded)

    // Print interaction tensor info
    const char* tensor_mode_names[] = {"2D Matrix", "3D Distance", "4D Velocity", "5D Angular"};
    if (custom_tensor_loaded) {
        printf("📊 Tensor Mode: %s - CUSTOM LOADED FROM FILE\n", tensor_mode_names[pl_tensor_mode]);
    } else {
        printf("📊 Tensor Mode: %s - PREDATOR-PREY DYNAMICS\n", tensor_mode_names[pl_tensor_mode]);
        printf("   🔴 Red = Predator (chases green)\n");
        printf("   🟢 Green = Prey (flees from red)\n");
        printf("   🔵 Blue = Neutral (weak clustering)\n");
    }

    if (pl_tensor_mode == TENSOR_MODE_DISTANCE) {
        printf("   Distance bins: [0,40) [40,80) [80,120) [120,160) [160,200)\n");
        for (int i = 0; i < pl_num_types; i++) {
            for (int j = 0; j < pl_num_types; j++) {
                printf("   Type %d→%d: ", i, j);
                for (int b = 0; b < pl_num_distance_bins; b++) {
                    printf("%6.2f ", pl_tensor_3d[i][j][b]);
                }
                printf("\n");
            }
        }
    } else if (pl_tensor_mode == TENSOR_MODE_VELOCITY) {
        printf("   Distance bins: [0,40) [40,80) [80,120) [120,160) [160,200)\n");
        printf("   Velocity bins: [0,50) [50,150) [150+) px/s\n");
        printf("   Red→Green chase strength varies with prey velocity (slow=1.5x, fast=0.6x)\n");
        printf("   Green→Red flee strength varies with predator velocity (slow=0.8x, fast=1.5x)\n");
    } else if (pl_tensor_mode == TENSOR_MODE_ANGULAR) {
        printf("   Distance bins: [0,40) [40,80) [80,120) [120,160) [160,200)\n");
        printf("   Velocity bins: [0,50) [50,150) [150+) px/s\n");
        printf("   Angle bins: 8 sectors of 22.5° (vision cones)\n");
        printf("   Red has forward vision cone (can't chase prey behind)\n");
        printf("   Green flees omnidirectionally (stronger when running away)\n");
    } else {
        printf("   Simple 2D Matrix (%dx%d):\n", pl_num_types, pl_num_types);
        for (int i = 0; i < pl_num_types; i++) {
            printf("   ");
            for (int j = 0; j < pl_num_types; j++) {
                printf("%6.2f ", pl_interaction_matrix[i][j]);
            }
            printf("\n");
        }
    }

    printf("🎨 Particle types: %d\n", pl_num_types);
    printf("🔄 Cyclic boundaries: %s\n", pl_cyclic_boundaries ? "ON" : "OFF");
    printf("📏 Interaction radius (rmax): %.4f\n", pl_rmax);
    printf("💨 Friction coefficient: %.4f\n", pl_friction);
    printf("⚡ Force multiplier: %.4f\n", pl_force_multiplier);
    printf("✅ Particle-life initialized with %d particles\n", n);
}

// ##CHRIS: Load custom tensor from file
// File format:
//   # Comments start with #
//   MODE <matrix|distance|velocity|angular>
//   TYPES <num_types>
//   DISTANCE_BINS <num_bins> <edge0> <edge1> ...
//   VELOCITY_BINS <num_bins> <edge0> <edge1> ...
//   ANGLE_BINS <num_bins> <edge0> <edge1> ...
//   TENSOR <type_i> <type_j> [<dist_bin>] [<vel_bin>] [<angle_bin>] <value>
static bool pl_load_tensor_file(const char *filename) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "❌ Error: Cannot open tensor file '%s'\n", filename);
        return false;
    }

    printf("📂 Loading tensor from file: %s\n", filename);
    char line[512];
    int line_num = 0;

    while (fgets(line, sizeof(line), fp)) {
        line_num++;

        // Strip newline
        size_t len = strlen(line);
        if (len > 0 && line[len-1] == '\n') line[len-1] = '\0';

        // Skip empty lines and comments
        if (line[0] == '\0' || line[0] == '#') continue;

        // Parse command
        char cmd[64];
        if (sscanf(line, "%63s", cmd) != 1) continue;

        if (strcmp(cmd, "MODE") == 0) {
            char mode_str[32];
            if (sscanf(line, "MODE %31s", mode_str) != 1) {
                fprintf(stderr, "❌ Line %d: Invalid MODE syntax\n", line_num);
                fclose(fp);
                return false;
            }
            if (strcmp(mode_str, "matrix") == 0) pl_tensor_mode = TENSOR_MODE_MATRIX;
            else if (strcmp(mode_str, "distance") == 0) pl_tensor_mode = TENSOR_MODE_DISTANCE;
            else if (strcmp(mode_str, "velocity") == 0) pl_tensor_mode = TENSOR_MODE_VELOCITY;
            else if (strcmp(mode_str, "angular") == 0) pl_tensor_mode = TENSOR_MODE_ANGULAR;
            else {
                fprintf(stderr, "❌ Line %d: Unknown mode '%s'\n", line_num, mode_str);
                fclose(fp);
                return false;
            }
            printf("   Mode: %s\n", mode_str);

        } else if (strcmp(cmd, "TYPES") == 0) {
            if (sscanf(line, "TYPES %d", &pl_num_types) != 1 || pl_num_types < 1 || pl_num_types > MAX_PARTICLE_TYPES) {
                fprintf(stderr, "❌ Line %d: Invalid TYPES (must be 1-%d)\n", line_num, MAX_PARTICLE_TYPES);
                fclose(fp);
                return false;
            }
            printf("   Types: %d\n", pl_num_types);

        } else if (strcmp(cmd, "DISTANCE_BINS") == 0) {
            int num_bins;
            char *ptr = line + strlen("DISTANCE_BINS");
            if (sscanf(ptr, "%d", &num_bins) != 1 || num_bins < 1 || num_bins > MAX_DISTANCE_BINS) {
                fprintf(stderr, "❌ Line %d: Invalid DISTANCE_BINS count\n", line_num);
                fclose(fp);
                return false;
            }
            pl_num_distance_bins = num_bins;

            // Parse bin edges
            ptr = strchr(ptr, ' ');
            for (int i = 0; i <= num_bins; i++) {
                float edge;
                if (!ptr || sscanf(ptr, "%f", &edge) != 1) {
                    fprintf(stderr, "❌ Line %d: Missing distance bin edge %d\n", line_num, i);
                    fclose(fp);
                    return false;
                }
                pl_distance_bin_edges[i] = edge;
                ptr = strchr(ptr + 1, ' ');
            }
            printf("   Distance bins: %d\n", num_bins);

        } else if (strcmp(cmd, "TENSOR") == 0) {
            // Parse tensor entry based on current mode
            int type_i, type_j;
            char *ptr = line + strlen("TENSOR");

            // Parse type indices and find where parsing stopped
            int chars_read;
            if (sscanf(ptr, "%d %d%n", &type_i, &type_j, &chars_read) != 2) {
                fprintf(stderr, "❌ Line %d: Invalid TENSOR syntax\n", line_num);
                fclose(fp);
                return false;
            }
            ptr += chars_read;  // Now ptr points right after the two indices

            if (type_i < 0 || type_i >= pl_num_types || type_j < 0 || type_j >= pl_num_types) {
                fprintf(stderr, "❌ Line %d: Type indices out of range\n", line_num);
                fclose(fp);
                return false;
            }

            if (pl_tensor_mode == TENSOR_MODE_MATRIX) {
                float value;
                if (sscanf(ptr, "%f", &value) != 1) {
                    fprintf(stderr, "❌ Line %d: Missing value for 2D tensor\n", line_num);
                    fclose(fp);
                    return false;
                }
                pl_interaction_matrix[type_i][type_j] = value;

            } else if (pl_tensor_mode == TENSOR_MODE_DISTANCE) {
                int dist_bin;
                float value;
                if (sscanf(ptr, "%d %f", &dist_bin, &value) != 2) {
                    fprintf(stderr, "❌ Line %d: Missing dist_bin or value for 3D tensor\n", line_num);
                    fclose(fp);
                    return false;
                }
                if (dist_bin < 0 || dist_bin >= pl_num_distance_bins) {
                    fprintf(stderr, "❌ Line %d: Distance bin index out of range\n", line_num);
                    fclose(fp);
                    return false;
                }
                pl_tensor_3d[type_i][type_j][dist_bin] = value;

            } else if (pl_tensor_mode == TENSOR_MODE_VELOCITY) {
                int dist_bin, vel_bin;
                float value;
                if (sscanf(ptr, "%d %d %f", &dist_bin, &vel_bin, &value) != 3) {
                    fprintf(stderr, "❌ Line %d: Missing bins or value for 4D tensor\n", line_num);
                    fclose(fp);
                    return false;
                }
                pl_tensor_4d[type_i][type_j][dist_bin][vel_bin] = value;

            } else if (pl_tensor_mode == TENSOR_MODE_ANGULAR) {
                int dist_bin, vel_bin, angle_bin;
                float value;
                if (sscanf(ptr, "%d %d %d %f", &dist_bin, &vel_bin, &angle_bin, &value) != 4) {
                    fprintf(stderr, "❌ Line %d: Missing bins or value for 5D tensor\n", line_num);
                    fclose(fp);
                    return false;
                }
                pl_tensor_5d[type_i][type_j][dist_bin][vel_bin][angle_bin] = value;
            }
        }
    }

    fclose(fp);
    printf("✅ Tensor file loaded successfully\n");
    return true;
}

// Clean up particle-life data
static void pl_cleanup(void) {
    if (pl_particle_types) {
        free(pl_particle_types);
        pl_particle_types = NULL;
    }
    pl_enabled = 0;
}

//============================================================================
// ##CHRIS: END PARTICLE-LIFE INITIALIZATION AND UTILITIES
//============================================================================

// Function to check and handle boundary collisions
// newtons secon law
// ##CHRIS: MODIFIED - Added heat bath support for all 4 outer walls
void handle_boundary_collision(int i) {
    float gas_temp = compute_measured_temperature_from_ke();  // ##CHRIS: Get current gas temperature
    const float Ri = (Radius[i] > 0.0) ? (float)Radius[i] : (float)PARTICLE_RADIUS;
    const bool sz_unstick = (cli_experiment_preset == EXPERIMENT_PRESET_SZILARD_ENGINE &&
                             sz_interactive_enabled && cli_szilard_unstick);
    bool hit_left = false;
    bool hit_right = false;

    // Left wall (x = XW1, normal points right)
    if (X[i] - Ri < XW1) {
        hit_left = true;
        X[i] = XW1 + Ri;
        if (heatbath_enabled && fabsf(gas_temp - heatbath_temperature) > 0.0001f * heatbath_temperature) {
            // ##CHRIS: Heat bath active - select thermal mode
            if (thermal_wall_mode == 1) {
                thermal_wall_bounce(&Vx[i], &Vy[i], 1.0f, 0.0f, heatbath_temperature);
            } else if (thermal_wall_mode == 2) {
                adaptive_thermal_wall_bounce(&Vx[i], &Vy[i], 1.0f, 0.0f, gas_temp, heatbath_temperature);
            } else {
                gradual_damping_bounce(&Vx[i], &Vy[i], gas_temp, heatbath_temperature);
            }
        } else {
            // Normal elastic reflection
            Vx[i] = -Vx[i];
        }
    }

    // Right wall (x = XW2, normal points left)
    if (X[i] + Ri > XW2) {
        hit_right = true;
        X[i] = XW2 - Ri;
        if (heatbath_enabled && fabsf(gas_temp - heatbath_temperature) > 0.0001f * heatbath_temperature) {
            // ##CHRIS: Heat bath active - select thermal mode
            if (thermal_wall_mode == 1) {
                thermal_wall_bounce(&Vx[i], &Vy[i], -1.0f, 0.0f, heatbath_temperature);
            } else if (thermal_wall_mode == 2) {
                adaptive_thermal_wall_bounce(&Vx[i], &Vy[i], -1.0f, 0.0f, gas_temp, heatbath_temperature);
            } else {
                gradual_damping_bounce(&Vx[i], &Vy[i], gas_temp, heatbath_temperature);
            }
        } else {
            // Normal elastic reflection
            Vx[i] = -Vx[i];
        }
    }

    // Bottom wall (y = YW1, normal points up)
    if (Y[i] - Ri < YW1) {
        Y[i] = YW1 + Ri;
        if (heatbath_enabled && fabsf(gas_temp - heatbath_temperature) > 0.0001f * heatbath_temperature) {
            // ##CHRIS: Heat bath active - select thermal mode
            if (thermal_wall_mode == 1) {
                thermal_wall_bounce(&Vx[i], &Vy[i], 0.0f, 1.0f, heatbath_temperature);
            } else if (thermal_wall_mode == 2) {
                adaptive_thermal_wall_bounce(&Vx[i], &Vy[i], 0.0f, 1.0f, gas_temp, heatbath_temperature);
            } else {
                gradual_damping_bounce(&Vx[i], &Vy[i], gas_temp, heatbath_temperature);
            }
        } else {
            // Normal elastic reflection
            Vy[i] = -Vy[i];
        }
    }

    // Top wall (y = YW2, normal points down)
    if (Y[i] + Ri > YW2) {
        Y[i] = YW2 - Ri;
        if (heatbath_enabled && fabsf(gas_temp - heatbath_temperature) > 0.0001f * heatbath_temperature) {
            // ##CHRIS: Heat bath active - select thermal mode
            if (thermal_wall_mode == 1) {
                thermal_wall_bounce(&Vx[i], &Vy[i], 0.0f, -1.0f, heatbath_temperature);
            } else if (thermal_wall_mode == 2) {
                adaptive_thermal_wall_bounce(&Vx[i], &Vy[i], 0.0f, -1.0f, gas_temp, heatbath_temperature);
            } else {
                gradual_damping_bounce(&Vx[i], &Vy[i], gas_temp, heatbath_temperature);
            }
        } else {
            // Normal elastic reflection
            Vy[i] = -Vy[i];
        }
    }

    // Szilard-only: if a particle keeps getting pinned to an outer wall with ~zero Vx,
    // redirect a small fraction of its speed into x while conserving |v| (KE).
    if (sz_unstick) {
        int dir = hit_left ? +1 : (hit_right ? -1 : 0); // push away from wall
        if (dir != 0) {
            if (fabs(Vx[i]) < 1e-6) sz_unstick_hit_count[i]++; else sz_unstick_hit_count[i] = 0;
            if (sz_unstick_hit_count[i] >= cli_szilard_unstick_hits) {
                const double s2 = Vx[i] * Vx[i] + Vy[i] * Vy[i];
                if (s2 > 1e-18) {
                    const double s = sqrt(s2);
                    double vx_mag = fmax(1e-4, (double)cli_szilard_unstick_vx_frac * s);
                    vx_mag = fmin(vx_mag, 0.95 * s);
                    const double vx_new = (double)dir * vx_mag;
                    const double vy_sign = (Vy[i] >= 0.0) ? 1.0 : -1.0;
                    const double vy_new = vy_sign * sqrt(fmax(0.0, s * s - vx_new * vx_new));
                    Vx[i] = vx_new;
                    Vy[i] = vy_new;
                }
                // Ensure we are not re-penetrating the wall after the nudge.
                const float eps = 1e-6f;
                if (hit_left) X[i] = XW1 + Ri + eps;
                if (hit_right) X[i] = XW2 - Ri - eps;
                sz_unstick_hit_count[i] = 0;
            }
        } else {
            sz_unstick_hit_count[i] = 0;
        }
    }
}
// ##CHRIS: END MODIFIED

//============================================================================
// ##CHRIS: Andersen thermostat call (after wall collisions)
// This provides bulk thermalization to prevent "cold core" problem
// where slow particles in center rarely hit walls
//============================================================================
static inline void apply_andersen_if_enabled(double dt) {
    if (andersen_enabled && heatbath_enabled) {
        andersen_thermostat_step(dt);
    }
}
//============================================================================



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

// CCD-style piston collision using start/end positions within the substep.
// Uses X_old[i] (start-of-substep) and current X[i] (after wall CCD) to detect crossings.
static inline void handle_piston_collisions_ccd_substep(int i, float sub_dt,
                                                        float pL0, float pR0,
                                                        float vL, float vR) {
    const float eps = 1e-6f;
    const float R = Radius[i];
    const float x0 = X_old[i];
    const float x1 = X[i];
    const double m = (double)PARTICLE_MASS;

    // Left piston (particle hits from the right): contact when (x - R) == pL
    {
        const float pL1 = pL0 + vL * sub_dt;
        const float r0 = (x0 - R) - pL0;
        const float r1 = (x1 - R) - pL1;
        if (r0 > 0.0f && r1 <= 0.0f) {
            float denom = r0 - r1;
            float f = (fabsf(denom) > 1e-12f) ? (r0 / denom) : 0.0f;
            if (f < 0.0f) f = 0.0f; if (f > 1.0f) f = 1.0f;
            float t = f * sub_dt;
            float pC = pL0 + vL * t;
            double u = (double)Vx[i];
            double v = 2.0 * (double)vL - u;
            double dE = 0.5 * m * (v * v - u * u);
            piston_work_left += dE;
            piston_hits_left++;
            Vx[i] = (float)v;
            X[i] = pC + R + eps + Vx[i] * (sub_dt - t);
        } else if (r0 <= 0.0f) {
            // Already overlapping at start: push out and reflect if moving into piston
            float vrel = Vx[i] - vL;
            if (vrel < 0.0f) {
                double u = (double)Vx[i];
                double v = 2.0 * (double)vL - u;
                double dE = 0.5 * m * (v * v - u * u);
                piston_work_left += dE;
                piston_hits_left++;
                Vx[i] = (float)v;
            }
            X[i] = pL0 + R + eps;
        }
    }

    // Right piston (particle hits from the left): contact when (x + R) == pR
    {
        const float pR1 = pR0 + vR * sub_dt;
        const float r0 = (x0 + R) - pR0;
        const float r1 = (x1 + R) - pR1;
        if (r0 < 0.0f && r1 >= 0.0f) {
            float denom = r0 - r1;
            float f = (fabsf(denom) > 1e-12f) ? (r0 / denom) : 0.0f;
            if (f < 0.0f) f = 0.0f; if (f > 1.0f) f = 1.0f;
            float t = f * sub_dt;
            float pC = pR0 + vR * t;
            double u = (double)Vx[i];
            double v = 2.0 * (double)vR - u;
            double dE = 0.5 * m * (v * v - u * u);
            piston_work_right += dE;
            piston_hits_right++;
            Vx[i] = (float)v;
            X[i] = pC - R - eps + Vx[i] * (sub_dt - t);
        } else if (r0 >= 0.0f) {
            // Already overlapping at start: push out and reflect if moving into piston
            float vrel = Vx[i] - vR;
            if (vrel > 0.0f) {
                double u = (double)Vx[i];
                double v = 2.0 * (double)vR - u;
                double dE = 0.5 * m * (v * v - u * u);
                piston_work_right += dE;
                piston_hits_right++;
                Vx[i] = (float)v;
            }
            X[i] = pR0 - R - eps;
        }
    }
}



//////// WALL FUNCTIONS /////////

void draw_wall() {
    if (!wall_enabled) return;  // Don't draw if the wall is disabled
       if (isnan(wall_x) || wall_x < XW1 || wall_x > XW2) {
        printf("🚨 Invalid wall_x = %f\n", wall_x);
            exit(1);
        }
    const bool sz_mode = (cli_experiment_preset == EXPERIMENT_PRESET_SZILARD_ENGINE && sz_interactive_enabled);
    float vis_thick = (wall_thickness_visual_runtime > 0.0f) ? wall_thickness_visual_runtime : WALL_THICKNESS;
    // In Szilard mode use a thin visual wall for both red/blue barriers.
    float draw_thick = sz_mode ? fmaxf(2.0f, vis_thick * 0.25f) : vis_thick;
    SDL_Rect wall_rect = { (int)(wall_x - (draw_thick/2)), 0, (int)draw_thick, SIM_HEIGHT };  // Thin vertical wall

    if (sz_mode) {
        // Szilard visual language:
        // - x-gated walls use species colors (red/blue)
        // - y-gated walls use memory colors (green/yellow)
        if (sz_perm_use_memory) {
            if (sz_perm_species == 0) SDL_SetRenderDrawColor(_renderer, 80, 220, 120, 255);   // y=0 / green
            else                      SDL_SetRenderDrawColor(_renderer, 255, 220, 80, 255);   // y=1 / yellow
        } else {
            if (sz_interactive_phase == 3) SDL_SetRenderDrawColor(_renderer, 255, 80, 80, 255);
            else SDL_SetRenderDrawColor(_renderer, 90, 140, 255, 255);
        }
    } else {
        SDL_SetRenderDrawColor(_renderer, 255, 255, 255, 255);  // default white
    }

    SDL_RenderFillRect(_renderer, &wall_rect);

    // Fixed semipermeable marker at L0 while active.
    if (sz_mode && sz_fixed_gate_active) {
        float gate_wf = fmaxf(2.0f, draw_thick);
        int gate_w = (int)gate_wf;
        float marker_x = sz_interactive_target_px;
        if (sz_fixed_gate_use_memory) {
            if (sz_fixed_gate_species == 0) SDL_SetRenderDrawColor(_renderer, 80, 220, 120, 255);   // y=0 / green
            else                            SDL_SetRenderDrawColor(_renderer, 255, 220, 80, 255);   // y=1 / yellow
        } else {
            if (sz_fixed_gate_species == 0) SDL_SetRenderDrawColor(_renderer, 255, 80, 80, 255);
            else SDL_SetRenderDrawColor(_renderer, 90, 140, 255, 255);
        }
        if (fabsf(wall_x - marker_x) < draw_thick * 0.75f) {
            // Overlap state: draw colored side rails around moving wall so both stay visible.
            int half = (int)fmaxf(1.0f, draw_thick * 0.5f);
            SDL_Rect left_rail  = { (int)(marker_x - half - 1), 0, 1, SIM_HEIGHT };
            SDL_Rect right_rail = { (int)(marker_x + half), 0, 1, SIM_HEIGHT };
            SDL_RenderFillRect(_renderer, &left_rail);
            SDL_RenderFillRect(_renderer, &right_rail);
        } else {
            SDL_Rect gate_rect = { (int)(marker_x - gate_w / 2), 0, gate_w, SIM_HEIGHT };
            SDL_RenderFillRect(_renderer, &gate_rect);
        }
    }

    for (int i = 0; i < extra_wall_count; ++i) {
        float wx = extra_wall_positions[i];
        float vth = (wall_thickness_visual_runtime > 0.0f) ? wall_thickness_visual_runtime : WALL_THICKNESS;
        SDL_Rect extra_rect = { (int)(wx - (vth/2)), 0, (int)vth, SIM_HEIGHT };
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
            // Use σ-time units for logging/analysis consistency.
            // `simulation_time` is advanced in the main loop as `fixed_dt_runtime / PIXELS_PER_SIGMA`.
            wall_release_time = simulation_time;
            primary_wall_release_time = simulation_time;
            if (!cli_quiet) {
                printf("🔔 Wall released at step %d (t = %.3f)\n", steps_elapsed, wall_release_time);
            }
        }
        return;
    }

    // === Normal Wall Motion Phase ===
    wall_x_old = wall_x;  // Store previous wall position

    if (!time_mode_walls_integrated_this_step) {
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
    }

    if (!cli_quiet && (steps_elapsed % 100 == 0)) {
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
        time_mode_walls_integrated_this_step = false;
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

        if (!time_mode_walls_integrated_this_step) {
            extra_wall_old_positions[j] = extra_wall_positions[j];
            extra_wall_positions[j] += extra_wall_velocity[j] * dt;
        }

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
    time_mode_walls_integrated_this_step = false;
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
In ccd_wall_step(i, dt, wall_x0, wall_vx_live, prev_side, wall_impulse_x_accum, wall_mass, force_locked):
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
                          const float wall_mass, // wall mass for finite-mass response (absolute units)
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
    const bool locked_wall = force_locked || (wall_hold_enabled && !wall_is_released) || (wall_mass <= 0.0f); // held before release or explicitly locked

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
                const float M = wall_mass;
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
                        *wall_vx_live -= Jp / M;
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
            const float M = wall_mass;            // wall mass
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
                    *wall_vx_live -= Jp / M;
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
            const float M = wall_mass;
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
                    *wall_vx_live -= Jp / M;
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
void update_particles_with_substepping(float dt, int* out_left, int* out_right,
                                       double current_sim_time, bool is_wall_released, double wall_released_at) {
    if (dt <= 0.0f) return;
    time_mode_walls_integrated_this_step = false;

    int n_active = particles_active > 0 ? particles_active : 0;

    int substeps = 1;
    #if SUBSTEPPING
        substeps = SUBSTEPS;
    #endif
    if (cli_substeps_override > 0) {
        substeps = cli_substeps_override;
    } else if (cli_substeps_auto) {
        // Adaptive substepping heuristic (fast + robust):
        // Require relative displacement per substep to be small compared to particle diameter.
        //   2*vmax*dt/substeps <= dx_frac*(2R)
        // so:
        //   substeps >= (2*vmax*dt) / (dx_frac*(2R))
        double vmax2 = 0.0;
        for (int i = 0; i < n_active; ++i) {
            const double vx = Vx[i];
            const double vy = Vy[i];
            const double v2 = vx * vx + vy * vy;
            if (v2 > vmax2) vmax2 = v2;
        }
        // Include moving walls/pistons (x-direction only)
        {
            const double vw = (double)vx_wall;
            const double v2 = vw * vw;
            if (v2 > vmax2) vmax2 = v2;
        }
        if (extra_wall_count > 0 && extra_wall_velocity) {
            for (int w = 0; w < extra_wall_count; ++w) {
                const double vw = (double)extra_wall_velocity[w];
                const double v2 = vw * vw;
                if (v2 > vmax2) vmax2 = v2;
            }
        }
        {
            const double vL = (double)vx_piston_left;
            const double v2 = vL * vL;
            if (v2 > vmax2) vmax2 = v2;
        }
        {
            const double vR = (double)vx_piston_right;
            const double v2 = vR * vR;
            if (v2 > vmax2) vmax2 = v2;
        }

        const double vmax = (vmax2 > 0.0) ? sqrt(vmax2) : 0.0;
        const double diameter_px = 2.0 * (double)PARTICLE_RADIUS;
        double dx_allow_px = (double)cli_substeps_auto_dx_frac * diameter_px;
        if (!isfinite(dx_allow_px) || dx_allow_px <= 1e-12) dx_allow_px = 1e-12;

        const double dx_rel_px = 2.0 * vmax * (double)dt;
        int desired = (dx_rel_px > 0.0) ? (int)ceil(dx_rel_px / dx_allow_px) : cli_substeps_auto_min;
        if (desired < cli_substeps_auto_min) desired = cli_substeps_auto_min;
        if (desired > cli_substeps_auto_max) desired = cli_substeps_auto_max;
        if (desired < 1) desired = 1;
        substeps = desired;
    }
    // Expose for HUD/diagnostics (e.g. energy_transfer_wall_positions.csv).
    g_last_substeps_used = substeps;
    float sub_dt = dt / (float)substeps;

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

    // Substep-local wall positions (to account for moving boundaries within the fixed step)
    float wall_x_sub = wall_x_old;
    float *extra_wall_sub = NULL;
    if (extra_wall_count > 0) {
        extra_wall_sub = (float *)malloc((size_t)extra_wall_count * sizeof(float));
        if (!extra_wall_sub) {
            printf("❌ Failed to allocate extra_wall_sub buffer.\n");
            exit(1);
        }
        for (int j = 0; j < extra_wall_count; ++j) {
            extra_wall_sub[j] = extra_wall_old_positions[j];
        }
    }

    for (int step = 0; step < substeps; ++step) {
        // ##CHRIS: Track substep time for substep-level CSV logging
        double substep_time = g_substep_base_time + step * (sub_dt / PIXELS_PER_SIGMA);

        const float pL0_sub = piston_left_x + vx_piston_left * (step * sub_dt);
        const float pR0_sub = piston_right_x + vx_piston_right * (step * sub_dt);
        const float pL1_sub = pL0_sub + vx_piston_left * sub_dt;
        const float pR1_sub = pR0_sub + vx_piston_right * sub_dt;
        if (wall_enabled &&
            energy_measurement.spring_enabled &&
            num_internal_walls > 0 &&
            !(wall_hold_enabled && !wall_is_released) &&
            wall_mass_runtime > 0.0f) {
            const float disp = wall_x_sub - energy_measurement.equilibrium_position;
            const float accel = (-energy_measurement.spring_constant * disp) / wall_mass_runtime;
            vx_wall += accel * sub_dt;
        }
        const float w_x0_sub = wall_x_sub;   // substep pose (moving boundary)

        for (int i = 0; i < n_active; i++) {
            X_old[i] = X[i];

            if (wall_enabled) {
                if (X_old[i] < w_x0_sub - WALL_THICKNESS / 2)      Side_old[i] = -1;
                else if (X_old[i] > w_x0_sub + WALL_THICKNESS / 2) Side_old[i] =  1;
                else                                             Side_old[i] =  0;
            } else {
                Side_old[i] = 0;
            }

            const int prev_side = Side_old[i];   // -1 left, 0 inside, +1 right

            // CCD handles ballistic advance + wall interaction (now updating wall velocity instantly).
            if (cli_experiment_preset == EXPERIMENT_PRESET_SZILARD_ENGINE &&
                sz_species &&
                sz_perm_wall_index == primary_wall_index) {
                bool apply_gate = false;
                int gate_target_side = 0; // 0=left, 1=right
                const int gate_label = szilard_particle_gate_label(i, sz_perm_use_memory);
                if (sz_perm_dual_active) {
                    // Dual semipermeable gates at L0:
                    // species 0 always filtered left, species 1 always filtered right.
                    if (gate_label == 0 || gate_label == 1) {
                        apply_gate = true;
                        gate_target_side = (gate_label == 0) ? 0 : 1;
                    }
                } else if (sz_perm_active) {
                    // Single active gate (used while divider is moving).
                    gate_target_side = sz_perm_target_side;
                    if (gate_label == sz_perm_species || sz_partial_obs_enabled) {
                        apply_gate = true;
                    }
                }
                if (apply_gate) {
                    // Semipermeable one-way gate:
                    //  - target_side=left: allow particles on right/inside to pass left, block left->right
                    //  - target_side=right: allow particles on left/inside to pass right, block right->left
                    const float R = Radius[i];
                    const float wL0 = w_x0_sub - 0.5f * WALL_THICKNESS;
                    const float wR0 = w_x0_sub + 0.5f * WALL_THICKNESS;
                    const float L0p = wL0 - R;
                    const float R0p = wR0 + R;

                    bool allow_pass = false;
                    if (gate_target_side == 0) {
                        // keep species on left: allow crossing from right->left (or inside)
                        if (gate_label == sz_perm_species && prev_side >= 0) allow_pass = true;
                    } else {
                        // keep species on right: allow crossing from left->right (or inside)
                        if (gate_label == sz_perm_species && prev_side <= 0) allow_pass = true;
                    }

                    if (!allow_pass && sz_partial_obs_enabled && !sz_perm_dual_active &&
                        szilard_particle_hot_enough_for_gate_leak(i)) {
                        allow_pass = true;
                    }

                    if (allow_pass) {
                        X[i] += Vx[i] * sub_dt;
                        Y[i] += Vy[i] * sub_dt;
                        // If we land inside the wall thickness, nudge through to avoid "sticking"
                        if (X[i] > L0p && X[i] < R0p) {
                            const float eps = fmaxf(1e-7f, 1e-6f * fmaxf(R, 1.f));
                            if (gate_target_side == 0) X[i] = L0p - eps;
                            else X[i] = R0p + eps;
                        }
                    } else {
                        bool sz_lock = (cli_experiment_preset == EXPERIMENT_PRESET_SZILARD_ENGINE && sz_interactive_enabled);
                        ccd_wall_step(i, sub_dt, w_x0_sub, &vx_wall, prev_side, &wall_impulse_x_accum, wall_mass_runtime, sz_lock);
                    }
                } else {
                    bool sz_lock = (cli_experiment_preset == EXPERIMENT_PRESET_SZILARD_ENGINE && sz_interactive_enabled);
                    ccd_wall_step(i, sub_dt, w_x0_sub, &vx_wall, prev_side, &wall_impulse_x_accum, wall_mass_runtime, sz_lock);
                }
            } else {
                bool sz_lock = (cli_experiment_preset == EXPERIMENT_PRESET_SZILARD_ENGINE && sz_interactive_enabled);
                ccd_wall_step(i, sub_dt, w_x0_sub, &vx_wall, prev_side, &wall_impulse_x_accum, wall_mass_runtime, sz_lock);
            }

            if (wall_enabled && extra_wall_count > 0) {
                for (int w = 0; w < extra_wall_count; ++w) {
                    float wx = extra_wall_sub ? extra_wall_sub[w] : extra_wall_old_positions[w];
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
                    float m_wall = (extra_wall_masses && w < extra_wall_count) ? extra_wall_masses[w] : wall_mass_runtime;
                    ccd_wall_step(i, sub_dt, wx, &extra_wall_velocity[w], extra_prev_side, NULL, m_wall, hold_active);
                }
            }

            handle_piston_collisions_ccd_substep(i, sub_dt, pL0_sub, pR0_sub, vx_piston_left, vx_piston_right);
            handle_boundary_collision(i);

            // Szilard fixed species gate at L0.
            if (cli_experiment_preset == EXPERIMENT_PRESET_SZILARD_ENGINE &&
                sz_interactive_enabled &&
                sz_fixed_gate_active &&
                szilard_particle_gate_label(i, sz_fixed_gate_use_memory) == sz_fixed_gate_species) {
                float gate_vis_thick = (wall_thickness_visual_runtime > 0.0f) ? wall_thickness_visual_runtime : WALL_THICKNESS;
                float gate_thick = fmaxf(2.0f, gate_vis_thick * 0.25f);
                float gate_x = sz_interactive_target_px;
                float left_face_center = gate_x - 0.5f * gate_thick - Radius[i];
                float right_face_center = gate_x + 0.5f * gate_thick + Radius[i];
                const float eps = fmaxf(1e-7f, 1e-6f * fmaxf(Radius[i], 1.f));

                if (sz_fixed_gate_block_dir == 0) {
                    // Block left->right crossings
                    if (X_old[i] <= left_face_center && X[i] > left_face_center) {
                        if (!(sz_partial_obs_enabled && szilard_particle_hot_enough_for_gate_leak(i))) {
                            X[i] = left_face_center - eps;
                            if (Vx[i] > 0.0f) Vx[i] = -Vx[i];
                        }
                    }
                } else {
                    // Block right->left crossings
                    if (X_old[i] >= right_face_center && X[i] < right_face_center) {
                        if (!(sz_partial_obs_enabled && szilard_particle_hot_enough_for_gate_leak(i))) {
                            X[i] = right_face_center + eps;
                            if (Vx[i] < 0.0f) Vx[i] = -Vx[i];
                        }
                    }
                }
            }
        }

        // Advance substep-local wall positions using (possibly updated) velocities.
        wall_x_sub += vx_wall * sub_dt;
        if (extra_wall_sub && extra_wall_velocity) {
            for (int j = 0; j < extra_wall_count; ++j) {
                extra_wall_sub[j] += extra_wall_velocity[j] * sub_dt;
            }
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
        if (!cli_no_pp_collisions) {
            if (!cli_pp_toi) {
                for (int i = 0; i < n_active; ++i) {
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
            } else {
                const int ring = cli_pp_toi_neighbor_ring < 0 ? 0 : cli_pp_toi_neighbor_ring;
                for (int iter = 0; iter < cli_pp_toi_max_events; ++iter) {
                    int hits = 0;
                    for (int i = 0; i < n_active; ++i) {
                        int cell_index = particle_cell[i];
                        int cx = cell_index % grid_cols;
                        int cy = cell_index / grid_cols;

                        for (int dy = -ring; dy <= ring; ++dy) {
                            int ny = cy + dy;
                            if (ny < 0 || ny >= grid_rows) continue;
                            for (int dx = -ring; dx <= ring; ++dx) {
                                int nx = cx + dx;
                                if (nx < 0 || nx >= grid_cols) continue;
                                int neighbor_idx = ny * grid_cols + nx;
                                for (int j = cell_head[neighbor_idx]; j != -1; j = cell_next[j]) {
                                    if (j <= i) continue;

                                    // Relative state at end of substep
                                    double rx1 = (double)(X[i] - X[j]);
                                    double ry1 = (double)(Y[i] - Y[j]);
                                    double r2_1 = rx1 * rx1 + ry1 * ry1;
                                    double min_d = 2.0 * (double)PARTICLE_RADIUS;
                                    double min_d2 = min_d * min_d;
                                    if (r2_1 >= min_d2 || r2_1 <= 0.0) continue;

                                    double dist = sqrt(r2_1);
                                    double overlap = min_d - dist;
                                    if (overlap > missed_threshold) {
                                        missed_collision_events++;
                                        if (overlap > worst_penetration_observed) worst_penetration_observed = overlap;
                                    }

                                    double vx = (double)(Vx[i] - Vx[j]);
                                    double vy = (double)(Vy[i] - Vy[j]);
                                    const double dt_step = (double)sub_dt;

                                    // Relative position at start of substep
                                    double rx0 = rx1 - vx * dt_step;
                                    double ry0 = ry1 - vy * dt_step;
                                    double r2_0 = rx0 * rx0 + ry0 * ry0;

                                    // Solve |r0 + v*t| = d for t in [0, dt]
                                    double t_hit = 0.0;
                                    bool   has_hit = false;
                                    if (r2_0 <= min_d2) {
                                        // Already overlapping at substep start
                                        t_hit = 0.0;
                                        has_hit = true;
                                    } else {
                                        double a = vx * vx + vy * vy;
                                        if (a > 1e-12) {
                                            double b = 2.0 * (rx0 * vx + ry0 * vy);
                                            double c = r2_0 - min_d2;
                                            double disc = b * b - 4.0 * a * c;
                                            if (disc >= 0.0) {
                                                double sqrtD = sqrt(disc);
                                                double t0 = (-b - sqrtD) / (2.0 * a);
                                                if (t0 >= 0.0 && t0 <= dt_step) {
                                                    t_hit = t0;
                                                    has_hit = true;
                                                }
                                            }
                                        }
                                    }
                                    if (!has_hit) continue;

                                    // Contact normal at impact
                                    double rx_c = rx0 + vx * t_hit;
                                    double ry_c = ry0 + vy * t_hit;
                                    double r2_c = rx_c * rx_c + ry_c * ry_c;
                                    if (r2_c <= 1e-18) continue;
                                    double inv_r = 1.0 / sqrt(r2_c);
                                    double nxp = rx_c * inv_r;
                                    double nyp = ry_c * inv_r;

                                    double vxi0 = (double)Vx[i];
                                    double vyi0 = (double)Vy[i];
                                    double vxj0 = (double)Vx[j];
                                    double vyj0 = (double)Vy[j];
                                    double vn = vx * nxp + vy * nyp;

                                    // Elastic collision along normal only if approaching
                                    double vxi1 = vxi0, vyi1 = vyi0;
                                    double vxj1 = vxj0, vyj1 = vyj0;
                                    if (vn < 0.0) {
                                        vxi1 = vxi0 - vn * nxp;
                                        vyi1 = vyi0 - vn * nyp;
                                        vxj1 = vxj0 + vn * nxp;
                                        vyj1 = vyj0 + vn * nyp;
                                    }

                                    // Rewind to impact, then advance with new velocities for remaining time
                                    double t_remain = dt_step - t_hit;
                                    double xi_hit = (double)X[i] - vxi0 * t_remain;
                                    double yi_hit = (double)Y[i] - vyi0 * t_remain;
                                    double xj_hit = (double)X[j] - vxj0 * t_remain;
                                    double yj_hit = (double)Y[j] - vyj0 * t_remain;
                                    X[i] = (float)(xi_hit + vxi1 * t_remain);
                                    Y[i] = (float)(yi_hit + vyi1 * t_remain);
                                    X[j] = (float)(xj_hit + vxj1 * t_remain);
                                    Y[j] = (float)(yj_hit + vyj1 * t_remain);

                                    Vx[i] = (float)vxi1;
                                    Vy[i] = (float)vyi1;
                                    Vx[j] = (float)vxj1;
                                    Vy[j] = (float)vyj1;

                                    // Residual overlap correction (tiny)
                                    double dist_after = sqrt((double)((X[i]-X[j])*(X[i]-X[j]) + (Y[i]-Y[j])*(Y[i]-Y[j])));
                                    double overlap2 = min_d - dist_after;
                                    if (overlap2 > 0.0) {
                                        double push = 0.5 * overlap2;
                                        X[i] += (float)(nxp * push);
                                        Y[i] += (float)(nyp * push);
                                        X[j] -= (float)(nxp * push);
                                        Y[j] -= (float)(nyp * push);
                                    }

                                    hits++;
                                }
                            }
                        }
                    }
                    if (hits == 0) break;
                }
            }
        }

        // Szilard: PP overlap resolution can push particles across barriers after the CCD/gate handling.
        // Re-enforce divider (moving wall) and fixed L0 gate constraints once per substep to prevent
        // apparent "teleporting" through semipermeable/solid walls.
        if (cli_experiment_preset == EXPERIMENT_PRESET_SZILARD_ENGINE && sz_interactive_enabled && sz_species && wall_enabled &&
            sz_perm_wall_index == primary_wall_index) {
            const float w_x1_sub = wall_x_sub; // wall pose at end of this substep (already advanced)
            const float wL0_0 = w_x0_sub - 0.5f * WALL_THICKNESS;
            const float wR0_0 = w_x0_sub + 0.5f * WALL_THICKNESS;
            const float wL0_1 = w_x1_sub - 0.5f * WALL_THICKNESS;
            const float wR0_1 = w_x1_sub + 0.5f * WALL_THICKNESS;

            for (int i = 0; i < n_active; ++i) {
                const float R = Radius[i];
                const float L0p = wL0_1 - R;
                const float R0p = wR0_1 + R;
                const float eps = fmaxf(1e-7f, 1e-6f * fmaxf(R, 1.f));

                int prev_side_local = 0;
                if (X_old[i] < wL0_0) prev_side_local = -1;
                else if (X_old[i] > wR0_0) prev_side_local = 1;

                bool apply_gate = false;
                int gate_target_side = 0;
                const int gate_label = szilard_particle_gate_label(i, sz_perm_use_memory);
                if (sz_perm_dual_active) {
                    if (gate_label == 0 || gate_label == 1) {
                        apply_gate = true;
                        gate_target_side = (gate_label == 0) ? 0 : 1;
                    }
                } else if (sz_perm_active && gate_label == sz_perm_species) {
                    apply_gate = true;
                    gate_target_side = sz_perm_target_side;
                }

                // If particle is inside the wall thickness (numerical), snap it out to avoid sticking.
                if (X[i] > L0p && X[i] < R0p) {
                    if (apply_gate) {
                        X[i] = (gate_target_side == 0) ? (L0p - eps) : (R0p + eps);
                    } else {
                        int face = prev_side_local;
                        if (face == 0) {
                            face = (fabsf(X_old[i] - L0p) <= fabsf(R0p - X_old[i])) ? -1 : +1;
                        }
                        X[i] = (face == -1) ? (L0p - eps) : (R0p + eps);
                    }
                }

                if (apply_gate) {
                    // One-way: block movement *away* from target side.
                    if (gate_target_side == 0) {
                        // Block left->right crossings for particles that were on the left.
                        if (prev_side_local == -1 && X[i] > L0p) {
                            X[i] = L0p - eps;
                            if (Vx[i] > 0.0f) Vx[i] = -Vx[i];
                        }
                    } else {
                        // Block right->left crossings for particles that were on the right.
                        if (prev_side_local == 1 && X[i] < R0p) {
                            X[i] = R0p + eps;
                            if (Vx[i] < 0.0f) Vx[i] = -Vx[i];
                        }
                    }
                } else {
                    // Solid wall: prevent crossings in either direction.
                    if (prev_side_local == -1 && X[i] > L0p) {
                        X[i] = L0p - eps;
                        if (Vx[i] > 0.0f) Vx[i] = -Vx[i];
                    } else if (prev_side_local == 1 && X[i] < R0p) {
                        X[i] = R0p + eps;
                        if (Vx[i] < 0.0f) Vx[i] = -Vx[i];
                    }
                }

                // Also re-enforce the fixed L0 gate if active.
                if (sz_fixed_gate_active &&
                    szilard_particle_gate_label(i, sz_fixed_gate_use_memory) == sz_fixed_gate_species) {
                    float gate_vis_thick = (wall_thickness_visual_runtime > 0.0f) ? wall_thickness_visual_runtime : WALL_THICKNESS;
                    float gate_thick = fmaxf(2.0f, gate_vis_thick * 0.25f);
                    float gate_x = sz_interactive_target_px;
                    float left_face_center = gate_x - 0.5f * gate_thick - R;
                    float right_face_center = gate_x + 0.5f * gate_thick + R;

                    if (sz_fixed_gate_block_dir == 0) {
                        if (X_old[i] <= left_face_center && X[i] > left_face_center) {
                            X[i] = left_face_center - eps;
                            if (Vx[i] > 0.0f) Vx[i] = -Vx[i];
                        }
                    } else {
                        if (X_old[i] >= right_face_center && X[i] < right_face_center) {
                            X[i] = right_face_center + eps;
                            if (Vx[i] < 0.0f) Vx[i] = -Vx[i];
                        }
                    }
                }

                // If the left pocket is collapsing during H/J/O, any particle that is
                // supposed to pass to the right should not get numerically trapped there.
                if ((sz_interactive_phase == 3 || sz_interactive_phase == 4 || sz_interactive_phase == 5) &&
                    vx_wall < 0.0f) {
                    int moving_target_side = -1;
                    const float left_face_center = wL0_1 - R;
                    const float pocket_w = left_face_center - ((float)XW1 + R + eps);
                    if (pocket_w <= (2.5f * R + 2.0f * eps) &&
                        szilard_particle_moving_gate_target_side(i, &moving_target_side) &&
                        moving_target_side == 1 &&
                        X[i] <= left_face_center + eps)
                    {
                        X[i] = R0p + eps;
                    }
                }

                // If the right piston closes onto the fixed L0 gate, any particle that is
                // allowed to move right->left through that gate should be released rather
                // than getting crushed in the tiny numerical pocket.
                if (sz_fixed_gate_active && vx_piston_right < 0.0f) {
                    float gate_vis_thick = (wall_thickness_visual_runtime > 0.0f) ? wall_thickness_visual_runtime : WALL_THICKNESS;
                    float gate_thick = fmaxf(2.0f, gate_vis_thick * 0.25f);
                    float gate_x = sz_interactive_target_px;
                    float gate_right_face = gate_x + 0.5f * gate_thick + R;
                    float gate_left_face = gate_x - 0.5f * gate_thick - R - eps;
                    float pocket_w = (pR1_sub - R - eps) - gate_right_face;
                    if (pocket_w <= (2.5f * R + 2.0f * eps) &&
                        X[i] >= gate_right_face - eps &&
                        szilard_particle_can_cross_fixed_gate_right_to_left(i))
                    {
                        X[i] = gate_left_face;
                    }
                }

                // Final piston-face enforcement after PP overlap correction and gate repair.
                if (X[i] - R <= pL1_sub) {
                    float vrel = Vx[i] - vx_piston_left;
                    if (vrel < 0.0f) Vx[i] = 2.0f * vx_piston_left - Vx[i];
                    X[i] = pL1_sub + R + eps;
                }
                if (X[i] + R >= pR1_sub) {
                    float vrel = Vx[i] - vx_piston_right;
                    if (vrel > 0.0f) Vx[i] = 2.0f * vx_piston_right - Vx[i];
                    X[i] = pR1_sub - R - eps;
                }
            }
        }

        // ##CHRIS: Substep-level CSV logging when output_dt <= 1.0
        // TIME/RK4 interpretation:
        //   - cli_output_dt <= 1.0 → "max resolution" logging at every substep
        //   - cli_output_dt  > 1.0 → coarse logging handled in the main loop
        if (g_substep_log && is_wall_released && cli_output_dt <= 1.0f) {
            // Count particles left/right of wall
            int left = 0, right = 0;
            for (int i = 0; i < n_active; i++) {
                if (X[i] < wall_x - WALL_THICKNESS / 2) left++;
                else if (X[i] > wall_x + WALL_THICKNESS / 2) right++;
            }

            double wall_x_sigma = wall_x / PIXELS_PER_SIGMA;
            double disp = wall_x_sigma - (XW1 + XW2) / (2.0 * PIXELS_PER_SIGMA);
            double time_since_release = substep_time - wall_released_at;

            fprintf(g_substep_log, "%.6f, %.6f, %.6f, %d, %d\n",
                    time_since_release, wall_x_sigma, disp, left, right);
            fflush(g_substep_log);
        }
    }

    if (wall_enabled) {
        wall_x = wall_x_sub;
        if (all_wall_positions && primary_wall_index >= 0 && primary_wall_index < num_internal_walls) {
            all_wall_positions[primary_wall_index] = wall_x;
        }
        time_mode_walls_integrated_this_step = true;
    }
    if (extra_wall_sub) {
        for (int j = 0; j < extra_wall_count; ++j) {
            extra_wall_positions[j] = extra_wall_sub[j];
        }
        time_mode_walls_integrated_this_step = true;
        free(extra_wall_sub);
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
                int locked_wall = (wall_hold_enabled && !wall_is_released) || (wall_mass_runtime <= 0.0f);
                double v1, v2;
                if (locked_wall) {
                    /* infinite-mass mirror in wall frame */
                    v1 = 2.0 * u2 - u1;
                    v2 = u2;
                } else {
                    double m = (double)PARTICLE_MASS;
                    double M = (double)wall_mass_runtime;
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
                    /* finite/infinite mass along x with wall's instantaneous vx (if provided) */
                    double u1 = Vx[i];
                    double u2 = (extra_wall_velocity ? extra_wall_velocity[w] : 0.0f);
                    double v1, v2;
                    int locked = 0;
                    if (extra_wall_hold_enabled && extra_wall_hold_enabled[w]) {
                        locked = 1;
                        if (extra_wall_is_released) locked = locked && !extra_wall_is_released[w];
                    }
                    double M = (extra_wall_masses && w < extra_wall_count) ? (double)extra_wall_masses[w]
                                                                          : (double)wall_mass_runtime;
                    if (locked || M <= 0.0) {
                        v1 = 2.0 * u2 - u1;
                        v2 = u2;
                    } else {
                        double m = (double)PARTICLE_MASS;
                        v1 = ((m - M) * u1 + 2.0 * M * u2) / (m + M);
                        v2 = ((M - m) * u2 + 2.0 * m * u1) / (m + M);
                    }
                    Vx[i] = (float)v1;
                    if (!locked && extra_wall_velocity) extra_wall_velocity[w] = (float)v2;
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

static inline double kinetic_energy_walls_total(void) {
    double ke = 0.0;
    // Primary wall (spring anchor)
    ke += 0.5 * (double)wall_mass_runtime * (double)vx_wall * (double)vx_wall;
    // Extra walls (if any)
    if (extra_wall_count > 0 && extra_wall_velocity) {
        for (int w = 0; w < extra_wall_count; ++w) {
            double M = (extra_wall_masses && w < extra_wall_count) ? (double)extra_wall_masses[w]
                                                                  : (double)wall_mass_runtime;
            double v = (double)extra_wall_velocity[w];
            ke += 0.5 * M * v * v;
        }
    }
    return ke;
}

double kinetic_energy_total_system() {
    double ke_particles = kinetic_energy(); // existing function (particles only)
    return ke_particles + kinetic_energy_walls_total();
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
    double ke_wall = kinetic_energy_walls_total();
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

// --- Szilard helpers ---
static void szilard_assign_species(int n_active) {
    if (n_active <= 0) return;
    if (!sz_species) sz_species = (int*)malloc((size_t)NUM_PARTICLES * sizeof(int));
    if (!sz_memory)  sz_memory  = (int*)malloc((size_t)NUM_PARTICLES * sizeof(int));
    if (!sz_species || !sz_memory) {
        printf("❌ Failed to allocate Szilard species/memory arrays\n");
        exit(1);
    }
    int nA = (int)llround((double)n_active * (double)cli_szilard_species_ratio);
    if (nA < 0) nA = 0; if (nA > n_active) nA = n_active;
    for (int i = 0; i < n_active; ++i) {
        sz_species[i] = (i < nA) ? 0 : 1;
        sz_memory[i] = -1;
    }
    // simple shuffle to randomize species distribution
    for (int i = n_active - 1; i > 0; --i) {
        int j = rand() % (i + 1);
        int tmp = sz_species[i]; sz_species[i] = sz_species[j]; sz_species[j] = tmp;
    }
}

static void szilard_count_species(double split_x,
                                  int *leftA, int *leftB, int *rightA, int *rightB,
                                  int *memA, int *memB) {
    if (leftA) *leftA = 0;
    if (leftB) *leftB = 0;
    if (rightA) *rightA = 0;
    if (rightB) *rightB = 0;
    if (memA) *memA = 0;
    if (memB) *memB = 0;
    int n = particles_active > 0 ? particles_active : 0;
    for (int i = 0; i < n; ++i) {
        int sp = sz_species ? sz_species[i] : 0;
        if (sz_memory && sz_memory[i] == 0 && memA) (*memA)++;
        if (sz_memory && sz_memory[i] == 1 && memB) (*memB)++;
        if (X[i] < split_x) {
            if (sp == 0 && leftA) (*leftA)++;
            else if (sp == 1 && leftB) (*leftB)++;
        } else {
            if (sp == 0 && rightA) (*rightA)++;
            else if (sp == 1 && rightB) (*rightB)++;
        }
    }
}

static void szilard_species_side_counts_at_split(double split_x,
                                                 int *nL0, int *nL1, int *nR0, int *nR1) {
    szilard_count_species(split_x, nL0, nL1, nR0, nR1, NULL, NULL);
}

static int szilard_count_side_for_perfect_sep(int i) {
    const float gate_x = sz_interactive_target_px;
    const float margin = 0.5f * WALL_THICKNESS + PARTICLE_RADIUS + 1e-3f;
    const float x = (float)X[i];

    if (x <= gate_x - margin) return 0;
    if (x >= gate_x + margin) return 1;

    /* Inside the membrane-adjacent ambiguity band: during G/hold, count particles
       by their intended target side so slow quasi-static runs are not marked as
       "wrong" just because a center is grazing the barrier. */
    if (sz_interactive_enabled &&
        cli_experiment_preset == EXPERIMENT_PRESET_SZILARD_ENGINE &&
        sz_interactive_phase == 1) {
        const int sp = (sz_species && sz_species[i] != 0) ? 1 : 0;
        return sp; /* species0->left(0), species1->right(1) */
    }

    return (x < gate_x) ? 0 : 1;
}

static void szilard_memory_counts(int *mem0, int *mem1, int *mem_unknown) {
    if (mem0) *mem0 = 0;
    if (mem1) *mem1 = 0;
    if (mem_unknown) *mem_unknown = 0;
    int n = particles_active > 0 ? particles_active : 0;
    for (int i = 0; i < n; ++i) {
        int m = (sz_memory ? sz_memory[i] : -1);
        if (m == 0) {
            if (mem0) (*mem0)++;
        } else if (m == 1) {
            if (mem1) (*mem1)++;
        } else {
            if (mem_unknown) (*mem_unknown)++;
        }
    }
}

static void szilard_xy_agreement_counts(int *xy_match, int *xy_mismatch, int *xy_unknown) {
    if (xy_match) *xy_match = 0;
    if (xy_mismatch) *xy_mismatch = 0;
    if (xy_unknown) *xy_unknown = 0;
    if (!sz_species) return;
    int n = particles_active > 0 ? particles_active : 0;
    for (int i = 0; i < n; ++i) {
        const int sp = (sz_species[i] != 0) ? 1 : 0;
        const int mem = (sz_memory ? sz_memory[i] : -1);
        if (mem != 0 && mem != 1) {
            if (xy_unknown) (*xy_unknown)++;
        } else if (sp == mem) {
            if (xy_match) (*xy_match)++;
        } else {
            if (xy_mismatch) (*xy_mismatch)++;
        }
    }
}

static double szilard_xy_mismatch_fraction(void) {
    int xy_match = 0, xy_mismatch = 0, xy_unknown = 0;
    szilard_xy_agreement_counts(&xy_match, &xy_mismatch, &xy_unknown);
    const int known = xy_match + xy_mismatch;
    if (known <= 0) return 0.0;
    return (double)xy_mismatch / (double)known;
}

static int szilard_memory_ready(void) {
    int mem0 = 0, mem1 = 0, memU = 0;
    szilard_memory_counts(&mem0, &mem1, &memU);
    return (particles_active > 0 && memU == 0);
}

static void szilard_measure_memory_from_species(void) {
    if (!sz_memory || !sz_species) return;
    for (int i = 0; i < particles_active; ++i) {
        sz_memory[i] = sz_species[i];
    }
    if (!cli_quiet) {
        printf("🧠 Szilard: measurement stored y=x for all particles (m)\n");
    }
}

static int szilard_particle_gate_label(int i, int use_memory) {
    const int sp = sz_species ? ((sz_species[i] != 0) ? 1 : 0) : -1;
    if (!use_memory) return sp;
    if (!sz_memory) return sp;
    const int mem = sz_memory[i];
    if (mem != 0 && mem != 1) return sp;

    /* Partially observable branch:
       once x has drifted away from stored y, the gate ceases to act on a
       perfectly addressable memory label and instead sees the current x-state. */
    if (sz_partial_obs_enabled && sp >= 0 && sp != mem) return sp;
    return mem;
}

static int szilard_particle_moving_gate_target_side(int i, int *target_side_out) {
    if (!target_side_out) return 0;
    *target_side_out = -1;
    if (i < 0 || i >= particles_active) return 0;

    const int gate_label = szilard_particle_gate_label(i, sz_perm_use_memory);
    if (sz_perm_dual_active) {
        if (gate_label == 0) {
            *target_side_out = 0;
            return 1;
        }
        if (gate_label == 1) {
            *target_side_out = 1;
            return 1;
        }
        return 0;
    }

    if (sz_perm_active && gate_label == sz_perm_species) {
        *target_side_out = (sz_perm_target_side != 0) ? 1 : 0;
        return 1;
    }
    return 0;
}

static int szilard_particle_can_cross_fixed_gate_right_to_left(int i) {
    if (!sz_fixed_gate_active) return 1;
    if (i < 0 || i >= particles_active) return 0;

    const int gate_label = szilard_particle_gate_label(i, sz_fixed_gate_use_memory);
    if (gate_label != sz_fixed_gate_species) return 1; /* transparent for non-blocked class */
    if (sz_fixed_gate_block_dir == 0) return 1;        /* blocks only L->R */
    if (sz_partial_obs_enabled && szilard_particle_hot_enough_for_gate_leak(i)) return 1;
    return 0;
}

static void szilard_refresh_gate_labels(void) {
    const int n = particles_active > 0 ? particles_active : 0;
    if (n <= 0) return;
    if (!sz_gate_labels) sz_gate_labels = (int*)malloc((size_t)NUM_PARTICLES * sizeof(int));
    if (!sz_gate_labels) {
        printf("❌ Failed to allocate Szilard gate-label array\n");
        exit(1);
    }
    for (int i = 0; i < n; ++i) {
        sz_gate_labels[i] = szilard_particle_gate_label(i, 1);
    }
}

static void szilard_speed_stats(double *mean_speed, double *std_speed) {
    if (mean_speed) *mean_speed = 0.0;
    if (std_speed) *std_speed = 0.0;
    const int n = particles_active > 0 ? particles_active : 0;
    if (n <= 0) return;

    double mean = 0.0;
    double m2 = 0.0;
    int count = 0;
    for (int i = 0; i < n; ++i) {
        const double vx = (double)Vx[i];
        const double vy = (double)Vy[i];
        const double speed = sqrt(vx * vx + vy * vy);
        count++;
        const double delta = speed - mean;
        mean += delta / (double)count;
        m2 += delta * (speed - mean);
    }
    if (mean_speed) *mean_speed = mean;
    if (std_speed) *std_speed = (count >= 2) ? sqrt(m2 / (double)(count - 1)) : 0.0;
}

static int szilard_particle_hot_enough_for_gate_leak(int i) {
    if (!sz_partial_obs_enabled) return 0;
    if (i < 0 || i >= particles_active) return 0;
    double mean_speed = 0.0, std_speed = 0.0;
    szilard_speed_stats(&mean_speed, &std_speed);
    const double ref = (mean_speed > 1e-12) ? mean_speed : fmax(std_speed, 1e-12);
    const double vx = (double)Vx[i];
    const double vy = (double)Vy[i];
    const double speed = sqrt(vx * vx + vy * vy);
    return (speed >= sz_gate_hot_speed_mult * ref);
}

static double szilard_memory_gate_information_bits_at_split(float split_x) {
    int memGateL0=0,memGateL1=0,memGateLU=0,memGateR0=0,memGateR1=0,memGateRU=0;
    szilard_memory_side_counts_at_split(split_x,
                                        &memGateL0,&memGateL1,&memGateLU,
                                        &memGateR0,&memGateR1,&memGateRU);
    return szilard_mutual_information_3state_nats(
        memGateL0,memGateL1,memGateLU,memGateR0,memGateR1,memGateRU
    ) / log(2.0);
}

static double szilard_partial_uncertainty_fraction_now(void) {
    const double xy_mismatch = szilard_xy_mismatch_fraction();
    double gate_uncertainty = 0.0;
    if (sz_partial_obs_ref_imem_gate_bits > 1e-12) {
        const double i_mem_gate_bits = szilard_memory_gate_information_bits_at_split(wall_x);
        gate_uncertainty = 1.0 - (i_mem_gate_bits / sz_partial_obs_ref_imem_gate_bits);
        if (gate_uncertainty < 0.0) gate_uncertainty = 0.0;
        if (gate_uncertainty > 1.0) gate_uncertainty = 1.0;
    }
    return (xy_mismatch > gate_uncertainty) ? xy_mismatch : gate_uncertainty;
}

static void szilard_update_partial_uncertainty_peak(void) {
    if (!sz_partial_obs_enabled) return;
    const double u = szilard_partial_uncertainty_fraction_now();
    if (u > sz_partial_obs_peak_uncertainty_frac) {
        sz_partial_obs_peak_uncertainty_frac = u;
    }
}

static void szilard_apply_species_flips(double dt_sigma) {
    if (!sz_xflip_enabled || !sz_species || particles_active <= 0) return;
    if (!(dt_sigma > 0.0) || !(sz_xflip_rate > 0.0)) return;
    double p_flip = sz_xflip_rate * dt_sigma;
    if (p_flip > 0.5) p_flip = 0.5;
    for (int i = 0; i < particles_active; ++i) {
        double u = (double)rand() / (double)RAND_MAX;
        if (u < p_flip) {
            sz_species[i] = 1 - sz_species[i];
        }
    }
}

static void szilard_apply_memory_erase_min_heat(void) {
    if (!sz_memory) return;
    int known = 0;
    for (int i = 0; i < particles_active; ++i) {
        if (sz_memory[i] == 0 || sz_memory[i] == 1) known++;
        sz_memory[i] = -1;
    }
    if (known <= 0) return;
    const double T = heatbath_enabled ? (double)heatbath_temperature : (double)temperature_runtime;
    const double qmin = (double)known * kB_effective() * T * log(2.0);
    sz_memory_erase_heat_min += qmin;
    if (!cli_quiet) {
        printf("♻️  Szilard: erased %d memory bits, Q_erase,min += %.6g\n", known, qmin);
    }
}

static void szilard_manual_erase_memory(void) {
    if (!sz_interactive_enabled || cli_experiment_preset != EXPERIMENT_PRESET_SZILARD_ENGINE) return;

    if (sz_interactive_phase == 7) {
        if (!cli_quiet) {
            printf("Szilard: memory already erased in this branch.\n");
        }
        return;
    }

    if (sz_interactive_phase != 6) {
        if (!cli_quiet) {
            printf("Szilard: erase is only available after H/J/O finishes. Wait for recombination to stop first.\n");
        }
        return;
    }

    if (!szilard_memory_ready()) {
        if (!cli_quiet) {
            printf("Szilard: no defined memory remains to erase.\n");
        }
        sz_interactive_phase = 7;
        sz_log_reset_requested = 1;
        return;
    }

    szilard_apply_memory_erase_min_heat();
    sz_interactive_phase = 7;
    sz_log_reset_requested = 1;
    szilard_sync_edmd_from_globals();
    if (!cli_quiet) {
        printf("♻️  Szilard: manual memory erase complete (E). The erase jump is now explicit in the log.\n");
    }
}

static void szilard_memory_side_counts(int *mL0, int *mL1, int *mLU,
                                       int *mR0, int *mR1, int *mRU) {
    const float gate_x = sz_interactive_target_px; /* L0 in pixels */
    if (mL0) *mL0 = 0;
    if (mL1) *mL1 = 0;
    if (mLU) *mLU = 0;
    if (mR0) *mR0 = 0;
    if (mR1) *mR1 = 0;
    if (mRU) *mRU = 0;

    int n = particles_active > 0 ? particles_active : 0;
    for (int i = 0; i < n; ++i) {
        const int side = (X[i] < gate_x) ? 0 : 1; /* 0=left, 1=right */
        const int m = (sz_memory ? sz_memory[i] : -1);
        if (side == 0) {
            if (m == 0) { if (mL0) (*mL0)++; }
            else if (m == 1) { if (mL1) (*mL1)++; }
            else { if (mLU) (*mLU)++; }
        } else {
            if (m == 0) { if (mR0) (*mR0)++; }
            else if (m == 1) { if (mR1) (*mR1)++; }
            else { if (mRU) (*mRU)++; }
        }
    }
}

static void szilard_memory_side_counts_at_split(float split_x,
                                                int *mL0, int *mL1, int *mLU,
                                                int *mR0, int *mR1, int *mRU) {
    if (mL0) *mL0 = 0;
    if (mL1) *mL1 = 0;
    if (mLU) *mLU = 0;
    if (mR0) *mR0 = 0;
    if (mR1) *mR1 = 0;
    if (mRU) *mRU = 0;
    int n = particles_active > 0 ? particles_active : 0;
    for (int i = 0; i < n; ++i) {
        const int side = (X[i] < split_x) ? 0 : 1; /* 0=left, 1=right */
        const int m = (sz_memory ? sz_memory[i] : -1);
        if (side == 0) {
            if (m == 0) { if (mL0) (*mL0)++; }
            else if (m == 1) { if (mL1) (*mL1)++; }
            else { if (mLU) (*mLU)++; }
        } else {
            if (m == 0) { if (mR0) (*mR0)++; }
            else if (m == 1) { if (mR1) (*mR1)++; }
            else { if (mRU) (*mRU)++; }
        }
    }
}

// --- Szilard interactive helpers (SDL mode) ---
static void szilard_prepare_interactive_config(void) {
    // One moving divider (B') that sweeps from x=0 to x=L0
    cli_requested_walls = 1;
    cli_num_walls = 1;
    preset_custom_positions = true;
    preset_wall_count = 1;
    preset_wall_fraction[0] = 0.5f; // seed using mid-divider (L0)

    // Use random seeding so particles don't start lined up on the left wall
    cli_seeding_mode = SEEDING_RANDOM;

    // Override particle count if provided
    if (cli_szilard_n > 0) cli_override_particles = cli_szilard_n;
    int N = (cli_override_particles > 0) ? cli_override_particles : particles_active;
    if (N < 2) N = 2;
    cli_override_particles = N;

    // Seed all particles in left half (N,0)
    if (cli_particles_box_counts) { free(cli_particles_box_counts); cli_particles_box_counts = NULL; }
    cli_particles_box_counts = (int*)malloc(2 * sizeof(int));
    cli_particles_box_counts_count = 2;
    cli_particles_box_counts[0] = N;
    cli_particles_box_counts[1] = 0;
    preset_custom_absolute_counts = true;
    cli_left_empty = false;
}

static void szilard_refresh_interactive_speed(void) {
    const double sep_duration = (cli_szilard_sep_duration > 1e-9) ? cli_szilard_sep_duration : 1.0;
    const double sep_dist_px = (double)(sz_interactive_target_px - sz_interactive_home_px);
    sz_interactive_vx_px = (float)(sep_dist_px / (sep_duration * PIXELS_PER_SIGMA));
    if (sz_interactive_vx_px < 0.0f) sz_interactive_vx_px = 0.0f;
}

static void szilard_apply_interactive_velocity_current_phase(void) {
    if (!sz_interactive_enabled) return;

    if (sz_interactive_phase == 1) {
        vx_wall = sz_interactive_vx_px;
        vx_piston_right = sz_interactive_vx_px;
    } else if (sz_interactive_phase == 3 || sz_interactive_phase == 4 || sz_interactive_phase == 5) {
        vx_wall = -sz_interactive_vx_px;
        vx_piston_right = -sz_interactive_vx_px;
    } else {
        vx_wall = 0.0f;
        vx_piston_right = 0.0f;
    }

    if (cli_experiment_preset == EXPERIMENT_PRESET_SZILARD_ENGINE && sim_mode == MODE_EDMD && g_edmd) {
        szilard_sync_edmd_from_globals();
    }
}

static void szilard_setup_interactive_state(void) {
    sz_interactive_enabled = 1;
    sz_interactive_phase = 0;
    g_sz_branch_snapshot.valid = 0;
    sz_log_index = 0;
    sz_log_branch_id = 0;
    snprintf(sz_log_branch_label, sizeof(sz_log_branch_label), "Eq");
    sz_sep_perfect_since = -1.0;
    sz_xflip_enabled = 0;
    sz_xeq_ready_since = -1.0;
    sz_memory_erase_heat_min = 0.0;
    sz_partial_obs_enabled = 0;
    sz_partial_obs_ref_imem_gate_bits = 1.0;
    sz_partial_obs_peak_uncertainty_frac = 0.0;
    sz_branch_ready_reason = 0;
    sz_pending_recombine_mode = 0;
    szilard_reset_visible_load_state(0.0);
    if (sim_mode == MODE_EDMD_HYBRID) {
        printf("⚠️ Szilard interactive in EDMD-HYBRID is not supported (semipermeable walls live in EDMD core).\n");
        printf("   Use: --mode=edmd for Szilard.\n");
        sim_mode = MODE_TIME;
    }

    // Divider B' starts at the left boundary and moves to L0 during separation.
    // Right piston B starts at L0 and moves to 2*L0.
    // Note: update_wall() enforces a +1px margin when bouncing off container boundaries,
    // so our "home" must match that clamp or the H protocol may never reach phase=0.
    sz_interactive_home_px = (float)XW1 + 0.5f * WALL_THICKNESS + 1.0f;
    sz_interactive_target_px = (float)XW1 + (float)(L0_UNITS * PIXELS_PER_SIGMA);
    sz_interactive_piston_home_px = sz_interactive_target_px;
    sz_interactive_piston_target_px = (float)XW1 + (float)(2.0f * L0_UNITS * PIXELS_PER_SIGMA);
    szilard_refresh_interactive_speed();

    piston_left_x = (float)XW1 - 1000.0f; // keep left piston far away in Szilard mode
    piston_right_x = sz_interactive_piston_home_px;
    vx_piston_left = 0.0f;
    vx_piston_right = 0.0f;
    piston_step_active = false;

    // Divider disabled initially (becomes semipermeable at separation)
    wall_enabled = 0;
    wall_hold_enabled = true;
    wall_is_released = false;
    wall_hold_steps = INT_MAX / 2;
    sz_perm_active = 0;
    sz_perm_species = 0;
    sz_perm_use_memory = 0;
    sz_perm_wall_index = primary_wall_index;
    sz_perm_dual_active = 0;
    sz_fixed_gate_active = 0;
    sz_fixed_gate_species = 0;
    sz_fixed_gate_block_dir = 0;
    sz_fixed_gate_use_memory = 0;
    memset(sz_unstick_hit_count, 0, sizeof(sz_unstick_hit_count));

    // Keep divider at L0 during equilibration; we will move it to the left boundary on separation start
    wall_x = sz_interactive_target_px;
    wall_x_old = wall_x;
    vx_wall = 0.0f;
    if (all_wall_positions) {
        all_wall_positions[primary_wall_index] = wall_x;
    }

    szilard_assign_species(particles_active);
    szilard_clamp_particles_left_half();
}

static void szilard_clamp_particles_left_half(void) {
    const int n = (particles_active > 0) ? particles_active : 0;
    if (n <= 0) return;
    const float R = PARTICLE_RADIUS;
    const float eps = 1.0f;
    const float x_min = (float)XW1 + R + eps;
    const float x_max = (float)XW1 + (float)(L0_UNITS * PIXELS_PER_SIGMA) - R - eps;
    const float y_min = (float)YW1 + R + eps;
    const float y_max = (float)YW2 - R - eps;
    if (!(x_max > x_min) || !(y_max > y_min)) return;

    for (int i = 0; i < n; ++i) {
        const float xi = (float)X[i];
        const float yi = (float)Y[i];
        bool in_bounds = isfinite(xi) && isfinite(yi) && xi >= x_min && xi <= x_max && yi >= y_min && yi <= y_max;
        bool overlaps = false;
        if (in_bounds) {
            for (int j = 0; j < i; ++j) {
                const float dx = xi - (float)X[j];
                const float dy = yi - (float)Y[j];
                const float d2 = dx * dx + dy * dy;
                const float min_d = 2.0f * R;
                if (d2 < (min_d * min_d)) { overlaps = true; break; }
            }
        }
        if (in_bounds && !overlaps) continue;

        for (int it = 0; it < 2000; ++it) {
            const float rx = (float)rand() / (float)RAND_MAX;
            const float ry = (float)rand() / (float)RAND_MAX;
            const float x_try = x_min + rx * (x_max - x_min);
            const float y_try = y_min + ry * (y_max - y_min);
            bool ok = true;
            for (int j = 0; j < i; ++j) {
                const float dx = x_try - (float)X[j];
                const float dy = y_try - (float)Y[j];
                const float d2 = dx * dx + dy * dy;
                const float min_d = 2.0f * R;
                if (d2 < (min_d * min_d)) { ok = false; break; }
            }
            if (ok) { X[i] = x_try; Y[i] = y_try; break; }
        }
    }
}

static void szilard_save_branch_snapshot(void) {
    g_sz_branch_snapshot.valid = 1;
    g_sz_branch_snapshot.particles_active = particles_active;
    g_sz_branch_snapshot.steps_elapsed = steps_elapsed;
    g_sz_branch_snapshot.simulation_started = simulation_started;
    g_sz_branch_snapshot.sz_interactive_phase = (sz_interactive_phase == 1) ? 2 : sz_interactive_phase;
    g_sz_branch_snapshot.wall_enabled = wall_enabled ? 1 : 0;
    g_sz_branch_snapshot.wall_hold_enabled = wall_hold_enabled ? 1 : 0;
    g_sz_branch_snapshot.wall_is_released = wall_is_released ? 1 : 0;
    g_sz_branch_snapshot.wall_hold_steps = wall_hold_steps;
    g_sz_branch_snapshot.sz_perm_active = sz_perm_active;
    g_sz_branch_snapshot.sz_perm_species = sz_perm_species;
    g_sz_branch_snapshot.sz_perm_use_memory = sz_perm_use_memory;
    g_sz_branch_snapshot.sz_perm_wall_index = sz_perm_wall_index;
    g_sz_branch_snapshot.sz_perm_target_side = sz_perm_target_side;
    g_sz_branch_snapshot.sz_perm_dual_active = sz_perm_dual_active;
    g_sz_branch_snapshot.sz_fixed_gate_active = sz_fixed_gate_active;
    g_sz_branch_snapshot.sz_fixed_gate_species = sz_fixed_gate_species;
    g_sz_branch_snapshot.sz_fixed_gate_block_dir = sz_fixed_gate_block_dir;
    g_sz_branch_snapshot.sz_fixed_gate_use_memory = sz_fixed_gate_use_memory;
    g_sz_branch_snapshot.sz_xflip_enabled = sz_xflip_enabled;
    g_sz_branch_snapshot.sz_partial_obs_enabled = sz_partial_obs_enabled;
    g_sz_branch_snapshot.sz_branch_ready_reason = sz_branch_ready_reason;
    g_sz_branch_snapshot.sz_load_enabled = sz_load_enabled;
    g_sz_branch_snapshot.simulation_time = simulation_time;
    g_sz_branch_snapshot.wall_release_time = wall_release_time;
    g_sz_branch_snapshot.sz_xeq_ready_since = sz_xeq_ready_since;
    g_sz_branch_snapshot.sz_memory_erase_heat_min = sz_memory_erase_heat_min;
    g_sz_branch_snapshot.sz_load_eta = sz_load_eta;
    g_sz_branch_snapshot.sz_load_k = sz_load_k;
    g_sz_branch_snapshot.sz_load_fx_start = sz_load_fx_start;
    g_sz_branch_snapshot.sz_load_fx_drop = sz_load_fx_drop;
    g_sz_branch_snapshot.sz_load_z = sz_load_z;
    g_sz_branch_snapshot.sz_load_energy = sz_load_energy;
    g_sz_branch_snapshot.sz_load_work = sz_load_work;
    g_sz_branch_snapshot.sz_load_qdiss = sz_load_qdiss;
    g_sz_branch_snapshot.wall_x = wall_x;
    g_sz_branch_snapshot.wall_x_old = wall_x_old;
    g_sz_branch_snapshot.vx_wall = vx_wall;
    g_sz_branch_snapshot.piston_left_x = piston_left_x;
    g_sz_branch_snapshot.piston_right_x = piston_right_x;
    g_sz_branch_snapshot.vx_piston_left = vx_piston_left;
    g_sz_branch_snapshot.vx_piston_right = vx_piston_right;

    if (particles_active > 0) {
        memcpy(g_sz_branch_snapshot.X, X, (size_t)particles_active * sizeof(X[0]));
        memcpy(g_sz_branch_snapshot.Y, Y, (size_t)particles_active * sizeof(Y[0]));
        memcpy(g_sz_branch_snapshot.Vx, Vx, (size_t)particles_active * sizeof(Vx[0]));
        memcpy(g_sz_branch_snapshot.Vy, Vy, (size_t)particles_active * sizeof(Vy[0]));
        if (sz_species) memcpy(g_sz_branch_snapshot.sz_species, sz_species, (size_t)particles_active * sizeof(int));
        if (sz_memory)  memcpy(g_sz_branch_snapshot.sz_memory,  sz_memory,  (size_t)particles_active * sizeof(int));
    }
}

static int szilard_gate_preferred_push_side(const EDMD_Params* prm, int divider_index, int species, int *target_side_out) {
    if (!prm || !target_side_out) return 0;
    if (divider_index < 0 || divider_index >= EDMD_MAX_DIVIDERS) return 0;

    const int mode = prm->divider_gate_mode[divider_index];
    int target_side = -1;

    if (mode == 1) {
        target_side = (prm->divider_gate_target_side[divider_index] != 0) ? 1 : 0;
        if (species != prm->divider_gate_species[divider_index]) {
            target_side = 1 - target_side; /* blocked species stays on non-target side */
        }
    } else if (mode == 2) {
        if (species == 0) target_side = 0;
        else if (species == 1) target_side = 1;
        else return 0;
    } else if (mode == 3) {
        /* Blocked species stays on the target side; others are transparent. */
        if (species != prm->divider_gate_species[divider_index]) return 0;
        target_side = (prm->divider_gate_target_side[divider_index] != 0) ? 1 : 0;
    } else if (mode == 4) {
        target_side = (prm->divider_gate_target_side[divider_index] != 0) ? 1 : 0;
        if (species != prm->divider_gate_species[divider_index]) {
            target_side = 1 - target_side; /* blocked-unless-hot species stays on non-target side */
        }
    } else {
        return 0;
    }

    *target_side_out = target_side;
    return 1;
}

static int szilard_repair_edmd_particles_for_pistons_and_pockets(void) {
    if (!g_edmd) return 0;
    if (sim_mode != MODE_EDMD) return 0;
    if (cli_experiment_preset != EXPERIMENT_PRESET_SZILARD_ENGINE) return 0;
    if (!sz_interactive_enabled) return 0;

    EDMD_Particle* P = (EDMD_Particle*)edmd_particles(g_edmd);
    const EDMD_Params* ep = edmd_params(g_edmd);
    if (!P || !ep) return 0;

    const double R = ep->radius;
    const double eps = fmax(1e-9, 1e-6 * fmax(R, 1.0));
    int changed = 0;

    double x_min = R;
    double x_max = ep->boxW - R;
    if (ep->has_pistonL) x_min = fmax(x_min, ep->pistonL_x + R + eps);
    if (ep->has_pistonR) x_max = fmin(x_max, ep->pistonR_x - R - eps);
    if (x_max < x_min) {
        const double mid = 0.5 * (x_min + x_max);
        x_min = mid;
        x_max = mid;
    }

    for (int i = 0; i < ep->N; ++i) {
        if (P[i].x < x_min) {
            if (ep->has_pistonL && P[i].vx < ep->pistonL_vx) {
                P[i].vx = 2.0 * ep->pistonL_vx - P[i].vx;
            }
            P[i].x = x_min;
            P[i].coll_count++;
            changed = 1;
        } else if (P[i].x > x_max) {
            if (ep->has_pistonR && P[i].vx > ep->pistonR_vx) {
                P[i].vx = 2.0 * ep->pistonR_vx - P[i].vx;
            }
            P[i].x = x_max;
            P[i].coll_count++;
            changed = 1;
        }

        if (P[i].y < R) {
            if (P[i].vy < 0.0) P[i].vy = -P[i].vy;
            P[i].y = R;
            P[i].coll_count++;
            changed = 1;
        } else if (P[i].y > ep->boxH - R) {
            if (P[i].vy > 0.0) P[i].vy = -P[i].vy;
            P[i].y = ep->boxH - R;
            P[i].coll_count++;
            changed = 1;
        }
    }

    if (ep->divider_count > 0 && ep->divider_thickness[0] > 0.0 &&
        (sz_interactive_phase == 3 || sz_interactive_phase == 4 || sz_interactive_phase == 5) &&
        ep->divider_vx[0] < 0.0)
    {
        const double left_face = ep->divider_x[0] - 0.5 * ep->divider_thickness[0] - R;
        const double pocket_w = left_face - x_min;
        if (pocket_w <= (2.5 * R + 2.0 * eps)) {
            const double right_face = ep->divider_x[0] + 0.5 * ep->divider_thickness[0] + R + eps;
            for (int i = 0; i < ep->N; ++i) {
                int target_side = -1;
                if (!szilard_particle_moving_gate_target_side(i, &target_side) || target_side != 1) continue;
                if (P[i].x <= left_face + eps) {
                    P[i].x = fmin(right_face, x_max);
                    P[i].coll_count++;
                    changed = 1;
                }
            }
        }
    }

    if (ep->divider_count > 1 && ep->divider_thickness[1] > 0.0 &&
        sz_fixed_gate_active && ep->has_pistonR && ep->pistonR_vx < 0.0)
    {
        const double gate_x = ep->divider_x[1];
        const double gate_right_face = gate_x + 0.5 * ep->divider_thickness[1] + R;
        const double gate_left_face = gate_x - 0.5 * ep->divider_thickness[1] - R - eps;
        const double pocket_w = (ep->pistonR_x - R - eps) - gate_right_face;
        if (pocket_w <= (2.5 * R + 2.0 * eps)) {
            for (int i = 0; i < ep->N; ++i) {
                if (!szilard_particle_can_cross_fixed_gate_right_to_left(i)) continue;
                if (P[i].x >= gate_right_face - eps) {
                    P[i].x = fmax(gate_left_face, x_min);
                    P[i].coll_count++;
                    changed = 1;
                }
            }
        }
    }

    return changed;
}

static void szilard_rebuild_edmd_from_current_globals(void) {
    if (sim_mode != MODE_EDMD || cli_experiment_preset != EXPERIMENT_PRESET_SZILARD_ENGINE) return;

    if (g_edmd) { edmd_destroy(g_edmd); g_edmd = NULL; }

    EDMD_Params prm = {0};
    prm.boxW = (double)(XW2 - XW1);
    prm.boxH = (double)(YW2 - YW1);
    prm.radius = (double)PARTICLE_RADIUS;
    prm.N = particles_active;
    prm.cell_size = 0.0;
    prm.divider_count = 2;
    prm.heatbath_enabled = heatbath_enabled;
    prm.heatbath_temperature = (double)heatbath_temperature;
    prm.thermal_wall_mode = thermal_wall_mode;
    prm.mb_overshoot_factor = (double)mb_overshoot_factor;
    prm.stability_window_percent = (double)stability_window_percent;
    prm.particle_mass = (double)PARTICLE_MASS;
    prm.kB = (double)kB_effective();
    prm.species = sz_species;

    g_edmd = edmd_create(&prm);
    if (!g_edmd) {
        fprintf(stderr, "Szilard EDMD recreate failed.\n");
        exit(1);
    }

    EDMD_Particle* P = (EDMD_Particle*)edmd_particles(g_edmd);
    for (int i = 0; i < particles_active; ++i) {
        P[i].x = (double)(X[i] - XW1);
        P[i].y = (double)(Y[i] - YW1);
        P[i].vx = (double)Vx[i];
        P[i].vy = (double)Vy[i];
        P[i].coll_count = 0;
    }

    szilard_sync_edmd_from_globals();
    edmd_reset_work(g_edmd);
}

static void szilard_restore_branch_snapshot(void) {
    if (!g_sz_branch_snapshot.valid) return;

    particles_active = g_sz_branch_snapshot.particles_active;
    steps_elapsed = g_sz_branch_snapshot.steps_elapsed;
    simulation_started = g_sz_branch_snapshot.simulation_started;
    sz_interactive_phase = g_sz_branch_snapshot.sz_interactive_phase;
    wall_enabled = g_sz_branch_snapshot.wall_enabled;
    wall_hold_enabled = g_sz_branch_snapshot.wall_hold_enabled;
    wall_is_released = g_sz_branch_snapshot.wall_is_released;
    wall_hold_steps = g_sz_branch_snapshot.wall_hold_steps;
    sz_perm_active = g_sz_branch_snapshot.sz_perm_active;
    sz_perm_species = g_sz_branch_snapshot.sz_perm_species;
    sz_perm_use_memory = g_sz_branch_snapshot.sz_perm_use_memory;
    sz_perm_wall_index = g_sz_branch_snapshot.sz_perm_wall_index;
    sz_perm_target_side = g_sz_branch_snapshot.sz_perm_target_side;
    sz_perm_dual_active = g_sz_branch_snapshot.sz_perm_dual_active;
    sz_fixed_gate_active = g_sz_branch_snapshot.sz_fixed_gate_active;
    sz_fixed_gate_species = g_sz_branch_snapshot.sz_fixed_gate_species;
    sz_fixed_gate_block_dir = g_sz_branch_snapshot.sz_fixed_gate_block_dir;
    sz_fixed_gate_use_memory = g_sz_branch_snapshot.sz_fixed_gate_use_memory;
    sz_xflip_enabled = g_sz_branch_snapshot.sz_xflip_enabled;
    sz_partial_obs_enabled = g_sz_branch_snapshot.sz_partial_obs_enabled;
    sz_branch_ready_reason = g_sz_branch_snapshot.sz_branch_ready_reason;
    sz_load_enabled = g_sz_branch_snapshot.sz_load_enabled;
    sz_partial_obs_ref_imem_gate_bits = 1.0;
    sz_partial_obs_peak_uncertainty_frac = 0.0;
    simulation_time = g_sz_branch_snapshot.simulation_time;
    wall_release_time = g_sz_branch_snapshot.wall_release_time;
    sz_xeq_ready_since = g_sz_branch_snapshot.sz_xeq_ready_since;
    sz_memory_erase_heat_min = g_sz_branch_snapshot.sz_memory_erase_heat_min;
    sz_load_eta = g_sz_branch_snapshot.sz_load_eta;
    sz_load_k = g_sz_branch_snapshot.sz_load_k;
    sz_load_fx_start = g_sz_branch_snapshot.sz_load_fx_start;
    sz_load_fx_drop = g_sz_branch_snapshot.sz_load_fx_drop;
    sz_load_z = g_sz_branch_snapshot.sz_load_z;
    sz_load_energy = g_sz_branch_snapshot.sz_load_energy;
    sz_load_work = g_sz_branch_snapshot.sz_load_work;
    sz_load_qdiss = g_sz_branch_snapshot.sz_load_qdiss;
    wall_x = g_sz_branch_snapshot.wall_x;
    wall_x_old = g_sz_branch_snapshot.wall_x_old;
    vx_wall = g_sz_branch_snapshot.vx_wall;
    piston_left_x = g_sz_branch_snapshot.piston_left_x;
    piston_right_x = g_sz_branch_snapshot.piston_right_x;
    vx_piston_left = g_sz_branch_snapshot.vx_piston_left;
    vx_piston_right = g_sz_branch_snapshot.vx_piston_right;

    if (particles_active > 0) {
        memcpy(X, g_sz_branch_snapshot.X, (size_t)particles_active * sizeof(X[0]));
        memcpy(Y, g_sz_branch_snapshot.Y, (size_t)particles_active * sizeof(Y[0]));
        memcpy(Vx, g_sz_branch_snapshot.Vx, (size_t)particles_active * sizeof(Vx[0]));
        memcpy(Vy, g_sz_branch_snapshot.Vy, (size_t)particles_active * sizeof(Vy[0]));
        if (sz_species) memcpy(sz_species, g_sz_branch_snapshot.sz_species, (size_t)particles_active * sizeof(int));
        if (sz_memory)  memcpy(sz_memory,  g_sz_branch_snapshot.sz_memory,  (size_t)particles_active * sizeof(int));
    }

    if (all_wall_positions) {
        all_wall_positions[primary_wall_index] = wall_x;
    }

    piston_work_left = 0.0;
    piston_work_right = 0.0;
    piston_work_baseline = 0.0;
    piston_work_delta_max = 0.0;
    piston_push_tracking = false;
    snprintf(sz_log_branch_label, sizeof(sz_log_branch_label), "restore");
    sz_log_reset_requested = 1;

    szilard_rebuild_edmd_from_current_globals();
}

static inline double szilard_safe_log(double x) { return (x > 0.0) ? log(x) : 0.0; }

static double szilard_mutual_information_nats(int nL0, int nL1, int nR0, int nR1) {
    const double N = (double)(nL0 + nL1 + nR0 + nR1);
    if (N <= 0.0) return 0.0;

    const double pL0 = (double)nL0 / N;
    const double pL1 = (double)nL1 / N;
    const double pR0 = (double)nR0 / N;
    const double pR1 = (double)nR1 / N;
    const double pL  = pL0 + pL1;
    const double pR  = pR0 + pR1;
    const double p0  = pL0 + pR0;
    const double p1  = pL1 + pR1;

    double I = 0.0;
    if (pL0 > 0.0 && pL > 0.0 && p0 > 0.0) I += pL0 * (szilard_safe_log(pL0) - szilard_safe_log(pL) - szilard_safe_log(p0));
    if (pL1 > 0.0 && pL > 0.0 && p1 > 0.0) I += pL1 * (szilard_safe_log(pL1) - szilard_safe_log(pL) - szilard_safe_log(p1));
    if (pR0 > 0.0 && pR > 0.0 && p0 > 0.0) I += pR0 * (szilard_safe_log(pR0) - szilard_safe_log(pR) - szilard_safe_log(p0));
    if (pR1 > 0.0 && pR > 0.0 && p1 > 0.0) I += pR1 * (szilard_safe_log(pR1) - szilard_safe_log(pR) - szilard_safe_log(p1));
    return I;
}

static double szilard_mutual_information_3state_nats(int nL0, int nL1, int nLU,
                                                     int nR0, int nR1, int nRU) {
    const double N = (double)(nL0 + nL1 + nLU + nR0 + nR1 + nRU);
    if (!(N > 0.0)) return 0.0;

    const double pL  = (double)(nL0 + nL1 + nLU) / N;
    const double pR  = (double)(nR0 + nR1 + nRU) / N;

    const double p0  = (double)(nL0 + nR0) / N;
    const double p1  = (double)(nL1 + nR1) / N;
    const double pU  = (double)(nLU + nRU) / N;

    const double pL0 = (double)nL0 / N;
    const double pL1 = (double)nL1 / N;
    const double pLU = (double)nLU / N;
    const double pR0 = (double)nR0 / N;
    const double pR1 = (double)nR1 / N;
    const double pRU = (double)nRU / N;

    double I = 0.0;
    if (pL0 > 0.0 && pL > 0.0 && p0 > 0.0) I += pL0 * (szilard_safe_log(pL0) - szilard_safe_log(pL) - szilard_safe_log(p0));
    if (pL1 > 0.0 && pL > 0.0 && p1 > 0.0) I += pL1 * (szilard_safe_log(pL1) - szilard_safe_log(pL) - szilard_safe_log(p1));
    if (pLU > 0.0 && pL > 0.0 && pU > 0.0) I += pLU * (szilard_safe_log(pLU) - szilard_safe_log(pL) - szilard_safe_log(pU));
    if (pR0 > 0.0 && pR > 0.0 && p0 > 0.0) I += pR0 * (szilard_safe_log(pR0) - szilard_safe_log(pR) - szilard_safe_log(p0));
    if (pR1 > 0.0 && pR > 0.0 && p1 > 0.0) I += pR1 * (szilard_safe_log(pR1) - szilard_safe_log(pR) - szilard_safe_log(p1));
    if (pRU > 0.0 && pR > 0.0 && pU > 0.0) I += pRU * (szilard_safe_log(pRU) - szilard_safe_log(pR) - szilard_safe_log(pU));
    return I;
}

static double szilard_side_purity(int n0, int n1) {
    const int nt = n0 + n1;
    if (nt <= 0) return 0.0;
    return (double)((n0 > n1) ? n0 : n1) / (double)nt;
}

static double szilard_separation_score(int nL0, int nL1, int nR0, int nR1) {
    const int nt = nL0 + nL1 + nR0 + nR1;
    if (nt <= 0) return 0.0;
    return (double)(((nL0 > nL1) ? nL0 : nL1) + ((nR0 > nR1) ? nR0 : nR1)) / (double)nt;
}

static const char* szilard_phase_label(int phase) {
    switch (phase) {
        case 1: return "G";
        case 2: return "Xeq";
        case 3: return "H";
        case 4: return "J";
        case 5: return "PO";
        case 6: return "E?";
        case 7: return "E";
        default: return "Eq";
    }
}

static const char* szilard_particle_view_mode_label(int mode) {
    switch (mode) {
        case SZ_VIEW_SPECIES: return "species";
        case SZ_VIEW_MEMORY: return "memory";
        case SZ_VIEW_OVERLAY: return "overlay";
        default: return "unknown";
    }
}

static const char* szilard_branch_ready_reason_label(int reason) {
    switch (reason) {
        case 1: return "full Xeq";
        case 2: return "5 mismatches";
        default: return "waiting";
    }
}

static double szilard_reference_kbt(void) {
    const double Tref = heatbath_enabled ? (double)heatbath_temperature : (double)temperature_runtime;
    return (double)kB_effective() * Tref;
}

static double szilard_available_free_energy_from_nats(double info_nats, int total_count) {
    if (!(info_nats > 0.0) || total_count <= 0) return 0.0;
    return (double)total_count * szilard_reference_kbt() * info_nats;
}

static double szilard_current_visible_free_energy_from_counts(int nL0, int nL1, int nR0, int nR1) {
    const double info_nats = szilard_mutual_information_nats(nL0, nL1, nR0, nR1);
    return szilard_available_free_energy_from_nats(info_nats, nL0 + nL1 + nR0 + nR1);
}

static void szilard_reset_visible_load_state(double fx_start) {
    sz_load_fx_start = (fx_start > 0.0) ? fx_start : 0.0;
    sz_load_fx_drop = 0.0;
    sz_load_z = 0.0;
    sz_load_energy = 0.0;
    sz_load_work = 0.0;
    sz_load_qdiss = 0.0;
}

static void szilard_update_visible_load_state(double fx_now) {
    if (!sz_load_enabled) {
        szilard_reset_visible_load_state(0.0);
        return;
    }

    double drop = sz_load_fx_start - fx_now;
    if (!(drop > 0.0)) drop = 0.0;
    sz_load_fx_drop = drop;

    const double eta = (sz_load_eta < 0.0) ? 0.0 : ((sz_load_eta > 1.0) ? 1.0 : sz_load_eta);
    sz_load_work = eta * drop;
    sz_load_energy = sz_load_work;
    sz_load_qdiss = drop - sz_load_work;
    if (sz_load_qdiss < 0.0) sz_load_qdiss = 0.0;

    if (sz_load_k > 1e-12 && sz_load_energy > 0.0) {
        sz_load_z = sqrt(2.0 * sz_load_energy / sz_load_k);
    } else {
        sz_load_z = 0.0;
    }
}

static void szilard_species_side_counts(int *nL0, int *nL1, int *nR0, int *nR1) {
    if (nL0) *nL0 = 0;
    if (nL1) *nL1 = 0;
    if (nR0) *nR0 = 0;
    if (nR1) *nR1 = 0;
    if (!sz_species) return;

    for (int i = 0; i < particles_active; ++i) {
        const int side = szilard_count_side_for_perfect_sep(i); /* 0=left, 1=right */
        const int sp = (sz_species[i] != 0) ? 1 : 0;
        if (side == 0) {
            if (sp == 0) (*nL0)++; else (*nL1)++;
        } else {
            if (sp == 0) (*nR0)++; else (*nR1)++;
        }
    }
}

static void szilard_sync_edmd_from_globals(void) {
    if (!g_edmd) return;
    if (sim_mode != MODE_EDMD) return;
    if (cli_experiment_preset != EXPERIMENT_PRESET_SZILARD_ENGINE) return;
    if (!sz_interactive_enabled) return;

    const double gate_thick = fmax(1.0, 0.25 * (double)WALL_THICKNESS);
    double cx[2] = {
        (double)wall_x - (double)XW1,
        (double)sz_interactive_target_px - (double)XW1
    };
    const int suppress_fixed_gate_overlap =
        (wall_enabled && sz_fixed_gate_active &&
         fabs(cx[0] - cx[1]) <= (0.5 * ((double)WALL_THICKNESS + gate_thick) + 2.0 * (double)PARTICLE_RADIUS + 1e-6));
    double th[2] = {
        wall_enabled ? (double)WALL_THICKNESS : 0.0,
        (wall_enabled && sz_fixed_gate_active && !sz_perm_dual_active && !suppress_fixed_gate_overlap) ? gate_thick : 0.0
    };
    double mass[2] = { 0.0, 0.0 }; /* prescribed motion */
    double vx[2] = { (double)vx_wall, 0.0 };
    double k[2] = { 0.0, 0.0 };
    double xeq[2] = { 0.0, 0.0 };

    edmd_config_dividers(g_edmd, 2, cx, th);
    edmd_set_divider_motions(g_edmd, 2, mass, vx);
    edmd_set_divider_springs(g_edmd, 2, k, xeq);

    /* Disable left piston entirely in Szilard; keep right piston synced. */
    const double xR = (double)piston_right_x - (double)XW1;
    edmd_config_pistons(g_edmd, 0, 0.0, 0.0, 0.0, 1, xR, (double)vx_piston_right, 0.0);

    EDMD_Params* prm = (EDMD_Params*)edmd_params(g_edmd);
    if (sz_partial_obs_enabled || sz_perm_use_memory || sz_fixed_gate_use_memory) {
        szilard_refresh_gate_labels();
        prm->species = sz_gate_labels;
    } else {
        prm->species = sz_species;
    }

    /* ##CHRIS: PP collisions on by default; disable during G if --szilard-no-pp was passed */
    prm->pp_collisions_enabled = (cli_szilard_no_pp_during_g && sz_interactive_phase == 1) ? 0 : 1;

    for (int d = 0; d < EDMD_MAX_DIVIDERS; ++d) {
        prm->divider_gate_mode[d] = 0;
        prm->divider_gate_species[d] = 0;
        prm->divider_gate_target_side[d] = 0;
        prm->divider_gate_hot_ratio[d] = 0.0;
        prm->divider_gate_speed_ref[d] = 0.0;
    }

    double mean_speed = 0.0, std_speed = 0.0;
    szilard_speed_stats(&mean_speed, &std_speed);
    const double speed_ref = (mean_speed > 1e-12) ? mean_speed : fmax(std_speed, 1e-12);

    /* Moving divider gate (index 0). */
    if (wall_enabled) {
        if (sz_perm_dual_active) {
            prm->divider_gate_mode[0] = 2;
        } else if (sz_perm_active) {
            prm->divider_gate_mode[0] = 1;
            prm->divider_gate_species[0] = sz_perm_species;
            prm->divider_gate_target_side[0] = (sz_perm_target_side != 0) ? 1 : 0;
            if (sz_partial_obs_enabled) {
                prm->divider_gate_hot_ratio[0] = sz_gate_hot_speed_mult;
                prm->divider_gate_speed_ref[0] = speed_ref;
            }
        }
    }

    /* Fixed L0 gate (index 1), only physical when th[1] > 0. */
    if (th[1] > 0.0 && sz_fixed_gate_active) {
        if (sz_interactive_phase == 1 && !sz_fixed_gate_use_memory) {
            /* During G, keep the original dual one-way sorter at L0 in EDMD. */
            prm->divider_gate_mode[1] = 2;
        } else {
            /* Single blocked-class fixed gate, optionally driven by memory labels. */
            prm->divider_gate_mode[1] = 3;
            prm->divider_gate_species[1] = sz_fixed_gate_species;
            prm->divider_gate_target_side[1] = (sz_fixed_gate_block_dir == 0) ? 0 : 1;
            if (sz_partial_obs_enabled) {
                prm->divider_gate_hot_ratio[1] = sz_gate_hot_speed_mult;
                prm->divider_gate_speed_ref[1] = speed_ref;
            }
        }
    }

    /* Push any particles out of newly-enabled divider slabs without reflecting (Szilard permeability). */
    if (wall_enabled) {
        EDMD_Particle* P = (EDMD_Particle*)edmd_particles(g_edmd);
        const EDMD_Params* ep = edmd_params(g_edmd);
        for (int d = 0; d < ep->divider_count && d < 2; ++d) {
            if (ep->divider_thickness[d] <= 0.0) continue;
            const double L = ep->divider_x[d] - 0.5 * ep->divider_thickness[d];
            const double Rf = ep->divider_x[d] + 0.5 * ep->divider_thickness[d];
            const double Rpart = ep->radius;
            const double left_face = L - Rpart;
            const double right_face = Rf + Rpart;
            for (int i = 0; i < ep->N; ++i) {
                if (P[i].x >= left_face && P[i].x <= right_face) {
                    int target_side = -1;
                    const int species = (ep->species ? ep->species[i] : 0);
                    if (szilard_gate_preferred_push_side(ep, d, species, &target_side)) {
                        if (target_side == 0) P[i].x = left_face - 1e-9;
                        else                  P[i].x = right_face + 1e-9;
                    } else {
                        const double gapL = fabs(P[i].x - left_face);
                        const double gapR = fabs(right_face - P[i].x);
                        if (gapL < gapR) P[i].x = left_face - 1e-9;
                        else             P[i].x = right_face + 1e-9;
                    }
                    P[i].coll_count++;
                }
            }
        }
    }

    // The slab push-out above is not an EDMD event and can introduce PP overlaps.
    // Also repair particles that numerically ended beyond the moving right piston or
    // in a vanishing pocket where the current gate should let them escape.
    (void)szilard_repair_edmd_particles_for_pistons_and_pockets();

    // Relax overlaps (Szilard-only) and rebuild the event schedule.
    edmd_relax_pp_overlaps_szilard(g_edmd);
    edmd_reschedule_all(g_edmd);
}

static void szilard_begin_separation(void) {
    if (!sz_interactive_enabled) return;
    if (sz_interactive_phase == 6) {
        if (!cli_quiet) {
            printf("Szilard: this branch is waiting for manual erase. Press 'e' before starting a new cycle.\n");
        }
        return;
    }
    if (sz_interactive_phase != 0 && sz_interactive_phase != 7) {
        if (!cli_quiet) {
            printf("Szilard: G can only start from equilibrium or after explicit erase.\n");
        }
        return;
    }
    if (!szilard_memory_ready()) {
        if (!cli_quiet) {
            printf("Szilard: measure first with 'm' so y is defined before G.\n");
        }
        return;
    }

    // Move divider to the left boundary before sweeping right
    wall_x = sz_interactive_home_px;
    wall_x_old = wall_x;
    if (all_wall_positions) {
        all_wall_positions[primary_wall_index] = wall_x;
    }
    wall_enabled = 1;
    wall_hold_enabled = false;  // allow divider to sweep right
    wall_is_released = true;
    wall_hold_steps = 0;
    sz_perm_active = 1;
    sz_perm_dual_active = 0;
    sz_perm_species = 0;        // red passes (remains left as divider sweeps right)
    sz_perm_use_memory = 0;
    sz_perm_target_side = 0;    // keep permeable species on left
    sz_perm_wall_index = primary_wall_index;
    sz_fixed_gate_active = 1;
    sz_fixed_gate_species = 0;   // red gate
    sz_fixed_gate_block_dir = 0; // block red L->R at L0
    sz_fixed_gate_use_memory = 0;
    vx_wall = sz_interactive_vx_px;
    vx_piston_right = sz_interactive_vx_px;
    sz_interactive_phase = 1;
    snprintf(sz_log_branch_label, sizeof(sz_log_branch_label), "G");
    sz_log_reset_requested = 1;
    sz_sep_perfect_since = -1.0;
    sz_xflip_enabled = 0;
    sz_xeq_ready_since = -1.0;
    sz_partial_obs_enabled = 0;
    sz_partial_obs_ref_imem_gate_bits = 1.0;
    sz_partial_obs_peak_uncertainty_frac = 0.0;
    sz_branch_ready_reason = 0;
    sz_pending_recombine_mode = 0;
    szilard_reset_visible_load_state(0.0);
    g_sz_branch_snapshot.valid = 0;
    simulation_started = 1;
    szilard_sync_edmd_from_globals();
    if (!cli_quiet) printf("▶ Szilard: separation started (G)\n");
}

static void szilard_begin_recombine(void) {
    if (!sz_interactive_enabled) return;
    if (!sz_stage1_baseline_mode && sz_interactive_phase == 2 && !g_sz_branch_snapshot.valid) {
        sz_pending_recombine_mode = 1;
        if (!cli_quiet) {
            printf("Szilard: H queued. Waiting for the saved post-Xeq branch point to be ready.\n");
        }
        return;
    }
    if (sz_interactive_phase != 2 || !g_sz_branch_snapshot.valid) {
        if (!cli_quiet) {
            printf("Szilard: H needs the saved post-Xeq branch point. Run G, wait for Xeq, or restore with 'n'.\n");
        }
        return;
    }
    // Recombine protocol (H):
    // Divider sweeps left (L0 -> 0). The already-separated red gas on the left should pass
    // into the blue side (left->right through the moving wall) so the species remix while the
    // divider returns home.
    sz_perm_active = 1;
    sz_perm_dual_active = 0;
    sz_perm_species = 0;        // memory-0 passes through moving wall
    sz_perm_use_memory = 1;
    sz_perm_target_side = 1;    // allow pass toward right
    sz_perm_wall_index = primary_wall_index;
    sz_fixed_gate_active = 1;
    sz_fixed_gate_species = 1;   // memory-1 crosses the fixed L0 gate
    sz_fixed_gate_block_dir = 0; // block left->right, allow right->left into the final compartment
    sz_fixed_gate_use_memory = 1;
    wall_enabled = 1;
    wall_hold_enabled = false;
    wall_is_released = true;
    wall_hold_steps = 0;
    vx_wall = -sz_interactive_vx_px;
    vx_piston_right = -sz_interactive_vx_px;
    sz_xflip_enabled = 0;
    sz_xeq_ready_since = -1.0;
    sz_partial_obs_enabled = 0;
    sz_partial_obs_ref_imem_gate_bits = 1.0;
    sz_partial_obs_peak_uncertainty_frac = 0.0;
    sz_pending_recombine_mode = 0;
    sz_interactive_phase = 3;
    sz_log_branch_id += 1;
    snprintf(sz_log_branch_label, sizeof(sz_log_branch_label), "H");
    sz_log_reset_requested = 1;
    simulation_started = 1;
    szilard_sync_edmd_from_globals();
    if (!cli_quiet) printf("◀ Szilard: memory-gated recombine started (H)\n");
}

static void szilard_begin_recombine_wrong(void) {
    if (!sz_interactive_enabled) return;
    if (!sz_stage1_baseline_mode && sz_interactive_phase == 2 && !g_sz_branch_snapshot.valid) {
        sz_pending_recombine_mode = 2;
        if (!cli_quiet) {
            printf("Szilard: J queued. Waiting for the saved post-Xeq branch point to be ready.\n");
        }
        return;
    }
    if (sz_interactive_phase != 2 || !g_sz_branch_snapshot.valid) {
        if (!cli_quiet) {
            printf("Szilard: J needs the saved post-Xeq branch point. Run G, wait for Xeq, or restore with 'n'.\n");
        }
        return;
    }
    // "Wrong wall" recombine protocol (J): intentionally choose the opposite permeability
    // for the moving gate, but now with respect to the stored memory bit y.
    sz_perm_active = 1;
    sz_perm_dual_active = 0;
    sz_perm_species = 1;        // WRONG: memory-1 passes through moving wall
    sz_perm_use_memory = 1;
    sz_perm_target_side = 1;    // allow pass toward right
    sz_perm_wall_index = primary_wall_index;
    sz_fixed_gate_active = 1;
    sz_fixed_gate_species = 0;   // complementary memory class crosses the fixed L0 gate
    sz_fixed_gate_block_dir = 0; // block left->right, allow right->left into the final compartment
    sz_fixed_gate_use_memory = 1;
    wall_enabled = 1;
    wall_hold_enabled = false;
    wall_is_released = true;
    wall_hold_steps = 0;
    vx_wall = -sz_interactive_vx_px;
    vx_piston_right = -sz_interactive_vx_px;
    sz_xflip_enabled = 0;
    sz_xeq_ready_since = -1.0;
    sz_partial_obs_enabled = 0;
    sz_partial_obs_ref_imem_gate_bits = 1.0;
    sz_partial_obs_peak_uncertainty_frac = 0.0;
    sz_pending_recombine_mode = 0;
    sz_interactive_phase = 4;
    sz_log_branch_id += 1;
    snprintf(sz_log_branch_label, sizeof(sz_log_branch_label), "J");
    sz_log_reset_requested = 1;
    simulation_started = 1;
    szilard_sync_edmd_from_globals();
    if (!cli_quiet) printf("◀ Szilard: wrong memory-gated recombine started (J)\n");
}

static void szilard_begin_recombine_partial(void) {
    if (!sz_interactive_enabled) return;
    if (sz_stage1_baseline_mode) {
        if (!cli_quiet) {
            printf("Szilard: O is disabled in the Stage-1 baseline. Finish the clean hidden-memory result first.\n");
        }
        return;
    }
    if (sz_interactive_phase == 2 && !g_sz_branch_snapshot.valid) {
        sz_pending_recombine_mode = 3;
        if (!cli_quiet) {
            printf("Szilard: O queued. Waiting for the saved post-Xeq branch point to be ready.\n");
        }
        return;
    }
    if (sz_interactive_phase != 2 || !g_sz_branch_snapshot.valid) {
        if (!cli_quiet) {
            printf("Szilard: O needs the saved post-Xeq branch point. Run G, wait for Xeq, or restore with 'n'.\n");
        }
        return;
    }
    /* Partially observable recombine:
       same geometry as H, but once x drifts away from stored y the gates no
       longer act on a perfectly readable memory label. */
    sz_perm_active = 1;
    sz_perm_dual_active = 0;
    sz_perm_species = 0;        // memory-0 / branch-intended class toward the right
    sz_perm_use_memory = 1;
    sz_perm_target_side = 1;
    sz_perm_wall_index = primary_wall_index;
    sz_fixed_gate_active = 1;
    sz_fixed_gate_species = 1;   // memory-1 into the final left compartment
    sz_fixed_gate_block_dir = 0;
    sz_fixed_gate_use_memory = 1;
    wall_enabled = 1;
    wall_hold_enabled = false;
    wall_is_released = true;
    wall_hold_steps = 0;
    vx_wall = -sz_interactive_vx_px;
    vx_piston_right = -sz_interactive_vx_px;
    sz_xflip_enabled = 1;
    sz_xeq_ready_since = -1.0;
    sz_partial_obs_enabled = 1;
    sz_partial_obs_ref_imem_gate_bits = szilard_memory_gate_information_bits_at_split(wall_x);
    if (sz_partial_obs_ref_imem_gate_bits < 1e-9) sz_partial_obs_ref_imem_gate_bits = 1.0;
    sz_partial_obs_peak_uncertainty_frac = 0.0;
    sz_pending_recombine_mode = 0;
    sz_interactive_phase = 5;
    sz_log_branch_id += 1;
    snprintf(sz_log_branch_label, sizeof(sz_log_branch_label), "O");
    sz_log_reset_requested = 1;
    simulation_started = 1;
    szilard_sync_edmd_from_globals();
    if (!cli_quiet) {
        printf("◀ Szilard: partially observable memory-gated recombine started (O)\n");
    }
}

static void szilard_update_interactive_poststep(void) {
    if (!sz_interactive_enabled) return;
    const float eps = 0.5f;
    if (sz_interactive_phase == 1) { // separating
        bool piston_done = false;
        bool wall_done = false;
        if (piston_right_x >= sz_interactive_piston_target_px - eps) {
            piston_right_x = sz_interactive_piston_target_px;
            vx_piston_right = 0.0f;
            piston_done = true;
        }
        if (wall_x >= sz_interactive_target_px - eps) {
            wall_x = sz_interactive_target_px;
            vx_wall = 0.0f;
            wall_done = true;
        }
        if (piston_done && wall_done) {
            wall_hold_enabled = true;
            wall_is_released = false;
            wall_hold_steps = INT_MAX / 2;
            vx_wall = 0.0f;
            vx_piston_right = 0.0f;
            /* End-of-G hold should use a single stationary dual one-way membrane at L0:
               species0 may only go left, species1 may only go right.
               This lets the last stragglers finish sorting without allowing wrong-direction leaks. */
            sz_perm_active = 0;
            sz_perm_dual_active = 1;
            sz_perm_use_memory = 0;
            sz_fixed_gate_active = 0;
            szilard_sync_edmd_from_globals();

            int nL0=0,nL1=0,nR0=0,nR1=0;
            szilard_species_side_counts(&nL0,&nL1,&nR0,&nR1);
            const int perfect = (nL1 == 0 && nR0 == 0 && nL0 > 0 && nR1 > 0);
            if (perfect) {
                if (sz_sep_perfect_since < 0.0) {
                    sz_sep_perfect_since = simulation_time;
                    if (!cli_quiet) {
                        printf("📌 Szilard: perfect x-separation reached; waiting %.1f sigma-time for stability\n",
                               sz_sep_perfect_hold_time);
                    }
                }
                if ((simulation_time - sz_sep_perfect_since) >= sz_sep_perfect_hold_time) {
                    sz_interactive_phase = 2;
                    wall_x = sz_interactive_target_px;
                    wall_x_old = wall_x;
                    wall_enabled = 1;
                    wall_hold_enabled = true;
                    wall_is_released = false;
                    wall_hold_steps = INT_MAX / 2;
                    vx_wall = 0.0f;
                    piston_right_x = sz_interactive_piston_target_px;
                    vx_piston_right = 0.0f;
                    sz_perm_active = 0;
                    sz_perm_dual_active = 0;
                    sz_perm_use_memory = 0;
                    sz_fixed_gate_active = 0;
                    sz_fixed_gate_use_memory = 0;
                    sz_xflip_enabled = 1;
                    sz_xeq_ready_since = -1.0;
                    sz_partial_obs_enabled = 0;
                    sz_partial_obs_ref_imem_gate_bits = 1.0;
                    sz_partial_obs_peak_uncertainty_frac = 0.0;
                    szilard_reset_visible_load_state(
                        szilard_current_visible_free_energy_from_counts(nL0, nL1, nR0, nR1)
                    );
                    g_sz_branch_snapshot.valid = 0;
                    /* ##CHRIS: re-enable PP collisions now that G is done (was disabled by --szilard-no-pp) */
                    if (cli_szilard_no_pp_during_g && g_edmd) {
                        EDMD_Params* prm2 = (EDMD_Params*)edmd_params(g_edmd);
                        prm2->pp_collisions_enabled = 1;
                        edmd_reschedule_all(g_edmd);
                    }
                    szilard_sync_edmd_from_globals();
                    if (!cli_quiet) {
                        if (sz_stage1_baseline_mode) {
                            printf("⏸️  Szilard: entering Xeq. x may now flip while y stays fixed; waiting for full visible-state equilibration.\n");
                        } else {
                            printf("⏸️  Szilard: entering Xeq. x may now flip while y stays fixed; waiting for species-side information to decay.\n");
                        }
                    }
                }
            } else {
                sz_sep_perfect_since = -1.0;
            }
        }
    } else if (sz_interactive_phase == 2) {
        if (!g_sz_branch_snapshot.valid) {
            int nL0=0,nL1=0,nR0=0,nR1=0;
            int xy_match=0, xy_mismatch=0, xy_unknown=0;
            szilard_species_side_counts(&nL0,&nL1,&nR0,&nR1);
            szilard_xy_agreement_counts(&xy_match, &xy_mismatch, &xy_unknown);
            const double I_nats = szilard_mutual_information_nats(nL0, nL1, nR0, nR1);
            const double I_bits = I_nats / log(2.0);
            const double F_x_now = szilard_available_free_energy_from_nats(I_nats, nL0 + nL1 + nR0 + nR1);
            const int left_count = nL0 + nL1;
            const int right_count = nR0 + nR1;
            const int info_ready = sz_stage1_baseline_mode
                                 ? (fabs(I_nats) <= 1e-12)
                                 : (I_bits <= sz_xeq_info_threshold_bits);
            const int mismatch_ready = (!sz_stage1_baseline_mode &&
                                        xy_mismatch >= sz_branch_ready_min_xy_mismatch);
            const int ready_reason = info_ready ? 1 : (mismatch_ready ? 2 : 0);
            const double ready_hold_target = sz_stage1_baseline_mode ? 0.0 : sz_xeq_ready_hold_time;

            szilard_update_visible_load_state(F_x_now);

            if (left_count > 0 && right_count > 0 && (info_ready || mismatch_ready)) {
                if (sz_xeq_ready_since < 0.0) {
                    sz_xeq_ready_since = simulation_time;
                    sz_branch_ready_reason = ready_reason;
                    if (!cli_quiet) {
                        if (mismatch_ready && !info_ready) {
                            printf("🧪 Szilard: partial x/y decorrelation ready (%d mismatches); holding %.1f sigma-time to define the branch point.\n",
                                   xy_mismatch, sz_xeq_ready_hold_time);
                        } else if (sz_stage1_baseline_mode) {
                            printf("🧪 Szilard: x-side information reached exact zero; saving the branch point now.\n");
                        } else {
                            printf("🧪 Szilard: x-side information low enough (%.3f bit); holding %.1f sigma-time to define the branch point.\n",
                                   I_bits, sz_xeq_ready_hold_time);
                        }
                    }
                }
                if ((simulation_time - sz_xeq_ready_since) >= ready_hold_target) {
                    sz_branch_ready_reason = ready_reason;
                    sz_xflip_enabled = 0;
                    szilard_save_branch_snapshot();
                    if (!cli_quiet) {
                        if (sz_stage1_baseline_mode) {
                            printf("🧠 Szilard: post-measurement branch point ready. Press H for correct y-gated recombination, J for sharp wrong y-gated recombination, n to restore later.\n");
                        } else {
                            printf("🧠 Szilard: post-measurement branch point ready. Press H for correct y-gated recombination, J for sharp wrong y-gated recombination, O for partially observable remixing, n to restore later.\n");
                        }
                    }
                    if (sz_pending_recombine_mode != 0) {
                        const int pending = sz_pending_recombine_mode;
                        sz_pending_recombine_mode = 0;
                        if (!cli_quiet) {
                            if (pending == 1) printf("↪ Szilard: launching queued H branch\n");
                            else if (pending == 2) printf("↪ Szilard: launching queued J branch\n");
                            else if (pending == 3 && !sz_stage1_baseline_mode) printf("↪ Szilard: launching queued O branch\n");
                        }
                        if (pending == 1) szilard_begin_recombine();
                        else if (pending == 2) szilard_begin_recombine_wrong();
                        else if (pending == 3 && !sz_stage1_baseline_mode) szilard_begin_recombine_partial();
                    }
                }
            } else {
                sz_xeq_ready_since = -1.0;
                sz_branch_ready_reason = 0;
            }
        }
    } else if (sz_interactive_phase == 3 || sz_interactive_phase == 4 || sz_interactive_phase == 5) { // recombining
        float piston_stop_px = sz_interactive_piston_home_px;
        float wall_stop_px = sz_interactive_home_px;
        if (sz_interactive_phase == 5) {
            szilard_update_partial_uncertainty_peak();
            double keep_frac = sz_partial_obs_peak_uncertainty_frac;
            if (keep_frac < 0.0) keep_frac = 0.0;
            if (keep_frac > 0.95) keep_frac = 0.95;
            piston_stop_px = sz_interactive_piston_home_px +
                             (float)(keep_frac * (double)(sz_interactive_piston_target_px - sz_interactive_piston_home_px));
            wall_stop_px = sz_interactive_home_px +
                           (float)(keep_frac * (double)(sz_interactive_target_px - sz_interactive_home_px));
        }
        bool piston_done = false;
        bool wall_done = false;
        if (piston_right_x <= piston_stop_px + eps) {
            piston_right_x = piston_stop_px;
            vx_piston_right = 0.0f;
            piston_done = true;
        }
        if (wall_x <= wall_stop_px + eps) {
            wall_x = wall_stop_px;
            vx_wall = 0.0f;
            wall_done = true;
        }
        if (piston_done && wall_done) {
            const int was_wrong = (sz_interactive_phase == 4);
            const int was_partial = (sz_interactive_phase == 5);
            wall_enabled = 0;
            wall_hold_enabled = true;
            wall_is_released = false;
            wall_hold_steps = INT_MAX / 2;
            sz_perm_active = 0;
            sz_perm_dual_active = 0;
            sz_perm_use_memory = 0;
            sz_fixed_gate_active = 0;
            sz_fixed_gate_use_memory = 0;
            sz_xflip_enabled = 0;
            sz_xeq_ready_since = -1.0;
            sz_partial_obs_enabled = 0;
            sz_interactive_phase = 6;
            sz_log_reset_requested = 1;
            szilard_sync_edmd_from_globals();
            if (!cli_quiet) {
                if (was_wrong) printf("✅ Szilard: wrong recombination complete; press 'e' to erase memory\n");
                else if (was_partial) printf("✅ Szilard: partially observable recombination complete; press 'e' to erase memory\n");
                else           printf("✅ Szilard: recombination complete; press 'e' to erase memory\n");
            }
        }
    }
}

// Distribution logging (x/y/speed histograms) for energy_transfer runs
static void log_distributions(double t_abs, double t_rel, int released) {
    if (!dist_log || dist_bins_runtime <= 0 || !dist_hist_x || !dist_hist_y || !dist_hist_v) return;
    if (dist_box_w_sigma <= 0.0 || dist_box_h_sigma <= 0.0 || dist_vmax_runtime <= 0.0) return;

    memset(dist_hist_x, 0, (size_t)dist_bins_runtime * sizeof(int));
    memset(dist_hist_y, 0, (size_t)dist_bins_runtime * sizeof(int));
    memset(dist_hist_v, 0, (size_t)dist_bins_runtime * sizeof(int));

    const double inv_dx = (double)dist_bins_runtime / dist_box_w_sigma;
    const double inv_dy = (double)dist_bins_runtime / dist_box_h_sigma;
    const double inv_dv = (double)dist_bins_runtime / dist_vmax_runtime;

    for (int i = 0; i < particles_active; ++i) {
        double x_sigma = ((double)X[i] - (double)XW1) / (double)PIXELS_PER_SIGMA;
        double y_sigma = ((double)Y[i] - (double)YW1) / (double)PIXELS_PER_SIGMA;
        if (x_sigma < 0.0) x_sigma = 0.0;
        if (y_sigma < 0.0) y_sigma = 0.0;
        if (x_sigma > dist_box_w_sigma) x_sigma = dist_box_w_sigma;
        if (y_sigma > dist_box_h_sigma) y_sigma = dist_box_h_sigma;

        int bx = (int)floor(x_sigma * inv_dx);
        int by = (int)floor(y_sigma * inv_dy);
        if (bx < 0) bx = 0; if (bx >= dist_bins_runtime) bx = dist_bins_runtime - 1;
        if (by < 0) by = 0; if (by >= dist_bins_runtime) by = dist_bins_runtime - 1;
        dist_hist_x[bx] += 1;
        dist_hist_y[by] += 1;

        const double v = sqrt((double)Vx[i] * (double)Vx[i] + (double)Vy[i] * (double)Vy[i]);
        double v_clamp = v;
        if (v_clamp < 0.0) v_clamp = 0.0;
        if (v_clamp > dist_vmax_runtime) v_clamp = dist_vmax_runtime;
        int bv = (int)floor(v_clamp * inv_dv);
        if (bv < 0) bv = 0; if (bv >= dist_bins_runtime) bv = dist_bins_runtime - 1;
        dist_hist_v[bv] += 1;
    }

    fprintf(dist_log,
            "%.6f, %.6f, %d, %d, %.6f, %.6f, %.6f, %.6f, %.6f, \"",
            t_abs, t_rel, released, dist_bins_runtime,
            0.0, dist_box_w_sigma, 0.0, dist_box_h_sigma, dist_vmax_runtime);
    for (int i = 0; i < dist_bins_runtime; ++i) {
        fprintf(dist_log, "%s%d", i ? ";" : "", dist_hist_x[i]);
    }
    fprintf(dist_log, "\", \"");
    for (int i = 0; i < dist_bins_runtime; ++i) {
        fprintf(dist_log, "%s%d", i ? ";" : "", dist_hist_y[i]);
    }
    fprintf(dist_log, "\", \"");
    for (int i = 0; i < dist_bins_runtime; ++i) {
        fprintf(dist_log, "%s%d", i ? ";" : "", dist_hist_v[i]);
    }
    fprintf(dist_log, "\"\n");
    fflush(dist_log);
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

            case SDLK_w:  // ##CHRIS: Toggle wall with 'w' or 'W' (syncs with EDMD)
                wall_enabled = !wall_enabled;

                // ##CHRIS: Sync wall state with EDMD's divider parameter
                if (g_edmd) {
                    if (wall_enabled) {
                        const EDMD_Params* ep = edmd_params(g_edmd);
                        double xs[EDMD_MAX_DIVIDERS];
                        double th[EDMD_MAX_DIVIDERS];
                        int dcount = num_internal_walls;
                        if (dcount > EDMD_MAX_DIVIDERS) dcount = EDMD_MAX_DIVIDERS;
                        for (int d = 0; d < dcount; ++d) {
                            xs[d] = ep->divider_x[d];
                            th[d] = (double)WALL_THICKNESS;
                        }
                        edmd_config_dividers(g_edmd, dcount, xs, th);
                    } else {
                        edmd_config_dividers(g_edmd, 0, NULL, NULL);
                    }

                    // Rebuild collision schedule so EDMD knows about the change
                    edmd_reschedule_all(g_edmd);

                    printf("Wall %s (EDMD divider synced)\n", wall_enabled ? "Enabled" : "Disabled");
                } else {
                    printf("Wall %s\n", wall_enabled ? "Enabled" : "Disabled");
                }
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
                printf("🔔 Wall manually released by key 'r'.\n");
            }
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
            case SDLK_g: { // Szilard: start separation
                if (cli_experiment_preset == EXPERIMENT_PRESET_SZILARD_ENGINE) {
                    szilard_begin_separation();
                }
            } break;
            case SDLK_h: { // Szilard: recombine
                if (cli_experiment_preset == EXPERIMENT_PRESET_SZILARD_ENGINE) {
                    szilard_begin_recombine();
                }
            } break;
            case SDLK_j: { // Szilard: recombine with wrong moving gate (A/B test)
                if (cli_experiment_preset == EXPERIMENT_PRESET_SZILARD_ENGINE) {
                    szilard_begin_recombine_wrong();
                }
            } break;
            case SDLK_o: { // Szilard: partially observable memory-gated recombine
                if (cli_experiment_preset == EXPERIMENT_PRESET_SZILARD_ENGINE) {
                    szilard_begin_recombine_partial();
                }
            } break;
            case SDLK_n: { // Szilard: restore saved post-G branch snapshot
                if (cli_experiment_preset == EXPERIMENT_PRESET_SZILARD_ENGINE && sz_interactive_enabled) {
                    if (g_sz_branch_snapshot.valid) {
                        sz_branch_restore_requested = 1;
                        if (!cli_quiet) {
                            printf("↩️  Szilard: restore to saved post-Xeq branch requested (n)\n");
                        }
                    } else if (!cli_quiet) {
                        printf("Szilard: no saved post-Xeq branch yet. Measure, run G, and wait for Xeq first.\n");
                    }
                }
            } break;

            // ##CHRIS: Heat bath toggle with 'b' key
            case SDLK_b: {
                heatbath_enabled = !heatbath_enabled;
                if (heatbath_enabled) {
                    // Capture current gas temperature as heat bath target
                    heatbath_temperature = compute_measured_temperature_from_ke();
                    // If using EDMD mode, update EDMD parameters
                    if (g_edmd) {
                        ((EDMD_Params*)edmd_params(g_edmd))->heatbath_enabled = 1;
                        ((EDMD_Params*)edmd_params(g_edmd))->heatbath_temperature = (double)heatbath_temperature;
                    }
                    printf("🔥 Heat bath ENABLED at T = %.4f\n", heatbath_temperature);
                } else {
                    // Disable heat bath
                    if (g_edmd) {
                        ((EDMD_Params*)edmd_params(g_edmd))->heatbath_enabled = 0;
                    }
                    printf("❄️  Heat bath DISABLED\n");
                }
            } break;

            // ##CHRIS: Andersen thermostat toggle with 'k' key (NOT 'a' - that's for piston!)
            case SDLK_k: {
                andersen_enabled = !andersen_enabled;
                if (andersen_enabled) {
                    printf("🔵 Andersen thermostat ENABLED (freq=%.3f)\n", andersen_collision_freq);
                } else {
                    printf("⚪ Andersen thermostat DISABLED\n");
                }
            } break;

            // ##CHRIS: Andersen frequency controls with ',' and '.'
            // In interactive Szilard mode, reuse ',' as a slower gate/protocol shortcut because it is
            // easier to hit on non-US keyboards. Outside Szilard, keep the original Andersen control.
            case SDLK_COMMA: {
                if (cli_experiment_preset == EXPERIMENT_PRESET_SZILARD_ENGINE && sz_interactive_enabled) {
                    cli_szilard_sep_duration *= 1.5;
                    if (cli_szilard_sep_duration > 10000.0) cli_szilard_sep_duration = 10000.0;
                    szilard_refresh_interactive_speed();
                    szilard_apply_interactive_velocity_current_phase();
                    if (!cli_quiet) {
                        printf("Szilard sep-duration: %.3f sigma-time  vx=%.4f px/step (; faster, , slower)\n",
                               cli_szilard_sep_duration, (double)sz_interactive_vx_px);
                    }
                } else {
                    if (andersen_collision_freq > 0.001) {
                        andersen_collision_freq -= 0.01;
                        if (andersen_collision_freq < 0.001) andersen_collision_freq = 0.001;
                        if (!cli_quiet) {
                            printf("Andersen freq: %.3f\n", andersen_collision_freq);
                        }
                    }
                }
            } break;

            case SDLK_PERIOD: {  // '.' key - Andersen frequency, or Szilard faster shortcut fallback
                if (cli_experiment_preset == EXPERIMENT_PRESET_SZILARD_ENGINE && sz_interactive_enabled) {
                    cli_szilard_sep_duration /= 1.5;
                    if (cli_szilard_sep_duration < 0.25) cli_szilard_sep_duration = 0.25;
                    szilard_refresh_interactive_speed();
                    szilard_apply_interactive_velocity_current_phase();
                    if (!cli_quiet) {
                        printf("Szilard sep-duration: %.3f sigma-time  vx=%.4f px/step (;/. faster, , slower)\n",
                               cli_szilard_sep_duration, (double)sz_interactive_vx_px);
                    }
                } else if (andersen_collision_freq < 2.0) {
                    andersen_collision_freq += 0.01;
                    if (andersen_collision_freq > 2.0) andersen_collision_freq = 2.0;
                    if (!cli_quiet) {
                        printf("Andersen freq: %.3f\n", andersen_collision_freq);
                    }
                }
            } break;

            // Time scale controls: '-' halves, '=' or '+' doubles, '0' resets
            case SDLK_MINUS:
                time_scale_runtime = fmaxf(0.1f, time_scale_runtime * 0.5f);
                if (!cli_quiet) {
                    printf("Time scale set to %.2fx\n", (double)time_scale_runtime);
                }
                break;
            case SDLK_EQUALS:
            case SDLK_PLUS:  // Also handle '+' key (Shift + '=')
                time_scale_runtime = fminf(10000.0f, time_scale_runtime * 2.0f);  // ##CHRIS: Increased max from 100× to 10000×
                if (!cli_quiet) {
                    printf("Time scale set to %.2fx\n", (double)time_scale_runtime);
                }
                break;
            case SDLK_0:
                time_scale_runtime = 1.0f;
                if (!cli_quiet) {
                    printf("Time scale reset to 1.00x\n");
                }
                break;

		            case SDLK_v: { // GUI: toggle distribution panels
		                gui_dist_mode = (gui_dist_mode + 1) % 4; // 0->1->2->3->0
		                if (!cli_quiet) {
		                    const char *name = (gui_dist_mode == 0) ? "off"
		                                     : (gui_dist_mode == 1) ? "velocity (global)"
		                                     : (gui_dist_mode == 2) ? "velocity (per-segment)"
		                                     : "density+velocity (per-segment)";
				                    printf("📊 Distributions: %s (v)\n", name);
				                    if (gui_dist_mode == 2 || gui_dist_mode == 3) {
				                        printf("    legend: vx=grey, vy=yellow; mode 3 stacks: x-density(global), [x-density(per-seg) with x], speed |v| (green), vx/vy; density overlays: baseline(u), curve(y); exact(z); markers(l)\n");
				                    }
				                }
				            } break;

                    case SDLK_z: { // GUI: toggle EXACT histogram mode (no smoothing/fit overlays)
                        gui_hist_exact = !gui_hist_exact;
                        if (gui_hist_exact) {
                            // In exact mode, suppress overlays that create "extra bins"/lines.
                            gui_vel_fit_overlay = 0;
                            gui_density_uniform_overlay = 0;
                            gui_density_curve_overlay = 0;
                        }
                        if (!cli_quiet) {
                            printf("📊 Hist exact: %s (z)\n", gui_hist_exact ? "ON (raw bins only)" : "off");
                        }
                    } break;

                    case SDLK_l: { // GUI: toggle marker lines in distribution panels
                        gui_dist_markers = !gui_dist_markers;
                        if (!cli_quiet) {
                            printf("📍 Dist markers: %s (l)\n", gui_dist_markers ? "ON" : "off");
                        }
                    } break;

		            case SDLK_f: { // GUI: toggle fixed velocity axis range (vmax)
		                gui_vel_range_freeze = !gui_vel_range_freeze;
		                if (gui_vel_range_freeze) {
		                    gui_vel_range_vmax_fixed = gui_compute_vmax_fixed();
		                }
		                if (!cli_quiet) {
		                    if (gui_vel_range_freeze) {
		                        printf("📏 Vel range: FIXED vmax=%.3f (f)\n", (double)gui_vel_range_vmax_fixed);
		                    } else {
		                        printf("📏 Vel range: AUTO (f)\n");
		                    }
		                }
		            } break;

                    case SDLK_c: { // GUI: toggle fitted distribution overlay on velocity panels
                        gui_vel_fit_overlay = !gui_vel_fit_overlay;
                        if (!cli_quiet) {
                            printf("📈 Vel fit overlay: %s (c)\n", gui_vel_fit_overlay ? "ON" : "off");
                        }
                    } break;

                    case SDLK_u: { // GUI: toggle uniform baseline overlay on density panels
                        gui_density_uniform_overlay = !gui_density_uniform_overlay;
                        if (!cli_quiet) {
                            printf("📏 Density baseline: %s (u)\n", gui_density_uniform_overlay ? "ON" : "off");
                        }
                    } break;

                    case SDLK_y: { // Szilard: cycle particle view, otherwise GUI density-curve toggle
                        if (cli_experiment_preset == EXPERIMENT_PRESET_SZILARD_ENGINE && sz_interactive_enabled) {
                            if (sz_view_mode == SZ_VIEW_OVERLAY) {
                                sz_view_mode = SZ_VIEW_SPECIES;
                            } else if (sz_view_mode == SZ_VIEW_SPECIES) {
                                sz_view_mode = SZ_VIEW_MEMORY;
                            } else {
                                sz_view_mode = SZ_VIEW_SPECIES;
                            }
                            if (!cli_quiet) {
                                printf("Szilard particle view: %s (y toggles species<->memory; overlay via CLI)\n",
                                       szilard_particle_view_mode_label((int)sz_view_mode));
                            }
                        } else {
                            gui_density_curve_overlay = !gui_density_curve_overlay;
                            if (!cli_quiet) {
                                printf("📈 Global x-density curve: %s (y)\n", gui_density_curve_overlay ? "ON" : "off");
                            }
                        }
                    } break;

                    case SDLK_x: { // GUI: toggle per-compartment x-density panel (stacked mode)
                        gui_show_x_density_by_segment = !gui_show_x_density_by_segment;
                        if (!cli_quiet) {
                            printf("📊 x-density per-compartment: %s (x)\n", gui_show_x_density_by_segment ? "ON" : "off");
                        }
                    } break;

                    case SDLK_SEMICOLON: { // Szilard: faster protocol (shorter duration)
                        if (cli_experiment_preset == EXPERIMENT_PRESET_SZILARD_ENGINE && sz_interactive_enabled) {
                            cli_szilard_sep_duration /= 1.5;
                            if (cli_szilard_sep_duration < 0.25) cli_szilard_sep_duration = 0.25;
                            szilard_refresh_interactive_speed();
                            szilard_apply_interactive_velocity_current_phase();
                            if (!cli_quiet) {
                                printf("Szilard sep-duration: %.3f sigma-time  vx=%.4f px/step (; faster, , slower)\n",
                                       cli_szilard_sep_duration, (double)sz_interactive_vx_px);
                            }
                        }
                    } break;

                    case SDLK_QUOTE: { // Szilard: slower protocol (longer duration)
                        if (cli_experiment_preset == EXPERIMENT_PRESET_SZILARD_ENGINE && sz_interactive_enabled) {
                            cli_szilard_sep_duration *= 1.5;
                            if (cli_szilard_sep_duration > 10000.0) cli_szilard_sep_duration = 10000.0;
                            szilard_refresh_interactive_speed();
                            szilard_apply_interactive_velocity_current_phase();
                            if (!cli_quiet) {
                                printf("Szilard sep-duration: %.3f sigma-time  vx=%.4f px/step (; faster, , slower)\n",
                                       cli_szilard_sep_duration, (double)sz_interactive_vx_px);
                            }
                        }
                    } break;

                    // GUI window height controls (pure UI): PgUp/PgDn changes the bottom plot area size.
                    case SDLK_PAGEUP: {
                        gui_adjust_window_height(+gui_window_height_step);
                    } break;
                    case SDLK_PAGEDOWN: {
                        gui_adjust_window_height(-gui_window_height_step);
                    } break;
	        }
	    }
	}


/// Constants for acceleration INIT piston
const float acceleration = 25;  // Increased acceleration
const float max_velocity = 1000000.0;   // Increased max speed
const float deceleration = 25;

static void configure_right_piston_protocol_from_cli(void) {
    // Protocol parameters: speed is in σ/σ-time (numerically same as px/pixel-time).
    // Prefer right-piston-specific overrides when provided.
    if (cli_piston_right_step_speed_set && cli_piston_right_step_speed > 0.0f) {
        piston_protocol_max_speed = cli_piston_right_step_speed;
    } else if (cli_piston_right_vmax_set && cli_piston_right_vmax > 0.0f) {
        piston_protocol_max_speed = cli_piston_right_vmax;
    } else {
        piston_protocol_max_speed = (cli_piston_speed_set && cli_piston_speed > 0.0f) ? cli_piston_speed : 1.0f;
    }

    if (cli_piston_right_duration_set && cli_piston_right_duration > 0.0f) {
        piston_protocol_duration = cli_piston_right_duration;
    } else {
        piston_protocol_duration = (cli_piston_duration_set && cli_piston_duration > 0.0f) ? cli_piston_duration : 2.0f;
    }

    piston_protocol_v0 = (cli_piston_right_v0_set && cli_piston_right_v0 >= 0.0f) ? cli_piston_right_v0 : 0.0f;
    piston_protocol_linear_use_gradient = (cli_piston_right_gradient_set);
    piston_protocol_linear_gradient = cli_piston_right_gradient_set ? cli_piston_right_gradient : 0.0f;
    piston_protocol_sigmoid_steepness = (cli_piston_right_sigmoid_steepness_set && cli_piston_right_sigmoid_steepness > 0.0f)
        ? cli_piston_right_sigmoid_steepness
        : 0.0f;
    piston_protocol_start_time = -1.0f;

    // Initialize piston speed consistent with protocol at t=0.
    switch (cli_protocol) {
        case PROTOCOL_STEP:
            piston_step_speed = piston_step_direction * fabsf(piston_protocol_max_speed);
            break;
        case PROTOCOL_LINEAR:
            piston_step_speed = piston_protocol_linear_use_gradient ? (piston_step_direction * fmaxf(0.0f, piston_protocol_v0)) : 0.0f;
            break;
        case PROTOCOL_SIGMOIDAL:
            piston_step_speed = piston_step_direction * fmaxf(0.0f, piston_protocol_v0);
            break;
        default:
            piston_step_speed = 0.0f;
            break;
    }
}

static void start_right_piston_step_to_target_px(float target_px) {
    float min_target = piston_left_x + 25.0f;
    float max_target = (float)XW2;
    float target = target_px;
    if (target < min_target) target = min_target;
    if (target > max_target) target = max_target;
    piston_step_target = target;

    if (fabsf(piston_step_target - piston_right_x) <= 1e-6f) {
        piston_step_active = false;
        vx_piston_right = 0.0f;
        piston_protocol_start_time = -1.0f;
        return;
    }

    piston_step_direction = (piston_step_target >= piston_right_x) ? 1.0f : -1.0f;
    configure_right_piston_protocol_from_cli();
    vx_piston_right = piston_step_speed;
    piston_step_active = true;
    simulation_started = 1;

    if (cli_enable_energy_measurement || cli_experiment_preset == EXPERIMENT_PRESET_ENERGY_TRANSFER) {
        initialize_energy_measurement();
    }
    start_piston_push_tracking();

    if (!cli_quiet) {
        printf("➡️  Piston step triggered: target %.2f px, direction=%s, protocol: %s\n",
               piston_step_target,
               piston_step_direction > 0.0f ? "expand" : "compress",
               cli_protocol == PROTOCOL_SIGMOIDAL ? "SIGMOIDAL" :
               cli_protocol == PROTOCOL_LINEAR ? "LINEAR" :
               cli_protocol == PROTOCOL_SINUSOIDAL ? "SINUSOIDAL" : "STEP");
    }
}

static void start_right_piston_step_to_target_sigma(float target_sigma) {
    start_right_piston_step_to_target_px((float)XW1 + (float)(target_sigma * (double)PIXELS_PER_SIGMA));
}

static void trigger_piston_step_with_protocol(void) {
    // Default: move right piston inward by 20% of box width.
    // Override with --piston-travel-sigma for experiments.
    float travel_sigma = 0.20f * ((float)SIM_WIDTH / (float)PIXELS_PER_SIGMA);
    if (cli_piston_right_travel_set && cli_piston_right_travel_sigma > 0.0f) {
        travel_sigma = cli_piston_right_travel_sigma;
    } else if (cli_piston_travel_set && cli_piston_travel_sigma > 0.0f) {
        travel_sigma = cli_piston_travel_sigma;
    }
    {
        float travel_px = travel_sigma * (float)PIXELS_PER_SIGMA;
        float target = (float)XW2 - travel_px;
        start_right_piston_step_to_target_px(target);
    }
}


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
            if (cli_experiment_preset == EXPERIMENT_PRESET_SZILARD_ENGINE) {
                szilard_begin_separation();
            } else {
                trigger_piston_step_with_protocol();
            }
        } break;
            
            case SDLK_e:
                if (cli_experiment_preset == EXPERIMENT_PRESET_SZILARD_ENGINE) {
                    szilard_manual_erase_memory();
                } else {
                    print_energy_measurement();
                }
                break;

            case SDLK_m: {
                if (cli_experiment_preset == EXPERIMENT_PRESET_SZILARD_ENGINE) {
                    szilard_measure_memory_from_species();
                } else if (pl_enabled) {
                    pl_tensor_mode = (TensorMode)((pl_tensor_mode + 1) % 4);
                    const char* mode_names[] = {"2D Matrix", "3D Distance", "4D Velocity", "5D Angular"};
                    printf("📊 Tensor mode switched to: %s\n", mode_names[pl_tensor_mode]);
                } else {
                    printf("⚠️  Particle-life not active (use --experiment=particlelife)\n");
                }
                break;
            }
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
        if (gui_mb_speed_view) {
            set_draw_color_speed_heat(((float)i + 0.5f) / (float)NUM_BINS, 220);
        } else {
            SDL_SetRenderDrawColor(_renderer, 255, 255, 0, 255);
        }
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

        if (gui_mb_speed_view) SDL_SetRenderDrawColor(_renderer, 230, 235, 230, 245);
        else SDL_SetRenderDrawColor(_renderer, 255, 255, 0, 255);
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

// Per-segment velocity component histograms (vx + vy) drawn in the bottom panel, aligned to compartments.
// Quick equilibrium check: vx and vy should both look ~Gaussian (same width) within each segment.
static inline float gui_smooth7(const float *a, int i, int n) {
    // Small symmetric kernel for visually smooth histograms at low N.
    // Weights sum to 1.0.
    const float w0 = 0.05f, w1 = 0.09f, w2 = 0.12f, w3 = 0.48f;
    int i0 = (i - 3 < 0) ? 0 : i - 3;
    int i1 = (i - 2 < 0) ? 0 : i - 2;
    int i2 = (i - 1 < 0) ? 0 : i - 1;
    int i3 = i;
    int i4 = (i + 1 >= n) ? (n - 1) : i + 1;
    int i5 = (i + 2 >= n) ? (n - 1) : i + 2;
    int i6 = (i + 3 >= n) ? (n - 1) : i + 3;
    return w0 * a[i0] + w1 * a[i1] + w2 * a[i2] + w3 * a[i3] + w2 * a[i4] + w1 * a[i5] + w0 * a[i6];
}

static inline int gui_panel_inner(int y0, int h, int *out_y, int *out_h) {
    // Reserve some space at the top for labels so bars don't overlap text and look less "cut off".
    const int top_pad = 28;
    const int bottom_pad = 4;
    int yy = y0 + top_pad;
    int hh = h - top_pad - bottom_pad;
    if (hh < 16) return 0;
    if (out_y) *out_y = yy;
    if (out_h) *out_h = hh;
    return 1;
}

static inline float gui_gauss_pdf(float v, float mu, float sigma) {
    if (!(sigma > 0.0f)) return 0.0f;
    const float z = (v - mu) / sigma;
    const float inv = 1.0f / (sigma * sqrtf(2.0f * (float)M_PI));
    return inv * expf(-0.5f * z * z);
}

// 2D speed distribution when vx,vy ~ N(0, sigma^2): Rayleigh(sigma)
static inline float gui_rayleigh_pdf(float v, float sigma) {
    if (!(sigma > 0.0f) || v < 0.0f) return 0.0f;
    const float s2 = sigma * sigma;
    return (v / s2) * expf(-(v * v) / (2.0f * s2));
}

static float gui_compute_vmax_fixed(void) {
    // Choose a fixed velocity axis range that avoids clipping current tails.
    // Units match Vx/Vy (simulation units).
    const int n = particles_active > 0 ? particles_active : 0;
    if (n <= 0) return 1.0f;
    double m2 = 0.0;
    double max_abs = 0.0;
    double max_sp = 0.0;
    for (int i = 0; i < n; ++i) {
        double vx = (double)Vx[i];
        double vy = (double)Vy[i];
        m2 += vx * vx + vy * vy;
        double ax = fabs(vx);
        double ay = fabs(vy);
        if (ax > max_abs) max_abs = ax;
        if (ay > max_abs) max_abs = ay;
        double sp = sqrt(vx * vx + vy * vy);
        if (sp > max_sp) max_sp = sp;
    }
    double rms = sqrt(m2 / (double)n);
    // Component sigma ≈ rms/sqrt(2). Use ~6σ for a visible tail.
    double v_from_rms = 6.0 * (rms / sqrt(2.0));
    double v = fmax(0.5, v_from_rms);
    v = fmax(v, 1.10 * max_abs);
    v = fmax(v, 1.10 * max_sp);
    if (v > (double)MAX_VELOCITY) v = (double)MAX_VELOCITY;
    return (float)v;
}

static void gui_mb_capture_speed_scale_now(void) {
    float v = gui_compute_vmax_fixed();
    if (!(v > 1e-6f) || !isfinite(v)) v = 1.0f;
    gui_mb_speed_vmax_ref = v;
}

static float gui_mb_speed_vmax_for_render(void) {
    if (gui_mb_speed_vmax_ref > 1e-6f && isfinite(gui_mb_speed_vmax_ref)) {
        return gui_mb_speed_vmax_ref;
    }
    if (cli_dist_vmax > 1e-6f && isfinite(cli_dist_vmax)) {
        return cli_dist_vmax;
    }
    if (dist_vmax_runtime > 1e-6 && isfinite(dist_vmax_runtime)) {
        return (float)dist_vmax_runtime;
    }

    // Deliberately do not auto-sample current particle speeds here.
    // Without a manual SET, use a fixed thermal reference scale instead of
    // MAX_VELOCITY; otherwise normal thermal particles all look dark blue.
    const double Tref = heatbath_enabled ? (double)heatbath_temperature : (double)temperature_runtime;
    const double kbt = (double)kB_effective() * fmax(1e-12, Tref);
    double v = 6.0 * sqrt(fmax(1e-12, 2.0 * kbt / (double)PARTICLE_MASS));
    if (!(v > 1e-6) || !isfinite(v)) v = 1.0;
    if (v > (double)MAX_VELOCITY) v = (double)MAX_VELOCITY;
    return (float)v;
}

// GUI helper: pick an "anchor" x-position (in pixels, global coords) to mark the effective left boundary.
// Requested behavior: choose whichever is furthest right among:
// - left piston plane (piston_left_x + 5)
// - spring wall / leftmost divider (approx: leftmost divider + thickness/2)
static float gui_compute_left_anchor_x_display(float pistonL_x_disp) {
    float anchor = pistonL_x_disp + 5.0f;
    if (num_internal_walls > 0) {
        float tmp_walls[MAX_DIVIDER_CAPACITY];
        get_sorted_wall_positions_display(tmp_walls);
        float leftmost_right_face = tmp_walls[0] + WALL_THICKNESS * 0.5f;
        if (leftmost_right_face > anchor) anchor = leftmost_right_face;
    }
    if (anchor < (float)XW1) anchor = (float)XW1;
    if (anchor > (float)XW2) anchor = (float)XW2;
    return anchor;
}

static inline int gui_clamp_int(int v, int lo, int hi) {
    if (v < lo) return lo;
    if (v > hi) return hi;
    return v;
}

// Choose x-density histogram bins so bin width roughly matches the particle diameter.
// This reduces jitter/aliasing when walls move and avoids bins so fine they are meaningless.
static inline int gui_bins_for_xdensity(float width_px, int max_bins) {
    float bin_px = 2.0f * (float)PARTICLE_RADIUS;
    if (!isfinite(bin_px) || bin_px < 2.0f) bin_px = 2.0f;
    int bins = (int)floorf(width_px / bin_px);
    bins = gui_clamp_int(bins, 16, max_bins);
    // Ensure we can draw at least ~1px per bin in tight compartments.
    if (width_px > 0.0f) {
        int cap = (int)floorf(width_px);
        if (cap >= 8 && bins > cap) bins = cap;
    }
    if (bins < 8) bins = 8;
    return bins;
}

// Choose velocity histogram bins based on sample count (so low-N doesn't look like sparse spikes).
static inline int gui_bins_for_velocity(int n_samples, float seg_w_px, int max_bins) {
    int n = (n_samples > 0) ? n_samples : 1;
    int bins = (int)floorf(6.0f * sqrtf((float)n)); // ~84 bins for n=200
    bins = gui_clamp_int(bins, 24, max_bins);
    // Rendering cap: avoid many bins mapping to the same x pixel in narrow segments.
    if (seg_w_px > 0.0f) {
        int cap = (int)floorf(seg_w_px);
        if (cap >= 8 && bins > cap) bins = cap;
    }
    if (bins < 8) bins = 8;
    return bins;
}

static void render_velocity_histograms_by_segment_panel(int y0, int h) {
    if (!segment_counts || segment_count <= 0 || segment_count > MAX_DIVIDER_CAPACITY + 1) return;
    if (particles_active <= 0) return;
    if (h < 20) return;

    SDL_SetRenderDrawBlendMode(_renderer, SDL_BLENDMODE_BLEND);
    sync_all_wall_positions();

    float left_bounds[MAX_DIVIDER_CAPACITY + 1];
    float right_bounds[MAX_DIVIDER_CAPACITY + 1];
    compute_segment_bounds_display(left_bounds, right_bounds);

    enum { GUI_HIST_MAX_BINS = 256, GUI_HIST_MAX_SEGS = MAX_DIVIDER_CAPACITY + 1 };
    const float alpha = 0.25f; // EMA smoothing for flicker reduction (time), bin smoothing below handles the rest
    const float vmax_hard = (float)MAX_VELOCITY;

    // Draw a faint background for the panel (clip to current piston window)
    const float pistonL_x_disp = piston_left_x + vx_piston_left * gui_render_dt_rem;
    const float pistonR_x_disp = piston_right_x + vx_piston_right * gui_render_dt_rem;
    const int panel_left = (int)fmaxf((float)XW1, pistonL_x_disp + 5.0f);
    const int panel_right = (int)fminf((float)XW2, pistonR_x_disp);
    const int panel_w = panel_right - panel_left;
    if (panel_w < 20) return;

    SDL_SetRenderDrawColor(_renderer, 20, 20, 20, 220);
    SDL_Rect bg = { panel_left, y0, panel_w, h };
    SDL_RenderFillRect(_renderer, &bg);

    // Draw vertical lines at wall boundaries
    if (gui_dist_markers) {
        SDL_SetRenderDrawColor(_renderer, 255, 255, 255, 40);
        if (num_internal_walls > 0) {
            float tmp_walls[MAX_DIVIDER_CAPACITY];
            get_sorted_wall_positions_display(tmp_walls);
            for (int w = 0; w < num_internal_walls; ++w) {
                int xw = (int)tmp_walls[w];
                if (xw < panel_left || xw > panel_right) continue;
                SDL_RenderDrawLine(_renderer, xw, y0, xw, y0 + h);
            }
        }
        // Piston boundary line (right)
        SDL_SetRenderDrawColor(_renderer, 0, 255, 0, 80);
        SDL_RenderDrawLine(_renderer, panel_right, y0, panel_right, y0 + h);

        // Requested: green anchor line for the effective left boundary (pistonL or spring/leftmost wall)
        {
            float ax = gui_compute_left_anchor_x_display(pistonL_x_disp);
            int xax = (int)ax;
            if (xax >= panel_left && xax <= panel_right) {
                SDL_SetRenderDrawColor(_renderer, 0, 255, 0, 50);
                SDL_RenderDrawLine(_renderer, xax, y0, xax, y0 + h);
            }
        }
    }

    int inner_y0 = y0, inner_h = h;
    if (!gui_panel_inner(y0, h, &inner_y0, &inner_h)) return;

    // EMA storage + per-segment auto scale.
    static float ema_vx[GUI_HIST_MAX_SEGS][GUI_HIST_MAX_BINS];
    static float ema_vy[GUI_HIST_MAX_SEGS][GUI_HIST_MAX_BINS];
    static float ema_vmax[GUI_HIST_MAX_SEGS];
    static int ema_inited = 0;
    static int ema_segs = 0;
    static int ema_bins_per_seg[GUI_HIST_MAX_SEGS];
    if (!ema_inited || ema_segs != segment_count) {
        memset(ema_vx, 0, sizeof(ema_vx));
        memset(ema_vy, 0, sizeof(ema_vy));
        for (int s = 0; s < GUI_HIST_MAX_SEGS; ++s) ema_vmax[s] = 0.0f;
        for (int s = 0; s < GUI_HIST_MAX_SEGS; ++s) ema_bins_per_seg[s] = 0;
        ema_inited = 1;
        ema_segs = segment_count;
    }

    // For each segment, compute and draw histograms
    for (int seg = 0; seg < segment_count; ++seg) {
        const float xl = left_bounds[seg];
        const float xr = right_bounds[seg];
        const float seg_w = xr - xl;
        if (seg_w < 8.0f) continue;

        int vx_hist[GUI_HIST_MAX_BINS];
        int vy_hist[GUI_HIST_MAX_BINS];

        // Auto-scale display range based on per-segment variance (so it doesn't look like a spike),
        // but avoid clipping fast tails after piston compression.
        double meanx = 0.0, meany = 0.0;
        double m2x = 0.0, m2y = 0.0;
        double max_abs = 0.0;
        int nseg = 0;
        for (int i = 0; i < particles_active; ++i) {
            int s = segment_index_for_position(X[i]);
            if (s != seg) continue;
            double vx = (double)Vx[i];
            double vy = (double)Vy[i];
            double ax = fabs(vx);
            double ay = fabs(vy);
            if (ax > max_abs) max_abs = ax;
            if (ay > max_abs) max_abs = ay;
            nseg++;
            double dx = vx - meanx;
            meanx += dx / (double)nseg;
            m2x += dx * (vx - meanx);
            double dy = vy - meany;
            meany += dy / (double)nseg;
            m2y += dy * (vy - meany);
        }
        float vmax_disp = vmax_hard;
        if (gui_vel_range_freeze && gui_vel_range_vmax_fixed > 0.0f) {
            vmax_disp = fminf(vmax_hard, gui_vel_range_vmax_fixed);
        }
        if (nseg >= 2) {
            double varx = m2x / (double)(nseg - 1);
            double vary = m2y / (double)(nseg - 1);
            double sig = sqrt(fmax(0.0, fmax(varx, vary)));
            double vnew = fmax(0.25, 6.0 * sig);
            vnew = fmax(vnew, 1.10 * max_abs);
            if (vnew > (double)vmax_hard) vnew = (double)vmax_hard;
            float vnewf = (float)vnew;
            if (!(gui_vel_range_freeze && gui_vel_range_vmax_fixed > 0.0f)) {
                if (ema_vmax[seg] <= 0.0f) ema_vmax[seg] = vnewf;
                else ema_vmax[seg] = 0.85f * ema_vmax[seg] + 0.15f * vnewf;
                vmax_disp = ema_vmax[seg];
            }
        } else if (ema_vmax[seg] > 0.0f) {
            if (!(gui_vel_range_freeze && gui_vel_range_vmax_fixed > 0.0f)) {
                vmax_disp = ema_vmax[seg];
            }
        }
        if (vmax_disp < 1e-6f) vmax_disp = 1e-6f;

        const int bins = gui_bins_for_velocity(nseg, seg_w, GUI_HIST_MAX_BINS);
        if (ema_bins_per_seg[seg] != bins) {
            for (int i = 0; i < GUI_HIST_MAX_BINS; ++i) { ema_vx[seg][i] = 0.0f; ema_vy[seg][i] = 0.0f; }
            ema_bins_per_seg[seg] = bins;
        }
        for (int i = 0; i < bins; ++i) { vx_hist[i] = 0; vy_hist[i] = 0; }

        for (int i = 0; i < particles_active; ++i) {
            int s = segment_index_for_position(X[i]);
            if (s != seg) continue;
            float vx = Vx[i];
            float vy = Vy[i];

            int b_vx = (int)(((vx + vmax_disp) / (2.0f * vmax_disp)) * (float)bins);
            int b_vy = (int)(((vy + vmax_disp) / (2.0f * vmax_disp)) * (float)bins);
            if (b_vx < 0) b_vx = 0; if (b_vx >= bins) b_vx = bins - 1;
            if (b_vy < 0) b_vy = 0; if (b_vy >= bins) b_vy = bins - 1;

            vx_hist[b_vx] += 1;
            vy_hist[b_vy] += 1;
        }

        float max_vx = 1.0f, max_vy = 1.0f;
        if (!gui_hist_exact) {
            // EMA smooth + track maxima
            for (int i = 0; i < bins; ++i) {
                float vvx = (float)vx_hist[i];
                float vvy = (float)vy_hist[i];
                if (ema_vx[seg][i] == 0.0f && ema_vy[seg][i] == 0.0f && (vvx > 0.0f || vvy > 0.0f)) {
                    ema_vx[seg][i] = vvx;
                    ema_vy[seg][i] = vvy;
                } else {
                    ema_vx[seg][i] = (1.0f - alpha) * ema_vx[seg][i] + alpha * vvx;
                    ema_vy[seg][i] = (1.0f - alpha) * ema_vy[seg][i] + alpha * vvy;
                }
                // Smooth for display scaling (reduces spiky appearance at low N).
                float sx = gui_smooth7(ema_vx[seg], i, bins);
                float sy = gui_smooth7(ema_vy[seg], i, bins);
                if (sx > max_vx) max_vx = sx;
                if (sy > max_vy) max_vy = sy;
            }
        } else {
            // Exact mode: raw bin counts only.
            for (int i = 0; i < bins; ++i) {
                float sx = (float)vx_hist[i];
                float sy = (float)vy_hist[i];
                if (sx > max_vx) max_vx = sx;
                if (sy > max_vy) max_vy = sy;
            }
        }

        // Draw vx histogram (grey)
        for (int i = 0; i < bins; ++i) {
            float sx = gui_hist_exact ? (float)vx_hist[i] : gui_smooth7(ema_vx[seg], i, bins);
            int bh = (int)((sx * (float)inner_h) / max_vx);
            SDL_SetRenderDrawColor(_renderer, 180, 180, 180, 120);
            int x0 = (int)floorf(xl + ((float)i / (float)bins) * seg_w);
            int x1 = (int)floorf(xl + ((float)(i + 1) / (float)bins) * seg_w);
            int bw = x1 - x0;
            if (bw < 1) bw = 1;
            if (x0 + bw > (int)xr) bw = (int)xr - x0;
            if (bw <= 0) continue;
            SDL_Rect r = { x0, inner_y0 + (inner_h - bh), bw, bh };
            SDL_RenderFillRect(_renderer, &r);
        }

        // Draw vy histogram (yellow) overlaid
        for (int i = 0; i < bins; ++i) {
            float sy = gui_hist_exact ? (float)vy_hist[i] : gui_smooth7(ema_vy[seg], i, bins);
            int bh = (int)((sy * (float)inner_h) / max_vy);
            SDL_SetRenderDrawColor(_renderer, 255, 220, 0, 140);
            int x0 = (int)floorf(xl + ((float)i / (float)bins) * seg_w);
            int x1 = (int)floorf(xl + ((float)(i + 1) / (float)bins) * seg_w);
            int bw = x1 - x0;
            if (bw < 1) bw = 1;
            if (x0 + bw > (int)xr) bw = (int)xr - x0;
            if (bw <= 0) continue;
            SDL_Rect r = { x0, inner_y0 + (inner_h - bh), bw, bh };
            SDL_RenderFillRect(_renderer, &r);
        }

        // Draw a reference line for vx=0 (middle of bins) in this segment
        if (gui_dist_markers) {
            int x_mid = (int)floorf(xl + 0.5f * seg_w);
            SDL_SetRenderDrawColor(_renderer, 255, 255, 255, 60);
            SDL_RenderDrawLine(_renderer, x_mid, inner_y0, x_mid, inner_y0 + inner_h);
        }

        // Draw mean/median markers for vx and vy (low-alpha so they don't dominate).
        // Median is computed from hist counts (no per-frame sorting needed).
        if (gui_dist_markers && nseg >= 2) {
            const float dv = (2.0f * vmax_disp) / (float)bins;

            // vx mean marker
            int x_mu_vx = (int)floorf(xl + (((float)((meanx + vmax_disp) / (2.0f * vmax_disp))) * seg_w));
            if (x_mu_vx >= (int)xl && x_mu_vx <= (int)xr) {
                SDL_SetRenderDrawColor(_renderer, 200, 200, 200, 90);
                SDL_RenderDrawLine(_renderer, x_mu_vx, inner_y0, x_mu_vx, inner_y0 + inner_h);
            }
            // vx median marker
            int cum = 0, med_bin = 0;
            for (int i = 0; i < bins; ++i) { cum += vx_hist[i]; if (cum * 2 >= nseg) { med_bin = i; break; } }
            float vx_med = -vmax_disp + ((float)med_bin + 0.5f) * dv;
            int x_med_vx = (int)floorf(xl + (((vx_med + vmax_disp) / (2.0f * vmax_disp)) * seg_w));
            if (x_med_vx >= (int)xl && x_med_vx <= (int)xr) {
                SDL_SetRenderDrawColor(_renderer, 255, 255, 255, 70);
                for (int yy = inner_y0; yy < inner_y0 + inner_h; yy += 10) {
                    SDL_RenderDrawLine(_renderer, x_med_vx, yy, x_med_vx, yy + 6);
                }
            }

            // vy mean marker
            int x_mu_vy = (int)floorf(xl + (((float)((meany + vmax_disp) / (2.0f * vmax_disp))) * seg_w));
            if (x_mu_vy >= (int)xl && x_mu_vy <= (int)xr) {
                SDL_SetRenderDrawColor(_renderer, 255, 240, 120, 90);
                SDL_RenderDrawLine(_renderer, x_mu_vy, inner_y0, x_mu_vy, inner_y0 + inner_h);
            }
            // vy median marker
            cum = 0; med_bin = 0;
            for (int i = 0; i < bins; ++i) { cum += vy_hist[i]; if (cum * 2 >= nseg) { med_bin = i; break; } }
            float vy_med = -vmax_disp + ((float)med_bin + 0.5f) * dv;
            int x_med_vy = (int)floorf(xl + (((vy_med + vmax_disp) / (2.0f * vmax_disp)) * seg_w));
            if (x_med_vy >= (int)xl && x_med_vy <= (int)xr) {
                SDL_SetRenderDrawColor(_renderer, 255, 220, 0, 70);
                for (int yy = inner_y0; yy < inner_y0 + inner_h; yy += 10) {
                    SDL_RenderDrawLine(_renderer, x_med_vy, yy, x_med_vy, yy + 6);
                }
            }

            // In-panel stats label (requested): show vx/vy distribution parameters inside each compartment region.
            if (font && seg_w >= 100.0f) {
                double sdx = 0.0, sdy = 0.0;
                if (nseg >= 2) {
                    sdx = sqrt(fmax(0.0, m2x / (double)(nseg - 1)));
                    sdy = sqrt(fmax(0.0, m2y / (double)(nseg - 1)));
                }
                SDL_Color ctxt = {230, 230, 230, 255};
                char sline1[160];
                char sline2[160];
                snprintf(sline1, sizeof(sline1), "vx mu=%.2f  med=%.2f  sd=%.2f", (double)meanx, (double)vx_med, sdx);
                snprintf(sline2, sizeof(sline2), "vy mu=%.2f  med=%.2f  sd=%.2f", (double)meany, (double)vy_med, sdy);
                TTF_Font *ff = font_small ? font_small : font;
                const int dy_stats = font_small ? TTF_FontLineSkip(ff) : 14;
                // Place below the panel title line to avoid collisions.
                draw_text_clipped(_renderer, ff, sline1, (int)xl + 4, y0 + 16, ctxt, (int)seg_w - 8);
                draw_text_clipped(_renderer, ff, sline2, (int)xl + 4, y0 + 16 + dy_stats, ctxt, (int)seg_w - 8);
            }
        }

        // Optional: overlay a fitted Gaussian on vx/vy (visual sanity check for MB-like behavior).
        // Uses the same (mu, sigma) estimated above (Welford), and scales to expected counts per bin.
        if (gui_vel_fit_overlay && !gui_hist_exact && nseg >= 10) {
            const float dv = (2.0f * vmax_disp) / (float)bins;
            const float mu_x = (float)meanx;
            const float mu_y = (float)meany;
            float sig_x = 0.0f, sig_y = 0.0f;
            if (nseg >= 2) {
                double varx = m2x / (double)(nseg - 1);
                double vary = m2y / (double)(nseg - 1);
                sig_x = (float)sqrt(fmax(0.0, varx));
                sig_y = (float)sqrt(fmax(0.0, vary));
            }

            // vx curve (light grey)
            if (sig_x > 0.0f) {
                SDL_SetRenderDrawColor(_renderer, 235, 235, 235, 200);
                int prev_x = -1, prev_y = -1;
                for (int i = 0; i < bins; ++i) {
                    float v = -vmax_disp + ((float)i + 0.5f) * dv;
                    float expected = (float)nseg * gui_gauss_pdf(v, mu_x, sig_x) * dv;
                    int bh = (int)((expected * (float)inner_h) / max_vx);
                    if (bh < 0) bh = 0;
                    if (bh > inner_h) bh = inner_h;
                    int cx = (int)floorf(xl + (((float)i + 0.5f) / (float)bins) * seg_w);
                    int cy = inner_y0 + (inner_h - bh);
                    if (prev_x >= 0) {
                        // dashed look: draw every other segment
                        if ((i % 2) == 0) SDL_RenderDrawLine(_renderer, prev_x, prev_y, cx, cy);
                    }
                    prev_x = cx; prev_y = cy;
                }
            }

            // vy curve (brighter yellow)
            if (sig_y > 0.0f) {
                SDL_SetRenderDrawColor(_renderer, 255, 245, 120, 210);
                int prev_x = -1, prev_y = -1;
                for (int i = 0; i < bins; ++i) {
                    float v = -vmax_disp + ((float)i + 0.5f) * dv;
                    float expected = (float)nseg * gui_gauss_pdf(v, mu_y, sig_y) * dv;
                    int bh = (int)((expected * (float)inner_h) / max_vy);
                    if (bh < 0) bh = 0;
                    if (bh > inner_h) bh = inner_h;
                    int cx = (int)floorf(xl + (((float)i + 0.5f) / (float)bins) * seg_w);
                    int cy = inner_y0 + (inner_h - bh);
                    if (prev_x >= 0) {
                        if ((i % 2) == 0) SDL_RenderDrawLine(_renderer, prev_x, prev_y, cx, cy);
                    }
                    prev_x = cx; prev_y = cy;
                }
            }
        }
    }

    // Small legend (top-left of panel)
    if (font) {
        SDL_Color txt = {220, 220, 220, 255};
        char label[128];
        if (gui_vel_range_freeze && gui_vel_range_vmax_fixed > 0.0f) {
            snprintf(label, sizeof(label), "vx/vy  FIXED vmax=%.2f (f)%s%s%s",
                     (double)gui_vel_range_vmax_fixed,
                     (gui_vel_fit_overlay && !gui_hist_exact) ? "  fit(c)" : "",
                     gui_hist_exact ? "  EXACT(z)" : "",
                     gui_dist_markers ? "" : "  markers(off)");
        } else {
            snprintf(label, sizeof(label), "vx/vy  AUTO (f)%s%s%s",
                     (gui_vel_fit_overlay && !gui_hist_exact) ? "  fit(c)" : "",
                     gui_hist_exact ? "  EXACT(z)" : "",
                     gui_dist_markers ? "" : "  markers(off)");
        }
        draw_text(_renderer, font, label, panel_left + 8, y0 + 4, txt);
        SDL_SetRenderDrawColor(_renderer, 180, 180, 180, 220);
        SDL_RenderDrawLine(_renderer, panel_left + 30, y0 + 12, panel_left + 60, y0 + 12);
        draw_text(_renderer, font, "vy (yellow)", panel_left + 180, y0 + 4, txt);
        SDL_SetRenderDrawColor(_renderer, 255, 220, 0, 220);
        SDL_RenderDrawLine(_renderer, panel_left + 170, y0 + 12, panel_left + 200, y0 + 12);

        // Mean/median legend (solid line = mean; dashed ticks = median)
        if (gui_dist_markers) {
            draw_text(_renderer, font, "mean|median", panel_left + 300, y0 + 4, txt);
            SDL_SetRenderDrawColor(_renderer, 235, 235, 235, 170);
            SDL_RenderDrawLine(_renderer, panel_left + 292, y0 + 6, panel_left + 292, y0 + 18);
            SDL_SetRenderDrawColor(_renderer, 255, 255, 255, 120);
            for (int yy = y0 + 6; yy < y0 + 18; yy += 5) SDL_RenderDrawLine(_renderer, panel_left + 284, yy, panel_left + 284, yy + 2);
        }
    }
}

static void render_velocity_histograms_by_segment(void) {
    const int y0 = YW1 + SIM_HEIGHT + 5;
    const int h = HIST_HEIGHT - 10;
    render_velocity_histograms_by_segment_panel(y0, h);
}

// Per-segment x-position (density) histograms drawn in the bottom panel, aligned to compartments.
// Useful to visually check if each compartment has roughly-uniform density once equilibrated.
static void render_density_histograms_by_segment_panel(int y0, int h) {
    if (!segment_counts || segment_count <= 0 || segment_count > MAX_DIVIDER_CAPACITY + 1) return;
    if (particles_active <= 0) return;
    if (h < 20) return;

    SDL_SetRenderDrawBlendMode(_renderer, SDL_BLENDMODE_BLEND);
    sync_all_wall_positions();

    float left_bounds[MAX_DIVIDER_CAPACITY + 1];
    float right_bounds[MAX_DIVIDER_CAPACITY + 1];
    compute_segment_bounds_display(left_bounds, right_bounds);

    enum { GUI_HIST_MAX_BINS = 256, GUI_HIST_MAX_SEGS = MAX_DIVIDER_CAPACITY + 1 };
    const float alpha = 0.25f;

    // Panel background (clip to current piston window)
    const float pistonL_x_disp = piston_left_x + vx_piston_left * gui_render_dt_rem;
    const float pistonR_x_disp = piston_right_x + vx_piston_right * gui_render_dt_rem;
    const int panel_left = (int)fmaxf((float)XW1, pistonL_x_disp + 5.0f);
    const int panel_right = (int)fminf((float)XW2, pistonR_x_disp);
    const int panel_w = panel_right - panel_left;
    if (panel_w < 20) return;

    SDL_SetRenderDrawColor(_renderer, 20, 20, 20, 220);
    SDL_Rect bg = { panel_left, y0, panel_w, h };
    SDL_RenderFillRect(_renderer, &bg);

    // Wall boundary lines
    if (gui_dist_markers) {
        SDL_SetRenderDrawColor(_renderer, 255, 255, 255, 40);
        if (num_internal_walls > 0) {
            float tmp_walls[MAX_DIVIDER_CAPACITY];
            get_sorted_wall_positions_display(tmp_walls);
            for (int w = 0; w < num_internal_walls; ++w) {
                int xw = (int)tmp_walls[w];
                if (xw < panel_left || xw > panel_right) continue;
                SDL_RenderDrawLine(_renderer, xw, y0, xw, y0 + h);
            }
        }
        // Piston boundary line (right)
        SDL_SetRenderDrawColor(_renderer, 0, 255, 0, 80);
        SDL_RenderDrawLine(_renderer, panel_right, y0, panel_right, y0 + h);

        // Requested: green anchor line for the effective left boundary
        {
            float ax = gui_compute_left_anchor_x_display(pistonL_x_disp);
            int xax = (int)ax;
            if (xax >= panel_left && xax <= panel_right) {
                SDL_SetRenderDrawColor(_renderer, 0, 255, 0, 50);
                SDL_RenderDrawLine(_renderer, xax, y0, xax, y0 + h);
            }
        }
    }

    int inner_y0 = y0, inner_h = h;
    if (!gui_panel_inner(y0, h, &inner_y0, &inner_h)) return;

    // EMA storage
    static float ema_den[GUI_HIST_MAX_SEGS][GUI_HIST_MAX_BINS];
    static float ema_den_s0[GUI_HIST_MAX_SEGS][GUI_HIST_MAX_BINS];
    static float ema_den_s1[GUI_HIST_MAX_SEGS][GUI_HIST_MAX_BINS];
    static int ema_inited = 0;
    static int ema_segs = 0;
    static int ema_bins_per_seg[GUI_HIST_MAX_SEGS];
    if (!ema_inited || ema_segs != segment_count) {
        memset(ema_den, 0, sizeof(ema_den));
        memset(ema_den_s0, 0, sizeof(ema_den_s0));
        memset(ema_den_s1, 0, sizeof(ema_den_s1));
        for (int s = 0; s < GUI_HIST_MAX_SEGS; ++s) ema_bins_per_seg[s] = 0;
        ema_inited = 1;
        ema_segs = segment_count;
    }

    const bool split_species = (cli_experiment_preset == EXPERIMENT_PRESET_SZILARD_ENGINE && sz_species);

    for (int seg = 0; seg < segment_count; ++seg) {
        const float xl = left_bounds[seg];
        const float xr = right_bounds[seg];
        const float seg_w = xr - xl;
        if (seg_w < 8.0f) continue;

        const int bins = gui_bins_for_xdensity(seg_w, GUI_HIST_MAX_BINS);
        if (ema_bins_per_seg[seg] != bins) {
            for (int i = 0; i < GUI_HIST_MAX_BINS; ++i) ema_den[seg][i] = 0.0f;
            for (int i = 0; i < GUI_HIST_MAX_BINS; ++i) ema_den_s0[seg][i] = 0.0f;
            for (int i = 0; i < GUI_HIST_MAX_BINS; ++i) ema_den_s1[seg][i] = 0.0f;
            ema_bins_per_seg[seg] = bins;
        }

        int hist[GUI_HIST_MAX_BINS];
        for (int i = 0; i < bins; ++i) hist[i] = 0;
        int hist0[GUI_HIST_MAX_BINS];
        int hist1[GUI_HIST_MAX_BINS];
        for (int i = 0; i < bins; ++i) { hist0[i] = 0; hist1[i] = 0; }

        int total = 0;
        for (int i = 0; i < particles_active; ++i) {
            int s = segment_index_for_position(X[i]);
            if (s != seg) continue;
            float u = (X[i] - xl) / (seg_w > 0.0f ? seg_w : 1.0f);
            int b = (int)(u * (float)bins);
            if (b < 0) b = 0;
            if (b >= bins) b = bins - 1;
            hist[b] += 1;
            if (split_species) {
                if (sz_species[i] == 0) hist0[b] += 1;
                else                    hist1[b] += 1;
            }
            total += 1;
        }
        if (total <= 0) continue;

        float maxc = 1.0f;
        if (!gui_hist_exact) {
            for (int i = 0; i < bins; ++i) {
                float v = (float)hist[i];
                if (ema_den[seg][i] == 0.0f && v > 0.0f) ema_den[seg][i] = v;
                else ema_den[seg][i] = (1.0f - alpha) * ema_den[seg][i] + alpha * v;
                float sv = gui_smooth7(ema_den[seg], i, bins);
                if (sv > maxc) maxc = sv;

                if (split_species) {
                    float v0 = (float)hist0[i];
                    float v1 = (float)hist1[i];
                    if (ema_den_s0[seg][i] == 0.0f && v0 > 0.0f) ema_den_s0[seg][i] = v0;
                    else ema_den_s0[seg][i] = (1.0f - alpha) * ema_den_s0[seg][i] + alpha * v0;
                    if (ema_den_s1[seg][i] == 0.0f && v1 > 0.0f) ema_den_s1[seg][i] = v1;
                    else ema_den_s1[seg][i] = (1.0f - alpha) * ema_den_s1[seg][i] + alpha * v1;
                }
            }
        } else {
            for (int i = 0; i < bins; ++i) {
                float sv = (float)hist[i];
                if (sv > maxc) maxc = sv;
            }
        }

        if (split_species) {
            // Species overlay: draw total (white) in the background, then species 0/1 on top (same baseline).
            for (int i = 0; i < bins; ++i) {
                float sv0 = gui_hist_exact ? (float)hist0[i] : gui_smooth7(ema_den_s0[seg], i, bins);
                float sv1 = gui_hist_exact ? (float)hist1[i] : gui_smooth7(ema_den_s1[seg], i, bins);
                float svt = sv0 + sv1;
                if (svt <= 0.0f) continue;

                int bht = (int)((svt * (float)inner_h) / maxc);
                int bh0 = (int)((sv0 * (float)inner_h) / maxc);
                int bh1 = (int)((sv1 * (float)inner_h) / maxc);

                int x0 = (int)floorf(xl + ((float)i / (float)bins) * seg_w);
                int x1 = (int)floorf(xl + ((float)(i + 1) / (float)bins) * seg_w);
                int bw = x1 - x0;
                if (bw < 1) bw = 1;
                if (x0 + bw > (int)xr) bw = (int)xr - x0;
                if (bw <= 0) continue;

                // Total (white) behind
                SDL_SetRenderDrawColor(_renderer, 255, 255, 255, 70);
                SDL_Rect rt = { x0, inner_y0 + (inner_h - bht), bw, bht };
                SDL_RenderFillRect(_renderer, &rt);

                // Species 0 (red-ish)
                SDL_SetRenderDrawColor(_renderer, 255, 80, 80, 150);
                SDL_Rect r0 = { x0, inner_y0 + (inner_h - bh0), bw, bh0 };
                SDL_RenderFillRect(_renderer, &r0);

                // Species 1 (blue-ish)
                SDL_SetRenderDrawColor(_renderer, 80, 120, 255, 150);
                SDL_Rect r1 = { x0, inner_y0 + (inner_h - bh1), bw, bh1 };
                SDL_RenderFillRect(_renderer, &r1);
            }
        } else {
            // Density bars: same blue->green->yellow->red palette as speed view.
            for (int i = 0; i < bins; ++i) {
                float sv = gui_hist_exact ? (float)hist[i] : gui_smooth7(ema_den[seg], i, bins);
                int bh = (int)((sv * (float)inner_h) / maxc);
                int x0 = (int)floorf(xl + ((float)i / (float)bins) * seg_w);
                int x1 = (int)floorf(xl + ((float)(i + 1) / (float)bins) * seg_w);
                int bw = x1 - x0;
                if (bw < 1) bw = 1;
                if (x0 + bw > (int)xr) bw = (int)xr - x0;
                if (bw <= 0) continue;
                set_draw_color_speed_heat(sv / fmaxf(maxc, 1e-6f), 170);
                SDL_Rect r = { x0, inner_y0 + (inner_h - bh), bw, bh };
                SDL_RenderFillRect(_renderer, &r);
            }
        }

        // Optional: overlay expected "uniform density" baseline for this segment (flat line).
        if (gui_density_uniform_overlay && !gui_hist_exact && total > 0) {
            float expected = (float)total / (float)bins;
            if (expected > maxc) expected = maxc;
            int y_line = inner_y0 + (inner_h - (int)((expected * (float)inner_h) / maxc));
            SDL_SetRenderDrawColor(_renderer, 255, 255, 255, 70);
            for (int x = (int)xl; x < (int)xr; x += 10) {
                int x1 = x + 6;
                if (x1 > (int)xr) x1 = (int)xr;
                SDL_RenderDrawLine(_renderer, x, y_line, x1, y_line);
            }
        }

        // Exact mode helper: show particle centers directly as a "rug" at the bottom of the panel.
        // This makes it unambiguous that each particle contributes to exactly one x position.
        if (gui_hist_exact) {
            SDL_SetRenderDrawColor(_renderer, 255, 255, 255, 120);
            const int tick_h = gui_clamp_int(inner_h / 6, 6, 18);
            const int y_tick0 = inner_y0 + inner_h - 1;
            for (int i = 0; i < particles_active; ++i) {
                int s = segment_index_for_position(X[i]);
                if (s != seg) continue;
                int px = (int)lroundf((float)X[i]);
                if (px < (int)xl || px > (int)xr) continue;
                SDL_RenderDrawLine(_renderer, px, y_tick0, px, y_tick0 - tick_h);
            }
        }
    }

    if (font) {
        SDL_Color txt = {220, 220, 220, 255};
        char title[128];
        snprintf(title, sizeof(title), "x-density%s%s%s",
                 split_species ? " (species-overlay)" : "",
                 gui_hist_exact ? "  EXACT(z)" : "",
                 gui_dist_markers ? "" : "  markers(off)");
        draw_text(_renderer, font, title, panel_left + 8, y0 + 4, txt);
    }
}

static void render_density_histograms_by_segment(void) {
    const int y0 = YW1 + SIM_HEIGHT + 5;
    const int h = HIST_HEIGHT - 10;
    render_density_histograms_by_segment_panel(y0, h);
}

// Per-segment speed magnitude |v| histogram (0..vmax) aligned to compartments.
static void render_speed_histograms_by_segment_panel(int y0, int h) {
    if (!segment_counts || segment_count <= 0 || segment_count > MAX_DIVIDER_CAPACITY + 1) return;
    if (particles_active <= 0) return;
    if (h < 20) return;

    SDL_SetRenderDrawBlendMode(_renderer, SDL_BLENDMODE_BLEND);
    sync_all_wall_positions();

    float left_bounds[MAX_DIVIDER_CAPACITY + 1];
    float right_bounds[MAX_DIVIDER_CAPACITY + 1];
    compute_segment_bounds_display(left_bounds, right_bounds);

    enum { GUI_HIST_MAX_BINS = 256, GUI_HIST_MAX_SEGS = MAX_DIVIDER_CAPACITY + 1 };
    const float alpha = 0.25f;
    const float vmax_hard = (float)MAX_VELOCITY;

    // Panel background (clip to piston window)
    const float pistonL_x_disp = piston_left_x + vx_piston_left * gui_render_dt_rem;
    const float pistonR_x_disp = piston_right_x + vx_piston_right * gui_render_dt_rem;
    const int panel_left = (int)fmaxf((float)XW1, pistonL_x_disp + 5.0f);
    const int panel_right = (int)fminf((float)XW2, pistonR_x_disp);
    const int panel_w = panel_right - panel_left;
    if (panel_w < 20) return;

    SDL_SetRenderDrawColor(_renderer, 20, 20, 20, 220);
    SDL_Rect bg = { panel_left, y0, panel_w, h };
    SDL_RenderFillRect(_renderer, &bg);

    // Wall boundaries
    if (gui_dist_markers) {
        SDL_SetRenderDrawColor(_renderer, 255, 255, 255, 40);
        if (num_internal_walls > 0) {
            float tmp_walls[MAX_DIVIDER_CAPACITY];
            get_sorted_wall_positions_display(tmp_walls);
            for (int w = 0; w < num_internal_walls; ++w) {
                int xw = (int)tmp_walls[w];
                if (xw < panel_left || xw > panel_right) continue;
                SDL_RenderDrawLine(_renderer, xw, y0, xw, y0 + h);
            }
        }
        // Piston boundary line (right)
        SDL_SetRenderDrawColor(_renderer, 0, 255, 0, 80);
        SDL_RenderDrawLine(_renderer, panel_right, y0, panel_right, y0 + h);

        // Requested: green anchor line for the effective left boundary
        {
            float ax = gui_compute_left_anchor_x_display(pistonL_x_disp);
            int xax = (int)ax;
            if (xax >= panel_left && xax <= panel_right) {
                SDL_SetRenderDrawColor(_renderer, 0, 255, 0, 50);
                SDL_RenderDrawLine(_renderer, xax, y0, xax, y0 + h);
            }
        }
    }

    int inner_y0 = y0, inner_h = h;
    if (!gui_panel_inner(y0, h, &inner_y0, &inner_h)) return;

    static float ema_sp[GUI_HIST_MAX_SEGS][GUI_HIST_MAX_BINS];
    static float ema_vmax[GUI_HIST_MAX_SEGS];
    static int ema_inited = 0;
    static int ema_segs = 0;
    static int ema_bins_per_seg[GUI_HIST_MAX_SEGS];
    if (!ema_inited || ema_segs != segment_count) {
        memset(ema_sp, 0, sizeof(ema_sp));
        for (int s = 0; s < GUI_HIST_MAX_SEGS; ++s) ema_vmax[s] = 0.0f;
        for (int s = 0; s < GUI_HIST_MAX_SEGS; ++s) ema_bins_per_seg[s] = 0;
        ema_inited = 1;
        ema_segs = segment_count;
    }

    for (int seg = 0; seg < segment_count; ++seg) {
        const float xl = left_bounds[seg];
        const float xr = right_bounds[seg];
        const float seg_w = xr - xl;
        if (seg_w < 8.0f) continue;

        // Pick a display vmax from RMS speed for this segment (auto zoom), but avoid clipping fast tails.
        double m2 = 0.0;
        double sum_sp = 0.0;
        double max_sp = 0.0;
        int nseg = 0;
        for (int i = 0; i < particles_active; ++i) {
            int s = segment_index_for_position(X[i]);
            if (s != seg) continue;
            double vx = (double)Vx[i];
            double vy = (double)Vy[i];
            m2 += vx * vx + vy * vy;
            double sp = sqrt(vx * vx + vy * vy);
            sum_sp += sp;
            if (sp > max_sp) max_sp = sp;
            nseg++;
        }
        float vmax_disp = vmax_hard;
        if (gui_mb_speed_view) {
            vmax_disp = fminf(vmax_hard, gui_mb_speed_vmax_for_render());
        } else if (gui_vel_range_freeze && gui_vel_range_vmax_fixed > 0.0f) {
            vmax_disp = fminf(vmax_hard, gui_vel_range_vmax_fixed);
        }
        if (nseg > 0) {
            double rms = sqrt(m2 / (double)nseg);
            double vnew = fmax(0.25, 6.0 * rms);
            // ensure we don't clip the current max speeds (add 10% headroom)
            vnew = fmax(vnew, 1.10 * max_sp);
            if (vnew > (double)vmax_hard) vnew = (double)vmax_hard;
            float vnewf = (float)vnew;
            if (!gui_mb_speed_view && !(gui_vel_range_freeze && gui_vel_range_vmax_fixed > 0.0f)) {
                if (ema_vmax[seg] <= 0.0f) ema_vmax[seg] = vnewf;
                else ema_vmax[seg] = 0.85f * ema_vmax[seg] + 0.15f * vnewf;
                vmax_disp = ema_vmax[seg];
            }
        } else if (ema_vmax[seg] > 0.0f) {
            if (!gui_mb_speed_view && !(gui_vel_range_freeze && gui_vel_range_vmax_fixed > 0.0f)) {
                vmax_disp = ema_vmax[seg];
            }
        }
        if (vmax_disp < 1e-6f) vmax_disp = 1e-6f;

        const int bins = gui_bins_for_velocity(nseg, seg_w, GUI_HIST_MAX_BINS);
        if (ema_bins_per_seg[seg] != bins) {
            for (int i = 0; i < GUI_HIST_MAX_BINS; ++i) ema_sp[seg][i] = 0.0f;
            ema_bins_per_seg[seg] = bins;
        }

        int hist[GUI_HIST_MAX_BINS];
        for (int i = 0; i < bins; ++i) hist[i] = 0;
        for (int i = 0; i < particles_active; ++i) {
            int s = segment_index_for_position(X[i]);
            if (s != seg) continue;
            float vx = Vx[i];
            float vy = Vy[i];
            float sp = sqrtf(vx * vx + vy * vy);
            int b = (int)((sp / vmax_disp) * (float)bins);
            if (b < 0) b = 0;
            if (b >= bins) b = bins - 1;
            hist[b] += 1;
        }

        float maxc = 1.0f;
        if (!gui_hist_exact) {
            for (int i = 0; i < bins; ++i) {
                float v = (float)hist[i];
                if (ema_sp[seg][i] == 0.0f && v > 0.0f) ema_sp[seg][i] = v;
                else ema_sp[seg][i] = (1.0f - alpha) * ema_sp[seg][i] + alpha * v;
                float sv = gui_smooth7(ema_sp[seg], i, bins);
                if (sv > maxc) maxc = sv;
            }
        } else {
            for (int i = 0; i < bins; ++i) {
                float sv = (float)hist[i];
                if (sv > maxc) maxc = sv;
            }
        }

        for (int i = 0; i < bins; ++i) {
            float sv = gui_hist_exact ? (float)hist[i] : gui_smooth7(ema_sp[seg], i, bins);
            int bh = (int)((sv * (float)inner_h) / maxc);
            int x0 = (int)floorf(xl + ((float)i / (float)bins) * seg_w);
            int x1 = (int)floorf(xl + ((float)(i + 1) / (float)bins) * seg_w);
            int bw = x1 - x0;
            if (bw < 1) bw = 1;
            if (x0 + bw > (int)xr) bw = (int)xr - x0;
            if (bw <= 0) continue;
            if (gui_mb_speed_view) {
                set_draw_color_speed_heat(((float)i + 0.5f) / (float)bins, 210);
            } else {
                SDL_SetRenderDrawColor(_renderer, 60, 220, 120, 160);
            }
            SDL_Rect r = { x0, inner_y0 + (inner_h - bh), bw, bh };
            SDL_RenderFillRect(_renderer, &r);
        }

        // Mean/median markers for |v| (requested: vertical marker lines + quick intuition for heating).
        if (gui_dist_markers && nseg > 0) {
            // mean(|v|)
            const float mu = (float)(sum_sp / (double)nseg);
            int x_mu = (int)floorf(xl + (fminf(fmaxf(mu / vmax_disp, 0.0f), 1.0f) * seg_w));
            if (x_mu >= (int)xl && x_mu <= (int)xr) {
                SDL_SetRenderDrawColor(_renderer, 240, 255, 240, 110);
                SDL_RenderDrawLine(_renderer, x_mu, inner_y0, x_mu, inner_y0 + inner_h);
            }

            // median(|v|) from histogram (no per-particle sort)
            int cum = 0, med_bin = 0;
            for (int i = 0; i < bins; ++i) { cum += hist[i]; if (cum * 2 >= nseg) { med_bin = i; break; } }
            int x_med = (int)floorf(xl + (((float)med_bin + 0.5f) / (float)bins) * seg_w);
            if (x_med >= (int)xl && x_med <= (int)xr) {
                SDL_SetRenderDrawColor(_renderer, 255, 255, 255, 70);
                for (int yy = inner_y0; yy < inner_y0 + inner_h; yy += 10) {
                    SDL_RenderDrawLine(_renderer, x_med, yy, x_med, yy + 6);
                }
            }

            // In-panel stats label (requested): show distribution parameters inside each compartment region.
            if (font && seg_w >= 70.0f) {
                const double mean_sq = (double)m2 / (double)nseg; // E[|v|^2]
                const double var = fmax(0.0, mean_sq - (double)mu * (double)mu);
                const double sd = sqrt(var);
                const float med = ((float)med_bin + 0.5f) * (vmax_disp / (float)bins);

                SDL_Color ctxt = {200, 245, 210, 255};
                TTF_Font *ff = font_small ? font_small : font;
                const int dy_stats = TTF_FontLineSkip(ff);

                // Use 3 short lines (mu/med/sd) so text doesn't get horizontally cut off
                // when compartments shrink during piston/wall motion.
                char l1[64], l2[64], l3[64];
                snprintf(l1, sizeof(l1), "mu=%.2f", (double)mu);
                snprintf(l2, sizeof(l2), "med=%.2f", (double)med);
                snprintf(l3, sizeof(l3), "sd=%.2f", sd);

                int tw1 = 0, th1 = 0, tw2 = 0, th2 = 0, tw3 = 0, th3 = 0;
                TTF_SizeText(ff, l1, &tw1, &th1);
                TTF_SizeText(ff, l2, &tw2, &th2);
                TTF_SizeText(ff, l3, &tw3, &th3);
                int bg_w = tw1;
                if (tw2 > bg_w) bg_w = tw2;
                if (tw3 > bg_w) bg_w = tw3;
                bg_w += 8;
                if (bg_w > (int)seg_w - 8) bg_w = (int)seg_w - 8;
                if (bg_w < 40) bg_w = 40;
                int bg_h = dy_stats * 3 + 4;

                int tx = (int)xl + 4;
                int ty = inner_y0 + 2; // inside bar area: doesn't require more top padding
                if (ty + bg_h > inner_y0 + inner_h) ty = inner_y0 + 2;

                SDL_SetRenderDrawColor(_renderer, 0, 0, 0, 120);
                SDL_Rect bg = { tx - 2, ty - 2, bg_w, bg_h };
                SDL_RenderFillRect(_renderer, &bg);

                draw_text_clipped(_renderer, ff, l1, tx, ty + 0 * dy_stats, ctxt, (int)seg_w - 8);
                draw_text_clipped(_renderer, ff, l2, tx, ty + 1 * dy_stats, ctxt, (int)seg_w - 8);
                draw_text_clipped(_renderer, ff, l3, tx, ty + 2 * dy_stats, ctxt, (int)seg_w - 8);
            }
        }

        // Optional: overlay a fitted Rayleigh curve on |v| (visual MB sanity check in 2D).
        // Rayleigh sigma is derived from rms speed: rms = sqrt(E[v^2]) = sqrt(2)*sigma.
        if ((gui_vel_fit_overlay || gui_mb_speed_view) && !gui_hist_exact && nseg >= 10) {
            const float dv = vmax_disp / (float)bins;
            float sigma = 0.0f;
            if (nseg > 0) {
                double rms = sqrt(m2 / (double)nseg);
                sigma = (float)(rms / sqrt(2.0));
            }
            if (sigma > 0.0f) {
                SDL_SetRenderDrawColor(_renderer, 225, 235, 225, gui_mb_speed_view ? 245 : 210);
                int prev_x = -1, prev_y = -1;
                for (int i = 0; i < bins; ++i) {
                    float v = ((float)i + 0.5f) * dv;
                    float expected = (float)nseg * gui_rayleigh_pdf(v, sigma) * dv;
                    int bh = (int)((expected * (float)inner_h) / maxc);
                    if (bh < 0) bh = 0;
                    if (bh > inner_h) bh = inner_h;
                    int cx = (int)floorf(xl + (((float)i + 0.5f) / (float)bins) * seg_w);
                    int cy = inner_y0 + (inner_h - bh);
                    if (prev_x >= 0) {
                        SDL_RenderDrawLine(_renderer, prev_x, prev_y, cx, cy);
                        if (gui_mb_speed_view) SDL_RenderDrawLine(_renderer, prev_x, prev_y + 1, cx, cy + 1);
                    }
                    prev_x = cx; prev_y = cy;
                }
            }
        }
    }

    if (font) {
        SDL_Color txt = {220, 220, 220, 255};
        char label[128];
        const char *view_label = gui_mb_speed_view ? "speed |v| (MB colors)" : "speed |v| (green)";
        const char *curve_label = (gui_mb_speed_view && !gui_hist_exact) ? "  MB curve"
                                : (gui_vel_fit_overlay && !gui_hist_exact) ? "  fit(c)"
                                : "";
        if (gui_vel_range_freeze && gui_vel_range_vmax_fixed > 0.0f) {
            snprintf(label, sizeof(label), "%s  FIXED vmax=%.2f (f)%s%s%s",
                     view_label,
                     (double)gui_vel_range_vmax_fixed,
                     curve_label,
                     gui_hist_exact ? "  EXACT(z)" : "",
                     gui_dist_markers ? "" : "  markers(off)");
        } else {
            snprintf(label, sizeof(label), "%s  AUTO (f)%s%s%s",
                     view_label,
                     curve_label,
                     gui_hist_exact ? "  EXACT(z)" : "",
                     gui_dist_markers ? "" : "  markers(off)");
        }
        draw_text(_renderer, font, label, panel_left + 8, y0 + 4, txt);

        // Mean/median legend (solid line = mean; dashed ticks = median)
        if (gui_dist_markers) {
            draw_text(_renderer, font, "mean|median", panel_left + 300, y0 + 4, txt);
            SDL_SetRenderDrawColor(_renderer, 220, 255, 220, 170);
            SDL_RenderDrawLine(_renderer, panel_left + 292, y0 + 6, panel_left + 292, y0 + 18);
            SDL_SetRenderDrawColor(_renderer, 255, 255, 255, 120);
            for (int yy = y0 + 6; yy < y0 + 18; yy += 5) SDL_RenderDrawLine(_renderer, panel_left + 284, yy, panel_left + 284, yy + 2);
        }
    }
}

// Global x-density histogram across the whole piston window (no per-segment separation).
// Useful as a "first view" to see overall compression/rarefaction across the chain.
static void render_density_histogram_global_panel(int y0, int h) {
    if (particles_active <= 0) return;
    if (h < 20) return;

    SDL_SetRenderDrawBlendMode(_renderer, SDL_BLENDMODE_BLEND);
    sync_all_wall_positions();

    enum { GUI_HIST_MAX_BINS = 512 };
    const float alpha = 0.25f;

    const float pistonL_x_disp = piston_left_x + vx_piston_left * gui_render_dt_rem;
    const float pistonR_x_disp = piston_right_x + vx_piston_right * gui_render_dt_rem;
    const int panel_left = (int)fmaxf((float)XW1, pistonL_x_disp + 5.0f);
    const int panel_right = (int)fminf((float)XW2, pistonR_x_disp);
    const int panel_w = panel_right - panel_left;
    if (panel_w < 20) return;
    const int bins = gui_bins_for_xdensity((float)panel_w, GUI_HIST_MAX_BINS);

    SDL_SetRenderDrawColor(_renderer, 20, 20, 20, 220);
    SDL_Rect bg = { panel_left, y0, panel_w, h };
    SDL_RenderFillRect(_renderer, &bg);

    if (gui_dist_markers) {
        // Wall boundaries (display positions)
        SDL_SetRenderDrawColor(_renderer, 255, 255, 255, 40);
        if (num_internal_walls > 0) {
            float tmp_walls[MAX_DIVIDER_CAPACITY];
            get_sorted_wall_positions_display(tmp_walls);
            for (int w = 0; w < num_internal_walls; ++w) {
                int xw = (int)tmp_walls[w];
                if (xw < panel_left || xw > panel_right) continue;
                SDL_RenderDrawLine(_renderer, xw, y0, xw, y0 + h);
            }
        }
        SDL_SetRenderDrawColor(_renderer, 0, 255, 0, 80);
        SDL_RenderDrawLine(_renderer, panel_right, y0, panel_right, y0 + h);

        // Requested: green anchor line for the effective left boundary
        {
            float ax = gui_compute_left_anchor_x_display(pistonL_x_disp);
            int xax = (int)ax;
            if (xax >= panel_left && xax <= panel_right) {
                SDL_SetRenderDrawColor(_renderer, 0, 255, 0, 50);
                SDL_RenderDrawLine(_renderer, xax, y0, xax, y0 + h);
            }
        }
    }

    int inner_y0 = y0, inner_h = h;
    if (!gui_panel_inner(y0, h, &inner_y0, &inner_h)) return;

    static float ema[GUI_HIST_MAX_BINS];
    static float ema0[GUI_HIST_MAX_BINS];
    static float ema1[GUI_HIST_MAX_BINS];
    static int ema_inited = 0;
    static int ema_bins = 0;
    if (!ema_inited || ema_bins != bins) {
        memset(ema, 0, sizeof(ema));
        memset(ema0, 0, sizeof(ema0));
        memset(ema1, 0, sizeof(ema1));
        ema_inited = 1;
        ema_bins = bins;
    }

    int hist[GUI_HIST_MAX_BINS];
    for (int i = 0; i < bins; ++i) hist[i] = 0;
    int hist0[GUI_HIST_MAX_BINS];
    int hist1[GUI_HIST_MAX_BINS];
    for (int i = 0; i < bins; ++i) { hist0[i] = 0; hist1[i] = 0; }

    const bool split_species = (cli_experiment_preset == EXPERIMENT_PRESET_SZILARD_ENGINE && sz_species);

    const float xl = (float)panel_left;
    const float xr = (float)panel_right;
    const float w = xr - xl;
    if (w <= 1.0f) return;

    int nwin = 0;
    for (int i = 0; i < particles_active; ++i) {
        float x = (float)X[i];
        if (x < xl || x >= xr) continue;
        int b = (int)(((x - xl) / w) * (float)bins);
        if (b < 0) b = 0;
        if (b >= bins) b = bins - 1;
        hist[b] += 1;
        if (split_species) {
            if (sz_species[i] == 0) hist0[b] += 1;
            else                    hist1[b] += 1;
        }
        nwin += 1;
    }

    float maxc = 1.0f;
    if (!gui_hist_exact) {
        for (int i = 0; i < bins; ++i) {
            float v = (float)hist[i];
            if (ema[i] == 0.0f && v > 0.0f) ema[i] = v;
            else ema[i] = (1.0f - alpha) * ema[i] + alpha * v;
            float sv = gui_smooth7(ema, i, bins);
            if (sv > maxc) maxc = sv;

            if (split_species) {
                float v0 = (float)hist0[i];
                float v1 = (float)hist1[i];
                if (ema0[i] == 0.0f && v0 > 0.0f) ema0[i] = v0;
                else ema0[i] = (1.0f - alpha) * ema0[i] + alpha * v0;
                if (ema1[i] == 0.0f && v1 > 0.0f) ema1[i] = v1;
                else ema1[i] = (1.0f - alpha) * ema1[i] + alpha * v1;
            }
        }
    } else {
        for (int i = 0; i < bins; ++i) {
            float sv = (float)hist[i];
            if (sv > maxc) maxc = sv;
        }
    }

    if (split_species) {
        // Species overlay: draw total (white) in the background, then species 0/1 on top (same baseline).
        for (int i = 0; i < bins; ++i) {
            float sv0 = gui_hist_exact ? (float)hist0[i] : gui_smooth7(ema0, i, bins);
            float sv1 = gui_hist_exact ? (float)hist1[i] : gui_smooth7(ema1, i, bins);
            float svt = sv0 + sv1;
            if (svt <= 0.0f) continue;
            int bht = (int)((svt * (float)inner_h) / maxc);
            int bh0 = (int)((sv0 * (float)inner_h) / maxc);
            int bh1 = (int)((sv1 * (float)inner_h) / maxc);
            int x0 = (int)floorf(xl + ((float)i / (float)bins) * w);
            int x1 = (int)floorf(xl + ((float)(i + 1) / (float)bins) * w);
            int bw = x1 - x0;
            if (bw < 1) bw = 1;
            if (x0 + bw > (int)xr) bw = (int)xr - x0;
            if (bw <= 0) continue;

            SDL_SetRenderDrawColor(_renderer, 255, 255, 255, 70);
            SDL_Rect rt = { x0, inner_y0 + (inner_h - bht), bw, bht };
            SDL_RenderFillRect(_renderer, &rt);

            SDL_SetRenderDrawColor(_renderer, 255, 80, 80, 150);
            SDL_Rect r0 = { x0, inner_y0 + (inner_h - bh0), bw, bh0 };
            SDL_RenderFillRect(_renderer, &r0);

            SDL_SetRenderDrawColor(_renderer, 80, 120, 255, 150);
            SDL_Rect r1 = { x0, inner_y0 + (inner_h - bh1), bw, bh1 };
            SDL_RenderFillRect(_renderer, &r1);
        }
    } else {
        // Global density: same blue->green->yellow->red palette as speed view.
        for (int i = 0; i < bins; ++i) {
            float sv = gui_hist_exact ? (float)hist[i] : gui_smooth7(ema, i, bins);
            int bh = (int)((sv * (float)inner_h) / maxc);
            int x0 = (int)floorf(xl + ((float)i / (float)bins) * w);
            int x1 = (int)floorf(xl + ((float)(i + 1) / (float)bins) * w);
            int bw = x1 - x0;
            if (bw < 1) bw = 1;
            if (x0 + bw > (int)xr) bw = (int)xr - x0;
            if (bw <= 0) continue;
            set_draw_color_speed_heat(sv / fmaxf(maxc, 1e-6f), 165);
            SDL_Rect r = { x0, inner_y0 + (inner_h - bh), bw, bh };
            SDL_RenderFillRect(_renderer, &r);
        }
    }

    // Optional: overlay expected "uniform density" baseline (flat line) across the whole piston window.
    // This makes it obvious when the system is not well-mixed (e.g., after compression).
    if (gui_density_uniform_overlay && !gui_hist_exact && nwin > 0) {
        float expected = (float)nwin / (float)bins;
        if (expected > maxc) expected = maxc;
        int y_line = inner_y0 + (inner_h - (int)((expected * (float)inner_h) / maxc));
        SDL_SetRenderDrawColor(_renderer, 255, 255, 255, 70);
        for (int x = (int)xl; x < (int)xr; x += 10) {
            int x1 = x + 6;
            if (x1 > (int)xr) x1 = (int)xr;
            SDL_RenderDrawLine(_renderer, x, y_line, x1, y_line);
        }
    }

    // Optional: overlay a smoothed density curve across the whole window (visual "fit", like the velocity panels).
    // This helps when transparent bars look like "just a line" in some lighting/monitor setups.
    if (gui_density_curve_overlay && !gui_hist_exact && !split_species) {
        SDL_SetRenderDrawColor(_renderer, 255, 200, 200, 220);
        int prev_x = -1, prev_y = -1;
        for (int i = 0; i < bins; ++i) {
            float sv = gui_smooth7(ema, i, bins);
            int bh = (int)((sv * (float)inner_h) / maxc);
            if (bh < 0) bh = 0;
            if (bh > inner_h) bh = inner_h;
            int cx = (int)floorf(xl + (((float)i + 0.5f) / (float)bins) * w);
            int cy = inner_y0 + (inner_h - bh);
            if (prev_x >= 0) {
                SDL_RenderDrawLine(_renderer, prev_x, prev_y, cx, cy);
                // fake thickness = 2px (helps on high-DPI monitors)
                SDL_RenderDrawLine(_renderer, prev_x, prev_y + 1, cx, cy + 1);
            }
            prev_x = cx; prev_y = cy;
        }
    }

    // Exact mode helper: show particle centers directly as a "rug" at the bottom of the panel.
    if (gui_hist_exact && nwin > 0) {
        const int tick_h = gui_clamp_int(inner_h / 6, 6, 18);
        const int y_tick0 = inner_y0 + inner_h - 1;
        for (int i = 0; i < particles_active; ++i) {
            float x = (float)X[i];
            if (x < xl || x > xr) continue;
            int px = (int)lroundf(x);
            if (px < (int)xl || px > (int)xr) continue;
            if (split_species) {
                if (sz_species[i] == 0) SDL_SetRenderDrawColor(_renderer, 255, 80, 80, 170);
                else                    SDL_SetRenderDrawColor(_renderer, 80, 120, 255, 170);
            } else {
                SDL_SetRenderDrawColor(_renderer, 255, 255, 255, 120);
            }
            SDL_RenderDrawLine(_renderer, px, y_tick0, px, y_tick0 - tick_h);
        }
    }

    if (font) {
        SDL_Color txt = {220, 220, 220, 255};
        char title[128];
        snprintf(title, sizeof(title), "x-density (global)%s%s%s%s%s",
                 split_species ? "  species" : "",
                 gui_hist_exact ? "  EXACT(z)" : "",
                 (gui_density_uniform_overlay && !gui_hist_exact) ? "  baseline(u)" : "",
                 (gui_density_curve_overlay && !gui_hist_exact && !split_species) ? "  curve(y)" : "",
                 gui_dist_markers ? "" : "  markers(off)");
        draw_text(_renderer, font, title, panel_left + 8, y0 + 4, txt);

        if (!split_species && nwin > 0) {
            const double expected = (double)nwin / (double)bins;
            const double contrast = (expected > 1e-12) ? fmax(0.0, (double)maxc / expected - 1.0) : 0.0;
            const double T_live = fmax(0.0, (double)compute_measured_temperature_from_ke());
            const double c_s = sqrt(fmax(0.0, 2.0 * (double)kB_effective() * T_live / (double)PARTICLE_MASS));
            char wave_label[160];
            snprintf(wave_label, sizeof(wave_label), "wave watch: rho_tt=cs^2 rho_xx  cs~%.2f  contrast=%.2f",
                     c_s, contrast);
            draw_text_clipped(_renderer, font_small ? font_small : font, wave_label,
                              panel_left + 8, y0 + 20, (SDL_Color){190, 205, 210, 255}, panel_w - 16);
        }
    }
}

// Combined distributions stacked vertically (requested row order):
// 1) global x-density across the whole piston window
// 2) per-segment x-density (aligned to compartments)
// 3) per-segment speed magnitude |v|
// 4) per-segment velocity components vx/vy
static void render_segment_distributions_stacked(void) {
    const int total_h = HIST_HEIGHT - 10;
    if (total_h < 70) {
        render_velocity_histograms_by_segment();
        return;
    }
    const int gap = 8;
    const int y_top = YW1 + SIM_HEIGHT + 5;

    if (total_h >= 170) {
        if (gui_show_x_density_by_segment) {
            // 4 panels (equal heights so velocity panels aren't huge vs density panels)
            const int gaps = gap * 3;
            const int base = total_h - gaps;
            const int h0 = base / 4;             // global density
            const int h1 = base / 4;             // per-seg density
            const int h2 = base / 4;             // speed
            const int h3 = base - h0 - h1 - h2;  // vx/vy (gets the remainder)
            const int y0 = y_top;
            const int y1 = y0 + h0 + gap;
            const int y2 = y1 + h1 + gap;
            const int y3 = y2 + h2 + gap;
            render_density_histogram_global_panel(y0, h0);
            render_density_histograms_by_segment_panel(y1, h1);
            render_speed_histograms_by_segment_panel(y2, h2);
            render_velocity_histograms_by_segment_panel(y3, h3);
        } else {
            // 3 panels: global density + speed + vx/vy (more compact; per-compartment density can be misleading)
            const int gaps = gap * 2;
            const int base = total_h - gaps;
            const int h0 = base / 3;             // global density
            const int h1 = base / 3;             // speed
            const int h2 = base - h0 - h1;       // vx/vy
            const int y0 = y_top;
            const int y1 = y0 + h0 + gap;
            const int y2 = y1 + h1 + gap;
            render_density_histogram_global_panel(y0, h0);
            render_speed_histograms_by_segment_panel(y1, h1);
            render_velocity_histograms_by_segment_panel(y2, h2);
        }
    } else if (total_h >= 110) {
        // 3 panels in tight space:
        // - If per-compartment x-density enabled: global density + per-compartment density + vx/vy
        // - Else: global density + speed + vx/vy
        const int gaps = gap * 2;
        const int h0 = (total_h - gaps) / 3;
        const int h1 = h0;
        const int h2 = total_h - gaps - h0 - h1;
        const int y0 = y_top;
        const int y1 = y0 + h0 + gap;
        const int y2 = y1 + h1 + gap;
        render_density_histogram_global_panel(y0, h0);
        if (gui_show_x_density_by_segment) {
            render_density_histograms_by_segment_panel(y1, h1);
        } else {
            render_speed_histograms_by_segment_panel(y1, h1);
        }
        render_velocity_histograms_by_segment_panel(y2, h2);
    } else {
        // 2 panels
        const int h1 = (total_h - gap) / 2;
        const int h2 = total_h - gap - h1;
        const int y1 = y_top;
        const int y2 = y1 + h1 + gap;
        render_density_histogram_global_panel(y1, h1);
        render_velocity_histograms_by_segment_panel(y2, h2);
    }
}










void draw_clear_screen() {
    SDL_SetRenderDrawColor(_renderer, 0, 0, 0, 255); // Black color
    SDL_RenderClear(_renderer);
}

void draw_simulation_boundary() {
    // ##CHRIS: Draw red walls when heat bath is enabled, yellow otherwise
    if (heatbath_enabled) {
        SDL_SetRenderDrawColor(_renderer, 255, 0, 0, 255); // Red color (heat bath active)
    } else {
        SDL_SetRenderDrawColor(_renderer, 255, 255, 0, 255); // Yellow color (normal)
    }
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

// Draw text truncated to fit max_w pixels (adds "..." when truncated).
static void draw_text_clipped(SDL_Renderer *renderer, TTF_Font *font, const char *text, int x, int y, SDL_Color color, int max_w) {
    if (!text) return;
    if (max_w <= 0) { draw_text(renderer, font, text, x, y, color); return; }

    int w = 0, h = 0;
    if (TTF_SizeText(font, text, &w, &h) != 0 || w <= max_w) {
        draw_text(renderer, font, text, x, y, color);
        return;
    }

    const char *ellipsis = "...";
    int w_ell = 0, h_ell = 0;
    if (TTF_SizeText(font, ellipsis, &w_ell, &h_ell) != 0) {
        draw_text(renderer, font, text, x, y, color);
        return;
    }
    if (w_ell >= max_w) { draw_text(renderer, font, ellipsis, x, y, color); return; }

    size_t n = strlen(text);
    char tmp[1024];
    if (n >= sizeof(tmp)) n = sizeof(tmp) - 1;

    size_t lo = 0, hi = n;
    while (lo + 1 < hi) {
        size_t mid = (lo + hi) / 2;
        memcpy(tmp, text, mid);
        tmp[mid] = '\0';
        strncat(tmp, ellipsis, sizeof(tmp) - strlen(tmp) - 1);
        int ww = 0, hh = 0;
        if (TTF_SizeText(font, tmp, &ww, &hh) == 0 && ww <= max_w) lo = mid;
        else hi = mid;
    }
    memcpy(tmp, text, lo);
    tmp[lo] = '\0';
    strncat(tmp, ellipsis, sizeof(tmp) - strlen(tmp) - 1);
    draw_text(renderer, font, tmp, x, y, color);
}

static double live_clamp(double value, double lo, double hi) {
    if (value < lo) return lo;
    if (value > hi) return hi;
    return value;
}

static double live_slider_get_value(const LiveSlider *slider) {
    if (!slider || !slider->target) return 0.0;
    if (slider->type == LIVE_PARAM_DOUBLE) return *((double *)slider->target);
    return (double)(*((float *)slider->target));
}

static void live_slider_set_value(int slider_index, double value) {
    if (slider_index < 0 || slider_index >= live_slider_count) return;
    LiveSlider *slider = &live_sliders[slider_index];
    value = live_clamp(value, slider->min_value, slider->max_value);
    if (slider->type == LIVE_PARAM_DOUBLE) *((double *)slider->target) = value;
    else *((float *)slider->target) = (float)value;

    // Keep values sane when the user drags aggressively.
    if (&time_scale_runtime == slider->target) {
        time_scale_runtime = fmaxf(0.01f, time_scale_runtime);
    } else if (&fixed_dt_runtime == slider->target) {
        fixed_dt_runtime = fmaxf(1e-6f, fixed_dt_runtime);
    } else if (&piston_protocol_max_speed == slider->target) {
        piston_protocol_max_speed = fmaxf(1e-6f, fabsf(piston_protocol_max_speed));
    } else if (&piston_protocol_duration == slider->target) {
        piston_protocol_duration = fmaxf(1e-6f, piston_protocol_duration);
    } else if (&heatbath_temperature == slider->target) {
        heatbath_temperature = fmaxf(1e-6f, heatbath_temperature);
        if (g_edmd) {
            ((EDMD_Params *)edmd_params(g_edmd))->heatbath_temperature = (double)heatbath_temperature;
        }
    } else if (&gui_particle_draw_min_px == slider->target) {
        gui_particle_draw_min_px = fmaxf(0.0f, gui_particle_draw_min_px);
    }
}

static void live_sync_edmd_spring_settings(void) {
    if (!g_edmd || num_internal_walls <= 0) return;

    EDMD_Params *ep = edmd_backend_params_mut(g_edmd);
    if (!ep || ep->divider_count <= 0) return;

    int primary_idx = primary_wall_index;
    if (primary_idx < 0 || primary_idx >= ep->divider_count) primary_idx = 0;

    if (energy_measurement.spring_enabled && energy_measurement.eff_output_mode == EFF_OUTPUT_SPRING) {
        const double xeq_px = (energy_measurement.equilibrium_position != 0.0f)
                                ? (double)energy_measurement.equilibrium_position
                                : ((double)XW1 + ep->divider_x[primary_idx]);
        ep->divider_k[primary_idx] = (double)energy_measurement.spring_constant;
        ep->divider_xeq[primary_idx] = xeq_px - (double)XW1;
    } else {
        ep->divider_k[primary_idx] = 0.0;
    }

    edmd_reschedule_all(g_edmd);
}

static void live_apply_spring_enabled(int enabled) {
    if (enabled) {
        if (!energy_measurement.eff_enabled) {
            cli_enable_energy_measurement = true;
            cli_eff_output_mode = EFF_OUTPUT_SPRING;
            initialize_energy_measurement();
        }
        energy_measurement.eff_enabled = true;
        energy_measurement.eff_output_mode = EFF_OUTPUT_SPRING;
        energy_measurement.spring_enabled = true;
    } else {
        energy_measurement.spring_enabled = false;
        energy_measurement.spring_force = 0.0f;
        energy_measurement.spring_energy = 0.0f;
    }

    live_sync_edmd_spring_settings();

    if (!cli_quiet) {
        printf("Live spring force: %s\n", energy_measurement.spring_enabled ? "ON" : "off");
    }
}

static double live_current_nominal_eta(void) {
    const double r_sigma = (double)PARTICLE_RADIUS / (double)PIXELS_PER_SIGMA;
    const double area = 2.0 * (double)L0_UNITS * (double)HEIGHT_UNITS;
    if (particles_active <= 0 || area <= 0.0) return 0.0;
    return (double)particles_active * M_PI * r_sigma * r_sigma / area;
}

static int live_scene_particle_max(void) {
    int max_n = NUM_PARTICLES;
    if (max_n > 20000) max_n = 20000; // SDL exploration cap; headless sweeps can still use CLI.
    if (max_n < 1) max_n = 1;
    return max_n;
}

static const char *live_scene_mode_label(LiveSceneMode mode) {
    switch (mode) {
        case LIVE_SCENE_MODE_ENERGY_TRANSFER: return "energy transfer";
        case LIVE_SCENE_MODE_SZILARD: return "Szilard sliding box";
        case LIVE_SCENE_MODE_SPEED_OF_SOUND: return "speed of sound";
        default: return "energy transfer";
    }
}

static LiveSceneMode live_scene_mode_from_preset(void) {
    if (cli_experiment_preset == EXPERIMENT_PRESET_SZILARD_ENGINE) {
        return LIVE_SCENE_MODE_SZILARD;
    }
    if (cli_experiment_preset == EXPERIMENT_PRESET_SPEED_OF_SOUND ||
        cli_experiment_preset == EXPERIMENT_PRESET_WALL_MID ||
        cli_experiment_preset == EXPERIMENT_PRESET_OFFCENTER_SINGLE) {
        return LIVE_SCENE_MODE_SPEED_OF_SOUND;
    }
    return LIVE_SCENE_MODE_ENERGY_TRANSFER;
}

static int live_scene_fixed_wall_count(void) {
    if (live_scene_mode == LIVE_SCENE_MODE_SZILARD) return 1;
    if (live_scene_mode == LIVE_SCENE_MODE_SPEED_OF_SOUND) return 1;
    return 0;
}

static int live_scene_spring_available(void) {
    return live_scene_mode == LIVE_SCENE_MODE_ENERGY_TRANSFER;
}

static const char *live_scene_engine_label(void) {
    switch (live_scene_engine) {
        case LIVE_SCENE_ENGINE_EDMD: return "EDMD";
        case LIVE_SCENE_ENGINE_HYBRID: return "HYBRID PP-EDMD";
        case LIVE_SCENE_ENGINE_TIME_TOI:
        default: return "TIME+TOI";
    }
}

static void live_scene_cycle_engine(int direction) {
    if (live_scene_mode == LIVE_SCENE_MODE_SZILARD) {
        live_scene_engine = LIVE_SCENE_ENGINE_EDMD;
        return;
    }
    if (direction >= 0) {
        if (live_scene_engine == LIVE_SCENE_ENGINE_EDMD) live_scene_engine = LIVE_SCENE_ENGINE_HYBRID;
        else if (live_scene_engine == LIVE_SCENE_ENGINE_HYBRID) live_scene_engine = LIVE_SCENE_ENGINE_TIME_TOI;
        else live_scene_engine = LIVE_SCENE_ENGINE_EDMD;
    } else {
        if (live_scene_engine == LIVE_SCENE_ENGINE_EDMD) live_scene_engine = LIVE_SCENE_ENGINE_TIME_TOI;
        else if (live_scene_engine == LIVE_SCENE_ENGINE_TIME_TOI) live_scene_engine = LIVE_SCENE_ENGINE_HYBRID;
        else live_scene_engine = LIVE_SCENE_ENGINE_EDMD;
    }
}

static double live_scene_spring_l0_min_sigma(void) {
    return 1.0;
}

static double live_scene_spring_l0_max_sigma(void) {
    return fmax(live_scene_spring_l0_min_sigma() + 0.5, 2.0 * (double)L0_UNITS - 0.5);
}

static double live_scene_clamped_spring_l0_sigma(void) {
    return live_clamp((double)live_scene_spring_l0_sigma,
                      live_scene_spring_l0_min_sigma(),
                      live_scene_spring_l0_max_sigma());
}

static int live_scene_param_visible(int index) {
    if (index == 0) return live_scene_fixed_wall_count() == 0;
    if (index >= 1 && index <= 3) return 1;
    if (index == 4) return live_scene_mode == LIVE_SCENE_MODE_ENERGY_TRANSFER;
    if (index >= 5 && index < LIVE_SCENE_ROW_COUNT) {
        return live_scene_mode == LIVE_SCENE_MODE_ENERGY_TRANSFER &&
               (index - 5) < live_scene_wall_count;
    }
    return 0;
}

static int live_scene_visible_param_count(void) {
    int count = 0;
    for (int i = 0; i < LIVE_SCENE_ROW_COUNT; ++i) {
        if (live_scene_param_visible(i)) count++;
    }
    return count;
}

static int live_scene_index_for_visible_row(int visible_row) {
    int row = 0;
    for (int i = 0; i < LIVE_SCENE_ROW_COUNT; ++i) {
        if (!live_scene_param_visible(i)) continue;
        if (row == visible_row) return i;
        row++;
    }
    return -1;
}

static void live_scene_apply_mode_constraints(void) {
    const int fixed_walls = live_scene_fixed_wall_count();
    if (fixed_walls > 0) {
        live_scene_wall_count = fixed_walls;
    } else {
        if (live_scene_wall_count < 1) live_scene_wall_count = 2;
        if (live_scene_wall_count > LIVE_SCENE_MAX_WALLS) live_scene_wall_count = LIVE_SCENE_MAX_WALLS;
    }
}

static int live_scene_segments_for_walls(int walls) {
    if (walls < 1) walls = 1;
    if (walls > LIVE_SCENE_MAX_WALLS) walls = LIVE_SCENE_MAX_WALLS;
    return walls + 1;
}

static double live_scene_wall_position_sigma(int walls, int wall_index) {
    if (walls < 1) walls = 1;
    if (walls > LIVE_SCENE_MAX_WALLS) walls = LIVE_SCENE_MAX_WALLS;
    if (wall_index < 0) wall_index = 0;
    if (wall_index >= walls) wall_index = walls - 1;

    if (live_scene_mode == LIVE_SCENE_MODE_ENERGY_TRANSFER) {
        const double spring_l0 = live_scene_clamped_spring_l0_sigma();
        const double right = 2.0 * (double)L0_UNITS;
        if (wall_index == 0) return spring_l0;
        const double span = fmax(1e-6, right - spring_l0);
        return spring_l0 + span * ((double)wall_index / (double)walls);
    }

    return 2.0 * (double)L0_UNITS * (double)(wall_index + 1) / (double)(walls + 1);
}

static double live_scene_segment_width_sigma(int walls, int segment) {
    const int segments = live_scene_segments_for_walls(walls);
    if (segment < 0) segment = 0;
    if (segment >= segments) segment = segments - 1;
    const double left = (segment == 0) ? 0.0 : live_scene_wall_position_sigma(walls, segment - 1);
    const double right = (segment < walls) ? live_scene_wall_position_sigma(walls, segment)
                                           : 2.0 * (double)L0_UNITS;
    return fmax(0.0, right - left);
}

static int live_scene_start_segment_for_walls(int walls) {
    const int segments = live_scene_segments_for_walls(walls);
    int start_seg = (live_scene_mode == LIVE_SCENE_MODE_ENERGY_TRANSFER) ? 1 : 0;
    if (start_seg >= segments) start_seg = 0;
    return start_seg;
}

static void live_scene_counts_for(int walls, int n, int *counts, int max_counts) {
    if (!counts || max_counts <= 0) return;
    const int segments = live_scene_segments_for_walls(walls);
    for (int s = 0; s < max_counts; ++s) counts[s] = 0;
    if (n < 0) n = 0;

    const int start_seg = live_scene_start_segment_for_walls(walls);
    int active_segments = segments - start_seg;
    if (active_segments < 1) active_segments = 1;
    for (int s = start_seg; s < segments && s < max_counts; ++s) {
        const int idx = s - start_seg;
        int count = n / active_segments;
        if (idx < (n % active_segments)) count++;
        counts[s] = count;
    }
}

static double live_scene_radius_for_eta_value(int walls, int n, double eta) {
    eta = live_clamp(eta, 0.02, 0.75);
    if (n < 1) n = 1;

    int counts[LIVE_SCENE_MAX_WALLS + 1];
    const int segments = live_scene_segments_for_walls(walls);
    live_scene_counts_for(walls, n, counts, (int)(sizeof(counts) / sizeof(counts[0])));

    double r_min = DBL_MAX;
    for (int s = 0; s < segments; ++s) {
        const double width_sigma = live_scene_segment_width_sigma(walls, s);
        const double area_sigma = width_sigma * (double)HEIGHT_UNITS;
        if (counts[s] <= 0 || area_sigma <= 0.0) continue;
        const double r_sq = eta * area_sigma / ((double)counts[s] * M_PI);
        if (r_sq > 0.0 && isfinite(r_sq)) {
            const double r = sqrt(r_sq);
            if (r < r_min) r_min = r;
        }
    }

    if (r_min == DBL_MAX) {
        const double area = 2.0 * (double)L0_UNITS * (double)HEIGHT_UNITS;
        r_min = sqrt(fmax(0.0, eta * area / ((double)n * M_PI)));
    }
    if (!isfinite(r_min) || r_min <= 0.0) r_min = 0.01;
    return r_min;
}

static double live_scene_eta_for_radius_value(int walls, int n, double radius_sigma) {
    if (n < 1) n = 1;
    if (!isfinite(radius_sigma) || radius_sigma <= 0.0) radius_sigma = 0.01;

    int counts[LIVE_SCENE_MAX_WALLS + 1];
    const int segments = live_scene_segments_for_walls(walls);
    live_scene_counts_for(walls, n, counts, (int)(sizeof(counts) / sizeof(counts[0])));

    double eta_max = 0.0;
    for (int s = 0; s < segments; ++s) {
        const double width_sigma = live_scene_segment_width_sigma(walls, s);
        const double area_sigma = width_sigma * (double)HEIGHT_UNITS;
        if (counts[s] <= 0 || area_sigma <= 0.0) continue;
        const double eta = (double)counts[s] * M_PI * radius_sigma * radius_sigma / area_sigma;
        if (isfinite(eta) && eta > eta_max) eta_max = eta;
    }
    if (eta_max <= 0.0) {
        const double area = 2.0 * (double)L0_UNITS * (double)HEIGHT_UNITS;
        eta_max = (double)n * M_PI * radius_sigma * radius_sigma / area;
    }
    return eta_max;
}

static void live_scene_update_radius_from_eta(void) {
    live_scene_radius_sigma = (float)live_scene_radius_for_eta_value(live_scene_wall_count,
                                                                     live_scene_particle_count,
                                                                     live_scene_eta);
}

static void live_scene_update_eta_from_radius(void) {
    live_scene_eta = (float)live_scene_eta_for_radius_value(live_scene_wall_count,
                                                           live_scene_particle_count,
                                                           live_scene_radius_sigma);
    if (live_scene_eta > 0.75f) {
        live_scene_eta = 0.75f;
        live_scene_update_radius_from_eta();
    }
    if (live_scene_eta < 0.0f) live_scene_eta = 0.0f;
}

static void live_scene_set_mode(LiveSceneMode mode) {
    if (mode < 0 || mode >= LIVE_SCENE_MODE_COUNT) mode = LIVE_SCENE_MODE_ENERGY_TRANSFER;
    LiveSceneMode old_mode = live_scene_mode;
    live_scene_mode = mode;
    live_scene_mode_dropdown_open = 0;

    if (live_scene_mode == LIVE_SCENE_MODE_ENERGY_TRANSFER &&
        old_mode != LIVE_SCENE_MODE_ENERGY_TRANSFER &&
        live_scene_wall_count < 2) {
        live_scene_wall_count = 2;
        live_scene_spring_l0_sigma = (float)L0_UNITS;
    }
    if (live_scene_mode == LIVE_SCENE_MODE_SZILARD) live_scene_engine = LIVE_SCENE_ENGINE_EDMD;
    live_scene_apply_mode_constraints();
    live_scene_update_radius_from_eta();
}

static void live_scene_sync_from_runtime(void) {
    live_scene_mode = live_scene_mode_from_preset();
    if (sim_mode == MODE_EDMD) live_scene_engine = LIVE_SCENE_ENGINE_EDMD;
    else if (sim_mode == MODE_EDMD_HYBRID) live_scene_engine = LIVE_SCENE_ENGINE_HYBRID;
    else live_scene_engine = LIVE_SCENE_ENGINE_TIME_TOI;
    live_scene_wall_count = num_internal_walls;
    if (live_scene_wall_count < 1) live_scene_wall_count = 1;
    if (live_scene_wall_count > LIVE_SCENE_MAX_WALLS) live_scene_wall_count = LIVE_SCENE_MAX_WALLS;
    live_scene_apply_mode_constraints();

    live_scene_particle_count = particles_active;
    if (cli_override_particles > 0) live_scene_particle_count = cli_override_particles;
    if (live_scene_particle_count < 1) live_scene_particle_count = 1;
    if (live_scene_particle_count > live_scene_particle_max()) live_scene_particle_count = live_scene_particle_max();

    const int radius_is_primary_control = (cli_particle_radius_set && cli_particle_radius_sigma > 0.0f);
    if (radius_is_primary_control) {
        live_scene_radius_sigma = cli_particle_radius_sigma;
        live_scene_update_eta_from_radius();
    } else if (cli_packing_fraction_set && cli_packing_fraction > 0.0f) {
        live_scene_eta = cli_packing_fraction;
        if (live_scene_eta > 0.75f) live_scene_eta = 0.75f;
        live_scene_update_radius_from_eta();
    } else {
        // Default for the reset panel: choose the radius that gives eta=0.2.
        // This is more useful for exploration than the compile-time fallback radius.
        (void)live_current_nominal_eta();
        live_scene_eta = 0.20f;
        live_scene_update_radius_from_eta();
    }
    if (live_scene_eta < 0.02f) live_scene_eta = 0.02f;
    if (live_scene_eta > 0.75f) live_scene_eta = 0.75f;

    if (live_scene_mode == LIVE_SCENE_MODE_ENERGY_TRANSFER) {
        if (cli_spring_eq_set) {
            live_scene_spring_l0_sigma = cli_spring_eq_sigma;
        } else if (num_internal_walls > 0 && all_wall_positions) {
            int idx = leftmost_wall_index;
            if (idx < 0 || idx >= num_internal_walls) idx = 0;
            live_scene_spring_l0_sigma = (float)(((double)all_wall_positions[idx] - (double)XW1) /
                                                 (double)PIXELS_PER_SIGMA);
        } else {
            live_scene_spring_l0_sigma = (float)L0_UNITS;
        }
        live_scene_spring_l0_sigma = (float)live_scene_clamped_spring_l0_sigma();
        if (radius_is_primary_control) live_scene_update_eta_from_radius();
        else live_scene_update_radius_from_eta();
    }

    for (int i = 0; i < LIVE_SCENE_MAX_WALLS; ++i) {
        float mf = 200.0f;
        if (cli_wall_masses_set && i < cli_wall_masses_count && cli_wall_masses[i] > 0.0f) {
            mf = cli_wall_masses[i];
        } else if (i < num_internal_walls && all_wall_masses && all_wall_masses[i] > 0.0f) {
            mf = all_wall_masses[i] / (float)PARTICLE_MASS;
        }
        if (mf < 2.0f) mf = 2.0f;
        if (mf > 16000.0f) mf = 16000.0f;
        live_scene_wall_mass[i] = mf;
    }
}

static double live_scene_get_value(int index) {
    if (index == 0) return (double)live_scene_wall_count;
    if (index == 1) return (double)live_scene_particle_count;
    if (index == 2) return (double)live_scene_radius_sigma;
    if (index == 3) return (double)live_scene_eta;
    if (index == 4) return (double)live_scene_spring_l0_sigma;
    if (index >= 5 && index < LIVE_SCENE_ROW_COUNT) return (double)live_scene_wall_mass[index - 5];
    return 0.0;
}

static void live_scene_set_value(int index, double value) {
    if (index == 0) {
        if (live_scene_fixed_wall_count() > 0) return;
        int walls = (int)lround(live_clamp(value, 1.0, (double)LIVE_SCENE_MAX_WALLS));
        if (walls < 1) walls = 1;
        if (walls > LIVE_SCENE_MAX_WALLS) walls = LIVE_SCENE_MAX_WALLS;
        live_scene_wall_count = walls;
        live_scene_update_radius_from_eta();
    } else if (index == 1) {
        int n = (int)lround(live_clamp(value, 1.0, (double)live_scene_particle_max()));
        if (n < 1) n = 1;
        live_scene_particle_count = n;
        live_scene_update_radius_from_eta();
    } else if (index == 2) {
        const double r_max = live_scene_radius_for_eta_value(live_scene_wall_count, live_scene_particle_count, 0.75);
        live_scene_radius_sigma = (float)live_clamp(value, 0.005, fmax(0.006, r_max));
        live_scene_update_eta_from_radius();
    } else if (index == 3) {
        live_scene_eta = (float)live_clamp(value, 0.02, 0.75);
        live_scene_update_radius_from_eta();
    } else if (index == 4) {
        live_scene_spring_l0_sigma = (float)live_clamp(value,
                                                      live_scene_spring_l0_min_sigma(),
                                                      live_scene_spring_l0_max_sigma());
        live_scene_update_radius_from_eta();
    } else if (index >= 5 && index < LIVE_SCENE_ROW_COUNT) {
        live_scene_wall_mass[index - 5] = (float)live_clamp(value, 2.0, 16000.0);
    }
}

static const char *live_scene_label(int index) {
    static char label[16];
    if (index == 0) return "walls";
    if (index == 1) return "particles";
    if (index == 2) return "radius σ";
    if (index == 3) return "eta";
    if (index == 4) return "spring L0";
    if (index >= 5 && index < LIVE_SCENE_ROW_COUNT) {
        snprintf(label, sizeof(label), "M%d", index - 4);
        return label;
    }
    return "?";
}

static int live_scene_log_scale(int index) {
    return (index == 1 || index == 2 || (index >= 5 && index < LIVE_SCENE_ROW_COUNT));
}

static double live_scene_min_value(int index) {
    if (index == 0) return 1.0;
    if (index == 1) return 1.0;
    if (index == 2) return 0.005;
    if (index == 3) return 0.02;
    if (index == 4) return live_scene_spring_l0_min_sigma();
    return 2.0;
}

static double live_scene_max_value(int index) {
    if (index == 0) return (double)LIVE_SCENE_MAX_WALLS;
    if (index == 1) return (double)live_scene_particle_max();
    if (index == 2) {
        const double r_max = live_scene_radius_for_eta_value(live_scene_wall_count, live_scene_particle_count, 0.75);
        return fmax(0.006, r_max);
    }
    if (index == 3) return 0.75;
    if (index == 4) return live_scene_spring_l0_max_sigma();
    return 16000.0;
}

static double live_scene_value_to_norm(int index, double value) {
    const double lo = live_scene_min_value(index);
    const double hi = live_scene_max_value(index);
    value = live_clamp(value, lo, hi);
    if (live_scene_log_scale(index)) {
        return live_clamp((log(value) - log(lo)) / (log(hi) - log(lo)), 0.0, 1.0);
    }
    return live_clamp((value - lo) / (hi - lo), 0.0, 1.0);
}

static double live_scene_norm_to_value(int index, double norm) {
    const double lo = live_scene_min_value(index);
    const double hi = live_scene_max_value(index);
    norm = live_clamp(norm, 0.0, 1.0);
    if (live_scene_log_scale(index)) return exp(log(lo) + norm * (log(hi) - log(lo)));
    return lo + norm * (hi - lo);
}

static void live_scene_layout(int *x, int *y, int *w, int *h) {
    const int panel_w = 420;
    const int visible_rows = live_scene_visible_param_count();
    const int panel_h = 134 + visible_rows * 30 + 72;
    int panel_x = MAX_X - panel_w - 12;
    if (panel_x < XW1 + 12) panel_x = XW1 + 12;
    int panel_y = YW1 + 12;
    if (live_controls_visible) {
        int cx = 0, cy = 0, cw = 0, ch = 0;
        live_controls_layout(&cx, &cy, &cw, &ch);
        panel_y = cy + ch + 10;
    }
    if (panel_y < YW1 + 10) panel_y = YW1 + 10;
    if (panel_y + panel_h > MAX_Y - 6) panel_y = MAX_Y - panel_h - 6;
    if (panel_y < 6) panel_y = 6;
    if (x) *x = panel_x;
    if (y) *y = panel_y;
    if (w) *w = panel_w;
    if (h) *h = panel_h;
}

static SDL_Rect live_scene_mode_rect(int panel_x, int panel_y, int panel_w) {
    SDL_Rect r;
    r.x = panel_x + 116;
    r.y = panel_y + 62;
    r.w = panel_w - 206;
    r.h = 22;
    return r;
}

static SDL_Rect live_scene_mode_option_rect(int panel_x, int panel_y, int panel_w, int option_index) {
    SDL_Rect r = live_scene_mode_rect(panel_x, panel_y, panel_w);
    r.y += 24 + option_index * 24;
    r.h = 22;
    return r;
}

static SDL_Rect live_scene_engine_rect(int panel_x, int panel_y, int panel_w) {
    SDL_Rect r;
    r.x = panel_x + 116;
    r.y = panel_y + 92;
    r.w = panel_w - 206;
    r.h = 22;
    return r;
}

static SDL_Rect live_scene_track_rect(int panel_x, int panel_y, int panel_w, int visible_row_index) {
    SDL_Rect r;
    r.x = panel_x + 116;
    r.y = panel_y + 130 + visible_row_index * 30;
    r.w = panel_w - 206;
    r.h = 8;
    return r;
}

static SDL_Rect live_scene_apply_rect(int panel_x, int panel_y, int panel_w) {
    SDL_Rect r;
    r.x = panel_x + 14;
    r.y = panel_y + 134 + live_scene_visible_param_count() * 30;
    r.w = panel_w - 28;
    r.h = 28;
    return r;
}

static void live_scene_set_from_mouse(int index, int visible_row, int mouse_x) {
    int px = 0, py = 0, pw = 0, ph = 0;
    live_scene_layout(&px, &py, &pw, &ph);
    SDL_Rect track = live_scene_track_rect(px, py, pw, visible_row);
    double norm = (track.w > 0) ? ((double)(mouse_x - track.x) / (double)track.w) : 0.0;
    live_scene_set_value(index, live_scene_norm_to_value(index, norm));
}

static void live_scene_prepare_rebuild(void) {
    live_scene_apply_mode_constraints();
    if (live_scene_mode == LIVE_SCENE_MODE_SZILARD) live_scene_engine = LIVE_SCENE_ENGINE_EDMD;
    if (live_scene_engine == LIVE_SCENE_ENGINE_EDMD) {
        sim_mode = MODE_EDMD;
    } else if (live_scene_engine == LIVE_SCENE_ENGINE_HYBRID) {
        sim_mode = MODE_EDMD_HYBRID;
    } else {
        sim_mode = MODE_TIME;
    }
    if (sim_mode == MODE_TIME) {
        // Realtime high-N TIME mode is only useful if it uses TOI + adaptive
        // substeps. Otherwise the first pressure wave creates deep overlaps
        // that later look like nonphysical clamped particle sheets.
        cli_pp_toi = 1;
        if (cli_pp_toi_max_events < 8) cli_pp_toi_max_events = 8;
        if (cli_pp_toi_neighbor_ring < 1) cli_pp_toi_neighbor_ring = 1;
        cli_substeps_override = -1;
        cli_substeps_auto = true;
        if (cli_substeps_auto_min < 4) cli_substeps_auto_min = 4;
        if (cli_substeps_auto_max < 240) cli_substeps_auto_max = 240;
        if (cli_substeps_auto_dx_frac > 0.12f) cli_substeps_auto_dx_frac = 0.12f;
    }

    int walls = live_scene_wall_count;
    if (walls < 1) walls = 1;
    if (walls > LIVE_SCENE_MAX_WALLS) walls = LIVE_SCENE_MAX_WALLS;

    cli_left_empty = false;
    cli_enable_energy_measurement = false;
    cli_auto_piston_step = false;
    cli_eff_output_mode = EFF_OUTPUT_SPRING;
    preset_custom_positions = true;
    preset_custom_fractions = false;
    preset_custom_absolute_counts = false;
    memset(preset_wall_fraction, 0, sizeof(preset_wall_fraction));
    memset(preset_segment_fraction, 0, sizeof(preset_segment_fraction));

    if (live_scene_mode == LIVE_SCENE_MODE_ENERGY_TRANSFER) {
        cli_experiment_preset = EXPERIMENT_PRESET_ENERGY_TRANSFER;
        enable_speed_of_sound_experiments = true;
        cli_left_empty = true;
        cli_enable_energy_measurement = true;
        cli_eff_output_mode = EFF_OUTPUT_SPRING;
        live_scene_spring_l0_sigma = (float)live_scene_clamped_spring_l0_sigma();
        live_apply_spring_enabled(1);
    } else if (live_scene_mode == LIVE_SCENE_MODE_SZILARD) {
        sim_mode = MODE_EDMD;
        cli_experiment_preset = EXPERIMENT_PRESET_SZILARD_ENGINE;
        enable_speed_of_sound_experiments = true;
        walls = 1;
        live_scene_wall_count = 1;
        cli_seeding_mode = SEEDING_RANDOM;
        cli_szilard_perm_swap = 1;
        live_apply_spring_enabled(0);
    } else if (live_scene_mode == LIVE_SCENE_MODE_SPEED_OF_SOUND) {
        cli_experiment_preset = EXPERIMENT_PRESET_SPEED_OF_SOUND;
        enable_speed_of_sound_experiments = true;
        walls = 1;
        live_scene_wall_count = 1;
        cli_seeding_mode = SEEDING_GRID;
        live_apply_spring_enabled(0);
    }

    int n_particles = live_scene_particle_count;
    if (n_particles < 1) n_particles = 1;
    if (n_particles > live_scene_particle_max()) n_particles = live_scene_particle_max();

    float radius_sigma = live_scene_radius_sigma;
    const float radius_max = (float)live_scene_radius_for_eta_value(walls, n_particles, 0.75);
    if (radius_sigma < 0.005f) radius_sigma = 0.005f;
    if (radius_sigma > radius_max) radius_sigma = radius_max;
    float eta = (float)live_scene_eta_for_radius_value(walls, n_particles, radius_sigma);
    if (eta < 0.02f) eta = 0.02f;
    if (eta > 0.75f) {
        eta = 0.75f;
        radius_sigma = (float)live_scene_radius_for_eta_value(walls, n_particles, eta);
    }
    live_scene_particle_count = n_particles;
    live_scene_radius_sigma = radius_sigma;
    live_scene_eta = eta;

    // Keep the live view honest for high-N eta checks. A render-only radius
    // larger than the collision radius makes eta=0.2 look like a much denser
    // gas. Still allow at least a subpixel/1px mark so tiny particles remain visible.
    const float real_radius_px = radius_sigma * (float)PIXELS_PER_SIGMA;
    if (gui_particle_draw_min_px > real_radius_px) {
        gui_particle_draw_min_px = fmaxf(0.5f, real_radius_px);
    }

    cli_num_walls = walls;
    cli_requested_walls = walls;
    num_internal_walls = walls;
    if (live_scene_mode == LIVE_SCENE_MODE_ENERGY_TRANSFER) {
        cli_spring_eq_set = 1;
        cli_spring_eq_sigma = live_scene_spring_l0_sigma;
    } else {
        cli_spring_eq_set = 0;
        cli_spring_eq_sigma = 0.0f;
    }
    cli_override_particles = n_particles;
    particles_active = n_particles;
    preset_custom_positions = true;
    preset_wall_count = walls;
    memset(preset_wall_fraction, 0, sizeof(preset_wall_fraction));
    memset(cli_wall_positions, 0, sizeof(cli_wall_positions));
    for (int w = 0; w < walls; ++w) {
        const float pos_sigma = (float)live_scene_wall_position_sigma(walls, w);
        const float frac = pos_sigma / (2.0f * L0_UNITS);
        preset_wall_fraction[w] = frac;
        cli_wall_positions[w] = pos_sigma;
    }

    cli_wall_masses_set = 1;
    cli_wall_masses_count = walls;
    for (int w = 0; w < 5; ++w) {
        float mf = live_scene_wall_mass[w];
        if (mf < 2.0f) mf = 2.0f;
        if (mf > 16000.0f) mf = 16000.0f;
        cli_wall_masses[w] = mf;
    }

    // Apply the real collision radius directly. The eta slider only computes this radius.
    cli_packing_fraction_set = false;
    cli_packing_fraction = eta;
    cli_particle_radius_set = 1;
    cli_particle_radius_sigma = radius_sigma;
    cli_particle_diameter_set = 0;
    cli_particle_diameter_sigma = 0.0f;
    particle_scale_runtime = radius_sigma / PARTICLE_RADIUS_UNIT;
    if (cli_override_wall_thickness_sigma <= 0.0f) {
        wall_thickness_runtime = 0.0f; // recompute from the new real collision radius.
    }

    const int segments = walls + 1;
    free(cli_particles_box_counts);
    cli_particles_box_counts = (int *)malloc((size_t)segments * sizeof(int));
    if (!cli_particles_box_counts) {
        fprintf(stderr, "Failed to allocate scene particle-box counts.\n");
        exit(EXIT_FAILURE);
    }
    cli_particles_box_counts_count = (size_t)segments;
    preset_custom_absolute_counts = true;
    preset_custom_fractions = false;
    cli_particles_box_left = -1;
    cli_particles_box_right = -1;

    live_scene_counts_for(walls, n_particles, cli_particles_box_counts, segments);

    live_scene_rebuild_requested = 1;
    if (!cli_quiet) {
        printf("Scene reset queued: mode=%s engine=%s walls=%d springL0=%.3fσ N=%d r=%.4fσ eta=%.3f masses=",
               live_scene_mode_label(live_scene_mode), live_scene_engine_label(),
               walls, (double)live_scene_spring_l0_sigma,
               n_particles, (double)radius_sigma, (double)eta);
        for (int w = 0; w < walls; ++w) printf("%s%.3g", w ? "," : "", (double)cli_wall_masses[w]);
        printf("\n");
    }
}

static void live_rebuild_edmd_from_current_state(void) {
    if (!(sim_mode == MODE_EDMD || sim_mode == MODE_EDMD_HYBRID)) {
        if (g_edmd) {
            edmd_destroy(g_edmd);
            g_edmd = NULL;
        }
        return;
    }
    if (g_edmd) {
        edmd_destroy(g_edmd);
        g_edmd = NULL;
    }

    EDMD_Params prm = {0};
    prm.boxW = (double)(XW2 - XW1);
    prm.boxH = (double)(YW2 - YW1);
    prm.radius = (double)PARTICLE_RADIUS;
    prm.N = particles_active;
    prm.cell_size = 0.0;

    int dcount = (num_internal_walls > 0) ? num_internal_walls : 0;
    if (dcount > EDMD_MAX_DIVIDERS) dcount = EDMD_MAX_DIVIDERS;
    prm.divider_count = dcount;
    int primary_idx = primary_wall_index;
    if (primary_idx < 0 || primary_idx >= dcount) primary_idx = 0;
    int extra_idx = 0;
    for (int w = 0; w < dcount; ++w) {
        double pos_px = all_wall_positions ? (double)all_wall_positions[w] : (double)wall_x;
        prm.divider_x[w] = pos_px - (double)XW1;
        prm.divider_thickness[w] = (double)WALL_THICKNESS;
        double mf = (double)(WALL_MASS / PARTICLE_MASS);
        if (cli_wall_masses_set && w < cli_wall_masses_count && cli_wall_masses[w] > 0.0f) {
            mf = (double)cli_wall_masses[w];
        }
        prm.divider_mass[w] = mf;
        prm.divider_vx[w] = 0.0;
        if (w == primary_idx) {
            prm.divider_vx[w] = (double)vx_wall;
        } else if (extra_wall_velocity && extra_idx < extra_wall_count) {
            prm.divider_vx[w] = (double)extra_wall_velocity[extra_idx++];
        }
    }

    prm.heatbath_enabled = heatbath_enabled;
    prm.heatbath_temperature = (double)heatbath_temperature;
    prm.thermal_wall_mode = thermal_wall_mode;
    prm.mb_overshoot_factor = (double)mb_overshoot_factor;
    prm.stability_window_percent = (double)stability_window_percent;
    prm.particle_mass = (double)PARTICLE_MASS;
    prm.kB = (double)kB_effective();
    prm.pp_collisions_enabled = 1;

    if (energy_measurement.spring_enabled && prm.divider_count > 0) {
        int pidx = primary_idx;
        const double xeq_px = (energy_measurement.equilibrium_position != 0.0f)
                                ? (double)energy_measurement.equilibrium_position
                                : (double)wall_x;
        prm.divider_k[pidx] = (double)energy_measurement.spring_constant;
        prm.divider_xeq[pidx] = xeq_px - (double)XW1;
    }

    g_edmd = edmd_create(&prm);
    if (!g_edmd) {
        fprintf(stderr, "EDMD create failed during live scene rebuild.\n");
        exit(EXIT_FAILURE);
    }
    edmd_config_pistons(g_edmd,
                        1, (double)((piston_left_x + 5.0f) - XW1), (double)vx_piston_left, 0.0,
                        1, (double)(piston_right_x - XW1), (double)vx_piston_right, 0.0);

    EDMD_Particle *P = (EDMD_Particle *)edmd_particles(g_edmd);
    for (int i = 0; i < particles_active; ++i) {
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

    if (sim_mode == MODE_EDMD) edmd_reschedule_all(g_edmd);
    else edmd_reschedule_all_pp_only(g_edmd);
    edmd_reset_work(g_edmd);
}

static double live_slider_value_to_norm(const LiveSlider *slider, double value) {
    if (!slider) return 0.0;
    value = live_clamp(value, slider->min_value, slider->max_value);
    if (slider->logarithmic) {
        double lo = fmax(slider->min_value, 1e-12);
        double hi = fmax(slider->max_value, lo * 1.000001);
        return live_clamp((log(value) - log(lo)) / (log(hi) - log(lo)), 0.0, 1.0);
    }
    double span = slider->max_value - slider->min_value;
    if (fabs(span) < 1e-12) return 0.0;
    return live_clamp((value - slider->min_value) / span, 0.0, 1.0);
}

static double live_slider_norm_to_value(const LiveSlider *slider, double norm) {
    if (!slider) return 0.0;
    norm = live_clamp(norm, 0.0, 1.0);
    if (slider->logarithmic) {
        double lo = fmax(slider->min_value, 1e-12);
        double hi = fmax(slider->max_value, lo * 1.000001);
        return exp(log(lo) + norm * (log(hi) - log(lo)));
    }
    return slider->min_value + norm * (slider->max_value - slider->min_value);
}

static void live_controls_layout(int *x, int *y, int *w, int *h) {
    const int panel_w = 420;
    const int panel_h = 58 + live_slider_count * 34 + 136;
    int panel_x = MAX_X - panel_w - 12;
    if (panel_x < XW1 + 12) panel_x = XW1 + 12;
    int panel_y = YW1 + 12;
    if (panel_y < YW1 + 10) panel_y = YW1 + 10;
    if (panel_y + panel_h > MAX_Y - 6) panel_y = MAX_Y - panel_h - 6;
    if (panel_y < 6) panel_y = 6;
    if (x) *x = panel_x;
    if (y) *y = panel_y;
    if (w) *w = panel_w;
    if (h) *h = panel_h;
}

static SDL_Rect live_slider_track_rect(int panel_x, int panel_y, int panel_w, int slider_index) {
    SDL_Rect r;
    r.x = panel_x + 116;
    r.y = panel_y + 49 + slider_index * 34;
    r.w = panel_w - 206;
    r.h = 8;
    return r;
}

static SDL_Rect live_spring_toggle_rect(int panel_x, int panel_y, int panel_w) {
    SDL_Rect r;
    r.x = panel_x + panel_w - 92;
    r.y = panel_y + 45 + live_slider_count * 34;
    r.w = 74;
    r.h = 22;
    return r;
}

static SDL_Rect live_mb_toggle_rect(int panel_x, int panel_y, int panel_w) {
    SDL_Rect r;
    r.x = panel_x + panel_w - 92;
    r.y = panel_y + 45 + (live_slider_count + 1) * 34;
    r.w = 74;
    r.h = 22;
    return r;
}

static SDL_Rect live_mb_scale_rect(int panel_x, int panel_y, int panel_w) {
    SDL_Rect r;
    r.x = panel_x + panel_w - 92;
    r.y = panel_y + 45 + (live_slider_count + 2) * 34;
    r.w = 74;
    r.h = 22;
    return r;
}

static SDL_Rect live_density_curve_toggle_rect(int panel_x, int panel_y, int panel_w) {
    SDL_Rect r;
    r.x = panel_x + panel_w - 92;
    r.y = panel_y + 45 + (live_slider_count + 3) * 34;
    r.w = 74;
    r.h = 22;
    return r;
}

static int live_point_in_rect(int x, int y, SDL_Rect r) {
    return (x >= r.x && x < r.x + r.w && y >= r.y && y < r.y + r.h);
}

static void live_controls_set_from_mouse(int slider_index, int mouse_x) {
    int px = 0, py = 0, pw = 0, ph = 0;
    live_controls_layout(&px, &py, &pw, &ph);
    SDL_Rect track = live_slider_track_rect(px, py, pw, slider_index);
    double norm = (track.w > 0) ? ((double)(mouse_x - track.x) / (double)track.w) : 0.0;
    double value = live_slider_norm_to_value(&live_sliders[slider_index], norm);
    live_slider_set_value(slider_index, value);
}

static int live_controls_runtime_panel_handle_event(const SDL_Event *e) {
    if (!e || !live_controls_visible) return 0;

    if (e->type == SDL_MOUSEBUTTONDOWN &&
        (e->button.button == SDL_BUTTON_LEFT || e->button.button == SDL_BUTTON_RIGHT)) {
        int mx = 0, my = 0;
        gui_mouse_to_canvas_xy(e->button.x, e->button.y, &mx, &my);
        int px = 0, py = 0, pw = 0, ph = 0;
        live_controls_layout(&px, &py, &pw, &ph);
        SDL_Rect panel = { px, py, pw, ph };
        if (!live_point_in_rect(mx, my, panel)) return 0;

        for (int i = 0; i < live_slider_count; ++i) {
            SDL_Rect row = { px + 8, py + 38 + i * 34, pw - 16, 30 };
            if (live_point_in_rect(mx, my, row)) {
                if (e->button.button == SDL_BUTTON_RIGHT) {
                    live_controls_drag_index = -1;
                    live_slider_set_value(i, live_sliders[i].default_value);
                } else {
                    live_controls_drag_index = i;
                    live_controls_set_from_mouse(i, mx);
                }
                return 1;
            }
        }

        SDL_Rect spring_row = { px + 8, py + 38 + live_slider_count * 34, pw - 16, 30 };
        if (live_point_in_rect(mx, my, spring_row)) {
            live_controls_drag_index = -1;
            if (live_scene_spring_available()) {
                live_apply_spring_enabled(!energy_measurement.spring_enabled);
            } else if (!cli_quiet) {
                printf("Live spring force: unavailable for %s\n", live_scene_mode_label(live_scene_mode));
            }
            return 1;
        }

        SDL_Rect mb_row = { px + 8, py + 38 + (live_slider_count + 1) * 34, pw - 16, 30 };
        if (live_point_in_rect(mx, my, mb_row)) {
            live_controls_drag_index = -1;
            gui_mb_speed_view = !gui_mb_speed_view;
            if (gui_mb_speed_view) {
                if (gui_dist_mode == 0) gui_dist_mode = 3;
            }
            if (!cli_quiet) {
                printf("🌈 MB speed view: %s\n", gui_mb_speed_view ? "ON" : "off");
            }
            return 1;
        }

        SDL_Rect mb_scale_row = { px + 8, py + 38 + (live_slider_count + 2) * 34, pw - 16, 30 };
        if (live_point_in_rect(mx, my, mb_scale_row)) {
            live_controls_drag_index = -1;
            gui_mb_capture_speed_scale_now();
            if (!cli_quiet) {
                printf("🌈 MB speed scale set: vmax=%.3f\n", (double)gui_mb_speed_vmax_ref);
            }
            return 1;
        }

        SDL_Rect den_curve_row = { px + 8, py + 38 + (live_slider_count + 3) * 34, pw - 16, 30 };
        if (live_point_in_rect(mx, my, den_curve_row)) {
            live_controls_drag_index = -1;
            gui_density_curve_overlay = !gui_density_curve_overlay;
            if (gui_density_curve_overlay && gui_dist_mode == 0) gui_dist_mode = 3;
            if (!cli_quiet) {
                printf("📈 Global x-density curve: %s\n", gui_density_curve_overlay ? "ON" : "off");
            }
            return 1;
        }

        return 1;
    }

    if (e->type == SDL_MOUSEBUTTONUP && e->button.button == SDL_BUTTON_LEFT) {
        if (live_controls_drag_index >= 0) {
            live_controls_drag_index = -1;
            return 1;
        }
    }

    if (e->type == SDL_MOUSEMOTION && live_controls_drag_index >= 0) {
        int mx = 0, my = 0;
        gui_mouse_to_canvas_xy(e->motion.x, e->motion.y, &mx, &my);
        live_controls_set_from_mouse(live_controls_drag_index, mx);
        return 1;
    }

    return 0;
}

static int live_controls_handle_event(const SDL_Event *e) {
    if (!e) return 0;

    if (e->type == SDL_KEYDOWN && e->key.keysym.sym == SDLK_F2) {
        live_controls_visible = !live_controls_visible;
        live_controls_drag_index = -1;
        if (!cli_quiet) {
            printf("Live SDL sliders: %s (F2)\n", live_controls_visible ? "ON" : "off");
        }
        return 1;
    }

    if (e->type == SDL_KEYDOWN && e->key.keysym.sym == SDLK_F3) {
        live_scene_visible = !live_scene_visible;
        live_scene_drag_index = -1;
        live_scene_drag_visible_row = -1;
        live_scene_mode_dropdown_open = 0;
        if (live_scene_visible) live_scene_sync_from_runtime();
        if (!cli_quiet) {
            printf("Live scene reset panel: %s (F3)\n", live_scene_visible ? "ON" : "off");
        }
        return 1;
    }

    // Give the blue F2 runtime panel first click priority. For short simulation
    // boxes the orange F3 panel can be clamped upward and overlap it.
    if (live_controls_runtime_panel_handle_event(e)) return 1;

    if (live_scene_visible) {
        if (e->type == SDL_MOUSEBUTTONDOWN &&
            (e->button.button == SDL_BUTTON_LEFT || e->button.button == SDL_BUTTON_RIGHT)) {
            int mx = 0, my = 0;
            gui_mouse_to_canvas_xy(e->button.x, e->button.y, &mx, &my);
            int sx = 0, sy = 0, sw = 0, sh = 0;
            live_scene_layout(&sx, &sy, &sw, &sh);
            SDL_Rect panel = { sx, sy, sw, sh };
            if (live_point_in_rect(mx, my, panel)) {
                SDL_Rect mode_row = { sx + 8, sy + 58, sw - 16, 27 };
                if (live_scene_mode_dropdown_open) {
                    for (int opt = 0; opt < LIVE_SCENE_MODE_COUNT; ++opt) {
                        SDL_Rect opt_rect = live_scene_mode_option_rect(sx, sy, sw, opt);
                        if (live_point_in_rect(mx, my, opt_rect)) {
                            live_scene_drag_index = -1;
                            live_scene_drag_visible_row = -1;
                            live_scene_set_mode((LiveSceneMode)opt);
                            return 1;
                        }
                    }
                    if (!live_point_in_rect(mx, my, mode_row)) {
                        live_scene_mode_dropdown_open = 0;
                    }
                }
                if (live_point_in_rect(mx, my, mode_row)) {
                    live_scene_drag_index = -1;
                    live_scene_drag_visible_row = -1;
                    if (e->button.button == SDL_BUTTON_RIGHT) {
                        int next = ((int)live_scene_mode + 1) % LIVE_SCENE_MODE_COUNT;
                        live_scene_set_mode((LiveSceneMode)next);
                    } else {
                        live_scene_mode_dropdown_open = !live_scene_mode_dropdown_open;
                    }
                    return 1;
                }
                SDL_Rect engine_row = { sx + 8, sy + 88, sw - 16, 27 };
                if (live_point_in_rect(mx, my, engine_row)) {
                    live_scene_drag_index = -1;
                    live_scene_drag_visible_row = -1;
                    live_scene_mode_dropdown_open = 0;
                    if (live_scene_mode == LIVE_SCENE_MODE_SZILARD) {
                        live_scene_engine = LIVE_SCENE_ENGINE_EDMD;
                        if (!cli_quiet) printf("Scene engine: Szilard uses EDMD for semipermeable gates.\n");
                    } else {
                        live_scene_cycle_engine(e->button.button == SDL_BUTTON_RIGHT ? -1 : 1);
                        if (!cli_quiet) printf("Scene engine queued: %s\n", live_scene_engine_label());
                    }
                    return 1;
                }
                SDL_Rect apply = live_scene_apply_rect(sx, sy, sw);
                if (live_point_in_rect(mx, my, apply)) {
                    live_scene_drag_index = -1;
                    live_scene_drag_visible_row = -1;
                    live_scene_mode_dropdown_open = 0;
                    live_scene_prepare_rebuild();
                    return 1;
                }
                const int visible_rows = live_scene_visible_param_count();
                for (int vr = 0; vr < visible_rows; ++vr) {
                    const int i = live_scene_index_for_visible_row(vr);
                    if (i < 0) continue;
                    SDL_Rect row = { sx + 8, sy + 118 + vr * 30, sw - 16, 27 };
                    if (live_point_in_rect(mx, my, row)) {
                        if (e->button.button == SDL_BUTTON_RIGHT) {
                            if (i == 0) live_scene_set_value(i, 2.0);
                            else if (i == 1) live_scene_set_value(i, (double)DEFAULT_PARTICLES_ACTIVE);
                            else if (i == 2) { live_scene_set_value(3, 0.20); }
                            else if (i == 3) live_scene_set_value(i, 0.20);
                            else if (i == 4) live_scene_set_value(i, (double)L0_UNITS);
                            else live_scene_set_value(i, 200.0);
                            live_scene_drag_index = -1;
                            live_scene_drag_visible_row = -1;
                        } else {
                            live_scene_drag_index = i;
                            live_scene_drag_visible_row = vr;
                            live_scene_set_from_mouse(i, vr, mx);
                        }
                        return 1;
                    }
                }
                return 1;
            }
        }
        if (e->type == SDL_MOUSEBUTTONUP && e->button.button == SDL_BUTTON_LEFT) {
            if (live_scene_drag_index >= 0) {
                live_scene_drag_index = -1;
                live_scene_drag_visible_row = -1;
                return 1;
            }
        }
        if (e->type == SDL_MOUSEMOTION && live_scene_drag_index >= 0) {
            int mx = 0, my = 0;
            gui_mouse_to_canvas_xy(e->motion.x, e->motion.y, &mx, &my);
            live_scene_set_from_mouse(live_scene_drag_index, live_scene_drag_visible_row, mx);
            return 1;
        }
    }

    return 0;
}

static void render_live_controls(SDL_Renderer *renderer, TTF_Font *font) {
    if (!renderer || !font) return;

    SDL_Color white = {245, 245, 245, 255};
    SDL_Color muted = {180, 190, 190, 255};
    SDL_Color cyan = {0, 210, 255, 255};
    SDL_Color green = {80, 230, 150, 255};

    if (!live_controls_visible) {
        draw_text(renderer, font, "F2: live sliders", XW1 + 10, YW2 - 28, muted);
    } else {
        int px = 0, py = 0, pw = 0, ph = 0;
        live_controls_layout(&px, &py, &pw, &ph);

        SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_BLEND);
        SDL_Rect panel = { px, py, pw, ph };
        SDL_SetRenderDrawColor(renderer, 6, 12, 14, 218);
        SDL_RenderFillRect(renderer, &panel);
        SDL_SetRenderDrawColor(renderer, 0, 210, 255, 180);
        SDL_RenderDrawRect(renderer, &panel);

        draw_text(renderer, font, "Live SDL controls (F2 hide)", px + 10, py + 8, cyan);
        draw_text(renderer, font, "drag: runtime only, right-click reset", px + 10, py + 28, muted);

        for (int i = 0; i < live_slider_count; ++i) {
            const LiveSlider *slider = &live_sliders[i];
            const int row_y = py + 43 + i * 34;
            SDL_Rect track = live_slider_track_rect(px, py, pw, i);
            double value = live_slider_get_value(slider);
            double norm = live_slider_value_to_norm(slider, value);
            int knob_x = track.x + (int)lround(norm * (double)track.w);

            char value_text[96];
            if (fabs(value) >= 1000.0 || fabs(value) < 0.01) {
                snprintf(value_text, sizeof(value_text), "%.3g%s", value, slider->unit ? slider->unit : "");
            } else {
                snprintf(value_text, sizeof(value_text), "%.3f%s", value, slider->unit ? slider->unit : "");
            }

            draw_text_clipped(renderer, font, slider->label, px + 10, row_y, white, 96);
            draw_text_clipped(renderer, font, value_text, px + pw - 78, row_y, white, 64);

            SDL_SetRenderDrawColor(renderer, 70, 80, 84, 230);
            SDL_RenderFillRect(renderer, &track);

            SDL_Rect fill = track;
            fill.w = knob_x - track.x;
            if (fill.w < 0) fill.w = 0;
            if (fill.w > track.w) fill.w = track.w;
            SDL_SetRenderDrawColor(renderer, 80, 230, 150, 230);
            SDL_RenderFillRect(renderer, &fill);

            SDL_Rect knob = { knob_x - 5, track.y - 5, 10, 18 };
            SDL_SetRenderDrawColor(renderer, 245, 245, 245, 245);
            SDL_RenderFillRect(renderer, &knob);

            if (i == live_controls_drag_index) {
                SDL_SetRenderDrawColor(renderer, green.r, green.g, green.b, 255);
                SDL_RenderDrawRect(renderer, &knob);
            }
        }

        {
            const int row_y = py + 43 + live_slider_count * 34;
            SDL_Rect sw = live_spring_toggle_rect(px, py, pw);
            const int spring_on = energy_measurement.spring_enabled ? 1 : 0;
            const int spring_available = live_scene_spring_available();
            draw_text_clipped(renderer, font,
                              spring_available ? "spring force" : "spring force n/a",
                              px + 10, row_y, spring_available ? white : muted, 150);
            SDL_SetRenderDrawColor(renderer,
                                   spring_available ? (spring_on ? 60 : 78) : 55,
                                   spring_available ? (spring_on ? 190 : 80) : 55,
                                   spring_available ? (spring_on ? 110 : 84) : 55,
                                   spring_available ? 235 : 150);
            SDL_RenderFillRect(renderer, &sw);
            SDL_SetRenderDrawColor(renderer, 235, 235, 235, 240);
            SDL_RenderDrawRect(renderer, &sw);
            draw_text(renderer, font,
                      spring_available ? (spring_on ? "ON" : "off") : "n/a",
                      sw.x + 19, sw.y + 2, spring_available ? white : muted);
        }

        {
            const int row_y = py + 43 + (live_slider_count + 1) * 34;
            SDL_Rect mb = live_mb_toggle_rect(px, py, pw);
            const int mb_on = gui_mb_speed_view ? 1 : 0;
            draw_text_clipped(renderer, font, "MB speed view", px + 10, row_y, white, 150);
            SDL_SetRenderDrawColor(renderer,
                                   mb_on ? 35 : 78,
                                   mb_on ? 145 : 80,
                                   mb_on ? 220 : 84,
                                   235);
            SDL_RenderFillRect(renderer, &mb);
            SDL_SetRenderDrawColor(renderer, 235, 235, 235, 240);
            SDL_RenderDrawRect(renderer, &mb);
            draw_text(renderer, font, mb_on ? "ON" : "off", mb.x + 21, mb.y + 2, white);
        }

        {
            const int row_y = py + 43 + (live_slider_count + 2) * 34;
            SDL_Rect btn = live_mb_scale_rect(px, py, pw);
            char scale_label[96];
            if (gui_mb_speed_vmax_ref > 1e-6f) {
                snprintf(scale_label, sizeof(scale_label), "MB scale  vmax=%.2f", (double)gui_mb_speed_vmax_ref);
            } else {
                snprintf(scale_label, sizeof(scale_label), "MB scale thermal %.2f",
                         (double)gui_mb_speed_vmax_for_render());
            }
            draw_text_clipped(renderer, font, scale_label, px + 10, row_y, white, 210);
            SDL_SetRenderDrawColor(renderer, 80, 90, 120, 235);
            SDL_RenderFillRect(renderer, &btn);
            SDL_SetRenderDrawColor(renderer, 235, 235, 235, 240);
            SDL_RenderDrawRect(renderer, &btn);
            draw_text(renderer, font, "SET", btn.x + 17, btn.y + 2, white);
        }

        {
            const int row_y = py + 43 + (live_slider_count + 3) * 34;
            SDL_Rect dc = live_density_curve_toggle_rect(px, py, pw);
            const int curve_on = gui_density_curve_overlay ? 1 : 0;
            draw_text_clipped(renderer, font, "x-density curve", px + 10, row_y, white, 150);
            SDL_SetRenderDrawColor(renderer,
                                   curve_on ? 190 : 78,
                                   curve_on ? 100 : 80,
                                   curve_on ? 50 : 84,
                                   235);
            SDL_RenderFillRect(renderer, &dc);
            SDL_SetRenderDrawColor(renderer, 235, 235, 235, 240);
            SDL_RenderDrawRect(renderer, &dc);
            draw_text(renderer, font, curve_on ? "ON" : "off", dc.x + 21, dc.y + 2, white);
        }
    }

    if (!live_scene_visible) {
        draw_text(renderer, font, "F3: scene reset", XW1 + 10, YW2 - 50, muted);
        return;
    }

    int sx = 0, sy = 0, sw = 0, sh = 0;
    live_scene_layout(&sx, &sy, &sw, &sh);
    SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_BLEND);
    SDL_Rect scene_panel = { sx, sy, sw, sh };
    SDL_SetRenderDrawColor(renderer, 18, 12, 6, 224);
    SDL_RenderFillRect(renderer, &scene_panel);
    SDL_SetRenderDrawColor(renderer, 255, 150, 55, 190);
    SDL_RenderDrawRect(renderer, &scene_panel);

    draw_text(renderer, font, "Scene reset controls (F3 hide)", sx + 10, sy + 8, (SDL_Color){255, 175, 80, 255});
    draw_text(renderer, font, "Apply restarts physics; eta capped at 0.75", sx + 10, sy + 28, muted);
    draw_text(renderer, font, "N/radius/eta are coupled real physics", sx + 10, sy + 48, muted);

    {
        const int row_y = sy + 64;
        SDL_Rect mode_box = live_scene_mode_rect(sx, sy, sw);
        SDL_SetRenderDrawColor(renderer, 82, 48, 24, 235);
        SDL_RenderFillRect(renderer, &mode_box);
        SDL_SetRenderDrawColor(renderer, 255, 180, 90, 245);
        SDL_RenderDrawRect(renderer, &mode_box);
        draw_text_clipped(renderer, font, "mode", sx + 10, row_y, white, 96);
        draw_text_clipped(renderer, font, live_scene_mode_label(live_scene_mode), mode_box.x + 6, row_y, white, mode_box.w - 26);
        draw_text(renderer, font, live_scene_mode_dropdown_open ? "^" : "v", mode_box.x + mode_box.w - 18, row_y, white);
    }

    {
        const int row_y = sy + 94;
        SDL_Rect engine_box = live_scene_engine_rect(sx, sy, sw);
        const int locked = (live_scene_mode == LIVE_SCENE_MODE_SZILARD);
        const int is_edmd = live_scene_engine == LIVE_SCENE_ENGINE_EDMD;
        const int is_hybrid = live_scene_engine == LIVE_SCENE_ENGINE_HYBRID;
        SDL_SetRenderDrawColor(renderer,
                               is_edmd ? 54 : (is_hybrid ? 54 : 88),
                               is_edmd ? 78 : (is_hybrid ? 86 : 54),
                               is_edmd ? 126 : (is_hybrid ? 54 : 34),
                               235);
        SDL_RenderFillRect(renderer, &engine_box);
        SDL_SetRenderDrawColor(renderer, locked ? 170 : 255, locked ? 170 : 180, locked ? 170 : 90, 245);
        SDL_RenderDrawRect(renderer, &engine_box);
        draw_text_clipped(renderer, font, "engine", sx + 10, row_y, locked ? muted : white, 96);
        draw_text_clipped(renderer, font,
                          locked ? "EDMD (required)" : live_scene_engine_label(),
                          engine_box.x + 6, row_y, white, engine_box.w - 12);
    }

    const int visible_rows = live_scene_visible_param_count();
    for (int vr = 0; vr < visible_rows; ++vr) {
        const int i = live_scene_index_for_visible_row(vr);
        if (i < 0) continue;
        SDL_Color row_color = white;
        const int row_y = sy + 124 + vr * 30;
        SDL_Rect track = live_scene_track_rect(sx, sy, sw, vr);
        double value = live_scene_get_value(i);
        double norm = live_scene_value_to_norm(i, value);
        int knob_x = track.x + (int)lround(norm * (double)track.w);
        char value_text[96];

        if (i == 0) snprintf(value_text, sizeof(value_text), "%d", live_scene_wall_count);
        else if (i == 1) snprintf(value_text, sizeof(value_text), "%d", live_scene_particle_count);
        else if (i == 2) snprintf(value_text, sizeof(value_text), "%.4f", (double)live_scene_radius_sigma);
        else if (i == 3) snprintf(value_text, sizeof(value_text), "%.3f", (double)live_scene_eta);
        else if (i == 4) snprintf(value_text, sizeof(value_text), "%.2f", (double)live_scene_spring_l0_sigma);
        else snprintf(value_text, sizeof(value_text), "%.0f", (double)live_scene_wall_mass[i - 5]);

        draw_text_clipped(renderer, font, live_scene_label(i), sx + 10, row_y, row_color, 96);
        draw_text_clipped(renderer, font, value_text, sx + sw - 78, row_y, row_color, 64);

        SDL_SetRenderDrawColor(renderer, 76, 62, 50, 230);
        SDL_RenderFillRect(renderer, &track);
        SDL_Rect fill = track;
        fill.w = knob_x - track.x;
        if (fill.w < 0) fill.w = 0;
        if (fill.w > track.w) fill.w = track.w;
        SDL_SetRenderDrawColor(renderer, 255, 150, 55, 230);
        SDL_RenderFillRect(renderer, &fill);
        SDL_Rect knob = { knob_x - 5, track.y - 5, 10, 18 };
        SDL_SetRenderDrawColor(renderer, 245, 245, 245, 245);
        SDL_RenderFillRect(renderer, &knob);
        if (i == live_scene_drag_index) {
            SDL_SetRenderDrawColor(renderer, 255, 210, 100, 255);
            SDL_RenderDrawRect(renderer, &knob);
        }
    }

    if (live_scene_mode_dropdown_open) {
        for (int opt = 0; opt < LIVE_SCENE_MODE_COUNT; ++opt) {
            SDL_Rect opt_rect = live_scene_mode_option_rect(sx, sy, sw, opt);
            const int selected = (opt == (int)live_scene_mode);
            SDL_SetRenderDrawColor(renderer,
                                   selected ? 210 : 54,
                                   selected ? 95 : 34,
                                   selected ? 34 : 22,
                                   245);
            SDL_RenderFillRect(renderer, &opt_rect);
            SDL_SetRenderDrawColor(renderer, 255, 205, 145, 245);
            SDL_RenderDrawRect(renderer, &opt_rect);
            draw_text_clipped(renderer, font, live_scene_mode_label((LiveSceneMode)opt),
                              opt_rect.x + 7, opt_rect.y + 3, white, opt_rect.w - 12);
        }
    }

    SDL_Rect apply = live_scene_apply_rect(sx, sy, sw);
    SDL_SetRenderDrawColor(renderer, 190, 70, 30, 235);
    SDL_RenderFillRect(renderer, &apply);
    SDL_SetRenderDrawColor(renderer, 255, 220, 180, 255);
    SDL_RenderDrawRect(renderer, &apply);
    draw_text(renderer, font, "APPLY: rebuild scene now", apply.x + 76, apply.y + 5, white);
    {
        char note[160];
        snprintf(note, sizeof(note), "Apply %s/%s: L0=%.2fσ N=%d r=%.4fσ eta=%.3f",
                 live_scene_mode_label(live_scene_mode),
                 live_scene_engine_label(),
                 (double)live_scene_spring_l0_sigma,
                 live_scene_particle_count,
                 (double)live_scene_radius_sigma,
                 (double)live_scene_eta);
        draw_text_clipped(renderer, font, note, sx + 12, apply.y + 34, muted, sw - 24);
    }
}

static int hud_build_experiment_key_lines(char lines[][256], int max_lines) {
    int n = 0;
    if (!lines || max_lines <= 0) return 0;

    switch (cli_experiment_preset) {
        case EXPERIMENT_PRESET_ENERGY_TRANSFER:
            if (n < max_lines) snprintf(lines[n++], 256, "Keys: t=run protocol  e=energy  p=pause  q=quit  +/-/0 speed  F2=sliders");
            if (n < max_lines) snprintf(lines[n++], 256, "Pistons: a/d left  left/right right  s/up stop");
            break;
        case EXPERIMENT_PRESET_SZILARD_ENGINE:
            if (n < max_lines) snprintf(lines[n++], 256, "Szilard: g=separate  h=recombine  o=partial  n=restore");
            if (n < max_lines) snprintf(lines[n++], 256, "Szilard: e=erase  m=measure  y=view  ;/, gate speed  +/-/0 time  F2");
            break;
        case EXPERIMENT_PRESET_PARTICLELIFE:
            if (n < max_lines) snprintf(lines[n++], 256, "ParticleLife: m=tensor mode  p=pause  q=quit  +/-/0 speed  F2=sliders");
            break;
        default:
            if (n < max_lines) snprintf(lines[n++], 256, "Keys: t=protocol step  e=energy  p=pause  q=quit  +/-/0 speed  F2=sliders");
            if (n < max_lines) snprintf(lines[n++], 256, "Pistons: a/d left  left/right right  s/up stop");
            break;
    }

    return n;
}

// ---------- HUD: energy & temperature ----------
static inline void render_energy_hud(SDL_Renderer *r, TTF_Font *f) {
    // compute energies
    const double KEp = kinetic_energy();
    const double KEw = kinetic_energy_walls_total();
    const double KEt = KEp + KEw;

    // temperatures (use runtime active N)
    int n_active = (particles_active > 0) ? particles_active : 1;
    const int n_walls = (wall_enabled ? num_internal_walls : 0);
    const double Tp  = KEp / (n_active * K_B);                                    // particles only (2D)
    const double Tw  = (n_walls > 0) ? (KEw / (0.5 * K_B * (double)n_walls)) : 0; // walls: 1 DOF each
    const double Tt  = KEt / ((n_active + 0.5 * (double)n_walls) * K_B);          // total incl. walls

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
             (double)simulation_time, (double)fixed_dt_runtime, g_last_substeps_used, (double)time_scale_runtime);
    char line_scale[256];
    const double draw_r_px = fmax((double)PARTICLE_RADIUS, (double)gui_particle_draw_min_px);
    snprintf(line_scale, sizeof(line_scale),
             "box: %.1fσ x %.1fσ -> %dpx x %dpx  r: %.4fσ -> %.3fpx draw %.2fpx  scale: %.1f px/σ",
             (double)L0_UNITS * 2.0, (double)HEIGHT_UNITS,
             SIM_WIDTH, SIM_HEIGHT,
             (double)PARTICLE_RADIUS / (double)PIXELS_PER_SIGMA,
             (double)PARTICLE_RADIUS,
             draw_r_px,
             (double)PIXELS_PER_SIGMA);
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
    // Boundary work (EDMD): divider + pistons, cumulative
    double W_div = 0.0, W_L = piston_work_left, W_R = piston_work_right;
    if (g_edmd && (sim_mode == MODE_EDMD || sim_mode == MODE_EDMD_HYBRID)) {
        W_div = edmd_work_divider(g_edmd);
        W_L   = edmd_work_pistonL(g_edmd);
        W_R   = edmd_work_pistonR(g_edmd);
    }
    double W_total = W_L + W_R + W_div;
    char line_work[256];
    snprintf(line_work, sizeof(line_work), "Work: L=%.3e R=%.3e Divider=%.3e", W_L, W_R, W_div);
    char line_work_total[256];
    snprintf(line_work_total, sizeof(line_work_total), "Work_total: %.3e", W_total);

    // Szilard EDMD validity check: show max particle–particle overlap (pixels).
    char line_pp_ov[128] = {0};
    char line_sz_mode[128] = {0};
    char line_sz_branch[256] = {0};
    char line_sz_sep[256] = {0};
    char line_sz_info[256] = {0};
    char line_sz_work[256] = {0};
    if (g_edmd && sim_mode == MODE_EDMD && cli_experiment_preset == EXPERIMENT_PRESET_SZILARD_ENGINE) {
        const double ov = edmd_compute_max_pp_overlap(g_edmd);
        if (ov > 1e-6) snprintf(line_pp_ov, sizeof(line_pp_ov), "PP_overlap_max: %.3g px", ov);
        else           snprintf(line_pp_ov, sizeof(line_pp_ov), "PP_overlap_max: 0");

        int nL0 = 0, nL1 = 0, nR0 = 0, nR1 = 0;
        szilard_species_side_counts(&nL0, &nL1, &nR0, &nR1);
        const double purity_L = szilard_side_purity(nL0, nL1);
        const double purity_R = szilard_side_purity(nR0, nR1);
        const double purity_mean = 0.5 * (purity_L + purity_R);
        const double sep_score = szilard_separation_score(nL0, nL1, nR0, nR1);
        const double I_nats = szilard_mutual_information_nats(nL0, nL1, nR0, nR1);
        const double I_bits = I_nats / log(2.0);
        const double kBT_ref = szilard_reference_kbt();
        int memL0=0,memL1=0,memLU=0,memR0=0,memR1=0,memRU=0;
        szilard_memory_side_counts(&memL0,&memL1,&memLU,&memR0,&memR1,&memRU);
        const double I_mem_nats = szilard_mutual_information_3state_nats(memL0,memL1,memLU,memR0,memR1,memRU);
        const double I_mem_bits = I_mem_nats / log(2.0);
        const int n_species_total = nL0 + nL1 + nR0 + nR1;
        const int n_memory_total = memL0 + memL1 + memLU + memR0 + memR1 + memRU;
        const double F_x = szilard_available_free_energy_from_nats(I_nats, n_species_total);
        const double F_y = szilard_available_free_energy_from_nats(I_mem_nats, n_memory_total);
        int xy_match = 0, xy_mismatch = 0, xy_unknown = 0;
        szilard_xy_agreement_counts(&xy_match, &xy_mismatch, &xy_unknown);
        const double xy_mismatch_frac = szilard_xy_mismatch_fraction();
        const double partial_uncertainty = szilard_partial_uncertainty_fraction_now();
        double mean_speed = 0.0, std_speed = 0.0;
        szilard_speed_stats(&mean_speed, &std_speed);
        double hold_now = 0.0;
        if (sz_sep_perfect_since >= 0.0) {
            hold_now = simulation_time - sz_sep_perfect_since;
            if (hold_now < 0.0) hold_now = 0.0;
        }
        if (hold_now > sz_sep_perfect_hold_time) hold_now = sz_sep_perfect_hold_time;
        const double xeq_hold_target = sz_stage1_baseline_mode ? 0.0 : sz_xeq_ready_hold_time;
        double xeq_hold_now = 0.0;
        if (sz_xeq_ready_since >= 0.0) {
            xeq_hold_now = simulation_time - sz_xeq_ready_since;
            if (xeq_hold_now < 0.0) xeq_hold_now = 0.0;
        }
        if (xeq_hold_now > xeq_hold_target) xeq_hold_now = xeq_hold_target;
        const int mem_ready = szilard_memory_ready();

        snprintf(line_sz_mode, sizeof(line_sz_mode), "Szilard mode: %s  view=%s",
                 szilard_phase_label(sz_interactive_phase),
                 szilard_particle_view_mode_label((int)sz_view_mode));
        const char *mem_state = mem_ready ? "ready" : "unknown";
        if (sz_interactive_phase == 6) mem_state = "erase-ready (press e)";
        else if (sz_interactive_phase == 7) mem_state = "erased";
        snprintf(line_sz_branch, sizeof(line_sz_branch), "Memory: %s  Branch: %s  ReadyBy: %s  Q_erase,min=%.3e",
                 mem_state,
                 g_sz_branch_snapshot.valid ? "READY" : "not ready",
                 szilard_branch_ready_reason_label(sz_branch_ready_reason),
                 sz_memory_erase_heat_min);
        snprintf(line_sz_sep, sizeof(line_sz_sep), "Purity L/R/avg=%.3f/%.3f/%.3f  Sep=%.3f  Ghold=%.1f/%.1f  Xeq=%.1f/%.1f",
                 purity_L, purity_R, purity_mean, sep_score,
                 hold_now, sz_sep_perfect_hold_time,
                 xeq_hold_now, xeq_hold_target);
        snprintf(line_sz_info, sizeof(line_sz_info), "I(x;side)=%.3f bit  I(y;side)=%.3f bit  kBT=%.3f",
                 I_bits, I_mem_bits, kBT_ref);
        if (sz_stage1_baseline_mode) {
            snprintf(line_sz_work, sizeof(line_sz_work), "F_x=%.3e  F_y=%.3e  LoadE=%.3e  LoadW=%.3e  xflip=%s  x!=y=%d/%d (%.3f)",
                     F_x, F_y,
                     sz_load_energy, sz_load_work,
                     sz_xflip_enabled ? "on" : "off",
                     xy_mismatch, xy_match + xy_mismatch, xy_mismatch_frac);
        } else {
            snprintf(line_sz_work, sizeof(line_sz_work), "W_gateMove=%.3e  W_total=%.3e  gate=%s  xflip=%s  x!=y=%d/%d (%.3f)  U=%.3f/%.3f  v_hot>=%.3f",
                     g_edmd ? edmd_work_divider_i(g_edmd, 0) : 0.0, W_total,
                     sz_partial_obs_enabled ? "partial(y~x)" : (sz_perm_use_memory ? "y" : "x"),
                     sz_xflip_enabled ? "on" : "off",
                     xy_mismatch, xy_match + xy_mismatch, xy_mismatch_frac,
                     partial_uncertainty, sz_partial_obs_peak_uncertainty_frac,
                     sz_gate_hot_speed_mult * ((mean_speed > 1e-12) ? mean_speed : fmax(std_speed, 1e-12)));
        }
    }
    // Quick comparison helpers: spring max vs piston input work (efficiency proxy)
    const bool eff_on = energy_measurement.eff_enabled;
    const bool eff_is_spring = eff_on && (energy_measurement.eff_output_mode == EFF_OUTPUT_SPRING);
    const bool eff_is_wall_ke = eff_on && (energy_measurement.eff_output_mode == EFF_OUTPUT_WALL_KE);

    double springE_now  = energy_measurement.spring_energy;
    double springE_max  = energy_measurement.spring_energy_max;
    // In wall-KE mode, "springE_*" fields represent load KE which is already included in KEw.
    double mechE_now    = eff_is_spring ? (KEp + KEw + springE_now) : (KEp + KEw);
    double W_in = piston_push_tracking ? piston_work_delta_max : (piston_work_left + piston_work_right);
    // Windowed gain metric: peak spring energy within a fixed window after piston stops,
    // baseline-subtracted at piston stop (reduces "max so far" drift).
    double dSpring_win = 0.0;
    if (energy_measurement.eff_window_started) {
        dSpring_win = (double)energy_measurement.spring_energy_peak_window - (double)energy_measurement.spring_energy_at_stop;
        if (dSpring_win < 0.0) dSpring_win = 0.0;
    }
    double ratio_gain_win = (W_in > 0.0) ? (dSpring_win / W_in) : 0.0;
    // Keep the all-time peak ratio for reference/debugging.
    double ratio_spring_all = (W_in > 0.0) ? (springE_max / W_in) : 0.0;
    double ratio_mech   = (W_in > 0.0) ? (mechE_now / W_in) : 0.0;
    char line_compare[256];
    snprintf(line_compare, sizeof(line_compare), "Gain_win/W_in=%.3f  Δ%s=%.3e  MechE/W_in=%.3f",
             ratio_gain_win, eff_is_wall_ke ? "LoadKE" : "Spring", dSpring_win, ratio_mech);
    char line_compare2[256];
    snprintf(line_compare2, sizeof(line_compare2), "%s_max(all)/W_in=%.3f",
             eff_is_wall_ke ? "LoadKE" : "Spring", ratio_spring_all);
    char line_compare_defs1[256], line_compare_defs2[256];
    snprintf(line_compare_defs1, sizeof(line_compare_defs1), "Formulas: W_in=piston work  Gain_win=peak(E-E_stop) in Tw=%.3g",
             (double)cli_eff_window_sigma);
    snprintf(line_compare_defs2, sizeof(line_compare_defs2), "          MechE=KEp+KEw%s",
             eff_is_spring ? "+SpringE" : "");

    SDL_Color white = {255,255,255,255};

    // Prefer a right-side "log panel" when the window is wide enough.
    // This avoids cutting off HUD lines when we also draw distributions below the box.
    const int panel_x = gui_scene_screen_x((float)XW2) + 10;
    int hud_right = MAX_X - 10;
    if (live_controls_visible || live_scene_visible) {
        int ux = MAX_X;
        if (live_controls_visible) {
            int cx = 0, cy = 0, cw = 0, ch = 0;
            live_controls_layout(&cx, &cy, &cw, &ch);
            if (cw > 0 && cx < ux) ux = cx;
        }
        if (live_scene_visible) {
            int sx = 0, sy = 0, sw = 0, sh = 0;
            live_scene_layout(&sx, &sy, &sw, &sh);
            if (sw > 0 && sx < ux) ux = sx;
        }
        if (ux < MAX_X) hud_right = ux - 12;
    }
    const int panel_w = hud_right - panel_x;
    const int panel_y = YW1 + 10;
    const bool use_right_panel = (panel_w >= 260);

		    if (use_right_panel) {
		        int x = panel_x;
		        int y = panel_y;
                const int y_max = MAX_Y - 10;
                const int dy = 20;

		        // Compact HUD (single column)
			        draw_text_clipped(r, f, line_clock, x, y, white, panel_w); y += dy;
			        draw_text_clipped(r, f, line_scale, x, y, white, panel_w); y += dy;
			        draw_text_clipped(r, f, line_thermo, x, y, white, panel_w); y += dy;
                // Restore the "classic" thermo/KE block in the right HUD panel (it used to be in the bottom HUD).
                draw_text_clipped(r, f, line0, x, y, white, panel_w); y += dy;

		        char ke_line[256];
		        snprintf(ke_line, sizeof(ke_line), "KEp=%.3e  KEw=%.3e", KEp, KEw);
		        draw_text_clipped(r, f, ke_line, x, y, white, panel_w); y += dy;

        // Per-compartment summary (counts/KE/T) — useful when the text overlay inside compartments is hard to read.
                if (segment_counts && segment_ke && segment_temperature && segment_count > 0) {
                    float left_bounds[MAX_DIVIDER_CAPACITY + 1];
                    float right_bounds[MAX_DIVIDER_CAPACITY + 1];
                    if (segment_count <= MAX_DIVIDER_CAPACITY + 1) {
                        compute_segment_bounds(left_bounds, right_bounds);
                    } else {
                        for (int s = 0; s < MAX_DIVIDER_CAPACITY + 1; ++s) { left_bounds[s] = (float)XW1; right_bounds[s] = (float)XW2; }
                    }
                    const double r_sigma = (double)PARTICLE_RADIUS / (double)PIXELS_PER_SIGMA;
                    const double area_particle = M_PI * r_sigma * r_sigma;
                    for (int seg = 0; seg < segment_count && seg < 6; ++seg) {
                        if (y + dy > y_max) break;
                        const double seg_w_sigma = fmax(1e-9, (double)(right_bounds[seg] - left_bounds[seg]) / (double)PIXELS_PER_SIGMA);
                        const double area_seg = seg_w_sigma * (double)HEIGHT_UNITS;
                        const double eta_seg = (area_seg > 0.0) ? ((double)segment_counts[seg] * area_particle / area_seg) : 0.0;
                        char seg_line[256];
                        snprintf(seg_line, sizeof(seg_line), "Box %d: N=%d  KE=%.3e  T=%.2f  eta=%.3f",
                                 seg + 1,
                                 segment_counts[seg],
                                 (double)segment_ke[seg],
                                 (double)segment_temperature[seg],
                                 eta_seg);
                        draw_text_clipped(r, f, seg_line, x, y, white, panel_w);
                        y += dy;

                        // Optional: compact velocity distribution parameters (mean/std) for quick sanity checks.
                        if (segment_vx_mean && segment_vx_var && segment_speed_mean && segment_speed_var) {
                            if (y + dy > y_max) break;
                            const double sx = sqrt(fmax(0.0, segment_vx_var[seg]));
                            const double ss = sqrt(fmax(0.0, segment_speed_var[seg]));
                            // Approximate median(|v|) assuming Rayleigh (good MB sanity check in 2D).
                            // Rayleigh: Var = (4-π)/2 * σ², Median = σ*sqrt(2 ln 2).
                            double v_med = segment_speed_mean[seg];
                            if (segment_speed_var[seg] > 0.0) {
                                const double denom = (4.0 - M_PI);
                                if (denom > 1e-12) {
                                    double sigma_hat = sqrt(fmax(0.0, segment_speed_var[seg]) * 2.0 / denom);
                                    v_med = sigma_hat * sqrt(2.0 * log(2.0));
                                }
                            }
                            char vel_line[256];
                            snprintf(vel_line, sizeof(vel_line), "   vx(mu=%.2f sd=%.2f)  |v|(mu=%.2f med=%.2f sd=%.2f)",
                                     (double)segment_vx_mean[seg], sx,
                                     (double)segment_speed_mean[seg], v_med, ss);
                            draw_text_clipped(r, f, vel_line, x, y, white, panel_w);
                            y += dy;
                        }
                    }
                }

		        if (y + dy <= y_max) { draw_text_clipped(r, f, line6, x, y, white, panel_w); y += dy; } // piston work delta
                if (y + dy <= y_max) { draw_text_clipped(r, f, line_work_total, x, y, white, panel_w); y += dy; }
                if (y + dy <= y_max) { draw_text_clipped(r, f, line_compare, x, y, white, panel_w); y += dy; }
                if (y + dy <= y_max) { draw_text_clipped(r, f, line_compare2, x, y, white, panel_w); y += dy; }
                if (line_sz_mode[0] && y + dy <= y_max) { draw_text_clipped(r, f, line_sz_mode, x, y, white, panel_w); y += dy; }
                if (line_sz_branch[0] && y + dy <= y_max) { draw_text_clipped(r, f, line_sz_branch, x, y, white, panel_w); y += dy; }
                if (line_sz_sep[0] && y + dy <= y_max) { draw_text_clipped(r, f, line_sz_sep, x, y, white, panel_w); y += dy; }
                if (line_sz_info[0] && y + dy <= y_max) { draw_text_clipped(r, f, line_sz_info, x, y, white, panel_w); y += dy; }
                if (line_sz_work[0] && y + dy <= y_max) { draw_text_clipped(r, f, line_sz_work, x, y, white, panel_w); y += dy; }
                if (line_pp_ov[0] && y + dy <= y_max) { draw_text_clipped(r, f, line_pp_ov, x, y, white, panel_w); y += dy; }

	        // Spring energy summary
	        char spring_line[256];
	        if (eff_is_wall_ke) {
            snprintf(spring_line, sizeof(spring_line), "LoadKE: %.3e (max %.3e)",
                     (double)energy_measurement.spring_energy,
                     (double)energy_measurement.spring_energy_max);
        } else if (energy_measurement.spring_enabled) {
            snprintf(spring_line, sizeof(spring_line), "SpringE: %.3e (max %.3e)",
                     (double)energy_measurement.spring_energy,
                     (double)energy_measurement.spring_energy_max);
	        } else {
	            snprintf(spring_line, sizeof(spring_line), "SpringE: (disabled)");
	        }
		        if (y + dy <= y_max) { draw_text_clipped(r, f, spring_line, x, y, white, panel_w); y += dy; }

	        // Thermostat status
	        char hb_line[256], andersen_line[256];
	        if (heatbath_enabled) {
            snprintf(hb_line, sizeof(hb_line), "HB: ON (T=%.3f, mode=%d)",
                     (double)heatbath_temperature, thermal_wall_mode);
	        } else {
	            snprintf(hb_line, sizeof(hb_line), "HB: off");
	        }
		        if (y + dy <= y_max) { draw_text_clipped(r, f, hb_line, x, y, white, panel_w); y += dy; }
	        if (andersen_enabled) {
	            snprintf(andersen_line, sizeof(andersen_line), "Andersen: ON (ν=%.3f)", andersen_collision_freq);
	        } else {
	            snprintf(andersen_line, sizeof(andersen_line), "Andersen: off");
	        }
		        if (y + dy <= y_max) { draw_text_clipped(r, f, andersen_line, x, y, white, panel_w); y += dy; }

	        // Windowed metric info
        if (energy_measurement.eff_window_started) {
            char win_line[256];
            double t_stop_rel = (wall_release_time >= 0.0) ? (energy_measurement.eff_piston_stop_time - (float)wall_release_time) : NAN;
            snprintf(win_line, sizeof(win_line), "t_stop(rel)=%.3f  Tw=%.3g",
                     (double)t_stop_rel, (double)cli_eff_window_sigma);
		            if (y + dy <= y_max) { draw_text_clipped(r, f, win_line, x, y, white, panel_w); y += dy; }

	            // Peak timing diagnostics (helps choose Tw large enough to capture the true peak)
	            char peak_line[256];
	            double t_peak_rel = (wall_release_time >= 0.0) ? (energy_measurement.spring_energy_peak_window_time - (float)wall_release_time) : NAN;
	            double dt_peak = (isfinite(t_stop_rel) && isfinite(t_peak_rel)) ? (t_peak_rel - t_stop_rel) : NAN;
            snprintf(peak_line, sizeof(peak_line), "t_peak(rel)=%.3f  Δt_peak=%.3f",
                     (double)t_peak_rel, (double)dt_peak);
	            if (y + dy <= y_max) { draw_text_clipped(r, f, peak_line, x, y, white, panel_w); y += dy; }
	        }
            {
                char key_lines[3][256];
                int key_n = hud_build_experiment_key_lines(key_lines, 3);
                for (int i = 0; i < key_n; ++i) {
                    if (y + dy > y_max) break;
                    draw_text_clipped(r, f, key_lines[i], x, y, white, panel_w);
                    y += dy;
                }
            }
	        return;
	    }

    // Fallback: draw HUD under the simulation window (older layout)
    int x = XW1 + 5;
    int y = YW1 + SIM_HEIGHT + 5;
    draw_text(r, f, line0, x, y, white);  y += 20;
    draw_text(r, f, line_clock, x, y, white); y += 20;
    draw_text(r, f, line_scale, x, y, white); y += 20;
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
    {
        char key_lines[3][256];
        int key_n = hud_build_experiment_key_lines(key_lines, 3);
        for (int i = 0; i < key_n; ++i) {
            draw_text(r, f, key_lines[i], x, y, white);
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

    float packing_unitless = particles_active * particle_area_unitless / box_area_unitless;
    float packing_actual_bound = particles_active * particle_area_actual / box_area_scaled;
    float packing_actual_excl_wall = particles_active * particle_area_actual / (box_area_scaled - wall_buffer_area);

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
    // Create timestamped folder to preserve old data
    time_t now = time(NULL);
    struct tm *t = localtime(&now);
    char timestamp[64];
    snprintf(timestamp, sizeof(timestamp), "simulation_%02d_%02d_%02d_%02d_%02d_%02d",
             t->tm_mday, t->tm_mon + 1, t->tm_year % 100, t->tm_hour, t->tm_min, t->tm_sec);

    const bool is_edmd_like = (sim_mode == MODE_EDMD || sim_mode == MODE_EDMD_HYBRID);
    const char* engine_folder = is_edmd_like ? "EDMD" : "TIME";

    char mode0_folder[512], mode1_folder[512];
    snprintf(mode0_folder, sizeof(mode0_folder),
             "experiments_speed_of_sound/%s/mode0_real_units/%s", engine_folder, timestamp);
    if (cli_speed_sound_run_dir && cli_speed_sound_run_dir[0]) {
        resolve_cli_output_path(mode1_folder, sizeof(mode1_folder), cli_speed_sound_run_dir);
    } else {
        snprintf(mode1_folder, sizeof(mode1_folder),
                 "experiments_speed_of_sound/%s/mode1_normalized_units/%s", engine_folder, timestamp);
    }

    // Create directories (but DON'T delete old data!)
    mkdir_p(mode0_folder);
    mkdir_p(mode1_folder);

    // Always record how this simulation batch was generated.
    write_command_md_to_dir(mode0_folder, g_main_argc, g_main_argv);
    write_command_md_to_dir(mode1_folder, g_main_argc, g_main_argv);
    write_speed_of_sound_plot_commands_to_dir(mode0_folder);
    write_speed_of_sound_plot_commands_to_dir(mode1_folder);

    #if MOLECULAR_MODE
        printf("📁 Output folder: %s\n", mode1_folder);
    #else
        printf("📁 Output folder: %s\n", mode0_folder);
    #endif

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

                // Use timestamped folder
                #if MOLECULAR_MODE
                    const char* mode_folder = mode1_folder;
                #else
                    const char* mode_folder = mode0_folder;
                #endif

                char filename[512];
                int L0_int = (int)(L0 * 10);
                snprintf(filename, sizeof(filename), "%s/wall_x_positions_L0_%d_wallmassfactor_%d_run%d.csv",
                         mode_folder, L0_int, wall_mass_factor, r);
                FILE *wall_log = fopen(filename, "w");
                if (!wall_log) { printf("❌ Could not open log file.\n"); continue; }

                // Keep this CSV schema stable (analysis scripts assume these columns).
                // Do NOT write substep logs into the same file (it breaks CSV parsing).
                g_substep_log = NULL;
                g_substep_base_time = 0.0;
                // First 5 columns are consumed by analysis scripts; append metadata after that.
                // Wall_X is an absolute σ-coordinate in the simulation frame (it is NOT equal to L0 if XW1 != 0).
                // Displacement(σ) = Wall_X - Center_X(σ), so ~0 means "divider centered".
                fprintf(wall_log, "Time,Wall_X,Displacement(σ),Left_Count,Right_Count,L0,eta,Center_X(σ)\n");

                printf("🔬 Running: L0 = %.1f, M = %d*m, run = %d\n", L0, wall_mass_factor, r);
                printf("🔍 Initial wall_x = %.3f, vx_wall = %.6f\n", wall_x, vx_wall);
                fflush(stdout);

                int left_particles = 0, right_particles = 0;
                int recorded_steps = 0;
                missed_collision_events = 0;
                worst_penetration_observed = 0.0;
                const double center_x_sigma_const = ((double)XW1 + (double)XW2) / (2.0 * (double)PIXELS_PER_SIGMA);
                const double height_sigma_const = ((double)YW2 - (double)YW1) / (double)PIXELS_PER_SIGMA;
                const double eta_nominal_const = ((double)particles_active * M_PI * (double)PARTICLE_RADIUS * (double)PARTICLE_RADIUS)
                    / fmax(1e-12, (2.0 * (double)L0_UNITS * height_sigma_const));

                if (sim_mode == MODE_EDMD || sim_mode == MODE_EDMD_HYBRID) {
                    if (sim_mode == MODE_EDMD_HYBRID) {
                        printf("⚠️ speed_of_sound requires boundary events; switching EDMD-HYBRID -> EDMD.\n");
                        sim_mode = MODE_EDMD;
                    }

                    // Build EDMD state from the currently seeded particles/wall.
                    if (g_edmd) { edmd_backend_destroy(g_edmd); g_edmd = NULL; }
                    EDMD_Params prm = {0};
                    prm.boxW = (double)(XW2 - XW1);
                    prm.boxH = (double)(YW2 - YW1);
                    prm.radius = (double)PARTICLE_RADIUS;
                    prm.N = particles_active;
                    prm.cell_size = 0.0; // auto
                    prm.divider_count = 1;
                    prm.divider_x[0] = (double)(wall_x - (double)XW1);
                    prm.divider_thickness[0] = (double)WALL_THICKNESS;
                    prm.divider_mass[0] = 0.0; // held phase: infinite mass (fixed divider)
                    prm.divider_vx[0] = 0.0;
                    prm.heatbath_enabled = heatbath_enabled;
                    prm.heatbath_temperature = (double)heatbath_temperature;
                    prm.thermal_wall_mode = thermal_wall_mode;
                    prm.mb_overshoot_factor = (double)mb_overshoot_factor;
                    prm.stability_window_percent = (double)stability_window_percent;
                    prm.particle_mass = (double)PARTICLE_MASS;
                    prm.kB = (double)kB_effective();
                    prm.pp_collisions_enabled = cli_no_pp_collisions ? 0 : 1;
                    g_edmd = edmd_backend_create(&prm);
                    if (!g_edmd) { fprintf(stderr, "EDMD create failed.\n"); exit(1); }
                    // Load state
                    EDMD_Particle* P = (EDMD_Particle*)edmd_backend_particles(g_edmd);
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
                    edmd_backend_divider_resolve_overlaps(g_edmd);
                    edmd_backend_reschedule_all(g_edmd);

                    // Hold/equilibrate with the divider fixed for wall_hold_steps "ticks"
                    wall_is_released = false;
                    wall_hold_enabled = true;
                    for (int s = 0; s < wall_hold_steps; ++s) {
                        double t0 = edmd_backend_time(g_edmd);
                        edmd_backend_advance_to(g_edmd, t0 + (double)fixed_dt_runtime);
                    }

                    // Release: set finite wall mass
                    {
                        const double mass = (double)wall_mass_factor;
                        const double vx = 0.0;
                        edmd_backend_set_divider_motions(g_edmd, 1, &mass, &vx);
                        edmd_backend_reschedule_all(g_edmd);
                    }
                    wall_release_time = edmd_backend_time(g_edmd) / (double)PIXELS_PER_SIGMA;
                    wall_is_released = true;
                    wall_hold_enabled = false;
                    recorded_steps = 0;
                    if (!cli_quiet) {
                        printf("🔔 Wall released (EDMD) at t = %.3f\n", wall_release_time);
                    }

                    while (recorded_steps < num_steps) {
                        // Advance EDMD to the next sample time (sampling interval = fixed_dt_runtime).
                        double t0 = edmd_backend_time(g_edmd);
                        edmd_backend_advance_to(g_edmd, t0 + (double)fixed_dt_runtime);
                        edmd_backend_divider_resolve_overlaps(g_edmd);

                        const EDMD_Params* ep = edmd_backend_params(g_edmd);
                        const EDMD_Particle* Pnow = edmd_backend_particles(g_edmd);
                        const double div_x_px = ep->divider_x[0];
                        int leftN = 0, rightN = 0;
                        for (int i = 0; i < particles_active; ++i) {
                            if (Pnow[i].x < div_x_px) leftN++;
                            else rightN++;
                        }
                        left_particles = leftN;
                        right_particles = rightN;

                        // Convert to σ-units for logging (Time is σ-time after release).
                        const double t_sigma = edmd_backend_time(g_edmd) / (double)PIXELS_PER_SIGMA;
                        const double wall_x_sigma = ((double)XW1 + div_x_px) / (double)PIXELS_PER_SIGMA;
                        const double disp = wall_x_sigma - ((double)XW1 + (double)XW2) / (2.0 * (double)PIXELS_PER_SIGMA);
                        const double t_after = t_sigma - wall_release_time;
                        fprintf(wall_log, "%.6f,%.6f,%.6f,%d,%d,%.6f,%.6f,%.6f\n",
                                t_after, wall_x_sigma, disp, left_particles, right_particles,
                                (double)L0_UNITS, eta_nominal_const, center_x_sigma_const);

                        recorded_steps++;

                        // bail-out if wall didn't get any motion for a long time
                        if (fabs(ep->divider_vx[0]) < 1e-15 && recorded_steps > wall_hold_steps * 10) {
                            printf("⚠️ Wall not moving — exiting early.\n");
                            break;
                        }
                    }
                } else {
                    while (recorded_steps < num_steps) {
                        // --- START OF STEP ---
                        wall_x_old = wall_x;   // ← CCD uses wall pose from start of step
                        steps_elapsed++;       // ← tick BEFORE update_wall so hold/release can trigger

                        // IMPORTANT: update piston *velocity* (protocol) before particle CCD, so CCD uses the correct moving piston.
                        apply_piston_protocol(fixed_dt_runtime);

                        // physics
                        wall_impulse_x_accum = 0.f; // (re)start accumulator for this step
                        update_particles_with_substepping(fixed_dt_runtime, &left_particles, &right_particles,
                                                          simulation_time, wall_is_released, wall_release_time);

                        // Advance piston positions after collision handling for this step.
                        move_left_piston(fixed_dt_runtime);
                        move_right_piston(fixed_dt_runtime);
                        update_wall(fixed_dt_runtime, L0_UNITS);
                        update_extra_walls(fixed_dt_runtime, L0_UNITS);
                        sync_all_wall_positions();
                        update_wall_metrics();

                        // ##CHRIS: time bookkeeping - convert pixel-time to σ-time
                        simulation_time += fixed_dt_runtime / PIXELS_PER_SIGMA;

                        // release check using true sim time
                        if (steps_elapsed == wall_hold_steps && !wall_is_released) {
                            wall_release_time = simulation_time;
                            wall_is_released  = true;
                            recorded_steps    = 0;
                            if (!cli_quiet) {
                                printf("🔔 Wall released at step %d (t = %.3f)\n", steps_elapsed, wall_release_time);
                            }
                        }

                        if (!wall_is_released) continue;

                        float time_after_release = (float)(simulation_time - wall_release_time);
                        float wall_x_sigma = wall_x / PIXELS_PER_SIGMA;
                        float disp = wall_x_sigma - (XW1 + XW2) / (2.0f * PIXELS_PER_SIGMA);

                        fprintf(wall_log, "%.6f,%.6f,%.6f,%d,%d,%.6f,%.6f,%.6f\n",
                                time_after_release, wall_x_sigma, disp, left_particles, right_particles,
                                (double)L0_UNITS, eta_nominal_const, center_x_sigma_const);

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
                }

                fclose(wall_log);
                printf("✅ Finished run %d\n\n", r);
                fflush(stdout);
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
    mkdir_p(g_energy_transfer_dir);
    write_command_md_to_dir(g_energy_transfer_dir, g_main_argc, g_main_argv);

    if (!energy_measurement.eff_enabled) {
        initialize_energy_measurement();
    }

    wall_enabled = 1;
    wall_position_is_managed_externally = 0;

    int left_particles = 0, right_particles = 0;

    char trace_path[512];
    if (cli_energy_transfer_trace_set && cli_energy_transfer_trace_path[0] != '\0') {
        resolve_cli_output_path(trace_path, sizeof(trace_path), cli_energy_transfer_trace_path);
    } else {
        snprintf(trace_path, sizeof(trace_path), "%s/energy_transfer_trace.csv", g_energy_transfer_dir);
    }
    mkdir_p_for_file(trace_path);
    FILE *elog = fopen(trace_path, "w");
    if (!elog) { printf("❌ Could not open energy transfer log file.\n"); return; }
    // Keep the first 10 columns stable for older scripts; append debug columns for wall/piston kinematics.
    fprintf(elog,
            "Time,Wall_X,SpringF,SpringE,EnergyTransferred,Left_Count,Right_Count,PistonWork,PistonWorkDeltaMax,SpringE_Max,"
            "W0_x_sigma,W0_v,W1_x_sigma,W1_v,W2_x_sigma,W2_v,W3_x_sigma,W3_v,W4_x_sigma,W4_v,"
            "PistonR_x_sigma,PistonR_v,PistonL_x_sigma,PistonL_v,"
            "PistonStop_t_rel,SpringE_stop,SpringE_peak_win,SpringE_delta_peak_win,SegCounts,SegEtas\n");

    // Reset sim
    initialize_simulation_dimensions();
    initialize_simulation_parameters();
    initialize_simulation();
    steps_elapsed = 0;
    wall_is_released = false;
    simulation_time = 0.0f;
    wall_release_time = -1.0;
    wall_x_old = wall_x;
    missed_collision_events = 0; worst_penetration_observed = 0.0;
    piston_step_active = false;
    vx_piston_right = 0.0f;
    vx_piston_left = 0.0f;
    piston_work_left = 0.0;
    piston_work_right = 0.0;
    piston_push_tracking = false;
    piston_work_baseline = 0.0;
    piston_work_delta_max = 0.0;
    energy_measurement.spring_energy_max = 0.0f;
    energy_measurement.energy_transferred = 0.0f;
    energy_measurement.eff_window_pending_start = false;
    energy_measurement.eff_window_started = false;
    energy_measurement.eff_window_active = false;
    energy_measurement.eff_window_done = false;
    energy_measurement.eff_piston_stop_time = 0.0f;
    energy_measurement.spring_energy_at_stop = 0.0f;
    energy_measurement.spring_energy_peak_window = 0.0f;

    int recorded_steps = 0;
    int target_steps = num_steps > 0 ? num_steps : (int)(3000.0f / fmaxf(fixed_dt_runtime, 1e-9f));
    bool piston_auto_triggered = false;

    // EDMD work tracking (if enabled)
    double edmd_prev_work_L = 0.0, edmd_prev_work_R = 0.0, edmd_prev_work_div = 0.0;

    // Initialize EDMD state for energy_transfer experiment if requested
    if (sim_mode == MODE_EDMD || sim_mode == MODE_EDMD_HYBRID) {
        if (num_internal_walls > EDMD_MAX_DIVIDERS && sim_mode == MODE_EDMD) {
            printf("⚠️ %d walls exceed EDMD capacity (%d); using EDMD-HYBRID.\n",
                   num_internal_walls, EDMD_MAX_DIVIDERS);
            sim_mode = MODE_EDMD_HYBRID;
        }
        if (g_edmd) { edmd_destroy(g_edmd); g_edmd = NULL; }
        EDMD_Params prm = {0};
        prm.boxW = (double)(XW2 - XW1);
        prm.boxH = (double)(YW2 - YW1);
        prm.radius = (double)PARTICLE_RADIUS;
        prm.N = particles_active;
        prm.cell_size = 0.0; // auto
        {
            int dcount = (num_internal_walls > 0 ? num_internal_walls : 0);
            if (dcount > EDMD_MAX_DIVIDERS) dcount = EDMD_MAX_DIVIDERS;
            prm.divider_count = dcount;
            int primary_idx = primary_wall_index;
            if (primary_idx < 0 || primary_idx >= dcount) primary_idx = 0;
            int extra_idx = 0;
            for (int w = 0; w < dcount; ++w) {
                double pos_px = (all_wall_positions ? (double)all_wall_positions[w] : (double)wall_x);
                prm.divider_x[w] = pos_px - (double)XW1;
                prm.divider_thickness[w] = (double)WALL_THICKNESS;
                double mf = (double)(WALL_MASS / PARTICLE_MASS);
                if (cli_wall_masses_set && w < (int)(sizeof(cli_wall_masses) / sizeof(cli_wall_masses[0])) && cli_wall_masses[w] > 0.0f) {
                    mf = (double)cli_wall_masses[w];
                }
                prm.divider_mass[w] = mf;
                double vxw = 0.0;
                if (w == primary_idx) {
                    vxw = (double)vx_wall;
                } else if (extra_wall_velocity && extra_idx < extra_wall_count) {
                    vxw = (double)extra_wall_velocity[extra_idx];
                    extra_idx++;
                }
                prm.divider_vx[w] = vxw;
            }
        }
        prm.heatbath_enabled = heatbath_enabled;
        prm.heatbath_temperature = (double)heatbath_temperature;
        prm.thermal_wall_mode = thermal_wall_mode;
        prm.mb_overshoot_factor = (double)mb_overshoot_factor;
        prm.stability_window_percent = (double)stability_window_percent;
        prm.particle_mass = (double)PARTICLE_MASS;
        prm.kB = (double)kB_effective();
        // EDMD harmonic spring wall: if we're outputting spring energy, model the primary divider as a
        // true harmonic oscillator inside EDMD (no timestep impulses).
        if (energy_measurement.spring_enabled && prm.divider_count > 0) {
            int primary_idx = primary_wall_index;
            if (primary_idx < 0 || primary_idx >= prm.divider_count) primary_idx = 0;
            const double xeq_px = (energy_measurement.equilibrium_position != 0.0f)
                                    ? (double)energy_measurement.equilibrium_position
                                    : (double)wall_x;
            prm.divider_k[primary_idx] = (double)energy_measurement.spring_constant;
            prm.divider_xeq[primary_idx] = xeq_px - (double)XW1;
        }
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
        // Szilard: keep EDMD in a valid (overlap-free) state even if the seeding/clamp left tiny overlaps.
        edmd_relax_pp_overlaps_szilard(g_edmd);
        if (sim_mode == MODE_EDMD) {
            edmd_reschedule_all(g_edmd);
            edmd_reset_work(g_edmd);
            printf("[EDMD] initialized: N=%d  box=(%.1fσ,%.1fσ) -> (%.1fpx,%.1fpx)  R=%.4fσ -> %.3fpx  scale=%.1f px/σ\n",
                   prm.N,
                   prm.boxW / (double)PIXELS_PER_SIGMA,
                   prm.boxH / (double)PIXELS_PER_SIGMA,
                   prm.boxW, prm.boxH,
                   prm.radius / (double)PIXELS_PER_SIGMA,
                   prm.radius,
                   (double)PIXELS_PER_SIGMA);
        } else {
            edmd_reschedule_all_pp_only(g_edmd);
            edmd_reset_work(g_edmd);
            printf("[EDMD-HYBRID] initialized (PP-only): N=%d  box=(%.1fσ,%.1fσ) -> (%.1fpx,%.1fpx)  R=%.4fσ -> %.3fpx  scale=%.1f px/σ\n",
                   prm.N,
                   prm.boxW / (double)PIXELS_PER_SIGMA,
                   prm.boxH / (double)PIXELS_PER_SIGMA,
                   prm.boxW, prm.boxH,
                   prm.radius / (double)PIXELS_PER_SIGMA,
                   prm.radius,
                   (double)PIXELS_PER_SIGMA);
        }
    }

    // Tunneling detection: snapshot initial segment counts at wall release
    int *initial_segment_counts = NULL;
    bool tunneling_detected = false;

    while (recorded_steps < target_steps) {
        wall_x_old = wall_x;
        steps_elapsed++;

        wall_impulse_x_accum = 0.f;
        if (sim_mode == MODE_EDMD || sim_mode == MODE_EDMD_HYBRID) {
            // EDMD / EDMD-HYBRID path
            if (sim_mode == MODE_EDMD) {
                const EDMD_Params* ep_state = edmd_params(g_edmd);
                int dcount = ep_state->divider_count;
                if (!wall_enabled) dcount = 0;
                if (dcount > EDMD_MAX_DIVIDERS) dcount = EDMD_MAX_DIVIDERS;
                int primary_idx = primary_wall_index;
                if (primary_idx < 0 || primary_idx >= dcount) primary_idx = 0;

                double div_x[EDMD_MAX_DIVIDERS];
                double div_vx[EDMD_MAX_DIVIDERS];
                double div_mass[EDMD_MAX_DIVIDERS];
                double div_th[EDMD_MAX_DIVIDERS];
                int extra_idx = 0;

                for (int w = 0; w < dcount; ++w) {
                    div_x[w] = ep_state->divider_x[w];
                    div_vx[w] = ep_state->divider_vx[w];
                    div_th[w] = (double)WALL_THICKNESS;
                    double mf = (double)(WALL_MASS / PARTICLE_MASS);
                    if (cli_wall_masses_set && w < (int)(sizeof(cli_wall_masses) / sizeof(cli_wall_masses[0])) && cli_wall_masses[w] > 0.0f) {
                        mf = (double)cli_wall_masses[w];
                    }
                    // Hold walls (infinite mass) before release
                    bool hold_active = false;
                    if (w == primary_idx) {
                        hold_active = wall_hold_enabled && !wall_is_released;
                    } else if (extra_wall_hold_enabled && extra_wall_is_released) {
                        int j = (w < primary_idx) ? w : (w - 1);
                        if (j >= 0 && j < extra_wall_count) {
                            hold_active = extra_wall_hold_enabled[j] && !extra_wall_is_released[j];
                        }
                    }
                    if (hold_active) {
                        div_mass[w] = 0.0;
                        div_vx[w] = 0.0;
                    } else {
                        div_mass[w] = mf;
                        if (extra_wall_velocity && extra_idx < extra_wall_count && w != primary_idx) {
                            div_vx[w] = (double)extra_wall_velocity[extra_idx];
                            extra_idx++;
                        }
                    }
                }

                piston_left_x  = (float)(XW1 + ep_state->pistonL_x - 5.0);
                piston_right_x = (float)(XW1 + ep_state->pistonR_x);

                // In pure EDMD mode, the spring wall is handled by the EDMD core (harmonic divider)
                // when --eff-output=spring is used. Do not apply timestep impulses here.

                edmd_config_dividers(g_edmd, dcount, div_x, div_th);
                edmd_set_divider_motions(g_edmd, dcount, div_mass, div_vx);

                static double prev_vxL = 0.0, prev_vxR = 0.0;
                apply_piston_protocol(fixed_dt_runtime);
                if (vx_piston_left != prev_vxL || vx_piston_right != prev_vxR) {
                    const EDMD_Params* epc = edmd_params(g_edmd);
                    edmd_config_pistons(g_edmd,
                        1, epc->pistonL_x, (double)vx_piston_left, epc->pistonL_mass,
                        1, epc->pistonR_x, (double)vx_piston_right, epc->pistonR_mass);
                    prev_vxL = vx_piston_left; prev_vxR = vx_piston_right;
                }
                edmd_reschedule_all(g_edmd);

                double t0 = edmd_time(g_edmd);
                edmd_advance_to(g_edmd, t0 + (double)fixed_dt_runtime);
                edmd_divider_resolve_overlaps(g_edmd);
                const EDMD_Particle* P = edmd_particles(g_edmd);
                for (int i=0;i<particles_active;i++) {
                    X[i] = (double)XW1 + (double)P[i].x;
                    Y[i] = (double)YW1 + (double)P[i].y;
                    Vx[i] = (double)P[i].vx;
                    Vy[i] = (double)P[i].vy;
                }
                const EDMD_Params* ep = edmd_params(g_edmd);
                int dcount2 = ep->divider_count;
                if (!wall_enabled) dcount2 = 0;
                if (dcount2 > EDMD_MAX_DIVIDERS) dcount2 = EDMD_MAX_DIVIDERS;
                int primary_idx2 = primary_wall_index;
                if (primary_idx2 < 0 || primary_idx2 >= dcount2) primary_idx2 = 0;
                int extra_idx2 = 0;
                for (int w = 0; w < dcount2; ++w) {
                    if (w == primary_idx2) {
                        wall_x = (double)XW1 + ep->divider_x[w];
                        vx_wall = (float)ep->divider_vx[w];
                    } else if (extra_wall_positions && extra_wall_velocity && extra_idx2 < extra_wall_count) {
                        extra_wall_positions[extra_idx2] = (float)((double)XW1 + ep->divider_x[w]);
                        extra_wall_velocity[extra_idx2] = (float)ep->divider_vx[w];
                        extra_idx2++;
                    }
                }
                piston_left_x  = (float)(XW1 + ep->pistonL_x - 5.0);
                piston_right_x = (float)(XW1 + ep->pistonR_x);
                double wL = edmd_work_pistonL(g_edmd);
                double wR = edmd_work_pistonR(g_edmd);
                double wD = edmd_work_divider(g_edmd);
                piston_work_left  += (wL - edmd_prev_work_L);
                piston_work_right += (wR - edmd_prev_work_R);
                edmd_prev_work_L = wL;
                edmd_prev_work_R = wR;
                edmd_prev_work_div = wD;
                recompute_segment_stats_counts_and_temperature();
                if (segment_counts && segment_count >= 2) {
                    left_particles = segment_counts[0];
                    right_particles = segment_counts[1];
                }
            } else {
                // EDMD-HYBRID
                double t0 = edmd_time(g_edmd);
                edmd_advance_pp_only_to(g_edmd, t0 + (double)fixed_dt_runtime);
                const EDMD_Particle* P = edmd_particles(g_edmd);
                for (int i=0;i<particles_active;i++) {
                    X[i] = (double)XW1 + (double)P[i].x;
                    Y[i] = (double)YW1 + (double)P[i].y;
                    Vx[i] = (double)P[i].vx;
                    Vy[i] = (double)P[i].vy;
                }
                hybrid_faces_correction_pass(fixed_dt_runtime);
                EDMD_Particle* Pw = (EDMD_Particle*)edmd_particles(g_edmd);
                for (int i=0;i<particles_active;i++) {
                    Pw[i].x = (double)(X[i] - XW1);
                    Pw[i].y = (double)(Y[i] - YW1);
                    Pw[i].vx = (double)Vx[i];
                    Pw[i].vy = (double)Vy[i];
                }
                edmd_reschedule_all_pp_only(g_edmd);
                recompute_segment_stats_counts_and_temperature();
                if (segment_counts && segment_count >= 2) {
                    left_particles = segment_counts[0];
                    right_particles = segment_counts[1];
                }
            }
            // Sync walls + metrics after EDMD step
            sync_all_wall_positions();
            update_wall_metrics();
            update_energy_measurement(fixed_dt_runtime);
        } else {
            // TIME / RK4 path
            // IMPORTANT: update piston *velocity* (protocol) before particle CCD, so CCD uses the correct moving piston.
            apply_piston_protocol(fixed_dt_runtime);
            g_substep_base_time = simulation_time; // ##CHRIS: Update base time for substep logging
            update_particles_with_substepping(fixed_dt_runtime, &left_particles, &right_particles,
                                              simulation_time, wall_is_released, wall_release_time);

            // Advance piston positions after collision handling for this step.
            move_left_piston(fixed_dt_runtime);
            move_right_piston(fixed_dt_runtime);
            update_wall(fixed_dt_runtime, L0_UNITS);
            update_extra_walls(fixed_dt_runtime, L0_UNITS);
            sync_all_wall_positions();
            update_wall_metrics();
            update_energy_measurement(fixed_dt_runtime);
        }

        // ##CHRIS: Convert pixel-time to σ-time units (TIME mode fix)
        simulation_time += fixed_dt_runtime / PIXELS_PER_SIGMA;

        if (steps_elapsed == wall_hold_steps && !wall_is_released) {
            wall_release_time = simulation_time;
            wall_is_released  = true;
            recorded_steps    = 0;
            if (!cli_quiet) {
                printf("🔔 Wall released at step %d (t = %.3f)\n", steps_elapsed, wall_release_time);
            }
            if (extra_wall_hold_enabled && extra_wall_is_released) {
                for (int j = 0; j < extra_wall_count; ++j) {
                    extra_wall_is_released[j] = true;
                    extra_wall_hold_enabled[j] = false;
                    if (extra_wall_release_time) extra_wall_release_time[j] = simulation_time;
                }
            }
            // Snapshot per-segment particle counts at release for tunneling detection
            if (segment_counts && segment_count > 0 && !initial_segment_counts) {
                initial_segment_counts = (int*)malloc((size_t)segment_count * sizeof(int));
                if (initial_segment_counts) {
                    memcpy(initial_segment_counts, segment_counts, (size_t)segment_count * sizeof(int));
                    if (!cli_quiet) {
                        printf("  Tunneling check: initial segment counts [");
                        for (int s = 0; s < segment_count; s++)
                            printf("%s%d", s ? ", " : "", initial_segment_counts[s]);
                        printf("]\n");
                    }
                }
            }
        }

        if (!wall_is_released) continue;

        // Tunneling detection: compare current segment counts to initial
        if (initial_segment_counts && segment_counts && segment_count > 0) {
            for (int s = 0; s < segment_count; s++) {
                if (segment_counts[s] != initial_segment_counts[s]) {
                    fprintf(stderr, "TUNNELING DETECTED at step %d (t=%.6f): "
                            "segment %d count changed from %d to %d\n",
                            steps_elapsed, simulation_time,
                            s, initial_segment_counts[s], segment_counts[s]);
                    fprintf(stderr, "  Initial: [");
                    for (int j = 0; j < segment_count; j++)
                        fprintf(stderr, "%s%d", j ? ", " : "", initial_segment_counts[j]);
                    fprintf(stderr, "]  Current: [");
                    for (int j = 0; j < segment_count; j++)
                        fprintf(stderr, "%s%d", j ? ", " : "", segment_counts[j]);
                    fprintf(stderr, "]\n");
                    tunneling_detected = true;
                    break;
                }
            }
            if (tunneling_detected) break;  // Exit experiment loop
        }

        if (cli_auto_piston_step && !piston_auto_triggered) {
            trigger_piston_step_with_protocol();
            piston_auto_triggered = true;
        }

        // Track piston work after protocol trigger
        if (piston_push_tracking) {
            double Wp = piston_work_left + piston_work_right;
            double delta = Wp - piston_work_baseline;
            if (delta > piston_work_delta_max) {
                piston_work_delta_max = delta;
            }
        }
        maybe_stop_right_piston_on_work_target();

        double t = simulation_time - wall_release_time;
        double Wp = piston_work_left + piston_work_right;
        // Wall and piston positions in sigma units from the left boundary (XW1).
        // NOTE: wall_x / piston_*_x are in pixels; dividing by PIXELS_PER_SIGMA gives sigma.
        double wx_sigma[5] = {NAN, NAN, NAN, NAN, NAN};
        double wv[5] = {NAN, NAN, NAN, NAN, NAN};

        if (num_internal_walls >= 1) {
            wx_sigma[0] = ((double)wall_x - (double)XW1) / (double)PIXELS_PER_SIGMA;
            wv[0] = (double)vx_wall;
        }
        for (int j = 0; j < extra_wall_count && j < 4; ++j) {
            if (extra_wall_positions) {
                wx_sigma[j + 1] = ((double)extra_wall_positions[j] - (double)XW1) / (double)PIXELS_PER_SIGMA;
            }
            if (extra_wall_velocity) {
                wv[j + 1] = (double)extra_wall_velocity[j];
            }
        }

        const double pistonR_x_sigma = ((double)piston_right_x - (double)XW1) / (double)PIXELS_PER_SIGMA;
        const double pistonL_x_sigma = ((double)piston_left_x - (double)XW1) / (double)PIXELS_PER_SIGMA;
        const double pistonR_v = (double)vx_piston_right;
        const double pistonL_v = (double)vx_piston_left;

        // Segment count string (for debugging / validation; helps detect wrong seeding)
        char seg_buf[256] = {0};
        char seg_eta_buf[256] = {0};
        if (segment_counts && segment_count > 0) {
            size_t off = 0;
            for (int s = 0; s < segment_count && off + 32 < sizeof(seg_buf); ++s) {
                int n = snprintf(seg_buf + off, sizeof(seg_buf) - off, "%s%d",
                                 (s == 0) ? "" : ";", segment_counts[s]);
                if (n < 0) break;
                off += (size_t)n;
            }
            build_segment_eta_csv(seg_eta_buf, sizeof(seg_eta_buf));
        }

        int log_left_particles = left_particles;
        int log_right_particles = right_particles;
        if (segment_counts && segment_count > 0) {
            log_left_particles = segment_counts[0];
            log_right_particles = 0;
            for (int s = 1; s < segment_count; ++s) log_right_particles += segment_counts[s];
        }

        fprintf(elog,
                "%.6f, %.6f, %.6e, %.6e, %.6e, %d, %d, %.9e, %.9e, %.9e,"
                "%.6f, %.6e, %.6f, %.6e, %.6f, %.6e, %.6f, %.6e, %.6f, %.6e,"
                "%.6f, %.6e, %.6f, %.6e,"
                "%.6f, %.9e, %.9e, %.9e,",
                t,
                (double)wall_x,
                (double)energy_measurement.spring_force,
                (double)energy_measurement.spring_energy,
                (double)energy_measurement.energy_transferred,
                log_left_particles, log_right_particles,
                Wp,
                piston_work_delta_max,
                (double)energy_measurement.spring_energy_max,
                wx_sigma[0], wv[0],
                wx_sigma[1], wv[1],
                wx_sigma[2], wv[2],
                wx_sigma[3], wv[3],
                wx_sigma[4], wv[4],
                pistonR_x_sigma, pistonR_v,
                pistonL_x_sigma, pistonL_v,
                (double)(energy_measurement.eff_window_started ? (energy_measurement.eff_piston_stop_time - (float)wall_release_time) : NAN),
                (double)(energy_measurement.eff_window_started ? energy_measurement.spring_energy_at_stop : NAN),
                (double)(energy_measurement.eff_window_started ? energy_measurement.spring_energy_peak_window : NAN),
                (double)(energy_measurement.eff_window_started ? fmaxf(0.0f, energy_measurement.spring_energy_peak_window - energy_measurement.spring_energy_at_stop) : NAN));
        csv_put_escaped(elog, seg_buf[0] ? seg_buf : NULL);
        fputc(',', elog);
        csv_put_escaped(elog, seg_eta_buf[0] ? seg_eta_buf : NULL);
        fputc('\n', elog);
        recorded_steps++;

        if (cli_eff_stop_after_window && energy_measurement.eff_window_done) {
            if (!cli_quiet) {
                printf("⏱️ Efficiency window done (Tw=%.3g σ-time) — stopping early.\n",
                       (double)cli_eff_window_sigma);
            }
            break;
        }
    }

    free(initial_segment_counts);
    initial_segment_counts = NULL;

    fclose(elog);

    if (tunneling_detected) {
        fprintf(stderr, "ABORTING: particle teleportation through wall detected. "
                        "Results for this run are invalid.\n");
        exit(1);
    }

    if (!cli_quiet) {
        printf("✅ Energy transfer experiment done. Output → %s\n", trace_path);
    }

	    // Append per-run summary for quick parameter sweeps.
	    {
	        char summary_path[512];
	        if (cli_energy_transfer_summary_set && cli_energy_transfer_summary_path[0] != '\0') {
	            resolve_cli_output_path(summary_path, sizeof(summary_path), cli_energy_transfer_summary_path);
	        } else {
	            snprintf(summary_path, sizeof(summary_path), "%s/energy_transfer_runs_%dwalls.csv",
	                     g_energy_transfer_dir, num_internal_walls);
	        }

		            const char *header =
		                "timestamp,mode,sim_mode,steps_after_release,wall_hold_steps,dt_sigma,"
		                "L0,height,particles_total,radius_sigma,eta_input,eta_nominal,"
		                "eta_particles_region,"
		                "num_walls,wall_positions_cli,wall_mass_factors_cli,particles_boxes_cli,"
		                "wall_thickness_sigma,wall_thickness_vis_sigma,"
		                "kbt1,kBT,kB_effective,temperature_runtime,"
		                "spring_k,spring_eq_sigma,spring_eq_sigma_effective,"
		                "piston_auto,piston_protocol,piston_target_sigma,piston_speed_px,"
		                "piston_right_protocol_mode,piston_right_travel_sigma,piston_right_v_step_constant,"
		                "piston_right_v0,piston_right_v_gradient,piston_right_vmax,piston_right_duration,piston_right_sigmoid_steepness,"
		                "W_in_max,W_out_max,ratio,loss_fraction,total_loss,"
		                "springE_stop,springE_peak_win,W_out_peak_win,ratio_win,eff_window_T,piston_stop_t_rel,"
		                "trace_path,spring_peak_t_rel,spring_peak_dt,spring_peak_near_end,command,seed\n";

        ensure_csv_header_schema(summary_path, header);
        const bool need_header = (access(summary_path, F_OK) != 0);
        FILE *sf = fopen(summary_path, "a");
        if (sf) {
            if (need_header) {
                fprintf(sf, "%s", header);
            }

            // Timestamp (local time)
            time_t now = time(NULL);
            struct tm tm_now;
            localtime_r(&now, &tm_now);
            char ts[64];
            strftime(ts, sizeof(ts), "%Y-%m-%d %H:%M:%S", &tm_now);

            // Basic geometry + particle size
            const double L0 = (double)L0_UNITS;
            const double H  = (double)HEIGHT_UNITS;
            const double r_sigma = (double)PARTICLE_RADIUS / (double)PIXELS_PER_SIGMA;
            const double eta_nominal = (particles_active > 0)
                ? ((double)particles_active * M_PI * r_sigma * r_sigma) / (2.0 * L0 * H)
                : 0.0;
            const double eta_input = cli_packing_fraction_set ? (double)cli_packing_fraction : NAN;

            // Wall thicknesses
            const double t_sigma = fabs((double)wall_thickness_runtime / (double)PIXELS_PER_SIGMA);
            const double t_vis_sigma = fabs((double)wall_thickness_visual_runtime / (double)PIXELS_PER_SIGMA);

            // CLI list strings
            char pos_buf[256] = {0};
            char mass_buf[256] = {0};
            char boxes_buf[256] = {0};

            // wall positions
            if (cli_num_walls > 0) {
                size_t off = 0;
                for (int w = 0; w < cli_num_walls && off + 32 < sizeof(pos_buf); ++w) {
                    int n = snprintf(pos_buf + off, sizeof(pos_buf) - off, "%s%.6g",
                                     (w == 0) ? "" : ";", (double)cli_wall_positions[w]);
                    if (n < 0) break;
                    off += (size_t)n;
                }
            }
            // wall mass factors
            if (cli_num_walls > 0) {
                size_t off = 0;
                for (int w = 0; w < cli_num_walls && off + 32 < sizeof(mass_buf); ++w) {
                    float mf = cli_wall_masses_set ? cli_wall_masses[w] : (float)(wall_mass_runtime / PARTICLE_MASS);
                    int n = snprintf(mass_buf + off, sizeof(mass_buf) - off, "%s%.6g",
                                     (w == 0) ? "" : ";", (double)mf);
                    if (n < 0) break;
                    off += (size_t)n;
                }
            }
            // particles-boxes (absolute counts if provided)
            if (cli_particles_box_counts && cli_particles_box_counts_count > 0) {
                size_t off = 0;
                for (size_t s = 0; s < cli_particles_box_counts_count && off + 32 < sizeof(boxes_buf); ++s) {
                    int n = snprintf(boxes_buf + off, sizeof(boxes_buf) - off, "%s%d",
                                     (s == 0) ? "" : ";", cli_particles_box_counts[s]);
                    if (n < 0) break;
                    off += (size_t)n;
                }
            }

            const char *mode_name = (sim_mode == MODE_EDMD) ? "edmd"
                                 : (sim_mode == MODE_EDMD_HYBRID) ? "edmd-hybrid"
                                 : (sim_mode == MODE_RK4) ? "rk4"
                                 : "time";

            const char *proto_name = (cli_protocol == PROTOCOL_SIGMOIDAL) ? "sigmoidal"
                                 : (cli_protocol == PROTOCOL_LINEAR) ? "linear"
                                 : (cli_protocol == PROTOCOL_SINUSOIDAL) ? "sinusoidal"
                                 : (cli_protocol == PROTOCOL_OPTIMAL) ? "optimal"
                                 : "step";

            const double dt_sigma = (double)fixed_dt_runtime / (double)PIXELS_PER_SIGMA;
            const double spring_eq_sigma_eff = ((double)energy_measurement.equilibrium_position - (double)XW1) / (double)PIXELS_PER_SIGMA;
            const double piston_target_sigma = ((double)piston_step_target - (double)XW1) / (double)PIXELS_PER_SIGMA;

            // "Particle region" packing fraction: eta computed using only segments that contain particles.
            double particle_region_width = 0.0;
            if (cli_particles_box_counts && cli_particles_box_counts_count == (size_t)(num_internal_walls + 1) && cli_num_walls == num_internal_walls) {
                double left = 0.0;
                for (int s = 0; s < num_internal_walls + 1; ++s) {
                    double right = (s < num_internal_walls) ? (double)cli_wall_positions[s] : (2.0 * L0);
                    double w = right - left;
                    if (w < 0.0) w = 0.0;
                    if (cli_particles_box_counts[s] > 0) particle_region_width += w;
                    left = right;
                }
            }
            const double eta_particles_region = (particles_active > 0 && particle_region_width > 0.0)
                ? ((double)particles_active * M_PI * r_sigma * r_sigma) / (particle_region_width * H)
                : NAN;

	            const double W_in_max = piston_work_delta_max;
	            const double W_out_max = (double)energy_measurement.spring_energy_max;
	            const double ratio = (W_in_max > 0.0) ? (W_out_max / W_in_max) : NAN;
	            const double loss_fraction = (W_in_max > 0.0) ? ((W_in_max - W_out_max) / W_in_max) : NAN;
	            const double total_loss = W_in_max - W_out_max;

		            const double springE_stop = energy_measurement.eff_window_started ? (double)energy_measurement.spring_energy_at_stop : NAN;
		            const double springE_peak_win = energy_measurement.eff_window_started ? (double)energy_measurement.spring_energy_peak_window : NAN;
	            double W_out_peak_win = NAN;
	            if (energy_measurement.eff_window_started) {
	                W_out_peak_win = springE_peak_win - springE_stop;
	                if (W_out_peak_win < 0.0) W_out_peak_win = 0.0;
	            }
		            const double ratio_win = (W_in_max > 0.0 && isfinite(W_out_peak_win)) ? (W_out_peak_win / W_in_max) : NAN;
		            const double piston_stop_t_rel = energy_measurement.eff_window_started ? (double)(energy_measurement.eff_piston_stop_time - (float)wall_release_time) : NAN;
		            const double spring_peak_t_rel = energy_measurement.eff_window_started
		                ? (double)(energy_measurement.spring_energy_peak_window_time - (float)wall_release_time)
		                : NAN;
		            const double spring_peak_dt = (isfinite(piston_stop_t_rel) && isfinite(spring_peak_t_rel))
		                ? (spring_peak_t_rel - piston_stop_t_rel)
		                : NAN;
		            const int spring_peak_near_end = (energy_measurement.eff_window_started &&
		                                             (cli_eff_window_sigma > 0.0f) &&
		                                             isfinite(spring_peak_dt) &&
		                                             (spring_peak_dt >= 0.90 * (double)cli_eff_window_sigma))
		                ? 1
		                : 0;

	            const double piston_right_travel_sigma = (cli_piston_right_travel_set ? (double)cli_piston_right_travel_sigma
	                                                : (cli_piston_travel_set ? (double)cli_piston_travel_sigma
	                                                : NAN));
	            const double piston_right_v_step = cli_piston_right_step_speed_set ? (double)cli_piston_right_step_speed : NAN;
	            const double piston_right_v0 = cli_piston_right_v0_set ? (double)cli_piston_right_v0 : NAN;
	            const double piston_right_v_gradient = cli_piston_right_gradient_set ? (double)cli_piston_right_gradient : NAN;
	            const double piston_right_vmax = cli_piston_right_vmax_set ? (double)cli_piston_right_vmax
	                                         : (cli_piston_speed_set ? (double)cli_piston_speed : NAN);
	            const double piston_right_duration = cli_piston_right_duration_set ? (double)cli_piston_right_duration
	                                             : (cli_piston_duration_set ? (double)cli_piston_duration : NAN);
	            const double piston_right_sigmoid_steep = cli_piston_right_sigmoid_steepness_set ? (double)cli_piston_right_sigmoid_steepness : NAN;
	            const char *piston_right_mode_str = (cli_piston_right_mode_set && cli_piston_right_mode[0] != '\0') ? cli_piston_right_mode : proto_name;

		            fprintf(sf,
		                    "%s,%s,%s,%d,%d,%.9g,"
		                    "%.6g,%.6g,%d,%.9g,%.9g,%.9g,"
		                    "%.9g,"
	                    "%d,%s,%s,%s,"
	                    "%.9g,%.9g,"
	                    "%d,%.9g,%.9g,%.9g,"
	                    "%.9g,%.9g,%.9g,"
	                    "%d,%s,%.9g,%.9g,"
			                    "%s,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,"
			                    "%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%s,",
	                    ts,
	                    "energy_transfer",
	                    mode_name,
	                    target_steps,
	                    wall_hold_steps,
	                    dt_sigma,
	                    L0, H, particles_active, r_sigma, eta_input, eta_nominal,
	                    eta_particles_region,
	                    num_internal_walls,
	                    pos_buf, mass_buf, boxes_buf,
	                    t_sigma, t_vis_sigma,
	                    cli_force_kbt_one ? 1 : 0,
	                    (double)kBT_effective(),
	                    (double)kB_effective(),
	                    (double)temperature_runtime,
	                    (double)energy_measurement.spring_constant,
	                    (double)(cli_spring_eq_set ? cli_spring_eq_sigma : NAN),
	                    spring_eq_sigma_eff,
	                    cli_auto_piston_step ? 1 : 0,
	                    proto_name,
	                    piston_target_sigma,
	                    (double)piston_step_speed,
	                    piston_right_mode_str,
	                    piston_right_travel_sigma,
	                    piston_right_v_step,
	                    piston_right_v0,
	                    piston_right_v_gradient,
	                    piston_right_vmax,
	                    piston_right_duration,
	                    piston_right_sigmoid_steep,
	                    W_in_max, W_out_max, ratio, loss_fraction, total_loss,
		                    springE_stop, springE_peak_win, W_out_peak_win, ratio_win,
		                    (double)cli_eff_window_sigma, piston_stop_t_rel,
		                    trace_path);
		            fprintf(sf, "%.9g,%.9g,%d,", spring_peak_t_rel, spring_peak_dt, spring_peak_near_end);
		            csv_put_escaped(sf, cli_command_line ? cli_command_line : "");
		            fprintf(sf, ",%u\n", cli_seed);
            fclose(sf);
        }
    }
}

typedef struct {
    double mean_x_sigma;
    double std_x_sigma;
    double mean_x_norm;
    double std_x_norm;
    double mean_vx;
    double mean_abs_vx;
    double anisotropy;
    double near_piston_frac;
    double temperature;
    double kinetic_energy;
} SimpleBoxStateFeatures;

typedef struct {
    double time;
    double work_piston_total;
    double heat_bath_total;
    double internal_energy;
} SimpleBoxThermoSnapshot;

typedef struct {
    int ready_now;
    int layout_relaxed_now;
    double temp_rel_err;
    double anisotropy_abs;
    double mean_vx_abs;
    double init_pos_mi_bits;
    double pos_entropy_bits;
    double pos_entropy_norm;
} SimpleBoxReadyStatus;

typedef struct {
    int paused;
    int step_requested;
    int autorun_enabled;
    int autorun_steps_remaining;
    int show_help;
    int ready;
    int waiting_for_trigger;
    int finished;
    int ready_hold_count;
    int ready_hold_required;
    double temp_rel_err;
    double anisotropy_abs;
    double mean_vx_abs;
    double init_pos_mi_bits;
    double pos_entropy_norm;
    double current_target_sigma;
    double display_target_sigma;
} SimpleBoxPreviewUiState;

#define SIMPLE_BOX_LAYOUT_BINS_X 8
#define SIMPLE_BOX_LAYOUT_BINS_Y 4
#define SIMPLE_BOX_LAYOUT_LABELS (SIMPLE_BOX_LAYOUT_BINS_X * SIMPLE_BOX_LAYOUT_BINS_Y)
static unsigned short simple_box_init_layout_label[NUM_PARTICLES];
static double simple_box_init_layout_entropy_bits = 0.0;

static double simple_box_clampd(double x, double lo, double hi) {
    if (x < lo) return lo;
    if (x > hi) return hi;
    return x;
}

static double simple_box_rand_unit(void) {
    return (double)rand() / (double)RAND_MAX;
}

static double simple_box_rand_uniform(double lo, double hi) {
    return lo + (hi - lo) * simple_box_rand_unit();
}

static double simple_box_rand_gaussian(void) {
    double u1 = 0.0;
    double u2 = 0.0;
    do {
        u1 = simple_box_rand_unit();
    } while (u1 <= 1e-12);
    u2 = simple_box_rand_unit();
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

static double simple_box_bath_kbt(void) {
    const double Tbath = heatbath_enabled ? (double)heatbath_temperature : (double)temperature_runtime;
    return (double)kB_effective() * Tbath;
}

static void simple_box_measure_state_features(SimpleBoxStateFeatures *out) {
    memset(out, 0, sizeof(*out));
    const int n = (particles_active > 0) ? particles_active : 0;
    if (n <= 0) return;

    const double piston_sigma = ((double)piston_right_x - (double)XW1) / (double)PIXELS_PER_SIGMA;
    const double avail_sigma = fmax(1e-9, piston_sigma);
    const double near_cut_sigma = fmax(0.0, piston_sigma - 0.2 * avail_sigma);

    double sum_x = 0.0;
    double sum_x2 = 0.0;
    double sum_vx = 0.0;
    double sum_abs_vx = 0.0;
    double sum_vx2 = 0.0;
    double sum_vy2 = 0.0;
    int near_count = 0;

    for (int i = 0; i < n; ++i) {
        double x_sigma = ((double)X[i] - (double)XW1) / (double)PIXELS_PER_SIGMA;
        x_sigma = simple_box_clampd(x_sigma, 0.0, piston_sigma);
        sum_x += x_sigma;
        sum_x2 += x_sigma * x_sigma;
        sum_vx += (double)Vx[i];
        sum_abs_vx += fabs((double)Vx[i]);
        sum_vx2 += (double)Vx[i] * (double)Vx[i];
        sum_vy2 += (double)Vy[i] * (double)Vy[i];
        if (x_sigma >= near_cut_sigma) near_count++;
    }

    out->mean_x_sigma = sum_x / (double)n;
    {
        double var_x = (sum_x2 / (double)n) - out->mean_x_sigma * out->mean_x_sigma;
        if (var_x < 0.0) var_x = 0.0;
        out->std_x_sigma = sqrt(var_x);
    }
    out->mean_x_norm = out->mean_x_sigma / avail_sigma;
    out->std_x_norm = out->std_x_sigma / avail_sigma;
    out->mean_vx = sum_vx / (double)n;
    out->mean_abs_vx = sum_abs_vx / (double)n;
    {
        const double denom = sum_vx2 + sum_vy2;
        out->anisotropy = (denom > 1e-12) ? ((sum_vx2 - sum_vy2) / denom) : 0.0;
    }
    out->near_piston_frac = (double)near_count / (double)n;
    out->temperature = (double)compute_measured_temperature_from_ke();
    out->kinetic_energy = kinetic_energy();
}

static SimpleBoxThermoSnapshot simple_box_thermo_snapshot(void) {
    SimpleBoxThermoSnapshot snap;
    snap.time = simulation_time;
    snap.work_piston_total = piston_work_left + piston_work_right;
    snap.heat_bath_total = g_edmd ? edmd_heat_bath(g_edmd) : NAN;
    snap.internal_energy = kinetic_energy();
    return snap;
}

static void simple_box_write_trace_row(FILE *trace_f,
                                       long long *trace_index,
                                       int cycle_index,
                                       const char *phase_label,
                                       double target_sigma)
{
    SimpleBoxStateFeatures feat;
    double init_mi_bits = 0.0, pos_entropy_bits = 0.0, pos_entropy_norm = 0.0;
    simple_box_measure_state_features(&feat);
    simple_box_layout_metrics(&init_mi_bits, &pos_entropy_bits, &pos_entropy_norm);
    const SimpleBoxThermoSnapshot snap = simple_box_thermo_snapshot();
    fprintf(trace_f,
            "%lld,%d,%s,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g\n",
            ++(*trace_index),
            cycle_index,
            phase_label,
            snap.time,
            target_sigma,
            ((double)piston_right_x - (double)XW1) / (double)PIXELS_PER_SIGMA,
            (double)vx_piston_right,
            snap.work_piston_total,
            snap.heat_bath_total,
            feat.kinetic_energy,
            feat.temperature,
            feat.mean_x_sigma,
            feat.std_x_sigma,
            feat.mean_x_norm,
            feat.std_x_norm,
            feat.anisotropy,
            feat.near_piston_frac,
            init_mi_bits,
            pos_entropy_bits,
            pos_entropy_norm);
}

static void simple_box_write_state_row(FILE *state_f,
                                       long long *state_index,
                                       int cycle_index,
                                       double x_t_sigma,
                                       double x_next_sigma)
{
    SimpleBoxStateFeatures feat;
    const SimpleBoxThermoSnapshot snap_state = simple_box_thermo_snapshot();
    double init_mi_bits = 0.0, pos_entropy_bits = 0.0, pos_entropy_norm = 0.0;
    simple_box_measure_state_features(&feat);
    simple_box_layout_metrics(&init_mi_bits, &pos_entropy_bits, &pos_entropy_norm);

    fprintf(state_f,
            "%lld,%d,%.9g,%.9g,%.9g,%s,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g\n",
            (*state_index)++,
            cycle_index,
            snap_state.time,
            x_t_sigma,
            x_next_sigma,
            cli_simple_box_target_protocol,
            snap_state.work_piston_total,
            snap_state.heat_bath_total,
            snap_state.internal_energy,
            feat.temperature,
            feat.mean_x_sigma,
            feat.std_x_sigma,
            feat.mean_x_norm,
            feat.std_x_norm,
            feat.mean_vx,
            feat.mean_abs_vx,
            feat.anisotropy,
            feat.near_piston_frac,
            init_mi_bits,
            pos_entropy_bits,
            pos_entropy_norm);
}

static void simple_box_write_step_row(FILE *step_f,
                                      long long *step_index,
                                      double t_start,
                                      double t_drive_end,
                                      int relax_steps_done,
                                      double x_prev_sigma,
                                      double x_next_sigma,
                                      const SimpleBoxThermoSnapshot *snap_start,
                                      const SimpleBoxThermoSnapshot *snap_drive,
                                      const SimpleBoxThermoSnapshot *snap_end)
{
    const double kbt_bath = simple_box_bath_kbt();
    const double deltaF_eq = (kbt_bath > 0.0 && x_prev_sigma > 1e-12 && x_next_sigma > 1e-12)
        ? ((double)particles_active * kbt_bath * log(x_prev_sigma / x_next_sigma))
        : NAN;
    const double W_drive = snap_drive->work_piston_total - snap_start->work_piston_total;
    const double W_step = snap_end->work_piston_total - snap_start->work_piston_total;
    const double Q_drive = snap_drive->heat_bath_total - snap_start->heat_bath_total;
    const double Q_relax = snap_end->heat_bath_total - snap_drive->heat_bath_total;
    const double Q_step = snap_end->heat_bath_total - snap_start->heat_bath_total;
    const double deltaE = snap_end->internal_energy - snap_start->internal_energy;
    const double betaW = (kbt_bath > 1e-12) ? (W_step / kbt_bath) : NAN;
    const double exp_minus_beta_W = isfinite(betaW) ? exp(-betaW) : NAN;
    const double W_diss_eq = isfinite(deltaF_eq) ? (W_step - deltaF_eq) : NAN;

    fprintf(step_f,
            "%lld,%s,%.9g,%.9g,%.9g,%d,%.9g,%.9g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g\n",
            (*step_index)++,
            cli_simple_box_target_protocol,
            t_start,
            t_drive_end,
            snap_end->time,
            relax_steps_done,
            x_prev_sigma,
            x_next_sigma,
            snap_start->work_piston_total,
            snap_drive->work_piston_total,
            snap_end->work_piston_total,
            W_drive,
            W_step,
            snap_start->heat_bath_total,
            snap_drive->heat_bath_total,
            snap_end->heat_bath_total,
            Q_drive,
            Q_relax,
            Q_step,
            snap_start->internal_energy,
            snap_drive->internal_energy,
            snap_end->internal_energy,
            deltaE,
            kbt_bath,
            deltaF_eq,
            W_diss_eq,
            exp_minus_beta_W);
}

static int simple_box_ready_hold_required(void) {
    int hold = cli_simple_box_relax_steps / 8;
    if (hold < 4) hold = 4;
    if (hold > 24) hold = 24;
    return hold;
}

static int simple_box_layout_label_from_xy(double x_px, double y_px) {
    double fx = ((x_px - (double)XW1) / fmax(1.0, (double)(XW2 - XW1)));
    double fy = ((y_px - (double)YW1) / fmax(1.0, (double)(YW2 - YW1)));
    int bx, by;
    if (fx < 0.0) fx = 0.0;
    if (fx >= 1.0) fx = 0.999999;
    if (fy < 0.0) fy = 0.0;
    if (fy >= 1.0) fy = 0.999999;
    bx = (int)floor(fx * SIMPLE_BOX_LAYOUT_BINS_X);
    by = (int)floor(fy * SIMPLE_BOX_LAYOUT_BINS_Y);
    if (bx < 0) bx = 0;
    if (bx >= SIMPLE_BOX_LAYOUT_BINS_X) bx = SIMPLE_BOX_LAYOUT_BINS_X - 1;
    if (by < 0) by = 0;
    if (by >= SIMPLE_BOX_LAYOUT_BINS_Y) by = SIMPLE_BOX_LAYOUT_BINS_Y - 1;
    return by * SIMPLE_BOX_LAYOUT_BINS_X + bx;
}

static double simple_box_discrete_entropy_bits_from_counts(const int *counts, int nlabels, int total) {
    double H = 0.0;
    if (!counts || nlabels <= 0 || total <= 0) return 0.0;
    for (int i = 0; i < nlabels; ++i) {
        if (counts[i] > 0) {
            const double p = (double)counts[i] / (double)total;
            H -= p * (log(p) / log(2.0));
        }
    }
    return H;
}

static void simple_box_capture_initial_layout(void) {
    int counts[SIMPLE_BOX_LAYOUT_LABELS];
    memset(counts, 0, sizeof(counts));
    for (int i = 0; i < particles_active; ++i) {
        const int lab = simple_box_layout_label_from_xy((double)X[i], (double)Y[i]);
        simple_box_init_layout_label[i] = (unsigned short)lab;
        counts[lab] += 1;
    }
    simple_box_init_layout_entropy_bits = simple_box_discrete_entropy_bits_from_counts(
        counts, SIMPLE_BOX_LAYOUT_LABELS, particles_active
    );
}

static void simple_box_layout_metrics(double *mi_bits, double *H_bits, double *H_norm) {
    int counts0[SIMPLE_BOX_LAYOUT_LABELS];
    int counts1[SIMPLE_BOX_LAYOUT_LABELS];
    int joint[SIMPLE_BOX_LAYOUT_LABELS][SIMPLE_BOX_LAYOUT_LABELS];
    double I = 0.0;
    double H = 0.0;
    const int total = (particles_active > 0) ? particles_active : 0;

    memset(counts0, 0, sizeof(counts0));
    memset(counts1, 0, sizeof(counts1));
    memset(joint, 0, sizeof(joint));

    if (total <= 0) {
        if (mi_bits) *mi_bits = 0.0;
        if (H_bits) *H_bits = 0.0;
        if (H_norm) *H_norm = 0.0;
        return;
    }

    for (int i = 0; i < total; ++i) {
        const int a = (int)simple_box_init_layout_label[i];
        const int b = simple_box_layout_label_from_xy((double)X[i], (double)Y[i]);
        if (a >= 0 && a < SIMPLE_BOX_LAYOUT_LABELS && b >= 0 && b < SIMPLE_BOX_LAYOUT_LABELS) {
            counts0[a] += 1;
            counts1[b] += 1;
            joint[a][b] += 1;
        }
    }

    H = simple_box_discrete_entropy_bits_from_counts(counts1, SIMPLE_BOX_LAYOUT_LABELS, total);
    for (int a = 0; a < SIMPLE_BOX_LAYOUT_LABELS; ++a) {
        for (int b = 0; b < SIMPLE_BOX_LAYOUT_LABELS; ++b) {
            if (joint[a][b] > 0 && counts0[a] > 0 && counts1[b] > 0) {
                const double p_ab = (double)joint[a][b] / (double)total;
                const double p_a = (double)counts0[a] / (double)total;
                const double p_b = (double)counts1[b] / (double)total;
                I += p_ab * (log(p_ab / (p_a * p_b)) / log(2.0));
            }
        }
    }

    if (mi_bits) *mi_bits = I;
    if (H_bits) *H_bits = H;
    if (H_norm) {
        const double Hmax = log((double)SIMPLE_BOX_LAYOUT_LABELS) / log(2.0);
        *H_norm = (Hmax > 1e-12) ? (H / Hmax) : 0.0;
    }
}

static SimpleBoxReadyStatus simple_box_ready_status(void) {
    SimpleBoxReadyStatus rs;
    SimpleBoxStateFeatures feat;
    const double temp_tol = fmax(0.12, 3.0 * (double)stability_window_percent);
    const double anisotropy_tol = 0.20;
    const double kbt_bath = simple_box_bath_kbt();
    const double v_scale = sqrt(fmax(1e-12, kbt_bath / (double)PARTICLE_MASS));
    const double mean_vx_tol = 0.35 * v_scale;
    const double target_T = heatbath_enabled ? (double)heatbath_temperature : (double)temperature_runtime;
    const double init_mi_tol = fmax(1.50, 0.50 * simple_box_init_layout_entropy_bits);
    const double pos_entropy_norm_tol = 0.70;

    simple_box_measure_state_features(&feat);
    simple_box_layout_metrics(&rs.init_pos_mi_bits, &rs.pos_entropy_bits, &rs.pos_entropy_norm);
    rs.temp_rel_err = (target_T > 1e-12) ? fabs(feat.temperature - target_T) / target_T : 0.0;
    rs.anisotropy_abs = fabs(feat.anisotropy);
    rs.mean_vx_abs = fabs(feat.mean_vx);
    rs.layout_relaxed_now = (rs.init_pos_mi_bits <= init_mi_tol &&
                             rs.pos_entropy_norm >= pos_entropy_norm_tol);
    rs.ready_now = (!piston_step_active &&
                    fabs((double)vx_piston_right) < 1e-9 &&
                    rs.temp_rel_err <= temp_tol &&
                    rs.anisotropy_abs <= anisotropy_tol &&
                    rs.mean_vx_abs <= mean_vx_tol);
    return rs;
}

static double simple_box_pick_next_target_sigma(int step_index,
                                                double current_sigma,
                                                double min_sigma,
                                                double max_sigma,
                                                double *corr_dx_sigma)
{
    const double span = fmax(1e-9, max_sigma - min_sigma);

    if (strcmp(cli_simple_box_target_protocol, "iid") == 0) {
        return simple_box_rand_uniform(min_sigma, max_sigma);
    }

    if (strcmp(cli_simple_box_target_protocol, "periodic") == 0) {
        const int period = (cli_simple_box_period < 2) ? 2 : cli_simple_box_period;
        const double phase = fmod((double)(step_index + 1), (double)period) / (double)period;
        const double tri = (phase < 0.5) ? (1.0 - 2.0 * phase) : (2.0 * phase - 1.0);
        return min_sigma + tri * span;
    }

    {
        double dx = cli_simple_box_corr_alpha * (*corr_dx_sigma);
        dx += cli_simple_box_noise_sigma * simple_box_rand_gaussian();
        if (fabs(dx) < 1e-6) {
            const double kick = fmax(0.25, 0.12 * span);
            dx = (simple_box_rand_unit() < 0.5 ? -kick : kick);
        }

        double x_next = current_sigma + dx;
        while (x_next < min_sigma || x_next > max_sigma) {
            if (x_next < min_sigma) {
                x_next = min_sigma + (min_sigma - x_next);
                dx = -dx;
            }
            if (x_next > max_sigma) {
                x_next = max_sigma - (x_next - max_sigma);
                dx = -dx;
            }
        }
        *corr_dx_sigma = dx;
        return simple_box_clampd(x_next, min_sigma, max_sigma);
    }
}

static void simple_box_initialize_edmd_state(void) {
    if (g_edmd) {
        edmd_destroy(g_edmd);
        g_edmd = NULL;
    }

    EDMD_Params prm = {0};
    prm.boxW = (double)(XW2 - XW1);
    prm.boxH = (double)(YW2 - YW1);
    prm.radius = (double)PARTICLE_RADIUS;
    prm.N = particles_active;
    prm.cell_size = 0.0;
    prm.divider_count = 0;
    prm.heatbath_enabled = heatbath_enabled;
    prm.heatbath_temperature = (double)heatbath_temperature;
    prm.thermal_wall_mode = thermal_wall_mode;
    prm.mb_overshoot_factor = (double)mb_overshoot_factor;
    prm.stability_window_percent = (double)stability_window_percent;
    prm.particle_mass = (double)PARTICLE_MASS;
    prm.kB = (double)kB_effective();
    /* Simple prediction box: use ideal-gas limit by default.
       Keep EDMD for exact wall/piston event handling, but disable particle-particle collisions. */
    prm.pp_collisions_enabled = 0;

    g_edmd = edmd_create(&prm);
    if (!g_edmd) {
        fprintf(stderr, "EDMD create failed for simple prediction box.\n");
        exit(EXIT_FAILURE);
    }

    edmd_config_dividers(g_edmd, 0, NULL, NULL);
    edmd_config_pistons(g_edmd,
                        1, (double)((piston_left_x + 5.0f) - XW1), (double)vx_piston_left, 0.0,
                        1, (double)(piston_right_x - XW1), (double)vx_piston_right, 0.0);

    EDMD_Particle *P = (EDMD_Particle *)edmd_particles(g_edmd);
    for (int i = 0; i < particles_active; ++i) {
        double R = (double)PARTICLE_RADIUS;
        double bx = (double)(X[i] - XW1);
        double by = (double)(Y[i] - YW1);
        if (bx < R) bx = R;
        if (bx > prm.boxW - R) bx = prm.boxW - R;
        if (by < R) by = R;
        if (by > prm.boxH - R) by = prm.boxH - R;
        P[i].x = bx;
        P[i].y = by;
        P[i].vx = (double)Vx[i];
        P[i].vy = (double)Vy[i];
        P[i].coll_count = 0;
    }

    edmd_relax_pp_overlaps_szilard(g_edmd);
    edmd_reschedule_all(g_edmd);
    edmd_reset_work(g_edmd);
}

static void simple_box_run_edmd_step(double *prev_vxL,
                                     double *prev_vxR,
                                     double *prev_work_L,
                                     double *prev_work_R)
{
    const double dt_sigma = (double)fixed_dt_runtime / (double)PIXELS_PER_SIGMA;
    apply_piston_protocol(fixed_dt_runtime);

    if ((double)vx_piston_left != *prev_vxL || (double)vx_piston_right != *prev_vxR) {
        const EDMD_Params *epc = edmd_params(g_edmd);
        edmd_config_pistons(g_edmd,
                            1, epc->pistonL_x, (double)vx_piston_left, epc->pistonL_mass,
                            1, epc->pistonR_x, (double)vx_piston_right, epc->pistonR_mass);
        *prev_vxL = (double)vx_piston_left;
        *prev_vxR = (double)vx_piston_right;
    }

    edmd_reschedule_all(g_edmd);
    {
        const double t0 = edmd_time(g_edmd);
        edmd_advance_to(g_edmd, t0 + (double)fixed_dt_runtime);
    }
    edmd_divider_resolve_overlaps(g_edmd);

    {
        const EDMD_Particle *P = edmd_particles(g_edmd);
        for (int i = 0; i < particles_active; ++i) {
            X[i] = (double)XW1 + (double)P[i].x;
            Y[i] = (double)YW1 + (double)P[i].y;
            Vx[i] = (double)P[i].vx;
            Vy[i] = (double)P[i].vy;
        }
    }

    {
        const EDMD_Params *ep = edmd_params(g_edmd);
        piston_left_x = (float)(XW1 + ep->pistonL_x - 5.0);
        piston_right_x = (float)(XW1 + ep->pistonR_x);
    }

    {
        const double wL = edmd_work_pistonL(g_edmd);
        const double wR = edmd_work_pistonR(g_edmd);
        piston_work_left += (wL - *prev_work_L);
        piston_work_right += (wR - *prev_work_R);
        *prev_work_L = wL;
        *prev_work_R = wR;
    }

    if (piston_step_active &&
        ((piston_step_direction < 0.0f && piston_right_x <= piston_step_target) ||
         (piston_step_direction > 0.0f && piston_right_x >= piston_step_target))) {
        const EDMD_Params *ep = edmd_params(g_edmd);
        piston_right_x = piston_step_target;
        vx_piston_right = 0.0f;
        piston_step_active = false;
        piston_step_direction = -1.0f;
        piston_protocol_start_time = -1.0f;
        edmd_config_pistons(g_edmd,
                            1, ep->pistonL_x, (double)vx_piston_left, ep->pistonL_mass,
                            1, (double)(piston_step_target - (float)XW1), 0.0, ep->pistonR_mass);
        *prev_vxR = 0.0;
    }

    recompute_segment_stats_counts_and_temperature();
    simulation_time += dt_sigma;

    if (piston_push_tracking) {
        const double Wp = piston_work_left + piston_work_right;
        const double delta = Wp - piston_work_baseline;
        if (delta > piston_work_delta_max) {
            piston_work_delta_max = delta;
        }
    }
    maybe_stop_right_piston_on_work_target();
}

static void simple_box_write_plot_commands(const char *run_dir) {
    char path[1200];
    snprintf(path, sizeof(path), "%s/01_PLOT_COMMANDS.md", run_dir);
    FILE *f = fopen(path, "w");
    if (!f) return;
    fprintf(f, "# Plot commands\n\n");
    fprintf(f, "Run from the `hspist3/` directory:\n\n");
    fprintf(f, "```sh\n");
    fprintf(f, "cd hspist3\n");
    fprintf(f, "python3 plot_simple_box_prediction.py --trace %s/simple_box_trace.csv --states %s/simple_box_state.csv --steps %s/simple_box_step.csv --out %s/plots\n",
            run_dir, run_dir, run_dir, run_dir);
    fprintf(f, "```\n");
    fclose(f);
}

static void simple_box_build_run_dir(char *run_dir, size_t run_dir_sz) {
    if (!run_dir || run_dir_sz == 0) return;
    run_dir[0] = '\0';

    if (cli_simple_box_run_dir && cli_simple_box_run_dir[0]) {
        resolve_cli_output_path(run_dir, run_dir_sz, cli_simple_box_run_dir);
        return;
    }

    {
        time_t now = time(NULL);
        struct tm *t = localtime(&now);
        char stamp[64];
        if (t) {
            snprintf(stamp, sizeof(stamp), "simulation_%02d_%02d_%02d_%02d_%02d",
                     t->tm_mday, t->tm_mon + 1, t->tm_year % 100, t->tm_hour, t->tm_min);
        } else {
            snprintf(stamp, sizeof(stamp), "simulation_%u", (unsigned int)now);
        }

        if (cli_simple_box_run_label && cli_simple_box_run_label[0]) {
            snprintf(run_dir, run_dir_sz, "%s/%s_%s", g_simple_box_dir, cli_simple_box_run_label, stamp);
        } else {
            snprintf(run_dir, run_dir_sz, "%s/%s", g_simple_box_dir, stamp);
        }
    }
}

static int simple_box_preview_poll_events(SimpleBoxPreviewUiState *ui) {
    SDL_Event e;
    while (SDL_PollEvent(&e)) {
        if (e.type == SDL_QUIT) return 0;
        if (live_controls_handle_event(&e)) continue;
        if (e.type == SDL_KEYDOWN) {
            switch (e.key.keysym.sym) {
                case SDLK_ESCAPE:
                case SDLK_q:
                    return 0;
                case SDLK_SPACE:
                case SDLK_p:
                    if (ui) ui->paused = !ui->paused;
                    break;
                case SDLK_t:
                    if (ui && ui->waiting_for_trigger && ui->ready && !ui->finished) {
                        ui->step_requested = 1;
                    }
                    break;
                case SDLK_a:
                    if (ui && !ui->finished) {
                        ui->autorun_enabled = !ui->autorun_enabled;
                    }
                    break;
                case SDLK_g:
                    if (ui && !ui->finished) {
                        if (ui->autorun_steps_remaining < 1000000 - 20) {
                            ui->autorun_steps_remaining += 20;
                        } else {
                            ui->autorun_steps_remaining = 1000000;
                        }
                    }
                    break;
                case SDLK_h:
                    if (ui) ui->show_help = !ui->show_help;
                    break;
                case SDLK_MINUS:
                case SDLK_KP_MINUS:
                    time_scale_runtime = fmaxf(0.1f, time_scale_runtime * 0.5f);
                    break;
                case SDLK_EQUALS:
                case SDLK_PLUS:
                case SDLK_KP_PLUS:
                    time_scale_runtime = fminf(256.0f, time_scale_runtime * 2.0f);
                    break;
                case SDLK_0:
                    time_scale_runtime = 1.0f;
                    break;
                default:
                    break;
            }
        }
    }
    return 1;
}

static int simple_box_gui_substeps_per_tick(void) {
    const double speed = (double)fmaxf(0.1f, time_scale_runtime);
    if (speed <= 1.0) return 1;
    if (speed >= 256.0) return 256;
    return (int)lround(speed);
}

static int simple_box_gui_effective_render_every(void) {
    const int base = (cli_simple_box_render_every > 0) ? cli_simple_box_render_every : 1;
    const int substeps = simple_box_gui_substeps_per_tick();
    const long long scaled = (long long)base * (long long)substeps;
    if (scaled <= 1LL) return 1;
    if (scaled >= 1000000LL) return 1000000;
    return (int)scaled;
}

static void simple_box_render_preview(const char *run_dir,
                                      int cycle_index,
                                      const char *phase_label,
                                      double target_sigma,
                                      const SimpleBoxPreviewUiState *ui)
{
    if (cli_headless) return;

    draw_clear_screen();
    gui_begin_scene_view();
    draw_coordinate_system(_renderer);
    draw_simulation_boundary();
    render_particles();
    render_pistons();
    gui_end_scene_view();
    render_energy_hud(_renderer, font);
    render_live_controls(_renderer, font_small ? font_small : font);

    {
        SDL_Color white = {255, 255, 255, 255};
        SDL_Color cyan = {90, 220, 255, 255};
        SDL_Color yellow = {255, 230, 120, 255};
        SDL_Color green = {110, 255, 140, 255};
        SDL_Color red = {255, 120, 120, 255};
        char line1[256];
        char line2[320];
        char line3[320];
        char line4[384];
        char line5[384];
        const char *base = run_dir;
        const char *slash = run_dir ? strrchr(run_dir, '/') : NULL;
        const int paused_now = ui ? ui->paused : 0;
        const int ready_now = ui ? ui->ready : 0;
        const int waiting_now = ui ? ui->waiting_for_trigger : 0;
        const int finished_now = ui ? ui->finished : 0;
        const int autorun_now = ui ? ui->autorun_enabled : 0;
        const int queued_steps = ui ? ui->autorun_steps_remaining : 0;
        if (slash && slash[1]) base = slash + 1;

        snprintf(line1, sizeof(line1), "Simple Box | run=%s | cycle=%d/%d | phase=%s%s",
                 base, cycle_index, cli_simple_box_cycles, phase_label ? phase_label : "?",
                 paused_now ? " | PAUSED" : "");
        snprintf(line2, sizeof(line2), "t=%.3f  piston=%.3f sigma  target=%.3f sigma  T=%.3f  protocol=%s",
                 (double)simulation_time,
                 ((double)piston_right_x - (double)XW1) / (double)PIXELS_PER_SIGMA,
                 target_sigma,
                 (double)compute_measured_temperature_from_ke(),
                 cli_simple_box_target_protocol);
        if (finished_now) {
            snprintf(line3, sizeof(line3), "READY final state saved | press q/esc to close");
        } else if (waiting_now && ready_now) {
            snprintf(line3, sizeof(line3), "READY | t=1 step  a=%s  g=+20 steps  queued=%d | hold=%d/%d  I0=%.2fb  Hxy=%.2f",
                     autorun_now ? "AUTO OFF" : "AUTO ON",
                     queued_steps,
                     ui ? ui->ready_hold_count : 0,
                     ui ? ui->ready_hold_required : 0,
                     ui ? ui->init_pos_mi_bits : NAN,
                     ui ? ui->pos_entropy_norm : NAN);
        } else {
            snprintf(line3, sizeof(line3), "Busy | auto=%s queued=%d | dT/T=%.3f  anis=%.3f  |mean vx|=%.3f  I0=%.2fb  Hxy=%.2f",
                     autorun_now ? "on" : "off",
                     queued_steps,
                     ui ? ui->temp_rel_err : NAN,
                     ui ? ui->anisotropy_abs : NAN,
                     ui ? ui->mean_vx_abs : NAN,
                     ui ? ui->init_pos_mi_bits : NAN,
                     ui ? ui->pos_entropy_norm : NAN);
        }
        snprintf(line4, sizeof(line4), "Keys: t=1 step  a=toggle autorun  g=queue 20 steps  +/-/0 speed  space/p pause  h help  q/esc quit | x%.2f",
                 (double)time_scale_runtime);
        snprintf(line5, sizeof(line5), "Set --simple-box-cycles=20 if you want a 20-step protocol end-to-end.");

        draw_text(_renderer, font, line1, XW1 + 8, 8, cyan);
        draw_text(_renderer, font, line2, XW1 + 8, 30, white);
        draw_text(_renderer, font_small ? font_small : font, line3, XW1 + 8, 52, ready_now ? green : yellow);
        if (ui && ui->show_help) {
            draw_text(_renderer, font_small ? font_small : font, line4, XW1 + 8, 74, yellow);
            draw_text(_renderer, font_small ? font_small : font, line5, XW1 + 8, 96, white);
            draw_text(_renderer, font_small ? font_small : font,
                      "Simple box uses local controls here; other experiments keep their own HUD key legend.",
                      XW1 + 8, 118, white);
        } else {
            draw_text(_renderer, font_small ? font_small : font, "Press h to toggle key help", XW1 + 8, 74, red);
        }
    }

    SDL_RenderPresent(_renderer);
    if (cli_simple_box_gui_delay_ms > 0) {
        double speed = (double)time_scale_runtime;
        double scaled_delay = (speed > 1e-9) ? ((double)cli_simple_box_gui_delay_ms / speed)
                                             : (double)cli_simple_box_gui_delay_ms;
        if (scaled_delay < 0.0) scaled_delay = 0.0;
        SDL_Delay((Uint32)llround(scaled_delay));
    }
}

static void run_simple_prediction_box_experiment(void) {
    mkdir_p(g_simple_box_dir);
    write_command_md_to_dir(g_simple_box_dir, g_main_argc, g_main_argv);

    if (sim_mode != MODE_EDMD) {
        printf("⚠️ Simple prediction box uses EDMD for clean piston-work / heat-bath accounting. Switching to EDMD.\n");
        sim_mode = MODE_EDMD;
    }

    if (!heatbath_enabled) {
        heatbath_enabled = 1;
        heatbath_temperature = (temperature_runtime > 0.0f) ? temperature_runtime : 1.0f;
    }

    cli_requested_walls = 0;
    cli_num_walls = 0;
    num_internal_walls = 0;
    cli_left_empty = false;
    wall_enabled = 0;

    initialize_simulation_dimensions();
    initialize_time_scale();
    initialize_simulation_parameters();
    update_dt_runtime_for_temperature();
    initialize_simulation();

    piston_left_x = (float)XW1 - 5.0f;
    piston_right_x = (float)XW2;
    vx_piston_left = 0.0f;
    vx_piston_right = 0.0f;
    piston_step_active = false;
    piston_step_speed = 0.0f;
    piston_step_direction = -1.0f;
    piston_protocol_start_time = -1.0f;
    piston_work_left = 0.0;
    piston_work_right = 0.0;
    piston_push_tracking = false;
    piston_work_baseline = 0.0;
    piston_work_delta_max = 0.0;
    simulation_time = 0.0;
    steps_elapsed = 0;
    wall_release_time = -1.0;
    missed_collision_events = 0;
    worst_penetration_observed = 0.0;
    {
        char run_dir[1024];
        int running = 1;
        int preview_step_counter = 0;
        int preview_cycle = -1;
        const char *preview_phase = "init";
        double preview_target_sigma = NAN;

        simple_box_build_run_dir(run_dir, sizeof(run_dir));
        mkdir_p(run_dir);
        write_command_md_to_dir(run_dir, g_main_argc, g_main_argv);
        write_latest_txt(g_simple_box_dir, run_dir);
        simple_box_write_plot_commands(run_dir);
        if (!cli_quiet) {
            printf("📁 Simple-box run dir: %s\n", run_dir);
        }

        {
            char rp[1200];
            snprintf(rp, sizeof(rp), "%s/run_params.json", run_dir);
            write_run_params_json_to(rp);
        }

        {
            char trace_path[1200];
            char state_path[1200];
            char step_path[1200];
            snprintf(trace_path, sizeof(trace_path), "%s/simple_box_trace.csv", run_dir);
            snprintf(state_path, sizeof(state_path), "%s/simple_box_state.csv", run_dir);
            snprintf(step_path, sizeof(step_path), "%s/simple_box_step.csv", run_dir);

            FILE *trace_f = fopen(trace_path, "w");
            FILE *state_f = fopen(state_path, "w");
            FILE *step_f = fopen(step_path, "w");
            if (!trace_f || !state_f || !step_f) {
                fprintf(stderr, "❌ Failed to open simple-box output CSVs in %s\n", run_dir);
                if (trace_f) fclose(trace_f);
                if (state_f) fclose(state_f);
                if (step_f) fclose(step_f);
                return;
            }

            fprintf(trace_f,
                    "LogIndex,Cycle,Phase,Time,X_target_sigma,Piston_x_sigma,Piston_v_sigma_per_time,"
                    "W_piston_total,Q_bath_total,E_internal,Temperature,MeanX_sigma,StdX_sigma,"
                    "MeanX_norm,StdX_norm,Anisotropy,NearPistonFrac,InitLayoutMI_bits,PosEntropy_bits,PosEntropyNorm\n");
            fprintf(state_f,
                    "StateIndex,Cycle,Time,X_t_sigma,X_next_sigma,ProtocolKind,W_piston_total,Q_bath_total,E_internal,"
                    "Temperature,MeanX_sigma,StdX_sigma,MeanX_norm,StdX_norm,MeanVx,MeanAbsVx,Anisotropy,NearPistonFrac,"
                    "InitLayoutMI_bits,PosEntropy_bits,PosEntropyNorm\n");
            fprintf(step_f,
                    "StepIndex,ProtocolKind,TimeStart,TimeDriveEnd,TimeEnd,RelaxSteps,X_prev_sigma,X_next_sigma,"
                    "W_piston_start,W_piston_drive_end,W_piston_end,W_drive,W_step_total,"
                    "Q_bath_start,Q_bath_drive_end,Q_bath_end,Q_drive,Q_relax,Q_step_total,"
                    "E_start,E_drive_end,E_end,DeltaE_step,kBT_bath,DeltaF_eq,W_diss_eq,exp_minus_beta_W\n");

            simple_box_initialize_edmd_state();

            double prev_vxL = 0.0;
            double prev_vxR = 0.0;
            double prev_work_L = 0.0;
            double prev_work_R = 0.0;
            long long trace_index = 0;
            long long state_index = 0;
            long long step_index = 0;
            double next_trace_time = 0.0;
            double corr_dx_sigma = 0.0;
            const double box_width_sigma = ((double)XW2 - (double)XW1) / (double)PIXELS_PER_SIGMA;
            double min_target_sigma = cli_simple_box_target_min_sigma;
            double max_target_sigma = cli_simple_box_target_max_sigma;
            double first_protocol_target_sigma = cli_simple_box_initial_target_sigma;
            const double particle_radius_sigma = (double)PARTICLE_RADIUS / (double)PIXELS_PER_SIGMA;
            const double hard_min_sigma = fmax(2.0 * particle_radius_sigma + 0.25, 0.5);
            const double full_equil_target_sigma = box_width_sigma;
            const int initial_equil_step_cap = cli_simple_box_equil_steps + ((cli_simple_box_equil_steps > 0) ? (8 * cli_simple_box_equil_steps) : 2000);
            SimpleBoxPreviewUiState ui;
            memset(&ui, 0, sizeof(ui));
            ui.show_help = 1;
            ui.ready_hold_required = simple_box_ready_hold_required();

            if (!(max_target_sigma > 0.0)) {
                const double safe_default_max = fmin(box_width_sigma - (2.0 * particle_radius_sigma + 0.5),
                                                     0.75 * box_width_sigma);
                max_target_sigma = safe_default_max;
            } else if (max_target_sigma > box_width_sigma) {
                max_target_sigma = box_width_sigma;
            }
            if (!(min_target_sigma >= 0.0)) {
                min_target_sigma = fmax(hard_min_sigma, 0.35 * max_target_sigma);
            }
            if (min_target_sigma < hard_min_sigma) min_target_sigma = hard_min_sigma;
            if (min_target_sigma > max_target_sigma - 0.25) {
                min_target_sigma = fmax(hard_min_sigma, max_target_sigma - 1.0);
            }
            if (first_protocol_target_sigma > 0.0) {
                first_protocol_target_sigma = simple_box_clampd(first_protocol_target_sigma, min_target_sigma, max_target_sigma);
            }
            simple_box_capture_initial_layout();

            if (cli_headless) {
                double current_target_sigma = full_equil_target_sigma;
                int eq_steps_done = 0;
                int ready_hold_count = 0;

                {
                    const double current_sigma = ((double)piston_right_x - (double)XW1) / (double)PIXELS_PER_SIGMA;
                    if (fabs(full_equil_target_sigma - current_sigma) > 1e-9) {
                        start_right_piston_step_to_target_sigma(full_equil_target_sigma);
                        while (running && piston_step_active) {
                            simple_box_run_edmd_step(&prev_vxL, &prev_vxR, &prev_work_L, &prev_work_R);
                            if (cli_simple_box_log_dt <= 0.0 || simulation_time + 1e-12 >= next_trace_time) {
                                simple_box_write_trace_row(trace_f, &trace_index, -1, "pre_equil_drive", full_equil_target_sigma);
                                if (cli_simple_box_log_dt > 0.0) next_trace_time = simulation_time + cli_simple_box_log_dt;
                            }
                        }
                    }
                }

                while (running) {
                    SimpleBoxReadyStatus rs;
                    simple_box_run_edmd_step(&prev_vxL, &prev_vxR, &prev_work_L, &prev_work_R);
                    eq_steps_done++;
                    rs = simple_box_ready_status();
                    ready_hold_count = (rs.ready_now && rs.layout_relaxed_now) ? (ready_hold_count + 1) : 0;
                    if (cli_simple_box_log_dt <= 0.0 || simulation_time + 1e-12 >= next_trace_time) {
                        simple_box_write_trace_row(trace_f, &trace_index, -1, "equil", full_equil_target_sigma);
                        if (cli_simple_box_log_dt > 0.0) next_trace_time = simulation_time + cli_simple_box_log_dt;
                    }
                    if (eq_steps_done >= cli_simple_box_equil_steps &&
                        ready_hold_count >= ui.ready_hold_required) {
                        break;
                    }
                    if (eq_steps_done >= initial_equil_step_cap) {
                        break;
                    }
                }

                for (int cycle = 0; running && cycle <= cli_simple_box_cycles; ++cycle) {
                    const double protocol_base_sigma = simple_box_clampd(current_target_sigma, min_target_sigma, max_target_sigma);
                    const double x_next_sigma = (cycle < cli_simple_box_cycles)
                        ? ((cycle == 0 && first_protocol_target_sigma > 0.0)
                            ? first_protocol_target_sigma
                            : simple_box_pick_next_target_sigma(cycle, protocol_base_sigma,
                                                                min_target_sigma, max_target_sigma,
                                                                &corr_dx_sigma))
                        : NAN;

                    simple_box_write_state_row(state_f, &state_index, cycle, current_target_sigma, x_next_sigma);
                    if (cycle >= cli_simple_box_cycles) break;

                    {
                        const SimpleBoxThermoSnapshot snap_start = simple_box_thermo_snapshot();
                        const double t_start = simulation_time;
                        start_right_piston_step_to_target_sigma(x_next_sigma);
                        {
                            const double dt_sigma = (double)fixed_dt_runtime / (double)PIXELS_PER_SIGMA;
                            const double vmax = fmax(1e-6, fabs((double)piston_protocol_max_speed));
                            int drive_limit = 100;
                            if (dt_sigma > 1e-12) {
                                drive_limit = 100 + (int)ceil((box_width_sigma / vmax) / dt_sigma * 20.0);
                            }
                            for (int it = 0; running && piston_step_active && it < drive_limit; ++it) {
                                simple_box_run_edmd_step(&prev_vxL, &prev_vxR, &prev_work_L, &prev_work_R);
                                if (cli_simple_box_log_dt <= 0.0 || simulation_time + 1e-12 >= next_trace_time) {
                                    simple_box_write_trace_row(trace_f, &trace_index, cycle, "drive", x_next_sigma);
                                    if (cli_simple_box_log_dt > 0.0) next_trace_time = simulation_time + cli_simple_box_log_dt;
                                }
                            }
                        }

                        {
                            const SimpleBoxThermoSnapshot snap_drive = simple_box_thermo_snapshot();
                            const double t_drive_end = simulation_time;
                            for (int relax = 0; running && relax < cli_simple_box_relax_steps; ++relax) {
                                simple_box_run_edmd_step(&prev_vxL, &prev_vxR, &prev_work_L, &prev_work_R);
                                if (cli_simple_box_log_dt <= 0.0 || simulation_time + 1e-12 >= next_trace_time) {
                                    simple_box_write_trace_row(trace_f, &trace_index, cycle, "relax", x_next_sigma);
                                    if (cli_simple_box_log_dt > 0.0) next_trace_time = simulation_time + cli_simple_box_log_dt;
                                }
                            }
                            {
                                const SimpleBoxThermoSnapshot snap_end = simple_box_thermo_snapshot();
                                simple_box_write_step_row(step_f, &step_index, t_start, t_drive_end,
                                                          cli_simple_box_relax_steps,
                                                          current_target_sigma, x_next_sigma,
                                                          &snap_start, &snap_drive, &snap_end);
                            }
                        }

                        current_target_sigma = x_next_sigma;
                    }
                }
            } else {
                enum {
                    SB_PHASE_PRE_EQUIL_DRIVE = 0,
                    SB_PHASE_EQUIL,
                    SB_PHASE_WAIT_TRIGGER,
                    SB_PHASE_DRIVE,
                    SB_PHASE_RELAX,
                    SB_PHASE_DONE
                };
                int sb_phase = SB_PHASE_EQUIL;
                int cycle = 0;
                int eq_steps_done = 0;
                int relax_steps_done = 0;
                int ready_hold_count = 0;
                int final_state_written = 0;
                double current_target_sigma = full_equil_target_sigma;
                double pending_target_sigma = NAN;
                const double current_sigma = ((double)piston_right_x - (double)XW1) / (double)PIXELS_PER_SIGMA;
                SimpleBoxThermoSnapshot snap_start = {0}, snap_drive = {0};
                double t_start = 0.0;
                double t_drive_end = 0.0;

                if (fabs(full_equil_target_sigma - current_sigma) > 1e-9) {
                    start_right_piston_step_to_target_sigma(full_equil_target_sigma);
                    sb_phase = SB_PHASE_PRE_EQUIL_DRIVE;
                }

                while (running) {
                    SimpleBoxReadyStatus rs = simple_box_ready_status();
                    const int effective_render_every = simple_box_gui_effective_render_every();
                    ui.ready = rs.ready_now;
                    ui.temp_rel_err = rs.temp_rel_err;
                    ui.anisotropy_abs = rs.anisotropy_abs;
                    ui.mean_vx_abs = rs.mean_vx_abs;
                    ui.init_pos_mi_bits = rs.init_pos_mi_bits;
                    ui.pos_entropy_norm = rs.pos_entropy_norm;
                    ui.ready_hold_count = ready_hold_count;
                    ui.waiting_for_trigger = (sb_phase == SB_PHASE_WAIT_TRIGGER);
                    ui.finished = (sb_phase == SB_PHASE_DONE);
                    ui.current_target_sigma = current_target_sigma;
                    ui.display_target_sigma = isfinite(pending_target_sigma) ? pending_target_sigma : current_target_sigma;
                    preview_cycle = cycle;

                    switch (sb_phase) {
                        case SB_PHASE_PRE_EQUIL_DRIVE:
                            preview_phase = "pre_equil_drive";
                            preview_target_sigma = full_equil_target_sigma;
                            break;
                        case SB_PHASE_EQUIL:
                            preview_phase = "equil";
                            preview_target_sigma = current_target_sigma;
                            break;
                        case SB_PHASE_WAIT_TRIGGER:
                            preview_phase = "ready";
                            preview_target_sigma = current_target_sigma;
                            break;
                        case SB_PHASE_DRIVE:
                            preview_phase = "drive";
                            preview_target_sigma = pending_target_sigma;
                            break;
                        case SB_PHASE_RELAX:
                            preview_phase = "relax";
                            preview_target_sigma = pending_target_sigma;
                            break;
                        default:
                            preview_phase = "done";
                            preview_target_sigma = current_target_sigma;
                            break;
                    }

                    if (!simple_box_preview_poll_events(&ui)) { running = 0; break; }
                    while (running && ui.paused) {
                        simple_box_render_preview(run_dir, preview_cycle, preview_phase, preview_target_sigma, &ui);
                        if (!simple_box_preview_poll_events(&ui)) { running = 0; break; }
                    }
                    if (!running) break;

                    if (sb_phase == SB_PHASE_WAIT_TRIGGER || sb_phase == SB_PHASE_DONE) {
                        simple_box_render_preview(run_dir, preview_cycle, preview_phase, preview_target_sigma, &ui);
                    }

                    switch (sb_phase) {
                        case SB_PHASE_PRE_EQUIL_DRIVE:
                            for (int sub = 0; running && sub < simple_box_gui_substeps_per_tick(); ++sub) {
                                simple_box_run_edmd_step(&prev_vxL, &prev_vxR, &prev_work_L, &prev_work_R);
                                if (cli_simple_box_log_dt <= 0.0 || simulation_time + 1e-12 >= next_trace_time) {
                                    simple_box_write_trace_row(trace_f, &trace_index, -1, "pre_equil_drive", full_equil_target_sigma);
                                    if (cli_simple_box_log_dt > 0.0) next_trace_time = simulation_time + cli_simple_box_log_dt;
                                }
                                if ((preview_step_counter++ % effective_render_every) == 0) {
                                    simple_box_render_preview(run_dir, preview_cycle, preview_phase, preview_target_sigma, &ui);
                                }
                                if (!piston_step_active) {
                                    sb_phase = SB_PHASE_EQUIL;
                                    eq_steps_done = 0;
                                    ready_hold_count = 0;
                                    break;
                                }
                            }
                            break;

                        case SB_PHASE_EQUIL:
                            for (int sub = 0; running && sub < simple_box_gui_substeps_per_tick(); ++sub) {
                                simple_box_run_edmd_step(&prev_vxL, &prev_vxR, &prev_work_L, &prev_work_R);
                                eq_steps_done++;
                                rs = simple_box_ready_status();
                                ready_hold_count = (rs.ready_now && rs.layout_relaxed_now) ? (ready_hold_count + 1) : 0;
                                if (cli_simple_box_log_dt <= 0.0 || simulation_time + 1e-12 >= next_trace_time) {
                                    simple_box_write_trace_row(trace_f, &trace_index, -1, "equil", current_target_sigma);
                                    if (cli_simple_box_log_dt > 0.0) next_trace_time = simulation_time + cli_simple_box_log_dt;
                                }
                                if ((preview_step_counter++ % effective_render_every) == 0) {
                                    ui.ready = rs.ready_now;
                                    ui.temp_rel_err = rs.temp_rel_err;
                                    ui.anisotropy_abs = rs.anisotropy_abs;
                                    ui.mean_vx_abs = rs.mean_vx_abs;
                                    ui.init_pos_mi_bits = rs.init_pos_mi_bits;
                                    ui.pos_entropy_norm = rs.pos_entropy_norm;
                                    ui.ready_hold_count = ready_hold_count;
                                    simple_box_render_preview(run_dir, preview_cycle, preview_phase, preview_target_sigma, &ui);
                                }
                                if (eq_steps_done >= cli_simple_box_equil_steps &&
                                    ready_hold_count >= ui.ready_hold_required) {
                                    sb_phase = SB_PHASE_WAIT_TRIGGER;
                                    break;
                                }
                                if (eq_steps_done >= initial_equil_step_cap) {
                                    sb_phase = SB_PHASE_WAIT_TRIGGER;
                                    break;
                                }
                            }
                            break;

                        case SB_PHASE_WAIT_TRIGGER:
                            if (cycle >= cli_simple_box_cycles) {
                                if (!final_state_written) {
                                    simple_box_write_state_row(state_f, &state_index, cycle, current_target_sigma, NAN);
                                    final_state_written = 1;
                                }
                                sb_phase = SB_PHASE_DONE;
                                ui.finished = 1;
                                ui.autorun_enabled = 0;
                                ui.autorun_steps_remaining = 0;
                            } else {
                                if (!ui.step_requested &&
                                    (ui.autorun_enabled || ui.autorun_steps_remaining > 0)) {
                                    ui.step_requested = 1;
                                }
                                if (ui.step_requested) {
                                    const double protocol_base_sigma = simple_box_clampd(current_target_sigma, min_target_sigma, max_target_sigma);
                                    ui.step_requested = 0;
                                    if (ui.autorun_steps_remaining > 0) {
                                        ui.autorun_steps_remaining--;
                                    }
                                    pending_target_sigma = (cycle == 0 && first_protocol_target_sigma > 0.0)
                                        ? first_protocol_target_sigma
                                        : simple_box_pick_next_target_sigma(cycle, protocol_base_sigma,
                                                                            min_target_sigma, max_target_sigma,
                                                                            &corr_dx_sigma);
                                    simple_box_write_state_row(state_f, &state_index, cycle, current_target_sigma, pending_target_sigma);
                                    snap_start = simple_box_thermo_snapshot();
                                    t_start = simulation_time;
                                    start_right_piston_step_to_target_sigma(pending_target_sigma);
                                    sb_phase = SB_PHASE_DRIVE;
                                    ready_hold_count = 0;
                                }
                            }
                            break;

                        case SB_PHASE_DRIVE:
                            for (int sub = 0; running && sub < simple_box_gui_substeps_per_tick(); ++sub) {
                                simple_box_run_edmd_step(&prev_vxL, &prev_vxR, &prev_work_L, &prev_work_R);
                                if (cli_simple_box_log_dt <= 0.0 || simulation_time + 1e-12 >= next_trace_time) {
                                    simple_box_write_trace_row(trace_f, &trace_index, cycle, "drive", pending_target_sigma);
                                    if (cli_simple_box_log_dt > 0.0) next_trace_time = simulation_time + cli_simple_box_log_dt;
                                }
                                if ((preview_step_counter++ % effective_render_every) == 0) {
                                    simple_box_render_preview(run_dir, preview_cycle, preview_phase, preview_target_sigma, &ui);
                                }
                                if (!piston_step_active) {
                                    snap_drive = simple_box_thermo_snapshot();
                                    t_drive_end = simulation_time;
                                    relax_steps_done = 0;
                                    ready_hold_count = 0;
                                    sb_phase = SB_PHASE_RELAX;
                                    break;
                                }
                            }
                            break;

                        case SB_PHASE_RELAX:
                            for (int sub = 0; running && sub < simple_box_gui_substeps_per_tick(); ++sub) {
                                simple_box_run_edmd_step(&prev_vxL, &prev_vxR, &prev_work_L, &prev_work_R);
                                relax_steps_done++;
                                rs = simple_box_ready_status();
                                ready_hold_count = rs.ready_now ? (ready_hold_count + 1) : 0;
                                if (cli_simple_box_log_dt <= 0.0 || simulation_time + 1e-12 >= next_trace_time) {
                                    simple_box_write_trace_row(trace_f, &trace_index, cycle, "relax", pending_target_sigma);
                                    if (cli_simple_box_log_dt > 0.0) next_trace_time = simulation_time + cli_simple_box_log_dt;
                                }
                                if ((preview_step_counter++ % effective_render_every) == 0) {
                                    ui.ready = rs.ready_now;
                                    ui.temp_rel_err = rs.temp_rel_err;
                                    ui.anisotropy_abs = rs.anisotropy_abs;
                                    ui.mean_vx_abs = rs.mean_vx_abs;
                                    ui.init_pos_mi_bits = rs.init_pos_mi_bits;
                                    ui.pos_entropy_norm = rs.pos_entropy_norm;
                                    ui.ready_hold_count = ready_hold_count;
                                    simple_box_render_preview(run_dir, preview_cycle, preview_phase, preview_target_sigma, &ui);
                                }
                                if (relax_steps_done >= cli_simple_box_relax_steps &&
                                    ready_hold_count >= ui.ready_hold_required) {
                                    {
                                        const SimpleBoxThermoSnapshot snap_end = simple_box_thermo_snapshot();
                                        simple_box_write_step_row(step_f, &step_index, t_start, t_drive_end,
                                                                  relax_steps_done,
                                                                  current_target_sigma, pending_target_sigma,
                                                                  &snap_start, &snap_drive, &snap_end);
                                    }
                                    current_target_sigma = pending_target_sigma;
                                    pending_target_sigma = NAN;
                                    cycle++;
                                    ready_hold_count = 0;
                                    sb_phase = SB_PHASE_WAIT_TRIGGER;
                                    break;
                                }
                            }
                            break;

                        case SB_PHASE_DONE:
                            ui.step_requested = 0;
                            break;
                    }
                }
            }

            if (!cli_headless && (running || ui.finished)) {
                preview_phase = ui.finished ? "done" : "ready";
                simple_box_render_preview(run_dir, preview_cycle, preview_phase, preview_target_sigma, &ui);
            }

            fclose(trace_f);
            fclose(state_f);
            fclose(step_f);
        }
    }
}

// ---------------------------------------------------------------------------
// Szilard engine particle separation experiment
// ---------------------------------------------------------------------------
static void run_szilard_engine_experiment(void) {
    mkdir_p(g_szilard_dir);
    write_command_md_to_dir(g_szilard_dir, g_main_argc, g_main_argv);

    // Force time mode for now (semi-permeable walls not implemented in EDMD yet)
    if (sim_mode == MODE_EDMD || sim_mode == MODE_EDMD_HYBRID) {
        printf("⚠️ Szilard experiment currently runs in TIME mode (semipermeable walls not in EDMD yet).\n");
        sim_mode = MODE_TIME;
    }

    // Configure a single divider at L0 (midpoint), but start disabled
    cli_requested_walls = 1;
    cli_num_walls = 1;
    preset_custom_positions = true;
    preset_wall_count = 1;
    preset_wall_fraction[0] = 0.5f;

    // Override particle count if provided
    if (cli_szilard_n > 0) {
        cli_override_particles = cli_szilard_n;
    }

    // Seed all particles in left half: counts {N,0}
    int N = (cli_override_particles > 0 ? cli_override_particles : particles_active);
    if (N < 1) N = 2;
    cli_override_particles = N;
    if (cli_particles_box_counts) { free(cli_particles_box_counts); cli_particles_box_counts = NULL; }
    cli_particles_box_counts = (int*)malloc(2 * sizeof(int));
    cli_particles_box_counts_count = 2;
    cli_particles_box_counts[0] = N;
    cli_particles_box_counts[1] = 0;
    preset_custom_absolute_counts = true;
    cli_left_empty = false;

    // Initialize system
    initialize_simulation_dimensions();
    initialize_time_scale();
    initialize_simulation_parameters();
    update_dt_runtime_for_temperature();
    initialize_simulation();

    // Pistons: left fixed outside, right starts at L0 (closing right half)
    piston_left_x = (float)XW1 - 5.0f;
    piston_right_x = (float)XW1 + (float)(L0_UNITS * PIXELS_PER_SIGMA);
    vx_piston_left = 0.0f;
    vx_piston_right = 0.0f;

    // Divider disabled initially
    wall_enabled = 0;
    sz_perm_active = 0;
    sz_perm_species = 0;
    sz_perm_wall_index = primary_wall_index;
    sz_perm_dual_active = 0;

    szilard_assign_species(particles_active);

    // Logging
    time_t now = time(NULL);
    struct tm tm_now;
    localtime_r(&now, &tm_now);
    char stamp[32];
    strftime(stamp, sizeof(stamp), "%Y%m%d_%H%M%S", &tm_now);
    char trace_path[512];
    snprintf(trace_path, sizeof(trace_path), "%s/szilard_trace_%s.csv", g_szilard_dir, stamp);
    FILE *szlog = fopen(trace_path, "w");
    if (!szlog) { printf("❌ Could not open Szilard trace file.\n"); return; }
    fprintf(szlog, "Time,Phase,PistonR_x_sigma,Divider_x_sigma,PermSpecies,"
                   "Left_A,Left_B,Right_A,Right_B,MemA,MemB\n");

    // Phase control
    enum { SZ_EQUIL=0, SZ_SEPARATE=1, SZ_RECOMBINE=2, SZ_DONE=3 };
    int phase = SZ_EQUIL;

    const double dt_sigma = fixed_dt_runtime / PIXELS_PER_SIGMA;
    const double sep_duration = (cli_szilard_sep_duration > 1e-9) ? cli_szilard_sep_duration : 1.0;
    const double sep_dist_px = (double)L0_UNITS * (double)PIXELS_PER_SIGMA;
    const double sep_vx_px = sep_dist_px / (sep_duration * PIXELS_PER_SIGMA); // see time scaling in sim

    int steps_total = (num_steps > 0) ? num_steps : (cli_szilard_equil_steps + (int)(2 * sep_duration / dt_sigma) + 1000);

    simulation_time = 0.0;
    steps_elapsed = 0;

    while (phase != SZ_DONE && steps_elapsed < steps_total) {
        steps_elapsed++;

        // Update physics
        wall_x_old = wall_x;
        wall_impulse_x_accum = 0.f;
        g_substep_base_time = simulation_time;
        update_particles_with_substepping(fixed_dt_runtime, NULL, NULL, simulation_time, false, -1.0);
        update_wall_metrics();

        // Phase logic
        if (phase == SZ_EQUIL) {
            if (steps_elapsed >= cli_szilard_equil_steps) {
                // Measurement: store memory bit = species
                for (int i = 0; i < particles_active; ++i) {
                    if (sz_memory) sz_memory[i] = sz_species ? sz_species[i] : 0;
                }
                // Enable divider and permeability
                wall_enabled = 1;
                sz_perm_active = 1;
                sz_perm_dual_active = 0;
                sz_perm_species = 0; // species A passes
                phase = SZ_SEPARATE;
            }
        } else if (phase == SZ_SEPARATE) {
            vx_piston_right = (float)sep_vx_px; // move right
            piston_right_x += vx_piston_right * fixed_dt_runtime;
            if (piston_right_x >= (float)XW1 + (float)(2.0f * L0_UNITS * PIXELS_PER_SIGMA)) {
                piston_right_x = (float)XW1 + (float)(2.0f * L0_UNITS * PIXELS_PER_SIGMA);
                vx_piston_right = 0.0f;
                if (cli_szilard_perm_swap) {
                    sz_perm_species = 1 - sz_perm_species;
                }
                phase = SZ_RECOMBINE;
            }
        } else if (phase == SZ_RECOMBINE) {
            vx_piston_right = (float)(-sep_vx_px); // move left
            piston_right_x += vx_piston_right * fixed_dt_runtime;
            if (piston_right_x <= (float)XW1 + (float)(L0_UNITS * PIXELS_PER_SIGMA)) {
                piston_right_x = (float)XW1 + (float)(L0_UNITS * PIXELS_PER_SIGMA);
                vx_piston_right = 0.0f;
                wall_enabled = 0;
                sz_perm_active = 0;
                sz_perm_dual_active = 0;
                phase = SZ_DONE;
            }
        }

        // Advance time
        simulation_time += dt_sigma;

        // Logging
        int leftA, leftB, rightA, rightB, memA, memB;
        double split_x = (double)XW1 + (double)(L0_UNITS * PIXELS_PER_SIGMA);
        szilard_count_species(split_x, &leftA, &leftB, &rightA, &rightB, &memA, &memB);
        double wall_sigma = wall_enabled ? ((double)wall_x - (double)XW1) / (double)PIXELS_PER_SIGMA : (double)L0_UNITS;
        fprintf(szlog, "%.6f,%d,%.6f,%.6f,%d,%d,%d,%d,%d,%d,%d\n",
                simulation_time, phase,
                ((double)piston_right_x - (double)XW1) / (double)PIXELS_PER_SIGMA,
                wall_sigma,
                sz_perm_dual_active ? -2 : (sz_perm_active ? sz_perm_species : -1),
                leftA, leftB, rightA, rightB, memA, memB);
    }

    fclose(szlog);
    printf("✅ Szilard experiment done. Output → %s\n", trace_path);
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
    if (cli_experiment_preset == EXPERIMENT_PRESET_SZILARD_ENGINE) {
        szilard_prepare_interactive_config();
        initialize_simulation_dimensions();
        initialize_time_scale();
        initialize_simulation_parameters();
        update_dt_runtime_for_temperature();
    }
    initialize_simulation();
    if (cli_experiment_preset == EXPERIMENT_PRESET_SZILARD_ENGINE) {
        szilard_setup_interactive_state();
        szilard_setup_timestamped_run_dir();
    }
    if (live_scene_visible) {
        live_scene_sync_from_runtime();
    }
    // If user requested pure EDMD with too many internal walls, fall back to EDMD-HYBRID.
    if (sim_mode == MODE_EDMD && num_internal_walls > EDMD_MAX_DIVIDERS) {
        printf("⚠️ %d walls exceed EDMD capacity (%d); using EDMD-HYBRID.\n",
               num_internal_walls, EDMD_MAX_DIVIDERS);
        sim_mode = MODE_EDMD_HYBRID;
    }
    // If EDMD or hybrid mode selected, initialize EDMD state from current particle arrays
    if (sim_mode == MODE_EDMD || sim_mode == MODE_EDMD_HYBRID) {
        if (g_edmd) { edmd_destroy(g_edmd); g_edmd = NULL; }
        EDMD_Params prm = {0};
        prm.boxW = (double)(XW2 - XW1);
        prm.boxH = (double)(YW2 - YW1);
        prm.radius = (double)PARTICLE_RADIUS;
        prm.N = particles_active;
        prm.cell_size = 0.0; // auto
        /* divider parameters (multi-wall capable) */
        {
            int dcount = (num_internal_walls > 0 ? num_internal_walls : 0);
            if (dcount > EDMD_MAX_DIVIDERS) dcount = EDMD_MAX_DIVIDERS;
            prm.divider_count = dcount;
            int primary_idx = primary_wall_index;
            if (primary_idx < 0 || primary_idx >= dcount) primary_idx = 0;
            int extra_idx = 0;
            for (int w = 0; w < dcount; ++w) {
                double pos_px = (all_wall_positions ? (double)all_wall_positions[w] : (double)wall_x);
                prm.divider_x[w] = pos_px - (double)XW1;
                prm.divider_thickness[w] = (double)WALL_THICKNESS;
                double mf = (double)(WALL_MASS / PARTICLE_MASS);
                if (cli_wall_masses_set && w < (int)(sizeof(cli_wall_masses) / sizeof(cli_wall_masses[0])) && cli_wall_masses[w] > 0.0f) {
                    mf = (double)cli_wall_masses[w];
                }
                prm.divider_mass[w] = mf;
                double vxw = 0.0;
                if (w == primary_idx) {
                    vxw = (double)vx_wall;
                } else if (extra_wall_velocity && extra_idx < extra_wall_count) {
                    vxw = (double)extra_wall_velocity[extra_idx];
                    extra_idx++;
                }
                prm.divider_vx[w] = vxw;
            }
        }
        /* ##CHRIS: Pass heat bath parameters to EDMD */
        prm.heatbath_enabled = heatbath_enabled;
        prm.heatbath_temperature = (double)heatbath_temperature;
        prm.thermal_wall_mode = thermal_wall_mode;  /* Pass thermal wall mode! */
        prm.mb_overshoot_factor = (double)mb_overshoot_factor;
        prm.stability_window_percent = (double)stability_window_percent;
        prm.particle_mass = (double)PARTICLE_MASS;
        prm.kB = (double)kB_effective();
        prm.pp_collisions_enabled = 1; /* ##CHRIS: default on; --szilard-no-pp overrides via szilard_sync_edmd_from_globals */
        // EDMD harmonic spring wall (interactive): model the primary divider as a true harmonic oscillator
        // when we want spring output (use --eff-output=spring).
        if (energy_measurement.spring_enabled && prm.divider_count > 0) {
            int primary_idx = primary_wall_index;
            if (primary_idx < 0 || primary_idx >= prm.divider_count) primary_idx = 0;
            const double xeq_px = (energy_measurement.equilibrium_position != 0.0f)
                                    ? (double)energy_measurement.equilibrium_position
                                    : (double)wall_x;
            prm.divider_k[primary_idx] = (double)energy_measurement.spring_constant;
            prm.divider_xeq[primary_idx] = xeq_px - (double)XW1;
        }
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
            edmd_reset_work(g_edmd);
            printf("[EDMD] initialized: N=%d  box=(%.1fσ,%.1fσ) -> (%.1fpx,%.1fpx)  R=%.4fσ -> %.3fpx  scale=%.1f px/σ\n",
                   prm.N,
                   prm.boxW / (double)PIXELS_PER_SIGMA,
                   prm.boxH / (double)PIXELS_PER_SIGMA,
                   prm.boxW, prm.boxH,
                   prm.radius / (double)PIXELS_PER_SIGMA,
                   prm.radius,
                   (double)PIXELS_PER_SIGMA);
        } else {
            edmd_reschedule_all_pp_only(g_edmd);
            edmd_reset_work(g_edmd);
            printf("[EDMD-HYBRID] initialized (PP-only): N=%d  box=(%.1fσ,%.1fσ) -> (%.1fpx,%.1fpx)  R=%.4fσ -> %.3fpx  scale=%.1f px/σ\n",
                   prm.N,
                   prm.boxW / (double)PIXELS_PER_SIGMA,
                   prm.boxH / (double)PIXELS_PER_SIGMA,
                   prm.boxW, prm.boxH,
                   prm.radius / (double)PIXELS_PER_SIGMA,
                   prm.radius,
                   (double)PIXELS_PER_SIGMA);
        }
        if (cli_experiment_preset == EXPERIMENT_PRESET_SZILARD_ENGINE) {
            /* Configure Szilard semipermeable behavior in EDMD core (dividers/piston). */
            szilard_sync_edmd_from_globals();
        }
    }

    int running = 1;
    // Track EDMD boundary work (cumulative) to feed piston_work_* deltas
    double edmd_prev_work_L = 0.0, edmd_prev_work_R = 0.0;
    double edmd_prev_work_div = 0.0;
    if (live_fft) {
        fft_cfg = kiss_fftr_alloc(FFT_SIZE, 0, NULL, NULL);
    }

    SDL_Event e;

    char energy_log_path[1024] = "energy_log.csv";
    char wall_pos_path[1024] = "wall_position.csv";
    if (g_szilard_run_active) {
        snprintf(energy_log_path, sizeof(energy_log_path), "%s/energy_log.csv", g_szilard_run_dir);
        snprintf(wall_pos_path, sizeof(wall_pos_path), "%s/wall_position.csv", g_szilard_run_dir);
    }
    FILE *logFile = fopen(energy_log_path, "w");
    FILE *wall_log = fopen(wall_pos_path, "w");
    if (!logFile || !wall_log) { printf("❌ Failed to open log files.\n"); exit(1); }

    // Optional: Szilard information/thermo log (mutual information + boundary work/heat)
    FILE *sz_log = NULL;
    double next_sz_log_abs_time = 0.0;
    char sz_log_path[1024] = {0};
    if (cli_experiment_preset == EXPERIMENT_PRESET_SZILARD_ENGINE) {
        const char *base_sz_log = (cli_szilard_log_path && cli_szilard_log_path[0]) ? cli_szilard_log_path : "szilard_info.csv";
        if (g_szilard_run_active && base_sz_log[0] != '/') {
            snprintf(sz_log_path, sizeof(sz_log_path), "%s/%s", g_szilard_run_dir, base_sz_log);
        } else {
            snprintf(sz_log_path, sizeof(sz_log_path), "%s", base_sz_log);
        }
        mkdir_p_for_file(sz_log_path);
        sz_log = fopen(sz_log_path, "w");
        if (!sz_log) { printf("❌ Failed to open szilard log: %s\n", sz_log_path); exit(1); }
        fprintf(sz_log,
                "LogIndex,BranchID,BranchLabel,Time,Phase,N,"
                "GateMove_x_sigma,GateMove_vx_sigma,GateMove_th_sigma,GateMove_mode,GateMove_species,GateMove_usesMemory,GateMove_target,"
                "GateL0_x_sigma,GateL0_vx_sigma,GateL0_th_sigma,GateL0_mode,GateL0_species,GateL0_usesMemory,GateL0_target,"
                "GatePartialObs,"
                "PistonR_x_sigma,PistonR_vx_sigma,"
                "nL_s0,nL_s1,nR_s0,nR_s1,"
                "nL_gate_s0,nL_gate_s1,nR_gate_s0,nR_gate_s1,"
                "purity_L,purity_R,purity_mean,separation_score,"
                "purity_gate_L,purity_gate_R,purity_gate_mean,separation_gate_score,"
                "mem0,mem1,memU,"
                "xy_match,xy_mismatch,xy_unknown,xy_mismatch_frac,"
                "partial_uncertainty_frac,partial_peak_uncertainty_frac,hot_speed_threshold,flip_tail_threshold,mean_speed,std_speed,"
                "memL0,memL1,memLU,memR0,memR1,memRU,"
                "memGateL0,memGateL1,memGateLU,memGateR0,memGateR1,memGateRU,"
                "I_mem_nats,I_mem_bits,"
                "I_mem_gate_nats,I_mem_gate_bits,"
                "I_nats,I_bits,"
                "I_gate_nats,I_gate_bits,"
                "kBT,F_x,F_x_gate,F_y,F_y_gate,"
                "LoadFxStart,LoadFxDrop,LoadZ,LoadE,LoadW,LoadQdiss,"
                "W_gateMove,W_gateL0,W_dividers,"
                "W_pistonL,W_pistonR,W_pistons,W_total,"
                "Q_memErase_min,Q_bath,E_ext\n");
        fflush(sz_log);
        // Baseline row at t=0
        int nL0=0,nL1=0,nR0=0,nR1=0;
        szilard_species_side_counts(&nL0,&nL1,&nR0,&nR1);
        const double purity_L = szilard_side_purity(nL0, nL1);
        const double purity_R = szilard_side_purity(nR0, nR1);
        const double purity_mean = 0.5 * (purity_L + purity_R);
        const double sep_score = szilard_separation_score(nL0, nL1, nR0, nR1);
        int mem0=0,mem1=0,memU=0;
        szilard_memory_counts(&mem0,&mem1,&memU);
        int xy_match=0,xy_mismatch=0,xy_unknown=0;
        szilard_xy_agreement_counts(&xy_match,&xy_mismatch,&xy_unknown);
        const double xy_mismatch_frac = szilard_xy_mismatch_fraction();
        const double partial_uncertainty_frac = szilard_partial_uncertainty_fraction_now();
        double mean_speed = 0.0, std_speed = 0.0;
        szilard_speed_stats(&mean_speed, &std_speed);
        const double hot_speed_threshold = sz_stage1_baseline_mode
                                         ? NAN
                                         : sz_gate_hot_speed_mult * ((mean_speed > 1e-12) ? mean_speed : fmax(std_speed, 1e-12));
        const double flip_tail_threshold = sz_stage1_baseline_mode
                                         ? NAN
                                         : mean_speed + sz_xflip_tail_sigma * ((std_speed > 1e-9) ? std_speed : fmax(mean_speed * 0.25, 1e-9));
        int memL0=0,memL1=0,memLU=0,memR0=0,memR1=0,memRU=0;
        szilard_memory_side_counts(&memL0,&memL1,&memLU,&memR0,&memR1,&memRU);
        const double I_mem_nats = szilard_mutual_information_3state_nats(memL0,memL1,memLU,memR0,memR1,memRU);
        const double I_mem_bits = I_mem_nats / log(2.0);
        const double I_nats = szilard_mutual_information_nats(nL0,nL1,nR0,nR1);
        const double I_bits = I_nats / log(2.0);
        const double kBT_ref = szilard_reference_kbt();
        const EDMD_Params* ep0 = (g_edmd ? edmd_params(g_edmd) : NULL);
        const double gate0_x = (ep0 && ep0->divider_count > 0) ? ep0->divider_x[0] : ((double)wall_x - (double)XW1);
        const double gate0_vx = (ep0 && ep0->divider_count > 0) ? ep0->divider_vx[0] : (double)vx_wall;
        const double gate0_th = (ep0 && ep0->divider_count > 0) ? ep0->divider_thickness[0] : (wall_enabled ? (double)WALL_THICKNESS : 0.0);
        const int gate0_mode = (ep0 ? ep0->divider_gate_mode[0] : 0);
        const int gate0_sp   = (ep0 ? ep0->divider_gate_species[0] : 0);
        const int gate0_use_mem = sz_perm_use_memory;
        const int gate0_tgt  = (ep0 ? ep0->divider_gate_target_side[0] : 0);

        const double gate1_x = (ep0 && ep0->divider_count > 1) ? ep0->divider_x[1] : ((double)sz_interactive_target_px - (double)XW1);
        const double gate1_vx = (ep0 && ep0->divider_count > 1) ? ep0->divider_vx[1] : 0.0;
        const double gate1_th = (ep0 && ep0->divider_count > 1) ? ep0->divider_thickness[1] : 0.0;
        const int gate1_mode = (ep0 ? ep0->divider_gate_mode[1] : 0);
        const int gate1_sp   = (ep0 ? ep0->divider_gate_species[1] : 0);
        const int gate1_use_mem = sz_fixed_gate_use_memory;
        const int gate1_tgt  = (ep0 ? ep0->divider_gate_target_side[1] : 0);

        int nLg0=0,nLg1=0,nRg0=0,nRg1=0;
        szilard_species_side_counts_at_split((double)XW1 + gate0_x, &nLg0,&nLg1,&nRg0,&nRg1);
        const double purity_gate_L = szilard_side_purity(nLg0, nLg1);
        const double purity_gate_R = szilard_side_purity(nRg0, nRg1);
        const double purity_gate_mean = 0.5 * (purity_gate_L + purity_gate_R);
        const double sep_gate_score = szilard_separation_score(nLg0, nLg1, nRg0, nRg1);

        int memGateL0=0,memGateL1=0,memGateLU=0,memGateR0=0,memGateR1=0,memGateRU=0;
        szilard_memory_side_counts_at_split((float)((double)XW1 + gate0_x),
                                            &memGateL0,&memGateL1,&memGateLU,&memGateR0,&memGateR1,&memGateRU);
        const double I_mem_gate_nats = szilard_mutual_information_3state_nats(
            memGateL0,memGateL1,memGateLU,memGateR0,memGateR1,memGateRU
        );
        const double I_mem_gate_bits = I_mem_gate_nats / log(2.0);
        const double I_gate_nats = szilard_mutual_information_nats(nLg0,nLg1,nRg0,nRg1);
        const double I_gate_bits = I_gate_nats / log(2.0);
        const double F_x = szilard_available_free_energy_from_nats(I_nats, nL0 + nL1 + nR0 + nR1);
        const double F_x_gate = szilard_available_free_energy_from_nats(I_gate_nats, nLg0 + nLg1 + nRg0 + nRg1);
        const double F_y = szilard_available_free_energy_from_nats(I_mem_nats, memL0 + memL1 + memLU + memR0 + memR1 + memRU);
        const double F_y_gate = szilard_available_free_energy_from_nats(
            I_mem_gate_nats,
            memGateL0 + memGateL1 + memGateLU + memGateR0 + memGateR1 + memGateRU
        );

        const double pistR_x = (ep0 ? ep0->pistonR_x : ((double)piston_right_x - (double)XW1));
        const double pistR_vx = (ep0 ? ep0->pistonR_vx : (double)vx_piston_right);

        const double wG0 = (g_edmd ? edmd_work_divider_i(g_edmd, 0) : 0.0);
        const double wG1 = (g_edmd ? edmd_work_divider_i(g_edmd, 1) : 0.0);
        const double wD  = wG0 + wG1;
        const double wPL = (g_edmd ? edmd_work_pistonL(g_edmd) : 0.0);
        const double wPR = (g_edmd ? edmd_work_pistonR(g_edmd) : 0.0);
        const double wP  = wPL + wPR;
        const double wT  = wD + wP;
        const double qB  = (g_edmd ? edmd_heat_bath(g_edmd) : 0.0);
        const double eext = wT + qB;

        fprintf(sz_log,
                "%lld,%d,%s,%.6f,%d,%d,"
                "%.6f,%.6f,%.6f,%d,%d,%d,%d,"
                "%.6f,%.6f,%.6f,%d,%d,%d,%d,"
                "%d,"
                "%.6f,%.6f,"
                "%d,%d,%d,%d,"
                "%d,%d,%d,%d,"
                "%.12g,%.12g,%.12g,%.12g,"
                "%.12g,%.12g,%.12g,%.12g,"
                "%d,%d,%d,"
                "%d,%d,%d,%.12g,"
                "%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,"
                "%d,%d,%d,%d,%d,%d,"
                "%d,%d,%d,%d,%d,%d,"
                "%.12g,%.12g,"
                "%.12g,%.12g,"
                "%.12g,%.12g,"
                "%.12g,%.12g,"
                "%.12g,%.12g,%.12g,%.12g,%.12g,"
                "%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,"
                "%.12g,%.12g,%.12g,"
                "%.12g,%.12g,%.12g,%.12g,"
                "%.12g,%.12g,%.12g\n",
                ++sz_log_index, sz_log_branch_id, sz_log_branch_label, simulation_time, sz_interactive_phase, particles_active,
                gate0_x / (double)PIXELS_PER_SIGMA, gate0_vx / (double)PIXELS_PER_SIGMA, gate0_th / (double)PIXELS_PER_SIGMA, gate0_mode, gate0_sp, gate0_use_mem, gate0_tgt,
                gate1_x / (double)PIXELS_PER_SIGMA, gate1_vx / (double)PIXELS_PER_SIGMA, gate1_th / (double)PIXELS_PER_SIGMA, gate1_mode, gate1_sp, gate1_use_mem, gate1_tgt,
                sz_partial_obs_enabled,
                pistR_x / (double)PIXELS_PER_SIGMA, pistR_vx / (double)PIXELS_PER_SIGMA,
                nL0,nL1,nR0,nR1,
                nLg0,nLg1,nRg0,nRg1,
                purity_L,purity_R,purity_mean,sep_score,
                purity_gate_L,purity_gate_R,purity_gate_mean,sep_gate_score,
                mem0,mem1,memU,
                xy_match,xy_mismatch,xy_unknown,xy_mismatch_frac,
                partial_uncertainty_frac,sz_partial_obs_peak_uncertainty_frac,hot_speed_threshold,flip_tail_threshold,mean_speed,std_speed,
                memL0,memL1,memLU,memR0,memR1,memRU,
                memGateL0,memGateL1,memGateLU,memGateR0,memGateR1,memGateRU,
                I_mem_nats,I_mem_bits,
                I_mem_gate_nats,I_mem_gate_bits,
                I_nats,I_bits,
                I_gate_nats,I_gate_bits,
                kBT_ref,F_x,F_x_gate,F_y,F_y_gate,
                sz_load_fx_start,sz_load_fx_drop,sz_load_z,sz_load_energy,sz_load_work,sz_load_qdiss,
                wG0,wG1,wD,
                wPL,wPR,wP,wT,
                sz_memory_erase_heat_min,qB,eext);
        fflush(sz_log);
        next_sz_log_abs_time = (cli_szilard_log_dt > 0.0) ? (simulation_time + cli_szilard_log_dt) : 0.0;
    }

    // Optional: detailed kinematics log for energy_transfer investigations.
    FILE *et_wall_log = NULL;
    if (cli_experiment_preset == EXPERIMENT_PRESET_ENERGY_TRANSFER) {
        et_wall_log = fopen("energy_transfer_wall_positions.csv", "w");
        if (!et_wall_log) { printf("❌ Failed to open energy_transfer_wall_positions.csv\n"); exit(1); }
        fprintf(et_wall_log,
                "TimeAbs,TimeAfterRelease,Released,"
                "W0_x_sigma,W0_v,W1_x_sigma,W1_v,W2_x_sigma,W2_v,W3_x_sigma,W3_v,W4_x_sigma,W4_v,"
                "PistonR_x_sigma,PistonR_v,PistonL_x_sigma,PistonL_v,SegCounts,"
                "Substeps,KEp,KEw,SpringE_now,SpringE_max,PistonWorkDelta,PistonWorkMax,MechE,"
                "WallKEs,SegKEs,SegTs\n");
    }

    // Optional: distribution log (x/y/speed) for energy_transfer
    if (cli_dist_log && cli_experiment_preset == EXPERIMENT_PRESET_ENERGY_TRANSFER) {
        dist_bins_runtime = (cli_dist_bins > 0) ? cli_dist_bins : 50;
        if (dist_bins_runtime > 2000) dist_bins_runtime = 2000;
        dist_hist_x = (int*)calloc((size_t)dist_bins_runtime, sizeof(int));
        dist_hist_y = (int*)calloc((size_t)dist_bins_runtime, sizeof(int));
        dist_hist_v = (int*)calloc((size_t)dist_bins_runtime, sizeof(int));
        if (!dist_hist_x || !dist_hist_y || !dist_hist_v) {
            printf("❌ Failed to allocate distribution histograms\n");
            exit(1);
        }
        dist_box_w_sigma = ((double)XW2 - (double)XW1) / (double)PIXELS_PER_SIGMA;
        dist_box_h_sigma = ((double)YW2 - (double)YW1) / (double)PIXELS_PER_SIGMA;
        if (cli_dist_vmax > 0.0f) {
            dist_vmax_runtime = (double)cli_dist_vmax;
        } else {
            const double kbt = (double)kB_effective() * (double)temperature_runtime;
            const double v_rms = (kbt > 0.0) ? sqrt(2.0 * kbt / (double)PARTICLE_MASS) : 1.0;
            dist_vmax_runtime = 6.0 * v_rms;
            if (!(dist_vmax_runtime > 0.0)) dist_vmax_runtime = 1.0;
        }
        char dist_path[1024];
        build_dist_log_path(dist_path, sizeof(dist_path));
        dist_log = fopen(dist_path, "w");
        if (!dist_log) { printf("❌ Failed to open %s\n", dist_path); exit(1); }
        fprintf(dist_log, "TimeAbs,TimeAfterRelease,Released,Bins,XMin,XMax,YMin,YMax,VMax,HistX,HistY,HistV\n");
        dist_next_log_abs_time = 0.0;
    }

    // ##CHRIS: Initialize substep logging pointer
    g_substep_log = wall_log;
    g_substep_base_time = 0.0;

    // Write run parameters for Python analysis tooling
    if (g_szilard_run_active) {
        char rp[1024];
        char wrp[1024];
        snprintf(rp, sizeof(rp), "%s/run_params.json", g_szilard_run_dir);
        snprintf(wrp, sizeof(wrp), "%s/wall_position_run_params.json", g_szilard_run_dir);
        write_run_params_json_to(rp);
        write_run_params_json_to(wrp);
    } else {
        write_run_params_json();
        // Also write a wall_position-specific params file, so GUI + wall_position.csv analysis
        // can't accidentally pick up a newer run_params.json from a different batch run.
        write_run_params_json_to("wall_position_run_params.json");
    }

    if (log_packing_fraction) {
        log_packing_fractions(wall_log, "main");
    } else {
        fprintf(wall_log, "Time, Wall_X, Displacement(σ), Left_Count, Right_Count\n");
    }
    fprintf(logFile, "KE_particles,KE_wall,KE_total,T,S_v,F\n");

    // initial sim state
    steps_elapsed        = 0;
    if (cli_experiment_preset == EXPERIMENT_PRESET_SZILARD_ENGINE) {
        // Szilard starts with divider disabled and piston closing the right half
        wall_hold_enabled = true;
        wall_is_released  = false;
        wall_enabled      = 0;
        wall_hold_steps   = INT_MAX / 2;
        sz_perm_active    = 0;
        sz_perm_dual_active = 0;
    } else {
        wall_hold_enabled = !cli_auto_release_wall;  // Auto-release if flag set
        wall_is_released  = cli_auto_release_wall;   // Auto-release if flag set
        wall_enabled      = 1;               // ensure wall is drawn/active
    }
    wall_x_old           = wall_x;          // first CCD step sees correct "old" pose
    missed_collision_events = 0;
    worst_penetration_observed = 0.0;

    // ##CHRIS: Auto-start simulation in headless mode or with auto-release
    if (cli_headless || cli_auto_release_wall) {
        simulation_started = 1;
        printf("▶ Simulation auto-started (headless=%d, auto-release=%d)\n", cli_headless, cli_auto_release_wall);
    }

    // ##CHRIS: wall_release_time is set in the loop when wall is actually released
    // Don't pre-set it here, or the logging initialization code won't trigger

    // ##CHRIS: In headless mode, we don't use SDL timing or events
    Uint32 last_time = cli_headless ? 0 : SDL_GetTicks();
    float accumulator = 0.0f;
    // Logging cadence (relative to release time)
    double next_wall_log_rel_time = 0.0;
    // Logging cadence (absolute time) for energy_transfer_wall_positions.csv
    double next_et_log_abs_time = 0.0;

    while (running) {
        // ##CHRIS: Handle timing differently in headless vs interactive mode
        if (!cli_headless) {
            Uint32 current_time = SDL_GetTicks();
            float frame_time = (current_time - last_time) / 1000.0f;
            last_time = current_time;
            if (frame_time > 0.1f) frame_time = 0.1f;
            accumulator += frame_time * time_scale_runtime;

            // events
            while (SDL_PollEvent(&e)) {
                if (e.type == SDL_QUIT) running = 0;
                if (e.type == SDL_WINDOWEVENT &&
                    (e.window.event == SDL_WINDOWEVENT_SIZE_CHANGED || e.window.event == SDL_WINDOWEVENT_RESIZED)) {
                    // Pure UI resize; physics box remains SIM_WIDTH×SIM_HEIGHT.
                    // Use drawable pixels, not logical window points, so Retina/HiDPI
                    // layouts and mouse hit boxes use the same coordinate system.
                    gui_sync_renderer_canvas_size();
                }
                if (live_controls_handle_event(&e)) {
                    continue;
                }
                keysSimulation(e);
                keysPiston(e);
            }
        } else {
            // ##CHRIS: In headless mode, just accumulate one timestep per iteration (run as fast as possible)
            accumulator += fixed_dt_runtime;
        }

        if (live_scene_rebuild_requested) {
            const bool spring_was_enabled = energy_measurement.spring_enabled;

            live_scene_rebuild_requested = 0;
            reset_simulation_state();

            wall_enabled = (num_internal_walls > 0) ? 1 : 0;
            wall_hold_enabled = !cli_auto_release_wall;
            wall_is_released = cli_auto_release_wall;
            simulation_started = 1;

            if (cli_enable_energy_measurement || cli_experiment_preset == EXPERIMENT_PRESET_ENERGY_TRANSFER) {
                initialize_energy_measurement();
                energy_measurement.spring_enabled = spring_was_enabled;
                if (num_internal_walls > 0 && all_wall_positions) {
                    energy_measurement.equilibrium_position = all_wall_positions[leftmost_wall_index];
                    energy_measurement.spring_force = 0.0f;
                    energy_measurement.spring_energy = 0.0f;
                    energy_measurement.spring_energy_max = 0.0f;
                }
                if (!spring_was_enabled) {
                    energy_measurement.spring_force = 0.0f;
                    energy_measurement.spring_energy = 0.0f;
                }
            }

            live_rebuild_edmd_from_current_state();
            if (cli_experiment_preset == EXPERIMENT_PRESET_SZILARD_ENGINE) {
                szilard_sync_edmd_from_globals();
            }

            edmd_prev_work_L = 0.0;
            edmd_prev_work_R = 0.0;
            edmd_prev_work_div = 0.0;
            gui_mb_speed_vmax_ref = 0.0f;
            accumulator = 0.0f;
            if (!cli_quiet) {
                printf("Live scene rebuilt: walls=%d eta=%.3f r=%.4fσ\n",
                       num_internal_walls,
                       (double)cli_packing_fraction,
                       (double)PARTICLE_RADIUS / (double)PIXELS_PER_SIGMA);
            }
        }

        int left_particles = 0, right_particles = 0;

        while (accumulator >= fixed_dt_runtime) {
            if (simulation_started && !paused) {
                if (sz_branch_restore_requested) {
                    sz_branch_restore_requested = 0;
                    szilard_restore_branch_snapshot();
                    edmd_prev_work_L = 0.0;
                    edmd_prev_work_R = 0.0;
                    edmd_prev_work_div = 0.0;
                    accumulator = 0.0f; // restore visually without immediately advancing a step
                    if (!cli_quiet) {
                        printf("↩️  Szilard: restored saved post-G branch. Work/heat counters reset for fresh H/J comparison.\n");
                    }
                    break;
                }
                // --- START OF STEP ---
                wall_x_old = wall_x;
                steps_elapsed++;

                // GUI convenience: auto-release after the configured hold steps.
                // TIME mode releases inside update_wall(), but EDMD mode does not call update_wall(),
                // so we do it here in a solver-agnostic place.
                if (cli_auto_release_after_hold && wall_hold_enabled && !wall_is_released && steps_elapsed >= wall_hold_steps) {
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
                    if (!cli_quiet) {
                        printf("🚀 Auto-release-after-hold: released at step %d (t=%.3f)\n", steps_elapsed, simulation_time);
                    }
                }

                // reset accumulator at the *start* of the step
                wall_impulse_x_accum = 0.f;

                if (sim_mode == MODE_EDMD) {
                    const bool sz_mode = (cli_experiment_preset == EXPERIMENT_PRESET_SZILARD_ENGINE && sz_interactive_enabled);
                    if (sz_mode) {
                        /* In Szilard EDMD mode, do not overwrite divider/piston configuration each step.
                           The protocol (G/H) drives boundary velocities/positions, and semipermeable
                           behavior is implemented inside the EDMD core. */
                        double t0 = edmd_time(g_edmd);
                        edmd_advance_to(g_edmd, t0 + (double)fixed_dt_runtime);

                        const EDMD_Particle* P = edmd_particles(g_edmd);
                        for (int i = 0; i < particles_active; i++) {
                            X[i] = (double)XW1 + (double)P[i].x;
                            Y[i] = (double)YW1 + (double)P[i].y;
                            Vx[i] = (double)P[i].vx;
                            Vy[i] = (double)P[i].vy;
                        }

                        const EDMD_Params* ep = edmd_params(g_edmd);
                        piston_left_x  = (float)(XW1 + ep->pistonL_x - 5.0);
                        piston_right_x = (float)(XW1 + ep->pistonR_x);
                        if (ep->divider_count > 0) {
                            wall_x = (double)XW1 + ep->divider_x[0];
                            vx_wall = (float)ep->divider_vx[0];
                        }

                        /* Track piston work for HUD consistency (even though Szilard doesn't use it for efficiency). */
                        double wL = edmd_work_pistonL(g_edmd);
                        double wR = edmd_work_pistonR(g_edmd);
                        double wD = edmd_work_divider(g_edmd);
                        piston_work_left  += (wL - edmd_prev_work_L);
                        piston_work_right += (wR - edmd_prev_work_R);
                        edmd_prev_work_L = wL;
                        edmd_prev_work_R = wR;
                        edmd_prev_work_div = wD;

                        recompute_segment_stats_counts_and_temperature();
                        if (segment_counts && segment_count >= 2) {
                            left_particles = segment_counts[0];
                            right_particles = segment_counts[1];
                        }
                    } else {
                    // Current divider state from EDMD (multi-wall)
                    const EDMD_Params* ep_state = edmd_params(g_edmd);
                    int dcount = ep_state->divider_count;
                    if (!wall_enabled) dcount = 0;
                    if (dcount > EDMD_MAX_DIVIDERS) dcount = EDMD_MAX_DIVIDERS;
                    int primary_idx = primary_wall_index;
                    if (primary_idx < 0 || primary_idx >= dcount) primary_idx = 0;

                    double div_x[EDMD_MAX_DIVIDERS];
                    double div_vx[EDMD_MAX_DIVIDERS];
                    double div_mass[EDMD_MAX_DIVIDERS];
                    double div_th[EDMD_MAX_DIVIDERS];
                    int extra_idx = 0;

                    for (int w = 0; w < dcount; ++w) {
                        div_x[w] = ep_state->divider_x[w];
                        div_vx[w] = ep_state->divider_vx[w];
                        div_th[w] = (double)WALL_THICKNESS;
                        double mf = (double)(WALL_MASS / PARTICLE_MASS);
                        if (cli_wall_masses_set && w < (int)(sizeof(cli_wall_masses) / sizeof(cli_wall_masses[0])) && cli_wall_masses[w] > 0.0f) {
                            mf = (double)cli_wall_masses[w];
                        }
                        // Hold walls (infinite mass) before release (match time-mode semantics).
                        bool hold_active = false;
                        if (w == primary_idx) {
                            hold_active = wall_hold_enabled && !wall_is_released;
                        } else if (extra_wall_hold_enabled && extra_wall_is_released) {
                            int j = (w < primary_idx) ? w : (w - 1);
                            if (j >= 0 && j < extra_wall_count) {
                                hold_active = extra_wall_hold_enabled[j] && !extra_wall_is_released[j];
                            }
                        }
                        if (hold_active) {
                            div_mass[w] = 0.0;
                            div_vx[w] = 0.0;
                        } else {
                            div_mass[w] = mf;
                        }

                        if (w == primary_idx) {
                            wall_x = (double)XW1 + div_x[w];
                            vx_wall = div_vx[w];
                        } else if (extra_wall_positions && extra_wall_velocity && extra_idx < extra_wall_count) {
                            extra_wall_positions[extra_idx] = (float)((double)XW1 + div_x[w]);
                            extra_wall_velocity[extra_idx] = (float)div_vx[w];
                            extra_idx++;
                        }
                    }

                    piston_left_x  = (float)(XW1 + ep_state->pistonL_x - 5.0);
                    piston_right_x = (float)(XW1 + ep_state->pistonR_x);

                    // In pure EDMD mode, the spring wall is handled by the EDMD core (harmonic divider)
                    // when --eff-output=spring is used. Do not apply timestep impulses here.

                    // Configure divider motion for EDMD (positions from EDMD state; velocities possibly updated)
                    edmd_config_dividers(g_edmd, dcount, div_x, div_th);
                    edmd_set_divider_motions(g_edmd, dcount, div_mass, div_vx);

                    // Keep pistons in sync with interactive controls; if velocities changed, update
                    static double prev_vxL = 0.0, prev_vxR = 0.0;
                    // Apply protocol-driven velocities (do not move positions here; EDMD owns positions)
                    apply_piston_protocol(fixed_dt_runtime);
                    if (vx_piston_left != prev_vxL || vx_piston_right != prev_vxR) {
                        const EDMD_Params* epc = edmd_params(g_edmd);
                        edmd_config_pistons(g_edmd,
                            1, epc->pistonL_x, (double)vx_piston_left, epc->pistonL_mass,
                            1, epc->pistonR_x, (double)vx_piston_right, epc->pistonR_mass);
                        prev_vxL = vx_piston_left; prev_vxR = vx_piston_right;
                    }
                    // Rebuild events after divider/piston velocity changes
                    edmd_reschedule_all(g_edmd);

                    double t0 = edmd_time(g_edmd);
                    edmd_advance_to(g_edmd, t0 + (double)fixed_dt_runtime);
                    // After advance, fix any divider overlaps caused by moving boundary and rebuild schedule.
                    // Szilard semipermeable walls are handled inside EDMD; do not apply a hard-wall overlap reflect.
                    const bool sz_gate = (cli_experiment_preset == EXPERIMENT_PRESET_SZILARD_ENGINE &&
                                          (sz_perm_active || sz_perm_dual_active || sz_fixed_gate_active));
                    if (!sz_gate) edmd_divider_resolve_overlaps(g_edmd);
                    // Szilard: if any PP overlaps exist (rare, but possible after non-EDMD push-outs),
                    // relax them to keep EDMD in a valid state.
                    edmd_relax_pp_overlaps_szilard(g_edmd);
                    edmd_enforce_szilard_gate_sides(g_edmd);
                    if (szilard_repair_edmd_particles_for_pistons_and_pockets()) {
                        edmd_relax_pp_overlaps_szilard(g_edmd);
                        edmd_reschedule_all(g_edmd);
                    }
                    // copy back for rendering/diagnostics
                    const EDMD_Particle* P = edmd_particles(g_edmd);
                    for (int i=0;i<particles_active;i++) {
                        X[i] = (double)XW1 + (double)P[i].x;
                        Y[i] = (double)YW1 + (double)P[i].y;
                        Vx[i] = (double)P[i].vx;
                        Vy[i] = (double)P[i].vy;
                    }
                    const EDMD_Params* ep = edmd_params(g_edmd);
                    int dcount2 = ep->divider_count;
                    if (!wall_enabled) dcount2 = 0;
                    if (dcount2 > EDMD_MAX_DIVIDERS) dcount2 = EDMD_MAX_DIVIDERS;
                    int primary_idx2 = primary_wall_index;
                    if (primary_idx2 < 0 || primary_idx2 >= dcount2) primary_idx2 = 0;
                    int extra_idx2 = 0;
                    for (int w = 0; w < dcount2; ++w) {
                        if (w == primary_idx2) {
                            wall_x = (double)XW1 + ep->divider_x[w];
                            vx_wall = (float)ep->divider_vx[w];
                        } else if (extra_wall_positions && extra_wall_velocity && extra_idx2 < extra_wall_count) {
                            extra_wall_positions[extra_idx2] = (float)((double)XW1 + ep->divider_x[w]);
                            extra_wall_velocity[extra_idx2] = (float)ep->divider_vx[w];
                            extra_idx2++;
                        }
                    }
                    // Sync piston face positions for visualization (assume left piston width=5)
                    piston_left_x  = (float)(XW1 + ep->pistonL_x - 5.0);
                    piston_right_x = (float)(XW1 + ep->pistonR_x);
                    // Accumulate work done by pistons (EDMD boundary KE change)
                    double wL = edmd_work_pistonL(g_edmd);
                    double wR = edmd_work_pistonR(g_edmd);
                    double wD = edmd_work_divider(g_edmd);
                    piston_work_left  += (wL - edmd_prev_work_L);
                    piston_work_right += (wR - edmd_prev_work_R);
                    edmd_prev_work_L = wL;
                    edmd_prev_work_R = wR;
                    edmd_prev_work_div = wD;
                    /* rebuild per-segment stats */
                    recompute_segment_stats_counts_and_temperature();
                    // ##CHRIS: Update local particle counts from segment stats for CSV logging
                    if (segment_counts && segment_count >= 2) {
                        left_particles = segment_counts[0];
                        right_particles = segment_counts[1];
                    }
                    }
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
                    // ##CHRIS: Update local particle counts from segment stats for CSV logging
                    if (segment_counts && segment_count >= 2) {
                        left_particles = segment_counts[0];
                        right_particles = segment_counts[1];
                    }
                } else {
                    // ##CHRIS: Particle-life mode uses custom physics (no hard-disk collisions)
                    if (pl_enabled) {
                        extern double cli_pl_dt;  // CLI override for particle-life timestep
                        double pl_dt = (cli_pl_dt > 0) ? cli_pl_dt : (double)fixed_dt_runtime;
                        pl_update_physics(pl_dt);
                } else {
                    // time-driven (Euler ballistic with CCD)
                    // IMPORTANT: update piston *velocity* (protocol) before particle CCD, so CCD uses the correct moving piston.
                    apply_piston_protocol(fixed_dt_runtime);
                    g_substep_base_time = simulation_time; // ##CHRIS: Update base time for substep logging
                    update_particles_with_substepping(fixed_dt_runtime, &left_particles, &right_particles,
                                                          simulation_time, wall_is_released, wall_release_time);
                }
                }

                //============================================================================
                // ##CHRIS: Apply Andersen thermostat after particle updates
                // This provides bulk thermalization to prevent "cold core" problem
                //============================================================================
                if (andersen_enabled && heatbath_enabled) {
                    andersen_thermostat_step((double)fixed_dt_runtime);

                    // For EDMD modes, need to write back to EDMD state
                    if (sim_mode == MODE_EDMD || sim_mode == MODE_EDMD_HYBRID) {
                        EDMD_Particle* Pw = (EDMD_Particle*)edmd_particles(g_edmd);
                        for (int i = 0; i < particles_active; i++) {
                            Pw[i].vx = (double)Vx[i];
                            Pw[i].vy = (double)Vy[i];
                        }
                    }
                }
                //============================================================================

                // wall velocity is already updated inside ccd_wall_step (impulses applied immediately).
                // wall_impulse_x_accum is kept for diagnostics only.

                // Advance piston positions after collision handling for this step.
                if (sim_mode != MODE_EDMD) {
                    move_left_piston(fixed_dt_runtime);
                    move_right_piston(fixed_dt_runtime);
                }
                if (cli_experiment_preset == EXPERIMENT_PRESET_SZILARD_ENGINE && sz_interactive_enabled) {
                    const double dt_sigma = (double)fixed_dt_runtime / (double)PIXELS_PER_SIGMA;
                    szilard_apply_species_flips(dt_sigma);
                    if (sz_partial_obs_enabled || sz_perm_use_memory || sz_fixed_gate_use_memory) {
                        szilard_refresh_gate_labels();
                    }
                    if (g_edmd && sim_mode == MODE_EDMD) {
                        /* Keep EDMD gate geometry/modes live during Szilard motion.
                           This is required for the fixed L0 memory gate: at branch start the
                           moving divider overlaps L0 so the fixed gate is suppressed once,
                           then it must become physical as soon as the moving divider pulls away. */
                        szilard_sync_edmd_from_globals();
                    }
                    szilard_update_partial_uncertainty_peak();
                }
                if (cli_experiment_preset == EXPERIMENT_PRESET_SZILARD_ENGINE) {
                    szilard_update_interactive_poststep();
                }
                if (sim_mode != MODE_EDMD) {
                    update_wall(fixed_dt_runtime, L0_UNITS);
                    update_extra_walls(fixed_dt_runtime, L0_UNITS);
                }
                sync_all_wall_positions();
                update_wall_metrics();
                update_energy_measurement(fixed_dt_runtime);

                // ##CHRIS: EDMD uses pixel coordinates, so dt is in pixel-time units.
                // Convert to σ-time units for logging by dividing by PIXELS_PER_SIGMA.
                simulation_time += fixed_dt_runtime / PIXELS_PER_SIGMA;

                // Update piston push energy tracking (interactive mode)
                if (piston_push_tracking) {
                    double Wp = piston_work_left + piston_work_right;
                    double delta = Wp - piston_work_baseline;
                    if (delta > piston_work_delta_max) {
                        piston_work_delta_max = delta;
                    }
                }
                maybe_stop_right_piston_on_work_target();

                // Diagnostics (add piston work, hits)
                if (!cli_quiet && (steps_elapsed % 1000) == 0) {
                    double KEp = kinetic_energy();
                    double KEw = kinetic_energy_walls_total();
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

                if (!cli_quiet && (steps_elapsed % 500) == 0) {
                    double ke_p = kinetic_energy();
                    double ke_w = kinetic_energy_walls_total();
                    printf("[ENERGY] step=%d  KE_p=%.6e  KE_w=%.6e  KE_tot=%.6e  vx_wall=%.6g\n",
                        steps_elapsed, ke_p, ke_w, ke_p + ke_w, vx_wall);
                }

                if (!cli_quiet && (steps_elapsed % 1000) == 0) {
                    printf("KE=%.9f  vx_wall=%.6g  J_accum=%.6g\n",
                        kinetic_energy(), vx_wall, wall_impulse_x_accum);
                }

                // establish release time once
                if (wall_is_released && wall_release_time < 0.0) {
                    wall_release_time = simulation_time;
                    if (!cli_quiet) {
                        printf("✅ Wall released at t = %.3f\n", wall_release_time);
                    }
                    next_wall_log_rel_time = 0.0; // reset relative timeline for output cadence
                }

                // Log all wall + piston positions for diagnosing energy-transfer transients.
                if (et_wall_log) {
                    bool et_do_log = false;
                    if (cli_output_dt <= 1.0f) {
                        et_do_log = true;  // per fixed step
                    } else {
                        et_do_log = (simulation_time + 1e-12 >= next_et_log_abs_time);
                    }

                    if (et_do_log) {
                        const double t_abs = simulation_time;
                        const int released = wall_is_released ? 1 : 0;
                        const double t_rel = (released && wall_release_time >= 0.0) ? (simulation_time - wall_release_time) : NAN;

                        double wx_sigma[5] = {NAN, NAN, NAN, NAN, NAN};
                        double wv[5] = {NAN, NAN, NAN, NAN, NAN};

                        // Primary wall (spring wall)
                        wx_sigma[0] = (double)wall_x / (double)PIXELS_PER_SIGMA;
                        wv[0] = (double)vx_wall;

                        // Extra internal walls (if present)
                        for (int j = 0; j < extra_wall_count && j < 4; ++j) {
                            wx_sigma[j + 1] = (double)extra_wall_positions[j] / (double)PIXELS_PER_SIGMA;
                            wv[j + 1] = (double)extra_wall_velocity[j];
                        }

                        const double pistonR_x_sigma = (double)piston_right_x / (double)PIXELS_PER_SIGMA;
                        const double pistonL_x_sigma = (double)piston_left_x / (double)PIXELS_PER_SIGMA;
                        const double pistonR_v = (double)vx_piston_right;
                        const double pistonL_v = (double)vx_piston_left;

                        char seg_buf[256] = {0};
                        if (segment_counts && segment_count > 0) {
                            size_t off = 0;
                            for (int s = 0; s < segment_count && off + 16 < sizeof(seg_buf); ++s) {
                                int n = snprintf(seg_buf + off, sizeof(seg_buf) - off, "%s%d", s ? ";" : "", segment_counts[s]);
                                if (n < 0) break;
                                off += (size_t)n;
                            }
                        }

                        const int sub_used = g_last_substeps_used;

                        // Energetics (use kB_effective so temperature comparisons are meaningful under --kbt1)
                        const double KEp_now = kinetic_energy();
                        const double springE_now = (energy_measurement.spring_enabled ? (double)energy_measurement.spring_energy : 0.0);
                        const double springE_max = (energy_measurement.spring_enabled ? (double)energy_measurement.spring_energy_max : 0.0);

                        // Piston work delta since push started (same quantity shown on HUD)
                        double pistonW_delta = 0.0;
                        double pistonW_max = 0.0;
                        if (piston_push_tracking) {
                            double Wp_now = piston_work_left + piston_work_right;
                            pistonW_delta = Wp_now - piston_work_baseline;
                            pistonW_max = piston_work_delta_max;
                        }

                        // Wall kinetic energies (mass factors may differ per wall if provided via --wall-mass-factors)
                        double wallKE_sum = 0.0;
                        char wallke_buf[256] = {0};
                        {
                            size_t off = 0;
                            const double mf0 = (cli_wall_masses_set && num_internal_walls > 0 && cli_wall_masses[0] > 0.0f)
                                                 ? (double)cli_wall_masses[0]
                                                 : (double)(WALL_MASS / PARTICLE_MASS);
                            for (int w = 0; w < 5; ++w) {
                                double ke_w = NAN;
                                if (w == 0) {
                                    // primary internal wall
                                    if (num_internal_walls > 0) {
                                        const double mf = mf0;
                                        const double M = mf * (double)PARTICLE_MASS;
                                        ke_w = 0.5 * M * (double)vx_wall * (double)vx_wall;
                                    }
                                } else {
                                    // extra internal walls
                                    const int j = w - 1;
                                    if (j < extra_wall_count && extra_wall_velocity) {
                                        double mf = mf0;
                                        if (cli_wall_masses_set && w < (int)(sizeof(cli_wall_masses) / sizeof(cli_wall_masses[0])) && cli_wall_masses[w] > 0.0f) {
                                            mf = (double)cli_wall_masses[w];
                                        }
                                        const double M = mf * (double)PARTICLE_MASS;
                                        const double vw = (double)extra_wall_velocity[j];
                                        ke_w = 0.5 * M * vw * vw;
                                    }
                                }

                                if (isfinite(ke_w)) wallKE_sum += ke_w;

                                if (off + 32 >= sizeof(wallke_buf)) break;
                                int n = snprintf(wallke_buf + off, sizeof(wallke_buf) - off, "%s%.6e",
                                                 w ? ";" : "", isfinite(ke_w) ? ke_w : NAN);
                                if (n < 0) break;
                                off += (size_t)n;
                            }
                        }

                        // Per-segment KE + temperature for diagnosing transient compression/cooling.
                        char segke_buf[512] = {0};
                        char segt_buf[512] = {0};
                        if (segment_counts && segment_count > 0) {
                            size_t off_ke = 0, off_t = 0;
                            for (int s = 0; s < segment_count && off_ke + 32 < sizeof(segke_buf) && off_t + 32 < sizeof(segt_buf); ++s) {
                                double ke_val = segment_ke ? segment_ke[s] : 0.0;
                                double t_val = segment_temperature ? segment_temperature[s] : 0.0;
                                int nk = snprintf(segke_buf + off_ke, sizeof(segke_buf) - off_ke, "%s%.6e", s ? ";" : "", ke_val);
                                int nt = snprintf(segt_buf + off_t, sizeof(segt_buf) - off_t, "%s%.6f", s ? ";" : "", t_val);
                                if (nk < 0 || nt < 0) break;
                                off_ke += (size_t)nk;
                                off_t += (size_t)nt;
                            }
                        }

                        const double mechE_now = KEp_now + wallKE_sum + springE_now;

                        fprintf(et_wall_log,
                                "%.6f, %.6f, %d,"
                                "%.6f, %.6e, %.6f, %.6e, %.6f, %.6e, %.6f, %.6e, %.6f, %.6e,"
                                "%.6f, %.6e, %.6f, %.6e, \"%s\","
                                "%d, %.6e, %.6e, %.6e, %.6e, %.6e, %.6e, %.6e, \"%s\", \"%s\", \"%s\"\n",
                                t_abs, t_rel, released,
                                wx_sigma[0], wv[0],
                                wx_sigma[1], wv[1],
                                wx_sigma[2], wv[2],
                                wx_sigma[3], wv[3],
                                wx_sigma[4], wv[4],
                                pistonR_x_sigma, pistonR_v,
                                pistonL_x_sigma, pistonL_v,
                                seg_buf,
                                sub_used,
                                KEp_now,
                                wallKE_sum,
                                springE_now,
                                springE_max,
                                pistonW_delta,
                                pistonW_max,
                                mechE_now,
                                wallke_buf,
                                segke_buf,
                                segt_buf);
                        fflush(et_wall_log);

                        if (cli_output_dt > 1.0f) {
                            next_et_log_abs_time += (double)cli_output_dt;
                        }
                    }
                }

                // Log x/y/speed distributions (optional)
                if (dist_log) {
                    bool dist_do_log = false;
                    if (cli_dist_dt <= 0.0f) {
                        dist_do_log = true;
                    } else {
                        dist_do_log = (simulation_time + 1e-12 >= dist_next_log_abs_time);
                    }
                    if (dist_do_log) {
                        const int released = wall_is_released ? 1 : 0;
                        const double t_rel = (released && wall_release_time >= 0.0) ? (simulation_time - wall_release_time) : NAN;
                        log_distributions(simulation_time, t_rel, released);
                        if (cli_dist_dt > 0.0f) {
                            dist_next_log_abs_time += (double)cli_dist_dt;
                        }
                    }
                }

                // Szilard: log mutual information + work/heat at a fixed cadence (absolute time).
                // Do not gate this on "wall_is_released" because Szilard holds the divider
                // at different phases (and we still want a continuous trace).
                if (sz_log) {
                    if (sz_log_reset_requested) {
                        next_sz_log_abs_time = simulation_time;
                        sz_log_reset_requested = 0;
                    }
                    bool do_sz_log = false;
                    if (cli_szilard_log_dt <= 0.0) {
                        do_sz_log = true;
                    } else {
                        do_sz_log = (simulation_time + 1e-12 >= next_sz_log_abs_time);
                    }
                    if (do_sz_log) {
                        int nL0=0,nL1=0,nR0=0,nR1=0;
                        szilard_species_side_counts(&nL0,&nL1,&nR0,&nR1);
                        const double purity_L = szilard_side_purity(nL0, nL1);
                        const double purity_R = szilard_side_purity(nR0, nR1);
                        const double purity_mean = 0.5 * (purity_L + purity_R);
                        const double sep_score = szilard_separation_score(nL0, nL1, nR0, nR1);
                        int mem0=0,mem1=0,memU=0;
                        szilard_memory_counts(&mem0,&mem1,&memU);
                        int xy_match=0,xy_mismatch=0,xy_unknown=0;
                        szilard_xy_agreement_counts(&xy_match,&xy_mismatch,&xy_unknown);
                        const double xy_mismatch_frac = szilard_xy_mismatch_fraction();
                        const double partial_uncertainty_frac = szilard_partial_uncertainty_fraction_now();
                        double mean_speed = 0.0, std_speed = 0.0;
                        szilard_speed_stats(&mean_speed, &std_speed);
                        const double hot_speed_threshold = sz_stage1_baseline_mode
                                                         ? NAN
                                                         : sz_gate_hot_speed_mult * ((mean_speed > 1e-12) ? mean_speed : fmax(std_speed, 1e-12));
                        const double flip_tail_threshold = sz_stage1_baseline_mode
                                                         ? NAN
                                                         : mean_speed + sz_xflip_tail_sigma * ((std_speed > 1e-9) ? std_speed : fmax(mean_speed * 0.25, 1e-9));
                        int memL0=0,memL1=0,memLU=0,memR0=0,memR1=0,memRU=0;
                        szilard_memory_side_counts(&memL0,&memL1,&memLU,&memR0,&memR1,&memRU);
                        const double I_mem_nats = szilard_mutual_information_3state_nats(memL0,memL1,memLU,memR0,memR1,memRU);
                        const double I_mem_bits = I_mem_nats / log(2.0);
                        const double I_nats = szilard_mutual_information_nats(nL0,nL1,nR0,nR1);
                        const double I_bits = I_nats / log(2.0);
                        const double kBT_ref = szilard_reference_kbt();
                        const EDMD_Params* ep0 = (g_edmd ? edmd_params(g_edmd) : NULL);
                        const double gate0_x = (ep0 && ep0->divider_count > 0) ? ep0->divider_x[0] : ((double)wall_x - (double)XW1);
                        const double gate0_vx = (ep0 && ep0->divider_count > 0) ? ep0->divider_vx[0] : (double)vx_wall;
                        const double gate0_th = (ep0 && ep0->divider_count > 0) ? ep0->divider_thickness[0] : (wall_enabled ? (double)WALL_THICKNESS : 0.0);
                        const int gate0_mode = (ep0 ? ep0->divider_gate_mode[0] : 0);
                        const int gate0_sp   = (ep0 ? ep0->divider_gate_species[0] : 0);
                        const int gate0_use_mem = sz_perm_use_memory;
                        const int gate0_tgt  = (ep0 ? ep0->divider_gate_target_side[0] : 0);

                        const double gate1_x = (ep0 && ep0->divider_count > 1) ? ep0->divider_x[1] : ((double)sz_interactive_target_px - (double)XW1);
                        const double gate1_vx = (ep0 && ep0->divider_count > 1) ? ep0->divider_vx[1] : 0.0;
                        const double gate1_th = (ep0 && ep0->divider_count > 1) ? ep0->divider_thickness[1] : 0.0;
                        const int gate1_mode = (ep0 ? ep0->divider_gate_mode[1] : 0);
                        const int gate1_sp   = (ep0 ? ep0->divider_gate_species[1] : 0);
                        const int gate1_use_mem = sz_fixed_gate_use_memory;
                        const int gate1_tgt  = (ep0 ? ep0->divider_gate_target_side[1] : 0);

                        int nLg0=0,nLg1=0,nRg0=0,nRg1=0;
                        szilard_species_side_counts_at_split((double)XW1 + gate0_x, &nLg0,&nLg1,&nRg0,&nRg1);
                        const double purity_gate_L = szilard_side_purity(nLg0, nLg1);
                        const double purity_gate_R = szilard_side_purity(nRg0, nRg1);
                        const double purity_gate_mean = 0.5 * (purity_gate_L + purity_gate_R);
                        const double sep_gate_score = szilard_separation_score(nLg0, nLg1, nRg0, nRg1);

                        int memGateL0=0,memGateL1=0,memGateLU=0,memGateR0=0,memGateR1=0,memGateRU=0;
                        szilard_memory_side_counts_at_split((float)((double)XW1 + gate0_x),
                                                            &memGateL0,&memGateL1,&memGateLU,&memGateR0,&memGateR1,&memGateRU);
                        const double I_mem_gate_nats = szilard_mutual_information_3state_nats(
                            memGateL0,memGateL1,memGateLU,memGateR0,memGateR1,memGateRU
                        );
                        const double I_mem_gate_bits = I_mem_gate_nats / log(2.0);
                        const double I_gate_nats = szilard_mutual_information_nats(nLg0,nLg1,nRg0,nRg1);
                        const double I_gate_bits = I_gate_nats / log(2.0);
                        const double F_x = szilard_available_free_energy_from_nats(I_nats, nL0 + nL1 + nR0 + nR1);
                        const double F_x_gate = szilard_available_free_energy_from_nats(I_gate_nats, nLg0 + nLg1 + nRg0 + nRg1);
                        const double F_y = szilard_available_free_energy_from_nats(I_mem_nats, memL0 + memL1 + memLU + memR0 + memR1 + memRU);
                        const double F_y_gate = szilard_available_free_energy_from_nats(
                            I_mem_gate_nats,
                            memGateL0 + memGateL1 + memGateLU + memGateR0 + memGateR1 + memGateRU
                        );

                        const double pistR_x = (ep0 ? ep0->pistonR_x : ((double)piston_right_x - (double)XW1));
                        const double pistR_vx = (ep0 ? ep0->pistonR_vx : (double)vx_piston_right);

                        const double wG0 = (g_edmd ? edmd_work_divider_i(g_edmd, 0) : 0.0);
                        const double wG1 = (g_edmd ? edmd_work_divider_i(g_edmd, 1) : 0.0);
                        const double wD = wG0 + wG1;
                        const double wPL = (g_edmd ? edmd_work_pistonL(g_edmd) : 0.0);
                        const double wPR = (g_edmd ? edmd_work_pistonR(g_edmd) : 0.0);
                        const double wP  = wPL + wPR;
                        const double wT  = wD + wP;
                        const double qB = (g_edmd ? edmd_heat_bath(g_edmd) : 0.0);
                        const double eext = wT + qB;
                        fprintf(sz_log,
                                "%lld,%d,%s,%.6f,%d,%d,"
                                "%.6f,%.6f,%.6f,%d,%d,%d,%d,"
                                "%.6f,%.6f,%.6f,%d,%d,%d,%d,"
                                "%d,"
                                "%.6f,%.6f,"
                                "%d,%d,%d,%d,"
                                "%d,%d,%d,%d,"
                                "%.12g,%.12g,%.12g,%.12g,"
                                "%.12g,%.12g,%.12g,%.12g,"
                                "%d,%d,%d,"
                                "%d,%d,%d,%.12g,"
                                "%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,"
                                "%d,%d,%d,%d,%d,%d,"
                                "%d,%d,%d,%d,%d,%d,"
                                "%.12g,%.12g,"
                                "%.12g,%.12g,"
                                "%.12g,%.12g,"
                                "%.12g,%.12g,"
                                "%.12g,%.12g,%.12g,%.12g,%.12g,"
                                "%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,"
                                "%.12g,%.12g,%.12g,"
                                "%.12g,%.12g,%.12g,%.12g,"
                                "%.12g,%.12g,%.12g\n",
                                ++sz_log_index, sz_log_branch_id, sz_log_branch_label, simulation_time, sz_interactive_phase, particles_active,
                                gate0_x / (double)PIXELS_PER_SIGMA, gate0_vx / (double)PIXELS_PER_SIGMA, gate0_th / (double)PIXELS_PER_SIGMA, gate0_mode, gate0_sp, gate0_use_mem, gate0_tgt,
                                gate1_x / (double)PIXELS_PER_SIGMA, gate1_vx / (double)PIXELS_PER_SIGMA, gate1_th / (double)PIXELS_PER_SIGMA, gate1_mode, gate1_sp, gate1_use_mem, gate1_tgt,
                                sz_partial_obs_enabled,
                                pistR_x / (double)PIXELS_PER_SIGMA, pistR_vx / (double)PIXELS_PER_SIGMA,
                                nL0,nL1,nR0,nR1,
                                nLg0,nLg1,nRg0,nRg1,
                                purity_L,purity_R,purity_mean,sep_score,
                                purity_gate_L,purity_gate_R,purity_gate_mean,sep_gate_score,
                                mem0,mem1,memU,
                                xy_match,xy_mismatch,xy_unknown,xy_mismatch_frac,
                                partial_uncertainty_frac,sz_partial_obs_peak_uncertainty_frac,hot_speed_threshold,flip_tail_threshold,mean_speed,std_speed,
                                memL0,memL1,memLU,memR0,memR1,memRU,
                                memGateL0,memGateL1,memGateLU,memGateR0,memGateR1,memGateRU,
                                I_mem_nats,I_mem_bits,
                                I_mem_gate_nats,I_mem_gate_bits,
                                I_nats,I_bits,
                                I_gate_nats,I_gate_bits,
                                kBT_ref,F_x,F_x_gate,F_y,F_y_gate,
                                sz_load_fx_start,sz_load_fx_drop,sz_load_z,sz_load_energy,sz_load_work,sz_load_qdiss,
                                wG0,wG1,wD,
                                wPL,wPR,wP,wT,
                                sz_memory_erase_heat_min,qB,eext);
                        fflush(sz_log);
                        if (cli_szilard_log_dt > 0.0) next_sz_log_abs_time += cli_szilard_log_dt;
                    }
                }

                // logging after release
                if (wall_is_released) {
                    double time_after_release = simulation_time - wall_release_time;

                    // ##CHRIS: Unified CSV logging semantics tuned for your workflow
                    //
                    //  - TIME / RK4 modes:
                    //        cli_output_dt <= 1.0  → "max resolution" logging is handled
                    //                                inside update_particles_with_substepping()
                    //                                (per substep). We skip here to avoid
                    //                                duplicate rows.
                    //        cli_output_dt  > 1.0  → coarse logging at fixed σ-time intervals.
                    //
                    //  - EDMD / EDMD_HYBRID modes:
                    //        cli_output_dt <= 1.0  → "max resolution" logging: log once per
                    //                                EDMD step in this loop.
                    //        cli_output_dt  > 1.0  → coarse logging at fixed σ-time intervals.
                    bool do_log = false;
                    if (sim_mode == MODE_EDMD || sim_mode == MODE_EDMD_HYBRID) {
                        if (cli_output_dt <= 1.0f) {
                            // Per-step logging for EDMD / hybrid
                            do_log = true;
                        } else {
                            // Interval-based logging for EDMD / hybrid
                            do_log = (time_after_release + 1e-12 >= next_wall_log_rel_time);
                        }
                    } else {
                        // TIME / RK4: substep logger owns high-res (output_dt <= 1.0)
                        if (cli_output_dt > 1.0f) {
                            do_log = (time_after_release + 1e-12 >= next_wall_log_rel_time);
                        }
                    }
                    if (do_log) {
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

                            float packing_unitless     = particles_active * particle_area_unitless / box_area_unitless;
                            float packing_actual_bound = particles_active * particle_area_actual / box_area_scaled;
                            float packing_actual_excl  = particles_active * particle_area_actual / (box_area_scaled - wall_buffer_area);

                            fprintf(wall_log, "%.3f, %.3f, %.3f, %d, %d, %.6f, %.6f, %.6f\n",
                                    time_after_release, wall_x_sigma, disp,
                                    left_particles, right_particles,
                                    packing_unitless, packing_actual_bound, packing_actual_excl);
                        }
                        // ##CHRIS: Flush to ensure data is written immediately
                        fflush(wall_log);
                        // Advance target only for interval-based logging
                        if (cli_output_dt > 1.0f) {
                            next_wall_log_rel_time += (double)cli_output_dt;
                        }
                    }
                    // Optional auto-stop after requested time span (post release)
                    if (cli_time_span > 0.0f && time_after_release >= (double)cli_time_span) {
                        running = 0; // exit main loop cleanly
                    }
                }
            }
            accumulator -= fixed_dt_runtime;
        }

        // === RENDERING ===
        // ##CHRIS: Skip all rendering in headless mode
        if (!cli_headless) {
            // GUI-only: stash the leftover accumulator time for smooth boundary extrapolation in histograms.
            // (Physics state remains at simulation_time; this only affects visualization.)
            gui_render_dt_rem = accumulator;

            draw_clear_screen();
            gui_begin_scene_view();
            draw_coordinate_system(_renderer);
            draw_simulation_boundary();
            render_particles();
            render_pistons();
            if (wall_enabled) draw_wall();
            // Draw spring visualization on leftmost segment
            render_left_spring();

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
                float xl = left_bounds[seg];
                float xr = right_bounds[seg];
                int seg_px_w = (int)fmaxf(0.0f, (xr - xl) - 8.0f);
                if (seg_px_w < 40) continue; // too narrow to label cleanly

                SDL_Color color = segment_palette[seg % palette_size];

                // Choose a compact label when the compartment is narrow to avoid overlaps.
                // Also show a lightweight per-compartment temperature estimate for live tracking.
                double ke_val = (segment_ke) ? segment_ke[seg] : 0.0;
                double temp_val = 0.0;
                if (segment_temperature) {
                    temp_val = segment_temperature[seg];
                } else if (count > 0) {
                    temp_val = ke_val / (count * K_B);
                }
                // Live packing fraction per compartment (nominal geometric area).
                // Uses the current compartment width between boundaries (including piston motion).
                const double r_sigma = (double)PARTICLE_RADIUS / (double)PIXELS_PER_SIGMA;
                const double area_particle = M_PI * r_sigma * r_sigma;
                const double seg_w_sigma = fmax(1e-9, (double)(xr - xl) / (double)PIXELS_PER_SIGMA);
                const double area_seg = seg_w_sigma * (double)HEIGHT_UNITS;
                const double eta_seg = (area_seg > 0.0) ? ((double)count * area_particle / area_seg) : 0.0;

                char seg_buffer_1[96];
                char seg_buffer_2[64];
                if (seg_px_w < 140) {
                    snprintf(seg_buffer_1, sizeof(seg_buffer_1), "N=%d", count);
                } else {
                    snprintf(seg_buffer_1, sizeof(seg_buffer_1), "Box %d: N=%d", seg + 1, count);
                }
                if (seg_px_w >= 160) {
                    snprintf(seg_buffer_2, sizeof(seg_buffer_2), "T=%.2f  eta=%.3f", temp_val, eta_seg);
                } else {
                    snprintf(seg_buffer_2, sizeof(seg_buffer_2), "T=%.2f", temp_val);
                }

                int text_y = YW1 + 6;
                int text_x = (int)xl + 4;
                int tw = 0, th = 0;
                if (TTF_SizeText(font, seg_buffer_1, &tw, &th) == 0) {
                    int max_x = (int)xr - 4 - tw;
                    if (text_x > max_x) text_x = max_x;
                    if (text_x < (int)xl + 2) text_x = (int)xl + 2;
                }
                draw_text(_renderer, font, seg_buffer_1, text_x, text_y, color);

                int twT = 0, thT = 0;
                int xT = (int)xl + 4;
                if (TTF_SizeText(font, seg_buffer_2, &twT, &thT) == 0) {
                    int max_xT = (int)xr - 4 - twT;
                    if (xT > max_xT) xT = max_xT;
                    if (xT < (int)xl + 2) xT = (int)xl + 2;
                }
                draw_text(_renderer, font, seg_buffer_2, xT, text_y + 18, color);

                // Optional third line (KE) only when there is enough space.
                if (seg_px_w >= 220) {
                    // In stacked distribution mode we show distribution metrics inside the histogram panels
                    // (less clutter inside the simulation box).
                    if (gui_dist_mode != 3) {
                        char ke_buffer[80];
                        snprintf(ke_buffer, sizeof(ke_buffer), "KE=%.2e", ke_val);
                        int tw2 = 0, th2 = 0;
                        int x2 = (int)xl + 4;
                        if (TTF_SizeText(font, ke_buffer, &tw2, &th2) == 0) {
                            int max_x2 = (int)xr - 4 - tw2;
                            if (x2 > max_x2) x2 = max_x2;
                            if (x2 < (int)xl + 2) x2 = (int)xl + 2;
                        }
                        draw_text(_renderer, font, ke_buffer, x2, text_y + 36, color);
                    }
                }
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

		            if (gui_dist_mode == 1) {
		                render_velocity_histograms();
		            } else if (gui_dist_mode == 2) {
		                render_velocity_histograms_by_segment();
		            } else if (gui_dist_mode == 3) {
		                render_segment_distributions_stacked();
		            }
            gui_end_scene_view();
            render_energy_hud(_renderer, font);
            render_live_controls(_renderer, font_small ? font_small : font);
            float T_measured = compute_measured_temperature_from_ke();
            char temp_label[128];
            snprintf(temp_label, sizeof(temp_label), "T_measured: %.2f", T_measured);
            //draw_text(_renderer, font, temp_label, XW1 + 5, YW1 + HEIGHT_UNITS*PIXELS_PER_SIGMA, yellow);

            SDL_RenderPresent(_renderer);
            SDL_Delay(5);
        } // ##CHRIS: End of rendering block (skipped in headless mode)
    }

    fclose(logFile);
    fclose(wall_log);
    if (sz_log) fclose(sz_log);
    if (et_wall_log) fclose(et_wall_log);
    if (dist_log) fclose(dist_log);
    if (dist_hist_x) { free(dist_hist_x); dist_hist_x = NULL; }
    if (dist_hist_y) { free(dist_hist_y); dist_hist_y = NULL; }
    if (dist_hist_v) { free(dist_hist_v); dist_hist_v = NULL; }
    if (cli_dist_log_path) { free(cli_dist_log_path); cli_dist_log_path = NULL; }
}

/////////////// SINGLE-TEST MODE (HEADLESS) /////////////
// Run simulation headless (no SDL window) but output to main directory for quick testing
void run_single_test_headless() {
    printf("🔬 Single-test mode: Running headless simulation...\n");
    printf("   Output will be written to wall_position.csv and energy_log.csv\n");
    printf("   Run 'python3 wall_x_FFT.py' to analyze results\n\n");

    // Initialize simulation (same as interactive mode)
    initialize_simulation();

    // Initialize EDMD if needed (MODE_EDMD or MODE_EDMD_HYBRID)
    if (sim_mode == MODE_EDMD || sim_mode == MODE_EDMD_HYBRID) {
        if (g_edmd) { edmd_destroy(g_edmd); g_edmd = NULL; }
        EDMD_Params prm = {0};
        prm.boxW = (double)(XW2 - XW1);
        prm.boxH = (double)(YW2 - YW1);
        prm.radius = (double)PARTICLE_RADIUS;
        prm.N = particles_active;
        prm.cell_size = 0.0; // auto
        {
            int dcount = (num_internal_walls > 0 ? num_internal_walls : 0);
            if (dcount > EDMD_MAX_DIVIDERS) dcount = EDMD_MAX_DIVIDERS;
            prm.divider_count = dcount;
            int primary_idx = primary_wall_index;
            if (primary_idx < 0 || primary_idx >= dcount) primary_idx = 0;
            int extra_idx = 0;
            for (int w = 0; w < dcount; ++w) {
                double pos_px = (all_wall_positions ? (double)all_wall_positions[w] : (double)wall_x);
                prm.divider_x[w] = pos_px - (double)XW1;
                prm.divider_thickness[w] = (double)WALL_THICKNESS;
                double mf = (double)(WALL_MASS / PARTICLE_MASS);
                if (cli_wall_masses_set && w < (int)(sizeof(cli_wall_masses) / sizeof(cli_wall_masses[0])) && cli_wall_masses[w] > 0.0f) {
                    mf = (double)cli_wall_masses[w];
                }
                prm.divider_mass[w] = mf;
                double vxw = 0.0;
                if (w == primary_idx) {
                    vxw = (double)vx_wall;
                } else if (extra_wall_velocity && extra_idx < extra_wall_count) {
                    vxw = (double)extra_wall_velocity[extra_idx];
                    extra_idx++;
                }
                prm.divider_vx[w] = vxw;
            }
        }
        prm.heatbath_enabled = heatbath_enabled;
        prm.heatbath_temperature = (double)heatbath_temperature;
        prm.thermal_wall_mode = thermal_wall_mode;
        prm.mb_overshoot_factor = (double)mb_overshoot_factor;
        prm.stability_window_percent = (double)stability_window_percent;
        prm.particle_mass = (double)PARTICLE_MASS;
        prm.kB = (double)kB_effective();
        // EDMD harmonic spring wall (single-test/headless): enable only when spring output is requested.
        if (energy_measurement.spring_enabled && prm.divider_count > 0) {
            int primary_idx = primary_wall_index;
            if (primary_idx < 0 || primary_idx >= prm.divider_count) primary_idx = 0;
            const double xeq_px = (energy_measurement.equilibrium_position != 0.0f)
                                    ? (double)energy_measurement.equilibrium_position
                                    : (double)wall_x;
            prm.divider_k[primary_idx] = (double)energy_measurement.spring_constant;
            prm.divider_xeq[primary_idx] = xeq_px - (double)XW1;
        }

        g_edmd = edmd_create(&prm);
        if (!g_edmd) { fprintf(stderr, "EDMD create failed.\n"); exit(1); }

        // Load current positions/velocities into EDMD state
        EDMD_Particle* P = (EDMD_Particle*)edmd_particles(g_edmd);
        for (int i = 0; i < particles_active; i++) {
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

        edmd_reschedule_all(g_edmd);
        edmd_reset_work(g_edmd);
        printf("[EDMD] initialized: N=%d box=(%.1f,%.1f) R=%.3f\n", prm.N, prm.boxW, prm.boxH, prm.radius);
    }

    // Open output files
    FILE *logFile = fopen("energy_log.csv", "w");
    FILE *wall_log = fopen("wall_position.csv", "w");
    if (!logFile || !wall_log) {
        fprintf(stderr, "❌ Failed to open log files for writing\n");
        exit(1);
    }

    // Write headers
    fprintf(wall_log, "Time, Wall_X, Displacement(σ), Left_Count, Right_Count\n");
    fprintf(logFile, "KE_particles,KE_wall,KE_total,T,S_v,F\n");

    // Write run parameters for Python analysis
    write_run_params_json();

    // Set initial state
    steps_elapsed = 0;
    simulation_time = 0.0;
    wall_release_time = -1.0;
    wall_hold_enabled = !cli_auto_release_wall;
    wall_is_released = cli_auto_release_wall;
    wall_enabled = 1;
    wall_x_old = wall_x;
    missed_collision_events = 0;

    if (cli_auto_release_wall) {
        wall_release_time = simulation_time;
        printf("🚀 Wall auto-released at t=%.3f\n", wall_release_time);
    }

    // Determine max steps to run (default 3000 if not overridden)
    int max_steps = (cli_override_num_steps > 0) ? cli_override_num_steps : 3000;
    printf("Running for %d steps...\n", max_steps);
    printf("Fixed timestep: %.6f σ-time units\n", fixed_dt_runtime / PIXELS_PER_SIGMA);

    // Progress reporting
    int log_interval = (cli_output_dt <= 1.0f) ? 1 : 100;  // Log every step if high-res, else every 100
    int progress_interval = max_steps / 20;  // Print progress 20 times (≈5% increments)
    if (progress_interval < 1) progress_interval = 1;

    double next_wall_log_rel_time = 0.0;
    int left_particles = 0, right_particles = 0;

    // Main simulation loop
    while (steps_elapsed < max_steps) {
        // Advance simulation (EDMD mode)
        if (sim_mode == MODE_EDMD) {
            wall_x_old = wall_x;

            const EDMD_Params* ep_state = edmd_params(g_edmd);
            int dcount = ep_state->divider_count;
            if (!wall_enabled) dcount = 0;
            if (dcount > EDMD_MAX_DIVIDERS) dcount = EDMD_MAX_DIVIDERS;
            int primary_idx = primary_wall_index;
            if (primary_idx < 0 || primary_idx >= dcount) primary_idx = 0;

            double div_x[EDMD_MAX_DIVIDERS];
            double div_vx[EDMD_MAX_DIVIDERS];
            double div_mass[EDMD_MAX_DIVIDERS];
            double div_th[EDMD_MAX_DIVIDERS];
            int extra_idx = 0;

            for (int w = 0; w < dcount; ++w) {
                div_x[w] = ep_state->divider_x[w];
                div_vx[w] = ep_state->divider_vx[w];
                div_th[w] = (double)WALL_THICKNESS;
                double mf = (double)(WALL_MASS / PARTICLE_MASS);
                if (cli_wall_masses_set && w < (int)(sizeof(cli_wall_masses) / sizeof(cli_wall_masses[0])) && cli_wall_masses[w] > 0.0f) {
                    mf = (double)cli_wall_masses[w];
                }
                // Hold walls (infinite mass) before release.
                bool hold_active = false;
                if (w == primary_idx) {
                    hold_active = wall_hold_enabled && !wall_is_released;
                } else if (extra_wall_hold_enabled && extra_wall_is_released) {
                    int j = (w < primary_idx) ? w : (w - 1);
                    if (j >= 0 && j < extra_wall_count) {
                        hold_active = extra_wall_hold_enabled[j] && !extra_wall_is_released[j];
                    }
                }
                if (hold_active) {
                    div_mass[w] = 0.0;
                    div_vx[w] = 0.0;
                } else {
                    div_mass[w] = mf;
                }

                if (w == primary_idx) {
                    wall_x = (double)XW1 + div_x[w];
                    vx_wall = div_vx[w];
                } else if (extra_wall_positions && extra_wall_velocity && extra_idx < extra_wall_count) {
                    extra_wall_positions[extra_idx] = (float)((double)XW1 + div_x[w]);
                    extra_wall_velocity[extra_idx] = (float)div_vx[w];
                    extra_idx++;
                }
            }

            // In pure EDMD mode, the spring wall is handled by the EDMD core (harmonic divider)
            // when --eff-output=spring is used. Do not apply timestep impulses here.

            edmd_config_dividers(g_edmd, dcount, div_x, div_th);
            edmd_set_divider_motions(g_edmd, dcount, div_mass, div_vx);
            edmd_reschedule_all(g_edmd);

            double t0 = edmd_time(g_edmd);
            edmd_advance_to(g_edmd, t0 + (double)fixed_dt_runtime);
            edmd_divider_resolve_overlaps(g_edmd);

            // Copy particle state back
            const EDMD_Particle* P = edmd_particles(g_edmd);
            for (int i = 0; i < particles_active; i++) {
                X[i] = (double)XW1 + (double)P[i].x;
                Y[i] = (double)YW1 + (double)P[i].y;
                Vx[i] = (double)P[i].vx;
                Vy[i] = (double)P[i].vy;
            }

            const EDMD_Params* ep = edmd_params(g_edmd);
            int dcount2 = ep->divider_count;
            if (!wall_enabled) dcount2 = 0;
            if (dcount2 > EDMD_MAX_DIVIDERS) dcount2 = EDMD_MAX_DIVIDERS;
            int primary_idx2 = primary_wall_index;
            if (primary_idx2 < 0 || primary_idx2 >= dcount2) primary_idx2 = 0;
            int extra_idx2 = 0;
            for (int w = 0; w < dcount2; ++w) {
                if (w == primary_idx2) {
                    wall_x = (double)XW1 + ep->divider_x[w];
                    vx_wall = (float)ep->divider_vx[w];
                } else if (extra_wall_positions && extra_wall_velocity && extra_idx2 < extra_wall_count) {
                    extra_wall_positions[extra_idx2] = (float)((double)XW1 + ep->divider_x[w]);
                    extra_wall_velocity[extra_idx2] = (float)ep->divider_vx[w];
                    extra_idx2++;
                }
            }

            // Recompute segment statistics
            recompute_segment_stats_counts_and_temperature();
            if (segment_counts && segment_count >= 2) {
                left_particles = segment_counts[0];
                right_particles = segment_counts[1];
            }
        }

        // Update simulation time
        simulation_time += fixed_dt_runtime / PIXELS_PER_SIGMA;
        steps_elapsed++;

        // Establish release time once
        if (wall_is_released && wall_release_time < 0.0) {
            wall_release_time = simulation_time;
            if (!cli_quiet) {
                printf("✅ Wall released at t = %.3f\n", wall_release_time);
            }
            next_wall_log_rel_time = 0.0;
        }

        // Log data after wall is released
        if (wall_is_released) {
            double time_after_release = simulation_time - wall_release_time;

            // Determine if we should log this step
            bool do_log = false;
            if (cli_output_dt <= 1.0f) {
                do_log = (steps_elapsed % log_interval == 0);
            } else {
                do_log = (time_after_release + 1e-12 >= next_wall_log_rel_time);
            }

            if (do_log) {
                double wall_x_sigma = wall_x / PIXELS_PER_SIGMA;
                double disp = wall_x_sigma - (XW1 + XW2) / (2.0 * PIXELS_PER_SIGMA);

                fprintf(wall_log, "%.3f, %.3f, %.3f, %d, %d\n",
                        time_after_release, wall_x_sigma, disp,
                        left_particles, right_particles);
                fflush(wall_log);

                if (cli_output_dt > 1.0f) {
                    next_wall_log_rel_time += (double)cli_output_dt;
                }
            }
        }

        // Log energy periodically
        if (steps_elapsed % 100 == 0) {
            double ke_particles = kinetic_energy();
            double ke_wall = kinetic_energy_walls_total();
            double T_current = (particles_active > 0) ? ke_particles / (particles_active * kB_effective()) : 0.0;
            fprintf(logFile, "KE_particles: %.8e, KE_wall: %.8e, KE_total: %.8e, T: %.8e, S_v: %.8e, F: nan\n",
                    ke_particles, ke_wall, ke_particles + ke_wall, T_current, 0.0);
            fflush(logFile);
        }

        // Print progress
        if (steps_elapsed % progress_interval == 0) {
            float percent = (100.0f * steps_elapsed) / max_steps;
            printf("  Progress: %.1f%% (step %d/%d, t=%.2f)\n",
                   percent, steps_elapsed, max_steps, simulation_time);
        }

        // Auto-stop if time span exceeded
        if (cli_time_span > 0.0f && wall_is_released) {
            double time_after_release = simulation_time - wall_release_time;
            if (time_after_release >= (double)cli_time_span) {
                printf("⏱️  Time span limit reached (%.2f σ-time units)\n", cli_time_span);
                break;
            }
        }
    }

    printf("\n✅ Simulation complete!\n");
    printf("   Total steps: %d\n", steps_elapsed);
    printf("   Final time: %.2f σ-time units\n", simulation_time);
    printf("   Output written to: wall_position.csv, energy_log.csv\n");
    printf("\n   Next step: Run analysis with:\n");
    printf("   python3 wall_x_FFT.py\n\n");

    // Close files
    fclose(wall_log);
    fclose(logFile);
}











/////////////// MAIN FUNCTION /////////////
// Main Function
int main(int argc, char* argv[]) {

    g_main_argc = argc;
    g_main_argv = argv;

    parse_cli_options(argc, argv);
    init_output_dirs_from_argv0(argv[0]);
    // Record the invoking command in the current working directory as well
    // (useful when external scripts set cwd to a per-run output folder).
    write_command_md_to_dir(".", argc, argv);
    if (cli_seed_set) {
        srand(cli_seed);
    }
    initialize_simulation_dimensions();     // sets SIM_WIDTH, SIM_HEIGHT, XW2 etc.
    initialize_time_scale();
    initialize_simulation_parameters();     // sets wall and piston regarding simulation dimensions
    // Set dt runtime based on preset and temperature
    update_dt_runtime_for_temperature();
    // Initialize energy/efficiency measurement if enabled (gated inside).
    initialize_energy_measurement();
    print_simulation_info();  // Print simulation parameters

    // ##CHRIS: Skip SDL/TTF initialization when running headless
    if (!cli_headless) {
        initSDL();
        TTF_Init();
        //kiss_fftr_cfg fft_cfg = kiss_fftr_alloc(FFT_SIZE, 0, NULL, NULL);
        //initialize_simulation();

        font = TTF_OpenFont("Roboto-Regular.ttf", 18);
        if (!font) {
            printf("Failed to load font: %s\n", TTF_GetError());
            return 1;
        }

        // Smaller font for dense in-panel histogram stats (avoids cut-off on lower vertical resolutions).
        font_small = TTF_OpenFont("Roboto-Regular.ttf", 12);
        if (!font_small && !cli_quiet) {
            printf("⚠️  Failed to load small font (continuing): %s\n", TTF_GetError());
        }
    } else {
        if (!cli_quiet) {
            printf("Running in headless mode (no GUI, CSV output enabled)\n");
        }
    }

  
   // printf("wall_x = %f, wall_enabled = %d\n", wall_x, wall_enabled);

    // ##CHRIS: Initialize particle-life experiment if selected
    if (cli_experiment_preset == EXPERIMENT_PRESET_PARTICLELIFE) {
        pl_initialize();
    }

    // Single-test mode takes priority over experiment modes
    if (cli_single_test_mode) {
        run_single_test_headless();
    } else if (!cli_force_no_experiments && cli_experiment_preset == EXPERIMENT_PRESET_SIMPLE_BOX_PREDICTION) {
        run_simple_prediction_box_experiment();
    } else if (!cli_force_no_experiments && cli_experiment_preset == EXPERIMENT_PRESET_ENERGY_TRANSFER) {
        run_energy_transfer_experiment();
    } else if (!cli_force_no_experiments && cli_experiment_preset == EXPERIMENT_PRESET_SZILARD_ENGINE) {
        run_szilard_engine_experiment();
    } else if (enable_speed_of_sound_experiments && !cli_force_no_experiments) {
        run_speed_of_sound_experiments();
    } else {
        simulation_loop();  // SDL interactive
    }
    

    // Print the first few velocities for testing
    if (!cli_quiet) {
        for (int i = 0; i < 10; i++) {
            printf("Particle %d: Vx = %f, Vy = %f\n", i, Vx[i], Vy[i]);
        }
    }

    // ##CHRIS: Clean up SDL/TTF only if they were initialized
	    if (!cli_headless) {
	        SDL_DestroyRenderer(_renderer);
	        SDL_DestroyWindow(_window);
	        //free(fft_cfg);
            if (font_small) { TTF_CloseFont(font_small); font_small = NULL; }
	        TTF_CloseFont(font);
	        TTF_Quit();
	        SDL_Quit();
	    }
    free_internal_wall_arrays();
    free(cli_experiment_lengths); cli_experiment_lengths = NULL;
    free(cli_experiment_wall_masses); cli_experiment_wall_masses = NULL;
    if (cli_speed_sound_run_dir) { free(cli_speed_sound_run_dir); cli_speed_sound_run_dir = NULL; }
    if (cli_simple_box_run_dir) { free(cli_simple_box_run_dir); cli_simple_box_run_dir = NULL; }
    if (cli_simple_box_run_label) { free(cli_simple_box_run_label); cli_simple_box_run_label = NULL; }
    // ##CHRIS: Clean up particle-life data
    pl_cleanup();
    if (sz_species) { free(sz_species); sz_species = NULL; }
    if (sz_memory) { free(sz_memory); sz_memory = NULL; }
    if (sz_gate_labels) { free(sz_gate_labels); sz_gate_labels = NULL; }
    return 0;
}
