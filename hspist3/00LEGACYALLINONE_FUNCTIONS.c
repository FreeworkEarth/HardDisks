
#include <SDL2/SDL.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <SDL2/SDL_ttf.h>
#include "kissfft/kiss_fft.h"
#include "kissfft/kiss_fftr.h"

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
#define SUBSTEPPING 1 // BOOL
#define INIT_SCALING 0 // BOOL
#define INIT_SCALING_WALL_RELEASE_STABLE 0 // BOOL

// PARTICLE PARAMETERS
#define NUM_PARTICLES 100      // total number of particles
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
    #define K_B 1.0f
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
        #define SUBSTEPS 20
        #define WALL_HOLD_STEPS 10000
        #define PARTICLE_SCALE_PARAM 1.0f // scaling factor for particle size to make more space for particles



    #else
        #define FIXED_DT 0.0016f  // in normalized time units They don’t state the exact dt, but normalized simulations use dt = 10^{-3} often.
        #define PIXELS_PER_SIGMA 40.0f      // rendering scaleq
        #define TEMPERATURE 1.0f            //SET to 0 if want to ssimulte traversing wave
        #define SUBSTEPS 10
        #define WALL_HOLD_STEPS 10000
        #define PARTICLE_SCALE_PARAM 0.01f // scaling factor for particle size


    #endif
#else
    #if EXPERIMENT_MODE
        #define FIXED_DT 2e-4f
        #define SUBSTEPS 20
        #define WALL_HOLD_STEPS 10000
        #define PIXELS_PER_SIGMA 1.0f      // rendering scale
        #define PARTICLE_SCALE_PARAM 1.0f // scaling factor for particle size


    #else
        #define FIXED_DT 1.6e-3f
        #define SUBSTEPS 10
        #define WALL_HOLD_STEPS 5000
        #define PIXELS_PER_SIGMA 20.0f      // rendering scale
        #define PARTICLE_SCALE_PARAM 1.0f // scaling factor for particle size



    #endif
#endif


#define DIAMETER (2 * PARTICLE_RADIUS_UNIT)
#define PARTICLE_RADIUS PARTICLE_RADIUS_UNIT * PIXELS_PER_SIGMA// / PARTICLE_SCALE_PARAM)) // in meters

int num_steps = TIME_UNITS_SIMULATION / FIXED_DT;  // = 3,000,000 with 3000 units and 0.001 dt
bool wall_position_is_managed_externally = false;



// WALL PARAMETERS
#define WALL_MASS_FACTOR 50                              // choose factor relative to PARTICLE_MASS
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


// /////////   THICKNESSS
#if EXPERIMENT_MODE
    #define WALL_THICKNESS 1.0f
#else
    #define WALL_THICKNESS fmaxf(PARTICLE_RADIUS*2, 0.10f)  // always enough to prevent tunneling //WALL_THICKNESS ((PARTICLE_RADIUS * 3 / 1.0f) < 1.0f  ? 1.0f : (PARTICLE_RADIUS * 3 / 1.0f))
#endif

//////v   WEIGHTED WALL MASS
float wall_mass_runtime = PARTICLE_MASS * WALL_MASS_FACTOR;
#define WALL_MASS (wall_mass_runtime)

bool wall_hold_enabled = true;         // Enable or disable hold wall mode
bool wall_is_released = false;          // Has the user released the wall?
int wall_hold_steps = WALL_HOLD_STEPS;             // Hold the wall for this many steps
int steps_elapsed = 0;    




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
}


// Global variables for pistons
float piston_left_x = 0;
float piston_right_x;
float vx_piston_left = 0; // Slow movement for piston velocity
float vx_piston_right = 0; // Slow movement for piston velocity

// Wall parameters
float wall_x_old;
float wall_x;
float vx_wall = 0.0;  // Velocity of the wall
int wall_enabled = 0;  // 0 = Wall is disabled, 1 = Wall is enabled

void initialize_simulation_parameters() {
    piston_right_x = MAX_X;

    // Initialize wall in the center of the box
    //wall_x = (XW1 + XW2) / 2.0f;
    wall_x = XW1 + (L0_UNITS * PIXELS_PER_SIGMA);  // Convert L0 to pixel space
    wall_x_old = wall_x;
    vx_wall = 0.0f;

    piston_left_x = 0.0f;
    vx_piston_left = 0.0f;
    vx_piston_right = 0.0f;
}

// Global:
//FILE *tunnel_log = NULL;

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
float X[NUM_PARTICLES];
float Y[NUM_PARTICLES];
float X_old[NUM_PARTICLES];    // previous X positions

float V_init[NUM_PARTICLES];   // initial velocity
float Vx[NUM_PARTICLES];
float Vy[NUM_PARTICLES];
float Radius[NUM_PARTICLES];

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







void handle_left_to_right_collision(int i) {
    printf("💥 L→R Collision with particle %d | v_particle = %.4f, v_wall = %.4f → new_v_wall = %.4f\n",i, Vx[i], vx_wall, ((WALL_MASS - PARTICLE_MASS) * vx_wall + 2 * PARTICLE_MASS * Vx[i]) / (PARTICLE_MASS + WALL_MASS));
    float v_particle = Vx[i];
    float v_wall = vx_wall;
    float rel_vel = v_particle - v_wall;

    float wall_left = wall_x - WALL_THICKNESS/2;

    if (X[i] + Radius[i] >= wall_left && rel_vel > 0) {
        float new_v_particle = ((PARTICLE_MASS - WALL_MASS) * v_particle + 2 * WALL_MASS * v_wall) / (PARTICLE_MASS + WALL_MASS);
        float new_v_wall = ((WALL_MASS - PARTICLE_MASS) * v_wall + 2 * PARTICLE_MASS * v_particle) / (PARTICLE_MASS + WALL_MASS);

        Vx[i] = new_v_particle;
        vx_wall = new_v_wall;
        X_old[i] = X[i];
        X[i] = wall_left - Radius[i] - 0.01f;
    }
}



void handle_right_to_left_collision(int i) {
    printf("💥 R→L Collision with particle %d | v_particle = %.4f, v_wall = %.4f → new_v_wall = %.4f\n",i, Vx[i], vx_wall, ((WALL_MASS - PARTICLE_MASS) * vx_wall + 2 * PARTICLE_MASS * Vx[i]) / (PARTICLE_MASS + WALL_MASS));    float v_particle = Vx[i];
    float v_wall = vx_wall;
    float rel_vel = v_particle - v_wall;

    float wall_right = wall_x + WALL_THICKNESS/2;

    if (X[i] - Radius[i] <= wall_right && rel_vel < 0) {
        float new_v_particle = ((PARTICLE_MASS - WALL_MASS) * v_particle + 2 * WALL_MASS * v_wall) / (PARTICLE_MASS + WALL_MASS);
        float new_v_wall = ((WALL_MASS - PARTICLE_MASS) * v_wall + 2 * PARTICLE_MASS * v_particle) / (PARTICLE_MASS + WALL_MASS);

        Vx[i] = new_v_particle;
        vx_wall = new_v_wall;
        
        X_old[i] = X[i];
        X[i] = wall_right + Radius[i] + 0.01f;
    }
}





void check_wall_collisions(int i) {
    if (!wall_enabled) return;

    float wall_left_now = wall_x - WALL_THICKNESS / 2;
    float wall_right_now = wall_x + WALL_THICKNESS / 2;
    float wall_left_old = wall_x_old - WALL_THICKNESS / 2;
    float wall_right_old = wall_x_old + WALL_THICKNESS / 2;

    float particle_x = X[i];
    float particle_r = Radius[i];

    // --- Predictive check: Did wall move across particle's X in this frame? ---
    bool wall_swept_across_left = (particle_x + particle_r >= wall_left_now && particle_x + particle_r <= wall_left_old && vx_wall < 0);
    bool wall_swept_across_right = (particle_x - particle_r <= wall_right_now && particle_x - particle_r >= wall_right_old && vx_wall > 0);

    if (wall_swept_across_left) {

        handle_left_to_right_collision(i);
        return;
    }

    if (wall_swept_across_right) {
        handle_right_to_left_collision(i);
        return;
    }

    // Original checks (based on particle movement)
    float prev_left_edge = X_old[i] - Radius[i];
    float prev_right_edge = X_old[i] + Radius[i];
    float curr_left_edge = X[i] - Radius[i];
    float curr_right_edge = X[i] + Radius[i];

    if (prev_right_edge < wall_left_now && curr_right_edge >= wall_left_now && Vx[i] > 0) {
        handle_left_to_right_collision(i);
        return;
    }

    if (prev_left_edge > wall_right_now && curr_left_edge <= wall_right_now && Vx[i] < 0) {
        handle_right_to_left_collision(i);
        return;
    }

    // Fallback if stuck inside wall
    if (curr_right_edge >= wall_left_now && curr_left_edge < wall_right_now) {
        if (Vx[i] > 0) {
            handle_left_to_right_collision(i);
        } else if (Vx[i] < 0) {
            handle_right_to_left_collision(i);
        }
    }
}













////////// PARTICLE COLLISION FUNCTIONS ///////////




/////////// TUNNELING DETECTION AND FIXING /////////////
/*
Why does one-particle tunneling matter?

Let’s quantify it:

1. Momentum imbalance

A single particle crossing through the wall introduces a net change in momentum to the wall:
\Delta p = 2m \cdot v
If a particle tunnels from right to left, it gives the wall a kick leftward. The wall then keeps oscillating around a new center (shifted).

2. In hard-disk systems (mass = 1, velocities ~1):

Each particle’s impulse can be ≈ 1–3 units of momentum.
If your wall mass is M = 50, one event gives:
\Delta v_\text{wall} \approx \frac{2}{50} = 0.04
Even that small kick can slowly accumulate displacement or bias the FFT.

3. In energy spectrum:
	•	This doesn’t ruin the frequency of oscillation,
	•	But it skews the baseline, shifts the wall’s equilibrium, and masks symmetry-based phenomena (like reflection symmetry or balanced fluctuations).
    
Tunneling occurs when a particle moves so far in one timestep that it completely skips over the wall or another object. To avoid it:

\Delta x_{\text{particle}} = v \cdot \Delta t < \text{wall thickness}

So if you want to reduce tunneling, make sure:
\Delta t < \frac{\text{wall thickness}}{\text{max particle speed}}
*/




void detect_and_fix_tunnel(float *X, float *X_old, float *Vx, int *Side_old, int *Side_new, int num_particles, float wall_x, float wall_thickness) {
    float wall_left = wall_x - wall_thickness / 2;
    float wall_right = wall_x + wall_thickness / 2;

    FILE *tunnel_log = fopen("tunnel_events.csv", "a");  // Append mode to avoid overwriting
    if (!tunnel_log) {
        printf("Error opening tunnel_events.csv!\n");
        return;
    }

    for (int i = 0; i < num_particles; i++) {
        bool was_left = X_old[i] + Radius[i] < wall_left;
        bool now_right = X[i] - Radius[i] > wall_right;
        bool was_right = X_old[i] - Radius[i] > wall_right;
        bool now_left = X[i] + Radius[i] < wall_left;

        if ((was_left && now_right) || (was_right && now_left)) {
            printf("❌ FATAL: Particle %d tunneled through wall at step %d. X_old = %.4f → X = %.4f\n", i, steps_elapsed, X_old[i], X[i]);
            fprintf(tunnel_log, "%d, %d, %.6f, %.6f, %.6f\n", steps_elapsed, i, X[i], wall_x, Vx[i]);
            exit(1);  // 💣 Hard kill if simulation reaches unphysical state
        }
    }

    fclose(tunnel_log);
}