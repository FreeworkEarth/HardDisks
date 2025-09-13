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
        #define SUBSTEPS 20 // default, can be changed externally
        #define WALL_HOLD_STEPS 10000
        #define PARTICLE_SCALE_PARAM 1.0f // scaling factor for particle size to make more space for particles



    #else
        #define FIXED_DT 0.0016f  // in normalized time units They don’t state the exact dt, but normalized simulations use dt = 10^{-3} often.
        #define PIXELS_PER_SIGMA 40.0f      // rendering scaleq
        #define TEMPERATURE 1.0f            //SET to 0 if want to ssimulte traversing wave
        #define SUBSTEPS 10 // default, can be changed externally
        #define WALL_HOLD_STEPS 10000
        #define PARTICLE_SCALE_PARAM 0.01f // scaling factor for particle size


    #endif
#else
    #if EXPERIMENT_MODE
        #define FIXED_DT 2e-4f
        #define SUBSTEPS 20 // default, can be changed externally
        #define WALL_HOLD_STEPS 10000
        #define PIXELS_PER_SIGMA 1.0f      // rendering scale
        #define PARTICLE_SCALE_PARAM 1.0f // scaling factor for particle size


    #else
        #define FIXED_DT 1.6e-3f
        #define SUBSTEPS 10 // default, can be changed externally
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







///////////  SIMULALTION INFO ////////////
void print_simulation_info() {
    printf("\n========== Simulation Info ==========\n");
    printf("Mode: %s\n", MOLECULAR_MODE ? "MOLECULAR (Normalized units)" : "MACROSCOPIC (Real units)");
    printf("Experiment Mode: %s\n", EXPERIMENT_MODE ? "ON (fine dt/substeps)" : "OFF (default dt)");
    printf("Fixed Time Step (dt): %.2e\n", FIXED_DT);
    printf("TAU scaling factor: %.4f\n", TAU_time_scale_factor_to_molecular);
    printf("Wall hold mode: %s\n", wall_hold_enabled ? "Enabled" : "Disabled");
    printf("Wall hold steps: %d\n", wall_hold_steps);
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
    int left_particles = NUM_PARTICLES / 2;
    int right_particles = NUM_PARTICLES - left_particles;

    float wall_buffer = WALL_THICKNESS * 2;

    for (int i = 0; i < NUM_PARTICLES; i++) {
        int valid_position = 0;
        int attempts = 0;

        while (!valid_position) {
            if (i < left_particles) {
                X[i] = ((float)rand() / RAND_MAX) * (wall_x - wall_buffer - XW1) + XW1 + PARTICLE_RADIUS;
            } else {
                X[i] = ((float)rand() / RAND_MAX) * (XW2 - wall_x - wall_buffer) + wall_x + wall_buffer;
            }

            Y[i] = ((float)rand() / RAND_MAX) * (YW2 - YW1 - 2 * PARTICLE_RADIUS) + YW1 + PARTICLE_RADIUS;



            valid_position = 1;  // Assume OK first
            for (int j = 0; j < i; j++) { // Check against all already placed particles
                float dx = X[i] - X[j];
                float dy = Y[i] - Y[j];
                float dist_sq = dx * dx + dy * dy;
                if (dist_sq < (4.0f * PARTICLE_RADIUS * PARTICLE_RADIUS)) {
                    valid_position = 0;  // Overlaps another particle
                    break;
                }
            }

            attempts++;
            if (attempts > 10000) {
                printf("❌ Could not place particle %d safely after 10000 tries! Exiting.\n", i);
                exit(1);
            }
        }

        if (X[i] > wall_x - wall_buffer && X[i] < wall_x + wall_buffer) {
            printf("⚠️ Particle %d is too close to wall! x = %f\n", i, X[i]);
        }
        

        Vx[i] = maxwell_boltzmann_velocity_gaussians(TEMPERATURE);
        Vy[i] = maxwell_boltzmann_velocity_gaussians(TEMPERATURE);

        Radius[i] = PARTICLE_RADIUS;

        printf("✅ Particle %d initialized: X = %f, Y = %f, Vx = %f, Vy = %f\n", i, X[i], Y[i], Vx[i], Vy[i]);
    }
}





/// place particles on grid honeycomb, works
void initialize_simulation_honeycomb() {
    float wall_buffer = WALL_THICKNESS + PARTICLE_RADIUS * 0.20f;  // space around wall
    float small_gap = PARTICLE_RADIUS * 0.01f;  // just a tiny gap

    float left_width = wall_x - wall_buffer - XW1;
    float right_width = XW2 - wall_x - wall_buffer;
    float height = YW2 - YW1;

    int particles_left = NUM_PARTICLES / 2;
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





void initialize_simulation() {
    float wall_buffer = WALL_THICKNESS + DIAMETER * 0.25f; // space from wall
    float margin = DIAMETER * 1.25f;

    float left_width = wall_x - wall_buffer - XW1;
    float right_width = XW2 - wall_x - wall_buffer;
    float height = YW2 - YW1;

    int particles_left = NUM_PARTICLES / 2;
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
    float left_box_edge = 0;
    float piston_width = 20;
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
    float right_box_edge = MAX_X - 20; // Adjusted for piston width
    float left_limit = piston_left_x + 20; // Ensure pistons don't overlap

    piston_right_x += vx_piston_right * dt;

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

    // Left piston (infinite mass moving at vx_piston_left)
    if (X[i] - Radius[i] <= piston_left_x) {
        float vrel = Vx[i] - vx_piston_left;
        if (vrel < 0.f) {
            // reflect in piston frame: v' = 2*Vp - v
            Vx[i] = 2.f * vx_piston_left - Vx[i];
        }
        // place just outside, do NOT advance by dt here
        X[i] = piston_left_x + Radius[i] + eps;
    }

    // Right piston
    if (X[i] + Radius[i] >= piston_right_x) {
        float vrel = Vx[i] - vx_piston_right;
        if (vrel > 0.f) {
            Vx[i] = 2.f * vx_piston_right - Vx[i];
        }
        X[i] = piston_right_x - Radius[i] - eps;
    }
}



//////// WALL FUNCTIONS /////////

void draw_wall() {
    if (!wall_enabled) return;  // Don't draw if the wall is disabled
       if (isnan(wall_x) || wall_x < 0 || wall_x > SIM_WIDTH) {
        printf("🚨 Invalid wall_x = %f\n", wall_x);
            exit(1);
        }
    SDL_SetRenderDrawColor(_renderer, 255, 255, 255, 255);  // White color
    SDL_Rect wall_rect = { (wall_x - (WALL_THICKNESS/2)), 0, WALL_THICKNESS, SIM_HEIGHT };  // Thin vertical wall

 
    SDL_RenderFillRect(_renderer, &wall_rect);
}


// Function to update wall position and velocity
void update_wall(float dt, float L0) {
    //steps_elapsed++;  // Increment the step counter

    // Inside update_wall()
    if (wall_hold_enabled && !wall_is_released) {
        vx_wall = 0.0f;

        if (!wall_position_is_managed_externally) {
            wall_x = XW1 + (L0 * PIXELS_PER_SIGMA);  // Only reset if not overridden
        }

        if (steps_elapsed >= wall_hold_steps) {
            wall_is_released = true;
            wall_release_time = steps_elapsed * dt;
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
        vx_wall = 0.0f;  // Stop wall if hitting left boundary
    }
    if (wall_x > XW2 - 0.1 * L0) {
        printf("🛑 Wall hit right boundary, clamping and resetting vx_wall\n");

        wall_x = XW2 - 0.1 * L0;
        vx_wall = 0.0f;  // Stop wall if hitting right boundary
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

static void ccd_wall_step(int i, float dt,
                          const float wall_x0,     // frozen pose (center) at substep start
                          const float wall_vx0,    // frozen velocity for this substep
                          const int   prev_side,   // -1 (left), 0 (inside), +1 (right)
                          float *wall_impulse_x_accum)
{
    if (!wall_enabled || dt <= 0.f) return;

    // Particle state
    float xs  = X[i],  ys  = Y[i];
    float vxs = Vx[i], vys = Vy[i];
    const float R = Radius[i];

    // Wall faces at t=0 (pose frozen this substep)
    const float wL0 = wall_x0 - 0.5f * WALL_THICKNESS;
    const float wR0 = wall_x0 + 0.5f * WALL_THICKNESS;
    const float vwx = wall_vx0;

    // Inflate by particle radius
    const float L0p = wL0 - R;
    const float R0p = wR0 + R;

    // Tolerances & hysteresis
    const float eps_pos = fmaxf(1e-6f, 1e-4f * fmaxf(R, 1.f));
    const float eps_t   = 1e-12f;
    const float v_hys   = 1e-6f; // co-moving band

    float t_remain = dt;

    for (int iter = 0; iter < 6 && t_remain > eps_t; ++iter) {
        const float vrel = vxs - vwx;
        const bool inside = (xs > L0p - eps_pos) && (xs < R0p + eps_pos);

        // A) Co-moving band: no resolution; eject if inside and free fly
        if (fabsf(vrel) < v_hys) {
            if (inside) {
                int face = (prev_side != 0) ? prev_side
                          : (fabsf(xs - L0p) <= fabsf(R0p - xs) ? -1 : +1);
                xs = (face == -1) ? (L0p - eps_pos) : (R0p + eps_pos);
            }
            xs += vxs * t_remain;
            ys += vys * t_remain;
            t_remain = 0.f;
            break;
        }

        // B) Start-inside (not co-moving): eject and bounce once (finite mass wall in x)
        if (inside) {
            int face = (prev_side != 0) ? prev_side
                      : (fabsf(xs - L0p) <= fabsf(R0p - xs) ? -1 : +1);

            const float m = (float)PARTICLE_MASS;
            const float M = (float)WALL_MASS;
            const float u1 = vxs;  // pre-impact particle x-vel
            const float u2 = vwx;  // frozen wall x-vel

            const float v1 = ((m - M)*u1 + 2.f*M*u2) / (m + M);
            vxs = v1;

            // impulse on particle Δp = m*(v1 - u1), equal/opposite to wall
            if (wall_impulse_x_accum) {
                float Jp = m * (v1 - u1);
                *wall_impulse_x_accum -= Jp;
            }

            xs = (face == -1) ? (L0p - eps_pos) : (R0p + eps_pos);
            xs += vxs * (1e-9f);
            continue;
        }

        // C) Start-outside: compute earliest valid TOI with either face
        float thit = INFINITY; int face = 0;

        if (vrel > 0.f) {
            float tL = (L0p - xs) / vrel; // to left face
            if (tL >= -eps_t && tL <= t_remain + eps_t) { thit = tL; face = -1; }
        }
        if (vrel < 0.f) {
            float tR = (R0p - xs) / vrel; // to right face (vrel < 0)
            if (tR >= -eps_t && tR <= t_remain + eps_t && tR < thit) { thit = tR; face = +1; }
        }

        if (face == 0 || !isfinite(thit) || thit > t_remain - eps_t) {
            // no impact within remainder: free flight
            xs += vxs * t_remain;
            ys += vys * t_remain;
            t_remain = 0.f;
            break;
        }

        // Advance to impact
        xs += vxs * thit;
        ys += vys * thit;

        // Elastic exchange along x against finite-mass wall (kinematic this substep)
        {
            const float m = (float)PARTICLE_MASS;
            const float M = (float)WALL_MASS;
            const float u1 = vxs, u2 = vwx;
            const float v1 = ((m - M)*u1 + 2.f*M*u2) / (m + M);
            if (wall_impulse_x_accum) {
                float Jp = m * (v1 - u1);
                *wall_impulse_x_accum -= Jp;  // opposite on wall
            }
            vxs = v1;
        }

        // Place just outside the impacting face + tiny nudge
        xs = (face == -1) ? (L0p - eps_pos) : (R0p + eps_pos);
        xs += vxs * (1e-9f);

        t_remain -= thit;
    }

    // Commit
    X[i]  = xs;  Y[i]  = ys;
    Vx[i] = vxs; Vy[i] = vys;
}






void update_particles_with_substepping(float dt, int* out_left, int* out_right) {
    if (dt <= 0.0f) return;
    
    #if SUBSTEPPING
        float sub_dt = dt / SUBSTEPS;
    #endif

    int Side_old[NUM_PARTICLES];
    int Side_new[NUM_PARTICLES];
    left_count = right_count = 0;

    for (int step = 0; step < (int) SUBSTEPS; ++step) {
        for (int i = 0; i < NUM_PARTICLES; i++) {
            X_old[i] = X[i];

            // classify side (optional)
            if (wall_enabled) {
                if (X_old[i] < wall_x - WALL_THICKNESS / 2)      Side_old[i] = -1;
                else if (X_old[i] > wall_x + WALL_THICKNESS / 2) Side_old[i] =  1;
                else                                             Side_old[i] =  0;
            }

            const float w_x0_sub = wall_x_old;   // frozen pose for this fixed step
            const float w_vx_sub = vx_wall;      // frozen vel for this substep
            int prev_side = Side_old[i];         // -1 left, 0 inside, +1 right

            ccd_wall_step(i, sub_dt, w_x0_sub, w_vx_sub, prev_side, &wall_impulse_x_accum);

            // pistons & boundaries
            handle_piston_collisions(i, sub_dt);
            handle_boundary_collision(i);

           // simple pair collisions (equal mass, perfectly elastic)
            for (int j = i + 1; j < NUM_PARTICLES; j++) {
                float dx = X[i] - X[j];
                float dy = Y[i] - Y[j];
                float dist2 = dx*dx + dy*dy;
                float min_d = 2.f * PARTICLE_RADIUS;

                if (dist2 < min_d * min_d && dist2 > 0.f) {
                    float dist = sqrtf(dist2);
                    float nx = dx / dist;
                    float ny = dy / dist;

                    // relative vel along normal
                    float relx = Vx[i] - Vx[j];
                    float rely = Vy[i] - Vy[j];
                    float vn = relx * nx + rely * ny;
                    if (vn < 0.f) {
                        // equal masses, e=1 → subtract/add vn along n
                        Vx[i] -= vn * nx;
                        Vy[i] -= vn * ny;
                        Vx[j] += vn * nx;
                        Vy[j] += vn * ny;

                        // positional correction (split overlap)
                        float overlap = min_d - dist;
                        float push = 0.5f * overlap + 1e-6f; // tiny bias to avoid recontact
                        X[i] += nx * push;  Y[i] += ny * push;
                        X[j] -= nx * push;  Y[j] -= ny * push;
                    }
                }
            }
        }
    }

    if (wall_enabled) {
        left_count = right_count = 0;
        for (int i = 0; i < NUM_PARTICLES; i++) {
            if (X[i] < wall_x - WALL_THICKNESS / 2) left_count++;
            else if (X[i] > wall_x + WALL_THICKNESS / 2) right_count++;
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
                printf("🔔 Wall manually released by key 'r'.\n");}
                break;

            case SDLK_p:  // Pause/unpause simulation
            paused = !paused;
            printf("Simulation %s\n", paused ? "Paused" : "Running");
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
            if (piston_left_x > 0) {
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
            break;
        case SDLK_RIGHT:
            if (piston_right_x < MAX_X) {
                vx_piston_right += acceleration;
                if (vx_piston_right > max_velocity) vx_piston_right = max_velocity;
            }
            break;
        case SDLK_UP: // Stop right piston immediately
            vx_piston_right = 0;
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








/// SPEED OF SOUND EXPERIMENTS

void run_speed_of_sound_experiments() {
    system("mkdir -p experiments_speed_of_sound/mode0_real_units");
    system("mkdir -p experiments_speed_of_sound/mode1_normalized_units");
    system("rm -f experiments_speed_of_sound/mode0_real_units/*.csv");
    system("rm -f experiments_speed_of_sound/mode1_normalized_units/*.csv");

    wall_enabled = 1;                        // ensure wall exists for CCD
    wall_position_is_managed_externally = 1; // you control center via L0

    float lengths[] = {7.5f, 10.0f, 15.0f};
    int wall_mass_factors[] = {20, 50, 100, 200, 500, 1000};
    int num_repeats = 1;

    printf("Running each experiment for %d steps after wall release\n", num_steps);
    if (num_steps > 10000000) {
        printf("⚠️ WARNING: Simulation time is very long! Consider reducing the number of steps.\n");
    } else if (num_steps < 1000) {
        printf("⚠️ WARNING: Simulation time is very short, SET EMERGENCY STEPS!\n");
        num_steps = (int)(200000 / FIXED_DT);
    }

    for (int l = 0; l < (int)(sizeof(lengths)/sizeof(lengths[0])); ++l) {
        for (int m = 0; m < (int)(sizeof(wall_mass_factors)/sizeof(wall_mass_factors[0])); ++m) {
            for (int r = 0; r < num_repeats; ++r) {

                float L0 = lengths[l];
                wall_mass_runtime = PARTICLE_MASS * wall_mass_factors[m];
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
                        mode_folder, L0_int, wall_mass_factors[m], r);
                FILE *wall_log = fopen(filename, "w");
                if (!wall_log) { printf("❌ Could not open log file.\n"); continue; }

                if (log_packing_fraction) {
                    log_packing_fractions(wall_log, "main");
                } else {
                    fprintf(wall_log, "Time, Wall_X, Displacement(σ), Left_Count, Right_Count\n");
                }

                printf("🔬 Running: L0 = %.1f, M = %d*m, run = %d\n", L0, wall_mass_factors[m], r);
                printf("🔍 Initial wall_x = %.3f, vx_wall = %.6f\n", wall_x, vx_wall);

                int left_particles = 0, right_particles = 0;
                int recorded_steps = 0;

                while (recorded_steps < num_steps) {
                    // --- START OF STEP ---
                    wall_x_old = wall_x;   // ← CCD uses wall pose from start of step
                    steps_elapsed++;       // ← tick BEFORE update_wall so hold/release can trigger

                    // physics
                    wall_impulse_x_accum = 0.f; // (re)start accumulator for this step
                    update_particles_with_substepping(FIXED_DT, &left_particles, &right_particles);

                    // 🔹 integrate wall velocity from accumulated impulses
                    if (wall_is_released) {
                        vx_wall += wall_impulse_x_accum / WALL_MASS;   // Δv = J / M
                    }
                    // (optional) clear accumulator here if you prefer step-scoped use
                    // wall_impulse_x_accum = 0.f;

                    update_pistons(FIXED_DT);
                    update_wall(FIXED_DT, L0_UNITS);

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

                // convert accumulated impulse to wall Δv once per step (after all substeps)
                if (wall_is_released) {
                    vx_wall += wall_impulse_x_accum / WALL_MASS;   // ← crucial
                }

                // optional: zero it here if you prefer step-scoped
                // wall_impulse_x_accum = 0.f;

                update_pistons(FIXED_DT);
                update_wall(FIXED_DT, L0_UNITS);

                simulation_time += FIXED_DT;

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
        draw_text(_renderer, font, temp_label, XW1 + 5, YW1 + HEIGHT_UNITS*PIXELS_PER_SIGMA, yellow);

        SDL_RenderPresent(_renderer);
        SDL_Delay(5);
    }

    fclose(logFile);
    fclose(wall_log);
}










/////////////// MAIN FUNCTION /////////////
// Main Function
int main(int argc, char* argv[]) {

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


