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




//////////// SIMULATION SETTINGS ////////////
// Global simulation flags

bool enable_speed_of_sound_experiments = false;  
#define MOLECULAR_MODE 0       // 1 = normalized units; 0 = real units
#define EXPERIMENT_MODE 0       // 1 = tuned dt/substep values for better experiments

// PARTICLE PARAMETERS
#define NUM_PARTICLES 100      // total number of particles
#define PARTICLE_RADIUS 1      // normalized radius = 1 (works in both modes!)
#define DIAMETER (2 * PARTICLE_RADIUS)
#define METERS_PER_PIXEL 1e-9f // pixel size for mapping to meters (important for real mode)

// --- Physical constants ---
#if MOLECULAR_MODE
    #define PARTICLE_MASS 1.0f
    #define K_B 1.0f
    #define TEMPERATURE 1.0f

#else
    #define PARTICLE_MASS 1.67e-27f      // kg (proton mass)
    #define K_B 1.38e-23f                // J/K
    #define TEMPERATURE 300.0f           // Kelvin
#endif

// --- Time step and substeps ---
#if MOLECULAR_MODE
    #if EXPERIMENT_MODE
        #define FIXED_DT 1e-14f
        #define SUBSTEPS 10
        #define WALL_HOLD_STEPS 10000

    #else
        #define FIXED_DT 1e-14f
        #define SUBSTEPS 5
        #define WALL_HOLD_STEPS 5000

    #endif
#else
    #if EXPERIMENT_MODE
        #define FIXED_DT 2e-4f
        #define SUBSTEPS 20
        #define WALL_HOLD_STEPS 10000

    #else
        #define FIXED_DT 1.6e-3f
        #define SUBSTEPS 10
        #define WALL_HOLD_STEPS 5000

    #endif
#endif




// WALL PARAMETERS
#define WALL_MASS_FACTOR 20                               // choose factor relative to PARTICLE_MASS

float wall_mass_runtime = PARTICLE_MASS * WALL_MASS_FACTOR;
#define WALL_MASS (wall_mass_runtime)
#define WALL_THICKNESS (PARTICLE_RADIUS)                  // simple for now
bool wall_hold_enabled = true;         // Enable or disable hold wall mode
bool wall_is_released = false;          // Has the user released the wall?
int wall_hold_steps = WALL_HOLD_STEPS;             // Hold the wall for this many steps
int steps_elapsed = 0;    


// --- Derived Time Scale τ ---
float TAU_time_scale_factor_to_molecular = 1.0f;           // default (molecular mode)
void initialize_time_scale() {
#if !MOLECULAR_MODE
    TAU_time_scale_factor_to_molecular = sqrtf((PARTICLE_MASS * powf(PARTICLE_RADIUS * METERS_PER_PIXEL, 2)) / (K_B * TEMPERATURE));
#endif
}



/// SIMULATION WINDOW PARAMETERS

// Set these in main or before simulation starts
float L0_UNITS = 25.0f;             // Box half-length in diameters (paper uses 7.5, 10, ...)
float HEIGHT_UNITS = 10.0f;
float PIXELS_PER_DIAMETER = 20.0f;

int SIM_WIDTH, SIM_HEIGHT;
int XW1 = 50;
int XW2, YW1 = 0, YW2;
int MAX_X, MAX_Y;
int HIST_WIDTH, HIST_HEIGHT;

void initialize_simulation_dimensions() {
    SIM_WIDTH = (int)(2 * L0_UNITS * PIXELS_PER_DIAMETER);
    SIM_HEIGHT = (int)(HEIGHT_UNITS * PIXELS_PER_DIAMETER);

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
    wall_x = (XW1 + XW2) / 2.0f;

    wall_x_old = wall_x;
    vx_wall = 0.0f;

    piston_left_x = 0.0f;
    vx_piston_left = 0.0f;
    vx_piston_right = 0.0f;
}



 // VARIABLES TO AVOID TUNNELING die to time step and stationary particles
float X_old[NUM_PARTICLES];  // Previous x-position
int left_count = 0;
int right_count = 0;



// SIMULATION PARAMETERS
#define MIN_SPATIAL_WINDOW (MAX_X/10 * MAX_Y/10)/ NUM_PARTICLES// Prevent too small window
#define NUM_BINS 100 // NUM bins for histogram
#define MAX_VELOCITY 10000.0
#define SEED 42  // Seed for random number generation (ensures repeatable randomness)
#define SCALE_Vxy  1 // Mass of particle (e.g., hydrogen atom) in kg
#define Scale_Energy 1e27/1000






// Particle states
float X[NUM_PARTICLES];
float Y[NUM_PARTICLES];
float X_old[NUM_PARTICLES];    // previous X positions

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



//////// INITIALIZATION FUNCTIONS /////////
// Function to generate Maxwell-Boltzmann distributed velocity 1D (v_x is normal distributed!!!!!)
// vx or vy 
// f(v_x) = \sqrt{\frac{m}{2\pi k_B T}}  \exp\left(-\frac{m v_x^2}{2 k_B T}\right)
//  •	Mean:  \mu = 0  (because velocity components are symmetrically distributed around zero)
//	•	Standard deviation: \sigma = \sqrt{\frac{k_B T}{m}}
// ---> f(v_x) = \frac{1}{\sqrt{2\pi \sigma^2}} \exp\left(-\frac{v_x^2}{2\sigma^2}\right)
//Since  v_x  follows a normal distribution, we need to sample from a Gaussian (normal) distribution with mean  0  and standard deviation  \sigma .
//The Box-Muller transform allows to generate 1 normally distributed numbers from two uniform random numbers  u_1  and  u_2 , which are sampled from the uniform distribution  U(0,1) .
// complete speed v
// f(v) \propto v^2 \exp\left(-\frac{mv^2}{2k_BT}\right)
float maxwell_boltzmann_velocity(float temperature) {
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




void initialize_simulation() {
    int left_particles = NUM_PARTICLES / 2;
    int right_particles = NUM_PARTICLES - left_particles;

    float wall_buffer = WALL_THICKNESS + PARTICLE_RADIUS*2;

    for (int i = 0; i < NUM_PARTICLES; i++) {
        int valid_position = 0;
        // Randomly assign particles to left or right side of the wall
        while (!valid_position) {
            if (i < left_particles) {
                // Prevent left particles from spawning too close to the wall
                X[i] = ((float)rand() / RAND_MAX) * (wall_x - wall_buffer - XW1) + XW1 + PARTICLE_RADIUS;
            } else {
                // Prevent right particles from spawning too close to the wall
                X[i] = ((float)rand() / RAND_MAX) * (XW2 - wall_x - wall_buffer) + wall_x + wall_buffer;
            }

            Y[i] = ((float)rand() / RAND_MAX) * (YW2 - YW1 - 2 * PARTICLE_RADIUS) + YW1 + PARTICLE_RADIUS;

            valid_position = 1;  // You could later add overlap checks here too
        }
        if (X[i] > wall_x - wall_buffer && X[i] < wall_x + wall_buffer) {
        printf("⚠️ Particle %d is too close to wall! x = %f\n", i, X[i]);
        }

        Vx[i] = maxwell_boltzmann_velocity(TEMPERATURE);
        Vy[i] = maxwell_boltzmann_velocity(TEMPERATURE);

        Radius[i] = PARTICLE_RADIUS;

        printf("Particle %d initialized: X = %f, Y = %f, Vx = %f, Vy = %f\n", i, X[i], Y[i], Vx[i], Vy[i]);
    }
}












/////////////////// PISTON FUNCTIONS /////////////////////////

// Function to render pistons
void render_pistons() {
    SDL_SetRenderDrawColor(_renderer, 0, 255, 0, 0);
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














//////////////////// PARTICLE FUNCTIONS --  ///////////////////////////
// Function to render particles using drawpoint (optimized) 
void SDL_RenderFillCircle(SDL_Renderer* renderer, int centerX, int centerY, int radius) {
    for (int w = 0; w < radius * 2; w++) {
        for (int h = 0; h < radius * 2; h++) {
            int dx = radius - w;
            int dy = radius - h;
            if ((dx * dx + dy * dy) <= (radius * radius)) {
                SDL_RenderDrawPoint(renderer, centerX + dx, centerY + dy);
            }
        }
    }
}


void render_particles() {
        SDL_SetRenderDrawColor(_renderer, 255, 255, 255, 255);
        for (int i = 0; i < NUM_PARTICLES; i++) {
            SDL_RenderFillCircle(_renderer, (int)X[i], (int)Y[i], PARTICLE_RADIUS);
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





void handle_piston_collisions(int i) {
    // Check collision with the left piston
    if (X[i] - Radius[i] <= piston_left_x) {  
        float relative_velocity = Vx[i] - vx_piston_left ;
        
        if (relative_velocity < 0) {  // Only if the particle is moving toward the piston
            Vx[i] = -Vx[i] + 2 * vx_piston_left ;
        } else {
            Vx[i] = -Vx[i];  // Simple reflection if moving away
        }
        
        X[i] = piston_left_x + Radius[i] + Vx[i] * FIXED_DT;  // Adjust for timestep
    }

    // Check collision with the right piston
    if (X[i] + Radius[i] >= piston_right_x) {  
        float relative_velocity = Vx[i] - vx_piston_right ;

        if (relative_velocity > 0) {  // Only if the piston moves into the particle
            Vx[i] = -Vx[i] + 2 * vx_piston_right ;
        } else {
            Vx[i] = -Vx[i];  // Simple reflection if moving away
        }

        X[i] = piston_right_x - Radius[i] + Vx[i] * FIXED_DT;  // Adjust for timestep
    }
}




//////// WALL FUNCTIONS /////////




void draw_wall() {
    if (!wall_enabled) return;  // Don't draw if the wall is disabled
    SDL_SetRenderDrawColor(_renderer, 255, 255, 255, 255);  // White color
    SDL_Rect wall_rect = { (wall_x - (WALL_THICKNESS/2)), 0, WALL_THICKNESS, SIM_HEIGHT };  // Thin vertical wall
    SDL_RenderFillRect(_renderer, &wall_rect);
}


// Function to update wall position and velocity
void update_wall(float dt) {
    steps_elapsed++;

    // === Wall Hold Phase ===
    if (wall_hold_enabled && !wall_is_released) {
        vx_wall = 0.0f;
        wall_x = (XW1 + XW2) / 2.0f;  // Center wall between boundaries

        // Auto-release after a certain number of steps
        if (steps_elapsed >= wall_hold_steps) {
            wall_is_released = true;
            printf("🔔 Wall released automatically after %d steps.\n", steps_elapsed);
        }
        return;  // Exit early, skip wall motion update
    }

    // === Normal Wall Motion Phase ===
    wall_x_old = wall_x;  // Store previous wall position
    wall_x += vx_wall * dt;  // Move wall based on current velocity

    // === Enforce Wall Boundaries ===
    if (wall_x < XW1 + 10) {
        wall_x = XW1 + 10;
        vx_wall = 0.0f;  // Stop wall if hitting left boundary
    }
    if (wall_x > XW2 - 10) {
        wall_x = XW2 - 10;
        vx_wall = 0.0f;  // Stop wall if hitting right boundary
    }

    printf("wall_x = %f, wall_enabled = %d\n", wall_x, wall_enabled);
}




void handle_left_to_right_collision(int i) {
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
    float v_particle = Vx[i];
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



bool swept_particle_collision_index(int i, int j, float dt, float *collision_time) {
    float dx = X[j] - X[i];
    float dy = Y[j] - Y[i];
    float dvx = Vx[j] - Vx[i];
    float dvy = Vy[j] - Vy[i];

    float radius_sum = 2 * PARTICLE_RADIUS;
    float r2 = radius_sum * radius_sum;

    float A = dvx * dvx + dvy * dvy;
    float B = 2 * (dx * dvx + dy * dvy);
    float C = dx * dx + dy * dy - r2;

    float discriminant = B * B - 4 * A * C;

    if (A == 0.0f || discriminant < 0.0f) return false;

    float sqrt_discriminant = sqrtf(discriminant);
    float t1 = (-B - sqrt_discriminant) / (2 * A);
    float t2 = (-B + sqrt_discriminant) / (2 * A);

    if (t1 >= 0.0f && t1 <= dt) {
        *collision_time = t1;
        return true;
    } else if (t2 >= 0.0f && t2 <= dt) {
        *collision_time = t2;
        return true;
    }

    return false;
}

void resolve_particle_collision_index(int i, int j) {
    float dx = X[j] - X[i];
    float dy = Y[j] - Y[i];
    float dist2 = dx * dx + dy * dy;
    if (dist2 == 0.0f) return;

    float dvx = Vx[j] - Vx[i];
    float dvy = Vy[j] - Vy[i];

    float dot = (dvx * dx + dvy * dy) / dist2;

    float fx = dot * dx;
    float fy = dot * dy;

    Vx[i] += fx;
    Vy[i] += fy;
    Vx[j] -= fx;
    Vy[j] -= fy;
}




void update_particles_with_collisions_indexed(float dt) {
    float sub_dt = dt / SUBSTEPS;

    for (int step = 0; step < SUBSTEPS; step++) {
        // Move all particles
        for (int i = 0; i < NUM_PARTICLES; i++) {
            X[i] += Vx[i] * sub_dt;
            Y[i] += Vy[i] * sub_dt;
        }

        // Handle particle-particle collisions
        for (int i = 0; i < NUM_PARTICLES; i++) {
            for (int j = i + 1; j < NUM_PARTICLES; j++) {
                float collision_time;

                if (swept_particle_collision_index(i, j, sub_dt, &collision_time)) {
                    // Rewind
                    float t = collision_time;
                    X[i] -= Vx[i] * (sub_dt - t);
                    Y[i] -= Vy[i] * (sub_dt - t);
                    X[j] -= Vx[j] * (sub_dt - t);
                    Y[j] -= Vy[j] * (sub_dt - t);

                    // Resolve
                    resolve_particle_collision_index(i, j);

                    // Move forward
                    X[i] += Vx[i] * (sub_dt - t);
                    Y[i] += Vy[i] * (sub_dt - t);
                    X[j] += Vx[j] * (sub_dt - t);
                    Y[j] += Vy[j] * (sub_dt - t);
                }
            }
        }
    }
}








void detect_and_fix_tunnel(float *X, float *X_old, float *Vx, int *Side_old, int *Side_new, int num_particles, float wall_x, float wall_thickness) {

/*
Tunneling occurs when a particle moves so far in one timestep that it completely skips over the wall or another object. To avoid it:

\Delta x_{\text{particle}} = v \cdot \Delta t < \text{wall thickness}

So if you want to reduce tunneling, make sure:
\Delta t < \frac{\text{wall thickness}}{\text{max particle speed}}
*/

    float wall_left = wall_x - wall_thickness / 2;
    float wall_right = wall_x + wall_thickness / 2;
    for (int i = 0; i < num_particles; i++) {
        if (Side_old[i] == -1 && Side_new[i] == 1) {
            printf("\u26d4\ufe0f Particle %d tunneled LEFT \u279e RIGHT. Reverting position.\n", i);
            X[i] = X_old[i];
            Vx[i] *= -1;
        } else if (Side_old[i] == 1 && Side_new[i] == -1) {
            printf("\u26d4\ufe0f Particle %d tunneled RIGHT \u279e LEFT. Reverting position.\n", i);
            X[i] = X_old[i];
            Vx[i] *= -1;
        }
    }
}











void update_particles_substep(float dt) {
    int Side_old[NUM_PARTICLES];
    int Side_new[NUM_PARTICLES];
    left_count = 0;
    right_count = 0;

    for (int i = 0; i < NUM_PARTICLES; i++) {
        X_old[i] = X[i];

        float wall_left = wall_x - WALL_THICKNESS / 2;
        float wall_right = wall_x + WALL_THICKNESS / 2;

        if (wall_enabled) {
            Side_old[i] = (X_old[i] < wall_left) ? -1 : (X_old[i] > wall_right) ? 1 : 0;
        }

        X[i] += Vx[i] * dt;
        Y[i] += Vy[i] * dt;

        if (wall_enabled) {
            Side_new[i] = (X[i] < wall_left) ? -1 : (X[i] > wall_right) ? 1 : 0;
        }

        if (wall_enabled) {
            if (X[i] < wall_left) left_count++;
            else if (X[i] > wall_right) right_count++;
        } else {
            right_count++;
        }

        bool was_left = X_old[i] + Radius[i] < wall_left;
        bool now_right = X[i] - Radius[i] > wall_right;
        bool was_right = X_old[i] - Radius[i] > wall_right;
        bool now_left = X[i] + Radius[i] < wall_left;

        if ((was_left && now_right) || (was_right && now_left)) {
            printf("\u26a0\ufe0f Particle %d tunneled in substep! Reverting position.\n", i);
            X[i] = X_old[i];
            Vx[i] *= -1;
            paused = true;
            return;
        }

        check_wall_collisions(i);
        handle_piston_collisions(i);
        handle_boundary_collision(i);

        for (int j = i + 1; j < NUM_PARTICLES; j++) {
            float dx = X[i] - X[j];
            float dy = Y[i] - Y[j];
            float dist = sqrtf(dx * dx + dy * dy);
            float min_dist = PARTICLE_RADIUS * 2;

            if (dist < min_dist) {
                float nx = dx / dist;
                float ny = dy / dist;
                float relVelX = Vx[i] - Vx[j];
                float relVelY = Vy[i] - Vy[j];
                float v_rel_dot = relVelX * nx + relVelY * ny;

                if (v_rel_dot < 0) {
                    float J = 2 * v_rel_dot / (PARTICLE_MASS + PARTICLE_MASS);

                    Vx[i] -= J * PARTICLE_MASS * nx;
                    Vy[i] -= J * PARTICLE_MASS * ny;
                    Vx[j] += J * PARTICLE_MASS * nx;
                    Vy[j] += J * PARTICLE_MASS * ny;

                    float overlap = min_dist - dist;
                    X[i] += nx * overlap / 2;
                    Y[i] += ny * overlap / 2;
                    X[j] -= nx * overlap / 2;
                    Y[j] -= ny * overlap / 2;
                }
            }
        }
    }

    if (wall_enabled && (left_count != right_count)) {
        detect_and_fix_tunnel(X, X_old, Vx, Side_old, Side_new, NUM_PARTICLES, wall_x, WALL_THICKNESS);
        paused = true;
    }
}






void update_particles(float dt) {
    int substeps = 4;
    float sub_dt = dt / substeps;
    for (int s = 0; s < substeps; s++) {
        if (!paused) {
            update_particles_substep(sub_dt);
        } else {
            break;
        }
    }
}




/// SECOND TRY no sub or time stepping





// Function to update particles with a fixed time step no substep
void update_particles2(float dt) {
    left_count = 0;
    right_count = 0;
    int Side_old[NUM_PARTICLES];
    int Side_new[NUM_PARTICLES];

    for (int i = 0; i < NUM_PARTICLES; i++) {
        X_old[i] = X[i];

        // Count before position update
        if (wall_enabled) {
            if (X[i] < wall_x - WALL_THICKNESS / 2) {
                left_count++;
            } else if (X[i] > wall_x + WALL_THICKNESS / 2) {
                right_count++;
            }
        } else {
            right_count++;  // All to the right if no wall
        }

        // Track side before update
        if (wall_enabled) {
            if (X_old[i] < wall_x - WALL_THICKNESS / 2)
                Side_old[i] = -1;
            else if (X_old[i] > wall_x + WALL_THICKNESS / 2)
                Side_old[i] = 1;
            else
                Side_old[i] = 0;
        }

        // Update positions
        X[i] += Vx[i] * dt;
        Y[i] += Vy[i] * dt;

        // Track side after update
        if (wall_enabled) {
            if (X[i] < wall_x - WALL_THICKNESS / 2)
                Side_new[i] = -1;
            else if (X[i] > wall_x + WALL_THICKNESS / 2)
                Side_new[i] = 1;
            else
                Side_new[i] = 0;

            // Tunnel detection based on position and radius
            float wall_left = wall_x - WALL_THICKNESS / 2;
            float wall_right = wall_x + WALL_THICKNESS / 2;

            bool was_left = X_old[i] + Radius[i] < wall_left;
            bool now_right = X[i] - Radius[i] > wall_right;
            bool was_right = X_old[i] - Radius[i] > wall_right;
            bool now_left = X[i] + Radius[i] < wall_left;

            if ((was_left && now_right) || (was_right && now_left)) {
                printf("⚠️ Particle %d tunneled through the wall! Reverting position.\n", i);
                X[i] = X_old[i];
                Vx[i] *= -1;
                //paused = true;
                return;
            }
        }

        // Boundary and wall handling
        check_wall_collisions(i);
        handle_piston_collisions(i);
        handle_boundary_collision(i);

        // Inter-particle collisions
        for (int j = i + 1; j < NUM_PARTICLES; j++) {
            float dx = X[i] - X[j];
            float dy = Y[i] - Y[j];
            float dist = sqrt(dx * dx + dy * dy);
            float min_dist = PARTICLE_RADIUS * 2;

            if (dist < min_dist) {
                float nx = dx / dist;
                float ny = dy / dist;
                float relVelX = Vx[i] - Vx[j];
                float relVelY = Vy[i] - Vy[j];
                float v_rel_dotProduct = (relVelX * nx + relVelY * ny);

                if (v_rel_dotProduct < 0) {
                    float J_impulse = 2 * v_rel_dotProduct / (PARTICLE_MASS + PARTICLE_MASS);

                    Vx[i] -= J_impulse * PARTICLE_MASS * nx;
                    Vy[i] -= J_impulse * PARTICLE_MASS * ny;
                    Vx[j] += J_impulse * PARTICLE_MASS * nx;
                    Vy[j] += J_impulse * PARTICLE_MASS * ny;

                    float overlap = min_dist - dist;
                    X[i] += nx * overlap / 2;
                    Y[i] += ny * overlap / 2;
                    X[j] -= nx * overlap / 2;
                    Y[j] -= ny * overlap / 2;
                }
            }
        }
    }

    // Redundant tunnel fix — optional if inline check already pauses
    if (wall_enabled && (left_count != right_count)) {
        detect_and_fix_tunnel(X, X_old, Vx, Side_old, Side_new, NUM_PARTICLES, wall_x, WALL_THICKNESS);
        //paused = true;
    }
}

// third with substepping 
void update_particles_with_substepping1(float dt) {
    // Define substeps
    float sub_dt = dt / SUBSTEPS; // Time step per substep

    // Variables for tracking state before and after movement
    int Side_old[NUM_PARTICLES];
    int Side_new[NUM_PARTICLES];

    // Track how many particles are on each side of the wall (reset only once per timestep)
    left_count = 0;
    right_count = 0;

    // Count particles' initial positions before any substeps
    if (wall_enabled) {
        for (int i = 0; i < NUM_PARTICLES; i++) {
            if (X[i] < wall_x - WALL_THICKNESS / 2) {
                left_count++;
            } else if (X[i] > wall_x + WALL_THICKNESS / 2) {
                right_count++;
            }
        }
    } else {
        right_count = NUM_PARTICLES;  // All particles are to the right if no wall
    }

    for (int step = 0; step < SUBSTEPS; step++) {
        for (int i = 0; i < NUM_PARTICLES; i++) {
            X_old[i] = X[i];  // Save current positions

            // Track side before update
            if (wall_enabled) {
                if (X_old[i] < wall_x - WALL_THICKNESS / 2)
                    Side_old[i] = -1;
                else if (X_old[i] > wall_x + WALL_THICKNESS / 2)
                    Side_old[i] = 1;
                else
                    Side_old[i] = 0;
            }

            // Move particle in this substep
            X[i] += Vx[i] * sub_dt;
            Y[i] += Vy[i] * sub_dt;

            // Track side after update
            if (wall_enabled) {
                if (X[i] < wall_x - WALL_THICKNESS / 2)
                    Side_new[i] = -1;
                else if (X[i] > wall_x + WALL_THICKNESS / 2)
                    Side_new[i] = 1;
                else
                    Side_new[i] = 0;

                // Tunnel detection for wall
                float wall_left = wall_x - WALL_THICKNESS / 2;
                float wall_right = wall_x + WALL_THICKNESS / 2;

                bool was_left = X_old[i] + Radius[i] < wall_left;
                bool now_right = X[i] - Radius[i] > wall_right;
                bool was_right = X_old[i] - Radius[i] > wall_right;
                bool now_left = X[i] + Radius[i] < wall_left;

                if ((was_left && now_right) || (was_right && now_left)) {
                    printf("⚠️ Particle %d tunneled through the wall! Reverting position.\n", i);
                    X[i] = X_old[i];  // Revert position to last valid state
                    Vx[i] *= -1;  // Reverse velocity
                    return; // Optional: Stop execution for this timestep
                }
            }

            // Boundary and wall handling (collision with walls or boundaries)
            check_wall_collisions(i);
            handle_piston_collisions(i);
            handle_boundary_collision(i);

            // Particle-particle collisions
            for (int j = i + 1; j < NUM_PARTICLES; j++) {
                float dx = X[i] - X[j];
                float dy = Y[i] - Y[j];
                float dist = sqrt(dx * dx + dy * dy);
                float min_dist = PARTICLE_RADIUS * 2;

                if (dist < min_dist) {
                    // Normalized direction
                    float nx = dx / dist;
                    float ny = dy / dist;

                    // Relative velocity
                    float relVelX = Vx[i] - Vx[j];
                    float relVelY = Vy[i] - Vy[j];

                    // Dot product for collision response
                    float v_rel_dotProduct = (relVelX * nx + relVelY * ny);

                    if (v_rel_dotProduct < 0) {
                        float J_impulse = 2 * v_rel_dotProduct / (PARTICLE_MASS + PARTICLE_MASS);

                        // Apply impulse to velocities
                        Vx[i] -= J_impulse * PARTICLE_MASS * nx;
                        Vy[i] -= J_impulse * PARTICLE_MASS * ny;
                        Vx[j] += J_impulse * PARTICLE_MASS * nx;
                        Vy[j] += J_impulse * PARTICLE_MASS * ny;

                        // Apply overlap correction
                        float overlap = min_dist - dist;
                        X[i] += nx * overlap / 2;
                        Y[i] += ny * overlap / 2;
                        X[j] -= nx * overlap / 2;
                        Y[j] -= ny * overlap / 2;
                    }
                }
            }
        }
    }

    // Tunnel fix for edge cases (optional, as done at the end of substeps)
    if (wall_enabled && (left_count != right_count)) {
        detect_and_fix_tunnel(X, X_old, Vx, Side_old, Side_new, NUM_PARTICLES, wall_x, WALL_THICKNESS);
    }
}





//// 4th with substepping and variable time step

void update_particles_with_substepping2(float dt) {
    // Define substeps
    float sub_dt = dt / SUBSTEPS; // Time step per substep

    // Variables for tracking state before and after movement
    int Side_old[NUM_PARTICLES];
    int Side_new[NUM_PARTICLES];

    // Track how many particles are on each side of the wall
    left_count = 0;
    right_count = 0;

    if (wall_enabled) {
        // Count particles' initial positions before any substeps
        for (int i = 0; i < NUM_PARTICLES; i++) {
            if (X[i] < wall_x - WALL_THICKNESS / 2) {
                left_count++;
            } else if (X[i] > wall_x + WALL_THICKNESS / 2) {
                right_count++;
            }
        }
    } else {
        right_count = NUM_PARTICLES;  // All particles are to the right if no wall
    }

    // Loop over substeps
    for (int step = 0; step < SUBSTEPS; step++) {
        for (int i = 0; i < NUM_PARTICLES; i++) {
            X_old[i] = X[i];  // Save current positions

            // Track side before update
            if (wall_enabled) {
                if (X_old[i] < wall_x - WALL_THICKNESS / 2)
                    Side_old[i] = -1;
                else if (X_old[i] > wall_x + WALL_THICKNESS / 2)
                    Side_old[i] = 1;
                else
                    Side_old[i] = 0;
            }

            // Calculate adaptive time step based on velocity (simplified for X direction)
            float velocity = fabs(Vx[i]); // You can use the max velocity in the x or y direction
            float max_move = velocity * sub_dt; // Max distance the particle can move in one substep

            // Safety factor to avoid too large movements in one substep
            if (max_move > WALL_THICKNESS) {
                sub_dt = WALL_THICKNESS / velocity;  // Adjust time step to prevent tunneling
            }

            // Move particle in this substep
            X[i] += Vx[i] * sub_dt;
            Y[i] += Vy[i] * sub_dt;

            // Track side after update
            if (wall_enabled) {
                if (X[i] < wall_x - WALL_THICKNESS / 2)
                    Side_new[i] = -1;
                else if (X[i] > wall_x + WALL_THICKNESS / 2)
                    Side_new[i] = 1;
                else
                    Side_new[i] = 0;

                // Continuous collision detection (check if particle crosses the wall during the substep)
                float wall_left = wall_x - WALL_THICKNESS / 2;
                float wall_right = wall_x + WALL_THICKNESS / 2;

                bool was_left = X_old[i] + Radius[i] < wall_left;
                bool now_right = X[i] - Radius[i] > wall_right;
                bool was_right = X_old[i] - Radius[i] > wall_right;
                bool now_left = X[i] + Radius[i] < wall_left;

                if ((was_left && now_right) || (was_right && now_left)) {
                    // Revert position if tunneling occurs
                    printf("⚠️ Particle %d tunneled through the wall! Reverting position.\n", i);
                    X[i] = X_old[i];  // Revert position to last valid state
                    Vx[i] *= -1;  // Reverse velocity
                    return; // Optional: Stop execution for this timestep
                }
            }

            // Boundary and wall handling (collision with walls or boundaries)
            check_wall_collisions(i);
            handle_piston_collisions(i);
            handle_boundary_collision(i);

            // Particle-particle collisions
            for (int j = i + 1; j < NUM_PARTICLES; j++) {
                float dx = X[i] - X[j];
                float dy = Y[i] - Y[j];
                float dist = sqrt(dx * dx + dy * dy);
                float min_dist = PARTICLE_RADIUS * 2;

                if (dist < min_dist) {
                    // Normalized direction
                    float nx = dx / dist;
                    float ny = dy / dist;

                    // Relative velocity
                    float relVelX = Vx[i] - Vx[j];
                    float relVelY = Vy[i] - Vy[j];

                    // Dot product for collision response
                    float v_rel_dotProduct = (relVelX * nx + relVelY * ny);

                    if (v_rel_dotProduct < 0) {
                        float J_impulse = 2 * v_rel_dotProduct / (PARTICLE_MASS + PARTICLE_MASS);

                        // Apply impulse to velocities
                        Vx[i] -= J_impulse * PARTICLE_MASS * nx;
                        Vy[i] -= J_impulse * PARTICLE_MASS * ny;
                        Vx[j] += J_impulse * PARTICLE_MASS * nx;
                        Vy[j] += J_impulse * PARTICLE_MASS * ny;

                        // Apply overlap correction
                        float overlap = min_dist - dist;
                        X[i] += nx * overlap / 2;
                        Y[i] += ny * overlap / 2;
                        X[j] -= nx * overlap / 2;
                        Y[j] -= ny * overlap / 2;
                    }
                }
            }
        }
    }

    // Tunnel fix for edge cases (optional, as done at the end of substeps)
    if (wall_enabled && (left_count != right_count)) {
        detect_and_fix_tunnel(X, X_old, Vx, Side_old, Side_new, NUM_PARTICLES, wall_x, WALL_THICKNESS);
    }
}





//// 5th with substepping and variable time step including window reducing comp complexity
void update_particles_with_substepping(float dt,int* out_left, int* out_right) {
    // Define substeps
    float sub_dt = dt / SUBSTEPS; // Base time step per substep

    // Variables for tracking state before and after movement
    int Side_old[NUM_PARTICLES];
    int Side_new[NUM_PARTICLES];

    // Track how many particles are on each side of the wall
    left_count = 0;
    right_count = 0;

    if (wall_enabled) {
        // Count particles' initial positions before any substeps
        for (int i = 0; i < NUM_PARTICLES; i++) {
            if (X[i] < wall_x - WALL_THICKNESS / 2) {
                left_count++;
            } else if (X[i] > wall_x + WALL_THICKNESS / 2) {
                right_count++;
            }
        }
    } else {
        right_count = NUM_PARTICLES;  // All particles are to the right if no wall
    }

    // Loop over substeps
    for (int step = 0; step < SUBSTEPS; step++) {
        for (int i = 0; i < NUM_PARTICLES; i++) {
            X_old[i] = X[i];  // Save current positions

            // Track side before update
            if (wall_enabled) {
                if (X_old[i] < wall_x - WALL_THICKNESS / 2)
                    Side_old[i] = -1;
                else if (X_old[i] > wall_x + WALL_THICKNESS / 2)
                    Side_old[i] = 1;
                else
                    Side_old[i] = 0;
            }

            // --- VELOCITY-BASED SPATIAL WINDOW FOR COLLISIONS ---
            float velocity = sqrt(Vx[i] * Vx[i] + Vy[i] * Vy[i]);  // Magnitude of velocity
            float spatial_window = velocity * sub_dt;              // Effective influence range
            if (spatial_window < MIN_SPATIAL_WINDOW) spatial_window = MIN_SPATIAL_WINDOW;

            // --- ADAPTIVE TIMESTEP TO PREVENT TUNNELING ---
            float max_move = fabs(Vx[i]) * sub_dt;
            if (max_move > WALL_THICKNESS) {
                sub_dt = WALL_THICKNESS / fabs(Vx[i]);
            }

            // Move particle
            X[i] += Vx[i] * sub_dt;
            Y[i] += Vy[i] * sub_dt;

            // Track side after update
            if (wall_enabled) {
                if (X[i] < wall_x - WALL_THICKNESS / 2)
                    Side_new[i] = -1;
                else if (X[i] > wall_x + WALL_THICKNESS / 2)
                    Side_new[i] = 1;
                else
                    Side_new[i] = 0;

                // Tunnel detection logic
                float wall_left = wall_x - WALL_THICKNESS / 2;
                float wall_right = wall_x + WALL_THICKNESS / 2;

                bool was_left = X_old[i] + Radius[i] < wall_left;
                bool now_right = X[i] - Radius[i] > wall_right;
                bool was_right = X_old[i] - Radius[i] > wall_right;
                bool now_left = X[i] + Radius[i] < wall_left;

                if ((was_left && now_right) || (was_right && now_left)) {
                    printf("⚠️ Particle %d tunneled through the wall! Reverting position.\n", i);
                    X[i] = X_old[i];
                    Vx[i] *= -1;
                    return;
                }

                // Return the final counts
                if (out_left) *out_left = left_count;
                if (out_right) *out_right = right_count;
            }

            // Boundary + wall
            check_wall_collisions(i);
            handle_piston_collisions(i);
            handle_boundary_collision(i);

            // Particle-particle collisions with spatial window filtering
            for (int j = i + 1; j < NUM_PARTICLES; j++) {
                float dx = X[i] - X[j];
                float dy = Y[i] - Y[j];
                float dist_sq = dx * dx + dy * dy;

                // Skip if particles are outside each other's influence radius
                if (dist_sq > spatial_window * spatial_window) continue;

                float dist = sqrt(dist_sq);
                float min_dist = PARTICLE_RADIUS * 2;

                if (dist < min_dist) {
                    // Normalized direction
                    float nx = dx / dist;
                    float ny = dy / dist;

                    // Relative velocity
                    float relVelX = Vx[i] - Vx[j];
                    float relVelY = Vy[i] - Vy[j];
                    float v_rel_dot = relVelX * nx + relVelY * ny;

                    if (v_rel_dot < 0) {
                        float J_impulse = 2 * v_rel_dot / (PARTICLE_MASS + PARTICLE_MASS);

                        Vx[i] -= J_impulse * PARTICLE_MASS * nx;
                        Vy[i] -= J_impulse * PARTICLE_MASS * ny;
                        Vx[j] += J_impulse * PARTICLE_MASS * nx;
                        Vy[j] += J_impulse * PARTICLE_MASS * ny;

                        float overlap = min_dist - dist;
                        X[i] += nx * overlap / 2;
                        Y[i] += ny * overlap / 2;
                        X[j] -= nx * overlap / 2;
                        Y[j] -= ny * overlap / 2;
                    }
                }
            }
        }
    }

    // Final check for wall tunnel
    if (wall_enabled && (left_count != right_count)) {
        detect_and_fix_tunnel(X, X_old, Vx, Side_old, Side_new, NUM_PARTICLES, wall_x, WALL_THICKNESS);
    }
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

// Expected kinetic energy for ideal gas: (3/2) N k_B T
double kinetic_energy_expected() {
    return (3.0 / 2.0) * NUM_PARTICLES * K_B * TEMPERATURE;  // In Joules
}

double average_kinetic_energy_per_particle() {
    return kinetic_energy() / NUM_PARTICLES;
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




void log_energy(FILE *fptr) {
    float ke = kinetic_energy();
    float ke_expected = kinetic_energy_expected();
    float avg_ke_pp = average_kinetic_energy_per_particle();
    float T = temperature();


    // Compute kinetic and position entropy
    double S_kinetic = entropy_kinetic();
    double S_position = entropy_position();
    double S_total = S_kinetic + S_position;
    
    printf("Kinetic Entropy: %e J/K\n", S_kinetic);
    printf("Position Entropy: %e J/K\n", S_position);
    printf("Total Entropy: %e J/K\n", S_total);

    compute_velocity_histogram();
    float S_v = entropy();

    compute_position_histogram();
    float S_p = position_entropy();


    float F = total_free_energy(T, S_total);

    fprintf(fptr, "Kinetic Energy: %.8e J, Temperature: %.8e K, Entropy velocities: %.8e J/K, Total Free Energy: %.8e J\n", ke, T, S_v, F);
    // Print the calculated values to the console
    printf("Kinetic Energy: %.6e J, Kinetic Energy expected: %.6e , AVG Kinetic Energy per particle: %.6e, Temperature: %.6e K, Entropy velocity: %.6e J/K, Entropy position: %.6e J/K, Total Free Energy: %.6e J\n", ke,ke_expected, avg_ke_pp, T, S_v, S_p, F);
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


void render_velocity_histograms() {
    int binsX[NUM_BINS] = {0}; // Histogram for Vx
    int binsY[NUM_BINS] = {0}; // Histogram for Vy

    // Calculate velocity bins
    for (int i = 0; i < NUM_PARTICLES; i++) {
        int binX = (int)((Vx[i] + MAX_VELOCITY) / (2 * MAX_VELOCITY) * NUM_BINS);
        int binY = (int)((Vy[i] + MAX_VELOCITY) / (2 * MAX_VELOCITY) * NUM_BINS);
        
        if (binX >= 0 && binX < NUM_BINS) binsX[binX]++;
        if (binY >= 0 && binY < NUM_BINS) binsY[binY]++;
    }

    // Find max bin count for normalization
    int max_binX = 1, max_binY = 1;
    for (int i = 0; i < NUM_BINS; i++) {
        if (binsX[i] > max_binX) max_binX = binsX[i];
        if (binsY[i] > max_binY) max_binY = binsY[i];
    }

    // Render Vx histogram at the bottom
    for (int i = 0; i < NUM_BINS; i++) {
        int height = (binsX[i] * HIST_HEIGHT) / max_binX;
        SDL_SetRenderDrawColor(_renderer, 255, 255, 255, 255); // White bars
        SDL_Rect rect = {XW1 + i * (SIM_WIDTH / NUM_BINS), YW1 + SIM_HEIGHT + (HIST_HEIGHT - height), SIM_WIDTH / NUM_BINS, height};
        SDL_RenderFillRect(_renderer, &rect);
    }

    // Render Vy histogram on the right
    for (int i = 0; i < NUM_BINS; i++) {
        int width = (binsY[i] * HIST_WIDTH) / max_binY;
        SDL_SetRenderDrawColor(_renderer, 255, 255, 255, 255); // White bars
        SDL_Rect rect = {XW1 + SIM_WIDTH, YW1 + i * (SIM_HEIGHT / NUM_BINS), width, SIM_HEIGHT / NUM_BINS};
        SDL_RenderFillRect(_renderer, &rect);
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





////////////////////////  EXPERIMENT FUNCTIONS //////////////////////////////////
void run_speed_of_sound_experiments() {
    system("mkdir -p experiments_speed_of_sound/mode0_real_units");
    system("mkdir -p experiments_speed_of_sound/mode1_normalized_units");
    system("rm -f experiments_speed_of_sound/mode0_real_units/*.csv");
    system("rm -f experiments_speed_of_sound/mode1_normalized_units/*.csv");
    wall_enabled = true;

    /// PARAMETER VARIATIONS
    float lengths[] = {7.5f, 10.0f, 15.0f};  // Box lengths
    int wall_mass_factors[] = {20, 100, 500};  // Wall mass relative to particle mass
    int num_repeats = 1;

    float total_simulation_time = 3000.0f;  // target time after wall release
    int num_steps = (int)(total_simulation_time / (FIXED_DT / (MOLECULAR_MODE ? 1.0f : TAU_time_scale_factor_to_molecular)));
    printf("Running each experiment for %d steps after wall release\n", num_steps);

    for (int l = 0; l < sizeof(lengths) / sizeof(lengths[0]); ++l) {
        for (int m = 0; m < sizeof(wall_mass_factors) / sizeof(wall_mass_factors[0]); ++m) {
            for (int r = 0; r < num_repeats; ++r) {
                float L0 = lengths[l];
                float wall_mass = PARTICLE_MASS * wall_mass_factors[m];
                wall_mass_runtime = wall_mass;

                // --- Initialize wall hold behavior ---
                wall_hold_enabled = true;
                wall_is_released = false;
                steps_elapsed = 0;

                wall_x = L0 * PIXELS_PER_DIAMETER;
                wall_x_old = wall_x;
                vx_wall = 0.0f;

                initialize_simulation();

                #if MOLECULAR_MODE
                    const char* mode_folder = "experiments_speed_of_sound/mode1_normalized_units/";
                #else
                    const char* mode_folder = "experiments_speed_of_sound/mode0_real_units/";
                #endif

                char filename[512];
                int L0_int = (int)(L0 * 10);  // e.g., 7.5 → 75
                sprintf(filename, "%swall_x_positions_L0_%d_wallmassfactor_%d_run%d.csv",
                        mode_folder, L0_int, wall_mass_factors[m], r);

                FILE *wall_log = fopen(filename, "w");
                fprintf(wall_log, "Time, Wall_X, Left_Count, Right_Count\n");
                printf("Running: L0 = %.1f, M = %d * PARTICLE_MASS, run = %d\n", L0, wall_mass_factors[m], r);

                int left_particles = 0;
                int right_particles = 0;

                int recorded_steps = 0;   // This counts only after release

                while (recorded_steps < num_steps) {
                    update_particles_with_substepping(FIXED_DT, &left_particles, &right_particles);
                    update_wall(FIXED_DT);
                    check_wall_collisions(FIXED_DT);

                    // Only log after the wall is released
                    if (wall_is_released) {
                        fprintf(wall_log, "%f, %.10e, %d, %d\n",
                                recorded_steps * FIXED_DT,
                                wall_x * METERS_PER_PIXEL,
                                left_particles,
                                right_particles);
                        recorded_steps++;
                    }

                    if (fabs(vx_wall) < 1e-12 && recorded_steps > 100) {
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
    int running = 1;
    if(live_fft == true){
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
    fprintf(wall_log, "Time, Wall_X\n");
    // 🔥 INIT Wall Hold Behavior for visual mode
    wall_hold_enabled = true;  // or false if you don't want wall hold in visual
    wall_is_released = false;
    steps_elapsed = 0;

    Uint32 last_time = SDL_GetTicks();
    float accumulator = 0.0f;
    
    while (running) {
        Uint32 current_time = SDL_GetTicks();
        float frame_time = (current_time - last_time) / 1000.0f;  // Convert ms to seconds
        last_time = current_time;

        // Prevent spiral of death (if frame_time is too large, limit it)
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

        // Process physics updates at fixed intervals
        while (accumulator >= FIXED_DT) {
            if (simulation_started == 1 && !paused) {
                //update_particles_with_collisions_indexed(FIXED_DT);
                update_particles_with_substepping(FIXED_DT, &left_particles, &right_particles);
                //update_particles_with_substepping_time_windowed_collisions(FIXED_DT);   
                //update_particles(FIXED_DT);
                update_pistons(FIXED_DT);
                //check_wall_collisions();  // Ensure wall-particle collisions
                update_wall(FIXED_DT);  // Move the wall
                check_wall_collisions(FIXED_DT);  // Optional but necessary if wall should affect particles

                // LIVE FFT 
                if (live_fft == true){
                    // Fill buffers
                    in_wall[sample_index] = wall_x;
                    in_energy[sample_index] = kinetic_energy();

                    sample_index = (sample_index + 1) % FFT_SIZE;

                    // Once buffer full
                    if (sample_index == 0) {
                        kiss_fftr(fft_cfg, in_wall, out_wall);
                        kiss_fftr(fft_cfg, in_energy, out_energy);

                        for (int i = 0; i < FFT_SIZE/2 + 1; i++) {
                            wall_fft_magnitude[i] = sqrtf(out_wall[i].r * out_wall[i].r + out_wall[i].i * out_wall[i].i);
                            energy_fft_magnitude[i] = sqrtf(out_energy[i].r * out_energy[i].r + out_energy[i].i * out_energy[i].i);
                        }
                    }
                 }

                // Log wall position
                //fprintf(wall_log, "%f, %f\n", current_time / 1000.0, wall_x); // saved in pixels
                fprintf(wall_log, "%f, %.10e\n", current_time / 1000.0, wall_x * METERS_PER_PIXEL); // in meters no pixels
                log_energy(logFile);
            }
            accumulator -= FIXED_DT;
        }

        // DRAWING AND RENDERING
        draw_clear_screen();
        draw_coordinate_system(_renderer);
        draw_simulation_boundary();
        render_velocity_histograms();
        render_particles();
        render_pistons();
        render_fft_histogram_bottom(_renderer, wall_fft_magnitude, NUM_BINS);
        //draw_fft_bars(_renderer, energy_fft_magnitude, FFT_SIZE / 2 + 1);

        if (wall_enabled) {  // Check if wall should be drawn
            draw_wall();
        }
        // Draw text for counter of particles on each side of the wall
        SDL_Color red = {255, 0, 0, 255};
        SDL_Color blue = {0, 0, 255, 255};

        char buffer[64];

        // Left counter (top-left corner near simulation region)
        sprintf(buffer, "Left: %d", left_count);
        draw_text(_renderer, font, buffer, XW1 + 5, YW1 + 10, red);

        // Right counter (top-right corner near simulation region)
        sprintf(buffer, "Right: %d", right_count);
        draw_text(_renderer, font, buffer, XW2 - 100, YW1 + 10, blue);

        if (wall_hold_enabled && !wall_is_released) {
            SDL_Color yellow = {255, 255, 0, 255};
            char buffer[64];
            sprintf(buffer, "Wall held: %d steps left", wall_hold_steps - steps_elapsed);
            draw_text(_renderer, font, buffer, XW1 + 5, YW1 + 50, yellow);
        }
        SDL_RenderPresent(_renderer);
        SDL_Delay(10);  // Reduce lag while keeping real-time updates
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

    font = TTF_OpenFont("Roboto-Regular.ttf", 18);
    if (!font) {
        printf("Failed to load font: %s\n", TTF_GetError());
        return 1;
    }

    initialize_simulation();
    // Print the first few velocities for testing
    for (int i = 0; i < 10; i++) {
        printf("Particle %d: Vx = %f, Vy = %f\n", i, Vx[i], Vy[i]);
        
    }
   // printf("wall_x = %f, wall_enabled = %d\n", wall_x, wall_enabled);
    if (enable_speed_of_sound_experiments) {
        run_speed_of_sound_experiments();
    } else {
        simulation_loop();  // your SDL interactive version
    }
    SDL_DestroyRenderer(_renderer);
    SDL_DestroyWindow(_window);
    //free(fft_cfg);
    TTF_CloseFont(font);
    TTF_Quit();
    SDL_Quit();
    return 0;
}