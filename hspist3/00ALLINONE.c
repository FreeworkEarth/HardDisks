#include <SDL2/SDL.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>



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


// Constants and Definitions
#define PARTICLE_RADIUS 2  // Particle radius (used for collision detection and rendering)
#define NUM_PARTICLES 200  // Total number of particles in the simulation
#define NUM_BINS 100
#define MAX_VELOCITY 10000.0
#define K_B 1.38e-23
#define TEMPERATURE 30  // Adjust as needed


#define MAX_X 800  // Maximum width of the simulation window (in pixels)
#define MAX_Y 600  // Maximum height of the simulation window (in pixels)

#define SIM_SCALE 0.85  // Scaling factor for simulation area
#define XW1 50  // X-coordinate left boundary for the simulation region
#define XW2 ((int)(MAX_X * SIM_SCALE )) // X-coordinate right boundary for the simulation region
#define YW1 0  // Y-coordinate top boundary for the simulation region
#define YW2 ((int)(MAX_Y * SIM_SCALE))  // Y-coordinate bottom boundary for the simulation region

#define SIM_WIDTH   (XW2 - XW1)  // Simulation width
#define SIM_HEIGHT  (YW2 - YW1)  // Simulation height

#define HIST_HEIGHT  (MAX_Y - SIM_HEIGHT)  // Vx histogram space
#define HIST_WIDTH   (MAX_X - SIM_WIDTH)   // Vy histogram space



#define SEED 42  // Seed for random number generation (ensures repeatable randomness)
#define MASS 1.67e-27    // Mass of particle (e.g., hydrogen atom) in kg

#define SCALE_Vxy 0.001     // Mass of particle (e.g., hydrogen atom) in kg
#define Scale_Energy 1e27/1000

// Particle states
float X[NUM_PARTICLES];
float Y[NUM_PARTICLES];
float Vx[NUM_PARTICLES];
float Vy[NUM_PARTICLES];
float Radius[NUM_PARTICLES];

int velocity_histogram[NUM_BINS] = {0};
int energy_histogram[NUM_BINS] = {0};
int simulation_started = 0;  // Flag to track simulation start


// Global Variables
SDL_Window* _window = NULL;
SDL_Renderer* _renderer = NULL;
float piston_left_x = 0;
float piston_right_x = 800.0;
float vx_piston_left = 0.1; // Slow movement for piston velocity
float vx_piston_right = 0.1; // Slow movement for piston velocity




// f(v_x) = \sqrt{\frac{m}{2\pi k_B T}} \exp\left(-\frac{m v_x^2}{2 k_B T}\right)
// f(v) \propto v^2 \exp\left(-\frac{mv^2}{2k_BT}\right)

// Function to generate Maxwell-Boltzmann distributed velocities
float maxwell_boltzmann_velocity(float temperature) {
    float sigma = sqrt(K_B * temperature / MASS);  // Standard deviation (velocity scale)
    printf("sigma = %f\n", sigma);
    float u1, u2, z;

    // Box-Muller transform: Generate Gaussian-distributed value
    do {
        u1 = (float)rand() / RAND_MAX;
        u2 = (float)rand() / RAND_MAX;
        z = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);  // Standard normal random variable
        printf("u1 = %f, u2 = %f, z = %f\n", u1, u2, z);
    } while (isnan(z));  // Ensure numerical stability

    return sigma * z ;  // Scale the Gaussian value by sigma to match Maxwell-Boltzmann distribution
}


// Function to initialize and correct missing definitions
void initialize_simulation() {
    for (int i = 0; i < NUM_PARTICLES; i++) {
        X[i] = rand() % (MAX_X / 6);  // Left part
        Y[i] = rand() % MAX_Y;
        Vx[i] = maxwell_boltzmann_velocity(TEMPERATURE);
        Vy[i] = maxwell_boltzmann_velocity(TEMPERATURE);
        
        //Vx[i] = (rand() % 2 == 0 ? 1 : -1) * (rand() % 5);
        //Vy[i] = (rand() % 2 == 0 ? 1 : -1) * (rand() % 5);  // Added Vy (vertical velocity)#
        printf("Random Vx = %f, Vy = %f\n", Vx[i], Vy[i]);
        Radius[i] = PARTICLE_RADIUS;
    }
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




/////////////////// PISTON FUNCTIONS /////////////////////////

// Function to render pistons
void render_pistons() {
    SDL_SetRenderDrawColor(_renderer, 0, 255, 0, 0);
    SDL_Rect left_piston = { (int)piston_left_x , 0, 5, YW2 };
    SDL_RenderFillRect(_renderer, &left_piston);
    SDL_Rect right_piston = { (int)piston_right_x , 0, 5, YW2 };
    SDL_RenderFillRect(_renderer, &right_piston);
}



void move_left_piston() {
    float left_box_edge = 0;
    float piston_width = 20;
    float right_limit = piston_right_x - piston_width; // Ensure pistons don't overlap

    piston_left_x += vx_piston_left;
    
    if (piston_left_x < left_box_edge) {
        piston_left_x = left_box_edge;
        vx_piston_left = 0;
    }
    if (piston_left_x > right_limit) {
        piston_left_x = right_limit;
        vx_piston_left = 0;
    }
}

void move_right_piston() {
    float right_box_edge = 800 - 20; // Adjusted for piston width
    float left_limit = piston_left_x + 20; // Ensure pistons don't overlap

    piston_right_x += vx_piston_right;
    
    if (piston_right_x > right_box_edge) {
        piston_right_x = right_box_edge;
        vx_piston_right = 0;
    }
    if (piston_right_x < left_limit) {
        piston_right_x = left_limit;
        vx_piston_right = 0;
    }
}

void update_pistons() {
    move_left_piston();
    move_right_piston();
}














//////////////////// PARTICLE FUNCTIONS -- PHYSICS ////////////////////



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
        SDL_SetRenderDrawColor(_renderer, 255, 0, 0, 255);
        for (int i = 0; i < NUM_PARTICLES; i++) {
            SDL_RenderFillCircle(_renderer, (int)X[i], (int)Y[i], PARTICLE_RADIUS);
            }
    }


void handle_piston_collisions() {
    for (int i = 0; i < NUM_PARTICLES; i++) {
        if (X[i] - Radius[i] <= piston_left_x) {  
            Vx[i] = -Vx[i];
            X[i] = piston_left_x + Radius[i];
        }
        if (X[i] + Radius[i] >= piston_right_x) {  
            Vx[i] = -Vx[i];
            X[i] = piston_right_x - Radius[i];
        }
    }
}


// Function to check and handle boundary collisions
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


// Update particle positions and handle collisions (boundary and inter-particle)
void update_particles() {
    // Loop through each particle
    for (int i = 0; i < NUM_PARTICLES; i++) {
        // Update the position of the particle based on its velocity
        X[i] += Vx[i] * SCALE_Vxy;
        Y[i] += Vy[i] * SCALE_Vxy;
           
        // Debugging: Print velocity and position
        printf("Particle %d: Vx = %.2f, Vy = %.2f, X = %.2f, Y = %.2f\n", i, Vx[i], Vy[i], X[i], Y[i]);
    
        // Handle boundary collisions
        handle_boundary_collision(i);


        // Check for collisions between this particle and all other particles
        for (int j = i + 1; j < NUM_PARTICLES; j++) {
            // Calculate the distance between particles i and j
            float dx = X[i] - X[j];  // X-distance between particles
            float dy = Y[i] - Y[j];  // Y-distance between particles
            float dist = sqrt(dx * dx + dy * dy);  // Euclidean distance between particles i and j
            float min_dist = PARTICLE_RADIUS * 2;  // Minimum distance before particles collide (twice the radius)

            // If particles are closer than the minimum distance, they are colliding
            if (dist < min_dist) {
                // Calculate the normal vector between the two particles (unit vector in the direction of collision)
                float nx = dx / dist;  // Normal X-component
                float ny = dy / dist;  // Normal Y-component

                // Calculate the relative velocity between the particles along the normal direction
                float relVelX = Vx[i] - Vx[j];  // Relative velocity in the X direction
                float relVelY = Vy[i] - Vy[j];  // Relative velocity in the Y direction
                float dotProduct = (relVelX * nx + relVelY * ny);  // Dot product of the relative velocity and normal vector

                // If the particles are moving towards each other (negative dot product)
                if (dotProduct < 0) {
                    // Assume equal mass for both particles (can be adjusted if necessary)
                    float mass_i = MASS;  // Mass of particle i
                    float mass_j = MASS;  // Mass of particle j

                    // Calculate impulse (change in velocity due to collision)
                    float impulse = 2 * dotProduct / (mass_i + mass_j);  // Impulse formula for elastic collision

                    // Update velocities for both particles based on impulse and normal direction
                    Vx[i] -= impulse * mass_j * nx;  // Update velocity of particle i in X direction
                    Vy[i] -= impulse * mass_j * ny;  // Update velocity of particle i in Y direction
                    Vx[j] += impulse * mass_i * nx;  // Update velocity of particle j in X direction
                    Vy[j] += impulse * mass_i * ny;  // Update velocity of particle j in Y direction

                    // Push particles apart to avoid overlap (separate them by half the overlap distance)
                    float overlap = min_dist - dist;  // Calculate the amount of overlap
                    X[i] += nx * overlap / 2;  // Push particle i apart along the normal
                    Y[i] += ny * overlap / 2;  // Push particle i apart along the normal
                    X[j] -= nx * overlap / 2;  // Push particle j apart along the normal
                    Y[j] -= ny * overlap / 2;  // Push particle j apart along the normal
                }
            }
        }
    }
}












///////////////////////////   ENERGY CALCULATIONS   ///////////////////////
float kinetic_energy() {
    float ke_total = 0.0;
    for (int i = 0; i < NUM_PARTICLES; i++) {
        // Use unscaled velocities for KE calculation
        float particle_ke = 0.5 * MASS * (Vx[i] * Vx[i] + Vy[i] * Vy[i])  * Scale_Energy;
        printf("Particle %d KE: %.6f\n", i, particle_ke);  // Print kinetic energy for each particle
        ke_total += particle_ke;
    }
    return ke_total;
}

float temperature() {
    return (2.0 / 3.0) * kinetic_energy() / (NUM_PARTICLES * K_B * Scale_Energy) ; 
}




///////////////////////////   ENTROPY CALCULATION   ///////////////////////
float entropy() {
    float S = 0.0;
    int total_counts = 0;
    
    for (int i = 0; i < NUM_BINS; i++) {
        total_counts += velocity_histogram[i];
    }
    
    for (int i = 0; i < NUM_BINS; i++) {
        if (velocity_histogram[i] > 0) {
            float p_i = (float)velocity_histogram[i] / total_counts;
            S -= p_i * log(p_i) * Scale_Energy ;  // Shannon entropy
        }
    }
    return S * K_B;  // Multiply by Boltzmann constant for physical units
}


float total_free_energy(float T, float S) {
    float ke = kinetic_energy();  // Internal energy
    float F = ke - T * S;  // Free energy
    return F;
}


///////////////////////////   PISTON FORCE TRACKING   ///////////////////////
float piston_force(float dt, float *momentum_transfer) {
    return *momentum_transfer / dt;
}

///////////////////////////   LOGGING ENERGY AND TEMPERATURE   ///////////////////////
void log_energy(FILE *fptr) {
    float ke = kinetic_energy();
    float T = temperature();
    float S = entropy();
    float F = total_free_energy(T, S);
    fprintf(fptr, "%f, %f, %f, %f\n", ke, T, S, F);
}










////////// KEYBOARD INPUT HANDLING //////////
void keysSimulation(SDL_Event e) {
    if (e.type == SDL_KEYDOWN && e.key.keysym.sym == SDLK_SPACE) {
        simulation_started = 1;  // Start the simulation
    }
}

/// Constants for acceleration INIT
const float acceleration = 0.5;  // Increased acceleration
const float max_velocity = 2.0;   // Increased max speed
const float deceleration = 0.5;


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


///////////////////// SIMULATION FUNCTIONS /////


int frame_count = 0;
// Updated simulation loop with missing corrections
void simulation_loop() {
    int running = 1;
    SDL_Event e;
    FILE *logFile = fopen("energy_log.csv", "w");  // Log file for energy data
    if (!logFile) {
        printf("Error: Unable to open log file!\n");
        return;
    }
    fprintf(logFile, "Kinetic Energy, Temperature, Entropy, Total Free Energy\n");

    while (running) {
        // Poll events
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT) {
                running = 0;
            }
            keysSimulation(e);
            keysPiston(e);
        }

        
        if (simulation_started == 1) {  // Corrected line
        // Update particle positions and handle collisions
        update_particles();
        handle_piston_collisions();
        update_pistons();
        

        // Update velocity histogram (needed for entropy calculation)
        //update_histogram();


        // Calculate energy values
        float ke = kinetic_energy();      // Calculate kinetic energy
        float T = temperature();         // Calculate temperature from kinetic energy
        float S = entropy();             // Calculate entropy from the velocity distribution
        float fe = total_free_energy(T, S); // Calculate total free energy

        // Print the calculated values to the console
        printf("Kinetic Energy: %.2f, Temperature: %.2f, Entropy: %.2f, Total Free Energy: %.2f\n", ke, T, S, fe);

        // Log the values to the CSV file
        fprintf(logFile, "%f, %f, %f, %f\n", ke, T, S, fe);

        }




        // DRAWING AND RENDERINGS
        draw_clear_screen();
        draw_simulation_boundary();
        render_velocity_histograms();
        render_particles();
        render_pistons();

        // Update the screen
        // Display everything
        SDL_RenderPresent(_renderer);
        SDL_Delay(16); // Frame delay to control speed (about 60 FPS)


    
    
      
    }
    
    // Close log file when simulation ends
    fclose(logFile);

}



/////////////// MAIN FUNCTION /////////////


// Main Function
int main(int argc, char* argv[]) {
    initSDL();
    initialize_simulation();

    // Print the first few velocities for testing
    for (int i = 0; i < 10; i++) {
        printf("Particle %d: Vx = %f, Vy = %f\n", i, Vx[i], Vy[i]);
    }
    simulation_loop();

    SDL_DestroyRenderer(_renderer);
    SDL_DestroyWindow(_window);
    SDL_Quit();
    return 0;
}