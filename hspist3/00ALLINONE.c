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
#define PARTICLE_RADIUS 20  // Particle radius (used for collision detection and rendering)
#define NUM_PARTICLES 20  // Total number of particles in the simulation

#define MAX_X 800  // Maximum width of the simulation window (in pixels)
#define MAX_Y 600  // Maximum height of the simulation window (in pixels)
#define XW1 100  // X-coordinate left boundary for the simulation region
#define XW2 700  // X-coordinate right boundary for the simulation region
#define YW1 100  // Y-coordinate top boundary for the simulation region
#define YW2 500  // Y-coordinate bottom boundary for the simulation region
#define SEED 42  // Seed for random number generation (ensures repeatable randomness)
#define MASS 1.0  // Mass of each particle (used in energy calculations, e.g., kinetic energy)
#define SPRING_CONSTANT 10.0  // Spring constant for the pistons (affects potential energy in the system)

// Particle states
float X[NUM_PARTICLES];
float Y[NUM_PARTICLES];
float Vx[NUM_PARTICLES];
float Vy[NUM_PARTICLES];
float Radius[NUM_PARTICLES];

// Global Variables
SDL_Window* _window = NULL;
SDL_Renderer* _renderer = NULL;
float piston_left_x = 100.0;
float piston_right_x = 700.0;
float vx_piston_left = 0.5; // Slow movement for piston velocity
float vx_piston_right = 0.5; // Slow movement for piston velocity





// Function to handle piston collisions
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



// Function to move pistons
void move_pistons() {
    if (piston_right_x - piston_left_x > 200.0) { // Avoid pistons getting too close
        piston_left_x += vx_piston_left;
        piston_right_x -= vx_piston_right;
    }
}



// Function to render particles
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



// Function to render pistons
void render_pistons() {
    SDL_SetRenderDrawColor(_renderer, 0, 255, 0, 255);
    SDL_Rect left_piston = { (int)piston_left_x - 10, 0, 20, MAX_Y };
    SDL_RenderFillRect(_renderer, &left_piston);
    SDL_Rect right_piston = { (int)piston_right_x - 10, 0, 20, MAX_Y };
    SDL_RenderFillRect(_renderer, &right_piston);
}





// Energy calculations
float kinetic_energy() {
    float ke_total = 0.0;
    for (int i = 0; i < NUM_PARTICLES; i++) {
        ke_total += 0.5 * MASS * Vx[i] * Vx[i];
    }
    return ke_total;
}

float potential_energy() {
    float pe_total = 0.0;
    for (int i = 0; i < NUM_PARTICLES; i++) {
        float distance_left = X[i] - piston_left_x;
        float distance_right = X[i] - piston_right_x;
        pe_total += 0.5 * SPRING_CONSTANT * (distance_left * distance_left + distance_right * distance_right);
    }
    return pe_total;
}

float total_free_energy() {
    float ke = kinetic_energy();
    float pe = potential_energy();
    
    float T = 1.0;  // Arbitrary temperature
    float S = 0.0;  // Ignoring entropy for now

    float F = (ke + pe) - T * S;
    return F;
}





// Function to Initialize SDL
void initSDL() {
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        printf("SDL could not initialize! SDL_Error: %s\n", SDL_GetError());
        exit(1);
    }
    _window = SDL_CreateWindow("Physics Simulation", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, MAX_X, MAX_Y, SDL_WINDOW_SHOWN);
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



// Function to initialize and correct missing definitions
void initialize_simulation() {
    for (int i = 0; i < NUM_PARTICLES; i++) {
        X[i] = rand() % (int)(piston_right_x - piston_left_x) + piston_left_x;
        Y[i] = rand() % MAX_Y;
        Vx[i] = (rand() % 2 == 0 ? 1 : -1) * (rand() % 5);
        Vy[i] = (rand() % 2 == 0 ? 1 : -1) * (rand() % 5);  // Added Vy (vertical velocity)
        Radius[i] = PARTICLE_RADIUS;
    }
}



// Fixed and complete update_particles function
void update_particles(float left_box_edge, float right_box_edge, float top_box_edge, float bottom_box_edge) {
    // Loop through all particles to update positions
    for (int i = 0; i < NUM_PARTICLES; i++) {
        // Update the positions of the particles
        X[i] += Vx[i];
        Y[i] += Vy[i];

        // Check particle collisions with the walls (left and right boundaries)
        if (X[i] - Radius[i] < left_box_edge) {
            X[i] = left_box_edge + Radius[i];  // Keep particle at the boundary
            Vx[i] = -Vx[i];                   // Reverse horizontal velocity (bounce)
        }
        if (X[i] + Radius[i] > right_box_edge) {
            X[i] = right_box_edge - Radius[i]; // Keep particle at the boundary
            Vx[i] = -Vx[i];                   // Reverse horizontal velocity (bounce)
        }
        if (Y[i] - Radius[i] < top_box_edge) {
            Y[i] = top_box_edge + Radius[i];   // Keep particle at the boundary
            Vy[i] = -Vy[i];                   // Reverse vertical velocity (bounce)
        }
        if (Y[i] + Radius[i] > bottom_box_edge) {
            Y[i] = bottom_box_edge - Radius[i]; // Keep particle at the boundary
            Vy[i] = -Vy[i];                    // Reverse vertical velocity (bounce)
        }

        // Particle-particle collision detection and response
        for (int j = i + 1; j < NUM_PARTICLES; j++) {
            float dx = X[i] - X[j];
            float dy = Y[i] - Y[j];
            float dist = sqrt(dx * dx + dy * dy);
            float min_dist = Radius[i] + Radius[j];

            if (dist < min_dist) { // Collision detected
                // Calculate the normal and tangent components of the velocities
                float normal_x = dx / dist;
                float normal_y = dy / dist;
                float relative_velocity_x = Vx[i] - Vx[j];
                float relative_velocity_y = Vy[i] - Vy[j];
                float dot_product = relative_velocity_x * normal_x + relative_velocity_y * normal_y;

                // Only resolve the collision if the particles are moving towards each other
                if (dot_product < 0) {
                    // Update velocities based on elastic collision
                    float impulse = 2 * dot_product / (Radius[i] + Radius[j]);

                    // Apply the impulse to the particles' velocities
                    Vx[i] -= impulse * Radius[j] * normal_x;
                    Vy[i] -= impulse * Radius[j] * normal_y;
                    Vx[j] += impulse * Radius[i] * normal_x;
                    Vy[j] += impulse * Radius[i] * normal_y;

                    // Move particles to prevent overlap (adjust their positions)
                    float overlap = min_dist - dist;
                    X[i] += overlap * (Radius[i] / (Radius[i] + Radius[j])) * normal_x;
                    Y[i] += overlap * (Radius[i] / (Radius[i] + Radius[j])) * normal_y;
                    X[j] -= overlap * (Radius[j] / (Radius[i] + Radius[j])) * normal_x;
                    Y[j] -= overlap * (Radius[j] / (Radius[i] + Radius[j])) * normal_y;
                }
            }
        }
    }
}



// Updated simulation loop with missing corrections
void simulation_loop() {
    int running = 1;
    SDL_Event e;

    // Define boundaries for pistons and the box
    float left_box_edge = 0;
    float right_box_edge = 600;
    float top_box_edge = 0;
    float bottom_box_edge = MAX_Y;
    float piston_width = 20; // Assume pistons have a width of 20 units

    while (running) {
        // Poll events
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT) {
                running = 0;
            }

            // Handle key presses to control piston movement
            if (e.type == SDL_KEYDOWN) {
                switch (e.key.keysym.sym) {
                    case SDLK_UP:    // Increase left piston velocity
                        vx_piston_left += 1.0;
                        break;
                    case SDLK_DOWN:  // Decrease left piston velocity
                        vx_piston_left -= 1.0;
                        break;
                    case SDLK_RIGHT: // Increase right piston velocity
                        vx_piston_right += 1.0;
                        break;
                    case SDLK_LEFT:  // Decrease right piston velocity
                        vx_piston_right -= 1.0;
                        break;
                    case SDLK_w:     // Move left piston faster
                        piston_left_x += 5.0;
                        break;
                    case SDLK_s:     // Move left piston slower
                        piston_left_x -= 5.0;
                        break;
                    case SDLK_d:     // Move right piston faster
                        piston_right_x += 5.0;
                        break;
                    case SDLK_a:     // Move right piston slower
                        piston_right_x -= 5.0;
                        break;
                }
            }
        }

        // Move pistons
        piston_left_x += vx_piston_left;
        piston_right_x += vx_piston_right;

        // Constrain pistons to stay within the box boundaries
        if (piston_left_x < left_box_edge) {
            piston_left_x = left_box_edge;
        }
        if (piston_left_x > right_box_edge - piston_width) {
            piston_left_x = right_box_edge - piston_width;
        }

        if (piston_right_x < left_box_edge + piston_width) {
            piston_right_x = left_box_edge + piston_width;
        }
        if (piston_right_x > right_box_edge) {
            piston_right_x = right_box_edge;
        }

        // Update particles based on new positions
        update_particles(left_box_edge, right_box_edge, top_box_edge, bottom_box_edge);

        // Handle piston collisions with particles
        handle_piston_collisions();

        // Render simulation
        SDL_SetRenderDrawColor(_renderer, 0, 0, 0, 255); // Clear screen with black
        SDL_RenderClear(_renderer);

        // Render particles (red)
        SDL_SetRenderDrawColor(_renderer, 255, 0, 0, 255);
        for (int i = 0; i < NUM_PARTICLES; i++) {
            SDL_RenderFillCircle(_renderer, (int)X[i], (int)Y[i], PARTICLE_RADIUS);
        }

        // Render pistons (green)
        render_pistons();

        // Update the screen
        SDL_RenderPresent(_renderer);
        SDL_Delay(16); // Frame delay to control speed (about 60 FPS)
    }
}






// Main Function
int main(int argc, char* argv[]) {
    initSDL();
    initialize_simulation();
    simulation_loop();

    SDL_DestroyRenderer(_renderer);
    SDL_DestroyWindow(_window);
    SDL_Quit();
    return 0;
}