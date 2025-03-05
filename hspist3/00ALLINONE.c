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
#define PARTICLE_RADIUS 10  // Particle radius (used for collision detection and rendering)
#define NUM_PARTICLES 30  // Total number of particles in the simulation
#define NUM_BINS 100 
#define MAX_VELOCITY 10000.0
#define K_B 1.38e-23
#define TEMPERATURE 300   // in K Adjust as needed

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
#define SCALE_Vxy 0.08  // Mass of particle (e.g., hydrogen atom) in kg
const float FIXED_DT = 0.0016f;  // 16ms per physics update (matches 60 FPS)

#define Scale_Energy 1e27/1000

// Particle states
float X[NUM_PARTICLES];
float Y[NUM_PARTICLES];
float Vx[NUM_PARTICLES];
float Vy[NUM_PARTICLES];
float Radius[NUM_PARTICLES];

int velocity_histogram[NUM_BINS] = {0};
int position_histogram[NUM_BINS * NUM_BINS] = {0};
int energy_histogram[NUM_BINS] = {0};
int simulation_started = 0;  // Flag to track simulation start


// Global Variables
SDL_Window* _window = NULL;
SDL_Renderer* _renderer = NULL;
float piston_left_x = 0;
float piston_right_x = MAX_X;
float vx_piston_left = 0; // Slow movement for piston velocity
float vx_piston_right = 0; // Slow movement for piston velocity







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
    float sigma = sqrt(K_B * temperature / MASS);  // Standard deviation (velocity scale)
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
    float sigma = sqrt(K_B * temperature / MASS);  // Compute sigma
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





// ONLY 2D - Function to initialize and correct missing definitions
void initialize_simulation() {
    for (int i = 0; i < NUM_PARTICLES; i++) {
        // TODO:check overlapping conditions
        X[i] = rand() % (MAX_X / 2);  // Left part
        Y[i] = rand() % MAX_Y;
        Vx[i] = maxwell_boltzmann_velocity(TEMPERATURE);
        Vy[i] = maxwell_boltzmann_velocity(TEMPERATURE);
        
        //Vx[i] = (rand() % 2 == 0 ? 1 : -1) * (rand() % 5);
        //Vy[i] = (rand() % 2 == 0 ? 1 : -1) * (rand() % 5);  // Added Vy (vertical velocity)#
        printf("Random Vx = %f, Vy = %f\n", Vx[i], Vy[i]);
        Radius[i] = PARTICLE_RADIUS;
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
    float right_box_edge = 800 - 20; // Adjusted for piston width
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














//////////////////// PARTICLE FUNCTIONS --  ////////////////////

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
        float relative_velocity = Vx[i] - vx_piston_left * SCALE_Vxy;
        
        if (relative_velocity < 0) {  // Only if the particle is moving toward the piston
            Vx[i] = -Vx[i] + 2 * vx_piston_left * SCALE_Vxy;
        } else {
            Vx[i] = -Vx[i];  // Simple reflection if moving away
        }

        X[i] = piston_left_x + Radius[i];  // Prevent overlap
    }

    // Check collision with the right piston
    if (X[i] + Radius[i] >= piston_right_x) {  
        float relative_velocity = Vx[i] - vx_piston_right * SCALE_Vxy;

        if (relative_velocity > 0) {  // Only if the piston moves into the particle
            Vx[i] = -Vx[i] + 2 * vx_piston_right * SCALE_Vxy;
        } else {
            Vx[i] = -Vx[i];  // Simple reflection if moving away
        }

        X[i] = piston_right_x - Radius[i];  // Prevent overlap
    }
}




///// ELASTIC COLLISION BETWEEN PARTICLES -- HARD SPHERES ////////
// Update particle positions and handle collisions (boundary and inter-particle)
/* 
basic equations elastic collision:
\mathbf{v_1{\prime}} = \mathbf{v_1} - \frac{2 m_2}{m_1 + m_2} \cdot \frac{(\mathbf{v_1} - \mathbf{v_2}) \cdot (\mathbf{r_1} - \mathbf{r_2})}{|\mathbf{r_1} - \mathbf{r_2}|^2} (\mathbf{r_1} - \mathbf{r_2})
\mathbf{v_2{\prime}} = \mathbf{v_2} - \frac{2 m_1}{m_1 + m_2} \cdot \frac{(\mathbf{v_2} - \mathbf{v_1}) \cdot (\mathbf{r_2} - \mathbf{r_1})}{|\mathbf{r_2} - \mathbf{r_1}|^2} (\mathbf{r_2} - \mathbf{r_1})
*/
void update_particles(float dt) {
    // Loop through each particle
    for (int i = 0; i < NUM_PARTICLES; i++) {
        /*
        •	X[i]  and  Y[i]  are the x and y coordinates of particle  i .
	    •	Vx[i]  and  Vy[i]  are the velocities in the  x  and  y  directions.
	    •	SCALE_Vxy is a factor that adjusts the velocity magnitude for the visualisation. */
        // Update the position of the particle based on its velocity
        X[i] += Vx[i] * dt * SCALE_Vxy;
        Y[i] += Vy[i] * dt * SCALE_Vxy;

    
           
        // Debugging: Print velocity and position
        printf("Particle %d: Vx = %.2f, Vy = %.2f, X = %.2f, Y = %.2f\n", i, Vx[i], Vy[i], X[i], Y[i]);
    
        // Handle boundary collisions - direct reflection at a wall
        handle_piston_collisions(i);
        handle_boundary_collision(i);


        // Check for collisions BETWEEN the current particle i  and all other particles j
        for (int j = i + 1; j < NUM_PARTICLES; j++) {

            // distance between particles i and j
            float dx = X[i] - X[j];  // X-distance between particles
            float dy = Y[i] - Y[j];  // Y-distance between particles
            float dist = sqrt(dx * dx + dy * dy);  // Euclidean distance between particles i and j
            float min_dist = PARTICLE_RADIUS * 2;  // Minimum distance before particles collide (twice the radius)

            // If particles are closer than the minimum distance, they are colliding
            if (dist < min_dist) {

                // normal vector between the two particles (unit vector in the direction!!!! of collision)
                // 	•	The unit vector (nx, ny) represents the collision direction is normalized--> :n_x^2 + n_y^2 = 1
                float nx = dx / dist;  // Normal X-component
                float ny = dy / dist;  // Normal Y-component

                // relative velocity between the particles along the normal direction:  v_{\text{rel}, n} = (\mathbf{v}_i - \mathbf{v}_j) \cdot \mathbf{n}
                /*	•	 v_{\text{rel}, n}  is the relative velocity along the normal.
	                •	 \mathbf{n} = (n_x, n_y)  is the normal direction.
                */
                float relVelX = Vx[i] - Vx[j];  // Relative velocity in the X direction to the Collision Axis
                float relVelY = Vy[i] - Vy[j];  // Relative velocity in the Y direction to the Collision Axis
                float v_rel_dotProduct = (relVelX * nx + relVelY * ny);  //, compute its projection onto the collision normal: Dot product of the relative velocity and normal vector

                // If the particles are moving towards each other (negative dot product)
                if (v_rel_dotProduct < 0) {
                    // Assume equal mass for both particles for homogenous gas (can be adjusted if necessary)
                    float mass_i = MASS;  // Mass of particle i
                    float mass_j = MASS;  // Mass of particle j

                    /* 
                    elastic collision, we use the impulse equation with impulse magnitude J  --- NOT momentum P
                         J_impulse = \frac{2 v_{\text{rel}, n}}{m_1 + m_2} : with J = immpulse magnitude.
                        impulse (change in velocity due to collision) is calculated using the relative velocity along the normal direction and the masses of the particles.
                         formula comes from conservation of momentum and kinetic energy in an elastic collision.
                    */
                    float J_impulse = 2 * v_rel_dotProduct / (mass_i + mass_j);  // Impulse formula for elastic collision

                    // Update velocities for both particles based on impulse and normal direction
                    /*
                        \mathbf{v}_i{\prime} = \mathbf{v}_i - J m_2 \mathbf{n}
                        \mathbf{v}_j{\prime} = \mathbf{v}_j + J m_1 \mathbf{n}
                            •	 \mathbf{v}_i{\prime}  and  \mathbf{v}_j{\prime}  are the new velocities after the collision.
                            •	This follows from Newton’s third law, where equal and opposite forces act on the particles.
                    */
                    Vx[i] -= J_impulse * mass_j * nx ;  // Update velocity of particle i in X direction
                    Vy[i] -= J_impulse * mass_j * ny ;  // Update velocity of particle i in Y direction
                    Vx[j] += J_impulse * mass_i * nx ;  // Update velocity of particle j in X direction
                    Vy[j] += J_impulse * mass_i * ny ;  // Update velocity of particle j in Y direction

                    // Push particles apart to avoid overlap (separate them by half the overlap distance)
                        /*
                        X[i] = X[i] + \frac{n_x \cdot \text{overlap}}{2}
                        X[j] = X[j] - \frac{n_x \cdot \text{overlap}}{2}
                            •	This separates the particles along the normal direction.
                            •	Dividing by 2 ensures that both particles move equally apart.
                        */
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
double kinetic_energy() {
    double ke_total = 0.0;
    for (int i = 0; i < NUM_PARTICLES; i++) {
        double particle_ke = 0.5 * MASS * (Vx[i] * Vx[i] + Vy[i] * Vy[i]);  // In Joules
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
            double f_v = exp(-speed_square / (2 * K_B * TEMPERATURE / MASS));
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
void keysSimulation(SDL_Event e) {
    if (e.type == SDL_KEYDOWN && e.key.keysym.sym == SDLK_SPACE) {
        simulation_started = 1;  // Start the simulation
    }
}

/// Constants for acceleration INIT piston
const float acceleration = 25;  // Increased acceleration
const float max_velocity = 500.0;   // Increased max speed
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




///////////////////// SIMULATION FUNCTIONS //////////////////////

void simulation_loop() {
    int running = 1;
    SDL_Event e;
    FILE *logFile = fopen("energy_log.csv", "w");
    if (!logFile) {
        printf("Error: Unable to open log file!\n");
        return;
    }
    fprintf(logFile, "Kinetic Energy, Temperature, Entropy, Total Free Energy\n");

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

        // Process physics updates at fixed intervals
        while (accumulator >= FIXED_DT) {
            if (simulation_started == 1) {
                update_particles(FIXED_DT);
                update_pistons(FIXED_DT);
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

        SDL_RenderPresent(_renderer);
        SDL_Delay(0.01);  // Reduce lag while keeping real-time updates
    }

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