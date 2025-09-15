#include <SDL2/SDL.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <SDL2/SDL_ttf.h>


// PARTICLE PARAMETERS
#define NUM_PARTICLES 1 // total number of particles

// --- Physical constants ---
#define PARTICLE_RADIUS 1e-9 // 1 nm radius
#define PARTICLE_MASS 1.67e-27     // kg (proton mass)
#define K_B 1.38e-23f              // J/K
#define TEMPERATURE 300.0f         // Kelvin
#define METERS_PER_PIXEL 1e-9 // pixel size for mapping to meters (important for real mode)




/////////////// MAIN FUNCTION /////////////
// Main Function
int main(int argc, char* argv[]) {

    initialize_simulation_dimensions();     // sets SIM_WIDTH, SIM_HEIGHT, XW2 etc.
    initialize_time_scale();
    initialize_simulation_parameters();     // sets wall and piston regarding simulation dimensions
    print_simulation_info();  // Print simulation parameters

    initSDL();
    TTF_Init();
    
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