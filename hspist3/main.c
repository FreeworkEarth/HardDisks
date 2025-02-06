#include "globals.h"
#include "physicssimulation.h"
#include "display.h"
#include "piston.h"
#include <SDL2/SDL.h>

int main(int argc, char* argv[]) {
    SDL_Init(SDL_INIT_VIDEO);
    _window = SDL_CreateWindow("Hard Sphere Gas", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, SIZEX, SIZEY, SDL_WINDOW_SHOWN);
    _renderer = SDL_CreateRenderer(_window, -1, SDL_RENDERER_ACCELERATED);

    while (1) {
        updateParticles();
        handleCollisions();
        movePiston();
        renderSimulation();
        SDL_Delay(16); // ~60 FPS
    }

    SDL_DestroyRenderer(_renderer);
    SDL_DestroyWindow(_window);
    SDL_Quit();
    return 0;
}