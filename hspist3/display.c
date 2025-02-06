#include "globals.h"
#include "display.h"
#include <SDL2/SDL.h>

// Function to render simulation
void renderSimulation() {
    SDL_SetRenderDrawColor(_renderer, 0, 0, 0, 255);
    SDL_RenderClear(_renderer);

    SDL_SetRenderDrawColor(_renderer, 255, 255, 255, 255);
    for (int i = 0; i < nspheres; i++) {
        SDL_Rect rect = { (int)X[i], (int)Y[i], 5, 5 };
        SDL_RenderFillRect(_renderer, &rect);
    }

    SDL_RenderPresent(_renderer);
}



void Display(SDL_Window* _window, SDL_Renderer* _renderer)
{
    // Initialization
    if (ncount == 0) {
        initialize_simulation();
    }

    // Simulation update loop
    for (long i = 0; i < slowctl; i++)
    {
        update_particles();
        handle_collisions();
        if (!(ncount % NSORT)) {
            sort_and_update_neighbors();
        }

        dobumpsa(X, Y);
        ++ncount;
        total_time = ncount * DELT;

        pproto(); // Presumably some debug output

        update_histograms();
        log_data_if_needed();
    }

    // Rendering
    render_display(_renderer);

    if (!(ncount % keyslowctl)) slowctl = keyslowctl;

    handle_fullscreen(_window);
}