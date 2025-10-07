#ifndef UI_CONTROL_H
#define UI_CONTROL_H

#include <SDL2/SDL.h>
#include "state.h"

#ifdef __cplusplus
extern "C" {
#endif

void handle_simulation_keys(const SDL_Event *e);
void handle_piston_keys(const SDL_Event *e);

#ifdef __cplusplus
}
#endif

#endif /* UI_CONTROL_H */
