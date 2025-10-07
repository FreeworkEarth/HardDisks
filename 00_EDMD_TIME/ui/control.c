#include "control.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static const float acceleration = 25.0f;
static const float max_velocity = 1000000.0f;
static const float deceleration = 25.0f;

void handle_simulation_keys(const SDL_Event *e) {
    if (e->type == SDL_KEYDOWN) {
        switch (e->key.keysym.sym) {
            case SDLK_SPACE:
                simulation_started = 1;
                break;

            case SDLK_q:
                exit(0);
                break;

            case SDLK_w:
                wall_enabled = !wall_enabled;
                printf("Wall %s\n", wall_enabled ? "Enabled" : "Disabled");
                break;

            case SDLK_r:
                if (wall_hold_enabled && !wall_is_released) {
                    wall_is_released = true;
                    printf("🔔 Wall manually released by key 'r'.\n");
                }
                break;

            case SDLK_p:
                paused = !paused;
                printf("Simulation %s\n", paused ? "Paused" : "Running");
                break;
        }
    }
}

void handle_piston_keys(const SDL_Event *e) {
    if (e->type == SDL_KEYDOWN) {
        switch (e->key.keysym.sym) {
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
            case SDLK_s:
                vx_piston_left = 0;
                break;

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
            case SDLK_UP:
                vx_piston_right = 0;
                break;
        }
    } else if (e->type == SDL_KEYUP) {
        switch (e->key.keysym.sym) {
            case SDLK_a:
            case SDLK_d:
                if (vx_piston_left > 0) vx_piston_left -= deceleration;
                if (vx_piston_left < 0) vx_piston_left += deceleration;
                if (fabsf(vx_piston_left) < deceleration) vx_piston_left = 0;
                break;

            case SDLK_LEFT:
            case SDLK_RIGHT:
                if (vx_piston_right > 0) vx_piston_right -= deceleration;
                if (vx_piston_right < 0) vx_piston_right += deceleration;
                if (fabsf(vx_piston_right) < deceleration) vx_piston_right = 0;
                break;
        }
    }
}
