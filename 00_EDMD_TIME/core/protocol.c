#include "protocol.h"

#include <math.h>
#include <stdio.h>

void move_left_piston(float dt) {
    float left_box_edge = 0.0f;
    float piston_width = 20.0f;
    float right_limit = piston_right_x - piston_width;

    piston_left_x += vx_piston_left * dt;

    if (piston_left_x < left_box_edge) {
        piston_left_x = left_box_edge;
        vx_piston_left = 0.0f;
    }
    if (piston_left_x > right_limit) {
        piston_left_x = right_limit;
        vx_piston_left = 0.0f;
    }
}

void move_right_piston(float dt) {
    float right_box_edge = MAX_X - 20.0f;
    float left_limit = piston_left_x + 20.0f;

    piston_right_x += vx_piston_right * dt;

    if (piston_right_x > right_box_edge) {
        piston_right_x = right_box_edge;
        vx_piston_right = 0.0f;
    }
    if (piston_right_x < left_limit) {
        piston_right_x = left_limit;
        vx_piston_right = 0.0f;
    }
}

void update_pistons(float dt) {
    move_left_piston(dt);
    move_right_piston(dt);
}

void update_wall(float dt, float L0_units) {
    if (wall_hold_enabled && !wall_is_released) {
        vx_wall = 0.0f;

        if (!wall_position_is_managed_externally) {
            wall_x = XW1 + (L0_units * PIXELS_PER_SIGMA);
        }

        if (steps_elapsed >= wall_hold_steps) {
            wall_is_released = true;
            wall_release_time = steps_elapsed * dt;
            printf("🔔 Wall released at step %d (t = %.3f)\n", steps_elapsed, wall_release_time);
        }
        return;
    }

    wall_x_old = wall_x;
    wall_x += vx_wall * dt;

    if (steps_elapsed % 100 == 0) {
        printf("📊 Step %d | wall_x = %.6f | vx_wall = %.8f\n", steps_elapsed, wall_x, vx_wall);
        printf("wall_enabled = %d, wall_released = %d, left = %d, right = %d\n",
               wall_enabled, wall_is_released, left_count, right_count);
    }

    if (wall_x < XW1 + 0.1f * L0_units) {
        printf("🛑 Wall hit left boundary, clamping and resetting vx_wall\n");
        wall_x = XW1 + 0.1f * L0_units;
        vx_wall = -vx_wall;
    }
    if (wall_x > XW2 - 0.1f * L0_units) {
        printf("🛑 Wall hit right boundary, clamping and resetting vx_wall\n");
        wall_x = XW2 - 0.1f * L0_units;
        vx_wall = -vx_wall;
    }
}
