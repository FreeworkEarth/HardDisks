#include "step_core.h"
#include "protocol.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static void handle_boundary_collision(int i);
static void handle_piston_collisions(int i, float sub_dt);
static void ccd_wall_step(int idx,
                          float dt,
                          float wall_x0,
                          float *wall_vx_live,
                          int prev_side,
                          float *wall_impulse_x_accum);

void remove_drift(void) {
    float vx_sum = 0.0f, vy_sum = 0.0f;

    for (int i = 0; i < NUM_PARTICLES; i++) {
        vx_sum += (float)Vx[i];
        vy_sum += (float)Vy[i];
    }

    float vx_mean = vx_sum / NUM_PARTICLES;
    float vy_mean = vy_sum / NUM_PARTICLES;

    for (int i = 0; i < NUM_PARTICLES; i++) {
        Vx[i] -= vx_mean;
        Vy[i] -= vy_mean;
    }

    printf("✅ Drift removed: mean Vx = %.4f, mean Vy = %.4f\n", vx_mean, vy_mean);
}

double kinetic_energy(void) {
    double ke_total = 0.0;
    for (int i = 0; i < NUM_PARTICLES; i++) {
        double particle_ke = 0.5 * PARTICLE_MASS * (Vx[i] * Vx[i] + Vy[i] * Vy[i]);
        ke_total += particle_ke;
    }
    return ke_total;
}

float compute_measured_temperature_from_ke(void) {
    double ke_total = kinetic_energy();
    return (float)(ke_total / (NUM_PARTICLES * K_B));
}

void maxwell_boltzmann_2D(float temperature, float *vx, float *vy) {
    float sigma = sqrtf(K_B * temperature / PARTICLE_MASS);
    float u1, u2, z0, z1;

    do {
        u1 = (float)rand() / RAND_MAX;
        u2 = (float)rand() / RAND_MAX;
    } while (u1 <= 1e-10f);

    z0 = sqrtf(-2.0f * logf(u1)) * cosf(2.0f * (float)M_PI * u2);
    z1 = sqrtf(-2.0f * logf(u1)) * sinf(2.0f * (float)M_PI * u2);

    *vx = sigma * z0;
    *vy = sigma * z1;
}

float maxwell_boltzmann_velocity_gaussians(float temperature) {
    float sigma = sqrtf(K_B * temperature / PARTICLE_MASS);
    printf("sigma = %f\n", sigma);
    float u1, u2, z;

    do {
        u1 = (float)rand() / RAND_MAX;
        u2 = (float)rand() / RAND_MAX;
        z = sqrtf(-2.0f * logf(u1)) * cosf(2.0f * (float)M_PI * u2);
        printf("u1 = %f, u2 = %f, z = %f\n", u1, u2, z);
    } while (isnan(z));
    return sigma * z;
}

float sample_gaussian_velocity(float temperature) {
    float sigma = sqrtf(K_B * temperature / PARTICLE_MASS);
    float u1, u2, z;

    do {
        u1 = (float)rand() / RAND_MAX;
        u2 = (float)rand() / RAND_MAX;
        z = sqrtf(-2.0f * logf(u1)) * cosf(2.0f * (float)M_PI * u2);
    } while (isnan(z));

    return sigma * z;
}

void maxwell_boltzmann_velocity_ARRAY_ALL(float temperature) {
    float sigma = sqrtf(K_B * temperature / PARTICLE_MASS);

    for (int i = 0; i < NUM_PARTICLES; i += 2) {
        float u1 = (float)rand() / RAND_MAX;
        float u2 = (float)rand() / RAND_MAX;

        float z0 = sqrtf(-2.0f * logf(u1)) * cosf(2.0f * (float)M_PI * u2);
        float z1 = sqrtf(-2.0f * logf(u1)) * sinf(2.0f * (float)M_PI * u2);

        Vx[i] = sigma * z0;
        Vy[i] = sigma * z1;

        if (i + 1 < NUM_PARTICLES) {
            float u3 = (float)rand() / RAND_MAX;
            float u4 = (float)rand() / RAND_MAX;
            float z2 = sqrtf(-2.0f * logf(u3)) * cosf(2.0f * (float)M_PI * u4);
            float z3 = sqrtf(-2.0f * logf(u3)) * sinf(2.0f * (float)M_PI * u4);
            Vx[i + 1] = sigma * z2;
            Vy[i + 1] = sigma * z3;
        }
    }
}

float maxwell_boltzmann_speed3D(float temperature) {
    float sigma = sqrtf(K_B * temperature / PARTICLE_MASS);
    printf("sigma = %f\n", sigma);

    float u1, u2, u3, z1, z2, z3;

    do {
        u1 = (float)rand() / RAND_MAX;
        u2 = (float)rand() / RAND_MAX;
        u3 = (float)rand() / RAND_MAX;

        z1 = sqrtf(-2.0f * logf(u1)) * cosf(2.0f * (float)M_PI * u2);
        z2 = sqrtf(-2.0f * logf(u1)) * sinf(2.0f * (float)M_PI * u2);
        z3 = sqrtf(-2.0f * logf(u3)) * cosf(2.0f * (float)M_PI * u3);

        printf("u1 = %f, u2 = %f, u3 = %f, z1 = %f, z2 = %f, z3 = %f\n", u1, u2, u3, z1, z2, z3);

    } while (isnan(z1) || isnan(z2) || isnan(z3));

    float v = sigma * sqrtf(z1 * z1 + z2 * z2 + z3 * z3);

    float r1 = (float)rand() / RAND_MAX;
    float r2 = (float)rand() / RAND_MAX;
    float theta = acosf(2.0f * r1 - 1.0f);
    float phi = 2.0f * (float)M_PI * r2;

    float vx = v * sinf(theta) * cosf(phi);
    float vy = v * sinf(theta) * sinf(phi);
    float vz = v * cosf(theta);

    (void)vx;
    (void)vy;
    (void)vz;
    return v;
}

void initialize_simulation_random(void) {
    int left_particles = NUM_PARTICLES / 2;
    int right_particles = NUM_PARTICLES - left_particles;

    float wall_buffer = WALL_THICKNESS * 2.0f;

    for (int i = 0; i < NUM_PARTICLES; i++) {
        int valid_position = 0;
        int attempts = 0;

        while (!valid_position) {
            if (i < left_particles) {
                X[i] = ((float)rand() / RAND_MAX) * (wall_x - wall_buffer - XW1) + XW1 + PARTICLE_RADIUS;
            } else {
                X[i] = ((float)rand() / RAND_MAX) * (XW2 - wall_x - wall_buffer) + wall_x + wall_buffer;
            }

            Y[i] = ((float)rand() / RAND_MAX) * (YW2 - YW1 - 2 * PARTICLE_RADIUS) + YW1 + PARTICLE_RADIUS;

            valid_position = 1;
            for (int j = 0; j < i; j++) {
                float dx = (float)(X[i] - X[j]);
                float dy = (float)(Y[i] - Y[j]);
                float dist_sq = dx * dx + dy * dy;
                if (dist_sq < (4.0f * PARTICLE_RADIUS * PARTICLE_RADIUS)) {
                    valid_position = 0;
                    break;
                }
            }

            attempts++;
            if (attempts > 10000) {
                printf("❌ Could not place particle %d safely after 10000 tries! Exiting.\n", i);
                exit(1);
            }
        }

        if (X[i] > wall_x - wall_buffer && X[i] < wall_x + wall_buffer) {
            printf("⚠️ Particle %d is too close to wall! x = %f\n", i, X[i]);
        }

        Vx[i] = maxwell_boltzmann_velocity_gaussians(TEMPERATURE);
        Vy[i] = maxwell_boltzmann_velocity_gaussians(TEMPERATURE);

        Radius[i] = PARTICLE_RADIUS;

        printf("✅ Particle %d initialized: X = %f, Y = %f, Vx = %f, Vy = %f\n", i, X[i], Y[i], Vx[i], Vy[i]);
    }
}

void initialize_simulation_honeycomb(void) {
    float wall_buffer = WALL_THICKNESS + PARTICLE_RADIUS * 0.20f;
    float left_width = wall_x - wall_buffer - XW1;
    float right_width = XW2 - wall_x - wall_buffer;
    float height = YW2 - YW1;

    int particles_left = NUM_PARTICLES / 2;
    int particles_right = NUM_PARTICLES - particles_left;

    int num_rows_left = (int)sqrtf((float)particles_left);
    int num_cols_left = (particles_left + num_rows_left - 1) / num_rows_left;

    int num_rows_right = (int)sqrtf((float)particles_right);
    int num_cols_right = (particles_right + num_rows_right - 1) / num_rows_right;

    float horizontal_spacing_left = left_width / num_cols_left;
    float vertical_spacing_left = height / num_rows_left;

    float horizontal_spacing_right = right_width / num_cols_right;
    float vertical_spacing_right = height / num_rows_right;

    int idx = 0;

    for (int row = 0; row < num_rows_left; row++) {
        for (int col = 0; col < num_cols_left; col++) {
            if (idx >= particles_left) break;

            float x = XW1 + (col + 0.5f) * horizontal_spacing_left;
            float y = YW1 + (row + 0.5f) * vertical_spacing_left;

            X[idx] = x;
            Y[idx] = y;
            Vx[idx] = maxwell_boltzmann_velocity_gaussians(TEMPERATURE);
            Vy[idx] = maxwell_boltzmann_velocity_gaussians(TEMPERATURE);
            Radius[idx] = PARTICLE_RADIUS;
            idx++;
        }
    }

    for (int row = 0; row < num_rows_right; row++) {
        for (int col = 0; col < num_cols_right; col++) {
            if (idx >= NUM_PARTICLES) break;

            float x = wall_x + wall_buffer + (col + 0.5f) * horizontal_spacing_right;
            float y = YW1 + (row + 0.5f) * vertical_spacing_right;

            X[idx] = x;
            Y[idx] = y;
            Vx[idx] = maxwell_boltzmann_velocity_gaussians(TEMPERATURE);
            Vy[idx] = maxwell_boltzmann_velocity_gaussians(TEMPERATURE);
            Radius[idx] = PARTICLE_RADIUS;
            idx++;
        }
    }

    printf("✅ Finished placing all %d particles evenly across box!\n", NUM_PARTICLES);
    printf("🔍 Initial particle velocities (first 10):\n");
    for (int i = 0; i < 10; i++) {
        printf("Particle %2d: Vx = %.6f, Vy = %.6f\n", i, Vx[i], Vy[i]);
    }
}

void initialize_simulation(void) {
    float wall_buffer = WALL_THICKNESS + DIAMETER * 0.25f;
    float margin = DIAMETER * 1.25f;

    float left_width = wall_x - wall_buffer - XW1;
    float right_width = XW2 - wall_x - wall_buffer;
    float height = YW2 - YW1;

    int particles_left = NUM_PARTICLES / 2;
    int particles_right = NUM_PARTICLES - particles_left;

    int num_rows_left = (int)sqrtf((float)particles_left);
    int num_cols_left = (particles_left + num_rows_left - 1) / num_rows_left;

    int num_rows_right = (int)sqrtf((float)particles_right);
    int num_cols_right = (particles_right + num_rows_right - 1) / num_rows_right;

    float spacing_left_x = left_width / num_cols_left;
    float spacing_left_y = height / num_rows_left;

    float spacing_right_x = right_width / num_cols_right;
    float spacing_right_y = height / num_rows_right;

    int idx = 0;

    for (int row = 0; row < num_rows_left; row++) {
        for (int col = 0; col < num_cols_left; col++) {
            if (idx >= particles_left) break;
            X[idx] = XW1 + margin + col * spacing_left_x;
            Y[idx] = YW1 + margin + row * spacing_left_y;

            float vx, vy;
            maxwell_boltzmann_2D(TEMPERATURE, &vx, &vy);
            Vx[idx] = vx;
            Vy[idx] = vy;
            V_init[idx] = sqrtf(vx * vx + vy * vy);
            Radius[idx] = PARTICLE_RADIUS;
            idx++;
        }
    }

    for (int row = 0; row < num_rows_right; row++) {
        for (int col = 0; col < num_cols_right; col++) {
            if (idx >= NUM_PARTICLES) break;
            X[idx] = XW2 - margin - col * spacing_right_x;
            Y[idx] = YW1 + margin + row * spacing_right_y;
            float vx, vy;
            maxwell_boltzmann_2D(TEMPERATURE, &vx, &vy);
            Vx[idx] = vx;
            Vy[idx] = vy;
            V_init[idx] = sqrtf(vx * vx + vy * vy);
            Radius[idx] = PARTICLE_RADIUS;
            idx++;
        }
    }

    printf("✅ Particles placed with buffer. Wall_x = %.3f\n", wall_x);
    fflush(stdout);
    double actual_ke = kinetic_energy();
    double target_ke = NUM_PARTICLES * K_B * TEMPERATURE;
    double scale = sqrt(target_ke / actual_ke);
    printf("✅ Scaling factor for velocities: %.4f\n", scale);

    for (int i = 0; i < NUM_PARTICLES; i++) {
        Vx[i] *= scale;
        Vy[i] *= scale;
    }
}

static void handle_boundary_collision(int i) {
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

static void handle_piston_collisions(int i, float sub_dt) {
    (void)sub_dt;
    const float eps = 1e-6f;

    if (X[i] - Radius[i] <= piston_left_x) {
        float vrel = (float)Vx[i] - vx_piston_left;
        if (vrel < 0.f) {
            double u = Vx[i];
            double v = 2.f * vx_piston_left - Vx[i];
            double dE = 0.5 * PARTICLE_MASS * (v * v - u * u);
            piston_work_left += dE;
            piston_hits_left++;
            Vx[i] = v;
        }
        X[i] = piston_left_x + Radius[i] + eps;
    }

    if (X[i] + Radius[i] >= piston_right_x) {
        float vrel = (float)Vx[i] - vx_piston_right;
        if (vrel > 0.f) {
            double u = Vx[i];
            double v = 2.f * vx_piston_right - Vx[i];
            double dE = 0.5 * PARTICLE_MASS * (v * v - u * u);
            piston_work_right += dE;
            piston_hits_right++;
            Vx[i] = v;
        }
        X[i] = piston_right_x - Radius[i] - eps;
    }
}

static void ccd_wall_step(int idx,
                          float dt,
                          float wall_x0,
                          float *wall_vx_live,
                          int prev_side,
                          float *wall_impulse_accum) {
    if (!wall_enabled) {
        X[idx] += Vx[idx] * dt;
        Y[idx] += Vy[idx] * dt;
        return;
    }

    const bool locked_wall = (wall_hold_enabled && !wall_is_released);

    float xs  = (float)X[idx];
    float ys  = (float)Y[idx];
    float vxs = (float)Vx[idx];
    float vys = (float)Vy[idx];
    const float R = (float)Radius[idx];

    const float wL0 = wall_x0 - 0.5f * WALL_THICKNESS;
    const float wR0 = wall_x0 + 0.5f * WALL_THICKNESS;
    float vwx = wall_vx_live ? *wall_vx_live : vx_wall;

    const float L0p = wL0 - R;
    const float R0p = wR0 + R;

    const float eps_pos = fmaxf(1e-8f, 1e-6f * fmaxf(R, 1.f));
    const float eps_t   = 1e-12f;
    const float v_hys   = 1e-8f;

    float t_remain = dt;

    for (int iter = 0; iter < 6 && t_remain > eps_t; ++iter) {
        float vrel = vxs - vwx;

        const bool inside = (xs > L0p - eps_pos) && (xs < R0p + eps_pos);

        if (fabsf(vrel) < v_hys) {
            if (inside) {
                int face = (prev_side != 0)
                         ? prev_side
                         : (fabsf(xs - L0p) <= fabsf(R0p - xs) ? -1 : +1);
                xs = (face == -1) ? (L0p - eps_pos) : (R0p + eps_pos);
            }
            xs += vxs * t_remain;
            ys += vys * t_remain;
            t_remain = 0.f;
            break;
        }

        if (inside) {
            int face = (prev_side != 0)
                     ? prev_side
                     : (fabsf(xs - L0p) <= fabsf(R0p - xs) ? -1 : +1);

            const float m = (float)PARTICLE_MASS;
            const float M = (float)WALL_MASS;
            const float u1 = vxs;
            const float u2 = vwx;

            float v1;
            if (locked_wall) {
                v1 = 2.f * u2 - u1;
            } else {
                v1 = ((m - M) * u1 + 2.f * M * u2) / (m + M);
                float Jp = m * (v1 - u1);
                if (wall_impulse_accum) {
                    *wall_impulse_accum -= Jp;
                }
                if (wall_vx_live) {
                    *wall_vx_live -= Jp / WALL_MASS;
                    vwx = *wall_vx_live;
                }
            }
            vxs = v1;

            xs = (face == -1) ? (L0p - eps_pos) : (R0p + eps_pos);
            xs += vxs * 1e-12f;
            continue;
        }

        float thit = INFINITY;
        int face = 0;

        if (vrel > 0.f) {
            float tL = (L0p - xs) / vrel;
            if (tL >= -eps_t && tL <= t_remain + eps_t) { thit = tL; face = -1; }
        }
        if (vrel < 0.f) {
            float tR = (R0p - xs) / vrel;
            if (tR >= -eps_t && tR <= t_remain + eps_t && tR < thit) { thit = tR; face = +1; }
        }

        if (face == 0 || !isfinite(thit) || thit > t_remain - eps_t) {
            xs += vxs * t_remain;
            ys += vys * t_remain;
            t_remain = 0.f;
            break;
        }

        xs += vxs * thit;
        ys += vys * thit;

        {
            const float m = (float)PARTICLE_MASS;
            const float M = (float)WALL_MASS;
            const float u1 = vxs;
            const float u2 = vwx;

            float v1;
            if (locked_wall) {
                v1 = 2.f * u2 - u1;
            } else {
                v1 = ((m - M) * u1 + 2.f * M * u2) / (m + M);
                float Jp = m * (v1 - u1);
                if (wall_impulse_accum) {
                    *wall_impulse_accum -= Jp;
                }
                if (wall_vx_live) {
                    *wall_vx_live -= Jp / WALL_MASS;
                    vwx = *wall_vx_live;
                }
            }
            vxs = v1;
        }

        xs = (face == -1) ? (L0p - eps_pos) : (R0p + eps_pos);
        xs += vxs * 1e-12f;

        t_remain -= thit;
    }

    X[idx]  = xs;
    Y[idx]  = ys;
    Vx[idx] = vxs;
    Vy[idx] = vys;
}

/*
 * update_particles_with_substepping
 * ----------------------------------
 * Integrates every disk across the given time step using optional substeps,
 * running CCD against the divider wall, piston impacts, boundary reflections,
 * and neighbour-limited pair collisions. Any overlap larger than a small
 * tolerance is recorded so we can keep track of potential missed collisions.
 */
void update_particles_with_substepping(float dt, int *out_left, int *out_right) {
    if (dt <= 0.0f) return;

    float sub_dt = dt;

#if SUBSTEPPING
    sub_dt = dt / SUBSTEPS;
#endif

    int Side_old[NUM_PARTICLES];
    left_count = right_count = 0;

    const float cell_size = fmaxf(4.0f * PARTICLE_RADIUS, 1.0f);
    const float missed_threshold = PARTICLE_RADIUS * 0.05f;
    const float domain_width  = (float)(XW2 - XW1);
    const float domain_height = (float)(YW2 - YW1);
    int grid_cols = (int)ceilf(domain_width  / cell_size);
    int grid_rows = (int)ceilf(domain_height / cell_size);
    if (grid_cols < 1) grid_cols = 1;
    if (grid_rows < 1) grid_rows = 1;
    int max_cells = grid_cols * grid_rows;

    int *cell_head = (int*)malloc((size_t)max_cells * sizeof(int));
    int *cell_next = (int*)malloc((size_t)NUM_PARTICLES * sizeof(int));
    int *particle_cell = (int*)malloc((size_t)NUM_PARTICLES * sizeof(int));
    if (!cell_head || !cell_next || !particle_cell) {
        printf("❌ Failed to allocate collision grid.\n");
        exit(1);
    }

    for (int step = 0; step < (int)SUBSTEPS; ++step) {
        for (int i = 0; i < NUM_PARTICLES; i++) {
            X_old[i] = X[i];

            if (wall_enabled) {
                if (X_old[i] < wall_x - WALL_THICKNESS / 2)      Side_old[i] = -1;
                else if (X_old[i] > wall_x + WALL_THICKNESS / 2) Side_old[i] =  1;
                else                                             Side_old[i] =  0;
            } else {
                Side_old[i] = 0;
            }

            const float w_x0_sub = wall_x_old;
            const int prev_side = Side_old[i];

            ccd_wall_step(i, sub_dt, w_x0_sub, &vx_wall, prev_side, &wall_impulse_x_accum);

            handle_piston_collisions(i, sub_dt);
            handle_boundary_collision(i);
        }

        for (int c = 0; c < max_cells; ++c) cell_head[c] = -1;

        for (int i = 0; i < NUM_PARTICLES; ++i) {
            float x_local = (float)(X[i] - XW1);
            float y_local = (float)(Y[i] - YW1);
            int cx = (int)floorf(x_local / cell_size);
            int cy = (int)floorf(y_local / cell_size);
            if (cx < 0) cx = 0; else if (cx >= grid_cols) cx = grid_cols - 1;
            if (cy < 0) cy = 0; else if (cy >= grid_rows) cy = grid_rows - 1;
            int cell_index = cy * grid_cols + cx;
            particle_cell[i] = cell_index;
            cell_next[i] = cell_head[cell_index];
            cell_head[cell_index] = i;
        }

        for (int i = 0; i < NUM_PARTICLES; ++i) {
            int cell_index = particle_cell[i];
            int cx = cell_index % grid_cols;
            int cy = cell_index / grid_cols;

            for (int dy = -1; dy <= 1; ++dy) {
                int ny = cy + dy;
                if (ny < 0 || ny >= grid_rows) continue;
                for (int dx = -1; dx <= 1; ++dx) {
                    int nx = cx + dx;
                    if (nx < 0 || nx >= grid_cols) continue;
                    int neighbor_idx = ny * grid_cols + nx;
                    for (int j = cell_head[neighbor_idx]; j != -1; j = cell_next[j]) {
                        if (j <= i) continue;

                        float dxp = (float)(X[i] - X[j]);
                        float dyp = (float)(Y[i] - Y[j]);
                        float dist2 = dxp * dxp + dyp * dyp;
                        float min_d = 2.f * PARTICLE_RADIUS;

                        if (dist2 < min_d * min_d && dist2 > 0.f) {
                            float dist = sqrtf(dist2);
                            float nxp = dxp / dist;
                            float nyp = dyp / dist;

                            float relx = (float)(Vx[i] - Vx[j]);
                            float rely = (float)(Vy[i] - Vy[j]);
                            float vn = relx * nxp + rely * nyp;
                            if (vn < 0.f) {
                                Vx[i] -= vn * nxp;
                                Vy[i] -= vn * nyp;
                                Vx[j] += vn * nxp;
                                Vy[j] += vn * nyp;

                                float overlap = min_d - dist;
                                if (overlap > missed_threshold) {
                                    missed_collision_events++;
                                    if (overlap > worst_penetration_observed) {
                                        worst_penetration_observed = overlap;
                                    }
                                }
                                float push = 0.5f * overlap;
                                X[i] += nxp * push;  Y[i] += nyp * push;
                                X[j] -= nxp * push;  Y[j] -= nyp * push;
                            }
                        }
                    }
                }
            }
        }
    }

    free(cell_head);
    free(cell_next);
    free(particle_cell);

    if (wall_enabled) {
        left_count = right_count = 0;
        for (int i = 0; i < NUM_PARTICLES; i++) {
            if (X[i] < wall_x - WALL_THICKNESS / 2) left_count++;
            else if (X[i] > wall_x + WALL_THICKNESS / 2) right_count++;
        }
    }

    if (out_left)  *out_left  = left_count;
    if (out_right) *out_right = right_count;
}

double kinetic_energy_expected(void) {
    return (3.0 / 2.0) * NUM_PARTICLES * K_B * TEMPERATURE;
}

double average_kinetic_energy_per_particle(void) {
    return kinetic_energy() / NUM_PARTICLES;
}

double kinetic_energy_total_system(void) {
    double ke_particles = kinetic_energy();
    double ke_wall = 0.5 * WALL_MASS * vx_wall * vx_wall;
    return ke_particles + ke_wall;
}

double temperature(void) {
    return (2.0 / 2.0) * kinetic_energy() / (NUM_PARTICLES * K_B);
}

void compute_velocity_histogram(void) {
    for (int i = 0; i < NUM_BINS; i++) velocity_histogram[i] = 0;

    for (int i = 0; i < NUM_PARTICLES; i++) {
        float speed = (float)sqrt(Vx[i] * Vx[i] + Vy[i] * Vy[i]);
        int bin = (int)(speed / MAX_VELOCITY * NUM_BINS);
        if (bin >= 0 && bin < NUM_BINS) velocity_histogram[bin]++;
    }
}

void compute_position_histogram(void) {
    for (int i = 0; i < NUM_BINS * NUM_BINS; i++) position_histogram[i] = 0;

    for (int i = 0; i < NUM_PARTICLES; i++) {
        int bin_x = (int)(X[i] / MAX_X * NUM_BINS);
        int bin_y = (int)(Y[i] / MAX_Y * NUM_BINS);
        int bin = bin_x + bin_y * NUM_BINS;

        if (bin >= 0 && bin < NUM_BINS * NUM_BINS) position_histogram[bin]++;
    }
}

double entropy_kinetic(void) {
    double S_kinetic = 0.0;

    for (int i = 0; i < NUM_PARTICLES; i++) {
        double speed_square = Vx[i] * Vx[i] + Vy[i] * Vy[i];
        if (speed_square > 0) {
            double f_v = exp(-speed_square / (2 * K_B * TEMPERATURE / PARTICLE_MASS));
            S_kinetic -= f_v * log(f_v);
        }
    }
    return S_kinetic * K_B;
}

double entropy_position(void) {
    return K_B * log((XW2 - XW1) * (YW2 - YW1));
}

float entropy(void) {
    double S = 0.0;
    int total_counts = 0;

    for (int i = 0; i < NUM_BINS; i++) {
        total_counts += velocity_histogram[i];
    }

    if (total_counts == 0) return 0.0f;

    for (int i = 0; i < NUM_BINS; i++) {
        if (velocity_histogram[i] > 0) {
            double p_i = (double)velocity_histogram[i] / total_counts;
            S -= p_i * log(p_i);
        }
    }
    return (float)(S * K_B);
}

float position_entropy(void) {
    double S = 0.0;
    int total_counts = 0;

    for (int i = 0; i < NUM_BINS * NUM_BINS; i++) {
        total_counts += position_histogram[i];
    }

    if (total_counts == 0) return 0.0f;

    for (int i = 0; i < NUM_BINS * NUM_BINS; i++) {
        if (position_histogram[i] > 0) {
            double p_i = (double)position_histogram[i] / total_counts;
            S -= p_i * log(p_i);
        }
    }
    return (float)(S * K_B);
}

double total_free_energy(double T, double S) {
    double ke = kinetic_energy();
    double F = ke - T * S;
    return F;
}
