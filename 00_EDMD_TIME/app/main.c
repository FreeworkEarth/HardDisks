#include "config_cli.h"
#include "state.h"
#include "step_core.h"
#include "protocol.h"
#include "io_log.h"
#include "ui/viz.h"
#include "ui/control.h"

#include <SDL2/SDL.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void run_speed_of_sound_experiments(void);
static void simulation_loop(void);

int main(int argc, char *argv[]) {
    config_cli_parse(argc, argv);

    initialize_simulation_dimensions();
    initialize_time_scale();
    initialize_simulation_parameters();
    print_simulation_info();

    if (!viz_init(MAX_X, MAX_Y)) {
        fprintf(stderr, "Failed to initialise SDL/TTF.\n");
        return 1;
    }

    if (enable_speed_of_sound_experiments) {
        run_speed_of_sound_experiments();
    } else {
        simulation_loop();
    }

    viz_shutdown();
    return 0;
}

static void simulation_loop(void) {
    initialize_simulation();

    int running = 1;
    if (live_fft) {
        fft_cfg = kiss_fftr_alloc(FFT_SIZE, 0, NULL, NULL);
    }

    SDL_Event e;

    FILE *logFile = fopen("energy_log.csv", "w");
    FILE *wall_log = fopen("wall_position.csv", "w");
    if (!logFile || !wall_log) {
        printf("❌ Failed to open log files.\n");
        exit(1);
    }

    if (log_packing_fraction) {
        log_packing_fractions(wall_log, "main");
    } else {
        fprintf(wall_log, "Time, Wall_X, Displacement(σ), Left_Count, Right_Count\n");
    }
    fprintf(logFile, "KE_particles,KE_wall,KE_total,T,S_v,F\n");

    steps_elapsed        = 0;
    wall_hold_enabled    = true;
    wall_is_released     = false;
    wall_enabled         = 1;
    wall_x_old           = wall_x;
    missed_collision_events = 0;
    worst_penetration_observed = 0.0;

    Uint32 last_time = SDL_GetTicks();
    float accumulator = 0.0f;

    while (running) {
        Uint32 current_time = SDL_GetTicks();
        float frame_time = (current_time - last_time) / 1000.0f;
        last_time = current_time;
        if (frame_time > 0.1f) frame_time = 0.1f;
        accumulator += frame_time;

        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT) running = 0;
            handle_simulation_keys(&e);
            handle_piston_keys(&e);
        }

        int left_particles = 0, right_particles = 0;

        while (accumulator >= FIXED_DT) {
            if (simulation_started && !paused) {
                wall_x_old = wall_x;
                steps_elapsed++;

                wall_impulse_x_accum = 0.f;

                update_particles_with_substepping(FIXED_DT, &left_particles, &right_particles);

                // wall velocity already adjusted inside ccd_wall_step; accumulator is kept for diagnostics.

                update_pistons(FIXED_DT);
                update_wall(FIXED_DT, L0_UNITS);

                simulation_time += FIXED_DT;

                if ((steps_elapsed % 1000) == 0) {
                    double KEp = kinetic_energy();
                    double KEw = 0.5 * WALL_MASS * vx_wall * vx_wall;
                    double Wp  = piston_work_left + piston_work_right;
                    printf("[E] step=%d  KEp=%.6e  KEw=%.6e  KE=%.6e  Pwork=%.6e  (KE-Pwork)=%.6e  hits(L,R)=(%ld,%ld)\n",
                           steps_elapsed, KEp, KEw, KEp + KEw, Wp, (KEp + KEw) - Wp,
                           piston_hits_left, piston_hits_right);
                    if (missed_collision_events > 0) {
                        printf("[CCD] step=%d  potential_misses=%ld  worst_penetration=%.6g\n",
                               steps_elapsed, missed_collision_events, worst_penetration_observed);
                    }
                }

                if ((steps_elapsed % 100) == 0) {
                    log_energy(logFile);
                    fflush(logFile);
                }

                if ((steps_elapsed % 500) == 0) {
                    double ke_p = kinetic_energy();
                    double ke_w = 0.5 * WALL_MASS * vx_wall * vx_wall;
                    printf("[ENERGY] step=%d  KE_p=%.6e  KE_w=%.6e  KE_tot=%.6e  vx_wall=%.6g\n",
                           steps_elapsed, ke_p, ke_w, ke_p + ke_w, vx_wall);
                }

                if ((steps_elapsed % 1000) == 0) {
                    printf("KE=%.9f  vx_wall=%.6g  J_accum=%.6g\n",
                           kinetic_energy(), vx_wall, wall_impulse_x_accum);
                }

                if (wall_is_released && wall_release_time < 0.0) {
                    wall_release_time = simulation_time;
                    printf("✅ Wall released at t = %.3f\n", wall_release_time);
                }

                if (wall_is_released) {
                    double time_after_release = simulation_time - wall_release_time;
                    double wall_x_sigma = wall_x / PIXELS_PER_SIGMA;
                    double disp = wall_x_sigma - (XW1 + XW2) / (2.0 * PIXELS_PER_SIGMA);

                    if (!log_packing_fraction) {
                        fprintf(wall_log, "%.3f, %.3f, %.3f, %d, %d\n",
                                time_after_release, wall_x_sigma, disp,
                                left_particles, right_particles);
                    } else {
                        float particle_area_unitless = (float)(M_PI * PARTICLE_RADIUS_UNIT * PARTICLE_RADIUS_UNIT);
                        float particle_area_actual   = (float)(M_PI * PARTICLE_RADIUS * PARTICLE_RADIUS);
                        float box_area_unitless      = 2.0f * L0_UNITS * HEIGHT_UNITS;

                        float L0_scaled      = L0_UNITS * PIXELS_PER_SIGMA;
                        float height_scaled  = HEIGHT_UNITS * PIXELS_PER_SIGMA;
                        float box_area_scaled= 2.0f * L0_scaled * height_scaled;
                        float wall_area      = 0.5f * WALL_THICKNESS * height_scaled;
                        float wall_buffer_area = 2.0f * wall_area;

                        float packing_unitless     = NUM_PARTICLES * particle_area_unitless / box_area_unitless;
                        float packing_actual_bound = NUM_PARTICLES * particle_area_actual / box_area_scaled;
                        float packing_actual_excl  = NUM_PARTICLES * particle_area_actual / (box_area_scaled - wall_buffer_area);

                        fprintf(wall_log, "%.3f, %.3f, %.3f, %d, %d, %.6f, %.6f, %.6f\n",
                                time_after_release, wall_x_sigma, disp,
                                left_particles, right_particles,
                                packing_unitless, packing_actual_bound, packing_actual_excl);
                    }
                }
            }
            accumulator -= FIXED_DT;
        }

        draw_clear_screen();
        draw_coordinate_system();
        draw_simulation_boundary();
        render_particles();
        render_pistons();
        if (wall_enabled) draw_wall();

        render_energy_hud();

        SDL_Color red = {255, 0, 0, 255};
        SDL_Color blue = {0, 0, 255, 255};
        SDL_Color yellow = {255, 255, 0, 255};

        char buffer[64];
        snprintf(buffer, sizeof(buffer), "Left: %d", left_particles);
        draw_text(_renderer, font, buffer, XW1 + 5, YW1 + 10, red);
        snprintf(buffer, sizeof(buffer), "Right: %d", right_particles);
        draw_text(_renderer, font, buffer, XW2 - 100, YW1 + 10, blue);

        if (wall_hold_enabled && !wall_is_released) {
            snprintf(buffer, sizeof(buffer), "Wall held: %d steps left", wall_hold_steps - steps_elapsed);
            draw_text(_renderer, font, buffer, XW1 + 5, YW1 + 50, yellow);
        }

        render_velocity_histograms();

        SDL_RenderPresent(_renderer);
        SDL_Delay(5);
    }

    fclose(logFile);
    fclose(wall_log);

    if (fft_cfg) {
        free(fft_cfg);
        fft_cfg = NULL;
    }
}

static void run_speed_of_sound_experiments(void) {
    system("mkdir -p experiments_speed_of_sound/mode0_real_units");
    system("mkdir -p experiments_speed_of_sound/mode1_normalized_units");
    system("rm -f experiments_speed_of_sound/mode0_real_units/*.csv");
    system("rm -f experiments_speed_of_sound/mode1_normalized_units/*.csv");

    wall_enabled = 1;
    wall_position_is_managed_externally = 1;

    float lengths[] = {7.5f, 10.0f, 15.0f};
    int wall_mass_factors[] = {20, 50, 100, 200, 500, 1000};
    int num_repeats = 1;

    if (num_steps > 10000000) {
        printf("⚠️ WARNING: Simulation time is very long! Consider reducing the number of steps.\n");
    } else if (num_steps < 1000) {
        printf("⚠️ WARNING: Simulation time is very short, SET EMERGENCY STEPS!\n");
        num_steps = (int)(200000 / FIXED_DT);
    }

    for (size_t l = 0; l < sizeof(lengths) / sizeof(lengths[0]); ++l) {
        for (size_t m = 0; m < sizeof(wall_mass_factors) / sizeof(wall_mass_factors[0]); ++m) {
            for (int r = 0; r < num_repeats; ++r) {
                float L0 = lengths[l];
                wall_mass_runtime = PARTICLE_MASS * wall_mass_factors[m];
                L0_UNITS = L0;
                initialize_simulation_dimensions();

                wall_x     = XW1 + (L0 * PIXELS_PER_SIGMA);
                wall_x_old = wall_x;
                vx_wall    = 0.0f;

                wall_hold_enabled = true;
                wall_is_released  = false;
                steps_elapsed     = 0;
                simulation_time   = 0.0;

                initialize_simulation();

#if MOLECULAR_MODE
                const char* mode_folder = "experiments_speed_of_sound/mode1_normalized_units/";
#else
                const char* mode_folder = "experiments_speed_of_sound/mode0_real_units/";
#endif

                char filename[512];
                int L0_int = (int)(L0 * 10);
                sprintf(filename, "%swall_x_positions_L0_%d_wallmassfactor_%d_run%d.csv",
                        mode_folder, L0_int, wall_mass_factors[m], r);
                FILE *wall_log = fopen(filename, "w");
                if (!wall_log) { printf("❌ Could not open log file.\n"); continue; }

                if (log_packing_fraction) {
                    log_packing_fractions(wall_log, "main");
                } else {
                    fprintf(wall_log, "Time, Wall_X, Displacement(σ), Left_Count, Right_Count\n");
                }

                printf("🔬 Running: L0 = %.1f, M = %d*m, run = %d\n", L0, wall_mass_factors[m], r);
                printf("🔍 Initial wall_x = %.3f, vx_wall = %.6f\n", wall_x, vx_wall);

                int left_particles = 0, right_particles = 0;
                int recorded_steps = 0;
                missed_collision_events = 0;
                worst_penetration_observed = 0.0;

                while (recorded_steps < num_steps) {
                    wall_x_old = wall_x;
                    steps_elapsed++;

                    wall_impulse_x_accum = 0.f;
                    update_particles_with_substepping(FIXED_DT, &left_particles, &right_particles);

                    // wall velocity already updated inside ccd_wall_step; accumulator is diagnostic only.

                    update_pistons(FIXED_DT);
                    update_wall(FIXED_DT, L0_UNITS);

                    simulation_time += FIXED_DT;

                    if (!wall_is_released && steps_elapsed >= wall_hold_steps) {
                        wall_is_released = true;
                        wall_release_time = simulation_time;
                    }

                    if (wall_is_released) {
                        double time_after_release = simulation_time - wall_release_time;
                        double wall_x_sigma = wall_x / PIXELS_PER_SIGMA;
                        double disp = wall_x_sigma - (XW1 + XW2) / (2.0 * PIXELS_PER_SIGMA);

                    if (!log_packing_fraction) {
                        fprintf(wall_log, "%.3f, %.3f, %.3f, %d, %d\n",
                                time_after_release, wall_x_sigma, disp,
                                left_particles, right_particles);
                    } else {
                            float particle_area_unitless = (float)(M_PI * PARTICLE_RADIUS_UNIT * PARTICLE_RADIUS_UNIT);
                            float particle_area_actual   = (float)(M_PI * PARTICLE_RADIUS * PARTICLE_RADIUS);
                            float box_area_unitless      = 2.0f * L0_UNITS * HEIGHT_UNITS;

                            float L0_scaled      = L0_UNITS * PIXELS_PER_SIGMA;
                            float height_scaled  = HEIGHT_UNITS * PIXELS_PER_SIGMA;
                            float box_area_scaled= 2.0f * L0_scaled * height_scaled;
                            float wall_area      = 0.5f * WALL_THICKNESS * height_scaled;
                            float wall_buffer_area = 2.0f * wall_area;

                            float packing_unitless     = NUM_PARTICLES * particle_area_unitless / box_area_unitless;
                            float packing_actual_bound = NUM_PARTICLES * particle_area_actual / box_area_scaled;
                            float packing_actual_excl  = NUM_PARTICLES * particle_area_actual / (box_area_scaled - wall_buffer_area);

                            fprintf(wall_log, "%.3f, %.3f, %.3f, %d, %d, %.6f, %.6f, %.6f\n",
                                    time_after_release, wall_x_sigma, disp,
                                    left_particles, right_particles,
                                    packing_unitless, packing_actual_bound, packing_actual_excl);
                        }
                    }

                    recorded_steps++;

                    if ((recorded_steps % 1000) == 0 && missed_collision_events > 0) {
                        printf("[CCD][exp] step=%d  potential_misses=%ld  worst_penetration=%.6g\n",
                               recorded_steps, missed_collision_events, worst_penetration_observed);
                    }
                }

                fclose(wall_log);
            }
        }
    }
}
