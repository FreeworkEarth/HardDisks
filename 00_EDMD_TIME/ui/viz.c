#include "viz.h"
#include "step_core.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

bool viz_init(int width, int height) {
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        printf("SDL could not initialize! SDL_Error: %s\n", SDL_GetError());
        return false;
    }

    if (TTF_Init() < 0) {
        printf("TTF could not initialize! SDL_ttf Error: %s\n", TTF_GetError());
        return false;
    }

    _window = SDL_CreateWindow("Hard Spheres Simulation",
                               SDL_WINDOWPOS_UNDEFINED,
                               SDL_WINDOWPOS_UNDEFINED,
                               width,
                               height,
                               SDL_WINDOW_SHOWN);
    if (!_window) {
        printf("Window could not be created! SDL_Error: %s\n", SDL_GetError());
        return false;
    }

    _renderer = SDL_CreateRenderer(_window, -1, SDL_RENDERER_ACCELERATED);
    if (!_renderer) {
        printf("Renderer could not be created! SDL_Error: %s\n", SDL_GetError());
        return false;
    }

    font = TTF_OpenFont("Roboto-Regular.ttf", 18);
    if (!font) {
        printf("Failed to load font: %s\n", TTF_GetError());
        return false;
    }

    return true;
}

void viz_shutdown(void) {
    if (font) {
        TTF_CloseFont(font);
        font = NULL;
    }
    if (_renderer) {
        SDL_DestroyRenderer(_renderer);
        _renderer = NULL;
    }
    if (_window) {
        SDL_DestroyWindow(_window);
        _window = NULL;
    }
    TTF_Quit();
    SDL_Quit();
}

static void SDL_RenderFillCircle(SDL_Renderer *renderer, int centerX, int centerY, int radius_pixels) {
    for (int w = 0; w < radius_pixels * 2; w++) {
        for (int h = 0; h < radius_pixels * 2; h++) {
            int dx = radius_pixels - w;
            int dy = radius_pixels - h;
            if ((dx * dx + dy * dy) <= (radius_pixels * radius_pixels)) {
                SDL_RenderDrawPoint(renderer, centerX + dx, centerY + dy);
            }
        }
    }
}

void render_particles(void) {
    SDL_SetRenderDrawColor(_renderer, 255, 255, 255, 255);

    for (int i = 0; i < NUM_PARTICLES; i++) {
        int pixel_x = (int)(X[i]);
        int pixel_y = (int)(Y[i]);
        int radius_pixels = (int)(PARTICLE_RADIUS);

        SDL_RenderFillCircle(_renderer, pixel_x, pixel_y, radius_pixels);
    }
}

void render_pistons(void) {
    SDL_SetRenderDrawColor(_renderer, 0, 255, 0, 150);
    SDL_Rect left_piston = { (int)piston_left_x , 0, 5, YW2 };
    SDL_RenderFillRect(_renderer, &left_piston);
    SDL_Rect right_piston = { (int)piston_right_x , 0, 5, YW2 };
    SDL_RenderFillRect(_renderer, &right_piston);
}

void draw_wall(void) {
    if (!wall_enabled) return;
    if (isnan(wall_x) || wall_x < 0 || wall_x > SIM_WIDTH) {
        printf("🚨 Invalid wall_x = %f\n", wall_x);
        exit(1);
    }
    SDL_SetRenderDrawColor(_renderer, 255, 255, 255, 255);
    SDL_Rect wall_rect = { (int)(wall_x - (WALL_THICKNESS/2)), 0, (int)WALL_THICKNESS, SIM_HEIGHT };
    SDL_RenderFillRect(_renderer, &wall_rect);
}

void draw_clear_screen(void) {
    SDL_SetRenderDrawColor(_renderer, 0, 0, 0, 255);
    SDL_RenderClear(_renderer);
}

void draw_simulation_boundary(void) {
    SDL_SetRenderDrawColor(_renderer, 255, 255, 0, 255);
    SDL_Rect sim_border = {XW1, YW1, SIM_WIDTH, SIM_HEIGHT};
    SDL_RenderDrawRect(_renderer, &sim_border);
}

void draw_coordinate_system(void) {
    SDL_SetRenderDrawColor(_renderer, 255, 0, 0, 255);

    int origin_x = 20;
    int origin_y = 20;
    int axis_length = 50;

    SDL_RenderDrawLine(_renderer, origin_x, origin_y, origin_x + axis_length, origin_y);
    SDL_RenderDrawLine(_renderer, origin_x + axis_length, origin_y, origin_x + axis_length - 5, origin_y - 5);
    SDL_RenderDrawLine(_renderer, origin_x + axis_length, origin_y, origin_x + axis_length - 5, origin_y + 5);

    SDL_RenderDrawLine(_renderer, origin_x, origin_y, origin_x, origin_y + axis_length);
    SDL_RenderDrawLine(_renderer, origin_x, origin_y + axis_length, origin_x - 5, origin_y + axis_length - 5);
    SDL_RenderDrawLine(_renderer, origin_x, origin_y + axis_length, origin_x + 5, origin_y + axis_length - 5);
}

void draw_text(SDL_Renderer *renderer, TTF_Font *font_ref, const char *text, int x, int y, SDL_Color color) {
    SDL_Surface *text_surface = TTF_RenderText_Blended(font_ref, text, color);
    SDL_Texture *text_texture = SDL_CreateTextureFromSurface(renderer, text_surface);
    SDL_Rect dst = { x, y, text_surface->w, text_surface->h };
    SDL_RenderCopy(renderer, text_texture, NULL, &dst);
    SDL_FreeSurface(text_surface);
    SDL_DestroyTexture(text_texture);
}

static void draw_maxwell_boltzmann_curve(float temperature, int max_count) {
    float sigma = sqrtf(K_B * temperature / PARTICLE_MASS);
    float norm_factor = NUM_PARTICLES * (SIM_WIDTH / NUM_BINS);

    for (int i = 0; i < NUM_BINS; i++) {
        float v = MAX_VELOCITY * i / NUM_BINS;
        float f_v = (v * v) * expf(-PARTICLE_MASS * v * v / (2.0f * K_B * temperature));
        float f_scaled = f_v * norm_factor;

        int height = (int)((f_scaled * HIST_HEIGHT) / max_count);

        SDL_SetRenderDrawColor(_renderer, 255, 255, 0, 255);
        SDL_RenderDrawPoint(_renderer, XW1 + i * (SIM_WIDTH / NUM_BINS), YW1 + SIM_HEIGHT + (HIST_HEIGHT - height));
    }
}

static void draw_thick_line_vertical(SDL_Renderer *renderer, int x, int y1, int y2, int thickness) {
    if (y1 > y2) {
        int temp = y1;
        y1 = y2;
        y2 = temp;
    }
    SDL_Rect rect = {x - thickness / 2, y1, thickness, y2 - y1 + 1};
    SDL_RenderFillRect(renderer, &rect);
}

void render_velocity_histograms(void) {
    SDL_SetRenderDrawBlendMode(_renderer, SDL_BLENDMODE_BLEND);

    float measured_T = compute_measured_temperature_from_ke();

    float center_x = XW1 + L0_UNITS * PIXELS_PER_SIGMA;
    int bin_width = SIM_WIDTH / NUM_BINS;

    int binsX[NUM_BINS] = {0};
    int binsY[NUM_BINS] = {0};
    int binsSpeed[NUM_BINS] = {0};

    int max_binX = 1, max_binY = 1, max_binSpeed = 1;

    for (int i = 0; i < NUM_PARTICLES; i++) {
        int binX = (int)((Vx[i] + MAX_VELOCITY) / (2 * MAX_VELOCITY) * NUM_BINS);
        int binY = (int)((Vy[i] + MAX_VELOCITY) / (2 * MAX_VELOCITY) * NUM_BINS);
        float speed = (float)sqrt(Vx[i] * Vx[i] + Vy[i] * Vy[i]);
        int binV = (int)(speed / MAX_VELOCITY * NUM_BINS);

        if (binX >= 0 && binX < NUM_BINS) binsX[binX]++;
        if (binY >= 0 && binY < NUM_BINS) binsY[binY]++;
        if (binV >= 0 && binV < NUM_BINS) binsSpeed[binV]++;
    }

    for (int i = 0; i < NUM_BINS; i++) {
        if (binsX[i] > max_binX) max_binX = binsX[i];
        if (binsY[i] > max_binY) max_binY = binsY[i];
        if (binsSpeed[i] > max_binSpeed) max_binSpeed = binsSpeed[i];
    }

    for (int i = 0; i < NUM_BINS; i++) {
        int height = (binsX[i] * HIST_HEIGHT) / max_binX;
        SDL_SetRenderDrawColor(_renderer, 255, 255, 255, 150);
        SDL_Rect rect = {
            (int)(center_x + (i - NUM_BINS / 2) * bin_width),
            YW1 + SIM_HEIGHT + (HIST_HEIGHT - height),
            bin_width,
            height
        };
        SDL_RenderFillRect(_renderer, &rect);
    }

    for (int i = 0; i < NUM_BINS; i++) {
        int height = (binsSpeed[i] * HIST_HEIGHT) / max_binSpeed;
        SDL_SetRenderDrawColor(_renderer, 255, 255, 0, 255);
        SDL_Rect rect = {
            (int)(center_x + (i - NUM_BINS / 2) * bin_width),
            YW1 + SIM_HEIGHT + (HIST_HEIGHT - height),
            bin_width,
            height
        };
        SDL_RenderFillRect(_renderer, &rect);
    }

    for (int i = 0; i < NUM_BINS; i++) {
        int width = (binsY[i] * HIST_WIDTH) / max_binY;
        SDL_SetRenderDrawColor(_renderer, 100, 200, 255, 180);
        SDL_Rect rect = {
            XW1 + SIM_WIDTH - HIST_WIDTH,
            YW1 + i * (SIM_HEIGHT / NUM_BINS),
            width,
            SIM_HEIGHT / NUM_BINS
        };
        SDL_RenderFillRect(_renderer, &rect);
    }

    float normalization = NUM_PARTICLES * bin_width * 50.0f;
    int max_count = max_binSpeed ? max_binSpeed : 1;
    draw_maxwell_boltzmann_curve(measured_T, max_count);
    (void)normalization;
}

void render_energy_hud(void) {
    const double KEp = kinetic_energy();
    const double KEw = 0.5 * WALL_MASS * vx_wall * vx_wall;
    const double KEt = KEp + KEw;

    const double Tp  = KEp / (NUM_PARTICLES * K_B);
    const double Tw  = KEw / (0.5 * K_B);
    const double Tt  = KEt / ((NUM_PARTICLES + 0.5) * K_B);

    char line0[256], line1[256], line2[256], line3[256], line4[256], line5[256];
    snprintf(line0, sizeof(line0), "T_total: %.3f", Tt);
    snprintf(line1, sizeof(line1), "KE_tot: %.6e", KEt);
    snprintf(line2, sizeof(line2), "KE_particles: %.6e", KEp);
    snprintf(line3, sizeof(line3), "KE_wall: %.6e", KEw);
    snprintf(line4, sizeof(line4), "T_particles: %.3f", Tp);
    snprintf(line5, sizeof(line5), "T_wall: %.3f", Tw);

    SDL_Color white = {255,255,255,255};
    int x = XW1 + 1000;
    int y = YW1 + HEIGHT_UNITS * PIXELS_PER_SIGMA;

    draw_text(_renderer, font, line0, x, y, white);  y += 20;
    draw_text(_renderer, font, line1, x, y, white);  y += 20;
    draw_text(_renderer, font, line2, x, y, white);  y += 20;
    draw_text(_renderer, font, line3, x, y, white);  y += 20;
    draw_text(_renderer, font, line4, x, y, white);  y += 20;
    draw_text(_renderer, font, line5, x, y, white);
}

void render_fft_histogram_bottom(SDL_Renderer *renderer, float *fft_data, int fft_len) {
    float max_val = 1.0f;
    for (int i = 0; i < fft_len; ++i) {
        if (fft_data[i] > max_val) max_val = fft_data[i];
    }

    for (int i = 0; i < fft_len; ++i) {
        float normalized = fft_data[i] / max_val;
        int height = (int)(normalized * HIST_HEIGHT);

        SDL_SetRenderDrawColor(renderer, 255, 0, 0, 200);
        SDL_Rect rect = {
            XW1 + i * (SIM_WIDTH / fft_len),
            YW1 + SIM_HEIGHT + (HIST_HEIGHT - height),
            SIM_WIDTH / fft_len,
            height
        };
        SDL_RenderFillRect(renderer, &rect);
    }
}
