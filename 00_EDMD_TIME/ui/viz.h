#ifndef UI_VIZ_H
#define UI_VIZ_H

#include <stdbool.h>
#include "state.h"

#ifdef __cplusplus
extern "C" {
#endif

bool viz_init(int width, int height);
void viz_shutdown(void);

void draw_clear_screen(void);
void draw_coordinate_system(void);
void draw_simulation_boundary(void);
void render_particles(void);
void render_pistons(void);
void draw_wall(void);
void render_velocity_histograms(void);
void render_energy_hud(void);
void render_fft_histogram_bottom(SDL_Renderer *renderer, float *fft_data, int fft_len);
void draw_text(SDL_Renderer *renderer, TTF_Font *font, const char *text, int x, int y, SDL_Color color);

#ifdef __cplusplus
}
#endif

#endif /* UI_VIZ_H */
