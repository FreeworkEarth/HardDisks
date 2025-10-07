/*
 * GUI Configuration Interface for Hard Disk Simulation
 * Can be compiled with Emscripten for web deployment
 */

#include <SDL2/SDL.h>
#include <SDL2/SDL_ttf.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#define GUI_WIDTH 800
#define GUI_HEIGHT 600
#define BUTTON_WIDTH 150
#define BUTTON_HEIGHT 40
#define MARGIN 20

typedef struct {
    int x, y, w, h;
    char text[64];
    bool active;
    int id;
} Button;

typedef struct {
    int x, y, w, h;
    char text[64];
    int value;
    int min_val, max_val;
    int id;
} Slider;

typedef struct {
    int x, y, w, h;
    char text[64];
    char value[32];
    int id;
    bool active;
} TextInput;

typedef struct {
    SDL_Window* window;
    SDL_Renderer* renderer;
    TTF_Font* font;
    
    // GUI Elements
    Button buttons[10];
    int num_buttons;
    
    Slider sliders[10];
    int num_sliders;
    
    TextInput inputs[10];
    int num_inputs;
    
    // Configuration
    int particles;
    float l0_units;
    float height_units;
    int num_walls;
    float wall_masses[5];
    float wall_positions[5];
    int protocol;
    int repeats;
    int steps;
    
    bool running;
    bool config_changed;
} GUIConfig;

static GUIConfig gui;

// Colors
#define COLOR_BG_R 50
#define COLOR_BG_G 50
#define COLOR_BG_B 50
#define COLOR_BG_A 255

#define COLOR_BUTTON_R 100
#define COLOR_BUTTON_G 150
#define COLOR_BUTTON_B 200
#define COLOR_BUTTON_A 255

#define COLOR_BUTTON_HOVER_R 120
#define COLOR_BUTTON_HOVER_G 170
#define COLOR_BUTTON_HOVER_B 220
#define COLOR_BUTTON_HOVER_A 255

#define COLOR_TEXT_R 255
#define COLOR_TEXT_G 255
#define COLOR_TEXT_B 255
#define COLOR_TEXT_A 255

void gui_init_ui();
void gui_generate_command();
void gui_save_config();
void gui_load_config();

void gui_init() {
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        printf("SDL could not initialize! SDL_Error: %s\n", SDL_GetError());
        return;
    }
    
    if (TTF_Init() == -1) {
        printf("TTF could not initialize! TTF_Error: %s\n", TTF_GetError());
        return;
    }
    
    gui.window = SDL_CreateWindow("Hard Disk Simulation - Configuration",
                                 SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                                 GUI_WIDTH, GUI_HEIGHT, SDL_WINDOW_SHOWN);
    if (!gui.window) {
        printf("Window could not be created! SDL_Error: %s\n", SDL_GetError());
        return;
    }
    
    gui.renderer = SDL_CreateRenderer(gui.window, -1, SDL_RENDERER_ACCELERATED);
    if (!gui.renderer) {
        printf("Renderer could not be created! SDL_Error: %s\n", SDL_GetError());
        return;
    }
    
    // Load font (same as main simulation)
    gui.font = TTF_OpenFont("Roboto-Regular.ttf", 18);
    if (!gui.font) {
        printf("Failed to load font: %s\n", TTF_GetError());
        // Try fallback fonts
        gui.font = TTF_OpenFont("/System/Library/Fonts/Arial.ttf", 18);
        if (!gui.font) {
            gui.font = TTF_OpenFont("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 18);
        }
        if (!gui.font) {
            printf("Could not load any font! TTF_Error: %s\n", TTF_GetError());
        }
    }
    
    // Initialize default values
    gui.particles = 1000;
    gui.l0_units = 20.0f;
    gui.height_units = 10.0f;
    gui.num_walls = 1;
    gui.wall_masses[0] = 200.0f;
    gui.wall_positions[0] = 0.0f;
    gui.protocol = 1; // sigmoidal
    gui.repeats = 1;
    gui.steps = 3000;
    
    gui.running = true;
    gui.config_changed = false;
    
    gui_init_ui();
}

void gui_init_ui() {
    int y = MARGIN;
    
    // Title
    gui.num_buttons = 0;
    gui.num_sliders = 0;
    gui.num_inputs = 0;
    
    // Particles slider
    gui.sliders[gui.num_sliders] = (Slider){MARGIN, y, 200, 30, "Particles", gui.particles, 100, 5000, 0};
    gui.num_sliders++;
    y += 50;
    
    // L0 units slider
    gui.sliders[gui.num_sliders] = (Slider){MARGIN, y, 200, 30, "L0 Units", (int)gui.l0_units, 5, 50, 1};
    gui.num_sliders++;
    y += 50;
    
    // Height units slider
    gui.sliders[gui.num_sliders] = (Slider){MARGIN, y, 200, 30, "Height Units", (int)gui.height_units, 5, 30, 2};
    gui.num_sliders++;
    y += 50;
    
    // Number of walls slider
    gui.sliders[gui.num_sliders] = (Slider){MARGIN, y, 200, 30, "Number of Walls", gui.num_walls, 1, 5, 3};
    gui.num_sliders++;
    y += 50;
    
    // Wall mass factor slider
    gui.sliders[gui.num_sliders] = (Slider){MARGIN, y, 200, 30, "Wall Mass Factor", (int)gui.wall_masses[0], 10, 1000, 4};
    gui.num_sliders++;
    y += 50;
    
    // Protocol buttons
    y += 20;
    gui.buttons[gui.num_buttons] = (Button){MARGIN, y, BUTTON_WIDTH, BUTTON_HEIGHT, "Step Protocol", false, 0};
    gui.num_buttons++;
    y += 50;
    
    gui.buttons[gui.num_buttons] = (Button){MARGIN, y, BUTTON_WIDTH, BUTTON_HEIGHT, "Sigmoidal Protocol", true, 1};
    gui.num_buttons++;
    y += 50;
    
    gui.buttons[gui.num_buttons] = (Button){MARGIN, y, BUTTON_WIDTH, BUTTON_HEIGHT, "Linear Protocol", false, 2};
    gui.num_buttons++;
    y += 50;
    
    // Action buttons
    y += 30;
    gui.buttons[gui.num_buttons] = (Button){MARGIN, y, BUTTON_WIDTH, BUTTON_HEIGHT, "Run Simulation", false, 10};
    gui.num_buttons++;
    
    gui.buttons[gui.num_buttons] = (Button){MARGIN + BUTTON_WIDTH + 20, y, BUTTON_WIDTH, BUTTON_HEIGHT, "Save Config", false, 11};
    gui.num_buttons++;
    
    gui.buttons[gui.num_buttons] = (Button){MARGIN + 2*(BUTTON_WIDTH + 20), y, BUTTON_WIDTH, BUTTON_HEIGHT, "Load Config", false, 12};
    gui.num_buttons++;
    
    y += 60;
    gui.buttons[gui.num_buttons] = (Button){MARGIN, y, BUTTON_WIDTH, BUTTON_HEIGHT, "Exit", false, 13};
    gui.num_buttons++;
}

void gui_draw_button(Button* btn) {
    if (btn->active) {
        SDL_SetRenderDrawColor(gui.renderer, COLOR_BUTTON_HOVER_R, COLOR_BUTTON_HOVER_G, COLOR_BUTTON_HOVER_B, COLOR_BUTTON_HOVER_A);
    } else {
        SDL_SetRenderDrawColor(gui.renderer, COLOR_BUTTON_R, COLOR_BUTTON_G, COLOR_BUTTON_B, COLOR_BUTTON_A);
    }
    
    SDL_Rect rect = {btn->x, btn->y, btn->w, btn->h};
    SDL_RenderFillRect(gui.renderer, &rect);
    
    // Draw border
    SDL_SetRenderDrawColor(gui.renderer, 255, 255, 255, 255);
    SDL_RenderDrawRect(gui.renderer, &rect);
    
    // Draw text (same approach as main simulation)
    if (gui.font) {
        SDL_Color text_color = {COLOR_TEXT_R, COLOR_TEXT_G, COLOR_TEXT_B, COLOR_TEXT_A};
        SDL_Surface* text_surface = TTF_RenderText_Blended(gui.font, btn->text, text_color);
        if (text_surface) {
            SDL_Texture* text_texture = SDL_CreateTextureFromSurface(gui.renderer, text_surface);
            if (text_texture) {
                SDL_Rect text_rect = {btn->x + 10, btn->y + 10, text_surface->w, text_surface->h};
                SDL_RenderCopy(gui.renderer, text_texture, NULL, &text_rect);
                SDL_DestroyTexture(text_texture);
            }
            SDL_FreeSurface(text_surface);
        }
    }
}

void gui_draw_slider(Slider* slider) {
    // Draw label
    if (gui.font) {
        SDL_Color text_color = {COLOR_TEXT_R, COLOR_TEXT_G, COLOR_TEXT_B, COLOR_TEXT_A};
        SDL_Surface* text_surface = TTF_RenderText_Blended(gui.font, slider->text, text_color);
        if (text_surface) {
            SDL_Texture* text_texture = SDL_CreateTextureFromSurface(gui.renderer, text_surface);
            if (text_texture) {
                SDL_Rect text_rect = {slider->x, slider->y - 20, text_surface->w, text_surface->h};
                SDL_RenderCopy(gui.renderer, text_texture, NULL, &text_rect);
                SDL_DestroyTexture(text_texture);
            }
            SDL_FreeSurface(text_surface);
        }
    }
    
    // Draw slider track
    SDL_SetRenderDrawColor(gui.renderer, 100, 100, 100, 255);
    SDL_Rect track = {slider->x, slider->y + 10, slider->w, 10};
    SDL_RenderFillRect(gui.renderer, &track);
    
    // Draw slider handle
    float ratio = (float)(slider->value - slider->min_val) / (slider->max_val - slider->min_val);
    int handle_x = slider->x + (int)(ratio * (slider->w - 20));
    
    SDL_SetRenderDrawColor(gui.renderer, 200, 200, 200, 255);
    SDL_Rect handle = {handle_x, slider->y, 20, 30};
    SDL_RenderFillRect(gui.renderer, &handle);
    
    // Draw value
    char value_str[32];
    snprintf(value_str, sizeof(value_str), "%d", slider->value);
    if (gui.font) {
        SDL_Color text_color = {COLOR_TEXT_R, COLOR_TEXT_G, COLOR_TEXT_B, COLOR_TEXT_A};
        SDL_Surface* text_surface = TTF_RenderText_Blended(gui.font, value_str, text_color);
        if (text_surface) {
            SDL_Texture* text_texture = SDL_CreateTextureFromSurface(gui.renderer, text_surface);
            if (text_texture) {
                SDL_Rect text_rect = {slider->x + slider->w + 10, slider->y, text_surface->w, text_surface->h};
                SDL_RenderCopy(gui.renderer, text_texture, NULL, &text_rect);
                SDL_DestroyTexture(text_texture);
            }
            SDL_FreeSurface(text_surface);
        }
    }
}

void gui_render() {
    SDL_SetRenderDrawColor(gui.renderer, 50, 50, 50, 255);
    SDL_RenderClear(gui.renderer);
    
    // Draw title
    if (gui.font) {
        SDL_Color text_color = {COLOR_TEXT_R, COLOR_TEXT_G, COLOR_TEXT_B, COLOR_TEXT_A};
        SDL_Surface* text_surface = TTF_RenderText_Blended(gui.font, "Hard Disk Simulation Configuration", text_color);
        if (text_surface) {
            SDL_Texture* text_texture = SDL_CreateTextureFromSurface(gui.renderer, text_surface);
            if (text_texture) {
                SDL_Rect text_rect = {MARGIN, 10, text_surface->w, text_surface->h};
                SDL_RenderCopy(gui.renderer, text_texture, NULL, &text_rect);
                SDL_DestroyTexture(text_texture);
            }
            SDL_FreeSurface(text_surface);
        }
    }
    
    // Draw sliders
    for (int i = 0; i < gui.num_sliders; i++) {
        gui_draw_slider(&gui.sliders[i]);
    }
    
    // Draw buttons
    for (int i = 0; i < gui.num_buttons; i++) {
        gui_draw_button(&gui.buttons[i]);
    }
    
    SDL_RenderPresent(gui.renderer);
}

bool gui_handle_click(int x, int y) {
    // Check button clicks
    for (int i = 0; i < gui.num_buttons; i++) {
        Button* btn = &gui.buttons[i];
        if (x >= btn->x && x <= btn->x + btn->w && y >= btn->y && y <= btn->y + btn->h) {
            if (btn->id < 10) {
                // Protocol selection
                for (int j = 0; j < gui.num_buttons; j++) {
                    if (gui.buttons[j].id < 10) {
                        gui.buttons[j].active = false;
                    }
                }
                btn->active = true;
                gui.protocol = btn->id;
            } else {
                // Action buttons
                switch (btn->id) {
                    case 10: // Run Simulation
                        gui_generate_command();
                        return true; // Exit GUI
                    case 11: // Save Config
                        gui_save_config();
                        break;
                    case 12: // Load Config
                        gui_load_config();
                        break;
                    case 13: // Exit
                        gui.running = false;
                        return true;
                }
            }
            return false;
        }
    }
    
    // Check slider clicks
    for (int i = 0; i < gui.num_sliders; i++) {
        Slider* slider = &gui.sliders[i];
        if (x >= slider->x && x <= slider->x + slider->w && y >= slider->y && y <= slider->y + slider->h) {
            float ratio = (float)(x - slider->x) / slider->w;
            slider->value = slider->min_val + (int)(ratio * (slider->max_val - slider->min_val));
            
            // Update corresponding values
            switch (slider->id) {
                case 0: gui.particles = slider->value; break;
                case 1: gui.l0_units = (float)slider->value; break;
                case 2: gui.height_units = (float)slider->value; break;
                case 3: gui.num_walls = slider->value; break;
                case 4: gui.wall_masses[0] = (float)slider->value; break;
            }
            return false;
        }
    }
    
    return false;
}

void gui_generate_command() {
    printf("Generated command:\n");
    printf("./00ALLINONE --no-experiments");
    printf(" --particles=%d", gui.particles);
    printf(" --l0=%.1f", gui.l0_units);
    printf(" --height=%.1f", gui.height_units);
    printf(" --num-walls=%d", gui.num_walls);
    printf(" --wall-mass-factors=%.1f", gui.wall_masses[0]);
    printf(" --protocol=");
    switch (gui.protocol) {
        case 0: printf("step"); break;
        case 1: printf("sigmoidal"); break;
        case 2: printf("linear"); break;
    }
    printf("\n");
}

void gui_save_config() {
    printf("Saving configuration...\n");
    // Implementation for saving config to file
}

void gui_load_config() {
    printf("Loading configuration...\n");
    // Implementation for loading config from file
}

void gui_cleanup() {
    if (gui.font) TTF_CloseFont(gui.font);
    if (gui.renderer) SDL_DestroyRenderer(gui.renderer);
    if (gui.window) SDL_DestroyWindow(gui.window);
    TTF_Quit();
    SDL_Quit();
}

void gui_run() {
    gui_init();
    
    SDL_Event e;
    while (gui.running) {
        while (SDL_PollEvent(&e)) {
            switch (e.type) {
                case SDL_QUIT:
                    gui.running = false;
                    break;
                case SDL_MOUSEBUTTONDOWN:
                    if (e.button.button == SDL_BUTTON_LEFT) {
                        if (gui_handle_click(e.button.x, e.button.y)) {
                            gui_cleanup();
                            return;
                        }
                    }
                    break;
            }
        }
        
        gui_render();
        SDL_Delay(16); // ~60 FPS
    }
    
    gui_cleanup();
}

int main(int argc, char* argv[]) {
    printf("Hard Disk Simulation - Configuration GUI\n");
    printf("This GUI can be compiled with Emscripten for web deployment\n\n");
    
    gui_run();
    return 0;
}
