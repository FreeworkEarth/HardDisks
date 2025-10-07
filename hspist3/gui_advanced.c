/*
 * Advanced GUI Configuration Interface for Hard Disk Simulation
 * Matches Python GUI functionality with tabs and comprehensive controls
 * Can be compiled with Emscripten for web deployment
 */

#include <SDL2/SDL.h>
#include <SDL2/SDL_ttf.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#define GUI_WIDTH 1000
#define GUI_HEIGHT 700
#define BUTTON_WIDTH 150
#define BUTTON_HEIGHT 40
#define MARGIN 20
#define TAB_HEIGHT 40

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
    int x, y, w, h;
    char text[128];
    bool active;
    int id;
} Tab;

typedef struct {
    SDL_Window* window;
    SDL_Renderer* renderer;
    TTF_Font* font;
    TTF_Font* font_large;
    
    // GUI Elements
    Button buttons[20];
    int num_buttons;
    
    Slider sliders[20];
    int num_sliders;
    
    TextInput inputs[20];
    int num_inputs;
    
    Tab tabs[10];
    int num_tabs;
    int active_tab;
    
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
    int experiment_type;
    
    bool running;
    bool config_changed;
} GUIConfig;

static GUIConfig gui;

// Colors
#define COLOR_BG_R 30
#define COLOR_BG_G 30
#define COLOR_BG_B 30
#define COLOR_BG_A 255

#define COLOR_TAB_R 60
#define COLOR_TAB_G 60
#define COLOR_TAB_B 60
#define COLOR_TAB_A 255

#define COLOR_TAB_ACTIVE_R 100
#define COLOR_TAB_ACTIVE_G 150
#define COLOR_TAB_ACTIVE_B 200
#define COLOR_TAB_ACTIVE_A 255

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
    
    gui.window = SDL_CreateWindow("Hard Disk Simulation - Advanced Configuration",
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
    
    // Load fonts (same as main simulation)
    gui.font = TTF_OpenFont("Roboto-Regular.ttf", 16);
    if (!gui.font) {
        gui.font = TTF_OpenFont("/System/Library/Fonts/Arial.ttf", 16);
    }
    if (!gui.font) {
        gui.font = TTF_OpenFont("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 16);
    }
    
    gui.font_large = TTF_OpenFont("Roboto-Regular.ttf", 20);
    if (!gui.font_large) {
        gui.font_large = TTF_OpenFont("/System/Library/Fonts/Arial.ttf", 20);
    }
    if (!gui.font_large) {
        gui.font_large = gui.font;
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
    gui.experiment_type = 0; // none
    
    gui.running = true;
    gui.config_changed = false;
    gui.active_tab = 0;
    
    gui_init_ui();
}

void gui_init_ui() {
    // Initialize tabs
    gui.num_tabs = 0;
    gui.tabs[gui.num_tabs] = (Tab){MARGIN, MARGIN, 120, TAB_HEIGHT, "Experiment", true, 0};
    gui.num_tabs++;
    gui.tabs[gui.num_tabs] = (Tab){MARGIN + 130, MARGIN, 120, TAB_HEIGHT, "Walls", false, 1};
    gui.num_tabs++;
    gui.tabs[gui.num_tabs] = (Tab){MARGIN + 260, MARGIN, 120, TAB_HEIGHT, "Protocol", false, 2};
    gui.num_tabs++;
    gui.tabs[gui.num_tabs] = (Tab){MARGIN + 390, MARGIN, 120, TAB_HEIGHT, "Simulation", false, 3};
    gui.num_tabs++;
    
    // Initialize other elements
    gui.num_buttons = 0;
    gui.num_sliders = 0;
    gui.num_inputs = 0;
}

void gui_draw_tab(Tab* tab) {
    if (tab->active) {
        SDL_SetRenderDrawColor(gui.renderer, COLOR_TAB_ACTIVE_R, COLOR_TAB_ACTIVE_G, COLOR_TAB_ACTIVE_B, COLOR_TAB_ACTIVE_A);
    } else {
        SDL_SetRenderDrawColor(gui.renderer, COLOR_TAB_R, COLOR_TAB_G, COLOR_TAB_B, COLOR_TAB_A);
    }
    
    SDL_Rect rect = {tab->x, tab->y, tab->w, tab->h};
    SDL_RenderFillRect(gui.renderer, &rect);
    
    // Draw border
    SDL_SetRenderDrawColor(gui.renderer, 255, 255, 255, 255);
    SDL_RenderDrawRect(gui.renderer, &rect);
    
    // Draw text
    if (gui.font) {
        SDL_Color text_color = {COLOR_TEXT_R, COLOR_TEXT_G, COLOR_TEXT_B, COLOR_TEXT_A};
        SDL_Surface* text_surface = TTF_RenderText_Blended(gui.font, tab->text, text_color);
        if (text_surface) {
            SDL_Texture* text_texture = SDL_CreateTextureFromSurface(gui.renderer, text_surface);
            if (text_texture) {
                SDL_Rect text_rect = {tab->x + 10, tab->y + 10, text_surface->w, text_surface->h};
                SDL_RenderCopy(gui.renderer, text_texture, NULL, &text_rect);
                SDL_DestroyTexture(text_texture);
            }
            SDL_FreeSurface(text_surface);
        }
    }
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
    
    // Draw text
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

void gui_render_tab_content() {
    int y = MARGIN + TAB_HEIGHT + 20;
    
    switch (gui.active_tab) {
        case 0: // Experiment
            gui.num_buttons = 0;
            gui.num_sliders = 0;
            
            // Experiment type buttons
            gui.buttons[gui.num_buttons] = (Button){MARGIN, y, BUTTON_WIDTH, BUTTON_HEIGHT, "None", gui.experiment_type == 0, 0};
            gui.num_buttons++;
            y += 50;
            
            gui.buttons[gui.num_buttons] = (Button){MARGIN, y, BUTTON_WIDTH, BUTTON_HEIGHT, "Speed of Sound", gui.experiment_type == 1, 1};
            gui.num_buttons++;
            y += 50;
            
            gui.buttons[gui.num_buttons] = (Button){MARGIN, y, BUTTON_WIDTH, BUTTON_HEIGHT, "Energy Transfer", gui.experiment_type == 2, 2};
            gui.num_buttons++;
            y += 50;
            
            gui.buttons[gui.num_buttons] = (Button){MARGIN, y, BUTTON_WIDTH, BUTTON_HEIGHT, "Multi-Wall", gui.experiment_type == 3, 3};
            gui.num_buttons++;
            y += 50;
            
            gui.buttons[gui.num_buttons] = (Button){MARGIN, y, BUTTON_WIDTH, BUTTON_HEIGHT, "Szilard Engine", gui.experiment_type == 4, 4};
            gui.num_buttons++;
            break;
            
        case 1: // Walls
            gui.num_buttons = 0;
            gui.num_sliders = 0;
            
            // Number of walls slider
            gui.sliders[gui.num_sliders] = (Slider){MARGIN, y, 200, 30, "Number of Walls", gui.num_walls, 1, 5, 0};
            gui.num_sliders++;
            y += 50;
            
            // Wall mass factor slider
            gui.sliders[gui.num_sliders] = (Slider){MARGIN, y, 200, 30, "Wall Mass Factor", (int)gui.wall_masses[0], 10, 1000, 1};
            gui.num_sliders++;
            y += 50;
            
            // Wall position input (simplified for now)
            gui.sliders[gui.num_sliders] = (Slider){MARGIN, y, 200, 30, "Wall Position", (int)(gui.wall_positions[0] * 100), 0, 100, 2};
            gui.num_sliders++;
            break;
            
        case 2: // Protocol
            gui.num_buttons = 0;
            gui.num_sliders = 0;
            
            // Protocol buttons
            gui.buttons[gui.num_buttons] = (Button){MARGIN, y, BUTTON_WIDTH, BUTTON_HEIGHT, "Step", gui.protocol == 0, 10};
            gui.num_buttons++;
            y += 50;
            
            gui.buttons[gui.num_buttons] = (Button){MARGIN, y, BUTTON_WIDTH, BUTTON_HEIGHT, "Sigmoidal", gui.protocol == 1, 11};
            gui.num_buttons++;
            y += 50;
            
            gui.buttons[gui.num_buttons] = (Button){MARGIN, y, BUTTON_WIDTH, BUTTON_HEIGHT, "Linear", gui.protocol == 2, 12};
            gui.num_buttons++;
            y += 50;
            
            gui.buttons[gui.num_buttons] = (Button){MARGIN, y, BUTTON_WIDTH, BUTTON_HEIGHT, "Sinusoidal", gui.protocol == 3, 13};
            gui.num_buttons++;
            break;
            
        case 3: // Simulation
            gui.num_buttons = 0;
            gui.num_sliders = 0;
            
            // Simulation parameters
            gui.sliders[gui.num_sliders] = (Slider){MARGIN, y, 200, 30, "Particles", gui.particles, 100, 5000, 3};
            gui.num_sliders++;
            y += 50;
            
            gui.sliders[gui.num_sliders] = (Slider){MARGIN, y, 200, 30, "L0 Units", (int)gui.l0_units, 5, 50, 4};
            gui.num_sliders++;
            y += 50;
            
            gui.sliders[gui.num_sliders] = (Slider){MARGIN, y, 200, 30, "Height Units", (int)gui.height_units, 5, 30, 5};
            gui.num_sliders++;
            y += 50;
            
            gui.sliders[gui.num_sliders] = (Slider){MARGIN, y, 200, 30, "Steps", gui.steps, 1000, 10000, 6};
            gui.num_sliders++;
            break;
    }
    
    // Action buttons (always visible)
    y = GUI_HEIGHT - 100;
    gui.buttons[gui.num_buttons] = (Button){MARGIN, y, BUTTON_WIDTH, BUTTON_HEIGHT, "Run Simulation", false, 20};
    gui.num_buttons++;
    
    gui.buttons[gui.num_buttons] = (Button){MARGIN + BUTTON_WIDTH + 20, y, BUTTON_WIDTH, BUTTON_HEIGHT, "Save Config", false, 21};
    gui.num_buttons++;
    
    gui.buttons[gui.num_buttons] = (Button){MARGIN + 2*(BUTTON_WIDTH + 20), y, BUTTON_WIDTH, BUTTON_HEIGHT, "Load Config", false, 22};
    gui.num_buttons++;
    
    gui.buttons[gui.num_buttons] = (Button){MARGIN + 3*(BUTTON_WIDTH + 20), y, BUTTON_WIDTH, BUTTON_HEIGHT, "Exit", false, 23};
    gui.num_buttons++;
}

void gui_render() {
    SDL_SetRenderDrawColor(gui.renderer, COLOR_BG_R, COLOR_BG_G, COLOR_BG_B, COLOR_BG_A);
    SDL_RenderClear(gui.renderer);
    
    // Draw title
    if (gui.font_large) {
        SDL_Color text_color = {COLOR_TEXT_R, COLOR_TEXT_G, COLOR_TEXT_B, COLOR_TEXT_A};
        SDL_Surface* text_surface = TTF_RenderText_Blended(gui.font_large, "Hard Disk Simulation - Advanced Configuration", text_color);
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
    
    // Draw tabs
    for (int i = 0; i < gui.num_tabs; i++) {
        gui_draw_tab(&gui.tabs[i]);
    }
    
    // Draw tab content
    gui_render_tab_content();
    
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
    // Check tab clicks
    for (int i = 0; i < gui.num_tabs; i++) {
        Tab* tab = &gui.tabs[i];
        if (x >= tab->x && x <= tab->x + tab->w && y >= tab->y && y <= tab->y + tab->h) {
            // Deactivate all tabs
            for (int j = 0; j < gui.num_tabs; j++) {
                gui.tabs[j].active = false;
            }
            // Activate clicked tab
            tab->active = true;
            gui.active_tab = tab->id;
            return false;
        }
    }
    
    // Check button clicks
    for (int i = 0; i < gui.num_buttons; i++) {
        Button* btn = &gui.buttons[i];
        if (x >= btn->x && x <= btn->x + btn->w && y >= btn->y && y <= btn->y + btn->h) {
            if (btn->id < 10) {
                // Experiment selection
                for (int j = 0; j < gui.num_buttons; j++) {
                    if (gui.buttons[j].id < 10) {
                        gui.buttons[j].active = false;
                    }
                }
                btn->active = true;
                gui.experiment_type = btn->id;
            } else if (btn->id < 20) {
                // Protocol selection
                for (int j = 0; j < gui.num_buttons; j++) {
                    if (gui.buttons[j].id >= 10 && gui.buttons[j].id < 20) {
                        gui.buttons[j].active = false;
                    }
                }
                btn->active = true;
                gui.protocol = btn->id - 10;
            } else {
                // Action buttons
                switch (btn->id) {
                    case 20: // Run Simulation
                        gui_generate_command();
                        return true; // Exit GUI
                    case 21: // Save Config
                        gui_save_config();
                        break;
                    case 22: // Load Config
                        gui_load_config();
                        break;
                    case 23: // Exit
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
                case 0: gui.num_walls = slider->value; break;
                case 1: gui.wall_masses[0] = (float)slider->value; break;
                case 2: gui.wall_positions[0] = (float)slider->value / 100.0f; break;
                case 3: gui.particles = slider->value; break;
                case 4: gui.l0_units = (float)slider->value; break;
                case 5: gui.height_units = (float)slider->value; break;
                case 6: gui.steps = slider->value; break;
            }
            return false;
        }
    }
    
    return false;
}

void gui_generate_command() {
    printf("Generated command:\n");
    printf("./00ALLINONE");
    
    if (gui.experiment_type == 0) {
        printf(" --no-experiments");
    } else {
        const char* experiments[] = {"", "speed_of_sound", "energy_transfer", "multi_wall", "szilard_engine"};
        printf(" --experiment=%s", experiments[gui.experiment_type]);
    }
    
    printf(" --particles=%d", gui.particles);
    printf(" --l0=%.1f", gui.l0_units);
    printf(" --height=%.1f", gui.height_units);
    printf(" --num-walls=%d", gui.num_walls);
    printf(" --wall-mass-factors=%.1f", gui.wall_masses[0]);
    printf(" --wall-positions=%.2f", gui.wall_positions[0]);
    
    const char* protocols[] = {"step", "sigmoidal", "linear", "sinusoidal"};
    printf(" --protocol=%s", protocols[gui.protocol]);
    
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
    if (gui.font_large) TTF_CloseFont(gui.font_large);
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
    printf("Hard Disk Simulation - Advanced Configuration GUI\n");
    printf("This GUI matches the Python GUI functionality\n\n");
    
    gui_run();
    return 0;
}
