// keyboard (and mouse) stuff
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <SDL2/SDL.h>
#include "glgraph.h"
#include "rstuff.h"
#include "globals.h"


//#define REAL float
//#define XFAC (REAL)SIZEX/SIZEY
// big box outline...
//#define XW1 -1.0*XFAC

extern int firstex;			// signal to videoUpdateFrame()
extern int gl_keys[512]; // for recording key presses


// Add SDL2 initialization and event handling code here
void handleEvents() {
    SDL_Event event;
    while (SDL_PollEvent(&event)) {
        switch (event.type) {
            case SDL_QUIT:
                exit(0);
                break;
            case SDL_KEYDOWN:
                // Handle key down events
                break;
            case SDL_KEYUP:
                // Handle key up events
                break;
            case SDL_MOUSEBUTTONDOWN:
                // Handle mouse button down events
                break;
            case SDL_MOUSEBUTTONUP:
                // Handle mouse button up events
                break;
            case SDL_MOUSEMOTION:
                // Handle mouse motion events
                break;
            default:
                break;
        }
    }
}




#define SDL_SCANCODE_LEFTBRACKET SDL_SCANCODE_LEFTBRACKET
#define SDL_SCANCODE_RIGHTBRACKET SDL_SCANCODE_RIGHTBRACKET
#define SDL_SCANCODE_COMMA SDL_SCANCODE_COMMA
#define SDL_SCANCODE_PERIOD SDL_SCANCODE_PERIOD



void keyCB(SDL_Event event) {
    if (event.type == SDL_KEYDOWN) {
        SDL_Keycode key = event.key.keysym.sym;
        SDL_Keymod mod = SDL_GetModState();

        if (mod & KMOD_SHIFT) gl_shift = 1;
        if (mod & KMOD_CTRL) gl_ctl = 1;
        if (mod & KMOD_ALT) gl_opt = 1;

        gl_keys[key] = 1; // mark key pressed

        switch (key) {
            long i;
            case SDLK_ESCAPE: // escape
                printf("\nrob shaw  may 2011\n");
                printf("rob@haptek.com\n\n");
                exit(0);
            case 'b': // for adjusting radii
                radctl ^= 1;
                printf("radctl %d\n", radctl);
                break;
            case 'B':
                dispctl ^= 1;
                printf("dispctl %d\n", dispctl);
                break;
            case 'c': // reverse all velocities!
                for (i = 0; i < nspheres; i++) {
                    VX[i] = -VX[i];
                    VY[i] = -VY[i];
                }
                vpist1 = -vpist1;
                printf("reverse velocities!\n");
                break;
            case 'C': // set square boundaries
                break;
            case 'd': // change display mode
                dispctl += 1;
                if (dispctl > 3) dispctl = 0;
                printf("dispctl %d\n", dispctl);
                break;
            case 'D': // change second display mode
                Dispctl += 1;
                if (Dispctl > 3) Dispctl = 0;
                printf("Dispctl %d\n", Dispctl);
                break;
            case 'e': // reset disk energies to one
                reseteng();
                break;
            case 'E': // double all velocities
                for (i = 0; i < nspheres; i++) {
                    VX[i] *= 2;
                    VY[i] *= 2;
                }
                break;
            case 'f':
                gl_fullscreen ^= 1;
                if (gl_fullscreen) {
                    SDL_SetWindowFullscreen(_window, SDL_WINDOW_FULLSCREEN);
                } else {
                    SDL_SetWindowFullscreen(_window, 0);
                    SDL_SetWindowPosition(_window, gl_left, gl_top);
                    SDL_SetWindowSize(_window, gl_ix, gl_iy);
                }
                break;
            case 'r':
                gl_flip ^= 1;
                break;
            case 'o': {
                float dt;
                static long oldncount = 0;
                extern long ncount;
                gettimeofday(&gl_tm1, NULL);
                dt = gl_tm1.tv_sec - gl_tm.tv_sec + .000001 * (gl_tm1.tv_usec - gl_tm.tv_usec);
                gl_tm.tv_sec = gl_tm1.tv_sec;
                gl_tm.tv_usec = gl_tm1.tv_usec;
                printf("elapsed time %.2f\t", dt);
                printf("frames per second %.2f\n", gl_framecount / dt);
                gl_framecount = 0;
                printf("iterations per second %.2f\n", (ncount - oldncount) / dt);
                oldncount = ncount;
                break;
            }
            case 'h': // turn on histogramming
                ++histctl;
                if (histctl > 3) histctl = 0;
                printf("histctl %ld\n", histctl);
                break;
            case 'H': // clear histogram
                for (i = 0; i < NENTS; i++) Histo1[i] = Histowall[i] = HistoW[i] = 0;
                break;
            case 'i': // toggle black/white
                gl_bwctl ^= 1;
                if (gl_bwctl) {
                    g_backcolor(BLACK);
                    g_color(WHITE);
                } else {
                    g_backcolor(WHITE);
                    g_color(BLACK);
                }
                break;
            case 'p': // print parameters
                gl_print ^= 1;
                if (gl_print == 1) {
                    printf("xpos %.2f  ypos %.2f  zpos %.2f\n", gl_xpos, gl_ypos, gl_zpos);
                    printf("xrot %.2f  yrot %.2f  zrot %.2f\n", gl_xrot, gl_yrot, gl_zrot);
                    printf("printing on\n");
                } else
                    printf("printing off\n");
                break;
            case ' ': // start off pistons...
                pistctl = 1;
                ncount = 1;
                if (pistctl) printf("right piston\n");
                else printf("left piston\n");
                break;
            case '/': // print help
                printhelp();
                break;
            case 'N': // reset
                ncount = 0;
                break;
            case 'R':
                gl_xpos = gl_ypos = gl_zpos = 0;
                gl_xrot = gl_yrot = gl_zrot = 0;
                break;
            case 'a':
                Slowctl *= 10;
                if (Slowctl > 1000) Slowctl = 1;
                printf("Slowctl %ld\n", Slowctl);
                break;
            case 'A':
                Slowctl /= 10;
                if (Slowctl == 0) Slowctl = 1;
                printf("Slowctl %ld\n", Slowctl);
                break;
            case 's':
                if (gl_opt) keyslowctl /= 2;
                else keyslowctl *= 10;
                if (keyslowctl > 10000) keyslowctl = 1;
                if (keyslowctl == 0) keyslowctl = 1;
                printf("slowctl %ld\n", keyslowctl);
                break;
            case 19: // ctl s
                keyslowctl *= 2;
                printf("slowctl %ld\n", keyslowctl);
                break;
            case 'S':
                keyslowctl /= 10;
                if (keyslowctl == 0) keyslowctl = 1;
                printf("slowctl %ld\n", keyslowctl);
                break;
            case 't': // trail control
                ++trailctl;
                color(WHITE);
                clear();
                if (trailctl > 3) trailctl = 0;
                printf("trailctl %ld\n", trailctl);
                break;
            case 'T': // erase trails
                color(WHITE);
                clear();
                break;
            case 'n': // reset pistons
                if (gl_shift) {
                    printf("tweak one!\n");
                }
                ncount = 0;
                break;
            case 'V':
                gl_sync ^= 1;
                break;
        }
    }
}




void specialKeyPressed(SDL_Keycode key) {
    REAL rad;

    gl_keys[key + 256] = 1; // mark key pressed
    const Uint8 *state = SDL_GetKeyboardState(NULL);
    gl_shift = state[SDL_SCANCODE_LSHIFT] || state[SDL_SCANCODE_RSHIFT];
    gl_ctl = state[SDL_SCANCODE_LCTRL] || state[SDL_SCANCODE_RCTRL];
    gl_opt = state[SDL_SCANCODE_LALT] || state[SDL_SCANCODE_RALT];

    if (gl_keys[SDLK_UP + 256]) { // up arrow
        rad = gl_shift ? rad2 : rad1;
        if (gl_opt)
            addisks1(10, rad);
        else if (gl_ctl)
            addisks1(100, rad);
        else
            addisks1(1, rad);
    }

    if (gl_keys[SDLK_DOWN + 256]) { // down arrow
        rad = gl_shift ? rad2 : rad1;
        if (gl_opt)
            subdisks1(10, rad);
        else if (gl_ctl)
            subdisks1(100, rad);
        else
            subdisks1(1, rad);
    }

    if (gl_keys[SDLK_LEFT + 256]) { // left arrow
        if (!gl_opt && !radctl) { // not a display rotation or radius change
            if (gl_ctl) { // move right-hand piston
                vpist2 = gl_shift ? -5 * dxpist / (slowctl * DELT) : -dxpist / (slowctl * DELT);
            } else { // move left-hand piston
                vpist1 = gl_shift ? -5 * dxpist / (slowctl * DELT) : -dxpist / (slowctl * DELT);
            }
        }
        if (!gl_opt && radctl) { // a radius change, not a display rotation
            int i;
            if (gl_shift) {
                rad1 -= .01;
                for (i = 0; i < nspheres; i++) Radius[i] -= .01;
            } else {
                rad1 -= .001;
                for (i = 0; i < nspheres; i++) Radius[i] -= .001;
            }
            printf("rad1 %f\n", rad1);
        }
    }

    if (gl_keys[SDLK_RIGHT + 256]) { // right arrow
        if (!gl_opt && !radctl) { // not a display rotation or radius change
            if (gl_ctl) { // move right-hand piston
                vpist2 = gl_shift ? 5 * dxpist / (slowctl * DELT) : dxpist / (slowctl * DELT);
            } else { // move left-hand piston
                vpist1 = gl_shift ? 5 * dxpist / (slowctl * DELT) : dxpist / (slowctl * DELT);
            }
        }
        if (!gl_opt && radctl) { // a radius change, not a display rotation
            int i;
            if (gl_shift) {
                rad1 += .01;
                for (i = 0; i < nspheres; i++) Radius[i] += .01;
            } else {
                rad1 += .001;
                for (i = 0; i < nspheres; i++) Radius[i] += .001;
            }
            printf("rad1 %f\n", rad1);
        }
    }
}




void specialKeyUp(SDL_Keycode key) {
    gl_keys[key + 256] = 0; // mark key released
    gl_shift = gl_ctl = gl_opt = 0; // turn these on and off every keystroke
    if (key == SDLK_LEFT) { // left arrow
        vpist1 = 0;
        vpist2 = 0;
    }
    if (key == SDLK_RIGHT) { // right arrow
        vpist1 = 0;
        vpist2 = 0;
    }
}








void clickCB(int button, int state, int x, int y) {
    printf("%d %d %d %d\n", button, state, x, y);
    if (!state)
        SDL_SetCursor(SDL_CreateSystemCursor(SDL_SYSTEM_CURSOR_HAND));
    else
        SDL_SetCursor(SDL_CreateSystemCursor(SDL_SYSTEM_CURSOR_NO));

    gl_imx = x; gl_imy = y;
    gl_butnum = button; gl_butstate = state;
    gl_mousex = (float)x / (float)gl_xsize;
    gl_mousey = (float)y / (float)gl_ysize;
    printf("%f %f\n", gl_mousex, gl_mousey);
}

void mouseCB(int x, int y) {
    printf("%d %d\n", x, y);
}

void updatekeys() {
    if (gl_opt) {
        if (gl_keys[SDL_SCANCODE_LEFT]) // left arrow
            gl_yrot -= 1 * gl_delt;
        if (gl_keys[SDL_SCANCODE_RIGHT]) // right arrow
            gl_yrot += 1 * gl_delt;
        if (gl_keys[SDL_SCANCODE_UP]) // up arrow
            gl_zpos += 10 * gl_delt;
        if (gl_keys[SDL_SCANCODE_DOWN]) // down arrow
            gl_zpos -= 10 * gl_delt;
        if (gl_keys[SDL_SCANCODE_LEFTBRACKET]) // left
            gl_xpos -= 10 * gl_delt;
        if (gl_keys[SDL_SCANCODE_RIGHTBRACKET]) // right
            gl_xpos += 10 * gl_delt;
        if (gl_keys[SDL_SCANCODE_COMMA]) // down
            gl_ypos += 10 * gl_delt;
        if (gl_keys[SDL_SCANCODE_PERIOD]) // up
            gl_ypos -= 10 * gl_delt;
        if (gl_keys[SDL_SCANCODE_COMMA]) // x axis rotations
            gl_xrot -= 1 * gl_delt;
        if (gl_keys[SDL_SCANCODE_PERIOD])
            gl_xrot += 1 * gl_delt;
        if (gl_keys[SDL_SCANCODE_COMMA]) // z axis rotations
            gl_zrot -= 1 * gl_delt;
        if (gl_keys[SDL_SCANCODE_PERIOD])
            gl_zrot += 1 * gl_delt;
    
    } else {
        if (gl_keys[SDL_SCANCODE_LEFT]) { // left arrow
            if (gl_ctl) { // right piston
                if (xpist2 < XW1) {
                    xpist2 = XW1;
                    vpist2 = 0;
                }
            } else { // left piston
                if (xpist1 < XW1) {
                    xpist1 = XW1;
                    vpist1 = 0;
                }
            }
        }
        if (gl_keys[SDL_SCANCODE_RIGHT]) { // right arrow
            if (gl_ctl) { // right piston
                if (xpist2 > XW2) {
                    xpist2 = XW2;
                    vpist2 = 0;
                }
            } else { // left piston
                if (xpist1 > XW2) {
                    xpist1 = XW2;
                    vpist1 = 0;
                }
            }
        }
    }
}

void printhelp() {
    printf("\nleft, right arrow keys - move left piston (shift increases speed)\n");
    printf("control left, right arrow keys - move right piston (shift increases speed)\n");
    printf("n - reset piston to left side\n");
    printf("option n - reset piston to center\n");
    printf("N - reset disks to right side\n");
    printf("up, down arrow keys - add, subtract disks\n");
    printf("\tshift for smaller disks, ctl for 100 disks at a time\n");
    printf("h - change histogrammer mode\n");
    printf("H - erase history histogram\n");
    printf("=,- - increase, decrease disk radius 1 (opt for greater change)\n");
    printf("+,_ - increase, decrease disk radius 2 (opt for greater change)\n");
    printf("t - turn on trails\n");
    printf("space - erase trails\n");
    printf("d - toggle display\n");
    printf("s, S - change time steps per display draw\n");
    printf("a, A - change screen displays per display draw\n");
    printf("c - reverse velocities\n");
    printf("\n3, 4 - decrease, increase timer interval\n");
    printf("shift 3, 4 - steps of 10 milliseconds\n");
    printf("option up, down arrow keys - zoom\n");
    printf("option left, right arrow keys - y axis rotation\n");
    printf("option [ , ] keys - move left, right\n");
    printf("option { , } keys - move down, up\n");
    printf("option , , . keys - x axis rotation\n");
    printf("option < , > keys - z axis rotation\n");
    printf("f - toggle full screen\n");
    printf("i - toggle black, white inversion\n");
    printf("r - toggle mirror reflection\n");
    printf("p - print parameters\n");
    printf("R - reset position\n");
    printf("o - print frames per second\n");
    printf("/ - print help\n");
    printf("<esc> - exit\n\n");
}