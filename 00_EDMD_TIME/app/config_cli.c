#include "config_cli.h"
#include "state.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void config_cli_parse(int argc, char **argv) {
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--no-wall-hold") == 0) {
            wall_hold_enabled = false;
        } else if (strcmp(argv[i], "--wall") == 0) {
            wall_enabled = 1;
        } else if (strcmp(argv[i], "--speed-of-sound") == 0) {
            enable_speed_of_sound_experiments = true;
        } else if (strcmp(argv[i], "--help") == 0) {
            printf("Options:\n");
            printf("  --wall             Enable the movable wall at startup\n");
            printf("  --no-wall-hold     Release the wall immediately\n");
            printf("  --speed-of-sound   Run the speed of sound experiment batch\n");
            exit(0);
        }
    }
}
