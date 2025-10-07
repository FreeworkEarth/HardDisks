#ifndef CORE_PROTOCOL_H
#define CORE_PROTOCOL_H

#include "state.h"

#ifdef __cplusplus
extern "C" {
#endif

void move_left_piston(float dt);
void move_right_piston(float dt);
void update_pistons(float dt);
void update_wall(float dt, float L0_units);

#ifdef __cplusplus
}
#endif

#endif /* CORE_PROTOCOL_H */
