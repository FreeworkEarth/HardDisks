#ifndef CORE_IO_LOG_H
#define CORE_IO_LOG_H

#include <stdio.h>
#include "state.h"
#include "step_core.h"

#ifdef __cplusplus
extern "C" {
#endif

void log_energy(FILE *fptr);
void log_packing_fractions(FILE *log, const char *label);

#ifdef __cplusplus
}
#endif

#endif /* CORE_IO_LOG_H */
