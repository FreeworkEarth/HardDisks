#ifndef EXPERIMENT_VALIDATION_H
#define EXPERIMENT_VALIDATION_H

#include <stddef.h>

typedef struct {
    int particle_count;
    const double *x;
    const double *y;
    const double *vx;
    const double *vy;
    double particle_radius;

    double x_min_face;
    double x_max_face;
    double y_min;
    double y_max;

    int wall_count;
    const double *wall_x;
    const double *wall_vx;
    double wall_half_thickness;

    double kinetic_energy;
    double spring_energy;
    double piston_work_left;
    double piston_work_right;
    long forced_advance_events;
} ExperimentValidationState;

typedef struct {
    int failed;
    char phase[32];
    char reason[64];
    char detail[512];
    int step;
    double time_sigma;
    int particle_a;
    int particle_b;
    int wall_a;
    int wall_b;
    double value;
    double limit;
} ExperimentValidationFailure;

typedef struct {
    int initial_particle_count;
    int *initial_particle_segments;
    int *cell_head;
    int *cell_next;
    int cell_capacity;
    int next_capacity;
    ExperimentValidationFailure failure;
} ExperimentValidator;

void experiment_validator_init(ExperimentValidator *validator);
void experiment_validator_destroy(ExperimentValidator *validator);

/* Snapshot each particle's compartment. Solid-wall experiments can later use
   check_particle_identity=1 to catch crossings even when two opposite
   crossings leave the aggregate compartment counts unchanged. */
int experiment_validator_snapshot(ExperimentValidator *validator,
                                  const ExperimentValidationState *state,
                                  const char *phase,
                                  int step,
                                  double time_sigma);

/* Returns 1 while valid and 0 after the first invariant violation. */
int experiment_validator_check(ExperimentValidator *validator,
                               const ExperimentValidationState *state,
                               const char *phase,
                               int step,
                               double time_sigma,
                               int check_particle_identity);

const ExperimentValidationFailure *experiment_validator_failure(
    const ExperimentValidator *validator);

/* Record an experiment-specific invariant using the same first-failure
   semantics and failure ledger fields as the generic checks. */
int experiment_validator_record_failure(ExperimentValidator *validator,
                                        const char *phase,
                                        const char *reason,
                                        const char *detail,
                                        int step,
                                        double time_sigma,
                                        int particle_a,
                                        int particle_b,
                                        int wall_a,
                                        int wall_b,
                                        double value,
                                        double limit);

#endif
