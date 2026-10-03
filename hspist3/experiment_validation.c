#include "experiment_validation.h"

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int segment_for_x(const ExperimentValidationState *state, double x) {
    int segment = 0;
    for (int w = 0; w < state->wall_count; ++w) {
        if (x > state->wall_x[w]) segment++;
        else break;
    }
    return segment;
}

static int fail(ExperimentValidator *validator,
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
                double limit) {
    ExperimentValidationFailure *f = &validator->failure;
    if (f->failed) return 0;
    f->failed = 1;
    snprintf(f->phase, sizeof(f->phase), "%s", phase ? phase : "unknown");
    snprintf(f->reason, sizeof(f->reason), "%s", reason ? reason : "unknown");
    snprintf(f->detail, sizeof(f->detail), "%s", detail ? detail : "");
    f->step = step;
    f->time_sigma = time_sigma;
    f->particle_a = particle_a;
    f->particle_b = particle_b;
    f->wall_a = wall_a;
    f->wall_b = wall_b;
    f->value = value;
    f->limit = limit;
    return 0;
}

void experiment_validator_init(ExperimentValidator *validator) {
    if (!validator) return;
    memset(validator, 0, sizeof(*validator));
    validator->failure.particle_a = -1;
    validator->failure.particle_b = -1;
    validator->failure.wall_a = -1;
    validator->failure.wall_b = -1;
    validator->failure.value = NAN;
    validator->failure.limit = NAN;
}

void experiment_validator_destroy(ExperimentValidator *validator) {
    if (!validator) return;
    free(validator->initial_particle_segments);
    free(validator->cell_head);
    free(validator->cell_next);
    validator->initial_particle_segments = NULL;
    validator->cell_head = NULL;
    validator->cell_next = NULL;
    validator->initial_particle_count = 0;
    validator->cell_capacity = 0;
    validator->next_capacity = 0;
}

static int state_shape_is_valid(ExperimentValidator *validator,
                                const ExperimentValidationState *state,
                                const char *phase,
                                int step,
                                double time_sigma) {
    if (!state || state->particle_count < 0 ||
        (state->particle_count > 0 &&
         (!state->x || !state->y || !state->vx || !state->vy)) ||
        state->wall_count < 0 ||
        (state->wall_count > 0 && (!state->wall_x || !state->wall_vx))) {
        return fail(validator, phase, "invalid_audit_input",
                    "Experiment validator received missing arrays or negative counts.",
                    step, time_sigma, -1, -1, -1, -1, NAN, NAN);
    }
    return 1;
}

int experiment_validator_snapshot(ExperimentValidator *validator,
                                  const ExperimentValidationState *state,
                                  const char *phase,
                                  int step,
                                  double time_sigma) {
    if (!validator || validator->failure.failed) return 0;
    if (!state_shape_is_valid(validator, state, phase, step, time_sigma)) return 0;
    free(validator->initial_particle_segments);
    validator->initial_particle_segments = NULL;
    validator->initial_particle_count = state->particle_count;
    if (state->particle_count > 0) {
        validator->initial_particle_segments =
            (int *)malloc((size_t)state->particle_count * sizeof(int));
        if (!validator->initial_particle_segments) {
            return fail(validator, phase, "audit_allocation_failed",
                        "Could not allocate per-particle compartment identities.",
                        step, time_sigma, -1, -1, -1, -1,
                        (double)state->particle_count, NAN);
        }
        for (int i = 0; i < state->particle_count; ++i) {
            validator->initial_particle_segments[i] = segment_for_x(state, state->x[i]);
        }
    }
    return 1;
}

static int check_overlaps(ExperimentValidator *validator,
                          const ExperimentValidationState *state,
                          const char *phase,
                          int step,
                          double time_sigma) {
    if (state->particle_count <= 1) return 1;
    const double diameter = 2.0 * state->particle_radius;
    const double box_w = fmax(1.0, state->x_max_face - state->x_min_face);
    const double box_h = fmax(1.0, state->y_max - state->y_min);
    const double target_cells = fmax(64.0, 4.0 * (double)state->particle_count);
    const double cell_size = fmax(diameter, sqrt((box_w * box_h) / target_cells));
    const int gw = (int)fmax(1.0, ceil(box_w / cell_size));
    const int gh = (int)fmax(1.0, ceil(box_h / cell_size));
    if (gw > INT_MAX / gh) {
        return fail(validator, phase, "audit_grid_invalid",
                    "Overlap-audit grid dimensions overflowed.",
                    step, time_sigma, -1, -1, -1, -1, (double)gw, (double)gh);
    }
    const int cell_count = gw * gh;
    if (cell_count > validator->cell_capacity) {
        int *new_head = (int *)realloc(validator->cell_head,
                                      (size_t)cell_count * sizeof(int));
        if (!new_head) {
            return fail(validator, phase, "audit_allocation_failed",
                        "Could not allocate overlap-audit cell heads.",
                        step, time_sigma, -1, -1, -1, -1,
                        (double)cell_count, NAN);
        }
        validator->cell_head = new_head;
        validator->cell_capacity = cell_count;
    }
    if (state->particle_count > validator->next_capacity) {
        int *new_next = (int *)realloc(validator->cell_next,
                                      (size_t)state->particle_count * sizeof(int));
        if (!new_next) {
            return fail(validator, phase, "audit_allocation_failed",
                        "Could not allocate overlap-audit particle links.",
                        step, time_sigma, -1, -1, -1, -1,
                        (double)state->particle_count, NAN);
        }
        validator->cell_next = new_next;
        validator->next_capacity = state->particle_count;
    }
    for (int c = 0; c < cell_count; ++c) validator->cell_head[c] = -1;

    const double tolerance = fmax(1e-7, 1e-6 * diameter);
    for (int i = 0; i < state->particle_count; ++i) {
        int cx = (int)floor((state->x[i] - state->x_min_face) / cell_size);
        int cy = (int)floor((state->y[i] - state->y_min) / cell_size);
        if (cx < 0) cx = 0; else if (cx >= gw) cx = gw - 1;
        if (cy < 0) cy = 0; else if (cy >= gh) cy = gh - 1;
        for (int dy = -1; dy <= 1; ++dy) {
            const int ny = cy + dy;
            if (ny < 0 || ny >= gh) continue;
            for (int dx = -1; dx <= 1; ++dx) {
                const int nx = cx + dx;
                if (nx < 0 || nx >= gw) continue;
                for (int j = validator->cell_head[ny * gw + nx];
                     j >= 0; j = validator->cell_next[j]) {
                    const double ddx = state->x[i] - state->x[j];
                    const double ddy = state->y[i] - state->y[j];
                    const double distance2 = ddx * ddx + ddy * ddy;
                    if (distance2 < diameter * diameter) {
                        const double penetration =
                            diameter - sqrt(fmax(0.0, distance2));
                        if (penetration > tolerance) {
                            const double distance = sqrt(fmax(0.0, distance2));
                            const double rel_normal_velocity = (distance > 0.0)
                                ? (((state->vx[i] - state->vx[j]) * ddx +
                                    (state->vy[i] - state->vy[j]) * ddy) / distance)
                                : NAN;
                            char detail[256];
                            snprintf(detail, sizeof(detail),
                                     "Particles %d and %d overlap by %.9g px (tolerance %.9g px, relative normal velocity %.9g px/time).",
                                     i, j, penetration, tolerance, rel_normal_velocity);
                            return fail(validator, phase, "particle_particle_overlap", detail,
                                        step, time_sigma, i, j, -1, -1,
                                        penetration, tolerance);
                        }
                    }
                }
            }
        }
        const int cell = cy * gw + cx;
        validator->cell_next[i] = validator->cell_head[cell];
        validator->cell_head[cell] = i;
    }
    return 1;
}

int experiment_validator_check(ExperimentValidator *validator,
                               const ExperimentValidationState *state,
                               const char *phase,
                               int step,
                               double time_sigma,
                               int check_particle_identity) {
    if (!validator || validator->failure.failed) return 0;
    if (!state_shape_is_valid(validator, state, phase, step, time_sigma)) return 0;
    const double radius = state->particle_radius;
    const double tolerance = fmax(1e-6, 1e-6 * fmax(1.0, radius));

    if (!(radius > 0.0) || !isfinite(radius)) {
        return fail(validator, phase, "invalid_particle_radius",
                    "Particle radius is non-finite or non-positive.",
                    step, time_sigma, -1, -1, -1, -1, radius, 0.0);
    }
    if (!isfinite(time_sigma)) {
        return fail(validator, phase, "nonfinite_time",
                    "Simulation time became NaN or infinite.",
                    step, time_sigma, -1, -1, -1, -1, time_sigma, NAN);
    }
    if (state->particle_count != validator->initial_particle_count) {
        char detail[192];
        snprintf(detail, sizeof(detail), "Particle count changed from %d to %d.",
                 validator->initial_particle_count, state->particle_count);
        return fail(validator, phase, "particle_count_changed", detail,
                    step, time_sigma, -1, -1, -1, -1,
                    (double)state->particle_count,
                    (double)validator->initial_particle_count);
    }
    if (!(state->x_min_face < state->x_max_face) ||
        !(state->y_min < state->y_max)) {
        return fail(validator, phase, "boundary_order_invalid",
                    "Piston or outer boundary faces touched or crossed.",
                    step, time_sigma, -1, -1, -1, -1,
                    state->x_max_face - state->x_min_face, tolerance);
    }

    for (int w = 0; w < state->wall_count; ++w) {
        const double wx = state->wall_x[w];
        const double wv = state->wall_vx[w];
        if (!isfinite(wx) || !isfinite(wv)) {
            char detail[192];
            snprintf(detail, sizeof(detail), "Wall %d has a non-finite position or velocity.", w);
            return fail(validator, phase, "nonfinite_wall_state", detail,
                        step, time_sigma, -1, -1, w, -1, wx, NAN);
        }
        const double left_face = wx - state->wall_half_thickness;
        const double right_face = wx + state->wall_half_thickness;
        if (left_face <= state->x_min_face + tolerance ||
            right_face >= state->x_max_face - tolerance) {
            char detail[256];
            snprintf(detail, sizeof(detail),
                     "Wall %d contacted a piston/outer boundary at x=%.9g px.", w, wx);
            return fail(validator, phase, "wall_boundary_contact", detail,
                        step, time_sigma, -1, -1, w, -1,
                        fmin(left_face - state->x_min_face,
                             state->x_max_face - right_face), tolerance);
        }
        if (w > 0) {
            const double previous_right =
                state->wall_x[w - 1] + state->wall_half_thickness;
            const double face_gap = left_face - previous_right;
            if (face_gap <= tolerance) {
                char detail[256];
                snprintf(detail, sizeof(detail),
                         "Walls %d and %d touched, intersected, or changed order (face gap %.9g px).",
                         w - 1, w, face_gap);
                return fail(validator, phase, "wall_wall_contact", detail,
                            step, time_sigma, -1, -1, w - 1, w,
                            face_gap, tolerance);
            }
        }
    }

    for (int i = 0; i < state->particle_count; ++i) {
        const double x = state->x[i];
        const double y = state->y[i];
        const double vx = state->vx[i];
        const double vy = state->vy[i];
        if (!isfinite(x) || !isfinite(y) || !isfinite(vx) || !isfinite(vy)) {
            char detail[192];
            snprintf(detail, sizeof(detail), "Particle %d has a non-finite position or velocity.", i);
            return fail(validator, phase, "nonfinite_particle_state", detail,
                        step, time_sigma, i, -1, -1, -1, x, NAN);
        }
        if (x - radius < state->x_min_face - tolerance ||
            x + radius > state->x_max_face + tolerance ||
            y - radius < state->y_min - tolerance ||
            y + radius > state->y_max + tolerance) {
            char detail[256];
            snprintf(detail, sizeof(detail),
                     "Particle %d escaped a piston/outer boundary at (%.9g, %.9g) px.", i, x, y);
            return fail(validator, phase, "particle_boundary_escape", detail,
                        step, time_sigma, i, -1, -1, -1, x, NAN);
        }
        for (int w = 0; w < state->wall_count; ++w) {
            const double penetration = radius + state->wall_half_thickness -
                                       fabs(x - state->wall_x[w]);
            if (penetration > tolerance) {
                char detail[256];
                snprintf(detail, sizeof(detail),
                         "Particle %d penetrated wall %d by %.9g px "
                         "(particle x=%.17g vx=%.17g, wall x=%.17g vx=%.17g, "
                         "radius=%.17g, half-thickness=%.17g).",
                         i, w, penetration, x, vx,
                         state->wall_x[w], state->wall_vx[w],
                         radius, state->wall_half_thickness);
                return fail(validator, phase, "particle_wall_penetration", detail,
                            step, time_sigma, i, -1, w, -1,
                            penetration, tolerance);
            }
        }
        if (check_particle_identity && validator->initial_particle_segments) {
            const int current_segment = segment_for_x(state, x);
            if (current_segment != validator->initial_particle_segments[i]) {
                char detail[256];
                snprintf(detail, sizeof(detail),
                         "Particle %d changed compartment from %d to %d; aggregate counts can miss swaps.",
                         i, validator->initial_particle_segments[i], current_segment);
                return fail(validator, phase, "particle_crossed_wall", detail,
                            step, time_sigma, i, -1, -1, -1,
                            (double)current_segment,
                            (double)validator->initial_particle_segments[i]);
            }
        }
    }

    if (!isfinite(state->kinetic_energy) || state->kinetic_energy < 0.0 ||
        !isfinite(state->spring_energy) || state->spring_energy < -tolerance ||
        !isfinite(state->piston_work_left) ||
        !isfinite(state->piston_work_right)) {
        return fail(validator, phase, "nonfinite_energy_state",
                    "Kinetic energy, spring energy, or piston work became invalid.",
                    step, time_sigma, -1, -1, -1, -1,
                    state->kinetic_energy, 0.0);
    }
    if (state->forced_advance_events > 0) {
        char detail[256];
        snprintf(detail, sizeof(detail),
                 "EDMD hit an event avalanche/stagnation limit and forced %ld free-flight advance(s).",
                 state->forced_advance_events);
        return fail(validator, phase, "edmd_forced_advance", detail,
                    step, time_sigma, -1, -1, -1, -1,
                    (double)state->forced_advance_events, 0.0);
    }
    return check_overlaps(validator, state, phase, step, time_sigma);
}

const ExperimentValidationFailure *experiment_validator_failure(
    const ExperimentValidator *validator) {
    return validator ? &validator->failure : NULL;
}

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
                                        double limit) {
    if (!validator) return 0;
    return fail(validator, phase, reason, detail, step, time_sigma,
                particle_a, particle_b, wall_a, wall_b, value, limit);
}
