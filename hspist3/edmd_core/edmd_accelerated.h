#ifndef EDMD_ACCELERATED_API_H
#define EDMD_ACCELERATED_API_H

/* Accelerated EDMD backend API.
   Same types as edmd.h, but all functions are exported as edmd_acc_* so it can
   be linked together with the default edmd.c backend. */

#include "edmd.h"

EDMD*  edmd_acc_create(const EDMD_Params* prm);
void   edmd_acc_destroy(EDMD* S);
void   edmd_acc_init_random_gas(EDMD* S, unsigned long long seed);

double                edmd_acc_time(const EDMD* S);
const EDMD_Params*    edmd_acc_params(const EDMD* S);
const EDMD_Particle*  edmd_acc_particles(const EDMD* S);
int                   edmd_acc_count(const EDMD* S);

double edmd_acc_advance_to(EDMD* S, double t_target);
double edmd_acc_advance_pp_only_to(EDMD* S, double t_target);

void   edmd_acc_reschedule_all(EDMD* S);
void   edmd_acc_reschedule_all_pp_only(EDMD* S);

void   edmd_acc_config_dividers(EDMD* S, int count, const double* cx, const double* thickness);
void   edmd_acc_set_divider_motions(EDMD* S, int count, const double* mass, const double* vx);
void   edmd_acc_config_divider(EDMD* S, int enabled, double cx, double thickness);
void   edmd_acc_set_divider_motion(EDMD* S, double mass, double vx);

void   edmd_acc_config_pistons(EDMD* S,
                               int hasL, double xL, double vxL, double mL,
                               int hasR, double xR, double vxR, double mR);

double edmd_acc_total_kinetic_energy(const EDMD* S, double mass);
double edmd_acc_gas_temperature(const EDMD* S, double mass, double kB);

double edmd_acc_work_divider(const EDMD* S);
double edmd_acc_work_divider_i(const EDMD* S, int divider_index);
double edmd_acc_work_pistonL(const EDMD* S);
double edmd_acc_work_pistonR(const EDMD* S);
double edmd_acc_heat_bath(const EDMD* S);
void   edmd_acc_reset_work(EDMD* S);

void   edmd_acc_divider_resolve_overlaps(EDMD* S);

#endif /* EDMD_ACCELERATED_API_H */
