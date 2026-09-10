/* ##CHRIS: equilibrium pressure validation for the reference EDMD core.
 *
 * Measures the hard-disk compressibility factor Z = P/(rho kB T) with THREE
 * estimators that share no arithmetic, so agreement between them is evidence
 * rather than construction:
 *
 *   Z_pair   = 1 + W/(2 KE dt),  W = sum m|dv_n| sigma over PAIR collisions
 *   Z_wall_x = (I_L + I_R) W / (2 dt KE)      direct momentum flux, x walls
 *   Z_wall_y = (I_B + I_T) H / (2 dt KE)      direct momentum flux, y walls
 *
 * Wall impulses are deliberately NOT folded into W: the pair virial is already
 * the complete collisional virial, and adding them would double-count pressure.
 *
 * The pilot found the wall estimators sit above the pair virial by an amount
 * that shrinks roughly as perimeter/area (~1/sqrt(N)) -- consistent with a
 * hard-wall boundary finite-size contribution. This runner exists to test that
 * with seeds and block uncertainties instead of single trajectories.
 *
 * Geometry: square hard-wall box, no piston, no divider, no heat bath, no
 * driving. Box scales with sqrt(N) at fixed eta.
 *
 * Usage:
 *   calibrate:  pressure_validation <eta> <N> <seed> 0 0 0 - - calibrate
 *                 -> prints the largest chunk that ran completely clean
 *   produce:    pressure_validation <eta> <N> <seed> <blocks> <block_dt> <equil_t>
 *                                   <traj_csv> <block_csv> <chunk>
 *
 * Calibration is separated from production so it can be amortised over the
 * scientific seeds of one (eta,N): running a full disposable calibration before
 * EVERY seed cost ~280 simulated time units per trajectory, which made the dense
 * states a multi-day job. The driver calibrates on two disposable seeds and gives
 * production 0.8*min(safe_A, safe_B). Caching is a speed optimisation ONLY -- the
 * all-or-nothing health rule on production trajectories is unchanged.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "edmd_core/edmd.h"

#define PSI6_CUTOFF 1.4   /* neighbour shell in sigma, matches 00ALLINONE.c */

/* Returns |<psi6>| and, via out_local, <|psi6_i|>. The two differ when the box
   holds several ordered domains at different orientations: the global average
   cancels while local six-fold order is still strong, so reporting only the
   first can hide real structure. */
static double psi6_global(const EDMD* S, double* out_local){
    const EDMD_Particle* P = edmd_particles(S);
    const int N = edmd_params(S)->N;
    const double sigma = 2.0 * edmd_params(S)->radius;
    const double rc = PSI6_CUTOFF * sigma, rc2 = rc*rc;
    double sre=0.0, sim=0.0, slocal=0.0; int counted=0;
    for(int i=0;i<N;i++){
        double re=0.0, im=0.0; int nb=0;
        for(int j=0;j<N;j++){
            if(j==i) continue;
            double dx=P[j].x-P[i].x, dy=P[j].y-P[i].y;
            double d2=dx*dx+dy*dy;
            if(d2>rc2 || d2<=0.0) continue;
            double th=atan2(dy,dx);
            re+=cos(6.0*th); im+=sin(6.0*th); nb++;
        }
        if(nb>0){
            const double lr=re/nb, li=im/nb;
            sre+=lr; sim+=li;
            slocal += sqrt(lr*lr + li*li);   /* magnitude first: no cancellation */
            counted++;
        }
    }
    if(counted==0){ if(out_local) *out_local=NAN; return NAN; }
    sre/=counted; sim/=counted;
    if(out_local) *out_local = slocal/counted;
    return sqrt(sre*sre+sim*sim);
}

/* Advance in fixed steps. One large edmd_advance_to() can exhaust the per-call
   event budget and trigger a forced advance, which silently corrupts both the
   dynamics and the measurement window: the pilot saw forced_advance ~ 40 and
   ~4500 overlap repairs at eta = 0.5, N = 1600 with a 1.0 chunk.
   The chunk is FIXED for a production trajectory. It is chosen beforehand on a
   throwaway instance, so no adaptation ever runs inside a scientific run. */
static double g_chunk = 1.0;
static int seed_system(EDMD* S, double eta, unsigned long long seed, const char* mode);
static double g_cal_eta = 0.0;      /* eta/mode for calibration seeding */
static const char* g_cal_mode = "auto";

static void advance_fixed(EDMD* S, double dt_total){
    const double target = edmd_time(S) + dt_total;
    while(edmd_time(S) < target - 1e-12){
        double next = edmd_time(S) + g_chunk;
        if(next > target) next = target;
        edmd_advance_to(S, next);
    }
}

static long health_sum(const EDMD* S){
    return edmd_forced_advance_count(S) + edmd_overlap_repair_count(S)
         + edmd_clamp_repair_count(S)   + edmd_wall_overdue_count(S);
}

/* Pick a safe chunk on a DISPOSABLE system, then throw it away.
   The resolved pair-collision rate only guides the first guess: the internal
   per-call budget is also consumed by event-calendar work, so the rate is not a
   guarantee. The candidate is therefore verified directly -- any health event at
   all means the chunk is too coarse -- and a wide safety margin is used. */
/* Chunk calibration is a PERFORMANCE HEURISTIC, not part of the acceptance proof.
   Validity is decided solely by the production health counters: any forced
   advance / overlap repair / clamp repair / overdue wall anywhere in a scientific
   trajectory discards the whole trajectory. So calibration only has to answer
   "does this candidate look obviously too coarse?", and a short interval suffices
   -- the 560-unit production run is itself the long stress test. Verifying for
   280 simulated units before every 560-unit run was pure redundancy and roughly
   doubled the campaign cost. */
#define TARGET_EVENT_FRACTION 0.15      /* of the ~250k per-call budget */
#define ASSUMED_EVENT_BUDGET  250000.0
#define CALIB_MAX_ATTEMPTS    8
/* Returns a clean chunk, or -1.0 if no candidate survived (cell must not run). */
/* ##CHRIS 2026-09-08 -- calibrate on an EQUILIBRATED disposable instance.
   The previous version probed the event rate 0.02 units after seeding and
   verified for 5 units, both on a fresh lattice. That rate is lower than the
   equilibrated fluid's: in the 2026-09-07 dense campaign the avalanche warnings
   sat at t = 7..423 (median ~100), so a fresh-seed probe passed chunks 0.8 at
   N=900/1600 and 0.4 at N=1600 that then failed during equilibration -- 17
   discards, five empty cells. Now the disposable instance is equilibrated for
   CALIB_PREEQUIL_TIME first (at a conservative chunk, and that phase must itself
   be clean), the rate is probed on the equilibrated state, and each candidate is
   verified over CALIB_VERIFY_CHUNKS consecutive chunks. Still a heuristic:
   production validity is decided only by the production health counters. */
#define CALIB_PREEQUIL_TIME   100.0     /* simulated units before probing */
#define CALIB_VERIFY_CHUNKS   20        /* consecutive chunks per verification */
static double calibrate_chunk_disposable(const EDMD_Params* prm, unsigned long long seed){
    /* conservative chunk for the pre-equilibration: half the empirical 320/N
       rule, never above 0.05, so this phase does not itself blow the budget */
    double pre_chunk = 0.5 * 320.0 / (double)prm->N;
    if(pre_chunk > 0.05) pre_chunk = 0.05;

    double chunk = -1.0;
    for(int attempt=0; attempt<CALIB_MAX_ATTEMPTS; ++attempt){
        EDMD* T = edmd_create(prm);
        if(!T) return -1.0;
        if(!seed_system(T, g_cal_eta, seed, g_cal_mode)){ edmd_destroy(T); return -1.0; }

        /* (1) equilibrate the disposable instance; must be clean */
        g_chunk = pre_chunk;
        advance_fixed(T, CALIB_PREEQUIL_TIME);
        if(health_sum(T) != 0){
            fprintf(stderr,"calib: pre-equilibration itself unclean at chunk %g (N=%d)\n",
                    pre_chunk, prm->N);
            edmd_destroy(T); return -1.0;
        }

        /* (2) first attempt: probe the EQUILIBRATED event rate for the candidate */
        if(chunk < 0.0){
            const double t0=edmd_time(T); const long e0=edmd_virial_pair_events(T);
            g_chunk = pre_chunk;
            advance_fixed(T, 1.0);
            const double used=edmd_time(T)-t0; const long ev=edmd_virial_pair_events(T)-e0;
            chunk = (used>0.0 && ev>0)
                  ? (TARGET_EVENT_FRACTION*ASSUMED_EVENT_BUDGET)/((double)ev/used)
                  : 1.0;
            if(chunk > 1.0)  chunk = 1.0;
            if(chunk < 1e-5) chunk = 1e-5;
        }

        /* (3) verify the candidate over >= CALIB_VERIFY_CHUNKS consecutive chunks */
        g_chunk = chunk;
        advance_fixed(T, CALIB_VERIFY_CHUNKS * chunk);
        const long h = health_sum(T);
        edmd_destroy(T);
        if(h == 0) return chunk;
        chunk *= 0.5;
        if(chunk < 1e-5) break;
    }
    return -1.0;   /* calibration_failed */
}

/* Returns 1 on success. "auto" picks lattice where random insertion is not
   reliably terminating. */
static int seed_system(EDMD* S, double eta, unsigned long long seed, const char* mode){
    int use_lattice;
    if(strcmp(mode,"random")==0)       use_lattice = 0;
    else if(strcmp(mode,"lattice")==0) use_lattice = 1;
    else                               use_lattice = (eta >= 0.50);
    if(use_lattice) return edmd_init_lattice_gas(S, seed);
    edmd_init_random_gas(S, seed);
    return 1;
}

/* One block: advance block_dt, form Z from accumulator DIFFERENCES, write a row.
   Used for measurement blocks and, when requested, equilibration blocks. */
typedef struct { double prevT, prevW, prevI[4]; } BlockState;
static void block_state_init(BlockState* b, const EDMD* S){
    b->prevT=edmd_time(S); b->prevW=edmd_virial_accum(S);
    for(int w=0;w<4;w++) b->prevI[w]=edmd_wall_impulse(S,w);
}
static void run_block(EDMD* S, BlockState* b, double block_dt, int idx, int N,
                      double boxW, double boxH, double eta, unsigned long long seed,
                      FILE* fb, double* Zp_out, double* Zx_out, double* Zy_out,
                      double* T_out, double* ps_out, double* psl_out){
    advance_fixed(S, block_dt);
    const double t1 = edmd_time(S), dt = t1 - b->prevT;
    const double W  = edmd_virial_accum(S);
    double I[4]; for(int w=0;w<4;w++) I[w]=edmd_wall_impulse(S,w);
    const double T  = edmd_gas_temperature(S, 1.0, 1.0);
    const double KE = (double)N * T;
    const double dW = W - b->prevW;
    const double dI_x = (I[0]-b->prevI[0]) + (I[1]-b->prevI[1]);
    const double dI_y = (I[2]-b->prevI[2]) + (I[3]-b->prevI[3]);
    const double Zp = (dt>0 && KE>0) ? 1.0 + dW/(2.0*KE*dt) : NAN;
    const double Zx = (dt>0 && KE>0) ? dI_x * boxW / (2.0*dt*KE) : NAN;
    const double Zy = (dt>0 && KE>0) ? dI_y * boxH / (2.0*dt*KE) : NAN;
    double ps_loc=NAN; const double ps = psi6_global(S, &ps_loc);
    if(fb) fprintf(fb,"%.6g,%d,%llu,%d,%.6f,%.6f,%.8g,%.8g,%.8g,%.8g,%.6g,%.6g\n",
                   eta,N,seed,idx,b->prevT,t1,T,Zp,Zx,Zy,ps,ps_loc);
    b->prevW=W; for(int w=0;w<4;w++) b->prevI[w]=I[w]; b->prevT=t1;
    *Zp_out=Zp; *Zx_out=Zx; *Zy_out=Zy; *T_out=T; *ps_out=ps; *psl_out=ps_loc;
}

int main(int argc, char** argv){
    if(argc < 9){
        fprintf(stderr,"usage: %s eta N seed blocks block_dt equil_t traj_csv block_csv "
                       "[chunk|calibrate]\n",argv[0]);
        return 2;
    }
    const int calibrate_only = (argc>9 && strcmp(argv[9],"calibrate")==0);
    const double chunk_arg   = (argc>9 && !calibrate_only) ? atof(argv[9]) : 0.0;
    /* Seeding mode. edmd_init_random_gas() HANGS above eta ~ 0.55 (unbounded
       rejection sampling), so dense states must use the checked lattice seeder.
       Explicit "random"/"lattice" exists so the two can be compared at densities
       where both work -- if they disagree on Z, the seeder is biasing the result
       and the dense numbers could not be trusted. */
    const char* seed_mode = (argc>10) ? argv[10] : "auto";
    /* ##CHRIS 2026-09-08 -- TASK 4 (equilibration ladder). When set, the
       equilibration phase is ALSO written to the block CSV, one row per block_dt,
       with NEGATIVE block_index (-n_eq .. -1, chronological when sorted), same
       columns. Default off: production schema and behaviour unchanged. */
    const int equil_log = (argc>11) ? atoi(argv[11]) : 0;
    const double eta      = atof(argv[1]);
    const int    N        = atoi(argv[2]);
    const unsigned long long seed = strtoull(argv[3],NULL,10);
    const int    nblocks  = atoi(argv[4]);
    const double block_dt = atof(argv[5]);
    const double equil_t  = atof(argv[6]);
    const char*  traj_csv = argv[7];
    const char*  block_csv= argv[8];

    const double r = 0.5;
    const double A = N * M_PI * r * r / eta;
    const double L = sqrt(A);                 /* square box, scales as sqrt(N) */

    EDMD_Params p; memset(&p,0,sizeof p);
    p.boxW=L; p.boxH=L; p.radius=r; p.N=N;
    p.pp_collisions_enabled=1; p.particle_mass=1.0; p.kB=1.0;
    /* explicitly: no piston, no divider, no heat bath -- all zero from memset */

    g_cal_eta = eta; g_cal_mode = seed_mode;
    if(calibrate_only){
        /* Disposable calibration only. Prints the largest chunk that ran a full
           verification interval with zero health events, or FAILED. */
        const double c = calibrate_chunk_disposable(&p, seed);
        if(c <= 0.0){ printf("FAILED\n"); return 3; }
        printf("%.10g\n", c);
        return 0;
    }
    /* Production uses the chunk the driver cached for this (eta,N). Falling back
       to self-calibration keeps the binary usable standalone, but the campaign
       always supplies it. */
    if(chunk_arg > 0.0){
        g_chunk = chunk_arg;
    } else {
        const double c = calibrate_chunk_disposable(&p, seed);
        if(c <= 0.0){ fprintf(stderr,"calibration_failed eta=%g N=%d\n",eta,N); return 3; }
        g_chunk = c;
    }

    EDMD* S = edmd_create(&p);
    if(!S){ fprintf(stderr,"edmd_create failed (eta=%g N=%d)\n",eta,N); return 1; }
    if(!seed_system(S, eta, seed, seed_mode)){
        fprintf(stderr,"seeding_failed eta=%g N=%d mode=%s\n",eta,N,seed_mode);
        edmd_destroy(S); return 4;
    }

    FILE* fb = fopen(block_csv,"a");
    if(equil_log && block_dt > 0.0){
        /* equilibration written as negative-index blocks; accumulators reset so
           differences are meaningful. Health is still judged cumulatively below. */
        edmd_reset_virial(S);
        BlockState bs; block_state_init(&bs, S);
        const int n_eq = (int)ceil(equil_t / block_dt);
        double d1,d2,d3,d4,d5,d6;
        for(int k=0;k<n_eq;k++){
            const double this_dt = (k==n_eq-1) ? (equil_t - k*block_dt) : block_dt;
            if(this_dt <= 0.0) break;
            run_block(S,&bs,this_dt,-(n_eq-k),N,p.boxW,p.boxH,eta,seed,fb,&d1,&d2,&d3,&d4,&d5,&d6);
        }
    } else {
        advance_fixed(S, equil_t);
    }
    /* An accepted validation trajectory must be clean from INITIALIZATION
       onward, not merely during the measurement window: a forced advance in
       equilibration reaches the measurement through nonphysical dynamics, and
       for a validation result there is no reason to accept that ambiguity. */
    if(health_sum(S) != 0){
        fprintf(stderr,"DISCARD eta=%g N=%d seed=%llu: health during equilibration "
                       "(fa=%ld orep=%ld crep=%ld wovd=%ld) chunk=%g\n",
                eta,N,seed,edmd_forced_advance_count(S),edmd_overlap_repair_count(S),
                edmd_clamp_repair_count(S),edmd_wall_overdue_count(S),g_chunk);
        edmd_destroy(S);
        return 1;
    }

    /* Blocks are contiguous windows; per-block values come from DIFFERENCES of
       the running accumulators, so no core changes are needed and the blocks are
       genuinely independent samples of the same stationary process. */
    BlockState bs; 
    edmd_reset_virial(S);
    block_state_init(&bs, S);

    double sZp=0, sZx=0, sZy=0, sT=0, sPsi=0, sPsiL=0;
    double s2Zp=0, s2Zx=0, s2Zy=0;
    int nb_ok=0;
    long pair_ev_total=0;

    for(int b=0;b<nblocks;b++){
        double Zp,Zx,Zy,T,ps,ps_loc;
        run_block(S,&bs,block_dt,b,N,p.boxW,p.boxH,eta,seed,fb,&Zp,&Zx,&Zy,&T,&ps,&ps_loc);
        if(isfinite(Zp)&&isfinite(Zx)&&isfinite(Zy)){
            sZp+=Zp; sZx+=Zx; sZy+=Zy; sT+=T; sPsi+=ps; sPsiL+=ps_loc;
            s2Zp+=Zp*Zp; s2Zx+=Zx*Zx; s2Zy+=Zy*Zy; nb_ok++;
        }
    }
    if(fb) fclose(fb);
    pair_ev_total = edmd_virial_pair_events(S);

    /* Cumulative, i.e. from initialization: any health event anywhere in this
       trajectory disqualifies it. No adaptation, no rescue, no partial credit. */
    const long m_fa = edmd_forced_advance_count(S);
    const long m_or = edmd_overlap_repair_count(S);
    const long m_cr = edmd_clamp_repair_count(S);
    const long m_wo = edmd_wall_overdue_count(S);
    const long health = m_fa + m_or + m_cr + m_wo;
    /* A trajectory with no collisions at all is not a measurement -- it yields
       Z_pair = 1 and Z_wall = 0 and would otherwise sail through every check.
       Caught when a mis-seeded lattice left the event calendar empty. */
    const int  moved  = (pair_ev_total > 0) &&
                        (edmd_wall_events(S,0)+edmd_wall_events(S,1)
                        +edmd_wall_events(S,2)+edmd_wall_events(S,3) > 0);
    const int  valid  = (health==0 && nb_ok==nblocks && moved) ? 1 : 0;
    if(!moved) fprintf(stderr,"INERT eta=%g N=%d seed=%llu: no collisions recorded\n",
                       eta,N,seed);

    #define MEAN(s) ((nb_ok>0)?(s)/nb_ok:NAN)
    #define SEM(s,s2) ((nb_ok>1)? sqrt(fmax(0.0,((s2)/nb_ok-((s)/nb_ok)*((s)/nb_ok)))/(nb_ok-1)) : NAN)
    FILE* ft = fopen(traj_csv,"a");
    if(ft){
        fprintf(ft,"%.6g,%d,%llu,%.6f,%.6f,%.4f,%.4f,%ld,%ld,%ld,%ld,%ld,"
                   "%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%d,%ld,%ld,%ld,%ld,%d\n",
                eta,N,seed,p.boxW,p.boxH,equil_t,nblocks*block_dt,
                pair_ev_total,
                edmd_wall_events(S,0),edmd_wall_events(S,1),
                edmd_wall_events(S,2),edmd_wall_events(S,3),
                MEAN(sZp),SEM(sZp,s2Zp),MEAN(sZx),SEM(sZx,s2Zx),
                MEAN(sZy),SEM(sZy,s2Zy),MEAN(sT),MEAN(sPsi),MEAN(sPsiL),
                nb_ok,
                m_fa,m_or,m_cr,m_wo,valid);
        fclose(ft);
    }
    printf("PROD eta=%.4g N=%d seed=%llu Z_pair=%.5f+-%.5f Z_wx=%.5f Z_wy=%.5f blocks=%d chunk=%g health=%ld valid=%d\n",
           eta,N,seed,MEAN(sZp),SEM(sZp,s2Zp),MEAN(sZx),MEAN(sZy),nb_ok,g_chunk,health,valid);
    edmd_destroy(S);
    return valid?0:1;
}
