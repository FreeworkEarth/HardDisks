/* ##CHRIS: regression tests for edmd_init_lattice_gas().
 *
 * Each test here locks out a defect that actually occurred and that produced
 * plausible-looking output rather than an obvious error:
 *   - missing reschedule    -> inert system, Z_pair = 1.000, Z_wall = 0, valid=1
 *   - row-major corner fill -> geometrically legal, but the local density was
 *                              ~2.5x the requested eta with most of the box empty
 *   - placement at exactly x = R -> gap = 0 with the wall face, i.e. an overdue
 *                              wall collision at t = 0
 * "No overlaps" alone did not catch any of these, which is why coverage and
 * liveness are tested explicitly.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "edmd_core/edmd.h"

static int fails = 0;
static void ck(int cond, const char* name, const char* detail){
    printf("  [%s] %-56s %s\n", cond?"PASS":"FAIL", name, cond?"":detail);
    if(!cond) fails++;
}
static EDMD_Params mkp(double eta, int N, double r){
    EDMD_Params p; memset(&p,0,sizeof p);
    const double A = N*M_PI*r*r/eta, L = sqrt(A);
    p.boxW=L; p.boxH=L; p.radius=r; p.N=N;
    p.pp_collisions_enabled=1; p.particle_mass=1.0; p.kB=1.0;
    return p;
}
static double min_pair_sep(const EDMD* S){
    const EDMD_Particle* P = edmd_particles(S);
    const int N = edmd_params(S)->N;
    double m = 1e30;
    for(int i=0;i<N;i++) for(int j=i+1;j<N;j++){
        const double dx=P[i].x-P[j].x, dy=P[i].y-P[j].y;
        const double d=sqrt(dx*dx+dy*dy); if(d<m) m=d;
    }
    return m;
}
static double min_wall_gap(const EDMD* S){
    const EDMD_Particle* P = edmd_particles(S);
    const EDMD_Params* q = edmd_params(S);
    double m = 1e30;
    for(int i=0;i<q->N;i++){
        double g[4] = { P[i].x - q->radius, q->boxW - q->radius - P[i].x,
                        P[i].y - q->radius, q->boxH - q->radius - P[i].y };
        for(int k=0;k<4;k++) if(g[k]<m) m=g[k];
    }
    return m;
}
static long health(const EDMD* S){
    return edmd_forced_advance_count(S)+edmd_overlap_repair_count(S)
         + edmd_clamp_repair_count(S)+edmd_wall_overdue_count(S);
}

int main(void){
    printf("lattice seeder regression tests\n\n");

    /* 1. VALID GEOMETRY AT t=0, at densities where random insertion cannot finish */
    printf("1. valid geometry at t=0\n");
    {   const double etas[3]={0.60,0.69,0.72}; const int Ns[3]={400,900,1600};
        for(int k=0;k<3;k++){
            EDMD_Params p = mkp(etas[k],Ns[k],0.5);
            EDMD* S = edmd_create(&p);
            char nm[128];
            int ok = S ? edmd_init_lattice_gas(S,4242+k) : 0;
            snprintf(nm,sizeof nm,"eta=%.2f N=%d seeds successfully",etas[k],Ns[k]);
            ck(ok, nm, "initializer returned failure");
            if(!ok){ if(S) edmd_destroy(S); continue; }
            const EDMD_Particle* P = edmd_particles(S);
            int finite=1, inside=1;
            for(int i=0;i<p.N;i++){
                if(!isfinite(P[i].x)||!isfinite(P[i].y)||!isfinite(P[i].vx)||!isfinite(P[i].vy)) finite=0;
                if(P[i].x < p.radius || P[i].x > p.boxW-p.radius ||
                   P[i].y < p.radius || P[i].y > p.boxH-p.radius) inside=0;
            }
            snprintf(nm,sizeof nm,"eta=%.2f all coordinates finite",etas[k]);
            ck(finite,nm,"non-finite value");
            snprintf(nm,sizeof nm,"eta=%.2f all centers inside wall region",etas[k]);
            ck(inside,nm,"center outside allowed region");
            const double sep = min_pair_sep(S), gap = min_wall_gap(S);
            snprintf(nm,sizeof nm,"eta=%.2f no pair overlap (min sep %.4f > 2r)",etas[k],sep);
            ck(sep > 2.0*p.radius, nm, "overlap at t=0");
            snprintf(nm,sizeof nm,"eta=%.2f strictly off the walls (min gap %.2e)",etas[k],gap);
            ck(gap > 0.0, nm, "particle touching a wall");
            edmd_destroy(S);
        }
    }

    /* 2. EVENT CALENDAR IS LIVE -- guards the missing-reschedule bug */
    printf("\n2. event calendar is live after seeding\n");
    {   EDMD_Params p = mkp(0.60,400,0.5);
        EDMD* S = edmd_create(&p); edmd_init_lattice_gas(S,7);
        for(int i=0;i<10;i++) edmd_advance_to(S, edmd_time(S)+1.0);
        const long pe = edmd_virial_pair_events(S);
        long we=0; for(int w=0;w<4;w++) we += edmd_wall_events(S,w);
        ck(pe > 0, "pair collisions occur after a short advance", "system is inert");
        ck(we > 0, "wall collisions occur after a short advance", "system is inert");
        ck(health(S) == 0, "no health events during the short advance", "health nonzero");
        edmd_destroy(S); }

    /* 3. DETERMINISM */
    printf("\n3. determinism\n");
    {   EDMD_Params p = mkp(0.60,400,0.5);
        EDMD *A=edmd_create(&p), *B=edmd_create(&p), *C=edmd_create(&p);
        edmd_init_lattice_gas(A,999); edmd_init_lattice_gas(B,999); edmd_init_lattice_gas(C,1000);
        const EDMD_Particle *pa=edmd_particles(A), *pb=edmd_particles(B), *pc=edmd_particles(C);
        int same=1, diff=0;
        for(int i=0;i<p.N;i++){
            if(pa[i].x!=pb[i].x||pa[i].y!=pb[i].y||pa[i].vx!=pb[i].vx||pa[i].vy!=pb[i].vy) same=0;
            if(pa[i].vx!=pc[i].vx || pa[i].x!=pc[i].x) diff=1;
        }
        ck(same,"same seed reproduces identical state exactly","seed not deterministic");
        ck(diff,"different seed changes jitter/velocities","different seed gave identical state");
        edmd_destroy(A); edmd_destroy(B); edmd_destroy(C); }

    /* 4. COUNT AND BOX COVERAGE -- locks out corner packing */
    printf("\n4. particle count and box coverage\n");
    {   const double etas[2]={0.30,0.60};
        for(int k=0;k<2;k++){
            EDMD_Params p = mkp(etas[k],400,0.5);
            EDMD* S=edmd_create(&p); edmd_init_lattice_gas(S,55);
            const EDMD_Particle* P=edmd_particles(S);
            double xmn=1e30,xmx=-1e30,ymn=1e30,ymx=-1e30;
            for(int i=0;i<p.N;i++){
                if(P[i].x<xmn)xmn=P[i].x; if(P[i].x>xmx)xmx=P[i].x;
                if(P[i].y<ymn)ymn=P[i].y; if(P[i].y>ymx)ymx=P[i].y;
            }
            const double availW=p.boxW-2*p.radius, availH=p.boxH-2*p.radius;
            const double fx=(xmx-xmn)/availW, fy=(ymx-ymn)/availH;
            char nm[128];
            snprintf(nm,sizeof nm,"eta=%.2f spans x (%.2f of box, need >0.8)",etas[k],fx);
            ck(fx>0.8,nm,"packed into part of the box");
            snprintf(nm,sizeof nm,"eta=%.2f spans y (%.2f of box, need >0.8)",etas[k],fy);
            ck(fy>0.8,nm,"packed into part of the box");
            edmd_destroy(S);
        }
    }

    /* 5. WALL MARGIN -- locks out gap==0 at t=0 producing overdue wall events */
    printf("\n5. wall-margin regression\n");
    {   EDMD_Params p = mkp(0.69,900,0.5);
        EDMD* S=edmd_create(&p); edmd_init_lattice_gas(S,31337);
        ck(min_wall_gap(S) > 0.0,"min wall gap strictly positive at t=0","touching a wall");
        for(int i=0;i<8;i++) edmd_advance_to(S, edmd_time(S)+1.0);
        ck(edmd_wall_overdue_count(S)==0,"wall_overdue == 0 after advancing","overdue wall events");
        edmd_destroy(S); }

    /* 6. IMPOSSIBLE GEOMETRY -- must refuse, and refuse promptly */
    printf("\n6. impossible geometry is refused\n");
    {   EDMD_Params p; memset(&p,0,sizeof p);
        p.boxW=3.0; p.boxH=3.0; p.radius=0.5; p.N=500;   /* cannot possibly fit */
        p.pp_collisions_enabled=1; p.particle_mass=1.0; p.kB=1.0;
        EDMD* S=edmd_create(&p);
        const int ok = S ? edmd_init_lattice_gas(S,1) : 0;
        ck(ok==0,"returns 0 for an impossible box (no infinite retry)","did not refuse");
        if(S) edmd_destroy(S); }

    printf("\n%s  (%d failures)\n", fails?"TESTS FAILED":"ALL TESTS PASSED", fails);
    return fails?1:0;
}
