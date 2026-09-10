/* ##CHRIS: one dense smoke case -- reports the diagnostics the seeder tests ask
   for. psi6 is REPORTED, not asserted: a rectangular lattice and a hex fallback
   legitimately start with very different orientational order, and what matters
   scientifically is what equilibration does to it. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "edmd_core/edmd.h"
static double psi6(const EDMD* S, double* loc){
    const EDMD_Particle* P=edmd_particles(S); const int N=edmd_params(S)->N;
    const double sig=2.0*edmd_params(S)->radius, rc=1.4*sig, rc2=rc*rc;
    double sre=0,sim=0,sl=0; int cnt=0;
    for(int i=0;i<N;i++){ double re=0,im=0; int nb=0;
        for(int j=0;j<N;j++){ if(j==i) continue;
            double dx=P[j].x-P[i].x, dy=P[j].y-P[i].y, d2=dx*dx+dy*dy;
            if(d2>rc2||d2<=0) continue; double th=atan2(dy,dx);
            re+=cos(6*th); im+=sin(6*th); nb++; }
        if(nb>0){ double lr=re/nb, li=im/nb; sre+=lr; sim+=li; sl+=sqrt(lr*lr+li*li); cnt++; } }
    if(!cnt){ if(loc)*loc=NAN; return NAN; }
    sre/=cnt; sim/=cnt; if(loc)*loc=sl/cnt; return sqrt(sre*sre+sim*sim);
}
int main(int argc,char**argv){
    if(argc<3) return 2;
    double eta=atof(argv[1]); int N=atoi(argv[2]); double r=0.5;
    double A=N*M_PI*r*r/eta, L=sqrt(A);
    EDMD_Params p; memset(&p,0,sizeof p);
    p.boxW=L;p.boxH=L;p.radius=r;p.N=N;p.pp_collisions_enabled=1;
    p.particle_mass=1.0;p.kB=1.0;
    clock_t t0=clock();
    EDMD* S=edmd_create(&p); if(!S) return 1;
    if(!edmd_init_lattice_gas(S,20260907)){ printf("SEED_FAILED\n"); return 1; }
    double init_s=(double)(clock()-t0)/CLOCKS_PER_SEC;
    const EDMD_Particle* P=edmd_particles(S);
    double sep=1e30,gap=1e30;
    for(int i=0;i<N;i++){ for(int j=i+1;j<N;j++){ double dx=P[i].x-P[j].x,dy=P[i].y-P[j].y,d=sqrt(dx*dx+dy*dy); if(d<sep)sep=d; }
        double g[4]={P[i].x-r,L-r-P[i].x,P[i].y-r,L-r-P[i].y}; for(int k=0;k<4;k++) if(g[k]<gap) gap=g[k]; }
    double pl,pg=psi6(S,&pl);
    /* small chunks: the campaign calibrates these; 1.0 is known to blow the
       event budget at these densities. */
    for(int i=0;i<200;i++) edmd_advance_to(S, edmd_time(S)+0.05);
    long we=0; for(int w=0;w<4;w++) we+=edmd_wall_events(S,w);
    long h=edmd_forced_advance_count(S)+edmd_overlap_repair_count(S)
          +edmd_clamp_repair_count(S)+edmd_wall_overdue_count(S);
    printf("eta=%-5.2f N=%-6d %-9.3f %-11.4f %-10.2e %-9.4f %-9.4f %-8ld %-8ld %ld\n",
           eta,N,init_s,sep,gap,pg,pl,edmd_virial_pair_events(S),we,h);
    edmd_destroy(S); return 0;
}
