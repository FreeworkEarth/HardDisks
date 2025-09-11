#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "edmd_core/edmd.h"

/* Simple CLI holder */
typedef enum { ENG_EDMD=0, ENG_STEPS=1 } Engine;

typedef struct {
    Engine engine;
    int    N;
    double R, boxW, boxH;
    unsigned int seed;
    int    fps;
    double tmax;
} Cli;

static void set_defaults(Cli* c){
    c->engine = ENG_EDMD;
    c->N = 800;
    c->R = 0.05;
    c->boxW = 12.0;
    c->boxH = 6.0;
    c->seed = 42;
    c->fps  = 240;
    c->tmax = 5.0; /* seconds */
}

static void parse_cli(int argc, char** argv, Cli* c){
    set_defaults(c);
    for(int i=1;i<argc;i++){
        if(strncmp(argv[i],"--engine=",9)==0){
            const char* v = argv[i]+9;
            if(strcmp(v,"edmd")==0) c->engine=ENG_EDMD;
            else if(strcmp(v,"steps")==0) c->engine=ENG_STEPS;
        } else if(sscanf(argv[i],"--N=%d",&c->N)==1){
        } else if(sscanf(argv[i],"--radius=%lf",&c->R)==1){
        } else if(sscanf(argv[i],"--box=%lfx%lf",&c->boxW,&c->boxH)==2){
        } else if(sscanf(argv[i],"--seed=%u",&c->seed)==1){
        } else if(sscanf(argv[i],"--fps=%d",&c->fps)==1){
        } else if(sscanf(argv[i],"--tmax=%lf",&c->tmax)==1){
        } else {
            fprintf(stderr, "Unknown arg: %s\n", argv[i]);
        }
    }
}

static void run_edmd(const Cli* c){
    /* set up EDMD */
    EDMD_Params p = {
        .boxW=c->boxW, .boxH=c->boxH, .radius=c->R, .N=c->N,
        .max_events_hint=c->N*16,
        .cell_size=0.0
    };
    EDMD* S = edmd_create(&p);
    edmd_init_random_gas(S, c->seed);

    double dt = 1.0 / (double)c->fps;
    double t_print_next = 0.0;

    printf("# engine=edmd N=%d R=%.6f box=%.3fx%.3f fps=%d tmax=%.3f seed=%u\n",
           c->N, c->R, c->boxW, c->boxH, c->fps, c->tmax, c->seed);

    while(edmd_time(S) < c->tmax){
        double t_target = edmd_time(S) + dt;
        edmd_advance_to(S, t_target);

        /* print once per ~1s */
        if(edmd_time(S) >= t_print_next){
            double U = edmd_total_kinetic_energy(S, 1.0);
            printf("t=%.3f  U_kin=%.6f\n", edmd_time(S), U);
            fflush(stdout);
            t_print_next += 1.0;
        }
    }

    edmd_destroy(S);
}

static void run_steps_placeholder(const Cli* c){
    (void)c;
    fprintf(stderr,
        "[steps] placeholder: wire this to your 00ALLINONE.c entry point.\n"
        "         For now, use --engine=edmd to run the event-driven version.\n");
}

int main(int argc, char** argv){
    Cli c; parse_cli(argc, argv, &c);
    if(c.engine == ENG_EDMD) run_edmd(&c);
    else run_steps_placeholder(&c);
    return 0;
}