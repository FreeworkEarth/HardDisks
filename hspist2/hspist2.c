// hspist.c - april 29, 2011
// observe effect of piston motion on hard disk gas
// hspist_wall.c - collect data on wall density profile
// hspist_mb.c - collect maxwell-boltzmann statistics, reverse statistics
// rob shaw
#include <stdio.h>
#include <stdlib.h>		// for exit()
#include <math.h>
#include <string.h>
#include <SDL2/SDL.h>
#include "glgraph.h"
#include "rstuff.h"
#include "myrand4.h"
#include "globals.h"


#define FULLSCR 1			// open program in full screen



unsigned long int _seed = SEED;

REAL X[NSPHERES],Y[NSPHERES],VX[NSPHERES],VY[NSPHERES],DM2[NSPHERES];
REAL Radius[NSPHERES],LastBumpTime[NSPHERES],BumpInterval[NSPHERES];
unsigned long Histo[NENTS],Histo1[NENTS],Histo2[NENTS],Histowall[NENTS],HistoW[NENTS];
long Cells[NCELLS * CSIZE];
long DiskNbrs[NSPHERES * CSIZE];
unsigned long Vout[SIZEX * SIZEY];			// output display image

REAL xw1=XW1,xw2=XW2,yw1=YW1,yw2=YW2;
REAL rad1,rad2;			// disk radii
long nrad1,nrad2;		// number of disks with RAD1, RAD2...
long nspheres = NBEGIN;	// number of disks, must be smaller than NSPHERES
long seed = SEED;
long bumpnum = 0;
long bumptot = 0;
long wallbump = 0;
long rowl,nrows,ncells;	// neighborhood indices
long slowctl = 1;		// slow display flag...
long keyslowctl = 1;	// changed in kbd.c
long Slowctl = 1;		// screen refresh slow display...
int dispctl = 0;		// display control
int Dispctl = 0;		// second display control
int radctl = 0;			// for adjusting radii
long ncount = 0;		// iteration counter
long nruns = 0;
REAL total_time = 0;	// simulation time
REAL xpist1 = XW1;		// left piston x position
REAL vpist1 = 0;		// left piston x velocity
REAL xpist2 = XW2;		// right piston x position
REAL vpist2 = 0;		// right piston x velocity
REAL dxpist = .001;		// piston x step per display
REAL Ktot;				// total kinetic energy
REAL KE0,KE1,KE2;		// KE initial, after first motion, after second motion
long histctl = 0;		// display histograms
REAL vmax = 3.0;		// maximum for velocity histogram
REAL Phend;				// henderson approximation to hard disk pressure
long trailctl = 0;		// trail control
char *mypath;			// program path
char pathName[1024];	// path name constructed with loadpath()
long datacnt = 0;		// number of times histogram totals are summed
REAL binwidth;			// histogram bin width
long nbins;				// number of histogram bins used
int forward = 1;		// flag for trajectory direction
long ndone;				// number of ncounts per run
int pistctl = 0;		// piston protocol control...
REAL llo,lhi,rlo,rhi;	// left, right piston starting, ending positions
REAL lvel,rvel;			// left, right piston velocities
REAL pdelay;			// time delay for second (right) piston
REAL dt1,dt2;			// time for first, second piston motion
int modectl = 1;		// demo mode...
int ppctl = 0;			// program protocol control...
REAL dx;				// distance of piston stroke...

struct bump{
	long disknum;
	long iden;	// positive for second disk, negative for walls, force centers...
	REAL fcenx;	// force center coordinates...
	REAL fceny;
};
struct bump bumps[MAXBUMPS];




// Function definitions
void cleanup() {
    // Add your cleanup code here
    // For example, freeing allocated memory, closing files, etc.
    printf("Cleanup function called.\n");
}

void Display(SDL_Window* _window, SDL_Renderer* _renderer)
{
#if FULLSCR
    static int fadecount = 0;
#endif
    long i;

    // initialize
    if (ncount == 0)
    {
        xpist1 = XW1;
        xpist2 = rlo;
        initgauss(nrad1, nrad2, RAD1, 0);
        reseteng();
        pdelay = 0.26;
        lvel = 1;
        rvel = 1;
        KE0 = kengy();
        KE1 = KE2 = 0;
        pdelay = dt1 + 0.00;

        switch (modectl)
        {
        case 1:
            lvel = 0.1;
            rvel = 0.1;
            dt1 = (lhi - llo) / lvel;
            printf("dt1 %f\n", dt1);
            pdelay = dt1 + 0.30;
            break;
        case 2:
            lvel = 2;
            rvel = 0.1;
            dt1 = (lhi - llo) / lvel;
            printf("dt1 %f\n", dt1);
            pdelay = dt1 + 0.30;
            break;
        case 3:
            lvel = 2;
            rvel = 2;
            dt1 = (lhi - llo) / lvel;
            printf("dt1 %f\n", dt1);
            pdelay = dt1 + 1;
            break;
        case 4:
            lvel = 2;
            rvel = 2;
            dt1 = (lhi - llo) / lvel;
            printf("dt1 %f\n", dt1);
            pdelay = dt1 + 0.08;
            break;
        case 5:
            lvel = 2;
            rvel = 2;
            dt1 = (lhi - llo) / lvel;
            printf("dt1 %f\n", dt1);
            pdelay = 0;
            break;
        }
        keyslowctl = 2;
    }

    for (i = 0; i < slowctl; i++)
    {
#if CBNDYS
        updateCB(X, Y, VX, VY, DELT);
#else
        update4b(X, Y, VX, VY, DELT);
#endif
        bumpnum = 0;

        if (!(ncount % NSORT))
        {
            sortcells();
#if CBNDYS
            mkdisknbrlistCB();
#else
            mkdisknbrlist();
#endif
        }
#if CBNDYS
        bouncecheckCB();
#else
        bouncecheckd();
#endif

        bigbox4(X, Y);
        dobumpsa(X, Y);
        ++ncount;
        total_time = ncount * DELT;
        pproto();

        if (histctl == 1) ldiskhist(Histo1);
        if (histctl == 2) ldiskhist(Histo2);
        if (histctl == 3) ldspeedhist(Histo1);

#if FCTL
        if (ncount > NTRANS)
        {
            if (!(ncount % DATAINT))
            {
                ldwallhist(Histowall);
                ++datacnt;
                if (datacnt >= DATACNT)
                {
                    writedata();
                    exit(0);
                }
            }
        }
#endif
    }

    display();
    if (!(ncount % keyslowctl)) slowctl = keyslowctl;

#if FULLSCR
    if (!fadecount++)
        SDL_SetWindowFullscreen(_window, SDL_WINDOW_FULLSCREEN);
    if (fadecount == 4)
        fadein();
#endif

    SDL_RenderPresent(_renderer);
}








int main(int argc, char* argv[])
{
    long i;
    char* mypath = (char *)argv[0]; // for loadpath()
    loadpath(mypath, "data/");

    #if(FLOATRES == 32)
    printf("reals are floats\n");
    #elif (FLOATRES == 80)
    printf("reals are long doubles\n");
    #else
    printf("reals are doubles\n");
    #endif

    if (SDL_Init(SDL_INIT_VIDEO) < 0)
    {
        fprintf(stderr, "Could not initialize SDL: %s\n", SDL_GetError());
        return 1;
    }

    SDL_Window* _window = SDL_CreateWindow("Simulation", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, 800, 600, SDL_WINDOW_SHOWN);
    if (!_window)
    {
        fprintf(stderr, "Could not create window: %s\n", SDL_GetError());
        SDL_Quit();
        return 1;
    }

    SDL_Renderer* _renderer = SDL_CreateRenderer(_window, -1, SDL_RENDERER_ACCELERATED);
    if (!_renderer)
    {
        fprintf(stderr, "Could not create renderer: %s\n", SDL_GetError());
        SDL_DestroyWindow(_window);
        SDL_Quit();
        return 1;
    }

    #if FULLSCR
    SDL_SetWindowFullscreen(_window, SDL_WINDOW_FULLSCREEN);
    atexit(restorecolor);
    atexit(cleanup);
    fadeout(); // fade to black for smooth screen change...
    #endif

    scale(MFAC * -2, MFAC * 2, MFAC * -1, MFAC * 1);
    nrad1 = 800; nrad2 = 0;
    rad1 = RAD1; rad2 = RAD2;
    dx = .4;
    xpist1 = XW1;
    xpist2 = XW2 - dx;
    vpist1 = 0;
    llo = XW1;
    lhi = XW1 + dx;
    rlo = XW2 - dx;
    rhi = XW2;
    initgauss(nrad1, nrad2, RAD1, 0);
    rowl = ceil((XW2 - XW1) / LS);
    nrows = (YW2 - YW1) / LS;
    ncells = rowl * nrows;
    printf("rowl %ld  nrows %ld  ncells %ld\n", rowl, nrows, ncells);
    for (i = 0; i < NENTS; i++) Histo[i] = Histo1[i] = Histo2[i] = Histowall[i] = 0;

    #if FCTL
    initfile();
    writeheader();
    #endif

    while (1)
    {
        SDL_Event event;
        while (SDL_PollEvent(&event))
        {
            if (event.type == SDL_QUIT)
            {
                SDL_DestroyRenderer(_renderer);
                SDL_DestroyWindow(_window);
                SDL_Quit();
                return 0;
            }
        }

        Display(_window, _renderer);
    }

    return 0;
}








   


void loadpath(char *arg, char *name)
{
char *lastSeparator;

strcpy(pathName, arg);
lastSeparator = strrchr(pathName, '/');
if(lastSeparator == NULL)
	strcpy(pathName, name);
else
	strcpy(lastSeparator+1, name); 
printf("path is %s\n", pathName);
}

void display()
{
long i;
REAL ymax=0;

if(trailctl != 1)
	{
	color(WHITE);
	clear();
	}

if(trailctl == 0)
//if(0)
	{
for(i=0;i<nspheres;i++)
	{
	REAL rad;

	rad = Radius[i];
	if(rad == rad1) color(RED);
	else color(BLUE);
	circf(X[i],Y[i],rad);
		
#if CB1 || CBNDYS
	if((YW2 - Y[i]) < rad)		// for circular bndy conditions (top edge)
		{
		circf(X[i],Y[i]-2*YW2,rad);
		color(BLACK);
		circ(X[i],Y[i]-2*YW2,rad);
		}
	if((Y[i] - YW1) < rad)		// for circular bndy conditions (bottom edge)
		{
		circf(X[i],Y[i]+2*YW2,rad);
		color(BLACK);
		circ(X[i],Y[i]+2*YW2,rad);
		}
#endif
#if CBNDYS
	if((XW2 - X[i]) < rad)		// for circular bndy conditions (right edge)
		{
		color(RED);		// corner disks can stay black
		circf(X[i]-(XW2-xpist1),Y[i],rad);
		color(BLACK);
		circ(X[i]-(XW2-xpist1),Y[i],rad);
		}
	if((X[i] - xpist1) < rad)		// for circular bndy conditions (left edge)
		{
		color(RED);
		circf(X[i]+(XW2-xpist1),Y[i],rad);
		color(BLACK);
		circ(X[i]+(XW2-xpist1),Y[i],rad);
		}
#endif
	color(BLACK);
	circ(X[i],Y[i],rad);
	}
	}

if(trailctl == 1)
	{
	color(BLACK);
	for(i=0;i<nspheres;i++)
		dotat(X[i],Y[i]);
	}

if(trailctl == 2)
	{
	for(i=0;i<nspheres;i++)
		{
		float x1,y1,x2,y2;
		float fac = 0.05;
		REAL maxe;
		maxe = VX[i]*VX[i] + VY[i]*VY[i];
		if(maxe > 0.0)
			{
			color(RED);
			x1 = X[i] - fac * VX[i];
			y1 = Y[i] - fac * VY[i];
//			x2 = X[i] + fac * VX[i];
//			y2 = Y[i] + fac * VY[i];
			x2 = X[i];
			y2 = Y[i];
			move(x1,y1); pendn();
			move(x2,y2); penup();
			color(BLACK);
			circf(x2,y2,.015);
			}
		}
	}

if(trailctl == 3)
	{
	for(i=0;i<nspheres;i++)
		{
		float x1,y1,x2,y2;
		float fac = 0.05;
//		float radfac = .01;
		float radfac = .00005;
		REAL maxe,rad;
		maxe = VX[i]*VX[i] + VY[i]*VY[i];
//		rad  = radfac * sqrt(DM2[i]);
//		rad /= BumpInterval[i];
//		rad = radfac * DM2[i] / BumpInterval[i];
		rad = radfac / BumpInterval[i];
//		if(i==100) printf("%f\n",BumpInterval[i]);
		if(maxe > 0.0)
			{
			color(RED);
			x1 = X[i] - fac * VX[i];
			y1 = Y[i] - fac * VY[i];
//			x2 = X[i] + fac * VX[i];
//			y2 = Y[i] + fac * VY[i];
			x2 = X[i];
			y2 = Y[i];
			move(x1,y1); pendn();
			move(x2,y2); penup();
			color(BLACK);
			circf(x2,y2,rad);
			}
		}
	}

color(WHITE);					// mask off left piston side
rectf(xpist1-.10,YW1,xpist1,YW2);
color(BLACK);					// draw left piston
rectf(xpist1-.05,YW1,xpist1,YW2);
color(WHITE);					// mask off box (y direction)
rectf(MFAC*-2,MFAC*-1,MFAC*2,-1);
rectf(MFAC*-2,MFAC*1,MFAC*2,1);
//rectf(MFAC*2,-1,2,1);			// mask off right edge
rectf(xpist2,YW1,xpist2+.10,YW2);	// mask off right piston side
color(BLACK);					// draw right piston
rectf(xpist2,YW1,xpist2+.05,YW2);
//color(BLACK);					// draw box outline
move(XW1,YW1); pendn();
move(XW1,YW2);move(XW2,YW2);move(XW2,YW1);move(XW1,YW1);penup();

if(dispctl)		// draw neighborhoods
	{
	int nx,ny,i;
	float x,y;
	nx = (XW2-XW1)/LS;
	ny = (YW2-YW1)/LS;
	for(i=0;i<=nx;i++)
		{
		x = XW1 + i*LS;
		move(x,YW1); pendn(); move(x,YW2); penup();
		}
	for(i=0;i<=ny;i++)
		{
		y = YW1 + i*LS;
		move(XW1,y); pendn(); move(XW2,y); penup();
		}
	}


if(histctl)
	{
//	ldhist0(Histo);
//	plthist(Histo,ROWL,100);
	if(histctl==1)  ymax = copyhist(Histo1,Histo,ROWL);
	if(histctl==2)  ymax = copyclearhist(Histo2,Histo,ROWL);
	if(histctl==3)  ymax = copyhist(Histo1,Histo,ROWL);
//	if(histctl == 3) printf("ymax %f\n",ymax);
	if(histctl==4)  ymax = copyhist(Histowall,Histo,ROWL);
	plthist(Histo,ROWL,(5*(ymax/4)));
	}

drawnum(ncount,800,38);			// labels
//drawnum(nrad1,100,38);
//drawnum(nrad2,300,38);
Printf("\n%d disks",nrad1);
labit(50,38);
Ktot = kengy();					// kinetic energy

#if 0
#if FLOATRES == 80
Printf("\nKtot %.1Lf",Ktot);
#else
Printf("\nKtot %.1f",Ktot);
#endif
labit(250,480);
Phend = hender();				// henderson approximation to hard disk pressure
//Printf("\nPh %.1f",Phend);
//labit(100,480);
#if FLOATRES == 80
Printf("\nvpist1 %.2Lf",vpist1);
#else
Printf("\nvpist1 %.2f",vpist1);
#endif
labit(100,480);
#endif

Printf("\nKE %.1f",Ktot);
labit(40,480);
Printf("\nKE0 %.1f",KE0);
labit(320,480);
if(KE1 != 0.0)
	{
	Printf("\nKE1 %.1f",KE1);
	labit(520,480);
	}
if(KE2 != 0.0)
	{
	Printf("\nKE2 %.1f",KE2);
	labit(720,480);
//	labelcolor(RED);
//	Printf("\ndE %.1f",KE2-KE0);
//	labit(760,480);
//	labelcolor(BLUE);
	labelcolor(RED);
//	Printf("\nKE2 - KE0 = %.1f",KE2-KE0);
	Printf("\nDissipation = %.1f",KE2-KE0);
	labit(350,38);
	labelcolor(BLUE);
	}
else
	{
	labelcolor(RED);
	switch(modectl)
		{
		case 1:
			Printf("\nquasistatic");
		break;
		case 2:
			Printf("\nimpulse");
		break;
		case 3:
			Printf("\nimpulse, long delay");
		break;
		case 4:
			Printf("\nimpulse, good delay");
		break;
		case 5:
			Printf("\nimpulse, too short delay");
		break;
		}
	labit(350,38);
	labelcolor(BLUE);
	}
//Printf("\nbumps %d",bumpnum);
//labit(700,480);

}

// randomly distribute within piston volume, start with maxwell-boltzmann statistics,
// gaussian distributed random velocities, with unit mean, kT = 1/2 per dof
// don't worry about overlap
// we are using the globally defined random seed _seed
#define M_SQRT1_2   0.707106781186547524400844362104849039  /* 1/sqrt(2) */
void initgauss(long nr1, long nr2, REAL r1, REAL r2)
{
double xx,yy;
long i;

for(i=0;i<nr1;i++)
	{
	xx = f_rand();
//	xx = (2-xpist1)*xx + xpist1;	// from xpist1 to 2 (xpist1 can be negative)
	xx = (xpist2-xpist1)*xx + xpist1;	// from xpist1 to xpist2 (xpist1 can be negative)
	yy = f1_rand();				// from -1 to 1
	X[i] = xx; Y[i] = yy;
	Radius[i] = r1;
	gaussrand(&xx,&yy);
	VX[i] = xx; VY[i] = yy;
	}
nrad1 = nr1; rad1 = r1;
// X[0] = .01; X[1] = 1.99;
// Y[0] = -.2; Y[1] = .2;
//VX[0] = VX[1] = 0;
//VY[0] = -.5; VY[1] = 0;
//VY[1] = 0; VY[0] = .5;

for(i=nr1;i<nr2;i++)
	{
	xx = f_rand();
//	xx = (2-xpist1)*xx + xpist1;	// from xpist1 to 2 (xpist1 can be negative)
	xx = (xpist2-xpist1)*xx + xpist1;	// from xpist1 to xpist2 (xpist1 can be negative)
	yy = f1_rand();				// from -1 to 1
	X[i] = xx; Y[i] = yy;
	Radius[i] = r2;
	gaussrand(&xx,&yy);
	VX[i] = xx; VY[i] = yy;
	}
nrad2 = nr2; rad2 = r2;

nspheres = nrad1 + nrad2;
}

// distribute randomly over right side, don't worry about overlap
void initspheres9(long nr1,long nr2,long seed)		// number of rad1, rad2, random seed...
{
long i;
REAL theta;
REAL xx,yy,delx,dely;
register unsigned long _seed;

_seed = seed;
for(i=0;i<17*seed;i++)
	myrand();

delx = XW2 - MR - 2*rad1;
dely = YW2 - YW1 - 2*rad1;
for(i=0;i<nr1;i++)
	{
	xx = f_rand();
	xx *= delx;
	yy = f_rand();
	yy *= dely;
	X[i] = MR + xx + rad1;
	Y[i] = YW1 + yy + rad1;
	Radius[i] = rad1;
	theta = ang_rand();
//	VX[i] = 1.0 * cos(theta);
//	VY[i] = 1.0 * sin(theta);
	VX[i] = M_SQRT2 * cos(theta);	// correct normalizatiion for 1/2 kT per dof, i think
	VY[i] = M_SQRT2 * sin(theta);	// M_SQRT2 is in math.h
	}

delx = XW2 - MR - 2*rad2;
dely = YW2 - YW1 - 2*rad2;
for(i=nr1;i<nr1+nr2;i++)
	{
	xx = f_rand();
	xx *= delx;
	yy = f_rand();
	yy *= dely;
	X[i] = MR + xx + rad2;
	Y[i] = YW1 + yy + rad2;
	Radius[i] = rad2;
	theta = ang_rand();
	VX[i] = 1.0 * cos(theta);
	VY[i] = 1.0 * sin(theta);
	}

nrad1 = nr1; nrad2 = nr2;
nspheres = nr1 + nr2;
}

// time stepper...
void update4b(REAL *XP,REAL *YP,REAL *VXP,REAL *VYP,REAL dtt)
{
long i;
register REAL dt;

dt = dtt;

for(i=0;i<nspheres;i++)
	{
	XP[i] += dt*VXP[i];
	YP[i] += dt*VYP[i];
#if CB1
	if(YP[i] >= YW2)				// the equality case can happen!
		{
		YP[i] -= (YW2-YW1);
//		if(i == markctl) Xold = 100;		// mark beginning of trail
		}
	if(YP[i] < YW1)
		{
		YP[i] += (YW2-YW1);
//		if(i == markctl) Xold = 100;		// mark beginning of trail
		}
#endif
	}
xpist1 += dt * vpist1;					// update piston position
xpist2 += dt * vpist2;
}

// time stepper, circular boundary conditions...
// horizontal box is defined on the left by the piston position
void updateCB(REAL *XP,REAL *YP,REAL *VXP,REAL *VYP,REAL dtt)
{
long i;
register REAL dt;

dt = dtt;

for(i=0;i<nspheres;i++)
	{
	XP[i] += dt*VXP[i];
	YP[i] += dt*VYP[i];

	if(YP[i] >= YW2)				// the equality case can happen!
		YP[i] -= (YW2-YW1);
	if(YP[i] < YW1)
		YP[i] += (YW2-YW1);
	if(XP[i] >= XW2)
		XP[i] -= (XW2 - xpist1);
	if(XP[i] < xpist1)
		XP[i] += (XW2 - xpist1);
	}
xpist1 += dt * vpist1;					// update piston position
xpist2 += dt * vpist2;
}

// puts boundary collisions in bumps structure...
// accumulate back wall "pressure" in VtotR and VtotL
void bigbox2(REAL *XP, REAL *YP)
{
long i;
REAL rad;

for(i=0;i<nspheres;i++)
	{
	rad = Radius[i];

	if((XP[i]-rad <= XW1) && VX[i] < 0)
		{
		bumps[bumpnum].disknum = i;
		bumps[bumpnum].iden = -1;	// collision with left edge...
		++bumpnum;
		++bumptot;
//		VtotL += 2 * VX[i];			// the bounce adds to the pressure!						
		}
	if((XP[i]+rad >= XW2) && VX[i] > 0)
		{
		bumps[bumpnum].disknum = i;
		bumps[bumpnum].iden = -2;	// collision with right edge...
		++bumpnum;
		++bumptot;
//		VtotR += 2 * VX[i];			// the bounce adds to the pressure!						
		}
//#if !CB1
	if((YP[i]-rad <= YW1) && VY[i] < 0)
		{
		bumps[bumpnum].disknum = i;
		bumps[bumpnum].iden = -3;	// collision with bottom edge...
		++bumpnum;
		++bumptot;		
		}
	if((YP[i]+rad >= YW2) && VY[i] > 0)
		{
		bumps[bumpnum].disknum = i;
		bumps[bumpnum].iden = -4;	// collision with top edge...
		++bumpnum;
		++bumptot;
		}
//#endif
	}
}

// puts boundary collisions in bumps structure...
// left wall is moveable (piston!)
// also right wall!
void bigbox4(REAL *XP, REAL *YP)
{
long i;
REAL rad;

for(i=0;i<nspheres;i++)
	{
	rad = Radius[i];
	
	if((X[i]-rad <= xpist1) && (VX[i] - vpist1) < 0)
//	if((X[i]-rad <= XW1) && VX[i] < 0)
		{
		bumps[bumpnum].disknum = i;
//		bumps[bumpnum].iden = -1;	// collision with left edge...
		bumps[bumpnum].iden = -9;	// collision with left edge (moveable!)
		++bumpnum;
		++bumptot;
//		VtotL += 2 * VX[i];			// the bounce adds to the pressure!						
		}
	if((X[i]+rad >= xpist2) && (VX[i] - vpist2) > 0)
//	if((X[i]+rad >= XW2) && VX[i] > 0)
		{
		bumps[bumpnum].disknum = i;
//		bumps[bumpnum].iden = -2;	// collision with right edge...
		bumps[bumpnum].iden = -10;	// collision with right edge (moveable!)
		++bumpnum;
		++bumptot;
//		VtotR += 2 * VX[i];			// the bounce adds to the pressure!						
		}
#if !CB1
	if((Y[i]-rad <= YW1) && VY[i] < 0)
		{
		bumps[bumpnum].disknum = i;
		bumps[bumpnum].iden = -3;	// collision with bottom edge...
		++bumpnum;
		++bumptot;		
		}
	if((Y[i]+rad >= YW2) && VY[i] > 0)
		{
		bumps[bumpnum].disknum = i;
		bumps[bumpnum].iden = -4;	// collision with top edge...
		++bumpnum;
		++bumptot;
		}
#endif
	}
}

// do all collisions on bump list...
// variable radius version...
void dobumpsa(REAL *XP, REAL *YP)
{
long i,dnum,iden;
REAL xcen,ycen;

for(i=0;i<bumpnum;i++)
	{
	dnum = bumps[i].disknum;
	iden = bumps[i].iden;
	if(iden >= 0)
		{
//		dobouncea(XP,YP,dnum,iden);
//		dobounceCB(XP,YP,dnum,iden);
		dobounceCBa(XP,YP,dnum,iden);		// this one keeps momentum change...
		}
	else
		{	
		switch(iden)
			{
			case -1:						// collision with left edge
			case -2:						// collision with right edge
			case -5:						// collision with right membrane edge
			case -6:						// collision with left membrane edge
				VX[dnum] = -VX[dnum];
				break;
			case -3:						// collision with bottom edge
			case -4:						// collision with top edge
				VY[dnum] = -VY[dnum];
				break;
			case -8:						// collision with force center
				xcen = bumps[i].fcenx;
				ycen = bumps[i].fceny;
				dofcbouncea(XP,YP,dnum,xcen,ycen);
				break;
			case -9:						// collision with moving left edge
				VX[dnum] = 2*vpist1 - VX[dnum];
//				VX[dnum] = -VX[dnum];
				break;
			case -10:						// collision with moving right edge
				VX[dnum] = 2*vpist2 - VX[dnum];
			}
		}
	}
}

// perform bounce between argument disks...
// cleaned up a little...
void dobouncea(REAL *XP,REAL *YP,long d1,long d2)
{
	register REAL x0,y0,x,y,delx,dely,vx,vy,vx0,vy0;
	register REAL vcmx,vcmy;
	register REAL R,sino,coso,R2;
	register REAL v1x,v1y,v2x,v2y,A;
	x0 = XP[d1];
	y0 = YP[d1];
	vx0 = VX[d1];
	vy0 = VY[d1];
	
	x = XP[d2];
	y = YP[d2];
	vx = VX[d2];
	vy = VY[d2];
	delx = x - x0;
	dely = y - y0;
	R2 = delx*delx + dely*dely;

#if CB1 || CBNDYS
		{
		REAL R3,dely1;
		dely1 = dely + (YW2-YW1);
		R3 = delx*delx + dely1*dely1;
		if(R3 < R2)
			{
			R2 = R3;
			dely = dely1;
			}
		dely1 = dely - (YW2-YW1);
		R3 = delx*delx + dely1*dely1;
		if(R3 < R2)
			{
			R2 = R3;
			dely = dely1;
			}
		}
#endif
#if CBNDYS
		{
		REAL R3,delx1;
		delx1 = delx + (XW2-xpist1);
//		printf("%f %f\n",delx,delx1);
		R3 = delx1*delx1 + dely*dely;
		if(R3 < R2)
			{
			R2 = R3;
			delx = delx1;
			}
		delx1 = delx - (XW2-xpist1);
		R3 = delx1*delx1 + dely*dely;
		if(R3 < R2)
			{
			R2 = R3;
			delx = delx1;
			}
		}
#endif	
//printf("%f %f %f\n",delx,dely,sqrt(R2));
	// move into CM coordinates
	vcmx = .5 * (vx0 + vx);
	vcmy = .5 * (vy0 + vy);
	// first sphere	
	v1x = vx0 - vcmx;
	v1y = vy0 - vcmy;
	R = sqrt(R2);
	coso = delx/R;
	sino = dely/R;
	A = -(v1x*coso + v1y*sino);
	if(A < 0)	// make sure velocity is inward...
		{
		v2x = v1x + 2*A*coso;
		v2y = v1y + 2*A*sino;
		// back to non CM coords, change velocities...
		VX[d1] = v2x + vcmx;
		VY[d1] = v2y + vcmy;
		// second sphere
		VX[d2] = -v2x + vcmx;
		VY[d2] = -v2y + vcmy;
		}
}

// perform bounce between argument disks...
// cleaned up a little more...
// xpist1 seems to have to be on a cell boundary for this to work correctly...
void dobounceCB(REAL *XP,REAL *YP,long d1,long d2)
{
	register REAL x1,y1,x2,y2,delx,dely,vx1,vy1,vx2,vy2;
	register REAL vcmx,vcmy;
	register REAL R,sino,coso,R2;
	register REAL v1x,v1y,v2x,v2y,A;
	x1 = XP[d1];
	y1 = YP[d1];
	vx1 = VX[d1];
	vy1 = VY[d1];
	
	x2 = XP[d2];
	y2 = YP[d2];
	vx2 = VX[d2];
	vy2 = VY[d2];
	
	delx = x2-x1;
	dely = y2-y1;

#if CB1 || CBNDYS
	{
	REAL dely1;
	if(y2 > y1)
		{
		dely1 = dely - (YW2-YW1);
		if(dely+dely1 > 0)
			dely = dely1;
		}
	else
		{
		dely1 = dely + (YW2-YW1);
		if(dely+dely1 < 0)
			dely = dely1;
		}
	}
#endif
#if CBNDYS
	{
	REAL delx1;
	if(x2 > x1)
		{
		delx1 = delx - (XW2-xpist1);
		if(delx+delx1 > 0)
			delx = delx1;
		}
	else
		{
		delx1 = delx + (XW2-xpist1);
		if(delx+delx1 < 0)
			delx = delx1;
		}
	}
#endif

	R2 = delx*delx + dely*dely;
//printf("%f %f %f\n",delx,dely,sqrt(R2));
	// move into CM coordinates
	vcmx = .5 * (vx1 + vx2);
	vcmy = .5 * (vy1 + vy2);
	// first sphere	
	v1x = vx1 - vcmx;
	v1y = vy1 - vcmy;
	R = sqrt(R2);
	coso = delx/R;
	sino = dely/R;
	A = -(v1x*coso + v1y*sino);
	if(A < 0)	// make sure velocity is inward...
		{
		v2x = v1x + 2*A*coso;
		v2y = v1y + 2*A*sino;
		// back to non CM coords, change velocities...
		VX[d1] = v2x + vcmx;
		VY[d1] = v2y + vcmy;
		// second sphere
		VX[d2] = -v2x + vcmx;
		VY[d2] = -v2y + vcmy;
		}
}

// perform bounce between argument disks...
// cleaned up a little more...
// xpist1 seems to have to be on a cell boundary for this to work correctly...
// this one computes (square of) momentum change, places in DM2[]
// also fills in LastBumpTime[],BumpInterval[]
void dobounceCBa(REAL *XP,REAL *YP,long d1,long d2)
{
	register REAL x1,y1,x2,y2,delx,dely,vx1,vy1,vx2,vy2;
	register REAL vcmx,vcmy;
	register REAL R,sino,coso,R2;
	register REAL v1x,v1y,v2x,v2y,A;
	x1 = XP[d1];
	y1 = YP[d1];
	vx1 = VX[d1];
	vy1 = VY[d1];
	
	x2 = XP[d2];
	y2 = YP[d2];
	vx2 = VX[d2];
	vy2 = VY[d2];
	
	delx = x2-x1;
	dely = y2-y1;

#if CB1 || CBNDYS
	{
	REAL dely1;
	if(y2 > y1)
		{
		dely1 = dely - (YW2-YW1);
		if(dely+dely1 > 0)
			dely = dely1;
		}
	else
		{
		dely1 = dely + (YW2-YW1);
		if(dely+dely1 < 0)
			dely = dely1;
		}
	}
#endif
#if CBNDYS
	{
	REAL delx1;
	if(x2 > x1)
		{
		delx1 = delx - (XW2-xpist1);
		if(delx+delx1 > 0)
			delx = delx1;
		}
	else
		{
		delx1 = delx + (XW2-xpist1);
		if(delx+delx1 < 0)
			delx = delx1;
		}
	}
#endif

	R2 = delx*delx + dely*dely;
//printf("%f %f %f\n",delx,dely,sqrt(R2));
	// move into CM coordinates
	vcmx = .5 * (vx1 + vx2);
	vcmy = .5 * (vy1 + vy2);
	// first sphere	
	v1x = vx1 - vcmx;
	v1y = vy1 - vcmy;
	R = sqrt(R2);
	coso = delx/R;
	sino = dely/R;
	A = -(v1x*coso + v1y*sino);
	if(A < 0)	// make sure velocity is inward...
		{
		v2x = v1x + 2*A*coso;
		v2y = v1y + 2*A*sino;
		// back to non CM coords, change velocities...
		VX[d1] = v2x + vcmx;
		VY[d1] = v2y + vcmy;
		// second sphere
		VX[d2] = -v2x + vcmx;
		VY[d2] = -v2y + vcmy;
		// compute momentum changes
		REAL dmx,dmy;
		dmx = VX[d1] - vx1;
		dmy = VX[d1] - vy1;
//		DM2[d1] = dmx*dmx + dmy*dmy;
		DM2[d1] = -A;
		dmx = VX[d2] - vx2;
		dmy = VX[d2] - vy2;
//		DM2[d2] = dmx*dmx + dmy*dmy;
		DM2[d1] = -A;
		BumpInterval[d1] = total_time - LastBumpTime[d1];
		LastBumpTime[d1] = total_time;
		BumpInterval[d2] = total_time - LastBumpTime[d2];
		LastBumpTime[d2] = total_time;
		}
}

// perform bounce between argument disk and argument force center...
void dofcbouncea(REAL *XP,REAL *YP,long d1,REAL xcen,REAL ycen)
{
register REAL x0,y0,x,y,delx,dely,vx0,vy0;
register REAL R,sino,coso,R2,A,rsq;

x0 = XP[d1];
y0 = YP[d1];
vx0 = VX[d1];
vy0 = VY[d1];
rsq = Radius[d1];
rsq *= rsq;

x = xcen;
y = ycen;
delx = x - x0;
dely = y - y0;
R2 = delx*delx + dely*dely;
// already checked going into this routine...
if(R2 > rsq) return;
	
R = sqrt(R2);
coso = delx/R;
sino = dely/R;
A = -(vx0*coso + vy0*sino);

// change velocities...
if(A < 0)
	{
	VX[d1] = vx0 + 2*A*coso;
	VY[d1] = vy0 + 2*A*sino;
	}
}

// sort disks into cells, put counts and disk indices into Cells[]...
// there are ncells entries, each of size CSIZE.
// the first number in each entry is the number of disks in that cell,
// then follows a list of the disk indices, maximum length, CSIZE-1
void sortcells(void)
{
long i,ix,iy,index,nn;

for(i=0;i<ncells;i++)				// these are the neighbor counts
	Cells[CSIZE * i] = 0;
	
for(i=0;i<nspheres;i++)
	{
	ix = (X[i] - XW1) * oLS;
	iy = (Y[i] - YW1) * oLS;
	index = iy * rowl + ix;
	index *= CSIZE;
	
	nn = (Cells[index++])++;		// load and increment cell count...
//	if(nn > 25) printf("%d oops!\n",nn);
	Cells[index+nn] = i;			// install index in cell!
	}		
}

// this one runs through the list of disks in each cell,
// assembled by sortcells() and placed in Cells[],
// and makes a list of neighbors of each disk, which is
// put in DiskNbrs[].  the format is, number of neighbors,
// followed by a list of the neighbor indices, close packed.
// only neighbors with a higher index than the current disk
// are listed, as each pair only has to be checked once.
// the collision detection check is done in bouncecheckd()
void mkdisknbrlist()
{
long i,j,index,*DNcptr,*DNlptr,cnt,ix,iy,nn,ii;

DNcptr = DNlptr = DiskNbrs;
for(i=0;i<nspheres;i++)
	{
	ix = (X[i] - XW1) * oLS;
	iy = (Y[i] - YW1) * oLS;
	index = iy * rowl + ix;
	
	index *= CSIZE;
	nn = Cells[index++];				// center cell index count
	
	++DNlptr;						// jump past neighbor count entry
	cnt = 0;
	for(j=0;j<nn;j++)
		{
		ii = Cells[index++];
//		if(abs(ii) > 10000) printf("wha1? %ld\n",ii);
		if(i < ii)					// only record neighbors with a higher disk number
			{
			*DNlptr++ = ii;			// enter a neighbor onto the list!
			++cnt;					// increment neighbor count
			}
		}

	if(iy != 0)						// if not on lowest row, examine lower three cells
		{
		if((ix%rowl) != 0)			// if not all the way to the left, examine lower left
			{
			index = (iy-1)*rowl + ix - 1;
			index *= CSIZE;
			nn = Cells[index++];
			for(j=0;j<nn;j++)
				{
				ii = Cells[index++];
//				if(abs(ii) > 10000) printf("wha2? %ld\n",ii);
				if(i < ii)					// only record neighbors with a higher disk number
					{
					*DNlptr++ = ii;			// enter a neighbor onto the list!
					++cnt;
					}
				}
			}

		index = (iy-1)*rowl + ix;	// examine lower center
		index *= CSIZE;
		nn = Cells[index++];
		for(j=0;j<nn;j++)
			{
			ii = Cells[index++];
//			if(abs(ii) > 10000) printf("wha3? %ld\n",ii);
			if(i < ii)					// only record neighbors with a higher disk number
				{
				*DNlptr++ = ii;			// enter a neighbor onto the list!
				++cnt;					// increment neighbor count
				}
			}
		
		if((ix%rowl) != (rowl-1))	// if not all the way to the right, examine lower right
			{
			index = (iy-1)*rowl + ix + 1;
			index *= CSIZE;
			nn = Cells[index++];
			for(j=0;j<nn;j++)
				{
				ii = Cells[index++];
//				if(abs(ii) > 10000) printf("wha4? %ld\n",ii);
				if(i < ii)					// only record neighbors with a higher disk number
					{
					*DNlptr++ = ii;			// enter a neighbor onto the list!
					++cnt;
					}
				}
			}
		}
#if CB1
	else					// we are on lowest row, check highest row (circular boundary conditions)
		{
		if((ix%rowl) != 0)			// if not all the way to the left, examine top row left
			{
			index = (nrows-1)*rowl + ix - 1;
			index *= CSIZE;
			nn = Cells[index++];
			for(j=0;j<nn;j++)
				{
				ii = Cells[index++];
				if(i < ii)					// only record neighbors with a higher disk number
					{
					*DNlptr++ = ii;			// enter a neighbor onto the list!
					++cnt;
					}
				}
			}

		index = (nrows-1)*rowl + ix;	// examine top row center
		index *= CSIZE;
		nn = Cells[index++];
		for(j=0;j<nn;j++)
			{
			ii = Cells[index++];
			if(i < ii)					// only record neighbors with a higher disk number
				{
				*DNlptr++ = ii;			// enter a neighbor onto the list!
				++cnt;
				}
			}
		
		if((ix%rowl) != (rowl-1))	// if not all the way to the right, examine top row right
			{
			index = (nrows-1)*rowl + ix + 1;
			index *= CSIZE;
			nn = Cells[index++];
			for(j=0;j<nn;j++)
				{
				ii = Cells[index++];
				if(i < ii)					// only record neighbors with a higher disk number
					{
					*DNlptr++ = ii;			// enter a neighbor onto the list!
					++cnt;
					}
				}
			}
		}
#endif
	if((ix%rowl) != 0)			// if not all the way to the left, examine middle left
		{
		index = iy*rowl + ix - 1;
		index *= CSIZE;
		nn = Cells[index++];
		for(j=0;j<nn;j++)
			{
			ii = Cells[index++];
//			if(abs(ii) > 10000) printf("wha5? %ld\n",ii);
			if(i < ii)					// only record neighbors with a higher disk number
				{
				*DNlptr++ = ii;			// enter a neighbor onto the list!
				++cnt;
				}
			}
		}	

	if((ix%rowl) != (rowl-1))	// if not all the way to the right, examine middle right
		{
		index = iy*rowl + ix + 1;
		index *= CSIZE;
		nn = Cells[index++];
		for(j=0;j<nn;j++)
			{
			ii = Cells[index++];
//			if(abs(ii) > 10000) printf("wha6? %ld\n",ii);
			if(i < ii)					// only record neighbors with a higher disk number
				{
				*DNlptr++ = ii;			// enter a neighbor onto the list!
				++cnt;
				}
			}
		}

	if(iy != (nrows-1))				// if not on highest row, examine upper three cells
		{
		if((ix%rowl) != 0)			// if not all the way to the left, examine upper left
			{
			index = (iy+1)*rowl + ix - 1;
			index *= CSIZE;
			nn = Cells[index++];
			for(j=0;j<nn;j++)
				{
				ii = Cells[index++];
//				if(abs(ii) > 10000) printf("wha7? %ld\n",ii);
				if(i < ii)					// only record neighbors with a higher disk number
					{
					*DNlptr++ = ii;			// enter a neighbor onto the list!
					++cnt;
					}
				}
			}

		index = (iy+1)*rowl + ix;	// examine upper center
		index *= CSIZE;
		nn = Cells[index++];
		for(j=0;j<nn;j++)
			{
			ii = Cells[index++];
//			if(abs(ii) > 10000) printf("wha8? %ld\n",ii);
			if(i < ii)					// only record neighbors with a higher disk number
				{
				*DNlptr++ = ii;			// enter a neighbor onto the list!
				++cnt;
				}
			}
		
		if((ix%rowl) != (rowl-1))	// if not all the way to the right, examine upper right
			{
			index = (iy+1)*rowl + ix + 1;
			index *= CSIZE;
			nn = Cells[index++];
			for(j=0;j<nn;j++)
				{
				ii = Cells[index++];
//				if(abs(ii) > 10000) printf("wha9? %ld\n",ii);
				if(i < ii)					// only record neighbors with a higher disk number
					{
					*DNlptr++ = ii;			// enter a neighbor onto the list!
					++cnt;
					}
				}
			}
		}
#if CB1
	else						// we are on top row, check lowest row, circular boundary conditions
		{
		if((ix%rowl) != 0)			// if not all the way to the left, examine bottom row left
			{
			index = ix - 1;
			index *= CSIZE;
			nn = Cells[index++];
			for(j=0;j<nn;j++)
				{
				ii = Cells[index++];
				if(i < ii)					// only record neighbors with a higher disk number
					{
					*DNlptr++ = ii;			// enter a neighbor onto the list!
					++cnt;
					}
				}
			}

		index = ix;					// examine bottom row center
		index *= CSIZE;
		nn = Cells[index++];
		for(j=0;j<nn;j++)
			{
			ii = Cells[index++];
			if(i < ii)					// only record neighbors with a higher disk number
				{
				*DNlptr++ = ii;			// enter a neighbor onto the list!
				++cnt;
				}
			}
		
		if((ix%rowl) != (rowl-1))	// if not all the way to the right, examine bottom row right
			{
			index = ix + 1;
			index *= CSIZE;
			nn = Cells[index++];
			for(j=0;j<nn;j++)
				{
				ii = Cells[index++];
				if(i < ii)					// only record neighbors with a higher disk number
					{
					*DNlptr++ = ii;			// enter a neighbor onto the list!
					++cnt;
					}
				}
			}
		}
#endif
		*DNcptr++ = cnt;			// record total neighbor count
		DNcptr += cnt;				// increment count pointer past neighbor entries
	}
}

// this one runs through the list of disks in each cell,
// assembled by sortcells() and placed in Cells[],
// and makes a list of neighbors of each disk, which is
// put in DiskNbrs[].  the format is, number of neighbors,
// followed by a list of the neighbor indices, close packed.
// only neighbors with a higher index than the current disk
// are listed, as each pair only has to be checked once.
// this one uses circlular boundary conditions,
// in the y direction YW2 wraps around to YW1 and vice versa, and
// in the x direction XW2 wraps around to xpist1, and vice versa.
// we continue to use the cells esablished for the full box, but
// we'll attempt to wrap around at the x piston location...
// the collision detection check is done in bouncecheckCB()
void mkdisknbrlistCB()
{
long i,j,index,*DNcptr,*DNlptr,cnt,ix,iy,nn,ii,ixp;

ixp = (xpist1 - XW1) * oLS;			// x index of the first piston location...
//printf("%Lf %Lf %Lf %Le\n",xpist,xpist - XW1,(xpist - XW1) * oLS,(REAL)oLS-10.0);
DNcptr = DNlptr = DiskNbrs;

for(i=0;i<nspheres;i++)
	{
	ix = (X[i] - XW1) * oLS;
	iy = (Y[i] - YW1) * oLS;
	index = iy * rowl + ix;
	
	index *= CSIZE;
	nn = Cells[index++];			// center cell index count
	
	++DNlptr;						// jump past neighbor count entry
	cnt = 0;
	for(j=0;j<nn;j++)
		{
		ii = Cells[index++];
		if(i < ii)					// only record neighbors with a higher disk number
			{
			*DNlptr++ = ii;			// enter a neighbor onto the list!
			++cnt;					// increment neighbor count
			}
		}

// bottom three cells
	if(iy != 0)						// if not on lowest row, examine lower three cells
		{
		if((ix%rowl) != 0)			// if not all the way to the left, examine lower left
			{
			if(ix > ixp)			// also if not at piston cell...
				{
				index = (iy-1)*rowl + ix - 1;
				index *= CSIZE;
				nn = Cells[index++];
				for(j=0;j<nn;j++)
					{
					ii = Cells[index++];
					if(i < ii)					// only record neighbors with a higher disk number
						{
						*DNlptr++ = ii;			// enter a neighbor onto the list!
						++cnt;
						}
					}
				}
			else					// at piston wall cell! (look at right wall, one row down)
				{
				index = (iy-1)*rowl + rowl - 1;
				index *= CSIZE;
				nn = Cells[index++];
				for(j=0;j<nn;j++)
					{
					ii = Cells[index++];
					if(i < ii)					// only record neighbors with a higher disk number
						{
						*DNlptr++ = ii;			// enter a neighbor onto the list!
						++cnt;
						}
					}
				}
			}
		else						// all the way to the left, examine right wall
			{
			index = (iy-1)*rowl + rowl - 1;
			index *= CSIZE;
			nn = Cells[index++];
			for(j=0;j<nn;j++)
				{
				ii = Cells[index++];
				if(i < ii)					// only record neighbors with a higher disk number
					{
					*DNlptr++ = ii;			// enter a neighbor onto the list!
					++cnt;
					}
				}
			}

		index = (iy-1)*rowl + ix;	// examine lower center
		index *= CSIZE;
		nn = Cells[index++];
		for(j=0;j<nn;j++)
			{
			ii = Cells[index++];
			if(i < ii)					// only record neighbors with a higher disk number
				{
				*DNlptr++ = ii;			// enter a neighbor onto the list!
				++cnt;					// increment neighbor count
				}
			}
		
		if((ix%rowl) != (rowl-1))	// if not all the way to the right, examine lower right
			{
			index = (iy-1)*rowl + ix + 1;
			index *= CSIZE;
			nn = Cells[index++];
			for(j=0;j<nn;j++)
				{
				ii = Cells[index++];
				if(i < ii)					// only record neighbors with a higher disk number
					{
					*DNlptr++ = ii;			// enter a neighbor onto the list!
					++cnt;
					}
				}
			}
		else						// we are all the way to the right, wrap around to piston cell
			{
			index = (iy-1)*rowl +ixp;
			index *= CSIZE;
			nn = Cells[index++];
			for(j=0;j<nn;j++)
				{
				ii = Cells[index++];
				if(i < ii)					// only record neighbors with a higher disk number
					{
					*DNlptr++ = ii;			// enter a neighbor onto the list!
					++cnt;
					}
				}
			}
		}
#if CB1
	else					// we are on lowest row, check highest row (circular boundary conditions)
		{
		if((ix%rowl) != 0)			// if not all the way to the left, examine top row left
			{
			if(ix > ixp)			// also if not at piston cell...
				{
				index = (nrows-1)*rowl + ix - 1;
				index *= CSIZE;
				nn = Cells[index++];
				for(j=0;j<nn;j++)
					{
					ii = Cells[index++];
					if(i < ii)					// only record neighbors with a higher disk number
						{
						*DNlptr++ = ii;			// enter a neighbor onto the list!
						++cnt;
						}
					}
				}
			else					// at piston wall cell! (top row, x direction wraparound)
				{
				index = nrows*rowl - 1;
				index *= CSIZE;
				nn = Cells[index++];
				for(j=0;j<nn;j++)
					{
					ii = Cells[index++];
					if(i < ii)					// only record neighbors with a higher disk number
						{
						*DNlptr++ = ii;			// enter a neighbor onto the list!
						++cnt;
						}
					}
				}
			}
		else						// all the way to the left, examine right wall
			{
			index = nrows*rowl - 1;
			index *= CSIZE;
			nn = Cells[index++];
			for(j=0;j<nn;j++)
				{
				ii = Cells[index++];
				if(i < ii)					// only record neighbors with a higher disk number
					{
					*DNlptr++ = ii;			// enter a neighbor onto the list!
					++cnt;
					}
				}
			}

		index = (nrows-1)*rowl + ix;	// examine top row center
		index *= CSIZE;
		nn = Cells[index++];
		for(j=0;j<nn;j++)
			{
			ii = Cells[index++];
			if(i < ii)					// only record neighbors with a higher disk number
				{
				*DNlptr++ = ii;			// enter a neighbor onto the list!
				++cnt;
				}
			}
		
		if((ix%rowl) != (rowl-1))	// if not all the way to the right, examine top row right
			{
			index = (nrows-1)*rowl + ix + 1;
			index *= CSIZE;
			nn = Cells[index++];
			for(j=0;j<nn;j++)
				{
				ii = Cells[index++];
				if(i < ii)					// only record neighbors with a higher disk number
					{
					*DNlptr++ = ii;			// enter a neighbor onto the list!
					++cnt;
					}
				}
			}
		else						// all the way to the right, examine piston cell!
			{
			index = (nrows-1)*rowl + ixp;
			index *= CSIZE;
			nn = Cells[index++];
			for(j=0;j<nn;j++)
				{
				ii = Cells[index++];
				if(i < ii)					// only record neighbors with a higher disk number
					{
					*DNlptr++ = ii;			// enter a neighbor onto the list!
					++cnt;
					}
				}
			}
		}
#endif

// middle cells (center already taken care of)
	if((ix%rowl) != 0)			// if not all the way to the left, examine middle left
		{
		if(ix > ixp)			// also if not at piston cell...
			{
			index = iy*rowl + ix - 1;
			index *= CSIZE;
			nn = Cells[index++];
			for(j=0;j<nn;j++)
				{
				ii = Cells[index++];
				if(i < ii)					// only record neighbors with a higher disk number
					{
					*DNlptr++ = ii;			// enter a neighbor onto the list!
					++cnt;
					}
				}
			}
		else					// at piston cell, examine right wall
			{
			index = iy*rowl + rowl - 1;
			index *= CSIZE;
			nn = Cells[index++];
			for(j=0;j<nn;j++)
				{
				ii = Cells[index++];
				if(i < ii)					// only record neighbors with a higher disk number
					{
					*DNlptr++ = ii;			// enter a neighbor onto the list!
					++cnt;
					}
				}
			}
		}
	else						// all the way to the left, examine right wall
		{
		index = iy*rowl + rowl - 1;
		index *= CSIZE;
		nn = Cells[index++];
		for(j=0;j<nn;j++)
			{
			ii = Cells[index++];
			if(i < ii)					// only record neighbors with a higher disk number
				{
				*DNlptr++ = ii;			// enter a neighbor onto the list!
				++cnt;
				}
			}
		}
	if((ix%rowl) != (rowl-1))	// if not all the way to the right, examine middle right
		{
		index = iy*rowl + ix + 1;
		index *= CSIZE;
		nn = Cells[index++];
		for(j=0;j<nn;j++)
			{
			ii = Cells[index++];
			if(i < ii)					// only record neighbors with a higher disk number
				{
				*DNlptr++ = ii;			// enter a neighbor onto the list!
				++cnt;
				}
			}
		}
	else						// all the way to the right, examine piston cell!
		{
		index = iy*rowl + ixp;
		index *= CSIZE;
		nn = Cells[index++];
		for(j=0;j<nn;j++)
			{
			ii = Cells[index++];
			if(i < ii)					// only record neighbors with a higher disk number
				{
				*DNlptr++ = ii;			// enter a neighbor onto the list!
				++cnt;
				}
			}
		}







// top three cells
	if(iy != (nrows-1))				// if not on highest row, examine upper three cells
		{
		if((ix%rowl) != 0)			// if not all the way to the left, examine upper left
			{
			if(ix > ixp)			// also if not at piston cell...
				{
				index = (iy+1)*rowl + ix - 1;
				index *= CSIZE;
				nn = Cells[index++];
				for(j=0;j<nn;j++)
					{
					ii = Cells[index++];
					if(i < ii)					// only record neighbors with a higher disk number
						{
						*DNlptr++ = ii;			// enter a neighbor onto the list!
						++cnt;
						}
					}
				}
			else					// at piston wall cell! (look at right wall, one row up)
				{
				index = (iy+1)*rowl + rowl - 1;
				index *= CSIZE;
				nn = Cells[index++];
				for(j=0;j<nn;j++)
					{
					ii = Cells[index++];
					if(i < ii)					// only record neighbors with a higher disk number
						{
						*DNlptr++ = ii;			// enter a neighbor onto the list!
						++cnt;
						}
					}
				}
			}
		else						// all the way to the left, examine right wall
			{
			index = (iy+1)*rowl + rowl - 1;
			index *= CSIZE;
			nn = Cells[index++];
			for(j=0;j<nn;j++)
				{
				ii = Cells[index++];
				if(i < ii)					// only record neighbors with a higher disk number
					{
					*DNlptr++ = ii;			// enter a neighbor onto the list!
					++cnt;
					}
				}
			}

		index = (iy+1)*rowl + ix;	// examine upper center
		index *= CSIZE;
		nn = Cells[index++];
		for(j=0;j<nn;j++)
			{
			ii = Cells[index++];
			if(i < ii)					// only record neighbors with a higher disk number
				{
				*DNlptr++ = ii;			// enter a neighbor onto the list!
				++cnt;
				}
			}



		
		if((ix%rowl) != (rowl-1))	// if not all the way to the right, examine upper right
			{
			index = (iy+1)*rowl + ix + 1;
			index *= CSIZE;
			nn = Cells[index++];
			for(j=0;j<nn;j++)
				{
				ii = Cells[index++];
				if(i < ii)					// only record neighbors with a higher disk number
					{
					*DNlptr++ = ii;			// enter a neighbor onto the list!
					++cnt;
					}
				}
			}
		else						// we are all the way to the right, wrap around to piston cell
			{
			index = (iy+1)*rowl +ixp;
			index *= CSIZE;
			nn = Cells[index++];
			for(j=0;j<nn;j++)
				{
				ii = Cells[index++];
				if(i < ii)					// only record neighbors with a higher disk number
					{
					*DNlptr++ = ii;			// enter a neighbor onto the list!
					++cnt;
					}
				}
			}
		}
#if CB1
	else						// we are on top row, check lowest row, circular boundary conditions
		{
		if((ix%rowl) != 0)			// if not all the way to the left, examine bottom row left
			{
			if(ix > ixp)			// also if not at piston cell...
				{
				index = ix - 1;
				index *= CSIZE;
				nn = Cells[index++];
				for(j=0;j<nn;j++)
					{
					ii = Cells[index++];
					if(i < ii)					// only record neighbors with a higher disk number
						{
						*DNlptr++ = ii;			// enter a neighbor onto the list!
						++cnt;
						}
					}
				}
			else					// at piston wall cell! (bottom row, x direction wraparound)
				{
				index = rowl - 1;
				index *= CSIZE;
				nn = Cells[index++];
				for(j=0;j<nn;j++)
					{
					ii = Cells[index++];
					if(i < ii)					// only record neighbors with a higher disk number
						{
						*DNlptr++ = ii;			// enter a neighbor onto the list!
						++cnt;
						}
					}
				}
			}
		else						// all the way to the left, examine right wall
			{
			index = rowl - 1;
			index *= CSIZE;
			nn = Cells[index++];
			for(j=0;j<nn;j++)
				{
				ii = Cells[index++];
				if(i < ii)					// only record neighbors with a higher disk number
					{
					*DNlptr++ = ii;			// enter a neighbor onto the list!
					++cnt;
					}
				}
			}

		index = ix;					// examine bottom row center
		index *= CSIZE;
		nn = Cells[index++];
		for(j=0;j<nn;j++)
			{
			ii = Cells[index++];
			if(i < ii)					// only record neighbors with a higher disk number
				{
				*DNlptr++ = ii;			// enter a neighbor onto the list!
				++cnt;
				}
			}

		if((ix%rowl) != (rowl-1))	// if not all the way to the right, examine bottom row right
			{
			index = ix + 1;
			index *= CSIZE;
			nn = Cells[index++];
			for(j=0;j<nn;j++)
				{
				ii = Cells[index++];
				if(i < ii)					// only record neighbors with a higher disk number
					{
					*DNlptr++ = ii;			// enter a neighbor onto the list!
					++cnt;
					}
				}
			}
		else						// all the way to the right, examine piston cell!
			{
			index = ixp;
			index *= CSIZE;
			nn = Cells[index++];
			for(j=0;j<nn;j++)
				{
				ii = Cells[index++];
				if(i < ii)					// only record neighbors with a higher disk number
					{
					*DNlptr++ = ii;			// enter a neighbor onto the list!
					++cnt;
					}
				}
			}
		}
#endif
		*DNcptr++ = cnt;			// record total neighbor count
		DNcptr += cnt;				// increment count pointer past neighbor entries
	}
}

// check if balls overlapping, if so, put on bump list...
// variable radius version...
void bouncechecka(REAL *XP,REAL *YP)
{
long i,j;

for(i=0;i<nspheres-1;i++)
	{
	register REAL x0,y0,x,y,delx,dely;
	//register REAL vcmx,vcmy;
	register REAL R2;
	register REAL r1,r12sq;
	
	x0 = XP[i];
	y0 = YP[i];
	//vx0 = VX[i];
	//vy0 = VY[i];
	r1 = Radius[i];
	for(j=i+1;j<nspheres;j++)
		{
		x = XP[j];
		y = YP[j];
	
		//vx = VX[j];
		//vy = VY[j];
		r12sq = r1 + Radius[j];
		r12sq *= r12sq;		
		delx = x - x0;
		dely = y - y0;
		R2 = delx*delx + dely*dely;		
		if(R2 > r12sq) continue;
		bumps[bumpnum].disknum = i; // first disk number...
		bumps[bumpnum].iden = j;	// second disk number...
		++bumpnum;
		++bumptot;
//		if(r12sq-R2 > squish) ++overlap;
		}

	}
}

// check if balls overlapping, if so, put on bump list...
// variable radius version...
// checks vertical wraparound if CB1 defined
void bouncecheckaa(REAL *XP,REAL *YP)
{
long i,j;

for(i=0;i<nspheres-1;i++)
	{
	register REAL x1,y1,x2,y2,delx,dely;
	register REAL R2;
	register REAL r1,r12sq;
	
	x1 = XP[i];
	y1 = YP[i];
	r1 = Radius[i];
#if CB1		
	if((y1+2*r1 < YW2) && (y1-2*r1 > YW1))
		{
#endif
		for(j=i+1;j<nspheres;j++)
			{
			x2 = XP[j];
			y2 = YP[j];

			r12sq = r1 + Radius[j];
			r12sq *= r12sq;		
			delx = x2 - x1;
			dely = y2 - y1;
			R2 = delx*delx + dely*dely;		
			if(R2 > r12sq) continue;
			bumps[bumpnum].disknum = i; // first disk number...
			bumps[bumpnum].iden = j;	// second disk number...
			++bumpnum;
			++bumptot;
//			if(r12sq-R2 > squish) ++overlap;
			}
#if CB1
		}
	else				// close to upper or lower edge
		{
		for(j=i+1;j<nspheres;j++)
			{
			REAL R3;
			x2 = XP[j];
			y2 = YP[j];

			r12sq = r1 + Radius[j];
			r12sq *= r12sq;		
			delx = x2 - x1;
			dely = y2 - y1;
			R2 = delx*delx + dely*dely;
			dely += (YW2-YW1);
			R3 = delx*delx + dely*dely;
			if(R3 < R2) R2 = R3;
			dely -= 2*(YW2-YW1);
			R3 = delx*delx + dely*dely;
			if(R3 < R2) R2 = R3;

			if(R2 > r12sq) continue;
			bumps[bumpnum].disknum = i; // first disk number...
			bumps[bumpnum].iden = j;	// second disk number...
			++bumpnum;
			++bumptot;
//			if(r12sq-R2 > squish) ++overlap;
			}
		}
#endif
	}
}

// neighborhood version, uses DiskNbrs[]
// use with mkdisknbrlist()
// checks vertical wraparound if CB1 defined
void bouncecheckd()
{
long i,j,*DNptr,nnbrs;

DNptr = DiskNbrs;

for(i=0;i<nspheres;i++)
	{
	register REAL x1,y1,x2,y2,delx,dely;
	register REAL R2,r1,r12sq;

	nnbrs = *DNptr++;
	if(nnbrs)
		{
		x1 = X[i];
		y1 = Y[i];
		r1 = Radius[i];
#if CB1		
		if((y1+2*r1 < YW2) && (y1-2*r1 > YW1))
			{
#endif
			for(j=0;j<nnbrs;j++)
				{
				long jdex;
			
				jdex = *DNptr++;
				x2 = X[jdex];
				y2 = Y[jdex];
				r12sq = r1 + Radius[jdex];
				r12sq *= r12sq;
				delx = x2 - x1;
				dely = y2 - y1;
				R2 = delx*delx + dely*dely;
				if(R2 > r12sq)
					continue;					// no overlap
				bumps[bumpnum].disknum = i;		// first disk number...
				bumps[bumpnum].iden = jdex;		// second disk number...
				++bumpnum;
				++bumptot;
#if 0
				if(r12sq-R2 > squish)
					{
					Overlist[2*overlap] = i;
					Overlist[2*overlap+1] = jdex;
					++overlap;
					}
#endif
				}
#if CB1
			}
		else				// close to upper or lower edge
			{
			for(j=0;j<nnbrs;j++)
				{
				long jdex;
				REAL R3;
			
				jdex = *DNptr++;
				x2 = X[jdex];
				y2 = Y[jdex];
				r12sq = r1 + Radius[jdex];
				r12sq *= r12sq;
				delx = x2 - x1;
				dely = y2 - y1;
				R2 = delx*delx + dely*dely;
				
				dely += (YW2-YW1);
				R3 = delx*delx + dely*dely;
				if(R3 < R2) R2 = R3;
				dely -= 2*(YW2-YW1);
				R3 = delx*delx + dely*dely;
				if(R3 < R2) R2 = R3;
			
				if(R2 > r12sq)
					continue;					// no overlap
				bumps[bumpnum].disknum = i;		// first disk number...
				bumps[bumpnum].iden = jdex;		// second disk number...
				++bumpnum;
				++bumptot;
//				if(r12sq-R2 > squish) ++overlap;
				}
			}
#endif
		}
	}
}

// neighborhood version, uses DiskNbrs[]
// use with mkdisknbrlistCB()
// checks vertical wraparound if CB1 defined
// horizontal w.r.t. XW2 and piston position xpist1
void bouncecheckCB()
{
long i,j,*DNptr,nnbrs;

DNptr = DiskNbrs;

for(i=0;i<nspheres;i++)
	{
	register REAL x1,y1,x2,y2,delx,dely;
	register REAL R2,r1,r12sq;

	nnbrs = *DNptr++;
	if(nnbrs)
		{
		x1 = X[i];
		y1 = Y[i];
		r1 = Radius[i];
#if CB1		
		if((y1+2*r1 < YW2) && (y1-2*r1 > YW1) && (x1-2*r1 > xpist1) && (x1+2*r1 < XW2))
			{
#endif
			for(j=0;j<nnbrs;j++)
				{
				long jdex;
			
				jdex = *DNptr++;
				x2 = X[jdex];
				y2 = Y[jdex];
				r12sq = r1 + Radius[jdex];
				r12sq *= r12sq;
				delx = x2 - x1;
				dely = y2 - y1;
				R2 = delx*delx + dely*dely;
				if(R2 > r12sq)
					continue;					// no overlap
				bumps[bumpnum].disknum = i;		// first disk number...
				bumps[bumpnum].iden = jdex;		// second disk number...
				++bumpnum;
				++bumptot;
				}
#if CB1
			}
		else				// close to upper or lower edge
			{
			for(j=0;j<nnbrs;j++)
				{
				long jdex;
				REAL dely1;
			
				jdex = *DNptr++;
				x2 = X[jdex];
				y2 = Y[jdex];
				r12sq = r1 + Radius[jdex];
				r12sq *= r12sq;
				delx = x2 - x1;
				dely = y2 - y1;

				dely1 = dely + (YW2-YW1);
				if(fabs(dely1) < fabs(dely)) dely = dely1;		// compiler inserts fabsl()?
				dely1 = dely - (YW2-YW1);
				if(fabs(dely1) < fabs(dely)) dely = dely1;
#if CBNDYS
				{
				REAL delx1;
				delx1 = delx + (XW2-xpist1);
				if(fabs(delx1) < fabs(delx)) delx = delx1;
				delx1 = delx - (XW2-xpist1);
				if(fabs(delx1) < fabs(delx)) delx = delx1;
				}
#endif
				R2 = delx*delx + dely*dely;
//				printf("%f %f %f\n",delx,dely,R2);




//				R3 = delx*delx + dely1*dely1;
//				if(R3 < R2) R2 = R3;
//				dely1 = dely - (YW2-YW1);
//				R3 = delx*delx + dely1*dely1;
//				if(R3 < R2) R2 = R3;
//				delx += (XW2-xpist);
//				R3 = delx*delx + dely*dely;
//				if(R3 < R2) R2 = R3;
//				delx -= 2*(XW2-xpist);
//				R3 = delx*delx + dely*dely;
//				if(R3 < R2) R2 = R3;
			
				if(R2 > r12sq)
					continue;					// no overlap
				bumps[bumpnum].disknum = i;		// first disk number...
				bumps[bumpnum].iden = jdex;		// second disk number...
				++bumpnum;
				++bumptot;
//				if(r12sq-R2 > squish) ++overlap;
				}
			}
#endif
		}
	}
}

// add additional disks to simulation, restricted to the "backroom"
void addisks1(long nadd, REAL rad)
{
long i;
REAL theta,x,y;

for(i=nspheres;i<nspheres+nadd;i++)
	{
	x = f_rand();
	x *= (BACKROOM - 2*rad);
	X[i] = xw2 - rad - x;
	y = f_rand();
	y *= (yw2-yw1 - 2*rad);
	Y[i] = yw1 + rad + y;
	myrand();
	theta = ang_rand();
	VX[i] = M_SQRT2 * cos(theta);
	VY[i] = M_SQRT2 * sin(theta);
	Radius[i] = rad;
	}
nspheres += nadd;
if(rad == rad1) nrad1 += nadd;
if(rad == rad2) nrad2 += nadd;
}

// remove disks from simulation, if radius matches "rad".
void subdisks1(long nsub, REAL rad)
{
long i,n,nrem=0;

n = nspheres - 1;
for(i=n;i>=0;i--)
	{
	if(nrem == nsub) break;		// we're done!
	if(Radius[i] == rad)		// found one!
		{
		--nspheres;
		++nrem;
		if(i != nspheres)		// move down top entry
			{
			X[i] = X[nspheres];
			Y[i] = Y[nspheres];
			VX[i] = VX[nspheres];
			VY[i] = VY[nspheres];
			Radius[i] = Radius[nspheres];
			}
		}
	}

if(rad == rad1)
	{
	nrad1 -= nsub;
	if(nrad1 < 0) nrad1 = 0;
	}
if(rad == rad2)
	{
	nrad2 -= nsub;
	if(nrad2 < 0) nrad2 = 0;
	}
}

// returns kinetic energy of all spheres...
REAL kengy()
{
long i;
REAL total,V;

total = 0.0;
for(i=0;i<nspheres;i++)
	{
	V = VX[i]*VX[i] + VY[i]*VY[i];
	total += V;
	}
return(.5 * total);
}

#define HX1 00
//#define HX2 400
#define HY1 00
#define HY2 200
#define HBIN 1

void plthist(unsigned long H[],long nmax,long ymax)
{
long i,ix,iy;
REAL yscale;

yscale = (REAL)(HY2 - HY1)/(REAL)ymax;
color(WHITE);
//gl_rect(HX1,HY1,HX2,HY2);
brectf(HX1,HY1,HX1+nmax*HBIN,HY2);
color(BLACK);
//recti(HX1,HY1,HX2,HY2);
brect(HX1,HY1,HX1+nmax*HBIN,HY2);
ix = HX1;
for(i=0;i<nmax;i++)
	{
	iy = H[i]*yscale + .49;
	brect(ix,HY2,ix+HBIN,HY2-iy);
	ix += HBIN;
	}
}

// copy the first histogram to the second, clear the first,
// return the maximum bin entry.
unsigned long copyclearhist(unsigned long H1[],unsigned long H2[],long nents)
{
long i,nmax;

nmax = 0;
for(i=0;i<nents;i++)
	{
	H2[i] = H1[i];
	if(H1[i] > nmax) nmax = H1[i];
	H1[i] = 0;
	}
return nmax;
}

// copy the first histogram to the second, leave the first alone,
// return the maximum bin entry.
unsigned long copyhist(unsigned long H1[],unsigned long H2[],long nents)
{
long i;
unsigned long nmax;

nmax = 0;
for(i=0;i<nents;i++)
	{
	H2[i] = H1[i];
	if(H1[i] > nmax) nmax = H1[i];
	}
return nmax;
}

// load histogram, test 
void ldhist0(unsigned long H[])
{
long i;

for(i=0;i<ROWL;i++)
	H[i] = i;
}

// histogram of disks positions along x axis
void ldiskhist(unsigned long H[])
{
long n,ii;
REAL x,wid;

wid = (xw2-xw1)/ROWL;
for(n=0;n<nspheres;n++)
	{
	x = X[n];
	x -= xw1;
	ii = x/wid;
	++H[ii];
	}
}

// histogram of speeds (maxwell-boltzmann in equilibrium?)
void ldspeedhist(unsigned long H[])
{
long n,ii;
REAL vel,wid;

wid = vmax/ROWL;
for(n=0;n<nspheres;n++)
	{
	vel = VX[n]*VX[n] + VY[n]*VY[n];
	vel = sqrt(vel);
	ii = vel/wid;
//	printf("ii %d\n",ii);
	if(ii >= NENTS) ii = NENTS-1;
	++H[ii];
	}
}

// histogram of wall profile, wall at x = 0
void ldwallhist(unsigned long H[])
{
long n,ii;
REAL x,wid;

//wid = (xw2-xw1)/ROWL;
wid = (xw2-xpist1)/ROWL;
for(n=0;n<nspheres;n++)
	{
	x = X[n];
	ii = x/wid;
	++H[ii];
	}
}

// compute virial approximation for hard-disk pressure
// Santos et al, J Chem Phys 103, 4622 (1995)			
REAL henderson(REAL d)
{
REAL alphaH,pi=M_PI;

alphaH = 7.0/3.0 - 4*sqrt(3.0)/pi;	// constant for virial pressure calculation
//printf("%f\n",alphaH);
return (1 + alphaH*d*d)/((1-d)*(1-d));
}

// compute approximate pressure in the chamber
// using henderson expression
REAL hender()
{
REAL volume,density,pi=M_PI,frac,vfac;

volume = (xw2 - xpist1) * (yw2 - yw1);
density = nspheres/volume;		// number density
frac = nrad1*pi*rad1*rad1/volume + nrad2*pi*rad2*rad2/volume;
vfac = henderson(frac);
//return (vfac * density);
return vfac ;
}

// change radii that match the first argument to the second argument
void changerad(REAL orad, REAL nrad)
{
long i;

for(i=0;i<nspheres;i++)
	{
	if(Radius[i] == orad)
		Radius[i] = nrad;
	}
}

// reset energies so that all disks have 1/2 kT
void reseteng()
{
long i;
REAL theta,sum=0,Ktot,cfac;

// printf("reseteng! %ld\n",nspheres);

for(i=0;i<nspheres;i++)
	{
	theta = ang_rand();
	VX[i] = M_SQRT2 * cos(theta);	// correct normalization for 1/2 kT per dof, i think
	VY[i] = M_SQRT2 * sin(theta);	// M_SQRT2 is in math.h
//	VX1[i] = M_SQRT2 * cos(theta);	// correct normalization for 1/2 kT per dof, i think
//	VY1[i] = M_SQRT2 * sin(theta);	// M_SQRT2 is in math.h
	sum += VX[i]*VX[i] + VY[i]*VY[i];
	}
//printf("energy sum %Lf nspheres %ld\n",.5*sum,nspheres);
//Ktot = kengy();
//printf("Ktot %Lf\n",Ktot);

// now remove translational momenta...
// assume masses all the same
sum = 0;	// fix x momenta
for(i=0;i<nspheres;i++)
	sum += VX[i];
sum /= nspheres;
for(i=0;i<nspheres;i++)
	VX[i] -= sum;
sum = 0;	// fix y momenta
for(i=0;i<nspheres;i++)
	sum += VY[i];
sum /= nspheres;
for(i=0;i<nspheres;i++)
	VY[i] -= sum;
Ktot = kengy();		// reset total energy
printf("Ktot %f\n",Ktot);
cfac = nspheres/Ktot;
cfac = sqrt(cfac);
printf("cfac %f\n",cfac);
for(i=0;i<nspheres;i++)
	{
	VX[i] *= cfac;
	VY[i] *= cfac;
	}

	
}

// conducts piston movement protocols...
// pistctl = 0	nothing happening
// pistctl = 1	initiate protocol
// pistctl = 2	protocol in action
void pproto()
{
//static REAL dt1=0;
//static int ctl = 0;
static REAL time,starttime;
if((pistctl == 1) && !ppctl)
	{
	starttime = total_time;
	vpist1 = lvel;
	ppctl = 1;
	pistctl = 2;
	}
if(ppctl & 1)			// move left piston
	{
	time = total_time - starttime;
	xpist1 = time * lvel + llo;
	if(xpist1 > lhi)	// end of left piston stroke
		{
		xpist1 = lhi;
		vpist1 = 0;
		// save KE
		KE1 = kengy();
		ppctl ^= 1;
		}
	}
if(pistctl == 2)	// move right piston
	{
	time = total_time - starttime;
	if(time >= pdelay)
		{
		vpist2 = rvel;
		ppctl |= 2;
		pistctl = 0;
		}
	}
if(ppctl & 2)
	{
	if(xpist2 > rhi)	// end of right piston stroke
		{
		xpist2 = rhi;
		vpist2 = 0;
		// save KE
		KE2 = kengy();
		ppctl ^= 2;
		}
	}
}
