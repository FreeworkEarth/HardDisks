// hsfast2.c - hard disks, SDL2 and wasm
// include histogramming...
#include <stdio.h>
#include <stdlib.h>
#include "macsdl2.h"
#include "myrand4a.h"
#if __EMSCRIPTEN__
#include <emscripten/emscripten.h>
#endif
#include <SDL2/SDL.h>

#define REAL float
#define SIZEX 1024
#define SIZEY 512
#define XFAC (REAL)SIZEX/SIZEY
#define NSPHERES 20000	// maximum number disks...
#define DELT .002		// time step
#define CB1	0			// vertical wraparound
#define CB2 0			// horizontal wraparound
#define NSORT 10		// steps between cell sorts...
#define CSIZE 128		// cell data size, maximum neighbors plus one...
#define NCELLS 1024*4*4	// max number of neighborhoods (at least rowl*nrows)...
#define LS .1			// length of neighborhood edge... at least as big as biggest disk!
#define NENTS	4096	// max number of histogram entries...
#define ROWL	100		// number histogram entries
#define NTRANS	100		// transient length...

// big box outline...
#define XW1 -1.0*XFAC
#define XW2 1.0*XFAC
#define YW1 -1.0
#define YW2 1.0

#define MAXBUMPS 100000

REAL X[NSPHERES],Y[NSPHERES],VX[NSPHERES],VY[NSPHERES];
REAL Radius[NSPHERES];
unsigned int Histo[NENTS],Histo1[NENTS],Histo2[NENTS];
unsigned int copyclearhist(unsigned int H1[],unsigned int H2[],int nents);

int Cells[NCELLS * CSIZE];
int DiskNbrs[NSPHERES * CSIZE];
int rowl,nrows,ncells;	// neighborhood indices
REAL LSx,LSy;			// neighborhood size
REAL oLSx,oLSy;			// inverse of above, to speed calculations

static int quit = 0;
int parity = 0;
int ncount=0,nc=0,ntick=0;
REAL x,y,theta=0,inc=.05;
REAL rad = .1;
//REAL rad1 = .03;
REAL rad1 = .01;
int nspheres = 5000;	// number of disks, must be smaller than NSPHERES
REAL xw1,xw2,yw1,yw2;	// box sides
REAL xscale,yscale;		// display scaling parameters
int fullctl = 0;		// full screen control
int slowctl = 1;		// slow display flag...
int histctl = 0;		// display histograms
int memctl = 0;			// membrane control
int idealctl = 0;		// ideal gas control
REAL vmax = 3.0;		// maximum for velocity histogram
unsigned int _seed = 1;	// for random number generator...
int framecount = 0;
int pnum,pnumax = 7;	// for presets
int ntrans = 1000;		// for transient control

struct bump{
	int disknum;
	int iden;	// positive for second disk, negative for walls, force centers...
	REAL fcenx;	// force center coordinates...
	REAL fceny;
};
struct bump bumps[MAXBUMPS];
int bumpnum = 0;

void initspheres(int nr1,float rad1,int seed);
void update();
void bigbox();
void dobumps();
void bouncechecka();
void dobounce(int d1,int d2);
void sortcells(void);
void mkdisknbrlist();
void bouncecheckd();
extern void brect(int ix1,int iy1,int ix2,int iy2);
extern void brectf(int x1,int y1,int x2,int y2);
extern void brcircf(long x0,long y0,long rad);
void ldhist0(unsigned int H[]);
void ldiskhist(unsigned int H[]);
void ldiskhist_r(unsigned int H[]);
void ldspeedhist(unsigned int H[]);
void plthist(unsigned int H[],int nmax,int ymax);
unsigned int copyhist(unsigned int H1[],unsigned int H2[],int nents);
void gaussrand(double *g1, double *g2);
REAL kengy();
void preset(int);		// parameter presets

void display()
{
int i;
REAL ymax;
char Cbuf[256];

color(WHITE);
clear();
/*
color(BLUE);
circf(x,y,rad);
color(BLACK);
circ(x,y,rad);
*/
for(i=0;i<nspheres/2;i++)
	{
	color(RED);
	circf(X[i],Y[i],Radius[i]);
	color(BLACK);
	circ(X[i],Y[i],Radius[i]);
	}

for(i=nspheres/2;i<nspheres;i++)
	{
	color(GREEN);
	circf(X[i],Y[i],Radius[i]);
	color(BLACK);
	circ(X[i],Y[i],Radius[i]);
	}

color(BLACK);		// draw box
move(xw1,yw1);
pendn();
move(xw1,yw2);
move(xw2,yw2);
move(xw2,yw1);
move(xw1,yw1);
penup();

if(memctl)			// draw membrane
	{
	color(BLACK);
	move(0,yw1); pendn();
	move(0,yw2); penup();
	}

if(histctl)			// histograms
	{
//	ldhist0(Histo);
//	plthist(Histo,200,200);
	if(histctl==1)	ymax = copyhist(Histo1,Histo,ROWL);
	if(histctl==2)	ymax = copyclearhist(Histo2,Histo,ROWL);
	if(histctl==3)  ymax = copyhist(Histo1,Histo,ROWL);
	if(histctl==4)  ymax = copyclearhist(Histo2,Histo,ROWL);
	if(histctl==5)  ymax = copyclearhist(Histo2,Histo,ROWL);	// for red disks
//	if(histctl==3)  ymax = copyclearhist(Histo2,Histo,ROWL);
	plthist(Histo,ROWL,(5*(ymax/4)));
	}
// p(v) = v * exp(-x*x/2)
// maximum at v = 1, p(1) = exp(-1/2) = 0.6065306597126334
if(histctl == 3 || histctl == 4)	// draw maxwell distribution
	{
	float v,y,iy,yscale;
	yscale = (REAL)(400)*0.6065306597126334;
//	plthist(Histo,ROWL,(5*(ymax/4)));
	color(RED);
	
//	for(i=0;i<SIZEX;i++)
	for(i=0;i<4*ROWL;i++)
		{
		v = (i*vmax)/(float)(2.8*ROWL);
		y = v * exp(-(v*v)/2);
//		if(i==SIZEX-10)
//			printf("v %f  y %f\n",v,y);
		iy = -260*y + 200;
//		iy = -yscale*y + 200;
		brcircf(i,iy,1);
		}
	}

drawtexture();
//fontet();
drawnum(ncount,820,0);
//iPrintf(10,0,"%d disks",nspheres);
sprintf(Cbuf,"%d disks",nspheres);
label(Cbuf,10,0);
//label("disks",0,0);
//drawnum(nspheres,100,0);
//printf("boop!\n");
drawbuffer();
++framecount;
}

void main_tick() {
SDL_Event event;

while (SDL_PollEvent(&event)) {
	switch (event.type) {
		case SDL_MOUSEBUTTONDOWN:
		case SDL_FINGERDOWN: {
			memctl = 1;
			initspheres(nspheres,rad1,12);
			ncount = 0;
//			ntrans = 1000;
			for(int i=0;i<NENTS;i++) Histo[i] = Histo1[i] = Histo2[i] = 0;
		}
		case SDL_KEYDOWN: {
			switch (event.key.keysym.sym) {
				case SDLK_x: {
					quit = 1;
					break; }
				case SDLK_ESCAPE: {
					quit = 1;
					break; }
				case SDLK_n: {
					memctl = 1;
					initspheres(nspheres,rad1,12);
					ncount = 0;
					for(int i=0;i<NENTS;i++) Histo[i] = Histo1[i] = Histo2[i] = 0;
					if(SDL_GetModState() & KMOD_SHIFT) ntrans = -1;	// turn off auto membrane removal
					break; }
				case SDLK_f: {
					fullctl += 1;
					if(fullctl > 1) fullctl = 0;
					printf("fullctl %d\n",fullctl);
					fullscreen(fullctl);
					break; }
				case SDLK_h: {
					if(SDL_GetModState() & KMOD_SHIFT)	// clear histograms
						{
						int i;
						for(i=0;i<NENTS;i++) Histo1[i] = Histo2[i] = 0;
						break;
						}
					else
						{
						++histctl;
						if(histctl > 5) histctl = 0;
						for(int i=0;i<NENTS;i++) Histo1[i] = Histo2[i] = 0;
						printf("histctl %d\n",histctl);
						}
					break; }
				case SDLK_i: {
					idealctl ^= 1;
					printf("idealctl %d\n",idealctl);
					if(idealctl)
						initspheres(100000,rad1,12);
					else
						initspheres(5000,rad1,12);
					break; }
				case SDLK_k: {
					REAL KE=0;
					KE = kengy();
					printf("kinetic energy %f\n",KE);
					break; }
				case SDLK_m: {
					memctl ^= 1;
					printf("memctl %d\n",memctl);
					break; }
				case SDLK_o: {
					float dt;
					static int oldncount = 0;
					static int oldticks = 0;
					int newticks;
					//extern int ncount;
					newticks = SDL_GetTicks();
					dt = newticks - oldticks;
					dt *= .001;
					oldticks = newticks;
					printf("elapsed time %.2f\t",dt);
					printf("frames per second %.2f\n",framecount/dt);
					framecount = 0;
					printf("iterations per second %.2f\n",(ncount - oldncount)/dt);
					oldncount = ncount;
					break; }
				case SDLK_s: {
					if(SDL_GetModState() & KMOD_SHIFT)
						{
						if(SDL_GetModState() & KMOD_ALT) slowctl /= 2;
						else slowctl /= 10;
						if(slowctl == 0) slowctl = 1;
						}
					else
						{
						if(SDL_GetModState() & KMOD_ALT) slowctl *= 2;
						else slowctl *= 10;
						if(slowctl > 10000) slowctl = 1;
						}
					printf("slowctl %d\n",slowctl);
					break; }
				case SDLK_5: {
					++pnum;
					if(pnum > pnumax) pnum = 0;
					preset(pnum);
					break; }
				}
			}
		}
	}

for(int i = 0;i<slowctl;i++)
	{
update();
bumpnum = 0;					// this counts all collisions each time step
if(ncount == ntrans)
	memctl = 0;					// remove membrane after transient...
//bigbox();
if (!idealctl)
{
#if 1
	if(!(ncount%NSORT))	// sort disks into cells, neighborhood list...
		{
		sortcells();
		mkdisknbrlist();
		}
//		bouncechecke();			// some disks may have infinite mass
		bouncecheckd();			// neighborhood version, uses DiskNbrs[]
#else
	bouncechecka();
#endif
}	// end idealctl
bigbox();
dobumps();

if(histctl==1) ldiskhist(Histo1);		// accumulate histogram
if(histctl==2) ldiskhist(Histo2);		// single time step histogram
if(histctl==3) ldspeedhist(Histo1);		// accumulate velocity histogram
if(histctl==4) ldspeedhist(Histo2);		// accumulate velocity histogram (single step)
if(histctl==5) ldiskhist_r(Histo2);		// red disks single time step histogram
}

	display();
}

void main_loop()
{
#if __EMSCRIPTEN__
	    emscripten_set_main_loop(main_tick, -1, 1);
//		emscripten_set_main_loop(main_tick, 300, 1);
#else
    while (0 == quit){main_tick();}
#endif
}

int main()
{
rowl = (XW2 - XW1)/LS;
nrows = (YW2 - YW1)/LS;
ncells = rowl * nrows;
//printf("rowl %d  nrows %d  ncells %d\n",rowl,nrows,ncells);
for(int i=0;i<NENTS;i++) Histo[i] =  Histo1[i] = 0;
	
prefsize(1024,512);
initgraph("hi there!");
yscale = 1.2;
xscale = yscale * XFAC;
scale(-xscale,xscale,-yscale,yscale);
xw1 = XW1; xw2 = XW2;
yw1 = YW1; yw2 = YW2;
//printf("%f %f %f %f\n",xw1,xw2,yw1,yw2);
LSx = (xw2-xw1)/rowl;
LSy = (yw2-yw1)/nrows;
//printf("LSx %f  LSy %f\n",LSx,LSy);
oLSx = 1.0/LSx;
oLSy = 1.0/LSy;
initspheres(nspheres,rad1,12);
labelsize(30);
labelcolor(0,0,255,255);
initlabel();
main_loop();
exitgraph();
return 0;
}

void initspheres(int nr1,float rad1,int seed)		// number of rad1, radius, random seed...
{
int i;
REAL theta;
REAL xx,yy,delx,dely;
//register unsigned int _seed;

//_seed = seed;
for(i=0;i<17*seed;i++)
	myrand();

nspheres = nr1;

delx = xw2 - 2*rad1;
dely = yw2 - yw1 - 2*rad1;

for(i=0;i<nr1;i++)
	{
	xx = f_rand();
	xx /= 2.0;
	if(i < nr1/2)
		xx += .5;
	xx *= delx;
	yy = f_rand();
	yy *= dely;
	X[i] = xx + rad1;
	Y[i] = yw1 + yy + rad1;
	Radius[i] = rad1;
	theta = ang_rand();
	if(idealctl)
//	if(0)
		{
		double w1,w2;
		gaussrand(&w1,&w2);
		VX[i] = w1/M_SQRT2;
//		VY[i] = 0;
		VY[i] = w2/M_SQRT2;
		}
	else
		{
		VX[i] = 1.0 * cos(theta);
//		VY[i] = 0;
		VY[i] = 1.0 * sin(theta);
		}
	}
}

// time stepper...
// variable box size
void update()
{
int i;
register REAL dt;

dt = DELT;

for(i=0;i<nspheres;i++)
	{
	X[i] += dt*VX[i];
	Y[i] += dt*VY[i];
#if CB1
	if(Y[i] >= yw2)				// the equality case can happen!
		Y[i] -= (yw2-yw1);
	if(Y[i] < yw1)
		Y[i] += (yw2-yw1);
#endif
#if CB2
	if(X[i] >= xw2)				// the equality case can happen!
		X[i] -= (xw2-xw1);
	if(X[i] < xw1)
		X[i] += (xw2-yw1);
#endif
	}
x = sin(theta);
theta += inc;
y = 0;
	
++nc;
if(nc == 60)
	{
	++ntick;
//	printf("ntick %d\n",ntick);
	nc = 0;
	}
++ncount;
}

// puts boundary collisions in bumps structure...
// accumulate back wall "pressure" in VtotR and VtotL
// wall positions are variable (slows things a couple percent)
void bigbox()
{
int i,xxw1;
REAL rad;

if(memctl) xxw1 = 0;
else xxw1 = xw1;

for(i=0;i<nspheres;i++)
	{
	rad = Radius[i];
#if !CB2	
	if((X[i]-rad <= xxw1) && VX[i] < 0)
//	if((X[i]-rad <= XW1) && VX[i] < 0)
		{
		bumps[bumpnum].disknum = i;
		bumps[bumpnum].iden = -1;	// collision with left edge...
		++bumpnum;
//		++bumptot;
//		VtotL += 2 * VX[i];			// the bounce adds to the pressure!						
		}
	if((X[i]+rad >= xw2) && VX[i] > 0)
//	if((X[i]+rad >= XW2) && VX[i] > 0)
		{
		bumps[bumpnum].disknum = i;
		bumps[bumpnum].iden = -2;	// collision with right edge...
		++bumpnum;
//		++bumptot;
//		VtotR += 2 * VX[i];			// the bounce adds to the pressure!						
		}
#endif
#if !CB1
	if((Y[i]-rad <= yw1) && VY[i] < 0)
//	if((Y[i]-rad <= YW1) && VY[i] < 0)
		{
		bumps[bumpnum].disknum = i;
		bumps[bumpnum].iden = -3;	// collision with bottom edge...
		++bumpnum;
//		++bumptot;		
		}
	if((Y[i]+rad >= yw2) && VY[i] > 0)
//	if((Y[i]+rad >= YW2) && VY[i] > 0)
		{
		bumps[bumpnum].disknum = i;
		bumps[bumpnum].iden = -4;	// collision with top edge...
		++bumpnum;
//		++bumptot;
		}
#endif
	}
}

// do all collisions on bump list...
// variable radius version...
// type of collision is encoded in bumps[].iden
// from 0 up it's a collision with another disk with that index
// negative numbers, collisions with various obstacles
// -1 left edge, -2 right edge, -3 bottom edge, -4 top edge
// -5 right membrane edge, -6 left membrane edge
// -8 "force center", a point on the circumference of the disk
//    actually, force applied on the circumference where the line
//    from the force center to he disk center intersects it
void dobumps()
{
int i,dnum,iden;
REAL xcen,ycen;

for(i=0;i<bumpnum;i++)
	{
	dnum = bumps[i].disknum;
	iden = bumps[i].iden;
	if(iden >= 0)
		{
		dobounce(dnum,iden);
		}
	else
		{	
		switch(iden)
			{
			case -1:
			case -2:
			case -5:
			case -6:
				VX[dnum] = -VX[dnum];
				break;
			case -3:
			case -4:
				VY[dnum] = -VY[dnum];
				break;
			case -8:
				xcen = bumps[i].fcenx;
				ycen = bumps[i].fceny;
//				dofcbouncea(XP,YP,dnum,xcen,ycen);
				break;
			}
		}
	}
}

// check if balls overlapping, if so, put on bump list...
// variable radius version...
// checks each pair...
void bouncechecka()
{
int i,j;

for(i=0;i<nspheres-1;i++)
	{
	register REAL x0,y0,x,y,delx,dely;
	register REAL R2;
	register REAL r1,r12sq;
	
	x0 = X[i];
	y0 = Y[i];
	r1 = Radius[i];
	for(j=i+1;j<nspheres;j++)
		{
		x = X[j];
		y = Y[j];

		r12sq = r1 + Radius[j];
		r12sq *= r12sq;		
		delx = x - x0;
		dely = y - y0;
		R2 = delx*delx + dely*dely;		
		if(R2 > r12sq) continue;
		bumps[bumpnum].disknum = i; // first disk number...
		bumps[bumpnum].iden = j;	// second disk number...
		++bumpnum;
//		++bumptot;
//		if(r12sq-R2 > squish) ++overlap;
		}
	}
}

// perform bounce between argument disks...
// cleaned up a little...
// variable box size
void dobounce(int d1,int d2)
{
	register REAL x0,y0,x,y,delx,dely,vx,vy,vx0,vy0;
	register REAL vcmx,vcmy;
	register REAL R,sino,coso,R2;
	register REAL v1x,v1y,v2x,v2y,A;
	x0 = X[d1];
	y0 = Y[d1];
	vx0 = VX[d1];
	vy0 = VY[d1];
	
	x = X[d2];
	y = Y[d2];
	vx = VX[d2];
	vy = VY[d2];
	delx = x - x0;
	dely = y - y0;
	R2 = delx*delx + dely*dely;

#if CB1
		{
		REAL R3,dely1;
		dely1 = dely + (yw2-yw1);
		R3 = delx*delx + dely1*dely1;
		if(R3 < R2)
			{
			R2 = R3;
			dely = dely1;
			}
		dely1 = dely - (yw2-yw1);
		R3 = delx*delx + dely1*dely1;
		if(R3 < R2)
			{
			R2 = R3;
			dely = dely1;
			}
		}
#endif	

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

// perform bounce between argument disk and argument force center...
void dofcbouncea(REAL *XP,REAL *YP,int d1,REAL xcen,REAL ycen)
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
// don't use, so we can bounce against stationary disks
//if(R2 > rsq) return;
	
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
int i,ix,iy,index,nn;

for(i=0;i<ncells;i++)				// these are the neighbor counts
	Cells[CSIZE * i] = 0;
	
for(i=0;i<nspheres;i++)
	{
	ix = (X[i] - xw1) * oLSx;
	iy = (Y[i] - yw1) * oLSy;
	index = iy * rowl + ix;
	index *= CSIZE;
	
	nn = (Cells[index++])++;		// load and increment cell count...
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
int i,j,index,*DNcptr,*DNlptr,cnt,ix,iy,nn,ii;

DNcptr = DNlptr = DiskNbrs;
for(i=0;i<nspheres;i++)
	{
	ix = (X[i] - xw1) * oLSx;
	iy = (Y[i] - yw1) * oLSy;
	index = iy * rowl + ix;
	
	index *= CSIZE;
	nn = Cells[index++];				// center cell index count
	
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

// neighborhood version, uses DiskNbrs[]
// use with mkdisknbrlist()
// checks vertical wraparound if CB1 defined
// variable box size
void bouncecheckd()
{
int i,j,*DNptr,nnbrs;

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
		if((y1+2*r1 < yw2) && (y1-2*r1 > yw1))
			{
#endif
			for(j=0;j<nnbrs;j++)
				{
				int jdex;
			
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
	//			++bumptot;
	//			if(r12sq-R2 > squish)
	//				{
	//				Overlist[2*overlap] = i;
	//				Overlist[2*overlap+1] = jdex;
	//				++overlap;
	//				}
				}
#if CB1
			}
		else				// close to upper or lower edge
			{
			for(j=0;j<nnbrs;j++)
				{
				int jdex;
				REAL R3;
			
				jdex = *DNptr++;
				x2 = X[jdex];
				y2 = Y[jdex];
				r12sq = r1 + Radius[jdex];
				r12sq *= r12sq;
				delx = x2 - x1;
				dely = y2 - y1;
				R2 = delx*delx + dely*dely;
				
				dely += (yw2-yw1);
				R3 = delx*delx + dely*dely;
				if(R3 < R2) R2 = R3;
				dely -= 2*(yw2-yw1);
				R3 = delx*delx + dely*dely;
				if(R3 < R2) R2 = R3;
			
				if(R2 > r12sq)
					continue;					// no overlap
				bumps[bumpnum].disknum = i;		// first disk number...
				bumps[bumpnum].iden = jdex;		// second disk number...
				++bumpnum;
//				++bumptot;
//				if(r12sq-R2 > squish) ++overlap;
				}
			}
#endif
		}
	}
}

#define HX1 00
//#define HX2 400
#define HY1 00
#define HY2 200
#define HBIN 4

void plthist(unsigned int H[],int nmax,int ymax)
{
int i,ix,iy;
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
unsigned int copyclearhist(unsigned int H1[],unsigned int H2[],int nents)
{
int i,nmax;

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
unsigned int copyhist(unsigned int H1[],unsigned int H2[],int nents)
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
void ldhist0(unsigned int H[])
{
int i;

for(i=0;i<ROWL;i++)
	H[i] = i;
}

// histogram of disks positions along x axis
void ldiskhist(unsigned int H[])
{
int n,ii;
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

// histogram of disks positions along x axis
void ldiskhist_r(unsigned int H[])
{
int n,ii;
REAL x,wid;

wid = (xw2-xw1)/ROWL;
//for(n=0;n<nspheres/2;n++)
for(n=nspheres/2;n<nspheres;n++)
	{
	x = X[n];
	x -= xw1;
	ii = x/wid;
	++H[ii];
	}
}

// histogram of speeds (maxwell-boltzmann in equilibrium?)
void ldspeedhist(unsigned int H[])
{
int n,ii;
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

// sets up parameter presets...
#if __EMSCRIPTEN__
EMSCRIPTEN_KEEPALIVE 
#endif
void preset(int n)
{
printf("preset new? %d\n",n);
switch(n) {
	case 0:
		idealctl = 0;
		slowctl = 1;
		histctl = 0;
		ntrans = -1;
		nspheres = 5000;
		rad1 = .01;
		break;
	case 1:
		idealctl = 0;
		slowctl = 10;
		histctl = 2;
		ntrans = 500;
		memctl = 1;
		nspheres = 5000;
		rad1 = .01;
		break;
	case 2:
		idealctl = 0;
		slowctl = 10;
		histctl = 0;
		ntrans = 500;
		memctl = 1;
		nspheres = 5000;
		rad1 = .01;	
		break;
	case 3:
		idealctl = 0;
		slowctl = 10;
		histctl = 0;
		ntrans = 500;
		memctl = 1;
		nspheres = 3000;
		rad1 = .015;
		break;
	case 4:
		idealctl = 0;
		slowctl = 10;
		ntrans = 500;
		memctl = 1;
		nspheres = 3000;
		rad1 = .015;
		histctl = 2;
		break;
	case 5:
		idealctl = 0;
		slowctl = 10;
		ntrans = 500;
		memctl = 1;
		nspheres = 3000;
		rad1 = .015;
		histctl = 4;
		break;
	case 6:
		idealctl = 1;
		slowctl = 10;
		ntrans = 1000;
		memctl = 1;
		nspheres = 10000;
		rad1 = .015;
		histctl = 5;
		break;
	case 7:
		idealctl = 0;
		slowctl = 10;
		ntrans = 1000;
		memctl = 1;
		nspheres = 3000;
		rad1 = .015;
		histctl = 5;
		break;
	}
}
