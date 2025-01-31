// vbase, basic video manipulation program
// digitizer dimensions for isight, (320 x 240), (352 x 288), (640 x 480)
// digitizer dimensions for sony camcorder, (720 x 480)
#define FLOATRES 64
#if FLOATRES == 32
	#define REAL float
#elif FLOATRES == 80
	#define REAL long double
#else
	#define REAL double		// (FLOATRES 64)
#endif
//#define SIZEX	640		// texture and Vout[] dimensions
//#define SIZEY	480
#define SIZEX	960
#define SIZEY	480
#define XFAC	(REAL)SIZEX/SIZEY
#define MFAC 1.2			// display margin factor
// big box outline...
#define XW1		(REAL)-1.0*XFAC
#define XW2		1.0*XFAC
#define YW1		-1.0
#define YW2		1.0
// membrane outline...
#define ML -0.1
#define MR 0.1
#define DELT	.001		// simulation time step
#define NSPHERES 30000		// maximum number spheres...
#define RAD1	.04			// initial disk radii
#define RAD2	.03			// RAD1 is bigger if disks different sizes...
#define NENTS	4096		// max number of histogram entries...
#define ROWL	1000		// number histogram entries
#define SEED	101			// random number seed
#define NBEGIN 2000
#define MAXBUMPS 50000		// maximum size of "bumps" list
#define NSORT 10			// steps between cell sorts...
#define LS .2				// length of neighborhood edge... at least as big as biggest disk!
#define oLS (1.0/LS)		// to hopefully speed up computation (need parentheses!  can round off)
#define CSIZE 332			// cell data size, maximum neighbors plus one...
#define NCELLS 1024*4*4		// max number of neighborhoods (at least rowl*nrows)...
#define BACKROOM 1			// area on right hand side to replace 'recycled' disks
#define CB1		1			// circular boundary conditions in vertical direction
#define CBNDYS	0			// circular boundary conditions in both directions (piston is left boundary)
#define FCTL	0			// file control, 1 to write to file, 0 otherwise...
#define NRUNS	101			// number of runs in datafile...
#define DATAFILE	"walltet.dat"
#define REMARK		"wall density profile"
#define NTRANS	10000		// transient wait
#define DATAINT	100			// number of time steps between each histogram sweep
#define DATACNT	10000		// length of data run
#define NWAIT	1000		// wait between data points
#define REVSTEP .0001		// time reversal step
#define REVMAX 0.2			// max reversal time
#define NSTROKES 10000000	// number piston strokes in data run
#define VPIST 0.1			// piston velocity

// prototypes..
void Display();
void display();
void fadeout();
void fadein();
void restorecolor();
void initspheres9(long r1,long r2,long seed);
void initgauss(long nr1, long nr2, REAL r1, REAL r2);
void update4b(REAL *XP,REAL *YP,REAL *VXP,REAL *VYP,REAL delt);
void update5(REAL *XP,REAL *YP,REAL *VXP,REAL *VYP,REAL dtt);
void updateCB(REAL *XP,REAL *YP,REAL *VXP,REAL *VYP,REAL dtt);
void bigbox2(REAL *XP, REAL *YP);
void bigbox4(REAL *XP, REAL *YP);
void dobumpsa(REAL *XP, REAL *YP);
void dobouncea(REAL *XP,REAL *YP,long d1,long d2);
void dobounceCB(REAL *XP,REAL *YP,long d1,long d2);
void dobounceCBa(REAL *XP,REAL *YP,long d1,long d2);
void dofcbouncea(REAL *XP,REAL *YP,long d1,REAL xcen,REAL ycen);
void bouncechecka(REAL *XP,REAL *YP);
void bouncecheckaa(REAL *XP,REAL *YP);
void bouncecheckd();
void bouncecheckCB();
void bouncecheckCB();
void sortcells(void);
void mkdisknbrlist();
void mkdisknbrlistCB();
void addisks1(long nadd, REAL rad);
void subdisks1(long nsub, REAL rad);
void ldhist0(unsigned long H[]);
void ldiskhist(unsigned long H[]);
void ldspeedhist(unsigned long H[]);
void ldwallhist(unsigned long H[]);
void plthist(unsigned long H[],long nmax,long ymax);
unsigned long copyclearhist(unsigned long H1[],unsigned long H2[],long nents);
unsigned long copyhist(unsigned long H1[],unsigned long H2[],long nents);
void changerad(REAL orad, REAL nrad);
void initfile(void);
void writeheader();
void writedata();
void loadpath(char *arg, char *name);
REAL kengy();
REAL hender();
void gaussrand(double *g1, double *g2);
void reseteng();
void pproto();

// globals..
extern long slowctl;
extern long keyslowctl;			// changed in kbd.c
extern long Slowctl;			// screen refresh slow display...
extern int dispctl;				// display control
extern int Dispctl;				// second display control
extern long nspheres;			// total number of disks
extern long nrad1,nrad2;		// number of disks with RAD1, RAD2...
extern long seed;
extern long ncount;				// iteration counter
extern REAL xpist1;				// left piston x position
extern REAL xpist2;				// right piston x position
extern REAL vpist1;				// left piston x velocity
extern REAL vpist2;				// right piston x velocity
extern REAL dxpist;				// piston x step per display
extern REAL rad1,rad2;			// disk radii
extern long nrad1,nrad2;		// number of disks with RAD1, RAD2...
extern long histctl;			// display histograms
extern unsigned long Histo[],Histo1[],Histo2[],Histowall[],HistoW[];
extern long trailctl;			// trail control
extern char *mypath;			// program path
extern char pathName[];			// path name constructed with loadpath()
extern REAL X[],Y[],VX[],VY[],Radius[];
extern REAL binwidth;			// histogram bin width
extern long nbins;				// number of histogram bins used
extern long wallbump;
extern int forward;				// flag for trajectory direction
extern int pistctl;
extern REAL llo,lhi,rlo,rhi;	// left, right piston starting, ending positions
extern REAL lvel,rvel;			// left, right piston velocities
extern int radctl;				// for adjusting radii
extern int modectl;				// demo mode...
