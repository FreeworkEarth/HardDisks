// boxtest.c - example code, disks bouncing in a box
#include <stdio.h>
#include "macsdl2.h"
#include "myrand4a.h"
#if __EMSCRIPTEN__
#include <emscripten/emscripten.h>
#endif
#include <SDL2/SDL.h>


#define REAL float		// WebAssembly likes floats
#define SIZEX 1024		// size of display screen
#define SIZEY 512
#define XFAC (REAL)SIZEX/SIZEY
#define NDISKS 10000	// maximum number disks...
#define DELT .0001		// time step

REAL X[NDISKS],Y[NDISKS],VX[NDISKS],VY[NDISKS]; // position and velocity arrays for disks
REAL rightwall_X, topwall_Y, rightwall_VX, topwall_VY; // initialize right and top walls

static REAL pi = 3.141592653589;
static int quit = 0;
int ncount = 0;			// total number of updates
REAL xscale,yscale;		// display scaling parameters
unsigned int _seed = 1;	// for random number generator...
int framecount = 0;
int ndisks = 100;		// number of disks displayed
int slowctl = 10;		// number of updates per display
static REAL bx1= -1,bx2 = 1;	// box size
static REAL by1= -1,by2 = 1;
int heatbath = 0; // boolean value for heat on/off
int equilibrium_pressure = 0; //brings outside pressure to inside-ish

REAL energy;

//physical parameters
REAL radius = .05;		// radius of disks
REAL mass = 1;
REAL wallmass = 100;
REAL heatbathtemp = 0.66;

//histogram
REAL maxsize = 0;
#define numbins 24
int histogram_global[numbins];



//constants
REAL k_B = 1.38e-23;

REAL sample_gaussian(REAL mean, REAL dev) // uses box muller transform from wikipedia
{
	REAL rand1 = f_rand();
	REAL rand2 = f_rand();

	REAL r = sqrt(-2.0 * log(rand1));
    REAL theta = 2.0 * pi * rand2;

    REAL z1 = r * cos(theta);
    REAL z2 = r * sin(theta);

	return mean + dev * z2;
}

// REAL sample_boltzmann(REAL mean, REAL dev) // uses box muller transform from wikipedia
// {
// 	//
// 	break;
// }

REAL getVavg() {
	REAL vsum = 0;
	for(int i=0;i<ndisks;i++) {
		vsum += fabs(VY[i]);
	}
	return (vsum/ndisks);
}

REAL getVrms() {
	REAL vsquaredsum = 0;
	for(int i=0;i<ndisks;i++) {
		vsquaredsum += powf(VX[i],2) + powf(VY[i],2);
	}
	return (vsquaredsum/ndisks);
}

REAL getEnergy() {
	REAL vrms = getVrms();
	return 0.5*mass*vrms;
}

REAL getPistonFrequency() {
	REAL v_avg = getVavg();
	REAL frequency = ndisks*v_avg/(2*(topwall_Y-by1)); //per second
	REAL freq = frequency/DELT;
	return freq;
}

REAL getHeatBathFrequency() {
	REAL v_avg = getVavg();
	REAL bottom_frequency = ndisks*v_avg/(2*(topwall_Y-by1)); //per second
	REAL left_frequency = ndisks*v_avg/(2*(rightwall_X-bx1)); //per second
	REAL freq = (bottom_frequency+left_frequency)/DELT;
	return freq;
}

REAL getHeatBathTemp() {
	REAL v_avg = getVavg();
	REAL bottom_frequency = ndisks*v_avg/(2*(topwall_Y-by1)); //per second
	REAL left_frequency = ndisks*v_avg/(2*(rightwall_X-bx1)); //per second
	REAL freq = (bottom_frequency+left_frequency)/DELT;
	return freq;
}

void thermalize(int toggle, REAL temp){
	heatbath = toggle;
	heatbathtemp = temp;
}

void initialize(int nd,float rad)
{
for(int i=0;i<123;i++) myrand();	// scramble random number generator
for(int i=0;i<ndisks;i++)			// initial positions
	{
	float x,y; //declare x and y for a given particle
	x = f1_rand();		// random between -1 and 1
	x *= (1 - radius);	// keep entire disk inside box; -0.95 < x < 0.95
	X[i] = x; 		// sets x position for ith disk in array
	y = f1_rand();	// random between -1 and 1
	y *= (1 - radius);	// keep entire disk inside box; -0.95 < y < 0.95
	Y[i] = y;			// sets y position for ith disk in array
	}
for(int i=0;i<ndisks;i++)			// initial velocities
	{
//	float theta;
	//float theta = ang_rand();	// ang_rad() provides random angle from 0 to 2pi
	//VX[i] = 1.0 * cos(theta);	// x component of v, |v| = 1
	//VY[i] = 1.0 * sin(theta);	// y component of v, |v| = 1
	REAL random = f1_rand();
	int sign;
	if (random > 0) {
		sign = 1;
	} else {
		sign = -1;
	}

	VX[i] = f1_rand();
	VY[i] = f1_rand();
	//VY[i] = sign*sqrtf(1 - powf(VX[i],2));
	}

//init walls
	rightwall_X = bx2; 
	topwall_Y = by2;
}

void histogram(int coord, int nbins, REAL minvalue, REAL maxvalue) {
	REAL delta = (maxvalue-minvalue)/nbins; //hist width
	int bins[nbins];  // declares the buckets
	memset(bins, 0, sizeof(bins)); // makes the buckets empty
	memset(histogram_global, 0, sizeof(histogram_global)); // makes the buckets empty

	REAL values[ndisks]; // this is the array to be sampled from
	switch (coord)
	{
	case 0:
		memcpy(&values, &X, sizeof(values));
		break;
	case 1:
		memcpy(&values, &Y, sizeof(values));
		break;
	case 2:
		memcpy(&values, &VX, sizeof(values));
		break;
	case 3:
		memcpy(&values, &VY, sizeof(values));
		break;
	default:
		break;
	}

	for(int i=0;i<ndisks;i++) {
		for (float val=minvalue; val<maxvalue; val+=delta) {
			if (values[i] > val && values[i] < val + delta) {
				bins[(int)((val-minvalue)/delta)] += 1;
				histogram_global[(int)((val-minvalue)/delta)] += 1;
			}
		}
	}

	// for (int i = 0; i < nbins; i++) {
	// 	//
	// 	printf(" bin %d is %d \n", i, bins[i]);
	// }
}

REAL getmagnitude(REAL x, REAL y){
	return sqrtf(powf(x, 2) + powf(x, 2));
}

void update() // updates position using velocity and euler method
{

// for (int i = 0; i < ndisks; i++){
// 	printf("disk %d has velocity %f \n", i, VX[i]);
// }


register REAL dt=DELT; // stores time step size in cpu register

//update wall positions
topwall_Y += dt*topwall_VY; 
rightwall_X += dt*rightwall_VX; 
if(topwall_Y < (by1 + 0.4)) topwall_VY = 0;
if(rightwall_X < (bx1 + 0.4)) rightwall_VX = 0;
if(topwall_Y > by2) {
	topwall_Y = by2;
	topwall_VY = 0;
}
if(rightwall_X > bx2) rightwall_VX = 0;

for(int i=0;i<ndisks;i++)
	{
		
	// perform bounces...
	if(X[i] < bx1+radius) {
		X[i] = bx1+radius;	
		// REAL angle = atan(VY[i]/VX[i]);
		// angle += 0.05*f1_rand()*angle;

		// REAL magnitude = getmagnitude(VX[i], VY[i]);
		if (heatbath != 0) {
			VX[i] = sample_gaussian(heatbathtemp, 0.05*heatbathtemp);
		} else {
			VX[i] = -VX[i];
		}
		REAL original_mag = sqrtf(VX[i]*VX[i] + VY[i]*VY[i]);
		VX[i] += sample_gaussian(0, 0.1*original_mag);
		VY[i] += sample_gaussian(0, 0.1*original_mag);
		REAL new_mag = sqrtf(VX[i]*VX[i] + VY[i]*VY[i]);
		REAL scale = original_mag / new_mag;
    	VX[i] *= scale;
    	VY[i] *= scale;

	}
	if(X[i] > bx2-radius || X[i] > rightwall_X-radius) {
		X[i] = rightwall_X-radius;
		if (rightwall_VX != 0) {
			VX[i] -= 2*rightwall_VX;
		}
		VX[i] = -VX[i];	

		REAL original_mag = sqrtf(VX[i]*VX[i] + VY[i]*VY[i]);
		VX[i] += sample_gaussian(0, 0.1*original_mag);
		VY[i] += sample_gaussian(0, 0.1*original_mag);
		REAL new_mag = sqrtf(VX[i]*VX[i] + VY[i]*VY[i]);
		REAL scale = original_mag / new_mag;
    	VX[i] *= scale;
    	VY[i] *= scale;
	}
	if(Y[i] < by1+radius) {
		Y[i] = by1+radius;
		if (heatbath != 0) {
			VY[i] = sample_gaussian(heatbathtemp, 0.05*heatbathtemp);
		} else {
			VY[i] = -VY[i];
		}

		REAL original_mag = sqrtf(VX[i]*VX[i] + VY[i]*VY[i]);
		VX[i] += sample_gaussian(0, 0.1*original_mag);
		VY[i] += sample_gaussian(0, 0.1*original_mag);
		REAL new_mag = sqrtf(VX[i]*VX[i] + VY[i]*VY[i]);
		REAL scale = original_mag / new_mag;
    	VX[i] *= scale;
    	VY[i] *= scale;
		
	} 
	if(Y[i] > by2-radius || Y[i] > topwall_Y-radius) {
		Y[i] = topwall_Y-radius;
		 if (topwall_VY != 0) {
		 	VY[i] -= 2*topwall_VY;
		 }
		 VY[i] = -VY[i];

		REAL original_mag = sqrtf(VX[i]*VX[i] + VY[i]*VY[i]);
		VX[i] += sample_gaussian(0, 0.1*original_mag);
		VY[i] += sample_gaussian(0, 0.1*original_mag);
		REAL new_mag = sqrtf(VX[i]*VX[i] + VY[i]*VY[i]);
		REAL scale = original_mag / new_mag;
    	VX[i] *= scale;
    	VY[i] *= scale;


		
	}
	X[i] += dt*VX[i]; // increments ith particle's x position 
	Y[i] += dt*VY[i]; // increments ith particle's y position 
	}
++ncount; // counts how many time steps have elapsed 
}

void display()
{
color(WHITE);			// set brush color to white
clear();				// clears Plane1[] with selected color (white)
color(BLACK);			// set brush color to black
// draw surrounding box
move(bx1,by1); pendn(); move(bx2,by1); move(bx2,by2); move(bx1,by2); move(bx1,by1); penup();
if (heatbath) {
	color(RED);
	move(bx1,by1); pendn(); move(bx2,by1); penup();
	move(bx1,by1); pendn(); move(bx1,by2); penup();
}
// draw disks
for(int i=0;i<ndisks;i++)
	{
	color(WHITE);
	circf(X[i],Y[i],radius);
	color(BLACK);
	circ(X[i],Y[i],radius);
	}
// draw lines
move(rightwall_X,by1); pendn(); move(rightwall_X,by2); penup();
move(bx1,topwall_Y); pendn(); move(bx2,topwall_Y); penup();


histogram(0, numbins, -1.1, 1.1);
for (int i = 0; i < numbins; i++) {
		//
		if (- 1 + 1.5*(histogram_global[i]/ndisks) > 1.1) {
			printf("histogram bar too high");
		} else {
			move(1.5+ i*0.03, -1); pendn(); 
			move(1.5 + i*0.03, - 1 + 1.5*(histogram_global[i])/ndisks); 
			move(1.5 + (i+1)*0.03, -1 + 1.5*(histogram_global[i])/ndisks); 
			move(1.5 + (i+1)*0.03, -1); penup();
		}
	}
	
histogram(1, numbins, -1.1, 1.1);
for (int i = 0; i < numbins; i++) {
		//
		if (- 1 + 1.5*(histogram_global[i]/ndisks) > 1.1) {
			printf("histogram bar too high");
		} else {
			move(1.5+ i*0.03, -0.42); pendn(); 
			move(1.5 + i*0.03, - 0.42 + 1.5*(histogram_global[i])/ndisks); 
			move(1.5 + (i+1)*0.03, - 0.42 + 1.5*(histogram_global[i])/ndisks); 
			move(1.5 + (i+1)*0.03, -0.42); penup();
		}
	}

histogram(2, numbins, -1.5, 1.5);
for (int i = 0; i < numbins; i++) {
		//
		if (- 1 + 1.5*(histogram_global[i]/ndisks) > 1.1) {
			printf("histogram bar too high");
		} else {
			move(1.5+ i*0.03, 0.16); pendn(); 
			move(1.5 + i*0.03, 0.16 + 1.5*(histogram_global[i])/ndisks); 
			move(1.5 + (i+1)*0.03, 0.16 + 1.5*(histogram_global[i])/ndisks); 
			move(1.5 + (i+1)*0.03, 0.16); penup();
		}
	}

histogram(3, numbins, -1.5, 1.5);
for (int i = 0; i < numbins; i++) {
		//
		if (- 1 + 1.5*(histogram_global[i]/ndisks) > 1.1) {
			printf("histogram bar too high");
		} else {
			move(1.5+ i*0.03, 0.74); pendn(); 
			move(1.5 + i*0.03, 0.74 + 1.5*(histogram_global[i])/ndisks); 
			move(1.5 + (i+1)*0.03, 0.74 + 1.5*(histogram_global[i])/ndisks); 
			move(1.5 + (i+1)*0.03, 0.74); penup();
		}
	}

drawtexture();			// copies Plane1[] to texture
drawnum(ncount,820,0);	// draws on texture
iPrintf(10,0,"%d disks in the box",ndisks);
iPrintf(10,100,"n - change ndisks");
iPrintf(10,150,"s - change speed");
iPrintf(10,200,"x - exit");

iPrintf(770,445,"X");
iPrintf(770,315,"Y");
iPrintf(770,185,"Vx");
iPrintf(770,55,"Vy");


if (framecount % 60 == 0) {
	energy = getEnergy();
}
iPrintf(10,250,"%f Average Energy", energy);
if (heatbath == 0) {
	iPrintf(1,300,"Thermal Bath is off");
} else {
	iPrintf(1,300,"Thermal Bath is on");
}
drawbuffer();			// copies texture to display
++framecount;
}

void main_tick() {
SDL_Event event;

while (SDL_PollEvent(&event)) {
	switch (event.type) {
		case SDL_KEYDOWN: {
			switch (event.key.keysym.sym) {
				case SDLK_ESCAPE:
				case SDLK_x: {
					quit = 1;
					break; }
				case SDLK_n: {
					ndisks *= 10;
					if(ndisks > 10000) ndisks = 10;
					initialize(ndisks,radius);
					}
				case SDLK_o: {
					float dt;
					static int oldncount = 0, oldticks = 0;
					int newticks = SDL_GetTicks();
					dt = .001 * (newticks - oldticks);
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
				case SDLK_UP: {
					REAL vel = getVavg();
					if (topwall_VY < 0) {
						topwall_VY = 0;
					} else {
						topwall_VY = fabs(vel)/1000;
					}
					break;
				}
				case SDLK_DOWN: {
					REAL vel = getVavg();
					if (topwall_VY > 0) {
						topwall_VY = 0;
					} else {
						topwall_VY = -fabs(vel)/1000;
					}
					break;
				}
				case SDLK_RIGHT: {
					REAL vel = getVavg();
					if (rightwall_VX < 0) {
						rightwall_VX = 0;
					} else {
						rightwall_VX = fabs(vel)/1000;
					}
					break;
				}
				case SDLK_LEFT: {
					REAL vel = getVavg();
					if (rightwall_VX > 0) {
						rightwall_VX = 0;
					} else {
						rightwall_VX = -fabs(vel)/1000;
					}
					break;
				}
				case SDLK_h: {
					if (heatbath) {
						heatbath = 0;
					} else {
						heatbath = 1;
					}
					break;
				}
				default: {
					rightwall_VX = 0;
					topwall_VY = 0;
            	break; 
				}
				}
			}
		}
	}
for(int i = 0;i<slowctl;i++)
	update();
display();
}

void main_loop()
{
#if __EMSCRIPTEN__
	emscripten_set_main_loop(main_tick, -1, 1);
#else
    while (0 == quit){main_tick();}
#endif
}

int main(int argc, char *argv[])
{
if(argc == 2 && argv[1] != NULL) ndisks = atoi(argv[1]);
prefsize(SIZEX,SIZEY);
initgraph("hi there!");

yscale = 1.2;
xscale = yscale * XFAC;
scale(-xscale,xscale,-yscale,yscale);
initialize(ndisks,radius);
labelsize(30);
labelcolor(0,0,255,255);
initlabel();
main_loop();
exitgraph();
return 0;
}
