// keyboard (and mouse) stuff
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <GLUT/glut.h>
#include <OpenGL/OpenGL.h>	// for CG stuff
#include "glgraph.h"
#include "rstuff.h"

//#define REAL float
//#define XFAC (REAL)SIZEX/SIZEY
// big box outline...
//#define XW1 -1.0*XFAC

extern int firstex;			// signal to videoUpdateFrame()

void keyCB(unsigned char key, int x, int y)
{
if(glutGetModifiers() & 1) gl_shift = 1;
if(glutGetModifiers() & 2) gl_ctl = 1;
if(glutGetModifiers() & 4) gl_opt = 1;
//printf("gl_ctl %d\n",gl_ctl);

gl_keys[key] = 1;	// mark key pressed
//printf("key down %d\n",key);

switch(key)
	{
	long i;
	case 27:		// escape
		printf("\nrob shaw  may 2011\n");
		printf("rob@haptek.com\n\n");
		exit(0);
	case 'b':	// for adjusting radii
		radctl ^= 1;
		printf("radctl %d\n",radctl);
		break;
	case 'B':
		dispctl ^= 1;
		printf("dispctl %d\n",dispctl);
		break;
	case 'c':		// reverse all velocities!
		for(i=0;i<nspheres;i++)
			{
			VX[i] = -VX[i];
			VY[i] = -VY[i];
			}
		vpist1 = -vpist1;
		printf("reverse velocities!\n");
		break;
	case 'C':		// set square boundaries
		break;
	case 'd':		// change display mode
		dispctl += 1;
		if(dispctl > 3) dispctl = 0;
		printf("dispctl %d\n",dispctl);
		break;
	case 'D':		// change second display mode
		Dispctl += 1;
		if(Dispctl > 3) Dispctl = 0;
		printf("Dispctl %d\n",Dispctl);
		break;
	case 'e':		// reset disk energies to one
		reseteng();
		break;
	case 'E':		// double all velocities
		{
		int i;
		for(i=0;i<nspheres;i++)
			{VX[i] *= 2; VY[i] *= 2;}
		break;
		}
	case 'f':
		{
		gl_fullscreen ^= 1;
		if(gl_fullscreen)
			glutFullScreen();
		else
			{
//			printf("%d %d\n",gl_xsize,gl_ysize);
			glutPositionWindow(gl_left, gl_top);
			glutReshapeWindow(gl_ix, gl_iy);
			}
		break;
		}
	case 'r':
		{
		gl_flip ^= 1;
		break;
		}
	case 'o':
		{
		float dt;
		static long oldncount = 0;
		extern long ncount;
//		printf("gl_tm.tv_sec %d\n",gl_tm.tv_sec);
		gettimeofday(&gl_tm1,NULL);
		dt = gl_tm1.tv_sec - gl_tm.tv_sec + .000001 * (gl_tm1.tv_usec - gl_tm.tv_usec);
		gl_tm.tv_sec = gl_tm1.tv_sec;
		gl_tm.tv_usec = gl_tm1.tv_usec;
		printf("elapsed time %.2f\t",dt);
		printf("frames per second %.2f\n",gl_framecount/dt);
		gl_framecount = 0;
		printf("iterations per second %.2f\n",(ncount - oldncount)/dt);
		oldncount = ncount;
		break;
		}
#if 0
	case 'v':		// toggle video frame rate counting
		{
		gl_vframerate ^= 1;
		if(gl_vframerate)
			{
			printf("video frame rate\n");
			gl_delt *= 8;	// speed up buttons
			}
		else
			{
			printf("timer frame rate\n");
			gl_delt /= 8;
			}
		break;
		}
#endif
	case 'h':		// turn on histogramming
		++histctl;
		if(histctl > 3) histctl = 0;
		printf	("histctl %ld\n",histctl);
		break;
	case 'H':		// clear histogram
		{
		long i;
		for(i=0;i<NENTS;i++) Histo1[i] = Histowall[i] = HistoW[i] = 0;
		break;
		}
	case 'i':		// toggle black/white
		{
		gl_bwctl ^= 1;
		if(gl_bwctl)
			{
			g_backcolor(BLACK);
			g_color(WHITE);
			}
		else
			{
			g_backcolor(WHITE);
			g_color(BLACK);
			}
		break;
		}
	case 'p':		// print parameters
		{
		gl_print ^= 1;
		if(gl_print == 1)
			{
//			printf("amp %.2f  nlines %ld\n",gl_amp,gl_lines);
			printf("xpos %.2f  ypos %.2f  zpos %.2f\n",gl_xpos,gl_ypos,gl_zpos);
			printf("xrot %.2f  yrot %.2f  zrot %.2f\n",gl_xrot,gl_yrot,gl_zrot);
			printf("printing on\n");
			}
		else
			printf("printing off\n");
		break;
		}
//	case 'P':		// switch pistons...
	case ' ':		// start off pistons...
//		pistctl ^= 1;
		pistctl = 1;
		ncount = 1;
		if(pistctl) printf("right piston\n");
		else printf("left piston\n");
		break;
	case '/':		// print help
		{
		printhelp();
		break;
		}
#if 0
	case 'n':		// reset piston
		if(gl_opt)
			{
			xpist1 = 0;
			xpist2 = XW2;
			}
		else
			{
			xpist1 = XW1;
			xpist2 = XW2;
			}
		break;
#endif
//extern long mdex;
	case 'N':		// reset
//		initspheres9(nrad1,nrad2,seed);		// random throughout right side
		ncount = 0;
//		mdex = 0;
//		copytank();
		break;
	case 'm':		// change mode
//		++mode;
//		if(mode > 5) mode = 1;
		break;
	case 'M':		// toggle absorber
//		Mode ^= 1;
		break;
	case 'R':
		gl_xpos = gl_ypos = gl_zpos = 0;
		gl_xrot = gl_yrot = gl_zrot = 0;
		break;
	case 'a':
		Slowctl *= 10;
		if(Slowctl > 1000) Slowctl = 1;
		printf("Slowctl %ld\n",Slowctl);
		break;
	case 'A':
		Slowctl /= 10;
		if(Slowctl == 0) Slowctl = 1;
		printf("Slowctl %ld\n",Slowctl);
		break;
	case 's':
		if(gl_opt) keyslowctl /= 2;
		else keyslowctl *= 10;
		if(keyslowctl > 10000) keyslowctl = 1;
		if(keyslowctl == 0) keyslowctl = 1;
		printf("slowctl %ld\n",keyslowctl);
		break;
	case 19:		// ctl s
		keyslowctl *= 2;
		printf("slowctl %ld\n",keyslowctl);
		break;
	case 'S':
		keyslowctl /= 10;
		if(keyslowctl == 0) keyslowctl = 1;
		printf("slowctl %ld\n",keyslowctl);
		break;
	case 't':		// trail control
		++trailctl;
		color(WHITE);
		clear();
		if(trailctl > 3) trailctl = 0;
		printf("trailctl %ld\n",trailctl);
		break;
	case 'T':		// erase trails
		color(WHITE);
		clear();
		break;
//	case ' ':		// (erase trails)  perturb single disk
	case 'n':		// reset pistons
//		target(2);
//		color(WHITE);
//		clear();
//		target(1);
//		X1[0] += .0000001;
//		X1[0] += 1e-19;			// minimum perturbation for long double		(64 bit mantissa)
//		X1[0] += 1e-16;			// minimum perturbation for double			(52 bit mantissa)
//		X1[0] += 1e-7;			// minimum perturbation for float			(23 bit mantissa)
		if(gl_shift)
			{
			printf("tweak one!\n");
//			X1[0] += 1e-18;
			}
//		else
//			tweakem(tweakval,2);
		ncount = 0;
//		mdex = 0;
		break;
#if 0
	case 'T':		// tank views
		++tankview;
		if(tankview > 3) tankview = 1;
		if(tankview == 1) printf("first tank\n");
		if(tankview == 2) printf("second tank\n");
		if(tankview == 3) printf("both tanks\n");
		break;
#endif
	case 'V':
		gl_sync ^= 1;
		CGLSetParameter(CGLGetCurrentContext(), kCGLCPSwapInterval, (GLint*)&gl_sync);
		if(gl_sync) printf("screen sync on\n");
		else printf("screen sync off\n");
		break;
	REAL oldrad;
	case '=':			// increase radius 1
		oldrad = rad1;
		if(gl_opt) rad1 += .01;
		else rad1 += .001;
		changerad(oldrad,rad1);
#if FLOATRES == 80
		printf("rad1 %Lf\n",rad1);
#else
		printf("rad1 %f\n",rad1);
#endif
		break;
	case '-':			// decrease radius 1
		oldrad = rad1;
		if(gl_opt) rad1 -= .01;
		else rad1 -= .001;
		changerad(oldrad,rad1);
#if FLOATRES == 80
		printf("rad1 %Lf\n",rad1);
#else
		printf("rad1 %f\n",rad1);
#endif
		break;
	case '+':			// increase radius 2
		oldrad = rad2;
		if(gl_opt) rad2 += .01;
		else rad2 += .001;
		changerad(oldrad,rad2);
#if FLOATRES == 80
		printf("rad2 %Lf\n",rad2);
#else
		printf("rad2 %f\n",rad1);
#endif
		break;
	case '_':			// decrease radius 2
		oldrad = rad2;
		if(gl_opt) rad2 -= .01;
		else rad2 -= .001;
		changerad(oldrad,rad2);
#if FLOATRES == 80
		printf("rad2 %Lf\n",rad2);
#else
		printf("rad2 %f\n",rad1);
#endif
		break;
	case '!':			// shift 1 (to back through presets)
		break;
#if 0
	case '1':
		firstex = 0;	// for image reset, glTexImage2D()
//		++texctl;
//		if(texctl > 1)
//			texctl = 0;
		break;
	case '2':
		++gl_overlay;
		if(gl_overlay > 2)
			gl_overlay = 0;
		break;
	case '3':
		--gl_timer;
		if(gl_timer < 0) gl_timer = 0;
		printf("timer %d  max fps %.1f\n",gl_timer,1000.0/gl_timer);
		break;
	case '#':
		gl_timer -= 10;
		if(gl_timer < 0) gl_timer = 0;
		printf("timer %d  max fps %.1f\n",gl_timer,1000.0/gl_timer);
		break;
	case '4':
		++gl_timer;
		printf("timer %d  max fps %.1f\n",gl_timer,1000.0/gl_timer);
		break;
	case '$':
		gl_timer += 10;
		printf("timer %d  max fps %.1f\n",gl_timer,1000.0/gl_timer);
		break;
#endif
	case '1':
		modectl = 1;
		ncount = 0;
		break;
	case '2':
		modectl = 2;
		ncount = 0;
		break;
	case '3':
		modectl = 3;
		ncount = 0;
		break;
	case '4':
		modectl = 4;
		ncount = 0;
		break;
	case '5':
		modectl = 5;
		ncount = 0;
		break;
	case '6':
		break;
	case '7':
		break;
	case '8':
		break;
	case '9':
//		meanfac *= 2;
//		printf("meanfac %Lf\n",meanfac);
		break;
	case '0':
//		meanfac /= 2;
//		printf("meanfac %Lf\n",meanfac);
		break;
	}
}

void keyUpCB(unsigned char key, int x, int y)
{
	gl_keys[key] = 0;			// mark key released
	if(key >= 97)
		gl_keys[key-32] = 0;	// also mark shifted key released
	if(key >= 65 && key <= 93)
		gl_keys[key+32] = 0;
	if(key == 60)				// special case '<'
		gl_keys[44] = 0;
	if(key == 62)				// special case '>'
		gl_keys[46] = 0;
	if(key == 44)				// special case ','
		gl_keys[60] = 0;
	if(key == 46)				// special case '.'
		gl_keys[62] = 0;
	gl_shift = gl_ctl = gl_opt = 0;	// turn these on and off every keystroke
//	printf("key %d\n",key);
//	gl_count = 100 / (2*gl_lines) + 20;	// set delay until key repeat...
}

void specialKeyPressed(int key, int x, int y)
{
REAL rad;

gl_keys[key+256] = 1;	// mark key pressed
if(glutGetModifiers() & 1) gl_shift = 1;
if(glutGetModifiers() & 2) gl_ctl = 1;
if(glutGetModifiers() & 4) gl_opt = 1;

if(gl_keys[101+256])	// up arrow
	{
	if(gl_shift) rad = rad2;
	else rad = rad1;
	if(gl_opt)
		addisks1(10,rad);
	else if(gl_ctl)
		addisks1(100,rad);
	else
		addisks1(1,rad);
	}

if(gl_keys[103+256])	// down arrow
	{
	if(gl_shift) rad = rad2;
	else rad = rad1;
	if(gl_opt)
		subdisks1(10,rad);
	else if(gl_ctl)
		subdisks1(100,rad);
	else
		subdisks1(1,rad);
//	resetflg = 1;		// to call reseteng()
	}

if(gl_keys[100+256])	// left arrow
	{
	if(!gl_opt && !radctl)			// not a display rotation or radius change
		{
		if(gl_ctl)		// move right-hand piston
			{
			if(gl_shift)
				vpist2 = -5 * dxpist / (slowctl * DELT);
			else
				vpist2 = -dxpist / (slowctl * DELT);
			}
		else			// move left-hand piston
			{
			if(gl_shift)
				vpist1 = -5 * dxpist / (slowctl * DELT);
			else
				vpist1 = -dxpist / (slowctl * DELT);
			}
		}
	if(!gl_opt && radctl)			// a radius change, not a display rotation
		{
		int i;
		if(gl_shift)
			{
			rad1 -= .01;
			for(i=0;i<nspheres;i++) Radius[i] -= .01;
			}
		else
			{
			rad1 -= .001;
			for(i=0;i<nspheres;i++) Radius[i] -= .001;
			}
		printf("rad1 %f\n",rad1);
		}
#if FLOATRES == 80
//	printf("vpist1 %Lf\n",vpist1);
//	printf("vpist2 %Lf\n",vpist2);
#else
//	printf("vpist1 %f\n",vpist1);
//	printf("vpist2 %f\n",vpist2);
#endif
	}

if(gl_keys[102+256])	// right arrow
	{
	if(!gl_opt && !radctl)			// not a display rotation or radius change
		{
		if(gl_ctl)		// move right-hand piston
			{
			if(gl_shift)
				vpist2 = 5 * dxpist / (slowctl * DELT);
			else
				vpist2 = dxpist / (slowctl * DELT);
			}
		else			// move left-hand piston
			{
			if(gl_shift)
				vpist1 = 5 * dxpist / (slowctl * DELT);
			else
				vpist1 = dxpist / (slowctl * DELT);
			}
		}
	if(!gl_opt && radctl)			// a radius change, not a display rotation
		{
		int i;
		if(gl_shift)
			{
			rad1 += .01;
			for(i=0;i<nspheres;i++) Radius[i] += .01;
			}
		else
			{
			rad1 += .001;
			for(i=0;i<nspheres;i++) Radius[i] += .001;
			}
		printf("rad1 %f\n",rad1);
		}
#if FLOATRES == 80
//	printf("vpist1 %Lf\n",vpist1);
//	printf("vpist2 %Lf\n",vpist2);
#else
//	printf("vpist1 %f\n",vpist1);
//	printf("vpist2 %f\n",vpist2);
#endif
	}

//	printf("special key %d\n",key);
}

void specialKeyUp(int key, int x, int y)
{
	gl_keys[key+256] = 0;	// mark key released
gl_shift = gl_ctl = gl_opt = 0;	// turn these on and off every keystroke
if(key == 100)				// left arrow
	{vpist1 = 0; vpist2 = 0;}
if(key == 102)				// right arrow
	{vpist1 = 0; vpist2 = 0;}
}

// button = 0 default
// button = 1 for alt click
// button = 2 for ctl click
void clickCB(int button, int state, int x, int y)	// mouse click
{
printf("%d %d %d %d\n",button,state,x,y);
if(!state)
	glutSetCursor(GLUT_CURSOR_INFO);	// little hand icon
else
	glutSetCursor(GLUT_CURSOR_NONE);
	
gl_imx = x; gl_imy = y;
gl_butnum = button; gl_butstate = state;
gl_mousex = (float)x/(float)gl_xsize;
gl_mousey = (float)y/(float)gl_ysize;
printf("%f %f\n",gl_mousex,gl_mousey);
}

void mouseCB(int x, int y)	// mouse motion
{
printf("%d %d\n",x,y);
}

//void joyCB(unsigned int mask, int x, int y, int z)
//{
//printf("%x, %d %d %d\n",mask,x,y,z);
//}

void updatekeys()
{
if(gl_opt) {
if(gl_keys[100+256])	// left arrow
	gl_yrot -= 1 * gl_delt;
if(gl_keys[102+256])	// right arrow
	gl_yrot += 1 * gl_delt;
if(gl_keys[101+256])	// up arrow
	gl_zpos += 10 * gl_delt;
if(gl_keys[103+256])	// down arrow
	gl_zpos -= 10 * gl_delt;
if(gl_keys['['])				// left
	gl_xpos -= 10 * gl_delt;
if(gl_keys[']'])				// right
	gl_xpos += 10 * gl_delt;
if(gl_keys['{'])				// down
	gl_ypos += 10 * gl_delt;
if(gl_keys['}'])				// up
	gl_ypos -= 10 * gl_delt;
if(gl_keys['.'])				// x axis rotations
	gl_xrot -= 1 * gl_delt;
if(gl_keys[','])
	gl_xrot += 1 * gl_delt;
if(gl_keys['<'])				// z axis rotations
	gl_zrot -= 1 * gl_delt;
if(gl_keys['>'])
	gl_zrot += 1 * gl_delt; }
else {
	if(gl_keys[100+256])	// left arrow
		{
		if(gl_ctl)			// right piston
			{
			if(xpist2 < XW1)
				{
				xpist2 = XW1;
				vpist2 = 0;
				}
			}
		else				// left piston
			{
			if(xpist1 < XW1)
				{
				xpist1 = XW1;
				vpist1 = 0;
				}
			}
		}
	if(gl_keys[102+256])	// right arrow
		{
		if(gl_ctl)			// right piston
			{
			if(xpist2 > XW2)
				{
				xpist2 = XW2;
				vpist2 = 0;
				}
			}
		else				// left piston
			{
			if(xpist1 > XW2)
				{
				xpist1 = XW2;
				vpist1 = 0;
				}
			}
		}
	}
}

void printhelp()
{
//printf("\n1 - toggle Vptr[], Vout[] display\n");
//printf("2 - video, overlay, video and overlay\n");
printf("\nleft, right arrow keys - move left piston (shift increases speed)\n");
printf("control left, right arrow keys - move right piston (shift increases speed)\n");
printf("n - reset piston to left side\n");
printf("option n - reset piston to center\n");
printf("N - reset disks to right side\n");
printf("up, down arrow keys - add, subtract disks\n");
printf("\tshift for smaller disks, ctl for 100 disks at a time\n");
printf("h - change histogrammer mode\n");
printf("H - erase history histogram\n");
printf("=,- - increase, decrease disk radius 1 (opt for greater change)\n");
printf("+,_ - increase, decrease disk radius 2 (opt for greater change)\n");
printf("t - turn on trails\n");
printf("space - erase trails\n");
printf("d - toggle display\n");
printf("s, S - change time steps per display draw\n");
printf("a, A - change screen displays per display draw\n");
printf("c - reverse velocities\n");
printf("\n3, 4 - decrease, increase timer interval\n");
printf("shift 3, 4 - steps of 10 milliseconds\n");
printf("option up, down arrow keys - zoom\n");
printf("option left, right arrow keys - y axis rotation\n");
printf("option [ , ] keys - move left, right\n");
printf("option { , } keys - move down, up\n");
printf("option , , . keys - x axis rotation\n");
printf("option < , > keys - z axis rotation\n");
//printf("9 , 0 keys - vary number of lines\n");
////printf("- , = keys - vary line displacement amplitude\n");
printf("f - toggle full screen\n");
printf("i - toggle black, white inversion\n");
printf("r - toggle mirror reflection\n");
printf("p - print parameters\n");
printf("R - reset position\n");
printf("o - print frames per second\n");
//printf("v - sync to video frame rate\n");
//printf("s - sync to screen refresh\n");
printf("/ - print help\n");
printf("<esc> - exit\n\n");
}
