// macsdl2.c - a version for SDL2, hopefully will also run in browser windows!
// february 14, 2019

#include <stdio.h>
#include <stdlib.h>
#include "macsdl2.h"

#if __EMSCRIPTEN__
#include <emscripten/emscripten.h>
#endif
#include <SDL2/SDL.h>
#include "SDL_FontCache.h"

struct Rect	gl_screenRect;
int	gl_ix = DEFAULTX,gl_iy = DEFAULTY,gl_BPL,gl_SPL;
float	gl_oleft,gl_oright,gl_obottom,gl_otop;
float	gl_oxfact,gl_oyfact;
int	gl_ixhold,gl_iyhold;
int	gl_ixminclip,gl_ixmaxclip,gl_iyminclip,gl_iymaxclip;
unsigned int Plane1[DEFAULTX*DEFAULTY],Plane2[DEFAULTX*DEFAULTY];
unsigned int gl_color;
char	*gl_DBPtr,*gl_DBmaxPtr,*gl_DBminPtr;
FC_Font* font;
int fontSize = 20;				// default font size
int fontColor[] = {0,0,0,255};
int fontR = 0, fontG = 0, fontB = 0, fontA = 255;
const char fontName[] = "Crimson-Roman.ttf";
static int texturedrawn = 0;

SDL_Window *gl_window = NULL;
SDL_Renderer *gl_renderer = NULL;
SDL_Texture *gl_texture = NULL;

int gl_pen = 0;

void initgraph(const char *gloop)
{
SDL_Init(SDL_INIT_VIDEO);

gl_window = SDL_CreateWindow(
	gloop,
//	SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
0,0,
	gl_ix,gl_iy,
//	SDL_WINDOW_SHOWN | SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE);
	SDL_WINDOW_SHOWN | SDL_WINDOW_OPENGL );
	
//SDL_SetWindowBordered(gl_window,SDL_FALSE);	// works, but no title bar to grab
//SDL_SetWindowFullscreen(gl_window,SDL_WINDOW_FULLSCREEN);

//gl_renderer = SDL_CreateRenderer(gl_window, -1, SDL_RENDERER_ACCELERATED);
gl_renderer = SDL_CreateRenderer(gl_window, -1, SDL_RENDERER_PRESENTVSYNC);

gl_texture = SDL_CreateTexture(
	gl_renderer,SDL_PIXELFORMAT_ARGB8888,SDL_TEXTUREACCESS_STREAMING,gl_ix,gl_iy);

gl_screenRect.left = gl_screenRect.top = 0;
gl_screenRect.right = gl_ix;
gl_screenRect.bottom = gl_iy;

gl_BPL = gl_ix * 4;		// assume 32 bit pixels
gl_SPL = gl_BPL/4;
gl_DBPtr = (char*) Plane1;
gl_DBminPtr = gl_DBPtr;
gl_DBmaxPtr = gl_DBPtr + gl_BPL*gl_iy;
gl_ixminclip = 0;
gl_ixmaxclip = gl_ix-1;
gl_iyminclip = 0;
gl_iymaxclip = gl_iy-1;
}

void prefsize(int ix,int iy)
{
gl_ix=ix;
//gl_ix = 2*(gl_ix/2);
gl_iy=iy;
}

// scale into full screen
void scale(float xlo,float xhi,float ylo,float yhi)
{
float rr,ll,tt,bb,ixhi,ixlo,iyhi,iylo;

gl_oleft=xlo; gl_oright=xhi;
gl_obottom=ylo; gl_otop=yhi;

rr = gl_oright; ll = gl_oleft; tt = gl_otop; bb = gl_obottom;
ixhi = gl_screenRect.right - 1; ixlo = gl_screenRect.left;	// adjust for
iyhi = gl_screenRect.bottom - 1; iylo = gl_screenRect.top;	// mac coordinates
gl_oxfact = (ixhi-ixlo)/(rr-ll);
gl_oyfact = (iyhi-iylo)/(bb-tt);
}

void drawtexture()
{
SDL_UpdateTexture(gl_texture, NULL, Plane1, gl_ix * sizeof (Uint32));
SDL_RenderCopy(gl_renderer, gl_texture, NULL, NULL);
texturedrawn = 1;	
}

void drawbuffer()
{
if(texturedrawn == 0)
	{
	SDL_UpdateTexture(gl_texture, NULL, Plane1, gl_ix * sizeof (Uint32));
	SDL_RenderCopy(gl_renderer, gl_texture, NULL, NULL);
	}
//FC_Draw(font, gl_renderer, 0, 0, "This is %s.\n It works.", "example text");
SDL_RenderPresent(gl_renderer);
texturedrawn = 0;
}

/* 32 */
void color(unsigned int value)
{
unsigned int lval=0;
switch(value)
	{
	case BLACK:
		lval = 0x00000000;
		break;
	case WHITE:
		lval = 0x00ffffff;
		break;
	case RED:
		lval = 0x00ff0000;
		break;
	case GREEN:
		lval = 0x0000ff00;
		break;
	case BLUE:
		lval = 0x000000ff;
		break;
	case YELLOW:
		lval = 0x00ffff00;	// R + G
		break;
	case MAGENTA:
		lval = 0x00ff00ff;	// R + B
		break;
	case CYAN:
		lval = 0x0000ffff;	// B + G
		break;
	default:
		lval = value;
		break;
	}
gl_color = lval;
}

/* 32 */
unsigned int RGBcolor(float r,float g,float b)
{
unsigned int ir,ig,ib;

ir = 256 * r;
if(ir == 256) ir = 255;
ig = 256 * g;
if(ig == 256) ig = 255;
ib = 256 * b;
if(ib == 256) ib = 255;
ir &= 255;
ig &= 255;
ib &= 255;
gl_color = ir<<16 | ig<<8 | ib;
return(gl_color);
}

/* 32 */
void setrgbcolor(float r,float g,float b)
{
gl_color = RGBcolor(r,g,b);
}

/* 32 */
void setgray(float a)
{
gl_color = RGBcolor(a,a,a);
}

/* 32 */
unsigned int setcolor(int ir,int ig,int ib)
{
ir &= 255;
ig &= 255;
ib &= 255;
gl_color = ir<<16 | ig<<8 | ib;
return(gl_color);
}

void clear(void)
{
gl_setbuf(gl_DBPtr,gl_color);
}

/* set buffer to given int value... */
/* assumes gl_ix by gl_iy size */

/*  32  */
void gl_setbuf(char *Buf,unsigned int value)
{
register unsigned int i,icnt,lval;
unsigned int L[2];
register double *Bptr,fval;

lval = value;					// make 2 copies
L[0] = lval; L[1] = lval;
fval = *((double *)L);

//icnt = (gl_ix*gl_iy)>>2;		// 4 shorts to a double
icnt = (gl_BPL*gl_iy)>>3;		// 8 bytes to a double

Bptr = (double*)Buf;
		
for(i=0;i<icnt;i++)
	{
	*Bptr = fval;
	++Bptr;
	}
}

/*  32  */
//#define macro_dotati(x,y) *(gl_DBPtr + y*gl_BPL + x) = gl_color
#if CLIP
#define macro_dotati(x,y) dotati(x,y)
#else
#define macro_dotati(x,y) *((unsigned int*)gl_DBPtr + y*gl_SPL + x) = gl_color
#endif
//#define macro_dotati(x,y) dotati(x,y)

/* bresenham line-drawing algorithm */
void bresh(int x1,int y1,int x2,int y2)
{
int dx,dy,d,ink1,ink2,x,y,oct;

if(x1>x2)
        {
        x=x1;x1=x2;x2=x;
        y=y1;y1=y2;y2=y;
        }
dx = x2 - x1;
dy = y2 - y1;
x = x1; y = y1;
macro_dotati(x,y);

/* octant dispatch table... */

if(y2>y1)
        {
        if(dx>dy)
                oct=1;
        else
                oct=2;
        }
else
        {
        if(dx>-dy)
                oct=3;
        else
                oct=4;
        }

switch(oct){

case(1):
d = (dy<<1) - dx;
ink1 = (dy<<1);
ink2 = ((dy-dx)<<1);
while(x<x2)
         {
         ++x;
         if(d <= 0 )
                 d += ink1;
         else
                 {
                 d += ink2;
                 ++y;
                 }
         macro_dotati(x,y);
         }
break;

case(3):
d = (-dy<<1) - dx;
ink1 = (-dy<<1);
ink2 = ((-dy-dx)<<1);
while(x<x2)
         {
         ++x;
         if(d <= 0 )
                 d += ink1;
         else
                 {
                 d += ink2;
                 --y;
                 }
         macro_dotati(x,y);
         }
break;

case(2):
d = (dx<<1) - dy;
ink1 = (dx<<1);
ink2 = ((dx-dy)<<1);
while(y<y2)
         {
         ++y;
         if(d <= 0 )
                 d += ink1;
         else
                 {
                 d += ink2;
                 ++x;
                 }
         macro_dotati(x,y);
         }
break;

case(4):
d = (dx<<1) + dy;
ink1 = (dx<<1);
ink2 = ((dx+dy)<<1);
while(y>y2)
         {
         --y;
         if(d <= 0 )
                 d += ink1;
         else
                 {
                 d += ink2;
                 ++x;
                 }
         macro_dotati(x,y);
         }
}
}

/*  32  */
void dotati(int x,int y)
{
register unsigned int *lPtr;

if(x<gl_ixminclip)return;
if(x>gl_ixmaxclip)return;
if(y<gl_iyminclip)return;
if(y>gl_iymaxclip)return;
lPtr = (unsigned int *)gl_DBPtr + y*gl_SPL + x;
#if CLIP
if((lPtr>=(unsigned int*)gl_DBminPtr)&&(lPtr<(unsigned int*)gl_DBmaxPtr))
#endif
	*lPtr = gl_color;
}

/* bresenham circle drawing routine */
/*  ..see foley and van dam, pg 87  */

#define dfn dotati
#define cpts(x, y)\
dfn(x+x0,y+y0);dfn(y+x0,x+y0);dfn(-x+x0,y+y0);dfn(y+x0,-x+y0);\
dfn(x+x0,-y+y0);dfn(-y+x0,x+y0);dfn(-x+x0,-y+y0);dfn(-y+x0,-x+y0)

void brcirc(int x0,int y0,int rad)
{
int x,y,d,ink1,ink2;

x=0;
y=rad;
d=1-rad;
ink1=3;
ink2=5-(rad<<1);
cpts(x,y);

while(y>x)
        {
        ++x;
        ink1+=2;
        if(d<0)
        		{
        		d+=ink1;
                ink2+=2;
                }
        else
                {
                d+=ink2;
                ink2+=4;
                --y;
                }
        cpts(x,y);
        }
}

/* bresenham filled circle drawing routine */

void brcircf(int x0,int y0,int rad)
{
int x,y,d,ink1,ink2;

x=0;
y=rad;
d=1-rad;
ink1=3;
ink2=5-(rad<<1);

while(y>x)
        {
        ink1+=2;
        if(d<0)
        		{
        		d+=ink1;
                ink2+=2;

                }
        else
                {
                d+=ink2;
                ink2+=4;
                gl_scan(x0-x,x0+x,y0+y);
                gl_scan(x0-x,x0+x,y0-y);
                --y;
                }
        gl_scan(x0-y,x0+y,y0-x);
        gl_scan(x0-y,x0+y,y0+x);
        ++x;
        }
gl_scan(x0-x,x0+x,y0+y);
gl_scan(x0-x,x0+x,y0-y);
}
/* bresenham ellipse drawing routine */
/*  ..see foley and van dam, pg 90  */

#define epts(x, y)\
dfn(x+x0,y+y0);dfn(-x+x0,y+y0);dfn(x+x0,-y+y0);dfn(-x+x0,-y+y0)

void brellip(int x0,int y0,int a, int b)
{
int x,y;
double d1,d2;

x=0;
y=b;
d1 = b*b - a*a*b + .25*a*a;
epts(x,y);

while(a*a*(y-.5) > b*b*(x+1))
	{
	if(d1 < 0)
		d1 += b*b*(2*x+3);
	else
		{
		d1 += b*b*(2*x+3) + a*a*(-2*y+2);
		y--;
		}
	x++;
	epts(x,y);
	}
	
d2 = b*b*(x+.5)*(x+.5) + a*a*(y-1)*(y-1) - a*a*b*b;
while(y > 0)
	{
	if(d2 < 0)
		{
		d2 += b*b*(2*x+2) + a*a*(-2*y+3);
		x++;
		}
	else
		d2 += a*a*(-2*y+3);
	y--;
	epts(x,y);
	}
}

/* bresenham ellipse drawing routine */

void brellipf(int x0,int y0,int a, int b)
{
int x,y;
double d1,d2;

x=0;
y=b;
d1 = b*b - a*a*b + .25*a*a;
//epts(x,y);

while(a*a*(y-.5) > b*b*(x+1))
	{
	if(d1 < 0)
		d1 += b*b*(2*x+3);
	else
		{
		d1 += b*b*(2*x+3) + a*a*(-2*y+2);
		y--;
		}
	x++;
	gl_scan(x0-x,x0+x,y0-y);
	gl_scan(x0-x,x0+x,y0+y);
	//epts(x,y);
	}
	
d2 = b*b*(x+.5)*(x+.5) + a*a*(y-1)*(y-1) - a*a*b*b;
while(y > 0)
	{
	if(d2 < 0)
		{
		d2 += b*b*(2*x+2) + a*a*(-2*y+3);
		x++;
		}
	else
		d2 += a*a*(-2*y+3);
	y--;
	gl_scan(x0-x,x0+x,y0-y);
	gl_scan(x0-x,x0+x,y0+y);
	//epts(x,y);
	}
}

void origgl_scan(int x1,int x2,int y)		// horizontal line-fill routine...
{
while(x1<=x2)
	{
	dotati(x1,y);
	++x1;
	}
}

/*  32  */
void gl_scan(int x1,int x2,int y)		// horizontal line-fill routine...
											//  fills x1 to x2 inclusive
{
register unsigned int *Ptr1,*Ptr2;

if(y<gl_iyminclip)return;
if(y>gl_iymaxclip)return;
if(x2<gl_ixminclip)return;
if(x1>gl_ixmaxclip)return;
if(x2-x1 < 0) return;
if(x1<gl_ixminclip) x1 = gl_ixminclip;
if(x2>gl_ixmaxclip) x2 = gl_ixmaxclip;

Ptr1 = (unsigned int *)gl_DBPtr + y*gl_SPL + x1;
Ptr2 = (unsigned int *)gl_DBPtr + y*gl_SPL + x2;

while(Ptr1 <= Ptr2)	// whole-word writes
	*Ptr1++ = gl_color;
}

void circ(float x,float y,float rad)		// floating point circle routine, radius in x-scaled units...
//float x,y,rad;
{
int ix,iy,irad;

x -= gl_oleft;
x *= gl_oxfact;
ix=x;
y -= gl_otop;
y *= gl_oyfact;
iy=y;
rad *= gl_oxfact;
irad=rad;
brcirc(ix,iy,irad);
}

void circf(float x,float y,float rad)	// floating point filled circle routine, radius in x-scaled units...
{
int ix,iy,irad;

x -= gl_oleft;
x *= gl_oxfact;
ix=x;
y -= gl_otop;
y *= gl_oyfact;
iy=y;
rad *= gl_oxfact;
irad=rad;
brcircf(ix,iy,irad);
}

void ellipse(float x,float y,float a,float b)		// floating point ellipse routine
{
int ix,iy,ia,ib;

x -= gl_oleft;
x *= gl_oxfact;
ix=x;
y -= gl_otop;
y *= gl_oyfact;
iy=y;
a *= gl_oxfact;
b *= gl_oyfact;
ia=a; ib= -b;
brellip(ix,iy,ia,ib);
}

void ellipsef(float x,float y,float a,float b)	// floating point filled ellipse routine
{
int ix,iy,ia,ib;

x -= gl_oleft;
x *= gl_oxfact;
ix=x;
y -= gl_otop;
y *= gl_oyfact;
iy=y;
a *= gl_oxfact;
b *= gl_oyfact;
ia=a; ib= -b;
brellipf(ix,iy,ia,ib);
}

void brect(int ix1,int iy1,int ix2,int iy2)		// pixel units rectangle drawing routine
{
gl_scan(ix1,ix2,iy1);
gl_scan(ix1,ix2,iy2);
bresh(ix2,iy1,ix2,iy2);
bresh(ix1,iy1,ix1,iy2);
}

void brectf(int x1,int y1,int x2,int y2)
{
while(y1<=y2)
	{
	gl_scan(x1,x2,y1);
	++y1;
	}
}

void rect(float x1,float y1,float x2,float y2)		// floating point rectangle drawing routine
{
int ix1,iy1,ix2,iy2,ihold;

x1 -= gl_oleft;
x1 *= gl_oxfact;
ix1=x1;
y1 -= gl_otop;
y1 *= gl_oyfact;
iy1=y1;
x2 -= gl_oleft;
x2 *= gl_oxfact;
ix2=x2;
y2 -= gl_otop;
y2 *= gl_oyfact;
iy2=y2;
if(ix1>ix2) {ihold=ix1;ix1=ix2;ix2=ihold;}

gl_scan(ix1,ix2,iy1);
gl_scan(ix1,ix2,iy2);
bresh(ix2,iy1,ix2,iy2);
bresh(ix1,iy1,ix1,iy2);
//printf("%d %d %d %d\n",ix1,iy1,ix2,iy2);
}

void recti(int ix1,int iy1,int ix2,int iy2)		// pixel units rectangle drawing routine
{
gl_scan(ix1,ix2,iy1);
gl_scan(ix1,ix2,iy2);
bresh(ix2,iy1,ix2,iy2);
bresh(ix1,iy1,ix1,iy2);
}

void rectf(float x1,float y1,float x2,float y2)		// floating point rectangle fill routine
{
int ix1,iy1,ix2,iy2,ihold;

x1 -= gl_oleft;
x1 *= gl_oxfact;
ix1=x1;
y1 -= gl_otop;
y1 *= gl_oyfact;
iy1=y1;
x2 -= gl_oleft;
x2 *= gl_oxfact;
ix2=x2;
y2 -= gl_otop;
y2 *= gl_oyfact;
iy2=y2;

if(ix1>ix2) {ihold=ix1;ix1=ix2;ix2=ihold;}
if(iy1>iy2) {ihold=iy1;iy1=iy2;iy2=ihold;}
#if TILERECTF
--ix2; --iy2;
#endif
//printf("%d %d %d %d wha?\n",ix1,iy1,ix2,iy2);
while(iy1<=iy2)
	{
	gl_scan(ix1,ix2,iy1);
	++iy1;
	}
}

void gl_rect(int x1,int y1,int x2,int y2)
{
while(y1<=y2)
	{
	gl_scan(x1,x2,y1);
	++y1;
	}
}

/*  *graph-type calls, penup(), pendn(), move(x,y), dotat(x,y)...  */
void penup()
{
gl_pen = 0;
}

void pendn()
{
gl_pen = 1;
}

void move(float x, float y)
{
 register int ix,iy;
 
	x -= gl_oleft;
	x *= gl_oxfact;
	ix = x;
	y -= gl_otop;
	y *= gl_oyfact;
	iy = y;
	
	if(gl_pen)
		bresh(gl_ixhold,gl_iyhold,ix,iy);
	gl_ixhold = ix;
	gl_iyhold = iy;
}

void movei(int ix, int iy)
{
	if(gl_pen)
		bresh(gl_ixhold,gl_iyhold,ix,iy);
	gl_ixhold = ix;
	gl_iyhold = iy;
}

void dotat(float x, float y)
{
 register int ix,iy;
 
	x -= gl_oleft;
	x *= gl_oxfact;
	ix = x;
	y -= gl_otop;
	y *= gl_oyfact;
	iy = y;

	dotati(ix,iy);
}

void linef(float x1, float y1, float x2, float y2)
{
 register int ix1,iy1,ix2,iy2;
 
	x1 -= gl_oleft;
	x1 *= gl_oxfact;
	ix1 = x1;
	y1 -= gl_otop;
	y1 *= gl_oyfact;
	iy1 = y1;
	x2 -= gl_oleft;
	x2 *= gl_oxfact;
	ix2 = x2;
	y2 -= gl_otop;
	y2 *= gl_oyfact;
	iy2 = y2;

	bresh(ix1,iy1,ix2,iy2);
}

void fullscreen(int val)
{
if(val == 0) SDL_SetWindowFullscreen(gl_window,0);
if(val == 1) SDL_SetWindowFullscreen(gl_window,SDL_WINDOW_FULLSCREEN);
if(val == 2) SDL_SetWindowFullscreen(gl_window,SDL_WINDOW_FULLSCREEN_DESKTOP);
}

void exitgraph()
{
SDL_DestroyTexture(gl_texture);
SDL_DestroyRenderer(gl_renderer);
SDL_DestroyWindow(gl_window);
SDL_Quit();
exit(0);
}

static int		initlabelflag = 0;	// so that initlabel() isn't called twice...

// call this one before initlabel(), otherwise default size is 20
void labelsize(int sz)
{
fontSize = sz;
}

// call this one before initlabel(), otherwise default color is black (0,0,0,255)
void labelcolor(int r, int g, int b, int a)
{
fontR = r; fontG = g; fontB = b; fontA = a;
}

void initlabel()
{
font = FC_CreateFont();
FC_LoadFont(font, gl_renderer, fontName, fontSize, FC_MakeColor(fontR,fontG,fontB,fontA), 	TTF_STYLE_NORMAL);
}

void fontet()
{
FC_Draw(font, gl_renderer, 0, 0, "This is %s.\n It works. %d", "example text",123);
}

void drawnum(int n, int xpos, int ypos)
{
FC_Draw(font, gl_renderer, xpos, ypos, "%d", n);
}

void label(char* cname, int xpos, int ypos)
{
FC_Draw(font, gl_renderer, xpos, ypos, "%s", cname);
}

#include <stdarg.h>
void iPrintf(int xpos, int ypos, char *format,...)
{
char Cbuf[256];
va_list args;

va_start(args,format);
vsprintf(Cbuf,format,args);
FC_Draw(font, gl_renderer, xpos, ypos, "%s", Cbuf);
}

