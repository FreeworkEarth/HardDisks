#define DEFAULTX	1024	// default window x size (pixels)
#define DEFAULTY	512		// default window y size (pixels)

struct Rect	{		// original mac definition
	short top;
	short left;
	short bottom;
	short right;
};

#define BLACK                   0
#define WHITE                   1
#define RED						2
#define GREEN					3
#define BLUE					4
#define YELLOW					5
#define MAGENTA					6
#define CYAN					7

void initgraph(const char *gloop);
void prefsize(int ix,int iy);
void scale(float xlo,float xhi,float ylo,float yhi);
void drawtexture();
void drawbuffer();
void clear();
void gl_setbuf(char *Buf,unsigned int value);
void gl_scan(int x1,int x2,int y);
void penup();
void pendn();
void move(float x, float y);
void movei(int ix, int iy);
void dotat(float x, float y);
void linef(float x1, float y1, float x2, float y2);
void fullscreen(int val);
void exitgraph();
void initlabel();
void labelsize(int sz);
void labelcolor(int r, int g, int b, int a);
void fontet();
void drawnum(int n, int xpos, int ypos);
void label(char* cname, int xpos, int ypos);
#include <stdarg.h>
void iPrintf(int xpos, int ypos, char *format,...);

void			color(unsigned int);
unsigned int  	RGBcolor(float,float,float);
void			setrgbcolor(float,float,float);
void			setgray(float);
unsigned int	setcolor(int,int,int);

void 	circ( float, float, float);
void	circf(float x,float y,float rad);
void	ellipse(float x,float y,float a,float b);
void	ellipsef(float x,float y,float a,float b);
void	rect(float x1,float y1,float x2,float y2);
void	recti(int ix1,int iy1,int ix2,int iy2);
void	rectf( float, float, float, float);
void	gl_rect(int x1,int y1,int x2,int y2);

extern char *gl_DBPtr;
