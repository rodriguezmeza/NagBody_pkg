#ifndef __OPENGLPRIMITIVES_H__
#define __OPENGLPRIMITIVES_H__

#include <GL/gl.h>			// Header File For The OpenGL32 Library
#include <GL/glu.h>			// Header File For The GLu32 Library
//#include <GL/glaux.h>		// Header File For The Glaux Library
#include <GL/glut.h>

inline void GLSphere(const double x,const double y,const double z,const double radius,const int slices,const int stacks)
{
  static GLUquadricObj* q=gluNewQuadric();

  glPushMatrix();
    glTranslated(x,y,z);
    gluSphere(q,radius,slices,stacks);
  glPopMatrix();
}

inline void GLColorb(const unsigned char R,const unsigned char G,const unsigned char B)
{
  glColor3b(R,G,B);
}

inline void GLColor(const double R,const double G,const double B,const double alpha=1.0)
{
  glColor4d(R,G,B,alpha);
}

inline void GLPointSize(const double s)
{
  glPointSize((float)s);
}

inline void GLLine(const double x1,const double y1,const double z1,const double x2,const double y2,const double z2)
{
  glBegin(GL_LINES);
    glVertex3d( x1, y1, z1);
	  glVertex3d( x2, y2, z2);
  glEnd();
}

inline void GLBox(const double x1,const double y1,const double z1,const double x2,const double y2,const double z2)
{
  glBegin(GL_LINES);
    glVertex3d( x1, y1, z1);
	  glVertex3d( x2, y1, z1);
	  glVertex3d( x2, y1, z1);
	  glVertex3d( x2, y2, z1);
	  glVertex3d( x2, y2, z1);
	  glVertex3d( x1, y2, z1);
	  glVertex3d( x1, y2, z1);
	  glVertex3d( x1, y1, z1);

	  glVertex3d( x1, y1, z2);
	  glVertex3d( x2, y1, z2);
	  glVertex3d( x2, y1, z2);
	  glVertex3d( x2, y2, z2);
	  glVertex3d( x2, y2, z2);
	  glVertex3d( x1, y2, z2);
	  glVertex3d( x1, y2, z2);
	  glVertex3d( x1, y1, z2);

	  glVertex3d( x1, y1, z1);
	  glVertex3d( x1, y1, z2);
	  glVertex3d( x2, y2, z1);
	  glVertex3d( x2, y2, z2);
	  glVertex3d( x1, y2, z1);
	  glVertex3d( x1, y2, z2);
	  glVertex3d( x2, y1, z1);
	  glVertex3d( x2, y1, z2);
	glEnd();
}

inline void GLBoxFill(const double x1,const double y1,const double z1,const double x2,const double y2,const double z2)
{
  glBegin(GL_QUADS);
    glVertex3d( x1, y1, z1);
	  glVertex3d( x2, y1, z1);
	  glVertex3d( x2, y1, z1);
	  glVertex3d( x2, y2, z1);
	  glVertex3d( x2, y2, z1);
	  glVertex3d( x1, y2, z1);
	  glVertex3d( x1, y2, z1);
	  glVertex3d( x1, y1, z1);

	  glVertex3d( x1, y1, z2);
	  glVertex3d( x2, y1, z2);
	  glVertex3d( x2, y1, z2);
	  glVertex3d( x2, y2, z2);
	  glVertex3d( x2, y2, z2);
	  glVertex3d( x1, y2, z2);
	  glVertex3d( x1, y2, z2);
	  glVertex3d( x1, y1, z2);

	  glVertex3d( x1, y1, z1);
	  glVertex3d( x1, y1, z2);
	  glVertex3d( x2, y2, z1);
	  glVertex3d( x2, y2, z2);
	  glVertex3d( x1, y2, z1);
	  glVertex3d( x1, y2, z2);
	  glVertex3d( x2, y1, z1);
	  glVertex3d( x2, y1, z2);
	glEnd();
}

/*
inline void GLBoxFill(const float x1,const float y1,const float z1,const float x2,const float y2,const float z2)
{
   glPushMatrix();
   glTranslatef((x2+x1)/2,(y2+y1)/2,(z2+z1)/2);
   auxSolidBox(x2-x1,y2-y1,z2-z1);
   glPopMatrix();
}
*/

inline void GLPlane(const float x1,const float y1,const float z1,const float x2,const float y2,const float z2)
{
	glBegin(GL_QUADS);
		glVertex3f(x1, y1, z2);
		glVertex3f(x1, y1, z1);

		glVertex3f(x2, y2, z1);
		glVertex3f(x2, y2, z2);
	glEnd();
}

inline void GLText(const double x,const double y,const double z,const char* text,const double scale)
{
	const int n=strlen(text);
	glListBase(1000);
	glPushMatrix();
	glTranslated(x, y, z) ;
	glScaled(scale,scale,scale);
    glCallLists (n, GL_UNSIGNED_BYTE, text);
	glPopMatrix();
	glFlush();
}

inline void GLAxis(const double scale,const double ticks,const bool grey)
{
  const double scalefont=scale/8;
  GLPointSize(1);
  int sub=6;

  if (grey) GLColor(1,1,1);

  if (!grey) GLColor(1,0,0);
  GLText(scale*1,0,0,"X",scalefont);
  GLLine(0,0,0,scale,0,0);
  if(ticks>0) for(double i=0;i<scale;i+=ticks) GLLine(i,-ticks/sub,0,i,ticks/sub,0);

  if (!grey) GLColor(0,1,0);
  GLText(0,scale*1,0,"Y",scalefont);
  GLLine(0,0,0,0,scale,0);
  if (ticks>0) for(double i=0;i<scale;i+=ticks) GLLine(-ticks/sub,i,0,ticks/sub,i,0);
  if (!grey) GLColor(0,0,1);
  GLText(0,0,scale*1,"Z",scalefont);
  GLLine(0,0,0,0,0,scale);
  if (ticks>0) for(double i=0;i<scale;i+=ticks) GLLine(0,-ticks/sub,i,0,ticks/sub,i);
}

inline void GLPoint(const double x,const double y,const double z)
{
  glBegin(GL_POINTS);
    glVertex3d( x, y, z);
  glEnd();
}

inline void GLPoint(const float x,const float y,const float z)
{
  glBegin(GL_POINTS);
    glVertex3f( x, y, z);
  glEnd();
}

#ifdef _WIN32_

inline bool GLSetup(HDC pData)
{
	HDC hDC;
	HFONT hFont;
	GLYPHMETRICSFLOAT agmf[256]; // Throw away
	LOGFONT logfont;

	logfont.lfHeight = -1;
	logfont.lfWidth = 0;
	logfont.lfEscapement = 0;
	logfont.lfOrientation = 0;
	logfont.lfWeight = FW_BOLD;
	logfont.lfItalic = FALSE;
	logfont.lfUnderline = FALSE;
	logfont.lfStrikeOut = FALSE;
	logfont.lfCharSet = ANSI_CHARSET;
	logfont.lfOutPrecision = OUT_DEFAULT_PRECIS;
	logfont.lfClipPrecision = CLIP_DEFAULT_PRECIS;
	logfont.lfQuality = DRAFT_QUALITY;
	logfont.lfPitchAndFamily = DEFAULT_PITCH;
	strcpy(logfont.lfFaceName,"Arial");

	hDC = pData;
	hFont = CreateFontIndirect(&logfont);
	SelectObject (hDC, hFont);

	//create display lists for glyphs 0 through 255 with 0.1 extrusion
	// and default deviation. The display list numbering starts at 1000
	// (it could be any number).
	wglUseFontOutlines(hDC, 0, 255, 1000, 0.0f, 0.01f, WGL_FONT_LINES, agmf);

	DeleteObject(hFont);
  return true;
}
#endif

inline void GLString(const char *text,const float x,const float y,const float newlinestep=0.037)
{
    //glRasterPos2d(x, y);
    //glCallLists(strlen(text), GL_BYTE, (GLbyte *) text);
	float nl=0;
	const int n=strlen(text);
	glRasterPos2f(x,y);
	for(int i=0;i<n;++i) {
		const char c=text[i];
		if (c=='\n') {nl += newlinestep; glRasterPos2f(x,y-nl); }
		else glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,c);
	}
}

inline void GLBegin2D(int width, int height)
{
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(0.0F, (GLfloat) width, 0.0F, (GLfloat) height, -1.0F, 1.0F);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
}

inline void GLEnd2D(void)
{
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
}

#endif // __OPENGLPRIMITIVES_H__
