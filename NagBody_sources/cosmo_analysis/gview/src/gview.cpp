#include "structs.h"
#include "exception.h"
#include "utils.h"
#include "openglprimitives.h"
#include "ver.h"

static int Window = 0;
static settings s(string(getenv("HOME")) + "/.gview.ini");
static bool g_mouse_left_button=false;
static int  g_key_modifiers=0;
bool DEBUG_PRINT=false;

bool RenderGadget(const settings& s);
void SetSnapshot(const int dir,settings& s);
void Timer(int);
int Usage();
string Version();

void Display( void )
{
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	glPushMatrix();
	glTranslatef( s.tx, s.ty, 0  );
	glRotatef(s.rotx, 1, 0, 0);
	glRotatef(s.roty, 0, 1, 0);
	glScalef (s.zoom,s.zoom,s.zoom);

	if (RenderGadget(s) && s.timer<=0) glutSetWindowTitle(("gview :: " + s.snapshot).c_str());

	if (s.showtext){
		GLBegin2D(0,0);
		GLColor(255,55,55);
		assert( s.rs.p0.size()%3==0 );
		string msg="model=" + s.snapshot;
		if (s.rs.snapshots.size()>1) msg +=" [/" + suffix(s.rs.snapshots.size()-1) +"] ";
		msg += " sz=" + tostring(s.rs.p0.size()/3);

		GLString(msg.c_str(),-.99,-.99);

		const io_header_1 x=s.rs.h0;
		msg = "Particles: ";
		for(int i=0;i<6;++i) msg += tostring(x.npart[i]) + " ";
		msg += "\nMasses: ";
		    for(int i=0;i<6;++i) msg += tostring(x.mass[i]) + " ";
		    msg += "\nTotals: ";
		for(int i=0;i<6;++i) msg += tostring(x.npartTotal[i]) + " ";

		msg +=  "\na="         + tostring(s.rs.h0.time) + '\n'
		     +  "z="           + tostring(s.rs.h0.redshift) + '\n'
		     +  "bbox="        + tostring(s.rs.h0.BoxSize) + '\n'
			 +  "Omega0="      + tostring(s.rs.h0.Omega0) + '\n'
		     +  "OmegaLambda=" + tostring(s.rs.h0.OmegaLambda) + '\n'
		     +  "h="           + tostring(s.rs.h0.HubbleParam);

		msg += "\n\nview: ";
		if (s.animate)          msg += "\n  animate";
		if (s.migrateparticles) msg += "\n  migrateparticles";
		if (s.debug)            msg += "\n  debug";
		if (s.gridsize>1)       msg += "\n  gridsize="  + tostring(s.gridsize);
 		if (s.timer>0)          msg += "\n  timer="     + tostring(s.timer);
		if (s.subsample>1)      msg += "\n  subsample=" + tostring(s.subsample);
		if (s.subsnapstep>0)    msg += "\n  subsnapstep="  + tostring(s.subsnapstep);

		const string dbgmsg=s.rs.msg();
		if (dbgmsg.size()>0) msg += "\n\n\n** debug message:\n" + dbgmsg;

		GLString(msg.c_str(),-.99,.97);
		GLEnd2D();
	}

	glPopMatrix();
	glutSwapBuffers();
	if(s.timer>0) glutTimerFunc(static_cast<int>(s.timer),Timer,0);
}

void Exit()
{
	glutDestroyWindow(Window);
	exit(0);
}

void Timer( int )
{
	//sleep(10);
	if (s.animate){
		const float drot=0.3*s.timer/30; //XXX /100
		const float dzoom=0.005;
		static float x=0;

		s.rotx += 2*drot;
		s.roty += drot;

		s.zoom += sin(x)/5000;
		x +=dzoom;
	}
	s.incsubsnap();
	Display();
}

void Redraw( int )
{
   if (s.timer<=0) glutPostRedisplay();
}

void Reshape( int width, int height )
{
	float ratio = (float) width / (float) height;
	s.wx = width;
	s.wy = height;

	glViewport( 0, 0, width, height );
	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();
	glFrustum( -ratio, ratio, -1.0, 1.0, 5.0, 30.0 );
	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity();
	glTranslatef( 0.0, 0.0, -20.0 );
}

void Key( unsigned char key, int x, int y )
{
	const float tx=1.0 * s.wx / 1000;
	const float ty=1.0 * s.wy / 1000;
	assert( s.rs.p0.size()%3==0 );

	switch (key) {
		case 27: Exit(); break;
		case 'q': Exit(); break;
		case 'd': s.tx   += tx; break;
		case 'a': s.tx   -= tx; break;
		case 'w': s.ty   += ty; break;
		case 's': s.ty   -= ty; break;
		case 'x': s.zoom /= 1.05; break;
		case 'z': s.zoom *= 1.05; break;
		case 'r': s.reset(); break;
		case '1': s.axis  = !s.axis; break;
		case '2': s.drawbbox = !s.drawbbox; break;
		case '3': s.drawpoints = !s.drawpoints; break;
		case '4': s.smoothlens = !s.smoothlens; break;
		case '5': s.migrateparticles = !s.migrateparticles; break;
		case '6': s.colorpoints = !s.colorpoints; break;
		case '%': s.animate = !s.animate; break;
		case '!': s.debug= !s.debug; break;
		case 'i': s.showtext = !s.showtext; break;

		case 'm': {s.subsnapstep *= 1.05; if (s.subsnapstep>1.0) s.subsnapstep=1.0; break;}
		case 'M': {s.subsnapstep /= 1.05; break;}
		case 'n': s.setsnapshot( 1); break;
		case 'N': s.setsnapshot(-1); break;
		case 'p': s.setsnapshot(-1); break;
		case 't': {if (s.timer<=0) {s.timer==0 ? s.timer=10 : s.timer=-s.timer;} else s.timer++; if (s.timer>0) glutPostRedisplay(); return;}
		case 'T': {--s.timer; if (s.timer>0) glutPostRedisplay(); return;}
		case '&': {s.timer = -s.timer; if (s.timer>0) glutPostRedisplay(); return;}
		case 'u': ++s.subsample; break;
		case 'U': {--s.subsample; if (s.subsample<=0) s.subsample=1; break;}

		case 'h': Usage(); break;
		case 'g': ++s.gridsize; break;
		case 'G': s.gridsize=max(s.gridsize-1,0); break;
		case 'j': {--s.follow;    if(s.follow<-1) s.follow=s.rs.p0.size()/3-1; break;}
		case 'J': {s.follow-=100; if(s.follow<-1) s.follow=s.rs.p0.size()/3-1; break;}
		case 'k': {++s.follow;    if(static_cast<size_t>(s.follow)*3>=s.rs.p0.size()) s.follow=-1; break;}
		case 'K': {s.follow+=100; if(static_cast<size_t>(s.follow)*3>=s.rs.p0.size()) s.follow=-1; break;}
		case 'l':  s.follow=-1; break;
	}
	Redraw(0);
}

void SpecialKey(int key,int x,int y )
{
	const float dr=10;
	switch (key) {
		case GLUT_KEY_RIGHT: s.roty += dr; break;
		case GLUT_KEY_LEFT:  s.roty -= dr; break;
		case GLUT_KEY_UP:    s.rotx -= dr; break;
		case GLUT_KEY_DOWN:  s.rotx += dr; break;
	}
	g_key_modifiers=glutGetModifiers();
	Redraw(0);
}

void MouseMove(int x, int y )
{
	if (g_mouse_left_button){
		s.rotx = y - s.wy / 2;
		s.roty = x - s.wx / 2;
		Redraw(0);
	}
}

void MouseClick(int button, int state,int x, int y )
{
	if (button==GLUT_LEFT_BUTTON && state==1) g_mouse_left_button=!g_mouse_left_button;
}

string Version() {return string("VERSION: ") + tostring(VERSION) + "." + tostring(VERSION_REV);}

string Config()
{
	std::string s;
	#ifdef NDEBUG
		s+="NDEBUG";
	#else
		#ifdef PROFILE
			s+="PROFILE";
		#else
			s+="DEBUG";
		#endif
	#endif
	if      (sizeof(void*)==4) s+=" 32BIT";
	else if (sizeof(void*)==8) s+=" 64BIT";
	else                       s+=" XXBIT";

	const long one= 1;
	const int big=!(*(reinterpret_cast<const char *>(&one)));
	if (big) s+=" BIGENDIAN";
	else     s+=" LITENDIAN";

	return s;
}

int Usage()
{
	cerr << "gview  - displays a 3D view of a series of gadget snapshot files.\n";
	cerr << "  " << Version() <<  " " << Config() <<  " \n";
	cerr << "  usage: gview [-d | -n] [snapshot file]\n\n";
	cerr << "  option: -d = debug\n";
	cerr << "  option: -n = no reodering\n";
	cerr << "  keyboard: \n";
	cerr << "    a or left arrow: translate left\n";
	cerr << "    d or right arrow: translate right\n";
	cerr << "    s or down arrow: translate down\n";
	cerr << "    w or up arrow: translate up\n";
	cerr << "    z: zoom in\n";
	cerr << "    x: zoom out\n";
	cerr << "    g: show grid and increase grid refinement\n";
	cerr << "    G: decrease grid refinement\n";
	cerr << "    h: show this help string\n";
	cerr << "    i: show snapshot info text on screen\n";
	cerr << "    n: load next snapshot\n";
	cerr << "    p or N: load previous snapshot\n";
	cerr << "    q or ESC: quits program\n";
	cerr << "    r: resets view\n";
	cerr << "    t: increase animation rate\n";
	cerr << "    T: decrease animation rate\n";
	cerr << "    u: increase subsample particle rate (for big files)\n";
	cerr << "    U: decrease subasmple particle rate (for smaller files)\n";
	cerr << "    m: encrease migration stepsize\n";
	cerr << "    M: decrease migration stepsize\n";
	cerr << "    k or K: follow next (k: +1, K:+100), particle\n";
	cerr << "    j or J: follow previous (j: +1, J:+100), particle\n";
	cerr << "    l: follow no particle\n";
	cerr << "    1: enable/disable axis\n";
	cerr << "    2: enable/disable bounding box\n";
	cerr << "    3: enable/disable drawing of points\n";
	cerr << "    4: enable/disable drawing of smoothing lengths\n";
	cerr << "    5: enable/disable animate particle migration (snapshotfiles must include id tags)\n";
	cerr << "    6: enable/disable particle colors\n";
	cerr << "    !: enable/disable debug informations\n";
	cerr << "    %: enable/disable animation\n";
	cerr << "    &: autoplay, rotate and zoom animation (start/stops)\n";
	cerr << "See also gview(1) manual page\n";
	cerr << endl;
	return -1;
}

int main(int argc, char *argv[])
{
	try{
		int a=0;
		s.debug=false;
		s.reorder=true;

		if (argc>1 && (string(argv[1])=="-?" || string(argv[1])=="--?" || string(argv[1])=="-h" || string(argv[1])=="--h")) return Usage();
		else if (argc==3 && string(argv[1])=="-d") {s.debug=true; a=2;}
		else if (argc==3 && string(argv[2])=="-d") {s.debug=true; a=1;}
		else if (argc==2 && string(argv[1])=="-d") {s.debug=true; a=0;}
		else if (argc==2 && string(argv[1])!="-d") {a=1;}
		else if (argc==3 && string(argv[1])=="-n") {s.reorder=false; a=2;}
		else if (argc==3 && string(argv[2])=="-n") {s.reorder=false; a=1;}
		else if (argc==2 && string(argv[1])=="-n") {s.reorder=false; a=0;}
		else if (argc==2 && string(argv[1])!="-n") {a=1;}
		else if (argc!=1) return Usage();

		if (s.debug){
		    DEBUG_PRINT=true;
		    cerr << "## debug output on" << endl;
		}
		s.setsnapshot(a>0 ? argv[a] : "" );

		glutInitWindowPosition(s.wox,s.woy);
		glutInitWindowSize(s.wx,s.wy);
		glutInit( &argc, argv );
		glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH | GLUT_MULTISAMPLE);

		Window = glutCreateWindow(argv[0]);
		if (!Window) {
			cerr << "Error, couldn't open window\n";
			exit(1);
		}

		glutVisibilityFunc( Redraw );
		glutReshapeFunc   ( Reshape );
		glutKeyboardFunc  ( Key );
		glutSpecialFunc   ( SpecialKey );
		glutDisplayFunc   ( Display );
		glutPassiveMotionFunc( MouseMove );
		glutMouseFunc        ( MouseClick );

		glutSetWindowTitle("gview");
		glutShowWindow();

		glutMainLoop();
	}
	catch(const char*   m)   {cout.flush(); cerr << "caught exception chars: "  << m;}
	catch(const string& m)   {cout.flush(); cerr << "caught exception string: " << m;}
	catch(...) 				 {cout.flush(); cerr << "caught unknown exception";}
	cerr << "...aborting\n";
	abort();
}
