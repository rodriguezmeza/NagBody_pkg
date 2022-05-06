#include <string>
#include "structs.h"
#include "exception.h"
#include "utils.h"
#include "rfile.h"
#include "configfile.h"
#include "triple.h"
#include "array.h"
#include "openglprimitives.h"

const static float kpc2mpc=1E-3; // kpc 2 Mpc

settings::settings(const string& file)
{
	defaultvalues();
	filename=file;

	if(FileExists(file)){
		Configfile c(file);
		if (c.hasEntry("zoom"))       zoom=c.Get<float>("zoom");
		if (c.hasEntry("axis"))       axis=c.Get<bool>("axis");
		if (c.hasEntry("smoothlens")) smoothlens=c.Get<bool>("smoothlens");
		if (c.hasEntry("migrateparticles")) migrateparticles=c.Get<bool>("migrateparticles");
		if (c.hasEntry("drawpoints")) drawpoints=c.Get<bool>("drawpoints");
		if (c.hasEntry("colorpoints"))colorpoints=c.Get<bool>("colorpoints");
		if (c.hasEntry("animate"))    animate=c.Get<bool>("animate");
		if (c.hasEntry("showtext"))   showtext=c.Get<bool>("showtext");
		if (c.hasEntry("reorder"))    debug=c.Get<bool>("reorder");
		if (c.hasEntry("debug"))      debug=c.Get<bool>("debug");
		if (c.hasEntry("drawbbox"))   drawbbox=c.Get<bool>("drawbbox");
		if (c.hasEntry("windowx"))    wx  =c.Get<int>("windowx");
		if (c.hasEntry("windowy"))    wy  =c.Get<int>("windowy");
		if (c.hasEntry("windowox"))   wox =c.Get<int>("windowox");
		if (c.hasEntry("windowoy"))   woy =c.Get<int>("windowoy");
		if (c.hasEntry("translatex")) tx  =c.Get<float>("translatex");
		if (c.hasEntry("translatey")) ty  =c.Get<float>("translatey");
		if (c.hasEntry("rotx"))       rotx=c.Get<float>("rotx");
		if (c.hasEntry("roty"))       roty=c.Get<float>("roty");
		if (c.hasEntry("timer"))      timer=c.Get<float>("timer");
		if (c.hasEntry("snapshot"))   snapshot=c.Get<string>("snapshot");
		if (c.hasEntry("subsample"))  subsample=c.Get<int>("subsample");
		if (c.hasEntry("gridsize"))   gridsize=c.Get<int>("gridsize");
		if (c.hasEntry("follow"))     follow=c.Get<int>("follow");
		if (c.hasEntry("subsnap"))    subsnap=c.Get<float>("subsnap");
		if (c.hasEntry("subsnapstep"))subsnapstep=c.Get<float>("subsnapstep");
	}
	else{
	    cerr << "File <" << file << "> does not exist, please copy it from installdir/gview.ini\n";
	}

	assert( wx>=0 && wy>=0 && wox>=0 && woy>=0);
	assert( timer>=0 );
	assert( subsample>=0 );
	assert( gridsize>=0 );
	assert( subsnap>=0 && subsnap<1.0 );
	assert( subsnapstep>=0 && subsnapstep<1.0 );
}

void settings::setsnapshot(const string& snapshotfile)
{
	if (snapshotfile!="") snapshot=snapshotfile;
	const size_t t=rs.snapshots.size();
	rs.snapshots=FindSnapshots(snapshot);
	if (rs.snapshots.size()>t){
		snap=FindSnapNumber(snapshot,rs.snapshots);
		rs.smoothlens.resize(rs.snapshots.size());
	}
}

void Renorm(vector<float>& y,const float bbox,const float factor=1.0)
{
	assert(bbox>0);
	const float off=bbox/2;

	for(size_t i=0;i<y.size();++i) {
		y[i] -= off;
		y[i] *= factor;
	}
}

string LoadSnapshot(const string& filename,io_header_1& h,vector<float>& p,vector<size_t>& id,const bool reorder,const bool debug)
{
	pair<io_header_1,vector<size_t> > l=Load_sph_data(filename,p,reorder,debug);
	h=l.first;
	id=l.second;
	Renorm(p,h.BoxSize,kpc2mpc);
	return filename;
}

inline float SubPeriodic(const float& p,const float& q,const float& bbox)
{
	const float x0=p-q;
	const bool sign=x0>=float(0);
	const float x1=sign ?  x0 : -x0;
	const float x2=bbox-x1;
	return x1<x2 ? x0 : (sign ? -x2 : x2);
}

// inline void WrapToInsideBox_1D(float& p,const float& bbox)
// {
//
// 	//assert( bbox > float(0) && numeric_limits<T>::is_integer==false && numeric_limits<T>::is_exact==false );
// 	//if      (p<float(0)) {p = bbox-float(fmod(-p/float(1),bbox/float(1)));}
// 	//else if (p>bbox) {p = float(fmod(p/float(1),bbox/float(1)));}
// 	//assert( p>=float(0) && p<bbox );
// 	p += bbox/2;
// 	p = fmod(p,bbox);
// 	if(!( p>=0 && p<bbox )) cerr << "bad point p=" << p;
// 	p -= bbox/2;
// }

inline void Interpol(const vector<float>& p0,const vector<float>& p1,const size_t& i,const float& bbox,const float& subsnap,float *px)
{
	assert( i<p0.size() && i<p1.size() );
	assert( subsnap>=0 && subsnap<1 );
	assert( bbox>0 );

	for(size_t j=0;j<3;++j){
		px[j]=SubPeriodic(p0[i+j],p1[i+j],bbox)*subsnap + p1[i+j];
		//WrapToInsideBox_1D(px[j],bbox); // not needed, produces ok visual output without
	}
}

int PointType(const settings& s,const size_t i,const size_t sz)
{
	assert( i<sz );
	unsigned int j=0,n=0;
	do{
 		assert( j<6 );
 		n += 3*s.rs.h0.npart[j];
 	}
	while(i>n && ++j<6);
	assert( j>=0 && j<6 );
	return j;
}

void ColorPoint(const settings& s,const size_t i,const size_t sz,const float alpha=1.0)
{
	if (!s.colorpoints) return;
	const int j=PointType(s,i,sz);

	if      (j==0) GLColor(255,255,255,alpha);
	else if (j==1) GLColor(155,155,0,alpha);
	else if (j==2) GLColor(0,100,100,alpha);
	else if (j==3) GLColor(0,255,0,alpha);
	else if (j==4) GLColor(0,0,255,alpha);
	else if (j==5) GLColor(255,0,0,alpha);
	else GLColor(1,1,1,alpha);
}

bool RenderGadget(const settings& s)
{
	static bool inRender=false;
	if (inRender) return false;
	inRender=true;

	bool changed=false;

	if (s.snapshot!=s.rs.file0){
		if (s.rs.snapshots.size()>0){
		    if (s.rs.p0.size()>0){
				s.rs.file1=s.rs.file0;
				s.rs.h1   =s.rs.h0;
				swap(s.rs.p1 ,s.rs.p0);
				swap(s.rs.id1,s.rs.id0);
			} else {
				const int m=s.getindex(s.snap-1);
				s.rs.file1=LoadSnapshot(s.rs.snapshots[m],s.rs.h1,s.rs.p1,s.rs.id1,s.reorder,s.debug);
			}
		}
		s.rs.file0=LoadSnapshot(s.snapshot,s.rs.h0,s.rs.p0,s.rs.id0,s.reorder,s.debug);
		s.rs.smoothlens.resize(0);
		changed=true;
 	}

    glEnable(GL_LIGHTING);
    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);

    const float point_alpha=0.5;
    GLPointSize(1.0);
    GLColor(255,255,255,point_alpha);
    glDisable(GL_LIGHTING);

    const vector<float>& p0=s.rs.p0;
    const vector<float>& p1=s.rs.p1;
    const size_t step=3*s.subsample;
    const float bbox=s.rs.h0.BoxSize;
    const float bboxrenorm=bbox*kpc2mpc;
    const float offset    =bboxrenorm/2;

    // draw points
    assert( s.subsample>0 );
    if (s.drawpoints){
		glBegin(GL_POINTS);
		if(s.migrateparticles==false || p1.size()==0){
			for(size_t i=0;i<p0.size();i+=step){
				//GLPoint(p[i],p[i+1],p[i+2]);
				ColorPoint(s,i,p0.size(),point_alpha);
				glVertex3f(p0[i],p0[i+1],p0[i+2]);
			}
		} else {
			assert( p0.size()==p1.size() );
			float px[3];
			for(size_t i=0;i<p0.size();i+=step){
				Interpol(p0,p1,i,bboxrenorm,s.subsnap,px);
				ColorPoint(s,i,p0.size(),point_alpha);
				glVertex3f(px[0],px[1],px[2]);
			}
		}
		glEnd();
	}

	// Draw smoothinglengths
	static bool smoothingexe=true;
	if (s.smoothlens && smoothingexe){
		if (s.rs.smoothlens.size()==0) {
			// load smoothlengths
			const string sfile=s.rs.file0 + ".sml";
			if (smoothingexe && (!FileExists(sfile) || isFileNewer(s.rs.file0,sfile)) ){
				cout << "** Generating smoothinglengths for file <" << s.rs.file0 << ">..." << endl;
				const string cmd="gfilter -s -ns " + s.rs.file0;
				const int n=system( cmd.c_str() );
				if (n!=0) {
					s.rs +=  "\nError running command: '" + cmd + "'";
					smoothingexe=false;
				}
				else cout << "done" << endl;
			}
			if (smoothingexe && FileExists(sfile)) s.rs.smoothlens=loadvector<float>(sfile);
			else {
				s.rs += "\nCould not open or create smoothing file <" + sfile + ">";
				smoothingexe=false;
			}

			vector<float>& l=s.rs.smoothlens;
			const size_t sz=l.size();

			if (sz==0) s.rs += "\nSmoothinglengths of size zero";
 			else {
				float ml=l[0];
				for(size_t i=0;i<sz;++i){
					if (ml<l[i]) ml=l[i];
				}
				l.resize(2*sz);
				for(size_t i=0;i<sz;++i){
					const float beta=asin(max(0.,1.0-l[i]/ml))/90+.0001;
					l[sz+i]=beta;
					l[i] *= kpc2mpc*0.5 + 0.001;
				}
			}
		}

		const vector<float>& l=s.rs.smoothlens;
		if (l.size()>0){
			const size_t offset=l.size()/2;
			assert( offset*2 == l.size() );
			if ( l.size()*3!=2*p0.size() ) throw_("bad smoothing size");

			for(size_t i=0,j=0;i<p0.size();i+=step,++j){
				const float smooth=l[j];
				const float beta=l[j+offset];
				if (s.colorpoints) ColorPoint(s,i,p0.size(),beta);
				else GLColor(1,1,1,beta);
				GLSphere(p0[i],p0[i+1],p0[i+2],smooth,4,3);
			}
		}
	} // end Draw smoothinglengths

	// Draw density
	// 	if (false) {
	// 		static array3d<float> a;
	// 		const string dfile=s.rs.file0 + ".den";
	// 		if (FileExists(dfile) || !isFileNewer(dfile,s.rs.file0)) warn_("** snapshotfile newer that density file");
	// 		if (FileExists(dfile) && a.size()==0) a=array3d<float>::Readascii(dfile);
	// 		if (!a.IsSquare()) throw_("bad density file, not square");
	//
	// 		const float dx=bboxrenorm/a.sizex();
	// 		vector<float> r(3);
	//
	// 		for(array3d<float>::const_iterator itt=a.begin();itt!=a.end();++itt){
	// 			const float d=*itt+1;
	// 			const t_idx idx=itt.toidx();
	// 			r[0]=idx[0]*dx;
	// 			r[1]=idx[1]*dx;
	// 			r[2]=idx[2]*dx;
	// 			Renorm(r,bbox,kpc2mpc);
	// 			GLColor(d,d,d,0.3);
	// 			GLBoxFill(r[0],r[1],r[2],r[0]+dx,r[1]+dx,r[2]+dx);
	// 		}
	// 	}

	// draw axis
	if (s.axis)	{
		glDisable(GL_LIGHTING);
		GLAxis(25,5,false);
	}

	// draw bounding box
    if (s.drawbbox) {
		GLColor(0.7,0.7,0.7);
		GLBox(-offset,-offset,-offset,+offset,+offset,+offset);
	}

	// draw grids
	if (s.gridsize>1){
		const float x0=-offset;
		const float dx=(-2*x0)/s.gridsize;
		GLColor(0.6,0.6,0.6,0.25);

		for(int k=0;k<s.gridsize;++k)
		for(int j=0;j<s.gridsize;++j)
		for(int i=0;i<s.gridsize;++i){
			const float x=x0+i*dx;
			const float y=x0+j*dx;
			const float z=x0+k*dx;
			GLBox(x,y,z,x+dx,y+dx,z+dx);
		}
	}

	// draw follow particle sphere and text
	if (s.follow>=0 && static_cast<size_t>(s.follow)*3<p0.size() ){
		const size_t i=3*s.follow;
		float px[3];

		if(s.migrateparticles==false || p1.size()==0){
			px[0]=p0[i];
			px[1]=p0[i+1];
			px[2]=p0[i+2];
		} else {
			assert( p0.size()==p1.size() );
			Interpol(p0,p1,i,bboxrenorm,s.subsnap,px);
		}

 		GLColor(255,255,255,0.4);
		GLPointSize(4.0);
		GLPoint(px[0],px[1],px[2]);

		if(s.colorpoints) GLColor(255,255,255,0.15);
		else              GLColor(0,255,255,0.15);

		const float r=bbox/45000;
		GLSphere(px[0],px[1],px[2],r,20,20);

		string msg="Follow particle: n=" + tostring(s.follow) ;
		const int t=PointType(s,i,p0.size());
		msg += "    type=";
		switch(t){
			case 0 : msg += "gas\n";  break;
			case 1 : msg += "halo\n"; break;
			case 2 : msg += "disk\n"; break;
			case 3 : msg += "bulge\n"; break;
			case 4 : msg += "star\n";  break;
			case 5 : msg += "bndry\n"; break;
			default : msg += "??\n"; break;
		}
		if (t>=0 && t<6){
			const float m=s.rs.h0.mass[t];
			msg += "    mass=" + tostring(m) + " msun10E10\n";
		}
		const float f=offset;
		msg += "    p=( " + tostring(p0[i]+f) + " ; " + tostring(p0[i+1]+f) + " ; " + tostring(p0[i+2]+f) + " ) h^{-1} Mpc \n";

		// add text direct from snapshot file
		static bool hasGhead=true;
		if(hasGhead && (s.rs.follow_id!=s.follow || changed))
		try{
			s.rs.follow_id=s.follow;
			assert( s.rs.id0.size()==0 || s.rs.follow_id<static_cast<int>(s.rs.id0.size()) );
			const size_t id=s.rs.id0.size()>0 ? s.rs.id0[s.rs.follow_id] : s.rs.follow_id;
			string msg2=System("ghead " + s.rs.file0 + " -p " + tostring(id),true,true);

			size_t n=msg2.find("pos=(");
			if (n!=string::npos)  msg2=msg2.substr(0,n) + "\n    "  + msg2.substr(n,-1);
			n=msg2.find("vel=(");
			if (n!=string::npos)  msg2=msg2.substr(0,n) + "\n    " + msg2.substr(n,-1);
			s.rs.follow_particle_msg=msg2;

			istringstream is(msg2.substr(n+10,-1));
			for(int i=0;i<3;++i){
				s.rs.follow_vx[i]=0;
				char c;
				is >> s.rs.follow_vx[i] >> c;
				s.rs.follow_vx[i] /= (bbox/1000);
			}
		}
		catch(...){hasGhead=false; s.rs.follow_particle_msg="";}
		if (hasGhead) {
			msg += "\n" + s.rs.follow_particle_msg;
			GLColor(255,255,255,0.4);
			GLPointSize(4.0);
			GLLine(px[0],px[1],px[2],px[0]+s.rs.follow_vx[0],px[1]+s.rs.follow_vx[1],px[2]+s.rs.follow_vx[2]);
		}

		GLBegin2D(0,0);
			GLColor(255,55,55);
			GLString(msg.c_str(),0.3,.97);
		GLEnd2D();
	}

	inRender=false;
	return changed;
}
