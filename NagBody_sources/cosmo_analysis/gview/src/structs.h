#ifndef __STRUCTS_H__
#define __STRUCTS_H__

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
#include <cmath>
#include <cassert>
#include <cstring>
#include <cstdlib>

using namespace std;

struct io_header_1
{
	int      npart[6];
	double   mass[6];
	double   time;
	double   redshift;
	int      flag_sfr;
	int      flag_feedback;
	int      npartTotal[6];
	int      flag_cooling;
	int      num_files;
	double   BoxSize;
	double   Omega0;
	double   OmegaLambda;
	double   HubbleParam;
	char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  // fills to 256 Bytes
};

struct render_settings
{
	vector<float>  p0,p1;
	vector<size_t> id0,id1;
	vector<float> smoothlens;
	io_header_1 h0,h1;
	string file0,file1;

	vector<string> snapshots;
	string message;

	const char* msg() {string s=message; message=""; return s.c_str();}
	void operator+=(const string& s){message +=s;}

	string follow_particle_msg;
	int    follow_id;
	float  follow_vx[3];

	render_settings() : follow_id(-1) {}
};

struct settings
{
	bool axis,drawbbox,smoothlens,migrateparticles,drawpoints,colorpoints,animate,showtext,reorder,debug;
	int wx,wy,wox,woy;
	float tx,ty;
	float rotx,roty;
	float zoom;
	float timer;
	string filename,snapshot,snapshot_next;
	int snap,subsample,gridsize,follow;
	float subsnap,subsnapstep;

	mutable render_settings rs;

	settings(const string& filename);
	~settings(){ofstream s(filename.c_str()); s << *this;}

	void setsnapshot(const string& snapshotfile);
	void reset(){string f=filename,s=snapshot; defaultvalues(); filename=f;snapshot=s;}

	void defaultvalues()
	{
		// default values
		axis=drawbbox=drawpoints=animate=true;
		colorpoints=migrateparticles=smoothlens=showtext=reorder=debug=false;
		wx=400;
		wy=300;
		wox=woy=20;
		tx=ty=rotx=roty=0.0;
		zoom=0.09;
		timer=0;
		filename=snapshot=snapshot_next="";
		snap=subsample=1;
		gridsize=0;
		follow=-1;
		subsnap=0;
		subsnapstep=0.05;
	}

	void setsnapshot(const int dir)
	{
		if (rs.snapshots.size()==0) return;
		// 		if (rs.snapshots!=FindSnapshots(snapshot)){
		// 			const int s=snap;
		// 			setsnapshot(snapshot);
		// 			if (rs.snapshots.size()<s) snap=s;
		// 			else snap=0;
		// 		}

		snap += dir;
		snapshot=rs.snapshots[getindex(snap)];
	}

	int getindex(const int i) const
	{
		if (i<0) return getindex(-i);
		assert( i>=0 );
		const size_t sz=rs.snapshots.size();
		const int n=i/sz;
		int m=i%(sz);
		if (n%2!=0) m=sz-1-m;
		assert( m>=0 && m<static_cast<int>(sz) );
		return m;
	}

	void incsubsnap()
	{
		subsnap += subsnapstep;
		if (subsnap>=1.0) {
			subsnap=0;
			setsnapshot(1);
		}
	}

	friend ostream& operator<<(ostream& s,const settings& x)
	{
		s << "% Configfile_begin\n";
		s << "gview_settings:\n";
		s << "  axis " << x.axis << "\n";
		s << "  drawbbox " << x.drawbbox << "\n";
		s << "  smoothlens " << x.smoothlens << "\n";
		s << "  migrateparticles " << x.migrateparticles << "\n";
		s << "  drawpoints " << x.drawpoints << "\n";
		s << "  colorpoints " << x.colorpoints << "\n";
		s << "  animate " << x.animate << "\n";
		s << "  showtext " << x.showtext << "\n";
		s << "  reorder " << x.reorder << "\n";
		s << "  debug " << x.debug << "\n";
		s << "  windowx " << x.wx << "\n";
		s << "  windowy " << x.wy << "\n";
		s << "  windowox " << x.wox << "\n";
		s << "  windowoy " << x.woy << "\n";
		s << "  translatex " << x.tx << "\n";
		s << "  translatey " << x.ty << "\n";
		s << "  rotx " << x.rotx << "\n";
		s << "  roty " << x.roty << "\n";
		s << "  zoom " << x.zoom << "\n";
		s << "  timer " << x.timer << "\n";
		s << "  snapshot " << x.snapshot << "\n";
		s << "  subsample " << x.subsample << "\n";
		s << "  gridsize " << x.gridsize << "\n";
		s << "  follow " << x.follow << "\n";
		s << "  subsnap " << x.subsnap << "\n";
		s << "  subsnapstep " << x.subsnapstep << "\n";
		s << "end" << "\n";
		s << "% Configfile_end\n";
		return s;
	}
};

#endif // __STRUCTS_H__
