template<class T>
inline T Readtyp(ifstream& s){
	T x;
	s.read(reinterpret_cast<char*>(&x),sizeof(x));
	if(!s) throw_("bad stream");
	return x;
}

template<class T>
inline void Writetyp(ofstream& o,const T& x){
	o.write(reinterpret_cast<const char*>(&x),sizeof(x));
	if(!o) throw_( "bad stream");
}

inline string Readstring(ifstream& s){
	char c=0;
	string t;
	do{
		c=Readtyp<char>(s);
		if( c!=0) t+=c;
	}
	while (c!=0);
	return t;
}

void readtag(ifstream& s,const int tag)
{
	int t=0;
	assert(sizeof(t)==4);
	if (!s) throw_("bad stream");
	s.read(reinterpret_cast<char*>(&t),sizeof(t));
	if (!s) throw_("bad stream");
	// cerr << "tag=" << tag << " t=" << t;
	if (tag!=t) throw_(string("bad tag in stream, found <") + tostring(t) + "> expected <" + tostring(tag) + ">");
}

template<class T>
void readbinary(istream& s,const size_t n,T& data)
{
	if(!s) throw_("bad stream");
	const size_t readsize=n*sizeof(T);
	s.read(reinterpret_cast<char*>(&data),readsize);
	if(!s) throw_("bad stream");
}

template<class T>
void readvector(istream& s,vector<T>& v){readbinary(s,v.size(),v.front());}

inline string Readline(istream& is)
{
	char buff[16*1024];
	is.getline(buff,16*1024);
	return string(buff);
}

template<class T>
vector<T> loadvector(const string& file)
{
	ifstream s(file.c_str());
	Readline(s);
	vector<T> v;
	while(true){
		T x;
		s >> x;
		if (!s) return v;
		v.push_back(x);
	}
}

pair<io_header_1,vector<size_t> > Load_sph_data(const string& file,vector<float>& v,const bool reorder,const bool debug)
{
	if (debug) cerr << "## Load_sph_data <" << file << ">" << endl;

	io_header_1 header1;
	vector<size_t> r;
	memset(&header1,sizeof(header1),0);

	assert( sizeof(float)==4 );
	assert( sizeof(int)==4 );
	assert( sizeof(io_header_1)==256 );

	ifstream s(file.c_str(),ios::binary);
	if(!s) throw string("file <") + file + ">" + " does not exsist";

	if (debug) cerr << "## reading header...[begintags]";
	readtag(s,256);
	if (debug) cerr << "[data]";
	readbinary(s,1,header1);
	if (debug) cerr << "[endtags]";
	readtag(s,256);
	if (debug) cerr << "...done" << endl;

	size_t NumPart=0;
	for(int k=0; k<6; k++) NumPart += header1.npart[k];

	if (debug) cerr << "## reading position data...(" << NumPart << ")";
	if (debug) cerr << "[begintags]";
	const int disksize=NumPart*3*sizeof(float);
	readtag(s,disksize);
	v.resize(NumPart*3);
	if (debug) cerr << "[data]";
	readvector(s,v);

// temp test cluster data
//  const float offset=495000;
//  for(size_t j=0;j<v.size();++j){
// 	if (j<10) cout << "p[" << j << "]=" << v[j] << endl;
// 	if      (v[j]>offset) v[j]-= offset;
// 	//else if (v[j]<0)      v[j]+= offset/2;
//  }

	if (debug) cerr << "[endtags]";
	readtag(s,disksize);
	if (debug) cerr << "...done" << endl;

	if (reorder){
		bool sortneeded=false;
		bool errorloadingids=false;

		try{
			if (debug) cerr << "## reorder...";

			if (debug) cerr << "[skip vel data]";
			readtag(s,disksize); // velo beg
			size_t pos=s.tellg();
			s.seekg(pos+disksize,ios::beg);
			readtag(s,disksize); // velo end

			if (debug) cerr << "[read ids]";
			readtag(s,NumPart*sizeof(int)); // id beg
			r.resize(NumPart);
			#ifndef NDEBUG
				for(size_t i=0;i<NumPart;++i) r[i]=-1;
			#endif
			size_t idmax=0,idmin=NumPart;
			for(size_t i=0;i<NumPart;++i) {
				size_t id=0;
				readbinary(s,1,id);
				if ( id<1 ) throw_("ID is less than 1, is ID convention [1;N] or [0;N-1] ?");
				idmax=max(idmax,id);
				idmin=min(idmin,id);
				--id; // gadget ids in range [1;N] ?
				assert( r[i]==static_cast<size_t>(-1) );
				r[i]=id;
				if (id!=i) sortneeded=true;
			}
			readtag(s,NumPart*sizeof(int)); // id end
			if (idmax!=NumPart || idmin!=1) throw_("bad ID array, not in [1;N] format");
		}
		catch (const Exception& e){
			errorloadingids=true;
		}

		if (sortneeded && !errorloadingids){
			if (debug) cerr << "[sorting]";

			vector<float> w(NumPart*3);
			assert( v.size()==w.size() &&  w.size()==3*r.size() );

			for(size_t i=0;i<NumPart;++i) {
				int from=i*3;
				int to  =r[i]*3;
				for(int j=0;j<3;++j) {
					if (from<0 || to<0){
						cerr << "** f=" << from << " t=" << to <<endl;
					}
					assert( from>=0 && to>=0 );
					if (!(static_cast<size_t>(to)<w.size() && static_cast<size_t>(from)<v.size())) throw_("id out of range, id=" + tostring(r[i]) + " to=" + tostring(to) + " from=" + tostring(from) + " size=" + tostring(v.size()) );
					w[to++]=v[from++];
				}
			}
			swap(v,w);
		}
		else if (debug && !errorloadingids){
			cerr << "[no sorting needed]";
			r.clear();
		}
		else if (errorloadingids){
			cerr << "[error loading ids]";
			r.clear();
		}
		if (debug) cerr << "...done" << endl;
	}

	if (debug) cerr << "## done" << endl;
	return make_pair(header1,r);
}
