#ifndef __UTILS_H__
#define __UTILS_H__

extern bool DEBUG_PRINT;

inline string removerems(const string& s,const string rem)
{
    const size_t n=s.find_first_of(rem);
	return s.substr(0,n);
}

template<typename T>
inline string tostring(const T& x)
{
	ostringstream os;
	os << x;
	return os.str();
}

template<typename T>
inline T totype(const string& s)
{
	istringstream is(s);
	T x;
	is >> x;
	return x;
}

inline string suffix(const int n)
{
	assert(n>=0 && n<999);
	if (n<=9) return "00" + tostring(n);
	else if (n<=99) return "0" + tostring(n);
	else return tostring(n);
}

inline string strip(const string& s,const char ch=' ')
{
	const size_t n=s.find_first_not_of(ch);
	const size_t m=s.find_last_not_of (ch);

	if (n==string::npos || m==string::npos) return "";
	return s.substr(n,m-n+1);
}

inline bool FileExists(const std::string& f)
{
	ifstream s(f.c_str());
	return (!s)==false;
}

inline const vector<string> FindSnapshots(const string& filename)
{
	vector<string> v;
	const size_t n=filename.find_last_of('.');
	const string sfx = n==string::npos ? "" : filename.substr(n,filename.size());
	const string base= n==string::npos ? filename : filename.substr(0,n); 

	if (DEBUG_PRINT) cerr << "## Looking for file pattern:  base=<" << base << ">  sfx=<" << sfx << ">\n";

	
	if (base.size()<5 || base[base.size()-4]!='_') return v;
	const string base_x=base.substr(0,base.size()-4);
	for(int i=0;i<999;++i){
		const string f=base_x+"_"+suffix(i)+sfx;
		if (FileExists(f)) {
		    if (DEBUG_PRINT) cerr << "## found file <" << f << ">\n";
		    v.push_back(f);
		}
		else return v;
	}
	return v;
}

inline size_t FindSnapNumber(const string filename,const vector<string>& v)
{
	for(size_t i=0;i<v.size();++i) if (v[i]==filename) return i;
	return 0;
}

inline size_t GetThreadId()
{
	#ifdef WIN32
		assert( sizeof(size_t)==sizeof(void*) );
		return reinterpret_cast<size_t>(GetCurrentThreadId());
	#else
		// may be replaced by return 0; if phtread not found!
		assert( sizeof(size_t)==sizeof(pthread_t) );
		const pthread_t p=pthread_self();
		size_t q=0;
		memcpy(&q,&p,sizeof(q));
		return q;
	#endif
}

inline string System(const string& cmd,const bool throwexception=true,const bool captureoutput=false)
{
	if (!captureoutput){
		const int n=system(cmd.c_str());
		if (n!=0 && throwexception) throw_(string("system command failed with code=") + tostring(n) + " cmd=<" + cmd + ">");
		return "";
	} else {
		static size_t n=0;
		const string file("/tmp/tempfile." + tostring(GetThreadId()) + "."  + tostring(n++) ); // , or some random number, or process id
		ifstream s1(file.c_str());
		if(s1) {
			s1.close();
			System(("rm " + file).c_str(),true,false);
		}
		System(cmd + " > " + file,throwexception,false);

		string t;
		char buff[16*1024];
		ifstream s2(file.c_str());
		while(s2) {
			s2.getline(buff,16*1024);
			if (s2) t += buff;
		}
		System(("rm " + file).c_str(),true,false);
		return t;
	}
}

inline size_t FileTime(const string& file)
{
	if (!FileExists(file)) throw_("file does not exist: <" + file + ">");
	const string t=System("date -r " + file + " +%s",true,true); // seconds since 1970-01-01 00:00:00 UTC
	return totype<size_t>(t);
}

inline bool isFileNewer(const string& file0,const string& file1)
{
	return FileTime(file0)>FileTime(file1);
}

// dummy funs for sync purpose only
inline string replace(const string&,const string&,const string&){throw "unimplemented";}
template<class T> void Readbin(istream&,const T&) {throw "unimplemented";}



#endif // __UTILS_H__
