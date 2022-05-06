#ifndef __ARRAY_H__
#define __ARRAY_H__

template<class T>
class array3d
{
	private:
		size_t  m_nx,m_ny,m_nz,m_nt;
		size_t* m_ly;
		size_t* m_lyz;
		T*      m_d;

		inline bool CheckIdx(const size_t nx,const size_t ny,const size_t nz) const
		{
			if (nx<m_nx && ny<m_ny && nz<m_nz) return true;
			else return false;
		}

		inline size_t Idx(const size_t nx,const size_t ny,const size_t nz) const
		{
			assert( CheckIdx(nx,ny,nz) );
			const size_t idx=nx+m_ly[ny]+m_lyz[nz];
			assert( idx<m_nt);
			return idx;
		}

		inline const T& Get(const size_t nx,const size_t ny,const size_t nz) const {assert(m_d); return m_d[Idx(nx,ny,nz)];}
		inline       T& Get(const size_t nx,const size_t ny,const size_t nz)       {assert(m_d); return m_d[Idx(nx,ny,nz)];}

		t_idx Findminmax(const bool findmin) const
		{
			if (size()==0) throw_("Empty array");
			T mx=(*this)(0,0,0);
			t_idx mx_idx(0,0,0);
			for(size_t k=0;k<sizez();++k)
			for(size_t j=0;j<sizey();++j)
			for(size_t i=0;i<sizex();++i) {
				const T& v=(*this)(i,j,k);
				if (findmin ? mx>v : mx<v){
					mx=v;
					mx_idx=t_idx(i,j,k);
				}
			}
			return mx_idx;
		}

		void Delete()
		{
			if (m_d!=0)   {assert( m_ly!=0 && m_lyz!=0); delete[] m_d;}
			if (m_ly!=0)  delete[] m_ly;
			if (m_lyz!=0) delete[] m_lyz;
			m_d=0; m_ly=0; m_lyz=0;
			m_nx=m_ny=m_nz=m_nt=0;
		}

		void Init(const size_t nx,const size_t ny,const size_t nz,const bool allocate=true)
		{
			m_nx=nx,m_ny=ny,m_nz=nz,m_nt=nx*ny*nz;

			assert( m_nx<static_cast<size_t>(1<<31) );
			assert( m_ny<static_cast<size_t>(1<<31) );
			assert( m_nz<static_cast<size_t>(1<<31) );

			const int mem=static_cast<int>(1.0*m_nx*m_ny*m_nz*sizeof(T)/(1024*1024)+.5);
			if (mem>200) warn_(string("array3d<T> requiring a large memory footprint, memory at least=") + tostring(mem) + " Mb");

			assert(m_d==0);
			m_d=allocate ? new T[nx*ny*nz] : 0;
			m_ly=new size_t[ny];
			m_lyz=new size_t[nz];

			const size_t nxy=m_nx*m_ny;
			for(size_t i=0;i<m_ny;++i) m_ly [i]=m_nx*i;
			for(size_t i=0;i<m_nz;++i) m_lyz[i]=nxy*i;
		}

	public:
	array3d() 													: m_d(0) {Init(0,0,0);}
	array3d(const t_idx& sz)								 	: m_d(0) {Init(sz[0],sz[1],sz[2]);}
	array3d(const size_t nx,const size_t ny,const size_t nz,const bool allocate=true) : m_d(0) {Init(nx,ny,nz,allocate);}
	array3d(const array3d<T>& a) : m_d(0) {Init(a.sizex(),a.sizey(),a.sizez()); this->operator=(a);}
	~array3d() {Delete();}

	void operator=(const array3d<T>& a){
		resize(a.sizex(),a.sizey(),a.sizez());
		for(size_t t=0;t<size();++t) (*this)[t]=a[t];
	}

	void resize(const t_idx& sz){
		if(sz==sizet()) return;
		Delete();
		Init(sz[0],sz[1],sz[2]);
	}

	void resize(const size_t nx,const size_t ny,const size_t nz){resize(t_idx(nx,ny,nz));}
	void clear() {resize(t_idx(0,0,0));}

	const T& operator()(const size_t nx,const size_t ny,const size_t nz) const {return Get(nx,ny,nz);}
	      T& operator()(const size_t nx,const size_t ny,const size_t nz)       {return Get(nx,ny,nz);}

	const T& operator[](const size_t k) const {assert(k<m_nt && m_d); return m_d[k];}
	      T& operator[](const size_t k)       {assert(k<m_nt && m_d); return m_d[k];}

	const T& operator[](const t_idx& k) const {return (*this)(k[0],k[1],k[2]);}
	      T& operator[](const t_idx& k)       {return (*this)(k[0],k[1],k[2]);}

	void operator= (const T& x) {assert(m_d); for(size_t i=0;i<m_nt;++i) m_d[i]=x;}

	void operator*=(const T& x) {assert(m_d); for(size_t i=0;i<m_nt;++i) m_d[i]*=x;}
	void operator/=(const T& x) {assert(m_d); for(size_t i=0;i<m_nt;++i) m_d[i]/=x;}
	void operator+=(const T& x) {assert(m_d); for(size_t i=0;i<m_nt;++i) m_d[i]+=x;}
	void operator-=(const T& x) {assert(m_d); for(size_t i=0;i<m_nt;++i) m_d[i]-=x;}

	size_t sizex() const {return m_nx;}
	size_t sizey() const {return m_ny;}
	size_t sizez() const {return m_nz;}
	size_t size () const {return m_nt;}
	t_idx  sizet() const {return t_idx(sizex(),sizey(),sizez());}
	t_idx  toidx(const size_t idx) const {return t_idx::toidx(idx,m_nx,m_ny,m_nz);}
	T      sum  () const {T x=T(0); for(size_t i=0;i<m_nt;++i) x+= (*this)[i]; return x;}

	vector<vector<T> > slicez(const size_t k) const
	{
		vector<vector<T> > s;
		for(size_t j=0;j<sizey();++j){
			vector<T> ss;
			for(size_t i=0;i<sizex();++i) ss.push_back((*this)(i,j,k));
			s.push_back(ss);
		}
		return s;
	}

	vector<T> slice(const t_idx& begx,const t_idx& endx) const
	{
		vector<T> s;
		// 		const Triplerange<t_idx> r(beg,end);
		// 		for(Triplerange<t_idx>::const_iterator itt=r.begin();itt!=r.end();++itt){
		// 			const t_idx idx=itt.index();
		// 			s.push_back((*this)[idx]);
		// 		}
		for(size_t k=begx[2];k<endx[2];++k){
			for(size_t j=begx[1];j<endx[1];++j){
				for(size_t i=begx[0];i<endx[0];++i){
					s.push_back((*this)(i,j,k));
				}
			}
		}

		return s;
	}

	void transpose_slice(const size_t k)
	{
		for(size_t j=0;j<sizey();++j){
			for(size_t i=0;i<sizex();++i){
				swap((*this)(i,j,k),(*this)(j,i,k));
			}
		}
	}

	bool IndexOk(const t_idx& idx) const
	{
		if (idx[0]<sizex() && idx[1]<sizey() && idx[2]<sizez()) return true;
		else return false;
	}

	bool isSquare() const {return (sizex()==sizey() && sizex()==sizez());}

	t_idx FindMin() const {return Findminmax(true);}
	t_idx FindMax() const {return Findminmax(false);}

	friend ostream& operator<<(ostream& s,const array3d<T>& v)
	{
		s << "% AR3D(ascii) " << static_cast<unsigned long>(v.sizex()) << " " << static_cast<unsigned long>(v.sizey()) << " " << static_cast<unsigned long>(v.sizez()) << "\n";
		for(size_t k=0;k<v.sizez();++k){
			for(size_t j=0;j<v.sizey();++j){
				for(size_t i=0;i<v.sizex();++i){
					s << v(i,j,k)/T(1) << " ";
				}
				s << "\n";
			}
		}
		s << "% END\n";
		return s;
	}

	void Writeascii(const string& filename,const string comment="") const
	{
		ofstream s(filename.c_str());
		if (!s) throw_("bad write");
		if (comment.size()>0) s << "% " << replace(comment,"\n","\n%") << "\n";
		s << *this;
	}

	static array3d<T> Readascii(const string& filename)
	{
		ifstream s(filename.c_str());
		if (!s) throw_("bad read");

		string t1,t2;
		char c;

		// skip initial  comments
		while(s){
		    s >> t1 >> t2;
			while(t2=="%"){
				t1=t2;
				s >> t2;
			}
		    if (t1=="%" && t2=="AR3D(ascii)") break;
			Readline(s);
		}
		if (t2!="AR3D(ascii)") throw_("bad tag, not in array3d format");

		size_t x,y,z;
		x=y=z=0;
		s >> x >> y >> z;

		// skip more comments
		while(s){
			s >> c;
		    if (c!='%') { s.putback(c); break;}
			Readline(s);
		}

		array3d<T> a(x,y,z);
		for(size_t i=0;i<a.size();++i) {
			s >> a[i];
			if (!s) throw_("bad or truncated stream in array3d format");
		}
		s >> t1 >> t1;
		if (t1!="END") throw_("bad tag, array3d array truncated");
		return a;
	}

	static array3d<T> Readmatlab(const string& filename)
	{
		int columns=0,rows=0;
		{
			ifstream s(filename.c_str());
			if (!s) throw_("bad read, file <" + filename +  "> not found");

			// skip more comments
			char c;
			while(s){
				s >> c;
				if (c!='%') { s.putback(c); break;}
				Readline(s);
			}

			istringstream is(Readline(s));
			while(is) {T x; is >> x; if (is) {++columns; rows=1;}}
			while(s)  {Readline(s); if (s) ++rows;}
		}

		array3d<T> a(columns,rows,1);
		ifstream s(filename.c_str());

		// skip more comments
		char c;
		while(s){
			s >> c;
		    if (c!='%') { s.putback(c); break;}
			Readline(s);
		}

		for(size_t i=0;i<a.size();++i) {
			s >> a[i];
			if (!s) throw_("bad read" + tostring(i) );
		}
		return a;
	}

	void Writebin(const string& filename) const
	{
		ofstream s(filename.c_str());
		if (!s) throw_("bad write");
		Writetyp(s,"AR3D");
		Writetyp(s,sizex());
		Writetyp(s,sizey());
		Writetyp(s,sizez());
		assert(m_d!=0 || (m_d==0 && size()==0));
		s.write(reinterpret_cast<const char*>(m_d),static_cast<streamsize>(size()*sizeof(T)));
		if(!s) throw_( "bad write");
		Writetyp(s,"END");
	}

	static array3d<T> Readbin(const string& filename)
	{
		ifstream s(filename.c_str());
		if (!s) throw_("bad read");
		string t=Readstring(s);
		if (t!="AR3D") throw_("bad tag, not in array3d format");

		const size_t x=Readtyp<size_t>(s);
		const size_t y=Readtyp<size_t>(s);
		const size_t z=Readtyp<size_t>(s);
		array3d<T> a(x,y,z);

		char* c=reinterpret_cast<char*>(&a(0,0,0));
		s.read(c,static_cast<streamsize>(a.size()*sizeof(T)));
		if(!s) throw_( "bad Read");

		t=Readstring(s);
		if (t!="END") throw_("bad tag, array3d array truncated");

		return a;
	}

	static array3d<T> ReadFortran(ifstream& s,const Triple<size_t>& sz)
	{
		// Single dim 3D data:
		//   C:       X(0,0,0) X(1,0,0) ... X(szi-1,0,0) X(0,1,0)  ...  X(szi-1,szj-1,szk-1)
		//   Fortran: X(1,1,1)X(1,1,2)... X(1,1,szk) X(1,2,1) ... X(szi,szj,szk)
		// Triple 3D data:
		//   Fortran: handled by multiple files here: <file>x, <file>y, <file>z
		array3d<T> a(sz);

		for(size_t k=0;k<sz[2];++k){
			const unsigned int tag1=Readtyp<unsigned int>(s);
			for(size_t j=0;j<sz[1];++j){
				const vector<T> x=::Readbin<T>(s,sz[0]);
				for(size_t i=0;i<sz[1];++i) a(i,j,k)=x[i];
			}
			const unsigned int tag2=Readtyp<unsigned int>(s);
			if (tag1!=tag2) throw_("fortran stream error, bad header begin tag"); // 4 byte Fortran header begin
		}
		return a;
	}

		  array3d<T>& begin()       {return *this;}
		  size_t      end  ()       {return size();}
	const array3d<T>& begin() const {return *this;}
	const size_t      end  () const {return size();}

	struct iterator
	{
		private:
			size_t m_n;
			array3d<T>& m_a;
		public:
			iterator(array3d<T>& a) : m_n(0), m_a(a) {};

			void  operator=  (const size_t x) {m_n=x; assert(m_n<m_a.size());}
			void  operator++ ()       {assert(m_n<m_a.size()); ++m_n;}
		   	T&    operator*  ()       {assert(m_n<m_a.size()); return m_a[m_n];}
			const t_idx toidx() const {assert(m_n<m_a.size()); return m_a.toidx(m_n);}

			bool operator==(const size_t x) const {return m_n==x;}
			bool operator!=(const size_t x) const {return m_n!=x;}
	};

	struct const_iterator
	{
		private:
			size_t m_n;
			const array3d<T>& m_a;
		public:
			const_iterator(const array3d<T>& a) : m_n(0), m_a(a) {};

			void     operator= (const size_t x) {m_n=x; assert(m_n<m_a.size());}
			void     operator++()       {assert(m_n<m_a.size()); ++m_n;}
			const T& operator* () const {assert(m_n<m_a.size()); return m_a[m_n];}
			const t_idx  toidx () const {assert(m_n<m_a.size()); return m_a.toidx(m_n);}

			bool operator==(const size_t x) const {return m_n==x;}
			bool operator!=(const size_t x) const {return m_n!=x;}
	};

};

#endif // __ARRAY_H__
