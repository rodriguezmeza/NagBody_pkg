#ifndef __TRIPLE_H__
#define __TRIPLE_H__

template<class T>
class Triple
{
  private:
    T m_t[3];

  public:
    Triple() {}
    Triple(const T& t0,const T& t1,const T& t2)      {m_t[0]=t0,m_t[1]=t1,m_t[2]=t2;}
	Triple(const T& t)                               {m_t[0]=m_t[1]=m_t[2]=t;}
	Triple(const T  t[3])                            {m_t[0]=t[0];m_t[1]=t[1];m_t[2]=t[2];}
	Triple(const Triple<T>& org)                     {m_t[0]=org.m_t[0];m_t[1]=org.m_t[1];m_t[2]=org.m_t[2];}

    void      operator= (const Triple<T>& org)       {m_t[0]=org.m_t[0]; m_t[1]=org.m_t[1]; m_t[2]=org.m_t[2];}

	const T&  operator[](const size_t i) const       {assert(i<size()); return m_t[i];}
    T&        operator[](const size_t i)             {assert(i<size()); return m_t[i];}

    template<class R> Triple<T> operator* (const R& rhs) const {Triple<T> x=*this; x*=rhs; return x;}
    template<class R> Triple<T> operator+ (const R& rhs) const {Triple<T> x=*this; x+=rhs; return x;}
    template<class R> Triple<T> operator- (const R& rhs) const {Triple<T> x=*this; x-=rhs; return x;}
    template<class R> Triple<T> operator/ (const R& rhs) const {Triple<T> x=*this; x/=rhs; return x;}
    template<class R> void      operator*=(const R& rhs)       {for(size_t i=0;i<size();++i) m_t[i] *= rhs;}
    template<class R> void      operator+=(const R& rhs)       {for(size_t i=0;i<size();++i) m_t[i] += rhs;}
    template<class R> void      operator-=(const R& rhs)       {for(size_t i=0;i<size();++i) m_t[i] -= rhs;}
    template<class R> void      operator/=(const R& rhs)       {for(size_t i=0;i<size();++i) m_t[i] /= rhs;}

	Triple<T> operator* (const T& rhs)         const {Triple<T> x=*this; x*=rhs; return x;}
    Triple<T> operator+ (const T& rhs)         const {Triple<T> x=*this; x+=rhs; return x;}
    Triple<T> operator- (const T& rhs)         const {Triple<T> x=*this; x-=rhs; return x;}
    Triple<T> operator/ (const T& rhs)         const {Triple<T> x=*this; x/=rhs; return x;}
    void      operator*=(const T& rhs)               {for(size_t i=0;i<size();++i) m_t[i] *= rhs;}
    void      operator+=(const T& rhs)               {for(size_t i=0;i<size();++i) m_t[i] += rhs;}
    void      operator-=(const T& rhs)               {for(size_t i=0;i<size();++i) m_t[i] -= rhs;}
    void      operator/=(const T& rhs)               {for(size_t i=0;i<size();++i) m_t[i] /= rhs;}

	Triple<T> operator* (const Triple<T>& rhs) const {Triple<T> x=*this; x*=rhs; return x;}
    Triple<T> operator+ (const Triple<T>& rhs) const {Triple<T> x=*this; x+=rhs; return x;}
    Triple<T> operator- (const Triple<T>& rhs) const {Triple<T> x=*this; x-=rhs; return x;}
    Triple<T> operator/ (const Triple<T>& rhs) const {Triple<T> x=*this; x/=rhs; return x;}
    Triple<T> operator- ()                     const {Triple<T> x=*this; for(size_t i=0;i<3;++i) {x.m_t[i]=-m_t[i];} return x;}
    void      operator*=(const Triple<T>& rhs)       {for(size_t i=0;i<3;++i) {m_t[i]*=rhs.m_t[i];}}
    void      operator+=(const Triple<T>& rhs)       {for(size_t i=0;i<3;++i) {m_t[i]+=rhs.m_t[i];}}
    void      operator-=(const Triple<T>& rhs)       {for(size_t i=0;i<3;++i) {m_t[i]-=rhs.m_t[i];}}
    void      operator/=(const Triple<T>& rhs)       {for(size_t i=0;i<3;++i) {m_t[i]/=rhs.m_t[i];}}

	bool      operator==(const Triple<T>& x)  const  {return m_t[0]==x.m_t[0] && m_t[1]==x.m_t[1] && m_t[2]==x.m_t[2];}
    bool      operator!=(const Triple<T>& x)  const  {return !operator==(x);}

	static size_t size() {return 3;}

	const T& first () const {return m_t[0];}
	const T& second() const {return m_t[1];}
   	const T& third () const {return m_t[2];}
	T&       first ()       {return m_t[0];}
	T&       second()       {return m_t[1];}
   	T&       third ()       {return m_t[2];}

    const T& X() const {return m_t[0];}
    const T& Y() const {return m_t[1];}
    const T& Z() const {return m_t[2];}
    T&       X()       {return m_t[0];}
    T&       Y()       {return m_t[1];}
    T&       Z()       {return m_t[2];}

	#ifdef max
      #undef max
    #endif
    #ifdef min
      #undef min
    #endif

    void min(const Triple<T>& y) {for(size_t i=0;i<3;++i) if (m_t[i]>y.m_t[i]) m_t[i]=y.m_t[i];}
    void max(const Triple<T>& y) {for(size_t i=0;i<3;++i) if (m_t[i]<y.m_t[i]) m_t[i]=y.m_t[i];}

	void limitlo (const Triple<T>& lim) {for(size_t i=0;i<size();++i) if (m_t[i]< lim[i]) m_t[i]=lim[i];}
    void limithi (const Triple<T>& lim) {for(size_t i=0;i<size();++i) if (m_t[i]>=lim[i]) m_t[i]=lim[i]-1;}

	template<class R>
	void dist2(const Triple<T>& rhs,R& d) const
	{
		d=R(0);
      	for (size_t i=0;i<3;++i) {
        	const T x=m_t[i]-rhs.m_t[i];
        	d+=x*x;
      	}
    }

	template<class R>
    void length2(R& d) const
    {
      	d=R(0);
      	for (size_t i=0;i<3;++i) {
        	d += m_t[i]*m_t[i];
		}
    }

	template<class R> R dist2	(const Triple<T>& rhs) 	const {R d; dist2(rhs,d);return d;}
	template<class R> R length2	() 						const {R d; length2(d);  return d;}
    template<class R> T dist 	(const Triple<T>& rhs) 	const {return sqrt(dist2<R>(rhs));}
	template<class R> T length	()						const {return sqrt(length2<R>());}

	template<class R,class S>
	S dot(const Triple<R>& rhs) const
    {
      S d(0);
      for (size_t i=0;i<3;++i) {
        d += (*this)[i]*rhs[i];
      }
      return d;
    }

	template<class R,class S>
    Triple<S> cross(const Triple<R>& rhs) const
    {
      Triple<S> c;
	  c[0] = (*this)[2] * rhs[1] - (*this)[1] * rhs[2];
	  c[1] = (*this)[0] * rhs[2] - (*this)[2] * rhs[0];
	  c[2] = (*this)[1] * rhs[0] - (*this)[0] * rhs[1];
      return c;
    }

	T 		  dot  (const Triple<T>& rhs) const {return (*this).dot<T,T>(rhs);}
	Triple<T> cross(const Triple<T>& rhs) const {return (*this).cross<T,T>(rhs);}

    Triple<T> projection(const Triple<T>& b) const
    {
       const T d=this->dot(b);  // a dot b
       const T l=b.dot(b);   // length: |b|^2
       if (l==0) throw_("null lenght in projection");
       Triple<T> ab(b);
       ab *= d/l;
       return ab;
    }

    T lengthfromline(const Triple<T>& b,const Triple<T>& p) const
    {
      const Triple<T> linea=p-*this;
      const Triple<T> lineb=b-*this;

      const Triple<T> lineab=linea.projection(lineb);    // a projected onto b = ab
      const Triple<T> linec =lineab-linea;
      const T d=linec.dot(linec);
      assert(d>=0);
      const T distance=sqrt(d);
      assert(distance>=0);
      return distance;
    }

    Triple<T> topolar() // polar: X=r, Y=theta, Z=Z
    {
      const T theta = X!=0 ? atan(Y/X) : 0;
      return Triple<T>(dist(Triple<T>(0,0,0)),theta,Z);
    }

    Triple<T> frompolar()
    {
      return Triple<T>(X*cos(Y),X*sin(Y),Z);
    }

    Triple<T> tospherical() // spherical: X=r, Y=theta, Z=phi
    {
      const T r     = dist(Triple<T>(0,0,0));
      const T theta = X!=0 ? atan(Y/X) : 0;
      const T phi   = r!=0 ? acos(Z/r) : 0;
      return Triple<T>(r,theta,phi);
    }

    Triple<T> fromspherical()
    {
      return Triple<T> (X*sin(Z)*cos(Y),X*sin(Z)*sin(Y),X*cos(Z));
    }

	// indexing functions

	const T toidx(const T& sx,const T& sy,const T& sz) const
	{
		assert( numeric_limits<T>::is_integer==true );
		assert( numeric_limits<T>::is_signed==false );
		const T n=m_t[0]+(m_t[1]+m_t[2]*sz)*sy;
		assert( m_t[0]<sx && m_t[1]<sy && m_t[2]<sz && n<sx*sy*sz );
		return n;
	}

	static const Triple<T> toidx(const T& n,const T& sx,const T& sy,const T& sz)
	{
		assert( numeric_limits<T>::is_integer==true );
		assert( numeric_limits<T>::is_signed==false );
		assert( n<sx*sy*sz );

		const T nx=n % sx;
		const T ny=(n / sx) % sy;
		const T nz=n/(sx*sy);
		const Triple<T> idx(nx,ny,nz);
		assert( idx.toidx(sx,sy,sz)==n );
		return idx;
	}

	const               T  toidx(const T& s)            const {return toidx(s,s,s);}
	static const Triple<T> toidx(const T& n,const T& s)       {return toidx(n,s,s,s);}

	bool isinrange(const Triple<T>& idx1,const Triple<T>& idx2) const
	{
		assert( numeric_limits<T>::is_integer==true );
		assert( numeric_limits<T>::is_signed==false );

		for(unsigned int i=0;i<3;++i) {
			assert( idx1[i]!=idx2[i] );
			if (idx1[i]<idx2[i]){
				if (!(idx1[i]<=m_t[i] && m_t[i]<=idx2[i])) return false;
			} else {
				if (  idx2[i]< m_t[i] && m_t[i]< idx1[i] ) return false;
			}
		}
		return true;
	}

	template<class R> Triple<R> periodic(const R& period) const
	{
		assert( period>R(1) );
		Triple<R> idx_fixed_periodic;

		if (numeric_limits<T>::is_signed==false){
			for(size_t i=0;i<3;++i) idx_fixed_periodic[i] = (*this)[i] % period;
		}
		else{
			// T should be signed, R unsigned, both integers
			assert( numeric_limits<T>::is_integer==true );
			assert( numeric_limits<R>::is_integer==true );
			assert( numeric_limits<T>::is_signed==true );
			assert( numeric_limits<R>::is_signed==false );

			for(size_t i=0;i<3;++i){
				const bool ispos=(*this)[i]>=0;
				R x=ispos ? (*this)[i] : -(*this)[i];

				if (x>=period) x = x % period;
				assert( x<period );
				idx_fixed_periodic[i]= ispos ? x : x==0 ? 0 : period-x;
				assert( idx_fixed_periodic[i]>=0 && idx_fixed_periodic[i]<period );
			}
		}
		return idx_fixed_periodic;
	}

    friend ostream& operator<<(ostream& s,const Triple<T>& x){
       return s << "(" << x.X() << ";" << x.Y() << ";" << x.Z() << ")";
    }

	friend istream& operator>>(istream& s,Triple<T>& x){
		char c[4];
		s >> c[0] >> x.X() >> c[1] >> x.Y() >> c[2] >> x.Z() >> c[3];
		if (!s || c[0]!='(' || c[1]!=';' ||  c[2]!=';' ||  c[3]!=')') throw_("bad format in Triple<T> stream");
		return s;
    }
};

template<class T>
class Triplerange
{
	private:
	Triple<T> m_beg;
	Triple<T> m_end;
	T m_sz;

	public:
	Triplerange(const Triple<T>& begx,const Triple<T>& endx,const T sz=T(0))
	: m_beg(begx), m_end(endx), m_sz(sz)
	{
		assert( numeric_limits<T>::is_integer==true );
		assert( m_sz!=0 || (m_beg[0]<m_end[0] && m_beg[1]<m_end[1] && m_beg[2]<m_end[2]) );
		assert( m_sz==0 || (m_beg[0]<m_sz && m_beg[1]<m_sz && m_beg[2]<m_sz && m_end[0]<m_sz && m_end[1]<m_sz && m_end[2]<m_sz) );
	}

	const T            period() const {return m_sz;}
	const Triplerange& begin()  const {return *this;}
	const Triple<T>&   end  ()  const {return m_end;}
	const Triple<T>&   first()  const {return m_beg;}

    friend ostream& operator<<(ostream& s,const Triplerange<T>& x)
    {
       return s << "Triplerange: {" << x.m_beg << ";" << x.m_end << ";" << x.m_sz << "}";
	}

	class const_iterator
 	{
		private:
		Triplerange m_r;
		mutable Triple<T> m_curr;

		public:
		const_iterator(const Triplerange<T>& r)
		: m_r(r), m_curr(m_r.first())
		{
			assert( numeric_limits<T>::is_integer==true );
		}

		bool operator==(const Triple<T>& itt) const {return m_curr==itt;}
		bool operator!=(const Triple<T>& itt) const {return !((*this)==itt);}
		Triple<T>      operator*    () const {return index();}
		Triple<T>      index        () const {
			assert( m_curr[0]<m_r.period() || m_curr[0]<m_r.m_end[0]);
			assert( m_curr[1]<m_r.period() || m_curr[1]<m_r.m_end[1]);
			assert( m_curr[2]<m_r.period() || m_curr[2]<m_r.m_end[2]);
		    return m_curr;
	    }

		void operator++() const
		{
			++m_curr[0];
			if (m_r.period()!=0) m_curr[0] %= m_r.period();
			if (m_curr[0]==m_r.m_end[0]) {
				m_curr[0]=m_r.m_beg[0];
				++m_curr[1];
				if (m_r.period()!=0) m_curr[1] %= m_r.period();
				if (m_curr[1]==m_r.m_end[1]) {
					m_curr[1]=m_r.m_beg[1];
					++m_curr[2];
					if (m_r.period()!=0) m_curr[2] %= m_r.period();
					if (m_curr[2]==m_r.m_end[2]) {m_curr=m_r.m_end; return;}
				}
			}
			assert( m_curr[0]<m_r.period() || m_curr[0]<m_r.m_end[0]);
			assert( m_curr[1]<m_r.period() || m_curr[1]<m_r.m_end[1]);
			assert( m_curr[2]<m_r.period() || m_curr[2]<m_r.m_end[2]);
		}
	};
};

typedef Triple<float>  pointf;
typedef Triple<double> pointd;
typedef Triple<size_t> t_idx;

inline pointf tof(const pointd& p) {pointf t; t[0]=static_cast<float>(p[0]),t[1]=static_cast<float>(p[1]),t[2]=static_cast<float>(p[2]); return t;}
inline pointd tod(const pointf& p) {pointd t; t[0]=p[0],t[1]=p[1],t[2]=p[2]; return t;}

#endif // __TRIPLE_H__
