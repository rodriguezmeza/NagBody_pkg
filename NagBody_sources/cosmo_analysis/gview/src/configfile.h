#ifndef __CONFIGFILE_H__
#define __CONFIGFILE_H__

class Configfile {
	private:
	typedef vector<string> t_value;
	t_value m_text;
	map<string,t_value> m_vals;
	string m_rem;

	void Load(istream &s)
	{
		if (!s) throw_("bad stream in configfile input operator");
		char buff[32*1024];
		bool ini=true;
		while(s.getline(buff,32*1024)) {
			const string b(strip(buff));
			if (ini && b!=m_rem + " Configfile_begin") throw_("bad stream in configfile input operator, missing begin tag");
			if (b==m_rem + " Configfile_end") return;
			const string t(strip(removerems(b,m_rem)));
			if (t.size()>0){
				m_text.push_back(t);
				const size_t n=t.find_first_of("=");
				istringstream s2(t.substr(0,n) + " " + (n==string::npos ? "" : t.substr(n+1)));
				string v,k;
				s2 >> k;
				while ((s2 >> v)) m_vals[k].push_back(v);
			}
			ini=false;
		}
		throw_("bad stream in configfile input operator, missing end tag");
	}

	public:
	Configfile() : m_rem("%") {}

	Configfile(const string& filename,const string rem="%") : m_rem(rem)
	{
		ifstream s(filename.c_str());
		if (!s) throw_("File <" + filename + "> does not exist");
		Load(s);
	}

	Configfile(istream& s,const string rem="%") : m_rem(rem)
	{
		if (!s) throw_("stream is invalid");
		Load(s);
	}

	size_t         size()                        const {return m_text.size();}
	const string&  operator[](const size_t n)    const {assert(n<m_text.size()); return m_text[n];}
	const t_value& operator()(const string& key) const {return Get(key);}
	bool           hasEntry  (const string& key) const {map<string,t_value >::const_iterator itt=m_vals.find(key); return itt!=m_vals.end();}
	void           Save(const string& filename) {ofstream s(filename.c_str()); s << *this;}

	const t_value& Get(const string& key) const
	{
		map<string,t_value >::const_iterator itt=m_vals.find(key);
		if (itt==m_vals.end()) throw_("No such entry, <" + key + ">, in configfile");
		return itt->second;
	}

	template<class T> const T Get(const string& key,const bool fullline=false)
	{
		const t_value& v=Get(key);
		assert( v.size()>0 );
		string s=v[0];
		if (fullline) for(size_t i=1;i<v.size();++i) s += " " + v[i];
		T t=totype<T>(s);
		return t;
	}

	template<class T> void Set(const string& key,const T& v)
	{
		t_value val;
		val.push_back(tostring(v));
		m_vals[key]=val;
		m_text.push_back(key + " " + val[0]);
	}

	friend ostream& operator<<(ostream& s,const Configfile& x)
	{
		s << x.m_rem << " Configfile_begin\n";
		for(map<string,t_value>::const_iterator itt1=x.m_vals.begin();itt1!=x.m_vals.end();++itt1){
			s << "\t" << itt1->first << " = ";
			for(t_value::const_iterator itt2=itt1->second.begin();itt2!=itt1->second.end();++itt2) s << *itt2 << " ";
			s << "\n";
		}
		s << x.m_rem << " Configfile_end\n";
		if (!s) throw_("bad stream in configfile output operator");
		return s;
	}

	friend istream& operator>>(istream& s,Configfile& x){x.Load(s); return s;}
};
#endif // __CONFIGFILE_H__
