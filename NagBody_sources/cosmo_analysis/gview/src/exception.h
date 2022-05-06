#ifndef __EXCEPTION_H__
#define __EXCEPTION_H__

class Exception{
private:
	const string m_msg;
	const string m_file;
	const int m_line;

public:
	Exception(const string msg,const string file,const int line) : m_msg(msg), m_file(file), m_line(line) {}
	Exception(const char*  msg,const string file,const int line) : m_msg(msg), m_file(file), m_line(line) {}

	inline static string FormatCompilerMsg(const string& file,const int line,const bool warnonly=false)
	{
		#ifdef WIN32
			return file + "(" + Exception::tostring(line) + ") " + (warnonly ? "warning " : "error ");
		#else
			return file + ":" + Exception::tostring(line) + ": " + (warnonly ? "warning: " : "error: ");
		#endif
	}

	inline static void throw_fun(const string& msg,const string& file,const int line)
	{
		cout.flush();
		const string msg2=Exception::FormatCompilerMsg(file,line) + "throwing exception: " + msg;
		cerr << msg2 << "\n";
		cerr.flush();
		throw Exception(msg,file,line);
	}

	inline string Msg() const
	{
		return FormatCompilerMsg(m_file,m_line) + "Exception: " + m_msg;
	}

	friend ostream& operator<<(ostream& os,const Exception& e)
	{
		return os << e.Msg();
	}

	template<typename T>
	static string tostring(const T& x)
	{
		ostringstream os;
		os << x;
		return os.str();
	}
};

#define throw_(msg) Exception::throw_fun(msg, __FILE__, __LINE__)
#define warn_(msg) cerr << Exception::FormatCompilerMsg(__FILE__, __LINE__,true) << msg << endl

#define CATCH_ALL\
	catch(const char*   s)   {cout.flush(); cerr << Exception::FormatCompilerMsg(__FILE__, __LINE__) << "caught exception chars: " << s;}\
	catch(const string& s)   {cout.flush(); cerr << Exception::FormatCompilerMsg(__FILE__, __LINE__) << "caught exception string: " << s;}\
	catch(const Exception& s){cout.flush(); cerr << Exception::FormatCompilerMsg(__FILE__, __LINE__) << "caught Exception class: "  << s;}\
	catch(...) 				 {cout.flush(); cerr << Exception::FormatCompilerMsg(__FILE__, __LINE__) << "caught unknown exception";}\
	cerr << "...aborting" << endl;

#endif // __EXCEPTION_H__
