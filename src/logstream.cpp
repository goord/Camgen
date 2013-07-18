//
// This file is part of the CAMGEN library.
// Copyright (C) 2010 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <ctime>
#include <cstdlib>
#include <config.h>
#include <Camgen/logstream.h>

namespace Camgen
{
    end_msg::end_msg(){}

    end_msg endlog=end_msg();

    logstream::logstream():enable_level(log_level::message),level(log_level::message),streaming(true),writing(false),prompting(true),is(std::cin),os(std::cout),enabled(true){}

    logstream::logstream(const std::string& file):enable_level(log_level::message),level(log_level::message),streaming(false),writing(true),prompting(false),is(std::cin),os(std::cout),fs(file.c_str()),enabled(true)
    {
	init_logfile();
    }

    logstream::logstream(const std::string& file,std::istream& is_,std::ostream& os_):enable_level(log_level::message),level(log_level::message),streaming(true),writing(true),prompting(true),is(is_),os(os_),fs(file.c_str()),enabled(true)
    {
	init_logfile();
    }

    logstream::~logstream()
    {
	close();
    }

    logstream& logstream::open(const std::string& file)
    {
	if(fs.is_open())
	{
	    fs.close();
	}
	fs.open(file.c_str());
	init_logfile();
	return *this;
    }

    logstream& logstream::close()
    {
	if(fs.is_open())
	{
	    fs.close();
	}
	return *this;
    }

    logstream& logstream::operator()(log_level::type level_)
    {
	level=level_;
	switch(level)
	{
	    case log_level::message:
		(*this)<<"Message: ";
		break;
	    case log_level::warning:
		(*this)<<"Warning: ";
		break;
	    case log_level::error:
		(*this)<<"Error: ";
		break;
	    case log_level::abort:
		(*this)<<"Fatal error:";
		break;
	}
	return *this;
    }

    logstream& logstream::assert(bool condition)
    {
	enabled=condition;
	return *this;
    }

    bool logstream::logging() const
    {
	return enabled and level>=enable_level;
    }

    logstream& logstream::end_message()
    {
	if(streaming)
	{
	    os<<std::endl;
	}
	if(writing and fs.is_open())
	{
	    fs<<std::endl;
	}
	if(level==log_level::abort)
	{
	    close();
	    std::exit(EXIT_FAILURE);
	}
	if(prompting)
	{
	    std::cerr<<"continue (y(yes)/n(no)/c(stop prompting)) ? ";
	    char c;
	    is>>c;
	    if(c=='n')
	    {
		close();
		std::exit(EXIT_FAILURE);
	    }
	    if(c=='c')
	    {
		prompting=false;
	    }
	}
	level=log_level::message;
	return *this;
    }

    void logstream::init_logfile()
    {
	if(!writing)
	{
	    return;
	}
	if(fs.is_open())
	{
	    std::time_t rawtime;
	    struct tm* timeinfo;
	    std::time(&rawtime);
	    timeinfo=std::localtime(&rawtime);
	    fs<<"This is a log file by "<<PACKAGE_STRING<<", created on "<<std::asctime(timeinfo)<<std::endl;
	}
	else
	{
	    os<<"Error opening log file.....no file streaming will be performed."<<std::endl;
	}
    }

    logstream log("Camgen.log",std::cin,std::cout);

    logstream& operator << (logstream& ls,std::ostream& (*fp)(std::ostream&))
    {
	if(!ls.logging())
	{
	    return ls;
	}
	if(ls.streaming)
	{
	    fp(ls.os);
	}
	if(ls.writing and ls.fs.is_open())
	{
	    fp(ls.fs);
	}
	return ls;
    }

    logstream& operator << (logstream& ls,const end_msg& end)
    {
	if(!ls.logging())
	{
	    return ls;
	}
	return ls.end_message();
    }
}

